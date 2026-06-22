#!/usr/bin/env Rscript
# tests/integration_smoke.R
# ==============================================================================
# END-TO-END INTEGRATION SMOKE TEST
# ------------------------------------------------------------------------------
# Proves the "add a new cohort by adding a config block — no code changes"
# contract: builds a tiny SYNTHETIC longitudinal cohort (FACS-only) + a minimal
# config that points at it, runs the full pipeline (steps 01 → 04 → 06 incl. the
# clinical+immune added-value layer) via main.R with CONFIG_PATH, and asserts the
# run completes and emits the expected artifacts.
#
# It checks STRUCTURE, not numbers (those depend on RNG/data). It writes only to
# a tempdir output_root, so the real results/ and results/latest are untouched.
#
# Run:  Rscript tests/integration_smoke.R
# Exit: 0 on pass, non-zero on failure.
# ==============================================================================

suppressMessages({
  library(here)
  library(yaml)
  library(writexl)
})

set.seed(2026)

cat("[SMOKE] Building synthetic cohort fixtures...\n")

n_per_group <- 20
markers <- c("KI67NAIVE", "CD28KI67", "KI67CD4", "EFF", "EM", "NAIVE", "PD1", "CD3tot")
work    <- tempfile("smoke_run_"); dir.create(work)
t0_path <- file.path(work, "synth_T0.xlsx")
t1_path <- file.path(work, "synth_T1.xlsx")

# Two groups: responders (RP, code 1) start HIGHER on the KI67 markers and
# CONTRACT more by T1 — a clean Timepoint×Group signal so the LMM gate is found
# and the downstream classifier/added-value layer have something to model.
make_frame <- function(timepoint) {
  ids   <- sprintf("SYN_%03d", seq_len(2 * n_per_group))
  group <- rep(c(1L, 0L), each = n_per_group)            # 1 = RP, 0 = SD_PD
  df <- data.frame(Patient_ID = ids, stringsAsFactors = FALSE)
  for (m in markers) {
    base <- 10 + 8 * group                               # responders higher at baseline
    drop <- if (timepoint == "T1") 6 * group else 0      # responders contract by T1
    df[[m]] <- pmax(0.5, base - drop + rnorm(length(ids), 0, 3))
  }
  # Per-patient outcome + clinical model vars + survival (constant across timepoints)
  df[["resp_code"]] <- group                             # 1 = RP, 0 = SD_PD
  df[["clin_a"]]    <- 30 * group + rnorm(length(ids), 0, 15)   # PD-L1-like
  df[["clin_b"]]    <- 3 - 1.5 * group + rnorm(length(ids), 0, 1) # NLR-like
  df[["OS"]]        <- round(runif(length(ids), 3, 40), 1)
  df[["PFS"]]       <- round(runif(length(ids), 2, 30), 1)
  df[["death"]]     <- rbinom(length(ids), 1, 0.5)
  df
}

writexl::write_xlsx(make_frame("T0"), t0_path)
writexl::write_xlsx(make_frame("T1"), t1_path)

# --- Build a minimal config from the real one (inherits qc/transform/stats/etc.),
#     replacing experiments with one synthetic block and lowering all cost knobs.
cfg <- yaml::read_yaml(here("config/global_params.yml"))

cfg$output_root           <- file.path(work, "out")
cfg$input_file            <- t0_path
cfg$input_file_t0         <- t0_path
cfg$input_file_t1         <- t1_path
cfg$run_machine_learning  <- TRUE
cfg$run_clinical_utility  <- TRUE
cfg$run_publication_figures <- FALSE          # skip figure rendering for speed/deps

# Cost knobs → tiny (structure test, not statistics)
cfg$lmm_bootstrap$enabled         <- FALSE
cfg$multivariate$n_repeat_cv      <- 5
cfg$multivariate$permutation$n_perm <- 19
cfg$machine_learning$n_perm       <- 19
cfg$machine_learning$clinical_utility$n_boot <- 30

cfg$experiments <- list(
  Synthetic = list(
    is_longitudinal = TRUE,
    run_network     = FALSE,
    input_file      = t0_path,
    input_file_t0   = t0_path,
    input_file_t1   = t1_path,
    features = list(facs = as.list(markers), soluble = list()),
    clinical = list(
      target_column       = "resp_code",
      responder_label     = "RP",
      non_responder_label = "SD_PD",
      mapping             = list(RP = list(1), SD_PD = list(0)),
      benchmark_column    = "clin_a",
      benchmark_label     = "synthetic benchmark"
    ),
    machine_learning = list(
      clinical_model = list(
        vars            = list(A = "clin_a", B = "clin_b"),
        impute          = "median",
        calibration_ridge = TRUE,
        specificity_null  = list(n_random = 20),
        export_locked_model = TRUE,
        survival        = list(os = "OS", pfs = "PFS", death = "death")
      )
    )
  )
)

cfg_path <- file.path(work, "smoke_config.yml")
yaml::write_yaml(cfg, cfg_path)

cat(sprintf("[SMOKE] Running pipeline via main.R (CONFIG_PATH=%s)...\n", cfg_path))

rscript <- file.path(R.home("bin"), "Rscript")
status <- system2(
  rscript, args = here("main.R"),
  env = c(sprintf("CONFIG_PATH=%s", cfg_path),
          "EXPERIMENTS=Synthetic",
          "RUN_LABEL=smoke"),
  stdout = file.path(work, "run.log"), stderr = file.path(work, "run.log")
)

fail <- function(msg) {
  cat(sprintf("\n[SMOKE] FAIL: %s\n", msg))
  cat("---- last 25 log lines ----\n")
  lg <- file.path(work, "run.log")
  if (file.exists(lg)) cat(paste(tail(readLines(lg, warn = FALSE), 25), collapse = "\n"), "\n")
  quit(status = 1)
}

if (status != 0) fail(sprintf("main.R exited with status %d", status))

# --- Assert expected artifacts exist (structure, not numbers) ---
run_id  <- readLines(file.path(cfg$output_root, "latest_run.txt"), warn = FALSE)[1]
exp_dir <- file.path(cfg$output_root, run_id, "Synthetic")

expected <- c(
  proc = list.files(file.path(exp_dir, "01_data_processing"), pattern = "\\.rds$", full.names = TRUE)[1],
  lmm  = file.path(exp_dir, "04_longitudinal_analysis", "Machine_Metrics_LMM_Synthetic.json"),
  ml   = file.path(exp_dir, "06_machine_learning", "Machine_Metrics_ML_Synthetic.json")
)

for (nm in names(expected)) {
  p <- expected[[nm]]
  if (is.na(p) || !file.exists(p)) fail(sprintf("missing expected artifact '%s' (%s)", nm, p))
}

# ML JSON must be valid and carry the added-value node (the new-cohort contract)
ml <- tryCatch(jsonlite::read_json(expected[["ml"]], simplifyVector = TRUE),
               error = function(e) fail(sprintf("ML JSON unreadable: %s", e$message)))
if (is.null(ml$clinical_immune_added_value))
  fail("ML JSON missing clinical_immune_added_value node (added-value layer did not run)")

cat("\n[SMOKE] PASS — pipeline ran end-to-end on a synthetic config and emitted:\n")
cat(sprintf("  - %s\n", expected[["proc"]]))
cat(sprintf("  - %s\n", expected[["lmm"]]))
cat(sprintf("  - %s (added-value node present)\n", expected[["ml"]]))
quit(status = 0)
