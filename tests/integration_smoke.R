#!/usr/bin/env Rscript
# tests/integration_smoke.R
# ==============================================================================
# END-TO-END INTEGRATION SMOKE TEST
# ------------------------------------------------------------------------------
# Proves the "add a new cohort by adding a config block — no code changes"
# contract: builds a tiny SYNTHETIC longitudinal cohort (FACS-only) + a minimal
# config that points at it, runs the full pipeline (steps 01 → 04 → 06 incl. the
# clinical+immune added-value layer, the selection-aware LRT, the patient
# bootstrap and the LRT leave-one-out, then step 07 post hoc and the publication
# collector) via main.R with CONFIG_PATH + RUN_NAME, and asserts the run completes,
# logs no [FATAL ERROR] and emits the expected artifacts.
#
# It checks STRUCTURE, not numbers (those depend on RNG/data). It writes only to
# a tempdir output_root, so the real results/ are untouched.
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

# Two groups: responders (RP) start HIGHER on the two KI67 markers and CONTRACT more
# by T1 — a clean Timepoint×Group signal on those two only, so the LMM gate is the
# two-marker gate (as in the real data) and the downstream layers have something to
# model. The other markers carry no group signal.
make_frame <- function(timepoint) {
  ids   <- sprintf("SYN_%03d", seq_len(2 * n_per_group))
  group <- rep(c(1L, 0L), each = n_per_group)            # 1 = RP, 0 = SD_PD
  df <- data.frame(Patient_ID = ids, stringsAsFactors = FALSE)
  signal <- c("KI67NAIVE", "CD28KI67")
  for (m in markers) {
    g    <- if (m %in% signal) group else 0L
    base <- 10 + 8 * g                                   # responders higher at baseline
    drop <- if (timepoint == "T1") 6 * g else 0          # responders contract by T1
    df[[m]] <- pmax(0.5, base - drop + rnorm(length(ids), 0, 3))
  }
  # Per-patient outcome + clinical model vars + survival (constant across timepoints).
  # Best-response codes as in the real data: 2 = PR (RP), 3 = SD, 4 = PD.
  df[["resp_code"]] <- ifelse(group == 1L, 2L, rep(c(3L, 4L), length.out = length(ids)))
  df[["clin_a"]]    <- 30 * group + rnorm(length(ids), 0, 15)   # PD-L1-like
  df[["clin_b"]]    <- 3 - 1.5 * group + rnorm(length(ids), 0, 1) # NLR-like
  df[["OS"]]        <- round(runif(length(ids), 3, 40), 1)
  df[["PFS"]]       <- round(runif(length(ids), 2, 30), 1)
  df[["death"]]     <- rbinom(length(ids), 1, 0.5)
  df
}

writexl::write_xlsx(make_frame("T0"), t0_path)
writexl::write_xlsx(make_frame("T1"), t1_path)

# Therapy sheets for the step-07 therapy block (same layout as data/therapy_fix/).
th_paths <- file.path(work, sprintf("synth_therapy_%s.xlsx", c("T0", "T1")))
for (p in th_paths) {
  ids <- sprintf("SYN_%03d", seq_len(2 * n_per_group))
  writexl::write_xlsx(data.frame(
    Patient_ID = ids,
    `Terapia_Pembro_1_Pembro_CHT_2_ipi/nivo_CHT_3;_Atezo_CHT_4` = rep(c(1L, 2L), length.out = length(ids)),
    Terapia_regime_IOmono_1_chemoIO_2 = rep(c(1L, 2L), length.out = length(ids)),
    Terapia_radio_tra_T0_T1_annotata = 0L, check.names = FALSE), p)
}

# --- Build a minimal config from the real one (inherits qc/transform/stats/etc.),
#     replacing experiments with one synthetic block and lowering all cost knobs.
cfg <- yaml::read_yaml(here("config/global_params.yml"))

cfg$output_root           <- file.path(work, "out")
cfg$input_file            <- t0_path
cfg$input_file_t0         <- t0_path
cfg$input_file_t1         <- t1_path
cfg$run_machine_learning  <- TRUE
cfg$run_clinical_utility  <- TRUE
cfg$run_publication_figures <- TRUE           # the collector needs the figure set

# Cost knobs → tiny (structure test, not statistics)
cfg$lmm_bootstrap$enabled         <- FALSE
cfg$multivariate$n_repeat_cv      <- 5
cfg$multivariate$permutation$n_perm <- 19
cfg$machine_learning$n_perm       <- 19
cfg$machine_learning$clinical_utility$n_boot <- 30
cfg$machine_learning$clinical_immune_added_value <- list(
  selection_aware_lrt = list(enabled = TRUE, B = 20, topk = "gate_size", cores = 2, seed = 2026),
  patient_bootstrap   = list(enabled = TRUE, B = 20, reps_in_boot = 2, cores = 2, seed = 2026),
  lrt_leave_one_out   = list(enabled = TRUE))
cfg$post_hoc <- list(
  grouping       = list(enabled = TRUE, subcodes = list(PR = 2, SD = 3, PD = 4)),
  duplicate_rows = list(enabled = TRUE, drop_t1 = list("SYN_001", "SYN_021")),
  therapy        = list(enabled = TRUE, file_t0 = th_paths[1], file_t1 = th_paths[2],
                        tag = "synthetic", B = 100, seed = 56),
  km             = list(enabled = TRUE, landmark_days = 21))

cfg$experiments <- list(
  Synthetic = list(
    is_longitudinal = TRUE,
    input_file      = t0_path,
    input_file_t0   = t0_path,
    input_file_t1   = t1_path,
    features = list(facs = as.list(markers), soluble = list()),
    clinical = list(
      target_column       = "resp_code",
      responder_label     = "RP",
      non_responder_label = "SD_PD",
      mapping             = list(RP = list(2), SD_PD = list(3, 4)),
      sensitivity_covariates = list("clin_a"),
      benchmark_column    = "clin_a",
      benchmark_label     = "synthetic benchmark"
    ),
    machine_learning = list(
      clinical_model = list(
        vars            = list(A = "clin_a", B = "clin_b"),
        impute          = "median",
        composite_timepoint = "delta",
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
          "RUN_NAME=smoke"),
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
run_dir <- file.path(cfg$output_root, "smoke")
exp_dir <- file.path(run_dir, "Synthetic")

# main.R exits 1 on a failed pass; the log grep also names the failing log.
logs <- c(file.path(work, "run.log"), list.files(exp_dir, pattern = "^log_.*\\.txt$", full.names = TRUE),
          list.files(run_dir, pattern = "^log_.*\\.txt$", full.names = TRUE))
fatal <- logs[vapply(logs, function(f) any(grepl("FATAL ERROR", readLines(f, warn = FALSE), fixed = TRUE)), logical(1))]
if (length(fatal)) fail(sprintf("[FATAL ERROR] in %s", paste(basename(fatal), collapse = ", ")))

ml_dir <- file.path(exp_dir, "06_machine_learning")
ph_dir <- file.path(exp_dir, "07_post_hoc")
expected <- c(
  proc     = list.files(file.path(exp_dir, "01_data_processing"), pattern = "\\.rds$", full.names = TRUE)[1],
  lmm      = file.path(exp_dir, "04_longitudinal_analysis", "Machine_Metrics_LMM_Synthetic.json"),
  ml       = file.path(ml_dir, "Machine_Metrics_ML_Synthetic.json"),
  sel_rds  = file.path(ml_dir, "SelectionAware_LRT_Synthetic.rds"),
  sel_csv  = file.path(ml_dir, "SelectionAware_LRT_Synthetic.csv"),
  boot_csv = file.path(ml_dir, "PatientBootstrap_AUC_Synthetic.csv"),
  loo_csv  = file.path(ml_dir, "LRT_LeaveOneOut_Synthetic.csv"),
  grouping = file.path(ph_dir, "grouping", "grouping_lmm_all8_unadjusted.csv"),
  dup      = file.path(ph_dir, "duplicate_rows", "duplicate_rows_scenarios.csv"),
  therapy  = file.path(ph_dir, "therapy", "therapy_tests.csv"),
  km       = file.path(ph_dir, "km", "km_tests.csv"),
  manifest = file.path(run_dir, "run_manifest.yml"),
  index    = file.path(run_dir, "publication", "INDEX.md")
)

for (nm in names(expected)) {
  p <- expected[[nm]]
  if (is.na(p) || !file.exists(p) || file.size(p) == 0) fail(sprintf("missing expected artifact '%s' (%s)", nm, p))
}

# ML JSON must be valid and carry the added-value node (the new-cohort contract)
# and the three inference nodes.
ml <- tryCatch(jsonlite::read_json(expected[["ml"]], simplifyVector = TRUE),
               error = function(e) fail(sprintf("ML JSON unreadable: %s", e$message)))
if (is.null(ml$clinical_immune_added_value))
  fail("ML JSON missing clinical_immune_added_value node (added-value layer did not run)")
for (node in c("selection_aware_lrt", "patient_bootstrap", "lrt_leave_one_out"))
  if (is.null(ml[[node]])) fail(sprintf("ML JSON missing the %s node", node))

# Publication set: every mapped main / supplementary / table entry was collected,
# named by its paper label (e.g. main/Fig_1_Biology.pdf).
source(here("R/utils_publication.R"))   # pub_dest_stem()
pub <- file.path(run_dir, "publication")
for (role in c("main", "supplementary", "tables")) {
  sub <- c(main = "main", supplementary = "supplementary", tables = "supplementary/tables")[[role]]
  for (e in cfg$publication[[role]]) {
    hit <- list.files(file.path(pub, sub), pattern = sprintf("^%s(_p[0-9]+)?\\.", pub_dest_stem(e)))
    if (!length(hit)) fail(sprintf("publication/%s: %s (%s) not collected", sub, pub_dest_stem(e), e$stem))
  }
}

# Manifest: frame Ns and the config md5 are recorded.
man <- yaml::read_yaml(expected[["manifest"]])
if (is.null(man$config_md5) || is.null(man$experiments$Synthetic$n_frames$delta))
  fail("run_manifest.yml lacks config_md5 or the frame Ns")
if (!"failed_passes" %in% names(man) || length(man$failed_passes))
  fail("run_manifest.yml: failed_passes missing or not empty")

cat("\n[SMOKE] PASS — pipeline ran end-to-end on a synthetic config and emitted:\n")
for (nm in names(expected)) cat(sprintf("  - %s\n", expected[[nm]]))
quit(status = 0)
