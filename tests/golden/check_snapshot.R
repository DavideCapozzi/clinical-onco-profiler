#!/usr/bin/env Rscript
# tests/golden/check_snapshot.R
# ==============================================================================
# ANTI-DRIFT SNAPSHOT (adjustment #1). Re-runs the CURRENT code for the two
# experiments whose numbers are publishable / quoted in CLAUDE.md —
# BestResponse_2v3_4 (v2, discovery-73 split) and HNSCC_Response — and ASSERTS
# their headline metrics (perm p, nested-LOO, FDR count, gate, AUC) against the
# frozen tests/golden/expected_metrics.yml.
#
# Together with replicate_original.R this guarantees no refactor silently moves a
# publishable number.
#
# Run:  Rscript tests/golden/check_snapshot.R
# Pass = exit 0; Fail = exit 1. Heavy (~12-15 min). NOTE: this runs main.R and so
# updates the results/latest pointer (run id suffixed `_golden_snapshot`).
# ==============================================================================

suppressPackageStartupMessages({ library(yaml); library(here); library(jsonlite) })
source(here("R/utils_io.R"))
source(here("R/utils_metrics.R"))
source(here("tests/golden/assert_metrics.R"))

expected    <- yaml::read_yaml(here("tests/golden/expected_metrics.yml"))
base_config <- yaml::read_yaml(here("config/global_params.yml"))
EXPS <- c("BestResponse_2v3_4", "HNSCC_Response")

# --- Re-run current code for exactly these experiments ------------------------
rscript <- file.path(R.home("bin"), "Rscript")
message(sprintf("[GOLDEN] running main.R for: %s", paste(EXPS, collapse = ", ")))
rc <- system2(rscript, args = here("main.R"),
              env = c(sprintf("EXPERIMENTS=%s", paste(EXPS, collapse = ",")),
                      "RUN_LABEL=golden_snapshot"),
              stdout = "", stderr = "")
if (rc != 0) { message(sprintf("[GOLDEN] main.R exited %d", rc)); quit(status = 1, save = "no") }

# --- Assert each experiment against the frozen snapshot -----------------------
run_root <- resolve_run_root(base_config$output_root)
message(sprintf("\n================ GOLDEN SNAPSHOT: %s ================", basename(run_root)))

all_fails <- character(0)
for (exp in EXPS) {
  if (is.null(expected[[exp]])) {
    message(sprintf("  [SKIP] %s — no frozen block in expected_metrics.yml", exp)); next
  }
  m <- extract_run_metrics(file.path(run_root, exp), exp)
  all_fails <- c(all_fails, report_assertions(exp, m, expected[[exp]]))
}

if (length(all_fails) > 0) {
  message(sprintf("\n[GOLDEN] FAILED — %d publishable metric(s) drifted.", length(all_fails)))
  quit(status = 1, save = "no")
}
message("\n[GOLDEN] PASSED — v2 + HNSCC headline metrics intact.")
quit(status = 0, save = "no")
