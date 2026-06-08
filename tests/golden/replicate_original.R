#!/usr/bin/env Rscript
# tests/golden/replicate_original.R
# ==============================================================================
# SACRED GUARDRAIL. Runs the CURRENT code on the ORIGINAL anonimi cohort (77 raw)
# with NO validation split — the exact data conditions of oldresults (0.716) —
# and ASSERTS the headline metrics against tests/golden/expected_metrics.yml.
#
# This is the regression test that protects every refactor: if a change moves the
# replicated 0.7158 / perm-p 0.001 / nested-LOO 0.7542 / gate, this exits non-zero.
#
# Run:  Rscript tests/golden/replicate_original.R
# Pass = exit 0; Fail = exit 1 (prints every mismatch).
# Heavy (~8-9 min): full Step 01 -> 04 -> 06 with permutation test.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(yaml); library(here)
  library(openxlsx); library(jsonlite); library(readxl)
})
options(crayon.enabled = FALSE)
list.files(here("R"), pattern = "\\.R$", full.names = TRUE) %>% purrr::walk(source)
source(here("tests/golden/assert_metrics.R"))

expected <- yaml::read_yaml(here("tests/golden/expected_metrics.yml"))
EXP_NAME <- "BestResponse_2v3_4"

base_config <- yaml::read_yaml(here("config/global_params.yml"))
base_config$lmm_bootstrap$enabled <- FALSE        # match oldresults conditions
exp_cfg <- base_config$experiments[[EXP_NAME]]

# Force ORIGINAL files, NO split (replicate oldresults conditions exactly).
exp_cfg$validation_split <- NULL
exp_cfg$input_file    <- "data/Dati_NSCLC_standardizzati_anonimi_T0.xlsx"
exp_cfg$input_file_t0 <- "data/Dati_NSCLC_standardizzati_anonimi_T0.xlsx"
exp_cfg$input_file_t1 <- "data/Dati_NSCLC_standardizzati_anonimi_T1.xlsx"

run_root <- file.path(base_config$output_root,
                      paste0("golden_replicate_", format(Sys.time(), "%Y%m%d_%H%M%S")))
exp_out  <- file.path(run_root, EXP_NAME)
dir.create(exp_out, recursive = TRUE, showWarnings = FALSE)
message(sprintf("[GOLDEN] replicate_original -> %s", run_root))

# PASS 1: Step 01 standard
config <- base_config
config$project_name <- EXP_NAME; config$run_mode <- "standard"
config$clinical <- exp_cfg$clinical; config$features <- exp_cfg$features
config$is_longitudinal <- FALSE
config$input_file <- exp_cfg$input_file
config$output_root <- exp_out
source(here("src/01_data_processing.R"), echo = FALSE, local = FALSE)

# PASS 2: Step 01 long + Step 04 gate
config <- base_config
config$project_name <- EXP_NAME; config$run_mode <- "longitudinal"
config$clinical <- exp_cfg$clinical; config$features <- exp_cfg$features
config$is_longitudinal <- TRUE; config$qc$remove_outliers <- FALSE
config$input_file_t0 <- exp_cfg$input_file_t0
config$input_file_t1 <- exp_cfg$input_file_t1
config$output_root <- exp_out
source(here("src/01_data_processing.R"), echo = FALSE, local = FALSE)
source(here("src/04_longitudinal_lmm.R"), echo = FALSE, local = FALSE)

# PASS 3: Step 06 ML
config <- base_config
config$project_name <- EXP_NAME; config$run_mode <- "machine_learning"
config$clinical <- exp_cfg$clinical; config$features <- exp_cfg$features
config$input_file <- exp_cfg$input_file; config$input_file_t0 <- exp_cfg$input_file_t0
if (!is.null(exp_cfg$machine_learning))
  config$machine_learning <- modifyList(
    if (!is.null(config$machine_learning)) config$machine_learning else list(),
    exp_cfg$machine_learning)
config$output_root <- exp_out
source(here("src/06_machine_learning.R"), echo = FALSE, local = FALSE)

# ---- ASSERT against the frozen snapshot --------------------------------------
m <- extract_run_metrics(exp_out, EXP_NAME)
message("\n================ GOLDEN: original 77, no split ================")
fails <- report_assertions("original_cohort", m, expected$original_cohort)

if (length(fails) > 0) {
  message("\n[GOLDEN] FAILED — the sacred replication has drifted. Refactor is NOT safe.")
  quit(status = 1, save = "no")
}
message("\n[GOLDEN] PASSED — replication intact (0.7158 / 0.001 / 0.7542 / gate).")
quit(status = 0, save = "no")
