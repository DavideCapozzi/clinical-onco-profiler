# sweep_n_discovery.R
# ==============================================================================
# Sensitivity sweep: how does the BestResponse_2v3_4 (NSCLC v2) ML result move
# as the temporal discovery/validation split point (n_discovery) changes?
#
# For each n_discovery value it replicates the BestResponse passes that feed the
# LMM-gated Step 06 classifier:
#   Pass 0 (split) -> Pass 1 (Step 01 standard) -> Pass 2 (Step 01 long + Step 04)
#   -> Pass 3 (Step 06 ML)
# Steps 02/03/05 are skipped (not inputs to the lmm-gated Step 06).
#
# Each point writes to results/<sweep_id>/nd_<NN>/BestResponse_2v3_4/...
# A tidy summary is printed and written to results/<sweep_id>/SWEEP_SUMMARY.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
  library(here)
  library(openxlsx)
  library(jsonlite)
  library(readxl)
})
options(crayon.enabled = FALSE)

list.files(here("R"), pattern = "\\.R$", full.names = TRUE) %>% purrr::walk(source)

base_config <- yaml::read_yaml(here("config/global_params.yml"))
# Disable the Step 04 cluster bootstrap CI: it costs ~12 min/point and only feeds
# the forest-plot CIs, not the FDR+LOO gate that Step 06 consumes. Irrelevant to
# the AUC sensitivity this sweep measures.
base_config$lmm_bootstrap$enabled <- FALSE
EXP_NAME    <- "BestResponse_2v3_4"
exp_base    <- base_config$experiments[[EXP_NAME]]

# Discovery sizes to probe. QC pool = 85 -> valid range 1..84.
# `full` (NA) disables the split entirely (uses the whole v2 cohort).
ND_VALUES <- c(60, 65, 70, 73, 77, 81, 84, NA)

sweep_id   <- paste0("sweep_n_discovery_", format(Sys.time(), "%Y%m%d_%H%M%S"))
sweep_root <- file.path(base_config$output_root, sweep_id)
dir.create(sweep_root, recursive = TRUE, showWarnings = FALSE)
message(sprintf("[SWEEP] root: %s", sweep_root))

run_one <- function(nd) {
  is_full   <- is.na(nd)
  label     <- if (is_full) "full" else sprintf("%02d", nd)
  run_root  <- file.path(sweep_root, paste0("nd_", label))
  exp_out   <- file.path(run_root, EXP_NAME)
  dir.create(exp_out, recursive = TRUE, showWarnings = FALSE)
  message(sprintf("\n================ n_discovery = %s ================", label))

  exp_cfg <- exp_base
  n_disc_real <- NA_integer_; n_valid_real <- NA_integer_

  # ---- PASS 0: cohort split (skipped when full) ----------------------------
  if (!is_full) {
    split_config <- base_config
    if (!is.null(exp_cfg$clinical)) split_config$clinical <- exp_cfg$clinical
    if (!is.null(exp_cfg$features)) split_config$features <- exp_cfg$features

    split_cfg <- exp_cfg$validation_split
    split_cfg$n_discovery <- nd

    t0_path <- exp_cfg$input_file_t0
    t1_path <- exp_cfg$input_file_t1

    split <- compute_cohort_split(readxl::read_excel(t0_path), split_config, split_cfg)
    n_disc_real  <- sum(split$report$assignment == "discovery")
    n_valid_real <- sum(split$report$assignment == "validation")

    split_dir <- file.path(exp_out, "00_cohort_split")
    dir.create(split_dir, recursive = TRUE, showWarnings = FALSE)
    paths_t0 <- write_split_excels(t0_path, "T0", split, split_dir, EXP_NAME)
    write_split_excels(t1_path, "T1", split, split_dir, EXP_NAME)

    exp_cfg$input_file    <- paths_t0$discovery
    exp_cfg$input_file_t0 <- paths_t0$discovery
    exp_cfg$input_file_t1 <- file.path(split_dir, sprintf("%s_discovery_T1.xlsx", EXP_NAME))
  }

  # ---- PASS 1: Step 01 standard only --------------------------------------
  config <<- base_config
  config$project_name <<- EXP_NAME
  config$run_mode     <<- "standard"
  if (!is.null(exp_cfg$clinical)) config$clinical <<- exp_cfg$clinical
  if (!is.null(exp_cfg$features)) config$features <<- exp_cfg$features
  config$is_longitudinal <<- FALSE
  if (!is.null(exp_cfg$input_file)) config$input_file <<- exp_cfg$input_file
  config$output_root <<- exp_out
  source(here("src/01_data_processing.R"), echo = FALSE, local = FALSE)

  # ---- PASS 2: Step 01 longitudinal + Step 04 LMM gate --------------------
  config <<- base_config
  config$project_name <<- EXP_NAME
  config$run_mode     <<- "longitudinal"
  if (!is.null(exp_cfg$clinical)) config$clinical <<- exp_cfg$clinical
  if (!is.null(exp_cfg$features)) config$features <<- exp_cfg$features
  config$is_longitudinal <<- TRUE
  config$qc$remove_outliers <<- FALSE
  config$input_file_t0 <<- exp_cfg$input_file_t0
  config$input_file_t1 <<- exp_cfg$input_file_t1
  config$output_root <<- exp_out
  source(here("src/01_data_processing.R"), echo = FALSE, local = FALSE)
  source(here("src/04_longitudinal_lmm.R"), echo = FALSE, local = FALSE)

  # ---- PASS 3: Step 06 ML -------------------------------------------------
  config <<- base_config
  config$project_name <<- EXP_NAME
  config$run_mode     <<- "machine_learning"
  if (!is.null(exp_cfg$clinical)) config$clinical <<- exp_cfg$clinical
  if (!is.null(exp_cfg$features)) config$features <<- exp_cfg$features
  if (!is.null(exp_cfg$input_file))    config$input_file    <<- exp_cfg$input_file
  if (!is.null(exp_cfg$input_file_t0)) config$input_file_t0 <<- exp_cfg$input_file_t0
  if (!is.null(exp_cfg$machine_learning)) {
    config$machine_learning <<- modifyList(
      if (!is.null(config$machine_learning)) config$machine_learning else list(),
      exp_cfg$machine_learning)
  }
  config$output_root <<- exp_out
  source(here("src/06_machine_learning.R"), echo = FALSE, local = FALSE)

  # ---- Harvest metrics ----------------------------------------------------
  # Step 06 is SKIPPED (no ML json) when the LMM gate is empty — a meaningful
  # outcome, not an error. Read each json defensively.
  ml_json  <- file.path(exp_out, "06_machine_learning",
                        sprintf("Machine_Metrics_ML_%s.json", EXP_NAME))
  lmm_json <- file.path(exp_out, "04_longitudinal_analysis",
                        sprintf("Machine_Metrics_LMM_%s.json", EXP_NAME))
  ml  <- if (file.exists(ml_json))  jsonlite::fromJSON(ml_json,  simplifyVector = FALSE) else NULL
  lmm <- if (file.exists(lmm_json)) jsonlite::fromJSON(lmm_json, simplifyVector = FALSE) else NULL

  gate_mk <- if (is.null(ml)) "EMPTY (gate skipped)" else tryCatch(
    paste(sapply(ml$nested_loocv_validation$gate_stability, function(x) x$Marker), collapse = "+"),
    error = function(e) NA_character_)
  sig_lmm <- tryCatch({
    fr <- lmm$full_results
    mk <- sapply(fr, function(x) x$Marker)
    fd <- sapply(fr, function(x) x$FDR_Interaction)
    paste(mk[fd < 0.05], collapse = "+")
  }, error = function(e) NA_character_)

  tibble::tibble(
    n_discovery   = label,
    n_disc_pool   = n_disc_real,
    n_valid_pool  = n_valid_real,
    n_samples_ml  = if (is.null(ml)) NA else (ml$n_samples %||% NA),
    n_features    = if (is.null(ml)) NA else (ml$n_features_total %||% NA),
    primary       = if (is.null(ml)) NA else (ml$primary_method %||% NA),
    svm_auc       = if (is.null(ml)) NA else (ml$svm_rbf$metrics$auc %||% NA),
    svm_perm_p    = if (is.null(ml)) NA else (ml$permutation_test$svm_rbf$p_value %||% NA),
    en_auc        = if (is.null(ml)) NA else (ml$elastic_net$metrics$auc %||% NA),
    en_perm_p     = if (is.null(ml)) NA else (ml$permutation_test$elastic_net$p_value %||% NA),
    nested_loo    = if (is.null(ml)) NA else (ml$nested_loocv_validation$auc %||% NA),
    lmm_n_sig     = if (is.null(lmm)) NA else (lmm$significant_features_fdr %||% NA),
    lmm_sig_mk    = sig_lmm,
    gate_markers  = gate_mk
  )
}

results <- purrr::map_dfr(ND_VALUES, function(nd) {
  tryCatch(run_one(nd), error = function(e) {
    message(sprintf("[SWEEP] n_discovery=%s FAILED: %s", nd, e$message))
    tibble::tibble(n_discovery = if (is.na(nd)) "full" else sprintf("%02d", nd),
                   svm_auc = NA, svm_perm_p = NA)
  })
})

out_csv <- file.path(sweep_root, "SWEEP_SUMMARY.csv")
write.csv(results, out_csv, row.names = FALSE)
message("\n================ SWEEP SUMMARY ================")
print(as.data.frame(results))
message(sprintf("\n[SWEEP] summary written to %s", out_csv))
