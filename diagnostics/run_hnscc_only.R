# diagnostics/run_hnscc_only.R
# Run HNSCC_Response experiment only (Pass 1 standard + Pass 3 ML).
# Does not touch any NSCLC experiment results.

suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
  library(here)
  library(openxlsx)
})

options(crayon.enabled = FALSE)

base_config_path <- here("config/global_params.yml")
base_config      <- yaml::read_yaml(base_config_path)

message("\n>>> LOADING MODULES <<<")
list.files(here("R"), pattern = "\\.R$", full.names = TRUE) %>% purrr::walk(source)
message("[System] Modules loaded.")

# Generator runner: create and publish a fresh timestamped run, mirroring main.R.
run_id   <- make_run_id(base_config$run_label)
run_root <- file.path(base_config$output_root, run_id)
if (!dir.exists(run_root)) dir.create(run_root, recursive = TRUE)
message(sprintf("[System] Run ID: %s  ->  %s", run_id, run_root))

exp_name <- "HNSCC_Response"
exp_cfg  <- base_config$experiments[[exp_name]]

if (is.null(exp_cfg))
  stop(sprintf("[FATAL] Experiment '%s' not found in config.", exp_name))

run_machine_learning <- isTRUE(base_config$run_machine_learning)

message(sprintf("\n========================================================"))
message(sprintf("EXPERIMENT: %s", exp_name))
message(sprintf("========================================================\n"))

# ── PASS 1: Standard (Steps 01 → 02 → 03) ────────────────────────────────────
{
  config              <- base_config
  config$project_name <- exp_name
  config$run_mode     <- "standard"

  if (!is.null(exp_cfg$clinical)) config$clinical <- exp_cfg$clinical
  if (!is.null(exp_cfg$features)) {
    config$features <- exp_cfg$features
    if (length(unlist(config$features$facs)) == 0)
      stop(sprintf("[CONFIG] Experiment '%s': features.facs is empty", exp_name))
    if (length(unlist(config$features$soluble)) == 0)
      stop(sprintf("[CONFIG] Experiment '%s': features.soluble is empty", exp_name))
  }
  config$is_longitudinal <- FALSE
  if (!is.null(exp_cfg$input_file)) config$input_file <- exp_cfg$input_file

  config$output_root <- file.path(run_root, config$project_name)
  if (!dir.exists(config$output_root)) dir.create(config$output_root, recursive = TRUE)

  log_file <- file.path(config$output_root,
                        paste0("log_standard_hnscc_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))

  message(sprintf("--- [PASS 1] STANDARD PIPELINE: %s ---", config$project_name))
  cat(sprintf("=== STANDARD PIPELINE STARTED: %s ===\n", Sys.time()), file = log_file)

  tryCatch({
    withCallingHandlers({
      message(">>> PHASE 1: DATA PROCESSING <<<")
      source(here("src/01_data_processing.R"), echo = FALSE, local = FALSE)

      message("\n>>> PHASE 2: VISUALIZATION <<<")
      source(here("src/02_visualization.R"), echo = FALSE, local = FALSE)

      message("\n>>> PHASE 3: STATISTICAL ANALYSIS <<<")
      source(here("src/03_statistical_analysis.R"), echo = FALSE, local = FALSE)

      message(sprintf("\n[SUCCESS] Standard pipeline completed for: %s", exp_name))
    },
    message = function(m) cat(conditionMessage(m), file = log_file, append = TRUE, sep = "\n"),
    warning = function(w) cat(paste0("WARNING: ", conditionMessage(w)),
                              file = log_file, append = TRUE, sep = "\n"))
  }, error = function(e) {
    err_msg <- paste0("\n[FATAL ERROR] Standard pipeline stopped: ", e$message)
    cat(err_msg, file = log_file, append = TRUE, sep = "\n")
    message(err_msg)
  })
}

# ── PASS 3: Machine Learning (Step 06) ───────────────────────────────────────
if (run_machine_learning) {
  config              <- base_config
  config$project_name <- exp_name
  config$run_mode     <- "machine_learning"

  if (!is.null(exp_cfg$clinical))         config$clinical        <- exp_cfg$clinical
  if (!is.null(exp_cfg$features))         config$features        <- exp_cfg$features
  if (!is.null(exp_cfg$input_file))       config$input_file      <- exp_cfg$input_file
  if (!is.null(exp_cfg$input_file_t0))    config$input_file_t0   <- exp_cfg$input_file_t0
  if (!is.null(exp_cfg$machine_learning)) {
    config$machine_learning <- modifyList(
      if (!is.null(config$machine_learning)) config$machine_learning else list(),
      exp_cfg$machine_learning
    )
  }

  config$output_root <- file.path(run_root, config$project_name)
  if (!dir.exists(config$output_root)) dir.create(config$output_root, recursive = TRUE)

  log_file_ml <- file.path(config$output_root,
                            paste0("log_ml_hnscc_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))

  message(sprintf("\n--- [PASS 3] MACHINE LEARNING PIPELINE: %s ---", config$project_name))
  cat(sprintf("=== ML PIPELINE STARTED: %s ===\n", Sys.time()), file = log_file_ml)

  tryCatch({
    withCallingHandlers({
      message(">>> PHASE 6: MACHINE LEARNING CLASSIFICATION <<<")
      source(here("src/06_machine_learning.R"), echo = FALSE, local = FALSE)

      message(sprintf("\n[SUCCESS] ML pipeline completed for: %s", exp_name))
    },
    message = function(m) cat(conditionMessage(m), file = log_file_ml, append = TRUE, sep = "\n"),
    warning = function(w) cat(paste0("WARNING: ", conditionMessage(w)),
                              file = log_file_ml, append = TRUE, sep = "\n"))
  }, error = function(e) {
    err_msg <- paste0("\n[FATAL ERROR] ML pipeline stopped: ", e$message)
    cat(err_msg, file = log_file_ml, append = TRUE, sep = "\n")
    message(err_msg)
  })
}

# Publish provenance + the `latest` pointer for follow-up runners.
run_passes <- list()
run_passes[[exp_name]] <- c("standard", if (run_machine_learning) "machine_learning")
write_run_manifest(file.path(run_root, "run_manifest.yml"), run_id, base_config, run_passes)
update_latest_pointer(base_config$output_root, run_id)

message(sprintf("\n=== HNSCC_Response RUN COMPLETE (run: %s) ===", run_id))
