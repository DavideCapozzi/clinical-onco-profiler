# main.R
# ==============================================================================
# MAIN PIPELINE ORCHESTRATOR: clinical-onco-profiler
# Description: End-to-end execution for Hybrid Data (FACS + Solubles)
#              Goal: Identify drivers of clinical response via sPLS-DA & LMM.
#              Supports nested multi-experiment and flat single-cohort configs.
# ==============================================================================

# 1. Environment Setup
rm(list = ls())
graphics.off()

suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
  library(here)
  library(openxlsx)
})

options(crayon.enabled = FALSE)

# Load BASE Configuration 
base_config_path <- here("config/global_params.yml")
if (!file.exists(base_config_path)) stop("[FATAL] Base Config file not found at expected path.")
base_config <- yaml::read_yaml(base_config_path)

# --- Configuration Router (Nested vs Flat Support) ---
experiments_list <- list()

if (!is.null(base_config$experiments) && length(base_config$experiments) > 0) {
  # Mode: Nested Multi-Experiment
  experiments_list <- base_config$experiments
  message("[System] Detected nested multi-experiment configuration.")
} else {
  # Mode: Flat Single-Cohort Configuration (Fallback)
  exp_name <- if (!is.null(base_config$project_name)) base_config$project_name else "Single_Cohort_Run"
  experiments_list[[exp_name]] <- base_config
  message("[System] Detected flat single-cohort configuration.")
}

# Optional experiment filter: EXPERIMENTS="A,B" restricts this run to the named
# experiments (e.g. for the golden snapshot test). Unset -> run all (default).
exp_filter <- Sys.getenv("EXPERIMENTS", unset = "")
if (nzchar(exp_filter)) {
  wanted <- trimws(strsplit(exp_filter, ",")[[1]])
  missing <- setdiff(wanted, names(experiments_list))
  if (length(missing) > 0)
    stop(sprintf("[System] EXPERIMENTS names not found in config: %s", paste(missing, collapse = ", ")))
  experiments_list <- experiments_list[wanted]
  message(sprintf("[System] EXPERIMENTS filter active -> %s", paste(wanted, collapse = ", ")))
}

# --- Module Loading ---
message("\n>>> LOADING MODULES <<<")
list.files(here("R"), pattern = "\\.R$", full.names = TRUE) %>% purrr::walk(source)
message("[System] Modules loaded successfully.")
validate_config(base_config)

# Isolate this invocation in one immutable, timestamped run directory. Every
# pass routes its output_root under run_root, so all cross-step artifacts of a
# single main.R call stay self-contained and reproducible.
run_id   <- make_run_id(base_config$run_label)
run_root <- file.path(base_config$output_root, run_id)
if (!dir.exists(run_root)) dir.create(run_root, recursive = TRUE)
message(sprintf("[System] Run ID: %s  ->  %s", run_id, run_root))

run_passes <- list()  # experiment -> passes executed (for the run manifest)

# 2. Pipeline Execution Loop
# ------------------------------------------------------------------------------
for (exp_name in names(experiments_list)) {

  exp_cfg <- experiments_list[[exp_name]]

  # Inherit boolean flags with safety fallbacks
  run_standard        <- if (!is.null(exp_cfg$run_standard))         as.logical(exp_cfg$run_standard)         else TRUE
  run_longitudinal    <- if (!is.null(exp_cfg$is_longitudinal))      as.logical(exp_cfg$is_longitudinal)      else FALSE
  run_network         <- if (!is.null(exp_cfg$run_network))          as.logical(exp_cfg$run_network)          else FALSE
  # ML is a global flag (base_config), not per-experiment — Step 06 self-gates on LOO-robust feature availability
  run_machine_learning <- if (!is.null(base_config$run_machine_learning)) as.logical(base_config$run_machine_learning) else FALSE

  run_passes[[exp_name]] <- c(
    if (run_standard)         "standard",
    if (run_network)          "network",
    if (run_longitudinal)     "longitudinal",
    if (run_machine_learning) "machine_learning"
  )

  # ============================================================================
  # PASS 0: COHORT SPLIT (discovery / validation)
  # Runs once per experiment when validation_split is enabled. Materializes the
  # discovery + held-out validation Excel snapshots, then redirects every
  # downstream pass to the discovery files. Fatal on failure by design: a broken
  # split must never silently fall back to the full cohort.
  # ============================================================================
  split_cfg <- exp_cfg$validation_split
  if (!is.null(split_cfg) && isTRUE(split_cfg$enabled)) {
    message(sprintf("\n--- [PASS 0] COHORT SPLIT: %s ---", exp_name))

    split_config <- base_config
    if (!is.null(exp_cfg$clinical)) split_config$clinical <- exp_cfg$clinical
    if (!is.null(exp_cfg$features)) split_config$features <- exp_cfg$features

    t0_path <- if (!is.null(exp_cfg$input_file_t0)) exp_cfg$input_file_t0 else base_config$input_file_t0
    t1_path <- if (!is.null(exp_cfg$input_file_t1)) exp_cfg$input_file_t1 else base_config$input_file_t1
    if (!file.exists(t0_path)) stop(sprintf("[SPLIT] T0 input not found: %s", t0_path))

    split <- compute_cohort_split(readxl::read_excel(t0_path), split_config, split_cfg)

    split_dir <- file.path(run_root, exp_name, "00_cohort_split")
    if (!dir.exists(split_dir)) dir.create(split_dir, recursive = TRUE)

    paths_t0 <- write_split_excels(t0_path, "T0", split, split_dir, exp_name)
    if (run_longitudinal && !is.null(t1_path) && file.exists(t1_path))
      write_split_excels(t1_path, "T1", split, split_dir, exp_name)

    write.csv(split$report,
              file.path(split_dir, sprintf("split_manifest_%s.csv", exp_name)),
              row.names = FALSE)

    # Redirect every downstream pass to the discovery snapshot.
    exp_cfg$input_file    <- paths_t0$discovery
    exp_cfg$input_file_t0 <- paths_t0$discovery
    if (run_longitudinal)
      exp_cfg$input_file_t1 <- file.path(split_dir, sprintf("%s_discovery_T1.xlsx", exp_name))

    run_passes[[exp_name]] <- c("cohort_split", run_passes[[exp_name]])
  }

  message(sprintf("\n========================================================"))
  message(sprintf("STARTING EXPERIMENT BLOCK: %s", exp_name))
  message(sprintf("========================================================\n"))
  
  # ============================================================================
  # PASS 1: STANDARD CROSS-SECTIONAL PIPELINE (01, 02, 03[, 05])
  # Note: step 05 (network) runs here, before pass 2 step 04 (longitudinal).
  #       The two steps belong to independent analytical passes and different data modalities.
  # ============================================================================
  if (run_standard) {
    # Isolate configuration state for this pass
    config <- base_config
    config$project_name <- exp_name 
    config$run_mode <- "standard"   
    
    # Overwrite clinical logic if defined specifically in the experiment block
    if (!is.null(exp_cfg$clinical)) config$clinical <- exp_cfg$clinical

    # Per-experiment feature panel override (e.g. HNSCC has a different marker set).
    # soluble may legitimately be empty (FACS-only panels, e.g. the NSCLC v2 cohort).
    if (!is.null(exp_cfg$features)) {
      config$features <- exp_cfg$features
      if (length(unlist(config$features$facs)) == 0)
        stop(sprintf("[CONFIG] Experiment '%s': features.facs is empty", exp_name))
    }

    config$is_longitudinal <- FALSE

    if (!is.null(exp_cfg$input_file)) config$input_file <- exp_cfg$input_file
    
    config$output_root <- file.path(run_root, config$project_name)
    if (!dir.exists(config$output_root)) dir.create(config$output_root, recursive = TRUE)
    
    log_file_std <- file.path(config$output_root, paste0("log_standard_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
    
    message(sprintf("\n--- [PASS 1] STANDARD PIPELINE: %s ---", config$project_name))
    cat(sprintf("=== STANDARD PIPELINE STARTED: %s ===\n", Sys.time()), file = log_file_std)
    
    tryCatch({
      withCallingHandlers({
        message(">>> RUNNING PHASE 1: DATA PROCESSING <<<")
        source(here("src/01_data_processing.R"), echo = FALSE, local = FALSE)

        message("\n>>> RUNNING PHASE 2: VISUALIZATION <<<")
        source(here("src/02_visualization.R"), echo = FALSE, local = FALSE)

        message("\n>>> RUNNING PHASE 3: STATISTICAL ANALYSIS & REPORTING <<<")
        source(here("src/03_statistical_analysis.R"), echo = FALSE, local = FALSE)

        if (run_network) {
          message("\n>>> RUNNING PHASE 5: DIFFERENTIAL NETWORK ANALYSIS <<<")
          source(here("src/05_network_analysis.R"), echo = FALSE, local = FALSE)
        }

        message(sprintf("\n[SUCCESS] Standard Pipeline completed for: %s", exp_name))
      }, message = function(m) cat(conditionMessage(m), file = log_file_std, append = TRUE, sep = "\n"),
      warning = function(w) cat(paste0("WARNING: ", conditionMessage(w)), file = log_file_std, append = TRUE, sep = "\n"))
    }, error = function(e) {
      err_msg <- paste0("\n[FATAL ERROR] Standard Pipeline stopped: ", e$message)
      cat(err_msg, file = log_file_std, append = TRUE, sep = "\n")
      message(err_msg)
    })
  }
  
  # ============================================================================
  # PASS 2: LONGITUDINAL PIPELINE (01 Joint, 04 LMM)
  # ============================================================================
  if (run_longitudinal) {
    # Isolate configuration state for this pass
    config <- base_config
    config$project_name <- exp_name     
    config$run_mode <- "longitudinal"   
    
    if (!is.null(exp_cfg$clinical)) config$clinical <- exp_cfg$clinical

    # Forward-compat: allow longitudinal experiments with custom feature panels
    if (!is.null(exp_cfg$features)) config$features <- exp_cfg$features

    config$is_longitudinal <- TRUE

    # Critical Structural Rule: Disable Multivariate Outliers to prevent temporal censoring
    config$qc$remove_outliers <- FALSE
    
    if (!is.null(exp_cfg$input_file_t0)) config$input_file_t0 <- exp_cfg$input_file_t0
    if (!is.null(exp_cfg$input_file_t1)) config$input_file_t1 <- exp_cfg$input_file_t1
    
    config$output_root <- file.path(run_root, config$project_name)
    if (!dir.exists(config$output_root)) dir.create(config$output_root, recursive = TRUE)
    
    log_file_long <- file.path(config$output_root, paste0("log_longitudinal_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
    
    message(sprintf("\n--- [PASS 2] LONGITUDINAL PIPELINE: %s ---", config$project_name))
    cat(sprintf("=== LONGITUDINAL PIPELINE STARTED: %s ===\n", Sys.time()), file = log_file_long)
    
    tryCatch({
      withCallingHandlers({
        message(">>> RUNNING PHASE 1: JOINT DATA PROCESSING <<<")
        source(here("src/01_data_processing.R"), echo = FALSE, local = FALSE)
        
        message("\n>>> RUNNING PHASE 4: LONGITUDINAL LMM ANALYSIS <<<")
        source(here("src/04_longitudinal_lmm.R"), echo = FALSE, local = FALSE)
        
        message(sprintf("\n[SUCCESS] Longitudinal Pipeline completed for: %s", exp_name))
      }, message = function(m) cat(conditionMessage(m), file = log_file_long, append = TRUE, sep = "\n"),
      warning = function(w) cat(paste0("WARNING: ", conditionMessage(w)), file = log_file_long, append = TRUE, sep = "\n"))
    }, error = function(e) {
      err_msg <- paste0("\n[FATAL ERROR] Longitudinal Pipeline stopped: ", e$message)
      cat(err_msg, file = log_file_long, append = TRUE, sep = "\n")
      message(err_msg)
    })
  }
  # ============================================================================
  # PASS 3: MACHINE LEARNING PIPELINE (06)
  # Reads LOO-robust features from Step 04 JSON. Degrades gracefully when
  # no features survive the FDR + LOO gate — no per-experiment flag needed.
  # ============================================================================
  if (run_machine_learning) {
    config <- base_config
    config$project_name <- exp_name
    config$run_mode     <- "machine_learning"

    if (!is.null(exp_cfg$clinical)) config$clinical <- exp_cfg$clinical

    # Per-experiment overrides needed for multi-dataset support
    if (!is.null(exp_cfg$features))      config$features      <- exp_cfg$features
    if (!is.null(exp_cfg$input_file))    config$input_file    <- exp_cfg$input_file
    if (!is.null(exp_cfg$input_file_t0)) config$input_file_t0 <- exp_cfg$input_file_t0
    if (!is.null(exp_cfg$machine_learning)) {
      config$machine_learning <- modifyList(
        if (!is.null(config$machine_learning)) config$machine_learning else list(),
        exp_cfg$machine_learning
      )
    }

    config$output_root <- file.path(run_root, config$project_name)
    if (!dir.exists(config$output_root)) dir.create(config$output_root, recursive = TRUE)

    log_file_ml <- file.path(config$output_root, paste0("log_machine_learning_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))

    message(sprintf("\n--- [PASS 3] MACHINE LEARNING PIPELINE: %s ---", config$project_name))
    cat(sprintf("=== MACHINE LEARNING PIPELINE STARTED: %s ===\n", Sys.time()), file = log_file_ml)

    tryCatch({
      withCallingHandlers({
        message(">>> RUNNING PHASE 6: MACHINE LEARNING CLASSIFICATION <<<")
        source(here("src/06_machine_learning.R"), echo = FALSE, local = FALSE)

        message(sprintf("\n[SUCCESS] Machine Learning Pipeline completed for: %s", exp_name))
      }, message = function(m) cat(conditionMessage(m), file = log_file_ml, append = TRUE, sep = "\n"),
      warning = function(w) cat(paste0("WARNING: ", conditionMessage(w)), file = log_file_ml, append = TRUE, sep = "\n"))
    }, error = function(e) {
      err_msg <- paste0("\n[FATAL ERROR] Machine Learning Pipeline stopped: ", e$message)
      cat(err_msg, file = log_file_ml, append = TRUE, sep = "\n")
      message(err_msg)
    })
  }
}

# Publish provenance + the `latest` pointer once all experiments have run.
write_run_manifest(file.path(run_root, "run_manifest.yml"), run_id, base_config, run_passes)
update_latest_pointer(base_config$output_root, run_id)

options(crayon.enabled = TRUE)
message(sprintf("\n=== ALL EXPERIMENT BLOCKS COMPLETED (run: %s) ===", run_id))