# main.R
# ==============================================================================
# Pipeline orchestrator: clinical-onco-profiler
# Runs, for each experiment of the config: pass 1 (steps 01-03, cross-sectional),
# pass 2 (01 + 04, longitudinal LMM), pass 3 (06, machine learning / added value),
# pass 4 (07, post hoc), then collects the publication set and writes the run
# manifest. Usage (from the repository root):
#   RUN_NAME=<run> CONFIG_PATH=config/<config>.yml Rscript main.R
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

# Configuration: CONFIG_PATH overrides the default (config/0final.yml).
.config_override <- Sys.getenv("CONFIG_PATH")
base_config_path <- if (nzchar(.config_override)) .config_override else here("config/0final.yml")
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
# experiments. Unset -> run all (default).
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
check_figure_export(base_config)   # figures.export needs poppler: fail now, not after the run

# NLR switch: whether the clinical model includes NLR (config `include_nlr`; the
# INCLUDE_NLR environment variable overrides it). A timestamped run dir (no RUN_NAME)
# is labelled with_nlr / no_nlr. Applied in step 06; reaches every pass via base_config.
include_nlr <- resolve_include_nlr(base_config)   # R/utils_io.R — shared with prepare_pass_config()
base_config$include_nlr <- include_nlr
.variant <- if (include_nlr) "with_nlr" else "no_nlr"
base_config$run_label <- if (!is.null(base_config$run_label) && nzchar(base_config$run_label))
                           paste(base_config$run_label, .variant, sep = "_") else .variant
message(sprintf("[System] include_nlr = %s (clinical model %s NLR)",
                include_nlr, if (include_nlr) "includes" else "excludes"))

# Isolate this invocation in one run directory: <output_root>/<RUN_NAME> when
# RUN_NAME is set (refuses an existing dir unless FORCE=1), else a timestamped id.
# Every pass routes its output_root under run_root, so all cross-step artifacts of
# a single main.R call stay self-contained and reproducible.
.run     <- resolve_run_dir(base_config$output_root, base_config$run_label)
run_id   <- .run$run_id
run_root <- .run$run_root
if (!dir.exists(run_root)) dir.create(run_root, recursive = TRUE)
message(sprintf("[System] Run ID: %s  ->  %s", run_id, run_root))

run_passes <- list()  # experiment -> passes executed (for the run manifest)
failed_passes <- character(0)  # "<experiment>: <pass>" for every pass that hit a FATAL ERROR

# 2. Pipeline Execution Loop
# ------------------------------------------------------------------------------
for (exp_name in names(experiments_list)) {

  exp_cfg <- experiments_list[[exp_name]]

  # Skip experiments disabled in the config (absence of the key means enabled).
  if (!is.null(exp_cfg$enabled) && !isTRUE(exp_cfg$enabled)) {
    message(sprintf("[System] Skipping disabled experiment: %s", exp_name))
    next
  }

  # Inherit boolean flags with safety fallbacks
  run_standard        <- if (!is.null(exp_cfg$run_standard))         as.logical(exp_cfg$run_standard)         else TRUE
  run_longitudinal    <- if (!is.null(exp_cfg$is_longitudinal))      as.logical(exp_cfg$is_longitudinal)      else FALSE
  # ML is a global flag (base_config), not per-experiment — Step 06 self-gates on LOO-robust feature availability
  run_machine_learning <- if (!is.null(base_config$run_machine_learning)) as.logical(base_config$run_machine_learning) else FALSE
  # Post hoc (step 07) needs steps 04 and 06; it self-skips (logged) when their artifacts are missing.
  run_post_hoc <- run_longitudinal && run_machine_learning && !isFALSE(base_config$run_post_hoc)

  run_passes[[exp_name]] <- c(
    if (run_standard)         "standard",
    if (run_longitudinal)     "longitudinal",
    if (run_machine_learning) "machine_learning",
    if (run_post_hoc)         "post_hoc"
  )

  message(sprintf("\n========================================================"))
  message(sprintf("STARTING EXPERIMENT BLOCK: %s", exp_name))
  message(sprintf("========================================================\n"))
  
  # ============================================================================
  # PASS 1: STANDARD CROSS-SECTIONAL PIPELINE (01, 02, 03)
  # ============================================================================
  if (run_standard) {
    # Isolate configuration state for this pass (see prepare_pass_config in utils_io.R)
    config <- prepare_pass_config(base_config, exp_cfg, exp_name, "standard", run_root)

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

        message(sprintf("\n[SUCCESS] Standard Pipeline completed for: %s", exp_name))
      }, message = function(m) cat(conditionMessage(m), file = log_file_std, append = TRUE, sep = "\n"),
      warning = function(w) cat(paste0("WARNING: ", conditionMessage(w)), file = log_file_std, append = TRUE, sep = "\n"))
    }, error = function(e) {
      err_msg <- paste0("\n[FATAL ERROR] Standard Pipeline stopped: ", e$message)
      failed_passes <<- c(failed_passes, paste0(exp_name, ": standard"))
      cat(err_msg, file = log_file_std, append = TRUE, sep = "\n")
      message(err_msg)
    })
  }
  
  # ============================================================================
  # PASS 2: LONGITUDINAL PIPELINE (01 Joint, 04 LMM)
  # ============================================================================
  if (run_longitudinal) {
    # Isolate configuration state for this pass (see prepare_pass_config in utils_io.R)
    config <- prepare_pass_config(base_config, exp_cfg, exp_name, "longitudinal", run_root)

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
      failed_passes <<- c(failed_passes, paste0(exp_name, ": longitudinal"))
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
    # Isolate configuration state for this pass (see prepare_pass_config in utils_io.R)
    config <- prepare_pass_config(base_config, exp_cfg, exp_name, "machine_learning", run_root)

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
      failed_passes <<- c(failed_passes, paste0(exp_name, ": machine_learning"))
      cat(err_msg, file = log_file_ml, append = TRUE, sep = "\n")
      message(err_msg)
    })
  }

  # ============================================================================
  # PASS 4: POST HOC (07) — grouping, duplicated rows, therapy, KM on the
  # finished step-04 / step-06 artifacts. Self-skips when they are missing.
  # ============================================================================
  if (run_post_hoc) {
    config <- prepare_pass_config(base_config, exp_cfg, exp_name, "post_hoc", run_root)

    log_file_ph <- file.path(config$output_root, paste0("log_post_hoc_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))

    message(sprintf("\n--- [PASS 4] POST HOC: %s ---", config$project_name))
    cat(sprintf("=== POST HOC STARTED: %s ===\n", Sys.time()), file = log_file_ph)

    tryCatch({
      withCallingHandlers({
        message(">>> RUNNING PHASE 7: POST HOC ANALYSES <<<")
        source(here("src/07_post_hoc.R"), echo = FALSE, local = FALSE)

        message(sprintf("\n[SUCCESS] Post hoc completed for: %s", exp_name))
      }, message = function(m) cat(conditionMessage(m), file = log_file_ph, append = TRUE, sep = "\n"),
      warning = function(w) cat(paste0("WARNING: ", conditionMessage(w)), file = log_file_ph, append = TRUE, sep = "\n"))
    }, error = function(e) {
      err_msg <- paste0("\n[FATAL ERROR] Post Hoc Pipeline stopped: ", e$message)
      failed_passes <<- c(failed_passes, paste0(exp_name, ": post_hoc"))
      cat(err_msg, file = log_file_ph, append = TRUE, sep = "\n")
      message(err_msg)
    })
  }
}

# Publication set: copy the mapped figures / tables into <run>/publication/ (+ INDEX.md).
# A mapped main or supplementary file that is missing is FATAL (logged, like a failed pass).
if (!is.null(base_config$publication)) {
  log_file_pub <- file.path(run_root, paste0("log_publication_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
  cat(sprintf("=== PUBLICATION COLLECTOR STARTED: %s ===\n", Sys.time()), file = log_file_pub)
  tryCatch(
    withCallingHandlers(collect_publication(run_root, base_config, experiments = names(run_passes)),
      message = function(m) cat(conditionMessage(m), file = log_file_pub, append = TRUE, sep = "\n"),
      warning = function(w) cat(paste0("WARNING: ", conditionMessage(w)), file = log_file_pub, append = TRUE, sep = "\n")),
    error = function(e) {
      err_msg <- paste0("\n[FATAL ERROR] Publication collector stopped: ", e$message)
      failed_passes <<- c(failed_passes, "publication")
      cat(err_msg, file = log_file_pub, append = TRUE, sep = "\n")
      message(err_msg)
    })
}

# Publish provenance once all experiments have run.
write_run_manifest(file.path(run_root, "run_manifest.yml"), run_id, base_config, run_passes,
                   config_path = base_config_path, experiments_cfg = experiments_list,
                   failed_passes = failed_passes)

options(crayon.enabled = TRUE)

# A failed pass is logged and the remaining passes still run, but the process must
# not exit 0: stop() at top level makes Rscript exit with status 1.
if (length(failed_passes) > 0)
  stop(sprintf("[FATAL ERROR] %d failed pass(es) in run %s: %s",
               length(failed_passes), run_id, paste(failed_passes, collapse = "; ")), call. = FALSE)
message(sprintf("\n=== ALL EXPERIMENT BLOCKS COMPLETED (run: %s) ===", run_id))