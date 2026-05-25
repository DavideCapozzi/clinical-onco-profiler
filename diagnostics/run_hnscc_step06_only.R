# diagnostics/run_hnscc_step06_only.R
# Re-runs only Step 06 (ML) for HNSCC_Response.
# Requires Step 01 RDS to already exist (run run_hnscc_only.R first).

suppressPackageStartupMessages({
  library(tidyverse); library(yaml); library(here); library(openxlsx)
})
options(crayon.enabled = FALSE)

base_config <- yaml::read_yaml(here("config/global_params.yml"))
list.files(here("R"), pattern = "\\.R$", full.names = TRUE) |> purrr::walk(source)

exp_name <- "HNSCC_Response"
exp_cfg  <- base_config$experiments[[exp_name]]
if (is.null(exp_cfg)) stop(sprintf("[FATAL] Experiment '%s' not found.", exp_name))

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
config$output_root <- file.path(base_config$output_root, exp_name)

message(sprintf("\n--- HNSCC Step 06 only (updated CPS parser + stratified analysis) ---"))

tryCatch({
  source(here("src/06_machine_learning.R"), echo = FALSE, local = FALSE)
  message("\n[SUCCESS] Step 06 completed.")
}, error = function(e) {
  message(sprintf("\n[FATAL] Step 06 failed: %s", e$message))
})
