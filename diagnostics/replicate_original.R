# replicate_original.R
# ==============================================================================
# Definitive replication test: run the CURRENT code on the ORIGINAL anonimi
# cohort (77 raw) with NO validation split — i.e. the exact data conditions of
# oldresults (0.716). Isolates "cohort + no-split" as the cause of the drop.
# Steps: 01 standard -> 01 long + 04 (LMM gate) -> 06 ML. Bootstrap CI off.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(yaml); library(here)
  library(openxlsx); library(jsonlite); library(readxl)
})
options(crayon.enabled = FALSE)
list.files(here("R"), pattern = "\\.R$", full.names = TRUE) %>% purrr::walk(source)

base_config <- yaml::read_yaml(here("config/global_params.yml"))
base_config$lmm_bootstrap$enabled <- FALSE
EXP_NAME <- "BestResponse_2v3_4"
exp_cfg  <- base_config$experiments[[EXP_NAME]]

# Force ORIGINAL files, NO split (replicate oldresults conditions exactly)
exp_cfg$validation_split <- NULL
exp_cfg$input_file    <- "data/Dati_NSCLC_standardizzati_anonimi_T0.xlsx"
exp_cfg$input_file_t0 <- "data/Dati_NSCLC_standardizzati_anonimi_T0.xlsx"
exp_cfg$input_file_t1 <- "data/Dati_NSCLC_standardizzati_anonimi_T1.xlsx"

run_root <- file.path(base_config$output_root,
                      paste0("replicate_original_", format(Sys.time(), "%Y%m%d_%H%M%S")))
exp_out  <- file.path(run_root, EXP_NAME)
dir.create(exp_out, recursive = TRUE, showWarnings = FALSE)
message(sprintf("[REPLICATE] root: %s", run_root))

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
  config$machine_learning <- modifyList(config$machine_learning %||% list(), exp_cfg$machine_learning)
config$output_root <- exp_out
source(here("src/06_machine_learning.R"), echo = FALSE, local = FALSE)

# Harvest
ml <- jsonlite::fromJSON(file.path(exp_out, "06_machine_learning",
        sprintf("Machine_Metrics_ML_%s.json", EXP_NAME)), simplifyVector = FALSE)
gate <- paste(sapply(ml$nested_loocv_validation$gate_stability, function(x) x$Marker), collapse = "+")
message("\n================ REPLICATION RESULT (original 77, no split) ================")
message(sprintf("n_samples   : %s", ml$n_samples))
message(sprintf("n_features  : %s", ml$n_features_total))
message(sprintf("primary     : %s", ml$primary_method))
message(sprintf("SVM AUC     : %.4f  (perm p = %s)", ml$svm_rbf$metrics$auc, ml$permutation_test$svm_rbf$p_value))
message(sprintf("EN  AUC     : %.4f  (perm p = %s)", ml$elastic_net$metrics$auc, ml$permutation_test$elastic_net$p_value))
message(sprintf("nested LOO  : %.4f", ml$nested_loocv_validation$auc))
message(sprintf("gate        : %s", gate))
message("\nEXPECTED (oldresults): SVM AUC 0.7158, perm p 0.001, nested LOO 0.7542")
