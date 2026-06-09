# R/modules_ml.R
# ==============================================================================
# MACHINE LEARNING MODULE — aggregator
# Description: Nested-LOOCV classification on LOO-robust biomarkers from LMM.
#              Primary method: Elastic Net Logistic Regression (glmnet).
#              Secondary method: SVM with RBF Kernel (e1071).
# Dependencies: dplyr, ggplot2, tidyr, glmnet, e1071, pROC
#
# The implementation lives in four themed files, split by responsibility:
#   modules_ml_gate.R     — feature-gate selection + gate signal decomposition
#   modules_ml_cv.R       — matrix building, nested-LOOCV (glmnet/SVM), metrics, permutation
#   modules_ml_utility.R  — clinical benchmark, stratified/combined models, clinical-utility layer
#   modules_ml_plots.R    — all ggplot/PDF rendering
#
# This file exists so `source("R/modules_ml.R")` (src/06_machine_learning.R) keeps
# loading the whole module. When the pipeline loader globs R/*.R the four files are
# also picked up directly; re-sourcing them here is idempotent (function defs only).
# ==============================================================================

local({
  ml_files <- c("modules_ml_gate.R", "modules_ml_cv.R",
                "modules_ml_utility.R", "modules_ml_plots.R")
  dir <- if (requireNamespace("here", quietly = TRUE)) here::here("R") else "R"
  for (f in ml_files) source(file.path(dir, f))
})
