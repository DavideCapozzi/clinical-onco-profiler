# utils_metrics.R
# ==============================================================================
# Single source of truth for reading HEADLINE metrics out of a completed run.
#
# JSON paths and field locations are defined here ONCE so that every consumer
# reads them the same way:
#   - the golden regression tests (tests/golden/*)        [Fase 1]
#   - compare_runs() / publishability_verdict()           [Fase 3]
#   - the /post-change skill (via a thin R runner)        [Fase 3]
#
# Fase 1 scope = extract_run_metrics() only (minimal: just the headline numbers
# the golden tests assert on). Fase 3 adds compare_runs() + publishability_verdict().
# ==============================================================================

# Internal: first non-empty (non-NULL, non-zero-length) of a, else b.
# Deliberately NOT exported as `%||%` to avoid clobbering rlang's operator (which
# only coalesces on NULL) when all of R/ is sourced into a pipeline run — here we
# also need to fall back on zero-length values (e.g. an empty gate list()).
.or_else <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# Per-experiment JSON locations — the ONLY place these paths are spelled out.
# `exp_root` is the directory that holds the 0X_* step folders for one experiment
# (i.e. <run_root>/<experiment>).
metrics_json_paths <- function(exp_root, experiment) {
  list(
    ml     = file.path(exp_root, "06_machine_learning",     sprintf("Machine_Metrics_ML_%s.json",  experiment)),
    lmm    = file.path(exp_root, "04_longitudinal_analysis", sprintf("Machine_Metrics_LMM_%s.json", experiment)),
    splsda = file.path(exp_root, "03_statistical_analysis",  sprintf("Machine_Metrics_%s.json",     experiment))
  )
}

# Extract the headline metrics for ONE experiment from its run directory.
# Returns a normalized list; nodes are NULL when the corresponding step did not
# run (e.g. no Step 04 for a cross-sectional/univariate experiment).
extract_run_metrics <- function(exp_root, experiment = basename(exp_root)) {
  if (!dir.exists(exp_root))
    stop(sprintf("[metrics] experiment dir not found: %s", exp_root))

  paths <- metrics_json_paths(exp_root, experiment)
  out <- list(experiment = experiment, exp_root = exp_root)

  # --- Step 06: machine learning -------------------------------------------
  if (file.exists(paths$ml)) {
    ml <- jsonlite::fromJSON(paths$ml, simplifyVector = FALSE)

    # Gate markers live in nested_loocv_validation$gate_stability (lmm path) or
    # top-level gate_stability (univariate path, already fully nested).
    gate_src <- .or_else(ml$nested_loocv_validation$gate_stability, ml$gate_stability)
    gate <- if (!is.null(gate_src))
      vapply(gate_src, function(x) .or_else(x$Marker, NA_character_), character(1)) else character(0)

    # nested-LOO AUC: explicit node for the lmm path; for the univariate path the
    # primary SVM AUC is itself fully nested by construction.
    nested_loo <- .or_else(ml$nested_loocv_validation$auc, ml$svm_rbf$metrics$auc)

    out$ml <- list(
      gate_method    = .or_else(ml$gate_method, "lmm"),
      n_samples      = ml$n_samples,
      n_features     = .or_else(ml$n_features_total, ml$n_features_post_filter),
      primary_method = ml$primary_method,
      svm_auc        = ml$svm_rbf$metrics$auc,
      svm_perm_p     = ml$permutation_test$svm_rbf$p_value,
      en_auc         = ml$elastic_net$metrics$auc,
      en_perm_p      = ml$permutation_test$elastic_net$p_value,
      nested_loo_auc = nested_loo,
      gate           = gate
    )
  }

  # --- Step 04: longitudinal LMM (FDR count) -------------------------------
  if (file.exists(paths$lmm)) {
    lmm <- jsonlite::fromJSON(paths$lmm, simplifyVector = FALSE)
    sf  <- lmm$significant_features_fdr
    # `significant_features_fdr` may be serialized either as a scalar count or as
    # a list of marker records — normalize both to an integer count.
    n_sig <- if (is.null(sf)) NA_integer_
             else if (length(sf) == 1 && is.numeric(sf[[1]])) as.integer(sf[[1]])
             else length(sf)
    out$lmm <- list(n_sig_fdr = n_sig, n_paired = lmm$n_paired)
  }

  out
}
