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

# Per-experiment JSON locations — the ONLY place these file names are spelled out.
# `exp_root` is the directory that holds the 0X_* step folders for one experiment
# (i.e. <run_root>/<experiment>). Step directory names come from .STEP_DIRS in
# utils_io.R via step_output_path() (sourced before this file), so the folder
# names live in exactly one place across the codebase.
metrics_json_paths <- function(exp_root, experiment) {
  cfg <- list(output_root = exp_root)
  list(
    ml     = step_output_path(cfg, 6, sprintf("Machine_Metrics_ML_%s",  experiment), "json"),
    lmm    = step_output_path(cfg, 4, sprintf("Machine_Metrics_LMM_%s", experiment), "json"),
    splsda = step_output_path(cfg, 3, sprintf("Machine_Metrics_%s",     experiment), "json")
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

# Flatten an extract_run_metrics() result ($ml + $lmm) into one named list, so the
# headline fields can be addressed uniformly by assertions, comparisons and the
# publishability verdict. Single source of truth for that flattening (the golden
# assertion helper delegates here).
flatten_run_metrics <- function(m) {
  out <- if (!is.null(m$ml)) m$ml else list()
  if (!is.null(m$lmm)) out$n_sig_fdr <- m$lmm$n_sig_fdr
  out
}

# List the experiment subdirectories of a run root (those that actually contain at
# least one 0X_* step folder). Skips bookkeeping entries like 00_cohort_split.
list_run_experiments <- function(run_root) {
  if (!dir.exists(run_root)) stop(sprintf("[metrics] run root not found: %s", run_root))
  subs <- list.dirs(run_root, recursive = FALSE, full.names = FALSE)
  subs[vapply(subs, function(s)
    any(dir.exists(file.path(run_root, s, .STEP_DIRS))), logical(1))]
}

# ==============================================================================
# PUBLISHABILITY VERDICT
# Encodes the pre-specified go/no-go thresholds (CLAUDE.md / the /post-change
# skill) in ONE place so every consumer reads the same rule:
#   - permutation p < 0.05 on the PRIMARY method (SVM-RBF unless EN was promoted)
#   - nested-LOO AUC > 0.6
#   - (LMM-gated path only) FDR-significant feature count > 0
# ==============================================================================

#' @param m Output of \code{extract_run_metrics()} for one experiment.
#' @param thresholds Named list overriding the defaults (perm_p, nested_loo_auc,
#'   n_sig_fdr).
#' @return A list with \code{experiment}, \code{primary_method},
#'   \code{publishable} (logical, NA when ML did not run) and a \code{criteria}
#'   list (each: metric, value, op, threshold, pass).
publishability_verdict <- function(m,
                                   thresholds = list(perm_p = 0.05,
                                                     nested_loo_auc = 0.6,
                                                     n_sig_fdr = 0)) {
  if (is.null(m$ml))
    return(list(experiment = m$experiment, primary_method = NA_character_,
                publishable = NA, reason = "no Step 06 metrics (ML did not run)",
                criteria = list()))

  ml      <- m$ml
  primary <- .or_else(ml$primary_method, "SVM-RBF")
  is_svm  <- grepl("svm", primary, ignore.case = TRUE)
  perm_p  <- if (is_svm) ml$svm_perm_p else ml$en_perm_p

  mk <- function(metric, value, op, thr, pass)
    list(metric = metric, value = value, op = op, threshold = thr, pass = isTRUE(pass))

  criteria <- list(
    perm_p         = mk(sprintf("perm p (%s)", primary), perm_p, "<",
                        thresholds$perm_p, isTRUE(perm_p < thresholds$perm_p)),
    nested_loo_auc = mk("nested-LOO AUC", ml$nested_loo_auc, ">",
                        thresholds$nested_loo_auc,
                        isTRUE(ml$nested_loo_auc > thresholds$nested_loo_auc))
  )

  # FDR criterion is meaningful only for the LMM-gated (longitudinal) path; the
  # univariate cross-sectional path has no Step 04 and reports n_sig_fdr = NA.
  if (identical(.or_else(ml$gate_method, "lmm"), "lmm") && !is.null(m$lmm)) {
    criteria$n_sig_fdr <- mk("FDR-significant features", m$lmm$n_sig_fdr, ">",
                             thresholds$n_sig_fdr,
                             isTRUE(m$lmm$n_sig_fdr > thresholds$n_sig_fdr))
  }

  list(
    experiment     = m$experiment,
    primary_method = primary,
    publishable    = all(vapply(criteria, function(c) isTRUE(c$pass), logical(1))),
    criteria       = criteria
  )
}

# ==============================================================================
# RUN-TO-RUN COMPARISON
# Delta of headline metrics between two runs, for refactor regression checks and
# the /post-change workflow. Reads both runs through extract_run_metrics(), so the
# compared numbers are exactly those the golden tests assert on.
# ==============================================================================

#' @param run_base,run_new Run roots (dirs holding the per-experiment subdirs).
#' @param experiments Optional character vector; default = experiments present in
#'   BOTH runs.
#' @param tol Numeric tolerance for flagging a change (default 5e-4, matching the
#'   golden AUC tolerance).
#' @return A data.frame: experiment, metric, base, new, delta, changed.
compare_runs <- function(run_base, run_new, experiments = NULL, tol = 5e-4) {
  exps <- .or_else(experiments,
                   intersect(list_run_experiments(run_base),
                             list_run_experiments(run_new)))
  if (length(exps) == 0) return(data.frame())

  num_fields <- c("svm_auc", "svm_perm_p", "en_auc", "en_perm_p",
                  "nested_loo_auc", "n_samples", "n_features", "n_sig_fdr")
  chr_fields <- c("primary_method", "gate_method")

  fmt <- function(x) if (is.null(x) || length(x) == 0) "—"
                     else if (is.numeric(x)) format(x) else as.character(x)
  add <- function(rows, exp, metric, a, b, delta = NA_real_, changed = NA) {
    rbind(rows, data.frame(
      experiment = exp, metric = metric, base = fmt(a), new = fmt(b),
      delta = if (is.na(delta)) "" else sprintf("%+.4g", delta),
      changed = changed, stringsAsFactors = FALSE))
  }

  rows <- data.frame()
  for (exp in exps) {
    fa <- flatten_run_metrics(extract_run_metrics(file.path(run_base, exp), exp))
    fb <- flatten_run_metrics(extract_run_metrics(file.path(run_new, exp), exp))

    for (f in num_fields) {
      a <- fa[[f]]; b <- fb[[f]]
      if (is.null(a) && is.null(b)) next
      delta   <- if (is.numeric(a) && is.numeric(b)) b - a else NA_real_
      changed <- !isTRUE(all.equal(a, b, tolerance = tol))
      rows <- add(rows, exp, f, a, b, delta, changed)
    }
    for (f in chr_fields) {
      a <- fa[[f]]; b <- fb[[f]]
      if (is.null(a) && is.null(b)) next
      rows <- add(rows, exp, f, a, b, changed = !isTRUE(identical(a, b)))
    }
    ga <- .or_else(fa$gate, character(0)); gb <- .or_else(fb$gate, character(0))
    if (length(ga) > 0 || length(gb) > 0)
      rows <- add(rows, exp, "gate",
                  paste(sort(ga), collapse = ","), paste(sort(gb), collapse = ","),
                  changed = !setequal(ga, gb))
  }
  rows
}
