# tests/golden/assert_metrics.R
# ==============================================================================
# Shared assertion helper for the golden regression tests.
# Compares a metrics list (from extract_run_metrics()) against a frozen
# `expected` block (from expected_metrics.yml) and collects every mismatch.
# Returns the vector of failure strings (empty == pass) so callers can decide
# how to exit. Numeric fields use `tol`; perm-p is exact (seed-fixed pipeline);
# the gate is compared as an unordered set.
# ==============================================================================

# Compare one extracted experiment against its expected block.
# `got` is the $ml (+ optional $lmm) portion flattened by the caller.
assert_block <- function(label, got, expected, tol_auc = 5e-4) {
  fails <- character(0)
  chk_num <- function(field, exact = FALSE) {
    if (is.null(expected[[field]])) return(invisible())
    g <- got[[field]]
    if (is.null(g) || length(g) == 0 || is.na(g)) {
      fails[[length(fails) + 1]] <<- sprintf("%s/%s: missing (expected %s)", label, field, expected[[field]])
      return(invisible())
    }
    ok <- if (exact) isTRUE(g == expected[[field]]) else isTRUE(abs(g - expected[[field]]) <= tol_auc)
    if (!ok)
      fails[[length(fails) + 1]] <<- sprintf("%s/%s: got %s, expected %s%s",
        label, field, format(g), format(expected[[field]]),
        if (exact) "" else sprintf(" (tol %s)", tol_auc))
  }
  chk_chr <- function(field) {
    if (is.null(expected[[field]])) return(invisible())
    g <- got[[field]]
    if (!isTRUE(g == expected[[field]]))
      fails[[length(fails) + 1]] <<- sprintf("%s/%s: got '%s', expected '%s'", label, field, g, expected[[field]])
  }

  chk_num("svm_auc");        chk_num("en_auc");        chk_num("nested_loo_auc")
  chk_num("svm_perm_p", exact = TRUE)
  chk_num("n_samples", exact = TRUE); chk_num("n_features", exact = TRUE)
  chk_num("n_sig_fdr", exact = TRUE)
  chk_chr("primary_method"); chk_chr("gate_method")

  # Use [["gate"]] (not $gate): `$` does partial matching and would silently
  # resolve to `gate_method` for experiments that omit an explicit gate (e.g. the
  # univariate HNSCC block). `<-` (not `<<-`): `fails` is local to this function.
  if (!is.null(expected[["gate"]])) {
    g <- got$gate
    if (!setequal(g, expected[["gate"]]))
      fails[[length(fails) + 1]] <- sprintf("%s/gate: got {%s}, expected {%s}",
        label, paste(sort(g), collapse = ","), paste(sort(unlist(expected[["gate"]])), collapse = ","))
  }
  fails
}

# Flatten an extract_run_metrics() result for assert_block(). Delegates to the
# single source of truth in R/utils_metrics.R (sourced before this file in both
# golden runners), so the flattening rule lives in exactly one place.
flatten_metrics <- function(m) flatten_run_metrics(m)

# Top-level: given extracted metrics + expected block, print and return result.
report_assertions <- function(label, m, expected, tol_auc = 5e-4) {
  fails <- assert_block(label, flatten_metrics(m), expected, tol_auc)
  if (length(fails) == 0) {
    message(sprintf("  [PASS] %s — all headline metrics match frozen snapshot", label))
  } else {
    message(sprintf("  [FAIL] %s — %d mismatch(es):", label, length(fails)))
    for (f in fails) message(sprintf("    - %s", f))
  }
  fails
}
