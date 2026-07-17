# R/modules_ml_utility.R
# ==============================================================================
# MACHINE LEARNING MODULE — UTILITY
# Split from the former monolithic R/modules_ml.R (sourced via that aggregator).
# Dependencies: dplyr, ggplot2, tidyr, glmnet, e1071, pROC
# ==============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)

#' @title Run Clinical Benchmark Comparison
#' @description
#' Loads a clinical biomarker column from the raw input Excel, merges by Patient_ID,
#' computes its univariate AUC on available cases, and performs a DeLong test against
#' the primary model's out-of-fold predicted probabilities. Returns NULL gracefully when
#' the column is absent, has too many NAs, or the input file cannot be read.
#'
#' @param DATA Named list. Processed data object from Step 01 (for Patient_ID and Group).
#' @param primary_probs Numeric vector. Out-of-fold predicted probabilities from the primary model.
#' @param y_primary Factor. True class labels (same order as primary_probs).
#' @param positive_label String. The positive class label.
#' @param input_file String. Path to the raw input Excel (config$input_file_t0).
#' @param benchmark_col String. Column name for the clinical biomarker.
#' @param benchmark_label String. Human-readable label for reporting.
#' @param min_n Integer. Minimum valid cases required to proceed (default 10).
#' @return A named list: label, auc, auc_ci, n_valid, n_na, delong_p, predicted_probs,
#'   y_true, positive_label. Returns NULL on failure.
#' @export
run_clinical_benchmark <- function(DATA, primary_probs, y_primary, positive_label,
                                   input_file, benchmark_col, benchmark_label = NULL,
                                   min_n = 10L) {
  if (!requireNamespace("pROC",    quietly = TRUE)) return(NULL)
  if (!requireNamespace("readxl",  quietly = TRUE)) return(NULL)
  if (!requireNamespace("dplyr",   quietly = TRUE)) return(NULL)

  if (is.null(benchmark_label)) benchmark_label <- benchmark_col

  if (!file.exists(input_file)) {
    message(sprintf("   [ML Benchmark] Input file not found: %s", input_file))
    return(NULL)
  }

  df_raw <- tryCatch(
    readxl::read_excel(input_file, sheet = 1),
    error = function(e) { message(sprintf("   [ML Benchmark] Cannot read %s: %s", input_file, e$message)); NULL }
  )
  if (is.null(df_raw)) return(NULL)

  if (!benchmark_col %in% names(df_raw)) {
    message(sprintf("   [ML Benchmark] Column '%s' not found in %s", benchmark_col, basename(input_file)))
    return(NULL)
  }

  # Merge: raw file → processed metadata by Patient_ID
  df_bench <- df_raw %>%
    dplyr::select(Patient_ID, dplyr::all_of(benchmark_col)) %>%
    dplyr::mutate(Patient_ID = as.character(Patient_ID))

  df_meta <- DATA$metadata %>%
    dplyr::mutate(Patient_ID = as.character(Patient_ID)) %>%
    dplyr::left_join(df_bench, by = "Patient_ID")

  bench_vals <- as.numeric(df_meta[[benchmark_col]])
  n_na    <- sum(is.na(bench_vals))
  n_valid <- sum(!is.na(bench_vals))

  message(sprintf("   [ML Benchmark] '%s': %d valid / %d total (NA=%d)",
                  benchmark_label, n_valid, nrow(df_meta), n_na))

  if (n_valid < min_n) {
    message(sprintf("   [ML Benchmark] Insufficient valid cases (%d < %d). Skipping.", n_valid, min_n))
    return(NULL)
  }

  valid_mask <- !is.na(bench_vals)
  y_bench    <- y_primary[valid_mask]
  y_bin      <- as.integer(y_bench == positive_label)
  x_bench    <- bench_vals[valid_mask]

  # Univariate AUC (direction auto-selected for maximum AUC)
  roc_fwd <- tryCatch(pROC::roc(y_bin, x_bench, direction = "<", quiet = TRUE), error = function(e) NULL)
  roc_rev <- tryCatch(pROC::roc(y_bin, x_bench, direction = ">", quiet = TRUE), error = function(e) NULL)
  auc_fwd <- if (!is.null(roc_fwd)) as.numeric(pROC::auc(roc_fwd)) else 0
  auc_rev <- if (!is.null(roc_rev)) as.numeric(pROC::auc(roc_rev)) else 0
  roc_best  <- if (auc_fwd >= auc_rev) roc_fwd else roc_rev
  bench_auc <- max(auc_fwd, auc_rev)

  ci_bench <- tryCatch(
    as.numeric(pROC::ci.auc(roc_best, method = "delong")),
    error = function(e) c(NA_real_, bench_auc, NA_real_)
  )

  # DeLong test: primary model (subset) vs benchmark
  probs_subset <- primary_probs[valid_mask]
  y_bin_sub    <- as.integer(y_primary[valid_mask] == positive_label)
  roc_primary  <- tryCatch(
    pROC::roc(y_bin_sub, probs_subset, direction = "<", quiet = TRUE),
    error = function(e) NULL
  )
  primary_auc_subset <- if (!is.null(roc_primary)) as.numeric(pROC::auc(roc_primary)) else NA_real_

  delong_result <- if (!is.null(roc_primary) && !is.null(roc_best)) {
    tryCatch(pROC::roc.test(roc_primary, roc_best, method = "delong"), error = function(e) NULL)
  } else NULL

  delong_p <- if (!is.null(delong_result)) round(delong_result$p.value, 5) else NA_real_

  message(sprintf(
    "   [ML Benchmark] %s: AUC=%.3f [%.3f-%.3f] | Primary(n=%d)=%.3f | DeLong p=%.4f",
    benchmark_label, bench_auc, ci_bench[1], ci_bench[3],
    n_valid, primary_auc_subset, if (is.na(delong_p)) 0 else delong_p
  ))

  list(
    label                 = benchmark_label,
    column                = benchmark_col,
    auc                   = bench_auc,
    auc_ci                = ci_bench,
    n_valid               = n_valid,
    n_na                  = n_na,
    primary_auc_on_subset = primary_auc_subset,
    delong_p              = delong_p,
    direction             = if (auc_fwd >= auc_rev) "<" else ">",
    predicted_probs       = x_bench,
    y_true                = y_bench,
    positive_label        = positive_label
  )
}

#' @title Benchmark-Stratified Subgroup Analysis
#' @description
#' Stratifies the cohort by clinical-standard categorical bins of the benchmark
#' biomarker (e.g. PD-L1 TPS at <1%, 1-49%, >=50%) and quantifies:
#'   - Frequency table of bins vs outcome (Fisher exact + Cochran-Armitage trend)
#'   - Binary cut performance ("benchmark high" vs "benchmark low") as classifier
#'   - Primary-model AUC restricted to the "benchmark low" subgroup, which is
#'     the clinically actionable population where the standard biomarker is
#'     insufficient for treatment decision
#'
#' @param DATA Step 01 standard DATA object (with metadata$Patient_ID).
#' @param X_main Numeric matrix of main-effect features (n_patients x n_main).
#' @param y Factor outcome aligned with X_main.
#' @param positive_label Positive class label.
#' @param input_file Path to raw clinical Excel (T0) with benchmark column.
#' @param benchmark_col Column name in raw Excel.
#' @param benchmark_label Optional human-readable label.
#' @param bin_breaks Numeric vector of bin breakpoints (default c(0.99, 49.99)
#'   for PD-L1 TPS clinical strata <1%, 1-49%, >=50%).
#' @param bin_labels Character vector of bin labels (default
#'   c("neg(<1%)", "low(1-49%)", "high(>=50%)")).
#' @param high_threshold Numeric. Cutoff defining "benchmark high" for the
#'   binary clinical-cut classifier (default 50).
#' @param C_grid,gamma_grid,inner_folds,seed SVM tuning hyperparameters.
#' @param min_subgroup_n Minimum patients per subgroup to attempt SVM fit
#'   (default 20). Smaller subgroups are reported with raw counts only.
#' @return A list with crosstab, fisher_p, ca_trend_p, binary_cut metrics, and
#'   per-subgroup nested-LOOCV SVM AUC (when subgroup_n >= min_subgroup_n).
#'   Returns NULL if benchmark column is unavailable.
#' @export
run_pdl1_stratified <- function(DATA, X_main, y, positive_label,
                                input_file, benchmark_col, benchmark_label = NULL,
                                bin_breaks     = c(0.99, 49.99),
                                bin_labels     = c("neg(<1%)", "low(1-49%)", "high(>=50%)"),
                                high_threshold = 50,
                                C_grid         = c(0.01, 0.1, 1, 10, 100),
                                gamma_grid     = c(0.01, 0.1, 1, 10),
                                inner_folds    = 5L,
                                seed           = 2026,
                                min_subgroup_n = 15L) {
  if (!requireNamespace("readxl", quietly = TRUE)) return(NULL)
  if (!requireNamespace("dplyr",  quietly = TRUE)) return(NULL)
  if (!requireNamespace("pROC",   quietly = TRUE)) return(NULL)

  if (is.null(benchmark_label)) benchmark_label <- benchmark_col
  if (!file.exists(input_file)) {
    message(sprintf("   [ML Stratified] Input file not found: %s", input_file))
    return(NULL)
  }

  df_raw <- tryCatch(readxl::read_excel(input_file, sheet = 1),
                     error = function(e) NULL)
  if (is.null(df_raw) || !benchmark_col %in% names(df_raw)) {
    message(sprintf("   [ML Stratified] Column '%s' not available — skipping.", benchmark_col))
    return(NULL)
  }

  df_bench <- df_raw %>%
    dplyr::select(Patient_ID, dplyr::all_of(benchmark_col)) %>%
    dplyr::mutate(Patient_ID = as.character(Patient_ID))

  df_meta <- DATA$metadata %>%
    dplyr::mutate(Patient_ID = as.character(Patient_ID)) %>%
    dplyr::left_join(df_bench, by = "Patient_ID")

  bench_vals <- as.numeric(df_meta[[benchmark_col]])
  n_valid <- sum(!is.na(bench_vals))
  if (n_valid < min_subgroup_n) {
    message(sprintf("   [ML Stratified] Insufficient valid cases (%d). Skipping.", n_valid))
    return(NULL)
  }

  # ---- Categorical bins (3 levels) ----
  bins <- cut(bench_vals,
              breaks = c(-Inf, bin_breaks, Inf),
              labels = bin_labels, include.lowest = TRUE)

  y_chr  <- as.character(y)
  pos    <- positive_label
  neg    <- setdiff(unique(y_chr), pos)
  tab3   <- table(bins, factor(y_chr, levels = c(neg, pos)), useNA = "no")
  fish3  <- tryCatch(fisher.test(tab3)$p.value, error = function(e) NA_real_)

  ca_p <- NA_real_
  if (nrow(tab3) >= 2) {
    rs <- rowSums(tab3)
    if (all(rs > 0)) {
      ca_x <- tab3[, pos]
      ca_n <- rs
      ca_p <- tryCatch(prop.trend.test(x = ca_x, n = ca_n)$p.value,
                       error = function(e) NA_real_)
    }
  }

  rr_per_bin <- as.numeric(prop.table(tab3, margin = 1)[, pos]) * 100
  bin_summary <- data.frame(
    Bin            = rownames(tab3),
    N              = as.numeric(rowSums(tab3)),
    N_Responder    = as.numeric(tab3[, pos]),
    Response_Rate  = round(rr_per_bin, 1),
    stringsAsFactors = FALSE
  )

  # ---- Binary cut (clinical threshold) as a classifier ----
  high_mask <- !is.na(bench_vals) & bench_vals >= high_threshold
  low_mask  <- !is.na(bench_vals) & bench_vals <  high_threshold
  y_b <- y_chr[high_mask | low_mask]
  pred_high <- ifelse(bench_vals[high_mask | low_mask] >= high_threshold, pos, neg)
  sens_b <- mean(pred_high[y_b == pos] == pos)
  spec_b <- mean(pred_high[y_b == neg] == neg)
  bal_b  <- (sens_b + spec_b) / 2
  tab2   <- table(factor(pred_high, levels = c(neg, pos)),
                  factor(y_b, levels = c(neg, pos)))
  fish2  <- tryCatch(fisher.test(tab2), error = function(e) NULL)
  binary_cut <- list(
    threshold         = high_threshold,
    n                 = length(y_b),
    n_high            = sum(high_mask),
    n_low             = sum(low_mask),
    balanced_accuracy = round(bal_b, 4),
    sensitivity       = round(sens_b, 4),
    specificity       = round(spec_b, 4),
    fisher_p          = if (!is.null(fish2)) round(fish2$p.value, 4) else NA_real_,
    odds_ratio        = if (!is.null(fish2)) round(as.numeric(fish2$estimate), 4) else NA_real_
  )

  # ---- Per-subgroup nested-LOOCV SVM on the LMM-gate features ----
  run_subgroup_svm <- function(mask, label) {
    n_s <- sum(mask)
    if (n_s < min_subgroup_n) {
      return(list(n = n_s,
                  n_responder = sum(y[mask] == pos),
                  auc = NA_real_, auc_ci = c(NA_real_, NA_real_, NA_real_),
                  predicted_probs = NULL, y_true = NULL, benchmark_vals = NULL,
                  note = sprintf("n=%d below min_subgroup_n=%d — not fitted.",
                                 n_s, min_subgroup_n)))
    }
    y_sub <- y[mask]
    if (length(unique(y_sub)) < 2 || min(table(y_sub)) < 3) {
      return(list(n = n_s,
                  n_responder = sum(y_sub == pos),
                  auc = NA_real_, auc_ci = c(NA_real_, NA_real_, NA_real_),
                  predicted_probs = NULL, y_true = NULL, benchmark_vals = NULL,
                  note = "Degenerate class distribution — not fitted."))
    }
    X_sub <- X_main[mask, , drop = FALSE]
    res <- run_nested_loocv_svm(
      X = X_sub, y = factor(as.character(y_sub), levels = c(neg, pos)),
      C_grid = C_grid, gamma_grid = gamma_grid,
      inner_folds = inner_folds, seed = seed
    )
    metrics <- compute_classification_metrics(res$y_true, res$predicted_probs, pos)
    small_n <- n_s < 20L
    inverted <- !is.na(metrics$auc) && metrics$auc < 0.5
    notes <- sprintf("SVM-RBF nested-LOOCV on KI67 gate, %s subgroup.", label)
    if (small_n) {
      notes <- paste0(notes,
                      sprintf(" CAVEAT: n=%d is small (< 20); estimate is indicative only.",
                              n_s))
    }
    if (inverted) {
      notes <- paste0(notes,
                      " WARNING: AUC < 0.5 means the gate-features ordering inverts on this subgroup — interpret as no signal (or anti-correlation) on the subset, not negative information.")
    }
    list(
      n                 = n_s,
      n_responder       = sum(y_sub == pos),
      auc               = round(metrics$auc, 4),
      auc_ci            = round(metrics$auc_ci, 4),
      balanced_accuracy = round(metrics$balanced_accuracy, 4),
      sensitivity       = round(metrics$sensitivity, 4),
      specificity       = round(metrics$specificity, 4),
      small_n_caveat    = small_n,
      inverted_auc      = inverted,
      predicted_probs   = res$predicted_probs,
      y_true            = as.character(res$y_true),
      benchmark_vals    = as.numeric(bench_vals[mask]),
      note              = notes
    )
  }

  sub_low  <- run_subgroup_svm(low_mask,  sprintf("%s < %g", benchmark_label, high_threshold))
  sub_high <- run_subgroup_svm(high_mask, sprintf("%s >= %g", benchmark_label, high_threshold))

  list(
    label              = benchmark_label,
    column             = benchmark_col,
    n_valid            = n_valid,
    bin_breaks         = bin_breaks,
    bin_labels         = bin_labels,
    bin_crosstab       = bin_summary,
    fisher_3bins_p     = round(fish3, 4),
    ca_trend_p         = round(ca_p, 4),
    binary_cut         = binary_cut,
    subgroup_low       = sub_low,
    subgroup_high      = sub_high,
    note               = sprintf(
      "Stratified analysis on n=%d valid %s values. Categorical bins at %s. Binary clinical cut at %g. Subgroup SVM uses LMM gate features.",
      n_valid, benchmark_label,
      paste(bin_breaks, collapse = "/"), high_threshold
    )
  )
}

#' @title Combined Benchmark + Model Information Gain (IDI/NRI)
#' @description
#' Quantifies how much the LMM-gate features add over a continuous clinical
#' benchmark biomarker (e.g. PD-L1). Two complementary perspectives:
#'   1. Logistic (apparent/in-sample): benchmark alone vs benchmark + main
#'      features, scored with apparent (resubstitution) probabilities — the
#'      conventional setup for Pencina IDI / continuous NRI. Patient-level
#'      bootstrap refits both models inside each iteration to give an honest
#'      CI accounting for fit uncertainty. We deliberately avoid LOOCV here
#'      because univariate LOOCV on a weak baseline (e.g. PD-L1 with AUC ~ 0.6)
#'      is unstable: fold-to-fold sign flips of the slope can drive the
#'      "benchmark-alone" AUC below chance, which mechanically inflates IDI.
#'   2. SVM-RBF nested-LOOCV: same-subset comparison of model-with-benchmark
#'      vs model-without (informs whether the benchmark adds noise to the
#'      non-linear classifier).
#'
#' @param DATA Step 01 DATA object (for metadata).
#' @param X_main Numeric matrix of main-effect features (already z-scored at
#'   Step 01) aligned with DATA$metadata$Patient_ID.
#' @param y Factor outcome aligned with X_main.
#' @param positive_label Positive class label.
#' @param input_file Path to raw clinical Excel (T0).
#' @param benchmark_col Column name in raw Excel.
#' @param benchmark_label Optional label.
#' @param n_boot Number of patient-level bootstrap iterations for IDI/cNRI CI.
#' @param C_grid,gamma_grid,inner_folds SVM tuning.
#' @param seed RNG seed.
#' @return List with logistic, idi/nri, and SVM comparison. NULL if benchmark
#'   unavailable.
#' @export
run_combined_benchmark_model <- function(DATA, X_main, y, positive_label,
                                         input_file, benchmark_col,
                                         benchmark_label = NULL,
                                         n_boot          = 1000L,
                                         C_grid          = c(0.01, 0.1, 1, 10, 100),
                                         gamma_grid      = c(0.01, 0.1, 1, 10),
                                         inner_folds     = 5L,
                                         seed            = 2026) {
  if (!requireNamespace("readxl", quietly = TRUE)) return(NULL)
  if (!requireNamespace("dplyr",  quietly = TRUE)) return(NULL)
  if (!requireNamespace("pROC",   quietly = TRUE)) return(NULL)

  if (is.null(benchmark_label)) benchmark_label <- benchmark_col
  if (!file.exists(input_file)) return(NULL)

  df_raw <- tryCatch(readxl::read_excel(input_file, sheet = 1),
                     error = function(e) NULL)
  if (is.null(df_raw) || !benchmark_col %in% names(df_raw)) return(NULL)

  df_bench <- df_raw %>%
    dplyr::select(Patient_ID, dplyr::all_of(benchmark_col)) %>%
    dplyr::mutate(Patient_ID = as.character(Patient_ID))
  df_meta <- DATA$metadata %>%
    dplyr::mutate(Patient_ID = as.character(Patient_ID)) %>%
    dplyr::left_join(df_bench, by = "Patient_ID")

  bench_vals <- as.numeric(df_meta[[benchmark_col]])
  valid_mask <- !is.na(bench_vals)
  if (sum(valid_mask) < 20) return(NULL)

  pos <- positive_label
  neg <- setdiff(levels(y), pos)
  y_sub <- y[valid_mask]
  y_bin <- as.integer(y_sub == pos)
  x_bench <- bench_vals[valid_mask]
  X_sub  <- X_main[valid_mask, , drop = FALSE]

  # ---- Apparent (in-sample) logistic: benchmark alone vs benchmark + features
  # We use resubstitution predictions (the standard Pencina 2008 IDI setup).
  # CI honesty is delivered by the patient-level bootstrap below, which refits
  # both models on each resample.
  df_lr <- data.frame(y = y_bin, bench = x_bench, X_sub, check.names = FALSE)
  feat_names_lr <- make.names(colnames(X_sub), unique = TRUE)
  colnames(df_lr)[3:ncol(df_lr)] <- feat_names_lr

  formula_bench <- as.formula("y ~ bench")
  formula_comb  <- as.formula(paste("y ~ bench +", paste(feat_names_lr, collapse = " + ")))

  insample_glm <- function(formula_obj, df) {
    fit <- suppressWarnings(glm(formula_obj, data = df, family = binomial()))
    as.numeric(predict(fit, newdata = df, type = "response"))
  }

  p_bench <- insample_glm(formula_bench, df_lr)
  p_comb  <- insample_glm(formula_comb,  df_lr)

  roc_bench <- pROC::roc(y_bin, p_bench, direction = "<", quiet = TRUE)
  roc_comb  <- pROC::roc(y_bin, p_comb,  direction = "<", quiet = TRUE)
  auc_bench <- as.numeric(pROC::auc(roc_bench))
  auc_comb  <- as.numeric(pROC::auc(roc_comb))
  ci_bench  <- tryCatch(as.numeric(pROC::ci.auc(roc_bench, method = "delong")),
                        error = function(e) c(NA_real_, auc_bench, NA_real_))
  ci_comb   <- tryCatch(as.numeric(pROC::ci.auc(roc_comb,  method = "delong")),
                        error = function(e) c(NA_real_, auc_comb,  NA_real_))
  delong_p  <- tryCatch(pROC::roc.test(roc_bench, roc_comb,
                                       method = "delong", paired = TRUE)$p.value,
                        error = function(e) NA_real_)

  # ---- IDI / Continuous NRI ----
  idi_components <- function(p_old, p_new, y_bin) {
    evt <- y_bin == 1
    imp_evt <- mean(p_new[evt]) - mean(p_old[evt])
    red_nev <- mean(p_old[!evt]) - mean(p_new[!evt])
    list(IDI = imp_evt + red_nev, IDI_event = imp_evt, IDI_nonevent = red_nev)
  }
  cnri_components <- function(p_old, p_new, y_bin) {
    evt <- y_bin == 1
    up_evt <- mean((p_new - p_old)[evt] > 0) - mean((p_new - p_old)[evt] < 0)
    dn_nev <- mean((p_new - p_old)[!evt] < 0) - mean((p_new - p_old)[!evt] > 0)
    list(NRI = up_evt + dn_nev, NRI_event = up_evt, NRI_nonevent = dn_nev)
  }

  idi <- idi_components(p_bench, p_comb, y_bin)
  nri <- cnri_components(p_bench, p_comb, y_bin)

  # Patient-level bootstrap CI for IDI / cNRI. Refits both logistic models
  # inside each iteration so the CI reflects fit-uncertainty in addition to
  # sampling variability (resampling predicted-probability tuples alone would
  # understate variance with apparent predictions).
  set.seed(seed)
  n_y <- length(y_bin)
  idi_boot <- numeric(n_boot); nri_boot <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    idx <- sample(n_y, n_y, replace = TRUE)
    if (length(unique(y_bin[idx])) < 2) { idi_boot[b] <- NA_real_; nri_boot[b] <- NA_real_; next }
    df_b <- df_lr[idx, , drop = FALSE]
    pb <- tryCatch(insample_glm(formula_bench, df_b), error = function(e) NULL)
    pc <- tryCatch(insample_glm(formula_comb,  df_b), error = function(e) NULL)
    if (is.null(pb) || is.null(pc)) { idi_boot[b] <- NA_real_; nri_boot[b] <- NA_real_; next }
    yb <- y_bin[idx]
    idi_boot[b] <- idi_components(pb, pc, yb)$IDI
    nri_boot[b] <- cnri_components(pb, pc, yb)$NRI
  }
  idi_ci <- as.numeric(quantile(idi_boot, c(0.025, 0.975), na.rm = TRUE))
  nri_ci <- as.numeric(quantile(nri_boot, c(0.025, 0.975), na.rm = TRUE))
  idi_p  <- 2 * min(mean(idi_boot <= 0, na.rm = TRUE), mean(idi_boot > 0, na.rm = TRUE))
  nri_p  <- 2 * min(mean(nri_boot <= 0, na.rm = TRUE), mean(nri_boot > 0, na.rm = TRUE))

  # ---- SVM-RBF on the same subset: benchmark + features vs features only ----
  bench_z <- as.numeric(scale(x_bench))  # standardize benchmark to feature scale
  X_with  <- cbind(BENCHMARK = bench_z, X_sub)
  res_with <- run_nested_loocv_svm(X = X_with, y = y_sub,
                                   C_grid = C_grid, gamma_grid = gamma_grid,
                                   inner_folds = inner_folds, seed = seed)
  res_only <- run_nested_loocv_svm(X = X_sub, y = y_sub,
                                   C_grid = C_grid, gamma_grid = gamma_grid,
                                   inner_folds = inner_folds, seed = seed)
  m_with <- compute_classification_metrics(res_with$y_true, res_with$predicted_probs, pos)
  m_only <- compute_classification_metrics(res_only$y_true, res_only$predicted_probs, pos)
  roc_with <- tryCatch(pROC::roc(as.integer(res_with$y_true == pos),
                                 res_with$predicted_probs, direction = "<", quiet = TRUE),
                       error = function(e) NULL)
  roc_only <- tryCatch(pROC::roc(as.integer(res_only$y_true == pos),
                                 res_only$predicted_probs, direction = "<", quiet = TRUE),
                       error = function(e) NULL)
  svm_delong_p <- if (!is.null(roc_with) && !is.null(roc_only)) {
    tryCatch(pROC::roc.test(roc_only, roc_with, method = "delong", paired = TRUE)$p.value,
             error = function(e) NA_real_)
  } else NA_real_

  list(
    label   = benchmark_label,
    column  = benchmark_col,
    n_valid = sum(valid_mask),
    logistic = list(
      auc_benchmark           = round(auc_bench, 4),
      auc_benchmark_ci        = as.list(round(ci_bench, 4)),
      auc_combined            = round(auc_comb, 4),
      auc_combined_ci         = as.list(round(ci_comb, 4)),
      delta_auc               = round(auc_comb - auc_bench, 4),
      delong_p                = round(delong_p, 5)
    ),
    information_gain = list(
      IDI                = round(idi$IDI, 4),
      IDI_event          = round(idi$IDI_event, 4),
      IDI_nonevent       = round(idi$IDI_nonevent, 4),
      IDI_95CI           = as.list(round(idi_ci, 4)),
      IDI_bootstrap_p    = round(idi_p, 4),
      cNRI               = round(nri$NRI, 4),
      cNRI_event         = round(nri$NRI_event, 4),
      cNRI_nonevent      = round(nri$NRI_nonevent, 4),
      cNRI_95CI          = as.list(round(nri_ci, 4)),
      cNRI_bootstrap_p   = round(nri_p, 4),
      n_boot             = n_boot
    ),
    svm_comparison = list(
      auc_features_only        = round(m_only$auc, 4),
      auc_features_with_bench  = round(m_with$auc, 4),
      delta_auc                = round(m_with$auc - m_only$auc, 4),
      delong_p                 = round(svm_delong_p, 5),
      note = sprintf(
        "Same n=%d subset. Negative delta means adding %s to the SVM-RBF kernel inputs degrades AUC — interpretable as the benchmark contributing noise relative to the LMM gate.",
        sum(valid_mask), benchmark_label
      )
    ),
    note = sprintf(
      "Apparent (in-sample) logistic: %s alone vs %s + %d gate marker(s). IDI/cNRI estimated with patient-level bootstrap (B=%d) that refits both models per resample.",
      benchmark_label, benchmark_label, ncol(X_sub), n_boot
    ),
    plot_data = list(
      y_bin                      = y_bin,
      positive_label             = pos,
      p_bench_logistic           = p_bench,
      p_combined_logistic        = p_comb,
      p_svm_features_only        = res_only$predicted_probs,
      p_svm_features_with_bench  = res_with$predicted_probs,
      auc_svm_features_only      = round(m_only$auc, 4),
      auc_svm_features_with_bench= round(m_with$auc, 4)
    )
  )
}

# ==============================================================================
# CLINICAL-UTILITY LAYER (calibration + decision curve + optimism correction)
# ==============================================================================
#' Clinical-utility evaluation of a pre-specified composite biomarker.
#'
#' Builds a single pre-specified composite score (mean z of the supplied gate
#' markers) at T0, fits a 1-feature logistic, and evaluates it the way a
#' clinical-prediction manuscript (TRIPOD/REMARK) requires: discrimination
#' (apparent + leakage-free LOO AUC + permutation p), Harrell enhanced-bootstrap
#' OPTIMISM correction, CALIBRATION (intercept/slope/Brier, apparent + LOO) with
#' logistic RECALIBRATION, and a DECISION CURVE (net benefit vs treat-all/none).
#'
#' Dataset-agnostic: the caller supplies `gate_markers` and `gate_provenance`;
#' nothing is re-derived here. When the gate comes from the same T0 data
#' (`gate_provenance = "univariate-selected"`) the composite is NOT independent —
#' the returned `provenance_note` flags it as exploratory/conditional-on-selection
#' so the optimism≈0 claim is never overstated.
#'
#' @param DATA_T0       Step 01 standard RDS list (`hybrid_data_z`, `metadata`).
#' @param gate_markers  Character vector of marker names forming the composite.
#' @param resp_label    Positive (responder) class label, e.g. "RP".
#' @param gate_provenance "lmm-prespecified" or "univariate-selected".
#' @param n_boot,n_perm Bootstrap / permutation reps (default 2000 each).
#' @param min_n         Minimum sample size; below this the block is skipped.
#' @param seed          RNG seed (state saved/restored — no side effects).
#' @return Named list (scalars + `decision_curve`, `per_patient` data.frames), or
#'         `NULL` (with an explanatory message) when prerequisites are not met.
run_clinical_utility <- function(DATA_T0, gate_markers, resp_label,
                                 gate_provenance = c("lmm-prespecified",
                                                     "univariate-selected"),
                                 n_boot = 2000L, n_perm = 2000L,
                                 min_n = 20L, seed = 2026L) {

  gate_provenance <- match.arg(gate_provenance)

  # ── Reproducible RNG without disturbing the caller's stream ─────────────────
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  }
  set.seed(seed)

  # ── Assemble composite + outcome ────────────────────────────────────────────
  z      <- as.data.frame(DATA_T0$hybrid_data_z)
  avail  <- intersect(gate_markers, colnames(z))
  group  <- DATA_T0$metadata$Group
  if (length(avail) == 0) {
    message("[ML][CU] No gate markers present in T0 matrix — clinical utility skipped.")
    return(NULL)
  }
  neg_label <- setdiff(unique(group), resp_label)
  if (length(neg_label) != 1) {
    message("[ML][CU] Outcome is not binary — clinical utility skipped.")
    return(NULL)
  }
  Zc <- z[, avail, drop = FALSE]
  Zc[] <- lapply(Zc, function(x) { x <- as.numeric(x); x[is.na(x)] <- median(x, na.rm = TRUE); x })
  comp <- rowMeans(Zc)
  keep <- !is.na(comp) & group %in% c(resp_label, neg_label)
  comp <- comp[keep]; grp <- group[keep]
  yb   <- as.integer(grp == resp_label)
  y    <- factor(grp, levels = c(neg_label, resp_label))
  n    <- length(yb)
  if (n < min_n || length(unique(yb)) < 2) {
    message(sprintf("[ML][CU] n=%d (< min_n=%d) or single class — clinical utility skipped.",
                    n, min_n))
    return(NULL)
  }
  message(sprintf("[ML][CU] composite=%s | n=%d (%s=%d, %s=%d) | provenance=%s",
                  paste(avail, collapse = "+"), n, resp_label, sum(yb),
                  neg_label, sum(1 - yb), gate_provenance))

  # ── Local helpers (orientation fixed: positive class = resp_label) ──────────
  auc_pos <- function(yv, p) as.numeric(pROC::auc(pROC::roc(
    yv, p, levels = c(neg_label, resp_label), direction = "<", quiet = TRUE)))
  brier   <- function(yvb, p) mean((p - yvb)^2)
  calib_params <- function(yvb, p) {           # CITL + slope by logistic recalibration
    p  <- pmin(pmax(p, 1e-6), 1 - 1e-6); lp <- qlogis(p)
    slope <- tryCatch(coef(glm(yvb ~ lp, family = binomial()))[["lp"]],
                      error = function(e) NA_real_)
    citl  <- tryCatch(coef(glm(yvb ~ 1, family = binomial(), offset = lp))[[1]],
                      error = function(e) NA_real_)
    c(intercept = citl, slope = slope)
  }
  fac <- function(v) factor(ifelse(v == 1, resp_label, neg_label),
                            levels = c(neg_label, resp_label))

  # ── 1-feature logistic: apparent + leakage-free LOO probabilities ───────────
  df  <- data.frame(yb = yb, comp = comp)
  fit <- glm(yb ~ comp, data = df, family = binomial())
  p_app <- as.numeric(predict(fit, type = "response"))
  p_loo <- numeric(n)
  for (i in seq_len(n)) {
    fi <- suppressWarnings(glm(yb ~ comp, data = df[-i, , drop = FALSE], family = binomial()))
    p_loo[i] <- as.numeric(predict(fi, newdata = df[i, , drop = FALSE], type = "response"))
  }
  auc_app <- auc_pos(y, p_app); auc_loo <- auc_pos(y, p_loo)
  br_app  <- brier(yb, p_app);  br_loo  <- brier(yb, p_loo)
  cp_app  <- calib_params(yb, p_app); cp_loo <- calib_params(yb, p_loo)

  # permutation p for the LOO AUC (shuffle labels, full LOO refit)
  cnt <- 1L
  for (b in seq_len(n_perm)) {
    yp <- sample(yb); dfp <- data.frame(yb = yp, comp = comp); pl <- numeric(n)
    for (i in seq_len(n)) {
      fi <- suppressWarnings(glm(yb ~ comp, data = dfp[-i, ], family = binomial()))
      pl[i] <- as.numeric(predict(fi, newdata = dfp[i, ], type = "response"))
    }
    if (auc_pos(fac(yp), pl) >= auc_loo) cnt <- cnt + 1L
  }
  perm_p_loo <- cnt / (n_perm + 1)

  # ── Harrell enhanced bootstrap: optimism of AUC & calibration slope ─────────
  opt_auc <- opt_slope <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    idx <- sample(n, n, replace = TRUE); dfb <- df[idx, , drop = FALSE]
    fb  <- suppressWarnings(glm(yb ~ comp, data = dfb, family = binomial()))
    pb_boot <- as.numeric(predict(fb, newdata = dfb, type = "response"))
    pb_orig <- as.numeric(predict(fb, newdata = df,  type = "response"))
    a_boot  <- tryCatch(auc_pos(fac(dfb$yb), pb_boot), error = function(e) NA_real_)
    opt_auc[b]   <- a_boot - auc_pos(y, pb_orig)
    opt_slope[b] <- calib_params(dfb$yb, pb_boot)["slope"] - calib_params(yb, pb_orig)["slope"]
  }
  opt_auc_m   <- mean(opt_auc,   na.rm = TRUE)
  opt_slope_m <- mean(opt_slope, na.rm = TRUE)
  auc_corr    <- auc_app - opt_auc_m
  auc_corr_ci <- as.numeric(quantile(auc_app - opt_auc, c(.025, .975), na.rm = TRUE))
  slope_corr  <- cp_app[["slope"]] - opt_slope_m

  # ── Logistic recalibration of the LOO probabilities (shrinkage) ─────────────
  lp_loo  <- qlogis(pmin(pmax(p_loo, 1e-6), 1 - 1e-6))
  p_recal <- as.numeric(plogis(cp_loo[["intercept"]] + cp_loo[["slope"]] * lp_loo))
  cp_recal <- calib_params(yb, p_recal); br_recal <- brier(yb, p_recal)

  # ── Decision curve (net benefit on honest LOO probabilities) ────────────────
  prev    <- mean(yb)
  pt_grid <- seq(0.01, 0.60, by = 0.01)
  dc <- do.call(rbind, lapply(pt_grid, function(pt) {
    pos <- p_loo >= pt; w <- pt / (1 - pt)
    data.frame(threshold = pt,
               model      = sum(pos & yb == 1) / n - sum(pos & yb == 0) / n * w,
               treat_all  = prev - (1 - prev) * w,
               treat_none = 0)
  }))
  dom <- dc[dc$model > dc$treat_all & dc$model > dc$treat_none, , drop = FALSE]
  dom_range <- if (nrow(dom))
    sprintf("%.0f%%-%.0f%%", 100 * min(dom$threshold), 100 * max(dom$threshold)) else "none"

  # ── Youden cutpoint on the composite score ──────────────────────────────────
  roc_comp <- pROC::roc(y, comp, levels = c(neg_label, resp_label),
                        direction = "<", quiet = TRUE)
  yj_df <- as.data.frame(pROC::coords(roc_comp, "best", best.method = "youden",
                                      ret = c("threshold", "sensitivity", "specificity"),
                                      transpose = FALSE))
  yj <- as.numeric(yj_df[1, c("threshold", "sensitivity", "specificity")])

  provenance_note <- if (gate_provenance == "lmm-prespecified") {
    paste("Composite pre-specified from the Step 04 LMM gate (longitudinal data,",
          "disjoint from this T0 classifier): optimism correction is a genuine",
          "stability property, not an in-sample artefact.")
  } else {
    paste("Composite built from markers selected on these same T0 data",
          "(univariate gate): EXPLORATORY / conditional-on-selection — the",
          "optimism estimate is a lower bound and calibration is not independent.")
  }
  message(sprintf("[ML][CU] AUC app=%.3f LOO=%.3f (perm p=%.4f) | corrected=%.3f [%.3f-%.3f] | optimism=%.4f",
                  auc_app, auc_loo, perm_p_loo, auc_corr, auc_corr_ci[1], auc_corr_ci[2], opt_auc_m))
  message(sprintf("[ML][CU] calib LOO: CITL=%.2f slope=%.2f Brier=%.3f | DCA dominates over pt [%s]",
                  cp_loo[["intercept"]], cp_loo[["slope"]], br_loo, dom_range))

  list(
    composite_markers = avail,
    gate_provenance   = gate_provenance,
    provenance_note   = provenance_note,
    n                 = n, n_pos = sum(yb), n_neg = sum(1 - yb),
    prevalence        = prev,
    positive_label    = resp_label, negative_label = neg_label,
    discrimination = list(
      auc_apparent = auc_app, auc_loo = auc_loo, loo_perm_p = perm_p_loo,
      n_perm = n_perm, auc_optimism = opt_auc_m,
      auc_corrected = auc_corr, auc_corrected_ci = auc_corr_ci),
    calibration = list(
      apparent    = list(intercept = cp_app[["intercept"]], slope = cp_app[["slope"]], brier = br_app),
      loo         = list(intercept = cp_loo[["intercept"]], slope = cp_loo[["slope"]], brier = br_loo),
      recalibrated = list(intercept = cp_recal[["intercept"]], slope = cp_recal[["slope"]], brier = br_recal),
      slope_optimism = opt_slope_m, slope_corrected = slope_corr),
    youden = list(threshold = yj[1], sensitivity = yj[2], specificity = yj[3]),
    decision_curve = list(prevalence = prev, dominance_range = dom_range, curve = dc),
    per_patient = data.frame(
      Patient_ID  = as.character(DATA_T0$metadata$Patient_ID[keep]),
      True_Group  = as.character(grp),
      Composite   = round(comp, 4),
      Prob_LOO    = round(p_loo, 4),
      Prob_Recal  = round(p_recal, 4),
      stringsAsFactors = FALSE)
  )
}

#' JSON-serializable view of a clinical-utility result (drops the per-patient
#' table, which is written to Excel instead). Returns NULL passthrough.
clinical_utility_json <- function(cu) {
  if (is.null(cu)) return(NULL)
  cu$per_patient <- NULL
  cu
}

#' Flat Metric/Value summary of a clinical-utility result, for the Excel report.
write_clinical_utility_summary <- function(cu) {
  d <- cu$discrimination; cb <- cu$calibration; dc <- cu$decision_curve; yj <- cu$youden
  data.frame(
    Metric = c("Composite markers", "Gate provenance", "N (pos/neg)", "Prevalence",
               "AUC apparent", "AUC LOO", "LOO permutation p", "AUC optimism",
               "AUC corrected", "AUC corrected CI",
               "Calibration CITL (LOO)", "Calibration slope (LOO)", "Brier (LOO)",
               "Calibration slope (recalibrated)", "Brier (recalibrated)",
               "Decision-curve dominance range",
               "Youden threshold", "Youden sensitivity", "Youden specificity",
               "Provenance note"),
    Value = c(
      paste(cu$composite_markers, collapse = "+"), cu$gate_provenance,
      sprintf("%d (%d/%d)", cu$n, cu$n_pos, cu$n_neg), sprintf("%.3f", cu$prevalence),
      sprintf("%.3f", d$auc_apparent), sprintf("%.3f", d$auc_loo),
      sprintf("%.4f", d$loo_perm_p), sprintf("%.4f", d$auc_optimism),
      sprintf("%.3f", d$auc_corrected),
      sprintf("[%.3f, %.3f]", d$auc_corrected_ci[1], d$auc_corrected_ci[2]),
      sprintf("%.3f", cb$loo$intercept), sprintf("%.3f", cb$loo$slope),
      sprintf("%.3f", cb$loo$brier), sprintf("%.3f", cb$recalibrated$slope),
      sprintf("%.3f", cb$recalibrated$brier), dc$dominance_range,
      sprintf("%.3f", yj$threshold), sprintf("%.3f", yj$sensitivity),
      sprintf("%.3f", yj$specificity), cu$provenance_note),
    stringsAsFactors = FALSE)
}

# =============================================================================
# Predictor-timepoint adapter for the added-value layer
# -----------------------------------------------------------------------------
#' Re-express a processed cohort at a chosen predictor timepoint.
#'
#' The immune composite (and the random-marker specificity pool) is normally
#' evaluated on baseline T0 z-scores. For longitudinal cohorts the SAME gate can be
#' read at the on-treatment value (`"T1"`) or the on-treatment change
#' (`"delta"` = T1 − T0). This returns a `DATA_T0`-shaped object
#' (`hybrid_data_z` + `metadata`) for the requested timepoint so the entire
#' added-value machinery downstream is timepoint-agnostic and unchanged.
#'
#' `"T0"` returns `DATA_T0` untouched (canonical path, byte-for-byte). `"T1"`/`"delta"`
#' are built from the paired (T0+T1) patients of `DATA_LONG` using the joint-z
#' convention (matches run_gate_signal_decomposition: delta = z_t1 − z_t0). Returns
#' NULL when paired data is unavailable so the caller skips gracefully; an unknown
#' timepoint falls back to T0.
build_timepoint_composite_data <- function(DATA_T0, DATA_LONG, timepoint = "T0") {
  if (identical(timepoint, "T0")) return(DATA_T0)
  if (!timepoint %in% c("T1", "delta")) {
    message(sprintf("[ML][AV] unknown composite_timepoint '%s' — using T0.", timepoint))
    return(DATA_T0)
  }
  if (is.null(DATA_LONG) || is.null(DATA_LONG$hybrid_data_z)) {
    message(sprintf("[ML][AV] composite_timepoint='%s' needs longitudinal data — skipped.", timepoint))
    return(NULL)
  }
  META    <- c("Patient_ID", "Sample_ID", "Timepoint", "Group")
  markers <- setdiff(colnames(DATA_LONG$hybrid_data_z), META)
  z <- as.data.frame(DATA_LONG$hybrid_data_z[, markers, drop = FALSE])
  z$Patient_ID <- as.character(DATA_LONG$metadata$Patient_ID)
  z$Timepoint  <- DATA_LONG$metadata$Timepoint
  z$Group      <- DATA_LONG$metadata$Group
  z0 <- z[z$Timepoint == "T0", , drop = FALSE]
  z1 <- z[z$Timepoint == "T1", , drop = FALSE]
  paired <- intersect(z0$Patient_ID, z1$Patient_ID)
  if (length(paired) < 10) {
    message(sprintf("[ML][AV] only %d paired patients — composite_timepoint='%s' skipped.",
                    length(paired), timepoint))
    return(NULL)
  }
  z0 <- z0[match(paired, z0$Patient_ID), , drop = FALSE]
  z1 <- z1[match(paired, z1$Patient_ID), , drop = FALSE]
  M  <- if (identical(timepoint, "T1")) z1[, markers, drop = FALSE]
        else as.data.frame(as.matrix(z1[, markers]) - as.matrix(z0[, markers]))  # delta = z_t1 − z_t0
  list(hybrid_data_z  = M,
       metadata       = data.frame(Patient_ID = paired, Group = z0$Group,
                                   stringsAsFactors = FALSE),
       hybrid_markers = markers)
}

# =============================================================================
# Clinical + cytometric ADDED-VALUE layer
# -----------------------------------------------------------------------------
# Pre-specified question: does the immune composite (gate markers) add predictive
# value BEYOND a standard-of-care clinical model? Additive — never alters the
# gate/X/classifiers. Compares three nested logistic models (clinical-only,
# immune-only, clinical+immune) with leakage-free LOO + in-fold imputation,
# Harrell optimism on the combined model, IDI/cNRI (apparent + bootstrap CI),
# DeLong on LOO probabilities, a likelihood-ratio test of the immune increment
# (asymptotic + exact permutation + Firth penalized — the correct nested-model
# test, since DeLong on AUC is insensitive for nested models), calibration and a
# decision curve (clinical vs combined), and an OS/PFS Cox secondary (C + LRT).
#
# `clinical_vars` is a named character vector: names = model-safe labels,
# values = exact column names in `input_file` (e.g. c(NLR="neutrophils/lymphocytes")).
# =============================================================================

# Row label for the 1-df immune composite on nomogram panel B. Single source so the
# spec, the var_scale/var_convert lookups and the figure caption cannot drift apart.
COMPOSITE_DISPLAY <- "Immune composite"

#' Exact display label for one nomogram tick.
#'
#' The nomogram is fitted on the model's own scale (z of logit(cell fraction), or the
#' Δ of those z's). Both map BIJECTIVELY back to a quantity a clinician reads, so a tick
#' can be RELABELLED exactly — no refit, no approximation, model untouched:
#'   "fc"  -> exp(value * sigma). Since logit(p) = log(p/(1-p)) and p is a cell fraction
#'           of a parent gate, p/(1-p) is the positive:negative cell RATIO (for the KI67
#'           gate: proliferating:resting). So exp(Δlogit) is EXACTLY the fold change in
#'           that ratio. Centering cancels in a Δ (Δz*sigma == Δlogit to machine
#'           precision), hence sigma alone inverts it — no mu needed. For a COMPOSITE
#'           (mean_i Δlogit_i/sigma_i) pass sigma = 1: exp(composite) is then the exact
#'           weighted geometric-mean ratio fold change, weights w_i = 1/(k*sigma_i).
#'           NB this is the fold change of the RATIO, not of the fraction: the latter is
#'           FC = exp(Δlogit)*(1-p1)/(1-p0) and is NOT identified from Δlogit alone.
#'           For low-frequency subsets the two agree closely (r ~ 0.999) — say so in the
#'           caption, but never relabel the axis as the fraction fold change.
#'   "pct" -> plogis(value * sigma + mu) * 100, the cell fraction (%) itself. Only
#'           meaningful at a single timepoint (no fold change exists at baseline).
#' Falls back to the plain numeric value when no conversion applies.
#' @param v          numeric tick value(s), on the model scale.
#' @param conv       list(type="fc", sigma=) | list(type="pct", mu=, sigma=) | NULL.
#' @param scale_type the predictor's `var_scale` entry (used only by the fallback).
nomo_tick_label <- function(v, conv, scale_type) {
  if (!is.null(conv) && identical(conv$type, "fc")) {
    f <- exp(v * conv$sigma)
    return(ifelse(f >= 10, sprintf("%.1fx", f), sprintf("%.2fx", f)))
  }
  if (!is.null(conv) && identical(conv$type, "pct"))
    return(sprintf("%.1f%%", stats::plogis(v * conv$sigma + conv$mu) * 100))
  if (identical(scale_type, "z")) return(sprintf("%.1f", v))
  ifelse(abs(v - round(v)) < 1e-6, sprintf("%.0f", v), sprintf("%.1f", v))
}

#' Build a display-only nomogram spec from a logistic fit on individual predictors.
#'
#' Fits `glm(y ~ x1 + ... + xk, binomial)` on the supplied (already-imputed)
#' predictor matrix and returns everything the ggplot builder needs to draw a
#' classic points-based nomogram: per-predictor points axes, a total-points axis
#' and the total-points→predicted-probability mapping. This is an APPARENT model
#' display (fit on the full analytic sample) — a communication/deployment aid,
#' NOT a validated performance claim (those live in the added-value LOO layer).
#' The point scale follows the rms convention: the predictor whose linear-predictor
#' contribution spans the widest observed range gets 0–`points_max` points, and
#' every other predictor is scaled to the same points-per-logit unit.
#'
#' @param y01        integer 0/1 outcome (1 = positive/`resp_label`).
#' @param X          data.frame of numeric predictors (imputed; one column each).
#' @param var_scale  optional named vector mapping display name -> "z" | "fc" | "pct" |
#'                   "raw" (drives value-label formatting + axis colour).
#' @param var_convert optional named list mapping display name -> an EXACT display-scale
#'                   conversion (see `nomo_tick_label()`). Relabels ticks only: the fit,
#'                   the points and the tick POSITIONS are untouched.
#' @param points_max points assigned to the widest-range predictor (default 100).
build_nomogram_spec <- function(y01, X, resp_label = NULL, var_scale = NULL,
                                var_convert = NULL, points_max = 100) {
  X    <- as.data.frame(X)
  disp <- colnames(X)
  cn   <- make.names(disp, unique = TRUE)
  colnames(X) <- cn
  dat  <- data.frame(.y = as.integer(y01), X, check.names = FALSE)
  fit  <- suppressWarnings(glm(reformulate(cn, ".y"), data = dat, family = binomial()))
  b    <- coef(fit); b0 <- b[["(Intercept)"]]
  betas <- b[cn]; betas[!is.finite(betas)] <- 0

  # Per-predictor DISPLAY range + linear-predictor contribution range. For
  # continuous predictors the range is robustified to the 5th–95th percentile so a
  # single skewed outlier (common for Δ z-scores) does not dominate the shared point
  # scale and compress the other axes; discrete/ordinal predictors keep min–max.
  # This affects the drawn axis only — the glm fit/coefficients are unchanged.
  info <- lapply(seq_along(cn), function(j) {
    v  <- X[[cn[j]]]; uq <- unique(v[is.finite(v)])
    if (length(uq) <= 6) { lo <- min(v, na.rm = TRUE); hi <- max(v, na.rm = TRUE) }
    else { q <- quantile(v, c(.05, .95), na.rm = TRUE, names = FALSE); lo <- q[1]; hi <- q[2] }
    if (!is.finite(lo) || !is.finite(hi) || hi <= lo) {
      lo <- min(v, na.rm = TRUE); hi <- max(v, na.rm = TRUE) }
    bj <- betas[[j]]; c_lo <- bj * lo; c_hi <- bj * hi
    list(lo = lo, hi = hi, beta = bj,
         cmin = min(c_lo, c_hi), range = abs(bj) * (hi - lo))
  })
  ranges    <- vapply(info, function(z) z$range, numeric(1))
  max_range <- max(ranges); if (!is.finite(max_range) || max_range <= 0) max_range <- 1
  scale     <- points_max / max_range
  base_lp   <- b0 + sum(vapply(info, function(z) z$cmin, numeric(1)))

  # Per-predictor tick tables (value -> points), on the shared 0..points_max scale.
  preds <- lapply(seq_along(cn), function(j) {
    z     <- info[[j]]
    uq    <- sort(unique(X[[cn[j]]][is.finite(X[[cn[j]]])]))
    mp    <- z$range * scale                    # axis length in points
    # Ordinal/integer predictors (few distinct levels, e.g. ECOG PS 0/1/2): tick at
    # the actual observed levels — pretty() would invent meaningless half-steps.
    # For continuous predictors thin the ticks on SHORT (near-null) axes so labels
    # don't overlap: a barely-contributing predictor gets a short axis + few ticks.
    ntick <- if (mp < 8) 2L else if (mp < 25) 3L else 5L
    ticks <- if (length(uq) <= 6) uq else pretty(c(z$lo, z$hi), n = ntick)
    ticks <- ticks[ticks >= z$lo - 1e-9 & ticks <= z$hi + 1e-9]
    if (length(ticks) < 2) ticks <- unique(c(z$lo, z$hi))
    if (mp < 8 && length(ticks) > 2) ticks <- range(ticks)   # collapsed axis → endpoints only
    st <- if (!is.null(var_scale)   && disp[j] %in% names(var_scale))   var_scale[[disp[j]]]   else "raw"
    cv <- if (!is.null(var_convert) && disp[j] %in% names(var_convert)) var_convert[[disp[j]]] else NULL
    data.frame(display    = disp[j],
               value      = ticks,                            # model scale (positions ticks)
               label      = nomo_tick_label(ticks, cv, st),   # display scale (printed)
               points     = (z$beta * ticks - z$cmin) * scale,
               scale_type = st,
               max_points = z$range * scale,
               stringsAsFactors = FALSE)
  })
  names(preds) <- disp

  total_max <- sum(ranges) * scale
  # Total-points axis ticks + total-points -> probability mapping. A thinned prob
  # set (logit spacing compresses mid-probabilities → labels would collide near 0.5).
  probs   <- c(.05, .1, .2, .35, .5, .65, .8, .9, .95)
  prob_tp <- (qlogis(probs) - base_lp) * scale
  keep_p  <- prob_tp >= -1e-6 & prob_tp <= total_max + 1e-6
  prob_axis <- data.frame(prob = probs[keep_p],
                          total_points = pmin(pmax(prob_tp[keep_p], 0), total_max))
  tp_ticks  <- pretty(c(0, total_max), n = 6)
  tp_ticks  <- tp_ticks[tp_ticks >= 0 & tp_ticks <= total_max]

  list(
    predictors       = preds,
    points_max       = points_max,
    total_max_points = total_max,
    tp_ticks         = tp_ticks,
    prob_axis        = prob_axis,
    base_lp          = base_lp,
    scale            = scale,
    coefficients     = as.list(b),
    positive_label   = resp_label,
    auc_apparent     = tryCatch(as.numeric(pROC::auc(pROC::roc(
                          y01, as.numeric(predict(fit, type = "response")),
                          direction = "<", quiet = TRUE))), error = function(e) NA_real_),
    n                = length(y01))
}

run_clinical_immune_added_value <- function(DATA_T0, gate_markers, resp_label,
                                            input_file, clinical_vars,
                                            gate_provenance = "lmm-prespecified",
                                            surv_cols = c(os = "OS", pfs = "PFS",
                                                          death = "Alive_0/Dead_1"),
                                            n_boot = 2000L, n_perm = 2000L,
                                            ridge = TRUE,
                                            baseline_sensitivity_cols = NULL,
                                            n_random_null = 0L,
                                            lmm_interaction = NULL,
                                            composite_timepoint = "T0", DATA_LONG = NULL,
                                            baseline_exclude = NULL, comparator_label = NULL,
                                            nomogram_clinical_vars = NULL,
                                            min_n = 30L, seed = 2026L) {
  if (!requireNamespace("readxl", quietly = TRUE) ||
      !requireNamespace("pROC",   quietly = TRUE) ||
      !requireNamespace("dplyr",  quietly = TRUE)) return(NULL)
  if (is.null(input_file) || !file.exists(input_file)) {
    message("[ML][AV] input_file missing — added-value layer skipped."); return(NULL)
  }
  # Re-express the cohort at the requested predictor timepoint (T0 unchanged →
  # canonical path byte-for-byte; T1/delta built from paired longitudinal data).
  if (!identical(composite_timepoint, "T0")) {
    DATA_T0 <- build_timepoint_composite_data(DATA_T0, DATA_LONG, composite_timepoint)
    if (is.null(DATA_T0)) return(NULL)
  }
  clinical_vars <- unlist(clinical_vars)
  if (length(clinical_vars) == 0) return(NULL)
  if (is.null(names(clinical_vars))) names(clinical_vars) <- clinical_vars

  if (exists(".Random.seed", envir = .GlobalEnv)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  }
  set.seed(seed)

  # ── Immune composite (mean z of gate markers; median-imputed) ───────────────
  z     <- as.data.frame(DATA_T0$hybrid_data_z)
  avail <- intersect(gate_markers, colnames(z))
  if (length(avail) == 0) { message("[ML][AV] no gate markers in T0 — skipped."); return(NULL) }
  group <- DATA_T0$metadata$Group
  neg_label <- setdiff(unique(group), resp_label)
  if (length(neg_label) != 1) { message("[ML][AV] outcome not binary — skipped."); return(NULL) }
  Zc <- z[, avail, drop = FALSE]
  Zc[] <- lapply(Zc, function(x) { x <- as.numeric(x); x[is.na(x)] <- median(x, na.rm = TRUE); x })
  composite <- rowMeans(Zc)

  # ── Clinical design from the raw file, joined by Patient_ID ─────────────────
  raw <- tryCatch(readxl::read_excel(input_file, sheet = 1), error = function(e) NULL)
  if (is.null(raw)) return(NULL)
  names(raw) <- trimws(names(raw))
  missing_cols <- setdiff(unname(clinical_vars), names(raw))
  if (length(missing_cols)) {
    message(sprintf("[ML][AV] clinical columns absent (%s) — skipped.",
                    paste(missing_cols, collapse = ", "))); return(NULL)
  }
  raw$Patient_ID <- as.character(raw$Patient_ID)
  meta_ids <- as.character(DATA_T0$metadata$Patient_ID)
  ridx <- match(meta_ids, raw$Patient_ID)
  Cmat <- sapply(clinical_vars, function(col) suppressWarnings(as.numeric(raw[[col]][ridx])))
  Cmat <- matrix(Cmat, nrow = length(meta_ids),
                 dimnames = list(NULL, names(clinical_vars)))

  # ── Keep patients with a defined outcome (clinical NAs handled by imputation) ─
  keep <- group %in% c(resp_label, neg_label) & !is.na(composite)
  Cmat <- Cmat[keep, , drop = FALSE]; comp <- composite[keep]
  grp  <- group[keep]; yb <- as.integer(grp == resp_label)
  y    <- factor(grp, levels = c(neg_label, resp_label))
  pid  <- meta_ids[keep]
  n    <- length(yb)

  # ── Optional: drop var(s) from the FORMAL clinical baseline, keep FULL model as a
  #    display-only comparator curve. `baseline_exclude` (labels) => the formal
  #    clinical/combined models + every downstream test (LRT/Firth/IDI/ridge/optimism/
  #    specificity/DCA) are recomputed on the reduced baseline (Cmat is reduced IN PLACE
  #    so all downstream references follow automatically); `Cmat_full` retains the full
  #    clinical set only to draw the comparator curve. Empty exclude => byte-for-byte.
  cl_all         <- colnames(Cmat)
  excl           <- intersect(as.character(unlist(baseline_exclude)), cl_all)
  Cmat_full      <- Cmat
  has_comparator <- length(excl) > 0 && !is.null(comparator_label) &&
                    length(setdiff(cl_all, excl)) > 0
  if (has_comparator) {
    Cmat <- Cmat[, setdiff(cl_all, excl), drop = FALSE]
    message(sprintf("[ML][AV] FORMAL baseline reduced to {%s}; comparator '%s' = {%s}",
                    paste(colnames(Cmat), collapse = "+"), comparator_label,
                    paste(cl_all, collapse = "+")))
  }
  n_cc <- sum(complete.cases(Cmat))
  if (n < min_n || length(unique(yb)) < 2) {
    message(sprintf("[ML][AV] n=%d (<%d) or single class — skipped.", n, min_n)); return(NULL)
  }
  message(sprintf("[ML][AV] clinical={%s} + immune composite={%s} | n=%d (cc=%d, %s=%d/%s=%d)",
                  paste(names(clinical_vars), collapse = "+"), paste(avail, collapse = "+"),
                  n, n_cc, resp_label, sum(yb), neg_label, sum(1 - yb)))

  cl_names <- colnames(Cmat)
  # ── Helpers ─────────────────────────────────────────────────────────────────
  impute_median <- function(M, train_rows) {
    med <- apply(M[train_rows, , drop = FALSE], 2, median, na.rm = TRUE)
    med[!is.finite(med)] <- 0
    for (j in seq_len(ncol(M))) M[is.na(M[, j]), j] <- med[j]
    M
  }
  fit_predict <- function(Cm, cv, train_rows, test_rows, use_immune) {
    Ci <- impute_median(Cm, train_rows)
    dat <- data.frame(yb = yb, Ci, comp = cv, check.names = TRUE)
    # Derive names from Cm's own columns so any clinical column SUBSET works
    # (identical to the full Cmat call, where colnames(Cm) == cl_names).
    cn  <- make.names(colnames(Cm), unique = TRUE); colnames(dat)[2:(1 + length(cn))] <- cn
    rhs <- if (use_immune) c(cn, "comp") else cn
    fit <- suppressWarnings(glm(reformulate(rhs, "yb"), data = dat[train_rows, , drop = FALSE],
                                family = binomial()))
    as.numeric(predict(fit, newdata = dat[test_rows, , drop = FALSE], type = "response"))
  }
  auc_pos <- function(p, yv = y) as.numeric(pROC::auc(pROC::roc(
    yv, p, levels = c(neg_label, resp_label), direction = "<", quiet = TRUE)))
  calib_params <- function(yvb, p) {
    p <- pmin(pmax(p, 1e-6), 1 - 1e-6); lp <- qlogis(p)
    c(intercept = tryCatch(coef(glm(yvb ~ 1, family = binomial(), offset = lp))[[1]], error = function(e) NA_real_),
      slope     = tryCatch(coef(glm(yvb ~ lp, family = binomial()))[["lp"]],          error = function(e) NA_real_))
  }
  brier <- function(yvb, p) mean((p - yvb)^2)
  all_rows <- seq_len(n)

  # ── Apparent clinical & combined fits (reused for predictions + LRT/Firth) ──
  Ci_app  <- impute_median(Cmat, all_rows)
  cn_app  <- make.names(cl_names, unique = TRUE)
  dat_app <- data.frame(yb = yb, Ci_app, comp = comp, check.names = TRUE)
  colnames(dat_app)[2:(1 + length(cn_app))] <- cn_app
  fc_app  <- suppressWarnings(glm(reformulate(cn_app, "yb"), data = dat_app, family = binomial()))
  fk_app  <- suppressWarnings(glm(reformulate(c(cn_app, "comp"), "yb"), data = dat_app, family = binomial()))

  # ── Apparent + leakage-free LOO probabilities ───────────────────────────────
  pa_clin <- as.numeric(predict(fc_app, type = "response"))
  pa_comb <- as.numeric(predict(fk_app, type = "response"))
  pa_imm  <- { dat <- data.frame(yb = yb, comp = comp)
               f <- suppressWarnings(glm(yb ~ comp, data = dat, family = binomial()))
               as.numeric(predict(f, type = "response")) }
  pl_clin <- pl_comb <- pl_imm <- numeric(n)
  for (i in all_rows) {
    tr <- all_rows[-i]
    pl_clin[i] <- fit_predict(Cmat, comp, tr, i, FALSE)
    pl_comb[i] <- fit_predict(Cmat, comp, tr, i, TRUE)
    fi <- suppressWarnings(glm(yb ~ comp, data = data.frame(yb = yb, comp = comp)[tr, ], family = binomial()))
    pl_imm[i] <- as.numeric(predict(fi, newdata = data.frame(comp = comp[i]), type = "response"))
  }
  # ── Comparator = FULL clinical (incl. excluded var) + immune composite — the
  #    'everything' model (e.g. 'clinical + immune + NLR'); display-only ROC/DCA curve,
  #    does NOT enter the formal increment (which stays immune over the reduced baseline).
  pl_comp <- pa_comp <- NULL
  if (has_comparator) {
    cnc     <- make.names(colnames(Cmat_full), unique = TRUE)
    Ci_full <- impute_median(Cmat_full, all_rows)
    dfc     <- data.frame(yb = yb, Ci_full, comp = comp, check.names = TRUE)
    colnames(dfc)[2:(1 + length(cnc))] <- cnc
    fca     <- suppressWarnings(glm(reformulate(c(cnc, "comp"), "yb"), data = dfc, family = binomial()))
    pa_comp <- as.numeric(predict(fca, type = "response"))
    pl_comp <- numeric(n)
    for (i in all_rows) pl_comp[i] <- fit_predict(Cmat_full, comp, all_rows[-i], i, TRUE)
  }
  # ── Apparent clinical sub-model probabilities (PD-L1 / PS / PD-L1+PS vs clinical+immune),
  #    retained as the DESCRIPTIVE PD-L1/PS/PD-L1+PS benchmark in the _Patients Excel sheet
  #    (diag_28: report this split descriptively, NOT as the formal LOO baseline). The
  #    exploratory Figure2CLONE that used to plot these was removed 2026-07-01. NB: leakage-
  #    free LOO of these weak / missing-heavy single clinical predictors (PD-L1 ~20% missing)
  #    is numerically DEGENERATE — per-fold coefficient sign is unstable and the LOO ranking
  #    inverts (AUC→0), the same pathology the PD-L1 benchmark below avoids with complete-case
  #    fits — hence these are APPARENT (in-sample) only.
  #    Deterministic (no RNG) and BEFORE the bootstrap loops → primary numbers reproduce.
  #    Columns selected BY NAME for new-dataset versatility; absent column → all-NA (skipped).
  app_submodel <- function(cols) {
    if (!all(cols %in% cl_names)) return(rep(NA_real_, n))
    Ci <- impute_median(Cmat[, cols, drop = FALSE], all_rows)
    cn <- make.names(cols, unique = TRUE)
    d  <- data.frame(yb = yb, Ci, check.names = TRUE); colnames(d)[2:(1 + length(cn))] <- cn
    f  <- suppressWarnings(glm(reformulate(cn, "yb"), data = d, family = binomial()))
    as.numeric(predict(f, type = "response"))
  }
  submodel_app <- list(PDL1 = app_submodel("PD_L1"), PS = app_submodel("PS"),
                       PDL1PS = app_submodel(c("PD_L1", "PS")))
  # Benchmark (PD-L1 = 1st clinical var) ALONE — APPARENT complete-case univariate
  # logistic. PD-L1 is a single, externally pre-fixed biomarker (no selection, 1 param)
  # → apparent ≈ honest AUC; LOO of this near-null predictor is numerically pathological
  # (flips the AUC), and median-imputing the ~20% missing also inverts it → complete
  # cases only, NA elsewhere. Reported as the cross-sectional companion-biomarker AUC.
  bench_raw <- as.numeric(Cmat[, 1]); cc <- which(is.finite(bench_raw))
  pl_bench <- rep(NA_real_, n)
  if (length(cc) >= 10) {
    fb <- suppressWarnings(glm(yb ~ x, data = data.frame(yb = yb, x = bench_raw)[cc, , drop = FALSE], family = binomial()))
    pl_bench[cc] <- as.numeric(predict(fb, newdata = data.frame(x = bench_raw[cc]), type = "response"))
  }
  auc <- list(
    clinical = list(apparent = auc_pos(pa_clin), loo = auc_pos(pl_clin)),
    immune   = list(apparent = auc_pos(pa_imm),  loo = auc_pos(pl_imm)),
    combined = list(apparent = auc_pos(pa_comb), loo = auc_pos(pl_comb)))

  # ── DeLong (LOO): clinical vs combined ──────────────────────────────────────
  r_clin <- pROC::roc(y, pl_clin, levels = c(neg_label, resp_label), direction = "<", quiet = TRUE)
  r_comb <- pROC::roc(y, pl_comb, levels = c(neg_label, resp_label), direction = "<", quiet = TRUE)
  delong_p <- tryCatch(pROC::roc.test(r_clin, r_comb, method = "delong", paired = TRUE)$p.value,
                       error = function(e) NA_real_)

  # ── IDI / continuous NRI (apparent + patient-level bootstrap CI) ─────────────
  idi_fun <- function(po, pn) { e <- yb == 1; (mean(pn[e]) - mean(po[e])) + (mean(po[!e]) - mean(pn[!e])) }
  nri_fun <- function(po, pn) { e <- yb == 1; d <- pn - po
    (mean(d[e] > 0) - mean(d[e] < 0)) + (mean(d[!e] < 0) - mean(d[!e] > 0)) }
  idi_app <- idi_fun(pa_clin, pa_comb); nri_app <- nri_fun(pa_clin, pa_comb)
  idi_bs <- nri_bs <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    idx <- sample(n, n, replace = TRUE)
    Cb <- Cmat[idx, , drop = FALSE]; cb <- comp[idx]; ybb <- yb[idx]
    Ci <- impute_median(Cb, seq_len(n))
    cn <- make.names(cl_names, unique = TRUE)
    db <- data.frame(yb = ybb, Ci, comp = cb); colnames(db)[2:(1 + length(cn))] <- cn
    fcl <- suppressWarnings(glm(reformulate(cn, "yb"), data = db, family = binomial()))
    fco <- suppressWarnings(glm(reformulate(c(cn, "comp"), "yb"), data = db, family = binomial()))
    po <- as.numeric(predict(fcl, type = "response")); pn <- as.numeric(predict(fco, type = "response"))
    e <- ybb == 1
    idi_bs[b] <- (mean(pn[e]) - mean(po[e])) + (mean(po[!e]) - mean(pn[!e]))
    d <- pn - po
    nri_bs[b] <- (mean(d[e] > 0) - mean(d[e] < 0)) + (mean(d[!e] < 0) - mean(d[!e] > 0))
  }
  qci <- function(v) as.numeric(quantile(v, c(.025, .975), na.rm = TRUE))

  # ── Likelihood-ratio test of the immune increment (asymptotic + permutation + Firth) ──
  # The correct nested-model test for the pre-specified composite's added value: DeLong on
  # AUC is insensitive for nested models, so it can be ns while the increment is real. The
  # permutation null is exact (no chi-square reliance at low events-per-variable); the Firth
  # penalized profile-likelihood test is separation-robust. All on the apparent models.
  lrt_aov  <- tryCatch(anova(fc_app, fk_app, test = "Chisq"), error = function(e) NULL)
  p_col    <- if (!is.null(lrt_aov)) grep("^Pr\\(", colnames(lrt_aov), value = TRUE) else character(0)
  lrt_stat <- if (!is.null(lrt_aov)) lrt_aov$Deviance[2]                 else NA_real_
  lrt_p    <- if (length(p_col))    lrt_aov[[p_col[1]]][2]               else NA_real_
  lrt_perm_p <- NA_real_
  if (is.finite(lrt_stat) && n_perm > 0) {
    ln <- numeric(n_perm)
    for (b in seq_len(n_perm)) {
      dp <- dat_app; dp$yb <- sample(yb)
      fcp <- suppressWarnings(glm(reformulate(cn_app, "yb"),         data = dp, family = binomial()))
      fkp <- suppressWarnings(glm(reformulate(c(cn_app, "comp"), "yb"), data = dp, family = binomial()))
      ln[b] <- deviance(fcp) - deviance(fkp)
    }
    lrt_perm_p <- (1 + sum(ln >= lrt_stat)) / (n_perm + 1)
  }
  lrt_firth_p <- if (requireNamespace("logistf", quietly = TRUE))
    tryCatch({ ff <- logistf::logistf(reformulate(c(cn_app, "comp"), "yb"), data = dat_app)
               unname(ff$prob[match("comp", names(coef(ff)))]) }, error = function(e) NA_real_)
    else NA_real_

  # ── Harrell optimism of the combined-model AUC ──────────────────────────────
  # Self-contained (fit_predict() binds the global original-order yb, so it must
  # not be used on a bootstrap resample): fit on the resample, then evaluate that
  # same model apparent-on-resample vs on the original cohort.
  cn_opt  <- make.names(cl_names, unique = TRUE)
  opt_auc <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    idx <- sample(n, n, replace = TRUE)
    Ci_boot <- impute_median(Cmat[idx, , drop = FALSE], seq_len(n))
    db <- data.frame(yb = yb[idx], Ci_boot, comp = comp[idx]); colnames(db)[2:(1 + length(cn_opt))] <- cn_opt
    fb <- suppressWarnings(glm(reformulate(c(cn_opt, "comp"), "yb"), data = db, family = binomial()))
    a_boot <- tryCatch(auc_pos(as.numeric(predict(fb, type = "response")),
                               factor(grp[idx], levels = c(neg_label, resp_label))),
                       error = function(e) NA_real_)
    Cf <- impute_median(Cmat, idx); df0 <- data.frame(Cf, comp = comp); colnames(df0)[1:length(cn_opt)] <- cn_opt
    a_orig <- tryCatch(auc_pos(as.numeric(predict(fb, newdata = df0, type = "response"))),
                       error = function(e) NA_real_)
    opt_auc[b] <- a_boot - a_orig
  }
  opt_m    <- mean(opt_auc, na.rm = TRUE)
  auc_corr <- auc$combined[["apparent"]] - opt_m

  # ── Calibration + decision curve (combined vs clinical, on LOO probs) ────────
  cp_comb <- calib_params(yb, pl_comb)
  prev <- mean(yb); pt_grid <- seq(0.01, 0.60, by = 0.01)
  nb <- function(p) sapply(pt_grid, function(pt) {
    pos <- p >= pt; w <- pt / (1 - pt)
    sum(pos & yb == 1) / n - sum(pos & yb == 0) / n * w })
  # PD-L1 net benefit on its complete cases (own denominator; see Fig 2 caveat)
  nb_cc <- function(p) { idx <- which(is.finite(p)); m <- length(idx)
    if (m == 0) return(rep(NA_real_, length(pt_grid)))
    sapply(pt_grid, function(pt) { pos <- p[idx] >= pt; w <- pt / (1 - pt)
      sum(pos & yb[idx] == 1) / m - sum(pos & yb[idx] == 0) / m * w }) }
  dc <- data.frame(threshold = pt_grid,
                   clinical   = nb(pl_clin), combined = nb(pl_comb),
                   pdl1       = nb_cc(pl_bench),
                   treat_all  = prev - (1 - prev) * (pt_grid / (1 - pt_grid)),
                   treat_none = 0)
  if (has_comparator) dc$comparator <- nb(pl_comp)   # display-only comparator (full clinical + immune) net benefit

  # ── Survival secondary (OS/PFS): Cox clinical vs clinical+immune ─────────────
  survival <- NULL
  if (requireNamespace("survival", quietly = TRUE) &&
      all(surv_cols %in% names(raw))) {
    os  <- suppressWarnings(as.numeric(raw[[surv_cols[["os"]]]][ridx]))[keep]
    pfs <- suppressWarnings(as.numeric(raw[[surv_cols[["pfs"]]]][ridx]))[keep]
    dth <- suppressWarnings(as.integer(raw[[surv_cols[["death"]]]][ridx]))[keep]
    ev_os  <- dth
    ev_pfs <- as.integer(dth == 1 | (dth == 0 & pfs < os))
    Cimp <- impute_median(Cmat, all_rows); cn <- make.names(cl_names, unique = TRUE)
    cox_block <- function(time, event) {
      ok <- is.finite(time) & time > 0 & !is.na(event)
      # `time > 0` silently drops patients recorded at time 0 — some of whom DIED (a real event).
      # Report the exclusion instead of hiding it.
      n_drop_zero <- sum(is.finite(time) & time <= 0 & !is.na(event))
      n_drop_ev   <- sum(is.finite(time) & time <= 0 & !is.na(event) & event == 1)
      if (sum(ok) < 20 || length(unique(event[ok])) < 2) return(NULL)
      d <- data.frame(time = time[ok], event = event[ok], Cimp[ok, , drop = FALSE], comp = comp[ok])
      colnames(d)[3:(2 + length(cn))] <- cn
      fc <- tryCatch(survival::coxph(reformulate(cn, "survival::Surv(time, event)"), data = d), error = function(e) NULL)
      fk <- tryCatch(survival::coxph(reformulate(c(cn, "comp"), "survival::Surv(time, event)"), data = d), error = function(e) NULL)
      if (is.null(fc) || is.null(fk)) return(NULL)
      lrt <- tryCatch(anova(fc, fk, test = "Chisq")[["Pr(>|Chi|)"]][2], error = function(e) NA_real_)
      # Effect size for the immune term + a PH check. Previously computed and DISCARDED: a C-index
      # and an LRT p cannot tell a reader whether a null is "no effect" or "underpowered". The HR
      # CI can, and `min detectable HR` below turns the null into an informative result.
      sk <- summary(fk)
      hr <- tryCatch(unname(sk$conf.int["comp", c(1, 3, 4)]), error = function(e) rep(NA_real_, 3))
      hp <- tryCatch(unname(sk$coefficients["comp", ncol(sk$coefficients)]), error = function(e) NA_real_)
      zp <- tryCatch(survival::cox.zph(fk)$table["comp", "p"], error = function(e) NA_real_)
      # smallest per-unit HR detectable at 80% power with the observed events (Schoenfeld)
      v  <- stats::var(comp[ok]); ev_n <- sum(event[ok])
      mdh <- if (is.finite(v) && v > 0 && ev_n > 0)
               exp((stats::qnorm(0.975) + stats::qnorm(0.80)) / sqrt(ev_n * v)) else NA_real_
      list(n = sum(ok), events = ev_n,
           c_clinical = summary(fc)$concordance[["C"]],
           c_combined = sk$concordance[["C"]],
           lrt_p = lrt,
           immune_hr = hr[1], immune_hr_ci = hr[2:3], immune_p = hp,
           ph_zph_p = zp, min_detectable_hr_80 = mdh,
           n_excluded_time_zero = n_drop_zero, n_excluded_time_zero_events = n_drop_ev)
    }
    survival <- list(OS = cox_block(os, ev_os), PFS = cox_block(pfs, ev_pfs))
    # PFS here is a PROXY, not an observed PFS: this cohort has NO progression-date column, so the
    # event is derived as death | (alive & PFS<OS). It is also near-tautological with best response
    # (BR=PD => progression at the first scan), so it cannot independently validate a response
    # model. OS is the defensible survival endpoint. Flag it in the payload so no downstream
    # consumer can quietly treat it as a real PFS.
    if (!is.null(survival$PFS))
      survival$PFS$endpoint_caveat <- paste(
        "DERIVED PROXY: no progression date exists in this cohort; event = death |",
        "(alive & PFS < OS). Near-tautological with best response (BR=PD => progression at",
        "first scan). Descriptive only — OS is the primary survival endpoint.")
  }

  # ── Sensitivity layers (additive; placed AFTER all pre-existing RNG consumers so the
  #    primary numbers — optimism etc. — reproduce byte-for-byte). Shared helpers: ───────
  idi_y <- function(po, pn, yv) { e <- yv == 1
    (mean(pn[e]) - mean(po[e])) + (mean(po[!e]) - mean(pn[!e])) }
  # anova.glm names the p column "Pr(>Chi)" (coxph uses "Pr(>|Chi|)") — match robustly.
  lrt_chisq_p <- function(a) {
    if (is.null(a)) return(NA_real_)
    pc <- grep("^Pr\\(", colnames(a), value = TRUE)
    if (!length(pc)) return(NA_real_)
    v <- a[[pc[1]]][2]; if (length(v) != 1) NA_real_ else v
  }

  # ── C1: leakage-free LOO IDI ────────────────────────────────────────────────
  # idi_app is in-sample (apparent, Pencina 2008 setup) and optimistic; the honest
  # reclassification metric is the IDI on the cross-validated LOO probs (matches diag_23).
  idi_loo <- idi_y(pl_clin, pl_comb, yb)
  idi_loo_bs <- numeric(n_boot)
  for (b in seq_len(n_boot)) { ix <- sample(n, n, replace = TRUE)
    idi_loo_bs[b] <- idi_y(pl_clin[ix], pl_comb[ix], yb[ix]) }
  idi_loo_ci <- qci(idi_loo_bs)

  # ── 1b: ridge-penalized combined model (calibration sensitivity), LOO leakage-free ──
  # The unpenalized combined model over-fits at EPV~9 (calib slope ~0.66, predictions too
  # extreme). A ridge (L2, inner-CV lambda, leakage-free LOO) recalibrates (slope -> ~0.9)
  # WITHOUT losing discrimination — disarms "the IDI is inflated by mis-calibration". The
  # unpenalized model stays PRIMARY (it carries the LRT); this is a sensitivity. See diag_23.
  combined_ridge <- NULL
  if (isTRUE(ridge) && requireNamespace("glmnet", quietly = TRUE)) {
    pr <- rep(NA_real_, n)
    for (i in all_rows) {
      tr  <- all_rows[-i]
      Cit <- impute_median(Cmat, tr)
      Xtr <- as.matrix(cbind(Cit[tr, , drop = FALSE], comp = comp[tr]))
      Xte <- as.matrix(cbind(Cit[i,  , drop = FALSE], comp = comp[i]))
      cvf <- tryCatch(glmnet::cv.glmnet(Xtr, yb[tr], family = "binomial", alpha = 0, nfolds = 10),
                      error = function(e) NULL)
      if (!is.null(cvf)) pr[i] <- as.numeric(predict(cvf, newx = Xte, s = "lambda.min", type = "response"))
    }
    ok <- is.finite(pr)
    if (sum(ok) >= min_n) {
      cpr <- calib_params(yb[ok], pr[ok]); wok <- which(ok)
      ir  <- idi_y(pl_clin[ok], pr[ok], yb[ok])
      irb <- numeric(n_boot)
      for (b in seq_len(n_boot)) { ix <- sample(wok, length(wok), replace = TRUE)
        irb[b] <- idi_y(pl_clin[ix], pr[ix], yb[ix]) }
      combined_ridge <- list(
        method           = "ridge L2 (alpha=0); inner 10-fold cv.glmnet lambda.min; LOO outer",
        n_used           = sum(ok),
        auc_clinical_loo = auc_pos(pl_clin[ok], y[ok]),
        auc_combined_loo = auc_pos(pr[ok], y[ok]),
        delta_auc_loo    = unname(auc_pos(pr[ok], y[ok]) - auc_pos(pl_clin[ok], y[ok])),
        calib_slope_loo  = unname(cpr[["slope"]]), calib_citl_loo = unname(cpr[["intercept"]]),
        brier_loo        = brier(yb[ok], pr[ok]),
        idi_loo          = ir, idi_loo_ci = qci(irb))
    }
  }

  # ── 1c: increment vs progressively enriched clinical baselines (immune held fixed) ──
  # Pre-empts "your SoC baseline is weak / a proper clinical model would absorb the immune
  # signal". Each enriched baseline re-tests the SAME composite (LRT df=1, leakage-free LOO
  # AUCs, apparent IDI). See diag_24 (trap: improving DeLong/dAUC = baseline degradation).
  baseline_sensitivity <- NULL
  if (length(baseline_sensitivity_cols)) {
    inc_for <- function(Cm2) {
      cl2 <- colnames(Cm2); cn2 <- make.names(cl2, unique = TRUE)
      plc <- plk <- numeric(n)
      for (i in all_rows) {
        tr  <- all_rows[-i]; Cit <- impute_median(Cm2, tr)
        d <- data.frame(yb = yb, Cit, comp = comp, check.names = TRUE)
        colnames(d)[2:(1 + length(cn2))] <- cn2
        fcl <- suppressWarnings(glm(reformulate(cn2, "yb"),            data = d[tr, ], family = binomial()))
        fco <- suppressWarnings(glm(reformulate(c(cn2, "comp"), "yb"), data = d[tr, ], family = binomial()))
        plc[i] <- as.numeric(predict(fcl, newdata = d[i, ], type = "response"))
        plk[i] <- as.numeric(predict(fco, newdata = d[i, ], type = "response"))
      }
      Cia <- impute_median(Cm2, all_rows)
      da  <- data.frame(yb = yb, Cia, comp = comp, check.names = TRUE)
      colnames(da)[2:(1 + length(cn2))] <- cn2
      fca <- suppressWarnings(glm(reformulate(cn2, "yb"),            data = da, family = binomial()))
      fka <- suppressWarnings(glm(reformulate(c(cn2, "comp"), "yb"), data = da, family = binomial()))
      lp  <- lrt_chisq_p(tryCatch(anova(fca, fka, test = "Chisq"), error = function(e) NULL))
      fp  <- if (requireNamespace("logistf", quietly = TRUE))
        tryCatch({ ff <- logistf::logistf(reformulate(c(cn2, "comp"), "yb"), data = da)
                   unname(ff$prob[match("comp", names(coef(ff)))]) }, error = function(e) NA_real_) else NA_real_
      ia  <- idi_fun(as.numeric(predict(fca, type = "response")), as.numeric(predict(fka, type = "response")))
      list(k_clinical = length(cl2),
           clinical_loo_auc = auc_pos(plc), combined_loo_auc = auc_pos(plk),
           delta_auc_loo = unname(auc_pos(plk) - auc_pos(plc)),
           lrt_p = lp, lrt_firth_p = fp, idi = ia)
    }
    baseline_sensitivity <- list()
    for (lab in names(baseline_sensitivity_cols)) {
      extra <- intersect(baseline_sensitivity_cols[[lab]], names(raw))
      if (!length(extra)) { message(sprintf("[ML][AV] baseline '%s': extra columns absent — skipped.", lab)); next }
      ec  <- sapply(extra, function(col) suppressWarnings(as.numeric(raw[[col]][ridx])))
      ec  <- matrix(ec, nrow = length(meta_ids), dimnames = list(NULL, extra))[keep, , drop = FALSE]
      Cm2 <- cbind(Cmat, ec)
      baseline_sensitivity[[lab]] <- tryCatch(inc_for(Cm2), error = function(e) {
        message(sprintf("[ML][AV] baseline '%s' failed: %s", lab, e$message)); NULL })
    }
    if (!length(baseline_sensitivity)) baseline_sensitivity <- NULL
  }

  # ── 1d: specificity null — gate vs random k-marker immune composites over the FIXED SoC ──
  # "With 39 markers, isn't a composite that beats a weak baseline easy to find by chance?"
  # Random triplets scored for the SAME apparent increment; specificity p = how extreme the
  # gate is vs the null. Claim rests on PROVENANCE (LMM + published biology), not rarity (diag_25).
  specificity_null <- NULL
  if (n_random_null > 0L) {
    zk <- as.matrix(z[keep, , drop = FALSE])
    zk <- apply(zk, 2, function(x) suppressWarnings(as.numeric(x)))
    keep_col <- apply(zk, 2, function(x) sum(is.finite(x)) > 0 &&
                        stats::sd(x, na.rm = TRUE) > 0)
    zk <- zk[, keep_col, drop = FALSE]
    for (j in seq_len(ncol(zk))) { col <- zk[, j]; col[is.na(col)] <- median(col, na.rm = TRUE); zk[, j] <- col }
    k_gate <- length(avail)
    if (ncol(zk) >= k_gate) {
      Ci_all <- impute_median(Cmat, all_rows); cnf <- make.names(cl_names, unique = TRUE)
      base_df <- data.frame(yb = yb, Ci_all); colnames(base_df)[2:(1 + length(cnf))] <- cnf
      fc0 <- suppressWarnings(glm(reformulate(cnf, "yb"), data = base_df, family = binomial()))
      pc0 <- as.numeric(predict(fc0, type = "response")); auc_c_app <- auc_pos(pc0)
      gate_dauc_app <- unname(auc$combined[["apparent"]] - auc$clinical[["apparent"]])
      rl <- ri <- rd <- rep(NA_real_, n_random_null)
      for (b in seq_len(n_random_null)) {
        mk <- sample(ncol(zk), k_gate)
        d2 <- base_df; d2$comp <- rowMeans(zk[, mk, drop = FALSE])
        fk <- suppressWarnings(glm(reformulate(c(cnf, "comp"), "yb"), data = d2, family = binomial()))
        rl[b] <- lrt_chisq_p(tryCatch(anova(fc0, fk, test = "Chisq"), error = function(e) NULL))
        pk <- as.numeric(predict(fk, type = "response"))
        ri[b] <- idi_fun(pc0, pk); rd[b] <- auc_pos(pk) - auc_c_app
      }
      fin <- is.finite(rl)
      specificity_null <- list(
        n_random = n_random_null, k_markers = k_gate, n_pool = ncol(zk),
        gate_lrt_p = lrt_p, gate_idi = idi_app, gate_delta_auc_apparent = gate_dauc_app,
        spec_p_lrt  = (1 + sum(rl[fin] <= lrt_p))           / (sum(fin) + 1),
        spec_p_idi  = (1 + sum(ri[fin] >= idi_app))         / (sum(fin) + 1),
        spec_p_dauc = (1 + sum(rd[fin] >= gate_dauc_app))   / (sum(fin) + 1),
        frac_sig_lrt = mean(rl[fin] < 0.05),
        null_dauc_q  = as.list(round(quantile(rd[fin], c(.5, .9, .95, .99), na.rm = TRUE), 4)),
        note = paste("Apparent LRT/IDI/dAUC vs random k-marker composites; specificity rests on",
                     "provenance (Step-04 LMM + published Ki67 biology), not numerical rarity (diag_25)."))
    }
  }

  # ── Dynamics <-> baseline coupling (Blocco G; gate rationale / Limitations) ──
  # WHY the dynamically-selected gate predicts at baseline: across the marker panel,
  # is LMM Timepoint x Group interaction strength correlated with each marker's
  # standalone baseline (T0) discrimination? A positive correlation = convergent
  # validity (regression-to-the-mean: responders start higher AND contract more) and
  # also bounds the residual selection optimism we disclose. Fully DETERMINISTIC
  # (univariate AUC + analytic correlation; no sampling) -> consumes no RNG, so the
  # primary numbers reproduce byte-for-byte regardless of where this block sits.
  dynamics_baseline_coupling <- NULL
  if (!is.null(lmm_interaction) && is.data.frame(lmm_interaction) &&
      all(c("Marker", "T_Value_Interaction") %in% names(lmm_interaction))) {
    zp <- as.matrix(z[keep, , drop = FALSE])
    zp <- apply(zp, 2, function(x) suppressWarnings(as.numeric(x)))
    okc <- apply(zp, 2, function(x) sum(is.finite(x)) > 0 && stats::sd(x, na.rm = TRUE) > 0)
    zp <- zp[, okc, drop = FALSE]
    for (j in seq_len(ncol(zp))) { col <- zp[, j]; col[is.na(col)] <- median(col, na.rm = TRUE); zp[, j] <- col }
    # standalone T0 discrimination per panel marker (oriented univariate AUC, deterministic)
    t0_auc <- vapply(seq_len(ncol(zp)), function(j) auc_pos(zp[, j]), numeric(1))
    pm <- data.frame(Marker = colnames(zp), T0_AUC = t0_auc, stringsAsFactors = FALSE)
    li_cols <- intersect(c("Marker", "Estimate_Interaction", "T_Value_Interaction",
                           "P_Value_Interaction", "FDR_Interaction"), names(lmm_interaction))
    cp <- merge(pm, lmm_interaction[, li_cols, drop = FALSE], by = "Marker")
    cp$absT    <- abs(cp$T_Value_Interaction)
    cp$is_gate <- cp$Marker %in% avail
    cp <- cp[is.finite(cp$absT) & is.finite(cp$T0_AUC), , drop = FALSE]
    if (nrow(cp) >= 8L) {
      pe   <- suppressWarnings(stats::cor.test(cp$absT, cp$T0_AUC, method = "pearson"))
      sp   <- suppressWarnings(stats::cor.test(cp$absT, cp$T0_AUC, method = "spearman"))
      ng   <- cp[!cp$is_gate, , drop = FALSE]
      peng <- if (nrow(ng) >= 8L) suppressWarnings(stats::cor.test(ng$absT, ng$T0_AUC, method = "pearson")) else NULL
      spng <- if (nrow(ng) >= 8L) suppressWarnings(stats::cor.test(ng$absT, ng$T0_AUC, method = "spearman")) else NULL
      cp$rank_T0  <- rank(-cp$T0_AUC, ties.method = "min")
      cp$rank_int <- rank(-cp$absT,  ties.method = "min")
      top_t0 <- cp$Marker[order(cp$rank_T0)][1:min(3L, nrow(cp))]
      dynamics_baseline_coupling <- list(
        metric   = "abs_interaction_t_vs_univariate_T0_AUC",
        n_panel  = nrow(cp), n_gate = sum(cp$is_gate),
        pearson_r  = unname(pe$estimate), pearson_ci = unname(pe$conf.int), pearson_p = pe$p.value,
        spearman_rho = unname(sp$estimate), spearman_p = sp$p.value,
        pearson_r_excl_gate  = if (!is.null(peng)) unname(peng$estimate) else NA_real_,
        pearson_p_excl_gate  = if (!is.null(peng)) peng$p.value else NA_real_,
        spearman_rho_excl_gate = if (!is.null(spng)) unname(spng$estimate) else NA_real_,
        gate_rank_t0  = sort(cp$rank_T0[cp$is_gate]),
        gate_rank_int = sort(cp$rank_int[cp$is_gate]),
        top_t0_markers = top_t0,
        gate_is_top_t0_marker = isTRUE(cp$is_gate[order(cp$rank_T0)][1]),
        per_marker = cp[order(cp$rank_int),
                        c("Marker", "is_gate", "Estimate_Interaction", "P_Value_Interaction",
                          "FDR_Interaction", "absT", "T0_AUC", "rank_int", "rank_T0")],
        note = paste("Panel-level convergent validity: LMM Timepoint x Group interaction strength",
                     "vs standalone baseline (T0) AUC. A positive correlation explains why the",
                     "dynamically-selected gate also predicts at baseline (regression-to-the-mean)",
                     "and bounds the residual selection optimism. The gate is NOT the top-T0 set,",
                     "so selection is dynamics-driven, not a baseline scan. Mechanism /",
                     "Limitations only; NOT proof of gate optimality (Blocco G)."))
    }
  }

  provenance_note <- if (gate_provenance == "lmm-prespecified")
    paste("Immune composite pre-specified from the Step-04 LMM gate (disjoint from this T0",
          "model). Clinical baseline pre-specified from NSCLC-ICI literature; not data-selected.")
  else
    paste("Immune composite selected on these T0 data (univariate/screen) — added value is",
          "EXPLORATORY / conditional-on-selection; treat IDI and DeLong p as a lower bound.")

  # ── Locked model: a frozen, pre-specified scorer for prospective/temporal validation ──
  # Bundles everything needed to score NEW patients WITHOUT refitting: the apparent
  # combined-model coefficients (developed on this full cohort) PLUS the development-cohort
  # preprocessing constants, so a new patient's RAW markers map to a probability through the
  # exact same transform. The immune side reverses Step 01's pipeline: logit/log2 (frozen
  # epsilon, attached by the caller) -> z with FROZEN per-marker center/scale (recovered from
  # `hybrid_data_raw`, the pre-z transformed matrix) -> median-impute on the z scale. The
  # clinical side stores per-variable medians (raw scale). See diagnostics/run_validation.R.
  #
  # Constants MUST come from the same rows Step 01's scale() saw, because it discards
  # `scaled:center`/`scaled:scale` (modules_coda.R:276-277) and recomputing from raw is the
  # only recovery route. On a T1/delta timepoint `DATA_T0` was REPLACED above by
  # build_timepoint_composite_data() (no `hybrid_data_raw`) and the z's came from the
  # POOLED T0+T1 longitudinal matrix — so read the constants from the full, unrestricted
  # DATA_LONG. Using the T0 RDS's constants here would be a silent wrong answer, and
  # reading the destroyed DATA_T0 yields NaN/NA (a lock that scores every patient to a
  # constant instead of erroring).
  raw_src  <- if (identical(composite_timepoint, "T0")) DATA_T0 else DATA_LONG
  raw_df   <- as.data.frame(raw_src$hybrid_data_raw)
  z_center <- sapply(avail, function(m) mean(suppressWarnings(as.numeric(raw_df[[m]])),     na.rm = TRUE))
  z_scale  <- sapply(avail, function(m) stats::sd(suppressWarnings(as.numeric(raw_df[[m]])), na.rm = TRUE))
  z_imp    <- sapply(avail, function(m) median(suppressWarnings(as.numeric(z[[m]])),         na.rm = TRUE))
  # A non-scorable lock must never ship silently — fail at the source instead.
  bad_const <- !is.finite(z_center) | !is.finite(z_scale) | z_scale <= 0
  if (any(bad_const))
    stop(sprintf(paste("[ML][AV] locked model: non-finite z constants for {%s} (timepoint=%s).",
                       "The frozen scorer would be unusable — check %s$hybrid_data_raw."),
                 paste(avail[bad_const], collapse = ", "), composite_timepoint,
                 if (identical(composite_timepoint, "T0")) "DATA_T0" else "DATA_LONG"))
  clin_med <- apply(Cmat, 2, median, na.rm = TRUE); clin_med[!is.finite(clin_med)] <- 0
  locked_model <- list(
    schema_version    = 1L,
    positive_label    = resp_label, negative_label = neg_label,
    prevalence        = prev, threshold_default = round(prev, 4),
    gate_provenance   = gate_provenance, provenance_note = provenance_note,
    # Predictor timepoint the composite is read at. "delta"/"T1" => the scorer needs
    # PAIRED T0+T1 draws (z at each timepoint, then differenced); "T0" => a single draw.
    composite_timepoint = composite_timepoint,
    immune = list(
      markers  = avail,
      z_center = as.list(z_center), z_scale = as.list(z_scale), z_impute = as.list(z_imp)),
    clinical = list(
      vars          = as.list(clinical_vars),              # label -> Excel column
      coef_names    = as.list(setNames(cn_app, cl_names)), # label -> coefficient name (make.names)
      impute_median = as.list(clin_med)),                  # label -> median (raw clinical scale)
    coefficients = list(
      combined      = as.list(coef(fk_app)),               # (Intercept) + clinical make.names + comp
      clinical_only = as.list(coef(fc_app))),
    calibration_loo  = list(citl = unname(cp_comb[["intercept"]]), slope = unname(cp_comb[["slope"]])),
    development_auc  = list(clinical_loo      = auc$clinical[["loo"]],
                           combined_loo      = auc$combined[["loo"]],
                           combined_apparent = auc$combined[["apparent"]]))

  message(sprintf("[ML][AV] AUC clinical LOO=%.3f -> combined LOO=%.3f (DeLong p=%.4f; LRT p=%.4f, perm=%.4f) | IDI app=%+.3f [%.3f,%.3f] LOO=%+.3f [%.3f,%.3f] | optimism=%.3f corr=%.3f",
                  auc$clinical[["loo"]], auc$combined[["loo"]], delong_p, lrt_p, lrt_perm_p,
                  idi_app, qci(idi_bs)[1], qci(idi_bs)[2],
                  idi_loo, idi_loo_ci[1], idi_loo_ci[2], opt_m, auc_corr))
  if (!is.null(combined_ridge))
    message(sprintf("[ML][AV] ridge sensitivity: combined LOO AUC=%.3f, calib slope %.3f->%.3f, IDI-LOO=%+.3f [%.3f,%.3f]",
                    combined_ridge$auc_combined_loo, cp_comb[["slope"]], combined_ridge$calib_slope_loo,
                    combined_ridge$idi_loo, combined_ridge$idi_loo_ci[1], combined_ridge$idi_loo_ci[2]))
  if (!is.null(specificity_null))
    message(sprintf("[ML][AV] specificity null (N=%d): spec-p LRT=%.3f IDI=%.3f dAUC=%.3f | frac LRT<.05=%.3f",
                    specificity_null$n_random, specificity_null$spec_p_lrt, specificity_null$spec_p_idi,
                    specificity_null$spec_p_dauc, specificity_null$frac_sig_lrt))
  if (!is.null(dynamics_baseline_coupling)) {
    dbc <- dynamics_baseline_coupling
    message(sprintf("[ML][AV] dynamics<->baseline coupling (n=%d panel): Pearson r=%+.3f [%.3f,%.3f] p=%.4f (Spearman rho=%+.3f) | excl-gate r=%+.3f | gate T0-AUC rank %s/%d | gate-is-topT0=%s",
                    dbc$n_panel, dbc$pearson_r, dbc$pearson_ci[1], dbc$pearson_ci[2], dbc$pearson_p,
                    dbc$spearman_rho, dbc$pearson_r_excl_gate,
                    paste(dbc$gate_rank_t0, collapse = ","), dbc$n_panel, dbc$gate_is_top_t0_marker))
  }

  # ── Publication nomograms ───────────────────────────────────────────────────────
  #   A) the LMM gate markers as INDIVIDUAL predictors — a display-only per-marker read
  #      (apparent fit; richer than the 1-df composite, so more apparent optimism).
  #   B) the FORMAL model itself: the 1-df immune composite + this run's own clinical
  #      vars. Same predictors + same data as `fk_app`, so B reproduces the formal
  #      apparent fit exactly and its validated discrimination is the layer's LOO AUC
  #      (not an inflated apparent number). `include_nlr` therefore drives B per run.
  #
  #   Immune axes are RELABELLED onto their exact clinical scale (nomo_tick_label):
  #   Δ -> fold change in the positive:negative cell ratio, exp(value*sigma), with the
  #   composite at sigma=1; T0 -> the cell fraction (%) itself. Both are bijections of
  #   the model scale, so this is display-only: no refit, no approximation, tick
  #   positions and every fitted number unchanged. The T0 composite has no such reading
  #   (it averages standardized values across markers with different mu) -> left as an
  #   index. RNG-free, additive → primary numbers unaffected.
  nomo <- tryCatch({
    Zc_keep  <- as.data.frame(Zc[keep, , drop = FALSE])            # gate markers (z, imputed)
    clin_sel <- if (is.null(nomogram_clinical_vars)) cl_names
                else intersect(as.character(unlist(nomogram_clinical_vars)), cl_names)
    is_delta <- !identical(composite_timepoint, "T0")
    # Exact per-marker relabel constants — the SAME frozen z constants the locked model
    # ships (pooled longitudinal on Δ/T1, T0 RDS on T0), so figure and scorer agree.
    imm_conv <- setNames(lapply(colnames(Zc_keep), function(m) if (is_delta)
                           list(type = "fc",  sigma = unname(z_scale[[m]]))
                         else
                           list(type = "pct", mu = unname(z_center[[m]]),
                                              sigma = unname(z_scale[[m]]))),
                         colnames(Zc_keep))
    spec_imm <- build_nomogram_spec(
      yb, Zc_keep, resp_label = resp_label,
      var_scale   = setNames(rep(if (is_delta) "fc" else "pct", ncol(Zc_keep)), colnames(Zc_keep)),
      var_convert = imm_conv)
    spec_ci <- NULL
    if (length(clin_sel)) {
      # B == the formal model: composite + formal clinical, apparent-imputed (== fk_app).
      d_ci <- cbind(setNames(data.frame(comp), COMPOSITE_DISPLAY),
                    as.data.frame(Ci_app[, clin_sel, drop = FALSE]))
      spec_ci <- build_nomogram_spec(
        yb, d_ci, resp_label = resp_label,
        var_scale   = c(setNames(if (is_delta) "fc" else "z", COMPOSITE_DISPLAY),
                        setNames(rep("raw", length(clin_sel)), clin_sel)),
        # sigma = 1: the per-marker sigmas are already absorbed into the composite mean,
        # so exp(composite) is the exact weighted geometric-mean ratio fold change.
        var_convert = if (is_delta) setNames(list(list(type = "fc", sigma = 1)),
                                             COMPOSITE_DISPLAY) else NULL)
    }
    list(timepoint = composite_timepoint, positive_label = resp_label,
         gate_markers = avail, clinical_vars = clin_sel,
         composite_display = COMPOSITE_DISPLAY,
         # Panel B is the formal model → report its LEAKAGE-FREE LOO AUC (already
         # computed by this layer), never the apparent one it is drawn from.
         combined_loo_auc = auc$combined[["loo"]],
         immune = spec_imm, clinical_immune = spec_ci)
  }, error = function(e) { message(sprintf("[ML][AV] nomogram spec failed (non-fatal): %s", e$message)); NULL })

  # ── Repeated stratified k-fold CV (additive; the REPORTABLE discrimination) ───────────────
  # WHY THIS EXISTS. The LOO AUCs above are artifact-prone for TIED/DISCRETE predictors. Leaving
  # patient i out shifts the fit in a direction that depends on y_i, which breaks ties among the
  # many patients sharing a value ANTI-correlated with y_i, dragging AUC below 0.5. It is not
  # anti-prediction — it is a CV pathology, and it scales with the tie fraction (diag_41):
  #     immune  0% tied pairs -> LOO 0.771 vs 10-fold 0.778   (no deficit)
  #     PS     44% tied pairs -> LOO 0.424 vs 10-fold 0.551   (-0.13)
  #     PD-L1  14% tied pairs -> LOO 0.237 vs 10-fold 0.393   (-0.16)
  # Reporting a clinical baseline "AUC 0.309" invites the fair objection that the CV is broken,
  # and makes the increment over it uninterpretable. Repeated stratified k-fold has many patients
  # per held-out fold, so no single y_i can steer the tie-break.
  # NOTE: this changes NO test. The primary test for the increment is the LRT, which uses no
  # cross-validation at all. These are discrimination DESCRIPTORS. Folds are shared between the
  # clinical and combined models within each rep, so the delta is paired.
  # RNG-isolated (state saved/restored, local seed) => every pre-existing number reproduces
  # byte-for-byte.
  cv_kfold <- tryCatch({
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      .old_seed_cv <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", .old_seed_cv, envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed + 977L)
    k    <- min(10L, max(3L, floor(min(sum(yb), sum(1 - yb)) / 2)))
    reps <- 50L
    R <- matrix(NA_real_, reps, 3, dimnames = list(NULL, c("clinical", "immune", "combined")))
    for (r in seq_len(reps)) {
      fold <- integer(n)
      for (cl in unique(yb)) {                       # stratified fold assignment
        ix <- which(yb == cl)
        fold[ix] <- sample(rep(seq_len(k), length.out = length(ix)))
      }
      pc <- pk <- pi_ <- numeric(n)
      for (f in seq_len(k)) {
        te <- which(fold == f); tr <- which(fold != f)
        if (length(unique(yb[tr])) < 2) { pc[te] <- pk[te] <- pi_[te] <- mean(yb[tr]); next }
        pc[te]  <- fit_predict(Cmat, comp, tr, te, use_immune = FALSE)
        pk[te]  <- fit_predict(Cmat, comp, tr, te, use_immune = TRUE)
        di      <- data.frame(yb = yb, comp = comp)
        fi      <- suppressWarnings(glm(yb ~ comp, data = di[tr, , drop = FALSE], family = binomial()))
        pi_[te] <- as.numeric(predict(fi, newdata = di[te, , drop = FALSE], type = "response"))
      }
      R[r, ] <- c(auc_pos(pc), auc_pos(pi_), auc_pos(pk))
    }
    dlt <- R[, "combined"] - R[, "clinical"]
    list(method = sprintf("repeated stratified %d-fold CV, %d reps (paired folds)", k, reps),
         k = k, reps = reps,
         auc_clinical = unname(mean(R[, "clinical"])),
         auc_immune   = unname(mean(R[, "immune"])),
         auc_combined = unname(mean(R[, "combined"])),
         delta_auc    = unname(mean(dlt)),
         delta_auc_ci = unname(stats::quantile(dlt, c(.025, .975))),
         baseline_degenerate = unname(mean(R[, "clinical"]) < 0.5),
         note = paste("Reportable discrimination. LOO values in `auc` are retained for",
                      "continuity but are artifact-prone for tied/discrete clinical predictors",
                      "(see diag_41); a sub-chance LOO AUC is a CV pathology, not anti-prediction.",
                      "The increment's primary test is the LRT, which uses no CV."))
  }, error = function(e) {
    message(sprintf("[ML][AV] k-fold CV block failed (non-fatal): %s", e$message)); NULL
  })

  # Panel B of the nomogram IS the formal model, so it must be labelled with the
  # REPORTABLE discrimination (repeated stratified k-fold), not the LOO AUC: the latter
  # is tie-artifact-prone for the discrete clinical vars this model contains (diag_41),
  # and is the number the project retired. `combined_loo_auc` is retained for continuity
  # (older persisted publication_data rds carry only that field); the figure prefers
  # `combined_cv_auc` and falls back to LOO only when it is absent. Display-only.
  if (!is.null(nomo) && !is.null(cv_kfold)) {
    nomo$combined_cv_auc <- cv_kfold$auc_combined
    nomo$combined_cv_k   <- cv_kfold$k
  }

  list(
    clinical_vars   = as.list(clinical_vars),
    formal_vars     = as.list(colnames(Cmat)),   # vars actually in the FORMAL baseline
    nomogram        = nomo,
    cv_kfold        = cv_kfold,
    locked_model    = locked_model,
    composite_markers = avail,
    composite_timepoint = composite_timepoint,
    gate_provenance = gate_provenance, provenance_note = provenance_note,
    n = n, n_complete_case = n_cc, n_pos = sum(yb), n_neg = sum(1 - yb),
    prevalence = prev, positive_label = resp_label, negative_label = neg_label,
    auc = auc,
    comparator = if (has_comparator) list(
      label = comparator_label, vars = as.list(colnames(Cmat_full)),
      auc = list(apparent = auc_pos(pa_comp), loo = auc_pos(pl_comp)), n = n) else NULL,
    increment = list(
      delta_auc_loo = unname(auc$combined[["loo"]] - auc$clinical[["loo"]]),
      delong_p_loo  = delong_p,
      lrt_p = lrt_p, lrt_perm_p = lrt_perm_p, lrt_firth_p = lrt_firth_p,
      # `jsonlite::write_json(digits = 4)` rounds to 4 DECIMAL PLACES, so any p < 0.00005
      # serialises as literally `0` — which is what the JSON has been reporting for the
      # headline LRT (true value 3.787e-05). "p = 0" is not a reportable number. These
      # companions survive the rounding and are exactly invertible: p = 10^-lrt_neglog10_p.
      # Additive; the raw `lrt_p` fields above are untouched.
      lrt_neglog10_p       = if (is.finite(lrt_p)      && lrt_p      > 0) -log10(lrt_p)      else NA_real_,
      lrt_firth_neglog10_p = if (is.finite(lrt_firth_p) && lrt_firth_p > 0) -log10(lrt_firth_p) else NA_real_,
      lrt_p_sci            = if (is.finite(lrt_p)) format(lrt_p, digits = 4, scientific = TRUE) else NA_character_,
      idi = idi_app, idi_ci = qci(idi_bs),
      idi_loo = idi_loo, idi_loo_ci = idi_loo_ci,
      cnri = nri_app, cnri_ci = qci(nri_bs)),
    combined_optimism = list(auc_optimism = opt_m, auc_corrected = auc_corr,
                             calib_slope_loo = cp_comb[["slope"]],
                             calib_citl_loo = cp_comb[["intercept"]],
                             brier_loo = brier(yb, pl_comb)),
    combined_ridge       = combined_ridge,
    baseline_sensitivity = baseline_sensitivity,
    specificity_null     = specificity_null,
    dynamics_baseline_coupling = dynamics_baseline_coupling,
    decision_curve = list(prevalence = prev, curve = dc),
    survival = survival,
    per_patient = data.frame(
      Patient_ID = pid, True_Group = as.character(grp),
      Composite = round(comp, 4),
      Prob_Clinical = round(pl_clin, 4), Prob_Combined = round(pl_comb, 4),
      Prob_ClinicalComp = if (!is.null(pl_comp)) round(pl_comp, 4) else NA_real_,
      Prob_PDL1 = round(pl_bench, 4), PD_L1 = round(as.numeric(Cmat[, 1]), 4),
      Prob_PDL1_APP = round(submodel_app$PDL1, 4),
      Prob_PS_APP = round(submodel_app$PS, 4),
      Prob_PDL1PS_APP = round(submodel_app$PDL1PS, 4),
      Prob_Combined_APP = round(pa_comb, 4),
      stringsAsFactors = FALSE)
  )
}

#' JSON-serializable view of an added-value result (drops per-patient + curve).
clinical_immune_json <- function(av) {
  if (is.null(av)) return(NULL)
  av$per_patient <- NULL
  av$decision_curve$curve <- NULL
  av$locked_model <- NULL          # serialized separately via write_locked_model()
  av$nomogram     <- NULL          # figure/RDS artifact only (kept out of the JSON contract)
  av
}

#' Serialize a locked added-value model (RDS + JSON) for prospective/temporal validation.
#'
#' The locked model is the pre-specified scorer built by run_clinical_immune_added_value();
#' `extra` is merged at the top level by the caller (the transform recipe from `config` plus
#' development-run metadata) so the on-disk artifact is fully self-describing. RDS is the
#' source of truth (full numeric precision); JSON is a human-readable mirror.
write_locked_model <- function(lm, rds_path, json_path = NULL, extra = list()) {
  if (is.null(lm)) return(invisible(NULL))
  lm <- modifyList(lm, extra)
  saveRDS(lm, rds_path)
  if (!is.null(json_path) && requireNamespace("jsonlite", quietly = TRUE))
    jsonlite::write_json(lm, json_path, auto_unbox = TRUE, pretty = TRUE,
                         digits = 12, null = "null")
  invisible(lm)
}

#' Flat Metric/Value summary of an added-value result, for the Excel report.
write_clinical_immune_summary <- function(av) {
  if (is.null(av)) return(data.frame(Metric = character(), Value = character()))
  inc <- av$increment; sv <- av$survival
  surv_str <- function(s) if (is.null(s)) "n/a"
    else sprintf("C %.3f->%.3f (LRT p=%.3f)", s$c_clinical, s$c_combined, s$lrt_p)
  fmt_p <- function(p) if (is.null(p) || is.na(p)) "n/a" else sprintf("%.4f", p)
  rg <- av$combined_ridge; sn <- av$specificity_null
  base_metric <- c("Clinical baseline", "Immune composite", "Gate provenance",
               "N (complete-case)", "N pos/neg", "Prevalence",
               "AUC clinical (LOO)", "AUC immune (LOO)", "AUC combined (LOO)",
               "Delta AUC (combined-clinical, LOO)", "DeLong p (LOO)",
               "LRT p (asymptotic)", "LRT p (permutation)", "LRT p (Firth penalized)",
               "IDI apparent [95% CI]", "IDI LOO leakage-free [95% CI]", "cNRI [95% CI]",
               "Combined AUC optimism", "Combined AUC corrected",
               "Combined calib slope (LOO)", "Combined Brier (LOO)",
               "Survival OS (secondary)", "Survival PFS (secondary)",
               "Provenance note")
  base_value <- c(
      paste(names(av$clinical_vars), collapse = "+"),
      paste(av$composite_markers, collapse = "+"), av$gate_provenance,
      sprintf("%d (%d)", av$n, av$n_complete_case),
      sprintf("%d/%d", av$n_pos, av$n_neg), sprintf("%.3f", av$prevalence),
      sprintf("%.3f", av$auc$clinical[["loo"]]), sprintf("%.3f", av$auc$immune[["loo"]]),
      sprintf("%.3f", av$auc$combined[["loo"]]),
      sprintf("%+.3f", inc$delta_auc_loo), sprintf("%.4f", inc$delong_p_loo),
      fmt_p(inc$lrt_p), fmt_p(inc$lrt_perm_p), fmt_p(inc$lrt_firth_p),
      sprintf("%+.3f [%.3f, %.3f]", inc$idi, inc$idi_ci[1], inc$idi_ci[2]),
      if (!is.null(inc$idi_loo)) sprintf("%+.3f [%.3f, %.3f]", inc$idi_loo, inc$idi_loo_ci[1], inc$idi_loo_ci[2]) else "n/a",
      sprintf("%+.3f [%.3f, %.3f]", inc$cnri, inc$cnri_ci[1], inc$cnri_ci[2]),
      sprintf("%.4f", av$combined_optimism$auc_optimism),
      sprintf("%.3f", av$combined_optimism$auc_corrected),
      sprintf("%.3f", av$combined_optimism$calib_slope_loo),
      sprintf("%.3f", av$combined_optimism$brier_loo),
      surv_str(sv$OS), surv_str(sv$PFS), av$provenance_note)

  # ── Ridge calibration sensitivity (1b) ──────────────────────────────────────
  if (!is.null(rg)) {
    base_metric <- c(base_metric,
      "[Ridge] Combined AUC (LOO)", "[Ridge] Delta AUC (LOO)",
      "[Ridge] calib slope (LOO)", "[Ridge] IDI LOO [95% CI]")
    base_value <- c(base_value,
      sprintf("%.3f", rg$auc_combined_loo), sprintf("%+.3f", rg$delta_auc_loo),
      sprintf("%.3f", rg$calib_slope_loo),
      sprintf("%+.3f [%.3f, %.3f]", rg$idi_loo, rg$idi_loo_ci[1], rg$idi_loo_ci[2]))
  }
  # ── Baseline-specification sensitivity (1c) ─────────────────────────────────
  if (!is.null(av$baseline_sensitivity)) {
    for (lab in names(av$baseline_sensitivity)) {
      bs <- av$baseline_sensitivity[[lab]]; if (is.null(bs)) next
      base_metric <- c(base_metric, sprintf("[Baseline %s] clin/comb LOO AUC; LRT; IDI", lab))
      base_value  <- c(base_value, sprintf("%.3f -> %.3f (dAUC %+.3f); LRT %s (Firth %s); IDI %+.3f",
        bs$clinical_loo_auc, bs$combined_loo_auc, bs$delta_auc_loo,
        fmt_p(bs$lrt_p), fmt_p(bs$lrt_firth_p), bs$idi))
    }
  }
  # ── Specificity null (1d) ───────────────────────────────────────────────────
  if (!is.null(sn)) {
    base_metric <- c(base_metric,
      "[SpecNull] N random / pool", "[SpecNull] spec-p (LRT / IDI / dAUC)",
      "[SpecNull] frac random LRT<0.05")
    base_value <- c(base_value,
      sprintf("%d / %d", sn$n_random, sn$n_pool),
      sprintf("%.3f / %.3f / %.3f", sn$spec_p_lrt, sn$spec_p_idi, sn$spec_p_dauc),
      sprintf("%.3f", sn$frac_sig_lrt))
  }
  # ── Dynamics<->baseline coupling (Blocco G) ─────────────────────────────────
  cpl <- av$dynamics_baseline_coupling
  if (!is.null(cpl)) {
    base_metric <- c(base_metric,
      "[Coupling] panel n / gate n",
      "[Coupling] Pearson r [95% CI] (p)",
      "[Coupling] Spearman rho (p)",
      "[Coupling] Pearson r excl. gate (p)",
      "[Coupling] gate T0-AUC ranks / interaction ranks",
      "[Coupling] gate is top-T0 marker")
    base_value <- c(base_value,
      sprintf("%d / %d", cpl$n_panel, cpl$n_gate),
      sprintf("%+.3f [%.3f, %.3f] (%s)", cpl$pearson_r, cpl$pearson_ci[1], cpl$pearson_ci[2], fmt_p(cpl$pearson_p)),
      sprintf("%+.3f (%s)", cpl$spearman_rho, fmt_p(cpl$spearman_p)),
      sprintf("%+.3f (%s)", cpl$pearson_r_excl_gate, fmt_p(cpl$pearson_p_excl_gate)),
      sprintf("%s / %s", paste(cpl$gate_rank_t0, collapse = ","), paste(cpl$gate_rank_int, collapse = ",")),
      as.character(cpl$gate_is_top_t0_marker))
  }
  data.frame(Metric = base_metric, Value = base_value, stringsAsFactors = FALSE)
}

# ── Nested gate-selection validation (anti-circularity / selection-optimism) ──
#' Re-runs the LMM Time×Group gate selection INSIDE each outer LOO fold and
#' recomputes the clinical vs clinical+immune out-of-fold predictions, so the
#' gate-choice variability is priced in (the pre-specified-gate added-value layer
#' is honest only CONDITIONAL on the fixed gate). Answers the #1 reviewer attack
#' ("gate selected + evaluated in one cohort") with a direct test. Promoted from
#' diagnostics/diag_30. Selection = FDR<0.05 on the interaction, NO collinearity
#' filter (matches the composite, which averages the FDR+LOO markers; VIF≈1).
#' Additive / RNG-free → never alters primary numbers. Returns NULL on failure.
#' @return list(selection_rule, n_folds, n_markers, auc=list(prespec,nested),
#'   delta=c(prespec,nested), delong_nested_p, gate_freq (named over all markers),
#'   prespec_markers, prespec_recovery, extra_markers, median_k, k_range)
run_nested_gate_validation <- function(DATA_T0, DATA_LONG, input_file,
                                       clinical_vars, prespec_gate, resp_label,
                                       fdr_thr = 0.05, seed = 2026) {
  if (!all(vapply(c("lmerTest", "lme4", "pROC", "readxl"), requireNamespace,
                  logical(1), quietly = TRUE))) {
    warning("[ML][Nested] lmerTest/lme4/pROC/readxl unavailable — skipped."); return(NULL)
  }
  set.seed(seed)
  META <- c("Patient_ID", "Sample_ID", "Timepoint", "Group")
  z_t0 <- as.data.frame(DATA_T0$hybrid_data_z)
  pid  <- DATA_T0$metadata$Patient_ID
  y    <- as.integer(DATA_T0$metadata$Group == resp_label)
  df_long <- as.data.frame(DATA_LONG$hybrid_data_z)
  markers <- intersect(setdiff(colnames(z_t0), META), setdiff(colnames(df_long), META))
  prespec <- intersect(prespec_gate, markers)
  if (length(prespec) == 0 || !"Timepoint" %in% colnames(df_long)) return(NULL)

  clin_raw <- tryCatch(as.data.frame(readxl::read_excel(input_file)), error = function(e) NULL)
  if (is.null(clin_raw) || !all(c("Patient_ID", clinical_vars) %in% colnames(clin_raw))) {
    warning("[ML][Nested] clinical columns not found — skipped."); return(NULL)
  }
  clin <- clin_raw[match(pid, clin_raw$Patient_ID), clinical_vars, drop = FALSE]
  colnames(clin) <- names(clinical_vars)
  clin[] <- lapply(clin, function(v) suppressWarnings(as.numeric(v)))

  # FDR<0.05 Time×Group selection on a training subset (no collinearity filter)
  select_gate <- function(df_sub) {
    p <- vapply(markers, function(m) {
      v <- as.numeric(df_sub[[m]]); s <- sd(v, na.rm = TRUE); if (!is.finite(s) || s == 0) s <- 1
      d <- data.frame(V = v / s, Ti = factor(df_sub$Timepoint), G = factor(df_sub$Group),
                      ID = factor(df_sub$Patient_ID)); d <- d[complete.cases(d), ]
      if (nrow(d) < 10 || length(unique(d$G)) < 2) return(NA_real_)
      tryCatch({
        mod <- suppressMessages(suppressWarnings(lmerTest::lmer(
          V ~ Ti * G + (1 | ID), data = d, REML = TRUE,
          control = lme4::lmerControl(calc.derivs = FALSE))))
        ct <- summary(mod)$coefficients; i <- grep("Ti.*:G", rownames(ct))
        if (length(i) != 1) NA_real_ else ct[i, grep("Pr\\(>\\|t\\|\\)", colnames(ct))]
      }, error = function(e) NA_real_)
    }, numeric(1))
    fdr <- p.adjust(p, method = "BH")
    markers[!is.na(fdr) & fdr < fdr_thr]
  }
  fold_predict <- function(i, gate) {
    tr <- setdiff(seq_along(pid), i)
    cl <- clin
    for (cc in names(cl)) cl[is.na(cl[, cc]), cc] <- median(cl[tr, cc], na.rm = TRUE)
    dC <- data.frame(y = y, cl)
    fitC <- suppressWarnings(glm(y ~ ., data = dC, family = binomial(), subset = tr))
    pC <- as.numeric(predict(fitC, dC[i, , drop = FALSE], type = "response"))
    if (length(gate) == 0) return(c(clin = pC, comb = pC, k = 0))
    dK <- data.frame(y = y, cl, composite = rowMeans(z_t0[, gate, drop = FALSE]))
    fitK <- suppressWarnings(glm(y ~ ., data = dK, family = binomial(), subset = tr))
    c(clin = pC, comb = as.numeric(predict(fitK, dK[i, , drop = FALSE], type = "response")),
      k = length(gate))
  }
  auc1 <- function(p) as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE, direction = "<")))

  pre <- t(vapply(seq_along(pid), function(i) fold_predict(i, prespec), numeric(3)))
  message(sprintf("[ML][Nested] re-selecting gate inside %d LOO folds...", length(pid)))
  sel_count <- setNames(integer(length(markers)), markers); kvec <- integer(length(pid))
  nes <- matrix(NA_real_, length(pid), 2, dimnames = list(NULL, c("clin", "comb")))
  recov <- 0L
  for (i in seq_along(pid)) {
    gate <- select_gate(df_long[df_long$Patient_ID %in% pid[setdiff(seq_along(pid), i)], ])
    sel_count[gate] <- sel_count[gate] + 1L
    if (all(prespec %in% gate)) recov <- recov + 1L
    pr <- fold_predict(i, gate); nes[i, ] <- pr[c("clin", "comb")]; kvec[i] <- pr["k"]
  }
  dl <- tryCatch(pROC::roc.test(pROC::roc(y, nes[, "comb"], quiet = TRUE, direction = "<"),
                                pROC::roc(y, nes[, "clin"], quiet = TRUE, direction = "<"),
                                method = "delong")$p.value, error = function(e) NA_real_)
  extra <- setdiff(names(sel_count)[sel_count == length(pid)], prespec)
  list(
    selection_rule  = "FDR<0.05 Time×Group interaction, re-selected per outer LOO fold; no collinearity filter (matches the averaging composite)",
    n_folds = length(pid), n_markers = length(markers),
    auc = list(prespec = c(clinical = auc1(pre[, "clin"]), combined = auc1(pre[, "comb"])),
               nested  = c(clinical = auc1(nes[, "clin"]), combined = auc1(nes[, "comb"]))),
    delta = c(prespec = auc1(pre[, "comb"]) - auc1(pre[, "clin"]),
              nested  = auc1(nes[, "comb"]) - auc1(nes[, "clin"])),
    delong_nested_p = dl,
    gate_freq = sel_count,
    prespec_markers = prespec,
    prespec_recovery = recov / length(pid),
    extra_markers = extra,
    median_k = stats::median(kvec), k_range = range(kvec))
}

# Compact JSON node for the nested-validation result.
nested_validation_json <- function(nv) {
  if (is.null(nv)) return(NULL)
  list(selection_rule = nv$selection_rule, n_folds = nv$n_folds, n_markers = nv$n_markers,
       auc = nv$auc, delta = as.list(nv$delta), delong_nested_p = nv$delong_nested_p,
       prespec_markers = nv$prespec_markers, prespec_recovery = nv$prespec_recovery,
       extra_markers = nv$extra_markers, median_k = nv$median_k, k_range = nv$k_range,
       gate_freq = as.list(round(nv$gate_freq / nv$n_folds, 3)))
}

