# R/modules_ml_gate.R
# ==============================================================================
# MACHINE LEARNING MODULE — GATE
# Split from the former monolithic R/modules_ml.R (sourced via that aggregator).
# Dependencies: dplyr, ggplot2, tidyr, glmnet, e1071, pROC
# ==============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)

#' @title Load LMM LOO-Robust Features
#' @description
#' Reads the Step 04 JSON and extracts markers that survive both FDR correction
#' and Leave-One-Out sensitivity testing. Returns an empty list when no markers
#' qualify, enabling graceful step-level skip without pipeline interruption.
#'
#' @param lmm_json_path String. Path to Machine_Metrics_LMM_<name>.json from Step 04.
#' @param fdr_threshold Numeric. Maximum FDR_Interaction to pass the gate (default 0.05).
#' @param loo_threshold Numeric. Maximum Max_P_Value_LOO to pass the gate (default 0.05).
#' @return A named list: markers (character), n_robust (integer), df_robust (data.frame).
#' @export
load_lmm_robust_features <- function(lmm_json_path,
                                     fdr_threshold = 0.05,
                                     loo_threshold = 0.05) {
  empty <- list(markers = character(0), n_robust = 0L, df_robust = data.frame(),
                full_results = data.frame())

  if (!file.exists(lmm_json_path)) {
    message(sprintf("   [ML] LMM JSON not found: %s. Skipping ML analysis.", lmm_json_path))
    return(empty)
  }

  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    message("   [ML] Package 'jsonlite' required to read LMM JSON. Skipping ML analysis.")
    return(empty)
  }

  metrics <- jsonlite::read_json(lmm_json_path, simplifyVector = TRUE)
  df      <- as.data.frame(metrics$full_results)

  if (nrow(df) == 0 || !"FDR_Interaction" %in% names(df)) {
    message("   [ML] LMM JSON contains no results. Skipping ML analysis.")
    return(empty)
  }

  # jsonlite drops all-NULL columns when simplifyVector = TRUE; restore as NA
  if (!"Max_P_Value_LOO" %in% names(df)) df$Max_P_Value_LOO <- NA_real_

  df_robust <- df %>%
    dplyr::filter(
      !is.na(FDR_Interaction),
      !is.na(Max_P_Value_LOO),
      FDR_Interaction < fdr_threshold,
      Max_P_Value_LOO < loo_threshold
    ) %>%
    dplyr::arrange(FDR_Interaction)

  if (nrow(df_robust) == 0) {
    message(sprintf(
      "   [ML] No markers survive LOO gate (FDR < %.2f AND LOO p < %.2f). Skipping ML analysis.",
      fdr_threshold, loo_threshold
    ))
    return(empty)
  }

  message(sprintf("   [ML] %d LOO-robust marker(s) selected for ML feature gate:", nrow(df_robust)))
  for (i in seq_len(nrow(df_robust))) {
    message(sprintf("         -> %s  (FDR=%.4f | LOO_max_p=%.4f)",
                    df_robust$Marker[i],
                    df_robust$FDR_Interaction[i],
                    df_robust$Max_P_Value_LOO[i]))
  }

  # full_results = the whole per-marker interaction table (all panel markers), used by
  # the added-value layer's dynamics<->baseline coupling analysis (Blocco G). Carried
  # alongside df_robust (the gate) so the JSON is parsed once.
  list(markers = df_robust$Marker, n_robust = nrow(df_robust), df_robust = df_robust,
       full_results = df)
}

#' @title LMM Gate Selection for a Single LOO Fold
#' @description
#' Fits per-marker LMMs on a longitudinal training subset (n-1 patients), applies
#' BH correction, and runs LOO sensitivity on FDR-passing markers to produce the
#' same two-stage gate used by Step 04, but trained on fold-local data only.
#'
#' @param df_lng_train Data.frame. Longitudinal data (meta + markers) for n-1 patients.
#' @param markers Character vector. All candidate marker names.
#' @param fdr_thr Numeric. FDR threshold (default 0.05).
#' @param loo_thr Numeric. Max LOO p-value threshold (default 0.05).
#' @return Character vector of markers that pass both gates. May be empty.
select_gate_for_fold <- function(df_lng_train, markers, fdr_thr = 0.05, loo_thr = 0.05) {
  results <- lapply(markers, function(mk) {
    vals   <- as.numeric(df_lng_train[[mk]])
    val_sd <- sd(vals, na.rm = TRUE)
    if (is.na(val_sd) || val_sd == 0) val_sd <- 1
    df_m <- data.frame(
      Value = vals / val_sd,
      Time  = as.factor(df_lng_train$Timepoint),
      Group = as.factor(df_lng_train$Group),
      ID    = as.factor(df_lng_train$Patient_ID)
    )
    df_m <- df_m[complete.cases(df_m), ]
    if (nrow(df_m) < 10 || length(unique(df_m$Group)) < 2)
      return(data.frame(Marker = mk, P = NA_real_))
    tryCatch({
      mod <- suppressMessages(suppressWarnings(
        lmerTest::lmer(Value ~ Time * Group + (1 | ID), data = df_m,
                       REML = TRUE, control = lme4::lmerControl(calc.derivs = FALSE))
      ))
      ct  <- summary(mod)$coefficients
      idx <- grep("Time.*:Group", rownames(ct))
      if (length(idx) != 1) return(data.frame(Marker = mk, P = NA_real_))
      p_col <- grep("Pr\\(>\\|t\\|\\)", colnames(ct))
      data.frame(Marker = mk, P = ct[idx, p_col])
    }, error = function(e) data.frame(Marker = mk, P = NA_real_))
  })

  res_df     <- do.call(rbind, results)
  res_df$FDR <- p.adjust(res_df$P, method = "BH")
  fdr_pass   <- res_df$Marker[!is.na(res_df$FDR) & res_df$FDR < fdr_thr]
  if (length(fdr_pass) == 0) return(character(0))

  loo_pass <- character(0)
  for (mk in fdr_pass) {
    unique_ids <- unique(df_lng_train$Patient_ID)
    max_p <- 0
    for (drop_id in unique_ids) {
      sub2  <- df_lng_train[df_lng_train$Patient_ID != drop_id, ]
      vals2 <- as.numeric(sub2[[mk]])
      sd2   <- sd(vals2, na.rm = TRUE)
      if (is.na(sd2) || sd2 == 0) sd2 <- 1
      df_m2 <- data.frame(
        Value = vals2 / sd2,
        Time  = as.factor(sub2$Timepoint),
        Group = as.factor(sub2$Group),
        ID    = as.factor(sub2$Patient_ID)
      )
      df_m2 <- df_m2[complete.cases(df_m2), ]
      if (nrow(df_m2) < 8 || length(unique(df_m2$Group)) < 2) next
      p_loo <- tryCatch({
        mod2 <- suppressMessages(suppressWarnings(
          lmerTest::lmer(Value ~ Time * Group + (1 | ID), data = df_m2,
                         REML = TRUE, control = lme4::lmerControl(calc.derivs = FALSE))
        ))
        ct2  <- summary(mod2)$coefficients
        idx2 <- grep("Time.*:Group", rownames(ct2))
        if (length(idx2) != 1) NA_real_
        else ct2[idx2, grep("Pr\\(>\\|t\\|\\)", colnames(ct2))]
      }, error = function(e) NA_real_)
      if (!is.na(p_loo) && p_loo > max_p) max_p <- p_loo
    }
    if (max_p < loo_thr) loo_pass <- c(loo_pass, mk)
  }
  loo_pass
}

#' Select Univariate Feature Gate (within a single LOO fold)
#'
#' Applies Wilcoxon rank-sum test per marker on the training fold and returns
#' markers surviving BH-FDR correction. Falls back to top-k by raw p-value if
#' fewer than 2 markers pass. Used exclusively when gate_method = "univariate".
#'
#' The FDR threshold is intentionally relaxed (default 0.20) because selection
#' operates on n-1 patients with limited power. The outer LOO permutation test
#' is the primary inferential layer.
#'
#' @param X_train Numeric matrix. Training fold (n-1 rows x p cols).
#' @param y_train Factor. Training fold binary outcome.
#' @param positive_label String. Positive class label.
#' @param fdr_threshold Numeric. BH-FDR threshold (default 0.20).
#' @param fallback_k Integer. Top-k markers used if fewer than 2 pass FDR (default 3).
#' @return Named list: markers (character vector), used_fallback (logical),
#'   pvalues (named numeric vector of raw Wilcoxon p-values for all markers).
select_univariate_gate <- function(X_train, y_train, positive_label,
                                   fdr_threshold = 0.20,
                                   fallback_k    = 3L) {
  markers <- colnames(X_train)
  y_bin   <- as.integer(y_train == positive_label)
  stopifnot(length(markers) >= 1L, sum(y_bin) >= 2L, sum(1L - y_bin) >= 2L)

  pvals <- vapply(markers, function(m) {
    x <- as.numeric(X_train[, m])
    v <- var(x, na.rm = TRUE)
    if (!is.finite(v) || v == 0) return(1.0)
    tryCatch(
      wilcox.test(x[y_bin == 1L], x[y_bin == 0L], exact = FALSE)$p.value,
      error = function(e) 1.0
    )
  }, numeric(1L))

  fdr_vals      <- p.adjust(pvals, method = "BH")
  gate          <- markers[!is.na(fdr_vals) & fdr_vals < fdr_threshold]
  used_fallback <- length(gate) < 2L

  if (used_fallback) {
    k    <- min(as.integer(fallback_k), length(markers))
    gate <- markers[order(pvals)[seq_len(k)]]
  }

  list(markers = gate, used_fallback = used_fallback, pvalues = pvals)
}

#' @title Gate Signal Decomposition (T0 / T1 / Delta)
#' @description
#' For LMM-gated longitudinal experiments, decomposes the discriminative signal
#' of gate markers across three representations: T0 (baseline), T1 (post-
#' treatment), and Delta (T1−T0). Reports univariate AUC, Cohen's d, and
#' Wilcoxon p for each, plus a DeLong comparison (T0 vs Delta), and ranks the
#' gate markers in a pure T0 Wilcoxon scan over all markers to quantify the
#' feature-selection contribution of the LMM gate.
#'
#' Addresses reviewer critique RC1 (LMM feature selection leakage),
#' RC2 (baseline-prognostic vs dynamic-predictive signal), and
#' RC3 (whether LMM guidance adds over a naive T0 scan).
#'
#' @param DATA_T0      Named list. Step 01 standard processed data object (T0).
#' @param DATA_LONG    Named list. Step 01 longitudinal processed data object.
#' @param gate_markers Character vector. LMM-gated features (post-collinearity).
#' @param resp_label   Character. Positive class label (e.g., "RP").
#' @return Named list: marker_decomp (data.frame), delong_comparisons (data.frame),
#'   gate_t0_scan_rank (data.frame), n_paired (integer), scenario (character),
#'   note (character). NULL when paired n < 10 or data unavailable.
#' @export
run_gate_signal_decomposition <- function(DATA_T0, DATA_LONG,
                                          gate_markers, resp_label) {

  META_COLS <- c("Patient_ID", "Sample_ID", "Timepoint", "Group")

  # ── Build paired T0/T1 frames ───────────────────────────────────────────────
  z_long        <- as.data.frame(DATA_LONG$hybrid_data_z)
  z_long$Patient_ID <- DATA_LONG$metadata$Patient_ID
  z_long$Timepoint  <- DATA_LONG$metadata$Timepoint
  z_long$Group      <- DATA_LONG$metadata$Group

  z_t0 <- z_long[z_long$Timepoint == "T0", ]
  z_t1 <- z_long[z_long$Timepoint == "T1", ]

  paired_ids <- intersect(z_t0$Patient_ID, z_t1$Patient_ID)
  n_paired   <- length(paired_ids)

  if (n_paired < 10) {
    message("[ML][GSD] Fewer than 10 paired patients — gate signal decomposition skipped.")
    return(NULL)
  }

  z_t0p <- z_t0[z_t0$Patient_ID %in% paired_ids, ]
  z_t0p <- z_t0p[order(z_t0p$Patient_ID), ]
  z_t1p <- z_t1[z_t1$Patient_ID %in% paired_ids, ]
  z_t1p <- z_t1p[order(z_t1p$Patient_ID), ]

  avail <- intersect(gate_markers, colnames(z_t0p))
  if (length(avail) == 0) return(NULL)

  y_bin <- as.integer(z_t0p$Group == resp_label)

  # ── Helper: Wilcoxon + AUC (DeLong CI) + Cohen's d ─────────────────────────
  decomp_one <- function(x, y_b, tp_label) {
    cc <- !is.na(x) & !is.na(y_b)
    if (sum(cc) < 6 || length(unique(y_b[cc])) < 2) {
      return(data.frame(Timepoint = tp_label, AUC = NA_real_, AUC_CI_Lo = NA_real_,
                        AUC_CI_Hi = NA_real_, Wilcox_p = NA_real_,
                        Cohen_d = NA_real_, N = sum(cc)))
    }
    x_c <- x[cc]; y_c <- y_b[cc]
    wt  <- wilcox.test(x_c[y_c == 1], x_c[y_c == 0], exact = FALSE)
    roc <- pROC::roc(y_c, x_c, quiet = TRUE, direction = "auto")
    ci  <- as.numeric(pROC::ci.auc(roc, method = "delong"))
    g1  <- x_c[y_c == 1]; g0 <- x_c[y_c == 0]
    sp  <- sqrt(((length(g1)-1)*var(g1) + (length(g0)-1)*var(g0)) /
                  (length(g1) + length(g0) - 2))
    d   <- (mean(g1) - mean(g0)) / sp
    data.frame(Timepoint = tp_label,
               AUC       = round(as.numeric(pROC::auc(roc)), 3),
               AUC_CI_Lo = round(ci[1], 3),
               AUC_CI_Hi = round(ci[3], 3),
               Wilcox_p  = round(wt$p.value, 4),
               Cohen_d   = round(d, 3),
               N         = sum(cc))
  }

  # ── Per-marker decomposition ─────────────────────────────────────────────────
  decomp_rows  <- list()
  delong_rows  <- list()

  for (m in avail) {
    x_t0 <- z_t0p[[m]]
    x_t1 <- z_t1p[[m]]
    x_d  <- x_t1 - x_t0

    r_t0 <- decomp_one(x_t0, y_bin, "T0")
    r_t1 <- decomp_one(x_t1, y_bin, "T1")
    r_d  <- decomp_one(x_d,  y_bin, "Delta (T1-T0)")

    r_t0$Marker <- r_t1$Marker <- r_d$Marker <- m
    decomp_rows <- c(decomp_rows, list(r_t0, r_t1, r_d))

    # DeLong T0 vs Delta (correlated samples — same patients)
    cc_d  <- !is.na(x_t0) & !is.na(x_d)
    p_dl  <- if (sum(cc_d) >= 10) {
      roc_t0 <- pROC::roc(y_bin[cc_d], x_t0[cc_d], quiet = TRUE, direction = "auto")
      roc_d  <- pROC::roc(y_bin[cc_d], x_d[cc_d],  quiet = TRUE, direction = "auto")
      tryCatch(pROC::roc.test(roc_t0, roc_d, method = "delong")$p.value,
               error = function(e) NA_real_)
    } else NA_real_

    delong_rows[[length(delong_rows) + 1]] <- data.frame(
      Marker                = m,
      AUC_T0                = r_t0$AUC,
      AUC_Delta             = r_d$AUC,
      DeLong_T0_vs_Delta_p  = round(p_dl, 3)
    )
  }

  marker_decomp <- do.call(rbind, decomp_rows)[
    , c("Marker", "Timepoint", "AUC", "AUC_CI_Lo", "AUC_CI_Hi",
        "Wilcox_p", "Cohen_d", "N")]
  delong_comp   <- do.call(rbind, delong_rows)
  rownames(marker_decomp) <- rownames(delong_comp) <- NULL

  # ── Pure T0 Wilcoxon scan: rank gate markers among all T0 features ──────────
  z_std  <- as.data.frame(DATA_T0$hybrid_data_z)
  y_std  <- as.integer(DATA_T0$metadata$Group == resp_label)
  all_m  <- intersect(DATA_T0$hybrid_markers, colnames(z_std))

  scan_res <- lapply(all_m, function(m) {
    x  <- z_std[[m]]
    cc <- !is.na(x)
    if (sum(cc) < 6 || length(unique(x[cc])) < 2) return(NULL)
    wt <- wilcox.test(x[cc][y_std[cc] == 1], x[cc][y_std[cc] == 0], exact = FALSE)
    data.frame(Marker = m, T0_Wilcox_p = wt$p.value)
  })
  scan_df        <- do.call(rbind, Filter(Negate(is.null), scan_res))
  scan_df        <- scan_df[order(scan_df$T0_Wilcox_p), ]
  scan_df$T0_Rank <- seq_len(nrow(scan_df))
  scan_df$N_Total <- nrow(scan_df)
  rownames(scan_df) <- NULL

  gate_rank_df <- merge(
    data.frame(Marker = avail),
    scan_df[, c("Marker", "T0_Wilcox_p", "T0_Rank", "N_Total")],
    by = "Marker", all.x = TRUE
  )
  gate_rank_df <- gate_rank_df[order(gate_rank_df$T0_Rank), ]
  rownames(gate_rank_df) <- NULL

  # ── Scenario determination ───────────────────────────────────────────────────
  t0_aucs <- marker_decomp$AUC[marker_decomp$Timepoint == "T0" & !is.na(marker_decomp$AUC)]
  d_aucs  <- marker_decomp$AUC[marker_decomp$Timepoint == "Delta (T1-T0)" &
                                  !is.na(marker_decomp$AUC)]
  scenario <- if (length(t0_aucs) > 0 && length(d_aucs) > 0) {
    md <- mean(t0_aucs, na.rm = TRUE) - mean(d_aucs, na.rm = TRUE)
    if (abs(md) < 0.05) "A/Mixed: T0 and Delta carry equivalent signal"
    else if (md > 0)    "C: T0 baseline is the primary discriminator"
    else                "B: Delta (T1-T0) is the primary discriminator"
  } else NA_character_

  message(sprintf("   [ML][GSD] n_paired=%d | Scenario: %s", n_paired, scenario))

  list(
    marker_decomp      = marker_decomp,
    delong_comparisons = delong_comp,
    gate_t0_scan_rank  = gate_rank_df,
    n_paired           = n_paired,
    scenario           = scenario,
    note               = sprintf(
      "Gate signal decomposition on %d paired patients. %s. Gate markers ranked among %d T0 markers in pure Wilcoxon scan.",
      n_paired, scenario, nrow(scan_df)
    )
  )
}

