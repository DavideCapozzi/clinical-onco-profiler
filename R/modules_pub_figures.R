# R/modules_pub_figures.R
# ==============================================================================
# PUBLICATION FIGURE BUILDERS  (manuscript figure set)
# ------------------------------------------------------------------------------
# Each builder takes a LIVE pipeline result object (clin_addval, gate_decomp,
# stratified_result, LMM bootstrap frames, ...) and returns a ggplot/patchwork
# object styled via R/modules_pub_style.R. No frozen-file reads, no recompute:
# when the analysis changes, re-running the pipeline regenerates these figures.
#
# Narrative (rebalanced 2026-06-24 toward immunology): MAIN Fig1 biology
# (gate-marker cell frequencies T0/T1×group, log-y + LMM forest) / Fig2 added-value
# (LRT-led) / Fig3 PD-L1 context (promoted) / Fig4 gate-signal decomposition
# (T0/T1/Δ, a results figure). SUPP S1 CONSORT / S4 baseline-invariance /
# S5 specificity-null / S6 standalone / S7 coupling / S8 robustness / S9 calibration
# (calibration + robustness demoted from main). Supp renumbering (gaps S2/S3) is
# deferred to the caption doc pass.
# ==============================================================================

suppressMessages({ library(ggplot2); library(patchwork) })

# ── MAIN Fig 1 — LMM Ki67 forest (pre-specified gate) ─────────────────────────
#' @param boot_df data.frame: Marker, Median_Beta_Boot, CI_Lower_2.5, CI_Upper_97.5, Pct_FDR_Significant
#' @param obs_df  data.frame: Marker, Estimate_Interaction
pub_fig_lmm_forest <- function(boot_df, obs_df = NULL, markers = NULL) {
  if (is.null(boot_df) || !nrow(boot_df)) return(NULL)
  d <- boot_df
  if (!is.null(markers)) d <- d[d$Marker %in% markers, , drop = FALSE]
  if (!nrow(d)) return(NULL)
  if (!is.null(obs_df) && all(c("Marker", "Estimate_Interaction") %in% names(obs_df))) {
    d <- merge(d, obs_df[, c("Marker", "Estimate_Interaction")], by = "Marker", all.x = TRUE)
  } else d$Estimate_Interaction <- d$Median_Beta_Boot
  d$Marker <- factor(d$Marker, levels = d$Marker[order(d$Estimate_Interaction)])
  # Single point estimate (observed β) + bootstrap 95% CI — standard forest layout.
  # The bootstrap median overlaps the observed β whenever bias ≈ 0 (the normal case),
  # so it is omitted; bootstrap robustness is conveyed by the % FDR label instead.
  ggplot(d, aes(y = Marker)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = pub_palette[["ref"]]) +
    geom_errorbarh(aes(xmin = CI_Lower_2.5, xmax = CI_Upper_97.5),
                   height = 0.22, colour = "grey35") +
    geom_point(aes(x = Estimate_Interaction), shape = 19, size = 2.6, colour = pub_palette[["combined"]]) +
    geom_text(aes(x = CI_Upper_97.5, label = sprintf(" %.0f%%", Pct_FDR_Significant)),
              hjust = 0, size = 2.5, colour = "grey30") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.18))) +
    labs(x = "Time x Group interaction (β)", y = NULL) +
    coord_cartesian(clip = "off") +
    theme_publication() + theme(plot.margin = margin(6, 12, 6, 6))
}

# ── MAIN Fig 1 (biology) — gate-marker cell frequencies (T0/T1 × group) + forest ─
#' Panel A: per-patient raw cell frequencies of the immune-gate markers at baseline
#' (T0) vs on-treatment (T1), split by response group, with within-patient spaghetti
#' lines (responders start higher and contract more). Panel B: the LMM Time×Group
#' interaction forest restricted to the same gate markers. Falls back to forest-only
#' when per-patient values are unavailable.
#' @param gate_decomp result of run_gate_signal_decomposition() (uses $per_patient_values)
#' @param boot_df,obs_df LMM bootstrap-CI and observed-interaction frames (forest)
pub_fig_biology <- function(gate_decomp, boot_df, obs_df = NULL, resp_label = NULL) {
  ppv <- if (!is.null(gate_decomp)) gate_decomp$per_patient_values else NULL
  forest <- pub_fig_lmm_forest(boot_df, obs_df,
                               markers = if (!is.null(ppv)) unique(ppv$Marker) else NULL)
  yvar <- if (!is.null(ppv) && "Value_pct" %in% names(ppv)) "Value_pct" else "Value_raw"
  if (is.null(ppv) || !nrow(ppv) || all(is.na(ppv[[yvar]]))) return(forest)

  # Responder = primary blue (consistent with the positive/combined class elsewhere),
  # non-responder = vermillion; responder column first. Fall back to alpha order.
  grps <- unique(as.character(ppv$Group))
  if (!is.null(resp_label) && resp_label %in% grps) {
    ord <- c(resp_label, setdiff(grps, resp_label))
    ppv$Group <- factor(ppv$Group, levels = ord)
  } else ord <- sort(grps)
  grp_cols <- setNames(c(pub_palette[["combined"]], pub_palette[["unpen"]])[seq_along(ord)], ord)
  # log10 y (per-facet free) — flow-cytometry-standard; decompresses the bulk so the
  # boxplots stay readable despite the heavy right-skew (rare high-proliferation
  # non-responders), without hiding any data.
  pA <- ggplot(ppv, aes(Timepoint, .data[[yvar]])) +
    geom_line(aes(group = Patient_ID), colour = "grey75", linewidth = 0.3, alpha = 0.5) +
    geom_boxplot(aes(colour = Group), fill = NA, width = 0.5,
                 outlier.shape = NA, linewidth = 0.5) +
    geom_point(aes(colour = Group), size = 0.7, alpha = 0.45) +
    facet_grid(Marker ~ Group, scales = "free_y") +
    scale_y_log10() +
    scale_colour_manual(values = grp_cols, guide = "none") +
    labs(x = NULL, y = "Cell frequency (%, log scale)") +
    theme_publication() + theme(panel.grid = element_blank())
  if (is.null(forest)) return(pub_tag(pA))
  pub_tag((pA | forest) + patchwork::plot_layout(widths = c(1.9, 1.1)))
}

# ── MAIN Fig 2 — Added value (LRT-led): ROC + decision curve ──────────────────
pub_fig_added_value <- function(av) {
  if (is.null(av) || is.null(av$per_patient)) return(NULL)
  pp  <- av$per_patient; pos <- av$positive_label
  y   <- as.integer(pp$True_Group == pos)
  rc  <- pub_roc_df(pp$Prob_Clinical, y); ro <- pub_roc_df(pp$Prob_Combined, y)
  lc  <- sprintf("Clinical (AUC %.2f)", rc$auc)
  lk  <- sprintf("Clinical + immune (AUC %.2f)", ro$auc)
  roc <- rbind(data.frame(rc$df, M = lc), data.frame(ro$df, M = lk))
  roc$M <- factor(roc$M, levels = c(lc, lk))
  inc <- av$increment
  pA <- ggplot(roc, aes(FPR, TPR, colour = M, linetype = M)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = pub_palette[["ref"]]) +
    geom_path(linewidth = 0.8) +
    annotate("text", x = 0.97, y = 0.06, hjust = 1, vjust = 0, size = 2.7,
             label = sprintf("LRT p = %.3f (perm %.3f)", inc$lrt_p,
                             if (!is.null(inc$lrt_perm_p)) inc$lrt_perm_p else NA)) +
    scale_colour_manual(values = setNames(c(pub_palette[["clinical"]], pub_palette[["combined"]]), c(lc, lk))) +
    scale_linetype_manual(values = setNames(c(2, 1), c(lc, lk))) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(x = "1 - Specificity", y = "Sensitivity") +
    theme_publication() + theme(legend.direction = "vertical")

  cur <- av$decision_curve$curve
  nb <- rbind(
    data.frame(t = cur$threshold, nb = cur$clinical,   s = "Clinical"),
    data.frame(t = cur$threshold, nb = cur$combined,   s = "Clinical + immune"),
    data.frame(t = cur$threshold, nb = cur$treat_all,  s = "Treat all"),
    data.frame(t = cur$threshold, nb = cur$treat_none, s = "Treat none"))
  nb$s <- factor(nb$s, levels = c("Clinical", "Clinical + immune", "Treat all", "Treat none"))
  pB <- ggplot(nb, aes(t, nb, colour = s, linetype = s)) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(values = c("Clinical" = pub_palette[["clinical"]],
      "Clinical + immune" = pub_palette[["combined"]], "Treat all" = pub_palette[["treat_all"]],
      "Treat none" = pub_palette[["treat_none"]])) +
    scale_linetype_manual(values = c("Clinical" = 2, "Clinical + immune" = 1,
      "Treat all" = 1, "Treat none" = 3)) +
    coord_cartesian(ylim = c(min(-0.02, min(cur$treat_all)), max(c(cur$clinical, cur$combined)) + 0.02)) +
    labs(x = "Threshold probability", y = "Net benefit") +
    theme_publication() + theme(legend.direction = "vertical")
  pub_tag(pA | pB)
}

# ── MAIN Fig 3 — Calibration (unpenalized vs ridge) + IDI fragility ───────────
pub_fig_calibration_idi <- function(av) {
  if (is.null(av) || is.null(av$per_patient)) return(NULL)
  pp <- av$per_patient; pos <- av$positive_label
  y  <- as.integer(pp$True_Group == pos); p <- pp$Prob_Combined
  brk <- unique(quantile(p, seq(0, 1, 0.2), na.rm = TRUE))
  q   <- cut(p, breaks = brk, include.lowest = TRUE)
  bins <- do.call(rbind, lapply(split(seq_along(y), q), function(ix) {
    o <- mean(y[ix]); data.frame(pm = mean(p[ix]), obs = o, se = sqrt(o * (1 - o) / length(ix)))
  }))
  g   <- seq(min(p), max(p), length.out = 100)
  fit <- glm(y ~ qlogis(pub_clamp01(p)), family = binomial())
  cr  <- av$combined_optimism; rg <- av$combined_ridge
  line_u <- pub_invlogit(coef(fit)[1] + coef(fit)[2] * qlogis(pub_clamp01(g)))
  line_r <- pub_invlogit(rg$calib_citl_loo + rg$calib_slope_loo * qlogis(pub_clamp01(g)))
  cal <- rbind(
    data.frame(g, val = line_u, k = sprintf("Unpenalized (slope %.2f)", cr$calib_slope_loo)),
    data.frame(g, val = line_r, k = sprintf("Ridge (slope %.2f)", rg$calib_slope_loo)))
  pA <- ggplot() +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = pub_palette[["ref"]]) +
    geom_line(data = cal, aes(g, val, colour = k, linetype = k), linewidth = 0.8) +
    geom_point(data = bins, aes(pm, obs), size = 1.7, colour = pub_palette[["observed"]]) +
    geom_errorbar(data = bins, aes(x = pm, ymin = pmax(0, obs - se), ymax = pmin(1, obs + se)),
                  width = 0.012, colour = pub_palette[["observed"]]) +
    scale_colour_manual(values = setNames(c(pub_palette[["unpen"]], pub_palette[["ridge"]]), unique(cal$k))) +
    scale_linetype_manual(values = setNames(c(2, 1), unique(cal$k))) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(x = "Predicted probability", y = "Observed fraction") +
    theme_publication() + theme(legend.direction = "vertical")

  idi <- data.frame(
    model = c("Unpenalized\n(logistic)", "Ridge (L2)"),
    est   = c(av$increment$idi_loo, rg$idi_loo),
    lo    = c(av$increment$idi_loo_ci[1], rg$idi_loo_ci[1]),
    hi    = c(av$increment$idi_loo_ci[2], rg$idi_loo_ci[2]))
  idi$excl0 <- idi$lo > 0
  idi$model <- factor(idi$model, levels = rev(idi$model))
  pB <- ggplot(idi, aes(est, model, colour = excl0)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = pub_palette[["ref"]]) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18) +
    geom_point(size = 2.6) +
    geom_text(aes(label = sprintf("%+.3f [%.3f, %.3f]", est, lo, hi)),
              vjust = -1.0, size = 2.5, colour = "grey20") +
    scale_colour_manual(values = c("TRUE" = pub_palette[["ok"]], "FALSE" = pub_palette[["fragile"]]),
                        guide = "none") +
    labs(x = "IDI (leakage-free, LOO)", y = NULL) +
    coord_cartesian(clip = "off") +
    theme_publication() + theme(plot.margin = margin(6, 14, 6, 6))
  pub_tag(pA | pB)
}

# ── SUPP S1 — CONSORT / patient flow ──────────────────────────────────────────
#' @param cc named list: raw_n, raw_rp, raw_sdpd, n_miss, n_out, n, n_pos, n_neg,
#'                       n_paired, n_cc, surv_n, surv_os_ev, surv_pfs_ev, clin_label
pub_fig_consort <- function(cc) {
  if (is.null(cc)) return(NULL)
  box <- function(x, y, w, h, lab, fill = "white") list(
    r = data.frame(xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2, fill = fill),
    t = data.frame(x = x, y = y, lab = lab))
  seg <- function(x0, y0, x1, y1) data.frame(x = x0, y = y0, xe = x1, ye = y1)
  b1 <- box(3, 9, 4.6, 1.5, sprintf(
    "NSCLC anti-PD1 cohort (baseline T0)\nn = %d  ·  %d RP / %d SD-PD", cc$raw_n, cc$raw_rp, cc$raw_sdpd))
  bX <- box(6.4, 7.6, 3.4, 1.3, sprintf(
    "Excluded (n = %d)\n• High missingness (n = %d)\n• PCA outliers (n = %d)",
    cc$n_miss + cc$n_out, cc$n_miss, cc$n_out), fill = "grey95")
  b2 <- box(3, 6.2, 4.2, 1.3, sprintf(
    "Analytic cohort, T0\nn = %d  ·  %d RP / %d SD-PD", cc$n, cc$n_pos, cc$n_neg), fill = "#DCEAF5")
  b3 <- box(1.3, 3.3, 2.5, 1.5, sprintf("Paired T0+T1\n(LMM gate)\nn = %d", cc$n_paired))
  b4 <- box(3.6, 3.3, 2.7, 1.5, sprintf("Complete-case\nclinical model\n(%s)\nn = %d", cc$clin_label, cc$n_cc), fill = "#DCEAF5")
  b5 <- box(6.2, 3.3, 2.7, 1.5, sprintf("Survival follow-up\nn = %d\n(OS %d / PFS %d events)",
                                        cc$surv_n, cc$surv_os_ev, cc$surv_pfs_ev))
  rects <- do.call(rbind, lapply(list(b1,bX,b2,b3,b4,b5), `[[`, "r"))
  txts  <- do.call(rbind, lapply(list(b1,bX,b2,b3,b4,b5), `[[`, "t"))
  arr <- rbind(seg(3, 8.25, 3, 6.85), seg(3, 5.55, 3, 4.05),
               seg(3, 4.9, 1.3, 4.05), seg(3, 4.9, 6.2, 4.05),
               data.frame(x = 3, y = 7.6, xe = 4.7, ye = 7.6))
  ggplot() +
    geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = rects$fill, colour = "grey40", linewidth = 0.4) +
    geom_segment(data = arr, aes(x = x, y = y, xend = xe, yend = ye),
                 arrow = arrow(length = unit(1.6, "mm")), colour = "grey40", linewidth = 0.4) +
    geom_text(data = txts, aes(x, y, label = lab), size = 2.5, lineheight = 0.92) +
    coord_cartesian(xlim = c(-0.2, 8.2), ylim = c(2.3, 10)) +
    theme_void()
}

# ── SUPP S2 — PD-L1 context: strata bars + subgroup-AUC forest ────────────────
pub_fig_pdl1_context <- function(sr) {
  if (is.null(sr)) return(NULL)
  bins <- sr$bin_crosstab
  bins$Bin <- factor(bins$Bin, levels = sr$bin_labels)
  bins$lab <- sprintf("%.0f%%\n(N=%d, R=%d)", bins$Response_Rate, bins$N, bins$N_Responder)
  pA <- ggplot(bins, aes(Bin, Response_Rate)) +
    geom_col(fill = pub_palette[["combined"]], width = 0.62, alpha = 0.9) +
    geom_text(aes(label = lab), vjust = -0.3, size = 2.6, lineheight = 0.85) +
    geom_hline(yintercept = 50, linetype = "dashed", colour = pub_palette[["ref"]]) +
    scale_y_continuous(limits = c(0, max(bins$Response_Rate) + 20), expand = c(0, 0),
                       breaks = seq(0, 100, 25), labels = function(x) paste0(x, "%")) +
    labs(x = sr$label, y = "Responder rate") + theme_publication()
  sl <- sr$subgroup_low; sh <- sr$subgroup_high
  fo <- data.frame(
    grp = c(sprintf("PD-L1 < %g\n(n=%d)", sr$binary_cut$threshold, sl$n),
            sprintf("PD-L1 ≥ %g\n(n=%d)", sr$binary_cut$threshold, sh$n)),
    est = c(sl$auc, sh$auc),
    lo  = c(sl$auc_ci[1], sh$auc_ci[1]),
    hi  = c(sl$auc_ci[3], sh$auc_ci[3]))
  fo$grp <- factor(fo$grp, levels = rev(fo$grp))
  pB <- ggplot(fo, aes(est, grp)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = pub_palette[["ref"]]) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.16, colour = "grey35") +
    geom_point(size = 2.8, colour = pub_palette[["combined"]]) +
    geom_text(aes(label = sprintf("%.2f [%.2f, %.2f]", est, lo, hi)), vjust = -1.0, size = 2.5, colour = "grey20") +
    coord_cartesian(xlim = c(0, 1), clip = "off") +
    labs(x = "Gate-model AUC within PD-L1 subgroup", y = NULL) +
    theme_publication() + theme(plot.margin = margin(6, 12, 6, 6))
  pub_tag(pA | pB)
}

# ── SUPP S3 — Gate-signal decomposition (T0 / T1 / Delta) ─────────────────────
pub_fig_gate_signal <- function(gd) {
  if (is.null(gd) || is.null(gd$marker_decomp)) return(NULL)
  d <- gd$marker_decomp; d <- d[!is.na(d$AUC), ]
  d$Timepoint <- factor(d$Timepoint, levels = c("T0", "T1", "Delta (T1-T0)"))
  ggplot(d, aes(Timepoint, AUC, colour = Timepoint)) +
    geom_hline(yintercept = 0.5, linetype = "dotted", colour = pub_palette[["ref"]]) +
    geom_line(aes(group = Marker), colour = "grey70", linetype = "dashed") +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = AUC_CI_Lo, ymax = AUC_CI_Hi), width = 0.15) +
    facet_wrap(~ Marker) +
    scale_colour_manual(values = c("T0" = "#56B4E9", "T1" = "#009E73", "Delta (T1-T0)" = "#D55E00"),
                        guide = "none") +
    labs(x = NULL, y = "Univariate AUC (95% CI)") +
    theme_publication() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

# ── SUPP S4 — Baseline invariance of the increment ────────────────────────────
pub_fig_baseline_invariance <- function(av) {
  bsen <- av$baseline_sensitivity
  rows <- list(list(b = "Reference\nclinical model", lrt = av$increment$lrt_p, idi = av$increment$idi))
  for (nm in names(bsen)) rows[[length(rows) + 1]] <-
    list(b = sub("^\\+", "+ ", nm), lrt = bsen[[nm]]$lrt_p, idi = bsen[[nm]]$idi)
  d <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  d$b <- factor(d$b, levels = rev(d$b))
  ggplot(d, aes(lrt, b)) +
    geom_vline(xintercept = 0.05, linetype = "dashed", colour = pub_palette[["unpen"]]) +
    geom_point(size = 3, colour = pub_palette[["combined"]]) +
    geom_text(aes(label = sprintf("LRT p=%.4f  |  IDI %+.3f", lrt, idi)),
              hjust = -0.08, size = 2.6, colour = "grey20") +
    scale_x_log10(limits = c(0.004, 0.2)) +
    labs(x = "LRT p (log scale) - increment vs each baseline", y = NULL) +
    coord_cartesian(clip = "off") +
    theme_publication() + theme(plot.margin = margin(6, 60, 6, 6))
}

# ── SUPP S5 — Specificity vs random k-marker composites ───────────────────────
pub_fig_specificity_null <- function(av) {
  sn <- av$specificity_null; if (is.null(sn)) return(NULL)
  qd <- data.frame(pct = c(50, 90, 95, 99),
                   dauc = as.numeric(unlist(sn$null_dauc_q[c("50%","90%","95%","99%")])))
  ggplot(qd, aes(dauc, pct)) +
    geom_line(colour = "grey55") + geom_point(colour = "grey40", size = 2) +
    geom_vline(xintercept = sn$gate_delta_auc_apparent, colour = pub_palette[["combined"]], linewidth = 0.9) +
    annotate("text", x = sn$gate_delta_auc_apparent, y = 60, hjust = -0.05, size = 2.6,
             colour = pub_palette[["combined"]],
             label = sprintf("Ki67 gate ΔAUC=%.3f\nspec-p(LRT)=%.3f; %.0f%% random LRT<0.05",
                             sn$gate_delta_auc_apparent, sn$spec_p_lrt, 100 * sn$frac_sig_lrt)) +
    labs(x = sprintf("Null ΔAUC across %d random %d-marker composites", sn$n_random, sn$k_markers),
         y = "Null percentile") +
    theme_publication()
}

# ── SUPP S6 — Standalone classifier ROC (secondary; not the headline) ─────────
#' @param results_list named list of nested-LOOCV result objects (Elastic Net, SVM-RBF),
#'        each with $per_patient_predictions or out-of-fold probs + $positive_label,
#'        OR a data.frame of out-of-fold predictions.
#' @param preds_df data.frame: True_Group, Prob_ElasticNet, Prob_SVM_RBF
#' @param perm named list: en, svm permutation p-values; auc named list en, svm
pub_fig_standalone <- function(preds_df, positive_label, perm = NULL, auc = NULL) {
  if (is.null(preds_df)) return(NULL)
  y <- as.integer(preds_df$True_Group == positive_label)
  re <- pub_roc_df(preds_df$Prob_ElasticNet, y); rs <- pub_roc_df(preds_df$Prob_SVM_RBF, y)
  pe <- if (!is.null(perm$en))  sprintf(", perm p=%.3f", perm$en)  else ""
  ps <- if (!is.null(perm$svm)) sprintf(", perm p=%.3f", perm$svm) else ""
  le <- sprintf("Elastic Net (AUC %.2f%s)", re$auc, pe)
  ls <- sprintf("SVM-RBF (AUC %.2f%s)", rs$auc, ps)
  d <- rbind(data.frame(re$df, M = le), data.frame(rs$df, M = ls))
  d$M <- factor(d$M, levels = c(ls, le))
  ggplot(d, aes(FPR, TPR, colour = M)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = pub_palette[["ref"]]) +
    geom_path(linewidth = 0.8) +
    scale_colour_manual(values = setNames(c(pub_palette[["combined"]], pub_palette[["clinical"]]), c(ls, le))) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(x = "1 - Specificity", y = "Sensitivity") +
    theme_publication() + theme(legend.direction = "vertical")
}

# ── SUPP S7 — Dynamics<->baseline coupling (gate rationale; Blocco G) ─────────
#' Across the marker panel: LMM Timepoint×Group interaction strength (|t|) vs each
#' marker's standalone baseline (T0) AUC. Positive correlation = convergent validity
#' (why a dynamically-selected gate predicts at baseline). Gate markers highlighted;
#' the top T0-AUC markers are labelled to show the gate is NOT the baseline-max set.
pub_fig_coupling <- function(av) {
  cpl <- av$dynamics_baseline_coupling
  if (is.null(cpl) || is.null(cpl$per_marker) || !nrow(cpl$per_marker)) return(NULL)
  d <- cpl$per_marker
  d$grp <- ifelse(d$is_gate, "Immune gate", "Other panel marker")
  lab_top <- d$Marker[order(d$rank_T0)][1:min(2L, nrow(d))]   # top non-/gate T0 markers
  d$lab <- ifelse(d$is_gate | d$Marker %in% lab_top, d$Marker, "")
  yr <- range(d$T0_AUC); xr <- range(d$absT)
  ggplot(d, aes(absT, T0_AUC)) +
    geom_hline(yintercept = 0.5, linetype = "dotted", colour = pub_palette[["ref"]]) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                colour = "grey55", fill = "grey85", linewidth = 0.6) +
    geom_point(aes(colour = grp, size = grp)) +
    geom_text(aes(label = lab), vjust = -0.8, hjust = 0.5, size = 2.4, colour = "grey20") +
    annotate("text", x = xr[1], y = yr[2], hjust = 0, vjust = 1, size = 2.7, colour = "grey20",
             label = sprintf("Pearson r=%+.2f [%.2f, %.2f], p=%.3f\nexcl. gate r=%+.2f; Spearman rho=%+.2f",
                             cpl$pearson_r, cpl$pearson_ci[1], cpl$pearson_ci[2], cpl$pearson_p,
                             cpl$pearson_r_excl_gate, cpl$spearman_rho)) +
    scale_colour_manual(values = c("Immune gate" = pub_palette[["combined"]],
                                   "Other panel marker" = "grey60"), name = NULL) +
    scale_size_manual(values = c("Immune gate" = 2.6, "Other panel marker" = 1.7), guide = "none") +
    labs(x = "LMM Timepoint × Group interaction strength  |t|",
         y = "Standalone baseline (T0) AUC") +
    theme_publication()
}

# ── MAIN Fig 4 — Robustness / evidence-stability of the increment ─────────────
#' Defensibility centrepiece: (A) the LRT anchor stays significant across every
#' specification we tested; (B) the IDI estimator gradient is honestly fragile.
#' Consolidates the baseline-invariance, specificity-null and IDI panels.
pub_fig_robustness <- function(av) {
  if (is.null(av)) return(NULL)
  inc <- av$increment; bsen <- av$baseline_sensitivity; sn <- av$specificity_null
  # ---- Panel A: LRT p across specifications --------------------------------
  rowsA <- list(
    c("Asymptotic χ²",        inc$lrt_p),
    c("Permutation",          inc$lrt_perm_p),
    c("Firth penalized",      inc$lrt_firth_p))
  if (!is.null(bsen[["+smoking"]])) rowsA[[length(rowsA) + 1]] <- c("Baseline + smoking",      bsen[["+smoking"]]$lrt_p)
  if (!is.null(bsen[["+burden"]]))  rowsA[[length(rowsA) + 1]] <- c("Baseline + tumor burden", bsen[["+burden"]]$lrt_p)
  if (!is.null(sn))                 rowsA[[length(rowsA) + 1]] <- c("Specificity (vs random)", sn$spec_p_lrt)
  dA <- data.frame(spec = vapply(rowsA, `[`, "", 1),
                   p    = as.numeric(vapply(rowsA, `[`, "", 2)), stringsAsFactors = FALSE)
  dA <- dA[is.finite(dA$p), , drop = FALSE]
  dA$spec <- factor(dA$spec, levels = rev(dA$spec))
  pA <- ggplot(dA, aes(p, spec)) +
    geom_vline(xintercept = 0.05, linetype = "dashed", colour = pub_palette[["unpen"]]) +
    geom_point(size = 2.8, colour = pub_palette[["combined"]]) +
    geom_text(aes(label = sprintf("p=%.3f", p)), vjust = -1.0, size = 2.6, colour = "grey20") +
    scale_x_log10(limits = c(0.004, 0.2)) +
    labs(x = "Likelihood-ratio test p (log scale)", y = NULL) +
    coord_cartesian(clip = "off") +
    theme_publication() + theme(plot.margin = margin(6, 16, 6, 6))
  # ---- Panel B: IDI estimator gradient (apparent -> LOO -> ridge) ----------
  rowsB <- list(list("Apparent\n(in-sample)", inc$idi, inc$idi_ci),
                list("Leakage-free\n(LOO)",   inc$idi_loo, inc$idi_loo_ci))
  if (!is.null(av$combined_ridge))
    rowsB[[length(rowsB) + 1]] <- list("Ridge (L2)", av$combined_ridge$idi_loo, av$combined_ridge$idi_loo_ci)
  dB <- do.call(rbind, lapply(rowsB, function(r)
    data.frame(model = r[[1]], est = r[[2]], lo = r[[3]][1], hi = r[[3]][2], stringsAsFactors = FALSE)))
  dB$excl0 <- dB$lo > 0
  dB$model <- factor(dB$model, levels = rev(dB$model))
  pB <- ggplot(dB, aes(est, model, colour = excl0)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = pub_palette[["ref"]]) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18) +
    geom_point(size = 2.8) +
    geom_text(aes(label = sprintf("%+.3f [%.3f, %.3f]", est, lo, hi)),
              vjust = -1.0, size = 2.5, colour = "grey20") +
    scale_colour_manual(values = c("TRUE" = pub_palette[["ok"]], "FALSE" = pub_palette[["fragile"]]),
                        guide = "none") +
    labs(x = "IDI (95% CI)", y = NULL) +
    coord_cartesian(clip = "off") +
    theme_publication() + theme(plot.margin = margin(6, 14, 6, 6))
  pub_tag(pA | pB)
}

# ── render the whole figure set (single source of figure list + sizes) ────────
#' Build + save every manuscript figure from a bundle of live objects.
#' Used by the pipeline (Step 06) AND by the fast re-render script, so figure
#' sizes/styling live in exactly one place.
#' @param objs list: clin_addval, gate_decomp, stratified_result, lmm_boot
#'        (list boot/obs), consort, df_preds, positive_label, perm (en/svm).
pub_render_all <- function(objs, out_dir, project_name) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  pf <- function(n) file.path(out_dir, sprintf("%s_%s.pdf", n, project_name))
  av <- objs$clin_addval
  # MAIN — Fig 1 biology (cells + forest) / Fig 2 added-value / Fig 3 calibration /
  #        Fig 4 gate-signal decomposition (results). Robustness demoted to Supp.
  if (!is.null(objs$lmm_boot) || !is.null(objs$gate_decomp))
    save_pub_figure(pub_fig_biology(objs$gate_decomp, objs$lmm_boot$boot, objs$lmm_boot$obs,
                                    resp_label = objs$positive_label),
                    pf("Figure1_Biology"), PUB_W2, 120)
  if (!is.null(objs$gate_decomp))
    save_pub_figure(pub_fig_gate_signal(objs$gate_decomp),        pf("Figure4_GateSignal"),    PUB_W2, 110)
  if (!is.null(av)) {
    save_pub_figure(pub_fig_added_value(av),         pf("Figure2_AddedValue"),      PUB_W2, 120)
    # SUPP
    save_pub_figure(pub_fig_calibration_idi(av),     pf("FigureS9_Calibration_IDI"), PUB_W2, 120)
    save_pub_figure(pub_fig_robustness(av),          pf("FigureS8_Robustness"),     PUB_W2, 95)
    if (!is.null(objs$consort))
      save_pub_figure(pub_fig_consort(objs$consort), pf("FigureS1_CONSORT"),        PUB_W2, 120)
    save_pub_figure(pub_fig_baseline_invariance(av), pf("FigureS4_BaselineInvariance"), PUB_W2, 82)
    save_pub_figure(pub_fig_specificity_null(av),    pf("FigureS5_SpecificityNull"),    PUB_W1, 95)
    if (!is.null(av$dynamics_baseline_coupling))
      save_pub_figure(pub_fig_coupling(av),          pf("FigureS7_Coupling"),           PUB_W1, 95)
  }
  if (!is.null(objs$stratified_result))
    save_pub_figure(pub_fig_pdl1_context(objs$stratified_result), pf("Figure3_PDL1_context"), PUB_W2, 115)
  if (!is.null(objs$df_preds))
    save_pub_figure(pub_fig_standalone(objs$df_preds, objs$positive_label, perm = objs$perm),
                    pf("FigureS6_StandaloneClassifier"), PUB_W1, 120)
  invisible(out_dir)
}

# ── helper: CONSORT counts from live objects + same-run QC report ──────────────
pub_consort_counts <- function(config, av, gate_decomp) {
  # downstream counts from live objects
  cc <- list(
    n = av$n, n_pos = av$n_pos, n_neg = av$n_neg,
    n_cc = av$n_complete_case,
    n_paired = if (!is.null(gate_decomp)) gate_decomp$n_paired else NA_integer_,
    surv_n = av$survival$OS$n, surv_os_ev = av$survival$OS$events,
    surv_pfs_ev = av$survival$PFS$events,
    clin_label = gsub("_", "-", paste(names(av$clinical_vars), collapse = "+")),
    n_miss = NA_integer_, n_out = NA_integer_,
    raw_n = av$n, raw_rp = av$n_pos, raw_sdpd = av$n_neg)
  # upstream exclusions from the same-run Step-01 standard QC report (best-effort)
  qc <- tryCatch({
    f <- file.path(step_dir(config, 1, create = FALSE),
                   sprintf("QC_Filtering_Report_%s_standard.xlsx", config$project_name))
    if (!file.exists(f)) NULL else {
      dd <- as.data.frame(readxl::read_excel(f, sheet = "Details_Dropped", skip = 1))
      names(dd) <- c("Patient_ID", "NA_Percent", "Reason", "Original_Source")[seq_len(ncol(dd))]
      dd
    }
  }, error = function(e) NULL)
  if (!is.null(qc) && nrow(qc)) {
    cc$n_miss <- sum(grepl("missing", qc$Reason, ignore.case = TRUE))
    cc$n_out  <- sum(grepl("outlier", qc$Reason, ignore.case = TRUE))
    drp_rp    <- sum(qc$Original_Source == config$clinical$responder_label, na.rm = TRUE)
    drp_sd    <- nrow(qc) - drp_rp
    cc$raw_n   <- av$n + nrow(qc)
    cc$raw_rp  <- av$n_pos + drp_rp
    cc$raw_sdpd<- av$n_neg + drp_sd
  } else { cc$n_miss <- 0L; cc$n_out <- 0L }
  cc
}

# ── helper: read same-run Step-04 LMM bootstrap frames for the forest ─────────
pub_read_lmm_bootstrap <- function(config) {
  tryCatch({
    f <- file.path(step_dir(config, 4, create = FALSE),
                   sprintf("Longitudinal_LMM_Report_%s.xlsx", config$project_name))
    if (!file.exists(f)) return(NULL)
    sheets <- readxl::excel_sheets(f)
    if (!"LMM_Bootstrap_CI" %in% sheets) return(NULL)
    boot <- as.data.frame(readxl::read_excel(f, sheet = "LMM_Bootstrap_CI"))
    obs  <- if ("LMM_Interaction_Results" %in% sheets)
      as.data.frame(readxl::read_excel(f, sheet = "LMM_Interaction_Results"))[, c("Marker", "Estimate_Interaction")]
    else NULL
    list(boot = boot, obs = obs)
  }, error = function(e) NULL)
}
