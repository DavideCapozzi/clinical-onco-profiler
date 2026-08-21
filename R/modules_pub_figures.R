# R/modules_pub_figures.R
# ==============================================================================
# PUBLICATION FIGURE BUILDERS  (manuscript figure set)
# ------------------------------------------------------------------------------
# Each builder takes a LIVE pipeline result object (clin_addval, gate_decomp,
# stratified_result, LMM bootstrap frames, ...) and returns a ggplot/patchwork
# object styled via R/modules_pub_style.R. No frozen-file reads, no recompute:
# when the analysis changes, re-running the pipeline regenerates these figures.
#
# Narrative (rebalanced 2026-06-24 toward immunology, 3 main): MAIN Fig1 biology
# (gate-marker cell frequencies T0/T1×group, log-y + LMM forest) / Fig2 added-value
# (LRT-led) / Fig3 PD-L1 context. The gate-signal decomposition and the
# selection-aware circularity nulls are DEFENSE figures (statistical / invariance ⇒
# weak main-figure payoff) → SUPP (promote S9 to main only for a Trustworthy-AI
# venue). SUPP S1 CONSORT / S2 baseline-invariance / S3 specificity-null / S4
# standalone / S5 coupling / S6 robustness / S7 calibration / S8 gate-signal
# decomposition / S9 selection-aware circularity nulls (pub_fig_selection_aware:
# the Δ-headline anti-circularity evidence — replays the WHOLE gate selection under
# permuted labels; A = naive→selection-honest p ladder, B = increment/LRT null,
# C = nested-AUC null. Built OUT-OF-PIPELINE by manuscript/figures/render_selection_aware.R
# from the persisted diag_39b / diag_44 nulls — NOT pub_render_all) /
# S10 nomograms (pub_fig_nomogram, per timepoint; panel A = 3 gate markers individual, apparent
# display-only; panel B = immune composite + THIS run's clinical vars = the formal model, so
# include_nlr drives it — reported with its leakage-free LOO AUC. Immune axes are relabelled
# onto their exact clinical scale: Δ = fold change in the Ki67+:Ki67− ratio, T0 = % of parent
# gate; a bijection of the model scale, so nothing fitted changes — see nomo_tick_label()) /
# S11 classification performance (pub_fig_classification: confusion matrix + precision-recall from
# the combined model's leakage-free LOO probs at the pre-specified prevalence threshold; NOT Youden) /
# S12 raw fold-change of the gate markers (pub_fig_foldchange: clinician-facing, computed directly
# from raw cell % — responders halve Ki67+ subsets, non-responders rise ~20%; descriptive, formal
# test stays the Step-04 LMM; agrees with the logit-scale model at r=0.999 for these small fractions).
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
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    labs(x = "Time×Group β", y = NULL) +
    coord_cartesian(clip = "off") +
    theme_publication() + theme(plot.margin = margin(6, 14, 6, 6))
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
  # non-responder = vermillion; responder column first. Display-relabel group codes
  # (RP→PR, SD_PD→SD/PD) before faceting so strips/keys read in RECIST English.
  ppv$Group <- pub_relabel_group(ppv$Group)
  resp_disp <- if (!is.null(resp_label)) pub_relabel_group(resp_label) else NULL
  grps <- unique(as.character(ppv$Group))
  if (!is.null(resp_disp) && resp_disp %in% grps) {
    ord <- c(resp_disp, setdiff(grps, resp_disp))
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
    labs(x = NULL, y = "Frequency (% of parent population)") +
    theme_publication() + theme(panel.grid = element_blank())
  if (is.null(forest)) return(pub_tag(pA))
  pub_tag((pA | forest) + patchwork::plot_layout(widths = c(1.9, 1.1)))
}

# ── MAIN Fig 2 — Added value (LRT-led): ROC + decision curve ──────────────────
pub_fig_added_value <- function(av) {
  if (is.null(av) || is.null(av$per_patient)) return(NULL)
  pp  <- av$per_patient; pos <- av$positive_label
  y   <- as.integer(pp$True_Group == pos)
  inc <- av$increment
  PDL1_RED <- "#B2182B"; COMP_COL <- "#CC79A7"      # PD-L1 dark red; comparator reddish-purple
  clab <- if (!is.null(av$comparator$label)) av$comparator$label else "Clinical + immune + NLR"
  # Name the base clinical model by its constituents (e.g. "Clinical (PD-L1+PS)") so the
  # reader sees what 'clinical' contains; downstream curves build on it verbally.
  clin_str <- if (!is.null(av$formal_vars) && length(av$formal_vars))
                gsub("_", "-", paste(unlist(av$formal_vars), collapse = "+")) else NULL
  clin_nm  <- if (!is.null(clin_str)) sprintf("Clinical (%s)", clin_str) else "Clinical"
  nfin <- function(v) sum(is.finite(v))
  # honest p formatting (avoid 'p = 0.000' on asymptotic underflow)
  op_p  <- function(p) if (is.null(p) || !is.finite(p)) "= NA" else if (p < 0.001) "< 0.001" else sprintf("= %.3f", p)
  val_p <- function(p) if (is.null(p) || !is.finite(p)) "NA"   else if (p < 1e-4)  "< 1e-4"  else sprintf("%.4f", p)
  pdl1_ok  <- "Prob_PDL1" %in% names(pp) && any(is.finite(pp$Prob_PDL1))
  imm_ok   <- "Prob_Immune" %in% names(pp) && any(is.finite(pp$Prob_Immune))
  comp_ok  <- "Prob_ClinicalComp" %in% names(pp) && any(is.finite(pp$Prob_ClinicalComp))
  # NB every Prob_* column here is a LEAVE-ONE-OUT prediction (pl_* upstream), so
  # every curve — and every AUC in the legend — is LOO. Labelled as such because
  # LOO is NOT the reportable discrimination: for tied/discrete clinical vars it is
  # tie-artifact-prone (PS has 44% tied pairs) and under-reports the clinical arm,
  # which would INFLATE the visually-read increment. The reportable repeated
  # stratified 10-fold pair is annotated on the panel (see cv_lab below).
  rc  <- pub_roc_df(pp$Prob_Clinical, y); ro <- pub_roc_df(pp$Prob_Combined, y)
  lc  <- sprintf("%s (LOO AUC %.2f, n=%d)", clin_nm, rc$auc, nfin(pp$Prob_Clinical))
  lk  <- sprintf("Clinical + immune (LOO AUC %.2f, n=%d)", ro$auc, nfin(pp$Prob_Combined))
  roc_l <- list(data.frame(rc$df, M = lc), data.frame(ro$df, M = lk))
  lev   <- c(lc, lk); cols <- c(pub_palette[["clinical"]], pub_palette[["combined"]]); lts <- c(2, 1)
  if (comp_ok) {                                # display-only 'clinical + NLR' comparator
    rcmp <- pub_roc_df(pp$Prob_ClinicalComp, y)
    lm_  <- sprintf("%s (LOO AUC %.2f, n=%d)", clab, rcmp$auc, nfin(pp$Prob_ClinicalComp))
    roc_l <- c(roc_l, list(data.frame(rcmp$df, M = lm_)))
    lev <- c(lev, lm_); cols <- c(cols, COMP_COL); lts <- c(lts, 5)
  }
  if (pdl1_ok) {
    fin <- is.finite(pp$Prob_PDL1)              # complete-case only (PD-L1 missingness)
    rp <- pub_roc_df(pp$Prob_PDL1[fin], y[fin])
    lp <- sprintf("PD-L1 alone (LOO AUC %.2f, n=%d)", rp$auc, sum(fin))
    roc_l <- c(list(data.frame(rp$df, M = lp)), roc_l)
    lev <- c(lp, lev); cols <- c(PDL1_RED, cols); lts <- c(4, lts)
  }
  # Immune composite ALONE. Shown deliberately: out-of-sample it can sit slightly ABOVE the
  # combined model when the clinical baseline is near-chance (the clinical coefficients are
  # near-null, so they cost more in variance than they buy in fit). Omitting the arm while
  # its AUC sits in the tables is the kind of gap a reader is entitled to read as concealment;
  # the honest move is to plot it and quantify the gap (cv_kfold$delta_auc_vs_immune —
  # quantified in the Figure 2 caption; see the note at the panel-A annotate() below).
  if (imm_ok) {
    ri  <- pub_roc_df(pp$Prob_Immune, y)
    li  <- sprintf("Immune composite alone (LOO AUC %.2f, n=%d)", ri$auc, nfin(pp$Prob_Immune))
    roc_l <- c(roc_l, list(data.frame(ri$df, M = li)))
    lev <- c(lev, li); cols <- c(cols, pub_palette[["immune"]]); lts <- c(lts, 6)
  }
  roc <- do.call(rbind, roc_l); roc$M <- factor(roc$M, levels = lev)
  # Reportable discrimination = repeated stratified k-fold; the curves above are LOO.
  # Stated on the panel so the reader never reads the increment off the LOO curves.
  # Guarded: runs / persisted rds predating the cv_kfold node simply omit the line.
  cv      <- av$cv_kfold
  cv_clin <- suppressWarnings(as.numeric(cv$auc_clinical)[1])
  cv_comb <- suppressWarnings(as.numeric(cv$auc_combined)[1])
  cv_imm  <- suppressWarnings(as.numeric(cv$auc_immune)[1])
  cv_k    <- suppressWarnings(as.numeric(cv$k)[1])
  # Immune-alone is stated next to the pair, not hidden behind it: when the clinical arm is
  # near-chance the combined model can fall slightly BELOW it, and the paired delta (same
  # folds) is what tells the reader whether that ordering is signal or noise.
  cv_imm_lab <- if (length(cv_imm) && isTRUE(is.finite(cv_imm)))
    sprintf(", immune alone %.2f", cv_imm) else ""
  cv_lab  <- if (length(cv_clin) && length(cv_comb) &&
                 isTRUE(is.finite(cv_clin)) && isTRUE(is.finite(cv_comb)))
    sprintf("Reportable %d-fold CV AUC: %.2f → %.2f%s (curves = LOO)\n",
            if (isTRUE(is.finite(cv_k))) as.integer(cv_k) else 10L,
            cv_clin, cv_comb, cv_imm_lab) else ""
  # The paired combined−immune interval (cv_kfold$delta_auc_vs_immune) is NOT annotated
  # on the panel: the on-plot text block ran to three lines and the two curves already sit
  # visibly on top of each other, so the panel makes the point without the number. It is
  # NOT dropped — it travels in the Figure 2 caption (manuscript/figure_captions.md), the
  # JSON and the Excel. Do not remove the immune-alone CURVE to match: the arm being
  # visible is what keeps this a presentation choice rather than a concealment (guardrail 2).
  pA <- ggplot(roc, aes(FPR, TPR, colour = M, linetype = M)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = pub_palette[["ref"]]) +
    geom_path(linewidth = 0.8) +
    annotate("text", x = 0.97, y = 0.06, hjust = 1, vjust = 0, size = 2.7,
             label = sprintf("%sLRT p %s (perm %s)", cv_lab,
                             op_p(inc$lrt_p), val_p(inc$lrt_perm_p))) +
    scale_colour_manual(values = setNames(cols, lev)) +
    scale_linetype_manual(values = setNames(lts, lev)) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(x = "1 - Specificity", y = "Sensitivity") +
    theme_publication() + theme(legend.direction = "vertical")

  cur <- av$decision_curve$curve
  nb_l <- list(data.frame(t = cur$threshold, nb = cur$clinical, s = clin_nm),
               data.frame(t = cur$threshold, nb = cur$combined, s = "Clinical + immune"))
  b_cols <- setNames(c(pub_palette[["clinical"]], pub_palette[["combined"]]), c(clin_nm, "Clinical + immune"))
  b_lts  <- setNames(c(2, 1), c(clin_nm, "Clinical + immune"))
  if (!is.null(cur$comparator)) {
    nb_l <- c(nb_l, list(data.frame(t = cur$threshold, nb = cur$comparator, s = clab)))
    b_cols <- c(b_cols, setNames(COMP_COL, clab)); b_lts <- c(b_lts, setNames(5, clab))
  }
  if (!is.null(cur$pdl1)) {
    nb_l <- c(nb_l, list(data.frame(t = cur$threshold, nb = cur$pdl1, s = "PD-L1 alone")))
    b_cols <- c(b_cols, "PD-L1 alone" = PDL1_RED); b_lts <- c(b_lts, "PD-L1 alone" = 4)
  }
  # Immune alone on the decision curve too — DCA is co-primary with the LRT, so the arm that
  # may dominate on net benefit has to be visible here, not only on the ROC.
  if (!is.null(cur$immune)) {
    nb_l <- c(nb_l, list(data.frame(t = cur$threshold, nb = cur$immune, s = "Immune alone")))
    b_cols <- c(b_cols, "Immune alone" = pub_palette[["immune"]])
    b_lts  <- c(b_lts, "Immune alone" = 6)
  }
  nb_l <- c(nb_l, list(data.frame(t = cur$threshold, nb = cur$treat_all,  s = "Treat all"),
                       data.frame(t = cur$threshold, nb = cur$treat_none, s = "Treat none")))
  b_cols <- c(b_cols, "Treat all" = pub_palette[["treat_all"]], "Treat none" = pub_palette[["treat_none"]])
  b_lts  <- c(b_lts, "Treat all" = 1, "Treat none" = 3)
  nb <- do.call(rbind, nb_l)
  nb$s <- factor(nb$s, levels = c(clin_nm, "Clinical + immune",
                                  if (!is.null(cur$comparator)) clab,
                                  if (!is.null(cur$pdl1)) "PD-L1 alone",
                                  if (!is.null(cur$immune)) "Immune alone",
                                  "Treat all", "Treat none"))
  pB <- ggplot(nb, aes(t, nb, colour = s, linetype = s)) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(values = b_cols) +
    scale_linetype_manual(values = b_lts) +
    # cur$immune belongs in the max(): it is the highest net-benefit arm whenever the
    # clinical baseline is near-chance, and would otherwise be clipped off the panel.
    coord_cartesian(ylim = c(min(-0.02, min(cur$treat_all)),
                             max(c(cur$clinical, cur$combined, cur$comparator,
                                   cur$immune)) + 0.02)) +
    # Net benefit is computed on the same LOO probabilities as panel A, so the
    # clinical reference arm carries the same tie artifact — disclosed on-panel
    # rather than left for the caption, since the curve is read directly.
    annotate("text", x = -Inf, y = -Inf, hjust = -0.04, vjust = -0.8, size = 2.3,
             colour = "grey30", label = "Net benefit from LOO probabilities") +
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
# Renders the COHORT LEDGER (R/utils_io.R): one row per real filtering point, recorded
# where the filtering happens. The previous version reconstructed the flow arithmetically
# from mismatched denominators (analytic-cohort exclusions added to the paired-subset n)
# and drew a starting cohort that never existed. Nothing is inferred here any more: the
# figure only renders rows, so the number of stages is data-driven and a cohort with no
# paired timepoint (cross-sectional) renders correctly with no code change.
#' @param led Cohort ledger data.frame: stage, n_in, n_out, n_dropped, reason, n_resp, n_nonresp.
#' @param terminal Optional named list of terminal ANNOTATIONS (not exclusions), e.g.
#'   list("Complete-case PD-L1+PS" = "n = 57", "Survival follow-up" = "n = 60 (41 OS events)").
#' @param title Optional label for the first box.
pub_fig_consort <- function(led, terminal = NULL, title = NULL) {
  if (is.null(led) || !nrow(led)) return(NULL)
  k    <- nrow(led)
  top  <- 10; step <- 2.2                       # one flow row per ledger stage
  ys   <- top - step * seq_len(k)               # y of each post-stage box
  grp  <- function(i) if (is.na(led$n_resp[i]) || is.na(led$n_nonresp[i])) "" else
    sprintf("  ·  %d / %d", led$n_resp[i], led$n_nonresp[i])
  nice <- function(s) gsub("_", " ", s)
  box <- function(x, y, w, h, lab, fill = "white") list(
    r = data.frame(xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2, fill = fill),
    t = data.frame(x = x, y = y, lab = lab))

  boxes <- list(box(3, top, 5.0, 1.2, sprintf("%s\nn = %d%s",
                    if (is.null(title)) "Screened cohort" else title,
                    led$n_in[1], if (k) "" else "")))
  arrs  <- list()
  for (i in seq_len(k)) {
    boxes[[length(boxes) + 1]] <- box(
      3, ys[i], 5.0, 1.2,
      sprintf("%s\nn = %d%s", nice(led$stage[i]), led$n_out[i], grp(i)),
      fill = if (i == k) "#DCEAF5" else "white")
    arrs[[length(arrs) + 1]] <- data.frame(x = 3, y = (if (i == 1) top else ys[i - 1]) - 0.6,
                                           xe = 3, ye = ys[i] + 0.6)
    if (led$n_dropped[i] > 0) {                 # side box only when patients actually left
      boxes[[length(boxes) + 1]] <- box(
        7.0, (if (i == 1) top else ys[i - 1]) - step / 2, 3.6, 1.15,
        sprintf("Excluded (n = %d)\n%s", led$n_dropped[i],
                if (is.na(led$reason[i])) "" else strwrap_lab(led$reason[i], 34)),
        fill = "grey95")
      arrs[[length(arrs) + 1]] <- data.frame(x = 3, y = (if (i == 1) top else ys[i - 1]) - step / 2,
                                             xe = 5.2, ye = (if (i == 1) top else ys[i - 1]) - step / 2)
    }
  }
  # Terminal annotations describe the FINAL cohort; they are not exclusions and must never
  # be drawn in the vertical flow (nobody is dropped for an imputed clinical value).
  if (length(terminal)) {
    tl <- paste(sprintf("%s: %s", names(terminal), unlist(terminal)), collapse = "\n")
    boxes[[length(boxes) + 1]] <- box(3, ys[k] - step, 5.0, 0.5 + 0.42 * length(terminal), tl)
    arrs[[length(arrs) + 1]] <- data.frame(x = 3, y = ys[k] - 0.6, xe = 3,
                                           ye = ys[k] - step + 0.35 + 0.21 * length(terminal))
  }
  rects <- do.call(rbind, lapply(boxes, `[[`, "r"))
  txts  <- do.call(rbind, lapply(boxes, `[[`, "t"))
  arr   <- do.call(rbind, arrs)
  ggplot() +
    geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = rects$fill, colour = "grey40", linewidth = 0.4) +
    geom_segment(data = arr, aes(x = x, y = y, xend = xe, yend = ye),
                 arrow = arrow(length = unit(1.6, "mm")), colour = "grey40", linewidth = 0.4) +
    geom_text(data = txts, aes(x, y, label = lab), size = 2.4, lineheight = 0.95) +
    coord_cartesian(xlim = c(0.2, 9.0),
                    ylim = c(min(rects$ymin) - 0.4, top + 0.9)) +
    theme_void()
}

# wrap a long exclusion reason onto <= `w`-char lines so side boxes stay inside the panel
strwrap_lab <- function(s, w = 34) paste(strwrap(s, width = w), collapse = "\n")

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
    # wrap strata ticks onto two lines so neg/low/high don't collide horizontally;
    # prettify ">=" → "≥" for display (data labels keep the ASCII form)
    scale_x_discrete(labels = function(v)
      sub("(", "\n(", gsub(">=", "≥", v, fixed = TRUE), fixed = TRUE)) +
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
    geom_text(aes(label = sprintf("%.2f [%.2f, %.2f]", est, lo, hi)),
              vjust = -1.5, size = 2.5, colour = "grey20") +
    scale_y_discrete(expand = expansion(add = c(0.6, 1.0))) +
    coord_cartesian(xlim = c(0, 1), clip = "off") +
    labs(x = "AUC within PD-L1 subgroup", y = NULL) +
    theme_publication() + theme(plot.margin = margin(6, 14, 6, 6))
  pub_tag(pA | pB)
}


# ── SUPP — Selection-aware circularity nulls (Δ headline) ─────────────────────
# Renders from the two persisted diagnostics nulls (NOT live pipeline objects,
# so it is driven by a standalone script, not pub_render_all):
#   nested_null = diag_39b_nested_null.rds  (fully-nested AUC permutation null)
#   lrt_null    = diag_44_selection_aware_lrt.rds (selection-aware LRT null)
# The circularity attack: the gate was chosen by a Time×Group interaction and
# Group IS best response, so markers were picked using the outcome and tested on
# the same outcome/patients. Both nulls REPLAY the whole selection under permuted
# labels, so the p-value pays for the selection step (not just the model fit).
# Three panels: A = the naive→selection-honest p ladder (the interpretive panel);
# B = the LRT (increment) null; C = the nested-AUC (discrimination) null. Fat
# tails (null maxima) are drawn explicitly — they are the honest caveat.
pub_fig_selection_aware <- function(nested_null, lrt_null) {
  if (is.null(nested_null) || is.null(lrt_null)) return(NULL)

  # ---- recompute every number from the rds so nothing can drift ----
  auc_perm <- nested_null$perm; auc_obs <- nested_null$observed
  p_auc    <- (sum(auc_perm >= auc_obs) + 1) / (length(auc_perm) + 1)
  auc_med  <- median(auc_perm); auc_max <- max(auc_perm)

  lr_null  <- lrt_null$lr_null; lr_obs <- lrt_null$lr_obs   # FDR rule = the pipeline's own rule
  p_lrt    <- (1 + sum(lr_null >= lr_obs)) / (length(lr_null) + 1)   # 0.0080 = the headline p
  lr_null3 <- lrt_null$lr_null3                            # top-3 rule = the non-degenerate companion
  p_lrt3   <- (1 + sum(lr_null3 >= lr_obs)) / (length(lr_null3) + 1) # 0.0100 (agrees with FDR)
  naive_p  <- pchisq(lr_obs, df = 1, lower.tail = FALSE)   # the gate-FIXED asymptotic LRT
  lr3_max  <- max(lr_null3)
  frac_empty_fdr <- mean(lr_null == 0)                     # FDR gate empty under permuted labels
  inflation <- signif(p_lrt / naive_p, 1)                  # ~200x

  pfmt <- function(p) ifelse(p < 1e-4, sprintf("%.1e", p), sprintf("%.4f", p))

  # ---- Panel A — the p-value ladder (naive vs selection-honest) ----
  lad <- data.frame(
    y     = c(3, 2, 1),   # Naive top, then the two honest tests
    label = c("Naive test\n(gate held fixed)",
              "Selection-honest\nincrement (LRT)",
              "Selection-honest\ndiscrimination (AUC)"),
    p     = c(naive_p, p_lrt, p_auc),
    kind  = c("Optimistic (leaks selection)", "Honest", "Honest"),
    stringsAsFactors = FALSE)
  pA <- ggplot(lad, aes(p, y, colour = kind)) +
    geom_vline(xintercept = 0.05, linetype = "dashed", colour = "grey50") +
    annotate("text", x = 0.05, y = 3.62, label = "p = 0.05", size = 2.5,
             colour = "grey40", vjust = 0, hjust = 0.5) +
    geom_point(size = 3.4) +
    geom_text(aes(label = pfmt(p)), vjust = -1.1, size = 2.7, colour = "grey20") +
    annotate("segment", x = naive_p * 1.4, xend = p_lrt / 1.4, y = 2.5, yend = 2.5,
             colour = "grey35", arrow = grid::arrow(length = unit(1.6, "mm"), ends = "both")) +
    annotate("text", x = sqrt(naive_p * p_lrt), y = 2.5, vjust = -0.7, size = 2.6,
             colour = "grey25", label = sprintf("selection inflates p ~%dx", inflation)) +
    scale_colour_manual(values = c("Optimistic (leaks selection)" = pub_palette[["unpen"]],
                                   "Honest" = pub_palette[["combined"]]), name = NULL) +
    scale_x_log10(breaks = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1),
                  labels = c("1e-5", "1e-4", "1e-3", "0.01", "0.1"),
                  limits = c(naive_p / 3, 0.3)) +
    scale_y_continuous(breaks = lad$y, labels = lad$label,
                       limits = c(0.6, 3.8)) +
    labs(x = "p-value (log scale) — lower = stronger evidence", y = NULL) +
    theme_publication() + theme(plot.margin = margin(6, 12, 4, 6),
                                legend.position = "bottom")

  # ---- Panel B — the LRT (increment) null ----
  # Shows the top-3-rule null (a full, legible distribution). The pipeline's own
  # FDR-rule null is degenerate — its gate is EMPTY in ~96% of permutations, so
  # there is usually nothing to add — which is why the FDR p (0.0080) is even
  # slightly smaller than the top-3 p (0.0100). Both p's annotated.
  pB <- ggplot(data.frame(lr = lr_null3), aes(lr)) +
    geom_histogram(bins = 34, fill = "grey82", colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = lr_obs, colour = pub_palette[["combined"]], linewidth = 1) +
    annotate("text", x = lr_obs, y = Inf, vjust = 1.4, hjust = 1.06, size = 2.5,
             colour = pub_palette[["combined"]],
             label = sprintf("observed %.1f\np = %s (FDR)\np = %s (top-3)",
                             lr_obs, pfmt(p_lrt), pfmt(p_lrt3))) +
    annotate("text", x = lr3_max, y = -Inf, vjust = -0.6, hjust = 1.0, size = 2.3,
             colour = "grey45", label = sprintf("null max %.1f", lr3_max)) +
    labs(x = "Increment likelihood-ratio statistic\n(permuted-label null)",
         y = "Permutations") +
    theme_publication()

  # ---- Panel C — the nested-AUC (discrimination) null ----
  pC <- ggplot(data.frame(auc = auc_perm), aes(auc)) +
    geom_histogram(bins = 30, fill = "grey82", colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = 0.5, linetype = "dotted", colour = "grey40") +
    geom_vline(xintercept = auc_obs, colour = pub_palette[["immune"]], linewidth = 1) +
    annotate("text", x = 0.49, y = -Inf, vjust = -0.7, hjust = 1, size = 2.3,
             colour = "grey40", label = sprintf("null median %.2f\n(= chance)", auc_med)) +
    annotate("text", x = auc_obs, y = Inf, vjust = 1.4, hjust = 1.08, size = 2.6,
             colour = pub_palette[["immune"]],
             label = sprintf("observed %.3f\np = %s", auc_obs, pfmt(p_auc))) +
    annotate("text", x = 0.84, y = -Inf, vjust = -0.7, hjust = 1, size = 2.3,
             colour = "grey45", label = sprintf("null max %.2f", auc_max)) +
    scale_x_continuous(limits = c(0, 0.85), breaks = c(0, 0.25, 0.5, 0.75)) +
    labs(x = "Nested cross-validated AUC\n(permuted-label null)",
         y = "Permutations") +
    theme_publication()

  pub_tag(pA / (pB | pC) + patchwork::plot_layout(heights = c(1, 1.15)))
}

# ── SUPP — Gate-signal decomposition (T0 / T1 / Delta) ────────────────────────
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

# ── helper: adaptive log10 p-value axis ───────────────────────────────────────
# T0 LRT p's sat in ~0.004-0.2, so the axis limits used to be hardcoded there.
# Under Δ-primary the increment p's are 10-100x smaller and the asymptotic χ² LRT
# can underflow to exactly 0 (log10(0) = -Inf → point silently dropped). Floor
# positions AND labels at PVAL_FLOOR (an honest "<1e-4" reporting floor for an
# underflowed asymptotic p) and derive the axis limits from the data, so every
# point stays on-axis and the p=0.05 reference line remains visible at either scale.
PVAL_FLOOR <- 1e-4
pub_pval_x   <- function(p) pmax(p, PVAL_FLOOR)                       # on-axis position (no -Inf)
pub_pval_fmt <- function(p) ifelse(!is.finite(p) | p < PVAL_FLOOR,   # honest label for underflow
                                   sprintf("<%g", PVAL_FLOOR), sprintf("%.4f", p))
pub_pval_log_scale <- function(p, upper = 0.2) {
  pv <- pub_pval_x(p[is.finite(p)])
  scale_x_log10(limits = c(min(pv) / 1.8, max(upper, max(pv)) * 1.3))
}

# ── SUPP S2 — Baseline invariance of the increment ────────────────────────────
pub_fig_baseline_invariance <- function(av) {
  bsen <- av$baseline_sensitivity
  rows <- list(list(b = "Reference\nclinical model", lrt = av$increment$lrt_p, idi = av$increment$idi))
  for (nm in names(bsen)) rows[[length(rows) + 1]] <-
    list(b = sub("^\\+", "+ ", nm), lrt = bsen[[nm]]$lrt_p, idi = bsen[[nm]]$idi)
  d <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  d$b <- factor(d$b, levels = rev(d$b))
  ggplot(d, aes(pub_pval_x(lrt), b)) +
    geom_vline(xintercept = 0.05, linetype = "dashed", colour = pub_palette[["unpen"]]) +
    geom_point(size = 3, colour = pub_palette[["combined"]]) +
    geom_text(aes(label = sprintf("LRT p=%s  |  IDI %+.3f", pub_pval_fmt(lrt), idi)),
              hjust = -0.08, size = 2.6, colour = "grey20") +
    pub_pval_log_scale(d$lrt) +
    labs(x = "LRT p (log scale) - increment vs each baseline", y = NULL) +
    coord_cartesian(clip = "off") +
    theme_publication() + theme(plot.margin = margin(6, 60, 6, 6))
}

# ── SUPP S3 — Specificity vs random k-marker composites ───────────────────────
pub_fig_specificity_null <- function(av) {
  sn <- av$specificity_null; if (is.null(sn)) return(NULL)
  qd <- data.frame(pct = c(50, 90, 95, 99),
                   dauc = as.numeric(unlist(sn$null_dauc_q[c("50%","90%","95%","99%")])))
  ggplot(qd, aes(dauc, pct)) +
    geom_line(colour = "grey55") + geom_point(colour = "grey40", size = 2) +
    geom_vline(xintercept = sn$gate_delta_auc_apparent, colour = pub_palette[["combined"]], linewidth = 0.9) +
    # gate ΔAUC lands beyond the 99% null quantile → the vline sits at the right
    # edge; annotate to the LEFT of the line (hjust=1) so the label reads into the
    # empty upper-left area instead of clipping off the panel.
    annotate("text", x = sn$gate_delta_auc_apparent, y = 60, hjust = 1.05, size = 2.6,
             colour = pub_palette[["combined"]],
             label = sprintf("Ki67 gate ΔAUC=%.3f\nspec-p(LRT)=%.3f; %.0f%% random LRT<0.05",
                             sn$gate_delta_auc_apparent, sn$spec_p_lrt, 100 * sn$frac_sig_lrt)) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.10))) +
    labs(x = sprintf("Null ΔAUC across %d random %d-marker composites", sn$n_random, sn$k_markers),
         y = "Null percentile") +
    theme_publication()
}

# ── SUPP S4 — Standalone classifier ROC (secondary; not the headline) ─────────
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

# ── SUPP S6 — Robustness / evidence-stability of the increment ────────────────
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
  # Loop over whatever baseline specifications the config defined, rather than naming two of
  # them: config-driven specs used to be silently dropped from this panel, which is precisely
  # the invariance claim the panel exists to make.
  bsen_lab <- c("+smoking" = "Baseline + smoking", "+burden" = "Baseline + tumor burden")
  for (bn in names(bsen)) {
    if (is.null(bsen[[bn]]$lrt_p)) next
    lab <- if (bn %in% names(bsen_lab)) unname(bsen_lab[bn]) else sprintf("Baseline %s", bn)
    rowsA[[length(rowsA) + 1]] <- c(lab, bsen[[bn]]$lrt_p)
  }
  if (!is.null(sn))                 rowsA[[length(rowsA) + 1]] <- c("Specificity (vs random)", sn$spec_p_lrt)
  dA <- data.frame(spec = vapply(rowsA, `[`, "", 1),
                   p    = as.numeric(vapply(rowsA, `[`, "", 2)), stringsAsFactors = FALSE)
  dA <- dA[is.finite(dA$p), , drop = FALSE]
  dA$spec <- factor(dA$spec, levels = rev(dA$spec))
  pA <- ggplot(dA, aes(pub_pval_x(p), spec)) +
    geom_vline(xintercept = 0.05, linetype = "dashed", colour = pub_palette[["unpen"]]) +
    geom_point(size = 2.8, colour = pub_palette[["combined"]]) +
    geom_text(aes(label = sprintf("p=%s", pub_pval_fmt(p))), vjust = -1.0, size = 2.6, colour = "grey20") +
    pub_pval_log_scale(dA$p) +
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
# ── SUPP — nomograms (A: display-only apparent fit; B: the run's formal model) ──
#' Draw one classic points-based nomogram from a `build_nomogram_spec()` object.
#' Rows top→bottom: Points ruler, one ruler per predictor, Total-points ruler,
#' Predicted-probability ruler. The upper block (Points + predictors) shares the
#' 0..points_max scale; the lower block (Total points + probability) shares the
#' 0..total_max scale. Both are drawn on a common [0,1] canvas.
pub_nomo_panel <- function(s, header = NULL) {
  if (is.null(s) || is.null(s$predictors) || !length(s$predictors)) return(NULL)
  FONT    <- "sans"                              # cairo → Helvetica/Arial at embed
  CHARC   <- "#222222"                           # all TEXT: high-contrast charcoal
  col_imm <- unname(pub_palette[["immune"]])     # gate-marker axis/tick colour (green)
  col_cln <- unname(pub_palette[["clinical"]])    # clinical-var axis/tick colour (orange)
  col_track <- "grey80"                          # full-width guide rail behind each axis
  band_fill <- "#e9edf1"                          # crisper zebra band vs white canvas
  LW_AXIS <- 0.9; LW_TICK <- 0.7; LW_TRACK <- 0.6   # heavier strokes → sharp at 300+ dpi
  # Tick labels are precomputed on the display scale by build_nomogram_spec() (all unit
  # math lives there — see nomo_tick_label). fmt_val is the fallback for specs persisted
  # before the `label` column existed.
  fmt_val <- function(v, type) {
    if (identical(type, "z")) sprintf("%.1f", v)
    else if (all(abs(v - round(v)) < 1e-6)) sprintf("%.0f", v)
    else sprintf("%.1f", v)
  }
  tick_labels <- function(d) {
    if (!is.null(d$label)) return(as.character(d$label))
    vapply(seq_len(nrow(d)), function(i) fmt_val(d$value[i], d$scale_type[i]), character(1))
  }
  preds   <- s$predictors
  npred   <- length(preds)
  n_row   <- npred + 3                         # Points + preds + Total points + Prob
  y_of    <- function(top_idx) n_row - top_idx + 1   # top row = highest y
  TICK    <- 0.14
  X0 <- -0.34; X1 <- 1.04                       # canvas bounds (labels | rulers)
  # Keep the far endpoint, then greedily drop interior LABELS closer than min_gap
  # (tick marks are all retained) — guarantees legible labels on collapsed axes.
  thin_lab <- function(xs, txt, min_gap = 0.055) {
    o <- order(xs); xs <- xs[o]; txt <- txt[o]; n <- length(xs)
    keep <- logical(n); keep[1] <- TRUE; last <- xs[1]
    if (n > 1) for (i in 2:n) if (xs[i] - last >= min_gap) { keep[i] <- TRUE; last <- xs[i] }
    keep[n] <- TRUE
    if (n > 2 && (xs[n] - xs[max(which(keep[-n]))]) < min_gap) keep[max(which(keep[-n]))] <- FALSE
    data.frame(x = xs[keep], label = txt[keep])
  }
  axes <- arrows <- tracks <- ticks <- labs <- rowlab <- bullets <- bands <- notes <- list()
  # Scale rulers (Points / Total points) get a directional arrowhead; predictor and
  # probability rulers are plain segments. Ticks hang BELOW the ruler, labels above.
  # Tick MARKS carry the axis colour; all TEXT is charcoal for maximum contrast.
  add_ruler  <- function(y, x0, x1, col, arrow = FALSE) {
    row <- data.frame(x = x0, xe = x1, y = y, ye = y, col = col)
    if (arrow) arrows[[length(arrows) + 1]] <<- row
    else       axes[[length(axes)   + 1]] <<- row
  }
  add_track  <- function(y) tracks[[length(tracks) + 1]] <<- data.frame(x = 0, xe = 1, y = y, ye = y)
  add_ticks  <- function(y, xs, txt, col) {
    ticks[[length(ticks) + 1]] <<- data.frame(x = xs, xe = xs, y = y - TICK, ye = y, col = col)
    tl <- thin_lab(xs, txt)
    if (any(nzchar(tl$label)))
      labs[[length(labs) + 1]] <<- data.frame(x = tl$x, y = y + 0.20, label = tl$label)
  }
  add_rowlab <- function(y, txt, bullet = NA) {
    rowlab[[length(rowlab) + 1]] <<- data.frame(x = -0.06, y = y, label = txt)
    if (!is.na(bullet)) bullets[[length(bullets) + 1]] <<- data.frame(x = -0.03, y = y, col = bullet)
  }
  add_note   <- function(x, y, txt) notes[[length(notes) + 1]] <<- data.frame(x = x, y = y, label = txt)

  # Row 1: Points ruler (0..points_max), scaled to [0,1]
  yP <- y_of(1)
  pt_ticks <- pretty(c(0, s$points_max), n = 6)
  pt_ticks <- pt_ticks[pt_ticks >= 0 & pt_ticks <= s$points_max]
  add_ruler(yP, 0, 1, CHARC, arrow = TRUE); add_ticks(yP, pt_ticks / s$points_max, sprintf("%.0f", pt_ticks), CHARC)
  add_rowlab(yP, "Points")

  # Predictor rulers (share the Points scale → x = points / points_max). A full-width
  # grey guide RAIL underlays every axis so short (low-weight) axes are not orphaned in
  # empty space — the coloured active segment still encodes the true points range, so
  # the mathematical weights are unchanged. Gate markers green, clinical vars orange;
  # a colour bullet keeps the group cue while the name stays charcoal.
  for (j in seq_len(npred)) {
    d   <- preds[[j]]; y <- y_of(1 + j)
    col <- if (d$scale_type[1] %in% c("z", "fc", "pct")) col_imm else col_cln
    if (j %% 2 == 1)
      bands[[length(bands) + 1]] <- data.frame(xmin = X0, xmax = X1, ymin = y - 0.5, ymax = y + 0.5)
    add_track(y)
    mp  <- max(d$points, na.rm = TRUE)
    xr  <- range(d$points, na.rm = TRUE) / s$points_max
    add_ruler(y, xr[1], xr[2], col)
    if (mp < 4) {                               # negligible axis → clean charcoal note
      add_ticks(y, xr, c("", ""), col)
      add_note(xr[2] + 0.03, y, "≈ 0 points (negligible)")
    } else {
      add_ticks(y, d$points / s$points_max, tick_labels(d), col)
    }
    add_rowlab(y, d$display[1], bullet = col)
  }

  # Total points ruler (0..total_max → [0,1])
  yT <- y_of(npred + 2); tmax <- s$total_max_points; if (!is.finite(tmax) || tmax <= 0) tmax <- 1
  add_ruler(yT, 0, 1, CHARC, arrow = TRUE); add_ticks(yT, s$tp_ticks / tmax, sprintf("%.0f", s$tp_ticks), CHARC)
  add_rowlab(yT, "Total points")

  # Predicted-probability ruler (placed by total-points position). A perceptually
  # uniform Viridis bar (colour-blind- and greyscale-safe) replaces the plain rule so
  # the probability scale reads at a glance and survives black-and-white printing.
  yQ <- y_of(npred + 3); pa <- s$prob_axis; gradbar <- NULL
  if (!is.null(pa) && nrow(pa)) {
    gx0 <- min(pa$total_points) / tmax; gx1 <- max(pa$total_points) / tmax
    ng  <- 128; gx <- seq(gx0, gx1, length.out = ng + 1)
    gpal <- grDevices::hcl.colors(ng, "Viridis")
    gradbar <- data.frame(xmin = gx[-(ng + 1)], xmax = gx[-1],
                          ymin = yQ - 0.10, ymax = yQ + 0.10, fill = gpal)
    add_ticks(yQ, pa$total_points / tmax, sprintf("%.2g", pa$prob), CHARC)
  } else add_ruler(yQ, 0, 1, CHARC)
  add_rowlab(yQ, "Predicted prob.")

  bands  <- if (length(bands)) do.call(rbind, bands) else NULL
  tracks <- do.call(rbind, tracks)
  axes   <- do.call(rbind, axes); arrows <- if (length(arrows)) do.call(rbind, arrows) else NULL
  ticks  <- do.call(rbind, ticks)
  labs   <- do.call(rbind, labs);  rowlab <- do.call(rbind, rowlab)
  bullets<- if (length(bullets)) do.call(rbind, bullets) else NULL
  notes  <- if (length(notes)) do.call(rbind, notes) else NULL

  p <- ggplot()
  if (!is.null(bands))
    p <- p + geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                       fill = band_fill, colour = NA)
  if (!is.null(gradbar))
    p <- p + geom_rect(data = gradbar, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
                       colour = NA) + scale_fill_identity()
  p <- p +
    geom_segment(data = tracks, aes(x = x, xend = xe, y = y, yend = ye),
                 colour = col_track, linewidth = LW_TRACK, lineend = "round") +
    geom_segment(data = axes,  aes(x = x, xend = xe, y = y, yend = ye, colour = col), linewidth = LW_AXIS, lineend = "round") +
    geom_segment(data = ticks, aes(x = x, xend = xe, y = y, yend = ye, colour = col), linewidth = LW_TICK, lineend = "round") +
    geom_text(data = labs,   aes(x = x, y = y, label = label), colour = CHARC, family = FONT, size = 2.35) +
    geom_text(data = rowlab, aes(x = x, y = y, label = label), colour = CHARC, family = FONT,
              hjust = 1, vjust = 0.5, size = 2.8, fontface = "bold")
  if (!is.null(bullets))
    p <- p + geom_point(data = bullets, aes(x = x, y = y, colour = col), shape = 15, size = 1.7)
  if (!is.null(arrows))
    p <- p + geom_segment(data = arrows, aes(x = x, xend = xe, y = y, yend = ye), colour = CHARC,
                          linewidth = LW_AXIS, lineend = "round",
                          arrow = arrow(length = unit(1.7, "mm"), type = "closed"))
  p <- p + scale_colour_identity() +
    coord_cartesian(xlim = c(X0 - 0.02, X1 + 0.03), ylim = c(0.25, n_row + 1.2), clip = "off") +
    theme_void(base_size = 9, base_family = FONT)
  if (!is.null(notes))
    p <- p + geom_text(data = notes, aes(x = x, y = y, label = label), colour = CHARC, family = FONT,
                       hjust = 0, vjust = 0.5, size = 2.05, fontface = "italic")
  if (!is.null(header))
    p <- p +
      annotate("rect", xmin = X0, xmax = X1, ymin = n_row + 0.55, ymax = n_row + 1.05,
               fill = "#dde4ea", colour = NA) +
      annotate("text", x = X0 + 0.01, y = n_row + 0.80, hjust = 0, family = FONT,
               label = header, fontface = "bold", size = 3.1, colour = CHARC)
  p
}

#' Two-panel nomogram figure: (A) LMM gate markers, (B) immune composite + clinical.
#' Panel A = per-marker display read (apparent fit on individual markers). Panel B = the
#' run's FORMAL model (1-df immune composite + that run's clinical vars, so `include_nlr`
#' drives it) — drawn from the apparent fit but reported with its leakage-free LOO AUC.
#' Immune axes are relabelled onto their exact clinical scale (see `nomo_tick_label()`):
#' a bijection of the model scale, so nothing fitted changes.
pub_fig_nomogram <- function(nomo) {
  if (is.null(nomo)) return(NULL)
  tp  <- if (!is.null(nomo$timepoint)) nomo$timepoint else "T0"
  tpl <- switch(tp, delta = "Δ (T1−T0)", T1 = "T1", "T0")
  # Δ  -> exp(Δlogit) = EXACT fold change in the positive:negative cell ratio.
  # T0 -> the cell fraction (%) itself. Clinical vars (panel B) stay on their raw scale.
  is_delta  <- identical(tp, "delta")
  imm_scale <- if (is_delta) "fold change, Ki67+:Ki67− ratio" else "% of parent gate"
  pA  <- pub_nomo_panel(nomo$immune,
           header = sprintf("A   Immune gate markers  [%s; %s]", tpl, imm_scale))
  clab <- if (!is.null(nomo$clinical_vars)) gsub("_", "-", paste(nomo$clinical_vars, collapse = " + ")) else "clinical"
  pB  <- if (!is.null(nomo$clinical_immune))
           pub_nomo_panel(nomo$clinical_immune,
             header = sprintf("B   Immune composite [%s] + %s  =  formal model", tpl, clab)) else NULL
  if (is.null(pA) && is.null(pB)) return(NULL)
  if (is.null(pB)) return(pA)
  if (is.null(pA)) return(pB)
  nA <- length(nomo$immune$predictors) + 3
  nB <- length(nomo$clinical_immune$predictors) + 3
  # Keep every caption line under ~110 characters: at size 7.6 on PUB_W2 (190 mm) a
  # longer line is silently clipped at the panel edge rather than wrapped.
  # Panel B is the formal model → report its REPORTABLE discrimination (repeated
  # stratified k-fold). LOO is the retired, tie-artifact-prone quantity (diag_41) and is
  # used only as a fallback for publication_data rds persisted before combined_cv_auc.
  cvauc <- nomo$combined_cv_auc
  loo   <- nomo$combined_loo_auc
  auc_txt <- if (!is.null(cvauc) && is.finite(cvauc))
    sprintf(" (%d-fold CV AUC %.3f)", if (!is.null(nomo$combined_cv_k)) nomo$combined_cv_k else 10L, cvauc)
  else if (!is.null(loo) && is.finite(loo)) sprintf(" (LOO AUC %.3f)", loo)
  else ""
  foot <- paste0(
    "Sum each predictor's Points → Total points → Predicted probability of response.\n",
    "A: apparent fit, per-marker display only. B: the formal model", auc_txt,
    " — axis positions from its apparent fit.\n",
    "Continuous axes span the central 90% (5th–95th percentile) of observed values.\n",
    if (is_delta)
      paste0("Immune axes = exp(Δlogit) = fold change in the Ki67+:Ki67− (proliferating:resting) ",
             "cell ratio — exact;\ncomposite = weighted geometric mean of the three. ")
    else
      "Immune axes = cell fraction (% of parent gate) — exact; composite is an index (z).\n",
    "Clinical axes on raw scale.")
  (pA / pB + patchwork::plot_layout(heights = c(nA, nB))) +
    patchwork::plot_annotation(
      caption = foot,
      theme = ggplot2::theme(plot.caption = ggplot2::element_text(
        size = 7.6, colour = "#222222", family = "sans", hjust = 0)))   # +20% & charcoal for scaled-down legibility
}

# ── SUPP — classification performance (confusion matrix + precision/recall) ────
#' Confusion matrix + precision–recall of the combined (clinical + immune) model,
#' from its LEAKAGE-FREE LOO probabilities (`per_patient$Prob_Combined`). The
#' operating point is the PRE-SPECIFIED disease prevalence (= locked_model
#' threshold_default); 0.5 is also tabulated. The threshold is deliberately NOT
#' optimised (Youden on n=49 would be optimistic). Panel A confusion matrix,
#' Panel B precision–recall curve (combined vs clinical, no-skill = prevalence).
# ── SUPP S13 — predictive-vs-prognostic dissociation ─────────────────────────
# Panel A: every predictor on BOTH endpoints, same rank scale, so the crossover is read
# directly. Panel B: the bootstrapped CONTRASTS — the actual test, because "significant
# on one endpoint, not the other" is not a dissociation (Gelman & Stern 2006).
#' @param ds `av$survival$dissociation` from run_dissociation_analysis().
pub_fig_dissociation <- function(ds) {
  if (is.null(ds) || !is.null(ds$skipped) || is.null(ds$cells)) return(NULL)
  lab <- function(v) if (identical(v, "immune")) "Immune composite" else gsub("_", "-", v)
  A <- do.call(rbind, lapply(ds$cells, function(c1) data.frame(
    var = lab(c1$variable),
    endpoint = c("Response (AUC)", "Survival (Harrell C)"),
    est = c(c1$response$estimate, c1$survival$estimate),
    lo  = c(c1$response$ci[1], c1$survival$ci[1]),
    hi  = c(c1$response$ci[2], c1$survival$ci[2]), stringsAsFactors = FALSE)))
  A <- A[is.finite(A$est), , drop = FALSE]
  if (!nrow(A)) return(NULL)
  A$var <- factor(A$var, levels = rev(unique(A$var)))
  pA <- ggplot(A, aes(est, var, colour = endpoint)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = pub_palette[["ref"]]) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.14,
                   position = position_dodge(width = 0.45), alpha = 0.75) +
    geom_point(size = 2.6, position = position_dodge(width = 0.45)) +
    scale_colour_manual(values = c("Response (AUC)" = unname(pub_palette[["combined"]]),
                                   "Survival (Harrell C)" = unname(pub_palette[["unpen"]]))) +
    labs(x = "P(correct ordering)   ·   0.5 = no information", y = NULL, colour = NULL,
         subtitle = sprintf("Same rank scale, same %d patients (%d events)", ds$n, ds$events)) +
    theme_publication() + theme(legend.position = "top")

  rows <- list(data.frame(spec = "Immune: response - survival",
                          est = ds$contrasts$single_dissociation$estimate,
                          lo = ds$contrasts$single_dissociation$ci[1],
                          hi = ds$contrasts$single_dissociation$ci[2], stringsAsFactors = FALSE))
  for (v in names(ds$contrasts$by_comparator)) {
    cv <- ds$contrasts$by_comparator[[v]]; star <- if (identical(v, ds$comparator)) " *" else ""
    rows <- c(rows, list(
      data.frame(spec = sprintf("Crossover vs %s%s", lab(v), star),
                 est = cv$crossover$estimate, lo = cv$crossover$ci[1], hi = cv$crossover$ci[2]),
      data.frame(spec = sprintf("Response: immune - %s", lab(v)),
                 est = cv$within_response$estimate, lo = cv$within_response$ci[1], hi = cv$within_response$ci[2]),
      data.frame(spec = sprintf("Survival: immune - %s", lab(v)),
                 est = cv$within_survival$estimate, lo = cv$within_survival$ci[1], hi = cv$within_survival$ci[2])))
  }
  B <- do.call(rbind, rows); B <- B[is.finite(B$est), , drop = FALSE]
  B$spec <- factor(B$spec, levels = rev(B$spec))
  B$excl <- is.finite(B$lo) & is.finite(B$hi) & (B$lo > 0 | B$hi < 0)
  pB <- ggplot(B, aes(est, spec, colour = excl)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = pub_palette[["ref"]]) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.16) +
    geom_point(size = 2.6) +
    scale_colour_manual(values = c(`TRUE` = unname(pub_palette[["ok"]]),
                                   `FALSE` = unname(pub_palette[["fragile"]])), guide = "none") +
    labs(x = "Difference in P(correct ordering)  [95% bootstrap CI]", y = NULL,
         subtitle = if (is.null(ds$comparator)) "* = pre-specified comparator" else
           sprintf("* = pre-specified comparator (%s)", lab(ds$comparator))) +
    theme_publication()
  pub_tag(pA / pB)
}

pub_fig_classification <- function(av) {
  if (is.null(av) || is.null(av$per_patient) ||
      !"Prob_Combined" %in% names(av$per_patient)) return(NULL)
  pp  <- av$per_patient
  pos <- av$positive_label; neg <- av$negative_label
  ok  <- is.finite(pp$Prob_Combined) & !is.na(pp$True_Group)
  y   <- as.integer(pp$True_Group[ok] == pos); p <- pp$Prob_Combined[ok]
  if (length(unique(y)) < 2) return(NULL)
  n <- length(y); prev <- mean(y)
  pc <- if ("Prob_Clinical" %in% names(pp)) pp$Prob_Clinical[ok] else NULL

  metr <- function(pr, thr) {
    pred <- as.integer(pr >= thr)
    TP <- sum(pred & y == 1); FP <- sum(pred & y == 0)
    FN <- sum(!pred & y == 1); TN <- sum(!pred & y == 0)
    rec  <- TP / (TP + FN); prc <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
    spec <- TN / (TN + FP); npv <- if ((TN + FN) > 0) TN / (TN + FN) else NA_real_
    f1   <- if (is.finite(prc) && (prc + rec) > 0) 2 * prc * rec / (prc + rec) else NA_real_
    list(TP = TP, FP = FP, FN = FN, TN = TN, rec = rec, prc = prc, spec = spec,
         npv = npv, f1 = f1, acc = (TP + TN) / n, bacc = (rec + spec) / 2)
  }
  pr_curve <- function(pr) {
    thr <- sort(unique(pr), decreasing = TRUE)
    d <- do.call(rbind, lapply(thr, function(t) { m <- metr(pr, t); data.frame(recall = m$rec, precision = m$prc) }))
    d <- d[is.finite(d$precision) & is.finite(d$recall), , drop = FALSE]
    d <- d[order(d$recall), , drop = FALSE]
    if (nrow(d)) d <- rbind(data.frame(recall = 0, precision = d$precision[1]), d)
    d
  }
  auprc <- function(d) if (nrow(d) < 2) NA_real_ else
    sum(diff(d$recall) * (head(d$precision, -1) + tail(d$precision, -1)) / 2)

  mP <- metr(p, prev); mH <- metr(p, 0.5)
  prc_comb <- pr_curve(p); ap_comb <- auprc(prc_comb)
  prc_clin <- if (!is.null(pc)) pr_curve(pc) else NULL
  ap_clin  <- if (!is.null(prc_clin)) auprc(prc_clin) else NA_real_

  # Panel A — confusion matrix at the prevalence operating point (LOO)
  posd <- pub_relabel_group(pos); negd <- pub_relabel_group(neg)
  cm <- data.frame(
    Actual = factor(c(posd, posd, negd, negd), levels = c(posd, negd)),
    Pred   = factor(c(posd, negd, posd, negd), levels = c(negd, posd)),
    n      = c(mP$TP, mP$FN, mP$FP, mP$TN),
    kind   = c("correct", "error", "error", "correct"))
  pA <- ggplot(cm, aes(Pred, Actual)) +
    geom_tile(aes(fill = kind), colour = "white", linewidth = 2) +
    geom_text(aes(label = n), size = 6, fontface = "bold", colour = "black") +
    scale_fill_manual(values = c(correct = "#A6DBC9", error = "grey85"), guide = "none") +
    scale_x_discrete(position = "top") +
    labs(x = "Predicted", y = "Actual") +
    coord_equal() + theme_publication() +
    theme(panel.grid = element_blank(), axis.text = element_text(size = 10))

  # Panel B — precision–recall curve (combined vs clinical; no-skill = prevalence)
  pB <- ggplot(prc_comb, aes(recall, precision)) +
    geom_hline(yintercept = prev, linetype = "dotted", colour = "grey50")
  if (!is.null(prc_clin))
    pB <- pB + geom_path(data = prc_clin, aes(recall, precision),
                         colour = pub_palette[["clinical"]], linewidth = 0.7, linetype = "dashed")
  pB <- pB +
    geom_path(colour = pub_palette[["combined"]], linewidth = 0.9) +
    annotate("text", x = 0.02, y = prev - 0.05, hjust = 0, size = 2.3, colour = "grey45",
             label = sprintf("no-skill = %.2f", prev)) +
    annotate("text", x = 0.98, y = 0.05, hjust = 1, size = 2.5, colour = pub_palette[["combined"]],
             label = sprintf("AUPRC %.2f", ap_comb)) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "Recall (sensitivity)", y = "Precision (PPV)") +
    theme_publication()

  # Panel C — metrics table (monospace so columns align) + provenance line
  hdr  <- sprintf("%-20s %5s %5s %5s %5s %5s %5s", "Operating point", "Sens", "Prec", "Spec", "NPV", "F1", "Acc")
  fmt  <- function(lab, m) sprintf("%-20s %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f",
                                   lab, m$rec, m$prc, m$spec, m$npv, m$f1, m$acc)
  tbl  <- paste(hdr, fmt(sprintf("prevalence (%.2f)", prev), mP), fmt("0.50", mH), sep = "\n")
  note <- sprintf("Combined (clinical+immune) model · leakage-free LOO probabilities (n=%d).\nAUPRC %.2f vs clinical %.2f (no-skill = prevalence %.2f) · threshold pre-specified, not optimised.",
                  n, ap_comb, ap_clin, prev)
  pC <- ggplot() +
    annotate("text", x = 0, y = 1, hjust = 0, vjust = 1, family = "mono", size = 2.6, label = tbl) +
    annotate("text", x = 0, y = 0.06, hjust = 0, vjust = 0, size = 2.2, colour = "grey40", label = note) +
    xlim(0, 1) + ylim(0, 1) + theme_void()

  ((pA | pB) / pC + patchwork::plot_layout(heights = c(3.2, 1.15))) +
    patchwork::plot_annotation(tag_levels = list(c("A", "B", ""))) &
    ggplot2::theme(plot.tag = ggplot2::element_text(size = 11, face = "bold"))
}

# ── SUPP — on-treatment fold-change of the gate markers (raw, clinician-facing) ─
#' Descriptive fold-change view of the immune-gate dynamic, computed DIRECTLY from
#' the raw cell fractions (`gate_decomp$per_patient_values$Value_pct`) — NOT derived
#' from the logit/z model, so it carries no transformation-equivalence assumption.
#' Panel A: per-marker log2 fold-change (T1/T0) by response group, with a fold-change
#' secondary axis and a reference line at "no change" (×1). Panel B: median %
#' change per group. The formal inferential test remains the Step-04 LMM
#' Time×Group interaction; the between-group p-values here are descriptive only.
pub_fig_foldchange <- function(gate_decomp, resp_label = NULL) {
  ppv <- if (!is.null(gate_decomp)) gate_decomp$per_patient_values else NULL
  if (is.null(ppv) || !nrow(ppv) || !"Value_pct" %in% names(ppv)) return(NULL)
  d  <- ppv[is.finite(ppv$Value_pct), c("Patient_ID", "Marker", "Timepoint", "Group", "Value_pct")]
  t0 <- d[d$Timepoint == "T0", c("Patient_ID", "Marker", "Group", "Value_pct")]; names(t0)[4] <- "pct0"
  t1 <- d[d$Timepoint == "T1", c("Patient_ID", "Marker", "Value_pct")];         names(t1)[3] <- "pct1"
  w  <- merge(t0, t1, by = c("Patient_ID", "Marker"))
  w  <- w[is.finite(w$pct0) & is.finite(w$pct1) & w$pct0 > 0, , drop = FALSE]
  if (!nrow(w)) return(NULL)
  w$log2FC <- log2(w$pct1 / w$pct0)
  w$pctchg <- 100 * (w$pct1 - w$pct0) / w$pct0

  # group colours/order matched to Fig 1 biology: responder = blue, non-resp = vermillion
  w$Group   <- pub_relabel_group(w$Group)
  resp_disp <- if (!is.null(resp_label)) pub_relabel_group(resp_label) else NULL
  grps <- unique(as.character(w$Group))
  ord  <- if (!is.null(resp_disp) && resp_disp %in% grps) c(resp_disp, setdiff(grps, resp_disp)) else sort(grps)
  w$Group  <- factor(w$Group, levels = ord)
  grp_cols <- setNames(c(pub_palette[["combined"]], pub_palette[["unpen"]])[seq_along(ord)], ord)
  w$Marker <- factor(w$Marker, levels = unique(w$Marker))

  # descriptive between-group Wilcoxon on log2FC (per marker)
  pl <- lapply(levels(w$Marker), function(mk) {
    s <- w[w$Marker == mk, ]
    p <- tryCatch(suppressWarnings(wilcox.test(log2FC ~ Group, data = s)$p.value), error = function(e) NA_real_)
    data.frame(Marker = mk, y = max(s$log2FC, na.rm = TRUE) + 0.4,
               label = if (is.finite(p)) sprintf("p=%.3f", p) else "")
  })
  plab <- do.call(rbind, pl)

  # Panel A — log2 fold-change distribution
  pA <- ggplot(w, aes(Marker, log2FC, colour = Group)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
    geom_boxplot(fill = NA, width = 0.6, outlier.shape = NA,
                 position = position_dodge(0.7), linewidth = 0.5) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7),
               size = 0.7, alpha = 0.45) +
    geom_text(data = plab, aes(Marker, y, label = label), inherit.aes = FALSE,
              size = 2.2, colour = "grey35") +
    scale_colour_manual(values = grp_cols) +
    scale_y_continuous(sec.axis = sec_axis(~ 2^., name = "fold-change (T1/T0)",
                                           breaks = c(0.125, 0.25, 0.5, 1, 2, 4))) +
    guides(colour = "none") +                       # shared legend comes from panel B
    labs(x = NULL, y = expression(log[2]~"fold-change (T1/T0)"),
         caption = "dashed line = no change (×1)") +
    theme_publication() +
    theme(plot.caption = element_text(size = 7, colour = "grey45", hjust = 0))

  # Panel B — median % change per group
  med <- aggregate(pctchg ~ Marker + Group, data = w, FUN = median)
  pB <- ggplot(med, aes(pctchg, Marker)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
    geom_line(aes(group = Marker), colour = "grey75", linewidth = 0.6) +
    geom_point(aes(colour = Group), size = 2.8) +
    geom_text(aes(colour = Group, label = sprintf("%+.0f%%", pctchg)),
              vjust = -1.1, size = 2.3, show.legend = FALSE) +
    scale_colour_manual(values = grp_cols) +
    scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
    labs(x = "median change T0→T1 (%)", y = NULL) +
    theme_publication()

  pub_tag((pA | pB) + patchwork::plot_layout(widths = c(1.5, 1), guides = "collect") &
          theme(legend.position = "bottom"))
}

pub_render_all <- function(objs, out_dir, project_name) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  pf <- function(n) file.path(out_dir, sprintf("%s_%s.pdf", n, project_name))
  av <- objs$clin_addval
  # MAIN (3) — Fig1 biology (cells + forest) / Fig2 added-value / Fig3 PD-L1 context.
  # SUPP (S1–S9) — S1 CONSORT, S2 baseline-inv, S3 specificity, S4 standalone,
  #        S5 coupling, S6 robustness, S7 calibration, S8 gate-signal decomposition,
  #        S9 nested-selection validation. Defense figures kept out of main.
  if (!is.null(objs$lmm_boot) || !is.null(objs$gate_decomp))
    save_pub_figure(pub_fig_biology(objs$gate_decomp, objs$lmm_boot$boot, objs$lmm_boot$obs,
                                    resp_label = objs$positive_label),
                    pf("Figure1_Biology"), PUB_W2, 120)
  if (!is.null(av)) {
    save_pub_figure(pub_fig_added_value(av),         pf("Figure2_AddedValue"),      PUB_W2, 120)
    # SUPP (S1–S9)
    if (!is.null(objs$consort))
      save_pub_figure(pub_fig_consort(objs$consort, terminal = objs$consort_terminal),
                                                     pf("FigureS1_CONSORT"),            PUB_W2, 140)
    save_pub_figure(pub_fig_baseline_invariance(av), pf("FigureS2_BaselineInvariance"), PUB_W2, 82)
    save_pub_figure(pub_fig_specificity_null(av),    pf("FigureS3_SpecificityNull"),    PUB_W1, 95)
    if (!is.null(av$dynamics_baseline_coupling))
      save_pub_figure(pub_fig_coupling(av),          pf("FigureS5_Coupling"),           PUB_W1, 95)
    save_pub_figure(pub_fig_robustness(av),          pf("FigureS6_Robustness"),         PUB_W2, 95)
    save_pub_figure(pub_fig_calibration_idi(av),     pf("FigureS7_Calibration_IDI"),    PUB_W2, 120)
    if (!is.null(av$nomogram))
      save_pub_figure(pub_fig_nomogram(av$nomogram), pf("FigureS10_Nomogram"),          PUB_W2, 172)
    save_pub_figure(pub_fig_classification(av),      pf("FigureS11_Classification"),    PUB_W2, 120)
    # S13 dissociation — a normal live-object figure, so unlike S9 it needs no
    # OUT_OF_PIPELINE protection in make_manuscript_figures.R.
    if (!is.null(av$survival$dissociation) && is.null(av$survival$dissociation$skipped))
      save_pub_figure(pub_fig_dissociation(av$survival$dissociation),
                                                     pf("FigureS13_Dissociation"),      PUB_W2, 150)
  }
  if (!is.null(objs$stratified_result))
    save_pub_figure(pub_fig_pdl1_context(objs$stratified_result), pf("Figure3_PDL1_context"), PUB_W2, 115)
  if (!is.null(objs$df_preds))
    save_pub_figure(pub_fig_standalone(objs$df_preds, objs$positive_label, perm = objs$perm),
                    pf("FigureS4_StandaloneClassifier"), PUB_W1, 120)
  # S8 gate-signal decomposition (demoted from main). S9 is NOT rendered here: the
  # selection-aware nulls are one-time computations, not per-run live objects — see
  # manuscript/figures/render_selection_aware.R. objs$nested_val is still computed
  # and lives in the JSON as the in-pipeline anti-circularity record.
  if (!is.null(objs$gate_decomp))
    save_pub_figure(pub_fig_gate_signal(objs$gate_decomp),         pf("FigureS8_GateSignal"),     PUB_W2, 110)
  # S12 raw fold-change of the gate markers (clinician-facing biology; computed from raw %)
  if (!is.null(objs$gate_decomp))
    save_pub_figure(pub_fig_foldchange(objs$gate_decomp, resp_label = objs$positive_label),
                    pf("FigureS12_FoldChange"), PUB_W2, 100)
  # Secondary-timepoint added-value artifacts (e.g. T0 reference when delta is primary,
  # or delta when T0 is primary). Builders are timepoint-agnostic → render per node so
  # the publication/ dir always reflects every analyzed timepoint.
  if (!is.null(objs$clin_addval_secondary)) {
    for (tp in names(objs$clin_addval_secondary)) {
      av_tp <- objs$clin_addval_secondary[[tp]]
      if (is.null(av_tp)) next
      save_pub_figure(pub_fig_added_value(av_tp),           pf(paste0("Figure2_AddedValue_", tp)), PUB_W2, 120)
      if (!is.null(av_tp$nomogram))
        save_pub_figure(pub_fig_nomogram(av_tp$nomogram),   pf(paste0("FigureS10_Nomogram_", tp)), PUB_W2, 172)
      save_pub_figure(pub_fig_classification(av_tp),        pf(paste0("FigureS11_Classification_", tp)), PUB_W2, 120)
    }
  }
  invisible(out_dir)
}

# ── helper: CONSORT inputs from the live cohort ledger ────────────────────────
# Replaces pub_consort_counts(), which reconstructed the flow as
#   raw_n = av$n + nrow(qc_dropped)
# i.e. it added exclusions counted from the ANALYTIC cohort to the PAIRED-subset n.
# Those are different denominators, and the result was a starting cohort that never
# existed (69 for a 94-row file whose analytic cohort is 82). The ledger is recorded
# at each filtering point, so there is nothing left to reconstruct: this helper only
# assembles the TERMINAL annotations, which are properties of the final cohort and
# explicitly not exclusions.
pub_consort_terminal <- function(av) {
  tl <- list()
  if (!is.null(av$n_complete_case))
    tl[[sprintf("Complete-case %s",
                gsub("_", "-", paste(names(av$clinical_vars), collapse = "+")))]] <-
      sprintf("n = %d (NAs median-imputed in-fold; nobody excluded)", av$n_complete_case)
  if (!is.null(av$survival$OS$n))
    tl[["Survival follow-up"]] <- sprintf("n = %d (%d OS / %s PFS events)",
      av$survival$OS$n, av$survival$OS$events,
      if (is.null(av$survival$PFS$events)) "-" else av$survival$PFS$events)
  if (!length(tl)) NULL else tl
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
