# R/modules_posthoc_grouping.R
# ==============================================================================
# RESPONSE GROUPING SENSITIVITY (step 07, grouping/; formerly diagnostics/diag_50)
# PR / SD / PD instead of PR vs SD/PD, on the step-04 LMM frame:
#   0. assert the two-group fit reproduces the step-04 LMM p and beta (every marker)
#   1. three-group LMM: omnibus Time x Group3 (F, 2 df, BH over markers) + pairwise
#      contrasts of the T0->T1 change (PR-SD, PR-PD, SD-PD), Holm within marker
#   2. formal grouping tests (ML, baseline always 3 groups):
#      (a) LRT change PR | SD+PD vs free (H0: SD and PD change alike)
#          LRT change PR+SD | PD vs free (H0: PR and SD change alike)
#      (c) PR | SD+PD vs PR+SD | PD: same k, compared by AIC
#      plus the disease-control gate (PR+SD vs PD) re-selected from scratch
#   3. 1-2 adjusted for the step-04 sensitivity covariates (PD-L1 + therapy, T0 file)
# Response sub-codes and covariates come from the RUN'S OWN input files.
# Supplementary consistency check: selection-biased for the gate markers.
# Dependencies: lmerTest, lme4, emmeans
# ==============================================================================

#' @title Response-grouping sensitivity of the LMM frame
#' @param DATA_LONG Step-01 longitudinal object.
#' @param t0_raw,t1_raw The run's T0 / T1 input sheets (data frames).
#' @param target_col Best-response column (config clinical$target_column).
#' @param subcodes Named codes, e.g. c(PR = 2, SD = 3, PD = 4).
#' @param responder_label Group label of PR in the step-01 metadata (e.g. "RP").
#' @param covs Covariates of the adjusted panel (config clinical$sensitivity_covariates).
#' @param step04_full Step-04 JSON full_results (Marker, P_Value_Interaction,
#'   Estimate_Interaction, N_Observations; 4 dp) for the reproduction gate.
#' @param report_markers Markers detailed in the log (the step-06 gate).
#' @param say sprintf-style logger.
#' @return List: unadjusted / adjusted per-marker tables, frame counts.
run_response_grouping_sensitivity <- function(DATA_LONG, t0_raw, t1_raw, target_col,
                                              subcodes = c(PR = 2, SD = 3, PD = 4),
                                              responder_label = "RP", covs = character(0),
                                              step04_full, report_markers, say = message) {
  for (pkg in c("lmerTest", "lme4", "emmeans"))
    if (!requireNamespace(pkg, quietly = TRUE)) stop(sprintf("[FATAL] grouping sensitivity needs package '%s'", pkg))
  emmeans::emm_options(lmer.df = "satterthwaite", lmerTest.limit = 1e5)
  sc <- as.character(unlist(subcodes)[c("PR", "SD", "PD")])
  if (anyNA(sc)) stop("[FATAL] post_hoc.grouping.subcodes must name PR, SD and PD")

  # ── Data ────────────────────────────────────────────────────────────────────
  D  <- DATA_LONG
  mk_all <- D$hybrid_markers
  df <- cbind(D$metadata, as.data.frame(D$hybrid_data_raw[, mk_all]))
  df$Timepoint  <- factor(df$Timepoint, levels = c("T0", "T1"))
  df$Patient_ID <- as.character(df$Patient_ID)

  code <- setNames(as.character(t0_raw[[target_col]]), t0_raw$Patient_ID)
  c1   <- setNames(as.character(t1_raw[[target_col]]), t1_raw$Patient_ID)
  both <- intersect(names(code), names(c1))
  n_mismatch <- sum(code[both] != c1[both])
  say("T0/T1 response code mismatches among patients in both files: %d", n_mismatch)
  if (n_mismatch != 0) stop("[FATAL] grouping: T0 and T1 response codes disagree")
  code <- c(code, c1[setdiff(names(c1), names(code))])

  df$Code <- unname(code[df$Patient_ID])
  if (!all(df$Code %in% sc)) stop("[FATAL] grouping: a response code of the frame is not PR/SD/PD")
  if (!all((df$Code == sc[1]) == (df$Group == responder_label)))
    stop("[FATAL] grouping: PR sub-code does not match the two-group responder label")
  df$G3  <- factor(setNames(c("PR", "SD", "PD"), sc)[df$Code], levels = c("PD", "SD", "PR"))
  df$ORR <- factor(ifelse(df$G3 == "PR", "PR", "SD_PD"), levels = c("SD_PD", "PR"))
  df$DCR <- factor(ifelse(df$G3 == "PD", "PD", "DC"), levels = c("PD", "DC"))
  df$t1    <- as.numeric(df$Timepoint == "T1")
  df$t1_PR <- df$t1 * (df$G3 == "PR")
  df$t1_SD <- df$t1 * (df$G3 == "SD")
  df$t1_PD <- df$t1 * (df$G3 == "PD")

  if (length(covs)) {
    miss <- setdiff(covs, names(t0_raw))
    if (length(miss)) stop(sprintf("[FATAL] grouping: covariate(s) not in the T0 input: %s", paste(miss, collapse = ", ")))
    cov_t0 <- t0_raw[, c("Patient_ID", covs)]
    df <- merge(df, cov_t0, by = "Patient_ID", all.x = TRUE, sort = FALSE)
  }

  tp <- split(as.character(df$Timepoint), df$Patient_ID)
  pat <- data.frame(Patient_ID = names(tp),
                    paired = vapply(tp, function(v) all(c("T0", "T1") %in% v), logical(1)),
                    only_T0 = vapply(tp, function(v) identical(v, "T0"), logical(1)),
                    only_T1 = vapply(tp, function(v) identical(v, "T1"), logical(1)))
  pat$G3 <- df$G3[match(pat$Patient_ID, df$Patient_ID)]
  say("frame: %d patients, %d observations, %d markers", nrow(pat), nrow(df), length(mk_all))
  counts <- list()
  for (g in levels(df$G3)) {
    p <- pat[pat$G3 == g, ]
    say("  %s: %d patients (%d paired, %d T0-only, %d T1-only), %d obs",
        g, nrow(p), sum(p$paired), sum(p$only_T0), sum(p$only_T1), sum(df$G3 == g))
    counts[[g]] <- c(patients = nrow(p), paired = sum(p$paired), only_T0 = sum(p$only_T0),
                     only_T1 = sum(p$only_T1), obs = sum(df$G3 == g))
  }

  # ── Model helpers (mirror fit_feature_lmm: Value / sd over the full marker) ─
  cov_terms <- function(covs) if (length(covs)) paste("+", paste(sprintf("`%s`", covs), collapse = " + ")) else ""

  marker_frame <- function(mk, covs) {
    v  <- as.numeric(df[[mk]])
    sdv <- sd(v, na.rm = TRUE)
    if (is.na(sdv) || sdv == 0) sdv <- 1
    d <- data.frame(Value = v / sdv, ID = factor(df$Patient_ID),
                    Time = df$Timepoint, ORR = df$ORR, DCR = df$DCR, G3 = df$G3,
                    t1 = df$t1, t1_PR = df$t1_PR, t1_SD = df$t1_SD, t1_PD = df$t1_PD)
    for (cv in covs) d[[cv]] <- df[[cv]]
    d <- d[complete.cases(d), ]
    list(d = d, sd = sdv)
  }

  fit <- function(f, d, reml) suppressMessages(suppressWarnings(
    lmerTest::lmer(as.formula(f), data = d, REML = reml,
                   control = lme4::lmerControl(calc.derivs = FALSE))))

  analyse_marker <- function(mk, covs = character(0)) {
    mf <- marker_frame(mk, covs); d <- mf$d; s <- mf$sd
    ct <- cov_terms(covs)

    # two-group REML (= step 04)
    m2 <- fit(paste("Value ~ Time * ORR", ct, "+ (1 | ID)"), d, TRUE)
    co2 <- summary(m2)$coefficients
    r2 <- grep("Time.*:ORR", rownames(co2))

    # two-group REML with disease control (PR+SD vs PD): the alternative gate
    m2d <- fit(paste("Value ~ Time * DCR", ct, "+ (1 | ID)"), d, TRUE)
    co2d <- summary(m2d)$coefficients
    r2d <- grep("Time.*:DCR", rownames(co2d))

    # three-group REML: omnibus + change contrasts
    m3 <- fit(paste("Value ~ Time * G3", ct, "+ (1 | ID)"), d, TRUE)
    av <- anova(m3)
    ch <- emmeans::contrast(emmeans::emmeans(m3, ~ Time | G3), "revpairwise")   # T1 - T0 per group
    chs <- summary(ch, infer = TRUE)
    pw <- summary(emmeans::contrast(ch, "revpairwise", by = NULL, adjust = "none"), infer = TRUE)
    pw$contrast <- gsub("\\(T1 - T0 |\\)", "", sub("^\\(T1 - T0 (\\w+)\\) - \\(T1 - T0 (\\w+)\\)$", "\\1-\\2", pw$contrast))
    pw$p_holm <- p.adjust(pw$p.value, "holm")
    grab <- function(lbl, col) { i <- which(pw$contrast == lbl); if (length(i) == 1) pw[[col]][i] else NA_real_ }

    # ML fits for the grouping tests (baseline 3 groups in every model)
    base <- paste("Value ~ G3 + t1", ct)
    mF   <- fit(paste(base, "+ t1_PR + t1_SD + (1 | ID)"), d, FALSE)   # change free in 3 groups
    mORR <- fit(paste(base, "+ t1_PR + (1 | ID)"), d, FALSE)           # SD change = PD change
    mDCR <- fit(paste(base, "+ t1_PD + (1 | ID)"), d, FALSE)           # PR change = SD change
    llF <- as.numeric(logLik(mF)); llO <- as.numeric(logLik(mORR)); llD <- as.numeric(logLik(mDCR))
    stopifnot(attr(logLik(mORR), "df") == attr(logLik(mDCR), "df"))
    lrt_sd_pd <- 2 * (llF - llO); lrt_pr_sd <- 2 * (llF - llD)
    dAIC <- AIC(mDCR) - AIC(mORR)                                       # > 0: PR | SD+PD fits better

    data.frame(
      Marker = mk, N_obs = nrow(d), N_pat = length(unique(d$ID)),
      p_2group = co2[r2, "Pr(>|t|)"], beta_2group = co2[r2, "Estimate"] * s,
      p_DCR = co2d[r2d, "Pr(>|t|)"], beta_DCR = co2d[r2d, "Estimate"] * s,
      F_omnibus = av["Time:G3", "F value"], df_den_omnibus = av["Time:G3", "DenDF"],
      p_omnibus = av["Time:G3", "Pr(>F)"],
      change_PR = chs$estimate[chs$G3 == "PR"] * s,
      change_SD = chs$estimate[chs$G3 == "SD"] * s,
      change_PD = chs$estimate[chs$G3 == "PD"] * s,
      PR_SD = grab("PR-SD", "estimate") * s, PR_SD_lo = grab("PR-SD", "lower.CL") * s,
      PR_SD_hi = grab("PR-SD", "upper.CL") * s, p_PR_SD = grab("PR-SD", "p.value"),
      PR_PD = grab("PR-PD", "estimate") * s, PR_PD_lo = grab("PR-PD", "lower.CL") * s,
      PR_PD_hi = grab("PR-PD", "upper.CL") * s, p_PR_PD = grab("PR-PD", "p.value"),
      SD_PD = grab("SD-PD", "estimate") * s, SD_PD_lo = grab("SD-PD", "lower.CL") * s,
      SD_PD_hi = grab("SD-PD", "upper.CL") * s, p_SD_PD = grab("SD-PD", "p.value"),
      p_holm_PR_SD = grab("PR-SD", "p_holm"), p_holm_PR_PD = grab("PR-PD", "p_holm"),
      p_holm_SD_PD = grab("SD-PD", "p_holm"),
      LRT_SDeqPD = lrt_sd_pd, p_LRT_SDeqPD = pchisq(lrt_sd_pd, 1, lower.tail = FALSE),
      LRT_PReqSD = lrt_pr_sd, p_LRT_PReqSD = pchisq(lrt_pr_sd, 1, lower.tail = FALSE),
      dAIC_DCR_minus_ORR = dAIC, evidence_ratio_ORR = exp(dAIC / 2),
      singular_3group = lme4::isSingular(m3),
      stringsAsFactors = FALSE)
  }

  run_panel <- function(covs) {
    res <- do.call(rbind, lapply(mk_all, analyse_marker, covs = covs))
    res$FDR_2group  <- p.adjust(res$p_2group, "BH")
    res$FDR_DCR     <- p.adjust(res$p_DCR, "BH")
    res$FDR_omnibus <- p.adjust(res$p_omnibus, "BH")
    res[order(res$p_omnibus), ]
  }

  # ── 0. Reproduction gate ───────────────────────────────────────────────────
  ref <- step04_full[, c("Marker", "P_Value_Interaction", "Estimate_Interaction", "N_Observations")]
  say("\n== Unadjusted ==")
  res <- run_panel(character(0))
  chk <- merge(ref, res[, c("Marker", "p_2group", "beta_2group", "N_obs")], by = "Marker")
  # the step-04 JSON stores p and beta rounded to 4 decimals
  same_p    <- all(round(chk$p_2group, 4) == chk$P_Value_Interaction)
  same_beta <- all(round(chk$beta_2group, 4) == chk$Estimate_Interaction)
  same_n    <- all(chk$N_obs == chk$N_Observations)
  say("[0] reproduction vs step-04 (4 dp as stored): %d/%d markers, p equal: %s, beta equal: %s, N_obs equal: %s",
      nrow(chk), length(mk_all), same_p, same_beta, same_n)
  if (!(nrow(chk) == length(mk_all) && same_p && same_beta && same_n))
    stop("[FATAL] grouping: the two-group fit does not reproduce step 04")

  report <- function(res, tag) {
    say("\n-- [%s] FDR markers of the primary two-group gate --", tag)
    for (mk in report_markers) {
      r <- res[res$Marker == mk, ]
      say("%s (n obs %d, pts %d)", mk, r$N_obs, r$N_pat)
      say("  change T1-T0 (logit): PR %+.3f | SD %+.3f | PD %+.3f", r$change_PR, r$change_SD, r$change_PD)
      say("  omnibus Time x Group3: F = %.2f (2, %.1f), p = %.4g, FDR(%d) = %.4g",
          r$F_omnibus, r$df_den_omnibus, r$p_omnibus, nrow(res), r$FDR_omnibus)
      say("  PR-SD %+.3f [%+.3f, %+.3f] p %.4g (Holm %.4g)", r$PR_SD, r$PR_SD_lo, r$PR_SD_hi, r$p_PR_SD, r$p_holm_PR_SD)
      say("  PR-PD %+.3f [%+.3f, %+.3f] p %.4g (Holm %.4g)", r$PR_PD, r$PR_PD_lo, r$PR_PD_hi, r$p_PR_PD, r$p_holm_PR_PD)
      say("  SD-PD %+.3f [%+.3f, %+.3f] p %.4g (Holm %.4g)", r$SD_PD, r$SD_PD_lo, r$SD_PD_hi, r$p_SD_PD, r$p_holm_SD_PD)
      say("  (a) LRT SD=PD: chi2 %.2f, p %.4g | LRT PR=SD: chi2 %.2f, p %.4g",
          r$LRT_SDeqPD, r$p_LRT_SDeqPD, r$LRT_PReqSD, r$p_LRT_PReqSD)
      say("  (c) AIC(PR+SD|PD) - AIC(PR|SD+PD) = %+.2f (evidence ratio for PR|SD+PD %.1f)",
          r$dAIC_DCR_minus_ORR, r$evidence_ratio_ORR)
    }
    say("\n-- [%s] all %d markers --", tag, nrow(res))
    say("omnibus FDR < 0.05: %s", paste(res$Marker[res$FDR_omnibus < 0.05], collapse = ", "))
    say("omnibus nominal p < 0.05: %d (%s)", sum(res$p_omnibus < 0.05),
        paste(res$Marker[res$p_omnibus < 0.05], collapse = ", "))
    say("two-group FDR < 0.05: %s", paste(res$Marker[res$FDR_2group < 0.05], collapse = ", "))
    say("(d) markers where PR|SD+PD fits better than PR+SD|PD (dAIC > 0): %d / %d",
        sum(res$dAIC_DCR_minus_ORR > 0), nrow(res))
    say("    of which dAIC > 2: %d | PR+SD|PD better by > 2: %d",
        sum(res$dAIC_DCR_minus_ORR > 2), sum(res$dAIC_DCR_minus_ORR < -2))
    say("    LRT SD=PD rejected (p < 0.05): %d | LRT PR=SD rejected: %d",
        sum(res$p_LRT_SDeqPD < 0.05), sum(res$p_LRT_PReqSD < 0.05))
    say("disease-control gate (PR+SD vs PD) FDR < 0.05: %s",
        if (any(res$FDR_DCR < 0.05)) paste(res$Marker[res$FDR_DCR < 0.05], collapse = ", ") else "none")
    say("    min FDR %.3f | nominal p < 0.05: %d (objective-response gate: %d)",
        min(res$FDR_DCR), sum(res$p_DCR < 0.05), sum(res$p_2group < 0.05))
    top <- head(res[order(res$p_DCR), ], 5)
    say("    top-5 disease-control markers (p, FDR, beta | objective-response p | dAIC):")
    for (i in seq_len(nrow(top)))
      say("      %-10s p %.4f FDR %.3f beta %+.2f | ORR p %.4f | dAIC %+.2f",
          top$Marker[i], top$p_DCR[i], top$FDR_DCR[i], top$beta_DCR[i], top$p_2group[i],
          top$dAIC_DCR_minus_ORR[i])
    nom <- res[res$p_omnibus < 0.05, ]
    if (nrow(nom)) {
      say("    among omnibus p < 0.05 markers (%d): PR|SD+PD better in %d", nrow(nom), sum(nom$dAIC_DCR_minus_ORR > 0))
      for (i in seq_len(nrow(nom)))
        say("      %-10s omni p %.4f | PR %+.2f SD %+.2f PD %+.2f | dAIC %+.2f | p SD=PD %.3f, p PR=SD %.3f",
            nom$Marker[i], nom$p_omnibus[i], nom$change_PR[i], nom$change_SD[i], nom$change_PD[i],
            nom$dAIC_DCR_minus_ORR[i], nom$p_LRT_SDeqPD[i], nom$p_LRT_PReqSD[i])
    }
  }
  report(res, "unadjusted")

  res_adj <- NULL
  if (length(covs)) {
    say("\n== Adjusted for %s (numeric, T0 file) ==", paste(covs, collapse = " + "))
    res_adj <- run_panel(covs)
    report(res_adj, "adjusted")
  }
  say("\nR %s | lme4 %s | lmerTest %s | emmeans %s", getRversion(), packageVersion("lme4"),
      packageVersion("lmerTest"), packageVersion("emmeans"))
  list(unadjusted = res, adjusted = res_adj, n_markers = length(mk_all),
       n_patients = nrow(pat), n_obs = nrow(df), counts = counts, covariates = covs)
}
