# R/modules_therapy_sensitivity.R
# ==============================================================================
# THERAPY SENSITIVITY MODULE — is the response-associated T0->T1 change of the
# LMM gate markers explained by concomitant chemotherapy?
# Every patient receives an immune checkpoint inhibitor (ICI); the comparison is
# ICI monotherapy vs chemo-immunotherapy. Gate markers are FIXED (selected on the
# whole cohort), so every p here is a sensitivity / exploratory value.
#
# Pure functions: data frames in, R objects out. No paths, no file I/O — the
# driver (src/07_post_hoc.R, therapy/) owns all paths. Sensitivity analysis only:
# its numbers never feed the step-04 / step-06 results.
#
# Tests per marker (numbering as in diagnostics/diag_56, which this replaces):
#   [1] LMM V ~ Time*Group + Time*Rx + (1|ID); heterogeneity Time:Group:Rx
#       (Wald + 1-df LRT, ML)
#   [2] per-stratum LMM; paired delta exact Wilcoxon and ANCOVA (delta ~ Group + T0)
#   [3] paired delta: Group permuted within Rx (B), OLS delta ~ Group + Rx,
#       ANCOVA delta ~ Group + T0 + Rx, heterogeneity Group x Rx with T0 (F test)
#   [4] [1] and [3] without patients with radiotherapy between T0 and T1
#   [5] baseline-adjusted mean change by Group x Rx: lm(delta ~ Group*Rx + T0)
#       predicted at the mean T0 (the figure)
# Values are the step-04 scale (hybrid_data_raw, logit); LMMs are fitted on
# Value / sd(Value) and back-scaled, as in fit_feature_lmm().
# Dependencies: lmerTest, lme4, ggplot2 (+ modules_pub_style.R for the figure)
# ==============================================================================

RX_LEVELS <- c("ChemoIO", "IOmono")
RX_LABELS <- c(ChemoIO = "Chemo-immunotherapy", IOmono = "ICI monotherapy")

#' @title Per-patient regimen table for the LMM frame
#' @description One row per patient of `df_long`. The regimen is read from the sheet
#'   of the patient's FIRST sample in the frame (T0, or T1 for T1-only patients such
#'   as NSCLC_64, whose T0 and T1 rows are two persons).
#' @param sheets named list (T0, T1) of data frames with Patient_ID, code, regime
#'   (1 = ICI monotherapy, 2 = chemo-immunotherapy, NA = unknown), radio (1 = annotated
#'   radiotherapy between T0 and T1)
#' @param df_long long frame with Patient_ID and Timepoint
therapy_regimen_table <- function(sheets, df_long) {
  ids   <- unique(as.character(df_long$Patient_ID))
  first <- tapply(as.character(df_long$Timepoint), as.character(df_long$Patient_ID), min)[ids]
  get   <- function(v) mapply(function(id, tp) { s <- sheets[[tp]]; s[[v]][match(id, s$Patient_ID)] },
                              ids, first, USE.NAMES = FALSE)
  data.frame(Patient_ID = ids, first_sheet = unname(first), code = get("code"),
             Rx = factor(c("IOmono", "ChemoIO")[get("regime")], levels = RX_LEVELS),   # regime code 1 = ICI mono, 2 = chemo-IO
             radio = get("radio"), stringsAsFactors = FALSE)
}

# LMM on Value / sd, back-scaled (same scaling as fit_feature_lmm)
.ts_fit <- function(d, m, rhs, reml = TRUE) {
  s <- sd(d[[m]], na.rm = TRUE); d$V <- d[[m]] / s
  f <- suppressMessages(suppressWarnings(lmerTest::lmer(
         as.formula(paste("V ~", rhs, "+ (1 | Patient_ID)")), data = d, REML = reml,
         control = lme4::lmerControl(calc.derivs = FALSE))))
  ct <- summary(f)$coefficients
  ct[, c("Estimate", "Std. Error")] <- ct[, c("Estimate", "Std. Error")] * s
  list(fit = f, ct = ct, n_pts = length(unique(d$Patient_ID)), n_obs = nrow(d))
}
.ts_term <- function(F, t) { x <- F$ct[t, ]; c(beta = x[["Estimate"]], se = x[["Std. Error"]], p = x[["Pr(>|t|)"]]) }
.ts_lrt  <- function(d, m) anova(.ts_fit(d, m, "Timepoint * Group * Rx - Timepoint:Group:Rx", reml = FALSE)$fit,
                                 .ts_fit(d, m, "Timepoint * Group * Rx", reml = FALSE)$fit)
.TG <- "TimepointT1:GroupRP"

# Group permuted within Rx strata; statistic = stratum-size-weighted mean PR - SD/PD delta.
# The seed is set on every call (as diag_56 did), so each call is reproducible on its own.
.ts_strat_perm <- function(w, B, seed) {
  st <- function(g) sum(sapply(levels(w$Rx), function(r) {
    i <- w$Rx == r; (mean(w$delta[i & g == "RP"]) - mean(w$delta[i & g == "SD_PD"])) * sum(i) })) / nrow(w)
  RNGkind("Mersenne-Twister"); set.seed(seed)
  obs <- st(w$Group)
  idx <- split(seq_len(nrow(w)), w$Rx)
  null <- replicate(B, { g <- w$Group; for (i in idx) g[i] <- g[i][sample.int(length(i))]; st(g) })
  c(obs = obs, p = (1 + sum(abs(null) >= abs(obs))) / (B + 1))
}

# paired delta frame (one row per patient with both samples)
.ts_delta <- function(dm, m) {
  w0 <- dm[dm$Timepoint == "T0", c("Patient_ID", "Group", "Rx", "radio", m)]
  w1 <- dm[dm$Timepoint == "T1", c("Patient_ID", m)]
  w  <- merge(w0, w1, by = "Patient_ID", suffixes = c("_T0", "_T1"))
  w$T0 <- w[[paste0(m, "_T0")]]; w$delta <- w[[paste0(m, "_T1")]] - w$T0
  w
}

#' @title Baseline-adjusted mean change by response x regimen
#' @description lm(delta ~ Group * Rx + T0) predicted at the mean T0 of `w`, with 95%
#'   CI. Per-patient adjusted changes delta - b_T0 * (T0 - mean T0) average exactly to
#'   the cell estimates (cell dummies + common T0 slope), so they can be drawn behind.
therapy_adjusted_change <- function(w) {
  f  <- lm(delta ~ Group * Rx + T0, w)
  nd <- expand.grid(Group = levels(w$Group), Rx = levels(w$Rx))
  nd$T0 <- mean(w$T0)
  pr <- predict(f, nd, interval = "confidence")
  cell <- interaction(w$Group, w$Rx)
  key  <- interaction(nd$Group, nd$Rx)
  out <- data.frame(nd[, c("Group", "Rx")], fit = pr[, "fit"], lwr = pr[, "lwr"], upr = pr[, "upr"],
                    n = as.vector(table(cell)[as.character(key)]),
                    raw_mean = as.vector(tapply(w$delta, cell, mean)[as.character(key)]),
                    T0_mean  = as.vector(tapply(w$T0, cell, mean)[as.character(key)]))
  pts <- data.frame(Patient_ID = w$Patient_ID, Group = w$Group, Rx = w$Rx,
                    adj_delta = w$delta - coef(f)[["T0"]] * (w$T0 - mean(w$T0)))
  list(cells = out, points = pts, slope_T0 = coef(f)[["T0"]], T0_ref = mean(w$T0))
}

#' @title Therapy sensitivity of the LMM gate
#' @param df_long long frame: Patient_ID, Timepoint (T0/T1), Group (SD_PD/RP), markers
#' @param gate character vector of fixed gate markers
#' @param regimen output of therapy_regimen_table()
#' @param B permutations; @param seed permutation seed
#' @param delta_ids optional Patient_IDs of the Δ frame (complete-case option): the paired-Δ
#'   analyses ([2] Wilcoxon / ANCOVA, [3], [4] permutation, adjusted change) use only these;
#'   the LMM rows stay on the full frame. NULL = all paired patients.
#' @return list(regimen, counts, markers = per-marker list, tests = tidy table)
run_therapy_sensitivity <- function(df_long, gate, regimen, B = 10000, seed = 56, delta_ids = NULL) {
  df <- df_long
  df$Patient_ID <- as.character(df$Patient_ID)
  df$Timepoint  <- factor(df$Timepoint, levels = c("T0", "T1"))
  df$Group      <- factor(df$Group, levels = c("SD_PD", "RP"))
  df$Rx    <- regimen$Rx[match(df$Patient_ID, regimen$Patient_ID)]
  df$radio <- regimen$radio[match(df$Patient_ID, regimen$Patient_ID)]
  dk <- df[!is.na(df$Rx), ]
  rad <- unique(dk$Patient_ID[dk$radio %in% 1])

  tests <- list(); res <- list()
  add <- function(m, test, n_pts, n_obs, n_paired, est = NA, se = NA, p = NA)
    tests[[length(tests) + 1]] <<- data.frame(marker = m, test = test, n_pts = n_pts, n_obs = n_obs,
                                              n_paired = n_paired, estimate = unname(est), se = unname(se),
                                              p = unname(p))
  for (m in gate) {
    dm <- dk[!is.na(dk[[m]]), ]
    dall <- df[!is.na(df[[m]]), ]
    # [1]
    F0 <- .ts_fit(dm, m, "Timepoint * Group")
    F1 <- .ts_fit(dm, m, "Timepoint * Group + Timepoint * Rx")
    F2 <- .ts_fit(dm, m, "Timepoint * Group * Rx")
    lrt <- .ts_lrt(dm, m)
    a0 <- .ts_term(F0, .TG); a1 <- .ts_term(F1, .TG); tr <- .ts_term(F1, "TimepointT1:RxIOmono")
    h  <- .ts_term(F2, "TimepointT1:GroupRP:RxIOmono")
    w_all <- .ts_delta(dm, m)
    w  <- if (is.null(delta_ids)) w_all else w_all[w_all$Patient_ID %in% delta_ids, , drop = FALSE]
    np <- nrow(w); np_lmm <- nrow(w_all)
    Fa <- .ts_fit(dall, m, "Timepoint * Group"); aa <- .ts_term(Fa, .TG)
    add(m, "reference LMM, all patients", Fa$n_pts, Fa$n_obs, sum(table(dall$Patient_ID) == 2), aa[1], aa[2], aa[3])
    add(m, "reference LMM, regimen known", F0$n_pts, F0$n_obs, np_lmm, a0[1], a0[2], a0[3])
    add(m, "[1] + Time x regimen", F1$n_pts, F1$n_obs, np_lmm, a1[1], a1[2], a1[3])
    add(m, "[1] Time x regimen (ICI mono vs chemo-IO change)", F1$n_pts, F1$n_obs, np_lmm, tr[1], tr[2], tr[3])
    add(m, "heterogeneity Time x response x regimen, Wald", F2$n_pts, F2$n_obs, np_lmm, h[1], h[2], h[3])
    add(m, "heterogeneity Time x response x regimen, LRT 1 df", F2$n_pts, F2$n_obs, np_lmm,
        lrt$Chisq[2], NA, lrt$`Pr(>Chisq)`[2])
    # [2] strata
    strat <- list()
    for (r in levels(dm$Rx)) {
      Fs <- .ts_fit(dm[dm$Rx == r, ], m, "Timepoint * Group"); a <- .ts_term(Fs, .TG); tt <- .ts_term(Fs, "TimepointT1")
      s  <- w[w$Rx == r, ]
      wx <- suppressWarnings(wilcox.test(delta ~ Group, s, exact = TRUE))$p.value
      an <- summary(lm(delta ~ Group + T0, s))$coefficients["GroupRP", ]
      strat[[r]] <- data.frame(marker = m, stratum = r, n_pts = Fs$n_pts, n_obs = Fs$n_obs, n_paired = nrow(s),
                               n_paired_PR = sum(s$Group == "RP"), beta = a[1], se = a[2], p = a[3],
                               sdpd_change = tt[1], sdpd_change_se = tt[2], p_wilcox = wx,
                               ancova_beta = an[1], p_ancova = an[4], row.names = NULL)
      s_all <- w_all[w_all$Rx == r, ]   # LMM rows report the full frame's pairs
      add(m, sprintf("[2] stratum %s: Time x response", RX_LABELS[[r]]), Fs$n_pts, Fs$n_obs, nrow(s_all), a[1], a[2], a[3])
      add(m, sprintf("[2] stratum %s: SD/PD change", RX_LABELS[[r]]), Fs$n_pts, Fs$n_obs,
          sum(s_all$Group == "SD_PD"), tt[1], tt[2], tt[3])
    }
    sall <- data.frame(marker = m, stratum = "all", n_pts = Fa$n_pts, n_obs = Fa$n_obs,
                       n_paired = sum(table(dall$Patient_ID) == 2), n_paired_PR = NA, beta = aa[1], se = aa[2],
                       p = aa[3], sdpd_change = .ts_term(Fa, "TimepointT1")[1], sdpd_change_se = NA,
                       p_wilcox = NA, ancova_beta = NA, p_ancova = NA, row.names = NULL)
    strata <- rbind(sall, do.call(rbind, strat))
    # [3]
    pr  <- .ts_strat_perm(w, B, seed)
    ols <- summary(lm(delta ~ Group + Rx, w))$coefficients
    anc <- summary(lm(delta ~ Group + T0 + Rx, w))$coefficients
    hT0 <- anova(lm(delta ~ Group + Rx + T0, w), lm(delta ~ Group * Rx + T0, w))
    add(m, sprintf("[3] response permuted within regimen (B = %d)", B), NA, NA, np, pr[["obs"]], NA, pr[["p"]])
    add(m, "[3] OLS delta ~ response + regimen", NA, NA, np, ols["GroupRP", 1], ols["GroupRP", 2], ols["GroupRP", 4])
    add(m, "[3] ANCOVA delta ~ response + T0 + regimen", NA, NA, np, anc["GroupRP", 1], anc["GroupRP", 2], anc["GroupRP", 4])
    add(m, "heterogeneity response x regimen on delta, T0-adjusted (F)", NA, NA, np, hT0$F[2], NA, hT0$`Pr(>F)`[2])
    # [4]
    no_rad <- NULL
    if (length(rad) > 0) {
      dr <- dm[!dm$Patient_ID %in% rad, ]; wr <- w[!w$Patient_ID %in% rad, ]
      nr_lmm <- sum(!w_all$Patient_ID %in% rad)
      F4 <- .ts_fit(dr, m, "Timepoint * Group + Timepoint * Rx"); a4 <- .ts_term(F4, .TG)
      l4 <- .ts_lrt(dr, m); p4 <- .ts_strat_perm(wr, B, seed)
      add(m, "[4] no radiotherapy T0-T1: + Time x regimen", F4$n_pts, F4$n_obs, nr_lmm, a4[1], a4[2], a4[3])
      add(m, "[4] no radiotherapy T0-T1: permutation", NA, NA, nrow(wr), p4[["obs"]], NA, p4[["p"]])
      add(m, "[4] no radiotherapy T0-T1: heterogeneity LRT", F4$n_pts, F4$n_obs, nr_lmm,
          l4$Chisq[2], NA, l4$`Pr(>Chisq)`[2])
      no_rad <- list(excluded = rad, adj = a4, lrt = l4, perm = p4)
    }
    res[[m]] <- list(adj = rbind(unadjusted = a0, time_x_rx_adjusted = a1), time_x_rx = tr,
                     heterogeneity_wald = h, heterogeneity_lrt = lrt, strata = strata, delta = w,
                     perm = c(pr, B = B), ols = ols, ancova = anc, heterogeneity_delta_T0 = hT0,
                     no_radio = no_rad, adjusted_change = therapy_adjusted_change(w))
  }
  pt <- df[!duplicated(df$Patient_ID), ]
  paired <- names(which(table(df$Patient_ID) == 2))
  counts <- list(n_pts = nrow(pt), n_obs = nrow(df), n_paired = length(paired),
                 rx_by_group = table(Rx = pt$Rx, Group = pt$Group, useNA = "ifany"),
                 rx_by_group_paired = table(Rx = pt$Rx[pt$Patient_ID %in% paired],
                                            Group = pt$Group[pt$Patient_ID %in% paired], useNA = "ifany"),
                 no_regimen = pt$Patient_ID[is.na(pt$Rx)], radio = rad)
  out <- list(regimen = regimen, counts = counts, markers = res, tests = do.call(rbind, tests),
              B = B, seed = seed)
  if (!is.null(delta_ids)) out$delta_ids <- delta_ids
  out
}

#' @title Figure: baseline-adjusted T0->T1 change by response within regimen
#' @description Points = adjusted cell means with 95% CI (therapy_adjusted_change);
#'   faint dots = per-patient baseline-adjusted changes. Facet per gate marker. No n or
#'   p is printed on the panels (user, 2026-10-27): both go in the caption.
#'   Requires modules_pub_style.R (theme_publication, pub_palette, pub_relabel_*).
#' @param res output of run_therapy_sensitivity()
#' @param adjusted FALSE draws raw (unadjusted) mean changes instead (supplementary)
pub_fig_therapy_change <- function(res, adjusted = TRUE) {
  ms <- names(res$markers)
  lab_m <- setNames(pub_relabel_marker(ms), ms)
  cells <- do.call(rbind, lapply(ms, function(m) {
    a <- res$markers[[m]]$adjusted_change$cells
    if (!adjusted) {
      w <- res$markers[[m]]$delta
      ci <- t(sapply(seq_len(nrow(a)), function(i) {
        x <- w$delta[w$Group == a$Group[i] & w$Rx == a$Rx[i]]
        mean(x) + c(0, -1, 1) * qt(0.975, length(x) - 1) * sd(x) / sqrt(length(x)) }))
      a$fit <- ci[, 1]; a$lwr <- ci[, 2]; a$upr <- ci[, 3]
    }
    cbind(marker = m, a)
  }))
  pts <- do.call(rbind, lapply(ms, function(m) {
    k <- res$markers[[m]]
    p <- k$adjusted_change$points
    if (!adjusted) p$adj_delta <- k$delta$delta[match(p$Patient_ID, k$delta$Patient_ID)]
    cbind(marker = m, p)
  }))
  for (d in c("cells", "pts")) {
    x <- get(d); x$marker <- factor(lab_m[x$marker], levels = lab_m)
    if (!is.null(x$Rx))    x$Rx <- factor(RX_LABELS[as.character(x$Rx)], levels = RX_LABELS)
    if (!is.null(x$Group)) x$Group <- factor(pub_relabel_group(x$Group), levels = c("PR", "SD/PD"))
    assign(d, x)
  }
  cols <- c("PR" = pub_palette[["combined"]], "SD/PD" = pub_palette[["unpen"]])
  dg <- position_dodge(width = 0.7)
  ggplot(cells, aes(Rx, fit, colour = Group)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
    geom_point(data = pts, aes(Rx, adj_delta, colour = Group), inherit.aes = FALSE, alpha = 0.3, size = 0.9,
               position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.7, seed = 1)) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.15, linewidth = 0.6, position = dg) +
    geom_point(size = 2.4, position = dg) +
    facet_wrap(~ marker, nrow = 1) +
    scale_colour_manual(values = cols) +
    scale_x_discrete(labels = c("Chemo-immunotherapy" = "Chemo-\nimmunotherapy", "ICI monotherapy" = "ICI\nmonotherapy")) +
    scale_y_continuous(expand = expansion(mult = 0.05)) +
    labs(x = NULL, y = if (adjusted) "T0→T1 change, baseline-adjusted\n(logit)"
                       else "T0→T1 change, unadjusted\n(logit)") +
    theme_publication(9) + theme(panel.grid.major.x = element_blank())
}
