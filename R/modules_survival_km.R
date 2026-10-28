# R/modules_survival_km.R
# ==============================================================================
# KAPLAN-MEIER MODULE — OS / PFS by median split of the Δ immune composite
# (post hoc, descriptive). The composite is the step-06 added-value composite
# (mean Δz of the fixed LMM gate markers); the split is at its median in the frame.
#
# Pure functions: vectors in, R objects out. No paths, no file I/O — the driver
# (src/07_post_hoc.R, km/) owns all paths. Post hoc and descriptive.
#
# Landmark: Δ is only known at T1, so time runs from the landmark (T1) and patients
# with OS <= landmark are excluded (recorded). With OS in whole months (min 1) and a
# 21-day landmark nobody is excluded and every test equals the un-landmarked one up to
# a constant time shift.
#
# Tests:
#   primary   log-rank Δ > median vs Δ <= median (survdiff), Cox HR [95% CI], cox.zph
#   secondary the same stratified by best response (post hoc; does the split separate
#             OS within PR and within SD/PD?) — reported in the tests table only
# Dependencies: survival, ggplot2, patchwork (+ modules_pub_style.R for the figure)
# ==============================================================================

KM_ARMS <- c(low = "Δ ≤ median", high = "Δ > median")
KM_COLS <- setNames(c("#2166AC", "#B2182B"), unname(KM_ARMS))   # blue / red (user, 2026-10-27)

.km_logrank_p <- function(sd) pchisq(sd$chisq, df = length(sd$n) - 1L, lower.tail = FALSE)

.km_cox_row <- function(fit, term = "highTRUE") {
  s <- summary(fit)
  c(hr = s$conf.int[term, 1], hr_lo = s$conf.int[term, 3], hr_hi = s$conf.int[term, 4],
    wald_p = s$coefficients[term, "Pr(>|z|)"])
}

#' @title Survival by median split of a continuous score, with a landmark
#' @param time OS or PFS in months from the start of therapy
#' @param event 1 = event (death; death or progression for PFS), 0 = censored
#' @param x score to split (here the Δ composite)
#' @param ids patient IDs (anonymised)
#' @param group best-response group (used only for the secondary stratified tests)
#' @param landmark_months landmark in months (T1); 0 = no landmark
#' @param risk_times times (months from the landmark) for the number-at-risk table
#'   (default every 6 months up to the first multiple of 6 past the last follow-up)
#' @return list(patients, tests, km, medians, risk, split_by_group, median_x, n_excluded, ...)
km_median_split <- function(time, event, x, ids, group, landmark_months = 0,
                            risk_times = NULL) {
  stopifnot(length(time) == length(x), length(event) == length(x), length(ids) == length(x),
            length(group) == length(x))
  ok_data <- is.finite(time) & !is.na(event) & is.finite(x)
  ok_lm   <- ok_data & time > landmark_months
  excluded <- data.frame(Patient_ID = ids[ok_data & !ok_lm], time_start = time[ok_data & !ok_lm],
                         event = event[ok_data & !ok_lm])
  d <- data.frame(Patient_ID = ids[ok_lm], group = group[ok_lm], x = x[ok_lm],
                  time_start = time[ok_lm], time = time[ok_lm] - landmark_months,
                  event = as.integer(event[ok_lm]), stringsAsFactors = FALSE)
  med <- median(d$x)
  d$high <- d$x > med
  d$arm  <- factor(ifelse(d$high, KM_ARMS[["high"]], KM_ARMS[["low"]]), levels = unname(KM_ARMS))

  S   <- survival::Surv(d$time, d$event)
  sd0 <- survival::survdiff(S ~ d$high)
  sds <- survival::survdiff(S ~ d$high + survival::strata(d$group))
  c0  <- survival::coxph(S ~ high, data = d)
  cs  <- survival::coxph(S ~ high + survival::strata(group), data = d)
  zph <- tryCatch(survival::cox.zph(c0)$table["high", "p"], error = function(e) NA_real_)
  arm_n <- function() c(n_low = sum(!d$high), events_low = sum(d$event[!d$high]),
                         n_high = sum(d$high), events_high = sum(d$event[d$high]))
  tests <- rbind(
    data.frame(analysis = "primary", test = "log-rank, Δ > median vs Δ ≤ median",
               t(arm_n()), chisq = sd0$chisq, p = .km_logrank_p(sd0),
               t(.km_cox_row(c0)), cox_zph_p = zph, check.names = FALSE),
    data.frame(analysis = "secondary (post hoc)", test = "log-rank stratified by best response",
               t(arm_n()), chisq = sds$chisq, p = .km_logrank_p(sds),
               t(.km_cox_row(cs)), cox_zph_p = NA_real_, check.names = FALSE))
  rownames(tests) <- NULL

  fit <- survival::survfit(S ~ arm, data = d)
  sf  <- summary(fit, censored = TRUE)
  arm_of <- function(strata) sub("^arm=", "", as.character(strata))
  km <- data.frame(arm = arm_of(sf$strata), time = sf$time, n_risk = sf$n.risk, n_event = sf$n.event,
                   n_censor = sf$n.censor, surv = sf$surv, lower = sf$lower, upper = sf$upper)
  tab <- summary(fit)$table
  medians <- data.frame(arm = arm_of(rownames(tab)), n = tab[, "records"], events = tab[, "events"],
                        median_months = tab[, "median"], lcl = tab[, "0.95LCL"], ucl = tab[, "0.95UCL"],
                        row.names = NULL)
  if (is.null(risk_times)) risk_times <- seq(0, 6 * ceiling(max(d$time) / 6), by = 6)
  rs <- summary(fit, times = risk_times, extend = TRUE)
  risk <- data.frame(arm = arm_of(rs$strata), time = rs$time, n_risk = rs$n.risk)

  list(patients = d, tests = tests, km = km, medians = medians, risk = risk,
       split_by_group = table(arm = d$arm, group = d$group), median_x = med,
       landmark_months = landmark_months, n_excluded = nrow(excluded), excluded = excluded,
       fit = fit, cox = c0, cox_strat = cs)
}

#' @title Figure panel: Kaplan-Meier curves by median split, log-rank p in the corner
#' @description Step curves from the landmark with censor ticks, red/blue, no grid, primary
#'   log-rank p top-right, colour legend at the bottom (user, 2026-10-27). Requires
#'   modules_pub_style.R (theme_publication).
#' @param res output of km_median_split()
#' @param y_lab y-axis title (endpoint name)
#' @param x_max last month on the axis (default: first multiple of 6 past the last follow-up)
pub_fig_km_median <- function(res, y_lab = "Overall survival", x_max = max(res$risk$time)) {
  km <- res$km
  # step curves need the (0, 1) start of each arm
  curve <- rbind(data.frame(arm = unname(KM_ARMS), time = 0, surv = 1), km[, c("arm", "time", "surv")])
  curve <- curve[order(curve$arm, curve$time), ]
  cens  <- km[km$n_censor > 0, ]
  for (v in c("curve", "cens")) { z <- get(v); z$arm <- factor(z$arm, levels = unname(KM_ARMS)); assign(v, z) }
  p_lr <- res$tests$p[res$tests$analysis == "primary"]
  lab  <- if (p_lr < 0.001) "Log-rank p < 0.001" else sprintf("Log-rank p = %.3f", p_lr)
  ggplot(curve, aes(time, surv, colour = arm)) +
    geom_step(linewidth = 0.7) +
    geom_point(data = cens, aes(time, surv, colour = arm), shape = 3, size = 1.6, stroke = 0.6,
               inherit.aes = FALSE) +
    annotate("text", x = x_max, y = 1, label = lab, hjust = 1, vjust = 1, size = 3.2) +
    scale_colour_manual(values = KM_COLS) +
    scale_x_continuous(breaks = seq(0, x_max, by = 6), limits = c(0, x_max),
                       expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(limits = c(0, 1), labels = function(b) paste0(b * 100, "%"),
                       expand = expansion(mult = c(0.01, 0.02))) +
    labs(x = "Months from T1 (landmark)", y = y_lab) +
    theme_publication(9) + theme(panel.grid = element_blank())
}

#' @title Figure: one Kaplan-Meier panel per endpoint (A, B, ...), shared bottom legend
#' @param res_list named list of km_median_split() results, in panel order
#' @param y_labs y-axis titles, same order
pub_fig_km_panels <- function(res_list, y_labs) {
  ps <- Map(pub_fig_km_median, res_list, y_labs)
  pub_tag(patchwork::wrap_plots(ps, nrow = 1) + patchwork::plot_layout(guides = "collect")) &
    theme(legend.position = "bottom")
}
