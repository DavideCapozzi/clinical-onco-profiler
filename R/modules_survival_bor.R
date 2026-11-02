# R/modules_survival_bor.R
# ==============================================================================
# SURVIVAL BY BEST OVERALL RESPONSE AND BY Δ IMMUNE SCORE CUT-OFF (post hoc, descriptive)
# Frame: the step-06 Δ added-value frame (ClinImmune_Patients), survival from the
# experiment's T0 input, PFS event = the step-06 derived proxy.
#   1. follow-up   median follow-up, reverse Kaplan-Meier on OS (death = censoring)
#   2. BOR         OS / PFS by PR / SD / PD from the start of therapy: KM medians
#                  [95% CI], log-rank across the three groups and PR vs SD/PD, Cox HR
#                  PR vs SD/PD. Repeated with a landmark (default 3 months) as the
#                  guarantee-time sensitivity: BOR is only known after restaging.
#   3. score       OS / PFS by the Δ composite dichotomised at fixed (0 = decreased vs
#                  stable/increased) and quantile cut-offs, from the T1 landmark of
#                  post_hoc.km (as Figure 4); continuous Cox; the same within regimen.
#
# Standalone: tools/survival_bor.R --run <run> (no re-run; writes into the run's
# 07_post_hoc/bor_survival/). Not wired into src/07_post_hoc.R yet; to integrate, add
# inside its post hoc block (after 7d):
#   bs_cfg <- ph_cfg$bor_survival
#   if (!is.null(bs_cfg) && !isFALSE(bs_cfg$enabled))
#     bor_survival_run(config, file.path(out_root, "bor_survival"), by = "src/07_post_hoc.R (bor_survival)")
# with source("R/modules_survival_bor.R") at the top and, in the config,
#   post_hoc: {bor_survival: {enabled: true, landmark_months: 3}}
#
# Asserts (FATAL otherwise): frame n / groups = step-06 JSON; OS / PFS n and events =
# JSON survival; Cox OS PR vs SD/PD = JSON mediation response_to_survival_hr; Cox
# clinical + composite = JSON immune_hr; median split = post_hoc km/km_tests.csv (when
# present). No random numbers (no seed).
# Dependencies: survival, readxl, jsonlite, ggplot2, patchwork (+ modules_pub_style.R,
# utils_io.R for step_output_path / step_dir_name).
# ==============================================================================

BOR_LEVELS <- c("PR", "SD", "PD")
BOR_COLS   <- c(PR = "#0072B2", SD = "#E69F00", PD = "#D55E00")   # Okabe-Ito
BOR_REGIMENS <- c("1" = "ICI monotherapy", "2" = "chemo-immunotherapy")

.bor_lr_p <- function(sd) pchisq(sd$chisq, df = length(sd$n) - 1L, lower.tail = FALSE)

.bor_cox_row <- function(fit) {
  s <- summary(fit)
  zph <- tryCatch(survival::cox.zph(fit)$table[1, "p"], error = function(e) NA_real_)
  data.frame(hr = s$conf.int[1, 1], hr_lo = s$conf.int[1, 3], hr_hi = s$conf.int[1, 4],
             wald_p = s$coefficients[1, "Pr(>|z|)"], cox_zph_p = zph)
}

.bor_medians <- function(fit, prefix) {
  tab <- summary(fit)$table
  if (is.null(dim(tab))) tab <- t(tab)
  data.frame(group = sub(paste0("^", prefix, "="), "", rownames(tab)), n = tab[, "records"],
             events = tab[, "events"], median_months = tab[, "median"], lcl = tab[, "0.95LCL"],
             ucl = tab[, "0.95UCL"], row.names = NULL)
}

#' @title Median follow-up by reverse Kaplan-Meier (events and censoring swapped)
#' @param time OS in months; @param event 1 = death
bor_followup <- function(time, event) {
  ok <- is.finite(time) & !is.na(event)
  fit <- survival::survfit(survival::Surv(time[ok], 1L - event[ok]) ~ 1)
  tab <- summary(fit)$table
  cens <- time[ok][event[ok] == 0]
  data.frame(method = "reverse Kaplan-Meier on OS", n = sum(ok), n_alive = length(cens),
             median_months = tab[["median"]], lcl = tab[["0.95LCL"]], ucl = tab[["0.95UCL"]],
             min_followup_alive = if (length(cens)) min(cens) else NA_real_,
             max_followup_alive = if (length(cens)) max(cens) else NA_real_)
}

#' @title Survival by best overall response, optionally from a landmark
#' @param time months from the start of therapy; @param event 1 = event
#' @param group factor PR / SD / PD; @param landmark_months 0 = none (patients with
#'   time <= landmark are excluded and time restarts at the landmark)
#' @return list(patients, medians, tests, risk, km, excluded, landmark_months)
bor_group_survival <- function(time, event, group, landmark_months = 0, risk_times = NULL) {
  ok <- is.finite(time) & !is.na(event) & !is.na(group)
  k  <- ok & time > landmark_months
  d  <- data.frame(group = droplevels(factor(group[k], levels = BOR_LEVELS)),
                   time = time[k] - landmark_months, event = as.integer(event[k]))
  d$pr <- d$group == "PR"
  S  <- survival::Surv(d$time, d$event)
  sd3 <- survival::survdiff(S ~ d$group)
  sd2 <- survival::survdiff(S ~ d$pr)
  cx  <- survival::coxph(S ~ pr, data = d)
  base <- data.frame(n = nrow(d), events = sum(d$event))
  tests <- rbind(
    cbind(test = sprintf("log-rank %s (df %d)", paste(levels(d$group), collapse = " / "), nlevels(d$group) - 1L),
          base, chisq = sd3$chisq, p = .bor_lr_p(sd3),
          data.frame(hr = NA_real_, hr_lo = NA_real_, hr_hi = NA_real_, wald_p = NA_real_, cox_zph_p = NA_real_)),
    cbind(test = "log-rank + Cox, PR vs SD/PD (HR PR / SD+PD)", base, chisq = sd2$chisq,
          p = .bor_lr_p(sd2), .bor_cox_row(cx)))
  fit <- survival::survfit(S ~ group, data = d)
  sf  <- summary(fit, censored = TRUE)
  km  <- data.frame(group = sub("^group=", "", sf$strata), time = sf$time, n_risk = sf$n.risk,
                    n_event = sf$n.event, n_censor = sf$n.censor, surv = sf$surv)
  if (is.null(risk_times)) risk_times <- seq(0, 6 * ceiling(max(d$time) / 6), by = 6)
  rs <- summary(fit, times = risk_times, extend = TRUE)
  excl <- data.frame(group = group[ok & !k], time = time[ok & !k], event = event[ok & !k])
  list(patients = d, medians = .bor_medians(fit, "group"), tests = tests,
       risk = data.frame(group = sub("^group=", "", rs$strata), time = rs$time, n_risk = rs$n.risk),
       km = km, excluded = excl, landmark_months = landmark_months, cox = cx)
}

#' @title Survival by a dichotomised score: "low" = x <= cut (for cut 0: decreased)
#' @return list(test = one-row data.frame, medians); HR = low vs high
bor_score_cut <- function(time, event, x, cut, label) {
  ok  <- is.finite(time) & !is.na(event) & is.finite(x)
  low <- x[ok] <= cut; S <- survival::Surv(time[ok], event[ok])
  arm <- factor(ifelse(low, "low", "high"), levels = c("low", "high"))
  sd  <- survival::survdiff(S ~ arm)
  cx  <- survival::coxph(S ~ low)
  md  <- .bor_medians(survival::survfit(S ~ arm), "arm")
  list(test = cbind(cutoff = label, cut_value = cut, n_low = sum(low), events_low = sum(event[ok][low]),
                    n_high = sum(!low), events_high = sum(event[ok][!low]),
                    median_low = md$median_months[md$group == "low"], median_high = md$median_months[md$group == "high"],
                    logrank_p = .bor_lr_p(sd), .bor_cox_row(cx)),
       medians = cbind(cutoff = label, md))
}

#' @title Cox HR per unit of a continuous score (unadjusted)
bor_score_continuous <- function(time, event, x) {
  ok <- is.finite(time) & !is.na(event) & is.finite(x)
  .bor_cox_row(survival::coxph(survival::Surv(time[ok], event[ok]) ~ x[ok]))
}

#' @title Score vs survival within strata (treatment regimen)
#' @return data.frame: per stratum (continuous and low-vs-high at cut), stratified Cox,
#'   interaction LRT (score x stratum)
bor_score_by_stratum <- function(time, event, x, cut, stratum) {
  ok <- is.finite(time) & !is.na(event) & is.finite(x) & !is.na(stratum)
  d  <- data.frame(time = time[ok], event = event[ok], x = x[ok], low = x[ok] <= cut, s = factor(stratum[ok]))
  row <- function(lbl, term, fit, dd) cbind(stratum = lbl, model = term, n = nrow(dd), events = sum(dd$event),
                                            .bor_cox_row(fit))
  out <- list()
  for (lv in levels(d$s)) {
    dd <- d[d$s == lv, ]
    out[[length(out) + 1]] <- row(lv, "continuous (per unit)", survival::coxph(survival::Surv(time, event) ~ x, data = dd), dd)
    if (length(unique(dd$low)) == 2)
      out[[length(out) + 1]] <- row(lv, sprintf("x <= %g vs > %g", cut, cut),
                                    survival::coxph(survival::Surv(time, event) ~ low, data = dd), dd)
  }
  out[[length(out) + 1]] <- row("stratified by regimen", "continuous (per unit)",
                                survival::coxph(survival::Surv(time, event) ~ x + survival::strata(s), data = d), d)
  out[[length(out) + 1]] <- row("stratified by regimen", sprintf("x <= %g vs > %g", cut, cut),
                                survival::coxph(survival::Surv(time, event) ~ low + survival::strata(s), data = d), d)
  res <- do.call(rbind, out)
  f0 <- survival::coxph(survival::Surv(time, event) ~ x + s, data = d)
  f1 <- survival::coxph(survival::Surv(time, event) ~ x * s, data = d)
  attr(res, "interaction_p") <- anova(f0, f1)[2, "Pr(>|Chi|)"]   # match on name (coxph uses |Chi|)
  res
}

#' @title Figure: KM by best overall response, one panel per endpoint (A, B)
#' @param res_list named list of bor_group_survival() results; @param y_labs y-axis titles
pub_fig_bor_km <- function(res_list, y_labs) {
  panel <- function(res, y_lab) {
    km <- res$km; lv <- levels(res$patients$group)
    curve <- rbind(data.frame(group = lv, time = 0, surv = 1), km[, c("group", "time", "surv")])
    curve$group <- factor(curve$group, levels = lv)
    cens <- km[km$n_censor > 0, ]; cens$group <- factor(cens$group, levels = lv)
    p3 <- res$tests$p[1]
    lab <- if (p3 < 0.001) "Log-rank p < 0.001" else sprintf("Log-rank p = %.3f", p3)
    x_max <- max(res$risk$time)
    ggplot(curve, aes(time, surv, colour = group)) +
      geom_step(linewidth = 0.7) +
      geom_point(data = cens, aes(time, surv, colour = group), shape = 3, size = 1.6, stroke = 0.6,
                 inherit.aes = FALSE, show.legend = FALSE) +   # PD has no censoring in PFS: legends would not merge
      annotate("text", x = x_max, y = 1, label = lab, hjust = 1, vjust = 1, size = 3.2) +
      scale_colour_manual(values = BOR_COLS[lv]) +
      scale_x_continuous(breaks = seq(0, x_max, by = 6), limits = c(0, x_max),
                         expand = expansion(mult = c(0.02, 0.02))) +
      scale_y_continuous(limits = c(0, 1), labels = function(b) paste0(b * 100, "%"),
                         expand = expansion(mult = c(0.01, 0.02))) +
      labs(x = if (res$landmark_months > 0) sprintf("Months from landmark (%g months)", res$landmark_months)
               else "Months from start of therapy", y = y_lab) +
      theme_publication(9) + theme(panel.grid = element_blank())
  }
  ps <- Map(panel, res_list, y_labs)
  pub_tag(patchwork::wrap_plots(ps, nrow = 1) + patchwork::plot_layout(guides = "collect")) &
    theme(legend.position = "bottom")
}

#' @title Driver: run the whole analysis on a finished run and write it to out_dir
#' @param config post hoc pass config (prepare_pass_config(..., "post_hoc", run_root))
#' @param out_dir output directory (<run>/<exp>/07_post_hoc/bor_survival)
#' @param by caller recorded in the manifest
#' @return invisible list of results
bor_survival_run <- function(config, out_dir, by = "tools/survival_bor.R") {
  fatal <- function(...) stop(sprintf("[FATAL] bor_survival: %s", sprintf(...)), call. = FALSE)
  exp  <- config$project_name
  ph   <- if (!is.null(config$post_hoc)) config$post_hoc else list()
  bcfg <- if (!is.null(ph$bor_survival)) ph$bor_survival else list()
  lm_bor   <- if (!is.null(bcfg$landmark_months)) as.numeric(bcfg$landmark_months) else 3
  lm_days  <- if (!is.null(ph$km$landmark_days)) as.numeric(ph$km$landmark_days) else 21
  lm_score <- lm_days / 30.4375                                  # as post_hoc km
  q_cuts   <- if (!is.null(bcfg$score_quantiles)) as.numeric(unlist(bcfg$score_quantiles)) else c(1/3, 1/2, 2/3)
  if (!is.finite(lm_bor) || lm_bor < 0) fatal("post_hoc.bor_survival.landmark_months must be >= 0")

  ml_json <- step_output_path(config, 6, sprintf("Machine_Metrics_ML_%s", exp), "json")
  ml_xlsx <- step_output_path(config, 6, sprintf("ML_Classification_Report_%s", exp), "xlsx")
  if (!all(file.exists(c(ml_json, ml_xlsx)))) fatal("step-06 artifacts missing in %s", dirname(ml_json))
  av <- jsonlite::fromJSON(ml_json, simplifyVector = FALSE)$clinical_immune_added_value
  if (is.null(av)) fatal("step 06 wrote no clinical_immune_added_value node")
  if (!identical(av$composite_timepoint, "delta")) fatal("added-value composite is not Δ")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  con <- file(file.path(out_dir, "console.txt"), open = "wt"); on.exit(close(con))
  say  <- function(...) { t <- sprintf(...); message(t); cat(t, "\n", sep = "", file = con) }
  show <- function(x) invisible(lapply(capture.output(print(x, row.names = FALSE)), function(l) {
    message("  ", l); cat("  ", l, "\n", sep = "", file = con) }))
  rnd <- function(df) { for (k in names(df)) if (is.numeric(df[[k]])) df[[k]] <- signif(df[[k]], 4); df }

  # Frame (= step 06)
  pts <- as.data.frame(readxl::read_excel(ml_xlsx, sheet = "ClinImmune_Patients"))
  ids <- as.character(pts$Patient_ID)
  if (length(ids) != av$n || anyDuplicated(ids)) fatal("ClinImmune_Patients has %d IDs, JSON n = %d", length(ids), av$n)
  in_t0 <- config$input_file_t0
  s0 <- as.data.frame(readxl::read_excel(in_t0, sheet = 1)); names(s0) <- trimws(names(s0))
  s0$Patient_ID <- as.character(s0$Patient_ID)
  i0 <- match(ids, s0$Patient_ID); if (anyNA(i0)) fatal("Δ-frame IDs missing from %s", in_t0)
  sv <- config$machine_learning$clinical_model$survival
  SURV <- c(os = if (!is.null(sv$os)) sv$os else "OS", pfs = if (!is.null(sv$pfs)) sv$pfs else "PFS",
            dead = if (!is.null(sv$death)) sv$death else "Alive_0/Dead_1")
  tcol <- config$clinical$target_column
  if (!all(c(SURV, tcol) %in% names(s0))) fatal("columns missing in %s: %s", in_t0,
                                                paste(setdiff(c(SURV, tcol), names(s0)), collapse = ", "))
  num <- function(k) suppressWarnings(as.numeric(s0[[k]][i0]))
  os <- num(SURV[["os"]]); pfs <- num(SURV[["pfs"]]); dth <- as.integer(num(SURV[["dead"]]))
  pfs_ev <- as.integer(dth == 1 | (dth == 0 & pfs < os))        # step-06 derived proxy
  sub <- ph$grouping$subcodes; if (is.null(sub)) sub <- list(PR = 2, SD = 3, PD = 4)
  code <- num(tcol)
  bor  <- factor(names(sub)[match(code, unlist(sub))], levels = BOR_LEVELS)
  if (anyNA(bor)) fatal("best-response codes outside %s for %s", paste(unlist(sub), collapse = "/"),
                        paste(ids[is.na(bor)], collapse = ", "))
  if (!identical(bor == "PR", pts$True_Group == av$positive_label)) fatal("PR / SD+PD differs from step-06 True_Group")
  say("bor_survival  %s | exp %s | by %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), exp, by)
  say("frame: step-06 Δ added-value frame, n = %d (%s) = JSON n %d (%d %s / %d %s): PASS", length(ids),
      paste(sprintf("%s %d", BOR_LEVELS, table(bor)), collapse = " / "), av$n, av$n_pos, av$positive_label,
      av$n_neg, av$negative_label)
  EP <- list(OS = list(time = os, event = dth, label = "Overall survival",
                       def = sprintf("%s (months from start of therapy), event = %s == 1", SURV[["os"]], SURV[["dead"]])),
             PFS = list(time = pfs, event = pfs_ev, label = "Progression-free survival",
                        def = sprintf("%s (months), event = death | (alive & PFS < OS): derived proxy, as step 06", SURV[["pfs"]])))
  for (e in names(EP)) {
    J <- av$survival[[e]]; ok <- is.finite(EP[[e]]$time) & EP[[e]]$time > 0 & !is.na(EP[[e]]$event)
    if (sum(ok) != J$n || sum(EP[[e]]$event[ok]) != J$events)
      fatal("%s n/events %d/%d != JSON %d/%d", e, sum(ok), sum(EP[[e]]$event[ok]), J$n, J$events)
    say("%s: %s; n = %d, %d events = JSON: PASS", e, EP[[e]]$def, J$n, J$events)
  }

  # 1. Follow-up
  fu <- bor_followup(os, dth)
  say("\n==== 1. median follow-up (reverse Kaplan-Meier on OS): %.1f months [95%% CI %.1f, %s]; %d alive, follow-up %g-%g months",
      fu$median_months, fu$lcl, if (is.na(fu$ucl)) "NR" else sprintf("%.1f", fu$ucl), fu$n_alive,
      fu$min_followup_alive, fu$max_followup_alive)

  # 2. BOR
  BOR <- list()
  for (e in names(EP)) for (lm in unique(c(0, lm_bor))) {
    key <- sprintf("%s_landmark%g", e, lm)
    BOR[[key]] <- r <- bor_group_survival(EP[[e]]$time, EP[[e]]$event, bor, landmark_months = lm)
    say("\n==== 2. %s by best overall response, %s", e,
        if (lm == 0) "from the start of therapy" else sprintf("landmark %g months (guarantee-time sensitivity)", lm))
    if (lm > 0) say("excluded (%s <= %g months): %d (%s)", e, lm, nrow(r$excluded),
                    paste(sprintf("%s %d", BOR_LEVELS, table(factor(r$excluded$group, levels = BOR_LEVELS))), collapse = " / "))
    say("medians (months) [95%% CI, log]:"); show(rnd(r$medians))
    say("tests:"); show(rnd(r$tests))
  }
  hr_os <- BOR$OS_landmark0$tests$hr[2]; J_hr <- av$survival$dissociation$mediation$response_to_survival_hr
  if (!is.null(J_hr)) {
    if (abs(hr_os - J_hr) > 5e-5 + 1e-9) fatal("Cox OS PR vs SD/PD HR %.4f != JSON response_to_survival_hr %.4f", hr_os, J_hr)
    say("\ncheck: Cox OS PR vs SD/PD HR %.4f = JSON dissociation.mediation.response_to_survival_hr %.4f: PASS", hr_os, J_hr)
  } else say("\ncheck: JSON response_to_survival_hr absent (mediation off): not checked")

  # 3. Score (Δ composite) from the T1 landmark, as Figure 4
  x <- pts$Composite
  cuts <- c(c("0 (decreased vs stable/increased)" = 0),
            setNames(stats::quantile(x, q_cuts, names = FALSE), sprintf("quantile %.3g", q_cuts)))
  SC <- list(); SCm <- list(); SR <- list(); inter <- c()
  # continuous / adjusted rows in the cut-off table layout (arm columns NA)
  sc_row <- function(e, label, cox) data.frame(endpoint = e, cutoff = label, cut_value = NA_real_,
    n_low = NA_integer_, events_low = NA_integer_, n_high = NA_integer_, events_high = NA_integer_,
    median_low = NA_real_, median_high = NA_real_, logrank_p = NA_real_,
    cox[, c("hr", "hr_lo", "hr_hi", "wald_p", "cox_zph_p")])
  for (e in names(EP)) {
    if (isTRUE(any(EP[[e]]$time <= lm_score))) fatal("%s <= T1 landmark for some patients: exclusions not handled", e)
    t_lm <- EP[[e]]$time - lm_score
    cr <- lapply(names(cuts), function(nm) bor_score_cut(t_lm, EP[[e]]$event, x, cuts[[nm]], nm))
    SC[[e]]  <- rbind(cbind(endpoint = e, do.call(rbind, lapply(cr, `[[`, "test"))),
                      sc_row(e, "continuous (per unit)", bor_score_continuous(t_lm, EP[[e]]$event, x)))
    SCm[[e]] <- cbind(endpoint = e, do.call(rbind, lapply(cr, `[[`, "medians")))
  }
  # check: median split = post_hoc km (HR there is high vs low -> reciprocal)
  km_csv <- file.path(step_dir(config, 7), "km", "km_tests.csv")
  if (file.exists(km_csv)) {
    kt <- utils::read.csv(km_csv, check.names = FALSE); kt <- kt[kt$analysis == "primary", ]
    for (e in names(EP)) {
      m <- SC[[e]][SC[[e]]$cutoff == "quantile 0.5", ]; k <- kt[kt$endpoint == e, ]
      if (abs(m$hr - 1 / k$hr) > 1e-6 || abs(m$logrank_p - k$p) > 1e-8)
        fatal("%s median split HR %.6f / p %.6g != km 1/HR %.6f / p %.6g", e, m$hr, m$logrank_p, 1 / k$hr, k$p)
    }
    say("\ncheck: median split = post_hoc km/km_tests.csv (HR reciprocal, log-rank p) for OS and PFS: PASS")
  } else say("\ncheck: km/km_tests.csv absent: median split not cross-checked")
  # check: step-06 Cox (clinical vars + composite, no landmark) reproduced
  cl <- unlist(av$clinical_vars)
  C  <- as.data.frame(lapply(cl, function(col) {
    v <- if (col %in% names(pts)) pts[[col]] else s0[[col]][i0]
    v <- suppressWarnings(as.numeric(v)); v[is.na(v)] <- median(v, na.rm = TRUE); v }))
  names(C) <- make.names(cl, unique = TRUE)
  for (e in names(EP)) {
    fk <- survival::coxph(survival::Surv(EP[[e]]$time, EP[[e]]$event) ~ ., data = cbind(C, comp = x))
    hr <- unname(exp(coef(fk))["comp"]); J <- av$survival[[e]]
    if (abs(hr - J$immune_hr) > 5e-4) fatal("%s Cox HR %.4f != JSON immune_hr %.4f", e, hr, J$immune_hr)
    s  <- summary(fk)
    SC[[e]] <- rbind(SC[[e]], sc_row(e, sprintf("continuous, adjusted for %s (step 06, no landmark)", paste(cl, collapse = " + ")),
                                     data.frame(hr = hr, hr_lo = s$conf.int["comp", 3], hr_hi = s$conf.int["comp", 4],
                                                wald_p = s$coefficients["comp", "Pr(>|z|)"], cox_zph_p = NA_real_)))
    say("check: %s Cox %s + composite HR %.4f = JSON immune_hr %.4f: PASS", e, paste(cl, collapse = " + "), hr, J$immune_hr)
  }
  # regimen (post_hoc.therapy files, parser v3.1)
  reg <- rep(NA_character_, length(ids)); th <- ph$therapy
  if (!is.null(th$file_t0) && file.exists(th$file_t0)) {
    rx <- as.data.frame(readxl::read_excel(th$file_t0)); rx$Patient_ID <- as.character(rx$Patient_ID)
    reg <- unname(BOR_REGIMENS[as.character(rx$Terapia_regime_IOmono_1_chemoIO_2[match(ids, rx$Patient_ID)])])
    for (e in names(EP)) {
      sr <- bor_score_by_stratum(EP[[e]]$time - lm_score, EP[[e]]$event, x, 0, reg)
      inter[[e]] <- attr(sr, "interaction_p"); SR[[e]] <- cbind(endpoint = e, sr)
    }
  }
  for (e in names(EP)) {
    say("\n==== 3. %s by Δ immune score (composite), from the T1 landmark (%g days); HR = score <= cut vs > cut", e, lm_days)
    show(rnd(SC[[e]][, -1]))
    say("medians (months from the landmark):"); show(rnd(SCm[[e]][, -1]))
    if (length(SR)) {
      say("by regimen (%s; no regimen: %s); HR per unit / <= 0 vs > 0:", basename(th$file_t0),
          paste(ids[is.na(reg)], collapse = ", "))
      show(rnd(SR[[e]][, -1])); say("interaction score x regimen (LRT, continuous): p = %.3g", inter[[e]])
    }
  }

  # Paste-ready values
  f1 <- function(v) ifelse(is.na(v), "NR", sprintf("%.1f", v))
  med_txt <- function(r) paste(sprintf("%s %s (95%% CI %s-%s)", r$medians$group, f1(r$medians$median_months),
                                       f1(r$medians$lcl), f1(r$medians$ucl)), collapse = "; ")
  p_txt <- function(p) if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)
  say("\n==== paragraph values (Δ frame, n = %d)", length(ids))
  say("follow-up %s months (95%% CI %s-%s); %d PFS events (proxy), %d deaths", f1(fu$median_months), f1(fu$lcl),
      f1(fu$ucl), sum(pfs_ev), sum(dth))
  for (e in names(EP)) {
    r <- BOR[[sprintf("%s_landmark0", e)]]; t2 <- r$tests[2, ]
    say("%s by BOR: %s; log-rank 3 groups %s; HR PR vs SD/PD %.2f (95%% CI %.2f-%.2f)", e, med_txt(r),
        p_txt(r$tests$p[1]), t2$hr, t2$hr_lo, t2$hr_hi)
    if (lm_bor > 0) { rl <- BOR[[sprintf("%s_landmark%g", e, lm_bor)]]; tl <- rl$tests[2, ]
      say("  landmark %g months: n %d (%d excluded), HR PR vs SD/PD %.2f (95%% CI %.2f-%.2f), %s", lm_bor,
          tl$n, nrow(rl$excluded), tl$hr, tl$hr_lo, tl$hr_hi, p_txt(tl$wald_p)) }
    s0r <- SC[[e]][SC[[e]]$cut_value %in% 0 & !is.na(SC[[e]]$n_low), ][1, ]
    say("%s score decreased (n %d) vs stable/increased (n %d): median %s vs %s months; HR %.2f (95%% CI %.2f-%.2f), p = %.2f",
        e, s0r$n_low, s0r$n_high, f1(s0r$median_low), f1(s0r$median_high), s0r$hr, s0r$hr_lo, s0r$hr_hi, s0r$wald_p)
  }

  # Write
  bind <- function(k) do.call(rbind, lapply(names(BOR), function(n) cbind(analysis = n, BOR[[n]][[k]])))
  utils::write.csv(fu, file.path(out_dir, "bor_followup.csv"), row.names = FALSE)
  utils::write.csv(bind("medians"), file.path(out_dir, "bor_medians.csv"), row.names = FALSE)
  utils::write.csv(bind("tests"), file.path(out_dir, "bor_tests.csv"), row.names = FALSE)
  utils::write.csv(bind("risk"), file.path(out_dir, "bor_risk_table.csv"), row.names = FALSE)
  utils::write.csv(do.call(rbind, SC), file.path(out_dir, "score_cutoff_tests.csv"), row.names = FALSE)
  utils::write.csv(do.call(rbind, SCm), file.path(out_dir, "score_cutoff_medians.csv"), row.names = FALSE)
  if (length(SR)) utils::write.csv(do.call(rbind, SR), file.path(out_dir, "score_by_regimen.csv"), row.names = FALSE)
  utils::write.csv(data.frame(Patient_ID = ids, BOR = bor, regimen = reg, composite = x, OS = os, death = dth,
                              PFS = pfs, pfs_event = pfs_ev),
                   file.path(out_dir, "bor_patients.csv"), row.names = FALSE)
  saveRDS(list(followup = fu, bor = BOR, score = SC, score_medians = SCm, regimen = SR), file.path(out_dir, "bor_survival.rds"))
  save_pub_figure(pub_fig_bor_km(BOR[c("OS_landmark0", "PFS_landmark0")], c(EP$OS$label, EP$PFS$label)),
                  file.path(out_dir, "FigureBOR_KM_OS_PFS"), PUB_W2, 95, config = config)
  sha <- tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE), error = function(e) NA_character_)
  dirty <- tryCatch(length(system2("git", c("status", "--porcelain", "R/modules_survival_bor.R", "tools/survival_bor.R"),
                                   stdout = TRUE, stderr = FALSE)) > 0, error = function(e) NA)
  jsonlite::write_json(list(
    step = by, module = "R/modules_survival_bor.R", date = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    experiment = exp, git_sha = if (length(sha)) sha[1] else NA_character_, module_uncommitted = dirty,
    module_md5 = if (file.exists("R/modules_survival_bor.R")) unname(tools::md5sum("R/modules_survival_bor.R")) else NA,
    input_t0 = in_t0, input_t0_md5 = unname(tools::md5sum(in_t0)),
    therapy_file = if (!is.null(th$file_t0)) th$file_t0 else NA,
    frame = list(n = length(ids), PR = sum(bor == "PR"), SD = sum(bor == "SD"), PD = sum(bor == "PD")),
    endpoints = lapply(EP, `[[`, "def"), events = list(OS = sum(dth), PFS = sum(pfs_ev)),
    bor_landmark_months = lm_bor, score_landmark_days = lm_days, score_cuts = as.list(cuts),
    km_ci = "log (survfit default, as post_hoc km)",
    checks = list(frame_and_events_vs_step06 = TRUE, response_to_survival_hr = !is.null(J_hr),
                  median_split_vs_km = file.exists(km_csv), step06_cox_immune_hr = TRUE)),
    file.path(out_dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
  message(sprintf("   [PostHoc] bor_survival -> %s", out_dir))
  invisible(list(followup = fu, bor = BOR, score = SC, regimen = SR))
}
