# src/07_post_hoc.R
# ==============================================================================
# STEP 07: POST HOC ANALYSES (on the finished step-04 / step-06 results)
# Description: Sensitivity and descriptive analyses that re-use the run's own
#              frames and its fixed gate; no marker re-selection, no new model.
#   grouping/        PR / SD / PD response grouping of the LMM frame
#   duplicate_rows/  drop the T1 rows listed in post_hoc.duplicate_rows.drop_t1
#   therapy/         ICI monotherapy vs chemo-immunotherapy (gate change)
#   km/              Kaplan-Meier OS / PFS by median split of the Δ composite
# Every analysis asserts that it reproduces its step-04 / step-06 inputs (FATAL
# otherwise). Self-skips, with a logged reason, when those artifacts are missing.
# Config block: post_hoc.{grouping, duplicate_rows, therapy, km}.
# ==============================================================================

source("R/utils_io.R")
source("R/utils_figures.R")
source("R/modules_longitudinal.R")
source("R/modules_ml.R")
source("R/modules_pub_style.R")
source("R/modules_posthoc_grouping.R")
source("R/modules_posthoc_duplicates.R")
source("R/modules_therapy_sensitivity.R")
source("R/modules_survival_km.R")

message("\n=== PIPELINE STEP 7: POST HOC ANALYSES ===")

if (!exists("config")) stop("[FATAL] Global configuration object not detected in environment.")
ph_cfg <- if (!is.null(config$post_hoc)) config$post_hoc else list()
exp    <- config$project_name

rds_std  <- step_output_path(config, 1, sprintf("data_processed_%s_standard", exp), "rds")
rds_long <- step_output_path(config, 1, sprintf("data_processed_%s_longitudinal", exp), "rds")
lmm_json <- step_output_path(config, 4, sprintf("Machine_Metrics_LMM_%s", exp), "json")
lmm_xlsx <- step_output_path(config, 4, sprintf("Longitudinal_LMM_Report_%s", exp), "xlsx")
ml_json  <- step_output_path(config, 6, sprintf("Machine_Metrics_ML_%s", exp), "json")
ml_xlsx  <- step_output_path(config, 6, sprintf("ML_Classification_Report_%s", exp), "xlsx")
needed   <- c(rds_std, rds_long, lmm_json, lmm_xlsx, ml_json, ml_xlsx)
ml_av    <- if (file.exists(ml_json))
              jsonlite::fromJSON(ml_json, simplifyVector = FALSE)$clinical_immune_added_value else NULL

if (!all(file.exists(needed))) {
  message(sprintf("[PostHoc] SKIPPED for '%s': missing step-01/04/06 artifact(s): %s", exp,
                  paste(basename(needed[!file.exists(needed)]), collapse = ", ")))
} else if (is.null(ml_av)) {
  message(sprintf("[PostHoc] SKIPPED for '%s': step 06 wrote no clinical_immune_added_value node (empty gate?)", exp))
} else {

  out_root <- step_dir(config, 7, create = TRUE)
  DATA      <- readRDS(rds_std)
  DATA_LONG <- readRDS(rds_long)
  J4   <- jsonlite::fromJSON(lmm_json)
  gate <- as.character(unlist(ml_av$composite_markers))
  seed <- if (!is.null(config$stats$seed)) as.integer(config$stats$seed) else 2026L
  ml_cfg <- if (!is.null(config$machine_learning)) config$machine_learning else list()
  # Clinical model exactly as step 06 used it (NLR toggle applied there; recorded in the JSON).
  clin_vars <- unlist(ml_av$clinical_vars)
  read_input <- function(f) {
    if (is.null(f) || !file.exists(f)) stop(sprintf("[FATAL] post hoc: input file not found: %s", f))
    x <- as.data.frame(readxl::read_excel(f, sheet = 1)); names(x) <- trimws(names(x))
    x$Patient_ID <- as.character(x$Patient_ID); x
  }
  # Each analysis writes a console.txt next to its tables.
  open_log <- function(dir) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    con <- file(file.path(dir, "console.txt"), open = "wt")
    list(con = con,
         say  = function(...) { t <- sprintf(...); message(t); cat(t, "\n", sep = "", file = con) },
         show = function(x) invisible(lapply(capture.output(print(x)), function(l) {
                  message("  ", l); cat("  ", l, "\n", sep = "", file = con) })))
  }
  run_sha <- tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
                      error = function(e) NA_character_)

  # --- 7a. Response grouping (PR / SD / PD) ---------------------------------------------
  g_cfg <- ph_cfg$grouping
  if (!isFALSE(g_cfg$enabled)) {
    d <- file.path(out_root, "grouping"); L <- open_log(d)
    L$say("response grouping sensitivity  %s | exp %s | LMM frame of step 04", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), exp)
    L$say("inputs (run's own): %s | %s", config$input_file_t0, config$input_file_t1)
    sub <- if (!is.null(g_cfg$subcodes)) unlist(g_cfg$subcodes) else c(PR = 2, SD = 3, PD = 4)
    covs <- as.character(unlist(config$clinical$sensitivity_covariates))
    gr <- run_response_grouping_sensitivity(
      DATA_LONG, read_input(config$input_file_t0), read_input(config$input_file_t1),
      target_col = config$clinical$target_column, subcodes = sub,
      responder_label = config$clinical$responder_label, covs = covs,
      step04_full = J4$full_results, report_markers = gate, say = L$say)
    n_mk <- gr$n_markers
    utils::write.csv(gr$unadjusted, file.path(d, sprintf("grouping_lmm_all%d_unadjusted.csv", n_mk)), row.names = FALSE)
    if (!is.null(gr$adjusted))
      utils::write.csv(gr$adjusted, file.path(d, sprintf("grouping_lmm_all%d_adj_PDL1_therapy.csv", n_mk)), row.names = FALSE)
    saveRDS(gr, file.path(d, "grouping_sensitivity.rds"))
    close(L$con)
    message(sprintf("   [PostHoc] grouping -> %s", d))
  }

  # --- 7b. Duplicated T1 rows ---------------------------------------------------------
  dr_cfg <- ph_cfg$duplicate_rows
  if (!isFALSE(dr_cfg$enabled) && length(unlist(dr_cfg$drop_t1))) {
    d <- file.path(out_root, "duplicate_rows"); L <- open_log(d)
    ids <- as.character(unlist(dr_cfg$drop_t1))
    L$say("duplicated-row sensitivity  %s | exp %s | drop T1 of: %s | gate fixed: %s",
          format(Sys.time(), "%Y-%m-%d %H:%M:%S"), exp, paste(ids, collapse = ", "), paste(gate, collapse = ", "))
    dup <- run_duplicate_row_sensitivity(DATA, DATA_LONG, gate, ids, config, clin_vars, seed = seed, say = L$say)
    if (dup$lmm_patients[1] != J4$n_patients || dup$lmm_obs[1] != J4$n_observations ||
        dup$delta_n[1] != ml_av$n)
      stop(sprintf("[FATAL] duplicate rows: 'as is' frame %d/%d, delta %d != step 04 %d/%d, step 06 %d",
                   dup$lmm_patients[1], dup$lmm_obs[1], dup$delta_n[1], J4$n_patients, J4$n_observations, ml_av$n))
    utils::write.csv(dup, file.path(d, "duplicate_rows_scenarios.csv"), row.names = FALSE)
    close(L$con)
    message(sprintf("   [PostHoc] duplicate rows -> %s", d))
  }

  # --- 7c. Therapy (ICI monotherapy vs chemo-immunotherapy) ----------------------------
  th_cfg <- ph_cfg$therapy
  if (!isFALSE(th_cfg$enabled) && !is.null(th_cfg$file_t0)) {
    d <- file.path(out_root, "therapy"); L <- open_log(d)
    th_files <- c(T0 = th_cfg$file_t0, T1 = th_cfg$file_t1)
    if (!all(file.exists(th_files))) stop(sprintf("[FATAL] therapy: file(s) not found: %s",
                                                  paste(th_files[!file.exists(th_files)], collapse = ", ")))
    B_th   <- if (!is.null(th_cfg$B)) as.integer(th_cfg$B) else 10000L
    seed_t <- if (!is.null(th_cfg$seed)) as.integer(th_cfg$seed) else 56L
    tag    <- if (!is.null(th_cfg$tag)) as.character(th_cfg$tag) else NA_character_
    sheets <- lapply(th_files, function(f) {
      x <- readxl::read_excel(f)
      data.frame(Patient_ID = as.character(x$Patient_ID), code = x[[grep("^Terapia_Pembro", names(x))]],
                 regime = x$Terapia_regime_IOmono_1_chemoIO_2, radio = x$Terapia_radio_tra_T0_T1_annotata)
    })
    df <- cbind(DATA_LONG$metadata, as.data.frame(DATA_LONG$hybrid_data_raw[, DATA_LONG$hybrid_markers]))
    df$Patient_ID <- as.character(df$Patient_ID)
    fr <- as.data.frame(readxl::read_excel(lmm_xlsx, "LMM_Interaction_Results"))   # full precision
    gate_th <- fr$Marker[fr$FDR_Interaction < 0.05]
    L$say("therapy sensitivity  %s | exp %s | therapy data %s (tag %s) | B %d | seed %d",
          format(Sys.time(), "%Y-%m-%d %H:%M:%S"), exp, dirname(th_files[[1]]), tag, B_th, seed_t)
    L$say("gate (step-04 FDR < 0.05, fixed): %s", paste(gate_th, collapse = ", "))

    # reproduction of step 04: counts and gate beta / p
    n_pts <- length(unique(df$Patient_ID)); n_obs <- nrow(df); n_pair <- sum(table(df$Patient_ID) == 2)
    L$say("frame %d patients / %d obs / %d paired (step 04: %d / %d / %d)", n_pts, n_obs, n_pair,
          J4$n_patients, J4$n_observations, J4$n_paired)
    if (!(n_pts == J4$n_patients && n_obs == J4$n_observations && n_pair == J4$n_paired))
      stop("[FATAL] therapy: LMM frame counts differ from step 04")
    for (m in gate_th) {
      r <- fit_feature_lmm(df, m); ref <- fr[fr$Marker == m, ]
      L$say("%s base model: beta %.4f p %.5f (step 04: %.4f, %.5f)", m, r$Estimate_Interaction,
            r$P_Value_Interaction, ref$Estimate_Interaction, ref$P_Value_Interaction)
      if (!(abs(r$Estimate_Interaction - ref$Estimate_Interaction) < 1e-6 &&
            abs(r$P_Value_Interaction - ref$P_Value_Interaction) < 1e-8))
        stop(sprintf("[FATAL] therapy: base model does not reproduce step 04 for %s", m))
    }

    reg <- therapy_regimen_table(sheets, df)
    res <- run_therapy_sensitivity(df, gate_th, reg, B = B_th, seed = seed_t)
    L$say("T1-only patients (regimen from the T1 sheet): %s", paste(reg$Patient_ID[reg$first_sheet == "T1"], collapse = ", "))
    L$say("without regimen: %s", paste(res$counts$no_regimen, collapse = ", "))
    L$say("radiotherapy between T0 and T1 (annotated): %s", paste(res$counts$radio, collapse = ", "))
    L$say("regimen x response, all patients:"); L$show(res$counts$rx_by_group)
    L$say("regimen x response, paired:");        L$show(res$counts$rx_by_group_paired)
    tt <- res$tests
    L$say("\ntests (estimate: beta logit, PR - SD/PD; LRT chi2; F):")
    L$show(transform(tt, estimate = signif(estimate, 4), se = signif(se, 3), p = signif(p, 4)))
    adj <- do.call(rbind, lapply(names(res$markers), function(m) cbind(marker = m, res$markers[[m]]$adjusted_change$cells)))
    L$say("\nbaseline-adjusted mean change (at mean T0) [95%% CI], raw mean, mean T0:")
    L$show(transform(adj, fit = round(fit, 3), lwr = round(lwr, 3), upr = round(upr, 3),
                     raw_mean = round(raw_mean, 3), T0_mean = round(T0_mean, 3)))

    saveRDS(res, file.path(d, "therapy_sensitivity.rds"))
    utils::write.csv(tt, file.path(d, "therapy_tests.csv"), row.names = FALSE)
    utils::write.csv(adj, file.path(d, "therapy_adjusted_change.csv"), row.names = FALSE)
    save_pub_figure(pub_fig_therapy_change(res, adjusted = TRUE),
                    file.path(d, "FigureT_TherapyChange"), PUB_W1 + 20, 85, config = config)
    save_pub_figure(pub_fig_therapy_change(res, adjusted = FALSE),
                    file.path(d, "FigureT_TherapyChange_unadjusted"), PUB_W1 + 20, 85, config = config)
    jsonlite::write_json(list(
      step = "src/07_post_hoc.R (therapy)", date = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      experiment = exp, tag = tag, git_sha = if (length(run_sha)) run_sha[1] else NA_character_,
      therapy_files = as.list(th_files), therapy_md5 = as.list(setNames(unname(tools::md5sum(th_files)), names(th_files))),
      parser = "data/therapy_fix/database_parser_nsclc_v3_therapy.R",
      B = B_th, seed = res$seed, gate = gate_th,
      frame = list(n_pts = n_pts, n_obs = n_obs, n_paired = n_pair),
      regimen_known = list(n_pts = sum(!is.na(reg$Rx)), no_regimen = res$counts$no_regimen,
                           radio = res$counts$radio),
      step04_reproduced = TRUE), file.path(d, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
    close(L$con)
    message(sprintf("   [PostHoc] therapy -> %s", d))
  }

  # --- 7d. Kaplan-Meier by median split of the Δ composite -----------------------------
  km_cfg <- ph_cfg$km
  if (!isFALSE(km_cfg$enabled)) {
    d <- file.path(out_root, "km"); L <- open_log(d)
    fatal <- function(...) stop(sprintf("[FATAL] km: %s", sprintf(...)), call. = FALSE)
    lm_days <- if (!is.null(km_cfg$landmark_days)) as.numeric(km_cfg$landmark_days) else 21
    if (!is.finite(lm_days) || lm_days < 0) fatal("post_hoc.km.landmark_days must be >= 0")
    lm_months <- lm_days / 30.4375
    sv <- ml_cfg$clinical_model$survival
    SURV <- c(os = if (!is.null(sv$os)) sv$os else "OS", pfs = if (!is.null(sv$pfs)) sv$pfs else "PFS",
              dead = if (!is.null(sv$death)) sv$death else "Alive_0/Dead_1")
    av <- ml_av
    if (!identical(av$composite_timepoint, "delta")) fatal("added-value composite is not Δ")
    pts <- as.data.frame(readxl::read_excel(ml_xlsx, sheet = "ClinImmune_Patients"))
    ids <- as.character(pts$Patient_ID)
    L$say("km median split  %s | exp %s | landmark %g days (%.3f months)",
          format(Sys.time(), "%Y-%m-%d %H:%M:%S"), exp, lm_days, lm_months)
    L$say("composite: mean Δz of the fixed gate (%s), as step 06", paste(gate, collapse = ", "))

    # Check 1: the frame
    if (length(ids) != av$n || anyDuplicated(ids)) fatal("ClinImmune_Patients has %d IDs, JSON n = %d", length(ids), av$n)
    if (sum(pts$True_Group == av$positive_label) != av$n_pos || sum(pts$True_Group != av$positive_label) != av$n_neg)
      fatal("group counts differ from JSON n_pos/n_neg")
    L$say("frame: Δ, n = %d (%d %s / %d %s) = JSON: PASS", length(ids), av$n_pos, av$positive_label,
          av$n_neg, av$negative_label)

    # Survival from the experiment's T0 input (= step 06's source)
    in_t0 <- config$input_file_t0; in_t1 <- config$input_file_t1
    rd <- function(f) { x <- read_input(f)
      if (!all(SURV %in% names(x))) fatal("survival columns missing in %s", f); x }
    s0 <- rd(in_t0); i0 <- match(ids, s0$Patient_ID)
    if (anyNA(i0)) fatal("Δ-frame IDs missing from %s", in_t0)
    os  <- suppressWarnings(as.numeric(s0[[SURV[["os"]]]][i0]))
    pfs <- suppressWarnings(as.numeric(s0[[SURV[["pfs"]]]][i0]))
    dth <- suppressWarnings(as.integer(s0[[SURV[["dead"]]]][i0]))
    # PFS event exactly as step 06 (no progression date exists: DERIVED PROXY)
    EP <- list(OS  = list(time = os,  event = dth, label = "Overall survival", json = av$survival$OS,
                          def = "OS (months from start of therapy), event = Alive_0/Dead_1 == 1"),
               PFS = list(time = pfs, event = as.integer(dth == 1 | (dth == 0 & pfs < os)),
                          label = "Progression-free survival", json = av$survival$PFS,
                          def = "PFS (months from start of therapy), event = death | (alive & PFS < OS): derived proxy, as step 06"))

    # Check 2: step-06 OS/PFS Cox reproduced (clinical model + composite). Clinical
    # columns come from the ClinImmune_Patients sheet when it has them, else from the
    # T0 input (the sheet does not carry every clinical variable, e.g. NLR).
    cl <- unlist(av$clinical_vars)
    C  <- as.data.frame(lapply(cl, function(col) {
      x <- if (col %in% names(pts)) pts[[col]] else s0[[col]][i0]
      x <- suppressWarnings(as.numeric(x)); x[is.na(x)] <- median(x, na.rm = TRUE); x }))
    names(C) <- make.names(cl, unique = TRUE)
    for (e in names(EP)) {
      t_ <- EP[[e]]$time; ev_ <- EP[[e]]$event; Jx <- EP[[e]]$json
      ok <- is.finite(t_) & t_ > 0 & !is.na(ev_)
      if (sum(ok) != Jx$n || sum(ev_[ok]) != Jx$events)
        fatal("%s n/events %d/%d != JSON %d/%d", e, sum(ok), sum(ev_[ok]), Jx$n, Jx$events)
      fk <- survival::coxph(survival::Surv(t_, ev_) ~ ., data = cbind(C, comp = pts$Composite)[ok, , drop = FALSE])
      hr <- unname(exp(coef(fk))["comp"])
      if (abs(hr - Jx$immune_hr) > 5e-4) fatal("%s Cox HR %.4f != JSON immune_hr %.4f", e, hr, Jx$immune_hr)
      L$say("%s (T0 file %s): n = %d, %d events = JSON; Cox %s + composite HR %.4f = JSON %.4f: PASS",
            e, basename(in_t0), sum(ok), sum(ev_[ok]), paste(cl, collapse = " + "), hr, Jx$immune_hr)
    }

    # Data check: survival fields T0 vs T1 file (reported, not corrected)
    disc <- NULL
    if (!is.null(in_t1) && file.exists(in_t1)) {
      s1 <- rd(in_t1); i1 <- match(ids, s1$Patient_ID)
      num <- function(s, i, k) suppressWarnings(as.numeric(s[[SURV[[k]]]][i]))
      differs <- function(a, b) !(is.na(a) & is.na(b)) & (is.na(a) | is.na(b) | a != b)
      dd <- !is.na(i1) & (differs(num(s0, i0, "os"), num(s1, i1, "os")) |
                          differs(num(s0, i0, "pfs"), num(s1, i1, "pfs")) |
                          differs(num(s0, i0, "dead"), num(s1, i1, "dead")))
      disc <- data.frame(Patient_ID = ids, Group = pts$True_Group,
                         OS_T0 = num(s0, i0, "os"), OS_T1 = num(s1, i1, "os"),
                         PFS_T0 = num(s0, i0, "pfs"), PFS_T1 = num(s1, i1, "pfs"),
                         Dead_T0 = num(s0, i0, "dead"), Dead_T1 = num(s1, i1, "dead"))[dd, , drop = FALSE]
      rownames(disc) <- NULL
      L$say("\nOS/PFS/status differ between the T0 and T1 files for %d of %d Δ-frame patients (analysis uses T0, as step 06):",
            nrow(disc), length(ids))
      if (nrow(disc)) L$show(disc)
    }

    # Analysis, per endpoint
    KM <- list()
    for (e in names(EP)) {
      res <- km_median_split(EP[[e]]$time, EP[[e]]$event, pts$Composite, ids, pts$True_Group,
                             landmark_months = lm_months)
      KM[[e]] <- res; Jx <- EP[[e]]$json
      L$say("\n==== %s: %s", e, EP[[e]]$def)
      L$say("landmark %.3f months: %d patient(s) excluded (%s <= landmark)%s", lm_months, res$n_excluded, e,
            if (res$n_excluded) paste0(": ", paste(res$excluded$Patient_ID, collapse = ", ")) else "")
      L$say("median Δ composite %.4f; arms: %s", res$median_x,
            paste(sprintf("%s n = %d", names(table(res$patients$arm)), table(res$patients$arm)), collapse = ", "))
      L$say("split x best response:"); L$show(res$split_by_group)
      L$say("\nmedian %s from the landmark (months) [95%% CI]:", e)
      L$show(transform(res$medians, median_months = round(median_months, 2), lcl = round(lcl, 2), ucl = round(ucl, 2)))
      L$say("\ntests (HR = Δ > median vs Δ ≤ median):")
      tt <- res$tests
      L$show(transform(tt[, c("analysis", "test", "n_low", "events_low", "n_high", "events_high", "chisq", "p",
                              "hr", "hr_lo", "hr_hi", "cox_zph_p")],
                       chisq = round(chisq, 3), p = signif(p, 3), hr = round(hr, 3), hr_lo = round(hr_lo, 3),
                       hr_hi = round(hr_hi, 3), cox_zph_p = signif(cox_zph_p, 3)))
      L$say("\nnumber at risk:"); L$show(reshape(res$risk, idvar = "arm", timevar = "time", direction = "wide"))
      pr <- tt[tt$analysis == "primary", ]
      L$say("\ncaption values %s: Δ frame, n = %d, %d events; split at the median Δ composite (%.3f); log-rank p = %.3f; HR %.2f [%.2f, %.2f] (Δ > vs ≤ median)",
            e, nrow(res$patients), sum(res$patients$event), res$median_x, pr$p, pr$hr, pr$hr_lo, pr$hr_hi)
      L$say(paste("notes %s: post hoc and descriptive; underpowered (min detectable HR %.2f per unit at 80%% power,",
                  "step 06); a median split loses power vs the continuous Cox (step 06 HR %.2f, p %.2f). The stratified",
                  "row is secondary/post hoc: no 'prognostic marker' or 'predictive vs prognostic' wording. With",
                  "n_excluded = 0 the %g-day landmark only shifts time: tests equal the un-landmarked analysis."),
            e, Jx$min_detectable_hr_80, Jx$immune_hr, Jx$immune_p, lm_days)
      if (e == "PFS") L$say("PFS caveat (step 06): %s", Jx$endpoint_caveat)
    }
    if (!identical(KM$OS$patients$arm, KM$PFS$patients$arm)) fatal("OS and PFS arms differ (landmark exclusions?)")

    bind <- function(k) do.call(rbind, lapply(names(KM), function(e) cbind(endpoint = e, KM[[e]][[k]])))
    saveRDS(KM, file.path(d, "km_median.rds"))
    utils::write.csv(bind("tests"), file.path(d, "km_tests.csv"), row.names = FALSE)
    utils::write.csv(bind("medians"), file.path(d, "km_medians.csv"), row.names = FALSE)
    utils::write.csv(bind("risk"), file.path(d, "km_risk_table.csv"), row.names = FALSE)
    utils::write.csv(bind("patients")[, c("endpoint", "Patient_ID", "group", "x", "arm", "time_start", "time", "event")],
                     file.path(d, "km_split_patients.csv"), row.names = FALSE)
    if (!is.null(disc)) utils::write.csv(disc, file.path(d, "os_discrepancies_T0_vs_T1.csv"), row.names = FALSE)
    save_pub_figure(pub_fig_km_panels(KM, sapply(EP, `[[`, "label")),
                    file.path(d, "FigureKM_DeltaComposite_OS_PFS"), PUB_W2, 95, config = config)
    jsonlite::write_json(list(
      step = "src/07_post_hoc.R (km)", date = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      experiment = exp, git_sha = if (length(run_sha)) run_sha[1] else NA_character_,
      input_t0 = in_t0, input_t0_md5 = unname(tools::md5sum(in_t0)),
      endpoints = lapply(EP, `[[`, "def"),
      score = "Δ composite (ClinImmune_Patients$Composite)", gate = gate,
      landmark_days = lm_days, landmark_months = lm_months,
      n_excluded_landmark = lapply(KM, `[[`, "n_excluded"),
      median_score = KM$OS$median_x, n = nrow(KM$OS$patients),
      events = lapply(KM, function(r) sum(r$patients$event)),
      os_discrepancies_T0_vs_T1 = if (is.null(disc)) NA else nrow(disc),
      step06_cox_reproduced = TRUE), file.path(d, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
    close(L$con)
    message(sprintf("   [PostHoc] km -> %s", d))
  }
}

message("=== STEP 7 COMPLETE ===\n")
