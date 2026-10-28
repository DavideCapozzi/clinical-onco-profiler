#!/usr/bin/env Rscript
# tools/km_median_split.R
# ==============================================================================
# Kaplan-Meier OS and PFS (panels A, B) by median split of the Δ immune composite, on an
# EXISTING run (no re-analysis, no marker re-selection; post hoc, descriptive).
# Logic in R/modules_survival_km.R; this driver owns the paths.
#
# Usage:
#   Rscript tools/km_median_split.R --run <run> [--exp A,B] [--landmark_days 21]
#                                   [--out DIR] [--preview DIR] [--force]
#
#   --run            run root, absolute or under results/ (REQUIRED: no default run)
#   --landmark_days  T1 landmark in days from the start of therapy (default 21 = after
#                    the 1st cycle; user, 2026-10-27). OS is in months: days / 30.4375
#   --out            output dir (default <run>/<exp>/post_hoc/km_median)
#   --force          overwrite an existing output dir (default: refuse)
#
# Reads the step-06 ML JSON (frame counts, composite timepoint, OS block) and the
# ClinImmune_Patients sheet of the ML xlsx (per-patient Δ composite), plus OS, PFS and
# Alive_0/Dead_1 from the experiment's T0 input file (as step 06; PFS event is the
# step-06 derived proxy). Aborts (non-zero exit) unless the frame and the step-06 OS and
# PFS Cox (PD-L1 + composite) are reproduced.
# Never writes into the run's step dirs or publication/.
# Outputs: FigureKM_DeltaComposite_OS_PFS.pdf, km_tests.csv, km_medians.csv,
# km_risk_table.csv, km_split_patients.csv, os_discrepancies_T0_vs_T1.csv,
# km_median.rds, console.txt, manifest.json
# ==============================================================================

suppressPackageStartupMessages({ library(here); library(jsonlite); library(readxl); library(survival) })
source(here("R/utils_io.R"))
source(here("R/utils_figures.R"))
source(here("R/utils_metrics.R"))            # list_run_experiments()
source(here("R/modules_pub_style.R"))
source(here("R/modules_survival_km.R"))

parse_args <- function(args) {
  out <- list(); i <- 1
  while (i <= length(args)) {
    a <- args[[i]]
    if (startsWith(a, "--")) {
      key <- sub("^--", "", a)
      val <- if (i + 1 <= length(args) && !startsWith(args[[i + 1]], "--")) args[[i + 1]] else ""
      out[[key]] <- val
      i <- i + (if (nzchar(val)) 2 else 1)
    } else i <- i + 1
  }
  out
}
resolve_run_arg <- function(x) {
  if (dir.exists(x)) return(x)
  cand <- file.path("results", x)
  if (dir.exists(cand)) return(cand)
  stop(sprintf("[km_median_split] run not found: %s", x))
}
git_out <- function(args) tryCatch(system2("git", args, stdout = TRUE, stderr = FALSE), error = function(e) NA_character_)
`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && !nzchar(a))) b else a
fatal <- function(...) stop(sprintf("[km_median_split] FATAL: %s", sprintf(...)), call. = FALSE)

opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(opt$run) || !nzchar(opt$run)) stop("[km_median_split] --run is required")
run  <- resolve_run_arg(opt$run)
lm_days <- if (!is.null(opt$landmark_days)) as.numeric(opt$landmark_days) else 21
if (!is.finite(lm_days) || lm_days < 0) stop("[km_median_split] --landmark_days must be >= 0")
lm_months <- lm_days / 30.4375
exps <- if (!is.null(opt$exp)) trimws(strsplit(opt$exp, ",")[[1]]) else list_run_experiments(run)
if (!is.null(opt$preview)) options(pub_fig.preview_dir = opt$preview)
run_manifest <- file.path(run, "run_manifest.yml")
run_sha <- if (file.exists(run_manifest)) sub("^git_sha: *", "", grep("^git_sha:", readLines(run_manifest), value = TRUE)[1]) else NA
SURV <- c(os = "OS", pfs = "PFS", dead = "Alive_0/Dead_1")

for (exp in exps) {
  cfg  <- list(output_root = file.path(run, exp), project_name = exp)   # enough for step_output_path()
  ml_json <- step_output_path(cfg, 6, sprintf("Machine_Metrics_ML_%s", exp), "json")
  ml_xlsx <- step_output_path(cfg, 6, sprintf("ML_Classification_Report_%s", exp), "xlsx")
  lon_rds <- step_output_path(cfg, 1, sprintf("data_processed_%s_longitudinal", exp), "rds")
  if (!all(file.exists(c(ml_json, ml_xlsx, lon_rds)))) { message(sprintf("[skip] %s: no step-01/06 artifacts", exp)); next }
  out_dir <- if (!is.null(opt$out)) { if (length(exps) > 1) file.path(opt$out, exp) else opt$out }
             else file.path(run, exp, "post_hoc", "km_median")
  if (dir.exists(out_dir) && length(list.files(out_dir)) && is.null(opt$force))
    stop(sprintf("[km_median_split] %s exists; use --force", out_dir))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  con <- file(file.path(out_dir, "console.txt"), open = "wt")
  say <- function(...) { t <- sprintf(...); cat(t, "\n"); cat(t, "\n", file = con) }
  show <- function(x) invisible(lapply(capture.output(print(x)), function(l) say("  %s", l)))

  av <- fromJSON(ml_json, simplifyVector = FALSE)$clinical_immune_added_value
  if (!identical(av$composite_timepoint, "delta")) fatal("%s: added-value composite is not Δ", exp)
  gate <- unlist(av$composite_markers)
  pts  <- as.data.frame(read_excel(ml_xlsx, sheet = "ClinImmune_Patients"))
  ids  <- as.character(pts$Patient_ID)
  say("km_median_split  %s | run %s (code %s) | exp %s | landmark %g days (%.3f months)",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"), run, substr(run_sha, 1, 7), exp, lm_days, lm_months)
  say("composite: mean Δz of the fixed gate (%s), as step 06", paste(gate, collapse = ", "))

  # ---- Check 1: the frame ---------------------------------------------------------
  if (length(ids) != av$n || anyDuplicated(ids)) fatal("ClinImmune_Patients has %d IDs, JSON n = %d", length(ids), av$n)
  if (sum(pts$True_Group == av$positive_label) != av$n_pos || sum(pts$True_Group != av$positive_label) != av$n_neg)
    fatal("group counts differ from JSON n_pos/n_neg")
  say("frame: Δ, n = %d (%d %s / %d %s) = JSON: PASS", length(ids), av$n_pos, av$positive_label,
      av$n_neg, av$negative_label)

  # ---- Survival from the experiment's input files (T0 = step 06's source) ----------
  lon   <- readRDS(lon_rds)
  ecfg  <- lon$config$experiments[[exp]]
  in_t0 <- ecfg$input_file %||% lon$config$input_file
  in_t1 <- ecfg$input_file_t1 %||% lon$config$input_file_t1
  rd <- function(f) { x <- as.data.frame(read_excel(f, sheet = 1)); names(x) <- trimws(names(x))
    if (!all(SURV %in% names(x))) fatal("survival columns missing in %s", f)
    x$Patient_ID <- as.character(x$Patient_ID); x }
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

  # ---- Check 2: step-06 OS/PFS Cox reproduced (clinical model + composite) --------
  cl <- unlist(av$clinical_vars)
  C  <- as.data.frame(lapply(pts[cl], function(x) { x <- as.numeric(x); x[is.na(x)] <- median(x, na.rm = TRUE); x }))
  names(C) <- make.names(names(C), unique = TRUE)
  for (e in names(EP)) {
    t_ <- EP[[e]]$time; ev_ <- EP[[e]]$event; J <- EP[[e]]$json
    ok <- is.finite(t_) & t_ > 0 & !is.na(ev_)
    if (sum(ok) != J$n || sum(ev_[ok]) != J$events)
      fatal("%s n/events %d/%d != JSON %d/%d", e, sum(ok), sum(ev_[ok]), J$n, J$events)
    fk <- coxph(Surv(t_, ev_) ~ ., data = cbind(C, comp = pts$Composite)[ok, , drop = FALSE])
    hr <- unname(exp(coef(fk))["comp"])
    if (abs(hr - J$immune_hr) > 5e-4) fatal("%s Cox HR %.4f != JSON immune_hr %.4f", e, hr, J$immune_hr)
    say("%s (T0 file %s): n = %d, %d events = JSON; Cox %s + composite HR %.4f = JSON %.4f: PASS",
        e, basename(in_t0), sum(ok), sum(ev_[ok]), paste(cl, collapse = " + "), hr, J$immune_hr)
  }

  # ---- Data check: survival fields T0 vs T1 file (reported, not corrected) --------
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
    say("\nOS/PFS/status differ between the T0 and T1 files for %d of %d Δ-frame patients (analysis uses T0, as step 06):",
        nrow(disc), length(ids))
    if (nrow(disc)) show(disc)
  }

  # ---- Analysis, per endpoint -------------------------------------------------------
  R <- list()
  for (e in names(EP)) {
    res <- km_median_split(EP[[e]]$time, EP[[e]]$event, pts$Composite, ids, pts$True_Group,
                           landmark_months = lm_months)
    R[[e]] <- res; J <- EP[[e]]$json
    say("\n==== %s: %s", e, EP[[e]]$def)
    say("landmark %.3f months: %d patient(s) excluded (%s <= landmark)%s", lm_months, res$n_excluded, e,
        if (res$n_excluded) paste0(": ", paste(res$excluded$Patient_ID, collapse = ", ")) else "")
    say("median Δ composite %.4f; arms: %s", res$median_x,
        paste(sprintf("%s n = %d", names(table(res$patients$arm)), table(res$patients$arm)), collapse = ", "))
    say("split x best response:"); show(res$split_by_group)
    say("\nmedian %s from the landmark (months) [95%% CI]:", e)
    show(transform(res$medians, median_months = round(median_months, 2), lcl = round(lcl, 2), ucl = round(ucl, 2)))
    say("\ntests (HR = Δ > median vs Δ ≤ median):")
    tt <- res$tests
    show(transform(tt[, c("analysis", "test", "n_low", "events_low", "n_high", "events_high", "chisq", "p",
                          "hr", "hr_lo", "hr_hi", "cox_zph_p")],
                   chisq = round(chisq, 3), p = signif(p, 3), hr = round(hr, 3), hr_lo = round(hr_lo, 3),
                   hr_hi = round(hr_hi, 3), cox_zph_p = signif(cox_zph_p, 3)))
    say("\nnumber at risk:"); show(reshape(res$risk, idvar = "arm", timevar = "time", direction = "wide"))
    pr <- tt[tt$analysis == "primary", ]
    say("\ncaption values %s: Δ frame, n = %d, %d events; split at the median Δ composite (%.3f); log-rank p = %.3f; HR %.2f [%.2f, %.2f] (Δ > vs ≤ median)",
        e, nrow(res$patients), sum(res$patients$event), res$median_x, pr$p, pr$hr, pr$hr_lo, pr$hr_hi)
    say(paste("notes %s: post hoc and descriptive; underpowered (min detectable HR %.2f per unit at 80%% power,",
              "step 06); a median split loses power vs the continuous Cox (step 06 HR %.2f, p %.2f). The stratified",
              "row is secondary/post hoc: no 'prognostic marker' or 'predictive vs prognostic' wording. With",
              "n_excluded = 0 the %g-day landmark only shifts time: tests equal the un-landmarked analysis."),
        e, J$min_detectable_hr_80, J$immune_hr, J$immune_p, lm_days)
    if (e == "PFS") say("PFS caveat (step 06): %s", J$endpoint_caveat)
  }
  if (!identical(R$OS$patients$arm, R$PFS$patients$arm)) fatal("OS and PFS arms differ (landmark exclusions?)")

  bind <- function(k) do.call(rbind, lapply(names(R), function(e) cbind(endpoint = e, R[[e]][[k]])))
  saveRDS(R, file.path(out_dir, "km_median.rds"))
  write.csv(bind("tests"), file.path(out_dir, "km_tests.csv"), row.names = FALSE)
  write.csv(bind("medians"), file.path(out_dir, "km_medians.csv"), row.names = FALSE)
  write.csv(bind("risk"), file.path(out_dir, "km_risk_table.csv"), row.names = FALSE)
  write.csv(bind("patients")[, c("endpoint", "Patient_ID", "group", "x", "arm", "time_start", "time", "event")],
            file.path(out_dir, "km_split_patients.csv"), row.names = FALSE)
  if (!is.null(disc)) write.csv(disc, file.path(out_dir, "os_discrepancies_T0_vs_T1.csv"), row.names = FALSE)
  save_pub_figure(pub_fig_km_panels(R, sapply(EP, `[[`, "label")),
                  file.path(out_dir, "FigureKM_DeltaComposite_OS_PFS.pdf"), PUB_W2, 95)

  write_json(list(
    tool = "tools/km_median_split.R", date = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    run = run, run_git_sha = run_sha, experiment = exp,
    git_sha = git_out(c("rev-parse", "HEAD")),
    git_dirty = length(git_out(c("status", "--porcelain", "--", "R/modules_survival_km.R",
                                 "tools/km_median_split.R"))) > 0,
    input_t0 = in_t0, input_t0_md5 = unname(tools::md5sum(in_t0)),
    endpoints = lapply(EP, `[[`, "def"),
    score = "Δ composite (ClinImmune_Patients$Composite)", gate = gate,
    landmark_days = lm_days, landmark_months = lm_months,
    n_excluded_landmark = lapply(R, `[[`, "n_excluded"),
    median_score = R$OS$median_x, n = nrow(R$OS$patients),
    events = lapply(R, function(r) sum(r$patients$event)),
    os_discrepancies_T0_vs_T1 = if (is.null(disc)) NA else nrow(disc),
    step06_cox_reproduced = TRUE), file.path(out_dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
  say("\nwritten: %s", out_dir)
  close(con)
}
