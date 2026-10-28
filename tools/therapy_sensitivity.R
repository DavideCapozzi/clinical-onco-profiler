#!/usr/bin/env Rscript
# tools/therapy_sensitivity.R
# ==============================================================================
# Therapy sensitivity of the LMM gate on an EXISTING run (no re-analysis, no marker
# re-selection): is the response-associated change explained by concomitant
# chemotherapy? Logic in R/modules_therapy_sensitivity.R; this driver owns the paths.
#
# Usage:
#   Rscript tools/therapy_sensitivity.R --run <run> --tag <tag> [--therapy DIR] [--exp A,B]
#                                       [--B 10000] [--out DIR] [--preview DIR] [--force]
#
#   --run      run root, absolute or under results/ (REQUIRED: no default run)
#   --tag      label of the therapy data version, e.g. v3.1 (REQUIRED; names the output dir)
#   --therapy  dir with Dati_NSCLC_standardizzati_latest_anonimi_therapyfix_{T0,T1}.xlsx
#              (default data/therapy_fix, written by database_parser_nsclc_v3_therapy.R)
#   --out      output dir (default <run>/<exp>/post_hoc/therapy_<tag>)
#   --force    overwrite an existing output dir (default: refuse)
#
# Reads the run's step-01 longitudinal rds and step-04 LMM JSON (counts) + xlsx report
# (gate, beta, p at full precision; the JSON rounds to 4 dp) through step_dir();
# never writes into the run's step dirs or publication/. Aborts (non-zero exit) unless
# the base model reproduces the step-04 gate (beta, p) and the frame counts.
# Outputs: therapy_sensitivity.rds, therapy_tests.csv, therapy_adjusted_change.csv,
# console.txt, manifest.json, FigureT_TherapyChange{,_unadjusted}.pdf
# ==============================================================================

suppressPackageStartupMessages({ library(here); library(jsonlite); library(readxl); library(lmerTest) })
source(here("R/utils_io.R"))
source(here("R/utils_figures.R"))
source(here("R/utils_metrics.R"))            # list_run_experiments()
source(here("R/modules_longitudinal.R"))     # fit_feature_lmm()
source(here("R/modules_pub_style.R"))
source(here("R/modules_therapy_sensitivity.R"))

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
  stop(sprintf("[therapy_sensitivity] run not found: %s", x))
}
git_out <- function(args) tryCatch(system2("git", args, stdout = TRUE, stderr = FALSE), error = function(e) NA_character_)

opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(opt$run) || !nzchar(opt$run)) stop("[therapy_sensitivity] --run is required")
if (is.null(opt$tag) || !nzchar(opt$tag)) stop("[therapy_sensitivity] --tag is required")
run   <- resolve_run_arg(opt$run)
thdir <- if (!is.null(opt$therapy)) opt$therapy else "data/therapy_fix"
B     <- if (!is.null(opt$B)) as.integer(opt$B) else 10000L
exps  <- if (!is.null(opt$exp)) trimws(strsplit(opt$exp, ",")[[1]]) else list_run_experiments(run)
if (!is.null(opt$preview)) options(pub_fig.preview_dir = opt$preview)

th_files <- setNames(file.path(thdir, sprintf("Dati_NSCLC_standardizzati_latest_anonimi_therapyfix_%s.xlsx",
                                              c("T0", "T1"))), c("T0", "T1"))
stopifnot(all(file.exists(th_files)))
sheets <- lapply(th_files, function(f) {
  x <- read_excel(f)
  data.frame(Patient_ID = as.character(x$Patient_ID), code = x[[grep("^Terapia_Pembro", names(x))]],
             regime = x$Terapia_regime_IOmono_1_chemoIO_2, radio = x$Terapia_radio_tra_T0_T1_annotata)
})
run_manifest <- file.path(run, "run_manifest.yml")
run_sha <- if (file.exists(run_manifest)) sub("^git_sha: *", "", grep("^git_sha:", readLines(run_manifest), value = TRUE)[1]) else NA

for (exp in exps) {
  cfg  <- list(output_root = file.path(run, exp), project_name = exp)   # enough for step_dir()
  out_dir <- if (!is.null(opt$out)) { if (length(exps) > 1) file.path(opt$out, exp) else opt$out }
             else file.path(run, exp, "post_hoc", paste0("therapy_", opt$tag))
  if (dir.exists(out_dir) && length(list.files(out_dir)) && is.null(opt$force))
    stop(sprintf("[therapy_sensitivity] %s exists; use a new --tag or --force", out_dir))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  con <- file(file.path(out_dir, "console.txt"), open = "wt")
  say <- function(...) { t <- sprintf(...); cat(t, "\n"); cat(t, "\n", file = con) }
  show <- function(x) invisible(lapply(capture.output(print(x)), function(l) say("  %s", l)))

  rds  <- file.path(step_dir(cfg, 1), sprintf("data_processed_%s_longitudinal.rds", exp))
  lmmj <- file.path(step_dir(cfg, 4), sprintf("Machine_Metrics_LMM_%s.json", exp))
  lmmx <- file.path(step_dir(cfg, 4), sprintf("Longitudinal_LMM_Report_%s.xlsx", exp))
  if (!all(file.exists(c(rds, lmmj, lmmx)))) { message(sprintf("[skip] %s: no step-01/04 artifacts", exp)); close(con); next }
  DL <- readRDS(rds)
  df <- cbind(DL$metadata, as.data.frame(DL$hybrid_data_raw[, DL$hybrid_markers]))
  df$Patient_ID <- as.character(df$Patient_ID)
  J  <- fromJSON(lmmj)                                                   # counts
  fr <- as.data.frame(read_excel(lmmx, "LMM_Interaction_Results"))     # full precision (JSON rounds to 4 dp)
  gate <- fr$Marker[fr$FDR_Interaction < 0.05]
  say("therapy_sensitivity  %s | run %s (code %s) | exp %s | therapy %s (tag %s) | B %d",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"), run, substr(run_sha, 1, 7), exp, thdir, opt$tag, B)
  say("gate (step-04 FDR < 0.05, fixed): %s", paste(gate, collapse = ", "))

  # ---- reproduction of step 04: counts and gate beta / p ------------------------
  n_pts <- length(unique(df$Patient_ID)); n_obs <- nrow(df); n_pair <- sum(table(df$Patient_ID) == 2)
  say("frame %d patients / %d obs / %d paired (step 04: %d / %d / %d)", n_pts, n_obs, n_pair,
      J$n_patients, J$n_observations, J$n_paired)
  stopifnot(n_pts == J$n_patients, n_obs == J$n_observations, n_pair == J$n_paired)
  for (m in gate) {
    r <- fit_feature_lmm(df, m); ref <- fr[fr$Marker == m, ]
    say("%s base model: beta %.4f p %.5f (step 04: %.4f, %.5f)", m, r$Estimate_Interaction,
        r$P_Value_Interaction, ref$Estimate_Interaction, ref$P_Value_Interaction)
    stopifnot(abs(r$Estimate_Interaction - ref$Estimate_Interaction) < 1e-6,
              abs(r$P_Value_Interaction - ref$P_Value_Interaction) < 1e-8)
  }

  # ---- analysis -------------------------------------------------------------------
  reg <- therapy_regimen_table(sheets, df)
  res <- run_therapy_sensitivity(df, gate, reg, B = B, seed = 56)
  say("T1-only patients (regimen from the T1 sheet): %s", paste(reg$Patient_ID[reg$first_sheet == "T1"], collapse = ", "))
  say("without regimen: %s", paste(res$counts$no_regimen, collapse = ", "))
  say("radiotherapy between T0 and T1 (annotated): %s", paste(res$counts$radio, collapse = ", "))
  say("regimen x response, all patients:"); show(res$counts$rx_by_group)
  say("regimen x response, paired:");        show(res$counts$rx_by_group_paired)
  tt <- res$tests
  say("\ntests (estimate: beta logit, PR - SD/PD; LRT chi2; F):")
  show(transform(tt, estimate = signif(estimate, 4), se = signif(se, 3), p = signif(p, 4)))
  adj <- do.call(rbind, lapply(names(res$markers), function(m) cbind(marker = m, res$markers[[m]]$adjusted_change$cells)))
  say("\nbaseline-adjusted mean change (at mean T0) [95%% CI], raw mean, mean T0:")
  show(transform(adj, fit = round(fit, 3), lwr = round(lwr, 3), upr = round(upr, 3),
                 raw_mean = round(raw_mean, 3), T0_mean = round(T0_mean, 3)))

  saveRDS(res, file.path(out_dir, "therapy_sensitivity.rds"))
  write.csv(tt, file.path(out_dir, "therapy_tests.csv"), row.names = FALSE)
  write.csv(adj, file.path(out_dir, "therapy_adjusted_change.csv"), row.names = FALSE)
  save_pub_figure(pub_fig_therapy_change(res, adjusted = TRUE),
                  file.path(out_dir, "FigureT_TherapyChange.pdf"), PUB_W1 + 20, 85)
  save_pub_figure(pub_fig_therapy_change(res, adjusted = FALSE),
                  file.path(out_dir, "FigureT_TherapyChange_unadjusted.pdf"), PUB_W1 + 20, 85)

  md5 <- unname(tools::md5sum(th_files))
  write_json(list(
    tool = "tools/therapy_sensitivity.R", date = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    run = run, run_git_sha = run_sha, experiment = exp, tag = opt$tag,
    git_sha = git_out(c("rev-parse", "HEAD")),
    git_dirty = length(git_out(c("status", "--porcelain", "--", "R/modules_therapy_sensitivity.R",
                                 "tools/therapy_sensitivity.R"))) > 0,
    therapy_files = as.list(setNames(th_files, c("T0", "T1"))), therapy_md5 = as.list(setNames(md5, c("T0", "T1"))),
    parser = "data/therapy_fix/database_parser_nsclc_v3_therapy.R",
    B = B, seed = res$seed, gate = gate,
    frame = list(n_pts = n_pts, n_obs = n_obs, n_paired = n_pair),
    regimen_known = list(n_pts = sum(!is.na(reg$Rx)), no_regimen = res$counts$no_regimen,
                         radio = res$counts$radio),
    step04_reproduced = TRUE), file.path(out_dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
  say("\nwritten: %s", out_dir)
  close(con)
}
