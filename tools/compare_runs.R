#!/usr/bin/env Rscript
# tools/compare_runs.R
# ==============================================================================
# Side-by-side table of the numbers each manuscript claim rests on, for several runs
# (step-04 xlsx, step-06 JSON, the step-06 patient-bootstrap and selection-aware rds;
# for runs that predate their integration, the old diagnostics/ rds). Read-only.
# Missing pieces print "—". The gate is read from each run's step-06 JSON.
# Usage: Rscript tools/compare_runs.R label=run_dir [label=run_dir ...] [--exp EXP] [--out file.md]
# ==============================================================================
suppressMessages({ library(readxl); library(jsonlite) })
args <- commandArgs(TRUE)
out_file <- if ("--out" %in% args) args[which(args == "--out") + 1] else NULL
EXP  <- if ("--exp" %in% args) args[which(args == "--exp") + 1] else "BestResponse_2v3_4"
args <- args[grepl("=", args)]
if (!length(args)) stop("[compare_runs] give at least one label=run_dir")
RUNS <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
for (r in RUNS) if (!dir.exists(r)) stop(sprintf("[compare_runs] run not found: %s", r))

f3 <- function(x, d = 3) if (length(x) == 0 || all(is.na(x))) "—" else formatC(x, format = "f", digits = d)
fp <- function(x) if (length(x) == 0 || is.na(x)) "—" else if (x < 1e-3) formatC(x, format = "e", digits = 2) else formatC(x, format = "f", digits = 4)
ci <- function(v, d = 3) if (length(v) < 2 || anyNA(v)) "—" else sprintf("[%s, %s]", f3(v[1], d), f3(v[2], d))
pair <- function(a, b, f = f3) sprintf("%s / %s", f(a), f(b))
sheet <- function(x, s) tryCatch(as.data.frame(read_excel(x, sheet = s)), error = function(e) NULL)
row_of <- function(df, mk, col) if (is.null(df) || !col %in% names(df) || !mk %in% df$Marker) NA else df[[col]][df$Marker == mk][1]

collect <- function(rd) {
  m <- list()
  jf <- Sys.glob(file.path(rd, EXP, "06_*", sprintf("Machine_Metrics_ML_%s.json", EXP)))
  J  <- if (length(jf)) fromJSON(jf) else NULL
  GATE <- if (!is.null(J$clinical_immune_added_value$composite_markers))
            as.character(J$clinical_immune_added_value$composite_markers) else c(NA_character_, NA_character_)
  if (length(GATE) < 2) GATE <- c(GATE, NA_character_)
  m$gate      <- paste(GATE[!is.na(GATE)], collapse = " / ")
  x4 <- Sys.glob(file.path(rd, EXP, "04_*", sprintf("Longitudinal_LMM_Report_%s.xlsx", EXP)))
  L  <- sheet(x4, "LMM_Interaction_Results")
  m$n_fdr     <- if (is.null(L)) NA else sum(L$FDR_Interaction < 0.05, na.rm = TRUE)
  m$beta      <- pair(row_of(L, GATE[1], "Estimate_Interaction"), row_of(L, GATE[2], "Estimate_Interaction"))
  m$p         <- pair(row_of(L, GATE[1], "P_Value_Interaction"), row_of(L, GATE[2], "P_Value_Interaction"), fp)
  m$fdr       <- pair(row_of(L, GATE[1], "FDR_Interaction"), row_of(L, GATE[2], "FDR_Interaction"), fp)
  nxt <- if (is.null(L)) NULL else L[!L$Marker %in% GATE, ][1, ]
  m$next_fdr  <- if (is.null(nxt)) "—" else sprintf("%s %s", nxt$Marker, fp(nxt$FDR_Interaction))
  m$loo       <- pair(row_of(L, GATE[1], "Max_P_Value_LOO"), row_of(L, GATE[2], "Max_P_Value_LOO"), fp)
  R <- sheet(x4, "LMM_RankANCOVA_Confirm")
  m$rank      <- pair(row_of(R, GATE[1], "RankANCOVA_Perm_P"), row_of(R, GATE[2], "RankANCOVA_Perm_P"), fp)
  P <- sheet(x4, "LMM_Paired_Delta_Sensitivity")
  m$paired    <- pair(row_of(P, GATE[1], "P_Value"), row_of(P, GATE[2], "P_Value"), fp)
  B <- sheet(x4, "LMM_Bootstrap_CI")
  m$boot      <- pair(row_of(B, GATE[1], "Pct_FDR_Significant"), row_of(B, GATE[2], "Pct_FDR_Significant"), function(v) f3(v, 1))
  C <- sheet(x4, "LMM_Covariate_Sensitivity")
  m$cov_fdr   <- pair(row_of(C, GATE[1], "FDR_Adj"), row_of(C, GATE[2], "FDR_Adj"), fp)

  A  <- J$clinical_immune_added_value
  cv <- A$cv_kfold
  m$auc_cv    <- if (is.null(cv)) "—" else sprintf("%s / %s / %s", f3(cv$auc_clinical), f3(cv$auc_immune), f3(cv$auc_combined))
  m$dauc      <- if (is.null(cv)) "—" else f3(cv$delta_auc)
  lp <- if (is.null(A)) NA else suppressWarnings(as.numeric(A$increment$lrt_p_sci))   # lrt_p itself is rounded
  m$lrt       <- if (is.na(lp)) "—" else sprintf("%s (p %s)", f3(qchisq(lp, 1, lower.tail = FALSE), 2), fp(lp))
  m$idi       <- if (is.null(A)) "—" else sprintf("%s %s", f3(A$increment$idi), ci(A$increment$idi_ci))
  os <- A$survival$OS
  m$os        <- if (is.null(os)) "—" else sprintf("HR %s %s, p %s", f3(os$immune_hr, 2), ci(os$immune_hr_ci, 2), f3(os$immune_p, 2))
  cu <- J$clinical_utility$discrimination
  m$t0_util   <- if (is.null(cu)) "—" else sprintf("app %s, LOO %s, perm p %s", f3(cu$auc_apparent), f3(cu$auc_loo), f3(cu$loo_perm_p))
  m$classif   <- if (is.null(J)) "—" else sprintf("EN %s / SVM %s / nested %s", f3(J$elastic_net$metrics$auc), f3(J$svm_rbf$metrics$auc), f3(J$nested_loocv_validation$auc))

  first <- function(p) { p <- p[file.exists(p)]; if (length(p)) p[1] else NA_character_ }
  pb <- first(c(Sys.glob(file.path(rd, EXP, "06_*", sprintf("PatientBootstrap_AUC_%s.rds", EXP))),
                file.path(rd, "diagnostics", "diag_45_bootstrap_delta_auc.rds"),
                file.path(rd, "diagnostics", "top2", "diag_45_bootstrap_delta_auc.rds")))
  g <- if (!is.na(pb)) { r <- readRDS(pb)$results; if (!is.null(r$gate)) r$gate else r$gate_2marker } else NULL
  m$ci_arms   <- if (is.null(g)) "—" else sprintf("%s / %s / %s", ci(g$bootA$auc$clinical), ci(g$bootA$auc$immune), ci(g$bootA$auc$combined))
  m$dauc_ci   <- if (is.null(g)) "—" else ci(g$bootA$delta_auc)
  m$dauc_imm  <- if (is.null(g)) "—" else sprintf("%s %s", f3(g$point$delta_auc_vs_immune), ci(g$bootA$delta_auc_vs_immune))
  sa <- first(c(Sys.glob(file.path(rd, EXP, "06_*", sprintf("SelectionAware_LRT_%s.rds", EXP))),
                file.path(rd, "diagnostics", "top2", "diag_44_selection_aware_lrt.rds")))
  s <- if (!is.na(sa)) readRDS(sa) else NULL
  nul <- if (is.null(s)) NULL else if (!is.null(s$lr_null_topk)) s$lr_null_topk else s$lr_null3
  m$sel_p     <- if (is.null(s)) "—" else { k <- sum(nul >= s$lr_obs)
    sprintf("%s (%d/%d)", f3((1 + k) / (length(nul) + 1), 4), k, length(nul)) }
  m$null95    <- if (is.null(s)) "—" else sprintf("%s vs %s", f3(quantile(nul, .95), 2), f3(s$lr_obs, 2))
  L1 <- J$lrt_leave_one_out
  m$lrt_loo   <- if (is.null(L1)) "—" else sprintf("%s (−%s) to %s; %s below null 95th pct",
                   f3(L1$lr_min, 2), L1$lr_min_dropped, f3(L1$lr_max, 2), L1$n_below_null_q95)
  N1 <- J$clinical_immune_added_value
  m$frames    <- sprintf("T0 %s · LMM %s/%s · Δ %s", if (is.null(J$clinical_utility$n)) "—" else J$clinical_utility$n,
                   if (length(x4)) tryCatch(fromJSON(Sys.glob(file.path(rd, EXP, "04_*", sprintf("Machine_Metrics_LMM_%s.json", EXP))))$n_patients, error = function(e) "—") else "—",
                   if (length(x4)) tryCatch(fromJSON(Sys.glob(file.path(rd, EXP, "04_*", sprintf("Machine_Metrics_LMM_%s.json", EXP))))$n_observations, error = function(e) "—") else "—",
                   if (is.null(N1$n)) "—" else N1$n)
  m
}

LAB <- c(frames = "frames: T0 · LMM patients/obs · Δ", gate = "gate (step-06 composite)",
         n_fdr = "Gate: markers at FDR < 0.05 (LMM frame)", beta = "β (gate markers, in order)", p = "p", fdr = "FDR",
         next_fdr = "next marker (FDR)", loo = "LOO max p", rank = "rank ANCOVA perm p", paired = "paired Δ OLS p",
         boot = "bootstrap % FDR < 0.05 (B 500)", cov_fdr = "FDR adjusted PD-L1 + therapy",
         auc_cv = "AUC clinical / immune / combined (Δ, CV)", ci_arms = "per-arm 95% CI (patient bootstrap)",
         dauc = "ΔAUC combined − clinical", dauc_ci = "ΔAUC 95% CI (patient bootstrap)",
         dauc_imm = "ΔAUC combined − immune [CI]", lrt = "LRT χ² (naive p)", sel_p = "selection-aware p, top-K (k ≥ obs / B)",
         null95 = "null 95th pct vs observed LR", lrt_loo = "fixed-gate LRT leave-one-out",
         idi = "IDI [CI]", os = "OS: immune HR [CI], p",
         t0_util = "T0 composite (T0 frame)", classif = "classifiers: EN / SVM / nested LOOCV AUC")
M <- lapply(RUNS, collect)
lines <- c(sprintf("| | %s |", paste(names(RUNS), collapse = " | ")),
           sprintf("|---|%s", paste(rep("---|", length(RUNS)), collapse = "")))
for (k in names(LAB)) lines <- c(lines, sprintf("| %s | %s |", LAB[[k]],
                                                paste(vapply(M, function(m) as.character(m[[k]]), ""), collapse = " | ")))
cat(lines, sep = "\n")
if (!is.null(out_file)) writeLines(c(sprintf("Runs: %s", paste(sprintf("%s = `%s`", names(RUNS), RUNS), collapse = "; ")), "", lines), out_file)
