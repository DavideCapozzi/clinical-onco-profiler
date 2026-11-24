#!/usr/bin/env Rscript
# tests/test_added_value_pairs.R
# ==============================================================================
# Unit test, on a SYNTHETIC added-value frame, of the "clinical var + immune score" arms
# (R/modules_ml_utility.R: av_cv_kfold(pairs, curves), av_attach_pair_curves();
# R/modules_ml_bootstrap_ci.R: patient_bootstrap_table(); R/modules_pub_figures.R:
# pub_fig_pair_roc()):
#   - adding the pair arms and the curves changes no other number (same folds, no RNG draw)
#     and leaves the caller's RNG state untouched;
#   - each averaged ROC curve has the area of its mean CV AUC;
#   - av_attach_pair_curves() adds the arms only when every existing arm is reproduced, and
#     stops on any mismatch (curve, AUC, seed);
#   - the bootstrap table carries the pair rows; the figure has one panel per clinical var.
#
# Run:  Rscript tests/test_added_value_pairs.R      (from the repository root, ~45 s)
# Exit: 0 on pass, 1 on failure.
# ==============================================================================
suppressPackageStartupMessages({
  library(here)
  source(here("R/utils_io.R")); source(here("R/utils_figures.R"))
  source(here("R/modules_ml_utility.R")); source(here("R/modules_ml_bootstrap_ci.R"))
  source(here("R/modules_pub_style.R")); source(here("R/modules_pub_figures.R"))
})

fails <- character(0)
check <- function(ok, what) if (!isTRUE(ok)) fails <<- c(fails, what)
errs  <- function(expr) inherits(tryCatch({ force(expr); NULL }, error = function(e) e), "error")

# Synthetic Δ frame: 50 patients (23 / 27), an informative composite, a weak and a null
# clinical variable, two missing values in the first (median-imputed in-fold).
set.seed(11)
n  <- 50; yb <- c(rep(1L, 23), rep(0L, 27))
comp <- -0.9 * yb + rnorm(n)
Cmat <- cbind(PD_L1 = round(runif(n, 0, 90)), NLR = 4 - 0.8 * yb + rnorm(n, 0, 2))
Cmat[c(3, 40), "PD_L1"] <- NA
A <- list(comp = comp, Cmat = Cmat, yb = yb, y = factor(ifelse(yb == 1, "RP", "SD_PD"), levels = c("SD_PD", "RP")),
          neg_label = "SD_PD", resp_label = "RP")
seed <- 2026L

# --- 1. pairs / curves change nothing else; RNG untouched ------------------------------
set.seed(99); r_ref <- runif(3)
set.seed(99)
cv0 <- av_cv_kfold(A, reps = 20L, seed = seed + 977L, single = TRUE)
cv1 <- av_cv_kfold(A, reps = 20L, seed = seed + 977L, single = TRUE, pairs = TRUE, curves = TRUE)
check(identical(runif(3), r_ref), "caller's RNG state restored")
check(identical(cv0$R, cv1$R) && identical(cv0$Rs, cv1$Rs), "per-repetition AUCs of the existing arms unchanged")
check(identical(cv0[c("auc_clinical", "auc_immune", "auc_combined", "auc_single", "delta_auc", "delta_auc_vs_immune")],
                cv1[c("auc_clinical", "auc_immune", "auc_combined", "auc_single", "delta_auc", "delta_auc_vs_immune")]),
      "summary AUCs unchanged")
check(identical(names(cv1$auc_pair), c("PD_L1", "NLR")) && is.null(cv0$auc_pair), "one pair arm per clinical var, only on request")
check(is.null(av_cv_kfold(within(A, Cmat <- Cmat[, "NLR", drop = FALSE]), reps = 5L, pairs = TRUE)$auc_pair),
      "one clinical var: no pair arm (it is the combined model)")
# a pair arm is the 2-variable logistic model: with the composite alone informative, var + immune
# tracks the immune arm, not the var alone
check(abs(cv1$auc_pair$PD_L1 - cv1$auc_immune) < abs(cv1$auc_pair$PD_L1 - cv1$auc_single$PD_L1),
      "pair arm follows the immune composite")

# --- 2. curves: grid, columns, area = mean AUC -----------------------------------------
rc <- cv1$curves$roc; nb <- cv1$curves$net_benefit
arms <- c("clinical", "immune", "combined", "single_PD_L1", "single_NLR", "pair_PD_L1", "pair_NLR")
check(identical(names(rc), c("fpr", arms)) && identical(names(nb), c("threshold", arms)), "curve columns")
check(nrow(rc) == 501 && nrow(nb) == 60, "curve grids (FPR step 0.002, thresholds 0.01-0.60)")
trap <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
area <- vapply(arms, function(a) trap(rc$fpr, rc[[a]]), 0)
mauc <- c(cv1$auc_clinical, cv1$auc_immune, cv1$auc_combined, unlist(cv1$auc_single), unlist(cv1$auc_pair))
check(all(abs(area - mauc) < 0.005), sprintf("ROC area = mean CV AUC (max diff %.4f)", max(abs(area - mauc))))

# --- 3. av_attach_pair_curves() ---------------------------------------------------------
cvs <- av_cv_kfold(A, seed = seed + 977L, single = TRUE, curves = TRUE)     # the step-06 arms
old <- c("clinical", "immune", "combined", "single_PD_L1", "single_NLR")
av <- list(roc_cv = cvs$curves$roc[, c("fpr", old)],
           decision_curve = list(curve_cv = data.frame(cvs$curves$net_benefit[, c("threshold", old)],
                                                       treat_all = 0.4, treat_none = 0)),
           cv_kfold = cvs[c("auc_clinical", "auc_immune", "auc_combined", "auc_single")])
av2 <- suppressMessages(av_attach_pair_curves(av, A, seed))
check(all(c("pair_PD_L1", "pair_NLR") %in% names(av2$roc_cv)) &&
      all(c("pair_PD_L1", "pair_NLR") %in% names(av2$decision_curve$curve_cv)), "pair curves attached")
check(identical(av2$roc_cv[, c("fpr", old)], av$roc_cv) && identical(av2$cv_kfold[names(av$cv_kfold)], av$cv_kfold) &&
      identical(av2$decision_curve$curve_cv[, names(av$decision_curve$curve_cv)], av$decision_curve$curve_cv),
      "existing arms untouched")
check(identical(names(av2$cv_kfold$auc_pair), c("PD_L1", "NLR")), "auc_pair recorded")
bad <- av; bad$roc_cv$combined <- bad$roc_cv$combined + 1e-6
check(errs(suppressMessages(av_attach_pair_curves(bad, A, seed))), "tampered ROC curve: FATAL")
bad <- av; bad$cv_kfold$auc_single$NLR <- bad$cv_kfold$auc_single$NLR + 1e-6
check(errs(suppressMessages(av_attach_pair_curves(bad, A, seed))), "tampered single AUC: FATAL")
check(errs(suppressMessages(av_attach_pair_curves(av, A, seed + 1L))), "other seed (other folds): FATAL")
check(identical(av_attach_pair_curves(list(roc_cv = NULL), A, seed), list(roc_cv = NULL)), "no CV curves: unchanged")

# --- 4. bootstrap table rows ------------------------------------------------------------
pb <- list(B = 20L, seed = seed, results = list(gate = list(
  markers = c("CD28KI67", "KI67NAIVE"), n = 50L, n_pos = 23L,
  point = list(auc_clinical = .57, auc_immune = .75, auc_combined = .76, delta_auc = .19, delta_auc_vs_immune = .01,
               auc_single = list(PD_L1 = .34, NLR = .62), auc_pair = list(PD_L1 = .73, NLR = .79),
               delta_pair_vs_single = list(PD_L1 = .39, NLR = .16)),
  mc_spread = list(delta_auc = c(.1, .3)),
  bootA = list(n_valid = 20L, delta_auc = c(.04, .35), delta_auc_vs_immune = c(-.05, .13), p_two_sided = .02,
               auc = list(clinical = c(.39, .77), immune = c(.61, .87), combined = c(.64, .90)),
               auc_single = list(PD_L1 = c(.27, .62), NLR = c(.43, .78)),
               auc_pair = list(PD_L1 = c(.60, .87), NLR = c(.66, .90)),
               delta_pair_vs_single = list(PD_L1 = c(.12, .51), NLR = c(.03, .33))),
  delong = list(delta_auc = .2, delta_auc_ci = c(.01, .4)))))
tb <- patient_bootstrap_table(pb)
row <- function(q) tb[tb$quantity == q, ]
check(nrow(tb) == 7 + 2 + 4, sprintf("bootstrap table: 13 rows (%d)", nrow(tb)))
check(isTRUE(all.equal(unlist(row("AUC NLR + immune")[, c("point", "ci_low", "ci_high")]), c(point = .79, ci_low = .66, ci_high = .90))),
      "bootstrap table: AUC NLR + immune")
check(isTRUE(all.equal(unlist(row("dAUC PD_L1 + immune - PD_L1 alone")[, c("point", "ci_low", "ci_high")]),
                       c(point = .39, ci_low = .12, ci_high = .51))), "bootstrap table: dAUC pair - single")
pb$results$gate$point$auc_pair <- NULL
check(nrow(patient_bootstrap_table(pb)) == 9, "bootstrap table: no pair rows for a result without them")

# --- 5. figure: one ROC panel per clinical var --------------------------------------------
p <- pub_fig_pair_roc(av2)
check(inherits(p, "patchwork") && length(p) == 2, "figure: two panels")
lab1 <- levels(p[[1]]$data$M); lab2 <- levels(p[[2]]$data$M)
check(identical(lab1, c("TPS", "TPS + immune score")) && identical(lab2, c("NLR", "NLR + immune score")),
      "figure: TPS panel then NLR panel, display names")
b1 <- ggplot2::ggplot_build(p[[1]])$data[[2]]
check(setequal(unique(b1$colour), c("#B2182B", pub_palette[["combined"]])), "figure: TPS dark red, + immune score blue")
check(isTRUE(all.equal(p[[2]]$data$TPR[p[[2]]$data$M == "NLR + immune score"], av2$roc_cv$pair_NLR)),
      "figure: curves are the CV curves")
check(is.null(pub_fig_pair_roc(av)), "figure: NULL without pair arms")
work <- tempfile("pair_roc_"); dir.create(work)
f <- save_pub_figure(p, file.path(work, "FigurePairROC_test"), PUB_W2, 110)
check(!is.null(f) && file.exists(f) && file.size(f) > 0, "figure saved")
unlink(work, recursive = TRUE)

if (length(fails)) {
  cat(sprintf("[FAIL] %d check(s):\n", length(fails))); cat(sprintf("  - %s\n", fails), sep = "")
  quit(status = 1)
}
cat("[OK] clinical var + immune score arms: every check passed.\n")
