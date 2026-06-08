# multilevel_t0_utility.R
# ==============================================================================
# Publishability / clinical-utility diagnostic for BestResponse_2v3_4.
#
# GOAL (per user): the end product must be something measurable AT T0 that is
# clinically useful (pre-treatment prediction of RP vs SD_PD). This script tests
# whether richer longitudinal methods recover signal that the current Step 06
# nested-LOOCV (SVM AUC=0.56, perm p=0.20) misses, and — crucially — whether any
# of that signal is available from BASELINE alone.
#
# Three nested questions, increasing clinical value:
#   (D) Multilevel sPLS-DA   - does TREATMENT DYNAMICS separate groups? (needs T0+T1)
#   (C) Delta (T1-T0)        - does EARLY ON-TREATMENT CHANGE predict?  (needs T0+T1)
#   (B) T0-only per-marker   - which BASELINE markers carry signal?     (T0 only)  <-- the goal
#   (A) T0 composite score   - a single pre-treatment biomarker a clinician could use
#
# Note: for a 2-timepoint design, one-factor multilevel sPLS-DA discriminates on
# the within-patient direction == the Delta. So (D) and (C) should agree; we run
# both and cross-check (D = the named, citable method; C = the transparent one).
#
# Runs on the EXISTING discovery-73 longitudinal RDS (honest cohort the pipeline
# trains on). To probe the optimistic ceiling instead, point RDS_PATH at an
# original/full-cohort run. Writes JSON + PDFs to results/diag_multilevel_t0_<ts>/.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(here); library(jsonlite)
  library(mixOmics); library(pROC); library(glmnet); library(ggplot2)
})
options(crayon.enabled = FALSE)
set.seed(2026)

EXP   <- "BestResponse_2v3_4"
GATE  <- c("KI67NAIVE", "CD28KI67", "KI67CD4")   # LMM FDR-robust markers (Step 04)
N_PERM_FAST <- 2000   # for cheap LOO-AUC permutations
N_PERM_ML   <- 99     # for expensive multilevel BER permutations

RDS_PATH <- Sys.getenv("RDS_PATH", unset = here(
  "results/latest", EXP, "01_data_processing",
  sprintf("data_processed_%s_longitudinal.rds", EXP)))

stopifnot(file.exists(RDS_PATH))
d <- readRDS(RDS_PATH)
markers <- d$markers
meta    <- as.data.frame(d$metadata)
Z       <- d$hybrid_data_z[, markers, drop = FALSE]   # transformed + z-scored
Z       <- as.data.frame(lapply(Z, as.numeric))
# residual NAs (30 cells / 4 rows: markers fully missing in a few samples, not
# imputed by the longitudinal median step) -> median-impute for this diagnostic
Z[]     <- lapply(Z, function(x) { x[is.na(x)] <- median(x, na.rm = TRUE); x })
stopifnot(nrow(Z) == nrow(meta))
GATE <- intersect(GATE, markers)

out_dir <- here(paste0("results/diag_multilevel_t0_", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message(sprintf("[diag] RDS: %s", RDS_PATH))
message(sprintf("[diag] out: %s", out_dir))
message(sprintf("[diag] %d markers, %d rows (T0=%d, T1=%d), groups: %s",
                length(markers), nrow(meta), sum(meta$Timepoint == "T0"),
                sum(meta$Timepoint == "T1"),
                paste(names(table(meta$Group)), table(meta$Group), collapse = " ")))

# ---- helpers ----------------------------------------------------------------
auc_pos <- function(y, score) {  # AUC with RP as positive class; orientation-fixed
  as.numeric(pROC::auc(pROC::roc(y, score, levels = c("SD_PD", "RP"),
                                 direction = "<", quiet = TRUE)))
}
cohen_d <- function(x, g) {      # RP - SD_PD
  a <- x[g == "RP"]; b <- x[g == "SD_PD"]; a <- a[!is.na(a)]; b <- b[!is.na(b)]
  s <- sqrt(((length(a)-1)*var(a) + (length(b)-1)*var(b)) / (length(a)+length(b)-2))
  (mean(a) - mean(b)) / s
}
auc_orient <- function(y, score) {  # orientation-agnostic discrimination
  a <- auc_pos(y, score); max(a, 1 - a)
}
# LOO AUC for a ridge-logistic on matrix X (rows=patients) vs binary y
loo_auc_ridge <- function(X, y) {
  X <- as.matrix(X); n <- nrow(X); pred <- numeric(n)
  yb <- as.integer(y == "RP")
  for (i in seq_len(n)) {
    fit <- glmnet(X[-i, , drop = FALSE], yb[-i], family = "binomial",
                  alpha = 0, lambda = exp(seq(2, -4, length.out = 50)))
    # pick lambda by inner 5-fold CV on training rows
    cvf <- tryCatch(cv.glmnet(X[-i, , drop = FALSE], yb[-i], family = "binomial",
                              alpha = 0, nfolds = 5), error = function(e) NULL)
    lam <- if (!is.null(cvf)) cvf$lambda.min else median(fit$lambda)
    pred[i] <- predict(fit, X[i, , drop = FALSE], s = lam, type = "response")
  }
  auc_pos(y, pred)
}
perm_p_loo <- function(X, y, observed, n_perm, fun) {
  cnt <- 1
  for (p in seq_len(n_perm)) {
    yp <- sample(y)
    a  <- tryCatch(fun(X, yp), error = function(e) NA_real_)
    if (!is.na(a) && a >= observed) cnt <- cnt + 1
  }
  cnt / (n_perm + 1)
}

# =============================================================================
# (B) + (A)  T0-ONLY  — the clinically actionable layer
# =============================================================================
t0     <- meta$Timepoint == "T0"
Zt0    <- Z[t0, , drop = FALSE]
yt0    <- meta$Group[t0]
message(sprintf("\n[T0-only] n=%d  (RP=%d, SD_PD=%d)",
                length(yt0), sum(yt0 == "RP"), sum(yt0 == "SD_PD")))

# (B) per-marker baseline signal
t0_tab <- map_dfr(markers, function(m) {
  x <- Zt0[[m]]
  tibble(Marker = m,
         d_RPvsSDPD = cohen_d(x, yt0),
         auc_T0     = auc_pos(yt0, x),
         wilcox_p   = wilcox.test(x ~ yt0)$p.value)
}) %>% mutate(FDR = p.adjust(wilcox_p, "BH")) %>% arrange(wilcox_p)
write_csv(t0_tab, file.path(out_dir, "T0_per_marker.csv"))
message("[T0-only] top baseline markers (by Wilcoxon p):")
print(head(t0_tab, 8))

# (A) pre-specified T0 composite Ki67 score = mean z of GATE markers
comp_t0 <- rowMeans(Zt0[, GATE, drop = FALSE])
auc_comp <- auc_pos(yt0, comp_t0)
roc_comp <- pROC::roc(yt0, comp_t0, levels = c("SD_PD","RP"), direction = "<", quiet = TRUE)
ci_comp  <- as.numeric(pROC::ci.auc(roc_comp, method = "bootstrap", boot.n = 2000))
yj       <- pROC::coords(roc_comp, "best", best.method = "youden",
                         ret = c("threshold","sensitivity","specificity"))
# permutation p for the fixed composite (shuffle labels)
perm_comp <- mean(c(TRUE, replicate(N_PERM_FAST,
                auc_pos(sample(yt0), comp_t0) >= auc_comp)))
message(sprintf("[T0-only] composite Ki67 score: AUC=%.3f [%.3f-%.3f], perm p=%.3f",
                auc_comp, ci_comp[1], ci_comp[3], perm_comp))
message(sprintf("[T0-only] Youden cutpoint=%.3f -> sens=%.2f spec=%.2f",
                yj[[1]], yj[[2]], yj[[3]]))

# T0 composite as an honest LOO classifier (single feature, ridge-free: score itself)
auc_comp_loo <- auc_comp  # fixed pre-specified score => LOO == apparent for a 1-D score

# =============================================================================
# (C)  DELTA (T1 - T0)  — early on-treatment dynamics
# =============================================================================
paired_ids <- meta %>% count(Patient_ID) %>% filter(n == 2) %>% pull(Patient_ID)
Zd <- map_dfr(paired_ids, function(pid) {
  r0 <- which(meta$Patient_ID == pid & meta$Timepoint == "T0")
  r1 <- which(meta$Patient_ID == pid & meta$Timepoint == "T1")
  as_tibble(Z[r1, , drop = FALSE] - Z[r0, , drop = FALSE])
})
yd <- meta$Group[match(paired_ids, meta$Patient_ID)]
message(sprintf("\n[Delta] n_paired=%d  (RP=%d, SD_PD=%d)",
                length(yd), sum(yd == "RP"), sum(yd == "SD_PD")))

# per-marker Delta signal
d_tab <- map_dfr(markers, function(m) {
  x <- Zd[[m]]
  tibble(Marker = m, d_delta = cohen_d(x, yd),
         auc_delta = auc_pos(yd, x), wilcox_p = wilcox.test(x ~ yd)$p.value)
}) %>% mutate(FDR = p.adjust(wilcox_p, "BH")) %>% arrange(wilcox_p)
write_csv(d_tab, file.path(out_dir, "Delta_per_marker.csv"))
message("[Delta] top markers (by Wilcoxon p):"); print(head(d_tab, 8))

# composite Delta Ki67 score (orientation-agnostic: RP contracts Ki67 MORE, so
# the raw RP-high AUC is <0.5; report discrimination magnitude + direction)
comp_d   <- rowMeans(Zd[, GATE, drop = FALSE])
auc_cd_r <- auc_pos(yd, comp_d)
auc_cd   <- max(auc_cd_r, 1 - auc_cd_r)
dir_cd   <- ifelse(auc_cd_r < 0.5, "RP_decreases", "RP_increases")
perm_cd  <- mean(c(TRUE, replicate(N_PERM_FAST, auc_orient(sample(yd), comp_d) >= auc_cd)))
message(sprintf("[Delta] composite Ki67 Delta score: AUC=%.3f (%s), perm p=%.3f",
                auc_cd, dir_cd, perm_cd))

# multivariate Delta classifier (ridge LOO over GATE markers)
auc_d_ml  <- loo_auc_ridge(Zd[, GATE, drop = FALSE], yd)
perm_d_ml <- perm_p_loo(Zd[, GATE, drop = FALSE], yd, auc_d_ml, N_PERM_ML, loo_auc_ridge)
message(sprintf("[Delta] ridge LOO (gate markers): AUC=%.3f, perm p=%.3f", auc_d_ml, perm_d_ml))

# =============================================================================
# (D)  MULTILEVEL sPLS-DA  (mixOmics, named/citable method)
# =============================================================================
# NOTE: mixOmics' native `multilevel=` argument is buggy in this version
# ('x must be an array of at least two dimensions'). We reproduce the method
# explicitly: withinVariation() decomposes the paired matrix into the
# within-patient component, then a sparse sPLS-DA discriminates on it. This is
# exactly what multilevel sPLS-DA does internally. Because each patient's two
# within-rows are mirror images, the patient-level signal == the Delta direction,
# so the LEAKAGE-FREE significance is the Delta ridge LOO above (perm p=0.010);
# here we report the citable variable selection + perf AUC point estimate.
multilevel_res <- tryCatch({
  pr_idx <- meta$Patient_ID %in% paired_ids
  Xp <- as.matrix(Z[pr_idx, , drop = FALSE]); Yp <- factor(meta$Group[pr_idx])
  design <- data.frame(sample = as.numeric(factor(meta$Patient_ID[pr_idx])))
  Xw  <- withinVariation(X = Xp, design = design)   # within-patient component
  kx  <- 5
  fit <- splsda(Xw, Yp, ncomp = 2, keepX = rep(kx, 2))
  sel <- selectVar(fit, comp = 1)$name
  pf  <- suppressWarnings(perf(fit, validation = "Mfold", folds = 5,
                               nrepeat = 50, auc = TRUE, progressBar = FALSE))
  auc_raw <- as.numeric(pf$auc$comp1["AUC.mean"])
  list(method = "withinVariation + sparse sPLS-DA (multilevel equivalent)",
       keepX = kx, auc_oriented = max(auc_raw, 1 - auc_raw),
       ber = as.numeric(pf$error.rate$BER[1, "max.dist"]),
       selected = head(sel, 10),
       significance_note = sprintf("Patient-level significance = Delta ridge LOO perm p=%.3f", perm_d_ml))
}, error = function(e) { message("[multilevel] FAILED: ", e$message); NULL })

if (!is.null(multilevel_res)) {
  message(sprintf("\n[multilevel] withinVariation+sPLS-DA keepX=%d  AUC(oriented)=%.3f  selected: %s",
                  multilevel_res$keepX, multilevel_res$auc_oriented,
                  paste(multilevel_res$selected, collapse = ", ")))
}

# =============================================================================
# BRIDGE: which dynamics-selected markers ALSO have T0 baseline signal?
# =============================================================================
tryCatch({
  bridge <- t0_tab %>%
    select(Marker, auc_T0, d_T0 = d_RPvsSDPD, p_T0 = wilcox_p) %>%
    left_join(d_tab %>% select(Marker, auc_delta, d_delta, p_delta = wilcox_p), by = "Marker") %>%
    mutate(in_gate = Marker %in% GATE) %>%
    arrange(p_T0)
  write_csv(bridge, file.path(out_dir, "bridge_T0_vs_delta.csv"))
}, error = function(e) message("[bridge] ", e$message))

# =============================================================================
# PLOTS
# =============================================================================
tryCatch({
  topm <- head(t0_tab, 15) %>% mutate(Marker = fct_reorder(Marker, d_RPvsSDPD))
  p1 <- ggplot(topm, aes(d_RPvsSDPD, Marker, fill = d_RPvsSDPD > 0)) +
    geom_col() + geom_vline(xintercept = 0) +
    scale_fill_manual(values = c(`TRUE` = "#2E8B57", `FALSE` = "#B2182B"), guide = "none") +
    labs(title = "T0 baseline effect size (Cohen's d, RP - SD_PD)",
         subtitle = sprintf("%s discovery cohort, n=%d", EXP, length(yt0)),
         x = "Cohen's d", y = NULL) + theme_bw()
  ggsave(file.path(out_dir, "T0_effect_sizes.pdf"), p1, width = 7, height = 6)

  pdf(file.path(out_dir, "T0_composite_ROC.pdf"), width = 6, height = 6)
  plot(roc_comp, main = sprintf("T0 composite Ki67 score — AUC=%.3f [%.3f-%.3f]",
                                auc_comp, ci_comp[1], ci_comp[3]), col = "#2166AC", lwd = 2)
  dev.off()
}, error = function(e) message("[plot] ", e$message))

# =============================================================================
# SUMMARY JSON + verdict
# =============================================================================
summary <- list(
  rds = RDS_PATH, n_T0 = length(yt0), n_paired = length(yd), gate = GATE,
  T0_only = list(
    composite_auc = auc_comp, composite_ci = ci_comp[c(1,3)], composite_perm_p = perm_comp,
    youden = list(threshold = yj[[1]], sens = yj[[2]], spec = yj[[3]]),
    best_marker = t0_tab$Marker[1], best_marker_auc = t0_tab$auc_T0[1],
    best_marker_p = t0_tab$wilcox_p[1], n_FDR_sig = sum(t0_tab$FDR < 0.05)),
  delta = list(composite_auc = auc_cd, composite_direction = dir_cd,
               composite_perm_p = perm_cd,
               ridge_loo_auc = auc_d_ml, ridge_loo_perm_p = perm_d_ml,
               n_FDR_sig = sum(d_tab$FDR < 0.05)),
  multilevel = multilevel_res,
  step06_reference = list(svm_auc = 0.56, perm_p = 0.198)
)
write_json(summary, file.path(out_dir, "summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = 4)

message("\n=================== VERDICT ===================")
message(sprintf("T0 composite (clinical target): AUC=%.3f, perm p=%.3f  [%s]",
                auc_comp, perm_comp, ifelse(perm_comp < 0.05, "SIGNIFICANT", "n.s.")))
message(sprintf("T0 best single marker         : %s AUC=%.3f, FDR-sig markers=%d",
                t0_tab$Marker[1], t0_tab$auc_T0[1], sum(t0_tab$FDR < 0.05)))
message(sprintf("Delta composite               : AUC=%.3f, perm p=%.3f", auc_cd, perm_cd))
message(sprintf("Delta ridge LOO               : AUC=%.3f, perm p=%.3f", auc_d_ml, perm_d_ml))
if (!is.null(multilevel_res))
  message(sprintf("Multilevel sPLS-DA            : AUC(oriented)=%.3f (sig via Delta LOO)",
                  multilevel_res$auc_oriented))
message(sprintf("Step 06 reference (T0 SVM)     : AUC=0.56, perm p=0.20"))
message(sprintf("\nOutputs -> %s", out_dir))
