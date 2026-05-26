# =============================================================================
# diag_14_t0_delta_signal_decomposition.R
#
# PURPOSE: Decompose the Step 06 ML signal into T0 baseline vs longitudinal
#          delta components. Determines the correct scientific narrative and
#          checks for feature-selection leakage via the LMM gate.
#
# REVIEWER CRITIQUE CHECKS:
#   RC1 – "LMM selected features using T1; T0 AUC may be inflated by selection bias"
#   RC2 – "Is the signal baseline-prognostic or dynamic-predictive?"
#   RC3 – "What AUC would you get without LMM guidance — just top T0 markers?"
#   RC4 – "Does baseline KI67 level predict the magnitude of the reset?"
#   RC5 – "Is the Timepoint×Group interaction driven by T0 differences, T1
#          differences, or the change itself?"
#   RC6 – "Are delta-based features sufficient to replace the LMM pipeline?"
#
# TESTS:
#   1.  Descriptive: T0, T1, delta medians ± IQR by group (RP / SD+PD)
#   2.  Wilcoxon + AUC: T0, T1, delta for KI67NAIVE / KI67CD4 / CD28KI67
#   3.  DeLong pairwise: AUC(T0) vs AUC(delta) vs AUC(T1) per marker
#   4.  Spaghetti plot: individual T0→T1 trajectories coloured by group
#   5.  Correlation: T0 level vs delta magnitude (reset potential hypothesis)
#   6.  Pure T0 Wilcoxon scan: which markers would be selected without LMM?
#        Do KI67NAIVE / KI67CD4 appear in top-k?
#   7.  LOO SVM on delta features (parallel to existing T0 SVM)
#   8.  LOO SVM on T0-only univariate gate (no LMM) — leakage upper bound
#   9.  LOO SVM on T1 features — upper bound of post-treatment signal
#   10. Leakage sensitivity: compare AUC(LMM-gated T0) vs AUC(T0-gate naive)
#   11. Summary verdict table
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(pROC)
  library(ggplot2)
  library(gridExtra)
  library(e1071)
  library(readxl)
})

# Source pipeline modules to use the exact same SVM implementation
suppressMessages({
  source(here::here("R/modules_ml.R"))
})

here_root <- function(...) {
  root <- tryCatch(here::here(), error = function(e) getwd())
  file.path(root, ...)
}

set.seed(2026)
N_PERM   <- 999
GATE_MARKERS <- c("KI67NAIVE", "KI67CD4")   # Step 06 LMM gate
ALL_KI67     <- c("KI67NAIVE", "KI67CD4", "CD28KI67")  # all LMM-sig

cat("==========================================================\n")
cat("diag_14 — T0 / Delta / T1 Signal Decomposition\n")
cat("==========================================================\n\n")

# ── 1. LOAD DATA ─────────────────────────────────────────────────────────────

rds_std  <- readRDS(here_root("results/BestResponse_2v3_4/01_data_processing",
                               "data_processed_BestResponse_2v3_4_standard.rds"))
rds_long <- readRDS(here_root("results/BestResponse_2v3_4/01_data_processing",
                               "data_processed_BestResponse_2v3_4_longitudinal.rds"))

meta_std  <- rds_std$metadata
meta_long <- rds_long$metadata

z_std  <- as.data.frame(rds_std$hybrid_data_z)
z_long <- as.data.frame(rds_long$hybrid_data_z)

# Patient_ID is a column, Sample_ID encodes Patient_ID_T0 / _T1
z_long$Patient_ID <- meta_long$Patient_ID
z_long$Timepoint  <- meta_long$Timepoint
z_long$Group      <- meta_long$Group

z_std$Patient_ID  <- meta_std$Patient_ID
z_std$Group       <- meta_std$Group

cat("Standard RDS: n =", nrow(z_std), " patients\n")
cat("Longitudinal RDS: n_rows =", nrow(z_long),
    " (T0:", sum(z_long$Timepoint == "T0"),
    " T1:", sum(z_long$Timepoint == "T1"), ")\n\n")

# Build paired dataset (patients with both T0 and T1)
z_t0 <- z_long %>% filter(Timepoint == "T0")
z_t1 <- z_long %>% filter(Timepoint == "T1")

paired_ids <- intersect(z_t0$Patient_ID, z_t1$Patient_ID)
cat("Paired patients (T0+T1):", length(paired_ids), "\n\n")

z_t0p <- z_t0 %>% filter(Patient_ID %in% paired_ids) %>%
  arrange(Patient_ID)
z_t1p <- z_t1 %>% filter(Patient_ID %in% paired_ids) %>%
  arrange(Patient_ID)

# Delta = T1 - T0 (on paired patients)
markers_avail <- intersect(c(GATE_MARKERS, ALL_KI67),
                            intersect(names(z_t0p), names(z_t1p)))

delta_df <- z_t0p[, c("Patient_ID", "Group")]
for (m in markers_avail) {
  delta_df[[m]] <- z_t1p[[m]] - z_t0p[[m]]
}

# ── 2. HELPER FUNCTIONS ───────────────────────────────────────────────────────

wilcox_auc <- function(x, y_bin) {
  # y_bin: 1=RP, 0=NR; x = continuous marker
  complete <- !is.na(x) & !is.na(y_bin)
  x <- x[complete]; y_bin <- y_bin[complete]
  wt  <- wilcox.test(x[y_bin == 1], x[y_bin == 0], exact = FALSE)
  roc <- pROC::roc(y_bin, x, quiet = TRUE)
  auc_val <- as.numeric(pROC::auc(roc))
  ci_val  <- as.numeric(pROC::ci.auc(roc, method = "delong"))
  list(p = wt$p.value, auc = auc_val,
       ci_lo = ci_val[1], ci_hi = ci_val[3],
       n = sum(complete))
}

cohen_d <- function(x, y_bin) {
  g1 <- x[y_bin == 1]; g0 <- x[y_bin == 0]
  sp <- sqrt(((length(g1)-1)*var(g1) + (length(g0)-1)*var(g0)) /
               (length(g1)+length(g0)-2))
  (mean(g1) - mean(g0)) / sp
}

perm_auc <- function(x, y_bin, n_perm = N_PERM, seed = 42) {
  set.seed(seed)
  obs <- as.numeric(pROC::auc(pROC::roc(y_bin, x, quiet = TRUE)))
  null <- replicate(n_perm, {
    as.numeric(pROC::auc(pROC::roc(sample(y_bin), x, quiet = TRUE)))
  })
  mean(null >= obs)
}

svm_loo_auc_from_res <- function(res) {
  probs <- res$predicted_probs
  y_bin <- as.integer(res$y_true == res$positive_label)
  as.numeric(pROC::auc(pROC::roc(y_bin, probs, quiet = TRUE, direction = "<")))
}

run_svm_loo <- function(X, y_factor, label = "") {
  cc <- complete.cases(X)
  X  <- X[cc, , drop = FALSE]
  y  <- y_factor[cc]
  message(sprintf("   [diag_14] %s: n=%d", label, nrow(X)))
  run_nested_loocv_svm(X, y)
}

# ── 3. DESCRIPTIVE STATISTICS ─────────────────────────────────────────────────

cat("── 3. Descriptive: T0 / T1 / Delta medians by group ──────────────────\n")

desc_rows <- list()
for (m in markers_avail) {
  y_t0 <- as.integer(z_t0p$Group == "RP")
  y_t1 <- as.integer(z_t1p$Group == "RP")
  y_d  <- as.integer(delta_df$Group == "RP")

  for (tp in c("T0","T1","Delta")) {
    x <- if (tp == "T0") z_t0p[[m]] else
         if (tp == "T1") z_t1p[[m]] else delta_df[[m]]
    y <- if (tp == "Delta") y_d else y_t0
    rp  <- x[y == 1]; nr  <- x[y == 0]
    desc_rows[[length(desc_rows)+1]] <- data.frame(
      Marker = m, Timepoint = tp,
      Median_RP  = round(median(rp, na.rm=TRUE), 3),
      Median_NR  = round(median(nr, na.rm=TRUE), 3),
      Diff       = round(median(rp, na.rm=TRUE) - median(nr, na.rm=TRUE), 3)
    )
  }
}
desc_tbl <- do.call(rbind, desc_rows)
print(desc_tbl, row.names = FALSE)
cat("\n")

# ── 4. AUC / WILCOXON DECOMPOSITION ──────────────────────────────────────────

cat("── 4. AUC / Wilcoxon decomposition: T0 vs T1 vs Delta ────────────────\n")

auc_rows <- list()
for (m in markers_avail) {
  y_paired <- as.integer(z_t0p$Group == "RP")

  r_t0 <- wilcox_auc(z_t0p[[m]], y_paired)
  r_t1 <- wilcox_auc(z_t1p[[m]], y_paired)
  r_d  <- wilcox_auc(delta_df[[m]], y_paired)
  d_t0 <- cohen_d(z_t0p[[m]], y_paired)
  d_d  <- cohen_d(delta_df[[m]], y_paired)

  for (tp in c("T0","T1","Delta")) {
    r <- if (tp=="T0") r_t0 else if (tp=="T1") r_t1 else r_d
    d <- if (tp=="T0") d_t0 else if (tp=="Delta") d_d else NA
    auc_rows[[length(auc_rows)+1]] <- data.frame(
      Marker = m, Timepoint = tp,
      AUC = round(r$auc, 3),
      CI_lo = round(r$ci_lo, 3),
      CI_hi = round(r$ci_hi, 3),
      Wilcox_p = round(r$p, 4),
      Cohen_d = round(d, 3),
      n = r$n
    )
  }
}
auc_tbl <- do.call(rbind, auc_rows)
print(auc_tbl, row.names = FALSE)
cat("\n")

# ── 5. DELONG: T0 vs DELTA vs T1 per marker ──────────────────────────────────

cat("── 5. DeLong pairwise: AUC(T0) vs AUC(Delta) vs AUC(T1) ─────────────\n")

for (m in markers_avail) {
  y <- as.integer(z_t0p$Group == "RP")
  roc_t0    <- pROC::roc(y, z_t0p[[m]], quiet = TRUE, direction = "auto")
  roc_t1    <- pROC::roc(y, z_t1p[[m]], quiet = TRUE, direction = "auto")
  roc_delta <- pROC::roc(y, delta_df[[m]], quiet = TRUE, direction = "auto")

  p_t0_d  <- pROC::roc.test(roc_t0, roc_delta, method = "delong")$p.value
  p_t0_t1 <- pROC::roc.test(roc_t0, roc_t1,    method = "delong")$p.value
  p_d_t1  <- pROC::roc.test(roc_delta, roc_t1,  method = "delong")$p.value

  cat(sprintf("  %s: T0 vs Delta p=%.3f | T0 vs T1 p=%.3f | Delta vs T1 p=%.3f\n",
              m, p_t0_d, p_t0_t1, p_d_t1))
}
cat("\n")

# ── 6. CORRELATION: T0 LEVEL vs DELTA MAGNITUDE ──────────────────────────────

cat("── 6. Correlation: T0 level vs delta magnitude (reset potential) ──────\n")
cat("   Hypothesis: high T0 KI67 → bigger drop in RP (reset potential)\n\n")

for (m in markers_avail) {
  y <- as.integer(z_t0p$Group == "RP")
  # All patients: T0 vs delta
  ct_all <- cor.test(z_t0p[[m]], delta_df[[m]], method = "spearman")
  # RP only
  ct_rp  <- cor.test(z_t0p[[m]][y==1], delta_df[[m]][y==1], method = "spearman")
  # NR only
  ct_nr  <- cor.test(z_t0p[[m]][y==0], delta_df[[m]][y==0], method = "spearman")

  cat(sprintf("  %s: rho(all)=%.3f p=%.3f | rho(RP)=%.3f p=%.3f | rho(NR)=%.3f p=%.3f\n",
              m, ct_all$estimate, ct_all$p.value,
              ct_rp$estimate, ct_rp$p.value,
              ct_nr$estimate, ct_nr$p.value))
}
cat("\n")

# ── 7. PURE T0 WILCOXON SCAN (no LMM guidance) ───────────────────────────────

cat("── 7. Pure T0 Wilcoxon scan — which markers survive without LMM? ──────\n")
cat("   RC1/RC3 check: would KI67NAIVE/KI67CD4 be selected from T0 alone?\n\n")

all_markers <- intersect(rds_std$hybrid_markers, names(z_std))
all_markers <- setdiff(all_markers, c("Patient_ID","Group","Timepoint","Sample_ID"))

y_std <- as.integer(z_std$Group == "RP")

t0_scan <- lapply(all_markers, function(m) {
  x <- z_std[[m]]
  if (all(is.na(x)) || length(unique(x[!is.na(x)])) < 2) return(NULL)
  wt  <- wilcox.test(x[y_std==1], x[y_std==0], exact = FALSE)
  roc <- pROC::roc(y_std, x, quiet = TRUE)
  data.frame(Marker = m, p_t0 = wt$p.value,
             AUC_t0 = as.numeric(pROC::auc(roc)))
})
t0_scan <- do.call(rbind, Filter(Negate(is.null), t0_scan))
t0_scan$FDR_t0 <- p.adjust(t0_scan$p_t0, method = "BH")
t0_scan <- t0_scan[order(t0_scan$p_t0), ]

cat("  Top 10 T0 markers by Wilcoxon p (all", nrow(t0_scan),"markers):\n")
print(head(t0_scan, 10), row.names = FALSE)

gate_rank_naive <- which(t0_scan$Marker == "KI67NAIVE")
gate_rank_cd4   <- which(t0_scan$Marker == "KI67CD4")
cat(sprintf("\n  KI67NAIVE rank in T0 scan: %d / %d (p=%.4f, FDR=%.3f, AUC=%.3f)\n",
            gate_rank_naive, nrow(t0_scan),
            t0_scan$p_t0[gate_rank_naive],
            t0_scan$FDR_t0[gate_rank_naive],
            t0_scan$AUC_t0[gate_rank_naive]))
cat(sprintf("  KI67CD4 rank in T0 scan:   %d / %d (p=%.4f, FDR=%.3f, AUC=%.3f)\n\n",
            gate_rank_cd4, nrow(t0_scan),
            t0_scan$p_t0[gate_rank_cd4],
            t0_scan$FDR_t0[gate_rank_cd4],
            t0_scan$AUC_t0[gate_rank_cd4]))

# ── 8. LOO SVM: T0 FEATURES (current pipeline — exact replication) ───────────

cat("── 8. LOO SVM AUC: T0 features (LMM gate, nested inner CV) ───────────\n")
cat("   Uses run_nested_loocv_svm() — identical to Step 06 pipeline\n\n")

X_t0  <- as.matrix(z_std[, GATE_MARKERS])
y_sv  <- factor(ifelse(z_std$Group == "RP", "RP", "NR"), levels = c("NR","RP"))
res_t0 <- run_svm_loo(X_t0, y_sv, "T0 LMM gate")
auc_t0_svm <- svm_loo_auc_from_res(res_t0)
cat(sprintf("  LOO SVM AUC (T0, n=%d): %.3f\n\n", nrow(X_t0), auc_t0_svm))

# ── 9. LOO SVM: DELTA FEATURES ────────────────────────────────────────────────

cat("── 9. LOO SVM AUC: Delta features (paired patients only) ─────────────\n")

X_delta <- as.matrix(delta_df[, GATE_MARKERS])
y_delta <- factor(ifelse(delta_df$Group == "RP", "RP", "NR"), levels = c("NR","RP"))
res_delta <- run_svm_loo(X_delta, y_delta, "Delta LMM gate")
auc_delta_svm <- svm_loo_auc_from_res(res_delta)
cat(sprintf("  LOO SVM AUC (Delta, n=%d): %.3f\n\n", nrow(X_delta), auc_delta_svm))

# ── 10. LOO SVM: T1 FEATURES ─────────────────────────────────────────────────

cat("── 10. LOO SVM AUC: T1 features (paired patients only) ───────────────\n")

X_t1v <- as.matrix(z_t1p[, GATE_MARKERS])
y_t1v <- factor(ifelse(z_t1p$Group == "RP", "RP", "NR"), levels = c("NR","RP"))
res_t1 <- run_svm_loo(X_t1v, y_t1v, "T1 LMM gate")
auc_t1_svm <- svm_loo_auc_from_res(res_t1)
cat(sprintf("  LOO SVM AUC (T1, n=%d): %.3f\n\n", nrow(X_t1v), auc_t1_svm))

# ── 11. LOO SVM: PURE T0 NAIVE GATE (top-2 T0 markers, no LMM) ──────────────

cat("── 11. LOO SVM: T0-naive gate (top-2 by Wilcoxon T0, no LMM) ─────────\n")
cat("   RC3 leakage check: what AUC without LMM feature selection?\n\n")

top2_t0 <- as.character(head(t0_scan$Marker, 2))
cat("  Top-2 T0 markers selected:", paste(top2_t0, collapse=", "), "\n")

X_naive <- as.matrix(z_std[, top2_t0, drop = FALSE])
res_naive <- run_svm_loo(X_naive, y_sv, "T0 naive top-2")
auc_naive_svm <- svm_loo_auc_from_res(res_naive)
cat(sprintf("  LOO SVM AUC (T0 naive top-2, n=%d): %.3f\n\n",
            nrow(X_naive), auc_naive_svm))

# Top-3 with KI67NAIVE included
top3_ki67 <- c("KI67NAIVE", "KI67CD4", "CD28KI67")
X_ki67_3 <- as.matrix(z_std[, top3_ki67, drop = FALSE])
res_ki67_3 <- run_svm_loo(X_ki67_3, y_sv, "T0 all-KI67 top-3")
auc_ki67_3_svm <- svm_loo_auc_from_res(res_ki67_3)
cat(sprintf("  LOO SVM AUC (T0 KI67 trio, n=%d): %.3f\n\n",
            nrow(X_ki67_3), auc_ki67_3_svm))

# ── 12. PERMUTATION TEST FOR DELTA SVM ───────────────────────────────────────

cat("── 12. Permutation p-value: Delta SVM ────────────────────────────────\n")
cat("   (uses 199 perms for speed — increase N_PERM for final report)\n\n")

N_PERM_FAST <- 199L
set.seed(2026)
auc_delta_null <- replicate(N_PERM_FAST, {
  y_perm <- factor(sample(as.character(y_delta)), levels = levels(y_delta))
  tryCatch(
    svm_loo_auc_from_res(run_nested_loocv_svm(X_delta, y_perm)),
    error = function(e) NA_real_
  )
})
auc_delta_null <- auc_delta_null[!is.na(auc_delta_null)]
perm_p_delta <- mean(auc_delta_null >= auc_delta_svm)
cat(sprintf("  Delta SVM AUC=%.3f, perm p=%.3f (n_perm=%d)\n\n",
            auc_delta_svm, perm_p_delta, length(auc_delta_null)))

# ── 13. SPAGHETTI PLOTS ───────────────────────────────────────────────────────

cat("── 13. Generating spaghetti plots (T0→T1 trajectories) ───────────────\n")

spaghetti_data <- bind_rows(
  z_t0p %>% select(Patient_ID, Group, all_of(markers_avail)) %>%
    mutate(Timepoint = "T0"),
  z_t1p %>% select(Patient_ID, Group, all_of(markers_avail)) %>%
    mutate(Timepoint = "T1")
) %>% pivot_longer(cols = all_of(markers_avail),
                   names_to = "Marker", values_to = "Z_score")

spaghetti_data$Timepoint <- factor(spaghetti_data$Timepoint, levels = c("T0","T1"))
spaghetti_data$Group_col <- ifelse(spaghetti_data$Group == "RP",
                                    "#2166ac", "#d6604d")

p_spag <- ggplot(spaghetti_data,
                 aes(x = Timepoint, y = Z_score,
                     group = Patient_ID, colour = Group)) +
  geom_line(alpha = 0.4, linewidth = 0.5) +
  geom_point(size = 1.2, alpha = 0.7) +
  stat_summary(aes(group = Group), fun = mean, geom = "line",
               linewidth = 1.5, linetype = "dashed") +
  scale_colour_manual(values = c("RP" = "#2166ac", "SD_PD" = "#d6604d")) +
  facet_wrap(~Marker, scales = "free_y") +
  labs(title = "diag_14: T0→T1 trajectories by group",
       subtitle = "Dashed = group mean; blue=RP, red=SD/PD",
       y = "Z-score", x = NULL) +
  theme_bw(base_size = 11)

ggsave(here_root("diagnostics/diag_14_spaghetti.pdf"),
       p_spag, width = 10, height = 4)
cat("  Saved: diagnostics/diag_14_spaghetti.pdf\n")

# Boxplot: T0, T1, delta per marker per group
box_data <- bind_rows(
  z_t0p %>% select(Patient_ID, Group, all_of(markers_avail)) %>%
    mutate(Timepoint = "T0"),
  z_t1p %>% select(Patient_ID, Group, all_of(markers_avail)) %>%
    mutate(Timepoint = "T1"),
  delta_df %>% select(Patient_ID, Group, all_of(markers_avail)) %>%
    mutate(Timepoint = "Delta (T1-T0)")
) %>% pivot_longer(cols = all_of(markers_avail),
                   names_to = "Marker", values_to = "Z_score")

box_data$Timepoint <- factor(box_data$Timepoint,
                              levels = c("T0","T1","Delta (T1-T0)"))

p_box <- ggplot(box_data,
                aes(x = Timepoint, y = Z_score, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  scale_fill_manual(values = c("RP" = "#2166ac", "SD_PD" = "#d6604d")) +
  facet_wrap(~Marker, scales = "free_y") +
  labs(title = "diag_14: T0 / T1 / Delta distributions by group",
       y = "Z-score", x = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(here_root("diagnostics/diag_14_boxplots.pdf"),
       p_box, width = 10, height = 4)
cat("  Saved: diagnostics/diag_14_boxplots.pdf\n\n")

# Scatter: T0 vs delta coloured by group
scatter_plots <- lapply(markers_avail, function(m) {
  df_sc <- data.frame(T0 = z_t0p[[m]], Delta = delta_df[[m]],
                      Group = z_t0p$Group)
  ggplot(df_sc, aes(T0, Delta, colour = Group)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, formula = y~x, linewidth = 0.8) +
    scale_colour_manual(values = c("RP"="#2166ac","SD_PD"="#d6604d")) +
    labs(title = m, x = "T0 z-score", y = "Delta (T1-T0) z-score") +
    theme_bw(base_size = 10)
})
p_scatter <- do.call(gridExtra::grid.arrange,
                     c(scatter_plots, list(ncol = length(markers_avail))))
ggsave(here_root("diagnostics/diag_14_scatter_t0_vs_delta.pdf"),
       p_scatter, width = 4 * length(markers_avail), height = 4)
cat("  Saved: diagnostics/diag_14_scatter_t0_vs_delta.pdf\n\n")

# ── 14. SUMMARY VERDICT ───────────────────────────────────────────────────────

cat("==========================================================\n")
cat("SUMMARY VERDICT\n")
cat("==========================================================\n\n")

cat(sprintf("  LOO SVM AUC — T0    (LMM gate, n=%d): %.3f\n",
            nrow(X_t0), auc_t0_svm))
cat(sprintf("  LOO SVM AUC — T1    (LMM gate, n=%d): %.3f\n",
            nrow(X_t1v), auc_t1_svm))
cat(sprintf("  LOO SVM AUC — Delta (LMM gate, n=%d): %.3f  [perm p=%.3f]\n",
            nrow(X_delta), auc_delta_svm, perm_p_delta))
cat(sprintf("  LOO SVM AUC — T0 naive top-2:          %.3f\n",
            auc_naive_svm))
cat(sprintf("  LOO SVM AUC — T0 KI67 trio (3 feat):   %.3f\n\n",
            auc_ki67_3_svm))

auc_t0_num    <- auc_t0_svm
auc_delta_num <- auc_delta_svm
auc_t1_num    <- auc_t1_svm
auc_naive_num <- auc_naive_svm

cat("  Narrative determination:\n")

if (abs(auc_t0_num - auc_delta_num) < 0.05) {
  cat("  SCENARIO A/MIXED: T0 and Delta give similar AUC.\n")
  cat("  -> T0 baseline already captures the discriminative signal.\n")
  cat("  -> Narrative: LMM identifies biological mechanism (reset); T0 prediction\n")
  cat("     sufficient for clinical application. Both findings complement each other.\n")
} else if (auc_delta_num > auc_t0_num + 0.05) {
  cat("  SCENARIO B: Delta AUC substantially exceeds T0 AUC.\n")
  cat("  -> Dynamic change is the primary signal; T0 is weaker.\n")
  cat("  -> Narrative: KI67 reset IS the finding; T1 required for full signal.\n")
  cat("  -> Step 06 should use delta features for the primary ML result.\n")
} else {
  cat("  SCENARIO C: T0 AUC >= Delta AUC.\n")
  cat("  -> Baseline KI67 is the primary discriminator.\n")
  cat("  -> LMM confirms dynamic biology but T0 alone sufficient for prediction.\n")
}

cat("\n  Feature-selection leakage check:\n")
lmm_lift <- auc_t0_num - auc_naive_num
if (abs(lmm_lift) < 0.05) {
  cat(sprintf("  -> LMM gate lift over naive T0 gate: ΔAUC=%.3f (marginal)\n", lmm_lift))
  cat("  -> LMM leakage unlikely to explain T0 AUC — signal appears genuine at T0.\n")
} else if (lmm_lift > 0.05) {
  cat(sprintf("  -> LMM gate lift over naive T0 gate: ΔAUC=+%.3f\n", lmm_lift))
  cat("  -> LMM feature selection contributes to T0 AUC. Acknowledge in paper;\n")
  cat("     delta-based SVM is the cleanest validation.\n")
} else {
  cat(sprintf("  -> Naive T0 gate outperforms LMM gate by ΔAUC=%.3f\n", -lmm_lift))
  cat("  -> T0 signal is independent of LMM guidance — strong baseline marker.\n")
}

ki67_rank <- gate_rank_naive
if (!is.null(ki67_rank) && ki67_rank <= 5) {
  cat(sprintf("\n  KI67NAIVE ranks #%d in pure T0 scan — would be selected without LMM.\n",
              ki67_rank))
} else if (!is.null(ki67_rank)) {
  cat(sprintf("\n  KI67NAIVE ranks #%d in pure T0 scan — LMM guidance needed for selection.\n",
              ki67_rank))
}

cat("\n==========================================================\n")
cat("diag_14 complete.\n")
cat("Output files: diag_14_spaghetti.pdf, diag_14_boxplots.pdf,\n")
cat("              diag_14_scatter_t0_vs_delta.pdf\n")
cat("==========================================================\n")
