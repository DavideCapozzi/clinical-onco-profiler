# diagnostics/diag_hnscc_02_signal_and_models.R
# Diagnostic battery for HNSCC_Response:
#   1. Signal characterisation: permutation of max-AUC, Cohen's d + bootstrap
#   2. Power analysis: n needed for 80% power given observed effect sizes
#   3. Alternative models (LOO-unbiased): Ridge, Random Forest, PCA-Logistic, forced 2-marker
#   4. Unsupervised structure: PCA biplot, correlation heatmap of top markers
#
# All multivariate models use strict LOO (no data leakage).
# Univariate permutation test: max-AUC across 30 markers under label permutation.

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(pROC)
  library(glmnet)
  library(e1071)
  library(ranger)
  library(ggplot2)
  library(patchwork)
  library(factoextra)
})

options(crayon.enabled = FALSE)
set.seed(2026)

out_dir <- here("results/HNSCC_Response/diagnostics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── 0. Load Step 01 output ────────────────────────────────────────────────────
rds_path <- here("results/HNSCC_Response/01_data_processing/data_processed_HNSCC_Response_standard.rds")
if (!file.exists(rds_path)) stop("[FATAL] Step 01 RDS not found. Run Pass 1 first.")
DATA <- readRDS(rds_path)

markers <- DATA$hybrid_markers
X_full  <- as.matrix(DATA$hybrid_data_z[, markers, drop = FALSE])
y       <- DATA$metadata$Group   # "R" / "NR"
pos_lab <- "R"
y_bin   <- as.integer(y == pos_lab)

n  <- nrow(X_full)
p  <- ncol(X_full)
nR <- sum(y_bin)
nN <- n - nR

cat(sprintf("\n[Data] n=%d  p=%d  R=%d  NR=%d\n", n, p, nR, nN))

# Helper: flip AUC < 0.5 to max(AUC, 1-AUC)
max_auc <- function(roc_obj) max(as.numeric(pROC::auc(roc_obj)), 1 - as.numeric(pROC::auc(roc_obj)))

# Helper: safe Wilcoxon AUC for a single vector
marker_auc <- function(x, y_bin) {
  r1 <- tryCatch(pROC::roc(y_bin, x, direction = "<", quiet = TRUE), error = function(e) NULL)
  r2 <- tryCatch(pROC::roc(y_bin, x, direction = ">", quiet = TRUE), error = function(e) NULL)
  a1 <- if (!is.null(r1)) as.numeric(pROC::auc(r1)) else 0
  a2 <- if (!is.null(r2)) as.numeric(pROC::auc(r2)) else 0
  max(a1, a2)
}

# ── 1. SIGNAL CHARACTERISATION ────────────────────────────────────────────────
cat("\n=== [1] SIGNAL CHARACTERISATION ===\n")

# 1a. Full-data Wilcoxon + Cohen's d (bootstrap 5000 resamples)
B_boot <- 5000
sig_res <- purrr::map_dfr(markers, function(mk) {
  x   <- X_full[, mk]
  xR  <- x[y_bin == 1]; xN <- x[y_bin == 0]
  wt  <- tryCatch(wilcox.test(xR, xN, exact = FALSE), error = function(e) list(p.value = NA))
  d_obs <- (mean(xR) - mean(xN)) / sd(x)
  auc_v <- marker_auc(x, y_bin)

  # Bootstrap CI for Cohen's d
  d_boot <- replicate(B_boot, {
    bR <- sample(xR, length(xR), replace = TRUE)
    bN <- sample(xN, length(xN), replace = TRUE)
    (mean(bR) - mean(bN)) / sd(c(bR, bN))
  })
  data.frame(
    Marker    = mk,
    d         = round(d_obs, 3),
    d_lo      = round(quantile(d_boot, 0.025), 3),
    d_hi      = round(quantile(d_boot, 0.975), 3),
    raw_p     = wt$p.value,
    fdr       = NA_real_,
    full_auc  = round(auc_v, 3),
    stringsAsFactors = FALSE
  )
})
sig_res$fdr <- p.adjust(sig_res$raw_p, method = "BH")
sig_res <- sig_res[order(sig_res$raw_p), ]

cat("\nTop 10 markers by Wilcoxon raw p (full data, n=52):\n")
print(head(sig_res[, c("Marker","d","d_lo","d_hi","raw_p","fdr","full_auc")], 10),
      row.names = FALSE, digits = 3)

# 1b. Permutation of max single-marker AUC under label permutation
cat("\n--- Max-AUC permutation test (n_perm=1999) ---\n")
N_PERM <- 1999
obs_max_auc <- max(sig_res$full_auc)
perm_max <- replicate(N_PERM, {
  y_perm <- sample(y_bin)
  max(vapply(markers, function(mk) marker_auc(X_full[, mk], y_perm), numeric(1)))
})
p_perm_max <- mean(perm_max >= obs_max_auc)
cat(sprintf("Observed best single-marker AUC (EFF): %.3f\n", obs_max_auc))
cat(sprintf("Permutation p (max-AUC test, %d perms): %.4f\n", N_PERM, p_perm_max))
cat(sprintf("95th pctile of null max-AUC: %.3f\n", quantile(perm_max, 0.95)))
cat(sprintf("99th pctile of null max-AUC: %.3f\n", quantile(perm_max, 0.99)))

# ── 2. POWER ANALYSIS ────────────────────────────────────────────────────────
cat("\n=== [2] POWER ANALYSIS ===\n")
# Wilcoxon rank-sum power via normal approximation (pwr::pwr.norm.test on Cohen's d)
# Cohen's d → AUC relationship: AUC ≈ pnorm(d/sqrt(2))
# For each top marker, compute n needed for 80% power (two-sided alpha=0.05)
top5 <- head(sig_res, 5)
cat("\nN per group needed for 80% power (alpha=0.05, two-sided Wilcoxon approx.):\n")
for (i in seq_len(nrow(top5))) {
  d_abs <- abs(top5$d[i])
  if (d_abs < 0.01) { cat(sprintf("  %-20s  d=%.3f  →  n=Inf (no effect)\n", top5$Marker[i], d_abs)); next }
  # approx formula: n = 2*(z_alpha/2 + z_beta)^2 / d^2, z_beta=0.842 (80%), z_alpha=1.96
  n_needed <- ceiling(2 * ((1.96 + 0.842)^2) / d_abs^2)
  cat(sprintf("  %-20s  d=%.3f  CI[%.3f, %.3f]  →  n=%d per group (total=%d)\n",
              top5$Marker[i], d_abs, top5$d_lo[i], top5$d_hi[i],
              n_needed, 2 * n_needed))
}
cat(sprintf("\n  Current study: nR=%d, nNR=%d (total=%d)\n", nR, nN, n))

# ── 3. ALTERNATIVE MODELS (LOO) ───────────────────────────────────────────────
cat("\n=== [3] ALTERNATIVE MODELS (strict LOO) ===\n")

# Shared inner CV utility: k-fold inner on training fold
inner_auc_glmnet <- function(X_tr, y_tr, alpha_val, n_lambda = 50, inner_k = 5) {
  if (length(unique(y_tr)) < 2) return(list(lambda = 0.01, auc = 0.5))
  fit <- tryCatch(
    glmnet::cv.glmnet(X_tr, y_tr, family = "binomial", alpha = alpha_val,
                      nfolds = inner_k, type.measure = "auc", nlambda = n_lambda),
    error = function(e) NULL
  )
  if (is.null(fit)) return(list(lambda = 0.01, auc = 0.5))
  list(lambda = fit$lambda.min, auc = max(fit$cvm, na.rm = TRUE))
}

# --- 3a. Ridge Regression (alpha=0) LOO ---
cat("\n--- 3a. Ridge Regression (alpha=0) LOO ---\n")
ridge_probs <- numeric(n)
for (i in seq_len(n)) {
  X_tr <- X_full[-i, , drop = FALSE]; y_tr <- as.integer(y[-i] == pos_lab)
  X_te <- X_full[i,  , drop = FALSE]
  res  <- inner_auc_glmnet(X_tr, y_tr, alpha_val = 0)
  fit  <- glmnet::glmnet(X_tr, y_tr, family = "binomial", alpha = 0, lambda = res$lambda)
  ridge_probs[i] <- as.numeric(predict(fit, newx = X_te, type = "response"))
}
roc_ridge   <- pROC::roc(y_bin, ridge_probs, direction = "<", quiet = TRUE)
auc_ridge   <- as.numeric(pROC::auc(roc_ridge))
ci_ridge    <- as.numeric(pROC::ci.auc(roc_ridge, method = "delong"))
perm_ridge  <- replicate(999, {
  as.numeric(pROC::auc(pROC::roc(sample(y_bin), ridge_probs, quiet = TRUE)))
})
p_ridge <- mean(perm_ridge >= auc_ridge)
cat(sprintf("  Ridge LOO AUC = %.3f [%.3f–%.3f]  perm p = %.3f\n",
            auc_ridge, ci_ridge[1], ci_ridge[3], p_ridge))

# --- 3b. Random Forest (ranger) LOO with OOB importance ---
cat("\n--- 3b. Random Forest (ranger) LOO ---\n")
rf_probs <- numeric(n)
for (i in seq_len(n)) {
  X_tr  <- X_full[-i, , drop = FALSE]
  y_tr  <- factor(y[-i], levels = c("NR","R"))
  X_te  <- X_full[i,  , drop = FALSE]
  df_tr <- as.data.frame(X_tr); df_tr$outcome <- y_tr
  fit_rf <- ranger::ranger(outcome ~ ., data = df_tr,
                            num.trees = 500, probability = TRUE,
                            class.weights = c(NR = nR / n, R = nN / n),
                            seed = 2026 + i, verbose = FALSE)
  df_te <- as.data.frame(X_te)
  rf_probs[i] <- predict(fit_rf, data = df_te)$predictions[, "R"]
}
roc_rf  <- pROC::roc(y_bin, rf_probs, direction = "<", quiet = TRUE)
auc_rf  <- as.numeric(pROC::auc(roc_rf))
ci_rf   <- as.numeric(pROC::ci.auc(roc_rf, method = "delong"))
perm_rf <- replicate(999, {
  as.numeric(pROC::auc(pROC::roc(sample(y_bin), rf_probs, quiet = TRUE)))
})
p_rf <- mean(perm_rf >= auc_rf)
cat(sprintf("  RF LOO AUC = %.3f [%.3f–%.3f]  perm p = %.3f\n",
            auc_rf, ci_rf[1], ci_rf[3], p_rf))

# RF variable importance (full model)
df_all    <- as.data.frame(X_full); df_all$outcome <- factor(y, levels = c("NR","R"))
fit_rf_all <- ranger::ranger(outcome ~ ., data = df_all,
                              num.trees = 1000, probability = TRUE,
                              importance = "impurity_corrected",
                              class.weights = c(NR = nR / n, R = nN / n),
                              seed = 2026, verbose = FALSE)
imp <- sort(fit_rf_all$variable.importance, decreasing = TRUE)
cat("\n  RF variable importance (top 10, impurity-corrected):\n")
print(round(head(imp, 10), 4))

# --- 3c. PCA-Logistic LOO (reduce to PC1-3 then logistic) ---
cat("\n--- 3c. PCA-Logistic LOO (PC1-3) ---\n")
pca_probs <- numeric(n)
for (i in seq_len(n)) {
  X_tr  <- X_full[-i, , drop = FALSE]
  X_te  <- X_full[i,  , drop = FALSE]
  y_tr  <- y_bin[-i]
  # Fit PCA on training fold only
  pca_tr <- prcomp(X_tr, center = TRUE, scale. = FALSE)  # already z-scored
  n_pc   <- min(3L, ncol(pca_tr$rotation))
  Xpc_tr <- pca_tr$x[, seq_len(n_pc), drop = FALSE]
  pc_names <- paste0("PC", seq_len(n_pc))
  colnames(Xpc_tr) <- pc_names
  te_scores <- as.numeric(predict(pca_tr, newdata = X_te))
  Xpc_te <- as.data.frame(matrix(te_scores[seq_len(n_pc)], nrow = 1,
                                   dimnames = list(NULL, pc_names)))
  df_pc  <- as.data.frame(Xpc_tr); df_pc$y <- y_tr
  fit_lr <- glm(y ~ ., data = df_pc, family = binomial, control = list(maxit = 100))
  pca_probs[i] <- predict(fit_lr, newdata = Xpc_te, type = "response")
}
roc_pca  <- pROC::roc(y_bin, pca_probs, direction = "<", quiet = TRUE)
auc_pca  <- as.numeric(pROC::auc(roc_pca))
ci_pca   <- as.numeric(pROC::ci.auc(roc_pca, method = "delong"))
perm_pca <- replicate(999, {
  as.numeric(pROC::auc(pROC::roc(sample(y_bin), pca_probs, quiet = TRUE)))
})
p_pca <- mean(perm_pca >= auc_pca)
cat(sprintf("  PCA-Logistic LOO AUC = %.3f [%.3f–%.3f]  perm p = %.3f\n",
            auc_pca, ci_pca[1], ci_pca[3], p_pca))

# --- 3d. Forced 2-marker model: EFF + sBTLA (top two by LOO-BalAcc) ---
cat("\n--- 3d. Forced 2-marker (EFF + sBTLA) Logistic LOO ---\n")
forced_markers <- c("EFF", "sBTLA")
forced_probs   <- numeric(n)
for (i in seq_len(n)) {
  X_tr <- X_full[-i, forced_markers, drop = FALSE]
  X_te <- X_full[i,  forced_markers, drop = FALSE]
  y_tr <- y_bin[-i]
  df_tr <- as.data.frame(cbind(X_tr, y = y_tr))
  fit_f <- tryCatch(glm(y ~ ., data = df_tr, family = binomial, control = list(maxit = 100)),
                    error = function(e) NULL)
  if (is.null(fit_f)) { forced_probs[i] <- 0.5; next }
  forced_probs[i] <- predict(fit_f, newdata = as.data.frame(X_te), type = "response")
}
roc_f2  <- pROC::roc(y_bin, forced_probs, direction = "<", quiet = TRUE)
auc_f2  <- as.numeric(pROC::auc(roc_f2))
ci_f2   <- as.numeric(pROC::ci.auc(roc_f2, method = "delong"))
perm_f2 <- replicate(999, {
  as.numeric(pROC::auc(pROC::roc(sample(y_bin), forced_probs, quiet = TRUE)))
})
p_f2 <- mean(perm_f2 >= auc_f2)
cat(sprintf("  2-marker LOO AUC = %.3f [%.3f–%.3f]  perm p = %.3f\n",
            auc_f2, ci_f2[1], ci_f2[3], p_f2))

# Also try EFF alone
eff_probs <- numeric(n)
for (i in seq_len(n)) {
  x_tr <- X_full[-i, "EFF"]; y_tr <- y_bin[-i]
  x_te <- X_full[i,  "EFF"]
  df_tr <- data.frame(x = x_tr, y = y_tr)
  fit_e <- tryCatch(glm(y ~ x, data = df_tr, family = binomial, control = list(maxit = 100)),
                    error = function(e) NULL)
  eff_probs[i] <- if (!is.null(fit_e)) predict(fit_e, newdata = data.frame(x = x_te), type = "response") else 0.5
}
roc_eff  <- pROC::roc(y_bin, eff_probs, direction = "<", quiet = TRUE)
auc_eff  <- as.numeric(pROC::auc(roc_eff))
ci_eff   <- as.numeric(pROC::ci.auc(roc_eff, method = "delong"))
perm_eff <- replicate(999, as.numeric(pROC::auc(pROC::roc(sample(y_bin), eff_probs, quiet=TRUE))))
p_eff    <- mean(perm_eff >= auc_eff)
cat(sprintf("  EFF alone LOO AUC = %.3f [%.3f–%.3f]  perm p = %.3f\n",
            auc_eff, ci_eff[1], ci_eff[3], p_eff))

# ── 4. SUMMARY TABLE ─────────────────────────────────────────────────────────
cat("\n=== SUMMARY: All models (LOOCV AUC) ===\n")
summary_df <- data.frame(
  Model      = c("EN (original)","SVM-RBF (original)","Ridge (alpha=0)","Random Forest",
                 "PCA-Logistic (PC1-3)","2-marker (EFF+sBTLA)","EFF alone"),
  LOO_AUC    = round(c(0.277, 0.486, auc_ridge, auc_rf, auc_pca, auc_f2, auc_eff), 3),
  CI_Lo      = round(c(0.129, 0.326, ci_ridge[1], ci_rf[1], ci_pca[1], ci_f2[1], ci_eff[1]), 3),
  CI_Hi      = round(c(0.426, 0.647, ci_ridge[3], ci_rf[3], ci_pca[3], ci_f2[3], ci_eff[3]), 3),
  Perm_p     = round(c(0.996, 0.567, p_ridge, p_rf, p_pca, p_f2, p_eff), 3),
  stringsAsFactors = FALSE
)
print(summary_df, row.names = FALSE)

# ── 5. PLOTS ──────────────────────────────────────────────────────────────────
cat("\n=== [5] GENERATING PLOTS ===\n")

# 5a. Volcano: Cohen's d vs -log10(raw p)
sig_res$sig <- sig_res$fdr < 0.20
p_volcano <- ggplot(sig_res, aes(x = d, y = -log10(raw_p),
                                  colour = sig, label = Marker)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = c(-0.3, 0.3),  linetype = "dashed", colour = "grey80") +
  geom_point(size = 2.5) +
  ggrepel::geom_text_repel(data = head(sig_res, 8), size = 3, max.overlaps = 12) +
  scale_colour_manual(values = c("FALSE" = "grey60", "TRUE" = "#B2182B"),
                      labels = c("FDR≥0.20","FDR<0.20"), name = "") +
  labs(title = "HNSCC: Marker effect sizes (full data, n=52)",
       x = "Cohen's d (R − NR)", y = "−log₁₀(Wilcoxon raw p)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

# 5b. Permutation null distribution of max single-marker AUC
perm_df <- data.frame(max_auc = perm_max)
p_perm <- ggplot(perm_df, aes(x = max_auc)) +
  geom_histogram(bins = 50, fill = "grey70", colour = "white") +
  geom_vline(xintercept = obs_max_auc, colour = "#B2182B", linewidth = 1.2) +
  annotate("text", x = obs_max_auc + 0.005, y = Inf,
           label = sprintf("EFF AUC=%.3f\np=%.3f", obs_max_auc, p_perm_max),
           vjust = 1.5, hjust = 0, colour = "#B2182B", size = 3.5) +
  labs(title = "Null distribution: max single-marker AUC",
       subtitle = sprintf("1999 label permutations, n=52, p=30"),
       x = "Best single-marker AUC (permuted labels)", y = "Count") +
  theme_bw(base_size = 12)

# 5c. Model comparison forest-style
summary_df_plot <- summary_df
summary_df_plot$Model <- factor(summary_df_plot$Model, levels = rev(summary_df_plot$Model))
summary_df_plot$sig   <- summary_df_plot$Perm_p < 0.05
p_models <- ggplot(summary_df_plot, aes(x = LOO_AUC, y = Model, colour = sig)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = CI_Lo, xmax = CI_Hi), height = 0.25, linewidth = 0.8) +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("%.3f (p=%.3f)", LOO_AUC, Perm_p)),
            hjust = -0.15, size = 3, colour = "black") +
  scale_colour_manual(values = c("FALSE" = "grey50", "TRUE" = "#2E8B57"),
                      labels = c("perm p≥0.05","perm p<0.05"), name = "") +
  scale_x_continuous(limits = c(0.1, 0.95)) +
  labs(title = "HNSCC: LOO AUC comparison across models",
       x = "LOO AUC [95% CI DeLong]", y = "") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

# 5d. PCA biplot (top markers)
pca_full <- prcomp(X_full, center = TRUE, scale. = FALSE)
pc_var   <- round(100 * pca_full$sdev^2 / sum(pca_full$sdev^2), 1)
pca_df   <- data.frame(PC1 = pca_full$x[,1], PC2 = pca_full$x[,2], Group = y)
p_pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = Group, shape = Group)) +
  stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.10, colour = NA) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = c("R" = "#2E8B57", "NR" = "#B2182B")) +
  scale_fill_manual(values   = c("R" = "#2E8B57", "NR" = "#B2182B")) +
  labs(title = "HNSCC: PCA biplot (PC1 × PC2)",
       x = sprintf("PC1 (%.1f%%)", pc_var[1]),
       y = sprintf("PC2 (%.1f%%)", pc_var[2])) +
  theme_bw(base_size = 12)

# 5e. Top-marker boxplots (top 6 by raw p)
top6 <- head(sig_res$Marker, 6)
box_base <- as.data.frame(X_full[, top6, drop = FALSE])
box_base$Group <- y
box_df <- tidyr::pivot_longer(box_base, cols = all_of(top6),
                               names_to = "Marker", values_to = "Value")
box_df$Group <- factor(box_df$Group, levels = c("R", "NR"))
p_box <- ggplot(box_df, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(outlier.size = 1.2, linewidth = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  scale_fill_manual(values = c("R" = "#2E8B57", "NR" = "#B2182B")) +
  facet_wrap(~ Marker, scales = "free_y", nrow = 2) +
  labs(title = "HNSCC: Top 6 markers by Wilcoxon p (z-scored)",
       x = "", y = "Z-score") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", strip.background = element_blank())

# Save all plots to one PDF
pdf(file.path(out_dir, "HNSCC_Diagnostics_Signal_Models.pdf"), width = 12, height = 16)
print(p_volcano)
print(p_perm)
print(p_models)
print(p_pca_plot)
print(p_box)
dev.off()
cat(sprintf("[Output] PDF: %s/HNSCC_Diagnostics_Signal_Models.pdf\n", out_dir))

# ── 6. RF IMPORTANCE PLOT ─────────────────────────────────────────────────────
imp_df <- data.frame(Marker = names(imp), Importance = as.numeric(imp))
imp_df <- imp_df[order(-imp_df$Importance), ][1:15, ]
imp_df$Marker <- factor(imp_df$Marker, levels = rev(imp_df$Marker))
p_imp <- ggplot(imp_df, aes(x = Importance, y = Marker,
                             fill = Importance > 0)) +
  geom_col() +
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "#1F77B4")) +
  geom_vline(xintercept = 0, colour = "black") +
  labs(title = "RF: Variable importance (impurity-corrected, full data)",
       x = "Impurity-corrected importance", y = "") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

pdf(file.path(out_dir, "HNSCC_RF_Importance.pdf"), width = 8, height = 6)
print(p_imp)
dev.off()

# ── 7. SAVE TEXT SUMMARY ─────────────────────────────────────────────────────
sink(file.path(out_dir, "HNSCC_Diagnostics_Summary.txt"))
cat("=== HNSCC_Response Diagnostic Summary ===\n")
cat(sprintf("Date: %s\n\n", Sys.time()))
cat(sprintf("Data: n=%d  p=%d  R=%d  NR=%d\n\n", n, p, nR, nN))

cat("--- Signal: Full-data Wilcoxon + Cohen's d (top 10) ---\n")
print(head(sig_res[, c("Marker","d","d_lo","d_hi","raw_p","fdr","full_auc")], 10),
      row.names = FALSE, digits = 3)

cat(sprintf("\n--- Max-AUC permutation test ---\n"))
cat(sprintf("Observed best AUC: %.3f  Perm p: %.4f  (n_perm=%d)\n",
            obs_max_auc, p_perm_max, N_PERM))
cat(sprintf("95th pctile null: %.3f | 99th pctile: %.3f\n\n",
            quantile(perm_max, 0.95), quantile(perm_max, 0.99)))

cat("--- Power analysis (n per group for 80% power, alpha=0.05) ---\n")
for (i in seq_len(nrow(top5))) {
  d_abs <- abs(top5$d[i])
  n_needed <- if (d_abs > 0.01) ceiling(2 * (1.96 + 0.842)^2 / d_abs^2) else Inf
  cat(sprintf("  %-20s  d=%.3f [%.3f,%.3f]  n=%d/group (total=%d)\n",
              top5$Marker[i], d_abs, top5$d_lo[i], top5$d_hi[i],
              n_needed, 2*n_needed))
}

cat("\n--- Model comparison (LOO AUC) ---\n")
print(summary_df, row.names = FALSE)

cat("\n--- RF variable importance (top 15, impurity-corrected) ---\n")
print(round(head(sort(fit_rf_all$variable.importance, decreasing=TRUE), 15), 4))
sink()

cat(sprintf("[Output] Text: %s/HNSCC_Diagnostics_Summary.txt\n", out_dir))
cat("\n=== DIAGNOSTICS COMPLETE ===\n")
