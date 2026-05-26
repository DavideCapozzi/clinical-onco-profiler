# =============================================================================
# diag_hnscc_03_immune_score_vs_cps.R
#
# PURPOSE: Pre-specified immune score analysis — can EFF (and composites)
#          discriminate R from NR better than CPS alone?
#
# DESIGN:
#   - Primary score: EFF_z (pre-specified: 100% fold selection in prior CV)
#   - Secondary: composite rank score rank(EFF) + rank(PD1) (non-parametric)
#   - "Range immunitario": tertile analysis with Cochran-Armitage trend test
#   - Combined model: LOO logistic(EFF_z + CPS_z) vs CPS alone
#   - NRI/IDI: combined vs CPS, and EFF alone vs CPS
#   - All AUC comparisons via DeLong (correlated samples, same patients)
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(pROC)
  library(ggplot2)
  library(gridExtra)
})

here_root <- function(...) {
  root <- tryCatch(here::here(), error = function(e) getwd())
  file.path(root, ...)
}

# ── 1. LOAD DATA ─────────────────────────────────────────────────────────────

rds  <- readRDS(here_root("results/HNSCC_Response/01_data_processing",
                           "data_processed_HNSCC_Response_standard.rds"))
z_df <- as.data.frame(rds$hybrid_data_z)           # includes Patient_ID col
meta <- rds$metadata

raw  <- readxl::read_excel(here_root("data/Dati_HNSCC_standardizzati_anonimi.xlsx"))
names(raw)[1]  <- "Patient_ID"
names(raw)[7]  <- "Therapy"
names(raw)[8]  <- "CPS"
names(raw)[21] <- "Response"
raw$CPS <- suppressWarnings(as.numeric(raw$CPS))

# Merge: keep only QC-passed patients
df <- merge(meta, raw[, c("Patient_ID","CPS","Response","Therapy")],
            by = "Patient_ID", all.x = TRUE)
df$y <- as.integer(df$Group == "R")   # 1=Responder, 0=Non-responder

# Extract marker z-scores
get_z <- function(marker) {
  stopifnot(marker %in% names(z_df))
  setNames(z_df[[marker]], z_df$Patient_ID)
}

eff_z   <- get_z("EFF")
pd1_z   <- get_z("PD1")
sbtla_z <- get_z("sBTLA")

# Align to df row order
df$EFF_z   <- eff_z[df$Patient_ID]
df$PD1_z   <- pd1_z[df$Patient_ID]
df$sBTLA_z <- sbtla_z[df$Patient_ID]
df$CPS_z   <- scale(df$CPS)[, 1]

# ── 2. COMPOSITE SCORES ───────────────────────────────────────────────────────

# Rank composite (non-parametric): higher = more immune-activated
df$rank_EFF_PD1   <- rank(df$EFF_z) + rank(df$PD1_z)
df$rank_EFF_sBTLA <- rank(df$EFF_z) + rank(df$sBTLA_z)
df$mean_EFF_PD1   <- (df$EFF_z + df$PD1_z) / 2

n      <- nrow(df)
n_R    <- sum(df$y == 1)
n_NR   <- sum(df$y == 0)

cat(sprintf("\n=== DATASET: n=%d (R=%d, NR=%d) ===\n", n, n_R, n_NR))

# ── 3. ROC + DeLong COMPARISONS ──────────────────────────────────────────────

# CPS: empirically lower CPS -> R in this cohort (AUC with higher->R = 0.449 < 0.5).
# Negate so direction="<" is consistent across all scores (higher = more likely R).
df$neg_CPS_z <- -df$CPS_z

roc_cps      <- roc(df$y, df$neg_CPS_z,       quiet = TRUE, direction = "<")
roc_eff      <- roc(df$y, df$EFF_z,            quiet = TRUE, direction = "<")
roc_pd1      <- roc(df$y, df$PD1_z,            quiet = TRUE, direction = "<")
roc_rank_ep  <- roc(df$y, df$rank_EFF_PD1,     quiet = TRUE, direction = "<")
roc_rank_es  <- roc(df$y, df$rank_EFF_sBTLA,   quiet = TRUE, direction = "<")
roc_mean_ep  <- roc(df$y, df$mean_EFF_PD1,     quiet = TRUE, direction = "<")

auc_ci <- function(r) {
  ci_v <- ci.auc(r, method = "bootstrap", boot.n = 2000)
  sprintf("%.3f [%.3f–%.3f]", ci_v[2], ci_v[1], ci_v[3])
}

cat("\n=== ROC AUC (bootstrap 95% CI) ===\n")
cat("CPS alone (best dir):", auc_ci(roc_cps), "\n")
cat("EFF_z (primary)    :", auc_ci(roc_eff), "\n")
cat("PD1_z              :", auc_ci(roc_pd1), "\n")
cat("rank(EFF)+rank(PD1):", auc_ci(roc_rank_ep), "\n")
cat("rank(EFF)+rank(sBTLA):", auc_ci(roc_rank_es), "\n")
cat("mean_z(EFF,PD1)    :", auc_ci(roc_mean_ep), "\n")

cat("\n=== DeLong tests vs CPS ===\n")
delong_test <- function(roc1, roc2, label) {
  t <- roc.test(roc1, roc2, method = "delong")
  cat(sprintf("  %-30s AUC_diff=%.3f  p=%.4f\n",
              label, as.numeric(auc(roc2)) - as.numeric(auc(roc1)), t$p.value))
}
delong_test(roc_cps, roc_eff,     "EFF vs CPS")
delong_test(roc_cps, roc_rank_ep, "rank(EFF+PD1) vs CPS")
delong_test(roc_cps, roc_rank_es, "rank(EFF+sBTLA) vs CPS")
delong_test(roc_cps, roc_mean_ep, "mean_z(EFF+PD1) vs CPS")

# Also EFF vs rank composites to choose best
cat("\n=== DeLong: composites vs EFF alone ===\n")
delong_test(roc_eff, roc_rank_ep, "rank(EFF+PD1) vs EFF")
delong_test(roc_eff, roc_rank_es, "rank(EFF+sBTLA) vs EFF")
delong_test(roc_eff, roc_mean_ep, "mean_z(EFF+PD1) vs EFF")

# ── 4. WILCOXON TESTS ─────────────────────────────────────────────────────────

cat("\n=== Wilcoxon R vs NR ===\n")
wtest <- function(x, label) {
  w <- wilcox.test(x[df$y==1], x[df$y==0], exact = FALSE)
  r_med <- median(x[df$y==1], na.rm=TRUE)
  nr_med <- median(x[df$y==0], na.rm=TRUE)
  cat(sprintf("  %-30s median(R)=%.2f  median(NR)=%.2f  p=%.4f\n",
              label, r_med, nr_med, w$p.value))
}
wtest(df$EFF_z,          "EFF_z")
wtest(df$rank_EFF_PD1,   "rank(EFF+PD1)")
wtest(df$rank_EFF_sBTLA, "rank(EFF+sBTLA)")
wtest(df$mean_EFF_PD1,   "mean_z(EFF+PD1)")
wtest(df$CPS_z,          "CPS_z")

# ── 5. RANGE IMMUNITARIO — TERTILE ANALYSIS ───────────────────────────────────

cat("\n=== RANGE IMMUNITARIO: EFF tertiles ===\n")

df$EFF_tertile <- with(df, cut(EFF_z,
  breaks = quantile(EFF_z, probs = c(0, 1/3, 2/3, 1)),
  labels = c("EFF-Low", "EFF-Mid", "EFF-High"),
  include.lowest = TRUE))

tbl_tert <- df %>%
  group_by(EFF_tertile) %>%
  summarise(
    n         = n(),
    n_R       = sum(y == 1),
    n_NR      = sum(y == 0),
    resp_rate = round(100 * mean(y == 1), 1),
    .groups   = "drop"
  )
print(as.data.frame(tbl_tert))

# Cochran-Armitage trend test
ca_scores <- c("EFF-Low"=1, "EFF-Mid"=2, "EFF-High"=3)
tbl_matrix <- matrix(
  c(tbl_tert$n_R, tbl_tert$n_NR),
  nrow = 2, byrow = TRUE
)
# prop.trend.test expects: x = successes, n = total, score = group scores
ca_test <- prop.trend.test(
  x = tbl_tert$n_R,
  n = tbl_tert$n,
  score = 1:3
)
cat(sprintf("\nCochran-Armitage trend test: chi2=%.3f, df=1, p=%.4f\n",
            ca_test$statistic, ca_test$p.value))

# Also for rank(EFF+PD1) composite
df$rank_tertile <- cut(df$rank_EFF_PD1,
  breaks = quantile(df$rank_EFF_PD1, probs = c(0, 1/3, 2/3, 1)),
  labels = c("Score-Low", "Score-Mid", "Score-High"),
  include.lowest = TRUE)

tbl_rank <- df %>%
  group_by(rank_tertile) %>%
  summarise(n=n(), n_R=sum(y==1), resp_rate=round(100*mean(y==1),1), .groups="drop")

cat("\n=== RANGE IMMUNITARIO: rank(EFF+PD1) tertiles ===\n")
print(as.data.frame(tbl_rank))

ca_rank <- prop.trend.test(x=tbl_rank$n_R, n=tbl_rank$n, score=1:3)
cat(sprintf("CA trend (rank composite): chi2=%.3f, p=%.4f\n",
            ca_rank$statistic, ca_rank$p.value))

# ── 6. COMBINED MODEL: EFF + CPS (LOO logistic) ──────────────────────────────

cat("\n=== COMBINED MODEL: LOO logistic (EFF_z + CPS_z) ===\n")

# LOO: fit on n-1, predict on held-out
loo_prob_combined <- numeric(n)
loo_prob_cps_only <- numeric(n)
loo_prob_eff_only <- numeric(n)

for (i in seq_len(n)) {
  train <- df[-i, ]
  test  <- df[i,  ]

  m_comb <- glm(y ~ EFF_z + neg_CPS_z, data = train, family = binomial())
  m_cps  <- glm(y ~ neg_CPS_z,         data = train, family = binomial())
  m_eff  <- glm(y ~ EFF_z,             data = train, family = binomial())

  loo_prob_combined[i] <- predict(m_comb, newdata=test, type="response")
  loo_prob_cps_only[i] <- predict(m_cps,  newdata=test, type="response")
  loo_prob_eff_only[i] <- predict(m_eff,  newdata=test, type="response")
}

roc_loo_comb <- roc(df$y, loo_prob_combined, quiet=TRUE)
roc_loo_cps  <- roc(df$y, loo_prob_cps_only, quiet=TRUE)
roc_loo_eff  <- roc(df$y, loo_prob_eff_only, quiet=TRUE)

cat("LOO AUC — CPS only    :", round(auc(roc_loo_cps),  3), "\n")
cat("LOO AUC — EFF only    :", round(auc(roc_loo_eff),  3), "\n")
cat("LOO AUC — EFF + CPS   :", round(auc(roc_loo_comb), 3), "\n")

cat("\nDeLong (LOO): EFF+CPS vs CPS\n")
dl_comb_vs_cps <- roc.test(roc_loo_cps, roc_loo_comb, method="delong")
cat(sprintf("  delta=%.3f  p=%.4f\n",
    auc(roc_loo_comb) - auc(roc_loo_cps), dl_comb_vs_cps$p.value))

cat("DeLong (LOO): EFF vs CPS\n")
dl_eff_vs_cps <- roc.test(roc_loo_cps, roc_loo_eff, method="delong")
cat(sprintf("  delta=%.3f  p=%.4f\n",
    auc(roc_loo_eff) - auc(roc_loo_cps), dl_eff_vs_cps$p.value))

# ── 7. NRI / IDI (bootstrap) ─────────────────────────────────────────────────

cat("\n=== NRI / IDI: combined vs CPS (bootstrap B=2000) ===\n")

set.seed(42)
B <- 2000
idi_vals <- numeric(B)
nri_vals <- numeric(B)

for (b in seq_len(B)) {
  idx <- sample(n, n, replace=TRUE)
  d_b <- df[idx, ]

  m_new <- glm(y ~ EFF_z + neg_CPS_z, data=d_b, family=binomial())
  m_ref <- glm(y ~ neg_CPS_z,         data=d_b, family=binomial())

  p_new <- predict(m_new, type="response")
  p_ref <- predict(m_ref, type="response")

  idi_b <- mean(p_new[d_b$y==1] - p_ref[d_b$y==1]) -
           mean(p_new[d_b$y==0] - p_ref[d_b$y==0])
  idi_vals[b] <- idi_b

  # continuous NRI
  up_event   <- mean((p_new > p_ref)[d_b$y==1])
  down_event <- mean((p_new < p_ref)[d_b$y==1])
  up_non     <- mean((p_new > p_ref)[d_b$y==0])
  down_non   <- mean((p_new < p_ref)[d_b$y==0])
  nri_vals[b] <- (up_event - down_event) - (up_non - down_non)
}

idi_obs <- mean(loo_prob_combined[df$y==1] - loo_prob_cps_only[df$y==1]) -
           mean(loo_prob_combined[df$y==0] - loo_prob_cps_only[df$y==0])
nri_obs <- (mean((loo_prob_combined > loo_prob_cps_only)[df$y==1]) -
             mean((loo_prob_combined < loo_prob_cps_only)[df$y==1])) -
           (mean((loo_prob_combined > loo_prob_cps_only)[df$y==0]) -
             mean((loo_prob_combined < loo_prob_cps_only)[df$y==0]))

idi_ci <- quantile(idi_vals, c(0.025, 0.975))
nri_ci <- quantile(nri_vals, c(0.025, 0.975))
idi_p  <- 2 * min(mean(idi_vals <= 0), mean(idi_vals >= 0))
nri_p  <- 2 * min(mean(nri_vals <= 0), mean(nri_vals >= 0))

cat(sprintf("IDI = %.3f [%.3f–%.3f]  p=%.4f\n",
            idi_obs, idi_ci[1], idi_ci[2], idi_p))
cat(sprintf("cNRI= %.3f [%.3f–%.3f]  p=%.4f\n",
            nri_obs, nri_ci[1], nri_ci[2], nri_p))

# ── 8. OPTIMAL CUT-POINT TABLE ────────────────────────────────────────────────

cat("\n=== OPTIMAL CUT-POINTS (Youden) ===\n")

get_cutpoint <- function(roc_obj, label) {
  co <- coords(roc_obj, "best", best.method="youden",
               ret=c("threshold","sensitivity","specificity","ppv","npv"),
               transpose=FALSE)
  cat(sprintf("  %-30s cut=%.2f  Sn=%.2f  Sp=%.2f  PPV=%.2f  NPV=%.2f\n",
              label,
              co$threshold[1], co$sensitivity[1],
              co$specificity[1], co$ppv[1], co$npv[1]))
}

get_cutpoint(roc_cps,     "CPS alone")
get_cutpoint(roc_eff,     "EFF_z (full data)")
get_cutpoint(roc_rank_ep, "rank(EFF+PD1)")

# CPS ≥20 cut-point (clinical standard)
cps20 <- df %>%
  summarise(
    label      = "CPS >= 20 (clinical)",
    Sn         = mean((CPS >= 20)[y == 1]),
    Sp         = mean((CPS < 20) [y == 0]),
    PPV        = mean(y[CPS >= 20] == 1),
    NPV        = mean(y[CPS < 20] == 0)
  )
cat(sprintf("  %-30s cut=20     Sn=%.2f  Sp=%.2f  PPV=%.2f  NPV=%.2f\n",
            "CPS>=20 (clinical)", cps20$Sn, cps20$Sp, cps20$PPV, cps20$NPV))

# ── 9. 2×2 GRID: EFF × CPS ────────────────────────────────────────────────────

cat("\n=== 2x2 GRID: EFF (median split) x CPS (>=20) ===\n")

eff_med <- median(df$EFF_z)
df$EFF_high <- df$EFF_z >= eff_med
df$CPS_high <- df$CPS >= 20

grid <- df %>%
  mutate(category = case_when(
    EFF_high & CPS_high   ~ "EFF-high + CPS-high",
    EFF_high & !CPS_high  ~ "EFF-high + CPS-low",
    !EFF_high & CPS_high  ~ "EFF-low + CPS-high",
    TRUE                  ~ "EFF-low + CPS-low"
  )) %>%
  group_by(category) %>%
  summarise(n=n(), n_R=sum(y==1),
            resp_rate=round(100*mean(y==1),1), .groups="drop") %>%
  arrange(desc(resp_rate))
print(as.data.frame(grid))

# Fisher: EFF-high vs EFF-low regardless of CPS
fisher_eff <- fisher.test(table(df$EFF_high, df$y))
cat(sprintf("\nFisher EFF (median split): OR=%.2f  p=%.4f\n",
            fisher_eff$estimate, fisher_eff$p.value))

fisher_cps <- fisher.test(table(df$CPS_high, df$y))
cat(sprintf("Fisher CPS >= 20:          OR=%.2f  p=%.4f\n",
            fisher_cps$estimate, fisher_cps$p.value))

# ── 10. PLOTS ─────────────────────────────────────────────────────────────────

out_dir <- here_root("results/HNSCC_Response/diagnostics")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# --- Plot 1: ROC curves comparison ---
pdf(file.path(out_dir, "HNSCC_ImmuneScore_ROC.pdf"), width=6, height=6)

plot(roc_cps, col="#888888", lwd=2, main="ROC Curves: Immune Score vs CPS",
     xlab="1 - Specificity", ylab="Sensitivity")
plot(roc_eff,     col="#E41A1C", lwd=2.5, add=TRUE)
plot(roc_rank_ep, col="#377EB8", lwd=2,   add=TRUE, lty=2)
plot(roc_loo_comb,col="#4DAF4A", lwd=2.5, add=TRUE, lty=1)

legend("bottomright", bty="n", lwd=c(2,2.5,2,2.5), lty=c(1,1,2,1),
       col=c("#888888","#E41A1C","#377EB8","#4DAF4A"),
       legend=c(
         sprintf("CPS alone (AUC=%.3f)", auc(roc_cps)),
         sprintf("EFF (AUC=%.3f)", auc(roc_eff)),
         sprintf("rank(EFF+PD1) (AUC=%.3f)", auc(roc_rank_ep)),
         sprintf("EFF+CPS LOO (AUC=%.3f)", auc(roc_loo_comb))
       ), cex=0.82)
dev.off()

# --- Plot 2: Tertile bar chart ---
tbl_plot <- rbind(
  tbl_tert %>% mutate(score="EFF alone",
    category=as.character(EFF_tertile)) %>%
    select(score, category, n, n_R, resp_rate),
  tbl_rank %>% mutate(score="rank(EFF+PD1)",
    category=as.character(rank_tertile)) %>%
    select(score, category, n, n_R, resp_rate)
)
tbl_plot$category <- factor(tbl_plot$category,
  levels=c("EFF-Low","EFF-Mid","EFF-High",
           "Score-Low","Score-Mid","Score-High"))

p_tert <- ggplot(tbl_plot, aes(x=category, y=resp_rate, fill=score)) +
  geom_col(position="dodge", width=0.7, alpha=0.85) +
  geom_text(aes(label=sprintf("%.0f%%\n(n=%d)", resp_rate, n)),
            position=position_dodge(width=0.7), vjust=-0.3, size=3) +
  scale_fill_manual(values=c("EFF alone"="#E41A1C","rank(EFF+PD1)"="#377EB8")) +
  labs(title="Response rate by immune tertile ('Range Immunitario')",
       x=NULL, y="Response rate (%)", fill="Score") +
  theme_bw(base_size=11) +
  coord_cartesian(ylim=c(0,80)) +
  theme(legend.position="top")

pdf(file.path(out_dir, "HNSCC_ImmuneRange_Tertiles.pdf"), width=7, height=5)
print(p_tert)
dev.off()

# --- Plot 3: EFF boxplot R vs NR ---
p_box <- ggplot(df, aes(x=Group, y=EFF_z, fill=Group)) +
  geom_boxplot(outlier.shape=21, alpha=0.8, width=0.5) +
  geom_jitter(width=0.15, alpha=0.5, size=1.5) +
  scale_fill_manual(values=c("NR"="#4393C3","R"="#D6604D")) +
  labs(title="EFF (z-score) by clinical response",
       subtitle=sprintf("Wilcoxon p=%.4f  |  AUC=%.3f",
                        wilcox.test(EFF_z~Group,data=df,exact=FALSE)$p.value,
                        as.numeric(auc(roc_eff))),
       x=NULL, y="EFF z-score") +
  theme_bw(base_size=11) + theme(legend.position="none")

# --- Plot 4: 2x2 grid heatmap ---
grid_wide <- df %>%
  mutate(
    EFF_cat = ifelse(EFF_high, "EFF-high", "EFF-low"),
    CPS_cat = ifelse(CPS_high, "CPS>=20", "CPS<20")
  ) %>%
  group_by(EFF_cat, CPS_cat) %>%
  summarise(n=n(), n_R=sum(y==1),
            resp_rate=round(100*mean(y==1),1), .groups="drop")

p_grid <- ggplot(grid_wide, aes(x=CPS_cat, y=EFF_cat, fill=resp_rate)) +
  geom_tile(color="white", linewidth=1) +
  geom_text(aes(label=sprintf("%.0f%%\nn=%d", resp_rate, n)),
            size=4.5, fontface="bold") +
  scale_fill_gradient2(low="#4393C3", mid="#FFFFFF", high="#D6604D",
                       midpoint=42, name="Resp.\nrate (%)") +
  labs(title="Response rate: EFF (median) × CPS threshold",
       x="CPS category", y="EFF category") +
  theme_bw(base_size=11)

pdf(file.path(out_dir, "HNSCC_ImmuneScore_Panels.pdf"), width=10, height=5)
grid.arrange(p_box, p_grid, ncol=2)
dev.off()

cat("\n=== DONE — plots saved to results/HNSCC_Response/diagnostics/ ===\n")
