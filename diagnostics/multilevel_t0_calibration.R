# multilevel_t0_calibration.R
# ==============================================================================
# Clinical-utility layer on top of the T0 composite Ki67 score (the pre-specified
# pre-treatment biomarker that IS significant: AUC=0.66, perm p=0.0085 - see
# multilevel_t0_utility.R). Discrimination alone (AUC) is not enough for a
# clinical-prediction manuscript; reviewers (TRIPOD) demand:
#   (1) CALIBRATION          - do predicted probabilities match observed risk?
#   (2) DECISION CURVE / NB  - is the model clinically useful vs treat-all/none?
#   (3) OPTIMISM CORRECTION  - Harrell enhanced-bootstrap apparent->corrected.
#
# Design choice that makes this honest and defensible:
#   - The CLINICAL MODEL is a 1-feature logistic: y ~ composite_T0, where
#     composite_T0 = mean z of the PRE-SPECIFIED gate (Ki67 subsets from Step 04).
#     No in-fold feature selection => minimal optimism (the whole point: contrast
#     with the nested-LOOCV ML whose selection variance kills the signal).
#   - All "honest" curves use LEAKAGE-FREE LOO predicted probabilities.
#   - Optimism is quantified two independent ways (LOO and bootstrap) so we can
#     report a corrected AUC/slope with a CI.
#
# Runs on the discovery-73 longitudinal RDS (same input as multilevel_t0_utility).
# Writes JSON + PDFs to results/diag_t0_calibration_<ts>/.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(here); library(jsonlite); library(pROC); library(ggplot2)
})
options(crayon.enabled = FALSE)
set.seed(2026)

EXP    <- "BestResponse_2v3_4"
GATE   <- c("KI67NAIVE", "CD28KI67", "KI67CD4")   # LMM FDR-robust gate (Step 04)
N_BOOT <- 2000                                    # enhanced-bootstrap reps
N_PERM <- 2000

RDS_PATH <- Sys.getenv("RDS_PATH", unset = here(
  "results/latest", EXP, "01_data_processing",
  sprintf("data_processed_%s_longitudinal.rds", EXP)))
stopifnot(file.exists(RDS_PATH))

d       <- readRDS(RDS_PATH)
markers <- d$markers
meta    <- as.data.frame(d$metadata)
Z       <- as.data.frame(lapply(d$hybrid_data_z[, markers, drop = FALSE], as.numeric))
Z[]     <- lapply(Z, function(x) { x[is.na(x)] <- median(x, na.rm = TRUE); x })
GATE    <- intersect(GATE, markers)

out_dir <- here(paste0("results/diag_t0_calibration_", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- T0 slice + composite ---------------------------------------------------
t0   <- meta$Timepoint == "T0"
Zt0  <- Z[t0, , drop = FALSE]
y    <- factor(meta$Group[t0], levels = c("SD_PD", "RP"))   # RP = positive (event)
yb   <- as.integer(y == "RP")
comp <- rowMeans(Zt0[, GATE, drop = FALSE])
n    <- length(yb)
message(sprintf("[calib] RDS: %s", RDS_PATH))
message(sprintf("[calib] n=%d  (RP=%d, SD_PD=%d)  gate=%s", n, sum(yb), sum(1-yb),
                paste(GATE, collapse = "+")))

# ---- helpers ----------------------------------------------------------------
auc_pos <- function(yv, p) as.numeric(pROC::auc(pROC::roc(
  yv, p, levels = c("SD_PD", "RP"), direction = "<", quiet = TRUE)))
brier   <- function(yvb, p) mean((p - yvb)^2)
# calibration intercept (CITL) & slope by logistic recalibration on the linear predictor
calib_params <- function(yvb, p) {
  p   <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  lp  <- qlogis(p)
  slope <- tryCatch(coef(glm(yvb ~ lp, family = binomial()))[["lp"]], error = function(e) NA_real_)
  citl  <- tryCatch(coef(glm(yvb ~ 1, family = binomial(), offset = lp))[[1]], error = function(e) NA_real_)
  c(intercept = citl, slope = slope)
}

# ---- 1-feature logistic: probabilities (apparent + LOO) ---------------------
df_full <- data.frame(yb = yb, comp = comp)
fit_app <- glm(yb ~ comp, data = df_full, family = binomial())
p_app   <- as.numeric(predict(fit_app, type = "response"))

p_loo <- numeric(n)
for (i in seq_len(n)) {
  fi <- glm(yb ~ comp, data = df_full[-i, , drop = FALSE], family = binomial())
  p_loo[i] <- as.numeric(predict(fi, newdata = df_full[i, , drop = FALSE], type = "response"))
}

auc_app <- auc_pos(y, p_app);  auc_loo <- auc_pos(y, p_loo)
br_app  <- brier(yb, p_app);   br_loo  <- brier(yb, p_loo)
cp_app  <- calib_params(yb, p_app); cp_loo <- calib_params(yb, p_loo)
message(sprintf("[calib] apparent : AUC=%.3f Brier=%.3f  CITL=%.2f slope=%.2f",
                auc_app, br_app, cp_app["intercept"], cp_app["slope"]))
message(sprintf("[calib] LOO       : AUC=%.3f Brier=%.3f  CITL=%.2f slope=%.2f",
                auc_loo, br_loo, cp_loo["intercept"], cp_loo["slope"]))

# permutation p for LOO AUC (shuffle labels, full LOO refit)
perm_auc_loo <- function() {
  yp <- sample(yb); pl <- numeric(n)
  dfp <- data.frame(yb = yp, comp = comp)
  for (i in seq_len(n)) {
    fi <- suppressWarnings(glm(yb ~ comp, data = dfp[-i, ], family = binomial()))
    pl[i] <- as.numeric(predict(fi, newdata = dfp[i, ], type = "response"))
  }
  auc_pos(factor(ifelse(yp == 1, "RP", "SD_PD"), levels = c("SD_PD","RP")), pl)
}
cnt <- 1; for (b in seq_len(N_PERM)) if (perm_auc_loo() >= auc_loo) cnt <- cnt + 1
perm_p_loo <- cnt / (N_PERM + 1)
message(sprintf("[calib] LOO AUC perm p (n=%d) = %.4f", N_PERM, perm_p_loo))

# =============================================================================
# (3) HARRELL ENHANCED BOOTSTRAP - optimism of AUC & calibration slope
# =============================================================================
opt_auc <- opt_slope <- numeric(N_BOOT)
for (b in seq_len(N_BOOT)) {
  idx <- sample(n, n, replace = TRUE)
  dfb <- df_full[idx, , drop = FALSE]
  fb  <- suppressWarnings(glm(yb ~ comp, data = dfb, family = binomial()))
  # performance on bootstrap sample (boot) and on original sample (orig)
  pb_boot <- as.numeric(predict(fb, newdata = dfb, type = "response"))
  pb_orig <- as.numeric(predict(fb, newdata = df_full, type = "response"))
  yb_boot <- dfb$yb
  a_boot <- tryCatch(auc_pos(factor(ifelse(yb_boot==1,"RP","SD_PD"), levels=c("SD_PD","RP")), pb_boot),
                     error = function(e) NA_real_)
  a_orig <- auc_pos(y, pb_orig)
  s_boot <- calib_params(yb_boot, pb_boot)["slope"]
  s_orig <- calib_params(yb,      pb_orig)["slope"]
  opt_auc[b]   <- a_boot - a_orig
  opt_slope[b] <- s_boot - s_orig
}
opt_auc_m   <- mean(opt_auc,   na.rm = TRUE)
opt_slope_m <- mean(opt_slope, na.rm = TRUE)
auc_corr    <- auc_app - opt_auc_m
slope_corr  <- cp_app["slope"] - opt_slope_m
# bootstrap CI on the corrected AUC (apparent - per-rep optimism)
auc_corr_ci <- quantile(auc_app - opt_auc, c(.025, .975), na.rm = TRUE)
message(sprintf("[boot] optimism AUC=%.4f  -> corrected AUC=%.3f [%.3f-%.3f]",
                opt_auc_m, auc_corr, auc_corr_ci[1], auc_corr_ci[2]))
message(sprintf("[boot] optimism slope=%.4f -> corrected slope=%.3f",
                opt_slope_m, slope_corr))

# =============================================================================
# (2) DECISION CURVE ANALYSIS  (net benefit, on honest LOO probabilities)
# =============================================================================
prev <- mean(yb)
pt_grid <- seq(0.01, 0.60, by = 0.01)
nb <- map_dfr(pt_grid, function(pt) {
  pos  <- p_loo >= pt
  TP   <- sum(pos & yb == 1); FP <- sum(pos & yb == 0)
  w    <- pt / (1 - pt)
  nb_model <- TP/n - FP/n * w
  nb_all   <- prev - (1 - prev) * w     # treat-all
  tibble(threshold = pt, model = nb_model, treat_all = nb_all, treat_none = 0)
})
nb_long <- nb %>% pivot_longer(c(model, treat_all, treat_none),
                               names_to = "strategy", values_to = "net_benefit")
write_csv(nb, file.path(out_dir, "decision_curve.csv"))
# range of thresholds where the model dominates BOTH defaults
dom <- nb %>% filter(model > treat_all & model > treat_none)
dom_range <- if (nrow(dom)) sprintf("%.0f%%-%.0f%%", 100*min(dom$threshold), 100*max(dom$threshold)) else "none"
message(sprintf("[dca] model dominates treat-all & treat-none over pt in [%s]", dom_range))

# =============================================================================
# PLOTS
# =============================================================================
# calibration plot (LOO): loess smooth + quantile bins
calib_df <- tibble(p = p_loo, y = yb)
bins <- calib_df %>% mutate(bin = ntile(p, 4)) %>% group_by(bin) %>%
  summarise(p_mean = mean(p), obs = mean(y), nseg = n(),
            se = sqrt(obs*(1-obs)/nseg), .groups = "drop")
pcal <- ggplot(calib_df, aes(p, y)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_smooth(method = "loess", se = TRUE, colour = "#2166AC", fill = "#D1E5F0", span = 1) +
  geom_point(data = bins, aes(p_mean, obs), inherit.aes = FALSE, size = 2.5) +
  geom_errorbar(data = bins, aes(x = p_mean, ymin = pmax(0, obs-se), ymax = pmin(1, obs+se)),
                inherit.aes = FALSE, width = 0.01) +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  labs(title = "Calibration - T0 Ki67 composite (LOO)",
       subtitle = sprintf("CITL=%.2f  slope=%.2f  Brier=%.3f  (n=%d)",
                          cp_loo["intercept"], cp_loo["slope"], br_loo, n),
       x = "Predicted probability (RP)", y = "Observed fraction RP") + theme_bw()
ggsave(file.path(out_dir, "calibration_LOO.pdf"), pcal, width = 6, height = 6)

pdca <- ggplot(nb_long, aes(threshold, net_benefit, colour = strategy)) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c(model = "#2166AC", treat_all = "#B2182B", treat_none = "grey40")) +
  coord_cartesian(ylim = c(min(-0.02, min(nb$treat_all)), max(nb$model) + 0.02)) +
  labs(title = "Decision curve - T0 Ki67 composite (LOO probabilities)",
       subtitle = sprintf("model dominates defaults over pt in [%s]; prevalence=%.2f", dom_range, prev),
       x = "Threshold probability", y = "Net benefit", colour = NULL) +
  theme_bw() + theme(legend.position = "top")
ggsave(file.path(out_dir, "decision_curve.pdf"), pdca, width = 6.5, height = 5.5)

# =============================================================================
# SUMMARY JSON + verdict
# =============================================================================
summary <- list(
  rds = RDS_PATH, n = n, n_RP = sum(yb), n_SDPD = sum(1-yb), gate = GATE, prevalence = prev,
  discrimination = list(auc_apparent = auc_app, auc_loo = auc_loo, loo_perm_p = perm_p_loo,
                        auc_optimism = opt_auc_m, auc_corrected = auc_corr,
                        auc_corrected_ci = as.numeric(auc_corr_ci)),
  calibration = list(
    apparent = list(intercept = cp_app[["intercept"]], slope = cp_app[["slope"]], brier = br_app),
    loo      = list(intercept = cp_loo[["intercept"]], slope = cp_loo[["slope"]], brier = br_loo),
    slope_optimism = opt_slope_m, slope_corrected = as.numeric(slope_corr)),
  decision_curve = list(prevalence = prev, dominance_range = dom_range,
                        nb_model_at_prev = nb$model[which.min(abs(pt_grid - prev))]),
  interpretation = "Pre-specified 1-feature logistic on Ki67 composite: minimal optimism by design, in contrast to nested-LOOCV ML (AUC 0.56)."
)
write_json(summary, file.path(out_dir, "summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = 4)

message("\n=================== VERDICT ===================")
message(sprintf("AUC apparent=%.3f | LOO=%.3f (perm p=%.4f) | boot-corrected=%.3f [%.3f-%.3f]",
                auc_app, auc_loo, perm_p_loo, auc_corr, auc_corr_ci[1], auc_corr_ci[2]))
message(sprintf("Calibration LOO: CITL=%.2f slope=%.2f Brier=%.3f", cp_loo["intercept"], cp_loo["slope"], br_loo))
message(sprintf("Optimism: AUC=%.4f slope=%.4f (small => pre-specified score is stable)", opt_auc_m, opt_slope_m))
message(sprintf("DCA: net benefit > treat-all/none over pt in [%s]", dom_range))
message(sprintf("\nOutputs -> %s", out_dir))
