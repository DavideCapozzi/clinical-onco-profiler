# R/modules_longitudinal.R
# ==============================================================================
# LONGITUDINAL LMM MODULE
# Description: Wrapper functions for Linear Mixed Models handling time & group.
# Dependencies: lmerTest, dplyr, ggplot2, ggrepel
# ==============================================================================

library(dplyr)
library(ggplot2)
library(ggrepel)

#' @title Marginal and Conditional R² for LMM (Nakagawa & Schielzeth 2013)
#' @description Computes variance-partition R² without requiring external packages.
#' @param mod A fitted lmerMod object.
#' @return Named list with R2m (marginal, fixed effects only) and R2c (conditional).
r2_nakagawa <- function(mod) {
  tryCatch({
    sigma2_f <- var(predict(mod, re.form = NA))
    vc       <- as.data.frame(lme4::VarCorr(mod))
    sigma2_r <- sum(vc$vcov[vc$grp != "Residual"])
    sigma2_e <- sigma(mod)^2
    total    <- sigma2_f + sigma2_r + sigma2_e
    list(R2m = round(sigma2_f / total, 4),
         R2c = round((sigma2_f + sigma2_r) / total, 4))
  }, error = function(e) list(R2m = NA_real_, R2c = NA_real_))
}

#' @title Run Linear Mixed Model on a single feature
#' @description Fits LMM with an interaction term between Time and Clinical Group.
#' Includes internal standard deviation scaling and Marginal Time Effect extraction.
#' @param data_long Dataframe in long format.
#' @param feature String. Name of the column containing feature values.
#' @param group_col String. Name of the column containing clinical labels.
#' @param time_col String. Name of the column containing timepoints.
#' @param id_col String. Name of the column containing Patient IDs.
#' @param covariates Vector of strings. Names of clinical covariates to include.
#' @return A one-row dataframe with model statistics.
fit_feature_lmm <- function(data_long, feature, group_col = "Group", 
                            time_col = "Timepoint", id_col = "Patient_ID",
                            covariates = NULL) {
  
  if (!requireNamespace("lmerTest", quietly = TRUE)) {
    stop("Package 'lmerTest' is required for p-value calculation in LMM.")
  }
  
  raw_vals <- as.numeric(data_long[[feature]])
  val_sd <- sd(raw_vals, na.rm = TRUE)
  if (is.na(val_sd) || val_sd == 0) val_sd <- 1 # Safety fallback for zero variance
  
  df_list <- list(
    Value = raw_vals / val_sd, # Internal scaling for optimizer convergence
    Time = as.factor(data_long[[time_col]]),
    Group = as.factor(data_long[[group_col]]),
    ID = as.factor(data_long[[id_col]])
  )
  
  valid_covs <- c()
  if (!is.null(covariates) && length(covariates) > 0) {
    valid_covs <- intersect(covariates, colnames(data_long))
    for (cov in valid_covs) {
      df_list[[cov]] <- data_long[[cov]]
    }
  }
  
  df_model <- as.data.frame(df_list, check.names = FALSE)
  df_model <- df_model[complete.cases(df_model), ]
  
  result <- data.frame(
    Marker = feature,
    Estimate_Interaction = NA,
    T_Value_Interaction = NA,
    Std_Error = NA,
    P_Value_Interaction = NA,
    Estimate_Time_Main = NA,
    P_Value_Time_Main = NA,
    Model_Converged = FALSE,
    Is_Singular = NA,
    N_Observations = nrow(df_model),
    R2m = NA_real_,
    R2c = NA_real_
  )
  
  if (nrow(df_model) < 10) return(result)
  
  tryCatch({
    # --- MODEL 1: Interaction Model (Primary Objective) ---
    formula_int <- "Value ~ Time * Group"
    if (length(valid_covs) > 0) formula_int <- paste(formula_int, "+", paste(sprintf("`%s`", valid_covs), collapse = " + "))
    formula_int <- paste(formula_int, "+ (1 | ID)")
    
    mod_int <- suppressMessages(suppressWarnings(
      lmerTest::lmer(as.formula(formula_int), data = df_model, 
                     REML = TRUE, control = lme4::lmerControl(calc.derivs = FALSE))
    ))
    
    result$Is_Singular <- lme4::isSingular(mod_int)
    
    coef_table_int <- summary(mod_int)$coefficients
    interaction_idx <- grep("Time.*:Group", rownames(coef_table_int))
    
    if (length(interaction_idx) == 1) {
      # Reverse scaling to restore native hybrid logit/log2 effect sizes
      result$Estimate_Interaction <- coef_table_int[interaction_idx, "Estimate"] * val_sd
      result$Std_Error <- coef_table_int[interaction_idx, "Std. Error"] * val_sd
      
      t_col <- grep("t value", colnames(coef_table_int))
      if (length(t_col) == 1) result$T_Value_Interaction <- coef_table_int[interaction_idx, t_col]
      
      p_col <- grep("Pr\\(>\\|t\\|\\)", colnames(coef_table_int))
      if (length(p_col) == 1) {
        result$P_Value_Interaction <- coef_table_int[interaction_idx, p_col]
        result$Model_Converged <- TRUE
      }
    }

    r2 <- r2_nakagawa(mod_int)
    result$R2m <- r2$R2m
    result$R2c <- r2$R2c
    
    # --- MODEL 2: Marginal Time Model (Positive Control) ---
    formula_time <- "Value ~ Time"
    if (length(valid_covs) > 0) formula_time <- paste(formula_time, "+", paste(sprintf("`%s`", valid_covs), collapse = " + "))
    formula_time <- paste(formula_time, "+ (1 | ID)")
    
    mod_time <- suppressMessages(suppressWarnings(
      lmerTest::lmer(as.formula(formula_time), data = df_model, 
                     REML = TRUE, control = lme4::lmerControl(calc.derivs = FALSE))
    ))
    
    coef_table_time <- summary(mod_time)$coefficients
    time_idx <- grep("Time", rownames(coef_table_time))
    
    if (length(time_idx) == 1) {
      result$Estimate_Time_Main <- coef_table_time[time_idx, "Estimate"]
      p_col <- grep("Pr\\(>\\|t\\|\\)", colnames(coef_table_time))
      if (length(p_col) == 1) result$P_Value_Time_Main <- coef_table_time[time_idx, p_col]
    }
    
  }, error = function(e) {
    # Silent fail on calculation error, returning NA defaults
  })
  
  return(result)
}

#' @title Run Leave-One-Out (LOO) Sensitivity Analysis
#' @description Iteratively drops one patient at a time to test interaction robustness.
#' @return Numeric. The maximum (worst-case) interaction P-value found across all iterations.
run_loo_sensitivity <- function(data_long, feature, group_col = "Group", 
                                time_col = "Timepoint", id_col = "Patient_ID",
                                covariates = NULL) {
  
  unique_ids <- unique(data_long[[id_col]])
  max_p_val <- 0
  
  for (drop_id in unique_ids) {
    # Create LOO dataset
    df_subset <- data_long[data_long[[id_col]] != drop_id, ]
    
    # Fit model silently
    res <- fit_feature_lmm(data_long = df_subset, feature = feature, 
                           group_col = group_col, time_col = time_col, 
                           id_col = id_col, covariates = covariates)
    
    # Track the worst p-value
    current_p <- res$P_Value_Interaction
    if (!is.na(current_p) && current_p > max_p_val) {
      max_p_val <- current_p
    }
  }
  
  # max_p_val stays 0 only when every LOO fold returned NA (all models failed).
  # Return NA so the gate filter (!is.na) correctly excludes the marker.
  return(if (max_p_val == 0) NA else max_p_val)
}

#' @title Split-Group LMM Supplementary Analysis
#' @description
#' Supplementary analysis for FDR-significant markers when the non-responder
#' label collapses multiple raw clinical codes (e.g. "SD_PD" = {SD=3, PD=4}).
#' Refits the LMM with the non-responder split back into its sub-levels,
#' keeping the responder side as a single class. The reference is the FIRST
#' sub-level (so betas are vs the most adverse outcome).
#'
#' Why this is worth running: the dichotomization RP vs (SD union PD) trades
#' biological texture for power. If responders show one direction and SD vs
#' PD show different directions, the dichotomized interaction is dominated
#' by the responder contrast — a 3-way diagnostic clarifies which sub-group
#' is actually driving the LMM hit.
#'
#' @param data_long Long-format dataframe with Patient_ID, Timepoint, Group,
#'   and the marker columns.
#' @param features Character vector of marker names to test.
#' @param patient_subgroup Named character vector mapping Patient_ID -> raw
#'   sub-level label (e.g. "SD"/"PD"/"RP"). Patients absent from this vector
#'   are excluded.
#' @param ref_level Character. Reference level for the 3-way Group factor
#'   (defaults to the first observed level).
#' @return Data.frame with one row per (marker x non-reference level):
#'   Marker, Level, Estimate_Interaction, Std_Error, T_Value, P_Value,
#'   N_Observations, Is_Singular. NULL on failure.
run_splitgroup_lmm <- function(data_long, features, patient_subgroup,
                               ref_level = NULL) {
  if (!requireNamespace("lmerTest", quietly = TRUE)) {
    warning("[Split-Group LMM] lmerTest not installed; skipping.")
    return(NULL)
  }
  if (length(features) == 0 || length(patient_subgroup) == 0) return(NULL)

  d <- data_long
  d$GroupSplit <- unname(patient_subgroup[as.character(d$Patient_ID)])
  d <- d[!is.na(d$GroupSplit), , drop = FALSE]
  if (nrow(d) == 0) return(NULL)

  levels_obs <- unique(d$GroupSplit)
  if (is.null(ref_level) || !(ref_level %in% levels_obs)) {
    ref_level <- levels_obs[1]
  }
  d$GroupSplit <- factor(d$GroupSplit,
                         levels = c(ref_level, setdiff(levels_obs, ref_level)))

  rows <- list()
  for (mk in features) {
    sub <- d[, c("Patient_ID", "Timepoint", "GroupSplit", mk), drop = FALSE]
    names(sub) <- c("ID", "Time", "Group", "Value")
    sub <- sub[complete.cases(sub), ]
    if (nrow(sub) < 10) next

    val_sd <- sd(sub$Value, na.rm = TRUE)
    if (is.na(val_sd) || val_sd == 0) val_sd <- 1
    sub$Value <- sub$Value / val_sd
    sub$Time  <- factor(sub$Time, levels = sort(unique(as.character(sub$Time))))

    fit <- tryCatch(
      suppressMessages(suppressWarnings(
        lmerTest::lmer(Value ~ Time * Group + (1 | ID), data = sub,
                       REML = TRUE,
                       control = lme4::lmerControl(calc.derivs = FALSE))
      )),
      error = function(e) NULL
    )
    if (is.null(fit)) next

    co <- summary(fit)$coefficients
    is_sing <- isTRUE(lme4::isSingular(fit))
    int_rows <- grep("^Time.*:Group", rownames(co), value = TRUE)
    for (rn in int_rows) {
      lvl <- sub("^Time.*:Group", "", rn)
      r <- co[rn, , drop = FALSE]
      rows[[length(rows) + 1]] <- data.frame(
        Marker               = mk,
        Level                = lvl,
        Reference            = ref_level,
        Estimate_Interaction = round(unname(r[1, "Estimate"]) * val_sd, 4),
        Std_Error            = round(unname(r[1, "Std. Error"]) * val_sd, 4),
        T_Value              = round(unname(r[1, "t value"]), 4),
        P_Value              = round(unname(r[1, "Pr(>|t|)"]), 6),
        N_Observations       = nrow(sub),
        Is_Singular          = is_sing,
        stringsAsFactors     = FALSE
      )
    }
  }
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}


#' @title Paired-Only Delta Sensitivity Analysis
#' @description
#' Sensitivity check for FDR-significant LMM interaction findings. Restricts
#' the cohort to patients with BOTH T0 and T1 measurements for the marker and
#' fits a plain OLS on the within-patient delta:
#'   lm(delta ~ Group), where delta_i = value_i(T1) - value_i(T0)
#'
#' This sidesteps the random-intercept variance estimation that drives LMM
#' singular fits when many patients contribute only one timepoint. Under
#' balanced paired data the OLS-on-delta slope equals the LMM Time x Group
#' interaction term, so agreement between the two estimators is the relevant
#' robustness signal — divergence flags an identification problem in the LMM.
#'
#' Effect sizes are returned in the native (un-scaled) hybrid scale of the
#' marker, matching the LMM Estimate_Interaction column.
#'
#' @param data_long Dataframe in long format with at least Patient_ID,
#'   Timepoint (with levels including "T0" and "T1"), Group, and the feature
#'   columns.
#' @param features Character vector of marker names to test.
#' @param group_col,time_col,id_col Column names in data_long.
#' @return Data.frame with one row per marker: Marker, Estimate_Delta,
#'   Std_Error, T_Value, P_Value, FDR, N_Pairs, plus the two timepoint
#'   contributions used (n_only_T0, n_only_T1) for transparency.
run_paired_only_sensitivity <- function(data_long, features,
                                        group_col = "Group",
                                        time_col  = "Timepoint",
                                        id_col    = "Patient_ID") {
  if (length(features) == 0) return(NULL)

  pid_vec  <- as.character(data_long[[id_col]])
  time_vec <- as.character(data_long[[time_col]])

  tp_per_pid <- split(time_vec, pid_vec)
  has_t0 <- vapply(tp_per_pid, function(v) "T0" %in% v, logical(1))
  has_t1 <- vapply(tp_per_pid, function(v) "T1" %in% v, logical(1))
  paired_pids <- names(tp_per_pid)[has_t0 & has_t1]

  n_only_t0 <- sum(has_t0 & !has_t1)
  n_only_t1 <- sum(has_t1 & !has_t0)

  if (length(paired_pids) < 5) {
    warning(sprintf("[Paired Sensitivity] Only %d paired patient(s); skipping.",
                    length(paired_pids)))
    return(NULL)
  }

  d_paired <- data_long[pid_vec %in% paired_pids, , drop = FALSE]

  rows <- lapply(features, function(mk) {
    sub <- d_paired[, c(id_col, time_col, group_col, mk), drop = FALSE]
    names(sub) <- c("ID", "Time", "Group", "Value")
    sub <- sub[complete.cases(sub), ]
    if (nrow(sub) < 2) {
      return(data.frame(Marker = mk, Estimate_Delta = NA_real_, Std_Error = NA_real_,
                        T_Value = NA_real_, P_Value = NA_real_, N_Pairs = 0L,
                        stringsAsFactors = FALSE))
    }
    # Pivot to wide. tapply over a factor ID column returns all factor levels
    # (NA for those without rows in this stratum), which would over-report
    # pairing — drop NAs before intersecting.
    sub$ID <- as.character(sub$ID)
    val_t0 <- tapply(sub$Value[sub$Time == "T0"], sub$ID[sub$Time == "T0"], mean)
    val_t1 <- tapply(sub$Value[sub$Time == "T1"], sub$ID[sub$Time == "T1"], mean)
    val_t0 <- val_t0[!is.na(val_t0)]
    val_t1 <- val_t1[!is.na(val_t1)]
    grp    <- tapply(as.character(sub$Group), sub$ID, function(g) g[1])
    ids    <- intersect(names(val_t0), names(val_t1))
    if (length(ids) < 5) {
      return(data.frame(Marker = mk, Estimate_Delta = NA_real_, Std_Error = NA_real_,
                        T_Value = NA_real_, P_Value = NA_real_,
                        N_Pairs = length(ids), stringsAsFactors = FALSE))
    }
    delta <- val_t1[ids] - val_t0[ids]
    # Preserve the level ordering of the input Group factor (which has the
    # non-responder label set as the reference in Step 04). With reference =
    # non-responder, the slope here matches the sign convention of the LMM
    # Time:Group interaction term — i.e. negative slope = responders contract
    # more than non-responders.
    if (is.factor(data_long[[group_col]])) {
      g <- factor(grp[ids], levels = levels(data_long[[group_col]]))
    } else {
      g <- factor(grp[ids])
    }
    fit <- tryCatch(lm(delta ~ g), error = function(e) NULL)
    if (is.null(fit) || nrow(summary(fit)$coefficients) < 2) {
      return(data.frame(Marker = mk, Estimate_Delta = NA_real_, Std_Error = NA_real_,
                        T_Value = NA_real_, P_Value = NA_real_,
                        N_Pairs = length(ids), stringsAsFactors = FALSE))
    }
    co <- summary(fit)$coefficients[2, , drop = FALSE]
    data.frame(
      Marker         = mk,
      Estimate_Delta = round(co[1, "Estimate"], 4),
      Std_Error      = round(co[1, "Std. Error"], 4),
      T_Value        = round(co[1, "t value"], 4),
      P_Value        = round(co[1, "Pr(>|t|)"], 6),
      N_Pairs        = length(ids),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$FDR <- p.adjust(out$P_Value, method = "BH")
  attr(out, "n_paired_total") <- length(paired_pids)
  attr(out, "n_only_T0")      <- as.integer(n_only_t0)
  attr(out, "n_only_T1")      <- as.integer(n_only_t1)
  out
}


#' @title Rank ANCOVA Confirmatory Robustness Test ("LMM selects -> rank ANCOVA confirms")
#' @description
#' Distribution-free, baseline-adjusted, transform-invariant confirmation of the
#' FDR-significant LMM interaction findings. The LMM remains the gate SELECTOR
#' (published-dynamics provenance + nested-fold validation + uses all patients);
#' this is the CONFIRMATORY test, NOT a re-selection. For each target marker, on
#' the paired (T0+T1) subset, it computes:
#'   - Quade rank ANCOVA: lm(rank(T1) ~ rank(T0) + Group), group-term estimate/p
#'   - permutation p (patient-level Group label swap) for exact small-n inference
#'   - assumption-justification metrics documenting WHY a robust test is used:
#'       * LMM residual Shapiro-Wilk p (normality)
#'       * residual heteroscedasticity p (lm(resid^2 ~ fitted))
#'       * ranova random-intercept LRT p (is the random effect even supported?)
#'       * max Cook's distance from the parametric ANCOVA (influential observation)
#'
#' Rank-based => invariant to the logit/log2 transform by construction;
#' baseline-adjusted (rank(T0)) => robust to baseline imbalance / regression-to-
#' the-mean that biases change-scores. RNG is self-contained (seed saved/restored)
#' so primary pipeline results reproduce byte-for-byte.
#'
#' @param data_long Long-format df (Patient_ID, Timepoint with T0/T1, Group, markers).
#' @param features Character vector of markers to confirm (FDR-significant subset).
#' @param group_col,time_col,id_col Column names in data_long. Group must be a
#'   factor with the non-responder label as its reference (matching Step 04), so
#'   the group term sign matches the LMM interaction convention.
#' @param n_perm Permutation iterations (default 2000).
#' @param seed RNG seed (default 2026).
#' @return Data.frame, one row per marker (Marker, N_Pairs, RankANCOVA_Estimate,
#'   RankANCOVA_P, RankANCOVA_Perm_P, Resid_Shapiro_P, Hetero_P, Ranova_RE_P,
#'   Max_Cooks_D), with n_perm/seed attributes. NULL if no marker has paired data.
run_rank_ancova_confirmation <- function(data_long, features,
                                         group_col = "Group", time_col = "Timepoint",
                                         id_col = "Patient_ID", n_perm = 2000L, seed = 2026L) {
  if (length(features) == 0) return(NULL)

  # Self-contained RNG: restore the global stream on exit so nothing downstream shifts.
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  }
  set.seed(seed)

  glev <- levels(as.factor(data_long[[group_col]]))   # reference (non-responder) first

  build_wide <- function(mk) {
    sub <- data.frame(ID    = as.character(data_long[[id_col]]),
                      Time  = as.character(data_long[[time_col]]),
                      Group = as.character(data_long[[group_col]]),
                      Value = as.numeric(data_long[[mk]]), stringsAsFactors = FALSE)
    sub <- sub[complete.cases(sub), ]
    t0 <- tapply(sub$Value[sub$Time == "T0"], sub$ID[sub$Time == "T0"], mean)
    t1 <- tapply(sub$Value[sub$Time == "T1"], sub$ID[sub$Time == "T1"], mean)
    t0 <- t0[!is.na(t0)]; t1 <- t1[!is.na(t1)]
    grp <- tapply(sub$Group, sub$ID, function(g) g[1])
    ids <- intersect(names(t0), names(t1))
    if (length(ids) < 5) return(NULL)
    data.frame(T0 = t0[ids], T1 = t1[ids],
               g  = factor(grp[ids], levels = glev),
               rT0 = rank(t0[ids]), rT1 = rank(t1[ids]))
  }
  grp_term <- function(fit, what) {
    co <- summary(fit)$coefficients
    i  <- grep("^g", rownames(co))
    if (length(i) != 1) return(NA_real_)
    if (what == "est") return(unname(co[i, "Estimate"]))
    pc <- grep("Pr\\(>\\|t\\|\\)", colnames(co))
    if (length(pc) != 1) return(NA_real_)
    unname(co[i, pc])
  }

  rows <- lapply(features, function(mk) {
    na_row <- data.frame(Marker = mk, N_Pairs = 0L, RankANCOVA_Estimate = NA_real_,
                         RankANCOVA_P = NA_real_, RankANCOVA_Perm_P = NA_real_,
                         Resid_Shapiro_P = NA_real_, Hetero_P = NA_real_,
                         Ranova_RE_P = NA_real_, Max_Cooks_D = NA_real_,
                         stringsAsFactors = FALSE)
    w <- build_wide(mk)
    if (is.null(w)) return(na_row)

    f_rank <- lm(rT1 ~ rT0 + g, data = w)
    obs    <- grp_term(f_rank, "est")
    null   <- replicate(n_perm, {
      w2 <- w; w2$g <- factor(sample(as.character(w$g)), levels = glev)
      grp_term(lm(rT1 ~ rT0 + g, data = w2), "est")
    })
    perm_p <- (sum(abs(null) >= abs(obs), na.rm = TRUE) + 1) / (sum(!is.na(null)) + 1)

    cookmax <- tryCatch(suppressWarnings(max(cooks.distance(lm(T1 ~ T0 + g, data = w)), na.rm = TRUE)),
                        error = function(e) NA_real_)

    # Refit the LMM only to extract residual / random-effect justification metrics
    # (the primary pipeline discards the model object).
    sh <- NA_real_; het <- NA_real_; reP <- NA_real_
    md <- data.frame(Value = as.numeric(data_long[[mk]]),
                     Time  = as.factor(data_long[[time_col]]),
                     Group = as.factor(data_long[[group_col]]),
                     ID    = as.factor(data_long[[id_col]]))
    md <- md[complete.cases(md), ]
    s <- sd(md$Value); if (is.na(s) || s == 0) s <- 1; md$Value <- md$Value / s
    m <- tryCatch(suppressWarnings(suppressMessages(
           lmerTest::lmer(Value ~ Time * Group + (1 | ID), data = md, REML = TRUE,
                          control = lme4::lmerControl(calc.derivs = FALSE)))),
         error = function(e) NULL)
    if (!is.null(m)) {
      r    <- tryCatch(resid(m, type = "pearson"), error = function(e) NULL)
      fitv <- tryCatch(fitted(m), error = function(e) NULL)
      if (!is.null(r))               sh  <- tryCatch(shapiro.test(r)$p.value, error = function(e) NA_real_)
      if (!is.null(r) && !is.null(fitv)) het <- tryCatch(summary(lm(I(r^2) ~ fitv))$coefficients[2, 4], error = function(e) NA_real_)
      rv <- tryCatch(lmerTest::ranova(m), error = function(e) NULL)
      if (!is.null(rv) && "Pr(>Chisq)" %in% names(rv)) reP <- rv$`Pr(>Chisq)`[2]
    }

    data.frame(Marker = mk, N_Pairs = nrow(w),
               RankANCOVA_Estimate = round(obs, 4),
               RankANCOVA_P        = signif(grp_term(f_rank, "p"), 4),
               RankANCOVA_Perm_P   = signif(perm_p, 4),
               Resid_Shapiro_P     = signif(sh, 4),
               Hetero_P            = signif(het, 4),
               Ranova_RE_P         = signif(reP, 4),
               Max_Cooks_D         = round(cookmax, 3),
               stringsAsFactors    = FALSE)
  })
  out <- do.call(rbind, rows)
  attr(out, "n_perm") <- as.integer(n_perm)
  attr(out, "seed")   <- as.integer(seed)
  out
}


#' @title Cluster Bootstrap CI for LMM Interaction Betas
#' @description
#' Patient-level cluster bootstrap (resample patient IDs with replacement,
#' preserving the within-patient T0/T1 pairing). For each bootstrap iteration
#' refits the full LMM panel and re-applies BH-FDR, then reports for the
#' supplied target markers:
#'   - bootstrap median estimate
#'   - 2.5% / 97.5% bootstrap quantiles (95% CI)
#'   - two-sided proportion-based bootstrap p-value
#'   - fraction of resamples in which the marker still passes FDR < fdr_threshold
#'
#' @param data_long Dataframe in long format (one row per (Patient_ID, Timepoint)).
#' @param all_markers Character vector of all markers fitted in the primary LMM
#'   panel (needed for BH-FDR re-computation inside each bootstrap iteration).
#' @param target_markers Character vector of markers to report results for
#'   (typically the FDR-significant subset from the primary run).
#' @param group_col,time_col,id_col Column names in data_long.
#' @param covariates Optional clinical covariates passed to fit_feature_lmm.
#' @param n_boot Integer. Number of bootstrap iterations (default 500).
#' @param fdr_threshold Numeric. FDR cutoff for the per-iteration significance
#'   tally (default 0.05).
#' @param seed Integer. RNG seed for reproducibility (default 2026).
#' @param progress_message Logical. Print every 100 iterations (default TRUE).
#' @return A list with:
#'   - summary_df: one row per target marker (Marker, Median_Beta_Boot,
#'                 CI_Lower_2.5, CI_Upper_97.5, Bootstrap_P, Pct_FDR_Significant,
#'                 N_Valid_Iterations)
#'   - beta_matrix: B x length(target_markers) matrix of bootstrap betas
#'   - fdr_matrix:  B x length(target_markers) matrix of bootstrap FDR values
#'   - n_boot, seed, fdr_threshold
run_lmm_bootstrap_ci <- function(data_long,
                                 all_markers,
                                 target_markers,
                                 group_col       = "Group",
                                 time_col        = "Timepoint",
                                 id_col          = "Patient_ID",
                                 covariates      = NULL,
                                 n_boot          = 500L,
                                 fdr_threshold   = 0.05,
                                 seed            = 2026L,
                                 progress_message = TRUE) {

  target_markers <- intersect(target_markers, all_markers)
  if (length(target_markers) == 0) {
    return(list(summary_df = data.frame(), beta_matrix = NULL, fdr_matrix = NULL,
                n_boot = 0L, seed = seed, fdr_threshold = fdr_threshold))
  }

  unique_pids <- unique(data_long[[id_col]])
  n_pat       <- length(unique_pids)

  beta_mat <- matrix(NA_real_, nrow = n_boot, ncol = length(target_markers),
                     dimnames = list(NULL, target_markers))
  fdr_mat  <- matrix(NA_real_, nrow = n_boot, ncol = length(target_markers),
                     dimnames = list(NULL, target_markers))

  set.seed(seed)
  t0 <- Sys.time()

  for (b in seq_len(n_boot)) {
    sampled_pids <- sample(unique_pids, n_pat, replace = TRUE)

    rows_list <- lapply(seq_along(sampled_pids), function(j) {
      pid <- sampled_pids[j]
      sub <- data_long[data_long[[id_col]] == pid, , drop = FALSE]
      sub[[id_col]] <- paste0(pid, "_b", j)  # unique pseudo-ID per draw
      sub
    })
    data_boot <- do.call(rbind, rows_list)

    res_list <- lapply(all_markers, function(mk) {
      fit_feature_lmm(data_long = data_boot, feature = mk,
                      group_col = group_col, time_col = time_col,
                      id_col = id_col, covariates = covariates)
    })
    res_df <- do.call(rbind, res_list)
    res_df$FDR <- p.adjust(res_df$P_Value_Interaction, method = "BH")

    for (m in target_markers) {
      row <- res_df[res_df$Marker == m, , drop = FALSE]
      if (nrow(row) == 1 && !is.na(row$Estimate_Interaction)) {
        beta_mat[b, m] <- row$Estimate_Interaction
        fdr_mat[b, m]  <- row$FDR
      }
    }

    if (progress_message && b %% 100 == 0) {
      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      eta     <- elapsed / b * (n_boot - b)
      message(sprintf("      [Bootstrap] iter %d/%d  [%.0fs elapsed, ~%.0fs remaining]",
                      b, n_boot, elapsed, eta))
    }
  }

  summary_rows <- lapply(target_markers, function(m) {
    b_vec <- beta_mat[, m]; b_vec <- b_vec[!is.na(b_vec)]
    f_vec <- fdr_mat[, m];  f_vec <- f_vec[!is.na(f_vec)]
    if (length(b_vec) == 0) {
      return(data.frame(Marker = m, Median_Beta_Boot = NA_real_,
                        CI_Lower_2.5 = NA_real_, CI_Upper_97.5 = NA_real_,
                        Bootstrap_P = NA_real_, Pct_FDR_Significant = NA_real_,
                        N_Valid_Iterations = 0L, stringsAsFactors = FALSE))
    }
    ci      <- quantile(b_vec, c(0.025, 0.975))
    one_sided <- mean(b_vec >= 0)
    two_sided <- 2 * min(one_sided, 1 - one_sided)
    pct_fdr <- mean(f_vec < fdr_threshold) * 100
    data.frame(
      Marker             = m,
      Median_Beta_Boot   = round(median(b_vec), 4),
      CI_Lower_2.5       = round(ci[1], 4),
      CI_Upper_97.5      = round(ci[2], 4),
      Bootstrap_P        = round(two_sided, 4),
      Pct_FDR_Significant = round(pct_fdr, 1),
      N_Valid_Iterations = length(b_vec),
      stringsAsFactors   = FALSE
    )
  })
  summary_df <- do.call(rbind, summary_rows)

  list(
    summary_df    = summary_df,
    beta_matrix   = beta_mat,
    fdr_matrix    = fdr_mat,
    n_boot        = n_boot,
    seed          = seed,
    fdr_threshold = fdr_threshold
  )
}


#' @title Plot Longitudinal Trajectories (Spaghetti + Boxplot)
#' @description Visualizes patient trajectories over time, split by group.
#' @param data_long Dataframe in long format.
#' @param feature String. Marker name.
#' @param group_col String.
#' @param time_col String.
#' @param id_col String.
#' @param colors Named vector of colors.
#' @param p_val Optional numeric. P-value or FDR to display in subtitle.
#' @return A ggplot object.
plot_lmm_trajectories <- function(data_long, feature, group_col = "Group", 
                                  time_col = "Timepoint", id_col = "Patient_ID",
                                  colors = NULL, p_val = NULL) {
  
  require(ggplot2)
  
  # Ensure NAs are dropped for pure visualization
  plot_df <- data_long[!is.na(data_long[[feature]]), ]
  
  sub_title <- "Patient trajectories over time"
  if (!is.null(p_val)) {
    metric_name <- if (p_val < 0.05) "FDR" else "Interaction P-Value"
    sub_title <- sprintf("%s: %.4f", metric_name, p_val)
  }
  
  p <- ggplot(plot_df, aes(x = .data[[time_col]], y = .data[[feature]], fill = .data[[group_col]])) +
    
    # Base Boxplot distribution
    geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.4) +
    
    # Spaghetti Lines (Connecting individual patients)
    geom_line(aes(group = .data[[id_col]], color = .data[[group_col]]), alpha = 0.3, linewidth = 0.6) +
    
    # Individual Points
    geom_point(aes(fill = .data[[group_col]]), shape = 21, size = 2.5, color = "white", stroke = 0.3) +
    
    # Split by clinical group
    facet_wrap(as.formula(paste("~", group_col))) +
    
    # Aesthetics scaling
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    
    labs(
      title = paste("Trajectory:", feature),
      subtitle = sub_title,
      x = "Timepoint", 
      y = "Expression Level (Hybrid Scale)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "gray95"),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40")
    )
  
  return(p)
}

#' @title Volcano Plot for LMM Results
#' @description Plots the t-statistic vs -log10(FDR) of the interaction terms.
#' @param results_df Dataframe output from LMM loop.
#' @param title Plot title.
#' @return A ggplot object.
plot_lmm_volcano <- function(results_df, title = "Longitudinal Interaction (LMM)") {
  
  plot_df <- results_df %>% filter(!is.na(P_Value_Interaction), !is.na(T_Value_Interaction))
  if (nrow(plot_df) == 0) return(NULL)
  
  plot_df$logP <- -log10(plot_df$P_Value_Interaction)
  plot_df$logFDR <- -log10(plot_df$FDR_Interaction)
  
  plot_df$Significance <- "Not Significant"
  plot_df$Significance[plot_df$FDR_Interaction < 0.05 & plot_df$T_Value_Interaction > 0] <- "Positive Interaction"
  plot_df$Significance[plot_df$FDR_Interaction < 0.05 & plot_df$T_Value_Interaction < 0] <- "Negative Interaction"
  
  if (all(plot_df$Significance == "Not Significant")) {
    plot_df$Significance[plot_df$P_Value_Interaction < 0.05 & plot_df$T_Value_Interaction > 0] <- "Positive (Raw P < 0.05)"
    plot_df$Significance[plot_df$P_Value_Interaction < 0.05 & plot_df$T_Value_Interaction < 0] <- "Negative (Raw P < 0.05)"
  }
  
  colors <- c(
    "Not Significant" = "gray70",
    "Positive Interaction" = "#B2182B",
    "Negative Interaction" = "#2166AC",
    "Positive (Raw P < 0.05)" = "#F4A582",
    "Negative (Raw P < 0.05)" = "#92C5DE"
  )
  
  p <- ggplot(plot_df, aes(x = T_Value_Interaction, y = logFDR, color = Significance)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkred") +
    geom_point(size = 3, alpha = 0.8) +
    geom_text_repel(
      data = subset(plot_df, Significance != "Not Significant" | P_Value_Interaction < 0.01),
      aes(label = Marker), size = 3, show.legend = FALSE, max.overlaps = 20
    ) +
    scale_color_manual(values = colors) +
    labs(
      title = title,
      subtitle = "X: t-statistic of Time:Group Interaction | Y: -log10(FDR)",
      x = "LMM Interaction (t-statistic)",
      y = "-log10(FDR)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )

  return(p)
}


#' @title Forest Plot of Bootstrap LMM Beta Estimates
#' @description
#' Renders a horizontal forest plot of the bootstrap-derived interaction beta
#' estimates with 95% bootstrap CIs. Each row is one target marker; the
#' observed (point) estimate from the primary LMM is plotted as a coloured
#' point, the bootstrap CI as a horizontal segment, and the bootstrap median
#' as a small dashed vertical tick. A reference line at beta=0 marks the null.
#'
#' @param boot_summary Dataframe from run_lmm_bootstrap_ci()$summary_df. Must
#'   contain Marker, Median_Beta_Boot, CI_Lower_2.5, CI_Upper_97.5,
#'   Pct_FDR_Significant.
#' @param observed_df Optional dataframe with columns Marker and
#'   Estimate_Interaction (primary LMM point estimate). If provided, plotted
#'   alongside the bootstrap median.
#' @param title Plot title.
#' @return A ggplot object, or NULL if no rows.
plot_lmm_forest <- function(boot_summary, observed_df = NULL,
                            title = "Bootstrap 95% CI of LMM Interaction Betas") {
  if (is.null(boot_summary) || nrow(boot_summary) == 0) return(NULL)
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)

  plot_df <- boot_summary
  if (!is.null(observed_df) && all(c("Marker", "Estimate_Interaction") %in% colnames(observed_df))) {
    plot_df <- merge(plot_df,
                     observed_df[, c("Marker", "Estimate_Interaction")],
                     by = "Marker", all.x = TRUE)
  } else {
    plot_df$Estimate_Interaction <- plot_df$Median_Beta_Boot
  }

  # Order by observed effect size (ascending — most negative at top)
  plot_df$Marker <- factor(plot_df$Marker,
                           levels = plot_df$Marker[order(plot_df$Estimate_Interaction)])

  plot_df$Label <- sprintf("%s  (%.0f%% FDR<0.05)", plot_df$Marker, plot_df$Pct_FDR_Significant)

  p <- ggplot(plot_df, aes(y = Marker)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = CI_Lower_2.5, xmax = CI_Upper_97.5),
                   height = 0.25, color = "gray30") +
    geom_point(aes(x = Median_Beta_Boot), shape = 4, size = 2.5, color = "gray40") +
    geom_point(aes(x = Estimate_Interaction), shape = 19, size = 3.5, color = "#B2182B") +
    geom_text(aes(x = CI_Upper_97.5,
                  label = sprintf("  %.0f%% FDR", Pct_FDR_Significant)),
              hjust = 0, size = 3.2, color = "gray30") +
    labs(
      title    = title,
      subtitle = "Red dot: observed beta | X: bootstrap median | Bar: 95% bootstrap CI",
      x        = "Time x Group interaction (beta, hybrid scale)",
      y        = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
      panel.grid.minor = element_blank(),
      axis.text.y   = element_text(size = 11)
    ) +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(8, 60, 8, 8))

  return(p)
}