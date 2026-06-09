# src/04_longitudinal_lmm.R
# ==============================================================================
# STEP 04: LONGITUDINAL LMM ANALYSIS
# Description: Fits Linear Mixed Models to evaluate Time x Clinical Interaction.
#              Includes LOO Sensitivity Analysis to protect against outlier bias.
# ==============================================================================

source("R/utils_io.R")
source("R/modules_longitudinal.R")
source("R/modules_viz.R")

message("\n=== PIPELINE STEP 4: LONGITUDINAL ANALYSIS (LMM) ===")

# 1. Configuration & Data Parsing
# ------------------------------------------------------------------------------
if (!exists("config")) stop("[FATAL] Global configuration object not detected in environment.")

input_rds <- step_output_path(config, 1, sprintf("data_processed_%s_longitudinal", config$project_name), "rds")
if (!file.exists(input_rds)) stop(sprintf("[FATAL] Step 01 longitudinal processed dataset not found at: %s", input_rds))

out_dir <- step_dir(config, 4, create = TRUE)

message(sprintf("[Data] Loading Joint Pharmacodynamic Dataset: %s", basename(input_rds)))
DATA <- readRDS(input_rds)

df_long <- DATA$hybrid_data_raw
meta_long <- DATA$metadata

df_model <- cbind(meta_long, as.data.frame(df_long[, DATA$hybrid_markers]))

df_model$Timepoint <- factor(df_model$Timepoint, levels = c("T0", "T1"))
df_model$Patient_ID <- as.factor(df_model$Patient_ID)

ref_level <- config$clinical$non_responder_label
avail_levels <- unique(as.character(df_model$Group))

if (ref_level %in% avail_levels) {
  df_model$Group <- factor(df_model$Group, levels = c(ref_level, setdiff(avail_levels, ref_level)))
} else {
  df_model$Group <- factor(df_model$Group)
}

n_t0 <- sum(df_model$Timepoint == "T0")
n_t1 <- sum(df_model$Timepoint == "T1")

tp_per_pid <- split(as.character(df_model$Timepoint), as.character(df_model$Patient_ID))
n_paired    <- sum(vapply(tp_per_pid, function(v) ("T0" %in% v) && ("T1" %in% v), logical(1)))
n_only_t0   <- sum(vapply(tp_per_pid, function(v) ("T0" %in% v) && !("T1" %in% v), logical(1)))
n_only_t1   <- sum(vapply(tp_per_pid, function(v) !("T0" %in% v) && ("T1" %in% v), logical(1)))

message(sprintf("   [Data] Longitudinal dimensions established: %d matrix observations (T0: %d, T1: %d)",
                nrow(df_model), n_t0, n_t1))
message(sprintf("   [Data] Patient pairing: %d paired (T0+T1) | %d T0-only | %d T1-only",
                n_paired, n_only_t0, n_only_t1))
if (n_t0 != n_t1) {
  message(sprintf(
    "   [Data] Note: %d patient(s) have only one timepoint. LMM is valid under MAR (Missing At Random); verify that dropout is not response-correlated. A paired-only delta sensitivity is reported below for FDR-significant markers.",
    abs(n_t0 - n_t1)
  ))
}

covariates_list <- if (!is.null(config$clinical$covariates)) unlist(config$clinical$covariates) else NULL

# --- SAFETY CHECK: Prevent Silent Cohort Decimation via Covariates ---
safe_max_na_covariates <- if(!is.null(config$qc$max_na_covariates)) config$qc$max_na_covariates else 0.10

if (!is.null(covariates_list) && length(intersect(covariates_list, colnames(df_model))) > 0) {
  valid_covs <- intersect(covariates_list, colnames(df_model))
  na_rows <- sum(!complete.cases(df_model[, valid_covs, drop = FALSE]))
  pct_na <- na_rows / nrow(df_model)
  
  if (pct_na > safe_max_na_covariates) {
    stop(sprintf("[FATAL] Security Gate Triggered: Covariates induce %.1f%% missingness, exceeding the safety threshold (%.1f%%). Aborting to prevent silent cohort decimation.", pct_na * 100, safe_max_na_covariates * 100))
  } else if (pct_na > 0) {
    message(sprintf("   [Data] Notice: Covariate inclusion will drop %d observations (%.1f%% of cohort).", na_rows, pct_na * 100))
  }
}

colors_viz <- get_clinical_colors(config)

# 2. Linear Mixed Models Computation
# ------------------------------------------------------------------------------
message("\n[Stats] Initiating Linear Mixed Models (LMM) for longitudinal variance extraction...")

results_list <- list()
pb <- txtProgressBar(min = 0, max = length(DATA$hybrid_markers), style = 3)

for (i in seq_along(DATA$hybrid_markers)) {
  mk <- DATA$hybrid_markers[i]
  res <- fit_feature_lmm(data_long = df_model, feature = mk, group_col = "Group", 
                         time_col = "Timepoint", id_col = "Patient_ID", 
                         covariates = covariates_list)
  results_list[[mk]] <- res
  setTxtProgressBar(pb, i)
}
close(pb)

df_results <- do.call(rbind, results_list)
rownames(df_results) <- NULL

# 3. Multi-Testing Correction (FDR)
# ------------------------------------------------------------------------------
df_results <- df_results %>%
  dplyr::filter(!is.na(P_Value_Interaction)) %>%
  dplyr::mutate(FDR_Interaction = p.adjust(P_Value_Interaction, method = "BH")) %>%
  dplyr::arrange(P_Value_Interaction)

n_sig_raw <- sum(df_results$P_Value_Interaction < 0.05, na.rm = TRUE)
n_sig_fdr <- sum(df_results$FDR_Interaction < 0.05, na.rm = TRUE)

message(sprintf("\n[Stats] LMM integration successful. Identified %d raw significant markers, %d survive FDR adjustment.", 
                n_sig_raw, n_sig_fdr))

# 4. LOO Sensitivity Security Check (FDR subset only)
# ------------------------------------------------------------------------------
df_results$Max_P_Value_LOO <- NA

if (n_sig_fdr > 0) {
  message("   [Stats] Executing Leave-One-Out (LOO) Sensitivity protocol on topological drivers...")
  sig_markers <- df_results$Marker[which(df_results$FDR_Interaction < 0.05)]

  for (mk in sig_markers) {
    max_p <- run_loo_sensitivity(data_long = df_model, feature = mk, group_col = "Group",
                                 time_col = "Timepoint", id_col = "Patient_ID",
                                 covariates = covariates_list)

    df_results$Max_P_Value_LOO[df_results$Marker == mk] <- max_p

    if (!is.na(max_p) && max_p < 0.05) {
      message(sprintf("      -> %s: Structurally Robust (Max P-Value = %.4f)", mk, max_p))
    } else {
      warning(sprintf("      -> %s: OUTLIER BIAS DETECTED (Max P-Value spikes to %.4f upon LOO)", mk, max_p))
    }
  }

  # Flag LMM singular fits among the FDR-significant set so reviewers can see
  # them at a glance. Singular RE => random-intercept variance estimated as 0
  # (LMM degenerates to OLS-with-clustering). Fixed-effect inference remains
  # valid but the paired-only sensitivity below is the relevant cross-check.
  sing_mask <- df_results$Marker %in% sig_markers & isTRUE(df_results$Is_Singular)
  sing_markers <- df_results$Marker[df_results$Marker %in% sig_markers &
                                    !is.na(df_results$Is_Singular) &
                                    df_results$Is_Singular]
  if (length(sing_markers) > 0) {
    message(sprintf(
      "   [Stats] LMM singular random-effect fit on %d/%d FDR-significant marker(s): %s",
      length(sing_markers), length(sig_markers), paste(sing_markers, collapse = ", ")
    ))
  }
}

# 4b. Paired-Only Delta Sensitivity (FDR-significant subset only)
# OLS on within-patient delta on the subset of patients with both T0 and T1.
# Sidesteps random-intercept identification; under balanced pairing the slope
# equals the LMM Time x Group interaction term.
# ------------------------------------------------------------------------------
paired_results <- NULL
if (n_sig_fdr > 0) {
  sig_markers <- df_results$Marker[which(df_results$FDR_Interaction < 0.05)]
  message(sprintf(
    "\n[Stats] Paired-only delta sensitivity (n_paired=%d) on %d FDR-significant marker(s)...",
    n_paired, length(sig_markers)
  ))
  paired_results <- tryCatch(
    run_paired_only_sensitivity(
      data_long = df_model, features = sig_markers,
      group_col = "Group", time_col = "Timepoint", id_col = "Patient_ID"
    ),
    error = function(e) {
      warning(sprintf("[Stats] Paired-only sensitivity failed: %s", e$message))
      NULL
    }
  )
  if (!is.null(paired_results) && nrow(paired_results) > 0) {
    for (i in seq_len(nrow(paired_results))) {
      r <- paired_results[i, ]
      message(sprintf("      -> %s: delta=%.3f (SE=%.3f, p=%.4f, FDR=%.4f, n_pairs=%d)",
                      r$Marker, r$Estimate_Delta, r$Std_Error, r$P_Value, r$FDR, r$N_Pairs))
    }
  }
}

# 5. Covariate Sensitivity Analysis (optional, config-driven)
# ------------------------------------------------------------------------------
sensitivity_results <- NULL
sensitivity_covariates <- config$clinical$sensitivity_covariates

if (!is.null(sensitivity_covariates) && length(sensitivity_covariates) > 0) {
  sig_markers_sens <- df_results$Marker[!is.na(df_results$FDR_Interaction) &
                                          df_results$FDR_Interaction < 0.05]
  if (length(sig_markers_sens) > 0) {
    message(sprintf("\n[Stats] Running covariate sensitivity analysis on %d FDR-significant marker(s)...",
                    length(sig_markers_sens)))
    message(sprintf("   [Stats] Adjustment covariates: %s", paste(sensitivity_covariates, collapse = ", ")))

    if (!file.exists(config$input_file_t0)) {
      warning(sprintf("[Stats] Covariate sensitivity skipped: raw T0 file not found at '%s'.", config$input_file_t0))
    } else tryCatch({
      df_raw_t0 <- readxl::read_excel(config$input_file_t0)
      avail_covs <- intersect(sensitivity_covariates, colnames(df_raw_t0))
      if (length(avail_covs) == 0) stop("None of the sensitivity covariates found in raw T0 file.")

      df_clin_sens <- df_raw_t0[, c("Patient_ID", avail_covs), drop = FALSE]
      df_model_adj <- dplyr::left_join(df_model, df_clin_sens, by = "Patient_ID")

      results_sens <- purrr::map_dfr(sig_markers_sens, function(mk) {
        fit_feature_lmm(data_long = df_model_adj, feature = mk,
                        group_col = "Group", time_col = "Timepoint", id_col = "Patient_ID",
                        covariates = avail_covs)
      })
      results_sens$FDR_Adj <- p.adjust(results_sens$P_Value_Interaction, method = "BH")
      results_sens$Covariates_Used <- paste(avail_covs, collapse = "; ")
      sensitivity_results <- results_sens

      for (i in seq_len(nrow(results_sens))) {
        mk_s <- results_sens$Marker[i]
        p_s  <- results_sens$P_Value_Interaction[i]
        fdr_s <- results_sens$FDR_Adj[i]
        n_s  <- results_sens$N_Observations[i]
        message(sprintf("      -> %s: adj. p=%.4f  FDR=%.4f  n=%d%s",
                        mk_s, p_s, fdr_s, n_s,
                        if (!is.na(fdr_s) && fdr_s < 0.05) "  [sig]" else ""))
      }
    }, error = function(e) {
      warning(sprintf("[Stats] Covariate sensitivity analysis failed: %s", e$message))
    })
  } else {
    message("[Stats] No FDR-significant markers to adjust — skipping sensitivity analysis.")
  }
}

# 5b. Split-Group Supplementary LMM (optional, config-driven)
# When the non-responder label collapses several raw clinical codes (e.g.
# SD_PD = {SD=3, PD=4}), this block re-fits the LMM with the non-responder
# side split into its sub-levels, so reviewers can see whether the primary
# interaction is graded across response severity or driven by one sub-group.
# Gated by clinical.split_nonresponder: true in config.
# ------------------------------------------------------------------------------
split_results <- NULL
split_cfg <- isTRUE(config$clinical$split_nonresponder)
if (split_cfg && n_sig_fdr > 0) {
  raw_codes <- config$clinical$mapping[[config$clinical$non_responder_label]]
  resp_codes <- config$clinical$mapping[[config$clinical$responder_label]]
  tgt_col   <- config$clinical$target_column
  raw_path  <- if (!is.null(config$input_file_t0)) config$input_file_t0 else config$input_file
  sig_markers_split <- df_results$Marker[!is.na(df_results$FDR_Interaction) &
                                          df_results$FDR_Interaction < 0.05]

  if (length(raw_codes) >= 2 && file.exists(raw_path)) {
    tryCatch({
      df_raw_t0 <- readxl::read_excel(raw_path)
      if (!(tgt_col %in% colnames(df_raw_t0))) {
        warning(sprintf("[Stats] Split-group LMM: column '%s' missing in %s — skipping.",
                        tgt_col, basename(raw_path)))
      } else {
        pid_vec  <- as.character(df_raw_t0$Patient_ID)
        code_vec <- as.integer(df_raw_t0[[tgt_col]])
        # Map every raw code to a label string (sub-level)
        lab_vec <- as.character(code_vec)
        # Patients outside responder/non-responder universe are excluded
        keep <- code_vec %in% c(raw_codes, resp_codes)
        pid_vec <- pid_vec[keep]
        lab_vec <- lab_vec[keep]
        patient_subgroup <- setNames(lab_vec, pid_vec)
        # Reference: clinically "most adverse" sub-code. For outcome scales
        # where higher integer = worse (e.g. RECIST 1=RC, 2=RP, 3=SD, 4=PD),
        # max == PD. The user can override this via config$clinical$split_reference_code.
        ref_lvl <- if (!is.null(config$clinical$split_reference_code)) {
          as.character(config$clinical$split_reference_code)
        } else {
          as.character(max(as.integer(raw_codes)))
        }

        message(sprintf(
          "\n[Stats] Split-group supplementary LMM (reference code=%s) on %d FDR-significant marker(s)...",
          ref_lvl, length(sig_markers_split)
        ))
        split_results <- run_splitgroup_lmm(
          data_long       = df_model,
          features        = sig_markers_split,
          patient_subgroup = patient_subgroup,
          ref_level       = ref_lvl
        )
        if (!is.null(split_results) && nrow(split_results) > 0) {
          for (i in seq_len(nrow(split_results))) {
            r <- split_results[i, ]
            message(sprintf("      -> %s [Group %s vs %s]: beta=%+.3f (SE=%.3f, p=%.4f, n=%d%s)",
                            r$Marker, r$Level, r$Reference, r$Estimate_Interaction,
                            r$Std_Error, r$P_Value, r$N_Observations,
                            if (r$Is_Singular) ", singular" else ""))
          }
        }
      }
    }, error = function(e) {
      warning(sprintf("[Stats] Split-group LMM failed: %s", e$message))
    })
  } else if (length(raw_codes) < 2) {
    message("[Stats] split_nonresponder: true but non-responder maps to a single code — nothing to split.")
  }
}

# 6. Bootstrap CI on FDR-significant Interaction Betas (optional, config-driven)
# Patient-level cluster bootstrap that re-runs the LMM panel + BH-FDR on each
# resample. Reports bootstrap 95% CI and the fraction of resamples in which
# each target marker still survives FDR < fdr_threshold.
# ------------------------------------------------------------------------------
bootstrap_results <- NULL
boot_cfg <- config$lmm_bootstrap

if (!is.null(boot_cfg) && isTRUE(boot_cfg$enabled) && n_sig_fdr > 0) {
  n_boot_iter <- if (!is.null(boot_cfg$n_boot)) as.integer(boot_cfg$n_boot) else 500L
  boot_seed   <- if (!is.null(boot_cfg$seed))   as.integer(boot_cfg$seed)   else 2026L

  target_markers <- df_results$Marker[!is.na(df_results$FDR_Interaction) &
                                        df_results$FDR_Interaction < 0.05]
  message(sprintf("\n[Stats] Running patient-level cluster bootstrap (n_boot=%d) on %d FDR-significant marker(s)...",
                  n_boot_iter, length(target_markers)))

  bootstrap_results <- tryCatch(
    run_lmm_bootstrap_ci(
      data_long      = df_model,
      all_markers    = DATA$hybrid_markers,
      target_markers = target_markers,
      group_col      = "Group",
      time_col       = "Timepoint",
      id_col         = "Patient_ID",
      covariates     = covariates_list,
      n_boot         = n_boot_iter,
      fdr_threshold  = 0.05,
      seed           = boot_seed,
      progress_message = TRUE
    ),
    error = function(e) {
      warning(sprintf("[Stats] LMM bootstrap failed: %s", e$message))
      NULL
    }
  )

  if (!is.null(bootstrap_results) && nrow(bootstrap_results$summary_df) > 0) {
    for (i in seq_len(nrow(bootstrap_results$summary_df))) {
      r <- bootstrap_results$summary_df[i, ]
      message(sprintf("      -> %s: median=%.3f [%.3f, %.3f]  %%FDR<0.05=%.1f%%  (n_valid=%d)",
                      r$Marker, r$Median_Beta_Boot, r$CI_Lower_2.5, r$CI_Upper_97.5,
                      r$Pct_FDR_Significant, r$N_Valid_Iterations))
    }
  }
}

# 7. Output Serialization
# ------------------------------------------------------------------------------
json_path <- file.path(out_dir, sprintf("Machine_Metrics_LMM_%s.json", config$project_name))
sig_singular_markers <- if (n_sig_fdr > 0) {
  df_results$Marker[!is.na(df_results$FDR_Interaction) &
                    df_results$FDR_Interaction < 0.05 &
                    !is.na(df_results$Is_Singular) &
                    df_results$Is_Singular]
} else character(0)

machine_output <- list(
  project_name = config$project_name,
  clinical_target = config$clinical$target_column,
  model_type = "LMM_Interaction",
  n_observations = nrow(df_model),
  n_patients = length(unique(df_model$Patient_ID)),
  n_paired = n_paired,
  n_only_T0 = n_only_t0,
  n_only_T1 = n_only_t1,
  significant_features_fdr = n_sig_fdr,
  primary_markers_singular = sig_singular_markers,
  full_results = df_results,
  split_group_lmm = if (!is.null(split_results) && nrow(split_results) > 0) {
    list(
      method  = "LMM Time x Group with non-responder split into raw clinical sub-levels",
      results = split_results,
      note    = "Use to inspect whether the primary interaction is graded vs reference (most adverse code) or driven by one sub-group; n smaller than the dichotomized model."
    )
  } else NULL,
  paired_sensitivity = if (!is.null(paired_results) && nrow(paired_results) > 0) {
    list(
      method         = "OLS on within-patient delta (T1 - T0), paired-only subset",
      n_paired_total = n_paired,
      n_only_T0      = n_only_t0,
      n_only_T1      = n_only_t1,
      results        = paired_results,
      note           = "Under balanced pairing OLS slope == LMM Time:Group interaction term. Sign + p-value agreement with the LMM table is the relevant robustness signal."
    )
  } else NULL,
  covariate_sensitivity = if (!is.null(sensitivity_results)) sensitivity_results else NULL,
  bootstrap_ci = if (!is.null(bootstrap_results)) {
    list(
      method        = "Patient-level cluster bootstrap with per-iteration BH-FDR re-computation",
      n_boot        = bootstrap_results$n_boot,
      seed          = bootstrap_results$seed,
      fdr_threshold = bootstrap_results$fdr_threshold,
      summary       = bootstrap_results$summary_df,
      note          = sprintf(
        "Cluster bootstrap (B=%d) resamples patient IDs with replacement, preserving T0/T1 pairing within patient. Beta CI excludes 0 + %%FDR<0.05 reported per marker.",
        bootstrap_results$n_boot
      )
    )
  } else NULL
)

if (requireNamespace("jsonlite", quietly = TRUE)) {
  jsonlite::write_json(machine_output, json_path, pretty = TRUE, auto_unbox = TRUE)
}

excel_path <- file.path(out_dir, sprintf("Longitudinal_LMM_Report_%s.xlsx", config$project_name))
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "LMM_Interaction_Results")
openxlsx::writeData(wb, "LMM_Interaction_Results", df_results)

sig_style <- openxlsx::createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
p_col_idx <- which(names(df_results) == "P_Value_Interaction")

if (length(p_col_idx) > 0) {
  openxlsx::conditionalFormatting(wb, "LMM_Interaction_Results", cols = p_col_idx,
                                  rows = 2:(nrow(df_results)+1), rule = "< 0.05", style = sig_style)
}

if (!is.null(sensitivity_results) && nrow(sensitivity_results) > 0) {
  openxlsx::addWorksheet(wb, "LMM_Covariate_Sensitivity")
  openxlsx::writeData(wb, "LMM_Covariate_Sensitivity", sensitivity_results)
}

if (!is.null(paired_results) && nrow(paired_results) > 0) {
  openxlsx::addWorksheet(wb, "LMM_Paired_Delta_Sensitivity")
  openxlsx::writeData(wb, "LMM_Paired_Delta_Sensitivity", paired_results)
}

if (!is.null(split_results) && nrow(split_results) > 0) {
  openxlsx::addWorksheet(wb, "LMM_Split_Group_Supplementary")
  openxlsx::writeData(wb, "LMM_Split_Group_Supplementary", split_results)
}

if (!is.null(bootstrap_results) && nrow(bootstrap_results$summary_df) > 0) {
  openxlsx::addWorksheet(wb, "LMM_Bootstrap_CI")
  openxlsx::writeData(wb, "LMM_Bootstrap_CI", bootstrap_results$summary_df)
}

openxlsx::saveWorkbook(wb, excel_path, overwrite = TRUE)
message(sprintf("   [Output] Mathematical models matrix saved: %s", basename(excel_path)))

# 7. Volcano Plot & Trajectory Rendering
# ------------------------------------------------------------------------------
plot_path <- file.path(out_dir, sprintf("Volcano_LMM_%s.pdf", config$project_name))
pdf(plot_path, width = 9, height = 7)
tryCatch({
  p_volcano <- plot_lmm_volcano(df_results, title = sprintf("LMM Pharmacodynamics: Time x %s", config$clinical$target_column))
  if (!is.null(p_volcano)) print(p_volcano)
}, error = function(e) warning(paste("Volcano plot rendering failed:", e$message)))
dev.off()

message("\n[Viz] Generating Patient Trajectory vectors...")

top_df <- df_results %>% dplyr::filter(FDR_Interaction < 0.05)
sig_type <- "FDR"

if (nrow(top_df) == 0) {
  top_df <- df_results %>% dplyr::arrange(P_Value_Interaction) %>% head(4)
  sig_type <- "RAW"
  message("   [Viz] Null set for FDR boundaries. Defaulting plot limits to top 4 absolute P-values.")
}

if (nrow(top_df) > 0) {
  traj_path <- file.path(out_dir, sprintf("Trajectories_LMM_%s.pdf", config$project_name))
  pdf(traj_path, width = 8, height = 6)
  
  for (i in 1:nrow(top_df)) {
    mk <- top_df$Marker[i]
    pval_disp <- if (sig_type == "FDR") top_df$FDR_Interaction[i] else top_df$P_Value_Interaction[i]
    
    tryCatch({
      p_traj <- plot_lmm_trajectories(
        data_long = df_model,
        feature = mk,
        group_col = "Group",
        time_col = "Timepoint",
        id_col = "Patient_ID",
        colors = colors_viz,
        p_val = pval_disp
      )
      print(p_traj)
    }, error = function(e) warning(sprintf("Trajectory matrix error on marker %s: %s", mk, e$message)))
  }
  
  dev.off()
  message(sprintf("   [Output] Temporal trajectory visualizations exported: %s", basename(traj_path)))
}

# Forest plot: bootstrap 95% CI of interaction betas (only if bootstrap ran)
if (!is.null(bootstrap_results) && nrow(bootstrap_results$summary_df) > 0) {
  forest_path <- file.path(out_dir, sprintf("Forest_Bootstrap_LMM_%s.pdf", config$project_name))
  pdf(forest_path, width = 9, height = max(3, 1.0 + 0.5 * nrow(bootstrap_results$summary_df)))
  tryCatch({
    p_forest <- plot_lmm_forest(
      boot_summary = bootstrap_results$summary_df,
      observed_df  = df_results[, c("Marker", "Estimate_Interaction")],
      title        = sprintf("Bootstrap 95%% CI of LMM Betas (B=%d) - %s",
                             bootstrap_results$n_boot, config$project_name)
    )
    if (!is.null(p_forest)) print(p_forest)
  }, error = function(e) warning(sprintf("Forest plot failed: %s", e$message)))
  dev.off()
  message(sprintf("   [Output] Bootstrap forest plot exported: %s", basename(forest_path)))
}

message("=== STEP 4 COMPLETE ===\n")