# R/utils_io.R
# ==============================================================================
# INFRASTRUCTURE UTILITIES
# Description: Config loading, Logging setup, and basic I/O helpers.
# ==============================================================================

library(yaml)
library(readxl)
library(dplyr)
library(openxlsx)

#' @title Load Global Configuration
#' @description Loads the YAML config and validates critical paths.
#' @param config_path Path to the YAML file.
#' @return A list containing configuration parameters.
load_config <- function(config_path = "config/global_params.yml") {
  if (!file.exists(config_path)) {
    stop(sprintf("Configuration file not found at: %s", config_path))
  }
  
  config <- yaml::read_yaml(config_path)
  
  # If project_name is defined, nest it inside output_root
  if (!is.null(config$project_name) && config$project_name != "") {
    config$output_root <- file.path(config$output_root, config$project_name)
  }
  
  message(sprintf("[System] Configuration loaded from %s", config_path))

  return(config)
}

# ==============================================================================
# RUN MANAGEMENT
# A single main.R invocation = one immutable, timestamped run directory under
# `output_root` (e.g. results/20260605_143000/). Every step derives its paths
# from config$output_root, so isolating the run at the root keeps all cross-step
# artifacts of one invocation self-contained and reproducible. A `latest`
# pointer (symlink + portable text file) lets follow-up tooling resolve the most
# recent run without hardcoding a timestamp.
# ==============================================================================

#' @title Build a Run Identifier
#' @description Timestamp-based run id, optionally suffixed with a human label
#'   (config `run_label` or the RUN_LABEL environment variable, env wins).
#' @param label Optional character label appended after the timestamp.
#' @return Character run id, e.g. "20260605_143000_disc73".
make_run_id <- function(label = NULL) {
  env_label <- Sys.getenv("RUN_LABEL", unset = "")
  if (nzchar(env_label)) label <- env_label

  id <- format(Sys.time(), "%Y%m%d_%H%M%S")
  if (!is.null(label) && nzchar(label)) {
    safe <- gsub("[^A-Za-z0-9._-]+", "-", label)
    id <- paste0(id, "_", safe)
  }
  id
}

#' @title Publish the `latest` Run Pointer
#' @description Records `run_id` as the most recent run via a portable text file
#'   (`latest_run.txt`, the source of truth) and a best-effort relative symlink
#'   (`latest`) for convenient path navigation. The symlink is non-fatal: some
#'   filesystems do not support it, in which case the text pointer still works.
#' @param results_base Base output directory (config$output_root, e.g. "results").
#' @param run_id The run id to publish.
#' @return Invisibly, the pointer file path.
update_latest_pointer <- function(results_base, run_id) {
  if (!dir.exists(results_base)) dir.create(results_base, recursive = TRUE)

  pointer <- file.path(results_base, "latest_run.txt")
  writeLines(run_id, pointer)

  symlink <- file.path(results_base, "latest")
  if (file.exists(symlink) || nzchar(Sys.readlink(symlink))) unlink(symlink, force = TRUE)
  ok <- tryCatch(file.symlink(run_id, symlink),
                 warning = function(w) FALSE, error = function(e) FALSE)
  if (!isTRUE(ok))
    message("[RUN] 'latest' symlink not created (filesystem may not support it); using latest_run.txt.")

  invisible(pointer)
}

#' @title Resolve the Run Root for Follow-up Tooling
#' @description Returns the directory of the run that single-step re-runners
#'   (diagnostics/) should read from. Resolution order: RUN_ID env override,
#'   then the `latest_run.txt` pointer, then the `latest` symlink.
#' @param results_base Base output directory (default "results").
#' @return Absolute-or-relative path to the resolved run root.
resolve_run_root <- function(results_base = "results") {
  env_id <- Sys.getenv("RUN_ID", unset = "")
  if (nzchar(env_id)) {
    run_root <- file.path(results_base, env_id)
    if (!dir.exists(run_root))
      stop(sprintf("[RUN] RUN_ID '%s' set but %s does not exist.", env_id, run_root))
    return(run_root)
  }

  pointer <- file.path(results_base, "latest_run.txt")
  if (file.exists(pointer)) {
    run_id   <- trimws(readLines(pointer, warn = FALSE)[1])
    run_root <- file.path(results_base, run_id)
    if (dir.exists(run_root)) return(run_root)
  }

  symlink <- file.path(results_base, "latest")
  target  <- Sys.readlink(symlink)
  if (nzchar(target)) {
    run_root <- if (grepl("^(/|[A-Za-z]:)", target)) target else file.path(results_base, target)
    if (dir.exists(run_root)) return(run_root)
  }

  stop(sprintf("[RUN] No resolvable run under '%s'. Run main.R first, or set RUN_ID.", results_base))
}

#' @title Write the Run Provenance Manifest
#' @description Serializes a self-describing record (git state, R version, seed,
#'   experiments + passes executed) to the run root, making each results folder
#'   traceable to the exact code and parameters that produced it.
#' @param path Output path for the YAML manifest.
#' @param run_id The run id.
#' @param config The base configuration object.
#' @param experiments_meta Named list mapping experiment name -> character vector
#'   of passes executed.
#' @return Invisibly, the manifest path.
write_run_manifest <- function(path, run_id, config, experiments_meta) {
  git_out <- function(args) tryCatch(
    system2("git", args, stdout = TRUE, stderr = FALSE),
    error = function(e) character(0), warning = function(w) character(0))

  sha    <- git_out(c("rev-parse", "HEAD"))
  branch <- git_out(c("rev-parse", "--abbrev-ref", "HEAD"))
  dirty  <- length(git_out(c("status", "--porcelain"))) > 0

  manifest <- list(
    run_id      = run_id,
    timestamp   = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    git_sha     = if (length(sha))    sha[1]    else NA_character_,
    git_branch  = if (length(branch)) branch[1] else NA_character_,
    git_dirty   = dirty,
    r_version   = R.version.string,
    seed        = if (!is.null(config$stats$seed)) as.integer(config$stats$seed) else NULL,
    run_label   = if (!is.null(config$run_label)) config$run_label else NULL,
    experiments = experiments_meta
  )

  yaml::write_yaml(manifest, path)
  message(sprintf("[RUN] Manifest written: %s", path))
  invisible(path)
}

#' @title Per-Step Output Directory Names
#' @description Single source of truth for the canonical step subdirectory names
#'   (indexed by step number 1..6). Every step derives its paths from these, so a
#'   rename happens in exactly one place.
.STEP_DIRS <- c(
  "01_data_processing",
  "02_visualization",
  "03_statistical_analysis",
  "04_longitudinal_analysis",
  "05_network_analysis",
  "06_machine_learning"
)

#' @title Resolve a Step's Output Directory
#' @param config Pass config (uses \code{config$output_root}).
#' @param step Integer step number (1..6).
#' @param create If TRUE, create the directory (recursively) when missing.
#' @return The step directory path under the run's experiment root.
step_dir <- function(config, step, create = FALSE) {
  d <- file.path(config$output_root, .STEP_DIRS[[step]])
  if (isTRUE(create) && !dir.exists(d)) dir.create(d, recursive = TRUE)
  d
}

#' @title Build a Path to a Step Artifact
#' @description Centralizes the \code{file.path(config$output_root, "0N_...", ...)}
#'   pattern repeated across \code{src/}. The step subdirectory name lives only in
#'   \code{.STEP_DIRS}.
#' @param config Pass config (uses \code{config$output_root}).
#' @param step Integer step number (1..6) owning the artifact.
#' @param name File name without extension (or full name if \code{ext} is NULL).
#' @param ext Optional extension (no dot), e.g. "rds", "json".
#' @return The full artifact path.
step_output_path <- function(config, step, name, ext = NULL) {
  fname <- if (!is.null(ext) && nzchar(ext)) paste0(name, ".", ext) else name
  file.path(step_dir(config, step), fname)
}


# ── Cohort ledger — the auditable patient-flow record ─────────────────────────
# WHY THIS EXISTS. Nothing used to record *why* a patient left the cohort, so the
# CONSORT figure reconstructed the flow by arithmetic — and every reconstruction was
# wrong (the rendered figure said 69 patients, the docs said 91 raw; the input file has
# 94). A ledger row is appended at each filtering point instead, and `ledger_validate()`
# makes a discontinuity FATAL rather than something a reader discovers in review.
# It is also how the silent outcome-eligibility leak became visible: `unmapped_values`
# records the target-column values that matched no mapping entry.

#' @title Start an empty cohort ledger
#' @return A zero-row data.frame with the ledger schema.
ledger_new <- function() {
  data.frame(stage = character(), n_in = integer(), n_out = integer(),
             n_dropped = integer(), reason = character(),
             n_resp = integer(), n_nonresp = integer(),
             dropped_ids = character(), unmapped_values = character(),
             stringsAsFactors = FALSE)
}

#' @title Append one filtering point to the ledger
#' @param led Ledger data.frame (\code{ledger_new()} or a previous result).
#' @param stage Short machine-readable stage name, e.g. "outcome_eligibility".
#' @param n_in,n_out Cohort size entering / leaving this stage.
#' @param reason Human-readable exclusion reason (shown on the CONSORT panel).
#' @param by_group Optional named counts of the surviving cohort per outcome class.
#' @param dropped_ids Optional character vector of the patient IDs removed here.
#' @param unmapped_values Optional raw values that matched no outcome mapping.
ledger_add <- function(led, stage, n_in, n_out, reason = NA_character_,
                       by_group = NULL, dropped_ids = NULL, unmapped_values = NULL) {
  cl <- function(x) if (is.null(x) || !length(x)) NA_character_ else paste(x, collapse = ", ")
  rbind(led, data.frame(
    stage = as.character(stage), n_in = as.integer(n_in), n_out = as.integer(n_out),
    n_dropped = as.integer(n_in) - as.integer(n_out), reason = as.character(reason),
    n_resp    = if (!is.null(by_group) && length(by_group) >= 2) as.integer(by_group[[2]]) else NA_integer_,
    n_nonresp = if (!is.null(by_group) && length(by_group) >= 1) as.integer(by_group[[1]]) else NA_integer_,
    dropped_ids = cl(dropped_ids), unmapped_values = cl(unique(unmapped_values)),
    stringsAsFactors = FALSE))
}

#' @title Validate that a cohort ledger reconciles
#' @description Two invariants: each row's arithmetic closes, and consecutive stages
#'   chain (\code{n_out[i] == n_in[i+1]}). A break means the flow diagram would be
#'   unreconcilable with the data, so this is FATAL by design — the same posture the
#'   Pass-0 cohort split takes.
#' @param led Ledger data.frame.
#' @param context Optional label included in the error message.
#' @return The ledger, invisibly, when valid.
ledger_validate <- function(led, context = "") {
  if (is.null(led) || !nrow(led)) return(invisible(led))
  bad <- which(led$n_in - led$n_dropped != led$n_out)
  if (length(bad))
    stop(sprintf("[FATAL] Cohort ledger %sdoes not reconcile at stage(s) %s: n_in - n_dropped != n_out.",
                 if (nzchar(context)) paste0("(", context, ") ") else "",
                 paste(led$stage[bad], collapse = ", ")))
  if (nrow(led) > 1) {
    brk <- which(led$n_out[-nrow(led)] != led$n_in[-1])
    if (length(brk))
      stop(sprintf("[FATAL] Cohort ledger %sis discontinuous: stage '%s' ends at n=%d but '%s' starts at n=%d.",
                   if (nzchar(context)) paste0("(", context, ") ") else "",
                   led$stage[brk[1]], led$n_out[brk[1]],
                   led$stage[brk[1] + 1], led$n_in[brk[1] + 1]))
  }
  invisible(led)
}

#' @title Save Quality Control Report to Excel
#' @description Generates a multi-sheet Excel report with filtering statistics, group breakdowns, and categorized marker lists.
#' @param qc_list A list containing summary stats, dataframes of dropped items, group mappings, and optionally imputed_details.
#' @param out_path Path to save the .xlsx file.
#' @param config Configuration object (optional) containing 'qc_reporting' settings for marker categorization.
#' @param ledger Optional cohort ledger (\code{ledger_new()}/\code{ledger_add()}). When supplied it
#'   is written as its own sheet and overrides the Summary's reconstructed grand totals.
save_qc_report <- function(qc_list, out_path, config = NULL, ledger = NULL) {
  
  wb <- createWorkbook()
  
  # --- Sheet 1: Summary with Group Breakdown ---
  addWorksheet(wb, "Summary")
  
  # 1. Data Preparation & Calculation
  
  # Extract Base Counts per Subgroup
  final_counts <- if(!is.null(qc_list$breakdown_final)) as.numeric(qc_list$breakdown_final) else numeric(0)
  if(!is.null(qc_list$breakdown_final)) names(final_counts) <- names(qc_list$breakdown_final)
  
  # Dropped QC Counts
  drop_qc_counts <- numeric(0)
  if (!is.null(qc_list$dropped_rows_detail) && nrow(qc_list$dropped_rows_detail) > 0) {
    tbl <- table(qc_list$dropped_rows_detail$Original_Source)
    drop_qc_counts <- as.numeric(tbl)
    names(drop_qc_counts) <- names(tbl)
  }
  
  # Dropped A Priori Counts
  drop_pre_counts <- numeric(0)
  if (!is.null(qc_list$dropped_samples_apriori) && nrow(qc_list$dropped_samples_apriori) > 0) {
    tbl <- table(qc_list$dropped_samples_apriori$Original_Source)
    drop_pre_counts <- as.numeric(tbl)
    names(drop_pre_counts) <- names(tbl)
  }
  
  # Calculate Derived Metrics per Subgroup
  all_subgroups <- unique(c(names(final_counts), names(drop_qc_counts), names(drop_pre_counts)))
  all_subgroups <- sort(all_subgroups)
  
  subgroup_stats <- list()
  for(sg in all_subgroups) {
    v_fin <- if(sg %in% names(final_counts)) final_counts[[sg]] else 0
    v_qc  <- if(sg %in% names(drop_qc_counts)) drop_qc_counts[[sg]] else 0
    v_pre <- if(sg %in% names(drop_pre_counts)) drop_pre_counts[[sg]] else 0
    # `Initial` is RECONSTRUCTED as final + dropped. That is only correct if `breakdown_final`
    # is the true end state — it is not: run_qc_pipeline() captures it after the missingness
    # filter, while the PCA outliers are appended to dropped_rows_detail afterwards in
    # src/01. The result double-counts and this sheet reported Initial 96 / Final 87 for a
    # 94-row file whose analytic cohort is 82. The ledger (see `ledger_*`) is the authority
    # when present; the reconstruction survives only as a fallback for callers without one.
    v_init <- v_fin + v_qc + v_pre

    subgroup_stats[[sg]] <- list(Initial = v_init, DropPre = v_pre, DropQC = v_qc, Final = v_fin)
  }
  
  # Helper to aggregate metrics for a list of subgroups
  aggregate_metrics <- function(target_subgroups) {
    sum_init <- 0
    sum_drop <- 0
    sum_fin  <- 0
    
    for(sg in target_subgroups) {
      if (sg %in% names(subgroup_stats)) {
        s <- subgroup_stats[[sg]]
        sum_init <- sum_init + s$Initial
        sum_drop <- sum_drop - (s$DropPre + s$DropQC) # Sum and ensure negative for display
        sum_fin  <- sum_fin + s$Final
      }
    }
    return(c(sum_init, sum_drop, sum_fin))
  }
  
  # Grand Totals
  all_subs_list <- names(subgroup_stats)
  grand_totals <- aggregate_metrics(all_subs_list)
  total_init <- grand_totals[1]
  total_dropped_combined <- grand_totals[2]
  total_final <- grand_totals[3]
  total_pre <- sum(sapply(subgroup_stats, function(x) x$DropPre))
  total_qc  <- sum(sapply(subgroup_stats, function(x) x$DropQC))

  # Ledger overrides the reconstructed grand totals: it is recorded AT each filtering
  # point rather than inferred after the fact, so it cannot disagree with the data.
  # NB the per-subgroup columns necessarily describe only the patients whose outcome COULD
  # be mapped — a patient dropped for an unmappable outcome has no subgroup by definition.
  # So the total is set to the ELIGIBLE cohort (keeping this sheet internally consistent)
  # and the screened-vs-eligible step is reported separately below. The full flow lives in
  # the Cohort_Ledger sheet, which is the authority.
  n_screened <- NA_integer_
  n_unmappable <- 0L
  if (!is.null(ledger) && nrow(ledger)) {
    elig <- ledger$stage == "outcome_eligibility"
    n_screened   <- ledger$n_in[1]
    n_unmappable <- sum(ledger$n_dropped[elig])
    total_init   <- if (any(elig)) ledger$n_out[which(elig)[1]] else ledger$n_in[1]
    total_final  <- ledger$n_out[nrow(ledger)]
    total_pre    <- 0L
    total_qc     <- total_init - total_final
    total_dropped_combined <- -(total_pre + total_qc)
  }
  
  # 2. Build Tables
  
  # --- Table 1: Detailed Dropped Metrics (Subgroup Level) ---
  metrics_detailed <- c("Initial Samples", "Dropped Samples (A Priori)", "Dropped Samples (QC)", "Final Samples")
  totals_vector_det <- c(total_init, -total_pre, -total_qc, total_final)

  df_detailed <- data.frame(Metric = metrics_detailed, Total = totals_vector_det, stringsAsFactors = FALSE)
  for (sg in all_subgroups) {
    s <- subgroup_stats[[sg]]
    df_detailed[[sg]] <- c(s$Initial, -s$DropPre, -s$DropQC, s$Final)
  }
  # Screened -> eligible, prepended so the sheet starts where the input file starts. These
  # patients have no subgroup (their outcome could not be mapped), hence Total-only.
  if (is.finite(n_screened) && n_unmappable > 0) {
    pre <- data.frame(Metric = c("Screened (input file)", "Excluded: outcome not mappable"),
                      Total = c(n_screened, -n_unmappable), stringsAsFactors = FALSE)
    for (sg in all_subgroups) pre[[sg]] <- c(NA_real_, NA_real_)
    df_detailed <- rbind(pre, df_detailed)
  }
  
  # --- Table 2: By Group Dropping Metrics (Parent Level) ---
  if(!is.null(qc_list$group_mapping)) {
    mapping <- qc_list$group_mapping
    unique_parents <- sort(unique(mapping$Group))
  } else {
    mapping <- data.frame(Subgroup = all_subgroups, Group = all_subgroups, stringsAsFactors=FALSE)
    unique_parents <- all_subgroups
  }
  
  metrics_summary <- c("Initial", "Dropped Samples", "Final")
  df_parent <- data.frame(Metric = metrics_summary, Total = c(total_init, total_dropped_combined, total_final), stringsAsFactors = FALSE)
  
  for (parent in unique_parents) {
    subs_in_parent <- mapping$Subgroup[mapping$Group == parent]
    col_name <- paste0(parent, "_tot")
    df_parent[[col_name]] <- aggregate_metrics(subs_in_parent)
  }
  
  # --- Table 3: By Clinical Status (EP vs LS) ---
  ep_subgroups <- grep("_EP$", all_subgroups, value = TRUE)
  ls_subgroups <- grep("_LS$", all_subgroups, value = TRUE)
  
  # Define logic for controls (assuming "Healthy" or "Control" in name)
  is_control_func <- function(g_name) grepl("Healthy|Control", g_name, ignore.case = TRUE)
  ctrl_subgroups <- all_subgroups[sapply(all_subgroups, is_control_func)]
  
  df_clinical <- data.frame(Metric = metrics_summary, Total = c(total_init, total_dropped_combined, total_final), stringsAsFactors = FALSE)
  
  # Column: Healthy_Donors_tot (or similar control group)
  if(length(ctrl_subgroups) > 0) {
    # We use the parent name from mapping if possible for cleaner header, otherwise "Control_tot"
    ctrl_header <- "Healthy_Donors_tot" 
    df_clinical[[ctrl_header]] <- aggregate_metrics(ctrl_subgroups)
  }
  
  # Column: EP_tot
  if(length(ep_subgroups) > 0) {
    df_clinical[["EP_tot"]] <- aggregate_metrics(ep_subgroups)
  }
  
  # Column: LS_tot
  if(length(ls_subgroups) > 0) {
    df_clinical[["LS_tot"]] <- aggregate_metrics(ls_subgroups)
  }
  
  # --- Table 4: By Group Dropping Metrics (Aggregated Control vs Case) ---
  ctrl_parents <- unique_parents[sapply(unique_parents, is_control_func)]
  case_parents <- unique_parents[!sapply(unique_parents, is_control_func)]
  
  df_macro <- data.frame(Metric = metrics_summary, Total = c(total_init, total_dropped_combined, total_final), stringsAsFactors = FALSE)
  
  # Aggregated Control Column
  if(length(ctrl_parents) > 0) {
    all_ctrl_subs <- mapping$Subgroup[mapping$Group %in% ctrl_parents]
    # Apply suffix to each parent individually before collapsing
    col_name <- paste0(paste0(ctrl_parents, "_tot"), collapse = " + ")
    df_macro[[col_name]] <- aggregate_metrics(all_ctrl_subs)
  }
  
  # Aggregated Case Column
  if(length(case_parents) > 0) {
    all_case_subs <- mapping$Subgroup[mapping$Group %in% case_parents]
    # Apply suffix to each parent individually before collapsing
    col_name <- paste0(paste0(case_parents, "_tot"), collapse = " + ")
    df_macro[[col_name]] <- aggregate_metrics(all_case_subs)
  }
  
  # 3. Writing to Excel
  curr_row <- 1
  
  # Section 1: By Group Dropping Metrics (Parent Level)
  writeData(wb, "Summary", "BY GROUP DROPPING METRICS:", startRow = curr_row)
  addStyle(wb, "Summary", createStyle(textDecoration = "bold"), rows = curr_row, cols = 1)
  curr_row <- curr_row + 1
  
  writeData(wb, "Summary", df_parent, startRow = curr_row)
  addStyle(wb, "Summary", createStyle(textDecoration = "bold"), rows = curr_row, cols = 1:ncol(df_parent), gridExpand = TRUE)
  
  curr_row <- curr_row + nrow(df_parent) + 2
  
  # Section 2: By Clinical Status (EP vs LS)
  writeData(wb, "Summary", "BY CLINICAL STATUS METRICS (EP vs LS):", startRow = curr_row)
  addStyle(wb, "Summary", createStyle(textDecoration = "bold"), rows = curr_row, cols = 1)
  curr_row <- curr_row + 1
  
  writeData(wb, "Summary", df_clinical, startRow = curr_row)
  addStyle(wb, "Summary", createStyle(textDecoration = "bold"), rows = curr_row, cols = 1:ncol(df_clinical), gridExpand = TRUE)
  
  curr_row <- curr_row + nrow(df_clinical) + 2
  
  # Section 3: By Group Dropping Metrics (Aggregated)
  writeData(wb, "Summary", "BY GROUP DROPPING METRICS (AGGREGATED):", startRow = curr_row)
  addStyle(wb, "Summary", createStyle(textDecoration = "bold"), rows = curr_row, cols = 1)
  curr_row <- curr_row + 1
  
  writeData(wb, "Summary", df_macro, startRow = curr_row)
  addStyle(wb, "Summary", createStyle(textDecoration = "bold"), rows = curr_row, cols = 1:ncol(df_macro), gridExpand = TRUE)
  
  curr_row <- curr_row + nrow(df_macro) + 2
  
  # Section 4: Detailed Dropped Metrics
  writeData(wb, "Summary", "DETAILED DROPPED METRICS:", startRow = curr_row)
  addStyle(wb, "Summary", createStyle(textDecoration = "bold"), rows = curr_row, cols = 1)
  curr_row <- curr_row + 1
  
  writeData(wb, "Summary", df_detailed, startRow = curr_row)
  addStyle(wb, "Summary", createStyle(textDecoration = "bold"), rows = curr_row, cols = 1:ncol(df_detailed), gridExpand = TRUE)
  
  curr_row <- curr_row + nrow(df_detailed) + 2 
  
  # Section 5: Markers Metrics
  n_markers_apriori <- if(!is.null(qc_list$dropped_markers_apriori)) nrow(qc_list$dropped_markers_apriori) else 0
  totals_markers <- c(qc_list$n_col_init,
                      n_markers_apriori,
                      qc_list$n_col_dropped,
                      qc_list$n_col_zerovar,
                      qc_list$n_col_final)
  metrics_markers <- c("Initial Markers", "Dropped Markers (A Priori)", "Dropped Markers (High NA)", "Dropped Markers (Zero Var)", "Final Markers")
  df_markers <- data.frame(Metric = metrics_markers, Total = totals_markers, stringsAsFactors = FALSE)
  
  writeData(wb, "Summary", "MARKERS METRICS:", startRow = curr_row)
  addStyle(wb, "Summary", createStyle(textDecoration = "bold"), rows = curr_row, cols = 1)
  curr_row <- curr_row + 1
  writeData(wb, "Summary", df_markers, startRow = curr_row)
  addStyle(wb, "Summary", createStyle(textDecoration = "bold"), rows = curr_row, cols = 1:ncol(df_markers), gridExpand = TRUE)
  curr_row <- curr_row + nrow(df_markers) + 2 
  
  # Section 6: Final Markers List
  if (!is.null(qc_list$final_markers_names) && length(qc_list$final_markers_names) > 0) {
    final_mks <- qc_list$final_markers_names
    formatted_df <- NULL
    
    if (!is.null(config) && !is.null(config$qc_reporting$marker_categories)) {
      # Categorized Logic
      cats <- config$qc_reporting$marker_categories
      df_list <- list()
      categorized_markers <- c()
      
      for (cat_name in names(cats)) {
        expected <- cats[[cat_name]]
        present <- intersect(expected, final_mks)
        categorized_markers <- c(categorized_markers, present)
        # Force column to have at least one valid entry or empty string to preserve structure
        df_list[[cat_name]] <- if (length(present) > 0) present else character(0)
      }
      
      # Catch-all for markers not in categories
      others <- setdiff(final_mks, categorized_markers)
      if (length(others) > 0) df_list[["Others"]] <- others
      
      # Pad with empty strings to make a data.frame
      max_len <- max(sapply(df_list, length))
      if (max_len > 0) {
        formatted_df <- data.frame(lapply(df_list, function(x) {
          c(x, rep("", max_len - length(x)))
        }), stringsAsFactors = FALSE)
      }
    } else {
      # Default 5-column Grid Logic
      n_per_row <- 5
      n_needed <- ceiling(length(final_mks) / n_per_row) * n_per_row
      mks_padded <- c(final_mks, rep("", n_needed - length(final_mks)))
      formatted_df <- as.data.frame(matrix(mks_padded, ncol = n_per_row, byrow = TRUE))
      colnames(formatted_df) <- paste0("Col_", 1:n_per_row)
    }
    
    if (!is.null(formatted_df)) {
      writeData(wb, "Summary", "FINAL MARKERS LIST:", startRow = curr_row)
      addStyle(wb, "Summary", createStyle(textDecoration = "bold"), rows = curr_row, cols = 1)
      # Write data without column names (lineage headers)
      writeData(wb, "Summary", formatted_df, startRow = curr_row + 1, colNames = FALSE)
    }
  }
  
  # --- Sheet 2: Details Dropped ---
  addWorksheet(wb, "Details_Dropped")
  curr_row_det <- 1
  dropped_patients_all <- data.frame()
  if (!is.null(qc_list$dropped_samples_apriori) && nrow(qc_list$dropped_samples_apriori) > 0) {
    dropped_patients_all <- rbind(dropped_patients_all, qc_list$dropped_samples_apriori)
  }
  if (!is.null(qc_list$dropped_rows_detail) && nrow(qc_list$dropped_rows_detail) > 0) {
    dropped_patients_all <- rbind(dropped_patients_all, qc_list$dropped_rows_detail)
  }
  
  if (nrow(dropped_patients_all) > 0) {
    # Ensure nice column ordering
    cols_order <- intersect(c("Patient_ID", "NA_Percent", "Reason", "Original_Source"), names(dropped_patients_all))
    cols_order <- c(cols_order, setdiff(names(dropped_patients_all), cols_order))
    dropped_patients_all <- dropped_patients_all[, cols_order, drop=FALSE]
    
    writeData(wb, "Details_Dropped", "Dropped Samples:", startRow = curr_row_det)
    addStyle(wb, "Details_Dropped", createStyle(textDecoration = "bold"), rows = curr_row_det, cols = 1)
    writeData(wb, "Details_Dropped", dropped_patients_all, startRow = curr_row_det + 1)
    curr_row_det <- curr_row_det + nrow(dropped_patients_all) + 3
  }
  
  if (!is.null(qc_list$dropped_markers_apriori) && nrow(qc_list$dropped_markers_apriori) > 0) {
    writeData(wb, "Details_Dropped", "Markers Excluded (A Priori):", startRow = curr_row_det)
    addStyle(wb, "Details_Dropped", createStyle(textDecoration = "bold"), rows = curr_row_det, cols = 1)
    writeData(wb, "Details_Dropped", qc_list$dropped_markers_apriori, startRow = curr_row_det + 1)
    curr_row_det <- curr_row_det + nrow(qc_list$dropped_markers_apriori) + 3
  }
  
  if (nrow(qc_list$dropped_cols_detail) > 0) {
    writeData(wb, "Details_Dropped", "Dropped Markers (QC):", startRow = curr_row_det)
    addStyle(wb, "Details_Dropped", createStyle(textDecoration = "bold"), rows = curr_row_det, cols = 1)
    writeData(wb, "Details_Dropped", qc_list$dropped_cols_detail, startRow = curr_row_det + 1)
  }
  
  # --- Sheet 3: Details Imputed (New addition) ---
  if (!is.null(qc_list$imputed_details) && nrow(qc_list$imputed_details) > 0) {
    addWorksheet(wb, "Details_Imputed")
    curr_row_imp <- 1
    
    writeData(wb, "Details_Imputed", "Samples requiring Imputation:", startRow = curr_row_imp)
    addStyle(wb, "Details_Imputed", createStyle(textDecoration = "bold"), rows = curr_row_imp, cols = 1)
    
    # Write Main Table
    writeData(wb, "Details_Imputed", qc_list$imputed_details, startRow = curr_row_imp + 1)
    
    # Calculate Total Row Position (Header + Data + 2 blank lines)
    total_row_idx <- curr_row_imp + 1 + nrow(qc_list$imputed_details) + 2
    
    # Write Total Summary
    total_msg <- paste("Total imputed samples:", nrow(qc_list$imputed_details))
    writeData(wb, "Details_Imputed", total_msg, startRow = total_row_idx)
    addStyle(wb, "Details_Imputed", createStyle(textDecoration = "bold"), rows = total_row_idx, cols = 1)
  }
  
  # --- Sheet: Cohort_Ledger (authoritative patient flow, one row per filtering point) ---
  if (!is.null(ledger) && nrow(ledger)) {
    addWorksheet(wb, "Cohort_Ledger")
    writeData(wb, "Cohort_Ledger",
              "Patient flow recorded AT each filtering point (not reconstructed).", startRow = 1)
    addStyle(wb, "Cohort_Ledger", createStyle(textDecoration = "bold"), rows = 1, cols = 1)
    writeData(wb, "Cohort_Ledger", ledger, startRow = 3)
    setColWidths(wb, "Cohort_Ledger", cols = seq_len(ncol(ledger)), widths = "auto")
  }

  saveWorkbook(wb, out_path, overwrite = TRUE)
  message(sprintf("[IO] QC Report saved to: %s", out_path))
}

#' @title Safe Excel Sheet Name
#' @description Truncates a string to comply with Excel's 31-character limit for sheet names.
#' @param name String to truncate.
#' @param max_len Integer. Maximum character length (default 31).
#' @return Truncated string safe for openxlsx.
safe_sheet_name <- function(name, max_len = 31) {
  return(substr(name, 1, max_len))
}


#' @title Validate Pipeline Configuration
#' @description
#' Pre-flight check that the loaded YAML config has the required structural keys
#' before any expensive computation begins. Checks invariants only (not every
#' parameter leaf), so the function stays stable as optional keys are added.
#'
#' Stops with an actionable message if a required key is missing or empty.
#'
#' @param cfg Named list. The result of \code{yaml::read_yaml()}.
#' @return Invisibly TRUE if all checks pass.
validate_config <- function(cfg) {
  if (is.null(cfg$output_root))
    stop("[CONFIG] Missing required key: output_root")

  if (is.null(cfg$features$facs) || length(cfg$features$facs) == 0)
    stop("[CONFIG] features.facs is empty or missing — define the FACS marker panel")

  if (is.null(cfg$features$soluble) || length(cfg$features$soluble) == 0)
    stop("[CONFIG] features.soluble is empty or missing — define the soluble marker panel")

  if (!is.null(cfg$experiments) && length(cfg$experiments) > 0) {
    for (nm in names(cfg$experiments)) {
      exp <- cfg$experiments[[nm]]
      # Disabled experiments (enabled: false) are skipped by main.R; their
      # columns need not exist, so don't validate them.
      if (!is.null(exp$enabled) && !isTRUE(exp$enabled)) next
      if (is.null(exp$clinical$target_column))
        stop(sprintf("[CONFIG] Experiment '%s' missing clinical.target_column", nm))
      validate_experiment_columns(cfg, exp, nm)
    }
  } else {
    if (is.null(cfg$clinical$target_column))
      stop("[CONFIG] Flat config missing clinical.target_column")
  }

  invisible(TRUE)
}

#' @title Read only the column names of an Excel file
#' @description Header-only read used by pre-flight config validation. Returns an
#'   empty character vector (never errors) when the file is absent or unreadable,
#'   so validation degrades to a no-op rather than crashing.
#' @param path Excel file path (resolved relative to the project root via here()).
#' @return Character vector of column names, or character(0) on any failure.
read_excel_header <- function(path) {
  if (is.null(path)) return(character(0))
  abs <- if (file.exists(path)) path else here::here(path)
  if (!file.exists(abs)) return(character(0))
  tryCatch(names(readxl::read_excel(abs, n_max = 0)),
           error = function(e) character(0))
}

#' @title Validate that an experiment's referenced columns exist in its input
#' @description Pre-flight, per-experiment column check against the T0 input
#'   header. Stops on a missing \code{clinical.target_column} or any missing
#'   \code{machine_learning.clinical_model.vars} column (the standard-of-care
#'   predictors); warns on missing feature markers (the pipeline already
#'   tolerates these with an in-step warning). Skipped silently when the input
#'   file is absent/unreadable (file existence is enforced later, per pass).
#' @param cfg The full config (for global fallbacks).
#' @param exp The per-experiment block.
#' @param nm Experiment name (for messages).
#' @return Invisibly TRUE.
validate_experiment_columns <- function(cfg, exp, nm) {
  t0 <- if (!is.null(exp$input_file_t0)) exp$input_file_t0 else cfg$input_file_t0
  hdr <- read_excel_header(t0)
  if (length(hdr) == 0) return(invisible(TRUE))  # can't read header → skip

  tc <- exp$clinical$target_column
  if (!is.null(tc) && !(tc %in% hdr))
    stop(sprintf("[CONFIG] Experiment '%s': clinical.target_column '%s' not found in %s",
                 nm, tc, t0))

  cm_vars <- exp$machine_learning$clinical_model$vars
  if (!is.null(cm_vars)) {
    cols <- unlist(cm_vars, use.names = FALSE)
    miss <- setdiff(cols, hdr)
    if (length(miss) > 0)
      stop(sprintf("[CONFIG] Experiment '%s': clinical_model.vars column(s) not found in %s: %s",
                   nm, t0, paste(miss, collapse = ", ")))
  }

  facs <- if (!is.null(exp$features$facs)) exp$features$facs else cfg$features$facs
  sol  <- if (!is.null(exp$features$soluble)) exp$features$soluble else cfg$features$soluble
  markers <- unique(c(unlist(facs, use.names = FALSE), unlist(sol, use.names = FALSE)))
  miss_mk <- setdiff(markers, hdr)
  if (length(miss_mk) > 0)
    warning(sprintf("[CONFIG] Experiment '%s': %d feature marker(s) not in %s: %s",
                    nm, length(miss_mk), t0, paste(head(miss_mk, 10), collapse = ", ")))

  invisible(TRUE)
}

#' @title Build an Isolated Per-Pass Configuration
#' @description
#' Encapsulates the null-guarded config assembly repeated once per pass in
#' \code{main.R} (and the single-step diagnostics runners). Starts from a copy of
#' \code{base_config}, layers the per-experiment overrides relevant to the given
#' pass, points \code{output_root} at the experiment's run subdirectory (creating
#' it), and returns the isolated config. Mutations never leak back into
#' \code{base_config} or sibling passes.
#'
#' Pass-specific behaviour (mirrors the original inline blocks exactly):
#'   - \strong{standard}: forces \code{is_longitudinal = FALSE}; applies
#'     \code{input_file}; asserts non-empty \code{features.facs}.
#'   - \strong{longitudinal}: forces \code{is_longitudinal = TRUE} and
#'     \code{qc$remove_outliers = FALSE} (paired LMM needs both timepoints);
#'     applies \code{input_file_t0}/\code{input_file_t1}.
#'   - \strong{machine_learning}: applies \code{input_file}/\code{input_file_t0}
#'     and merges \code{machine_learning} overrides via \code{modifyList}.
#'
#' @param base_config The global configuration (single source of truth).
#' @param exp_cfg The per-experiment override block.
#' @param exp_name Experiment name (becomes \code{project_name}).
#' @param pass_mode One of "standard", "longitudinal", "machine_learning".
#' @param run_root The run root directory; \code{output_root} becomes
#'   \code{run_root/exp_name}.
#' @return The isolated, pass-ready config list.
#' Resolve the NLR toggle: env INCLUDE_NLR overrides the config default.
#'
#' Single source of truth for the with-/without-NLR variants. `main.R` needs it early
#' (it labels the run dir with_nlr/no_nlr), and `prepare_pass_config()` re-applies it so
#' EVERY consumer — including the `diagnostics/` single-step runners, which read the YAML
#' directly and never saw main.R's resolution — honours the env var. Idempotent: main.R
#' writes the resolved value back into base_config, so re-resolving yields the same answer.
#' Without this the runners silently ran the WITH-NLR model regardless of INCLUDE_NLR, and
#' would happily write it into a `*_no_nlr` run dir — a mislabeled artifact with no error.
#' @return TRUE/FALSE
resolve_include_nlr <- function(base_config) {
  e <- Sys.getenv("INCLUDE_NLR", unset = NA_character_)
  if (!is.na(e) && nzchar(e)) tolower(e) %in% c("1", "true", "yes", "t")
  else !isFALSE(base_config$include_nlr)
}

prepare_pass_config <- function(base_config, exp_cfg, exp_name, pass_mode, run_root) {
  valid_modes <- c("standard", "longitudinal", "machine_learning")
  if (!pass_mode %in% valid_modes)
    stop(sprintf("[CONFIG] Unknown pass_mode '%s' (expected one of: %s)",
                 pass_mode, paste(valid_modes, collapse = ", ")))

  config <- base_config
  config$project_name <- exp_name
  config$run_mode     <- pass_mode
  # Re-apply the NLR toggle here so a runner that never called main.R's resolution
  # cannot silently produce a with-NLR model under INCLUDE_NLR=false (Step 06 reads
  # config$include_nlr). No-op when the caller already resolved it.
  config$include_nlr  <- resolve_include_nlr(base_config)

  # Clinical override (all passes)
  if (!is.null(exp_cfg$clinical)) config$clinical <- exp_cfg$clinical

  # Feature panel override (all passes). soluble may legitimately be empty
  # (FACS-only panels); the standard pass additionally asserts non-empty FACS.
  if (!is.null(exp_cfg$features)) {
    config$features <- exp_cfg$features
    if (pass_mode == "standard" && length(unlist(config$features$facs)) == 0)
      stop(sprintf("[CONFIG] Experiment '%s': features.facs is empty", exp_name))
  }

  if (pass_mode == "standard") {
    config$is_longitudinal <- FALSE
    if (!is.null(exp_cfg$input_file)) config$input_file <- exp_cfg$input_file

  } else if (pass_mode == "longitudinal") {
    config$is_longitudinal <- TRUE
    # Critical structural rule: disable multivariate outliers to prevent temporal censoring.
    config$qc$remove_outliers <- FALSE
    if (!is.null(exp_cfg$input_file_t0)) config$input_file_t0 <- exp_cfg$input_file_t0
    if (!is.null(exp_cfg$input_file_t1)) config$input_file_t1 <- exp_cfg$input_file_t1

  } else if (pass_mode == "machine_learning") {
    if (!is.null(exp_cfg$input_file))    config$input_file    <- exp_cfg$input_file
    if (!is.null(exp_cfg$input_file_t0)) config$input_file_t0 <- exp_cfg$input_file_t0
    # Step 06 reads T1 via the Step 01 longitudinal RDS, not this path — but the value is
    # frozen into locked_model$development so a Δ scorer can name its own T1 cohort file.
    # Without the override it would silently inherit the top-level default (a DIFFERENT
    # extraction), and anonymized IDs are not stable across extractions.
    if (!is.null(exp_cfg$input_file_t1)) config$input_file_t1 <- exp_cfg$input_file_t1
    if (!is.null(exp_cfg$machine_learning)) {
      config$machine_learning <- modifyList(
        if (!is.null(config$machine_learning)) config$machine_learning else list(),
        exp_cfg$machine_learning
      )
    }
  }

  config$output_root <- file.path(run_root, config$project_name)
  if (!dir.exists(config$output_root)) dir.create(config$output_root, recursive = TRUE)

  config
}