# R/modules_split.R
# ==============================================================================
# COHORT SPLIT MODULE (Discovery / Validation)
# Description: Pass-independent temporal partition of the baseline (T0) cohort
#   into a discovery set (the earliest n_discovery patients passing missingness
#   QC) and a held-out validation set (the remainder). The split is defined once
#   on the T0 QC-A pool — the basic-missingness filter that is identical across
#   the standard and longitudinal passes — so every downstream pass filters to
#   the same discovery patients and the validation set is never seen in training.
# ==============================================================================

#' @title Compute Discovery/Validation Cohort Split
#' @description Mirrors Step 01's clinical filter + basic-missingness QC on the T0
#'   matrix to define a stable QC pool, then ranks survivors by an enrollment-date
#'   column and assigns the earliest \code{n_discovery} to discovery and the rest
#'   to validation. QC here uses \code{remove_outliers = FALSE} (the pass-invariant
#'   QC-A stage); multivariate outlier removal happens later, inside the discovery
#'   run only, and never reassigns a patient to validation.
#' @param raw_t0 Data frame of the raw T0 cohort (first column = patient id).
#' @param config Experiment-resolved configuration (clinical, features, qc).
#' @param split_cfg The \code{validation_split} block: n_discovery, order_by, order_format.
#' @return list(discovery_ids, validation_ids, report); report is a data frame with
#'   Patient_ID, Group, enroll_date, rank, assignment.
compute_cohort_split <- function(raw_t0, config, split_cfg) {
  clin_col  <- config$clinical$target_column
  resp_lbl  <- config$clinical$responder_label
  nresp_lbl <- config$clinical$non_responder_label
  order_by  <- split_cfg$order_by
  order_fmt <- if (!is.null(split_cfg$order_format)) split_cfg$order_format else "%d/%m/%Y"
  n_disc    <- as.integer(split_cfg$n_discovery)

  if (is.null(order_by))  stop("[SPLIT] validation_split.order_by is required.")
  if (is.na(n_disc))      stop("[SPLIT] validation_split.n_discovery must be an integer.")

  raw <- as.data.frame(raw_t0, stringsAsFactors = FALSE)
  colnames(raw)[1] <- "Patient_ID"
  raw$Patient_ID <- as.character(raw$Patient_ID)

  if (!clin_col %in% colnames(raw))
    stop(sprintf("[SPLIT] target_column '%s' not found in T0 data.", clin_col))
  if (!order_by %in% colnames(raw))
    stop(sprintf("[SPLIT] order_by column '%s' not found in T0 data.", order_by))

  # 1. Clinical filter + group mapping (mirror of Step 01) -------------------
  if (!is.null(config$clinical$mapping)) {
    mapped_resp  <- as.character(unlist(config$clinical$mapping[[resp_lbl]]))
    mapped_nresp <- as.character(unlist(config$clinical$mapping[[nresp_lbl]]))
    raw$Group <- dplyr::case_when(
      as.character(raw[[clin_col]]) %in% mapped_resp  ~ resp_lbl,
      as.character(raw[[clin_col]]) %in% mapped_nresp ~ nresp_lbl,
      TRUE ~ NA_character_
    )
  } else {
    raw$Group <- ifelse(
      as.character(raw[[clin_col]]) %in% c(resp_lbl, nresp_lbl),
      as.character(raw[[clin_col]]), NA_character_
    )
  }
  raw <- raw[!is.na(raw$Group), , drop = FALSE]
  if (nrow(raw) == 0) stop("[SPLIT] No patients with valid clinical labels in T0.")

  # 2. Basic-missingness QC-A on the effective panel (remove_outliers = FALSE) -
  facs_feats <- unlist(config$features$facs)
  sol_feats  <- unlist(config$features$soluble)
  markers    <- intersect(unique(c(facs_feats, sol_feats)), colnames(raw))
  if (length(markers) == 0) stop("[SPLIT] No configured markers found in T0 data.")

  mat <- as.matrix(raw[, markers, drop = FALSE])
  rownames(mat) <- raw$Patient_ID

  qc_cfg <- config$qc
  qc_cfg$remove_outliers <- FALSE
  qc <- run_qc_pipeline(
    mat_raw          = mat,
    metadata         = data.frame(Patient_ID = raw$Patient_ID, Group = raw$Group,
                                  stringsAsFactors = FALSE),
    qc_config        = qc_cfg,
    facs_features    = intersect(facs_feats, markers),
    soluble_features = intersect(sol_feats, markers),
    stratification_col = "Group"
  )
  pool <- raw[match(rownames(qc$data), raw$Patient_ID), , drop = FALSE]

  # 3. Temporal rank + split -------------------------------------------------
  enroll <- as.Date(pool[[order_by]], format = order_fmt)
  if (any(is.na(enroll)))
    stop(sprintf("[SPLIT] %d enrollment date(s) failed to parse with format '%s' in column '%s'.",
                 sum(is.na(enroll)), order_fmt, order_by))
  if (n_disc < 1 || n_disc >= nrow(pool))
    stop(sprintf("[SPLIT] n_discovery=%d invalid for QC pool of %d patients (need 1..%d).",
                 n_disc, nrow(pool), nrow(pool) - 1L))

  ord    <- order(enroll, pool$Patient_ID)  # ties broken by id for determinism
  pool   <- pool[ord, , drop = FALSE]
  enroll <- enroll[ord]
  assignment <- c(rep("discovery", n_disc), rep("validation", nrow(pool) - n_disc))

  report <- data.frame(
    Patient_ID  = pool$Patient_ID,
    Group       = pool$Group,
    enroll_date = format(enroll, "%Y-%m-%d"),
    rank        = seq_len(nrow(pool)),
    assignment  = assignment,
    stringsAsFactors = FALSE
  )

  message(sprintf(
    "[SPLIT] QC pool=%d | discovery=%d (%s) | validation=%d (%s)",
    nrow(pool), n_disc,
    paste(sprintf("%s:%d", names(table(report$Group[assignment == "discovery"])),
                  table(report$Group[assignment == "discovery"])), collapse = " "),
    nrow(pool) - n_disc,
    paste(sprintf("%s:%d", names(table(report$Group[assignment == "validation"])),
                  table(report$Group[assignment == "validation"])), collapse = " ")
  ))

  list(
    discovery_ids  = report$Patient_ID[assignment == "discovery"],
    validation_ids = report$Patient_ID[assignment == "validation"],
    report         = report
  )
}

#' @title Materialize Discovery/Validation Excel Files
#' @description Writes row-subset Excel snapshots of a source cohort file, split by
#'   the discovery/validation patient ids. All original columns are preserved so the
#'   downstream pipeline reads the discovery file exactly as it would the full file.
#' @param src_path Path to the source Excel (T0 or T1).
#' @param timepoint Label "T0"/"T1" used in the output filename.
#' @param split Result of \code{compute_cohort_split()}.
#' @param out_dir Output directory.
#' @param exp_name Experiment name used in output filenames.
#' @return Invisibly, a named list with the discovery and validation file paths.
write_split_excels <- function(src_path, timepoint, split, out_dir, exp_name) {
  raw <- readxl::read_excel(src_path)
  ids <- as.character(raw[[1]])

  disc <- raw[ids %in% split$discovery_ids, , drop = FALSE]
  val  <- raw[ids %in% split$validation_ids, , drop = FALSE]

  disc_path <- file.path(out_dir, sprintf("%s_discovery_%s.xlsx",  exp_name, timepoint))
  val_path  <- file.path(out_dir, sprintf("%s_validation_%s.xlsx", exp_name, timepoint))
  openxlsx::write.xlsx(disc, disc_path, overwrite = TRUE)
  openxlsx::write.xlsx(val,  val_path,  overwrite = TRUE)

  message(sprintf("[SPLIT]   %s: discovery %d rows -> %s | validation %d rows -> %s",
                  timepoint, nrow(disc), basename(disc_path), nrow(val), basename(val_path)))
  invisible(list(discovery = disc_path, validation = val_path))
}
