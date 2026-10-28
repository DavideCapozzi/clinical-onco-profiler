# R/modules_posthoc_duplicates.R
# ==============================================================================
# DUPLICATED-ROW SENSITIVITY (step 07, duplicate_rows/; formerly diagnostics/diag_55)
# In the final inputs the T1 rows of NSCLC_10 (PR) and NSCLC_12 (SD) are identical on
# every FACS column (two different persons), so at most one of them is that patient's
# sample. The dataset is NOT changed; this drops the T1 sample of each listed ID alone
# and of all of them together, and repeats
#   [1] the step-04 gate (same model as fit_feature_lmm, BH over the panel)
#   [2] the added-value layer with the gate held fixed (CV AUC, LRT), Z recomputed on the
#       reduced longitudinal frame exactly as step 01 fast mode does (scale() of the raw logit).
# A dropped T1 sample makes that patient T0-only, so it leaves the paired delta frame.
# The selection-aware p is not repeated (it needs a full null per scenario).
# Dependencies: modules_longitudinal.R (fit_feature_lmm), modules_ml_utility.R (av_*)
# ==============================================================================

#' @title Duplicated-T1-row sensitivity of the gate and the added-value layer
#' @param DATA_T0,DATA_LONG Step-01 standard / longitudinal objects.
#' @param gate The step-06 gate (composite markers), held fixed.
#' @param drop_t1 Patient IDs whose T1 row is dropped (config post_hoc.duplicate_rows.drop_t1).
#' @param config Pass config (clinical labels, input_file_t0).
#' @param clinical_vars Clinical-model columns (after the NLR toggle).
#' @param seed stats.seed (cv_kfold uses seed + 977, as step 06).
#' @param say sprintf-style logger.
#' @return data.frame, one row per scenario (and per gate marker in long columns).
run_duplicate_row_sensitivity <- function(DATA_T0, DATA_LONG, gate, drop_t1, config,
                                          clinical_vars, seed = 2026L, say = message) {
  DL <- DATA_LONG
  resp_label <- config$clinical$responder_label
  R <- av_context(DATA_T0, DATA_LONG, resp_label, config$input_file_t0, clinical_vars)

  mk <- DL$hybrid_markers
  raw <- as.matrix(DL$hybrid_data_raw[, mk])
  ok_z <- max(abs(scale(raw) - as.matrix(DL$hybrid_data_z[, mk])), na.rm = TRUE)
  say("fast-mode Z = scale(raw logit) on the full frame: max |diff| %.1e", ok_z)
  if (!(ok_z < 1e-10)) stop("[FATAL] duplicate rows: longitudinal Z is not scale(raw logit); Z cannot be recomputed on a reduced frame")
  rows <- lapply(drop_t1, function(id) which(DL$metadata$Patient_ID == id & DL$metadata$Timepoint == "T1"))
  names(rows) <- drop_t1
  miss <- names(rows)[lengths(rows) != 1]
  if (length(miss)) stop(sprintf("[FATAL] duplicate rows: no unique T1 row for %s", paste(miss, collapse = ", ")))
  if (length(drop_t1) == 2)
    say("T1 rows of %s and %s identical on all %d markers (raw logit): %s", drop_t1[1], drop_t1[2], length(mk),
        isTRUE(all.equal(unname(raw[rows[[1]], ]), unname(raw[rows[[2]], ]))))

  drop_rows <- function(r) {
    if (!length(r)) return(DL)
    keep <- setdiff(seq_len(nrow(raw)), r)
    D <- DL
    D$metadata <- DL$metadata[keep, , drop = FALSE]
    D$hybrid_data_raw <- DL$hybrid_data_raw[keep, , drop = FALSE]
    z <- scale(raw[keep, , drop = FALSE]); attr(z, "scaled:center") <- NULL; attr(z, "scaled:scale") <- NULL
    D$hybrid_data_z <- z
    D
  }
  # Scenarios: as is, each listed ID alone, all of them together (IDs per scenario).
  SC_ids <- c(list("as is" = character(0)),
              setNames(as.list(drop_t1), sprintf("%s T1 dropped", drop_t1)))
  if (length(drop_t1) > 1) SC_ids[[if (length(drop_t1) == 2) "both T1 dropped" else "all T1 dropped"]] <- drop_t1
  SC <- lapply(SC_ids, function(ids) unlist(rows[ids], use.names = FALSE))

  out <- list()
  for (nm in names(SC)) {
    D <- drop_rows(SC[[nm]])
    df <- cbind(D$metadata, as.data.frame(D$hybrid_data_raw[, mk]))
    df$Timepoint <- factor(df$Timepoint, levels = c("T0", "T1"))
    df$Group <- factor(df$Group, levels = c(config$clinical$non_responder_label, resp_label))
    df$Patient_ID <- as.factor(df$Patient_ID)
    P <- do.call(rbind, lapply(mk, function(m) fit_feature_lmm(df, m)))
    P$FDR <- p.adjust(P$P_Value_Interaction, "BH"); P <- P[order(P$P_Value_Interaction), ]
    R2 <- R; R2$DATA_LONG <- D
    A <- av_assemble(R2, gate, "delta"); cv <- av_cv_kfold(A, seed = seed + 977L); lr <- av_lrt(A)
    g <- P[P$Marker %in% gate, ]
    k <- length(gate)
    say("\n%s: LMM %d pts / %d obs | FDR<0.05: %d (%s) | top-%d: %s", nm, length(unique(df$Patient_ID)), nrow(df),
        sum(P$FDR < 0.05), paste(P$Marker[P$FDR < 0.05], collapse = ","), k, paste(P$Marker[seq_len(k)], collapse = ","))
    for (m in gate) say("    %-10s beta %+.3f  p %.5f  FDR %.4f", m, g$Estimate_Interaction[g$Marker == m],
                        g$P_Value_Interaction[g$Marker == m], g$FDR[g$Marker == m])
    say("    next: %s FDR %.4f", P$Marker[k + 1], P$FDR[k + 1])
    say("    delta frame n %d (%d PR / %d SD/PD) | CV AUC clinical %.3f immune %.3f combined %.3f | dAUC %+.3f | LRT chi2 %.2f (naive p %.2e)",
        length(A$yb), sum(A$yb), sum(1 - A$yb), cv$auc_clinical, cv$auc_immune, cv$auc_combined, cv$delta_auc, lr$stat, lr$p)
    row <- data.frame(scenario = nm, dropped_t1 = paste(SC_ids[[nm]], collapse = ","),
                      lmm_patients = length(unique(df$Patient_ID)), lmm_obs = nrow(df),
                      n_fdr = sum(P$FDR < 0.05), fdr_markers = paste(P$Marker[P$FDR < 0.05], collapse = ","),
                      topk_markers = paste(P$Marker[seq_len(k)], collapse = ","),
                      next_marker = P$Marker[k + 1], next_FDR = P$FDR[k + 1],
                      delta_n = length(A$yb), delta_pos = sum(A$yb), delta_neg = sum(1 - A$yb),
                      cv_auc_clinical = cv$auc_clinical, cv_auc_immune = cv$auc_immune,
                      cv_auc_combined = cv$auc_combined, cv_delta_auc = cv$delta_auc,
                      lrt_chi2 = lr$stat, lrt_naive_p = lr$p, stringsAsFactors = FALSE)
    for (m in gate) {
      row[[paste0(m, "_beta")]] <- g$Estimate_Interaction[g$Marker == m]
      row[[paste0(m, "_p")]]    <- g$P_Value_Interaction[g$Marker == m]
      row[[paste0(m, "_FDR")]]  <- g$FDR[g$Marker == m]
    }
    out[[nm]] <- row
  }
  res <- do.call(rbind, out); rownames(res) <- NULL
  res
}
