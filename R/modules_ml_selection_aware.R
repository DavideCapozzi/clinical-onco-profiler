# R/modules_ml_selection_aware.R
# ==============================================================================
# SELECTION-AWARE LRT — permutation null for the immune increment that replays
# the gate SELECTION, not only the coefficient (formerly diagnostics/diag_44).
#
# The step-06 gate-fixed permutation p (increment$lrt_perm_p) shuffles labels while
# holding the pre-specified composite constant, so it prices the COEFFICIENT, not
# the SELECTION. But the gate was chosen by a Time x Group interaction where Group IS
# best response: markers picked using the outcome, then tested on the same outcome
# and patients. This replays the WHOLE procedure under permuted labels: shuffle Group
# at the PATIENT level (both timepoints move together), refit every marker's LMM,
# re-select the gate, rebuild the composite, refit both glms, record the LR statistic.
#
# Two selection rules on the SAME permuted fit:
#   top-K - the K markers with the smallest interaction p, K = observed gate size.
#           Never empty. THE reported rule (the one the Methods describe).
#   FDR   - BH FDR < 0.05 (the pipeline's rule). Degenerate: mostly empty gates, an
#           empty gate records LR = 0. Kept for provenance only; not the reported p.
#
# All B permutation label vectors are drawn in the parent BEFORE the parallel loop,
# so the result does not depend on the core count. Correctness gates (FATAL): the
# observed FDR gate and the top-K gate must equal the step-06 gate, and the observed
# LR must equal the step-06 LRT statistic.
# ==============================================================================

#' @title Selection-aware permutation LRT of the immune increment
#' @param DATA_T0 Step-01 standard object; @param DATA_LONG Step-01 longitudinal object.
#' @param av The step-06 added-value object (\code{clin_addval}).
#' @param config Pass config (clinical labels; \code{clinical$covariates} must be NULL).
#' @param clinical_vars Named clinical-model columns (after the NLR toggle).
#' @param B Number of permutations.
#' @param topk "gate_size" (default) or an integer K for the fixed-size rule.
#' @param cores Worker processes for \code{parallel::mclapply} (results do not depend on it).
#' @param seed Seed for the permutation draws.
#' @param fdr_thr FDR threshold of the pipeline's gate rule.
#' @param log_file Optional progress log path.
#' @return List: observed statistic and gates, both nulls (full vectors), p-values,
#'   Monte-Carlo CI, null 95th percentile and provenance.
run_selection_aware_lrt <- function(DATA_T0, DATA_LONG, av, config, clinical_vars,
                                    B = 10000L, topk = "gate_size", cores = 1L,
                                    seed = 2026L, fdr_thr = 0.05, log_file = NULL) {
  for (pkg in c("lme4", "lmerTest", "parallel"))
    if (!requireNamespace(pkg, quietly = TRUE)) stop(sprintf("[FATAL] selection-aware LRT needs package '%s'", pkg))
  log_msg <- function(...) {
    m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(...))
    message("   [SelAware] ", m)
    if (!is.null(log_file)) cat(m, "\n", file = log_file, append = TRUE)
  }
  # Leave the caller's RNG stream untouched.
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    .old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", .old_seed, envir = .GlobalEnv), add = TRUE)
  }

  R <- av_context(DATA_T0, DATA_LONG, config$clinical$responder_label,
                  config$input_file_t0, clinical_vars)
  GATE_PUB <- as.character(av$composite_markers)
  TOPK <- if (is.null(topk) || identical(as.character(topk), "gate_size")) length(GATE_PUB) else as.integer(topk)
  B <- as.integer(B); cores <- max(1L, as.integer(cores))

  # ── Longitudinal frame, exactly as src/04_longitudinal_lmm.R ──────────────────
  DL      <- DATA_LONG
  markers <- DL$hybrid_markers
  df_model <- cbind(DL$metadata, as.data.frame(DL$hybrid_data_raw[, markers]))
  df_model$Timepoint  <- factor(df_model$Timepoint, levels = c("T0", "T1"))
  df_model$Patient_ID <- as.factor(df_model$Patient_ID)
  ref_level <- config$clinical$non_responder_label
  avail     <- unique(as.character(df_model$Group))
  df_model$Group <- factor(df_model$Group, levels = c(ref_level, setdiff(avail, ref_level)))
  # The step-04 primary sweep is unadjusted only when clinical$covariates is NULL
  # (sensitivity_covariates is a separate block); the null must refit the same model.
  if (!is.null(config$clinical$covariates))
    stop("[FATAL] selection-aware LRT assumes the unadjusted step-04 sweep (clinical$covariates must be NULL)")

  pat_ids  <- levels(df_model$Patient_ID)
  pat_grp  <- as.character(df_model$Group[match(pat_ids, as.character(df_model$Patient_ID))])
  log_msg("%d markers | %d patients (%d obs) | B=%d | top-K K=%d | cores=%d",
          length(markers), length(pat_ids), nrow(df_model), B, TOPK, cores)

  # ── Fixed pieces: the delta z-matrix and the clinical design never change ─────
  A_pub <- av_assemble(R, GATE_PUB, "delta")
  Dd    <- build_timepoint_composite_data(R$DATA, R$DATA_LONG, "delta")
  Zd    <- as.data.frame(Dd$hybrid_data_z)[, markers, drop = FALSE]
  Zd[]  <- lapply(Zd, function(x) { x <- as.numeric(x); x[is.na(x)] <- median(x, na.rm = TRUE); x })
  pid_d <- as.character(Dd$metadata$Patient_ID)
  # Complete-case option: Zd is imputed on the full paired frame, then rows leave (as A_pub).
  cc_d  <- av_cc_keep(pid_d, R$DATA_LONG)
  if (!is.null(cc_d)) { Zd <- Zd[cc_d, , drop = FALSE]; pid_d <- pid_d[cc_d] }

  raw <- as.data.frame(readxl::read_excel(R$input_file, sheet = 1))
  names(raw) <- trimws(names(raw)); raw$Patient_ID <- as.character(raw$Patient_ID)
  Cd <- sapply(R$clinical_vars, function(col) suppressWarnings(as.numeric(raw[[col]][match(pid_d, raw$Patient_ID)])))
  Cd <- matrix(Cd, nrow = length(pid_d), dimnames = list(NULL, names(R$clinical_vars)))
  if (!identical(nrow(Cd), length(A_pub$yb)) || !identical(pid_d, A_pub$pid))
    stop("[FATAL] selection-aware LRT: delta frame and added-value frame disagree")

  # ── One lean interaction p per marker (same model as fit_feature_lmm's model 1) ─
  lmm_pvals <- function(dfm) {
    vapply(markers, function(mk) {
      v <- as.numeric(dfm[[mk]]); s <- sd(v, na.rm = TRUE)
      if (!is.finite(s) || s == 0) s <- 1
      d <- data.frame(V = v / s, Ti = dfm$Timepoint, G = dfm$Group, ID = dfm$Patient_ID)
      d <- d[complete.cases(d), ]
      if (nrow(d) < 10 || length(unique(d$G)) < 2) return(NA_real_)
      tryCatch({
        m <- suppressMessages(suppressWarnings(lmerTest::lmer(
          V ~ Ti * G + (1 | ID), data = d, REML = TRUE,
          control = lme4::lmerControl(calc.derivs = FALSE))))
        ct <- summary(m)$coefficients; i <- grep("Ti.*:G", rownames(ct))
        if (length(i) != 1) NA_real_ else ct[i, grep("Pr\\(>\\|t\\|\\)", colnames(ct))]
      }, error = function(e) NA_real_)
    }, numeric(1))
  }

  # LR statistic of the immune increment for a given gate and label vector.
  lr_for <- function(gate, yb) {
    if (length(gate) == 0) return(0)                 # empty gate: nothing to add
    comp <- rowMeans(Zd[, gate, drop = FALSE])
    Ci   <- av_impute_median(Cd, seq_len(nrow(Cd)))
    cn   <- make.names(colnames(Cd), unique = TRUE)
    dat  <- data.frame(yb = yb, Ci, comp = comp, check.names = TRUE)
    colnames(dat)[2:(1 + length(cn))] <- cn
    fc <- suppressWarnings(glm(reformulate(cn, "yb"), data = dat, family = binomial()))
    fk <- suppressWarnings(glm(reformulate(c(cn, "comp"), "yb"), data = dat, family = binomial()))
    d  <- deviance(fc) - deviance(fk)
    if (!is.finite(d) || d < 0) 0 else d
  }

  # ── OBSERVED ─────────────────────────────────────────────────────────────────
  p_obs    <- lmm_pvals(df_model)
  fdr_obs  <- p.adjust(p_obs, method = "BH")
  gate_obs <- markers[!is.na(fdr_obs) & fdr_obs < fdr_thr]
  topk_obs <- markers[order(p_obs, na.last = NA)][seq_len(min(TOPK, sum(!is.na(p_obs))))]
  lr_obs   <- lr_for(gate_obs, A_pub$yb)
  log_msg("observed gate (FDR<%.2f): %s | top-%d: %s", fdr_thr,
          paste(gate_obs, collapse = " + "), TOPK, paste(topk_obs, collapse = " + "))
  if (!setequal(gate_obs, GATE_PUB))
    stop(sprintf("[FATAL] selection-aware LRT: observed gate {%s} != step-06 gate {%s}",
                 paste(gate_obs, collapse = ","), paste(GATE_PUB, collapse = ",")))
  if (TOPK == length(GATE_PUB) && !setequal(topk_obs, GATE_PUB))
    stop(sprintf("[FATAL] selection-aware LRT: observed top-%d {%s} != step-06 gate {%s}", TOPK,
                 paste(topk_obs, collapse = ","), paste(GATE_PUB, collapse = ",")))
  lr_pub <- av_lrt(A_pub)$stat
  if (abs(lr_obs - lr_pub) > 1e-6)
    stop(sprintf("[FATAL] selection-aware LRT: observed LR %.6f != step-06 LRT statistic %.6f", lr_obs, lr_pub))
  log_msg("correctness gates passed (gate and LR = %.4f reproduce step 06)", lr_obs)

  # ── NULL ─────────────────────────────────────────────────────────────────────
  # All permutations drawn once in the parent (Mersenne-Twister, as in the pipeline)
  # so the result does not depend on how mclapply seeds or schedules its workers.
  .old_kind <- RNGkind()[1]
  on.exit(suppressWarnings(RNGkind(.old_kind)), add = TRUE)
  suppressWarnings(RNGkind("Mersenne-Twister"))
  set.seed(seed)
  perms <- lapply(seq_len(B), function(b) sample(pat_grp))
  map_d <- match(pid_d, pat_ids)            # paired patients -> position in pat_ids
  resp  <- R$resp_label

  t0 <- Sys.time()
  out <- parallel::mclapply(seq_len(B), function(b) {
    g_perm <- perms[[b]]
    dfp <- df_model
    dfp$Group <- factor(g_perm[match(as.character(dfp$Patient_ID), pat_ids)],
                        levels = levels(df_model$Group))
    yb_perm <- as.integer(g_perm[map_d] == resp)
    if (length(unique(yb_perm)) < 2) return(c(NA_real_, NA_real_, NA_real_))
    p   <- lmm_pvals(dfp)
    fdr <- p.adjust(p, method = "BH")
    g_fdr  <- markers[!is.na(fdr) & fdr < fdr_thr]
    g_topk <- markers[order(p, na.last = NA)][seq_len(min(TOPK, sum(!is.na(p))))]
    c(lr_for(g_fdr, yb_perm), lr_for(g_topk, yb_perm), length(g_fdr))
  }, mc.cores = cores)

  bad <- vapply(out, function(x) !is.numeric(x) || length(x) != 3, logical(1))
  if (any(bad))
    stop(sprintf("[FATAL] selection-aware LRT: %d permutation worker(s) failed", sum(bad)))
  M <- do.call(rbind, out)
  keep <- is.finite(M[, 1]) & is.finite(M[, 2])
  M <- M[keep, , drop = FALSE]
  log_msg("null done in %.1f min | %d/%d valid permutations",
          as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(M), B)

  pp <- function(nul) (1 + sum(nul >= lr_obs)) / (length(nul) + 1)
  p_topk <- pp(M[, 2]); p_fdr <- pp(M[, 1])
  mc_half <- 1.96 * sqrt(p_topk * (1 - p_topk) / nrow(M))
  naive_p <- pchisq(lr_obs, df = 1, lower.tail = FALSE)
  list(
    method = paste("selection-aware permutation LRT: patient-level label shuffle, every",
                   "marker's Time x Group LMM refitted, gate re-selected, composite rebuilt,",
                   "nested logistic LRT refitted"),
    rule = sprintf("top-%d by Time x Group interaction p (K = observed gate size)", TOPK),
    topk = TOPK, B = B, B_valid = nrow(M), seed = as.integer(seed),
    lr_obs = lr_obs,
    p_selection_aware = p_topk,
    n_null_ge_obs = as.integer(sum(M[, 2] >= lr_obs)),
    p_mc_ci = c(max(0, p_topk - mc_half), min(1, p_topk + mc_half)),
    null_q95 = unname(stats::quantile(M[, 2], 0.95)),
    null_max = max(M[, 2]),
    naive_lrt_p = naive_p,
    inflation_vs_naive = p_topk / naive_p,
    gate_fixed_perm_p = as.numeric(av$increment$lrt_perm_p),
    fdr_rule = list(p = p_fdr, frac_empty_gate = mean(M[, 3] == 0),
                    note = "FDR<0.05 re-selection: degenerate (mostly empty gates, LR = 0); provenance only, not the reported p"),
    gate_obs = gate_obs, topk_obs = topk_obs,
    n_markers = length(markers), n_patients_lmm = length(pat_ids),
    n_obs_lmm = nrow(df_model), n_paired = length(pid_d),
    lr_null_topk = M[, 2], lr_null_fdr = M[, 1], k_null_fdr = M[, 3],
    note = paste("Report p_selection_aware (top-K rule). Never quote naive_lrt_p.",
                 "gate_fixed_perm_p prices the coefficient only."))
}

#' @title JSON-ready summary of a selection-aware LRT result (drops the null vectors)
selection_aware_json <- function(sa) {
  if (is.null(sa)) return(NULL)
  sa[setdiff(names(sa), c("lr_null_topk", "lr_null_fdr", "k_null_fdr"))]
}

#' @title One-row summary table of a selection-aware LRT result (csv)
selection_aware_table <- function(sa) {
  if (is.null(sa)) return(NULL)
  data.frame(rule = sa$rule, K = sa$topk, B = sa$B, B_valid = sa$B_valid, seed = sa$seed,
             LR_observed = sa$lr_obs, n_null_ge_observed = sa$n_null_ge_obs,
             p_selection_aware = sa$p_selection_aware,
             p_mc_ci_low = sa$p_mc_ci[1], p_mc_ci_high = sa$p_mc_ci[2],
             null_q95 = sa$null_q95, null_max = sa$null_max,
             gate_fixed_perm_p = sa$gate_fixed_perm_p,
             gate = paste(sa$gate_obs, collapse = "+"),
             n_patients_lmm = sa$n_patients_lmm, n_obs_lmm = sa$n_obs_lmm, n_paired = sa$n_paired,
             stringsAsFactors = FALSE)
}
