# R/modules_ml_bootstrap_ci.R
# ==============================================================================
# PATIENT-BOOTSTRAP CIs for the per-arm CV AUCs and the ΔAUC increment
# (formerly diagnostics/diag_45).
#
# cv_kfold$delta_auc_ci in the step-06 JSON is quantile() over the 50 REPETITION-level
# deltas of repeated CV on one fixed set of patients: Monte-Carlo spread of the fold
# assignment, whose width -> 0 as reps -> Inf. It contains NO patient-level sampling
# uncertainty and is NOT a confidence interval. This module gives the reportable one.
#
#   (A) full-procedure patient bootstrap - resample the patients, re-run the ENTIRE
#       repeated k-fold CV inside each resample. Prices patient sampling AND fold noise
#       AND coefficient instability. PRIMARY (the reported CIs).
#       Caveat: resamples contain duplicate patients, so a duplicate can sit in both
#       train and test within a resample; the per-arm AUCs of a resample are inflated,
#       but the deltas compare two arms on identical folds, so the bias largely cancels.
#   (B) conditional patient bootstrap - resample patients over the FIXED out-of-fold
#       probabilities. Conditions on the fitted models. Narrower by design.
#   (C) DeLong paired CI on the same out-of-fold probabilities.
#
# Resample indices are drawn ONCE in the parent under a fixed seed and each resample's
# CV seeds itself (seed + b), so the result does not depend on the core count.
# ==============================================================================

#' @title Out-of-fold probabilities per arm on the original frame (for B and C)
av_oof_probs <- function(A, reps = 50L, seed = 2026L + 977L) {
  yb <- A$yb; n <- length(yb); set.seed(seed)
  k <- min(10L, max(3L, floor(min(sum(yb), sum(1 - yb)) / 2)))
  acc <- matrix(0, n, 3, dimnames = list(NULL, c("clinical", "immune", "combined")))
  for (r in seq_len(reps)) {
    fold <- integer(n)
    for (cl in unique(yb)) { ix <- which(yb == cl)
      fold[ix] <- sample(rep(seq_len(k), length.out = length(ix))) }
    for (f in seq_len(k)) {
      te <- which(fold == f); tr <- which(fold != f)
      if (length(unique(yb[tr])) < 2) next
      acc[te, "clinical"] <- acc[te, "clinical"] + av_fit_predict(A$Cmat, A$comp, yb, tr, te, FALSE)
      acc[te, "combined"] <- acc[te, "combined"] + av_fit_predict(A$Cmat, A$comp, yb, tr, te, TRUE)
      di <- data.frame(yb = yb, comp = A$comp)
      fi <- suppressWarnings(glm(yb ~ comp, data = di[tr, , drop = FALSE], family = binomial()))
      acc[te, "immune"] <- acc[te, "immune"] +
        as.numeric(predict(fi, newdata = di[te, , drop = FALSE], type = "response"))
    }
  }
  acc / reps      # mean out-of-fold probability per patient across repetitions
}

#' @title Patient-bootstrap CIs of the CV AUCs for one marker set
#' @param R Context from \code{av_context()}; @param markers Marker set (composite).
#' @param timepoint Composite timepoint (the step-06 one).
#' @param B Bootstrap resamples; @param reps_in_boot CV repetitions inside each resample.
#' @param seed Base seed (CV point estimate seed + 977, resamples seed + 101,
#'   resample b's CV seed + b).
#' @param cores Worker processes (results do not depend on it).
#' @param single Also bootstrap each clinical var alone (same folds inside each resample,
#'   no extra RNG draws, so every other number is unchanged); only when > 1 clinical var.
av_bootstrap_marker_set <- function(R, markers, timepoint, B, reps_in_boot, seed, cores,
                                    single = FALSE) {
  A <- av_assemble(R, markers, timepoint)
  cv_full <- av_cv_kfold(A, seed = seed + 977L, single = single)   # step-06-style point estimates
  sv <- names(cv_full$auc_single)
  nA <- 5L + length(sv)
  P  <- av_oof_probs(A, seed = seed + 977L)
  ap <- function(p) av_auc_pos(p, A$y, A$neg_label, A$resp_label)

  set.seed(seed + 101L)
  idx_list <- lapply(seq_len(B), function(b) sample(length(A$yb), replace = TRUE))

  # (A) full-procedure patient bootstrap -------------------------------------
  bootA <- parallel::mclapply(seq_len(B), function(b) {
    idx <- idx_list[[b]]
    if (length(unique(A$yb[idx])) < 2) return(rep(NA_real_, nA))
    Ab <- list(comp = A$comp[idx], Cmat = A$Cmat[idx, , drop = FALSE],
               yb = A$yb[idx], y = A$y[idx],
               neg_label = A$neg_label, resp_label = A$resp_label)
    cvb <- tryCatch(av_cv_kfold(Ab, reps = reps_in_boot, seed = seed + b, single = length(sv) > 0),
                    error = function(e) NULL)
    if (is.null(cvb)) return(rep(NA_real_, nA))
    c(cvb$delta_auc, cvb$delta_auc_vs_immune, cvb$auc_clinical, cvb$auc_immune, cvb$auc_combined,
      vapply(sv, function(v) { a <- cvb$auc_single[[v]]; if (is.null(a)) NA_real_ else a }, numeric(1)))
  }, mc.cores = cores)
  if (any(vapply(bootA, function(x) !is.numeric(x) || length(x) != nA, logical(1))))
    stop("[FATAL] patient bootstrap: a worker failed (full-procedure resamples)")
  MA <- do.call(rbind, bootA)
  MA <- MA[complete.cases(MA[, 1:5, drop = FALSE]), , drop = FALSE]   # validity = the 5 original arms

  # (B) conditional patient bootstrap over fixed out-of-fold probs ------------
  bootB <- parallel::mclapply(seq_len(B), function(b) {
    idx <- idx_list[[b]]          # same resamples as (A), so the two are comparable
    yv <- A$y[idx]
    if (length(unique(yv)) < 2) return(c(NA_real_, NA_real_))
    a <- function(col) av_auc_pos(P[idx, col], yv, A$neg_label, A$resp_label)
    cb <- a("combined")
    c(cb - a("clinical"), cb - a("immune"))
  }, mc.cores = cores)
  if (any(vapply(bootB, function(x) !is.numeric(x) || length(x) != 2, logical(1))))
    stop("[FATAL] patient bootstrap: a worker failed (conditional resamples)")
  MB <- do.call(rbind, bootB); MB <- MB[complete.cases(MB), , drop = FALSE]

  # (C) DeLong paired CI ------------------------------------------------------
  rc <- pROC::roc(A$y, P[, "clinical"], levels = c(A$neg_label, A$resp_label), direction = "<", quiet = TRUE)
  rk <- pROC::roc(A$y, P[, "combined"], levels = c(A$neg_label, A$resp_label), direction = "<", quiet = TRUE)
  ri <- pROC::roc(A$y, P[, "immune"],   levels = c(A$neg_label, A$resp_label), direction = "<", quiet = TRUE)
  dl_ci <- function(r1, r2) {
    t <- pROC::roc.test(r1, r2, method = "delong", paired = TRUE)
    d <- as.numeric(t$estimate[1] - t$estimate[2])
    se <- abs(d) / abs(stats::qnorm(t$p.value / 2)); if (!is.finite(se)) se <- NA_real_
    list(delta = d, ci = c(d - 1.96 * se, d + 1.96 * se), p = t$p.value)
  }
  dC_clin <- dl_ci(rk, rc); dC_imm <- dl_ci(rk, ri)

  qci <- function(v) unname(quantile(v, c(.025, .975), na.rm = TRUE))
  lr  <- av_lrt(A)

  out <- list(
    markers = markers, n = length(A$yb), n_pos = sum(A$yb),
    point = list(auc_clinical = cv_full$auc_clinical, auc_immune = cv_full$auc_immune,
                 auc_combined = cv_full$auc_combined,
                 delta_auc = cv_full$delta_auc, delta_auc_vs_immune = cv_full$delta_auc_vs_immune),
    mc_spread = list(delta_auc = qci(cv_full$R[, "combined"] - cv_full$R[, "clinical"]),
                     delta_auc_vs_immune = qci(cv_full$R[, "combined"] - cv_full$R[, "immune"])),
    bootA = list(n_valid = nrow(MA),
                 delta_auc = qci(MA[, 1]), delta_auc_vs_immune = qci(MA[, 2]),
                 median_delta = median(MA[, 1], na.rm = TRUE),
                 auc = list(clinical = qci(MA[, 3]), immune = qci(MA[, 4]), combined = qci(MA[, 5])),
                 p_two_sided = 2 * min(mean(MA[, 1] <= 0), mean(MA[, 1] >= 0))),
    bootB = list(n_valid = nrow(MB),
                 delta_auc = qci(MB[, 1]), delta_auc_vs_immune = qci(MB[, 2])),
    delong = list(delta_auc = dC_clin$delta, delta_auc_ci = dC_clin$ci, delta_auc_p = dC_clin$p,
                  delta_vs_immune = dC_imm$delta, delta_vs_immune_ci = dC_imm$ci,
                  delta_vs_immune_p = dC_imm$p),
    lrt = list(stat = lr$stat, p = lr$p),
    auc_oof = c(clinical = ap(P[, "clinical"]), immune = ap(P[, "immune"]), combined = ap(P[, "combined"]))
  )
  if (length(sv)) {                     # each clinical var alone (keys absent otherwise)
    out$point$auc_single <- cv_full$auc_single
    out$bootA$auc_single <- setNames(lapply(seq_along(sv), function(j) qci(MA[, 5L + j])), sv)
  }
  out
}

#' @title Patient-bootstrap AUC CIs for the gate composite and each gate marker alone
#' @param DATA_T0,DATA_LONG Step-01 objects; @param av The step-06 \code{clin_addval}.
#' @param config Pass config; @param clinical_vars Clinical-model columns (after the NLR toggle).
#' @param B Resamples (2000 for the reported CIs); @param reps_in_boot CV reps per resample.
#' @param seed Base seed; @param cores Worker processes.
#' @return List: \code{results} (per marker set: "gate" then one entry per marker),
#'   B, reps, seed, the step-06 Monte-Carlo spread for contrast, and a note.
run_patient_bootstrap_auc <- function(DATA_T0, DATA_LONG, av, config, clinical_vars,
                                      B = 2000L, reps_in_boot = 10L, seed = 2026L, cores = 1L) {
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    .old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", .old_seed, envir = .GlobalEnv), add = TRUE)
  }
  .old_kind <- RNGkind()[1]
  on.exit(suppressWarnings(RNGkind(.old_kind)), add = TRUE)
  suppressWarnings(RNGkind("Mersenne-Twister"))

  R <- av_context(DATA_T0, DATA_LONG, config$clinical$responder_label,
                  config$input_file_t0, clinical_vars)
  gate <- as.character(av$composite_markers)
  sets <- c(list(gate = gate), setNames(as.list(gate), gate))
  if (length(gate) == 1) sets <- sets[1]
  t0 <- Sys.time()
  res <- lapply(names(sets), function(nm) {
    message(sprintf("   [PatientBoot] %s (%s): B=%d", nm, paste(sets[[nm]], collapse = "+"), B))
    av_bootstrap_marker_set(R, sets[[nm]], av$composite_timepoint, as.integer(B),
                            as.integer(reps_in_boot), as.integer(seed), max(1L, as.integer(cores)),
                            single = identical(nm, "gate"))   # clinical vars alone: gate set only
  })
  names(res) <- names(sets)
  message(sprintf("   [PatientBoot] done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  list(results = res, B = as.integer(B), reps_in_bootstrap = as.integer(reps_in_boot),
       seed = as.integer(seed), composite_timepoint = av$composite_timepoint,
       step06_mc_spread_delta_auc = as.numeric(av$cv_kfold$delta_auc_ci),
       note = paste("cv_kfold$delta_auc_ci in the JSON is CV Monte-Carlo spread across 50",
                    "repetitions on fixed patients, NOT a confidence interval. bootA (the",
                    "full-procedure patient bootstrap) is the reportable interval."))
}

#' @title Long summary table of a patient-bootstrap result (csv)
patient_bootstrap_table <- function(pb) {
  if (is.null(pb)) return(NULL)
  full <- "patient bootstrap (full procedure)"
  out <- do.call(rbind, lapply(names(pb$results), function(nm) {
    r <- pb$results[[nm]]; a <- r$bootA
    d <- data.frame(
      marker_set = nm, markers = paste(r$markers, collapse = "+"), n = r$n, n_pos = r$n_pos,
      quantity = c("AUC clinical", "AUC immune", "AUC combined", "dAUC combined - clinical",
                   "dAUC combined - immune", "dAUC combined - clinical", "dAUC combined - clinical"),
      point    = c(r$point$auc_clinical, r$point$auc_immune, r$point$auc_combined, r$point$delta_auc,
                   r$point$delta_auc_vs_immune, r$point$delta_auc, r$delong$delta_auc),
      ci_low   = c(a$auc$clinical[1], a$auc$immune[1], a$auc$combined[1], a$delta_auc[1],
                   a$delta_auc_vs_immune[1], r$mc_spread$delta_auc[1], r$delong$delta_auc_ci[1]),
      ci_high  = c(a$auc$clinical[2], a$auc$immune[2], a$auc$combined[2], a$delta_auc[2],
                   a$delta_auc_vs_immune[2], r$mc_spread$delta_auc[2], r$delong$delta_auc_ci[2]),
      interval = c(rep(full, 5), "CV Monte-Carlo spread (NOT a CI)",
                   "DeLong paired (out-of-fold probabilities)"),
      boot_p_two_sided = c(NA, NA, NA, a$p_two_sided, NA, NA, NA),
      n_valid_resamples = c(rep(a$n_valid, 5), NA, NA),
      stringsAsFactors = FALSE)
    sv <- names(r$point$auc_single)
    if (length(sv)) d <- rbind(d, data.frame(   # each clinical var alone (same CV, same resamples)
      marker_set = nm, markers = paste(r$markers, collapse = "+"), n = r$n, n_pos = r$n_pos,
      quantity = sprintf("AUC %s alone", sv), point = unlist(r$point$auc_single[sv]),
      ci_low = vapply(sv, function(v) a$auc_single[[v]][1], numeric(1)),
      ci_high = vapply(sv, function(v) a$auc_single[[v]][2], numeric(1)),
      interval = full, boot_p_two_sided = NA, n_valid_resamples = a$n_valid,
      stringsAsFactors = FALSE, row.names = NULL))
    d
  }))
  out$B <- pb$B; out$seed <- pb$seed
  out
}
