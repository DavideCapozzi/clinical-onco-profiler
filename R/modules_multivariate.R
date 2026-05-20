# R/modules_multivariate.R
# ==============================================================================
# MULTIVARIATE ANALYSIS MODULE (sPLS-DA)
# Description: Wrappers for Sparse PLS-DA with automatic feature tuning.
# Dependencies: mixOmics, dplyr
# ==============================================================================

library(dplyr)
library(mixOmics)

#' @title Run Sparse PLS-DA (sPLS-DA) with Tuning
#' @description 
#' Fits a supervised sparse PLS-DA model.
#' Performs automatic tuning to select the optimal number of variables (keepX)
#' for each component using cross-validation.
#' 
#' @param data_z Matrix of Z-scored data (Samples x Features).
#' @param metadata Dataframe containing metadata (must match data_z rows).
#' @param group_col Name of the column in metadata to use as the class/Y.
#' @param n_comp Number of components to compute (default: 2).
#' @param validation_method CV method ("Mfold" or "loo"). 
#' @param folds Number of folds for Mfold CV (default: 5).
#' @param n_repeat Number of repetitions for Cross-Validation (default: 50).
#' @return A list containing the final optimized model, performance metrics, tuning results, and explicit validation method.
#' @export
run_splsda_model <- function(data_z, metadata, group_col = "Group", n_comp = 2, 
                             validation_method = "Mfold", folds = 5, n_repeat = 50) {
  
  requireNamespace("mixOmics", quietly = TRUE)
  
  Y <- as.factor(metadata[[group_col]])
  X <- as.matrix(data_z)
  
  if (nrow(X) != length(Y)) stop("Dimension mismatch between Data matrix and Group labels.")
  
  min_samples_per_class <- min(table(Y))
  if (validation_method == "Mfold" && min_samples_per_class < folds) {
    message(sprintf("   [sPLS-DA] Warning: Smallest class has %d samples. Adjusting CV folds to %d.", 
                    min_samples_per_class, min_samples_per_class))
    folds <- min_samples_per_class
  }
  
  if (validation_method == "Mfold" && folds < 3) {
    stop("Minimum 3 samples per class required for M-fold cross-validation. Analysis aborted.")
  }
  
  message(sprintf("   [sPLS-DA] Tuning optimal features (keepX) for %d components...", n_comp))
  
  n_vars <- ncol(X)
  if (n_vars < 5) {
    list_keepX <- c(n_vars)
  } else {
    list_keepX <- seq(5, min(30, n_vars), by = 5)
    max_val <- min(30, n_vars)
    if (!max_val %in% list_keepX) list_keepX <- c(list_keepX, max_val)
    if (n_vars <= 30 && !n_vars %in% list_keepX) list_keepX <- c(list_keepX, n_vars)
  }
  list_keepX <- sort(unique(list_keepX))
  
  tune_splsda <- mixOmics::tune.splsda(
    X = X, Y = Y, ncomp = n_comp, test.keepX = list_keepX, 
    validation = validation_method, folds = folds, 
    dist = "max.dist", progressBar = FALSE, nrepeat = n_repeat  
  )
  
  choice_keepX <- tune_splsda$choice.keepX
  message(sprintf("      -> Optimal keepX selected: Comp1=%d, Comp2=%d", 
                  choice_keepX[1], choice_keepX[2]))
  
  message("   [sPLS-DA] Fitting final sparse model with optimized parameters...")
  final_model <- mixOmics::splsda(X, Y, ncomp = n_comp, keepX = choice_keepX)
  
  perf_splsda <- tryCatch({
    mixOmics::perf(final_model, validation = validation_method, folds = folds, 
                   progressBar = FALSE, nrepeat = n_repeat, auc = TRUE) 
  }, error = function(e) {
    message(paste("      [WARN] Performance evaluation encountered an error:", e$message))
    return(NULL) 
  })
  
  # Ensure explicit propagation of the validation method to avoid "Unknown" downstream bugs
  return(list(
    model = final_model, 
    performance = perf_splsda, 
    tuning = tune_splsda,
    explicit_validation_method = validation_method
  ))
}

#' @title Extract PLS-DA Performance Metrics
#' @description Returns BER (Balanced Error Rate) from the performance object and safely extracts the validation strategy.
#' @param pls_result The list returned by run_splsda_model.
#' @return A dataframe containing Component, Overall_BER, and Validation_Method.
#' @export
extract_plsda_performance <- function(pls_result) {
  perf <- pls_result$performance
  n_comp <- pls_result$model$ncomp
  
  # Structural fallback if mixOmics performance fails
  if (is.null(perf)) {
    return(data.frame(
      Component = 1:n_comp,
      Overall_BER = NA,
      Validation_Method = "EVALUATION_FAILED"
    ))
  }
  
  ber_vals <- rep(NA, n_comp)
  tryCatch({
    if (!is.null(perf$error.rate$BER)) {
      if (is.matrix(perf$error.rate$BER)) {
        raw_ber <- perf$error.rate$BER[, "max.dist"]
        if(length(raw_ber) >= n_comp) ber_vals <- raw_ber[1:n_comp]
      } else if (is.vector(perf$error.rate$BER)) {
        ber_vals[1] <- perf$error.rate$BER["max.dist"]
      }
    }
  }, error = function(e) {
    warning("      [WARN] BER extraction encountered a matrix structure issue. Outputting NAs.")
  })
  
  # Safely extract explicitly propagated validation method
  val_method <- if(!is.null(pls_result$explicit_validation_method)) pls_result$explicit_validation_method else "Unknown"
  
  df_perf <- data.frame(
    Component = 1:n_comp,
    Overall_BER = ber_vals,
    Validation_Method = val_method 
  )
  
  return(df_perf)
}

#' @title Permutation Test on sPLS-DA Balanced Error Rate
#' @description
#' Empirical null for the per-component BER under a shuffled outcome. Refits
#' the sPLS-DA at the same keepX as the observed model (does NOT re-tune —
#' keepX tuning under shuffled labels is dominated by noise and doubles the
#' compute budget without changing the inference) and measures the per-fold
#' BER on each shuffle. The observed BER (computed at full CV repeats) is
#' compared against the permuted distribution; we report a one-sided p-value
#' (lower BER = better than null).
#'
#' @param data_z Z-scored data matrix (samples x features), same as the one
#'   passed to run_splsda_model().
#' @param metadata Metadata data.frame aligned with data_z rows.
#' @param keepX Integer vector of per-component sparsity (typically
#'   pls_result$model$keepX), used in every permutation refit.
#' @param observed_ber Numeric vector of observed BER per component (from
#'   extract_plsda_performance(pls_result)$Overall_BER). Length == length(keepX).
#' @param group_col Metadata column with the outcome (default "Group").
#' @param validation_method "Mfold" or "loo".
#' @param folds Folds per CV pass under permutation (default 5).
#' @param n_repeat CV repeats per permutation (default 1 — single-pass is
#'   sufficient because we average across permutations).
#' @param n_perm Permutations (default 200).
#' @param seed RNG seed (default 2026).
#' @return Data.frame with one row per component: Component, Observed_BER,
#'   Perm_Mean, Perm_SD, Perm_P_Value, N_Valid_Perm. Plus attribute "n_perm".
run_splsda_permutation <- function(data_z, metadata, keepX, observed_ber,
                                   group_col       = "Group",
                                   validation_method = "Mfold",
                                   folds           = 5,
                                   n_repeat        = 1,
                                   n_perm          = 200,
                                   seed            = 2026) {
  requireNamespace("mixOmics", quietly = TRUE)

  Y_obs <- as.factor(metadata[[group_col]])
  X     <- as.matrix(data_z)
  n_comp <- length(keepX)
  if (length(observed_ber) != n_comp) {
    stop("[Permutation] observed_ber length must match length(keepX).")
  }

  # Folds may already have been auto-shrunk in run_splsda_model when the
  # smallest class is small; apply the same clamp here for consistency.
  min_class_n <- min(table(Y_obs))
  if (validation_method == "Mfold" && min_class_n < folds) folds <- min_class_n
  if (validation_method == "Mfold" && folds < 3) {
    warning("[Permutation] Smallest class < 3 samples; permutation test skipped.")
    return(NULL)
  }

  perm_ber <- matrix(NA_real_, nrow = n_perm, ncol = n_comp)
  set.seed(seed)
  for (b in seq_len(n_perm)) {
    Y_shuf <- sample(Y_obs)
    fit <- tryCatch(
      mixOmics::splsda(X, Y_shuf, ncomp = n_comp, keepX = keepX),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    pr <- tryCatch(
      mixOmics::perf(fit, validation = validation_method, folds = folds,
                     progressBar = FALSE, nrepeat = n_repeat, auc = FALSE),
      error = function(e) NULL
    )
    if (is.null(pr) || is.null(pr$error.rate$BER)) next
    ber_row <- pr$error.rate$BER
    if (is.matrix(ber_row)) {
      v <- as.numeric(ber_row[, "max.dist"])
    } else {
      v <- rep(NA_real_, n_comp)
      v[1] <- as.numeric(ber_row["max.dist"])
    }
    if (length(v) >= n_comp) perm_ber[b, ] <- v[seq_len(n_comp)]
  }

  res <- data.frame(
    Component   = seq_len(n_comp),
    Observed_BER = round(observed_ber, 4),
    stringsAsFactors = FALSE
  )
  res$Perm_Mean <- round(colMeans(perm_ber, na.rm = TRUE), 4)
  res$Perm_SD   <- round(apply(perm_ber, 2, sd, na.rm = TRUE), 4)
  res$N_Valid_Perm <- as.integer(colSums(!is.na(perm_ber)))

  # One-sided permutation p-value: fraction of permuted BERs <= observed,
  # with the standard (k+1)/(N+1) small-sample adjustment (Phipson & Smyth).
  res$Perm_P_Value <- vapply(seq_len(n_comp), function(j) {
    v <- perm_ber[, j]; v <- v[!is.na(v)]
    if (length(v) == 0) return(NA_real_)
    round((sum(v <= observed_ber[j]) + 1) / (length(v) + 1), 4)
  }, numeric(1))

  attr(res, "n_perm") <- as.integer(n_perm)
  res
}


#' @title Extract Discriminant Features (Sparse Loadings)
#' @description 
#' Extracts loading weights from the sPLS-DA model. 
#' Filters out variables with zero weight (not selected by the model).
#' 
#' @param pls_result The list returned by run_splsda_model.
#' @return A dataframe of non-zero loadings.
extract_plsda_loadings <- function(pls_result) {
  model <- pls_result$model
  
  # Check if loadings exist
  if(is.null(model$loadings$X)) return(data.frame())
  
  # Convert loadings matrix to dataframe
  loadings_raw <- as.data.frame(model$loadings$X)
  loadings_raw$Marker <- rownames(loadings_raw)
  
  # Reshape to long format or keep wide? 
  # Strategy: Create a clean table with Comp1 and Comp2 weights
  
  df_clean <- loadings_raw %>%
    dplyr::select(Marker, everything()) %>%
    dplyr::rename(Comp1_Weight = comp1)
  
  if(model$ncomp > 1) {
    df_clean <- df_clean %>% dplyr::rename(Comp2_Weight = comp2)
  }
  
  # Calculate Importance (Absolute weight on Comp1)
  df_clean$Importance <- abs(df_clean$Comp1_Weight)
  
  # Filter: In sPLS-DA, we only care about non-zero weights
  # We keep rows where at least one component has a non-zero weight
  df_clean <- df_clean %>%
    filter(Importance > 0 | (if("Comp2_Weight" %in% names(.)) abs(Comp2_Weight) > 0 else FALSE)) %>%
    arrange(desc(Importance)) %>%
    mutate(
      Direction = ifelse(Comp1_Weight > 0, "Positive_Assoc", "Negative_Assoc")
    )
  
  return(df_clean)
}