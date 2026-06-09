# R/modules_ml_cv.R
# ==============================================================================
# MACHINE LEARNING MODULE — CV
# Split from the former monolithic R/modules_ml.R (sourced via that aggregator).
# Dependencies: dplyr, ggplot2, tidyr, glmnet, e1071, pROC
# ==============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)

#' @title Filter Collinear Features
#' @description
#' Iteratively removes near-collinear markers from a feature matrix. At each
#' step, the pair with the highest absolute Pearson correlation is identified;
#' the marker with the higher MEAN absolute correlation to all other remaining
#' markers is dropped (the more globally redundant one). Repeats until no pair
#' exceeds the threshold. This correctly handles chains of collinearity without
#' relying on univariate AUC, which can be noisy in small samples.
#'
#' @param X_main Numeric matrix. Main-effect feature matrix (samples x features).
#' @param cor_threshold Numeric. Drop one of a pair when |r| exceeds this value (default 0.85).
#' @return A named list: X (filtered matrix), kept (character), dropped (character).
#' @export
filter_collinear_features <- function(X_main, cor_threshold = 0.85) {
  markers <- colnames(X_main)
  dropped <- character(0)

  repeat {
    if (length(markers) <= 1) break

    cm <- abs(cor(X_main[, markers, drop = FALSE], method = "pearson"))
    diag(cm) <- 0
    max_r <- max(cm, na.rm = TRUE)
    if (!is.finite(max_r) || max_r < cor_threshold) break

    # Identify the most correlated pair
    idx    <- which(cm == max_r, arr.ind = TRUE)[1, ]
    m1     <- markers[idx[1]]
    m2     <- markers[idx[2]]
    others <- markers[!markers %in% c(m1, m2)]

    if (length(others) == 0) {
      # Only two markers left and both are collinear with each other.
      # With no other markers to compute mean |r| against, break the tie
      # by convention (drop m2, keep m1).
      dropped <- c(dropped, m2)
      markers <- m1
      break
    }

    # Drop whichever has higher mean |r| to all other markers (more redundant)
    mean_r_m1 <- mean(abs(cm[m1, others]), na.rm = TRUE)
    mean_r_m2 <- mean(abs(cm[m2, others]), na.rm = TRUE)
    to_drop   <- if (mean_r_m1 >= mean_r_m2) m1 else m2

    dropped <- c(dropped, to_drop)
    markers <- markers[markers != to_drop]
  }

  list(X = X_main[, markers, drop = FALSE], kept = markers, dropped = dropped)
}

#' @title Build ML Feature Matrix
#' @description
#' Extracts the Z-scored feature columns for the LOO-robust markers from the
#' Step 01 processed data object. Applies a collinearity filter before adding
#' interaction terms to prevent near-duplicate features from degrading the
#' classifier. Optionally appends pairwise interaction terms guided by the
#' network co-activation topology.
#'
#' @param DATA Named list. Processed data object from Step 01 (standard mode).
#' @param robust_markers Character vector. Markers selected by the LMM LOO gate.
#' @param include_interactions Logical. Append all pairwise products (default TRUE).
#' @param cor_threshold Numeric. Collinearity filter threshold (default 0.85).
#'   Set to 1 to disable filtering.
#' @return A named list: X (matrix), y (factor), feature_names (character), n_main (integer).
#' @export
build_ml_matrix <- function(DATA, robust_markers, include_interactions = TRUE,
                            cor_threshold = 0.85) {
  # Resolve marker availability against the transformed Z-score matrix
  all_z_cols  <- DATA$hybrid_markers
  available   <- intersect(robust_markers, all_z_cols)

  if (length(available) < length(robust_markers)) {
    missing <- setdiff(robust_markers, all_z_cols)
    warning(sprintf("[ML] %d marker(s) not found in Step 01 Z-score matrix: %s",
                    length(missing), paste(missing, collapse = ", ")))
  }
  if (length(available) == 0) stop("[ML][FATAL] No robust markers present in data matrix.")

  X <- as.matrix(DATA$hybrid_data_z[, available, drop = FALSE])
  mode(X) <- "numeric"
  rownames(X) <- make.unique(as.character(DATA$metadata$Patient_ID))

  # Collinearity filter: drop the more redundant of any near-duplicate pair
  if (length(available) >= 2 && is.finite(cor_threshold) && cor_threshold < 1) {
    filt <- filter_collinear_features(X, cor_threshold = cor_threshold)
    if (length(filt$dropped) > 0) {
      message(sprintf(
        "   [ML] Collinearity filter (|r| > %.2f): removed %s (kept %s)",
        cor_threshold,
        paste(filt$dropped, collapse = ", "),
        paste(filt$kept,    collapse = ", ")
      ))
    }
    X         <- filt$X
    available <- filt$kept
  }

  n_main <- length(available)

  if (include_interactions && n_main >= 2) {
    combos <- utils::combn(available, 2, simplify = FALSE)
    for (pair in combos) {
      col_name        <- paste(pair[1], pair[2], sep = "_x_")
      X               <- cbind(X, as.numeric(X[, pair[1]]) * as.numeric(X[, pair[2]]))
      colnames(X)[ncol(X)] <- col_name
    }
    message(sprintf("   [ML] Feature matrix: %d main effects + %d interaction term(s) = %d total",
                    n_main, length(combos), ncol(X)))
  } else {
    message(sprintf("   [ML] Feature matrix: %d main effect(s) (interactions disabled or p < 2)", n_main))
  }

  y <- DATA$metadata$Group  # Factor: levels = c(non_responder, responder)

  list(X = X, y = y, feature_names = colnames(X), n_main = n_main)
}

#' @title Run Nested-LOOCV with Elastic Net Logistic Regression
#' @description
#' Outer loop: Leave-One-Out (n folds). Inner loop: k-fold cross-validation via
#' glmnet::cv.glmnet to jointly select alpha (Ridge/Elastic Net/LASSO) and lambda.
#' No feature re-selection occurs inside the loop; the feature set is frozen prior
#' to LOOCV, ensuring no data leakage from the classification data.
#'
#' @param X Numeric matrix. Feature matrix (samples x features).
#' @param y Factor. Binary outcome (two levels; positive class = second level).
#' @param alpha_grid Numeric vector. Alpha values to try in inner CV (default c(0, 0.5, 1)).
#' @param n_lambda Integer. Number of lambda values on the regularization path (default 100).
#' @param inner_folds Integer. k for inner k-fold CV (default 5).
#' @param seed Integer. Random seed for reproducibility.
#' @return A named list: method, predicted_probs, y_true, positive_label, coef_matrix, feature_names.
#' @export
run_nested_loocv_glmnet <- function(X, y,
                                    alpha_grid  = c(0, 0.5, 1),
                                    n_lambda    = 100L,
                                    inner_folds = 5L,
                                    seed        = 2026) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("[ML] Package 'glmnet' required. Add r-glmnet to env/environment.yml and rebuild.")
  }

  n           <- nrow(X)
  pos_label   <- levels(y)[2]
  y_bin       <- as.integer(y) - 1L  # 0 = reference, 1 = positive class

  pred_probs  <- numeric(n)
  coef_list   <- vector("list", n)
  coef_names  <- c("(Intercept)", colnames(X))

  message(sprintf("   [Elastic Net] Outer LOO: %d folds | Inner CV: %d folds | alpha grid: {%s}",
                  n, inner_folds, paste(alpha_grid, collapse = ", ")))

  for (i in seq_len(n)) {
    X_train <- X[-i, , drop = FALSE]
    y_train <- y_bin[-i]
    X_test  <- X[i,  , drop = FALSE]

    n_min_class  <- min(sum(y_train == 0L), sum(y_train == 1L))
    safe_folds   <- max(2L, min(as.integer(inner_folds), n_min_class))

    # Inner CV: select best (alpha, lambda) by maximum inner AUC
    best_alpha      <- alpha_grid[1]
    best_lambda     <- NULL
    best_inner_auc  <- -Inf

    for (a in alpha_grid) {
      set.seed(seed + i)
      cv_fit <- tryCatch(
        glmnet::cv.glmnet(
          x            = X_train,
          y            = as.numeric(y_train),
          family       = "binomial",
          alpha        = a,
          nfolds       = safe_folds,
          type.measure = "auc",
          nlambda      = as.integer(n_lambda)
        ),
        error = function(e) NULL
      )

      if (!is.null(cv_fit)) {
        inner_auc <- suppressWarnings(max(cv_fit$cvm, na.rm = TRUE))
        if (is.finite(inner_auc) && inner_auc > best_inner_auc) {
          best_inner_auc <- inner_auc
          best_alpha     <- a
          best_lambda    <- cv_fit$lambda.min
        }
      }
    }

    # Fit final model on full training fold with best hyperparameters.
    # glmnet with a fixed lambda is deterministic — no set.seed needed here.
    final_fit <- tryCatch(
      glmnet::glmnet(X_train, as.numeric(y_train), family = "binomial",
                     alpha = best_alpha, lambda = best_lambda),
      error = function(e) NULL
    )

    if (!is.null(final_fit) && !is.null(best_lambda)) {
      pred_probs[i] <- tryCatch(
        as.numeric(predict(final_fit, newx = X_test, s = best_lambda, type = "response")),
        error = function(e) 0.5
      )
      coef_list[[i]] <- tryCatch(
        as.numeric(coef(final_fit, s = best_lambda)),
        error = function(e) rep(NA_real_, length(coef_names))
      )
    } else {
      pred_probs[i]  <- 0.5
      coef_list[[i]] <- rep(NA_real_, length(coef_names))
    }
  }

  coef_matrix           <- do.call(rbind, coef_list)
  colnames(coef_matrix) <- coef_names

  list(
    method        = "Elastic Net Logistic Regression (Nested-LOOCV)",
    predicted_probs = pred_probs,
    y_true        = y,
    positive_label = pos_label,
    coef_matrix   = coef_matrix,
    feature_names = colnames(X)
  )
}

#' @title Run Nested-LOOCV with SVM (RBF Kernel)
#' @description
#' Outer loop: Leave-One-Out (n folds). Inner loop: grid search over (C, gamma)
#' via e1071::tune() with k-fold CV. Probability estimates obtained via Platt
#' scaling (probability = TRUE in e1071::svm). Serves as a non-linear comparison
#' arm to validate whether a linear decision boundary is sufficient.
#'
#' @param X Numeric matrix. Feature matrix (samples x features).
#' @param y Factor. Binary outcome (two levels; positive class = second level).
#' @param C_grid Numeric vector. Cost values for inner grid search.
#' @param gamma_grid Numeric vector. Gamma values for inner grid search.
#' @param inner_folds Integer. k for inner k-fold CV (default 5).
#' @param seed Integer. Random seed for reproducibility.
#' @return A named list: method, predicted_probs, y_true, positive_label, coef_matrix (NULL), feature_names.
#' @export
run_nested_loocv_svm <- function(X, y,
                                 C_grid      = c(0.01, 0.1, 1, 10, 100),
                                 gamma_grid  = c(0.01, 0.1, 1, 10),
                                 inner_folds = 5L,
                                 seed        = 2026) {
  if (!requireNamespace("e1071", quietly = TRUE)) {
    stop("[ML] Package 'e1071' required. Add r-e1071 to env/environment.yml and rebuild.")
  }

  n          <- nrow(X)
  pos_label  <- levels(y)[2]
  pred_probs <- numeric(n)

  message(sprintf("   [SVM-RBF] Outer LOO: %d folds | Inner CV: %d folds | Grid: C(%d) x gamma(%d)",
                  n, inner_folds, length(C_grid), length(gamma_grid)))

  for (i in seq_len(n)) {
    X_train <- X[-i, , drop = FALSE]
    y_train <- y[-i]
    X_test  <- X[i,  , drop = FALSE]

    n_min_class <- min(table(y_train))
    safe_folds  <- max(2L, min(as.integer(inner_folds), n_min_class))

    set.seed(seed + i)
    tune_res <- tryCatch(
      suppressWarnings(e1071::tune(
        e1071::svm,
        train.x    = X_train,
        train.y    = y_train,
        kernel     = "radial",
        ranges     = list(cost = C_grid, gamma = gamma_grid),
        tunecontrol = e1071::tune.control(sampling = "cross", cross = safe_folds)
      )),
      error = function(e) NULL
    )

    if (!is.null(tune_res)) {
      bp <- tune_res$best.parameters

      set.seed(seed + i)
      final_svm <- tryCatch(
        e1071::svm(x = X_train, y = y_train, kernel = "radial",
                   cost = bp$cost, gamma = bp$gamma, probability = TRUE),
        error = function(e) NULL
      )

      if (!is.null(final_svm)) {
        pred_obj <- tryCatch(
          predict(final_svm, newdata = X_test, probability = TRUE),
          error = function(e) NULL
        )

        if (!is.null(pred_obj)) {
          prob_mat       <- attr(pred_obj, "probabilities")
          pred_probs[i]  <- if (!is.null(prob_mat) && pos_label %in% colnames(prob_mat)) {
            as.numeric(prob_mat[1, pos_label])
          } else 0.5
        } else {
          pred_probs[i] <- 0.5
        }
      } else {
        pred_probs[i] <- 0.5
      }
    } else {
      pred_probs[i] <- 0.5
    }
  }

  list(
    method         = "SVM with RBF Kernel (Nested-LOOCV)",
    predicted_probs = pred_probs,
    y_true         = y,
    positive_label = pos_label,
    coef_matrix    = NULL,
    feature_names  = colnames(X)
  )
}

#' @title Compute Classification Performance Metrics
#' @description
#' Computes AUC with 95% CI (DeLong method), balanced accuracy, sensitivity,
#' and specificity at the Youden-J optimal threshold from out-of-fold predictions.
#'
#' @param y_true Factor. True class labels.
#' @param y_pred_prob Numeric vector. Predicted probabilities for the positive class.
#' @param positive_label String. The positive class label (second factor level).
#' @return A named list: auc, auc_ci (length-3 numeric: lower, estimate, upper),
#'   balanced_accuracy, ber, sensitivity, specificity, threshold.
#' @export
compute_classification_metrics <- function(y_true, y_pred_prob, positive_label) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("[ML] Package 'pROC' required. Add r-proc to env/environment.yml and rebuild.")
  }

  y_bin   <- as.integer(y_true == positive_label)
  na_mask <- !is.na(y_pred_prob)

  if (sum(na_mask) < length(y_bin)) {
    warning(sprintf("[ML] %d NA predictions replaced with 0.5 before metric computation.",
                    sum(!na_mask)))
    y_pred_prob[!na_mask] <- 0.5
  }

  roc_obj <- tryCatch(
    pROC::roc(response = y_bin, predictor = y_pred_prob, direction = "<", quiet = TRUE),
    error = function(e) {
      warning(paste("[ML] pROC::roc() failed:", e$message))
      return(NULL)
    }
  )

  if (is.null(roc_obj)) {
    return(list(auc = NA_real_, auc_ci = c(NA, NA, NA),
                balanced_accuracy = NA_real_, ber = NA_real_,
                sensitivity = NA_real_, specificity = NA_real_, threshold = 0.5))
  }

  auc_val <- as.numeric(pROC::auc(roc_obj))
  auc_ci  <- tryCatch(
    as.numeric(pROC::ci.auc(roc_obj, method = "delong")),
    error = function(e) c(NA_real_, auc_val, NA_real_)
  )

  # Youden-J optimal threshold
  coords_df <- tryCatch(
    pROC::coords(roc_obj, x = "best", best.method = "youden",
                 ret = c("threshold", "sensitivity", "specificity"),
                 transpose = FALSE),
    error = function(e) NULL
  )

  if (!is.null(coords_df) && nrow(coords_df) > 0) {
    opt_thresh <- as.numeric(coords_df$threshold[1])
    opt_sens   <- as.numeric(coords_df$sensitivity[1])
    opt_spec   <- as.numeric(coords_df$specificity[1])
  } else {
    opt_thresh <- 0.5
    opt_sens   <- NA_real_
    opt_spec   <- NA_real_
  }

  bal_acc <- if (!is.na(opt_sens) && !is.na(opt_spec)) (opt_sens + opt_spec) / 2 else NA_real_

  list(
    auc               = auc_val,
    auc_ci            = auc_ci,
    balanced_accuracy = bal_acc,
    ber               = if (!is.na(bal_acc)) 1 - bal_acc else NA_real_,
    sensitivity       = opt_sens,
    specificity       = opt_spec,
    threshold         = opt_thresh
  )
}

#' @title Permutation AUC Test
#' @description
#' Empirical p-value for the observed AUC by permuting class labels against
#' fixed out-of-fold predicted probabilities. No ML models are re-fitted.
#' Uses the inclusive formula p = (B_extreme + 1) / (n_perm + 1).
#'
#' @param y_true Factor. True class labels (same order as predicted_probs).
#' @param predicted_probs Numeric vector. Out-of-fold predicted probabilities.
#' @param positive_label String. The positive class label.
#' @param n_perm Integer. Number of label permutations (default 2000).
#' @param seed Integer. Random seed for reproducibility.
#' @return A named list: observed_auc, p_value, n_perm, method_name.
#' @export
run_permutation_auc_test <- function(y_true, predicted_probs, positive_label,
                                     n_perm = 2000L, seed = 2026) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("[ML] Package 'pROC' required for permutation AUC test.")
  }

  y_bin    <- as.integer(y_true == positive_label)
  roc_obs  <- tryCatch(
    pROC::roc(response = y_bin, predictor = predicted_probs, direction = "<", quiet = TRUE),
    error = function(e) NULL
  )
  if (is.null(roc_obs)) return(list(observed_auc = NA_real_, p_value = NA_real_, n_perm = n_perm))

  obs_auc <- as.numeric(pROC::auc(roc_obs))
  n       <- length(y_bin)

  set.seed(seed)
  perm_aucs <- vapply(seq_len(n_perm), function(b) {
    y_perm <- sample(y_bin, size = n, replace = FALSE)
    roc_p  <- tryCatch(
      pROC::roc(response = y_perm, predictor = predicted_probs, direction = "<", quiet = TRUE),
      error = function(e) NULL
    )
    if (is.null(roc_p)) return(NA_real_)
    as.numeric(pROC::auc(roc_p))
  }, numeric(1))

  perm_aucs <- perm_aucs[!is.na(perm_aucs)]
  p_val     <- (sum(perm_aucs >= obs_auc) + 1L) / (length(perm_aucs) + 1L)

  list(observed_auc = obs_auc, p_value = round(p_val, 5), n_perm = n_perm)
}

#' @title Univariate AUC (Parameter-Free)
#' @description
#' Computes AUC with 95% CI (DeLong) for each marker independently using the
#' Z-score as predictor on all n samples. Direction is auto-selected (max AUC).
#' No hyperparameter tuning and no threshold selection occur, so no LOOCV is
#' needed — the rank-based AUC estimate is unbiased by construction.
#'
#' @param X_main Numeric matrix. Main-effect columns only (no interaction terms).
#' @param y Factor. Binary outcome.
#' @param positive_label String. The positive class label.
#' @return A data.frame with columns: Marker, AUC, CI_Lower, CI_Upper.
#' @export
run_univariate_auc <- function(X_main, y, positive_label) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("[ML] Package 'pROC' required for univariate AUC.")
  }

  y_bin   <- as.integer(y == positive_label)
  markers <- colnames(X_main)

  results <- purrr::map_dfr(markers, function(mk) {
    x_vec <- as.numeric(X_main[, mk])

    roc_fwd <- tryCatch(
      pROC::roc(response = y_bin, predictor = x_vec, direction = "<", quiet = TRUE),
      error = function(e) NULL
    )
    roc_rev <- tryCatch(
      pROC::roc(response = y_bin, predictor = x_vec, direction = ">", quiet = TRUE),
      error = function(e) NULL
    )

    auc_fwd <- if (!is.null(roc_fwd)) as.numeric(pROC::auc(roc_fwd)) else 0
    auc_rev <- if (!is.null(roc_rev)) as.numeric(pROC::auc(roc_rev)) else 0
    roc_obj <- if (auc_fwd >= auc_rev) roc_fwd else roc_rev
    auc_val <- max(auc_fwd, auc_rev)

    ci <- tryCatch(
      as.numeric(pROC::ci.auc(roc_obj, method = "delong")),
      error = function(e) c(NA_real_, auc_val, NA_real_)
    )

    data.frame(
      Marker   = mk,
      AUC      = round(auc_val, 4),
      CI_Lower = round(ci[1],   4),
      CI_Upper = round(ci[3],   4),
      stringsAsFactors = FALSE
    )
  })

  results
}

#' @title Univariate LOO Threshold Classification
#' @description
#' For each marker, runs outer LOO: selects the Youden-J optimal threshold on
#' the training fold ROC and applies it blindly to the held-out sample.
#' Threshold is the only learned parameter; LOO guarantees unbiased evaluation.
#'
#' @param X_main Numeric matrix. Main-effect columns only (no interaction terms).
#' @param y Factor. Binary outcome.
#' @param positive_label String. The positive class label.
#' @param seed Integer. Random seed (for reproducibility of any ties).
#' @return A data.frame with columns: Marker, AUC, Balanced_Accuracy, BER,
#'   Sensitivity, Specificity.
#' @export
run_univariate_loo_threshold <- function(X_main, y, positive_label, seed = 2026) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("[ML] Package 'pROC' required for LOO threshold analysis.")
  }

  y_bin   <- as.integer(y == positive_label)
  markers <- colnames(X_main)
  n       <- nrow(X_main)

  results <- purrr::map_dfr(markers, function(mk) {
    x_vec    <- as.numeric(X_main[, mk])
    preds    <- numeric(n)

    for (i in seq_len(n)) {
      x_train <- x_vec[-i]
      y_train <- y_bin[-i]

      roc_train <- tryCatch(
        pROC::roc(response = y_train, predictor = x_train, quiet = TRUE),
        error = function(e) NULL
      )

      if (is.null(roc_train)) { preds[i] <- 0.5; next }

      coords_df <- tryCatch(
        pROC::coords(roc_train, x = "best", best.method = "youden",
                     ret = "threshold", transpose = FALSE),
        error = function(e) NULL
      )
      thresh <- if (!is.null(coords_df) && nrow(coords_df) > 0) {
        as.numeric(coords_df$threshold[1])
      } else {
        median(x_train)
      }

      preds[i] <- if (x_vec[i] >= thresh) 1L else 0L
    }

    metrics <- compute_classification_metrics(y, as.numeric(preds), positive_label)

    data.frame(
      Marker            = mk,
      AUC               = round(metrics$auc,               4),
      Balanced_Accuracy = round(metrics$balanced_accuracy, 4),
      BER               = round(metrics$ber,               4),
      Sensitivity       = round(metrics$sensitivity,       4),
      Specificity       = round(metrics$specificity,       4),
      stringsAsFactors  = FALSE
    )
  })

  results
}

#' @title Fully-Nested LOO Validation for SVM-RBF Classifier
#' @description
#' For each outer LOO fold, LMM-based feature selection (FDR + LOO sensitivity +
#' collinearity filter) is re-run on the n-1 training patients using the
#' longitudinal dataset, then an SVM-RBF is trained on T0 cross-sectional data
#' and tested on the held-out patient. This provides an unbiased AUC estimate
#' and quantifies the stability of the feature gate across folds.
#'
#' @param DATA_T0         Named list. Step 01 standard (T0) processed data object.
#' @param DATA_LONG       Named list. Step 01 longitudinal processed data object.
#' @param fdr_thr         Numeric. FDR threshold for LMM gate (default 0.05).
#' @param loo_thr         Numeric. LOO p-value threshold for LMM gate (default 0.05).
#' @param cor_threshold   Numeric. Collinearity filter threshold (default 0.85).
#' @param C_grid          Numeric vector. SVM cost values for inner CV.
#' @param gamma_grid      Numeric vector. SVM gamma values for inner CV.
#' @param inner_folds     Integer. Number of inner CV folds (default 5).
#' @param seed            Integer. Random seed for reproducibility.
#' @return Named list: metrics, gate_stability (data.frame), n_empty_folds, predicted_probs.
#' @export
run_nested_loocv_svm_validated <- function(DATA_T0, DATA_LONG,
                                           fdr_thr        = 0.05,
                                           loo_thr        = 0.05,
                                           cor_threshold  = 0.85,
                                           C_grid         = c(0.01, 0.1, 1, 10, 100),
                                           gamma_grid     = c(0.01, 0.1, 1, 10),
                                           inner_folds    = 5L,
                                           seed           = 2026L,
                                           positive_label = NULL) {

  META_COLS   <- c("Patient_ID", "Sample_ID", "Timepoint", "Group")
  std_markers <- setdiff(colnames(DATA_T0$hybrid_data_z),   META_COLS)
  lng_markers <- setdiff(colnames(DATA_LONG$hybrid_data_z), META_COLS)
  shared_mkrs <- intersect(std_markers, lng_markers)

  data_z_t0 <- as.matrix(DATA_T0$hybrid_data_z[, shared_mkrs, drop = FALSE])
  meta_t0   <- DATA_T0$metadata

  # Preserve existing factor level ordering (non-responder first, responder second)
  grp_fac   <- if (is.factor(meta_t0$Group)) meta_t0$Group
               else factor(meta_t0$Group)
  grp_lvls  <- levels(grp_fac)
  pos_label <- if (!is.null(positive_label)) positive_label else grp_lvls[2]
  neg_label <- setdiff(grp_lvls, pos_label)[1]

  pid_t0 <- meta_t0$Patient_ID
  n_t0   <- nrow(meta_t0)

  data_z_lng <- as.data.frame(DATA_LONG$hybrid_data_z[, shared_mkrs, drop = FALSE])
  meta_lng   <- DATA_LONG$metadata

  probs         <- numeric(n_t0)
  gate_per_fold <- vector("list", n_t0)

  for (i in seq_len(n_t0)) {
    pid_test  <- pid_t0[i]
    idx_train <- setdiff(seq_len(n_t0), i)
    X_train   <- data_z_t0[idx_train, , drop = FALSE]
    # Use group labels directly so SVM probability columns match pos_label / neg_label
    y_train   <- as.character(grp_fac[idx_train])
    X_test    <- data_z_t0[i, , drop = FALSE]

    lng_mask     <- meta_lng$Patient_ID != pid_test
    df_lng_train <- cbind(as.data.frame(meta_lng[lng_mask, ]),
                          as.data.frame(data_z_lng[lng_mask, , drop = FALSE]))

    gate <- select_gate_for_fold(df_lng_train, shared_mkrs, fdr_thr, loo_thr)

    if (length(gate) > 1) {
      filt <- filter_collinear_features(X_train[, gate, drop = FALSE], cor_threshold)
      gate <- filt$kept
    }
    gate_per_fold[[i]] <- gate

    if (length(gate) == 0) { probs[i] <- 0.5; next }

    X_tr   <- X_train[, gate, drop = FALSE]
    X_te   <- X_test[,  gate, drop = FALSE]
    mu     <- colMeans(X_tr, na.rm = TRUE)
    sds    <- apply(X_tr, 2, sd, na.rm = TRUE); sds[sds == 0] <- 1
    X_tr_s <- scale(X_tr, center = mu, scale = sds)
    X_te_s <- scale(X_te, center = mu, scale = sds)

    y_tr_fac   <- factor(y_train, levels = grp_lvls)
    safe_folds <- max(2L, min(inner_folds, min(table(y_train))))
    set.seed(seed + i)
    tune_res <- tryCatch(
      e1071::tune(e1071::svm, train.x = X_tr_s, train.y = y_tr_fac,
                  kernel = "radial", probability = TRUE,
                  ranges      = list(cost = C_grid, gamma = gamma_grid),
                  tunecontrol = e1071::tune.control(sampling = "cross", cross = safe_folds)),
      error = function(e) NULL
    )
    if (is.null(tune_res)) { probs[i] <- 0.5; next }

    fit <- tryCatch(
      e1071::svm(X_tr_s, y_tr_fac, kernel = "radial",
                 cost  = tune_res$best.parameters$cost,
                 gamma = tune_res$best.parameters$gamma,
                 probability = TRUE),
      error = function(e) NULL
    )
    if (is.null(fit)) { probs[i] <- 0.5; next }

    pred_obj <- predict(fit, X_te_s, probability = TRUE)
    pm       <- attr(pred_obj, "probabilities")
    probs[i] <- if (pos_label %in% colnames(pm)) pm[, pos_label] else 0.5
  }

  metrics <- compute_classification_metrics(grp_fac, probs, pos_label)

  all_gates <- unlist(gate_per_fold)
  freq_tab  <- sort(table(all_gates), decreasing = TRUE)
  gate_stab <- data.frame(
    Marker         = names(freq_tab),
    Folds_Selected = as.integer(freq_tab),
    Total_Folds    = n_t0,
    Pct            = round(100 * as.integer(freq_tab) / n_t0, 1),
    stringsAsFactors = FALSE
  )

  list(metrics        = metrics,
       gate_stability  = gate_stab,
       n_empty_folds   = sum(sapply(gate_per_fold, length) == 0L),
       predicted_probs = probs)
}

#' Run Fully-Nested LOOCV with Univariate Gate
#'
#' Outer loop: Leave-One-Out (n folds). Inside each fold:
#'   1. Optional per-fold collinearity filter on X_train (skipped when cor_threshold >= 1).
#'   2. Wilcoxon + BH gate on (filtered) X_train vs y_train.
#'   3. Elastic Net (inner k-fold CV for alpha/lambda) on gate features.
#'   4. SVM-RBF (inner k-fold CV for C/gamma) on gate features.
#'   5. Predict on held-out patient with both models.
#'
#' X_all should be globally collinearity-filtered before calling this function.
#' Pass cor_threshold = 1.0 (default) to skip the per-fold filter.
#'
#' Gate stability is computed post-loop: percentage of folds each marker is selected.
#'
#' @param X_all Numeric matrix. Full feature matrix (n x p), pre-collinearity-filtered.
#' @param y Factor. Binary outcome (two levels; positive class = second level by convention).
#' @param positive_label String. Positive class label.
#' @param fdr_threshold Numeric. BH-FDR gate threshold (default 0.20).
#' @param fallback_k Integer. Fallback gate size if fewer than 2 pass FDR (default 3).
#' @param cor_threshold Numeric. Per-fold collinearity threshold (default 1.0, i.e. skip).
#' @param alpha_grid Numeric vector. Elastic Net alpha values for inner CV.
#' @param n_lambda Integer. Lambda path length for glmnet.
#' @param inner_folds Integer. k for inner CV (Elastic Net and SVM).
#' @param C_grid Numeric vector. SVM cost values.
#' @param gamma_grid Numeric vector. SVM gamma values.
#' @param seed Integer. Random seed for reproducibility.
#' @return Named list:
#'   glmnet_probs, svm_probs — out-of-fold predicted probabilities (length n);
#'   y_true — factor of true labels;
#'   positive_label — character;
#'   coef_matrix — NULL (gate varies per fold; cross-fold comparison not meaningful);
#'   gate_log — list of length n, each entry: (markers, used_fallback, pvalues);
#'   n_folds_fallback — count of folds that used the top-k fallback;
#'   gate_stability — data.frame (Marker, Pct_Folds_Selected, Median_Raw_P).
run_nested_loocv_univariate_gate <- function(X_all,
                                              y,
                                              positive_label,
                                              fdr_threshold = 0.20,
                                              fallback_k    = 3L,
                                              cor_threshold = 1.0,
                                              alpha_grid    = c(0, 0.5, 1),
                                              n_lambda      = 100L,
                                              inner_folds   = 5L,
                                              C_grid        = c(0.01, 0.1, 1, 10, 100),
                                              gamma_grid    = c(0.01, 0.1, 1, 10),
                                              seed          = 2026L) {

  if (!requireNamespace("glmnet", quietly = TRUE))
    stop("[ML] Package 'glmnet' required for run_nested_loocv_univariate_gate().")
  if (!requireNamespace("e1071", quietly = TRUE))
    stop("[ML] Package 'e1071' required for run_nested_loocv_univariate_gate().")

  n           <- nrow(X_all)
  pos_label   <- positive_label
  y_bin_full  <- as.integer(y == pos_label)
  all_markers <- colnames(X_all)

  glmnet_probs <- numeric(n)
  svm_probs    <- numeric(n)
  gate_log     <- vector("list", n)

  message(sprintf(
    "   [Univariate Gate] Outer LOO: %d folds | FDR<%.2f | fallback_k=%d | inner_folds=%d",
    n, fdr_threshold, as.integer(fallback_k), as.integer(inner_folds)
  ))

  for (i in seq_len(n)) {

    X_train_full <- X_all[-i, , drop = FALSE]
    y_train_fac  <- y[-i]
    y_train_bin  <- y_bin_full[-i]
    X_test_full  <- X_all[i, , drop = FALSE]

    # Step 1 — Per-fold collinearity filter (only when explicitly requested)
    if (cor_threshold < 1.0 && ncol(X_train_full) >= 2L) {
      filt      <- filter_collinear_features(X_train_full, cor_threshold)
      X_tr_filt <- filt$X
    } else {
      X_tr_filt <- X_train_full
    }

    # Step 2 — Wilcoxon + BH gate on training fold
    gate_i    <- select_univariate_gate(
      X_train        = X_tr_filt,
      y_train        = y_train_fac,
      positive_label = pos_label,
      fdr_threshold  = fdr_threshold,
      fallback_k     = fallback_k
    )
    gate_log[[i]] <- gate_i

    gate_cols <- gate_i$markers
    if (length(gate_cols) == 0L) {
      glmnet_probs[i] <- 0.5
      svm_probs[i]    <- 0.5
      next
    }

    X_tr_gate <- X_tr_filt[, gate_cols, drop = FALSE]
    X_te_gate <- X_test_full[, gate_cols, drop = FALSE]

    n_min_class  <- min(sum(y_train_bin == 0L), sum(y_train_bin == 1L))
    safe_folds_k <- max(2L, min(as.integer(inner_folds), n_min_class))

    # Step 3 — Elastic Net: inner k-fold CV over alpha_grid x lambda path
    best_alpha     <- alpha_grid[1L]
    best_lambda    <- NULL
    best_inner_auc <- -Inf

    for (a in alpha_grid) {
      set.seed(seed + i)
      cv_fit <- tryCatch(
        glmnet::cv.glmnet(
          x            = X_tr_gate,
          y            = as.numeric(y_train_bin),
          family       = "binomial",
          alpha        = a,
          nfolds       = safe_folds_k,
          type.measure = "auc",
          nlambda      = as.integer(n_lambda)
        ),
        error = function(e) NULL
      )
      if (!is.null(cv_fit)) {
        inner_auc <- suppressWarnings(max(cv_fit$cvm, na.rm = TRUE))
        if (is.finite(inner_auc) && inner_auc > best_inner_auc) {
          best_inner_auc <- inner_auc
          best_alpha     <- a
          best_lambda    <- cv_fit$lambda.min
        }
      }
    }

    final_en <- tryCatch(
      glmnet::glmnet(X_tr_gate, as.numeric(y_train_bin),
                     family = "binomial",
                     alpha  = best_alpha,
                     lambda = best_lambda),
      error = function(e) NULL
    )

    glmnet_probs[i] <- if (!is.null(final_en) && !is.null(best_lambda)) {
      tryCatch(
        as.numeric(predict(final_en, newx = X_te_gate,
                           s = best_lambda, type = "response")),
        error = function(e) 0.5
      )
    } else 0.5

    # Step 4 — SVM-RBF: grid search over C_grid x gamma_grid
    set.seed(seed + i)
    tune_res <- tryCatch(
      suppressWarnings(e1071::tune(
        e1071::svm,
        train.x     = X_tr_gate,
        train.y     = y_train_fac,
        kernel      = "radial",
        ranges      = list(cost = C_grid, gamma = gamma_grid),
        tunecontrol = e1071::tune.control(sampling = "cross",
                                          cross    = safe_folds_k)
      )),
      error = function(e) NULL
    )

    svm_probs[i] <- if (!is.null(tune_res)) {
      bp <- tune_res$best.parameters
      set.seed(seed + i)
      final_svm <- tryCatch(
        e1071::svm(x = X_tr_gate, y = y_train_fac,
                   kernel      = "radial",
                   cost        = bp$cost,
                   gamma       = bp$gamma,
                   probability = TRUE),
        error = function(e) NULL
      )
      if (!is.null(final_svm)) {
        pred_obj <- tryCatch(
          predict(final_svm, newdata = X_te_gate, probability = TRUE),
          error = function(e) NULL
        )
        if (!is.null(pred_obj)) {
          pm <- attr(pred_obj, "probabilities")
          if (!is.null(pm) && pos_label %in% colnames(pm)) {
            as.numeric(pm[1L, pos_label])
          } else 0.5
        } else 0.5
      } else 0.5
    } else 0.5
  }

  # Post-loop: gate stability aggregation
  stability_pct <- vapply(all_markers, function(m) {
    mean(vapply(gate_log, function(g) m %in% g$markers, logical(1L)))
  }, numeric(1L))

  median_pvals <- vapply(all_markers, function(m) {
    pv <- vapply(gate_log, function(g) {
      if (m %in% names(g$pvalues)) g$pvalues[[m]] else NA_real_
    }, numeric(1L))
    median(pv, na.rm = TRUE)
  }, numeric(1L))

  n_folds_fallback <- sum(vapply(gate_log, function(g) isTRUE(g$used_fallback), logical(1L)))

  gate_stability_df <- data.frame(
    Marker             = all_markers,
    Pct_Folds_Selected = round(stability_pct * 100, 1),
    Median_Raw_P       = round(median_pvals, 4),
    stringsAsFactors   = FALSE
  )
  gate_stability_df <- gate_stability_df[order(-gate_stability_df$Pct_Folds_Selected), ]
  rownames(gate_stability_df) <- NULL

  list(
    glmnet_probs     = glmnet_probs,
    svm_probs        = svm_probs,
    y_true           = y,
    positive_label   = pos_label,
    coef_matrix      = NULL,          # variable gate per fold — cross-fold comparison not meaningful
    gate_log         = gate_log,
    n_folds_fallback = n_folds_fallback,
    gate_stability   = gate_stability_df
  )
}


# ==============================================================================
# Gate Signal Decomposition
# ==============================================================================

