# R/modules_ml_plots.R
# ==============================================================================
# MACHINE LEARNING MODULE — PLOTS
# Split from the former monolithic R/modules_ml.R (sourced via that aggregator).
# Dependencies: dplyr, ggplot2, tidyr, glmnet, e1071, pROC
# ==============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)

#' @title Plot Nested-LOOCV ROC Curves
#' @description
#' Overlays ROC curves for multiple methods on a single ggplot with AUC annotated
#' in the legend. ML model curves are drawn as solid lines; optional clinical
#' benchmark curves are drawn as dashed lines at a distinct colour to signal their
#' different nature (univariate biomarker, potentially on a smaller subset).
#'
#' @param results_list Named list of ML result objects (each with predicted_probs,
#'   y_true, positive_label, method).
#' @param colors_viz Named character vector of clinical group colours (optional).
#' @param title String. Plot title.
#' @param benchmark_list Optional list of benchmark result objects from
#'   run_clinical_benchmark(). Each element must contain predicted_probs, y_true,
#'   positive_label, label, n_valid, and direction fields.
#' @return A ggplot object, or NULL on failure.
#' @export
plot_ml_roc <- function(results_list,
                        colors_viz     = NULL,
                        title          = "Nested-LOOCV ROC Curves",
                        benchmark_list = NULL) {
  if (!requireNamespace("pROC", quietly = TRUE)) return(NULL)

  method_colors    <- c("Elastic Net" = "#2171B5", "SVM-RBF" = "#CB181D")
  benchmark_colors <- c("#6A3D9A", "#FF7F00", "#33A02C")

  # ── ML model ROC data (solid curves) ────────────────────────────────────────
  roc_data <- purrr::map_dfr(names(results_list), function(nm) {
    res   <- results_list[[nm]]
    y_bin <- as.integer(res$y_true == res$positive_label)
    tryCatch({
      roc_obj <- pROC::roc(response = y_bin, predictor = res$predicted_probs,
                           direction = "<", quiet = TRUE)
      auc_val <- round(as.numeric(pROC::auc(roc_obj)), 3)
      data.frame(
        FPR    = 1 - roc_obj$specificities,
        TPR    = roc_obj$sensitivities,
        Method = sprintf("%s\n(AUC = %.3f)", nm, auc_val),
        stringsAsFactors = FALSE
      )
    }, error = function(e) data.frame())
  })

  if (nrow(roc_data) == 0) return(NULL)

  # ── Benchmark ROC data (dashed curves, possibly on a subset) ────────────────
  bench_data <- data.frame()
  if (!is.null(benchmark_list)) {
    bench_data <- purrr::map_dfr(seq_along(benchmark_list), function(bi) {
      br    <- benchmark_list[[bi]]
      y_bin <- as.integer(br$y_true == br$positive_label)
      x_vec <- as.numeric(br$predicted_probs)
      dir   <- if (!is.null(br$direction)) br$direction else "<"
      tryCatch({
        roc_obj <- pROC::roc(response = y_bin, predictor = x_vec,
                             direction = dir, quiet = TRUE)
        auc_val <- round(as.numeric(pROC::auc(roc_obj)), 3)
        data.frame(
          FPR    = 1 - roc_obj$specificities,
          TPR    = roc_obj$sensitivities,
          Method = sprintf("%s\n(AUC = %.3f, n=%d)", br$label, auc_val, br$n_valid),
          stringsAsFactors = FALSE
        )
      }, error = function(e) data.frame())
    })
  }

  # ── Colour map ───────────────────────────────────────────────────────────────
  # ML methods: keyed by short name prefix; benchmarks: cycle through palette.
  color_map <- stats::setNames(
    sapply(unique(roc_data$Method), function(m) {
      key <- sub("\n.*$", "", m)
      if (key %in% names(method_colors)) method_colors[[key]] else "#636363"
    }),
    unique(roc_data$Method)
  )
  if (nrow(bench_data) > 0) {
    bench_methods <- unique(bench_data$Method)
    color_map[bench_methods] <- benchmark_colors[
      ((seq_along(bench_methods) - 1L) %% length(benchmark_colors)) + 1L
    ]
  }

  # ── Build plot: two separate geom_line layers preserve linetype clarity ──────
  p <- ggplot2::ggplot(mapping = ggplot2::aes(x = FPR, y = TPR, color = Method)) +
    ggplot2::geom_line(data = roc_data, linewidth = 1.2, na.rm = TRUE)

  if (nrow(bench_data) > 0) {
    p <- p + ggplot2::geom_line(data = bench_data, linewidth = 1.0,
                                linetype = "dashed", na.rm = TRUE)
  }

  p <- p +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                         color = "grey55", linewidth = 0.5) +
    ggplot2::scale_color_manual(values = color_map, name = NULL) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0),
                                 breaks = seq(0, 1, 0.2)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0),
                                 breaks = seq(0, 1, 0.2)) +
    ggplot2::labs(
      title    = title,
      x        = "1 - Specificity (FPR)",
      y        = "Sensitivity (TPR)",
      color    = NULL,
      caption  = if (nrow(bench_data) > 0) "Dashed: clinical benchmark (univariate, available-cases subset)" else NULL
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position  = "bottom",
      legend.direction = "vertical",
      plot.title       = ggplot2::element_text(size = 11, face = "bold"),
      plot.caption     = ggplot2::element_text(size = 8, color = "grey45")
    )

  p
}

#' @title Plot ML Predicted Probabilities by Clinical Group
#' @description
#' Jitter strip plot of out-of-fold predicted probabilities faceted by method.
#' A dashed reference line at 0.5 indicates the neutral decision boundary.
#'
#' @param predictions_df Data.frame with columns: Group, Prob_ElasticNet, Prob_SVM_RBF.
#' @param colors_viz Named character vector mapping group labels to hex colours.
#' @param title String. Plot title.
#' @return A ggplot object, or NULL on failure.
#' @export
plot_ml_predictions <- function(predictions_df,
                                colors_viz = NULL,
                                title      = "Predicted Probabilities by Clinical Group") {
  tryCatch({
    df_long <- predictions_df %>%
      dplyr::select(Group, dplyr::starts_with("Prob_")) %>%
      tidyr::pivot_longer(
        cols      = dplyr::starts_with("Prob_"),
        names_to  = "Method",
        values_to = "Probability"
      ) %>%
      dplyr::mutate(
        Method = dplyr::recode(Method,
          "Prob_ElasticNet" = "Elastic Net",
          "Prob_SVM_RBF"    = "SVM-RBF"
        ),
        Group = factor(Group)
      )

    p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Group, y = Probability, color = Group)) +
      ggplot2::geom_jitter(width = 0.18, alpha = 0.80, size = 2.2, na.rm = TRUE) +
      ggplot2::stat_summary(fun = median, geom = "crossbar", width = 0.4,
                            color = "grey25", linewidth = 0.5, fatten = 1.5) +
      ggplot2::facet_wrap(~Method) +
      ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed",
                          color = "grey45", linewidth = 0.6) +
      ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
      ggplot2::labs(title = title, x = NULL,
                    y = sprintf("Predicted P(Positive Class)")) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(legend.position = "none",
                     plot.title = ggplot2::element_text(size = 11, face = "bold"))

    if (!is.null(colors_viz) && length(colors_viz) > 0) {
      avail <- intersect(levels(df_long$Group), names(colors_viz))
      if (length(avail) == length(levels(df_long$Group))) {
        p <- p + ggplot2::scale_color_manual(values = colors_viz)
      }
    }

    p
  }, error = function(e) {
    warning(paste("[ML] plot_ml_predictions failed:", e$message))
    return(NULL)
  })
}

#' @title Plot Univariate ROC Curves
#' @description
#' Overlaid ROC curves for each marker from the univariate AUC analysis.
#' Uses the same visual style as plot_ml_roc() (bw theme, diagonal reference,
#' AUC in legend). Reconstructs ROC objects from the full dataset Z-scores.
#'
#' @param X_main Numeric matrix. Main-effect columns (same used in run_univariate_auc).
#' @param y Factor. Binary outcome.
#' @param positive_label String. The positive class label.
#' @param univariate_df Data.frame returned by run_univariate_auc() (for AUC labels).
#' @param colors_viz Named character vector (optional, for marker colours).
#' @param title String. Plot title.
#' @return A ggplot object, or NULL on failure.
#' @export
plot_univariate_roc <- function(X_main, y, positive_label, univariate_df,
                                colors_viz = NULL,
                                title      = "Univariate ROC Curves") {
  if (!requireNamespace("pROC", quietly = TRUE)) return(NULL)

  y_bin   <- as.integer(y == positive_label)
  markers <- colnames(X_main)

  roc_data <- purrr::map_dfr(markers, function(mk) {
    x_vec <- as.numeric(X_main[, mk])

    roc_fwd <- tryCatch(pROC::roc(response = y_bin, predictor = x_vec, direction = "<", quiet = TRUE), error = function(e) NULL)
    roc_rev <- tryCatch(pROC::roc(response = y_bin, predictor = x_vec, direction = ">", quiet = TRUE), error = function(e) NULL)
    auc_fwd <- if (!is.null(roc_fwd)) as.numeric(pROC::auc(roc_fwd)) else 0
    auc_rev <- if (!is.null(roc_rev)) as.numeric(pROC::auc(roc_rev)) else 0
    roc_obj <- if (auc_fwd >= auc_rev) roc_fwd else roc_rev
    auc_val <- max(auc_fwd, auc_rev)

    if (is.null(roc_obj)) return(data.frame())

    data.frame(
      FPR    = 1 - roc_obj$specificities,
      TPR    = roc_obj$sensitivities,
      Marker = sprintf("%s (AUC = %.3f)", mk, auc_val),
      stringsAsFactors = FALSE
    )
  })

  if (nrow(roc_data) == 0) return(NULL)

  p <- ggplot2::ggplot(roc_data, ggplot2::aes(x = FPR, y = TPR, color = Marker)) +
    ggplot2::geom_line(linewidth = 1.1, na.rm = TRUE) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         color = "grey55", linewidth = 0.6) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0),
                                 breaks = seq(0, 1, 0.2)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0),
                                 breaks = seq(0, 1, 0.2)) +
    ggplot2::labs(title = title, x = "1 - Specificity (FPR)",
                  y = "Sensitivity (TPR)", color = NULL) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position  = "bottom",
      legend.direction = "vertical",
      plot.title       = ggplot2::element_text(size = 11, face = "bold")
    )

  p
}


# ==============================================================================
# NESTED LOO VALIDATION — fully-nested LMM feature selection inside outer fold
# ==============================================================================

#' @title Plot Benchmark-Stratified Subgroup Analysis (2-panel)
#' @description
#' Renders a publication-ready 2-panel figure for run_pdl1_stratified() output:
#'   - Panel A: bar chart of response rate by benchmark bin (with N per bin
#'     labelled on bars) plus the chi-squared trend test p-value
#'   - Panel B: ROC curve of the model restricted to the "benchmark-low"
#'     subgroup (typically the clinically actionable population) with AUC and
#'     95% CI annotated. A single point at the benchmark binary-cut performance
#'     in the same subgroup is overlaid as a reference, when defined.
#'
#' Uses patchwork to lay out the two panels side by side. Pass the same
#' colors_viz used elsewhere to keep palette consistency.
#'
#' @param stratified_result List returned by run_pdl1_stratified().
#' @param colors_viz Named vector of clinical colours (optional; if missing,
#'   responder/non-responder default to standard journal palette).
#' @param positive_label Positive class label (e.g. "RP"). Defaults to inferred.
#' @param title Optional overall title for the figure.
#' @return A patchwork object combining the two panels, or NULL on failure.
#' @export
plot_benchmark_stratified <- function(stratified_result,
                                      colors_viz     = NULL,
                                      positive_label = NULL,
                                      title          = NULL) {
  if (is.null(stratified_result)) return(NULL)
  if (!requireNamespace("ggplot2",   quietly = TRUE)) return(NULL)
  if (!requireNamespace("patchwork", quietly = TRUE)) return(NULL)
  if (!requireNamespace("pROC",      quietly = TRUE)) return(NULL)

  pos_lbl <- if (!is.null(positive_label)) positive_label else "Responder"
  resp_color  <- if (!is.null(colors_viz) && pos_lbl %in% names(colors_viz)) colors_viz[[pos_lbl]] else "#2E8B57"
  nonresp_color <- if (!is.null(colors_viz)) {
    cn <- setdiff(names(colors_viz), pos_lbl)
    if (length(cn) > 0) colors_viz[[cn[1]]] else "#B2182B"
  } else "#B2182B"

  # Short label for titles — strip trailing "expression (%)" / "(%)" noise
  bench_short <- sub("\\s*(expression\\s*\\(%\\)|\\(%\\)|\\(.*\\))\\s*$", "",
                     stratified_result$label)
  bench_short <- trimws(bench_short)
  if (bench_short == "") bench_short <- stratified_result$label

  # ---- Panel A: Response rate bar chart ----
  bin_df <- stratified_result$bin_crosstab
  bin_df$Bin <- factor(bin_df$Bin, levels = stratified_result$bin_labels)
  bin_df$Label <- sprintf("%.1f%%\n(N=%d, R=%d)",
                          bin_df$Response_Rate, bin_df$N, bin_df$N_Responder)

  pA <- ggplot2::ggplot(bin_df, ggplot2::aes(x = Bin, y = Response_Rate)) +
    ggplot2::geom_col(fill = resp_color, alpha = 0.85, width = 0.65) +
    ggplot2::geom_text(ggplot2::aes(label = Label),
                       vjust = -0.4, size = 4.2, lineheight = 0.85) +
    ggplot2::geom_hline(yintercept = 50, linetype = "dashed", color = "grey60") +
    ggplot2::scale_y_continuous(limits = c(0, max(bin_df$Response_Rate) + 18),
                                 expand = c(0, 0),
                                 breaks = seq(0, 100, 25),
                                 labels = function(x) paste0(x, "%")) +
    ggplot2::labs(
      title = sprintf("Response rate by %s strata", bench_short),
      subtitle = sprintf("3-bin Fisher p = %.3f | Cochran-Armitage trend p = %.3f | N total = %d",
                         stratified_result$fisher_3bins_p,
                         stratified_result$ca_trend_p,
                         stratified_result$n_valid),
      x = stratified_result$label, y = "Responder rate"
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "gray35", size = 11),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.text.x   = ggplot2::element_text(face = "bold")
    )

  # ---- Panel B: ROC in benchmark-low subgroup ----
  sl <- stratified_result$subgroup_low
  pB <- NULL
  if (!is.null(sl$predicted_probs) && !is.na(sl$auc)) {
    y_bin   <- as.integer(sl$y_true == pos_lbl)
    roc_obj <- tryCatch(pROC::roc(y_bin, sl$predicted_probs, direction = "<", quiet = TRUE),
                        error = function(e) NULL)
    if (!is.null(roc_obj)) {
      roc_df <- data.frame(
        FPR = 1 - roc_obj$specificities,
        TPR = roc_obj$sensitivities
      )

      # Binary-cut reference point (TPS>=threshold within the full cohort)
      bc <- stratified_result$binary_cut
      bc_x <- 1 - bc$specificity; bc_y <- bc$sensitivity

      auc_lab <- sprintf("KI67-gate model (subgroup)\nAUC = %.3f [%.3f-%.3f]",
                         sl$auc, sl$auc_ci[1], sl$auc_ci[3])

      bench_threshold <- bc$threshold
      bench_lab <- sprintf("%s >=%g cut\n(full cohort)", bench_short, bench_threshold)

      pB <- ggplot2::ggplot(roc_df, ggplot2::aes(x = FPR, y = TPR)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                              color = "grey55", linewidth = 0.6) +
        ggplot2::geom_line(color = resp_color, linewidth = 1.3) +
        ggplot2::annotate("point",  x = bc_x, y = bc_y, color = nonresp_color,
                          shape = 18, size = 4.5) +
        ggplot2::annotate("text",   x = bc_x + 0.04, y = bc_y - 0.02,
                          label = bench_lab, hjust = 0, size = 3.4,
                          color = nonresp_color, lineheight = 0.9) +
        ggplot2::annotate("label",  x = 0.98, y = 0.06, label = auc_lab,
                          hjust = 1, size = 3.8, color = resp_color, fill = "white",
                          label.size = 0.4, fontface = "bold", lineheight = 0.9) +
        ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0.005, 0)) +
        ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0.005, 0)) +
        ggplot2::labs(
          title = sprintf("Model in %s < %g subgroup",
                          bench_short, bc$threshold),
          subtitle = sprintf("n = %d (responders=%d) | clinical context where %s alone is non-discriminative",
                             sl$n, sl$n_responder, bench_short),
          x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)"
        ) +
        ggplot2::theme_bw(base_size = 13) +
        ggplot2::theme(
          plot.title    = ggplot2::element_text(face = "bold"),
          plot.subtitle = ggplot2::element_text(color = "gray35", size = 11),
          panel.grid.minor = ggplot2::element_blank()
        )
    }
  }

  if (is.null(pB)) {
    out <- pA
  } else {
    out <- patchwork::wrap_plots(pA, pB, ncol = 2, widths = c(1, 1))
  }
  if (!is.null(title)) {
    out <- out + patchwork::plot_annotation(
      title = title,
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 14))
    )
  }
  out
}

#' @title Plot Combined-Model Information Gain (2-panel)
#' @description
#' Publication figure for run_combined_benchmark_model(). Panel A overlays the
#' three ROC curves of interest on the same n_valid subset (benchmark alone,
#' benchmark + LMM gate combined via logistic LOOCV, and the LMM-gate SVM-RBF
#' nested-LOOCV). Panel B is a forest-style plot of the three information-gain
#' point estimates (IDI, cNRI, ΔAUC) with their bootstrap 95% CIs, mapping
#' directly onto the Pencina reclassification framework expected by reviewers.
#'
#' @param combined_result List returned by run_combined_benchmark_model().
#' @param title Optional overall title.
#' @return A patchwork object, or NULL on failure.
#' @export
plot_combined_information_gain <- function(combined_result, title = NULL) {
  if (is.null(combined_result) || is.null(combined_result$plot_data)) return(NULL)
  if (!requireNamespace("ggplot2",   quietly = TRUE)) return(NULL)
  if (!requireNamespace("patchwork", quietly = TRUE)) return(NULL)
  if (!requireNamespace("pROC",      quietly = TRUE)) return(NULL)

  pd <- combined_result$plot_data
  lg <- combined_result$logistic
  ig <- combined_result$information_gain
  sv <- combined_result$svm_comparison

  y_bin <- pd$y_bin
  bench_label <- combined_result$label
  bench_short <- sub("\\s*(expression\\s*\\(%\\)|\\(%\\)|\\(.*\\))\\s*$", "", bench_label)
  bench_short <- trimws(bench_short)
  if (bench_short == "") bench_short <- bench_label

  build_roc <- function(p_vec, label, auc_val, auc_ci_lo, auc_ci_hi) {
    roc_obj <- tryCatch(pROC::roc(y_bin, p_vec, direction = "<", quiet = TRUE),
                        error = function(e) NULL)
    if (is.null(roc_obj)) return(NULL)
    data.frame(
      FPR = 1 - roc_obj$specificities,
      TPR = roc_obj$sensitivities,
      Model = sprintf("%s (AUC=%.3f [%.3f-%.3f])",
                      label, auc_val, auc_ci_lo, auc_ci_hi),
      stringsAsFactors = FALSE
    )
  }

  # SVM CIs not stored explicitly; reconstruct quickly for labelling
  roc_only <- pROC::roc(y_bin, pd$p_svm_features_only,       direction = "<", quiet = TRUE)
  ci_only  <- tryCatch(as.numeric(pROC::ci.auc(roc_only, method = "delong")),
                       error = function(e) c(NA, pd$auc_svm_features_only, NA))
  roc_with <- pROC::roc(y_bin, pd$p_svm_features_with_bench, direction = "<", quiet = TRUE)
  ci_with  <- tryCatch(as.numeric(pROC::ci.auc(roc_with, method = "delong")),
                       error = function(e) c(NA, pd$auc_svm_features_with_bench, NA))

  rocs <- rbind(
    build_roc(pd$p_bench_logistic,    sprintf("%s alone (logistic)", bench_short),
              lg$auc_benchmark,
              lg$auc_benchmark_ci[[1]], lg$auc_benchmark_ci[[3]]),
    build_roc(pd$p_combined_logistic, sprintf("%s + KI67 gate (logistic)", bench_short),
              lg$auc_combined,
              lg$auc_combined_ci[[1]], lg$auc_combined_ci[[3]]),
    build_roc(pd$p_svm_features_only, "KI67 gate alone (SVM-RBF)",
              pd$auc_svm_features_only, ci_only[1], ci_only[3])
  )

  model_colors <- c("#7A5CA6", "#2E8B57", "#B2182B")
  names(model_colors) <- unique(rocs$Model)

  pA <- ggplot2::ggplot(rocs, ggplot2::aes(x = FPR, y = TPR, color = Model)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                          color = "grey55", linewidth = 0.6) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_color_manual(values = model_colors) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0.005, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0.005, 0)) +
    ggplot2::labs(
      title    = "ROC: benchmark, combined, and gate-only",
      subtitle = sprintf("Same cohort subset (n = %d, non-NA %s)",
                         combined_result$n_valid, bench_short),
      x        = "1 - Specificity (FPR)", y = "Sensitivity (TPR)",
      color    = NULL
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "gray35", size = 11),
      legend.position = "bottom",
      legend.direction = "vertical",
      legend.text   = ggplot2::element_text(size = 9.5),
      legend.key.height = ggplot2::unit(8, "pt"),
      panel.grid.minor = ggplot2::element_blank()
    )

  # ---- Panel B: Information-gain forest ----
  forest_df <- data.frame(
    Metric = c(
      sprintf("IDI\n(bootstrap p = %.3f)",      ig$IDI_bootstrap_p),
      sprintf("Continuous NRI\n(bootstrap p = %.3f)", ig$cNRI_bootstrap_p),
      sprintf("Delta AUC (logistic)\n(DeLong p = %.3f)", lg$delong_p),
      sprintf("Delta AUC (SVM combined)\n(DeLong p = %.3f)", sv$delong_p)
    ),
    Estimate = c(ig$IDI, ig$cNRI, lg$delta_auc, sv$delta_auc),
    Lower    = c(ig$IDI_95CI[[1]], ig$cNRI_95CI[[1]], NA_real_, NA_real_),
    Upper    = c(ig$IDI_95CI[[2]], ig$cNRI_95CI[[2]], NA_real_, NA_real_),
    stringsAsFactors = FALSE
  )
  forest_df$Metric <- factor(forest_df$Metric, levels = rev(forest_df$Metric))
  forest_df$Significant <- ifelse(
    !is.na(forest_df$Lower) & forest_df$Lower > 0, "Positive (CI > 0)",
    ifelse(!is.na(forest_df$Lower) & forest_df$Upper < 0, "Negative (CI < 0)",
           ifelse(forest_df$Estimate > 0, "Positive (point)", "Negative (point)"))
  )
  sig_colors <- c(
    "Positive (CI > 0)"  = "#2E8B57",
    "Negative (CI < 0)"  = "#B2182B",
    "Positive (point)"   = "#7AC480",
    "Negative (point)"   = "#E08983"
  )

  pB <- ggplot2::ggplot(forest_df,
                        ggplot2::aes(x = Estimate, y = Metric, color = Significant)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = Lower, xmax = Upper),
                             height = 0.2, linewidth = 0.9, na.rm = TRUE) +
    ggplot2::geom_point(size = 4.5) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.3f", Estimate)),
                        hjust = -0.25, vjust = -0.7, size = 4, fontface = "bold",
                        color = "black") +
    ggplot2::scale_color_manual(values = sig_colors, na.translate = FALSE) +
    ggplot2::labs(
      title    = sprintf("Information gain (KI67 gate vs %s)", bench_short),
      subtitle = sprintf("Point estimates with 95%% bootstrap CIs (B=%d) where applicable",
                         ig$n_boot),
      x        = "Estimate (positive = KI67 gate adds value)",
      y        = NULL,
      color    = NULL
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "gray35", size = 11),
      legend.position = "bottom",
      axis.text.y   = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = ggplot2::margin(8, 50, 8, 8))

  out <- patchwork::wrap_plots(pA, pB, ncol = 2, widths = c(1.05, 1))
  if (!is.null(title)) {
    out <- out + patchwork::plot_annotation(
      title = title,
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 14))
    )
  }
  out
}

# ==============================================================================
# UNIVARIATE GATE PATH — Cross-sectional feature selection (no T1 required)
# ==============================================================================

