# R/modules_viz.R
# ==============================================================================
# VISUALIZATION MODULE
# Description: Custom themes, PCA wrappers, and Plotting functions.
# Dependencies: ggplot2, FactoMineR, factoextra
# ==============================================================================

library(ggplot2)
library(FactoMineR)
library(factoextra)
library(ggrepel)
library(igraph)
library(tidygraph)
library(ggraph)
library(ComplexHeatmap)
library(circlize)
library(patchwork)

# ---- Shared viz colour constants --------------------------------------------
# Single definition of recurring hex literals (RdBu diverging ramp + semantic
# accent colours), referenced across viz functions in place of copy-pasted hex.
# These are the internal-QC palette; manuscript figures use pub_palette
# (Okabe-Ito) in modules_pub_style.R.
VIZ_BLUE      <- "#2166AC"  # blue: positive-class default / structural connector
VIZ_GREY_MID  <- "#F7F7F7"  # neutral midpoint of the diverging ramp
VIZ_RED       <- "#B2182B"  # firebrick: negative class / treat-all / master regulator
VIZ_RED_LT    <- "#D6604D"  # light red: solo driver
VIZ_DIVERGING <- c(VIZ_BLUE, VIZ_GREY_MID, VIZ_RED)  # RdBu ramp for z-score heatmaps

#' @title Standard Project Theme
#' @description A consistent ggplot2 theme for publication-quality figures.
#' @return A ggplot theme object.
theme_coda <- function() {
  theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "gray95")
    )
}


#' @title Plot Merged Raw Distribution with Highlights & Stats
#' @description 
#' Visualizes raw marker distribution.
#' Features: Violin (density), Jitter (points), Mean (Dashed Red), Median (Solid Black).
#' @param data_df Dataframe containing 'Group', metadata columns, and the marker.
#' @param marker_name Name of the marker column to plot.
#' @param colors Named vector of colors for the Groups.
#' @param highlight_pattern String pattern to search in 'Original_Source' (default: "_LS").
#' @param highlight_color Color for the border of highlighted points (default: "#FFD700").
#' @return A ggplot object.
plot_raw_distribution_merged <- function(data_df, marker_name, colors, 
                                         highlight_pattern = "_LS", 
                                         highlight_color = "#FFD700") {
  
  require(ggplot2)
  require(dplyr)
  
  # 1. Validation
  if (!marker_name %in% names(data_df)) return(NULL)
  
  # 2. Prepare Data & Highlighting Logic
  src_col <- names(data_df)[grepl("Original|Source", names(data_df), ignore.case = TRUE)][1]
  
  plot_data <- data_df
  plot_data$Value <- plot_data[[marker_name]]
  
  # Define Highlight Attributes defaults
  plot_data$Highlight_Type <- "Standard"
  plot_data$Stroke_Size <- 0.2
  
  if (!is.na(src_col)) {
    is_target <- grepl(highlight_pattern, plot_data[[src_col]])
    if (any(is_target)) {
      plot_data$Highlight_Type[is_target] <- "Target"
      plot_data$Stroke_Size[is_target] <- 1.5
    }
  }
  
  # 3. Calculate Stats for Subtitle
  n_na <- sum(is.na(plot_data$Value))
  pct_na <- round((n_na / nrow(plot_data)) * 100, 1)
  
  # 4. Explicitly remove NAs 
  plot_data_clean <- plot_data %>% filter(!is.na(Value))
  
  # 5. Plotting
  p <- ggplot(plot_data_clean, aes(x = Group, y = Value, fill = Group)) +
    
    # Layer 1: Violin
    geom_violin(alpha = 0.4, trim = FALSE, color = NA, scale = "width") +
    
    # Layer 2: Jittered Points
    geom_jitter(aes(color = Highlight_Type, stroke = Stroke_Size), 
                width = 0.2, size = 2.5, shape = 21, alpha = 0.8) +
    
    # Layer 3: Median (Solid Black Line)
    stat_summary(fun = median, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
                 width = 0.5, linewidth = 0.8, color = "black") +
    
    # Layer 4: Mean (Dashed Red Line)
    stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
                 width = 0.5, linewidth = 0.8, color = "darkred", linetype = "dashed") +
    
    # Colors
    scale_fill_manual(values = colors) +
    scale_color_manual(values = c("Standard" = "white", "Target" = highlight_color), 
                       guide = "none") +
    scale_continuous_identity(aesthetics = "stroke") +
    
    # Styling
    labs(
      title = paste("Raw Distribution:", marker_name),
      subtitle = sprintf("Stats: Black Line = Median | Red Dashed = Mean | Highlight: '%s'", 
                         highlight_pattern),
      y = "Raw Value",
      x = NULL
    ) +
    theme_coda() +
    theme(
      legend.position = "none", # Legend redundant as X axis labels exist
      axis.text.x = element_text(face = "bold", size = 11)
    )
  
  return(p)
}

#' @title Save Distribution Report (Multi-Page PDF)
#' @description 
#' Iterates through all markers and saves their raw distribution plots 
#' into a single PDF file. Encapsulates the looping logic.
#' 
#' @param data_df Merged dataframe with metadata and raw values.
#' @param markers Vector of marker names to plot.
#' @param file_path Output path for the PDF.
#' @param colors Named vector of colors for groups.
#' @param hl_pattern String pattern for highlighting (e.g. "_LS").
viz_save_distribution_report <- function(data_df, markers, file_path, colors, hl_pattern = "") {
  
  if (length(markers) == 0) {
    warning("[Viz] No markers provided for distribution report.")
    return(NULL)
  }
  
  message(sprintf("   [Viz] Saving distribution plots to: %s", basename(file_path)))
  
  pdf(file_path, width = 8, height = 6)
  # Ensure device is closed even if the loop crashes
  on.exit(try(dev.off(), silent = TRUE), add = TRUE)
  
  for (marker in markers) {
    if (marker %in% names(data_df)) {
      # Try-catch ensures one bad plot doesn't crash the whole PDF generation
      tryCatch({
        p <- plot_raw_distribution_merged(
          data_df = data_df, 
          marker_name = marker,
          colors = colors,     
          highlight_pattern = hl_pattern,   
          highlight_color = "#FFD700"
        )
        if (!is.null(p)) print(p)
      }, error = function(e) {
        warning(sprintf("      [WARN] Failed to plot marker '%s': %s", marker, e$message))
      })
    }
  }
  
  dev.off()
}

#' @title Custom PCA Biplot with Multi-Pattern Highlighting
#' @description Plots Individuals with flexible highlighting based on metadata substrings.
#' @param pca_res Result from run_coda_pca().
#' @param metadata Dataframe containing 'Patient_ID', 'Group' and the subgroup column.
#' @param colors Named vector of colors for main Groups (fill).
#' @param dims Integer vector of length 2 indicating which PCs to plot.
#' @param show_labels Logical. If TRUE, adds Patient_ID labels.
#' @param highlight_patterns Named vector (Pattern -> Color). E.g., c("_LS" = "#FFD700").
#' @param find_col_keyword Keyword to identify the metadata column to search.
#' @return A ggplot object.
plot_pca_custom <- function(pca_res, metadata, colors, dims = c(1, 2), 
                            show_labels = FALSE, 
                            highlight_patterns = NULL,
                            find_col_keyword = "Original|Source") {
  
  # 1. Setup Coordinates & Metadata
  if (length(dims) != 2) stop("dims parameter must be a vector of length 2")
  ind_coords <- as.data.frame(pca_res$ind$coord[, dims])
  colnames(ind_coords) <- c("X_Coord", "Y_Coord")
  
  if (nrow(ind_coords) != nrow(metadata)) stop("Metadata/PCA dimension mismatch.")
  plot_data <- cbind(metadata, ind_coords)
  
  # 2. Logic for Highlights
  # Default state: No Highlight
  plot_data$Highlight_Type <- "Standard"
  plot_data$Stroke_Size <- 0.5 
  
  # Identify the column to search
  target_col <- names(metadata)[grepl(find_col_keyword, names(metadata), ignore.case = TRUE)][1]
  
  has_highlights <- !is.null(highlight_patterns) && length(highlight_patterns) > 0 && !is.na(target_col)
  
  if (has_highlights) {
    target_values <- as.character(plot_data[[target_col]])
    
    for (pattern in names(highlight_patterns)) {
      matches <- grepl(pattern, target_values)
      # Only overwrite if currently Standard (prevents overwriting higher priority matches)
      to_update <- matches & (plot_data$Highlight_Type == "Standard")
      
      if (any(to_update)) {
        plot_data$Highlight_Type[to_update] <- pattern 
        plot_data$Stroke_Size[to_update] <- 2.0  
      }
    }
  }
  
  # Prepare Color Scale for Borders
  border_colors <- c("Standard" = "white")
  if (has_highlights) {
    border_colors <- c(highlight_patterns, border_colors)
  }
  
  # 3. Variance Explained
  eig_val <- get_eigenvalue(pca_res)
  var_x <- round(eig_val[dims[1], 2], 1)
  var_y <- round(eig_val[dims[2], 2], 1)
  
  # 4. Plotting
  p <- ggplot(plot_data, aes(x = X_Coord, y = Y_Coord, fill = Group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray80") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray80")
  
  # Identify groups where ellipse calculation is mathematically possible
  valid_groups <- c()
  unique_grps <- unique(as.character(plot_data$Group))
  
  for (g in unique_grps) {
    sub_dat <- plot_data[plot_data$Group == g, c("X_Coord", "Y_Coord")]
    
    # Attempt covariance calculation as proxy for stat_ellipse feasibility
    is_computable <- tryCatch({
      if (nrow(sub_dat) < 3) stop("Insufficient points")
      # cov() will fail or return NA/Inf if variance is problematic
      cov_mat <- cov(sub_dat)
      if (any(is.na(cov_mat))) stop("NA in covariance")
      TRUE
    }, error = function(e) { FALSE })
    
    if (is_computable) valid_groups <- c(valid_groups, g)
  }
  
  # Only add the layer if we have valid groups
  if (length(valid_groups) > 0) {
    p <- p + stat_ellipse(
      data = plot_data[plot_data$Group %in% valid_groups, ],
      geom = "polygon", 
      alpha = 0.1, 
      show.legend = FALSE, 
      level = 0.95
    )
  }
  
  p <- p +
    # Points 
    geom_point(aes(color = Highlight_Type, stroke = Stroke_Size), 
               size = 3.5, shape = 21) +
    
    # Main Fill Colors (Clinical Group)
    scale_fill_manual(values = colors, guide = guide_legend(override.aes = list(shape = 21))) +
    
    # Border Colors (Highlights)
    scale_color_manual(name = "Subgroup", values = border_colors) +
    scale_continuous_identity(aesthetics = "stroke") +
    
    labs(
      x = sprintf("PC%d (%s%%)", dims[1], var_x),
      y = sprintf("PC%d (%s%%)", dims[2], var_y),
      title = sprintf("PCA (PC%d vs PC%d)", dims[1], dims[2]),
      subtitle = if(has_highlights) "Colored borders indicate specific subgroups" else NULL
    ) +
    theme_coda()
  
  # Hide legend if no highlights exist
  if (!has_highlights) {
    p <- p + guides(color = "none")
  }
  
  if (show_labels) {
    p <- p + geom_text_repel(aes(label = Patient_ID), size = 3, show.legend = FALSE, max.overlaps = 40)
  }
  
  return(p)
}


#' @title Plot PCA Variance Dashboard
#' @description 
#' Visualizes explained variance per component with bars, curve, labels, 
#' and a side-list of cumulative variance values.
#' All labels are strictly in English.
#' 
#' @param pca_res Result object from FactoMineR::PCA.
#' @param n_list Integer. Number of top components to list on the side (default 10).
#' @return A ggplot object.
plot_pca_variance_dashboard <- function(pca_res, n_list = 10) {
  
  # 1. Extract Eigenvalues/Variance
  eig_df <- as.data.frame(pca_res$eig)
  # Standardize column names from FactoMineR
  colnames(eig_df) <- c("eigenvalue", "variance_percent", "cumulative_variance_percent")
  
  # Create PC identifiers (PC1, PC2...) and ensure factor order
  eig_df$PC <- factor(rownames(eig_df), levels = rownames(eig_df))
  eig_df$PC_Num <- 1:nrow(eig_df)
  
  # Limit to dimensions present
  n_pcs <- nrow(eig_df)
  # Dynamic limit for the plot (show max 15 bars for readability, or all if low dim)
  n_plot <- min(n_pcs, 15) 
  plot_data <- eig_df[1:n_plot, ]
  
  # 2. Prepare Side List Text
  n_list_actual <- min(n_pcs, n_list)
  list_data <- eig_df[1:n_list_actual, ]
  
  # Formatting the text table
  table_text <- paste0(
    "PC", list_data$PC_Num, ": ", 
    sprintf("%.1f%%", list_data$variance_percent), 
    " (Cum: ", sprintf("%.1f%%", list_data$cumulative_variance_percent), ")"
  )
  final_label <- paste(table_text, collapse = "\n")
  header_label <- paste0("Top ", n_list_actual, " Components\n(Var / Cumulative)")
  
  # 3. Plotting
  # We use a secondary axis scaling factor if needed, but for simplicity
  # since both variance and cumulative are %, we plot variance on Y.
  # The "Curve" usually represents the scree (variance), not cumulative, 
  # to match the histogram profile.
  
  p <- ggplot(plot_data, aes(x = PC, y = variance_percent)) +
    
    # A. Histograms (Bars)
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7, width = 0.7) +
    
    # B. The Curve (Scree line connecting bars)
    geom_line(aes(group = 1), color = "darkred", linewidth = 1, linetype = "dashed") +
    geom_point(color = "darkred", size = 2) + # Point size remains 'size'
    
    # C. Labels above histograms
    geom_text(aes(label = sprintf("%.1f%%", variance_percent)), 
              vjust = -0.5, size = 3.5, fontface = "bold") +
    
    # D. The Side List (Annotation)
    # We create a text annotation on the right. 
    # Logic: Place it at x = n_plot + 0.5, y = max_variance.
    annotate("text", x = n_plot + 0.6, y = max(plot_data$variance_percent) * 0.9, 
             label = header_label, hjust = 0, vjust = 1, fontface = "bold", size = 4, color = "gray20") +
    annotate("text", x = n_plot + 0.6, y = max(plot_data$variance_percent) * 0.8, 
             label = final_label, hjust = 0, vjust = 1, size = 3.5, color = "gray30", lineheight = 1.2) +
    
    # Scales & Expansion
    scale_y_continuous(limits = c(0, max(plot_data$variance_percent) * 1.15), 
                       expand = c(0, 0)) +
    # Expand X axis to make room for the list on the right
    scale_x_discrete(expand = expansion(add = c(0.6, 4))) + 
    
    # Theme & Labels
    labs(
      title = "PCA Scree Plot & Variance Explained",
      subtitle = "Explained Variance per Principal Component",
      x = "Principal Component",
      y = "Explained Variance (%)"
    ) +
    theme_coda() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  return(p)
}

#' @title Plot Stratification Heatmap (ComplexHeatmap)
#' @description 
#' Generates a clustered heatmap using Z-scored data.
#' @param mat_z A Z-scored numeric matrix (Samples x Markers).
#' @param metadata Dataframe containing patient metadata (aligned with mat_z rows).
#' @param annotation_colors_list A named list of color vectors for annotations.
#' @param title Plot title.
#' @return A ComplexHeatmap object (drawn).
plot_stratification_heatmap <- function(mat_z, metadata, annotation_colors_list, title = "Stratification") {
  
  # 1. Prepare Data
  mat_plot <- t(mat_z)
  
  if (!all(rownames(mat_z) == metadata$Patient_ID)) {
    stop("Mismatch between matrix rownames and metadata Patient_ID")
  }
  
  # Ensure metadata is a pure data.frame (not tibble) for HeatmapAnnotation
  metadata_clean <- as.data.frame(metadata)
  
  # 2. Setup Annotations
  ha <- HeatmapAnnotation(
    df = metadata_clean,
    col = annotation_colors_list,
    simple_anno_size = unit(0.3, "cm"),
    show_annotation_name = TRUE
  )
  
  # 3. Define Color Map for Z-Scores
  col_fun <- colorRamp2(c(-2, 0, 2), VIZ_DIVERGING)
  
  # 4. Draw Heatmap
  hm <- Heatmap(
    mat_plot,
    name = "Z-Score",
    col = col_fun,
    
    # Clustering
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    
    # Annotations
    top_annotation = ha,
    
    # Aesthetics
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 6),
    show_column_names = FALSE, 
    row_dend_width = unit(2, "cm"),
    column_dend_height = unit(2, "cm"),
    
    column_title = paste0(title, " (n=", ncol(mat_plot), ")"),
    row_title = "Markers"
  )
  
  return(hm)
}


#' @title Plot Signature Boxplots (Top Drivers)
#' @description 
#' Helper function to visualize the distribution of the top discriminating markers.
#' Draws boxplots of Z-scores for the features with highest loadings.
#' 
#' @param data_matrix Numeric matrix of data used in PLS (usually Z-scored).
#' @param group_factor Factor vector defining groups.
#' @param loadings_df Dataframe of loadings from extract_plsda_loadings.
#' @param comp Integer. Which component to visualize?
#' @param n_top Integer. How many top markers to plot?
#' @param colors Named vector of colors.
#' @return A ggplot object.
viz_plot_signature_boxplots <- function(data_matrix, group_factor, loadings_df, 
                                        comp = 1, n_top = 6, colors) {
  
  col_name <- paste0("Comp", comp, "_Weight")
  if (!col_name %in% names(loadings_df)) return(NULL)
  
  # Filter top absolute weights
  top_mks <- loadings_df %>%
    arrange(desc(abs(!!sym(col_name)))) %>%
    head(n_top) %>%
    pull(Marker)
  
  if (length(top_mks) == 0) return(NULL)
  
  # Prepare long dataframe for plotting
  plot_df <- as.data.frame(data_matrix[, top_mks, drop=FALSE])
  plot_df$Group <- group_factor
  
  long_df <- plot_df %>%
    pivot_longer(-Group, names_to = "Marker", values_to = "Z_Score") %>%
    mutate(Marker = factor(Marker, levels = top_mks)) # Preserve order
  
  p <- ggplot(long_df, aes(x = Group, y = Z_Score, fill = Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    facet_wrap(~Marker, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = colors) +
    labs(
      title = sprintf("Top %d Drivers Distribution (Comp %d)", n_top, comp),
      subtitle = "Visual validation of PLS-DA weights (Z-Scores)",
      y = "Standardized Abundance (Z)", x = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", strip.background = element_rect(fill="gray95"))
  
  return(p)
}

#' @title Plot sPLS-DA Signature Heatmap
#' @description 
#' Custom heatmap for sPLS-DA markers with patients on rows and markers on columns.
#' @param mat_z Numeric matrix of Z-scored data (Patients x Markers).
#' @param metadata Dataframe with patient metadata.
#' @param group_col Column name defining the group.
#' @param colors Named vector of group colors.
#' @param title Plot title.
#' @return A ComplexHeatmap object.
viz_plot_splsda_heatmap <- function(mat_z, metadata, group_col, colors, title = "Global Signature Heatmap") {
  requireNamespace("ComplexHeatmap", quietly = TRUE)
  requireNamespace("circlize", quietly = TRUE)
  requireNamespace("grid", quietly = TRUE)
  
  if (!all(rownames(mat_z) == metadata$Patient_ID)) {
    stop("Mismatch between matrix rownames and metadata Patient_ID")
  }
  
  # Setup Row Annotations (Patients)
  meta_df <- as.data.frame(metadata[, group_col, drop = FALSE])
  anno_colors <- list()
  anno_colors[[group_col]] <- colors
  
  ra <- ComplexHeatmap::rowAnnotation(
    df = meta_df,
    col = anno_colors,
    show_annotation_name = FALSE,
    simple_anno_size = grid::unit(0.4, "cm")
  )
  
  # Standardized Z-Score color mapping
  col_fun <- circlize::colorRamp2(c(-3, 0, 3), VIZ_DIVERGING)
  
  hm <- ComplexHeatmap::Heatmap(
    mat_z, 
    name = "Z-Score",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "ward.D2",
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "ward.D2",
    left_annotation = ra,
    row_names_gp = grid::gpar(fontsize = 7),
    column_names_gp = grid::gpar(fontsize = 8),
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_title = title,
    row_title = "Patients",
    column_title_side = "top"
  )
  
  return(hm)
}

#' @title Report sPLS-DA Visualization 
#' @description 
#' Generates a comprehensive PDF report for sPLS-DA results.
#' Includes robust error handling for PDF device closing and ggplot fallback.
#' 
#' @param pls_res The result object from run_splsda_model.
#' @param drivers_df Dataframe of extracted loadings.
#' @param metadata_viz Dataframe for plotting.
#' @param colors_viz Named vector of colors for groups.
#' @param out_path Path to save the PDF.
#' @param group_col The name of the metadata column to use for grouping/coloring (default: "Group").
#' @param n_top_boxplots Number of boxplots to plot. 
viz_report_plsda <- function(pls_res, drivers_df, metadata_viz, colors_viz, out_path, group_col = "Group", n_top_boxplots = 9) {
  
  if (is.null(pls_res$model)) return(NULL)
  
  requireNamespace("mixOmics", quietly = TRUE)
  requireNamespace("ComplexHeatmap", quietly = TRUE)
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(ggrepel)
  
  message(sprintf("   [Viz] Generating sPLS-DA graphical report: %s", basename(out_path)))
  
  if (!group_col %in% names(metadata_viz)) {
    stop(sprintf("Column '%s' not found in metadata for sPLS-DA visualization.", group_col))
  }
  
  pdf(out_path, width = 11, height = 8)
  on.exit(try(dev.off(), silent = TRUE), add = TRUE) 
  
  plot_group_factor <- as.factor(metadata_viz[[group_col]])
  levels_present <- levels(plot_group_factor)
  plot_colors <- colors_viz[levels_present]
  
  if (any(is.na(plot_colors))) {
    missing_grps <- levels_present[is.na(plot_colors)]
    warning(paste("Missing colors for:", paste(missing_grps, collapse=", ")))
    plot_colors[is.na(plot_colors)] <- "gray50"
    names(plot_colors) <- levels_present
  }
  
  n_comps <- pls_res$model$ncomp
  
  # Helper for clean ggplot error messages
  plot_error_msg <- function(msg) {
    ggplot() + 
      annotate("text", x = 0, y = 0, label = msg, color = "darkred", size = 5, fontface = "bold") + 
      theme_void()
  }
  
  # 1. Native ggplot2 Biplot
  tryCatch({
    variates_df <- as.data.frame(pls_res$model$variates$X[, 1:2, drop = FALSE])
    colnames(variates_df) <- c("PC1", "PC2")
    variates_df$Group <- plot_group_factor
    variates_df$Patient_ID <- rownames(pls_res$model$variates$X)
    
    var_pct <- round(pls_res$model$prop_expl_var$X[1:2] * 100, 1)
    
    p_biplot <- ggplot(variates_df, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
      stat_ellipse(geom = "polygon", alpha = 0.1, level = 0.95, show.legend = FALSE) +
      geom_point(size = 3.5, shape = 21, color = "white", stroke = 0.5) +
      geom_text_repel(aes(label = Patient_ID), size = 3, show.legend = FALSE, max.overlaps = 20) +
      scale_fill_manual(values = plot_colors) +
      scale_color_manual(values = plot_colors) +
      labs(
        title = "sPLS-DA: Model Separation",
        subtitle = paste("Grouping by:", group_col),
        x = sprintf("X-variate 1: %s%% expl. var", var_pct[1]),
        y = sprintf("X-variate 2: %s%% expl. var", var_pct[2])
      ) +
      theme_coda()
    
    print(p_biplot)
  }, error = function(e) {
    print(plot_error_msg(paste("Biplot Error:", e$message)))
  })
  
  # 2. Loadings & Boxplots
  for (i in 1:n_comps) {
    comp_col <- paste0("Comp", i, "_Weight")
    if (!comp_col %in% names(drivers_df)) next
    
    variates <- pls_res$model$variates$X[, i]
    stat_groups <- plot_group_factor
    group_means <- tapply(variates, stat_groups, mean)
    
    pos_group_name <- names(group_means)[which.max(group_means)]
    neg_group_name <- names(group_means)[which.min(group_means)]
    
    label_sub <- sprintf("Direction: Positive -> %s | Negative -> %s", pos_group_name, neg_group_name)
    
    current_fill_colors <- c()
    if(pos_group_name %in% names(plot_colors)) current_fill_colors[[pos_group_name]] <- plot_colors[[pos_group_name]]
    if(neg_group_name %in% names(plot_colors)) current_fill_colors[[neg_group_name]] <- plot_colors[[neg_group_name]]
    if(is.null(current_fill_colors[[pos_group_name]])) current_fill_colors[[pos_group_name]] <- "gray"
    if(is.null(current_fill_colors[[neg_group_name]])) current_fill_colors[[neg_group_name]] <- "gray"
    
    df_comp <- drivers_df[abs(drivers_df[[comp_col]]) > 0, ]
    
    if (nrow(df_comp) > 0) {
      df_comp$Association <- ifelse(df_comp[[comp_col]] > 0, pos_group_name, neg_group_name)
      
      p_load <- ggplot(df_comp, aes(x = reorder(Marker, abs(.data[[comp_col]])), 
                                    y = .data[[comp_col]], fill = Association)) +
        geom_bar(stat = "identity", width = 0.7) +
        coord_flip() +
        scale_fill_manual(values = current_fill_colors) +
        labs(title = sprintf("sPLS-DA Loadings (Component %d)", i), 
             subtitle = label_sub, x = "Marker", y = "Weight Contribution", fill = "Associated Group") +
        theme_coda() + theme(legend.position = "bottom")
      print(p_load)
      
      if(exists("viz_plot_signature_boxplots")) {
        p_box <- viz_plot_signature_boxplots(
          data_matrix = pls_res$model$X, 
          group_factor = plot_group_factor,
          loadings_df = df_comp, 
          comp = i, 
          n_top = n_top_boxplots, 
          colors = plot_colors
        )
        if (!is.null(p_box)) print(p_box)
      }
    }
  }
  
  # 3. Transposed Signature Heatmap
  if (nrow(drivers_df) > 1) {
    tryCatch({
      mks <- unique(drivers_df$Marker)
      mat_sub <- pls_res$model$X[, mks, drop = FALSE]
      
      meta_hm <- metadata_viz[match(rownames(mat_sub), metadata_viz$Patient_ID), ]
      
      hm <- viz_plot_splsda_heatmap(
        mat_z = mat_sub,
        metadata = meta_hm,
        group_col = group_col,
        colors = plot_colors,
        title = "sPLS-DA Signature Heatmap"
      )
      
      ComplexHeatmap::draw(hm, merge_legend = TRUE)
    }, error = function(e) {
      print(plot_error_msg(paste("Heatmap Error:", e$message)))
    })
  }
}


#' @title Extract Clinical Colors
#' @description Safely extracts colors for responder and non-responder labels from the global configuration.
#' Provides robust hex color fallbacks if configuration is missing.
#' @param config Global configuration list object.
#' @return A named character vector of hex color codes.
#' @export
get_clinical_colors <- function(config) {
  resp_lbl <- config$clinical$responder_label
  nresp_lbl <- config$clinical$non_responder_label
  
  colors_viz <- c()
  
  if (!is.null(config$colors$groups[[resp_lbl]])) {
    colors_viz[resp_lbl] <- config$colors$groups[[resp_lbl]]
  } else {
    colors_viz[resp_lbl] <- "#2E8B57" # Default SeaGreen
    warning(sprintf("[Viz] Missing color configuration for '%s'. Applying default green.", resp_lbl))
  }
  
  if (!is.null(config$colors$groups[[nresp_lbl]])) {
    colors_viz[nresp_lbl] <- config$colors$groups[[nresp_lbl]]
  } else {
    colors_viz[nresp_lbl] <- VIZ_RED # Default Firebrick
    warning(sprintf("[Viz] Missing color configuration for '%s'. Applying default red.", nresp_lbl))
  }

  return(colors_viz)
}

# ==============================================================================
# NETWORK VISUALIZATION FUNCTIONS
# ==============================================================================

#' @title Plot Network Structure (ggraph)
#' @description Force-directed network layout colored by edge sign.
#' @param adj_mat Adjacency matrix (0/1).
#' @param weight_mat Partial correlation matrix.
#' @param title Plot title.
#' @param layout_type ggraph layout algorithm.
#' @param min_cor Minimum absolute weight to render edge.
#' @return ggplot object.
plot_network_structure <- function(adj_mat, weight_mat, title = "Network",
                                   layout_type = "nicely", min_cor = 0) {

  requireNamespace("igraph", quietly = TRUE)
  requireNamespace("tidygraph", quietly = TRUE)
  requireNamespace("ggraph", quietly = TRUE)

  g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE)

  igraph::E(g)$weight_raw <- NA
  igraph::E(g)$sign <- NA
  igraph::E(g)$weight <- NA

  el <- igraph::as_data_frame(g, what = "edges")

  if (nrow(el) > 0) {
    weights <- numeric(nrow(el))
    signs   <- character(nrow(el))
    keep_edge <- logical(nrow(el))

    for (k in 1:nrow(el)) {
      w <- weight_mat[el[k, 1], el[k, 2]]
      if (abs(w) >= min_cor) {
        keep_edge[k] <- TRUE
        weights[k] <- abs(w)
        signs[k]   <- ifelse(w > 0, "Positive", "Negative")
      }
    }

    g <- igraph::delete_edges(g, igraph::E(g)[!keep_edge])

    if (igraph::ecount(g) > 0) {
      igraph::E(g)$weight <- weights[keep_edge]
      igraph::E(g)$sign   <- signs[keep_edge]
    }
  }

  igraph::V(g)$degree <- igraph::degree(g)
  tg <- tidygraph::as_tbl_graph(g)

  p <- ggraph::ggraph(tg, layout = layout_type) +
    ggraph::geom_edge_link(ggplot2::aes(width = weight, color = sign), alpha = 0.6) +
    ggraph::scale_edge_width(range = c(0.2, 1.5), guide = "none") +
    ggraph::scale_edge_color_manual(values = c("Positive" = "#00BFC4", "Negative" = "#F8766D")) +
    ggraph::geom_node_point(ggplot2::aes(size = degree), color = "gray20", fill = "white", shape = 21) +
    ggraph::geom_node_text(ggplot2::aes(label = name), repel = TRUE, size = 3, max.overlaps = 20) +
    ggplot2::scale_size_continuous(range = c(2, 8)) +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("Nodes: %d | Edges: %d (Filter: |rho| > %.2f)",
                         igraph::vcount(g), igraph::ecount(g), min_cor),
      edge_color = "Association"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "gray50"),
      legend.position = "bottom"
    )

  return(p)
}

#' @title Plot Partial Correlation Density
#' @description Density of edge weights with threshold annotation.
#' @param pcor_mat Numeric partial correlation matrix.
#' @param adj_mat Optional adjacency matrix (0/1) to report stable edge count.
#' @param threshold Magnitude threshold used for filtering.
#' @param group_label Label for the group.
#' @return ggplot object.
viz_plot_edge_density <- function(pcor_mat, adj_mat = NULL, threshold = 0.15, group_label = "") {

  vals <- pcor_mat[upper.tri(pcor_mat)]
  n_total <- length(vals)
  df_plot <- data.frame(Value = vals)

  if (!is.null(adj_mat)) {
    n_stable <- sum(adj_mat[upper.tri(adj_mat)])
    sub_text <- sprintf("Threshold: |rho| > %.4f | Final Stable Edges: %d", threshold, n_stable)
  } else {
    n_kept <- sum(abs(vals) >= threshold)
    pct_kept <- round((n_kept / n_total) * 100, 1)
    sub_text <- sprintf("Threshold: |rho| > %.4f (Keeps %.1f%% of raw edges)", threshold, pct_kept)
  }

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Value)) +
    ggplot2::geom_density(fill = "steelblue", alpha = 0.3)

  if (threshold > 0) {
    p <- p + ggplot2::geom_vline(xintercept = c(-threshold, threshold),
                                 linetype = "dashed", color = "red", linewidth = 0.8)
  }

  p <- p +
    ggplot2::labs(
      title = paste("Edge Weight Distribution:", group_label),
      subtitle = sub_text,
      x = "Partial Correlation (Shrinkage)",
      y = "Density"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 11, hjust = 0.5, color = "gray40")
    ) +
    ggplot2::xlim(-1, 1)

  return(p)
}

#' @title Plot Hub-Driver Quadrant
#' @description sPLS-DA Importance vs Network Topology scatter with role classification.
#' @param hub_driver_df Output from integrate_hub_drivers().
#' @param y_label Label for Y axis ("Degree" or "Betweenness").
#' @param title_suffix String to append to plot title.
#' @return ggplot object.
plot_hub_driver_quadrant <- function(hub_driver_df, y_label = "Degree", title_suffix = "") {

  if (is.null(hub_driver_df) || nrow(hub_driver_df) == 0) return(NULL)

  y_col_name <- if ("Topology_Metric_Value" %in% names(hub_driver_df)) "Topology_Metric_Value" else "Degree"

  x_mid <- median(hub_driver_df$Importance, na.rm = TRUE)
  y_mid <- median(hub_driver_df[[y_col_name]], na.rm = TRUE)

  p <- ggplot2::ggplot(hub_driver_df, ggplot2::aes(x = Importance, y = .data[[y_col_name]], fill = Role)) +
    ggplot2::geom_vline(xintercept = x_mid, linetype = "dashed", color = "gray60") +
    ggplot2::geom_hline(yintercept = y_mid, linetype = "dashed", color = "gray60") +
    ggplot2::geom_point(size = 4, shape = 21, alpha = 0.8) +
    ggrepel::geom_text_repel(ggplot2::aes(label = Marker), size = 3.5, max.overlaps = 20) +
    ggplot2::scale_fill_manual(values = c(
      "Master_Regulator"    = VIZ_RED,
      "Solo_Driver"         = VIZ_RED_LT,
      "Structural_Connector" = VIZ_BLUE,
      "Background"          = "gray80"
    )) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks()) +
    ggplot2::labs(
      title    = paste("Hub-Driver Analysis", title_suffix),
      subtitle = paste("sPLS-DA Importance vs", y_label),
      x        = "Statistical Importance (sPLS-DA)",
      y        = paste0("Topological Centrality (", y_label, ")"),
      caption  = "Quadrants defined by median values"
    ) +
    theme_coda() +
    ggplot2::theme(legend.position = "bottom")

  return(p)
}



#' @title Clinical-Utility Figure (Calibration + Decision Curve)
#' @description Paper-ready two-panel figure for the output of
#'   \code{run_clinical_utility()}. Panel A: leakage-free (LOO) calibration with
#'   quartile bins and the logistic-recalibrated curve overlaid (shows the
#'   shrinkage). Panel B: decision curve (net benefit vs treat-all / treat-none)
#'   with the model-dominance range shaded. Styled with \code{theme_coda()}.
#' @param clin_util  List returned by \code{run_clinical_utility()}.
#' @param colors     Named vector of group colours (from \code{get_clinical_colors}).
#' @param title_prefix Optional string prepended to the figure title (e.g. cohort).
#' @return A patchwork object (A | B). ASCII-only titles (PDF font safe).
viz_plot_clinical_utility <- function(clin_util, colors, title_prefix = "") {
  if (is.null(clin_util)) return(NULL)
  pos    <- clin_util$positive_label
  col_mod <- if (!is.null(colors[[pos]])) colors[[pos]] else VIZ_BLUE
  pp     <- clin_util$per_patient
  yb     <- as.integer(pp$True_Group == pos)
  cl     <- clin_util$calibration; dc <- clin_util$decision_curve

  # ── Panel A: calibration (LOO raw vs recalibrated) ──────────────────────────
  calib_long <- rbind(
    data.frame(pred = pp$Prob_LOO,   obs = yb, kind = "LOO (raw)"),
    data.frame(pred = pp$Prob_Recal, obs = yb, kind = "Recalibrated")
  )
  bins <- do.call(rbind, lapply(split(seq_along(yb), dplyr::ntile(pp$Prob_LOO, 4)),
    function(ix) {
      o <- mean(yb[ix])
      data.frame(p_mean = mean(pp$Prob_LOO[ix]), obs = o,
                 se = sqrt(o * (1 - o) / length(ix)))
    }))
  pA <- ggplot(calib_long, aes(pred, obs, colour = kind, linetype = kind)) +
    geom_abline(slope = 1, intercept = 0, linetype = 3, colour = "grey50") +
    geom_smooth(method = "loess", se = FALSE, span = 1, linewidth = 0.9) +
    geom_point(data = bins, aes(p_mean, obs), inherit.aes = FALSE, size = 2.4) +
    geom_errorbar(data = bins,
                  aes(x = p_mean, ymin = pmax(0, obs - se), ymax = pmin(1, obs + se)),
                  inherit.aes = FALSE, width = 0.012) +
    scale_colour_manual(values = c("LOO (raw)" = col_mod, "Recalibrated" = "grey40")) +
    scale_linetype_manual(values = c("LOO (raw)" = 1, "Recalibrated" = 2)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = "Calibration (LOO)",
         subtitle = sprintf("CITL=%.2f  slope=%.2f  Brier=%.3f  (recal. slope=%.2f)",
                            cl$loo$intercept, cl$loo$slope, cl$loo$brier, cl$recalibrated$slope),
         x = "Predicted probability", y = "Observed fraction", colour = NULL, linetype = NULL) +
    theme_coda()

  # ── Panel B: decision curve ─────────────────────────────────────────────────
  cur <- dc$curve
  dom <- cur[cur$model > cur$treat_all & cur$model > cur$treat_none, , drop = FALSE]
  nb_long <- rbind(
    data.frame(threshold = cur$threshold, net_benefit = cur$model,      strategy = "Model"),
    data.frame(threshold = cur$threshold, net_benefit = cur$treat_all,  strategy = "Treat all"),
    data.frame(threshold = cur$threshold, net_benefit = cur$treat_none, strategy = "Treat none")
  )
  pB <- ggplot(nb_long, aes(threshold, net_benefit, colour = strategy))
  if (nrow(dom) > 0)
    pB <- pB + annotate("rect", xmin = min(dom$threshold), xmax = max(dom$threshold),
                        ymin = -Inf, ymax = Inf, alpha = 0.08, fill = col_mod)
  pB <- pB +
    geom_line(linewidth = 0.9) +
    scale_colour_manual(values = c("Model" = col_mod, "Treat all" = VIZ_RED,
                                   "Treat none" = "grey45")) +
    coord_cartesian(ylim = c(min(-0.02, min(cur$treat_all)), max(cur$model) + 0.02)) +
    labs(title = "Decision curve",
         subtitle = sprintf("Net benefit > defaults over pt [%s]; prevalence=%.2f",
                            dc$dominance_range, dc$prevalence),
         x = "Threshold probability", y = "Net benefit", colour = NULL) +
    theme_coda()

  ttl <- paste0(if (nzchar(title_prefix)) paste0(title_prefix, " - ") else "",
                "Clinical utility of pre-specified composite (",
                paste(clin_util$composite_markers, collapse = "+"), ")")
  sub <- sprintf("AUC corrected=%.3f [%.3f-%.3f] | optimism=%.4f | gate: %s",
                 clin_util$discrimination$auc_corrected,
                 clin_util$discrimination$auc_corrected_ci[1],
                 clin_util$discrimination$auc_corrected_ci[2],
                 clin_util$discrimination$auc_optimism, clin_util$gate_provenance)

  patchwork::wrap_plots(pA, pB, nrow = 1) +
    patchwork::plot_annotation(title = ttl, subtitle = sub, tag_levels = "A",
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5),
                    plot.subtitle = element_text(hjust = 0.5, colour = "gray40")))
}

#' @title Plot clinical + cytometric added-value (decision curve + LOO AUC)
#' @description Two panels: (A) decision curve comparing the standard-of-care
#'   clinical model vs clinical+immune (vs treat-all / treat-none), and (B) a bar
#'   of leakage-free LOO AUCs (clinical / immune / combined). Subtitle carries the
#'   increment (ΔAUC, DeLong p, IDI with 95% CI).
#' @param av Result of run_clinical_immune_added_value().
viz_plot_added_value <- function(av, colors, title_prefix = "") {
  if (is.null(av)) return(NULL)
  pos     <- av$positive_label
  col_comb <- if (!is.null(colors[[pos]])) colors[[pos]] else VIZ_BLUE
  cur <- av$decision_curve$curve; inc <- av$increment

  # ── Panel A: decision curve (clinical vs combined) ──────────────────────────
  nb_long <- rbind(
    data.frame(threshold = cur$threshold, net_benefit = cur$clinical,   strategy = "Clinical only"),
    data.frame(threshold = cur$threshold, net_benefit = cur$combined,   strategy = "Clinical + immune"),
    data.frame(threshold = cur$threshold, net_benefit = cur$treat_all,  strategy = "Treat all"),
    data.frame(threshold = cur$threshold, net_benefit = cur$treat_none, strategy = "Treat none")
  )
  pA <- ggplot(nb_long, aes(threshold, net_benefit, colour = strategy, linetype = strategy)) +
    geom_line(linewidth = 0.9) +
    scale_colour_manual(values = c("Clinical only" = "#E08214", "Clinical + immune" = col_comb,
                                   "Treat all" = VIZ_RED, "Treat none" = "grey45")) +
    scale_linetype_manual(values = c("Clinical only" = 2, "Clinical + immune" = 1,
                                     "Treat all" = 3, "Treat none" = 3)) +
    coord_cartesian(ylim = c(min(-0.02, min(cur$treat_all)),
                             max(c(cur$clinical, cur$combined)) + 0.02)) +
    labs(title = "Decision curve", x = "Threshold probability", y = "Net benefit",
         colour = NULL, linetype = NULL,
         subtitle = sprintf("prevalence=%.2f", av$decision_curve$prevalence)) +
    theme_coda()

  # ── Panel B: LOO AUC bars ───────────────────────────────────────────────────
  auc_df <- data.frame(
    model = factor(c("Clinical", "Immune", "Clinical+immune"),
                   levels = c("Clinical", "Immune", "Clinical+immune")),
    auc   = c(av$auc$clinical[["loo"]], av$auc$immune[["loo"]], av$auc$combined[["loo"]]))
  pB <- ggplot(auc_df, aes(model, auc, fill = model)) +
    geom_col(width = 0.6) +
    geom_hline(yintercept = 0.5, linetype = 3, colour = "grey50") +
    geom_text(aes(label = sprintf("%.3f", auc)), vjust = -0.4, size = 3.4) +
    scale_fill_manual(values = c("Clinical" = "#E08214", "Immune" = "#80CDC1",
                                 "Clinical+immune" = col_comb), guide = "none") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = "Discrimination (leakage-free LOO AUC)", x = NULL, y = "AUC") +
    theme_coda()

  ttl <- paste0(if (nzchar(title_prefix)) paste0(title_prefix, " - ") else "",
                "Added value of immune profiling over clinical (",
                paste(names(av$clinical_vars), collapse = "+"), ")")
  # LRT of the increment = the correct primary test for nested models (DeLong on
  # ΔAUC is expected to be ns); lead with it, guard if absent (older runs / no logistf).
  lrt_txt <- if (!is.null(inc$lrt_p)) {
    pp <- if (!is.null(inc$lrt_perm_p)) sprintf(" (perm %.3f)", inc$lrt_perm_p) else ""
    sprintf("LRT p=%.4f%s   ", inc$lrt_p, pp)
  } else ""
  sub <- sprintf("%sDelta-AUC(LOO)=%+.3f  DeLong p=%.3f  IDI=%+.3f [%.3f, %.3f]  | %s",
                 lrt_txt, inc$delta_auc_loo, inc$delong_p_loo, inc$idi,
                 inc$idi_ci[1], inc$idi_ci[2], av$gate_provenance)

  patchwork::wrap_plots(pA, pB, nrow = 1) +
    patchwork::plot_annotation(title = ttl, subtitle = sub, tag_levels = "A",
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5),
                    plot.subtitle = element_text(hjust = 0.5, colour = "gray40")))
}