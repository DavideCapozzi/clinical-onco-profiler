# R/modules_pub_style.R
# ==============================================================================
# PUBLICATION STYLE MODULE  (single source of truth for manuscript figures)
# ------------------------------------------------------------------------------
# Frontiers in Immunology compliance: cairo_pdf (fonts embedded), RGB, column
# widths 85 mm / 180 mm, >=8 pt text at final size, NO embedded titles/subtitles
# (stats live in the figure captions), bold A/B panel tags, Okabe-Ito colour-
# blind-safe palette (no red/green collision).
#
# These helpers are consumed by R/modules_pub_figures.R and are intentionally
# decoupled from the pipeline's internal-QC theme (theme_coda) so the QC figures
# and the publication figures can evolve independently.
# ==============================================================================

suppressMessages({ library(ggplot2); library(patchwork) })

# Canvas widths (mm). Frontiers print maxima are 85 (single) / 180 (double); we
# render a touch larger for a comfortable, legible working draft — final
# typesetting rescales to the column grid without loss (vector PDF).
PUB_W1 <- 120   # ~1.4x single column
PUB_W2 <- 190   # ~full width

# Okabe-Ito colour-blind-safe palette, mapped to manuscript roles.
pub_palette <- c(
  clinical   = "#E69F00",  # orange
  immune     = "#009E73",  # bluish green
  combined   = "#0072B2",  # blue
  resp       = "#0072B2",  # positive class
  observed   = "#000000",
  ridge      = "#0072B2",
  unpen      = "#D55E00",  # vermillion
  treat_all  = "grey55",
  treat_none = "grey20",
  ref        = "grey60",
  ok         = "#0072B2",  # CI excludes the null
  fragile    = "grey55"    # CI includes the null
)

#' Publication ggplot theme (no title/subtitle; captions are external).
#' base_size chosen so axis text stays >= 8 pt at the final printed mm size.
theme_publication <- function(base = 11) {
  theme_bw(base_size = base, base_family = "sans") +
    theme(
      plot.title       = element_blank(),
      plot.subtitle    = element_blank(),
      axis.title       = element_text(size = base, face = "bold"),
      axis.text        = element_text(size = base - 1, colour = "black"),
      legend.title     = element_blank(),
      legend.text      = element_text(size = base - 1),
      legend.key.size  = unit(3.4, "mm"),
      legend.position  = "bottom",
      legend.margin    = margin(0, 0, 0, 0),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", colour = "grey80"),
      strip.text       = element_text(size = base - 1, face = "bold")
    )
}

#' Apply bold A/B... panel tags to a patchwork composition.
pub_tag <- function(p) p + patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 11, face = "bold"))

#' Save a publication figure as embedded-font vector PDF at an exact mm size.
save_pub_figure <- function(plot, path, width_mm, height_mm) {
  if (is.null(plot)) return(invisible(NULL))
  ggsave(path, plot, device = grDevices::cairo_pdf,
         width = width_mm, height = height_mm, units = "mm")
  message(sprintf("   [PubFig] %s  (%.0f x %.0f mm)", basename(path), width_mm, height_mm))
  invisible(path)
}

# ---- small numeric helpers ---------------------------------------------------
pub_invlogit <- function(x) 1 / (1 + exp(-x))
pub_clamp01  <- function(p, e = 1e-4) pmin(pmax(p, e), 1 - e)

# Display-only relabel of group codes for figures (config labels are unchanged).
# Named-lookup with passthrough — NEVER blanket gsub("_","/") (would corrupt other
# experiments, e.g. Tox_Yes). "RP" is the Italian "Risposta Parziale" → English RECIST
# "PR"; "SD_PD" underscore is a code artefact → "SD/PD". Extend the map as needed.
PUB_GROUP_LABELS <- c("RP" = "PR", "SD_PD" = "SD/PD")
pub_relabel_group <- function(x) {
  x  <- as.character(x)
  hit <- x %in% names(PUB_GROUP_LABELS)
  x[hit] <- PUB_GROUP_LABELS[x[hit]]
  x
}

#' ROC data.frame + AUC from a probability vector and a 0/1 positive indicator.
pub_roc_df <- function(prob, y_pos) {
  if (!requireNamespace("pROC", quietly = TRUE)) return(NULL)
  r <- pROC::roc(as.integer(y_pos), prob, direction = "<", quiet = TRUE)
  list(df = data.frame(FPR = 1 - r$specificities, TPR = r$sensitivities),
       auc = as.numeric(pROC::auc(r)))
}
