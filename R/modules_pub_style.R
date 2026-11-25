# R/modules_pub_style.R
# ==============================================================================
# PUBLICATION STYLE MODULE  (single source of truth for manuscript figures)
# ------------------------------------------------------------------------------
# Figure conventions: cairo_pdf (fonts embedded), RGB, column
# widths 85 mm / 180 mm, >=8 pt text at final size, NO embedded titles/subtitles
# (stats live in the figure captions), bold A/B panel tags, Okabe-Ito colour-
# blind-safe palette (no red/green collision).
#
# These helpers are consumed by R/modules_pub_figures.R and are intentionally
# decoupled from the pipeline's internal-QC theme (theme_coda) so the QC figures
# and the publication figures can evolve independently.
# ==============================================================================

suppressMessages({ library(ggplot2); library(patchwork) })

# Canvas widths (mm). Typical print maxima are 85 (single) / 180 (double); we
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

#' A P / FDR value as the text writes it: "< 0.001" below 0.001, else 2 significant
#' digits (0.0235 -> 0.023, 0.238 -> 0.24, 0.9 -> 0.90).
pub_fmt_p <- function(p) {
  ifelse(p < 0.001, "< 0.001", formatC(signif(p, 2), digits = 2, format = "fg", flag = "#"))
}
#' Panel label for a P / FDR value: "<name> < 0.001" or "<name> = x" (pub_fmt_p()).
pub_lab_p <- function(p, name = "P") {
  ifelse(p < 0.001, sprintf("%s < 0.001", name), sprintf("%s = %s", name, pub_fmt_p(p)))
}

#' Save a publication figure at an exact mm size through save_figure() (formats and
#' dpi from config$figures; default embedded-font cairo PDF only).
#' @param path Output path; a format extension (e.g. ".pdf") is stripped to the stem.
#' @param config Pass config or its \code{figures} block (NULL = defaults).
#' @return Invisibly, the primary output file, or NULL when \code{plot} is NULL.
save_pub_figure <- function(plot, path, width_mm, height_mm, config = NULL) {
  if (is.null(plot)) return(invisible(NULL))
  files <- save_figure(plot, path, width_mm, height_mm, units = "mm", config = config)
  message(sprintf("   [PubFig] %s  (%.0f x %.0f mm)", basename(files[1]), width_mm, height_mm))
  invisible(files[1])
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

# Display-only marker names, matching the manuscript text (config/data codes unchanged).
# Same named-lookup-with-passthrough contract as pub_relabel_group(): unmapped markers
# keep their code. Clinical variables only need the underscore code artefact removed.
# The whole FACS panel is mapped (full-panel forest): each antigen of the code gets ⁺, the
# memory-subset codes (CM, EFF, EM, EMRA) stay as written, naïve is lowercase as in the text
# (antigen⁺ naïve, with a space), "tot" = the whole
# population. KI657EMRA is the database's spelling of the Ki67⁺EMRA column.
PUB_MARKER_LABELS <- c(
  "KI67NAIVE" = "Ki67⁺ naïve", "CD28KI67" = "CD28⁺Ki67⁺",
  "CD3tot" = "CD3⁺", "CD4tot" = "CD4⁺", "CD8tot" = "CD8⁺",
  "NAIVE" = "naïve", "CM" = "CM", "EFF" = "EFF", "EM" = "EM", "EMRA" = "EMRA",
  "CD28" = "CD28⁺", "CD28CD4" = "CD28⁺CD4⁺", "CD28CD8" = "CD28⁺CD8⁺", "CD28NAIVE" = "CD28⁺ naïve",
  "CD28CM" = "CD28⁺CM", "CD28EFF" = "CD28⁺EFF", "CD28EMRA" = "CD28⁺EMRA",
  "CD137tot" = "CD137⁺", "CD137CD4" = "CD137⁺CD4⁺", "CD137CD8" = "CD137⁺CD8⁺",
  "CD137KI67" = "CD137⁺Ki67⁺", "CD137PD1" = "CD137⁺PD1⁺", "CD137NAIVE" = "CD137⁺ naïve",
  "CD137CM" = "CD137⁺CM", "CD137EFF" = "CD137⁺EFF", "CD137EMRA" = "CD137⁺EMRA",
  "KI67" = "Ki67⁺", "KI67CD4" = "Ki67⁺CD4⁺", "KI67CD8" = "Ki67⁺CD8⁺",
  "KI67CM" = "Ki67⁺CM", "KI67EFF" = "Ki67⁺EFF", "KI657EMRA" = "Ki67⁺EMRA",
  "PD1" = "PD1⁺", "PD1CD4" = "PD1⁺CD4⁺", "PD1CD8" = "PD1⁺CD8⁺", "PD1KI67" = "PD1⁺Ki67⁺",
  "PD1NAIVE" = "PD1⁺ naïve", "PD1CM" = "PD1⁺CM", "PD1EFF" = "PD1⁺EFF", "PD1EMRA" = "PD1⁺EMRA")
pub_relabel_marker <- function(x) {
  x  <- as.character(x)
  hit <- x %in% names(PUB_MARKER_LABELS)
  x[hit] <- PUB_MARKER_LABELS[x[hit]]
  x
}
# PD-L1 is shown as TPS and the immune composite as "Immune score".
PUB_VAR_LABELS <- c("PD_L1" = "TPS", "PD-L1" = "TPS", "Immune composite" = "Immune score")
pub_relabel_var <- function(x) {
  x  <- as.character(x)
  hit <- x %in% names(PUB_VAR_LABELS)
  x[hit] <- PUB_VAR_LABELS[x[hit]]
  gsub("_", "-", x, fixed = TRUE)
}

#' ROC data.frame + AUC from a probability vector and a 0/1 positive indicator.
pub_roc_df <- function(prob, y_pos) {
  if (!requireNamespace("pROC", quietly = TRUE)) return(NULL)
  r <- pROC::roc(as.integer(y_pos), prob, direction = "<", quiet = TRUE)
  list(df = data.frame(FPR = 1 - r$specificities, TPR = r$sensitivities),
       auc = as.numeric(pROC::auc(r)))
}
