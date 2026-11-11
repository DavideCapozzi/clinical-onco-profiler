# R/utils_figures.R
# ==============================================================================
# FIGURE EXPORT — single helper for every figure the pipeline writes
# Description: save_figure() renders a ggplot / patchwork, a list of plots
#              (multipage) or a draw callback (ComplexHeatmap, base graphics)
#              to every format in config$figures$formats at config$figures$dpi.
#   pdf  = grDevices::cairo_pdf (embedded fonts, vector, multipage in one file)
#   png  = ragg::agg_png
#   tiff = ragg::agg_tiff(compression = "lzw")
#   svg  = grDevices::svg (cairo)
# Raster / svg exports of a multipage figure are written one file per page
# (<stem>_p01.png, <stem>_p02.png, ...); the pdf keeps all pages in one file.
# ==============================================================================

FIGURE_FORMATS <- c("pdf", "png", "tiff", "svg")

#' @title Resolve the Figure Export Settings
#' @param config A pass config (uses \code{config$figures}), a \code{figures}
#'   block itself, or NULL.
#' @return List with \code{formats} (character) and \code{dpi} (numeric).
figure_settings <- function(config = NULL) {
  fc <- if (!is.null(config$figures)) config$figures
        else if (!is.null(config$formats) || !is.null(config$dpi)) config
        else list()
  formats <- if (!is.null(fc$formats)) tolower(as.character(unlist(fc$formats))) else "pdf"
  bad <- setdiff(formats, FIGURE_FORMATS)
  if (length(bad))
    stop(sprintf("[FATAL] figures.formats: unsupported format(s) %s (allowed: %s)",
                 paste(bad, collapse = ", "), paste(FIGURE_FORMATS, collapse = ", ")))
  if (!length(formats)) formats <- "pdf"
  dpi <- if (!is.null(fc$dpi)) as.numeric(fc$dpi) else 300
  list(formats = unique(formats), dpi = dpi)
}

#' @title Save a Figure in Every Configured Format
#' @param plot_or_drawfn A ggplot/patchwork object, a list of them (one page
#'   each), or a function with no arguments that draws on the open device (it may
#'   draw several pages).
#' @param path_stem Output path WITHOUT extension (a trailing format extension is
#'   stripped).
#' @param width,height Figure size in \code{units}.
#' @param units "mm" (default) or "in".
#' @param config Pass config or \code{figures} block (see \code{figure_settings}).
#' @param multipage For a draw callback: TRUE when it draws several pages (raster
#'   and svg exports then get a page suffix). Ignored for plots and lists.
#' @return Invisibly, the character vector of files written (the first is the
#'   primary format's file), or NULL when there was nothing to draw.
save_figure <- function(plot_or_drawfn, path_stem, width, height, units = "mm",
                        config = NULL, multipage = NULL) {
  if (is.null(plot_or_drawfn)) return(invisible(NULL))
  units <- match.arg(units, c("mm", "in"))
  w_in <- if (units == "mm") width / 25.4 else width
  h_in <- if (units == "mm") height / 25.4 else height
  fs   <- figure_settings(config)
  path_stem <- sub(sprintf("\\.(%s)$", paste(FIGURE_FORMATS, collapse = "|")), "", path_stem)
  dir.create(dirname(path_stem), showWarnings = FALSE, recursive = TRUE)

  is_fn   <- is.function(plot_or_drawfn)
  is_list <- !is_fn && is.list(plot_or_drawfn) && !inherits(plot_or_drawfn, "gg")
  pages   <- if (is_list) Filter(Negate(is.null), plot_or_drawfn) else NULL
  if (is_list && !length(pages)) return(invisible(NULL))
  multi <- if (is_fn) isTRUE(multipage) else (is_list && length(pages) > 1)

  draw <- function() {
    if (is_fn) plot_or_drawfn()
    else if (is_list) for (p in pages) print(p)
    else print(plot_or_drawfn)
  }

  written <- character(0)
  for (fmt in fs$formats) {
    file <- if (fmt != "pdf" && multi) sprintf("%s_p%%02d.%s", path_stem, fmt)
            else sprintf("%s.%s", path_stem, fmt)
    switch(fmt,
      pdf  = grDevices::cairo_pdf(file, width = w_in, height = h_in, onefile = TRUE),
      png  = ragg::agg_png(file, width = w_in, height = h_in, units = "in", res = fs$dpi),
      tiff = ragg::agg_tiff(file, width = w_in, height = h_in, units = "in", res = fs$dpi,
                            compression = "lzw"),
      svg  = grDevices::svg(file, width = w_in, height = h_in, onefile = !multi))
    dev_id <- grDevices::dev.cur()
    tryCatch(draw(), finally = if (dev_id > 1 && dev_id %in% grDevices::dev.list())
                                 grDevices::dev.off(dev_id))
    written <- c(written, file)
  }
  invisible(written)
}

# ==============================================================================
# PUBLICATION EXPORT — raster / vector copies of the collected publication pdfs
# config figures.export: {formats: [png], dpi: 600, roles: [main, supplementary]}.
# collect_publication() converts every pdf of those roles into <role>/<format>/ with
# poppler's pdftocairo (env/environment.yml), so each copy is the vetted pdf itself
# (same layout and embedded fonts), not a second rendering. White background; tiff is
# LZW-compressed; a multipage pdf gives one file per page (<stem>_p01.<ext>, ...).
# ==============================================================================

EXPORT_FORMATS <- c(png = "png", tiff = "tif", jpeg = "jpg", svg = "svg", eps = "eps")
EXPORT_ROLES   <- c("main", "supplementary", "additional")

#' @title Resolve the Publication Export Settings
#' @param config Base config (uses \code{config$figures$export}).
#' @return NULL when no export is configured, else list(formats, dpi, roles).
figure_export_settings <- function(config) {
  ex <- config$figures$export
  if (is.null(ex) || !length(ex$formats)) return(NULL)
  formats <- unique(tolower(as.character(unlist(ex$formats))))
  bad <- setdiff(formats, names(EXPORT_FORMATS))
  if (length(bad))
    stop(sprintf("[FATAL] figures.export.formats: unsupported %s (allowed: %s)",
                 paste(bad, collapse = ", "), paste(names(EXPORT_FORMATS), collapse = ", ")))
  native <- intersect(formats, figure_settings(config)$formats)
  if (length(native))
    stop(sprintf("[FATAL] figures.export.formats: %s already rendered natively (figures.formats)",
                 paste(native, collapse = ", ")))
  dpi <- if (!is.null(ex$dpi)) as.numeric(ex$dpi) else 600
  if (length(dpi) != 1 || is.na(dpi) || dpi <= 0) stop("[FATAL] figures.export.dpi must be a positive number")
  roles <- if (!is.null(ex$roles)) as.character(unlist(ex$roles)) else c("main", "supplementary")
  bad <- setdiff(roles, EXPORT_ROLES)
  if (length(bad))
    stop(sprintf("[FATAL] figures.export.roles: unknown %s (allowed: %s)",
                 paste(bad, collapse = ", "), paste(EXPORT_ROLES, collapse = ", ")))
  list(formats = formats, dpi = dpi, roles = roles)
}

#' @title Pre-flight: the export converter is available when an export is configured
#' @return Invisibly TRUE; stops (FATAL) when pdftocairo / pdfinfo are missing.
check_figure_export <- function(config) {
  if (is.null(figure_export_settings(config))) return(invisible(TRUE))
  miss <- c("pdftocairo", "pdfinfo")[!nzchar(Sys.which(c("pdftocairo", "pdfinfo")))]
  if (length(miss))
    stop(sprintf("[FATAL] figures.export needs poppler (%s not found): micromamba install -c conda-forge poppler",
                 paste(miss, collapse = ", ")))
  invisible(TRUE)
}

#' @title Number of pages of a pdf (poppler pdfinfo)
pdf_page_count <- function(pdf) {
  out <- suppressWarnings(system2("pdfinfo", shQuote(pdf), stdout = TRUE, stderr = TRUE))
  n <- suppressWarnings(as.integer(sub("^Pages:\\s+", "", grep("^Pages:", out, value = TRUE))))
  if (!length(n) || is.na(n)) stop(sprintf("[FATAL] pdfinfo could not read %s", pdf))
  n
}

#' @title Convert one pdf into one export format
#' @param pdf Source pdf.
#' @param out_dir Output directory (created).
#' @param format One of \code{names(EXPORT_FORMATS)}.
#' @param dpi Resolution of the raster formats (svg / eps are vector: dpi ignored).
#' @return Character vector of the files written (one per page).
export_pdf <- function(pdf, out_dir, format, dpi = 600) {
  format <- match.arg(format, names(EXPORT_FORMATS))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  stem   <- file.path(out_dir, sub("\\.pdf$", "", basename(pdf), ignore.case = TRUE))
  n      <- pdf_page_count(pdf)
  raster <- format %in% c("png", "tiff", "jpeg")
  opts   <- switch(format,
    png  = c("-png", "-r", dpi),
    tiff = c("-tiff", "-tiffcompression", "lzw", "-r", dpi),
    jpeg = c("-jpeg", "-jpegopt", "quality=95", "-r", dpi),
    svg  = "-svg",
    eps  = "-eps")
  files <- character(0)
  for (p in seq_len(n)) {
    out <- if (n > 1) sprintf("%s_p%02d", stem, p) else stem
    target <- if (raster) out else paste0(out, ".", EXPORT_FORMATS[[format]])
    args <- c(opts, "-f", p, "-l", p, if (raster) "-singlefile", shQuote(pdf), shQuote(target))
    msg <- suppressWarnings(system2("pdftocairo", args, stdout = TRUE, stderr = TRUE))
    file <- paste0(out, ".", EXPORT_FORMATS[[format]])
    if (!is.null(attr(msg, "status")) || !file.exists(file))
      stop(sprintf("[FATAL] pdftocairo failed on %s (%s): %s", basename(pdf), format, paste(msg, collapse = " ")))
    files <- c(files, file)
  }
  files
}
