#!/usr/bin/env Rscript
# tests/test_figure_export.R
# ==============================================================================
# Unit test of the publication export (R/utils_figures.R export_pdf(),
# figure_export_settings(); R/utils_publication.R collect_publication()) on
# SYNTHETIC pdfs: every export format, page suffixes of a multipage pdf, raster size
# and dpi tag, the <role>/<format>/ layout of a collected run, INDEX rows, roles left
# out of the export, and the config validation.
#
# Run:  Rscript tests/test_figure_export.R      (from the repository root, ~5 s; needs poppler)
# Exit: 0 on pass, 1 on failure.
# ==============================================================================
suppressPackageStartupMessages({
  library(here); library(jsonlite)
  source(here("R/utils_io.R")); source(here("R/utils_figures.R")); source(here("R/utils_publication.R"))
})

fails <- character(0)
check <- function(ok, what) if (!isTRUE(ok)) fails <<- c(fails, what)

make_pdf <- function(path, pages = 1, w_mm = 100, h_mm = 50) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::cairo_pdf(path, width = w_mm / 25.4, height = h_mm / 25.4, onefile = TRUE)
  for (p in seq_len(pages)) { plot.new(); text(0.5, 0.5, sprintf("page %d", p)) }
  invisible(grDevices::dev.off())
  invisible(path)
}
png_info <- function(f) {   # width, height, pHYs pixels per metre (x)
  b <- readBin(f, "raw", file.info(f)$size); i <- 9; out <- list()
  while (i < length(b)) {
    n <- sum(as.integer(b[i:(i + 3)]) * 256^(3:0)); t <- rawToChar(b[(i + 4):(i + 7)])
    d <- if (n) b[(i + 8):(i + 7 + n)] else raw(0)
    if (t == "IHDR") { out$w <- sum(as.integer(d[1:4]) * 256^(3:0)); out$h <- sum(as.integer(d[5:8]) * 256^(3:0)) }
    if (t == "pHYs") out$ppm <- sum(as.integer(d[1:4]) * 256^(3:0))
    if (t == "IDAT") break
    i <- i + 12 + n
  }
  out
}

check_figure_export(list(figures = list(export = list(formats = "png"))))   # stops if poppler is missing
work <- tempfile("fig_export_"); dir.create(work)

# --- 1. export_pdf(): every format, single and multipage --------------------------
one <- make_pdf(file.path(work, "src", "One.pdf"))
two <- make_pdf(file.path(work, "src", "Two.pdf"), pages = 2)
for (fmt in names(EXPORT_FORMATS)) {
  ext <- EXPORT_FORMATS[[fmt]]
  f1 <- export_pdf(one, file.path(work, "out", fmt), fmt, dpi = 300)
  check(identical(basename(f1), sprintf("One.%s", ext)), sprintf("%s single-page name: %s", fmt, paste(basename(f1), collapse = ",")))
  f2 <- export_pdf(two, file.path(work, "out", fmt), fmt, dpi = 300)
  check(identical(basename(f2), sprintf("Two_p%02d.%s", 1:2, ext)), sprintf("%s multipage names: %s", fmt, paste(basename(f2), collapse = ",")))
  check(all(file.size(c(f1, f2)) > 0), sprintf("%s: empty file", fmt))
}
# cairo_pdf rounds the page to whole points (<= 1 pt = 300/72 px here); the pHYs tag is truncated
pi1 <- png_info(file.path(work, "out", "png", "One.png"))
tol <- 300 / 72 + 1
check(abs(pi1$w - 100 / 25.4 * 300) <= tol && abs(pi1$h - 50 / 25.4 * 300) <= tol,
      sprintf("png size %d x %d at 300 dpi", pi1$w, pi1$h))
check(abs(pi1$ppm - 300 / 0.0254) <= 1, sprintf("png dpi tag %s", pi1$ppm))

# --- 2. collect_publication(): <role>/<format>/ layout and INDEX ------------------
run <- file.path(work, "run"); exp <- "Synthetic"
make_pdf(file.path(run, exp, "06_machine_learning", "publication", sprintf("Figure1_Biology_%s.pdf", exp)))
make_pdf(file.path(run, exp, "07_post_hoc", "therapy", "FigureT_TherapyChange.pdf"), pages = 2)
make_pdf(file.path(run, exp, "06_machine_learning", "publication", sprintf("FigureS1_CONSORT_%s.pdf", exp)))
cfg <- list(
  figures = list(formats = "pdf", dpi = 300, export = list(formats = c("png", "tiff"), dpi = 150,
                                                           roles = c("main", "supplementary"))),
  publication = list(
    main          = list(list(stem = "Figure1_Biology", from = "06_machine_learning/publication", label = "Figure 1")),
    supplementary = list(list(stem = "FigureT_TherapyChange", from = "07_post_hoc/therapy", label = "Figure S2")),
    additional    = list(list(stem = "FigureS1_CONSORT", from = "06_machine_learning/publication"))))
idx <- suppressMessages(collect_publication(run, cfg))
pub <- file.path(run, "publication")
want <- c("main/Fig_1_Biology.pdf", "main/png/Fig_1_Biology.png", "main/tiff/Fig_1_Biology.tif",
          "supplementary/Fig_S2_TherapyChange.pdf",
          "supplementary/png/Fig_S2_TherapyChange_p01.png", "supplementary/png/Fig_S2_TherapyChange_p02.png",
          "supplementary/tiff/Fig_S2_TherapyChange_p01.tif", "additional/CONSORT.pdf", "INDEX.md")
for (f in want) check(file.exists(file.path(pub, f)), sprintf("collected file missing: %s", f))
check(!dir.exists(file.path(pub, "additional", "png")), "additional/ exported although not in export roles")
check(setequal(idx$file, setdiff(c(want, "supplementary/tiff/Fig_S2_TherapyChange_p02.tif"), "INDEX.md")),
      sprintf("INDEX rows: %s", paste(idx$file, collapse = ", ")))
check(abs(png_info(file.path(pub, "main/png/Fig_1_Biology.png"))$ppm - 150 / 0.0254) <= 1,
      "collected png does not carry the configured dpi")
check(any(grepl("main/png/Fig_1_Biology.png", readLines(file.path(pub, "INDEX.md")), fixed = TRUE)),
      "INDEX.md lacks the exported png")

# --- 3. config validation ----------------------------------------------------------
err <- function(x) inherits(try(figure_export_settings(x), silent = TRUE), "try-error")
check(is.null(figure_export_settings(list(figures = list(formats = "pdf")))), "no export block must mean no export")
check(err(list(figures = list(export = list(formats = "bmp")))), "unsupported format accepted")
check(err(list(figures = list(formats = c("pdf", "png"), export = list(formats = "png")))),
      "export format already rendered natively accepted")
check(err(list(figures = list(export = list(formats = "png", roles = "tables")))), "unknown role accepted")
check(err(list(figures = list(export = list(formats = "png", dpi = 0)))), "dpi 0 accepted")

unlink(work, recursive = TRUE)
if (length(fails)) {
  cat(sprintf("[FAIL] %s\n", fails), sep = "")
  quit(status = 1)
}
cat("[PASS] test_figure_export: formats, page suffixes, dpi, <role>/<format>/ layout, INDEX, validation\n")
