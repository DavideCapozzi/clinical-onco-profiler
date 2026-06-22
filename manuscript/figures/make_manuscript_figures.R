#!/usr/bin/env Rscript
# ==============================================================================
# manuscript/figures/make_manuscript_figures.R  (FAST RE-RENDER, no re-analysis)
# ------------------------------------------------------------------------------
# The manuscript figures are produced by the pipeline (Step 06) from the live
# result objects, via R/modules_pub_style.R + R/modules_pub_figures.R. Step 06
# also persists those objects to
#   <run>/<experiment>/06_machine_learning/publication/publication_data_<exp>.rds
#
# This script RE-RENDERS the whole set from that bundle using the SAME builders
# (pub_render_all) — so you can iterate on figure sizes/styling in seconds
# WITHOUT re-running the (slow) analysis. Output also copied to
# manuscript/figures/output/ for manuscript assembly. One source of plotting
# logic + sizes (R/modules_pub_*); this is just a thin driver.
#
# USAGE:
#   Rscript manuscript/figures/make_manuscript_figures.R
#   RUN_DIR=/abs/.../<run>/BestResponse_2v3_4 Rscript manuscript/figures/make_manuscript_figures.R
# ==============================================================================

suppressMessages({ library(here); library(jsonlite) })
repo <- tryCatch(here::here(), error = function(e) getwd())
invisible(lapply(list.files(file.path(repo, "R"), pattern = "modules_pub_.*\\.R$",
                            full.names = TRUE), source))

EXP <- Sys.getenv("EXPERIMENT", "BestResponse_2v3_4")
run_dir <- Sys.getenv("RUN_DIR", "")
if (!nzchar(run_dir)) {
  latest <- file.path(repo, "results", "latest")
  if (!dir.exists(latest)) stop("No results/latest; set RUN_DIR=")
  run_dir <- file.path(normalizePath(latest), EXP)
}
pub_dir <- file.path(run_dir, "06_machine_learning", "publication")
rds <- file.path(pub_dir, sprintf("publication_data_%s.rds", EXP))
if (!file.exists(rds))
  stop(sprintf("No %s\nRun the pipeline / Step 06 with run_publication_figures: true first.", rds))

objs <- readRDS(rds)
out  <- file.path(repo, "manuscript", "figures", "output")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
unlink(file.path(out, "*.pdf"))

# re-render at the sizes defined in R/modules_pub_figures.R (single source)
pub_render_all(objs, pub_dir, EXP)          # refresh the canonical run set
pub_render_all(objs, out,     EXP)          # + manuscript assembly copy

cat(sprintf("\nRe-rendered manuscript figures from:\n  %s\nto:\n  %s\n", rds, out))
for (f in sort(basename(list.files(out, pattern = "\\.pdf$")))) cat("   -", f, "\n")
cat("\nCaptions: manuscript/figure_captions.md\n")
