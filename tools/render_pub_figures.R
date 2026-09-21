#!/usr/bin/env Rscript
# tools/render_pub_figures.R
# ==============================================================================
# Re-render the publication figure set of an EXISTING run from its persisted
# `publication_data_<exp>.rds` — no analysis is re-run, so every number drawn is
# the run's own by construction. Use it after a figure-code (presentation) change.
#
# Usage:
#   Rscript tools/render_pub_figures.R [--run <run>] [--exp A,B] [--out DIR] [--preview DIR]
#
#   --run      run root (default: latest, via resolve_run_root); absolute or under results/
#   --exp      comma-separated experiments (default: all in the run)
#   --out      output dir (default: <run>/<exp>/06_machine_learning/publication_corrected;
#              with several experiments, <out>/<exp>)
#   --preview  also write PNG previews there (ragg)
#
# The run's own publication/ dir is never touched. Exits non-zero if any figure the
# renderer attempted came back NULL (save_pub_figure() would otherwise no-op silently).
# ==============================================================================

suppressPackageStartupMessages({ library(here); library(jsonlite) })
source(here("R/utils_io.R"))
source(here("R/utils_metrics.R"))        # list_run_experiments()
source(here("R/modules_pub_style.R"))
source(here("R/modules_pub_figures.R"))

parse_args <- function(args) {
  out <- list(); i <- 1
  while (i <= length(args)) {
    a <- args[[i]]
    if (startsWith(a, "--")) {
      key <- sub("^--", "", a)
      val <- if (i + 1 <= length(args) && !startsWith(args[[i + 1]], "--")) args[[i + 1]] else ""
      out[[key]] <- val
      i <- i + (if (nzchar(val)) 2 else 1)
    } else i <- i + 1
  }
  out
}
resolve_run_arg <- function(x) {
  if (dir.exists(x)) return(x)
  cand <- file.path("results", x)
  if (dir.exists(cand)) return(cand)
  stop(sprintf("[render_pub_figures] run not found: %s", x))
}

opt  <- parse_args(commandArgs(trailingOnly = TRUE))
run  <- if (!is.null(opt$run)) resolve_run_arg(opt$run) else resolve_run_root()
exps <- if (!is.null(opt$exp)) trimws(strsplit(opt$exp, ",")[[1]]) else list_run_experiments(run)
if (!is.null(opt$preview)) options(pub_fig.preview_dir = opt$preview)

failed <- character(0)
for (exp in exps) {
  exp_root <- file.path(run, exp)
  cfg      <- list(output_root = exp_root, project_name = exp)   # enough for step_dir()
  rds      <- file.path(step_dir(cfg, 6), "publication",
                        sprintf("publication_data_%s.rds", exp))
  if (!file.exists(rds)) { message(sprintf("[skip] %s: no %s", exp, basename(rds))); next }
  objs <- readRDS(rds)
  # Runs persisted before lmm_frame was plumbed into Step 06: read it from the same
  # run's Step-04 JSON, exactly as the pipeline now does.
  if (is.null(objs$lmm_frame)) objs$lmm_frame <- pub_read_lmm_frame(cfg)

  out <- if (is.null(opt$out)) file.path(step_dir(cfg, 6), "publication_corrected")
         else if (length(exps) > 1) file.path(opt$out, exp) else opt$out
  cat(sprintf("\n=== %s -> %s ===\n", exp, out))
  log <- pub_render_all(objs, out, exp)
  print(log, row.names = FALSE)
  failed <- c(failed, sprintf("%s/%s", exp, log$figure[!log$saved]))
}
if (length(failed)) {
  cat(sprintf("\n[FAIL] figure(s) returned NULL: %s\n", paste(failed, collapse = ", ")))
  quit(status = 1)
}
cat("\n[OK] every attempted figure was written.\n")
