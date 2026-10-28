#!/usr/bin/env Rscript
# tools/render_pub_figures.R
# ==============================================================================
# Re-render the publication figure set of an EXISTING run from its persisted
# `publication_data_<exp>.rds` — no analysis is re-run, so every number drawn is
# the run's own by construction. Use it after a figure-code (presentation) change.
#
# Usage:
#   Rscript tools/render_pub_figures.R --run <run> [--exp A,B] [--config FILE] [--out DIR] [--preview DIR]
#
#   --run      run root (REQUIRED); absolute or under results/
#   --exp      comma-separated experiments (default: all in the run)
#   --config   config to use (default: the run manifest's config_path; its md5 must match
#              the manifest's config_md5, else pass --config explicitly)
#   --out      render into DIR only (with several experiments, <out>/<exp>) and do NOT touch
#              the run. Default: re-render into <run>/<exp>/06_machine_learning/publication
#              and refresh <run>/publication/ through collect_publication().
#   --preview  also write PNG previews there (ragg)
#
# Formats and dpi come from the config's figures block. Exits non-zero if any figure the
# renderer attempted came back NULL (save_pub_figure() would otherwise no-op silently).
# ==============================================================================

suppressPackageStartupMessages({ library(here); library(jsonlite) })
source(here("R/utils_io.R"))
source(here("R/utils_figures.R"))
source(here("R/utils_metrics.R"))        # list_run_experiments()
source(here("R/modules_pub_style.R"))
source(here("R/modules_pub_figures.R"))
source(here("R/utils_publication.R"))     # collect_publication()

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
if (is.null(opt$run) || !nzchar(opt$run)) stop("[render_pub_figures] --run <run> is required (there is no latest-run default).")
run  <- resolve_run_arg(opt$run)
exps <- if (!is.null(opt$exp)) trimws(strsplit(opt$exp, ",")[[1]]) else list_run_experiments(run)
if (!is.null(opt$preview)) options(pub_fig.preview_dir = opt$preview)

# The run's config: formats/dpi for the figures and the publication mapping for the collector.
cfg_path <- opt$config
if (is.null(cfg_path) || !nzchar(cfg_path)) {
  man <- file.path(run, "run_manifest.yml")
  if (!file.exists(man)) stop("[render_pub_figures] no run_manifest.yml in the run: pass --config")
  m <- yaml::read_yaml(man)
  cfg_path <- m$config_path
  if (is.null(cfg_path) || !file.exists(cfg_path))
    stop("[render_pub_figures] the manifest's config_path is missing or unreadable: pass --config")
  if (!identical(unname(tools::md5sum(cfg_path)), m$config_md5))
    stop(sprintf("[render_pub_figures] %s changed since the run (md5 differs from the manifest): pass --config", cfg_path))
}
run_cfg <- yaml::read_yaml(cfg_path)

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

  out <- if (is.null(opt$out)) file.path(step_dir(cfg, 6), "publication")
         else if (length(exps) > 1) file.path(opt$out, exp) else opt$out
  cat(sprintf("\n=== %s -> %s ===\n", exp, out))
  log <- pub_render_all(objs, out, exp, fig_cfg = run_cfg)
  print(log, row.names = FALSE)
  failed <- c(failed, sprintf("%s/%s", exp, log$figure[!log$saved]))
}
if (length(failed)) {
  cat(sprintf("\n[FAIL] figure(s) returned NULL: %s\n", paste(failed, collapse = ", ")))
  quit(status = 1)
}
cat("\n[OK] every attempted figure was written.\n")
if (is.null(opt$out)) {
  collect_publication(run, run_cfg)
  cat(sprintf("[OK] %s refreshed through collect_publication().\n", file.path(run, "publication")))
}
