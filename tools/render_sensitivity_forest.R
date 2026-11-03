#!/usr/bin/env Rscript
# tools/render_sensitivity_forest.R
# ==============================================================================
# LMM Time×Group forest of a SENSITIVITY run (e.g. nlr_v2_cc50: gate selected on its own
# LMM frame), for another run's publication set (the final run's Figure S1). Nothing is
# re-run: the forest is the sensitivity run's Step-04 report, drawn with the same builder
# as Figure 1 panel A (pub_fig_lmm_forest), restricted to that run's own gate.
#
# The inputs are FROZEN in a snapshot (byte copies of the Step-04 report + LMM JSON,
# SOURCE.yml with the gate and reference group, MD5SUMS), so the receiving run no longer
# depends on the sensitivity run existing. Step 07 draws the figure itself when the
# config has post_hoc.sensitivity_lmm.snapshot (as config/final_gate84_nlr50.yml does).
#
# Usage:
#   Rscript tools/render_sensitivity_forest.R --from <sensitivity run> --snapshot <dir>
#       freeze the inputs (once; an existing snapshot is never overwritten)
#   Rscript tools/render_sensitivity_forest.R --run <run> --snapshot <dir> [--config FILE]
#       draw into a finished run without re-running it (= what step 07 does), then
#       re-collect with tools/render_pub_figures.R --run <run>
#
#   --config  config for figure formats / dpi (default: --run's manifest config_path)
#
# Writes <run>/<exp>/07_post_hoc/sensitivity_lmm/FigureSens_LMMForest.<ext> + source.txt.
# ==============================================================================

suppressPackageStartupMessages({ library(here); library(jsonlite) })
source(here("R/utils_io.R"))
source(here("R/utils_figures.R"))
source(here("R/utils_metrics.R"))        # list_run_experiments()
source(here("R/modules_pub_style.R"))
source(here("R/modules_pub_figures.R"))

args <- commandArgs(trailingOnly = TRUE)
opt  <- list()
for (i in seq_along(args)) if (startsWith(args[[i]], "--") && i < length(args)) opt[[sub("^--", "", args[[i]])]] <- args[[i + 1]]
resolve_run <- function(x) {
  if (is.null(x) || !nzchar(x)) return(NULL)
  for (cand in c(x, file.path("results", x))) if (dir.exists(cand)) return(cand)
  stop(sprintf("[render_sensitivity_forest] run not found: %s", x))
}
snap <- opt$snapshot
if (is.null(snap) || !nzchar(snap)) stop("[render_sensitivity_forest] --snapshot is required.")
if (!is.null(opt$from) == !is.null(opt$run)) stop("[render_sensitivity_forest] give exactly one of --from (freeze) or --run (draw).")

if (!is.null(opt$from)) {
  from <- resolve_run(opt$from)
  exps <- list_run_experiments(from)
  if (length(exps) != 1) stop(sprintf("[render_sensitivity_forest] %s has %d experiments; one snapshot holds one", from, length(exps)))
  pub_snapshot_sensitivity_lmm(from, exps, snap)
  cat(sprintf("[OK] %s: inputs of the %s LMM forest frozen in %s (MD5SUMS)\n", exps, basename(normalizePath(from)), snap))
} else {
  run <- resolve_run(opt$run)
  cfg_path <- if (!is.null(opt$config)) opt$config else yaml::read_yaml(file.path(run, "run_manifest.yml"))$config_path
  fig_cfg  <- yaml::read_yaml(cfg_path)
  exp <- yaml::read_yaml(file.path(snap, "SOURCE.yml"))$experiment
  if (!exp %in% list_run_experiments(run)) stop(sprintf("[render_sensitivity_forest] %s has no experiment %s", run, exp))
  out_dir <- file.path(run, exp, step_dir_name(7), "sensitivity_lmm")
  r <- pub_draw_sensitivity_lmm(snap, exp, out_dir, fig_cfg = fig_cfg, by = "tools/render_sensitivity_forest.R")
  cat(sprintf("[OK] %s: forest of %s (LMM %d / %d) -> %s\n", exp, r$source_run,
              r$frame$n_patients, r$frame$n_observations, out_dir))
}
