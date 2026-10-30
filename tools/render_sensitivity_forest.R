#!/usr/bin/env Rscript
# tools/render_sensitivity_forest.R
# ==============================================================================
# Draw the LMM Time×Group forest of a SENSITIVITY run (e.g. nlr_v2_cc50: gate selected
# on its own LMM frame) into an EXISTING run, so the final run's publication set can
# carry it as a supplementary figure. Nothing is re-run: the forest is read from the
# sensitivity run's Step-04 report (pub_read_lmm_bootstrap) and drawn with the same
# builder as Figure 1 panel A (pub_fig_lmm_forest), restricted to that run's own gate.
#
# Usage:
#   Rscript tools/render_sensitivity_forest.R --run <run> --from <sensitivity run> [--config FILE]
#
#   --run     run that receives the figure (REQUIRED)
#   --from    run whose Step-04 LMM is drawn (REQUIRED)
#   --config  config for figure formats / dpi (default: --run's manifest config_path)
#
# Writes <run>/<exp>/07_post_hoc/sensitivity_lmm/FigureSens_LMMForest.<ext> + source.txt
# (source report, md5, frame). Then re-collect with tools/render_pub_figures.R --run <run>.
# After a full re-run of <run> this must be run again BEFORE the publication collector
# finds the file: until then collect_publication() stops on the missing supplementary.
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
run  <- resolve_run(opt$run)
from <- resolve_run(opt$from)
if (is.null(run) || is.null(from)) stop("[render_sensitivity_forest] --run and --from are both required.")

cfg_path <- if (!is.null(opt$config)) opt$config else yaml::read_yaml(file.path(run, "run_manifest.yml"))$config_path
fig_cfg  <- yaml::read_yaml(cfg_path)

for (exp in intersect(list_run_experiments(run), list_run_experiments(from))) {
  from_cfg <- list(output_root = file.path(from, exp), project_name = exp)
  lmm <- pub_read_lmm_bootstrap(from_cfg)
  frm <- pub_read_lmm_frame(from_cfg)
  rds <- file.path(step_dir(from_cfg, 6, create = FALSE), "publication", sprintf("publication_data_%s.rds", exp))
  if (is.null(lmm) || is.null(frm) || !file.exists(rds)) stop(sprintf("[render_sensitivity_forest] %s: Step-04 report, LMM JSON or Step-06 rds missing in %s", exp, from))
  objs    <- readRDS(rds)
  markers <- unique(objs$gate_decomp$per_patient_values$Marker)     # that run's own gate
  ref     <- pub_relabel_group(objs$clin_addval$negative_label)
  p <- pub_fig_lmm_forest(lmm$boot, lmm$obs, markers = markers, ref_label = ref)
  if (is.null(p)) stop(sprintf("[render_sensitivity_forest] %s: forest came back NULL", exp))

  out_dir <- file.path(run, exp, step_dir_name(7), "sensitivity_lmm")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  save_pub_figure(p, file.path(out_dir, "FigureSens_LMMForest"), PUB_W1, 95, config = fig_cfg)
  src <- file.path(step_dir(from_cfg, 4, create = FALSE), sprintf("Longitudinal_LMM_Report_%s.xlsx", exp))
  writeLines(c(
    sprintf("written: %s by tools/render_sensitivity_forest.R", format(Sys.time(), "%Y-%m-%d %H:%M")),
    sprintf("source: %s (md5 %s)", src, unname(tools::md5sum(src))),
    sprintf("frame: LMM %d patients / %d obs (%d paired) of %s", frm$n_patients, frm$n_observations,
            frm$n_paired, basename(normalizePath(from))),
    sprintf("markers: %s; reference group: %s", paste(markers, collapse = ", "), ref)),
    file.path(out_dir, "source.txt"))
  cat(sprintf("[OK] %s: forest of %s (LMM %d / %d) -> %s\n", exp, basename(normalizePath(from)),
              frm$n_patients, frm$n_observations, out_dir))
}
