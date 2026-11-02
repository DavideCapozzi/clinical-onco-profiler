#!/usr/bin/env Rscript
# tools/survival_bor.R
# ==============================================================================
# OS / PFS by best overall response and by Δ immune score cut-off on a FINISHED run
# (R/modules_survival_bor.R, bor_survival_run()). Nothing is re-run: reads the run's
# step-06 frame and JSON, its step-07 km tests and the experiment's T0 input.
#
# Usage:
#   Rscript tools/survival_bor.R --run <run dir> [--config FILE] [--out DIR]
#     --config  default: the run's manifest config_path (its md5 is compared with the
#               manifest's and a mismatch is reported, not fatal)
#     --out     default: <run>/<exp>/07_post_hoc/bor_survival (one experiment only)
#
# Same pass config as step 07 (prepare_pass_config(..., "post_hoc")), so wiring the
# driver into src/07_post_hoc.R later gives identical output.
# ==============================================================================

suppressPackageStartupMessages({ library(here); library(jsonlite) })
setwd(here())
source("R/utils_io.R")
source("R/utils_figures.R")
source("R/utils_metrics.R")        # list_run_experiments()
source("R/modules_pub_style.R")
source("R/modules_survival_bor.R")

args <- commandArgs(trailingOnly = TRUE)
opt  <- list()
for (i in seq_along(args)) if (startsWith(args[[i]], "--") && i < length(args)) opt[[sub("^--", "", args[[i]])]] <- args[[i + 1]]
if (is.null(opt$run)) stop("[survival_bor] --run is required.")
run <- opt$run
if (!dir.exists(run) && dir.exists(file.path("results", run))) run <- file.path("results", run)
if (!dir.exists(run)) stop(sprintf("[survival_bor] run not found: %s", opt$run))
man <- file.path(run, "run_manifest.yml")
if (!file.exists(man)) stop(sprintf("[survival_bor] %s has no run_manifest.yml (unfinished run?)", run))
mf  <- yaml::read_yaml(man)
cfg_path <- if (!is.null(opt$config)) opt$config else mf$config_path
if (!file.exists(cfg_path)) stop(sprintf("[survival_bor] config not found: %s", cfg_path))
md5 <- unname(tools::md5sum(cfg_path))
if (!identical(md5, mf$config_md5))
  message(sprintf("[survival_bor] NOTE: %s md5 %s differs from the run's (%s): config edited after the run",
                  cfg_path, md5, mf$config_md5))

base <- yaml::read_yaml(cfg_path)
exps <- if (!is.null(base$experiments)) base$experiments else setNames(list(base), base$project_name)
exps <- exps[names(exps) %in% list_run_experiments(run)]
if (!length(exps)) stop(sprintf("[survival_bor] no experiment of %s found in %s", cfg_path, run))
if (!is.null(opt$out) && length(exps) > 1) stop("[survival_bor] --out needs a single-experiment run")
for (exp in names(exps)) {
  config <- prepare_pass_config(base, exps[[exp]], exp, "post_hoc", run)
  out <- if (!is.null(opt$out)) opt$out else file.path(run, exp, step_dir_name(7), "bor_survival")
  bor_survival_run(config, out, by = sprintf("tools/survival_bor.R --run %s (config %s, md5 %s)", run, cfg_path, md5))
  cat(sprintf("[OK] %s -> %s\n", exp, out))
}
