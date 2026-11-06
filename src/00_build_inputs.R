# src/00_build_inputs.R
# ==============================================================================
# Builds the analysis input files from the curated database.
#
#   Rscript src/00_build_inputs.R [--config config/data_build.yml]
#                                 [--subset <file with one patient ID per line>]
#                                 [--out <output dir>]
#
# Writes <out>/<prefix>_<sheet>.xlsx (analysis inputs), <out>/<prefix>_therapy_<sheet>.xlsx
# (therapy tables for step 07), build_log.txt and MD5SUMS. Run from the repository root.
# Rules: config/data_build.yml; documentation: docs/DATA_PROCESSING.md.
# ==============================================================================
suppressPackageStartupMessages({ library(yaml); library(readxl); library(openxlsx) })
source("R/modules_data_build.R")

args <- commandArgs(trailingOnly = TRUE)
opt <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i)) default else args[i + 1] }
cfg_path <- opt("--config", "config/data_build.yml")
cfg <- yaml::read_yaml(cfg_path)
out_dir <- opt("--out", cfg$output_dir)
subset_file <- opt("--subset")
subset <- if (!is.null(subset_file)) { x <- trimws(readLines(subset_file)); x[nzchar(x)] } else NULL

if (!file.exists(cfg$source)) stop(sprintf("[FATAL] curated database not found: %s", cfg$source))
if (dir.exists(out_dir) && length(list.files(out_dir)) && Sys.getenv("FORCE") != "1")
  stop(sprintf("[FATAL] %s is not empty (set FORCE=1 to overwrite)", out_dir))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_lines <- character()
say <- function(fmt, ...) { line <- sprintf(fmt, ...); message(line); log_lines <<- c(log_lines, line) }
say("build inputs  %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("config  %s (md5 %s)", cfg_path, unname(tools::md5sum(cfg_path)))
say("source  %s (md5 %s)", cfg$source, unname(tools::md5sum(cfg$source)))
if (!is.null(subset)) say("subset  %s (%d IDs)", subset_file, length(subset))

res <- db_build_inputs(cfg, subset = subset, say = say)

files <- character()
for (nm in names(res)) {
  f <- file.path(out_dir, sprintf("%s_%s.xlsx", cfg$output_prefix, nm))
  openxlsx::write.xlsx(res[[nm]], f, overwrite = TRUE)
  files <- c(files, f)
  say("written %s (%d x %d)", f, nrow(res[[nm]]), ncol(res[[nm]]))
}
writeLines(log_lines, file.path(out_dir, "build_log.txt"))
md5 <- tools::md5sum(files)
writeLines(sprintf("%s  %s", md5, basename(names(md5))), file.path(out_dir, "MD5SUMS"))
