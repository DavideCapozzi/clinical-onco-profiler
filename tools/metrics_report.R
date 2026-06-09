#!/usr/bin/env Rscript
# tools/metrics_report.R
# ==============================================================================
# Headline-metrics report for the /post-change workflow (and manual use):
#   - a publishability VERDICT for each experiment in a run, and
#   - an optional DELTA table versus a baseline run.
#
# All numbers and thresholds come from R/utils_metrics.R — the same single source
# of truth the golden tests assert on — so this script does no JSON path or
# threshold logic of its own (it only parses args, calls the helpers, and prints).
#
# Usage:
#   Rscript tools/metrics_report.R [--new <run>] [--base <run>] [--exp A,B]
#
#   --new   run root to evaluate         (default: latest, via resolve_run_root)
#   --base  baseline run root to diff    (omit -> verdict only, no delta table)
#   --exp   comma-separated experiments  (default: all in the run[s])
#
# Run roots may be absolute or relative to results/. Exit code 0 always (this is a
# report, not a gate — use tests/golden/* for pass/fail gating).
# ==============================================================================

suppressPackageStartupMessages({ library(here); library(jsonlite) })
source(here("R/utils_io.R"))
source(here("R/utils_metrics.R"))

# --- minimal --flag value arg parser -----------------------------------------
parse_args <- function(args) {
  out <- list()
  i <- 1
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

# Resolve a run argument to an existing directory (absolute, or under results/).
resolve_run_arg <- function(x) {
  if (dir.exists(x)) return(x)
  cand <- file.path("results", x)
  if (dir.exists(cand)) return(cand)
  stop(sprintf("[metrics_report] run not found: %s", x))
}

opt  <- parse_args(commandArgs(trailingOnly = TRUE))
new  <- if (!is.null(opt$new)) resolve_run_arg(opt$new) else resolve_run_root()
base <- if (!is.null(opt$base)) resolve_run_arg(opt$base) else NULL
exps <- if (!is.null(opt$exp)) trimws(strsplit(opt$exp, ",")[[1]]) else list_run_experiments(new)

# --- publishability verdicts --------------------------------------------------
cat(sprintf("\n=== PUBLISHABILITY — run %s ===\n", basename(new)))
for (exp in exps) {
  m <- extract_run_metrics(file.path(new, exp), exp)
  v <- publishability_verdict(m)
  if (is.na(v$publishable)) {
    cat(sprintf("\n%-22s  [N/A] %s\n", exp, v$reason)); next
  }
  cat(sprintf("\n%-22s  %s  (primary: %s)\n",
              exp, if (v$publishable) "PUBLISHABLE" else "NOT publishable",
              v$primary_method))
  for (c in v$criteria)
    cat(sprintf("   [%s] %-26s %s %s %s\n",
                if (isTRUE(c$pass)) "PASS" else "FAIL",
                c$metric, format(c$value), c$op, format(c$threshold)))
}

# --- delta table vs baseline --------------------------------------------------
if (!is.null(base)) {
  cat(sprintf("\n=== DELTA vs baseline %s ===\n", basename(base)))
  cmp <- compare_runs(base, new, experiments = if (!is.null(opt$exp)) exps else NULL)
  if (nrow(cmp) == 0) {
    cat("(no shared experiments to compare)\n")
  } else {
    cmp$flag <- ifelse(cmp$changed, "** CHANGED **", "=")
    # Fixed-width printer with truncation: long set-valued cells (e.g. the gate)
    # would otherwise make print.data.frame wrap the whole table vertically.
    cols   <- c("experiment", "metric", "base", "new", "delta", "flag")
    widths <- c(experiment = 20, metric = 15, base = 26, new = 26, delta = 8, flag = 13)
    trunc  <- function(s, n) { s <- as.character(s)
      ifelse(nchar(s) > n, paste0(substr(s, 1, n - 1), "…"), s) }
    pad    <- function(s, n) formatC(trunc(s, n), width = -n, flag = " ")
    line   <- function(vals) cat(paste(mapply(pad, vals, widths[cols]), collapse = " "), "\n")
    line(cols)
    line(strrep("-", widths[cols]))
    for (i in seq_len(nrow(cmp))) line(unlist(cmp[i, cols]))
    cat(sprintf("\n%d metric(s) changed beyond tolerance.\n", sum(cmp$changed, na.rm = TRUE)))
  }
}
