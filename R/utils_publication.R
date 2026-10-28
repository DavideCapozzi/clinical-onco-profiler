# R/utils_publication.R
# ==============================================================================
# PUBLICATION COLLECTOR — copies the paper's figures and tables from the step
# directories into <run>/publication/{main, supplementary, supplementary/tables,
# additional}/ and writes publication/INDEX.md (one row per file).
# The mapping lives in config `publication:` (role -> list of entries); each entry:
#   stem      file stem (the step output is <stem>_<EXP>.<ext> or <stem>.<ext>)
#   from      directory under <run>/<EXP>/ (e.g. 06_machine_learning/publication)
#   type      figure (default; formats from config figures.formats) | xlsx | csv
#   label     paper label (provisional; numbering may change); names the copy:
#             "Figure 2" -> Fig_2_<name>, "Figure S1" -> Fig_S1_<name>, "Table S1" -> Table_S1_<name>
#   name      descriptive part of the copied file name (default: the stem without its
#             legacy Figure<N>_ prefix); entries without a label are named <name> only
#   producer  producing function (for INDEX.md)
#   frame     cohort-frame text; {key} placeholders are filled from the run's JSONs
#             (keys of collect_frame_ns() plus t0_pdl1_valid)
# The collector COPIES, so the step directories stay complete (step file names unchanged;
# INDEX.md maps every copy to its source). A missing main or supplementary file is FATAL;
# a missing additional file is a warning.
# ==============================================================================

PUBLICATION_ROLES <- c(main = "main", supplementary = "supplementary",
                       tables = "supplementary/tables", additional = "additional")

#' @title Fill {key} placeholders of a frame description
.pub_fill_frame <- function(txt, ns) {
  if (is.null(txt) || !nzchar(txt)) return("")
  for (k in names(ns)) txt <- gsub(sprintf("{%s}", k), as.character(ns[[k]]), txt, fixed = TRUE)
  txt
}

#' @title File-name stem of a collected copy (without extension / page suffix)
#' @param e A mapping entry (uses \code{stem}, \code{label}, \code{name}).
#' @return E.g. \code{"Fig_2_AddedValue_CVcurves"}; \code{"Coupling"} for an unlabelled entry.
pub_dest_stem <- function(e) {
  nm <- if (!is.null(e$name)) e$name else sub("^Figure[A-Z]*[0-9]*[a-z]?_", "", e$stem)
  if (is.null(e$label) || !nzchar(e$label)) return(nm)
  paste(gsub("\\s+", "_", sub("^Figure", "Fig", trimws(e$label))), nm, sep = "_")
}

#' @title Name of the copy of one step file: dest stem + the file's page suffix / extension
.pub_dest_name <- function(f, e, exp) {
  b <- basename(f)
  for (src in c(sprintf("%s_%s", e$stem, exp), e$stem))
    if (startsWith(b, src)) return(paste0(pub_dest_stem(e), substring(b, nchar(src) + 1)))
  stop(sprintf("[FATAL] publication: %s does not start with stem %s", b, e$stem))
}

#' @title Files a mapping entry resolves to in a step directory
.pub_entry_files <- function(src_dir, stem, exp, type, formats) {
  exts <- switch(type, figure = formats, xlsx = "xlsx", csv = "csv",
                 stop(sprintf("[FATAL] publication: unknown entry type '%s'", type)))
  unlist(lapply(exts, function(ext) {
    cand <- file.path(src_dir, c(sprintf("%s_%s.%s", stem, exp, ext), sprintf("%s.%s", stem, ext)))
    hit <- cand[file.exists(cand)]
    if (length(hit)) return(hit[1])
    # multipage raster / svg export: <stem>[_<exp>]_pNN.<ext>
    pages <- list.files(src_dir, pattern = sprintf("^%s(_%s)?_p[0-9]+\\.%s$", stem, exp, ext), full.names = TRUE)
    if (length(pages)) pages else NA_character_
  }))
}

#' @title Collect the publication set of a finished run
#' @param run_dir Run root (\code{<output_root>/<run name>}).
#' @param config Base config (uses \code{publication} and \code{figures}).
#' @param experiments Experiments to collect (default: all with step folders).
#' @return Invisibly, the INDEX data.frame.
collect_publication <- function(run_dir, config, experiments = NULL) {
  map <- config$publication
  if (is.null(map) || !length(map)) {
    message("[Publication] no `publication:` mapping in config — nothing collected."); return(invisible(NULL))
  }
  bad <- setdiff(names(map), names(PUBLICATION_ROLES))
  if (length(bad)) stop(sprintf("[FATAL] publication: unknown role(s) %s", paste(bad, collapse = ", ")))
  formats <- figure_settings(config)$formats
  exps <- if (!is.null(experiments)) experiments else list_run_experiments(run_dir)
  if (!length(exps)) stop(sprintf("[FATAL] publication: no experiment directories under %s", run_dir))

  pub_root <- file.path(run_dir, "publication")
  if (dir.exists(pub_root)) unlink(pub_root, recursive = TRUE)   # always a fresh copy of the step outputs
  rows <- list(); missing_req <- character(0)
  for (exp in exps) {
    exp_root <- file.path(run_dir, exp)
    base <- if (length(exps) > 1) file.path(pub_root, exp) else pub_root
    ns <- collect_frame_ns(exp_root, exp)
    j6 <- file.path(exp_root, step_dir_name(6), sprintf("Machine_Metrics_ML_%s.json", exp))
    if (file.exists(j6)) {
      nv <- tryCatch(jsonlite::fromJSON(j6, simplifyVector = FALSE)$clinical_benchmark$n_valid, error = function(e) NULL)
      if (!is.null(nv)) ns$t0_pdl1_valid <- as.integer(nv)
    }
    for (role in names(PUBLICATION_ROLES)) {
      entries <- map[[role]]
      if (!length(entries)) next
      dest <- file.path(base, PUBLICATION_ROLES[[role]])
      dir.create(dest, recursive = TRUE, showWarnings = FALSE)
      for (e in entries) {
        type <- if (!is.null(e$type)) e$type else "figure"
        files <- .pub_entry_files(file.path(exp_root, e$from), e$stem, exp, type, formats)
        found <- files[!is.na(files)]
        if (!length(found) || anyNA(files)) {
          msg <- sprintf("%s/%s (%s) missing in %s", role, e$stem, type, file.path(exp, e$from))
          if (role == "additional") { warning(sprintf("[Publication] %s", msg)); if (!length(found)) next }
          else { missing_req <- c(missing_req, msg); next }
        }
        out <- vapply(found, .pub_dest_name, character(1), e = e, exp = exp, USE.NAMES = FALSE)
        clash <- file.exists(file.path(dest, out))
        if (any(clash)) stop(sprintf("[FATAL] publication: two entries map to %s", paste(out[clash], collapse = ", ")))
        ok <- file.copy(found, file.path(dest, out), copy.date = TRUE)
        if (!all(ok)) stop(sprintf("[FATAL] publication: copy failed for %s", paste(found[!ok], collapse = ", ")))
        for (i in seq_along(found))
          rows[[length(rows) + 1]] <- data.frame(
            experiment = exp, role = PUBLICATION_ROLES[[role]],
            label = if (!is.null(e$label)) e$label else "—",
            file = if (length(exps) > 1) file.path(exp, PUBLICATION_ROLES[[role]], out[i])
                   else file.path(PUBLICATION_ROLES[[role]], out[i]),
            source = file.path(exp, e$from, basename(found[i])),
            producer = if (!is.null(e$producer)) e$producer else "—",
            frame = .pub_fill_frame(e$frame, ns), stringsAsFactors = FALSE)
      }
    }
  }
  if (length(missing_req))
    stop(sprintf("[FATAL] publication: mapped main/supplementary file(s) missing: %s",
                 paste(missing_req, collapse = "; ")))

  idx <- do.call(rbind, rows)
  esc <- function(x) gsub("|", "\\|", x, fixed = TRUE)
  lines <- c(
    sprintf("# Publication set — %s", basename(normalizePath(run_dir))),
    "",
    sprintf("Collected %s by `collect_publication()` (R/utils_publication.R) from the step directories (copies;",
            format(Sys.time(), "%Y-%m-%d %H:%M")),
    "the step outputs stay in place). Paper labels are provisional: numbering may still change.",
    sprintf("Figure formats: %s. Frames are named with the n of this run.", paste(formats, collapse = ", ")),
    "",
    "| role | paper label | file | source (step output) | producing function | n / frame |",
    "|---|---|---|---|---|---|",
    sprintf("| %s | %s | `%s` | `%s` | %s | %s |", esc(idx$role), esc(idx$label), esc(idx$file),
            esc(idx$source), esc(idx$producer), esc(idx$frame)),
    "")
  writeLines(lines, file.path(pub_root, "INDEX.md"))
  message(sprintf("[Publication] %d file(s) collected -> %s (INDEX.md)", nrow(idx), pub_root))
  invisible(idx)
}
