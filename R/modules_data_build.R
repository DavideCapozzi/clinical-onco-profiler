# R/modules_data_build.R
# ==============================================================================
# Data build: curated database -> analysis input files (src/00_build_inputs.R).
#
# Every rule is generic and driven by config/data_build.yml. Patient-level
# content (corrections, supplements, exclusions) is data and lives in the
# curated database; nothing patient-specific is written here.
#
# Order of operations per sheet (the baseline sheet first):
#   1. read as text; drop rows without a patient ID
#   2. rename / clean / alias the headers; split off exclusion and therapy columns
#   3. dates, derived binary variables, metastatic-site count, PD-L1 parsing
#   4. numeric coercion of every non-text column (logged)
#   5. row filter on missing values, duplicate records, baseline back-fill
#   6. exclusions, optional subset, output schema
# ==============================================================================

#' Clean a source header into a pipeline column name
#' (drops "%", turns spaces and + - ( , ) = into "_", collapses repeats).
db_clean_name <- function(x) {
  x <- gsub("%", "", x)
  x <- gsub("[+[:space:](,)=-]+", "_", x)
  x <- gsub("_+", "_", x)
  sub("_$", "", x)
}

#' Normalised person columns (lower case, single spaces; blank = NA), for record matching only.
db_person <- function(df, cols) {
  as.data.frame(lapply(setNames(cols, cols), function(cn) {
    v <- tolower(gsub("[[:space:]]+", " ", trimws(df[[cn]]))); v[!is.na(v) & !nzchar(v)] <- NA; v
  }), check.names = FALSE, stringsAsFactors = FALSE)
}

#' Two records belong to the same person when the first person column (family name) is equal
#' and every other person column is equal where both records have it (a blank is not a mismatch).
db_same_person <- function(a, b) {
  same <- !is.na(a[[1]]) & !is.na(b[[1]]) & a[[1]] == b[[1]]
  for (cn in names(a)[-1]) same <- same & (is.na(a[[cn]]) | is.na(b[[cn]]) | a[[cn]] == b[[cn]])
  same
}

#' Binary flag from a free-text column: zero tokens -> "0", blank or "na" -> NA, any other text -> "1".
db_binary_flag <- function(x, zero) {
  low <- tolower(trimws(x))
  ifelse(is.na(x) | low == "na", NA_character_, ifelse(low %in% zero, "0", "1"))
}

#' Number of distinct metastatic-site codes (integers in [min_code, max_code]) in a free-text cell.
db_count_codes <- function(x, min_code = 0L, max_code = 10L) {
  vapply(x, function(v) {
    if (is.na(v) || trimws(v) == "") return(NA_integer_)
    codes <- suppressWarnings(as.integer(regmatches(v, gregexpr("\\b\\d+\\b", v))[[1]]))
    codes <- codes[!is.na(codes) & codes >= min_code & codes <= max_code]
    if (length(codes) == 0L) NA_integer_ else length(unique(codes))
  }, integer(1), USE.NAMES = FALSE)
}

#' PD-L1 TPS free text -> percentage.
#' Plain number, "N:" category prefix and "%" are stripped; "X-Y" or "X/Y" -> midpoint;
#' "<X" -> X/2; ">X" -> (X+100)/2; anything else -> NA.
db_parse_pdl1 <- function(x) {
  vapply(as.character(x), function(s) {
    s <- trimws(s)
    if (is.na(s) || !nzchar(s)) return(NA_real_)
    s <- trimws(gsub("%", "", sub("^\\s*\\d+\\s*:\\s*", "", s)))
    v <- suppressWarnings(as.numeric(s)); if (!is.na(v)) return(v)
    m <- regmatches(s, regexpr("^<\\s*(\\d+(?:\\.\\d+)?)", s, perl = TRUE))
    if (length(m)) return(as.numeric(sub("^<\\s*", "", m)) / 2)
    m <- regmatches(s, regexpr("^>\\s*(\\d+(?:\\.\\d+)?)", s, perl = TRUE))
    if (length(m)) return((as.numeric(sub("^>\\s*", "", m)) + 100) / 2)
    m <- regmatches(s, regexpr("^(\\d+(?:\\.\\d+)?)\\s*[-/]\\s*(\\d+(?:\\.\\d+)?)", s, perl = TRUE))
    if (length(m)) return(mean(as.numeric(strsplit(m, "\\s*[-/]\\s*")[[1L]])))
    NA_real_
  }, numeric(1), USE.NAMES = FALSE)
}

#' Columns of `nms` matched by `patterns` (entries starting with "^" are regular expressions).
db_match_cols <- function(nms, patterns) {
  unique(unlist(lapply(patterns, function(p)
    if (startsWith(p, "^")) grep(p, nms, value = TRUE) else intersect(p, nms))))
}

#' Read and standardise one sheet of the curated database (steps 1-4).
#' @return list(data, person, excluded, therapy, log)
db_prepare_sheet <- function(path, sheet, cfg, say) {
  cc <- cfg$columns
  raw <- as.data.frame(readxl::read_excel(path, sheet = sheet, col_types = "text", .name_repair = "minimal"),
                       check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(cc$id %in% names(raw), all(cc$person %in% names(raw)), cc$exclusion %in% names(raw))
  keep <- !is.na(raw[[cc$id]]) & nzchar(trimws(raw[[cc$id]]))
  say("[%s] %d patient rows (%d rows without an ID ignored)", sheet, sum(keep), sum(!keep))
  raw <- raw[keep, , drop = FALSE]

  person   <- db_person(raw, cc$person)
  excluded <- data.frame(id = raw[[cc$id]], flag = as.numeric(raw[[cc$exclusion]]),
                         reason = raw[[cc$exclusion_reason]], stringsAsFactors = FALSE)
  th_map  <- cfg$therapy[[sheet]]
  therapy <- if (length(th_map)) {
    stopifnot(all(names(th_map) %in% names(raw)))
    t <- data.frame(raw[[cc$id]], lapply(raw[names(th_map)], as.numeric), check.names = FALSE)
    setNames(t, c(cc$id, unname(unlist(th_map))))
  } else NULL
  df <- raw[, setdiff(names(raw), c(cc$person, cc$exclusion, cc$exclusion_reason, names(th_map))), drop = FALSE]

  # headers: explicit renames, then cleaning, then aliases
  nm <- names(df)
  for (from in names(cc$rename)) nm[nm == from] <- cc$rename[[from]]
  nm <- ifelse(nm == cc$id, nm, db_clean_name(nm))
  for (from in names(cc$aliases)) if (from %in% nm && !cc$aliases[[from]] %in% nm) nm[nm == from] <- cc$aliases[[from]]
  stopifnot(!anyDuplicated(nm[!nm %in% unlist(lapply(cfg$derived, `[[`, "name"))]))
  names(df) <- nm

  # dates (Excel serial numbers -> dd/mm/yyyy; other text kept)
  for (cn in intersect(cc$dates, names(df))) {
    v <- df[[cn]]; num <- !is.na(v) & grepl("^[0-9]+$", v)
    v[num] <- format(openxlsx::convertToDate(as.numeric(v[num])), "%d/%m/%Y")
    df[[cn]] <- v
  }
  # derived binary variables; a supplied column of the same name fills undetermined cells
  for (d in cfg$derived) {
    src <- setdiff(db_match_cols(names(df), d$from), d$name); stopifnot(length(src) == 1)
    flag <- db_binary_flag(df[[src]], unlist(d$zero))
    if (d$name %in% names(df)) {
      n_fill <- sum(is.na(flag) & !is.na(df[[d$name]]))
      flag <- ifelse(is.na(flag), df[[d$name]], flag)
      say("[%s] %s: derived from %s; %d undetermined cell(s) taken from the supplied column", sheet, d$name, src, n_fill)
      df[[d$name]] <- NULL
    }
    df[[d$name]] <- flag
  }
  mc  <- cfg$mets_count
  src <- db_match_cols(names(df), mc$from); stopifnot(length(src) == 1)
  df[[mc$name]] <- db_count_codes(df[[src]], mc$min_code, mc$max_code)
  # PD-L1
  if (cc$pdl1 %in% names(df)) {
    p <- db_parse_pdl1(df[[cc$pdl1]])
    say("[%s] %s: %d free-text value(s) parsed, %d missing", sheet, cc$pdl1,
        sum(!is.na(p) & is.na(suppressWarnings(as.numeric(df[[cc$pdl1]])))), sum(is.na(p)))
    df[[cc$pdl1]] <- p
  }
  # numeric coercion of every non-text column
  text_cols <- c(cc$id, db_match_cols(names(df), c(unlist(cc$text), unlist(cc$dates))),
                 unlist(lapply(cfg$derived, `[[`, "name")))
  for (cn in setdiff(names(df), c(text_cols, cc$pdl1, mc$name))) {
    v <- df[[cn]]; n <- suppressWarnings(as.numeric(v))
    lost <- !is.na(v) & is.na(n)
    if (any(lost)) say("[%s] %s: %d non-numeric cell(s) set to missing", sheet, cn, sum(lost))
    df[[cn]] <- n
  }
  list(data = df, person = person, excluded = excluded, therapy = therapy)
}

#' Build the analysis input files from the curated database.
#' @param cfg parsed config/data_build.yml
#' @param subset optional character vector of patient IDs to keep
#' @param say logger: function(fmt, ...)
#' @return list of data frames: <sheet> (analysis) and therapy_<sheet>
db_build_inputs <- function(cfg, subset = NULL, say = function(fmt, ...) message(sprintf(fmt, ...))) {
  cc <- cfg$columns; id <- cc$id
  sheets <- unlist(cfg$sheets); base <- sheets[1]
  P <- lapply(setNames(sheets, sheets), function(s) db_prepare_sheet(cfg$source, s, cfg, say))
  pl <- unlist(cfg$patient_level)
  baseline <- P[[base]]$data[, c(id, intersect(pl, names(P[[base]]$data))), drop = FALSE]
  base_person <- P[[base]]$person

  out <- list()
  for (s in sheets) {
    df <- P[[s]]$data; person <- P[[s]]$person
    schema <- unlist(cfg$output_columns[[s]])
    miss <- setdiff(schema, names(df))
    if (length(miss)) stop(sprintf("[%s] output columns not found: %s", s, paste(miss, collapse = ", ")))
    # rows with too many missing values (share over the output columns)
    na_share <- rowMeans(is.na(df[, schema, drop = FALSE]))
    drop <- na_share > cfg$qc$max_na_row
    if (any(drop)) say("[%s] removed %d row(s) with > %.0f%% missing values: %s", s, sum(drop),
                       100 * cfg$qc$max_na_row, paste(df[[id]][drop], collapse = ", "))
    df <- df[!drop, , drop = FALSE]; person <- person[!drop, , drop = FALSE]
    # duplicate records of the same person (all person columns equal): keep the first
    key <- do.call(paste, c(lapply(person, function(v) ifelse(is.na(v), "", v)), sep = "\r"))
    dup <- duplicated(key) & !is.na(person[[1]])
    if (any(dup)) say("[%s] removed %d duplicate record(s) of the same person: %s", s, sum(dup),
                      paste(df[[id]][dup], collapse = ", "))
    df <- df[!dup, , drop = FALSE]; person <- person[!dup, , drop = FALSE]
    # patient-level attributes: blank follow-up cells take the baseline value of the same person
    if (s != base) for (cn in intersect(pl, names(df))) {
      j <- match(df[[id]], baseline[[id]])
      same <- !is.na(j) & db_same_person(base_person[j, , drop = FALSE], person)
      src <- baseline[[cn]][j]
      fill <- is.na(df[[cn]]) & !is.na(src) & same
      if (any(fill)) { df[[cn]][fill] <- src[fill]
        say("[%s] %s: %d blank cell(s) taken from %s: %s", s, cn, sum(fill), base, paste(df[[id]][fill], collapse = ", ")) }
      blocked <- is.na(df[[cn]]) & !is.na(src) & !same
      if (any(blocked)) say("[%s] %s: not filled, %s row is another person: %s", s, cn, base,
                            paste(df[[id]][blocked], collapse = ", "))
    }
    # exclusions
    ex <- P[[s]]$excluded; ex_ids <- ex$id[ex$flag %in% 1]
    hit <- df[[id]] %in% ex_ids
    if (any(hit)) say("[%s] excluded %d patient(s) flagged in the curated database: %s", s, sum(hit),
                      paste(df[[id]][hit], collapse = ", "))
    df <- df[!hit, , drop = FALSE]
    if (!is.null(subset)) {
      df <- df[df[[id]] %in% subset, , drop = FALSE]
      say("[%s] subset: %d of %d listed patients present", s, nrow(df), length(subset))
    }
    extra <- setdiff(names(df), schema)
    say("[%s] %d rows x %d columns; not in the output schema: %s", s, nrow(df), length(schema),
        if (length(extra)) paste(extra, collapse = ", ") else "none")
    out[[s]] <- df[, schema, drop = FALSE]
    if (!is.null(P[[s]]$therapy)) {
      th <- P[[s]]$therapy
      out[[paste0("therapy_", s)]] <- th[match(df[[id]], th[[id]]), , drop = FALSE]
    }
  }
  out
}
