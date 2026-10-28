#!/usr/bin/env Rscript
# tools/data_integrity_audit.R
# ==============================================================================
# Data-integrity audit, raw DB -> parser -> excl6 -> the run's inputs (formerly
# diagnostics/diag_53). Read-only on the data; standalone (needs the raw DB).
#
#   [1] raw DB vs standardized files (cell by cell), pre_labelfix vs labelfix
#   [2] standardized -> excl6 -> the run's inputs -> step-01 raw_matrix
#   [3] compositional identities, every family, every row of T0 / T1 / Progressione
#   [4] duplicate / near-duplicate rows and family blocks, within and across sheets
#   [5] plausibility: range, zeros, extreme T0->T1 changes (> 3 robust SD)
#   [6] clinical consistency: best response vs 3/6-month response, PFS, OS, dates, therapy
#   [7] KI67 total in the rows where it disagrees with both decompositions
#
# The raw DB holds names. They are used IN MEMORY ONLY to align rows across sheets;
# nothing written here contains a name, and the output still goes to the git-ignored
# private/ dir. IDs: NSCLC_k, k = Excel row - 1 of T0/T1. Progressione rows get the
# T1 ID of the same full name.
#
# Usage:
#   Rscript tools/data_integrity_audit.R --run <run> --raw <raw DB xlsx>
#           [--allowed-edits CSV] [--excl6 PATTERN] [--out DIR]
#
#   --run            run root (REQUIRED): its step-01 frames and, via run_manifest.yml,
#                    its T0 / T1 input files ("the run's inputs" below)
#   --raw            raw clinical DB (REQUIRED; e.g. data/DB_NSCLC_T0-T1-PD_latest.xlsx)
#   --allowed-edits  CSV (Patient_ID, Timepoint, Marker[, Value]) of the cells deliberately
#                    changed between the excl6 files and the run's inputs (primary run:
#                    data/0final/primary_PIregated_noPS/PI_regated_T1_cells.csv)
#   --excl6          sprintf pattern of the excl6 files
#                    (default data/0final/_upstream_excl6/Dati_NSCLC_excl6_%s.xlsx)
#   --out            output dir (default private/<run name>/data_integrity_audit)
# ==============================================================================
suppressMessages({ library(here); library(readxl); library(stringr) })
setwd(here::here())

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
opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(opt$run) || !nzchar(opt$run)) stop("[data_integrity_audit] --run is required")
if (is.null(opt$raw) || !nzchar(opt$raw)) stop("[data_integrity_audit] --raw is required")
RUN_DIR <- opt$run
if (!dir.exists(RUN_DIR)) stop(sprintf("[data_integrity_audit] run not found: %s", RUN_DIR))
RAW     <- opt$raw
if (!file.exists(RAW)) stop(sprintf("[data_integrity_audit] raw DB not found: %s", RAW))
STD     <- "data/Dati_NSCLC_standardizzati_latest_anonimi_%s.xlsx"
PRE     <- "data/Dati_NSCLC_standardizzati_latest_anonimi_%s_pre_labelfix.xlsx"
EXCL6   <- if (!is.null(opt$excl6)) opt$excl6 else "data/0final/_upstream_excl6/Dati_NSCLC_excl6_%s.xlsx"
OUT_DIR <- if (!is.null(opt$out)) opt$out else file.path("private", basename(normalizePath(RUN_DIR)), "data_integrity_audit")
EXP     <- "BestResponse_2v3_4"
# The run's own inputs, from its manifest.
.man <- yaml::read_yaml(file.path(RUN_DIR, "run_manifest.yml"))
.inp <- vapply(.man$experiments[[EXP]]$inputs, `[[`, "", "path")
RUN_IN <- c(T0 = .inp[grepl("T0\\.xlsx$", .inp)][1], T1 = .inp[grepl("T1\\.xlsx$", .inp)][1])
if (anyNA(RUN_IN) || !all(file.exists(RUN_IN))) stop("[data_integrity_audit] the run's T0/T1 inputs (manifest) are not readable")
RECON   <- function(s) RUN_IN[[s]]
TGT     <- "Best_response_rate_RC_1_RP_2_SD_3_PD_4"
GATE    <- c("KI67NAIVE", "CD28KI67")
EXCL    <- paste0("NSCLC_", c(27, 39, 47, 59, 63, 71))
# Cells deliberately changed between the excl6 files and the run's inputs.
RECON_CELLS <- list()
if (!is.null(opt[["allowed-edits"]])) {
  ae <- read.csv(opt[["allowed-edits"]], stringsAsFactors = FALSE)
  for (s in unique(ae$Timepoint))
    RECON_CELLS[[s]] <- list(ids = unique(ae$Patient_ID[ae$Timepoint == s]), cols = unique(ae$Marker[ae$Timepoint == s]))
}

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
log_con <- file(file.path(OUT_DIR, "diag_53_log.txt"), open = "wt")
say <- function(...) { txt <- sprintf(...); cat(txt, "\n", sep = ""); cat(txt, "\n", sep = "", file = log_con) }
say("diag_53_data_integrity  %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("raw DB: %s | frames from: %s", RAW, RUN_DIR)

# Parser helpers, taken from the parser itself (function definitions only, no side effects)
local({
  ex <- parse("data/database_parser_nsclc_v2.R")
  for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) &&
                    identical(e[[3]][[1]], as.name("function"))) eval(e, globalenv())
})
stopifnot(exists("clean_marker_name"), exists("parse_pdl1_value"), exists("extract_binary_from_raw"),
          exists("extract_n_mets_sites"))

num  <- function(x) suppressWarnings(as.numeric(as.character(x)))
fmt  <- function(x) ifelse(is.na(x), "NA", formatC(signif(x, 5), format = "fg", digits = 5))
nkey <- function(x) tolower(gsub("\\s+", " ", trimws(as.character(x))))

# ── findings accumulator ─────────────────────────────────────────────────────
FIND <- list()
add_finding <- function(id, sheet, excel_row, columns, values, type, check, severity, note = "") {
  FIND[[length(FIND) + 1]] <<- data.frame(ID = id, sheet = sheet, excel_row = excel_row, columns = columns,
    values = values, type = type, check = check, severity = severity, note = note, stringsAsFactors = FALSE)
}

# ── raw sheets, canonical column names, IDs ──────────────────────────────────
rd_raw <- function(s) { d <- as.data.frame(read_excel(RAW, sheet = s, col_types = "text")); names(d) <- trimws(names(d)); d }
RAWS <- list(T0 = rd_raw("T0"), T1 = rd_raw("T1"), PG = rd_raw("Progressione"))
CANON <- c(CD3 = "CD3tot", CD4 = "CD4tot", CD8 = "CD8tot", CD137 = "CD137tot", EM = "EMx", EMRA = "EMx",
           `Treg tot` = "Treg_tot", `Treg active` = "Treg_active", `Treg resting` = "Treg_resting",
           `Treg non suppr` = "Treg_non_suppr")
canon_names <- function(n) { m <- n %in% names(CANON); n[m] <- CANON[n[m]]; n }
MK <- c("CD3tot", "CD4tot", "CD8tot", "CM", "EFF", "EMx", "NAIVE",
        "CD28", "CD28CD4", "CD28CD8", "CD28KI67", "CD28CM", "CD28EFF", "CD28EMRA", "CD28NAIVE",
        "CD137tot", "CD137CD4", "CD137CD8", "CD137KI67", "CD137PD1", "CD137CM", "CD137EFF", "CD137EMRA", "CD137NAIVE",
        "KI67", "KI67CD4", "KI67CD8", "KI67CM", "KI67EFF", "KI657EMRA", "KI67NAIVE",
        "PD1", "PD1CD4", "PD1CD8", "PD1KI67", "PD1CM", "PD1EFF", "PD1EMRA", "PD1NAIVE",
        "Treg_tot", "Treg_active", "Treg_resting", "Treg_non_suppr")
ANALYSED <- setdiff(MK, c("Treg_tot", "Treg_active", "Treg_resting", "Treg_non_suppr"))   # 39 (EMx = EM at T0, EMRA at T1)
NUMS <- lapply(RAWS, function(d) {
  m <- d; names(m) <- canon_names(names(m))
  as.data.frame(lapply(m[, intersect(MK, names(m))], num), check.names = FALSE)
})
k0 <- nkey(paste(RAWS$T0[[1]], RAWS$T0[[2]])); k1 <- nkey(paste(RAWS$T1[[1]], RAWS$T1[[2]]))
kp <- nkey(paste(RAWS$PG[[1]], RAWS$PG[[2]]))
pg_id <- match(kp, k1)
pg_fam <- match(nkey(RAWS$PG[[1]]), nkey(RAWS$T1[[1]]))
pg_id[is.na(pg_id) & !is.na(RAWS$PG[[1]])] <- pg_fam[is.na(pg_id) & !is.na(RAWS$PG[[1]])]
IDS <- list(T0 = paste0("NSCLC_", seq_len(nrow(RAWS$T0))), T1 = paste0("NSCLC_", seq_len(nrow(RAWS$T1))),
            PG = ifelse(is.na(pg_id), NA, paste0("NSCLC_", pg_id)))
EXROW <- lapply(RAWS, function(d) seq_len(nrow(d)) + 1L)

# ── analysis frames (the audited run) ──────────────────────────────────────
DS <- readRDS(Sys.glob(file.path(RUN_DIR, EXP, "01_*", sprintf("data_processed_%s_standard.rds", EXP))))
DL <- readRDS(Sys.glob(file.path(RUN_DIR, EXP, "01_*", sprintf("data_processed_%s_longitudinal.rds", EXP))))
F_T0  <- as.character(DS$metadata$Patient_ID)
F_LMM <- unique(as.character(DL$metadata$Patient_ID))
LMM_OBS <- paste(DL$metadata$Patient_ID, DL$metadata$Timepoint)
F_DEL <- names(which(table(DL$metadata$Patient_ID) == 2))
GRP <- setNames(as.character(DL$metadata$Group), DL$metadata$Patient_ID)
GRP <- c(GRP, setNames(as.character(DS$metadata$Group), DS$metadata$Patient_ID))
GRP <- GRP[!duplicated(names(GRP))]
say("frames: T0 %d | LMM %d pts / %d obs | delta %d (%d RP / %d SD_PD)", length(F_T0), length(F_LMM),
    nrow(DL$metadata), length(F_DEL), sum(GRP[F_DEL] == "RP"), sum(GRP[F_DEL] == "SD_PD"))
br0 <- num(RAWS$T0[["Best response rate (RC=1, RP=2, SD=3, PD=4)"]])
BR3 <- setNames(c(NA, "PR", "SD", "PD")[match(br0, 1:4)], IDS$T0)       # raw 3-level best response (T0 sheet)

in_frames <- function(id, sheet) {
  tp <- c(T0 = "T0", T1 = "T1", PG = NA)[sheet]
  c(T0_76 = sheet == "T0" && id %in% F_T0,
    LMM = !is.na(tp) && paste(id, tp) %in% LMM_OBS,
    delta = !is.na(tp) && id %in% F_DEL)
}
severity_of <- function(id, sheet, cols, base = "media") {
  fr <- in_frames(id, sheet)
  if (!any(fr)) return("bassa")
  if (any(cols %in% GATE) && (fr["LMM"] || fr["delta"])) return("alta")
  if (any(cols %in% ANALYSED)) return(base)
  "bassa"
}

# ==============================================================================
# [1] raw DB vs standardized files
# ==============================================================================
say("\n[1] RAW DB -> STANDARDIZED FILES")
say("rows: raw T0 %d, T1 %d, Progressione %d (not parsed)", nrow(RAWS$T0), nrow(RAWS$T1), nrow(RAWS$PG))
pos_mis <- which(k0 != k1)
fam_mis <- which(nkey(RAWS$T0[[1]]) != nkey(RAWS$T1[[1]]))
say("T0/T1 positional alignment: full name differs at %d positions, family name at %d (%s)",
    length(pos_mis), length(fam_mis), paste(fam_mis, collapse = ", "))
first_na <- pos_mis[is.na(RAWS$T1[pos_mis, 2]) | is.na(RAWS$T0[pos_mis, 2])]
say("  first name blank in one sheet: %s", paste(first_na, collapse = ", "))
dcols <- c(5, 6, 11)   # T0 diagnosis / metastatic / therapy start; T1 = same - 1
same_dates <- function(i) all(mapply(identical, RAWS$T0[i, dcols], RAWS$T1[i, dcols - 1]))
for (i in setdiff(pos_mis, first_na)) {
  sd <- same_dates(i); t1_has <- sum(!is.na(NUMS$T1[i, ])) > 0
  say("  NSCLC_%d: %s differs; diagnosis/mets/start dates %s; T1 markers %s", i,
      if (i %in% fam_mis) "family name" else "first name", if (sd) "identical" else "DIFFER",
      if (t1_has) "present" else "absent")
  if (!sd && i == 64)
    add_finding(sprintf("NSCLC_%d", i), "T0/T1", sprintf("T0 %d / T1 %d", i + 1, i + 1), "Family/First Name",
      "T0 row is a copy of row 63 (NSCLC_62); T1 row is another person", "allineamento T0/T1: persona diversa", "1", "info",
      "noto e gestito: la riga T0 è scartata come duplicato, NSCLC_64 entra solo con T1")
  else if (!sd)
    add_finding(sprintf("NSCLC_%d", i), "T0/T1", sprintf("T0 %d / T1 %d", i + 1, i + 1), "Family/First Name",
      "different person at the same position", "allineamento T0/T1: persona diversa", "1",
      if (t1_has && paste0("NSCLC_", i) %in% F_LMM) "alta" else "bassa",
      "ID posizionali: la riga T1 non è lo stesso paziente della riga T0")
  else if (t1_has)
    add_finding(sprintf("NSCLC_%d", i), "T0/T1", sprintf("T0 %d / T1 %d", i + 1, i + 1), "Family/First Name",
      sprintf("%s name differs, dates identical", if (i %in% fam_mis) "family" else "first"),
      "variante del nome", "1", "info", "stesso paziente (date identiche); grafia del nome da uniformare")
}
fam0 <- nkey(RAWS$T0[[1]])
dupfam <- which(duplicated(fam0) | duplicated(fam0, fromLast = TRUE))
for (g in split(dupfam, fam0[dupfam])) {
  same_person <- length(unique(k0[g])) == 1
  exact_dup <- length(unique(as.character(RAWS$T0[g, 1]))) == 1
  say("  T0 rows sharing a family name: NSCLC_%s -> %s; parser key (exact family string) %s",
      paste(g, collapse = "/"), if (same_person) "same full name" else "DIFFERENT persons",
      if (exact_dup) "collides (later row dropped)" else "does not collide (case/space differs)")
  if (!same_person)
    add_finding(paste0("NSCLC_", g, collapse = "/"), "T0", paste(g + 1, collapse = "/"), "Family Name",
      if (exact_dup) "exact family-name collision" else "same family name, differs only in case/spaces",
      "dedup parser su solo cognome", "1", "bassa",
      if (exact_dup) "il parser scarta il secondo paziente come duplicato" else
        "oggi salvo per caso: il dedup del parser confronta solo la stringa del cognome")
}

# raw -> standardized, column by column (parser renaming rules)
raw_to_parser <- function(d) {
  orig <- names(d)
  resp_idx <- grep("risposta|response", tolower(orig))[1]
  names(d)[resp_idx] <- "Clinical_Response"
  fi <- 3:ncol(d); names(d)[fi] <- sapply(names(d)[fi], clean_marker_name, USE.NAMES = FALSE)
  al <- c(CD4 = "CD4tot", CD8 = "CD8tot", CD137 = "CD137tot")
  for (o in names(al)) if (o %in% names(d) && !al[o] %in% names(d)) names(d)[names(d) == o] <- al[o]
  d
}
cmp_rows <- list()
for (s in c("T0", "T1")) {
  r  <- raw_to_parser(RAWS[[s]])
  st <- as.data.frame(read_excel(sprintf(STD, s)))
  idx <- as.integer(sub("NSCLC_", "", st$Patient_ID))
  date_cols <- setdiff(grep("date|data|follow_up", names(st), ignore.case = TRUE, value = TRUE),
                       "Date_of_birth_at_tumor_metastasi")
  common <- setdiff(intersect(names(st), names(r)), "Patient_ID")
  n_cells <- 0; n_diff <- 0; n_fill <- 0
  for (cn in common) {
    rv <- r[[cn]][idx]; sv <- st[[cn]]
    if (cn %in% date_cols) {
      rv <- ifelse(grepl("^[0-9]+$", rv), format(openxlsx::convertToDate(num(rv)), "%d/%m/%Y"), rv)
      same <- (is.na(rv) & is.na(sv)) | (!is.na(rv) & !is.na(sv) & rv == as.character(sv))
    } else if (cn == "PD_L1") {
      same <- (is.na(parse_pdl1_value(rv)) & is.na(sv)) | abs(parse_pdl1_value(rv) - num(sv)) < 1e-9
    } else if (is.numeric(sv)) {
      same <- (is.na(num(rv)) & is.na(sv)) | abs(num(rv) - sv) < 1e-9
    } else {
      same <- (is.na(rv) & is.na(sv)) | (!is.na(rv) & !is.na(sv) & trimws(rv) == trimws(as.character(sv)))
    }
    same[is.na(same)] <- FALSE
    n_cells <- n_cells + length(same)
    for (i in which(!same)) {
      filled <- is.na(rv[i]) || (is.numeric(sv) && is.na(num(rv[i])) && !is.na(sv[i]))
      if (filled) n_fill <- n_fill + 1 else n_diff <- n_diff + 1
      cmp_rows[[length(cmp_rows) + 1]] <- data.frame(sheet = s, ID = st$Patient_ID[i], column = cn,
        raw = as.character(r[[cn]][idx[i]]), std = as.character(sv[i]),
        kind = if (filled) "raw blank/non-numeric -> std filled" else "DIFFERENT", stringsAsFactors = FALSE)
    }
  }
  mk_common <- intersect(common, ANALYSED)
  say("%s: %d std rows (raw IDs missing: %s); %d common columns (%d analysed markers), %d cells compared: %d identical, %d filled, %d different",
      s, nrow(st), paste(setdiff(seq_len(nrow(RAWS[[s]])), idx), collapse = ", "), length(common),
      length(mk_common), n_cells, n_cells - n_fill - n_diff, n_fill, n_diff)
}
CMP <- do.call(rbind, cmp_rows)
if (!is.null(CMP)) {
  write.csv(CMP, file.path(OUT_DIR, "raw_vs_std_differences.csv"), row.names = FALSE)
  agg <- aggregate(ID ~ sheet + column + kind, CMP, length)
  for (k in seq_len(nrow(agg))) say("  %s %-45s %-38s %d cells", agg$sheet[k], substr(agg$column[k], 1, 45), agg$kind[k], agg$ID[k])
  mk_diff <- CMP[CMP$column %in% c(ANALYSED, "EM", "EMRA") & CMP$kind == "DIFFERENT", ]
  say("  analysed-marker cells that differ raw vs std: %d", nrow(mk_diff))
  for (k in seq_len(nrow(mk_diff)))
    add_finding(mk_diff$ID[k], mk_diff$sheet[k], as.integer(sub("NSCLC_", "", mk_diff$ID[k])) + 1, mk_diff$column[k],
      sprintf("raw %s / std %s", mk_diff$raw[k], mk_diff$std[k]), "DB grezzo != file standardizzato", "1",
      severity_of(mk_diff$ID[k], mk_diff$sheet[k], mk_diff$column[k]))
  br_diff <- CMP[CMP$column == TGT, ]
  for (k in seq_len(nrow(br_diff)))
    add_finding(br_diff$ID[k], br_diff$sheet[k], as.integer(sub("NSCLC_", "", br_diff$ID[k])) + 1, TGT,
      sprintf("raw %s / std %s", br_diff$raw[k], br_diff$std[k]), "best response riempita dal parser", "1",
      if (any(in_frames(br_diff$ID[k], br_diff$sheet[k]))) "media" else "bassa",
      if (br_diff$sheet[k] == "T0") "da file clinico (non più presente in data/): fonte da confermare" else "back-fill T0 -> T1 (fix 2026-09-14)")
}
# derived clinical columns recomputed with the parser's own functions (T0)
r0 <- raw_to_parser(RAWS$T0); s0 <- as.data.frame(read_excel(sprintf(STD, "T0"))); i0 <- as.integer(sub("NSCLC_", "", s0$Patient_ID))
comorb <- extract_binary_from_raw(r0[[grep("comorbid", names(r0), ignore.case = TRUE)[1]]], c("0"))[i0]
toss   <- extract_binary_from_raw(r0[[grep("tossicit", names(r0), ignore.case = TRUE)[1]]], c("0", "no"))[i0]
nmets  <- extract_n_mets_sites(r0[[grep("sites.*of.*mets", names(r0), ignore.case = TRUE)[1]]])[i0]
say("derived (T0): comorbidity %d/%d equal (%d filled where raw-derived NA), toxicity %d/%d (%d filled), n_mets %d/%d",
    sum(comorb == as.character(s0[["comorbidità_No_0_Si_1"]]), na.rm = TRUE), nrow(s0),
    sum(is.na(comorb) & !is.na(s0[["comorbidità_No_0_Si_1"]])),
    sum(toss == as.character(s0[["tossicità_da_immuno_No_0_Si_1"]]), na.rm = TRUE), nrow(s0),
    sum(is.na(toss) & !is.na(s0[["tossicità_da_immuno_No_0_Si_1"]])),
    sum(nmets == s0$n_mets_sites | (is.na(nmets) & is.na(s0$n_mets_sites)), na.rm = TRUE), nrow(s0))

# pre_labelfix vs labelfix
cell_diff <- function(a, b, key = "Patient_ID") {
  a <- as.data.frame(a); b <- as.data.frame(b)
  out <- list()
  ids <- union(a[[key]], b[[key]])
  only_a <- setdiff(a[[key]], b[[key]]); only_b <- setdiff(b[[key]], a[[key]])
  cols <- union(names(a), names(b))
  for (cn in setdiff(intersect(names(a), names(b)), key)) {
    ia <- match(intersect(a[[key]], b[[key]]), a[[key]]); ib <- match(a[[key]][ia], b[[key]])
    va <- a[[cn]][ia]; vb <- b[[cn]][ib]
    eq <- (is.na(va) & is.na(vb)) | (!is.na(va) & !is.na(vb) &
          (if (is.numeric(va) && is.numeric(vb)) abs(va - vb) < 1e-9 else as.character(va) == as.character(vb)))
    eq[is.na(eq)] <- FALSE
    for (i in which(!eq)) out[[length(out) + 1]] <- data.frame(ID = a[[key]][ia[i]], column = cn,
      old = as.character(va[i]), new = as.character(vb[i]), stringsAsFactors = FALSE)
  }
  list(diff = do.call(rbind, out), only_a = only_a, only_b = only_b,
       cols_only_a = setdiff(names(a), names(b)), cols_only_b = setdiff(names(b), names(a)))
}
for (s in c("T0", "T1")) {
  cd <- cell_diff(read_excel(sprintf(PRE, s)), read_excel(sprintf(STD, s)))
  say("pre_labelfix -> labelfix %s: %d cells changed%s; rows only-old %d, only-new %d; cols only-old %d, only-new %d",
      s, if (is.null(cd$diff)) 0 else nrow(cd$diff),
      if (is.null(cd$diff)) "" else sprintf(" (%s)", paste(sprintf("%s %s: %s -> %s", cd$diff$ID, substr(cd$diff$column, 1, 18), cd$diff$old, cd$diff$new), collapse = "; ")),
      length(cd$only_a), length(cd$only_b), length(cd$cols_only_a), length(cd$cols_only_b))
}

# ==============================================================================
# [2] standardized -> excl6 -> the run's inputs -> step-01 raw_matrix
# ==============================================================================
say("\n[2] STANDARDIZED -> EXCL6 -> RUN INPUTS -> RUN")
for (s in c("T0", "T1")) {
  st <- as.data.frame(read_excel(sprintf(STD, s))); ex <- as.data.frame(read_excel(sprintf(EXCL6, s)))
  cd <- cell_diff(st, ex)
  say("std -> excl6 %s: rows %d -> %d, removed {%s}, added {%s}, cells changed %d, cols dropped %d",
      s, nrow(st), nrow(ex), paste(cd$only_a, collapse = ","), paste(cd$only_b, collapse = ","),
      if (is.null(cd$diff)) 0 else nrow(cd$diff), length(cd$cols_only_a))
  stopifnot(setequal(cd$only_a, EXCL), length(cd$only_b) == 0)
  rc <- as.data.frame(read_excel(RECON(s)))
  cd2 <- cell_diff(ex, rc)
  expct <- RECON_CELLS[[s]]
  exp_keys <- if (is.null(expct)) character(0) else as.vector(outer(expct$ids, expct$cols, paste))
  got_keys <- if (is.null(cd2$diff)) character(0) else paste(cd2$diff$ID, cd2$diff$column)
  say("excl6 -> run inputs %s: %d cells changed; expected %d; unexpected %d; expected-but-unchanged %d",
      s, length(got_keys), length(exp_keys), length(setdiff(got_keys, exp_keys)), length(setdiff(exp_keys, got_keys)))
  if (!is.null(cd2$diff)) for (k in seq_len(nrow(cd2$diff)))
    say("    %s %-10s %s -> %s", cd2$diff$ID[k], cd2$diff$column[k], cd2$diff$old[k], cd2$diff$new[k])
}
# the run read exactly its inputs
rm_chk <- function(D, s) {
  rc <- as.data.frame(read_excel(RECON(s)))
  meta <- D$metadata; tp <- if (s == "T0") c("T0", "Baseline") else "T1"
  rows <- which(meta$Timepoint %in% tp)
  if (!length(rows)) return(NA)
  m <- D$raw_matrix[rows, , drop = FALSE]
  r <- match(meta$Patient_ID[rows], rc$Patient_ID)
  cols <- intersect(colnames(m), names(rc))
  max(abs(m[, cols] - as.matrix(rc[r, cols])), na.rm = TRUE)
}
say("step-01 raw_matrix vs the run's inputs, max |diff|: standard T0 %.3g, longitudinal T0 %.3g, T1 %.3g",
    rm_chk(DS, "T0"), rm_chk(DL, "T0"), rm_chk(DL, "T1"))

# ==============================================================================
# [3] compositional identities
# ==============================================================================
say("\n[3] COMPOSITIONAL IDENTITIES (raw DB, every row with marker data)")
FAM <- list(KI67 = c(tot = "KI67", m = c("KI67NAIVE", "KI67CM", "KI67EFF", "KI657EMRA"), c = c("KI67CD4", "KI67CD8")),
            CD137 = c(tot = "CD137tot", m = c("CD137NAIVE", "CD137CM", "CD137EFF", "CD137EMRA"), c = c("CD137CD4", "CD137CD8")),
            PD1 = c(tot = "PD1", m = c("PD1NAIVE", "PD1CM", "PD1EFF", "PD1EMRA"), c = c("PD1CD4", "PD1CD8")),
            CD28 = c(tot = "CD28", m = c("CD28NAIVE", "CD28CM", "CD28EFF", "CD28EMRA"), c = c("CD28CD4", "CD28CD8")))
# tolerances, set from the healthy rows (memory ratio IQR 0.999-1.001; CD4+CD8 runs below the total by the
# CD4-CD8- / CD4+CD8+ cells, 5th pct ~0.85)
ok_abs <- 0.05; ok_rel <- 0.01; rnd_abs <- 0.2; rnd_rel <- 0.05
agree <- function(a, b) abs(a - b) <= pmax(ok_abs, ok_rel * pmax(abs(a), abs(b)))
near  <- function(a, b) abs(a - b) <= pmax(rnd_abs, rnd_rel * pmax(abs(a), abs(b)))
ID_ROWS <- list()
for (s in names(NUMS)) {
  d <- NUMS[[s]]; has <- rowSums(!is.na(d)) > 0
  for (f in names(FAM)) {
    v <- FAM[[f]]; tot <- d[[v["tot"]]]
    mem <- rowSums(d[, v[grep("^m", names(v))]]); cd <- rowSums(d[, v[grep("^c", names(v))]])
    cls <- rep("ok", nrow(d))
    mem_ok <- agree(mem, tot); cd_ok <- cd <= tot * (1 + ok_rel) + ok_abs & cd >= 0.70 * tot
    for (i in which(has & !is.na(tot))) {
      if (is.na(mem[i]) || is.na(cd[i])) { cls[i] <- "NA in sottoinsiemi"; next }
      if (mem_ok[i] && cd_ok[i]) next
      comp_m <- abs(mem[i] - 100) <= 2; comp_c <- abs(cd[i] - 100) <= 3
      cls[i] <- if ((comp_m || comp_c) && tot[i] < 50) "sottoinsiemi in composizione (somma ~100)"
        else if (near(mem[i], tot[i]) && cd_ok[i]) "arrotondamento"
        else if (!mem_ok[i] && agree(mem[i], cd[i])) "totale anomalo (scomposizioni concordi)"
        else if (!mem_ok[i] && !cd_ok[i] && !near(mem[i], cd[i])) "non spiegabile (totale e scomposizioni tutti diversi)"
        else if (!mem_ok[i] && cd_ok[i]) "scomposizione memoria anomala"
        else if (mem_ok[i] && !cd_ok[i] && cd[i] < 0.70 * tot[i]) "CD4+CD8 << totale"
        else if (mem_ok[i] && !cd_ok[i]) "CD4+CD8 > totale"
        else "non spiegabile"
    }
    ID_ROWS[[length(ID_ROWS) + 1]] <- data.frame(sheet = s, ID = IDS[[s]], excel_row = EXROW[[s]], family = f,
      total = tot, sum_memory = mem, sum_cd4cd8 = cd, class = cls, stringsAsFactors = FALSE)[has & !is.na(tot), ]
  }
  # phenotype: NAIVE+CM+EFF+EM(RA) = 100; CD4tot + CD8tot <= ~100 (% of CD3)
  ph <- rowSums(d[, c("NAIVE", "CM", "EFF", "EMx")]); c48 <- d$CD4tot + d$CD8tot
  cls <- ifelse(abs(ph - 100) <= 1, "ok", ifelse(abs(ph - 100) <= 3, "arrotondamento", "non spiegabile"))
  ID_ROWS[[length(ID_ROWS) + 1]] <- data.frame(sheet = s, ID = IDS[[s]], excel_row = EXROW[[s]], family = "MEMORY_PHENOTYPE",
    total = 100, sum_memory = ph, sum_cd4cd8 = NA, class = cls, stringsAsFactors = FALSE)[has & !is.na(ph), ]
  c34 <- d$CD3tot + d$CD4tot
  cls <- ifelse(abs(c34 - 100) < 0.06, "CD4tot = 100 - CD3tot", ifelse(c48 > 101, "CD4+CD8 > 100", ifelse(c48 < 80, "CD4+CD8 < 80", "ok")))
  ID_ROWS[[length(ID_ROWS) + 1]] <- data.frame(sheet = s, ID = IDS[[s]], excel_row = EXROW[[s]], family = "CD4CD8_of_CD3",
    total = d$CD3tot, sum_memory = NA, sum_cd4cd8 = c48, class = cls, stringsAsFactors = FALSE)[has & !is.na(c48), ]
}
IDR <- do.call(rbind, ID_ROWS)
write.csv(IDR, file.path(OUT_DIR, "identity_rows.csv"), row.names = FALSE)
tab <- table(paste(IDR$family), IDR$class)
say("identity classes over all sheets (rows):")
for (f in rownames(tab)) say("  %-17s %s", f, paste(sprintf("%s=%d", colnames(tab), tab[f, ])[tab[f, ] > 0], collapse = " | "))
for (s in names(NUMS)) say("  %s: CD3tot + CD4tot = 100 (+-0.05) in %d / %d rows; CD4tot + CD8tot = 100 in %d", s,
    sum(abs(NUMS[[s]]$CD3tot + NUMS[[s]]$CD4tot - 100) < 0.06, na.rm = TRUE), sum(!is.na(NUMS[[s]]$CD3tot + NUMS[[s]]$CD4tot)),
    sum(abs(NUMS[[s]]$CD4tot + NUMS[[s]]$CD8tot - 100) < 0.06, na.rm = TRUE))
for (s in names(NUMS)) for (f in names(FAM)) {
  x <- IDR[IDR$sheet == s & IDR$family == f & IDR$class == "ok", ]
  say("  %s %-5s healthy rows %d: memory/total median %.4f IQR %.4f-%.4f; (CD4+CD8)/total median %.3f, 5th pct %.3f",
      s, f, nrow(x), median(x$sum_memory / x$total), quantile(x$sum_memory / x$total, .25), quantile(x$sum_memory / x$total, .75),
      median(x$sum_cd4cd8 / x$total), quantile(x$sum_cd4cd8 / x$total, .05))
}
EXC <- IDR[!IDR$class %in% c("ok", "arrotondamento"), ]
fam_cols <- function(f) if (f %in% names(FAM)) unname(FAM[[f]]) else if (f == "MEMORY_PHENOTYPE") c("NAIVE", "CM", "EFF", "EM/EMRA") else c("CD4tot", "CD8tot", "CD3tot")
for (k in seq_len(nrow(EXC))) {
  e <- EXC[k, ]
  if (is.na(e$ID)) next
  cols <- fam_cols(e$family)
  # the columns the class puts under suspicion (they set the severity)
  sus <- if (!e$family %in% names(FAM)) cols
         else if (grepl("^totale anomalo", e$class)) cols[1]
         else if (grepl("CD4\\+CD8", e$class)) cols[6:7]
         else if (grepl("memoria", e$class)) cols[2:5]
         else cols
  handled <- e$sheet %in% names(RECON_CELLS) && e$ID %in% RECON_CELLS[[e$sheet]]$ids && any(cols %in% RECON_CELLS[[e$sheet]]$cols)
  sev <- if (e$sheet == "PG") "bassa" else if (e$class == "NA in sottoinsiemi") "info" else severity_of(e$ID, e$sheet, sus)
  if (e$class == "CD4+CD8 << totale" && sev == "alta") sev <- "media"   # may be CD4-CD8- cells, not an error
  if (e$family == "CD4CD8_of_CD3") cols <- if (grepl("CD3tot", e$class)) c("CD4tot", "CD3tot") else c("CD4tot", "CD8tot")
  add_finding(e$ID, e$sheet, e$excel_row, paste(cols, collapse = ","),
    sprintf("total %s | sum memory %s | CD4+CD8 %s", fmt(e$total), fmt(e$sum_memory), fmt(e$sum_cd4cd8)),
    e$class, sprintf("3:%s", e$family), if (handled) "info" else sev,
    if (handled) "già corretto negli input della run" else "")
}
say("identity exceptions (excluding ok/rounding): %d rows", nrow(EXC))
for (k in seq_len(nrow(EXC))) say("  %-3s %-9s row %3d %-17s total %-8s mem %-8s cd %-8s %s", EXC$sheet[k], EXC$ID[k], EXC$excel_row[k],
                                  EXC$family[k], fmt(EXC$total[k]), fmt(EXC$sum_memory[k]), fmt(EXC$sum_cd4cd8[k]), EXC$class[k])

# double-positive bounds: a cell positive for X AND Y cannot exceed the X+ (or Y+) cells
BOUNDS <- list(c("CD28KI67", "KI67"), c("CD28KI67", "CD28"), c("PD1KI67", "KI67"), c("PD1KI67", "PD1"),
               c("CD137KI67", "KI67"), c("CD137KI67", "CD137tot"), c("CD137PD1", "CD137tot"), c("CD137PD1", "PD1"),
               c("KI67NAIVE", "NAIVE"), c("KI67CM", "CM"), c("KI67EFF", "EFF"), c("KI657EMRA", "EMx"),
               c("KI67CD4", "CD4tot"), c("KI67CD8", "CD8tot"), c("CD28CD4", "CD4tot"), c("CD28CD8", "CD8tot"))
say("double-positive bounds (value <= bound + 0.05):")
BND <- list()
for (b in BOUNDS) {
  viol <- 0; n <- 0
  for (s in names(NUMS)) {
    d <- NUMS[[s]]; x <- d[[b[1]]]; y <- d[[b[2]]]
    # KI67 bound also tested against the KI67 decomposition (the stored KI67 total is itself suspect in some rows)
    ydec <- if (b[2] == "KI67") rowSums(d[, FAM$KI67[2:5]]) else y
    ok <- !is.na(x) & !is.na(y); n <- n + sum(ok)
    w <- which(ok & x > y + 0.05)
    viol <- viol + length(w)
    for (i in w) BND[[length(BND) + 1]] <- data.frame(sheet = s, ID = IDS[[s]][i], excel_row = EXROW[[s]][i],
      value_col = b[1], bound_col = b[2], value = x[i], bound = y[i], bound_decomposition = ydec[i],
      exceeds_decomposition_too = x[i] > ydec[i] + 0.05, stringsAsFactors = FALSE)
  }
  say("  %-9s <= %-8s violated in %3d / %3d rows (%.1f%%)", b[1], b[2], viol, n, 100 * viol / n)
}
BND <- do.call(rbind, BND)
write.csv(BND, file.path(OUT_DIR, "bound_violations.csv"), row.names = FALSE)
C34 <- IDR[IDR$family == "CD4CD8_of_CD3" & IDR$class == "CD4tot = 100 - CD3tot", ]
if (!is.null(BND)) for (k in seq_len(nrow(BND))) {
  e <- BND[k, ]; if (is.na(e$ID)) next
  if (e$bound_col == "CD4tot" && paste(e$sheet, e$ID) %in% paste(C34$sheet, C34$ID)) {
    add_finding(e$ID, e$sheet, e$excel_row, paste(e$value_col, e$bound_col, sep = " vs "),
      sprintf("%s %s > %s %s", e$value_col, fmt(e$value), e$bound_col, fmt(e$bound)),
      "conseguenza di CD4tot = 100 - CD3tot", "3:bounds", "info", "conferma che è CD4tot a essere sbagliato")
    next
  }
  tot_sus <- e$bound_col == "KI67" && !e$exceeds_decomposition_too
  sev <- if (e$sheet == "PG") "bassa" else severity_of(e$ID, e$sheet, if (tot_sus) e$bound_col else e$value_col)
  handled <- e$sheet %in% names(RECON_CELLS) && e$ID %in% RECON_CELLS[[e$sheet]]$ids && e$value_col %in% RECON_CELLS[[e$sheet]]$cols
  add_finding(e$ID, e$sheet, e$excel_row, paste(e$value_col, e$bound_col, sep = " vs "),
    sprintf("%s %s > %s %s (decomposition %s)", e$value_col, fmt(e$value), e$bound_col, fmt(e$bound), fmt(e$bound_decomposition)),
    if (e$bound_col == "KI67" && !e$exceeds_decomposition_too) "doppio positivo > totale (totale sospetto)" else "doppio positivo > singolo positivo",
    "3:bounds", if (handled) "info" else sev, if (handled) "già corretto negli input della run" else "")
}

# ==============================================================================
# [4] duplicates
# ==============================================================================
say("\n[4] DUPLICATES (marker columns, all sheets)")
ALL <- do.call(rbind, lapply(names(NUMS), function(s) {
  d <- NUMS[[s]][, intersect(MK, names(NUMS[[s]]))]
  for (m in setdiff(MK, names(d))) d[[m]] <- NA
  data.frame(sheet = s, ID = IDS[[s]], excel_row = EXROW[[s]], d[, MK], check.names = FALSE, stringsAsFactors = FALSE)
}))
ALL <- ALL[rowSums(!is.na(ALL[, MK])) >= 10, ]
BLOCKS <- list(phenotype = MK[1:7], CD28 = MK[8:15], CD137 = MK[16:24], KI67 = MK[25:31], PD1 = MK[32:39], Treg = MK[40:43])
DUP <- list()
M <- as.matrix(ALL[, MK])
for (a in seq_len(nrow(M) - 1)) for (b in (a + 1):nrow(M)) {
  both <- !is.na(M[a, ]) & !is.na(M[b, ])
  if (sum(both) < 5) next
  eq <- both & abs(M[a, ] - M[b, ]) < 1e-9 & M[a, ] != 0
  frac <- sum(eq) / sum(both)
  blk <- names(BLOCKS)[vapply(BLOCKS, function(cc) { k <- cc[both[cc]]; length(k) >= 4 && all(eq[k]) }, logical(1))]
  if (frac >= 0.5 || length(blk) > 0)
    DUP[[length(DUP) + 1]] <- data.frame(sheet_a = ALL$sheet[a], ID_a = ALL$ID[a], row_a = ALL$excel_row[a],
      sheet_b = ALL$sheet[b], ID_b = ALL$ID[b], row_b = ALL$excel_row[b], n_compared = sum(both),
      n_equal_nonzero = sum(eq), frac_equal = round(frac, 3), identical_blocks = paste(blk, collapse = ","),
      stringsAsFactors = FALSE)
}
DUP <- do.call(rbind, DUP)
write.csv(DUP, file.path(OUT_DIR, "duplicates.csv"), row.names = FALSE)
say("pairs with >= 50%% identical non-zero values or an identical family block: %d", if (is.null(DUP)) 0 else nrow(DUP))
if (!is.null(DUP)) for (k in seq_len(nrow(DUP))) {
  e <- DUP[k, ]
  same_pt <- !is.na(e$ID_a) && !is.na(e$ID_b) && e$ID_a == e$ID_b
  kind <- if (same_pt && setequal(c(e$sheet_a, e$sheet_b), c("T1", "PG"))) "Progressione = T1 (stesso paziente)"
          else if (setequal(c(e$sheet_a, e$sheet_b), c("T1", "PG")) && e$frac_equal >= 0.9) "Progressione = T1 di un ALTRO paziente"
          else if (same_pt) "stesso paziente, T0 = T1" else "pazienti diversi"
  say("  %s %s r%d <-> %s %s r%d: %d/%d equal (%.2f) blocks [%s] -> %s", e$sheet_a, e$ID_a, e$row_a, e$sheet_b, e$ID_b,
      e$row_b, e$n_equal_nonzero, e$n_compared, e$frac_equal, e$identical_blocks, kind)
  ids <- unique(na.omit(c(e$ID_a, e$ID_b)))
  known <- setequal(ids, c("NSCLC_62", "NSCLC_64")) && e$sheet_a == "T0" && e$sheet_b == "T0"
  sev <- if (known) "info"
         else if (kind == "Progressione = T1 (stesso paziente)") "bassa"
         else if (kind == "Progressione = T1 di un ALTRO paziente") "media"
         else if ("PG" %in% c(e$sheet_a, e$sheet_b) || identical(e$identical_blocks, "Treg")) "bassa"
         else if (any(vapply(seq_along(ids), function(j) any(in_frames(ids[j], c(e$sheet_a, e$sheet_b)[j])), logical(1)))) "alta" else "bassa"
  add_finding(paste(c(e$ID_a, e$ID_b), collapse = " / "), paste(e$sheet_a, e$sheet_b, sep = "/"),
    paste(e$row_a, e$row_b, sep = "/"), if (nzchar(e$identical_blocks)) e$identical_blocks else "markers",
    sprintf("%d/%d non-zero values identical", e$n_equal_nonzero, e$n_compared), paste("duplicato:", kind), "4", sev)
}

# can a T1 profile be matched to its own T0? (would expose row shifts/swaps)
fp_mk <- c("CD3tot", "CD8tot", "CM", "EFF", "NAIVE", "CD28", "CD28CD4", "CD28CD8", "CD28CM", "CD28EFF", "CD28EMRA", "CD28NAIVE",
           "PD1", "PD1CD4", "PD1CD8", "PD1CM", "PD1EFF", "PD1EMRA", "PD1NAIVE", "CD137tot", "CD28KI67")
lgt <- function(p) { p <- pmin(pmax(p / 100, 1e-3), 1 - 1e-3); log(p / (1 - p)) }
FA <- sapply(fp_mk, function(m) lgt(NUMS$T0[[m]])); FB <- sapply(fp_mk, function(m) lgt(NUMS$T1[[m]]))
ctr <- apply(rbind(FA, FB), 2, median, na.rm = TRUE); scl <- apply(rbind(FA, FB), 2, mad, na.rm = TRUE)
FA <- sweep(sweep(FA, 2, ctr), 2, scl, "/"); FB <- sweep(sweep(FB, 2, ctr), 2, scl, "/")
pair <- which(rowSums(!is.na(FA)) >= 15 & rowSums(!is.na(FB)) >= 15)
own_rank <- vapply(pair, function(i) {
  dd <- vapply(which(rowSums(!is.na(FA)) >= 15), function(j) sqrt(mean((FB[i, ] - FA[j, ])^2, na.rm = TRUE)), numeric(1))
  match(i, which(rowSums(!is.na(FA)) >= 15)[order(dd)])
}, numeric(1))
say("fingerprint T1 -> T0 (%d markers, %d patients): own T0 nearest in %d, in top 3 in %d -> %s", length(fp_mk), length(pair),
    sum(own_rank == 1), sum(own_rank <= 3), if (sum(own_rank == 1) < length(pair) / 2) "NOT informative (no individual signature); row shifts cannot be detected this way" else "informative")

# ==============================================================================
# [5] plausibility
# ==============================================================================
say("\n[5] PLAUSIBILITY")
for (s in names(NUMS)) {
  d <- NUMS[[s]]
  oor <- which(as.matrix(d) < 0 | as.matrix(d) > 100, arr.ind = TRUE)
  say("%s: out of [0,100]: %d cells%s", s, nrow(oor),
      if (nrow(oor)) paste0(" (", paste(sprintf("%s %s=%s", IDS[[s]][oor[, 1]], names(d)[oor[, 2]], d[oor]), collapse = "; "), ")") else "")
  for (k in seq_len(nrow(oor)))
    add_finding(IDS[[s]][oor[k, 1]], s, EXROW[[s]][oor[k, 1]], names(d)[oor[k, 2]], fmt(d[oor][k]), "fuori range [0,100]", "5",
                severity_of(IDS[[s]][oor[k, 1]], s, names(d)[oor[k, 2]]))
  z <- colSums(d == 0, na.rm = TRUE); z <- z[z > 0]
  say("%s: zeros per column: %s", s, if (length(z)) paste(sprintf("%s=%d", names(z), z), collapse = " ") else "none")
}
# extreme T0 -> T1 change on the pipeline's own scale (step-01 hybrid_data_raw = logit, LOD, no Z)
H <- cbind(DL$metadata[, c("Patient_ID", "Timepoint")], as.data.frame(DL$hybrid_data_raw[, DL$hybrid_markers]))
h0 <- H[H$Timepoint == "T0", ]; h1 <- H[H$Timepoint == "T1", ]
pp <- intersect(h0$Patient_ID, h1$Patient_ID)
DLT <- as.matrix(h1[match(pp, h1$Patient_ID), DL$hybrid_markers]) - as.matrix(h0[match(pp, h0$Patient_ID), DL$hybrid_markers])
rownames(DLT) <- pp
EXT <- list()
REC0 <- as.data.frame(read_excel(RECON("T0"))); REC1 <- as.data.frame(read_excel(RECON("T1")))
for (m in colnames(DLT)) {
  x <- DLT[, m]; s_ <- mad(x, na.rm = TRUE); if (!is.finite(s_) || s_ == 0) s_ <- IQR(x, na.rm = TRUE) / 1.349
  rz <- (x - median(x, na.rm = TRUE)) / s_
  for (i in which(abs(rz) > 3)) {
    id <- pp[i]
    mm <- if (m == "EM" || m == "EMRA") "EMx" else m
    # the run's inputs carry the corrected values; the raw DB still has the wrong ones
    v0 <- REC0[[m]][REC0$Patient_ID == id]; v1 <- REC1[[m]][REC1$Patient_ID == id]
    fam_exc <- IDR[IDR$ID %in% id & IDR$sheet %in% c("T0", "T1") & !IDR$class %in% c("ok", "arrotondamento") &
                   vapply(IDR$family, function(f) mm %in% fam_cols(f) || (f == "KI67" && mm == "KI67"), logical(1)), ]
    lod <- (!is.na(v0) && v0 == 0) || (!is.na(v1) && v1 == 0)
    blk_m <- names(BLOCKS)[vapply(BLOCKS, function(cc) mm %in% cc, logical(1))]
    dup <- !is.null(DUP) && any((DUP$ID_a %in% id | DUP$ID_b %in% id) & DUP$sheet_a != "PG" & DUP$sheet_b != "PG" &
                                !(DUP$ID_a %in% "NSCLC_62" & DUP$ID_b %in% "NSCLC_64") &
                                (DUP$frac_equal >= 0.5 | grepl(blk_m, DUP$identical_blocks, fixed = TRUE)))
    expl <- if (nrow(fam_exc)) sprintf("eccezione di identità (%s %s: %s)", fam_exc$sheet[1], fam_exc$family[1], fam_exc$class[1])
            else if (lod) "zero a un timepoint (LOD: zero -> min/2 prima del logit)"
            else if (dup) "riga duplicata"
            else "nessun errore di dati evidente"
    EXT[[length(EXT) + 1]] <- data.frame(ID = id, group = unname(GRP[id]), br3 = unname(BR3[id]), marker = m, T0 = v0, T1 = v1,
      delta_logit = round(x[i], 3), robust_z = round(rz[i], 2), explanation = expl, stringsAsFactors = FALSE)
  }
}
EXT <- do.call(rbind, EXT)
write.csv(EXT, file.path(OUT_DIR, "extreme_changes.csv"), row.names = FALSE)
say("paired patients %d x %d markers: %d changes with |robust z| > 3 (%.1f%%) in %d patients",
    nrow(DLT), ncol(DLT), nrow(EXT), 100 * nrow(EXT) / sum(!is.na(DLT)), length(unique(EXT$ID)))
say("  by explanation: %s", paste(sprintf("%s=%d", names(table(EXT$explanation)), table(EXT$explanation)), collapse = " | "))
say("  by group: %s", paste(sprintf("%s=%d", names(table(EXT$group)), table(EXT$group)), collapse = " "))
g <- EXT[EXT$marker %in% GATE, ]
say("  gate markers: %d extreme changes", nrow(g))
for (k in seq_len(nrow(g))) say("    %s (%s/%s) %s T0 %s -> T1 %s, delta %.2f, z %.1f: %s", g$ID[k], g$group[k], g$br3[k], g$marker[k],
                                fmt(g$T0[k]), fmt(g$T1[k]), g$delta_logit[k], g$robust_z[k], g$explanation[k])
pt <- sort(table(EXT$ID), decreasing = TRUE)
say("  patients with >= 4 extreme changes: %s", paste(sprintf("%s(%s,%d)", names(pt)[pt >= 4], GRP[names(pt)[pt >= 4]], pt[pt >= 4]), collapse = " "))
for (id in names(pt)[pt >= 4]) {
  e <- EXT[EXT$ID == id, ]
  add_finding(id, "T0/T1", sprintf("T0 %d / T1 %d", as.integer(sub("NSCLC_", "", id)) + 1, as.integer(sub("NSCLC_", "", id)) + 1),
    paste(e$marker, collapse = ","), sprintf("%d markers with |robust z| > 3 (max %.1f)", nrow(e), max(abs(e$robust_z))),
    "variazione T0->T1 estrema su molti marcatori", "5", if (any(e$marker %in% GATE)) "alta" else "media",
    paste(unique(e$explanation), collapse = "; "))
}
for (k in seq_len(nrow(g))) if (!g$ID[k] %in% names(pt)[pt >= 4])
  add_finding(g$ID[k], "T0/T1", sprintf("T0 %d / T1 %d", as.integer(sub("NSCLC_", "", g$ID[k])) + 1, as.integer(sub("NSCLC_", "", g$ID[k])) + 1),
    g$marker[k], sprintf("T0 %s -> T1 %s, z %.1f", fmt(g$T0[k]), fmt(g$T1[k]), g$robust_z[k]),
    "variazione T0->T1 estrema (marcatore del gate)", "5", "alta", g$explanation[k])

# ==============================================================================
# [6] clinical consistency
# ==============================================================================
say("\n[6] CLINICAL CONSISTENCY (T0 sheet = patient-level reference)")
C0 <- RAWS$T0
cc <- function(pat, d = C0) d[[grep(pat, names(d), fixed = TRUE)[1]]]
CL <- data.frame(ID = IDS$T0, excel_row = EXROW$T0, br = br0, br_raw = cc("Best response rate", C0),
  r3 = cc("(3 mesi)"), r6 = cc("(6 mesi)"), pfs = num(cc("PFS")), os = num(cc("OS")), dead = num(cc("Alive=0")),
  ther = cc("Terapia"), d_diag = num(cc("Date of  tumor")), d_mets = num(cc("Date of metastatic")),
  d_start = num(cc("Data inizio")), d_end = num(cc("Data fine")), d_fu = num(cc("ultimo follow")), stringsAsFactors = FALSE)
CL <- CL[!is.na(C0[[1]]), ]
fr <- function(id) any(c(id %in% F_T0, id %in% F_LMM))
say("best response x 3-month response:"); print_tab <- function(t) for (r in seq_len(nrow(t))) say("   %-16s %s", rownames(t)[r], paste(sprintf("%s=%d", colnames(t), t[r, ]), collapse = "  "))
print_tab(table(CL$br_raw, CL$r3, useNA = "ifany")); say("best response x 6-month response:"); print_tab(table(CL$br_raw, CL$r6, useNA = "ifany"))
say("PFS (months) by best response: %s", paste(tapply(CL$pfs, CL$br, function(x) sprintf("[%s]", paste(quantile(x, c(0, .5, 1), na.rm = TRUE), collapse = "/"))), collapse = " "))
mon <- function(a, b) (b - a) / 30.4375
CL$os_dates <- mon(CL$d_start, CL$d_fu); CL$ther_dur <- mon(CL$d_start, CL$d_end)
clin_flag <- function(i, cols, vals, type, note = "", sev_in = "media") {
  id <- CL$ID[i]
  add_finding(id, "T0", CL$excel_row[i], cols, vals, type, "6", if (fr(id)) sev_in else "bassa", note)
}
for (i in seq_len(nrow(CL))) {
  b <- CL$br[i]
  if (is.na(b) && !is.na(CL$br_raw[i])) clin_flag(i, "Best response", CL$br_raw[i], "codice best response non numerico",
                                                  "escluso dal parser/pipeline (outcome_eligibility)")
  if (!is.na(b) && b == 2 && !is.na(CL$r3[i]) && CL$r3[i] == "NR" && !is.na(CL$r6[i]) && CL$r6[i] == "NR")
    clin_flag(i, "Best response, Risposta 3/6 mesi", sprintf("BR 2, 3m NR, 6m NR, PFS %s", fmt(CL$pfs[i])), "PR ma NR a 3 e 6 mesi", "significato di R/NR da chiarire")
  if (!is.na(b) && b == 4 && !is.na(CL$r3[i]) && CL$r3[i] == "R")
    clin_flag(i, "Best response, Risposta 3 mesi", sprintf("BR 4, 3m R, 6m %s, PFS %s", CL$r6[i], fmt(CL$pfs[i])), "PD ma R a 3 mesi", "significato di R/NR da chiarire")
  if (!is.na(b) && b == 2 && !is.na(CL$pfs[i]) && CL$pfs[i] < 3)
    clin_flag(i, "Best response, PFS", sprintf("BR 2, PFS %s", fmt(CL$pfs[i])), "PR con PFS < 3 mesi", sev_in = "alta")
  if (!is.na(b) && b == 4 && !is.na(CL$pfs[i]) && CL$pfs[i] > 6)
    clin_flag(i, "Best response, PFS", sprintf("BR 4, PFS %s", fmt(CL$pfs[i])), "PD con PFS > 6 mesi", sev_in = "alta")
  if (!is.na(CL$pfs[i]) && !is.na(CL$os[i]) && CL$pfs[i] > CL$os[i] + 0.5)
    clin_flag(i, "PFS, OS", sprintf("PFS %s > OS %s", fmt(CL$pfs[i]), fmt(CL$os[i])), "PFS > OS")
  if (!is.na(CL$d_start[i]) && !is.na(CL$d_end[i]) && CL$d_end[i] < CL$d_start[i])
    clin_flag(i, "Data inizio/fine terapia", sprintf("fine - inizio = %d giorni", CL$d_end[i] - CL$d_start[i]), "fine terapia prima dell'inizio")
  if (!is.na(CL$d_diag[i]) && !is.na(CL$d_mets[i]) && CL$d_mets[i] < CL$d_diag[i])
    clin_flag(i, "Date diagnosi/metastasi", sprintf("metastasi - diagnosi = %d giorni", CL$d_mets[i] - CL$d_diag[i]), "metastasi prima della diagnosi")
  if (!is.na(CL$d_start[i]) && !is.na(CL$d_mets[i]) && CL$d_start[i] < CL$d_mets[i])
    clin_flag(i, "Date metastasi/inizio terapia", sprintf("inizio - metastasi = %d giorni", CL$d_start[i] - CL$d_mets[i]), "terapia iniziata prima della data di metastasi")
  if (!is.na(CL$os[i]) && !is.na(CL$os_dates[i]) && abs(CL$os[i] - CL$os_dates[i]) > 3)
    clin_flag(i, "OS, Data inizio, ultimo follow up", sprintf("OS %s vs follow-up - inizio %.1f mesi", fmt(CL$os[i]), CL$os_dates[i]), "OS non coerente con le date")
  if (!is.na(CL$pfs[i]) && !is.na(CL$os_dates[i]) && CL$pfs[i] > CL$os_dates[i] + 3)
    clin_flag(i, "PFS, date", sprintf("PFS %s > follow-up - inizio %.1f", fmt(CL$pfs[i]), CL$os_dates[i]), "PFS oltre l'ultimo follow-up")
  th <- CL$ther[i]
  if (!is.na(th) && !th %in% c("1", "2", "3", "4"))
    clin_flag(i, "Terapia", th, "codice terapia fuori legenda (1-4)",
              "entra come covariata numerica (sensitivity_covariates) se numerico; NA se testo")
}
say("therapy codes (T0 sheet): %s", paste(sprintf("'%s'=%d", names(table(CL$ther, useNA = "ifany")), table(CL$ther, useNA = "ifany")), collapse = " "))
say("therapy duration vs PFS (months), PD only: median duration %.1f, PFS %.1f; patients on therapy > PFS + 6 months: %d",
    median(CL$ther_dur[CL$br %in% 4], na.rm = TRUE), median(CL$pfs[CL$br %in% 4], na.rm = TRUE),
    sum(CL$ther_dur > CL$pfs + 6, na.rm = TRUE))
# patient-level values: T0 sheet vs T1 sheet (aligned positions) and Progressione vs T1
PL <- c("Sex", "PD-L1", "Terapia", "Data inizio", "PFS", "OS", "Best response", "(3 mesi)", "(6 mesi)", "PS prima")
plv <- function(d, p) as.character(d[[grep(p, names(d), fixed = TRUE)[1]]])
same_pt <- which(nkey(RAWS$T0[[1]]) == nkey(RAWS$T1[[1]]) | !(seq_len(nrow(RAWS$T0)) %in% fam_mis))
same_pt <- setdiff(same_pt, 64)
disc <- list()
for (p in PL) {
  a <- plv(RAWS$T0, p)[same_pt]; b <- plv(RAWS$T1, p)[same_pt]
  w <- which(!is.na(a) & !is.na(b) & trimws(a) != trimws(b))
  if (length(w)) disc[[p]] <- sprintf("NSCLC_%d %s/%s", same_pt[w], a[w], b[w])
  for (j in w) add_finding(sprintf("NSCLC_%d", same_pt[j]), "T0/T1", sprintf("T0 %d / T1 %d", same_pt[j] + 1, same_pt[j] + 1), p,
    sprintf("T0 %s / T1 %s", a[j], b[j]), "valore paziente diverso tra fogli T0 e T1", "6",
    if (fr(sprintf("NSCLC_%d", same_pt[j])) && p %in% c("Best response", "PD-L1", "Terapia")) "media" else "bassa",
    "la pipeline legge PD-L1, terapia, PFS e OS dal file T0; il foglio T1 conta solo per la best response di T1")
}
say("T0 vs T1 patient-level disagreements: %s", if (length(disc)) paste(sprintf("%s: %s", names(disc), sapply(disc, paste, collapse = ", ")), collapse = " | ") else "none")
pgm <- which(!is.na(pg_id))
for (p in PL) {
  a <- plv(RAWS$PG, p)[pgm]; b <- plv(RAWS$T1, p)[pg_id[pgm]]
  w <- which(!is.na(a) & !is.na(b) & trimws(a) != trimws(b))
  if (length(w)) say("Progressione vs T1 %s: %s", p, paste(sprintf("%s %s/%s", IDS$PG[pgm][w], a[w], b[w]), collapse = ", "))
}

# ==============================================================================
# [7] KI67 total in the rows that disagree with both decompositions
# ==============================================================================
say("\n[7] TOTALS that disagree with both decompositions (KI67 first, then the other families)")
K <- IDR[IDR$family %in% names(FAM) & grepl("^totale anomalo|^non spiegabile \\(totale", IDR$class), ]
K <- K[order(K$family != "KI67", K$family, K$sheet), ]
say("%d rows: %s", nrow(K), paste(sprintf("%s %s", K$sheet, K$ID), collapse = ", "))
cand_par <- c("CD3tot", "CD4tot", "CD8tot", "CD28", "PD1", "CD137tot", "NAIVE", "CM", "EFF", "EMx")
K7 <- list()
for (k in seq_len(nrow(K))) {
  s <- K$sheet[k]; i <- match(K$excel_row[k], EXROW[[s]]); d <- NUMS[[s]]
  tot <- K$total[k]; dec <- K$sum_memory[k]; fam <- K$family[k]; tcol <- FAM[[fam]][["tot"]]
  # (a) does the stored total equal another column of the same row?
  rawrow <- vapply(RAWS[[s]][i, ], num, numeric(1))
  eq_col <- names(rawrow)[!is.na(rawrow) & abs(rawrow - tot) < 0.006 & !names(rawrow) %in% c(tcol, "CD137")]
  # (b) does it equal the decomposition (or the total) of another row / sheet? (row shift or swap)
  dec_all <- do.call(rbind, lapply(names(NUMS), function(ss) data.frame(sheet = ss, ID = IDS[[ss]], tot = NUMS[[ss]][[tcol]],
                     dec = rowSums(NUMS[[ss]][, FAM[[fam]][2:5]]), stringsAsFactors = FALSE)))
  m_dec <- dec_all[!is.na(dec_all$dec) & abs(dec_all$dec - tot) < 0.006, ]
  m_tot <- dec_all[!is.na(dec_all$tot) & abs(dec_all$tot - dec) < 0.006 & !(dec_all$sheet == s & dec_all$ID == K$ID[k]), ]
  # (c) is the ratio total/decomposition a parent-population ratio?
  ratio <- tot / dec
  pr <- vapply(cand_par, function(p) d[[p]][i] / 100, numeric(1))
  best <- names(pr)[which.min(abs(log(pmax(pr, 1e-6) / ratio)))]
  best_inv <- names(pr)[which.min(abs(log(pmax(1 / pr, 1e-6) / ratio)))]
  other <- if (s == "T0") NUMS$T1[[tcol]][i] else if (s == "T1") NUMS$T0[[tcol]][i] else NA
  decimal <- abs(ratio - 10) < 0.6 || abs(ratio - 0.1) < 0.006
  K7[[k]] <- data.frame(family = fam, sheet = s, ID = K$ID[k], excel_row = K$excel_row[k], total = tot, decomposition_memory = dec,
    decomposition_cd4cd8 = K$sum_cd4cd8[k], ratio_total_over_dec = round(ratio, 3),
    same_value_in_columns = paste(eq_col, collapse = ","),
    equals_decomposition_of = paste(sprintf("%s %s", m_dec$sheet, m_dec$ID), collapse = ","),
    decomposition_equals_total_of = paste(sprintf("%s %s", m_tot$sheet, m_tot$ID), collapse = ","),
    closest_parent_ratio = sprintf("%s/100 = %.3f", best, pr[best]), closest_inverse = sprintf("100/%s = %.3f", best_inv, 1 / pr[best_inv]),
    decimal_point_error = decimal, total_other_timepoint = other, CD28KI67 = d$CD28KI67[i], PD1KI67 = d$PD1KI67[i], CD137KI67 = d$CD137KI67[i],
    group = unname(GRP[K$ID[k]]), stringsAsFactors = FALSE)
  say("  %-5s %-3s %-9s row %3d total %-6s dec %-6s/%-6s ratio %-6.3f%s | same-row equal: {%s} | = dec of {%s} | dec = total of {%s} | %s, %s | other tp %s | CD28KI67 %s",
      fam, s, K$ID[k], K$excel_row[k], fmt(tot), fmt(dec), fmt(K$sum_cd4cd8[k]), ratio, if (decimal) " DECIMAL?" else "", paste(eq_col, collapse = ","),
      K7[[k]]$equals_decomposition_of, K7[[k]]$decomposition_equals_total_of, K7[[k]]$closest_parent_ratio,
      K7[[k]]$closest_inverse, fmt(other), fmt(d$CD28KI67[i]))
}
K7 <- do.call(rbind, K7)
write.csv(K7, file.path(OUT_DIR, "ki67_total_rows.csv"), row.names = FALSE)
kk <- K7[K7$family == "KI67", ]
if (nrow(kk)) say("KI67: ratio total/decomposition median %.3f, range %.3f-%.3f (%d of %d below 1); correlation with CD3tot/100 r = %.2f; CD28KI67 > total in %d, > decomposition in %d",
    median(kk$ratio_total_over_dec), min(kk$ratio_total_over_dec), max(kk$ratio_total_over_dec), sum(kk$ratio_total_over_dec < 1), nrow(kk),
    suppressWarnings(cor(kk$ratio_total_over_dec, mapply(function(s, r) NUMS[[s]]$CD3tot[match(r, EXROW[[s]])] / 100, kk$sheet, kk$excel_row))),
    sum(kk$CD28KI67 > kk$total + 0.05, na.rm = TRUE), sum(kk$CD28KI67 > kk$decomposition_memory + 0.05, na.rm = TRUE))

# ==============================================================================
# findings table
# ==============================================================================
FD <- do.call(rbind, FIND)
fr_cols <- t(mapply(function(id, sh) {
  ids <- trimws(strsplit(id, "/")[[1]]); shs <- trimws(strsplit(sh, "/")[[1]])
  if (length(shs) == 1) shs <- rep(shs, length(ids))
  if (length(ids) == 1 && length(shs) > 1) ids <- rep(ids, length(shs))
  ff <- mapply(function(a, b) if (b == "PG" || is.na(a)) c(FALSE, FALSE, FALSE) else in_frames(a, b), ids, shs)
  apply(matrix(ff, nrow = 3), 1, any)
}, FD$ID, FD$sheet))
FD$in_T0_76 <- fr_cols[, 1]; FD$in_LMM <- fr_cols[, 2]; FD$in_delta <- fr_cols[, 3]
FD$gate_marker <- vapply(FD$columns, function(x) any(vapply(GATE, grepl, logical(1), x = x, fixed = TRUE)), logical(1))
FD$group <- vapply(FD$ID, function(x) paste(unique(na.omit(GRP[trimws(strsplit(x, "/")[[1]])])), collapse = "/"), character(1))
FD <- FD[order(match(FD$severity, c("alta", "media", "bassa", "info")), FD$check, FD$ID), ]
write.csv(FD, file.path(OUT_DIR, "findings.csv"), row.names = FALSE)
say("\nFINDINGS: %d rows -> %s", nrow(FD), file.path(OUT_DIR, "findings.csv"))
tb <- table(FD$type, FD$severity)
for (r in rownames(tb)) say("  %-48s %s", substr(r, 1, 48), paste(sprintf("%s=%d", colnames(tb), tb[r, ])[tb[r, ] > 0], collapse = " "))
close(log_con)
