#!/usr/bin/env Rscript
# tests/test_build_inputs.R
# ==============================================================================
# Unit test of the step-00 data build (R/modules_data_build.R) on a tiny SYNTHETIC
# curated workbook (placeholder names, invented values). Each row exercises one rule:
# header cleaning and aliases, dates, PD-L1 parsing, numeric coercion, derived
# flags (with a supplied fallback column), metastatic-site count, the missing-value
# row filter, duplicate records, the same-person back-fill, exclusions, the subset
# option and the therapy tables.
#
# Run:  Rscript tests/test_build_inputs.R      (from the repository root, ~5 s)
# Exit: 0 on pass, 1 on failure.
# ==============================================================================
suppressPackageStartupMessages({ library(openxlsx); library(readxl) })
source("R/modules_data_build.R")

work <- tempfile("build_inputs_"); dir.create(work)
src  <- file.path(work, "curated.xlsx")

# --- Synthetic curated workbook ------------------------------------------------
# T0: P1 normal · P2 non-numeric marker cell, first name blank on T1 · P3 excluded ·
#     P4 mostly empty (row filter) · P5 duplicate record of P1 (names differ only in
#     case / spaces) · P6 undetermined comorbidity filled by the supplied column.
#     The last row has no Patient_ID (ignored).
T0 <- data.frame(
  `Family Name`      = c("Fam01", "Fam02", "Fam03", "Fam04", "  FAM01 ", "Fam06", NA),
  `First Name`       = c("Giv01", "Giv02", "Giv03", "Giv04", "giv01", "Giv06", NA),
  Patient_ID         = c("P1", "P2", "P3", "P4", "P5", "P6", NA),
  `Diag date`        = c(45000, 44000, 44500, NA, 45000, 43000, NA),
  PD_L1              = c("50", "<1", "1-49", ">50", "50", "1/49", NA),
  `Best response`    = c("2", "3", "4", NA, "2", "4", NA),
  CD4                = c("40", "35", "30", NA, "40", "20", NA),
  `KI67 %`           = c("10", "abc", "12", NA, "10", "8", NA),
  `Comorbidity (text)` = c("no", "COPD", "0", NA, "no", "na", NA),
  comorb_flag        = c(NA, NA, NA, NA, NA, "1", NA),
  `Sites of mets`    = c("1, 6", "2", "0", NA, "1, 6", "1, 1, 12", NA),
  `Therapy code`     = c("1", "2", "3", "4", "1", "2", NA),
  Excluded           = c("0", "0", "1", "0", "0", "0", NA),
  Exclusion_reason   = c(NA, NA, "synthetic reason", NA, NA, NA, NA),
  check.names = FALSE, stringsAsFactors = FALSE)
# T1: P1 best response blank (filled from T0) · P2 first name blank (still the same
#     person: filled) · P3 excluded · P6 row belongs to another person (not filled).
T1 <- data.frame(
  `Family Name`      = c("Fam01", "Fam02", "Fam03", "Fam09"),
  `First Name`       = c("Giv01", NA, "Giv03", "Giv09"),
  Patient_ID         = c("P1", "P2", "P3", "P6"),
  PD_L1              = c("2: 30%", "<1", "1-49", "1/49"),
  `Best response`    = c(NA, NA, "4", NA),
  CD4tot             = c("41", "36", "31", "21"),
  `KI67 %`           = c("5", "9", "11", "7"),
  `Comorbidity (text)` = c("no", "COPD", "0", "no"),
  `Sites of mets`    = c("1, 6", "2", "0", "3"),
  `Therapy code`     = c("1", "2", "3", "2"),
  Excluded           = c("0", "0", "1", "0"),
  Exclusion_reason   = c(NA, NA, "synthetic reason", NA),
  check.names = FALSE, stringsAsFactors = FALSE)
wb <- createWorkbook()
for (s in c("T0", "T1")) { addWorksheet(wb, s); writeData(wb, s, get(s)) }
saveWorkbook(wb, src)

cfg <- list(
  source = src, sheets = list("T0", "T1"),
  columns = list(id = "Patient_ID", person = list("Family Name", "First Name"),
                 exclusion = "Excluded", exclusion_reason = "Exclusion_reason",
                 rename = list(), aliases = list(CD4 = "CD4tot"),
                 dates = list("Diag_date"), text = list("Best_response", "^Sites_of_mets"),
                 pdl1 = "PD_L1"),
  derived = list(list(name = "comorb_flag", from = "^Comorbidity", zero = list("0", "no"))),
  mets_count = list(name = "n_mets", from = "^Sites_of_mets", min_code = 0, max_code = 10),
  patient_level = list("Best_response"),
  qc = list(max_na_row = 0.5),
  therapy = list(T0 = list(`Therapy code` = "th_code"), T1 = list(`Therapy code` = "th_code")),
  output_columns = list(
    T0 = list("Patient_ID", "Diag_date", "PD_L1", "Best_response", "CD4tot", "KI67", "comorb_flag", "n_mets"),
    T1 = list("Patient_ID", "PD_L1", "Best_response", "CD4tot", "KI67")))

quiet <- function(fmt, ...) invisible(NULL)
res <- db_build_inputs(cfg, say = quiet)
sub <- db_build_inputs(cfg, subset = c("P1", "P6"), say = quiet)

# --- Checks ----------------------------------------------------------------------
fails <- 0L
check <- function(label, ok) {
  ok <- isTRUE(ok)
  cat(sprintf("[%s] %s\n", if (ok) "PASS" else "FAIL", label))
  if (!ok) fails <<- fails + 1L
}
t0 <- res$T0; t1 <- res$T1
row <- function(df, id) df[df$Patient_ID == id, , drop = FALSE]

check("outputs: T0, T1 and the two therapy tables",
      setequal(names(res), c("T0", "T1", "therapy_T0", "therapy_T1")))
check("T0 rows: excluded, mostly-empty and duplicate records removed; ID-less row ignored",
      identical(t0$Patient_ID, c("P1", "P2", "P6")))
check("T1 rows: excluded patient removed", identical(t1$Patient_ID, c("P1", "P2", "P6")))
check("output schema and column order", identical(names(t0), unlist(cfg$output_columns$T0)) &&
        identical(names(t1), unlist(cfg$output_columns$T1)))
check("no person column in any output",
      !any(unlist(lapply(res, names)) %in% c("Family Name", "First Name", "Family_Name", "First_Name")))
check("header cleaning ('KI67 %' -> KI67) and alias (CD4 -> CD4tot)",
      identical(t0$CD4tot, c(40, 35, 20)) && identical(t1$KI67, c(5, 9, 7)))
check("Excel serial date -> dd/mm/yyyy", identical(row(t0, "P1")$Diag_date, "15/03/2023"))
check("PD-L1: plain, <X -> X/2, X/Y -> midpoint",
      identical(t0$PD_L1, c(50, 0.5, 25)))
check("PD-L1: 'N:' category prefix and % stripped", identical(row(t1, "P1")$PD_L1, 30))
check("non-numeric marker cell -> missing", is.na(row(t0, "P2")$KI67))
check("derived flag: 'no' -> 0, other text -> 1, 'na' -> supplied column",
      identical(t0$comorb_flag, c("0", "1", "1")))
check("metastatic sites: distinct codes in range", identical(t0$n_mets, c(2L, 1L, 1L)))
check("back-fill from T0 for the same person (blank first name is not a mismatch)",
      identical(row(t1, "P1")$Best_response, "2") && identical(row(t1, "P2")$Best_response, "3"))
check("no back-fill when the T1 row is another person", is.na(row(t1, "P6")$Best_response))
check("therapy table aligned with the analysis rows",
      identical(res$therapy_T0$Patient_ID, t0$Patient_ID) && identical(res$therapy_T0$th_code, c(1, 2, 2)))
check("subset keeps only the listed patients",
      identical(sub$T0$Patient_ID, c("P1", "P6")) && identical(sub$T1$Patient_ID, c("P1", "P6")))

unlink(work, recursive = TRUE)
if (fails > 0L) { cat(sprintf("\n%d check(s) failed\n", fails)); quit(status = 1) }
cat("\nAll checks passed\n")
