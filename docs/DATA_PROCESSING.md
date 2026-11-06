# Data processing

How the study database becomes the analysis input files, and which data rules the pipeline applies after that. This
document lists rules and categories only; patient-level details (identifiers, values, reasons) are in the data
package, which is not public.

## 1. Study database

One Excel workbook with one row per patient and sample time: sheet `T0` (before the first ICI cycle) and sheet `T1`
(after the first cycle). Columns: demographics, clinical history, PD-L1 TPS, therapy, best overall response,
survival, haematology, and the PBMC flow-cytometry panel (percentages of the parent population). A third sheet
(`Progressione`) is not analysed.

## 2. Curated database (`data/curated/DB_NSCLC_curated.xlsx`)

The curated database is the study database with the information the analysis needs and the study database does not
hold. Every cell not listed in its change log is identical to the study database. Each change is logged with sheet,
cell, old value, new value, source and rule. Categories:

| Category | What it adds |
|---|---|
| Pseudonymous identifier | `Patient_ID` on both sheets (assigned by row order; a patient's T0 and T1 rows are checked to belong to the same person before any value is copied between them) |
| Clinical supplement | Best overall response where the database cell is blank; toxicity and comorbidity flags; T0 haematology (total WBC, neutrophils, lymphocytes, NLR) and NOX2. Values are never overwritten |
| Haematology extract | Missing haematology values of a few patients, from a later extract of the same database (never overwrites) |
| Re-gated samples | The Ki67-subset cells of three T1 samples, re-gated by the investigators after the original gating had used the wrong parent population. The original gating files were not available |
| Investigator decisions | Best overall response recorded as clinical progression without radiological assessment, coded as progressive disease |
| Therapy coding | Therapy code, regimen class (ICI monotherapy / chemo-immunotherapy) and a radiotherapy annotation, from the free-text therapy column by a fixed rule (below) |
| Exclusions | `Excluded` (1 / 0) and `Exclusion_reason`, set by investigator decision |

### Therapy coding rule

1. The legend in the column header applies first: 1 pembrolizumab; 2 pembrolizumab + chemotherapy;
   3 ipilimumab + nivolumab + chemotherapy; 4 atezolizumab + chemotherapy. A code followed by a note keeps its
   meaning; a note about radiotherapy between T0 and T1 sets the radiotherapy annotation.
2. Codes outside the header legend are read with the investigators' detailed legend (pembrolizumab alone → 1;
   nivolumab + ipilimumab + carboplatin + pemetrexed → 3).
3. Free text is read literally (drug names → code where the header legend has one; otherwise only the regimen class).
4. Empty cell: no regimen, unless the investigators supplied it.

Regimen class: 1 = ICI monotherapy (any checkpoint inhibitor alone), 2 = chemo-immunotherapy. The radiotherapy
annotation is 1 when written in the cell and 0 otherwise (0 is not a confirmed absence).

## 3. Step 00: analysis inputs (`src/00_build_inputs.R`)

Rules are in `config/data_build.yml`; code in `R/modules_data_build.R`. Per sheet, the baseline sheet first:

1. Read every cell as text; rows without a `Patient_ID` are ignored.
2. Person columns (names), exclusion columns and curated therapy columns are split off; names are used only to match
   records and are never written.
3. Headers: explicit renames, then cleaning (`%` dropped; spaces and `+ - ( , ) =` → `_`), then aliases (the same
   marker named differently on the two sheets).
4. Dates: Excel serial numbers → `dd/mm/yyyy`.
5. Derived binary variables from free text (comorbidity, immune-related toxicity): `0` / `no` → 0, blank or `na` →
   missing, any other text → 1. Where the free text leaves a cell undetermined, a supplied column of the same name
   fills it.
6. Number of metastatic sites = number of distinct site codes (0–10) in the free-text cell.
7. PD-L1 TPS free text → percentage: plain number as is; range `X-Y` or `X/Y` → midpoint; `<X` → X/2;
   `>X` → (X + 100)/2; anything else → missing.
8. Every other non-text column is coerced to numeric; text cells that become missing are counted in the build log.
9. Rows with more than 50% missing values over the output columns are removed.
10. Duplicate records of the same person on a sheet: the first is kept.
11. Patient-level attributes (best overall response): a blank T1 cell takes the T0 value of the same patient, only if
    both rows belong to the same person.
12. Rows flagged `Excluded = 1` are removed; with `--subset`, only the listed patients are kept.
13. The output columns are written in the order of the config (`output_columns`, the data dictionary), with a build
    log and checksums. The curated therapy columns are written to separate tables used only by step 07.

## 4. Rules inside the pipeline

Set in the run config (`config/0final.yml`; `config/nlr_v2_cc50.yml` for the sensitivity run).

- **Outcome.** Best overall response: partial response (code 2) vs stable or progressive disease (codes 3 and 4).
  No patient in the input files has a complete response (code 1).
- **Panel.** 39 flow-cytometry markers at T0. The effector-memory column is named `EM` on the T0 sheet and `EMRA` on
  the T1 sheet, so the longitudinal panel has 38 markers.
- **Step 01 QC and transformation.** Samples with more than 40% missing markers and markers with more than 20%
  missing values are removed; within each response group, samples whose Mahalanobis distance on the first principal
  components (at most 5) exceeds the 99% chi-square quantile are removed as outliers. Percentages are bounded away from 0
  and 100 (epsilon in the config), logit-transformed, imputed by Bayesian PCA and z-scored. Soluble markers would be
  log2-transformed (none in the final panel).
- **Frames.** T0 cross-sectional frame: every T0 sample that passes QC. Added-value frame: patients with T0 and T1
  samples after QC and a pre-therapy NLR (`clinical_model.complete_case_vars`). LMM frame: every T0 sample plus the
  T1 samples of the added-value frame (`longitudinal.t1_restrict_to_delta_frame`).
- **Clinical model.** PD-L1 TPS + NLR; NLR is not imputed (the frame requires it); missing PD-L1 values are
  median-imputed inside each cross-validation fold.
- **Survival.** Overall survival and the survival status are read from the T0 sheet. Progression-free survival is a
  proxy derived in step 06 (no progression date is recorded).
- **Step 07 duplicated-row sensitivity.** The patients in `data/curated/subsets/duplicate_t1_rows.txt` have identical
  T1 rows; step 07 re-fits the LMM dropping each in turn and then both.

## 5. Known data anomalies (kept as recorded)

These are reported, not corrected:

- Two patients with identical T1 rows (step 07 sensitivity above).
- One T0 sample whose Ki67 and CD137 memory-subset values follow the same composition pattern as the re-gated T1
  samples.
- One T1 sample with CD28⁺Ki67⁺ larger than total Ki67⁺.
- Fourteen samples in which the Ki67 subsets do not sum to the Ki67 total.
- CD4 total computed as 100 − CD3 total by a formula in part of the T0 sheet.
- One T0 record that duplicates another patient's record (removed by step 00, rule 10).
- Overall survival recorded differently on the T0 and T1 sheets for four patients (the analysis reads T0).
- The sampling window of the pre-therapy NLR is not recorded; radiotherapy between T0 and T1 is known only where
  annotated.
