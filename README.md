# clinical-onco-profiler

An R-based, config-driven pipeline for immunological profiling of oncology cohorts. Designed for hybrid immunophenotyping data (FACS flow cytometry + soluble biomarkers), it performs:

- **Cross-sectional analysis** (Steps 01–03): QC → visualization → sPLS-DA multivariate modeling
- **Longitudinal analysis** (Steps 01 + 04): joint T0/T1 processing → Linear Mixed Models (LMM) with Leave-One-Out sensitivity and covariate adjustment
- **Machine learning and added value** (Step 06): nested-LOOCV Elastic Net + SVM-RBF gated on LMM LOO-robust features; clinical + immune added-value layer with the selection-aware LRT (label permutation with full gate re-selection), patient-bootstrap CIs for the cross-validated AUCs (including each clinical variable alone when the clinical model has more than one, drawn in Figure 2b from the same CV folds), and the fixed-gate LRT leave-one-out
- **Post hoc sensitivity** (Step 07): response grouping (PR / SD / PD), duplicated-row sensitivity, therapy regimen, Kaplan-Meier median split; descriptive, gate fixed

Applied to NSCLC patients treated with immune checkpoint inhibitors (ICI) to identify response biomarkers. Converging evidence across three analytical layers (sPLS-DA, LMM, ML) points to a **KI67 proliferation axis that moves with early radiological response**. The reported result is the *incremental* value of that axis over a PD-L1 baseline — it is a single marker axis selected by the LMM, with family-level corroboration, **not a validated multi-marker signature and not a standalone classifier**.

---

## Requirements

- [conda](https://docs.conda.io) or [micromamba](https://mamba.readthedocs.io)
- R 4.3.x (managed via the conda environment)

---

## Installation

```bash
conda env create -f env/environment.yml
conda activate clinical-onco-profiler
```

Or with micromamba:

```bash
micromamba env create -f env/environment.yml
micromamba activate clinical-onco-profiler
```

The environment pins R 4.x with Bioconductor packages (`mixOmics`, `pcaMethods`, `speckle`, `ALDEx2`) and key CRAN packages (`lmerTest`, `lme4`, `glmnet`, `e1071`, `pROC`, `igraph`). `env/environment.yml` is the human-editable spec; for an exactly reproducible build use the fully-resolved lock:

```bash
micromamba create -f env/environment.lock.yml   # exact versions + build strings
```

---

## Usage

```bash
# Canonical runs (~100 and ~90 min; the selection-aware LRT uses 15 cores)
RUN_NAME=primary_PIregated_noPS CONFIG_PATH=config/global_params.yml Rscript main.R
RUN_NAME=nlr_v2_cc50            CONFIG_PATH=config/nlr_v2_cc50.yml    Rscript main.R

# Any other config — no repo edits; restrict to named experiments
RUN_NAME=my_run CONFIG_PATH=/path/to/alt_config.yml EXPERIMENTS=BestResponse_2v3_4 Rscript main.R
```

The run directory is `<output_root>/<RUN_NAME>` (without `RUN_NAME`, a timestamped name). `main.R` refuses to
write into an existing run directory unless `FORCE=1`. Each experiment runs Pass 1 standard (`01→02→03`), Pass 2
longitudinal (`01→04`), Pass 3 ML (`06`) and Pass 4 post hoc (`07`); then `publication/` is collected and
`run_manifest.yml` (git sha, config and input md5, seeds, N per frame, failed passes) is written.

A failed pass logs `[FATAL ERROR]` and the remaining passes still run; at the end `main.R` exits with
status 1 and lists the failed passes. `grep -l "FATAL ERROR" <run>/*/log_*.txt` shows which log to read.

The pipeline auto-detects whether the config is **nested multi-experiment** (has an `experiments:` key) or **flat single-cohort**. At startup `validate_config()` checks that each enabled experiment's outcome and clinical-model columns exist in its input Excel, so a misconfigured cohort fails fast. Set `enabled: false` on an experiment block to skip it without deleting it.

### Tests

```bash
Rscript tests/integration_smoke.R              # end-to-end smoke test on a synthetic cohort (~10 min)
```

### Tools (all read an existing run, named explicitly)

```bash
Rscript tools/metrics_report.R --new <run>        # headline metrics + publishability verdict
Rscript tools/render_pub_figures.R --run <run>    # re-render the figures from publication_data (~10 s), re-collect publication/
Rscript tools/compare_runs.R a=<run> b=<run> [--out f.md]   # side-by-side comparison of runs
Rscript tools/data_integrity_audit.R --run <run> --raw <raw DB>   # raw DB → inputs → step-01 audit (output in private/)
```

`tools/patient_inclusion_map.R` and `tools/patient_delta_export.R` (local only) write identifiable files to
`private/<run name>/`, which is gitignored.

---

## Project structure

```
config/                    Pipeline configuration (YAML)
  global_params.yml        Primary run config
  nlr_v2_cc50.yml          NLR run config (clinical model PD-L1 + NLR)

data/                      Input Excel files (not tracked in git)
  0final/                  Frozen inputs of the canonical runs, with MD5SUMS

docs/
  CHANGELOG.md             Analysis history, refactor record, old → new names

env/
  environment.yml          Conda environment specification (human-editable)
  environment.lock.yml     Fully-resolved lock (exact versions + builds)

R/                         Reusable module library (sourced by main.R at startup)
  modules_qc.R             PCA outlier detection, QC pipeline
  modules_coda.R           Logit/log2 transformation, BPCA imputation
  modules_multivariate.R   sPLS-DA fitting and extraction
  modules_longitudinal.R   LMM fitting, LOO sensitivity, R² (Nakagawa 2013)
  modules_split.R          Discovery/validation cohort split
  modules_ml.R             Aggregator — sources the four modules_ml_* files
  modules_ml_gate.R        LMM feature gate, per-fold gate selection
  modules_ml_cv.R          Nested-LOOCV classifiers, collinearity filter, metrics
  modules_ml_utility.R     Clinical benchmark, added-value layer, utility, dissociation, LRT leave-one-out
  modules_ml_selection_aware.R  Selection-aware LRT (label permutation with gate re-selection)
  modules_ml_bootstrap_ci.R     Patient-bootstrap CIs for the CV AUCs
  modules_posthoc_grouping.R    Step 07: response grouping sensitivity
  modules_posthoc_duplicates.R  Step 07: duplicated-row sensitivity
  modules_therapy_sensitivity.R Step 07: therapy regimen sensitivity
  modules_survival_km.R         Step 07: Kaplan-Meier median split
  modules_pub_figures.R    Manuscript figure set (pub_fig_*, pub_render_all)
  modules_pub_style.R      Publication theme, palette, figure saving
  modules_viz.R            ggplot rendering functions
  utils_io.R               Config loading, run dirs, cohort ledger, manifests
  utils_metrics.R          Headline metric extraction, publishability verdict, run compare
  utils_figures.R          save_figure(): every figure export (formats / dpi from config)
  utils_publication.R      collect_publication(): publication/ folder and INDEX.md

tools/                     Stand-alone CLIs on an existing run (see Tools above)

src/                       Pipeline step scripts (sourced in order by main.R)
  01_data_processing.R
  02_visualization.R
  03_statistical_analysis.R
  04_longitudinal_lmm.R
  06_machine_learning.R
  07_post_hoc.R

tests/                     Test scripts
  integration_smoke.R      End-to-end smoke test on a synthetic cohort

results/                   All pipeline outputs (not tracked in git)
  0final/                  The canonical runs, README.md and VERIFICATION.md
  _archive/                Superseded runs kept for provenance (never quoted)
  <run name>/
    run_manifest.yml       git sha, config + input md5, seeds, N per frame, timestamp
    publication/           Collected figures and tables (see below) + INDEX.md
    <experiment_name>/
      01_data_processing/  Processed .rds objects, cohort ledger
      02_visualization/    Distribution and PCA PDFs
      03_statistical_analysis/ sPLS-DA results, JSON drivers
      04_longitudinal_analysis/ LMM tables (Excel + JSON), volcano/trajectory PDFs
      06_machine_learning/ Classifier and added-value results (Excel, JSON, csv)
        publication/       Manuscript figures + publication_data.rds
      07_post_hoc/{grouping,duplicate_rows,therapy,km}/

private/                   Identifiable exports (names), local only
```

### Publication layout

`collect_publication()` copies files (it never re-renders) into `<run>/publication/`, following the `publication:`
block of the config: `main/`, `supplementary/`, `supplementary/tables/` (LMM and ML xlsx reports) and `additional/`
(other figures and the selection-aware LRT, patient-bootstrap and LRT leave-one-out csv summaries). Copies are named
by their paper label (`main/Fig_2_AddedValue_CVcurves.pdf`, `supplementary/Fig_S1_CONSORT.pdf`,
`supplementary/tables/Table_S1_Longitudinal_LMM_Report.xlsx`; unlabelled files keep only the descriptive part); the
step outputs keep their original names. `INDEX.md` has one row per file: role, paper label (provisional), file,
source step file, producing function, frame and n. A missing main or supplementary file is a `[FATAL ERROR]`.

---

## Configuration

All pipeline parameters are controlled by `config/global_params.yml`. Key sections:

| Section | Description |
|---------|-------------|
| `experiments` | Named experiment blocks (override clinical mapping, input files, flags) |
| `features.facs` / `features.soluble` | Marker panel definitions |
| `clinical.target_column` | Column encoding the response group |
| `clinical.mapping` | Maps raw integer codes to responder/non-responder labels |
| `qc` | Missingness thresholds, outlier detection settings |
| `machine_learning` | FDR/LOO thresholds, collinearity filter, CV grids |
| `machine_learning.clinical_immune_added_value.selection_aware_lrt` | B (10000), top-K rule, cores, seed |
| `machine_learning.clinical_immune_added_value.patient_bootstrap` | B (2000), CV repetitions per resample, cores, seed |
| `post_hoc` | Step 07 blocks: `grouping`, `duplicate_rows`, `therapy` (regimen files, B, seed), `km` (landmark) |
| `figures` | Export formats (`pdf`, `png`, `tiff`, `svg`) and dpi |
| `publication` | Mapping of figures and tables into `publication/` with labels and frames |
| `stats.seed` | Global seed; B values and seeds behind reported numbers must not change |

See `config/global_params.yml` for inline documentation of each parameter.

---

## Input data format

Input files are Excel workbooks (`.xlsx`) with one row per patient, one column per marker. Required columns: `Patient_ID` plus all marker columns listed under `features.facs` and `features.soluble` in the config. Longitudinal analyses expect separate T0 and T1 files with matching `Patient_ID` values.

Input files are **not tracked in git** (the `data/` directory is gitignored) to protect patient data.

---

## Citation

If you use this pipeline, please cite

---

## License

This project is released for academic use. Contact the authors for other uses.
