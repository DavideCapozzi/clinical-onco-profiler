# clinical-onco-profiler

An R-based, config-driven pipeline for immunological profiling of oncology cohorts. Designed for hybrid immunophenotyping data (FACS flow cytometry + soluble biomarkers), it performs:

- **Cross-sectional analysis** (Steps 01–03): QC → visualization → sPLS-DA multivariate modeling
- **Longitudinal analysis** (Steps 01 + 04): joint T0/T1 processing → Linear Mixed Models (LMM) with Leave-One-Out sensitivity and covariate adjustment
- **Differential network analysis** (Step 05, optional): bootstrap partial correlation baselines per clinical group → permutation-tested differential overlay → Cytoscape-ready CSV exports
- **Machine learning classification** (Step 06): nested-LOOCV Elastic Net + SVM-RBF classifier gated on LMM LOO-robust features

Applied to NSCLC patients treated with immune checkpoint inhibitors (ICI) to identify response biomarkers. Converging evidence across three analytical layers (sPLS-DA, LMM, ML) points to a **KI67 proliferation axis that moves with early radiological response**. The reported result is the *incremental* value of that axis over a PD-L1 + ECOG-PS baseline — it is a single pre-specified marker with family-level corroboration, **not a validated multi-marker signature and not a standalone classifier**.

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
# Full pipeline (reads config/global_params.yml)
Rscript main.R

# With a custom config (e.g. a new cohort) — no repo edits
CONFIG_PATH=/path/to/alt_config.yml Rscript main.R

# Restrict a run to named experiments
EXPERIMENTS=BestResponse_2v3_4 Rscript main.R
```

The pipeline auto-detects whether the config is **nested multi-experiment** (has an `experiments:` key) or **flat single-cohort**, then runs each experiment through two independent passes (standard cross-sectional + longitudinal), followed by an optional ML pass. At startup `validate_config()` checks that each enabled experiment's outcome and clinical-model columns exist in its input Excel, so a misconfigured cohort fails fast. Set `enabled: false` on an experiment block to skip it without deleting it.

### Tests

```bash
Rscript tests/integration_smoke.R           # end-to-end smoke test on a synthetic cohort (01→04→06)
Rscript tools/metrics_report.R --new <run_id>  # headline metrics + publishability verdict for a run
```

---

## Project structure

```
config/                    Pipeline configuration (YAML)
  global_params.yml        Multi-experiment config (single source of truth)

data/                      Input Excel files (not tracked in git — add your own)

diagnostics/               Supplementary stand-alone analyses (gitignored — local only)
  _av_common.R             Rebuilds the added-value objects outside the pipeline
  diag_00_verify_harness.R Asserts the harness reproduces the published run — run this first

env/
  environment.yml          Conda environment specification (human-editable)
  environment.lock.yml     Fully-resolved lock (exact versions + builds)

R/                         Reusable module library (sourced by main.R at startup)
  modules_qc.R             PCA outlier detection, QC pipeline
  modules_coda.R           Logit/log2 transformation, BPCA imputation
  modules_multivariate.R   sPLS-DA fitting and extraction
  modules_longitudinal.R   LMM fitting, LOO sensitivity, R² (Nakagawa 2013)
  modules_network.R        Bootstrap partial correlations, differential overlay
  modules_split.R          Discovery/validation cohort split
  modules_ml.R             Aggregator — sources the four modules_ml_* files
  modules_ml_gate.R        LMM feature gate, per-fold gate selection
  modules_ml_cv.R          Nested-LOOCV classifiers, collinearity filter, metrics
  modules_ml_utility.R     Clinical benchmark, added-value layer, utility, dissociation
  modules_ml_plots.R       ML ggplot/PDF rendering
  modules_pub_figures.R    Manuscript figure set (pub_fig_*, pub_render_all)
  modules_pub_style.R      Publication theme, palette, figure saving
  modules_viz.R            ggplot rendering functions
  utils_io.R               Config loading, run dirs, cohort ledger, manifests
  utils_metrics.R          Headline metric extraction, publishability verdict, run compare

tools/
  metrics_report.R         Publishability verdict / run-to-run comparison CLI

src/                       Pipeline step scripts (sourced in order by main.R)
  01_data_processing.R
  02_visualization.R
  03_statistical_analysis.R
  04_longitudinal_lmm.R
  05_network_analysis.R
  06_machine_learning.R

tests/                     Test scripts
  integration_smoke.R      End-to-end smoke test on a synthetic cohort

results/                   All pipeline outputs (not tracked in git)
  latest_run.txt           Source of truth for "the last run"
  <run_id>/                One timestamped dir per invocation, e.g. 20260822_173722_no_nlr
    run_manifest.yml       run_id, git_sha, git_dirty, r_version, seed, passes
    <experiment_name>/
      01_data_processing/  Processed .rds objects
      02_visualization/    Distribution and PCA PDFs
      03_statistical_analysis/ sPLS-DA results, JSON drivers
      04_longitudinal_analysis/ LMM tables (Excel + JSON), volcano/trajectory PDFs
      05_network_analysis/ Cytoscape CSVs, network PDFs
      06_machine_learning/ Classifier results (Excel, PDFs, JSON)
        publication/       Manuscript-ready PDFs + publication_data.rds
```

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
| `network` | Bootstrap/permutation counts, edge thresholds, Cytoscape export filters |
| `machine_learning` | FDR/LOO thresholds, collinearity filter, CV grids |

See `config/global_params.yml` for inline documentation of each parameter.

---

## Input data format

Input files are Excel workbooks (`.xlsx`) with one row per patient, one column per marker. Required columns: `Patient_ID` plus all marker columns listed under `features.facs` and `features.soluble` in the config. Longitudinal analyses expect separate T0 and T1 files with matching `Patient_ID` values.

Input files are **not tracked in git** (the `data/` directory is gitignored) to protect patient data.

---

## Supplementary diagnostics

Scripts in `diagnostics/` are stand-alone analyses that extend or validate the main pipeline results. They read processed `.rds` outputs from `results/` and can be run independently. The directory is **gitignored**, so its contents vary by working copy.

`diagnostics/_av_common.R` rebuilds the added-value objects (composite, clinical design matrix, outcome) outside the pipeline, so a diagnostic can vary one thing at a time without re-implementing the composite or the cross-validation.

```bash
Rscript diagnostics/diag_00_verify_harness.R   # ALWAYS run first
```

`diag_00_verify_harness.R` asserts that the harness reproduces the published cross-validated AUCs and likelihood-ratio test **exactly**. If it fails, every downstream diagnostic is invalid — fix the assembly, never adjust the target.

These scripts require the main pipeline to have been run first (so that `results/` outputs exist).

---

## Citation

If you use this pipeline, please cite

---

## License

This project is released for academic use. Contact the authors for other uses.
