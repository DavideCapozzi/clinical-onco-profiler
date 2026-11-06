# clinical-onco-profiler

Config-driven R pipeline for immune profiling of peripheral blood (PBMC flow cytometry) in a cohort of patients with
non-small cell lung cancer (NSCLC) treated with immune checkpoint inhibitors (ICI). Blood is sampled before therapy
(T0) and after the first cycle (T1). The pipeline asks whether early changes in circulating T-cell subsets track the
best overall response (partial response, PR, vs stable or progressive disease, SD/PD) and whether they add to a
pre-specified clinical model (PD-L1 TPS + pre-therapy neutrophil-to-lymphocyte ratio, NLR).

The analysis is exploratory, single-arm and has no external validation cohort.

## What the pipeline does

| Step | Script | Content |
|---|---|---|
| 00 | `src/00_build_inputs.R` | Curated database → analysis input files (rules in `config/data_build.yml`; see [docs/DATA_PROCESSING.md](docs/DATA_PROCESSING.md)) |
| 01 | `src/01_data_processing.R` | QC (row / column missingness, PCA outliers), detection-limit clamp, logit (flow cytometry %) / log2 (soluble), BPCA imputation, z-scores |
| 02 | `src/02_visualization.R` | Distributions and PCA (QC figures) |
| 03 | `src/03_statistical_analysis.R` | Baseline (T0) sPLS-DA with repeated CV and a label-permutation null of the balanced error rate |
| 04 | `src/04_longitudinal_lmm.R` | Linear mixed models per marker (`marker ~ Time × Group + (1 \| patient)`), FDR, leave-one-patient-out stability, patient-level cluster bootstrap, covariate and response-subgroup sensitivity |
| 06 | `src/06_machine_learning.R` | Marker gate from step 04 (FDR and leave-one-out rule, collinearity filter); immune score = mean z of the gate markers' change; added value over the clinical model: repeated stratified k-fold CV AUC, patient-bootstrap CIs, likelihood-ratio test (LRT) with a **selection-aware permutation null** (the gate is re-selected in every permutation), LRT leave-one-out, calibration and decision curves; exploratory nested-LOOCV classifiers (elastic net, SVM) |
| 07 | `src/07_post_hoc.R` | Post hoc, gate fixed: response grouping (PR / SD / PD), duplicated-row sensitivity, therapy regimen (ICI monotherapy vs chemo-immunotherapy), survival (Kaplan-Meier and Cox by best response and by immune-score change), LMM forest of the sensitivity run |

Step numbering has a gap (05) on purpose: a network step was removed before release.

`main.R` runs, for each experiment of the config: Pass 1 standard (`01 → 02 → 03`), Pass 2 longitudinal
(`01 → 04`), Pass 3 (`06`), Pass 4 post hoc (`07`); then it collects the figures and tables into `publication/` and
writes `run_manifest.yml`.

## Requirements

- [micromamba](https://mamba.readthedocs.io) or conda.
- R 4.3.3 and the packages pinned in `env/`. For an exact rebuild use the lock file (linux-64):

```bash
micromamba create -n clinical-onco-profiler -f env/environment.lock.yml
micromamba activate clinical-onco-profiler
```

`env/environment.yml` lists the direct dependencies with their versions (any platform). `logistf` is optional and
not in the reference environment; without it the Firth LRT is reported as missing.

- 15 CPU cores for the selection-aware LRT and the patient bootstrap (`cores:` in the config; the results do not
  depend on the number of cores). Run one analysis at a time.
- Run every command from the repository root.

## Data availability

The patient-level data are not public: they contain clinical information on identifiable patients. They are
available from the corresponding author on reasonable request, subject to approval by the ethics committee.

The data package has this layout (all under `data/`, which is never committed):

```
data/
  curated/
    DB_NSCLC_curated.xlsx        study database plus the curated supplements (input of step 00)
    subsets/nlr_v2_cc50_ids.txt  patient list of the sensitivity run
    subsets/duplicate_t1_rows.txt  patients with identical T1 rows (step 07 sensitivity)
  snapshots/fig_s1_nlr_v2_cc50/  frozen inputs of Supplementary Figure S1 (optional: can be rebuilt, see below)
  inputs/                        written by step 00
```

What the curated database adds to the study database, and every rule applied to it, is described in
[docs/DATA_PROCESSING.md](docs/DATA_PROCESSING.md).

## Reproducing the results

Two runs are reported: the final analysis (`config/0final.yml`) and its sensitivity analysis on the 50 patients of
the added-value frame (`config/nlr_v2_cc50.yml`). Supplementary Figure S1 of the final run is the LMM forest of the
sensitivity run, so the sensitivity run comes first.

```bash
# 1. Analysis inputs from the curated database
Rscript src/00_build_inputs.R                                   # -> data/inputs/cohort/
Rscript src/00_build_inputs.R --subset data/curated/subsets/nlr_v2_cc50_ids.txt \
                              --out data/inputs/nlr_v2_cc50     # -> data/inputs/nlr_v2_cc50/

# 2. Sensitivity run (~90 min)
RUN_NAME=nlr_v2_cc50 CONFIG_PATH=config/nlr_v2_cc50.yml Rscript main.R

# 3. Freeze the inputs of Figure S1 (skip if data/snapshots/fig_s1_nlr_v2_cc50/ is in the data package;
#    an existing snapshot is never overwritten)
Rscript tools/render_sensitivity_forest.R --from results/0final/nlr_v2_cc50 \
                                          --snapshot data/snapshots/fig_s1_nlr_v2_cc50

# 4. Final run (~110 min)
RUN_NAME=0final CONFIG_PATH=config/0final.yml Rscript main.R
```

Every random draw has a fixed seed in the config (global `stats.seed`, and a seed and number of resamples per
procedure), so the same inputs give the same numbers. A rebuilt snapshot holds the same values as the distributed one,
but its file checksums differ (the copied Excel report carries its own write time).

Expected frames of the final run (`run_manifest.yml`, `experiments: … n_frames`): T0 cross-sectional frame 78 patients; LMM 84
patients / 134 observations; added-value frame 50 paired patients (23 PR / 27 SD/PD). Selection-aware LRT p = 0.024
(236 of 10000 permutations).

### Environment variables of `main.R`

| Variable | Effect |
|---|---|
| `CONFIG_PATH` | Config file (default `config/0final.yml`) |
| `RUN_NAME` | Run directory `<output_root>/<RUN_NAME>`; an existing directory is refused unless `FORCE=1` (then it is replaced) |
| `RUN_LABEL` | Without `RUN_NAME`: suffix of the timestamped run directory |
| `EXPERIMENTS` | Comma-separated experiment names: run only these |
| `INCLUDE_NLR` | `true` / `false`: overrides `include_nlr` (clinical model with or without NLR) |

`src/00_build_inputs.R` also refuses to write into a non-empty output directory unless `FORCE=1`.

A failed pass logs `[FATAL ERROR]`; the other passes still run, then `main.R` exits with status 1 and
`run_manifest.yml` lists the pass under `failed_passes`.

### Tests and tools

```bash
Rscript tests/test_build_inputs.R          # step-00 rules on a synthetic curated workbook (~5 s)
Rscript tests/integration_smoke.R          # end-to-end run on a synthetic cohort (~10 min)
Rscript tools/render_pub_figures.R --run <run dir>         # re-draw the figures of a finished run (~10 s)
Rscript tools/compare_runs.R a=<run dir> b=<run dir> [--out file.md]   # headline numbers side by side
Rscript tools/render_sensitivity_forest.R --run <run dir> --snapshot <dir>   # draw Figure S1 into a finished run
```

Tools only read finished runs; `render_pub_figures.R` re-draws from the persisted figure data and never re-runs an
analysis.

## Outputs

```
results/0final/<RUN_NAME>/
  run_manifest.yml          git commit, config and input checksums, data-file checksums, seeds, package versions,
                            N per frame, failed passes
  publication/              figures and tables of the article, collected by label, with INDEX.md
    main/                   Figures 1–3
    supplementary/          Figures S1–S3; tables/ Tables S1–S2
    additional/             other figures and csv summaries (not in the article)
  BestResponse_2v3_4/       one directory per experiment
    01_data_processing/     processed data (.rds), cohort ledger, QC report
    02_visualization/       QC figures
    03_statistical_analysis/  sPLS-DA
    04_longitudinal_analysis/ LMM report (xlsx) and metrics (json)
    06_machine_learning/    added-value and classifier results (json, csv, xlsx), publication figures
    07_post_hoc/            grouping/ duplicate_rows/ therapy/ km/ bor_survival/ sensitivity_lmm/
    log_*.txt               one log per pass
```

`publication/INDEX.md` gives, for every file, its label, the step output it was copied from, the function that drew
it and the analysis frame with its n. The headline numbers are in
`06_machine_learning/Machine_Metrics_ML_<experiment>.json` (nodes `selection_aware_lrt`, `patient_bootstrap`,
`lrt_leave_one_out`, `clinical_immune_added_value`).

## Repository layout

```
main.R                 pipeline driver
config/                0final.yml (final run), nlr_v2_cc50.yml (sensitivity run), data_build.yml (step 00 rules)
src/                   step scripts 00–07
R/                     function modules (sourced by main.R)
tools/                 command-line tools on finished runs
tests/                 synthetic tests (step 00 rules; end-to-end run)
env/                   conda environment and lock file
docs/                  DATA_PROCESSING.md, CHANGELOG.md
```

The config files document every parameter inline. The values behind the reported numbers (number of resamples and
seeds) must not be changed if the results are to be reproduced.

## License

This project is released for academic use. Contact the authors for other uses.
