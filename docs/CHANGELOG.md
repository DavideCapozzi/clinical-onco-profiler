# Changelog

## 1.0.0 (2026-11-06)

First public release: the code behind the results of the accompanying article (final run `config/0final.yml`,
sensitivity run `config/nlr_v2_cc50.yml`).

- Step 00 (`src/00_build_inputs.R`, `R/modules_data_build.R`, `config/data_build.yml`): builds the analysis inputs
  from the curated database with generic, config-driven rules; no patient identifier or value in the code
  (`docs/DATA_PROCESSING.md`).
- Steps 01–04: QC, logit / log2 transformation, BPCA imputation, sPLS-DA, per-marker linear mixed models with
  leave-one-out, cluster-bootstrap and covariate sensitivity.
- Step 06: LMM marker gate, immune score, added value over the clinical model (CV AUC with patient-bootstrap CIs,
  selection-aware LRT, LRT leave-one-out, calibration, decision curves), exploratory classifiers.
- Step 07: post hoc analyses (response grouping, duplicated rows, therapy regimen, survival, sensitivity-run LMM forest
  from a frozen snapshot).
- `RUN_NAME` / `FORCE` run directories, `run_manifest.yml` (git commit, config, input and data-file checksums, seeds,
  package versions, frames), `publication/` collector with `INDEX.md`, exit status 1 on any failed pass.
- Tools on finished runs: `render_pub_figures.R`, `render_sensitivity_forest.R`, `compare_runs.R`.
- Synthetic tests: `tests/test_build_inputs.R` (step-00 rules), `tests/integration_smoke.R` (end-to-end run).
- Environment: `env/environment.yml` (direct dependencies) and `env/environment.lock.yml` (R 4.3.3, linux-64).
