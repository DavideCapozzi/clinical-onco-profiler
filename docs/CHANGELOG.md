# Changelog

## Unreleased

- Publication figures in other formats: `figures.export` (`formats`, `dpi`, `roles`) converts the collected pdfs of
  `publication/main` and `publication/supplementary` with poppler's `pdftocairo` into `<role>/<format>/` (png, tiff,
  jpeg, svg or eps; `config/0final.yml`: png at 600 dpi). The pdfs are unchanged; each copy is listed in `INDEX.md`.
  Any non-pdf format rendered natively (`figures.formats`) is also collected into `<role>/<format>/`.
- `tools/render_pub_figures.R --collect-only` rebuilds `publication/` of a finished run without re-drawing; the
  `--preview` option is removed (superseded by the export).
- `main.R` and the collector stop early when an export is configured and poppler is missing.
- Environment: `poppler=26.02.0` added (installed with the other packages frozen; no other package changed).
- Test: `tests/test_figure_export.R`; `tests/integration_smoke.R` checks the exported copies.

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
