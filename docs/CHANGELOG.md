# Changelog

## Unreleased

- Whole-panel LMM forest (`FigureS13_LMMForest_AllMarkers`, `pub_fig_lmm_forest_all()`): every marker's Time×Group
  β with a model-based 95% CI (Satterthwaite t, recovered exactly from the stored t and P, so each CI excludes 0
  exactly when P < 0.05) and its FDR; the gate markers highlighted. Collected as `additional` in both configs.
  Drawn from the Step-04 report: `tools/render_pub_figures.R` adds it to an existing run without re-running.
  `pub_read_lmm_report()` now also returns the panel table and reads a report without the bootstrap sheet.
  `PUB_MARKER_LABELS` covers the whole FACS panel.
- Each clinical variable + immune score: the repeated stratified CV of Step 06 also scores each clinical
  variable of the formal model together with the immune composite (TPS + immune score, NLR + immune score) on the
  same folds, with no extra random draw, so every existing AUC, curve and CI is unchanged (`av_cv_kfold(pairs,
  curves)`; `av_attach_pair_curves()` adds the arms only after reproducing every published arm exactly, FATAL
  otherwise). JSON `cv_kfold$auc_pair`; patient-bootstrap CIs (`AUC <var> + immune`, `dAUC <var> + immune - <var>
  alone`) in `PatientBootstrap_AUC.csv`. `tools/render_pub_figures.R` rebuilds the arms for a run that predates them.
  Test: `tests/test_added_value_pairs.R`.
- Paper set (`config/0final.yml`): Figure 2 = CV ROC (A) + decision curve (B) + nomogram of the combined model (C)
  in one figure (`Figure2c_AddedValue_Nomogram`, `pub_fig_added_value_nomogram()`; the former Figure 3 is panel C);
  Figure S3 = ROC of TPS and NLR, each alone vs + immune score (`FigurePairROC_ClinicalPlusImmune`,
  `pub_fig_pair_roc()`); the survival figure becomes S4; Table S1 = every panel marker's LMM Time×Group β, P and FDR
  (`Table_LMM_AllMarkers`, `pub_table_lmm_all()`), the two reports become Tables S2–S3. The separate CV-curve and
  nomogram figures are kept as `additional`.
- The whole-panel cell-frequency figure (S14) is removed with its helpers (too large for the supplement; the table
  above replaces it). Test: `tests/test_panel_figures.R` (S13 and the table).
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
