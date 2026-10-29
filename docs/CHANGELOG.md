# Changelog

Analysis history of the NSCLC response study and the record of the 2026-10-28 refactor. The final run is
`results/0final/final_gate84_nlr50` (2026-10-29); `results/0final/primary_PIregated_noPS` and
`results/0final/nlr_v2_cc50` are its sensitivity analyses. Every number from an earlier run is superseded and is not
quoted anywhere.

## Analysis history

- **2026-08-22**: first full run on the parsed database (pre-label-fix).
- **2026-09-10**: audit found response labels misassigned for a group of patients by the parser.
- **2026-09-14**: label-fixed data files. Analysis cohort fixed: patients with an undefined T1 response target removed
  at both timepoints (`data/0final/_upstream_excl6/make_excl6_inputs.R`). This with-PS run is archived:
  `results/_archive/20260914_135629_EXPLORATORY_excl6_no_nlr/`.
- **2026-10-19/20**: ECOG PS removed from every analysis (PI decision); clinical model = PD-L1 TPS alone. This
  no-PS run is archived: `results/_archive/no_PS_20261020_112735_EXPLORATORY_excl6_no_PS_no_nlr/`.
- **2026-10-22**: selection-aware LRT fixed to the top-K rule (K = observed gate size, the rule the Methods describe).
- **2026-10-23**: T1 Ki67 data error. Three patients had the six Ki67-subset columns on the wrong denominator; the PI
  re-gated those samples and only those 21 cells were replaced (`data/0final/primary_PIregated_noPS/`). This is the
  final dataset.
- **2026-10-26/27**: therapy regimen parsed from the free-text field (parser v3.1, mixed code legend; open with the PI);
  therapy sensitivity (ICI monotherapy vs chemo-immunotherapy).
- **2026-10-27**: Kaplan-Meier median split of the Δ composite (post hoc, descriptive).
- **2026-10-27**: NLR v2. The whole pipeline on the 50 Δ-frame patients with a pre-therapy NLR, clinical model
  PD-L1 + NLR (`data/0final/nlr_v2_cc50/`). It replaces the NLR v1 test (gate selected on the larger LMM frame).
- **2026-10-28**: refactor (below). Both canonical runs re-run on the refactored code and verified equal to the
  2026-10-23 / 2026-10-27 runs (`results/0final/VERIFICATION.md`: 136/136 and 129/129 checks pass; every stochastic
  quantity bit-identical).
- **2026-10-28**: NLR run re-run on `6fab3a4` (each clinical variable alone in CV, with CI, and in Figure 2b): verified
  129/129 against the reference; versus the previous NLR run only the new keys and rows were added (that run was
  then deleted). NLR alone
  0.624 [0.426, 0.784], PD-L1 alone 0.339 [0.274, 0.619] (n = 50). The primary run was not re-run: with one clinical
  variable the change adds nothing.
- **2026-10-29**: final design. Gate selected on all longitudinal data (LMM 84 / 144, as in the primary); added value
  tested on the 50 paired patients with a pre-therapy NLR, clinical model PD-L1 + NLR, NLR not imputed. Rationale:
  selection uses every sample; the clinical-vs-immune comparison is made only where the clinical model is complete.
  This is the NLR v1 design (first run 2026-10-27, selection-aware p 0.022), now in the pipeline as the config option
  `clinical_model.complete_case_vars` (default off; off = byte-identical to before) and run end-to-end as
  `results/0final/final_gate84_nlr50` (`config/final_gate84_nlr50.yml`, inputs `data/0final/final_gate84_nlr50/`).
  Verified A1–A11 (`results/0final/VERIFICATION.md`): steps 01–04 identical to the primary; CV 0.570 / 0.751 / 0.764;
  LR 12.98; selection-aware p 0.022 (223/10000, same count as v1); LOO 0/50 below the null q95. The primary (Δ 60,
  PD-L1 alone, p 0.044) and NLR v2 (whole pipeline on the 50, p 0.037) are reported as sensitivity analyses; the
  three are nested analyses of the same patients, never compared as improvements.

## Refactor 2026-10-28 (branch `refactor/final-results`, tag `pre-refactor-20261028` = state before)

- Inputs frozen in `data/0final/` (md5 in each `MD5SUMS`); `config/global_params.yml` = primary config,
  `config/nlr_v2_cc50.yml` = NLR config.
- `RUN_NAME` names the run dir (refuses an existing dir unless `FORCE=1`); the `latest` pointer is gone; tools take
  `--run`; every run writes `run_manifest.yml` (git sha, config and input md5, seeds, N per frame).
- Removed: step 05 (networks), the HNSCC / Comorbidity / Toxicity experiments, the univariate path, 8 legacy step-06
  figures, the dead selection-aware builder, `R/modules_ml_plots.R`, `R/modules_network.R`, `diagnostics/`.
- Figures exported by one helper, `save_figure()` (formats and dpi from config `figures:`).
- Step 06 now computes the selection-aware LRT, the patient-bootstrap CIs and the LRT leave-one-out.
- New step 07 (`src/07_post_hoc.R`): grouping, duplicated T1 rows, therapy, KM.
- `publication/` collector with `INDEX.md`, driven by config `publication:`.
- Identifiable exports go to `private/<run name>/` (gitignored).

### Old name → new name

| old | new |
|---|---|
| `diagnostics/diag_44_selection_aware_lrt.R` | `run_selection_aware_lrt()`, `R/modules_ml_selection_aware.R` (step 06) |
| `diagnostics/diag_45_bootstrap_delta_auc.R` | `run_patient_bootstrap_auc()`, `R/modules_ml_bootstrap_ci.R` (step 06) |
| fixed-gate LRT leave-one-out (computed inline) | `run_lrt_leave_one_out()`, `R/modules_ml_utility.R` (step 06) |
| `diagnostics/diag_50_response_grouping.R` | `run_response_grouping_sensitivity()`, `R/modules_posthoc_grouping.R` (step 07 `grouping/`) |
| `diagnostics/diag_55_duplicate_t1_sensitivity.R` | `run_duplicate_row_sensitivity()`, `R/modules_posthoc_duplicates.R` (step 07 `duplicate_rows/`) |
| `diagnostics/diag_56_therapy_stratification.R`, `tools/therapy_sensitivity.R` | `R/modules_therapy_sensitivity.R` via step 07 `therapy/` |
| `tools/km_median_split.R` | `R/modules_survival_km.R` via step 07 `km/` |
| `diagnostics/diag_53_data_integrity.R` | `tools/data_integrity_audit.R --run --raw` |
| `diagnostics/diag_54_compare_runs.R` | `tools/compare_runs.R` |
| `diagnostics/_av_common.R` (harness) | `av_*` helpers in `R/modules_ml_utility.R` |
| `diag_00`, `diag_46`–`49`, `diag_51`, `diag_51b`, `diag_57`, `run_d44_d45.sh`, `run_nlr*.sh` | removed (superseded or set aside) |
| `results/final_PIregated/<run>` | `results/0final/primary_PIregated_noPS` |
| `results/nlr_v2/<run>` | `results/0final/nlr_v2_cc50` |
| `<run>/BestResponse_2v3_4/post_hoc/{km_median,therapy_v3.1}` | `<run>/BestResponse_2v3_4/07_post_hoc/{km,therapy}` |
| `<run>/BestResponse_2v3_4/06_machine_learning/publication/` (all figures) | also collected into `<run>/publication/{main,supplementary,additional}` |

Earlier one-off diagnostics (diag_09–14, the HNSCC prototypes, the multilevel T0 utility prototype) were already
retired before the refactor; their manifest was `diagnostics/archive/README.md`
(`git show pre-refactor-20261028:diagnostics/archive/README.md`).

### Archive and deletion record

- Backup of `results/`, `diagnostics/`, `data/`, `config/`, `manuscript/` and `CLAUDE.md` before any change:
  `~/backups/clinical-onco-profiler_pre-refactor_20261028.tar.gz`, 1316 entries, sha256
  `42f8cf4766b537a805dda555fbd8b4b8b210d56dcd30a855b59b94fc048c0edd`.
- Archived (`results/_archive/`, with a README): the 2026-10-20 no-PS run (with its `_inputs/` and README) and the
  2026-09-14 with-PS run.
- Run READMEs, comparison tables and the original configs and console logs of the 2026-10-23 and 2026-10-27 runs:
  `docs/archive/run_history/` (local only).
- Deleted from `results/`: `20260822_173722_no_nlr`, `20260914_134721_labelfix_no_nlr`,
  `20260914_163809_EXPLORATORY_76only`, `BestResponse_2v3_4`, `Comorbidity_1v0`, `Toxicity_1v0`, the rest of
  `no_PS/` (two 2026-10-19 runs, `_inputs/`, README, pointers), `with_PS/`, `ki67_sensitivity/`,
  `ki67_reconstructed/`, `ki67_PI_regated/`, `grouping_sensitivity/`, `final_PIregated/`, `nlr_v2/`,
  `data_integrity_audit/`, `latest`, `latest_run.txt`. Deleted: `diagnostics/` (including `diagnostics/out/`).
  Everything deleted is in the backup above.

## After the canonical runs (2026-10-28, no effect on any number)

- `main.R` exits with status 1 when any pass or the publication collector logs a `[FATAL ERROR]` (the other passes
  still run first) and `run_manifest.yml` lists them under `failed_passes` (empty list = clean run). Before this, a
  failed run exited 0 and only the log grep caught it.
- Header and inline comments of `config/global_params.yml` and `config/nlr_v2_cc50.yml` rewritten: superseded numbers
  (with-PS era AUCs and naive LRT p), ECOG PS wording and pointers to deleted diagnostics removed. Apart from the
  publication `name:` key below, the parsed YAML is identical (every analysis setting unchanged), but the file md5
  changed: `global_params.yml` b5293bcd… → 11e967e4…, `nlr_v2_cc50.yml` 4f09a9d9… → 56365076…. The `config_md5` in
  the two canonical `run_manifest.yml` files is the pre-edit md5; the old files are `git show 03fcb81:config/<file>`.
- `publication/` copies are named by paper label (`Fig_1_Biology.pdf`, `Fig_S1_CONSORT.pdf`, `Table_S1_….xlsx`;
  unlabelled files: descriptive part only, without the legacy `FigureN_` prefix and experiment suffix). The step
  outputs keep their names; `INDEX.md` gained a source column. Both canonical `publication/` folders were re-collected
  (same 25 files, byte-identical contents). Config entries accept an optional `name:` (used for
  `AddedValue_LOOcurves`).
- Each clinical variable alone (only when the clinical model has more than one, i.e. the NLR run): 1-variable
  logistic scored on the same repeated 10-fold CV folds as the three arms (no extra random draws), with its
  patient-bootstrap CI (gate set, same resamples). New: JSON `cv_kfold.auc_single`, rows `AUC <var> alone` in
  `PatientBootstrap_AUC_*.csv`, `single_<var>` curves in Figure 2b (ROC and net benefit). These CV curves replace the
  LOO "PD-L1 alone" curve that Figure 2b of the NLR run used to draw next to CV curves. With a one-variable clinical
  model (primary run) nothing is added.
