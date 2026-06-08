# diagnostics/archive

One-off diagnostic scripts that have served their purpose. They are **kept on
disk, not deleted**: their conclusions are already recorded in the agent memory
(`memory/project_diagnostics_09_13.md`, `memory/project_hnscc_extension.md`,
`memory/project_supplementary_analyses.md`) and CLAUDE.md.

> **Note:** the `diagnostics/` tree is intentionally **not tracked in git** for
> now (see `.gitignore`). This README is the only tracked file under
> `diagnostics/` and serves as a manifest of what exists locally. The scripts
> themselves live in each developer's working copy.

| Script | What it explored | Verdict (see memory/) |
|---|---|---|
| `diag_09_gate_composite.R` | Single KI67 composite vs multi-marker gate | composite NO |
| `diag_10_lmm_bootstrap.R` | Bootstrap CIs on LMM betas | bootstrap LMM SÌ (integrated into Step 04) |
| `diag_11_perfold_zscore.R` | Per-fold z-scoring (leakage check for SVM) | per-fold z NO |
| `diag_12_pdl1_stratified.R` | PD-L1-stratified analysis | PD-L1 stratified SÌ |
| `diag_13_combined_pdl1_ki67.R` | IDI/NRI for combined PD-L1+KI67 model | IDI/NRI SÌ |
| `diag_14_t0_delta_signal_decomposition.R` | T0/T1/Delta gate signal decomposition | integrated into Step 06 |
| `diag_hnscc_02_signal_and_models.R` | HNSCC signal + candidate models | superseded by HNSCC_Response experiment |
| `diag_hnscc_03_immune_score_vs_cps.R` | Immune score vs CPS benchmark | superseded by Step 06 CPS benchmark |
| `diag_hnscc_04_stratified_fixed_gate.R` | HNSCC stratified fixed-gate analysis | superseded by HNSCC_Response experiment |
| `multilevel_t0_utility.R` | Multilevel T0 clinical-utility prototype | superseded by `multilevel_t0_calibration.R` |

To resurrect one, `git mv` it back to `diagnostics/`.
