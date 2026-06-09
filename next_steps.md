# Next steps — refactoring pipeline (handoff per nuova chat)

Roadmap completa: `~/.claude/plans/voglio-che-scrivi-una-breezy-peacock.md`.
Le 3 modifiche concordate in revisione + le decisioni di questa sessione sono qui sotto.

## Stato attuale (fatto)

- **Fase 0 — pulizia git** ✅
  - `diagnostics/` **non più tracciato** (solo `diagnostics/archive/README.md` resta, come manifest).
    Gli script restano su disco; orfani spostati in `diagnostics/archive/`.
  - Untrack artefatti rigenerabili: `Rplots.pdf`, `*.rds`, `*.pdf`, `diag_*_output.txt`.
    `.gitignore` aggiornato (`*.rds`, `*.pdf`, `Rplots.pdf`, `diagnostics/*` con eccezione README).
  - `run_step06_only.R` **tenuto** (è il runner NSCLC Step06-only, NON un duplicato dell'HNSCC).
- **Fase 1 — guardrail golden** ✅
  - `R/utils_metrics.R`: `extract_run_metrics()` — unica fonte dei path JSON + metriche headline (Step 03/04/06).
  - `tests/golden/` (TRACKED): `replicate_original.R` (sacro), `check_snapshot.R` (v2+HNSCC),
    `assert_metrics.R`, `expected_metrics.yml` (valori congelati), `README.md`.
  - `main.R`: filtro `EXPERIMENTS="A,B"` (unset → tutti) per scoping run.
  - Dry-run asserzioni: **PASS** su original_cohort + BestResponse_2v3_4 + HNSCC_Response.
- Commit: dato il messaggio oneline; **committa l'utente manualmente**.

## Vincoli / decisioni (rispettare nelle prossime fasi)

- **Guardrail sacro**: la replica `replicate_original` deve restare `SVM AUC 0.7158 / perm p 0.001 /
  nested-LOO 0.7542 / gate KI67NAIVE+KI67CD4+CD28KI67`. Eseguire `tests/golden/replicate_original.R`
  (≈8-9 min) dopo ogni refactor; per i numeri pubblicabili v2/HNSCC `tests/golden/check_snapshot.R` (≈12-15 min).
- **Commit**: quando l'utente lo chiede, fornire **un solo messaggio oneline** che riassume TUTTE le
  modifiche correnti (no commit per-fase, no etichette "Fase N"). **Mai** eseguire `git add/commit/reset/push`.
  (Vedi skill `/post-change` Fase 2 + memoria `feedback_commits_manual.md`.)
- **CLAUDE.md** è untracked (`.gitignore`) → non entra nei commit; aggiornarlo comunque su disco (compatto).
- **Generalità**: scrivere il codice dataset-agnostico salvo forzature pesanti (richiesta esplicita utente).
- Rscript: `/home/laboratorio/micromamba/envs/clinical-onco-profiler/bin/Rscript`.

## Da fare

### Fase 2 — Compattazione via helper condivisi ✅
- `prepare_pass_config(base_config, exp_cfg, exp_name, pass_mode, run_root)` in `R/utils_io.R`:
  incapsula il blocco null-guarded di `main.R` (Pass 1/2/3) → main.R −55 righe nette.
- `step_dir()` / `step_output_path()` + `.STEP_DIRS` in `R/utils_io.R`: SSOT dei nomi cartelle
  step; sostituiti i ~19 `file.path(config$output_root, "0X_...", ...)` in tutti i `src/01–06`.
- Runner `diagnostics/` generalizzati: `run_step06_only.R` e `run_step04_only.R` ora parametrici
  via env `EXPERIMENT` (default BestResponse_2v3_4); `run_hnscc_step06_only.R` → wrapper di 8 righe;
  `run_hnscc_only.R` su `prepare_pass_config`.
- Wrapper opzionali `save_pdf()`/`save_workbook()`: **saltati** (erano marcati opzionali).
- **Validato**: `replicate_original.R` PASS (0.7158/0.001/0.7542/gate); `check_snapshot.R` PASS
  (BestResponse_2v3_4 + HNSCC_Response, numeri invariati).
- **Bug-fix collaterale** (tracked, scoperto dal golden): `tests/golden/assert_metrics.R` —
  `expected$gate` faceva partial-match su `gate_method` (→ `expected[["gate"]]`), e riga 44 usava
  `<<-` invece di `<-` nel corpo della funzione. Nessun valore atteso cambiato.

### Fase 3 — Robustezza skill: logica in codice
- Ampliare `R/utils_metrics.R` con:
  - `compare_runs(run_a, run_b)` → tabella delta tra due run.
  - `publishability_verdict(metrics)` → verdetto (perm p<0.05, nested-LOO>0.6, n_sig_fdr>0) + flag.
- Aggiornare `/post-change` (SKILL.md) perché chiami questi helper (via piccolo runner R) invece di
  jq/estrazioni ad-hoc → la skill smette di hardcodare path JSON e soglie. Stessa fonte del golden test.
- Unificare `metrics_json_paths()` (in `utils_metrics.R`) su `.STEP_DIRS`/`step_output_path()` —
  rimandato dalla Fase 2 per non sconfinare di fase: ora i nomi cartelle step sono definiti in 2 posti.

### Fase 4 — Spezzare i god file (protetta dal golden test)
- `R/modules_ml.R` (2623 righe) → split per responsabilità: `modules_ml_gate.R`, `modules_ml_cv.R`,
  `modules_ml_utility.R`, `modules_ml_plots.R` + un file che li sorgia. Nessun cambio comportamento.
- Idem in piccolo per `src/06_machine_learning.R` se resta troppo denso.

## Fuori scope ora (deprioritizzato)
`renv.lock`, riconciliazione `.gitignore`↔CLAUDE.md, unificazione README/CLAUDE.md, CI.
