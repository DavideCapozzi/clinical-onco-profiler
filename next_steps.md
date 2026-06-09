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

### Fase 3 — Robustezza skill: logica in codice ✅
- `R/utils_metrics.R` ampliato: `publishability_verdict(m)` (perm p<0.05 sul metodo primario,
  nested-LOO>0.6, n_sig_fdr>0 solo path lmm), `compare_runs(base, new)` (tabella delta via
  `extract_run_metrics`), `flatten_run_metrics()` (SSOT del flatten), `list_run_experiments()`.
- `tools/metrics_report.R` (tracked, nuovo): CLI `--new/--base/--exp` → verdetto + delta table;
  solo I/O+printing, tutta la logica/soglie negli helper.
- `metrics_json_paths()` unificato su `.STEP_DIRS`/`step_output_path()` (nomi cartelle step in 1 posto).
- `tests/golden/assert_metrics.R`: `flatten_metrics()` delega a `flatten_run_metrics()` (dedup).
- `/post-change` (SKILL.md): fasi 1b+1c ora un solo comando `tools/metrics_report.R` (niente jq/sprintf);
  tabella runner aggiornata ai runner `EXPERIMENT`-parametrici.
- **Validato**: golden assertions re-check PASS (no rerun pipeline); CLI testato (base+new, verdict-only,
  default-latest) — delta 0 metric changed.

### Fase 4 — Spezzare i god file (protetta dal golden test) ✅
- `R/modules_ml.R` (2623 righe) → split per responsabilità: `modules_ml_gate.R` (4 fn),
  `modules_ml_cv.R` (10), `modules_ml_utility.R` (6), `modules_ml_plots.R` (5).
  `modules_ml.R` resta come **aggregatore** (sorgia i 4 via `here()`) così
  `src/06: source("R/modules_ml.R")` funziona invariato; sotto il glob `R/*.R` i 4 file
  sono caricati direttamente (re-source idempotente).
- **Validato**: deparse di tutte e 25 le funzioni **identico** pre/post split; `replicate_original`
  PASS (lmm) + `check_snapshot` PASS (univariate HNSCC + v2). Nessun cambio comportamento.
- `src/06_machine_learning.R` (1415 righe): **non splittato** — è uno script procedurale sequenziale
  (non una libreria di funzioni); spezzarlo in sub-script che condividono variabili sarebbe fragile e
  di basso valore. Deferred.

## Fuori scope ora (deprioritizzato)
`renv.lock`, riconciliazione `.gitignore`↔CLAUDE.md, unificazione README/CLAUDE.md, CI.
