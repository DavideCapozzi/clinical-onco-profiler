# tests/golden — regression guardrail

These tests freeze the **publishable headline numbers** of the pipeline and fail
(exit ≠ 0) if a code change moves them. They exist so that the refactoring work
(helper extraction, god-file splits) can proceed safely: *no refactor is allowed
to change a number that could end up in a manuscript.*

The pipeline is seed-fixed (`seed: 2026`), so every metric — including the
permutation p-values — is deterministic and exactly assertable.

## What is protected

| Test | Cohort / experiment | Asserted (frozen in `expected_metrics.yml`) |
|---|---|---|
| `replicate_original.R` | NSCLC original 77, **no split** (oldresults conditions) | SVM AUC **0.7158**, perm p **0.001**, nested-LOO **0.7542**, gate, n_sig_fdr — the *sacred* numbers |
| `check_snapshot.R` | `BestResponse_2v3_4` (v2, discovery-73) + `HNSCC_Response` | perm p, nested-LOO, FDR count, gate, AUC of the current honest results |

## Running

```bash
RSCRIPT=/home/laboratorio/micromamba/envs/clinical-onco-profiler/bin/Rscript
$RSCRIPT tests/golden/replicate_original.R   # ~8-9 min  — the sacred guardrail
$RSCRIPT tests/golden/check_snapshot.R       # ~12-15 min — anti-drift snapshot
```

`check_snapshot.R` runs `main.R` scoped to the two experiments (via the
`EXPERIMENTS` env filter) and therefore **updates `results/latest`** (run id
suffixed `_golden_snapshot`).

## Design (single source of truth)

- JSON paths + field locations live **only** in `R/utils_metrics.R`
  (`extract_run_metrics()`), shared by these tests and — from Fase 3 — by the
  `/post-change` skill. No JSON path is hardcoded twice.
- Comparison logic lives **only** in `assert_metrics.R`.
- Expected values live **only** in `expected_metrics.yml`.

## Re-freezing the snapshot (intentional changes only)

If a change is *meant* to move a number (e.g. a new cohort, a deliberate method
change), update `expected_metrics.yml` in the **same commit**, with the run id /
git SHA that produced the new values noted in the commit message. Never edit it
to make a red test green without understanding why it moved.
