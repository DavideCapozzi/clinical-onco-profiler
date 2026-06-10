# Lavorare con Claude Code su questo progetto — note pratiche

Appunti ricavati dalla sessione di refactoring (Fasi 2–4). Obiettivo: ridurre i
tempi morti e gli errori ricorrenti nelle prossime conversazioni.

> Questo file è meta-documentazione personale. Valuta se aggiungerlo a `.gitignore`
> (come `next_steps.md` / `CLAUDE.md`) se non vuoi versionarlo.

---

## 1. Il collo di bottiglia #1: i golden test sono lenti — usa una validazione "a scaglioni"

In questa sessione gran parte del tempo è stato attesa: `replicate_original.R`
(~8–9 min) e `check_snapshot.R` (~12–15 min), rieseguiti per ogni fase. Spesso
**non servivano entrambi**.

**Validazione dal più economico al più costoso — fermati al primo che fallisce:**

| # | Check | Tempo | Quando basta |
|---|---|---|---|
| 1 | `parse()` dei file modificati | secondi | sempre, primo gate |
| 2 | **deparse-equivalence** vs HEAD | secondi | refactor "puro" (spostamento/dedup di codice) |
| 3 | scenari di caricamento (`source` aggregatore + glob `R/*.R`) | secondi | se cambia *come* i file vengono sorgiati |
| 4 | `tools/metrics_report.R --base <run>` su un run esistente | secondi | confronto metriche senza rieseguire la pipeline |
| 5 | `replicate_original.R` (lmm path) | ~9 min | logica ML core / prima del commit |
| 6 | `check_snapshot.R` (univariate HNSCC + v2) | ~15 min | tocchi il path univariate o `main.R` |

**Regola pratica:** per un refactor che NON cambia la logica (es. split di
`modules_ml.R`), il check #2 (deparse identico su tutte le funzioni) è una **prova
più forte** del golden test per i corpi delle funzioni. Lancia il golden solo per
ciò che il deparse non copre (il *caricamento*) e come conferma finale.

**Cosa puoi dirmi per risparmiare ~25 min:** all'inizio di un refactor di' una
frase tipo *"è un refactor puro: usa deparse-equivalence, golden solo come
conferma finale"*. Eviti che io rilanci entrambi i golden ad ogni fase.

### Da aggiungere al repo (una volta)
- `tests/check_refactor_equivalence.R`: confronta il `deparse` di ogni funzione tra
  `git show HEAD:<file>` e la working copy. Riutilizzabile per ogni refactor.
- **Fast smoke test**: una variante dei golden con `n_perm`/`n_boot` ridotti
  (es. 99 / 100) che dà un pass/fail strutturale in ~1–2 min prima del guardrail pieno.

---

## 2. Bug che è costato tempo e come prevenirlo (R-specifici)

Durante la sessione è emerso un bug latente in `tests/golden/assert_metrics.R` che
ha bloccato `check_snapshot` a metà. Due trappole R da mettere nel CLAUDE.md come
**vincoli**, così le evito a priori:

- **`$` fa partial matching sulle liste.** `expected$gate` ha matchato silenziosamente
  `gate_method` (perché la chiave `gate` non esisteva). → **Usare sempre `[["..."]]`**
  per leggere config/liste attese, mai `$`.
- **`x[[i]] <<- v` nel corpo di una funzione (non in una closure) dà
  `object 'x' not found`**: `<<-` salta il frame corrente e cerca nei parent. → usare
  `<-` per gli accumulatori locali; `<<-` solo dentro closure annidate.
- **`vapply(..., integer(1))`** fallisce se la funzione ritorna `double`
  (es. `i + 1`): usare `as.integer()`.

### Prevenzione consigliata
Aggiungere uno **unit test del test-harness** (senza pipeline, fixture statiche):
verifica `assert_block` / `report_assertions` / `publishability_verdict` con casi
pass/fail finti. Cattura questi bug in secondi invece che dentro una run da 15 min.

---

## 3. Igiene git da decidere a monte

- **File di handoff/scratch** (`next_steps.md`): decidere **subito** se tracciato o no.
  In questa sessione è stato committato per errore e poi rimosso con `git rm --cached`.
  → metterli in `.gitignore` al momento della creazione.
- **`CLAUDE.md` è gitignored** in questo repo: comodo (non sporca i commit) ma
  **non è condiviso né versionato** — chi clona non lo riceve, e le sue modifiche non
  sono nella history. Tienilo a mente se mai collabori o cambi macchina.
- **Policy commit già nota** (e rispettata): i commit li fai tu a mano; io do solo il
  messaggio oneline. Eseguo comandi git che toccano indice/history **solo** se lo
  chiedi esplicitamente in quel momento (es. `git rm --cached`).

---

## 4. Skill: cosa usare e cosa aggiungere

**Già utili e usate bene:**
- `/post-change` — ora aggganciata a `tools/metrics_report.R` (verdetto + delta da
  un'unica fonte). Solido. Usalo dopo ogni modifica alla pipeline.

**Da provare / aggiungere:**
- **`/fewer-permission-prompts`** — scansiona i transcript e pre-autorizza i comandi
  read-only ricorrenti (in questa sessione molti `grep`/`wc`/`sed` di sola lettura).
  Riduce i prompt di permesso → meno interruzioni.
- **`/code-review`** (o `/code-review ultra`) sul diff prima del commit, per i refactor
  grossi: caccia bug che i golden (solo numeri headline) non vedono.
- **Skill custom "refactor-verify"**: incapsula la validazione a scaglioni della §1
  (parse → deparse-equivalence → load check). Te la posso creare se vuoi.
- **`update-config`** per hook automatici: es. un hook che lancia il `parse()`-check
  sui file `R/*.R` salvati, o un reminder allo `Stop`.

---

## 5. Prompting: frasi che fanno risparmiare round-trip

- **Dichiara la natura della modifica**: *"refactor puro, numeri invariati"* vs
  *"cambio di metodo, i numeri possono cambiare"* → scelgo la validazione giusta.
- **Fissa il livello di validazione**: *"solo deparse-check, niente golden ora"* oppure
  *"lancia entrambi i golden prima di chiudere"*.
- **Batch delle decisioni a monte**: tracked/untracked dei nuovi file, dove mettere
  i nuovi script (`tools/`, `R/`, `tests/`), se aggiornare CLAUDE.md/memory.
- **Per task lunghi**: ok lasciarmi backgroundare i golden e proseguire; chiedimi un
  *"riassunto del diff finora"* quando vuoi un checkpoint senza aspettare.
- **Quando dici "prosegui"**: se intendi una fase precisa di un piano, nominala — evita
  che io debba inferire l'ambito.

---

## 6. Cose da aggiungere a CLAUDE.md (per evitare errori futuri)

Sintesi degli elementi non ovvi emersi, candidati a entrare in `## Key design constraints`:

1. **Caricamento moduli**: tutti gli `R/*.R` sono sorgiati nel *global env* via
   `list.files(here("R")) |> walk(source)`; `src/06` fa anche `source("R/modules_ml.R")`
   (aggregatore). Nuovi file in `R/` si auto-caricano; il doppio source è idempotente
   (solo def di funzioni).
2. **Trappole R**: `$` partial-match → usare `[[ ]]`; `<<-` solo in closure (vedi §2).
3. **SSOT dei path/metriche**: nomi cartelle step in `.STEP_DIRS` (`utils_io.R`); path
   JSON + soglie + flatten in `utils_metrics.R`; non duplicarli altrove.
4. **Validazione refactor**: deparse-equivalence prima dei golden (vedi §1).

---

## TL;DR (le 3 cose che pesano di più)
1. **Validazione a scaglioni**: deparse-equivalence (secondi) prima dei golden (minuti).
   Dimmi se è un "refactor puro" e salti ~25 min di rerun.
2. **Metti in CLAUDE.md le due trappole R** (`$` partial-match, `<<-`) + uno unit test
   del test-harness: evitano bug latenti scoperti solo dentro run lente.
3. **`/fewer-permission-prompts`** + decidere tracked/untracked dei file scratch a monte.
