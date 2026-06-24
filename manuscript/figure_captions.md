# Figure captions — NSCLC clinical+immune added-value manuscript

**Target:** Frontiers in Immunology (Cancer Immunity & Immunotherapy). Exploratory /
hypothesis-generating. All numbers from canonical run `20260623_114344_canonical`
(BestResponse_2v3_4; n=79, 34 PR / 45 SD/PD; n_paired=50; complete-case for the clinical
model = 55; PD-L1 valid = 63). Figures rendered by the pipeline (Step 06) from live
objects via `R/modules_pub_style.R` + `R/modules_pub_figures.R` (cairo_pdf, fonts
embedded; Okabe–Ito colour-blind-safe palette; panel labels A/B; 120/190 mm @ 11 pt).
Re-render without re-analysis: `Rscript manuscript/figures/make_manuscript_figures.R`.

**Abbreviations** (define on first figure use): PR, partial response; SD, stable disease;
PD, progressive disease; LMM, linear mixed model; LOO, leave-one-out; LRT, likelihood-ratio
test; IDI, integrated discrimination improvement; DCA, decision-curve analysis; AUC, area
under the ROC curve; CI, confidence interval; NLR, neutrophil-to-lymphocyte ratio; PS,
ECOG performance status; TPS, tumour proportion score.

> Figure set rebalanced toward immunology (2026-06-24): **3 MAIN figures** = Fig 1 biology,
> Fig 2 added value, Fig 3 PD-L1 context. Defense/robustness material (gate-signal
> decomposition, calibration, robustness, nested-selection validation) is Supplementary
> — invariance/statistical figures have little main-figure payoff and read as over-defense.
> Captions carry every statistic; panels carry only minimal in-plot anchors. Cite values to
> the precision below (canonical = pipeline JSON).

---

## Main figures

**Figure 1. Pre-specified immune gate: longitudinal Ki67⁺ T-cell dynamics.**
The immune composite is the mean z-score of three Ki67-expressing T-cell subsets
(CD28KI67, KI67NAIVE, KI67CD4). **(A)** Per-patient cell frequency (% of the relevant
parent population; **log₁₀ scale**) of each subset at baseline (T0) versus on-treatment
(T1) in the n=50 paired patients, split by best response — partial response (PR, blue)
versus stable/progressive disease (SD/PD, vermillion); grey lines connect each patient's
two timepoints. Responders (PR) start higher and contract after anti-PD1. **(B)** Linear
mixed-model Time×Group interaction coefficient (β; point = observed estimate, bar = 95 %
bootstrap CI, B = 500 resamples) for the three subsets: CD28KI67 β = −0.90 [−1.35, −0.46];
KI67NAIVE β = −1.30 [−2.09, −0.57]; KI67CD4 β = −1.05 [−1.72, −0.47] (all FDR = 0.038, LOO
worst-case p ≤ 0.006); right-hand labels = % of bootstrap resamples retaining FDR < 0.05.
Negative β = greater post-immunotherapy contraction in responders. The composite was
pre-specified from these longitudinal dynamics — **independently of the baseline T0→response
association used downstream** — and anchored by published biology (Mazzaschi et al. 2024).
Double column.

**Figure 2. Incremental value of the immune composite over a pre-specified clinical
reference model.** Leave-one-out (LOO), leakage-free. **(A)** ROC for the clinical
reference model (PD-L1 + NLR + PS; orange dashed, LOO AUC 0.605) versus clinical + immune
(blue, LOO AUC 0.675). The **primary test of the increment is the likelihood-ratio test**
(**LRT p = 0.0098**; exact permutation p = 0.0125; Firth-penalized p = 0.0127) — the correct
nested-model test; the DeLong test on ΔAUC is non-significant (p = 0.143), **as expected for
nested models** where AUC is insensitive to a hierarchically-ordered predictor.
ΔAUC(LOO) = +0.070. **(B)** Decision-curve analysis (net benefit versus threshold
probability; prevalence 0.43): clinical + immune yields the highest net benefit across the
clinically relevant threshold range, supporting the LRT as a co-primary axis of evidence.
Supportive: leakage-free LOO IDI = **+0.066 [95 % CI 0.005, 0.132]** (apparent in-sample
+0.075 [0.004, 0.214]; see Fig. S7 for its fragility under penalization). Clinical NAs
in-fold median-imputed (n = 79; complete-case n = 55 is a sensitivity analysis). Double column.

**Figure 3. PD-L1 context.** Positioning the immune gate against PD-L1, the FDA companion
biomarker, in this cohort. **(A)** Best-response rate by PD-L1 (TPS %) stratum — negative
< 1 % (35 %; N = 26, R = 9), low 1–49 % (52 %; N = 21, R = 11), high ≥ 50 % (56 %; N = 16,
R = 9); dashed line = 50 %. PD-L1 alone is weakly discriminative here (3-bin Fisher p =
0.348; Cochran–Armitage trend p = 0.146; n_valid = 63). **(B)** Immune-gate-model AUC
computed within PD-L1 subgroups (dashed line = 0.5): PD-L1 < 50 % (n = 47) AUC 0.626 [0.464,
0.788]; PD-L1 ≥ 50 % (n = 16) AUC 0.254 [0, 0.532]. The high-expressor estimate is
**underpowered (n = 16) and inverts below 0.5 — interpret as no usable signal on that
subset, not negative information**; not interpretable. Double column.

---

## Supplementary figures

**Figure S1. CONSORT / patient flow.** Sapienza–Policlinico Umberto I NSCLC cohort
(anti-PD1, non-oncogene-addicted): n = 89 baseline (T0; 37 PR / 52 SD-PD) → 10 exclusions
(4 high-missingness > 40 %; 6 PCA outliers) → analytic cohort n = 79 (34 PR / 45 SD-PD), of
which n = 50 have paired T0+T1 (LMM gate), n = 55 are complete-case for the clinical
reference model (PD-L1 + NLR + PS), and n = 75 have survival follow-up (54 OS / 64 PFS
events). Double column.

**Figure S2. Invariance of the increment to clinical-baseline specification.** LRT p of the
immune increment (log scale; dashed line = 0.05) against the pre-specified reference model
and two enriched baselines: reference PD-L1+NLR+PS (LRT 0.0098, IDI +0.066); + smoking (LRT
0.0102, IDI +0.074); + tumour burden (LRT 0.0085, IDI +0.080). The increment is robust to
baseline specification. (Apparent clinical-alone AUC *decreases* under enrichment — baseline
overfitting on n = 79, not evidence of immune strength; see Limitations.) Double column.

**Figure S3. Specificity versus random marker composites.** Null ΔAUC percentiles from 1000
random 3-marker composites drawn from the 39-marker panel (50th/90th/95th/99th pct = 0.008 /
0.044 / 0.054 / 0.073); blue line = observed Ki67-gate ΔAUC (0.052, ≈ 93rd percentile).
Specificity-permutation p(LRT) = 0.048, but **18.5 % of random triplets reach LRT < 0.05** —
therefore specificity rests on **provenance** (the independent Step-04 LMM + published Ki67
biology), **not numerical rarity**. Single column.

**Figure S4. Standalone immune classifier (secondary; not the headline).** Nested-LOOCV ROC
for standalone classifiers on the gate features: SVM-RBF (AUC 0.611, permutation p = 0.049)
and Elastic Net (AUC 0.415, permutation p = 0.893). The standalone classifier is weak on the
full v2 cohort and motivated the pivot to the incremental-value framework; retained as the
honest backstory only and **cut from the main text**. Single column.

**Figure S5. Dynamics↔baseline coupling (gate rationale).** Across the 39-marker panel, each
marker's Step-04 LMM Time×Group interaction strength (|t|) versus its standalone baseline
(T0) AUC. Immune-gate markers (blue) versus other panel markers (grey); dotted line = 0.5;
grey line = linear fit (95 % CI). Interaction strength correlates with baseline discrimination
(Pearson r = +0.65 [95 % CI 0.42, 0.80], p < 0.001; Spearman ρ = +0.53; holds excluding the
gate, r = +0.50, p = 0.002) — selecting markers on longitudinal dynamics enriches for
baseline-discriminating markers (regression-to-the-mean: responders start higher and contract
more), explaining why the dynamically-selected gate also predicts at baseline. The gate ranks
1–3 on interaction strength but only 3–5 on baseline AUC, and the two highest-T0-AUC markers
(PD1KI67, KI67CM) are **not** in the gate — so the gate is **not** a baseline-AUC-maximizing
selection. Convergent validity / gate rationale; bounds (does not eliminate) the residual
selection optimism; **not** evidence of optimality. Single column.

**Figure S6. Robustness and evidence-stability of the increment.** **(A)** LRT p for the
immune increment across every specification tested (log scale; dashed line = 0.05):
asymptotic χ² (0.0098), exact permutation (0.0125), Firth-penalized (0.0127), enriched
baselines + smoking (0.0102) and + tumour burden (0.0085), and the specificity permutation
versus 1000 random 3-marker composites (0.048) — all significant. **(B)** The IDI is honestly
fragile: apparent in-sample +0.075 [0.004, 0.214] and leakage-free LOO +0.066 [0.005, 0.132]
exclude 0 (blue), but under ridge (L2) penalization +0.034 [−0.024, 0.091] **includes 0**
(grey). Hence the LRT (and net benefit, Fig 2B) — not the IDI — is the primary axis of
evidence. Double column.

**Figure S7. Calibration and the fragility of the discrimination-improvement metric.**
**(A)** LOO calibration of the clinical+immune model. Points = observed event fraction within
predicted-probability quintiles (bars = SE); dashed orange = unpenalized logistic
(slope 0.66, under-calibrated; Brier 0.224); solid blue = ridge-penalized (L2)
(slope 0.87, near-diagonal; Brier 0.221). Discrimination is essentially unchanged by ridge
(combined LOO AUC 0.675 → 0.674): the mis-calibration reflects mild overfitting of the
unpenalized fit, not absence of signal. **(B)** Leakage-free LOO IDI: unpenalized +0.066
[0.005, 0.132] (excludes 0) versus ridge +0.034 [−0.024, 0.091] (**includes 0**); the IDI is
fragile and conditional on the fitted model's calibration. *Note: ridge per-patient
probabilities are not stored in the frozen run; the ridge calibration curve is drawn from its
fitted slope/intercept (standard calibration-line representation).* Double column.

**Figure S8. Gate-signal decomposition (T0 / T1 / Δ).** Univariate LOO AUC (95 % DeLong CI)
of each gate marker at baseline (T0), on-treatment (T1), and within-patient change
(Δ = T1 − T0) in the n = 50 paired subset (dotted line = 0.5):

| Marker | T0 | T1 | Δ (T1−T0) |
|---|---|---|---|
| CD28KI67 | 0.70 [0.55, 0.85] | 0.64 [0.49, 0.80] | 0.76 [0.62, 0.90] |
| KI67NAIVE | 0.71 [0.56, 0.86] | 0.66 [0.50, 0.81] | 0.78 [0.65, 0.91] |
| KI67CD4 | 0.68 [0.52, 0.83] | 0.63 [0.47, 0.78] | 0.77 [0.65, 0.90] |

Each subset discriminates responders at baseline (T0 AUC 0.68–0.71) — signal present
pre-treatment and usable prospectively. The within-patient change Δ has the highest point
estimate (0.76–0.78), **but is not significantly different from T0** (DeLong T0 vs Δ
p = 0.15–0.57); i.e. baseline and dynamics carry **statistically comparable** signal. The
deployed added-value composite uses **T0 only** for baseline prospective applicability
(Δ requires a second on-treatment sample). Double column.

**Figure S9. The increment survives nested gate re-selection (selection-optimism test).**
The pre-specified-gate added-value (Fig 2) is honest *conditional on the fixed gate*; here the
LMM Time×Group gate is **re-selected inside every outer LOO fold** so the gate-choice
variability is priced in (diag_30, in-pipeline `run_nested_gate_validation`). **(A)** The
increment is invariant: clinical→clinical+immune LOO AUC is 0.605→0.675 (ΔAUC +0.070) for the
pre-specified gate and 0.605→0.676 (+0.071) when the gate is re-selected per fold (DeLong
p = 0.146, n.s., as expected for nested models). **(B)** Selection frequency across folds: the
three pre-specified markers (CD28KI67, KI67CD4, KI67NAIVE; blue) are **recovered in 100 % of
folds** (+ KI67CM, a collinear Ki67 marker, also 100 %; grey); the other 35 panel markers are
never selected. The immune increment is therefore **not an artifact of fixing the gate** — it
is a robustly re-selected Ki67-proliferation module. (The nested rule also picks the high-T0
marker KI67CM → consistent with the dynamics↔baseline coupling, Fig S5; this is convergent
validity, *not* independence from baseline.) Double column.

---

## Cross-cutting statements (Methods / Limitations)

- **Tier:** exploratory / hypothesis-generating; single cohort with internal bootstrap
  optimism correction only (optimism-corrected combined AUC 0.708; apparent 0.748,
  optimism 0.041). **External validation in an independent cohort (≥ 40–50 events, complete
  PD-L1 + NLR + FACS) is required** before any clinical claim.
- **Survival is null:** OS Cox c-index 0.675 → 0.672 (LRT p = 0.245); PFS 0.638 → 0.649
  (LRT p = 0.702). The immune signature does not improve prognostication beyond response.
- **Prior art / overlap:** this is the predictive-modeling arm of an already-published cohort
  (Gelibter 2024; Mazzaschi 2024; Cancers 2025; Tuosto/Gelibter/Napoletano 2026); the only
  unpublished contribution is the formal incremental-value ML framework. Disclose all four
  prior slices.
- **Disclosure:** n_paired = 50; singular random-effect fits flagged for Ki67-subset LMMs;
  the gate is pre-specified but not re-selected per fold (mild residual selection optimism,
  conditional on the gate); PD-L1 ≥ 50 % subgroup (Fig 3B) underpowered.
