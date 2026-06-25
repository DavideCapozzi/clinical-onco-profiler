# Figure captions — NSCLC clinical+immune added-value manuscript

**Target:** Frontiers in Immunology (Cancer Immunity & Immunotherapy). Exploratory /
hypothesis-generating. All numbers from canonical run `20260624_173705_newdata_pdl1fix`
(BestResponse_2v3_4, **value-filled cohort**; n=82, 36 PR / 46 SD/PD; n_paired=49;
complete-case for the clinical model = 62; PD-L1 valid = 79). Figures rendered by the
pipeline (Step 06) from live objects via `R/modules_pub_style.R` +
`R/modules_pub_figures.R` (cairo_pdf, fonts embedded; Okabe–Ito colour-blind-safe palette;
panel labels A/B; 120/190 mm @ 11 pt). Re-render without re-analysis:
`Rscript manuscript/figures/make_manuscript_figures.R`.

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
(T1) in the n=49 paired patients, split by best response — partial response (PR, blue)
versus stable/progressive disease (SD/PD, vermillion); grey lines connect each patient's
two timepoints. Responders (PR) start higher and contract after anti-PD1. **(B)** Linear
mixed-model Time×Group interaction coefficient (β; point = observed estimate, bar = 95 %
bootstrap CI, B = 500 resamples) for the three subsets: CD28KI67 β = −0.90 [−1.37, −0.45];
KI67NAIVE β = −1.32 [−2.07, −0.57]; KI67CD4 β = −1.04 [−1.66, −0.50] (all FDR = 0.038, LOO
worst-case p ≤ 0.006); right-hand labels = % of bootstrap resamples retaining FDR < 0.05.
Negative β = greater post-immunotherapy contraction in responders. The composite was
pre-specified from these longitudinal dynamics — **independently of the baseline T0→response
association used downstream** — and anchored by published biology (Mazzaschi et al. 2024).
Double column.

**Figure 2. Incremental value of the immune composite over a pre-specified clinical
reference model.** Leave-one-out (LOO), leakage-free. **(A)** ROC for the clinical
reference model (PD-L1 + NLR + PS; orange dashed, LOO AUC 0.560) versus clinical + immune
(blue, LOO AUC 0.639), with **PD-L1 alone (red dot-dash, AUC 0.52, n = 79)** as the
companion-biomarker reference. Notably **PD-L1 alone is non-discriminative (0.52) and the
full clinical model adds little (0.56)** — NLR + PS contribute marginally — whereas the
immune composite lifts discrimination to 0.639. The **primary test of the increment is the
likelihood-ratio test** (**LRT p = 0.0075**; exact permutation p = 0.009; Firth-penalized
p = 0.0098) — the correct nested-model test; the DeLong test on ΔAUC is non-significant
(p = 0.116), **as expected for nested models** where AUC is insensitive to a
hierarchically-ordered predictor. ΔAUC(LOO) = +0.079. **(B)** Decision-curve analysis (net
benefit versus threshold probability; prevalence 0.44): clinical + immune yields the highest
net benefit across the clinically relevant threshold range, supporting the LRT as a
co-primary axis of evidence. Supportive: leakage-free LOO IDI = **+0.069 [95 % CI 0.006,
0.136]** (apparent in-sample +0.077 [0.008, 0.210]; see Fig. S6 for its fragility under
penalization). Clinical NAs in-fold median-imputed (n = 82; complete-case n = 62 is a
sensitivity analysis). *PD-L1 alone is the **apparent** univariate benchmark on its
**complete cases (n = 79)** — a single externally-fixed biomarker for which LOO is
uninformative (near-null logistic); in (B) its decision curve uses the n = 79
complete-case denominator and is a reference, not strictly comparable to the n = 82 LOO
curves.* Double column.

**Figure 3. PD-L1 context.** Positioning the immune gate against PD-L1, the FDA companion
biomarker, in this cohort. **(A)** Best-response rate by PD-L1 (TPS %) stratum — negative
< 1 % (34 %; N = 32, R = 11), low 1–49 % (50 %; N = 26, R = 13), high ≥ 50 % (48 %; N = 21,
R = 10); dashed line = 50 %. PD-L1 alone is weakly discriminative here (3-bin Fisher p =
0.442; Cochran–Armitage trend p = 0.293; n_valid = 79). **(B)** Immune-gate-model AUC
computed within PD-L1 subgroups (dashed line = 0.5): PD-L1 < 50 % (n = 58) AUC 0.615 [0.466,
0.764]; PD-L1 ≥ 50 % (n = 21) AUC 0.236 [0.022, 0.451]. The high-expressor estimate is
**underpowered (n = 21) and inverts below 0.5 — interpret as no usable signal on that
subset, not negative information**; not interpretable. Double column.

---

## Supplementary figures

**Figure S1. CONSORT / patient flow.** Sapienza–Policlinico Umberto I NSCLC cohort
(anti-PD1, non-oncogene-addicted): n = 91 baseline (T0; 38 PR / 53 SD-PD) → 9 exclusions
(4 high-missingness > 40 %; 5 PCA outliers) → analytic cohort n = 82 (36 PR / 46 SD-PD), of
which n = 49 have paired T0+T1 (LMM gate), n = 62 are complete-case for the clinical
reference model (PD-L1 + NLR + PS), and n = 81 have survival follow-up (59 OS / 71 PFS
events). Double column.

**Figure S2. Invariance of the increment to clinical-baseline specification.** LRT p of the
immune increment (log scale; dashed line = 0.05) against the pre-specified reference model
and two enriched baselines: reference PD-L1+NLR+PS (LRT 0.0075, IDI +0.069); + smoking (LRT
0.0076, IDI +0.077); + tumour burden (LRT 0.0052, IDI +0.087). The increment is robust to
baseline specification. (Apparent clinical-alone LOO AUC does not improve under enrichment
— + smoking 0.521, + burden 0.565 versus reference 0.560 — i.e. no available factor
strengthens the baseline; this is invariance, not "survives a stronger model"; see
Limitations.) Double column.

**Figure S3. Specificity versus random marker composites.** Null ΔAUC percentiles from 1000
random 3-marker composites drawn from the 39-marker panel (50th/90th/95th/99th pct = 0.007 /
0.045 / 0.058 / 0.076); blue line = observed Ki67-gate ΔAUC (0.058, ≈ 94th percentile).
Specificity-permutation p(LRT) = 0.05, but **19.8 % of random triplets reach LRT < 0.05** —
therefore specificity rests on **provenance** (the independent Step-04 LMM + published Ki67
biology), **not numerical rarity**. Single column.

**Figure S4. Standalone immune classifier (secondary; not the headline).** Nested-LOOCV ROC
for standalone classifiers on the gate features: SVM-RBF (AUC 0.525, permutation p = 0.342)
and Elastic Net (AUC 0.393, permutation p = 0.954). The standalone classifier is null on the
value-filled cohort and motivated the pivot to the incremental-value framework; retained as
the honest backstory only and **cut from the main text**. Single column.

**Figure S5. Dynamics↔baseline coupling (gate rationale).** Across the 39-marker panel, each
marker's Step-04 LMM Time×Group interaction strength (|t|) versus its standalone baseline
(T0) AUC. Immune-gate markers (blue) versus other panel markers (grey); dotted line = 0.5;
grey line = linear fit (95 % CI). Interaction strength correlates with baseline discrimination
(Pearson r = +0.59 [95 % CI 0.34, 0.77], p < 0.001; Spearman ρ = +0.39, p = 0.017; holds
excluding the gate, r = +0.38, p = 0.026) — selecting markers on longitudinal dynamics
enriches for baseline-discriminating markers (regression-to-the-mean: responders start higher
and contract more), explaining why the dynamically-selected gate also predicts at baseline.
The gate ranks 1–3 on interaction strength but only 3–5 on baseline AUC, and the two
highest-T0-AUC markers (KI67CM, PD1KI67) are **not** in the gate — so the gate is **not** a
baseline-AUC-maximizing selection. Convergent validity / gate rationale; bounds (does not
eliminate) the residual selection optimism; **not** evidence of optimality. Single column.

**Figure S6. Robustness and evidence-stability of the increment.** **(A)** LRT p for the
immune increment across every specification tested (log scale; dashed line = 0.05):
asymptotic χ² (0.0075), exact permutation (0.009), Firth-penalized (0.0098), enriched
baselines + smoking (0.0076) and + tumour burden (0.0052), and the specificity permutation
versus 1000 random 3-marker composites (0.05) — all significant. **(B)** The IDI is honestly
fragile: apparent in-sample +0.077 [0.008, 0.210] and leakage-free LOO +0.069 [0.006, 0.136]
exclude 0 (blue), but under ridge (L2) penalization +0.030 [−0.026, 0.083] **includes 0**
(grey). Hence the LRT (and net benefit, Fig 2B) — not the IDI — is the primary axis of
evidence. Double column.

**Figure S7. Calibration and the fragility of the discrimination-improvement metric.**
**(A)** LOO calibration of the clinical+immune model. Points = observed event fraction within
predicted-probability quintiles (bars = SE); dashed orange = unpenalized logistic
(slope 0.59, under-calibrated; Brier 0.235); solid blue = ridge-penalized (L2)
(slope 0.67, better-calibrated; Brier 0.235). Discrimination is only modestly reduced by
ridge (combined LOO AUC 0.639 → 0.620): the mis-calibration reflects mild overfitting of the
unpenalized fit, not absence of signal. **(B)** Leakage-free LOO IDI: unpenalized +0.069
[0.006, 0.136] (excludes 0) versus ridge +0.030 [−0.026, 0.083] (**includes 0**); the IDI is
fragile and conditional on the fitted model's calibration. *Note: ridge per-patient
probabilities are not stored in the frozen run; the ridge calibration curve is drawn from its
fitted slope/intercept (standard calibration-line representation).* Double column.

**Figure S8. Gate-signal decomposition (T0 / T1 / Δ).** Univariate LOO AUC (95 % DeLong CI)
of each gate marker at baseline (T0), on-treatment (T1), and within-patient change
(Δ = T1 − T0) in the n = 49 paired subset (dotted line = 0.5):

| Marker | T0 | T1 | Δ (T1−T0) |
|---|---|---|---|
| CD28KI67 | 0.71 [0.56, 0.86] | 0.64 [0.49, 0.80] | 0.76 [0.62, 0.90] |
| KI67NAIVE | 0.73 [0.59, 0.88] | 0.66 [0.50, 0.81] | 0.78 [0.65, 0.91] |
| KI67CD4 | 0.69 [0.53, 0.84] | 0.63 [0.47, 0.78] | 0.77 [0.65, 0.90] |

Each subset discriminates responders at baseline (T0 AUC 0.69–0.73) — signal present
pre-treatment and usable prospectively. The within-patient change Δ has the highest point
estimate (0.76–0.78), **but is not significantly different from T0** (DeLong T0 vs Δ
p = 0.15–0.57); i.e. baseline and dynamics carry **statistically comparable** signal. The
deployed added-value composite uses **T0 only** for baseline prospective applicability
(Δ requires a second on-treatment sample). Double column.

**Figure S9. The increment survives nested gate re-selection (selection-optimism test).**
The pre-specified-gate added-value (Fig 2) is honest *conditional on the fixed gate*; here the
LMM Time×Group gate is **re-selected inside every outer LOO fold** so the gate-choice
variability is priced in (diag_30, in-pipeline `run_nested_gate_validation`). **(A)** The
increment is invariant: clinical→clinical+immune LOO AUC is 0.560→0.639 (ΔAUC +0.079) for the
pre-specified gate and 0.560→0.641 (+0.081) when the gate is re-selected per fold (DeLong
p = 0.109, n.s., as expected for nested models). **(B)** Selection frequency across folds: the
three pre-specified markers (CD28KI67, KI67CD4, KI67NAIVE; blue) are **recovered in 100 % of
folds** (+ KI67CM, a collinear Ki67 marker, also 100 %; grey); the other 35 panel markers are
never selected. The immune increment is therefore **not an artifact of fixing the gate** — it
is a robustly re-selected Ki67-proliferation module. (The nested rule also picks the high-T0
marker KI67CM → consistent with the dynamics↔baseline coupling, Fig S5; this is convergent
validity, *not* independence from baseline.) Double column.

---

## Cross-cutting statements (Methods / Limitations)

- **Tier:** exploratory / hypothesis-generating; single cohort with internal bootstrap
  optimism correction only (optimism-corrected combined AUC 0.688; apparent 0.732,
  optimism 0.044). **External validation in an independent cohort (≥ 40–50 events, complete
  PD-L1 + NLR + FACS) is required** before any clinical claim.
- **Survival is null:** OS Cox c-index 0.671 → 0.663 (LRT p = 0.321); PFS 0.633 → 0.645
  (LRT p = 0.546). The immune signature does not improve prognostication beyond response.
- **Prior art / overlap:** this is the predictive-modeling arm of an already-published cohort
  (Gelibter 2024; Mazzaschi 2024; Cancers 2025; Tuosto/Gelibter/Napoletano 2026); the only
  unpublished contribution is the formal incremental-value ML framework. Disclose all four
  prior slices.
- **Disclosure:** n_paired = 49; singular random-effect fits flagged for Ki67-subset LMMs;
  the gate is pre-specified but not re-selected per fold (mild residual selection optimism,
  conditional on the gate); PD-L1 ≥ 50 % subgroup (Fig 3B) underpowered.
