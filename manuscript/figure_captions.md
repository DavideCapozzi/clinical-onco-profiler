# Figure captions — NSCLC clinical+immune added-value manuscript

**Target:** Frontiers in Immunology (Cancer Immunity & Immunotherapy). Exploratory /
hypothesis-generating. All numbers from canonical run `20260622_121316_canonical`
(BestResponse_2v3_4; n=79, 34 RP / 45 SD-PD; n_paired=50; n complete-case for the
clinical model = 55). Figures rendered by `manuscript/figures/make_manuscript_figures.R`
(cairo_pdf, fonts embedded; Okabe–Ito colour-blind-safe palette; panel labels A/B).

> Captions carry every statistic; panels carry none beyond minimal in-plot anchors.
> Cite values to the precision below (canonical = pipeline JSON).

---

## Main figures

**Figure 1. Pre-specified immune gate from longitudinal T-cell dynamics.**
Bootstrap-resampled Time×Group interaction coefficients (β) from the Step-04 linear
mixed model for the three Ki67-expressing T-cell subsets that define the immune
composite (CD28KI67, KI67CD4, KI67NAIVE). Filled blue points = observed β; grey ×
= bootstrap median; horizontal bars = 95% bootstrap CI (B bootstrap resamples);
right-hand labels = fraction of LOO refits retaining FDR<0.05. Negative β = greater
post-immunotherapy contraction of the subset in responders. The composite was
pre-specified from these longitudinal dynamics (Timepoint×Group), **independently of
the baseline T0→response association used downstream**, and is anchored by published
biology (Mazzaschi et al. 2024). Single column.

**Figure 2. Incremental value of the immune composite over a pre-specified clinical
reference model.** Leave-one-out (LOO), leakage-free. **(A)** ROC for the clinical
reference model (PD-L1 + NLR + PS; orange dashed, LOO AUC 0.605) vs clinical+immune
(blue, LOO AUC 0.675). The primary test of the increment is the likelihood-ratio test
(**LRT p = 0.0098**; permutation p = 0.0125; Firth-penalized p = 0.0127) — the correct
nested-model test; the DeLong test on ΔAUC is non-significant (p = 0.143), **as expected
for nested models** where AUC is insensitive to a hierarchically-ordered predictor.
ΔAUC(LOO) = +0.070. **(B)** Decision-curve analysis (net benefit vs threshold
probability; prevalence 0.43): the clinical+immune model yields the highest net benefit
across the clinically relevant threshold range, supporting the LRT as a co-primary axis
of evidence. Supportive: IDI (leakage-free LOO) = **+0.066 [95% CI 0.005, 0.132]**
(see Fig. 3 for its fragility). Double column.

**Figure 3. Calibration and the fragility of the discrimination-improvement metric.**
**(A)** LOO calibration of the clinical+immune model. Points = observed event fraction
within predicted-probability quintiles (bars = SE); dashed orange = fitted calibration
line of the unpenalized logistic model (slope 0.66, under-calibrated; Brier 0.224);
solid blue = ridge-penalized (L2) model (slope 0.87, near-diagonal; Brier 0.221).
Discrimination is essentially unchanged by ridge (combined LOO AUC 0.675 → 0.674):
the mis-calibration reflects mild overfitting of the unpenalized fit, not absence of
signal. **(B)** Leakage-free LOO IDI: unpenalized +0.066 [0.005, 0.132] (excludes 0,
blue) vs ridge +0.034 [−0.024, 0.091] (**includes 0**, grey). The IDI is therefore
fragile and conditional on the fitted model's calibration; the LRT and net-benefit are
the robust axes. Double column.
*Note: ridge per-patient probabilities are not stored in the frozen run; the ridge
calibration curve is drawn from its frozen fitted slope/intercept (standard
calibration-line representation).*

**Figure 4. Robustness and evidence-stability of the incremental value.**
Defensibility summary. **(A)** Likelihood-ratio test p for the immune increment across
every specification tested (log scale; dashed line = α=0.05): asymptotic χ² (p=0.0098),
exact permutation (0.0125), Firth-penalized (0.0127), enriched baselines +smoking
(0.0102) and +tumor burden (0.0085), and the specificity permutation vs 1,000 random
3-marker composites (spec-p=0.048). The LRT anchor remains significant under all — the
increment is not a single fragile p. **(B)** The IDI is honestly fragile: apparent
in-sample +0.075 [0.004, 0.214] and leakage-free LOO +0.066 [0.005, 0.132] exclude 0
(blue), but under ridge penalization +0.034 [−0.024, 0.091] includes 0 (grey). Hence
the LRT (and net benefit, Fig 2B) — not the IDI — is the primary axis of evidence.
Specificity rests on provenance (independent LMM gate + published Ki67 biology), not
numerical rarity (18.5% of random triplets reach LRT<0.05). Double column.

*(Main set: Fig 1 forest, Fig 2 added-value, Fig 3 calibration, Fig 4 robustness.
The PD-L1/benchmark figure is Supplementary — see Figure S2. Fig 4 consolidates the
material otherwise split across the baseline-invariance / specificity / IDI panels.)*

---

## Supplementary figures

**Figure S1. CONSORT / patient-flow.** Sapienza–Policlinico Umberto I NSCLC cohort
(anti-PD1, non-oncogene-addicted): n=89 baseline (T0; 37 RP / 52 SD-PD) → 10 exclusions
(4 high-missingness >40%; 6 PCA outliers) → analytic cohort n=79 (34 RP / 45 SD-PD),
from which n=50 have paired T0+T1 (LMM gate), n=55 are complete-case for the clinical
reference model (PD-L1+NLR+PS), and n=75 have survival follow-up (54 OS / 64 PFS events).
Double column.

**Figure S2. PD-L1 context.** **(A)** Best-response rate by PD-L1 (TPS%) stratum
(neg <1%, low 1–49%, high ≥50%); labels show N and responders per stratum. Categorical
association is non-significant (3-bin Fisher p = 0.348; Cochran–Armitage trend
p = 0.146), i.e. PD-L1 alone is weakly discriminative in this cohort. **(B)** Gate-model
AUC within PD-L1 subgroups (dashed line = AUC 0.5): PD-L1<50% (n=47) AUC 0.626
[0.464, 0.788]; PD-L1≥50% (n=16) AUC 0.254 [0.000, 0.532]. The high-expressor estimate
is **underpowered (n=16)** and not interpretable. Double column.

**Figure S3. Gate-signal decomposition (T0 / T1 / Δ).** Univariate LOO AUC (95% DeLong
CI) of each gate marker at baseline (T0), on-treatment (T1), and within-patient change
(Δ = T1−T0) in the n=50 paired subset (dotted line = AUC 0.5). Establishes that baseline
(T0) signal drives the classifier while Δ adds little.

**Figure S4. Invariance of the increment to clinical-baseline specification.** LRT p of
the immune increment (log scale; dashed line = 0.05) against the pre-specified reference
model and two enriched baselines: reference PD-L1+NLR+PS (LRT 0.0098, IDI +0.066);
+smoking (LRT 0.0102, IDI +0.074); +tumor burden (LRT 0.0085, IDI +0.080). The increment
is robust to baseline specification. (Apparent clinical-alone AUC decreases under
enrichment — baseline overfitting, not evidence of immune strength; see Limitations.)

**Figure S5. Specificity vs random marker composites.** Null ΔAUC percentiles from 1000
random 3-marker composites drawn from the 39-marker panel; blue line = observed Ki67-gate
ΔAUC (0.052, ≈93rd percentile). Specificity-permutation p(LRT) = 0.048, but 18.5% of
random triplets reach LRT<0.05 — therefore specificity rests on **provenance** (the
independent Step-04 LMM + published Ki67 biology), **not numerical rarity**.

**Figure S6. Standalone immune classifier (secondary; not the headline).** Nested-LOOCV
ROC for the standalone classifiers on the gate features: SVM-RBF (AUC 0.611, permutation
p ≈ 0.05) and Elastic Net (AUC 0.415). The standalone classifier is weak on the full v2
cohort and motivated the pivot to the incremental-value framework; it is retained as the
honest backstory only and is **cut from the main text**.

---

## Cross-cutting statements (Methods / Limitations)

- **Tier:** exploratory / hypothesis-generating; single cohort with internal bootstrap
  optimism correction only. **External validation in an independent cohort
  (≥40–50 events, complete PD-L1+NLR+FACS) is required** before any clinical claim.
- **Survival is null:** OS Cox c-index 0.675→0.672 (LRT p=0.245); PFS 0.638→0.649
  (LRT p=0.702). The immune signature does not improve prognostication beyond response.
- **Prior art / overlap:** this is the predictive-modeling arm of an already-published
  cohort (Gelibter 2024; Mazzaschi 2024; Cancers 2025; Tuosto/Gelibter/Napoletano 2026);
  the only unpublished contribution is the formal incremental-value ML framework. Disclose
  all four prior slices.
- **Disclosure:** n_paired=50; singular random-effect fits flagged for Ki67-subset LMMs;
  the gate is pre-specified but not re-selected per fold (mild residual selection optimism,
  conditional on the gate).
