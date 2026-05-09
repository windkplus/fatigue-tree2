# Paper #2 Outline — Probabilistic Eccentric Decay and Wind Fatigue of Urban Street Trees

> Working title: *Probabilistic Eccentric Decay Geometry: Monte Carlo Confidence Intervals on the Wind-Fatigue Life of Urban Street Trees*
>
> Target journal (primary): **Trees - Structure and Function** (Springer)
> Target journal (alternate): **Engineering Failure Analysis** (Elsevier)
> Companion to paper #1 (Trees / IJF dual submission, 2026).
> Author: Junsuk Kang (Seoul National University)
> Status: Outline draft — ready for full-text writing once Phase C2 (sensitivity) closes.

---

## Headline Claims (one-sentence summaries)

1. The concentric hollow-cylinder decay assumption used in Mattheck (1994) and in our paper #1 *systematically over-estimates* fatigue life relative to the geometrically realistic eccentric/asymmetric decay pockets observed in field tomography.
2. Under a Beta(2, 2) prior on cavity eccentricity and a uniform prior on cavity orientation, the *median* fatigue life at the Mattheck threshold $r_d/R = 0.7$ is **25 % of the concentric estimate** (lower-bound S-N), and the *5 th-percentile* tail is shorter by a factor of 47.
3. Under a uniform cavity-orientation prior, the wind-direction distribution (isotropic, Seoul rose, Jeju rose) has < 0.3 % effect on the life ratio — i.e. *cavity geometry, not wind climate, is the dominant uncertainty source*.
4. A risk-based inspection schedule keyed to the 5 th-percentile fatigue-life lower bound (rather than the deterministic concentric estimate) should be adopted; for ginkgo at $r_d/R = 0.5$ in Seoul this lowers the threshold from "25 yr" (concentric) to "$\leq$ 5 yr for 21 % of trees" (eccentric MC).

---

## 1. Introduction (~ 1,000 words)

### 1.1 Motivation
- Urban tree windthrow risk; brief recap of paper #1's framework and findings.
- The decay-accelerated low-cycle fatigue paradigm: stress amplification from $r_d/R^4$ singularity is the dominant lever.

### 1.2 Limits of the concentric assumption
- The classical Mattheck (1994) hollow-cylinder formula $\sigma/\sigma_0 = 1/(1-(r_d/R)^4)$ is the de facto standard.
- Sonic tomography surveys (Brazee 2011; Wang 2007; Goh 2018) show real cavities are *eccentric*, *non-circular*, and frequently *multi-pocketed*.
- Mattheck and Breloer (1994) themselves note that asymmetric defects can amplify peak stress by 20-40 % at the same $r_d/R$, but their treatment is illustrative rather than probabilistic.
- Paper #1 (Section 4.5, fifth limitation) explicitly identifies this as a non-conservatism that warrants probabilistic treatment.

### 1.3 Research questions (carry over from `docs/research_plan.md`)
1. How does eccentric decay geometry change fatigue life relative to the concentric assumption?
2. Can deterministic point-estimates of fatigue life be replaced by 95 % confidence intervals derived from observable cavity-geometry priors?
3. Does this enable a risk-based inspection schedule?

### 1.4 Contributions
- A vectorisable, analytically grounded framework for arbitrary eccentric / elliptical cavities (extends Mattheck);
- A 30-scenario Monte Carlo sweep (5 decay levels × 3 wind distributions × 2 S-N curves, $n = 10\,000$ each) anchored to the calibrated absolute lives of paper #1;
- Demonstration that wind-direction distributions decouple from the life ratio under uniform cavity priors;
- A policy-ready decay-ratio × percentile risk table for Korean street-tree management.

---

## 2. Methods (~ 1,500 words)

### 2.1 Eccentric cross-section model
- Outer disk radius $R$ (intact), single decay cavity of arbitrary geometry.
- Closed-form moments-of-inertia for circular and elliptical cavities (Section 2.1.1, Section 2.1.2).
- Composite-area / parallel-axis derivation; net centroid shift; verification against Cartesian grid integration.
- Reference: Boresi & Schmidt (2003) Advanced Mechanics of Materials.

### 2.2 Stress Amplification Factor (SAF) under wind bending
- Define $\text{SAF}(\theta_w) = S_0 / S(\theta_w)$ where $S_0 = \pi R^3 / 4$ is the sound section modulus and $S(\theta_w)$ is the eccentric section modulus for a wind direction $\theta_w$.
- General formula derived in Eq. (4-6) of `stress_amplification.py` docstring.
- *Continuity check*: $\text{SAF}(\theta_w) \to 1/(1-(r_d/R)^4)$ as eccentricity $\to 0$ for any $\theta_w$ (Theorem A.1, Appendix).

### 2.3 Monte Carlo sampling
- **Cavity-geometry priors (Phase B confirmed)**:
  - Normalised eccentricity $e/(R - r_d) \sim \text{Beta}(2, 2)$ (symmetric, weakly informative).
  - Cavity offset angle $\sim U[0, 2\pi)$.
- **Wind-direction priors**: isotropic, Seoul 8-bin empirical rose, Jeju 8-bin empirical rose (Section 2.5).
- **Joint**: cavity orientation independent of wind direction (defended in Section 4.x).
- **Sample size**: $n = 10\,000$ per scenario, seed = 42 (reproducible).
- **Sweep**: $r_d/R \in \{0.3, 0.5, 0.6, 0.7, 0.8\}$.

### 2.4 Damage-multiplier formulation
- Scaling argument: under multiplicative stress amplification and a power-law S-N curve $\sigma \propto N^{-b}$, Miner damage scales as $\text{SAF}^{1/b}$.
- Hence the per-sample fatigue life is

  $$T_{\text{ecc}}(\text{sample}) = T_{\text{baseline}} \,\Big/\, \mathbb{E}_{\theta_w}\!\bigl[\text{SAF}(\theta_w; \text{sample})^{1/b}\bigr].$$

- Reported as the *life ratio* $\rho = T_{\text{ecc}} / T_{\text{conc}}$ which is invariant to $T_{\text{baseline}}$ and lets us anchor to paper #1.

### 2.5 Anchoring to paper #1 absolute lives
- Direct re-use of paper #1 Table 2 values (Seoul, $H = 8\;\text{m}$, $\text{DBH} = 20\;\text{cm}$, ginkgo).
- $T_{\text{ecc, sample}} = T_{\text{conc, paper #1}} \times \rho_{\text{sample}}$ — exact under the multiplicative-stress assumption.
- Treats paper #1's full pipeline (DAF, rainflow, Weibull) as a calibrated "black-box" baseline.

### 2.6 Software & reproducibility
- All code: `https://github.com/.../fatigue-tree2` (release at submission).
- 4-stage regression test (`scripts/00_regression_test_concentric.py`) verifies the framework reproduces Mattheck (1994) to machine precision in the concentric limit.
- Seed = 42, NumPy 2.x, Python 3.11+.

---

## 3. Results (~ 1,500 words)

### 3.1 Concentric-limit verification
- Figure 1 (proposed): regression test results — 7 decay levels × 4 tests, all passing.
- Headline: $|\text{SAF}_{\text{model}} - \text{SAF}_{\text{Mattheck}}| < 10^{-15}$ (machine precision).

### 3.2 Eccentricity-driven SAF amplification (single-sample analysis)
- **Figure 2 (proposed)**: 6-panel SAF demo (`output/figures/eccentric_demo.png`).
  - (a) Cross-section schematics for concentric, eccentric, elliptical cases.
  - (b) Polar SAF($\theta_w$) for the same cases.
  - (c) Worst-case SAF vs. eccentricity at four decay levels.
  - (d) SAF vs. decay ratio: concentric / eccentric ($e_{\text{norm}} = 0.5$) / elliptical 2:1.
  - (e) Net-centroid shift vs. eccentricity.
  - (f) Failure-direction selection (worst $\theta_w$ vs. eccentricity).
- Headline: at $r_d/R = 0.7$ with the cavity tangent to the outer wall ($e_{\text{norm}} = 0.95$), SAF is **2.16 ×** the Mattheck value.

### 3.3 Monte Carlo life-ratio distributions
- **Figure 3 (proposed)**: 5-panel CDF (`fig_mc_cdf.png`).
- **Figure 4 (proposed)**: percentile bands vs. $r_d/R$ (`fig_mc_percentile.png`).
- **Figure 5 (proposed)**: life-shortening heatmap (`fig_mc_heatmap.png`).
- **Headline figure (likely Fig. 6)**: histogram at $r_d/R = 0.7$ (`fig_mc_concentric_vs_ecc.png`).

#### 3.3.1 Lower-bound S-N (b = 0.125)
| $r_d/R$ | $\rho_{p5}$ | $\rho_{p50}$ | $\rho_{p95}$ | $1/\rho_{p5}$ |
|---------|------------|-------------|-------------|---------------|
| 0.3     | 0.341      | 0.672       | 0.942       | 2.9           |
| 0.5     | 0.083      | 0.413       | 0.865       | 12            |
| 0.6     | 0.042      | 0.322       | 0.826       | 24            |
| 0.7     | **0.021**  | **0.250**   | 0.778       | **47**        |
| 0.8     | 0.013      | 0.205       | 0.753       | 78            |

#### 3.3.2 Mean S-N (b = 0.10)
| $r_d/R$ | $\rho_{p5}$ | $\rho_{p50}$ | $\rho_{p95}$ | $1/\rho_{p5}$ |
|---------|------------|-------------|-------------|---------------|
| 0.7     | **0.007**  | 0.167       | 0.729       | **143**       |
| 0.8     | 0.004      | 0.129       | 0.699       | 250           |

### 3.4 Wind-direction distribution decouples
- The three wind distributions (isotropic, Seoul, Jeju) produce indistinguishable life-ratio CDFs (max $|\Delta \rho_{p50}| < 0.003$, see Table B.1).
- **Theorem (Section 4.x)**: under independent Uniform priors on cavity orientation, the convolution with the wind-direction distribution is itself uniform.
- Practical implication: site-specific wind-rose data is *not* required for the eccentric-decay correction; cavity-geometry priors are.

### 3.5 Absolute-life anchoring (paper #1 Table 2)
- **Figure 7 (proposed)**: absolute-year distribution at $r_d/R = 0.7$ (`fig_absolute_distribution_rd07.png`).
- **Figure 8 (proposed)**: percentile bands of absolute years vs. $r_d/R$ (`fig_absolute_percentile.png`).
- **Figure 9 (proposed)**: fraction of trees with life $< T_{\text{threshold}}$ vs. $r_d/R$ (`fig_inspection_threshold.png`).

#### Key numbers (Seoul, $H = 8\;\text{m}$, DBH = 20 cm, ginkgo, lower-bound S-N)

| $r_d/R$ | concentric (paper #1) | eccentric $T_{p5}$ | $T_{p50}$ | $T_{p95}$ | % below 5 yr |
|---------|----------------------|--------------------|-----------|-----------|--------------|
| 0.50    | 25.5 yr              | 2.1 yr             | 10.5 yr   | 22.1 yr   | **21 %**     |
| 0.60    | 9.5 yr               | 0.4 yr             | 3.1 yr    | 7.8 yr    | **72 %**     |
| 0.70    | 1.7 yr               | 0.04 yr            | **0.42 yr** | 1.3 yr  | **100 %**    |

---

## 4. Discussion (~ 1,500 words)

### 4.1 Why eccentricity matters more than wind direction
- Independence assumption + uniform cavity orientation $\Rightarrow$ wind-direction integration self-averages.
- Robustness across regions; only the absolute baseline $T_{\text{conc}}$ varies (Seoul vs. Jeju).
- This is a falsifiable prediction: if cavity orientation correlates with wind (e.g. lee-side decay), it will break.

### 4.2 Comparison with Mattheck and Breloer (1994)
- M & B's "20-40 % stress increase" is consistent with our $\rho_{p50}$ ($\approx 0.25 - 0.41$) at moderate decay.
- Our result extends to the *tail*: $\rho_{p5} \approx 0.02$ at the Mattheck threshold has not been quantified before.
- Their $t/R = 0.3$ threshold remains an excellent indicator for the *median* sample, but is non-conservative for the bottom tail.

### 4.3 Inspection-cycle implications
- Under the conservative ("plan for the bottom 5 %") philosophy, a 5-year inspection cycle is non-conservative for $r_d/R \geq 0.5$ in Seoul (Lower-bound S-N).
- Recommendation: inspection interval $\Delta t \leq T_{p5}$, not $\Delta t \leq T_{\text{conc}}$.
- Full risk matrix for adoption by Korea Forest Service / municipal forestry departments (Table 4).

### 4.4 Compatibility with paper #1
- Concentric values reproduced exactly in the limit;
- Absolute-life anchoring uses paper #1 verbatim;
- Eccentric correction acts purely as a multiplicative stress factor — physically the same operator paper #1 used for sound stress.

### 4.5 Limitations of the present work
- (i) Decay is treated as a single circular or elliptical cavity. Multiple pockets, irregular boundaries, and channelled decay (Schwarze 2007) require polygon-based extension.
- (ii) Beta(2, 2) prior is weakly informative; Phase C2 will demonstrate sensitivity to prior choice.
- (iii) Cavity-geometry data from sonic tomography is sparse — proposed campaign with NIFoS (cf. paper #1 future-work item 4).
- (iv) Independence between cavity orientation and wind has not been validated empirically; it is the single largest unverified assumption.
- (v) Time evolution of eccentricity not modelled (paper #3 candidate).

### 4.6 Future work bridges to paper #1's research roadmap
- This paper closes paper #1 future-work item 4 (`docs/future_work.md`).
- Outputs become learning targets for the GNN surrogate (item 6) — paper #3 candidate.
- Coupling with reversed-loading S-N (item 1) would remove the residual upper-bound bias.

---

## 5. Conclusions (≈ 200 words, 5 numbered)

1. The concentric-decay assumption is *systematically non-conservative* in fatigue-life prediction; the median eccentric tree at the Mattheck threshold has only 25 % of the concentric life.
2. The bottom-5 % tail is shorter by a factor of 47 (lower-bound S-N) or 143 (mean S-N), invalidating deterministic single-point estimates for risk-based management.
3. Wind-direction distributions are immaterial to the life-ratio correction under uniform cavity-orientation priors.
4. Risk-based inspection cycles should be keyed to the 5 th-percentile lower bound rather than the concentric point estimate.
5. Cavity-geometry tomography is the single most valuable empirical input to make these CIs site-specific.

---

## Figures (provisional list)

| # | File | Caption (placeholder) |
|---|------|-----------------------|
| 1 | `regression_test_summary.png` (to generate) | Concentric-limit verification |
| 2 | `eccentric_demo.png` ✓ | Single-sample SAF behaviour |
| 3 | `fig_mc_cdf.png` ✓ | MC life-ratio CDF |
| 4 | `fig_mc_percentile.png` ✓ | Percentile bands vs decay ratio |
| 5 | `fig_mc_heatmap.png` ✓ | Life-shortening heatmap |
| 6 | `fig_mc_concentric_vs_ecc.png` ✓ | Distribution at Mattheck threshold |
| 7 | `fig_absolute_distribution_rd07.png` ✓ | Absolute-year distribution at $r_d/R=0.7$ |
| 8 | `fig_absolute_percentile.png` ✓ | Absolute-year percentile bands |
| 9 | `fig_inspection_threshold.png` ✓ | Fraction below inspection thresholds |

## Tables

| # | Source | Content |
|---|--------|---------|
| 1 | `mc_baseline_summary.csv` ✓ | 30-scenario MC summary |
| 2 | `mc_policy_table.csv` ✓ | Life-shortening factors |
| 3 | `absolute_life_summary.csv` ✓ | Anchored absolute years |
| 4 | (to draft) | Inspection-cycle policy matrix |

---

## Outstanding to-do before submission

| # | Item | Priority | Status |
|---|------|----------|--------|
| 1 | Phase C2 sensitivity (Beta prior, aspect ratio) | High | Pending |
| 2 | Generate Fig. 1 from regression test output | Med | Pending |
| 3 | Polish Tables 1 and 4 for paper format | Med | Pending |
| 4 | Full text drafting | High | Pending |
| 5 | English language polish / external review | Med | Pending |
| 6 | Cover letter (Trees journal) | Low | Pending |

---

## References (anticipated, ~ 30)

To be sourced on full-text drafting; minimum core set:

1. Mattheck, C. & Breloer, H. (1994). The Body Language of Trees. HMSO.
2. Boresi, A.P. & Schmidt, R.J. (2003). Advanced Mechanics of Materials, 6th ed.
3. **Kang, J.** (2026). Fatigue Life Assessment of Urban Street Trees under Wind Loading. *Trees / Int. J. Fatigue* (paper #1).
4. Wang, X. et al. (2007). Acoustic tomography for decay detection.
5. Brazee, N.J. et al. (2011). Field assessment of internal decay using sonic tomography.
6. Schwarze, F.W.M.R. (2007). Wood decay under the microscope. *Fungal Biology Reviews* 21:133-170.
7. James, K.R., Haritos, N. & Ades, P.K. (2006). Mechanical stability of trees under dynamic loads. *American Journal of Botany* 93:1522-1530.
8. Yang, C. et al. (2025). Reversed-loading fatigue of wood. *Int. J. Fatigue* 194:108807.
9. KMA. *Climatological Normals of Korea* (2020).
10. Ministry of Land, Infrastructure and Transport (2019). Korean street tree selection guidelines.
11. - 30: TBD on full draft.

---

## Estimated paper length

- Body: ~ 5,000 words
- Figures: 9
- Tables: 4
- References: ~ 30
- Total page count (Trees format): ~ 18 pages

Within Trees-Springer normal-article limit. Submission ready when C2 sensitivity is closed.
