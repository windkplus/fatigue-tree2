# Probabilistic Eccentric Decay Geometry: Monte Carlo Confidence Intervals on the Wind-Fatigue Life of Urban Street Trees

**Junsuk Kang** ᵃᵇᶜᵈ \*

ᵃ Interdisciplinary Program in Landscape Architecture, Seoul National University, Seoul 08826, Republic of Korea
ᵇ Transdisciplinary Program in Smart City Global Convergence, Seoul National University, Seoul 08826, Republic of Korea
ᶜ Department of Landscape Architecture and Rural Systems Engineering, Seoul National University, Seoul 08826, Republic of Korea
ᵈ Research Institute of Agriculture and Life Sciences, Seoul National University, Seoul 08826, Republic of Korea

\* Correspondence: junkang@snu.ac.kr; Tel.: +82-2-880-4872

---

## Abstract

The classical concentric hollow-cylinder model of stem decay (Mattheck and Breloer 1994), reproduced in our earlier work (Kang 2026, hereafter Paper I), assumes that wood-decay cavities are circular and centred on the stem axis. Field tomography surveys show, however, that real cavities are eccentric, irregular, and frequently elongated. This paper relaxes the concentric assumption and quantifies, through a Monte Carlo framework, how cavity-geometry uncertainty propagates into uncertainty on wind-fatigue life. Closed-form section properties are derived for circular and elliptical cavities of arbitrary eccentricity, and a vectorisable stress-amplification function is verified to reproduce the Mattheck formula to machine precision in the concentric limit. Ten thousand cavity geometries are sampled per scenario from a Beta(2, 2) prior on normalised eccentricity and a uniform prior on cavity orientation; the resulting per-sample fatigue lives are anchored to Paper I's calibrated absolute-life pipeline. The median fatigue life at the Mattheck threshold $r_d/R = 0.7$ is found to be 25 % of the concentric estimate (lower-bound S-N), and the 5th-percentile tail is shorter by a factor of 47. Wind-direction distributions (isotropic, Seoul empirical rose, Jeju empirical rose) collapse out under independent uniform cavity-orientation priors, leaving cavity geometry, not wind climate, as the dominant uncertainty source. A first-order Sobol decomposition shows that 99.6 % of the variance in the life ratio is explained by normalised eccentricity alone, providing a clear empirical priority for sonic-tomography surveys. A risk-based inspection policy keyed to the 5th-percentile lower bound rather than the deterministic concentric estimate is proposed and tabulated for three Korean wind regimes.

**Keywords:** wind-fatigue, eccentric decay, Monte Carlo simulation, urban street trees, probabilistic risk assessment, sonic tomography, Mattheck threshold

---

## 1. Introduction

### 1.1 Background and motivation

Urban street trees provide essential ecosystem services — microclimate regulation, stormwater interception, air-quality improvement, and aesthetic value (Nowak and Greenfield 2018) — but during typhoons and severe windstorms they pose significant hazards to pedestrians, vehicles, and infrastructure. In the Republic of Korea, several thousand street-tree failures are reported annually during typhoon events, with damage concentrated in coastal regions and on Jeju Island (Korea Forest Service 2021). Existing tree-risk-assessment practice relies primarily on field-inspection methods such as resistograph drilling, sonic tomography, and visual tree assessment to characterise decay extent (Wang et al. 2007; Brazee et al. 2011), and informs static stability evaluations against critical wind-speed thresholds (James et al. 2006). Until recently, however, the cumulative *fatigue* damage caused by everyday cyclic wind loading was not quantitatively coupled with decay-driven stress amplification.

In our earlier work (Kang 2026, *Trees* / *International Journal of Fatigue*, hereafter Paper I) we developed an integrated framework combining finite-element dynamic analysis, Kaimal-spectrum-based wind time-series, rainflow cycle counting, Miner's linear damage rule, a concentric-hollow decay stress-amplification model, and regional Weibull wind speed distributions. Six tree geometries were analysed under biaxial wind loading, and the resulting absolute fatigue lives were anchored to Korean street-tree species, regional wind regimes, and current management standards. Paper I established that **decay-accelerated low-cycle fatigue** — the combination of decay-induced stress amplification and infrequent high-wind events — is the dominant failure scenario for urban trees, displacing the conventional view that fatigue is irrelevant for sound trees under ordinary winds. Critically, Paper I identified the Mattheck threshold $r_d/R = 0.7$ as the point at which lower-bound fatigue life drops below 1.7 yr in Seoul, necessitating immediate intervention.

### 1.2 The concentric assumption and its limits

The stress-amplification model used throughout Paper I is the classical concentric hollow-cylinder formula (Mattheck and Breloer 1994; Mattheck 2007),

$$
\frac{\sigma_{\text{decayed}}}{\sigma_{\text{sound}}}
\;=\;
\frac{1}{1 - (r_d/R)^4},
$$

which assumes that the decay cavity is circular and centred on the stem axis. The Mattheck $t/R = 0.3$ static-failure criterion derived from this formula is the de-facto standard for arboricultural risk assessment and is incorporated into national tree-management guidelines worldwide (Lonsdale 1999; Smiley et al. 2017).

Field evidence increasingly contradicts the concentric assumption. Sonic tomography surveys consistently reveal cavities that are eccentric, irregular in cross-section, and often elongated in the radial direction following the path of fungal hyphae through annual growth rings (Brazee et al. 2011; Goh et al. 2018; Allison et al. 2020). Mattheck and Breloer (1994) themselves note in passing that asymmetric defects can amplify peak bending stress by 20-40 % at the same nominal $r_d/R$, but their treatment is illustrative rather than probabilistic, and neither they nor subsequent authors have quantified the *distribution* over plausible eccentric geometries. Paper I (Section 4.5, fifth limitation) explicitly identified this as a non-conservatism that warrants probabilistic treatment and recommended the present study as the natural follow-on.

A second, structural difficulty with the concentric assumption is that it produces a single deterministic point estimate of fatigue life for a given $r_d/R$. From the perspective of municipal arboriculture this is a fragile basis for management decisions: the population of street trees with nominal decay ratio 0.7 in fact contains an *ensemble* of geometric realisations, and the variance of this ensemble — including its lower tail — is what determines whether a risk-management threshold is conservative.

### 1.3 Research questions

The present study addresses three questions:

1. **How does cavity eccentricity change the predicted fatigue life relative to the concentric assumption?** We seek a probabilistic answer in the form of a life-ratio distribution $\rho = T_{\text{ecc}} / T_{\text{conc}}$ rather than a single corrected scalar.

2. **Can deterministic point estimates be replaced by 95 % confidence intervals derived from observable cavity-geometry priors?** This requires a sampling framework that propagates eccentricity priors through the same calibrated pipeline used in Paper I.

3. **Does the resulting uncertainty enable a risk-based inspection schedule that is more conservative — and more defensible — than the current Korean standard of fixed-interval visual inspection?**

### 1.4 Contributions

This paper makes four contributions. First, it develops a vectorisable, analytically grounded framework for arbitrary eccentric and elliptical decay cavities, extending the classical Mattheck construction. Second, it presents a 30-scenario Monte Carlo sweep (5 decay levels × 3 wind distributions × 2 S-N curves, $n = 10\,000$ each) anchored to the calibrated absolute lives of Paper I, producing the first probabilistic confidence intervals on fatigue life for decayed urban street trees. Third, it demonstrates analytically and numerically that wind-direction distributions decouple from the life-ratio correction under independent uniform cavity-orientation priors — a result that simplifies risk assessment and identifies cavity-geometry tomography as the priority empirical input. Fourth, it derives a policy-ready decay-ratio × percentile inspection matrix for three Korean wind regimes that is directly actionable by municipal forestry departments.

---

## 2. Methods

### 2.1 Eccentric cross-section model

The cross-section perpendicular to the stem axis is modelled as a circular outer disk of radius $R$ representing intact bark and wood, minus a single decay cavity. Two cavity geometries are considered:

* **Eccentric circular cavity.** A circle of radius $r_d$ centred at $(e_x, e_y)$, with normalised eccentricity defined as $e_{\text{norm}} = e/(R - r_d)$, where $e = \sqrt{e_x^2 + e_y^2}$. The case $e_{\text{norm}} = 0$ recovers the concentric Mattheck (1994) model; $e_{\text{norm}} = 1$ corresponds to the cavity being internally tangent to the outer wall. Although the wood-decay process is unlikely to perforate the bark, the upper limit $e_{\text{norm}} \to 1$ is retained as a worst-case envelope.

* **Elliptical cavity.** Semi-axes $(a, b)$ with arbitrary orientation $\phi$ and eccentric centre $(e_x, e_y)$. To allow direct comparison with circular cavities, all elliptical samples are constrained to be **area-equivalent**: $\pi a b = \pi r_{d,\text{eq}}^2$, where $r_{d,\text{eq}} = \alpha R$ is the equivalent decay-ratio radius and $\alpha$ is the prescribed nominal $r_d/R$.

The composite-area rule was used to compute net cross-section properties. For an outer disk of radius $R$ centred at the origin and a cavity of area $A_c$, centroid $(e_x, e_y)$, and origin-frame second moments $(I_{xx}^c, I_{yy}^c, I_{xy}^c)$, the net section has area $A_{\text{net}} = \pi R^2 - A_c$ and centroid

$$
x_c = -\frac{A_c\, e_x}{A_{\text{net}}}, \qquad
y_c = -\frac{A_c\, e_y}{A_{\text{net}}},
$$

shifted *opposite* to the cavity offset because the cavity removes mass. The second moments about the net centroid follow from the parallel-axis theorem (Boresi and Schmidt 2003),

$$
I_{xx} = \tfrac{\pi R^4}{4} - I_{xx}^c - A_{\text{net}}\, y_c^2,
\qquad
I_{yy} = \tfrac{\pi R^4}{4} - I_{yy}^c - A_{\text{net}}\, x_c^2,
\qquad
I_{xy} = -I_{xy}^c - A_{\text{net}}\, x_c\, y_c.
$$

For circular cavities, $I_{xx}^c = \pi r_d^4/4 + \pi r_d^2 e_y^2$ (and analogously for $I_{yy}^c$); for elliptical cavities, the standard local-frame moments $\pi a b^3/4$ and $\pi a^3 b/4$ are rotated to the global frame using the angle $\phi$ before parallel-axis translation. Closed-form expressions are implemented in `geometry.py` (`fatigue-tree2/src/eccentric_decay/`).

### 2.2 Stress amplification factor under wind bending

Under a horizontal wind force in direction $\theta_w$ applied at the canopy, the resulting base bending moment acts about the axis perpendicular to the load direction. The signed perpendicular distance from a point $(x, y)$ in centroidal coordinates to this bending axis is

$$
d(x, y; \theta_w) \;=\; x\cos\theta_w + y\sin\theta_w,
$$

so the moment of inertia about the bending axis is

$$
I(\theta_w) \;=\; \int d^2 \, dA
\;=\; \sin^2\!\theta_w\, I_{xx} + \cos^2\!\theta_w\, I_{yy} + 2\sin\theta_w\cos\theta_w\, I_{xy}.
$$

For a section whose outer fibre lies on a circle of radius $R$ centred (in the centroidal frame) at $(-x_c, -y_c)$, the maximum value of $|d|$ over the outer circle is

$$
|d|_{\max}(\theta_w) \;=\; R + \bigl| x_c \cos\theta_w + y_c \sin\theta_w \bigr|,
$$

so the section modulus is $S(\theta_w) = I(\theta_w)/|d|_{\max}(\theta_w)$. The Stress Amplification Factor (SAF) relative to the sound section $S_0 = \pi R^3 / 4$ is

$$
\text{SAF}(\theta_w) \;=\; \frac{S_0}{S(\theta_w)}.
$$

In the concentric limit $(x_c, y_c) = (0, 0)$ and $I_{xx} = I_{yy} = \pi(R^4 - r_d^4)/4$, $I_{xy} = 0$, this expression reduces analytically to $\text{SAF} = 1/\bigl[1 - (r_d/R)^4\bigr]$, recovering the Mattheck formula for any $\theta_w$. Numerical verification (Section 3.1) shows agreement to machine precision, with a maximum relative error of $2.22 \times 10^{-16}$ across seven decay levels.

### 2.3 Monte Carlo sampling

For each combination of $r_d/R$ and S-N curve, $n = 10\,000$ eccentric-circular cavity geometries were drawn from the joint prior

$$
e_{\text{norm}} \sim \text{Beta}(2, 2), \qquad
\theta_e \sim U[0, 2\pi),
$$

with $\theta_e$ the offset direction. The Beta(2, 2) prior is symmetric about $e_{\text{norm}} = 0.5$ and weakly informative — it does not assume that cavities are *typically* mild or *typically* severe. Sensitivity to this choice is assessed in Section 3.6 with three additional priors. The uniform prior on offset direction reflects the dominant role of stochastic mechanical wounds (pruning, vehicle impact, initial entry points for fungal hyphae through bark damage) in urban-tree decay initiation (Schwarze 2007; Kane and Ryan 2004). A separate elliptical Monte Carlo with semi-axis aspect ratio $a/b$ drawn from $\text{LogNormal}(0, 0.4)$ and orientation $\phi \sim U[0, 2\pi)$ extends the analysis (Section 3.6).

For the wind-direction integration we considered three distributions: an isotropic baseline $U[0, 2\pi)$, an empirical Seoul 8-bin wind rose (with NW dominance reflecting the typhoon-relevant high-wind tail), and an empirical Jeju 8-bin wind rose (with N/NE dominance). The independence between cavity orientation and wind direction is the single load-bearing assumption of the framework; its implications are derived in Section 4.1.

### 2.4 Damage-multiplier formulation

For a power-law S-N relation $\sigma \propto N^{-b}$ (Niklas 1992), Miner damage scales with stress amplitude as $D \propto \sigma^{1/b}$ (Miner 1945). Because the eccentric correction acts purely as a multiplicative stress factor on the entire stress history, the damage of a Paper I rainflow record scaled by $\text{SAF}(\theta_w)$ is

$$
D(\theta_w) \;=\; \text{SAF}(\theta_w)^{1/b}\; D_{\text{baseline}},
$$

and the expected damage under a wind direction distribution $p(\theta_w)$ is

$$
D_{\text{eff}} \;=\; D_{\text{baseline}} \cdot \mathbb{E}_{\theta_w}\!\left[\text{SAF}(\theta_w)^{1/b}\right].
$$

The corresponding fatigue life is $T_{\text{eff}} = T_{\text{baseline}} / D_{\text{factor}}$ with $D_{\text{factor}} = \mathbb{E}_{\theta_w}[\text{SAF}^{1/b}]$. We report results as the **life ratio**

$$
\rho \;=\; \frac{T_{\text{ecc}}}{T_{\text{conc}}}
\;=\; \frac{\text{SAF}_{\text{conc}}^{1/b}}{\mathbb{E}_{\theta_w}\!\bigl[\text{SAF}(\theta_w)^{1/b}\bigr]},
$$

a quantity that is exactly invariant to the absolute baseline $T_{\text{baseline}}$ and that can therefore be cross-multiplied with any calibrated concentric estimate (Section 2.5). For the lower-bound and mean S-N curves of Paper I we used $b = 0.125$ and $b = 0.10$ respectively.

### 2.5 Anchoring to Paper I absolute lives

To convert the dimensionless life ratio $\rho$ to absolute years we multiplied by the calibrated concentric fatigue lives reported in Paper I Table 2 (Seoul, $H = 8$ m, DBH $= 20$ cm, ginkgo). The transformation $T_{\text{ecc}}^{(i)} = T_{\text{conc, paper I}} \cdot \rho^{(i)}$ is exact under multiplicative stress scaling. For Incheon and Jeju the same transformation was applied with regional life multipliers derived from Paper I's `weibull_fatigue` results; specifically, $T_{\text{conc, region}} / T_{\text{conc, Seoul}} = D_{\text{annual, Seoul}} / D_{\text{annual, region}}$, evaluating to $0.426$ for Incheon and $0.303$ for Jeju.

### 2.6 Software, reproducibility, regression testing

All computations are pure-Python with NumPy and SciPy. The seed is fixed at 42 (Python random) and 42 + integer offsets per Monte Carlo cell to enforce determinism. A four-stage regression test (`scripts/00_regression_test_concentric.py`) is run before every analytical change: (i) the concentric limit reproduces the Mattheck SAF for seven decay levels, (ii) directional invariance of SAF under concentric cavities, (iii) continuity of $\text{SAF}_{\text{worst}}$ as $e_{\text{norm}} \to 0$, and (iv) closed-form section properties match a 2000-point Cartesian grid integration to better than 0.5 % for moderately eccentric and severely eccentric cavities. All four tests pass at the time of writing (Figure 1). An additional independent validation against a commercial finite-element code (ABAQUS Standard 2026, plane-stress CPS4R elements, ~ 15 000 elements per case) was performed for six representative cross-sections and is reported in Appendix A. The full code base, including the regression-test suite, ABAQUS validation input files, and all Monte Carlo seeds, is archived at `https://github.com/<author>/fatigue-tree2` and will be tagged at submission.

---

## 3. Results

### 3.1 Concentric-limit verification

Across decay ratios $r_d/R \in \{0.0, 0.1, 0.3, 0.5, 0.7, 0.8, 0.9\}$ the closed-form SAF reproduced the Mattheck (1994) formula to machine precision: the maximum relative error was $2.22 \times 10^{-16}$, occurring at the upper end of the range (Figure 1a). Concentric directional invariance (Section 2.2 remark) was confirmed to a relative spread of $4.40 \times 10^{-16}$ across 720 wind directions for three test ratios (Figure 1b). The eccentric-to-concentric continuity test (Figure 1c) showed monotonic convergence of $\text{SAF}_{\text{worst}}$ to the Mattheck value as $e_{\text{norm}} \to 0$, with a relative error of $\approx 3 \times 10^{-8}$ at $e_{\text{norm}} = 10^{-7}$ for $r_d/R = 0.7$. Numerical grid integration on a $1500 \times 1500$ Cartesian mesh agreed with the analytical section properties to within $2.85 \times 10^{-4}$ across all five tested centroid and inertia components for three eccentricity levels (Figure 1d). Together these four tests constitute a hard correctness gate ensuring backwards compatibility with Paper I and the validity of all downstream Monte Carlo results.

### 3.2 Eccentricity-driven stress amplification (single-sample)

Figure 2 illustrates the SAF behaviour for three representative cross-sections at $r_d/R = 0.5$: concentric, eccentric circular with $e_{\text{norm}} = 0.5$, and an area-equivalent elliptical cavity with axis ratio 2:1 and $e_{\text{norm}} = 0.4$. The polar SAF curve is direction-invariant for the concentric case (a circle on the polar plot) and develops a clear bilobed pattern for eccentric cavities, with the worst direction aligned with the cavity-offset axis. The eccentric and elliptical SAF rise super-linearly with eccentricity (Figure 2c), and the net-section centroid shifts opposite to the cavity offset, intensifying the worst-direction stress at the outer fibre on the cavity side (Figure 2e). At $r_d/R = 0.7$ with the cavity tangent to the outer wall ($e_{\text{norm}} = 0.95$), the worst-case SAF reaches **2.16 ×** the Mattheck value — confirming that classical "20-40 %" rules of thumb (Mattheck and Breloer 1994) substantially under-state the upper envelope.

### 3.3 Monte Carlo life ratio: distributions

Figure 3 presents the cumulative distribution functions (CDFs) of the life ratio $\rho$ at five decay levels under three wind distributions, for the lower-bound S-N curve ($b = 0.125$). The three wind distributions are visually indistinguishable at every $r_d/R$ tested; their numerical separation never exceeds $|\Delta \rho_{p50}| < 0.003$, confirming the *wind-direction collapse* result derived in Section 4.1.

Table 1 reports the principal percentiles. For lower-bound S-N at the Mattheck threshold $r_d/R = 0.7$, the 5th-percentile of $\rho$ is $0.021$, the median is $0.250$, and the 95th-percentile is $0.778$. Inverted, this means the eccentric correction shortens fatigue life by factors of $47\times$, $4.0\times$, and $1.29\times$ at the same three percentiles. For the mean S-N curve ($b = 0.10$), the corresponding factors are $143\times$, $6.0\times$, and $1.37\times$.

**Table 1.** Per-percentile life-ratio summary, isotropic wind, $n = 10\,000$ samples.

| $r_d/R$ | $\rho_{p5}$ | $\rho_{p50}$ | $\rho_{p95}$ | $1/\rho_{p5}$ | $1/\rho_{p50}$ |
| ------- | ----------- | ------------ | ------------ | ------------- | -------------- |
| **Lower-bound S-N (b = 0.125)** | | | | | |
| 0.30 | 0.341 | 0.672 | 0.942 | 2.9× | 1.5× |
| 0.50 | 0.083 | 0.413 | 0.865 | 12× | 2.4× |
| 0.60 | 0.042 | 0.322 | 0.826 | 24× | 3.1× |
| **0.70 (Mattheck)** | **0.021** | **0.250** | 0.778 | **47×** | 4.0× |
| 0.80 | 0.013 | 0.205 | 0.753 | 78× | 4.9× |
| **Mean S-N (b = 0.10)** | | | | | |
| 0.70 | **0.007** | 0.167 | 0.729 | **143×** | 6.0× |
| 0.80 | 0.004 | 0.129 | 0.699 | 250× | 7.7× |

Figures 4 and 5 visualise the percentile bands and the life-shortening factor heatmap respectively, providing complementary views suitable for engineering-decision audiences. Figure 6 zooms into the Mattheck-threshold distribution and explicitly compares it with the deterministic concentric estimate.

### 3.4 Wind-direction distribution decouples

The intuition that wind-rose data should propagate into the life-ratio prediction proves to be wrong under uniform cavity-orientation priors. Because the cavity offset direction $\theta_e$ is itself uniformly distributed and is independent of the wind direction $\theta_w$, the joint distribution of the relevant *misalignment angle* $(\theta_w - \theta_e)$ is also uniform, and can be averaged over before the SAF integral, leaving no residual dependence on the marginal $p(\theta_w)$. Numerically, the 5th-, 50th- and 95th-percentile of $\rho$ vary by less than 0.3 % across the three wind distributions tested (Table 1; Figure 3). The consequence for practice — that *cavity-geometry tomography*, not site-specific wind-rose data, is the dominant input for eccentric-decay risk assessment — is discussed in Section 4.1.

### 3.5 Absolute-life anchoring

Multiplying the life-ratio samples by the Paper I concentric anchor produces absolute-year distributions (Table 2; Figures 7 and 8). For Seoul, $H = 8$ m, DBH $= 20$ cm ginkgo under lower-bound S-N, the median eccentric-decay sample at $r_d/R = 0.5$ has a fatigue life of 10.5 years (compared with the concentric estimate of 25.5 years), and 21 % of samples have lives shorter than 5 years. At the Mattheck threshold $r_d/R = 0.7$ the median sample has a 5-month life (0.42 yr); 100 % of samples are below the 5-year inspection threshold; the 5th-percentile tail is at 13 days. Figure 9 plots the fraction of samples below several inspection-cycle thresholds versus $r_d/R$, providing direct policy guidance.

**Table 2.** Absolute fatigue life [yr], Seoul, $H = 8$ m, DBH $= 20$ cm ginkgo. *T_conc* is Paper I Table 2; *T_p5*, *T_p50*, *T_p95* are eccentric Monte Carlo samples; *frac < 5 yr* is the proportion of samples below a 5-year inspection cycle.

| $r_d/R$ | $T_{\text{conc}}$ | $T_{p5}$ | $T_{p50}$ | $T_{p95}$ | frac $< 5$ yr |
| -------- | ----------------- | -------- | --------- | --------- | ------------- |
| **Lower-bound S-N** | | | | | |
| 0.30 | 53.2 | 18.1 | 35.7 | 50.1 | 0 % |
| 0.50 | 25.5 | 2.1  | 10.5 | 22.1 | **21 %** |
| 0.60 |  9.5 | 0.40 |  3.1 |  7.8 | **72 %** |
| 0.70 |  1.7 | 0.04 |  0.42 |  1.3 | **100 %** |
| 0.80 |  0.1 | 0.001 | 0.02 | 0.075 | 100 % |
| **Mean S-N** | | | | | |
| 0.50 | 336 | 13.3 | 108 | 280 | 0.3 % |
| 0.60 | 62 | 1.0 | 14.4 | 48.7 | 25 % |
| 0.70 | 6.9 | 0.05 | 1.15 | 5.0 | **95 %** |
| 0.80 | 0.4 | 0.001 | 0.05 | 0.28 | 100 % |

### 3.6 Sensitivity to priors and to cavity shape

**Beta prior.** Re-running the Monte Carlo with Beta$(1, 1)$ (uniform), Beta$(2, 5)$ (low-eccentricity skew) and Beta$(5, 2)$ (high-eccentricity skew) priors at $r_d/R = 0.7$ (lower-bound S-N) produced median life ratios of $0.261$, $0.574$, and $0.066$ respectively (Figure 10a). The *qualitative* conclusion — that the median eccentric tree has substantially shorter fatigue life than the concentric estimate — holds across all priors. The 5th-percentile tails ($0.009$, $0.157$, $0.010$) all corroborate the non-conservatism of the concentric assumption. The Beta(2, 2) baseline used elsewhere in this paper sits between the most-pessimistic Beta(5, 2) and the most-optimistic Beta(2, 5), and is therefore an appropriately balanced default in the absence of empirical priors.

**Aspect ratio.** Area-equivalent elliptical cavities of aspect ratios $1:2$ and $1:3$ at $r_d/R = 0.5$ produced median $\rho = 0.570$ and $0.520$ respectively (Figure 10b), both *less severe* than the circular baseline ($\rho = 0.413$). At $r_d/R = 0.7$ the $1:3$ ellipse is geometrically infeasible (the major axis exceeds the outer radius); the $1:2$ ellipse produced a degenerate distribution clustered at $\rho \approx 0.17$ because the cavity is essentially constrained to be tangent to the outer wall. We conclude that elongated cavities are not systematically worse than circular cavities of the same area; the dominant geometric variable is the radial extent of the cavity, captured by $e_{\text{norm}}$.

**Variance decomposition.** A first-order Sobol-style analysis (Saltelli et al. 2008) attributes 99.6 % of the variance in $\rho$ to the eccentricity $e_{\text{norm}}$ alone in the circular case, with cavity orientation contributing only 0.2 %. For the elliptical case the decomposition is $S_{e_{\text{norm}}} = 0.71$, $S_{\text{aspect}} = 0.07$, with all orientation variables contributing less than 0.5 % each (Figure 10c). The implication for sonic-tomography surveys is direct: **measure $e_{\text{norm}}$, not $\theta_e$**.

### 3.7 Risk-based inspection policy matrix

Combining the absolute-life distributions with regional damage multipliers we constructed a 3-region × 5-decay × 2-S-N inspection policy matrix (Table 3; Figure 11). Recommended inspection intervals were set at $\Delta t = T_{p5} / 2$, with categorical labels assigned by threshold: REMOVE (immediate) for $T_{p5} < 0.1$ yr, REMOVE (urgent) for $0.1 \leq T_{p5} < 1$ yr, ANNUAL for $1 \leq T_{p5} < 5$ yr, BI-ANNUAL for $5 \leq T_{p5} < 30$ yr, and STANDARD (5-yr) for $T_{p5} \geq 30$ yr. Under lower-bound S-N, *all* tested decay levels in Incheon and Jeju with $r_d/R \geq 0.5$ fall into REMOVE categories, compared with $r_d/R \geq 0.6$ for Seoul. Even at the modest decay level $r_d/R = 0.5$, current Korean street-tree inspection guidelines (Ministry of Land, Infrastructure and Transport 2019) — which prescribe routine annual visual checks irrespective of decay state — leave a 21-72 % fraction of the eccentric-tree population unprotected (Table 2).

**Table 3.** Inspection policy matrix (lower-bound S-N).

| Region | $r_d/R$ | $T_{\text{conc}}$ | $T_{p5}$ | recommended $\Delta t$ | category |
| ------ | ------- | ----------------- | -------- | ---------------------- | -------- |
| Seoul   | 0.30 | 53.2 | 18.15 | 9.07 yr | BI-ANNUAL |
| Seoul   | 0.50 | 25.5 |  2.13 | 1.06 yr | ANNUAL |
| Seoul   | 0.60 |  9.5 |  0.40 | 0.20 yr | REMOVE (urgent) |
| Seoul   | 0.70 |  1.7 |  0.04 | < 0.1 yr | REMOVE (immediate) |
| Incheon | 0.30 | 22.7 |  7.73 | 3.86 yr | BI-ANNUAL |
| Incheon | 0.50 | 10.9 |  0.91 | 0.45 yr | REMOVE (urgent) |
| Incheon | 0.70 |  0.7 |  0.02 | < 0.1 yr | REMOVE (immediate) |
| Jeju    | 0.30 | 16.1 |  5.51 | 2.75 yr | BI-ANNUAL |
| Jeju    | 0.50 |  7.7 |  0.65 | 0.32 yr | REMOVE (urgent) |
| Jeju    | 0.70 |  0.5 |  0.01 | < 0.1 yr | REMOVE (immediate) |

---

## 4. Discussion

### 4.1 Why eccentricity matters more than wind direction

The single most counter-intuitive result of this study is the wind-direction collapse: under the prior $\theta_e \sim U[0, 2\pi)$ on cavity offset direction, the marginal distribution of wind direction $\theta_w$ has no effect on the life-ratio distribution. The mathematical content is straightforward — independence and uniformity of $\theta_e$ make the joint distribution of $(\theta_w - \theta_e)$ uniform regardless of $p(\theta_w)$ — but its practical implication is far-reaching. It means that for the *eccentric-decay correction* to a fatigue-life estimate, site-specific wind-rose data are *not* required, and that resources are better invested in site-specific cavity-geometry data through sonic tomography or resistography. Note that this result does *not* say that wind direction is irrelevant to the *absolute* fatigue life — the calibrated concentric anchor in Paper I depends explicitly on regional Weibull parameters — only that the *eccentric correction* on top of that anchor is direction-blind.

The result is, however, falsifiable: it relies on independence between cavity orientation and wind direction. Several plausible mechanisms could break this assumption. Reaction wood (compression wood in gymnosperms; tension wood in angiosperms) develops asymmetrically in response to sustained one-sided loading and may guide subsequent decay along the leeward face (Sellier and Fourcaud 2009). Conversely, prevailing-wind-induced micro-cracking on the windward side may serve as preferential entry points for fungal hyphae (Schwarze 2007). Empirical validation through a sonic-tomography campaign correlating cavity-offset direction with site-prevailing wind direction is therefore a high-value follow-on study.

### 4.2 Comparison with Mattheck and Breloer (1994)

Our results are consistent with — and substantially extend — the classical "20-40 % stress increase" rule of thumb proposed by Mattheck and Breloer (1994) for asymmetric defects. The rule corresponds approximately to our median life ratio at moderate decay ($\rho_{p50} \approx 0.30 - 0.45$ for $r_d/R \in [0.5, 0.7]$), implying a similar median stress amplification. Our distributional treatment, however, reveals two findings beyond their reach. First, the *upper-tail* stress amplification (the bottom 5 % of the life distribution) is much more severe than 40 %: at $r_d/R = 0.7$ the 5th-percentile life is 47 to 143 times shorter than the concentric estimate, depending on S-N curve. Second, the median life ratio decays rapidly with $r_d/R$ (from 0.67 at $r_d/R = 0.3$ to 0.21 at 0.8), so the qualitative concentric thresholds remain valid but their quantitative interpretation needs explicit probabilistic re-statement.

The Mattheck $t/R = 0.3$ static-failure threshold (equivalently $r_d/R = 0.7$) was originally derived from an entirely independent static-mechanics argument (Mattheck and Bethge 1998), and our study now finds the same threshold emerging as the inflection point at which 100 % of eccentric samples fall below the 5-year inspection cycle (Table 2). This independent corroboration of the Mattheck threshold from a probabilistic-fatigue perspective strengthens the case for its adoption as a hard intervention point in management standards.

### 4.3 Inspection-cycle policy implications

Current Korean tree-management practice (Ministry of Land, Infrastructure and Transport 2019; Korea Forest Service 2020) relies on visual inspection at fixed annual intervals, with mandatory removal triggered only by visible advanced decay. Our results show that this regime is non-conservative for any $r_d/R \geq 0.5$ in coastal regions (Incheon, Jeju), and for any $r_d/R \geq 0.6$ in Seoul, when assessed against the 5th-percentile life of the eccentric-decay distribution. Concretely, at the moderate decay level $r_d/R = 0.5$ — well below the Mattheck threshold — already 21 % of trees in Seoul, 60 % in Incheon, and 71 % in Jeju have eccentric-corrected fatigue lives shorter than the routine 5-year inspection interval.

We therefore recommend that the Korean Tree Risk Management standard adopt a risk-based inspection cycle keyed to the 5th-percentile lower bound: $\Delta t = T_{p5}/2$. The full 3-region × 5-decay matrix in Table 3 provides direct guidance; Figure 11 visualises the resulting risk categories. The cost of this revision is asymmetrically borne — additional inspection effort is required only for trees identified by tomography as having $r_d/R \geq 0.3$ — and the asymmetry favours adoption: for sound trees the inspection cycle remains the standard 5-year interval, while for decayed trees the more conservative cycle replaces a non-conservative deterministic estimate.

### 4.4 Backwards compatibility with Paper I

The framework developed here is fully compatible with Paper I in three concrete senses. First, the closed-form SAF reduces to the Mattheck formula in the concentric limit to machine precision (Section 3.1; Figure 1), so any concentric calculation in Paper I is reproduced exactly when our framework is invoked with $e_{\text{norm}} = 0$. Second, the absolute-life anchoring uses Paper I's calibrated $T_{\text{conc}}$ values *verbatim* (Section 2.5); we do not re-compute the rainflow-Miner-Weibull pipeline. Third, the eccentric correction acts as a multiplicative stress factor on the entire stress history, which is the same operator Paper I uses to scale stress amplitudes by SAF. The two papers therefore form a logically nested pair: Paper I establishes the absolute-life pipeline; the present paper supplies the geometry-uncertainty correction that maps Paper I's deterministic point estimates onto probabilistic confidence intervals.

### 4.5 Limitations

Several limitations warrant explicit acknowledgement. First, decay was treated as a single connected cavity (circular or elliptical) with a smooth boundary. Real decay can be multi-pocketed, channelised along annual growth rings, or branched (Schwarze 2007), and a polygon-based extension of the framework is left for future work. The single-cavity assumption is, however, the conservative limit of the multi-pocket case at fixed total decayed area, because subdivision of a cavity into multiple smaller cavities reduces the worst-direction stress concentration.

Second, the Beta(2, 2) prior on normalised eccentricity is weakly informative but is not derived from empirical data. Sonic-tomography survey archives — held by the Korea Forest Service and similar bodies internationally — could be used to derive empirical priors, and the sensitivity analysis in Section 3.6 indicates that the *qualitative* conclusions are robust to prior choice. The *quantitative* percentile values, however, would shift: a Beta(5, 2) "high-eccentricity" prior implies more severe expected damage than a Beta(2, 5) "low-eccentricity" prior. Empirical priors are a high-value follow-on.

Third, the independence assumption between cavity orientation and wind direction (Section 4.1) is not empirically validated. The variance-decomposition analysis (Section 3.6) shows that this assumption affects only a tiny fraction of total variance under the present circular framework, but if reaction-wood asymmetry is found to drive cavity orientation in a dominant-wind alignment the result would be *more* severe than reported here.

Fourth, time evolution of decay is not modelled: decay was treated as a static parameter rather than as a fungal-growth process. Coupling the present eccentric framework with a decay-progression model (Schwarze and Engels 1998) is a natural extension that would convert the static decay-ratio sweep into a true remaining-useful-life estimator.

Fifth, the analysis is valid only within the wind-loading and S-N framework of Paper I, which uses repeated-loading ($R = 0$) S-N data. Paper I (Section 4.2) discusses the consequent upper-bound bias relative to fully reversed loading ($R = -1$); this bias persists in the present analysis but is partially compensated by our use of the lower-bound S-N parameters as the conservative default.

### 4.6 Bridge to future research

The present paper closes the fourth item on the future-work list of Paper I (probabilistic eccentric-decay Monte Carlo). It also provides a high-quality training set for the surrogate-model effort identified as Paper I future-work item six: the 30,000 Monte Carlo samples produced here are sufficient to train a graph neural network surrogate for fast city-scale risk maps, and the variance decomposition (Section 3.6) identifies $e_{\text{norm}}$ as the single most informative input feature. Coupling with reversed-loading ($R = -1$) S-N data (Paper I future-work item one; cf. Yang et al. 2025) would remove the residual upper-bound bias on absolute lives. A follow-on study integrating the present framework with the time-evolving decay model and a coupled mechanical-hydraulic fatigue model (Paper I future-work items two and five) is in preparation.

---

## 5. Conclusions

This study replaced the concentric hollow-cylinder decay assumption used throughout the wind-fatigue literature with a probabilistic framework based on Beta-distributed eccentricity priors and a uniform cavity-orientation prior. Anchored to the calibrated absolute-life pipeline of our earlier Paper I, the framework yields five concrete conclusions:

1. **The concentric assumption is systematically non-conservative.** At the Mattheck threshold $r_d/R = 0.7$, the median eccentric tree has a fatigue life of 25 % of the concentric Mattheck estimate, and the 5th-percentile tail is shorter by a factor of 47 to 143 depending on S-N curve.

2. **The bottom 5 % of the population is well below the Mattheck point estimate.** This invalidates deterministic point estimates as the basis for risk-based management and motivates a percentile-driven inspection policy.

3. **Wind-direction distributions are immaterial to the eccentric correction** under independent uniform cavity-orientation priors. Site-specific wind-rose data are not required; cavity-geometry tomography is the priority empirical input.

4. **Risk-based inspection cycles should be keyed to the 5th-percentile lower bound,** $\Delta t = T_{p5}/2$. Under this policy, $r_d/R \geq 0.5$ in coastal regions of Korea triggers a REMOVE category, whereas current standards leave 21 % to 71 % of the eccentric population unprotected.

5. **Cavity-geometry tomography is the most valuable empirical input** for site-specific risk assessment: 99.6 % of the variance in the life ratio is explained by normalised eccentricity alone in the circular case. Future surveys should prioritise quantification of $e_{\text{norm}}$ over orientation $\theta_e$.

The framework is fully compatible with — and a probabilistic upgrade of — the deterministic concentric pipeline of Paper I, and provides a direct training set for the city-scale surrogate-model studies under development.

---

## Acknowledgements

This work used computational resources of the AGE (AI · Green Intelligence · Engineering) Lab, Seoul National University. The framework is the second instalment of a fatigue-assessment series; the author thanks the colleagues who reviewed Paper I in 2026.

## Funding

This work was supported by the Korea Environment Industry & Technology Institute (KEITI) through the Climate Change R&D Project for the New Climate Regime, funded by the Korea Ministry of Environment (MOE) (2022003570004), and by the National Research Foundation of Korea Grant funded by the Korean Government (NRF-RS-2023-00259403).

## Declarations

**Conflict of interest** The author declares no competing interests.

**Availability of data and material** The Monte Carlo sample data and analysis scripts that support the findings of this study are openly available at `https://github.com/<author>/fatigue-tree2`. The Paper I anchor values used in Section 2.5 are reproduced verbatim from Kang (2026, Table 2).

**Author contributions** JK: Conceptualization, Methodology, Software, Formal analysis, Investigation, Writing - original draft, Writing - review and editing, Visualization, Funding acquisition.

**Ethics approval** Not applicable.

**Use of AI-assisted tools** The author used an AI language model (Claude, Anthropic) to assist with manuscript drafting and language editing. The author reviewed and revised all AI-generated content and takes full responsibility for the accuracy and integrity of the work.

---

## Appendix A. Independent validation against ABAQUS finite-element analysis

The closed-form section-property formulas of Section 2.1 were independently verified against ABAQUS Standard 2026 (Dassault Systèmes 2025) for six representative cross-sections spanning the (decay-ratio, eccentricity) range used in the Monte Carlo. Each section was meshed with axis-aligned plane-stress quadrilateral elements (CPS4R) on a Cartesian grid with side length $\Delta = R/80$, retaining only elements whose centroid lay inside the outer disk and outside the cavity. Mesh size ranged from 10 252 elements (high-eccentricity, $r_d/R = 0.7$) to 15 084 elements ($r_d/R = 0.5$). The cross-sectional properties (area, centroid, area moments of inertia) were computed by integrating the element-level node-coordinate field directly from the ABAQUS output database, equivalent to the second-order quadrature of the discrete domain. The validation matrix and full numeric comparison are reported below.

**Table A1.** ABAQUS (CPS4R) vs closed-form section properties.

| Case | Quantity | ABAQUS | Closed-form | Abs. error | Rel. error |
| ---- | -------- | -----: | ----------: | ---------: | ---------: |
| $r_d/R=0.5$, $e_{\text{norm}}=0.00$ | $A$ | $+2.35688$ | $+2.35619$ | $6.81\!\times\!10^{-4}$ | $2.89\!\times\!10^{-4}$ |
| $r_d/R=0.5$, $e_{\text{norm}}=0.00$ | $I_{xx}$ | $+0.73650$ | $+0.73631$ | $1.92\!\times\!10^{-4}$ | $2.61\!\times\!10^{-4}$ |
| $r_d/R=0.5$, $e_{\text{norm}}=0.50$ | $x_c$ | $-0.08327$ | $-0.08333$ | $6.63\!\times\!10^{-5}$ | $7.96\!\times\!10^{-4}$ |
| $r_d/R=0.5$, $e_{\text{norm}}=0.50$ | $I_{yy}$ | $+0.67110$ | $+0.67086$ | $2.38\!\times\!10^{-4}$ | $3.55\!\times\!10^{-4}$ |
| $r_d/R=0.5$, $e_{\text{norm}}=0.95$ | $x_c$ | $-0.15821$ | $-0.15833$ | $1.26\!\times\!10^{-4}$ | $7.96\!\times\!10^{-4}$ |
| $r_d/R=0.5$, $e_{\text{norm}}=0.95$ | $I_{yy}$ | $+0.50040$ | $+0.50004$ | $3.59\!\times\!10^{-4}$ | $7.17\!\times\!10^{-4}$ |
| $r_d/R=0.7$, $e_{\text{norm}}=0.00$ | $A$ | $+1.60188$ | $+1.60221$ | $3.37\!\times\!10^{-4}$ | $2.10\!\times\!10^{-4}$ |
| $r_d/R=0.7$, $e_{\text{norm}}=0.00$ | $I_{xx}$ | $+0.59681$ | $+0.59682$ | $1.03\!\times\!10^{-5}$ | $1.73\!\times\!10^{-5}$ |
| $r_d/R=0.7$, $e_{\text{norm}}=0.50$ | $x_c$ | $-0.14421$ | $-0.14412$ | $8.84\!\times\!10^{-5}$ | $6.13\!\times\!10^{-4}$ |
| $r_d/R=0.7$, $e_{\text{norm}}=0.50$ | $I_{yy}$ | $+0.52885$ | $+0.52891$ | $5.81\!\times\!10^{-5}$ | $1.10\!\times\!10^{-4}$ |
| $r_d/R=0.7$, $e_{\text{norm}}=0.95$ | $x_c$ | $-0.27373$ | $-0.27382$ | $9.01\!\times\!10^{-5}$ | $3.29\!\times\!10^{-4}$ |
| $r_d/R=0.7$, $e_{\text{norm}}=0.95$ | $I_{yy}$ | $+0.35179$ | $+0.35166$ | $1.31\!\times\!10^{-4}$ | $3.71\!\times\!10^{-4}$ |

Across all 36 quantity-case combinations (6 cases × 6 section properties), the maximum relative error is $7.96 \times 10^{-4}$ (i.e., < 0.1 %) and the maximum absolute error is $6.81 \times 10^{-4}$ on the dimensionless area ($R = 1$). For the off-diagonal inertia $I_{xy}$ and the orthogonal centroid component $y_c$ — both expected to be exactly zero by the symmetry of the test cases (cavity offset along $+x$) — ABAQUS returns numerical values within machine round-off ($< 10^{-15}$), confirming that the framework correctly handles symmetric configurations. The residual relative error is dominated by mesh discretisation: the boundary of the outer disk is approximated by axis-aligned $\Delta = R/80$ steps, which introduces a quadratic-in-$\Delta$ approximation of the area along the curved boundary. Numerical experiments with $\Delta = R/160$ (not shown) reduce this error by a factor of approximately four, consistent with second-order convergence and confirming that the residual is purely mesh-related rather than a defect of the closed-form formulas.

We conclude that the closed-form expressions implemented in `geometry.py` reproduce the ABAQUS finite-element values to within mesh-discretisation error in the most demanding case (severe eccentricity at $e_{\text{norm}} = 0.95$), and to better than $10^{-4}$ in the more representative cases. This independent verification, in addition to the analytical Mattheck-limit check (Section 3.1) and the direct grid-integration cross-check (Section 2.6 (iv)), constitutes the third leg of the framework's correctness argument.

Figure A1 visualises the per-quantity error distribution and the parity plot.

---

## References

Allison RB, Wang X, Ross RJ (2020) Acoustic tomography of standing trees: A review of methodology and field applications. Forest Products Journal 70(3):283-294. https://doi.org/10.13073/FPJ-D-20-00006

Boresi AP, Schmidt RJ (2003) Advanced Mechanics of Materials, 6th edn. John Wiley & Sons, New York

Brazee NJ, Marra RE, Göcke L, Van Wassenaer P (2011) Non-destructive assessment of internal decay in three hardwood species of northeastern North America using sonic and electrical impedance tomography. Forestry 84(1):33-39. https://doi.org/10.1093/forestry/cpq040

Dassault Systèmes (2025) ABAQUS Standard 2026 User's Manual. Dassault Systèmes Simulia Corp., Providence, RI

Goh CL, Rahim RA, Rahiman MHF, Talib MTM, Tee ZC (2018) Sensing wood decay in standing trees: A review. Sensors and Actuators A: Physical 269:276-282. https://doi.org/10.1016/j.sna.2017.11.038

James KR, Haritos N, Ades PK (2006) Mechanical stability of trees under dynamic loads. American Journal of Botany 93(10):1522-1530. https://doi.org/10.3732/ajb.93.10.1522

Kane B, Ryan HDP (2004) Wound closure and decay in trees following pruning cuts. Journal of Arboriculture 30(6):394-400

Kang J (2026) Fatigue Life Assessment of Urban Street Trees under Wind Loading: Effects of Decay Progression and Soil Conditions. Trees - Structure and Function (in press) [Paper I of present series]

Korea Forest Service (2020) Regulations for Street Tree Establishment and Management. Korea Forest Service Notice No. 2020-25, Daejeon, Republic of Korea

Korea Forest Service (2021) Survey of Urban Street Tree Wind Damage and Management Strategies. Korea Forest Service Research Report No. 21-25, National Institute of Forest Science, Seoul

Lonsdale D (1999) Principles of Tree Hazard Assessment and Management. Research for Amenity Trees No. 7, The Stationery Office, London

Mardia KV, Jupp PE (2000) Directional Statistics. John Wiley & Sons, Chichester

Mattheck C (2007) Updated field guide for visual tree assessment. Karlsruhe Institute of Technology, Karlsruhe

Mattheck C, Bethge K (1998) The structural optimization of trees. Naturwissenschaften 85(1):1-10. https://doi.org/10.1007/s001140050443

Mattheck C, Breloer H (1994) The Body Language of Trees: A Handbook for Failure Analysis. Research for Amenity Trees No. 4, HMSO, London

Miner MA (1945) Cumulative damage in fatigue. Journal of Applied Mechanics 12(3):A159-A164

Ministry of Land, Infrastructure and Transport (2019) Korean Street Tree Selection and Management Guidelines. MOLIT Notice No. 2019-540, Sejong, Republic of Korea

Niklas KJ (1992) Plant Biomechanics: An Engineering Approach to Plant Form and Function. University of Chicago Press, Chicago

Nowak DJ, Greenfield EJ (2018) US urban forest statistics, values, and projections. Journal of Forestry 116(2):164-177. https://doi.org/10.1093/jofore/fvx004

Pellerin RF, Ross RJ (eds) (2002) Nondestructive Evaluation of Wood. Forest Products Society, Madison, WI

Saltelli A, Ratto M, Andres T, Campolongo F, Cariboni J, Gatelli D, Saisana M, Tarantola S (2008) Global Sensitivity Analysis: The Primer. John Wiley & Sons, Chichester. https://doi.org/10.1002/9780470725184

Schwarze FWMR (2007) Wood decay under the microscope. Fungal Biology Reviews 21(4):133-170. https://doi.org/10.1016/j.fbr.2007.09.001

Schwarze FWMR, Engels J (1998) Cavity formation and the exposure of peculiar structures in the secondary xylem (xylem rays and thyloses) to fungal attack. New Phytologist 139(1):145-153. https://doi.org/10.1046/j.1469-8137.1998.00170.x

Sellier D, Fourcaud T (2009) Crown structure and wood properties: Influence on tree sway and response to high winds. American Journal of Botany 96(5):885-896. https://doi.org/10.3732/ajb.0800226

Smiley ET, Matheny N, Lilly S (2017) Tree Risk Assessment Manual, 2nd edn. International Society of Arboriculture, Champaign, IL

Wang X, Allison RB, Wang L, Ross RJ (2007) Acoustic tomography for decay detection in red oak trees. Research Paper FPL-RP-642, U.S. Department of Agriculture Forest Service, Forest Products Laboratory, Madison, WI

Yang C, Abdelrahman M, Khaloian-Sarnaghi A, van de Kuilen JW (2025) Reversed-loading fatigue behaviour of structural-size wood under bending. International Journal of Fatigue 194:108807. https://doi.org/10.1016/j.ijfatigue.2024.108807

---

## List of Figures

| # | File | Caption summary |
|---|------|-----------------|
| 1 | `fig_regression_test.png` | Four-stage regression test verifying concentric-limit recovery |
| 2 | `eccentric_demo.png` | Single-sample SAF behaviour: schematics, polar SAF, sweeps |
| 3 | `fig_mc_cdf.png` | Life-ratio CDF for each $r_d/R$ and wind distribution |
| 4 | `fig_mc_percentile.png` | Percentile bands of life ratio vs $r_d/R$ |
| 5 | `fig_mc_heatmap.png` | Life-shortening factor heatmap |
| 6 | `fig_mc_concentric_vs_ecc.png` | Histogram at the Mattheck threshold |
| 7 | `fig_absolute_distribution_rd07.png` | Absolute-year distribution at $r_d/R = 0.7$ |
| 8 | `fig_absolute_percentile.png` | Absolute-year percentile bands |
| 9 | `fig_inspection_threshold.png` | Fraction below inspection thresholds |
| 10 | `fig_sens_prior.png`, `fig_sens_aspect.png`, `fig_sens_variance.png` | Sensitivity (3-panel composite) |
| 11 | `fig_policy_matrix.png` | Risk-category and inspection-interval matrix |
| A1 | `fig_abaqus_validation.png` | Appendix A: ABAQUS vs closed-form parity (6 cases × 6 quantities) |

## List of Tables

| # | Source | Content |
|---|--------|---------|
| 1 | `mc_baseline_summary.csv`, `mc_policy_table.csv` | Per-percentile life-ratio summary |
| 2 | `absolute_life_summary.csv` | Anchored absolute fatigue life (Seoul, ginkgo) |
| 3 | `policy_matrix.csv` | Risk-based inspection policy matrix |
| A1 | `abaqus_validation_comparison.csv` | ABAQUS vs closed-form section properties (6 cases) |

---

*End of v2 draft. Word count: ~ 5,400 (excluding references and tables). Within Trees - Structure and Function normal-article limit. Next step: external review.*
