# Probabilistic Eccentric Decay Geometry: Monte Carlo Confidence Intervals on the Wind-Fatigue Life of Urban Street Trees

**Junsuk Kang** ᵃᵇᶜᵈ \*

ᵃ Interdisciplinary Program in Landscape Architecture, Seoul National University, Seoul 08826, Republic of Korea
ᵇ Transdisciplinary Program in Smart City Global Convergence, Seoul National University, Seoul 08826, Republic of Korea
ᶜ Department of Landscape Architecture and Rural Systems Engineering, Seoul National University, Seoul 08826, Republic of Korea
ᵈ Research Institute of Agriculture and Life Sciences, Seoul National University, Seoul 08826, Republic of Korea

\* Correspondence: junkang@snu.ac.kr; Tel.: +82-2-880-4872

> **Status:** Working draft v1 (2026-05-09). Methods §2 and Results §3 are written in full. Introduction §1 and Discussion §4 follow the outline in `paper2_outline.md` and will be expanded once the present text is reviewed.

---

## Abstract (placeholder, ≤ 250 words)

The classical concentric hollow-cylinder model of stem decay (Mattheck 1994), reproduced in the companion paper (Kang 2026, hereafter Paper I), assumes that wood-decay cavities are circular and centred on the stem axis. Field tomography surveys consistently show, however, that real cavities are eccentric, irregular, and often elongated. This paper relaxes the concentric assumption and quantifies, through a Monte Carlo framework, how cavity-geometry uncertainty propagates into uncertainty on wind-fatigue life. Closed-form section properties are derived for circular and elliptical cavities of arbitrary eccentricity, and a vectorisable stress-amplification function is verified to reproduce the Mattheck formula to machine precision in the concentric limit. Ten thousand cavity geometries are sampled per scenario from a Beta(2, 2) prior on normalised eccentricity and a uniform prior on cavity orientation; the resulting per-sample fatigue lives are anchored to Paper I's calibrated absolute-life pipeline. The median fatigue life at the Mattheck threshold $r_d/R = 0.7$ is found to be 25 % of the concentric estimate (lower-bound S-N), and the 5th-percentile tail is shorter by a factor of 47. Wind-direction distributions (isotropic, Seoul empirical rose, Jeju empirical rose) collapse out under independent uniform cavity-orientation priors. A risk-based inspection policy keyed to the 5th-percentile lower bound rather than the deterministic concentric estimate is proposed and tabulated for three Korean wind regimes.

**Keywords:** wind-fatigue, eccentric decay, Monte Carlo, urban street trees, probabilistic risk assessment, Mattheck threshold

---

## 1. Introduction (outline only)

*See `paper2_outline.md §1` for the outline. Full text to be drafted in v2.*

Briefly: motivation (urban tree windthrow risk; recap of Paper I); the limits of the concentric assumption; sonic-tomography evidence for eccentric decay; research questions; contributions.

---

## 2. Methods

### 2.1 Eccentric cross-section model

The cross-section perpendicular to the stem axis is modelled as a circular outer disk of radius $R$ representing intact bark and wood, minus a single decay cavity. Two cavity geometries are considered:

* **Eccentric circular cavity.** A circle of radius $r_d$ centred at $(e_x, e_y)$, with normalised eccentricity defined as $e_{\text{norm}} = e/(R - r_d)$, where $e = \sqrt{e_x^2 + e_y^2}$. The case $e_{\text{norm}} = 0$ recovers the concentric Mattheck (1994) model; $e_{\text{norm}} = 1$ corresponds to the cavity being internally tangent to the outer wall. Although physically the wood-decay process is unlikely to perforate the bark, the upper limit $e_{\text{norm}} \to 1$ is retained as a worst-case envelope.

* **Elliptical cavity.** Semi-axes $(a, b)$ with arbitrary orientation $\phi$ and eccentric centre $(e_x, e_y)$. To allow direct comparison with circular cavities, all elliptical samples are constrained to be **area-equivalent**: $\pi a b = \pi r_{d,\text{eq}}^2$, where $r_{d,\text{eq}} = \alpha R$ is the equivalent decay-ratio radius and $\alpha$ is the prescribed nominal $r_d/R$.

The composite-area rule was used to compute net cross-section properties. For an outer disk of radius $R$ centred at the origin and a cavity of area $A_c$, centroid $(e_x, e_y)$, and origin-frame second moments $(I_{xx}^c, I_{yy}^c, I_{xy}^c)$, the net section has area $A_{\text{net}} = \pi R^2 - A_c$ and centroid

$$
x_c = -\frac{A_c\, e_x}{A_{\text{net}}}, \qquad
y_c = -\frac{A_c\, e_y}{A_{\text{net}}},
$$

shifted *opposite* to the cavity offset because the cavity removes mass. The second moments about the net centroid follow from the parallel-axis theorem,

$$
I_{xx} = \tfrac{\pi R^4}{4} - I_{xx}^c - A_{\text{net}}\, y_c^2, \qquad
I_{yy} = \tfrac{\pi R^4}{4} - I_{yy}^c - A_{\text{net}}\, x_c^2, \qquad
I_{xy} = -I_{xy}^c - A_{\text{net}}\, x_c\, y_c.
$$

For circular cavities, $I_{xx}^c = I_{yy}^c = \pi r_d^4 / 4 + \pi r_d^2 e_y^2$ (and analogously for $I_{yy}^c$); for elliptical cavities, the standard local-frame moments $\pi a b^3 / 4$ and $\pi a^3 b / 4$ are rotated to the global frame using the angle $\phi$ before parallel-axis translation. Closed-form expressions are implemented in `geometry.py` (`fatigue-tree2/src/eccentric_decay/`).

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

In the concentric limit $(x_c, y_c) = (0, 0)$ and $I_{xx} = I_{yy} = \pi(R^4 - r_d^4)/4$, $I_{xy} = 0$, this expression reduces analytically to $\text{SAF} = 1/\bigl[1 - (r_d/R)^4\bigr]$, recovering the Mattheck formula for any $\theta_w$. Numerical verification (Section 3.1) shows agreement to machine precision (relative error $< 10^{-15}$).

### 2.3 Monte Carlo sampling

For each combination of $r_d/R$ and S-N curve, $n = 10\,000$ eccentric-circular cavity geometries were drawn from the joint prior

$$
e_{\text{norm}} \sim \text{Beta}(2, 2), \qquad
\theta_e \sim U[0, 2\pi),
$$

with $\theta_e$ the offset direction. The Beta(2, 2) prior is symmetric about $e_{\text{norm}} = 0.5$ and weakly informative — it does not assume that cavities are *typically* mild or *typically* severe. Sensitivity to this choice is assessed in Section 3.6. The uniform prior on offset direction reflects the dominant role of stochastic mechanical wounds (pruning, vehicle impact) in urban-tree decay initiation. A separate elliptical Monte Carlo with semi-axis aspect ratio $a/b$ drawn from $\text{LogNormal}(0, 0.4)$ and orientation $\phi \sim U[0, 2\pi)$ extends the analysis (Section 3.6).

For the wind-direction integration we considered three distributions: an isotropic baseline $U[0, 2\pi)$, an empirical Seoul 8-bin wind rose, and an empirical Jeju 8-bin wind rose. The independence between cavity orientation and wind direction is the single load-bearing assumption of the framework; its consequences are discussed in Section 4.1.

### 2.4 Damage-multiplier formulation

For a power-law S-N relation $\sigma \propto N^{-b}$, Miner damage scales with stress amplitude as $D \propto \sigma^{1/b}$. Because the eccentric correction acts purely as a multiplicative stress factor on the entire stress history, the damage of a Paper I rainflow record scaled by $\text{SAF}(\theta_w)$ is

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

All computations are pure-Python with NumPy and SciPy. The seed is fixed at 42 (Python random) and 42 + integer offsets per Monte Carlo cell to enforce determinism. A four-stage regression test (`scripts/00_regression_test_concentric.py`) is run before every analytical change: (i) the concentric limit reproduces the Mattheck SAF for seven decay levels, (ii) directional invariance of SAF under concentric cavities, (iii) continuity of $\text{SAF}_{\text{worst}}$ as $e_{\text{norm}} \to 0$, and (iv) closed-form section properties match a 2000-point Cartesian grid integration to better than 0.5 % for moderately eccentric and severely eccentric cavities. All four tests pass at the time of writing. The full code base, including the regression-test suite and all Monte Carlo seeds, is archived at `https://github.com/<author>/fatigue-tree2` and will be tagged at submission.

---

## 3. Results

### 3.1 Concentric-limit verification

Across decay ratios $r_d/R \in \{0.0, 0.1, 0.3, 0.5, 0.7, 0.8, 0.9\}$ the closed-form SAF reproduced the Mattheck (1994) formula to machine precision: the maximum relative error was $1.69 \times 10^{-16}$, occurring at $r_d/R = 0.7$. Concentric directional invariance (Section 2.2 remark) was confirmed to a relative spread of $4.40 \times 10^{-16}$ across 720 wind directions. Numerical grid integration agreed with the analytical section properties to within 0.5 % for all six tested centroid and inertia components across three eccentricity levels. These results constitute a hard correctness gate ensuring backwards compatibility with Paper I.

### 3.2 Eccentricity-driven stress amplification (single-sample)

Figure 2 illustrates the SAF behaviour for three representative cross-sections at $r_d/R = 0.5$: concentric, eccentric circular with $e_{\text{norm}} = 0.5$, and an area-equivalent elliptical cavity with axis ratio 2:1 and $e_{\text{norm}} = 0.4$. The polar SAF curve is direction-invariant for the concentric case (a circle on the polar plot) and develops a clear bilobed pattern for eccentric cavities, with the worst direction aligned with the cavity-offset axis. At $r_d/R = 0.7$ with the cavity tangent to the outer wall ($e_{\text{norm}} = 0.95$), the worst-case SAF reaches **2.16 ×** the Mattheck value (Table A.1), confirming that classical "20-40 %" rules of thumb (Mattheck and Breloer 1994) substantially under-state the upper envelope.

### 3.3 Monte Carlo life ratio: distributions

Figure 3 presents the cumulative distribution functions (CDFs) of the life ratio $\rho$ at five decay levels under three wind distributions, for the lower-bound S-N curve ($b = 0.125$). The three wind distributions are visually indistinguishable; their numerical separation never exceeds $|\Delta \rho_{p50}| < 0.003$, confirming the *wind-direction collapse* result (Section 4.1).

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

### 3.4 Wind-direction distribution decouples

The intuition that wind-rose data should propagate into the life-ratio prediction proves to be wrong under uniform cavity-orientation priors. Because the cavity offset direction $\theta_e$ is itself uniformly distributed and is independent of the wind direction $\theta_w$, the joint distribution of the relevant *misalignment angle* $(\theta_w - \theta_e)$ is also uniform and can be averaged over before the SAF integral, leaving no residual dependence on the marginal $p(\theta_w)$. Numerically, the 5th-, 50th- and 95th-percentile of $\rho$ vary by less than 0.3 % across the three wind distributions tested (Table 1; Figure 3). The consequence for practice — that *cavity-geometry tomography*, not site-specific wind-rose data, is the dominant input for eccentric-decay risk assessment — is discussed in Section 4.1.

### 3.5 Absolute-life anchoring

Multiplying the life-ratio samples by the Paper I concentric anchor produces absolute-year distributions (Table 2; Figure 7). For Seoul, $H = 8$ m, DBH $= 20$ cm ginkgo under lower-bound S-N, the median eccentric-decay sample at $r_d/R = 0.5$ has a fatigue life of 10.5 years (compared with the concentric estimate of 25.5 years), and 21 % of samples have lives shorter than 5 years. At the Mattheck threshold $r_d/R = 0.7$ the median sample has a 5-month life (0.42 yr); 100 % of samples are below the 5-year inspection threshold; the 5th-percentile tail is at 13 days.

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

**Beta prior.** Re-running the Monte Carlo with Beta$(1, 1)$ (uniform), Beta$(2, 5)$ (low-eccentricity skew) and Beta$(5, 2)$ (high-eccentricity skew) priors at $r_d/R = 0.7$ (lower-bound S-N) produced median life ratios of $0.261$, $0.574$, and $0.066$ respectively (Figure 8a). The qualitative conclusion — that the median eccentric tree has substantially shorter fatigue life than the concentric estimate — holds across all priors. The 5th-percentile tails ($0.009$, $0.157$, $0.010$) all corroborate the non-conservatism of the concentric assumption.

**Aspect ratio.** Area-equivalent elliptical cavities of aspect ratios $1:2$ and $1:3$ at $r_d/R = 0.5$ produced median $\rho = 0.570$ and $0.520$ respectively (Figure 8b), both *less severe* than the circular baseline ($\rho = 0.413$). At $r_d/R = 0.7$ the $1:3$ ellipse is geometrically infeasible (the major axis exceeds the outer radius); the $1:2$ ellipse produced a degenerate distribution clustered at $\rho \approx 0.17$ because the cavity is essentially constrained to be tangent to the outer wall. We conclude that elongated cavities are not systematically worse than circular cavities of the same area; the dominant geometric variable is the radial extent of the cavity, captured by $e_{\text{norm}}$.

**Variance decomposition.** A first-order Sobol-style analysis (Section 2.3 protocol) attributes 99.6 % of the variance in $\rho$ to the eccentricity $e_{\text{norm}}$ alone in the circular case, with cavity orientation contributing only 0.2 %. For the elliptical case the decomposition is $S_{e_{\text{norm}}} = 0.71$, $S_{\text{aspect}} = 0.07$, with all orientation variables contributing less than 0.5 % each (Figure 8c). The implication for sonic-tomography surveys is direct: **measure $e_{\text{norm}}$, not $\theta_e$**.

### 3.7 Risk-based inspection policy matrix

Combining the absolute-life distributions with regional damage multipliers we constructed a 3-region × 5-decay × 2-S-N inspection policy matrix (Table 3; Figure 9). Recommended inspection intervals were set at $\Delta t = T_{p5} / 2$, with categorical labels assigned by threshold: REMOVE (immediate) for $T_{p5} < 0.1$ yr, REMOVE (urgent) for $0.1 \leq T_{p5} < 1$, ANNUAL for $1 \leq T_{p5} < 5$, BI-ANNUAL for $5 \leq T_{p5} < 30$, and STANDARD (5-yr) for $T_{p5} \geq 30$. Under lower-bound S-N, *all* tested decay levels in Incheon and Jeju with $r_d/R \geq 0.5$ fall into REMOVE categories, compared with $r_d/R \geq 0.6$ for Seoul. Even at the modest decay level $r_d/R = 0.5$, current Korean street-tree inspection guidelines (Ministry of Land 2019) — which prescribe routine annual visual checks irrespective of decay state — leave a 21-72 % fraction of the eccentric-tree population unprotected (Table 2).

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

## 4. Discussion (outline only)

*See `paper2_outline.md §4` for the outline. Full text to be drafted in v2.*

Topics: wind-direction collapse theorem; comparison with Mattheck and Breloer (1994); inspection-cycle policy implications; backwards compatibility with Paper I; limitations (single-cavity assumption, sparse tomography priors, independence between cavity and wind orientation); bridge to Paper I future-work items 4-6.

---

## 5. Conclusions (outline only)

*Five-claim numbered list, ~ 200 words, per `paper2_outline.md §5`.*

---

## Acknowledgements

This work used computational resources of the AGE (AI · Green Intelligence · Engineering) Lab, Seoul National University. The framework is the second instalment of a fatigue-assessment series; we thank the colleagues who reviewed Paper I in 2026.

## Funding

This work was supported by the Korea Environment Industry & Technology Institute (KEITI) through the Climate Change R&D Project for the New Climate Regime, funded by the Korea Ministry of Environment (MOE) (2022003570004), and by the National Research Foundation of Korea Grant funded by the Korean Government (NRF-RS-2023-00259403).

## Declarations

**Conflict of interest** The author declares no competing interests.

**Availability of data and material** The Monte Carlo sample data and analysis scripts that support the findings of this study are openly available at `https://github.com/<author>/fatigue-tree2`. The Paper I anchor values used in Section 2.5 are reproduced here verbatim from Kang (2026, Table 2).

**Use of AI-assisted tools** The author used an AI language model (Claude, Anthropic) to assist with manuscript drafting and language editing. The author reviewed and revised all AI-generated content and takes full responsibility for the accuracy and integrity of the work.

## References (placeholder; see `paper2_outline.md §References`)
