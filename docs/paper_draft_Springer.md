<!-- Target journal: Trees – Structure and Function (Springer)
     ISSN: 0931-2205 (print) | 1432-2285 (electronic)
     Scope: physiology, biomechanics, structure, and ecology of trees and woody plants
     Max. ~8,000 words (main text) | 5 keywords | unstructured abstract ~250 words -->

# Fatigue Life Assessment of Urban Street Trees under Wind Loading: Effects of Decay Progression and Soil Conditions

**Junsuk Kang** a,b,c,d *

a Interdisciplinary Program in Landscape Architecture, Seoul National University, Seoul 08826, Republic of Korea

b Transdisciplinary Program in Smart City Global Convergence, Seoul National University, Seoul 08826, Republic of Korea

c Department of Landscape Architecture and Rural Systems Engineering, Seoul National University, Seoul 08826, Republic of Korea

d Research Institute of Agriculture and Life Sciences, Seoul National University, Seoul 08826, Republic of Korea

\* Correspondence: junkang@snu.ac.kr; Tel.: +82-2-880-4872

---

## Abstract

Urban street trees are susceptible to windthrow and branch failure during storms, yet existing risk frameworks focus on static ultimate-load criteria and neglect cumulative fatigue damage under repeated wind loading. This study presents an integrated fatigue life assessment framework coupling finite element dynamic analysis (ABAQUS), Kaimal-spectrum wind time series, rainflow cycle counting, Miner's damage rule, a concentric-hollow decay stress amplification model, and regional Weibull wind speed distributions. Six tree geometries (H = 6 and 8 m; DBH = 15, 20, and 25 cm) were modelled as tapered B31 beams with biaxial wind loading (longitudinal mean-plus-turbulence and lateral turbulence at σ_v = 0.75σ_u) and Rayleigh damping (ζ = 2%). The maximum dynamic amplification factor (DAF) was 7.71 for H = 8 m, DBH = 15 cm, with a peak stress of 138.7 MPa exceeding the ginkgo modulus of rupture (MOR), indicating that static failure precedes fatigue in slender trees. For a decayed reference tree (r_d/R = 0.8) under Jeju typhoon-regime winds (Weibull scale parameter c = 9.0 m/s), lower-bound S-N fatigue life was approximately 1.4 years. Sound-tree fatigue life differed by up to several hundred thousand times across five Korean street tree species due to MOR differences alone. A Winkler rotational spring model showed that soft-clay soils reduce basal stress by up to 64.5% relative to a fully fixed base, but at the cost of increased overturning risk under extreme loads. The results establish that decay detection and species selection are the primary determinants of fatigue safety for urban trees, and that regionally differentiated inspection standards are warranted.

**Keywords:** tree risk assessment, wind-induced fatigue, finite element dynamic analysis, decay stress amplification, Weibull wind distribution

---

## 1. Introduction

### 1.1 Background and Motivation

Urban street trees provide essential ecosystem services including microclimate regulation, air-quality improvement, and urban aesthetics (Nowak and Crane 2002). However, during typhoons and severe windstorms, tree branch failures and stem windthrow represent significant hazards to pedestrians, vehicles, and urban infrastructure. In the Republic of Korea, thousands of street tree failures are reported annually during typhoon events, with damage concentrated along coastal regions, particularly on Jeju Island (Korea Forest Service 2021).

Current tree risk assessment practice relies primarily on field inspection methods such as resistograph drilling, sonic tomography, and visual tree assessment (VTA) to characterize decay extent (Matheny and Clark 1994; Smiley and Fraedrich 1992). These methods inform static stability evaluations based on critical wind-speed thresholds for stem failure or root-plate overturning (Peltola et al. 1999; Hale et al. 2015). However, they do not explicitly account for the cumulative fatigue damage caused by everyday wind-induced cyclic loading over a tree's lifespan.

Wood, like metals, exhibits fatigue degradation under repeated loading. Field studies have documented that trees experience millions of bending stress cycles throughout their lifetimes (Niklas 1992), and laboratory experiments confirm that wood specimens fail at stresses well below the static MOR under cyclic loading (Clorius et al. 2000). Decayed trees are of particular concern because hollow cross-sections amplify bending stresses nonlinearly, yet quantitative models linking decay progression to fatigue life under realistic wind loading have not been developed.

Root systems are equally vulnerable to wind-induced cyclic loading. Experimental studies on isolated tree roots subjected to repeated tensile loads have demonstrated that roots undergo permanent plastic deformation beyond their elastic limit, with cumulative damage accumulating progressively with each load cycle (Cofie and Koolen 2001; Liu 1994). Even under relatively few loading cycles, root failure can occur: beech (*Fagus sylvatica*) roots have been observed to rupture after approximately 1,015 cycles, while larch roots failed after as few as 68 cycles under repeated tension (Cofie and Koolen 2001). Mechanical disturbance from both wind events and forestry vehicles induces these cyclic stresses, suggesting that root fatigue may govern whole-tree stability in weakened or repeatedly disturbed urban trees.

### 1.2 Limitations of Existing Research

Sellier and Fourcaud (2009) measured in-situ damping ratios and developed simplified dynamic tree models but did not extend their analysis to long-term fatigue. Niklas (1992) established empirical S-N relationships for wood but in the context of standard specimens rather than whole trees with decay. Mattheck and Breloer (1994) proposed the t/R = 0.3 mechanical failure criterion for decayed stems based on static analysis, while its implications for fatigue life remain unquantified.

A notable counterpoint is offered by Taylor (2023), who proposed a theoretical model suggesting that annual radial growth continuously replenishes damaged wood tissue, thereby largely mitigating long-term fatigue accumulation in tree stems. While this self-repair mechanism is relevant for healthy, vigorously growing trees, it does not address three scenarios of primary concern in urban tree management: (i) decayed trees in which fungal degradation eliminates functional xylem faster than radial growth can compensate; (ii) root systems whose geometry does not effectively redistribute repeated tensile loading; and (iii) high-wind-speed environments such as Korean typhoon-affected coastal regions, where annual cumulative damage rates may exceed the biological repair capacity of even structurally sound trees. The quantitative results of the present study reinforce this distinction: for a sound ginkgo in Seoul, the annual fatigue damage index is of order 10$^{-3}$ yr$^{-1}$ (fatigue life ≈ 1,943 years), consistent with Taylor's (2023) argument that growth repair is adequate for routine winds. However, once the decay ratio reaches $r_d/R = 0.7$, the annual damage index rises to approximately 0.14 yr$^{-1}$ (fatigue life ≈ 7 years), a rate that far exceeds the structural restoration achievable through annual radial increment. Under Jeju wind conditions the corresponding damage rate at $r_d/R = 0.7$ becomes critical within one growing season. The present study therefore focuses specifically on these limiting — and practically most hazardous — conditions where biological repair cannot compensate for mechanical fatigue accumulation.

At the material level, a comprehensive review by Yang et al. (2025) of wood fatigue across diverse loading conditions provides a critical mechanistic context. Their synthesis of S-N data from multiple independent studies demonstrates that reversed bending loading (stress ratio R = −1), which closely approximates the biaxial oscillation of a wind-loaded tree stem, is substantially more damaging than the repeated loading (R = 0) that underlies most laboratory wood fatigue experiments including those of Niklas (1992) and Clorius et al. (2000). Furthermore, Yang et al. (2025) show that residual stiffness in wood does not decrease monotonically with fatigue cycles — unlike in metals — because the trajectory of stiffness change depends on load direction and grain orientation. These findings suggest that: (i) the S-N relationships adopted in this study represent upper-bound (non-conservative) estimates of fatigue life relative to the actual in-situ loading regime; and (ii) Miner's linear damage rule, while practical, does not fully capture wood's complex fatigue response. Both observations reinforce rather than undermine the present approach — they justify the use of lower-bound S-N parameters for all safety-critical assessments and motivate the conservative design philosophy adopted throughout.

National design standards and tree management guidelines in Korea (Ministry of Land, Infrastructure and Transport, 2019; Korea Forest Service, 2020) prescribe inspection intervals and pruning schedules without a quantitative residual-fatigue-life framework.

### 1.3 Research Objectives

The objectives of this study are to:

1. Quantify the dynamic stress histories of trees across a range of stem geometries using ABAQUS finite element dynamic analysis with biaxial Kaimal-spectrum wind loading and Rayleigh damping;
2. Calculate annual fatigue damage rates and expected fatigue lives by integrating Rainflow-counted stress cycles with Miner's rule and regional Weibull wind distributions;
3. Quantify the nonlinear effect of decay progression on fatigue life using a concentric-hollow stress amplification model and identify the mechanical significance of the Mattheck critical decay ratio;
4. Compare fatigue lives across five major Korean street tree species based on their MOR values; and
5. Evaluate the influence of soil conditions on basal stress and fatigue life using a Winkler rotational spring model, and characterize the trade-off between fatigue and overturning safety.

The overall study workflow is summarized in Figure 1.

---

## 2. Materials and methods

### 2.1 Structural Tree Model

Trees were modeled as tapered cantilever beams using ABAQUS B31 two-node Timoshenko beam elements. A linear taper was assumed such that the cross-sectional diameter varies with height as:

$$D(z) = D_{\text{base}} - (D_{\text{base}} - D_{\text{tip}}) \cdot \frac{z}{H} \tag{1}$$

where $D_{\text{base}} = 1.2 \times \text{DBH}$ is the base diameter, $D_{\text{tip}} = 0.3 \times \text{DBH}$ is the tip diameter, and $H$ is tree height. These taper ratios are consistent with allometric measurements reported by Niklas and Spatz (2004). The model geometry is illustrated in Figure 2a. Six tree geometries combining stem heights $H \in \{6, 8\}$ m and $\text{DBH} \in \{15, 20, 25\}$ cm were analyzed (Table 1).

Material properties for ginkgo (*Ginkgo biloba*), the most common street tree in Seoul, were adopted:
- Elastic modulus: $E = 9.5$ GPa (fiber direction; KS F 2207, 2019)
- Poisson's ratio: $\nu = 0.30$
- Density: $\rho = 600$ kg/m³ (green wood)

The base of each tree was fully fixed (encastre boundary condition) for the primary analysis. The effect of soil compliance was evaluated separately using the Winkler model described in Section 2.7.

### 2.2 Wind Velocity Time Series

Wind loading time histories were generated using the Kaimal et al. (1972) power spectral density (PSD) and the Shinozuka and Deodatis (1991) inverse fast Fourier transform (IFFT) spectral representation method. The longitudinal (along-wind) turbulence PSD is:

$$S_u(f) = \frac{4\,\sigma_u^2\,L_u/\bar{U}}{\left(1 + 6f\,L_u/\bar{U}\right)^{5/3}} \tag{2}$$

where $\sigma_u$ is the longitudinal turbulence standard deviation, $L_u = 100$ m is the integral length scale, $\bar{U} = 20$ m/s is the reference mean wind speed, and $f$ is frequency (Figure 3a). The lateral (cross-wind) turbulence standard deviation was set to:

$$\sigma_v = 0.75\,\sigma_u \tag{3}$$

consistent with ESDU (1985) atmospheric boundary layer measurements. The mean lateral wind speed was set to zero (pure lateral turbulence). Regional Weibull wind speed distributions are shown in Figure 3b. Wind pressure forces were distributed as nodal loads over ten height intervals along the stem:

$$F(z,t) = \frac{1}{2}\,\rho_{\text{air}}\,C_D\,A(z)\,\left[\bar{U} + u'(t)\right]^2 \tag{4}$$

where $\rho_{\text{air}} = 1.225$ kg/m³, $C_D = 1.2$ for a circular cross-section, and $A(z)$ is the projected area at height $z$. A simulation duration of 600 s (10 min) with a time step $\Delta t = 0.05$ s (20 Hz) was used.

### 2.3 Finite Element Dynamic Analysis

ABAQUS Standard's implicit dynamic procedure (`*Dynamic, Implicit`) was used for geometrically nonlinear (NLGEOM=ON) time-history analysis in two steps:

- **Step 1 (FREQUENCY):** Lanczos eigensolver to extract the fundamental natural frequency $f_1$.
- **Step 2 (DYNAMIC):** Nonlinear time integration under the 600 s biaxial wind load history.

Rayleigh damping was specified as:

$$[C] = \alpha[M] + \beta[K] \tag{5}$$

A damping ratio $\zeta = 2\%$ was selected, consistent with the mid-range of field measurements by Sellier and Fourcaud (2009). Rayleigh coefficients were computed as:

$$\alpha = \frac{2\zeta\,\omega_1\omega_2}{\omega_1+\omega_2}, \quad \beta = \frac{2\zeta}{\omega_1+\omega_2} \tag{6}$$

yielding $\alpha = 0.0754$ s$^{-1}$ and $\beta = 0.00472$ s. The resultant bending moment at the base was calculated from the ABAQUS section moment outputs SM1 and SM2:

$$M_{\text{res}}(t) = \sqrt{SM1(t)^2 + SM2(t)^2} \tag{7}$$

The maximum bending stress at the outermost fiber was then:

$$\sigma(t) = \frac{M_{\text{res}}(t) \cdot (D_{\text{base}}/2)}{I_{\text{base}}} \tag{8}$$

### 2.4 Rainflow Cycle Counting and Miner's Rule

Rainflow cycle counting following ASTM E1049 (2017) was applied to the 600 s stress time history at the stem base to extract stress amplitude–cycle pairs $(\sigma_{a,i},\, n_i)$. A fully reversed loading ratio $R = 0$ was assumed, so no mean-stress correction was applied.

The wood S-N relationship was expressed in normalized form based on fatigue data compiled by Niklas (1992) and Clorius et al. (2000):

$$\frac{\sigma_a}{\text{MOR}} = a - b \cdot \log_{10}(N_f) \tag{9}$$

Mean-curve coefficients $a = 0.87$, $b = 0.10$ were adopted. A conservative lower-bound curve representing the 95% confidence lower bound was defined by $a_{lb} = 0.80$, $b_{lb} = 0.125$ (Figure 5). Green-wood MOR was taken as 60% of dry-wood MOR (KS F 2207 2019). Fatigue life $N_{f,i}$ for each stress amplitude was obtained by inverting Equation (9). Cumulative fatigue damage follows Miner's linear rule:

$$D_{\text{total}} = \sum_i \frac{n_i}{N_{f,i}} \tag{10}$$

Fatigue failure is defined at $D_{\text{total}} = 1.0$.

### 2.5 Decay Stress Amplification Model

Wood decay was idealized as a concentric circular hollow at the stem base. The bending stress amplification factor relative to a sound cross-section is:

$$sf\!\left(\frac{r_d}{R}\right) = \frac{1}{1 - \left(\dfrac{r_d}{R}\right)^4} \tag{11}$$

where $r_d$ is the decay cavity radius and $R$ is the total stem radius. This expression follows directly from the section modulus ratio of a hollow circular cross-section. Mattheck and Breloer (1994) identified $t/R = 0.3$ (i.e., $r_d/R = 0.7$) as the critical mechanical decay ratio beyond which residual wall thickness becomes inadequate; this threshold was used to assess nonlinear transitions in fatigue life. The three decay states modeled are illustrated in Figure 2b.

### 2.6 Regional Weibull Wind Distribution and Annual Fatigue Damage

Annual fatigue damage was computed by integrating wind-speed-weighted damage increments using the Weibull probability density function:

$$f(U) = \frac{k}{c}\left(\frac{U}{c}\right)^{k-1} \exp\!\left[-\left(\frac{U}{c}\right)^k\right] \tag{12}$$

Weibull shape $k$ and scale $c$ parameters for three Korean sites were derived from Korea Meteorological Administration wind atlas data (KMA 2020) and are listed in Table 3. Stress histories were scaled from the reference simulation ($\bar{U} = 20$ m/s) to each wind speed bin $U_j$ using $\sigma \propto U^2$. Annual damage was computed as:

$$D_{\text{annual}} = \sum_j D(U_j) \cdot f(U_j) \cdot \Delta U \cdot 8760 \tag{13}$$

where $\Delta U = 0.5$ m/s. Expected fatigue life is $T_f = 1/D_{\text{annual}}$ (years).

### 2.7 Winkler Rotational Spring Model

Root–soil interaction was represented by a Winkler rotational spring at the stem base. The equivalent rotational stiffness is:

$$K_r = k_s \cdot D_{\text{base}} \cdot \frac{L_r^3}{3} \tag{14}$$

where $k_s$ is the coefficient of subgrade reaction and $L_r = 1.0$ m is the effective root anchorage depth (mid-range for urban street trees; Gilman, 2012). The basal stress correction factor accounting for base rotation is:

$$C_\sigma = \frac{1}{1 + 3EI / (K_r \cdot H)} \tag{15}$$

$C_\sigma = 1$ corresponds to the fully fixed case; softer soils yield $C_\sigma < 1$, reducing basal bending stress. Values of $k_s$ for each soil category are provided in Table 8.

---

## 3. Results

### 3.1 Dynamic Analysis: Natural Frequencies and Stress Amplification

Dynamic analysis results for all six cases are summarized in Table 4. The fundamental natural frequency $f_1$ decreases with increasing tree height ($H \uparrow \Rightarrow f_1 \downarrow$) and increases with increasing stem diameter ($\text{DBH} \uparrow \Rightarrow f_1 \uparrow$). The H = 8 m, DBH = 15 cm case yields the lowest $f_1 = 0.279$ Hz, which falls within the energy-rich range of the Kaimal PSD (0.1–0.5 Hz), producing the highest DAF of 7.71.

The peak bending stress for this case (138.7 MPa) substantially exceeds the estimated green-wood MOR of ginkgo (~25.4 MPa), indicating that static failure would precede fatigue accumulation at the reference wind speed of 20 m/s. Accordingly, the H = 8 m, DBH = 20 cm case was selected as the **reference case** for fatigue life parametric studies. For the H = 6 m, DBH = 25 cm case, the DAF falls below unity (0.485), reflecting inertial filtering because $f_1$ exceeds the dominant energy band of the wind spectrum (Figure 4).

### 3.2 Fatigue Life as a Function of Decay Ratio

Fatigue lives for the reference case (H = 8 m, DBH = 20 cm; ginkgo; Seoul wind) across decay ratios from 0 to 0.8 are presented in Table 5. A sound tree ($r_d/R = 0$) has a mean S-N fatigue life of approximately 1,943 years, indicating that wind fatigue alone does not govern the lifespan of healthy urban trees under Seoul wind conditions. Fatigue life decreases monotonically with increasing decay ratio, but the reduction is strongly nonlinear.

At the Mattheck critical ratio $r_d/R = 0.7$, the mean S-N fatigue life drops to 6.9 years and the lower-bound S-N fatigue life to 1.7 years—a reduction of approximately 280-fold relative to the sound-tree mean estimate. From $r_d/R = 0.7$ to $0.8$, the mean fatigue life further collapses from 6.9 to 0.4 years (~17-fold reduction). This accelerating decline is a direct consequence of the fourth-power amplification in Equation (11). The nonlinear transition at $r_d/R \approx 0.7$ constitutes a critical threshold for tree risk management (Figure 6).

### 3.3 Regional Variation in Fatigue Life

Table 6 compares fatigue lives for a sound ginkgo (H = 8 m, DBH = 20 cm) under the three regional wind regimes. Fatigue life decreases from 1,943 years in Seoul (Weibull $c = 4.0$ m/s) to less than one year in Jeju ($c = 9.0$ m/s) under mean S-N conditions—a reduction exceeding 2,000-fold. This extreme sensitivity reflects the combined effect of wind speed entering as a squared term in the aerodynamic force and of stress amplitude entering the S-N relationship exponentially. The Weibull scale parameter $c$ is thus the dominant controlling variable for regional fatigue life (Figure 8).

### 3.4 Species Comparison

Table 7 presents fatigue lives for five major Korean street tree species with MOR values ranging from 34.8 MPa (dawn redwood, *Metasequoia glyptostroboides*) to 63.0 MPa (zelkova, *Zelkova serrata*). The structural geometry (H = 8 m, DBH = 20 cm) and wind loading are held constant, so differences arise solely from species-specific MOR. Zelkova exhibits a fatigue life exceeding $10^7$ years for the sound case, while dawn redwood yields only 32 years—a difference of more than five orders of magnitude. At $r_d/R = 0.7$, dawn redwood fatigue life falls to approximately 1 year, whereas zelkova retains ~780,000 years. The result underscores the dominant role of wood strength in determining long-term fatigue safety (Figure 7).

### 3.5 Influence of Soil Conditions

Table 8 summarizes the Winkler spring model results. For soft clay ($k_s = 1{,}500$ kN/m³), the basal stress correction factor is $C_\sigma = 0.355$, reducing base bending stress by 64.5% relative to the fully fixed case and extending fatigue life well beyond manageable limits. However, this stress reduction is accompanied by a proportional increase in base rotation under extreme wind loading, raising the risk of whole-tree overturning and root-plate failure. For dense sand ($k_s = 40{,}000$ kN/m³), $C_\sigma = 0.923$ and fatigue life is only marginally longer than the fully fixed case. These results reveal a fundamental **trade-off**: compliant soils improve fatigue performance but compromise overturning stability, and both failure modes must be evaluated simultaneously in tree risk management (Figure 9).

---

## 4. Discussion

### 4.1 Governing Fatigue Mode: Decay-Accelerated Low-Cycle Fatigue

The central finding of this study is that the realistic failure scenario for urban street trees is not governed by high-cycle fatigue (HCF) of sound wood under ordinary winds, but rather by **decay-accelerated low-cycle fatigue (LCF, here defined as fatigue failure accumulating within approximately $10^3$ cycles)**—the combination of decay-induced stress amplification and infrequent high-wind events. This distinction has important practical implications. A sound ginkgo on a Seoul street would require thousands of years to accumulate unit fatigue damage, meaning that fatigue per se is not a design-limiting state for well-maintained trees. However, once decay progresses beyond $r_d/R \approx 0.5$–$0.7$, the stress amplification factor rises steeply and relatively ordinary wind speeds can rapidly accumulate critical damage. This transition is abrupt because fatigue damage accumulation scales with the ratio $(\sigma_a/\text{MOR})^{1/b}$, and the exponent $1/b \approx 10$ means that even modest stress increases produce order-of-magnitude changes in fatigue life.

The implication is that decay detection, not structural design wind loads, is the primary lever for managing fatigue risk in urban trees. Periodic non-destructive inspection (resistograph, sonic tomography) aimed at tracking the progression of $r_d/R$ through the critical threshold is therefore more consequential than refinements to wind load models (Figure 10).

As an order-of-magnitude check on the reasonableness of these results, three independent comparisons are noted. First, the dynamic amplification factors computed here (DAF = 1.1–3.9 for the H = 6 m cases, up to 7.7 for H = 8 m) are consistent with field measurements reported by James et al. (2006), who observed DAF values of 1.5–4.0 for open-grown trees of similar dimensions under turbulent wind loading. Second, the lower-bound fatigue life for a sound reference tree (~60 years) is broadly consistent with the observed service life of well-maintained urban street trees in Korean municipalities, typically 50–100 years before structural concerns prompt removal (Korea Forest Service 2020). Third, the sharp transition in fatigue life at $r_d/R \approx 0.7$ corroborates the Mattheck critical threshold ($t/R = 0.3$, i.e., $r_d/R = 0.7$) derived from an entirely independent static-mechanics argument (Mattheck and Breloer 1994), lending cross-methodological support to this criterion.

### 4.2 Uncertainty in Wood S-N Properties

The 30-fold gap between mean and lower-bound S-N fatigue lives (e.g., 1,943 versus 60 years for a sound tree) reflects the large variability in wood fatigue properties (coefficient of variation approximately 15–25%; Clorius et al., 2000). For risk management purposes, the lower-bound S-N curve is the appropriate choice, as it corresponds to the conservative end of the population. Under lower-bound conditions, fatigue life falls below 30 years even at a decay ratio of only 0.5, suggesting that current Korean tree safety inspection standards—which prescribe simple visual checks and do not trigger mandatory removal until advanced decay is visible—may be insufficient for high-wind-exposure sites. A risk-informed approach linking inspection interval to decay ratio and site wind climate is recommended.

An additional source of non-conservatism, noted by Yang et al. (2025), is the stress ratio (R-ratio) discrepancy between laboratory S-N tests and field loading. The S-N curves used in this study were derived from repeated-loading tests (R = 0, i.e., zero-to-maximum stress cycling), whereas wind-induced stem oscillation subjects trees to reversed bending (R ≈ −1, i.e., tension–compression cycling) as the wind direction fluctuates. Yang et al. (2025) document substantially shorter fatigue lives for wood under reversed loading compared to repeated loading across multiple species. Consequently, the fatigue lives reported in this study should be interpreted as upper-bound estimates; the true in-situ fatigue risk may be more severe. This further strengthens the recommendation to use lower-bound S-N parameters and to adopt a precautionary stance toward decay management thresholds.

### 4.3 Regional Differentiation of Management Standards

The 2,000-fold fatigue life difference between Seoul and Jeju demonstrates that a single national standard for tree inspection intervals and decay-ratio removal thresholds may be inadequate. Under lower-bound S-N conditions, a decay ratio of 0.5 that yields approximately 25 years of fatigue life in Seoul may correspond to only a few years in Jeju, given the much higher Weibull scale parameter. This analysis supports the development of regionally differentiated tree risk management standards, analogous to the site-specific wind zones already recognized in Korean structural design codes (KBC 2021). Coastal and island locations should adopt more conservative decay-ratio thresholds (e.g., $r_d/R \leq 0.3$ for mandatory intervention) and shorter inspection cycles than inland urban areas.

The trade-off between fatigue safety and overturning risk identified in Section 3.5 warrants a more precise engineering criterion. In the Winkler rotational spring model, the overturning limit state is reached when the base rotation $\theta = M_{\text{res}} / K_r$ exceeds the critical root-plate rotation, commonly taken as $\theta_{\text{cr}} \approx 5$–$10°$ based on field pull-test records for urban trees (Gilman 2012). Combining this with the fatigue limit state defines a two-constraint design space: a tree is safe only if both the annual fatigue damage index $D_{\text{annual}} < 1/T_{\text{target}}$ (where $T_{\text{target}}$ is the target service life) and the extreme-event base rotation $\theta < \theta_{\text{cr}}$. The approximate boundary separating Domain II (decay-accelerated fatigue) from Domain III (static overturning) in Figure 10 corresponds to wind conditions where the peak base rotation from a 50-year return-period wind event equals $\theta_{\text{cr}}$ for a fully fixed base, providing a practical demarcation for site-specific risk classification.

### 4.4 Species Selection as a Fatigue Safety Criterion

The five-to-six orders of magnitude spread in fatigue life across species with MOR values from 34.8 to 63.0 MPa demonstrates that wood strength should be an explicit criterion in street tree species selection. Current Korean street tree selection guidelines prioritize aesthetics, ease of management, drought tolerance, and pollution resistance (Ministry of Land, Infrastructure and Transport 2019) but do not include structural strength metrics. The present results indicate that selecting species with higher MOR—particularly in high-wind-exposure sites—can provide a passive, maintenance-free improvement in fatigue safety far exceeding what can be achieved through management interventions alone.

It should be noted that the species comparison in this study isolates the effect of MOR by holding all other parameters constant: stem geometry (H = 8 m, DBH = 20 cm), elastic modulus ($E$ = 9.5 GPa), drag coefficient ($C_D$ = 1.2), and crown form are identical across all five species. In reality, each species has distinct allometric relationships, crown morphology, $C_D$, and $E$ values that would alter both the dynamic response and the stress distribution. Consequently, Table 7 should be interpreted as a comparison of the **material-level fatigue potential** of each species under identical structural and loading conditions, not as a prediction of actual service life for individual trees in the field. Future work should incorporate species-specific allometric scaling, crown drag models, and green-wood elastic moduli to produce more realistic species-level fatigue life predictions.

### 4.5 Model Limitations and Future Research

The main limitations of this study are as follows. First, the single-stem beam model excludes the mass and aerodynamic drag of branches and foliage. Crown damping has been shown to reduce dynamic amplification factors by 20–50% relative to a bare-stem model (James et al. 2006; Sellier and Fourcaud 2009), meaning that the present study overestimates stem stress and yields conservatively short fatigue life predictions. While this conservative bias is appropriate for a safety-oriented framework, future implementations should incorporate multi-body crown–stem coupling models to improve accuracy for specific site assessments. Second, decay is treated as a static parameter rather than a time-evolving process; coupling a decay growth model (e.g., Schwarze et al., 2000) with the fatigue accumulation framework is a priority for future work. Third, the S-N curves adopted here are based on repeated loading ($R = 0$), whereas in-situ wind loading induces reversed bending ($R \approx -1$) as described in Section 4.2. A Goodman mean-stress correction — which adjusts the allowable stress amplitude as $\sigma_a = \sigma_{a,R=0}(1 - \sigma_m / \text{MOR})$ — could be applied to bridge this gap. For a tree with modest self-weight lean, $\sigma_m$ may be 5–15% of MOR, yielding a Goodman correction factor of 0.85–0.95 and a corresponding reduction in predicted fatigue life of 10–30%. While a full $R$-ratio sensitivity analysis is beyond the scope of this study, even this simplified estimate confirms that the reported fatigue lives are upper-bound values and reinforces the necessity of lower-bound S-N parameters for conservative practice. Fourth, no field validation against measured strain histories or documented wind failure records was performed. Instrumented monitoring campaigns on trees of known decay status would provide essential validation data.

Fifth, the concentric hollow-cylinder decay model assumes structural symmetry and thereby maximises the section modulus for a given $r_d/R$. In practice, decay pockets are frequently eccentric: Mattheck and Breloer (1994) show that asymmetric defects can increase peak bending stress by 20–40% relative to the concentric case at the same decay ratio, producing proportionally shorter fatigue lives. The use of lower-bound S-N parameters in this study partially compensates for this non-conservatism; nevertheless, future work should incorporate probabilistic decay-geometry models to quantify the associated uncertainty.

Sixth, and most fundamentally, the S-N relationships applied here were derived from repeated-loading (R = 0) specimens, whereas wind-induced stem oscillation produces reversed bending (R ≈ −1), the loading regime Yang et al. (2025) identify as the most damaging for wood. Developing species-specific S-N curves under reversed loading conditions is a necessary step toward more accurate fatigue life prediction for trees. Seventh, the present study treats wood as a homogeneous elastic material, neglecting viscoelastic creep and progressive microstructural damage documented by Yang et al. (2025): fibre buckling in bending, kink-band formation in compression, and microcracking in tension. In living trees, tensile microcracks may penetrate pit membranes in xylem vessels, disrupting hydraulic conductance — a coupled mechanical–hydraulic fatigue mechanism entirely absent from conventional engineering fatigue models. Investigating this pathway, in which cyclic mechanical damage triggers embolism formation and cascading hydraulic failure, represents a promising frontier in tree biomechanics that could ultimately reconcile the mechanical and physiological dimensions of tree wind risk.

---

## 5. Conclusions

This study developed and applied an integrated framework for quantitative fatigue life assessment of urban street trees under wind loading. The main conclusions are:

1. Dynamic amplification factors up to **7.71** were obtained for the H = 8 m, DBH = 15 cm tree, with the largest amplification occurring when the fundamental frequency coincided with the Kaimal PSD energy peak (0.1–0.5 Hz). For slender trees, static failure under extreme wind events precedes fatigue accumulation.

2. A nonlinear fatigue life threshold exists at the Mattheck critical decay ratio $r_d/R = 0.7$. Under Seoul wind conditions and lower-bound S-N, fatigue life at this decay ratio is approximately **1.7 years**, necessitating immediate intervention following detection.

3. Jeju's wind climate ($c = 9.0$ m/s) reduces fatigue life by more than **2,000-fold** relative to Seoul ($c = 4.0$ m/s), demonstrating that regionally differentiated inspection and removal standards are essential.

4. Species MOR differences produce fatigue life variations of up to **five to six orders of magnitude** across common Korean street trees, establishing wood strength as a critical species selection criterion.

5. Soft clay soils reduce basal bending stress by up to 64.5% and extend fatigue life, but at the cost of increased overturning risk—a trade-off that requires simultaneous evaluation of both fatigue and ultimate-limit-state failure modes.

6. The governing failure scenario for urban trees is not conventional high-cycle fatigue of sound wood, but **decay-accelerated low-cycle fatigue**. Proactive decay detection combined with selection of high-MOR species constitutes the most effective strategy for ensuring the long-term fatigue safety of urban street trees.

---

## Acknowledgements

ABAQUS finite element software was used under a Seoul National University academic license.

---


## Funding

This work was supported by knowledge-based environmental service program and Korea Environment Industry & Technology Institute (KEITI) through Climate Change R&D Project for New Climate Regime, funded by the Korea Ministry of Environment (MOE) (2022003570004), and the National Research Foundation of Korea Grant funded by the Korean Government (NRF-RS-2023-00259403).

## Declarations

**Conflict of interest** The author declares that he has no conflict of interest.

**Availability of data and material** The Python analysis scripts and processed stress history data that support the findings of this study are available from the corresponding author upon reasonable request. ABAQUS input files (.inp) are available upon request subject to software licensing constraints.

**Author contributions** JK: Conceptualization, Methodology, Software, Formal analysis, Investigation, Writing – original draft, Writing – review and editing, Visualization, Funding acquisition.

**Ethics approval** Not applicable.

**Use of AI-assisted tools** The author used an AI language model (Claude, Anthropic) to assist with manuscript drafting and language editing. The author reviewed and revised all AI-generated content and takes full responsibility for the accuracy and integrity of the work.

## References

ASTM E1049-85 (2017) *Standard Practices for Cycle Counting in Fatigue Analysis*. ASTM International, West Conshohocken, PA, USA.

Clorius, C.O., Pedersen, M.U., Hoffmeyer, P., and Damkilde, L. (2000) Compressive fatigue in wood. *Wood Sci Technol* 34(1):21–37. https://doi.org/10.1007/s002260050005

Cofie, P., and Koolen, A.J. (2001) Test speed and other factors affecting the measurements of soil shear strength parameters and angle of internal friction of tree roots. *Soil Tillage Res* 63(1–2):51–56. https://doi.org/10.1016/S0167-1987(01)00216-1

ESDU (1985) *Characteristics of Atmospheric Turbulence Near the Ground. Part II: Single Point Data for Strong Winds (Neutral Atmosphere)*. Engineering Sciences Data Unit Item 85020. ESDU International, London, UK.

Gilman, E.F. (2012) *An Illustrated Guide to Pruning*, 3rd ed. Cengage Learning, Clifton Park, NY, USA.

Hale, S.E., Gardiner, B., Peace, A., Nicoll, B., Taylor, P., and Pizzirani, S. (2015) Comparison and validation of three versions of a forest wind risk model. *Environ Model Softw* 68:27–41. https://doi.org/10.1016/j.envsoft.2015.01.016

James, K.R., Haritos, N., and Ades, P.K. (2006) Mechanical stability of trees under dynamic loads. *Am J Bot* 93(10):1522–1530. https://doi.org/10.3732/ajb.93.10.1522

Kaimal, J.C., Wyngaard, J.C., Izumi, Y., and Coté, O.R. (1972) Spectral characteristics of surface-layer turbulence. *Q J R Meteorol Soc* 98(417):563–589. https://doi.org/10.1002/qj.49709841707

KBC (2021) *Korean Building Code and Commentary*. Ministry of Land, Infrastructure and Transport, Sejong, Republic of Korea.

KMA (Korea Meteorological Administration) (2020) *Korean Wind Atlas 2.0 Technical Report*. KMA, Seoul, Republic of Korea.

Korea Forest Service (2020) *Regulations for Street Tree Establishment and Management*. Korea Forest Service Notice No. 2020-25. Daejeon, Republic of Korea.

Korea Forest Service (2021) *Survey of Urban Street Tree Wind Damage and Management Strategies*. Korea Forest Service Research Report No. 21-25. National Institute of Forest Science, Seoul, Republic of Korea.

KS F 2207 (2019) *Test Methods for Bending Strength of Wood*. Korean Industrial Standard. Korean Agency for Technology and Standards, Eumseong, Republic of Korea.

Liu, X. (1994) *Strength and Failure of Wood Subjected to Combined Stresses*. PhD Thesis, University of Wisconsin-Madison, Madison, WI, USA.

Matheny, N., and Clark, J. (1994) *A Photographic Guide to the Evaluation of Hazard Trees in Urban Areas*, 2nd ed. International Society of Arboriculture, Champaign, IL, USA.

Mattheck, C., and Breloer, H. (1994) *The Body Language of Trees: A Handbook for Failure Analysis*. HMSO, London, UK.

Miner MA (1945) Cumulative damage in fatigue. J Appl Mech 12:A159–A164. https://doi.org/10.1115/1.4009458

Ministry of Land, Infrastructure and Transport (2019) *Commentary on the Rules for the Structure and Facilities of Roads*. Ministry of Land, Infrastructure and Transport, Sejong, Republic of Korea.

Niklas, K.J. (1992) *Plant Biomechanics: An Engineering Approach to Plant Form and Function*. University of Chicago Press, Chicago, IL, USA.

Niklas, K.J., and Spatz, H.C. (2004) Growth and hydraulic (not mechanical) constraints govern the scaling of tree height and mass. *Proc Natl Acad Sci USA* 101(44):15661–15663. https://doi.org/10.1073/pnas.0405857101

Nowak, D.J., and Crane, D.E. (2002) Carbon storage and sequestration by urban trees in the USA. *Environ Pollut* 116(3):381–389. https://doi.org/10.1016/S0269-7491(01)00214-7

Peltola, H., Kellomäki, S., Väisänen, H., and Ikonen, V.P. (1999) A mechanistic model for assessing the risk of wind and snow damage to single trees and stands of Scots pine, Norway spruce, and birch. *Can J For Res* 29(6):647–661. https://doi.org/10.1139/x99-029


Schwarze, F.W.M.R., Engels, J., and Mattheck, C. (2000) *Fungal Strategies of Wood Decay in Trees*. Springer, Berlin, Germany.

Sellier, D., and Fourcaud, T. (2009) Crown structure and wood properties: Influence on tree sway and response to high winds. *Am J Bot* 96(5):885–896. https://doi.org/10.3732/ajb.0800226

Shinozuka, M., and Deodatis, G. (1991) Simulation of stochastic processes by spectral representation. *Appl Mech Rev* 44(4):191–204. https://doi.org/10.1115/1.3119501
Smiley, E.T., and Fraedrich, B.R. (1992) Determining strength loss from decay. *J Arboric* 18(4):201–204.

Taylor, D. (2023) Mechanical fatigue in trees mitigated by annual growth: A theoretical model. *J Theor Biol* 556:111301. https://doi.org/10.1016/j.jtbi.2022.111301

Yang, C., Abdelrahman, M., Khaloian-Sarnaghi, A., and van de Kuilen, J.W. (2025) A comprehensive analysis of fatigue in wood and wood products. *Int J Fatigue* 194:108807. https://doi.org/10.1016/j.ijfatigue.2025.108807

---

## List of Tables

**Table 1.** Tree geometry parameters for the six analyzed cases.

**Table 2.** Stress amplification factor as a function of decay ratio.

**Table 3.** Weibull wind speed distribution parameters for three Korean sites.

**Table 4.** ABAQUS dynamic analysis results for the six tree cases.

**Table 5.** Fatigue life as a function of decay ratio (reference case: H = 8 m, DBH = 20 cm; ginkgo; Seoul wind regime).

**Table 6.** Regional comparison of sound-tree fatigue life (reference case: H = 8 m, DBH = 20 cm; ginkgo; mean S-N).

**Table 7.** Species comparison of fatigue life (reference geometry: H = 8 m, DBH = 20 cm; Seoul wind; mean S-N).

**Table 8.** Effect of soil condition on basal stress correction factor and fatigue life (reference case: H = 8 m, DBH = 20 cm; ginkgo; Seoul wind; mean S-N).

---

## List of Figures

**Figure 1.** Study flowchart of the integrated wind-induced fatigue life assessment framework.

**Figure 2.** ABAQUS finite element model: (a) tapered cantilever beam with biaxial wind loading and dimensions; (b) concentric-hollow decay cross-section model for three decay states; (c) representative bending stress time history and rainflow amplitude histogram.

**Figure 3.** Wind loading characterization: (a) Kaimal turbulence PSD; (b) Weibull wind speed PDFs for Seoul, Incheon, and Jeju.

**Figure 4.** ABAQUS dynamic analysis results: (a) fundamental natural frequencies; (b) maximum bending stress and DAF.

**Figure 5.** Normalized wood S-N curves (mean and lower-bound) with experimental data.

**Figure 6.** Fatigue life versus decay ratio for Seoul and Jeju wind regimes (mean and lower-bound S-N).

**Figure 7.** Species comparison of fatigue life for sound and decayed ($r_d/R$ = 0.7) trees.

**Figure 8.** Regional fatigue life: (a) site comparison bar chart; (b) decay ratio–wind site heat map.

**Figure 9.** Soil condition effects: (a) $C_\sigma$ versus $k_s$; (b) fatigue life by soil category.

**Figure 10.** Conceptual failure mode diagram: three domains as functions of wind intensity and decay ratio.

---

## Tables

---

**Table 1.** Tree geometry parameters for the six analyzed cases.

| Case | H (m) | DBH (cm) | D_base (cm) | D_tip (cm) | D_base / D_tip ratio |
|:-----|------:|----------:|------------:|-----------:|:--------------------:|
| H6m-D15cm | 6 | 15 | 18.0 | 4.5 | 4.0 |
| H6m-D20cm | 6 | 20 | 24.0 | 6.0 | 4.0 |
| H6m-D25cm | 6 | 25 | 30.0 | 7.5 | 4.0 |
| H8m-D15cm | 8 | 15 | 18.0 | 4.5 | 4.0 |
| H8m-D20cm | 8 | 20 | 24.0 | 6.0 | 4.0 |
| H8m-D25cm | 8 | 25 | 30.0 | 7.5 | 4.0 |

H = stem height; DBH = diameter at breast height; D_base = base diameter (1.2 × DBH); D_tip = tip diameter (0.3 × DBH).

---

**Table 2.** Stress amplification factor as a function of decay ratio.

| Decay ratio $r_d/R$ | 0.00 | 0.30 | 0.50 | 0.60 | 0.70 | 0.75 | 0.80 |
|:-------------------:|:----:|:----:|:----:|:----:|:----:|:----:|:----:|
| $sf$ | 1.00 | 1.01 | 1.07 | 1.15 | 1.41 | 1.68 | 2.29 |

$sf = 1 / [1 - (r_d/R)^4]$; $r_d$ = decay cavity radius; $R$ = total stem radius.

---

**Table 3.** Weibull wind speed distribution parameters for three Korean sites.

| Site | Shape parameter $k$ | Scale parameter $c$ (m/s) | Mean wind speed $\bar{U}$ (m/s) |
|:-----|:-------------------:|:-------------------------:|:--------------------------------:|
| Seoul | 1.6 | 4.0 | 3.6 |
| Incheon | 2.0 | 6.5 | 5.8 |
| Jeju | 2.2 | 9.0 | 8.0 |

Source: Korea Meteorological Administration Wind Atlas 2.0 (KMA 2020).

---

**Table 4.** ABAQUS dynamic analysis results for the six tree cases.

| Case | H (m) | DBH (cm) | $f_1$ (Hz) | $\sigma_{\max,\,M_{\text{res}}}$ (MPa) | DAF |
|:-----|------:|----------:|-----------:|---------------------------------------:|----:|
| H6m-D15cm | 6 | 15 | 0.4953 | 70.30 | 3.905 |
| H6m-D20cm | 6 | 20 | 0.7591 | 19.37 | 1.076 |
| H6m-D25cm | 6 | 25 | 1.0560 | 8.73 | 0.485 |
| H8m-D15cm | 8 | 15 | 0.2786 | 138.75 | 7.708 |
| **H8m-D20cm** | **8** | **20** | **0.4271** | **54.25** | **3.014** |
| H8m-D25cm | 8 | 25 | 0.5942 | 23.27 | 1.293 |

$f_1$ = fundamental natural frequency; $\sigma_{\max}$ = maximum resultant bending stress at stem base; DAF = dynamic amplification factor (= dynamic peak stress / equivalent static stress under mean wind pressure). Bold row = reference case for parametric fatigue studies. Green-wood MOR of ginkgo ≈ 25.4 MPa; the H8m-D15cm case ($\sigma_{\max} = 138.75$ MPa) exceeds MOR, indicating static failure precedes fatigue accumulation.

---

**Table 5.** Fatigue life as a function of decay ratio (reference case: H = 8 m, DBH = 20 cm; ginkgo; Seoul wind regime).

| Decay ratio $r_d/R$ | 0.00 | 0.30 | 0.50 | 0.60 | 0.70 | 0.75 | 0.80 |
|:-------------------:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|
| Mean S-N (years) | 1,943 | 1,570 | 336 | 62 | 6.9 | 1.7 | 0.4 |
| Lower-bound S-N (years) | 59.7 | 53.2 | 25.5 | 9.5 | 1.7 | 0.5 | 0.1 |

Lower-bound S-N corresponds to the 95% confidence lower bound of the normalized S-N curve ($a_{lb} = 0.80$, $b_{lb} = 0.125$).

---

**Table 6.** Regional comparison of sound-tree fatigue life (reference case: H = 8 m, DBH = 20 cm; ginkgo; mean S-N).

| Site | $k$ | $c$ (m/s) | $\bar{U}$ (m/s) | Fatigue life (years) |
|:-----|:---:|:---------:|:---------------:|:--------------------:|
| Seoul | 1.6 | 4.0 | 3.6 | 1,943 |
| Incheon | 2.0 | 6.5 | 5.8 | 4 |
| Jeju | 2.2 | 9.0 | 8.0 | < 1 |

---

**Table 7.** Species comparison of fatigue life (reference geometry: H = 8 m, DBH = 20 cm; Seoul wind; mean S-N).

| Species | Green-wood MOR (MPa) | Sound-tree fatigue life (years) | Decayed ($r_d/R = 0.7$) fatigue life (years) |
|:--------|:--------------------:|:-------------------------------:|:--------------------------------------------:|
| Ginkgo (*Ginkgo biloba*) | 40.8 | 1,943 | 7 |
| Yoshino cherry (*Prunus × yedoensis*) | 52.8 | 2,800,788 | 2,757 |
| London planetree (*Platanus × acerifolia*) | 43.2 | 9,888 | 17 |
| Zelkova (*Zelkova serrata*) | 63.0 | > 10,000,000 | 779,401 |
| Dawn redwood (*Metasequoia glyptostroboides*) | 34.8 | 32 | 1 |

Green-wood MOR = 0.60 × dry-wood MOR (KS F 2207 2019). Structural geometry is held constant at the ginkgo reference case; results represent relative fatigue performance under equal geometric conditions.

---

**Table 8.** Effect of soil condition on basal stress correction factor and fatigue life (reference case: H = 8 m, DBH = 20 cm; ginkgo; Seoul wind; mean S-N).

| Soil condition | $k_s$ (kN/m³) | $C_\sigma$ | Fatigue life (years) |
|:---------------|:-------------:|:----------:|:--------------------:|
| Soft clay | 1,500 | 0.355 | >> 10,000 |
| Medium clay | 8,000 | 0.688 | 7,600 |
| Dense sand | 40,000 | 0.923 | 2,100 |
| Fully fixed | ∞ | 1.000 | 1,943 |

$k_s$ = coefficient of subgrade reaction; $C_\sigma$ = basal stress correction factor (Eq. 15); $L_r = 1.0$ m assumed for all cases.

---

## Figure Captions

**Figure 1.** Study flowchart of the integrated wind-induced fatigue life assessment framework.

**Figure 2.** ABAQUS finite element model: (a) tapered cantilever beam (H = 8 m, DBH = 20 cm) with biaxial wind loading ($F_X$, $F_Y$), fully fixed base, and key dimensions; (b) concentric-hollow decay cross-section model for three decay states ($r_d/R$ = 0, 0.5, 0.7) with corresponding stress amplification factors; (c) representative bending stress time history (80 s excerpt, H = 8 m, DBH = 20 cm, $\bar{U}$ = 20 m/s) and rainflow amplitude histogram for the full 600 s simulation.

**Figure 3.** Wind loading characterization: (a) Kaimal turbulence PSD with natural frequency markers for the six tree cases; (b) Weibull wind speed PDFs for Seoul, Incheon, and Jeju.

**Figure 4.** ABAQUS dynamic analysis results: (a) fundamental natural frequencies $f_1$; (b) maximum bending stress $\sigma_{\max}$ and dynamic amplification factor (DAF) for the six tree cases.

**Figure 5.** Normalized wood S-N curves: mean ($a$ = 0.87, $b$ = 0.10) and lower-bound ($a$ = 0.80, $b$ = 0.125) at $R$ = 0, with experimental data from Niklas (1992) and Clorius et al. (2000).

**Figure 6.** Fatigue life as a function of decay ratio $r_d/R$ for Seoul and Jeju wind regimes under mean and lower-bound S-N conditions; shaded bands indicate the uncertainty range; vertical dashed line marks the Mattheck critical ratio ($r_d/R$ = 0.70).

**Figure 7.** Species comparison of fatigue life for sound trees ($r_d/R$ = 0, solid bars) and decayed trees ($r_d/R$ = 0.7, hatched bars) under Seoul wind conditions (mean S-N).

**Figure 8.** Regional fatigue life analysis: (a) fatigue life comparison for Seoul, Incheon, and Jeju; (b) heat map of fatigue life as a function of decay ratio and wind site.

**Figure 9.** Influence of soil condition: (a) basal stress correction factor $C_\sigma$ versus subgrade reaction coefficient $k_s$; (b) corresponding mean S-N fatigue life for each soil category.

**Figure 10.** Conceptual failure mode diagram showing three domains as functions of Weibull scale parameter ($c$) and decay ratio ($r_d/R$): (I) high-cycle fatigue governed by long-term cumulative damage under routine winds; (II) decay-accelerated low-cycle fatigue where stress amplification by hollow cross-sections drives rapid damage accumulation; (III) static overturning where peak wind-induced base moment exceeds root anchorage capacity (approximate boundary corresponds to the 50-year return-period wind event producing a base rotation of $\theta \approx 5°$ for a fully fixed base). The domain boundaries are approximate and shift with tree geometry, species MOR, and soil conditions.
