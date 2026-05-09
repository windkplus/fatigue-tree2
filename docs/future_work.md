# Future Work: Wind-Induced Fatigue in Urban Street Trees

This document outlines research directions extending the present study on decay-accelerated, wind-induced fatigue of urban street trees. Each direction is motivated by an explicit limitation of the current framework or by an open question raised in the discussion.

---

## 1. Reversed-Loading (R ≈ −1) S-N Characterization for Tree-Relevant Species

**Motivation.** The S-N curves used in this study are derived from repeated-loading (R = 0) data, while wind-induced stem oscillation produces fully reversed bending (R ≈ −1). Yang et al. (2025) demonstrated that reversed loading is the most damaging regime for wood, yet published R = −1 fatigue data for the principal Korean street tree species are essentially absent.

**Proposed work.**
- Four-point reversed-bending tests on green-wood specimens of *Ginkgo biloba*, *Zelkova serrata*, *Platanus orientalis*, *Prunus serrulata*, and *Quercus palustris*.
- Comparative evaluation of Goodman, Gerber, and Soderberg mean-stress corrections against directly measured R = −1 data.
- Construction of a public S-N database for green-wood species under both repeated and reversed loading regimes.

**Expected impact.** Replaces upper-bound estimates with direct fatigue life predictions, removes the R-ratio non-conservatism, and provides the foundational dataset for tree-fatigue modelling generally.

---

## 2. Coupled Mechanical–Hydraulic Fatigue (Cyclic Embolism Pathway)

**Motivation.** Cyclic tensile microcracking documented by Yang et al. (2025) may compromise pit-membrane integrity in xylem vessels, lowering the threshold for embolism formation. This coupled mechanical–hydraulic fatigue mechanism is absent from conventional engineering fatigue models and represents the most distinctive feature of *living* tree fatigue.

**Proposed framework.**

```
σ_cycles ──► microcrack density ──► pit-membrane permeability ──► P50 shift ──► hydraulic failure
   │              │                       │                         │
   rainflow       Paris-law /              effective porosity        cavitation threshold
                  Weibull damage           model                     under drought
```

**Validation pathway.** Field instrumentation combining strain gauges with sap-flow sensors and stem psychrometers on trees of known decay/water-status during typhoon events. Collaboration with plant physiology groups (e.g., Korea Forest Research Institute) recommended.

**Expected impact.** First quantitative model linking storm-induced fatigue to drought-induced mortality, with implications for compound climate-event risk assessment.

---

## 3. Crown–Stem Multi-Body Dynamic Modelling

**Motivation.** The 20–50% reduction in DAF attributable to crown damping (James et al. 2006; Sellier and Fourcaud 2009) is currently an empirical correction. A bare-stem beam model—as used in this study—is conservative for fatigue but does not isolate the mechanism by which crowns dissipate energy.

**Proposed work.**
- Multi-body model with branched lumped-mass crown coupled to a tapered stem.
- Foliar drag formulated with Reynolds-dependent reconfiguration coefficient ψ(Re) following Vogel.
- Joint-compliance terms at major branch attachments.
- Calibration against published tree-shake experiments and on-site accelerometer arrays.

**Expected impact.** Replaces the conservative bare-stem assumption with physically grounded crown effects; quantifies the fatigue benefit of crown thinning vs. retention.

---

## 4. Probabilistic Eccentric Decay Geometry

**Motivation.** Real decay pockets are eccentric and irregular; the concentric hollow-cylinder assumption maximises section modulus and underestimates peak stress by 20–40% relative to asymmetric defects (Mattheck and Breloer 1994).

**Proposed work.**
- Treat the decayed cross-section as a random variable parameterised by centroid offset (eccentricity ratio), cavity aspect ratio, and circumferential wall-thickness distribution.
- Monte Carlo coupling with the existing rainflow–Miner pipeline to produce fatigue-life confidence intervals rather than point estimates.
- Calibration of geometric distributions from sonic tomography image archives.

**Expected impact.** Risk-based inspection thresholds keyed to *probabilistic* fatigue life rather than worst-case scenarios; directly actionable for municipal tree-management protocols.

---

## 5. Time-Evolving Decay + Fatigue Coupling (Remaining Useful Life)

**Motivation.** Fungal decay propagates at species-specific rates (typically 0.5–3 cm/yr radially for *Phellinus*, *Ganoderma*, and *Inonotus*). The present study treats decay as static; in reality, both the stress amplification factor and the fatigue accumulation rate evolve in time.

**Proposed formulation.**

$$D(t) = \int_0^t \dot{D}\left(r_d(\tau)\right)\, d\tau, \qquad r_d(\tau) = r_{d,0} + v_{\text{decay}}\,\tau$$

where the decay growth curve $r_d(\tau)$ is parameterised from forest-pathology literature.

**Expected impact.** Converts the static decay-ratio sweep into a Remaining Useful Life (RUL) estimator suitable for asset-management dashboards.

---

## 6. Graph Neural Network Surrogate for City-Scale Risk Assessment

**Motivation.** The full FEM–Kaimal–rainflow–Miner pipeline requires hours of compute per tree. City-scale fatigue assessment for a tree inventory of 10⁵–10⁶ specimens demands a surrogate model that retains physical fidelity while running in milliseconds per query.

**Proposed architecture.**
- Graph nodes: stem segments and branch junctions (taper, MOR, decay state).
- Global features: site Weibull (k, c), soil class, terrain roughness.
- Outputs: (DAF, peak stress, fatigue life) with Bayesian uncertainty intervals.
- Training data: 10⁴–10⁵ FEM realisations across the parameter envelope identified in this study.

**Expected impact.** Direct extension of the AGE Lab's CFD–GNN methodology to structural dynamics; enables real-time risk maps coupled to municipal tree inventories.

---

## 7. Non-Stationary Wind Regimes under Climate Change

**Motivation.** The Weibull (k, c) parameters used here assume stationarity. Under RCP 4.5/8.5 projections, typhoon frequency, intensity, and track latitude are all changing.

**Proposed work.**
- Time-varying Weibull parameters $(k(t), c(t))$ from CMIP6 / KMA downscaled wind projections.
- Annual-maxima extreme-value analysis (GEV) for tail behaviour.
- Compound-event scenarios (storm + drought + heatwave) for combined mechanical–hydraulic fatigue from §2.
- 2030 / 2050 / 2100 fatigue-life projections at major Korean cities.

**Expected impact.** First climate-projected fatigue-life atlas for Korean urban trees; informs species-selection guidelines for trees being planted *now* but expected to live through 2100.

---

## 8. Field Validation: Instrumented Tree Network

**Motivation.** No published study, including this one, validates wind-fatigue predictions against measured strain histories on living urban trees.

**Proposed campaign.**
- Cohort of 10–20 trees of varying species and known decay state (resistograph + sonic tomography baseline).
- Instrumentation:
  - Multi-axis strain gauges at base and DBH (rosette configuration).
  - Tri-axial accelerometers at crown centroid.
  - Anemometer at canopy height.
  - Optional: sap-flow sensors for §2 hydraulic coupling.
- 3–5 year monitoring window to capture ≥1 typhoon event per site.
- Coordination with Seoul Forest, National Institute of Forest Science.

**Expected impact.** Provides the strain-history validation dataset that has been the missing experimental anchor for tree-biomechanics fatigue models for two decades.

---

## 9. Adaptive Growth and Reaction-Wood Self-Repair

**Motivation.** Living trees redistribute material in response to mechanical loading (thigmomorphogenesis) and form reaction wood under sustained asymmetric loads. Conventional fatigue treats the cross-section as fixed; trees do not.

**Proposed direction.**
- Coupled growth–fatigue model with cross-section evolving as a function of recent peak-stress history.
- Empirical growth–stress relationships from wind-stake / pendulation experiments.
- Identifies whether mature urban trees converge to a "fatigue-resistant" geometry over decades.

**Expected impact.** Distinguishes engineering fatigue (passive degradation) from biological fatigue (active adaptation); reframes pruning and staking practice in terms of the adaptive response they provoke.

---

## 10. Multi-Scale Damage Mechanics: Cell Wall to Whole Tree

**Motivation.** Wood fatigue resistance is ultimately governed at the cell-wall scale by microfibril angle (MFA): high-MFA juvenile wood has lower stiffness but higher strain-to-failure.

**Proposed work.**
- MFA distributions from synchrotron x-ray diffraction across species and growth-ring positions.
- Cell-wall finite element models with appropriate yield criterion (Hill / Tsai-Wu).
- Hierarchical homogenisation linking cell-wall to clear-wood to defect-bearing-stem fatigue behaviour.

**Expected impact.** Mechanistic explanation for the order-of-magnitude fatigue-life differences observed across species in this study; guides material-level breeding/selection criteria for urban forestry.

---

## 11. Decision-Support Tool for Urban Foresters

**Motivation.** Research outputs that remain in journal articles do not influence management. A practitioner-accessible interface is the final translational step.

**Proposed deliverables.**
- Web tool: input (species, DBH, height, decay ratio, site) → output (fatigue life CI, recommended inspection interval, comparison against alternative species).
- Integration with municipal tree inventories (Seoul Tree Atlas; Statistics Korea forestry datasets).
- Risk-based inspection-interval recommendation engine grounded in §4 and §5 above.

**Expected impact.** Converts the fatigue framework from a research artefact into a deployed urban-forestry asset-management standard.

---

## Suggested Priority Order (2026–2029)

| Rank | Direction | Rationale | Lead time |
|------|-----------|-----------|-----------|
| 1 | §1 R = −1 S-N experiments | Foundational; unblocks every downstream prediction | 1–2 yr |
| 2 | §4 Probabilistic eccentric decay | Low-cost extension of existing pipeline | < 1 yr |
| 3 | §6 GNN surrogate | Leverages AGE Lab's existing GNN expertise | 1 yr |
| 4 | §8 Field validation network | Long lead time — must start early | 3–5 yr |
| 5 | §2 Mechanical–hydraulic coupling | Highest novelty; cross-disciplinary | 2–3 yr |
| 6 | §7 Climate-projected fatigue atlas | Builds on §1 + §4 | 2 yr after §1 |

---

## Cross-Cutting Themes

Three themes recur across the directions above and may merit dedicated framing in subsequent proposals:

1. **From upper-bound to expected-value.** The present study is deliberately conservative (R = 0 data, concentric decay, bare stem). Each future direction relaxes one of these conservatisms.
2. **From static to time-evolving.** Decay (§5), wind climate (§7), and tree geometry (§9) are all currently treated as static. A unified time-domain framework would bring tree fatigue assessment in line with civil-infrastructure asset management.
3. **From engineering fatigue to biological fatigue.** The coupling of mechanical damage with hydraulic failure (§2) and adaptive growth (§9) defines the unique research frontier that distinguishes tree fatigue from inert-material fatigue.
