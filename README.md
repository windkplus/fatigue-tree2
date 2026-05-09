# fatigue-tree2

> Probabilistic Eccentric Decay Geometry — Monte Carlo Confidence Intervals on the Wind-Fatigue Life of Urban Street Trees

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20097384.svg)](https://doi.org/10.5281/zenodo.20097384)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)

Companion repository to the manuscript *Probabilistic Eccentric Decay Geometry: Monte Carlo Confidence Intervals on the Wind-Fatigue Life of Urban Street Trees* (Kang, submitted to *Trees – Structure and Function*, 2026), and the second instalment of a planned two-part series. The first paper (Kang 2026, *Trees / Int. J. Fatigue*) established a deterministic absolute-life pipeline using the classical concentric hollow-cylinder decay model of Mattheck (1994); the present project relaxes that assumption and quantifies how cavity-geometry uncertainty propagates into the predicted fatigue life via Monte Carlo simulation.

## Headline result

At the Mattheck threshold $r_d/R = 0.7$, the *median* eccentric tree has a fatigue life of **25 % of the concentric estimate** (lower-bound S-N), and the 5th-percentile tail is shorter by a factor of **47** (lower-bound S-N) or **143** (mean S-N). The eccentric correction is robust across wind-direction distributions; cavity-geometry tomography, not site-specific wind data, is the dominant empirical uncertainty.

## Repository layout

```
fatigue-tree2/
├── src/eccentric_decay/    # core analytical framework (Phase A)
│   ├── geometry.py
│   ├── stress_amplification.py
│   ├── wind_directions.py
│   └── monte_carlo.py
├── src/core/               # modules inherited from Paper I (read-only)
├── scripts/                # reproducible analysis pipeline (00_…, 01_…, …)
├── data/
│   ├── abaqus_validation/  # ABAQUS .inp + extraction script (Appendix A)
│   ├── inherited/          # Paper I data used for anchoring (read-only)
│   └── generated/          # Monte Carlo outputs (CSV.gz)
├── output/
│   ├── figures/            # publication PNGs (Figs 1–11, A1)
│   └── tables/             # CSV summaries (mc_baseline_summary, …)
├── CITATION.cff
├── LICENSE
├── requirements.txt
└── README.md  ← this file
```

## Quick start

```powershell
# 1. Set up the environment
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt

# 2. Verify the framework against Mattheck (1994) — must pass before anything else
python scripts/00_regression_test_concentric.py

# 3. Reproduce the figures and tables of the paper, in order
python scripts/01_demo_eccentric_saf.py          # Fig. 2
python scripts/02_monte_carlo_baseline.py        # Tables 1; Figs 3–6
python scripts/03_visualize_mc.py
python scripts/04_anchor_to_paper1.py            # Table 2; Figs 7–9
python scripts/05_sensitivity.py                 # Fig. 10
python scripts/06_policy_matrix.py               # Table 3; Fig. 11
python scripts/07_regression_viz.py              # Fig. 1

# 4. Independent ABAQUS validation (requires ABAQUS Standard 2026; see Appendix A)
python scripts/08_abaqus_validation_setup.py     # generate 6 .inp files
#   then run data/abaqus_validation/run_validation.bat in the ABAQUS environment
python scripts/10_compare_validation.py          # Fig. A1

# 5. Build the manuscript and the cover letter
python build_paper2_docx.py
python build_cover_letter.py
```

## Methodology summary

| Phase | Output | Script(s) |
| ----- | ------ | --------- |
| A. Closed-form eccentric / elliptical SAF | `src/eccentric_decay/{geometry,stress_amplification}.py` | `00, 01, 07` |
| B. Monte Carlo over cavity geometry | `data/generated/mc_baseline_*.csv.gz` (300k samples) | `02, 03` |
| C1. Anchoring to Paper I absolute lives | `output/tables/absolute_life_summary.csv` | `04` |
| C2. Sensitivity (priors, aspect, Sobol) | `output/tables/sensitivity_*.csv` | `05` |
| Policy. Risk-based inspection matrix | `output/tables/policy_matrix.csv` | `06` |
| F. ABAQUS independent validation | `output/tables/abaqus_validation_comparison.csv` | `08, 10` |
| Paper. Trees-Springer DOCX + cover letter | `output/paper/paper2_Trees.docx`, `cover_letter_paper2.docx` | `build_*` |

## Reproducing the Monte Carlo

All Monte Carlo runs use a fixed seed (42) and per-cell offsets so the published numbers are bit-reproducible across machines.

```powershell
python scripts/00_regression_test_concentric.py   # PASS at machine precision
python scripts/02_monte_carlo_baseline.py         # ≈ 4 min on a laptop
```

The 300 k-sample CSV.gz files in `data/generated/` are checked into the repository as ground truth for downstream visualisations.

## Independent ABAQUS validation (Appendix A)

Six 2D plane-stress (CPS4R, ~ 15 000 elements per case) ABAQUS analyses verify the closed-form section properties to better than 0.1 % across all 36 case-quantity combinations. The input files, batch driver, and post-processing script are in `data/abaqus_validation/`. Reproducing this validation requires ABAQUS Standard 2026 (Dassault Systèmes) with a standard academic licence.

## How to cite

If you use this software or data, please cite both the paper and the archive:

```bibtex
@article{Kang2026EccentricDecay,
  author  = {Kang, Junsuk},
  title   = {Probabilistic Eccentric Decay Geometry: Monte Carlo Confidence Intervals on the Wind-Fatigue Life of Urban Street Trees},
  journal = {Trees - Structure and Function},
  year    = {2026},
  note    = {Paper II of a two-part series; Paper I: Fatigue Life Assessment of Urban Street Trees under Wind Loading}
}

@software{fatigue_tree2_2026,
  author       = {Kang, Junsuk},
  title        = {fatigue-tree2: Probabilistic eccentric-decay Monte Carlo for wind-fatigue of urban street trees},
  year         = 2026,
  version      = {v1.0.0},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.20097384},
  url          = {https://doi.org/10.5281/zenodo.20097384}
}
```

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20097384.svg)](https://doi.org/10.5281/zenodo.20097384)

## Companion paper repository

Paper I is archived at [github.com/windkplus/fatigue-tree](https://github.com/windkplus/fatigue-tree).

## Author and contact

**Junsuk Kang**, Department of Landscape Architecture and Rural Systems Engineering, Seoul National University. <junkang@snu.ac.kr>

## License

MIT (see [LICENSE](LICENSE)).
