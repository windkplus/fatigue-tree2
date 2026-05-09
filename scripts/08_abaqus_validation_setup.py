"""
08_abaqus_validation_setup.py
-----------------------------
Phase F-1 of paper #2: build six ABAQUS 2D plane-stress (CPS4R) input files
that mesh an eccentric circular cross-section. ABAQUS will be used as an
independent numerical witness for the closed-form section properties used
throughout the paper.

Validation matrix:
    r_d/R   :  0.5, 0.7
    e_norm  :  0.0, 0.5, 0.95
    -> 6 cases

Each .inp file:
  - Defines a 2D structured Cartesian mesh of CPS4R elements with side
    `dx = R / N_GRID` (default N_GRID = 80, giving ~6,000 elements per case),
    keeping only elements whose centroid lies inside the outer circle and
    outside the cavity.
  - Includes a single trivial static step with no external load. This forces
    ABAQUS to write an .odb whose part-instance carries the geometric mass
    properties (area, centroid, inertia) used by `09_extract_validation.py`.
  - Uses density = 1.0 and unit thickness so that the reported mass equals
    the cross-sectional area, the reported centroid equals the area centroid,
    and the reported inertia tensor (about the centroid) equals the area
    moments of inertia.

After all six jobs run, `09_extract_validation.py` queries each .odb via
ABAQUS Python and writes a CSV of the FEM-computed properties; the comparison
script `10_compare_validation.py` then juxtaposes them with the closed-form
predictions and generates the validation figure for the paper appendix.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from eccentric_decay.geometry import DecaySection  # noqa: E402

# ─────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────
R = 1.0                         # outer radius [m] (unit-radius)
N_GRID = 80                     # half-grid resolution in each axis
                                # actual elements ~ pi/4 * (2 * N_GRID)^2 ~ 20k
THICKNESS = 1.0                 # element thickness [m]
DENSITY = 1.0                   # so that mass = area, I_inertial = I_area

CASES = []
for r_dR in [0.5, 0.7]:
    for e_norm in [0.0, 0.5, 0.95]:
        CASES.append(dict(name=f"val_rd{int(r_dR*100):02d}_e{int(e_norm*100):03d}",
                           r_dR=r_dR, e_norm=e_norm))

OUT_DIR = ROOT / "data" / "abaqus_validation"
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ─────────────────────────────────────────────────────────────────────
# Mesh + inp writer
# ─────────────────────────────────────────────────────────────────────
def build_mesh(R: float, r_d: float, e_x: float, e_y: float, n: int):
    """
    Build a structured Cartesian mesh of CPS4R elements covering a circular
    disk of radius R minus an offset circular cavity (e_x, e_y, r_d).

    Returns
    -------
    nodes : (N, 2) array of node coordinates
    elems : (M, 4) array of node ids (1-indexed, ABAQUS convention)
    """
    side = 2.0 * R
    dx = side / (2 * n)               # uniform cell side length
    # candidate node grid (vertices of cells)
    xs = np.linspace(-R, R, 2 * n + 1)
    ys = np.linspace(-R, R, 2 * n + 1)
    # element centroids
    xc = 0.5 * (xs[:-1] + xs[1:])
    yc = 0.5 * (ys[:-1] + ys[1:])
    Xc, Yc = np.meshgrid(xc, yc, indexing="xy")
    inside_outer = (Xc ** 2 + Yc ** 2) <= R ** 2
    if r_d > 0:
        inside_cav = ((Xc - e_x) ** 2 + (Yc - e_y) ** 2) <= r_d ** 2
    else:
        inside_cav = np.zeros_like(inside_outer, dtype=bool)
    keep = inside_outer & ~inside_cav
    # build node ids only for nodes touched by kept elements
    used_node = np.zeros(len(xs) * len(ys), dtype=bool)
    elems = []
    for i in range(2 * n):                   # column index
        for j in range(2 * n):               # row index
            if not keep[j, i]:
                continue
            # CCW node ordering
            n0 = j * len(xs) + i              # bottom-left
            n1 = j * len(xs) + (i + 1)        # bottom-right
            n2 = (j + 1) * len(xs) + (i + 1)  # top-right
            n3 = (j + 1) * len(xs) + i        # top-left
            elems.append((n0, n1, n2, n3))
            used_node[[n0, n1, n2, n3]] = True
    elems = np.asarray(elems, dtype=np.int64)
    # compact node numbering
    new_id = -np.ones(len(xs) * len(ys), dtype=np.int64)
    new_id[used_node] = np.arange(used_node.sum())
    elems = new_id[elems]
    # build node coordinate array
    XX, YY = np.meshgrid(xs, ys, indexing="xy")
    flat_x = XX.ravel(order="C")
    flat_y = YY.ravel(order="C")
    # we ravelled with j (row) outer, i (col) inner, matching n_id = j*len(xs)+i? Check
    # Actually: meshgrid(xs, ys, indexing="xy") gives shape (len(ys), len(xs))
    # So XX[j, i] corresponds to xs[i], ys[j] -- correct.
    nodes = np.stack([flat_x[used_node], flat_y[used_node]], axis=-1)
    return nodes, elems + 1                   # ABAQUS 1-indexed


# ─────────────────────────────────────────────────────────────────────
# .inp template
# ─────────────────────────────────────────────────────────────────────
INP_TEMPLATE = """*Heading
 fatigue-tree2 paper #2 Appendix A: ABAQUS validation case {name}
 r_d/R = {r_dR:.2f}, e_norm = {e_norm:.2f}
 Generated by 08_abaqus_validation_setup.py
**
*Part, name=SECTION
*Node
{nodes}
*Element, type=CPS4R, elset=ALL
{elems}
*Solid Section, elset=ALL, material=WOOD
{thickness:.4f}
*End Part
**
*Material, name=WOOD
*Density
{density:.4f}
*Elastic
9.5e9, 0.3
**
*Assembly, name=ASSEMBLY
*Instance, name=SEC-1, part=SECTION
*End Instance
*End Assembly
**
*Step, name=trivial, perturbation=NO
*Static
1.0, 1.0, 1e-5, 1.0
*Output, field
*Element Output, elset=SEC-1.ALL
EVOL, IVOL, COORD
*End Step
"""


def write_node_block(nodes: np.ndarray) -> str:
    return "\n".join(
        f" {i+1}, {x:.10f}, {y:.10f}"
        for i, (x, y) in enumerate(nodes))


def write_element_block(elems: np.ndarray) -> str:
    return "\n".join(
        f" {i+1}, {n1}, {n2}, {n3}, {n4}"
        for i, (n1, n2, n3, n4) in enumerate(elems))


# ─────────────────────────────────────────────────────────────────────
# Loop over cases
# ─────────────────────────────────────────────────────────────────────
print(f"Generating {len(CASES)} ABAQUS 2D validation cases at {OUT_DIR}\n")
print(f"{'case':<24} {'r_d/R':>6} {'e_norm':>7} {'nodes':>8} {'elems':>8} "
      f"{'A_mesh':>10} {'A_exact':>10} {'rel.err':>8}")
print("-" * 96)

for case in CASES:
    name = case["name"]
    r_dR = case["r_dR"]
    e_norm = case["e_norm"]
    r_d = r_dR * R
    e_radial = e_norm * (R - r_d)
    e_x = e_radial             # offset along +x
    e_y = 0.0
    nodes, elems = build_mesh(R, r_d, e_x, e_y, n=N_GRID)
    A_mesh = (2.0 * R / (2 * N_GRID)) ** 2 * len(elems)
    sec = DecaySection.eccentric_circular(R, r_dR, ecc_norm=e_norm)
    A_exact = sec.section_properties()["A"]
    rel_A = abs(A_mesh - A_exact) / A_exact
    inp_text = INP_TEMPLATE.format(
        name=name, r_dR=r_dR, e_norm=e_norm,
        nodes=write_node_block(nodes),
        elems=write_element_block(elems),
        thickness=THICKNESS, density=DENSITY,
    )
    out_inp = OUT_DIR / f"{name}.inp"
    out_inp.write_text(inp_text, encoding="utf-8")
    print(f"{name:<24} {r_dR:>6.2f} {e_norm:>7.2f} "
          f"{len(nodes):>8d} {len(elems):>8d} "
          f"{A_mesh:>10.5f} {A_exact:>10.5f} {rel_A:>8.2e}")


# ─────────────────────────────────────────────────────────────────────
# Batch driver
# ─────────────────────────────────────────────────────────────────────
batch_path = OUT_DIR / "run_validation.bat"
extract_calls = "\n".join(
    f"call abaqus python extract_validation.py {c['name']}" for c in CASES)
job_calls = "\n".join(
    f"echo Running {c['name']}\ncall abaqus job={c['name']} cpus=1 "
    f"interactive ask_delete=off"
    for c in CASES)
batch = f"""@echo off
echo ============================================================
echo   ABAQUS 2D Validation Run (paper #2 Appendix A)
echo ============================================================
cd /d "{OUT_DIR.as_posix()}"

{job_calls}

echo.
echo ============================================================
echo   Extracting mass properties from each .odb
echo ============================================================

{extract_calls}

echo.
echo ============================================================
echo   Validation runs complete. See validation_results.csv
echo ============================================================
"""
batch_path.write_text(batch, encoding="utf-8")
print(f"\nBatch driver written: {batch_path}")
print(f"Run with:  cd {OUT_DIR}  &&  run_validation.bat\n")
