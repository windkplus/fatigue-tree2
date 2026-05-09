"""
wood_sn_plot.py
---------------
목재 S-N 곡선 시각화 및 수종 간 비교
wood_properties_db.py 의 데이터 활용
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from wood_properties_db import STATIC_PROPS, SN_PARAMS, MC_CORRECTION, sn_Nf

matplotlib.rcParams['font.family'] = 'Malgun Gothic'
matplotlib.rcParams['axes.unicode_minus'] = False

out_dir = __import__('pathlib').Path("output")
out_dir.mkdir(exist_ok=True)

N_plot = np.logspace(1, 8, 500)

fig, axes = plt.subplots(1, 3, figsize=(17, 6))
fig.suptitle("가로수 주요 수종 목재 S-N 곡선 비교\n"
             "(σ_a/MOR = a − b·log₁₀(N_f),  소건재 MC≈12%)",
             fontsize=11)

COLORS = ["steelblue", "hotpink", "darkorange", "forestgreen",
          "purple", "saddlebrown"]

# ── (a) R=0 편진: 기건 ────────────────────────────────────────────
ax = axes[0]
for (sp, params), color in zip(SN_PARAMS.items(), COLORS):
    MOR = STATIC_PROPS[sp]["MOR_MPa"]
    a, b = params["a_R0"], params["b_R0"]
    sigma_sn = MOR * (a - b * np.log10(N_plot))
    valid = (sigma_sn > 0) & (sigma_sn <= MOR * a)
    # 피로한계 수평선
    sigma_lim = MOR * params["sigma_limit_ratio"]
    ax.semilogx(N_plot[valid], sigma_sn[valid], lw=1.8,
                color=color, label=sp.split("(")[0].strip())
    ax.axhline(sigma_lim, color=color, lw=0.8, ls=":", alpha=0.6)
ax.set_xlabel("파단 반복수 N_f")
ax.set_ylabel("응력 진폭 σ_a [MPa]")
ax.set_title("(a) R=0 (편진, 풍하중 근사)\n기건재 기준")
ax.legend(fontsize=7, loc="upper right")
ax.grid(True, which="both", alpha=0.3)
ax.set_xlim([1e1, 1e8])
ax.set_ylim([0, 120])

# ── (b) R=-1 완전반전: 기건 ──────────────────────────────────────
ax = axes[1]
for (sp, params), color in zip(SN_PARAMS.items(), COLORS):
    MOR = STATIC_PROPS[sp]["MOR_MPa"]
    a, b = params["a_Rm1"], params["b_Rm1"]
    sigma_sn = MOR * (a - b * np.log10(N_plot))
    valid = (sigma_sn > 0) & (sigma_sn <= MOR * a)
    sigma_lim = MOR * params["sigma_limit_ratio"]
    ax.semilogx(N_plot[valid], sigma_sn[valid], lw=1.8,
                color=color, label=sp.split("(")[0].strip())
    ax.axhline(sigma_lim, color=color, lw=0.8, ls=":", alpha=0.6)
ax.set_xlabel("파단 반복수 N_f")
ax.set_ylabel("응력 진폭 σ_a [MPa]")
ax.set_title("(b) R=−1 (완전 반전)\n기건재 기준")
ax.legend(fontsize=7, loc="upper right")
ax.grid(True, which="both", alpha=0.3)
ax.set_xlim([1e1, 1e8])
ax.set_ylim([0, 120])

# ── (c) 생재 보정 효과 (R=0, 은행나무) ──────────────────────────
ax = axes[2]
sp = "은행나무 (Ginkgo biloba)"
params = SN_PARAMS[sp]
a, b = params["a_R0"], params["b_R0"]

MOR_dry   = STATIC_PROPS[sp]["MOR_MPa"]
MOR_green = MOR_dry * MC_CORRECTION["MOR"]

for MOR_val, label, ls, color in [
    (MOR_dry,   f"기건 (MOR={MOR_dry:.0f} MPa)",  "-",  "steelblue"),
    (MOR_green, f"생재 (MOR={MOR_green:.0f} MPa)", "--", "darkorange"),
]:
    sigma_sn = MOR_val * (a - b * np.log10(N_plot))
    valid = (sigma_sn > 0) & (sigma_sn <= MOR_val * a)
    ax.semilogx(N_plot[valid], sigma_sn[valid], lw=2.0,
                color=color, ls=ls, label=label)
    sigma_lim = MOR_val * params["sigma_limit_ratio"]
    ax.axhline(sigma_lim, color=color, lw=1.0, ls=":", alpha=0.7,
               label=f"피로한계 {sigma_lim:.0f} MPa")

# 현재 해석 응력 범위 표시
ax.axhspan(0, 18, alpha=0.08, color="green",
           label="현재 해석 응력범위 (0~18 MPa)")
ax.set_xlabel("파단 반복수 N_f")
ax.set_ylabel("응력 진폭 σ_a [MPa]")
ax.set_title(f"(c) 생재 vs 기건 비교\n({sp.split('(')[0].strip()}, R=0)")
ax.legend(fontsize=7.5, loc="upper right")
ax.grid(True, which="both", alpha=0.3)
ax.set_xlim([1e1, 1e8])
ax.set_ylim([0, 80])

plt.tight_layout()
fig_path = out_dir / "wood_sn_curves.png"
plt.savefig(fig_path, dpi=150, bbox_inches="tight")
plt.show()
print(f"S-N 곡선 그래프 저장: {fig_path}")
