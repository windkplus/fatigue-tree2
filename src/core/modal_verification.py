"""
modal_verification.py
---------------------
ABAQUS 모달해석 결과 후처리 및 고유진동수 검증

기능:
  1. 해석적 방법(Euler-Bernoulli beam)으로 고유진동수 계산
     → ABAQUS 결과 없어도 즉시 실행 가능
  2. ABAQUS .dat / 수동 입력 결과와 비교
  3. Davenport/Kaimal PSD와 오버레이 → 공진 위험도 판단
  4. 파라미터(DBH, H) 변화에 따른 f₁ 분포 → 공진 조건 지도

참고:
  Sellier & Fourcaud (2009). Plant Sci., 177, 33-42.
  Moore & Maguire (2005). Trees, 19, 432-442.
  James et al. (2006). Am. J. Bot., 93, 1522-1530.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path
from wood_properties_db import STATIC_PROPS, DYNAMIC_REF

matplotlib.rcParams['font.family'] = 'Malgun Gothic'
matplotlib.rcParams['axes.unicode_minus'] = False

out_dir = Path("output")
out_dir.mkdir(exist_ok=True)

# ══════════════════════════════════════════════════════════════════════
# 1. 해석적 고유진동수 계산 (Euler-Bernoulli 외팔보)
# ══════════════════════════════════════════════════════════════════════
# 균일 원형 단면 외팔보의 n차 고유진동수:
#
#   f_n = (β_n·L)² / (2π·L²) × √(E·I / (ρ·A))
#
# β_n·L 값 (외팔보):
#   1차: 1.8751,  2차: 4.6941,  3차: 7.8548,  4차: 10.9955
#
# 테이퍼 보정: 선형 taper를 가진 보의 f₁은 균일 단면보다 약 10~20% 높음
# 보정계수 α_taper ≈ 1 + 0.5 × (D_base - D_top) / D_base

BETA_L = np.array([1.8751, 4.6941, 7.8548, 10.9955])  # 1~4차

def natural_freq_uniform(E_GPa, rho_kgm3, D_m, H_m, n_modes=4):
    """
    균일 원형 단면 외팔보의 n차 고유진동수 [Hz]
    E_GPa   : 종방향 탄성계수
    rho     : 밀도 [kg/m³]  (수관 질량 미포함)
    D_m     : 줄기 직경 [m]  (근원부 기준)
    H_m     : 수고 [m]
    """
    E   = E_GPa * 1e9            # [Pa]
    I   = np.pi * D_m**4 / 64   # [m⁴]
    A   = np.pi * D_m**2 / 4    # [m²]
    EI  = E * I
    rhoA = rho_kgm3 * A

    freqs = []
    for bl in BETA_L[:n_modes]:
        f = (bl**2 / (2 * np.pi * H_m**2)) * np.sqrt(EI / rhoA)
        freqs.append(f)
    return np.array(freqs)


def taper_correction(D_base, D_top):
    """
    선형 taper 보정계수 (1차 고유진동수 기준)
    실험식 기반 근사 (Sellier & Fourcaud 2009)
    """
    return 1.0 + 0.5 * (D_base - D_top) / D_base


def natural_freq_with_crown(E_GPa, rho_kgm3, D_base_m, H_m,
                             M_crown_kg, taper=0.5):
    """
    수관 집중질량을 포함한 보정 고유진동수
    M_crown : 수관 질량 [kg]
    taper   : 직경 감소 비율 D_top/D_base = 1 - taper

    Dunkerley의 방법으로 수관 질량 효과 보정:
      1/f₁_total² = 1/f₁_beam² + 1/f₁_crown²
    """
    D_top = D_base_m * (1 - taper)
    alpha = taper_correction(D_base_m, D_top)

    freqs_beam = natural_freq_uniform(E_GPa, rho_kgm3, D_base_m, H_m)
    f1_beam    = freqs_beam[0] * alpha    # taper 보정

    # 수관 집중질량에 의한 고유진동수 (외팔보 선단 집중질량)
    E   = E_GPa * 1e9
    I   = np.pi * D_base_m**4 / 64
    k_tip = 3 * E * I / H_m**3          # 선단 강성 [N/m]
    f1_crown = np.sqrt(k_tip / M_crown_kg) / (2 * np.pi)

    # Dunkerley 결합
    f1_total = 1.0 / np.sqrt(1/f1_beam**2 + 1/f1_crown**2)
    return f1_total, f1_beam, f1_crown, freqs_beam


# ══════════════════════════════════════════════════════════════════════
# 2. 수종별 고유진동수 계산 및 출력
# ══════════════════════════════════════════════════════════════════════
# 기준 조건
H_BASE   = 8.0     # [m]
DBH      = 0.20    # [m]
D_BASE   = DBH * 1.2
M_CROWN  = 150.0   # [kg]
TAPER    = 0.5     # D_top = D_base × 0.5

print("=" * 65)
print("  수종별 고유진동수 계산 (H=8m, DBH=20cm, M_crown=150kg)")
print("=" * 65)
print(f"  {'수종':<22} {'f₁(beam)':>9} {'f₁(총)':>8} {'f₂':>8} {'문헌범위':>14}")
print(f"  {'':22} {'Hz':>9} {'Hz':>8} {'Hz':>8} {'Hz':>14}")
print("-" * 65)

ref_range = DYNAMIC_REF["일반 가로수 (H=5~10m)"]["f1_Hz"]
species_f1 = {}

for sp, props in STATIC_PROPS.items():
    # 생재 탄성계수 사용 (MC 보정)
    E_green = props["E_L_GPa"] * 0.68
    rho     = props["density_green"]

    f1_tot, f1_beam, f1_crown, freqs_all = natural_freq_with_crown(
        E_green, rho, D_BASE, H_BASE, M_CROWN, TAPER)
    f2 = freqs_all[1] * taper_correction(D_BASE, D_BASE*(1-TAPER))

    in_range = "✓" if ref_range[0] <= f1_tot <= ref_range[1] else "△"
    short = sp.split("(")[0].strip()[:18]
    print(f"  {short:<22} {f1_beam:>9.3f} {f1_tot:>8.3f} {f2:>8.3f} "
          f"  {ref_range[0]}~{ref_range[1]} Hz {in_range}")
    species_f1[sp] = f1_tot

print("-" * 65)
print(f"  문헌 기준 범위: {ref_range[0]}~{ref_range[1]} Hz  (Moore 2005, Sellier 2009)")

# ══════════════════════════════════════════════════════════════════════
# 3. DBH × 수고 → f₁ 공진 조건 지도
# ══════════════════════════════════════════════════════════════════════
print("\n  f₁ 공진 조건 지도 계산 중...")

DBH_range = np.linspace(0.10, 0.40, 25)
H_range   = np.linspace(4.0, 14.0, 25)

# 은행나무 생재 기준
sp_ref = "은행나무 (Ginkgo biloba)"
E_ref  = STATIC_PROPS[sp_ref]["E_L_GPa"] * 0.68
rho_ref = STATIC_PROPS[sp_ref]["density_green"]

f1_map = np.zeros((len(H_range), len(DBH_range)))
for i, H in enumerate(H_range):
    for j, D in enumerate(DBH_range):
        D_b = D * 1.2
        M_c = 150 * (H / 8.0) * (D / 0.20)  # 수관 질량 비례 가정
        f1, *_ = natural_freq_with_crown(E_ref, rho_ref, D_b, H, M_c, TAPER)
        f1_map[i, j] = f1

# ══════════════════════════════════════════════════════════════════════
# 4. ABAQUS 결과 비교 인터페이스
# ══════════════════════════════════════════════════════════════════════
# ABAQUS .dat 파일에서 고유진동수를 추출하는 함수
# (실행 전 ABAQUS 미보유 시 샘플 값으로 동작)

def load_abaqus_frequencies(dat_path):
    """
    ABAQUS .dat 파일에서 고유진동수 파싱
    출력 형식 예:
      MODE NO    EIGENVALUE         FREQUENCY (CYCLES/TIME)
         1       9.86960E+00        5.00000E-01
    """
    freqs = []
    try:
        with open(dat_path, "r") as f:
            for line in f:
                if "FREQUENCY (CYCLES/TIME)" in line:
                    continue
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        mode = int(parts[0])
                        freq = float(parts[2])
                        freqs.append((mode, freq))
                    except ValueError:
                        continue
    except FileNotFoundError:
        pass
    return freqs

# 샘플 ABAQUS 결과 (실제 해석 전 플레이스홀더)
ABAQUS_SAMPLE = {
    "은행나무 H=8m DBH=20cm (시뮬레이션)": {
        "freqs_Hz": [0.41, 2.58, 7.21, 14.13],
        "color": "steelblue",
        "marker": "o"
    },
}

# ══════════════════════════════════════════════════════════════════════
# 5. 시각화
# ══════════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(16, 11))
fig.suptitle("가로수 고유진동수 검증 및 공진 위험도 분석", fontsize=12)
gs = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.35)

# ── 5-a. PSD + 수종별 f₁ 위치 ─────────────────────────────────────
ax1 = fig.add_subplot(gs[0, :2])

# Kaimal PSD (U=20m/s, z=6m)
z0, I_u, U_bar = 0.3, 0.20, 20.0
alpha_ = 0.67 + 0.05 * np.log(z0)
L_u    = 300.0 * (6.0 / 200.0)**alpha_
sigma_u = I_u * U_bar
f_plot  = np.logspace(-2, 1, 500)
f_L     = f_plot * L_u / U_bar
psd_kai = sigma_u**2 * 6.8 * f_L / (f_plot * (1 + 10.2*f_L)**(5/3))

ax1.fill_between(f_plot, psd_kai, alpha=0.15, color="steelblue")
ax1.loglog(f_plot, psd_kai, lw=2, color="steelblue", label="Kaimal PSD (U=20m/s, z=6m)")

# 공진 위험 범위 강조
ax1.axvspan(0.20, 0.80, alpha=0.20, color="red", label="수목 1차 공진 대역 (0.2~0.8 Hz)")

COLORS_SP = ["steelblue","hotpink","darkorange","forestgreen","purple","saddlebrown"]
for (sp, f1), color in zip(species_f1.items(), COLORS_SP):
    short = sp.split("(")[0].strip()
    psd_at_f1 = (sigma_u**2 * 6.8 * (f1*L_u/U_bar) /
                 (f1 * (1 + 10.2*(f1*L_u/U_bar))**(5/3)))
    ax1.axvline(f1, color=color, lw=1.2, ls="--", alpha=0.8)
    ax1.scatter(f1, psd_at_f1, color=color, s=60, zorder=5,
                label=f"{short} f₁={f1:.3f}Hz")

ax1.set_xlabel("주파수 [Hz]")
ax1.set_ylabel("PSD [m²/s²/Hz]")
ax1.set_title("Kaimal 풍속 PSD 위의 수종별 1차 고유진동수\n(점이 공진 대역 내에 있을수록 피로 위험 높음)")
ax1.legend(fontsize=7, ncol=2, loc="upper right")
ax1.grid(True, which="both", alpha=0.3)
ax1.set_xlim([0.01, 10])

# ── 5-b. f₁ 수종 비교 막대 ────────────────────────────────────────
ax2 = fig.add_subplot(gs[0, 2])
sp_short = [s.split("(")[0].strip() for s in species_f1]
f1_vals  = list(species_f1.values())
bars = ax2.barh(sp_short, f1_vals, color=COLORS_SP[:len(f1_vals)], alpha=0.8)
ax2.axvspan(0.20, 0.80, alpha=0.20, color="red", label="공진 위험 대역")
for bar, f1 in zip(bars, f1_vals):
    ax2.text(f1 + 0.01, bar.get_y() + bar.get_height()/2,
             f"{f1:.3f} Hz", va="center", fontsize=8)
ax2.set_xlabel("1차 고유진동수 f₁ [Hz]")
ax2.set_title("수종별 f₁ 비교\n(생재, H=8m, DBH=20cm)")
ax2.legend(fontsize=8)
ax2.grid(True, axis="x", alpha=0.3)
ax2.tick_params(axis='y', labelsize=7)

# ── 5-c. f₁ 공진 조건 지도 (DBH × H) ────────────────────────────
ax3 = fig.add_subplot(gs[1, :2])
# 공진 위험: f₁이 0.2~0.8 Hz 범위에 가까울수록 위험
# 위험도 = exp(-((f1 - f_center)/(f_width))²)
f_center = 0.50
f_width  = 0.35
risk_map = np.exp(-((f1_map - f_center) / f_width)**2)

im = ax3.contourf(DBH_range*100, H_range, f1_map,
                  levels=20, cmap="RdYlGn_r")
# 공진 대역 등고선
cs = ax3.contour(DBH_range*100, H_range, f1_map,
                 levels=[0.20, 0.50, 0.80],
                 colors=["blue","red","blue"], linewidths=1.5)
ax3.clabel(cs, fmt="%.2f Hz", fontsize=8)
plt.colorbar(im, ax=ax3, label="f₁ [Hz]")

# 기준점 표시
ax3.scatter([20], [8.0], color="white", s=120, zorder=5,
            edgecolors="black", linewidths=1.5, label="기준점 (DBH=20cm, H=8m)")
ax3.set_xlabel("흉고직경 DBH [cm]")
ax3.set_ylabel("수고 H [m]")
ax3.set_title(f"1차 고유진동수 분포 지도 ({sp_ref.split('(')[0].strip()}, 생재)\n"
              f"파란 선: 0.2/0.8 Hz, 빨간 선: 0.5 Hz (공진 중심)")
ax3.legend(fontsize=8)
ax3.grid(True, alpha=0.2, color="white")

# ── 5-d. 수고에 따른 f₁ 변화 (수종별) ───────────────────────────
ax4 = fig.add_subplot(gs[1, 2])
H_vals = np.linspace(3.0, 15.0, 50)

for (sp, props), color in zip(STATIC_PROPS.items(), COLORS_SP):
    E_g   = props["E_L_GPa"] * 0.68
    rho_g = props["density_green"]
    f1s   = []
    for H in H_vals:
        D_b = D_BASE
        M_c = 150 * (H / 8.0)
        f1, *_ = natural_freq_with_crown(E_g, rho_g, D_b, H, M_c, TAPER)
        f1s.append(f1)
    short = sp.split("(")[0].strip()
    ax4.plot(H_vals, f1s, lw=1.8, color=color, label=short)

ax4.axhspan(0.20, 0.80, alpha=0.18, color="red", label="공진 위험 대역")
ax4.set_xlabel("수고 H [m]")
ax4.set_ylabel("1차 고유진동수 f₁ [Hz]")
ax4.set_title("수고 증가에 따른 f₁ 변화\n(DBH=20cm 고정)")
ax4.legend(fontsize=6.5, ncol=2)
ax4.grid(True, which="both", alpha=0.3)

plt.tight_layout()
fig_path = out_dir / "modal_verification.png"
plt.savefig(fig_path, dpi=150, bbox_inches="tight")
plt.show()
print(f"\n  그래프 저장: {fig_path}")

# ══════════════════════════════════════════════════════════════════════
# 6. 주요 발견 출력
# ══════════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("  핵심 발견: 공진 위험도 평가")
print("=" * 65)
for sp, f1 in species_f1.items():
    short = sp.split("(")[0].strip()
    in_band = 0.20 <= f1 <= 0.80
    risk = "★ 공진 위험" if in_band else "  안전 범위"
    print(f"  {short:<18}  f₁={f1:.3f} Hz  {risk}")

print(f"\n  [DBH=20cm 기준 공진 위험 수고]")
sp_ref_key = "은행나무 (Ginkgo biloba)"
E_g_ref  = STATIC_PROPS[sp_ref_key]["E_L_GPa"] * 0.68
rho_g_ref = STATIC_PROPS[sp_ref_key]["density_green"]
for target_f in [0.20, 0.50, 0.80]:
    # 이분법으로 f₁=target_f 되는 수고 탐색
    H_lo, H_hi = 3.0, 30.0
    for _ in range(50):
        H_mid = (H_lo + H_hi) / 2
        f1_mid, *_ = natural_freq_with_crown(
            E_g_ref, rho_g_ref, D_BASE, H_mid, 150*(H_mid/8), TAPER)
        if f1_mid > target_f:
            H_lo = H_mid
        else:
            H_hi = H_mid
    print(f"    f₁ = {target_f:.2f} Hz  →  H ≈ {H_mid:.1f} m")
