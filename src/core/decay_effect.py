"""
decay_effect.py
---------------
수목 부후(腐朽, Decay)가 근원부 단면 강성 및 피로수명에 미치는 영향 분석

부후 모델:
  - 중심부 공동(hollow) 모델: 반경 r_d 까지 부후로 탄성계수 0
  - 편심 공동 모델: 한쪽으로 치우친 부후 (응력집중 유발)
  - 잔존 강성 비율: η = I_residual / I_intact

실용 지표 — t/R 비율 (Shell thickness / Radius):
  t/R = (R - r_d) / R = 1 - r_d/R
  현장에서 저항추자(Resistograph)로 측정 가능
  임계값: t/R < 0.3 → 도복 위험 (Mattheck & Breloer 1994)

참고:
  Mattheck, C. & Breloer, H. (1994). The Body Language of Trees.
    HMSO, London.
  Kane, B. & Ryan, H.D.P. (2004). Wound closure and decay in trees.
    J. Arboric., 30(6), 394-400.
  Deflorio, G. et al. (2008). Dynamic properties of decayed wood.
    Tree Physiol., 28, 1043-1050.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import rainflow
import csv
from pathlib import Path
from wood_properties_db import STATIC_PROPS

matplotlib.rcParams['font.family'] = 'Malgun Gothic'
matplotlib.rcParams['axes.unicode_minus'] = False

out_dir = Path("output")
out_dir.mkdir(exist_ok=True)

# ══════════════════════════════════════════════════════════════════════
# 1. 단면 특성 계산 함수
# ══════════════════════════════════════════════════════════════════════

def section_intact(D_m):
    """건전 원형 단면: I, A, c"""
    R = D_m / 2
    I = np.pi * R**4 / 4
    A = np.pi * R**2
    return I, A, R

def section_hollow_concentric(D_m, decay_ratio):
    """
    동심원 공동(중심 부후) 단면
    decay_ratio = r_d / R  (부후 반경 / 전체 반경)
    I_hollow = π/4 × (R⁴ - r_d⁴)
    """
    R   = D_m / 2
    r_d = R * decay_ratio
    I   = np.pi * (R**4 - r_d**4) / 4
    A   = np.pi * (R**2 - r_d**2)
    c   = R          # 최외단까지 거리 변동 없음
    return I, A, c, r_d

def section_hollow_eccentric(D_m, decay_ratio, eccentricity=0.5):
    """
    편심 공동 단면 (한쪽 치우침)
    eccentricity: 공동 중심 이동량 / R (0=동심, 1=최대편심)
    → 응력집중계수 Kt = σ_max / σ_nominal 계산
    편심 공동의 단면 2차 모멘트는 평행축 정리로 계산
    """
    R   = D_m / 2
    r_d = R * decay_ratio
    e   = R * eccentricity * (1 - decay_ratio)  # 실제 편심 거리

    # 단면 2차 모멘트 (평행축 정리)
    I_full  = np.pi * R**4 / 4
    I_hole  = np.pi * r_d**4 / 4 + np.pi * r_d**2 * e**2
    I_resid = I_full - I_hole

    # 중립축 이동 (Ay = A1y1 - A2y2)
    A_full = np.pi * R**2
    A_hole = np.pi * r_d**2
    A_resid = A_full - A_hole
    y_na   = (A_full * 0 - A_hole * e) / A_resid  # 중립축 편이

    # 최대 응력 발생 위치: 부후 쪽 최외단
    c_max = R - y_na    # 부후 쪽 거리
    c_min = R + y_na    # 반대쪽 거리 (압축)

    # 응력집중계수 (편심 대비 동심)
    I_concentric = np.pi * (R**4 - r_d**4) / 4
    c_concentric = R
    Kt = (c_max / I_resid) / (c_concentric / I_concentric)

    return I_resid, A_resid, c_max, Kt


# ══════════════════════════════════════════════════════════════════════
# 2. 부후율별 단면 강성 비율 및 응력 증가율
# ══════════════════════════════════════════════════════════════════════

D_BASE  = 0.24    # [m] 근원부 직경 (DBH=20cm × 1.2)
I_int, A_int, R = section_intact(D_BASE)

decay_ratios = np.linspace(0, 0.90, 200)   # 부후 반경 / 전체 반경

# 동심 공동
I_hol   = np.array([section_hollow_concentric(D_BASE, dr)[0] for dr in decay_ratios])
eta_hol = I_hol / I_int      # 잔존 단면 강성 비율

# 응력 증가율 (F, L_arm, c 동일 → σ ∝ 1/I)
stress_amp_conc = I_int / I_hol   # 동심 공동

# 편심 공동 (eccentricity=0.5)
stress_amp_ecc  = np.array([
    section_hollow_eccentric(D_BASE, dr, eccentricity=0.5)[3]
    for dr in decay_ratios])
# 편심 공동의 절대 응력 증가율
I_ecc_arr = np.array([section_hollow_eccentric(D_BASE, dr, 0.5)[0] for dr in decay_ratios])
c_ecc_arr = np.array([section_hollow_eccentric(D_BASE, dr, 0.5)[2] for dr in decay_ratios])
stress_amp_ecc_abs = (c_ecc_arr / I_ecc_arr) / (R / I_int)

# ══════════════════════════════════════════════════════════════════════
# 3. Mattheck 임계 t/R 기준
# ══════════════════════════════════════════════════════════════════════
# t/R = (R - r_d) / R = 1 - decay_ratio
# 임계값: t/R = 0.3  → decay_ratio = 0.7
MATTHECK_CRITICAL = 0.70   # decay_ratio 기준

print("=" * 60)
print("  부후율별 잔존 단면 강성 및 응력 증가율")
print("=" * 60)
print(f"  {'부후율 r_d/R':>12}  {'t/R':>6}  {'η (강성비)':>10}  "
      f"{'응력배율(동심)':>12}  {'응력배율(편심)':>12}")
for dr in [0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
    idx = np.argmin(np.abs(decay_ratios - dr))
    tR  = 1 - dr
    eta = eta_hol[idx]
    sa_c = stress_amp_conc[idx]
    sa_e = stress_amp_ecc_abs[idx]
    flag = " ← Mattheck 임계" if abs(dr - 0.7) < 0.01 else ""
    print(f"  {dr:>12.1f}  {tR:>6.2f}  {eta:>10.3f}  "
          f"{sa_c:>12.2f}×  {sa_e:>12.2f}×{flag}")

# ══════════════════════════════════════════════════════════════════════
# 4. 부후율별 피로수명 감소 (Kaimal 기반 응력 적용)
# ══════════════════════════════════════════════════════════════════════
# 기존 wind_time_history.py 의 항력 시계열을 로드하여 사용

def load_drag_force(csv_path):
    t_list, F_list = [], []
    with open(csv_path, "r", encoding="utf-8") as f:
        for row in csv.reader(f):
            if not row or row[0].startswith("#") or row[0] == "time_s":
                continue
            t_list.append(float(row[0]))
            F_list.append(float(row[3]))
    return np.array(t_list), np.array(F_list)

def fatigue_life_decay(F_t, I_sec, c_sec, L_arm, MOR_MPa, a_sn, b_sn,
                       T_sim=600.0, T_year=100*3600):
    sigma_t = F_t * L_arm * c_sec / I_sec * 1e-6
    D = 0.0
    for rng, *_ in rainflow.extract_cycles(sigma_t):
        amp = rng / 2.0
        if amp <= 0:
            continue
        ratio = amp / MOR_MPa
        Nf = 1.0 if ratio >= a_sn else 10.0**((a_sn - ratio) / b_sn)
        D += 1.0 / Nf
    D_annual = D * (T_year / T_sim)
    return 1.0 / D_annual if D_annual > 0 else np.inf

# 항력 시계열 로드
csv_kai = out_dir / "wind_kaimal.csv"
if csv_kai.exists():
    _, F_arr = load_drag_force(csv_kai)
    L_ARM = 8.0 * (2/3)

    # S-N 파라미터 (은행나무, 생재 MOR 적용)
    sp_key  = "은행나무 (Ginkgo biloba)"
    props   = STATIC_PROPS[sp_key]
    MOR_g   = props["MOR_MPa"] * 0.60
    A_SN, B_SN = 0.97, 0.072

    print("\n  [부후율별 피로수명] 은행나무, 생재, Kaimal, 연 100h 강풍")
    print(f"  {'부후율':>8}  {'t/R':>5}  {'수명(동심)':>12}  {'수명(편심)':>12}  {'수명감소율(동심)':>14}")

    decay_study = [0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
    L_baseline  = None
    decay_lives_conc = []
    decay_lives_ecc  = []

    for dr in decay_study:
        I_c, _, c_c, _  = section_hollow_concentric(D_BASE, dr)
        I_e, _, c_e, _  = section_hollow_eccentric(D_BASE, dr, 0.5)

        L_conc = fatigue_life_decay(F_arr, I_c, c_c, L_ARM, MOR_g, A_SN, B_SN)
        L_ecc  = fatigue_life_decay(F_arr, I_e, c_e, L_ARM, MOR_g, A_SN, B_SN)

        if L_baseline is None:
            L_baseline = L_conc

        decay_lives_conc.append(min(L_conc, 1e8))
        decay_lives_ecc.append(min(L_ecc, 1e8))

        tR   = 1 - dr
        red  = (L_baseline - L_conc) / L_baseline * 100 if L_baseline < 1e7 else 0
        flag = " ★" if dr >= 0.7 else ""
        print(f"  {dr:>8.1f}  {tR:>5.2f}  {min(L_conc,1e8):>12.0f}  "
              f"{min(L_ecc,1e8):>12.0f}  {red:>12.1f}%{flag}")
else:
    decay_study = []
    decay_lives_conc = []
    decay_lives_ecc  = []
    F_arr = None
    print("  [경고] wind_kaimal.csv 없음 — wind_time_history.py 먼저 실행하세요")

# ══════════════════════════════════════════════════════════════════════
# 5. 시각화
# ══════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(2, 3, figsize=(16, 10))
fig.suptitle("부후(腐朽)가 가로수 구조 안전성에 미치는 영향\n"
             f"은행나무, H=8m, DBH=20cm, 생재(MOR={MOR_g:.0f}MPa)", fontsize=11)

# ── 5-a. 잔존 단면 강성 비율 ─────────────────────────────────────
ax = axes[0, 0]
ax.plot(decay_ratios, eta_hol, lw=2, color="steelblue", label="동심 공동 η = I/I₀")
I_full_arr = np.ones_like(decay_ratios)
ax.fill_between(decay_ratios, eta_hol, alpha=0.15, color="steelblue")
ax.axvline(MATTHECK_CRITICAL, color="red", lw=1.5, ls="--",
           label=f"Mattheck 임계 (r_d/R={MATTHECK_CRITICAL})")
ax.axhspan(0, 0.3, alpha=0.10, color="red")
ax.set_xlabel("부후 반경 비율 r_d/R")
ax.set_ylabel("잔존 강성 비율 η = I_부후 / I_건전")
ax.set_title("부후율에 따른 단면 강성 감소")
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
# 2차 x축 (t/R)
ax2 = ax.twiny()
ax2.set_xlim(ax.get_xlim())
ax2.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax2.set_xticklabels([f"{1-v:.1f}" for v in [0, 0.2, 0.4, 0.6, 0.8, 1.0]], fontsize=7)
ax2.set_xlabel("잔존 벽두께 비율 t/R", fontsize=8)

# ── 5-b. 응력 증가율 비교 ────────────────────────────────────────
ax = axes[0, 1]
ax.plot(decay_ratios, stress_amp_conc,     lw=2, color="steelblue",
        label="동심 공동")
ax.plot(decay_ratios, stress_amp_ecc_abs,  lw=2, color="crimson",
        label="편심 공동 (e=0.5R)")
ax.axvline(MATTHECK_CRITICAL, color="red", lw=1.5, ls="--",
           label="Mattheck 임계")
ax.axhline(1.0, color="gray", lw=0.8, ls=":")
ax.set_xlabel("부후 반경 비율 r_d/R")
ax.set_ylabel("응력 증가 배율 (건전 대비)")
ax.set_title("부후에 의한 근원부 응력 증폭")
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_ylim([0, 15])

# ── 5-c. 단면 형상 시각화 (3가지 부후 수준) ─────────────────────
ax = axes[0, 2]
ax.set_aspect("equal")
theta = np.linspace(0, 2*np.pi, 300)
R_plot = D_BASE / 2

# 건전 단면
ax.plot(R_plot*np.cos(theta), R_plot*np.sin(theta),
        "k-", lw=2, label="건전 단면")
ax.fill(R_plot*np.cos(theta), R_plot*np.sin(theta),
        color="burlywood", alpha=0.3)

# 부후 단면들 (동심)
for dr, color, label in [(0.4, "darkorange", "부후 40%"),
                          (0.7, "crimson",   "부후 70% (임계)")]:
    r_d = R_plot * dr
    ax.fill(r_d*np.cos(theta), r_d*np.sin(theta),
            color=color, alpha=0.5)
    ax.plot(r_d*np.cos(theta), r_d*np.sin(theta),
            ls="--", color=color, lw=1.5, label=label)

# Mattheck t/R = 0.3 최소 잔존 벽두께 화살표
t_min = R_plot * 0.3
ax.annotate("", xy=(R_plot*0.7, 0), xytext=(R_plot, 0),
            arrowprops=dict(arrowstyle="<->", color="red", lw=1.5))
ax.text(R_plot*0.83, 0.01, "t_min", color="red", fontsize=8, ha="center")

ax.set_xlim([-R_plot*1.3, R_plot*1.3])
ax.set_ylim([-R_plot*1.3, R_plot*1.3])
ax.set_title("부후 단면 형상 (동심 공동)")
ax.legend(fontsize=7, loc="upper right")
ax.grid(True, alpha=0.2)
ax.set_xlabel("x [m]")
ax.set_ylabel("y [m]")

# ── 5-d. 피로수명 vs 부후율 ──────────────────────────────────────
ax = axes[1, 0]
if decay_study:
    ax.semilogy(decay_study, decay_lives_conc, "o-", lw=2,
                color="steelblue", label="동심 공동", markersize=6)
    ax.semilogy(decay_study, decay_lives_ecc, "s--", lw=2,
                color="crimson", label="편심 공동", markersize=6)
    ax.axvline(MATTHECK_CRITICAL, color="red", lw=1.5, ls="--",
               label="Mattheck 임계")
    ax.set_xlabel("부후 반경 비율 r_d/R")
    ax.set_ylabel("피로수명 [년]")
    ax.set_title("부후율별 피로수명 감소\n(Kaimal, 연 100h 강풍, 은행나무 생재)")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)

# ── 5-e. 수종별 부후 임계 부후율 ─────────────────────────────────
ax = axes[1, 1]
# 파단 임계 부후율: σ_peak(U_10yr, 부후 단면) ≥ MOR
U_10YR = 12.3    # [m/s] 서울 10년 재현 풍속
RHO_AIR, CD = 1.225, 0.5
A_CROWN_REF = (8.0 * 0.4)**2 * np.pi / 4
F_10yr = 0.5 * RHO_AIR * CD * A_CROWN_REF * U_10YR**2
L_ARM_REF = 8.0 * (2/3)

COLORS_SP = ["steelblue","hotpink","darkorange","forestgreen","purple","saddlebrown"]
for (sp, props), color in zip(STATIC_PROPS.items(), COLORS_SP):
    MOR_g = props["MOR_MPa"] * 0.60
    # σ = F × L × c / I ≥ MOR → c/I ≥ MOR/(F×L)
    target = MOR_g * 1e6 / (F_10yr * L_ARM_REF)
    # 부후율 탐색
    crit_dr = None
    for dr_test in np.linspace(0, 0.95, 1000):
        I_c, _, c_c, _ = section_hollow_concentric(D_BASE, dr_test)
        if c_c / I_c >= target:
            crit_dr = dr_test
            break
    short = sp.split("(")[0].strip()
    if crit_dr is not None:
        ax.barh(short, crit_dr, color=color, alpha=0.8)
        ax.text(crit_dr + 0.01, short, f"{crit_dr:.2f}", va="center", fontsize=8)
    else:
        ax.barh(short, 0.95, color=color, alpha=0.4)
        ax.text(0.96, short, ">0.95", va="center", fontsize=8)

ax.axvline(MATTHECK_CRITICAL, color="red", lw=1.5, ls="--",
           label=f"Mattheck 임계 ({MATTHECK_CRITICAL})")
ax.set_xlabel("파단 임계 부후율 r_d/R")
ax.set_title(f"수종별 정적 파단 임계 부후율\n(U_10yr={U_10YR}m/s, 서울 기준)")
ax.legend(fontsize=8)
ax.grid(True, axis="x", alpha=0.3)
ax.set_xlim([0, 1.05])
ax.tick_params(axis='y', labelsize=7)

# ── 5-f. 종합 안전 평가 매트릭스 ─────────────────────────────────
ax = axes[1, 2]
# DBH별 Mattheck 임계 판정 + 피로 위험
DBH_range = np.array([0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40])
decay_levels = np.array([0.0, 0.3, 0.5, 0.7, 0.9])

# 위험도 점수: 응력 증가 × 피로 민감도 (단순 지수화)
risk_matrix = np.zeros((len(decay_levels), len(DBH_range)))
for i, dr in enumerate(decay_levels):
    for j, dbh in enumerate(DBH_range):
        D_b = dbh * 1.2
        I_c_, _, c_c_, _ = section_hollow_concentric(D_b, dr)
        I_i_, _, R_i_ = section_intact(D_b)
        stress_mult = (c_c_ / I_c_) / (R_i_ / I_i_)
        risk_matrix[i, j] = min(stress_mult, 10)

im = ax.imshow(risk_matrix, cmap="RdYlGn_r", aspect="auto",
               vmin=1, vmax=8, origin="lower")
ax.set_xticks(range(len(DBH_range)))
ax.set_xticklabels([f"{d*100:.0f}" for d in DBH_range])
ax.set_yticks(range(len(decay_levels)))
ax.set_yticklabels([f"{dr:.1f}" for dr in decay_levels])
ax.set_xlabel("흉고직경 DBH [cm]")
ax.set_ylabel("부후율 r_d/R")
ax.set_title("응력 증폭 배율 매트릭스\n(붉을수록 위험)")
plt.colorbar(im, ax=ax, label="응력 증폭 배율")
for i in range(len(decay_levels)):
    for j in range(len(DBH_range)):
        ax.text(j, i, f"{risk_matrix[i,j]:.1f}×",
                ha="center", va="center", fontsize=7,
                color="white" if risk_matrix[i,j] > 5 else "black")

plt.tight_layout()
fig_path = out_dir / "decay_effect.png"
plt.savefig(fig_path, dpi=150, bbox_inches="tight")
plt.show()
print(f"\n  그래프 저장: {fig_path}")
