"""
winkler_analysis.py
-------------------
Winkler 지반 스프링이 가로수 피로수명에 미치는 영향 분석
(ABAQUS 재해석 없이 Python으로 보정계수 적용)

원리:
  완전 고정(Clamped) 경계조건 → 과도하게 보수적
  실제 뿌리-토양 시스템은 유한 강성 → 근원부 회전 허용

Winkler 모델:
  토양 = 분산 스프링 k_s [kN/m/m]
  뿌리 회전 스프링 K_r = k_s * L_r^3 / 3  [kN.m/rad]
    L_r : 뿌리 영향 깊이 (= 0.5~1.5 * DBH)
    k_s : 토양 반력계수 (지반 종류별)

보정 방법:
  완전 고정 대비 유연 경계조건에서:
    1. 응력 보정: σ_winkler = σ_fixed * C_sigma(K_r)
    2. 고유진동수 보정: f1_winkler = f1_fixed * C_f(K_r)
  C_sigma, C_f 는 Euler-Bernoulli 보 이론으로 해석적 계산

참고:
  Spatz et al. (2007). Biomechanics of plant stems.
  Coutts (1983). Root architecture and tree stability.
  Rahardjo et al. (2014). Root anchorage of urban trees.
"""

import sys, os
PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
if PROJECT_DIR not in sys.path:
    sys.path.insert(0, PROJECT_DIR)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import rainflow
from pathlib import Path
from scipy.special import gamma as gamma_func
from wood_properties_db import STATIC_PROPS, SN_PARAMS, MC_CORRECTION

matplotlib.rcParams['font.family'] = 'Malgun Gothic'
matplotlib.rcParams['axes.unicode_minus'] = False

out_dir = Path("output")

# ======================================================================
# 1. 토양 종류별 Winkler 반력계수
# ======================================================================
# k_s [kN/m^3] : 수평 지반반력계수 (Broms 1964, IS 2911)
SOIL_TYPES = {
    "연약 점토 (도시 매립지)": {"k_s": 2000,   "color": "#e74c3c"},
    "중간 점토 (일반 도심)":   {"k_s": 8000,   "color": "#e67e22"},
    "단단한 점토/느슨한 모래": {"k_s": 25000,  "color": "#f1c40f"},
    "조밀한 모래/자갈":        {"k_s": 80000,  "color": "#2ecc71"},
    "완전 고정 (Clamped)":     {"k_s": np.inf, "color": "#2c3e50"},
}

# ======================================================================
# 2. 수목 파라미터 (H8m, DBH20cm 기준)
# ======================================================================
H      = 8.0          # [m]  수고
DBH    = 0.20         # [m]  흉고직경
D_BASE = DBH * 1.2    # [m]  근원부 직경
I_BASE = np.pi * D_BASE**4 / 64.0   # [m^4]
E_WOOD = 7.0e9 * 0.68                # [Pa]  생재 탄성계수
EI     = E_WOOD * I_BASE             # [N.m^2]

# 뿌리 영향 깊이
# 도시 가로수 뿌리 정착 깊이: 0.5~1.5 m (Coutts 1983, Rahardjo 2014)
# 기본값 1.0 m (중간값)
L_r = 1.0   # [m]

# ======================================================================
# 3. Winkler 회전 스프링 강성 계산
# ======================================================================
# 분산 스프링 모델: K_r = k_s * D_BASE * L_r^3 / 3
# (원통형 근원부, 수평 하중에 대한 등가 회전 강성)
def root_spring(k_s, D_base, L_r):
    """뿌리 회전 스프링 강성 K_r [N.m/rad]"""
    if np.isinf(k_s):
        return np.inf
    return k_s * 1000 * D_base * L_r**3 / 3.0   # k_s: kN->N

# ======================================================================
# 4. 응력 보정계수 및 고유진동수 보정계수
# ======================================================================
def stress_correction(K_r, EI, H):
    """
    Winkler 스프링 경계조건에서 응력 보정계수 C_sigma
    완전 고정 대비 근원부 모멘트 비율

    외팔보 + 근원부 회전 스프링:
      회전각 theta = M_base / K_r
      추가 변위  delta_extra = theta * H
      실제 모멘트 M_base = F*H / (1 + 3EI/(K_r*H))
      -> C_sigma = M_base(Winkler) / M_base(fixed) = 1 / (1 + 3EI/(K_r*H))

    (완전 고정: C_sigma = 1.0)
    """
    if np.isinf(K_r):
        return 1.0
    denom = 1.0 + 3.0 * EI / (K_r * H)
    return 1.0 / denom

def freq_correction(K_r, EI, H, M_total):
    """
    Winkler 스프링 경계조건에서 고유진동수 보정계수 C_f
    f1_winkler = f1_fixed * C_f

    등가 강성:
      k_eq = 3EI/H^3 + K_r/H^2  (수평 스프링 등가)
      or 회전 포함:  1/k_eff = H^3/(3EI) + H^2/K_r
      -> k_eff = 3EI*K_r / (3EI + K_r*H) * (1/H^2)
      f1 ~ sqrt(k_eff / M_total)
      C_f = f1_winkler / f1_fixed = sqrt(k_eff_winkler / k_eff_fixed)
    """
    if np.isinf(K_r):
        return 1.0
    k_fixed   = 3.0 * EI / H**3
    # 유연 경계: 등가 수평 강성
    k_winkler = 1.0 / (H**3/(3.0*EI) + H**2/K_r)
    return np.sqrt(k_winkler / k_fixed)

# ======================================================================
# 5. 응력 시계열 로드 및 Rainflow (H8m-DBH20cm 기준)
# ======================================================================
def load_stress_csv(csv_path):
    sigma = []
    with open(csv_path, 'r', encoding='utf-8-sig') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#') or s.startswith('time'):
                continue
            parts = s.split(',')
            try:
                sigma.append(abs(float(parts[-1])))
            except (ValueError, IndexError):
                continue
    return np.array(sigma)

def weibull_pdf(U, k, c):
    if U <= 0: return 0.0
    return (k/c) * (U/c)**(k-1) * np.exp(-(U/c)**k)

def stress_factor_decay(dr):
    return 1.0 / (1.0 - dr**4) if dr < 1.0 else np.inf

csv_ref = out_dir / "tree_H8m_DBH20cm_stress.csv"
sigma_abq = load_stress_csv(csv_ref)
cycles    = list(rainflow.extract_cycles(sigma_abq))
ranges_rf = np.array([x[0] for x in cycles])
counts_rf = np.array([x[2] for x in cycles])
sigma_max_fixed = sigma_abq.max()

U_REF  = 20.0
T_SIM  = 600.0
T_YEAR = 365.25 * 24 * 3600
U_BINS = np.arange(1.0, 35.0, 1.0)

SPECIES   = "은행나무 (Ginkgo biloba)"
MOR_GRN   = STATIC_PROPS[SPECIES]["MOR_MPa"] * MC_CORRECTION["MOR"]
A_SN, B_SN = SN_PARAMS[SPECIES]["a_R0"], SN_PARAMS[SPECIES]["b_R0"]

# Weibull 서울 도심
K_W, C_W = 1.6, 4.0

# 수관 질량 (고유진동수 계산용)
M_CROWN = 150.0   # [kg]
RHO_W   = 780.0
M_TRUNK = RHO_W * np.pi * (D_BASE/2)**2 * H * 0.33   # 등가 집중 질량(1/3 규칙)
M_TOTAL = M_CROWN + M_TRUNK

# 완전 고정 기준 f1 (ABAQUS 결과)
F1_FIXED = 0.4271   # [Hz] (summary.csv 에서)

# ======================================================================
# 6. 토양별 보정계수 계산
# ======================================================================
def annual_damage_winkler(ranges, counts, sf_decay, C_sigma,
                          MOR, a, b, k_w, c_w):
    """응력 보정계수 C_sigma 포함 연간 손상"""
    D_ann = 0.0
    for U in U_BINS:
        ws   = (U / U_REF)**2
        prob = weibull_pdf(U, k_w, c_w) * 1.0
        T_U  = prob * T_YEAR
        if T_U < 1.0: continue
        n_sim = T_U / T_SIM
        D = 0.0
        for rng, cnt in zip(ranges * ws, counts):
            amp = (rng * sf_decay * C_sigma) / 2.0
            if amp <= 0: continue
            ratio = amp / MOR
            Nf = 1.0 if ratio >= a else 10.0**((a-ratio)/b)
            D += cnt / Nf
        D_ann += D * n_sim
    return D_ann

print("=" * 65)
print("  Winkler 지반 스프링 민감도 분석")
print(f"  기준: H={H}m, DBH={DBH*100:.0f}cm, {SPECIES}")
print(f"  완전 고정 f1 = {F1_FIXED:.4f} Hz")
print(f"  완전 고정 sigma_max = {sigma_max_fixed:.2f} MPa")
print("=" * 65)

print(f"\n  {'토양 조건':22}  {'K_r [kN.m/rad]':>16}  "
      f"{'C_sigma':>8}  {'C_f':>6}  "
      f"{'f1 [Hz]':>8}  {'sigma_max [MPa]':>16}")
print("  " + "-" * 80)

soil_results = {}
for soil_name, sp in SOIL_TYPES.items():
    K_r    = root_spring(sp["k_s"], D_BASE, L_r)
    C_sig  = stress_correction(K_r, EI, H)
    C_f    = freq_correction(K_r, EI, H, M_TOTAL)
    f1_w   = F1_FIXED * C_f
    sig_w  = sigma_max_fixed * C_sig

    if np.isinf(K_r):
        Kr_str = "inf"
    else:
        Kr_str = f"{K_r/1000:.1f}"

    print(f"  {soil_name:22}  {Kr_str:>16}  "
          f"{C_sig:>8.4f}  {C_f:>6.4f}  "
          f"{f1_w:>8.4f}  {sig_w:>14.2f} MPa")

    soil_results[soil_name] = {
        "K_r": K_r, "C_sig": C_sig, "C_f": C_f,
        "f1": f1_w, "sigma_max": sig_w,
        "color": sp["color"]
    }

# ======================================================================
# 7. 피로수명 계산 (토양별 x 부후율)
# ======================================================================
DECAY_RATIOS = np.array([0.0, 0.3, 0.5, 0.6, 0.7, 0.75, 0.8])

print(f"\n  피로수명 [년] - 서울 도심, 평균 S-N")
print(f"  {'토양':22}", end="")
for dr in DECAY_RATIOS:
    print(f"  {f'dr={dr:.2f}':>10}", end="")
print()
print("  " + "-" * 78)

soil_lives = {}
for soil_name, sr in soil_results.items():
    soil_lives[soil_name] = []
    row = f"  {soil_name:22}"
    for dr in DECAY_RATIOS:
        sf_d  = stress_factor_decay(dr)
        D_ann = annual_damage_winkler(
            ranges_rf, counts_rf, sf_d, sr["C_sig"],
            MOR_GRN, A_SN, B_SN, K_W, C_W)
        L = 1.0 / D_ann if D_ann > 0 else np.inf
        soil_lives[soil_name].append(L)
        if L >= 1e7: s = ">1000만"
        elif L >= 1e4: s = f"{L:.0f}"
        elif L >= 1: s = f"{L:.1f}"
        else: s = f"{L:.2f}"
        row += f"  {s:>10}"
    print(row)

# ======================================================================
# 8. 시각화
# ======================================================================
fig, axes = plt.subplots(2, 3, figsize=(17, 11))
fig.suptitle(
    "Winkler 지반 스프링 민감도 분석\n"
    f"H={H}m, DBH={DBH*100:.0f}cm, {SPECIES}, 서울 도심",
    fontsize=11
)

# ── 8-a. 토양별 보정계수 ─────────────────────────────────────────
ax = axes[0, 0]
soil_names = list(soil_results.keys())
C_sigs = [sr["C_sig"] for sr in soil_results.values()]
C_fs   = [sr["C_f"]   for sr in soil_results.values()]
colors_soil = [sr["color"] for sr in soil_results.values()]
x = np.arange(len(soil_names))
b1 = ax.bar(x - 0.2, C_sigs, 0.35, label='응력 보정 C_sigma',
            color=colors_soil, alpha=0.9)
b2 = ax.bar(x + 0.2, C_fs,   0.35, label='진동수 보정 C_f',
            color=colors_soil, alpha=0.4)
ax.axhline(1.0, color='black', lw=1, ls='--')
for bar, v in zip(b1, C_sigs):
    ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.01,
            f'{v:.3f}', ha='center', va='bottom', fontsize=7)
ax.set_xticks(x)
ax.set_xticklabels([s.split('(')[0].strip() for s in soil_names],
                    rotation=25, fontsize=7, ha='right')
ax.set_ylabel('보정계수')
ax.set_title('토양별 응력/진동수 보정계수\n(1.0=완전 고정, <1.0=유연 경계)')
ax.legend(fontsize=8)
ax.grid(True, axis='y', alpha=0.3)

# ── 8-b. 토양별 f1 및 sigma_max ──────────────────────────────────
ax = axes[0, 1]
f1_vals  = [sr["f1"] for sr in soil_results.values()]
sig_vals = [sr["sigma_max"] for sr in soil_results.values()]
ax2 = ax.twinx()
bars = ax.bar(x, f1_vals, 0.4, color=colors_soil, alpha=0.7, label='f1 [Hz]')
ax2.plot(x, sig_vals, 'ko-', lw=2, markersize=8, label='sigma_max [MPa]')
ax2.axhline(MOR_GRN, color='red', lw=1.5, ls='--',
            label=f'MOR={MOR_GRN:.0f}MPa')
ax.set_xticks(x)
ax.set_xticklabels([s.split('(')[0].strip() for s in soil_names],
                    rotation=25, fontsize=7, ha='right')
ax.set_ylabel('f1 [Hz]', color='steelblue')
ax2.set_ylabel('sigma_max [MPa]', color='black')
ax.set_title('토양 조건별 f1 및 최대 응력')
lines1, labels1 = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax.legend(lines1+lines2, labels1+labels2, fontsize=7)
ax.grid(True, axis='y', alpha=0.3)

# ── 8-c. 토양별 부후율-수명 곡선 ─────────────────────────────────
ax = axes[0, 2]
for soil_name, sr in soil_results.items():
    lives = [min(L, 1e8) for L in soil_lives[soil_name]]
    lw = 2.5 if 'Clamped' in soil_name else 1.5
    ax.semilogy(DECAY_RATIOS, lives, 'o-', lw=lw,
                color=sr["color"], markersize=5,
                label=soil_name.split('(')[0].strip())
ax.axvline(0.7, color='red', lw=1.5, ls=':', label='Mattheck 임계')
ax.axhline(50,  color='gray', lw=1, ls='--', label='목표수명 50년')
ax.set_xlabel('부후율 r_d/R')
ax.set_ylabel('피로수명 [년]')
ax.set_title('토양 조건별 부후율-수명 곡선\n(서울 도심, 평균 S-N)')
ax.legend(fontsize=7)
ax.grid(True, which='both', alpha=0.3)

# ── 8-d. 뿌리 깊이 민감도 ────────────────────────────────────────
ax = axes[1, 0]
L_r_range = np.linspace(0.3, 2.0, 80)   # L_r [m]
soil_sel = {
    "연약 점토": 2000,
    "중간 점토": 8000,
    "조밀한 모래": 80000,
}
for soil_n, k_s_val in soil_sel.items():
    C_sigs_lr = []
    for lr_i in L_r_range:
        K_r_i = root_spring(k_s_val, D_BASE, lr_i)
        C_sigs_lr.append(stress_correction(K_r_i, EI, H))
    color = [sr["color"] for sn, sr in soil_results.items()
             if soil_n in sn][0]
    ax.plot(L_r_range, C_sigs_lr, lw=2, color=color, label=soil_n)
ax.axhline(1.0, color='black', lw=1, ls='--', label='완전 고정')
ax.axvline(L_r, color='gray', lw=1, ls=':', label=f'L_r={L_r}m (기준)')
ax.set_xlabel('뿌리 정착 깊이 L_r [m]')
ax.set_ylabel('응력 보정계수 C_sigma')
ax.set_title('뿌리 깊이에 따른 응력 보정계수')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_ylim([0, 1.1])

# ── 8-e. 건전목 수명 vs k_s (연속) ──────────────────────────────
ax = axes[1, 1]
k_s_range = np.logspace(2.5, 6, 100)   # 316 ~ 1,000,000 kN/m^3
L_healthy = []
for k_s_val in k_s_range:
    K_r_i  = root_spring(k_s_val, D_BASE, L_r)
    C_sig  = stress_correction(K_r_i, EI, H)
    D_ann  = annual_damage_winkler(
        ranges_rf, counts_rf, 1.0, C_sig,
        MOR_GRN, A_SN, B_SN, K_W, C_W)
    L = 1.0 / D_ann if D_ann > 0 else 1e9
    L_healthy.append(min(L, 1e9))

ax.semilogx(k_s_range, L_healthy, lw=2.5, color='steelblue')
# 토양 종류 마킹
for soil_name, sp in SOIL_TYPES.items():
    if np.isinf(sp["k_s"]): continue
    K_r_i = root_spring(sp["k_s"], D_BASE, L_r)
    C_sig = stress_correction(K_r_i, EI, H)
    D_ann = annual_damage_winkler(
        ranges_rf, counts_rf, 1.0, C_sig,
        MOR_GRN, A_SN, B_SN, K_W, C_W)
    L = min(1.0/D_ann if D_ann > 0 else 1e9, 1e9)
    ax.plot(sp["k_s"], L, 'o', markersize=10,
            color=sp["color"],
            label=soil_name.split('(')[0].strip())
ax.axhline(soil_lives["완전 고정 (Clamped)"][0],
           color='black', lw=1.5, ls='--', label='완전 고정')
ax.set_xlabel('지반반력계수 k_s [kN/m^3]')
ax.set_ylabel('피로수명 [년]')
ax.set_title('지반 강성 vs 건전목 피로수명\n(서울 도심, 평균 S-N)')
ax.legend(fontsize=7)
ax.grid(True, which='both', alpha=0.3)

# ── 8-f. 종합 안전 평가 (토양 x 부후율) ──────────────────────────
ax = axes[1, 2]
soil_labels_short = [s.split('(')[0].strip() for s in soil_names]
life_mat = np.zeros((len(SOIL_TYPES), len(DECAY_RATIOS)))
for i, (soil_name, _) in enumerate(soil_results.items()):
    for j, L in enumerate(soil_lives[soil_name]):
        life_mat[i, j] = np.log10(max(min(L, 1e8), 0.1))

im = ax.imshow(life_mat, cmap='RdYlGn', aspect='auto',
               vmin=0, vmax=8, origin='upper')
ax.set_xticks(range(len(DECAY_RATIOS)))
ax.set_xticklabels([f'{dr:.2f}' for dr in DECAY_RATIOS], fontsize=8)
ax.set_yticks(range(len(soil_labels_short)))
ax.set_yticklabels(soil_labels_short, fontsize=8)
ax.set_xlabel('부후율 r_d/R')
ax.set_title('피로수명 열지도 log10(년)\n(토양 x 부후율)')
cbar = plt.colorbar(im, ax=ax)
cbar.set_label('log10(년)')
cbar.set_ticks([0, 2, 4, 6, 8])
cbar.set_ticklabels(['<1', '100', '10k', '1M', '>1억'])
for i in range(len(SOIL_TYPES)):
    for j in range(len(DECAY_RATIOS)):
        L = 10**life_mat[i, j]
        if L < 1:     txt = f'{L:.1f}'
        elif L < 100: txt = f'{L:.0f}'
        elif L < 1e4: txt = f'{L/1000:.0f}k'
        else:         txt = f'{L/1e6:.0f}M'
        ax.text(j, i, txt, ha='center', va='center', fontsize=7,
                color='white' if life_mat[i, j] < 3 else 'black')

plt.tight_layout()
fig_path = out_dir / 'winkler_analysis.png'
plt.savefig(fig_path, dpi=150, bbox_inches='tight')
plt.show()
print(f'\n  그래프 저장: {fig_path}')

# ======================================================================
# 9. 결론 요약
# ======================================================================
print("\n" + "=" * 65)
print("  Winkler 분석 결론 요약")
print("=" * 65)
print(f"  기준 케이스: H={H}m, DBH={DBH*100:.0f}cm, 은행나무 생재")
print(f"  완전 고정 sigma_max = {sigma_max_fixed:.2f} MPa")
print()
for soil_name, sr in soil_results.items():
    L0 = soil_lives[soil_name][0]
    L7_idx = np.argmin(np.abs(DECAY_RATIOS - 0.7))
    L7 = soil_lives[soil_name][L7_idx]
    def fmt(L):
        if L >= 1e7: return ">1000만년"
        if L >= 1e4: return f"{L:.0f}년"
        return f"{L:.1f}년"
    print(f"  {soil_name}")
    print(f"    C_sigma={sr['C_sig']:.4f}  f1={sr['f1']:.4f}Hz  "
          f"sigma_max={sr['sigma_max']:.1f}MPa")
    print(f"    건전목={fmt(L0)}  /  부후0.7={fmt(L7)}")
    print()
print("  핵심 발견:")
print("  1. 연약 지반일수록 응력 감소 (경계조건 유연화)")
print("  2. 완전 고정 가정은 보수적 (안전측) 결과")
print("  3. 부후 영향이 지반 조건보다 수명에 더 크게 기여")
print("  4. 도심 식재지반 (중간 점토) 기준 완전 고정 대비")
C_ref = soil_results["중간 점토 (일반 도심)"]["C_sig"]
print(f"     응력 {(1-C_ref)*100:.1f}% 감소 -> 수명 보정 가능")
print("=" * 65)
