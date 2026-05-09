"""
decay_coupled_fatigue.py
------------------------
ABAQUS 동적해석 결과 + 부후(腐朽) + S-N 하한값 + Weibull 풍속분포
→ 현실적인 피로수명 매트릭스 도출

핵심 수치 (ABAQUS 해석):
  - σ_max = 20.734 MPa  (U=20 m/s, 10분 시뮬레이션)
  - DAF   = 1.130       (동적 증폭 계수, Python 단독 대비)

접근 방법:
  1. ABAQUS 응력 시계열을 부후율별로 증폭 (σ × stress_factor(dr))
  2. S-N 평균값 vs 하한값(95% 신뢰구간) 두 가지 적용
  3. Weibull 풍속 분포(서울/인천/제주)로 연간 손상 가중합산
  4. 부후율(0.0~0.9) × 지역 × S-N 신뢰수준 → 피로수명 매트릭스
"""

import sys
import os

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
if PROJECT_DIR not in sys.path:
    sys.path.insert(0, PROJECT_DIR)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import rainflow
import csv
from pathlib import Path
from scipy.special import gamma as gamma_func
from wood_properties_db import STATIC_PROPS, SN_PARAMS, MC_CORRECTION

matplotlib.rcParams['font.family'] = 'Malgun Gothic'
matplotlib.rcParams['axes.unicode_minus'] = False

out_dir = Path("output")
out_dir.mkdir(exist_ok=True)

# ══════════════════════════════════════════════════════════════════════
# 0. 수종 및 S-N 설정
# ══════════════════════════════════════════════════════════════════════
SPECIES  = "은행나무 (Ginkgo biloba)"
MOR_DRY  = STATIC_PROPS[SPECIES]["MOR_MPa"]       # 68 MPa (소건재)
MOR_GRN  = MOR_DRY * MC_CORRECTION["MOR"]          # 40.8 MPa (생재)

SN = SN_PARAMS[SPECIES]
# S-N 평균값 (현재 DB)
A_MEAN = SN["a_R0"]   # 0.97
B_MEAN = SN["b_R0"]   # 0.072

# S-N 하한값 — 95% 신뢰구간 (평균 - 1.645σ)
# 목재 S-N 산포: σ_a ≈ 5~10% of MOR (Tsai & Ansell 1990)
# → a 감소, b 증가로 보수적 표현
A_LB = A_MEAN * 0.92     # a 하한: 8% 감소
B_LB = B_MEAN * 1.25     # b 하한: 25% 증가 (더 빠른 손상)

print("=" * 65)
print("  [S-N 파라미터]")
print(f"  수종: {SPECIES}  (생재 MOR = {MOR_GRN:.1f} MPa)")
print(f"  평균값: a={A_MEAN}, b={B_MEAN}")
print(f"  하한값: a={A_LB:.3f}, b={B_LB:.3f}  (95% CI, 보수적)")
print("=" * 65)

# ══════════════════════════════════════════════════════════════════════
# 1. ABAQUS 응력 시계열 로드
# ══════════════════════════════════════════════════════════════════════
def load_abaqus_csv(csv_path):
    t_arr, sigma_arr = [], []
    with open(csv_path, 'r', encoding='utf-8-sig') as f:
        lines = f.readlines()
    header_idx = 0
    sigma_col = -1
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith('#') or not stripped:
            continue
        try:
            float(stripped.split(',')[0])
            header_idx = i
            break
        except ValueError:
            cols = [c.strip().lower() for c in stripped.split(',')]
            sigma_col = next((j for j, c in enumerate(cols)
                              if 'sigma' in c or 'stress' in c or 'bending' in c), -1)
            header_idx = i + 1
            break
    for line in lines[header_idx:]:
        stripped = line.strip()
        if not stripped or stripped.startswith('#'):
            continue
        parts = stripped.split(',')
        try:
            t = float(parts[0])
            if sigma_col == -1 or sigma_col >= len(parts):
                sigma_col = len(parts) - 1
            sigma = abs(float(parts[sigma_col]))
            t_arr.append(t)
            sigma_arr.append(sigma)
        except (ValueError, IndexError):
            continue
    return np.array(t_arr), np.array(sigma_arr)

csv_abaqus = out_dir / 'abaqus_stress_history.csv'
if csv_abaqus.exists():
    t_abq, sigma_abq = load_abaqus_csv(csv_abaqus)
    SIGMA_REF   = sigma_abq.max()     # U=20 m/s 기준 최대 응력
    U_REF       = 20.0                # [m/s] 기준 풍속
    T_SIM       = t_abq[-1] - t_abq[0]  # 시뮬레이션 시간
    DATA_SOURCE = "ABAQUS 동적해석"
    print(f"\n  [{DATA_SOURCE}] 로드 완료")
    print(f"  σ_max = {SIGMA_REF:.3f} MPa  /  T_sim = {T_SIM:.0f} s")
else:
    # ABAQUS CSV 없으면 Python 단독 결과로 대체
    from fatigue_analysis import load_drag_force, force_to_stress
    t_py, F_py = load_drag_force(out_dir / 'wind_kaimal.csv')
    D_B = 0.24
    I_B = np.pi * D_B**4 / 64
    C_B = D_B / 2
    sigma_abq   = force_to_stress(F_py, 8.0*(2/3), C_B, I_B)
    t_abq       = t_py
    SIGMA_REF   = sigma_abq.max()
    U_REF       = 20.0
    T_SIM       = 600.0
    DAF_FALLBACK = 1.130
    sigma_abq   = sigma_abq * DAF_FALLBACK   # DAF 보정
    SIGMA_REF   = sigma_abq.max()
    DATA_SOURCE = "Python 단독 + DAF=1.130 (ABAQUS 대체)"
    print(f"\n  [경고] ABAQUS CSV 없음 → {DATA_SOURCE}")
    print(f"  σ_max = {SIGMA_REF:.3f} MPa")

# ══════════════════════════════════════════════════════════════════════
# 2. 단면 특성 (부후율별 응력 증폭 계수)
# ══════════════════════════════════════════════════════════════════════
D_BASE = 0.24   # [m] 근원부 직경

def stress_factor_concentric(decay_ratio):
    """
    동심 부후 → 응력 증폭 계수 = (c/I)_decayed / (c/I)_intact
    건전 원형: I = πR⁴/4, c = R
    부후 중공: I = π(R⁴-r_d⁴)/4, c = R (최외단 동일)
    → factor = R⁴ / (R⁴ - r_d⁴) = 1/(1-dr⁴)
    """
    dr = decay_ratio
    if dr >= 1.0:
        return np.inf
    return 1.0 / (1.0 - dr**4)

# 부후율 스터디 포인트
DECAY_RATIOS = np.array([0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.80])

print("\n  [부후율별 응력 증폭]")
print(f"  {'r_d/R':>7}  {'t/R':>5}  {'응력증폭':>8}  {'σ_max증폭 [MPa]':>16}")
for dr in DECAY_RATIOS:
    sf  = stress_factor_concentric(dr)
    sig = SIGMA_REF * sf
    tR  = 1 - dr
    flag = " ← Mattheck 임계" if abs(dr - 0.7) < 0.01 else ""
    print(f"  {dr:>7.2f}  {tR:>5.2f}  {sf:>8.3f}x  {sig:>14.2f} MPa{flag}")

# ══════════════════════════════════════════════════════════════════════
# 3. Rainflow cycle 추출 (U_REF 기준)
# ══════════════════════════════════════════════════════════════════════
base_cycles = list(rainflow.extract_cycles(sigma_abq))
base_ranges  = np.array([c[0] for c in base_cycles])   # 응력 범위 [MPa]
base_counts  = np.array([c[2] for c in base_cycles])   # 사이클 수

print(f"\n  Rainflow (기준): {base_counts.sum():.0f} 사이클 / {T_SIM:.0f}s")

# ══════════════════════════════════════════════════════════════════════
# 4. Weibull 풍속 분포 설정
# ══════════════════════════════════════════════════════════════════════
WEIBULL_SITES = {
    "서울 도심":  {"k": 1.6, "c": 4.0,  "color": "steelblue"},
    "인천 해안":  {"k": 2.0, "c": 6.5,  "color": "darkorange"},
    "제주 해안":  {"k": 2.2, "c": 9.0,  "color": "crimson"},
}

def weibull_pdf(U, k, c):
    if U <= 0:
        return 0.0
    return (k / c) * (U / c)**(k-1) * np.exp(-(U/c)**k)

def weibull_mean(k, c):
    return c * gamma_func(1 + 1/k)

U_BINS   = np.arange(1.0, 35.0, 1.0)   # 대표 풍속 구간 중앙값 [m/s]
dU       = 1.0
T_YEAR   = 365.25 * 24 * 3600          # [s/년]

print("\n  [Weibull 지역별 파라미터]")
for site, p in WEIBULL_SITES.items():
    U_m = weibull_mean(p['k'], p['c'])
    print(f"  {site}: k={p['k']}, c={p['c']} → U_mean={U_m:.1f} m/s")

# ══════════════════════════════════════════════════════════════════════
# 5. 피로 손상 계산 함수
# ══════════════════════════════════════════════════════════════════════
def damage_per_sim(ranges, counts, stress_factor, MOR, a, b):
    """
    응력 범위를 stress_factor 배 증폭 후 Miner's damage 계산
    ranges : 기준 응력 범위 배열 [MPa]
    counts : 사이클 수 배열
    stress_factor : 부후에 의한 응력 증폭 배율
    반환: D (T_SIM 동안 누적 손상)
    """
    D = 0.0
    for rng, cnt in zip(ranges, counts):
        amp = (rng * stress_factor) / 2.0
        if amp <= 0:
            continue
        ratio = amp / MOR
        if ratio >= a:
            Nf = 1.0
        else:
            Nf = 10.0**((a - ratio) / b)
        D += cnt / Nf
    return D

def annual_damage_weibull(ranges, counts, stress_factor, MOR, a, b,
                          k_w, c_w, U_ref=20.0, T_sim=600.0):
    """
    Weibull 풍속 분포를 이용한 연간 누적 손상
    - 풍속 U에서 응력은 U_ref 기준 결과를 풍압비(U/U_ref)²로 스케일
    - 각 풍속 출현 시간 = pdf(U) × dU × T_year
    """
    D_annual = 0.0
    for U in U_BINS:
        # 풍압 비례: F ∝ U², σ ∝ U²  → ranges 스케일
        wind_scale = (U / U_ref)**2
        prob       = weibull_pdf(U, k_w, c_w) * dU
        T_U        = prob * T_YEAR    # [s/년] 해당 풍속 출현 시간
        if T_U < 1.0:
            continue
        n_sim = T_U / T_sim           # 해당 풍속 시뮬레이션 반복 횟수/년
        D_sim = damage_per_sim(ranges * wind_scale, counts,
                               stress_factor, MOR, a, b)
        D_annual += D_sim * n_sim
    return D_annual

# ══════════════════════════════════════════════════════════════════════
# 6. 메인 계산: 부후율 × 지역 × S-N 수준 → 피로수명 매트릭스
# ══════════════════════════════════════════════════════════════════════
SN_CASES = {
    "평균 S-N":  (A_MEAN, B_MEAN),
    "하한 S-N":  (A_LB,   B_LB),
}

results = {}   # results[site][sn_case][decay_ratio] = 수명(년)

print("\n" + "=" * 70)
print("  피로수명 매트릭스 (단위: 년)")
print("  수종: 은행나무 생재  /  ABAQUS DAF 포함  /  Weibull 풍속 분포")
print("=" * 70)

for site, wp in WEIBULL_SITES.items():
    results[site] = {}
    for sn_label, (a_sn, b_sn) in SN_CASES.items():
        lives = []
        for dr in DECAY_RATIOS:
            sf   = stress_factor_concentric(dr)
            D_ann = annual_damage_weibull(
                base_ranges, base_counts, sf,
                MOR_GRN, a_sn, b_sn,
                wp['k'], wp['c'], U_REF, T_SIM
            )
            life = 1.0 / D_ann if D_ann > 0 else np.inf
            lives.append(life)
        results[site][sn_label] = lives

# 표 출력
header = f"  {'부후율':>6}  {'t/R':>5}"
for site in WEIBULL_SITES:
    header += f"  {site[:4]+'평균':>10}  {site[:4]+'하한':>10}"
print(header)
print("  " + "-" * (len(header) - 2))

for i, dr in enumerate(DECAY_RATIOS):
    tR   = 1 - dr
    flag = " ←Mattheck" if abs(dr - 0.7) < 0.01 else ""
    row  = f"  {dr:>6.2f}  {tR:>5.2f}"
    for site in WEIBULL_SITES:
        for sn_label in ["평균 S-N", "하한 S-N"]:
            L = results[site][sn_label][i]
            if L >= 1e8:
                row += f"  {'> 1억년':>10}"
            elif L >= 1e6:
                row += f"  {L/1e6:>8.2f}백만"
            elif L >= 1e4:
                row += f"  {L:>10.0f}"
            else:
                row += f"  {L:>10.1f}"
    print(row + flag)
print("=" * 70)

# ══════════════════════════════════════════════════════════════════════
# 7. 핵심 시나리오 요약
# ══════════════════════════════════════════════════════════════════════
print("\n  [핵심 시나리오 요약 - 서울 도심 기준]")
print(f"  {'시나리오':35}  {'수명(년)':>12}  {'비고'}")
print("  " + "-" * 65)

scenarios = [
    ("건전목  + 평균 S-N", 0.0, "평균 S-N", "서울 도심"),
    ("건전목  + 하한 S-N", 0.0, "하한 S-N", "서울 도심"),
    ("부후 0.5 + 평균 S-N", 0.5, "평균 S-N", "서울 도심"),
    ("부후 0.5 + 하한 S-N", 0.5, "하한 S-N", "서울 도심"),
    ("부후 0.7 + 평균 S-N", 0.7, "평균 S-N", "서울 도심"),
    ("부후 0.7 + 하한 S-N", 0.7, "하한 S-N", "서울 도심"),
    ("부후 0.8 + 평균 S-N", 0.8, "평균 S-N", "서울 도심"),
    ("부후 0.8 + 하한 S-N", 0.8, "하한 S-N", "서울 도심"),
]
for label, dr, sn_label, site in scenarios:
    idx  = np.argmin(np.abs(DECAY_RATIOS - dr))
    L    = results[site][sn_label][idx]
    note = "★ Mattheck 임계" if dr >= 0.7 else ""
    if L >= 1e7:
        Lstr = f"> {L:.0e} 년"
    else:
        Lstr = f"{L:>12.0f} 년"
    print(f"  {label:35}  {Lstr:>15}  {note}")

# ══════════════════════════════════════════════════════════════════════
# 8. 시각화
# ══════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(2, 3, figsize=(17, 10))
fig.suptitle(
    f"부후 + S-N 하한값 + Weibull 조합 피로수명 평가\n"
    f"{SPECIES}, 생재 (MOR={MOR_GRN:.0f} MPa), ABAQUS DAF 포함",
    fontsize=11
)

# ── 8-a. 서울 도심: 부후율 vs 수명 (평균/하한 S-N 비교) ──────────
ax = axes[0, 0]
site = "서울 도심"
lives_mean = results[site]["평균 S-N"]
lives_lb   = results[site]["하한 S-N"]

# inf 처리 (1e8 상한)
lm = [min(L, 1e8) for L in lives_mean]
ll = [min(L, 1e8) for L in lives_lb]

ax.semilogy(DECAY_RATIOS, lm, 'o-', lw=2, color='steelblue',
            markersize=7, label='평균 S-N')
ax.semilogy(DECAY_RATIOS, ll, 's--', lw=2, color='crimson',
            markersize=7, label='하한 S-N (95% CI)')
ax.fill_between(DECAY_RATIOS, ll, lm, alpha=0.15, color='steelblue',
                label='불확실성 범위')
ax.axvline(0.7, color='red', lw=1.5, ls=':', label='Mattheck 임계 (0.7)')
ax.axhline(50, color='gray', lw=1, ls='--', label='목표수명 50년')
ax.set_xlabel('부후율 r_d/R')
ax.set_ylabel('피로수명 [년]')
ax.set_title(f'{site}: 부후율 vs 피로수명\n(S-N 평균 vs 하한값)')
ax.legend(fontsize=8)
ax.grid(True, which='both', alpha=0.3)

# ── 8-b. 3개 지역 비교 (하한 S-N 기준) ──────────────────────────
ax = axes[0, 1]
for site, wp in WEIBULL_SITES.items():
    ll = [min(L, 1e8) for L in results[site]["하한 S-N"]]
    ax.semilogy(DECAY_RATIOS, ll, 'o-', lw=2,
                color=wp['color'], markersize=6, label=site)
ax.axvline(0.7, color='red', lw=1.5, ls=':', label='Mattheck 임계')
ax.axhline(50, color='gray', lw=1, ls='--', label='목표수명 50년')
ax.set_xlabel('부후율 r_d/R')
ax.set_ylabel('피로수명 [년]')
ax.set_title('지역별 피로수명 비교\n(하한 S-N 기준)')
ax.legend(fontsize=8)
ax.grid(True, which='both', alpha=0.3)

# ── 8-c. 열지도: 부후율 × 지역 (하한 S-N) ──────────────────────
ax = axes[0, 2]
site_labels = list(WEIBULL_SITES.keys())
life_matrix = np.array([
    [min(results[site]["하한 S-N"][i], 1e8)
     for i, _ in enumerate(DECAY_RATIOS)]
    for site in site_labels
])
life_log = np.log10(np.clip(life_matrix, 1, 1e8))

im = ax.imshow(life_log, cmap='RdYlGn', aspect='auto', origin='lower',
               vmin=1, vmax=8)
ax.set_xticks(range(len(DECAY_RATIOS)))
ax.set_xticklabels([f'{dr:.2f}' for dr in DECAY_RATIOS], fontsize=7, rotation=45)
ax.set_yticks(range(len(site_labels)))
ax.set_yticklabels(site_labels, fontsize=8)
ax.set_xlabel('부후율 r_d/R')
ax.set_title('피로수명 열지도 log10(년)\n(하한 S-N, 초록=장수명)')
cbar = plt.colorbar(im, ax=ax)
cbar.set_label('log₁₀(수명 [년])')
cbar.set_ticks([1, 2, 3, 4, 5, 6, 7, 8])
cbar.set_ticklabels(['10','100','1k','10k','100k','1M','10M','>1억'])
for i in range(len(site_labels)):
    for j in range(len(DECAY_RATIOS)):
        L = life_matrix[i, j]
        txt = f'{L:.0f}' if L < 1000 else (f'{L/1000:.0f}k' if L < 1e6 else f'{L/1e6:.1f}M')
        ax.text(j, i, txt, ha='center', va='center', fontsize=6,
                color='white' if life_log[i,j] < 4 else 'black')

# ── 8-d. Weibull 풍속 PDF ────────────────────────────────────────
ax = axes[1, 0]
U_plot = np.linspace(0, 35, 200)
for site, wp in WEIBULL_SITES.items():
    pdf_vals = [weibull_pdf(U, wp['k'], wp['c']) for U in U_plot]
    U_m = weibull_mean(wp['k'], wp['c'])
    ax.plot(U_plot, pdf_vals, lw=2, color=wp['color'],
            label=f"{site} (k={wp['k']}, c={wp['c']}, U_m={U_m:.1f}m/s)")
ax.axvline(U_REF, color='black', lw=1.5, ls='--', label=f'U_ref={U_REF}m/s (ABAQUS)')
ax.set_xlabel('풍속 U [m/s]')
ax.set_ylabel('확률밀도 f(U)')
ax.set_title('Weibull 풍속 확률밀도 분포')
ax.legend(fontsize=7)
ax.grid(True, alpha=0.3)

# ── 8-e. 응력 증폭 vs 부후율 ─────────────────────────────────────
ax = axes[1, 1]
dr_fine = np.linspace(0, 0.92, 300)
sf_fine = [stress_factor_concentric(dr) for dr in dr_fine]
ax.plot(dr_fine, sf_fine, lw=2.5, color='darkorange')
ax.axvline(0.7, color='red', lw=1.5, ls='--', label='Mattheck 임계 (r_d/R=0.7)')
for dr_mark in [0.0, 0.5, 0.7, 0.8]:
    sf_m = stress_factor_concentric(dr_mark)
    sig_m = SIGMA_REF * sf_m
    ax.annotate(f'  {dr_mark:.1f}: ×{sf_m:.2f}\n  σ={sig_m:.1f}MPa',
                xy=(dr_mark, sf_m), fontsize=7,
                arrowprops=dict(arrowstyle='->', color='gray'),
                xytext=(dr_mark + 0.05, sf_m + 0.5 + dr_mark * 3))
ax.set_xlabel('부후율 r_d/R')
ax.set_ylabel('응력 증폭 배율')
ax.set_title(f'부후에 의한 근원부 응력 증폭\n(기준 σ_max={SIGMA_REF:.1f} MPa, ABAQUS)')
ax.set_ylim([0, min(max(sf_fine), 15)])
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# ── 8-f. 종합 시나리오 수명 막대 ─────────────────────────────────
ax = axes[1, 2]
bar_labels = []
bar_lives  = []
bar_colors = []

key_scenarios = [
    (0.0, "서울 도심", "평균 S-N", "steelblue",   "건전+평균\n서울"),
    (0.0, "서울 도심", "하한 S-N", "lightblue",   "건전+하한\n서울"),
    (0.7, "서울 도심", "평균 S-N", "darkorange",  "부후0.7+평균\n서울"),
    (0.7, "서울 도심", "하한 S-N", "orangered",   "부후0.7+하한\n서울"),
    (0.8, "서울 도심", "평균 S-N", "crimson",     "부후0.8+평균\n서울"),
    (0.8, "서울 도심", "하한 S-N", "darkred",     "부후0.8+하한\n서울"),
    (0.7, "제주 해안", "하한 S-N", "purple",      "부후0.7+하한\n제주"),
    (0.8, "제주 해안", "하한 S-N", "black",       "부후0.8+하한\n제주"),
]
for dr, site, sn_label, color, label in key_scenarios:
    idx = np.argmin(np.abs(DECAY_RATIOS - dr))
    L = min(results[site][sn_label][idx], 1e8)
    bar_labels.append(label)
    bar_lives.append(max(L, 0.1))
    bar_colors.append(color)

bars = ax.bar(range(len(bar_labels)), bar_lives,
              color=bar_colors, alpha=0.85, width=0.6)
for bar, L in zip(bars, bar_lives):
    txt = f'{L:.0f}년' if L < 1000 else (f'{L/1000:.0f}k년' if L < 1e6 else f'{L/1e6:.1f}M년')
    ax.text(bar.get_x() + bar.get_width()/2,
            bar.get_height() * 1.1, txt,
            ha='center', va='bottom', fontsize=7, rotation=30)
ax.set_yscale('log')
ax.axhline(50, color='gray', lw=1.5, ls='--', label='목표수명 50년')
ax.set_xticks(range(len(bar_labels)))
ax.set_xticklabels(bar_labels, fontsize=7)
ax.set_ylabel('피로수명 [년]')
ax.set_title('주요 시나리오 피로수명 비교\n(부후율 × 지역 × S-N 수준)')
ax.legend(fontsize=8)
ax.grid(True, axis='y', alpha=0.3)

plt.tight_layout()
fig_path = out_dir / 'decay_coupled_fatigue.png'
plt.savefig(fig_path, dpi=150, bbox_inches='tight')
plt.show()
print(f'\n  그래프 저장: {fig_path}')

# ══════════════════════════════════════════════════════════════════════
# 9. 결론 요약 출력
# ══════════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("  [최종 결론 요약]")
print("=" * 65)
print(f"  기준 조건: {SPECIES}, 생재 MOR={MOR_GRN:.0f}MPa")
print(f"  ABAQUS DAF 포함 σ_max={SIGMA_REF:.2f}MPa  ({DATA_SOURCE})")
print()

for site in ["서울 도심", "인천 해안", "제주 해안"]:
    print(f"  [{site}]")
    for sn_label in ["평균 S-N", "하한 S-N"]:
        print(f"    {sn_label}:")
        for dr in [0.0, 0.5, 0.7, 0.8]:
            idx = np.argmin(np.abs(DECAY_RATIOS - dr))
            L   = results[site][sn_label][idx]
            sf  = stress_factor_concentric(dr)
            sig = SIGMA_REF * sf
            flag = " ★★" if dr >= 0.7 else ""
            if L >= 1e7:
                Lstr = f"> 1,000만년"
            elif L >= 1e4:
                Lstr = f"{L:>12.0f}년"
            else:
                Lstr = f"{L:>12.1f}년"
            print(f"      r_d/R={dr:.1f} (σ={sig:5.1f}MPa): {Lstr}{flag}")
    print()

print("  핵심 발견:")
print("  1. 건전목: 수백만~수억년 → 고주기 피로는 지배 모드 아님")
print("  2. 부후 0.7(Mattheck 임계) + 하한 S-N → 수명 급격 감소")
print("  3. 부후 0.8 + 하한 S-N + 강풍지역 → 수십~수백년 도달 가능")
print("  4. 실제 도복 경로: 부후에 의한 강성 저하가 선행 필수 조건")
print("=" * 65)
