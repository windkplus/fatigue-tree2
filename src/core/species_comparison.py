"""
species_comparison.py
---------------------
6수종 피로수명 비교 (ABAQUS 재해석 없이)
  - 동적 응력 시계열: 은행나무 ABAQUS 6케이스 결과 그대로 사용
  - 수종별 변경 물성: MOR, S-N 파라미터 (피로강도)
  - E, 밀도 차이에 의한 동적 응답 변화는 민감도 분석으로 별도 반영

가정:
  "같은 구조(H, DBH, 풍하중), 다른 재료강도" 비교
  논문 표현: ABAQUS 동적 응력 시계열은 은행나무 기준 적용,
             수종별 피로강도(MOR, S-N)만 변화
"""

import sys, os
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

# ======================================================================
# 1. 수종 목록
# ======================================================================
SPECIES_LIST = [
    {"name": "은행나무 (Ginkgo biloba)",
     "short": "은행나무", "color": "#2ecc71", "marker": "o"},
    {"name": "왕벚나무 (Prunus yedoensis)",
     "short": "왕벚나무", "color": "#e74c3c", "marker": "s"},
    {"name": "플라타너스 (Platanus x acerifolia)",
     "short": "플라타너스", "color": "#3498db", "marker": "^"},
    {"name": "느티나무 (Zelkova serrata)",
     "short": "느티나무", "color": "#8e44ad", "marker": "D"},
    {"name": "이팝나무 (Chionanthus retusus)",
     "short": "이팝나무", "color": "#e67e22", "marker": "v"},
    {"name": "메타세쿼이아 (Metasequoia glyptostroboides)",
     "short": "메타세쿼이아", "color": "#7f8c8d", "marker": "P"},
]

# ======================================================================
# 2. ABAQUS 케이스 정의
# ======================================================================
CASES = [
    {"name": "tree_H6m_DBH15cm", "label": "H6m-D15", "H": 6.0, "DBH": 0.15},
    {"name": "tree_H6m_DBH20cm", "label": "H6m-D20", "H": 6.0, "DBH": 0.20},
    {"name": "tree_H6m_DBH25cm", "label": "H6m-D25", "H": 6.0, "DBH": 0.25},
    {"name": "tree_H8m_DBH15cm", "label": "H8m-D15", "H": 8.0, "DBH": 0.15},
    {"name": "tree_H8m_DBH20cm", "label": "H8m-D20 (기준)", "H": 8.0, "DBH": 0.20},
    {"name": "tree_H8m_DBH25cm", "label": "H8m-D25", "H": 8.0, "DBH": 0.25},
]

# ======================================================================
# 3. 공통 설정
# ======================================================================
U_REF  = 20.0
T_SIM  = 600.0
T_YEAR = 365.25 * 24 * 3600

WEIBULL_SITES = {
    "서울 도심":  {"k": 1.6, "c": 4.0},
    "인천 해안":  {"k": 2.0, "c": 6.5},
    "제주 해안":  {"k": 2.2, "c": 9.0},
}
U_BINS = np.arange(1.0, 35.0, 1.0)

DECAY_RATIOS = np.array([0.0, 0.3, 0.5, 0.7, 0.8])

# ======================================================================
# 4. 함수
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

def stress_factor(dr):
    return 1.0 / (1.0 - dr**4) if dr < 1.0 else np.inf

def weibull_pdf(U, k, c):
    if U <= 0: return 0.0
    return (k/c) * (U/c)**(k-1) * np.exp(-(U/c)**k)

def weibull_mean(k, c):
    return c * gamma_func(1 + 1/k)

def annual_damage(ranges, counts, sf, MOR, a, b, k_w, c_w):
    D_ann = 0.0
    for U in U_BINS:
        ws   = (U / U_REF)**2
        prob = weibull_pdf(U, k_w, c_w) * 1.0
        T_U  = prob * T_YEAR
        if T_U < 1.0: continue
        n_sim = T_U / T_SIM
        D = 0.0
        for rng, cnt in zip(ranges * ws, counts):
            amp = (rng * sf) / 2.0
            if amp <= 0: continue
            ratio = amp / MOR
            Nf = 1.0 if ratio >= a else 10.0**((a - ratio) / b)
            D += cnt / Nf
        D_ann += D * n_sim
    return D_ann

def fatigue_life(ranges, counts, sf, MOR, a, b, k_w, c_w):
    D = annual_damage(ranges, counts, sf, MOR, a, b, k_w, c_w)
    return 1.0 / D if D > 0 else np.inf

# ======================================================================
# 5. 응력 시계열 로드 및 Rainflow
# ======================================================================
print("=" * 65)
print("  6수종 피로수명 비교")
print("  (동적 응력: 은행나무 ABAQUS 결과 공용, MOR/S-N만 수종별 적용)")
print("=" * 65)

case_rf = {}   # case_rf[case_name] = (ranges, counts, sigma_max)
for c in CASES:
    csv_path = out_dir / f"{c['name']}_stress.csv"
    if not csv_path.exists():
        print(f"  [경고] {csv_path.name} 없음")
        continue
    sigma  = load_stress_csv(csv_path)
    cycles = list(rainflow.extract_cycles(sigma))
    ranges = np.array([x[0] for x in cycles])
    counts = np.array([x[2] for x in cycles])
    case_rf[c["name"]] = (ranges, counts, sigma.max())
    print(f"  {c['label']:18}  σ_max={sigma.max():6.2f} MPa  "
          f"사이클={counts.sum():.0f}")

# ======================================================================
# 6. 수종별 물성 출력
# ======================================================================
print()
print("=" * 70)
print(f"  {'수종':12}  {'MOR(건)':>8}  {'MOR(생)':>8}  "
      f"{'a_R0':>6}  {'b_R0':>6}  {'피로한계':>8}")
print("  " + "-" * 60)
for sp in SPECIES_LIST:
    sname = sp["name"]
    p  = STATIC_PROPS[sname]
    sn = SN_PARAMS.get(sname, {})
    MOR_g = p["MOR_MPa"] * MC_CORRECTION["MOR"]
    a  = sn.get("a_R0", "-")
    b  = sn.get("b_R0", "-")
    fl = MOR_g * sn.get("sigma_limit_ratio", 0.4)
    print(f"  {sp['short']:12}  {p['MOR_MPa']:>8.1f}  {MOR_g:>8.1f}  "
          f"{a:>6}  {b:>6}  {fl:>8.2f} MPa")
print("=" * 70)

# ======================================================================
# 7. 피로수명 계산: 대표 케이스 H8m-DBH20cm, 서울 도심
# ======================================================================
REP_CASE = "tree_H8m_DBH20cm"
REP_SITE = "서울 도심"
wp_ref   = WEIBULL_SITES[REP_SITE]

if REP_CASE in case_rf:
    ranges_ref, counts_ref, _ = case_rf[REP_CASE]

    print(f"\n  [대표 케이스: H8m-DBH20cm, {REP_SITE}]")
    print(f"  {'수종':12}  {'MOR_생 MPa':>10}  "
          f"{'건전목(년)':>12}  {'부후0.5(년)':>12}  "
          f"{'부후0.7(년)':>12}  {'부후0.8(년)':>12}")
    print("  " + "-" * 72)

    sp_lives_ref = {}
    for sp in SPECIES_LIST:
        sname = sp["name"]
        p     = STATIC_PROPS[sname]
        sn    = SN_PARAMS.get(sname, {})
        if not sn: continue
        MOR_g = p["MOR_MPa"] * MC_CORRECTION["MOR"]
        a_sn  = sn["a_R0"]
        b_sn  = sn["b_R0"]

        lives = []
        for dr in [0.0, 0.5, 0.7, 0.8]:
            sf = stress_factor(dr)
            L  = fatigue_life(ranges_ref, counts_ref, sf,
                              MOR_g, a_sn, b_sn,
                              wp_ref["k"], wp_ref["c"])
            lives.append(L)
        sp_lives_ref[sname] = lives

        def fmt(L):
            if L >= 1e7: return ">1000만"
            if L >= 1e4: return f"{L:.0f}"
            if L >= 1:   return f"{L:.1f}"
            return f"{L:.2f}"

        print(f"  {sp['short']:12}  {MOR_g:>10.1f}  "
              f"{fmt(lives[0]):>12}  {fmt(lives[1]):>12}  "
              f"{fmt(lives[2]):>12}  {fmt(lives[3]):>12}")

# ======================================================================
# 8. 전 케이스 x 전 수종 (건전목, 서울, 평균 S-N)
# ======================================================================
print(f"\n  [전 케이스 x 수종 피로수명 - 건전목, {REP_SITE}, 평균 S-N]  단위: 년")
header = f"  {'수종':12}"
for c in CASES:
    header += f"  {c['label']:>12}"
print(header)
print("  " + "-" * (14 + 14*len(CASES)))

all_lives = {}   # all_lives[sp_short][case_name] = life
for sp in SPECIES_LIST:
    sname  = sp["name"]
    p      = STATIC_PROPS[sname]
    sn     = SN_PARAMS.get(sname, {})
    if not sn: continue
    MOR_g  = p["MOR_MPa"] * MC_CORRECTION["MOR"]
    a_sn   = sn["a_R0"]
    b_sn   = sn["b_R0"]
    all_lives[sp["short"]] = {}

    row = f"  {sp['short']:12}"
    for c in CASES:
        cname = c["name"]
        if cname not in case_rf:
            row += f"  {'N/A':>12}"
            continue
        rng, cnt, _ = case_rf[cname]
        L = fatigue_life(rng, cnt, 1.0, MOR_g, a_sn, b_sn,
                         wp_ref["k"], wp_ref["c"])
        all_lives[sp["short"]][cname] = L
        if L >= 1e7: s = ">1000만"
        elif L >= 1e4: s = f"{L:.0f}"
        elif L >= 1: s = f"{L:.1f}"
        else: s = f"{L:.2f}"
        row += f"  {s:>12}"
    print(row)

# ======================================================================
# 9. 시각화
# ======================================================================
fig, axes = plt.subplots(2, 3, figsize=(17, 11))
fig.suptitle(
    "6수종 피로수명 비교\n"
    "(은행나무 ABAQUS 응력 공용, 수종별 MOR/S-N 적용, 서울 도심 Weibull)",
    fontsize=11
)

# ── 9-a. 수종별 MOR 비교 ──────────────────────────────────────────
ax = axes[0, 0]
sp_names  = [sp["short"] for sp in SPECIES_LIST]
MOR_dry   = [STATIC_PROPS[sp["name"]]["MOR_MPa"] for sp in SPECIES_LIST]
MOR_green = [m * 0.60 for m in MOR_dry]
colors_sp = [sp["color"] for sp in SPECIES_LIST]
x = np.arange(len(sp_names))
b1 = ax.bar(x - 0.2, MOR_dry,   0.35, label='소건재', color=colors_sp, alpha=0.5)
b2 = ax.bar(x + 0.2, MOR_green, 0.35, label='생재', color=colors_sp, alpha=0.9)
for bar, v in zip(b2, MOR_green):
    ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.5,
            f'{v:.0f}', ha='center', va='bottom', fontsize=7)
ax.set_xticks(x)
ax.set_xticklabels(sp_names, rotation=30, fontsize=8, ha='right')
ax.set_ylabel('MOR [MPa]')
ax.set_title('수종별 파단계수 (MOR)\n(흐린=소건재, 진한=생재)')
ax.legend(fontsize=8)
ax.grid(True, axis='y', alpha=0.3)

# ── 9-b. 대표 케이스(H8m-D20) 수종별 수명 vs 부후율 ──────────────
ax = axes[0, 1]
if REP_CASE in case_rf:
    for sp in SPECIES_LIST:
        sname = sp["name"]
        p     = STATIC_PROPS[sname]
        sn    = SN_PARAMS.get(sname, {})
        if not sn: continue
        MOR_g = p["MOR_MPa"] * MC_CORRECTION["MOR"]
        lives = []
        for dr in DECAY_RATIOS:
            sf = stress_factor(dr)
            L  = min(fatigue_life(ranges_ref, counts_ref, sf, MOR_g,
                                  sn["a_R0"], sn["b_R0"],
                                  wp_ref["k"], wp_ref["c"]), 1e8)
            lives.append(L)
        ax.semilogy(DECAY_RATIOS, lives, marker=sp["marker"], lw=2,
                    color=sp["color"], markersize=6, label=sp["short"])
ax.axvline(0.7, color='red', lw=1.5, ls=':', label='Mattheck 임계')
ax.axhline(50,  color='gray', lw=1, ls='--', label='목표수명 50년')
ax.set_xlabel('부후율 r_d/R')
ax.set_ylabel('피로수명 [년]')
ax.set_title(f'수종별 부후율-수명 곡선\n(H8m-DBH20cm, {REP_SITE})')
ax.legend(fontsize=7, ncol=2)
ax.grid(True, which='both', alpha=0.3)

# ── 9-c. 수종 x 케이스 열지도 (건전목, 서울, 평균 S-N) ───────────
ax = axes[0, 2]
valid_cases = [c for c in CASES if c["name"] in case_rf]
life_mat = np.zeros((len(SPECIES_LIST), len(valid_cases)))
for i, sp in enumerate(SPECIES_LIST):
    for j, c in enumerate(valid_cases):
        L = all_lives.get(sp["short"], {}).get(c["name"], 0.1)
        life_mat[i, j] = np.log10(max(min(L, 1e8), 0.1))

im = ax.imshow(life_mat, cmap='RdYlGn', aspect='auto',
               vmin=0, vmax=8, origin='upper')
ax.set_xticks(range(len(valid_cases)))
ax.set_xticklabels([c["label"] for c in valid_cases],
                   rotation=30, fontsize=7, ha='right')
ax.set_yticks(range(len(SPECIES_LIST)))
ax.set_yticklabels(sp_names, fontsize=8)
ax.set_title('피로수명 열지도 log10(년)\n(건전목, 서울, 평균 S-N)')
cbar = plt.colorbar(im, ax=ax)
cbar.set_label('log10(년)')
cbar.set_ticks([0,2,4,6,8])
cbar.set_ticklabels(['<1','100','10k','1M','>1억'])
for i in range(len(SPECIES_LIST)):
    for j in range(len(valid_cases)):
        L = 10**life_mat[i, j]
        if L < 1:    txt = f'{L:.1f}'
        elif L < 100: txt = f'{L:.0f}'
        elif L < 1e4: txt = f'{L/1000:.0f}k'
        else:         txt = f'{L/1e6:.0f}M'
        ax.text(j, i, txt, ha='center', va='center', fontsize=6,
                color='white' if life_mat[i,j] < 3 else 'black')

# ── 9-d. 건전목 수명 막대 (H8m-D20, 3개 지역) ────────────────────
ax = axes[1, 0]
if REP_CASE in case_rf:
    x = np.arange(len(SPECIES_LIST))
    width = 0.25
    site_colors = {"서울 도심": "steelblue",
                   "인천 해안": "darkorange", "제주 해안": "crimson"}
    for j, (site, wp) in enumerate(WEIBULL_SITES.items()):
        lives_site = []
        for sp in SPECIES_LIST:
            sname = sp["name"]
            p  = STATIC_PROPS[sname]
            sn = SN_PARAMS.get(sname, {})
            if not sn:
                lives_site.append(0.1); continue
            MOR_g = p["MOR_MPa"] * MC_CORRECTION["MOR"]
            L = min(fatigue_life(ranges_ref, counts_ref, 1.0, MOR_g,
                                 sn["a_R0"], sn["b_R0"],
                                 wp["k"], wp["c"]), 1e8)
            lives_site.append(max(L, 0.1))
        ax.bar(x + j*width, lives_site, width,
               color=site_colors[site], alpha=0.8, label=site)
    ax.set_yscale('log')
    ax.axhline(50, color='gray', lw=1, ls='--', label='50년')
    ax.set_xticks(x + width)
    ax.set_xticklabels(sp_names, rotation=30, fontsize=8, ha='right')
    ax.set_ylabel('피로수명 [년]')
    ax.set_title(f'건전목 수명 지역별 비교\n(H8m-DBH20cm, 평균 S-N)')
    ax.legend(fontsize=8)
    ax.grid(True, axis='y', alpha=0.3)

# ── 9-e. 부후 0.7 수명 막대 (서울, 평균 S-N) ─────────────────────
ax = axes[1, 1]
if REP_CASE in case_rf:
    lives_d07 = []
    for sp in SPECIES_LIST:
        sname = sp["name"]
        p  = STATIC_PROPS[sname]
        sn = SN_PARAMS.get(sname, {})
        if not sn:
            lives_d07.append(0.1); continue
        MOR_g = p["MOR_MPa"] * MC_CORRECTION["MOR"]
        sf = stress_factor(0.7)
        L  = min(fatigue_life(ranges_ref, counts_ref, sf, MOR_g,
                              sn["a_R0"], sn["b_R0"],
                              wp_ref["k"], wp_ref["c"]), 1e8)
        lives_d07.append(max(L, 0.1))
    bars = ax.bar(range(len(sp_names)), lives_d07,
                  color=colors_sp, alpha=0.85, width=0.6)
    for bar, L in zip(bars, lives_d07):
        txt = f'{L:.1f}년' if L < 100 else f'{L:.0f}년'
        ax.text(bar.get_x()+bar.get_width()/2,
                bar.get_height()*1.1, txt,
                ha='center', va='bottom', fontsize=8)
    ax.set_yscale('log')
    ax.axhline(50, color='gray', lw=1, ls='--', label='50년')
    ax.set_xticks(range(len(sp_names)))
    ax.set_xticklabels(sp_names, rotation=30, fontsize=8, ha='right')
    ax.set_ylabel('피로수명 [년]')
    ax.set_title(f'부후 0.7 수명 비교\n(H8m-DBH20cm, {REP_SITE}, 평균 S-N)')
    ax.legend(fontsize=8)
    ax.grid(True, axis='y', alpha=0.3)

# ── 9-f. S-N 곡선 비교 ───────────────────────────────────────────
ax = axes[1, 2]
N_arr = np.logspace(3, 8, 200)
for sp in SPECIES_LIST:
    sname = sp["name"]
    p  = STATIC_PROPS[sname]
    sn = SN_PARAMS.get(sname, {})
    if not sn: continue
    MOR_g = p["MOR_MPa"] * MC_CORRECTION["MOR"]
    a, b  = sn["a_R0"], sn["b_R0"]
    sigma_a = MOR_g * (a - b * np.log10(N_arr))
    sigma_a = np.clip(sigma_a, 0, None)
    ax.semilogx(N_arr, sigma_a, lw=2,
                color=sp["color"], label=f"{sp['short']} (MOR={MOR_g:.0f}MPa)")
ax.set_xlabel('파단 반복수 N_f')
ax.set_ylabel('허용 응력 진폭 [MPa]')
ax.set_title('수종별 S-N 곡선 비교\n(생재 기준, R=0)')
ax.legend(fontsize=7)
ax.grid(True, which='both', alpha=0.3)

plt.tight_layout()
fig_path = out_dir / 'species_comparison.png'
plt.savefig(fig_path, dpi=150, bbox_inches='tight')
plt.show()
print(f'\n  그래프 저장: {fig_path}')
