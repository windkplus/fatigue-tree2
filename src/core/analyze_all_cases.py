"""
analyze_all_cases.py
--------------------
ABAQUS 6케이스 결과 통합 피로해석
  - 6케이스 응력 시계열 로드
  - Weibull 풍속 분포 기반 연간 손상 계산
  - 부후율별 응력 증폭
  - 피로수명 매트릭스 (케이스 x 부후율 x 지역)
  - 종합 시각화

입력: output/tree_H*m_DBH*cm_stress.csv  (6개)
      output/abaqus_cases_summary.csv
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
# 1. 설정
# ======================================================================
SPECIES  = "은행나무 (Ginkgo biloba)"
MOR_DRY  = STATIC_PROPS[SPECIES]["MOR_MPa"]       # 68 MPa
MOR_GRN  = MOR_DRY * MC_CORRECTION["MOR"]          # 40.8 MPa
SN       = SN_PARAMS[SPECIES]
A_MEAN, B_MEAN = SN["a_R0"], SN["b_R0"]            # 0.97, 0.072
A_LB,   B_LB   = A_MEAN * 0.92, B_MEAN * 1.25      # 하한값

U_REF   = 20.0          # [m/s] ABAQUS 기준 풍속
T_SIM   = 600.0         # [s]   시뮬레이션 시간

DECAY_RATIOS = np.array([0.0, 0.3, 0.5, 0.6, 0.7, 0.75, 0.80])

WEIBULL_SITES = {
    "서울 도심":  {"k": 1.6, "c": 4.0,  "color": "steelblue"},
    "인천 해안":  {"k": 2.0, "c": 6.5,  "color": "darkorange"},
    "제주 해안":  {"k": 2.2, "c": 9.0,  "color": "crimson"},
}

U_BINS  = np.arange(1.0, 35.0, 1.0)
T_YEAR  = 365.25 * 24 * 3600

SN_CASES = {"평균 S-N": (A_MEAN, B_MEAN), "하한 S-N": (A_LB, B_LB)}

# ======================================================================
# 2. 케이스 정의
# ======================================================================
CASES = [
    {"name": "tree_H6m_DBH15cm", "H": 6.0, "DBH": 0.15,
     "label": "H6m-DBH15", "color": "#e74c3c", "ls": "--"},
    {"name": "tree_H6m_DBH20cm", "H": 6.0, "DBH": 0.20,
     "label": "H6m-DBH20", "color": "#e67e22", "ls": "--"},
    {"name": "tree_H6m_DBH25cm", "H": 6.0, "DBH": 0.25,
     "label": "H6m-DBH25", "color": "#f1c40f", "ls": "--"},
    {"name": "tree_H8m_DBH15cm", "H": 8.0, "DBH": 0.15,
     "label": "H8m-DBH15", "color": "#8e44ad", "ls": "-"},
    {"name": "tree_H8m_DBH20cm", "H": 8.0, "DBH": 0.20,
     "label": "H8m-DBH20 (기준)", "color": "#2980b9", "ls": "-"},
    {"name": "tree_H8m_DBH25cm", "H": 8.0, "DBH": 0.25,
     "label": "H8m-DBH25", "color": "#27ae60", "ls": "-"},
]

# ======================================================================
# 3. 함수 정의
# ======================================================================
def load_stress_csv(csv_path):
    """응력 시계열 CSV 로드 -> sigma_bending_MPa 배열"""
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


def stress_factor(decay_ratio):
    """동심 부후 응력 증폭: 1/(1-dr^4)"""
    dr = decay_ratio
    return 1.0 / (1.0 - dr**4) if dr < 1.0 else np.inf


def weibull_pdf(U, k, c):
    if U <= 0: return 0.0
    return (k/c) * (U/c)**(k-1) * np.exp(-(U/c)**k)


def weibull_mean(k, c):
    return c * gamma_func(1 + 1/k)


def damage_per_sim(ranges, counts, sf, MOR, a, b):
    """Miner's rule 손상 계산 (응력 증폭 sf 적용)"""
    D = 0.0
    for rng, cnt in zip(ranges, counts):
        amp = (rng * sf) / 2.0
        if amp <= 0: continue
        ratio = amp / MOR
        Nf = 1.0 if ratio >= a else 10.0**((a - ratio) / b)
        D += cnt / Nf
    return D


def annual_damage_weibull(ranges, counts, sf, MOR, a, b, k_w, c_w):
    """Weibull 가중 연간 손상"""
    D_ann = 0.0
    for U in U_BINS:
        wind_scale = (U / U_REF)**2
        prob = weibull_pdf(U, k_w, c_w) * 1.0
        T_U  = prob * T_YEAR
        if T_U < 1.0: continue
        n_sim = T_U / T_SIM
        D_sim = damage_per_sim(ranges * wind_scale, counts, sf, MOR, a, b)
        D_ann += D_sim * n_sim
    return D_ann

# ======================================================================
# 4. 응력 시계열 로드 및 Rainflow
# ======================================================================
print("=" * 65)
print(f"  통합 피로해석: {SPECIES} 생재 (MOR={MOR_GRN:.0f} MPa)")
print("=" * 65)

case_data = {}
for c in CASES:
    csv_path = out_dir / f"{c['name']}_stress.csv"
    if not csv_path.exists():
        print(f"  [경고] {csv_path.name} 없음 - 건너뜀")
        continue
    sigma = load_stress_csv(csv_path)
    cycles = list(rainflow.extract_cycles(sigma))
    ranges = np.array([x[0] for x in cycles])
    counts = np.array([x[2] for x in cycles])
    sigma_max = sigma.max()
    case_data[c['name']] = {
        "sigma": sigma, "ranges": ranges, "counts": counts,
        "sigma_max": sigma_max, "n_cycles": counts.sum(),
        **c
    }
    print(f"  {c['label']:20}  σ_max={sigma_max:6.2f} MPa  "
          f"사이클={counts.sum():.0f}")

print()

# ======================================================================
# 5. 피로수명 계산 (케이스 x 부후율 x 지역 x S-N)
# ======================================================================
# results[case_name][site][sn_label][decay_idx] = 수명(년)
results = {}

for cname, cd in case_data.items():
    results[cname] = {}
    for site, wp in WEIBULL_SITES.items():
        results[cname][site] = {}
        for sn_label, (a_sn, b_sn) in SN_CASES.items():
            lives = []
            for dr in DECAY_RATIOS:
                sf = stress_factor(dr)
                D_ann = annual_damage_weibull(
                    cd["ranges"], cd["counts"], sf,
                    MOR_GRN, a_sn, b_sn, wp["k"], wp["c"]
                )
                life = 1.0 / D_ann if D_ann > 0 else np.inf
                lives.append(life)
            results[cname][site][sn_label] = lives

# ======================================================================
# 6. 결과 테이블 출력 (서울, 하한 S-N 기준)
# ======================================================================
site_ref = "서울 도심"
sn_ref   = "하한 S-N"

print("=" * 70)
print(f"  피로수명 [년] - {site_ref}, {sn_ref}")
print(f"  {'케이스':20}", end="")
for dr in DECAY_RATIOS:
    print(f"  {f'dr={dr:.2f}':>10}", end="")
print()
print("  " + "-" * 68)

for c in CASES:
    cname = c["name"]
    if cname not in results: continue
    print(f"  {c['label']:20}", end="")
    for i, dr in enumerate(DECAY_RATIOS):
        L = results[cname][site_ref][sn_ref][i]
        if L >= 1e7:
            s = ">1000만"
        elif L >= 1e4:
            s = f"{L:.0f}"
        elif L >= 1:
            s = f"{L:.1f}"
        else:
            s = f"{L:.2f}"
        print(f"  {s:>10}", end="")
    print()
print("=" * 70)

# ======================================================================
# 7. 시각화
# ======================================================================
fig, axes = plt.subplots(2, 3, figsize=(17, 11))
fig.suptitle(
    f"ABAQUS 6케이스 통합 피로해석\n"
    f"{SPECIES}, 생재 MOR={MOR_GRN:.0f}MPa, Weibull 풍속분포, 다축하중",
    fontsize=11
)

# ── 7-a. σ_max 비교 막대 ────────────────────────────────────────────
ax = axes[0, 0]
case_labels = [cd["label"] for cd in case_data.values()]
sigma_maxes = [cd["sigma_max"] for cd in case_data.values()]
colors_bar  = [cd["color"] for cd in case_data.values()]
bars = ax.bar(range(len(case_labels)), sigma_maxes,
              color=colors_bar, alpha=0.8, width=0.6)
ax.axhline(MOR_GRN, color='red', lw=2, ls='--',
           label=f'생재 MOR={MOR_GRN:.0f}MPa')
ax.axhline(MOR_GRN * SN["sigma_limit_ratio"], color='orange', lw=1.5, ls=':',
           label=f'피로한계={MOR_GRN*SN["sigma_limit_ratio"]:.0f}MPa')
for bar, s in zip(bars, sigma_maxes):
    ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+1,
            f'{s:.1f}', ha='center', va='bottom', fontsize=8)
ax.set_xticks(range(len(case_labels)))
ax.set_xticklabels(case_labels, rotation=30, fontsize=7, ha='right')
ax.set_ylabel('σ_max [MPa]  (U=20m/s)')
ax.set_title('케이스별 최대 응력 (U=20m/s)\n빨간선=MOR, 주황선=피로한계')
ax.legend(fontsize=8)
ax.grid(True, axis='y', alpha=0.3)

# ── 7-b. f1 비교 ────────────────────────────────────────────────────
ax = axes[0, 1]
summary_path = out_dir / 'abaqus_cases_summary.csv'
f1_vals, case_names_sum = [], []
if summary_path.exists():
    with open(summary_path) as f:
        reader = csv.DictReader(f)
        seen = set()
        for row in reader:
            n = row['ODB']
            if n not in seen:
                f1_vals.append(float(row['f1_Hz']))
                case_names_sum.append(n.replace('tree_','').replace('_',' '))
                seen.add(n)

if f1_vals:
    colors_f1 = [cd["color"] for cd in list(case_data.values())[:len(f1_vals)]]
    bars2 = ax.bar(range(len(case_names_sum)), f1_vals,
                   color=colors_f1, alpha=0.8, width=0.6)
    ax.axhspan(0.2, 0.8, alpha=0.1, color='green', label='수목 f1 문헌 범위')
    for bar, f in zip(bars2, f1_vals):
        ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.01,
                f'{f:.3f}', ha='center', va='bottom', fontsize=8)
    ax.set_xticks(range(len(case_names_sum)))
    ax.set_xticklabels(case_names_sum, rotation=30, fontsize=7, ha='right')
ax.set_ylabel('1차 고유진동수 f1 [Hz]')
ax.set_title('케이스별 고유진동수\n(높을수록 강성 높음)')
ax.legend(fontsize=8)
ax.grid(True, axis='y', alpha=0.3)

# ── 7-c. Weibull PDF ────────────────────────────────────────────────
ax = axes[0, 2]
U_plot = np.linspace(0, 35, 300)
for site, wp in WEIBULL_SITES.items():
    pdf_v = [weibull_pdf(U, wp['k'], wp['c']) for U in U_plot]
    U_m   = weibull_mean(wp['k'], wp['c'])
    ax.plot(U_plot, pdf_v, lw=2, color=wp['color'],
            label=f"{site} (U_m={U_m:.1f}m/s)")
ax.axvline(U_REF, color='black', lw=1.5, ls='--', label=f'U_ref={U_REF}m/s')
ax.set_xlabel('풍속 [m/s]')
ax.set_ylabel('확률밀도 f(U)')
ax.set_title('Weibull 풍속 분포\n(U=20m/s는 극히 희귀)')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# ── 7-d. 부후율 vs 수명 (서울, 하한 S-N, 전 케이스) ──────────────
ax = axes[1, 0]
for c in CASES:
    cname = c["name"]
    if cname not in results: continue
    lives = [min(results[cname][site_ref][sn_ref][i], 1e8)
             for i in range(len(DECAY_RATIOS))]
    ax.semilogy(DECAY_RATIOS, lives, marker='o', lw=2,
                color=c["color"], ls=c["ls"], markersize=5,
                label=c["label"])
ax.axvline(0.7, color='red', lw=1.5, ls=':', label='Mattheck 임계')
ax.axhline(50, color='gray', lw=1, ls='--', label='목표수명 50년')
ax.set_xlabel('부후율 r_d/R')
ax.set_ylabel('피로수명 [년]')
ax.set_title(f'부후율 vs 피로수명\n({site_ref}, {sn_ref})')
ax.legend(fontsize=7, ncol=2)
ax.grid(True, which='both', alpha=0.3)

# ── 7-e. 지역별 비교 (건전목, 하한 S-N) ────────────────────────────
ax = axes[1, 1]
x = np.arange(len(CASES))
width = 0.25
site_list = list(WEIBULL_SITES.keys())
for j, (site, wp) in enumerate(WEIBULL_SITES.items()):
    lives_healthy = []
    for c in CASES:
        cname = c["name"]
        if cname not in results:
            lives_healthy.append(0.1)
            continue
        L = min(results[cname][site]["하한 S-N"][0], 1e8)
        lives_healthy.append(max(L, 0.1))
    ax.bar(x + j*width, lives_healthy, width,
           color=wp["color"], alpha=0.8, label=site)

ax.set_yscale('log')
ax.axhline(50, color='gray', lw=1, ls='--', label='50년')
ax.set_xticks(x + width)
ax.set_xticklabels([cd["label"] for cd in case_data.values()],
                   rotation=30, fontsize=7, ha='right')
ax.set_ylabel('피로수명 [년]')
ax.set_title('건전목 피로수명 비교\n(하한 S-N, 지역별)')
ax.legend(fontsize=8)
ax.grid(True, axis='y', alpha=0.3)

# ── 7-f. 열지도: 케이스 x 부후율 (서울, 하한 S-N) ──────────────────
ax = axes[1, 2]
life_mat = np.zeros((len(CASES), len(DECAY_RATIOS)))
valid_cases = [c for c in CASES if c["name"] in results]
for i, c in enumerate(valid_cases):
    for j in range(len(DECAY_RATIOS)):
        L = results[c["name"]][site_ref][sn_ref][j]
        life_mat[i, j] = np.log10(max(min(L, 1e8), 0.1))

im = ax.imshow(life_mat, cmap='RdYlGn', aspect='auto',
               vmin=0, vmax=8, origin='upper')
ax.set_xticks(range(len(DECAY_RATIOS)))
ax.set_xticklabels([f'{dr:.2f}' for dr in DECAY_RATIOS], fontsize=8)
ax.set_yticks(range(len(valid_cases)))
ax.set_yticklabels([c["label"] for c in valid_cases], fontsize=8)
ax.set_xlabel('부후율 r_d/R')
ax.set_title(f'피로수명 열지도 log10(년)\n({site_ref}, 하한 S-N)')
cbar = plt.colorbar(im, ax=ax)
cbar.set_label('log₁₀(수명 [년])')
cbar.set_ticks([0,2,4,6,8])
cbar.set_ticklabels(['<1년','100년','10k년','1M년','>1억년'])
for i in range(len(valid_cases)):
    for j in range(len(DECAY_RATIOS)):
        L = 10**life_mat[i, j]
        if L < 1:    txt = f'{L:.1f}'
        elif L < 100: txt = f'{L:.0f}'
        elif L < 1e4: txt = f'{L/1000:.0f}k'
        elif L < 1e6: txt = f'{L/1000:.0f}k'
        else:         txt = f'{L/1e6:.0f}M'
        ax.text(j, i, txt, ha='center', va='center', fontsize=6,
                color='white' if life_mat[i,j] < 3 else 'black')

plt.tight_layout()
fig_path = out_dir / 'analyze_all_cases.png'
plt.savefig(fig_path, dpi=150, bbox_inches='tight')
plt.show()
print(f'\n  그래프 저장: {fig_path}')

# ======================================================================
# 8. 최종 요약
# ======================================================================
print("\n" + "=" * 65)
print("  최종 결론 요약")
print("=" * 65)
print(f"  수종: {SPECIES}  생재 MOR={MOR_GRN:.0f}MPa")
print()
for c in CASES:
    cname = c["name"]
    if cname not in results: continue
    cd = case_data[cname]
    print(f"  [{c['label']}]  f1={cd.get('f1_Hz','-')}  sigma_max={cd['sigma_max']:.1f}MPa")
    for site in ["서울 도심", "제주 해안"]:
        for sn_label in ["평균 S-N", "하한 S-N"]:
            L0 = results[cname][site][sn_label][0]   # 건전목
            L7 = results[cname][site][sn_label][
                list(DECAY_RATIOS).index(0.7) if 0.7 in DECAY_RATIOS
                else np.argmin(np.abs(DECAY_RATIOS-0.7))]
            def fmt(L):
                if L>=1e7: return ">1000만년"
                if L>=1e4: return f"{L:.0f}년"
                return f"{L:.1f}년"
            print(f"    {site} {sn_label}: 건전={fmt(L0)}  부후0.7={fmt(L7)}")
    print()
print("=" * 65)
