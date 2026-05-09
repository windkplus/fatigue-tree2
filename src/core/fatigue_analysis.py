"""
fatigue_analysis.py
-------------------
가로수 줄기 근원부 피로수명 평가
  1. 풍하중 시간이력 → 보 굽힘 이론으로 응력 시계열 변환
  2. Rainflow counting (ASTM E1049)
  3. 목재 S-N 곡선 (문헌 다수 수종)
  4. Miner's rule 누적손상 D, 피로수명 예측

참고 문헌:
  Bohannan, B. (1966). Effect of size on bending strength of wood members.
    USDA Forest Service Research Paper FPL 56.
  Tsai, K.T. & Ansell, M.P. (1990). The fatigue properties of wood in flexure.
    J. Mater. Sci., 25, 865-878.
  Kohara, M. & Okuyama, T. (1992). Bending fatigue of small clear wood specimens.
    Mokuzai Gakkaishi, 38(6), 529-536.
  EN 1991-1-4 (Eurocode 1, 2005). Wind actions.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import csv
import rainflow
from pathlib import Path

matplotlib.rcParams['font.family'] = 'Malgun Gothic'
matplotlib.rcParams['axes.unicode_minus'] = False

# ══════════════════════════════════════════════════════════════════════
# 1. 수목 형상 파라미터
# ══════════════════════════════════════════════════════════════════════
H_tree  = 8.0        # [m]   수고 (전체 높이)
DBH     = 0.20       # [m]   흉고직경 (직경, 1.2 m 높이)
D_base  = DBH * 1.2  # [m]   근원부 직경 (흉고의 약 1.2배 가정)

# 근원부 단면 2차 모멘트 (원형 단면)
I_base  = np.pi * D_base**4 / 64.0   # [m⁴]
c_base  = D_base / 2.0               # [m]  중립축 ~ 최외단 거리

# 등가 모멘트 팔 길이: 수관 합력 작용점 ≈ 수고의 2/3
L_arm   = H_tree * (2.0 / 3.0)       # [m]

print("=" * 55)
print("  수목 형상 파라미터")
print("=" * 55)
print(f"  수고       H        = {H_tree:.1f} m")
print(f"  흉고직경   DBH      = {DBH*100:.1f} cm")
print(f"  근원부직경 D_base   = {D_base*100:.2f} cm")
print(f"  단면 2차 모멘트 I   = {I_base*1e6:.4f} cm⁴  ({I_base:.4e} m⁴)")
print(f"  모멘트 팔  L_arm    = {L_arm:.2f} m")

# ══════════════════════════════════════════════════════════════════════
# 2. 목재 S-N 곡선 데이터베이스 (문헌 기반)
# ══════════════════════════════════════════════════════════════════════
# 형식: S_a/MOR = a - b * log10(N_f)   →   N_f = 10^((a - S_a/MOR) / b)
# MOR: Modulus of Rupture (파단 계수, 정적 굽힘 강도) [MPa]
# a  : N=1 절편 (≈ 1.0)
# b  : 기울기 (클수록 피로 저항 낮음)
# 응력비 R ≈ 0 (편진, 풍하중은 단방향)
#
# 주의: 실제 연구에서는 대상 수종 실측값으로 교체 필요
WOOD_SN = {
    "은행나무 (Ginkgo biloba)": {
        "MOR_MPa": 68.0,   # Panshin & de Zeeuw (1980) 참고
        "a": 0.97,
        "b": 0.072,
        "color": "steelblue",
        "source": "Kohara & Okuyama (1992) 유추"
    },
    "벚나무 (Prunus serrulata)": {
        "MOR_MPa": 88.0,
        "a": 0.98,
        "b": 0.065,
        "color": "hotpink",
        "source": "Tsai & Ansell (1990) 유추"
    },
    "플라타너스 (Platanus x acerifolia)": {
        "MOR_MPa": 72.0,
        "a": 0.97,
        "b": 0.070,
        "color": "darkorange",
        "source": "Bohannan (1966) 유추"
    },
    "느티나무 (Zelkova serrata)": {
        "MOR_MPa": 95.0,
        "a": 0.98,
        "b": 0.060,
        "color": "forestgreen",
        "source": "Kohara & Okuyama (1992) 유추"
    },
}

def sn_cycles_to_failure(sigma_a_MPa, MOR_MPa, a, b):
    """
    S-N 곡선 기반 파단 반복수 계산
    N_f = 10^( (a - σ_a/MOR) / b )
    σ_a > a*MOR (정적 파단 초과) → N_f = 1 반환
    σ_a ≤ 0                       → N_f = inf 반환
    """
    if sigma_a_MPa <= 0:
        return np.inf
    ratio = sigma_a_MPa / MOR_MPa
    if ratio >= a:
        return 1.0
    return 10.0 ** ((a - ratio) / b)

# ══════════════════════════════════════════════════════════════════════
# 3. 풍하중 → 응력 시계열 변환
# ══════════════════════════════════════════════════════════════════════
def load_drag_force(csv_path):
    """output/wind_*.csv 에서 항력 시계열 로드"""
    t_list, F_list = [], []
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row or row[0].startswith("#") or row[0] == "time_s":
                continue
            t_list.append(float(row[0]))
            F_list.append(float(row[3]))   # F_drag_N 열
    return np.array(t_list), np.array(F_list)


def force_to_stress(F_t, L_arm, c, I):
    """
    단순 외팔보 가정: 근원부 굽힘 응력
        M(t) = F(t) * L_arm
        σ(t) = M(t) * c / I    [Pa → MPa]
    """
    M_t     = F_t * L_arm           # [N·m]
    sigma_t = M_t * c / I           # [Pa]
    return sigma_t * 1e-6            # [MPa]


# 두 스펙트럼 모델 모두 처리
out_dir   = Path("output")
datasets  = {}

for model, fname in [("Davenport", "wind_davenport.csv"),
                     ("Kaimal",    "wind_kaimal.csv")]:
    csv_path = out_dir / fname
    if not csv_path.exists():
        print(f"  [경고] {csv_path} 없음 — wind_time_history.py 먼저 실행하세요.")
        continue
    t_arr, F_arr = load_drag_force(csv_path)
    sigma_arr    = force_to_stress(F_arr, L_arm, c_base, I_base)
    datasets[model] = {"t": t_arr, "F": F_arr, "sigma": sigma_arr}
    print(f"\n  [{model}] 응력 통계")
    print(f"    σ 최대 = {sigma_arr.max():.2f} MPa")
    print(f"    σ 평균 = {sigma_arr.mean():.2f} MPa")
    print(f"    σ std  = {sigma_arr.std():.2f} MPa")

if not datasets:
    raise SystemExit("데이터 없음: wind_time_history.py 를 먼저 실행하세요.")

# ══════════════════════════════════════════════════════════════════════
# 4. Rainflow Counting (ASTM E1049)
# ══════════════════════════════════════════════════════════════════════
def run_rainflow(sigma_t):
    """
    rainflow 패키지로 사이클 추출.
    반환: [(range_MPa, mean_MPa, count, i_start, i_end), ...]
    amplitude = range / 2
    """
    cycles = list(rainflow.extract_cycles(sigma_t))
    # cycles: (range, mean, count, i_start, i_end)
    return cycles


print("\n" + "=" * 55)
print("  Rainflow Counting 결과")
print("=" * 55)

rf_results = {}
for model, ds in datasets.items():
    cycles = run_rainflow(ds["sigma"])
    rf_results[model] = cycles
    ranges   = np.array([c[0] for c in cycles])
    means    = np.array([c[1] for c in cycles])
    counts   = np.array([c[2] for c in cycles])
    amps     = ranges / 2.0
    print(f"\n  [{model}]")
    print(f"    총 사이클 수     = {counts.sum():.1f}")
    print(f"    응력 진폭 최대   = {amps.max():.3f} MPa")
    print(f"    응력 진폭 평균   = {np.average(amps, weights=counts):.3f} MPa")

# ══════════════════════════════════════════════════════════════════════
# 5. Miner's Rule 누적 손상 계산
# ══════════════════════════════════════════════════════════════════════
def miners_rule(cycles, MOR_MPa, a, b, T_sim=600.0, T_year=3.15e7):
    """
    Miner's rule:  D = Σ (n_i / N_fi)

    T_sim  : 시뮬레이션 시간 [s]  (기본 600 s = 10분)
    T_year : 연간 강풍 노출 시간 [s]
             (예: 연 100시간 강풍 → 3.6e5 s)
             기본값 3.15e7 = 1년 전체 (보수적 상한)

    반환: D_sim (10분), D_annual (연간), L_fatigue (피로수명, 년)
    """
    D_sim = 0.0
    for rng, mean, count, *_ in cycles:
        amp = rng / 2.0
        if amp <= 0:
            continue
        Nf  = sn_cycles_to_failure(amp, MOR_MPa, a, b)
        D_sim += count / Nf

    scale      = T_year / T_sim
    D_annual   = D_sim * scale
    L_fatigue  = 1.0 / D_annual if D_annual > 0 else np.inf
    return D_sim, D_annual, L_fatigue


print("\n" + "=" * 55)
print("  Miner's Rule 피로수명 평가")
print(f"  (수목 형상: H={H_tree}m, DBH={DBH*100:.0f}cm)")
print("=" * 55)

# 연간 강풍 노출 시간 가정 (논문에서 감도 분석 필요)
T_YEAR_EXPOSURE = 100 * 3600.0   # [s] = 연 100시간 강풍 노출 (보수적)

damage_table = {}   # {(model, species): (D_sim, D_annual, L_yr)}

for model, cycles in rf_results.items():
    print(f"\n  ── {model} ──")
    damage_table[model] = {}
    for species, params in WOOD_SN.items():
        D_sim, D_ann, L_yr = miners_rule(
            cycles,
            MOR_MPa = params["MOR_MPa"],
            a       = params["a"],
            b       = params["b"],
            T_sim   = 600.0,
            T_year  = T_YEAR_EXPOSURE
        )
        damage_table[model][species] = (D_sim, D_ann, L_yr)
        print(f"  {species[:18]:18s}  "
              f"D_10min={D_sim:.2e}  "
              f"D_annual={D_ann:.4f}  "
              f"수명={L_yr:.1f} 년")

# ══════════════════════════════════════════════════════════════════════
# 6. 결과 저장 (CSV)
# ══════════════════════════════════════════════════════════════════════
result_path = out_dir / "fatigue_results.csv"
with open(result_path, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["model", "species", "MOR_MPa",
                     "D_10min", "D_annual", "fatigue_life_yr"])
    for model in damage_table:
        for species, (D_sim, D_ann, L_yr) in damage_table[model].items():
            MOR = WOOD_SN[species]["MOR_MPa"]
            writer.writerow([model, species, MOR,
                             f"{D_sim:.4e}", f"{D_ann:.6f}", f"{L_yr:.2f}"])
print(f"\n  결과 저장: {result_path}")

# ══════════════════════════════════════════════════════════════════════
# 7. 시각화
# ══════════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(15, 11))
fig.suptitle(
    f"가로수 근원부 피로수명 평가\n"
    f"H={H_tree}m, DBH={DBH*100:.0f}cm, U_bar=20m/s, I_u=0.20  "
    f"(연간 강풍 노출 {T_YEAR_EXPOSURE/3600:.0f}h 가정)",
    fontsize=11)

gs = fig.add_gridspec(2, 3, hspace=0.45, wspace=0.35)

# ── 7-a. 응력 시계열 (처음 60초) ─────────────────────────────────────
ax1 = fig.add_subplot(gs[0, :2])
colors_model = {"Davenport": "steelblue", "Kaimal": "darkorange"}
for model, ds in datasets.items():
    mask = ds["t"] <= 60
    ax1.plot(ds["t"][mask], ds["sigma"][mask],
             lw=0.7, color=colors_model[model], label=model, alpha=0.85)
ax1.set_xlabel("시간 [s]")
ax1.set_ylabel("근원부 굽힘 응력 [MPa]")
ax1.set_title("근원부 응력 시계열 (처음 60 s)")
ax1.legend(fontsize=8)
ax1.grid(True, alpha=0.3)

# ── 7-b. Rainflow 진폭 히스토그램 ────────────────────────────────────
ax2 = fig.add_subplot(gs[0, 2])
bins = np.linspace(0, max(
    max(c[0]/2 for c in rf_results["Davenport"]),
    max(c[0]/2 for c in rf_results["Kaimal"])
) * 1.05, 30)

for model, cycles in rf_results.items():
    amps   = np.array([c[0]/2 for c in cycles])
    counts = np.array([c[2]   for c in cycles])
    hist, edges = np.histogram(amps, bins=bins, weights=counts)
    centers = (edges[:-1] + edges[1:]) / 2
    ax2.bar(centers, hist, width=np.diff(edges),
            color=colors_model[model], alpha=0.5, label=model)
ax2.set_xlabel("응력 진폭 [MPa]")
ax2.set_ylabel("사이클 수")
ax2.set_title("Rainflow 사이클 분포")
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.3)

# ── 7-c. S-N 곡선 + Rainflow 사이클 오버레이 ─────────────────────────
ax3 = fig.add_subplot(gs[1, :2])
N_plot = np.logspace(1, 9, 500)
for species, params in WOOD_SN.items():
    MOR, a, b = params["MOR_MPa"], params["a"], params["b"]
    # S-N 역산: σ_a = MOR * (a - b*log10(N_f))
    sigma_sn = MOR * (a - b * np.log10(N_plot))
    valid = sigma_sn > 0
    ax3.semilogx(N_plot[valid], sigma_sn[valid],
                 lw=1.5, color=params["color"], label=f"{species} (MOR={MOR}MPa)")

# Rainflow 사이클을 S-N 위에 점으로 표시 (Kaimal 기준)
if "Kaimal" in rf_results:
    for rng, mean, count, *_ in rf_results["Kaimal"]:
        amp = rng / 2.0
        if amp < 0.001:
            continue
        # 해당 진폭에서의 파단 반복수 (느티나무 기준)
        p = WOOD_SN["느티나무 (Zelkova serrata)"]
        Nf = sn_cycles_to_failure(amp, p["MOR_MPa"], p["a"], p["b"])
        if np.isfinite(Nf) and Nf < 1e9:
            ax3.scatter(Nf, amp, s=count*2, color="gray", alpha=0.3, zorder=3)

ax3.set_xlabel("파단 반복수 N_f")
ax3.set_ylabel("응력 진폭 σ_a [MPa]")
ax3.set_title("목재 S-N 곡선  (점: Kaimal Rainflow 사이클, 크기∝횟수)")
ax3.legend(fontsize=7, loc="upper right")
ax3.grid(True, which="both", alpha=0.3)
ax3.set_xlim([1e1, 1e9])

# ── 7-d. 피로수명 막대그래프 비교 ────────────────────────────────────
ax4 = fig.add_subplot(gs[1, 2])
species_labels  = [s.split("(")[0].strip() for s in WOOD_SN.keys()]
x = np.arange(len(WOOD_SN))
width = 0.35

for i, model in enumerate(["Davenport", "Kaimal"]):
    lives = [damage_table[model][sp][2] for sp in WOOD_SN]
    # 무한대 수명은 표시를 위해 상한 clamp
    lives_plot = [min(L, 500) for L in lives]
    bars = ax4.bar(x + (i - 0.5) * width, lives_plot, width,
                   color=colors_model[model], alpha=0.8, label=model)
    # 실제 수치 표시
    for bar, L in zip(bars, lives):
        label_txt = f"{L:.0f}" if L < 500 else ">500"
        ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 3,
                 label_txt, ha="center", va="bottom", fontsize=7)

ax4.set_xticks(x)
ax4.set_xticklabels(species_labels, rotation=25, ha="right", fontsize=7)
ax4.set_ylabel("피로수명 [년]")
ax4.set_title(f"수종별 피로수명 비교\n(연간 강풍 노출 {T_YEAR_EXPOSURE/3600:.0f}h)")
ax4.legend(fontsize=8)
ax4.grid(True, axis="y", alpha=0.3)

fig_path = out_dir / "fatigue_analysis.png"
plt.savefig(fig_path, dpi=150, bbox_inches="tight")
plt.show()
print(f"  그래프 저장: {fig_path}")
