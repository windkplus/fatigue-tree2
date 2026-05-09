"""
gumbel_extreme_wind.py
----------------------
Gumbel 분포 기반 극한 풍속(재현주기 풍속) 해석
→ 단기 극한 하중에 의한 정적 파단 위험도 평가
→ 피로 해석과 병행하여 "어느 파괴 모드가 지배적인가" 판단

배경:
  가로수 도복 사고는 두 가지 메커니즘 중 하나:
    A. 반복 피로 누적 (Miner's rule)   ← 앞선 분석
    B. 극한 돌풍에 의한 즉각 파단      ← 이 스크립트

  재현주기 T_R 의 설계풍속 U_R:
    U_R = U_n - (1/α) × ln(-ln(1 - 1/T_R))   (Gumbel)
    또는  P(U > U_R) = 1 - exp(-exp(-α(U - u_n)))

참고:
  Gumbel, E.J. (1958). Statistics of Extremes. Columbia Univ. Press.
  Cook, N.J. (1985). The Designer's Guide to Wind Loading. BRE/Butterworths.
  KDS 41 10 15 (2022). 건축구조기준 — 풍하중.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import csv
from pathlib import Path
from scipy.stats import gumbel_r
from wood_properties_db import STATIC_PROPS, SN_PARAMS, get_props

matplotlib.rcParams['font.family'] = 'Malgun Gothic'
matplotlib.rcParams['axes.unicode_minus'] = False

out_dir = Path("output")
out_dir.mkdir(exist_ok=True)

# ══════════════════════════════════════════════════════════════════════
# 1. 한국 기상청 연최대풍속 대표 데이터 (지역별)
#    출처: 기상청 기후통계 (1991~2020 평년값 기반 유추)
#    실제 연구에서는 대상 지점 실측 데이터 사용 권장
# ══════════════════════════════════════════════════════════════════════
WIND_SITES = {
    "서울 (도심)": {
        "annual_max_ms": [9.2, 10.1, 8.7, 11.3, 9.8, 12.5, 10.4, 8.9,
                          11.8, 10.7, 9.4, 13.2, 10.0, 9.6, 11.1, 12.0,
                          9.9, 10.8, 11.5, 9.3, 10.6, 12.1, 9.7, 11.9,
                          10.2, 9.5, 12.8, 10.9, 11.4, 10.3],
        "color": "steelblue"
    },
    "인천 (해안)": {
        "annual_max_ms": [13.5, 15.2, 12.8, 16.7, 14.3, 17.9, 15.1, 13.2,
                          16.4, 15.8, 14.0, 18.6, 14.8, 14.2, 16.0, 17.5,
                          14.6, 15.9, 16.8, 13.8, 15.4, 17.2, 14.4, 17.1,
                          15.0, 14.1, 18.2, 16.1, 16.5, 15.3],
        "color": "darkorange"
    },
    "제주 (강풍)": {
        "annual_max_ms": [22.0, 25.5, 20.1, 28.3, 23.7, 31.2, 24.8, 21.6,
                          27.5, 26.1, 22.9, 33.4, 24.0, 23.3, 26.8, 29.7,
                          24.5, 27.0, 28.6, 22.5, 25.8, 30.1, 23.8, 29.4,
                          25.2, 23.1, 32.0, 27.3, 28.0, 26.4],
        "color": "crimson"
    },
}

# ══════════════════════════════════════════════════════════════════════
# 2. Gumbel 분포 피팅 및 재현주기 풍속 계산
# ══════════════════════════════════════════════════════════════════════
RETURN_PERIODS = [2, 5, 10, 25, 50, 100, 200, 500]  # [년]

def fit_gumbel(data):
    """Gumbel 최대값 분포 피팅 → (위치모수 u, 척도모수 α)"""
    loc, scale = gumbel_r.fit(data)   # scipy: loc=u, scale=1/α
    return loc, scale

def return_period_wind(loc, scale, T_R):
    """재현주기 T_R 에 대한 설계풍속 [m/s]"""
    return gumbel_r.ppf(1 - 1/T_R, loc=loc, scale=scale)

def exceedance_prob(loc, scale, U):
    """풍속 U 의 연간 초과확률"""
    return 1 - gumbel_r.cdf(U, loc=loc, scale=scale)

print("=" * 60)
print("  Gumbel 분포 피팅 및 재현주기 풍속")
print("=" * 60)

gumbel_params = {}
for site, data in WIND_SITES.items():
    arr = np.array(data["annual_max_ms"])
    loc, scale = fit_gumbel(arr)
    gumbel_params[site] = {"loc": loc, "scale": scale,
                            "data": arr, "color": data["color"]}
    print(f"\n  [{site}]")
    print(f"    Gumbel 위치모수 u    = {loc:.2f} m/s")
    print(f"    Gumbel 척도모수 α⁻¹ = {scale:.2f} m/s")
    print(f"    데이터 평균 / 표준편차 = {arr.mean():.1f} / {arr.std():.1f} m/s")
    print(f"    {'T_R':>6}  {'U_R':>8}")
    for T_R in [10, 50, 100, 500]:
        U_R = return_period_wind(loc, scale, T_R)
        print(f"    {T_R:>6}년  {U_R:>7.1f} m/s")

# ══════════════════════════════════════════════════════════════════════
# 3. 극한 풍속에 의한 정적 파단 위험 평가
#    σ_max = 0.5 × ρ × Cd × A × U² × L_arm × c / I
#    파단 조건: σ_max ≥ MOR × 안전율(=1.0 극한)
# ══════════════════════════════════════════════════════════════════════

# 수목 파라미터 (기준)
H_TREE   = 8.0
DBH      = 0.20
D_BASE   = DBH * 1.2
I_SEC    = np.pi * D_BASE**4 / 64
C        = D_BASE / 2
L_ARM    = H_TREE * (2/3)
A_CROWN  = (H_TREE * 0.4)**2 * np.pi / 4
RHO_AIR  = 1.225
CD       = 0.5

def peak_stress(U_peak):
    """극한 풍속에서의 근원부 최대 굽힘 응력 [MPa]"""
    F = 0.5 * RHO_AIR * CD * A_CROWN * U_peak**2
    return F * L_ARM * C / I_SEC * 1e-6

# 파단 풍속 역산: σ = MOR → U_fail = sqrt(2×MOR×I/(ρ×Cd×A×L×c))
def failure_wind(MOR_MPa):
    """정적 파단을 일으키는 임계 풍속 [m/s]"""
    sigma_Pa = MOR_MPa * 1e6
    F_crit   = sigma_Pa * I_SEC / (L_ARM * C)
    U_crit   = np.sqrt(2 * F_crit / (RHO_AIR * CD * A_CROWN))
    return U_crit

print("\n" + "=" * 60)
print("  수종별 정적 파단 임계 풍속 (H=8m, DBH=20cm, 기건)")
print("=" * 60)
failure_winds = {}
for sp, props in STATIC_PROPS.items():
    MOR_dry   = props["MOR_MPa"]
    MOR_green = MOR_dry * 0.60   # 생재 보정
    U_dry   = failure_wind(MOR_dry)
    U_green = failure_wind(MOR_green)
    short = sp.split("(")[0].strip()
    print(f"  {short:<18}  기건 {U_dry:>5.1f} m/s  /  생재 {U_green:>5.1f} m/s")
    failure_winds[sp] = {"dry": U_dry, "green": U_green}

# 재현주기와 파단 위험 교차 분석
print("\n" + "=" * 60)
print("  파단 임계 풍속의 재현주기 (생재 기준)")
print("=" * 60)
print(f"  {'수종':<18}  ", end="")
for site in WIND_SITES:
    print(f"{site.split('(')[0].strip():>12}", end="")
print()
print("-" * 60)

for sp in STATIC_PROPS:
    U_fail = failure_winds[sp]["green"]
    short  = sp.split("(")[0].strip()
    print(f"  {short:<18}  ", end="")
    for site, params in gumbel_params.items():
        p_exceed = exceedance_prob(params["loc"], params["scale"], U_fail)
        if p_exceed > 0:
            T_R_fail = 1.0 / p_exceed
            txt = f"{T_R_fail:>10.0f}년"
        else:
            txt = "     >10000년"
        print(f"{txt:>12}", end="")
    print()

# ══════════════════════════════════════════════════════════════════════
# 4. 피로 vs 극한 파괴: 지배 모드 비교
#    피로수명(Weibull 기반) vs 극한 파괴 재현주기
# ══════════════════════════════════════════════════════════════════════
# 앞선 weibull_fatigue.py 결과 (서울 도심 기준, 은행나무)
FATIGUE_LIFE_REF = {
    "서울 (도심)":  755884,
    "인천 (해안)":  321837,
    "제주 (강풍)":  229340,
}

# ══════════════════════════════════════════════════════════════════════
# 5. 결과 CSV 저장
# ══════════════════════════════════════════════════════════════════════
csv_path = out_dir / "gumbel_results.csv"
with open(csv_path, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["site", "T_R_yr", "U_R_ms"])
    for site, params in gumbel_params.items():
        for T_R in RETURN_PERIODS:
            U_R = return_period_wind(params["loc"], params["scale"], T_R)
            writer.writerow([site, T_R, f"{U_R:.2f}"])
print(f"\n  결과 저장: {csv_path}")

# ══════════════════════════════════════════════════════════════════════
# 6. 시각화
# ══════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(2, 3, figsize=(16, 10))
fig.suptitle("극한 풍속 해석 (Gumbel 분포) 및 파괴 모드 비교\n"
             f"H={H_TREE}m, DBH={DBH*100:.0f}cm", fontsize=11)

# ── 6-a. Gumbel 확률지 ────────────────────────────────────────────
ax = axes[0, 0]
for site, params in gumbel_params.items():
    data_sorted = np.sort(params["data"])
    n = len(data_sorted)
    # Gumbel 확률지 변환: y = -ln(-ln(F))
    F_emp = (np.arange(1, n+1) - 0.44) / (n + 0.12)  # Gringorten 공식
    y_emp = -np.log(-np.log(F_emp))
    loc, scale = params["loc"], params["scale"]
    U_fit = np.linspace(data_sorted.min()*0.9, data_sorted.max()*1.2, 100)
    y_fit = (U_fit - loc) / scale
    ax.scatter(y_emp, data_sorted, s=20, color=params["color"], alpha=0.7)
    ax.plot(y_fit, U_fit, lw=1.5, color=params["color"],
            label=site.split("(")[0].strip())
ax.set_xlabel("Gumbel 환산변량 y = -ln(-ln(F))")
ax.set_ylabel("연최대풍속 [m/s]")
ax.set_title("Gumbel 확률지\n(점: 실측, 선: 피팅)")
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# ── 6-b. 재현주기-풍속 곡선 ──────────────────────────────────────
ax = axes[0, 1]
T_range = np.logspace(0, 3.5, 200)
for site, params in gumbel_params.items():
    U_curve = [return_period_wind(params["loc"], params["scale"], T)
               for T in T_range]
    ax.semilogx(T_range, U_curve, lw=2, color=params["color"],
                label=site.split("(")[0].strip())
# 수종별 파단 풍속 수평선
for sp, winds in failure_winds.items():
    short = sp.split("(")[0].strip()
    ax.axhline(winds["green"], lw=0.8, ls=":", color="gray", alpha=0.6)
    ax.text(1.1, winds["green"]+0.3, f"{short}(생재)", fontsize=6, color="gray")
ax.set_xlabel("재현주기 T_R [년]")
ax.set_ylabel("설계풍속 U_R [m/s]")
ax.set_title("재현주기-설계풍속 곡선\n(수평 점선: 수종별 정적 파단 풍속)")
ax.legend(fontsize=8)
ax.grid(True, which="both", alpha=0.3)
ax.set_xlim([1, 2000])

# ── 6-c. 극한 응력 vs MOR ─────────────────────────────────────────
ax = axes[0, 2]
U_range_plot = np.linspace(5, 60, 200)
sigma_curve  = np.array([peak_stress(U) for U in U_range_plot])
ax.plot(U_range_plot, sigma_curve, lw=2, color="black", label="근원부 최대 응력")

COLORS_SP = ["steelblue","hotpink","darkorange","forestgreen","purple","saddlebrown"]
for (sp, props), color in zip(STATIC_PROPS.items(), COLORS_SP):
    short = sp.split("(")[0].strip()
    MOR_g = props["MOR_MPa"] * 0.60
    ax.axhline(MOR_g, lw=1.2, ls="--", color=color, alpha=0.8,
               label=f"{short} MOR_생재={MOR_g:.0f}MPa")
ax.set_xlabel("풍속 U [m/s]")
ax.set_ylabel("응력 [MPa]")
ax.set_title("극한 풍속-응력 관계\n(수평선: 수종별 생재 MOR)")
ax.legend(fontsize=6.5, loc="upper left")
ax.grid(True, alpha=0.3)
ax.set_xlim([5, 60])

# ── 6-d. DBH별 파단 임계 풍속 ────────────────────────────────────
ax = axes[1, 0]
DBH_range_plot = np.linspace(0.10, 0.40, 100)
for (sp, props), color in zip(STATIC_PROPS.items(), COLORS_SP):
    MOR_g = props["MOR_MPa"] * 0.60
    U_fails = []
    for d in DBH_range_plot:
        D_b = d * 1.2
        I_  = np.pi * D_b**4 / 64
        c_  = D_b / 2
        F_c = MOR_g * 1e6 * I_ / (L_ARM * c_)
        U_f = np.sqrt(2 * F_c / (RHO_AIR * CD * A_CROWN))
        U_fails.append(U_f)
    short = sp.split("(")[0].strip()
    ax.plot(DBH_range_plot*100, U_fails, lw=1.8, color=color, label=short)

# 재현주기 100년 풍속 수평선
for site, params in gumbel_params.items():
    U100 = return_period_wind(params["loc"], params["scale"], 100)
    ax.axhline(U100, lw=1.0, ls=":", color=params["color"],
               label=f"{site.split('(')[0].strip()} U_100={U100:.0f}m/s")
ax.set_xlabel("흉고직경 DBH [cm]")
ax.set_ylabel("파단 임계 풍속 [m/s]")
ax.set_title("DBH별 정적 파단 임계 풍속\n(점선: 지역별 100년 재현 풍속)")
ax.legend(fontsize=6, ncol=2)
ax.grid(True, alpha=0.3)

# ── 6-e. 피로 vs 극한: 지배 파괴 모드 ───────────────────────────
ax = axes[1, 1]
sites_short = [s.split("(")[0].strip() for s in WIND_SITES]
x = np.arange(len(sites_short))
w = 0.35

# 피로수명 (은행나무 기준)
fatigue_lives = [FATIGUE_LIFE_REF[s] for s in WIND_SITES]
# 극한 파단 재현주기 (은행나무, 생재)
sp_ref = "은행나무 (Ginkgo biloba)"
U_fail_ginkgo = failure_winds[sp_ref]["green"]
extreme_Trs = []
for site, params in gumbel_params.items():
    p = exceedance_prob(params["loc"], params["scale"], U_fail_ginkgo)
    extreme_Trs.append(min(1/p if p > 0 else 1e7, 1e7))

b1 = ax.bar(x - w/2, fatigue_lives, w, color="steelblue",
            alpha=0.8, label="피로수명 [년]")
b2 = ax.bar(x + w/2, extreme_Trs, w, color="crimson",
            alpha=0.8, label="극한 파단 재현주기 [년]")
ax.set_yscale("log")
ax.set_xticks(x)
ax.set_xticklabels(sites_short, fontsize=9)
ax.set_ylabel("년 [log scale]")
ax.set_title(f"피로 수명 vs 극한 파단 재현주기\n({sp_ref.split('(')[0].strip()}, 생재)")
ax.legend(fontsize=8)
ax.grid(True, axis="y", alpha=0.3)

# 지배 모드 텍스트
for i, (fl, et) in enumerate(zip(fatigue_lives, extreme_Trs)):
    mode = "극한 지배" if et < fl else "피로 지배"
    ax.text(i, max(fl, et)*1.5, mode, ha="center", fontsize=8,
            color="crimson" if mode=="극한 지배" else "steelblue")

# ── 6-f. 풍속 초과확률 곡선 ──────────────────────────────────────
ax = axes[1, 2]
U_range2 = np.linspace(5, 60, 300)
for site, params in gumbel_params.items():
    p_exceed = [exceedance_prob(params["loc"], params["scale"], U)
                for U in U_range2]
    ax.semilogy(U_range2, p_exceed, lw=2, color=params["color"],
                label=site.split("(")[0].strip())
for (sp, winds), color in zip(failure_winds.items(), COLORS_SP):
    ax.axvline(winds["green"], lw=0.8, ls=":", color=color, alpha=0.6)
ax.set_xlabel("풍속 U [m/s]")
ax.set_ylabel("연간 초과확률")
ax.set_title("풍속별 연간 초과확률\n(수직 점선: 수종별 생재 파단 풍속)")
ax.legend(fontsize=8)
ax.grid(True, which="both", alpha=0.3)
ax.set_xlim([5, 60])
ax.set_ylim([1e-4, 1])

plt.tight_layout()
fig_path = out_dir / "gumbel_extreme_wind.png"
plt.savefig(fig_path, dpi=150, bbox_inches="tight")
plt.show()
print(f"  그래프 저장: {fig_path}")
