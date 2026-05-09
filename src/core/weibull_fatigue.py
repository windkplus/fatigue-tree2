"""
weibull_fatigue.py
------------------
Weibull 풍속 분포를 이용한 장기 피로수명 평가
 - 연간 풍속 통계를 Weibull 분포로 표현
 - 풍속 구간별 Rainflow + Miner's rule 적용
 - 전체 연간 누적손상 D_annual = Σ D(U_i) × P(U_i) × T_year

앞선 코드(param_study.py)의 단점:
  "강풍(20 m/s) 100시간/년" 이라는 단순 가정 대신
  실제 풍속 출현 빈도(Weibull 분포)를 반영해 현실적인 수명 예측

참고:
  Weibull, W. (1951). A statistical distribution function of wide applicability.
    J. Appl. Mech., 18, 293-297.
  Manwell et al. (2009). Wind Energy Explained. Wiley. Ch.2.
  IEC 61400-1 (2019). Wind energy — Wind turbine design requirements.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import rainflow
import csv
from pathlib import Path
from scipy.stats import weibull_min
from scipy.special import gamma as gamma_func

matplotlib.rcParams['font.family'] = 'Malgun Gothic'
matplotlib.rcParams['axes.unicode_minus'] = False

out_dir = Path("output")
out_dir.mkdir(exist_ok=True)

# ══════════════════════════════════════════════════════════════════════
# 1. Weibull 파라미터 — 한국 도시 지역 대표값
#    f(U) = (k/c) * (U/c)^(k-1) * exp(-(U/c)^k)
#    평균풍속: U_mean = c * Γ(1 + 1/k)
# ══════════════════════════════════════════════════════════════════════
WEIBULL_SITES = {
    "서울 도심 (약풍)":   {"k": 1.6, "c": 4.0,  "color": "steelblue"},
    "인천 해안 (중풍)":   {"k": 2.0, "c": 6.5,  "color": "darkorange"},
    "제주 해안 (강풍)":   {"k": 2.2, "c": 9.0,  "color": "crimson"},
}

def weibull_pdf(U, k, c):
    """Weibull PDF f(U)"""
    return (k / c) * (U / c)**(k - 1) * np.exp(-(U / c)**k)

def weibull_mean(k, c):
    return c * gamma_func(1 + 1/k)

# ══════════════════════════════════════════════════════════════════════
# 2. 풍속 구간 설정
#    U_min ~ U_max 를 N_bins 구간으로 분할
#    각 구간 중앙값을 대표 풍속으로 사용
# ══════════════════════════════════════════════════════════════════════
U_bins  = np.arange(0.5, 35.5, 1.0)   # [m/s] 구간 중앙값: 0.5, 1.5, ..., 34.5
dU      = 1.0                           # [m/s] 구간 폭
T_year  = 365.25 * 24 * 3600           # [s]   1년

# 수목 파라미터 (기준 조건)
H_TREE  = 8.0
DBH     = 0.20
Z0      = 0.3
I_U     = 0.20   # 난류강도 (도시 대표값)

# S-N 파라미터 (은행나무)
MOR, A_SN, B_SN = 68.0, 0.97, 0.072

# ══════════════════════════════════════════════════════════════════════
# 3. 공통 함수 (param_study.py에서 재사용)
# ══════════════════════════════════════════════════════════════════════
def kaimal_psd(f, U_bar, z, z0, I_u):
    alpha   = 0.67 + 0.05 * np.log(z0)
    L_u     = 300.0 * (z / 200.0) ** alpha
    sigma_u = I_u * U_bar
    f_L     = f * L_u / U_bar
    return sigma_u**2 * 6.8 * f_L / (f * (1.0 + 10.2 * f_L)**(5.0/3.0))

def generate_turbulence(freq, psd, N, I_u, U_bar, seed=42):
    rng   = np.random.default_rng(seed)
    df    = freq[1] - freq[0]
    amp   = np.sqrt(2.0 * psd * df)
    phase = rng.uniform(0, 2*np.pi, len(freq))
    spec  = amp * np.exp(1j * phase)
    spec[0] = 0.0
    u_t   = np.fft.irfft(spec, n=N)
    u_t  *= (I_u * U_bar) / np.std(u_t)
    return u_t

def stress_at_wind(U_bar, H_tree=H_TREE, DBH=DBH, z0=Z0, I_u=I_U,
                   T_total=600.0, dt=0.05):
    """주어진 U_bar에서 600초 응력 시계열 생성"""
    N     = int(T_total / dt)
    freq  = np.fft.rfftfreq(N, d=dt)
    freq[0] = freq[1] * 1e-6
    z     = H_tree * 0.5

    psd   = kaimal_psd(freq, U_bar, z, z0, I_u)
    u_t   = generate_turbulence(freq, psd, N, I_u, U_bar)

    U_t   = np.maximum(U_bar + u_t, 0.0)
    A_crown = (H_tree * 0.4)**2 * np.pi / 4
    F_t   = 0.5 * 1.225 * 0.5 * A_crown * U_t**2

    D_base = DBH * 1.2
    I_sec  = np.pi * D_base**4 / 64.0
    c      = D_base / 2.0
    L_arm  = H_tree * (2.0 / 3.0)
    return F_t * L_arm * c / I_sec * 1e-6   # [MPa]

def damage_per_600s(sigma_t, MOR, a, b):
    """10분 시뮬레이션에서의 손상 D"""
    D = 0.0
    for rng, *_ in rainflow.extract_cycles(sigma_t):
        amp = rng / 2.0
        if amp <= 0:
            continue
        ratio = amp / MOR
        Nf = 1.0 if ratio >= a else 10.0**((a - ratio) / b)
        D += 1.0 / Nf
    return D

# ══════════════════════════════════════════════════════════════════════
# 4. 풍속 구간별 손상 계산 (캐시 활용)
#    U < 5 m/s: 피로 무시 (응력 < 1 MPa)
#    U > 30 m/s: 출현 빈도 극히 낮음 → 개별 계산
# ══════════════════════════════════════════════════════════════════════
U_THRESHOLD = 5.0   # [m/s] 이하는 피로 무시

print("풍속 구간별 손상 계산 중 (U=5~35 m/s)...")
damage_cache = {}   # U_bin → D_600s
for U in U_bins:
    if U < U_THRESHOLD:
        damage_cache[U] = 0.0
        continue
    sig = stress_at_wind(U)
    damage_cache[U] = damage_per_600s(sig, MOR, A_SN, B_SN)
    print(f"  U={U:.1f} m/s  D_600s={damage_cache[U]:.3e}", end="\r")
print("\n  계산 완료.")

# ══════════════════════════════════════════════════════════════════════
# 5. Weibull 가중 연간 누적손상
#    D_annual = Σ_i  D_600s(U_i) * f(U_i)*dU * T_year / 600
# ══════════════════════════════════════════════════════════════════════
T_SIM = 600.0   # [s]

results = {}
for site, params in WEIBULL_SITES.items():
    k, c = params["k"], params["c"]
    U_mean = weibull_mean(k, c)
    D_annual = 0.0
    contrib  = []   # 구간별 기여도 (시각화용)
    for U in U_bins:
        p_U   = weibull_pdf(U, k, c) * dU   # 해당 구간 출현 확률
        T_U   = p_U * T_year                 # 연간 해당 풍속 노출 시간 [s]
        n_sim = T_U / T_SIM                  # 10분 시뮬레이션 횟수
        D_U   = damage_cache[U] * n_sim      # 해당 풍속의 연간 손상
        D_annual += D_U
        contrib.append((U, p_U, T_U/3600, D_U))

    L_yr = 1.0 / D_annual if D_annual > 0 else np.inf
    results[site] = {
        "k": k, "c": c, "U_mean": U_mean,
        "D_annual": D_annual, "L_yr": L_yr,
        "contrib": contrib,
        "color": params["color"]
    }
    print(f"  [{site}]  U_mean={U_mean:.1f}m/s  D_annual={D_annual:.4e}  수명={L_yr:.1f} 년")

# ══════════════════════════════════════════════════════════════════════
# 6. 결과 저장
# ══════════════════════════════════════════════════════════════════════
csv_path = out_dir / "weibull_fatigue_results.csv"
with open(csv_path, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["site", "k", "c", "U_mean_ms",
                     "D_annual", "fatigue_life_yr"])
    for site, r in results.items():
        writer.writerow([site, r['k'], r['c'], f"{r['U_mean']:.2f}",
                         f"{r['D_annual']:.4e}", f"{r['L_yr']:.2f}"])
print(f"  결과 저장: {csv_path}")

# ══════════════════════════════════════════════════════════════════════
# 7. 시각화
# ══════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(2, 3, figsize=(16, 10))
fig.suptitle(
    f"Weibull 풍속 분포 기반 장기 피로수명 평가\n"
    f"H={H_TREE}m, DBH={DBH*100:.0f}cm, S-N: 은행나무(MOR={MOR}MPa), Kaimal 스펙트럼",
    fontsize=11)

U_plot = np.linspace(0, 35, 300)

# ── 행 0: Weibull PDF, 손상 기여도, 누적손상 ──────────────────────
# 7-a. Weibull PDF 비교
ax = axes[0, 0]
for site, r in results.items():
    pdf = weibull_pdf(U_plot, r['k'], r['c'])
    ax.plot(U_plot, pdf, lw=2, color=r['color'],
            label=f"{site}\n(k={r['k']}, c={r['c']}, U_m={r['U_mean']:.1f}m/s)")
ax.set_xlabel("풍속 U [m/s]")
ax.set_ylabel("확률밀도 f(U)")
ax.set_title("지역별 Weibull 풍속 분포")
ax.legend(fontsize=6.5)
ax.grid(True, alpha=0.3)

# 7-b. 풍속 구간별 손상 기여도 (절대값)
ax = axes[0, 1]
for site, r in results.items():
    Us     = [c[0] for c in r['contrib']]
    D_Uns  = [c[3] for c in r['contrib']]
    ax.bar(Us, D_Uns, width=0.85, color=r['color'],
           alpha=0.5, label=site.split("(")[0])
ax.set_xlabel("풍속 U [m/s]")
ax.set_ylabel("연간 손상 기여 D(U)")
ax.set_title("풍속 구간별 연간 손상 기여도")
ax.legend(fontsize=7)
ax.grid(True, axis='y', alpha=0.3)

# 7-c. 연간 누적손상 → 피로수명 막대 비교
ax = axes[0, 2]
sites_short = [s.split("(")[0].strip() for s in results.keys()]
L_vals = [min(r['L_yr'], 1e6) for r in results.values()]
colors = [r['color'] for r in results.values()]
bars = ax.bar(sites_short, L_vals, color=colors, alpha=0.8, width=0.5)
for bar, r in zip(bars, results.values()):
    txt = f"{r['L_yr']:.0f} 년" if r['L_yr'] < 1e6 else ">1M 년"
    ax.text(bar.get_x() + bar.get_width()/2,
            bar.get_height() * 1.02, txt,
            ha='center', va='bottom', fontsize=8)
ax.set_ylabel("피로수명 [년]")
ax.set_title("지역별 피로수명 비교")
ax.set_yscale('log')
ax.grid(True, axis='y', alpha=0.3)
ax.tick_params(axis='x', labelsize=7)

# ── 행 1: DBH 파라미터 스터디 (Weibull 통합) ──────────────────────
print("\nDBH 파라미터 스터디 (Weibull 통합) 계산 중...")
DBH_range = np.array([0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40])

# 7-d. 수고 영향 (지역별 선)
ax = axes[1, 0]
H_range = np.array([4, 6, 8, 10, 12])
for site, r in results.items():
    k, c = r['k'], r['c']
    L_Hs = []
    for H in H_range:
        # DBH=20cm 고정, 수고 변화
        cache_H = {}
        for U in U_bins:
            if U < U_THRESHOLD:
                cache_H[U] = 0.0
                continue
            sig = stress_at_wind(U, H_tree=H, DBH=0.20)
            cache_H[U] = damage_per_600s(sig, MOR, A_SN, B_SN)
        D_ann = sum(damage_cache_H * weibull_pdf(U, k, c) * dU * T_year / T_SIM
                    for U, damage_cache_H in cache_H.items())
        L_Hs.append(min(1/D_ann if D_ann > 0 else 1e8, 1e8))
    ax.semilogy(H_range, L_Hs, lw=2, marker='o', color=r['color'],
                label=site.split("(")[0].strip())
ax.set_xlabel("수고 H [m]")
ax.set_ylabel("피로수명 [년]")
ax.set_title("수고 H 영향 (DBH=20cm)")
ax.legend(fontsize=7)
ax.grid(True, which='both', alpha=0.3)

# 7-e. DBH 영향 (지역별 선)
ax = axes[1, 1]
for site, r in results.items():
    k, c = r['k'], r['c']
    L_Ds = []
    for D in DBH_range:
        cache_D = {}
        for U in U_bins:
            if U < U_THRESHOLD:
                cache_D[U] = 0.0
                continue
            sig = stress_at_wind(U, H_tree=8.0, DBH=D)
            cache_D[U] = damage_per_600s(sig, MOR, A_SN, B_SN)
        D_ann = sum(dmg * weibull_pdf(U, k, c) * dU * T_year / T_SIM
                    for U, dmg in cache_D.items())
        L_Ds.append(min(1/D_ann if D_ann > 0 else 1e8, 1e8))
    ax.semilogy(DBH_range*100, L_Ds, lw=2, marker='o', color=r['color'],
                label=site.split("(")[0].strip())
ax.set_xlabel("흉고직경 DBH [cm]")
ax.set_ylabel("피로수명 [년]")
ax.set_title("DBH 영향 (H=8m)")
ax.legend(fontsize=7)
ax.grid(True, which='both', alpha=0.3)

# 7-f. 풍속-손상 관계 (10분 D vs U)
ax = axes[1, 2]
U_plot_bins = [U for U in U_bins if U >= U_THRESHOLD]
D_plot      = [damage_cache[U] for U in U_plot_bins]
ax.semilogy(U_plot_bins, D_plot, lw=2, color='black',
            marker='.', markersize=4, label="D/10min (Kaimal)")
ax.set_xlabel("풍속 U [m/s]")
ax.set_ylabel("10분 손상 D_600s")
ax.set_title("풍속별 10분 손상량\n(H=8m, DBH=20cm, 은행나무)")
ax.grid(True, which='both', alpha=0.3)
ax.legend(fontsize=8)

plt.tight_layout()
fig_path = out_dir / "weibull_fatigue.png"
plt.savefig(fig_path, dpi=150, bbox_inches='tight')
plt.show()
print(f"\n  그래프 저장: {fig_path}")

# ══════════════════════════════════════════════════════════════════════
# 최종 요약 출력
# ══════════════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("  Weibull 통합 피로수명 최종 요약")
print("  (은행나무, H=8m, DBH=20cm)")
print("=" * 60)
for site, r in results.items():
    print(f"  {site}")
    print(f"    U_mean={r['U_mean']:.1f}m/s  D_annual={r['D_annual']:.3e}  "
          f"수명={r['L_yr']:.1f} 년")
