"""
wind_time_history.py
--------------------
Davenport(1961) 및 Kaimal(1972) 스펙트럼 기반 풍속 시간이력 생성 및 비교
- 평균풍속 U_bar에 난류 성분 u(t)를 중첩
- 항력 공식: F(t) = 0.5 * rho * Cd * A * U(t)^2
- ABAQUS *AMPLITUDE 입력용 파일 출력 (두 모델 각각)

참고 문헌:
  Davenport, A.G. (1961). The spectrum of horizontal gustiness near the
    ground in high winds. Q.J.R. Meteorol. Soc., 87, 194-211.
  Kaimal, J.C. et al. (1972). Spectral characteristics of surface-layer
    turbulence. Q.J.R. Meteorol. Soc., 98, 563-589.
  EN 1991-1-4 (Eurocode 1, 2005). Wind actions — Annex B (Kaimal form).
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import csv
from pathlib import Path
from scipy.signal import welch

# Windows 한글 폰트 설정
matplotlib.rcParams['font.family'] = 'Malgun Gothic'
matplotlib.rcParams['axes.unicode_minus'] = False

# ══════════════════════════════════════════════════════════════════════
# 1. 해석 파라미터
# ══════════════════════════════════════════════════════════════════════
T_total = 600.0        # [s]   총 시뮬레이션 시간 (기상 표준 10분)
dt      = 0.05         # [s]   시간 간격  → Nyquist f_max = 10 Hz
N       = int(T_total / dt)

# 풍환경 파라미터
U_bar = 20.0           # [m/s]  기준 높이 평균풍속 (설계풍속)
I_u   = 0.20           # [-]    난류강도 (도심 0.2~0.3 / 외곽 0.1~0.15)
z     = 6.0            # [m]    수목 유효높이 (수고 절반 기준)
z0    = 0.3            # [m]    지표 조도 길이 (도시 지표)

# 수목 물성 (항력 계산용)
rho_air  = 1.225       # [kg/m³] 공기 밀도
Cd       = 0.5         # [-]     항력계수 (수관 형상에 따라 0.3~0.8)
A_crown  = 12.0        # [m²]   수관 투영 면적

# ══════════════════════════════════════════════════════════════════════
# 2. 스펙트럼 함수 정의
# ══════════════════════════════════════════════════════════════════════

def davenport_psd(f, U_bar, L=1200.0):
    """
    Davenport (1961) 종방향 풍속 PSD

        S_u(f) = 4·K·U_bar²·x² / [f·(1 + x²)^(4/3)]
        x = f·L / U_bar   (무차원 주파수)
        K = 지표면 항력계수 (도시 ≈ 0.005)
        L = 길이스케일 [m] (Davenport 고정값 1200 m)

    단위: [m²/s²/Hz]
    """
    K = 0.005
    x = f * L / U_bar
    return 4.0 * K * U_bar**2 * x**2 / (f * (1.0 + x**2)**(4.0 / 3.0))


def kaimal_psd(f, U_bar, z, z0, I_u):
    """
    Kaimal et al. (1972) / Eurocode 1 (EN 1991-1-4 Annex B) 종방향 풍속 PSD

        S_u(f) = σ_u²·6.8·f_L / [f·(1 + 10.2·f_L)^(5/3)]
        f_L    = f·L_u(z) / U_bar   (무차원 주파수)
        σ_u    = I_u · U_bar         (난류 표준편차)
        L_u(z) = 300·(z/200)^α      (Eurocode 1 길이스케일)
        α      = 0.67 + 0.05·ln(z0)

    단위: [m²/s²/Hz]
    참고: 고주파 영역 -5/3 지수법칙(Kolmogorov) 엄밀 반영
    """
    alpha = 0.67 + 0.05 * np.log(z0)          # 조도 의존 지수
    L_u   = 300.0 * (z / 200.0) ** alpha      # [m] 적분 길이스케일
    sigma_u = I_u * U_bar                      # [m/s] 난류 표준편차
    f_L   = f * L_u / U_bar                   # 무차원 주파수
    return sigma_u**2 * 6.8 * f_L / (f * (1.0 + 10.2 * f_L)**(5.0 / 3.0))


# ══════════════════════════════════════════════════════════════════════
# 3. 주파수 배열
# ══════════════════════════════════════════════════════════════════════
freq     = np.fft.rfftfreq(N, d=dt)    # [Hz], 0 ~ f_Nyquist
freq[0]  = freq[1] * 1e-6             # f=0 특이점 방지

# ══════════════════════════════════════════════════════════════════════
# 4. 스펙트럼 → 난류 시간이력 (Shinozuka & Jan, 1972)
# ══════════════════════════════════════════════════════════════════════
def generate_turbulence(freq, psd, N, dt, I_u, U_bar, seed=42):
    """
    랜덤 위상각 부여 후 IFFT로 난류 시간이력 복원.
    σ_u = I_u·U_bar 로 표준편차 스케일링.
    """
    rng     = np.random.default_rng(seed)
    df      = freq[1] - freq[0]
    amp     = np.sqrt(2.0 * psd * df)
    phase   = rng.uniform(0, 2 * np.pi, len(freq))

    spectrum    = amp * np.exp(1j * phase)
    spectrum[0] = 0.0                          # DC 제거

    u_t = np.fft.irfft(spectrum, n=N)

    sigma_target = I_u * U_bar
    u_t *= sigma_target / np.std(u_t)         # σ 정규화
    return u_t


# 같은 seed → 위상 구조 동일, 스펙트럼 형태만 다름 (공정 비교)
psd_dav  = davenport_psd(freq, U_bar)
psd_kai  = kaimal_psd(freq, U_bar, z, z0, I_u)

u_dav    = generate_turbulence(freq, psd_dav, N, dt, I_u, U_bar, seed=42)
u_kai    = generate_turbulence(freq, psd_kai, N, dt, I_u, U_bar, seed=42)

# ══════════════════════════════════════════════════════════════════════
# 5. 풍속 및 항력 계산
# ══════════════════════════════════════════════════════════════════════
t     = np.arange(N) * dt

U_dav = np.maximum(U_bar + u_dav, 0.0)
U_kai = np.maximum(U_bar + u_kai, 0.0)

F_dav = 0.5 * rho_air * Cd * A_crown * U_dav**2
F_kai = 0.5 * rho_air * Cd * A_crown * U_kai**2

# ══════════════════════════════════════════════════════════════════════
# 6. 통계 비교 출력
# ══════════════════════════════════════════════════════════════════════
def stats(label, u_t, U_t, F_t, U_bar, I_u):
    print(f"\n  [ {label} ]")
    print(f"  σ_u 목표 / 실제  = {I_u*U_bar:.2f} / {np.std(u_t):.2f} m/s")
    print(f"  최대 풍속        = {U_t.max():.2f} m/s")
    print(f"  항력 평균 / 최대 = {F_t.mean():.1f} / {F_t.max():.1f} N")

print("=" * 55)
print("  풍속 시간이력 생성 결과 비교")
print(f"  U_bar={U_bar} m/s  I_u={I_u}  z={z} m  T={T_total:.0f} s")
print("=" * 55)
stats("Davenport (1961)",  u_dav, U_dav, F_dav, U_bar, I_u)
stats("Kaimal (1972) / Eurocode 1", u_kai, U_kai, F_kai, U_bar, I_u)
print("=" * 55)

# Kaimal 길이스케일 정보 출력
alpha_val = 0.67 + 0.05 * np.log(z0)
L_u_val   = 300.0 * (z / 200.0) ** alpha_val
print(f"\n  Kaimal 길이스케일 L_u(z={z}m, z0={z0}m) = {L_u_val:.1f} m")
print(f"  Davenport 길이스케일 L                  = 1200 m")
print(f"  → L_u 차이가 저주파 에너지 분포 차이를 유발\n")

# ══════════════════════════════════════════════════════════════════════
# 7. 파일 저장
# ══════════════════════════════════════════════════════════════════════
out_dir = Path("output")
out_dir.mkdir(exist_ok=True)


def save_csv(path, t, U_t, u_t, F_t, label):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([f"# {label}"])
        writer.writerow(["time_s", "U_ms", "u_turbulent_ms", "F_drag_N"])
        for i in range(len(t)):
            writer.writerow([f"{t[i]:.4f}", f"{U_t[i]:.4f}",
                             f"{u_t[i]:.4f}", f"{F_t[i]:.4f}"])
    print(f"  CSV 저장: {path}")


def save_abaqus_amplitude(path, t, F_t, label):
    F_ref = F_t.max()
    amp   = F_t / F_ref
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"** {label}\n")
        f.write(f"** Fref = {F_ref:.2f} N  (ABAQUS Cload 에 이 값 입력)\n")
        f.write("** Generated by wind_time_history.py\n")
        name = "WIND_DAV" if "Davenport" in label else "WIND_KAI"
        f.write(f"*Amplitude, name={name}, definition=TABULAR\n")
        for i in range(0, len(t), 4):
            pairs = [f"{t[i+j]:.4f}, {amp[i+j]:.6f}"
                     for j in range(4) if i + j < len(t)]
            f.write(", ".join(pairs) + "\n")
    print(f"  ABAQUS amplitude 저장: {path}  (Fref={F_ref:.1f} N)")
    return F_ref


save_csv(out_dir / "wind_davenport.csv", t, U_dav, u_dav, F_dav, "Davenport 1961")
save_csv(out_dir / "wind_kaimal.csv",    t, U_kai, u_kai, F_kai, "Kaimal 1972 / Eurocode 1")

Fref_dav = save_abaqus_amplitude(out_dir / "amplitude_davenport.inp", t, F_dav, "Davenport (1961)")
Fref_kai = save_abaqus_amplitude(out_dir / "amplitude_kaimal.inp",    t, F_kai, "Kaimal (1972) / Eurocode 1")

# ══════════════════════════════════════════════════════════════════════
# 8. 시각화
# ══════════════════════════════════════════════════════════════════════
mask      = t <= 60                        # 처음 60초 표시
f_theory  = np.logspace(-3, 1, 500)
psd_dav_th = davenport_psd(f_theory, U_bar)
psd_kai_th = kaimal_psd(f_theory, U_bar, z, z0, I_u)

_, psd_dav_w = welch(u_dav, fs=1/dt, nperseg=N//8)
f_w, psd_kai_w = welch(u_kai, fs=1/dt, nperseg=N//8)

fig, axes = plt.subplots(3, 1, figsize=(13, 10))
fig.suptitle(
    f"Davenport vs. Kaimal 풍속 시간이력 비교\n"
    f"U_bar={U_bar} m/s, I_u={I_u}, z={z} m, z0={z0} m",
    fontsize=11)

# ── 8-a. 풍속 시간이력 비교 ──────────────────────────────────────────
ax = axes[0]
ax.plot(t[mask], U_dav[mask], lw=0.7, color="steelblue",  label="Davenport")
ax.plot(t[mask], U_kai[mask], lw=0.7, color="darkorange", label="Kaimal", alpha=0.8)
ax.axhline(U_bar, color="black", lw=1.0, ls="--", label=f"U_bar={U_bar} m/s")
ax.set_ylabel("풍속 U(t) [m/s]")
ax.set_xlabel("시간 [s]")
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# ── 8-b. 항력 시간이력 비교 ──────────────────────────────────────────
ax = axes[1]
ax.plot(t[mask], F_dav[mask], lw=0.7, color="steelblue",  label=f"Davenport  (Fref={Fref_dav:.0f} N)")
ax.plot(t[mask], F_kai[mask], lw=0.7, color="darkorange", label=f"Kaimal     (Fref={Fref_kai:.0f} N)", alpha=0.8)
ax.set_ylabel("항력 F(t) [N]")
ax.set_xlabel("시간 [s]")
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# ── 8-c. PSD 비교 (이론 + Welch 검증) ───────────────────────────────
ax = axes[2]
ax.loglog(f_theory, psd_dav_th, lw=2.0, ls="--", color="steelblue",
          label="Davenport 이론")
ax.loglog(f_theory, psd_kai_th, lw=2.0, ls="--", color="darkorange",
          label="Kaimal 이론 (Eurocode 1)")
ax.loglog(f_w[1:], psd_dav_w[1:], lw=0.8, color="steelblue",  alpha=0.6,
          label="Davenport 생성 신호 (Welch)")
ax.loglog(f_w[1:], psd_kai_w[1:], lw=0.8, color="darkorange", alpha=0.6,
          label="Kaimal 생성 신호 (Welch)")

# 수목 1차 고유진동수 범위 표시 (0.2 ~ 0.8 Hz)
ax.axvspan(0.2, 0.8, alpha=0.12, color="green", label="수목 고유진동수 범위 (0.2~0.8 Hz)")
ax.set_xlabel("주파수 [Hz]")
ax.set_ylabel("PSD [m²/s²/Hz]")
ax.legend(fontsize=7, ncol=2)
ax.grid(True, which="both", alpha=0.3)
ax.set_title("PSD 비교: Davenport vs. Kaimal — 녹색 영역이 수목 공진 대역", fontsize=9)

plt.tight_layout()
fig_path = out_dir / "wind_spectrum_comparison.png"
plt.savefig(fig_path, dpi=150)
plt.show()
print(f"\n  그래프 저장 완료: {fig_path}")
