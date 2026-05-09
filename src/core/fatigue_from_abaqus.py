"""
fatigue_from_abaqus.py
----------------------
ABAQUS 동적 해석 결과(응력 시계열) → Rainflow + Miner's rule 피로수명 평가
extract_odb.py 실행 후 생성된 CSV를 입력으로 사용

사용 흐름:
    1. abaqus python extract_odb.py        (ABAQUS 환경)
    2. python fatigue_from_abaqus.py       (일반 Python 환경)

또는 방법 A (CAE 수동 내보내기) 사용 시:
    CSV 컬럼: time_s, sigma_bending_MPa  (최소 2열)
    파일 경로: output/abaqus_stress_history.csv
"""

import sys
import os

# 어느 폴더에서 실행하든 프로젝트 모듈 탐색
PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
if PROJECT_DIR not in sys.path:
    sys.path.insert(0, PROJECT_DIR)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import rainflow
import csv
from pathlib import Path
from wood_properties_db import STATIC_PROPS, SN_PARAMS

matplotlib.rcParams['font.family'] = 'Malgun Gothic'
matplotlib.rcParams['axes.unicode_minus'] = False

out_dir = Path("output")
out_dir.mkdir(exist_ok=True)

# ══════════════════════════════════════════════════════════════════════
# 1. ABAQUS 결과 CSV 로드
#    방법 A (CAE 수동 내보내기): 컬럼명이 다를 수 있으므로 자동 탐지
# ══════════════════════════════════════════════════════════════════════

def load_abaqus_csv(csv_path):
    """
    ABAQUS History Output CSV 로드
    extract_odb.py 출력 또는 CAE 수동 내보내기 모두 지원
    """
    t_arr, sigma_arr = [], []

    with open(csv_path, 'r', encoding='utf-8-sig') as f:
        lines = f.readlines()

    # 헤더 탐색 (# 주석 및 빈 줄 건너뜀)
    header_idx = 0
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith('#') or not stripped:
            continue
        # 숫자로 시작하면 데이터 행
        try:
            float(stripped.split(',')[0])
            header_idx = i
            break
        except ValueError:
            # 헤더 행 — 컬럼 위치 탐색
            cols = [c.strip().lower() for c in stripped.split(',')]
            t_col     = next((i for i, c in enumerate(cols)
                              if 'time' in c), 0)
            sigma_col = next((i for i, c in enumerate(cols)
                              if 'sigma' in c or 'stress' in c
                              or 'sf3' in c or 'bending' in c), -1)
            header_idx = i + 1
            break
    else:
        t_col, sigma_col = 0, -1
        header_idx = 0

    # 데이터 파싱
    for line in lines[header_idx:]:
        stripped = line.strip()
        if not stripped or stripped.startswith('#'):
            continue
        parts = stripped.split(',')
        try:
            t = float(parts[0])
            # sigma 컬럼 자동 선택: 마지막 열 (CAE 내보내기 형식 대응)
            if sigma_col == -1 or sigma_col >= len(parts):
                sigma_col = len(parts) - 1
            sigma = abs(float(parts[sigma_col]))
            t_arr.append(t)
            sigma_arr.append(sigma)
        except (ValueError, IndexError):
            continue

    return np.array(t_arr), np.array(sigma_arr)


csv_abaqus = out_dir / 'abaqus_stress_history.csv'

if not csv_abaqus.exists():
    print(f"[경고] {csv_abaqus} 없음")
    print("  방법 A: ABAQUS CAE → Result → History Output → SF3 → Save As CSV")
    print("  방법 B: abaqus python extract_odb.py 실행 후 재시도")
    print("\n  데모용 임시 데이터로 실행합니다...")

    # 데모: Python 단독 결과를 ABAQUS 결과로 가정
    from wind_time_history import *   # noqa
    import importlib, sys
    # fatigue_analysis.py 의 stress 계산 재사용
    from fatigue_analysis import force_to_stress, load_drag_force
    t_demo, F_demo = load_drag_force(out_dir / 'wind_kaimal.csv')
    D_BASE = 0.24
    I_BASE = np.pi * D_BASE**4 / 64
    C_BASE = D_BASE / 2
    L_ARM  = 8.0 * (2/3)
    sigma_demo = force_to_stress(F_demo, L_ARM, C_BASE, I_BASE)

    # DAF 1.0 (Python 단독) — 비교 기준
    t_abaqus     = t_demo
    sigma_abaqus = sigma_demo
    DAF_note     = "DAF=1.0 (데모, ABAQUS 결과 아님)"
else:
    t_abaqus, sigma_abaqus = load_abaqus_csv(csv_abaqus)
    DAF_note = "ABAQUS 동적 해석 결과"

print(f"\n  [{DAF_note}]")
print(f"  시간 범위: {t_abaqus[0]:.1f} ~ {t_abaqus[-1]:.1f} s  ({len(t_abaqus)} 스텝)")
print(f"  σ 최대: {sigma_abaqus.max():.3f} MPa")
print(f"  σ 평균: {sigma_abaqus.mean():.3f} MPa")
print(f"  σ std : {sigma_abaqus.std():.3f} MPa")

# ══════════════════════════════════════════════════════════════════════
# 2. Python 단독 결과 로드 (비교용)
# ══════════════════════════════════════════════════════════════════════
from fatigue_analysis import load_drag_force, force_to_stress

t_py, F_py = load_drag_force(out_dir / 'wind_kaimal.csv')
D_BASE = 0.24
I_BASE = np.pi * D_BASE**4 / 64
C_BASE = D_BASE / 2
L_ARM  = 8.0 * (2/3)
sigma_py = force_to_stress(F_py, L_ARM, C_BASE, I_BASE)

# 동적 증폭 계수 (DAF)
DAF = sigma_abaqus.max() / sigma_py.max() if sigma_py.max() > 0 else 1.0
print(f"\n  동적 증폭 계수 DAF = {DAF:.3f}")
print(f"  (DAF>1이면 ABAQUS가 공진 증폭 포착)")

# ══════════════════════════════════════════════════════════════════════
# 3. Rainflow counting + Miner's rule
# ══════════════════════════════════════════════════════════════════════
def run_fatigue(sigma_t, label):
    cycles = list(rainflow.extract_cycles(sigma_t))
    ranges  = np.array([c[0] for c in cycles])
    counts  = np.array([c[2] for c in cycles])
    amps    = ranges / 2.0
    print(f"\n  [{label}] Rainflow 결과")
    print(f"    총 사이클: {counts.sum():.0f}")
    print(f"    진폭 최대: {amps.max():.4f} MPa  /  평균: {np.average(amps,weights=counts):.4f} MPa")
    return cycles

cycles_abaqus = run_fatigue(sigma_abaqus, DAF_note)
cycles_py     = run_fatigue(sigma_py,     "Python 단독 (Kaimal)")

def miners_damage(cycles, MOR, a, b, T_sim=600.0, T_year=100*3600):
    D = 0.0
    for rng, *_ in cycles:
        amp = rng / 2.0
        if amp <= 0: continue
        ratio = amp / MOR
        Nf = 1.0 if ratio >= a else 10.0**((a - ratio) / b)
        D += 1.0 / Nf
    D_ann = D * (T_year / T_sim)
    return D, D_ann, (1/D_ann if D_ann > 0 else np.inf)

# 은행나무 기준
sp    = "은행나무 (Ginkgo biloba)"
MOR_g = STATIC_PROPS[sp]["MOR_MPa"] * 0.60
a_sn  = SN_PARAMS[sp]["a_R0"]
b_sn  = SN_PARAMS[sp]["b_R0"]

D_sim_abq, D_ann_abq, L_abq = miners_damage(cycles_abaqus, MOR_g, a_sn, b_sn)
D_sim_py,  D_ann_py,  L_py  = miners_damage(cycles_py,     MOR_g, a_sn, b_sn)

print("\n" + "="*55)
print(f"  피로수명 비교 ({sp.split('(')[0].strip()}, 생재, 연 100h 강풍)")
print("="*55)
print(f"  {'':20}  {'D/10min':>12}  {'D/year':>12}  {'수명(년)':>12}")
print(f"  {'Python 단독':20}  {D_sim_py:>12.3e}  {D_ann_py:>12.4f}  {min(L_py,1e8):>12.0f}")
print(f"  {'ABAQUS 동적':20}  {D_sim_abq:>12.3e}  {D_ann_abq:>12.4f}  {min(L_abq,1e8):>12.0f}")
if L_py < 1e8 and L_abq < 1e8:
    print(f"  수명 감소율 (DAF 효과): {(L_py - L_abq)/L_py*100:.1f}%")
print("="*55)

# ══════════════════════════════════════════════════════════════════════
# 4. 시각화
# ══════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(2, 2, figsize=(14, 9))
fig.suptitle(f"ABAQUS 동적 해석 결과 vs Python 단독 비교\n"
             f"H=8m, DBH=20cm, Kaimal, 은행나무 생재(MOR={MOR_g:.0f}MPa)",
             fontsize=11)

mask = t_abaqus <= 60

# ── 응력 시계열 비교 ──────────────────────────────────────────────
ax = axes[0, 0]
ax.plot(t_abaqus[mask], sigma_abaqus[mask], lw=0.8, color='crimson',
        label=f'ABAQUS 동적 (max={sigma_abaqus.max():.2f}MPa)')
mask_py = t_py <= 60
ax.plot(t_py[mask_py], sigma_py[mask_py], lw=0.8, color='steelblue',
        alpha=0.7, label=f'Python 단독 (max={sigma_py.max():.2f}MPa)')
ax.set_xlabel('시간 [s]')
ax.set_ylabel('근원부 굽힘 응력 [MPa]')
ax.set_title(f'응력 시계열 비교 (처음 60s)  /  DAF={DAF:.3f}')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# ── Rainflow 진폭 분포 비교 ───────────────────────────────────────
ax = axes[0, 1]
amps_abq = np.array([c[0]/2 for c in cycles_abaqus])
cnt_abq  = np.array([c[2]   for c in cycles_abaqus])
amps_py  = np.array([c[0]/2 for c in cycles_py])
cnt_py   = np.array([c[2]   for c in cycles_py])

bins = np.linspace(0, max(amps_abq.max(), amps_py.max()) * 1.05, 35)
for amps, cnts, color, label in [
    (amps_abq, cnt_abq, 'crimson',   'ABAQUS 동적'),
    (amps_py,  cnt_py,  'steelblue', 'Python 단독'),
]:
    hist, edges = np.histogram(amps, bins=bins, weights=cnts)
    centers = (edges[:-1] + edges[1:]) / 2
    ax.bar(centers, hist, width=np.diff(edges),
           color=color, alpha=0.5, label=label)
ax.set_xlabel('응력 진폭 [MPa]')
ax.set_ylabel('사이클 수')
ax.set_title('Rainflow 사이클 분포 비교')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# ── PSD 비교 (응력 시계열) ────────────────────────────────────────
ax = axes[1, 0]
from scipy.signal import welch
dt = float(t_abaqus[1] - t_abaqus[0])
f_w, psd_abq = welch(sigma_abaqus, fs=1/dt, nperseg=len(sigma_abaqus)//8)
f_w, psd_py  = welch(sigma_py,     fs=1/dt, nperseg=len(sigma_py)//8)
ax.loglog(f_w[1:], psd_abq[1:], lw=1.2, color='crimson',   label='ABAQUS 동적')
ax.loglog(f_w[1:], psd_py[1:],  lw=1.2, color='steelblue', label='Python 단독')
ax.axvspan(0.2, 0.8, alpha=0.12, color='green', label='f₁ 공진 대역')
ax.axvline(0.582, color='red', lw=1.5, ls='--', label='f₁=0.582Hz (ABAQUS)')
ax.set_xlabel('주파수 [Hz]')
ax.set_ylabel('응력 PSD [MPa²/Hz]')
ax.set_title('응력 PSD 비교 (공진 대역 에너지 증폭 확인)')
ax.legend(fontsize=7)
ax.grid(True, which='both', alpha=0.3)

# ── 피로수명 비교 막대 ────────────────────────────────────────────
ax = axes[1, 1]
labels_bar = ['Python\n단독', 'ABAQUS\n동적']
lives_bar  = [min(L_py, 1e8), min(L_abq, 1e8)]
colors_bar = ['steelblue', 'crimson']
bars = ax.bar(labels_bar, lives_bar, color=colors_bar, alpha=0.8, width=0.4)
for bar, L in zip(bars, [L_py, L_abq]):
    txt = f'{L:.0f}년' if L < 1e8 else '>1억년'
    ax.text(bar.get_x() + bar.get_width()/2,
            bar.get_height() * 1.05, txt,
            ha='center', va='bottom', fontsize=9)
ax.set_yscale('log')
ax.set_ylabel('피로수명 [년]')
ax.set_title(f'피로수명 비교\n(DAF={DAF:.3f}, 은행나무 생재)')
ax.grid(True, axis='y', alpha=0.3)

plt.tight_layout()
fig_path = out_dir / 'fatigue_abaqus_comparison.png'
plt.savefig(fig_path, dpi=150, bbox_inches='tight')
plt.show()
print(f'\n  그래프 저장: {fig_path}')
