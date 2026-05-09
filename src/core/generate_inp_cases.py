"""
generate_inp_cases.py
---------------------
6케이스 ABAQUS INP 파일 자동 생성
  케이스: H x DBH = {6, 8}m x {15, 20, 25}cm
  하중:   X방향(평균+난류) + Y방향(난류전용) Kaimal 시계열
  출력:   inp_cases/ 폴더에 tree_H{H}m_DBH{DBH}cm.inp 6개

사용법:
  python generate_inp_cases.py
  -> 생성된 *.inp 파일을 ABAQUS 실행 폴더로 복사 후
     abaqus job=tree_H8m_DBH20cm interactive
"""

import numpy as np
import os
from pathlib import Path

# ======================================================================
# 1. 케이스 정의
# ======================================================================
CASES = [
    {"H": 6.0, "DBH": 0.15},
    {"H": 6.0, "DBH": 0.20},
    {"H": 6.0, "DBH": 0.25},
    {"H": 8.0, "DBH": 0.15},
    {"H": 8.0, "DBH": 0.20},   # 기존 기준 케이스
    {"H": 8.0, "DBH": 0.25},
]

# ======================================================================
# 2. 공통 파라미터
# ======================================================================
# 목재 물성 (은행나무 생재)
E_WOOD   = 7.0e9 * 0.68       # [Pa]  생재 탄성계수 (4.76 GPa)
NU_WOOD  = 0.35                # [-]   포아송비
RHO_WOOD = 780.0               # [kg/m3] 생재 밀도

# 공기 / 항력
RHO_AIR  = 1.225
CD       = 0.5
# 수관 투영면적 기준: H=8m -> 12 m^2  (crown_diam ~ 0.49H)
# 스케일: A ~ H^2
A_CROWN_REF = 12.0
H_REF       = 8.0

# 풍속 / 난류
U_BAR   = 20.0                 # [m/s]  설계 기준풍속
I_U     = 0.20                 # [-]    종방향 난류강도
I_V_RAT = 0.75                 # [-]    횡/종 난류강도 비 (ESDU 83045)
Z0      = 0.3                  # [m]    지표 조도 길이

# 시뮬레이션 설정
DT      = 0.05                 # [s]
T_SIM   = 600.0                # [s]
N_STEPS = int(T_SIM / DT)      # = 12000

# 트렁크 테이퍼
# D_base = DBH x 1.2,  D_tip = DBH x 0.3
TAPER_BASE = 1.2
TAPER_TIP  = 0.3

# 수관 집중 질량 스케일 (H=8m, DBH=0.20m -> 150 kg)
M_CROWN_REF = 150.0

# 출력 폴더
out_dir = Path("inp_cases")
out_dir.mkdir(exist_ok=True)

# ======================================================================
# 3. 함수 정의
# ======================================================================

def kaimal_psd(f, U_bar, z, z0, sigma_u):
    alpha = 0.67 + 0.05 * np.log(z0)
    L_u   = 300.0 * (z / 200.0) ** alpha
    f_L   = f * L_u / U_bar
    return sigma_u**2 * 6.8 * f_L / (f * (1.0 + 10.2 * f_L)**(5.0/3.0))


def gen_turbulence(freq, psd, N, sigma_target, seed):
    """Shinozuka IFFT 난류 시계열 생성"""
    rng     = np.random.default_rng(seed)
    df      = freq[1] - freq[0]
    amp     = np.sqrt(2.0 * psd * df)
    phase   = rng.uniform(0, 2 * np.pi, len(freq))
    spectrum    = amp * np.exp(1j * phase)
    spectrum[0] = 0.0
    u_t = np.fft.irfft(spectrum, n=N)
    u_t *= sigma_target / (np.std(u_t) + 1e-12)
    return u_t


def generate_amplitudes(H, A_crown, seed_x=42, seed_y=123):
    """
    X방향 (평균+난류) 및 Y방향 (난류전용) Kaimal 진폭 배열 생성
    반환: t, amp_X, F_X_ref, amp_Y, F_Y_ref
    """
    freq    = np.fft.rfftfreq(N_STEPS, d=DT)
    freq[0] = freq[1] * 1e-6   # DC 특이점 방지

    z_eff   = H * 0.7           # 유효 높이 (수고의 70% ~ 수관 중심)
    sigma_u = I_U * U_BAR
    sigma_v = I_V_RAT * sigma_u  # 횡방향 난류

    psd_u = kaimal_psd(freq, U_BAR, z_eff, Z0, sigma_u)
    psd_v = kaimal_psd(freq, U_BAR, z_eff, Z0, sigma_v)

    u_t = gen_turbulence(freq, psd_u, N_STEPS, sigma_u, seed_x)
    v_t = gen_turbulence(freq, psd_v, N_STEPS, sigma_v, seed_y)

    # X방향: F_X(t) = 0.5*rho*Cd*A*(U_bar+u(t))^2
    U_t  = np.maximum(U_BAR + u_t, 0.0)
    F_Xt = 0.5 * RHO_AIR * CD * A_crown * U_t**2
    F_X_ref = float(F_Xt.max())
    amp_X   = F_Xt / F_X_ref

    # Y방향: F_Y(t) = rho*Cd*A*U_bar*v(t) [선형화, 평균=0]
    F_Yt = RHO_AIR * CD * A_crown * U_BAR * v_t
    F_Y_ref = float(np.max(np.abs(F_Yt)))
    if F_Y_ref > 0:
        amp_Y = F_Yt / F_Y_ref
    else:
        amp_Y = np.zeros_like(F_Yt)

    t_arr = np.arange(N_STEPS) * DT
    return t_arr, amp_X, F_X_ref, amp_Y, F_Y_ref


def amp_to_lines(t_arr, amp_arr, pairs_per_line=4):
    """진폭 배열 -> ABAQUS *Amplitude 텍스트 라인 리스트"""
    lines = []
    N = len(t_arr)
    for i in range(0, N, pairs_per_line):
        pairs = []
        for j in range(pairs_per_line):
            if i + j < N:
                pairs.append(f"{t_arr[i+j]:.4f}, {amp_arr[i+j]:+.6f}")
        lines.append(", ".join(pairs))
    return lines


def section_radius_at(D_base, D_tip, n_elem, i_elem):
    """i번째 요소(0-indexed) 중앙점 반지름 (선형 테이퍼)"""
    t = (i_elem + 0.5) / n_elem
    D = D_base + (D_tip - D_base) * t
    return D / 2.0


# ======================================================================
# 4. INP 파일 작성
# ======================================================================

def write_inp(case):
    H   = case["H"]
    DBH = case["DBH"]

    D_base = DBH * TAPER_BASE
    D_tip  = DBH * TAPER_TIP
    n_elem = 10
    n_node = n_elem + 1       # 11

    # 절점 Z 좌표 (Z축 방향, 수직)
    z_nodes = np.linspace(0, H, n_node)

    # 수관 투영면적 (H^2 스케일)
    A_crown = A_CROWN_REF * (H / H_REF)**2

    # 수관 집중 질량 (H, DBH 선형 스케일)
    M_crown = M_CROWN_REF * (H / H_REF) * (DBH / 0.20)

    # Kaimal 진폭 생성
    print(f"  [{H:.0f}m / DBH={DBH*100:.0f}cm] Kaimal 진폭 생성 중...")
    t_arr, amp_X, F_X_ref, amp_Y, F_Y_ref = generate_amplitudes(H, A_crown)

    amp_X_lines = amp_to_lines(t_arr, amp_X)
    amp_Y_lines = amp_to_lines(t_arr, amp_Y)

    fname = out_dir / f"tree_H{H:.0f}m_DBH{DBH*100:.0f}cm.inp"

    # ------------------------------------------------------------------
    lines = []
    def w(s=""):
        lines.append(s)

    w("*Heading")
    w(f" Tree Model: H={H:.1f}m, DBH={DBH*100:.0f}cm, Multi-axial Kaimal Wind")
    w(f" F_X_ref={F_X_ref:.2f} N  F_Y_ref={F_Y_ref:.2f} N")
    w(f" Generated by generate_inp_cases.py")
    w("**")
    w("** ============================================================")
    w("** NODES")
    w("** ============================================================")
    w("*Node")
    for i, z in enumerate(z_nodes, start=1):
        w(f" {i},  0., 0., {z:.6f}")

    w("**")
    w("** ============================================================")
    w("** ELEMENTS - TRUNK (B31 beam, Z-axis)")
    w("** ============================================================")
    w("*Element, type=B31, elset=TRUNK")
    for i in range(1, n_elem + 1):
        w(f" {i}, {i}, {i+1}")

    # Crown mass element (dummy element at tip node)
    w("**")
    w("*Element, type=MASS, elset=CROWN_MASS_ELEM")
    w(f" 101, {n_node}")

    w("**")
    w("** ============================================================")
    w("** ELEMENT SETS (1 per element for tapered sections)")
    w("** ============================================================")
    for i in range(1, n_elem + 1):
        w(f"*Elset, elset=E{i}")
        w(f" {i},")

    w("**")
    w("** ============================================================")
    w("** NODE SETS")
    w("** ============================================================")
    w("*Nset, nset=BASE")
    w(" 1,")
    w("*Nset, nset=TIP")
    w(f" {n_node},")

    w("**")
    w("** ============================================================")
    w("** BEAM SECTIONS (tapered, CIRC, n1=(1,0,0))")
    w("** ============================================================")
    for i in range(1, n_elem + 1):
        R_i = section_radius_at(D_base, D_tip, n_elem, i - 1)
        w(f"*Beam Section, elset=E{i}, material=WOOD, section=CIRC")
        w(f" {R_i:.6f},")
        w(" 1.0, 0.0, 0.0")

    w("**")
    w("** ============================================================")
    w("** MATERIAL")
    w("** ============================================================")
    w("*Material, name=WOOD")
    w("*Elastic")
    w(f" {E_WOOD:.4e}, {NU_WOOD}")
    w("*Density")
    w(f" {RHO_WOOD},")
    w("** Rayleigh damping: zeta=2% at f1 range (Sellier & Fourcaud 2009)")
    w("** alpha=0.0754, beta=0.00472  -> zeta~2% at 0.3~1.0 Hz")
    w("*Damping, alpha=0.0754, beta=0.00472")

    w("**")
    w("** ============================================================")
    w("** CROWN MASS")
    w("** ============================================================")
    w("*Mass, elset=CROWN_MASS_ELEM")
    w(f" {M_crown:.2f},")

    w("**")
    w("** ============================================================")
    w("** BOUNDARY CONDITIONS (base: fully fixed)")
    w("** ============================================================")
    w("*Boundary")
    w(" 1, 1, 6")

    w("**")
    w("** ============================================================")
    w("** AMPLITUDE - KAIMAL X (mean + turbulence, always >= 0)")
    w(f"** F_X_ref = {F_X_ref:.2f} N  (Cload 입력값)")
    w(f"** A_crown = {A_crown:.2f} m2, U_bar = {U_BAR} m/s")
    w("** ============================================================")
    w("*Amplitude, name=KAIMAL_X, definition=TABULAR")
    for ln in amp_X_lines:
        w(ln)

    w("**")
    w("** ============================================================")
    w("** AMPLITUDE - KAIMAL Y (lateral turbulence only, bidirectional)")
    w(f"** F_Y_ref = {F_Y_ref:.2f} N  (Cload 입력값)")
    w(f"** sigma_v = {I_V_RAT * I_U * U_BAR:.2f} m/s  (= {I_V_RAT} x sigma_u)")
    w("** ============================================================")
    w("*Amplitude, name=KAIMAL_Y, definition=TABULAR")
    for ln in amp_Y_lines:
        w(ln)

    w("**")
    w("** ============================================================")
    w("** STEP 1: MODAL ANALYSIS (고유진동수)")
    w("** ============================================================")
    w("*Step, name=FREQUENCY, perturbation")
    w("*Frequency, eigensolver=Lanczos, normalization=displacement")
    w(" 5,")
    w("*Output, field")
    w("*Node Output")
    w(" U,")
    w("*End Step")

    w("**")
    w("** ============================================================")
    w("** STEP 2: WIND DYNAMIC ANALYSIS")
    w("** ============================================================")
    w("*Step, name=WIND_DYNAMIC, nlgeom=NO, inc=15000")
    w("*Dynamic, direct")
    w(f" {DT}, {T_SIM:.1f}")
    w("**")
    w("** X-direction: mean wind + turbulence")
    w("*Cload, amplitude=KAIMAL_X")
    w(f" {n_node}, 1, {F_X_ref:.4f}")
    w("**")
    w("** Y-direction: lateral turbulence only")
    w("*Cload, amplitude=KAIMAL_Y")
    w(f" {n_node}, 2, {F_Y_ref:.4f}")
    w("**")
    w("*Output, field, frequency=100")
    w("*Node Output")
    w(" U, V, A")
    w("*Element Output, elset=TRUNK")
    w(" SF, SM")
    w("**")
    w("*Output, history, frequency=1")
    w("** Base element: bending moment SM1, SM2")
    w("*Element Output, elset=E1")
    w(" SF1, SF2, SF3, SM1, SM2")
    w("** Tip displacement")
    w("*Node Output, nset=TIP")
    w(" U1, U2")
    w("** Base reaction force")
    w("*Node Output, nset=BASE")
    w(" RF1, RF2")
    w("*End Step")

    # ------------------------------------------------------------------
    with open(fname, "w", encoding="ascii", errors="replace") as f:
        f.write("\n".join(lines) + "\n")

    n_lines = len(lines)
    print(f"  -> {fname.name}  ({n_lines} lines,  F_X_ref={F_X_ref:.1f}N,  M_crown={M_crown:.1f}kg)")
    return fname


# ======================================================================
# 5. 전체 케이스 생성
# ======================================================================
print("=" * 62)
print("  ABAQUS INP 파일 자동 생성")
print(f"  케이스: {len(CASES)}개  (H x DBH)")
print(f"  다축하중: X(Kaimal 평균+난류) + Y(횡방향 난류)")
print(f"  시뮬레이션: {T_SIM:.0f}s, dt={DT}s, {N_STEPS}스텝")
print("=" * 62)

generated = []
for case in CASES:
    f = write_inp(case)
    generated.append(f)

print()
print("=" * 62)
print("  생성 완료")
print("=" * 62)
for f in generated:
    print(f"  {f}")

print()
print("  ABAQUS 실행 방법:")
print("  1. inp_cases/ 폴더를 ABAQUS 실행 컴퓨터에 복사")
print("  2. 각 케이스 실행:")
for case in CASES:
    H   = case["H"]
    DBH = case["DBH"]
    name = f"tree_H{H:.0f}m_DBH{DBH*100:.0f}cm"
    print(f"     abaqus job={name} interactive")
print()
print("  3. 해석 완료 후 각 ODB에 대해:")
print("     abaqus python extract_odb.py -- {odb_name}")
print("     (또는 extract_odb.py 내부 ODB_PATH 수정 후 실행)")
