"""
extract_odb.py
--------------
ABAQUS ODB 파일에서 결과 추출 -> 피로해석 연계
  - 고유진동수 f1 (FREQUENCY step)
  - 선단 변위 최대값 (U1, U2 at tip node)
  - 근원부 굽힘응력 시계열 (SM1 + SM2 -> M_res)

실행:
    abaqus python extract_odb.py
    (ODB_PATH 아래에서 수정 후 실행, 또는 명령행 인수 사용)

    abaqus python extract_odb.py -- tree_H8m_DBH20cm.odb

ABAQUS 보 요소 단면력:
    SF1 = 축력 [N]
    SF2 = 전단력 2방향 [N]
    SF3 = 전단력 3방향 [N]
    SM1 = 굽힘모멘트 1축 [N.m]
    SM2 = 굽힘모멘트 2축 [N.m]  <- X방향 풍하중 지배 성분
    SM3 = 비틀림 [N.m]

다축하중 추가 후:
    M_res = sqrt(SM1^2 + SM2^2)  -> 합성 굽힘모멘트
    sigma_bending = M_res * c / I  [Pa -> MPa]
"""

from odbAccess import openOdb
import os
import csv
import math
import sys

# ======================================================================
# 1. 설정
# ======================================================================
# 명령행 인수로 ODB 파일명 지정 가능
# abaqus python extract_odb.py -- mymodel.odb
if len(sys.argv) > 1:
    ODB_PATH = sys.argv[-1]
    if not ODB_PATH.endswith('.odb'):
        ODB_PATH += '.odb'
else:
    ODB_PATH = 'abaqus_tree_model.odb'

STEP_FREQ = 'FREQUENCY'
STEP_DYN  = 'WIND_DYNAMIC'
OUT_DIR   = 'output'

# 근원부 단면 (Element 1)
# generate_inp_cases.py 와 동일한 계산 사용
# ODB 파일명에서 DBH 자동 파싱 시도
D_BASE = 0.24   # 기본값 (DBH=0.20m -> D_base=0.24m)
try:
    import re
    m_dbh = re.search(r'DBH(\d+)cm', ODB_PATH)
    if m_dbh:
        DBH_cm = float(m_dbh.group(1))
        D_BASE = (DBH_cm / 100.0) * 1.2
        print('ODB 파일명에서 DBH={}cm 파싱 -> D_base={:.4f}m'.format(DBH_cm, D_BASE))
except Exception:
    pass

I_BASE = math.pi * D_BASE**4 / 64.0
C_BASE = D_BASE / 2.0

print('=' * 60)
print('  ABAQUS ODB 결과 추출 (다축하중 대응)')
print('=' * 60)
print('ODB: {}'.format(ODB_PATH))
print('단면: D={} m, I={:.4e} m4, c={} m'.format(D_BASE, I_BASE, C_BASE))

if not os.path.exists(OUT_DIR):
    os.makedirs(OUT_DIR)

# ======================================================================
# 2. ODB 열기
# ======================================================================
odb  = openOdb(path=ODB_PATH, readOnly=True)

# ======================================================================
# 3. 고유진동수 추출 (FREQUENCY step)
# ======================================================================
f1_Hz   = None
mode_freqs = []

if STEP_FREQ in odb.steps:
    freq_step = odb.steps[STEP_FREQ]
    print('\n[ FREQUENCY Step - 고유진동수 ]')
    for frame in freq_step.frames:
        if frame.frameId == 0:
            continue   # 초기 프레임 건너뜀
        try:
            freq_val = frame.frequency    # [Hz]
            mode_freqs.append(freq_val)
            print('  Mode {}: f = {:.4f} Hz  (omega = {:.3f} rad/s)'.format(
                frame.frameId, freq_val, 2*math.pi*freq_val))
        except AttributeError:
            pass
    if mode_freqs:
        f1_Hz = mode_freqs[0]
        print('  -> 1차 고유진동수 f1 = {:.4f} Hz'.format(f1_Hz))
else:
    print('\n[경고] FREQUENCY step 없음 (기존 모델)')

# ======================================================================
# 4. 동적 해석 결과 추출 (WIND_DYNAMIC step)
# ======================================================================
if STEP_DYN not in odb.steps:
    print('[오류] WIND_DYNAMIC step 없음.')
    odb.close()
    raise SystemExit(1)

step = odb.steps[STEP_DYN]
print('\nStep: {}  |  프레임 수: {}'.format(STEP_DYN, len(step.frames)))

# ── History Regions 확인 ────────────────────────────────────────────
print('\n[ History Regions ]')
for rname in step.historyRegions.keys():
    region = step.historyRegions[rname]
    print('  {} : {}'.format(rname, list(region.historyOutputs.keys())))

# ── Element 1 Region 탐색 (SM1, SM2 우선) ───────────────────────────
elem_region = None
elem_rname  = None

for rname, region in step.historyRegions.items():
    if 'Element' in rname:
        outputs = list(region.historyOutputs.keys())
        if 'SM2' in outputs or 'SM1' in outputs:
            elem_region = region
            elem_rname  = rname
            break

if elem_region is None:
    for rname, region in step.historyRegions.items():
        if 'Element' in rname:
            elem_region = region
            elem_rname  = rname
            break

if elem_region is None:
    print('[오류] Element history region 없음.')
    odb.close()
    raise SystemExit(1)

print('\n선택 Region: {}'.format(elem_rname))
outputs = list(elem_region.historyOutputs.keys())
print('사용 가능 변수: {}'.format(outputs))

# ── SM1, SM2 추출 ───────────────────────────────────────────────────
sm2_data, sm1_data = [], []

if 'SM2' in outputs:
    sm2_data = list(elem_region.historyOutputs['SM2'].data)
    print('SM2 추출 완료: {} 스텝'.format(len(sm2_data)))

if 'SM1' in outputs:
    sm1_data = list(elem_region.historyOutputs['SM1'].data)
    print('SM1 추출 완료: {} 스텝'.format(len(sm1_data)))

if not sm2_data and not sm1_data:
    print('[경고] SM1/SM2 없음 -> SF3으로 대체 (부정확)')
    sm2_data = list(elem_region.historyOutputs.get(
        'SF3', list(elem_region.historyOutputs.values())[0]).data)

# 시간 배열은 SM2 기준 (없으면 SM1)
ref_data = sm2_data if sm2_data else sm1_data
t_arr    = [d[0] for d in ref_data]
n        = len(t_arr)

sm2_arr = [d[1] for d in sm2_data] if sm2_data else [0.0] * n
sm1_arr = [d[1] for d in sm1_data] if sm1_data else [0.0] * n

# 합성 굽힘모멘트 M_res = sqrt(SM1^2 + SM2^2)
m_res_arr    = [math.sqrt(sm1_arr[i]**2 + sm2_arr[i]**2) for i in range(n)]
sigma_arr    = [m * C_BASE / I_BASE * 1e-6 for m in m_res_arr]  # [MPa]
sigma_sm2    = [abs(sm2_arr[i]) * C_BASE / I_BASE * 1e-6 for i in range(n)]  # SM2 단독

sigma_max_res  = max(sigma_arr)
sigma_max_sm2  = max(sigma_sm2)
sigma_mean     = sum(sigma_arr) / n
sigma_std      = (sum((s - sigma_mean)**2 for s in sigma_arr) / n) ** 0.5

print('\n' + '=' * 60)
print('  근원부 굽힘 응력 통계')
print('=' * 60)
if sm1_data:
    sm1_max = max(abs(v) for v in sm1_arr)
    sm2_max = max(abs(v) for v in sm2_arr)
    print('  SM1 최대: {:.1f} N.m'.format(sm1_max))
    print('  SM2 최대: {:.1f} N.m'.format(sm2_max))
    print('  SM2/SM1 비: {:.2f}  (>1이면 SM2 지배)'.format(
        sm2_max / sm1_max if sm1_max > 0 else float('inf')))
print('  sigma (SM2 단독) 최대 = {:.3f} MPa'.format(sigma_max_sm2))
print('  sigma (M_res) 최대   = {:.3f} MPa'.format(sigma_max_res))
print('  sigma 평균 = {:.3f} MPa'.format(sigma_mean))
print('  sigma std  = {:.3f} MPa'.format(sigma_std))
if f1_Hz:
    print('  고유진동수 f1 = {:.4f} Hz'.format(f1_Hz))
print('=' * 60)

# DAF 계산 (SM2 단독 기준, Python 단독 ~18 MPa 대비)
SIGMA_PY_REF = 18.0
DAF_sm2 = sigma_max_sm2 / SIGMA_PY_REF if SIGMA_PY_REF > 0 else 1.0
DAF_res = sigma_max_res / SIGMA_PY_REF if SIGMA_PY_REF > 0 else 1.0
print('  DAF (SM2 단독) = {:.3f}  (Python 단독 {}MPa 대비)'.format(DAF_sm2, SIGMA_PY_REF))
print('  DAF (M_res)    = {:.3f}'.format(DAF_res))

# ======================================================================
# 5. 선단 변위 추출 (Tip node)
# ======================================================================
u1_tip_max = 0.0
u2_tip_max = 0.0

node_region = None
for rname, region in step.historyRegions.items():
    if 'Node' in rname:
        node_region = region
        break

if node_region:
    node_outputs = list(node_region.historyOutputs.keys())
    if 'U1' in node_outputs:
        u1_arr = [d[1] for d in node_region.historyOutputs['U1'].data]
        u1_tip_max = max(abs(v) for v in u1_arr)
        print('\n  선단부 X변위 (U1) 최대: {:.4f} m'.format(u1_tip_max))
    if 'U2' in node_outputs:
        u2_arr = [d[1] for d in node_region.historyOutputs['U2'].data]
        u2_tip_max = max(abs(v) for v in u2_arr)
        print('  선단부 Y변위 (U2) 최대: {:.4f} m'.format(u2_tip_max))
    if u1_tip_max > 0 or u2_tip_max > 0:
        u_res_max = math.sqrt(u1_tip_max**2 + u2_tip_max**2)
        print('  선단부 합성변위 최대:    {:.4f} m'.format(u_res_max))

# ======================================================================
# 6. CSV 저장
# ======================================================================
# 파일명: ODB 이름에서 .odb 제거 후 _stress.csv
base_name = os.path.splitext(os.path.basename(ODB_PATH))[0]
out_csv   = os.path.join(OUT_DIR, base_name + '_stress.csv')
# 기존 코드와의 호환성: abaqus_stress_history.csv 도 함께 저장
compat_csv = os.path.join(OUT_DIR, 'abaqus_stress_history.csv')

def write_csv(path, use_compat=False):
    with open(path, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['# ABAQUS Dynamic Analysis - Stress History'])
        writer.writerow(['# ODB: {}'.format(ODB_PATH)])
        if f1_Hz:
            writer.writerow(['# f1={:.4f} Hz  U1_tip_max={:.4f} m  U2_tip_max={:.4f} m'.format(
                f1_Hz, u1_tip_max, u2_tip_max)])
        writer.writerow(['# sigma_max_SM2={:.3f} MPa  sigma_max_Mres={:.3f} MPa'.format(
            sigma_max_sm2, sigma_max_res)])
        writer.writerow(['# DAF_SM2={:.3f}  DAF_Mres={:.3f}'.format(DAF_sm2, DAF_res)])
        writer.writerow(['time_s', 'SM1_Nm', 'SM2_Nm', 'Mres_Nm', 'sigma_bending_MPa'])
        for i in range(n):
            writer.writerow([
                '{:.4f}'.format(t_arr[i]),
                '{:.4f}'.format(sm1_arr[i]),
                '{:.4f}'.format(sm2_arr[i]),
                '{:.4f}'.format(m_res_arr[i]),
                '{:.6f}'.format(sigma_arr[i])
            ])

write_csv(out_csv)
write_csv(compat_csv)
print('\n  CSV 저장 완료: {}'.format(out_csv))
print('  CSV 저장 완료: {}  (호환용)'.format(compat_csv))
print('  -> python fatigue_from_abaqus.py 로 피로해석 연계 가능')

# ======================================================================
# 7. 요약 파일 저장 (케이스별 비교용)
# ======================================================================
summary_path = os.path.join(OUT_DIR, 'abaqus_cases_summary.csv')
summary_exists = os.path.exists(summary_path)

with open(summary_path, 'a') as f:
    writer = csv.writer(f)
    if not summary_exists:
        writer.writerow(['ODB', 'f1_Hz', 'sigma_max_SM2_MPa', 'sigma_max_Mres_MPa',
                         'DAF_SM2', 'DAF_Mres', 'U1_tip_max_m', 'U2_tip_max_m'])
    writer.writerow([
        base_name,
        '{:.4f}'.format(f1_Hz) if f1_Hz else 'N/A',
        '{:.3f}'.format(sigma_max_sm2),
        '{:.3f}'.format(sigma_max_res),
        '{:.3f}'.format(DAF_sm2),
        '{:.3f}'.format(DAF_res),
        '{:.4f}'.format(u1_tip_max),
        '{:.4f}'.format(u2_tip_max),
    ])

print('  케이스 요약 추가: {}'.format(summary_path))

# ======================================================================
# 8. ODB 닫기
# ======================================================================
odb.close()
print('  ODB 닫기 완료.')
print('\n  다음 케이스 실행 방법:')
print('  abaqus python extract_odb.py -- tree_H6m_DBH15cm.odb')
