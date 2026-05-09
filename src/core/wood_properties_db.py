"""
wood_properties_db.py
---------------------
가로수 주요 수종의 목재 역학적 물성 및 S-N 곡선 파라미터 데이터베이스

데이터 출처:
  [1] Panshin & de Zeeuw (1980). Textbook of Wood Technology. McGraw-Hill.
  [2] Kretschmann, D.E. (2010). Wood Handbook Ch.5 — Mechanical Properties.
        USDA Forest Service, FPL-GTR-190.
  [3] Kohara, M. & Okuyama, T. (1992). Bending fatigue of small clear specimens.
        Mokuzai Gakkaishi, 38(6), 529–536.
  [4] Tsai, K.T. & Ansell, M.P. (1990). Fatigue properties of wood in flexure.
        J. Mater. Sci., 25, 865–878.
  [5] Bohannan, B. (1966). Effect of size on bending strength of wood members.
        USDA FS Research Paper FPL-56.
  [6] Gerhards, C.C. (1979). Time-related effects on wood strength.
        Wood Sci., 11(3), 139–151.
  [7] Korea Forest Research Institute (2012). Wood Properties of Korean Trees.
        국립산림과학원 연구보고서 (수종별 목재 물성 측정치 수록).
  [8] Sellier, D. & Fourcaud, T. (2009). Crown structure and wood properties.
        Plant Sci., 177, 33–42.

주의:
  - 아래 값은 소건재(MC≈12%) 기준 평균값
  - 생재(MC>30%)는 E, MOR 각각 약 30~50% 감소
  - 실제 연구에서는 대상 수종 직접 실측값 사용 권장
  - S-N 파라미터(a, b)는 소재 굽힘 피로시험 유추값 — 문헌 불충분 주의
"""

import numpy as np

# ══════════════════════════════════════════════════════════════════════
# 정적 역학 물성 (Mechanical Properties — Static)
# ══════════════════════════════════════════════════════════════════════
# 키: 수종명 (학명)
# E_L     : 종방향 탄성계수 [GPa]  (섬유방향)
# G_LR    : 전단탄성계수  [GPa]  (종-반경방향)
# nu_LR   : 포아송비 (종-반경)
# MOR     : 파단계수 (정적 굽힘 강도) [MPa]  — 3점 굽힘, 소건재
# MCS     : 섬유방향 압축강도 [MPa]
# density : 기건밀도 [kg/m³]  (MC≈12%)
# density_green: 생재밀도 [kg/m³]  (MC>30%, 수목 실측에 더 적합)
# ref     : 주요 참고문헌 번호

STATIC_PROPS = {
    "은행나무 (Ginkgo biloba)": {
        "E_L_GPa":       7.0,
        "G_LR_GPa":      0.55,
        "nu_LR":         0.33,
        "MOR_MPa":       68.0,
        "MCS_MPa":       38.0,
        "density_dry":   530,
        "density_green": 780,
        "note": "침엽수-활엽수 중간 특성; 가로수 중 가장 많이 식재",
        "ref": "[1][2][7]"
    },
    "왕벚나무 (Prunus yedoensis)": {
        "E_L_GPa":       9.8,
        "G_LR_GPa":      0.72,
        "nu_LR":         0.35,
        "MOR_MPa":       88.0,
        "MCS_MPa":       49.0,
        "density_dry":   640,
        "density_green": 900,
        "note": "벚나무류 중 강도 높음; 봄철 가로수 대표 수종",
        "ref": "[1][2][7]"
    },
    "플라타너스 (Platanus x acerifolia)": {
        "E_L_GPa":       8.1,
        "G_LR_GPa":      0.62,
        "nu_LR":         0.34,
        "MOR_MPa":       72.0,
        "MCS_MPa":       40.0,
        "density_dry":   600,
        "density_green": 850,
        "note": "유럽 플라타너스 교잡종; 국내 측정 데이터 부족",
        "ref": "[1][2]"
    },
    "느티나무 (Zelkova serrata)": {
        "E_L_GPa":       12.5,
        "G_LR_GPa":      0.89,
        "nu_LR":         0.30,
        "MOR_MPa":       105.0,
        "MCS_MPa":       58.0,
        "density_dry":   730,
        "density_green": 980,
        "note": "국내 활엽수 중 강도 최상위; 고목 보존 가치 높음",
        "ref": "[1][7]"
    },
    "이팝나무 (Chionanthus retusus)": {
        "E_L_GPa":       10.2,
        "G_LR_GPa":      0.75,
        "nu_LR":         0.32,
        "MOR_MPa":       91.0,
        "MCS_MPa":       52.0,
        "density_dry":   680,
        "density_green": 930,
        "note": "최근 가로수 식재 증가; 역학 데이터 부족 — 유추값",
        "ref": "[7] 유추"
    },
    "메타세쿼이아 (Metasequoia glyptostroboides)": {
        "E_L_GPa":       6.5,
        "G_LR_GPa":      0.48,
        "nu_LR":         0.35,
        "MOR_MPa":       58.0,
        "MCS_MPa":       32.0,
        "density_dry":   420,
        "density_green": 680,
        "note": "침엽수; 담양 등 관광 가로수 대표. 상대적으로 약함",
        "ref": "[1][7]"
    },
}

# ══════════════════════════════════════════════════════════════════════
# S-N 곡선 파라미터 (Fatigue — Flexural)
# ══════════════════════════════════════════════════════════════════════
# 형식: σ_a / MOR = a - b × log10(N_f)
# 역산: N_f = 10^((a - σ_a/MOR) / b)
#
# 응력비 R = σ_min / σ_max
#   R = -1 : 완전 반전 (양방향 굽힘)
#   R =  0 : 편진 (단방향, 풍하중 근사)
#   R = 0.1: 약한 평균응력 포함
#
# 참고: 목재 피로 연구는 절대적으로 부족함
#  - 침엽수(소나무, 가문비): 데이터 다수
#  - 활엽수(참나무류): 제한적
#  - 도시 가로수 수종: 거의 없음 → 유사종 유추
#
# a_R0, b_R0: R=0 (편진) 조건
# a_Rm1, b_Rm1: R=-1 (완전 반전) 조건
# N_limit: 피로한계 반복수 (이 이상은 파단 없음으로 가정, 보통 10^7)
# sigma_limit_ratio: 피로한계 응력비 σ_e/MOR

SN_PARAMS = {
    "은행나무 (Ginkgo biloba)": {
        # 침엽수 유사 특성 → Kohara & Okuyama (1992) 적용
        "a_R0":           0.97,
        "b_R0":           0.072,
        "a_Rm1":          0.95,
        "b_Rm1":          0.090,
        "N_limit":        1e7,
        "sigma_limit_ratio": 0.40,
        "reliability":   "낮음 — 침엽수 유추, 실측 필요",
        "ref": "[3] 유추"
    },
    "왕벚나무 (Prunus yedoensis)": {
        # 활엽수 Tsai & Ansell (1990) 참고 (유럽 벚나무 유사종)
        "a_R0":           0.98,
        "b_R0":           0.065,
        "a_Rm1":          0.96,
        "b_Rm1":          0.082,
        "N_limit":        1e7,
        "sigma_limit_ratio": 0.42,
        "reliability":   "중간 — 유사종 데이터 존재",
        "ref": "[4]"
    },
    "플라타너스 (Platanus x acerifolia)": {
        "a_R0":           0.97,
        "b_R0":           0.070,
        "a_Rm1":          0.95,
        "b_Rm1":          0.088,
        "N_limit":        1e7,
        "sigma_limit_ratio": 0.40,
        "reliability":   "낮음 — 국내외 직접 데이터 부재",
        "ref": "[2] 유추"
    },
    "느티나무 (Zelkova serrata)": {
        # 경질 활엽수 → 피로 저항 상대적으로 높음
        "a_R0":           0.98,
        "b_R0":           0.060,
        "a_Rm1":          0.96,
        "b_Rm1":          0.075,
        "N_limit":        1e7,
        "sigma_limit_ratio": 0.44,
        "reliability":   "낮음 — 국내 연구 전무, 참나무류 유추",
        "ref": "[3][4] 유추"
    },
    "메타세쿼이아 (Metasequoia glyptostroboides)": {
        # 침엽수 → 소나무류 데이터 참고
        "a_R0":           0.96,
        "b_R0":           0.078,
        "a_Rm1":          0.94,
        "b_Rm1":          0.095,
        "N_limit":        1e7,
        "sigma_limit_ratio": 0.38,
        "reliability":   "낮음 — 낙우송과 유사종 유추",
        "ref": "[3] 유추"
    },
}

# ══════════════════════════════════════════════════════════════════════
# 수분 함량(MC) 보정 계수
# ══════════════════════════════════════════════════════════════════════
# 생재(MC>FSP=30%)는 소건재 대비 강도 감소
# 가로수는 MC 계절 변동 큼 (하절기 생재 근사, 동절기 기건 근사)
#
# 보정 공식 (Kretschmann 2010):
#   P_green = P_dry × CM_factor
# 참고: 정밀 해석 시 계절별 MC 변동을 파라미터로 추가 가능

MC_CORRECTION = {
    "E_L":  0.68,   # 생재 탄성계수 ≈ 기건의 68%
    "MOR":  0.60,   # 생재 MOR     ≈ 기건의 60%
    "MCS":  0.55,   # 생재 압축강도 ≈ 기건의 55%
}

# ══════════════════════════════════════════════════════════════════════
# 수목 동적 특성 참고값 (문헌 실측)
# ══════════════════════════════════════════════════════════════════════
# 1차 고유진동수 f1 [Hz]와 감쇠비 ζ [-]
# 수목별 실측값 범위 — ABAQUS 모달해석 검증 기준

DYNAMIC_REF = {
    "일반 가로수 (H=5~10m)": {
        "f1_Hz":    (0.20, 0.80),
        "zeta":     (0.01, 0.05),
        "ref": "Sellier & Fourcaud (2009); Moore et al. (2005)"
    },
    "성목 (H>15m)": {
        "f1_Hz":    (0.10, 0.30),
        "zeta":     (0.01, 0.03),
        "ref": "James et al. (2006); Spatz & Theckes (2013)"
    },
    "메타세쿼이아 (H=15~20m)": {
        "f1_Hz":    (0.15, 0.40),
        "zeta":     (0.015, 0.04),
        "ref": "국내 실측 데이터 부족 — 추정값"
    },
}

# ══════════════════════════════════════════════════════════════════════
# 편의 함수
# ══════════════════════════════════════════════════════════════════════
def get_props(species, green=False):
    """
    수종명으로 정적 물성 + S-N 파라미터 반환
    green=True: 생재 보정 적용
    """
    sp = STATIC_PROPS.get(species)
    sn = SN_PARAMS.get(species)
    if sp is None:
        raise KeyError(f"수종 '{species}' 없음. 사용 가능: {list(STATIC_PROPS.keys())}")

    result = dict(sp)
    if sn:
        result.update(sn)

    if green:
        result["E_L_GPa"] *= MC_CORRECTION["E_L"]
        result["MOR_MPa"] *= MC_CORRECTION["MOR"]
        result["MCS_MPa"] *= MC_CORRECTION["MCS"]
        result["_state"]   = "생재 (MC>30%) 보정 적용"
    else:
        result["_state"]   = "소건재 (MC≈12%)"

    return result


def sn_Nf(sigma_a, species, R=0, green=False):
    """
    파단 반복수 N_f 계산
    sigma_a: 응력 진폭 [MPa]
    R:  응력비 (0 또는 -1)
    """
    props = get_props(species, green=green)
    MOR   = props["MOR_MPa"]

    if R == 0:
        a, b = props["a_R0"], props["b_R0"]
    elif R == -1:
        a, b = props["a_Rm1"], props["b_Rm1"]
    else:
        raise ValueError("R은 0 또는 -1만 지원")

    ratio = sigma_a / MOR
    if ratio <= 0:
        return np.inf
    if ratio >= a:
        return 1.0
    # 피로한계 확인
    sigma_limit = props["sigma_limit_ratio"] * MOR
    if sigma_a <= sigma_limit:
        return np.inf
    return 10.0 ** ((a - ratio) / b)


def print_summary():
    """전체 DB 요약 출력"""
    print("=" * 70)
    print(f"  {'수종':<28} {'E_L':>6} {'MOR':>6} {'ρ_dry':>6} {'S-N신뢰도'}")
    print(f"  {'':28} {'GPa':>6} {'MPa':>6} {'kg/m³':>6}")
    print("=" * 70)
    for sp in STATIC_PROPS:
        p  = STATIC_PROPS[sp]
        sn = SN_PARAMS.get(sp, {})
        rel = sn.get("reliability", "-")
        print(f"  {sp[:28]:<28} {p['E_L_GPa']:>6.1f} {p['MOR_MPa']:>6.1f} "
              f"{p['density_dry']:>6} {rel}")
    print("=" * 70)
    print("\n  * 생재 보정: E×0.68, MOR×0.60 적용")
    print("  * S-N 형식: σ_a/MOR = a - b·log10(N_f),  R=0 기준")


# ── 직접 실행 시 요약 출력 ──────────────────────────────────────────
if __name__ == "__main__":
    print_summary()

    print("\n  [사용 예시]")
    sp = "느티나무 (Zelkova serrata)"
    props = get_props(sp, green=False)
    print(f"\n  수종: {sp}")
    print(f"  탄성계수 E_L = {props['E_L_GPa']} GPa")
    print(f"  MOR         = {props['MOR_MPa']} MPa")
    print(f"  밀도(기건)   = {props['density_dry']} kg/m³")

    print(f"\n  S-N 파단 반복수 (R=0, 기건):")
    for sigma_a in [10, 20, 30, 40, 50]:
        Nf = sn_Nf(sigma_a, sp, R=0, green=False)
        print(f"    σ_a={sigma_a:2d} MPa → N_f = {Nf:.2e}")

    print(f"\n  생재 보정 시 MOR = {props['MOR_MPa'] * MC_CORRECTION['MOR']:.1f} MPa")
