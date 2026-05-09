# Research Plan: Probabilistic Eccentric Decay Monte Carlo

## 1. 연구 배경 및 목적

1차 논문은 동심 hollow-cylinder 부후 가정을 사용하여 응력 증폭 계수를 다음과 같이 정의했다:

$$\text{SAF}_{\text{conc}}(r_d/R) = \frac{1}{1 - (r_d/R)^4}$$

그러나 실제 부후는 편심·비대칭이며 (Mattheck & Breloer 1994), 동일 $r_d/R$ 에서도 단면 형상에 따라 응력 분포가 20–40% 변동한다. 본 프로젝트는 이 불확실성을 정량화하여 결정론적 피로수명을 *분포*로 대체한다.

## 2. 단계별 계획

### Phase A. 편심 기하 모델 (1–2개월)

**A1. 편심 단면 파라미터화**
- 결정 변수: 편심률 $e/R$, 부후 종횡비 $\alpha = a/b$, 잔존벽 두께 분포 $t(\theta)$
- Mattheck & Breloer 도판을 정량 데이터로 추출
- 형태 라이브러리: 동심, 편심 원형, 타원형, 초승달형, 다중 공동

**A2. 단면 모듈러스 계산**
- 적분 기반: $S_{\text{net}}(\theta_{\text{wind}}) = I_{\text{net}}/c_{\text{max}}$
- 풍향 $\theta_{\text{wind}}$ 에 대한 sweep → 최악 방향 식별
- ABAQUS 부분 검증 (parametric inp 자동 생성, 5–10 케이스)

**산출**: `src/eccentric_decay/geometry.py`, `stress_amplification.py`

### Phase B. Monte Carlo 샘플링 (1개월)

**B1. 분포 가정 (literature-based priors)**
| 변수 | 분포 | 출처 |
|------|------|------|
| $e/R$ | Beta(2, 5) — 편심 쪽 skewed | Mattheck 도판 + sonic tomography 사례 |
| $\alpha$ | LogNormal(μ, σ) | 동상 |
| 풍향 | Uniform[0, 2π] | 무지향성 가정 |

**B2. MC 실행**
- $N = 10,000$ 샘플 / $r_d/R$ 수준 (5단계: 0.3, 0.5, 0.6, 0.7, 0.8)
- 각 샘플: SAF 계산 → rainflow 적용 → Miner 누적 → 피로수명
- 풍속은 1차 논문의 Weibull 그대로 (서울/인천/제주)

**산출**: `data/generated/mc_results.parquet`, fatigue life CDF

### Phase C. 신뢰구간 + 정책 기여 (1개월)

**C1. 점검 주기 환산**
- 95% CI 하한 → mandatory removal threshold
- 50% CI → 권고 점검 주기
- 1차 논문의 단일 임계값 ($r_d/R = 0.7$) 과 비교

**C2. 민감도 분석**
- 가장 영향력 큰 기하 변수 식별 (Sobol indices)
- 우선 측정해야 할 현장 변수 도출 → sonic tomography 측정 프로토콜 권고

### Phase D. 논문 작성 (2개월)

**구성**:
1. Introduction — 1차 논문 직접 인용, §4.5 한계 명시
2. Methods — Phase A–C
3. Results — fatigue life CDF, regional comparison, 정책 매트릭스
4. Discussion — 동심 가정의 보수성/비보수성 정량화, 측정 프로토콜
5. Conclusion

**목표 저널**: Trees–Structure and Function (1차 논문과 동일 → 시리즈 평가 가능)

## 3. 마일스톤

| 시점 | 목표 |
|------|------|
| Month 1 | Phase A1–A2 완료, 편심 기하 코드 검증 |
| Month 2 | Phase B 완료, MC 결과 1차 산출 |
| Month 3 | Phase C 완료, 정책 매트릭스 |
| Month 4–5 | 논문 초안 |
| Month 6 | 외부 리뷰 + 투고 |

## 4. 핵심 의존성

- 1차 논문 코드 (`src/core/`) 의 wind/fatigue/Weibull 모듈
- ABAQUS 2026 (검증 단계 한정)
- 1차 논문의 stress history (`data/inherited/stress_history/`) — 회귀 테스트용

## 5. 위험 요소

| 위험 | 완화 |
|------|------|
| 편심 기하 분포가 임의적 | Mattheck 도판 + 가능하면 sonic tomography 데이터 협력 (국립산림과학원) |
| ABAQUS 자동화가 병목 | beam 단면 모델은 Python 직접 적분으로 대체 가능 |
| MC 수렴 불충분 | 10⁴ → 10⁵ 확장 가능, importance sampling 검토 |
| 결과가 1차 논문과 모순 | e = 0 한계에서 정확히 일치하는지 회귀 테스트 (필수) |

## 6. 예상 산출물

1. **논문 1편** — Trees 또는 Engineering Failure Analysis
2. **공개 코드** — Github (1차 논문 후속 시리즈로 묶기)
3. **MC 결과 데이터셋** — Zenodo 또는 부록
4. **위험 점검 매트릭스** (1쪽 표) — 한국 가로수 관리 기관 배포 가능
