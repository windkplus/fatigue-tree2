# CLAUDE.md — fatigue-tree2 (Eccentric Decay Monte Carlo)

> 본 프로젝트는 `fatigue-tree` 1차 논문(Springer/Elsevier 동시 투고)의 **§4.5 다섯 번째 한계** — concentric hollow-cylinder 부후 가정 — 를 직접 확장하는 후속 연구입니다.
> 1차 논문 정보·데이터는 `docs/` 와 `data/inherited/` 에 동결 보관되어 있습니다.

---

## 프로젝트 개요

**제목 (작업명)**: Probabilistic Eccentric Decay Geometry: Monte Carlo Fatigue Life Confidence Intervals for Urban Street Trees

**연구 질문**:
1. 동심 가정 대비 편심·비대칭 부후가 피로수명에 미치는 영향을 정량화할 수 있는가?
2. 현장에서 실측 가능한 부후 기하 변수(편심률, 종횡비, 잔존벽 두께 분포)에 대해, 결정론적 피로수명 점추정을 95% 신뢰구간으로 대체할 수 있는가?
3. 이 신뢰구간이 위험 기반 점검 주기 결정에 어떻게 활용되는가?

**연구 가설**:
H1. 동심 가정은 평균적으로 응력을 20–40% 과소평가 (Mattheck & Breloer 1994).
H2. 편심·종횡비 분포의 분산이 피로수명 신뢰구간의 폭을 결정한다.
H3. lower-bound S-N + 동심 가정의 조합은 편심 분포의 50–95 percentile 영역에 해당한다.

---

## 프로젝트 구조

```
fatigue-tree2/
├── CLAUDE.md                       ← 본 파일 (프로젝트 가이드)
├── README.md                       ← 운영 매뉴얼
├── docs/
│   ├── future_work.md              ← 1차 논문에서 도출된 11개 미래연구 (master 사본)
│   ├── research_plan.md            ← 본 프로젝트 상세 계획
│   ├── paper_draft_Springer.md     ← 1차 논문 (참고용, 수정 금지)
│   └── paper_draft_Elsevier.md     ← 1차 논문 (참고용, 수정 금지)
├── src/
│   ├── core/                       ← fatigue-tree에서 상속한 모듈 (수정 시 신중)
│   │   ├── wind_time_history.py
│   │   ├── weibull_fatigue.py
│   │   ├── fatigue_analysis.py
│   │   ├── wood_properties_db.py
│   │   ├── decay_effect.py         ← 동심 모델 (참고용; 새 모델로 대체 예정)
│   │   └── ...
│   ├── eccentric_decay/            ← 본 프로젝트 신규 코드
│   │   ├── geometry.py             ← 편심 단면 기하 모델
│   │   ├── stress_amplification.py ← 비동심 응력 증폭 계수
│   │   ├── monte_carlo.py          ← MC 샘플링 + 피로수명 분포
│   │   └── visualization.py
│   └── utils/
├── data/
│   ├── inherited/                  ← 1차 논문 데이터 (read-only)
│   │   ├── abaqus_inp/
│   │   ├── stress_history/
│   │   ├── wind/
│   │   └── parameters/
│   └── generated/                  ← MC 출력
├── output/
│   ├── figures/
│   ├── tables/
│   └── paper/                      ← 2차 논문 원고
├── notebooks/exploratory/          ← 탐색용 (재현 파이프라인은 scripts/)
└── scripts/                        ← 01_*, 02_* 순차 파이프라인
```

---

## 작업 규칙 (이 프로젝트 한정)

### 데이터 무결성
- `data/inherited/` 는 **read-only**. 직접 수정 금지. 필요 시 `data/generated/` 에 복사 후 변형.
- `docs/paper_draft_*.md` 는 1차 논문의 동결본. 본 프로젝트에서 *절대* 수정하지 않음. 1차 논문 수정이 필요하면 원본(`C:/ClaudeCode_Projects/fatigue-tree/output/02_paper/`) 에서 작업.

### 코드 작성 규칙
- 1차 논문과 동일한 SEED = 42, 동일 단위계 (응력 MPa, 길이 m, 시간 s)
- `src/core/` 모듈 변경 시 사유를 commit 메시지에 명시
- 새 코드는 모두 `src/eccentric_decay/` 에 위치
- 스크립트 명명: `01_*.py`, `02_*.py` … (글로벌 CLAUDE.md 규칙 준수)

### 재현성
- ABAQUS 라이선스 의존 부분은 별도 분리. 가능한 경우 inp 파일 + 추출된 stress history 만으로 재현 가능하도록 설계.
- Monte Carlo 결과는 `data/generated/` 에 seed 함께 저장.

### 모델 검증
- 새 편심 모델은 *반드시* e = 0 (동심) 한계에서 1차 논문 결과와 일치해야 함 → 회귀 테스트 필수.

---

## 1차 논문과의 관계

| 항목 | 1차 논문 (`fatigue-tree`) | 본 프로젝트 (`fatigue-tree2`) |
|------|--------------------------|-----------------------------|
| 부후 모델 | concentric hollow cylinder | eccentric, irregular, probabilistic |
| 결과 형태 | deterministic point estimate | distribution + 95% CI |
| 응력 증폭 | $\sigma_a = \sigma_0 / (1-(r_d/R)^4)$ | 편심 단면 모듈러스 + Monte Carlo |
| 인용 | — | 1차 논문 §4.5 직접 인용 |

---

## 목표 저널 (2차 논문)

후보 (영향력 + 적합도 순):
1. **Trees – Structure and Function** (Springer) — 1차 논문과 같은 저널 → 시리즈로 평가 가능
2. **Engineering Failure Analysis** (Elsevier) — 확률론적 손상 평가에 적합
3. **Urban Forestry & Urban Greening** (Elsevier) — 정책·관리 함의 강조 시

---

## 갱신 이력

- 2026-05-08: 초기 작성. fatigue-tree에서 분기 후 #2 (편심 부후 MC) 시작.
