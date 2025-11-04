# GLaSS-SS: Group-Laplacian Structured Shrinkage with Stability Selection

[![R](https://img.shields.io/badge/R-276DC3?style=flat&logo=r&logoColor=white)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

고차원 데이터에서 그룹 구조와 상관관계를 동시에 고려한 안정적인 변수 선택 방법

## 📋 목차

- [개요](#개요)
- [주요 기능](#주요-기능)
- [설치 방법](#설치-방법)
- [빠른 시작](#빠른-시작)
- [파일 구조](#파일-구조)
- [방법론](#방법론)
- [성능 비교](#성능-비교)
- [시뮬레이션 연구](#시뮬레이션-연구)
- [참고 문헌](#참고-문헌)

## 개요

GLaSS-SS는 고차원 생물학적 데이터 분석을 위한 강력한 변수 선택 프레임워크입니다. 변수 간 **그룹 구조**와 **상관관계**를 동시에 활용하여 기존 방법보다 안정적이고 해석 가능한 결과를 제공합니다.

### 핵심 특징

- **적응적 페널티 최적화**: Group Lasso와 Laplacian 정규화를 혼합
- **안정성 선택**: Stability Selection으로 False Discovery 제어
- **병렬 처리**: 대규모 데이터셋에서 효율적인 연산
- **유연한 그룹 구조**: 다양한 그룹 정의 지원

## 주요 기능

### 1. 그래프 구조 학습
```r
# GLASSO로 최적 인접 행렬 생성
adj_matrix <- create_optimal_adjacency_glasso(X)

# 상관관계 기반 라플라시안 행렬 생성
L <- create_laplacian_matrix(X)
```

### 2. 적응적 페널티 최적화
```r
# GLaSS 모델 적합
result <- fit_adaptive_penalty_optim_noCV(
  X = X_train,
  y = y_train,
  L_sparse = L,
  groups = groups,
  alpha_seq = seq(0, 1, 0.25),  # Lasso-Laplacian 혼합 비율
  nlambda = 10
)
```

### 3. 안정성 선택
```r
# 병렬 안정성 선택
stabsel_result <- my_stabsel_parallel(
  x = X_train,
  y = y_train,
  fitfun = my_adaptive_fitfun,
  args.fitfun = list(groups = groups, L_sparse = L),
  cutoff = 0.6,
  pfer_values = c(1, 5, 10),
  B = 50,
  mc.cores = 10
)
```

## 설치 방법

### 필수 패키지 설치

```r
# CRAN 패키지
install.packages(c(
  "Matrix", "Rcpp", "glmnet", "grpreg",
  "huge", "stabs", "parallel",
  "tidyverse", "mvtnorm", "corrplot",
  "pROC", "randomForest", "ggplot2", "gridExtra"
))
```

### GLaSS-SS 다운로드

```bash
git clone https://github.com/yonghankwon0/GLaSS-SS.git
cd GLaSS-SS
```

## 빠른 시작

### 간단한 예제 실행

```r
# 예제 스크립트 실행
source("example_quick_start.R")
```

이 스크립트는 다음을 포함합니다:
1. **기본 사용법**: 시뮬레이션 데이터로 GLaSS-SS 실행
2. **방법 비교**: Lasso, Elastic Net과 성능 비교
3. **시각화**: 선택 확률 분포 확인

### 직접 실행하기

```r
# 1. 메서드 로드
source("glass_ss_methods.R")

# 2. 데이터 생성
set.seed(123)
data <- generate_group_data(n = 100, snr = 1, half = 0)

# 3. 훈련 데이터 준비
train_idx <- sample(1:nrow(data$X), 70)
X_train <- data$X[train_idx, ]
y_train <- data$y[train_idx]

# 4. 그룹 및 그래프 구조 정의
groups <- rep(1:3, each = 40)
L <- create_laplacian_matrix(X_train)

# 5. GLaSS-SS 실행
result <- my_stabsel_parallel(
  x = X_train,
  y = y_train,
  fitfun = my_adaptive_fitfun,
  args.fitfun = list(
    groups = groups,
    L_sparse = L,
    alpha_seq = seq(0, 1, 0.25),
    nlambda = 10
  ),
  cutoff = 0.6,
  pfer_values = c(5),
  sampling.type = "SS",
  B = 20,
  mc.cores = 4
)

# 6. 결과 확인
selected_vars <- result[["5"]]$selected
cat("선택된 변수:", length(selected_vars), "개\n")
```

## 파일 구조

```
GLaSS-SS/
├── glass_ss_methods.R          # 핵심 메서드 구현
│   ├── 그래프 구조 함수
│   ├── GLaSS 최적화 알고리즘
│   ├── 안정성 선택 통합
│   └── 기준 방법 (Lasso, Elastic Net, Group Lasso)
│
├── simulation_study.R          # 시뮬레이션 및 성능 평가
│   ├── 데이터 생성 함수
│   ├── 안정성 메트릭 계산
│   ├── 성능 평가 프레임워크
│   └── 시각화 함수
│
├── example_quick_start.R       # 빠른 시작 예제
│   ├── 기본 사용법
│   ├── 방법 비교
│   └── 결과 시각화
│
└── README.md                   # 프로젝트 문서
```

## 방법론

### GLaSS-SS 목적 함수

GLaSS-SS는 다음 목적 함수를 최소화합니다:

```
minimize: L(β) + λ·α·Σ√(|G_j|)||β_G_j||_2 + λ·(1-α)·β^T L β
```

여기서:
- `L(β)`: 로지스틱 손실 함수
- `α ∈ [0,1]`: Group Lasso와 Laplacian 정규화의 혼합 비율
- `G_j`: j번째 그룹의 변수 인덱스
- `L`: 라플라시안 행렬 (변수 간 상관관계 인코딩)

### 최적화 알고리즘

**Generalized Forward-Backward Splitting (GFBS)** 알고리즘 사용:

1. **초기화**: β⁰ = 0, 두 개의 보조 변수 z_g⁰, z_l⁰ 설정
2. **반복**:
   - Gradient step: 로지스틱 손실의 그래디언트 계산
   - Proximal operator (Group Lasso): 그룹별 soft-thresholding
   - Proximal operator (Laplacian): 2차 정규화 해결
   - 변수 업데이트 및 수렴 확인

### 안정성 선택

**Stability Selection** 프레임워크로 False Discovery 제어:

1. **서브샘플링**: 데이터를 B번 반복 서브샘플링 (보통 B=50-100)
2. **변수 선택**: 각 서브샘플에서 GLaSS 실행
3. **선택 확률 계산**: P_hat(j) = (변수 j가 선택된 횟수) / B
4. **최종 선택**: P_hat(j) > cutoff인 변수만 선택 (보통 cutoff=0.6)

**Per-Family Error Rate (PFER) 제어**:
```
E[FP] ≤ (q²)/(πcutoff - 0.5) ≤ PFER
```

## 성능 비교

### 비교 대상 방법

1. **Lasso + Stability Selection**: 기본 L1 정규화
2. **Elastic Net + Stability Selection**: L1 + L2 혼합
3. **Group Lasso + Stability Selection**: 그룹 구조만 활용
4. **Elastic Net (CV)**: Cross-validation 기반
5. **Group Lasso (CV)**: Cross-validation 기반

### 평가 지표

- **True Positive Rate (TPR)**: 실제 신호 변수의 탐지율
- **Positive Predictive Value (PPV)**: 선택된 변수 중 실제 신호 비율
- **F1 Score**: TPR과 PPV의 조화 평균
- **AUC**: Random Forest 모델의 예측 성능
- **Nogueira Stability**: 선택의 안정성 측정
- **Jaccard Stability**: 반복 간 선택 일관성

## 시뮬레이션 연구

### 시뮬레이션 시나리오

```r
# simulation_study.R 실행
source("simulation_study.R")
```

테스트 시나리오:
- **샘플 크기**: n ∈ {60, 120}
- **SNR**: Signal-to-Noise Ratio ∈ {1}
- **신호 패턴**:
  - Block signal (3 groups)
  - Sparse signal (30% active in groups)
- **그룹 구조**:
  - Well-specified: 6 groups × 20 variables
  - Misspecified: 3 groups × 40 variables
  - Fine-grained: 12 groups × 10 variables

### 결과 출력

시뮬레이션 실행 시 생성되는 파일:
```
simulation_results_YYYY-MM-DD_HH-MM-SS/
├── *.txt                       # 상세 결과 로그
├── sim_result_*.rds            # R 객체 저장
└── all_plots_*.pdf             # 성능 비교 그래프
```

## 참고 문헌

### 이론적 배경

1. **Stability Selection**
   Meinshausen, N., & Bühlmann, P. (2010). Stability selection. *Journal of the Royal Statistical Society: Series B*, 72(4), 417-473.

2. **Group Lasso**
   Yuan, M., & Lin, Y. (2006). Model selection and estimation in regression with grouped variables. *Journal of the Royal Statistical Society: Series B*, 68(1), 49-67.

3. **Graph-Guided Fused Lasso**
   Kim, S., & Xing, E. P. (2009). Statistical estimation of correlated genome associations to a quantitative trait network. *PLoS Genetics*, 5(8), e1000587.

4. **Forward-Backward Splitting**
   Combettes, P. L., & Pesquet, J. C. (2011). Proximal splitting methods in signal processing. In *Fixed-point algorithms for inverse problems in science and engineering* (pp. 185-212). Springer.

### 관련 패키지

- **stabs**: Stability selection implementation
- **glmnet**: Lasso and Elastic Net
- **grpreg**: Group regularization
- **huge**: High-dimensional graph estimation

## 문의 및 기여

이슈 및 제안사항은 [GitHub Issues](https://github.com/yonghankwon0/GLaSS-SS/issues)에 등록해주세요.

## 라이선스

이 프로젝트는 MIT 라이선스 하에 배포됩니다.
