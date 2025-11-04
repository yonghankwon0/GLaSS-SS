#=============================================================
# GLaSS-SS Quick Start Example
# 간단한 사용 예제
#=============================================================

# 1. 필요한 파일 로드
source("glass_ss_methods.R")

#=============================================================
# 예제 1: 기본 사용법 - 간단한 데이터로 변수 선택
#=============================================================

cat("\n=== 예제 1: 기본 GLaSS-SS 변수 선택 ===\n")

# 데이터 생성 (n=100 샘플, SNR=1)
set.seed(123)
data <- generate_group_data(n = 100, snr = 1, half = 0)

# 데이터 분할 (훈련/테스트)
n <- nrow(data$X)
train_idx <- sample(1:n, size = floor(0.7 * n))
X_train <- data$X[train_idx, ]
y_train <- data$y[train_idx]

# 그룹 구조 정의 (3개 그룹, 각 40개 변수)
groups <- rep(1:3, each = 40)

# 라플라시안 행렬 생성
L <- create_laplacian_matrix(X_train)

cat("데이터 크기:", nrow(X_train), "x", ncol(X_train), "\n")
cat("응답 변수 분포:", table(y_train), "\n")
cat("실제 신호 변수 개수:", sum(data$beta_true != 0), "\n\n")

# GLaSS-SS 실행 (빠른 실행을 위해 단일 PFER 값 사용)
cat("GLaSS-SS 실행 중...\n")
result <- my_stabsel_parallel(
  x = X_train,
  y = y_train,
  fitfun = my_adaptive_fitfun,
  args.fitfun = list(
    groups = groups,
    L_sparse = L,
    alpha_seq = seq(0, 1, 0.25),  # Lasso와 Laplacian의 혼합 비율
    nlambda = 10
  ),
  cutoff = 0.6,           # 선택 임계값
  pfer_values = c(5),     # Per-Family Error Rate
  sampling.type = "SS",   # Stability Selection
  B = 20,                 # 서브샘플링 횟수 (빠른 실행)
  mc.cores = 4,           # 병렬 처리 코어 수
  verbose = TRUE
)

# 결과 확인
selected_vars <- result[["5"]]$selected
cat("\n선택된 변수:", length(selected_vars), "개\n")
cat("변수 인덱스:", sort(selected_vars), "\n")

# 진양성(True Positive) 확인
true_signals <- which(data$beta_true != 0)
tp <- length(intersect(selected_vars, true_signals))
fp <- length(setdiff(selected_vars, true_signals))
cat("\n성능:\n")
cat("  True Positives:", tp, "/", length(true_signals), "\n")
cat("  False Positives:", fp, "\n")
cat("  TPR (민감도):", round(tp / length(true_signals), 3), "\n")

#=============================================================
# 예제 2: 다른 방법들과 비교
#=============================================================

cat("\n\n=== 예제 2: 기준 방법과 비교 ===\n")

# Lasso with Stability Selection
cat("\nLasso 실행 중...\n")
lasso_result <- stabsel(
  x = X_train,
  y = y_train,
  fitfun = glmnet.lasso,
  cutoff = 0.6,
  PFER = 5,
  sampling.type = "SS",
  B = 20
)

# Elastic Net with CV
cat("Elastic Net (CV) 실행 중...\n")
library(glmnet)
set.seed(123)
en_cv <- cv.glmnet(X_train, y_train, family = "binomial", alpha = 0.5)
en_coef <- coef(en_cv, s = "lambda.min")
en_selected <- which(as.numeric(en_coef[-1]) != 0)

# 결과 비교
cat("\n=== 변수 선택 결과 비교 ===\n")
cat("GLaSS-SS:      ", length(selected_vars), "개 변수 선택\n")
cat("Lasso (SS):    ", length(lasso_result$selected), "개 변수 선택\n")
cat("Elastic Net:   ", length(en_selected), "개 변수 선택\n")

# TPR 계산
glass_tpr <- length(intersect(selected_vars, true_signals)) / length(true_signals)
lasso_tpr <- length(intersect(lasso_result$selected, true_signals)) / length(true_signals)
en_tpr <- length(intersect(en_selected, true_signals)) / length(true_signals)

cat("\n=== True Positive Rate (TPR) 비교 ===\n")
cat("GLaSS-SS:      ", round(glass_tpr, 3), "\n")
cat("Lasso (SS):    ", round(lasso_tpr, 3), "\n")
cat("Elastic Net:   ", round(en_tpr, 3), "\n")

#=============================================================
# 예제 3: 선택 확률 시각화
#=============================================================

cat("\n\n=== 예제 3: 선택 확률 시각화 ===\n")

# 선택 확률 추출
selection_probs <- result[["5"]]$phat
top_vars <- order(selection_probs, decreasing = TRUE)[1:20]

# 간단한 텍스트 기반 시각화
cat("\nTop 20 변수의 선택 확률:\n")
for (i in 1:20) {
  var_idx <- top_vars[i]
  prob <- selection_probs[var_idx]
  is_true_signal <- var_idx %in% true_signals

  # 막대 그래프 (텍스트 버전)
  bar_length <- round(prob * 40)
  bar <- paste(rep("█", bar_length), collapse = "")

  cat(sprintf("변수 %3d: %s %.3f %s\n",
              var_idx,
              bar,
              prob,
              ifelse(is_true_signal, "✓ (실제 신호)", "")))
}

cat("\n실행 완료! 더 자세한 분석은 simulation_study.R을 참고하세요.\n")
