#=============================================================
# GLaSS-SS Quick Start Example
# Simple Usage Examples
#=============================================================

# 1. Load required files
source("glass_ss_methods.R")

#=============================================================
# Example 1: Basic Usage - Variable Selection with Simple Data
#=============================================================

cat("\n=== Example 1: Basic GLaSS-SS Variable Selection ===\n")

# Generate data (n=100 samples, SNR=1)
set.seed(123)
data <- generate_group_data(n = 100, snr = 1, half = 0)

# Split data (train/test)
n <- nrow(data$X)
train_idx <- sample(1:n, size = floor(0.7 * n))
X_train <- data$X[train_idx, ]
y_train <- data$y[train_idx]

# Define group structure (3 groups, 40 variables each)
groups <- rep(1:3, each = 40)

# Create Laplacian matrix
L <- create_laplacian_matrix(X_train)

cat("Data size:", nrow(X_train), "x", ncol(X_train), "\n")
cat("Response distribution:", table(y_train), "\n")
cat("True signal variables:", sum(data$beta_true != 0), "\n\n")

# Run GLaSS-SS (using single PFER value for quick execution)
cat("Running GLaSS-SS...\n")
result <- my_stabsel_parallel(
  x = X_train,
  y = y_train,
  fitfun = my_adaptive_fitfun,
  args.fitfun = list(
    groups = groups,
    L_sparse = L,
    alpha_seq = seq(0, 1, 0.25),  # Mixing ratio of Lasso and Laplacian
    nlambda = 10
  ),
  cutoff = 0.6,           # Selection threshold
  pfer_values = c(5),     # Per-Family Error Rate
  sampling.type = "SS",   # Stability Selection
  B = 20,                 # Number of subsampling iterations (for quick run)
  mc.cores = 4,           # Number of cores for parallel processing
  verbose = TRUE
)

# Check results
selected_vars <- result[["5"]]$selected
cat("\nSelected variables:", length(selected_vars), "\n")
cat("Variable indices:", sort(selected_vars), "\n")

# Check True Positives
true_signals <- which(data$beta_true != 0)
tp <- length(intersect(selected_vars, true_signals))
fp <- length(setdiff(selected_vars, true_signals))
cat("\nPerformance:\n")
cat("  True Positives:", tp, "/", length(true_signals), "\n")
cat("  False Positives:", fp, "\n")
cat("  TPR (Sensitivity):", round(tp / length(true_signals), 3), "\n")

#=============================================================
# Example 2: Comparison with Other Methods
#=============================================================

cat("\n\n=== Example 2: Comparison with Baseline Methods ===\n")

# Lasso with Stability Selection
cat("\nRunning Lasso...\n")
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
cat("Running Elastic Net (CV)...\n")
library(glmnet)
set.seed(123)
en_cv <- cv.glmnet(X_train, y_train, family = "binomial", alpha = 0.5)
en_coef <- coef(en_cv, s = "lambda.min")
en_selected <- which(as.numeric(en_coef[-1]) != 0)

# Compare results
cat("\n=== Variable Selection Results Comparison ===\n")
cat("GLaSS-SS:      ", length(selected_vars), "variables selected\n")
cat("Lasso (SS):    ", length(lasso_result$selected), "variables selected\n")
cat("Elastic Net:   ", length(en_selected), "variables selected\n")

# Calculate TPR
glass_tpr <- length(intersect(selected_vars, true_signals)) / length(true_signals)
lasso_tpr <- length(intersect(lasso_result$selected, true_signals)) / length(true_signals)
en_tpr <- length(intersect(en_selected, true_signals)) / length(true_signals)

cat("\n=== True Positive Rate (TPR) Comparison ===\n")
cat("GLaSS-SS:      ", round(glass_tpr, 3), "\n")
cat("Lasso (SS):    ", round(lasso_tpr, 3), "\n")
cat("Elastic Net:   ", round(en_tpr, 3), "\n")

#=============================================================
# Example 3: Selection Probability Visualization
#=============================================================

cat("\n\n=== Example 3: Selection Probability Visualization ===\n")

# Extract selection probabilities
selection_probs <- result[["5"]]$phat
top_vars <- order(selection_probs, decreasing = TRUE)[1:20]

# Simple text-based visualization
cat("\nSelection probabilities for top 20 variables:\n")
for (i in 1:20) {
  var_idx <- top_vars[i]
  prob <- selection_probs[var_idx]
  is_true_signal <- var_idx %in% true_signals

  # Bar chart (text version)
  bar_length <- round(prob * 40)
  bar <- paste(rep("█", bar_length), collapse = "")

  cat(sprintf("Var %3d: %s %.3f %s\n",
              var_idx,
              bar,
              prob,
              ifelse(is_true_signal, "✓ (true signal)", "")))
}

cat("\nExecution complete! See simulation_study.R for more detailed analysis.\n")
