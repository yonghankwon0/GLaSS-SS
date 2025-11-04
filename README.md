# GLaSS-SS: Group-Laplacian Structured Shrinkage with Stability Selection

[![R](https://img.shields.io/badge/R-276DC3?style=flat&logo=r&logoColor=white)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

A robust variable selection method that simultaneously leverages group structure and correlation patterns in high-dimensional data

## 📋 Table of Contents

- [Overview](#overview)
- [Key Features](#key-features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [File Structure](#file-structure)
- [Methodology](#methodology)
- [Performance Comparison](#performance-comparison)
- [Simulation Studies](#simulation-studies)
- [References](#references)

## Overview

GLaSS-SS is a powerful variable selection framework for high-dimensional biological data analysis. It simultaneously leverages **group structure** and **correlation patterns** among variables to provide more stable and interpretable results compared to existing methods.

### Core Features

- **Adaptive Penalty Optimization**: Combines Group Lasso and Laplacian regularization
- **Stability Selection**: Controls False Discovery Rate through Stability Selection
- **Parallel Processing**: Efficient computation for large-scale datasets
- **Flexible Group Structure**: Supports various group definitions

## Key Features

### 1. Graph Structure Learning
```r
# Create optimal adjacency matrix using GLASSO
adj_matrix <- create_optimal_adjacency_glasso(X)

# Generate Laplacian matrix based on correlation structure
L <- create_laplacian_matrix(X)
```

### 2. Adaptive Penalty Optimization
```r
# Fit GLaSS model
result <- fit_adaptive_penalty_optim_noCV(
  X = X_train,
  y = y_train,
  L_sparse = L,
  groups = groups,
  alpha_seq = seq(0, 1, 0.25),  # Lasso-Laplacian mixing ratio
  nlambda = 10
)
```

### 3. Stability Selection
```r
# Parallel stability selection
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

## Installation

### Install Required Packages

```r
# CRAN packages
install.packages(c(
  "Matrix", "Rcpp", "glmnet", "grpreg",
  "huge", "stabs", "parallel",
  "tidyverse", "mvtnorm", "corrplot",
  "pROC", "randomForest", "ggplot2", "gridExtra"
))
```

### Download GLaSS-SS

```bash
git clone https://github.com/yonghankwon0/GLaSS-SS.git
cd GLaSS-SS
```

## Quick Start

### Run Simple Example

```r
# Run example script
source("example_quick_start.R")
```

This script includes:
1. **Basic Usage**: Run GLaSS-SS on simulated data
2. **Method Comparison**: Compare performance with Lasso and Elastic Net
3. **Visualization**: Check selection probability distribution

### Manual Execution

```r
# 1. Load methods
source("glass_ss_methods.R")

# 2. Generate data
set.seed(123)
data <- generate_group_data(n = 100, snr = 1, half = 0)

# 3. Prepare training data
train_idx <- sample(1:nrow(data$X), 70)
X_train <- data$X[train_idx, ]
y_train <- data$y[train_idx]

# 4. Define group and graph structure
groups <- rep(1:3, each = 40)
L <- create_laplacian_matrix(X_train)

# 5. Run GLaSS-SS
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

# 6. Check results
selected_vars <- result[["5"]]$selected
cat("Selected variables:", length(selected_vars), "\n")
```

## File Structure

```
GLaSS-SS/
├── glass_ss_methods.R          # Core method implementation
│   ├── Graph structure functions
│   ├── GLaSS optimization algorithm
│   ├── Stability selection integration
│   └── Baseline methods (Lasso, Elastic Net, Group Lasso)
│
├── simulation_study.R          # Simulation and performance evaluation
│   ├── Data generation functions
│   ├── Stability metrics calculation
│   ├── Performance evaluation framework
│   └── Visualization functions
│
├── example_quick_start.R       # Quick start example
│   ├── Basic usage
│   ├── Method comparison
│   └── Result visualization
│
└── README.md                   # Project documentation
```

## Methodology

### GLaSS-SS Objective Function

GLaSS-SS minimizes the following objective function:

```
minimize: L(β) + λ·α·Σ√(|G_j|)||β_G_j||_2 + λ·(1-α)·β^T L β
```

Where:
- `L(β)`: Logistic loss function
- `α ∈ [0,1]`: Mixing parameter between Group Lasso and Laplacian regularization
- `G_j`: Variable indices in group j
- `L`: Laplacian matrix (encoding correlation structure among variables)

### Optimization Algorithm

Uses **Generalized Forward-Backward Splitting (GFBS)** algorithm:

1. **Initialization**: Set β⁰ = 0, two auxiliary variables z_g⁰, z_l⁰
2. **Iteration**:
   - Gradient step: Compute gradient of logistic loss
   - Proximal operator (Group Lasso): Group-wise soft-thresholding
   - Proximal operator (Laplacian): Solve quadratic regularization
   - Update variables and check convergence

### Stability Selection

**Stability Selection** framework for controlling False Discovery:

1. **Subsampling**: Repeatedly subsample data B times (typically B=50-100)
2. **Variable Selection**: Run GLaSS on each subsample
3. **Selection Probability**: P_hat(j) = (# times variable j selected) / B
4. **Final Selection**: Select variables with P_hat(j) > cutoff (typically cutoff=0.6)

**Per-Family Error Rate (PFER) Control**:
```
E[FP] ≤ (q²)/(πcutoff - 0.5) ≤ PFER
```

## Performance Comparison

### Baseline Methods

1. **Lasso + Stability Selection**: Basic L1 regularization
2. **Elastic Net + Stability Selection**: L1 + L2 mixture
3. **Group Lasso + Stability Selection**: Leverages group structure only
4. **Elastic Net (CV)**: Cross-validation based
5. **Group Lasso (CV)**: Cross-validation based

### Evaluation Metrics

- **True Positive Rate (TPR)**: Detection rate of true signal variables
- **Positive Predictive Value (PPV)**: Proportion of true signals among selected variables
- **F1 Score**: Harmonic mean of TPR and PPV
- **AUC**: Prediction performance of Random Forest model
- **Nogueira Stability**: Measure of selection stability
- **Jaccard Stability**: Consistency of selections across repetitions

## Simulation Studies

### Simulation Scenarios

```r
# Run simulation_study.R
source("simulation_study.R")
```

Test scenarios:
- **Sample Size**: n ∈ {60, 120}
- **SNR**: Signal-to-Noise Ratio ∈ {1}
- **Signal Patterns**:
  - Block signal (3 groups)
  - Sparse signal (30% active in groups)
- **Group Structures**:
  - Well-specified: 6 groups × 20 variables
  - Misspecified: 3 groups × 40 variables
  - Fine-grained: 12 groups × 10 variables

### Output Files

Files generated when running simulations:
```
simulation_results_YYYY-MM-DD_HH-MM-SS/
├── *.txt                       # Detailed result logs
├── sim_result_*.rds            # Saved R objects
└── all_plots_*.pdf             # Performance comparison plots
```

## References

### Theoretical Background

1. **Stability Selection**
   Meinshausen, N., & Bühlmann, P. (2010). Stability selection. *Journal of the Royal Statistical Society: Series B*, 72(4), 417-473.

2. **Group Lasso**
   Yuan, M., & Lin, Y. (2006). Model selection and estimation in regression with grouped variables. *Journal of the Royal Statistical Society: Series B*, 68(1), 49-67.

3. **Graph-Guided Fused Lasso**
   Kim, S., & Xing, E. P. (2009). Statistical estimation of correlated genome associations to a quantitative trait network. *PLoS Genetics*, 5(8), e1000587.

4. **Forward-Backward Splitting**
   Combettes, P. L., & Pesquet, J. C. (2011). Proximal splitting methods in signal processing. In *Fixed-point algorithms for inverse problems in science and engineering* (pp. 185-212). Springer.

### Related Packages

- **stabs**: Stability selection implementation
- **glmnet**: Lasso and Elastic Net
- **grpreg**: Group regularization
- **huge**: High-dimensional graph estimation

## Contributing

Please submit issues and suggestions to [GitHub Issues](https://github.com/yonghankwon0/GLaSS-SS/issues).

## License

This project is distributed under the MIT License.
