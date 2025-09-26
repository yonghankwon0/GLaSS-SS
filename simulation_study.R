#=============================================================
# Simulation Study for GLaSS-SS
# Performance Evaluation and Comparison
#=============================================================

# Source the methods file
source("glass_ss_methods.R")

#=============================================================
# Additional Required Packages for Simulation
#=============================================================
library(tidyverse)
library(mvtnorm)
library(corrplot)
library(pROC)
library(randomForest)
library(ggplot2)
library(gridExtra)

#=============================================================
# 1. Stability Metrics
#=============================================================

#' Calculate Nogueira stability measure with confidence intervals
#' @param Z Binary selection matrix (rows: repetitions, cols: variables)
#' @return List with stability, variance, and confidence intervals
calculate_stability <- function(Z) {
  M <- nrow(Z)
  d <- ncol(Z)
  k_bar <- mean(rowSums(Z))
  p_hat <- colMeans(Z)
  s2 <- (M/(M-1)) * p_hat * (1 - p_hat)
  denominator <- (k_bar/d) * (1 - k_bar/d)
  stability <- 1 - sum(s2)/(d * denominator)
  
  # Calculate variance estimate for confidence intervals
  phi_i <- sapply(1:M, function(i) {
    ki <- sum(Z[i,])
    numerator <- sum(Z[i,] * p_hat)/d - ki * k_bar/d^2 +
      stability^2 * (2*k_bar*ki/d^2 - ki/d - k_bar/d + 1)
    numerator / denominator
  })
  phi_bar <- mean(phi_i)
  var_est <- 4/M^2 * sum((phi_i - phi_bar)^2)
  ci <- stability + c(-1,1) * qnorm(0.975) * sqrt(var_est)
  
  list(
    stability = stability,
    variance = var_est,
    ci_lower = ci[1],
    ci_upper = ci[2]
  )
}

#' Calculate Jaccard stability index
#' @param selection_matrix Binary selection matrix
#' @return List with mean Jaccard index and confidence intervals
calculate_jaccard_stability <- function(selection_matrix) {
  n_reps <- nrow(selection_matrix)
  
  if (n_reps < 2) {
    return(list(stability = NA, ci_lower = NA, ci_upper = NA, sd = NA))
  }
  
  jaccard_scores <- numeric()
  
  # Calculate pairwise Jaccard similarities
  for (i in 1:(n_reps-1)) {
    for (j in (i+1):n_reps) {
      set1 <- which(selection_matrix[i, ] == 1)
      set2 <- which(selection_matrix[j, ] == 1)
      
      if (length(set1) == 0 && length(set2) == 0) {
        jaccard <- 1  # Both empty sets
      } else {
        intersection <- length(intersect(set1, set2))
        union <- length(union(set1, set2))
        jaccard <- intersection / union
      }
      
      jaccard_scores <- c(jaccard_scores, jaccard)
    }
  }
  
  # Calculate mean and confidence interval
  mean_jaccard <- mean(jaccard_scores, na.rm = TRUE)
  sd_jaccard <- sd(jaccard_scores, na.rm = TRUE)
  
  if (length(jaccard_scores) > 1) {
    se <- sd_jaccard / sqrt(length(jaccard_scores))
    ci_lower <- mean_jaccard - 1.96 * se
    ci_upper <- mean_jaccard + 1.96 * se
  } else {
    ci_lower <- NA
    ci_upper <- NA
  }
  
  return(list(
    stability = mean_jaccard,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    sd = sd_jaccard
  ))
}

#' Find maximum valid PFER (Per-Family Error Rate) value
#' @param data_matrix Input data matrix
#' @param start_pfer Starting PFER value to test
#' @return Maximum valid PFER value
find_max_valid_pfer <- function(data_matrix, start_pfer = 1) {
  max_valid <- 1
  for (pfer in start_pfer:1000) {
    tryCatch({
      res <- stabsel_parameters(p = ncol(data_matrix),
                                cutoff = 0.6,
                                PFER = pfer)
      cat(sprintf("Testing PFER = %d ... Valid\n", pfer))
      max_valid <- pfer
    }, error = function(e) {
      return(max_valid)
    })
  }
  return(max_valid)
}

#=============================================================
# 2. Data Generation
#=============================================================

#' Generate simulated group-structured data
#' @param n Number of observations
#' @param snr Signal-to-noise ratio
#' @param half Signal pattern type:
#'   0: 3 groups (40 vars each) - block signal
#'   1: 6 groups (20 vars each) - block signal
#'   2: 3 groups with opposite signs within groups
#'   3: 6 groups with sparse signals (30% active)
#' @return List with X, y, true beta, covariance structure, and SNR
generate_group_data <- function(n, snr, half = 0) {
  p <- 120           # Total number of predictors
  k <- 3             # Number of major blocks
  group_size <- p / k
  
  # Function to create complex within-group covariance structure
  create_complex_sigma <- function(size) {
    sigma_mat <- matrix(0, nrow = size, ncol = size)
    
    # 1. Exponential decay correlation
    for(i in 1:size) {
      for(j in 1:size) {
        sigma_mat[i, j] <- 0.95 ^ abs(i - j)
      }
    }
    
    # 2. Add periodic patterns
    for(i in 1:size) {
      for(j in 1:size) {
        if(abs(i - j) %% 3 == 0) {
          sigma_mat[i, j] <- sigma_mat[i, j] + 0.15 * cos(2 * pi * (i + j) / size)
        }
      }
    }
    
    # 3. Add cluster structure (4 clusters within each group)
    cluster_size <- size / 4
    for(c in 1:4) {
      idx <- (((c - 1) * cluster_size) + 1):(c * cluster_size)
      for(i in idx) {
        for(j in idx) {
          if(i != j) {
            sigma_mat[i, j] <- sigma_mat[i, j] + 0.2
          }
        }
      }
    }
    
    # 4. Add random strong correlations
    random_pairs <- matrix(sample(1:size, size = floor(size/3) * 2, replace = FALSE), ncol = 2)
    for(i in 1:nrow(random_pairs)) {
      sigma_mat[random_pairs[i, 1], random_pairs[i, 2]] <- sigma_mat[random_pairs[i, 1], random_pairs[i, 2]] + 0.3
      sigma_mat[random_pairs[i, 2], random_pairs[i, 1]] <- sigma_mat[random_pairs[i, 2], random_pairs[i, 1]] + 0.3
    }
    
    # Symmetrize and ensure positive definiteness
    sigma_mat <- (sigma_mat + t(sigma_mat)) / 2
    eigen_decomp <- eigen(sigma_mat)
    min_eigen <- min(eigen_decomp$values)
    if(min_eigen < 0) {
      sigma_mat <- sigma_mat + (-min_eigen + 0.01) * diag(size)
    }
    
    # Normalize to correlation matrix
    D <- diag(1/sqrt(diag(sigma_mat)))
    sigma_mat <- D %*% sigma_mat %*% D
    return(sigma_mat)
  }
  
  # Create block-diagonal covariance matrix
  Sigma_blocks <- lapply(1:k, function(i) create_complex_sigma(group_size))
  Sigma <- Matrix::bdiag(Sigma_blocks)
  Sigma <- as.matrix(Sigma)
  
  # Add cross-group correlations
  for(i in 1:k) {
    for(j in 1:k) {
      if(i != j) {
        block_i <- ((i - 1) * group_size + 1):(i * group_size)
        block_j <- ((j - 1) * group_size + 1):(j * group_size)
        cross_corr <- matrix(runif(group_size^2, -0.1, 0.1), group_size, group_size)
        Sigma[block_i, block_j] <- cross_corr
        Sigma[block_j, block_i] <- t(cross_corr)
      }
    }
  }
  
  # Special case for half==2: negative correlations between sub-groups
  if(half == 2) {
    for(g in c(1, 3)) {
      block_idx <- ((g - 1) * group_size + 1):(g * group_size)
      sub1 <- block_idx[1:(length(block_idx)/2)]
      sub2 <- block_idx[(length(block_idx)/2 + 1):length(block_idx)]
      for(i in sub1) {
        for(j in sub2) {
          Sigma[i, j] <- -abs(Sigma[i, j])
          Sigma[j, i] <- Sigma[i, j]
        }
      }
    }
  }
  
  # Final symmetrization and positive definiteness check
  Sigma <- (Sigma + t(Sigma)) / 2
  eigen_decomp <- eigen(Sigma)
  min_eigen <- min(eigen_decomp$values)
  if(min_eigen < 0) {
    Sigma <- Sigma + (-min_eigen + 0.01) * diag(p)
  }
  Sigma <- Sigma * 3
  diag(Sigma) <- 1
  
  # Generate predictor matrix
  X <- mvtnorm::rmvnorm(n, sigma = Sigma)
  
  # Generate true coefficients based on pattern type
  if (half == 0) {
    # Block pattern: 3 groups (40 vars each)
    beta_dagger <- c(
      rnorm(40, 1.5, 0.5),  # Group 1: signal
      rnorm(40, 0, 0),      # Group 2: noise
      rnorm(40, 3, 0.5)     # Group 3: signal
    )
  } else if (half == 1) {
    # Block pattern: 6 groups (20 vars each)
    beta_dagger <- c(
      rnorm(20, 1.5, 0.5),  # Group 1: signal
      rnorm(20, 0, 0),      # Group 2: noise
      rnorm(20, 0, 0.5),    # Group 3: weak signal
      rnorm(20, 0, 0),      # Group 4: noise
      rnorm(20, 3, 0.5),    # Group 5: signal
      rnorm(20, 0, 0)       # Group 6: noise
    )
  } else if (half == 2) {
    # Opposite signs within signal groups
    beta_group1 <- c(rnorm(20, 1.5, 0.5), rnorm(20, -1.5, 0.5))
    beta_group2 <- rep(0, 40)
    beta_group3 <- c(rnorm(20, 3, 0.5), rnorm(20, -3, 0.5))
    beta_dagger <- c(beta_group1, beta_group2, beta_group3)
  } else if (half == 3) {
    # Sparse signals: only 30% of signal group variables are active
    create_sparse_signal_group <- function(size, mean, sd, signal_prop = 0.3) {
      n_signal <- round(size * signal_prop)
      signal_idx <- sample(1:size, n_signal)
      beta <- rep(0, size)
      beta[signal_idx] <- rnorm(n_signal, mean, sd)
      return(beta)
    }
    
    beta_dagger <- c(
      create_sparse_signal_group(20, 1.5, 0.5),
      rep(0, 20),
      create_sparse_signal_group(20, 0, 0.5),
      rep(0, 20),
      create_sparse_signal_group(20, 3, 0.5),
      rep(0, 20)
    )
  }
  
  # Calculate signal variance and noise level
  signal_variance <- as.numeric(t(beta_dagger) %*% Sigma %*% beta_dagger)
  sigma_final <- sqrt(signal_variance / snr)
  
  # Generate response
  epsilon <- rnorm(n, mean = 0, sd = sigma_final)
  linear_pred <- X %*% beta_dagger + epsilon
  probs <- 1 / (1 + exp(-linear_pred))
  y <- rbinom(n, size = 1, prob = probs)
  
  # Calculate realized SNR
  realized_snr <- signal_variance / (sigma_final^2)
  
  return(list(
    X = X, 
    y = y, 
    beta_true = beta_dagger, 
    Sigma = Sigma, 
    sigma = sigma_final, 
    realized_snr = realized_snr
  ))
}

#=============================================================
# 3. Performance Evaluation Framework
#=============================================================

#' Comprehensive performance evaluation across methods
#' @param groups Group membership vector
#' @param repeats Number of simulation repetitions
#' @param output_filename File for detailed results
#' @param data_generator Function to generate data
#' @param data_generator_args Arguments for data generator
#' @param output_dir Output directory
#' @return List of all results across PFER values
calculate_metrics_repeated_multi_pfer <- function(groups, repeats = 1, output_filename,
                                                  data_generator = generate_group_data,
                                                  data_generator_args = list(n = 60, snr = 1),
                                                  output_dir = ".") {
  set.seed(123)
  
  # Initialize output directory and file
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  detailed_output_path <- file.path(output_dir, output_filename)
  
  result_file <- file(detailed_output_path, "w")
  cat("Detailed Analysis Results\n", file = result_file)
  cat("=======================\n\n", file = result_file)
  
  # Get data dimensions
  initial_data <- do.call(data_generator, data_generator_args)
  total_vars <- ncol(initial_data$X)
  var_names <- if (!is.null(colnames(initial_data$X))) {
    colnames(initial_data$X)
  } else {
    paste0("V", 1:total_vars)
  }
  
  # Validate group vector
  if (length(groups) != total_vars) {
    stop(sprintf("Length of groups (%d) must match number of variables (%d)", 
                 length(groups), total_vars))
  }
  
  # Log basic information
  cat("Data Information:\n", file = result_file)
  cat(sprintf("Total variables: %d\n", total_vars), file = result_file)
  cat(sprintf("Number of groups: %d\n", length(unique(groups))), file = result_file)
  cat(sprintf("Sample size: %d\n\n", nrow(initial_data$X)), file = result_file)
  
  # Determine PFER values
  max_pfer <- find_max_valid_pfer(initial_data$X, start_pfer = 1)
  pfer_values <- c(1, max_pfer/4, (max_pfer/4)*2, (max_pfer/4)*3, max_pfer)
  rm(initial_data)
  gc()
  
  # Initialize result storage
  all_results <- list()
  for(pfer in pfer_values) {
    pfer_key <- as.character(pfer)
    all_results[[pfer_key]] <- list(
      adaptive_results = list(),
      lasso_results = list(),
      elasticnet_results = list(),
      grplasso_results = list(),
      elastic_results = list(),
      grpreg_results = list()
    )
  }
  
  # Initialize selection matrices for stability calculation
  selection_matrices <- list()
  pfer_keys <- as.character(pfer_values)
  for(pfer_key in pfer_keys) {
    selection_matrices[[pfer_key]] <- list(
      adaptive = matrix(0, nrow = repeats, ncol = total_vars),
      lasso = matrix(0, nrow = repeats, ncol = total_vars),
      elasticnet = matrix(0, nrow = repeats, ncol = total_vars),
      grplasso = matrix(0, nrow = repeats, ncol = total_vars),
      elastic = matrix(0, nrow = repeats, ncol = total_vars),
      grpreg = matrix(0, nrow = repeats, ncol = total_vars)
    )
  }
  
  # Main simulation loop
  for(rep in 1:repeats) {
    cat(sprintf("\nStarting repetition %d of %d\n", rep, repeats))
    cat(sprintf("\nRepetition %d\n", rep), file = result_file)
    cat("----------------------\n", file = result_file)
    
    # Generate data
    set.seed(123 + rep)
    generated_data <- do.call(data_generator, data_generator_args)
    x <- generated_data$X
    y <- generated_data$y
    colnames(x) <- var_names
    
    # Train-test split
    n <- nrow(x)
    train_idx <- sample(1:n, size = floor(0.7 * n))
    test_idx  <- setdiff(1:n, train_idx)
    
    x_train <- x[train_idx, ]
    y_train <- y[train_idx]
    x_test  <- x[test_idx, ]
    y_test  <- y[test_idx]
    
    # Create Laplacian matrix
    L <- create_laplacian_matrix(x_train)
    
    sampling.type <- "SS"
    
    # 1. GLaSS-SS (Adaptive Stability Selection)
    set.seed(123)
    stabsel_results <- my_stabsel_parallel(
      x = x_train, y = y_train,
      fitfun = my_adaptive_fitfun,
      args.fitfun = list(
        groups = groups,
        L_sparse = L,
        alpha_seq = seq(0, 1, 0.25),
        nlambda = 10
      ),
      cutoff = 0.6,
      pfer_values = pfer_values,
      sampling.type = "SS",
      B = 50,
      mc.cores = 20,
      verbose = FALSE
    )
    
    # 2. Lasso with Stability Selection
    set.seed(123)
    lasso_results <- lapply(pfer_values, function(pfer) {
      stabsel(
        x = x_train,
        y = y_train,
        fitfun = glmnet.lasso,
        cutoff = 0.6,
        PFER = pfer,
        sampling.type = "SS",
        assumption = "unimodal",
        mc.preschedule = TRUE,
        B = ifelse(sampling.type == "MB", 100, 50)
      )
    })
    names(lasso_results) <- as.character(pfer_values)
    
    # 3. Elastic Net with Stability Selection
    set.seed(123)
    elasticnet_results <- lapply(pfer_values, function(pfer) {
      stabsel(
        x = x_train,
        y = y_train,
        fitfun = glmnet.elasticnet,
        cutoff = 0.6,
        PFER = pfer,
        sampling.type = "SS",
        assumption = "unimodal",
        mc.preschedule = TRUE,
        B = ifelse(sampling.type == "MB", 100, 50)
      )
    })
    names(elasticnet_results) <- as.character(pfer_values)
    
    # 4. Group Lasso with Stability Selection
    set.seed(123)
    grplasso_results <- lapply(pfer_values, function(pfer) {
      stabsel(
        x = x_train,
        y = y_train,
        fitfun = grpreg.grouplasso,
        args.fitfun = list(group = groups),
        cutoff = 0.6,
        PFER = pfer,
        sampling.type = "SS",
        assumption = "unimodal",
        mc.preschedule = TRUE,
        B = ifelse(sampling.type == "MB", 100, 50)
      )
    })
    names(grplasso_results) <- as.character(pfer_values)
    
    # 5. Elastic Net with Cross-Validation
    alpha_candidates <- seq(0, 1, by = 0.2)
    best_alpha <- 0
    best_auc <- -Inf
    best_cv_fit <- NULL
    
    for(a in alpha_candidates) {
      set.seed(123)
      temp_cv <- cv.glmnet(x_train, y_train, family = "binomial", alpha = a, type.measure = "auc")
      current_max_auc <- max(temp_cv$cvm)
      if(current_max_auc > best_auc) {
        best_auc <- current_max_auc
        best_alpha <- a
        best_cv_fit <- temp_cv
      }
    }
    
    coefs_elastic <- coef(best_cv_fit, s = "lambda.min")
    nonzero_indices <- which(as.numeric(coefs_elastic) != 0)
    selected_elastic <- if(length(nonzero_indices) > 1) {
      nonzero_indices[-1] - 1
    } else {
      integer(0)
    }
    
    # 6. Group Lasso with Cross-Validation
    require(grpreg)
    set.seed(123)
    grpreg_cv <- cv.grpreg(x_train, y_train, group = groups, family = "binomial", penalty = "grLasso")
    
    coefs_grpreg <- coef(grpreg_cv)
    nonzero_indices_grpreg <- which(as.numeric(coefs_grpreg) != 0)
    selected_grpreg <- if(length(nonzero_indices_grpreg) > 1) {
      nonzero_indices_grpreg[-1] - 1
    } else {
      integer(0)
    }
    
    # Process results for each PFER value
    for(pfer in pfer_values) {
      pfer_key <- as.character(pfer)
      current_adaptive <- stabsel_results[[pfer_key]]
      current_lasso <- lasso_results[[pfer_key]]
      current_elasticnet <- elasticnet_results[[pfer_key]]
      current_grplasso <- grplasso_results[[pfer_key]]
      
      # Identify true positives
      true_pos <- which(generated_data$beta_true != 0)
      
      # Calculate false positives
      false_pos_adaptive <- setdiff(current_adaptive$selected, true_pos)
      false_pos_lasso  <- setdiff(current_lasso$selected, true_pos)
      false_pos_elasticnet  <- setdiff(current_elasticnet$selected, true_pos)
      false_pos_grplasso  <- setdiff(current_grplasso$selected, true_pos)
      false_pos_elastic  <- setdiff(selected_elastic, true_pos)
      false_pos_grpreg   <- setdiff(selected_grpreg, true_pos)
      
      # Function to calculate performance metrics
      calculate_metrics <- function(selected_vars, true_vars) {
        if(length(true_vars) > 0) {
          tpr <- length(intersect(selected_vars, true_vars)) / length(true_vars)
        } else {
          tpr <- NA
        }
        
        if(length(selected_vars) > 0) {
          ppv <- length(intersect(selected_vars, true_vars)) / length(selected_vars)
        } else {
          ppv <- NA
        }
        
        if(!is.na(tpr) && !is.na(ppv) && (tpr + ppv) > 0) {
          f1 <- 2 * tpr * ppv / (tpr + ppv)
        } else {
          f1 <- NA
        }
        
        return(list(tpr = tpr, ppv = ppv, f1 = f1))
      }
      
      # Calculate metrics for each method
      metrics_adaptive <- calculate_metrics(current_adaptive$selected, true_pos)
      metrics_lasso <- calculate_metrics(current_lasso$selected, true_pos)
      metrics_elasticnet <- calculate_metrics(current_elasticnet$selected, true_pos)
      metrics_grplasso <- calculate_metrics(current_grplasso$selected, true_pos)
      metrics_elastic <- calculate_metrics(selected_elastic, true_pos)
      metrics_grpreg <- calculate_metrics(selected_grpreg, true_pos)
      
      # Function to calculate AUC using Random Forest
      calculate_auc <- function(selected_vars, x_train, y_train, x_test, y_test) {
        if(length(selected_vars) > 0) {
          train_data <- data.frame(x_train[, selected_vars, drop = FALSE], y = factor(y_train))
          test_data <- data.frame(x_test[, selected_vars, drop = FALSE], y = factor(y_test))
          formula_str <- as.formula(paste("y ~", paste(colnames(train_data)[-ncol(train_data)], collapse = " + ")))
          
          set.seed(5223)
          model <- randomForest(formula_str, data = train_data, ntree = 500)
          pred_prob <- predict(model, newdata = test_data, type = "prob")[,2]
          auc_val <- auc(pROC::roc(test_data$y, pred_prob, direction = "<"))
          return(auc_val)
        } else {
          return(0.5)
        }
      }
      
      # Calculate AUC for each method
      auc_adaptive <- calculate_auc(current_adaptive$selected, x_train, y_train, x_test, y_test)
      auc_lasso <- calculate_auc(current_lasso$selected, x_train, y_train, x_test, y_test)
      auc_elasticnet <- calculate_auc(current_elasticnet$selected, x_train, y_train, x_test, y_test)
      auc_grplasso <- calculate_auc(current_grplasso$selected, x_train, y_train, x_test, y_test)
      auc_elastic <- calculate_auc(selected_elastic, x_train, y_train, x_test, y_test)
      auc_grpreg <- calculate_auc(selected_grpreg, x_train, y_train, x_test, y_test)
      
      # Store results (code continues with storing results and logging)

      all_results[[pfer_key]]$adaptive_results[[rep]] <- list(
        selected = current_adaptive$selected,
        fp = false_pos_adaptive,
        fp_count = length(false_pos_adaptive),
        auc = auc_adaptive,
        tpr = metrics_adaptive$tpr,
        ppv = metrics_adaptive$ppv,
        f1 = metrics_adaptive$f1,
        selection_record = current_adaptive$selection_record,
        params = current_adaptive$params
      )
      
      all_results[[pfer_key]]$lasso_results[[rep]] <- list(
        selected = current_lasso$selected,
        fp = false_pos_lasso,
        fp_count = length(false_pos_lasso),
        auc = auc_lasso,
        tpr = metrics_lasso$tpr,
        ppv = metrics_lasso$ppv,
        f1 = metrics_lasso$f1,
        selection_record = current_lasso$selected
      )
      
      all_results[[pfer_key]]$elasticnet_results[[rep]] <- list(
        selected = current_elasticnet$selected,
        fp = false_pos_elasticnet,
        fp_count = length(false_pos_elasticnet),
        auc = auc_elasticnet,
        tpr = metrics_elasticnet$tpr,
        ppv = metrics_elasticnet$ppv,
        f1 = metrics_elasticnet$f1,
        selection_record = current_elasticnet$selected
      )
      
      all_results[[pfer_key]]$grplasso_results[[rep]] <- list(
        selected = current_grplasso$selected,
        fp = false_pos_grplasso,
        fp_count = length(false_pos_grplasso),
        auc = auc_grplasso,
        tpr = metrics_grplasso$tpr,
        ppv = metrics_grplasso$ppv,
        f1 = metrics_grplasso$f1,
        selection_record = current_grplasso$selected
      )
      
      all_results[[pfer_key]]$elastic_results[[rep]] <- list(
        selected = selected_elastic,
        fp = false_pos_elastic,
        fp_count = length(false_pos_elastic),
        auc = auc_elastic,
        tpr = metrics_elastic$tpr,
        ppv = metrics_elastic$ppv,
        f1 = metrics_elastic$f1,
        best_alpha = best_alpha
      )
      
      all_results[[pfer_key]]$grpreg_results[[rep]] <- list(
        selected = selected_grpreg,
        fp = false_pos_grpreg,
        fp_count = length(false_pos_grpreg),
        auc = auc_grpreg,
        tpr = metrics_grpreg$tpr,
        ppv = metrics_grpreg$ppv,
        f1 = metrics_grpreg$f1,
        lambda_min = grpreg_cv$lambda.min
      )
      
      # log print
      cat(sprintf("\nPFER = %.2f Results:\n", pfer), file = result_file)
      cat(sprintf("Selected Lambda (mean ± sd): %.4f ± %.4f, Selected Alpha (mean ± sd): %.4f ± %.4f\n", 
                  current_adaptive$params_avg$lambda_mean, 
                  current_adaptive$params_avg$lambda_sd, 
                  current_adaptive$params_avg$alpha_mean,  
                  current_adaptive$params_avg$alpha_sd), 
          file = result_file)
      
      cat("AUC:\n", file = result_file)
      cat(sprintf("  GLaSS:     %.4f\n", auc_adaptive), file = result_file)
      cat(sprintf("  Lasso:     %.4f\n", auc_lasso), file = result_file)
      cat(sprintf("  ElasticNet:%.4f\n", auc_elasticnet), file = result_file)
      cat(sprintf("  GrpLasso:  %.4f\n", auc_grplasso), file = result_file)
      cat(sprintf("  EN(CV):    %.4f\n", auc_elastic), file = result_file)
      cat(sprintf("  GrpReg(CV):%.4f\n", auc_grpreg), file = result_file)
      
      cat("TPR, PPV and F1:\n", file = result_file)
      cat(sprintf("  GLaSS:     TPR = %.4f, PPV = %.4f, F1 = %.4f\n", 
                  metrics_adaptive$tpr, metrics_adaptive$ppv, metrics_adaptive$f1), file = result_file)
      cat(sprintf("  Lasso:     TPR = %.4f, PPV = %.4f, F1 = %.4f\n", 
                  metrics_lasso$tpr, metrics_lasso$ppv, metrics_lasso$f1), file = result_file)
      cat(sprintf("  ElasticNet:TPR = %.4f, PPV = %.4f, F1 = %.4f\n", 
                  metrics_elasticnet$tpr, metrics_elasticnet$ppv, metrics_elasticnet$f1), file = result_file)
      cat(sprintf("  GrpLasso:  TPR = %.4f, PPV = %.4f, F1 = %.4f\n", 
                  metrics_grplasso$tpr, metrics_grplasso$ppv, metrics_grplasso$f1), file = result_file)
      cat(sprintf("  EN(CV):    TPR = %.4f, PPV = %.4f, F1 = %.4f\n", 
                  metrics_elastic$tpr, metrics_elastic$ppv, metrics_elastic$f1), file = result_file)
      cat(sprintf("  GrpReg(CV):TPR = %.4f, PPV = %.4f, F1 = %.4f\n", 
                  metrics_grpreg$tpr, metrics_grpreg$ppv, metrics_grpreg$f1), file = result_file)
      
      # 선택 기록(Selection Matrix) 업데이트
      selection_matrices[[pfer_key]]$adaptive[rep, current_adaptive$selected] <- 1
      selection_matrices[[pfer_key]]$lasso[rep, current_lasso$selected] <- 1
      selection_matrices[[pfer_key]]$elasticnet[rep, current_elasticnet$selected] <- 1
      selection_matrices[[pfer_key]]$grplasso[rep, current_grplasso$selected] <- 1
      selection_matrices[[pfer_key]]$elastic[rep, selected_elastic] <- 1
      selection_matrices[[pfer_key]]$grpreg[rep, selected_grpreg] <- 1
      
      # Calculate cumulative stability (if rep > 1)
      if (rep > 1) {
        current_adaptive_matrix <- selection_matrices[[pfer_key]]$adaptive[1:rep, ]
        current_lasso_matrix <- selection_matrices[[pfer_key]]$lasso[1:rep, ]
        current_elasticnet_matrix <- selection_matrices[[pfer_key]]$elasticnet[1:rep, ]
        current_grplasso_matrix <- selection_matrices[[pfer_key]]$grplasso[1:rep, ]
        current_elastic_matrix <- selection_matrices[[pfer_key]]$elastic[1:rep, ]
        current_grpreg_matrix <- selection_matrices[[pfer_key]]$grpreg[1:rep, ]
        
        current_stability_adaptive <- calculate_fold_stability(current_adaptive_matrix)
        current_stability_lasso <- calculate_fold_stability(current_lasso_matrix)
        current_stability_elasticnet <- calculate_fold_stability(current_elasticnet_matrix)
        current_stability_grplasso <- calculate_fold_stability(current_grplasso_matrix)
        current_stability_elastic <- calculate_fold_stability(current_elastic_matrix)
        current_stability_grpreg <- calculate_fold_stability(current_grpreg_matrix)
        
        cat(sprintf("\nCumulative Stability Metrics up to Repetition %d (PFER = %.2f):\n", 
                    rep, pfer), file = result_file)
        cat("  GLaSS:     Stability = ", sprintf("%.4f (%.4f, %.4f)", 
                                                 current_stability_adaptive$stability,
                                                 current_stability_adaptive$stability_ci_lower,
                                                 current_stability_adaptive$stability_ci_upper), 
            "\n", file = result_file)
        cat("  Lasso:     Stability = ", sprintf("%.4f (%.4f, %.4f)", 
                                                 current_stability_lasso$stability,
                                                 current_stability_lasso$stability_ci_lower,
                                                 current_stability_lasso$stability_ci_upper), 
            "\n", file = result_file)
        cat("  ElasticNet:Stability = ", sprintf("%.4f (%.4f, %.4f)", 
                                                 current_stability_elasticnet$stability,
                                                 current_stability_elasticnet$stability_ci_lower,
                                                 current_stability_elasticnet$stability_ci_upper), 
            "\n", file = result_file)
        cat("  GrpLasso:  Stability = ", sprintf("%.4f (%.4f, %.4f)", 
                                                 current_stability_grplasso$stability,
                                                 current_stability_grplasso$stability_ci_lower,
                                                 current_stability_grplasso$stability_ci_upper), 
            "\n", file = result_file)
        cat("  EN(CV):    Stability = ", sprintf("%.4f (%.4f, %.4f)", 
                                                 current_stability_elastic$stability,
                                                 current_stability_elastic$stability_ci_lower,
                                                 current_stability_elastic$stability_ci_upper), 
            "\n", file = result_file)
        cat("  GrpReg(CV):Stability = ", sprintf("%.4f (%.4f, %.4f)", 
                                                 current_stability_grpreg$stability,
                                                 current_stability_grpreg$stability_ci_lower,
                                                 current_stability_grpreg$stability_ci_upper), 
            "\n\n", file = result_file)
      }
      
      cat("Selected Variables:\n", file = result_file)
      cat(sprintf("  GLaSS (%d vars):     %s\n", 
                  length(current_adaptive$selected),
                  paste(sort(current_adaptive$selected), collapse = ", ")), file = result_file)
      cat(sprintf("  Lasso (%d vars):     %s\n", 
                  length(current_lasso$selected),
                  paste(sort(current_lasso$selected), collapse = ", ")), file = result_file)
      cat(sprintf("  ElasticNet (%d vars):%s\n", 
                  length(current_elasticnet$selected),
                  paste(sort(current_elasticnet$selected), collapse = ", ")), file = result_file)
      cat(sprintf("  GrpLasso (%d vars):  %s\n", 
                  length(current_grplasso$selected),
                  paste(sort(current_grplasso$selected), collapse = ", ")), file = result_file)
      cat(sprintf("  EN(CV) (%d vars):    %s\n", 
                  length(selected_elastic),
                  paste(sort(selected_elastic), collapse = ", ")), file = result_file)
      cat(sprintf("  GrpReg(CV) (%d vars):%s\n", 
                  length(selected_grpreg),
                  paste(sort(selected_grpreg), collapse = ", ")), file = result_file)
      cat("----------------------\n", file = result_file)
      
      gc()  # Clean memory
    }
  }
  
  # Write final summary
  cat("\n\nFINAL SUMMARY\n", file = result_file)
  cat("=============\n\n", file = result_file)
  
  for(pfer in pfer_values) {
    pfer_key <- as.character(pfer)
    cat(sprintf("PFER = %.2f\n", pfer), file = result_file)
    cat("-----------------\n", file = result_file)
    
    # Lambda and Alpha statistics (GLaSS only)
    rep_lambda_means <- numeric(length(all_results[[pfer_key]]$adaptive_results))
    rep_alpha_means <- numeric(length(all_results[[pfer_key]]$adaptive_results))
    
    for(i in seq_along(all_results[[pfer_key]]$adaptive_results)) {
      current_result <- all_results[[pfer_key]]$adaptive_results[[i]]
      if(!is.null(current_result) && !is.null(current_result$params)) {
        params_vec <- unlist(current_result$params)
        lambda_values <- params_vec[names(params_vec) == "lambda"]
        alpha_values <- params_vec[names(params_vec) == "alpha"]
        
        rep_lambda_means[i] <- mean(lambda_values, na.rm = TRUE)
        rep_alpha_means[i] <- mean(alpha_values, na.rm = TRUE)
      } else {
        rep_lambda_means[i] <- NA
        rep_alpha_means[i] <- NA
      }
    }
    
    valid_lambdas <- rep_lambda_means[!is.na(rep_lambda_means)]
    valid_alphas <- rep_alpha_means[!is.na(rep_alpha_means)]
    
    cat("Mean Parameters across all repetitions:\n", file = result_file)
    if(length(valid_lambdas) > 0) {
      cat(sprintf("  Lambda: %.4f (%.4f)\n", 
                  mean(valid_lambdas), 
                  sd(valid_lambdas)), file = result_file)
    } else {
      cat("  Lambda: No valid data\n", file = result_file)
    }
    
    if(length(valid_alphas) > 0) {
      cat(sprintf("  Alpha:  %.4f (%.4f)\n", 
                  mean(valid_alphas), 
                  sd(valid_alphas)), file = result_file)
    } else {
      cat("  Alpha: No valid data\n", file = result_file)
    }
    
    # Performance statistics for each method
    methods <- c("adaptive", "lasso", "elasticnet", "grplasso", "elastic", "grpreg")
    method_names <- c("Adaptive", "Regular", "ElasticNet", "GrpLasso", "Elastic", "GrpReg")
    
    # AUC statistics
    cat("\nMean AUC (SD):\n", file = result_file)
    for(i in seq_along(methods)) {
      method <- methods[i]
      method_name <- method_names[i]
      method_key <- paste0(method, "_results")
      
      if(method_key %in% names(all_results[[pfer_key]])) {
        aucs <- sapply(all_results[[pfer_key]][[method_key]], function(x) x$auc)
        cat(sprintf("  %s: %.4f (%.4f)\n", 
                    method_name, 
                    mean(aucs, na.rm = TRUE), 
                    sd(aucs, na.rm = TRUE)), file = result_file)
      }
    }
    
    # TPR statistics
    cat("\nMean TPR (SD):\n", file = result_file)
    for(i in seq_along(methods)) {
      method <- methods[i]
      method_name <- method_names[i]
      method_key <- paste0(method, "_results")
      
      if(method_key %in% names(all_results[[pfer_key]])) {
        tprs <- sapply(all_results[[pfer_key]][[method_key]], function(x) x$tpr)
        cat(sprintf("  %s: %.4f (%.4f)\n", 
                    method_name, 
                    mean(tprs, na.rm = TRUE), 
                    sd(tprs, na.rm = TRUE)), file = result_file)
      }
    }
    
    # PPV statistics
    cat("\nMean PPV (SD):\n", file = result_file)
    for(i in seq_along(methods)) {
      method <- methods[i]
      method_name <- method_names[i]
      method_key <- paste0(method, "_results")
      
      if(method_key %in% names(all_results[[pfer_key]])) {
        ppvs <- sapply(all_results[[pfer_key]][[method_key]], function(x) x$ppv)
        cat(sprintf("  %s: %.4f (%.4f)\n", 
                    method_name, 
                    mean(ppvs, na.rm = TRUE), 
                    sd(ppvs, na.rm = TRUE)), file = result_file)
      }
    }
    
    # F1 statistics
    cat("\nMean F1 (SD):\n", file = result_file)
    for(i in seq_along(methods)) {
      method <- methods[i]
      method_name <- method_names[i]
      method_key <- paste0(method, "_results")
      
      if(method_key %in% names(all_results[[pfer_key]])) {
        f1s <- sapply(all_results[[pfer_key]][[method_key]], function(x) x$f1)
        cat(sprintf("  %s: %.4f (%.4f)\n", 
                    method_name, 
                    mean(f1s, na.rm = TRUE), 
                    sd(f1s, na.rm = TRUE)), file = result_file)
      }
    }
    
    # Stability analysis - Calculate Nogueira Stability
    cat("\nStability Analysis:\n", file = result_file)
    
    # Generate selection matrix and calculate stability for each method
    n_reps <- length(all_results[[pfer_key]]$adaptive_results)
    
    # Find maximum variable index
    max_var_index <- 0
    for(i in seq_along(methods)) {
      method_key <- paste0(methods[i], "_results")
      if(method_key %in% names(all_results[[pfer_key]])) {
        for(rep in 1:n_reps) {
          selected_vars <- all_results[[pfer_key]][[method_key]][[rep]]$selected
          if(length(selected_vars) > 0) {
            max_var_index <- max(max_var_index, max(selected_vars))
          }
        }
      }
    }
    
    # Calculate stability for each method
    for(i in seq_along(methods)) {
      method <- methods[i]
      method_name <- method_names[i]
      method_key <- paste0(method, "_results")
      
      if(method_key %in% names(all_results[[pfer_key]])) {
        # Generate selection matrix
        selection_matrix <- matrix(0, nrow = n_reps, ncol = max_var_index)
        
        for(rep in 1:n_reps) {
          selected_vars <- all_results[[pfer_key]][[method_key]][[rep]]$selected
          if(length(selected_vars) > 0) {
            selection_matrix[rep, selected_vars] <- 1
          }
        }
        
        # Calculate Nogueira stability
        stability_result <- calculate_stability(selection_matrix)
        cat(sprintf("%s:\n", method_name), file = result_file)
        cat(sprintf("  Stability: %.4f (%.4f, %.4f)\n", 
                    stability_result$stability,
                    stability_result$ci_lower,
                    stability_result$ci_upper), file = result_file)
      }
    }
    
    # Jaccard stability analysis
    cat("\nJaccard Stability (SD):\n", file = result_file)
    for(i in seq_along(methods)) {
      method <- methods[i]
      method_name <- method_names[i]
      method_key <- paste0(method, "_results")
      
      if(method_key %in% names(all_results[[pfer_key]])) {
        # Generate selection matrix (recalculated although already computed above)
        selection_matrix <- matrix(0, nrow = n_reps, ncol = max_var_index)
        
        for(rep in 1:n_reps) {
          selected_vars <- all_results[[pfer_key]][[method_key]][[rep]]$selected
          if(length(selected_vars) > 0) {
            selection_matrix[rep, selected_vars] <- 1
          }
        }
        
        # Calculate Jaccard stability
        jaccard_result <- calculate_jaccard_stability(selection_matrix)
        cat(sprintf("%s:\n", method_name), file = result_file)
        cat(sprintf("  Jaccard: %.4f (%.4f)\n", 
                    jaccard_result$stability,
                    jaccard_result$sd), file = result_file)
      }
    }
    
    # Statistics for number of selected variables
    cat("\nNumber of Selected Variables:\n", file = result_file)
    for(i in seq_along(methods)) {
      method <- methods[i]
      method_name <- method_names[i]
      method_key <- paste0(method, "_results")
      
      if(method_key %in% names(all_results[[pfer_key]])) {
        n_vars <- sapply(all_results[[pfer_key]][[method_key]], 
                         function(x) length(x$selected))
        cat(sprintf("%s:\n", method_name), file = result_file)
        cat(sprintf("  Mean: %.2f (SD: %.2f)\n", 
                    mean(n_vars, na.rm = TRUE), 
                    sd(n_vars, na.rm = TRUE)), file = result_file)
      }
    }
    
    # Range of selected variables
    cat("\nRange of Selected Variables:\n", file = result_file)
    for(i in seq_along(methods)) {
      method <- methods[i]
      method_name <- method_names[i]
      method_key <- paste0(method, "_results")
      
      if(method_key %in% names(all_results[[pfer_key]])) {
        n_vars <- sapply(all_results[[pfer_key]][[method_key]], 
                         function(x) length(x$selected))
        cat(sprintf("%s: Min = %d, Max = %d\n", 
                    method_name, 
                    min(n_vars, na.rm = TRUE), 
                    max(n_vars, na.rm = TRUE)), file = result_file)
      }
    }
    
    cat("\n", file = result_file)
  }
  
  close(result_file)
  return(all_results)
}

#=============================================================
# 4. Visualization Functions
#=============================================================

#' Generate comprehensive plots from simulation results
#' @param sim_result Simulation results object from calculate_metrics_repeated_multi_pfer
#' @param scenario_key Identifier for scenario (used in filename)
#' @param output_dir Directory for saving plots
#' @return List containing plot data and ggplot objects
make_plots_from_sim_result <- function(sim_result, scenario_key = NULL, output_dir = ".") {
  library(ggplot2)
  library(gridExtra)
  
  # ============================================================
  # 1. Prepare data for plots
  # ============================================================
  plot_data <- data.frame()    # FP and TPR data
  auc_data <- data.frame()     # AUC values
  n_vars_data <- data.frame()  # Number of selected variables
  
  for (pfer_key in names(sim_result)) {
    pfer_value <- as.numeric(pfer_key)
    n_rep <- length(sim_result[[pfer_key]]$adaptive_results)
    
    for (rep_idx in seq_len(n_rep)) {
      # Extract results for GLaSS-SS
      fp_count <- sim_result[[pfer_key]]$adaptive_results[[rep_idx]]$fp_count
      tpr_val  <- sim_result[[pfer_key]]$adaptive_results[[rep_idx]]$tpr
      auc_val  <- sim_result[[pfer_key]]$adaptive_results[[rep_idx]]$auc
      n_vars   <- length(sim_result[[pfer_key]]$adaptive_results[[rep_idx]]$selected)
      
      plot_data <- rbind(plot_data, data.frame(
        method = "GLaSS stabsel",
        pfer   = as.character(pfer_value),
        FP     = fp_count,
        TPR    = tpr_val,
        stringsAsFactors = FALSE
      ))
      
      auc_data <- rbind(auc_data, data.frame(
        method = "GLaSS stabsel",
        pfer   = as.character(pfer_value),
        AUC    = auc_val,
        stringsAsFactors = FALSE
      ))
      
      n_vars_data <- rbind(n_vars_data, data.frame(
        method = "GLaSS stabsel",
        pfer   = as.character(pfer_value),
        n_vars = n_vars,
        stringsAsFactors = FALSE
      ))
      
      # Extract results for Lasso with stability selection
      fp_count <- sim_result[[pfer_key]]$lasso_results[[rep_idx]]$fp_count
      tpr_val  <- sim_result[[pfer_key]]$lasso_results[[rep_idx]]$tpr
      auc_val  <- sim_result[[pfer_key]]$lasso_results[[rep_idx]]$auc
      n_vars   <- length(sim_result[[pfer_key]]$lasso_results[[rep_idx]]$selected)
      
      plot_data <- rbind(plot_data, data.frame(
        method = "Lasso stabsel",
        pfer   = as.character(pfer_value),
        FP     = fp_count,
        TPR    = tpr_val,
        stringsAsFactors = FALSE
      ))
      
      auc_data <- rbind(auc_data, data.frame(
        method = "Lasso stabsel",
        pfer   = as.character(pfer_value),
        AUC    = auc_val,
        stringsAsFactors = FALSE
      ))
      
      n_vars_data <- rbind(n_vars_data, data.frame(
        method = "Lasso stabsel",
        pfer   = as.character(pfer_value),
        n_vars = n_vars,
        stringsAsFactors = FALSE
      ))
      
      # Extract results for ElasticNet with stability selection
      fp_count <- sim_result[[pfer_key]]$elasticnet_results[[rep_idx]]$fp_count
      tpr_val  <- sim_result[[pfer_key]]$elasticnet_results[[rep_idx]]$tpr
      auc_val  <- sim_result[[pfer_key]]$elasticnet_results[[rep_idx]]$auc
      n_vars   <- length(sim_result[[pfer_key]]$elasticnet_results[[rep_idx]]$selected)
      
      plot_data <- rbind(plot_data, data.frame(
        method = "EN stabsel",
        pfer   = as.character(pfer_value),
        FP     = fp_count,
        TPR    = tpr_val,
        stringsAsFactors = FALSE
      ))
      
      auc_data <- rbind(auc_data, data.frame(
        method = "EN stabsel",
        pfer   = as.character(pfer_value),
        AUC    = auc_val,
        stringsAsFactors = FALSE
      ))
      
      n_vars_data <- rbind(n_vars_data, data.frame(
        method = "EN stabsel",
        pfer   = as.character(pfer_value),
        n_vars = n_vars,
        stringsAsFactors = FALSE
      ))
      
      # Extract results for Group Lasso with stability selection
      fp_count <- sim_result[[pfer_key]]$grplasso_results[[rep_idx]]$fp_count
      tpr_val  <- sim_result[[pfer_key]]$grplasso_results[[rep_idx]]$tpr
      auc_val  <- sim_result[[pfer_key]]$grplasso_results[[rep_idx]]$auc
      n_vars   <- length(sim_result[[pfer_key]]$grplasso_results[[rep_idx]]$selected)
      
      plot_data <- rbind(plot_data, data.frame(
        method = "grLasso stabsel",
        pfer   = as.character(pfer_value),
        FP     = fp_count,
        TPR    = tpr_val,
        stringsAsFactors = FALSE
      ))
      
      auc_data <- rbind(auc_data, data.frame(
        method = "grLasso stabsel",
        pfer   = as.character(pfer_value),
        AUC    = auc_val,
        stringsAsFactors = FALSE
      ))
      
      n_vars_data <- rbind(n_vars_data, data.frame(
        method = "grLasso stabsel",
        pfer   = as.character(pfer_value),
        n_vars = n_vars,
        stringsAsFactors = FALSE
      ))
      
      # Extract results for Elastic Net (CV-based)
      # Note: PFER doesn't apply to CV methods, so we use a fixed label
      fp_count <- sim_result[[pfer_key]]$elastic_results[[rep_idx]]$fp_count
      tpr_val  <- sim_result[[pfer_key]]$elastic_results[[rep_idx]]$tpr
      auc_val  <- sim_result[[pfer_key]]$elastic_results[[rep_idx]]$auc
      n_vars   <- length(sim_result[[pfer_key]]$elastic_results[[rep_idx]]$selected)
      
      plot_data <- rbind(plot_data, data.frame(
        method = "EN(CV)",
        pfer   = "EN(CV)",
        FP     = fp_count,
        TPR    = tpr_val,
        stringsAsFactors = FALSE
      ))
      
      auc_data <- rbind(auc_data, data.frame(
        method = "EN(CV)",
        pfer   = "EN(CV)",
        AUC    = auc_val,
        stringsAsFactors = FALSE
      ))
      
      n_vars_data <- rbind(n_vars_data, data.frame(
        method = "EN(CV)",
        pfer   = "EN(CV)",
        n_vars = n_vars,
        stringsAsFactors = FALSE
      ))
      
      # Extract results for Group Lasso (CV-based)
      fp_count <- sim_result[[pfer_key]]$grpreg_results[[rep_idx]]$fp_count
      tpr_val  <- sim_result[[pfer_key]]$grpreg_results[[rep_idx]]$tpr
      auc_val  <- sim_result[[pfer_key]]$grpreg_results[[rep_idx]]$auc
      n_vars   <- length(sim_result[[pfer_key]]$grpreg_results[[rep_idx]]$selected)
      
      plot_data <- rbind(plot_data, data.frame(
        method = "grLasso(CV)",
        pfer   = "grLasso(CV)",
        FP     = fp_count,
        TPR    = tpr_val,
        stringsAsFactors = FALSE
      ))
      
      auc_data <- rbind(auc_data, data.frame(
        method = "grLasso(CV)",
        pfer   = "grLasso(CV)",
        AUC    = auc_val,
        stringsAsFactors = FALSE
      ))
      
      n_vars_data <- rbind(n_vars_data, data.frame(
        method = "grLasso(CV)",
        pfer   = "grLasso(CV)",
        n_vars = n_vars,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # ============================================================
  # 2. Set up factor levels for proper ordering
  # ============================================================
  
  # Identify stability selection methods
  stability_methods <- c("GLaSS stabsel", "Lasso stabsel", "EN stabsel", "grLasso stabsel")
  non_special_methods <- plot_data$method %in% stability_methods
  non_special_levels <- sort(unique(as.numeric(plot_data$pfer[non_special_methods])))
  
  # Map PFER values to meaningful labels (min, Q1, median, Q3, max)
  pfer_names <- c()
  for (i in seq_along(non_special_levels)) {
    if (non_special_levels[i] == 1) {
      pfer_names[i] <- "min"
    } else if (non_special_levels[i] == 9) {
      pfer_names[i] <- "Q1"
    } else if (non_special_levels[i] == 18) {
      pfer_names[i] <- "median"
    } else if (non_special_levels[i] == 27) {
      pfer_names[i] <- "Q3"
    } else if (non_special_levels[i] == 36) {
      pfer_names[i] <- "max"
    } else {
      pfer_names[i] <- as.character(non_special_levels[i])
    }
  }
  
  # Create mapping between original values and new names
  pfer_mapping <- setNames(pfer_names, as.character(non_special_levels))
  
  # Apply PFER value mapping to data frames
  for (i in seq_along(non_special_levels)) {
    old_value <- as.character(non_special_levels[i])
    new_value <- pfer_names[i]
    
    plot_data$pfer[plot_data$pfer == old_value] <- new_value
    auc_data$pfer[auc_data$pfer == old_value] <- new_value
    n_vars_data$pfer[n_vars_data$pfer == old_value] <- new_value
  }
  
  final_levels <- c(pfer_names, "EN(CV)", "grLasso(CV)")
  
  # Set factor levels to control legend order
  method_levels <- c("GLaSS stabsel", "Lasso stabsel", "EN stabsel", "grLasso stabsel", "EN(CV)", "grLasso(CV)")
  plot_data$method <- factor(plot_data$method, levels = method_levels)
  auc_data$method <- factor(auc_data$method, levels = method_levels)
  n_vars_data$method <- factor(n_vars_data$method, levels = method_levels)
  
  plot_data$pfer <- factor(plot_data$pfer, levels = final_levels)
  auc_data$pfer <- factor(auc_data$pfer, levels = final_levels)
  n_vars_data$pfer <- factor(n_vars_data$pfer, levels = final_levels)
  
  # Create common theme for all plots
  center_title_theme <- theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )
  
  # ============================================================
  # 3. Create False Positive plot
  # ============================================================
  p_fp <- ggplot(plot_data, aes(x = pfer, y = FP, fill = method)) +
    geom_boxplot(outlier.shape = 21, alpha = 0.7) +
    labs(
      title = "False Positive Distribution",
      x = "Methods",
      y = "False Positives",
      fill = "Method"
    ) +
    theme_bw(base_size = 18) +
    center_title_theme +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))  # 2-row horizontal legend
  
  # Add red horizontal lines to show PFER limits
  unique_limits <- data.frame(
    pfer_old = as.character(non_special_levels),
    pfer = pfer_names,
    x = as.numeric(factor(pfer_names, levels = pfer_names)),
    y = non_special_levels
  )
  
  p_fp <- p_fp + 
    geom_segment(data = unique_limits, 
                 mapping = aes(x = x - 0.4, xend = x + 0.4, y = y, yend = y),
                 inherit.aes = FALSE, color = "red", size = 1) +
    scale_y_continuous(limits = c(0, 60),
                       breaks = seq(0, 60, by = 20))
  
  # ============================================================
  # 4. Create True Positive Rate plot
  # ============================================================
  p_tpr <- ggplot(plot_data, aes(x = pfer, y = TPR, fill = method)) +
    geom_boxplot(outlier.shape = 21, alpha = 0.7) +
    labs(
      title = "True Positive Rate Distribution",
      x = "Methods",
      y = "True Positive Rate",
      fill = "Method"
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_bw(base_size = 18) +
    center_title_theme +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
  
  # ============================================================
  # 5. Create AUC plot
  # ============================================================
  p_auc <- ggplot(auc_data, aes(x = pfer, y = AUC, fill = method)) +
    geom_boxplot(outlier.shape = 21, alpha = 0.7) +
    labs(
      title = "AUC Distribution",
      x = "Methods",
      y = "AUC",
      fill = "Method"
    ) +
    scale_y_continuous(limits = c(0.5, 1)) +
    theme_bw(base_size = 18) +
    center_title_theme +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
  
  # ============================================================
  # 6. Create Number of Selected Variables plot
  # ============================================================
  p_nvars <- ggplot(n_vars_data, aes(x = pfer, y = n_vars, fill = method)) +
    geom_boxplot(outlier.shape = 21, alpha = 0.7) +
    labs(
      title = "Selected Variables Distribution",
      x = "Methods",
      y = "Number of Selected Variables",
      fill = "Method"
    ) +
    scale_y_continuous(limits = c(0, 120),
                       breaks = seq(0, 120, by = 20)) +
    theme_bw(base_size = 18) +
    center_title_theme +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
  
  # ============================================================
  # 7. Prepare data for stability plots
  # ============================================================
  stability_data <- data.frame()
  jaccard_data <- data.frame()
  
  for (pfer_key in names(sim_result)) {
    pfer_value <- as.numeric(pfer_key)
    
    # Apply PFER value mapping
    pfer_label <- as.character(pfer_value)
    if (pfer_label %in% names(pfer_mapping)) {
      pfer_label <- pfer_mapping[pfer_label]
    }
    
    # Calculate stability metrics from selection matrices
    n_reps <- length(sim_result[[pfer_key]]$adaptive_results)
    
    # Determine maximum variable index
    max_var_index <- 0
    for (rep in 1:n_reps) {
      max_var_index <- max(max_var_index, 
                           max(c(0, sim_result[[pfer_key]]$adaptive_results[[rep]]$selected)),
                           max(c(0, sim_result[[pfer_key]]$lasso_results[[rep]]$selected)),
                           max(c(0, sim_result[[pfer_key]]$elasticnet_results[[rep]]$selected)),
                           max(c(0, sim_result[[pfer_key]]$grplasso_results[[rep]]$selected)),
                           max(c(0, sim_result[[pfer_key]]$elastic_results[[rep]]$selected)),
                           max(c(0, sim_result[[pfer_key]]$grpreg_results[[rep]]$selected)))
    }
    
    # Initialize selection matrices
    adaptive_selection_matrix <- matrix(0, nrow = n_reps, ncol = max_var_index)
    lasso_selection_matrix <- matrix(0, nrow = n_reps, ncol = max_var_index)
    elasticnet_selection_matrix <- matrix(0, nrow = n_reps, ncol = max_var_index)
    grplasso_selection_matrix <- matrix(0, nrow = n_reps, ncol = max_var_index)
    elastic_selection_matrix <- matrix(0, nrow = n_reps, ncol = max_var_index)
    grpreg_selection_matrix <- matrix(0, nrow = n_reps, ncol = max_var_index)
    
    # Fill selection matrices
    for (rep in 1:n_reps) {
      if (length(sim_result[[pfer_key]]$adaptive_results[[rep]]$selected) > 0) {
        adaptive_selection_matrix[rep, sim_result[[pfer_key]]$adaptive_results[[rep]]$selected] <- 1
      }
      if (length(sim_result[[pfer_key]]$lasso_results[[rep]]$selected) > 0) {
        lasso_selection_matrix[rep, sim_result[[pfer_key]]$lasso_results[[rep]]$selected] <- 1
      }
      if (length(sim_result[[pfer_key]]$elasticnet_results[[rep]]$selected) > 0) {
        elasticnet_selection_matrix[rep, sim_result[[pfer_key]]$elasticnet_results[[rep]]$selected] <- 1
      }
      if (length(sim_result[[pfer_key]]$grplasso_results[[rep]]$selected) > 0) {
        grplasso_selection_matrix[rep, sim_result[[pfer_key]]$grplasso_results[[rep]]$selected] <- 1
      }
      if (length(sim_result[[pfer_key]]$elastic_results[[rep]]$selected) > 0) {
        elastic_selection_matrix[rep, sim_result[[pfer_key]]$elastic_results[[rep]]$selected] <- 1
      }
      if (length(sim_result[[pfer_key]]$grpreg_results[[rep]]$selected) > 0) {
        grpreg_selection_matrix[rep, sim_result[[pfer_key]]$grpreg_results[[rep]]$selected] <- 1
      }
    }
    
    # Calculate stability metrics for each method
    methods_list <- list(
      "GLaSS stabsel" = adaptive_selection_matrix,
      "Lasso stabsel" = lasso_selection_matrix,
      "EN stabsel" = elasticnet_selection_matrix,
      "grLasso stabsel" = grplasso_selection_matrix,
      "EN(CV)" = elastic_selection_matrix,
      "grLasso(CV)" = grpreg_selection_matrix
    )
    
    for (method_name in names(methods_list)) {
      selection_matrix <- methods_list[[method_name]]
      
      # Calculate Nogueira stability
      stability_result <- calculate_stability(selection_matrix)
      
      # Set appropriate PFER label (CV methods don't use PFER)
      current_pfer_label <- if (method_name %in% c("EN(CV)", "grLasso(CV)")) {
        method_name
      } else {
        pfer_label
      }
      
      stability_data <- rbind(stability_data, data.frame(
        method = method_name,
        pfer = current_pfer_label,
        stability = stability_result$stability,
        ci_lower = stability_result$ci_lower,
        ci_upper = stability_result$ci_upper,
        stringsAsFactors = FALSE
      ))
      
      # Calculate Jaccard stability
      jaccard_result <- calculate_jaccard_stability(selection_matrix)
      
      jaccard_data <- rbind(jaccard_data, data.frame(
        method = method_name,
        pfer = current_pfer_label,
        stability = jaccard_result$stability,
        ci_lower = jaccard_result$ci_lower,
        ci_upper = jaccard_result$ci_upper,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Set factor levels for stability data
  stability_data$method <- factor(stability_data$method, levels = method_levels)
  jaccard_data$method <- factor(jaccard_data$method, levels = method_levels)
  stability_data$pfer <- factor(stability_data$pfer, levels = final_levels)
  jaccard_data$pfer <- factor(jaccard_data$pfer, levels = final_levels)
  
  # ============================================================
  # 8. Create Nogueira Stability plot
  # ============================================================
  p_stability <- ggplot(stability_data, aes(x = pfer, y = stability, color = method, group = method)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = -Inf, fill = "gray90", alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                  position = position_dodge(width = 0.5),
                  width = 0.2) +
    labs(
      title = "Nogueira Stability",
      x = "Methods",
      y = "Nogueira Stability",
      color = "Method"
    ) +
    scale_y_continuous(limits = c(-0.3, 1), 
                       breaks = seq(-0.3, 1, by = 0.1),
                       labels = function(x) sprintf("%.1f", x)) +
    theme_bw(base_size = 18) +
    center_title_theme +
    theme(panel.grid.minor = element_blank()) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE))
  
  # ============================================================
  # 9. Create Jaccard Stability plot
  # ============================================================
  p_jaccard <- ggplot(jaccard_data, aes(x = pfer, y = stability, color = method, group = method)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    labs(
      title = "Jaccard Stability",
      x = "Methods",
      y = "Jaccard Stability",
      color = "Method"
    ) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
    theme_bw(base_size = 18) +
    center_title_theme +
    theme(panel.grid.minor = element_blank()) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE))
  
  # ============================================================
  # 10. Save combined plots if scenario_key is provided
  # ============================================================
  if (!is.null(scenario_key)) {
    # Function to extract legend from a plot
    get_legend <- function(myplot) {
      tmp <- ggplot_build(myplot)
      g <- ggplot_gtable(tmp)
      leg <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
      legend <- g$grobs[[leg]]
      return(legend)
    }
    
    # Configure legend for bottom placement with horizontal orientation
    p_fp <- p_fp + theme(legend.position = "bottom", legend.direction = "horizontal")
    common_legend <- get_legend(p_fp)
    
    # Remove individual legends from all plots
    p_fp <- p_fp + theme(legend.position = "none")
    p_tpr <- p_tpr + theme(legend.position = "none")
    p_auc <- p_auc + theme(legend.position = "none")
    p_nvars <- p_nvars + theme(legend.position = "none")
    p_stability <- p_stability + theme(legend.position = "none")
    p_jaccard <- p_jaccard + theme(legend.position = "none")
    
    # Convert plots to grobs to ensure theme application
    g1 <- ggplotGrob(p_fp)
    g2 <- ggplotGrob(p_tpr)
    g3 <- ggplotGrob(p_auc)
    g4 <- ggplotGrob(p_nvars)
    g5 <- ggplotGrob(p_stability)
    g6 <- ggplotGrob(p_jaccard)
    
    # Arrange all plots on a single page with common legend
    all_plots_filename <- file.path(output_dir, paste0("all_plots_", scenario_key, ".pdf"))
    pdf(all_plots_filename, width = 18, height = 20)
    grid.arrange(
      arrangeGrob(g1, g2, g3, g4, g5, g6, ncol = 2),
      common_legend,
      nrow = 2,
      heights = c(10, 1)  # 10:1 ratio for plots:legend
    )
    dev.off()
  }
  
  # ============================================================
  # 11. Return all plot data and objects
  # ============================================================
  return(list(
    plot_data = plot_data,
    auc_data = auc_data,
    n_vars_data = n_vars_data,
    stability_data = stability_data,
    jaccard_data = jaccard_data,
    fp_plot = p_fp,
    tpr_plot = p_tpr,
    auc_plot = p_auc,
    n_vars_plot = p_nvars,
    stability_plot = p_stability,
    jaccard_plot = p_jaccard
  ))
}

#=============================================================
# 5. Main Simulation Execution
#=============================================================

# Set random seed for reproducibility
set.seed(42)

# Define simulation parameters
n_vals   <- c(120, 60)  # Sample sizes
snr_vals <- c(1)        # Signal-to-noise ratios

# Define group structures to test
groups_list <- list(
  rep(1:6, each = 20),   # 6 groups, 20 variables each (well-specified)
  rep(1:3, each = 40),   # 3 groups, 40 variables each (misspecified)
  rep(1:12, each = 10),  # 12 groups, 10 variables each (finer groups)
  rep(1:1, each = 120),  # 1 group, all variables (completely misspecified)
  rep(1:120, each = 1)   # 120 groups, 1 variable each (no grouping)
)
groups_names <- c("6groups_20each", "3group_40each", "12group_10each", 
                  "1groups_120each", "120groups_1each")

# Signal pattern types
half_vals <- c(1, 3)
half_names <- c("block_group", "NULL", "sparse_group")

# Number of repetitions for each scenario
repeats <- 30

# Create output directory with timestamp
output_dir_1 <- format(Sys.time(), "simulation_results_%Y-%m-%d_%H-%M-%S")
if (!dir.exists(output_dir_1)) {
  dir.create(output_dir_1, recursive = TRUE)
}

# Create list of all parameter combinations
scenarios <- expand.grid(
  n_val = n_vals,
  group_idx = seq_along(groups_list),
  snr_val = snr_vals,
  half_vals = half_vals
)
scenarios <- scenarios[order(scenarios$snr_val, scenarios$n_val, decreasing = TRUE), ]
rownames(scenarios) <- NULL

# Add readable names to scenarios
scenarios$current_groups_name <- groups_names[scenarios$group_idx]
scenarios$current_half_names <- half_names[scenarios$half_vals]

# Execute simulations
for(i in 1:nrow(scenarios)) {
  # Extract parameters for current scenario
  n_val <- scenarios$n_val[i]
  snr_val <- scenarios$snr_val[i]
  current_groups <- groups_list[[scenarios$group_idx[i]]]
  current_groups_name <- scenarios$current_groups_name[i]
  half_val <- scenarios$half_vals[i]
  half_val_name <- scenarios$current_half_names[i]
  
  # Generate scenario identifier and output filename
  scenario_key <- sprintf("%s_n%d_snr%.1f_%s", 
                          half_val_name, n_val, snr_val, current_groups_name)
  output_filename <- sprintf("%s_n%d_snr%.1f_%s_B50_30rep.txt", 
                            half_val_name, n_val, snr_val, current_groups_name)
  
  cat(sprintf("Running simulation %d/%d: n = %d, snr = %.1f, groups = %s, pattern = %s\n", 
              i, nrow(scenarios), n_val, snr_val, current_groups_name, half_val_name))
  
  # Run simulation
  sim_result <- calculate_metrics_repeated_multi_pfer(
    groups = current_groups,
    repeats = repeats,
    data_generator_args = list(n = n_val, snr = snr_val, half = half_val),
    output_filename = output_filename,
    output_dir = output_dir_1
  )
  
  # Save simulation result
  sim_result_file <- file.path(output_dir_1, sprintf("sim_result_%s.rds", scenario_key))
  saveRDS(sim_result, sim_result_file)
  cat(sprintf("Saved simulation result to %s\n", sim_result_file))
  
  # Generate and save plots
  plot_results <- make_plots_from_sim_result(
    sim_result = sim_result,
    scenario_key = output_filename,
    output_dir = output_dir_1
  )
}

cat("\n=== Simulation study completed successfully ===\n")
cat(sprintf("Results saved in: %s\n", output_dir_1))
