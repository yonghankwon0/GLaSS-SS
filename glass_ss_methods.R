#=============================================================
# GLaSS-SS: Group-Laplacian Structured Shrinkage with Stability Selection
# Core Methodology Implementation
#=============================================================

#=============================================================
# Required Packages
#=============================================================
library(Matrix)
library(Rcpp)
library(glmnet)
library(grpreg)
library(huge)
library(stabs)
library(parallel)

#=============================================================
# 1. Graph Construction Functions
#=============================================================

#' Create optimal adjacency matrix using GLASSO with STARS selection
#' @param X Input data matrix (n x p)
#' @return List containing adjacency matrix and optimal lambda
create_optimal_adjacency_glasso <- function(X) {
  set.seed(123)
  huge_model <- huge(X, method = "glasso", nlambda = 100)
  cv_model <- huge.select(huge_model)
  
  precision_matrix <- cv_model$refit %>% as.matrix()
  adj_matrix <- abs(precision_matrix)
  diag(adj_matrix) <- 0
  
  return(list(adj_matrix = adj_matrix, optimal_lambda = cv_model$opt.lambda))
}

#' Create Laplacian matrix from correlation structure
#' @param X Input data matrix (n x p)
#' @return Sparse Laplacian matrix
create_laplacian_matrix <- function(X) {
  cor_matrix <- cor(X)
  cor_matrix <- (1 + cor_matrix)/2  # Scale to [0,1]
  diag(cor_matrix) <- 0
  degree_matrix <- diag(rowSums(cor_matrix))
  L <- degree_matrix - cor_matrix
  return(as(L, "dgCMatrix"))
}

#=============================================================
# 2. GLaSS Optimization Algorithm
#=============================================================

#' Fit GLaSS model with adaptive penalty optimization
#' Uses Generalized Forward-Backward splitting algorithm
#' 
#' @param X Design matrix (n x p)
#' @param y Response vector (binary)
#' @param L_sparse Laplacian matrix (sparse)
#' @param groups Group membership vector
#' @param lambda_seq Sequence of lambda values (optional)
#' @param alpha_seq Sequence of alpha values for mixing parameter
#' @param nlambda Number of lambda values to generate
#' @param n_cores Number of cores for parallel processing
#' @param tol Convergence tolerance
#' @param max_iter Maximum number of iterations
#' @param step_size Initial step size (auto-computed if NULL)
#' @param step_size_decay Step size decay factor
#' @param min_step_size Minimum step size
#' @param relaxation Relaxation parameter for proximal updates
#' @param early_stop_count Number of stalls before early stopping
#' @return List with optimal beta, deviance, and fitting details
fit_adaptive_penalty_optim_noCV <- function(
    X, y, L_sparse, groups,
    lambda_seq = NULL,
    alpha_seq = seq(0, 1, length.out = 5),
    nlambda   = 10,
    n_cores   = 1,
    tol       = 1e-6,
    max_iter  = 2000,
    step_size = NULL,
    step_size_decay = 1.0,
    min_step_size   = NULL,
    relaxation      = 1,
    early_stop_count = 5) {
  
  set.seed(123)
  require(Matrix)
  require(Rcpp)
  
  # Data dimensions and matrix conversion
  n <- nrow(X); p <- ncol(X)
  if (inherits(X, "sparseMatrix") && p * n < 1e7) {
    X_dense  <- as.matrix(X)
    use_dense <- TRUE
  } else {
    X_dense  <- X
    use_dense <- FALSE
  }
  
  # Pre-compute matrices for efficiency
  XtX <- if (use_dense) crossprod(X_dense) else crossprod(X)
  Xty <- if (use_dense) crossprod(X_dense, y) else crossprod(X, y)
  L_quad <- (L_sparse + t(L_sparse)) / 2  # Symmetrize Laplacian
  L_quad <- as.matrix(L_quad)
  
  # Process group structure
  unique_groups     <- unique(groups)
  group_indices     <- lapply(unique_groups, function(g) which(groups == g))
  n_groups          <- length(group_indices)
  group_sizes       <- vapply(group_indices, length, 1L)
  group_indices_array <- lapply(group_indices, function(idx) idx - 1L)  # 0-based for C++
  
  # C++ implementation of group lasso proximal operator
  Rcpp::cppFunction('
  NumericVector group_lasso_prox_weighted_cpp(NumericVector beta,
                                              NumericVector lambda_vec,
                                              IntegerVector group_sizes,
                                              List group_indices) {
    int n_groups = group_sizes.size();
    NumericVector result = clone(beta);
    for (int g = 0; g < n_groups; ++g) {
      IntegerVector idx = group_indices[g];
      double l2 = 0.0;
      for (int j = 0; j < idx.size(); ++j) {
        l2 += beta[idx[j]] * beta[idx[j]];
      }
      double norm = sqrt(l2);
      double thr  = lambda_vec[g];
      if (norm > thr) {
        double shrink = 1.0 - thr / norm;
        for (int j = 0; j < idx.size(); ++j) {
          result[idx[j]] = beta[idx[j]] * shrink;
        }
      } else {
        for (int j = 0; j < idx.size(); ++j) {
          result[idx[j]] = 0.0;
        }
      }
    }
    return result;
  }')
  
  # Define proximal operators
  laplacian_prox <- function(beta, lambda_param) {
    A <- diag(p) + lambda_param * L_quad
    solve(A, beta)
  }
  
  # Gradient computation for logistic loss
  logistic_gradient <- function(Xbeta) {
    z    <- pmin(pmax(Xbeta, -30), 30)  # Clip for numerical stability
    p_hat <- 1 / (1 + exp(-z))
    if (use_dense) -crossprod(X_dense, y - p_hat) / n else -crossprod(X, y - p_hat) / n
  }
  
  # Generate lambda sequence if not provided
  Xt1 <- if (use_dense) crossprod(X_dense, rep(1, n)) else crossprod(X, rep(1, n))
  Xty_c <- Xty - 0.5 * Xt1
  group_norms <- vapply(group_indices, function(idx) sqrt(sum(Xty_c[idx]^2))/n, 0)
  
  if (is.null(lambda_seq)) {
    lambda_seq_list <- list()
    for (alpha_val in alpha_seq) {
      lambda_max <- if (alpha_val>0) max(group_norms)/alpha_val else max(group_norms)*10
      lambda_min <- if (n < p) 0.01*lambda_max else 1e-4*lambda_max
      lambda_seq_list[[as.character(alpha_val)]] <- exp(seq(log(lambda_max), log(lambda_min), length.out = nlambda))
    }
  }
  
  # Compute step size using power iteration for Lipschitz constant
  if (is.null(step_size)) {
    power_iter <- function(A, k=15) {
      v <- rnorm(ncol(A)); v <- v / sqrt(sum(v^2))
      for (i in 1:k) v <- A %*% v / sqrt(sum((A %*% v)^2))
      as.numeric(crossprod(v, A %*% v))
    }
    L_lip    <- power_iter(as.matrix(XtX)/n) / 4
    beta_val <- 1 / L_lip
    step_size <- min(1.9 * beta_val, beta_val)
  }
  if (is.null(min_step_size)) min_step_size <- max(0.1*step_size, 1e-10)
  
  # Calculate null model deviance
  eps      <- 1e-10
  p_null   <- mean(y)
  null_dev <- -2 * sum(y*log(p_null+eps) + (1-y)*log(1-p_null+eps))
  
  # Create grid of lambda-alpha combinations
  if (is.null(lambda_seq)) {
    grid <- do.call(rbind, lapply(alpha_seq, function(a){
      data.frame(lambda=lambda_seq_list[[as.character(a)]], alpha=a)
    }))
  } else {
    grid <- expand.grid(lambda=lambda_seq, alpha=alpha_seq)
  }
  grid <- grid[order(grid$alpha, decreasing=TRUE),]
  G    <- nrow(grid)
  res  <- vector("list", G)
  
  # Main optimization loop over lambda-alpha grid
  for (g in seq_len(G)) {
    lambda <- grid$lambda[g]; alpha <- grid$alpha[g]
    
    # Warm start from previous solution
    beta_k <- if (g>1 && grid$alpha[g]==grid$alpha[g-1]) res[[g-1]]$beta else rep(0, p)
    omega_g <- 0.5; omega_l <- 0.5
    z_g <- z_l <- beta_k
    cs  <- step_size; cr <- relaxation
    stall <- 0; prev_d <- Inf
    
    # Proximal gradient iterations
    for (iter in 1:max_iter) {
      beta_old <- beta_k
      x_k      <- omega_g*z_g + omega_l*z_l
      Xx_k     <- if (use_dense) X_dense %*% x_k else X %*% x_k
      grad_k   <- logistic_gradient(Xx_k)
      
      # Group lasso proximal update
      if (alpha>0) {
        prox_arg_g <- 2*x_k - z_g - cs*grad_k
        lambda_vec <- (lambda*alpha*cs/omega_g) * sqrt(group_sizes)
        prox_g     <- group_lasso_prox_weighted_cpp(prox_arg_g, lambda_vec, group_sizes, group_indices_array)
        z_g        <- z_g + cr*(prox_g - x_k)
      }
      
      # Laplacian proximal update
      if (alpha<1) {
        prox_arg_l <- 2*x_k - z_l - cs*grad_k
        prox_l     <- laplacian_prox(prox_arg_l, lambda*(1-alpha)*cs/omega_l)
        z_l        <- z_l + cr*(prox_l - x_k)
      }
      
      beta_k <- omega_g*z_g + omega_l*z_l
      
      # Check convergence
      delta  <- sqrt(sum((beta_k-beta_old)^2)) / (sqrt(sum(beta_old^2))+1e-10)
      if (delta < tol) break
      
      # Early stopping check
      if (abs(delta-prev_d) < tol*0.01) {
        stall <- stall+1; if (stall>=early_stop_count) break
      } else stall <- 0
      prev_d <- delta
      
      # Step size decay
      if (step_size_decay<1) {
        cs <- max(cs*step_size_decay, min_step_size)
        cs <- min(cs, 1.9*beta_val)
      }
    }
    
    # Calculate model deviance
    Xb    <- if (use_dense) X_dense %*% beta_k else X %*% beta_k
    z     <- pmin(pmax(Xb, -30), 30)
    p_hat <- 1/(1+exp(-z))
    dev   <- -2 * sum(y*log(p_hat+eps) + (1-y)*log(1-p_hat+eps))
    frac  <- (null_dev - dev)/null_dev
    
    res[[g]] <- list(beta=beta_k, dev=dev, frac=frac, 
                     params=list(lambda=lambda,alpha=alpha),
                     iter=iter, step=cs, conv=(delta<tol))
    
    if (g%%5==0 || g==G) cat(sprintf("[grid %d/%d] λ=%.3e, α=%.2f, frac=%.4f\n",g,G,lambda,alpha,frac))
  }
  
  # Select best model based on fraction of deviance explained
  best <- which.max(vapply(res, function(r) r$frac, 0))
  best_res <- res[[best]]
  
  cat("\nBest params -> λ=",best_res$params$lambda,", α=",best_res$params$alpha,"\n")
  cat("Null dev  :", sprintf("%.4f", null_dev), "\n")
  cat("Model dev :", sprintf("%.4f", best_res$dev),  "\n")
  cat("Frac expl :", sprintf("%.4f", best_res$frac), "\n")
  
  list(beta=best_res$beta,
       model_deviance = best_res$dev,
       frac_explained = best_res$frac,
       params         = best_res$params,
       all_results    = res)
}

#=============================================================
# 3. Stability Selection Integration
#=============================================================

#' Adaptive fitting function for stability selection
#' @param X Design matrix
#' @param y Response vector
#' @param q Number of variables to select
#' @param groups Group membership vector
#' @param L_sparse Laplacian matrix
#' @param lambda_seq Lambda sequence
#' @param alpha_seq Alpha sequence for mixing parameter
#' @param nlambda Number of lambda values
#' @return List with selected variables and model parameters
my_adaptive_fitfun <- function(X, y, q = 10,
                               groups = NULL,
                               L_sparse = NULL,
                               lambda_seq = NULL,
                               alpha_seq = seq(0, 1, length.out = 5),
                               nlambda = 10,
                               ...) {
  p <- ncol(X)
  
  # Set defaults
  if (is.null(groups)) groups <- seq_len(p)
  if (is.null(L_sparse)) L_sparse <- Diagonal(p, 0)
  
  # Fit model
  res_optim <- fit_adaptive_penalty_optim_noCV(
    X = X, y = y,
    L_sparse = L_sparse, groups = groups,
    lambda_seq = lambda_seq,
    alpha_seq = alpha_seq,
    nlambda = nlambda
  )
  
  # Define deviance calculation function
  calculate_deviance <- function(X, y, beta, selected_vars = NULL) {
    if (!is.null(selected_vars)) {
      X <- X[, selected_vars, drop = FALSE]
      beta <- beta[selected_vars]
    }
    
    Xbeta <- X %*% beta
    z_safe <- pmin(pmax(Xbeta, -30), 30)
    p_hat <- 1 / (1 + exp(-z_safe))
    
    eps <- 1e-10
    deviance <- -2 * sum(y * log(p_hat + eps) + (1 - y) * log(1 - p_hat + eps))
    
    p_null <- mean(y)
    null_deviance <- -2 * sum(y * log(p_null + eps) + (1 - y) * log(1 - p_null + eps))
    
    frac_explained <- (null_deviance - deviance) / null_deviance
    
    return(list(
      deviance = deviance,
      null_deviance = null_deviance,
      frac_explained = frac_explained
    ))
  }
  
  # Find best model for each alpha
  alpha_best_models <- list()
  
  for (alpha_val in alpha_seq) {
    alpha_models <- Filter(function(model) model$params$alpha == alpha_val, res_optim$all_results)
    if (length(alpha_models) == 0) next
    
    # Sort by lambda (descending)
    alpha_models <- alpha_models[order(sapply(alpha_models, function(model) model$params$lambda), decreasing = TRUE)]
    
    selected_vars <- integer(0)
    
    for (model in alpha_models) {
      beta <- model$beta
      
      # Add non-zero coefficients
      nonzero_idx <- which(abs(beta) > 1e-8)
      ordered_idx <- nonzero_idx[order(abs(beta[nonzero_idx]), decreasing = TRUE)]
      selected_vars <- c(selected_vars, setdiff(ordered_idx, selected_vars))
      
      if (length(selected_vars) >= q) {
        # Select top q variables by coefficient magnitude
        abs_beta_selected <- abs(beta[selected_vars])
        top_indices <- order(abs_beta_selected, decreasing = TRUE)[1:q]
        selected_vars <- selected_vars[top_indices]
        
        sel_vec <- rep(FALSE, p)
        sel_vec[selected_vars] <- TRUE
        
        dev_result <- calculate_deviance(X, y, beta, selected_vars)
        
        alpha_best_models[[as.character(alpha_val)]] <- list(
          beta = beta[selected_vars],
          selected = sel_vec,
          params = model$params,
          deviance = dev_result$deviance,
          frac_explained = dev_result$frac_explained
        )
        break
      }
    }
  }
  
  # Select optimal model based on deviance
  best_deviances <- sapply(alpha_best_models, function(model) model$deviance)
  best_alpha_key <- names(alpha_best_models)[which.min(best_deviances)]
  best_model <- alpha_best_models[[best_alpha_key]]
  
  return(list(
    best_beta = best_model$beta,
    selected = best_model$selected,
    path = NULL,
    params = best_model$params,
    deviance = best_model$deviance,
    frac_explained = best_model$frac_explained,
    all_alpha_models = alpha_best_models
  ))
}

#' Parallel stability selection with multiple PFER values
#' @param x Design matrix
#' @param y Response vector
#' @param fitfun Fitting function to use
#' @param args.fitfun Arguments for fitting function
#' @param cutoff Selection threshold
#' @param pfer_values Vector of PFER values to evaluate
#' @param B Number of subsampling iterations
#' @param sampling.type Sampling type ("SS" or "MB")
#' @param assumption Assumption type for error control
#' @param verbose Print progress
#' @param mc.cores Number of cores for parallel processing
#' @return List of results for each PFER value
my_stabsel_parallel <- function(x, y,
                                fitfun = my_adaptive_fitfun,
                                args.fitfun = list(),
                                cutoff = 0.6,   
                                pfer_values = NULL,
                                B = NULL,
                                sampling.type = c("SS", "MB"),
                                assumption = c("unimodal"),
                                verbose = TRUE,
                                mc.cores = 10) {
  set.seed(123)
  
  # Input validation and conversion
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  storage.mode(y) <- "numeric"
  
  if (!is.null(args.fitfun$L_sparse)) {
    if (is(args.fitfun$L_sparse, "list")) {
      args.fitfun$L_sparse <- args.fitfun$L_sparse$adj_matrix
    }
    args.fitfun$L_sparse <- as(args.fitfun$L_sparse, "dgCMatrix")
  }
  
  sampling.type <- match.arg(sampling.type)
  assumption <- match.arg(assumption)
  
  p <- ncol(x)
  n <- nrow(x)
  
  if (is.null(B)) {
    B <- if (sampling.type == "SS") 50 else 100
  }
  
  # Calculate q values for each PFER
  param_objs <- lapply(pfer_values, function(pfer) {
    stabsel_parameters(
      p = p,
      cutoff = cutoff,
      PFER = pfer,
      B = B,
      sampling.type = sampling.type,
      assumption = if (sampling.type == "SS") assumption else "none",
      verbose = verbose
    )
  })
  
  # Use maximum q value for initial selection
  q_max <- max(sapply(param_objs, function(x) round(x$q)))
  
  # Generate subsamples
  fold_mat <- subsample_half(n, B = B, sampling.type = sampling.type)
  n_sub <- ncol(fold_mat)
  
  # Function to fit single subsample
  single_fit <- function(b) {
    idx <- which(fold_mat[, b] == 1)
    Xb <- x[idx, , drop = FALSE]
    yb <- y[idx]
    fit_args <- c(list(X = Xb, y = yb, q = q_max), args.fitfun)
    res_fit <- do.call(fitfun, fit_args)
    list(
      beta = as.numeric(res_fit$best_beta),
      selected = as.numeric(res_fit$selected),
      params = res_fit$params
    )
  }
  
  # Parallel processing setup
  n_cores <- min(mc.cores, parallel::detectCores() - 1)
  cl <- parallel::makeCluster(n_cores)
  
  # Load required libraries in cluster
  parallel::clusterEvalQ(cl, {
    library(Matrix)
    library(stats)
    library(MASS)
    library(huge)
    library(glmnet)
    library(grpreg)
    NULL
  })
  
  # Export necessary objects to cluster
  parallel::clusterExport(cl, 
                          c("x", "y", "fold_mat", "q_max", "fitfun", "args.fitfun",
                            "fit_adaptive_penalty_optim_noCV", "my_adaptive_fitfun"), 
                          envir = environment())
  
  if (verbose) {
    cat(sprintf("\nParallel processing with %d cores...\n", n_cores))
  }
  
  # Run parallel processing
  result_list <- tryCatch({
    parallel::parLapply(cl, seq_len(n_sub), single_fit)
  }, error = function(e) {
    cat("Error in parallel processing:", conditionMessage(e), "\n")
    stop(e)
  }, finally = {
    parallel::stopCluster(cl)
  })
  
  # Create selection matrix and parameter list
  selection_record <- matrix(0, nrow = p, ncol = n_sub)
  params_list <- vector("list", n_sub)
  
  for (b in seq_len(n_sub)) {
    selection_record[(result_list[[b]]$selected == 1) %>% which, b] <- result_list[[b]]$beta
    params_list[[b]] <- result_list[[b]]$params
  }
  
  # Calculate results for each PFER value
  pfer_results <- list()
  for (i in seq_along(pfer_values)) {
    pfer <- pfer_values[i]
    param_obj <- param_objs[[i]]
    q_use <- round(param_obj$q)
    
    # Select top q_use variables for each subsample
    selection_record_pfer <- matrix(0, nrow = p, ncol = n_sub)
    for (b in seq_len(n_sub)) {
      sorted_idx <- order(abs(selection_record[, b]), decreasing = TRUE)
      selection_record_pfer[sorted_idx[1:q_use], b] <- 1
    }
    
    # Calculate selection probability
    phat <- rowMeans(selection_record_pfer)
    selected_ids <- which(phat > param_obj$cutoff)
    
    # Extract parameter information
    params_list <- lapply(result_list, function(res) res$params)
    
    # Calculate average lambda and alpha
    avg_lambda <- mean(sapply(params_list, function(par) if (!is.null(par$lambda)) par$lambda else NA), na.rm = TRUE)
    avg_alpha  <- mean(sapply(params_list, function(par) if (!is.null(par$alpha)) par$alpha else NA), na.rm = TRUE)
    sd_lambda <- sd(sapply(params_list, function(par) if (!is.null(par$lambda)) par$lambda else NA), na.rm = TRUE)
    sd_alpha  <- sd(sapply(params_list, function(par) if (!is.null(par$alpha)) par$alpha else NA), na.rm = TRUE)
    
    pfer_results[[as.character(pfer)]] <- list(
      param_obj = param_obj,
      phat = phat,
      selected = selected_ids,
      max = max(phat),
      q = param_obj$q,
      cutoff = param_obj$cutoff,
      PFER = param_obj$PFER,
      specifiedPFER = param_obj$specifiedPFER,
      B = param_obj$B,
      sampling.type = param_obj$sampling.type,
      assumption = param_obj$assumption,
      n_subsamples = n_sub,
      selection_record = selection_record_pfer,
      params = params_list,
      params_avg = list(
        lambda_mean = avg_lambda,
        lambda_sd = sd_lambda,
        alpha_mean = avg_alpha,
        alpha_sd = sd_alpha
      )
    )
  }
  
  return(pfer_results)
}

#' Generate subsample indices for stability selection
#' @param n Number of observations
#' @param B Number of subsample pairs
#' @param sampling.type Type of sampling ("SS" or "MB")
#' @return Matrix of subsample indicators
subsample_half <- function(n, B=50, sampling.type=c("SS","MB")) {
  sampling.type <- match.arg(sampling.type)
  fold_mat <- matrix(0, nrow=n, ncol=ifelse(sampling.type=="SS", 2*B, B))
  
  cat(sprintf("[subsample_half] sampling.type='%s', B=%d => producing %d folds\n",
              sampling.type, B, ncol(fold_mat)))
  
  set.seed(123)
  if (sampling.type=="SS") {
    # Complementary pairs stability selection
    for (b in seq_len(B)) {
      chosen <- sample.int(n, size=floor(n/2))
      fold_mat[chosen, b] <- 1
      not_chosen <- setdiff(seq_len(n), chosen)
      fold_mat[not_chosen, b + B] <- 1
    }
  } else {
    # Standard subsampling
    for (b in seq_len(B)) {
      chosen <- sample.int(n, size=floor(n/2))
      fold_mat[chosen, b] <- 1
    }
  }
  fold_mat
}

#=============================================================
# 4. Baseline Methods for Comparison
#=============================================================

#' Group Lasso with grpreg package
grpreg.grouplasso <- function(x, y, q, type = c("conservative", "anticonservative"),
                              group = NULL, ...) {
  if (!requireNamespace("grpreg", quietly = TRUE))
    stop("Package ", sQuote("grpreg"), " needed but not available")
  
  if (is.data.frame(x)) {
    message("Note: ", sQuote("x"), " is coerced to a model matrix without intercept")
    x <- model.matrix(~. - 1, x)
  }
  
  if ("lambda" %in% names(list(...)))
    stop("It is not permitted to specify the penalty parameter ",
         sQuote("lambda"), " for group lasso when used with stability selection.")
  
  if (is.null(group))
    stop("A 'group' parameter must be provided with length equal to the number of columns in x")
  
  if (length(group) != ncol(x))
    stop("The 'group' vector must have length equal to the number of columns in x")
  
  type <- match.arg(type)
  
  # Fit group LASSO model
  if (type == "conservative") {
    fit <- suppressWarnings(grpreg::grpreg(x, y, group = group,
                                           penalty = "grLasso",
                                           dfmax = q, ...))
  } else {
    fit <- grpreg::grpreg(x, y, group = group,
                          penalty = "grLasso",
                          dfmax = q - 1, ...)
  }
  
  # Get non-zero coefficients from last lambda
  coef_matrix <- fit$beta
  last_coef <- coef_matrix[, ncol(coef_matrix)]
  selected <- last_coef[-1] != 0  # Remove intercept
  
  ret <- logical(ncol(x))
  ret[selected] <- TRUE
  names(ret) <- colnames(x)
  
  beta_path <- as.matrix(fit$beta[-1, ])
  sequence <- beta_path != 0
  
  return(list(selected = ret, path = sequence))
}

#' Elastic Net with glmnet package
glmnet.elasticnet <- function (x, y, q, type = c("conservative", "anticonservative"), 
                               ...) {
  if (!requireNamespace("glmnet", quietly = TRUE)) 
    stop("Package ", sQuote("glmnet"), " needed but not available")
  
  if (is.data.frame(x)) {
    message("Note: ", sQuote("x"), " is coerced to a model matrix without intercept")
    x <- model.matrix(~. - 1, x)
  }
  
  if ("lambda" %in% names(list(...))) 
    stop("It is not permitted to specify the penalty parameter ", 
         sQuote("lambda"), " for lasso when used with stability selection.")
  
  type <- match.arg(type)
  
  if (type == "conservative") 
    fit <- suppressWarnings(glmnet::glmnet(x, y, pmax = q, alpha = 0.5, ...))
  else
    fit <- glmnet::glmnet(x, y, dfmax = q - 1, alpha = 0.5, ...)
  
  selected <- predict(fit, type = "nonzero")
  selected <- selected[[length(selected)]]
  ret <- logical(ncol(x))
  ret[selected] <- TRUE
  names(ret) <- colnames(x)
  
  cf <- fit$beta
  sequence <- as.matrix(cf != 0)
  
  return(list(selected = ret, path = sequence))
}
