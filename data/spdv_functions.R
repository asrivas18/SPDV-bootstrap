# =============================================================================
# spdv_functions.R
# Core functions for Scalable Studentized Bootstrap Variance Inference
# Manuscript: "Scalable Studentized Bootstrap Variance Inference via Linear-Time Pairwise Variance Representation"
# Journal: Journal of Statistical Computation and Simulation (JSCS)
# =============================================================================

# ----------------------------------------------------------------------------
# Fast O(n) variance (SPDV = s² via algebraic identity)
# ----------------------------------------------------------------------------
spdv_fast <- function(x) {
  if (length(x) < 2) stop("Sample size must be at least 2")
  var(x)  # R's var() is unbiased: sum((x - mean(x))^2) / (n-1)
}

# ----------------------------------------------------------------------------
# Empirical influence-function standard error (factor 4 from 4τ²)
# ----------------------------------------------------------------------------
influence_se <- function(x) {
  n <- length(x)
  if (n < 2) stop("Sample size must be at least 2")
  
  xbar <- mean(x)
  s2   <- var(x)
  
  # First-order projection values: ψ̂_i = [(x_i - x̄)^2 - s²]/2
  psi_hat <- ((x - xbar)^2 - s2) / 2
  
  # Asymptotic SE = sqrt(4/n * Var(ψ̂))
  se_hat <- sqrt(4 / n * var(psi_hat))
  
  # Numerical safeguard
  if (se_hat < 1e-10 || !is.finite(se_hat)) {
    warning("Very small or non-finite SE detected — possible extreme concentration")
    se_hat <- max(se_hat, 1e-8)
  }
  
  return(se_hat)
}

# ----------------------------------------------------------------------------
# Studentized bootstrap-t confidence interval
# Returns list with variance estimate, CI, valid pivot count
# Supports parallel execution via foreach/doParallel
# ----------------------------------------------------------------------------
studentized_bootstrap <- function(x, B = 50000, alpha = 0.05, parallel = TRUE, n_cores = NULL) {
  n <- length(x)
  if (n < 2) stop("Sample size must be at least 2")
  
  theta_hat <- spdv_fast(x)          # s²
  se_hat    <- influence_se(x)       # Original-sample SE
  
  if (parallel) {
    if (is.null(n_cores)) n_cores <- max(1, parallel::detectCores() - 1)
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
    
    boot_results <- foreach(b = 1:B, .combine = rbind, .packages = "dplyr") %dopar% {
      xb <- sample(x, n, replace = TRUE)
      theta_b <- spdv_fast(xb)
      se_b    <- influence_se(xb)
      data.frame(theta_b = theta_b, se_b = se_b)
    }
  } else {
    boot_results <- data.frame(theta_b = numeric(B), se_b = numeric(B))
    set.seed(12345 + sample(1:1000, 1))  # Vary slightly per run if needed
    for (b in 1:B) {
      xb <- sample(x, n, replace = TRUE)
      boot_results$theta_b[b] <- spdv_fast(xb)
      boot_results$se_b[b]    <- influence_se(xb)
    }
  }
  
  # Valid pivots only (avoid division by zero / Inf)
  valid <- boot_results$se_b > 1e-10 & is.finite(boot_results$se_b)
  valid_count <- sum(valid)
  
  if (valid_count == 0) stop("No valid bootstrap SEs — check data for extreme concentration")
  
  t_star <- (boot_results$theta_b[valid] - theta_hat) / boot_results$se_b[valid]
  
  t_low  <- quantile(t_star, alpha/2, na.rm = TRUE)
  t_high <- quantile(t_star, 1 - alpha/2, na.rm = TRUE)
  
  ci_var <- c(theta_hat - t_high * se_hat, theta_hat - t_low * se_hat)
  
  list(
    n             = n,
    variance      = theta_hat,
    se            = se_hat,
    ci_variance   = ci_var,
    ci_root       = sqrt(2 * ci_var),
    valid_pivots  = valid_count,
    total_B       = B,
    skipped_pct   = round(100 * (B - valid_count) / B, 2)
  )
}