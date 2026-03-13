# ============================================================================
# bootstrap_functions.R
# Core functions for SPDV studentized bootstrap variance inference
# Manuscript: Scalable Studentized Bootstrap Inference for Variance via Pairwise Difference Representation
# Journal: Computational Statistics & Data Analysis (CSDA)
# ============================================================================

# ----------------------------------------------------------------------------
# Function: compute_spdv
# Computes Sample Pairwise Difference Variance (SPDV) — algebraically = var(x)
# ----------------------------------------------------------------------------
compute_spdv <- function(x) {
  n <- length(x)
  if (n < 2) return(NA)
  # Direct pairwise (O(n²) — for illustration only; not used in bootstrap)
  # mean(outer(x, x, "-")^2) / 2
  # Instead, use O(n) identity
  var(x)  # unbiased sample variance s² = SPDV exactly
}

# ----------------------------------------------------------------------------
# Function: compute_influence_se
# Computes asymptotic SE of variance using unscaled influence function
# ψ_i = (X_i - \bar{X})^2 - s²
# SE = sqrt( Var(ψ) / n )
# ----------------------------------------------------------------------------
compute_influence_se <- function(x) {
  n <- length(x)
  if (n < 2) return(NA)
  xbar <- mean(x)
  s2   <- var(x)
  psi  <- (x - xbar)^2 - s2
  se   <- sqrt(var(psi) / n)
  # Prevent negative/zero variance (rare numerical issue)
  se   <- pmax(se, 1e-10)
  se
}

# ----------------------------------------------------------------------------
# Function: studentized_bootstrap
# Main function: Runs studentized bootstrap for variance
# Returns: list with CI bounds, raw bootstrap variances, pivots, and SEs
# ----------------------------------------------------------------------------
studentized_bootstrap <- function(x, B = 5000, alpha = 0.05, parallel = TRUE, n_cores = NULL) {
  n <- length(x)
  if (n < 2) stop("Sample size must be at least 2")

  s2 <- var(x)  # Original estimate (SPDV = s²)
  se_orig <- compute_influence_se(x)  # Original-sample SE

  # Bootstrap loop
  if (parallel && n_cores > 1) {
    if (is.null(n_cores)) n_cores <- max(1, parallel::detectCores() - 1)
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    doRNG::registerDoRNG(12345)

    boot_res <- foreach(b = 1:B, .combine = rbind, .packages = c("stats")) %dopar% {
      xb <- sample(x, n, replace = TRUE)
      s2_b <- var(xb)
      se_b <- compute_influence_se(xb)
      t_b  <- (s2_b - s2) / se_b
      data.frame(s2_b = s2_b, se_b = se_b, t_b = t_b)
    }

    parallel::stopCluster(cl)
  } else {
    # Sequential fallback
    boot_res <- data.frame(s2_b = numeric(B), se_b = numeric(B), t_b = numeric(B))
    for (b in 1:B) {
      xb <- sample(x, n, replace = TRUE)
      s2_b <- var(xb)
      se_b <- compute_influence_se(xb)
      t_b  <- (s2_b - s2) / se_b
      boot_res[b, ] <- c(s2_b, se_b, t_b)
    }
  }

  # Percentile CI (for comparison)
  pctl_ci <- quantile(boot_res$s2_b, c(alpha/2, 1 - alpha/2), na.rm = TRUE)
  names(pctl_ci) <- c("lower", "upper")

  # Studentized CI (bootstrap-t inversion)
  valid <- is.finite(boot_res$t_b) & boot_res$se_b > 0
  t_star <- boot_res$t_b[valid]

  if (length(t_star) < 10) {
    warning("Too few valid studentized pivots (<10) — possible instability")
    stud_ci <- c(NA, NA)
  } else {
    t_low  <- quantile(t_star, alpha/2, na.rm = TRUE)
    t_high <- quantile(t_star, 1 - alpha/2, na.rm = TRUE)
    stud_ci <- c(s2 - t_high * se_orig, s2 - t_low * se_orig)
    names(stud_ci) <- c("lower", "upper")
  }

  # Return results
  list(
    original_s2    = s2,
    original_se    = se_orig,
    boot_variances = boot_res$s2_b,
    boot_pivots    = boot_res$t_b,
    pctl_ci        = pctl_ci,
    stud_ci        = stud_ci,
    valid_pivots   = sum(valid),
    skipped_frac   = mean(!valid)
  )
}

# ----------------------------------------------------------------------------
# Example usage (for testing)
# ----------------------------------------------------------------------------
if (FALSE) {
  set.seed(12345)
  x <- rnorm(50)
  res <- studentized_bootstrap(x, B = 1000, parallel = TRUE)
  print(res$pctl_ci)
  print(res$stud_ci)
}