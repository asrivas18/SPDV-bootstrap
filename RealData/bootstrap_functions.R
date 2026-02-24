
### 2. bootstrap_functions.R (Core — used by both applications)
# Core functions for scalable SPDV bootstrap

library(tidyverse)

# Fast O(n) variance (equivalent to SPDV via identity)
spdv_fast <- function(x) {
  n <- length(x)
  if (n < 2) return(NA)
  var(x)  # R's var() is unbiased: (n-1)/n adjustment handled automatically
}

# Studentized bootstrap with empirical fourth-moment SE
studentized_bootstrap <- function(x, B = 50000, conf = 0.95) {
  n <- length(x)
  if (n < 2) stop("Sample size too small")
  
  s2 <- var(x)
  alpha <- 1 - conf
  
  # Empirical influence-function SE
  psi_hat <- ((x - mean(x))^2 - s2) / 2
  se_hat <- sqrt(4 / n^2 * sum(psi_hat^2))
  if (se_hat == 0) se_hat <- 1e-8  # Avoid division by zero
  
  # Bootstrap replicates
  boot_s2 <- numeric(B)
  set.seed(12345)  # Fixed for reproducibility
  for (b in seq_len(B)) {
    xb <- sample(x, n, replace = TRUE)
    boot_s2[b] <- var(xb)
  }
  
  # Studentized pivots
  t_star <- (boot_s2 - s2) / se_hat
  
  # Quantiles
  q <- quantile(t_star, probs = c(alpha/2, 1 - alpha/2))
  
  # CI for variance
  ci_var <- c(s2 - q[2] * se_hat, s2 - q[1] * se_hat)
  
  # Root pairwise CI
  ci_root <- sqrt(2 * ci_var)
  
  list(
    n = n,
    variance = s2,
    root_pairwise = sqrt(2 * s2),
    ci_variance = ci_var,
    ci_root = ci_root
  )
}