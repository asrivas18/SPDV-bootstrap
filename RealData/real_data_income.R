# real_data_income.R (Economic Inequality — Illustrative)


library(tidyverse)

source("bootstrap_functions.R")

# Synthetic heavy-tailed income data (lognormal, mean ≈ $60k, heavy right tail)
set.seed(123)
n <- 10000
income <- exp(rnorm(n, mean = log(60000), sd = 0.8))  # Approximate US income distribution

cat("Synthetic income data: n =", n, "\n")
result <- studentized_bootstrap(income, B = 50000)

cat("Sample variance s² ≈", round(result$variance, 0), "\n")
cat("Root pairwise difference ≈", round(result$root_pairwise, 0), "\n")
cat("95% CI for variance:", round(result$ci_variance, 0), "\n")
cat("95% CI for root pairwise:", round(result$ci_root, 0), "\n")

# Illustrative interpretation
cat("\nInterpretation: On average, two randomly selected households differ by ≈ $",
    round(result$root_pairwise, 0), " in income (squared units via pairwise view).\n")