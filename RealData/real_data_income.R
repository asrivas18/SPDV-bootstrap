# ============================================================================
# CPS ASEC 2022 Household Income Extraction & Studentized Bootstrap Demo
# Manuscript: Scalable Studentized Bootstrap Variance Inference via Linear-Time Pairwise Variance Representation
# Journal: Journal of Statistical Computation and Simulation (JSCS)
# Creates: data/cps_income_2022.csv (~60,000 rows) + computes studentized CI
# ============================================================================

rm(list = ls())

# ---------------------------------------------------------------------------
# 1. Install & Load Packages
# ---------------------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(haven, dplyr, data.table)

# ---------------------------------------------------------------------------
# 2. Load CPS SAS Household File
# ---------------------------------------------------------------------------
file_path <- "cps_income.sas7bdat"  # Adjust path as needed

#if (!file.exists(file_path)) {
#  stop("CPS file not found. Download from:\n",
#       "https://www.census.gov/data/datasets/2022/demo/cps/cps-asec-2022.html\n",
#       "Look for hhpub22.sas7bdat or asecpub22.sas7bdat")
#}

cat("Reading CPS 2022 household file...\n")
cps <- read_sas(file_path)

# ---------------------------------------------------------------------------
# 3. Extract & Clean Household Income (HTOTVAL)
# ---------------------------------------------------------------------------
cps_income <- cps %>%
  transmute(
    Household_ID   = as.character(SERIAL),     # Unique household ID
    Total_HH_Income = as.numeric(HHINCOME),    # Total household income (main var)
    HH_Weight      = as.numeric(ASECWT)      # Household weight (no scaling needed)
  ) %>%
  filter(!is.na(Total_HH_Income) & Total_HH_Income >= 0)  # Positive income only

cat("Rows after cleaning:", nrow(cps_income), "\n")

# ---------------------------------------------------------------------------
# 4. Reproducible Subsampling (~60,000 rows, matches manuscript)
# ---------------------------------------------------------------------------
set.seed(2022)
target_n <- 60000

if (nrow(cps_income) > target_n) {
  cps_income <- cps_income %>% slice_sample(n = target_n)
  cat("Subsampled to", nrow(cps_income), "households.\n")
} else {
  cat("Using full cleaned dataset (", nrow(cps_income), "rows).\n")
}

# ---------------------------------------------------------------------------
# 5. Save Processed Dataset
# ---------------------------------------------------------------------------
dir.create("data", showWarnings = FALSE)
fwrite(cps_income, "data/cps_income_2022.csv")
saveRDS(cps_income, "data/cps_income_2022.rds")
cat("Saved: data/cps_income_2022.csv and .rds\n")

# Quick summary statistics (matches manuscript description)
cat("\nSummary of Household Income (USD):\n")
print(summary(cps_income$Total_HH_Income))
cat("Empirical kurtosis:", moments::kurtosis(cps_income$Total_HH_Income), "\n")

# ---------------------------------------------------------------------------
# 6. SPDV Variance Estimator (for reference)
# ---------------------------------------------------------------------------
spdv <- function(x) {
  mean(x^2) - mean(x)^2
}

spdv_est <- spdv(cps_income$Total_HH_Income)
cat("\nSPDV variance estimate:", spdv_est, "\n")

# ---------------------------------------------------------------------------
# 7. Correct Studentized Bootstrap-t Interval (matches manuscript method)
# ---------------------------------------------------------------------------
studentized_bootstrap_ci <- function(x, B = 50000, alpha = 0.05) {
  n <- length(x)
  theta_hat <- spdv(x)  # or var(x) — same

  # Original sample projection values
  xbar <- mean(x)
  psi_hat <- (x - xbar)^2 - theta_hat
  se_hat <- sqrt(4 * var(psi_hat) / n)  # Factor 4 from 4τ²

  # Bootstrap loop
  boot_theta <- numeric(B)
  boot_se    <- numeric(B)

  for (b in 1:B) {
    xb <- sample(x, n, replace = TRUE)
    boot_theta[b] <- spdv(xb)

    xbar_b <- mean(xb)
    psi_b  <- (xb - xbar_b)^2 - boot_theta[b]
    boot_se[b] <- sqrt(4 * var(psi_b) / n)
  }

  # Studentized pivots (only valid resamples)
  valid <- boot_se > 1e-10 & is.finite(boot_se)
  if (sum(valid) == 0) stop("All bootstrap SEs invalid — check data")

  t_star <- (boot_theta[valid] - theta_hat) / boot_se[valid]

  t_low  <- quantile(t_star, alpha/2, na.rm = TRUE)
  t_high <- quantile(t_star, 1 - alpha/2, na.rm = TRUE)

  ci_low  <- theta_hat - t_high * se_hat
  ci_high <- theta_hat - t_low  * se_hat

  c(lower = ci_low, upper = ci_high)
}

# Run demo (B=50,000 — will take ~10–30 min single-threaded)
set.seed(123)
system.time({
  ci <- studentized_bootstrap_ci(cps_income$Total_HH_Income, B = 50000)
})

cat("\n95% Studentized Bootstrap CI for variance:\n")
print(ci)
cat("\nNote: Manuscript reports approximate [1.9e9, 2.3e9] — values may vary slightly due to random seed and subsampling.\n")
