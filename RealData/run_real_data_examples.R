# =============================================================================
# run_real_data_examples.R
# Real-data applications for SPDV studentized bootstrap:
#   - BMI heterogeneity (Louisiana BRFSS 2022)
#   - Economic inequality (US CPS 2022 income)
#
# Requires:
#   - R/spdv_functions.R
#   - data/brfss_bmi_subset.rds
#   - data/cps_income_subset.rds
# =============================================================================
rm(list = ls())

# ---- Load core functions ----
source("R/spdv_functions.R")

# ---- Helper: pretty printing ----
print_ci <- function(label, s2, ci_var, ci_root = NULL, units = "") {
  cat("\n", label, "\n", sep = "")
  cat(" s^2 estimate: ", formatC(s2, digits = 4, format = "g"), " ", units, "\n", sep = "")
  cat(" 95% CI for variance: [", formatC(ci_var[1], digits = 4, format = "g"), ", ",
      formatC(ci_var[2], digits = 4, format = "g"), "] ", units, "\n", sep = "")
  if (!is.null(ci_root)) {
    cat(" 95% CI for root pairwise difference: [", formatC(ci_root[1], digits = 4, format = "g"), ", ",
        formatC(ci_root[2], digits = 4, format = "g"), "]", "\n", sep = "")
  }
}

# =============================================================================
# 1. BMI heterogeneity (Louisiana BRFSS 2022)
# =============================================================================
cat("============================================================\n")
cat("1. Public Health: BMI Heterogeneity (BRFSS Louisiana 2022)\n")
cat("============================================================\n")

brfss_path <- "data/brfss_bmi_subset.rds"
if (!file.exists(brfss_path)) {
  stop("File not found: ", brfss_path, "\nSee data/README-data.md for instructions.")
}

brfss <- readRDS(brfss_path)
# Expect either numeric vector or data.frame with 'bmi' column
if (is.data.frame(brfss)) {
  if (!"bmi" %in% names(brfss)) stop("brfss_bmi_subset.rds must contain column 'bmi'")
  bmi <- brfss$bmi
} else {
  bmi <- as.numeric(brfss)
}
bmi <- bmi[is.finite(bmi) & bmi > 0]  # Safety: positive finite BMI
n_bmi <- length(bmi)
cat("Sample size (complete cases): n =", n_bmi, "\n")

# SPDV = s^2
s2_bmi <- var(bmi)  # or compute_spdv(bmi) if implemented
root_pairwise_bmi <- sqrt(2 * s2_bmi)

set.seed(20220301)
system.time({
  res_bmi <- studentized_bootstrap(bmi, B = 50000, alpha = 0.05, parallel = TRUE)
}) -> time_bmi

ci_bmi_var  <- res_bmi$stud_ci
ci_bmi_root <- sqrt(2 * ci_bmi_var)

print_ci("BMI example (SPDV studentized bootstrap)",
         s2_bmi, ci_bmi_var, ci_bmi_root, units = "")

cat(" Runtime (wall-clock): ", round(time_bmi["elapsed"], 1), " seconds\n", sep = "")
cat(" Valid studentized pivots: ", res_bmi$valid_pivots, " / ", res_bmi$B, "\n", sep = "")

# =============================================================================
# 2. Economic inequality (US CPS 2022 income)
# =============================================================================
cat("\n============================================================\n")
cat("2. Economic Inequality: CPS 2022 Household Income\n")
cat("============================================================\n")

cps_path <- "data/cps_income_subset.rds"
if (!file.exists(cps_path)) {
  stop("File not found: ", cps_path, "\nSee data/README-data.md for instructions.")
}

cps <- readRDS(cps_path)
# Expect data.frame with 'income' column (USD)
if (is.data.frame(cps)) {
  if (!"income" %in% names(cps)) stop("cps_income_subset.rds must contain column 'income'")
  inc <- cps$income
} else {
  inc <- as.numeric(cps)
}
inc <- inc[is.finite(inc) & inc > 0]  # Positive finite income
n_inc <- length(inc)
cat("Sample size (processed): n =", n_inc, "\n")

s2_inc <- var(inc)  # or compute_spdv(inc)

set.seed(20220302)
system.time({
  res_inc <- studentized_bootstrap(inc, B = 50000, alpha = 0.05, parallel = TRUE)
}) -> time_inc

ci_inc_var <- res_inc$stud_ci

print_ci("CPS income example (SPDV studentized bootstrap)",
         s2_inc, ci_inc_var, ci_root = NULL, units = "USD^2")

cat(" Runtime (wall-clock): ", round(time_inc["elapsed"], 1), " seconds\n", sep = "")
cat(" Valid studentized pivots: ", res_inc$valid_pivots, " / ", res_inc$B, "\n", sep = "")

cat("\nAll real-data examples completed.\n")