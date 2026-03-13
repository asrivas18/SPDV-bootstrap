# =============================================================================
# 2. Economic Inequality: CPS 2022 Household Income (ALL 4 Methods)
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

# Point estimates using SPDV functions
source("R/spdv_functions.R")
s2_inc <- compute_spdv(inc)
kappa_inc <- kurtosis_hat(inc)
se_inc <- se_spdv(inc)

cat(sprintf("Point estimates: s²=%.2e, SE=%.2e, κ̂=%.2f\n", s2_inc, se_inc, kappa_inc))

# =============================================================================
# METHOD 1: CHI-SQUARE (classical)
# =============================================================================
set.seed(20220301)
chi_ci_inc <- chi_square_ci(inc)
chi_length_inc <- chi_ci_inc[2] - chi_ci_inc[1]

cat(sprintf("Chi-square 95%% CI: [%.2e, %.2e] (length=%.2e)\n", 
            chi_ci_inc[1], chi_ci_inc[2], chi_length_inc))

# =============================================================================
# METHOD 2: NORMAL APPROXIMATION
# =============================================================================
normal_ci_inc <- normal_ci(inc)
normal_length_inc <- normal_ci_inc[2] - normal_ci_inc[1]

cat(sprintf("Normal 95%% CI: [%.2e, %.2e] (length=%.2e)\n", 
            normal_ci_inc[1], normal_ci_inc[2], normal_length_inc))

# =============================================================================
# METHOD 3: PERCENTILE BOOTSTRAP
# =============================================================================
pct_ci_inc <- percentile_ci(inc, B = 5000)
pct_length_inc <- pct_ci_inc[2] - pct_ci_inc[1]

cat(sprintf("Percentile Boot 95%% CI: [%.2e, %.2e] (length=%.2e)\n", 
            pct_ci_inc[1], pct_ci_inc[2], pct_length_inc))

# =============================================================================
# METHOD 4: STUDENTIZED BOOTSTRAP (SPDV)
# =============================================================================
stud_res_inc <- studentized_bootstrap(inc, B = 5000, parallel = FALSE)
stud_ci_inc <- stud_res_inc$ci_95
stud_length_inc <- stud_ci_inc[2] - stud_ci_inc[1]
stud_runtime_inc <- round(stud_res_inc$runtime_s, 1)

cat(sprintf("Studentized Boot 95%% CI: [%.2e, %.2e] (length=%.2e)\n", 
            stud_ci_inc[1], stud_ci_inc[2], stud_length_inc))
cat(sprintf("Studentized runtime: %.1fs\n", stud_runtime_inc))

cat(rep("=", 60), "\n")


# =============================================================================
# SAVE RESULTS
# =============================================================================
cps_results <- data.frame(
  Method = c("Chi-square", "Normal", "Pctl Boot", "Stud Boot"),
  s2 = rep(s2_inc, 4),
  CI = c(sprintf("[%.2e, %.2e]", chi_ci_inc[1], chi_ci_inc[2]),
         sprintf("[%.2e, %.2e]", normal_ci_inc[1], normal_ci_inc[2]),
         sprintf("[%.2e, %.2e]", pct_ci_inc[1], pct_ci_inc[2]),
         sprintf("[%.2e, %.2e]", stud_ci_inc[1], stud_ci_inc[2])),
  Length = c(chi_length_inc, normal_length_inc, pct_length_inc, stud_length_inc)/1e9
)

write.csv(cps_results, "results/cps_table5.csv", row.names = FALSE)
cat("\n✓ CPS results saved: results/cps_table5.csv\n")
cat("🎉 COMPLETE CPS analysis + Table 5 ready!\n")
