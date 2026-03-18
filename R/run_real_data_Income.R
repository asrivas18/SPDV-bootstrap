# =============================================================================
# 2. CPS 2022 Household Income Analysis – Matches Table \ref{tab:real-data}
# =============================================================================
cat("------------------------------------------------------------\n")
cat("2. Economic Inequality: CPS 2022 Household Income\n")
cat("------------------------------------------------------------\n\n")

# Path to preprocessed subset (from process_cps_2022.R)
cps_path <- "data/cps_income_subset_60k.rds"

if (!file.exists(cps_path)) {
  stop("File not found: ", cps_path, "\nSee data/README-data.md for instructions.")
}

# Load preprocessed data
cps <- readRDS(cps_path)

# Extract income column
if (is.data.frame(cps)) {
  if (!"Total_HH_Income" %in% names(cps)) {
    stop("cps_income_subset_60k.rds must contain column 'Total_HH_Income'")
  }
  inc <- cps$Total_HH_Income
} else {
  inc <- as.numeric(cps)
}

# Safety filter: positive finite values
inc <- inc[is.finite(inc) & inc > 0]
n_inc <- length(inc)

cat(sprintf("Sample size (processed): n = %d\n", n_inc))

# Point estimates
s2_inc    <- compute_spdv(inc)
kappa_inc <- moments::kurtosis(inc)
se_inc    <- se_spdv(inc)

cat(sprintf("Point estimates: s^2 = %.2e, SE = %.2e, kurtosis_hat = %.2f\n",
            s2_inc, se_inc, kappa_inc))

# METHOD 1: Chi-square interval (assumes normality – for comparison only)
chi_ci_inc <- chi_square_ci(inc)
chi_length_inc <- chi_ci_inc[2] - chi_ci_inc[1]
cat(sprintf("Chi-square 95%% CI: [%.2e, %.2e] (length = %.2e)\n",
            chi_ci_inc[1], chi_ci_inc[2], chi_length_inc))

# METHOD 2: Normal approximation (using influence-function SE)
normal_ci_inc <- normal_ci(inc)
normal_length_inc <- normal_ci_inc[2] - normal_ci_inc[1]
cat(sprintf("Normal 95%% CI: [%.2e, %.2e] (length = %.2e)\n",
            normal_ci_inc[1], normal_ci_inc[2], normal_length_inc))

# METHOD 3: Percentile bootstrap (B=5000 for speed in demo)
set.seed(20220302)
pct_ci_inc <- percentile_ci(inc, B = 5000)
pct_length_inc <- pct_ci_inc[2] - pct_ci_inc[1]
cat(sprintf("Percentile Boot 95%% CI: [%.2e, %.2e] (length = %.2e)\n",
            pct_ci_inc[1], pct_ci_inc[2], pct_length_inc))

# METHOD 4: Studentized bootstrap (SPDV) – full B=50,000
cat("\nRunning studentized bootstrap (B=50,000)...\n")
system.time({
  stud_res_inc <- studentized_bootstrap(inc, B = 50000, parallel = TRUE)
}) -> time_stud

stud_ci_inc <- stud_res_inc$ci_variance
stud_length_inc <- stud_ci_inc[2] - stud_ci_inc[1]

cat(sprintf("Studentized Boot 95%% CI: [%.2e, %.2e] (length = %.2e)\n",
            stud_ci_inc[1], stud_ci_inc[2], stud_length_inc))
cat(sprintf("Runtime (wall-clock): %.1f seconds\n", time_stud["elapsed"]))
cat(sprintf("Valid studentized pivots: %d / %d (%.1f%% skipped)\n",
            stud_res_inc$valid_pivots, stud_res_inc$B, stud_res_inc$skipped_pct))

# =============================================================================
# Save results to CSV for manuscript Table \ref{tab:real-data} reproduction
# =============================================================================
results_table <- data.frame(
  Method = c("Chi-square", "Normal", "Percentile Boot", "Studentized Boot"),
  s2     = rep(round(s2_inc, 2), 4),
  CI_lower = c(chi_ci_inc[1], normal_ci_inc[1], pct_ci_inc[1], stud_ci_inc[1]),
  CI_upper = c(chi_ci_inc[2], normal_ci_inc[2], pct_ci_inc[2], stud_ci_inc[2]),
  Length = c(chi_length_inc, normal_length_inc, pct_length_inc, stud_length_inc),
  Kurtosis = c(NA, NA, NA, round(kappa_inc, 2))
)

# Round scientific notation for readability
results_table <- results_table %>%
  mutate(across(c(CI_lower, CI_upper, Length),
                ~ formatC(., format = "e", digits = 2)))

cat("\nResults summary for Table \\ref{tab:real-data}:\n")
print(results_table)

# Save to CSV
write.csv(results_table, "data/cps_income_results.csv", row.names = FALSE)
cat("\nResults saved to: data/cps_income_results.csv\n")