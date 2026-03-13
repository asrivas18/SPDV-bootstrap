## =============================================================================
## run_real_data_examples.R
## Real-data examples: BRFSS BMI + CPS Income (Table 7)
## =============================================================================

rm(list = ls())
cat("=== Real-data variance examples: BRFSS BMI + CPS CPS income ===\n\n")

## -----------------------------------------------------------------------------
## Setup
## -----------------------------------------------------------------------------
library(dplyr)

source("R/spdv_functions.R")

## =============================================================================
## 1. BRFSS BMI Analysis (Louisiana 2022)
## =============================================================================
cat("------------------------------------------------------------\n")
cat("1. BRFSS BMI: Louisiana 2022\n")
cat("------------------------------------------------------------\n")

brfss <- readRDS("data/brfss_bmi_subset.rds")
if (is.data.frame(brfss)) {
  bmi_col <- grep("bmi|BMI", names(brfss), ignore.case = TRUE, value = TRUE)
  bmi <- brfss[[bmi_col[1]]]
} else {
  bmi <- as.numeric(brfss)
}
bmi <- bmi[is.finite(bmi) & bmi > 0]
n_bmi <- length(bmi)

cat(sprintf("BRFSS BMI: n = %d observations\n", n_bmi))

## Point estimates
s2_bmi    <- compute_spdv(bmi)
kappa_bmi <- kurtosis_hat(bmi)
se_bmi    <- se_spdv(bmi)

cat(sprintf("Point estimates: s^2 = %.2f, SE = %.4f, kappa_hat = %.2f\n",
            s2_bmi, se_bmi, kappa_bmi))

## METHOD 1: Chi-square
set.seed(20220301)
chi_ci_bmi     <- chi_square_ci(bmi)
chi_length_bmi <- chi_ci_bmi[2] - chi_ci_bmi[1]
cat(sprintf("Chi-square 95%% CI: [%.2f, %.2f] (length = %.2f)\n",
            chi_ci_bmi[1], chi_ci_bmi[2], chi_length_bmi))

## METHOD 2: Normal approximation
normal_ci_bmi     <- normal_ci(bmi)
normal_length_bmi <- normal_ci_bmi[2] - normal_ci_bmi[1]
cat(sprintf("Normal 95%% CI: [%.2f, %.2f] (length = %.2f)\n",
            normal_ci_bmi[1], normal_ci_bmi[2], normal_length_bmi))

## METHOD 3: Percentile bootstrap
pct_ci_bmi     <- percentile_ci(bmi, B = 5000)
pct_length_bmi <- pct_ci_bmi[2] - pct_ci_bmi[1]
cat(sprintf("Percentile Boot 95%% CI: [%.2f, %.2f] (length = %.2f)\n",
            pct_ci_bmi[1], pct_ci_bmi[2], pct_length_bmi))

## METHOD 4: Studentized bootstrap (SPDV)
stud_res_bmi    <- studentized_bootstrap(bmi, B = 5000, parallel = FALSE)
stud_ci_bmi     <- stud_res_bmi$ci_95
stud_length_bmi <- stud_ci_bmi[2] - stud_ci_bmi[1]
stud_runtime_bmi <- round(stud_res_bmi$runtime_s, 1)
cat(sprintf("Studentized Boot 95%% CI: [%.2f, %.2f] (length = %.2f)\n",
            stud_ci_bmi[1], stud_ci_bmi[2], stud_length_bmi))
cat(sprintf("Studentized runtime: %.1fs\n", stud_runtime_bmi))

cat("\n", rep("=", 80), "\n\n")


## =============================================================================
## 2. CPS 2022 Household Income Analysis
## =============================================================================
cat("------------------------------------------------------------\n")
cat("2. Economic Inequality: CPS 2022 Household Income\n")
cat("------------------------------------------------------------\n")

cps_path <- "data/cps_income_subset.rds"
if (!file.exists(cps_path)) {
  stop("File not found: ", cps_path, "\nSee data/README-data.md for instructions.")
}

cps <- readRDS(cps_path)
## Expect data.frame with 'income' column (USD)
if (is.data.frame(cps)) {
  if (!"income" %in% names(cps))
    stop("cps_income_subset.rds must contain column 'income'")
  inc <- cps$income
} else {
  inc <- as.numeric(cps)
}

## Positive finite income
inc   <- inc[is.finite(inc) & inc > 0]
n_inc <- length(inc)
cat(sprintf("Sample size (processed) CPS income: n = %d\n", n_inc))

## Point estimates
s2_inc    <- compute_spdv(inc)
kappa_inc <- kurtosis_hat(inc)
se_inc    <- se_spdv(inc)
cat(sprintf("Point estimates: s^2 = %.2e, SE = %.2e, kappa_hat = %.2f\n",
            s2_inc, se_inc, kappa_inc))

## METHOD 1: Chi-square
set.seed(20220301)
chi_ci_inc     <- chi_square_ci(inc)
chi_length_inc <- chi_ci_inc[2] - chi_ci_inc[1]
cat(sprintf("Chi-square 95%% CI: [%.2e, %.2e] (length = %.2e)\n",
            chi_ci_inc[1], chi_ci_inc[2], chi_length_inc))

## METHOD 2: Normal approximation
normal_ci_inc     <- normal_ci(inc)
normal_length_inc <- normal_ci_inc[2] - normal_ci_inc[1]
cat(sprintf("Normal 95%% CI: [%.2e, %.2e] (length = %.2e)\n",
            normal_ci_inc[1], normal_ci_inc[2], normal_length_inc))

## METHOD 3: Percentile bootstrap
pct_ci_inc     <- percentile_ci(inc, B = 5000)
pct_length_inc <- pct_ci_inc[2] - pct_ci_inc[1]
cat(sprintf("Percentile Boot 95%% CI: [%.2e, %.2e] (length = %.2e)\n",
            pct_ci_inc[1], pct_ci_inc[2], pct_length_inc))

## METHOD 4: Studentized bootstrap (SPDV)
stud_res_inc     <- studentized_bootstrap(inc, B = 5000, parallel = FALSE)
stud_ci_inc      <- stud_res_inc$ci_95
stud_length_inc  <- stud_ci_inc[2] - stud_ci_inc[1]
stud_runtime_inc <- round(stud_res_inc$runtime_s, 1)
cat(sprintf("Studentized Boot 95%% CI: [%.2e, %.2e] (length = %.2e)\n",
            stud_ci_inc[1], stud_ci_inc[2], stud_length_inc))
cat(sprintf("Studentized runtime: %.1fs\n", stud_runtime_inc))

cat(rep("=", 60), "\n\n")

## Save CPS numeric results (Table 7 / supplement)
cps_results <- data.frame(
  Method = c("Chi-square", "Normal", "Pctl Boot", "Stud Boot"),
  s2     = rep(s2_inc, 4),
  CI     = c(sprintf("[%.2e, %.2e]", chi_ci_inc[1],    chi_ci_inc[2]),
             sprintf("[%.2e, %.2e]", normal_ci_inc[1], normal_ci_inc[2]),
             sprintf("[%.2e, %.2e]", pct_ci_inc[1],    pct_ci_inc[2]),
             sprintf("[%.2e, %.2e]", stud_ci_inc[1],   stud_ci_inc[2])),
  Length_1e9 = c(chi_length_inc,
                 normal_length_inc,
                 pct_length_inc,
                 stud_length_inc) / 1e9
)
dir.create("results", showWarnings = FALSE)
write.csv(cps_results, "results/table7_cps_income.csv", row.names = FALSE)
cat("✓ CPS results saved: results/table7_cps_income.csv\n\n")


## =============================================================================
## 3. LaTeX Table 7 (BRFSS BMI + CPS Income)
## =============================================================================
cat("✓ LaTeX Table 7 (Real Data, All Methods):\n\n")

cat("\\begin{table}[htbp]\n")
cat("\\centering\n")
cat("\\caption{Real-data variance inference comparison}\n")
cat("\\label{tab:real-data}\n")
cat("\\setlength{\\tabcolsep}{4pt}\n")
cat("\\resizebox{\\textwidth}{!}{%\n")
cat("\\begin{threeparttable}\n")
cat("\\begin{tabular}{lccccccc}\n")
cat("\\toprule\n")
cat(sprintf("& \\multicolumn{3}{c}{BRFSS BMI ($n=%d$)} & ",
            n_bmi))
cat(sprintf("\\multicolumn{3}{c}{CPS Income ($n=%d$)} \\\\\n", n_inc))
cat("\\cmidrule(lr){2-4} \\cmidrule(lr){5-7}\n")
cat("Method & $s^2$ & 95\\% CI & Length & $s^2$ & 95\\% CI ($\\times 10^9$) & Length ($\\times 10^9$) \\\\\n")
cat("\\midrule\n")

## Chi-square row
cat(sprintf("Chi-square & %.2f & [%.2f, %.2f] & %.2f & %.2f & [%.2f, %.2f] & %.3f \\\\\n",
            s2_bmi,
            round(chi_ci_bmi[1], 2), round(chi_ci_bmi[2], 2),
            round(chi_length_bmi, 2),
            s2_inc / 1e9,
            chi_ci_inc[1] / 1e9, chi_ci_inc[2] / 1e9,
            chi_length_inc / 1e9))

## Normal row
cat(sprintf("Normal & %.2f & [%.2f, %.2f] & %.2f & %.2f & [%.2f, %.2f] & %.3f \\\\\n",
            s2_bmi,
            round(normal_ci_bmi[1], 2), round(normal_ci_bmi[2], 2),
            round(normal_length_bmi, 2),
            s2_inc / 1e9,
            normal_ci_inc[1] / 1e9, normal_ci_inc[2] / 1e9,
            normal_length_inc / 1e9))

## Percentile bootstrap row
cat(sprintf("Pctl Boot & %.2f & [%.2f, %.2f] & %.2f & %.2f & [%.2f, %.2f] & %.3f \\\\\n",
            s2_bmi,
            round(pct_ci_bmi[1], 2), round(pct_ci_bmi[2], 2),
            round(pct_length_bmi, 2),
            s2_inc / 1e9,
            pct_ci_inc[1] / 1e9, pct_ci_inc[2] / 1e9,
            pct_length_inc / 1e9))

## Studentized bootstrap row (bold)
cat(sprintf("\\textbf{Stud Boot} & \\textbf{%.2f} & \\textbf{[%.2f, %.2f]} & \\textbf{%.2f} & ",
            s2_bmi,
            round(stud_ci_bmi[1], 2), round(stud_ci_bmi[2], 2),
            round(stud_length_bmi, 2)))
cat(sprintf("\\textbf{%.2f} & \\textbf{[%.2f, %.2f]} & \\textbf{%.3f} \\\\\n",
            s2_inc / 1e9,
            stud_ci_inc[1] / 1e9, stud_ci_inc[2] / 1e9,
            stud_length_inc / 1e9))

cat("\\midrule\n")
cat(sprintf("Kurtosis $\\hat{\\kappa}$ & %.2f & --- & --- & %.2f & --- & --- \\\\\n",
            kappa_bmi, kappa_inc))
cat(sprintf("Runtime ($B=5$k) & %.1fs & --- & --- & %.1fs & --- & --- \\\\\n",
            stud_runtime_bmi, stud_runtime_inc))
cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\begin{tablenotes}\n")
cat("\\small\n")
cat("\\item Note: $s^2$ is the sample variance (SPDV). Length = upper $-$ lower bounds.\n")
cat("\\item Pctl Boot: Percentile bootstrap; Stud Boot: Studentized bootstrap (SPDV).\n")
cat("\\item Data: BRFSS 2022 Louisiana ($n=")
cat(n_bmi)
cat("$), CPS 2022 ($n=")
cat(n_inc)
cat("$).\n")
cat("\\end{tablenotes}\n")
cat("\\end{threeparttable}}\n")
cat("\\end{table}\n\n")

## Optional: save BMI-only summary to CSV (for supplement)
bmi_results <- data.frame(
  Method     = c("Chi-square", "Normal", "Pctl Boot", "Stud Boot"),
  s2         = rep(s2_bmi, 4),
  CI         = c(sprintf("[%.2f, %.2f]", chi_ci_bmi[1],     chi_ci_bmi[2]),
                 sprintf("[%.2f, %.2f]", normal_ci_bmi[1],  normal_ci_bmi[2]),
                 sprintf("[%.2f, %.2f]", pct_ci_bmi[1],     pct_ci_bmi[2]),
                 sprintf("[%.2f, %.2f]", stud_ci_bmi[1],    stud_ci_bmi[2])),
  Length     = c(chi_length_bmi, normal_length_bmi,
                 pct_length_bmi, stud_length_bmi),
  stringsAsFactors = FALSE
)
write.csv(bmi_results, "results/table7_brfss_bmi.csv", row.names = FALSE)
cat("✓ BRFSS BMI results saved: results/table7_brfss_bmi.csv\n")

cat("\n🎉 COMPLETE: BRFSS BMI + CPS income analysis and LaTeX Table 7 generated.\n")
