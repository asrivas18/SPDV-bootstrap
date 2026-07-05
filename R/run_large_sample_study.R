# =============================================================================
# run_large_sample_study.R
# Large-sample behavior and scalability for SPDV studentized bootstrap
# n in {100, 500, 2000, 10,000, 100,000}, designs: Normal, t5, Lognormal
# Outputs: results/large_sample_variance_summary.csv
# =============================================================================

rm(list = ls())
set.seed(987654)

## You need spdv_functions.R with:
## - chi_square_ci(x)
## - normal_ci(x)
## - percentile_ci(x, B)
## - studentized_bootstrap(x, B, parallel = FALSE)
source("R/spdv_functions.R")

dir.create("results", showWarnings = FALSE)

true_var <- 1
alpha <- 0.05

## Large-n grid
n_vec    <- c(100, 500, 2000, 10000, 100000)
designs  <- c("normal", "t5", "lnorm")

## Monte Carlo and bootstrap settings
## Start modest; increase if runtime is acceptable
M_large  <- 500    # MC replications per (n, design); you can raise 200 to 500
B_large  <- 5000   # bootstrap resamples per replication; you can raise 2000 to 5000

cat("Large-sample study:\n")
cat(sprintf("  n values: %s\n", paste(n_vec, collapse = ", ")))
cat(sprintf("  designs:  %s\n", paste(designs, collapse = ", ")))
cat(sprintf("  M = %d, B = %d\n\n", M_large, B_large))

## -----------------------------------------------------------------------------
## 1. Unit-variance data generator (consistent with your main script)
## -----------------------------------------------------------------------------
generate_data <- function(n, dist) {
  if (dist == "normal") {
    rnorm(n)                        # N(0,1), var=1
  } else if (dist == "t5") {
    rt(n, df = 5) / sqrt(5/3)       # Var = 1
  } else if (dist == "lnorm") {
    raw <- rlnorm(n, meanlog = 0, sdlog = 1)
    true_sd <- sqrt(exp(1) * (exp(1) - 1))
    (raw - exp(0.5)) / true_sd      # centered, approx var=1
  } else {
    stop("Unknown dist: ", dist)
  }
}

## -----------------------------------------------------------------------------
## 2. Main loop over (n, design)
## -----------------------------------------------------------------------------
large_results <- data.frame()

for (design_name in designs) {
  cat("Design:", design_name, "\n")
  for (n in n_vec) {
    cat(sprintf("  n = %d\n", n))

    cov_counts <- c(chi_sq = 0, normal = 0, pctl = 0, stud = 0)
    width_sums <- c(chi_sq = 0, normal = 0, pctl = 0, stud = 0)
    rt_sums    <- c(chi_sq = 0, normal = 0, pctl = 0, stud = 0)

    t_start <- proc.time()[3]

    for (m in seq_len(M_large)) {
      x <- generate_data(n, design_name)

      ## CHI-SQUARE
      t0 <- proc.time()[3]
      chi_ci <- chi_square_ci(x)
      rt_sums["chi_sq"] <- rt_sums["chi_sq"] + (proc.time()[3] - t0)
      cov_counts["chi_sq"] <- cov_counts["chi_sq"] +
        as.integer(chi_ci[1] <= true_var && chi_ci[2] >= true_var)
      width_sums["chi_sq"] <- width_sums["chi_sq"] + diff(chi_ci)

      ## NORMAL
      t0 <- proc.time()[3]
      norm_ci <- normal_ci(x)
      rt_sums["normal"] <- rt_sums["normal"] + (proc.time()[3] - t0)
      cov_counts["normal"] <- cov_counts["normal"] +
        as.integer(norm_ci[1] <= true_var && norm_ci[2] >= true_var)
      width_sums["normal"] <- width_sums["normal"] + diff(norm_ci)

      ## PERCENTILE BOOTSTRAP
      t0 <- proc.time()[3]
      pct_ci <- percentile_ci(x, B = B_large)
      rt_sums["pctl"] <- rt_sums["pctl"] + (proc.time()[3] - t0)
      cov_counts["pctl"] <- cov_counts["pctl"] +
        as.integer(pct_ci[1] <= true_var && pct_ci[2] >= true_var)
      width_sums["pctl"] <- width_sums["pctl"] + diff(pct_ci)

      ## STUDENTIZED BOOTSTRAP
      t0 <- proc.time()[3]
      stud_res <- studentized_bootstrap(x, B = B_large, parallel = TRUE)
      rt_sums["stud"] <- rt_sums["stud"] + (proc.time()[3] - t0)
      stud_ci <- stud_res$ci_95
      cov_counts["stud"] <- cov_counts["stud"] +
        as.integer(stud_ci[1] <= true_var && stud_ci[2] >= true_var)
      width_sums["stud"] <- width_sums["stud"] + diff(stud_ci)
    }

    t_total <- proc.time()[3] - t_start

    cov <- cov_counts / M_large
    width <- width_sums / M_large
    rt_avg <- rt_sums / M_large

    large_results <- rbind(
      large_results,
      data.frame(
        design = design_name,
        n = n,
        M = M_large,
        B = B_large,
        cov_chi_sq = cov["chi_sq"],
        cov_normal = cov["normal"],
        cov_pctl   = cov["pctl"],
        cov_stud   = cov["stud"],
        width_chi_sq = width["chi_sq"],
        width_normal = width["normal"],
        width_pctl   = width["pctl"],
        width_stud   = width["stud"],
        rt_chi_sq = rt_avg["chi_sq"],
        rt_normal = rt_avg["normal"],
        rt_pctl   = rt_avg["pctl"],
        rt_stud   = rt_avg["stud"],
        rt_total  = t_total
      )
    )

    cat(sprintf("    coverage (studentized): %.3f\n", cov["stud"]))
    cat(sprintf("    avg width (studentized): %.4f\n", width["stud"]))
    cat(sprintf("    total runtime for (n=%d, %s): %.2f s\n\n", n, design_name, t_total))
  }
}

## -----------------------------------------------------------------------------
## 3. Save summary for tables/plots
## -----------------------------------------------------------------------------
outfile <- "results/large_sample_variance_summary.csv"
write.csv(large_results, outfile, row.names = FALSE)
cat("✓ Large-sample summary saved to", outfile, "\n")

## Simple preview
print(head(large_results))

# =============================================================================
# FIGURE 3: Large-sample coverage line plot
# Reads data directly from large_sample_variance_summary.csv
# =============================================================================

rm(list = ls())
library(ggplot2)
library(dplyr)
library(tidyr)

# Create target directory if it does not exist
dir.create("figures", showWarnings = FALSE)

# Input/output paths
summary_file <- "large_sample_variance_summary.csv"
if (!file.exists(summary_file)) {
  summary_file <- "results/large_sample_variance_summary.csv"
}

# -------------------------------------------------------------------------
# Read summary file and reshape coverage data
# -------------------------------------------------------------------------
large_df <- read.csv(summary_file, stringsAsFactors = FALSE)

plot_data <- large_df %>%
  tidyr::pivot_longer(
    cols      = c(cov_chi_sq, cov_normal, cov_pctl, cov_stud),
    names_to  = "method",
    values_to = "coverage"
  ) %>%
  dplyr::mutate(
    method = factor(
      method,
      levels = c("cov_chi_sq", "cov_normal", "cov_pctl", "cov_stud"),
      labels = c("Chi-square", "Normal approx.", "Percentile", "Studentized")
    ),
    design = factor(
      design,
      levels = c("normal", "t5", "lnorm"),
      labels = c("Normal", "t[5]", "Lognormal")
    )
  )

# -------------------------------------------------------------------------
# Base plot (JSCS-style grayscale)
# -------------------------------------------------------------------------
p_large <- ggplot(
  plot_data,
  aes(x = n, y = coverage * 100,
      linetype = method, shape = method, color = method)
) +
  geom_hline(yintercept = 95, linetype = "dashed",
             linewidth = 0.7, color = "black") +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.8) +
  facet_wrap(~ design, ncol = 3, labeller = label_parsed) +
  scale_color_grey(start = 0.2, end = 0.7) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash", "longdash")) +
  scale_x_log10(
    breaks = c(100, 500, 2000, 10000, 100000),
    labels = c("100", "500", "2,000", "10,000", "100,000")
  ) +
  labs(
    x = "Sample size (n) on log scale",
    y = "Empirical coverage (%)",
    linetype = "Method",
    shape    = "Method",
    color    = "Method"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position  = "bottom",
    legend.title     = element_text(size = 11),
    legend.text      = element_text(size = 10),
    strip.text       = element_text(size = 12),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

# -------------------------------------------------------------------------
# Save grayscale and color versions
# -------------------------------------------------------------------------
ggsave("figures/Figure3_LargeSampleCoverage_JSCS.pdf",
       p_large, width = 14, height = 4.8, dpi = 300)

ggsave("figures/Figure3_LargeSampleCoverage_color.pdf",
       p_large + scale_color_brewer(palette = "Set1"),
       width = 14, height = 4.8, dpi = 300)

cat("✓ Figure 3 saved (grayscale and color versions)\n")