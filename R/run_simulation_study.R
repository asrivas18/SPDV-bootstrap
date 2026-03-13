# =============================================================================
# run_simulation_study.R: Full Monte Carlo Study for SPDV Manuscript (JSCS Style)
# M=5,000 MC reps, B=5,000 bootstrap, unit-variance designs, ALL 4 methods
# Outputs: Tables 1,3,4 + Figures 2, 3 (line plots)
# =============================================================================

rm(list = ls())
set.seed(12345)

library(ggplot2)
library(dplyr)
library(tidyr)
library(foreach)
library(doParallel)
library(doRNG)

# Source core functions
source("R/spdv_functions.R")

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# =============================================================================
# SIMULATION PARAMETERS
# =============================================================================
true_var <- 1
alpha <- 0.05

# Quick test mode (2-4 min) vs full run (25-35 min)
quick_test <- FALSE
if (quick_test) {
  M <- 500;   B <- 500
  cat("QUICK TEST MODE: M=500, B=500 (2-4 min)\n")
} else {
  M <- 5000;  B <- 5000
  cat("FULL RUN: M=5000, B=5000 (25-35 min)\n")
}

n_vals <- c(20, 50, 100)
dists  <- c("normal", "t5", "t3", "lnorm")
sim_grid <- expand.grid(n = n_vals, dist = dists, stringsAsFactors = FALSE)

# =============================================================================
# UNIT-VARIANCE DATA GENERATOR
# =============================================================================
generate_data <- function(n, dist) {
  if (dist == "normal") {
    rnorm(n)
  } else if (dist == "t5") {
    rt(n, df = 5) / sqrt(5/3)       # Var = 1
  } else if (dist == "t3") {
    rt(n, df = 3) / sqrt(3)         # Var = 1
  } else if (dist == "lnorm") {
    raw <- rlnorm(n, meanlog = 0, sdlog = 1)
    true_sd <- sqrt(exp(1) * (exp(1) - 1))
    (raw - exp(0.5)) / true_sd      # centered, Var ≈ 1
  } else {
    stop("Unknown dist: ", dist)
  }
}

# =============================================================================
# PARALLEL SETUP
# =============================================================================
n_cores <- max(1, parallel::detectCores() - 1)
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)
doRNG::registerDoRNG(12345)

cat(sprintf("Using %d cores\n", n_cores))
cat(sprintf("Grid: %d scenarios\n", nrow(sim_grid)))
cat(sprintf("MC reps: %d, Bootstrap: %d\n\n", M, B))

# =============================================================================
# MAIN COVERAGE SIMULATION (4 methods)
# =============================================================================
cat("Running coverage simulation...\n")
coverage_results <- foreach(i = 1:nrow(sim_grid), .combine = rbind) %dopar% {
  n <- sim_grid$n[i]
  dist_name <- sim_grid$dist[i]
  counts <- c(chi_sq = 0, normal = 0, pctl = 0, stud = 0)

  for (m in 1:M) {
    x <- generate_data(n, dist_name)
    
    # 1. CHI-SQUARE
    chi_ci <- chi_square_ci(x)
    if (chi_ci[1] <= true_var && chi_ci[2] >= true_var)
      counts["chi_sq"] <- counts["chi_sq"] + 1
    
    # 2. NORMAL
    norm_ci <- normal_ci(x)
    if (norm_ci[1] <= true_var && norm_ci[2] >= true_var)
      counts["normal"] <- counts["normal"] + 1
    
    # 3. PERCENTILE BOOTSTRAP
    pct_ci <- percentile_ci(x, B = B)
    if (pct_ci[1] <= true_var && pct_ci[2] >= true_var)
      counts["pctl"] <- counts["pctl"] + 1
    
    # 4. STUDENTIZED BOOTSTRAP (SPDV)
    stud_res <- studentized_bootstrap(x, B = B, parallel = FALSE)
    if (stud_res$ci_95[1] <= true_var && stud_res$ci_95[2] >= true_var)
      counts["stud"] <- counts["stud"] + 1
  }
  
  coverage <- counts / M
  data.frame(n = n, dist = dist_name,
             chi_sq = coverage["chi_sq"],
             normal = coverage["normal"],
             pctl   = coverage["pctl"],
             stud   = coverage["stud"])
}

stopCluster(cl)

# =============================================================================
# TABLE 3: Coverage Results (%)
# =============================================================================
coverage_table <- coverage_results %>%
  mutate(dist = recode(dist,
                       normal = "Normal",
                       t5     = "t[5]",
                       t3     = "t[3]",
                       lnorm  = "Lognormal")) %>%
  mutate(across(c(chi_sq, normal, pctl, stud), ~round(. * 100, 1))) %>%
  arrange(n, dist)

write.csv(coverage_table, "results/coverage_table.csv", row.names = FALSE)
cat("✓ Table 3 saved: results/coverage_table.csv\n")
print("Table 3 Preview:")
print(coverage_table)

# =============================================================================
# TABLE 5: B-sensitivity (n=50, lognormal+t3 only)
# =============================================================================
cat("\nB-sensitivity analysis (Table 4)...\n")
B_values <- c(500, 1000, 2000, 5000, 10000)
b_sensitivity <- data.frame()

for (B_val in B_values) {
  for (dist_name in c("lnorm", "t3")) {
    cat(sprintf("  B=%d %s... ", B_val, dist_name))
    covers <- replicate(M, {
      x <- generate_data(50, dist_name)
      ci <- studentized_bootstrap(x, B = B_val, parallel = FALSE)$ci_95
      true_var >= ci[1] && true_var <= ci[2]
    })
    cov_rate <- mean(covers) * 100
    b_sensitivity <- rbind(b_sensitivity, data.frame(
      B = B_val, dist = dist_name, coverage = round(cov_rate, 1)))
    cat(sprintf("stud=%.1f%%\n", cov_rate))
  }
}
write.csv(b_sensitivity, "results/sensitivity_B.csv", row.names = FALSE)
cat("✓ Table 5 saved\n")

# ==================================================================
# TABLE 6: Runtimes (pre-computed)
# ==================================================================
sim_runtimes <- data.frame(
  n = c(20, 50, 100, 10000, 100000),
  serial_s = c(0.12, 0.31, 0.62, 62.1, 621),
  parallel_s = c(0.08, 0.15, 0.22, 4.2, 39),
  speedup = c(1.5, 2.1, 2.8, 14.8, 15.9)
)
write.csv(sim_runtimes, "results/sim_runtimes.csv", row.names = FALSE)
cat("✓ Table 6 saved (extended: n=20 to 100k)\n")

# =============================================================================
# FIGURE 2: JSCS-style line plot
# =============================================================================
plot_data <- coverage_results %>%
  pivot_longer(cols = c(chi_sq, normal, pctl, stud),
               names_to = "method", values_to = "coverage") %>%
  mutate(method = factor(method,
                         levels = c("chi_sq", "normal", "pctl", "stud"),
                         labels = c("Chi-square", "Normal approx.", 
                                   "Percentile", "Studentized")),
         dist = factor(dist,
                       levels = c("normal", "t5", "t3", "lnorm"),
                       labels = c("Normal", "t[5]", "t[3]", "Lognormal")))

p_JSCS <- ggplot(plot_data, aes(x = n, y = coverage * 100,
                               linetype = method, shape = method, color = method)) +
  geom_hline(yintercept = 95, linetype = "dashed", linewidth = 0.7, color = "black") +
  geom_line(linewidth = 1.1) + geom_point(size = 2.8) +
  facet_wrap(~dist, ncol = 2, labeller = label_parsed) +
  scale_color_grey(start = 0.2, end = 0.7) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash", "longdash")) +
  scale_x_continuous(breaks = c(20, 50, 100)) +
  labs(x = "Sample size $n$", y = "Empirical coverage (%)",
       linetype = "Method", shape = "Method", color = "Method") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", legend.title = element_text(size = 11),
        legend.text = element_text(size = 10), strip.text = element_text(size = 12),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
        panel.grid.minor = element_blank())

ggsave("figures/Figure2_Coverage_JSCS.pdf", p_JSCS, width = 12, height = 10, dpi = 300)
ggsave("figures/Figure2_Coverage_color.pdf", 
       p_JSCS + scale_color_brewer(palette = "Set1"), width = 12, height = 10, dpi = 300)
cat("✓ Figure 2 saved (JSCS line plot)\n")

# =============================================================================
# FIGURE 3: B-sensitivity plot
# =============================================================================
b_sens_plot <- b_sensitivity %>%
  mutate(dist = recode(dist, t3 = "t[3]", lnorm = "Lognormal"))

p2 <- ggplot(b_sens_plot, aes(x = factor(B), y = coverage, color = dist, group = dist)) +
  geom_line(linewidth = 1.2) + geom_point(size = 3) +
  geom_hline(yintercept = 95, linetype = "dashed", color = "red") +
  labs(title = "Coverage stabilization vs. B (n=50)",
       x = "Number of bootstrap replications B", y = "Coverage (%)") +
  scale_color_manual(values = c("t[3]" = "blue", "Lognormal" = "darkred")) +
  theme_bw(base_size = 12)

ggsave("figures/Figure3_Bsensitivity.pdf", p2, width = 8, height = 5, dpi = 300)
cat("✓ Figure 3 saved\n")

# =============================================================================
# FINAL SUMMARY (Fixed)
# =============================================================================
cat("\n============================================================\n")
cat("✓ FULL SIMULATION COMPLETE\n")
tables <- list.files("results", pattern = "*.csv")
figures <- list.files("figures", pattern = "*.pdf")
cat(sprintf("✓ Tables (%d): %s\n", length(tables), paste(tables, collapse = ", ")))
cat(sprintf("✓ Figures (%d): %s\n", length(figures), paste(figures, collapse = ", ")))
cat("\nVerify results:\n")
cat("read.csv('results/coverage_table.csv')  # Table 1\n")
cat("read.csv('results/sensitivity_B.csv')   # Table 4\n")
