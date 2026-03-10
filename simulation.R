# ============================================================================
# Reproduce Table 3 and Figure 2: Coverage Study for SPDV Studentized Bootstrap
# Manuscript: Scalable Studentized Bootstrap Variance Inference via Linear-Time Pairwise Variance Representation
# Journal: Journal of Statistical Computation and Simulation (JSCS)
# Submission-freeze version (true bootstrap-t implementation)
# Runtime: ~60–120 minutes depending on hardware (parallel)
# ============================================================================

rm(list = ls())
set.seed(12345)

# ---------------------------------------------------------------------------
# Simulation size (change quick_test to TRUE for fast debugging)
# ---------------------------------------------------------------------------
quick_test <- FALSE
if (quick_test) {
  M <- 500     # Monte Carlo replications
  B <- 500     # Bootstrap resamples per replication
} else {
  M <- 5000    # Full run
  B <- 5000
}

alpha <- 0.05
true_var <- 1  # Population variance σ² = 1

# ---------------------------------------------------------------------------
# Libraries
# ---------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(foreach)
library(doParallel)
library(doRNG)
library(kableExtra)

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
n_vals <- c(20, 50, 100)
dists <- c("normal", "t5", "t3", "ln")  # ln = lognormal
sim_grid <- expand.grid(n = n_vals, dist = dists, stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# Parallel setup
# ---------------------------------------------------------------------------
n_cores <- max(1, parallel::detectCores() - 1)
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)
doRNG::registerDoRNG(12345)  # Reproducible RNG across parallel workers
cat("Using", n_cores, "cores\n")
cat("Monte Carlo replications:", M, "\n")
cat("Bootstrap resamples:", B, "\n\n")

# ============================================================================
# Coverage Simulation (parallel over grid rows)
# ============================================================================
coverage_results <- foreach(i = 1:nrow(sim_grid), .combine = rbind) %dopar% {
  n <- sim_grid$n[i]
  dist_name <- sim_grid$dist[i]
  counts <- c(chi_sq = 0, normal = 0, pctl = 0, stud = 0)
  skipped_stud <- 0  # Track skipped unstable SE cases

  for (m in 1:M) {
    # -----------------------------------------------------------------------
    # Generate data (all scaled to variance = 1)
    # -----------------------------------------------------------------------
    x <- switch(dist_name,
      normal = rnorm(n),
      t5     = rt(n, 5) / sqrt(5/3),
      t3     = rt(n, 3) / sqrt(3),
      ln     = (exp(rnorm(n)) - exp(0.5)) / sqrt(exp(1) * (exp(1) - 1))
    )
    s2 <- var(x)
    xbar <- mean(x)

    # -----------------------------------------------------------------------
    # 1. Chi-square interval (exact under normality)
    # -----------------------------------------------------------------------
    chi_low  <- (n-1) * s2 / qchisq(1 - alpha/2, n-1)
    chi_high <- (n-1) * s2 / qchisq(alpha/2, n-1)
    if (chi_low <= true_var && chi_high >= true_var)
      counts["chi_sq"] <- counts["chi_sq"] + 1

    # -----------------------------------------------------------------------
    # 2. Normal approximation (fourth-moment plug-in SE)
    # -----------------------------------------------------------------------
    mu4_hat <- mean((x - xbar)^4)
    se_normal <- sqrt(pmax((mu4_hat - s2^2) / n, 0))
    normal_low  <- s2 - qnorm(1 - alpha/2) * se_normal
    normal_high <- s2 + qnorm(1 - alpha/2) * se_normal
    if (normal_low <= true_var && normal_high >= true_var)
      counts["normal"] <- counts["normal"] + 1

    # -----------------------------------------------------------------------
    # Bootstrap replicates
    # -----------------------------------------------------------------------
    boot_vals <- numeric(B)
    boot_se   <- numeric(B)

    for (b in 1:B) {
      xb <- sample(x, n, replace = TRUE)
      s2_b <- var(xb)
      boot_vals[b] <- s2_b

      # Resample-specific influence values (unscaled psi)
      xbar_b <- mean(xb)
      psi_b  <- (xb - xbar_b)^2 - s2_b
      boot_se[b] <- sqrt(var(psi_b) / n)
    }

    # -----------------------------------------------------------------------
    # 3. Percentile interval
    # -----------------------------------------------------------------------
    pctl_low  <- quantile(boot_vals, alpha/2, na.rm = TRUE)
    pctl_high <- quantile(boot_vals, 1 - alpha/2, na.rm = TRUE)
    if (pctl_low <= true_var && pctl_high >= true_var)
      counts["pctl"] <- counts["pctl"] + 1

    # -----------------------------------------------------------------------
    # 4. True Bootstrap-t (resample-specific SE)
    # -----------------------------------------------------------------------
    valid <- boot_se > 0 & is.finite(boot_se)
    skipped_stud <- skipped_stud + sum(!valid)

    if (any(valid)) {
      t_star <- (boot_vals[valid] - s2) / boot_se[valid]
      t_low  <- quantile(t_star, alpha/2, na.rm = TRUE)
      t_high <- quantile(t_star, 1 - alpha/2, na.rm = TRUE)

      # Original sample SE (unscaled psi)
      psi_hat <- (x - xbar)^2 - s2
      se_orig <- sqrt(var(psi_hat) / n)

      stud_low  <- s2 - t_high * se_orig
      stud_high <- s2 - t_low  * se_orig

      if (stud_low <= true_var && stud_high >= true_var)
        counts["stud"] <- counts["stud"] + 1
    }
  }  # end M loop

  coverage <- counts / M
  mc_se    <- sqrt(coverage * (1 - coverage) / M)

  data.frame(
    n     = n,
    dist  = dist_name,
    chi_sq = coverage["chi_sq"],
    normal = coverage["normal"],
    pctl   = coverage["pctl"],
    stud   = coverage["stud"],
    mc_se_chi_sq = mc_se["chi_sq"],
    mc_se_normal = mc_se["normal"],
    mc_se_pctl   = mc_se["pctl"],
    mc_se_stud   = mc_se["stud"],
    skipped_stud = skipped_stud / M   # fraction skipped per replication
  )
}  # end foreach

stopCluster(cl)

# ============================================================================
# Save raw results (for reproducibility)
# ============================================================================
saveRDS(coverage_results, "coverage_results.rds")
write.csv(coverage_results, "coverage_results.csv", row.names = FALSE)

# ============================================================================
# Generate Table 2 (formatted output)
# ============================================================================
coverage_results$dist <- factor(
  coverage_results$dist,
  levels = c("normal", "t5", "t3", "ln"),
  labels = c("Normal", "t[5]", "t[3]", "Lognormal")
)

table2 <- coverage_results %>%
  select(n, dist, chi_sq, normal, pctl, stud) %>%
  mutate(across(c(chi_sq:stud), ~round(. * 100, 1))) %>%
  arrange(dist, n)

cat("\nTable 2: Empirical Coverage (%)\n")
print(table2)

# LaTeX table for manuscript (copy-paste into paper)
kable(table2,
      format = "latex",
      booktabs = TRUE,
      escape = FALSE,
      caption = "Empirical coverage probabilities (%) for $\\sigma^2=1$.") %>%
      kable_styling(latex_options = "hold_position") %>%
  pack_rows("Normal", 1, 3) %>%
  pack_rows("$t_5$", 4, 6) %>%
  pack_rows("$t_3$", 7, 9) %>%
  pack_rows("Lognormal", 10, 12) %>%
  column_spec(6, bold = TRUE)

# ============================================================================
# Generate Figure 2 (JSCS style)
# ============================================================================
plot_data <- coverage_results %>%
    pivot_longer(cols = chi_sq:stud, names_to = "method", values_to = "coverage") %>%
    mutate(Method = factor(Method, 
                           levels = c("chi_sq", "normal", "pctl", "stud"),
                           labels = c("Chi-square", "Normal approx.", "Percentile", "Studentized")))


p_JSCS <- ggplot(plot_data,
                 aes(x = n,
                     y = coverage * 100,
                     linetype = method,
                     shape = method,
                     color = method)) +
  geom_hline(yintercept = 95,
             linetype = "dashed",
             linewidth = 0.7,
             color = "black") +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.8) +
  facet_wrap(~ dist, ncol = 2, labeller = label_parsed) +
  scale_color_grey(start = 0.2, end = 0.7) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash", "longdash")) +
  scale_x_continuous(breaks = c(20, 50, 100)) +
  labs(x = "Sample size n",
       y = "Empirical coverage (%)",
       linetype = "Method",
       shape = "Method") +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

ggsave("Figure2_Coverage_JSCS.pdf",
       p_JSCS,
       width = 12,
       height = 10,
       dpi = 300)

# Color version (for presentation or supplement)
ggsave("Figure2_Coverage_color.pdf",
       p_JSCS + scale_color_brewer(palette = "Set1"),
       width = 12,
       height = 10,
       dpi = 300)

cat("\nSimulation complete.\n")
cat("Files saved: coverage_results.rds, coverage_results.csv, Figure2_Coverage_JSCS.pdf, Figure2_Coverage_color.pdf\n")