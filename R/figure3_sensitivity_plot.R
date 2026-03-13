# =============================================================================
# Generate Figure \ref{fig:sensitivity-B}: Coverage vs. B
# Manuscript: Scalable Studentized Bootstrap Variance Inference...
# Saves: figures/sensitivity_B_plot.pdf
# =============================================================================

library(ggplot2)
library(dplyr)

# ---------------------------------------------------------------------------
# 1. Sensitivity data (exact match to Table \ref{tab:sensitivity-B})
# ---------------------------------------------------------------------------
sensitivity_data <- data.frame(
  dist = rep(c("ln", "t3"), each = 5),
  B = rep(c(500, 1000, 2000, 5000, 10000), times = 2),
  # stud = c(78.4, 81.2, 82.6, 83.3, 83.7,    # Lognormal n=50
  #         79.8, 81.9, 82.8, 83.6, 83.9)     # t_3 n=50
    stud = c(81.5, 82.5, 81.7, 81.1, 82.3,    # Lognormal n=50
           83, 83.4, 83.2, 82.9, 84.0)     # t_3 n=50

) %>%
  mutate(
    dist = factor(dist, levels = c("ln", "t3"),
                  labels = c("Lognormal", expression(t[3]))),
    Coverage = stud  # For cleaner plotting
  )

# ---------------------------------------------------------------------------
# 2. Create publication-quality plot
# ---------------------------------------------------------------------------
p_sensitivity <- ggplot(sensitivity_data, 
                        aes(x = B, y = Coverage, 
                            color = dist, shape = dist, linetype = dist)) +
  geom_hline(yintercept = 95, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.5) +
  scale_x_continuous(breaks = c(500, 1000, 2000, 5000, 10000),
                     labels = scales::comma_format()) +
  scale_y_continuous(limits = c(75, 90), breaks = seq(75, 90, by = 5),
                     name = "Empirical coverage (%)") +
  scale_color_manual(values = c("Lognormal" = "#E41A1C", "t[3]" = "#377EB8")) +
  scale_linetype_manual(values = c("Lognormal" = "solid", "t[3]" = "longdash")) +
  scale_shape_manual(values = c("Lognormal" = 16, "t[3]" = 17)) +
  labs(
    x = "Number of bootstrap replications B",
    color = "Distribution", shape = "Distribution", linetype = "Distribution",
    title = "Coverage stabilization vs. B (n = 50)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

# ---------------------------------------------------------------------------
# 3. Save both color and grayscale versions
# ---------------------------------------------------------------------------
if (!dir.exists("figures")) dir.create("figures")

ggsave("figures/sensitivity_B_plot2.pdf", p_sensitivity, 
       width = 9, height = 6, dpi = 300, device = "pdf")

# Grayscale for manuscript (JSCS preference)
p_grayscale <- p_sensitivity + 
  scale_color_grey(start = 0.2, end = 0.7) +
  guides(color = guide_legend(override.aes = list(shape = c(16, 17))))

ggsave("figures/sensitivity_B_plot_grayscale2.pdf", p_grayscale, 
       width = 9, height = 6, dpi = 300, device = "pdf")

cat("✓ Figures saved:\n")
cat("  figures/sensitivity_B_plot2.pdf (color)\n")
cat("  figures/sensitivity_B_plot_grayscale.pdf\n")
