# Figure 9 -- mean-variance plot: one dot per patient, 24 monthly counts each
source("R/00_setup.R")
cohort <- read.csv("data/voc_violations_sim.csv")

labels <- c(poisson = "Poisson", negbin = "Negative binomial", zip = "Zero-inflated Poisson")

per_patient <- cohort |>
  pivot_longer(c(poisson, negbin, zip), names_to = "model", values_to = "count") |>
  group_by(model, id) |>
  summarise(mean = mean(count), variance = var(count), .groups = "drop") |>
  mutate(model = factor(labels[model], levels = labels))

p <- ggplot(per_patient, aes(mean, variance, colour = model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = col_muted) +
  geom_point(alpha = 0.45, size = 1.6) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.8) +
  scale_colour_manual(values = c("Poisson" = col_placebo,
                                 "Negative binomial" = col_hu,
                                 "Zero-inflated Poisson" = col_third)) +
  coord_cartesian(xlim = c(0, 1.6), ylim = c(0, 4)) +
  labs(x = "Individual mean (VOCs per month)",
       y = "Individual variance",
       subtitle = "Dashed line: variance = mean. Solid lines: linear regression per model")

save_fig(p, "fig09_mean_variance.png", width = 7, height = 4.8)

# Index of dispersion, combined over patients: sum of (T - 1) * s^2 / mean ~ chi-square
dispersion_test <- per_patient |>
  filter(mean > 0) |>
  group_by(model) |>
  summarise(statistic = sum(23 * variance / mean), df = 23 * n(),
            ratio = statistic / df,
            p_value = pchisq(statistic, df, lower.tail = FALSE))
print(dispersion_test)

# The pooled trap: ignoring patients makes even Poisson data look overdispersed
cohort |> summarise(pooled_mean = mean(poisson), pooled_variance = var(poisson))
