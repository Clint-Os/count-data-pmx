# Figure 6 -- from one population value (theta) to individual mean counts (lambda_i)
source("R/00_setup.R")

theta <- 4.5 / 12
omega <- 0.6

people <- tibble(who = c("5th percentile", "Typical patient", "95th percentile"),
                 eta = qnorm(c(0.05, 0.50, 0.95), mean = 0, sd = omega)) |>
  mutate(lambda = theta * exp(eta),
         panel  = sprintf("%s\nlambda = %.3f", who, lambda),
         panel  = factor(panel, levels = panel))

dens <- tibble(lambda = seq(0.01, 1.6, by = 0.005)) |>
  mutate(d = dlnorm(lambda, meanlog = log(theta), sdlog = omega))

p_dens <- ggplot(dens, aes(lambda, d)) +
  geom_area(fill = col_placebo, alpha = 0.18) +
  geom_line(colour = col_placebo, linewidth = 0.7) +
  geom_segment(data = people,
               aes(x = lambda, xend = lambda, y = 0,
                   yend = dlnorm(lambda, log(theta), omega)),
               colour = col_hu, linewidth = 0.7) +
  geom_point(data = people, aes(lambda, dlnorm(lambda, log(theta), omega)),
             colour = col_hu, size = 2.4) +
  labs(title = "Between patients: lambda_i = theta * exp(eta_i)",
       subtitle = "theta = 0.375 per month, omega = 0.6; markers at the 5th, 50th, 95th percentiles",
       x = "Individual mean count, lambda_i (VOCs per month)", y = "Density")

pmf <- people |>
  crossing(n = 0:5) |>
  mutate(prob = dpois(n, lambda))

p_pmf <- ggplot(pmf, aes(n, prob)) +
  geom_col(fill = col_hu, width = 0.6) +
  facet_wrap(~ panel) +
  scale_x_continuous(breaks = 0:5) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Within a patient: Y_ij ~ Poisson(lambda_i)",
       x = "VOCs in a month, n", y = "P(Y = n)")

save_fig(p_dens / p_pmf, "fig06_mixed_effects.png", height = 6.5)
