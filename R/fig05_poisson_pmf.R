# Figure 5 -- one hazard, three interval lengths, three Poisson pmfs
source("R/00_setup.R")

h <- 4.5 / 12                                   # events per month
intervals <- tibble(label = c("1 month", "3 months", "12 months"),
                    dt    = c(1, 3, 12)) |>
  mutate(lambda = h * dt,
         panel  = sprintf("%s: lambda = %.3f", label, lambda),
         panel  = factor(panel, levels = panel))

pmf <- intervals |>
  crossing(n = 0:12) |>
  mutate(prob = dpois(n, lambda))

p <- ggplot(pmf, aes(n, prob)) +
  geom_col(fill = col_placebo, width = 0.6) +
  geom_vline(aes(xintercept = lambda), colour = col_hu, linetype = "dashed",
             linewidth = 0.5) +
  facet_wrap(~ panel) +
  scale_x_continuous(breaks = seq(0, 12, 2)) +
  labs(x = "Count n", y = "P(Y = n)",
       subtitle = "Dashed line: lambda, which is both the mean and the variance")

save_fig(p, "fig05_poisson_pmf.png", height = 3.8)
