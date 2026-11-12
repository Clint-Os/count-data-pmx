# Figure 4 -- the Poisson distribution as the limit of the binomial
source("R/00_setup.R")

lambda <- 4.5 / 4                      # expected VOCs in one quarter (3 months)
n_slots <- c(3, 12, 90)                # months, weeks, days in a quarter

binom <- crossing(n = n_slots, k = 0:6) |>
  mutate(prob  = dbinom(k, size = n, prob = lambda / n),
         panel = factor(sprintf("n = %d slots, p = %.4f", n, lambda / n),
                        levels = sprintf("n = %d slots, p = %.4f", n_slots, lambda / n_slots)))
pois <- binom |> mutate(prob = dpois(k, lambda))

p <- ggplot(binom, aes(k, prob)) +
  geom_col(aes(fill = "Binomial(n, p = lambda/n)"), width = 0.6) +
  geom_point(data = pois, aes(colour = "Poisson(lambda = 1.125)"), size = 2.4) +
  facet_wrap(~ panel) +
  scale_fill_manual(values = c("Binomial(n, p = lambda/n)" = col_placebo)) +
  scale_colour_manual(values = c("Poisson(lambda = 1.125)" = col_hu)) +
  scale_x_continuous(breaks = 0:6) +
  guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2)) +
  labs(x = "Number of events in the quarter", y = "Probability")

save_fig(p, "fig04_binomial_to_poisson.png", height = 3.8)
