# Figure 12 -- three models with the same mean of 0.375 VOCs per month
source("R/00_setup.R")
source("R/02_count_distributions.R")

lambda <- 4.5 / 12
pi0    <- 0.4
n      <- 0:4

pmf <- bind_rows(
  tibble(n, model = "Poisson",                         prob = dpois(n, lambda)),
  tibble(n, model = "Negative binomial (O = 1)",       prob = dnb(n, lambda, 1)),
  tibble(n, model = "Zero-inflated Poisson (pi = 0.4)", prob = dzip(n, lambda / (1 - pi0), pi0))
) |> mutate(model = factor(model, levels = unique(model)))

p <- ggplot(pmf, aes(n, prob, fill = model)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  scale_fill_manual(values = c(col_placebo, col_hu, col_third)) +
  labs(x = "VOCs in a month, n", y = "P(Y = n)",
       subtitle = "All three have mean 0.375. They differ in zeros and in the tail")

save_fig(p, "fig12_same_mean_different_shape.png", width = 7, height = 4.2)

k <- 0:200
pmf_full <- list(Poisson = dpois(k, lambda), NB = dnb(k, lambda, 1),
                 ZIP = dzip(k, lambda / (1 - pi0), pi0))
round(sapply(pmf_full, function(p) c(mean = sum(k * p),
                                     variance = sum((k - sum(k * p))^2 * p),
                                     p_zero = p[1], p_two_or_more = 1 - p[1] - p[2])), 3)
