# Figure 10 -- Poissonness plot (Hoaglin 1980)
source("R/00_setup.R")
source("R/02_count_distributions.R")
set.seed(11)

lambda <- 4.5          # VOCs per year, identical for everyone in this illustration
n_obs  <- 2000         # patient-years per model

sims <- bind_rows(
  tibble(model = "Poisson",                      count = rpois(n_obs, lambda)),
  tibble(model = "Negative binomial (O = 0.5)",  count = rnb(n_obs, lambda, 0.5)),
  tibble(model = "Double Poisson (upsilon = 3)",
         count = r_from_pmf(n_obs, function(n) ddp(n, lambda, 3)))
) |>
  mutate(model = factor(model, levels = unique(model)))

# count metameter: log(n! * f_n / N); a straight line in n if the data are Poisson
metameter <- sims |>
  count(model, count, name = "f_n") |>
  group_by(model) |>
  mutate(N = sum(f_n), phi = log(f_n) + lgamma(count + 1) - log(N)) |>
  filter(f_n >= 5) |>
  ungroup()

p <- ggplot(metameter, aes(count, phi)) +
  geom_abline(intercept = -lambda, slope = log(lambda), colour = col_muted,
              linetype = "dashed") +
  geom_point(colour = col_placebo, size = 2) +
  facet_wrap(~ model) +
  labs(x = "Count value n (VOCs per year)", y = "Count metameter",
       subtitle = "Dashed line: what a Poisson distribution with lambda = 4.5 predicts (intercept -4.5, slope log 4.5)")

save_fig(p, "fig10_poissonness.png", width = 9, height = 3.8)

metameter |> group_by(model) |>
  summarise(intercept = coef(lm(phi ~ count))[1], slope = coef(lm(phi ~ count))[2])
