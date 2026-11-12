# Figure 8 -- lambda as a function of time and treatment
source("R/00_setup.R")
trial <- read.csv("data/voc_trial_sim.csv") |>
  mutate(arm = factor(arm, levels = c("Placebo", "Hydroxyurea")))
set.seed(7)

theta <- 4.5 / 12; omega <- 0.6
p_max <- 0.15; k_p  <- log(2) / 3
e_max <- 0.50; k_on <- log(2) / 2

# model-based population mean: E[theta * exp(eta)] = theta * exp(omega^2 / 2)
model <- crossing(arm = levels(trial$arm), month = seq(0.5, 24, by = 0.1)) |>
  mutate(arm    = factor(arm, levels = levels(trial$arm)),
         f_time = p_max * (1 - exp(-k_p * month)),
         g_drug = if_else(arm == "Hydroxyurea", e_max * (1 - exp(-k_on * month)), 0),
         mean_count = theta * exp(omega^2 / 2) * (1 - f_time) * (1 - g_drug))

# 95% prediction interval for the arm mean: re-simulate the whole trial 500 times
design <- distinct(trial, id, arm) |> crossing(month = 1:24) |>
  mutate(f_time = p_max * (1 - exp(-k_p * month)),
         g_drug = if_else(arm == "Hydroxyurea", e_max * (1 - exp(-k_on * month)), 0))
n_id <- n_distinct(design$id)

one_trial <- function(rep) {
  eta <- rnorm(n_id, 0, omega)
  design |>
    mutate(count = rpois(n(), theta * exp(eta[id]) * (1 - f_time) * (1 - g_drug))) |>
    group_by(arm, month) |>
    summarise(mean_count = mean(count), .groups = "drop")
}
pi95 <- bind_rows(lapply(1:500, one_trial)) |>
  group_by(arm, month) |>
  summarise(lo = quantile(mean_count, 0.025), hi = quantile(mean_count, 0.975),
            .groups = "drop")

observed <- trial |>
  group_by(arm, month) |>
  summarise(mean_count = mean(count), .groups = "drop")

p <- ggplot(mapping = aes(month, colour = arm, fill = arm)) +
  geom_ribbon(data = pi95, aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(data = model, aes(y = mean_count), linewidth = 0.8) +
  geom_point(data = observed, aes(y = mean_count), size = 1.8) +
  scale_colour_manual(values = c(Placebo = col_placebo, Hydroxyurea = col_hu)) +
  scale_fill_manual(values = c(Placebo = col_placebo, Hydroxyurea = col_hu)) +
  scale_x_continuous(breaks = seq(0, 24, 6)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Time (months)", y = "Mean VOCs per patient per month",
       subtitle = "Dots: observed arm means. Lines: model mean. Bands: 95% prediction interval (500 simulated trials)")

save_fig(p, "fig08_mean_count_function.png", height = 4.4)
