# Figure 1 -- counts over time, and the same counts as a probability mass function
source("R/00_setup.R")
trial <- read.csv("data/voc_trial_sim.csv")

placebo <- filter(trial, arm == "Placebo")
set.seed(1)
shown   <- filter(placebo, id %in% sample(unique(id), 25))
y_max   <- max(placebo$count)

p_time <- ggplot(shown, aes(month, count, group = id)) +
  geom_line(colour = col_placebo, alpha = 0.25, linewidth = 0.4) +
  geom_point(colour = col_placebo, alpha = 0.35, size = 1.2,
             position = position_jitter(width = 0.15, height = 0, seed = 1)) +
  scale_x_continuous(breaks = seq(0, 24, 6)) +
  scale_y_continuous(breaks = 0:y_max, limits = c(-0.5, y_max + 0.5)) +
  labs(title = "Individual counts over time",
       subtitle = "25 placebo patients",
       x = "Time (months)", y = "VOCs per month")

pmf <- placebo |>
  count(count, name = "n") |>
  mutate(pct = 100 * n / sum(n))

p_pmf <- ggplot(pmf, aes(pct, count)) +
  geom_col(orientation = "y", fill = col_placebo, width = 0.6) +
  geom_text(aes(label = ifelse(pct < 1, "<1%", sprintf("%.0f%%", pct))), hjust = -0.15,
            size = 3.2, colour = col_ink) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  scale_y_continuous(breaks = 0:y_max, limits = c(-0.5, y_max + 0.5)) +
  labs(title = "Count frequency (pmf)",
       subtitle = "All placebo patient-months",
       x = "Frequency (%)", y = NULL)

save_fig(p_time + p_pmf + plot_layout(widths = c(2.2, 1)),
         "fig01_counts_and_pmf.png")
