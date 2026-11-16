# Figure 14 -- diagnosing dependence between consecutive counts
# Rate ratio: mean count in a month that follows a crisis month,
#             divided by the mean count in a month that follows a crisis-free month.
source("R/00_setup.R")
cohort <- read.csv("data/voc_violations_sim.csv")

lagged <- cohort |>
  pivot_longer(c(poisson, markov), names_to = "model", values_to = "count") |>
  mutate(model = factor(ifelse(model == "poisson", "Poisson (independent)", "Markov Poisson"),
                        levels = c("Poisson (independent)", "Markov Poisson"))) |>
  group_by(model, id) |>
  arrange(month, .by_group = TRUE) |>
  mutate(previous = lag(count)) |>
  ungroup() |>
  filter(!is.na(previous)) |>
  mutate(after_crisis = previous > 0)

# (a) naive: pool all patient-months
pooled <- lagged |>
  group_by(model) |>
  summarise(rate_ratio = mean(count[after_crisis]) / mean(count[!after_crisis])) |>
  mutate(method = "Pooled over patients")

# (b) within patients: Mantel-Haenszel rate ratio, each patient is a stratum
within <- lagged |>
  group_by(model, id) |>
  summarise(a  = sum(count[after_crisis]),  n1 = sum(after_crisis),
            b  = sum(count[!after_crisis]), n0 = sum(!after_crisis), .groups = "drop") |>
  group_by(model) |>
  summarise(rate_ratio = sum(a * n0 / (n1 + n0)) / sum(b * n1 / (n1 + n0))) |>
  mutate(method = "Within patients")

both <- bind_rows(pooled, within) |>
  mutate(method = factor(method, levels = c("Pooled over patients", "Within patients")))

p <- ggplot(both, aes(model, rate_ratio, fill = method)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = col_muted) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  geom_text(aes(label = sprintf("%.2f", rate_ratio)),
            position = position_dodge(width = 0.7), vjust = -0.5, size = 3.4, colour = col_ink) +
  scale_fill_manual(values = c("Pooled over patients" = col_placebo,
                               "Within patients" = col_hu)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Rate ratio",
       title = "VOC rate after a crisis month relative to after a crisis-free month",
       subtitle = "True values used in the simulation: 1 for the Poisson cohort, 2 for the Markov cohort")

save_fig(p, "fig14_dependence.png", width = 7.5, height = 4.6)
print(both)
