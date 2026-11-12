# Figure 7 -- the Poisson model and the constant-hazard model are the same model
source("R/00_setup.R")
set.seed(42)

h <- 4.5 / 12                                   # constant hazard, per month

# (a) survival curve S(t) = exp(-h t), and P(Y = 0) from the Poisson pmf
surv   <- tibble(t = seq(0, 12, by = 0.05)) |> mutate(S = exp(-h * t))
marks  <- tibble(t = c(1, 3, 12)) |>
  mutate(p0    = dpois(0, lambda = h * t),
         label = sprintf("lambda = %.3f\nP(Y = 0) = %.3f", h * t, p0))

p_surv <- ggplot(surv, aes(t, S)) +
  geom_line(colour = col_ink, linewidth = 0.6) +
  geom_point(data = marks, aes(t, p0), colour = col_hu, size = 2.8) +
  geom_text(data = marks, aes(t, p0, label = label), hjust = c(-0.12, -0.12, 1.1),
            vjust = c(0.2, -0.3, -0.6), size = 3.1, colour = col_ink, lineheight = 0.95) +
  scale_x_continuous(breaks = seq(0, 12, 3)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Survival curve vs Poisson P(Y = 0)",
       subtitle = "Line: S(t) = exp(-h t). Dots: Poisson pmf at n = 0",
       x = "Interval length (months)", y = "Probability of no event")

# (b) simulate event TIMES with exponential gaps, count per month, compare with dpois
n_months <- 20000
t_ev   <- cumsum(rexp(2 * h * n_months, rate = h))
t_ev   <- t_ev[t_ev <= n_months]
counts <- tabulate(ceiling(t_ev), nbins = n_months)

emp <- tibble(n = 0:4) |>
  mutate(simulated = sapply(n, function(k) mean(counts == k)),
         poisson   = dpois(n, h))

p_sim <- ggplot(emp, aes(n)) +
  geom_col(aes(y = simulated, fill = "Counted from event times"), width = 0.6) +
  geom_point(aes(y = poisson, colour = "Poisson pmf"), size = 2.6) +
  scale_fill_manual(values = c("Counted from event times" = col_placebo)) +
  scale_colour_manual(values = c("Poisson pmf" = col_hu)) +
  labs(title = "Event times, counted per month",
       subtitle = "20,000 simulated patient-months, lambda = 0.375",
       x = "VOCs in a month, n", y = "Proportion of months") +
  guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2)) +
  theme(legend.position = "bottom")

save_fig(p_surv + p_sim, "fig07_poisson_hazard_link.png", width = 9)
print(emp)
