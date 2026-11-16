# Figure 15 -- a non-homogeneous Poisson process: hazard, cumulative hazard, survival, counts
source("R/00_setup.R")
set.seed(15)

h0    <- 4.5 / 12                 # VOCs per month before treatment
e_max <- 0.50
k_on  <- log(2) / 2

h <- function(t) h0 * (1 - e_max * (1 - exp(-k_on * t)))                 # hazard
H <- function(t) h0 * ((1 - e_max) * t + e_max * (1 - exp(-k_on * t)) / k_on)  # its integral

grid <- tibble(t = seq(0, 24, by = 0.05)) |>
  mutate(hazard = h(t), cum_hazard = H(t), survival = exp(-H(t)),
         survival_untreated = exp(-h0 * t))

p_h <- ggplot(grid, aes(t, hazard)) + geom_line(colour = col_hu, linewidth = 0.8) +
  scale_y_continuous(limits = c(0, 0.4)) + scale_x_continuous(breaks = seq(0, 24, 6)) +
  labs(title = "Hazard h(t)", x = NULL, y = "VOCs per month")
p_H <- ggplot(grid, aes(t, cum_hazard)) + geom_line(colour = col_hu, linewidth = 0.8) +
  scale_x_continuous(breaks = seq(0, 24, 6)) +
  labs(title = "Cumulative hazard H(t)", x = NULL, y = "Expected VOCs so far")
p_S <- ggplot(grid, aes(t)) +
  geom_line(aes(y = survival_untreated, linetype = "Untreated (constant hazard)"),
            colour = col_muted, linewidth = 0.6) +
  geom_line(aes(y = survival, linetype = "On hydroxyurea"), colour = col_hu, linewidth = 0.8) +
  scale_linetype_manual(values = c("On hydroxyurea" = "solid",
                                   "Untreated (constant hazard)" = "dashed")) +
  scale_x_continuous(breaks = seq(0, 24, 6)) +
  labs(title = "Survival S(t) = exp(-H(t))", x = "Time since starting treatment (months)",
       y = "P(no VOC yet)") +
  guides(linetype = guide_legend(ncol = 1)) +
  theme(legend.position = "bottom")

# expected counts per 6-month interval, three ways
ends <- seq(6, 24, by = 6)
expected <- tibble(end = ends, start = ends - 6) |>
  mutate(`Exact: integral of h(t)`        = H(end) - H(start),
         `Hazard at interval start x 6`   = h(start) * 6,
         `Hazard at interval end x 6`     = h(end) * 6)

# check the integral by simulation: thinning a homogeneous process of rate h0
simulate_times <- function() {
  t <- cumsum(rexp(60, rate = h0)); t <- t[t <= 24]
  t[runif(length(t)) < h(t) / h0]
}
sim_counts <- replicate(20000, tabulate(ceiling(simulate_times() / 6), nbins = 4))
expected$simulated <- rowMeans(sim_counts)

bars <- expected |>
  select(-simulated) |>
  pivot_longer(-c(start, end), names_to = "method", values_to = "count") |>
  mutate(interval = factor(sprintf("%d-%d", start, end), levels = sprintf("%d-%d", ends - 6, ends)),
         method = factor(method, levels = c("Hazard at interval start x 6",
                                            "Exact: integral of h(t)",
                                            "Hazard at interval end x 6")))
p_N <- ggplot(bars, aes(interval, count, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  geom_point(data = mutate(expected, interval = factor(sprintf("%d-%d", start, end))),
             aes(interval, simulated, shape = "Simulated mean (20,000 patients)"),
             inherit.aes = FALSE, size = 2.2, colour = col_ink) +
  scale_fill_manual(values = c(col_third, col_hu, col_placebo)) +
  labs(title = "Expected VOCs per 6-month interval", x = "Interval (months)", y = "Expected count") +
  guides(fill = guide_legend(ncol = 1, order = 1), shape = guide_legend(order = 2)) +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(0, 0, 0, 0))

save_fig((p_h + p_H) / (p_S + p_N) + plot_layout(heights = c(1, 1.5)),
         "fig15_nonhomogeneous.png", width = 10, height = 8.5)
print(as.data.frame(expected))

# how much does "hazard at the end of each month" (used in Part 1) miss with monthly intervals?
m <- 1:24
c(exact_24_months = H(24), end_of_month_sum = sum(h(m)), start_of_month_sum = sum(h(m - 1)))
