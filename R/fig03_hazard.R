# Figure 3 -- constant vs non-constant hazard, and what interval counting keeps
source("R/00_setup.R")

h0 <- 4.5 / 12
h_const <- function(t) rep(h0, length(t))
h_vary  <- function(t) h0 * (1 + 0.8 * sin(2 * pi * (t - 1) / 12))   # a seasonal swing

# average hazard over each interval of a given width = expected count / width
interval_mean <- function(h, width, t_end = 24) {
  starts <- seq(0, t_end - width, by = width)
  tibble(start = starts, end = starts + width,
         h_bar = sapply(starts, function(s) integrate(h, s, s + width)$value / width))
}

grid <- tibble(t = seq(0, 24, by = 0.05))
curves <- bind_rows(
  grid |> mutate(panel = "Constant hazard", h = h_const(t)),
  grid |> mutate(panel = "Non-constant hazard", h = h_vary(t))
)
steps <- bind_rows(
  interval_mean(h_const, 1) |> mutate(panel = "Constant hazard",     what = "Monthly counts"),
  interval_mean(h_const, 6) |> mutate(panel = "Constant hazard",     what = "6-monthly counts"),
  interval_mean(h_vary,  1) |> mutate(panel = "Non-constant hazard", what = "Monthly counts"),
  interval_mean(h_vary,  6) |> mutate(panel = "Non-constant hazard", what = "6-monthly counts")
) |>
  mutate(what = factor(what, levels = c("Monthly counts", "6-monthly counts")))

p <- ggplot() +
  geom_segment(data = filter(steps, what == "6-monthly counts"),
               aes(x = start, xend = end, y = h_bar, yend = h_bar, colour = what),
               linewidth = 2.2, alpha = 0.6) +
  geom_segment(data = filter(steps, what == "Monthly counts"),
               aes(x = start, xend = end, y = h_bar, yend = h_bar, colour = what),
               linewidth = 1.1) +
  geom_line(data = curves, aes(t, h, linetype = "True hazard h(t)"),
            colour = col_ink, linewidth = 0.6) +
  facet_wrap(~ panel) +
  scale_colour_manual(values = c("Monthly counts" = col_placebo,
                                 "6-monthly counts" = col_hu),
                      breaks = c("Monthly counts", "6-monthly counts")) +
  scale_linetype_manual(values = c("True hazard h(t)" = "solid")) +
  scale_x_continuous(breaks = seq(0, 24, 6)) +
  scale_y_continuous(limits = c(0, 0.75)) +
  labs(x = "Time (months)", y = "Hazard (events per month)",
       subtitle = "Coloured steps: the average rate that a count per interval can see")

save_fig(p, "fig03_hazard.png", height = 4.2)
