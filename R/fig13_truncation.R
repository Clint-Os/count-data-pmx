# Figure 13 -- truncated count distributions and the regression to the mean they cause
source("R/00_setup.R")

truncate_left <- function(lambda, c_l, n = 0:200) {        # keep counts >= c_l
  p <- dpois(n, lambda) * (n >= c_l)
  p / sum(p)
}

# (a) zero truncation: a clinic register that only contains patients seen with a crisis
n_a <- 0:6
a <- bind_rows(
  tibble(n = n_a, what = "Poisson", prob = dpois(n_a, 1.125)),
  tibble(n = n_a, what = "Truncated", prob = truncate_left(1.125, 1)[n_a + 1])
)
# (b) left truncation at 3: enrolment requires at least 3 crises in the previous year
n_b <- 0:12
b <- bind_rows(
  tibble(n = n_b, what = "Poisson", prob = dpois(n_b, 4.5)),
  tibble(n = n_b, what = "Truncated", prob = truncate_left(4.5, 3)[n_b + 1])
)

bars <- function(d, title, subtitle, xlab) {
  ggplot(d, aes(n, prob, fill = what)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72) +
    scale_fill_manual(values = c(Poisson = col_placebo, Truncated = col_hu)) +
    scale_x_continuous(breaks = unique(d$n)[unique(d$n) %% 2 == 0]) +
    labs(title = title, subtitle = subtitle, x = xlab, y = "P(Y = n)")
}
p_a <- bars(a, "Zero truncation", "lambda = 1.125", "VOCs per quarter")
p_b <- bars(b, "Left truncation at 3", "lambda = 4.5", "VOCs in the screening year")

# (c) regression to the mean: expected screening-year count among those enrolled
rtm <- tibble(lambda = seq(1, 8, by = 0.1)) |>
  rowwise() |>
  mutate(`Screening year, enrolled patients` = sum((0:200) * truncate_left(lambda, 3)),
         `Any later year` = lambda) |>
  ungroup() |>
  pivot_longer(-lambda, names_to = "what", values_to = "expected")

p_c <- ggplot(rtm, aes(lambda, expected, colour = what)) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c("Screening year, enrolled patients" = col_hu,
                                 "Any later year" = col_placebo)) +
  scale_x_continuous(breaks = 1:8) +
  labs(title = "Selection inflates the baseline", subtitle = "Enrollment requires 3 or more",
       x = "Patient's true rate (VOCs per year)", y = "Expected VOCs per year") +
  guides(colour = guide_legend(ncol = 1))

save_fig(p_a + p_b + p_c + plot_layout(widths = c(1, 1.3, 1.2)),
         "fig13_truncation.png", width = 11, height = 4.4)

rtm |> filter(lambda %in% c(2.5, 4.5, 6)) |> pivot_wider(names_from = what, values_from = expected)
