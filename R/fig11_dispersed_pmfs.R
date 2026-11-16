# Figure 11 -- under- and overdispersed pmfs with the same mean as a Poisson
source("R/00_setup.R")
source("R/02_count_distributions.R")

lambda <- 4.5
n <- 0:16

under <- bind_rows(
  tibble(n, model = "Poisson",                          prob = dpois(n, lambda)),
  tibble(n, model = "Double Poisson (upsilon = 3)",     prob = ddp(n, lambda, 3)),
  tibble(n, model = "Generalized Poisson (delta = -0.5)", prob = dgp(n, lambda, -0.5))
) |> mutate(panel = "Underdispersion: variance < mean")

over <- bind_rows(
  tibble(n, model = "Poisson",                          prob = dpois(n, lambda)),
  tibble(n, model = "Negative binomial (O = 0.5)",      prob = dnb(n, lambda, 0.5)),
  tibble(n, model = "Generalized Poisson (delta = 0.4)", prob = dgp(n, lambda, 0.4))
) |> mutate(panel = "Overdispersion: variance > mean")

make_panel <- function(d, colours) {
  d$model <- factor(d$model, levels = names(colours))
  ggplot(d, aes(n, prob, colour = model, shape = model)) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 2) +
    scale_colour_manual(values = colours) +
    scale_x_continuous(breaks = seq(0, 16, 4)) +
    scale_y_continuous(limits = c(0, 0.36)) +
    labs(title = unique(d$panel), x = "VOCs per year, n", y = "P(Y = n)") +
    guides(colour = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))
}

p_under <- make_panel(under, c("Poisson" = col_ink,
                               "Double Poisson (upsilon = 3)" = col_placebo,
                               "Generalized Poisson (delta = -0.5)" = col_third))
p_over  <- make_panel(over,  c("Poisson" = col_ink,
                               "Negative binomial (O = 0.5)" = col_hu,
                               "Generalized Poisson (delta = 0.4)" = col_third))

save_fig(p_under + p_over, "fig11_dispersed_pmfs.png", width = 9, height = 4.8)

# mean and variance of each pmf
summ <- function(p, k = 0:200) { m <- sum(k * p); c(mean = m, variance = sum((k - m)^2 * p)) }
k <- 0:200
round(rbind(Poisson = summ(dpois(k, lambda)),
            DP_3    = summ(ddp(k, lambda, 3)),
            GP_neg  = summ(dgp(k, lambda, -0.5)),
            NB_0.5  = summ(dnb(k, lambda, 0.5)),
            GP_pos  = summ(dgp(k, lambda, 0.4))), 2)
