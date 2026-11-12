# 01_simulate_trial.R -- simulate the teaching trial used throughout the series
# A hypothetical 24-month, two-arm trial of hydroxyurea in sickle cell anaemia.
# Outcome: number of vaso-occlusive crises (VOCs) per patient per month.
source("R/00_setup.R")
set.seed(2026)

n_per_arm <- 60
months    <- 1:24

theta  <- 4.5 / 12       # typical VOC rate per month (4.5 per year, MSH placebo arm)
omega  <- 0.6            # SD of the between-patient random effect (log scale)
p_max  <- 0.15           # maximal "placebo" / regression-to-the-mean reduction
k_p    <- log(2) / 3     # placebo effect half-time: 3 months
e_max  <- 0.50           # maximal hydroxyurea reduction (rate ratio 0.50)
k_on   <- log(2) / 2     # drug effect onset half-time: 2 months

subjects <- tibble(
  id  = seq_len(2 * n_per_arm),
  arm = rep(c("Placebo", "Hydroxyurea"), each = n_per_arm),
  eta = rnorm(2 * n_per_arm, mean = 0, sd = omega)
)

trial <- subjects |>
  crossing(month = months) |>
  mutate(
    f_time = p_max * (1 - exp(-k_p * month)),                      # f(t)
    g_drug = if_else(arm == "Hydroxyurea",
                     e_max * (1 - exp(-k_on * month)), 0),         # g(t)
    lambda = theta * exp(eta) * (1 - f_time) * (1 - g_drug),       # mean count
    count  = rpois(n(), lambda)                                    # observation
  ) |>
  mutate(arm = factor(arm, levels = c("Placebo", "Hydroxyurea"))) |>
  arrange(id, month)

dir.create("data", showWarnings = FALSE)
write.csv(trial, "data/voc_trial_sim.csv", row.names = FALSE)

# Quick look: events per patient-year by arm
trial |>
  group_by(arm) |>
  summarise(patients = n_distinct(id),
            events   = sum(count),
            rate_per_year = sum(count) / (n() / 12))
