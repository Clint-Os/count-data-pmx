# 03_simulate_violations.R -- four cohorts that differ only in how counts arise
# 100 untreated patients, 24 monthly VOC counts each, no time trend.
# The SAME patients (same lambda_i) are used in every scenario, so any difference
# between scenarios comes from the within-patient process alone.
source("R/00_setup.R")
source("R/02_count_distributions.R")
set.seed(2027)

n_id   <- 100
months <- 1:24
theta  <- 4.5 / 12      # typical VOCs per month
omega  <- 0.6           # between-patient SD on the log scale

omicron <- 1            # negative binomial overdispersion
pi0     <- 0.4          # zero-inflation probability
gamma   <- 1            # Markov: the rate doubles in the month after a crisis month

patients <- tibble(id = seq_len(n_id),
                   lambda_i = theta * exp(rnorm(n_id, 0, omega)))

simulate_markov <- function(lambda_i) {
  y    <- integer(length(months))
  prev <- rpois(1, lambda_i)                      # the unobserved month before the study
  for (j in seq_along(months)) {
    y[j] <- rpois(1, lambda_i * (1 + gamma * (prev > 0)))
    prev <- y[j]
  }
  y
}

cohort <- patients |>
  crossing(month = months) |>
  arrange(id, month) |>
  mutate(poisson = rpois(n(), lambda_i),
         negbin  = rnb(n(), lambda_i, omicron),
         zip     = rzip(n(), lambda_i / (1 - pi0), pi0)) |>
  group_by(id) |>
  mutate(markov = simulate_markov(first(lambda_i))) |>
  ungroup()

write.csv(cohort, "data/voc_violations_sim.csv", row.names = FALSE)

cohort |>
  pivot_longer(c(poisson, negbin, zip, markov), names_to = "scenario", values_to = "count") |>
  group_by(scenario) |>
  summarise(mean = mean(count), variance = var(count), pct_zero = 100 * mean(count == 0),
            max = max(count))
