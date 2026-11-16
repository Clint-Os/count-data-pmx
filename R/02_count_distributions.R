# 02_count_distributions.R -- pmfs and random generators for the models in Part 2
# All are parameterized by the mean count `lambda` plus one dispersion parameter,
# following Plan (2014), so that each collapses to the Poisson model.

# Negative binomial: Var = lambda * (1 + omicron * lambda); Poisson as omicron -> 0
dnb <- function(n, lambda, omicron) dnbinom(n, mu = lambda, size = 1 / omicron)
rnb <- function(k, lambda, omicron) rnbinom(k, mu = lambda, size = 1 / omicron)

# Generalized Poisson (Consul & Jain 1973): Var = lambda / (1 - delta)^2; Poisson at delta = 0
# For delta < 0 the pmf is cut off where lambda*(1-delta) + n*delta <= 0, then renormalized.
dgp <- function(n, lambda, delta, n_max = 200) {
  raw <- function(k) {
    a <- lambda * (1 - delta)
    b <- a + k * delta
    ifelse(b > 0, exp(log(a) + (k - 1) * log(pmax(b, 1e-300)) - b - lgamma(k + 1)), 0)
  }
  raw(n) / sum(raw(0:n_max))
}

# Double Poisson (Efron 1986): mean ~ lambda, Var ~ lambda / upsilon; Poisson at upsilon = 1
# Normalized numerically; Efron's closed-form constant is an approximation to this.
ddp <- function(n, lambda, upsilon, n_max = 200) {
  raw <- function(k) {
    k_log_k <- ifelse(k == 0, 0, k * log(k))
    exp(0.5 * log(upsilon) - upsilon * lambda - k + k_log_k - lgamma(k + 1) +
          upsilon * k * (1 + log(lambda)) - upsilon * k_log_k)
  }
  raw(n) / sum(raw(0:n_max))
}

# Zero-inflated Poisson: a structural zero with probability pi0, otherwise Poisson(lambda_c)
dzip <- function(n, lambda_c, pi0) pi0 * (n == 0) + (1 - pi0) * dpois(n, lambda_c)
rzip <- function(k, lambda_c, pi0) ifelse(runif(k) < pi0, 0L, rpois(k, lambda_c))

# Draw from any pmf defined on 0..n_max
r_from_pmf <- function(k, pmf, n_max = 200) {
  sample(0:n_max, size = k, replace = TRUE, prob = pmf(0:n_max))
}
