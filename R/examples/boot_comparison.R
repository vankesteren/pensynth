# basic bootstrap comparison
library(tidyverse)
devtools::load_all()
NBOOT <- 1000
LAMBDA <- 1e-1
NDONOR <- 10
PHI <- 0
sim_ar <- function(n = 10, mu = 0, sd = 0.1, phi = 0.1) {
  x <- numeric(n)
  x[1] <- rnorm(1, mu, sd)
  for (i in 2:n)
    x[i] <- mu + phi*(x[i-1]-mu) + rnorm(1, 0, sqrt(sd^2 - sd^2 * phi^2))
  return(x)
}

simfun <- function() {
  if (PHI == 0)
    return(simulate_data(treatment_effect = 1, N_donor = NDONOR))
  dat <- simulate_data(treatment_effect = 1, N_donor = NDONOR, sd_resid_X1 = 0, sd_resid_Y1 = 0, sd_resid_Z1 = 0)
  resid_X1 <- sim_ar(n = nrow(dat$X1), mu = 0, sd = 0.1, phi = PHI)
  resid_Y1 <- sim_ar(n = nrow(dat$Y1), mu = 0, sd = 0.1, phi = PHI)
  resid_Z1 <- sim_ar(n = nrow(dat$Z1), mu = 0, sd = 0.1, phi = PHI)
  dat$X1 <- dat$X1 + resid_X1
  dat$Y1 <- dat$Y1 + resid_Y1
  dat$Z1 <- dat$Z1 + resid_Z1
  return(dat)
}

# first, true sampling variance
set.seed(45)
fit_samplingv <- lapply(1:NBOOT, \(j) {
  d <- simfun()
  pensynth(d$X1, d$X0, lambda = LAMBDA)
})

# then, the regular bootstrap
set.seed(45)
d <- simfun()
fit_bootstrap <- lapply(1:NBOOT, \(j) {
  i <- sample(ncol(d$X0), replace = TRUE)
  pensynth(d$X1, d$X0[,i], lambda = LAMBDA)
})

# then, the weighted bootstrap
set.seed(45)
d <- simfun()
fit_bayesboot <- lapply(1:NBOOT, \(j) {
  a <- rgamma(ncol(d$X0) + 1, 1, 1)
  pensynth(d$X1, d$X0, lambda = LAMBDA, iweights = a / sum(a))
})

# now, standard regression weights
set.seed(45)
d <- simfun()
f <- pensynth(d$X1, d$X0, lambda = LAMBDA)
rss <- sum((d$X1 - predict(f, d$X0))^2)
ses <- sqrt(diag(rss*solve(crossprod(d$X0) + LAMBDA * diag(ncol(d$X0)))))


# let's look at weights (NB: first 4 are true weights)
w_samplingv <- sapply(fit_samplingv, weights)
w_bootstrap <- sapply(fit_bootstrap, weights)
w_bayesboot <- sapply(fit_bayesboot, weights)
w_leastsqrs <- replicate(NBOOT, rnorm(length(ses), weights(f), ses))

boot_weights <-
  bind_rows(
    "samplingv" = as_tibble(w_samplingv) |> rownames_to_column(var = "donor") |> pivot_longer(-donor),
    "bootstrap" = as_tibble(w_bootstrap) |> rownames_to_column(var = "donor") |> pivot_longer(-donor),
    "bayesboot" = as_tibble(w_bayesboot) |> rownames_to_column(var = "donor") |> pivot_longer(-donor),
    "leastsqrs" = as_tibble(w_leastsqrs) |> rownames_to_column(var = "donor") |> pivot_longer(-donor),
    .id = "method"
  ) |>
  mutate(name = parse_number(name), across(-value, as_factor)) |>
  rename(iter = name)

boot_weights |>
  filter(as.integer(donor) < 13) |>
  ggplot(aes(x = log(value), fill = as.integer(donor) < 5)) +
  geom_histogram(bins = 40) +
  facet_grid(rows = vars(method), cols = vars(donor)) +
  labs(x = "log(weight)", fill = "True donor")


# let's look at actual effect estimate
set.seed(45)
ate_samplingv <- sapply(1:NBOOT, \(j) {
  d <- simfun()
  mean(d$Y1 - predict(pensynth(d$X1, d$X0, lambda = LAMBDA), newdata = d$Y0))
})
ate_bootstrap <- sapply(fit_bootstrap, \(x) mean(d$Y1 - predict(x, newdata = d$Y0)))
ate_bayesboot <- sapply(fit_bayesboot, \(x) mean(d$Y1 - predict(x, newdata = d$Y0)))
ate_leastsqrs <- apply(w_leastsqrs,2, \(w) mean(d$Y1 - d$Y0%*%w))

# conformal ate
zprd_fit <- d$Z1 - predict(f, d$Z0)
zprd_bootstrap <- sapply(fit_bootstrap, \(x) d$Z1 - predict(x, d$Z0))
zprd_bayesboot <- sapply(fit_bayesboot, \(x) d$Z1 - predict(x, d$Z0))
zprd_leastsqrs <- apply(w_leastsqrs, 2, \(w) d$Z1 - d$Z0%*%w)

cte_fit <- replicate(NBOOT, mean(d$Y1 - predict(f, d$Y0) + sample(zprd_fit, nrow(d$Y1), TRUE)))
cte_bootstrap <- sapply(fit_bootstrap, \(x) mean(d$Y1 - predict(x, d$Y0) + sample(zprd_bootstrap, nrow(d$Y1), TRUE)))
cte_bayesboot <- sapply(fit_bayesboot, \(x) mean(d$Y1 - predict(x, d$Y0) + sample(zprd_bayesboot, nrow(d$Y1), TRUE)))
cte_leastsqrs <- apply(w_leastsqrs, 2, \(w) mean(d$Y1 - d$Y0%*%w + sample(zprd_leastsqrs, nrow(d$Y1), TRUE)))

boot_ate <-
  bind_rows(
    "samplingv" = tibble(iter = 1:NBOOT, value = ate_samplingv),
    "bootstrap" = tibble(iter = 1:NBOOT, value = ate_bootstrap),
    "bayesboot" = tibble(iter = 1:NBOOT, value = ate_bayesboot),
    "leastsqrs" = tibble(iter = 1:NBOOT, value = ate_leastsqrs),
    "conformal" = tibble(iter = 1:NBOOT, value = cte_fit),
    "cbootstrp" = tibble(iter = 1:NBOOT, value = cte_bootstrap),
    "cbayesbot" = tibble(iter = 1:NBOOT, value = cte_bayesboot),
    "clestsqrs" = tibble(iter = 1:NBOOT, value = cte_leastsqrs),
    .id = "method"
  ) |> mutate(method = as_factor(method))

boot_ate |>
  ggplot(aes(x = value, fill = method)) +
  geom_histogram(bins = 60) +
  facet_grid(rows = vars(method)) +
  labs(x = "ATE")


boot_ate |> summarize(mean(value), sd(value), .by = "method")




