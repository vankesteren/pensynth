library(tidyverse)
library(pensynth)
library(pbapply)
library(parallel)

# create grid
grid <- expand_grid(
  iter = 1:2000,
  N_donors = round(exp(seq(log(20), log(1000), len = 20))),
  N_covariates = round(exp(seq(log(1), log(25), len = 8))),
  weights = c(FALSE, TRUE)
)

# create result function
inch <- function(i) {
  # get param settings
  Nd <- grid[i,][["N_donors"]]
  Nc <- grid[i,][["N_covariates"]]
  Wt <- grid[i,][["weights"]]

  # generate data
  X0 <- matrix(rnorm(Nc*Nd), nrow = Nc)
  if (Wt) {
    w <- runif(Nd)
    w[5:Nd] <- 0
    w <- w / sum(w)
    X1 <- X0 %*% w + rnorm(Nc, sd = 0.1)
  } else {
    X1 <- matrix(rnorm(Nc), nrow = Nc)
  }

  # compute outcome
  in_convex_hull(X1, X0)
}

# make cluster
clus <- makeCluster(10)
clusterEvalQ(clus, library(pensynth))
clusterExport(clus, c("grid", "inch"))

# perform simulation (takes about 8 minutes)
grid$inch <- pbvapply(1:nrow(grid), inch, logical(1), cl = clus)
stopCluster(clus)

# analyse
grid_res <-
  grid |>
  summarize(
    prob = mean(inch, na.rm = TRUE),
    nacount = sum(is.na(inch)),
    .by = c(N_donors, N_covariates, weights)
  )

# fit logistic regression model
fit <- glm(
  prob ~ log(N_donors) + log(N_covariates) + weights + log(N_donors):weights + log(N_covariates):weights,
  family = "quasibinomial", data = grid_res, weights = rep(2000, nrow(grid_res))
)


# add predictions to grid_res for plotting
grid_res$pred <- predict(fit, grid_res, type = "response")

# save for easy adjustment
write_rds(grid_res, "img/grid_res.rds")

# plot
grid_res |>
  mutate(w = factor(ifelse(weights, "Synthetic control weights", "Random draw"))) |>
  ggplot(aes(x = N_donors, y = prob, colour = as_factor(N_covariates))) +
  geom_line(aes(y = pred), alpha = 0.2) +
  geom_line() +
  geom_point() +
  scale_colour_brewer(type = "qual", palette = 2) +
  scale_x_log10() +
  theme_linedraw() +
  facet_wrap(vars(w)) +
  labs(
    x = "# donor units",
    y = "Probability of non-unique solution",
    colour = "# covariates",
    title = "Probability of non-unique solution in synthetic control method",
    subtitle = "Probability increases with donors and decreases with covariates"
  )

# save the plot
ggsave("img/nonunique.png", width = 10, height = 6)
