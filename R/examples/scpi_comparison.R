# scpi comparison of boot
library(scpi)


scpi_data <-
  scdata(
    df = scpi_germany,
    id.var = "country",
    time.var = "year",
    outcome.var = "gdp",
    period.pre = (1960:1990),
    period.post = (1991:2003),
    unit.tr = "West Germany",
    unit.co = setdiff(unique(scpi_germany$country), "West Germany"),
    constant = FALSE,
    cointegrated.data = TRUE
  )

X1 <- scpi_data$A
X0 <- scpi_data$B

Y1 <- scpi_data$Y.post
Y0 <- scpi_data$P


res_scpi <- scpi(scpi_data, w.constr = list(name = "simplex", Q = 1), sims = 2000)
pred_scpi <- cbind(
  res_scpi$est.results$Y.post.fit,
  res_scpi$inference.results$CI.all.gaussian[,1:2]
)

res_pensynth <- pensynth(X1, X0, lambda = 1e-5)
bb_pensynth <- sapply(1:2000, \(i) pensynth(X1, X0, lambda = 1e-5, boot_weight = TRUE)$w)
Yhat_pensynth <- sapply(1:2000, \(i) Y0%*%bb_pensynth[,i])
pred_pensynth <- cbind(
  Y0%*%res_pensynth$w,
  t(apply(Yhat_pensynth, 1, \(x) quantile(x, probs = c(0.025, 0.975))))
)
library(tidyverse)

unname(rbind(pred_scpi, as.matrix(pred_pensynth))) |>
  as_tibble() |>
  set_names(c("est", "lwr", "upr")) |>
  mutate(method = as_factor(rep(c("scpi", "pensynth"), each = 13)),
         time = rep(1991:2003, 2)) |>
  ggplot(aes(x = time, color = method)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = method), alpha = 0.3, color = "transparent") +
  geom_line(aes(y = est)) +
  geom_line(data = tibble(time = 1991:2003, est = Y1, method = "observed"), aes(y = est), colour = "black") +
  theme_minimal()
