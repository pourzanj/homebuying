library(rstan)

fit <- stan("minimal_bimodal_horseshoe.stan", data = list(y = 100, sigma = 30), chains = 10, iter = 2e3)
samples <- extract(fit)
samples$lambda %>% qplot() + scale_x_log10()

y <- 1e2
s <- 30
f <- function(l) dnorm(y, scale = sqrt(s^2 + l^2))

tibble(lambda = seq(0, 2e2, 0.01)) %>%
  mutate(cauchy_pdf = 2 * dcauchy(lambda),
         normal_pdf = dnorm(y, sd = sqrt(s^2 + lambda^2))) %>%
  mutate(p = cauchy_pdf * normal_pdf) %>%
  ggplot(aes(lambda, p)) +
  geom_line() +
  scale_x_log10()
