library(tidyverse)
library(tibbletime)
library(tidyquant)
library(cmdstanr)
library(posterior)
library(bayesplot)

vix <-
  tq_get(c("^VIX"),
         get = "stock.prices",
         from = "1990-01-01",
         to = "2021-12-20") %>%
  arrange(date) %>%
  mutate(r = 100 * log(close/lag(close))) %>%
  # if any returns are exactly zero this causes the lp to evaluate to zero
  mutate(r = ifelse(r == 0.0, 1e-4, r)) %>%
  na.omit() %>%
  mutate(logVIX = log(close/sqrt(252)))

spx <-
  tq_get(c("^GSPC"),
         get = "stock.prices",
         from = "1990-01-01",
         to = "2021-12-20") %>%
  arrange(date) %>%
  mutate(r = 100 * log(close/lag(close))) %>%
  # if any returns are exactly zero this causes the lp to evaluate to zero
  mutate(r = ifelse(r == 0.0, 1e-4, r)) %>%
  na.omit() %>%
  mutate(logVIX = vix$logVIX)

model <- cmdstan_model("stan/kaeck_vov.stan")

dat <- list(T = nrow(df %>% filter(t %% 1 == 0)), logVIX_Obs = df %>% filter(t %% 1 == 0) %>% pull(logVIX), H = 50)

fit <-
  model$sample(
    data = dat,
    seed = 1994,
    iter_warmup = 2e2,
    iter_sampling = 2e2,
    chains = 4,
    parallel_chains = 4,
    refresh = 1e0,
    max_treedepth = 6,
    adapt_delta = 0.8,
    show_messages = FALSE,
    init = replicate(4, list(kappa = 0.011, theta = 3.073,
                             kappa_nu = 0.11, theta_nu = 0.349*1e-2, sigma_nu = 0.183*1e-1,
                             rho = 0.183,
                             V0 = 0.35*1e-2,
                             W_raw = matrix(0.0, dat$T, dat$H-1),
                             W_nu = matrix(0.0, dat$T, dat$H-1)), simplify = FALSE)
  )


# Basic VoV model
model <- cmdstan_model("stan/vov.stan")
model <- cmdstan_model("stan/vov_spx_test.stan")

dat <- list(T = nrow(spx), r = spx$r)

fit <-
  model$sample(
    data = dat,
    seed = 1994,
    iter_warmup = 2e2,
    iter_sampling = 2e2,
    chains = 4,
    parallel_chains = 4,
    refresh = 1e0,
    max_treedepth = 8,
    adapt_delta = 0.8,
    show_messages = FALSE,
    init = replicate(4, list(rp = 0.1, theta = 0.9, tau = 0.5, rho = -0.7,
                             s0 = 0.9, phi = 0.99, sigma = 0.21, rho_g = 0.5), simplify = FALSE)
  )

g <- fit$summary("g") %>% mutate(date = spx$date)
u <- fit$summary("u") %>% mutate(date = spx$date)

h <- fit$summary("h") %>% mutate(date = spx$date)
v <- fit$summary("v") %>% mutate(date = spx$date)

v_raw <- fit$summary("v_raw") %>% mutate(date = spx$date)
u_raw <- fit$summary("u_raw") %>% mutate(date = spx$date)

y_rep <- fit$summary("y_rep") %>% mutate(date = vix$date)

u_raw %>%
  filter(year(date) == 2008) %>%
  ggplot(aes(date, median)) + geom_point() +
  geom_errorbar(aes(ymin = q5, ymax = q95))

h_rep %>%
  mutate(h = log(vix$close)) %>%
  filter(year(date) == 2008) %>%
  ggplot(aes(date, median)) +
  geom_line() +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.5) +
  geom_line(aes(date, h), color = "red")

h_rep %>%
  mutate(h = log(vix$close)) %>%
  filter(year(date) == 2008, month(date) >= 8) %>%
  ggplot(aes(date, h - median)) +
  geom_point() +
  geom_errorbar(aes(ymin = q5 - median, ymax = q95-median), alpha = 0.5)

u_raw_draws <- fit$draws("u_raw") %>% as_draws_matrix()

# Forward sim with same parameters
vix_holdout <-
  tq_get(c("^VIX"),
         get = "stock.prices",
         from = "2010-05-01",
         to = "2021-12-20") %>%
  arrange(date) %>%
  mutate(r = 100 * log(close/lag(close))) %>%
  # if any returns are exactly zero this causes the lp to evaluate to zero
  mutate(r = ifelse(r == 0.0, 1e-4, r)) %>%
  na.omit() 

theta <- 0.93
tau <- 0.32
rho <- 0.7
mu <- 2.62
phi <- 0.99
sigma <- 0.06

g_draws <- fit$draws("g[8054]") %>% as_draws_matrix()
g0_draws <- g_draws[,1]

sample_day <- function(g0, h0, N) {
  g <- g0
  h <- h0
  for(t in 1:N) {
    u <- tau*rnorm(1)
    g <- theta*g + u
    h <- mu + phi*(h-mu) + rho*(u/tau)*sigma*exp(g/2) + sigma*exp(g/2)*sqrt(1-rho^2)*rnorm(1)
  }
  
  h
}

sample_day(g0_draws[1], 3.437851, 1)

g30 <- map_dbl(g0_draws, sample_day, h0 = 3.071303, N = 5)
qplot(exp(g30))















