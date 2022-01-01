library(tidyverse)
library(tidyquant)
library(cmdstanr)

##############
##### Verify that v_raw wants to be big after observing -3.35 drop on 2020-02-24

spy <-
  tq_get(c("SPY"),
         get = "stock.prices",
         from = "2019-08-20",
         to = "2020-02-21") %>%
  arrange(date) %>%
  mutate(pct_return = 100 * (close - lag(close)) / lag(close)) %>%
  # if any returns are exactly zero this causes the lp to evaluate to zero
  mutate(pct_return = ifelse(pct_return == 0.0, 1e-4, pct_return)) %>%
  na.omit() 

model <- cmdstan_model("stan/sv.stan")
model <- cmdstan_model("stan/sv_rhs.stan")
model <- cmdstan_model("stan/sv_rhs_ar_tau.stan")
model <- cmdstan_model("stan/multiscale_sv.stan")
model <- cmdstan_model("stan/sv_arg.stan")

dat <- list(T = nrow(spy), y = spy$pct_return, tau0 = 5e-1, nu = 50, scale = 0.6, tau0_mu = 1e-1, nu_mu = 50, scale_mu = 0.6)
dat <- list(T = nrow(spy), y = spy$pct_return, K = 5)

fit <-
  model$sample(
    data = dat,
    seed = 1994,
    iter_warmup = 1e3,
    iter_sampling = 1e3,
    chains = 4,
    parallel_chains = 4,
    refresh = 2e1,
    max_treedepth = 12,
    adapt_delta = 0.9999,
    show_messages = FALSE
  )

fit$summary("v") %>%
  mutate(date = spy$date) %>%
  ggplot(aes(date, median)) +
  geom_point() +
  geom_errorbar(aes(ymin = q5, ymax = q95))

fit$draws("v[126]") %>% as_draws_df()%>% mutate(s = exp(`v[126]`/2))  %>% ggplot(aes(`v[126]`)) + geom_histogram()
fit$draws("y_rep[126]") %>% as_draws_df() %>% ggplot(aes(`y_rep[126]`)) + geom_histogram()

fit$draws("y_rep_ahead") %>% as_draws_df() %>% ggplot(aes(y_rep_ahead)) + geom_histogram(binwidth = 0.25)
fit$draws("y_rep_ahead") %>% as_draws_df() %>% as_tibble() %>% summarize(p = mean(y_rep_ahead <= -3.32))

fit$draws("y_rep_ahead") %>% as_draws_df() %>% as_tibble() %>% pull(y_rep_ahead) -> y_rep
mean(log(1+1*y_rep/100))
mean(log(1+2*y_rep/100))
mean(log(1+3*y_rep/100))

fit$summary("s") %>%
  mutate(date = spy$date) %>%
  ggplot(aes(date, median)) +
  geom_point() +
  geom_errorbar(aes(ymin = q5, ymax = q95))

fit$summary("log_lik") %>%
  mutate(date = spy$date) %>%
  ggplot(aes(date, median)) +
  geom_point() +
  geom_errorbar(aes(ymin = q5, ymax = q95))

fit$summary("y_rep") %>%
  mutate(date = spy$date, r= spy$pct_return) %>%
  ggplot(aes(date, median)) +
  geom_point() +
  geom_errorbar(aes(ymin = q5, ymax = q95)) +
  geom_point(aes(date, r), color = "red")


# insample pvals
pvals <-
  fit$draws("y_rep") %>%
  as_draws_df() %>%
  as_tibble() %>%
  select(-.chain, -.draw, -.iteration) %>%
  as.matrix() %>%
  t() %>%
  as_tibble() %>%
  mutate(date = spy$date, r = spy$pct_return) %>%
  pivot_longer(c(-date, -r)) %>%
  group_by(date) %>%
  summarize(p = mean(value <= r)) %>%
  mutate(z = qnorm(p))

pvals %>% ggplot(aes(date, z)) + geom_point() + geom_line()
pvals$z %>% qqnorm()
pvals$z %>% qqline()

par_draws <- fit$draws(c("sigma", "sigmaq", "log_lik[10]", "log_lik[50]", "log_lik[100]", "log_lik[126]"))
mcmc_pairs(par_draws)









### sim AR1

sim_ar1 <- function(N, phi, sigma) {
  v_raw <- rnorm(N)
  v <- sigma * v_raw
  h <- rep(0.0, N)
  
  h[1] = v[1] / sqrt(1 - phi^2);
  
  for (t in 2:N) {
    h[t] = phi * h[t-1] + v[t];
  }
  
  tibble(i = 1:N, v_raw, v, h)
}

sim_ar1(100, 0.8, 0.1) %>% ggplot(aes(i, h)) + geom_point() + geom_line()

### sim ARG
sim_arg <- function(N, nu, C, phi, d2) {
  h <- rep(0.0, N)
  z <- rep(0L, N)
  r <- rep(0L, N)
  h[1] <- rgamma(1, nu, scale = C / (1-phi))
  
  for (t in 2:N) {
    r[t] <- d2 + phi*h[t-1] / C
    z[t] <- rpois(1, r[t])
    # z[t] <- rgamma(1, r[t], 1)
    h[t] <- rgamma(1, nu + z[t], scale = C)
  }
  
  tibble(i = 1:N, z, r, h)
}

df <- sim_arg(1000, nu = 0.0, C = 0.1, phi = 0.99, d2 = 0.05)
df %>% ggplot(aes(i, sqrt(h))) + geom_point() + geom_line()

N <- 1e3; h0 <- 0.5; C <- 0.2; phi <- 0.9; d2 <- 0.2; tibble(z = rpois(N, d2 + phi*h0 / C), h = rgamma(N, z, 1/C)) %>% pull(h) %>% qplot()
N <- 1e3; h0 <- 0.5; phi <- 0.9; c2 <- 1^2; tau <- 1e0; tibble(l2 = rcauchy(N)^2, lt = sqrt(c2*l2 / (c2 + tau^2*l2)), s = tau*lt) %>% pull(s) %>% qplot()


### sim nonaffine AR1

sim_ar1na <- function(N, phi, sigma) {
  v_raw <- rnorm(N)
  v <- sigma * v_raw
  h <- rep(0.0, N)
  
  h[1] = v[1] / sqrt(1 - phi^2);
  
  for (t in 2:N) {
    h[t] = phi * h[t-1] + exp(h[t-1]) * v[t];
  }
  
  tibble(i = 1:N, v_raw, v, h)
}

sim_ar1na(100, 0.9, 0.3) %>% ggplot(aes(i, exp(h/2 - 1.3))) + geom_point() + geom_line()

### sim horseshoe nonaffine
rrhs <- function(N, tau, C) {
  l2 <- rcauchy(N)^2
  tau * sqrt(C^2*l2 / (C^2 + tau^2*l2))
}

sim_ar1hs <- function(N, phi, tau, C) {
  s <- rep(1.0, N)
  h <- rep(0.0, N)

  h[1] <- 0

  for (t in 2:N) {
    s[t] <- rrhs(1, tau + phi*h[t-1], C)
    h[t] <- abs(rnorm(1, 0, s[t]))
    # h[t] <- rnorm(1, 0, s[t])
    # h[t] <- rsn(1, 0, s[t], alpha = 1.0)
    # h[t] <- s[t]
  }
  
  tibble(i = 1:N, s, h)
}

df <- sim_ar1hs(230, phi = 2.5, tau = 1e-2, C = 10.0)
df %>% ggplot(aes(i, h)) + geom_point() + geom_line()
df %>% ggplot(aes(i, sqrt(h))) + geom_point() + geom_line()
df %>% ggplot(aes(i, log(h))) + geom_point() + geom_line()
df %>% ggplot(aes(i, sqrt(h)-sqrt(lag(h)))) + geom_point() + geom_line()

g %>% filter(year(date) == 2021) %>% mutate(g = median) %>% select(date, g) %>% mutate(g_rep = df$h) %>% pivot_longer(-date) %>% ggplot(aes(date, value, color = name)) + geom_point() + geom_line()
df %>% mutate(v = log(h) - lag(log(h))) %>% ggplot(aes(i, v)) + geom_point() + geom_line()

### sim horseshoe2
rhs <- function(l, tau, C) tau * sqrt(C^2*l^2 / (C^2 + tau^2*l^2))

sim_ar1hs2 <- function(N, tau, C, theta, rho, phi) {
  l_raw <- abs(rcauchy(N))
  v_raw <- rnorm(N) + rexp(N)
  l <- rep(0.0, N)
  v <- rep(0.0, N)
  h <- rep(0.0, N)
  
  h[1] <- 0
  
  for (t in 2:N) {
    l[t] <- rrhs(1, tau + theta*l[t-1], C)
    v[t] <- l[t] * v_raw[t]
    h[t] <- phi*h[t-1] + v[t]
  }
  
  tibble(i = 1:N, l, v, h)
}

df <- sim_ar1hs2(126, tau = 1e-3, C = 2, theta = 0.7, rho = 1, phi = 0.5)
df %>% ggplot(aes(i, h)) + geom_point() + geom_line()

####
sim_doublear1 <- function(N, phi, rho, theta, tau) {
  u_raw <- rnorm(N)
  u <- tau * u_raw
  v_raw <- rnorm(N)
  g <- rep(0.0, N)
  h <- rep(0.0, N)
  v <- rep(0.0, N)
  
  g[1] = u[1] / sqrt(1 - theta^2);
  h[1] = v[1] / sqrt(1 - phi^2);
  
  for (t in 2:N) {
    g[t] <- theta*g[t-1] + u[t]
    v[t] <- rho*u[t] + exp(g[t]/2) * v_raw[t]
    h[t] = phi * h[t-1] + v[t]
  }
  
  tibble(i = 1:N, u_raw, u, v_raw, v, g, h)
}



df <- sim_doublear1(250, phi = 0.95, rho = 1, theta = 0.7, tau = 0.1)
df  %>% ggplot(aes(i, h)) + geom_point() + geom_line()

### sim nonaffine AR1

sim_ar1 <- function(N, phi, sigma, K, theta) {
  v_raw <- rnorm(N)
  g <- rep(0.0, N)
  v <- rep(0.0, N)
  h <- rep(0.0, N)
  
  h[1] = v_raw[1] / sqrt(1 - phi^2);
  
  for (t in 2:N) {
    g[t] <- sigma * exp(theta*v[t-1]) / (K + exp(theta*v[t-1]))
    v[t] <- exp(g[t]*v_raw[t])*v_raw[t]
    h[t] = phi * h[t-1] + v[t];
  }
  
  tibble(i = 1:N, v_raw, g, v, h)
}

df <- sim_ar1(250, phi = 0.95, sigma = 0.4, K = 0.5, theta = 0.9)
df  %>% ggplot(aes(i, g)) + geom_point() + geom_line()

### sim normal ARG
sim_normarg <- function(N, sigma, phi) {
  v_raw <- rnorm(N)
  v <- rep(0.0, N)
  v_raw[1] <- rnorm(1)
  v[1] <- sigma*v_raw[1]
  
  for (t in 2:N) {
    v_raw[t] <- rnorm(1)
    v[t] <- sigma*rnorm(1, phi*v[t-1], exp(0.1*phi*v[t-1]+0.1*v_raw[t]))
  }
  
  tibble(i = 1:N, v)
}

df <- sim_normarg(252, sigma = 1, phi = 0.99)
df %>% ggplot(aes(i, v)) + geom_point() + geom_line()

              