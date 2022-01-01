library(tidyverse)
library(tidyquant)
library(cmdstanr)

spy <-
  tq_get(c("SPY"),
         get = "stock.prices",
         from = "2018-11-20",
         to = "2021-11-24") %>%
  arrange(date) %>%
  mutate(pct_return = 100 * (close - lag(close)) / lag(close)) %>%
  # if any returns are exactly zero this causes the lp to evaluate to zero
  mutate(pct_return = ifelse(pct_return == 0.0, 1e-4, pct_return)) %>%
  na.omit() 

vix <-
  tq_get(c("^VIX"),
         get = "stock.prices",
         from = "2015-08-20",
         to = "2020-02-21") %>%
  arrange(date) %>%
  mutate(pct_return = 100 * (close - lag(close)) / lag(close)) %>%
  # if any returns are exactly zero this causes the lp to evaluate to zero
  mutate(pct_return = ifelse(pct_return == 0.0, 1e-4, pct_return)) %>%
  mutate(log_r = 100 * log(close/lag(close))) %>%
  na.omit() 

model <- cmdstan_model("stan/sv.stan")
model_vix <- cmdstan_model("stan/sv_joint_vix.stan")

model <- cmdstan_model("stan/sv_wang.stan")

dat <- list(T = nrow(spy), y = spy$pct_return, x = vix$log_r)

dat <- list(T = nrow(spy), y = spy$pct_return)

fit <-
  model$sample(
    data = dat,
    seed = 1994,
    iter_warmup = 1e3,
    iter_sampling = 1e4,
    chains = 4,
    parallel_chains = 4,
    refresh = 1e2,
    max_treedepth = 10,
    adapt_delta = 0.8,
    show_messages = FALSE
  )

log_s <- fit$summary("log_s")

s <- fit$summary("lambda_h") %>% mutate(date = spy$date)
s %>% ggplot(aes(date, median)) + geom_line() + geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.2) + scale_y_log10()

v <- fit$draws("v") %>% as_draws_matrix()

log_s %>%
  bind_cols(vix) %>%
  ggplot(aes(median, log(adjusted)/2 -1.4)) +
  geom_point() +
  # geom_errorbar(aes(ymin = q5, ymax = q95))
  stat_function(fun=function(x) 0 + 0.1*x - 0.00*x^2, color = "orange")

fit_vix <-
  model_vix$sample(
    data = dat,
    seed = 1994,
    iter_warmup = 1e3,
    iter_sampling = 1e3,
    chains = 4,
    parallel_chains = 4,
    refresh = 2e1,
    max_treedepth = 10,
    adapt_delta = 0.8,
    show_messages = FALSE,
    init = replicate(4, list(mu = -2, a2 = 0.001, a3 = 0.0001), simplify = FALSE)
  )

s_vix <- fit_vix$summary("v") %>% mutate(date = vix$date) %>% mutate(model = "vix")

s_vix %>%
  ggplot(aes(date, median, color = model)) +
  geom_pointrange(aes(ymin = q5, ymax = q95), position = position_dodge(width = 0.5))

bind_rows(s, s_vix) %>%
  ggplot(aes(date, median, color = model)) +
  geom_pointrange(aes(ymin = q5, ymax = q95), position = position_dodge(width = 0.5))

# Features to add
# 1. RHS on v
# 2. skew on v via v = exp(a*v_raw)*v_raw
# 3. AR on v via tau
# 4. quadratic rho
# 5. long term moving mu
# 6. volatility premium




