library(tidyverse)
library(tibbletime)
library(tidyquant)
library(cmdstanr)
library(posterior)
library(bayesplot)

spx <-
  tq_get(c("^GSPC"),
         get = "stock.prices",
         from = "2000-01-01",
         to = "2021-12-27") %>%
  arrange(date) %>%
  mutate(y0 = 0.0) %>%
  mutate(y = 100 * log(close/lag(close))) %>%
  mutate(a = ifelse(lag(close) < low, 0.0, 100 * log(low/lag(close)))) %>%
  mutate(b = ifelse(lag(close) > high, 0.0, 100 * log(high/lag(close)))) %>%
  # if any returns are exactly zero this causes the lp to evaluate to zero
  na.omit()

model <- cmdstan_model("stan/molina.stan")
model <- cmdstan_model("stan/molina_sv.stan")

h0 <- 2*log(spx$b - spx$a)
v0 <- h0-lag(h0)

fit <-
  model$sample(
    data = list(T = nrow(spx), y0 = spx$y0, y = spx$y, a = spx$a, b = spx$b),
    seed = 1994,
    iter_warmup = 2e2,
    iter_sampling = 2e2,
    chains = 4,
    parallel_chains = 4,
    refresh = 1e0,
    # adapt_engaged = FALSE,
    # metric = "unit_e",
    # step_size = 1e-2,
    max_treedepth = 8,
    adapt_delta = 0.8,
    show_messages = FALSE,
    init = replicate(4, list(phi = 0.67, s0 = 0.73, sigma = 0.5, theta = 0.5, tau = 0.1,
                             rho = -0.6, gamma = 0.0, rhog = 0.0,
                             h = 2*log(spx$b - spx$a), u_raw = rep(0.0, nrow(spx)), nu = 10, U = rep(1.0, nrow(spx))), simplify = FALSE)
  )

# m <- fit$summary("m") %>% mutate(date = spx$date)
s <- fit$summary("s") %>% mutate(date = spx$date)
h <- fit$summary("h") %>% mutate(date = spx$date, y = spx$y)
v_raw <- fit$summary("v_raw") %>% mutate(date = spx$date)
eps <- fit$summary("eps") %>% mutate(date = spx$date)
y_rep <- fit$summary("y_rep") %>% mutate(date = spx$date)
s_rep <- fit$summary("s_rep") %>% mutate(date = spx$date)


# m %>% ggplot(aes(date, median)) + geom_point() + geom_line() + geom_errorbar(aes(ymin = q5, ymax = q95))
s %>% ggplot(aes(date, median)) + geom_point() + geom_line() + geom_errorbar(aes(ymin = q5, ymax = q95))
h %>% ggplot(aes(date, median)) + geom_point() + geom_line() + geom_errorbar(aes(ymin = q5, ymax = q95))
v_raw %>% ggplot(aes(date, median)) + geom_point() + geom_line() + geom_errorbar(aes(ymin = q5, ymax = q95))
eps %>% ggplot(aes(date, median)) + geom_point() + geom_line() + geom_errorbar(aes(ymin = q5, ymax = q95))

y_rep %>%
  mutate(y = spx$y) %>%
  mutate(o5 = y <= q5, o95 = y >= q95, o = o5 | o95) %>%
  filter(year(date) == 2021) %>%
  ggplot(aes(date, median)) +
  geom_line() +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.4) +
  geom_point(aes(date, y, color = o))

y_rep %>%
  mutate(y = spx$y) %>%
  mutate(o5 = y <= q5, o95 = y >= q95) %>%
  group_by(year(date)) %>%
  summarize(n = n(), o5 = sum(o5), o95 = sum(o95)) %>%
  mutate(p5 = pbinom(o5, n, 0.05), p95 = pbinom(o95, n, 0.05))

h <- fit$summary("h") %>% mutate(date = spx$date)
h_rep <- fit$summary("h_rep") %>% mutate(date = spx$date) %>% bind_cols(select(h, value = median, v5 = q5, v95 = q95)) %>% mutate(name = "h")

v_raw <- fit$summary("v_raw") %>% mutate(date = spx$date)
v_raw_rep <- fit$summary("v_raw_rep") %>% mutate(date = spx$date) %>% bind_cols(select(v_raw, value = median, v5 = q5, v95 = q95)) %>% mutate(name = "v_raw")

y_rep <- fit$summary("y_rep") %>% mutate(date = spx$date, value = spx$y) %>% mutate(name = "y")

y_rep %>%
  bind_rows(h_rep) %>%
  bind_rows(v_raw_rep) %>%
  mutate(o5 = value <= q5, o95 = value >= q95, o = o5 | o95) %>%
  filter(year(date) == 2021, month(date) <= 8) %>%
  ggplot(aes(date, median)) +
  geom_line() +
  geom_errorbar(aes(ymin = q5, ymax = q95)) +
  geom_linerange(aes(ymin = v5, ymax = v95, color = o)) +
  geom_point(aes(date, value, color = o)) +
  facet_grid(name ~ ., scales = "free")




# Get individual h for each day with no prior and see how it predicts next day
# return.
hb <-
  h %>%
  mutate(hb = cut(median, quantile(median, seq(0, 1, by = 0.2)), include.lowest = TRUE)) %>%
  mutate(hb = lag(hb))

hb %>%
  ggplot(aes(y)) +
  geom_histogram(binwidth = 0.2) +
  facet_wrap(hb ~ ., scales = "free")

model <- cmdstan_model("stan/generative_sgt.stan")
model <- cmdstan_model("stan/ast.stan")

hb1 <- hb %>% filter(hb == "[-5.09,-2.3]")

dat <-
  list(N = nrow(hb1),
       y = hb1$y,
       tau0 = 1e-2, nu =6)

fit <-
  model$sample(
    data = dat,
    seed = 1994,
    iter_warmup = 2e2,
    iter_sampling = 2e2,
    chains = 4,
    parallel_chains = 4,
    refresh = 2e1,
    # adapt_engaged = FALSE,
    # metric = "unit_e",
    # step_size = 1e-2,
    max_treedepth = 8,
    adapt_delta = 0.8,
    show_messages = FALSE,
    init = replicate(4, list(mu = 0, sigma = 0.2, lambda = 0, nu = 10, tau = 1e-3), simplify = FALSE)
  )

# Just do real SGT

model <- cmdstan_model("stan/sgt.stan")

fit_sgt <- function(df) {
  dat <- list(N = nrow(df), y = df$y, m_p = 1, s_p = 1, m_q = 15, s_q = 5)
  
  fit <-
    model$sample(
      data = dat,
      seed = 1994,
      iter_warmup = 2e2,
      iter_sampling = 2e2,
      chains = 4,
      parallel_chains = 4,
      refresh = 1e2,
      # adapt_engaged = FALSE,
      # metric = "unit_e",
      # step_size = 1e-2,
      max_treedepth = 8,
      adapt_delta = 0.8,
      show_messages = FALSE,
      init = replicate(4, list(mu = 0, s = 0.2, l = 0, p = 2, q = 10), simplify = FALSE)
    )
  
  fit
}

fits <-
  hb %>%
  group_split(hb) %>%
  map(fit_sgt)

fits %>%
  map(function(f) f$summary()) %>%
  bind_rows() %>%
  mutate(i = rep(1:5, each = 6)) %>%
  inner_join(hb %>% group_by(hb) %>% summarize(med = median(lag(median), na.rm = TRUE)) %>% mutate(i = row_number())) %>%
  ggplot(aes(med, median)) +
  geom_point() +
  geom_errorbar(aes(ymin = q5, ymax = q95)) +
  facet_wrap(variable ~ ., scales = "free")

draw_sgt_total_ret <- function(mu, s, l, p, q, N = 1106) {
  y <- rsgt(N, mu, s, l, p, q)
  prod(exp(y/100))^(252/N)
}

fits[[5]]$draws(c("mu", "s", "l", "p", "q")) %>%
  as_draws_df() %>%
  as_tibble() %>%
  select(-.chain, -.iteration, -.draw) %>%
  pmap_dbl(draw_sgt_total_ret) %>%
  qplot()


