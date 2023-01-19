library(tidyverse)
library(tibbletime)
library(tidyquant)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(moments)

spx <-
  tq_get(c("^GSPC"),
         get = "stock.prices",
         from = "2002-01-01",
         to = "2022-12-31") %>%
  mutate(y = 100*log(close/lag(close))) %>%
  na.omit() %>%
  mutate(week = isoweek(date)) %>%
  mutate(year = if_else(month(date) == 12 & week == 1, year(date)+1, year(date))) %>%
  mutate(day = wday(date, label = TRUE))

#############
# Initialization fit
#############

model <- cmdstan_model("stan/sv13_1.stan")

init <-
  replicate(4, list(mu = 0.2, phi = 0.98, sigma = 0.17,
                    alpha = c(0.1, 0, 0),
                    gamma = c(0, 0, 0),
                    beta = c(0.7, 0),
                    q = 10,
                    v0 = 0,
                    v = rep(0.0, nrow(spx))), simplify = FALSE)

fit <-
  model$sample(
    data = list(T = nrow(spx), y = spx$y,
                M = 2, y_oos = c(0,0),
                m_p = 2, s_p = 2,
                m_q = 15, s_q = 6),
    seed = 1994,
    iter_warmup = 5e2,
    iter_sampling = 5e2,
    chains = 4,
    parallel_chains = 4,
    refresh = 10,
    max_treedepth = 8,
    adapt_delta = 0.9,
    init = init
  )

#############
# Analyze
#############

fit$summary(c("sum_log_lik", "sum_log_lik_ahead"))
fit$summary(c("lp__", "mu_ret", "alpha", "beta", "mu", "sigma", "phi", "gamma", "nu"))

s <- fit$summary("s") %>% bind_cols(spx)
s %>% ggplot(aes(date, median)) + geom_line() + geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.2) + scale_y_log10() 

h <- fit$summary("h") %>% bind_cols(spx)
h %>% ggplot(aes(date, median)) + geom_line() + geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.2)  

# Simulate y_rep using rsgt since Stan doesn't have it natively
inv_logit <- function(x) exp(x)/(1+exp(x))
y_rep <-
  fit$draws(c("mu_ret", "mu", "phi", "sigma", "h", "alpha", "beta", "gamma", "nu")) %>%
  as_draws_df() %>%
  as_tibble() %>%
  select(-.chain, -.iteration) %>%
  # sample_n(1) %>%
  pivot_longer(c(-.draw,-mu_ret, -mu, -phi, -sigma,
                 -`nu[1]`, -`nu[2]`, -`nu[3]`,
                 -`alpha[1]`, -`alpha[2]`,
                 -`beta[1]`, -`beta[2]`, -`beta[3]`,
                 -`gamma[1]`, -`gamma[2]`, -`gamma[3]`)) %>%
  mutate(i = parse_number(name)) %>%
  inner_join(transmute(spx, date, y, i = row_number())) %>%
  arrange(.draw, date) %>%
  group_by(.draw) %>%
  mutate(h1 = lag(value)) %>%
  ungroup() %>%
  na.omit() %>%
  mutate(h_rep = phi*h1 + sigma*rnorm(nrow(.))) %>%
  mutate(m_rep = mu_ret + `alpha[1]`*h1 + `alpha[2]`*h1^2) %>%
  mutate(p_rep = exp(`beta[1]` + `beta[2]`*h1 + `beta[3]`*h1^2)) %>%
  mutate(q_rep = exp(`nu[1]` + `nu[2]`*h1 + `nu[3]`*h1^2)) %>%
  mutate(l_rep = 2*inv_logit(`gamma[1]` + `gamma[2]`*h1 + `gamma[3]`*h1^2) - 1) %>%
  mutate(s_rep = exp(0.5*(mu + h_rep))) %>%
  mutate(y_rep = rsgt(nrow(.), m_rep, s_rep, l_rep, p_rep, q_rep, mean.cent = TRUE, var.adj = TRUE))

# Compare density of y and y_rep for different posterior draws
y_rep %>%
  filter(.draw %in% sample(1:2000, 1)) %>%
  mutate(h1bin = cut(h1, quantile(h1, seq(0, 1, by = 0.25)), labels = 1:4)) %>%
  na.omit() %>%
  select(date, h1bin, y, y_rep) %>%
  pivot_longer(c(y, y_rep)) %>%
  ggplot(aes(value, color = name)) +
  geom_density() +
  facet_wrap(h1bin ~ ., scales = "free")

# PPC moments of y and y_rep stratified by h-1
y_rep %>%
  group_by(.draw) %>%
  mutate(h1bin = cut(h1, quantile(h1, seq(0, 1, by = 0.25)), labels = 1:4)) %>%
  na.omit() %>%
  ungroup() %>%
  group_by(.draw, h1bin) %>%
  summarize(m1_true = mean(y),
            m2_true = sd(y),
            m3_true = mean((y-mean(y))^3 / sd(y)^3),
            m4_true = mean((y-mean(y))^4 / sd(y)^4),
            m1_rep = mean(y_rep),
            m2_rep = sd(y_rep),
            m3_rep = mean((y_rep-mean(y_rep))^3 / sd(y_rep)^3),
            m4_rep = mean((y_rep-mean(y_rep))^4 / sd(y_rep)^4)) %>%
  ungroup() %>%
  pivot_longer(c(-.draw, -h1bin), names_sep = "_", names_to = c("moment", "rep")) %>%
  ggplot(aes(value, color = rep)) +
  geom_density() +
  facet_wrap(h1bin + moment ~ ., scales = "free")

# Get p-values
pvals <-
  y_rep %>%
  select(.draw, date, y_rep, y) %>%
  arrange(date, .draw) %>%
  group_by(date) %>%
  summarize(p = mean(y > y_rep)) %>%
  arrange(date)
  
pvals %>%
  # inner_join(select(h, date, h = median)) %>%
  # mutate(h1 = lag(h)) %>%
  # na.omit() %>%
  # mutate(h1bin = cut(h1, quantile(h1, seq(0, 1, by = 0.25)), labels = 1:4)) %>%
  # na.omit() %>%
  ggplot(aes(p)) +
  geom_histogram(breaks = seq(0, 1, by = 0.1))

# Compare draws of h and h_rep
hh_rep <-
  fit$draws(c("h", "h_rep")) %>%
  as_draws_df() %>%
  as_tibble() %>%
  # sample_n(4) %>%
  select(-.chain, -.iteration) %>%
  pivot_longer(-.draw) %>%
  mutate(i = parse_number(name)) %>%
  mutate(name = str_extract(name, "[a-z_]+")) %>%
  inner_join(transmute(spx, date, i = row_number()))

hh_rep %>%
  filter(.draw %in% sample(1:2000, 4)) %>%
  ggplot(aes(date, value, color = name)) +
  geom_line() +
  facet_grid(.draw ~ name)

hh_rep %>%
  filter(.draw %in% sample(1:2000, 4)) %>%
  ggplot(aes(value, color = name)) +
  geom_density() +
  # geom_histogram() +
  # facet_grid(name ~ .draw) +
  facet_grid(.draw ~ .)

hh_rep %>%
  group_by(name, .draw) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  summarize(m1_true = moment(h, order = 1, central = FALSE),
            m2_true = moment(h, order = 2, central = TRUE),
            m3_true = moment(h, order = 3, central = TRUE),
            m4_true = moment(h, order = 4, central = TRUE),
            m1_rep = moment(h_rep, order = 1, central = FALSE),
            m2_rep = moment(h_rep, order = 2, central = TRUE),
            m3_rep = moment(h_rep, order = 3, central = TRUE),
            m4_rep = moment(h_rep, order = 4, central = TRUE)) %>%
  ungroup() %>%
  pivot_longer(-.draw, names_sep = "_", names_to = c("moment", "rep")) %>%
  ggplot(aes(value, color = rep)) +
  # geom_histogram() +
  geom_density() +
  # facet_grid(laghbin ~ moment, scales = "free") +
  facet_wrap(moment ~ ., scales = "free")

vv_rep <-
  fit$draws(c("v", "v_rep")) %>%
  as_draws_df() %>%
  as_tibble() %>%
  # sample_n(4) %>%
  select(-.chain, -.iteration) %>%
  pivot_longer(-.draw) %>%
  mutate(i = parse_number(name)) %>%
  mutate(name = str_extract(name, "[a-z_]+")) %>%
  inner_join(transmute(spx, date, i = row_number()))

vv_rep %>%
  group_by(name, .draw) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  summarize(m1_true = moment(v, order = 1, central = FALSE),
            m2_true = moment(v, order = 2, central = TRUE),
            m3_true = moment(v, order = 3, central = TRUE),
            m4_true = moment(v, order = 4, central = TRUE),
            m1_rep = moment(v_rep, order = 1, central = FALSE),
            m2_rep = moment(v_rep, order = 2, central = TRUE),
            m3_rep = moment(v_rep, order = 3, central = TRUE),
            m4_rep = moment(v_rep, order = 4, central = TRUE)) %>%
  ungroup() %>%
  pivot_longer(-.draw, names_sep = "_", names_to = c("moment", "rep")) %>%
  ggplot(aes(value, color = rep)) +
  # geom_histogram() +
  geom_density() +
  # facet_grid(laghbin ~ moment, scales = "free") +
  facet_wrap(moment ~ ., scales = "free")

#############
# Compute optimal leverages
#############

Er <-
  y_rep %>%
  select(date, y_rep) %>%
  mutate(r = exp(y_rep/100)-1) %>%
  crossing(k = c(-3, -1, 0, 1, 3)) %>%
  mutate(kr = k*r) %>%
  group_by(date, k) %>%
  summarize(Er = mean(log(1 + kr))) %>%
  mutate(Er_year = 100*((exp(Er)^252-1))) %>%
  ungroup()

Er %>%
  filter(year(date) >= 2021) %>%
  ggplot(aes(date, Er_year, color = factor(k))) +
  geom_line()

opt <-
  Er %>%
  inner_join(select(spx, date, y)) %>%
  # filter(year(date) >= 2020) %>%
  group_by(date) %>%
  filter(Er == max(Er)) %>%
  ungroup() %>%
  mutate(r_realized = exp(y/100)-1) %>%
  mutate(kr_realized = k*r_realized) 

# Plot returns by optimal leverage
opt %>%
  ggplot(aes(k, kr_realized)) +
  geom_jitter()

# Compute average returns by optimal leverage
opt %>%
  group_by(k) %>%
  summarize(n = n(), avg_ret = prod(1 + kr_realized)^(252/n))

# Compute yearly returns
opt %>%
  mutate(year = year(date)) %>%
  group_by(year) %>%
  summarize(n = n(), tot_ret = prod(1 + kr_realized)) %>%
  print(n = nrow(.))

# Overall yearly average return
opt %>%
  mutate(year = year(date)) %>%
  summarize(n = n(), avg_yearly_ret = prod(1 + kr_realized)^(252/n))

#############
# Try breaking up over time to see if parameters still hold
#############

spx_split <-
  spx %>%
  mutate(yearbin = cut(year, quantile(year, seq(0, 1, by = 0.1)))) %>%
  na.omit() %>%
  group_split(yearbin)

fit_wrapper <- function(i) {
  model$sample(
    data = list(T = nrow(spx_split[[i]]), y = spx_split[[i]]$y,
                m_p = 2, s_p = 2,
                m_q = 15, s_q = 6),
    seed = 1994,
    iter_warmup = 5e2,
    iter_sampling = 5e2,
    chains = 4,
    parallel_chains = 4,
    refresh = 200,
    max_treedepth = 8,
    adapt_delta = 0.8,
    init = init
  )
}
  

fits <- map(1:length(spx_split), fit_wrapper)

summaries <-
  fits %>%
  map2_dfr(1:length(spx_split),
       function(fit, i)
         fit$summary(c("lp__", "mu_ret", "alpha", "beta", "mu", "sigma", "phi", "gamma", "nu")) %>%
                       mutate(i = i)
)

summaries %>%
  ggplot(aes(i, median)) +
  geom_pointrange(aes(ymin = q5, ymax = q95)) +
  facet_wrap(variable ~ ., scales = "free")









