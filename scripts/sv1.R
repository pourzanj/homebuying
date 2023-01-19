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

vix <-
  tq_get(c("^VIX"),
         get = "stock.prices",
         from = "2002-01-01",
         to = "2022-12-31") %>%
  select(date, vix = close)
  
#############
# Initialization fit
#############

model <- cmdstan_model("stan/sv1.stan")
model <- cmdstan_model("stan/sv13.stan")
model <- cmdstan_model("stan/sv18.stan")

init <-
  replicate(4, list(mu_ret = 0.1,
                    phi_raw = c(0.99, 0.9), sigma = c(0.1,0.3),
                    v1 = rep(0.0, nrow(spx)),
                    v2 = rep(0.0, nrow(spx))), simplify = FALSE)

fit <-
  model$sample(
    data = list(T = nrow(spx), y = spx$y,
                M = 2, y_ahead = c(0, 0)),
    seed = 1994,
    iter_warmup = 5e2,
    iter_sampling = 5e2,
    chains = 4,
    parallel_chains = 4,
    refresh = 10,
    max_treedepth = 10,
    adapt_delta = 0.8,
    init = init
  )

fit$summary(c("lp__", "mu_ret", "mu", "phi", "sigma"))

#############
# Fit full model using initialization
#############

# model <- cmdstan_model("stan/sv11.stan")
# 
# # Set intial values
# pars0 <-
#   fit$draws(c("mu_ret", "rho", "mu", "phi", "sigma")) %>%
#   as_draws_df() %>%
#   as_tibble() %>%
#   filter(.draw == 1821)
# 
# v0 <- as.vector(fit$draws(c("v"), format = "draws_matrix")[1821,])
# 
# init <-
#   replicate(4, list(mu_ret = pars0$mu_ret,
#                     mu = pars0$mu, phi = pars0$phi, sigma = pars0$sigma,
#                     rho = pars0$rho,
#                     b = c(0, 0, 0),
#                     v = v0), simplify = FALSE)
# 
# saveRDS(init, "init.RDS")
# init <- readRDS(init, "init.RDS")
# 
# # Set inverse metric
# inv_metric_orig <- diag(fit$inv_metric()[[1]])
# inv_metric_new <- c(inv_metric_orig, rep(0.2^2, 3))
# 
# saveRDS(inv_metric_new, "inv_metric.RDS")
# inv_metric_new <- readRDS("inv_metric.RDS")
# 
# 
# # Fit
# fit <-
#   model$sample(
#     data = list(T = nrow(spx), y = spx$y, K = 3, knots = c(-1.0, -0.5, 0)),
#     seed = 1994,
#     iter_warmup = 0e3,
#     iter_sampling = 1e3,
#     chains = 4,
#     parallel_chains = 4,
#     refresh = 20,
#     max_treedepth = 6,
#     adapt_delta = 0.8,
#     metric = "diag_e",
#     inv_metric = inv_metric_new,
#     adapt_engaged = FALSE,
#     step_size = 0.02,
#     init = init
#   )

#############
# Analyze
#############

fit$summary(c("sum_log_lik", "sum_log_lik_ahead"))
fit$summary(c("lp__", "mu_ret", "mu", "sigma", "phi", "l", "p", "gamma"))

s <- fit$summary("s") %>% bind_cols(spx)
s %>% ggplot(aes(date, median)) + geom_line() + geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.2) + scale_y_log10() 

h <- fit$summary("h") %>% bind_cols(spx)
h %>% ggplot(aes(date, median)) + geom_line() + geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.2)  


v <- fit$summary("v") %>% bind_cols(spx)
v %>% filter(year(date) == 2022) %>% ggplot(aes(date, median)) + geom_pointrange(aes(ymin = q5, ymax = q95))

fit$draws(c("mu_ret", "rho", "mu", "phi", "sigma", "sum_log_lik", "sum_log_lik_ahead")) %>% mcmc_pairs()

y_rep <- fit$summary("y_rep") %>% bind_cols(spx)

y_rep_draws <- fit$draws("y_rep") %>% as_draws_matrix()
p <- map_dbl(1:nrow(spx), function(i) mean(spx$y[i] >= y_rep_draws[,i]))

# Redraw new y_rep
temp <-
  fit$draws(c("mu_ret", "rho", "mu", "phi", "sigma", "alpha", "beta", "h")) %>%
  as_draws_df() %>%
  as_tibble() %>%
  select(-.chain, -.iteration) %>%
  sample_n(1) %>%
  pivot_longer(c(-.draw,-mu_ret, -rho, -mu, -phi, -sigma, -alpha, -beta)) %>%
  mutate(date = spx$date, y = spx$y, s = s$median)

temp %>%
  rename(h = value) %>%
  # mutate(mul = exp(0.4*lag(h)+-0.3*lag(h)^2)) %>%
  mutate(mul = exp(alpha*lag(h)+beta*pmax(-0.7, lag(h)))) %>%
  mutate(v_rep = rnorm(nrow(.))) %>%
  mutate(v_rep1 = mul*rnorm(nrow(.))) %>%
  mutate(h_rep = phi*lag(h) + sigma*v_rep) %>%
  mutate(h_rep1 = phi*lag(h) + sigma*v_rep1) %>%
  mutate(s_rep = exp((mu+h_rep)/2)) %>%
  mutate(s_rep1 = exp((mu+h_rep1)/2)) %>%
  mutate(eps = rnorm(nrow(.))) %>%
  mutate(y0 = mu_ret + rho*v_rep*s_rep + sqrt(1-rho^2)*s_rep*eps) %>%
  mutate(y1 = mu_ret + rho*v_rep1*s_rep1 + sqrt(1-rho^2)*s_rep1*eps) %>%
  select(date, s, y, y0, y1) %>%
  mutate(sbin = lag(cut(s, quantile(s, c(0, 0.25, 0.5, 0.75, 1))))) %>%
  na.omit() %>%
  pivot_longer(c(-date, -s, -sbin)) %>%
  ggplot(aes(value, color = name)) +
  geom_density() +
  facet_wrap(sbin ~ ., scales = "free")


s %>%
  select(date, s = median, y) %>%
  mutate(sbin = lag(cut(s, quantile(s, seq(0, 1, by = 0.25), na.rm = TRUE)))) %>%
  mutate(y1 = as.vector(y_rep_draws[sample(1:dim(y_rep_draws)[1], 1),])) %>%
  # mutate(y2 = as.vector(y_rep_draws[sample(1:dim(y_rep_draws)[1], 1),])) %>%
  # mutate(y3 = as.vector(y_rep_draws[sample(1:dim(y_rep_draws)[1], 1),])) %>%
  filter(row_number() > 1) %>%
  na.omit() %>%
  # mutate(y0 = temp$y0) %>%
  # filter(year(date) == 2022) %>%
  pivot_longer(c(y, y1)) %>%
  ggplot(aes(value, color = name)) +
  geom_density() +
  # geom_histogram(binwidth = 0.25) +
  facet_wrap(sbin ~ ., scales = "free") +
  scale_y_sqrt()

s %>%
  select(date, s = median, y) %>%
  mutate(sbin = lag(cut(s, quantile(s, seq(0, 1, by = 0.25), na.rm = TRUE)))) %>%
  mutate(p = p) %>%
  filter(row_number() > 1) %>%
  na.omit() %>%
  ggplot(aes(p)) +
  geom_histogram(breaks = seq(0, 1, by = 0.1)) +
  facet_wrap(sbin ~ .)

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
  ggplot(aes(date, value, color = name)) +
  geom_line() +
  facet_grid(.draw ~ name)

hh_rep %>%
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

# Get moments of v
v_draws <- fit$draws("v_raw") %>% as_draws_matrix()

moments <-
  v_draws %>%
  as.matrix() %>%
  t() %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(date = spx$date, s = h$median) %>%
  mutate(sbin = lag(cut(s, quantile(s, seq(0, 1, by = 0.1))))) %>%
  mutate(s = lag(s)) %>%
  na.omit() %>%
  pivot_longer(c(-date, -s, -sbin)) %>%
  group_by(sbin, name) %>%
  summarize(m1 = moment(value, order = 1, central = FALSE),
            m2 = moment(value, order = 2, central = TRUE),
            m3 = moment(value, order = 3, central = TRUE),
            m4 = moment(value, order = 4, central = TRUE))

normal_moments <-
  tibble(d = 1:(dim(v_draws)[2])) %>%
  crossing(sbin = unique(moments$sbin)) %>%
  mutate(value = rnorm(nrow(.))) %>%
  group_by(sbin) %>%
  summarize(m1 = moment(value, order = 1, central = FALSE),
            m2 = moment(value, order = 2, central = TRUE),
            m3 = moment(value, order = 3, central = TRUE),
            m4 = moment(value, order = 4, central = TRUE)) %>%
  pivot_longer(-sbin, names_to = "moment")

moments %>%
  ungroup() %>%
  pivot_longer(c(-sbin, -name), names_to = "moment") %>%
  ggplot(aes(value)) +
  geom_histogram() +
  facet_grid(sbin ~ moment, scales = "free") +
  geom_vline(aes(xintercept = value), color = "red", data = normal_moments)

# Plot moments of y and v versus rep
hvy_samples <-
  fit$draws(c("h", "v", "y_rep")) %>%
  as_draws_df() %>%
  as_tibble() %>%
  select(-.chain, -.iteration) %>%
  pivot_longer(-.draw) %>%
  mutate(var = str_extract(name, "[a-z]+")) %>%
  mutate(idx = parse_number(name)) %>%
  select(-name) %>%
  pivot_wider(names_from = var, values_from = value)

hvy_samples %>%
  rename(y_rep = y) %>%
  inner_join(transmute(spx, date, idx = row_number(), y)) %>%
  filter(idx > 1) %>%
  arrange(.draw, date) %>%
  group_by(.draw) %>%
  mutate(hbin = cut(h, quantile(h, seq(0, 1, by = 0.25)), labels = 1:4)) %>%
  mutate(laghbin = lag(hbin)) %>%
  na.omit() %>%
  ungroup() %>%
  group_by(.draw, laghbin) %>%
  summarize(m1_true = moment(y, order = 1, central = FALSE),
            m2_true = moment(y, order = 2, central = TRUE),
            m3_true = moment(y, order = 3, central = TRUE),
            m4_true = moment(y, order = 4, central = TRUE),
            m1_rep = moment(y_rep, order = 1, central = FALSE),
            m2_rep = moment(y_rep, order = 2, central = TRUE),
            m3_rep = moment(y_rep, order = 3, central = TRUE),
            m4_rep = moment(y_rep, order = 4, central = TRUE)) %>%
  ungroup() %>%
  pivot_longer(c(-.draw, -laghbin), names_sep = "_", names_to = c("moment", "rep")) %>%
  ggplot(aes(value, color = rep)) +
  # geom_histogram() +
  geom_density() +
  # facet_grid(laghbin ~ moment, scales = "free") +
  facet_wrap(laghbin + moment ~ ., scales = "free")

hvy_samples %>%
  rename(y_rep = y) %>%
  inner_join(transmute(spx, date, idx = row_number(), y)) %>%
  mutate(v_rep = rnorm(nrow(.))) %>%
  filter(idx > 1) %>%
  arrange(.draw, date) %>%
  group_by(.draw) %>%
  mutate(hbin = cut(h, quantile(h, seq(0, 1, by = 0.25)), labels = 1:4)) %>%
  mutate(laghbin = lag(hbin)) %>%
  na.omit() %>%
  ungroup() %>%
  group_by(.draw, laghbin) %>%
  summarize(m1_true = moment(v, order = 1, central = FALSE),
            m2_true = moment(v, order = 2, central = TRUE),
            m3_true = moment(v, order = 3, central = TRUE),
            m4_true = moment(v, order = 4, central = TRUE),
            m1_rep = moment(v_rep, order = 1, central = FALSE),
            m2_rep = moment(v_rep, order = 2, central = TRUE),
            m3_rep = moment(v_rep, order = 3, central = TRUE),
            m4_rep = moment(v_rep, order = 4, central = TRUE)) %>%
  ungroup() %>%
  pivot_longer(c(-.draw, -laghbin), names_sep = "_", names_to = c("moment", "rep")) %>%
  ggplot(aes(value, color = rep)) +
  # geom_histogram() +
  geom_density() +
  # facet_grid(laghbin ~ moment, scales = "free") +
  facet_wrap(laghbin + moment ~ ., scales = "free")

# Generate horsehose draws
N <- 1e6; c2 <- 1; tau <- 3; nu <- 5;
tibble(l2 = rcauchy(N)^2) %>%
  mutate(lt = sqrt((c2*l2) / (c2+tau^2*l2))) %>%
  mutate(v1 = rt(N, df = nu)*tau*lt) %>%
  # mutate(v1 = rnorm(N)*tau*lt) %>%
  mutate(v2 = rnorm(N)) %>%
  select(v1, v2) %>%
  mutate(i = row_number()) %>%
  pivot_longer(-i) %>%
  ggplot(aes(value, color = name)) +
  # geom_histogram() +
  geom_density() +
  xlim(-5, 5)
  e# facet_grid(name ~ .)
  # scale_y_sqrt()

N <- 1e5
tibble(i = 1:N,
       x = rnorm(N)) %>%
  mutate(y = 0.3*pmin(0, x - -0.67) + x) %>%
  pivot_longer(-i) %>%
  ggplot(aes(value)) +
  geom_histogram(binwidth = 0.1) +
  facet_grid(name ~ .)

###### Play with weeks LOO
my.loglikelihood <- function(z, mu, logsigma, logdf) {
  n     <- length(z)
  df    <- exp(logdf)
  sigma <- exp(logsigma)
  LL    <- -(n/2)*log(pi) + n*lgamma((df+1)/2) - n*lgamma(df/2) - (n/2)*log(df) -
    n*log(sigma) - ((df+1)/2)*sum(log(1 + (1/df)*((z-mu)/sigma)^2))
  LL
}

my.MLE <- function(z, df) {
  NEGLOGLIKE <- function(par) { -my.loglikelihood(z, par[1], par[2], log(df)) }
  PAR0   <- c(mean(z), log(sd(z)))
  OPTIM  <- optim(fn = NEGLOGLIKE, par = PAR0)
  PARHAT <- OPTIM$par
  MLE    <- c(PARHAT[1], exp(PARHAT[2]))
  MLE
}

my.MLE(x, df = 1)

sd5 <- rollify(sd, 5)
mabs5 <- rollify(function(x) mean(abs(x)), 5)

spx %>%
  select(year, week, day, y) %>%
  crossing(holdout = c("Mon", "Tue", "Wed", "Thu", "Fri")) %>%
  arrange(year, week, holdout, day) %>%
  group_by(year, week, holdout) %>%
  filter(any(day == holdout)) %>%
  filter(day != holdout) %>%
  # summarize(n = n(), m = mean(y), s = sd(y)) %>%
  summarize(n = n(), m = median(y), b = mean(abs(y-m))) %>%
  # summarize(n = n(), m = my.MLE(y, 1)[1], s = my.MLE(y, 1)[2]) %>%
  ungroup() %>%
  inner_join(select(spx, year, week, holdout = day, y)) %>%
  # mutate(p = pnorm(y, m, s)) %>%
  mutate(p = plaplace(y, m, b)) %>%
  # mutate(p = pt2(y, m, s, 1)) %>%
  pull(p) %>% qplot(binwidth = 0.01)

spx %>%
  mutate(s = lag(sd5(y)), b = lag(mabs5(y))) %>%
  na.omit() %>%
  mutate(x = rnorm(nrow(.), 0.1, s),
         xl = rlaplace(nrow(.), 0.1, b)) %>%
  mutate(sbin = cut(s, quantile(s, c(0.0, 0.25, 0.5, 0.75, 1.0)))) %>%
  mutate(sbin = lag(sbin)) %>%
  na.omit() %>%
  select(date, x, y, xl, sbin) %>%
  pivot_longer(c(x, y, xl)) %>%
  ggplot(aes(value, color = name)) +
  geom_histogram() +
  facet_grid(name ~ sbin, scales = "free") +
  scale_y_sqrt()

y <-
  spx %>%
  mutate(s = lag(sd5(y)), b = lag(mabs5(y))) %>%
  na.omit() %>%
  mutate(x = rnorm(nrow(.), 0.1, s),
         xl = rlaplace(nrow(.), 0.1, b)) %>%
  mutate(sbin = cut(s, quantile(s, c(0.0, 0.25, 0.5, 0.75, 1.0)))) %>%
  mutate(sbin = lag(sbin)) %>%
  na.omit() %>%
  select(date, x, y, xl, sbin) %>%
  # filter(sbin == "(0.053,0.517]") %>%
  # filter(sbin == "(1.25,9.59]") %>%
  group_split(sbin)

model <- cmdstan_model("stan/normal_mixture.stan")

init <-
  replicate(4, list(theta = c(0.5, 0.4, 0.1),
                    mu = c(0.1, 0.0, 0.0),
                    sigma = c(0.1, 0.5, 2)), simplify = FALSE)

fit_wrapper <- function(i) {
  model$sample(
    data = list(K = 3, N = nrow(y[[i]]), y = y[[i]]$y),
    seed = 1994,
    iter_warmup = 5e2,
    iter_sampling = 5e2,
    parallel_chains = 4,
    init = init
  )
}

fits %>%
  map(function(fit) fit$summary("sigma")) %>%
  map2_dfr(1:length(.), function(df, i) mutate(df, i = i)) %>%
  ggplot(aes(i, median)) +
  geom_pointrange(aes(ymin = q5, ymax = q95)) +
  facet_grid(variable ~ ., scales = "free")

fits %>%
  map(function(fit) fit$summary("sigma")) %>%
  map2_dfr(1:length(.), function(df, i) mutate(df, i = i)) %>%
  ggplot(aes(i, median, color = variable)) +
  geom_pointrange(aes(ymin = q5, ymax = q95)) +
  scale_y_log10() +
  geom_smooth(method = "lm", se = FALSE)

y_rep <-
  fit$draws("y_rep") %>%
  as_draws_df() %>%
  as_tibble() %>%
  sample_n(3) %>%
  select(-.chain, -.iteration) %>%
  pivot_longer(-.draw) %>%
  bind_rows(tibble(.draw = -1, name = "y", value = y))

y_rep %>%
  ggplot(aes(value)) +
  geom_histogram() +
  facet_wrap(.draw ~ .)

####################
model <- cmdstan_model("stan/sv5.stan")

K <- 3
init <-
  replicate(4, list(phi = 0.98,
                    sigma = 0.3,
                    a_theta = rep(0.0, K), b_theta = rep(0.0, K),
                    a_m = rep(0.0, K), b_m = rep(0.0, K),
                    a_s = 1:K,
                    v = rep(0.0, nrow(spx))), simplify = FALSE)

init <-
  replicate(4, list(phi = 0.98,
                    sigma = 0.3,
                    mu = c(0, 0, 0),
                    theta = c(0.45, 0.45, 0.1),
                    a_s = 1:K,
                    v = rep(0.0, nrow(spx))), simplify = FALSE)

fit <-
  model$sample(
    data = list(K = 3, T = nrow(spx), y = spx$y),
    seed = 1994,
    iter_warmup = 5e2,
    iter_sampling = 5e2,
    chains = 4,
    parallel_chains = 4,
    refresh = 10,
    max_treedepth = 5,
    adapt_delta = 0.5,
    init = init
  )

fit$summary() %>% print(n = 15)

fit$summary("h") %>%
  bind_cols(spx) %>%
  ggplot(aes(date, median)) +
  geom_line() +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.2)

s <- fit$summary("s")

s %>%
  mutate(i = as.integer(str_extract(variable, "(?<=\\[)[0-9]+"))) %>%
  mutate(k = as.integer(str_extract(variable, "(?<=,)[0-9]+"))) %>%
  select(i, k, median, q5, q95) %>%
  inner_join(transmute(spx, i = row_number(), date)) %>%
  ggplot(aes(date, median, color = factor(k))) +
  geom_line() +
  scale_y_log10()






