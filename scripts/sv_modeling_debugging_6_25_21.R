kNumCalendarTrainDays <- round((kNumTrainDays / 252) * 365) + 100
kTodaysDate <- Sys.Date()

spy <-
  tq_get(c("^GSPC"),
         get = "stock.prices",
         from = ymd("2019-06-01"),
         to = kTodaysDate) %>%
  arrange(date) %>%
  mutate(pct_return = 100 * (close - lag(close)) / lag(close)) %>%
  # if any returns are exactly zero this causes the lp to evaluate to zero
  mutate(pct_return = ifelse(pct_return == 0.0, 1e-4, pct_return)) %>%
  na.omit()

## Plot pattern of weekly movements
spy %>%
  mutate(day = wday(date), week = epiweek(date), year = year(date)) %>%
  group_by(year, week) %>% mutate(start_date = min(date)) %>%
  ungroup() %>%
  group_by(start_date) %>%
  filter(year(date) == 2021) %>%
  ggplot(aes(date, pct_return, color = factor(day-1))) +
  geom_point() +
  facet_wrap(start_date ~ ., scales = "free_x") + geom_hline(yintercept = 0)

## Plot pattern of monthly movements
spy %>%
  mutate(year = year(date), month = month(date)) %>%
  group_by(year, month) %>%
  mutate(start_date = min(date)) %>%
  ungroup() %>%
  group_by(start_date) %>%
  mutate(cumret = cumprod(1 + pct_return/100)) %>%
  mutate(r = 100*(prod(1 + pct_return/100)-1),
         min_r = 100*(min(cumret)-1),
         max_r = 100*(max(cumret)-1)) %>%
  ungroup() %>%
  group_by(year) %>%
  mutate(cumret = cumprod(1 + pct_return/100)) %>%
  mutate(r_year = 100*(prod(1 + pct_return/100)-1),
         min_r_year = 100*(min(cumret)-1),
         max_r_year = 100*(max(cumret)-1)) %>%
  ggplot(aes(start_date, pct_return)) +
  geom_jitter(width = 5.0) +
  facet_wrap(year ~ ., scales = "free") +
  geom_point(aes(start_date, r,color = r > 1)) +
  geom_errorbar(aes(ymin = min_r, ymax= max_r, color = r > 1))
  # geom_hline(aes(yintercept = r_year)) +
  # geom_hline(aes(yintercept = min_r_year), color = "red") +
  # geom_hline(aes(yintercept = max_r_year), color = "green")

spy %>%
  mutate(year = year(date), month = month(date)) %>%
  group_by(year, month) %>%
  mutate(start_date = min(date)) %>%
  ungroup() %>%
  group_by(start_date) %>%
  mutate(cumret = cumprod(1 + pct_return/100)) %>%
  summarize(r = 100*(prod(1 + pct_return/100)-1),
         min_r = 100*(min(cumret)-1),
         max_r = 100*(max(cumret)-1)) %>%
  mutate(year = year(start_date)) %>%
  ggplot(aes(r)) + 
  geom_histogram() +
  facet_wrap(year ~ .)

## Does volatility ever go up from a positive day?
sd5 <- rollify(sd, window = 5)
spy %>%
  mutate(s5 = lag(sd5(pct_return))) %>%
  mutate(s_next5 = lead(s5, 6)) %>%
  mutate(z = pct_return / s5) %>%
  mutate(year = year(date)) %>%
  mutate(week = epiweek(date)) %>%
  ggplot(aes(date, z)) +
  geom_point() + 
  facet_wrap(year ~ ., scales = "free_x")

## Make fake leverage data
rho <- -0.7
mu <- 0.45
phi <- 0.9
sigma <- 0.36

N <- nrow(spy)
v <- rnorm(N)
h <- v
h[1] <- v[1] / sqrt(1 - phi^2)
for(t in 2:N) h[t] <- phi * h[t-1] + sigma * v[t]
s <- exp((mu + h) / 2)
y <- rho * s * v + sqrt(1-rho^2) * s * rnorm(N)

sim_dat <- tibble(date = spy$date, v, h, s, y)

sim_dat %>% 
  pivot_longer(-date) %>%
  ggplot(aes(date, value)) +
  geom_line() +
  facet_grid(name ~ ., scales = "free")
  


## Fit Stan
dat <- list(T = nrow(spy), y = spy$pct_return)
dat <- list(T = nrow(sim_dat), y = sim_dat$y)

model <- cmdstan_model("stan/sv_pos_hs.stan")

fit <-
  model$sample(
    data = dat,
    seed = 1994,
    iter_warmup = 1e3,
    iter_sampling = 1e3,
    chains = 1,
    parallel_chains = 1,
    refresh = 1e1,
    max_treedepth = 10,
    adapt_delta = 0.8,
    show_messages = FALSE
  )

phi <- rep(0.95, N)
for(t in 2:N) phi[t] <- rbeta(1, phi[t-1] * 10, (1-phi[t-1])*10)
qplot(1:N, phi)
fit$summary("phi") %>% mutate(time_id = as.integer(str_extract(variable, "[0-9]+"))) %>% ggplot(aes(time_id, median)) + geom_point() + geom_errorbar(aes(ymin = q5, ymax = q95))

draws <- fit$draws()

mcmc_pairs(draws, pars = c("lp__", "rho", "phi", "sigma"))
mcmc_pairs(draws, pars = c("rho", "v[10]","v[50]", "v[100]", "v[200]"))


y_rep <-
  fit$draws("y_rep") %>%
  as_draws_df() %>%
  select(-.chain, -.iteration) %>%
  pivot_longer(-.draw) %>%
  mutate(time_id = as.integer(str_extract(name, "[0-9]+"))) %>%
  select(-name)

y_true <-
  spy %>%
  transmute(.draw = NA_integer_, value = pct_return, time_id = row_number())

y_rep %>%
  filter(.draw %in% c(1, 101, 201, 301, 401, 501, 601, 701, 801, 901, 1000)) %>%
  bind_rows(y_true) %>%
  ggplot(aes(value, color = is.na(.draw))) +
  geom_histogram() +
  facet_wrap(.draw ~ .) +
  scale_y_log10()

y_rep %>%
  bind_rows(y_true) %>%
  group_by(.draw) %>%
  summarize(f = mean(value < -10)) %>%
  arrange(f) %>%
  mutate(rank = row_number()) %>%
  filter(is.na(.draw))

y_rep %>%
  filter(.draw < 200) %>%
  bind_rows(y_true) %>%
  #filter(time_id <= 200) %>%
  group_by(.draw) %>%
  mutate(value = cumprod(1 + value/100)) %>%
  ggplot(aes(time_id, value, group = .draw, color = is.na(.draw))) +
  geom_line(alpha = 0.5)

y_rep %>%
  filter(.draw < 200) %>%
  bind_rows(y_true) %>%
  group_by(.draw) %>%
  mutate(value = cumprod(1 + value/100)) %>%
  filter(time_id == max(time_id)) %>%
  ggplot(aes(value, color = is.na(.draw))) +
  geom_histogram()
