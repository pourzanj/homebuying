# compare simulated vs actual returns on a monthly basis
months <-
  vol %>%
  inner_join(select(test, date, ret_spy)) %>%
  mutate(year = year(date), month = month(date))

month_summary <-
  months %>%
  group_by(month, year) %>%
  summarize(start = head(date, 1), end = tail(date, 1),
            n = n(),
            r = 100*exp(mean(log(1+ret_spy/100)))-100) %>%
  ungroup() %>%
  arrange(year, month)

sim_month_draws <-
  y_rep_holdout %>%
  t() %>%
  as_tibble() %>%
  bind_cols(select(months, year, month)) %>%
  pivot_longer(c(-year, -month)) %>%
  group_by(year, month, name) %>%
  summarize(r = 100*exp(mean(log(1+value/100)))-100) %>%
  ungroup()

sim_month_summary <-
  sim_month_draws %>%
  inner_join(select(month_summary, month, year, r_obs = r)) %>%
  group_by(year, month) %>%
  summarize(p = mean(r <= r_obs),
            q05 = quantile(r, 0.05),
            q50 = quantile(r, 0.5),
            q95 = quantile(r, 0.95)) %>%
  mutate(z = qnorm(p)) %>%
  ungroup()

sim_month_summary %>%
  ggplot(aes(month, q50)) +
  geom_point() +
  geom_errorbar(aes(ymin = q05, ymax = q95)) +
  geom_point(aes(month, r), color = "red", data = month_summary) +
  facet_wrap(year ~ .)

# Break up time into when the algo switches leverage regimes and compare
# simulated versus actual returns
y_rep_holdout <- readRDS("data/sv_with_leverage_01_13_20/y_rep_holdout.RDS")
vol <- readRDS("data/sv_with_leverage_01_13_20/vol.RDS")
test <- readRDS("data/sv_with_leverage_01_13_20/test.RDS")

y_rep_holdout <- readRDS("data/orig2021_10_19_21/y_rep_holdout.RDS")
vol <- readRDS("data/orig2021_10_19_21/vol.RDS")
test <- readRDS("data/orig2021_10_19_21/test.RDS")

switches <-
  vol %>%
  mutate(lev_switch = k != lag(k)) %>%
  na.omit() %>%
  mutate(switch_id = cumsum(lev_switch))

switch_summary <-
  switches %>%
  inner_join(select(test, date, ret_spy)) %>%
  group_by(switch_id) %>%
  summarize(start = head(date, 1), end = tail(date, 1),
            k = head(k, 1), n = n(),
            r = 100*exp(mean(log(1+ret_spy/100)))-100)

sim_switch_draws <-
  y_rep_holdout[,-1] %>%
  t() %>%
  as_tibble() %>%
  mutate(date = switches$date,
         switch_id = switches$switch_id) %>%
  pivot_longer(c(-date, -switch_id)) %>%
  group_by(switch_id, name) %>%
  summarize(r = 100*exp(mean(log(1+value/100)))-100) %>%
  ungroup()
  
sim_switch_summary <-
  sim_switch_draws %>%
  inner_join(select(switch_summary, switch_id, r_obs = r)) %>%
  group_by(switch_id) %>%
  summarize(p = mean(r <= r_obs),
            q05 = quantile(r, 0.05),
            q50 = quantile(r, 0.5),
            q95 = quantile(r, 0.95)) %>%
  mutate(z = qnorm(p))

switch_summary %>%
  inner_join(sim_switch_summary) %>%
  ggplot(aes(start, z, color = factor(k))) +
  geom_point()

switch_summary %>%
  inner_join(sim_switch_summary) %>%
  ggplot(aes(z)) +
  geom_histogram() +
  facet_grid(k ~ .)

switch_summary %>%
  inner_join(sim_switch_summary) %>%
  ggplot(aes(start, q50, color = factor(k))) +
  geom_point() +
  geom_errorbar(aes(ymin = q05, ymax = q95)) +
  geom_point(aes(start, r, color = NULL))

# Compare returns versus distribution of predicted returns at differing predicted
# values of expected returns
y_rep_holdout <- readRDS("data/vol_prem_10_16_21/y_rep_holdout.RDS")
vol <- readRDS("data/vol_prem_10_16_21/vol.RDS")
test <- readRDS("data/vol_prem_10_16_21/test.RDS")

y_rep_holdout <- readRDS("data/vol_prem_tighter_prior_laplace_10_17_21/y_rep_holdout.RDS")
test <- readRDS("data/vol_prem_tighter_prior_laplace_10_17_21/test.RDS")

expect_daily_returns <- apply(y_rep_holdout, 2, function(r) 100*exp(mean(log(1+r/100)))-100)

test %>%
  select(date, ret_spy) %>%
  mutate(expec = expect_daily_returns) %>%
  # mutate(expec_bin = cut(expec, quantile(expec, seq(0,1,by=0.1)))) %>%
  mutate(expec_bin = cut(expec, quantile(expec, c(0, 0.05, 0.1, 0.2, 0.5, 0.75, 1.0)))) %>%
  pull(expec_bin) -> expec_bin

t(y_rep_holdout) %>%
  as_tibble() %>%
  mutate(date = test$date, expec_bin = expec_bin) %>%
  pivot_longer(c(-date, -expec_bin)) %>%
  group_by(expec_bin, name) %>%
  summarize(avg_daily_ret = 100*exp(mean(log(1+value/100)))-100) %>%
  group_by(expec_bin) %>%
  summarize(q05 = quantile(avg_daily_ret, 0.05),
            q50 = quantile(avg_daily_ret, 0.5),
            q95 = quantile(avg_daily_ret, 0.95)) -> avg_daily_ret_quant

test %>%
  select(date, ret_spy) %>%
  mutate(expec = expect_daily_returns) %>%
  # mutate(expec_bin = cut(expec, quantile(expec, seq(0,1,by=0.1)))) %>%
  mutate(expec_bin = expec_bin) %>%
  group_by(expec_bin) %>%
  summarize(avg_daily_ret = 100*exp(mean(log(1+ret_spy/100)))-100) %>%
  na.omit() -> avg_daily_ret_obs

avg_daily_ret_quant %>%
  na.omit() %>%
  ggplot(aes(expec_bin, q50)) +
  geom_point() +
  geom_errorbar(aes(ymin = q05, ymax = q95)) +
  geom_point(aes(expec_bin, avg_daily_ret), color = "red", data = avg_daily_ret_obs)


# How parameters change over time on independent training sets
parameters %>% mutate(month = month(date)) %>% filter(month %in% c(1, 7)) %>% group_by(variable, month, year(date)) %>% filter(date == min(date)) %>% ungroup() %>% ggplot(aes(date, median)) + geom_point() + geom_errorbar(aes(ymin = q5, ymax = q95)) + facet_grid(variable ~ ., scales = "free")