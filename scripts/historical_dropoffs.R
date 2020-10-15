library(tidyverse)
library(tidyquant)
library(lubridate)
library(tibbletime)

sp500_daily <- tq_get("^GSPC", get = "stock.prices", from = "2010-01-01", to = "2020-09-26")

sp500 <-
  tq_get("^GSPC", get = "stock.prices", from = "1928-01-01", to = "2020-09-28") %>%
  mutate(weekday = weekdays(date)) %>%
  mutate(last_weekday = lag(weekday))

sp500_weekly <-
  sp500 %>%
  # filter(weekday == "Monday" |
  #          (weekday == "Tuesday" & last_weekday == "Friday") |
  #          (weekday == "Tuesday" & last_weekday == "Thursday") |
  #          (weekday == "Wednesday" & last_weekday == "Friday")) %>%
  mutate(year = year(date), week = week(date)) %>%
  group_by(week, year) %>%
  filter(date == min(date)) %>%
  ungroup() %>%
  mutate(days_since_last = date - lag(date)) %>%
  #select(date, adjusted) %>%
  mutate(ath = cummax(adjusted)) %>%
  mutate(at_ath = adjusted == ath) %>%
  mutate(return_from_ath = adjusted / ath) %>%
  mutate(is_low = return_from_ath < 0.90) %>%
  mutate(change = is_low - lag(is_low)) %>% 
  filter(date != "2001-01-08")

# Plot time series
sp500_weekly  %>%
  pivot_longer(-date) %>%
  ggplot(aes(date, value)) +
  geom_line() +
  facet_grid(name ~ ., scales = "free") +
  geom_hline(yintercept = 0.90, color = "red") +
  geom_vline(xintercept = as.Date("2001-03-19"), color = "red") +
  geom_vline(xintercept = as.Date("2006-10-16"), color = "red") +
  geom_vline(xintercept = as.Date("2008-01-22"), color = "red") + 
  geom_vline(xintercept = as.Date("2013-04-01"), color = "red")
  
# Get down periods how long they last and how low they go
down_periods <-
  sp500_weekly %>%
  filter(change == 1 | at_ath) %>%
  filter((change == 1 & lag(at_ath)) | (at_ath & lag(change == 1))) %>%
  mutate(start_date = date, end_date = lead(date)) %>%
  select(start_date, end_date, open, ath, return_from_ath) %>%
  filter(return_from_ath != 1) %>%
  mutate(crash_id = row_number()) %>%
  mutate(length = end_date - start_date)

# Use a complex join to apend on the lowest low during the crash
down_periods <-
  sp500_weekly %>%
  mutate(dummy = TRUE) %>%
  full_join(select(down_periods, crash_id, start_date, end_date) %>% mutate(dummy = TRUE)) %>%
  filter(start_date <= date & date <= end_date) %>%
  group_by(crash_id) %>%
  mutate(low = min(return_from_ath)) %>%
  filter(return_from_ath == low) %>%
  ungroup() %>%
  mutate(time_to_bottom = date - start_date,
         time_bottom_to_end = end_date - date) %>%
  select(crash_id, low, time_to_bottom, time_bottom_to_end) %>%
  inner_join(down_periods) %>%
  select(start_date, end_date, length, open, ath, return_from_ath, low, time_to_bottom, time_bottom_to_end)

sp500_weekly %>%
  select(date, open, ath, return_from_ath) %>%
  pivot_longer(-date) %>%
  ggplot() +
  geom_line(aes(date, value)) +
  facet_grid(name ~ ., scales = "free") +
  geom_vline(aes(xintercept = start_date), color = "red", data = down_periods) +
  geom_vline(aes(xintercept = end_date), color = "green", data = down_periods) +
  geom_rect(aes(xmin = start_date, xmax = end_date,
                ymin = 0, ymax = 1),
            alpha = 0.3, data = down_periods) +
  geom_hline(aes(yintercept = 0.95), color = "blue")

down_periods %>%
  select(length, return_from_ath, low, time_to_bottom, time_bottom_to_end) %>%
  mutate_all(as.numeric) %>%
  mutate(crash_id = row_number()) %>%
  pivot_longer(-crash_id) %>%
  ggplot(aes(crash_id, value)) +
  geom_line() +
  geom_point() +
  facet_grid(name ~ ., scales = "free")


down_periods %>% ggplot(aes(total_periods)) + geom_histogram(binwidth = 10)
down_periods$total_periods %>% quantile(probs = seq(0.1, 0.9, 0.1))

down_periods %>% ggplot(aes(low_rfath)) + geom_histogram()
down_periods$low_rfath %>% quantile(probs = seq(0.1, 0.9, 0.1))

down_periods %>%
  ggplot(aes(total_periods, low_rfath)) +
  geom_point() +
  scale_x_log10()

resevoir <- function(rfath) {
  N <- length(return_from_ath)
  debt <- rep(0.0, N)
  leftover <- rep(0.0, N)
  purchase <- rep(0.0, N)
  
  for(i in 2:N) {
    purchase[i] <- 100 * (1 - rfath)
    leftover[i] <- 100 - purchase[i]
    
    if(leftover[i] > 0.0) {
      # If there's leftovers pay off any outstanding debt then buy more
      if(res[i-1] < 0.0) {
        
      }
    } else {
      # If there's no leftovers we finance the purchase with debt
      
    }

  }
}

# Getting rolling returns
sp500_monthly <-
  sp500_weekly %>%
  mutate(year = year(date), month = month(date)) %>%
  group_by(year, month) %>%
  filter(date == min(date)) %>%
  ungroup()
  
rolling_prod <- rollify(prod, window = 13)
rolling_sum <- rollify(sum, window = 12)

ret <- function(open) tail(open, 1) / head(open, 1)
rolling_ret <- rollify(ret, window = 13)

inv_total <- function(monthly_return) {
  N <- length(monthly_return)
  r <- 100 * cumprod(head(rev(monthly_return), N-1))
  sum(r)
}
inv_total <- function(open) {
  N <- length(open)
  r <- 100 * open[N] / head(open, N-1)
  sum(r)
}
rolling_inv_total <- rollify(inv_total, window = 13)

rolling_returns <-
  sp500_monthly %>%
  # mutate(monthly_return = open / lag(open)) %>%
  # filter(!is.na(monthly_return)) %>%
  #mutate(lead_monthly_return = lead(monthly_return)) %>%
  mutate(past_year_return = rolling_ret(open)) %>%
  filter(!is.na(past_year_return)) %>%
  mutate(pv_1200_inv_last_year = past_year_return * 1200) %>%
  mutate(pv_100_inv_monthly_year = rolling_inv_total(open)) %>%
  filter(!is.na(pv_100_inv_monthly_year))

rolling_returns %>%
  select(date, pv_100_inv_monthly_year, pv_1000_inv_last_year) %>%
  pivot_longer(-date) %>%
  ggplot(aes(date, value, color = name)) +
  geom_line()
  

qplot(rolling_returns$pv_100_inv_monthly_year)
mean(rolling_returns$pv_100_inv_monthly_year < 1000)
median(rolling_returns$pv_100_inv_monthly_year)
qplot(rolling_returns$pv_100_inv_last_year * 10)
mean((rolling_returns$pv_100_inv_last_year * 10) < 1000)
median(rolling_returns$pv_100_inv_last_year * 10)

rolling_returns %>%
  ggplot(aes(date, pv_100_inv_monthly_year)) +
  geom_line() +
  geom_hline(yintercept = 1000, color = "red")

# conditional
weekly_cond <-
  sp500_weekly %>%
  mutate(returns = open / lag(open)) %>%
  mutate(log_returns = log(returns)) %>%
  # mutate(previous_bin = lag(cut(log_returns,
  #                               c(-Inf, quantile(log_returns, c(0.02, 0.1, 0.25, 0.5, 0.75, 0.9, 0.98),
  #                                                na.rm = TRUE), Inf)))) %>%
  mutate(previous_bin = lag(cut(log_returns,
                                c(-Inf, seq(-0.05, 0.05, by = 0.01), Inf)))) %>%
  filter(!is.na(previous_bin))
  
weekly_cond %>%
  ggplot(aes(log_returns)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.01) +
  facet_grid(previous_bin ~ .)

# Fit conditional in Stan
fit_sgt <- function(y) {
  fit <-
    stan(file = "sgt.stan",
         data = list(N = length(y), y = y,
                     m_p = 2, s_p = 100,
                     m_q = 2, s_q = 100),
         chains = 1, iter = 1e3, refresh = 1e2)
}

extract_samples <- function(fit) {
  s <- extract(fit)
  tibble(mu = s$mu, sigma = s$s, lambda = s$l, p = s$p, q = s$q) %>%
    mutate(sample_id = row_number())
}

fits <-
  weekly_cond %>%
  group_split(previous_bin) %>%
  map(function(df) df$log_returns) %>%
  map(fit_sgt)

groups <- 
  weekly_cond %>%
  group_split(previous_bin) %>%
  map(function(df) df$previous_bin[1])

# Compare SGT parameters for different value of lag of return
fits %>%
  map(extract_samples) %>%
  map2(groups, function(df, group) mutate(df, previous_bin = group)) %>%
  bind_rows() %>%
  pivot_longer(c(-sample_id, -previous_bin)) %>%
  group_by(previous_bin, name) %>%
  summarize(q10 = quantile(value, 0.1),
            q50 = quantile(value, 0.5),
            q90 = quantile(value, 0.9)) %>%
  ungroup() %>%
  ggplot() +
  geom_pointrange(aes(previous_bin, q50, ymin = q10, ymax = q90)) +
  facet_grid(name ~ ., scales = "free")

# Plot density posterior over samples
grid <-
  tibble(mu = s$mu, sigma = s$s, lambda = s$l, p = s$p, q = s$q) %>%
  sample_n(20) %>%
  crossing(x = seq(-0.2, 0.2, by = 0.001)) %>%
  mutate(density = pmap_dbl(., dsgt, mean.cent = FALSE, var.adj = FALSE)) %>%
  inner_join(mutate(distinct(., mu, sigma, lambda, p, q), sample_id = row_number()))

weekly_cond %>%
  filter(previous_bin == "(-Inf,-0.05]") %>%
  ggplot(aes(log_returns)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.003) +
  geom_line(aes(x, density, group = sample_id), data = grid)
  

# Testing sgt_lpdf
# expose_stan_functions(fit)
# sgt_lpdf(array(0.1), 0, 0.02, 0, 2.0, 0.5)
# dsgt(0.1, mu = 0, sigma = 0.02, lambda = 0, p = 2.0, q = 0.5,
#      mean.cent = FALSE, var.adj = FALSE, log = TRUE)
