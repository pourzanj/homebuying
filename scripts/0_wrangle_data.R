library(tidyverse)
library(lubridate)
library(tidyquant)
library(rstan)


# Functions

# Takes SPX df with date, pct_return and holdout_idx columns and fits up to
# desired holdout_idx
fit_ahead <- function(dat, max_holdout_idx = 0) {
  
  print(paste("Fitting up to holdout index", max_holdout_idx))
  
  # Filter out data past the holdout max
  holdout <-
    dat %>%
    filter(holdout_idx <= max_holdout_idx)
  
  # Fit
  fit <-
    stan("stan/stoch_vol.stan",
         data = list(T = nrow(holdout), y = holdout$pct_return),
         cores = 4, chains = 4, iter = 5e3,
         control = list(adapt_delta = 0.999))
  
  # Return
  fit
}

# Same as fit ahead but only gets certain parameters as tibble to save memory
fit_ahead_compressed <- function(dat, max_holdout_idx = 0) {
  
  print(paste("Fitting up to holdout index", max_holdout_idx))
  
  # Filter out data past the holdout max
  holdout <-
    dat %>%
    filter(holdout_idx <= max_holdout_idx)
  
  # Fit
  fit <-
    stan("stan/stoch_vol.stan",
         data = list(T = nrow(holdout), y = holdout$pct_return),
         cores = 4, chains = 4, iter = 5e3,
         control = list(adapt_delta = 0.999))
  
  # Return
  s <- extract(fit)
  tibble(max_holdout_idx,
         mu = s$mu, phi = s$phi, sigma = s$sigma,
         h_ahead = s$h_ahead, y_ahead = s$y_ahead) %>%
    mutate(sample_id = row_number())
}



# Filters the options quotes for a single day to include only those included in
# the VIX
vix_filter <- function(dat) {
  chain <-
    dat %>%
    mutate(price = (bid + ask) / 2) %>%
    select(strike, type, price) %>%
    pivot_wider(names_from = type, values_from = price) %>%
    mutate(diff = call - put)
  
  # Pull forward price F (don't account for interest since short term)
  F <-
    chain %>%
    filter(abs(diff) == min(abs(diff))) %>%
    pull(strike)
  
  K0 <- F
  
  # Get lowest included put strike
  lowest_inc_put <-
    dat %>%
    filter(type == "put") %>%
    arrange(type, strike) %>%
    select(strike, bid) %>%
    mutate(bid_zero = bid == 0) %>%
    mutate(lead_bid_zero = lead(bid_zero),
           lead2_bid_zero = lead(bid_zero, 2)) %>%
    mutate(poss_excl_point = bid_zero & lead_bid_zero & !lead2_bid_zero) %>%
    filter(poss_excl_point) %>%
    pull(strike) %>%
    max()
  
  # Highest included call
  highest_inc_call <-
    dat %>%
    filter(type == "call") %>%
    arrange(type, strike) %>%
    select(strike, bid) %>%
    mutate(bid_zero = bid == 0) %>%
    mutate(lag_bid_zero = lag(bid_zero),
           lag2_bid_zero = lag(bid_zero, 2)) %>%
    mutate(poss_excl_point = bid_zero & lag_bid_zero & !lag2_bid_zero) %>%
    filter(poss_excl_point) %>%
    pull(strike) %>%
    min()
  
  # Filter past lowest, highest, and anything in between with bid == 0
  # Finally filter out puts above K0 and calls below K0
  dat %>%
    filter(lowest_inc_put <= strike & strike <= highest_inc_call) %>%
    filter(bid != 0) %>%
    filter((type == "put" & strike <= K0) | (type == "call" & strike >= K0))
}

####################
# SPX Prices
####################
holdout_start_date <- ymd("2019-09-13")

spx_weekly <-
  tq_get("^GSPC",
         get = "stock.prices",
         from = "2010-09-07",
         to = "2020-08-28") %>%
  mutate(year = year(date), week = epiweek(date)) %>%
  group_by(week, year) %>%
  mutate(weekday = weekdays(date)) %>%
  filter(weekday == "Friday" | weekday == "Thursday") %>%
  filter(date == max(date)) %>%
  ungroup() %>%
  mutate(pct_return = 100*(close - lag(close)) / lag(close)) %>%
  select(date, weekday, close, pct_return) %>%
  mutate(holdout = date >= holdout_start_date) %>%
  group_by(holdout) %>%
  mutate(holdout_idx = row_number()) %>%
  ungroup() %>%
  mutate(holdout_idx = ifelse(holdout, holdout_idx, 0)) %>%
  na.omit()

write_rds(spx_weekly, "data/spx_weekly.rds")

####################
# Options
####################
spxw_opt <-
  read_csv("data/Batch__LX1d2sapX/SPXW_2019.csv") %>%
  bind_rows(read_csv("data/Batch__LX1d2sapX/SPXW_2020.csv")) %>%
  select(underlying, underlying_last,
         type, expiration, quotedate, strike,
         last, bid, ask, volume, openinterest) %>%
  mutate_at(c("expiration", "quotedate"), mdy) %>%
  mutate(expir_weekday = weekdays(expiration),
         quote_weekday = weekdays(quotedate)) %>%
  filter(expir_weekday == "Friday" & quote_weekday == "Friday") %>%
  mutate(dte = difftime(expiration, quotedate, units = "days")) %>%
  filter(dte == 7) %>%
  group_split(quotedate) %>%
  map_dfr(vix_filter) %>%
  mutate(mid_price = (bid + ask) / 2) %>%
  mutate(pct_strike = 100 * (strike - underlying_last) / underlying_last)

write_rds(spxw_opt, "data/spxw_opt.rds")

####################
# Fit Models
####################

# Fit just first week
fit <- fit_ahead(spx_weekly,  max_holdout_idx = 0)

s <- extract(fit)
y_ahead <- tibble(y_ahead = s$y_ahead)

write_rds(y_ahead, "data/y_ahead_09_06_19.rds")


# Fit all models
fits <- map_dfr(0:1, fit_ahead_compressed, dat = spx_weekly)
