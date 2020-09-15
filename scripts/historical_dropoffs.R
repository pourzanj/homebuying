library(tidyverse)
library(tidyquant)
library(lubridate)

sp500 <-
  tq_get("^GSPC", get = "stock.prices", from = "1928-01-01", to = "2020-09-14") %>%
  mutate(weekday = weekdays(date)) %>%
  mutate(last_weekday = lag(weekday))

sp500_weekly <-
  sp500 %>%
  filter(weekday == "Monday" |
           (weekday == "Tuesday" & last_weekday == "Friday") |
           (weekday == "Tuesday" & last_weekday == "Thursday") |
           (weekday == "Wednesday" & last_weekday == "Friday")) %>%
  mutate(days_since_last = date - lag(date)) %>%
  select(date, open) %>%
  mutate(ath = cummax(open)) %>%
  mutate(at_ath = open == ath) %>%
  mutate(return_from_ath = open / ath) %>%
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
  geom_hline(aes(yintercept = 0.9), color = "blue")

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