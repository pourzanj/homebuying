us_vs_inter <-
  tq_get(c("VTMGX", "VTSAX", "MCHI"), get = "stock.prices", from = "1928-01-01", to = "2020-10-06")

us_vs_inter %>%
  mutate(year = year(date)) %>%
  filter(year > 2011) %>%
  group_by(symbol, year) %>%
  filter(date == min(date)) %>%
  ungroup() %>%
  group_by(symbol) %>%
  mutate(start_adjusted = adjusted * ifelse(date == min(date), 1, 0)) %>%
  mutate(adjusted = adjusted / max(start_adjusted)) %>%
  ungroup() %>%
  select(-start_adjusted) %>%
  group_by(symbol) %>%
  arrange(symbol, date) %>%
  mutate(return = (adjusted - lag(adjusted)) / lag(adjusted)) %>%
  ungroup() %>%
  select(symbol, date, adjusted, return) %>%
  pivot_longer(c(-date, -symbol)) %>%
  ggplot(aes(date, value, color = symbol)) +
  geom_line() +
  geom_point() +
  facet_grid(name ~ ., scales = "free")
  
mchi <-
  tq_get(c("MCHI"), get = "stock.prices", from = "1928-01-01", to = "2020-10-06")

mchi %>%
  ggplot(aes(date, adjusted)) +
  geom_line()
