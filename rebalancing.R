data <-
  tq_get(c("SPY", "DLR", "EQIX", "BLDR"),
         get = "stock.prices",
         from = "2010-01-01", to = "2020-09-30")

data %>%
  filter(symbol == "SPY") %>%
  select(date, close, adjusted) %>%
  mutate(adjusted = adjusted / adjusted[1]) %>%
  pivot_longer(-date) %>%
  ggplot(aes(date, value, color = name)) +
  geom_line()


data %>%
  group_by(symbol) %>%
  mutate(close = close / close[1],
         adjusted = adjusted / adjusted[1]) %>%
  ungroup() %>%
  select(date, symbol, close, adjusted) %>%
  pivot_longer(c(-date, -symbol)) %>%
  ggplot(aes(date, value, color = symbol)) +
  geom_line() +
  facet_grid(name ~ .) + scale_y_log10()

# Tail and SPY
data <-
  tq_get(c("SPY", "TAIL"),
         get = "stock.prices",
         from = "2010-01-01", to = "2020-09-30")

# Plot overall
data %>%
  select(symbol, date, adjusted) %>%
  pivot_wider(names_from = "symbol", values_from = "adjusted") %>%
  filter(!is.na(TAIL)) %>%
  pivot_longer(-date) %>%
  group_by(name) %>%
  mutate(value = value / head(value, 1)) %>%
  ungroup() %>%
  ggplot(aes(date, value, color = name)) +
  geom_line() 

# Portfolio 1: Buy and hold 50%
data %>%
  select(symbol, date, adjusted) %>%
  pivot_wider(names_from = "symbol", values_from = "adjusted") %>%
  filter(!is.na(TAIL)) %>%
  mutate(portfolio = 0.5 * (SPY / head(SPY, 1)) + 0.5 * (TAIL / head(TAIL, 1))) %>%
  mutate(portfolio = portfolio / head(portfolio, 1)) %>%
  ggplot(aes(date, portfolio)) +
  geom_line()

# Portfolio 1: Buy and hold 50%
df <-
  data %>%
  # mutate(month = month(date), year = year(date)) %>%
  # group_by(month, year) %>%
  mutate(week = week(date), year = year(date)) %>%
  group_by(week, year) %>%
  filter(date == min(date)) %>%
  ungroup() %>%
  select(symbol, date, adjusted) %>%
  pivot_wider(names_from = "symbol", values_from = "adjusted") %>%
  filter(!is.na(TAIL)) %>%
  mutate_at(c("SPY", "TAIL"), function(x) x / head(x, 1))

rebalancing_portfolio <- function(date, s1, s2) {
  N <- length(s1)
  
  Pi <- rep(2.0, N)
  shares1 <- rep(1.0, N)
  shares2 <- rep(1.0, N)
  
  w1 <- rep(0.5, N)
  w2 <- rep(0.5, N)
  
  for(t in 2:N) {
    Pi[t] <- shares1[t-1] * s1[t] + shares2[t-1] * s2[t]
    w1[t] <- (shares1[t-1] * s1[t]) / Pi[t]
    w2[t] <- (shares2[t-1] * s2[t]) / Pi[t]
    
    if(w1[t] > 0.5) {
      excess_weight <- w1[t] - 0.5
      diff_dollars <- excess_weight * Pi[t]
      shares1[t] <- shares1[t-1] - diff_dollars / s1[t]
      shares2[t] <- shares2[t-1] + diff_dollars / s2[t]
    } else {
      excess_weight <- w2[t] - 0.5
      diff_dollars <- excess_weight * Pi[t]
      shares1[t] <- shares1[t-1] + diff_dollars / s1[t]
      shares2[t] <- shares2[t-1] - diff_dollars / s2[t]
    }
  }
  tibble(date, s1, s2, Pi, w1, w2, shares1, shares2)
}

rebalancing_portfolio(df$date, df$SPY, df$TAIL) %>%
  ggplot(aes(date, Pi)) +
  geom_line()

rebalancing_portfolio(df$date, df$SPY, df$TAIL) %>%
  mutate(returns = log(Pi / lag(Pi)), PiNorm = Pi / 2) %>%
  pivot_longer(-date) %>%
  ggplot(aes(date, value)) +
  geom_line() +
  geom_point() +
  facet_grid(name ~ ., scales = "free")
