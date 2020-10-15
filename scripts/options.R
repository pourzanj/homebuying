spx <-
  tq_get("^GSPC", get = "stock.prices", from = "1928-01-01", to = "2020-10-04")

spx_monthly <-
  spx %>%
  mutate(year = year(date), month = month(date)) %>%
  group_by(month, year) %>%
  filter(date == min(date)) %>%
  ungroup() %>%
  mutate(return = close / lag(close),
         pct_return = (close - lag(close)) / lag(close)) %>%
  mutate(log_return = log(close))

spx_monthly %>%
  select(date, close, pct_return) %>%
  pivot_longer(-date) %>%
  ggplot(aes(date, value)) +
  geom_line() +
  facet_grid(name ~ ., scales = "free")

spx_options <-
  read_csv("data/Sample_SPX_20151001_to_20151030.csv") %>%
  filter(underlying == "SPX") %>%
  mutate_at(c("expiration", "quotedate"), mdy) %>%
  mutate(avg_price = (ask + bid) / 2) %>%
  mutate(pct_strike = (strike - underlying_last) / underlying_last) %>%
  mutate(pct_avg_price = avg_price / underlying_last)

# Plot prices for one-month-out options on single day
spx %>%
  filter(date %in% ymd(c("2015-10-20", "2015-11-20")))

dat <-
  spx_options %>%
  filter(quotedate == "2015-10-20",
         expiration == "2015-11-20") %>%
  filter(bid > 0) %>%
  filter((type == "call" & pct_strike > 0) | (type == "put" & pct_strike < 0)) %>%
  select(type, pct_strike, pct_avg_price) %>%
  arrange(pct_strike) %>%
  mutate(width = pct_strike - lag(pct_strike),
         height = (pct_avg_price + lag(pct_avg_price)) / 2) %>%
  mutate(area = width * height) %>%
  filter(!is.na(area)) %>%
  mutate(dens = pct_avg_price / sum(area)) %>%
  mutate(rr = ((2089-2031) / 2031)) %>%
  mutate(itm = ifelse(type == "call",
                      pct_strike < rr,
                      pct_strike > rr)) %>%
  mutate(rp = ifelse(itm, abs(rr - pct_strike) - pct_avg_price, -pct_avg_price))

# Price vs strike
dat %>%
  ggplot(aes(pct_strike, pct_avg_price)) +
  geom_line() +
  geom_point(aes(color = type))

# Realized payout vs strike
dat %>%
  ggplot(aes(pct_strike, rp)) +
  geom_line() +
  geom_point(aes(color = type))

# Density vs historical 1 month density when preivous month was close to zero return
dat %>%
  ggplot() +
  geom_line(aes(pct_strike, dens)) +
  geom_point(aes(pct_strike, dens, color = type)) +
  geom_density(aes(pct_return), data = filter(spx_monthly, -0.01 < lag(pct_return) & lag(pct_return) < 0.01))

# Compare implied vs historical CDF
p_h <-
  filter(spx_monthly, -0.01 < lag(pct_return) & lag(pct_return) < 0.01) %>%
  arrange(pct_return) %>%
  mutate(dummy = 1) %>%
  mutate(p_h = cumsum(dummy) / sum(dummy))

h <- p_h$pct_return

dat %>%
  mutate(p_i = cumsum(pct_avg_price) / sum(pct_avg_price),
         p_h = map_dbl(pct_strike, function(p) mean(h < p))) %>%
  select(pct_strike, p_i, p_h) %>%
  pivot_longer(-pct_strike) %>%
  ggplot(aes(pct_strike, value, color = name)) +
  geom_line()

dat %>%
  mutate(p_i = cumsum(pct_avg_price) / sum(pct_avg_price),
         p_h = map_dbl(pct_strike, function(p) mean(h < p))) %>%
  mutate(p_d = p_i - p_h) %>%
  ggplot(aes(pct_strike, p_d)) +
  geom_line()

# Simulate process
f <- function(i) sum(t(rmultinom(12, 1, c(0.62, 0.32, 0.06))) %*% c(0.02224551, -0.02865269, -0.02865269))
1:1e3 %>% map_dbl(f) %>% qplot()

# Get todays
spy_last <- 333.87

spy_opt <-
  getOptionChain("SPY", Exp = "2020-11-06") %>%
  map2(names(.), function(df, name) mutate(df, type = name)) %>%
  bind_rows() %>%
  mutate(underlying_last = spy_last) %>%
  mutate(pct_strike = (Strike - spy_last) / spy_last) %>%
  filter((type == "calls" & pct_strike > 0) | (type == "puts" & pct_strike < 0)) %>%
  mutate(avg_price = (Bid + Ask) / 2) %>%
  mutate(pct_avg_price = avg_price / underlying_last) %>%
  arrange(pct_strike)


spy_opt %>%
  ggplot(aes(pct_strike, pct_avg_price)) +
  geom_line() +
  geom_point(aes(color = type))
