library(tidyverse)
library(tidyquant)
library(rstan)
library(sgt)

sp500_daily <-
  tq_get("^GSPC", get = "stock.prices", from = "2010-01-01", to = "2020-10-07") %>%
  # mutate(year = year(date), week = week(date)) %>%
  # group_by(year, week) %>%
  # filter(date == min(date)) %>%
  # ungroup() %>%
  mutate(pct_return = (adjusted - lag(adjusted)) / lag(adjusted)) %>%
  mutate(pct_return = (open - lag(open)) / lag(open)) %>%
  na.omit

sp500_daily %>% ggplot(aes(date, pct_return)) + geom_line()
sp500_daily %>% ggplot(aes(pct_return)) + geom_histogram()

r <- sp500_daily$pct_return

fit <-
  stan(file = "sgt.stan",
       data = list(N = length(r), y = r,
                   m_p = 2, s_p = 100,
                   m_q = 2, s_q = 100),
       chains = 1, iter = 1e3, refresh = 1e2)

s <- extract(fit)

m <- stan_model("sgt.stan")
ml <- optimizing(m, data = list(N = length(r), y = r,
                                  m_p = 2, s_p = 100,
                                  m_q = 2, s_q = 100))

mlp <- which(s$lp__ == max(s$lp__))

sp500_daily %>%
  ggplot(aes(pct_return)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.001) +
  geom_function(fun = function(x) dsgt(x, s$mu[mlp], s$s[mlp], s$l[mlp], s$p[mlp], s$q[mlp], mean.cent = FALSE, var.adj = FALSE))

sgtized <-
  tibble(date = sp500_daily$date,
       r,
       p = psgt(r, s$mu[mlp], s$s[mlp], s$l[mlp], s$p[mlp], s$q[mlp],
                mean.cent = FALSE, var.adj = FALSE)) %>%
  mutate(q = qnorm(p))

sgtized %>%
  ggplot(aes(lag(p, 1), r)) +
  geom_point()  +
  geom_smooth()

sgtized %>%
  ggplot(aes(lag(r, 1), r^2)) +
  geom_point()  +
  geom_smooth()

sgtized %>%
  ggplot(aes(date, r)) +
  geom_line() 

sgtized %>% ggplot(aes(p)) + geom_histogram()

sgtized %>%
  ggplot(aes(lag(r), r)) +
  geom_point()
sgtized