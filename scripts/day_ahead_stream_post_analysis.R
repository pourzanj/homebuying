kOutdir <- "data/backtest_sv13_1_1992_2022/"
y_oos <- readRDS(paste0(kOutdir, "y_oos.RDS"))
daily_summary <- readRDS(paste0(kOutdir, "daily_summary.RDS"))
par_summary <- readRDS(paste0(kOutdir, "par_summary.RDS"))

key <-
  daily_summary %>%
  transmute(idx = row_number(),
            date,
            y_true = y,
            r_true = exp(y_true/100),
            mean) %>%
  filter(!is.na(mean)) %>%
  select(-mean) %>%
  mutate(epoch_idx = ((row_number()-1) %/% 1) + 1)

y_oos_df <-
  y_oos %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(.draw = row_number()) %>%
  pivot_longer(-.draw) %>%
  mutate(idx = as.integer(parse_number(name))) %>%
  inner_join(key) %>%
  select(-name, -idx)

epoch_returns <-
  y_oos_df %>%
  group_by(epoch_idx, .draw) %>%
  summarize(date = min(date),
            R = prod(exp(value/100)),
            R_true = prod(exp(y_true/100))) %>%
  ungroup()

epoch_pvals <-
  epoch_returns %>%
  group_by(date) %>%
  summarize(p = mean(R_true > R)) %>%
  mutate(z = qnorm(p))

epoch_pvals %>%
  ggplot(aes(p)) +
  geom_histogram(breaks = seq(0, 1, by = 0.1))

epoch_pvals %>%
  ggplot(aes(date, p)) +
  geom_point()

epoch_pvals %>%
  mutate(year = year(date)) %>%
  # filter(year(date) >= 2021) %>%
  ggplot(aes(p)) +
  geom_histogram(breaks = seq(0, 1, by = 0.2)) +
  facet_wrap(year ~ .)

epoch_pvals %>%
  filter(year(date) >= 2021) %>%
  ggplot(aes(date, p)) +
  geom_point()

##################
N <- 252
sim_ar1 <- function(.draw, phi, sigma, nu, h0, N) {
  h <- rep(0, N)
  h[1] <- h0
  for(n in 2:N) {
    # h[n] <- phi*h[n-1] + sigma*rnorm(1)
    # h[n] <- phi*h[n-1] + sigma*sqrt((nu-2)/nu)*rt(1, nu)
    h[n] <- phi*h[n-1] + sigma*exp(nu*h[n-1])*rnorm(1)
  }
  
  tibble(.draw = .draw, i = 1:N, h)
}

phi <- 0.988
c <- 0.015
nu <- 1.539

N <- 2520
sim_ar1 <- function(h0) {
  z <- rep(NA_real_, N)
  h <- rep(0, N)
  h[1] <- h0
  for(n in 2:N) {
    z[n] <- rpois(1, phi*h[n-1]/c)
    h[n] <- rgamma(1, nu + z[n], scale = c)
  }
  
  tibble(i = 1:N, z, h, s = sqrt(h))
}

sim_ar1(0.36) %>% ggplot(aes(i, s)) + geom_line()
