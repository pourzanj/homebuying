library(tidyverse)
library(tidyquant)
library(lubridate)
library(tibbletime)

sp500 <-
  tq_get("^GSPC", get = "stock.prices", from = "1928-01-01", to = "2020-09-28") %>%
  mutate(weekday = weekdays(date)) %>%
  mutate(last_weekday = lag(weekday))

sp500_weekly <-
  sp500 %>%
  mutate(year = year(date), week = week(date)) %>%
  group_by(week, year) %>%
  filter(date == min(date)) %>%
  ungroup() %>%
  mutate(return = adjusted / lag(adjusted),
         pct_return = (adjusted - lag(adjusted)) / lag(adjusted)) %>%
  mutate(log_return = log(return)) %>%
  filter(date > "2010-01-08")

vix_weekly <-
  tq_get("VIX", get = "stock.prices", from = "1928-01-01", to = "2020-09-28") %>%
  mutate(year = year(date), week = week(date)) %>%
  group_by(week, year) %>%
  filter(date == min(date)) %>%
  ungroup() %>%
  mutate(return = adjusted / lag(adjusted),
         pct_return = (adjusted - lag(adjusted)) / lag(adjusted)) %>%
  mutate(log_return = log(return)) %>%
  filter(date > "2010-01-08")


# Fit Stan model
r <- sp500_weekly$pct_return

fit <-
  stan("garch_1_1.stan", data = list(T = length(r), r = r),
       chains = 1, iter = 1e3)

print(fit, pars = c("mu", "alpha0", "alpha1", "beta1", "sigma1"))

# Extract and plot simulated data
s <- extract(fit)

r_rep <-
  s$r_rep %>%
  rbind(r) %>%
  as_tibble() %>%
  mutate(source = c(rep("sim", 500), "true"),
         sample_id = c(1:500, NA)) %>%
  pivot_longer(c(-sample_id, -source)) %>%
  mutate(t = as.integer(str_remove(name, "V"))) %>%
  select(-name)

rand_idx <- sample(1:500, 4)
r_rep %>%
  filter(sample_id %in% rand_idx | source == "true") %>%
  ggplot(aes(t, value, color = source)) +
  geom_line() +
  facet_grid(sample_id ~ .)

r_rep %>%
  filter(sample_id %in% rand_idx | source == "true") %>%
  ggplot(aes(t, value^2, color = source)) +
  geom_line() +
  facet_grid(sample_id ~ .)

r_rep %>%
  filter(sample_id %in% rand_idx | source == "true") %>%
  ggplot(aes(t, value^3, color = source)) +
  geom_line() +
  facet_grid(sample_id ~ .)

r_rep %>%
  filter(sample_id %in% rand_idx | source == "true") %>%
  ggplot(aes(t, value^4, color = source)) +
  geom_line() +
  facet_grid(sample_id ~ .)

# Histogram
rand_idx <- sample(1:500, 15)
r_rep %>%
  filter(sample_id %in% rand_idx | source == "true") %>%
  ggplot(aes(value, fill = source)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.005, center = 0.0025, closed = "right") +
  facet_wrap(sample_id ~ .)

idx <- 538
r_rep %>%
  filter(t == idx) %>%
  filter(source == "sim") %>%
  ggplot(aes(value)) +
  geom_histogram() +
  geom_vline(xintercept = r[idx], color = "red")

mode_count <- function(x) {
  bins <- cut(x, seq(-0.5, 0.5, by = 0.005))
  max_count <-
    tibble(x, bin = bins) %>%
    group_by(bin) %>%
    tally() %>%
    filter(n == max(n)) %>%
    head(1) %>%
    pull(n)
  
  # normalize by total length to get density
  max_count / length(x)
}

quantile_dist <- function(x, start, end) {
  q <- quantile(x, probs = c(start, end))
  q[2] - q[1]
}

test_stats <-
  r_rep %>%
  group_by(source, sample_id) %>%
  summarize(mean = mean(value),
            sd = sd(value),
            expon2 = mean(value^2),
            expon3 = mean(value^3),
            expon4 = mean(value^4),
            sd_p = sd(pmax(value, 0)),
            sd_m = sd(pmin(value, 0)),
            qdist_25_75 = quantile_dist(value, 0.25, 0.75),
            qdist_1_9 = quantile_dist(value, 0.1, 0.9),
            qdist_1_5 = quantile_dist(value, 0.1, 0.5),
            qdist_5_9 = quantile_dist(value, 0.5, 0.9),
            max = max(value),
            min = min(value),
            num_gt_5 = sum(value > 0.05),
            num_lt_5 = sum(value < -0.05),
            num_gt_10 = sum(value > 0.10),
            num_lt_10 = sum(value < -0.10),
            mode_count = mode_count(value)) %>%
  ungroup() %>%
  pivot_longer(c(-source, -sample_id))

test_stats %>%
  filter(source == "sim") %>%
  ggplot(aes(value)) +
  geom_histogram() +
  geom_vline(aes(xintercept = value), color = "red", data = filter(test_stats, source == "true")) +
  facet_wrap(. ~ name, scales = "free")

#################
## Try SGT marginally on all returns
#################
fit <-
  stan("sgt.stan",
       data = list(N = length(r), y = r,
                   m_p = 2, s_p = 1e2,
                   m_q = 2, s_q = 1e2),
       chains = 1, iter = 1e3, refresh = 25)

fit

s <- extract(fit)

# See if we captured tails and pointiness
post_sgt_par_draws <-
  tibble(mu = s$mu, sigma = s$s, lambda = s$l, p = s$p, q = s$q) %>%
  mutate(sample_id = row_number()) %>%
  sample_n(15)

sgt_draws <-
  post_sgt_par_draws %>%
  select(-sample_id) %>%
  pmap(rsgt, n = 365, mean.cent = FALSE, var.adj = FALSE) %>%
  map2(post_sgt_par_draws$sample_id,
       function(x, id) tibble(sample_id = id, source = "sim", pct_return = x)) %>%
  bind_rows()

grid <-
  post_sgt_par_draws %>%
  select(-sample_id) %>%
  crossing(x = seq(-0.2, 0.2, by = 0.001)) %>%
  mutate(density = pmap_dbl(., dsgt, mean.cent = FALSE, var.adj = FALSE)) %>%
  inner_join(post_sgt_par_draws)

true_sim_combined <-
  sp500_weekly %>%
  mutate(sample_id = as.integer(NA), source = "true") %>%
  mutate(pct_return = (open - lag(open)) / lag(open)) %>%
  select(sample_id, source, pct_return) %>%
  bind_rows(sgt_draws)
  
true_sim_combined %>%
  ggplot(aes(pct_return)) +
  geom_histogram(aes(y = ..density.., fill = source), binwidth = 0.005, center = 0.0025, closed = "right") +
  geom_line(aes(x, density, group = sample_id), data = grid) +
  facet_wrap(sample_id ~ .)

# PPCs
post_sgt_par_draws <-
  tibble(mu = s$mu, sigma = s$s, lambda = s$l, p = s$p, q = s$q) %>%
  mutate(sample_id = row_number())

sgt_draws <-
  post_sgt_par_draws %>%
  select(-sample_id) %>%
  pmap(rsgt, n = 365, mean.cent = FALSE, var.adj = FALSE) %>%
  map2(post_sgt_par_draws$sample_id,
       function(x, id) tibble(sample_id = id, source = "sim", pct_return = x)) %>%
  bind_rows()

true_sim_combined <-
  sp500_weekly %>%
  mutate(sample_id = as.integer(NA), source = "true") %>%
  select(sample_id, source, pct_return) %>%
  bind_rows(sgt_draws)

test_stats <-
  true_sim_combined %>%
  select(sample_id, source, value = pct_return) %>%
  group_by(source, sample_id) %>%
  summarize(mean = mean(value),
            sd = sd(value),
            expon2 = mean(value^2),
            expon3 = mean(value^3),
            expon4 = mean(value^4),
            sd_p = sd(pmax(value, 0)),
            sd_m = sd(pmin(value, 0)),
            qdist_25_75 = quantile_dist(value, 0.25, 0.75),
            qdist_1_9 = quantile_dist(value, 0.1, 0.9),
            qdist_1_5 = quantile_dist(value, 0.1, 0.5),
            qdist_5_9 = quantile_dist(value, 0.5, 0.9),
            max = max(value),
            min = min(value),
            num_gt_5 = sum(value > 0.05),
            num_lt_5 = sum(value < -0.05),
            num_gt_10 = sum(value > 0.10),
            num_lt_10 = sum(value < -0.10),
            mode_count = mode_count(value)) %>%
  ungroup() %>%
  pivot_longer(c(-source, -sample_id))

test_stats %>%
  filter(source == "sim") %>%
  ggplot(aes(value)) +
  geom_histogram() +
  geom_vline(aes(xintercept = value), color = "red", data = filter(test_stats, source == "true")) +
  facet_wrap(. ~ name, scales = "free")

#################
# FIT GARCH SGT
#################
fit <-
  stan("garch_1_1_sgt.stan",
       data = list(T = length(r), r = r,
                   m_alpha0 = 0.0, s_alpha0 = 1e1,
                   m_p = 2, s_p = 1e2,
                   m_q = 2, s_q = 1e2),
       chains = 4, cores= 4, iter = 4e3, refresh = 100)

print(fit, pars = c("mu", "alpha0", "alpha1", "beta1", "sigma1", "l", "p", "q"))

s <- extract(fit)

# Generated quantities for GARCH SGT model since Stan doesn't have native SGT
sgt_garch_gen_quant <- function(s, r) {

  post_garch_sgt_par_draws <-
    tibble(mu = s$mu,
           alpha0 = s$alpha0, alpha1 = s$alpha1, beta1 = s$beta1, sigma1 = s$sigma1,
           lambda = s$l, p = s$p, q = s$q) %>%
    mutate(sample_id = row_number())
  
  # Function to get sigma_t and r_rep
  # Also sigma_hat which is what the sigma would be fit to the true time series
  get_sigma <- function(mu, alpha0, alpha1, beta1, sigma1, lambda, p, q, sample_id) {
    
    N <- length(r)
    sigma_rep <- rep(1.0, N)
    r_rep <- rep(0.0, N)
    
    sigma_hat <- rep(1.0, N)
    r_hat <- rep(0.0, N)
    
    sigma_rep[1] <- sigma1
    sigma_hat[1] <- sigma1
    r_rep[1] <- rsgt(1, mu, sigma_rep[1], lambda, p, q, mean.cent = FALSE, var.adj = FALSE)
    r_hat[1] <- rsgt(1, mu, sigma_hat[1], lambda, p, q, mean.cent = FALSE, var.adj = FALSE)
    
    for (t in 2:N) {
      sigma_rep[t] = sqrt(alpha0 + alpha1 * (r_rep[t-1] - mu)^2 + beta1 * sigma_rep[t-1]^2)
      r_rep[t] = rsgt(1, mu, sigma_rep[t], lambda, p, q, mean.cent = FALSE, var.adj = FALSE)
      
      sigma_hat[t] = sqrt(alpha0 + alpha1 * (r[t-1] - mu)^2 + beta1 * sigma_hat[t-1]^2)
      r_hat[t] = rsgt(1, mu, sigma_hat[t], lambda, p, q, mean.cent = FALSE, var.adj = FALSE)
    }
    
    tibble(sigma_rep, r_rep, sigma_hat, r_hat) %>%
      mutate(sample_id = sample_id) %>%
      mutate(time_id = row_number())
  }
  
  sgt_draws <-
    post_garch_sgt_par_draws %>%
    pmap(get_sigma) %>%
    bind_rows() %>%
    inner_join(post_garch_sgt_par_draws)
  
  sgt_draws
}

system.time(r_rep <- sgt_garch_gen_quant(s, r))

# Compare real data and simulate

# Time Series
rand_idx <- sample(1:500, 4)
r_rep %>%
  select(sample_id, time_id, r = r_rep) %>%
  mutate(source = "sim") %>%
  bind_rows(tibble(sample_id = as.integer(NA), time_id = 1:length(r), r = r, source = "true")) %>%
  filter(sample_id %in% rand_idx | source == "true") %>%
  ggplot(aes(time_id, r, color = source)) +
  geom_line() +
  facet_grid(sample_id ~ .)

# Histogram
rand_idx <- sample(1:500, 15)
r_rep %>%
  select(sample_id, time_id, r = r_rep) %>%
  mutate(source = "sim") %>%
  bind_rows(tibble(sample_id = as.integer(NA), time_id = 1:length(r), r = r, source = "true")) %>%
  filter(sample_id %in% rand_idx | source == "true") %>%
  ggplot(aes(r, fill = source)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.005, center = 0.0025, closed = "right") +
  facet_wrap(sample_id ~ .)

# Test Stats
test_stats <-
  r_rep %>%
  select(sample_id, time_id, value = r_rep) %>%
  mutate(source = "sim") %>%
  bind_rows(tibble(sample_id = as.integer(NA), time_id = 1:length(r), value = r, source = "true")) %>%
  group_by(source, sample_id) %>%
  summarize(mean = mean(value),
            sd = sd(value),
            expon2 = mean(value^2),
            expon3 = mean(value^3),
            expon4 = mean(value^4),
            sd_p = sd(pmax(value, 0)),
            sd_m = sd(pmin(value, 0)),
            qdist_25_75 = quantile_dist(value, 0.25, 0.75),
            qdist_1_9 = quantile_dist(value, 0.1, 0.9),
            qdist_1_5 = quantile_dist(value, 0.1, 0.5),
            qdist_5_9 = quantile_dist(value, 0.5, 0.9),
            max = max(value),
            min = min(value),
            num_gt_5 = sum(value > 0.05),
            num_lt_5 = sum(value < -0.05),
            num_gt_10 = sum(value > 0.10),
            num_lt_10 = sum(value < -0.10),
            mode_count = mode_count(value)) %>%
  ungroup() %>%
  pivot_longer(c(-source, -sample_id))

test_stats %>%
  filter(source == "sim") %>%
  ggplot(aes(value)) +
  geom_histogram() +
  geom_vline(aes(xintercept = value), color = "red", data = filter(test_stats, source == "true")) +
  facet_wrap(. ~ name, scales = "free")

# Predict true market 

# Compare real data and simulate

# Time Series
rand_idx <- sample(1:500, 4)
r_rep %>%
  select(sample_id, time_id, r = r_hat) %>%
  mutate(source = "sim") %>%
  bind_rows(tibble(sample_id = as.integer(NA), time_id = 1:length(r), r = r, source = "true")) %>%
  filter(sample_id %in% rand_idx | source == "true") %>%
  ggplot(aes(time_id, r, color = source)) +
  geom_line() +
  facet_grid(sample_id ~ .)

# Plot outliers
r_rep %>%
  group_by(time_id) %>%
  summarize(q10 = quantile(r_hat, 0.1), q50 = quantile(r_hat, 0.5), q90 = quantile(r_hat, 0.9)) %>%
  inner_join(mutate(sp500_weekly, time_id = row_number()) %>% select(time_id, pct_return)) %>%
  mutate(error = pct_return - q50,in_90_int = q10 < pct_return & pct_return < q90) %>%
  ggplot(aes(time_id, pct_return, color = in_90_int)) +
  geom_point()

r_rep %>%
  group_by(time_id) %>%
  summarize(q10 = quantile(r_hat, 0.1), q50 = quantile(r_hat, 0.5), q90 = quantile(r_hat, 0.9)) %>%
  inner_join(mutate(sp500_weekly, time_id = row_number()) %>% select(time_id, pct_return)) %>%
  mutate(error = pct_return - q50,in_90_int = q10 < pct_return & pct_return < q90) %>%
  ggplot(aes(time_id, pct_return, color = in_90_int)) +
  geom_point()

r_rep %>%
  inner_join(mutate(sp500_weekly, time_id = row_number()) %>% select(time_id, pct_return)) %>%
  group_by(time_id) %>%
  summarize(q10 = quantile(r_hat, 0.1), q50 = quantile(r_hat, 0.5), q90 = quantile(r_hat, 0.9),
            med_sigma_hat = median(sigma_hat), med_eps = median((pct_return - mu)^2)) %>%
  inner_join(mutate(sp500_weekly, time_id = row_number()) %>% select(time_id, pct_return)) %>%
  mutate(error = pct_return - q50,in_90_int = q10 < pct_return & pct_return < q90) %>%
  select(time_id, in_90_int, med_sigma_hat, med_eps, pct_return) %>%
  pivot_longer(c(-time_id, -in_90_int)) %>%
  ggplot(aes(time_id, value, color = in_90_int)) +
  geom_point() +
  facet_grid(name ~ ., scales = "free")

## Compare to Vix
vxx_weekly

r_rep %>%
  group_by(time_id) %>%
  summarize(q10 = quantile(r_hat, 0.1), q50 = quantile(r_hat, 0.5), q90 = quantile(r_hat, 0.9), s = sd(r_hat)) %>%
  inner_join(mutate(vxx_weekly, time_id = row_number())) %>%
  ggplot(aes(s, pct_return)) +
  geom_point()






