S <- 1e4
Q <- 188
s <- matrix(rnorm(S*Q), nrow = S)

s_sorted <- t(apply(s, 1, sort))

q <- t(apply(s_sorted, 2, quantile, probs = c(0.025, 0.5, 0.975)))

as_tibble(q)
  
diagnostics %>%
  filter(year(date) >= 2020, month(date) == 3) %>%
  mutate(z = rnorm(nrow(.))) %>%
  mutate(cum_z = cumsum(z^2)) %>% ggplot(aes(date, cum_z)) + geom_point()

diagnostics %>%
  filter(year(date) >= 2020, month(date) == 3) %>%
  arrange(z) %>%
  bind_cols(as_tibble(q)) %>%
  ggplot(aes(`50%`, z)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`))






t <- y_rep_holdout[1:5, 1:3]

y_rep_resampled <-
  t(y_rep_holdout) %>%
  apply(1, sample, size = 1e3, replace = TRUE) %>%
  t() %>%
  as_tibble() %>%
  mutate(date = test$date, y = test$ret_spy) %>%
  filter(year(date) >= 2010)

# y_rep_resampled %>%
#   tidyr::pivot_longer(-date) %>%
#   ggplot(aes(value)) +
#   geom_histogram() +
#   facet_wrap(name ~ .)
bins <- c(-100, -10, -5, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 5, 100)
true_bins <-
  y_rep_resampled %>%
  tidyr::pivot_longer(-date) %>%
  filter(name == "y") %>%
  mutate(bin = cut(value, bins)) %>%
  group_by(name, bin) %>%
  summarize(p = n()) %>%
  ungroup()

y_rep_resampled %>%
  tidyr::pivot_longer(-date) %>%
  filter(name != "y") %>%
  mutate(bin = cut(value, bins)) %>%
  group_by(name, bin) %>%
  summarize(p = n()) %>%
  ungroup() %>%
  group_by(bin) %>%
  summarize(q5 = quantile(p, 0.05),
            q50 = quantile(p, 0.5),
            q95 = quantile(p, 0.95)) %>%
  mutate(true_bins = true_bins$p) %>%
  mutate(bin_start = bins[1:(length(bins)-1)],
         bin_end = bins[-1]) %>%
  mutate(mid = (bin_end + bin_start) / 2) %>%
  ggplot() +
  geom_rect(aes(xmin = bin_start,
                xmax = bin_end,
                ymin = 0,
                ymax = q50), alpha = 0.5) +
  geom_errorbar(aes(mid, ymin = q5, ymax = q95)) +
  geom_point(aes(mid, true_bins), color = "red") +
  xlim(-10, 10)
