vix <-
  tq_get(c("^VIX"),
         get = "stock.prices",
         from = "2002-01-01",
         to = "2022-12-31") %>%
  select(date, vix = close) %>%
  mutate(vixd = vix/sqrt(252), vixd2 = vixd^2, h = log(vixd2)) %>%
  mutate(v = h-lag(h)) %>%
  na.omit()

model <- cmdstan_model("stan/vix_innovations.stan")

fit <-
  model$sample(
    data = list(T = nrow(vix), v = vix$v, h = vix$h,
                m_p = 2, s_p = 1,
                m_q = 15, s_q = 5),
    seed = 1994,
    iter_warmup = 5e2,
    iter_sampling = 5e2,
    chains = 4,
    parallel_chains = 4,
    refresh = 10,
    max_treedepth = 10,
    adapt_delta = 0.8
  )

N <- nrow(vix)
v_rep <- rep(0, N)
h_rep <- rep(0.592, N)
for(n in 2:N) {
  v_rep[n] <- rsgt(1, 0.00465 + -0.2*h_rep[n-1], exp(-2+0.188*h_rep[n-1]), 2*ilog(0.413+-0.115*h_rep[n-1])-1, 1.68, 3.18)
  h_rep[n] <- h_rep[n-1] + v_rep[n]
}

vix %>%
  mutate(v_rep = v_rep, h_rep = h_rep) %>%
  ggplot(aes(date, h_rep)) +
  geom_line()





