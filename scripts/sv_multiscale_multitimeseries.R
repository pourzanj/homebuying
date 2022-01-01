library(tidyverse)
library(tibbletime)
library(tidyquant)
library(cmdstanr)
library(posterior)
library(bayesplot)

tlt <-
  tq_get(c("TLT"),
         get = "stock.prices",
         from = "2018-01-01",
         to = "2021-12-16") %>%
  arrange(date) %>%
  mutate(r = 100 * log(close/lag(close))) %>%
  # if any returns are exactly zero this causes the lp to evaluate to zero
  mutate(r = ifelse(r == 0.0, 1e-4, r)) %>%
  na.omit() 

vvix <-
  tq_get(c("^VVIX"),
         get = "stock.prices",
         from = "2018-01-01",
         to = "2021-12-16") %>%
  arrange(date) %>%
  mutate(r = 100 * log(close/lag(close))) %>%
  # if any returns are exactly zero this causes the lp to evaluate to zero
  mutate(r = ifelse(r == 0.0, 1e-4, r)) %>%
  na.omit() 

vix <-
  tq_get(c("^VIX"),
         get = "stock.prices",
         from = "1990-01-01",
         to = "2021-12-20") %>%
  arrange(date) %>%
  mutate(r = 100 * log(close/lag(close))) %>%
  # if any returns are exactly zero this causes the lp to evaluate to zero
  mutate(r = ifelse(r == 0.0, 1e-4, r)) %>%
  na.omit() 

spx <-
  tq_get(c("^GSPC"),
         get = "stock.prices",
         from = "1990-01-01",
         to = "2021-12-20") %>%
  arrange(date) %>%
  mutate(r = 100 * log(close/lag(close)),
         hl = 100 * log(high/low),
         vol = 100 * log(volume / lag(volume))) %>%
  # if any returns are exactly zero this causes the lp to evaluate to zero
  mutate(r = ifelse(r == 0.0, 1e-4, r)) %>%
  na.omit() %>%
  mutate(vix = vix$r)


alizadeh <- function(r, s) {
  ret <- 0.0
  for(n in 1:1e2) {
    eps <- 8 * (-1)^(n-1) * n^2*r/(sqrt(2*pi)*s) * exp(-n^2*r^2 / (2*s^2))
    ret <- ret + eps
    if(abs(eps) <= 1e-7) break
  }
  
  ret
}

alizadeh(100*(log(4667)-log(4600)), 1e0)

molina <- function(a, b, y, y0, m, s, N = 100) {
  d1 <- function(n) {
    (y - y0 - 2*n*(b-a))^2 / (2*s^2)
  }
  
  d2 <- function(n) {
    (y + y0 - 2*a - 2*n*(b-a))^2 / (2*s^2)
  }
  
  k <- 1/(sqrt(2*pi)*s^3)*exp(-(m^2-2*m*(y-y0))/(2*s^2))
  
  tibble(n = 1:N) %>%
    mutate(neg = 4*n^2*(2*d1(-n)-1)*exp(-d1(-n)) - 4*-n*(-n-1)*(2*d2(-n)-1)*exp(-d2(-n)),
           pos = 4*n^2*(2*d1(n)-1)*exp(-d1(n)) - 4*n*(n-1)*(2*d2(n)-1)*exp(-d2(n))) %>%
    mutate(sum = neg + pos) %>%
    mutate(cumsum = cumsum(sum)) %>%
    mutate(p = k*cumsum) %>%
    pull(p) %>%
    tail(1)
}



molina <- function(a, b, y, y0, m, s) {
  
  d1 <- function(n) {
    (y - y0 - 2*n*(b-a))^2 / (2*s^2)
  }
  
  d2 <- function(n) {
    (y + y0 - 2*a - 2*n*(b-a))^2 / (2*s^2)
  }
  
  k <- 1/(sqrt(2*pi)*s^3)*exp(-(m^2-2*m*(y-y0))/(2*s^2))
  
  ret <- 0.0
  for(n in 1:1e2) {
    d1n <- d1(-n)
    d2n <- d2(-n)
    neg <- 4*n^2*(2*d1n-1)*exp(-d1n) - 4*-n*(-n-1)*(2*d2n-1)*exp(-d2n)
    
    d1p <- d1(n)
    d2p <- d2(n)
    pos <- 4*n^2*(2*d1p-1)*exp(-d1p) - 4*n*(n-1)*(2*d2p-1)*exp(-d2p)
    
    eps <- k*(neg+pos)
    
    ret <- ret + eps
    if(abs(eps) <= 1e-7) break
    if(n == 1e2) print(paste("eps:", eps))
  }
  
  ret
}

molina(a = 100 * log(4600), b = 100 * log(4667), y = 100 * log(4621),
       y0 = 100 * log(4652), m = 0, s = 1.1e0,
       N = 1e3) %>%
  ggplot(aes(n, p)) + geom_point()

a <- 100 * log(4600)
y0 <- 100 * log(4652)
b <- 100 * log(4667)

tibble(a = seq(100 * log(4500), 100 * log(4621), by = 0.01)) %>%
  mutate(p = map_dbl(a, function(a) molina(a = a,
                                           b = b,
                                           y = 100 * log(4621),
                                           y0 = y0,
                                           m = 0.0, s = 1.0))) %>%
  mutate(cdf = cumsum(p*0.01)) %>%
  ggplot(aes(exp(a/100), p)) +
  geom_line() +
  geom_vline(xintercept = exp(y0/100)) + 
  geom_vline(xintercept = exp(a/100)) +
  geom_vline(xintercept = exp(b/100))

molinag <- function(x) molina(x[1], x[3], x[2], y0, 0.0, 1.0)
grad(molinag, c(a,100 * log(4621),b))

grad(molinag, c(843.381,843.837,844.827))

# Try Moline CDF


# Test appendix

# Sample a, b, y in Stan
model <- cmdstan_model("stan/rmolina.stan")

fit <-
  model$sample(
    data = list(y0 = 100 * log(4652), m = 0.0, s = 1.0),
    seed = 1994,
    iter_warmup = 0e0,
    iter_sampling = 1e4,
    chains = 4,
    parallel_chains = 4,
    refresh = 1e2,
    metric = "unit_e",
    max_treedepth = 10,
    adapt_engaged = FALSE,
    adapt_delta = 0.8,
    step_size = 1e-1,
    show_messages = FALSE,
    init = replicate(4, list(a = 100 * log(4600), b = 100 * log(4667), y = 100 * log(4621)), simplify = FALSE)
  )



model <- cmdstan_model("stan/sv_multiscale_multitimeseries2.stan")

model <- cmdstan_model("stan/sv_multiscale_multitimeseries2_vov_clean.stan")

model <- cmdstan_model("stan/sv_multiscale_multitimeseries2_vov_clean_skew.stan")

model <- cmdstan_model("stan/sv_multiscale_multitimeseries2_vovov.stan")

# model <- cmdstan_model("stan/sv_multiscale_multitimeseries2_arg.stan")

# model <- cmdstan_model("stan/sv_multiscale_multitimeseries2_fouque.stan")


dat <- list(T = nrow(spx),
            y = spx %>% select(r, hl, vix) %>% mutate(hl = log(hl)) %>% as.matrix(),
            p0 = 1e-1, scale = sqrt(1.0), nu = 5, C2 = 1^2, p = 2)

fit <-
  model$sample(
    data = dat,
    seed = 1994,
    iter_warmup = 2e2,
    iter_sampling = 2e2,
    chains = 4,
    parallel_chains = 4,
    refresh = 1e0,
    max_treedepth = 10,
    adapt_delta = 0.8,
    show_messages = FALSE,
    init = replicate(4, list(mu_ret = 0.05, mu_hloffset = 0.25, mu = -3,
                             psi = 0.8, kappa = 0.2,
                             theta = 0.85, tau = 0.26, lambda =2.0, xi = -1, nu = 8,
                             phi = 0.99, sigma = 0.075,
                             sigma_obs = c(0.4, 7),
                             nu = 25,theta = 0.6, tau = 0.4, eta=0,
                             v_raw = rep(0.0, dat$T), u_raw = rep(0.0, dat$T), w_raw = rep(0.0, dat$T),
                             L = diag(5)), simplify = FALSE)
  )

fit$save_object("data/vovov_fit.RDS")

fit$draws(c("mu_ret", "mu_hloffset", "mu", "phi", "sigma", "R[1,4]","R[2,4]", "R[3,4]", "lp__")) %>% mcmc_pairs()

fit$summary("v_raw") %>% mutate(date= spx$date) %>% ggplot(aes(date, median))+ geom_point() + geom_errorbar(aes(ymin = q5, ymax = q95),alpha=0.9)

h <- fit$draws("h") %>% as_draws_df() %>% as_tibble() %>% mutate(date= spx$date)

h %>% select(-.draw) %>% pivot_longer(c(-.iteration, -.chain)) %>% mutate(i = as.integer(str_extract(name, "[0-9]+"))) %>% filter(.iteration == 55) %>% ggplot(aes(i, value, color = factor(.chain))) + geom_line()

key <-
  spx %>%
  select(date, r, hl, vix) %>%
  mutate(hl = log(hl)) %>%
  mutate(time_id = row_number()) %>%
  pivot_longer(c(-date, -time_id)) %>%
  inner_join(tibble(name = c("r", "hl", "vix"), var_id = c(1,2,3)))

y_rep <-
  fit$summary("y_rep") %>%
  mutate(time_id = as.integer(str_extract(variable, "(?<=\\[)[0-9]+")),
         var_id = as.integer(str_extract(variable, "(?<=,)[0-9]+")))

y_rep %>%
  inner_join(key) %>%
  ggplot(aes(date, median)) +
  geom_line() +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.5) +
  geom_line(aes(date, value), color = "red") +
  facet_grid(name ~ ., scales = "free")

f <- fit$summary("f") %>% mutate(date = spx$date)
h <- fit$summary("h") %>% mutate(date = spx$date)
g <- fit$summary("g") %>% mutate(date = spx$date)

f_rep <- fit$summary("f_rep") %>% mutate(date = spx$date) %>% mutate(name = "f", value = f$median)
h_rep <- fit$summary("h_rep") %>% mutate(date = spx$date) %>% mutate(name = "h", value = h$median)
g_rep <- fit$summary("g_rep") %>% mutate(date = spx$date) %>% mutate(name = "g", value = g$median)

w <- fit$summary("w") %>% mutate(date = spx$date)
u <- fit$summary("u") %>% mutate(date = spx$date)
v <- fit$summary("v") %>% mutate(date = spx$date)

w_rep <- fit$summary("w_rep") %>% mutate(date = spx$date) %>% mutate(name = "w", value = w$median)
u_rep <- fit$summary("u_rep") %>% mutate(date = spx$date) %>% mutate(name = "u", value = u$median)
v_rep <- fit$summary("v_rep") %>% mutate(date = spx$date) %>% mutate(name = "v", value = v$median)

w_raw <- fit$summary("w_raw") %>% mutate(date = spx$date)
u_raw <- fit$summary("u_raw") %>% mutate(date = spx$date)
v_raw <- fit$summary("v_raw") %>% mutate(date = spx$date)

w_raw_rep <- fit$summary("w_raw_rep") %>% mutate(date = spx$date) %>% mutate(name = "w_raw", value = w_raw$median)
u_raw_rep <- fit$summary("u_raw_rep") %>% mutate(date = spx$date) %>% mutate(name = "u_raw", value = u_raw$median)
v_raw_rep <- fit$summary("v_raw_rep") %>% mutate(date = spx$date) %>% mutate(name = "v_raw", value = v_raw$median)

y_rep %>%
  inner_join(key) %>%
  bind_rows(h_rep) %>%
  bind_rows(g_rep) %>%
  filter(year(date) == 2020) %>%
  ggplot(aes(date, median)) +
  geom_line() +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.5) +
  geom_line(aes(date, value), color = "red") +
  geom_point(aes(date, value), color = "red") +
  facet_grid(name ~ ., scales = "free")

y_rep %>%
  inner_join(key) %>%
  bind_rows(h_rep) %>%
  bind_rows(g_rep) %>%
  # bind_rows(f_rep) %>%
  # bind_rows(w_rep) %>%
  bind_rows(u_rep) %>%
  bind_rows(v_rep) %>%
  # bind_rows(w_raw_rep) %>%
  bind_rows(u_raw_rep) %>%
  bind_rows(v_raw_rep) %>%
  # bind_rows(select(spx, date, median = adjusted, value = adjusted) %>% mutate(name = "SPX")) %>% 
  # bind_rows(select(spx %>% mutate(drawdown = 100*log(adjusted/cummax(adjusted))), date, median = drawdown, value = drawdown) %>% mutate(name = "drawdown")) %>% 
  # bind_rows(select(tlt, date, median = r, value = r) %>% mutate(name = "TLT")) %>% 
  # bind_rows(select(spx, date, median = vol, value = vol) %>% mutate(name = "vol")) %>% 
  # inner_join(select(k, date, k)) %>%
  filter(year(date) == 2021) %>%
  mutate(o5 = value <= q5, o95 = value >= q95, o = o5 | o95) %>%
  ggplot(aes(date, median)) +
  geom_line() +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.2) +
  geom_errorbar(aes(date, median, ymin = q5, ymax = q95), data = u_raw %>% mutate(name = "u_raw") %>% filter(year(date) == 2021)) +
  geom_line(aes(date, value), color = "red") +
  geom_point(aes(date, value, color = o)) +
  facet_grid(name ~ ., scales = "free")


y_rep_draws <- fit$draws("y_rep[553,1]") %>% as_draws_df()
y_rep_draws %>% ggplot(aes(`y_rep[553,1]`)) + geom_histogram()



fit$summary("u_rep") %>% mutate(date= spx$date) %>% mutate(var= "u_rep") -> u_rep
fit$summary("u") %>% mutate(date= spx$date) %>% mutate(var= "u") -> u

bind_rows(u, u_rep) %>% filter(year(date) == 2021) %>% ggplot(aes(date, median,color = var, fill = var)) + geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.5) + geom_line() + geom_point()

# Generate draws of y
draws <-
  fit$draws(c("tau", "theta", "sigma", "g[530]",
              "R[4,5]", "R22_R23R32", "phi", "h[530]",
              "mu", "mu_ret", "R[1,5]","R12_R13R32","R_e_hg[1,1]")) %>%
  as_draws_df() %>%
  as_tibble() %>%
  sample_n(1)

N <- 1e5
tau <- draws$tau
theta <- draws$theta
sigma <- draws$sigma
g0 <- draws$`g[530]`
R45 <- draws$`R[4,5]`
R22_R23R32 <- draws$R22_R23R32
phi <- draws$phi
h0 <- draws$`h[530]`
mu <- draws$mu
mu_ret1 <- draws$`mu_ret[1]`
mu_ret2 <- draws$`mu_ret[2]`
mu_ret3 <- draws$`mu_ret[3]`
R15 <- draws$`R[1,5]`
R12_R13R32 <- draws$`R12_R13R32[1]`
R_e_hg <- draws$`R_e_hg[1,1]`


sim_y <-
  tibble(
    l_raw = abs(rcauchy(N)),
    u_raw = rnorm(N),
       u = tau*u_raw,
       g = theta*g0 + u,
       v_raw = rnorm(N),
       v = sigma*exp(g/2)* (R45*u_raw + sqrt(R22_R23R32)*v_raw),
       h = phi*h0 + v,
       s = exp((mu+h)/2),
       m = mu_ret1 + mu_ret2*s + mu_ret3*s^2 + s*(R15*u_raw + R12_R13R32*v_raw),
       y = rnorm(N, m, s*sqrt(R_e_hg)))

qplot(sim_y$y, binwidth = 0.1)


# outside test
get_next_outside <- function(df) {
  N <- nrow(df)
  nexto95 <- rep(0L, N)
  
  for(i in 1:N) {
    if(df$o5[i]) {
      ctr <- 0L
      for(j in i:N) {
        if(df$o95[j]) {
          nexto95[i] <- ctr
          break
        } else {
          ctr <- ctr +1
        }
      }
    }  
  }
  
  mutate(df, nexto95)
}

get_next_outside(o) %>% filter(o5) %>% pull(nexto95) %>% median()

sim_o <- function(i) {
  select(o, date) %>%
    mutate(u = runif(nrow(.))) %>%
    mutate(o5 = u <= 0.05, o95 = u >= 0.95) %>%
    mutate(i = i)
}

osims <- map(1:1000, sim_o)
nextos <- map(osims, get_next_outside)
med_nexto95_sims <- map_dbl(nextos, function(df) df %>% filter(o5) %>% pull(nexto95) %>% median())


####
# Fit g by itself
g <- fit$summary("g") %>% mutate(date= spx$date)

g %>% ggplot(aes(date, median)) + geom_point() + geom_line()


model <- cmdstan_model("stan/ar1.stan")

dat <- list(N = nrow(g), y = g$median, tau0 = 1e-2, scale = 0.8, nu = 5)

fit_g <-
  model$sample(
    data = dat,
    seed = 1994,
    iter_warmup = 2e2,
    iter_sampling = 2e2,
    chains = 4,
    parallel_chains = 4,
    refresh = 1e0,
    max_treedepth = 12,
    adapt_delta = 0.99,
    show_messages = FALSE,
    init = replicate(4, list(tau=1e-2, C2 = 0.1, theta = 0.1), simplify = FALSE)
  )

y_rep <- fit_g$summary("y_rep")

y_rep %>%
  mutate(date = g$date, y = g$median) %>%
  filter(year(date) == 2021, month(date) >= 6) %>%
  ggplot(aes(date, median)) +
  geom_point() +
  geom_line() +
  geom_point(aes(date, y), color = "red") +
  geom_errorbar(aes(ymin = q5, ymax = q95), alpha = 0.5)

y_rep %>%
  mutate(date = g$date, y = g$median) %>%
  mutate(lt = y <= q5,
         gt = y >= q95) %>%
  na.omit() %>%
  summarize(n = n(), lt = sum(lt), gt = sum(gt)) %>%
  mutate(plt = pbinom(lt, n, 0.05), pgt = pbinom(gt, n, 0.05))

lt <- fit_g$summary("lt")
lt %>%
  mutate(date = g$date, y = g$median) %>%
  ggplot(aes(date, median)) +
  geom_point() +
  geom_errorbar(aes(ymin = q5, ymax = q95), alpha = 0.5)

fit_g$draws("y_rep[538]") %>% as_draws_df -> y_rep538
y_rep538$`y_rep[538]` %>% qplot()


# Sim the g process from VOV with skew
N <- 996
nu <- 6
# tau <- 1e3
tau <- 0.18
lambda <- 0.73
tibble(c2 = (1/rgamma(N, nu/2, rate = nu/2)),
       l = abs(rcauchy(N))) %>%
  # mutate(lt = tau * sqrt((c2*l^2) / (c2 + tau^2 * l^2))) %>%
  mutate(lt = c2) %>%
  mutate(T0 = abs(rnorm(N)),
         T1 = rnorm(N),
         delta = lambda*lt / sqrt(1+lambda^2 * lt^2),
         u_raw = -0.78 + lt*(delta*T0 + sqrt(1-delta^2)*T1),
         u = tau*u_raw) %>%
  mutate(g_sim = Reduce(function(g0, u1) 0.9*g0 + u1, u, accumulate = TRUE)) %>%
  mutate(date = spx$date) %>%
  mutate(g_true = g$median) %>%
  select(date, g_sim, g_true) %>%
  pivot_longer(-date) %>%
  filter(year(date) == 2021) %>%
  ggplot(aes(date, value, color = name)) + 
  geom_line() + geom_point()

# Simulate affine SDE
dt <- 0.25
N <- 250*(1/dt)
g0 <- 0.3

df <- 
  tibble(t = dt*0:N,
         u_raw = c(g0,rnorm(N))) %>%
  mutate(g = accumulate(u_raw,
                        function(g0, u1) g0 + 0.1*(0.35-g0)*dt + sqrt(g0)*0.2*u1*sqrt(dt)))

df %>% filter(t %% 1 == 0) %>%
  ggplot(aes(t, g)) +
  geom_line() +
  geom_point()

df %>% filter(t %% 1 == 0) %>%
  ggplot(aes(t, g-lag(g))) +
  geom_line() +
  geom_point()

# Compare how step size affects daily distribution by simulate daily distribution with different step size
# For lower g0 close to zero e.g. 0.1 there is heavy right tail, but for higher g0 close to e.g. 0.8 the
# distribution is more symmetric with a mode below g0 (wanting to revert back to mean value). Note the
# skew doesn't appear if dt = 1.0.
N <- 2e4

dt <- 0.1
H <- (1/dt)

g0 <- 0.1
df_g0_01 <- 
  tibble(i = rep(1:N, each = H+1),
         t = rep(dt*(0:H), N),
         u = rnorm(N*(H+1))) %>%
  group_by(i) %>%
  summarize(g = reduce(u[2:(H+1)], function(g0, u1) g0 + 0.1*(0.35-g0)*dt + (g0)^(0.5)*0.2*u1*sqrt(dt), .init = g0)) %>%
  ungroup() %>%
  mutate(g0 = "0.1")

g0 <- 0.5
df_g0_05 <- 
  tibble(i = rep(1:N, each = H+1),
         t = rep(dt*(0:H), N),
         u = rnorm(N*(H+1))) %>%
  group_by(i) %>%
  summarize(g = reduce(u[2:(H+1)], function(g0, u1) g0 + 0.1*(0.35-g0)*dt + (g0)^(0.5)*0.2*u1*sqrt(dt), .init = g0)) %>%
  ungroup() %>%
  mutate(g0 = "0.5")

g0 <- 1.0
df_g0_10 <- 
  tibble(i = rep(1:N, each = H+1),
         t = rep(dt*(0:H), N),
         u = rnorm(N*(H+1))) %>%
  group_by(i) %>%
  summarize(g = reduce(u[2:(H+1)], function(g0, u1) g0 + 0.1*(0.35-g0)*dt + (g0)^(0.5)*0.2*u1*sqrt(dt), .init = g0)) %>%
  ungroup() %>%
  mutate(g0 = "1.0")

df_g0_01 %>%
  bind_rows(df_g0_05) %>%
  bind_rows(df_g0_10) %>%
  ggplot(aes(g)) +
  geom_histogram(binwidth = 0.025) +
  facet_grid(g0 ~ .)

# Sim VoV
dt <- 0.02
N <- 999*(1/dt)
V0 <- 0.3
logVIX0 <- 2.21

k <- 0.01
theta <- 3.0

kv <- 0.1
thetav <- 0.35
sigmav <- 0.18
rho <- 0.653

df <- 
  tibble(t = dt*0:N,
         Wv = c(V0,rnorm(N)),
         W = c(logVIX0, rho*Wv[2:(N+1)] + sqrt(1-rho^2)*rnorm(N))) %>%
  mutate(V = accumulate(Wv, function(V00, Wv1) V00 + kv*(thetav-V00)*dt + sigmav*sqrt(V00)*Wv1*sqrt(dt))) %>%
  mutate(Wdt = c(logVIX0,sqrt(V[1:N]/100)*W[2:(N+1)]*sqrt(dt))) %>%
  mutate(logVIX = accumulate(Wdt, function(logVIX00, Wdt1) logVIX00 + k*(theta-logVIX00)*dt + Wdt1))

df %>%
  filter(t %% 1 == 0) %>%
  select(t, V, logVIX) %>%
  mutate(logVIX_true = log(vix$close) %>% head(1e3)) %>%
  pivot_longer(-t) %>%
  ggplot(aes(t, value)) +
  geom_line() +
  # geom_point() +
  facet_grid(name ~ ., scales = "free")

df %>%
  filter(t %% 1 == 0) %>%
  mutate(logVIX_true = log(vix$close)) %>%
  mutate(r = logVIX_true - lag(logVIX_true)) %>%
  mutate(r_sim = logVIX - lag(logVIX)) %>%
  select(t, r, r_sim) %>%
  pivot_longer(-t) %>%
  ggplot(aes(value)) +
  geom_histogram() +
  facet_grid(name ~ ., scales = "free")


# Generate from closed from non-central chisquare of CIR
rcir <- function(n, r0, a, b, sigma) {
  c <- 2*a / ((1-exp(-a))*sigma^2)
  rchisq(n, 4*a*b / sigma^2, 2*c*r0*exp(-a)) / (2*c)
}











