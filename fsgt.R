# Fernandex and Steel Generalized T density
dfsgt <- function(x, v, nu, gamma, log = FALSE) {
  term1 <- log(2) - log(gamma + 1/gamma)
  term2 <- log(v) - log(2*nu^(1/v))
  term3 <- -lbeta(nu/2, 1/v)
  main_term <- -(nu/2 + 1/v) * log(1 + (1/nu) * abs(x)^v / gamma^(v*sign(x)))
  
  log_p <- term1 + term2 + term3 + main_term
  
  if(log == TRUE) {
    return(log_p)
  } else {
    return(exp(log_p))
  }
}

tibble(x = seq(-3, 3, by = 0.02)) %>% mutate(p = dfsgt(x, v = 1.0, nu = 100, gamma = 1.0)) %>% ggplot(aes(x, p)) + geom_line()

# Expected value
expec_fsgt <- function(v, nu, gamma) {
  (gamma - 1/gamma) * ((gamma(2/v)*gamma(nu/2-1/v)) / (gamma(1/v)*gamma(nu/2))) * nu^(1/v)
}

expec_fsgt(v = 2, nu = 2, gamma = 1.3)
integrate(function(x, v, nu, gamma) x*dfsgt(x, v, nu, gamma), -Inf, Inf, v = 2, nu = 2, gamma = 1.3)




# f adds mean and scale
f <- function(x, mu, lambda, v, nu, gamma, log = FALSE) {
  eps <- (x - mu) / exp(lambda) + expec_fsgt(v, nu, gamma)
  log_p <- -lambda + dfsgt(eps, v, nu, gamma, log = TRUE)
  
  if(log == TRUE) {
    return(log_p)
  } else {
    return(exp(log_p))
  }
}

tibble(x = seq(-3, 3, by = 0.1)) %>% mutate(p = f(x, mu = 0.5, lambda = 0.3, v = 2, nu = 10, gamma = 1.0)) %>% ggplot(aes(x, p)) + geom_line()

# CDF
pf <- function(x, mu, lambda, v, nu, gamma) {
  integrate(dfsgt, -Inf, x, mu = mu, lambda = lambda, v = v, nu = nu, gamma = gamma, log = FALSE)
}

# Random sample using inverse CDF method
rf <- function(x, mu, lambda, v, nu, gamma) {
  u <- runif(1)
  uniroot(function(x) pf(x) - u,
          c(-Inf, Inf),
          mu = mu, lambda = lambda, v = v, nu = nu, gamma = gamma, log = FALSE)
}

# Partial w.r.t. lambda
df_dlambda <- function(x, mu, lambda, v, nu, gamma) {
  eps <- (x - mu) / exp(lambda) + expec_fsgt(v, nu, gamma)
  mu_eps <- expec_fsgt(v, nu, gamma)
  (nu*v/2 + 1) * (1 - mu_eps / eps) * (abs(eps)^v / (abs(eps)^v + nu*gamma^(v*sign(eps)))) - 1
}

df_dlambda(0.2, mu = 0.5, lambda = 0.3, v = 2, nu = 10, gamma = 1.0)

e <- 1e-4
(f(0.2, mu = 0.5, lambda = 0.3 + e, v = 2, nu = 10, gamma = 1.0, log = TRUE) -
    f(0.2, mu = 0.5, lambda = 0.3, v = 2, nu = 10, gamma = 1.0, log = TRUE)) / e

# Test Stan fit
fit <-
  stan(file = "stan/harvey/fsgt.stan",
       data = list(N = length(spx$pct_return), y = spx$pct_return,
                   m_mu = 0, s_mu = 1,
                   m_lambda = 0, s_lambda = 5,
                   m_v = 1.5, s_v = 1,
                   m_nu = 5, s_nu = 2,
                   m_gamma = 1, s_gamma = 1),
       chains = 1, iter = 1e3)

s <- extract(fit)

samples <-
  tibble(mu = s$mu, lambda = s$lambda, v = s$v, nu = s$nu, gamma = s$gamma) %>%
  mutate(sample_id = row_number()) %>%
  sample_n(5)

samples %>%  
  crossing(x = seq(-12, 12, by = 0.025)) %>%
  select(-sample_id) %>%
  mutate(p = pmap_dbl(., f)) %>%
  inner_join(samples) %>%
  ggplot() +
  geom_histogram(aes(pct_return, y = ..density..), binwidth = 0.1, data = spx) +
  geom_line(aes(x, p, group = sample_id))

# Harvey Model 1A
generate_model_1a <- function(i,
                                  omega, phi, k1, k2, lambda1,
                                  mu, v, nu, gamma,
                                  date) {
  # Initialize output vectors
  N <- length(date)
  eps <- rep(0.0, N)
  u_t <- rep(0.0, N)
  h <- rep(0.0, N)
  lambda <- rep(0.0, N)
  Eeps <- rep(0.0, N)
  y <- rep(0.0, N)
  
  lambda[1] <- lambda1
  Eeps[1] <- expec_fsgt(v, nu, gamma)
  eps[1] <- rt(1, nu)
  y[1] <- (eps[1]-Eeps[1])*exp(lambda[1]) + mu
  u_t[1] = df_dlambda(y[1], mu, lambda[1], v, nu, gamma)
  h[1] = k1*u_t[1] + k2*sign(-eps[1])*(u_t[1]+1)
  
  for (t in 2:N) {
    # Set current scale from last scale
    lambda[t] <- omega*(1-phi) + phi*lambda[t-1] + h[t-1]
    
    # Sim current
    Eeps[t] <- expec_fsgt(v, nu, gamma)
    eps[t] <- rt(1, nu)
    y[t] <- (eps[t]-Eeps[t])*exp(lambda[t]) + mu
    
    # Setup for the next iteration
    u_t[t] <- df_dlambda(y[t], mu, lambda[t], v, nu, gamma)
    h[t] <- k1*u_t[t] + k2*sign(-eps[t])*(u_t[t]+1)
  }
  
  tibble(sample_id = i, date, eps = eps, u_t = u_t, h, lambda, y)
}

sim_1a <- generate_model_1a(i = 1,
                            omega = -0.3, phi = 0.96, k1 = 0.038, k2 = 0.02, lambda1 = -1,
                            mu = 0, v = 2, nu = 10, gamma = 1,
                            date = spx$date[1:4000])

spx[1:4000] %>% ggplot(aes(date, pct_return))  + geom_line()
sim_1a %>% ggplot(aes(date, y))  + geom_line()

init <-
  list(list(omega = -0.3, phi = 0.96, k1 = 0.038, k2 = 0.02, lambda1 = -1, mu = 0))

fit <-
  stan(file = "stan/harvey/model_1a.stan",
       data = list(T = nrow(sim_1a[1:1000,]), y = sim_1a$y[1:1000],
                   m_omega = 0, s_omega = 0.5,
                   a_phi = 20, b_phi = 2,
                   m_k1 = 0, s_k1 = 0.1,
                   m_k2 = 0, s_k2 = 0.1,
                   m_lambda1 = 0, s_lambda1 = 1,
                   m_mu = 0, s_mu = 0.5,
                   m_v = 2, s_v = 1,
                   m_nu = 10, s_nu = 5,
                   m_gamma = 1, s_gamma = 1),
       control = list(adapt_delta = 0.9, max_treedepth = 10),
       chains = 1, iter = 1e3, init = init, refresh = 10)

print(fit, pars = c("omega", "phi", "k1", "k2", "lambda1", "mu", "v", "nu", "gamma"))

s <- extract(fit)
