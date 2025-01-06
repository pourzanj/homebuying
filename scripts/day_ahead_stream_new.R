library(tidyverse)
library(tibbletime)
library(tidyquant)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(moments)
library(sgt)
library(loo)
library(gridExtra)

kOutdir <- "data/backtest_sv13_1_1992_2022/"
dir.create(kOutdir)

spx <-
  tq_get(c("^GSPC"),
         get = "stock.prices",
         from = "1992-01-01",
         to = "2022-12-31") %>%
  mutate(y = 100*log(close/lag(close))) %>%
  na.omit() %>%
  mutate(week = isoweek(date)) %>%
  mutate(year = if_else(month(date) == 12 & week == 1, year(date)+1, year(date))) %>%
  mutate(day = wday(date, label = TRUE))

#############
# Fit first
#############

model <- cmdstan_model("stan/sv13_1.stan")

M <- 100
L <- 252*10
past <- 1:L
oos <- L + 1
df_past <- spx[past, , drop = FALSE]
df_oos <- spx[oos:min(oos+M-1, nrow(spx)), , drop = FALSE]

init <-
  replicate(8, list(mu = 0.2, phi = 0.98, sigma = 0.17,
                    alpha = c(0.1, 0, 0),
                    gamma = c(0, 0, 0),
                    beta = c(0.7, 0),
                    q = 10,
                    v0 = 0,
                    v = rep(0.0, nrow(df_past)),
                    v_oos = rep(0.0, nrow(df_oos))), simplify = FALSE)

fit_past <-
  model$sample(
    data = list(T = nrow(df_past), y = df_past$y,
                M = nrow(df_oos), y_oos = df_oos$y,
                m_p = 2, s_p = 2,
                m_q = 15, s_q = 6),
    seed = 1994,
    iter_warmup = 200,
    iter_sampling = 200,
    chains = 8,
    parallel_chains = 8,
    refresh = 200,
    max_treedepth = 8,
    adapt_delta = 0.8,
    init = init,
    show_messages = FALSE
  )

#############
# LFO
#############

# some helper functions we'll use throughout

# more stable than log(sum(exp(x))) 
log_sum_exp <- function(x) {
  max_x <- max(x)  
  max_x + log(sum(exp(x - max_x)))
}

# more stable than log(mean(exp(x)))
log_mean_exp <- function(x) {
  log_sum_exp(x) - log(length(x))
}

# compute log of raw importance ratios
# sums over observations *not* over posterior samples
sum_log_ratios <- function(loglik, ids = NULL) {
  if (!is.null(ids)) loglik <- loglik[, ids, drop = FALSE]
  rowSums(loglik)
}

# for printing comparisons later
rbind_print <- function(...) {
  round(rbind(...), digits = 2)
}

# initialize the process for i = L
k_thres <- 0.7
y_oos <- matrix(NA, ncol = nrow(spx), nrow = 8*200)
loglik <- fit_past$draws("loglik_oos", format = "draws_matrix")

# Draw y_oos using original past data
y_oos0 <-
  fit_past$draws(c("m_oos", "s_oos", "l_oos", "p_oos", "q")) %>%
  as_draws_df() %>%
  as_tibble() %>%
  select(-.chain, -.iteration) %>%
  pivot_longer(c(-.draw, -q)) %>%
  mutate(i = as.integer(parse_number(name)),
         name = str_extract(name, "[a-z_]+")) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  mutate(y_oos= rsgt(nrow(.), m_oos, s_oos, l_oos, p_oos, q)) %>%
  select(.draw, i, y_oos) %>%
  pivot_wider(names_from = i, values_from = y_oos) %>%
  select(-.draw) %>%
  as.data.frame() %>%
  as.matrix()

y_oos[,L+1] <- y_oos0[,1]
# qplot(y_oos[,L+1]) # plot just to have a plot to turn dev.off() later

# Initialize plotting data
daily_summary <-
  spx %>%
  select(date, y) %>%
  mutate(r = exp(y/100)-1) %>%
  mutate(mean = NA_real_,
         sd = NA_real_,
         skew = NA_real_,
         kurtosis = NA_real_,
         E_kn3 = NA_real_,
         E_kn1 = NA_real_,
         E_k0 = NA_real_,
         E_kp1 = NA_real_,
         E_kp3 = NA_real_,
         k_vs = NA_real_,
         k_vss = NA_real_,
         r_1x = NA_real_,
         r_3x = NA_real_,
         r_vs = NA_real_,
         r_vss = NA_real_,
         R_1x = NA_real_,
         R_3x = NA_real_,
         R_vs = NA_real_,
         R_vss = NA_real_)

get_annualized_expected_leveraged_return <- function(y, k) {
  r <- exp(y/100)-1
  E_logr <- -Inf
  if(all(k*r > -1)) E_logr <- mean(log(1 + k*r))
  100 * (exp(E_logr)^252 - 1)
}

get_optimal_leverage <- function(E_r, k_levels) {
  max_idx <- unique(which(E_r == max(E_r)))
  k_levels[max_idx]
}

daily_summary$mean[L+1] <- mean(y_oos[,L+1])
daily_summary$sd[L+1] <- sd(y_oos[,L+1])
daily_summary$skew[L+1] <- mean((y_oos[,L+1]-mean(y_oos[,L+1]))^3) / sd(y_oos[,L+1])^3
daily_summary$kurtosis[L+1] <- mean((y_oos[,L+1]-mean(y_oos[,L+1]))^4) / sd(y_oos[,L+1])^4

daily_summary$E_kn3[L+1] <- get_annualized_expected_leveraged_return(y_oos[,L+1], -3)
daily_summary$E_kn1[L+1] <- get_annualized_expected_leveraged_return(y_oos[,L+1], -1)
daily_summary$E_k0[L+1] <- get_annualized_expected_leveraged_return(y_oos[,L+1], 0)
daily_summary$E_kp1[L+1] <- get_annualized_expected_leveraged_return(y_oos[,L+1], 1)
daily_summary$E_kp3[L+1] <- get_annualized_expected_leveraged_return(y_oos[,L+1], 3)


daily_summary$k_vs[L+1] <-
  get_optimal_leverage(c(daily_summary$E_k0[L+1],
                         daily_summary$E_kp1[L+1],
                         daily_summary$E_kp3[L+1]),
                       c(0, 1, 3))

daily_summary$k_vss[L+1] <-
  get_optimal_leverage(c(daily_summary$E_kn3[L+1],
                         daily_summary$E_kn1[L+1],
                         daily_summary$E_k0[L+1],
                         daily_summary$E_kp1[L+1],
                         daily_summary$E_kp3[L+1]),
                       c(-3, -1, 0, 1, 3))

r <- daily_summary$r[L+1]
daily_summary$r_1x[L+1] <- 1*r
daily_summary$r_3x[L+1] <- 3*r
daily_summary$r_vs[L+1] <- daily_summary$k_vs[L+1]*r
daily_summary$r_vss[L+1] <- daily_summary$k_vss[L+1]*r

daily_summary$R_1x[L+1] <- 1
daily_summary$R_3x[L+1] <- 1
daily_summary$R_vs[L+1] <- 1
daily_summary$R_vss[L+1] <- 1

# Initialize model parameter data
par_summary <- 
  fit_past$summary(c("alpha", "beta", "mu", "sigma", "phi", "gamma", "q")) %>%
  mutate(date = daily_summary$date[L+1])

# iterate over i > L
i_refit <- L
refits <- L
ks <- NULL

# y_oos <- readRDS(paste0(kOutdir, "y_oos.RDS"))
# daily_summary <- readRDS(paste0(kOutdir, "daily_summary.RDS"))
# par_summary <- readRDS(paste0(kOutdir, "par_summary.RDS"))

# i tracks absolute index (row) in spx. in each iteration we add the ith data
# point into the dataset to predict the i+1
for (i in (L + 1):(nrow(spx) - 1)) {
  
  logratio <- sum_log_ratios(loglik, 1:(i-i_refit))
  psis_obj <- suppressWarnings(psis(logratio))
  k <- pareto_k_values(psis_obj)
  ks <- c(ks, k)
  
  if (k > k_thres | (i-i_refit) > (M-2)) {
    
    # We have to refit model
    
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(paste("Refitting", i-(L+1)+1, "/", (nrow(spx)-1)-(L+1)+1))
    
    i_refit <- i
    refits <- c(refits, i)
    
    # Only take last L days and fit
    past <- (i-L+1):i
    oos <- i + 1
    df_past <- spx[past, , drop = FALSE]
    df_oos <- spx[oos:min(oos+M-1, nrow(spx)), , drop = FALSE]
    
    fit_past <-
      model$sample(
        data = list(T = nrow(df_past), y = df_past$y,
                    M = nrow(df_oos), y_oos = df_oos$y,
                    m_p = 2, s_p = 2,
                    m_q = 15, s_q = 6),
        seed = 1994,
        iter_warmup = 200,
        iter_sampling = 200,
        chains = 8,
        parallel_chains = 8,
        refresh = 200,
        max_treedepth = 8,
        adapt_delta = 0.8,
        init = init,
        show_messages = FALSE
      )
    
    # Get loglik
    loglik <- fit_past$draws("loglik_oos", format = "draws_matrix")
    
    # Draw y_oos using original past data
    y_oos0 <-
      fit_past$draws(c("m_oos", "s_oos", "l_oos", "p_oos", "q")) %>%
      as_draws_df() %>%
      as_tibble() %>%
      select(-.chain, -.iteration) %>%
      pivot_longer(c(-.draw, -q)) %>%
      mutate(i = as.integer(parse_number(name)),
             name = str_extract(name, "[a-z_]+")) %>%
      pivot_wider(names_from = name, values_from = value) %>%
      mutate(y_oos= rsgt(nrow(.), m_oos, s_oos, l_oos, p_oos, q)) %>%
      select(.draw, i, y_oos) %>%
      pivot_wider(names_from = i, values_from = y_oos) %>%
      select(-.draw) %>%
      as.data.frame() %>%
      as.matrix()
    
    # For fresh fit can just use posterior samples for first look ahead
    y_oos[,i+1] <- y_oos0[,1]
    
    #############
    # Update
    #############
    
    daily_summary$mean[i+1] <- mean(y_oos[,i+1])
    daily_summary$sd[i+1] <- sd(y_oos[,i+1])
    daily_summary$skew[i+1] <- mean((y_oos[,i+1]-mean(y_oos[,i+1]))^3) / sd(y_oos[,i+1])^3
    daily_summary$kurtosis[i+1] <- mean((y_oos[,i+1]-mean(y_oos[,i+1]))^4) / sd(y_oos[,i+1])^4
    
    daily_summary$E_kn3[i+1] <- get_annualized_expected_leveraged_return(y_oos[,i+1], -3)
    daily_summary$E_kn1[i+1] <- get_annualized_expected_leveraged_return(y_oos[,i+1], -1)
    daily_summary$E_k0[i+1] <- get_annualized_expected_leveraged_return(y_oos[,i+1], 0)
    daily_summary$E_kp1[i+1] <- get_annualized_expected_leveraged_return(y_oos[,i+1], 1)
    daily_summary$E_kp3[i+1] <- get_annualized_expected_leveraged_return(y_oos[,i+1], 3)
    
    
    daily_summary$k_vs[i+1] <-
      get_optimal_leverage(c(daily_summary$E_k0[i+1],
                             daily_summary$E_kp1[i+1],
                             daily_summary$E_kp3[i+1]),
                           c(0, 1, 3))
    
    daily_summary$k_vss[i+1] <-
      get_optimal_leverage(c(daily_summary$E_kn3[i+1],
                             daily_summary$E_kn1[i+1],
                             daily_summary$E_k0[i+1],
                             daily_summary$E_kp1[i+1],
                             daily_summary$E_kp3[i+1]),
                           c(-3, -1, 0, 1, 3))
    
    r <- daily_summary$r[i+1]
    daily_summary$r_1x[i+1] <- 1*r
    daily_summary$r_3x[i+1] <- 3*r
    daily_summary$r_vs[i+1] <- daily_summary$k_vs[i+1]*r
    daily_summary$r_vss[i+1] <- daily_summary$k_vss[i+1]*r
    
    daily_summary$R_1x[i+1] <- prod(1+daily_summary$r_1x[(L+1):(i+1)])
    daily_summary$R_3x[i+1] <- prod(1+daily_summary$r_3x[(L+1):(i+1)])
    daily_summary$R_vs[i+1] <- prod(1+daily_summary$r_vs[(L+1):(i+1)])
    daily_summary$R_vss[i+1] <- prod(1+daily_summary$r_vss[(L+1):(i+1)])
    
    # Update parameter summary
    par_summary_new <- 
      fit_past$summary(c("alpha", "beta", "mu", "sigma", "phi", "gamma", "q")) %>%
      mutate(date = daily_summary$date[i])
    
    par_summary <- bind_rows(par_summary, par_summary_new)
    
    #############
    # Plot
    #############
    
    # dev.off()
    
    # 1A Performance of methods
    p1a <-
      daily_summary %>%
      select(date, R_1x, R_3x, R_vs, R_vss) %>%
      na.omit() %>%
      pivot_longer(-date) %>%
      ggplot(aes(date, value, color = name)) +
      geom_line() +
      labs(title = "Performance")
    
    # 1B Daily market returns over time
    p1b <-
      daily_summary %>%
      select(date, r, R_1x) %>%
      na.omit() %>%
      ggplot(aes(date, 100*r)) +
      geom_line() +
      labs(title = "Daily Returns")
    
    # 1C Predicted annualized return over time
    p1c <-
      daily_summary %>%
      select(date, E_kn3, E_kn1, E_k0, E_kp1, E_kp3) %>%
      na.omit() %>%
      pivot_longer(-date) %>%
      mutate(value = if_else(value > 200, 200, value)) %>%
      ggplot(aes(date, value, color = name)) +
      geom_line() +
      labs(title = "Predicted Annualized Return")
    
    # 1d Optimal leverage across time
    p1d <-
      daily_summary %>%
      select(date, k_vss) %>%
      na.omit() %>%
      ggplot(aes(date, k_vss)) +
      geom_line() +
      labs(title = "Optimal Leverage")
    
    # 2 Moments
    p2 <-
      daily_summary %>%
      select(date, mean, sd, skew, kurtosis) %>%
      na.omit() %>%
      pivot_longer(-date) %>%
      ggplot(aes(date, value)) +
      geom_line() +
      facet_grid(name ~ ., scales = "free") +
      labs(title = "Daily Moments")
    
    # 3 Model parameters across times
    p3 <-
      par_summary %>%
      select(date, variable, median, q5, q95) %>%
      ggplot(aes(date, median)) +
      geom_line() +
      geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.5) +
      facet_grid(variable ~ ., scales = "free") +
      labs(title = "Posterior Marginal Distributions")
    
    # Plot all
    p_all <- arrangeGrob(arrangeGrob(p1a, p1b, p1c, p1d, ncol = 1), p2, p3, nrow = 1)
    ggsave(file=paste0(kOutdir, "dashboard.pdf"), plot = p_all, width = 36, height = 12)
    
    #############
    # Save
    #############
    
    # Save data
    saveRDS(y_oos, paste0(kOutdir, "y_oos.RDS"))
    saveRDS(daily_summary, paste0(kOutdir, "daily_summary.RDS"))
    saveRDS(par_summary, paste0(kOutdir, "par_summary.RDS"))
    
  } else {
    # We're fine to just recycle old y_rep samples
    w <- weights(psis_obj, log = FALSE, normalize = TRUE)[, 1]
    y_oos[,i+1] <- sample(y_oos0[,i-i_refit+1], replace = TRUE, prob = w)
    
    #############
    # Update
    #############
    
    daily_summary$mean[i+1] <- mean(y_oos[,i+1])
    daily_summary$sd[i+1] <- sd(y_oos[,i+1])
    daily_summary$skew[i+1] <- mean((y_oos[,i+1]-mean(y_oos[,i+1]))^3) / sd(y_oos[,i+1])^3
    daily_summary$kurtosis[i+1] <- mean((y_oos[,i+1]-mean(y_oos[,i+1]))^4) / sd(y_oos[,i+1])^4
    
    daily_summary$E_kn3[i+1] <- get_annualized_expected_leveraged_return(y_oos[,i+1], -3)
    daily_summary$E_kn1[i+1] <- get_annualized_expected_leveraged_return(y_oos[,i+1], -1)
    daily_summary$E_k0[i+1] <- get_annualized_expected_leveraged_return(y_oos[,i+1], 0)
    daily_summary$E_kp1[i+1] <- get_annualized_expected_leveraged_return(y_oos[,i+1], 1)
    daily_summary$E_kp3[i+1] <- get_annualized_expected_leveraged_return(y_oos[,i+1], 3)
    
    
    daily_summary$k_vs[i+1] <-
      get_optimal_leverage(c(daily_summary$E_k0[i+1],
                             daily_summary$E_kp1[i+1],
                             daily_summary$E_kp3[i+1]),
                           c(0, 1, 3))
    
    daily_summary$k_vss[i+1] <-
      get_optimal_leverage(c(daily_summary$E_kn3[i+1],
                             daily_summary$E_kn1[i+1],
                             daily_summary$E_k0[i+1],
                             daily_summary$E_kp1[i+1],
                             daily_summary$E_kp3[i+1]),
                           c(-3, -1, 0, 1, 3))
    
    r <- daily_summary$r[i+1]
    daily_summary$r_1x[i+1] <- 1*r
    daily_summary$r_3x[i+1] <- 3*r
    daily_summary$r_vs[i+1] <- daily_summary$k_vs[i+1]*r
    daily_summary$r_vss[i+1] <- daily_summary$k_vss[i+1]*r
    
    daily_summary$R_1x[i+1] <- prod(1+daily_summary$r_1x[(L+1):(i+1)])
    daily_summary$R_3x[i+1] <- prod(1+daily_summary$r_3x[(L+1):(i+1)])
    daily_summary$R_vs[i+1] <- prod(1+daily_summary$r_vs[(L+1):(i+1)])
    daily_summary$R_vss[i+1] <- prod(1+daily_summary$r_vss[(L+1):(i+1)])
  }
} 


















