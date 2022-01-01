library(tidyverse)
library(lubridate)
library(tidyquant)
library(cmdstanr)
library(posterior)
library(tibbletime)
library(plotly)
library(bayesplot)
library(gridExtra)

##############################
# Define constants
##############################

kStanModelFile <- "stan/joint_spx_tlt.stan"

kOutdir <- "data/joint_stock_bonds_studentt5_06_27_21/"

kTrainStartDate <- "2019-01-01"
kTestEndDate <- "2021-06-26"
kNumTrainDays <- 252

kNumStanSamples <- 1e3

kStanParNames <- c("lp__", "mus_ret", "mub_ret", "mus", "mub", "phis", "phib",
                   "sigmas", "sigmab",
                   "Rho[1,2]", "Rho[1,3]", "Rho[1,4]", "Rho[2,3]", "Rho[2,4]", "Rho[3,4]")
kNumStanPars <- length(kStanParNames)

# kLeverageLevels <-
#   expand.grid(p = c(0.25, 0.5, 0.75),
#               ks = c(0, 1, 3),
#               kb = c(0, 1, 3)) %>%
#   filter(!(ks == 0 & kb == 0)) %>%
#   bind_rows(expand.grid(p = 0, ks = 0, kb = c(1, 3))) %>%
#   bind_rows(expand.grid(p = 1, ks = c(0, 1, 3), kb = 0)) %>%
#   as_tibble()

kLeverageLevels <-
  tribble(~stocks, ~bonds,
          0, 0,
          0, 1,
          0, 2,
          0, 3,
          1, 0,
          1, 1,
          1, 2,
          2, 0,
          2, 1,
          3, 0) %>%
  mutate(cash = 3 - stocks - bonds)

##############################
# Define utility functions
##############################

# compute_expected_ret <- function(p, ks, kb) {
#   next_day_ret_draws %>%
#     mutate(total_ret = 1 + p*ks*ys_rep/100 + (1-p)*kb*yb_rep/100) %>%
#     mutate(log_total_ret = log(total_ret)) %>%
#     pull(log_total_ret) %>%
#     mean() %>%
#     exp() %>%
#     (function(x) 100*(x^252 -1))
# }

compute_expected_ret <- function(stocks, bonds, cash) {
  next_day_ret_draws %>%
    mutate(total_ret = 1 + stocks*ys_rep/100 + bonds*yb_rep/100) %>%
    mutate(log_total_ret = log(total_ret)) %>%
    pull(log_total_ret) %>%
    mean() %>%
    exp() %>%
    (function(x) 100*(x^252 -1))
}

##############################
# 1) Get data
##############################

spy <-
  tq_get(c("^GSPC", "TLT"),
         get  = "stock.prices",
         from = kTrainStartDate,
         to   = kTestEndDate) %>%
  group_by(symbol) %>%
  arrange(symbol, date) %>%
  mutate(pct_return = 100 * (adjusted - lag(adjusted)) / lag(adjusted)) %>%
  mutate(pct_return = ifelse(pct_return == 0.0, 1e-4, pct_return)) %>%
  ungroup() %>%
  select(date, symbol, pct_return) %>%
  pivot_wider(names_from = symbol, values_from = pct_return) %>%
  na.omit() %>%
  mutate(idx = row_number())

##############################
# 2) Initialize Data Objects
##############################

n_test_days <- nrow(spy) - kNumTrainDays

# Initialize data to fill up during sampling
y_rep_holdout <- array(0.0, c(2, kNumStanSamples * 4, n_test_days))

vol <-
  spy[(kNumTrainDays + 1):nrow(spy),] %>%
  select(idx, date) %>%
  mutate(m_stocks = 0, m_bonds = 0,
         s_stocks = 1, s_bonds = 1,
         k_stocks = 0, k_bonds = 0, k_cash = 0, k_vs = 0)

test <-
  spy[(kNumTrainDays + 1):nrow(spy),] %>%
  mutate(upro = 3 * `^GSPC`,
         tmf = 3* TLT,
         half = 1.5*`^GSPC` + 1.5*TLT,
         vs = 0,
         vsb = 0,
         cum_spx = cumprod(1 + `^GSPC` / 100),
         cum_tlt = cumprod(1 + TLT / 100),
         cum_upro = cumprod(1 + upro / 100),
         cum_tmf = cumprod(1 + tmf / 100),
         cum_half = cumprod(1 + half / 100),
         cum_vs = 1,
         cum_vsb = 1)

diagnostics <-
  spy[(kNumTrainDays + 1):nrow(spy),] %>%
  select(date) %>%
  mutate(p_stocks = 0.5, z_stocks = 0, cumx2_stocks = pchisq(cumsum(z_stocks^2), df = row_number())) %>%
  mutate(p_bonds = 0.5, z_bonds = 0, cumx2_bonds = pchisq(cumsum(z_bonds^2), df = row_number()))
  
parameters <-
  tibble(date = rep(spy$date[(kNumTrainDays + 1):nrow(spy)], each = kNumStanPars)) %>%
  mutate(variable = "", q5 = 1, median = 2, q95 = 3)

##############################
# 3) Main loop
##############################

model <- cmdstan_model(kStanModelFile)

# y_rep_holdout <- readRDS(paste0(kOutdir, "y_rep_holdout.RDS"))
# vol <- readRDS(paste0(kOutdir, "vol.RDS"))
# test <- readRDS(paste0(kOutdir, "test.RDS"))
# diagnostics <- readRDS(paste0(kOutdir, "diagnostics.RDS"))
# parameters <- readRDS(paste0(kOutdir, "parameters.RDS"))

for(i in 7:n_test_days) {
  
  # Set up train and test data
  train <- spy[i:(i + kNumTrainDays - 1),]
  
  # Set up Stan data
  dat <- list(T = nrow(train), ys = train$`^GSPC`, yb = train$TLT)
  
  # Fit model
  fit <-
    model$sample(
      data = dat,
      seed = 11112,
      iter_warmup = 1e3,
      iter_sampling = kNumStanSamples,
      chains = 4,
      parallel_chains = 4,
      refresh = 2e2,
      max_treedepth = 10,
      adapt_delta = 0.8,
      show_messages = FALSE
    )
  
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`")
  print(paste("Finished fitting model", i, "/", n_test_days))
  
  # Extract draws
  next_day_ret_draws <- fit$draws(c("ys_rep", "yb_rep")) %>% as_draws_df() %>% select(ys_rep, yb_rep)
  y_rep_holdout[,,i] <- t(as.matrix(next_day_ret_draws))
  
  # Get expected SD of posterior predictive draws
  vol$m_stocks[i] <- mean(next_day_ret_draws$ys_rep)
  vol$m_bonds[i] <- mean(next_day_ret_draws$yb_rep)
  vol$s_stocks[i] <- sd(next_day_ret_draws$ys_rep)
  vol$s_bonds[i] <- sd(next_day_ret_draws$yb_rep)
  
  # Set optimal k
  expec_ann_returns <- 
    kLeverageLevels %>%
    mutate(expec_ann_ret = pmap_dbl(., compute_expected_ret)) %>%
    mutate(expec_ann_ret = ifelse(is.nan(expec_ann_ret), NA_real_, expec_ann_ret)) %>%
    arrange(desc(expec_ann_ret))
  
  vol$k_stocks[i] <- expec_ann_returns$stocks[1]
  vol$k_bonds[i] <- expec_ann_returns$bonds[1]
  vol$k_cash[i] <- expec_ann_returns$cash[1]

  vol$k_vs[i] <- expec_ann_returns %>% filter(bonds == 0) %>% pull(stocks) %>% head(1)
  
  # Set test results
  test$vs[i] <- vol$k_vs[i] * test$`^GSPC`[i]
  test$vsb[i] <- vol$k_stocks[i]*test$`^GSPC`[i] + vol$k_bonds[i]*test$TLT[i]
  test <- mutate(test, cum_vs = cumprod(1 + vs / 100))
  test <- mutate(test, cum_vsb = cumprod(1 + vsb / 100))
  
  # Set diagnostics
  actual_return_stocks <- spy$`^GSPC`[i + kNumTrainDays]
  actual_return_bonds <- spy$TLT[i + kNumTrainDays]
  
  diagnostics$p_stocks[i] <- mean(y_rep_holdout[1,,i] < actual_return_stocks)
  diagnostics$z_stocks[i] <- qnorm(diagnostics$p_stocks[i])
  diagnostics <- mutate(diagnostics, cumx2_stocks = pchisq(cumsum(z_stocks^2), df = row_number()))
  
  diagnostics$p_bonds[i] <- mean(y_rep_holdout[2,,i] < actual_return_bonds)
  diagnostics$z_bonds[i] <- qnorm(diagnostics$p_bonds[i])
  diagnostics <- mutate(diagnostics, cumx2_bonds = pchisq(cumsum(z_bonds^2), df = row_number()))
  
  # Extract and set parameter summaries
  par_start_idx <- (i-1)*kNumStanPars + 1
  par_end_idx <- par_start_idx + kNumStanPars - 1
  
  parameters[par_start_idx:par_end_idx, 2:5] <-
    fit$summary(kStanParNames) %>%
    select(variable, q5, median, q95)
  
  if(i %% 1 == 0) {
    dev.off()
    
    # Make left column plots
    plot1a <-
      qplot(y_rep_holdout[1,,i], binwidth = 0.1) +
      geom_vline(xintercept = actual_return_stocks, color = "red") +
      xlab(paste("SPX Percent Return", spy$date[i + kNumTrainDays])) +
      labs(title = "SPX LFO Predicted Return Vs. Actual",
           subtitle = paste(paste("Pct. Ret.:", round(actual_return_stocks, 2)),
                            paste("P-val:", round(diagnostics$p_stocks[i], 4)),
                            paste("Z-score:", round(diagnostics$z_stocks[i], 2)),
                            paste("Opt. Lev:", vol$k_stocks[i]),
                            paste("Opt. Lev VS:", vol$k_vs[i]),
                            sep = " // "))
    
    plot1b <-
      qplot(y_rep_holdout[2,,i], binwidth = 0.1) +
      geom_vline(xintercept = actual_return_bonds, color = "red") +
      xlab(paste("TLT Percent Return", spy$date[i + kNumTrainDays])) +
      labs(title = "TLT LFO Predicted Return Vs. Actual",
           subtitle = paste(paste("Pct. Ret.:", round(actual_return_bonds, 2)),
                            paste("P-val:", round(diagnostics$p_bonds[i], 4)),
                            paste("Z-score:", round(diagnostics$z_bonds[i], 2)),
                            paste("Opt. Lev:", vol$k_bonds[i]),
                            sep = " // "))
    
    plot1c <-
      qplot(y_rep_holdout[1,,i], y_rep_holdout[2,,i]) +
      geom_vline(xintercept = actual_return_stocks, color = "red") +
      geom_hline(yintercept = actual_return_bonds, color = "red") +
      xlab("SPX") +
      ylab("TLT")
    
    # Make middle column plots
    plot2a <-
      vol %>%
      head(i) %>%
      select(date, k_stocks:k_vs) %>%
      pivot_longer(k_stocks:k_cash, names_to = c("var", "asset"), names_sep = "_") %>%
      mutate(asset = factor(asset, levels = c("stocks", "bonds", "cash"))) %>%
      ggplot() +
      geom_area(aes(date, value, fill = asset)) +
      geom_line(aes(date, k_vs)) +
      facet_grid(var ~ ., scales = "free") +
      labs(title = "Opt. Lev. (k) Over Time")
    
    plot2b <-
      vol %>%
      head(i) %>%
      select(date, m_stocks, m_bonds, s_stocks, s_bonds) %>%
      pivot_longer(-date, names_to = c("var", "asset"), names_sep = "_") %>%
      ggplot(aes(date, value, color = asset)) +
      geom_point() +
      geom_line() +
      facet_grid(var ~ ., scales = "free") +
      labs(title = "Mean & SD Over Time")
    
    plot2c <-
      diagnostics %>%
      head(i) %>%
      select(-p_stocks, -p_bonds) %>%
      pivot_longer(-date, names_to = c("var", "asset") , names_sep = "_") %>%
      ggplot(aes(date, value, color = asset)) +
      geom_point() +
      geom_line() +
      facet_grid(var + asset ~ ., scales = "free") +
      labs(title = "Z-Scores and Cum. X^2 Quantiles Over Time")
    
    # Make third column plots
    plot3a <-
      test %>%
      head(i) %>%
      select(date, upro, half, vs, vsb) %>%
      pivot_longer(-date) %>%
      ggplot(aes(value, fill = name)) +
      geom_histogram() +
      facet_grid(name ~ .) +
      labs(title = "Hist. of Daily Returns")
    
    plot3b <-
      test %>%
      head(i) %>%
      select(date, upro, half, vs, vsb, cum_upro, cum_half, cum_vs, cum_vsb) %>%
      pivot_longer(-date) %>%
      mutate(strat = str_extract(name, "(?<=_).*"),
             cum = str_detect(name, "cum")) %>%
      mutate(strat = ifelse(is.na(strat), name, strat)) %>%
      select(-name) %>%
      ggplot(aes(date, value, color = strat)) +
      geom_line() +
      facet_grid(cum ~ ., scales = "free") +
      labs(title = "Daily and Cum. Returns")
    
    plot3c <-
      test %>%
      head(i) %>%
      ggplot(aes(`^GSPC`, vsb)) +
      geom_abline(slope = 1, intercept = 0) +
      geom_abline(slope = 3, intercept = 0) +
      geom_point() +
      labs(title = "Daily Returns of VSB vs. SPY")
    
    # Make 4th column plots
    plot4a <-
      parameters %>%
      head(i*kNumStanPars) %>%
      ggplot(aes(date, median)) +
      geom_point() +
      geom_errorbar(aes(ymin = q5, ymax = q95)) +
      facet_grid(variable ~ ., scales = "free") +
      labs(title = "Model Parameter Posterior Intervals Over Time")
    
    # Plot all together
    grid.arrange(grid.arrange(plot1a, plot1b, plot1c, ncol = 1),
                 grid.arrange(plot2a, plot2b, plot2c, ncol = 1),
                 grid.arrange(plot3a, plot3b, plot3c, ncol = 1),
                 grid.arrange(plot4a, ncol = 1),
                 nrow = 1)
    
    # Save stanfit and intermediary data
    saveRDS(y_rep_holdout, paste0(kOutdir, "y_rep_holdout.RDS"))
    saveRDS(vol, paste0(kOutdir, "vol.RDS"))
    saveRDS(test, paste0(kOutdir, "test.RDS"))
    saveRDS(diagnostics, paste0(kOutdir, "diagnostics.RDS"))
    saveRDS(parameters, paste0(kOutdir, "parameters.RDS"))
  }
  
  
}


