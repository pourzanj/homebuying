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

kStanModelFile <- "stan/sv.stan"

kOutdir <- "data/sv_with_leverage_01_13_20/"

kTrainStartDate <- "1990-06-01"
kTestEndDate <- "2021-01-13"
kNumTrainDays <- 126

kNumStanSamples <- 1e3

kLeverageLevels <- c(0.0, 0.5, 1.0, 2.0, 3.0)

kStanParNames <- c("mu_ret", "rho", "mu", "phi", "sigma", "lp__")
kNumStanPars <- length(kStanParNames)

##############################
# Define utility functions
##############################

compute_leverage <- function(r, k = c(0.0, 0.5, 1.0, 2.0, 3.0)) {
  expec_log_return <- map_dbl(k, function(k) mean(log(1 + k*r/100)))
  #k[which(expec_log_return == max(expec_log_return))]
  100 * (exp(expec_log_return)^252 - 1)
}

##############################
# 1) Get data
##############################

spy <-
  tq_get(c("SPY"),
         get = "stock.prices",
         from = kTrainStartDate,
         to = kTestEndDate) %>%
  arrange(date) %>%
  mutate(pct_return = 100 * (close - lag(close)) / lag(close)) %>%
  # if any returns are exactly zero this causes the lp to evaluate to zero
  mutate(pct_return = ifelse(pct_return == 0.0, 1e-4, pct_return)) %>%
  na.omit() %>%
  mutate(idx = row_number())

##############################
# 2) Initialize Data Objects
##############################

n_test_days <- nrow(spy) - kNumTrainDays

# Initialize data to fill up during sampling
# y_rep_holdout <- matrix(0.0, nrow = kNumStanSamples * 4, ncol = n_test_days)
# 
# vol <-
#   spy[(kNumTrainDays + 1):nrow(spy),] %>%
#   select(idx, date) %>%
#   mutate(s = 1, k = 0)
# 
# test <-
#   spy[(kNumTrainDays + 1):nrow(spy),] %>%
#   select(idx, date, ret_spy = pct_return) %>%
#   mutate(ret_spy3x = 3 * ret_spy, ret_vs = 0,
#          cum_spy = cumprod(1 + ret_spy / 100),
#          cum_spy3x = cumprod(1 + ret_spy3x / 100),
#          cum_vs = 1)
# 
# diagnostics <-
#   spy[(kNumTrainDays + 1):nrow(spy),] %>%
#   select(date) %>%
#   mutate(p = 0.5, z = 0, cum_x2 = pchisq(cumsum(z^2), df = row_number()))
# 
# parameters <-
#   tibble(date = rep(spy$date[(kNumTrainDays + 1):nrow(spy)], each = kNumStanPars)) %>%
#   mutate(variable = "", q5 = 1, median = 2, q95 = 3)

##############################
# 3) Main loop
##############################

model <- cmdstan_model(kStanModelFile)

y_rep_holdout <- readRDS(paste0(kOutdir, "y_rep_holdout.RDS"))
vol <- readRDS(paste0(kOutdir, "vol.RDS"))
test <- readRDS(paste0(kOutdir, "test.RDS"))
diagnostics <- readRDS(paste0(kOutdir, "diagnostics.RDS"))
parameters <- readRDS(paste0(kOutdir, "parameters.RDS"))

for(i in 4821:n_test_days) {
  
  # Set up train and test data
  train <- spy[i:(i + kNumTrainDays - 1),]

  # Set up Stan data
  dat <-
    list(T = nrow(train),
         y = train$pct_return,
         sigma1 = sd(train$pct_return[1:5]),
         s = 0.7)
  
  # Fit model
  fit <-
    model$sample(
      data = dat,
      seed = 11112,
      iter_warmup = 1e3,
      iter_sampling = kNumStanSamples,
      chains = 4,
      parallel_chains = 4,
      refresh = 5e2,
      max_treedepth = 10,
      adapt_delta = 0.99,
      show_messages = FALSE
    )
  
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`")
  print(paste("Finished fitting model", i, "/", n_test_days))
  
  # Extract draws
  y_rep_holdout[, i] <-
    fit$draws(c("y_rep_ahead")) %>%
    as_draws_matrix() %>%
    as.vector()
  
  # Get expected SD of posterior predictive draws
  vol$s[i] <- sd(y_rep_holdout[, i])
  
  # Set optimal k
  expec_ann_returns <- compute_leverage(y_rep_holdout[, i], kLeverageLevels)
  expec_ann_returns[is.nan(expec_ann_returns)] <- NA_real_
  vol$k[i] <- kLeverageLevels[which(expec_ann_returns == max(expec_ann_returns, na.rm = TRUE))[1]]
  
  # Set test results
  test$ret_vs[i] <- vol$k[i] * test$ret_spy[i]
  test <- mutate(test, cum_vs = cumprod(1 + ret_vs / 100))
  
  # Set diagnostics
  actual_return <- spy$pct_return[i + kNumTrainDays]
  
  diagnostics$p[i] <- mean(y_rep_holdout[, i] < actual_return)
  diagnostics$z[i] <- qnorm(diagnostics$p[i])
  diagnostics <- mutate(diagnostics, cum_x2 = pchisq(cumsum(z^2), df = row_number()))
  
  # Extract and set parameter summaries
  par_start_idx <- (i-1)*kNumStanPars + 1
  par_end_idx <- par_start_idx + kNumStanPars - 1
  
  parameters[par_start_idx:par_end_idx, 2:5] <-
    fit$summary(kStanParNames) %>%
    select(variable, q5, median, q95)
  
  if(i %% 20 == 1) {
    dev.off()
    
    # Make left column plots
    plot1a <-
      qplot(y_rep_holdout[, i], binwidth = 0.1) +
      geom_vline(xintercept = actual_return, color = "red") +
      xlab(paste("SPX Percent Return", spy$date[i + kNumTrainDays])) +
      labs(title = "SPX LFO Predicted Return Vs. Actual",
           subtitle = paste(paste("Day", i, "/", n_test_days),
                            paste("Pct. Ret.:", round(actual_return, 2)),
                            paste("P-val:", round(diagnostics$p[i], 4)),
                            paste("Z-score:", round(diagnostics$z[i], 2)),
                            paste("Opt. Lev:", vol$k[i]),
                            sep = " // "))
    
    plot1b <-
      tibble(k = kLeverageLevels, expec_ann_returns) %>%
      ggplot(aes(k, expec_ann_returns)) +
      geom_point() +
      geom_line() +
      xlab("Leverave Level (k)") +
      ylab("Expected Annual Return") +
      labs(title = "Returs vs. Leverage Level")
    
    # Make middle column plots
    plot2a <-
      vol %>%
      head(i) %>%
      select(-idx) %>%
      pivot_longer(-date) %>%
      mutate(name = factor(name, levels = c("s", "k"))) %>%
      ggplot(aes(date, value)) +
      geom_point() +
      geom_line() +
      facet_grid(name ~ ., scales = "free") +
      labs(title = "SD and Opt. Lev. (k) Over Time")
    
    plot2b <-
      diagnostics %>%
      head(i) %>%
      select(-p) %>%
      pivot_longer(-date) %>%
      ggplot(aes(date, value)) +
      geom_point() +
      geom_line() +
      facet_grid(name ~ ., scales = "free") +
      labs(title = "Z-Scores and Cum. X^2 Quantiles Over Time")
    
    # Make third column plots
    plot3a <-
      test %>%
      head(i) %>%
      select(date, ret_spy, ret_spy3x, ret_vs) %>%
      pivot_longer(-date) %>%
      ggplot(aes(value, fill = name)) +
      geom_histogram() +
      facet_grid(name ~ .) +
      labs(title = "Hist. of Daily Returns")
    
    plot3b <-
      test %>%
      head(i) %>%
      select(-idx) %>%
      pivot_longer(-date) %>%
      mutate(strat = str_extract(name, "(?<=_).*"),
             var = str_extract(name, "^.*(?=(_))")) %>%
      filter(!(strat == "spy3x" && var == "ret")) %>%
      select(-name) %>%
      ggplot(aes(date, value, color = strat)) +
      geom_line() +
      facet_grid(var ~ ., scales = "free") +
      labs(title = "Daily and Cum. Returns")
    
    plot3c <-
      test %>%
      head(i) %>%
      ggplot(aes(ret_spy, ret_vs)) +
      geom_abline(slope = 1, intercept = 0) +
      geom_abline(slope = 3, intercept = 0) +
      geom_point() +
      labs(title = "Daily Returns of VS vs. SPY")
    
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
    grid.arrange(grid.arrange(plot1a, plot1b, ncol = 1),
                 grid.arrange(plot2a, plot2b, ncol = 1),
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


