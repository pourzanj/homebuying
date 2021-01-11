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

kTrainStartDate <- "2017-01-01"
kTestEndDate <- "2021-01-10"
kNumTrainDays <- 502

kNumStanSamples <- 4e3

kLeverageLevels <- c(0.0, 0.5, 1.0, 2.0, 3.0)

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

# Matrix of posterior predictive draws for each day
y_rep_holdout <- matrix(0.0, nrow = kNumStanSamples * 4, ncol = n_test_days)

# Results on test set
vol <-
  spy[(kNumTrainDays + 1):nrow(spy),] %>%
  select(idx, date) %>%
  mutate(s = 1, k = 0)

test <-
  spy[(kNumTrainDays + 1):nrow(spy),] %>%
  select(idx, date, ret_spy = pct_return) %>%
  mutate(ret_spy3x = 3 * ret_spy, ret_vs = 0,
         cum_spy = cumprod(1 + ret_spy / 100),
         cum_spy3x = cumprod(1 + ret_spy3x / 100),
         cum_vs = 1)

diagnostics <-
  spy[(kNumTrainDays + 1):nrow(spy),] %>%
  select(date) %>%
  mutate(p = 0.5, z = 0, cum_x2 = pchisq(cumsum(z^2), df = row_number()))

parameters <-
  tibble(date = rep(spy$date[(kNumTrainDays + 1):nrow(spy)], each = 9)) %>%
  mutate(variable = "", q5 = 1, median = 2, q95 = 3)

##############################
# 3) Main loop
##############################

model <- cmdstan_model("stan/egarch_random_effects.stan")

for(i in 1:n_test_days) {
  
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
      seed = 102115,
      iter_warmup = 1e3,
      iter_sampling = kNumStanSamples,
      chains = 4,
      parallel_chains = 4,
      refresh = 1e3,
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
  par_start_idx <- (i-1)*9 + 1
  par_end_idx <- par_start_idx + 9 - 1
  
  parameters[par_start_idx:par_end_idx, 2:5] <-
    fit$summary(c("mu", "alpha0", "alpha1", "alpha2", "beta1", "c", "rho", "lp__", "h_rep_ahead")) %>%
    select(variable, q5, median, q95)
  
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
    labs(title = "Returns vs. Leverage Level")
  
  plot2b <-
    diagnostics %>%
    head(i) %>%
    select(-p) %>%
    pivot_longer(-date) %>%
    ggplot(aes(date, value)) +
    geom_point() +
    geom_line() +
    facet_grid(name ~ ., scales = "free") +
    labs(title = "Predicted SD and Optimal Leverage (k) Over Time")
  
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
    head(i*8) %>%
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
  # fit$save_object(file = paste0("data/run_01_10_20/fits/fit_", i, ".rds"))
  saveRDS(y_rep_holdout, "data/run_01_10_20/y_rep_holdout.RDS")
  saveRDS(vol, "data/run_01_10_20/vol.RDS")
  saveRDS(test, "data/run_01_10_20/test.RDS")
  saveRDS(diagnostics, "data/run_01_10_20/diagnostics.RDS")
  saveRDS(parameters, "data/run_01_10_20/parameters.RDS")
}


