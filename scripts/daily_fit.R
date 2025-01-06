#!/usr/bin/env Rscript

library(dplyr)
library(purrr)
library(tidyquant)
library(cmdstanr)
library(posterior)
# library(twilio)
# library(gmailr)
library(discordr)

# Define Constants
kStanModelFile <- "homebuying/stan/sv.stan"
kNumTrainMonths <- 6
kNumTrainDays <- 126

kNumCalendarTrainDays <- round((kNumTrainDays / 252) * 365) + 100

kTodaysDate <- Sys.Date()

kNumStanSamples <- 5e3

kLeverageLevels <- c(0.0, 0.5, 1.0, 2.0, 3.0)

# Define helper function
compute_leverage <- function(r, k = c(0.0, 0.5, 1.0, 2.0, 3.0)) {
  expec_log_return <- map_dbl(k, function(k) mean(log(1 + k*r/100)))
  #k[which(expec_log_return == max(expec_log_return))]
  100 * (exp(expec_log_return)^252 - 1)
}

main <- function() {
  
  ##############################
  # 1) Get data
  ##############################
  
  # Get previous days
  spy <-
    tq_get(c("SPY"),
           get = "tiingo",
           from = kTodaysDate - kNumCalendarTrainDays,
           to = kTodaysDate) %>%
    arrange(date) %>%
    mutate(pct_return = 100 * (close - lag(close)) / lag(close)) %>%
    # if any returns are exactly zero this causes the lp to evaluate to zero
    mutate(pct_return = ifelse(pct_return == 0.0, 1e-4, pct_return)) %>%
    na.omit() %>%
    tail(kNumTrainDays)
  
  # Append on current price
  spy_today <-
    getQuote("SPY") %>%
    rename_with(~ tolower(gsub(".", "_", .x, fixed = TRUE))) %>%
    mutate(date = kTodaysDate) %>%
    rename(close = last) %>%
    select(-change, -`trade time`, -`% change`)
  
  spy <-
    bind_rows(spy, spy_today) %>%
    mutate(pct_return = 100 * (close - lag(close)) / lag(close)) %>%
    # if any returns are exactly zero this causes the lp to evaluate to zero
    mutate(pct_return = ifelse(pct_return == 0.0, 1e-4, pct_return)) %>%
    filter(!is.na(pct_return))
  
  ##############################
  # 2) Fit Model
  ##############################
  
  model <- cmdstan_model(kStanModelFile)
  
  dat <- list(T = nrow(spy), y = spy$pct_return)
  
  # Fit model
  fit <-
    model$sample(
      data = dat,
      seed = 1994,
      iter_warmup = 1e3,
      iter_sampling = kNumStanSamples,
      chains = 4,
      parallel_chains = 4,
      refresh = 1e4,
      max_treedepth = 10,
      adapt_delta = 0.99,
      show_messages = FALSE
    )
  
  ##############################
  # 3) Display Results
  ##############################
  
  y_rep_ahead <-
    fit$draws(c("y_rep_ahead")) %>%
    as_draws_matrix() %>%
    as.vector()
  
  # Set optimal k
  expec_ann_returns <- compute_leverage(y_rep_ahead, kLeverageLevels)
  expec_ann_returns[is.nan(expec_ann_returns)] <- NA_real_
  max_expec_ann_returns <- max(expec_ann_returns, na.rm = TRUE)
  opt_k <- kLeverageLevels[which(expec_ann_returns == max_expec_ann_returns)[1]]
  
  ##############################
  # 4) Message Update
  ##############################
  
  body <-
    paste("Return today:", round(tail(spy$pct_return, 1), 2),
          "\n",
          "Optimal leverage for tomorrow:", opt_k,
          "\n",
          "Expected annualized return:", round(max_expec_ann_returns, 2))

  # Send via discord
  conn_obj <-
    create_discord_connection(webhook = 'https://discord.com/api/webhooks/1134185904282550283/qQjKc9PC8R1fQFxOVShhwK9kqAlUW-v3d1h3yBlQt54BteUcLyX7w3U4qfxZHys5zORx',
                              username = 'Gravekimjani Daily',
                              set_default = TRUE)
  
  send_webhook_message(body)
  
  # Send via Twilio
  # tw_send_message(
  #   to = Sys.getenv("A_PN"),
  #   from = "+18722405834",
  #   body = body
  # )
  # 
  # tw_send_message(
  #   to = Sys.getenv("J_PN"),
  #   from = "+18722405834",
  #   body = body
  # )
  # 
  # tw_send_message(
  #   to = Sys.getenv("C_PN"),
  #   from = "+18722405834",
  #   body = body
  # )
  # 
  # tw_send_message(
  #   to = Sys.getenv("T_PN"),
  #   from = "+18722405834",
  #   body = body
  # )
  
  # Code to send email
  # gm_auth_configure(path = Sys.getenv("GMAILR_APP"))
  # gm_auth(email = TRUE, cache = Sys.getenv("GMAILR_CACHE"))
  # 
  # addresses <-
  #   c(Sys.getenv("A_E"),
  #     Sys.getenv("J_E"),
  #     Sys.getenv("N_E"),
  #     Sys.getenv("C_E"),
  #     Sys.getenv("T_E"))
  # 
  # safe_send <- function(address) {
  #   email <-
  #     gm_mime() %>%
  #     gm_to(address) %>%
  #     gm_from("gravkimjani@gmail.com") %>%
  #     gm_subject(paste("Optimal leverage for tomorrow:", opt_k)) %>%
  #     gm_text_body(body)
  #   
  #   safely(gm_send_message)(email)
  # }
  # 
  # map(addresses, safe_send)
}

if (sys.nframe() == 0) {
  main()
}