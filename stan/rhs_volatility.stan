data {
  int<lower=0> T;
  vector[T] r_t;                // Returns
  
  // Prior hyperparameters
  real<lower=0> tau0;           // Sparsity of shocks
  real<lower=0> c;              // Scale of slab of shocks
}
parameters {
  real mu;                      // Mean of returns
  vector[T] eps;                // Horseshoe shock to SD
  
  // Horseshoe prior parameters
  real<lower=0> tau;
  vector<lower=0>[T] lambda;
}
transformed parameters {
  // Set time varying volatilities
  vector[T] log_sigma;
  
  log_sigma[1] = eps[1];
  for (t in 2:T) {
    // Volatility this week is equal to volatility last week plus a horseshoe
    // shock that is usually zero, but occasionally some large non-zero change.
    // Thus for weeks at a time volatility effectively does not change.
    log_sigma[t] = log_sigma[t-1] + eps[t];
  }
  // print("~~~~~");
  // print("eps: ", eps[1:5]);
  // print("log_sigma: ", log_sigma[1:5]);
}
model {
  // Horseshoe prior on variance shocks (usually variance will remain unchanged
  // week-to-week). First eps statement controls spike, second statement is slab
  tau ~ cauchy(0, tau0);
  lambda ~ cauchy(0, 1);
  eps ~ normal(0, tau*lambda);
  eps ~ normal(0, c);
  
  // Likelihood of returns data
  r_t ~ normal(mu, exp(log_sigma));
}
generated quantities {
  real r_t_rep[T] = normal_rng(mu, exp(log_sigma));
}

