data {
  int<lower=0> T;
  vector[T] r_t;                // Returns
}
parameters {
  real mu;                      // Mean of returns
  real<lower=0> log_sigma1;     // SD for first observation
  vector[T-1] eps;              // Horseshoe shock to SD
  
  // Horseshoe prior parameters
  real<lower=0> tau;
  vector<lower=0>[T-1] lambda;
}
transformed parameters {
  // Set time varying volatilities (mostly constant with discontinuous jumps) 
  vector[T] log_sigma_t;
  
  log_sigma_t[1] = log_sigma1;
  for (t in 2:T)
    log_sigma_t[t] += log_sigma_t[t-1] + eps[t];
}
model {
  // Horseshoe prior on variance shocks (usually variance will remain unchanged week-to-week)
  tau ~ cauchy(0, 1);
  lambda ~ cauchy(0, 1);
  eps ~ normal(0, tau*lambda);
  
  // Likelihood of returns data
  r_t ~ normal(mu, exp(log_sigma_t));
}

