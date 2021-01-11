data {
  int<lower=0> T;
  vector[T] y;
  real<lower=0> sigma1;
  
  real<lower=0> s;
}
parameters {
  // GARCH parameters
  real mu;
  real alpha0;
  real alpha1;
  real alpha2;
  real beta1;
  
  // SV parameters
  real rho;
  real<lower=0> c;
  vector[T] v_raw;
}
transformed parameters {

  vector[T] v = c * v_raw;

  vector[T] h_expected;
  vector[T] sigma_expected;
  vector[T] h;      // log of volatility squared
  vector[T] sigma;  // volatility (SD of returns)
  vector[T] e;
  
  // Set for time 1
  h_expected[1] = log(sigma1^2);
  sigma_expected[1] = exp(0.5 * h_expected[1]);
  h[1] = h_expected[1];
  sigma[1] = exp(0.5 * h[1]);
  e[1] = y[1] / sigma_expected[1];
  
  for (t in 2:T) {

    h_expected[t] =
      alpha0 +
      alpha1 * e[t-1] +
      alpha2 * e[t-1]^2 / (1 + e[t-1]^2) +
      beta1 * h_expected[t-1];
      
    sigma_expected[t] = exp(0.5 * h_expected[t]);
    
    h[t] = h_expected[t] + v[t];
    sigma[t] = exp(0.5 * h[t]);
    
    e[t] = y[t] / sigma_expected[t];
  }
}
model {
  // Prior on parameters
  mu ~ normal(0, 0.2);
  rho ~ normal(0, 1);
  
  // SV
  c ~ exponential(1 / s);
  v_raw ~ std_normal();

  // Likelihood
  y ~ normal(mu + rho * v, sigma);
}
generated quantities {
  // Turn off for production
  // vector[T] v_rep = to_vector(normal_rng(rep_array(0, T), c));
  // vector[T] y_rep = to_vector(normal_rng(mu + rho * v_rep, exp(0.5 * (h_expected + v_rep))));

  real h_expected_ahead = 
    alpha0 +
      alpha1 * e[T] +
      alpha2 * e[T]^2 / (1 + e[T]^2) +
      beta1 * h_expected[T];
  real v_rep_ahead = normal_rng(0.0, c);
  real h_rep_ahead = h_expected_ahead + v_rep_ahead;
  real sigma_rep_ahead = exp(0.5 * h_rep_ahead);
  real y_rep_ahead = normal_rng(mu + rho * v_rep_ahead, sigma_rep_ahead);
}
