data {
  int<lower=0> T;   // # time points (equally spaced)
  vector[T] y;      // mean corrected return at time t
}
parameters {
  real mu_ret;
  real rho;

  // SV
  real mu;                     // mean log volatility
  real<lower=-1, upper=1> phi;  // persistence of volatility
  real<lower=0> sigma;         // white noise shock scale
  
  vector[T] v_raw;             // std log volatility time t
}
transformed parameters {
  vector[T] v = sigma * v_raw; 
  vector[T] h;
  vector[T] s;
  
  h[1] = v[1] / sqrt(1 - phi^2);
  
  for (t in 2:T) {
    h[t] = phi * h[t-1] + v[t];
  }
  
  s = exp((mu + h) / 2);
}
model {
  // prior
  mu_ret ~ normal(0, 0.2);
  rho ~ normal(-1, 1);

  (phi+1)/2.0 ~ beta(8, 2);
  sigma ~ exponential(10);
  mu ~ normal(0, 1);
  
  // likelihood
  v_raw ~ std_normal();
  y ~ normal(mu_ret + rho * v .* s, s);
  
}
generated quantities {
  real v_rep_ahead = sigma * normal_rng(0, 1);
  real h_rep_ahead = phi * h[T] + v_rep_ahead;
  real s_rep_ahead = exp((mu + h_rep_ahead) / 2);
  real y_rep_ahead = normal_rng(mu_ret + rho * v_rep_ahead  * s_rep_ahead, s_rep_ahead);
}
