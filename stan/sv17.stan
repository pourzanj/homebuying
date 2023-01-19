data {
  int<lower=0> T;   // # time points (equally spaced)
  vector[T] y;      // mean corrected return at time t
  
  int<lower=2> M;
  vector[M] y_ahead;
}

parameters {
  real mu_ret;

  real<lower=0> gamma;

  // SV
  real mu;                     // mean log volatility
  real<lower=-1, upper=1> phi;  // persistence of volatility
  real<lower=0> sigma;         // white noise shock scale
  
  vector[T] v;             // std log volatility time t
  vector[M] v_rep_ahead;
  
}
transformed parameters {
  vector[T] h;
  vector[T] s;

  h[1] = sigma*v[1] / sqrt(1 - phi^2);
  
  for (t in 2:T) {
    h[t] = phi*h[t-1] + sigma*v[t];
  }
  s = exp((mu + h) / gamma);
}
model {
  // prior
  // mu_ret ~ normal(0, 0.2);
  // rho ~ normal(-1, 1);

  // (phi+1)/2.0 ~ beta(8, 2);
  // sigma ~ exponential(10);
  // mu ~ normal(0, 1);
  
  gamma ~ normal(2, 1);
  
  // likelihood
  v ~ std_normal();
  y ~ normal(mu_ret, s);
  
  // for LOO
  v_rep_ahead ~ std_normal();
  
}
generated quantities {
  vector[T] v_rep = to_vector(normal_rng(rep_array(0.0, T), 1));
  vector[T] h_rep;
  vector[T] s_rep;
  vector[T] y_rep;

  vector[T] log_lik;
  real sum_log_lik;
  vector[T] log_lik_ahead;
  real sum_log_lik_ahead;
  
  vector[M] h_rep_ahead;
  vector[M] s_rep_ahead;
  vector[M] loglik_rep_ahead;

  h_rep[1] = sigma*v_rep[1] / sqrt(1 - phi^2);
  
  for(t in 2:T) {
    h_rep[t] = phi*h[t-1] + sigma*v_rep[t];
    s_rep[t] = exp((mu+h_rep[t])/gamma);
    y_rep[t] = normal_rng(mu_ret, s_rep[t]);
    log_lik[t] = normal_lpdf(y[t] | mu_ret, s[t]);
    log_lik_ahead[t] = normal_lpdf(y[t] | mu_ret, s_rep[t]);
  }
  sum_log_lik = sum(log_lik[2:T]);
  sum_log_lik_ahead = sum(log_lik_ahead[2:T]);
  
  // Compute rep_ahead and loglik_rep_ahead
  h_rep_ahead[1] = phi*h[T] + sigma*v_rep_ahead[1];
  s_rep_ahead[1] = exp((mu + h_rep_ahead[1]) / gamma);
  loglik_rep_ahead[1] = normal_lpdf(y_ahead[1] | mu_ret, s_rep_ahead[1]);
  for (t in 2:M) {
    h_rep_ahead[t] = phi*h_rep_ahead[t-1] + sigma*v_rep_ahead[t];
    s_rep_ahead[t] = exp((mu + h_rep_ahead[t]) / gamma);
    loglik_rep_ahead[t] = normal_lpdf(y_ahead[t] | mu_ret, s_rep_ahead[t]);
  }
}
