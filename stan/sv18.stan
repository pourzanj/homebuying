data {
  int<lower=0> T;   // # time points (equally spaced)
  vector[T] y;      // mean corrected return at time t
  
  int<lower=2> M;
  vector[M] y_ahead;
}

parameters {
  real mu_ret;

  // SV
  real mu;                     // mean log volatility
  real<lower=0, upper=1> phi_raw[2];  // persistence of volatility
  real<lower=0> sigma[2];         // white noise shock scale
  // positive_ordered[2] sigma;
  
  vector[T] v1;             // std log volatility time t
  vector[T] v2; 
  // vector[M] v_rep_ahead;
  
}
transformed parameters {
  real phi[2];
  vector[T] h1;
  vector[T] h2;
  vector[T] h;
  vector[T] s;

  // Set phi to be ordered. Biggest phi is first
  // phi[1] = -1 + (1 - -1) * phi_raw[1];
  // phi[2] = -1 + (phi[1] - -1) * phi_raw[2];
  phi[1] = phi_raw[1];
  phi[2] = phi_raw[1] * phi_raw[2];

  h1[1] = sigma[1]*v1[1] / sqrt(1 - phi[1]^2);
  h2[1] = sigma[2]*v2[1] / sqrt(1 - phi[2]^2);
  h[1] = h1[1] + h2[1];
  
  for (t in 2:T) {
    h1[t] = phi[1]*h1[t-1] + sigma[1]*v1[t];
    h2[t] = h2[t-1] + sigma[2]*v2[t];
    h[t] = h1[t] + h2[t];
  }
  s = exp((mu + h1 + h2) / 2);
}
model {
  // prior
  // mu_ret ~ normal(0, 0.2);
  // rho ~ normal(-1, 1);

  // (phi+1)/2.0 ~ beta(8, 2);
  // sigma ~ exponential(10);
  // mu ~ normal(0, 1);
  
  // likelihood
  v1 ~ std_normal();
  v2 ~ std_normal();
  y ~ normal(mu_ret, s);
  
  // for LOO
  // v_rep_ahead ~ std_normal();
  
}
generated quantities {
  // vector[T] v_rep = to_vector(normal_rng(rep_array(0.0, T), 1));
  // vector[T] h_rep;
  // vector[T] s_rep;
  // vector[T] y_rep;
  // 
  // vector[T] log_lik;
  // real sum_log_lik;
  // vector[T] log_lik_ahead;
  // real sum_log_lik_ahead;
  // 
  // vector[M] h_rep_ahead;
  // vector[M] s_rep_ahead;
  // vector[M] loglik_rep_ahead;
  // 
  // h_rep[1] = sigma*v_rep[1] / sqrt(1 - phi^2);
  // 
  // for(t in 2:T) {
  //   h_rep[t] = phi*h[t-1] + sigma*v_rep[t];
  //   s_rep[t] = exp((mu+h_rep[t])/gamma);
  //   y_rep[t] = normal_rng(mu_ret, s_rep[t]);
  //   log_lik[t] = normal_lpdf(y[t] | mu_ret, s[t]);
  //   log_lik_ahead[t] = normal_lpdf(y[t] | mu_ret, s_rep[t]);
  // }
  // sum_log_lik = sum(log_lik[2:T]);
  // sum_log_lik_ahead = sum(log_lik_ahead[2:T]);
  // 
  // // Compute rep_ahead and loglik_rep_ahead
  // h_rep_ahead[1] = phi*h[T] + sigma*v_rep_ahead[1];
  // s_rep_ahead[1] = exp((mu + h_rep_ahead[1]) / gamma);
  // loglik_rep_ahead[1] = normal_lpdf(y_ahead[1] | mu_ret, s_rep_ahead[1]);
  // for (t in 2:M) {
  //   h_rep_ahead[t] = phi*h_rep_ahead[t-1] + sigma*v_rep_ahead[t];
  //   s_rep_ahead[t] = exp((mu + h_rep_ahead[t]) / gamma);
  //   loglik_rep_ahead[t] = normal_lpdf(y_ahead[t] | mu_ret, s_rep_ahead[t]);
  // }
}
