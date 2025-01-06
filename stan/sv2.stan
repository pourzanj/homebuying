// Added time-varying rho
data {
  int<lower=0> T;   // # time points (equally spaced)
  vector[T] y;      // mean corrected return at time t
}
parameters {
  real mu_ret;
  // real rho;

  // SV
  real mu;                     // mean log volatility
  real<lower=-1, upper=1> phi;  // persistence of volatility
  real<lower=0> sigma;         // white noise shock scale
  
  // Rho
  real mu_rho;
  real<lower=-1, upper=1> phi_rho;
  real<lower=0> sigma_rho;
  
  vector[T] v;             // std log volatility time t
  vector[T] w;
  vector<lower=0>[T] wl;
}
transformed parameters {
  vector[T] h;
  vector[T] s;
  
  vector[T] f;
  vector[T] rho;

  h[1] = sigma*v[1] / sqrt(1 - phi^2);
  f[1] = sigma_rho*w[1] / sqrt(1 - phi_rho^2);
  
  for (t in 2:T) {
    h[t] = phi*h[t-1] + sigma*v[t];
    // f[t] = phi_rho*f[t-1] + sigma_rho*w[t];
    f[t] = f[t-1] + sigma_rho*w[t]*wl[t];
  }
  
  s = exp((mu + h) / 2.0);
  rho = -1.0 + 2.0*inv_logit(mu_rho + f);
}
model {
  // prior
  // mu_ret ~ normal(0, 0.2);
  // rho ~ normal(-1, 1);

  // (phi+1)/2.0 ~ beta(8, 2);
  // sigma ~ exponential(10);
  // mu ~ normal(0, 1);
  
  mu_rho ~ normal(0.0, 1.0);
  (phi_rho+1)/2.0 ~ beta(495, 5);
  // sigma_rho ~ lognormal(log(0.05), 1.0);
  sigma_rho ~ normal(0.0, 0.5);
  
  // likelihood
  v ~ std_normal();
  w ~ std_normal();
  wl ~ student_t(3, 0, 1);
  y ~ normal(mu_ret + rho .* v .* s, sqrt(1-rho^2) .* s);
  
}
generated quantities {
  vector[T] v_rep = to_vector(normal_rng(rep_array(0.0, T), 1));
  vector[T] h_rep;
  vector[T] s_rep;
  
  vector[T] w_rep = to_vector(normal_rng(rep_array(0.0, T), 1));
  vector[T] f_rep;
  vector[T] rho_rep;
  
  vector[T] y_rep;

  vector[T] log_lik;
  real sum_log_lik;

  h_rep[1] = sigma*v_rep[1] / sqrt(1 - phi^2);
  f_rep[1] = sigma_rho*w_rep[1] / sqrt(1 - phi_rho^2);

  for(t in 2:T) {
    h_rep[t] = phi*h[t-1] + sigma*v_rep[t];
    s_rep[t] = exp((mu+h_rep[t])/2.0);
    
    f_rep[t] = phi_rho*f[t-1] + sigma_rho*w_rep[t];
    rho_rep[t] = -1.0 + 2.0*inv_logit(mu_rho + f_rep[t]);
    
    y_rep[t] = normal_rng(mu_ret + rho_rep[t]*v_rep[t]*s_rep[t], sqrt(1-rho_rep[t]^2)*s_rep[t]);
    log_lik[t] = normal_lpdf(y[t] | mu_ret + rho_rep[t]*v_rep[t]*s_rep[t], sqrt(1-rho_rep[t]^2)*s_rep[t]);
  }
  sum_log_lik = sum(log_lik[2:T]);
}
