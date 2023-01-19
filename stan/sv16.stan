data {
  int<lower=0> T;   // # time points (equally spaced)
  vector[T] y;      // mean corrected return at time t
  
  real tau0;
  real nu;
  real scale;
}

transformed data {
  // real<lower=0> tau = tau0;
}

parameters {
  real mu_ret;

  real<lower=0> tau;
  real<lower=0> c2;

  // SV
  real mu0;
  // real mu;                     // mean log volatility
  real<lower=-1, upper=1> phi;  // persistence of volatility
  real<lower=0> sigma;         // white noise shock scale
  
  vector[T] v;             // std log volatility time t
  
  vector[T] u_raw;             // std log volatility time t
  vector<lower=0>[T] l;
}
transformed parameters {
  vector[T] l2 = l .* l;
  vector[T] lt = sqrt((c2 * l2) ./ (c2 + square(tau) * l2));
  vector[T] u = u_raw .* (tau * lt); 
  
  vector[T] mu;
  vector[T] h;
  vector[T] s;

  mu[1] = mu0 + u[1];
  h[1] = sigma*v[1] / sqrt(1 - phi^2);
  
  for (t in 2:T) {
    mu[t] = mu[t-1] + u[t];
    h[t] = (1-phi)*mu[t]  + phi*h[t-1] + sigma*v[t];
  }
  s = exp(h / 2.0);
}
model {
  // prior
  // mu_ret ~ normal(0, 0.2);
  // rho ~ normal(-1, 1);

  // (phi+1)/2.0 ~ beta(8, 2);
  // sigma ~ exponential(10);
  // mu ~ normal(0, 1);
  mu0 ~ normal(0, 1);
  tau ~ cauchy(0, tau0);
  c2 ~ inv_gamma(0.5 * nu, 0.5 * nu * scale^2);
  
  // likelihood
  u_raw ~ std_normal();
  l ~ cauchy(0, 1);
  
  v ~ std_normal();
  y ~ normal(mu_ret, s);
  
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
  //   s_rep[t] = exp((mu+h_rep[t])/2.0);
  //   y_rep[t] = normal_rng(mu_ret + rho*v_rep[t]*s_rep[t], sqrt(1-rho^2)*s_rep[t]);
  //   log_lik[t] = normal_lpdf(y[t] | mu_ret + rho*v[t]*s[t], sqrt(1-rho^2)*s[t]);
  //   log_lik_ahead[t] = normal_lpdf(y[t] | mu_ret + rho*v_rep[t]*s_rep[t], sqrt(1-rho^2)*s_rep[t]);
  // }
  // sum_log_lik = sum(log_lik[2:T]);
  // sum_log_lik_ahead = sum(log_lik_ahead[2:T]);
  // 
  // // Compute rep_ahead and loglik_rep_ahead
  // h_rep_ahead[1] = phi*h[T] + sigma*v_rep_ahead[1];
  // s_rep_ahead[1] = exp((mu + h_rep_ahead[1]) / 2.0);
  // loglik_rep_ahead[1] = normal_lpdf(y_ahead[1] | mu_ret + rho*v_rep_ahead[1]*s_rep_ahead[1], sqrt(1-rho^2)*s_rep_ahead[1]);
  // for (t in 2:M) {
  //   h_rep_ahead[t] = phi*h_rep_ahead[t-1] + sigma*v_rep_ahead[t];
  //   s_rep_ahead[t] = exp((mu + h_rep_ahead[t]) / 2.0);
  //   loglik_rep_ahead[t] = normal_lpdf(y_ahead[t] | mu_ret + rho*v_rep_ahead[t]*s_rep_ahead[t], sqrt(1-rho^2)*s_rep_ahead[t]);
  // }
}
