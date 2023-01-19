// Mixture distribution controlled by a latent AR1
data {
  int K;
  int<lower=0> T;   // # time points (equally spaced)
  vector[T] y;      // mean corrected return at time t
}

parameters {
  // AR1
  real<lower=-1, upper=1> phi;  // persistence of volatility
  real<lower=0> sigma;         // white noise shock scale
  
  // Relationship between likelihood parameters and latent AR1
  vector[K] a_theta; vector[K] b_theta;
  // simplex[K] theta;
  real mu[K];
  // vector[K] a_m; vector[K] b_m;
  ordered[K] a_s; // positive_ordered[K] b_s;
  
  // Random effects for latent AR1
  vector[T] v;             // std log volatility time t
}
transformed parameters {
  vector[T] h;
  
  matrix[T, K] theta_raw;
  simplex[K] theta[T];
  // matrix[T, K] m;
  // matrix[T, K] s;

  h[1] = sigma*v[1] / sqrt(1 - phi^2);
  for (t in 2:T) {
    h[t] = phi*h[t-1] + sigma*v[t];
  }
  
  for(k in 1:K) {
    theta_raw[,k] = a_theta[k] + b_theta[k]*h;
    // m[,k] = a_m[k] + b_m[k]*h;
    // s[,k] = exp((a_s[k] + h)/2.0);
  }
  for(t in 1:T) {
    theta[t] = softmax(theta_raw[t,]');
  }
}
model {
  // prior
  sigma ~ normal(0, 1);
  
  // theta[1] ~ normal(0.45, 0.05);
  // theta[2] ~ normal(0.45, 0.05);
  // theta[3] ~ normal(0.1, 0.01);
  
  mu ~ normal(0, 10);
  // mu[1] ~ normal(0.1, 1);
  // mu[2] ~ normal(0.0, 1);
  // mu[3] ~ normal(-0.4, 1);
  
  // a_s[1] ~ normal(2*-0.7, 0.1);
  // a_s[2] ~ normal(2*0.0, 0.1);
  // a_s[3] ~ normal(2*0.6, 0.1);
  
  a_theta ~ normal(0, 10); b_theta ~ normal(0, 10);
  // a_m ~ normal(0, 10); b_m ~ normal(0, 10);
  a_s ~ normal(0, 10);
  
  // likelihood
  v ~ std_normal();
  // y ~ normal(mu_ret + rho*v .* s, sqrt(1-rho^2)*s);
  
  for (t in 1:T) {
    // vector[K] lps = log(theta);
    vector[K] lps = log(theta[t]);
    for (k in 1:K) {
      // lps[k] += normal_lpdf(y[t] | m[t,k], s[t,k]);
      // lps[k] += normal_lpdf(y[t] | mu[k], s[t,k]);
      lps[k] += normal_lpdf(y[t] | mu[k], exp((a_s[k] + h)/2.0));
    }
    target += log_sum_exp(lps);
  }
  
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
}
