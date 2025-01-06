data {
  int<lower=0> T;   // # time points (equally spaced)
  vector[T] y;      // mean corrected return at time t
  
  int<lower=0> K;   // # number of AR time scales to model
}
transformed data {
  vector[K] ones = rep_vector(1.0, K);
}
parameters {
  real mu_ret;
  real rho;

  real mu;                             // mean log volatility
  vector<lower=0, upper=1>[K] phi_raw; // persistence of volatility
  positive_ordered[K] sigma;            // white noise shock scale
  
  matrix[T, K] v_raw;                      // std log volatility time t
}
transformed parameters {
  matrix[T, K] v = diag_post_multiply(v_raw, sigma);
  vector[K] phi;
  matrix[T, K] h;
  vector[T] m;
  vector[T] s;
  
  // Set phi to be ordered. Biggest phi is first
  phi[1] = -1 + (1 - -1) * phi_raw[1];
  for(k in 2:K) {
    phi[k] = -1 + (phi[k-1] - -1) * phi_raw[k];
  }
  
  // Set h similar to Stan manual of SV model
  h[1,] = v[1,] ./ sqrt(1 - phi .* phi)';     // rescale h[1]
  for (t in 2:T) {
    h[t,] = phi' .* h[t-1,] + v[t,];
  }
  
  s = exp((mu + h*ones) / 2);
}
model {
  // prior
  mu_ret ~ normal(0, 0.2);
  rho ~ normal(0, 1);

  // (phi+1)/2.0 ~ beta(8, 2);
  sigma ~ exponential(10);
  mu ~ normal(0, 1);
  
  // likelihood
  to_vector(v_raw) ~ std_normal();
  
  // only fastest biggest v affects skew
  y ~ normal(mu_ret + rho * v[,K] .* s, s);
}
generated quantities {
  vector[K] v_rep_ahead = sigma .* to_vector(normal_rng(rep_vector(0.0, K), 1));
  vector[K] h_rep_ahead = (phi' .* h[T,])' + v_rep_ahead;
  real s_rep_ahead = exp((mu + dot_product(h_rep_ahead, ones)) / 2);
  real y_rep_ahead = normal_rng(mu_ret + rho * v_rep_ahead[K]  * s_rep_ahead, s_rep_ahead);
  
  vector[T] y_rep;
  vector[T] log_lik;
  
  for(t in 1:T) {
   log_lik[t] = normal_lpdf(y[t] | mu_ret + rho * v[t,K] .* s[t], s[t]);
  }
  
  for(t in 1:T) {
    y_rep[t] = normal_rng(mu_ret + rho * v[t,K]  * s[t], s[t]);
  }
}
