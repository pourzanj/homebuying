data {
  int<lower=0> T;   // # time points (equally spaced)
  vector[T] y;      // mean corrected return at time t
  
  real tau0;
  real nu;
  real scale;
}

transformed data {
  real<lower=-1, upper=1> rho = 0.0;
}

parameters {
  real mu_ret;
  // real<lower=-1, upper=1> rho;
  real alpha[2];
  real beta[2];
  
  real tau[2];
  // real<lower=0> tau;
  real<lower=0> c2;

  // SV
  real mu;                     // mean log volatility
  real<lower=-1, upper=1> phi;  // persistence of volatility
  real<lower=0> sigma;         // white noise shock scale
  
  vector[T] v;             // std log volatility time t
  
  vector<lower=0>[T] l;
}
transformed parameters {
  vector[T] l2 = l .* l;
  // vector[T] lt = sqrt((c2 * l2) ./ (c2 + square(tau) * l2));
  vector[T] lt;
  
  vector[T] h;
  vector[T] m;
  vector[T] a;
  vector[T] s;
  vector[T] u;

  h[1] = sigma*v[1] / sqrt(1 - phi^2);
  m[1] = mu_ret;
  a[1] = 0.0;
  u[1] = 1;
  lt[1] = 1;
  
  for (t in 2:T) {
    // m[t] = mu_ret;
    m[t] = mu_ret + beta[1]*h[t-1] + beta[2]*(h[t-1]^2-1);
    // a[t] = 0;
    // a[t] = alpha[1] + alpha[2]*h[t-1] + alpha[3]*(h[t-1]^2-1);
    a[t] = alpha[1] + alpha[2]*h[t-1];
    h[t] = phi*h[t-1] + sigma*v[t];
    // u[t] = exp(tau[1] + tau[2]*h[t-1] + tau[3]*h[t-1]^2-1);
    u[t] = exp(tau[1] + tau[2]*h[t-1] + tau[3]*h[t-1]^2-1);
    // u[t] = tau;
    // lt[t] = sqrt((c2 * l2[t]) ./ (c2 + square(u[t]) * l2[t]));
    s[t] = exp((mu + h[t]) / 2.0);
    lt[t] = sqrt((s[t]^2 * l2[t]) ./ (s[t]^2 + square(u[t]) * l2[t]));
  }
  
  // s = exp((mu + h) / 2.0);
}
model {
  // prior
  // mu_ret ~ normal(0, 0.2);
  // rho ~ normal(-1, 1);

  // (phi+1)/2.0 ~ beta(8, 2);
  // sigma ~ exponential(10);
  // mu ~ normal(0, 1);
  
  // tau ~ cauchy(0, tau0);
  tau[1] ~ normal(2, 0.1);
  tau[2] ~ normal(0, 1);
  tau[3] ~ normal(0, 0.1);
  c2 ~ inv_gamma(0.5 * nu, 0.5 * nu * scale^2);
  
  // likelihood
  v ~ std_normal();
  l ~ cauchy(0, 1);
  // y ~ skew_normal(m, s .* (u .* lt), a);
  y ~ skew_normal(m, u .* lt, a);
  
}
generated quantities {
  vector[T] v_rep = to_vector(normal_rng(rep_array(0.0, T), 1));
  vector[T] l2_rep = to_vector(cauchy_rng(rep_array(0.0, T), 1))^2;
  vector[T] lt_rep;
  vector[T] h_rep;
  vector[T] s_rep;
  vector[T] y_rep;

  vector[T] log_lik;
  real sum_log_lik;
  vector[T] log_lik_ahead;
  real sum_log_lik_ahead;

  h_rep[1] = sigma*v_rep[1] / sqrt(1 - phi^2);
  lt_rep[1] = 1;

  for(t in 2:T) {
    h_rep[t] = phi*h[t-1] + sigma*v_rep[t];
    s_rep[t] = exp((mu+h_rep[t])/2.0);
    lt_rep[t] = sqrt((s_rep[t]^2 * l2_rep[t]) ./ (s_rep[t]^2 + square(u[t]) * l2_rep[t]));
    
    y_rep[t] = skew_normal_rng(m[t], u[t]*lt_rep[t], a[t]);
    log_lik[t] = skew_normal_lpdf(y[t] | m[t], u[t]*lt_rep[t], a[t]);
    log_lik_ahead[t] = skew_normal_lpdf(y[t] | m[t], u[t]*lt_rep[t], a[t]);
  }
  sum_log_lik = sum(log_lik[2:T]);
  sum_log_lik_ahead = sum(log_lik_ahead[2:T]);
}
