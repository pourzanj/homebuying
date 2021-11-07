// A latent volatility for every day that can only go up according to a half
// horseshoe and only goes down according to the daily decay. 
data {
  int<lower=0> T;   // # time points (equally spaced)
  vector[T] y;      // mean corrected return at time t
}
parameters {
  // real mu_ret;
  real rho;
  real mu_ret;

  real<lower=0> alpha;
  real<lower=0> beta;

  // SV
  real<lower=0, upper=1> phi[T];  // persistence of volatility
  real<lower=0> sigma;         // white noise shock scale
  
  vector<lower=0>[T] v;
  vector<lower=0>[T] l;  
}
transformed parameters {
  vector<lower=0>[T] vls = sigma * (l .* v);
  vector<lower=0>[T] s;
  
  s[1] = vls[1];
  for (t in 2:T) {
    s[t] = phi[t] * s[t-1] + vls[t];
  }
  
}
model {
  // prior
  // mu_ret ~ normal(0, 0.2);
  // rho ~ normal(-1, 1);

  // (phi+1)/2.0 ~ beta(8, 2);
  // sigma ~ exponential(10);
  // mu ~ normal(0, 1);
  
  alpha ~ normal(92, 50);
  beta ~ normal(8, 0.1);
  phi ~ beta(alpha, beta);
  
  // likelihood
  l ~ cauchy(0, 1);
  v ~ std_normal();
  y ~ normal(rho * (vls .* s) + mu_ret, sqrt(1-rho^2)*s);
  // for (t in 1:T) {
  //   matrix[2,2] S = [[square(s[t]), s[t]*rho],[s[t]*rho, 1]];
  //   [y[t], v[t]] ~ multi_normal([0, 0], S);
  // }
  
}
generated quantities {
  vector[T] v_rep = to_vector(fabs(normal_rng(rep_array(0, T), 1)));
  vector[T] l_rep = to_vector(fabs(cauchy_rng(rep_array(0, T), 1)));
  vector[T] vls_rep;
  vector[T] s_rep;
  vector[T] y_rep;
  for(t in 1:T) l_rep[t] = min([l_rep[t], 50]);
  vls_rep = sigma * (v_rep .* l_rep);
  
  s_rep[1] = vls_rep[1];
  for(t in 2:T) s_rep[t] = phi[t] * s_rep[t-1] + vls_rep[t];
  
  y_rep = (rho * (vls_rep .* s_rep) + mu_ret) + sqrt(1-rho^2) * (s_rep .* to_vector(normal_rng(rep_array(0, T), 1)));
}
