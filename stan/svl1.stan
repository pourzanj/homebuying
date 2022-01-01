data {
  int<lower=0> T;   // # time points (equally spaced)
  vector[T] y;      // mean corrected return at time t
}
parameters {
  real mu_ret;
  real<lower=-1, upper=1> rho;

  real<lower=0> alpha;
  real<lower=0> beta;

  // SV
  real mu;                     // mean log volatility
  vector<lower=0, upper=1>[T] phi;  // persistence of volatility
  real<lower=0> sigma;         // white noise shock scale
  
  vector[T] v;            
}
transformed parameters {
  vector[T] h;
  vector[T] s;
  
  h[1] = sigma * v[1] / sqrt(1 - phi[1]^2);
  
  for (t in 2:T) {
    h[t] = phi[t] * h[t-1] + sigma * v[t];
  }
  
  s = exp((mu + h) / 2);
}
model {
  // prior
  mu_ret ~ normal(0, 0.2);
  // rho ~ normal(-1, 1);

  // (phi+1)/2.0 ~ beta(8, 2);
  // sigma ~ exponential(10);
  // mu ~ normal(0, 1);
  
  alpha ~ normal(98, 1);
  beta ~ normal(2, 0.1);
  phi[1] ~ beta(alpha, beta);
  for (t in 1:T) {
    phi[t] ~ beta(phi[t-1] * 100, (1-phi[t-1]) * 100);
  }
  
  // likelihood
  for (t in 1:T) {
    matrix[2,2] S = [[square(s[t]), s[t]*rho],[s[t]*rho, 1]];
    [y[t], v[t]] ~ multi_normal([mu_ret, 0], S);
  }
  
}
generated quantities {
  vector[T] v_rep = to_vector(normal_rng(rep_array(0, T), 1));
  vector[T] h_rep = v_rep;
  vector[T] s_rep;
  vector[T] y_rep;
  
  h_rep[1] = sigma * v_rep[1] / sqrt(1 - phi[1]^2);
  for(t in 2:T) h_rep[t] = phi[t] * h_rep[t-1] + sigma * v_rep[t];
  s_rep = exp((mu + h_rep) / 2);
  y_rep = mu_ret + rho * (s_rep .* v_rep) + sqrt(1-rho^2) * (s_rep .* to_vector(normal_rng(rep_array(0, T), 1)));
}
