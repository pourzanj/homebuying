// Add integrated likelihood
functions {
  real intnormal_lpdf(real y,
                      real h1,
                      real mu_ret, real mu, real phi, real sigma, real rho, 
                      vector grid) {
    int G = size(grid);
    real dp = 1.0/G;
    vector[G] lp;
    for(g in 1:G) {
      real v = grid[g];
      real h = phi*h1 + sigma*v;
      real s = exp((mu + h) / 2.0);
      
      lp[g] = normal_lpdf(y | mu_ret + rho*v .* s, sqrt(1-rho^2)*s) + log(dp);
    }
    return log_sum_exp(lp);
  }
}
data {
  int<lower=0> T;   // # time points (equally spaced)
  vector[T] y;      // mean corrected return at time t
  int G;
  vector[G] grid;
}

transformed data {
  // real<lower=-1, upper=1> rho = 0.0;
}

parameters {
  real mu_ret;
  real<lower=-1, upper=1> rho;

  // SV
  real mu;                     // mean log volatility
  real<lower=-1, upper=1> phi;  // persistence of volatility
  real<lower=0> sigma;         // white noise shock scale
  
  vector[T] v;             // std log volatility time t
}
transformed parameters {
  vector[T] h;
  vector[T] s;

  h[1] = sigma*v[1] / sqrt(1 - phi^2);
  
  for (t in 2:T) {
    h[t] = phi*h[t-1] + sigma*v[t];
  }
  
  s = exp((mu + h) / 2.0);
}
model {
  // prior
  // mu_ret ~ normal(0, 0.2);
  // rho ~ normal(-1, 1);

  // (phi+1)/2.0 ~ beta(8, 2);
  // sigma ~ exponential(10);
  // mu ~ normal(0, 1);
  
  // likelihood
  v ~ std_normal();
  y ~ normal(mu_ret + rho*v .* s, sqrt(1-rho^2)*s);
  for(t in 2:T) {
    y[t] ~ intnormal(h[t-1], mu_ret, mu, phi, sigma, rho, grid);
  }
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

  h_rep[1] = sigma*v_rep[1] / sqrt(1 - phi^2);
  
  for(t in 2:T) {
    h_rep[t] = phi*h[t-1] + sigma*v_rep[t];
    s_rep[t] = exp((mu+h_rep[t])/2.0);
    y_rep[t] = normal_rng(mu_ret + rho*v_rep[t]*s_rep[t], sqrt(1-rho^2)*s_rep[t]);
    log_lik[t] = normal_lpdf(y[t] | mu_ret + rho*v[t]*s[t], sqrt(1-rho^2)*s[t]);
    log_lik_ahead[t] = normal_lpdf(y[t] | mu_ret + rho*v_rep[t]*s_rep[t], sqrt(1-rho^2)*s_rep[t]);
  }
  sum_log_lik = sum(log_lik[2:T]);
  sum_log_lik_ahead = sum(log_lik_ahead[2:T]);
}
