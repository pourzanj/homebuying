functions {
  real[] put_normal_weighted_payout(real t, real[] y,
                                    real[] theta,
                                    real[] x_r, int[] x_i) {

    // Unpack parameters
    real mu = theta[1];
    real sigma = theta[2];
    
    // Unpack data
    real K = x_r[1];
    real x0 = x_r[2];

    // Transform underlying price to return
    // real r = 100 * (t - x0)/x0;
    real r = 100 * (log(t) - log(x0));

    // Compute terms of integrand and return product
    real return_dens = exp(normal_lpdf(r | mu, sigma));
    // real trans_adj = 100/x0;
    real trans_adj = 100 / t;
    real option_payoff = K - t;

    return {return_dens * trans_adj * option_payoff};
  }
  
  real[] call_normal_weighted_payout(real t, real[] y,
                                     real[] theta,
                                     real[] x_r, int[] x_i) {

    // Unpack parameters
    real mu = theta[1];
    real sigma = theta[2];
    
    // Unpack data
    real K = x_r[1];
    real x0 = x_r[2];

    // Transform underlying price to return
    // real r = 100 * (t - x0)/x0;
    real r = 100 * (log(t) - log(x0));

    // Compute terms of integrand and return product
    real return_dens = exp(normal_lpdf(r | mu, sigma));
    // real trans_adj = 100/x0;
    real trans_adj = 100 / t;
    real option_payoff = t - K;

    return {return_dens * trans_adj * option_payoff};
  }
}
data {
  
  // Current price of underlying asset
  real<lower=0> s0;
  
  // Number of observed puts and calls
  int<lower=0> N_p;
  int<lower=0> N_c;
  
  // Strike prices for puts and calls
  vector<lower=0>[N_p] k_p;
  vector<lower=0>[N_c] k_c;
  
  // Observed prices of puts and calls
  vector<lower=0>[N_p] p_o;
  vector<lower=0>[N_c] c_o;
  
  // Bid ask spreads
  vector<lower=0>[N_p] p_spread;
  vector<lower=0>[N_c] c_spread;
}
transformed data {
  vector[N_p] log_p_o = log(p_o);
  vector[N_c] log_c_o = log(c_o);
  
  int x_i[0];
}
parameters {
  // Percent returns distribution
  real mu;
  real<lower=0> sigma;
  
  // Random effects
  real<lower=0> tau;
  vector[N_p] log_p_e;
  vector[N_c] log_c_e;

}
transformed parameters {
  // Theoretical prices of puts and calls assuming underlying
  // return follows Gaussian and markets are efficient
  vector[N_p] p_t;
  vector[N_c] c_t;
  
  // Compute expected value of puts under distribution for every strike k
  for(k in 1:N_p) {
    p_t[k] = integrate_ode_rk45(put_normal_weighted_payout,
                                {0.0},
                                1.0, {k_p[k]},
                                {mu, sigma},
                                {k_p[k], s0}, x_i)[1,1];
                           
  }
  
  // Compute expected value of calls under distribution for every strike k
  for(k in 1:N_c) {
    c_t[k] = integrate_ode_rk45(call_normal_weighted_payout,
                                {0.0},
                                k_c[k], {2*s0},
                                {mu, sigma},
                                {k_c[k], s0}, x_i)[1,1];
                           
  }
  
}
model {
  // Priors
  mu ~ normal(0, 1);
  sigma ~ normal(0, 10);
  tau ~ normal(0, 10);
  
  // Expected prices are theoretical price plus random noise to absorb
  // either mis-specification of the underlying return distribution or
  // non-efficient pricing
  log_p_e ~ normal(log(p_t), tau);
  log_c_e ~ normal(log(c_t), tau);
  
  // Observed prices are a noisy observation of expected prices but should
  // be within the bid ask spread
  log_p_o ~ normal(log_p_e, p_spread/4);
  log_c_o ~ normal(log_c_e, c_spread/4);
}
generated quantities {
  vector[N_p] log_p_error = log_p_e - log(p_t);
  vector[N_c] log_c_error = log_c_e - log(c_t);
  real<lower=0> p_o_rep[N_p] = exp(normal_rng(log_p_e, p_spread/4));
  real<lower=0> c_o_rep[N_c] = exp(normal_rng(log_c_e, c_spread/4));
}
