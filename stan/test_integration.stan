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
  
  // Put strike price
  real<lower=0> k_p[N_p];
  real<lower=0> k_c[N_c];
  
  real mu;
  real<lower=0> sigma;
  
  real<lower=0> a;
}
transformed data {
  int x_i[0];
}
parameters {
  real x;
}
transformed parameters {
  // Theoretical prices of puts and calls assuming underlying
  // return follows Gaussian and markets are efficient
  real p_t[N_p];
  real c_t[N_c];

  // print("mu: ", mu);
  // print("sigma: ", sigma);
  
  // Compute expected value of puts under distribution for every strike k
  for(k in 1:N_p) {
    p_t[k] = integrate_ode_rk45(put_normal_weighted_payout,
                                {0.0},
                                a, {k_p[k]},
                                {mu, sigma},
                                {k_p[k], s0}, x_i)[1,1];
                           
  }
  
  for(k in 1:N_c) {
    c_t[k] = integrate_ode_rk45(call_normal_weighted_payout,
                                {0.0},
                                k_c[k], {2*s0},
                                {mu, sigma},
                                {k_c[k], s0}, x_i)[1,1];
                           
  }
  
  
  
}
model {
  x ~ normal(0, 1);
}
