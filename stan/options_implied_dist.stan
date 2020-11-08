functions {
  real call_price_normal(real x,        // Function argument
                         real xc,       // Complement of function argument
                                        //  on the domain (defined later)
                         real[] theta,  // parameters
                         real[] x_r,    // data (real)
                         int[] x_i) {   // data (integer)
    // Unpack parameters
    real mu = theta[1];
    real sigma = theta[2];
    
    // Unpack data
    real K = x_r[1];
    real x0 = x_r[2];
  
    // Transform underlying price to return
    real r = 100 * (x - x0)/x0;
    
    // Compute terms of integrand and return product
    real return_dens = exp(normal_lpdf(r | mu, sigma));
    real trans_adj = 100/x0;
    real option_payoff = x - K;
    
    return return_dens * trans_adj * option_payoff;
  }
  
  real put_price_normal(real x,        // Function argument
                        real xc,       // Complement of function argument
                                       //  on the domain (defined later)
                        real[] theta,  // parameters
                        real[] x_r,    // data (real)
                        int[] x_i) {   // data (integer)
    // Unpack parameters
    real mu = theta[1];
    real sigma = theta[2];
    
    // Unpack data
    real K = x_r[1];
    real x0 = x_r[2];
  
    // Transform underlying price to return
    real r = 100 * (x - x0)/x0;
    
    // Compute terms of integrand and return product
    real return_dens = exp(normal_lpdf(r | mu, sigma));
    real trans_adj = 100/x0;
    real option_payoff = K - x;
    
    return return_dens * trans_adj * option_payoff;
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
}
transformed data {
  vector<lower=0>[N_p] log_p_o = log(p_o);
  vector<lower=0>[N_c] log_c_o = log(c_o);
  
  int x_i[0];
  real x_r[0];
}
parameters {
  // Percent returns distribution
  real mu;
  real<lower=0> sigma;
  
  // Observational model of options prices
  real<lower=0> sigma_o;
}
transformed parameters {
  // Expected priced of puts and calls
  vector<lower=0>[N_p] log_p_e;
  vector<lower=0>[N_p] log_c_e;
  
  real theta[2] = { mu, sigma };
  
  // Compute expected value of puts under distribution for every strike k
  for(k in 1:N_p) {
    // log_p_e[k] = integrate_1d(put_price_normal, 0, k_p[k], { mu, sigma }, x_r, x_i);
    log_p_e[k] = log(integrate_1d(call_price_normal,
                             1.0,
                             positive_infinity(),
                             theta, x_r, x_i));
  }
}
model {
  // Observed prices are a noisy observation of expected prices
  // under the distribution
  log_p_o ~ normal(log_p_e, sigma_o);
  log_c_o ~ normal(log_c_e, sigma_o);
}

