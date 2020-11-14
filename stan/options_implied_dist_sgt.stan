functions {
  
  real get_v(real lambda, real p, real q) {
    return q^(-1.0/p) * ((3.0*lambda^2+1)*exp(lbeta(3.0/p,q-2.0/p)-lbeta(1.0/p,q)) - 4*lambda^2*exp(lbeta(2.0/p,q-1.0/p)-lbeta(1.0/p,q))^2)^(-1.0/2.0);
  }
  
  real get_m(real sigma, real lambda, real p, real q, real v) {
    return 2*v*sigma*lambda*q^(1.0/p) * exp(lbeta(2.0/p,q-1.0/p)-lbeta(1.0/p,q));
  }
  
  real sgt_lpdf(real x, real mu, real s, real l, real p, real q, real v, real m)
  { // Skewed generalised t
    real lz1;
    real lz2;
    real r;
    real out;
    real acc;
    real vs;
    int sign;
    lz1 = lbeta(1.0/p,q);
    lz2 = lbeta(2.0/p,q-1.0/p);
    vs = v*s;

    out = 0;
    out += log(p);
    out -= log(2*v);
    out -= log(s);
    out -= (1 / p) * log(q);
    out -= lz1;

    r = x-mu+m;

    acc = 0.0;

    if(r < 0) sign = -1;
    else sign = 1;

    acc += log1p((fabs(r) / (vs*(l*sign + 1)))^p / q);

    out -= (1.0/p+q) * acc;
    
    return out;
  }
  
  real[] put_sgt_weighted_payout(real t, real[] y,
                                    real[] theta,
                                    real[] x_r, int[] x_i) {

    // Unpack parameters
    real mu = theta[1];
    real sigma = theta[2];
    real lambda = theta[3];
    real p = theta[4];
    real q = theta[5];
    real v = theta[6];
    real m = theta[7];
    
    // Unpack data
    real K = x_r[1];
    real x0 = x_r[2];

    // Transform underlying price to return
    // real r = 100 * (t - x0)/x0;
    real r = 100 * (log(t) - log(x0));

    // Compute terms of integrand and return product
    real return_dens = exp(sgt_lpdf(r | mu, sigma, lambda, p, q, v, m));
    // real trans_adj = 100/x0;
    real trans_adj = 100 / t;
    real option_payoff = K - t;

    return {return_dens * trans_adj * option_payoff};
  }
  
  real[] call_sgt_weighted_payout(real t, real[] y,
                                     real[] theta,
                                     real[] x_r, int[] x_i) {

    // Unpack parameters
    real mu = theta[1];
    real sigma = theta[2];
    real lambda = theta[3];
    real p = theta[4];
    real q = theta[5];
    real v = theta[6];
    real m = theta[7];
    
    // Unpack data
    real K = x_r[1];
    real x0 = x_r[2];

    // Transform underlying price to return
    // real r = 100 * (t - x0)/x0;
    real r = 100 * (log(t) - log(x0));

    // Compute terms of integrand and return product
    real return_dens = exp(sgt_lpdf(r | mu, sigma, lambda, p, q, v, m));
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
  
  // Hyper-parameters
  real m_lambda; real<lower=0> s_lambda;
  real m_p; real<lower=0> s_p;
  real m_q; real<lower=0> s_q;
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
  real<lower=-1, upper = 1> lambda;
  real<lower=1> p;
  real<lower=2> q;
  
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
  real v = get_v(lambda, p, q);
  real m = get_m(sigma, lambda, p, q, v);
  
  // Compute expected value of puts under distribution for every strike k
  for(k in 1:N_p) {
    p_t[k] = integrate_ode_rk45(put_sgt_weighted_payout,
                                {0.0},
                                1.0, {k_p[k]},
                                {mu, sigma, lambda, p, q, v, m},
                                {k_p[k], s0}, x_i)[1,1];
                                
    if((p_t[k] <= 0.0) || is_nan(p_t[k])) {
      print("p_t[k]: ", p_t[k], " k: ", k, " k_p[k]: ", k_p[k], " mu: " , mu, " sigma: " , sigma, " p: ", p, " q: ", q);
    }
                           
  }
  
  // Compute expected value of calls under distribution for every strike k
  for(k in 1:N_c) {
    c_t[k] = integrate_ode_rk45(call_sgt_weighted_payout,
                                {0.0},
                                k_c[k], {2*s0},
                                {mu, sigma, lambda, p, q, v, m},
                                {k_c[k], s0}, x_i)[1,1];
                                
    if((c_t[k] <= 0.0) || is_nan(c_t[k])) {
      print("c_t[k]: ", c_t[k], " k: ", k, " k_c[k]: ", k_c[k], " mu: " , mu, " sigma: " , sigma, " p: ", p, " q: ", q);
    }
                           
  }
}
model {
  // Priors
  mu ~ normal(0, 1);
  sigma ~ normal(0, 10);
  lambda ~ normal(m_lambda, s_lambda);
  p ~ normal(m_p, s_p);
  q ~ normal(m_q, s_q);
  
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
  real pq = p*q;
}
