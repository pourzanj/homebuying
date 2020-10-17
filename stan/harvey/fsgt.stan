// Generalized T (GT) with Skew added using Fernandez and Steel (1998) method
functions {
  
  real sign(real x) {
    if(x < 0.0) {
      return -1.0;
    } else {
      return 1.0;
    }
  }
  
  vector sign_vec(vector x) {
    int N = num_elements(x);
    vector[N] ret;
    for(n in 1:N) {
      if(x[n] < 0.0) {
        ret[n] = -1.0;
      } else {
        ret[n] = 1.0;
      }
    }
    return ret;
  }
  
  vector pow_vec(vector x, real y) {
    int N = num_elements(x);
    vector[N] ret;
    for(n in 1:N) {
      ret[n] = x[n]^y;
    }
    return ret;
  }
  
  vector pow_real_vec(real x, vector y) {
    int N = num_elements(y);
    vector[N] ret;
    for(n in 1:N) {
      ret[n] = x^y[n];
    }
    return ret;
  }
  
  real expec_fsgt(real v, real nu, real gamma) {
    return ((gamma - 1/gamma) * ((tgamma(2/v)*tgamma(nu/2-1/v)) / (tgamma(1/v)*tgamma(nu/2))) * nu^(1/v));
  }

  real fsgt_lpdf(vector x, real mu, vector lambda, real v, real nu, real gamma) {
    int N = num_elements(x);
    vector[N] eps = (x - mu) ./ exp(lambda) + expec_fsgt(v, nu, gamma);
    real term1 = log(2) - log(gamma + 1/gamma);
    real term2 = log(v) - log(2*nu^(1/v));
    real term3 = -lbeta(nu/2, 1/v);
    vector[N] main_term = -(nu/2 + 1/v) * log(1 + (1/nu) * pow_vec(fabs(eps), v) ./ pow_real_vec(gamma, v * sign_vec(eps)));

    vector[N] log_p = -lambda + term1 + term2 + term3 + main_term;

    return sum(log_p);
  }
  
}
data {
  int<lower=0> N;
  vector[N] y;
  
  real m_mu; real<lower=0> s_mu;
  real m_lambda; real<lower=0> s_lambda;
  real m_v; real<lower=0> s_v;
  real m_nu; real<lower=0> s_nu;
  real m_gamma; real<lower=0> s_gamma;
}
transformed data {
  // real mu = 0;
  // real<lower=1> v = 2;
  // real<lower=2> nu = 20;
  // real<lower=0> gamma = 1;
}
parameters {
  real mu;
  real lambda;
  real<lower=1> v;
  real<lower=2> nu;
  real<lower=0> gamma;
}
model {
  // prior
  mu ~ normal(m_mu, s_mu);
  lambda ~ normal(m_mu, s_mu);
  v ~ normal(m_v, s_v);
  nu ~ normal(m_nu, s_nu);
  gamma ~ normal(m_gamma, s_gamma);
  
  // likelihood
  y ~ fsgt(mu, rep_vector(lambda, N), v, nu, gamma);
}

