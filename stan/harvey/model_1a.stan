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
    // return x[4];
    // return log_p[4];
    // return -lambda[4];
    // return main_term[4];
  }
  
  real u(real x, real mu, real lambda, real v, real nu, real gamma) {
    real eps = (x - mu) / exp(lambda) + expec_fsgt(v, nu, gamma);
    real mu_eps = expec_fsgt(v, nu, gamma);
    return((nu*v/2 + 1) * (1 - mu_eps / eps) * (fabs(eps)^v / (fabs(eps)^v + nu*gamma^(v*sign(eps)))) - 1);
  }
}
data {
  int<lower=0> T;
  vector[T] y;
  
  real m_omega; real<lower=0> s_omega;
  real a_phi; real<lower=0> b_phi;
  real m_k1; real<lower=0> s_k1;
  real m_k2; real<lower=0> s_k2;
  
  real m_lambda1; real<lower=0> s_lambda1;
  
  real m_mu; real<lower=0> s_mu;
  real m_v; real<lower=0> s_v;
  real m_nu; real<lower=0> s_nu;
  real m_gamma; real<lower=0> s_gamma;
}
transformed data {
  // real omega = -0.3;
  // real<lower=0, upper=1> phi = 0.96;
  // real k1 = 0.038;
  // real k2 = 0.02;
  // real lambda1 = -1;
  // real mu = 0;
  // real<lower=1> v = 2;
  // real<lower=2> nu = 10;
  // real<lower=0> gamma = 1;
}
parameters {
  
  // Scale parameters
  real omega;
  real<lower=0, upper=1> phi;
  real k1;
  real k2;

  // Scale of very first observation
  real lambda1;
  
  // Other FSGT Parameters
  real mu;
  real<lower=1> v;
  real<lower=2> nu;
  real<lower=0> gamma;
}
transformed parameters {
  vector[T] Eeps;
  vector[T] eps;
  vector[T] u_t;
  vector[T] h;
  vector[T] lambda;
  
  lambda[1] = lambda1;
  Eeps[1] = expec_fsgt(v, nu, gamma);
  eps[1] = (y[1] - mu)/exp(lambda[1]) + Eeps[1];
  u_t[1] = u(y[1], mu, lambda[1], v, nu, gamma);
  h[1] = k1*u_t[1] + k2*sign(-eps[1])*(u_t[1]+1);
  
  for (t in 2:T) {
    lambda[t] = omega*(1-phi) + phi*lambda[t-1] + h[t-1];
    
    Eeps[t] = expec_fsgt(v, nu, gamma);
    eps[t] = (y[t] - mu)/exp(lambda[t]) + Eeps[t];
    u_t[t] = u(y[t], mu, lambda[t], v, nu, gamma);
    h[t] = k1*u_t[t] + k2*sign(-eps[t])*(u_t[t]+1);
  }
}
model {
  // prior
  omega ~ normal(m_omega, s_omega);
  phi ~ beta(a_phi, b_phi);
  k1 ~ normal(m_k1, s_k1);
  k2 ~ normal(m_k2, s_k2);

  lambda1 ~ normal(m_lambda1, s_lambda1);

  mu ~ normal(m_mu, s_mu);
  v ~ normal(m_v, s_v);
  nu ~ normal(m_nu, s_nu);
  gamma ~ normal(m_gamma, s_gamma);
  
  // likelihood
  y ~ fsgt(mu, lambda, v, nu, gamma);
  // target += lambda[4];
  // target += lambda[3];
  // target += h[4];
  // target += u_t[4];
  // target += sign(eps[3]);
}
generated quantities {
  // measure of heavy-tailness
  real vnu2 = v*nu / 2;
}
