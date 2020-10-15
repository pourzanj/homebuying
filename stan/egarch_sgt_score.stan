functions {
  real sgt_lpdf(vector x, real mu, vector s, real l, real p, real q)
  { // Skewed generalised t
    int N = dims(x)[1];
    real lz1;
    real lz2;
    real v;
    real m;
    vector[N] r;
    real out;
    real acc;
    vector[N] vs;
    int sign[N];
    lz1 = lbeta(1.0/p,q);
    lz2 = lbeta(2.0/p,q-1.0/p);
    v = 1;
    vs = v*s;

    out = 0;
    out += log(p);
    out -= log(2*v);
    out -= (1 / p) * log(q);
    out -= lz1;
    out = N*out;
    out -= sum(log(s));
    
    m = 0;
    r = x-mu;
    
    acc = 0.0;
    
    for (n in 1:N) {
      if(r[n] < 0) sign[n] = -1;
      else sign[n] = 1;
    }
    
    for (n in 1:N) {
      acc += log1p((fabs(r[n]) / (vs[n]*(l*sign[n] + 1)))^p / q);
    }
    
    out -= (1.0/p+q) * acc;
    
    return out;
  }
}
data {
  int<lower=0> T;
  vector[T] y;
  
  real m_p; real<lower=0> s_p;
  real m_q; real<lower=0> s_q;
}
parameters {
  
  // Scale parameters
  real omega;
  real phi;
  // real iota;
  real k1;
  real k2;
  // real k3;
  // real k4;

  // Scale of very first observation
  real lambda1;
  
  // Other SGT Parameters
  real mu;
  real<lower=-1, upper = 1> l;
  real<lower=1> p;
  real<lower=2> q;
}
transformed parameters {
  vector[T] h;
  vector[T] lambda;
  
  h[1] = 0;
  lambda[1] = lambda1;
  
  for (t in 2:T) {
    // h[t] = k1*y[t-1] + k2*y[t-1]^2 + k3*y[t-1]^3 + k4*y[t-1]^4;
    h[t] = k1*y[t-1] + k2*y[t-1]^2;
    lambda[t] = omega + phi*lambda[t-1] + h[t];
  }
}
model {
  // prior
  // k3 ~ normal(0, 1e-1);
  p ~ normal(m_p, s_p);
  q ~ normal(m_q, s_q);
  
  // likelihood
  y ~ sgt(mu, exp(lambda), l, p, q);
}
