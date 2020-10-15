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
    //v = q^(-1.0/p)*((3*l^2+1)*exp(lbeta(3.0/p,q-2.0/p)-lz1)-4*l^2*exp(lz2-lz1)^2)^(-0.5);
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
    
    // acc += sum(log1p((fabs(r) ./ (vs*(l*sign + 1)))^p / q));
    
    out -= (1.0/p+q) * acc;
    
    return out;
  }
}
data {
  int<lower=0> T;
  vector[T] r;
  
  real m_alpha0; real<lower=0> s_alpha0;
  real m_p; real<lower=0> s_p;
  real m_q; real<lower=0> s_q;
}
parameters {
  real mu;
  
  real<lower=0> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=(1-alpha1)> beta1;
  real<lower=0> sigma1;
  
  real<lower=-1, upper = 1> l;
  real<lower=1> p;
  real<lower=0.5> q;
}
transformed parameters {
  vector<lower=0>[T] sigma;
  sigma[1] = sigma1;
  for (t in 2:T)
    sigma[t] = sqrt(alpha0
                     + alpha1 * pow(r[t-1] - mu, 2)
                     + beta1 * pow(sigma[t-1], 2));
}
model {
  // prior
  alpha0 ~ normal(m_alpha0, s_alpha0);
  p ~ normal(m_p, s_p);
  q ~ normal(m_q, s_q);
  
  // likelihood
  // r ~ sgt(mu, sigma, l, p, q);
  r ~ sgt(mu, sigma, l, p, q);
}
