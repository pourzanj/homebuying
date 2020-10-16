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
  
  real sign(real x) {
    if(x < 0.0) {
      return -1.0;
    } else {
      return 1.0;
    }
  }
  
  real u(real y, real mu, real sigma, real l, real p, real q) {
    // Get star
    real num = fabs(y - mu)^p;
    real den = q*sigma^p * (l*sign(y-mu) + 1)^p;
    real star = num / den;
    
    real d_star;
    real df_dsigma;
  
    print("y: ", y);
    print("mu: ", mu);
    print("p: ", p);
    
    print("num: ", num);
    print("den: ", den);
    print("star: ", star);
  
    // Get d_star
    num = p * fabs(y - mu)^p;
    den = q*sigma^(p+1) * (l*sign(y-mu) + 1)^p;
    d_star = -num / den;
    
    print("d_star: ", d_star);
  
    // Return final
    df_dsigma = -(1/sigma) - (1/p + q) * (1 / (star+1)) * d_star;
    print("df_dsigma: ", df_dsigma);
    return(df_dsigma * sigma);
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
  real k1;
  real k2;

  // Scale of very first observation
  real lambda1;
  
  // Other SGT Parameters
  real mu;
  real<lower=-1, upper = 1> l;
  real<lower=1> p;
  real<lower=2> q;
}
transformed parameters {
  vector[T] u_t;
  vector[T] h;
  vector[T] lambda;
  
  u_t[1] = 0;
  h[1] = 0;
  lambda[1] = lambda1;
  
  for (t in 2:T) {
    u_t[t] = u(y[t-1], mu, exp(lambda[t-1]), l, p, q);
    h[t] = k1*u_t[t] + k2*sign(-y[t-1])*(u_t[t]+1);
    lambda[t] = omega + phi*lambda[t-1] + h[t];
  }
}
model {
  // prior
  p ~ normal(m_p, s_p);
  q ~ normal(m_q, s_q);
  
  // likelihood
  y ~ sgt(mu, exp(lambda), l, p, q);
}
