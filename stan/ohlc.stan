functions {
  real sgt_lpdf(vector x, real mu, real s, real l, real p, real q)
  { // Skewed generalised t
    int N = dims(x)[1];
    real lz1;
    real lz2;
    real v;
    real m;
    vector[N] r;
    real out;
    real acc;
    real vs;
    int sign[N];
    lz1 = lbeta(1.0/p,q);
    lz2 = lbeta(2.0/p,q-1.0/p);
    //v = q^(-1.0/p)*((3*l^2+1)*exp(lbeta(3.0/p,q-2.0/p)-lz1)-4*l^2*exp(lz2-lz1)^2)^(-0.5);
    v = 1;
    vs = v*s;

    out = 0;
    out += log(p);
    out -= log(2*v);
    out -= log(s);
    out -= (1 / p) * log(q);
    out -= lz1;
    out = N*out;
    
    m = 0;
    r = x-mu;
    
    acc = 0.0;
    
    for (n in 1:N) {
      if(r[n] < 0) sign[n] = -1;
      else sign[n] = 1;
    }
    
    for (n in 1:N) {
      acc += log1p((fabs(r[n]) / (vs*(l*sign[n] + 1)))^p / q);
    }
    
    // acc += sum(log1p((fabs(r) ./ (vs*(l*sign + 1)))^p / q));
    
    out -= (1.0/p+q) * acc;
    
    return out;
  }
}
data {
  int<lower=0> T;
  
  vector[T] r_co;
  vector[T] r_oc;
  vector[T] hl;
}
parameters {
  // Likelihood
  real m_co;
  real m_oc;
  real m_hl;
  
  real rho_co;
  real rho_oc;
  real rho_hl;
  
  real<lower=0> a;
  real<lower=0> b;
  
  // real<lower=-1, upper = 1> l;
  // real<lower=0.5> p;
  // real<lower=2> q;

  // SV
  // real mu;                      // mean log volatility
  vector[T] mu_b;
  vector[T] mu_l;
  real<lower=-1, upper=1> phi;  // persistence of volatility
  real<lower=0> sigma;          // white noise shock scale

  vector[T] v;
  // vector<lower=0>[T] v_raw2;
  // vector<lower=0>[T] l;
}
transformed parameters {
  // vector[T] v = v_raw1;
  // vector[T] v2 = v_raw2 .* l;
  vector[T] h;
  vector[T] s;
  
  vector[T] mu;
  
  h[1] = (sigma * v[1]) / sqrt(1 - phi^2);

  for (t in 2:T) {
    h[t] = phi * h[t-1] + sigma * v[t];
  }
  
  mu[1] = -3.0;
  for (t in 2:T) {
    mu[t] = phi * h[t-1] + sigma * v[t];
  }
  
  s = exp((mu + h) / 2);
}
model {
  // v_raw1 ~ normal(0, 1);
  // v_raw2 ~ normal(0, 1);
  // l ~ cauchy(0, 1);
  // p ~ double_exponential(0.5, 1);
  // q ~ normal(20, 5);
  // v ~ sgt(0, 1, 0, 2, 100);
  v ~ normal(0, 1);

  r_co ~ double_exponential(m_co + rho_co * v.* s, s);
  r_oc ~ double_exponential(m_oc + rho_oc * v.* s, a * s);
  hl ~ normal(m_hl + rho_hl * v.* s, b * s);
  // r_co ~ double_exponential(m_co + rho_co * v.* s, s * sqrt(1-rho_co^2));
  // r_oc ~ double_exponential(m_oc + rho_oc * v.* s, a * s * sqrt(1-rho_oc^2));
  // hl ~ normal(m_hl + rho_hl * v.* s, b * s * sqrt(1-rho_hl^2));
}
generated quantities {
  vector[T] v_rep = to_vector(normal_rng(rep_vector(0.0, T), rep_vector(1.0, T)));
  vector[T] h_rep;
  vector[T] s_rep;
  
  vector[T] rco_rep;
  vector[T] roc_rep; 
  vector[T] r_rep;
  
  h_rep[1] = h[1];
  for (t in 2:T) {
    h_rep[t] = phi * h[t-1] + sigma * v_rep[t];
  }
  s_rep = exp((mu + h_rep) / 2);
  
  rco_rep = to_vector(double_exponential_rng(m_co + rho_co * v_rep .* s_rep, s_rep));
  roc_rep = to_vector(double_exponential_rng(m_oc + rho_oc * v_rep .* s_rep, a * s_rep));
  r_rep = 100*(exp(rco_rep/100.0) .* exp(roc_rep/100.0)-1.0);
}
