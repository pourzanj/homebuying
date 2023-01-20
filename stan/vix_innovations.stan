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
    v = q^(-1.0/p)*((3*l^2+1)*exp(lbeta(3.0/p,q-2.0/p)-lz1)-4*l^2*exp(lz2-lz1)^2)^(-0.5);
    // v = 1;
    vs = v*s;

    out = 0;
    out += log(p);
    out -= log(2*v);
    out -= log(s);
    out -= (1 / p) * log(q);
    out -= lz1;
    out = N*out;

    // m = 0;
    m = 2*vs*l*q^(1.0/p) * exp(lz2-lz1);
    r = x-mu+m;

    acc = 0.0;

    for (n in 1:N) {
      if(r[n] < 0) sign[n] = -1;
      else sign[n] = 1;
    }

    for (n in 1:N) {
      acc += log1p((abs(r[n]) / (vs*(l*sign[n] + 1)))^p / q);
    }

    // acc += sum(log1p((fabs(r) ./ (vs*(l*sign + 1)))^p / q));

    out -= (1.0/p+q) * acc;

    return out;
  }
}

data {
  int<lower=0> T;   // # time points (equally spaced)
  vector[T] v;
  vector[T] h;
  
  real m_p; real<lower=0> s_p;
  real m_q; real<lower=0> s_q;
}

parameters {
  real beta[2];
  real sigma[2];
  real lambda[2];
  real<lower=1> p;
  real<lower=2> q;
  // real q[2];
}

model {
  p ~ normal(m_p, s_p);
  q ~ normal(m_q, s_q);
  
  // q[1] ~ normal(0, 1);
  // q[2] ~ normal(0, 1);
  
  for(t in 2:T) {
    // (exp(q[1] + q[2]*h[t-1]) + 2) ~ normal(0, 20);
    // to_vector([v[t]]) ~ sgt(beta[1] + beta[2]*h[t-1],
    //                         exp(sigma[1] + sigma[2]*h[t-1]),
    //                         2*inv_logit(lambda[1] + lambda[2]*h[t-1]) - 1,
    //                         exp(p[1] + p[2]*h[t-1]) + 1,
    //                         exp(q[1] + q[2]*h[t-1]) + 2);
        to_vector([v[t]]) ~ sgt(beta[1] + beta[2]*h[t-1],
                            exp(sigma[1] + sigma[2]*h[t-1]),
                            2*inv_logit(lambda[1] + lambda[2]*h[t-1]) - 1,
                            p,
                            q);
  }
}
