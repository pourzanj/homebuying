functions {
  // real sgt_lpdf(real x, real mu, real s, real l, real p, real q)
  // { // Skewed generalised t
  //   // int N = dims(x)[1];
  //   real lz1;
  //   real lz2;
  //   real v;
  //   real m;
  //   real r;
  //   real out;
  //   real acc;
  //   real vs;
  //   int sign;
  //   lz1 = lbeta(1.0/p,q);
  //   lz2 = lbeta(2.0/p,q-1.0/p);
  //   // v = q^(-1.0/p)*((3*l^2+1)*exp(lbeta(3.0/p,q-2.0/p)-lz1)-4*l^2*exp(lz2-lz1)^2)^(-0.5);
  //   // v = 1;
  //   vs = v*s;
  // 
  //   out = 0;
  //   out += log(p);
  //   out -= log(2*v);
  //   out -= log(s);
  //   out -= (1 / p) * log(q);
  //   out -= lz1;
  //   // out = N*out;
  //   out = out;
  //   
  //   m = 0;
  //   r = x-mu;
  //   
  //   acc = 0.0;
  //   
  //   // for (n in 1:N) {
  //   //   if(r[n] < 0) sign[n] = -1;
  //   //   else sign[n] = 1;
  //   // }
  //   
  //   if(r < 0) sign = -1;
  //   else sign = 1;
  //   
  //   // for (n in 1:N) {
  //   //   acc += log1p((fabs(r[n]) / (vs*(l*sign[n] + 1)))^p / q);
  //   // }
  //   
  //   acc += log1p((fabs(r) / (vs*(l*sign + 1)))^p / q);
  //   
  //   // acc += sum(log1p((fabs(r) ./ (vs*(l*sign + 1)))^p / q));
  //   
  //   out -= (1.0/p+q) * acc;
  //   
  //   return out;
  // }
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
      acc += log1p((fabs(r[n]) / (vs*(l*sign[n] + 1)))^p / q);
    }

    // acc += sum(log1p((fabs(r) ./ (vs*(l*sign + 1)))^p / q));

    out -= (1.0/p+q) * acc;

    return out;
  }
}
data {
  int<lower=0> T;   // # time points (equally spaced)
  vector[T] y;      // mean corrected return at time t
  
  real m_p; real<lower=0> s_p;
  real m_q; real<lower=0> s_q;
}

parameters {
  real mu_ret;
  real alpha[2];
  // real<lower=-1, upper=1> l;
  real beta[3];
  // real<lower=1> p;
  // real<lower=2> q;
  real nu[3];
  real gamma[3];

  // SV
  real mu;                     // mean log volatility
  real<lower=-1, upper=1> phi;  // persistence of volatility
  real<lower=0> sigma;         // white noise shock scale
  
  vector[T] v;             // std log volatility time t
}
transformed parameters {
  vector[T] h;
  vector[T] p;
  vector[T] l;
  vector[T] m;
  vector[T] s;
  vector[T] q;

  m[1] = mu_ret;
  p[1] = exp(beta[1]);
  l[1] = 2*inv_logit(gamma[1]) - 1;
  h[1] = sigma*v[1] / sqrt(1 - phi^2);
  q[1] = exp(nu[1]);
  
  for (t in 2:T) {
    h[t] = phi*h[t-1] + sigma*v[t];
    m[t] = mu_ret + alpha[1]*h[t-1] + alpha[2]*h[t-1]^2;
    l[t] = 2*inv_logit(gamma[1] + gamma[2]*h[t-1] + gamma[3]*h[t-1]^2) - 1;
    p[t] = exp(beta[1] + beta[2]*h[t-1] + beta[3]*h[t-1]^2);
    q[t] = exp(nu[1] + nu[2]*h[t-1] + nu[3]*h[t-1]^2);
  }
  
  s = exp(0.5*(mu + h));
  // s = exp(mu + gamma[1]*h + gamma[2]*h^3);
}
model {
  // prior
  // mu_ret ~ normal(0, 0.2);
  // rho ~ normal(-1, 1);

  // (phi+1)/2.0 ~ beta(8, 2);
  // sigma ~ exponential(10);
  // mu ~ normal(0, 1);
  
  // gamma ~ normal(2, 2);
  // gamma[1] ~ normal(0, 1);
  // gamma[2] ~ normal(0, 0.1);
  
  p ~ normal(m_p, s_p);
  q ~ normal(m_q, s_q);
  
  // likelihood
  v ~ std_normal();
  for(t in 1:T) {
    to_vector([y[t]]) ~ sgt(m[t], s[t], l[t], p[t], q[t]);
  }
  
  
}
generated quantities {
  vector[T] v_rep = to_vector(normal_rng(rep_array(0.0, T), 1));
  vector[T] h_rep;
  vector[T] h_ahead_rep;
  // vector[T] s_rep;
  // vector[T] y_rep;
  // 
  // vector[T] log_lik;
  // real sum_log_lik;
  // vector[T] log_lik_ahead;
  // real sum_log_lik_ahead;
  // 
  h_rep[1] = sigma*v_rep[1] / sqrt(1 - phi^2);
  // 
  for(t in 2:T) {
    h_rep[t] = phi*h_rep[t-1] + sigma*v_rep[t];
    h_ahead_rep[t] = phi*h[t-1] + sigma*v_rep[t];
    // y_rep[t] = normal_rng(mu_ret, s_rep[t]);
    // log_lik[t] = normal_lpdf(y[t] | mu_ret, s[t]);
    // log_lik_ahead[t] = normal_lpdf(y[t] | mu_ret, s_rep[t]);
  }
  // sum_log_lik = sum(log_lik[2:T]);
  // sum_log_lik_ahead = sum(log_lik_ahead[2:T]);
}
