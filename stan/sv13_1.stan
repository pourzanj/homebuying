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
  
  int<lower=2> M;
  vector[M] y_oos;
  
  // Priors on p and q since they're hard to identify
  real m_p; real<lower=0> s_p;
  real m_q; real<lower=0> s_q;
}

parameters {

  // SV
  real mu;                     // mean log volatility
  real<lower=-1, upper=1> phi;  // persistence of volatility
  real<lower=0> sigma;         // white noise shock scale
  
  // Return Shape
  real alpha[3];
  real gamma[3];
  real beta[2];
  real<lower=2> q;
  
  // Volatility innovations
  real v0;
  vector[T] v;             // std log volatility time t
  
  vector[M] v_oos;
}
transformed parameters {
  real h0;
  vector[T] h;
  
  vector[T] m;
  vector[T] s;
  vector[T] l;
  vector[T] p;
  
  h0 = sigma*v0 / sqrt(1 - phi^2);
  h[1] = phi*h0 + sigma*v[1];
  m[1] = alpha[1] + alpha[2]*h0 + alpha[3]*h0^2;
  s[1] = exp(0.5*(mu + h[1]));
  l[1] = 2*inv_logit(gamma[1] + gamma[2]*h0 + gamma[3]*h0^2) - 1;
  p[1] = exp(beta[1] + beta[2]*h0);

  for (t in 2:T) {
    h[t] = phi*h[t-1] + sigma*v[t];
    m[t] = alpha[1] + alpha[2]*h[t-1] + alpha[3]*h[t-1]^2;
    s[t] = exp(0.5*(mu + h[t]));
    l[t] = 2*inv_logit(gamma[1] + gamma[2]*h[t-1] + gamma[3]*h[t-1]^2) - 1;
    p[t] = exp(beta[1] + beta[2]*h[t-1]);
  }
}
model {
  // prior
  mu ~ normal(0, 2);
  sigma ~ normal(0, 1);
  
  alpha ~ normal(0, 0.2);
  gamma ~ normal(0, 1);
  beta ~ normal(0, 2);
  
  p ~ normal(m_p, s_p);
  q ~ normal(m_q, s_q);
  
  // likelihood
  v0 ~ std_normal();
  v ~ std_normal();
  v_oos ~ std_normal();
  for(t in 1:T) {
    to_vector([y[t]]) ~ sgt(m[t], s[t], l[t], p[t], q);
  }
}
generated quantities {
  vector[M] h_oos;
  vector[M] m_oos;
  vector[M] s_oos;
  vector[M] l_oos;
  vector[M] p_oos;
  vector[M] loglik_oos;

  h_oos[1] = phi*h[T] + sigma*v_oos[1];
  m_oos[1] = alpha[1] + alpha[2]*h[T] + alpha[3]*h[T]^2;
  s_oos[1] = exp(0.5*(mu + h_oos[1]));
  l_oos[1] = 2*inv_logit(gamma[1] + gamma[2]*h[T] + gamma[3]*h[T]^2) - 1;
  p_oos[1] = exp(beta[1] + beta[2]*h[T]);
  loglik_oos[1] = sgt_lpdf(to_vector([y_oos[1]]) | m_oos[1], s_oos[1], l_oos[1], p_oos[1], q);
  
  for (t in 2:M) {
    h_oos[t] = phi*h_oos[t-1] + sigma*v_oos[t];
    m_oos[t] = alpha[1] + alpha[2]*h_oos[t-1] + alpha[3]*h_oos[t-1]^2;
    s_oos[t] = exp(0.5*(mu + h_oos[t]));
    l_oos[t] = 2*inv_logit(gamma[1] + gamma[2]*h_oos[t-1] + gamma[3]*h_oos[t-1]^2) - 1;
    p_oos[t] = exp(beta[1] + beta[2]*h_oos[t-1]);
    loglik_oos[t] = sgt_lpdf(to_vector([y_oos[t]]) | m_oos[t], s_oos[t], l_oos[t], p_oos[t], q);
  }
}
