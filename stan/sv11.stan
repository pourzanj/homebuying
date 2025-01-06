functions {
  real spline1(real x, vector b, vector t) {
    int N = size(t);
    vector[N] f;
    
    real Bi;
    real Bi1;
    real omega_ik;
    real omega_ik1;
    
    Bi = 0;
    Bi1 = x < t[2];
    omega_ik = 0;
    omega_ik1 = (x-t[1]) / (t[2]-t[1]);
    f[1] = omega_ik*Bi + (1-omega_ik1)*Bi1;
    
    for(i in 2:(N-1)) {
      Bi = t[i-1] <= x && x < t[i];
      Bi1 = t[i] <= x && x < t[i+1];
      omega_ik = (x-t[i-1]) / (t[i]-t[i-1]);
      omega_ik1 = (x-t[i]) / (t[i+1]-t[i]);
      f[i] = omega_ik*Bi + (1-omega_ik1)*Bi1;
    }
    
    Bi = t[N-1] <= x;
    Bi1 = 0;
    omega_ik = (x-t[N-1]) / (t[N]-t[N-1]);
    omega_ik1 = 0;
    f[N] = omega_ik*Bi + (1-omega_ik1)*Bi1;
    
    return dot_product(b, f);
  }
}

data {
  int<lower=0> T;   // # time points (equally spaced)
  vector[T] y;      // mean corrected return at time t
  
  int K;
  vector[K] knots;
}

parameters {
  real mu_ret;
  real<lower=-1, upper=1> rho;

  // SV
  real mu;                     // mean log volatility
  real<lower=-1, upper=1> phi;  // persistence of volatility
  real<lower=0> sigma;         // white noise shock scale
  
  vector[T] v_raw;             // std log volatility time t
  
  vector[K] b;
}
transformed parameters {
  vector[T] mul;
  vector[T] v;
  vector[T] h;
  vector[T] s;

  mul[1] = 1;
  v[1] = v_raw[1];
  h[1] = sigma*v[1] / sqrt(1 - phi^2);
  
  for (t in 2:T) {
    mul[t] = exp(spline1(h[t-1], b, knots));
    v[t] = mul[t]*v_raw[t];
    h[t] = phi*h[t-1] + sigma*v[t];
  }
  
  s = exp((mu + h) / 2.0);
}
model {
  // prior
  // mu_ret ~ normal(0, 0.2);
  // rho ~ normal(-1, 1);
  // rho ~ normal(-0.75, 0.05);

  // (phi+1)/2.0 ~ beta(8, 2);
  // sigma ~ exponential(10);
  sigma ~ normal(0, 0.2);
  mu ~ normal(0, 0.5);
  
  b ~ normal(0, 0.5);
  
  // likelihood
  v_raw ~ std_normal();
  y ~ normal(mu_ret + rho*v .* s, sqrt(1-rho^2)*s);
  
}
generated quantities {
  vector[T] v_raw_rep = to_vector(normal_rng(rep_array(0.0, T), 1));
  vector[T] v_rep;
  vector[T] h_rep;
  vector[T] s_rep;
  vector[T] y_rep;

  vector[T] log_lik;
  real sum_log_lik;
  vector[T] log_lik_ahead;
  real sum_log_lik_ahead;
  
  real m1_v_raw = mean(v_raw);
  real m2_v_raw = mean((v_raw-m1_v_raw)^2);
  real m3_v_raw = mean((v_raw-m1_v_raw)^3);
  real m4_v_raw = mean((v_raw-m1_v_raw)^4);
  
  real m1_v_raw_rep = mean(v_raw_rep);
  real m2_v_raw_rep = mean((v_raw_rep-m1_v_raw_rep)^2);
  real m3_v_raw_rep = mean((v_raw_rep-m1_v_raw_rep)^3);
  real m4_v_raw_rep = mean((v_raw_rep-m1_v_raw_rep)^4);

  v_rep[1] = v_raw_rep[1];
  h_rep[1] = sigma*v_rep[1] / sqrt(1 - phi^2);

  for(t in 2:T) {
    v_rep[t] = exp(spline1(h[t-1], b, knots))*v_raw_rep[t];
    h_rep[t] = phi*h[t-1] + sigma*v_rep[t];
    s_rep[t] = exp((mu+h_rep[t])/2.0);
    y_rep[t] = normal_rng(mu_ret + rho*v_rep[t]*s_rep[t], sqrt(1-rho^2)*s_rep[t]);
    log_lik[t] = normal_lpdf(y[t] | mu_ret + rho*v[t]*s[t], sqrt(1-rho^2)*s[t]);
    log_lik_ahead[t] = normal_lpdf(y[t] | mu_ret + rho*v_rep[t]*s_rep[t], sqrt(1-rho^2)*s_rep[t]);
  }
  sum_log_lik = sum(log_lik[2:T]);
  sum_log_lik_ahead = sum(log_lik_ahead[2:T]);
  
}
