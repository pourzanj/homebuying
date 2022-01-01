data {
  int<lower=0> T;   // # time points (equally spaced)
  matrix[T,2] y;
}
parameters {
  real mu_ret;
  real<lower=0> mu_hloffset;

  real mu;                     // mean log volatility
  real<lower=-1, upper=1> phi;  // persistence of volatility
  real<lower=0> sigma;         // white noise shock scale
  
  cholesky_factor_corr[3] L;
  
  vector[T] v_raw;
}
transformed parameters {
  corr_matrix[3] R = multiply_lower_tri_self_transpose(L);
  cov_matrix[2] R11_R12R21 = R[1:2,1:2] - R[1:2,3]*R[3,1:2];
  vector[T] v = sigma * v_raw; 
  vector[T] h;
  matrix[T,2] s;
  matrix[T,2] m;
  
  h[1] = v[1] / sqrt(1 - phi^2);
  s[1,] = [exp((mu+h[1])/2), 1.0];
  m[1,] = [mu_ret, mu_hloffset + (mu+h[1])/2];
  for (t in 2:T) {
    h[t] = phi * h[t-1] + v[t];
    s[t,] = [exp((mu+h[t])/2), 1.0];
    m[t,] = [mu_ret, mu_hloffset + (mu+h[t])/2] + R[1:2,3]'*v_raw[t];
  }
  
}
model {
  // prior
  phi ~ uniform(-1, 1);
  sigma ~ exponential(1.0);
  mu ~ normal(0, 1);
  
  mu_ret ~ normal(0.1, 0.05);
  mu_hloffset ~ normal(0.25, 0.025);

  L ~ lkj_corr_cholesky(10.0);

  // likelihood
  v_raw ~ std_normal();
  for (t in 1:T) {
    y[t,] ~ multi_normal(m[t,], quad_form_diag(R11_R12R21, s[t,]'));
    // y[t,] ~ multi_normal([0,0], [[1,0],[0,1]]);
  }
}
