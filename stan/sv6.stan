data {
  int<lower=0> T;   // # time points (equally spaced)
  matrix[T,2] y;
}

transformed data {
  // real phi = 1;
}

parameters {
  real mu_ret;
  
  real mu_hloffset;
  // real<lower=0> tau;

  // SV
  real mu;                     // mean log volatility
  real<lower=-1, upper=1> phi;  // persistence of volatility
  real<lower=0> sigma;         // white noise shock scale
  
  cholesky_factor_corr[3] L;
  
  vector[T] v;             // std log volatility time t
}
transformed parameters {
  corr_matrix[3] R = multiply_lower_tri_self_transpose(L);
  cov_matrix[2] R11_R12R21 = R[1:2,1:2] - R[1:2,3]*R[3,1:2];
  
  vector[T] h;
  
  matrix[T,2] s;
  matrix[T,2] m;
  
  // h[1] = sigma*v[1] / sqrt(1 - phi^2);
  h[1] = sigma*v[1];
  s[1,] = [exp((mu+h[1])/2), 0.38];
  m[1,] = [mu_ret, mu_hloffset + (mu+h[1])/2];
  
  for (t in 2:T) {
    h[t] = phi * h[t-1] + sigma*v[t];
    
    s[t,] = [exp((mu+h[t])/2), 0.38];
    m[t,] = [mu_ret, mu_hloffset + (mu+h[t])/2] + R[1:2,3]'*v[t];
  }
  
}
model {
  // prior
  // mu_ret ~ normal(0, 0.2);
  // rho ~ normal(-1, 1);

  // (phi+1)/2.0 ~ beta(8, 2);
  // sigma ~ exponential(10);
  // mu ~ normal(0, 1);
  
  // likelihood
  v ~ std_normal();
  for (t in 1:T) {
    y[t,] ~ multi_normal(m[t,], quad_form_diag(R11_R12R21, s[t,]'));
  }
  
}
generated quantities {
  vector[T] v_rep = to_vector(normal_rng(rep_array(0.0, T), 1));
  vector[T] h_rep;

  matrix[T,2] s_rep;
  matrix[T,2] m_rep;

  matrix[T,2] y_rep;

  // h_rep[1] = sigma*v_rep[1] / sqrt(1 - phi^2);
  h_rep[1] = sigma*v_rep[1];
  s_rep[1,] = [exp((mu+h_rep[1])/2), 0.38];
  m_rep[1,] = [mu_ret, mu_hloffset + (mu+h[1])/2];

  for(t in 2:T) {
    h_rep[t] = phi*h[t-1] + sigma*v_rep[t];

    s_rep[t,] = [exp((mu+h_rep[t])/2), 0.38];
    m_rep[t,] = [mu_ret, mu_hloffset + (mu+h_rep[t])/2] + R[1:2,3]'*v_rep[t];
    
    y_rep[t,] = multi_normal_rng(m_rep[t,], quad_form_diag(R11_R12R21, s_rep[t,]'))';
  }

}
