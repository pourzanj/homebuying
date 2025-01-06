data {
  int<lower=0> T;   // # time points (equally spaced)
  matrix[T,3] y;
}
parameters {
  real mu_ret;
  real<lower=0> mu_hloffset;

  real mu;                     // mean log volatility
  vector<lower=0, upper=1>[2] phi_raw; 
  positive_ordered[2] sigma;         // white noise shock scale
  
  real<lower=0> sigma_obs[2];
  
  cholesky_factor_corr[5] L;
  
  matrix[T,2] v_raw;
}
transformed parameters {
  corr_matrix[5] R = multiply_lower_tri_self_transpose(L);
  cov_matrix[3] R11_R12R21 = R[1:3,1:3] - R[1:3,4:5]*R[4:5,1:3];

  vector[2] phi;
  matrix[T,2] v; 
  matrix[T,2] h;
  vector[T] H;
  matrix[T,3] s;
  matrix[T,3] m;
  
  // Set phi
  phi[1] = -1 + 2 * phi_raw[1];
  for(k in 2:2) {
    phi[k] = -1 + (phi[k-1] - -1) * phi_raw[k];
  }
  
  // Set v, h
  v[1,] = sigma' .* v_raw[1,];
  h[1,] = v[1,] ./ sqrt(1 - phi^2)';
  H[1] = sum(h[1,]);
  
  s[1,] = [exp((mu+H[1])/2), sigma_obs[1], sigma_obs[2]];
  m[1,] = [mu_ret, mu_hloffset + (mu+H[1])/2, (100/2)*sum(v[1,])];
  
  for (t in 2:T) {
    v[t,] = sigma' .* v_raw[t,];
    h[t,] = phi' .* h[t-1,] + v[t,];
    H[t] = sum(h[t,]);
    
    s[t,] = [exp((mu+H[t])/2), sigma_obs[1], sigma_obs[2]];
    m[t,] = [mu_ret,
             mu_hloffset + (mu+H[t])/2,
             (100/2)*(H[t]-H[t-1])] +
                (R[1:3,4:5]*v_raw[t,]')';
  }
  
}
model {
  // prior
  phi ~ uniform(-1, 1);
  
  sigma[1] ~ exponential(10.0);
  sigma[2] ~ exponential(1.0);
  
  mu ~ normal(0, 1);

  mu_ret ~ normal(0.1, 0.05);
  mu_hloffset ~ normal(0.25, 0.025);

  sigma_obs[1] ~ exponential(1.0);
  sigma_obs[2] ~ exponential(1e-2);

  L ~ lkj_corr_cholesky(10.0);
  // (R[1,4]+1)/2.0 ~ beta(3, 10);

  // likelihood
  to_vector(v_raw) ~ std_normal();
  for (t in 1:T) {
    y[t,] ~ multi_normal(m[t,], quad_form_diag(R11_R12R21, s[t,]'));
  }
}
generated quantities {
  // vector[T] v_raw_rep = to_vector(normal_rng(rep_vector(0.0, T), 1.0));
  // vector[T] v_rep; 
  // vector[T] h_rep;
  // matrix[T,3] s_rep;
  // matrix[T,3] m_rep;
  // matrix[T,3] y_rep;
  // 
  // vector[T] log_lik;
  // vector[T] log_lik_rep;
  // 
  // v_rep[1] = sigma*v_raw_rep[1];
  // h_rep[1] = v_rep[1] / sqrt(1 - phi^2);
  // s_rep[1,] = [exp((mu+h_rep[1])/2), sigma_obs[1], sigma_obs[2]];
  // m_rep[1,] = [mu_ret, mu_hloffset + (mu+h_rep[1])/2, (100/2)*((phi-1)*h_rep[1]+v_rep[1])];
  // y_rep[1,] = multi_normal_rng(m_rep[1,], quad_form_diag(R11_R12R21, s_rep[1,]'))';
  // for (t in 2:T) {
  //   v_rep[t] = sigma*exp(a*h[t-1]/2+b*v_raw_rep[t] + c*v_raw_rep[t-1])*v_raw_rep[t];
  //   h_rep[t] = phi * h[t-1] + v_rep[t];
  //   s_rep[t,] = [exp((mu+h_rep[t])/2), sigma_obs[1], sigma_obs[2]];
  //   m_rep[t,] = [mu_ret,
  //                mu_hloffset + (mu+h_rep[t])/2,
  //                (100/2)*((phi-1)*h[t-1]+v_rep[t])] +
  //                   R[1:3,4]'*exp(a*h[t-1]/2 + b*v_raw_rep[t] + c*v_raw[t-1])*v_raw_rep[t];
  //                   
  //   // y_rep[t,] = multi_normal_rng(m_rep[t,], quad_form_diag(R11_R12R21, s_rep[t,]'))';
  //   // 
  //   // log_lik[t] = multi_normal_lpdf(y[t,] | m[t,], quad_form_diag(R11_R12R21, s[t,]'));
  //   // log_lik_rep[t] = multi_normal_lpdf(y[t,] | m_rep[t,], quad_form_diag(R11_R12R21, s_rep[t,]'));
  //   
  //   y_rep[t,] = multi_normal_cholesky_rng(m_rep[t,], diag_pre_multiply(s_rep[t,], L11_L12L21))';
  //   
  //   log_lik[t] = multi_normal_cholesky_lpdf(y[t,] | m[t,], diag_pre_multiply(s[t,], L11_L12L21));
  //   log_lik_rep[t] = multi_normal_cholesky_lpdf(y[t,] | m_rep[t,], diag_pre_multiply(s_rep[t,], L11_L12L21));
  // }
}



