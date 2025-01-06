data {
  int<lower=0> T;   // # time points (equally spaced)
  matrix[T,3] y;
  
  real<lower=0> p0;
  real<lower=0> scale;
  real<lower=0> nu;
  // real<lower=0> C2;
  // real<lower=0> p;
}
transformed data {
  // vector[T] v_raw_base = rep_vector(0.0, T);
  vector<lower=0>[T] l = rep_vector(1e-5, T);
}
parameters {
  real mu_ret;
  real<lower=0> mu_hloffset;

  real<lower=0> mu;                     // mean log volatility
  real<lower=-1, upper=1> phi;  // persistence of volatility
  real<lower=0> sigma;         // white noise shock scale
  
  // real<lower=0> p;
  // real<lower=0> C2;
  // real<lower=0> theta;
  
  real<lower=0> sigma_obs[2];
  
  cholesky_factor_corr[4] L;
  
  vector<lower=0>[T] v_raw;
  // vector<lower=0>[T] l;
}
transformed parameters {
  corr_matrix[4] R = multiply_lower_tri_self_transpose(L);
  cov_matrix[3] R11_R12R21 = R[1:3,1:3] - R[1:3,4]*R[4,1:3];

  // vector<lower=0>[T] pT;
  // vector[T] v_raw_tail;
  // vector[T] v_raw;
  vector[T] v; 
  vector[T] h;
  matrix[T,3] s;
  matrix[T,3] m;
  
  // pT[1] = p;
  // v_raw_tail[1] = pT[1] * sqrt(C2*l[1]^2 ./ (C2 + pT[1]^2 * l[1]^2));
  v[1] = sigma*v_raw[1];
  h[1] = mu + v[1] / sqrt(1 - phi^2);
  
  s[1,] = [h[1], sigma_obs[1], sigma_obs[2]];
  m[1,] = [mu_ret, mu_hloffset + h[1], 100*(log(v[1]))];
  
  for (t in 2:T) {
    // pT[t] = (1-theta)*pT[t-1]+theta*v_raw_tail[t-1];
    // v_raw_tail[t] = pT[t] * sqrt(C2*l[t]^2 ./ (C2 + pT[t]^2 * l[t]^2));
    v[t] = sigma*v_raw[t];
    h[t] = mu + phi * h[t-1] + v[t];
    
    s[t,] = [h[t], sigma_obs[1], sigma_obs[2]];
    m[t,] = [mu_ret,
             mu_hloffset + h[t],
             100*((phi-1)*h[t-1]+v[t])] +
                R[1:3,4]'*v_raw[t];
  }
  
}
model {
  // prior
  phi ~ uniform(-1, 1);
  sigma ~ exponential(1.0);
  mu ~ normal(-5, 2.5);

  mu_ret ~ normal(0.1, 0.05);
  mu_hloffset ~ normal(0.25, 0.025);

  // p ~ cauchy(0, p0);
  // C2 ~ inv_gamma(0.5*nu, 0.5*nu*scale^2);

  sigma_obs[1] ~ exponential(1.0);
  sigma_obs[2] ~ exponential(1e-2);

  L ~ lkj_corr_cholesky(10.0);
  // (R[1,4]+1)/2.0 ~ beta(3, 10);

  // likelihood
  v_raw ~ std_normal();
  // l ~ cauchy(0,1);
  
  for (t in 1:T) {
    y[t,] ~ multi_normal(m[t,], quad_form_diag(R11_R12R21, s[t,]'));

  }
}
generated quantities {
  // vector[T] v_raw_base_rep = to_vector(normal_rng(rep_vector(0.0, T), 1.0));
  // // vector[T] v_raw_base_rep = rep_vector(0.0, T);
  // // vector[T] l_rep = to_vector(cauchy_rng(rep_vector(0.0, T), 1.0));
  // vector[T] l_rep = rep_vector(1e-5, T);
  // // vector<lower=0>[T] pT_rep;
  // 
  // vector[T] v_raw_tail_rep;
  // vector[T] v_raw_rep;
  // vector[T] v_rep;
  // 
  // vector[T] h_rep;
  // matrix[T,3] s_rep;
  // matrix[T,3] m_rep;
  // matrix[T,3] y_rep;
  // 
  // vector[T] log_lik;
  // vector[T] log_lik_rep;
  // 
  // // pT_rep[1] = p;
  // // v_raw_tail_rep[1] = pT_rep[1] * sqrt(C2*l_rep[1]^2 ./ (C2 + pT_rep[1]^2 * l_rep[1]^2));
  // v_raw_rep[1] = v_raw_base_rep[1];
  // v_rep[1] = sigma*v_raw_rep[1];
  // h_rep[1] = v_rep[1] / sqrt(1 - phi^2);
  // 
  // s_rep[1,] = [exp((mu+h_rep[1])/2), sigma_obs[1], sigma_obs[2]];
  // m_rep[1,] = [mu_ret, mu_hloffset + (mu+h_rep[1])/2, (100/2)*((phi-1)*h_rep[1]+v_rep[1])];
  // 
  // y_rep[1,] = multi_normal_rng(m_rep[1,], quad_form_diag(R11_R12R21, s_rep[1,]'))';
  // 
  // for (t in 2:T) {
  //   // pT_rep[t] = (1-theta)*pT[t-1]+theta*v_raw_tail[t-1];
  //   // v_raw_tail_rep[t] = pT_rep[t] * sqrt(C2*l_rep[t]^2 ./ (C2 + pT_rep[t]^2 * l_rep[t]^2));
  //   v_raw_rep[t] = v_raw_base_rep[t];
  //   v_rep[t] = sigma*v_raw_rep[t];
  //   
  //   h_rep[t] = phi * h[t-1] + v_rep[t];
  //   
  //   s_rep[t,] = [exp((mu+h_rep[t])/2), sigma_obs[1], sigma_obs[2]];
  //   m_rep[t,] = [mu_ret,
  //                mu_hloffset + (mu+h_rep[t])/2,
  //                (100/2)*((phi-1)*h[t-1]+v_rep[t])] +
  //                   R[1:3,4]'*v_raw_rep[t];
  // 
  //   y_rep[t,] = multi_normal_rng(m_rep[t,], quad_form_diag(R11_R12R21, s_rep[t,]'))';
  // 
  //   log_lik[t] = multi_normal_lpdf(y[t,] | m[t,], quad_form_diag(R11_R12R21, s[t,]'));
  //   log_lik_rep[t] = multi_normal_lpdf(y[t,] | m_rep[t,], quad_form_diag(R11_R12R21, s_rep[t,]'));
  // }
}



