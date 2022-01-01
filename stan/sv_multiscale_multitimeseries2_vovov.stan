data {
  int<lower=0> T;   // # time points (equally spaced)
  matrix[T,3] y;
}
parameters {
  real mu_ret[2];
  real<lower=0> mu_hloffset;

  real<lower=-1, upper=1> psi;
  real<lower=0> kappa;

  real<lower=-1, upper=1> theta;
  real<lower=0> tau;
  
  real mu;                     // mean log volatility
  real<lower=-1, upper=1> phi;  // persistence of volatility
  real<lower=0> sigma;         // white noise shock scale
  
  real<lower=0> sigma_obs[2];
  
  cholesky_factor_corr[6] L;
  
  vector[T] w_raw;
  vector[T] u_raw;
  vector[T] v_raw;
}
transformed parameters {
  // Corr matrices and sub-matrices
  matrix[6,6] R = multiply_lower_tri_self_transpose(L);
  // corr_matrix[5] Rf = R[1:5,1:5] - R[1:5,6]*inv(R[6,6])*R[6,1:5];
  // corr_matrix[4] Rg = Rf[1:4,1:4] - Rf[1:4,5]*inv(Rf[5,5])*Rf[5,1:4];
  // corr_matrix[3] Rh = Rf[1:3,1:3] - Rg[1:3,4]*inv(Rg[4,4])*Rg[4,1:3];
  matrix[5,5] Rf = R[1:5,1:5] - quad_form_sym([[inv(R[6,6])]], [R[6,1:5]]);
  matrix[4,4] Rg = Rf[1:4,1:4] - quad_form_sym([[inv(Rf[5,5])]], [Rf[5,1:4]]);
  matrix[3,3] Rh = Rg[1:3,1:3] - quad_form_sym([[inv(Rg[4,4])]], [Rg[4,1:3]]);

  // Declare time variables
  vector[T] w;
  vector[T] f;
  
  vector[T] u;
  vector[T] g;
  
  vector[T] v; 
  vector[T] h;
  
  // SD of observations
  vector[T] S;
  matrix[T,3] s;
  
  // Mean of observations (separate by caused by corr and just for the data)
  matrix[T,3] m_corr;
  matrix[T,3] m_data;
  matrix[T,3] m;
  
  // print("~~~~~~~~~~~~~~~Rh: ");
  // print(Rh[1,]);
  // print(Rh[2,]);
  // print(Rh[3,]);
  
  // Set first time
  w[1] = kappa*w_raw[1];
  f[1] = w[1] / sqrt(1-psi^2);
  
  u[1] = tau*exp(f[1]/2) * (R[5,6]*w_raw[1] + sqrt(Rf[5,5])*u_raw[1]);
  g[1] = u[1] / sqrt(1-theta^2);
  
  v[1] = sigma*exp(g[1]/2) * (R[4,6]*w_raw[1] + Rf[4,5]*inv_sqrt(Rf[5,5])*u_raw[1] + sqrt(Rg[4,4])*v_raw[1]);
  h[1] = v[1] / sqrt(1 - phi^2);
  
  S[1] = exp((mu+h[1])/2);
  s[1,] = [S[1], sigma_obs[1], sigma_obs[2]];
  m_corr[1,] = (R[1:3,6]*w_raw[1] + Rf[1:3,5]*inv_sqrt(Rf[5,5])*u_raw[1] + Rg[1:3,4]*inv_sqrt(Rg[4,4])*v_raw[1])';
  m_data[1,] = [mu_ret[1] + mu_ret[2]*S[1],
                mu_hloffset + (mu+h[1])/2,
                (100/2)*v[1]];
  m[1,] = m_data[1,] + s[1,] .* m_corr[1,];

  // Set rest
  for (t in 2:T) {
    
    w[t] = kappa*w_raw[t];
    f[t] = psi*f[t-1] + w[t];
  
    u[t] = tau*exp(f[t]/2) * (R[5,6]*w_raw[t] + sqrt(Rf[5,5])*u_raw[t]);
    g[t] = theta*g[t-1] + u[t];
  
    v[t] = sigma*exp(g[t]/2) * (R[4,6]*w_raw[t] + Rf[4,5]*inv_sqrt(Rf[5,5])*u_raw[t] + sqrt(Rg[4,4])*v_raw[t]);
    h[t] = phi*h[t-1] + v[t];
  
    S[t] = exp((mu+h[t])/2);
    s[t,] = [S[t], sigma_obs[1], sigma_obs[2]];
    m_corr[t,] = (R[1:3,6]*w_raw[t] + Rf[1:3,5]*inv_sqrt(Rf[5,5])*u_raw[t] + Rg[1:3,4]*inv_sqrt(Rg[4,4])*v_raw[t])';
    m_data[t,] = [mu_ret[1] + mu_ret[2]*S[t],
                  mu_hloffset + (mu+h[t])/2,
                  (100/2)*(h[t] - h[t-1])];
    m[t,] = m_data[t,] + s[t,] .* m_corr[t,];
  }
  
}
model {
  // Prior
  psi ~ uniform(-1, 1);
  kappa ~ exponential(1.0);
  
  theta ~ uniform(-1, 1);
  tau ~ exponential(1.0);
  
  phi ~ uniform(-1, 1);
  sigma ~ exponential(1.0);
  mu ~ normal(-1, 2);

  mu_ret ~ normal(0.1, 0.05);
  mu_hloffset ~ normal(0.25, 0.025);

  sigma_obs[1] ~ exponential(1.0);
  sigma_obs[2] ~ exponential(1e-2);

  L ~ lkj_corr_cholesky(10.0);

  // Likelihood
  w_raw ~ std_normal();
  u_raw ~ std_normal();
  v_raw ~ std_normal();

  for (t in 2:T) {
    y[t,] ~ multi_normal(m[t,], quad_form_diag(Rh, s[t,]'));
  }
}
generated quantities {
  vector[T] w_raw_rep = to_vector(normal_rng(rep_vector(0.0, T), 1.0));
  vector[T] u_raw_rep = to_vector(normal_rng(rep_vector(0.0, T), 1.0));
  vector[T] v_raw_rep = to_vector(normal_rng(rep_vector(0.0, T), 1.0));
 
  vector[T] w_rep;
  vector[T] f_rep;
  
  vector[T] u_rep;
  vector[T] g_rep;
  
  vector[T] v_rep;
  vector[T] h_rep;
  
  vector[T] S_rep;
  matrix[T,3] s_rep;
  matrix[T,3] m_corr_rep;
  matrix[T,3] m_data_rep;
  matrix[T,3] m_rep;
  matrix[T,3] y_rep;

  vector[T] log_lik;
  vector[T] log_lik_rep;

  w_rep[1] = kappa*w_raw_rep[1];
  f_rep[1] = w_rep[1] / sqrt(1-psi^2);
  
  u_rep[1] = tau*exp(f_rep[1]/2) * (R[5,6]*w_raw_rep[1] + sqrt(Rf[5,5])*u_raw_rep[1]);
  g_rep[1] = u_rep[1] / sqrt(1-theta^2);
  
  v_rep[1] = sigma*exp(g_rep[1]/2) * (R[4,6]*w_raw_rep[1] + Rf[4,5]*inv_sqrt(Rf[5,5])*u_raw_rep[1] + sqrt(Rg[4,4])*v_raw_rep[1]);
  h_rep[1] = v_rep[1] / sqrt(1 - phi^2);
  
  S_rep[1] = exp((mu+h_rep[1])/2);
  s_rep[1,] = [S_rep[1], sigma_obs[1], sigma_obs[2]];
  m_corr_rep[1,] = (R[1:3,6]*w_raw_rep[1] + Rf[1:3,5]*inv_sqrt(Rf[5,5])*u_raw_rep[1] + Rg[1:3,4]*inv_sqrt(Rg[4,4])*v_raw_rep[1])';
  m_data_rep[1,] = [mu_ret[1] + mu_ret[2]*S_rep[1],
                    mu_hloffset + (mu+h_rep[1])/2,
                    (100/2)*v_rep[1]];
  m_rep[1,] = m_data_rep[1,] + s_rep[1,] .* m_corr_rep[1,];
  
  y_rep[1,] = multi_normal_rng(m_rep[1,], quad_form_diag(Rh, s_rep[1,]'))';
  
  for (t in 2:T) {
    
    w_rep[t] = kappa*w_raw_rep[t];
    f_rep[t] = psi*f[t-1] + w_rep[t];
  
    u_rep[t] = tau*exp(f_rep[t]/2) * (R[5,6]*w_raw_rep[t] + sqrt(Rf[5,5])*u_raw_rep[t]);
    g_rep[t] = theta*g[t-1] + u_rep[t];
  
    v_rep[t] = sigma*exp(g_rep[t]/2) * (R[4,6]*w_raw_rep[t] + Rf[4,5]*inv_sqrt(Rf[5,5])*u_raw_rep[t] + sqrt(Rg[4,4])*v_raw_rep[t]);
    h_rep[t] = phi*h[t-1] + v_rep[t];
  
    S_rep[t] = exp((mu+h_rep[t])/2);
    s_rep[t,] = [S_rep[t], sigma_obs[1], sigma_obs[2]];
    m_corr_rep[t,] = (R[1:3,6]*w_raw_rep[t] + Rf[1:3,5]*inv_sqrt(Rf[5,5])*u_raw_rep[t] + Rg[1:3,4]*inv_sqrt(Rg[4,4])*v_raw_rep[t])';
    m_data_rep[t,] = [mu_ret[1] + mu_ret[2]*S_rep[t],
                      mu_hloffset + (mu+h_rep[t])/2,
                      (100/2)*(h_rep[t] - h[t-1])];
    m_rep[t,] = m_data_rep[t,] + s_rep[t,] .* m_corr_rep[t,];

    y_rep[t,] = multi_normal_rng(m_rep[t,], quad_form_diag(Rh, s_rep[t,]'))';

    log_lik[t] = multi_normal_lpdf(y[t,] | m[t,], quad_form_diag(Rh, s[t,]'));
    log_lik_rep[t] = multi_normal_lpdf(y[t,] | m_rep[t,], quad_form_diag(Rh, s_rep[t,]'));
  }
}



