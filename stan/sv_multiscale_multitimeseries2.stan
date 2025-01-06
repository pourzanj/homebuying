data {
  int<lower=0> T;   // # time points (equally spaced)
  matrix[T,3] y;
}
transformed data {
  // real mu_ret = 0.05;
  // real<lower=0> mu_hloffset = 0.25;

  // real mu = -0.75;                     // mean log volatility
  // real<lower=-1, upper=1> phi = 0.97;  // persistence of volatility
  // real<lower=0> sigma = 0.2;
  real a = 0;
  real b = 0;
  real c = 0;
  // real<lower=-1, upper=1> phi = 0.0;
}
parameters {
  real mu_ret;
  real<lower=0> mu_hloffset;
  
  // real a;
  // real b;
  // real c;
  // real alpha;
  real<lower=2> nu;

  real mu;                     // mean log volatility
  real<lower=-1, upper=1> phi;  // persistence of volatility
  real<lower=0> sigma;         // white noise shock scale
  
  real<lower=0> sigma_obs[2];
  
  cholesky_factor_corr[4] L;
  
  vector[T] v_raw;
  // vector<lower=0>[T] v_raw_exp;
}
transformed parameters {
  corr_matrix[4] R = multiply_lower_tri_self_transpose(L);
  cov_matrix[3] R11_R12R21 = R[1:3,1:3] - R[1:3,4]*R[4,1:3];
  matrix[3,3] L11_L12L21 = cholesky_decompose(R11_R12R21);
  
  vector[T] v; 
  vector[T] h;
  matrix[T,3] s;
  vector[T] S;
  matrix[T,3] m;
  
  v[1] = sigma*exp(b*v_raw[1])*(v_raw[1]);
  h[1] = v[1] / sqrt(1 - phi^2);
  s[1,] = [exp((mu+h[1])/2), sigma_obs[1], sigma_obs[2]];
  m[1,] = [mu_ret, mu_hloffset + (mu+h[1])/2, (100/2)*((phi-1)*h[1]+v[1])];
  for (t in 2:T) {
    v[t] = sigma*exp(a*h[t-1]/2 + b*v_raw[t] + c*v_raw[t-1])*(v_raw[t]);
    h[t] = phi * h[t-1] + v[t];
    S[t] = exp((mu+h[t])/2);
    s[t,] = [S[t], sigma_obs[1], sigma_obs[2]];
    m[t,] = [mu_ret,
             mu_hloffset + (mu+h[t])/2,
             (100/2)*(h[t]-h[t-1])] +
                s[t,] .* R[1:3,4]'*exp(a*h[t-1]/2 + b*v_raw[t] + c*v_raw[t-1])*(v_raw[t]);
  }
  
}
model {
  // prior
  phi ~ uniform(-1, 1);
  sigma ~ exponential(0.1);
  // mu ~ normal(-5, 2.5);
  mu ~ normal(0, 1);

  // a ~ normal(0, 1.0);
  // b ~ normal(0, 1.0);
  // c ~ normal(0, 1.0);

  mu_ret ~ normal(0.1, 0.05);
  mu_hloffset ~ normal(0.25, 0.025);

  sigma_obs[1] ~ exponential(1.0);
  sigma_obs[2] ~ exponential(1e-2);

  L ~ lkj_corr_cholesky(10.0);
  // (R[1,4]+1)/2.0 ~ beta(3, 10);

  // likelihood
  // v_raw ~ std_normal();
  nu ~ gamma(2, 0.1);
  v_raw ~ student_t(nu, 0.0, 1.0);
  // alpha ~ normal(0, 1.0);
  // v_raw ~ skew_normal(0, 1, alpha);
  for (t in 1:T) {
    y[t,] ~ multi_normal(m[t,], quad_form_diag(R11_R12R21, s[t,]'));
  }
}
generated quantities {
  vector[T] v_raw_rep = to_vector(normal_rng(rep_vector(0.0, T), 1.0));
  vector[T] v_rep;
  vector[T] h_rep;
  matrix[T,3] s_rep;
  vector[T] S_rep;
  matrix[T,3] m_rep;
  matrix[T,3] y_rep;

  vector[T] log_lik;
  vector[T] log_lik_rep;

  v_rep[1] = sigma*v_raw_rep[1];
  h_rep[1] = v_rep[1] / sqrt(1 - phi^2);
  
  s_rep[1,] = [exp((mu+h_rep[1])/2), sigma_obs[1], sigma_obs[2]];
  m_rep[1,] = [mu_ret, mu_hloffset + (mu+h_rep[1])/2, (100/2)*((phi-1)*h_rep[1]+v_rep[1])];
  y_rep[1,] = multi_normal_rng(m_rep[1,], quad_form_diag(R11_R12R21, s_rep[1,]'))';
  
  for (t in 2:T) {
    v_rep[t] = sigma*exp(a*h[t-1]/2+b*v_raw_rep[t] + c*v_raw_rep[t-1])*v_raw_rep[t];
    h_rep[t] = phi * h[t-1] + v_rep[t];
    S_rep[t] = exp((mu+h_rep[t])/2);
    s_rep[t,] = [S_rep[t], sigma_obs[1], sigma_obs[2]];
    m_rep[t,] = [mu_ret,
                 mu_hloffset + (mu+h_rep[t])/2,
                 (100/2)*(h_rep[t]-h[t-1])] +
                    s_rep[t,] .* R[1:3,4]'*exp(a*h[t-1]/2 + b*v_raw_rep[t] + c*v_raw[t-1])*v_raw_rep[t];

    y_rep[t,] = multi_normal_rng(m_rep[t,], quad_form_diag(R11_R12R21, s_rep[t,]'))';

    log_lik[t] = multi_normal_lpdf(y[t,] | m[t,], quad_form_diag(R11_R12R21, s[t,]'));
    log_lik_rep[t] = multi_normal_lpdf(y[t,] | m_rep[t,], quad_form_diag(R11_R12R21, s_rep[t,]'));
  }
}



