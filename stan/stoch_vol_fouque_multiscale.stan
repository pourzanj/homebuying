data {
  int<lower=0> T;   // # time points (equally spaced)
  vector[T] y;      // mean corrected return at time t
  
  int<lower=0> K;   // # number of AR time scales to model
}
transformed data {
  vector[K] ones = rep_vector(1.0, K);
}
parameters {
  real mu_ret;
  unit_vector[K] rho_unscaled;
  real<lower=0,upper=1> rho_length;
  
  real mu;                             // mean log volatility
  vector<lower=0, upper=1>[K] phi_raw; // persistence of volatility
  vector<lower=0>[K] sigma;            // white noise shock scale
  
  matrix[T, K] v;                      // std log volatility time t
}
transformed parameters {
  vector[K] rho = rho_length * rho_unscaled;
  vector[K] phi;
  matrix[T, K] h;
  vector[T] m;
  vector[T] s;
  
  // Set phi
  phi[1] = -1 + 2 * phi_raw[1];
  for(k in 2:K) {
    phi[k] = -1 + (phi[k-1] - -1) * phi_raw[k];
  }
  
  // Set h similar to Stan manual of SV model
  for(k in 1:K) {
    h[, k] = v[, k] * sigma[k];
  }
  
  h[1,] = h[1,] ./ sqrt(1 - phi .* phi)';     // rescale h[1]
  for (t in 2:T) {
    h[t,] += phi' .* h[t-1,];
  }
  
  // Set m and s
  {
    vector[T] s_marginal = exp((h * ones + mu) / 2);
    m = mu_ret + s_marginal .* (v * rho);
    s = s_marginal * sqrt(1 - rho' * rho);
  }
}
model {
  // prior
  phi ~ uniform(-1, 1);
  sigma ~ cauchy(0, 5);
  mu ~ cauchy(0, 10);

  // likelihood
  to_vector(v) ~ std_normal();
  y ~ normal(m, s);
}
generated quantities {
  // h five days ahead i.e. only using h[t-1]
  matrix[T, K] v_rep[5];
  matrix[T, K] h_rep[5];
  
  // mean and SD
  matrix[T, 5] m_rep;
  matrix[T, 5] s_rep;
  
  //standardized and actual return
  matrix[T, 5] eps_rep;
  matrix[T, 5] y_rep;
  
  // standarized error term of data
  vector[T] eps;
  
  /////////////////////////
  // First day sampled separately since there is no day before
  /////////////////////////
  
  // Draw h for first observation then 5 days into the future
  v_rep[1,1,] = to_vector(normal_rng(rep_array(1.0, K), 1))';
  h_rep[1,1,] = v_rep[1,1,] .* (sigma ./ sqrt(1 - phi .* phi))';
  for(l in 2:5) {
    v_rep[l,1,] = to_vector(normal_rng(rep_array(1.0, K), 1))';
    h_rep[l,1,] = phi' .* h_rep[l-1,1,] + sigma' .* v_rep[l,1,];
  }

  /////////////////////////
  // Rest of days 
  /////////////////////////
  
  // Draw h for rest of observations
  for (t in 2:T) {
    // Draw today based on yesterday, h[t-1]
    v_rep[1,t,] = to_vector(normal_rng(rep_array(1.0, K), 1))';
    h_rep[1,t,] = phi' .* h_rep[1,t,] + sigma' .* v_rep[1,t,];

    // Draw h 4 more days into the future based on yesterday
    for(l in 2:5) {
      v_rep[l,t,] = to_vector(normal_rng(rep_array(1.0, K), 1))';
      h_rep[l,t,] = phi' .* h_rep[l-1,t,] + sigma' .* v_rep[l,t,];
    }
  }
  
  // Draw actual returns based on h
  for(l in 1:5) {
    vector[T] s_marginal_rep = exp((h_rep[l,,] * ones + mu) / 2);
    m_rep[,l] = mu_ret + s_marginal_rep .* (v_rep[l,,] * rho);
    s_rep[,l] = s_marginal_rep * sqrt(1 - rho' * rho);
    
    eps_rep[,l] = to_vector(normal_rng(rep_vector(0, T), 1));
    y_rep[,l] = m_rep[,l] + s_rep[,l] .* eps_rep[,l];
  }
  
  /////////////////////////
  // Standardized actual data
  /////////////////////////
  eps = (y - m_rep[, 1]) ./ s_rep[, 1];
}
