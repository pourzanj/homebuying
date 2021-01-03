data {
  int<lower=0> T;   // # time points (equally spaced)
  vector[T] y;      // mean corrected return at time t
}
parameters {
  real mu_ret;
  real<lower=-1,upper=1> rho;
  
  real mu;                     // mean log volatility
  real<lower=-1,upper=1> phi;  // persistence of volatility
  real<lower=0> sigma;         // white noise shock scale
  
  vector[T] h_std;             // std log volatility time t
}
transformed parameters {
  // h_{t+1} ~ N(mu + phi*(h_t - mu) , sigma)
  // h_{t+1} ~ N(phi*h_t + (1-phi)*mu, sigma)
  vector[T] h = h_std * sigma;  // now h ~ normal(0, sigma)
  h[1] /= sqrt(1 - phi * phi);  // rescale h[1]
  h += mu;
  for (t in 2:T)
    h[t] += phi * (h[t-1] - mu);
}
model {
  // prior
  phi ~ uniform(-1, 1);
  sigma ~ cauchy(0, 5);
  mu ~ cauchy(0, 10);

  // likelihood
  h_std ~ std_normal();
  y ~ normal(mu_ret + rho * exp(h / 2) .* h_std, sqrt(1-rho^2) * exp(h / 2));
}
generated quantities {
  // h five days ahead i.e. only using h[t-1]
  matrix[T, 5] h_std_rep;
  matrix[T, 5] h_rep;
  
  // mean and SD
  matrix[T, 5] m_rep;
  matrix[T, 5] s_rep;
  
  //standardized and actual return
  matrix[T, 5] eps_rep;
  matrix[T, 5] y_rep;
  
  // standarized error term of data
  vector[T] eps;
  
  // Leave-forward-out one day ahead prediction
  // h five days ahead i.e. only using h[t-1]
  real h_std_lfo;
  real h_lfo;
  real m_lfo;
  real s_lfo;
  real eps_lfo;
  real y_lfo;
  
  /////////////////////////
  // First day sampled separately since there is no day before
  /////////////////////////
  
  // Draw h for first observation then 5 days into the future
  h_std_rep[1, 1] = normal_rng(0, 1);
  h_rep[1, 1] = mu + (sigma / sqrt(1 - phi * phi)) * h_std_rep[1, 1];
  for(l in 2:5) {
    h_std_rep[1, l] = normal_rng(0, 1);
    h_rep[1, l] = mu + phi*(h_rep[1, l-1]-mu) + sigma * h_std_rep[1, l];
  }
  
  m_rep[1,] = mu_ret + rho * exp(h_rep[1,] / 2) .* h_std_rep[1,];
  s_rep[1,] = sqrt(1-rho^2) * exp(h_rep[1,] / 2);
  
  eps_rep[1,] = to_vector(normal_rng(rep_vector(0, 5), 1))';
  y_rep[1,] = m_rep[1,] + s_rep[1,] .* eps_rep[1,];
  
  /////////////////////////
  // Rest of days 
  /////////////////////////
  
  // Draw h for rest of observations
  for (t in 2:T) {
    // Draw today based on yesterday, h[t-1]
    h_std_rep[t, 1] = normal_rng(0, 1);
    h_rep[t, 1] = mu + phi*(h[t-1]-mu) + sigma * h_std_rep[t, 1];
    
    // Draw h 4 more days into the future based on yesterday
    for(l in 2:5) {
      h_std_rep[t, l] = normal_rng(0, 1);
      h_rep[t, l] = mu + phi*(h_rep[t, l-1]-mu) + sigma * h_std_rep[t, l];
    }
    
    // Draw actual returns based on h
    m_rep[t,] = mu_ret + rho * exp(h_rep[t,] / 2) .* h_std_rep[t,];
    s_rep[t,] = sqrt(1-rho^2) * exp(h_rep[t,] / 2);
  
    eps_rep[t,] = to_vector(normal_rng(rep_vector(0, 5), 1))';
    y_rep[t,] = m_rep[t,] + s_rep[t,] .* eps_rep[t,];
  }
  
  /////////////////////////
  // Standardized actual data
  /////////////////////////
  eps = (y - m_rep[, 1]) ./ s_rep[, 1];
  
  /////////////////////////
  // LFO
  /////////////////////////
  h_std_lfo = normal_rng(0, 1);
  h_lfo = mu + phi*(h[T]-mu) + sigma * h_std_lfo;

  m_lfo = mu_ret + rho * exp(h_lfo / 2) .* h_std_lfo;
  s_lfo = sqrt(1-rho^2) * exp(h_lfo / 2);
  
  eps_lfo =  normal_rng(0, 1);
  y_lfo = m_lfo + s_lfo .* eps_lfo;
}
