data {
  int<lower=0> T;   // # time points (equally spaced)
  vector[T] y;      // mean corrected return at time t
}
parameters {
  real mu;                     // mean log volatility
  real<lower=-1,upper=1> phi;  // persistence of volatility
  real<lower=0> sigma;         // white noise shock scale
  
  real mu_ret;
  real<lower=-1,upper=1> rho;
  
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
  // h one day ahead i.e. only using h[t-1]
  vector[T] h_std_rep;
  vector[T] h_rep;
  
  // mean and SD of y one day ahead
  vector[T] m_rep;
  vector[T] s_rep;
  
  // y one day ahead only using today's info
  vector[T] eps_rep;
  vector[T] y_rep;
  
  // Draw h
  h_std_rep[1] = normal_rng(0, 1);
  h_rep[1] = mu + (sigma / sqrt(1 - phi * phi)) * h_std_rep[1];
  for (t in 2:T) {
    h_std_rep[t] = normal_rng(0, 1);
    h_rep[t] = mu + phi*(h[t-1]-mu) + sigma * h_std_rep[t]; 
  }
   
  // Set mean and SD based on h
  m_rep = mu_ret + rho * exp(h_rep / 2) .* h_std_rep;
  s_rep = sqrt(1-rho^2) * exp(h_rep / 2);
  
  // Set eps and y
  eps_rep = to_vector(normal_rng(rep_vector(0, T), 1));
  y_rep = m_rep + s_rep .* eps_rep;
}
