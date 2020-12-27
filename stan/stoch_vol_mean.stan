data {
  int<lower=0> T;   // # time points (equally spaced)
  vector[T] y;      // mean corrected return at time t
}
parameters {
  real<lower=0> a;
  real<upper=0> b;
  real<lower=0> c;
  
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
  y ~ normal(c + a * exp(b * h / 2), exp(h / 2));
}
generated quantities {
  real y_rep[T] = normal_rng(a * exp((b+1) * h / 2), exp(h / 2));
  vector[T] h_rep;
  real y_rep_rep[T];
  vector[T] m_rep;
  vector[T] s_rep; 
  // real h_ahead = mu + phi * (h[T] - mu) + normal_rng(0, sigma);
  // real y_ahead = normal_rng(0, exp(h_ahead / 2));
  
  h_rep[1] = normal_rng(mu, sigma / sqrt(1 - phi * phi));
  for (t in 2:T)
    h_rep[t] = mu + phi*(h[t-1]-mu) + normal_rng(0, sigma);
    
  m_rep = c + a * exp(b * h / 2);
  s_rep = exp(h_rep / 2);
  y_rep_rep = normal_rng(m_rep, s_rep);
}
