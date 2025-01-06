data {
  int<lower=0> T;   // # time points (equally spaced)
  vector[T] y;      // mean corrected return at time t
}

parameters {
  real mu_ret;
  
  real mu;                     // mean log volatility
  real<lower=-1, upper=1> phi; // persistence of volatility
  real<lower=0> sigma;         // white noise shock scale
  
  real<lower=5> nu;
  
  vector[T] v;  // std log volatility time t
  vector<lower=0>[T] l2;
}

transformed parameters {
  vector[T] h;

  h[1] = mu + (sigma/sqrt(1-phi^2))*(sqrt(nu-2)/nu)*sqrt(l2[1])*v[1];
  
  for (t in 2:T) {
    h[t] = mu + phi*(h[t-1]-mu) + sigma*(sqrt(nu-2)/nu)*sqrt(l2[t])*v[t];
  }
}

model {
  mu_ret ~ normal(0, 1);
  phi ~ uniform(-1, 1);
  sigma ~ cauchy(0, 5);
  mu ~ cauchy(0, 10);

  v ~ std_normal();
  l2 ~ inv_gamma(nu/2.0, nu/2.0);
  nu ~ normal(20, 5);
  // v ~ student_t(nu, 0.0, (nu-2)/nu);
  
  y ~ normal(mu_ret, exp(h / 2));
}
