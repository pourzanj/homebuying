data {
  int<lower=1> N;
  vector[N] y;
  real<lower=0> tau0;
  real<lower=0> nu;
}

parameters {
  real mu;
  real<lower=0> sigma;
  real lambda;
  // real<lower=0> nu;
  real<lower=0> tau;
  
  vector<lower=0>[N] c2;
  vector<lower=0>[N] l;
}

transformed parameters {
  vector<lower=0>[N] lt = tau * sqrt((c2 .* l^2) ./ (c2 + tau^2 * l^2));
}

model {
  c2 ~ inv_gamma(nu/2, nu/2);
  l ~ cauchy(0, 1);
  
  nu ~ gamma(2, 0.1);
  tau ~ cauchy(0, tau0);
  
  y ~ skew_normal(mu, sigma*lt, lambda*lt);
}