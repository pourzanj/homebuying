// Generalized T (GT) with Skew added using Fernandez and Steel (1998) method
data {
  int<lower=0> N;
  vector[N] y;
}
parameters {
  real<lower=0> v;
  real<lower=0> nu;
  real<lower=0> gamma;
}
model {
  y ~ normal(mu, sigma);
}

