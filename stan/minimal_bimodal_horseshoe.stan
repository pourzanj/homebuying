data {
  real y;
  real<lower=0.0> sigma;
}
parameters {
  real<lower=0.0> lambda;
}
model {
  lambda ~ cauchy(0, 1);
  y ~ normal(0, sqrt(sigma^2 + lambda^2));
}

