data {
  int<lower=0> T;
  real r[T];
}
parameters {
  real mu;
  real<lower=0> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=(1-alpha1)> beta1;
  real<lower=0> sigma1;
}
transformed parameters {
  real<lower=0> sigma[T];
  sigma[1] = sigma1;
  for (t in 2:T)
    sigma[t] = sqrt(alpha0
                     + alpha1 * pow(r[t-1] - mu, 2)
                     + beta1 * pow(sigma[t-1], 2));
}
model {
  r ~ normal(mu, sigma);
}
generated quantities {
  real r_rep[T];
  real<lower=0> sigma_rep[T];
  sigma_rep[1] = sigma1;
  r_rep[1] = normal_rng(mu, sigma_rep[1]);
  
  for (t in 2:T) {
    sigma_rep[t] = sqrt(alpha0
                        + alpha1 * pow(r_rep[t-1] - mu, 2)
                        + beta1 * pow(sigma_rep[t-1], 2));
                        
    r_rep[t] = normal_rng(mu, sigma_rep[t]);
    
  }
    
}
