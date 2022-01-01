functions {
  real ast_lpdf(vector y, real mu, real sigma, real alpha, real nu1, real nu2)
  { // Skewed generalised t
    int N = dims(y)[1];
    real Knu1 = exp(lgamma((nu1+1)/2) - 0.5*log(pi()*nu1) - lgamma(nu1/2));
    real Knu2 = exp(lgamma((nu2+1)/2) - 0.5*log(pi()*nu2) - lgamma(nu2/2));
    real out = 0.0;
    
    for (n in 1:N) {
      if(y[n] <= mu) {
        out += -log(sigma) + -(nu1+1)/2 * log1p((1/nu1)*square((y[n]-mu)/(2*alpha*sigma*Knu1)));
      } else {
        out += -log(sigma) + -(nu2+1)/2 * log1p((1/nu2)*square((y[n]-mu)/(2*(1-alpha)*sigma*Knu2)));
      }
    }
    
    return out;
  }
  
  real sign(real x) {
    if(x > 0) {
      return 1;
    } else if(x < 0){
      return -1;
    } else {
      return 0;
    }
  } 
}
data {
  int<lower=0> N;
  vector[N] y;
}
parameters {
  real mu;
  real<lower=0> sigma;
  real<lower=0, upper = 1> alpha;
  real<lower=2> nu1;
  real<lower=2> nu2;
}
model {
  // prior
  nu1 ~ gamma(2, 0.1);
  nu2 ~ gamma(2, 0.1);

  // likelihood
  y ~ ast(mu, sigma, alpha, nu1, nu2);
}

generated quantities {
  vector[N] U = to_vector(uniform_rng(rep_vector(0.0, N), 1.0));
  vector[N] T1 = to_vector(student_t_rng(nu1, rep_vector(0.0, N), 1.0));
  vector[N] T2 = to_vector(student_t_rng(nu2, rep_vector(0.0, N), 1.0));
  vector[N] y_rep;
  
  real r;
  
  real Knu1 = exp(lgamma((nu1+1)/2) - 0.5*log(pi()*nu1) - lgamma(nu1/2));
  real Knu2 = exp(lgamma((nu2+1)/2) - 0.5*log(pi()*nu2) - lgamma(nu2/2));
  real astar = alpha*Knu1 / (alpha*Knu1 + (1-alpha)*Knu2);
  
  for(n in 1:N) {
    y_rep[n] = mu + sigma*(astar*fabs(T1[n])*(sign(U[n]-alpha)-1) + (1-astar)*fabs(T2[n])*(sign(U[n]-alpha)+1));  
  }
  
  r = prod(exp(y_rep/100))^(252.0/N);
}

