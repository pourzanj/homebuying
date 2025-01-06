data {
  int<lower=0> T;    // # time points (equally spaced)
  vector[T] logVIX_Obs;
  
  int<lower=2> H; // inverse of stepsize. number of Euler steps per day
}
transformed data {
  real<lower=0> h = 1.0 / H; // step size
}
parameters {
  real<lower=-1, upper=1> kappa;
  real theta;
  real<lower=-1, upper=1> kappa_nu;
  real theta_nu;
  real<lower=0> sigma_nu;
  real<lower=-1, upper=1> rho;
  
  real<lower=0> V0;
  
  matrix[T,H-1] W_raw;
  matrix[T,H] W_nu;
}
transformed parameters {
  matrix[T,H-1] W;
  
  // print("~~~~~~~~~~~~~~~~~~~~~`");
  // print("kappa: ", kappa);
  // print("theta: ", theta);
  // print("kappa_nu: ", kappa_nu);
  // print("theta_nu: ", theta_nu);
  // print("sigma_nu: ", sigma_nu);
  // print("rho: ", rho);
  
  // Goth observed and unobserved intraday VIX levels. First column is observed
  // levels at close. Rest of columns are latent up until next day's close.
  matrix[T,H] logVIX;
  
  matrix[T,H] V;
  
  for(t in 1:T) {
    logVIX[t, 1] = logVIX_Obs[t];
    
    if(t != 1) {
      V[t,1] = V[t-1,H] + kappa_nu*(theta_nu-V[t-1,H])*h + sqrt(h)*sigma_nu*sqrt(V[t-1,H])*W_nu[t-1,H];
    } else {
      V[t,1] = V0;
    }
    
    for(i in 2:H) {
      W[t,i-1] = rho*W_nu[t,i-1] + (1-rho^2)*W_raw[t,i-1];
      
      logVIX[t, i] = logVIX[t,i-1] + kappa*(theta-logVIX[t,i-1])*h + sqrt(h)*sqrt(V[t,i-1])*W[t,i-1];
      V[t, i] = V[t,i-1] + kappa_nu*(theta_nu-V[t,i-1])*h + sqrt(h)*sigma_nu*sqrt(V[t,i-1])*W_nu[t,i-1];
    }
    // print("~~~~~~~");
    // print("W: ", W[t,]);
    // print("logVIX: ", logVIX[t,]);
    // print("V: ", V[t,]);
  }
}
model {
  // Prior
  
  
  // Likelihood
  to_vector(W_raw) ~ std_normal();
  to_vector(W_nu) ~ std_normal();
  
  for(t in 2:T) {
    logVIX_Obs[t] ~ normal(logVIX[t-1,H] + kappa*(theta-logVIX[t-1,H])*h, sqrt(h)*sqrt(V[t-1,H]));
  }
}