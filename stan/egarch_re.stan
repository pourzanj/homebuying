functions {
   /**
    * Return natural cubic spline basis function at a specific knot evaluated
    * at a set of x-values
    *
    * @param xs An array of reals of x-values for which the basis function should be evaluated
    * @param l A real denoting \lambda_j for the j'th basis function in Royston & Parmar
    * @param k A real denoting k_j, i.e. the position of the j'th knot
    * @param kmax A real denoting the position of the right boundary knot
    * @param kmin A real denoting the position of the left boundary knot
    * @return A vector corresponding to the values of the spline basis function evaluated at x
    */
    vector basis_function_v(real[] xs,  real l, real k, real kmax, real kmin) {
        vector[size(xs)] vs;
        for(i in 1:size(xs)) 
            vs[i] = pow(max({0, xs[i] - k}),3) - l*pow(max({0, xs[i] - kmin}),3)s - (1-l)*pow(max({0, xs[i] - kmax}),3);
        return vs;
}
data {
  int<lower=0> T;
  vector[T] r;
  real<lower=0> sigma1;
  
  real tau0;
  real nu_global;
  real nu;
  real s;
}
transformed data {
  // vector<lower=0>[T] lambda = rep_vector(1, T);
  // real<lower=0> tau = tau0;
}
parameters {
  
  
  // GARCH parameters
  real mu;
  real alpha0;
  real alpha1;
  real alpha2;
  real beta1;
  
  // Models negative correlation between returns and innovations in volatility
  real rho;
  
  // RHS
  real<lower=0> tau;
  real<lower=0> c_sq_aux;
  vector[T] v_aux;
  vector<lower=0>[T] lambda;
}
transformed parameters {
  
  vector[T] h_expected;
  vector[T] h;      // log of volatility squared
  vector[T] sigma;  // volatility (SD of returns)
  vector[T] e;
  matrix[T, 4] x;   // design matrix

  // Set RHS
  real c = s * sqrt(c_sq_aux);
  vector<lower=0>[T] lambda_tilde = sqrt(c^2 * square(lambda) ./ (c^2 + tau^2*square(lambda)));
  vector[T] v = v_aux .* (tau * lambda_tilde);
  
  // Set variance
  e[1] = 0;
  h_expected[1] = log(sigma1^2);
  h[1] = h_expected[1] + v[1];
  sigma[1] = exp(0.5 * h[1]);
  for (t in 2:T) {
    e[t] = (r[t-1] - mu - rho * v[t-1]) / sigma[t-1];
    
    x[t, 1] = 1;
    x[t, 2] = e[t];
    x[t, 3] = e[t]^2 / (1 + e[t]^2);
    x[t, 4] = h[t-1];

    h_expected[t] = x[t,] * [alpha0, alpha1, alpha2, beta1]';
    h[t] = h_expected[t] + v[t];
    sigma[t] = exp(0.5 * h[t]);
  }
  
  
}
model {
  // Prior on parameters
  mu ~ normal(0, 0.2);
  rho ~ normal(0, 1);
  
  // RHS
  v_aux ~ std_normal();
  c_sq_aux ~ inv_gamma(0.5 * nu, 0.5 * nu);
  tau ~ student_t(nu_global, 0, tau0);
  lambda ~ cauchy(0, 1);
  
  // Likelihood
  r ~ normal(mu + rho * v, sigma);
}
generated quantities {
  vector<lower=0>[T] lambda_rep = fabs(to_vector(cauchy_rng(rep_vector(0, T), 1)));
  vector<lower=0>[T] lambda_tilde_rep = sqrt(c^2 * square(lambda_rep) ./ (c^2 + tau^2*square(lambda_rep)));
  vector[T] v_rep = to_vector(normal_rng(rep_vector(0, T), 1)) .* (tau * lambda_tilde);
  vector[T] r_rep = to_vector(normal_rng(mu + rho * v_rep, sigma .* exp(0.5 * v_rep)));
}
