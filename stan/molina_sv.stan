functions {
  real molina_lpdf(vector x, real y0, real m, real s) {
    real a = x[1];
    real b = x[2];
    real y = x[3];
    
    real k = 1/(sqrt(2*pi())*s^3)*exp(-(m^2-2*m*(y-y0))/(2*s^2));
  
    real ret = 0.0;
    for(n in 1:100) {
      real d1n = (y - y0 - 2*(-n)*(b-a))^2 / (2*s^2);
      real d2n = (y + y0 - 2*a - 2*(-n)*(b-a))^2 / (2*s^2);
      real neg = 4*n^2*(2*d1n-1)*exp(-d1n) - 4*-n*(-n-1)*(2*d2n-1)*exp(-d2n);
    
      real d1p = (y - y0 - 2*n*(b-a))^2 / (2*s^2);
      real d2p = (y + y0 - 2*a - 2*n*(b-a))^2 / (2*s^2);
      real pos = 4*n^2*(2*d1p-1)*exp(-d1p) - 4*n*(n-1)*(2*d2p-1)*exp(-d2p);
      
      real eps = (neg+pos);
      
      ret = ret + eps;
      if(fabs(eps) <= 1e-6) break;
      // if(n == 1000 || is_nan(ret)) print("Max iter reached", " iter: ", n, " k: ", k, " eps: ", eps, " m: ", m, " s: ", s, "a, b, c: " , x, " ret: ", ret, " neg: ", neg, " pos: ", pos);
    }
    
    if(a > min([y, y0])) print("a > min(y, y0)");
    if(b < max([y, y0])) print("b < max(y, y0)");
    if(a > b) print("a > b");
    
    // print("s : ", s, " y: ", y, " ", log(1/(sqrt(2*pi())*s^3)), " ", -(m^2-2*m*(y-y0))/(2*s^2), " ", log(ret));
    return log(1/(sqrt(2*pi())*s^3)) + -(m^2-2*m*(y-y0))/(2*s^2) + log(ret);
  }
}
data {
  int<lower=1> T;
  vector[T] y0;
  vector[T] a;
  vector[T] b;
  vector[T] y;
}
transformed data {
  // real<lower=-1, upper=1> phi = 0.0;
  // real<lower=0> sigma = 4.0;
  // real<lower=-1, upper=1> rho = 0.0; 
}
parameters {
  real alpha;
  real beta;
  
  // real gamma;
  
  // real<lower=0> s0;                     // mean log volatility
  real<lower=-1, upper=1> phi;  // persistence of volatility
  real mu;
  real<lower=0> sigma;
  real<lower=2> nu;
  real<lower=2> nuy;
  
  // real<lower=-1, upper=1> theta;
  // real<lower=0> tau;
  // real<lower=-1, upper=1> rhog;
  real<lower=-1, upper=1> rho;
  
  // vector[T] v_raw;
  vector<lower=0>[T] U;
  vector[T] h;
  // vector<lower=0>[T] s;
  
  // vector[T] u_raw;
}
transformed parameters {
  vector[T] u;
  vector[T] g;
  vector[T] v_raw;
  vector[T] eps;
  // // vector[T] h;
  vector<lower=0>[T] s;

  // u[1] = tau*u_raw[1];
  // g[1] = u[1] / sqrt(1-phi^2);
  
  v_raw[1] = (h[1]-mu)*(sqrt(1-phi^2)/(sigma*sqrt(nu/(nu/2))));

  for (t in 2:T) {
    // u[t] = tau*u_raw[t];
    // g[t] = theta*g[t-1] + u[t];

    // v_raw[t] = h[t]-h[t-1];
    v_raw[t] = (h[t]-mu*(1-phi)-phi*h[t-1])/(sigma*sqrt(nu/(nu/2)));
    // v[t] = sigma*v_raw[t];
    // h[t] = phi*h[t-1] + v[t];
  }
  // 
  s = sqrt(U) .* exp(h/2);
  eps = (y - alpha - beta*s + rho*v_raw .* s) ./ (s*sqrt(1-rho^2));
}
model {
  // u_raw ~ std_normal();
  
  // h[1] ~ normal(0.0, sigma/sqrt(1-phi^2));
  h[1] ~ normal(mu, sigma/sqrt(1-phi^2));
  [a[1], b[1], y[1]]' ~ molina(y0[1], alpha + beta*s[1] + rho*v_raw[1]*s[1], s[1]*sqrt(1-rho^2)); 
  
  nu ~ gamma(2, 0.1);
  nuy ~ gamma(2, 0.1);
  U ~ inv_gamma(nuy/2, nuy/2);
  // U ~ cauchy(0, 1);
  
  for(t in 2:T) {
    // h[t] ~ normal(phi*h[t-1], sigma);
    h[t] ~ student_t(nu, mu*(1-phi) + phi*h[t-1], sigma);
    [a[t], b[t], y[t]]' ~ molina(y0[t], alpha + beta*s[t] + rho*v_raw[t]*s[t], s[t]*sqrt(1-rho^2)); 
  }
}
generated quantities {
  vector[T] v_raw_rep = to_vector(student_t_rng(nu, rep_vector(0.0, T), 1.0));
  vector<lower=0>[T] U_rep = to_vector(inv_gamma_rng(rep_vector(nuy/2, T), nuy/2));
  vector[T] h_rep;
  vector[T] s_rep;
  vector[T] y_rep;

  h_rep[1] = mu + sigma*v_raw_rep[1] / sqrt(1-phi^2);
  s_rep[1] = sqrt(U_rep[1])*exp(h_rep[1]/2);
  y_rep[1] = normal_rng(alpha + beta*s[1] + rho*v_raw_rep[1]*s_rep[1], s_rep[1]*sqrt(1-rho^2));

  for (t in 2:T) {
    h_rep[t] = mu*(1-phi) + phi*h[t-1] + sigma*v_raw_rep[t];
    s_rep[t] = sqrt(U_rep[t])*exp(h_rep[t]/2);
    y_rep[t] = normal_rng(alpha + beta*s_rep[t] + rho*v_raw_rep[t]*s_rep[t], s_rep[t]*sqrt(1-rho^2));
  }
}
