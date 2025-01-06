data {
  int<lower=0> T;   // # time points (equally spaced)
  vector[T] y;      // mean corrected return at time t
}

transformed data{
  real a = 3.2;
}

parameters {
  real mu_ret;
  
  real mu;                     // mean log volatility
  // real<lower=-1, upper=1> phi; // persistence of volatility
  real<lower=0> sigma;         // white noise shock scale
  // vector[T] h;                 // log volatility at time t
  
  // real a;
  real<lower=-1, upper=1> b;
  real<lower=0> c;
  
  // real v0;
  vector[T] v;  // std log volatility time t
  
  // real w0;
  vector[T] w;
}

transformed parameters {
  // real g0;
  vector[T] g;
  vector[T] phi;
  
  real h0;
  vector[T] h;
  
  // g0 = c*w0;
  // g[1] = a + b*(g0-a) + c*sqrt(1-b^2)*w[1];
  g[1] = a + c*w[1];
  phi[1] = inv_logit(g[1]);
  
  // h0 = sigma*v0;
  // h[1] = mu + phi[1]*(h0-mu) + sigma*sqrt(1-phi[1]^2)*v[1];
  h[1] = mu + sigma*v[1];
  
  for (t in 2:T) {
    g[t] = a + b*(g[t-1]-a) + c*sqrt(1-b^2)*w[t];
    phi[t] = inv_logit(g[t]);
    
    // h[t] = mu + phi[t]*(h[t-1]-mu) + sigma*sqrt(1-phi[t]^2)*v[t];
    h[t] = mu + phi[t]*(h[t-1]-mu) + sigma*sqrt(1-phi[t]^2)*v[t];
  }
}

model {
  mu_ret ~ normal(0, 1);
  phi ~ uniform(-1, 1);
  sigma ~ cauchy(0, 5);
  mu ~ cauchy(0, 10);
  
  // a ~ normal(3, 1);
  // b ~ beta(99, 1);
  // c ~ normal(0, 1);
  
  // v0 ~ std_normal();
  v ~ std_normal();
  
  // w0 ~ std_normal();
  w ~ std_normal();
  
  y ~ normal(mu_ret, exp(h / 2));
}
