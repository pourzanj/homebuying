data {
  int<lower=0> T;   // # time points (equally spaced)
  vector[T] ys;      // stock returns
  vector[T] yb;    // bond returns
}
parameters {
  real mus_ret;
  real mub_ret;
  corr_matrix[4] Rho;

  // SV
  real mus;                     // mean log volatility
  real mub;                     // mean log volatility
  real<lower=0, upper=1> phis;  // persistence of volatility
  real<lower=0, upper=1> phib;  // persistence of volatility
  real<lower=0> sigmas;         // white noise shock scale
  real<lower=0> sigmab;         // white noise shock scale
  
  vector[T] vs;
  vector[T] vb;
}
transformed parameters {
  vector[T] hs;
  vector[T] hb;
  vector[T] ss;
  vector[T] sb;
  
  hs[1] = sigmas * vs[1] / sqrt(1 - phis^2);
  hb[1] = sigmab * vb[1] / sqrt(1 - phib^2);
  
  for (t in 2:T) {
    hs[t] = phis * hs[t-1] + sigmas * vs[t];
    hb[t] = phib * hb[t-1] + sigmab * vb[t];
  }
  
  ss = exp((mus + hs) / 2);
  sb = exp((mub + hb) / 2);
}
model {
  // prior
  mus_ret ~ normal(0, 0.2);
  mub_ret ~ normal(0, 0.2);

  // likelihood
  for (t in 1:T) {
    matrix[4,4] S = quad_form_diag(Rho, [ss[t], sb[t], 1, 1]);
    [ys[t], yb[t], vs[t], vb[t]] ~ multi_normal([mus_ret, mub_ret, 0, 0], S);
  }
  
}
generated quantities {
  vector[4] x = multi_normal_rng([0, 0, 0, 0], Rho);
  real hs_rep = phis * hs[T] + sigmas * x[3];
  real hb_rep = phib * hb[T] + sigmab * x[4];
  real ss_rep = exp((mus + hs_rep) / 2);
  real sb_rep = exp((mub + hb_rep) / 2);
  real ys_rep = mus_ret + ss_rep*x[1];
  real yb_rep = mub_ret + sb_rep*x[2];
}
