functions {
  real eval_b_spline(real x, vector k, int order, int idx);
  real eval_b_spline(real x, vector k, int order, int idx) {
  
    real f;
    real w1 = 0.0;
    real w2 = 0.0; 
    
    if (order == 1) {
      f = (k[idx] <= x) && (x < k[idx + 1]);
    } else {
      
      if (k[idx] != k[idx + order - 1]) {
        w1 = (x - k[idx]) / (k[idx + order - 1] - k[idx]);
      }
      if (k[idx + 1] != k[idx + order]) {
        w2 = 1 - ((x - k[idx + 1]) / (k[idx + order] - k[idx + 1]));
      }
      
      f = w1 * eval_b_spline(x, k, order - 1, idx) +
          w2 * eval_b_spline(x, k, order - 1, idx + 1);
    }
    
    return(f);
  }
  
  vector eval_b_spline_basis(real x, vector k, int order) {

    int K = num_elements(k);
    int degree = order - 1;
    int N = K + degree - 1; // Number of basis functions
    
    vector[K+2*degree] k_extended = 
      append_row(append_row(rep_vector(k[1], degree), k), rep_vector(k[K], degree));
  
    vector[N-2] h;

    for(i in 2:(N-1)) {
      h[i-1] = eval_b_spline(x, k_extended, order, i);
    }
    
    return(h);
  }
}
data {
  int<lower=0> T;
  real y[T];
  real<lower=0> sigma1;
  
  int K;                 // num of knots
  vector[K] k;           // the sequence of knots
  int order;        
}
transformed data {
  int degree = order - 1;
  int H = K + degree - 1 - 2;          // Number of basis functions
  vector[T] ones = rep_vector(1.0, T);
}
parameters {
  real mu;
  
  real alpha0;
  vector[H] alpha;
  real beta;
}
transformed parameters {
  
  vector[T] h_hat;
  vector[T] h;        // log of volatility squared (variance)
  vector[T] sigma;    // volatility (SD of returns)
  vector[T-1] e;
  matrix[T, H+2] X;   // design matrix
  
  h_hat[1] = log(sigma1^2);
  h[1] = h_hat[1];
  sigma[1] = exp(0.5 * h[1]);
  
  for (t in 2:T) {
    e[t-1] = y[t-1] / sigma[t-1];
    
    X[t, 1] = 1.0;
    X[t, 2:(H+1)] = eval_b_spline_basis(e[t-1], k, order)';
    X[t, H+2] = h[t-1];

    h_hat[t] = X[t,] * append_row(append_row(alpha0, alpha), beta);
    h[t] = h_hat[t];
    sigma[t] = exp(0.5 * h[t]);
  }

}
model {
  y ~ normal(mu, sigma);
}
