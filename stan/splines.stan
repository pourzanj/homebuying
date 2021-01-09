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
  int N;             // number of data points
  real x[N];
  real y[N];
  
  int K;            // num of knots
  vector[K] k;  // the sequence of knots
  int order;        
}
transformed data {
  int degree = order - 1;
  int H = K + degree - 1 - 2; // Number of basis functions
  vector[N] ones = rep_vector(1.0, N);
}
parameters {
  real a0; 
  vector[H] a;
}
transformed parameters {
  vector[N] h;
  matrix[N, H+1] X;
  X[, 1] = ones;
  
  for(i in 1:N) {
    X[i, 2:(H+1)] = eval_b_spline_basis(x[i], k, order)';
    h[i] = X[i,] * append_row(a0, a);
  }
}
model {
  // Priors
  a0 ~ normal(0, 1);
  a ~ normal(0, 1);

  //Likelihood
  y ~ normal(h, 1e-1);
}
