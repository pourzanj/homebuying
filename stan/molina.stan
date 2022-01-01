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
      if(n == 1000 || is_nan(ret)) print("Max iter reached", " iter: ", n, " k: ", k, " eps: ", eps, " m: ", m, " s: ", s, "a, b, c: " , x, " ret: ", ret, " neg: ", neg, " pos: ", pos);
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
parameters {
  // vector[T] m;
  vector[T] h;
}
transformed parameters {
  vector<lower=0>[T] s;
  
  s = exp(h/2);
}
model {
  // m ~ normal(0.0, 1.0);
  
  for(t in 1:T) {
   [a[t], b[t], y[t]]' ~ molina(y0[t], 0.0, s[t]); 
  }
}
generated quantities {
}