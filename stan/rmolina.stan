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
      
      real eps = k*(neg+pos);
      
      ret = ret + eps;
      if(fabs(eps) <= 1e-7) break;
      if(n == 100) print("Max iter reached", " eps: ", eps);
    }
    
    if(a > min([y, y0])) print("a > min(y, y0)");
    if(b < max([y, y0])) print("b < max(y, y0)");
    if(a > b) print("a > b");
    
    
    // print("ret: ", ret, " [a,y,b] = ", [a, y, b]);
    return log(ret);
  }
}
data {
  real y0;
  real m;
  real<lower=0> s;
}
parameters {
  real<upper=y0> a;
  real<lower=y0> b;
  real<lower=a,upper=b> y;
}
transformed parameters {
  
}
model {
  [a, b, y]' ~ molina(y0, m, s);
}
generated quantities {
  real r = y - y0;
  real r_max = a - y0;
  real r_min = b - y0;
}