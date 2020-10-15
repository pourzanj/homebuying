ln_fsgt <- function(y, mu, sigma, lambda, q, p) {
  num <- abs(y - mu)^p
  den <- q*sigma^p * (lambda*sign(y-mu) + 1)^p
  
  log(p) - log(2*sigma) - (1/p)*log(q) - lbeta(1/p, q) - (1/p + q)*log(num/den + 1)
}

# Make sure log density matches R implementation
ln_fsgt(1, mu = 0.23, sigma = 0.5, lambda = -0.1, q = 5, p = 2)
dsgt(1, mu = 0.23, sigma = 0.5, lambda = -0.1, q = 5, p = 2, mean.cent = FALSE, var.adj = FALSE, log = TRUE)

# Check that star derivative is correct
star <- function(y, mu, sigma, lambda, q, p) {
  num <- abs(y - mu)^p
  den <- q*sigma^p * (lambda*sign(y-mu) + 1)^p
  num / den
}

d_star <- function(y, mu, sigma, lambda, q, p) {
  num <- p * abs(y - mu)^p
  den <- q*sigma^(p+1) * (lambda*sign(y-mu) + 1)^p
  -num / den
}

d_star(1, mu = 0.23, sigma = 0.5 + eps, lambda = -0.1, q = 5, p = 2)

eps <- 1e-5
(star(1, mu = 0.23, sigma = 0.5 + eps, lambda = -0.1, q = 5, p = 2) -
    star(1, mu = 0.23, sigma = 0.5, lambda = -0.1, q = 5, p = 2)) / eps

# Check rest of function
d_ln_fsgt <- function(y, mu, sigma, lambda, q, p) {
  # Get star
  num <- abs(y - mu)^p
  den <- q*sigma^p * (lambda*sign(y-mu) + 1)^p
  star <- num / den
  
  # Get d_star
  num <- p * abs(y - mu)^p
  den <- q*sigma^(p+1) * (lambda*sign(y-mu) + 1)^p
  d_star <- -num / den
  
  # Return final
  -(1/sigma) - (1/p + q) * (1 / (star+1)) * d_star
}

d_ln_fsgt(1, mu = 0.23, sigma = 0.5 + eps, lambda = -0.1, q = 5, p = 2)

eps <- 1e-5
(ln_fsgt(1, mu = 0.23, sigma = 0.5 + eps, lambda = -0.1, q = 5, p = 2) -
    ln_fsgt(1, mu = 0.23, sigma = 0.5, lambda = -0.1, q = 5, p = 2)) / eps
