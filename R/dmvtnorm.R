# Dichtefunktion der Multivariaten Trunkierten Normalverteilung mit Trunkierungsvektor lower and upper
#
# vgl. Horrace (2005) "Some Results on the Multivariate Truncated Normal Distribution"
#
# @param x Argumentenvektor der Dichte der Länge n
# @param mean  Mittelwertvektor der Länge n
# @param sigma Kovarianzmatrix (n x n)
# @param lower unterer Trunkierungsvektor (k x 1) mit lower <= x <= upper
# @param upper oberer Trunkierungsvektor (k x 1) mit lower <= x <= upper
dtmvnorm <- function(x, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), lower = rep(-Inf, length = length(mean)), upper = rep( Inf, length = length(mean)))
{
  require(mvtnorm)
  
  if (is.vector(x))
  {
    x <- matrix(x, ncol = length(x))
  }
  
  if (missing(mean)) {
    mean <- rep(0, length = ncol(x))
  }
  if (missing(sigma)) {
    sigma <- diag(ncol(x))
  }
  
  if (NCOL(x) != NCOL(sigma)) {
    stop("x and sigma have non-conforming size")
  }

  if (NROW(sigma) != NCOL(sigma)) {
    stop("sigma must be a square matrix")
  }

  if (length(mean) != NROW(sigma)) {
    stop("mean and sigma have non-conforming size")
  }
  
  # lower <= x <= upper
  if (any(x < lower | x > upper)) return(NA) 
    
  f <- dmvnorm(x, mean=mean, sigma=sigma) / pmvnorm(lower=lower, upper=upper, mean=mean, sigma=sigma)
  return(f)
}

#dtmvnorm(x=c(0,0), mean=c(0,0), sigma=diag(2))
#dmvnorm(x=c(0,0), mean=c(0,0), sigma=diag(2))


# Verteilungsfunktion der truncated multivariate normal distribution
#
# @param lower unterer Trunkierungsvektor (k x 1) mit lower <= x <= upper
# @param upper oberer Trunkierungsvektor (k x 1) mit lower <= x <= upper
ptmvnorm <- function(lowerx, upperx, mean=rep(0, length(lowerx)), sigma, lower = rep(-Inf, length = length(mean)), upper = rep( Inf, length = length(mean)), maxpts = 25000, abseps = 0.001, releps = 0)
{
  if (missing(mean)) {
		mean <- rep(0, length = ncol(lowerx))
  }
  if (missing(sigma)) {
		sigma <- diag(ncol(lowerx))
  }
	
  if (NCOL(lower) != NCOL(upper)) {
		stop("lower and upper have non-conforming size")
  }
	
  if (NROW(sigma) != NCOL(sigma)) {
		stop("sigma must be a square matrix")
  }
	
  if (length(mean) != NROW(sigma)) {
		stop("mean and sigma have non-conforming size")
  }
	
  if (any(lowerx>=upperx))
  {
	stop("lowerx must be smaller than or equal to upperx (lowerx<=upperx)")
  }	
  if (any(lower>=upper))
  {
	stop("lower must be smaller than or equal to upper (lower<=upper)")
  }
  
  f <- pmvnorm(lower=lowerx, upper=upperx, mean=mean, sigma=sigma, maxpts = maxpts, abseps = abseps, releps = releps) / 
	   pmvnorm(lower=lower, upper=upper, mean=mean, sigma=sigma, maxpts = maxpts, abseps = abseps, releps = releps)
  return(f)
}

