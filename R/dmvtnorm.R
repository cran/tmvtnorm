# Dichtefunktion der Multivariaten Trunkierten Normalverteilung mit Trunkierungsvektor lower and upper
#
# vgl. Horrace (2005) "Some Results on the Multivariate Truncated Normal Distribution"
#
# @param x Argumentenvektor der Dichte der Länge n oder Matrix (T x n) mit T Beobachtungen
# @param mean  Mittelwertvektor der Länge n
# @param sigma Kovarianzmatrix (n x n)
# @param lower unterer Trunkierungsvektor (n x 1) mit lower <= x <= upper
# @param upper oberer Trunkierungsvektor (n x 1) mit lower <= x <= upper
dtmvnorm <- function(x, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), lower = rep(-Inf, length = length(mean)), upper = rep( Inf, length = length(mean)), log = FALSE)
{
  # Check inputs
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
  if (NCOL(lower) != NCOL(upper)) {
	stop("lower and upper have non-conforming size")
  }
  if (any(lower>=upper)) {
	stop("lower must be smaller than or equal to upper (lower<=upper)")
  }
  
  # Anzahl der Beobachtungen
  T = nrow(x)
  
  # check for each row if in support region
  insidesupportregion <- logical(T)
  for (i in 1:T)
  {
    insidesupportregion[i] = all(x[i,] >= lower & x[i,] <= upper & !any(is.infinite(x)))
  }
  
  # density value for points outside the support region
  dv = if (log) { -Inf } else { 0 }
  
  f <- ifelse(insidesupportregion, dmvnorm(x, mean=mean, sigma=sigma, log=log) / pmvnorm(lower=lower, upper=upper, mean=mean, sigma=sigma), dv)
  return(f)
}

#dtmvnorm(x=c(0,0), mean=c(0,0), sigma=diag(2))
#dmvnorm(x=c(0,0), mean=c(0,0), sigma=diag(2))
#dtmvnorm(x=matrix(c(0,0,1,1),2,2, byrow=TRUE), mean=c(0,0), sigma=diag(2))
#dtmvnorm(x=matrix(c(0,0,1,1),2,2, byrow=TRUE), mean=c(0,0), sigma=diag(2), lower=c(-1,-1), upper=c(0.5, 0.5))
dtmvnorm(x=matrix(c(0,0,1,1),2,2, byrow=TRUE), mean=c(0,0), sigma=diag(2), lower=c(-1,-1), upper=c(0.5, 0.5), log=TRUE)
dtmvnorm(as.matrix(seq(-1,2, by=0.1), ncol=1), mean=c(0.5), sigma=as.matrix(1.2^2), lower=0)


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

