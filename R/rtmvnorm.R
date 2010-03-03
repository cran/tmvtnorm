# Sampling from Truncated multivariate Gaussian distribution using 
#
# a) Rejection sampling
# b) Gibbs sampler
# 
# Author: Stefan Wilhelm
#
# Literatur:
# (1) Kotecha (1999)
# (2) Geweke (1991) : 
# "Effcient simulation from the multivariate normal and Student-t distributions 
# subject to linear constraints and the evaluation of constraint probabilities"
###############################################################################

# Uses partly checks as in mvtnorm:::checkmvArgs!
checkTmvArgs <- function(mean, sigma, lower, upper)
{
	if (is.null(lower) || any(is.na(lower))) 
		stop(sQuote("lower"), " not specified or contains NA")
	if (is.null(upper) || any(is.na(upper))) 
		stop(sQuote("upper"), " not specified or contains NA")
	if (!is.numeric(mean) || !is.vector(mean)) 
		stop(sQuote("mean"), " is not a numeric vector")
	if (is.null(sigma) || any(is.na(sigma))) 
		stop(sQuote("sigma"), " not specified or contains NA")
	
	if (!is.matrix(sigma)) {
		sigma = as.matrix(sigma)
	}
	
	if (NCOL(lower) != NCOL(upper)) {
		stop("lower and upper have non-conforming size")
	}
	
	if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps))) {
		stop("sigma must be a symmetric matrix")
	}
	
	if (NROW(sigma) != NCOL(sigma)) {
		stop("sigma must be a square matrix")
	}
	
	if (length(mean) != NROW(sigma)) {
		stop("mean and sigma have non-conforming size")
	}
	
	if (det(sigma) <= 0) {
		stop("sigma must be positive definite")
	}
	
	if (length(lower) != length(mean) || length(upper) != length(mean)) {
		stop("mean, lower and upper must have the same length")
	}
	
	if (any(lower>=upper)) {
		stop("lower must be smaller than or equal to upper (lower<=upper)")
	}
	
	# checked arguments
	cargs <- list(mean=mean, sigma=sigma, lower=lower, upper=upper)
	return(cargs)
}

# Erzeugt eine Matrix X (n x k) mit Zufallsrealisationen aus einer Trunkierten Multivariaten Normalverteilung mit k Dimensionen
# ¸ber Rejection Sampling oder Gibbs Sampler aus einer Multivariaten Normalverteilung
#
# @param n Anzahl der Realisationen
# @param mean Mittelwertvektor (k x 1) der Normalverteilung
# @param sigma Kovarianzmatrix (k x k) der Normalverteilung
# @param lower unterer Trunkierungsvektor (k x 1) mit lower <= x <= upper
# @param upper oberer Trunkierungsvektor (k x 1) mit lower <= x <= upper
# @param algorithm c("rejection", "gibbs", "gibbsR")
rtmvnorm <- function(n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), lower = rep(-Inf, length = length(mean)), upper = rep( Inf, length = length(mean)), 
		algorithm=c("rejection", "gibbs", "gibbsR"), ...)
{
  algorithm <- match.arg(algorithm)
  
  # check of standard tmvtnorm arguments
  cargs <- checkTmvArgs(mean, sigma, lower, upper)
  mean  <- cargs$mean
  sigma <- cargs$sigma
  lower <- cargs$lower
  upper <- cargs$upper
  
  # check of additional arguments
  if (n < 1 || !is.numeric(n) || n != as.integer(n) || length(n) > 1) {
	  stop("n must be a integer scalar > 0")
  }
  
  if (algorithm == "rejection") {
	retval <- rtmvnorm.rejection(n, mean, sigma, lower, upper, ...)  
  } else if (algorithm == "gibbs") {
	retval <- rtmvnorm.gibbs.Fortran(n, mean, sigma, lower, upper, ...)  
  } else if (algorithm == "gibbsR") {
	retval <- rtmvnorm.gibbs(n, mean, sigma, lower, upper, ...)  
  }
  return(retval)
}

# Erzeugt eine Matrix X (n x k) mit Zufallsrealisationen aus einer Trunkierten Multivariaten Normalverteilung mit k Dimensionen
# ¸ber Rejection Sampling aus einer Multivariaten Normalverteilung
#
# @param n Anzahl der Realisationen
# @param mean Mittelwertvektor (k x 1) der Normalverteilung
# @param sigma Kovarianzmatrix (k x k) der Normalverteilung
# @param lower unterer Trunkierungsvektor (k x 1) mit lower <= x <= upper
# @param upper oberer Trunkierungsvektor (k x 1) mit lower <= x <= upper
rtmvnorm.rejection <- function(n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), lower = rep(-Inf, length = length(mean)), upper = rep( Inf, length = length(mean)))
{
  # No check of input parameters, checks are done in rtmvnorm()!
  
  # k = Dimension
  k <- length(mean)
  
  # Ergebnismatrix (n x k)
  Y <- matrix(NA, n, k)
  
  # Anzahl der noch zu ziehenden Samples
  numSamples <- n
  
  # Anzahl der akzeptierten Samples insgesamt
  numAcceptedSamplesTotal <- 0
  
  # Akzeptanzrate alpha aus der Multivariaten Normalverteilung bestimmen
  alpha <- pmvnorm(lower=lower, upper=upper, mean=mean, sigma=sigma)
  
  if (alpha <= 0.01) warning("Acceptance rate is very low and rejection sampling becomes inefficient. Consider using Gibbs sampling.")
  
  # Ziehe wiederholt aus der Multivariaten NV und schaue, wieviel Samples nach Trunkierung ¸brig bleiben
  while(numSamples > 0)
  {
    # Erzeuge N/alpha Samples aus einer multivariaten Normalverteilung: Wenn alpha zu niedrig ist, wird Rejection Sampling ineffizient und N/alpha zu groﬂ. Dann nur N erzeugen
    nproposals <- ifelse (numSamples/alpha > 1000000, numSamples, ceiling(max(numSamples/alpha,10)))
    X <- rmvnorm(nproposals, mean=mean, sigma=sigma)
    
    # Bestimme den Anteil der Samples nach Trunkierung
    # Bug: ind= rowSums(lower <= X & X <= upper) == k
    # wesentlich schneller als : ind=apply(X, 1, function(x) all(x >= lower & x<=upper))
    ind <- logical(nproposals)
    for (i in 1:nproposals)
    {
      ind[i] = all(X[i,] >= lower & X[i,] <= upper)
    } 

    # Anzahl der akzeptierten Samples in diesem Durchlauf
    numAcceptedSamples <- length(ind[ind==TRUE])
        
    # Wenn nix akzeptiert wurde, dann weitermachen
    if (length(numAcceptedSamples) == 0 || numAcceptedSamples == 0) next
    
    #cat("numSamplesAccepted=",numAcceptedSamples," numSamplesToDraw = ",numSamples,"\n")
	numNeededSamples <- min(numAcceptedSamples, numSamples)
	Y[(numAcceptedSamplesTotal+1):(numAcceptedSamplesTotal+numNeededSamples),] <- X[which(ind)[1:numNeededSamples],]
        
    # Anzahl der akzeptierten Samples insgesamt
    numAcceptedSamplesTotal <- numAcceptedSamplesTotal + numAcceptedSamples
    
    # Anzahl der verbliebenden Samples
    numSamples <- numSamples - numAcceptedSamples 
  }
  Y
}

# Gibbs Sampler for Truncated Univariate Normal Distribution
#
# Jayesh H. Kotecha and Petar M. Djuric (1999) : GIBBS SAMPLING APPROACH FOR GENERATION OF TRUNCATED MULTIVARIATE GAUSSIAN RANDOM VARIABLES
#
# Im univariaten Fall sind die erzeugten Samples unabh‰ngig, deswegen gibt es hier keine Chain im eigentlichen Sinn und auch keinen Startwert/Burn-in/Thinning. 
#
# @param n Anzahl der Realisationen
# @param mu
# @param sigma
# @param a unterer Trunkierungspunkt
# @param b oberer Trunkierungspunkt
rtnorm.gibbs <- function(n, mu=0, sigma=1, a=-Inf, b=Inf)
{
   # Draw from Gaussian Distribution	
   #x <- rnorm(n, mu, sigma)
   
   # Draw from Uni(0,1)
   F <- runif(n) 	
   
   #Phi(a) und Phi(b)
   Fa <- pnorm(a, mu, sigma)
   Fb <- pnorm(b, mu, sigma)
   
   # Truncated Normal Distribution, see equation (6), but F(x) ~ Uni(0,1), 
   # so we directly draw from Uni(0,1)...
   #y  <-  mu + sigma * qnorm(pnorm(x)*(Fb - Fa) + Fa)
   y  <-  mu + sigma * qnorm(F * (Fb - Fa) + Fa)	
   
   y
}

# Gibbs Sampler f¸r Truncated Multivariate Normal Distribution
#
# Jayesh H. Kotecha and Petar M. Djuric (1999) : GIBBS SAMPLING APPROACH FOR GENERATION OF TRUNCATED MULTIVARIATE GAUSSIAN RANDOM VARIABLES
#
#
# @param n Anzahl der Realisationen
# @param mean Mittelwertvektor (k x 1) der Normalverteilung
# @param sigma Kovarianzmatrix (k x k) der Normalverteilung
# @param lower unterer Trunkierungsvektor (k x 1) mit lower <= x <= upper
# @param upper oberer Trunkierungsvektor (k x 1) mit lower <= x <= upper
# @param burn.in number of burn-in samples to be discarded
# @param start start value for Gibbs sampling
# @param thinning
rtmvnorm.gibbs <- function(n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), lower = rep(-Inf, length = length(mean)), upper = rep( Inf, length = length(mean)), 
		burn.in.samples = 0, start.value = NULL, thinning = 1)
{
  # We check only additional arguments like "burn.in.samples", "start.value" and "thinning"
  
  if (thinning < 1 || !is.numeric(thinning) || length(thinning) > 1) {
	  stop("thinning must be a integer scalar > 0")
  }
	
  # dimension of X
  d <- length(mean)
  
  # number of burn-in samples
  S <- burn.in.samples
  if (!is.null(S)) {
	if (S < 0) stop("number of burn-in samples must be non-negative")   
  }
	
  # Take start value given by user or determine from lower and upper	
  if (!is.null(start.value)) {
    if (length(mean) != length(start.value)) stop("mean and start value have non-conforming size")
	if (any(start.value<lower || start.value>upper)) stop("start value is not inside support region") 
	x0 <- start.value 
  } else {
    # Start value from support region, may be lower or upper bound, if they are finite, 
	# if both are infinite, we take 0.
	x0  <- ifelse(is.finite(lower), lower, ifelse(is.finite(upper), upper, 0))
  }
  
  # Sample from univariate truncated normal distribution which is very fast.
  if (d == 1)
  {
    X <- rtnorm.gibbs(n, mu=mean[1], sigma=sigma[1,1], a=lower[1], b=upper[1])
    return(X)
  }
      
  # Ergebnismatrix (n x k)
  X <- matrix(NA, n, d)
  
  # Draw from Uni(0,1)
  U <- runif((S + n*thinning) * d)
  l <- 1
  
  # List of conditional standard deviations can be pre-calculated
  sd <- list(d)
  # List of t(Sigma_i) %*% solve(Sigma) term
  P  <- list(d)
  
  for(i in 1:d)
  {
    # Partitioning of Sigma
    Sigma    <- sigma[-i,-i] # (d-1) x (d-1)
    sigma_ii <- sigma[i,i]   # 1 x 1
    Sigma_i  <- sigma[i,-i]  # 1 x (d-1)
    
    P[[i]]   <- t(Sigma_i) %*% solve(Sigma)  # (1 x (d-1)) * ((d-1) x (d-1)) =  (1 x (d-1))
    sd[[i]]  <- sqrt(sigma_ii - P[[i]] %*% Sigma_i)  # (1 x (d-1)) * ((d-1) x 1)
  }
  
  x <- x0
  
  # Runn chain from index (1 - #burn-in-samples):(n*thinning) and only record samples from j >= 1
  # which discards the burn-in-samples
  for (j in (1-S):(n*thinning))
  {
    # For all dimensions
    for(i in 1:d)
    {
      # Berechnung von bedingtem Erwartungswert und bedingter Varianz:
      # bedingte Varianz h‰ngt nicht von x[-i] ab!
      mu_i  <- mean[i]    + P[[i]] %*% (x[-i] - mean[-i])
      
      # Transformation
	  F.tmp <- pnorm(c(lower[i], upper[i]), mu_i, sd[[i]])
	  Fa    <- F.tmp[1]
      Fb    <- F.tmp[2]
	  x[i]  <- mu_i + sd[[i]] * qnorm(U[l] * (Fb - Fa) + Fa)
      l     <- l + 1
    }
	
	if (j > 0) {
	  if (thinning == 1) {
	    # no thinning, take all samples	except for burn-in-period
	    X[j,] <- x
	  }
	  else if (j %% thinning == 0){
	    X[j %/% thinning,] <- x	
	  }
    }
  }
  return(X)
}

# Versuch, die Gibbs Sampler Methode schneller zu machen mittels kompiliertem Fortran-Code
rtmvnorm.gibbs.Fortran <- function(n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), lower = rep(-Inf, length = length(mean)), upper = rep( Inf, length = length(mean)), 
		burn.in.samples = 0, start.value = NULL, thinning = 1)
{
  # No checks of input arguments, checks are done in rtmvnorm()
  
  # dimension of X
  d <- length(mean)
  
  # number of burn-in samples
  S <- burn.in.samples
  if (!is.null(S)) {
	if (S < 0) stop("number of burn-in samples must be non-negative")   
  }
	
  # Take start value given by user or determine from lower and upper	
  if (!is.null(start.value)) {
    if (length(mean) != length(start.value)) stop("mean and start value have non-conforming size")
	if (any(start.value<lower || start.value>upper)) stop("start value is not inside support region") 
	x0 <- start.value 
  } else {
    # Start value from support region, may be lower or upper bound, if they are finite, 
	# if both are infinite, we take 0.
	x0  <- ifelse(is.finite(lower), lower, ifelse(is.finite(upper), upper, 0))
  }
  
  # Sample from univariate truncated normal distribution which is very fast.
  if (d == 1) {
    X <- rtnorm.gibbs(n, mu=mean[1], sigma=sigma[1,1], a=lower[1], b=upper[1])
    return(X)
  }
      
  # Ergebnismatrix (n x d)
  X <- matrix(0, n, d)
  
  # Call to Fortran subroutine
  ret <- .Fortran("rtmvnormgibbs",
                              n     = as.integer(n),
                              d     = as.integer(d),
                              mean  = as.double(mean),
                              sigma = as.double(sigma),
                              lower = as.double(lower), 
                              upper = as.double(upper),
                              x0    = as.double(x0),
							  burnin   = as.integer(burn.in.samples),
							  thinning = as.integer(thinning),
                              X     = as.double(X), 
                              NAOK=TRUE, PACKAGE="tmvtnorm")
  X <- matrix(ret$X, ncol=d, byrow=TRUE)
  return(X)
}
