# Sampling from Truncated multivariate Gaussian distribution using 
#
# a) Rejection sampling
# b) Gibbs sampler
# 
# Author: Stefan Wilhelm
#
# Literatur:
# (1) Kotecha (1999)
###############################################################################

# Erzeugt eine Matrix X (n x k) mit Zufallsrealisationen aus einer Trunkierten Multivariaten Normalverteilung mit k Dimensionen
# über Rejection Sampling oder Gibbs Sampler aus einer Multivariaten Normalverteilung
#
# @param n Anzahl der Realisationen
# @param mean Mittelwertvektor (k x 1) der Normalverteilung
# @param sigma Kovarianzmatrix (k x k) der Normalverteilung
# @param lower unterer Trunkierungsvektor (k x 1) mit lower <= x <= upper
# @param upper oberer Trunkierungsvektor (k x 1) mit lower <= x <= upper
rtmvnorm <- function(n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), lower = rep(-Inf, length = length(mean)), upper = rep( Inf, length = length(mean)), 
		algorithm=c("rejection", "gibbs"), ...)
{
  algorithm <- match.arg(algorithm)
  if (algorithm == "rejection")
  {
	retval <- rtmvnorm.rejection(n, mean, sigma, lower, upper, ...)  
  }
  else if (algorithm == "gibbs")
  {
	retval <- rtmvnorm.gibbs(n, mean, sigma, lower, upper, ...)  
  }
  return(retval)
}

# Erzeugt eine Matrix X (n x k) mit Zufallsrealisationen aus einer Trunkierten Multivariaten Normalverteilung mit k Dimensionen
# über Rejection Sampling aus einer Multivariaten Normalverteilung
#
# @param n Anzahl der Realisationen
# @param mean Mittelwertvektor (k x 1) der Normalverteilung
# @param sigma Kovarianzmatrix (k x k) der Normalverteilung
# @param lower unterer Trunkierungsvektor (k x 1) mit lower <= x <= upper
# @param upper oberer Trunkierungsvektor (k x 1) mit lower <= x <= upper
rtmvnorm.rejection <- function(n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), lower = rep(-Inf, length = length(mean)), upper = rep( Inf, length = length(mean)))
{
  if (missing(mean)) {
    mean <- rep(0, length = nrow(sigma))
  }
  if (missing(sigma)) {
    sigma <- diag(length(mean))
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
  
  if (any(lower>=upper))
  {
    stop("lower must be smaller than or equal to upper (lower<=upper)")
  }

  # k = Dimension
  k = length(mean)
  
  # Ergebnismatrix (n x k)
  Y = matrix(NA, n, k)
  
  # Anzahl der noch zu ziehenden Samples
  numSamples = n
  
  # Anzahl der akzeptierten Samples insgesamt
  numAcceptedSamplesTotal = 0
  
  # Akzeptanzrate alpha aus der Multivariaten Normalverteilung bestimmen
  alpha = pmvnorm(lower=lower, upper=upper, mean=mean, sigma=sigma)
  
  if (alpha <= 0.01) warning("Acceptance rate is very low and rejection sampling becomes inefficient. Consider using Gibbs sampling.")
  
  # Ziehe wiederholt aus der Multivariaten NV und schaue, wieviel Samples nach Trunkierung übrig bleiben
  while(numSamples > 0)
  {
    #cat("numSamplesAcceptedTotal=",numAcceptedSamplesTotal," numSamplesToDraw = ",numSamples,"\n")
    
    # Erzeuge N/alpha Samples aus einer multivariaten Normalverteilung: Wenn alpha zu niedrig ist, wird Rejection Sampling ineffizient und N/alpha zu groß. Dann nur N erzeugen
    nproposals = ifelse (numSamples/alpha > 1000000, numSamples, ceiling(max(numSamples/alpha,10)))
    X = rmvnorm(nproposals, mean=mean, sigma=sigma)
    
    # Bestimme den Anteil der Samples nach Trunkierung
    # Bug: ind= rowSums(lower <= X & X <= upper) == k
    # wesentlich schneller als : ind=apply(X, 1, function(x) all(x >= lower & x<=upper))
    ind <- logical(nproposals)
    for (i in 1:nproposals)
    {
      ind[i] = all(X[i,] >= lower & X[i,] <= upper)
    } 

    # akzeptierte Samples
    X.ind = X[ind,]
    
    # Anzahl der akzeptierten Samples in diesem Durchlauf
    if (k == 1) # im univariaten Fall
    {
      numAcceptedSamples = length(X.ind)
    }
    else
    {
      numAcceptedSamples = nrow(X.ind)
    }
    
    # Wenn nix akzeptiert wurde, dann weitermachen
    if (length(numAcceptedSamples) == 0 || numAcceptedSamples == 0) next
    
    #cat("numSamplesAccepted=",numAcceptedSamples," numSamplesToDraw = ",numSamples,"\n")
    if (numAcceptedSamples <= numSamples)
    {
      # Übernehmen aller akzeptierten Samples
      Y[(numAcceptedSamplesTotal+1):(numAcceptedSamplesTotal+numAcceptedSamples),] = X.ind
    }
    else
    {
      # Übernehmen nur der benötigten Samples
      Y[(numAcceptedSamplesTotal+1):(numAcceptedSamplesTotal+numSamples),] = X.ind[1:numSamples,]
    }
    
    # Anzahl der akzeptierten Samples insgesamt
    numAcceptedSamplesTotal = numAcceptedSamplesTotal + numAcceptedSamples
    
    # Anzahl der verbliebenden Samples
    numSamples = numSamples - numAcceptedSamples 
  }
  Y
}

# Gibbs Sampler for Truncated Univariate Normal Distribution
#
# Jayesh H. Kotecha and Petar M. Djuric (1999) : GIBBS SAMPLING APPROACH FOR GENERATION OF TRUNCATED MULTIVARIATE GAUSSIAN RANDOM VARIABLES
#
# @param n Anzahl der Realisationen
# @param mu
# @param sigma
# @param a unterer Trunkierungspunkt
# @param b oberer Trunkierungspunkt
rtnorm.gibbs <- function(n, mu=0, sigma=1, a=-Inf, b=Inf)
{
   # Draw from Gaussian Distribution	
   #x = rnorm(n, mu, sigma)
   
   # Draw from Uni(0,1)
   F = runif(n) 	
   
   #Phi(a) und Phi(b)
   Fa = pnorm(a, mu, sigma)
   Fb = pnorm(b, mu, sigma)
   
   # Truncated Normal Distribution, see equation (6), but F(x) ~ Uni(0,1), 
   # so we directly draw from Uni(0,1)...
   #y  =  mu + sigma * qnorm(pnorm(x)*(Fb - Fa) + Fa)
   y  =  mu + sigma * qnorm(F * (Fb - Fa) + Fa)	
   
   y
}

# Gibbs Sampler für Truncated Multivariate Normal Distribution
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
rtmvnorm.gibbs <- function(n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), lower = rep(-Inf, length = length(mean)), upper = rep( Inf, length = length(mean)), burn.in.samples = 0, start.value = NULL)
{
  if (missing(mean)) {
	mean <- rep(0, length = nrow(sigma))
  }
  if (missing(sigma)) {
	sigma <- diag(length(mean))
  }
	
  if (NCOL(lower) != NCOL(upper)) {
	stop("lower and upper have non-conforming size")
  }
	
  if (NROW(sigma) != NCOL(sigma)) 
  {
	stop("sigma must be a square matrix")
  }
	
  if (length(mean) != NROW(sigma)) 
  {
	stop("mean and sigma have non-conforming size")
  }
	
  if (any(lower>=upper))
  {
   stop("lower must be smaller than or equal to upper (lower<=upper)")
  }
	
  # dimension of X
  d <- length(mean)
  
  # number of burn-in samples
  S <- burn.in.samples
  if (!is.null(S))
  {
	if (S < 0) stop("number of burn-in samples must be non-negative")   
  }
	
  # Take start value given by user or determine from lower and upper	
  if (!is.null(start.value))
  {
    if (length(mean) != length(start.value)) stop("mean and start value have non-conforming size")
	if (any(start.value<lower || start.value>upper)) stop("start value is not inside support region") 
	x0 <- start.value 
  }
  else
  {
    # Start value from support region, may be lower or upper bound, if they are finite, 
	# if both are infinite, we take 0.
	x0  <- ifelse(is.finite(lower), lower, ifelse(is.finite(upper), upper, 0))
  }
  
  # actual number of samples to draw (including burn-in samples)
  n = n + S
  
  # Sample from univariate truncated normal distribution which is very fast.
  if (d == 1)
  {
    X = rtnorm.gibbs(n, mu=mean[1], sigma=sigma[1,1], a=lower[1], b=upper[1])
    return(X)
  }
      
  # Ergebnismatrix (n x k)
  X <- matrix(NA, n, d)
  
  # Draw from Uni(0,1)
  F = matrix(runif(n * d), n, d)
  
  # List of conditional standard deviations can be pre-calculated
  sd <- list(d)
  # List of t(Sigma_i) %*% solve(Sigma) term
  P  <- list(d)
  
  for(i in 1:d)
  {
    # Partitioning of Sigma
    Sigma    = sigma[-i,-i] # (d-1) x (d-1)
    sigma_ii = sigma[i,i]   # 1 x 1
    Sigma_i  = sigma[i,-i]  # (d-1) x 1
    
    P[[i]]   = t(Sigma_i) %*% solve(Sigma)
    sd[[i]]  = sqrt(sigma_ii - P[[i]] %*% Sigma_i)
  }
  
  x <- x0
  # For all samples
  for(j in 1:n)
  {
    # For all dimensions
    for(i in 1:d)
    {
      # Berechnung von bedingtem Erwartungswert und bedingter Varianz:
      # bedingte Varianz hängt nicht von x[-i] ab!
      mu_i     = mean[i]    + P[[i]] %*% (x[-i] - mean[-i])
      
      # Transformation
	  # TODO: evtl. nur ein Aufruf von pnorm...	
      F.tmp = pnorm(c(lower[i], upper[i]), mu_i, sd[[i]])
	  Fa    = F.tmp[1]
      Fb    = F.tmp[2]
	  x[i]  = mu_i + sd[[i]] * qnorm(F[j,i] * (Fb - Fa) + Fa)
    }
    X[j,] = x
  }
  # if number of burn-in samples has been given, discard burn-in samples
  if (S > 0)
  {
    return(X[-(1:S),])
  }
  else
  {
	return(X)
  }
}