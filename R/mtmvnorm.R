# Expectation and covariance matrix computation based on the Algorithm by Lee (1979), Lee (1983) and Leppard and Tallis (1989)
#
# Literatur:
# Amemiya (1973) : Regression Analysis When the Dependent Variable is Truncated Normal
# Amemiya (1974) : Multivariate Regression and Simultaneous Equations Models When the Dependent Variables Are Truncated Normal
# Lee (1979)     : On the first and second moments of the truncated multi-normal distribution and a simple estimator
# Lee (1983)     : The Determination of Moments of the Doubly Truncated Multivariate Tobit Model
# Leppard and Tallis (1989) : 

# Mean and Covariance of the truncated multivariate distribution (double truncation, general sigma, general mean)
#
# @param mean Mittelwertvektor (k x 1)
# @param sigma Kovarianzmatrix (k x k)
# @param lower untere Trunkierungspunkte (k x 1)
# @param upper obere Trunkierungspunkte (k x 1)
mtmvnorm <- function(mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), lower = rep(-Inf, length = length(mean)), upper = rep( Inf, length = length(mean)))
{
  N = length(mean)
  
  # Check input parameters
  if (NROW(sigma) != NCOL(sigma)) {
	stop("sigma must be a square matrix")
  }
	
  if (length(mean) != NROW(sigma)) {
	stop("mean and sigma have non-conforming size")
  }
  
  if (NCOL(lower) != NCOL(upper)) {
	stop("lower and upper have non-conforming size")
  }
	
  if (any(lower>=upper))
  {
	stop("lower must be smaller than or equal to upper (lower<=upper)")
  }

  # Truncated Mean
  TMEAN <- numeric(N)
  # Truncated Covariance matrix
  TVAR  <- matrix(NA, N, N)
		
  # Verschiebe die Integrationsgrenzen um -mean, damit der Mittelwert 0 wird 
  a     <- lower - mean
  b     <- upper - mean
  lower <- lower - mean
  upper <- upper - mean 
    	
  # eindimensionale Randdichte
  F_a <- numeric(N)
  F_b <- numeric(N)
    	
  # 1. Bestimme E[X_i]
  for (i in 1:N)
  {
  	sum = 0
  	for (q in 1:N)
  	{
  		F_a[q] = dtmvnorm.marginal(xn=a[q], n = q, mean=rep(0,N), sigma=sigma, lower=lower, upper=upper)
  		F_b[q] = dtmvnorm.marginal(xn=b[q], n = q, mean=rep(0,N), sigma=sigma, lower=lower, upper=upper)
  		sum = sum + sigma[i, q] * (F_a[q] - F_b[q])
  	}
  	TMEAN[i] = sum
  	# general mean case : TMEAN[i] = mean[i] + sum
  }
  
  # 2. Bestimme E[X_i, X_j]
  for (i in 1:N)
  {
  	for (j in 1:N)
  	{
  		sum = 0
  		for (q in 1:N)
  		{
			# Check if a[q] = -Inf or b[q]=+Inf, then F_a[q]=F_b[q]=0, but a[q] * F_a[q] = NaN  and b[q] * F_b[q] = NaN  
			F_a_q = a[q] * F_a[q]
			F_b_q = b[q] * F_b[q]
			if (is.infinite(a[q]))
			{
				F_a_q = 0 
			}
			if (is.infinite(b[q]))
			{
				F_b_q = 0 
			}
			
			sum  = sum + sigma[i,q] * sigma[j,q] * (sigma[q,q])^(-1) * (F_a_q - F_b_q)
  			sum2 = 0
  			for (s in 1:N)
  			{
  				if (s != q)
  				{
  					sum2 = sum2 + (sigma[j,s] - sigma[q,s] * sigma[j,q] * (sigma[q,q])^(-1)) * 
  							((dtmvnorm.marginal2(xq=a[q], xr=a[s], q=q, r=s,  mean=rep(0,N), sigma=sigma, lower=lower, upper=upper)  - 
                              dtmvnorm.marginal2(xq=b[q], xr=a[s], q=q, r=s,  mean=rep(0,N), sigma=sigma, lower=lower, upper=upper)) -
                             (dtmvnorm.marginal2(xq=a[q], xr=b[s], q=q, r=s,  mean=rep(0,N), sigma=sigma, lower=lower, upper=upper) - 
                              dtmvnorm.marginal2(xq=b[q], xr=b[s], q=q, r=s,  mean=rep(0,N), sigma=sigma, lower=lower, upper=upper)
                             ))
  				}
  			}
  			sum2 = sigma[i, q] * sum2
  			sum  = sum + sum2
  		}
  		TVAR[i, j] = sigma[i, j] + sum
  		#general mean case: TVAR[i, j] = mean[j] * TMEAN[i] + mean[i] * TMEAN[j] - mean[i] * mean[j] + sigma[i, j] + sum
  	}
  }
    	
  # 3. Bestimme Varianz Cov(X_i, X_j) = E[X_i, X_j] - E[X_i]*E[X_j] für (0, sigma)-case
  TVAR = TVAR - TMEAN %*% t(TMEAN)
    	
  # 4. Rückverschiebung um +mean für (mu, sigma)-case
  TMEAN = TMEAN + mean
		
  return(list(tmean=TMEAN, tvar=TVAR))
}

# Bestimmung von Erwartungswert und Kovarianzmatrix über numerische Integration und die eindimensionale Randdichte
# d.h. 
# E[X_i]       = \int_{a_i}^{b_i}{x_i * f(x_i) d{x_i}}
# Var[x_i]     = \int_{a_i}^{b_i}{(x_i-\mu_i)^2 * f(x_i) d{x_i}}
# Cov[x_i,x_j] = \int_{a_i}^{b_i}\int_{a_j}^{b_j}{(x_i-\mu_i)(x_j-\mu_j) * f(x_i,x_j) d{x_i}d{x_j}}
#
# Die Bestimmung von E[X_i] und Var[x_i]
# Die Bestimmung der Kovarianz Cov[x_i,x_j] benötigt die zweidimensionale Randdichte.
# 
#
# @param mean Mittelwertvektor (k x 1)
# @param sigma Kovarianzmatrix (k x k)
# @param lower, upper obere und untere Trunkierungspunkte (k x 1)
mtmvnorm.quadrature <- function(mean, sigma, lower, upper)
{
  k       = length(mean)
  
  # Bestimmung des Erwartungswerts/Varianz über numerische Integration
  expectation <- function(x, n=1)
  {
    x * dtmvnorm.marginal(x, n=n, mean=mean, sigma=sigma, lower=lower, upper=upper)
  }

  variance <- function(x, n=1)
  {
    (x - m.integration[n])^2 * dtmvnorm.marginal(x, n=n, mean=mean, sigma=sigma, lower=a, upper=b)
  }

  # Determine expectation from one-dimensional marginal distribution using integration
  # i=1..k
  m.integration<-numeric(k)
  for (i in 1:k)
  {
    m.integration[1]=integrate(expectation, lower[i], upper[i], n=i)$value 
  }
  
  # Determine variances from one-dimensional marginal distribution using integration
  # i=1..k
  v.integration<-numeric(k)
  for (i in 1:k)
  {
    v.integration[1]=integrate(variance, lower[i], upper[i], n=i)$value 
  }
  
  return(list(m=m.integration, v=v.integration))
}
