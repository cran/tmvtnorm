# Computation of the bivariate marginal density F_{q,r}(x_q, x_r) (q != r)
# of truncated multivariate normal distribution 
# following the works of Tallis (1961), Leppard and Tallis (1989)
#
# References:
# Tallis (1961): 
#   "The Moment Generating Function of the Truncated Multi-normal Distribution"
# Leppard and Tallis (1989): 
#   "Evaluation of the Mean and Covariance of the Truncated Multinormal"
# Manjunath B G and Stefan Wilhelm (2009): 
#   "Moments Calculation for the Doubly Truncated Multivariate Normal Distribution"
#
# (n-2) Integral, d.h. zweidimensionale Randdichte in Dimension q und r, 
# da (n-2) Dimensionen rausintegriert werden.
# vgl. Tallis (1961), S.224 und Code Leppard (1989), S.550
#
# f(xq=b[q], xr=b[r])
#
# Attention: Function is not vectorized at the moment!
#
# @param xq
# @param xr
# @param q index for dimension q
# @param r Index für Dimension r
# @param mean
# @param sigma
# @param lower
# @param upper
# @param log=FALSE
dtmvnorm.marginal2 <- function(xq, xr, q, r, mean=rep(0, nrow(sigma)), 
  sigma=diag(length(mean)), lower=rep(-Inf, length = length(mean)), 
	upper=rep( Inf, length = length(mean)), log=FALSE) {
	
  # dimensionality
  n <- nrow(sigma)
  
  # input checks
  if (n < 2) stop("Dimension n must be >= 2!")
  
  if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps))) {
    stop("sigma must be a symmetric matrix")
  }
  
  if (length(mean) != NROW(sigma)) {
    stop("mean and sigma have non-conforming size")
  }
  
  if (!(q %in% 1:n && r %in% 1:n)) {
    stop("Indexes q and r must be integers in 1:n")
  }
  
  if (q == r) {
    stop("Index q must be different than r!")
  }
  
  # Skalierungsfaktor der gestutzten Dichte (Anteil nach Trunkierung)
  alpha = pmvnorm(lower=lower, upper=upper, mean=mean, sigma=sigma)
  
  if (n == 2) {
    if (xq < lower[q] | xq > upper[q] | xr < lower[r] | xr > upper[r] | is.infinite(xq) | is.infinite(xr)) {
	  density = 0
    } else {
      sigma2 <- sigma[c(q,r),c(q,r)]
      density <- dmvnorm(x=c(xq, xr), mean=mean[c(q,r)], sigma=sigma2) / alpha
    }
	  if (log == TRUE) {
	    return(log(density))
	  } else {
	    return(density)
	  }
  }
  
  # standard deviation for normalisation
  SD <- sqrt(diag(sigma))
  
  # normalised bounds
  lower.normalised  <- (lower - mean) / SD
  upper.normalised  <- (upper - mean) / SD
  
  xq.normalised     <- (xq - mean[q]) / SD[q]
  xr.normalised     <- (xr - mean[r]) / SD[r]
  
  # Computing correlation matrix R from sigma (matrix (n x n)): 
  # R = D % sigma %*% D with diagonal matrix D as sqrt(sigma)
  D  <- matrix(0, n, n)
  diag(D) <- sqrt(diag(sigma))^(-1)
  R <- D %*% sigma %*% D
  
  # lower and upper integration bounds
  AQR <- numeric(n-2)                    
  BQR <- numeric(n-2)
  
  #
  # Determine n-2 x n-2 correlation matrix RQR
  #
  RQR <- matrix(NA, n-2, n-2)
  RINV <- solve(R)
  WW <- matrix(NA, n-2, n-2)
  M1 <- 0
  for (i in 1:n) {
    if (i != q && i != r) {
      M1 <- M1 + 1
      M2 <- 0
      for (j in 1:n) {
        if (j != q && j != r) {
          M2 <- M2 + 1
          WW[M1, M2] <- RINV[i,j]
        }
      }
    }
  }
  WW <- solve(WW[1:(n-2),1:(n-2)])
  for(i in 1:(n-2)) {
    for(j in 1:(n-2)) {
       RQR[i, j] <- WW[i, j] / sqrt(WW[i,i] * WW[j,j])
    }
  }
  
  #
  # Determine bounds of integration vector aqr
  #
  M2 <- 0  # counter = 1..(n-2)
  for (i in 1:n) {
    if (i != q && i != r) {
      M2 <- M2 + 1
      BSQR <- (R[q, i] - R[q, r] * R[r, i]) / (1 - R[q, r]^2)    
      BSRQ <- (R[r, i] - R[q, r] * R[q, i]) / (1 - R[q, r]^2)    
      RSRQ <- (1 - R[i, q]^2) * (1 - R[q, r]^2)
      RSRQ <- (R[i, r] - R[i, q] * R[q, r]) / sqrt(RSRQ)         # partial correlation coefficient R[r,i] given q
      
      # lower integration bound
      AQR[M2] <- (lower.normalised[i] - BSQR * xq.normalised - BSRQ * xr.normalised) / sqrt((1 - R[i, q]^2) * (1 - RSRQ^2))
      if (is.nan(AQR[M2])) {
        AQR[M2] <- -Inf
      }
      # upper integration bound
      BQR[M2] <- (upper.normalised[i] - BSQR * xq.normalised - BSRQ * xr.normalised) / sqrt((1 - R[i, q]^2) * (1 - RSRQ^2))
      if (is.nan(BQR[M2])) {
        BQR[M2] <- Inf
      }
    }
  }
  
  # Correlation matrix for r and q
  R2 <- matrix(c(    1,      R[q,r], 
                R[q,r],          1), 2, 2)
              
  sigma2 <- sigma[c(q,r),c(q,r)]            
              
  if (xq < lower[q] | 
      xq > upper[q] | 
      xr < lower[r] | 
      xr > upper[r] | is.infinite(xq) | is.infinite(xr)) {
    density = 0
  }
  else {
	  density <- dmvnorm(x=c(xq, xr), mean=mean[c(q,r)], sigma=sigma2) * pmvnorm(lower=AQR, upper=BQR, sigma=RQR) / alpha
  }
  if (log == TRUE) {
	  return(log(density))
  } else {
	  return(density)
  }
}


