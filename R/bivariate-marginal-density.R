# Berechnung der bivariaten marginal density einer truncated multivariate normal distribution 
# nach Tallis (1961) und Leppard and Tallis (1989)
#
# Literatur:
# Tallis (1961): "The Moment Generating Function of the Truncated Multi-normal Distribution"
# Leppard and Tallis (1989): "Evaluation of the Mean and Covariance of the Truncated Multinormal"

# (n-2) Integral, d.h. zweidimensionale Randdichte in Dimension q und r, da (n-2) Dimensionen rausintegriert werden.
# vgl. Tallis (1961), S.224 und Code Leppard (1989), S.550
#
# f(xq=b[q], xr=b[r])
#
# Achtung: Funktion ist momentan nicht vektorisiert!
#
# @param xq
# @param xr
# @param q Index für Dimension q
# @param r Index für Dimension r
# @param mean
# @param sigma
# @param lower
# @param upper
dtmvnorm.marginal2 <- function(xq, xr, q, r, mean=rep(0, nrow(sigma)), sigma=diag(length(mean)), lower=rep(-Inf, length = length(mean)), upper=rep( Inf, length = length(mean)))
{
  # Dimension
  n = nrow(sigma)
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
  
  if (n == 2)
  {
    if (xq < lower[q] | xq > upper[q] | xr < lower[r] | xr > upper[r] | is.infinite(xq) | is.infinite(xr))
    {
      return(0)
    }
    else
    {
      sigma2=sigma[c(q,r),c(q,r)]
      return (dmvnorm(x=c(xq, xr), mean=mean[c(q,r)], sigma=sigma2) / alpha)
    }
  }
  
  # Standardabweichungen bestimmen
  SD = sqrt(diag(sigma))
  
  # Normalisierte Bounds berechnen
  lower.normalised  = (lower - mean) / SD
  upper.normalised  = (upper - mean) / SD
  
  xq.normalised     = (xq - mean[q]) / SD[q]
  xr.normalised     = (xr - mean[r]) / SD[r]
  
  # Korrelationsmatrix R aus sigma berechnen (Matrix (n x n)) : R = D % sigma %*% D mit D Diagonalmatrix von sqrt(sigma)
  D       = matrix(0, n, n)
  diag(D) = sqrt(diag(sigma))^(-1)
  R       = D %*% sigma %*% D
  
  # Integrationsgrenzen
  AQR   = numeric(n-2)                    
  BQR   = numeric(n-2)
  
  #
  # Determine n-2 x n-2 correlation matrix RQR
  #
  RQR        = matrix(NA, n-2, n-2)
  RINV       = solve(R)
  WW         = matrix(NA, n-2, n-2)
  M1 = 0
  for (i in 1:n)
  {
    if (i != q && i != r)
    {
      M1 = M1 + 1
      M2 = 0
      for (j in 1:n)
      {
        if (j != q && j != r)
        {
          M2 = M2 + 1
          WW[M1, M2] = RINV[i,j]
        }
      }
    }
  }
  WW = solve(WW[1:(n-2),1:(n-2)])
  for(i in 1:(n-2))
  {
    for(j in 1:(n-2))
    {
       RQR[i, j] = WW[i, j] / sqrt(WW[i,i] * WW[j,j])
    }
  }
  
  #
  # Determine bounds of integration vector aqr
  #
  M2 = 0
  for (i in 1:n)
  {
    if (i != q && i != r)
    {
      M2 = M2 + 1
      BSQR = (R[q, i] - R[q, r] * R[r, i]) / (1 - R[q, r]^2)
      BSRQ = (R[r, i] - R[q, r] * R[q, i]) / (1 - R[q, r]^2)
      RSRQ = (1 - R[i, q]^2) * (1 - R[q, r]^2)
      RSRQ = (R[i, r] - R[i, q] * R[q, r]) / sqrt(RSRQ)
      
      # Untere Integrationsgrenze
      AQR[M2] = (lower.normalised[i] - BSQR * xq.normalised - BSRQ * xr.normalised) / sqrt((1 - R[i, q]^2)*(1 - RSRQ^2))
      if (is.nan(AQR[M2]))
      {
        AQR[M2]=-Inf
      }
      # obere Integrationsgrenze
      BQR[M2] = (upper.normalised[i] - BSQR * xq.normalised - BSRQ * xr.normalised) / sqrt((1 - R[i, q]^2)*(1 - RSRQ^2))
      if (is.nan(BQR[M2]))
      {
        BQR[M2]=Inf
      }
    }
  }
  
  # Korrelationsmatrix zwischen r und q
  R2=matrix(c(1,      R[q,r], 
              R[q,r], 1), 2, 2)
              
  sigma2=sigma[c(q,r),c(q,r)]            
              
              
  if (xq < lower[q] | xq > upper[q] | xr < lower[r] | xr > upper[r] | is.infinite(xq) | is.infinite(xr))
  {
    return(0)
  }
  else
  {
    dmvnorm(x=c(xq, xr), mean=mean[c(q,r)], sigma=sigma2) * pmvnorm(lower=AQR, upper=BQR, sigma=RQR) / alpha
  }
}


