# Verteilungsfunktion für die eindimensionale Randdichte f(xn) einer Truncated Multivariate Normal Distribution,
# vgl. Jack Cartinhour (1990) "One-dimensional marginal density functions of a truncated multivariate normal density function" für die Dichtefunktion
#
# @param xn Vektor der Länge l von Punkten, an dem die Verteilungsfunktion ausgewertet wird
# @param i Index  (1..n) dessen Randdichte berechnet werden soll
# @param mean (nx1) Mittelwertvektor
# @param sigma (nxn)-Kovarianzmatrix
# @param lower,upper Trunkierungsvektor lower <= x <= upper
ptmvnorm.marginal <- function(xn, n=1, mean=rep(0, nrow(sigma)), sigma=diag(length(mean)), lower=rep(-Inf, length = length(mean)), upper=rep( Inf, length = length(mean)))
{
   if (NROW(sigma) != NCOL(sigma)) {
     stop("sigma must be a square matrix")
   }

   if (length(mean) != NROW(sigma)) {
    stop("mean and sigma have non-conforming size")
   }
   
   if (det(as.matrix(sigma)) <= 0) {
    stop("sigma must be positive definite")
   }
	
   if (any(lower>=upper))
   {
     stop("lower must be smaller than or equal to upper (lower<=upper)")
   }
   
   # Anzahl der Dimensionen                
   k = length(mean)
   
   if (n < 1 || n > length(mean) || !is.numeric(n) || length(n) > 1 ||  !n %in% 1:length(mean))
   {
     stop("n must be a integer scalar in 1..length(mean)")
   }

   Fx     = numeric(length(xn))
   upper2 = upper
   alpha  = pmvnorm(lower = lower, upper = upper, mean = mean, sigma = sigma)
   for (i in 1:length(xn))
   {
    upper2[n] = xn[i]
    Fx[i]     = pmvnorm(lower=lower, upper=upper2, mean=mean, sigma=sigma)
   }
   return (Fx/alpha)
}
