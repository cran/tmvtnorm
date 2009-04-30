% --- Source file: tmvtnorm.Rd ---
\name{rtmvnorm}
\alias{rtmvnorm}
\title{Sampling Random Numbers From Truncated Multivariate Normal Distribution}
\description{
  This function generates random numbers
  from the truncated multivariate normal
  distribution with mean equal to \code{mean} and covariance matrix
  \code{sigma}, lower and upper truncation points \code{lower} and \code{upper}
  with either rejection sampling or Gibbs sampling.
}
\usage{
rtmvnorm(n, mean = rep(0, nrow(sigma)), 
  sigma = diag(length(mean)), 
  lower=rep(-Inf, length = length(mean)), 
  upper=rep( Inf, length = length(mean)),
  algorithm=c("rejection", "gibbs"),
  ...)
}

\arguments{
 \item{n}{Number of observations.}
 \item{mean}{Mean vector, default is \code{rep(0, length = ncol(x))}.}
 \item{sigma}{Covariance matrix, default is \code{diag(ncol(x))}.}
 \item{lower}{Vector of lower truncation points,\\ 
              default is \code{rep(-Inf, length = length(mean))}.}
 \item{upper}{Vector of upper truncation points,\\ 
              default is \code{rep( Inf, length = length(mean))}.}
 \item{algorithm}{Method used, possible methods are rejection sampling ("rejection", default) and Gibbs sampler ("gibbs").}
 \item{...}{additional parameters for Gibbs sampling, given to \code{rmvtnorm.gibbs}, such as burn.in.samples and start.value}               
}

\details{

The generation of random numbers from a truncated multivariate normal distribution is done 
using either rejection sampling or Gibbs sampler.

Rejection sampling is done from the standard multivariate normal distribution. 
So we use the function \code{\link[mvtnorm]{rmvnorm}} of the mvtnorm package.
In order to speed up the generation of N samples from the truncated distribution, 
we first calculate the acceptance rate alpha from the truncation points and then generate N/alpha samples iteratively 
until we have got N samples. This typically does not take more then 2-3 iterations.
Rejection sampling may be very inefficient when the support region is small (i.e. in higher dimensions) 
which results in very low acceptance rates alpha. In this case the Gibbs sampler is preferable.

The Gibbs Sampler samples from univariate conditional distributions, 
so all samples can be accepted except for a burn-in period. 
The number of burn-in samples to be discarded can be specified, as well as a start value.
If no start value is given we determine a start value from the support region
using either lower bound or upper bound if they are finite, or 0 otherwise.
}
 
\author{Stefan Wilhelm <Stefan.Wilhelm@financial.com>, Manjunath B G <bgmanjunath@gmail.com>}
\seealso{\code{\link{ptmvnorm}}, \code{\link[mvtnorm]{pmvnorm}}, \code{\link[mvtnorm]{rmvnorm}}, \code{\link[mvtnorm]{dmvnorm}}}

\references{ 

Johnson, N./Kotz, S. (1970). Distributions in Statistics: Continuous Multivariate Distributions
\emph{Wiley & Sons}, pp. 70--73

Horrace, W. (2005). Some Results on the Multivariate Truncated Normal Distribution. 
\emph{Journal of Multivariate Analysis}, \bold{94}, 209--221

Jayesh H. Kotecha and Petar M. Djuric (1999). Gibbs Sampling Approach For Generation of Truncated Multivariate Gaussian Random Variables
\emph{IEEE Computer Society}, 1757--1760
}

\examples{
dtmvnorm(x=c(0,0))
dtmvnorm(x=c(0,0), mean=c(1,1), upper=c(0,0))

###########################################
#
# Example 1: 
# rejection sampling        
#
############################################

sigma <- matrix(c(4,2,2,3), ncol=2)
x <- rtmvnorm(n=500, mean=c(1,2), sigma=sigma, upper=c(1,0))
plot(x, main="samples from truncated bivariate normal distribution",
  xlim=c(-6,6), ylim=c(-6,6), 
  xlab=expression(x[1]), ylab=expression(x[2]))
abline(v=1, lty=3, lwd=2, col="gray")
abline(h=0, lty=3, lwd=2, col="gray")

###########################################
#
# Example 2: 
# Gibbs sampler for 4 dimensions
#
############################################

C = matrix(0.8, 4, 4)
diag(C)=rep(1, 4)
lower = rep(-4, 4)
upper = rep(-1, 4)

# acceptance rate alpha
alpha = pmvnorm(lower=lower, upper=upper, mean=rep(0,4), sigma=C)
alpha

# Gibbs sampler
X1=rtmvnorm(n=20000, mean = rep(0,4), sigma=C, lower=lower, upper=upper, algorithm="gibbs", burn.in.samples=100)
# Rejection sampling
X2=rtmvnorm(n=5000, mean = rep(0,4), sigma=C, lower=lower, upper=upper)

colMeans(X1)
colMeans(X2)

plot(density(X1[,1]), col="red", lwd=2, main="Gibbs vs. Rejection")
lines(density(X2[,1]), col="blue", lwd=2)
legend("topleft",legend=c("Gibbs Sampling","Rejection Sampling"), col=c("red","blue"), lwd=2)

}
\keyword{distribution}
\keyword{multivariate}