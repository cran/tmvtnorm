#############################################
#
# Beispiel
#
#############################################

Sigma = matrix(c(1, 0.95,
                 0.95, 1), 2, 2)
mu = c(0,0)
X = rmvnorm(1000, mu, Sigma)

# Trunkiere x2 <= 0
X.trunc = X[X[,2]<0,]

par(mfrow=c(2,2))
plot(X.trunc)

# Randdichte für x1
plot(density(X.trunc[,1]))

# Berechnung der Randdichte für x1
x <- seq(-5, 5, by=0.01)
fx <- c()
for (i in 1:length(x))
{
  fx[i] = dtmvnorm.marginal(x[i], n=1, mean=mu, sigma=Sigma, lower=c(-Inf,-Inf), upper=c(Inf,0))
}

lines(x, fx, lwd=2, col="red")

plot(x, fx, type="l", lwd=2, col="red", main="Randdichte")

Rprof("prof.out")
for (i in 1:length(x))
{
  fx[i] = dtmvnorm.marginal(x[i], n=1, mean=mu, sigma=Sigma, lower=c(-Inf,-Inf), upper=c(Inf,0))
}
Rprof(NULL)
summaryRprof("prof.out")

# Vektorisiert
fx <- c()
Rprof("prof2.out")
fx = dtmvnorm.marginal(x, n=1, mean=mu, sigma=Sigma, lower=c(-Inf,-Inf), upper=c(Inf,0))
Rprof(NULL)
summaryRprof("prof2.out")
lines(x, fx, lwd=2, col="blue")


#############################################
#
# Beispiel 2 : trivariat
#
#############################################

Sigma = outer(1:3,1:3,pmin)
mu = c(0,0,0)

X = rmvnorm(2000, mu, Sigma)

# Trunkiere x2 <= 0, x3 <= 0
X.trunc = X[X[,2]<0 & X[,3]<0,]

par(mfrow=c(2,2))
plot(X.trunc)

# Randdichte für x1
plot(density(X.trunc[,1]))

# Berechnung der Randdichte für x1
x <- seq(-5, 5, by=0.01)
fx <- c()
for (i in 1:length(x))
{
  fx[i] = dtmvnorm.marginal(x[i], n=1, mean=mu, sigma=Sigma, lower=c(-Inf,-Inf,-Inf), upper=c(Inf,0,0))
}

lines(x, fx, lwd=2, col="red")
plot(x, fx, type="l", lwd=2, col="red", main="Randdichte")

Rprof("prof.out")
for (i in 1:length(x))
{
  fx[i] = dtmvnorm.marginal(x[i], n=1, mean=mu, sigma=Sigma, lower=c(-Inf,-Inf,-Inf), upper=c(Inf,0,0))
}
Rprof(NULL)
summaryRprof("prof.out")

fx<-c()
Rprof("prof.out")
fx = dtmvnorm.marginal(x, n=1, mean=mu, sigma=Sigma, lower=c(-Inf,-Inf,-Inf), upper=c(Inf,0,0))
Rprof(NULL)
summaryRprof("prof.out")
lines(x, fx, lwd=2, col="blue")


#############################################
#
# Beispiel 3 : bivariat, aber trunkiert in 2. Dimension
#
#############################################

i = 2

sigma=0.4

# n=3 Bewertungstage
T <- c(0.5, 0.7, 1)

# Strikes
X <- c(100, 100)

Sigma = sigma^2*outer(T[1:i],T[1:i],pmin)
mu = c(0,0)

X = rmvnorm(2000, mu, Sigma)

# Trunkiere x1 <= 0, x2 <= 0
X.trunc = X[X[,1]<0 & X[,2]<0,]

par(mfrow=c(2,2))
plot(X.trunc)

# Randdichte für x2
plot(density(X.trunc[,2]))

# Berechnung der Randdichte für x2
x <- seq(-5, 5, by=0.01)
fx<-c()
Rprof("prof.out")
fx = dtmvnorm.marginal(x, n=2, mean=mu, sigma=Sigma, lower=c(-Inf,-Inf), upper=c(0,0))
Rprof(NULL)
summaryRprof("prof.out")
lines(x, fx, lwd=2, col="blue")

# Testen, ob Dichte zu 1 integriert
#integrate(dtmvnorm.marginal, lower=-Inf, upper=0, n=2, mean=mu, sigma=Sigma, lower=c(-Inf,-Inf), upper=c(0,0))
