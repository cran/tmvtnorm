y = rtnorm.gibbs(1000, a=-3, b=0)
plot(ecdf(y))

# Example 3: Profiling of Gibbs sampling: 10000 samples
Rprof("rtnorm.gibbs.out")
X=rtnorm.gibbs(100000, a=-3, b=0)
Rprof(NULL)
summaryRprof("rtnorm.gibbs.out")

###############################

X=rtmvnorm.gibbs(n=10, mean = c(0,0,0), sigma = diag(3), lower = c(0,0,0))

# Akzeptanzrate alpha aus der Multivariaten Normalverteilung bestimmen
alpha = pmvnorm(lower=c(0,0,0), mean=c(0,0,0), sigma=diag(3))
alpha

# Example 4: Profiling of Gibbs sampling: d=3, n=10000 samples ~ 2.2 second
Rprof("rtmvnorm.gibbs.out")
X1=rtmvnorm.gibbs(n=10000, mean = c(0,0,0), sigma = diag(3), lower = c(0,0,0))
Rprof(NULL)
summaryRprof("rtmvnorm.gibbs.out")

# Example 5: Profiling of rejection sampling: d=3, n=10000 samples ~ 2.2 second
Rprof("rtmvnorm.out")
X2=rtmvnorm(n=10000, mean = c(0,0,0), diag(3), lower = c(0,0,0))
Rprof(NULL)
summaryRprof("rtmvnorm.out")

colMeans(X1)
colMeans(X2)

###############################

###############################
#
# 10-dimensionales Beispiel aus Kotecha
#
###############################
C = matrix(0.8, 10, 10)
diag(C)=rep(1, 10)
#lower = rep(-4, 10)
#upper = rep(-3, 10)

lower = rep(-1, 10)
upper = rep(+1, 10)


# Akzeptanzrate alpha aus der Multivariaten Normalverteilung bestimmen
alpha = pmvnorm(lower=lower, upper=upper, mean=rep(0,10), sigma=C)
alpha

# Example 6: Profiling of Gibbs sampling: d=10, n=10000 samples ~ 6.8 second
Rprof("rtmvnorm.gibbs.out")
X1=rtmvnorm.gibbs(n=10000, mean = rep(0,10), sigma=C, lower=lower, upper=upper)
Rprof(NULL)
summaryRprof("rtmvnorm.gibbs.out")

# Example 7: Profiling of rejection sampling: d=10, n=10000 samples ~ 2.2 second
Rprof("rtmvnorm.out")
X2=rtmvnorm(n=10000, mean = rep(0,10), sigma=C, lower = lower, upper=upper)
Rprof(NULL)
summaryRprof("rtmvnorm.out")

colMeans(X1)
colMeans(X2)

###############################
#
# 4-dimensionales Beispiel aus Kotecha (statt 10-Dimensionen)
#
###############################

# --> Akzeptanzrate ist so niedrig 1e-5, dass Rejection Sampling in die Knie geht,
# selbst bei wenigen Samples (n=500)
# --> Gibbs Sampling ist davon unbeeindruckt!

C = matrix(0.8, 4, 4)
diag(C)=rep(1, 4)
lower = rep(-4, 4)
upper = rep(-1, 4)

# Akzeptanzrate alpha aus der Multivariaten Normalverteilung bestimmen
alpha = pmvnorm(lower=lower, upper=upper, mean=rep(0,4), sigma=C)
alpha

# Example 6: Profiling of Gibbs sampling: d=10, n=10000 samples ~ 6.8 second
Rprof("rtmvnorm.gibbs.out")
X1=rtmvnorm.gibbs(n=5000, mean = rep(0,4), sigma=C, lower=lower, upper=upper)
Rprof(NULL)
summaryRprof("rtmvnorm.gibbs.out")

# Example 7: Profiling of rejection sampling: d=10, n=10000 samples ~ 2.2 second
Rprof("rtmvnorm.out")
X2=rtmvnorm(n=5000, mean = rep(0,4), sigma=C, lower=lower, upper=upper)
Rprof(NULL)
summaryRprof("rtmvnorm.out")

colMeans(X1[-(1:20),])
colMeans(X2)

plot(density(X1[-(1:20),1]), col="red", lwd=2, main="Gibbs vs. Rejection")
lines(density(X2[,1]), col="blue", lwd=2)
legend("topleft",legend=c("Gibbs Sampling","Rejection Sampling"), col=c("red","blue"), lwd=2)

a=acf(X1[,1])
a2=acf(X1[,2])
plot(a$lag[-1], a$acf[-1], type="l", col="red")
lines(a2$lag[-1], a2$acf[-1], col="blue")

# Entwicklung des MCMC-Algorithmus
plot(X1[,1], type="l")

plot(X1[,1],X1[,2])

############################################

X1 = rtmvnorm.gibbs(n=5000, mean = rep(0,4), sigma=C, lower=lower, upper=upper, burn.in=10, start.value=rep(-5,4))
X1 = rtmvnorm.gibbs(n=5000, mean = rep(0,4), sigma=C, lower=lower, upper=upper, burn.in=10, start.value=rep(-3.5,4))


X1 = rtmvnorm(n=5000, mean = rep(0,4), sigma=C, lower=lower, upper=upper, algorithm="gibbs", burn.in=10, start.value=rep(-3.5,4))
