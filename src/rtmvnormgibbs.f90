! Gibbs sampling from a truncated multinormal distribution
!
! see Kotecha et al. (1999):
! Kotecha, J. H. & Djuric, P. M.
! "Gibbs sampling approach for generation of truncated multivariate Gaussian random variables",
! IEEE Computer Society, IEEE Computer Society, 1999, 1757-1760
!
! Code written by Stefan Wilhelm <wilhelm@financial.com> as part of the R package tmvtnorm.
! (http://CRAN.R-project.org/package=tmvtnorm)
!
! To cite package tmvtnorm in publications use:
!
!  Stefan Wilhelm, Manjunath B G (2010). tmvtnorm: Truncated
!  Multivariate Normal Distribution. R package version 1.1-5.
!
! A BibTeX entry for LaTeX users is
!
!  @Manual{,
!    title = {{tmvtnorm}: Truncated Multivariate Normal Distribution},
!    author = {Stefan Wilhelm and Manjunath B G},
!    year = {2010},
!    note = {R package version 1.1-5},
!    url = {http://CRAN.R-project.org/package=tmvtnorm},
!  }
! 
!
! @param n number of random sample to generate by Gibbs sampling
! @param d dimension (d >= 2)
! @param mean mean vector of dimension d (d x 1)
! @param sigma covariance matrix (d x d)
! @param lower lower truncation points (d x 1)
! @param upper upper truncation points (d x 1)
! @param x0 Startvektor (d x 1)
! @param burnin Number of Burn-in samples to be discarded
! @param thinning thinning factor for thinning the Markov chain
! @return return value X --> vektor (n * d) --> can be coerced into a (n x d) matrix
subroutine rtmvnormgibbs(n, d, mean, sigma, lower, upper, x0, burnin, thinning, X)

IMPLICIT NONE

integer :: n, d, i, j, k, l, ind = 0, error, burnin, thinning

! subindex "-i"
integer, dimension(d-1) :: minus_i

double precision :: unifrnd, qnormr, pnormr, u, q, prob, Fa, Fb, mu_i, s2
double precision, dimension(n*d), INTENT(INOUT) :: X
double precision, dimension(d-1)     :: s3
double precision, dimension(d)       :: x0, xr, mean, lower, upper, sd

! Kovarianzmatrix sigma und Partitionen Sigma_i, sigma_ii und S
double precision, dimension(d, d)    :: sigma
double precision, dimension(d, d-1)  :: Sigma_i
double precision                     :: sigma_ii
double precision, dimension(d-1,d-1) :: S
! S_inv (d-1 x d-1) ist die Inverse von S
double precision, dimension(d-1,d-1) :: S_inv

! Liste von d mal 1 x (d-1) Matrizen = d x (d-1) Matrix
double precision, dimension(d, d-1)   :: P

! Deklarationen fürs Matrix-Invertieren mit LAPACK-Routinen (Dimension d-1)
double precision, dimension( d-1 )    :: work
! ipiv = pivot indices
integer, dimension( d-1 )             :: ipiv
! lda = leading dimension
integer                               :: m, lda, lwork, info

! initialise R random number generator
call rndstart()

m    =d-1
lda  =d-1
lwork=d-1
ind  = 0

! Partitioning of sigma
! sigma   = [ sigma_ii   Sigma_i     ]
! (d x d)   [ (1 x 1)    (1 x d-1)   ]
!           [ Sigma_i'   S           ]
!           [ (d-1 x 1)  (d-1 x d-1) ]
! List of conditional variances sd(i) can be precalculated
do i = 1,d
  ! subindex "-i"
  minus_i  = (/ (j, j=1,i-1), (j, j=i+1,d) /)

  S            = sigma(minus_i, minus_i) ! (d-1) x (d-1)
  sigma_ii     = sigma(i,i)              ! 1 x 1
  Sigma_i(i,:) = sigma(i, minus_i)       ! 1 x (d-1)

  ! Matrix S --> S_inv umkopieren
  do k=1,(d-1)
    do l=1,(d-1)
      S_inv(k,l)=S(k,l)
    end do
  end do

  ! Matrix invertieren
  ! LU-Faktorisierung (Dreieckszerlegung) der Matrix S_inv
  call dgetrf( m, m, S_inv, lda, ipiv, info )

  ! Inverse der LU-faktorisierten Matrix S_inv
  call dgetri( m, S_inv, lda, ipiv, work, lwork, info )
  P(i,:) = pack(matmul(Sigma_i(i,:), S_inv), .TRUE.)  ! (1 x d-1) %*% (d-1 x d-1) = (1 x d-1)
  s2 = 0
  do j = 1,d-1
    s2 = s2 + P(i,j) * Sigma_i(i,j)
  end do
  sd(i)        = sqrt(sigma(i,i) - s2) ! (1 x d-1) * (d-1 x 1) --> sd[[i]] ist (1,1)
end do

! start value
xr = x0

! Actual number of samples to create:
! #burn-in-samples + n * #thinning-factor

!For all samples n times the thinning factor
do j = 1,(burnin + n * thinning)

  ! For all dimensions
  do i = 1,d
    ! Berechnung von bedingtem Erwartungswert und bedingter Varianz:
    ! bedingte Varianz hängt nicht von x[-i] ab!
    ! subindex "-i"
    minus_i  = (/ (k, k=1,i-1), (k, k=i+1,d) /)

    ! mu_i     = mean(i)    + P[[i]] %*% (x(-i) - mean(-i))
    s3(1:(d-1))= xr(minus_i) - mean(minus_i)
    s2 = 0
    do k = 1,d-1
      s2 = s2 + P(i,k) * s3(k)
    end do
    mu_i       = mean(i) + s2

    Fa         = pnormr(lower(i), mu_i, sd(i), 1, 0)
    Fb         = pnormr(upper(i), mu_i, sd(i), 1, 0)
    u          = unifrnd()
    prob       = u * (Fb - Fa) + Fa
    q          = qnormr(prob, 0.0d0, 1.0d0, 1, 0)
    xr(i)      = mu_i + sd(i) * q

    ! Nur für j > burnin samples aufzeichnen, Default ist thinning = 1
    ! bei Thinning nur jedes x-te Element nehmen
    if (j > burnin .AND. mod(j - burnin,thinning) == 0) then
      ind        = ind + 1
      X(ind)     = xr(i)
    end if
  end do
end do

! reset R random number generator
call rndend()
end subroutine rtmvnormgibbs
