! Gibbs sampling from a truncated multinormal distribution
!
! see Kotecha et al. (1999):
! Kotecha, J. H. & Djuric, P. M.
! "Gibbs sampling approach for generation of truncated multivariate Gaussian random variables",
! IEEE Computer Society, IEEE Computer Society, 1999, 1757-1760
!
! @param n Anzahl der MC-Wiederholungen
! @param d Dimension (d >= 2)
! @param x0 Startvektor
! @return Rückgabewert X --> Vektor (n * d)
subroutine rtmvnormgibbs(n, d, mean, sigma, lower, upper, x0, X)

IMPLICIT NONE

integer :: n, d, i, j, k, l, ind = 0, error

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

xr = x0

!For all samples
do j = 1,n
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
    ind        = ind + 1
    X(ind)     = xr(i)
  end do
end do

! reset R random number generator
call rndend()
end subroutine rtmvnormgibbs
