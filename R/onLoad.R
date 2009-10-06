# On Load of the Package load shared object/DLL for Fortran code
.First.lib <- function(libname, pkgname) 
{
  library.dynam("tmvtnorm", pkgname)
}
