module universal_parameters
	implicit none

	! Totally universal constants
	integer, parameter :: nm  = 100001
	integer, parameter :: N   = 10001
	integer, parameter :: npr = (N-1) / 10
	integer, parameter :: nx  = 20

	real*8, parameter :: Rg    = 8.314d0  ! Universal gas constant -- Pa*m3/mol/K
	real*8, parameter :: tol   = 1.d-6    ! Tolerance
	real*8, parameter :: Prel  = 1.0d0    ! Relative pressure
	real*8, parameter :: pi    = dacos(-1.d0)
	real*8, parameter :: phs   = 1.013d5
	real*8, parameter :: visc_co2 = 1.1d-5 ! Pa*s
        ! real, parameter :: visc_sul =
	! real, parameter :: t	  = 1.d0



end module universal_parameters
