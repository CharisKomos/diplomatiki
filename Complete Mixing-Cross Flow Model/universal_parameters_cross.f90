module universal_parameters_cross
	implicit none

	! Totally universal constants
	integer, parameter :: nx  = 1d6 ! Number of nodes [-]
	integer, parameter :: Ntubes = 4000 ! Number of tubes per module [-]

	real*8, parameter :: Lf0 = 123.6d0 / 3600d0 ! System feed [Nm3/s]
  real*8, parameter :: xf0 = 0.293d0 ! Molar fraction of CO2 [-]

	real*8, parameter :: Rg    = 8.314d0  ! Universal gas constant [Pa*m3/mol/K]
	real*8, parameter :: tol   = 1.d-5    ! Tolerance [-]
	real*8, parameter :: Prel  = 1.0d0    ! Relative pressure [bar]
	real*8, parameter :: pi    = dacos(-1.d0)
	real*8, parameter :: phs   = 1.013d5
	real*8, parameter :: visc_co2 = 0.85d-5 ! [Pa*s]
	real*8, parameter :: Pat 	 = 103.4d-6 ! [ccSTP / cm2*s*cmHg]
	real*8, parameter :: Pbt 	 = 2.63d-6 ! [ccSTP / cm2*s*cmHg]
end module universal_parameters_cross
