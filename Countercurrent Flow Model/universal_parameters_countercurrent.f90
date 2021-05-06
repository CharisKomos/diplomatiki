module universal_parameters_countercurrent
  implicit none

  integer, parameter :: nx  = 10000! Number of nodes for the solution
  integer, parameter :: N   = 4000 ! Number of hollow fibers

  real*8, parameter :: Lf0  = 123.6d0 / 3600d0 ! System feed [Nm3/s]
  real*8, parameter :: xf0  = 0.293d0 ! Molar fraction of CO2 [-]

  real*8, parameter :: Rg   = 8.314d0  ! Universal gas constant -- Pa*m3/mol/K
	real*8, parameter :: T 		= 293.15d0 ! K = 20oC
  real*8, parameter :: Prel = 1.d0
  real*8, parameter :: Pat  = 103.4d-6 ! ccSTP / cm2*s*cmHg
  real*8, parameter :: Pat_mol = 346.39d-12 ! mol*m / m2*s*Pa
  real*8, parameter :: Pbt  = 2.63d-6 ! ccSTP / cm2*s*cmHg
  real*8, parameter :: Pbt_mol = 8.8105d-12 ! mol*m / m2*s*Pa
  real*8, parameter :: visc_ch4 = 1.1d-5 ! Pa*s
  real*8, parameter :: visc_co2 = 0.85d-5 ! Pa*s
  real*8, parameter :: pi   = dacos(-1.d0)
  real*8, parameter :: D0   = 420d-6 ! m
  real*8, parameter :: Dt   = 280d-6 ! m
  real*8, parameter :: Di   = 240d-6 ! m
  real*8, parameter :: d    = 0.5d-6 ! m - Random assumption
  real*8, parameter :: len  = 1.d0  ! 1.2954d0 ! m = 51inches
  real*8, parameter :: gc   = 1.d0 ! Kg*m / N*s2 -> Dimensionless
end module universal_parameters_countercurrent
