! Based on the Geankoplis model for the cocurrent flow model for assymetric membranes
! Pressure loss is calculated

program countercurrentflow
  use membrane_modules_countercurrent

  implicit none
  real*8 :: Lf, xf, Lf_first, xf_old
  real*8 :: c1, c2, c3
  real*8 :: start, finish
  real*8 :: theta1, theta2, theta3

  ! Declare new instances of class "module"
  type(module) :: module1
  type(module) :: module2
  type(module) :: module3

  module1%Phf = 1.d5*(10.d0 + Prel) ! Feed side pressure of module 1 [Pa]
  module1%Ph  = module1%Phf
  module1%Plf = 1.d5*(1.d0  + Prel) ! Permeate side pressure of module 1 [Pa]

  module2%Phf = 1.d5*(10d0 + Prel) ! Feed side pressure of module 2 [Pa]
  module2%Ph  = module2%Phf
  module2%Plf = 1.d5*(0d0 + Prel) ! Permeate side pressure of module 2 [Pa]

  module3%Phf = 1.d5*(1d0 + Prel) ! Feed side pressure of module 3 [Pa]
  module3%Ph  = module3%Phf
  module3%Plf = 1.d5*(0d0 + Prel) ! Permeate side pressure of module 3 [Pa]

  ! Start counting calculation time from this point
  start = 0d0
  Call CPU_TIME(start)

  ! Some first estimations of xo e.g. xo1 = xf1 * c1
  c1 = 0.1d0
  c2 = 0.04d0
  c3 = 0.9899d0

  ! Theta values for each module
  theta1 = 0.374269d0
  theta2 = 0.346064d0
  theta3 = 0.515625d0

  ! Prepare variables before loop
  Lf = Lf0
  xf = xf0
  xf_old = 0d0

  ! Do the calculations until steady state
  Do While(.TRUE.)
    Call configure(module1, 'Module 1', xf, theta1, c1*xf, Lf)
    Call solve(module1)
    Call configure(module2, 'Module 2', module1%x(nx), theta2, c2*module1%x(nx), module1%Lo)
    Call solve(module2)
    Call configure(module3, 'Module 3', module1%y(1), theta3, c3*module1%y(1), module1%Vp)
    Call solve(module3)

    ! Recycling - Recalculate Lf and xf
    xf_old = xf
    Lf = (Lf0 + module2%Vp + module3%Lo)
    xf = (Lf0*xf0 + module2%Vp*module2%y(1) + module3%Lo*module3%x(nx))/Lf

    ! Check for xf convergence
    If(abs(xf - xf_old) .lt. 1d-5) Then
      Exit
    Else
      Continue
    End If
  End Do

  ! Display calculation time
  Call CPU_TIME(finish)
  print *, 'Process ended in ',finish-start,'sec'

  ! Display results for the whole system
  Call displayResults(module2)
  Call displayResults(module3)
End Program countercurrentflow
