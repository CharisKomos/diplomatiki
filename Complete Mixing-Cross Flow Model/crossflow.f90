program crossflow
  use membrane_modules_cross

  implicit none

  real*8, dimension(nx) :: z, V, L
  real*8 :: Am_study, Lf, xf, xf_old
  real*8 :: ytot, thetat, theta1, theta2, theta3
  real*8 :: a1, b1, c1, start, finish
  integer :: i

  type(module) :: module1
  type(module) :: module2
  type(module) :: module3

  module1%ph(1) = 1.d5*(10.d0 + Prel)
  module1%pl = 1.d5*(1.d0 +  Prel)

  module2%ph(1) = 1.d5*(10.d0 + Prel)
  module2%pl = 1.d5*(0.d0 +  Prel)

  module3%ph(1) = 1.d5*(1.d0 + Prel)
  module3%pl = 1.d5*(0.d0 +  Prel)

  start = 0d0
  Call CPU_TIME(start)

  ! Theta values for each module
  theta1 = 0.374269d0
  theta2 = 0.346064d0
  theta3 = 0.515625d0

  ! Prepare variables before loop
  Lf = Lf0
  xf = xf0
  xf_old   = 0.d0

  ! Do the calculations until steady state
  Do While (.TRUE.)
    Call configure(module1,'Module 1', Lf, xf, theta1)
    Call solve(module1)
    Call configure(module2,'Module 2', module1%Lo, module1%xo, theta2)
    Call solve(module2)
    Call configure(module3,'Module 3', module1%Vp, module1%yp, theta3)
    Call solve(module3)

    ! Recycling - Recalculate Lf and xf
    xf_old = module1%xf0
    Lf = Lf0 + module2%Vp + module3%Lo
    xf = (xf0*Lf0 + module2%Vp*module2%yp + module3%Lo*module3%xo) / Lf

    ! Check for xf convergence
    If (abs(xf - xf_old) .lt. tol) Then
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
end program crossflow
