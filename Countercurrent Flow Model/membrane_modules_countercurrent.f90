module membrane_modules_countercurrent
  use, intrinsic :: IEEE_ARITHMETIC ! Import built-in arithmetic library
  use universal_parameters_countercurrent ! Import library with constant parameters

  implicit none

  ! Define the "module" class
  type, public :: module
    real*8, dimension(nx) :: x, y, dumpy, Ph, ytone, r
    real*8 :: Plf, Phf, Lf, Lo, Am, Vp
    real*8 :: purity, recovery, theta, astar
    real*8 :: Lf0, xf0
    character(LEN=8) :: mname
  end type

 contains
 ! --------------------------------------------------
 !                   SUBROUTINES
 ! --------------------------------------------------
  subroutine configure(mn, mname, xf, theta, xo, Lf)
    ! Get some values and configure the specified module instance
    ! mn = module instance
    ! mname = a name for the mn
    ! xf = molar fraction of co2 in the module's inlet at feed-side
    ! theta = theta value of module
    ! xo = molar fraction of co2 in the module's outlet at feed-side
    implicit none

    character(LEN=8) :: mname
    real*8, intent(in) :: xf, xo, Lf, theta
    integer :: i
    class(module) :: mn

    mn % mname      = mname
    mn % astar      = Pat/Pbt
    mn % x(1)       = xf
    mn % x(nx)      = xo
    mn % Lf         = Lf
    mn % theta      = theta

    Do i=1,nx
      mn%r(i) = mn%Plf/mn%Ph(i)
    End Do
  end subroutine configure

  subroutine solve(mn)
    implicit none
    class(module) :: mn
    real*8 :: dx, dz, xo_calc, error

    ! A trial and error simple algorithm is used here to reach the right xo for the specified theta
    Do While(.TRUE.)
      dx = (mn%x(nx)-mn%x(1))/nx ! Step of dydx integration

      Call integrateDyDxRK4(mn, -dx) ! Obtain y-x relation - Integration of dy/dx with Runke Kutta 4th order

      xo_calc = (mn%x(1) - mn%theta*mn%y(1))/(1d0 - mn%theta) ! Newly calculated xo

      mn%Vp = mn%theta*mn%Lf ! Permeate flowrate calculation
      mn%Lo = mn%Lf - mn%Vp ! Outlet at feed-side flowrate calculation

      ! Trial and error simple algorithm
      error = mn%x(nx) - xo_calc
      If (abs(error) .lt. 1d-3) Then
        Exit
      Else
        mn%x(nx) = mn%x(nx) - error/6d0 ! New guess
        Continue
      End If
    End Do

    Call integrateDAmDxRK4(mn, -dx) ! Calculate area - Integration of dAm/dx with Runke Kutta 4th order
    Call integrateDph2DxRK4(mn, -dx) ! Calculate pressure drop - Integration of dPh2/dx with Runke Kutta 4th order

    ! Calculate purity and recovery of CH4 and CO2 for module 2 and 3 respectively
    If(mn%mname .eq. 'Module 3') Then
      mn%purity = 100.d0 * mn%y(1)
      mn%recovery = 100.d0 * (mn%Vp*mn%y(1)) / (Lf0*xf0)
    Else
      mn%purity = 100.d0*(1d0 - mn%x(nx))
      mn%recovery = 100.d0*(mn%Lo*(1d0 - mn%x(nx))) / (Lf0*(1d0 - xf0))
    End If
  end subroutine solve

  subroutine integrateDyDxRK4(mn, h)
    ! Use Runke Kutta 4th order to calculate relation of y-x
    implicit none
    real*8 :: k1, k2, k3, k4, h
    real*8 :: a, b, c
    integer :: l
    class(module) :: mn

    Do l=nx,2,-1
      ! Calculation of y' for each step
      a = 1.d0 - mn%astar
      b = -1.d0 + mn%astar + 1.d0/mn%r(l) + mn%x(l)*(mn%astar - 1.d0)/mn%r(l)
      c = -(mn%astar*mn%x(l))/mn%r(l)

      mn%ytone(l) = (-b + sqrt(b**2 - 4*a*c))/(2*a)

      ! For the right end use function dydxo but for every other node, use dydx
      ! This is due to the fact that at the right end of the membrane,
      ! the differential equation is indeterminated.
      If(l.eq.nx) Then
        mn%y(l) = mn%ytone(l)
        k1 = h*dydxo(mn%x(l), mn%y(l), mn%astar, mn%r(l))
        k2 = h*dydxo(mn%x(l)+ 0.5d0*h, mn%y(l)+ 0.5d0*k1, mn%astar, mn%r(l))
        k3 = h*dydxo(mn%x(l)+ 0.5d0*h, mn%y(l)+ 0.5d0*k2, mn%astar, mn%r(l))
        k4 = h*dydxo(mn%x(l)+h, mn%y(l)+k3, mn%astar, mn%r(l))
      Else
        k1 = h*dydx(mn%x(l), mn%y(l), mn%astar, mn%r(l), mn%x(nx), mn%ytone(l))
        k2 = h*dydx(mn%x(l)+ 0.5d0*h,mn%y(l)+ 0.5d0*k1, mn%astar, mn%r(l), mn%x(nx), mn%ytone(l))
        k3 = h*dydx(mn%x(l)+ 0.5d0*h,mn%y(l)+ 0.5d0*k2, mn%astar, mn%r(l), mn%x(nx), mn%ytone(l))
        k4 = h*dydx(mn%x(l)+h,mn%y(l)+k3, mn%astar, mn%r(l), mn%x(nx), mn%ytone(l))
      End If
      ! Calculate next step y and x
      mn%y(l-1) = mn%y(l) + (1d0/6d0)*(k1 + 2*k2 + 2*k3 + k4)
      mn%x(l-1) = mn%x(l) + h
    End Do
  end subroutine integrateDyDxRK4

  subroutine integrateDAmDxRK4(mn, h)
    ! Use Runke Kutta 4th order to calculate membrane area
    implicit none
    real*8 :: k1, k2, k3, k4, h, Am
    integer :: k
    class(module) :: mn

    Am = 0.d0
    Do k=1,nx-1

      k1 = h*dAmdx(mn%x(k), mn%y(k), -mn%Lo, mn%Phf, mn%astar, mn%r(k), mn%x(nx))
      k2 = h*dAmdx(mn%x(k)+ 0.5d0*h, mn%y(k), -mn%Lo, mn%Phf, mn%astar, mn%r(k), mn%x(nx))
      k3 = h*dAmdx(mn%x(k)+ 0.5d0*h, mn%y(k), -mn%Lo, mn%Phf, mn%astar, mn%r(k), mn%x(nx))
      k4 = h*dAmdx(mn%x(k)+ h, mn%y(k), -mn%Lo, mn%Phf, mn%astar, mn%r(k), mn%x(nx))

      Am = Am + (1d0/6d0)*(k1 + 2*k2 + 2*k3 + k4)
      mn%x(k+1) = mn%x(k) - h
    End Do
    mn%Am = Am
  end subroutine integrateDAmDxRK4

  subroutine integrateDph2DxRK4(mn, h)
    ! Use of Runke Kutta 4th order to calculate the new pressure profile at feed side
    ! The function given in Geankoplis for pressure drop is dPh^2/dx
    implicit none
    real*8 :: k1,k2,k3,k4,h,L,visc_mix
    real*8, dimension(nx) :: ph_square
    integer :: f
    class(module) :: mn

    ph_square(1) = mn%Ph(1)**2
    Do f=1,nx-1
      L = (mn%Lo*(288.15d0/273.15d0)*(101325d0/mn%Ph(f))*1d6)*(mn%x(nx) - mn%y(f))/(mn%x(f) - mn%y(f))
      visc_mix = (mn%x(f)*(visc_co2 - visc_ch4) + visc_ch4)

      k1 = h*dph2dx(visc_mix, L, T, Di, N)

      ph_square(f+1) = ph_square(f) - k1
    End do
    Do f=1,nx
      mn%Ph(f) = sqrt(ph_square(f))
    End Do
  end subroutine integrateDph2DxRK4

  subroutine displayResults(mn)
    ! Display results of variables in mn object
    implicit none
    class(module) :: mn

    print *, '----------------------------------'
    print *, mn%mname
    print *, '----------------------------------'
    print *, 'theta = ',mn%theta
    print *, 'Lf = ',mn%Lf
    print *, 'Vp = ',mn%Vp
    print *, 'Lo = ',mn%Lo
    print *, 'yp = ',mn%y(1)
    print *, 'ytn =',mn%y(nx)
    print *, 'xf = ',mn%x(1)
    print *, 'xo = ',mn%x(nx)
    print *, ''
    print *, 'Purity = ',mn%purity
    print *, 'Recovery = ',mn%recovery
    print *, 'Am [m2] =',mn%Am*1d-4
    print *, '----------------------------------'
  end subroutine displayResults

! --------------------------------------------------
!                    FUNCTIONS
! --------------------------------------------------
  real*8 function dydxo(xo,yo,astar,r)
    ! Contains the differential equation for right end of membrane
    implicit none
    real*8, intent(in) :: xo,yo,astar,r

    dydxo = ((yo-xo)*(astar-(astar-1d0)*yo)) / ( (astar*(1d0-xo)*(xo-r*yo))- &
            xo*((1d0-xo) - r*(1d0-yo))-((yo-xo)*((astar-1d0)*(2d0*r*yo - &
             xo - r)-1d0) ))
    return
  end function dydxo

  real*8 function dydx(x,y,astar,r,xo,ytn)
    ! Contains the differential equation dy/dx
    implicit none
    real*8, intent(in) :: x,y,astar,r,xo,ytn
    real*8 :: cstm1, cstm2, a, b, c

    cstm1 = astar*(x - r*ytn)
    cstm2 = ((1d0 - x) - r*(1d0 - ytn))

    dydx = ((y-xo)/(x-xo))*(((1d0-y)*cstm1-y*cstm2)/((1d0-x)*cstm1-x*cstm2))
    return
  end function dydx

  real*8 function dAmdx(x, y, Lo, ph, astar, r, xo)
    ! Contains the differential equation dAm/dx
    implicit none
    real*8, intent(in) :: x, y, astar, r, xo, Lo, ph
    real*8 :: a, b, c, ytone, cstm1, cstm2

    a = 1d0 - astar
    b = -1d0 + astar + 1d0/r + x*(astar - 1.d0)/r
    c = -(astar*x)/r

    ytone = (-b + sqrt(b**2 - 4*a*c))/(2*a) ! [-]

    cstm1 = astar*(x - r*ytone) ! [-]
    cstm2 = ((1d0 - x) - r*(1d0 - ytone)) ! [-]

    dAmdx = ((Lo*(288.15d0/273.15d0)*(101325d0/ph)*1d6)/(Pbt*ph*0.00075d0))*( ((y-xo)/(x-y)) /( (1d0-x)*cstm1 - x*cstm2) )
    return
  end function dAmdx

  real*8 function dph2dx(visc_mix, L, T, Di, N)
    ! Contains the differential equation dPh^2/dx
    implicit none
    real*8 :: visc_mix,L,T,Di
    real*8 :: Ps, Ts
    integer :: N

    Ps =  101325d0 ! [Pa] - Because viscosity is in Pa.s
    Ts = 273.15d0 ! [K]

    dph2dx = (-256d0*visc_mix*L*T*Ps)/(pi*((Di*1d3)**4)*N*Ts)
    return
  end function
end module membrane_modules_countercurrent
