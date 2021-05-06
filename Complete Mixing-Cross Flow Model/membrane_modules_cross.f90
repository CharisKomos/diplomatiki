module membrane_modules_cross
  use, intrinsic :: IEEE_ARITHMETIC
  use universal_parameters_cross

  implicit none

  type, public :: module
    real*8, dimension(nx) :: z, x, V, L, ph
    real*8 :: Am_study,Lf,theta,vp_total0,vp_total_i
    real*8 :: xf,xo,pl,ph0,r,astar,yp,Vp,Lf0,xf0,max, maxpur,maxrec,Lo

    real*8 :: Am, Amt, x0f
    real*8 :: gamma2
    real*8 :: y1, xN, purity, recovery
    character(LEN=8) :: mname
  end type

contains
  ! --------------------------------------------------
  !                   SUBROUTINES
  ! --------------------------------------------------
  subroutine configure(mn, mname, Lf, xf, theta)
    implicit none

    character(LEN=8) :: mname
    real*8 :: temp, Am_study, Lf, xf, theta, dz
    integer :: i

    class(module) :: mn

    ! Initialize some properties of mn object
    mn % mname     = mname
    mn % Am_study  = Am_study
    mn % Lf        = Lf
    mn % Lf0       = Lf
    mn % xf        = xf
    mn % xf0       = xf
    mn % theta     = theta
    mn % astar     = Pat/Pbt
    mn % vp_total0 = mn%Lf * mn%theta

    ! Divide space to 1/nx nodes and calculate the Vpi for each node
    dz = 1 / nx
    mn%vp_total_i = mn%vp_total0 / nx

    ! Initialize first value of vectors z and Vp
    ! Then set the whole vectors in regard to step dz and possible dP
    mn%z(1) = 0.d0
    mn%V(1) = mn%vp_total_i
    If(nx .gt. 1) Then
      do i=2, nx
        mn%z(i) = mn%z(i-1) + dz
        mn%V(i) = mn%vp_total_i
      end do
    End if
  end subroutine configure

  subroutine solve(mn)
    implicit none
    class(module) :: mn
    real*8 :: vp_total, vpy_total,Am_total
    real*8 :: a1,b1,c1
    real*8, dimension(1) :: yp
    integer :: i

    vp_total = 0.d0
    vpy_total = 0.d0
    Am_total = 0.d0

    Do i=1, nx
      If(i .gt. 1) mn%ph(i) = mn%ph(i-1)

      ! Calculate gamma (or r) and theta for node i
      mn%r = mn%pl / mn%ph(i)
      mn%theta = mn%V(i) / mn%Lf

      ! Ensure that theta is not greater than 1
      If(mn%theta .gt. 1.d0) Then
        print *, 'theta > 1 at : ', i, mn%theta
        stop
      End If

      ! Coefficients for yp of node i calculation
      a1 = mn%theta+mn%r-mn%r*mn%theta - mn%astar*mn%theta &
            - mn%astar*mn%r + mn%astar*mn%r*mn%theta
      b1 = 1.d0 - mn%theta - mn%xf -mn%r+mn%r*mn%theta &
           + mn%astar*mn%theta+mn%astar*mn%r-mn%astar*mn%r*mn%theta &
           + mn%astar*mn%xf
      c1 = -mn%astar*mn%xf

      yp(1) = (-b1 + sqrt(b1*b1 - 4*a1*c1)) / (2*a1)

      ! Calculation of residue fraction at feed side
      mn%xo = (mn%xf - mn%theta*yp(1)) / (1.d0 - mn%theta)

      ! Calculation of active membrane area needed
      mn%Am = (mn%theta*(mn%Lf*(288.15d0/273.15d0)*(101325d0/mn%ph(i))*1d6)*yp(1))/(Pat*(mn%ph(i)*mn%xo - mn%pl*yp(1)))/0.00075d0 ! cm2
      Am_total = Am_total + mn%Am

      mn%V(i) = mn%Lf*mn%theta
      mn%Lo = mn%Lf*(1.d0 - mn%theta)
      mn%yp = yp(1)

      vp_total = vp_total + mn%V(i)
      vpy_total = vpy_total + mn%V(i)*yp(1)

      mn%xf = mn%xo ! Next node's xf will be equivalent to the xo at the outlet of previous
      mn%Lf = mn%Lo ! Next node's Lf will be equivalent to the Lo at the outlet of previous

    End Do

    mn%Vp = vp_total
    mn%yp = vpy_total / vp_total
    mn%theta = vp_total / mn%Lf0
    mn%Amt = Am_total*1d-4 !m2

    If(mn%mname .eq. 'Module 3') Then
      mn%purity = 100.d0 * mn%yp
      mn%recovery = 100.d0 * (mn%Vp*mn%yp) &
                        / (Lf0*xf0)
    Else
      mn%purity = 100.d0 * (1.d0 - mn%xo)
      mn%recovery = 100.d0 * (mn%Lo*(1.d0 - mn%xo)) &
                        / (Lf0*(1d0-xf0))
    End If
  end subroutine solve

  subroutine displayResults(mn)
    ! Display all the appropriate results for mn object
    implicit none
    class(module) :: mn

    print *, '----------------------------------'
    print *, mn%mname
    print *, '----------------------------------'
    print *, 'theta = ',mn%theta
    print *, 'Lf = ',mn%Lf
    print *, 'Vp = ',mn%Vp
    print *, 'Lo = ',mn%Lo
    print *, 'yp = ',mn%yp
    print *, 'xf = ',mn%xf0
    print *, 'xo = ',mn%xo
    print *, ''
    print *, 'Purity = ',mn%purity
    print *, 'Recovery = ',mn%recovery
    print *, 'Am =',mn%Amt
    print *, '----------------------------------'
  end subroutine displayResults

end module membrane_modules_cross
