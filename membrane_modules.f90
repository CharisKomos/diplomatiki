module membrane_modules
   use universal_parameters
	 use, intrinsic :: IEEE_ARITHMETIC
   implicit none

   type, public :: flow
      real*8 :: temperature, pressure, flowrate
      real*8, dimension(3) :: composition ! (CO2, CH4, other)
   end type

   type, public :: membrane_module
      character(LEN=8) :: mname
      integer :: Ntubes
      real*8 :: temps, temp, Diam, Diam4, Diam_OD, Am0
      real*8 :: Am, lm, Lfth, Lfts, Lf, xf, x0, x0f
      real*8 :: fha, fhb, Pl0, Pl, Pl00, Ph0, Ph00, L00, Pat, Pbt
      real*8 :: r, gamma2, astar, Kap1, Kap2, Kap22
      real*8 :: theta, y1, xN, purity, recovery, purity0, recovery0, Amt
   	  real*8 :: dpress, flr
      type(flow) :: inlet
      type(flow) :: outlet_feed
      type(flow) :: outlet_permeate
   end type

   !type, extends(membrane_module), public :: membrane_results
   !	real*8 :: theta, y1, xN, purity, recovery, purity0, recovery0, Amt
   !	real*8 :: dpress, flr
   !end type

  contains

   subroutine configure(mn, mname, temp, Am0, Diam, &
			 Diam_OD, Ntubes, Lfth, xf, x0)
        implicit none
        character(Len=10) :: mname
        real*8 :: Am0, temp, Diam, Diam_OD, Lfth, xf, x0
        integer :: Ntubes
        class(membrane_module) :: mn

        mn % mname = mname
        mn % temp = temp
        mn % temps = temp + 273.d0
        mn % Am0 = Am0
        mn % Am = Am0 / dfloat(Ntubes)
        mn % Diam = Diam
        mn % Diam4 = Diam**4
        mn % Diam_OD = Diam_OD
        mn % Ntubes = Ntubes
        mn % lm = mn % Am / (pi*Diam_OD)
        mn % Lfth = Lfth
        mn % Lfts = Lfth / 3600.d0
        mn % Lf = mn % Lfts / dfloat(Ntubes)
        mn % xf = xf
        mn % x0 = x0
        mn % x0f = 0.95 * mn % xf

        mn % inlet % composition(1) = mn % xf
        mn % inlet % composition(2) = 1.d0 - mn % xf
        mn % inlet % composition(3) = 0.d0
				mn % outlet_feed % composition(1) = mn % x0

		    !Call prompt_user_for_variables(mn)
				If(mn%mname == 'Module 1') Then
					mn%Ph00 = 10.d0
					mn%Pl00 = 2.d0
					mn%L00  = 0.1d0
					mn%Pat 	= 85.3d-6* (76.d0/1.d5) / 100.d0 ! Pat is now Nm3/s/Pa
					mn%Pbt 	= 3.2d-6* (76.d0/1.d5) / 100.d0
					mn%fha 	= 0.952d0
					mn%fhb	= 1.d0

				Else If(mn%mname == 'Module 2') Then
					mn%Ph00 = 10.d0
					mn%Pl00 = 0.d0
					mn%L00  = 0.1d0
					mn%Pat 	= 85.3d-6* (76.d0/1.d5) / 100.d0 ! Pat is now Nm3/s/Pa
					mn%Pbt 	= 3.2d-6* (76.d0/1.d5) / 100.d0
					mn%fha 	= 0.952d0
					mn%fhb	= 1.d0

				Else If(mn%mname == 'Module 3') Then
					mn%Ph00 = 2.d0
					mn%Pl00 = 0.d0
					mn%L00  = 0.015d0
					mn%Pat 	= 95.0d-6* (76.d0/1.d5) / 100.d0 ! Pat is now Nm3/s/Pa
					mn%Pbt 	= 2.8d-6* (76.d0/1.d5) / 100.d0
					mn%fha 	= 1.d0
					mn%fhb	= 1.d0

				Else
					print *, 'An error occured with module name'

				End If


			  mn % Ph0 = (1.d5)*(mn % Ph00 + Prel) ! Pa
			  mn % Pl0 = (1.d5)*(mn % Pl00 + Prel) ! Pa

			  mn % astar  = mn % Pat / mn % Pbt
			  mn % Pl     = mn % Pl0
			  mn % r      = mn % Pl / mn % Ph0
			  mn % gamma2 = mn % r

			  mn % Kap1 = pi*mn%Diam_OD*(mn%lm/mn%Lf)*mn%Pbt*mn%Ph0
			  mn % Kap2 = 128.d0 * visc_co2 * (phs * mn%temp / mn%temps) &
					* mn%Lf * mn%lm / (pi*(mn%Ph0 ** 2)) * mn%Diam4
			  mn % Kap22 = 128.d0 * visc_co2 * Rg * mn%temp &
					 * (mn%Lf*1000.d0/22.4d0) * mn%lm / (pi*(mn%Ph0 ** 2)) * mn%Diam4
   end subroutine configure


	 subroutine prompt_user_for_variables(mn)

			implicit none
			class(membrane_module) :: mn

			character(LEN=1) :: option
			real*8 :: optionA, optionB

			! Clear all and show the interface
			Call System('clear')
			print *, '-------------------------------------------'
			print *, 'Pressure configuration for module ',mn%mname
			print *, '-------------------------------------------'
			print *, 'Please choose an option for the pressure : '
			print *, '(a) Ph = 10bar'
			print *, '(b) Ph = 6bar'
			print *, '(c) Ph = 2bar'
			print *, '-------------------------------------------'
			print *, 'I choose : '
			read *, option

			If (option.eq.'a') Then
				mn%Ph00 = 10.d0
			Else If (option.eq.'b') Then
				mn%Ph00 = 6.d0
			Else
				mn%Ph00 = 2.d0
			End If

			Select Case (nint(mn%Ph00))
				Case (10)
					mn%Pat 	= 85.3d-6
					mn%Pbt 	= 3.2d-6
					mn%fha 	= 0.952d0
					mn%fhb	= 1.d0
					optionA = 0.1d0
					optionB = 0.08d0
				Case (6)
					mn%Pat 	= 86.0d-6
					mn%Pbt 	= 3.0d-6
					mn%fha 	= 0.9712d0
					mn%fhb	= 1.d0
					optionA = 0.08d0
					optionB = 0.5d0
				Case (2)
					mn%Pat 	= 95.0d-6
					mn%Pbt 	= 2.8d-6
					mn%fha 	= 1.d0
					mn%fhb	= 1.d0
					optionA = 0.015d0
					optionB = 0.003d0
			End Select

			mn%Pat = mn%Pat * (76.d0/1.d5) / 100.d0 ! Pat is now Nm3/s/Pa
			mn%Pbt = mn%Pbt * (76.d0/1.d5) / 100.d0 ! Pbt is now Nm3/s/Pa

			print *, 'These options are '&
				'recommended for the flowrate :'
			print *, '(a) L00 = ', optionA
			print *, '(b) L00 = ', optionB
			print *, '-------------------------------------------'
			print *, 'I choose : '
			read *,  option
			If (option.eq.'a') Then
				mn%L00 = optionA
			Else If (option.eq.'b') Then
				mn%L00 = optionB
			Else
				print *, 'Setting default data'&
					' L00 = ', optionA
				mn%L00 = optionA
			End If

			print *, '-------------------------------------------'
			print *, 'Please provide pressure value for the permeate side :'
			read *, mn%Pl00
			print *, 'Permeate side pressure is set to Pl00 = ', mn%Pl00,'bar'
			print *, '-------------------------------------------'
			pause
	 end subroutine prompt_user_for_variables


   subroutine solve(mn, z, dz)!, dx0)

	   implicit none

		 real*8, dimension(nm) :: x, y, z, V, L, ph, dph, dph1
		 real*8 L0, L1, L2, rBC, rr0, rr1, dz, del_r, correction, err !,dx0
		 real*8 err_dp
		 integer :: g, h, iterp, iter

	   Class(membrane_module) :: mn

		 ! You can make the parametric analysis based on the change of this variable

		  x(N) = mn % outlet_feed % composition(1)
	   	V(N) = 0.d0

	   	do g=1, N
	   	  ph(g)  = 1.d0
	   	  dph(g) = 0.d0
	   	end do

	   	L0 = mn % L00
	   	L1 = L0 + L0/1000.d0

	   	rBC = mn % inlet % composition(1)

	   	iterp = 0

		 outer_loop: Do While (.TRUE.)
		   	  iterp = 0
		   	  iterp = iterp + 1

			    iter = 0
		 inner_loop: Do While (.TRUE.)
			    iter = iter + 1

				  L(N) = L0
				  Call loop_through_L(0, mn,x,y,z,V,L,ph,dph,dph1,dz,rr0,rr1)
					L(N) = L1
					Call loop_through_L(1, mn,x,y,z,V,L,ph,dph,dph1,dz,rr0,rr1)

					del_r = rr1 - rr0
					correction = ((rBC - rr1)/del_r)*(L1 - L0)

					L2 = L1 + correction

					err = dabs(L2-L1)

					L0 = L1
					L1 = L2

					If (dabs(correction).gt.tol) Then
					   Cycle inner_loop
					Else
					   Exit inner_loop
					End If

		    End Do inner_loop

		    err_dp = 0.d0

		    do h=2, N
					If (IEEE_IS_FINITE(dph1(h) - dph(h))) Then
					  err_dp = err_dp + (dph1(h) - dph(h))**2
					  dph(h) = dph1(h)
					  ph(h) = ph(h-1) + dph(h)
					Else
						Exit outer_loop
					End If
		    end do

		    If (iterp.gt.1 .and. ph(N).le.mn%gamma2) stop

				err_dp = dsqrt(err_dp)/dfloat(N)

				If (err_dp.gt.tol) Then
					 Cycle outer_loop
				Else
					 Exit outer_loop
				End If

		 End Do outer_loop


       mn % Amt = dfloat(mn%Ntubes)*pi* mn%Diam_OD * mn%lm
       mn % dpress = (ph(1) - ph(N))*mn%ph0 / 1.d5
       mn % theta = V(1) / L(1)
       mn % recovery = 100.d0 * (V(1)*y(1))/(L(1)*x(1))
       mn % purity = 100.d0 * y(1)
       mn % recovery0 = 100.d0*(L(N)*(1.d0-x(N))) &
       					/(L(1)*(1.d0-x(1)))
       mn % purity0 = 100.d0*(1.d0 - x(N))
       mn % flr = L(1) * mn%Lfth ! This is the inlet flowrate
       mn % y1 = y(1)
       mn % xN = x(N)

			 mn % inlet % flowrate = mn % flr
			 mn % outlet_feed % flowrate = L(N)
			 mn % outlet_permeate % flowrate = V(1)

       mn % outlet_feed % composition(1) = mn % xN
			 mn % outlet_feed % composition(2) = 1.d0 - mn % xN
			 mn % outlet_feed % composition(3) = 0.d0
       mn % outlet_permeate % composition(1) = mn % y1
			 mn % outlet_permeate % composition(2) = 1.d0 - mn % y1
			 mn % outlet_permeate % composition(3) = 0.d0
   end subroutine solve


   subroutine loop_through_L(flag, mn,x,y,z,V,L,ph,dph,dph1,dz,rr0,rr1)
      implicit none
      Class(membrane_module) :: mn
      real*8, dimension(nm) :: x, y, z, V, L, ph, dph, dph1
      real*8 :: a,b,c,yt,dfr1,dfr2,fx,fy,fv,fl,Vav,Lav,rr0,rr1,phav,dz
      integer :: flag, f

      do f = N, 2, -1
         a = mn%gamma2 * ( 1.d0 - mn%astar )
         b = ph(f) * mn%fhb + (mn%astar*mn%fha - mn%fhb)*ph(f)*x(f) &
         		+ mn%gamma2 * ( mn%astar - 1.d0 )
         c = -mn%astar * ph(f) * mn%fha*x(f)

         yt= (-b + sqrt(b**2 - 4*a*c)) / (2.d0*a)

         If (f .eq. N) y(f) = yt

         dfr1 = ph(f) * mn%fha * x(f) - mn%gamma2*yt
         dfr2 = ph(f) * mn%fhb * (1.d0 - x(f)) - mn%gamma2 *(1.d0 - yt)

         fx = -(mn%Kap1)*(mn%astar * (1.d0-x(f))*dfr1 - x(f)*dfr2)
         fy = -(mn%Kap1)*(mn%astar * (1.d0-y(f))*dfr1 - y(f)*dfr2)
         fl = -mn%Kap1 * (mn%astar*dfr1 + dfr2)
         fv = fl

         V(f - 1) = V(f) - fv*dz
         L(f - 1) = L(f) - fl*dz

         Vav = (V(f-1) + V(f)) / 2.d0
         Lav = (L(f-1) + L(f)) / 2.d0

         x(f - 1) = x(f) - (fx/Lav)*dz
         y(f - 1) = y(f) - (fy/Vav)*dz

         If (flag .eq. 1) Then
            phav = (ph(f-1) + ph(f)) / 2.d0
            dph1(f) = - (mn%Kap2 * Lav / phav)*dz
            rr1 = x(1)
         Else
         	rr0 = x(1)
         End If

      end do
   end subroutine loop_through_L


	 ! When two streams get mixed up, calculate temperature, pressure and consistency
   subroutine mix_flows(mod1,mod2,mod3)
	  implicit none
		Class(membrane_module) :: mod1
		Class(membrane_module) :: mod2
		Class(membrane_module) :: mod3

		mod1%inlet%composition(1) = (mod1%inlet%flowrate*mod1%inlet%composition(1) &
					 + mod2%outlet_permeate%flowrate*mod2%outlet_permeate%composition(1) &
									+ mod3%outlet_feed%flowrate*mod3%outlet_feed%composition(1)) &
												/ (mod1%inlet%flowrate + mod2%outlet_permeate%flowrate &
														+ mod3%outlet_feed%flowrate)
   end subroutine mix_flows
end module membrane_modules
