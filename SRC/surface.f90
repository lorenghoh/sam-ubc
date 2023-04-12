subroutine surface()
	
use vars
use params
use microphysics, only: micro_field, index_water_vapor
use rad, only: swdsvisxy, swdsnirxy, swdsvisdxy, swdsnirdxy, lwdsxy, coszrsxy
implicit none
	
real t_s, q_s, w_s, r_s, z0_s, t_h, q_h, ta_h, u_h, v_h, u_h0
real taux0, tauy0, xlmo
real diag_ustar, coef, coef1
real tmp1, tmp2, tmpu(0:nx,ny), tmpv(nx,0:ny)
real zr(nx,ny),presr(nx,ny),ur(nx,ny),vr(nx,ny),tr(nx,ny),qvr(nx,ny),t_sfc(nx,ny)
integer i,j,k
integer tag(2),req(2),irank,count
real(8) buffer(2), buffer1(2)
real buffers(max(nx,ny),2)
logical flag(2)
real, external :: qsatw,qsati

call t_startf ('surface')


if(.not.SFC_FLX_FXD) then
  
  if(sstxy(1,1).le.-100.) then
     print*,'surface: sst is undefined. Quitting...'
     call task_abort()
  end if

  if(OCEAN.or.LAND.and.SLM) then

       do j=1,ny
         do i=1,nx

          if(landmask(i,j).eq.0) then

           u_h = 0.5*(u(i+1,j,1)+u(i,j,1)) + ug
           v_h = 0.5*(v(i,j+YES3D,1)+v(i,j,1)) + vg
           ta_h = tabs(i,j,1)
           q_h = qv(i,j,1)
           t_h = tabs(i,j,1)*prespot(1)
           t_s = (sstxy(i,j)+t00)*prespoti(1)
           q_s = salt_factor*qsatw(sstxy(i,j)+t00,presi(1))

           call oceflx(rho(1), u_h, v_h, ta_h, q_h, t_h, z(1)-zi(1), t_s, q_s, &
                        fluxt0, fluxq0, taux0, tauy0)

           fluxbu(i,j) = taux0/rho(1)
           fluxbv(i,j) = tauy0/rho(1)
           fluxbt(i,j) = fluxt0
           fluxbq(i,j) = fluxq0

          end if  

         end do
       end do

	
  end if ! OCEAN


  if(LAND) then

       if(SLM) then

          do j=1,ny
           do i=1,nx
             zr(i,j) = z(1)
             presr(i,j) = pres(1)
             ur(i,j) = 0.5*(u(i+1,j,1)+u(i,j,1)) + ug
             vr(i,j) = 0.5*(v(i,j+YES3D,1)+v(i,j,1)) + vg
             tr(i,j) = tabs(i,j,1)
             qvr(i,j) = qv(i,j,1)
             t_sfc(i,j) = sstxy(i,j)+t00
           end do
          end do

          call run_slm(ur, vr, precinst, qvr, tr, t_sfc, zr, presr, &
              swdsvisxy, swdsnirxy, swdsvisdxy, swdsnirdxy, lwdsxy, coszrsxy, &
              precinstsoil(:,:), fluxbu, fluxbv, fluxbq, fluxbt)
              sstxy(:,:) = t_sfc(:,:)-t00

        else

          do j=1,ny
           do i=1,nx

            if(landmask(i,j).eq.1) then

              u_h = 0.5*(u(i+1,j,1)+u(i,j,1)) + ug
              v_h = 0.5*(v(i,j+YES3D,1)+v(i,j,1)) + vg
              q_h = qv(i,j,1)
              t_h = tabs(i,j,1)*prespot(1)
              t_s = (sstxy(i,j)+t00)*prespoti(1)
              q_s = salt_factor*qsatw(sstxy(i,j)+t00,presi(1))
 
              call landflx(t_h, t_s, q_h, q_s, u_h, v_h, z(1)-zi(1), z0, &
                          fluxbt(i,j), fluxbq(i,j), fluxbu(i,j), fluxbv(i,j))
            end if

           end do
          end do

       end if ! SLM

  end if ! LAND

! 
! compute stresses at velocity positions:

  tmpu(1:nx,1:ny) = fluxbu(1:nx,1:ny)
  tmpv(1:nx,1:ny) = fluxbv(1:nx,1:ny)
  if(dompi) then
     call task_receive_float(buffers(:,1),max(nx,ny),req(1))
     call task_receive_float(buffers(:,2),max(nx,ny),req(2))
     call task_bsend_float(rankee,fluxbu(nx,1:ny),ny,133)
     call task_bsend_float(ranknn,fluxbv(1:nx,ny),nx,134)
     count = 2
     flag(1:2) = .false.
     do while(count.gt.0)
      do k=1,2
        if(.not.flag(k)) then
          call task_test(req(k),flag(k),irank,tag(k))
          if(flag(k)) then
            if(tag(k).eq.133) then
             tmpu(0,1:ny) = buffers(1:ny,k)
             count=count-1
            else if(tag(k).eq.134) then
             tmpv(1:nx,0) = buffers(1:nx,k)
             count=count-1
            else
              print*,'surface: wrong tag. Should be 133 or 134, got ',tag
              call task_abort
            end if
          end if
        end if
      end do
     end do
  else
     tmpu(0,1:ny) = fluxbu(nx,1:ny)
     tmpv(1:nx,0) = fluxbv(1:nx,ny)
  end if
  do j=1,ny
   do i=1,nx
    fluxbu(i,j) = 0.5*(tmpu(i,j)+tmpu(i-1,j)) 
   end do
  end do

  do j=1,ny
   do i=1,nx
    fluxbv(i,j) = 0.5*(tmpv(i,j)+tmpv(i,j-YES3D)) 
   end do
  end do

  call task_barrier()

end if! .not.SFC_FLX_FXD



if(SFC_FLX_FXD) then

  u_h0 = max(1.,sqrt((u0(1)+ug)**2+(v0(1)+vg)**2))

  if(.not.SFC_TAU_FXD) then
    if(OCEAN) z0 = 0.0001  ! for LAND z0 should be set in namelist (default z0=0.035)

    tau0 = diag_ustar(z(1),  &
                bet(1)*(fluxt0+epsv*(t0(1)-gamaz(1))*fluxq0),u_h0,z0)**2  

  end if ! .not.SFC_TAU_FXD

  fluxbu(:,:) = -(u(1:nx,1:ny,1)+ug)/u_h0*tau0
  fluxbv(:,:) = -(v(1:nx,1:ny,1)+vg)/u_h0*tau0

  fluxbt(:,:) = fluxt0
  fluxbq(:,:) = fluxq0

end if ! SFC_FLX_FXD

!
! Homogenize the surface scalar fluxes if needed for sensitivity studies
!
   if(dosfchomo) then

	fluxt0 = 0.
	fluxq0 = 0.
	do j=1,ny
         do i=1,nx
	   fluxt0 = fluxt0 + fluxbt(i,j)
	   fluxq0 = fluxq0 + fluxbq(i,j)
         end do
        end do
	fluxt0 = fluxt0 / float(nx*ny)
	fluxq0 = fluxq0 / float(nx*ny)
        if(dompi) then
            buffer(1) = fluxt0
            buffer(2) = fluxq0
            call task_sum_real8(buffer,buffer1,2)
	    fluxt0 = buffer1(1) /float(nsubdomains)
	    fluxq0 = buffer1(2) /float(nsubdomains)
        end if ! dompi
	fluxbt(:,:) = fluxt0
	fluxbq(:,:) = fluxq0

   end if

shf_xy(:,:) = shf_xy(:,:) + fluxbt(:,:) * dtfactor
lhf_xy(:,:) = lhf_xy(:,:) + fluxbq(:,:) * dtfactor

call t_stopf ('surface')

end




! ----------------------------------------------------------------------
!
! DISCLAIMER : this code appears to be correct but has not been
!              very thouroughly tested. If you do notice any
!              anomalous behaviour then please contact Andy and/or
!              Bjorn
!
! Function diag_ustar:  returns value of ustar using the below 
! similarity functions and a specified buoyancy flux (bflx) given in
! kinematic units
!
! phi_m (zeta > 0) =  (1 + am * zeta)
! phi_m (zeta < 0) =  (1 - bm * zeta)^(-1/4)
!
! where zeta = z/lmo and lmo = (theta_rev/g*vonk) * (ustar^2/tstar)
!
! Ref: Businger, 1973, Turbulent Transfer in the Atmospheric Surface 
! Layer, in Workshop on Micormeteorology, pages 67-100.
!
! Code writen March, 1999 by Bjorn Stevens
!
! Code corrected 8th June 1999 (obukhov length was wrong way up,
! so now used as reciprocal of obukhov length)

      real function diag_ustar(z,bflx,wnd,z0)

      implicit none
      real, parameter      :: vonk =  0.4   ! von Karmans constant
      real, parameter      :: g    = 9.81   ! gravitational acceleration
      real, parameter      :: am   =  4.8   !   "          "         "
      real, parameter      :: bm   = 19.3   !   "          "         "
      real, parameter      :: eps  = 1.e-10 ! non-zero, small number

      real, intent (in)    :: z             ! height where u locates
      real, intent (in)    :: bflx          ! surface buoyancy flux (m^2/s^3)
      real, intent (in)    :: wnd           ! wind speed at z
      real, intent (in)    :: z0            ! momentum roughness height

      integer :: iterate
      real    :: lnz, klnz, c1, x, psi1, zeta, rlmo, ustar

      lnz   = log(z/z0) 
      klnz  = vonk/lnz              
      c1    = 3.14159/2. - 3.*log(2.)

      ustar =  wnd*klnz
      if (bflx /= 0.0) then 
        do iterate=1,4
          rlmo   = -bflx * vonk/(ustar**3 + eps)   !reciprocal of
                                                   !obukhov length
          zeta  = z*rlmo
          if (zeta > 0.) then
            ustar =  vonk*wnd  /(lnz + am*zeta)
          else
            x     = sqrt( sqrt( 1.0 - bm*zeta ) )
            psi1  = 2.*log(1.0+x) + log(1.0+x*x) - 2.*atan(x) + c1
            ustar = wnd*vonk/(lnz - psi1)
          end if
        end do
      end if

      diag_ustar = ustar

      return
      end function diag_ustar
! ----------------------------------------------------------------------

