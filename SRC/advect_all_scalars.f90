subroutine advect_all_scalars()

  use vars
  use microphysics
  use entrainment
  use sgs
  use tracers
  use params, only: dotracers
  implicit none
  real dummy(nz)
  integer k


!---------------------------------------------------------
!      advection of scalars :

    ! UBC ENT
    if (doentrainment) then
      call entrainment_advective_fluxes()
    end if
    ! End UBC ENT

     call advect_scalar(t,tadv,twle,t2leadv,t2legrad,twleadv,.true.)

!
!    Advection of microphysics prognostics:
!

     do k = 1,nmicro_fields
       if(   k.eq.index_water_vapor             &! transport water-vapor variable no metter what
         .or. docloud.and.flag_precip(k).ne.1    & ! transport non-precipitation vars
         .or. doprecip.and.flag_precip(k).eq.1 ) then
          if(flag_advect(k).eq.1) &
             call advect_scalar(micro_field(:,:,:,k),mkadv(:,k),mkwle(:,k),dummy,dummy,dummy,.false.)
       end if
     end do

!
!    Advection of sgs prognostics:
!

     if(dosgs.and.advect_sgs) then
       do k = 1,nsgs_fields
           call advect_scalar(sgs_field(:,:,:,k),sgsadv(:,k),sgswle(:,k),dummy,dummy,dummy,.false.)
       end do
     end if


!
!   Precipitation fallout:
!
    precinst(:,:) = 0.
    if(doprecip) then

       total_water_prec = total_water_prec + total_water()

       call micro_precip_fall()

       total_water_prec = total_water_prec - total_water()


    end if

    ! UBC ENT
    if (doentrainment) then
      call entrainment_diffusive_fluxes()
      call calculate_entrainment_fluxes()
    end if
    ! End UBC ENT

 ! advection of tracers:

     if(dotracers) then

        do k = 1,ntracers
         call advect_scalar(tracer(:,:,:,k),tradv(:,k),trwle(:,k),dummy,dummy,dummy,.false.)
        end do

     end if

end subroutine advect_all_scalars
