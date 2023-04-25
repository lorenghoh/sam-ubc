module tetrahedral_entrainment

use grid
use tetrahedral_interpolation
use vars

CONTAINS

!--------------------------------

  subroutine calculate_entrainment_surface(n_values, &
                                           values, &
                                           mask, &
                                           u_area_frac, &
                                           v_area_frac, &
                                           w_area_frac, &
                                           volume)
    ! Calculates the cell volume and surface fraction that is within the
    ! entrainment surface.

    ! Input variables
    ! number of values used to define the location of the entrainment surface
    integer n_values
    ! variable values used to define the location of the entrainment surface
    real, dimension(n_values, dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) :: values
    ! mask which is 1 for grid point inside the entrainment surface
    logical mask(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    ! Output variables
    ! fraction of the x, y, and z cell surfaces inside the entrainment surface
    real u_area_frac(dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm)
    real v_area_frac(dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm)
    real w_area_frac(dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz )
    ! cell volume inside the entrainment surface
    real volume(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)

    ! Calculate the fraction of the u, v, and w cell surfaces which are inside
    ! the entrainment surface
    call calculate_surface(n_values, &
                           values, &
                           mask, &
                           u_area_frac, &
                           v_area_frac, &
                           w_area_frac)

    ! Calculate the volume of each grid cell that is within the entrainment
    ! surface
    call calculate_volume(n_values, &
                          values, &
                          mask, &
                          volume)
    
  end subroutine calculate_entrainment_surface

!--------------------------------

  subroutine calculate_fluxes_through_entrainment_surface(uflux, &
                                                          vflux, &
                                                          wflux, &
                                                          u_area_frac, &
                                                          v_area_frac, &
                                                          w_area_frac, &
                                                          volume, &
                                                          fluxes)
                                                          
    ! Using the divergence theorem, calculate fluxes through the 
    ! entrainment surface.
    ! Postive surface_flux means the surface is entraining.

    ! Input variables
    ! flux of tracer/mass through each grid cell wall in (kg/m3)(m/s)
    real uflux(dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm)
    real vflux(dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm)
    real wflux(dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz )
    ! fraction of the x, y, and z cell surfaces inside the entrainment surface
    real u_area_frac(dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm)
    real v_area_frac(dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm)
    real w_area_frac(dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz )
    ! cell volume inside the entrainment surface
    real volume(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)

    ! Output variables
    ! net flux of tracer/mass through the entrainment surface in kg/m3/s
    real fluxes(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)

    integer i, j, k
   
    ! Calculates surface flux using divergence theorem
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          fluxes(i, j, k) = &
                     uflux(i+1, j, k)*u_area_frac(i+1, j, k)/dx &
                   - uflux(i, j, k)*u_area_frac(i, j, k)/dx &
                   + vflux(i, j+1, k)*v_area_frac(i, j+1, k)/dy &
                   - vflux(i, j, k)*v_area_frac(i, j, k)/dy &
                   + wflux(i, j, k+1)*w_area_frac(i, j, k+1)/(dz*adz(k)) &
                   - wflux(i, j, k)*w_area_frac(i, j, k)/(dz*adz(k))&
              - volume(i, j, k)*((uflux(i+1, j, k) - uflux(i, j, k))/dx &
                               + (vflux(i, j+1, k) - vflux(i, j, k))/dy &
                               + (wflux(i, j, k+1) - wflux(i, j, k))/(dz*adz(k)))
        end do
      end do
    end do

  end subroutine calculate_fluxes_through_entrainment_surface

!--------------------------------

  subroutine calculate_E_D(volume, old_volume, &
                           tracer, old_tracer, &
                           tracer_fluxes, E_tracer, D_tracer)
    ! Calculate fluxes due to cloud volume changes
    ! Input variables
    ! cell volume inside the entrainment surface
    real volume(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    ! cell volume inside the entrainment surface during the previous time step
    real old_volume(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    ! cell tracer value
    real tracer(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    ! cell tracer value during the previous time step
    real old_tracer(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    ! Output variables
    ! Flux of tracer through the entrainment surface
    real tracer_fluxes(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    ! Tracer entrainment/detrainment rate
    real E_tracer(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real D_tracer(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)

    integer i, j, k
    real volume_flux
    real temp_tracer
 
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          temp_tracer = (tracer(i,j,k) + old_tracer(i, j, k))/2.
          ! calculate cloud tracer volume change/dt in kg/s/m3
          volume_flux = rho(k)*temp_tracer &
                       *(volume(i, j, k) - old_volume(i, j, k))/dtn

          ! add to surface flux
          tracer_fluxes(i, j, k) = tracer_fluxes(i, j, k) + volume_flux

          ! if > 0 add to E, set D=0, if < 0 add to D, set E=0
          if (temp_tracer > 0.) then
            E_tracer(i, j, k) = max(0., tracer_fluxes(i, j, k))
            D_tracer(i, j, k) = max(0., -tracer_fluxes(i, j, k))
          else
            E_tracer(i, j, k) = min(0., tracer_fluxes(i, j, k))
            D_tracer(i, j, k) = min(0., -tracer_fluxes(i, j, k))
          end if
        end do
      end do
    end do   

  end subroutine calculate_E_D
 
end module tetrahedral_entrainment
