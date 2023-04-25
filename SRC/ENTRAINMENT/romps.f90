module romps

use grid

contains

  subroutine calculate_activity_mask(active_cells, activity_mask, adjacent_n)
    implicit none
   
    real active_cells(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    logical activity_mask(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real adjacent_n(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)

    integer i, j, k

    ! Calculate activity changes
    do k=2,nzm-1
      do j=1,ny
        do i=1,nx
          ! Because adjacent_n is used multiple times for each entrianed
          ! variable, all changes to adjacent_n are performed here for
          ! computational efficiency.
          ! If the cell was inactive at the previous timestep, reset 
          ! adjacent_n to zero.
          if (.not. activity_mask(i, j, k)) adjacent_n(i,j,k) = 0.
          ! If an active cell is adjacent to an inactive cell, that
          ! cell is participating in entrianment or detrainment.
          activity_mask(i, j, k) = &
                  ((active_cells(i, j, k) .ne. active_cells(i+1, j, k)) &
              .or. (active_cells(i, j, k) .ne. active_cells(i-1, j, k)) &
              .or. (active_cells(i, j, k) .ne. active_cells(i, j+1, k)) &
              .or. (active_cells(i, j, k) .ne. active_cells(i, j-1, k)) &
              .or. (active_cells(i, j, k) .ne. active_cells(i, j, k+1)) &
              .or. (active_cells(i, j, k) .ne. active_cells(i, j, k-1)))
          ! If the cell is active, then it is adjacent to the entrainment
          ! surface in this time step.
          if (activity_mask(i, j, k)) then
            adjacent_n(i, j, k) = adjacent_n(i, j, k) + 1.
          end if
        end do
      end do
    end do

  end subroutine calculate_activity_mask

!--------------------------

  subroutine calculate_romps_fluxes(u_flux, v_flux, w_flux, &
                                    flux_activity, active, active_mask)
    implicit none

    real u_flux(dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm)
    real v_flux(dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm)
    real w_flux(dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz )
    real flux_activity(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real active(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    logical active_mask(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)

    integer i, j, k
 
    ! Given the tracer fluxes and the activities, calculate the net
    ! flux through the entrainment surface.
    do k=2,nzm-1
      do j=1,ny
        do i=1,nx
          if (active_mask(i,j,k)) then
            flux_activity(i, j, k) = flux_activity(i, j, k) &
                + max(0., u_flux(i+1, j,   k  ))*active(i,   j,   k  )/dx &
                + min(0., u_flux(i+1, j,   k  ))*active(i+1, j,   k  )/dx &
                - max(0., u_flux(i,   j,   k  ))*active(i-1, j,   k  )/dx &
                - min(0., u_flux(i,   j,   k  ))*active(i,   j,   k  )/dx &
                + max(0., v_flux(i,   j+1, k  ))*active(i,   j,   k  )/dy &
                + min(0., v_flux(i,   j+1, k  ))*active(i,   j+1, k  )/dy &
                - max(0., v_flux(i,   j,   k  ))*active(i,   j-1, k  )/dy &
                - min(0., v_flux(i,   j,   k  ))*active(i,   j,   k  )/dy &
                + max(0., w_flux(i,   j,   k+1))*active(i,   j,   k  )/(dz*adz(k)) &
                + min(0., w_flux(i,   j,   k+1))*active(i,   j,   k+1)/(dz*adz(k)) &
                - max(0., w_flux(i,   j,   k  ))*active(i,   j,   k-1)/(dz*adz(k)) &
                - min(0., w_flux(i,   j,   k  ))*active(i,   j,   k  )/(dz*adz(k))
          end if
        end do
      end do
    end do

  end subroutine calculate_romps_fluxes

!----------------------

  subroutine calculate_romps_forcing(force, force_activity, &
                                     active, active_mask)
    implicit none

    real force(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real force_activity(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real active(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    logical active_mask(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)

    integer i, j, k

    ! Subtract forcing terms (buoyancy, rainfall) from the 
    ! fluxes through the entrainment surface.
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          if (active_mask(i,j,k)) then
            force_activity(i,j,k) = force_activity(i,j,k) &
                                  - force(i,j,k)*active(i,j,k)
          end if
        end do
      end do
    end do

  end subroutine calculate_romps_forcing

!----------------------

  subroutine calculate_romps_activity(mask, old_mask, &
                                      n_values, values, old_values, &
                                      active_cells)
    implicit none

    logical mask(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    logical old_mask(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    integer n_values
    real, dimension(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, n_values) :: values
    real, dimension(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, n_values) :: old_values
    real active_cells(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
 
    integer i, j, k, n
    real x1, x2

    ! Calculate new activities
    do k=1,nzm
      do j=dimy1_s,dimy2_s
        do i=dimx1_s,dimx2_s
          if (mask(i, j, k) .neqv. old_mask(i, j, k)) then
            ! if active_cell has flipped, interpolate flip time 
            ! using old grid values
            if (old_mask(i, j, k)) then
              x1 = 1
              do n=1,n_values
                if (values(i, j, k, n) > 0.) then
                  x2 = 1.
                else
                  x2 = old_values(i, j, k, n) & 
                       / (old_values(i, j, k, n) - values(i, j, k, n))
                end if
                x1 = min(x1, x2)
              end do      
              active_cells(i, j, k) = x1
            else
              x1 = 0
              do n=1,n_values
                if (old_values(i, j, k, n) > 0.) then
                  x2 = 0.
                else
                  x2 = old_values(i, j, k, n) &
                       / (old_values(i, j, k, n) - values(i, j, k, n))
                end if
                x1 = max(x1, x2)
              end do
              active_cells(i, j, k) = 1. - x1
            end if      
          else
            if (mask(i, j, k)) then
              active_cells(i, j, k) = 1.
            else
              active_cells(i, j, k) = 0.
            end if
          end if
        end do
      end do
    end do    

  end subroutine calculate_romps_activity

!----------------------

  subroutine calculate_romps_entrainment(active_cells, active_cells_old, &
                                         tracer, tracer_old, &
                                         E_tracer, D_tracer, &
                                         adjacent_n, activity_mask, &
                                         tracer_mean, tracer_activity)
    use vars, only: rho
    implicit none

    real active_cells(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real active_cells_old(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real tracer(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real tracer_old(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real E_tracer(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real D_tracer(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real adjacent_n(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    logical activity_mask(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real tracer_mean(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real tracer_activity(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)

    integer i, j, k
    real d_rho_A_tracer_dt

    E_tracer = 0.
    D_tracer = 0.

    ! Calculate activity changes
    do k=1,nzm
      do j=1,ny
        do i=1,nx

          if (activity_mask(i,j,k)) then
  
            ! If cell is currently active, add tracer time derivatives to 
            ! the activity fluxes
            d_rho_A_tracer_dt = &
                rho(k)*active_cells(i, j, k)*tracer(i, j, k)/dtn &
              - rho(k)*active_cells_old(i, j, k)*tracer_old(i, j, k)/dtn
            tracer_activity(i, j, k) = &
                tracer_activity(i, j, k) + d_rho_A_tracer_dt
            ! Record the mean tracer value to help decide if the fluxes should
            ! be assigned to entrainment or detrainment.
            tracer_mean(i, j, k) = tracer_mean(i, j, k) + tracer_old(i, j, k)
  
          else
  
            ! Cell is currently inactive
            ! If this cell was active in the previous timestep, dump activities
            if (adjacent_n(i, j, k) > 0.) then
            
              if (tracer_mean(i, j, k) > 0.) then
                ! If the tracer has been, on average, positive, then positive activity
                ! implies entrainment, and negative activity implies detrainment
                E_tracer(i, j, k) = max(0., tracer_activity(i, j, k))
                D_tracer(i, j, k) = max(0., -tracer_activity(i, j, k))
              else
                ! If the tracer has been, on average, negative, then positive activity
                ! implies negative tracer detrainment, and negative activity implies 
                ! negative tracer entrainment
                E_tracer(i, j, k) = min(0., tracer_activity(i, j, k))
                D_tracer(i, j, k) = min(0., -tracer_activity(i, j, k))
              end if
            
              tracer_activity(i, j, k) = 0.
              tracer_mean(i, j, k) = 0.
            end if
  
          end if        

        end do
      end do
    end do    

!    if (masterproc) print *, "tracer_activity: ", tracer_activity(9, 22:26, 10)
!    if (masterproc) print *, "tracer_mean: ", tracer_mean(9, 22:26, 10)
!    if (masterproc) print *, "E_tracer: ", E_tracer(9, 22:26, 10)
!    if (masterproc) print *, "D_tracer: ", D_tracer(9, 22:26, 10)

  end subroutine calculate_romps_entrainment

end module romps
