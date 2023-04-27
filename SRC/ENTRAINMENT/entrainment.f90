module entrainment

use vars, only: u1, v1, w1, rho, rhow
use sgs, only: tkh, nsgs_fields, nsgs_fields_diag
use params
use grid
use tetrahedral_entrainment
use romps
use microphysics, only: nmicro_fields, micro_field

implicit none

real qt_minus_qsat(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real delta_thetav(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real w_on_rho(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)

real, external :: qsatw, qsati

!-----------------------

! (nx, ny, nz, n) Tracer fluxes
! n=1 is mass
! n=2 is specific humidity
! n=3 is liquid water moist static energy
! n=4 is vertical velocity
real u_adv_fluxes(dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm, 4)
real v_adv_fluxes(dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm, 4)
real w_adv_fluxes(dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz,  4)

real u_dif_fluxes(dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm, 4)
real v_dif_fluxes(dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm, 4)
real w_dif_fluxes(dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz,  4)

! Store the advective w fluxes from the previous 3 timesteps for
! use in the Adams-Bashforth timestepping calculation
real u_w_adv_flux_store(dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm, 3)
real v_w_adv_flux_store(dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm, 3)
real w_w_adv_flux_store(dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz , 3)

! Store the advective w fluxes from the previous 3 timesteps for
! use in the Adams-Bashforth timestepping calculation
real u_w_dif_flux_store(dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm, 3)
real v_w_dif_flux_store(dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm, 3)
real w_w_dif_flux_store(dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz , 3)

! Changes in tracer variables due to forcing terms, advection terms
! and the total variable change during the current time step.
real d_force(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 4)
real d_adv(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 4)
real d_dif(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 4)
real d_total(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 4)

! Area of the grid cell faces occupied by the cloud core
real u_area_core_frac(dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm)
real v_area_core_frac(dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm)
real w_area_core_frac(dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz )

! Fraction of the grid cell volumes occupied by the cloud core
real volume_core(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
! Fraction of the grid cell volumes occupied by the cloud core
! in the previous time step
real old_volume_core(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)

! Tracer fluxes through the cloud core surface
real fluxes_core(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 4)

! Tracer entrainment/detrainment through the cloud core surface
real E_tetra_core(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 4)
real D_tetra_core(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 4)


! Area of the grid cell faces occupied by the cloud
real u_area_cloud_frac(dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm)
real v_area_cloud_frac(dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm)
real w_area_cloud_frac(dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz )

! Fraction of the grid cell volumes occupied by the cloud
real volume_cloud(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
! Fraction of the grid cell volumes occupied by the cloud
! in the previous time step
real old_volume_cloud(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)

! Tracer fluxes through the cloud surface
real fluxes_cloud(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 4)

! Tracer entrainment/detrainment through the cloud surface
real E_tetra_cloud(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 4)
real D_tetra_cloud(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 4)


!----------------------------------------------------------
! Romps variables

real entrain_vars(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 4)
real entrain_vars_old(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 4)

real core_values(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 3)
logical core_mask(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real old_core_values(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 3)
logical old_core_mask(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)

real mean_romps_core(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 4)
real E_romps_core(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 4)
real D_romps_core(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 4)

real activity_core(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 4)

real adjacent_n_core(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real active_cells_core(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real active_cells_old_core(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
logical activity_mask_core(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)


real cloud_values(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 1)
logical cloud_mask(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real old_cloud_values(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 1)
logical old_cloud_mask(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)

real mean_romps_cloud(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 4)
real E_romps_cloud(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 4)
real D_romps_cloud(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 4)

real activity_cloud(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 4)

real adjacent_n_cloud(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real active_cells_cloud(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real active_cells_old_cloud(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
logical activity_mask_cloud(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)


!----------------------------------------------------------

CONTAINS

!=============================

  subroutine init_entrainment()
    entrain_vars(:, :, :, 1) = 1.

    call interpolate_surface()

    cloud_values(:, :, :, 1) = qt_minus_qsat
    cloud_mask = all(cloud_values > 0., dim=4)

    core_values(:, :, :, 1) = qt_minus_qsat
    core_values(:, :, :, 2) = delta_thetav
    core_values(:, :, :, 3) = w_on_rho
    core_mask = all(core_values > 0., dim=4)

    entrain_vars_old(:, :, :, 1) = 1.
    entrain_vars_old(:, :, :, 2) = micro_field(:, :, :, 1)
    entrain_vars_old(:, :, :, 3) = t
    entrain_vars_old(:, :, :, 4) = w_on_rho

    call calculate_advective_momentum_fluxes()

    u_w_adv_flux_store(:, :, :, nb) = u_w_adv_flux_store(:, :, :, na)
    u_w_adv_flux_store(:, :, :, nc) = u_w_adv_flux_store(:, :, :, na)
    v_w_adv_flux_store(:, :, :, nb) = v_w_adv_flux_store(:, :, :, na)
    v_w_adv_flux_store(:, :, :, nc) = v_w_adv_flux_store(:, :, :, na)
    w_w_adv_flux_store(:, :, :, nb) = w_w_adv_flux_store(:, :, :, na)
    w_w_adv_flux_store(:, :, :, nc) = w_w_adv_flux_store(:, :, :, na)

    call calculate_diffusive_momentum_fluxes()

    u_w_dif_flux_store(:, :, :, nb) = u_w_dif_flux_store(:, :, :, na)
    u_w_dif_flux_store(:, :, :, nc) = u_w_dif_flux_store(:, :, :, na)
    v_w_dif_flux_store(:, :, :, nb) = v_w_dif_flux_store(:, :, :, na)
    v_w_dif_flux_store(:, :, :, nc) = v_w_dif_flux_store(:, :, :, na)
    w_w_dif_flux_store(:, :, :, nb) = w_w_dif_flux_store(:, :, :, na)
    w_w_dif_flux_store(:, :, :, nc) = w_w_dif_flux_store(:, :, :, na)

    d_adv = 0.
    d_dif = 0.

    adjacent_n_cloud = 0.
    activity_cloud = 0.
    mean_romps_cloud = 0.
    active_cells_cloud = 0.
    activity_mask_cloud = .false.

    adjacent_n_core = 0.
    activity_core = 0.
    mean_romps_core = 0.
    active_cells_core = 0.
    activity_mask_core = .false.

    where (cloud_mask)
      active_cells_cloud= 1.
    end where

    where (core_mask)
      active_cells_core = 1.
    end where

  end subroutine init_entrainment

!======================

  subroutine entrainment_advective_fluxes()
    integer n, i, j, k

    ! determine which cells are active
    call calculate_activity_mask(active_cells_cloud, &
                                 activity_mask_cloud, &
                                 adjacent_n_cloud)

    call calculate_activity_mask(active_cells_core, &
                                 activity_mask_core, &
                                 adjacent_n_core)

    ! u, v, w are courant numbers at this point in the main loop

    ! Find the u-, v-, and w- fluxes of each tracer and the total advective
    ! change d_adv
    do n=1,3
      call tracer_advective_flux3D(entrain_vars(:, :, :, n), &
                                   u1, v1, w1, &
                                   rho, rhow, &
                                   u_adv_fluxes(:, :, :, n), &
                                   v_adv_fluxes(:, :, :, n), &
                                   w_adv_fluxes(:, :, :, n), &
                                   d_adv(:, :, :, n))
    end do

    ! Find the u-, v-, and w- fluxes of momentum accounting for the Adams-
    ! Bashforth scheme.

    u_adv_fluxes(:, :, :, 4) = (at*u_w_adv_flux_store(:, :, :, na) &
                              + bt*u_w_adv_flux_store(:, :, :, nb) &
                              + ct*u_w_adv_flux_store(:, :, :, nc))
    v_adv_fluxes(:, :, :, 4) = (at*v_w_adv_flux_store(:, :, :, na) &
                              + bt*v_w_adv_flux_store(:, :, :, nb) &
                              + ct*v_w_adv_flux_store(:, :, :, nc))
    w_adv_fluxes(:, :, :, 4) = (at*w_w_adv_flux_store(:, :, :, na) &
                              + bt*w_w_adv_flux_store(:, :, :, nb) &
                              + ct*w_w_adv_flux_store(:, :, :, nc))

    ! Find the net advective momentum change.
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          d_adv(i, j, k, 4) = -(u_adv_fluxes(i+1, j,   k,   4)/dx &
                              - u_adv_fluxes(i  , j,   k,   4)/dx &
                              + v_adv_fluxes(i  , j+1, k,   4)/dy &
                              - v_adv_fluxes(i  , j,   k,   4)/dy &
                              + w_adv_fluxes(i  , j,   k+1, 4)/(dz*adz(k)) &
                              - w_adv_fluxes(i  , j,   k,   4)/(dz*adz(k)))
        end do
      end do
    end do

  end subroutine entrainment_advective_fluxes

!======================

  subroutine entrainment_diffusive_fluxes()
    integer i, j, k
    ! u, v, w are courant numbers at this point in the main loop

    ! Add eddy diffusion of tracers to the u-, v-, and w- fluxes.
    call tracer_diffusive_flux3D(entrain_vars(:, :, :, 2), &
                                 fluxbq, fluxtq, &
                                 tkh, &
                                 rho, rhow, &
                                 u_dif_fluxes(:, :, :, 2), &
                                 v_dif_fluxes(:, :, :, 2), &
                                 w_dif_fluxes(:, :, :, 2), &
                                 d_dif(:, :, :, 2))

    call tracer_diffusive_flux3D(entrain_vars(:, :, :, 3), &
                                 fluxbt, fluxtt, &
                                 tkh, &
                                 rho, rhow, &
                                 u_dif_fluxes(:, :, :, 3), &
                                 v_dif_fluxes(:, :, :, 3), &
                                 w_dif_fluxes(:, :, :, 3), &
                                 d_dif(:, :, :, 3))

    ! Add eddy diffusion of momentum to the u-, v-, and w- fluxes.
    u_dif_fluxes(:, :, :, 4) = (at*u_w_dif_flux_store(:, :, :, na) &
                              + bt*u_w_dif_flux_store(:, :, :, nb) &
                              + ct*u_w_dif_flux_store(:, :, :, nc))
    v_dif_fluxes(:, :, :, 4) = (at*v_w_dif_flux_store(:, :, :, na) &
                              + bt*v_w_dif_flux_store(:, :, :, nb) &
                              + ct*v_w_dif_flux_store(:, :, :, nc))
    w_dif_fluxes(:, :, :, 4) = (at*w_w_dif_flux_store(:, :, :, na) &
                              + bt*w_w_dif_flux_store(:, :, :, nb) &
                              + ct*w_w_dif_flux_store(:, :, :, nc))

    ! Find the net diffusive momentum change.
    ! TODO: Remove dx/dy/dz
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          d_dif(i, j, k, 4) = -(u_dif_fluxes(i+1, j,   k,   4)/dx &
                              - u_dif_fluxes(i,   j,   k,   4)/dx &
                              + v_dif_fluxes(i,   j+1, k,   4)/dy &
                              - v_dif_fluxes(i,   j,   k,   4)/dy &
                              + w_dif_fluxes(i,   j,   k+1, 4)/(dz*adz(k)) &
                              - w_dif_fluxes(i,   j,   k,   4)/(dz*adz(k)))
        end do
      end do
    end do


  end subroutine entrainment_diffusive_fluxes

!======================

  subroutine calculate_entrainment_fluxes()
    integer n
    ! u, v, w are courant numbers at this point in the main loop

    do n=1,4
      call calculate_fluxes_through_entrainment_surface( &
              u_adv_fluxes(:, :, :, n) + u_dif_fluxes(:, :, :, n), &
              v_adv_fluxes(:, :, :, n) + v_dif_fluxes(:, :, :, n), &
              w_adv_fluxes(:, :, :, n) + w_dif_fluxes(:, :, :, n), &
              u_area_cloud_frac, &
              v_area_cloud_frac, &
              w_area_cloud_frac, &
              volume_cloud, &
              fluxes_cloud(:, :, :, n))

      call calculate_romps_fluxes(u_adv_fluxes(:, :, :, n), &
                                  v_adv_fluxes(:, :, :, n), &
                                  w_adv_fluxes(:, :, :, n), &
                                  activity_cloud(:, :, :, n), &
                                  active_cells_cloud, &
                                  activity_mask_cloud)

      call calculate_fluxes_through_entrainment_surface( &
              u_adv_fluxes(:, :, :, n) + u_dif_fluxes(:, :, :, n), &
              v_adv_fluxes(:, :, :, n) + v_dif_fluxes(:, :, :, n), &
              w_adv_fluxes(:, :, :, n) + w_dif_fluxes(:, :, :, n), &
              u_area_core_frac, &
              v_area_core_frac, &
              w_area_core_frac, &
              volume_core, &
              fluxes_core(:, :, :, n))

      call calculate_romps_fluxes(u_adv_fluxes(:, :, :, n), &
                                  v_adv_fluxes(:, :, :, n), &
                                  w_adv_fluxes(:, :, :, n), &
                                  activity_core(:, :, :, n), &
                                  active_cells_core, &
                                  activity_mask_core)
    end do

  end subroutine calculate_entrainment_fluxes

!======================

  subroutine calculate_entrainment_rates()
    integer n

    call interpolate_surface()
    call calc_flux_forcing()

    do n=1,4
      ! Cloud calculation
      call calculate_E_D(volume_cloud, &
                         old_volume_cloud, &
                         entrain_vars(:, :, :, n), &
                         entrain_vars_old(:, :, :, n), &
                         fluxes_cloud(:, :, :, n), &
                         E_tetra_cloud(:, :, :, n), &
                         D_tetra_cloud(:, :, :, n))

      call calculate_romps_entrainment(active_cells_cloud, &
                                        active_cells_old_cloud, &
                                        entrain_vars(:, :, :, n), &
                                        entrain_vars_old(:, :, :, n), &
                                        E_romps_cloud(:, :, :, n), &
                                        D_romps_cloud(:, :, :, n), &
                                        adjacent_n_cloud, &
                                        activity_mask_cloud, &
                                        mean_romps_cloud(:, :, :, n), &
                                        activity_cloud(:, :, :, n))

      ! Core calculation
      call calculate_E_D(volume_core, &
                         old_volume_core, &
                         entrain_vars(:, :, :, n), &
                         entrain_vars_old(:, :, :, n), &
                         fluxes_core(:, :, :, n), &
                         E_tetra_core(:, :, :, n), &
                         D_tetra_core(:, :, :, n))

      call calculate_romps_entrainment(active_cells_core, &
                                        active_cells_old_core, &
                                        entrain_vars(:, :, :, n), &
                                        entrain_vars_old(:, :, :, n), &
                                        E_romps_core(:, :, :, n), &
                                        D_romps_core(:, :, :, n), &
                                        adjacent_n_core, &
                                        activity_mask_core, &
                                        mean_romps_core(:, :, :, n), &
                                        activity_core(:, :, :, n))

      entrain_vars_old(:, :, :, n) = entrain_vars(:, :, :, n)
    end do

  end subroutine calculate_entrainment_rates

!======================

  subroutine calc_flux_forcing()

    integer n

    do n=1,4
      d_total(:, :, :, n) = + (entrain_vars(:, :, :, n) &
                            - entrain_vars_old(:, :, :, n)) / dtn
      d_force(:, :, :, n) = + d_total(:, :, :, n) &
                            - d_adv(:, :, :, n) &
                            - d_dif(:, :, :, n)
      call calculate_romps_forcing(d_force(:, :, :, n), &
                                    activity_cloud(:, :, :, n), &
                                    active_cells_cloud, &
                                    activity_mask_cloud)
      call calculate_romps_forcing(d_force(:, :, :, n), &
                                    activity_core(:, :, :, n), &
                                    active_cells_core, &
                                    activity_mask_core)
    end do

  end subroutine calc_flux_forcing

!===================================

  subroutine interpolate_surface()

    call calculate_variable_interpolants()

    old_cloud_values = cloud_values
    cloud_values(:, :, :, 1) = qt_minus_qsat
    old_cloud_mask = cloud_mask
    cloud_mask = all(cloud_values > 0., dim=4)
    old_volume_cloud = volume_cloud

    call calculate_entrainment_surface(1, cloud_values, &
                                       cloud_mask, &
                                       u_area_cloud_frac, &
                                       v_area_cloud_frac, &
                                       w_area_cloud_frac, &
                                       volume_cloud)

    active_cells_old_cloud = active_cells_cloud
    call calculate_romps_activity(cloud_mask, old_cloud_mask, &
                                  1, cloud_values, old_cloud_values, &
                                  active_cells_cloud)

    old_core_values = core_values
    core_values(:, :, :, 1) = qt_minus_qsat
    core_values(:, :, :, 2) = delta_thetav
    core_values(:, :, :, 3) = w_on_rho
    old_core_mask = core_mask
    core_mask = all(core_values > 0., dim=4)
    old_volume_core = volume_core

    call calculate_entrainment_surface(3, core_values, &
                                       core_mask, &
                                       u_area_core_frac, &
                                       v_area_core_frac, &
                                       w_area_core_frac, &
                                       volume_core)

    active_cells_old_core = active_cells_core
    call calculate_romps_activity(core_mask, old_core_mask, &
                                  3, core_values, old_core_values, &
                                  active_cells_core)

  end subroutine interpolate_surface

!=============================

  subroutine calculate_variable_interpolants()
    integer :: i, j, k
    integer :: ent_id
    real :: thetav(nx, ny, nzm)
    real :: thetav_profile(nzm)

    do k=1,nzm
      do j=1, ny
        do i=1, nx
          thetav(i, j, k) = tabs(i, j, k)*prespot(k)* &
                                 (1. + epsv*qv(i,j,k) &
                                 - qcl(i, j, k) - qci(i, j, k) &
                                 - qpl(i, j, k) - qpi(i, j, k))

          qt_minus_qsat(i, j, k) = qv(i, j, k) + qcl(i, j, k) + qci(i, j, k) &
                                 + qpl(i, j, k) + qpi(i, j, k) &
                                 - qsatw(tabs(i, j, k), pres(k))
        end do
      end do
    end do

    ! Calculate w on tracer points
    w_on_rho(dimx1_w:dimx2_w, dimy1_w:dimy2_w, :) = (w1(:, :, 2:) + w1(:, :, :nzm)) / 2

    call averageXY_MPI(thetav, 1, nx, 1, ny, nzm, thetav_profile)

    do k=1,nzm
      delta_thetav(1:nx, 1:ny, k) = thetav(:, :, k) - thetav_profile(k)
    end do

    ! Exchange adjacent thetav and qt_minus_qsat halos
    ! w comes with halo points so we don't need to call task_exchange
    ent_id = nsgs_fields + nsgs_fields_diag + nmicro_fields + ntracers

    call task_exchange(qt_minus_qsat, &
                       dimx1_s, dimx2_s, &
                       dimy1_s, dimy2_s, &
                       nzm, 1,1,1,1, ent_id+1)

    call task_exchange(delta_thetav, &
                       dimx1_s, dimx2_s, &
                       dimy1_s, dimy2_s, &
                       nzm, 1,1,1,1, ent_id+2)

    entrain_vars(:, :, :, 2) = micro_field(:, :, :, 1)
    entrain_vars(:, :, :, 3) = t
    entrain_vars(:, :, :, 4) = w_on_rho

  end subroutine calculate_variable_interpolants

!=============================

  subroutine ent_hbuf_init(namelist,deflist,unitlist,status,average_type,count, entcount)
    character(*) namelist(*), deflist(*), unitlist(*)
    integer status(*), average_type(*), count, entcount

    character*8 name
    character*80 longname
    character*10 units

    entcount = 0

    ! Cloud statistics
    name = 'ETETCLD'
    longname = 'Tetrahedral Cloud Mass Entrainment'
    units = 'kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'DTETCLD'
    longname = 'Tetrahedral Cloud Mass Detrainment'
    units = 'kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'EQTETCLD'
    longname = 'Tetrahedral Cloud Specific Humidity Entrainment'
    units = 'kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'DQTETCLD'
    longname = 'Tetrahedral Cloud Specific Humidity Detrainment'
    units = 'kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'ETTETCLD'
    longname = 'Tetrahedral Cloud LWMSE Entrainmentrainment'
    units = 'K kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'DTTETCLD'
    longname = 'Tetrahedral Cloud LWMSE Detrainment'
    units = 'K kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'EWTETCLD'
    longname = 'Tetrahedral Cloud Vertical Velocity Entrainment'
    units = 'm/s kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'DWTETCLD'
    longname = 'Tetrahedral Cloud Vertical Velocity Detrainment'
    units = 'm/s kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'VTETCLD'
    longname = 'Tetrahedral Cloud Volume'
    units = 'm3/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'MFTETCLD'
    longname = 'Tetrahedral Cloud Vertical Mass Flux'
    units = 'kg/m2/s'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'EROMCLD'
    longname = 'Romps Cloud Mass Entrainment'
    units = 'kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'DROMCLD'
    longname = 'Romps Cloud Mass Detrainment'
    units = 'kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'EQROMCLD'
    longname = 'Romps Cloud Specific Humidity Entrainment'
    units = 'kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'DQROMCLD'
    longname = 'Romps Cloud Specific Humidity Detrainment'
    units = 'kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'ETROMCLD'
    longname = 'Romps Cloud LWMSE Entrainment'
    units = 'K kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'DTROMCLD'
    longname = 'Romps Cloud LWMSE Detrainment'
    units = 'K kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'EWROMCLD'
    longname = 'Romps Cloud Vertical Velocity Entrainment'
    units = 'm/s kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'DWROMCLD'
    longname = 'Romps Cloud Vertical Velocity Detrainment'
    units = 'm/s kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    ! Core statistics
    name = 'ETETCOR'
    longname = 'Tetrahedral Core Mass Entrainment'
    units = 'kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'DTETCOR'
    longname = 'Tetrahedral Core Mass Detrainment'
    units = 'kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'EQTETCOR'
    longname = 'Tetrahedral Core Specific Humidity Entrainment'
    units = 'kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'DQTETCOR'
    longname = 'Tetrahedral Core Specific Humidity Detrainment'
    units = 'kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'ETTETCOR'
    longname = 'Tetrahedral Core LWMSE Entrainment'
    units = 'K kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'DTTETCOR'
    longname = 'Tetrahedral Core LWMSE Detrainment'
    units = 'K kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'EWTETCOR'
    longname = 'Tetrahedral Core Vertical Velocity Entrainment'
    units = 'm/s kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'DWTETCOR'
    longname = 'Tetrahedral Core Vertical Velocity Detrainment'
    units = 'm/s kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'VTETCOR'
    longname = 'Tetrahedral Core Volume'
    units = 'm3/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'MFTETCOR'
    longname = 'Tetrahedral Core Vertical Mass Flux'
    units = 'kg/m2/s'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'EROMCOR'
    longname = 'Romps Core Mass Entrainment'
    units = 'kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'DROMCOR'
    longname = 'Romps Core Mass Detrainment'
    units = 'kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'EQROMCOR'
    longname = 'Romps Core Specific Humidity Entrainment'
    units = 'kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'DQROMCOR'
    longname = 'Romps Core Specific Humidity Detrainment'
    units = 'kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'ETROMCOR'
    longname = 'Romps Core LWMSE Entrainment'
    units = 'K kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'DTROMCOR'
    longname = 'Romps Core LWMSE Detrainment'
    units = 'K kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'EWROMCOR'
    longname = 'Romps Core Vertical Velocity Entrainment'
    units = 'm/s kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    name = 'DWROMCOR'
    longname = 'Romps Core Vertical Velocity Detrainment'
    units = 'm/s kg/s/m3'
    call add_to_namelist(count, entcount, name, longname, units, 0)

    if(masterproc) then
        write(*,*) 'Added ', entcount, ' arrays to statistics for entrainment'
    end if

  end subroutine ent_hbuf_init

  end module entrainment

