
subroutine calculate_advective_momentum_fluxes()

use grid
use vars
use entrainment

implicit none

real, dimension(0:nx, 1-YES3D:ny, nzm) :: fuw, fvw, fww
real :: dx25, dy25, dz25
integer :: i, j, k, ib, jb, kb, ic, jc, kc, kcu

dx25 = 0.25 / dx
dy25 = 0.25 / dy
dz25 = 0.25 / dz

w_w_adv_flux_store(:, :, nz, na) = 0.

do k = 1, nzm
  kc = k+1
  kcu = min(kc, nzm)
  do j = 1, ny+1
    jc = j+1
    do i = 1, nx+1
      ic = i+1
      fuw(i, j, k) = dx25 * (u(ic, j, k)*rho(k)*adz(k) + u(ic, j, kcu)*rho(kcu)*adz(kcu)) * (w(i, j, kc) + w(ic, j, kc))

      fvw(i, j, k) = dy25 * (v(i, jc, k)*rho(k)*adz(k) + v(i, jc, kcu)*rho(kcu)*adz(kcu)) * (w(i, j, kc) + w(i, jc, kc))

      fww(i, j, k) = dz25 * (w(i, j, kc)*rhow(kc) + w(i, j, k)*rhow(k)) * (w(i, j, kc) + w(i, j, k))
    end do
  end do
end do

do k = 1, nzm
  kb = k-1
  kc = k+1
  kcu = min(kc, nzm)
  do j = 1, ny
    jb = j-1
    do i = 1, nx
      ib = i-1
      u_w_adv_flux_store(i, j, k, na) = (fuw(ib, j, kcu) + fuw(ib, j, k)) * dx / 2.
      v_w_adv_flux_store(i, j, k, na) = (fvw(i, jb, kcu) + fvw(i, jb, k)) * dy / 2.
      w_w_adv_flux_store(i, j, k, na) = (fww(i, j, k) + fww(i, j, kb)) * dz / 2.
    end do
  end do
end do

w_w_adv_flux_store(:, :, 1, na) = 0.

end subroutine calculate_advective_momentum_fluxes


