subroutine calculate_diffusive_momentum_fluxes()
	
!        momentum tendency due to SGS diffusion

use vars
use sgs
use entrainment

implicit none

real rdx2, rdy2, rdz2, rdz, rdx25, rdy25
real rdx251, rdy251, rdz25
real dxz, dyz

integer i, j, k, ic, ib, jb, jc, kc, kcu, kb
real tkx, tky, tkz, rhoi, iadz
real fuw(0:nx, 0:ny, nz), fvw(0:nx, 0:ny, nz), fww(0:nx, 0:ny, nz)

rdx2 = 1./(dx*dx)
rdy2 = 1./(dy*dy)

rdx25 = 0.25*rdx2
rdy25 = 0.25*rdy2

dxz = dx/dz
dyz = dy/dz

do k = 1,nzm
  kc = k + 1
  kcu = min(kc, nzm)
  dxz = dx/(dz*adzw(kc))
  dyz = dy/(dz*adzw(kc))
  rdx251 = rdx25  * grdf_x(k)
  rdy251 = rdy25  * grdf_y(k)
  do j = 1,ny
    jb = j-1
    do i = 0,nx
      ic = i+1
      tkx = rdx251*(tk(i, j, k) + tk(ic, j, k) + tk(ic, j, kcu) + tk(ic, j, k)) 	
      fuw(i, j, k) = -tkx*(w(ic, j, kc) - w(i, j, kc) + (u(ic, j, kcu) - u(ic, j, k))*dxz)
    end do
  end do

  do j = 0,ny
    jc = j+1
    do i = 1,nx
      ib = i-1
      tky = rdy251*(tk(i, j, k) + tk(i, jc, k) + tk(i, j, kcu) + tk(i, jc, kcu))
      fvw(i, j, k) = -tky*(w(i, jc, kc) - w(i, j, kc) + (v(i, jc, kcu) - v(i, jc, k))*dyz)
    end do
  end do
end do
 
!-------------------------
rdz = 1./dz

do k = 1,nzm
  kc = k+1
  iadz = 1./adz(k)
  rdz2 = rdz*rdz * grdf_z(k)
  rdz25 = 0.25*rdz2
  do j = 1,ny
    jb = j-1
    do i = 1,nx
      ib = i-1
      tkz = rdz2*tk(i, j, k)
      fww(i, j, kc) = -2.*tkz*(w(i, j, kc) - w(i, j, k))*rho(k)*iadz
    end do
  end do
end do

do k = 1,nzm
  kb = k-1
  kc = k + 1
  kcu = min(kc, nzm)
  do j = 1,ny
    jb = j-1
    do i = 1,nx
      ib = i-1
      u_w_dif_flux_store(i, j, k, na) = (fuw(ib, j, kcu) + fuw(ib, j, k)) / 2.
      v_w_dif_flux_store(i, j, k, na) = (fvw(i, jb, kcu) + fvw(i, jb, k)) / 2.
      w_w_dif_flux_store(i, j, k, na) = (fww(i, j, k) + fww(i, j, kb)) / 2.
    end do
  end do
end do

end subroutine calculate_diffusive_momentum_fluxes
