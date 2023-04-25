subroutine tracer_diffusive_flux3d(input_field, fluxb, fluxt, tkh, rho, rhow, uflux, vflux, wflux, d_dif)

use grid
use sgs, only: dimx1_d, dimx2_d, dimy1_d, dimy2_d, grdf_x, grdf_y, grdf_z
use params

implicit none

! input
real input_field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real tkh(dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm)
real fluxb(nx, ny)
real fluxt(nx, ny)
real rho(nzm)
real rhow(nz)
real uflux(dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm)
real vflux(dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm)
real wflux(dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz )
real d_dif(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)

! local
real field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! scalar
real flx(0:nx,0:ny,0:nzm)
real rdx2,rdy2,rdz2,rdz,rdx5,rdy5,rdz5,tmp
real dxy,dxz,dyx,dyz,dzx,dzy,tkx,tky,tkz,rhoi
integer i,j,k,ib,ic,jb,jc,kc,kb

if(.not.dosgs) return

field = input_field

rdx2 = 1./(dx*dx)
rdy2 = 1./(dy*dy)
rdz2 = 1./(dz*dz)
rdz = 1./dz
dxy = dx/dy
dxz = dx/dz
dyx = dy/dx
dyz = dy/dz
dzx = dz/dx
dzy = dz/dy

!-----------------------------------------
if(dowallx) then

  if(mod(rank,nsubdomains_x).eq.0) then
    do k=1,nzm
      do j=1,ny
        field(0,j,k) = field(1,j,k)
      end do
    end do
  end if
  if(mod(rank, nsubdomains_x) .eq. nsubdomains_x-1) then
    do k=1,nzm
      do j=1,ny
        field(nx+1, j, k) = field(nx, j, k)
      end do
    end do
  end if

end if

if(dowally) then

  if(rank .lt. nsubdomains_x) then
    do k=1,nzm
       do i=1,nx
         field(i, 1-YES3D, k) = field(i, 1, k)
       end do
    end do
  end if
  if (rank .gt. nsubdomains-nsubdomains_x-1) then
    do k=1,nzm
       do i=1,ny
         field(i, ny+YES3D, k) = field(i, ny, k)
       end do
    end do
  end if

  call task_rank_to_index(rank, ib, jb)
  if (jb.eq.0) then
    do k=1,nzm
      do i=1,nx
        field(i, 1-YES3D, k) = field(i, 1, k)
      end do
    end do
  end if
  if (jb.eq.nsubdomains_y-1) then
    do k=1,nzm
      do i=1,nx
        field(i, ny+YES3D, k) = field(i, ny, k)
      end do
    end do
  end if

end if

!-----------------------------------------


!  Horizontal diffusion:

do k=1,nzm
	
  rdx5 = 0.5*rdx2 * grdf_x(k)
  rdy5 = 0.5*rdy2 * grdf_y(k)

  do j=1,ny
    do i=0,nx
      ic = i+1
      tkx = rdx5*(tkh(i,j,k) + tkh(ic,j,k)) 	
      flx(i, j, k) = -tkx*(field(ic, j, k) - field(i, j, k))
      uflux(i+1, j, k) = flx(i, j, k)
    end do
  end do 

  do j=0,ny
    jc=j+1
    do i=1,nx
      tky = rdy5*(tkh(i,j,k) + tkh(i,jc,k)) 	
      flx(i,j,k) = -tky*(field(i,jc,k) - field(i,j,k))
      vflux(i, j+1, k) = flx(i, j, k)
    end do
  end do 

end do ! k


!  Vertical diffusion:

tmp = 1./adzw(nz)
do j=1,ny
  do i=1,nx	
    flx(i, j, 0) = fluxb(i, j)*rdz*rhow(1)
    flx(i, j, nzm) = fluxt(i, j)*rdz*tmp*rhow(nz)
  end do
end do

do k=1,nzm-1
  kc=k+1	
  rhoi = rhow(kc)/adzw(kc)
  rdz5 = 0.5*rdz2 * grdf_z(k)
  do j=1,ny
    do i=1,nx
      tkz = rdz5*(tkh(i,j,k) + tkh(i,j,kc))
      flx(i, j, k) = -tkz*(field(i,j,kc) - field(i,j,k))*rhoi
    
      wflux(i, j, k+1) = flx(i, j, k)
    end do 
  end do
end do

do k=1,nzm
  rhoi = 1./(adz(k)*rho(k))
  do j=1,ny
    do i=1,nx
      d_dif(i, j, k) = -(uflux(i+1, j, k) - uflux(i, j, k) &
                       + vflux(i, j+1, k) - vflux(i, j, k) &
                       + wflux(i, j, k+1) - wflux(i, j, k)) * rhoi * dtn
    end do
  end do
end do
              

end subroutine tracer_diffusive_flux3d
