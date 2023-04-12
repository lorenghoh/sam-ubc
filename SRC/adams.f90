
subroutine adams

!       Adams-Bashforth scheme

use vars
use params, only: dowallx, dowally

implicit none

real dtdx, dtdy, dtdz, rhox, rhoy, rhoz , a1, a2
integer i,j,k

call t_startf('adams')

u1(:,:,:) = u(:,:,:)
v1(:,:,:) = v(:,:,:)
w1(:,:,:) = w(:,:,:)

do k=1,nzm
  do j=1,ny
   do i=1,nx
     u(i,j,k) = u1(i,j,k)+dt3(na) &
              *(at*dudt(i,j,k,na)+bt*dudt(i,j,k,nb)+ct*dudt(i,j,k,nc))
     v(i,j,k) = v1(i,j,k)+dt3(na) &
              *(at*dvdt(i,j,k,na)+bt*dvdt(i,j,k,nb)+ct*dvdt(i,j,k,nc))
     w(i,j,k) = w1(i,j,k)+dt3(na) &
              *(at*dwdt(i,j,k,na)+bt*dwdt(i,j,k,nb)+ct*dwdt(i,j,k,nc))
   end do
  end do
end do

! Exchange boundaries:

call boundaries(1)

! compute time averaged velocties for second-order advection of scalars:

dtdx = dtn/dx
dtdy = dtn/dy
dtdz = dtn/dz
a1 = 0.5
a2 = 0.5
if(nstep.eq.1) then
 a1 = 1.
 a2 = 0.
end if

do k=1,nzm
  rhox = rho(k)*dtdx
  do j=dimy1_u,dimy2_u
   do i=dimx1_u,dimx2_u
     u1(i,j,k) = (a1*u(i,j,k)+a2*u1(i,j,k))*rhox
   end do
  end do
end do
do k=1,nzm
  rhoy = rho(k)*dtdy
  do j=dimy1_v,dimy2_v
   do i=dimx1_v,dimx2_v
     v1(i,j,k) = (a1*v(i,j,k)+a2*v1(i,j,k))*rhoy
   end do
  end do
end do
do k=1,nzm
  do j=1,ny
   do i=1,nx
     misc(i,j,k) = (a1*w(i,j,k)+a2*w1(i,j,k))
   end do
  end do
end do
do k=1,nzm
  rhoz = rhow(k)*dtdz
  do j=dimy1_w,dimy2_w
   do i=dimx1_w,dimx2_w
     w1(i,j,k) = (a1*w(i,j,k)+a2*w1(i,j,k))*rhoz
   end do
  end do
end do


call t_stopf ('adams')

end subroutine adams

	
