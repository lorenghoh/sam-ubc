subroutine kurant

use vars
use sgs, only: kurant_sgs
use params, only: ncycle_max

implicit none

integer i, j, k
real cfl,cflz,cflh,cflz1,cflh1,coef,idx,idy,idz
real buf1(4),buf2(4)

call t_startf ('kurant')

ncycle = 1
cfl = 0.
cflz = 0.
cflh = 0.
idy = YES3D*dt/dy
idx = dt/dx
do k = 1,nzm
 idz = dt/(dz*adzw(k))
 do j=1,ny
  do i=1,nx
   cflz1 = abs(w(i,j,k))*idz
   cflh1 = sqrt((u(i,j,k)*idx)**2+(v(i,j,k)*idy)**2)
   cfl = max(cfl,sqrt(cflh1**2+cflz1**2))
   cflh = max(cflh,cflh1)
   cflz = max(cflz,cflz1)
  end do
 end do
end do
w_max=max(w_max,maxval(w(1:nx,1:ny,1:nz)))
u_max=max(u_max,sqrt(maxval(u(1:nx,1:ny,1:nzm)**2+YES3D*v(1:nx,1:ny,1:nzm)**2)))

call kurant_sgs(cfl_sgs)

cfl_adv=cfl
cfl_advh=cflh
cfl_advz=cflz
if(dompi) then
  buf1(1)=cfl_adv
  buf1(2)=cfl_advh
  buf1(3)=cfl_advz
  buf1(4)=cfl_sgs
  call task_max_real(buf1,buf2,4)
  cfl_adv=buf2(1)
  cfl_advh=buf2(2)
  cfl_advz=buf2(3)
  cfl_sgs=buf2(4)
end if
ncycle = max(1,ceiling(cfl_adv/0.7),ceiling(cfl_sgs/0.1))
if(ncycle.gt.ncycle_max) then
   if(masterproc) print *,'the number of cycles exceeded ', ncycle_max
   call stepout(-1)
   call write_fields2D()
   call write_fields3D()
   call task_abort()
end if

call t_stopf ('kurant')

end subroutine kurant   


