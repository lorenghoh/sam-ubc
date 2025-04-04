subroutine tracer_advective_flux3D (f_in, &
                                    u, v, w, &
                                    rho, rhow, &
                                    uflux, vflux, wflux, &
                                    d_adv)
 	
!     positively definite monotonic advection with non-oscillatory option

use grid
use params

implicit none

real, intent(in) :: f_in(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)

real, intent(in) :: u(dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm)
real, intent(in) :: v(dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm)
real, intent(in) :: w(dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz)
real, intent(in) :: rho(nzm)
real, intent(in) :: rhow(nz)

real, intent(inout) :: uflux(dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm)
real, intent(inout) :: vflux(dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm)
real, intent(inout) :: wflux(dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz)
real, intent(inout) :: d_adv(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)

real :: f(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)

real mx(0:nxp1, 0:nyp1, nzm)
real mn(0:nxp1, 0:nyp1, nzm)
real uuu(-1:nxp3, -1:nyp2, nzm)
real vvv(-1:nxp2, -1:nyp3, nzm)
real www(-1:nxp2, -1:nyp2, nz)

real eps, dd
real iadz(nzm), irho(nzm), irhow(nzm)
integer i, j, k, ic, ib, jc, jb, kc, kb
logical nonos

real x1, x2, a, b, a1, a2, y
real andiff, across, pp, pn

andiff(x1, x2, a, b) = (abs(a) - a*a*b)*0.5*(x2 - x1)
across(x1, a1, a2) = 0.03125*a1*a2*x1
pp(y) =  max(0.,y)
pn(y) = -min(0.,y)

nonos = .true.
eps = 1.e-10

f = f_in

www(:, :, nz) = 0.

!-----------------------------------------

if(nonos) then
 do k=1,nzm
  kc=min(nzm,k+1)
  kb=max(1,k-1)
  do j=0,nyp1
    jb=j-1
    jc=j+1
    do i=0,nxp1
      ib=i-1
      ic=i+1
      mx(i,j,k) = max(f(ib,j,k), f(ic,j,k), f(i,jb,k), &
                    f(i,jc,k), f(i,j,kb), f(i,j,kc), f(i,j,k))
      mn(i,j,k) = min(f(ib,j,k), f(ic,j,k), f(i,jb,k), &
                    f(i,jc,k), f(i,j,kb), f(i,j,kc), f(i,j,k))	  
    end do
  end do
 end do
end if ! nonos

do k=1,nzm
  do j=-1,nyp2
    do i=-1,nxp3
    uuu(i,j,k) = max(0.,u(i,j,k))*f(i-1,j,k)+min(0.,u(i,j,k))*f(i,j,k)
    uflux(i,j,k) = uuu(i,j,k)
    end do
  end do
end do

do k=1,nzm
  do j=-1,nyp3
    do i=-1,nxp2
    vvv(i,j,k) = max(0.,v(i,j,k))*f(i,j-1,k)+min(0.,v(i,j,k))*f(i,j,k)
    vflux(i,j,k) = vvv(i,j,k)
    end do
  end do
end do

do k=1,nzm
  kb=max(1,k-1)
  do j=-1,nyp2
    do i=-1,nxp2
    www(i,j,k)=max(0.,w(i,j,k))*f(i,j,kb)+min(0.,w(i,j,k))*f(i,j,k)
    wflux(i,j,k) = www(i,j,k)
    end do
  end do
end do

do k=1,nzm
  irho(k) = 1./rho(k)
  iadz(k) = 1./adz(k)
  do j=-1,nyp2
    do i=-1,nxp2
      f(i,j,k) = f(i,j,k) -(uuu(i+1,j,k)-uuu(i,j,k)+vvv(i,j+1,k)-vvv(i,j,k) &
                        +(www(i,j,k+1)-www(i,j,k))*iadz(k))*irho(k)
    end do
  end do
end do 

do k=1,nzm
  kc=min(nzm,k+1)
  kb=max(1,k-1)
  dd=2./(kc-kb)/adz(k)
  do j=0,nyp1
    jb=j-1
    jc=j+1
    do i=0,nxp2
    ib=i-1
    uuu(i,j,k) = andiff(f(ib,j,k),f(i,j,k),u(i,j,k),irho(k)) &
              -(across(f(ib,jc,k)+f(i,jc,k)-f(ib,jb,k)-f(i,jb,k), &
                    u(i,j,k), v(ib,j,k)+v(ib,jc,k)+v(i,jc,k)+v(i,j,k)) &
              +across(dd*(f(ib,j,kc)+f(i,j,kc)-f(ib,j,kb)-f(i,j,kb)), &
                    u(i,j,k), w(ib,j,k)+w(ib,j,kc)+w(i,j,k)+w(i,j,kc))) *irho(k)
    end do
  end do
end do

do k=1,nzm
  kc=min(nzm,k+1)
  kb=max(1,k-1)
  dd=2./(kc-kb)/adz(k)
  do j=0,nyp2
    jb=j-1
    do i=0,nxp1
    ib=i-1
    ic=i+1
    vvv(i,j,k) = andiff(f(i,jb,k),f(i,j,k),v(i,j,k),irho(k)) &
              -(across(f(ic,jb,k)+f(ic,j,k)-f(ib,jb,k)-f(ib,j,k), &
                    v(i,j,k), u(i,jb,k)+u(i,j,k)+u(ic,j,k)+u(ic,jb,k)) &
                +across(dd*(f(i,jb,kc)+f(i,j,kc)-f(i,jb,kb)-f(i,j,kb)), & 
                    v(i,j,k), w(i,jb,k)+w(i,j,k)+w(i,j,kc)+w(i,jb,kc))) *irho(k)
    end do
  end do
end do

do k=1,nzm
  kb=max(1,k-1)
  irhow(k)=1./(rhow(k)*adz(k))
  do j=0,nyp1
    jb=j-1
    jc=j+1
    do i=0,nxp1
    ib=i-1
    ic=i+1
    www(i,j,k) = andiff(f(i,j,kb),f(i,j,k),w(i,j,k),irhow(k)) &
              -(across(f(ic,j,kb)+f(ic,j,k)-f(ib,j,kb)-f(ib,j,k), &
                  w(i,j,k), u(i,j,kb)+u(i,j,k)+u(ic,j,k)+u(ic,j,kb)) &
              +across(f(i,jc,k)+f(i,jc,kb)-f(i,jb,k)-f(i,jb,kb), &
                  w(i,j,k), v(i,j,kb)+v(i,jc,kb)+v(i,jc,k)+v(i,j,k)))*irhow(k)
    end do
  end do
end do

www(:,:,1) = 0.

!---------- non-osscilatory option ---------------

if(nonos) then
  do k=1,nzm
    kc=min(nzm,k+1)
    kb=max(1,k-1)
    do j=0,nyp1
    jb=j-1
    jc=j+1
    do i=0,nxp1
      ib=i-1
      ic=i+1
      mx(i,j,k) = max(f(ib,j,k),f(ic,j,k),f(i,jb,k), &
                    f(i,jc,k),f(i,j,kb),f(i,j,kc),f(i,j,k),mx(i,j,k))
      mn(i,j,k) = min(f(ib,j,k),f(ic,j,k),f(i,jb,k), &
                    f(i,jc,k),f(i,j,kb),f(i,j,kc),f(i,j,k),mn(i,j,k))
    end do
    end do
  end do

  do k=1,nzm
    kc=min(nzm,k+1)
    do j=0,nyp1
      jc=j+1
      do i=0,nxp1
        ic=i+1
        mx(i,j,k) = rho(k)*(mx(i,j,k)-f(i,j,k))/ &
                          (pn(uuu(ic,j,k)) + pp(uuu(i,j,k))+ &
                            pn(vvv(i,jc,k)) + pp(vvv(i,j,k))+ &
                  iadz(k)*(pn(www(i,j,kc)) + pp(www(i,j,k)))+eps)
        mn(i,j,k) = rho(k)*(f(i,j,k)-mn(i,j,k))/ &
                          (pp(uuu(ic,j,k)) + pn(uuu(i,j,k))+ &
                            pp(vvv(i,jc,k)) + pn(vvv(i,j,k))+ &
                  iadz(k)*(pp(www(i,j,kc)) + pn(www(i,j,k)))+eps)
      end do
    end do
  end do

  do k=1,nzm
    do j=1,ny
      do i=1,nxp1
        ib=i-1
        uuu(i,j,k) = pp(uuu(i,j,k))*min(1.,mx(i,j,k), mn(ib,j,k)) &
                - pn(uuu(i,j,k))*min(1.,mx(ib,j,k),mn(i,j,k))
        uflux(i, j, k) = uflux(i, j, k) + uuu(i, j, k)
      end do
    end do
  end do

  do k=1,nzm
    do j=1,nyp1
      jb=j-1
      do i=1,nx
        vvv(i,j,k) = pp(vvv(i,j,k))*min(1.,mx(i,j,k), mn(i,jb,k)) &
                - pn(vvv(i,j,k))*min(1.,mx(i,jb,k),mn(i,j,k))
        vflux(i, j, k) = vflux(i, j, k) + vvv(i, j, k)
      end do
    end do
  end do

  do k=1,nzm
    kb=max(1,k-1)
    do j=1,ny
      do i=1,nx
        www(i,j,k) = pp(www(i,j,k))*min(1.,mx(i,j,k), mn(i,j,kb)) &
                - pn(www(i,j,k))*min(1.,mx(i,j,kb),mn(i,j,k))
        wflux(i, j, k) = wflux(i, j, k) + www(i, j, k)
      end do
    end do
  end do
endif ! nonos

uflux = uflux * dx / dtn
vflux = vflux * dy / dtn
wflux = wflux * dz / dtn

do k=1,nzm
  kc = k+1
  do j=1,ny
    do i=1,nx
      f(i, j, k) = max(0., f(i,j,k) - (uuu(i+1,j,k) - uuu(i,j,k) + vvv(i,j+1,k) - vvv(i,j,k) + (www(i,j,k+1) - www(i,j,k))*iadz(k))*irho(k))

      d_adv(i, j, k) = (f(i, j, k) - f_in(i, j, k)) / dtn
    end do
  end do
end do

end subroutine tracer_advective_flux3D
