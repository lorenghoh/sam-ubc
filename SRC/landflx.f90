
! Monin-Obukhov Similarity for land
! Coded using description document  LSM4 in CESM  
! (C) Marat Khairoutdinov, Oct 2016

subroutine landflx(th, ts, qh, qs, uh, vh, h, z0, shf, lhf, taux, tauy)

use grid, only: dtn

implicit none

! Input:

real th   ! pot. temperature at height h, K
real ts   ! pot. Temperature at z0, K
real qh   ! vapor at height h, g/g
real qs   ! saturated vapor at z0, g/g
real uh   ! zonal wind at height h, m/s
real vh   ! merid wind at height h, m/s
real h    ! height h, m
real z0   ! surfase roughness, m

! Output:

real shf   ! sensible heat flux (K m/s)
real lhf   ! latent heat flux (m/s)
real taux  ! zonal surface stress (m2/s2)
real tauy  ! merid surface stress (m2/s2)

! Local

real r     ! bulk Richardson number
real pii, zodym, zodyh, vel
real ustar, tstar, qstar
real x, x0, xm0, xh0, xsi, xsim0, xsih0, xsi1, fm, fh, error
real psim1, psim2, psim3, psim4, xx, yy, fm0
real psih1, psih2, psih3, psih4
real z0h, zTh, zt0
real, parameter :: xsim = -1.574
real, parameter :: xsih = -0.465
real, parameter :: xm = sqrt(sqrt((1.-16.*xsim)))
real, parameter :: xh = sqrt(sqrt((1.-16.*xsih)))
real, parameter :: errormax = 0.01
integer, parameter :: nitermax = 10
integer niter

xx(yy) = sqrt(sqrt((1.-16.*yy)))
! unstable: -1.574<xsi<0
psim1(x,x0)=2.*log((1.+x)/(1+x0))+log((1.+x*x)/(1.+x0*x0))-2.*(atan(x)-atan(x0))
psih1(x,x0)=2.*log((1.+x*x)/(1.+x0*x0))
! very unstable:  xsi < -1.574
psim2(xsi,xsim0,xm0) = log(xsim/xsim0)-psim1(xm,xm0)+1.14*((-xsi)**0.3333-(-xsim)**0.3333)
psih2(xsi,xsih0,xh0) = log(xsih/xsih0)-psih1(xh,xh0)+0.8*((-xsi)**0.3333-(-xsih)**0.3333)
! stable: 0 < xsi < 1
psim3(xsi,xsim0) = -5.*(xsi-xsim0)
psih3(xsi,xsih0) = -5.*(xsi-xsih0)
! very stable: 0 < xsi < 1
psim4(xsi,xsim0) = log(xsi**5/xsim0)+5.*(1.-xsim0)+xsi-1.
psih4(xsi,xsih0) = log(xsi**5/xsih0)+5.*(1.-xsih0)+xsi-1.

vel = sqrt(max(0.5,uh**2+vh**2))
r=max(-10.,min(0.19,9.81*((th-ts)/ts*(1.+0.61*qh)+0.61*(qh-qs))*h/vel**2))

! initial guess:

z0h=z0/h
zt0 = 0.135*z0  ! roughness length for scalars
zTh=zt0/h   ! (assume ln(z0/z0t) = 2.)
zodym = log(1./z0h)
zodyh = log(1./zTh)

if(r.gt.0.) then
 xsi = r*zodym/(1.-5.*r)
else
 xsi = r*zodym
end if

niter = 0
error = 1000
do while (error.gt.errormax.and.niter.lt.nitermax)

xsi1 = xsi
niter = niter + 1
xsim0 = z0h*xsi
xsih0 = zTh*xsi 

if(xsi.lt.-0.01) then
  if(xsi.ge.xsim) then
    x = xx(xsi)
    x0 = xx(xsim0)
    fm = zodym-psim1(x,x0)
  else
    xm0 = xx(xsim0)
    fm = psim2(xsi,xsim0,xm0)
  end if
  if(xsi.ge.xsih) then
    x = xx(xsi)
    x0 = xx(xsih0)
    fh = zodyh-psih1(x,x0)
  else
    xh0 = xx(xsih0)
    fh = psih2(xsi,xsih0,xh0)
  end if
elseif(xsi.gt.0.01) then
  if(xsi.le.1.) then
    fm = zodym-psim3(xsi,xsim0)
    fh = zodyh-psih3(xsi,xsih0)
  else
    fm = psim4(xsi,xsim0)
    fh = psih4(xsi,xsih0)
  end if
else
  fm = zodym
  fh = zodyh
end if
xsi = r*fm*fm/fh
error = abs(xsi-xsi1)

end do

! limit fh and fh to avoid too large fluxes especially over large surface
! roughness. Basically, make the maximum slowdown of the velocity not bigger
! than ! 50% in one timestep.

fm0 = sqrt(0.4**2*vel*dtn/0.5/h)
fm = max(fm0,fm)
fh = max(fh,fh/fm*fm0)

shf = 0.4**2*vel/(fm*fh)*(ts-th)
lhf = 0.4**2*vel/(fm*fh)*(qs-qh)
taux=-0.4**2*vel/(fm*fm)*uh
tauy=-0.4**2*vel/(fm*fm)*vh

return
end
