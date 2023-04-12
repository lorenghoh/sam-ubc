
! Monin-Obukhov Similarity for land
! Coded using description document  LSM4 in CESM  
! (C) Marat Khairoutdinov, 2016

subroutine transfer_coef(pref, ts, th, qh, qs, uh, vh, h, z0, disp, iind,jind)

use slm_vars, only : rgas, cp, pres0,pii, epsv,DBL,xsi, mom_trans_coef, heat_trans_coef, ustar, tstar, RiB, r_a, vel_m
use grid, only: dtn
implicit none

! Input:

real pref  ! reference level pressure (mb)
real th   ! temperature at height h, K
real (kind=DBL) :: ts   ! temperature at z0, K
real qh   ! vapor at height h, g/g
real (kind=DBL) :: qs   ! saturated vapor at z0, g/g
real  uh   ! zonal wind at height h, m/s
real  vh   ! merid wind at height h, m/s
real  h    ! height h, m
real (kind=DBL) :: z0   ! surfase roughness, m
real (kind=DBL) :: disp   ! displacement height
INTEGER , INTENT(IN) :: iind, jind ! indices

! Local

real (kind=DBL) :: thp, tsp
real (kind=DBL) :: r     ! bulk Richardson number
real (kind=DBL) :: zodym, zodyh, vel, qss
real (kind=DBL) :: x, x0, xm0, xh0, xsim0, xsih0, xsi1, fm, fh, error
real (kind=DBL) :: psim1, psim2, psim3, psim4, xx, yy, fm0
real (kind=DBL) :: psih1, psih2, psih3, psih4
real (kind=DBL) :: z0h, zTh, zt0
real (kind=DBL), parameter :: xsim = -1.574_DBL
real (kind=DBL), parameter :: xsih = -0.465_DBL
real (kind=DBL), parameter :: xm = sqrt(sqrt((1.-16.*xsim)))
real (kind=DBL), parameter :: xh = sqrt(sqrt((1.-16.*xsih)))
real (kind=DBL), parameter :: errormax = 0.01_DBL
real (kind=DBL), parameter :: kk = 0.35_DBL
integer, parameter :: nitermax = 10
integer niter

xx(yy) = sqrt(sqrt((1.-16.*yy)))
! unstable: -1.574<xsi<0
psim1(x,x0)=2._DBL*log((1._DBL+x)/(1._DBL+x0))+log((1._DBL+x*x)/(1._DBL+x0*x0))-2._DBL*(atan(x)-atan(x0))
psih1(x,x0)=2._DBL*log((1._DBL+x*x)/(1._DBL+x0*x0))
! very unstable:  xsi < -1.574
psim2(xsi,xsim0,xm0) = log(xsim/xsim0)-psim1(xm,xm0)+1.14_DBL*((-xsi)**0.3333_DBL-(-xsim)**0.3333_DBL)
psih2(xsi,xsih0,xh0) = log(xsih/xsih0)-psih1(xh,xh0)+0.8_DBL*((-xsi)**0.3333_DBL-(-xsih)**0.3333_DBL)
! stable: 0 < xsi < 1
psim3(xsi,xsim0) = -5._DBL*(xsi-xsim0)
psih3(xsi,xsih0) = -5._DBL*(xsi-xsih0)
! very stable: 0 < xsi < 1
psim4(xsi,xsim0) = log(xsi**5/xsim0)+5._DBL*(1.-xsim0)+xsi-1._DBL
psih4(xsi,xsih0) = log(xsi**5/xsih0)+5._DBL*(1.-xsih0)+xsi-1._DBL

tsp = ts*(1000./pres0)**(rgas/cp)
thp = th*(1000./pres0)**(rgas/cp)

! Add additional velocity depending on the stratification
if((thp-tsp).ge.0._DBL) then
        vel = sqrt(DBLE(uh)**2+DBLE(vh)**2+0.1_DBL**2)
else
        vel = sqrt(DBLE(uh)**2+DBLE(vh)**2+1.0_DBL)
END IF
r=9.81_DBL/tsp*(thp*(1._DBL+epsv*DBLE(qh))-tsp*(1._DBL+epsv*qs))*(DBLE(h)-disp)/(vel**2) ! bulk richardson number
r = max(-10._DBL,min(r, 0.19_DBL))

! initial guess:

z0h=z0/(h-disp)
!zt0 = 0.135*z0  ! roughness length for scalars
! Zeng and Dickinson (1998)
!       zt0 = z0*exp(-1.*0.13*(ustar(i,j)*z0sfc/(1.5e-5))**0.45)
! Yang (2008)
zt0 = max(0.0001_DBL,(70._DBL*1.5e-5_DBL/ustar(iind,jind))&
        *exp(-7.2_DBL*sqrt(ustar(iind,jind))*((abs(tstar(iind,jind)))**0.25_DBL)))
zTh=zt0/(h-disp)   ! (assume ln(z0/z0t) = 2.)
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
! roughness. Basically, make the maximum slowdown of the velocity not bigger than
! 50% in one timestep.
fm0 = sqrt(kk**2*vel*dtn/0.5/h)
fm = max(fm0,fm)
fh = max(fh,fh/fm*fm0)


! drag coefficient C_D = k**2/fm**2
! heat transfer coefficient C_H = k**2/fm/fh
!vel(i,j) = sqrt(uh(i,j)**2 + vh(i,j)**2)
mom_trans_coef = kk**2/fm**2
heat_trans_coef = kk**2/fm/fh
ustar(iind,jind) = sqrt(mom_trans_coef)*vel
! set ustar > 0.2 to avoid too calm conditions at night for the turbulent transfer
ustar(iind,jind) = max(0.2_DBL,ustar(iind,jind))

! aerodynamic resistance between surface and reference level
r_a = fh/kk/ustar(iind,jind)

vel_m = vel
RiB = r

return
end
