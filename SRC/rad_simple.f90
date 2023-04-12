
subroutine rad_simple
 	
!	Simple Interactive Radiation
!  Coded by Ping Zhu for DycomsII LES intercomparison


use grid
use vars
use params
use rad, only: qrad
implicit none
	
real f0, xk, coef, coef1
integer i,k

!pzhu
real qq1,qq2,hfact,flux(nz),FTHRL(nzm) !bloss: flux(nz) rather than flux(nzm)
integer j,kk,itop,item3
!pzhu-end

!bloss
real deltaq(nzm), cpmassl(nzm), qzinf, qzeroz
real, parameter :: cp_spec = 1015. ! from DycomsII LES intercomparison specs


if(.not.dolongwave) return

coef=70.
coef1=22.
f0=3.75e-6
xk=85.

do k=1,nzm
radlwdn(k) =0.
radqrlw(k) =0.
enddo

do k = 1,nzm
   cpmassl(k) = cp_spec*rho(k)*dz*adz(k) ! thermal mass of model layer.
end do

do i=1,nx
do j=1,ny

   ! search for inversion height (highest level where qt=8 g/kg) 
   !  and compute optical depth increments.
   itop = 1
   deltaq = 0.
   qzinf = 0. ! holds accumulated optical depth between z and domain top
   qzeroz = 0.! holds accumulated optical depth between surface and z
   do k = 1,nzm
      ! optical depth only includes that due to liquid water
      if(qcl(i,j,k)+qci(i,j,k).gt.0.) deltaq(k) = xk*rho(k)*(qcl(i,j,k)+qci(i,j,k))*adz(k)*dz

      ! inversion height is top of highest layer w/q>8g/kg.
      if(qv(i,j,k)+qcl(i,j,k)+qci(i,j,k).gt.0.008) itop=max(itop,k+1) ! note zi(k+1) is inversion hgt

      ! accumulate optical depth in qzinf
      qzinf = qzinf + deltaq(k) 
   end do
   ! note: qzinf now holds total optical depth of column (due to cloud).
   !       qzeroz is initialized to zero since first level is at surface.

   ! compute net upward longwave flux up to inversion height.
   do k = 1,itop
      flux(k) = coef*exp(-qzinf) + coef1*exp(-qzeroz)
      qzinf  = qzinf  - deltaq(k)
      qzeroz = qzeroz + deltaq(k)
   end do

   ! compute net upward longwave flux above inversion height.
   ! this includes correction for clearsky fluxes which balances
   ! the prescribed subsidence heating above the inversion.
   do k = itop+1,nzm
      flux(k) = coef*exp(-qzinf) + coef1*exp(-qzeroz) &
           + cp_spec*rhow(k)*f0*(0.25*zi(k)+0.75*zi(itop)) &
                               *(zi(k)-zi(itop))**(1./3.)
      qzinf  = qzinf  - deltaq(k)
      qzeroz = qzeroz + deltaq(k)
   end do
   flux(nz) = coef*exp(-qzinf) + coef1*exp(-qzeroz) &
        + cp_spec*rhow(nz)*f0*(0.25*zi(nz)+0.75*zi(itop)) &
                             *(zi(nz)-zi(itop))**(1./3.)

   ! note that our flux differs from the specification in that it
   ! uses the local density (rather than the interface density) in
   ! computing the correction to the flux above the inversion. 
   ! Formulating things in this way ensures a rough balance between
   ! the subsidence heating and radiative cooling above the inversion.

   ! compute radiative heating as divergence of net upward lw flux.
   do k=1,nzm
      FTHRL(k)=-(flux(k+1)-flux(k))/cpmassl(k)
      t(i,j,k) = t(i,j,k) + FTHRL(k) * dtn ! add radiative heating to sli
      radlwdn(k) = radlwdn(k) + flux(k)   ! net lw flux for stats
      radqrlw(k) = radqrlw(k) + FTHRL(k)  ! net lw heating for stats
      qrad(i,j,k) = FTHRL(k) ! store radiative heating for 3d output/stepout
   enddo
   
end do
end do


end 




