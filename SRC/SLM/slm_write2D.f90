!=====================================================================
! Based on write_fields2D.f90 of SAM
! Write 2D variables into the .2Dcom output file
! Caution: 
!	nfields (number of 2D output fields) should be updated with the 
!	number of SLM 2D output fields
! 	For instance, in SAM/SRC/write_fields2D.f90, 
!	nfields = nfields + (31+(nsoil*2)+dosoilwnudging*nsoil+dosoiltnudging*nsoil) <= for the current version
!	
!=====================================================================
subroutine slm_write2D(nfields1, coef)

use slm_vars

implicit none

!!!!!!!!!!!!!!!!!!!!!!!CAUTION!!!!!!!!!!!!!!!!!!!!
!
! MAKE SURE THE NUMBER OF DIAGNOSTICS HERE MATCH THOSE IN
! THE CALLING ROUTINE, write_fields2D.f90
!
!!!!!!!!!!!!!!!!!!!!!!!CAUTION!!!!!!!!!!!!!!!!!!!!

! argument list
   integer, intent(inout) ::nfields1
   real, intent(in) :: coef
    
! local variables
   integer :: i, j, k
   real(4), dimension(nx,ny,1) :: tmp
   real(kind=DBL) :: coef_dble
   character *80 long_name
   character *8 name
   character *4 numchar ! convert soil layer number to character
   character *10 units
 
   k = 0   
   coef_dble = DBLE(coef)
!====================================================================
!! Surface momentum flux in x
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_taux_sfc(i,j)*coef_dble,4)
       s_taux_sfc(i,j)= 0.
     end do
   end do
  name='S_TAUX'
  long_name='SFC Momentum Flux X'
  units='m2/s2'
 call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
!====================================================================
!! Surface momentum flux in y
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_tauy_sfc(i,j)*coef_dble,4)
       s_tauy_sfc(i,j)= 0.
     end do
   end do
  name='S_TAUY'
  long_name='SFC Momentum Flux Y'
  units='m2/s2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
!====================================================================
! precipitation interception rate by canopy 
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_precip(i,j)*86400._DBL*coef_dble,4)
       s_precip(i,j)= 0._DBL
     end do
   end do
  name='S_PCNP'
  long_name='Precip. Interception Rate'
  units='mm/d'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
!====================================================================
! precipitation drainage rate from canopy
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_drain(i,j)*86400._DBL*coef_dble,4)
       s_drain(i,j)= 0._DBL
     end do
   end do
  name='S_DCNP'
  long_name='Precip. Drainage Rate'
  units='mm/d'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
!====================================================================
! precip rate at surface 
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_precip_sfc(i,j)*86400._DBL*coef_dble,4)
       s_precip_sfc(i,j)= 0._DBL
     end do
   end do
  name='S_PSFC'
  long_name='Precip. Rate at soil surface'
  units='mm/d'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
!====================================================================
! net  rad  on canopy 
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_net_rad(1,i,j)*coef_dble,4)
       s_net_rad(1,i,j)= 0._DBL
     end do
   end do
  name='S_NetRC'
  long_name='Net Rad on canopy'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
!====================================================================
! net  rad  on soil 
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_net_rad(2,i,j)*coef_dble,4)
       s_net_rad(2,i,j)= 0._DBL
     end do
   end do
  name='S_NetRS'
  long_name='Net Rad on soil'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
!====================================================================
!!roughness length
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=z0_sfc(i,j)
     end do
   end do
  name='S_Z0'
  long_name='Surface Roughness Length'
  units='m'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
!====================================================================
!!Vegetation height
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=ztop(i,j)
     end do
   end do
  name='S_ZTOP'
  long_name='Vegetation Height'
  units='m'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

!====================================================================
! Friction Velocity 
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_ustar(i,j)*coef_dble,4)
       s_ustar(i,j)= 0._DBL
     end do
   end do
  name='S_USTAR'
  long_name='Friction Velocity'
  units='m/s'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
!====================================================================
! IGBP Landtype
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(landtype(i,j),4)
     end do
   end do
  name='S_Landtype'
  long_name='Landtype'
  units=''
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

!====================================================================
! LAI
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_LAI(i,j)*coef_dble,4)
       s_LAI(i,j)= 0._DBL
     end do
   end do
  name='S_LAI'
  long_name='Leaf Area Index'
  units='m2/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

!====================================================================
! Basal Area Index
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=BAI(i,j)
     end do
   end do
  name='S_BAI'
  long_name='Basal Area Index'
  units='sq.f/ac'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

!====================================================================
! Canopy heat capacity
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_Cveg(i,j)*coef_dble,4)
       s_LAI(i,j)= 0._DBL
     end do
   end do
  name='S_Cveg'
  long_name='Canopy Heat Capacity'
  units='J/m2/K'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

!====================================================================
! ref. resistance, Ra 
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_r_a(i,j)*coef_dble,4)
       s_r_a(i,j)= 0._DBL
     end do
   end do
  name='S_Ra'
  long_name='Aerodynamic Resistance'
  units='s/m'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
!========================
! leaf boundary resistance 
!========================
  nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_r_b(i,j)*coef_dble,4)
       s_r_b(i,j)= 0._DBL
     end do
   end do
  where(tmp(:,:,1).eq.10000.) tmp(:,:,1) = tmp(:,:,1)-10000.
  name='S_Rb'
  long_name='Leaf boundary Resistance'
  units='s/m'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
!========================
! stomatal resistance 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_r_c(i,j)*coef_dble,4)
       s_r_c(i,j)= 0._DBL
     end do
   end do
  where(tmp(:,:,1).eq.10000.) tmp(:,:,1) = tmp(:,:,1)-10000.
  name='S_Rc'
  long_name='Stomatal Resistance'
  units='s/m'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
!========================
! undercanopy resistance 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_r_d(i,j)*coef_dble,4)
       s_r_d(i,j)= 0._DBL
     end do
   end do
  where(tmp(:,:,1).eq.10000.) tmp(:,:,1) = tmp(:,:,1)-10000.
  name='S_Rd'
  long_name='Undercanopy Resistance'
  units='s/m'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
!========================
! soil surface resistance 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_r_soil(i,j)*coef_dble,4)
       s_r_soil(i,j)= 0._DBL
     end do
   end do
  where(tmp(:,:,1).eq.10000.) tmp(:,:,1) = tmp(:,:,1)-10000.
  name='S_Rsoil'
  long_name='Soil Resistance'
  units='s/m'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
!========================
! canopy sensible heat flux 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_shf_canop(i,j)*coef_dble,4)
       s_shf_canop(i,j)= 0._DBL
     end do
   end do
  name='S_SHFC'
  long_name='canopy sensible heat flux'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
!========================
! soil sensible heat flux 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_shf_soil(i,j)*coef_dble,4)
       s_shf_soil(i,j)= 0._DBL
     end do
   end do
  name='S_SHFS'
  long_name='soil sensible heat flux'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
!========================
! canopy latent heat flux 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_lhf_canop(i,j)*coef_dble,4)
       s_lhf_canop(i,j)= 0._DBL
     end do
   end do
  name='S_LHFC'
  long_name='canopy latent heat flux'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
!========================
! soil latent heat flux 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_lhf_soil(i,j)*coef_dble,4)
       s_lhf_soil(i,j)= 0._DBL
     end do
   end do
  name='S_LHFS'
  long_name='soil latent heat flux'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
!========================
! wet latent heat flux 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_lhf_wet(i,j)*coef_dble,4)
       s_lhf_wet(i,j)= 0._DBL
     end do
   end do
  name='S_LHFCW'
  long_name='canopy wet latent heat flux'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
!========================
! dry latent heat flux 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_lhf_tr(i,j)*coef_dble,4)
       s_lhf_tr(i,j)= 0._DBL
     end do
   end do
  name='S_LHFCD'
  long_name='canopy dry latent heat flux'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
!========================
! ground heat flux 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_grflux0(i,j)*coef_dble,4)
       s_grflux0(i,j)= 0._DBL
     end do
   end do
  name='S_GRFLUX0'
  long_name='ground heat flux'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
!========================
! Precipitation infiltrate rate
!========================

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_precip_ref(i,j)*86400._DBL*coef_dble,4)
       s_precip_ref(i,j)= 0._DBL
     end do
   end do
  name='S_REF'
  long_name='Precipitation at reference level'
  units='mm/d'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
!========================
! Precipitation infiltrate rate 
!========================

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_precip_in(i,j)*86400._DBL*coef_dble,4)
       s_precip_in(i,j)= 0._DBL
     end do
   end do
  name='S_PIN'
  long_name='Precipitation infiltrate rate'
  units='mm/d'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
!========================
! surface run off rate 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_run_off_sfc(i,j)*86400._DBL*coef_dble,4)
     s_run_off_sfc(i,j)= 0._DBL
     end do
   end do
  name='S_ROFFS'
  long_name='surface runoff rate'
  units='mm/d'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
!========================
! leaf temperature 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_tcnp(i,j)*coef_dble,4)
       s_tcnp(i,j)= 0._DBL
     end do
   end do
  where(tmp(:,:,1).ne.0.) tmp(:,:,1) = tmp(:,:,1) - tfriz
  name='S_TCNP'
  long_name='canopy leaf temperature'
  units='C'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
!========================
! CAS temperature 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_tcas(i,j)*coef_dble,4)
       s_tcas(i,j)= 0._DBL
     end do
   end do
  where(tmp(:,:,1).ne.0.) tmp(:,:,1) = tmp(:,:,1) - tfriz
  name='S_TCAS'
  long_name='canopy air space temperature'
  units='C'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
!========================
! CAS q 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_qcas(i,j)*1000._DBL*coef_dble,4)
       s_qcas(i,j)= 0._DBL
     end do
   end do
  name='S_QCAS'
  long_name='canopy air space q'
  units='g/kg'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
!========================
! canopy water storage 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_mw(i,j)*coef_dble,4)
       s_mw(i,j)= 0._DBL
     end do
   end do
  name='S_MW'
  long_name='canopy water storage'
  units='kg/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
!========================
! drainage rate 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_drainage(i,j)*86400._DBL*coef_dble,4)
       s_drainage(i,j)= 0._DBL
     end do
   end do
	
  name='S_WDRAIN'
  long_name='water drainage rate from bottom soil layer'
  units='mm/d'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

!========================
! puddle water storage
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_mws(i,j)*coef_dble,4)
       s_mws(i,j)= 0._DBL
     end do
   end do
  name='S_MWS'
  long_name='puddle water storage'
  units='kg/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

!========================
! net upward SW to ref.level 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_net_swup(1,i,j),4)*coef
       s_net_swup(1,i,j)= 0.
     end do
   end do
  name='S_SUP'
  long_name='SW upward to ref.level'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

!========================
! net upward LW to ref. level
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_net_lwup(1,i,j),4)*coef
       s_net_lwup(1,i,j)= 0.
     end do
   end do
  name='S_LUP'
  long_name='LW upward to ref. level'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
!========================
! Clay Content
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=maxval(CLAY(:,i,j))
     end do
   end do
  name='S_CLAY'
  long_name='Clay Maximum  Content'
  units='%'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

!========================
! Sand Content
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=maxval(SAND(:,i,j))
     end do
   end do
  name='S_SAND'
  long_name='Sand  Maximum  Content'
  units='%'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

do k = 1,nsoil
write(numchar,'(i0)') k
!========================
! soil wetness 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_soilw(k,i,j)*coef_dble,4)
       s_soilw(k,i,j)= 0._DBL
     end do
   end do
  name='S_WS'//numchar
  long_name='soil wetness'//numchar
  units=' '
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
!========================
! soil temperature
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=real(s_soilt(k,i,j)*coef_dble,4)
       s_soilt(k,i,j)= 0._DBL
     end do
   end do
  where(tmp(:,:,1).ne.0.) tmp(:,:,1) = tmp(:,:,1) - tfriz
  name='S_TS'//numchar
  long_name='Soil temperature'//numchar
  units='C'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

end do ! k
end subroutine
