subroutine albedo(lchnk, ncol, landmask, landtype, coszrs, albedovis_v, albedonir_v,  &
                  albedovis_s, albedonir_s, soilw, phi_1, phi_2, LAI, asdir, aldir, asdif, aldif )
  !-----------------------------------------------------------------------
  ! Computes surface albedos over ocean for Slab Ocean Model (SOM)

  ! and the surface (added by Marat Khairoutdinov)

  ! Two spectral surface albedos for direct (dir) and diffuse (dif)
  ! incident radiation are calculated. The spectral intervals are:
  !   s (shortwave)  = 0.2-0.7 micro-meters
  !   l (longwave)   = 0.7-5.0 micro-meters
  !
  ! Uses knowledge of surface type to specify albedo, as follows:
  !
  ! Ocean           Uses solar zenith angle to compute albedo for direct
  !                 radiation; diffuse radiation values constant; albedo
  !                 independent of spectral interval and other physical
  !                 factors such as ocean surface wind speed.
  !
  ! For more details , see Briegleb, Bruce P., 1992: Delta-Eddington
  ! Approximation for Solar Radiation in the NCAR Community Climate Model,
  ! Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
  !
  !---------------------------Code history--------------------------------
  !
  ! Original version:         CCM1
  ! Standardized:             J. Rosinski, June 1992
  ! Reviewed:                 J. Kiehl, B. Briegleb, August 1992
  ! Rewritten for ocean only: J. Rosinski, May 1994
  ! Modified for SOM:         B. Briegleb  April 1995
  ! Reviewed:                 B. Briegleb  March 1996
  ! Reviewed:                 J. Kiehl     April 1996
  ! Reviewed:                 B. Briegleb  May   1996
  ! Added Land albedo	    M. Kairoutdinov Sep. 1999
  !
  !-----------------------------------------------------------------------
  !
  ! $Id: albedo.f90,v 1.1 2004/09/08 19:46:02 samcvs Exp $
  ! $Author: samcvs $
  !
  use shr_kind_mod, only: r4 => shr_kind_r4
  use ppgrid
  use params, only: SLM
  implicit none
  !------------------------------Arguments--------------------------------
  !
  ! Input arguments

  !
  integer lchnk,ncol

  real(r4) coszrs    ! Cosine solar zenith angle
  integer landmask 
  integer landtype
  real albedovis_v, albedovis_s  ! visible albedo of vegetation and soil
  real albedonir_v, albedonir_s  ! near IR albedo of land
  real soilw  ! soil wetness
  real phi_1, phi_2, LAI
  !
  ! Output arguments
  !
  real(r4) asdir     ! Srf alb for direct rad   0.2-0.7 micro-ms
  real(r4) aldir     ! Srf alb for direct rad   0.7-5.0 micro-ms
  real(r4) asdif     ! Srf alb for diffuse rad  0.2-0.7 micro-ms
  real(r4) aldif     ! Srf alb for diffuse rad  0.7-5.0 micro-ms
  !
  !---------------------------Local variables-----------------------------
  !
  real(r4), parameter :: adif = 0.06
  real ka, explai, wetfactor
  !
  !
  ! Initialize all ocean surface albedos to zero
  !
  asdir = 0.
  aldir = 0.
  asdif = 0.
  aldif = 0.

  if (coszrs.lt.0.) return

  if (landmask.eq.0) then
     ! Ice-free ocean albedos function of solar zenith angle only, and
     ! independent of spectral interval (Briegleb et al 1986):
        aldir = .026/(coszrs**1.7 + .065) + .15*(coszrs - 0.10)*(coszrs - 0.50)*(coszrs - 1.00)
        asdir = aldir
        aldif = adif
        asdif = adif

  else ! land

    if(SLM) then
      if(landtype.ne.15) then
        ka = phi_1/coszrs + phi_2
        wetfactor = 1.-0.5*soilw
        explai = exp(-(phi_1/coszrs+phi_2)*LAI)
        asdir = albedovis_v*(1.-explai)+albedovis_s*wetfactor*explai
        aldir = albedonir_v*(1.-explai)+albedonir_s*wetfactor*explai
        ! assume diffuse albedo is just albedo for coszrs=0 (Marat K): 
        explai = exp(-(phi_1+phi_2)*LAI)
        asdif = albedovis_v*(1.-explai)+albedovis_s*wetfactor*explai
        aldif = albedonir_v*(1.-explai)+albedonir_s*wetfactor*explai
      else
     ! albedo of ice
        aldir = 0.65
        asdir = 0.95
        aldif = 0.65
        asdif = 0.95
      end if
     else
     ! Albedos for land type I (Briegleb)
        asdir = 1.4 * 0.06 / ( 1. + 0.8 * coszrs)
        asdif = 1.2 * 0.06
        aldir = 1.4 * 0.24 / ( 1. + 0.8 * coszrs)
        aldif = 1.2 * 0.24
     end if 

  endif
      
  return
end subroutine albedo

