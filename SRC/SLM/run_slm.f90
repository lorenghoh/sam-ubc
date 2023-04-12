!=============================================================================================
! subroutine run_slm 
! Note:
!	main() for SLM
!	called from SAM/SRC/surface()
!
! History: March, 2016
!=============================================================================================

SUBROUTINE run_slm(ur, vr, precip_ref, qr, tr, ts, zr, pref, &
                   swdsvisxy, swdsnirxy, swdsvisdxy, swdsnirdxy, lwdsxy, coszrsxy, &
                   prsfc, flbu, flbv, flbq, flbt)
			
use slm_vars
use soil_proc, only : soil_water, soil_temperature

IMPLICIT NONE
!=============================================================================================
! Input : variables at reference level
!=============================================================================================
REAL , INTENT(IN) , DIMENSION(nx,ny) :: ur, & ! u-wind
				        vr, & ! v-wind
				        qr, & ! qv at reference level
				        tr, & ! t  at reference level
				        precip_ref, & ! precipitation rate from reference level
                                        swdsvisxy, & ! downward shortwave direct visible flux
                                        swdsnirxy, & ! downward shortwave direct near-IR flux
                                        swdsvisdxy, & ! downward shortwave diffuse visible flux
                                        swdsnirdxy, & ! downward shortwave idiffuse near-IR flux
                                        lwdsxy, & ! downward longwave flux
                                        coszrsxy, &  ! zenith angle
                                        zr, &  ! height of ref. level
					pref   ! pressure at ref. level
!==============================================================================================
! Output : 
!==============================================================================================
REAL , INTENT(OUT), DIMENSION(nx,ny) :: flbu, & ! surface u-momentum
				        flbv, & ! surface v-momentum
				        flbq, & ! surface vapor flux
				        flbt, & ! surface heat flux
					prsfc   ! precipitation rate reaching to soil surface

!====================================================================================
! INPUT & OUTPUT : surface temperature 
!====================================================================================
REAL , INTENT(INOUT), DIMENSION(nx,ny) :: ts  ! surface temperature

!====================================================================================
! Local Variables
!====================================================================================
REAL (KIND=DBL) :: t_sfc, q_sfc,tmp, cp_vege_tot, coef
INTEGER :: i,j,k
real, external :: qsatw,qsati

!====================================================================================
! note:
! Downward longwave (lwdsxy_slm), shortwave radiation (swdsxy_slm), 
! and cosine of zenith angle (coszrs_slm) are read from the equivalent SAM variables 
! at each time step.
! When the run restarts, restarts swdsxy_slm, lwdsxy_slm, coszrs_slm are overwritten 
! from the values stored in the LAND restart file. 
! Otherwise values will be zero as radiation is calculated after surface() in SAM.
!====================================================================================
if(flag_vege_init) then
	if(nrestart.eq.0) then  
		! initialize vegetation and soil parameters from input file
		call init_slm_vars(tr, ts, qr)  
		
		if(masterproc) write(*,*) 'initializing SLM ......'
	else ! for restart run
		
		! read LAND restart file
		call read_statement_slm()
		
		! read from restart file.
		nstep_restart = nstep
		if(masterproc) write(*,*) 'restart SLM time step at:', nstep_restart
	end if ! nrestart
        flag_vege_init = .false.
        call slm_printout
end if

!====================================================================================
! local variable initialization
cnp_mw_drip 	= 0._DBL ! dripping from canopy water storage when the storage exceeds mw_mx : mm/s
drain   	= 0._DBL ! drainige rate from canopy
dtn_dbl = dtn

DO j = 1,ny
DO i = 1,nx
 if(landmask(i,j).eq.1) then
!====================================================================================
! Note:
! Nudging soil wetness and soil temperature profile to initial over time period, tau_soil.
! soil_relax_hgt (in the range of 0 - 1) variable determines the effectiveness of nudging as set in the 
! initial soil sounding. 
!====================================================================================
  if(dosoilwnudging) then
        DO k = 1, nsoil
        soilw(k,i,j) = soilw(k,i,j)-(soilw(k,i,j)-soilw_obs(k,i,j))*dtn_dbl/tau_soil(i,j)*soil_relax_hgt(k,i,j)
        soilw_nudge(k,i,j) = (soilw(k,i,j)-soilw_obs(k,i,j))*soil_relax_hgt(k,i,j)/tau_soil(i,j)
        END DO ! k
  end if

  if(dosoiltnudging) then
        DO k = 1, nsoil
        soilt(k,i,j) = soilt(k,i,j)-(soilt(k,i,j)-soilt_obs(k,i,j))*dtn_dbl/tau_soil(i,j)*soil_relax_hgt(k,i,j)
        soilt_nudge(k,i,j) = (soilt(k,i,j)-soilt_obs(k,i,j))*soil_relax_hgt(k,i,j)/tau_soil(i,j)
        END DO ! k
  end if

!====================================================================================
! Note:
! precip 	: precipitation interception rate at canoppy
! precip_ref 	: precipiptation rate at the reference height, zr
! precip_extinc : extinction coefficient for rainfall penetration through 
!		  vegetation layer
! LAI 		: leaf area index
! For baresoil, precip = 0, as LAI = 0
!====================================================================================
   precip = DBLE(precip_ref(i,j))*(1.-exp(-1.*precip_extinc(i,j)*LAI(i,j)))
!====================================================================================
! Note:
! drain 	: drainage rate from canopy 
! cnp_mw_drip 	: dripping rate from canopy
! precip_sfc 	: precipitation reaching to soil surface
! mw 		: water storage of canopy [kg/m^2]
! mw_mx 	: maximum water storage of canopy
! mw_inc 	: (prognostic variable) increment of moisture storage on canopy
!====================================================================================
   if(mw(i,j).lt.mw_mx(i,j)) then  
 	drain = 0.0_DBL 
   else
        ! When water holding storage exceeds its maximum, 
	! no precipiataion is intercepted
        drain = precip

        ! Excess water storage gets drained from mw [kg/m2] : Caution unit change
        cnp_mw_drip = MAX(mw(i,j)-mw_mx(i,j),0._DBL)/dtn_dbl ! mm/s

        drain = drain + cnp_mw_drip 
   end if
                
   precip_sfc = DBLE(precip_ref(i,j))-precip+drain 
   mw_inc = dtn_dbl*(precip-drain)
!====================================================================================
! Note: 
! Calculate net radiation absorbed by canopy and soil surface  
!====================================================================================
   call radiative_fluxes(i,j)

!====================================================================================
! Note:
! Assign appropriate "surface level" values for each land type, for the calculation  
! of surface turbulent fluxes
! 
! For baresoil,   surface level = soil surface
! For vegetation, surface level = canopy level
!====================================================================================
! specific humidity at top soil
  if(soilt(1,i,j).gt.tfriz) then
   if(mws(i,j).gt.0.) then
     q_gr = qsatw(real(soilt(1,i,j),KIND(1.)),pres0)
   else
     q_gr = fh_calc(soilt(1,i,j),m_pot_sat(1,i,j), soilw(1,i,j), Bconst(1,i,j)) &
                    *qsatw(real(soilt(1,i,j),KIND(1.)),pres0)
   end if
  else
     q_gr = qsati(real(soilt(1,i,j),KIND(1.)),pres0)
  end if

  IF(vegetype(i,j).ne.0) THEN ! for vegetated surfaces
	t_sfc = t_cas(i,j)
	q_sfc = q_cas(i,j)
  ELSE			    ! for baresoil
	t_sfc = soilt(1,i,j)
	q_sfc = q_gr
  END IF
		
!====================================================================================
! Note:
! Calculate turbulent transfer coefficient between reference level and surface
! Input : pref  	reference level pressure
!	  tr    	reference level temperature
!	  qr    	reference level specific humidity
!	  ur,vr 	reference level zonal, meridional velocities
!	  zr		reference level height
!	  t_sfc 	surface level temperature
!	  q_sfc 	surface level specific humidity
!	  z0_sfc 	surface roughness length
!	  disp_hgt 	displacement height
!====================================================================================
  call transfer_coef(pref(i,j), t_sfc, tr(i,j), qr(i,j), q_sfc, ur(i,j), vr(i,j), zr(i,j), &
	z0_sfc(i,j),disp_hgt(i,j),i,j)

!====================================================================================
! Note:
! Calculate surface momentum fluxes
! Surface Stress : kg/m/s2
!====================================================================================
  taux_sfc = -1._DBL*mom_trans_coef*vel_m*DBLE(ur(i,j))*(pres0*100._DBL/287._DBL/DBLE(tr(i,j)))
  tauy_sfc = -1._DBL*mom_trans_coef*vel_m*DBLE(vr(i,j))*(pres0*100._DBL/287._DBL/DBLE(tr(i,j)))

!====================================================================================
! Note:
! Calculate aerodynamic resistances + stomatal resistance
!====================================================================================
  call resistances(i,j)

!====================================================================================
! Note:
! Sensible heat fluxes
!====================================================================================
  shf_canop(i,j) =(t_canop(i,j)-t_sfc)*rhow(1)*cp/r_b*vege_YES(i,j)
  shf_soil(i,j) = (soilt(1,i,j)-t_sfc)*rhow(1)*cp/r_d
  shf_air(i,j) = (t_sfc-DBLE(tr(i,j)))*rhow(1)*cp/r_a

  IF(vegetype(i,j).eq.0) THEN ! for baresoil
 	shf_soil(i,j) = shf_air(i,j)
  ELSE
        shf_air(i,j) = shf_canop(i,j)+shf_soil(i,j)
  END IF

!====================================================================================
! Note:
! Calculate temperature scale for z0hsfc 
!====================================================================================
 tstar(i,j) = -1._DBL*shf_air(i,j)/rhow(1)/cp/ustar(i,j)

!====================================================================================
! Note: 
! Calculate latent heat fluxes
! Input: q_sfc		surface level specific humidity
!	 qr		reference level specific humidity
!	 i,j		current index of the grids
!====================================================================================
  call vapor_fluxes(q_sfc,qr(i,j),i,j)
!====================================================================================
! Note:
! Calculate soil moisture increment
!====================================================================================
  call soil_water(i,j)

!====================================================================================
! Note:
! Calculate soil temperature increment
! grflux0:	Total heat flux on top of soil surface : top boundary condition
! grflux0 < 0 => heat is added into soil top
!====================================================================================
  grflux0 = -1.0_DBL*(net_rad(2) - shf_soil(i,j) - lhf_soil(i,j))  
  call soil_temperature(i,j)

!====================================================================================
! Note: 
! Update vegetation moisture storage
!====================================================================================
  mw(i,j) = mw(i,j) +  mw_inc

!====================================================================================
! Note:
! Update puddle water storage
!====================================================================================
  mws(i,j) = mws(i,j) +  mws_inc

!====================================================================================
! Note:
! Update vegetatation temperature
! cp_vege_tot [J/m2/K]	: total vegetation heat capacity 
!         = heat capacity from vegetation + heat capacity of water held on vegetation    
!====================================================================================
 ! canapy heat capasity: assume 0.001 m leaf thickness, basal area in sq.feet/acre=43560 m2/m2,
 ! 900 kg/m3 density of leaves and wood, 2800 J/kg/K specific heat capacity.
  cp_vege(i,j) = (LAI(i,j)*0.001+ztop(i,j)*BAI(i,j)/43560.)*900.*2800.
  cp_vege_tot = cp_vege(i,j) + mw(i,j)*1.e-3_DBL*cp_water

 ! max(1.e-3, vege_cp) is for the baresoil case where vege_cp = 0.(to prevent divide by zero)
 ! vege_YES keeps t_canop_inc zero for baresoil case
  t_canop_inc = dtn_dbl/max(1.e-3_DBL,cp_vege_tot)&
                    *(net_rad(1)-shf_canop(i,j)-lhf_canop(i,j))*vege_YES(i,j)
  t_canop(i,j) = t_canop(i,j) +  t_canop_inc
!====================================================================================
! Note:
! Compute diagnostic variables at canopy air space level
!====================================================================================

 ! Calculate heat conducatnces ; non-zero only for canopy land type
  cond_heat = 1._DBL/r_a + 1._DBL/r_b + 1._DBL/r_d
  cond_href = 1._DBL/r_a/cond_heat*vege_YES(i,j)
  cond_hcnp = 1._DBL/r_b/cond_heat*vege_YES(i,j)
  cond_hundercnp = 1._DBL/r_d/cond_heat*vege_YES(i,j)

 ! Calculate Vapor conductances
  cond_vapor = 1./r_a &
		+ wet_canop/(2.*r_b) &
		+ (1.-wet_canop)/(2.*r_b+r_c) &
		+ 1./(r_d+r_soil+r_litter)
  cond_vref  = 1./r_a/cond_vapor*vege_YES(i,j)
  cond_vcnp  = (wet_canop/(2.*r_b)/cond_vapor &
		  + (1.-wet_canop)/(2.*r_b+r_c)/cond_vapor)*vege_YES(i,j)
  cond_vundercnp = 1./(r_d+r_soil+r_litter)/cond_vapor*vege_YES(i,j)

! Canopy Air Spcae 
  IF(vegetype(i,j).eq.0) THEN ! baresoil
        cond_hcnp = 1._DBL
        cond_vundercnp = 1._DBL
  END IF
  t_cas(i,j) = DBLE(tr(i,j))*cond_href &
              + t_canop(i,j) *cond_hcnp &
              + soilt(1,i,j) *cond_hundercnp
  if(t_canop(i,j).ge.tfriz) then
    qsat_canop 	= qsatw(real(t_canop(i,j),KIND(1.)),pres0)
  else
    qsat_canop 	= qsati(real(t_canop(i,j),KIND(1.)),pres0)
  end if
  if(soilt(1,i,j).ge.tfriz) then
    q_gr 	= qsatw(real(soilt(1,i,j),KIND(1.)),pres0)&
		 *fh_calc(soilt(1,i,j),m_pot_sat(1,i,j),soilw(1,i,j),Bconst(1,i,j))
  else
    q_gr 	= qsati(real(soilt(1,i,j),KIND(1.)),pres0)
  end if
  q_cas(i,j) 	= DBLE(qr(i,j))*cond_vref &
             	+ qsat_canop*cond_vcnp &
             	+ q_gr*cond_vundercnp
!====================================================================================
! output variables
!====================================================================================
!  IF(vegetype(i,j).ne.0) then
!        t_skin(i,j) = t_canop(i,j)
!  else
!        t_skin(i,j) = soilt(1,i,j)
!  end if
! t_skin is computed from upwelling LW radiation in raditive_fluxes()

  ts(i,j)   = real(t_skin(i,j),4)
  flbu(i,j) = real(taux_sfc,4)!/rhow(1)  ! Surface Stress in x-direction 
  flbv(i,j) = real(tauy_sfc,4)!/rhow(1)  ! Surface Stress in y-direction
  flbq(i,j) = real(lhf_air(i,j),4)/(lcond*rhow(1)) ! latent heat flux [W/m2] -> [1.e-3m/s]
  flbt(i,j) = real(shf_air(i,j),4)/(cp*rhow(1))    ! sensible heat flux [W/m2] -> [Km/s]
  prsfc(i,j) = real(precip_sfc,4)

!====================================================================================
! Collect statistics
!====================================================================================

 ! collect 2D stat variables
  call collect_2D_stat_vars(i,j)

  s_precip_ref(i,j) = s_precip_ref(i,j) + precip_ref(i,j)*dtfactor_dbl

 end if ! landmask
END DO
END DO
 
cast = real(t_cas,4) 
canopt = real(t_canop,4) 
casq = real(q_cas,4) 

! Save land variables to restart file
if(mod(nstep,nstat*(1+nrestart_skip)).eq.0.or.nstep.eq.nstop.or.nelapse.eq.0) then

  call write_statement_slm() ! save restart file

end if

END SUBROUTINE run_slm


