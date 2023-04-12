!=====================================================================================================================
! subroutine vapor_fluxes
! Note: 
! computes latent heat fluxes from topsoil, canopy (for vegetated land), and combines them to return the total to SAM
!
! Input:
! qsfc: specific humidity at surface level 
!	(the surface that the atmosphere sees: for instance, topsoil for baresoil case, and canopy for vegetated land)
! qr:   specific humidity at the reference level (input from the lowest level of SAM) 
! i,j : indices for the grid location
!	(Indices are passed to access other variables of the same grid point, in case landmask includes ocean)
!
! History : Jungmin Lee, March, 2016
!=====================================================================================================================

 SUBROUTINE vapor_fluxes(qsfc,qr,i,j)
	
 use slm_vars, only : DBL,nsoil, vegetype, soilw, w_s_FC, q_gr, pii, soilt, r_a, r_d, r_litter, r_b, &
		      w_s_WP,evapo_s, vege_YES, rhow, wet_canop, qsat_canop, t_canop, pres0, mw, &
		     mw_mx, evapo_wet,rootF, evapo_dry, dtn_dbl, lhf_soil, lhf_canop, lhf_air, lcond, r_c, r_soil, mw_inc
 IMPLICIT NONE

 REAL (KIND=DBL), INTENT(IN) :: qsfc 
 REAL, INTENT(IN) :: qr 	
 INTEGER,INTENT(IN) :: i,j
        
 INTEGER :: k
 REAL (KIND=DBL) :: totalR_soil,& ! total aerodynamic resistance over soil surface
		   qref_tmp, & 
		   soil_diff      ! soil diffusion factor for the topsoil evaporation from soil pore space to the overlying air
real, external :: qsatw
	
!=====================================================================================================================
! Evaporation from the Top soil layer 
!=====================================================================================================================

!=====================================================================================================================
! New reference level, which is the end point for water vapor turbulent flux, is defined over the vegetated land 
! 	for baresoil soil => reference level q stays to be qr
! 	for vegetated => vegetation layer q, qsfc, become the new reference level q
!=====================================================================================================================
 IF(vegetype(i,j).eq.0) THEN ! for baresoil
	qref_tmp = DBLE(qr)
 ELSE
	qref_tmp = qsfc
 END IF

!=====================================================================================================================
! Soil diffusion factor (soil_diff) follows "Sakaguchi and Zeng, 2009", but the snow factor is discarded in SLM adoptation
! Maximum evaporation is assumed when dew formation and/or soil top layer wetness is greater than field capacity
!=====================================================================================================================
 IF((soilw(1,i,j).ge.w_s_FC(1,i,j)).or.(qref_tmp.gt.q_gr)) THEN
	soil_diff = 1._DBL
 ELSE
       	soil_diff = 0.25_DBL*((1.0_DBL-cos(pii*max(0.01,soilw(1,i,j))/w_s_FC(1,i,j)))**2)
 END IF	
	
!=====================================================================================================================
! total resistance for the topsoil evaporation = moisture diffusion factor + aerodynamic resistance
! r_soil base value is applied to prevent too active evaporation and roughly represents the ground litter resistance
! r_litter at this version is not in use, contains zero
!=====================================================================================================================
 IF(vegetype(i,j).eq.0) then ! for baresoil
	r_soil = min(10000._DBL,max(100._DBL,r_a*(1.0_DBL/soil_diff-1.0_DBL))) 
	totalR_soil = r_soil + r_a + r_litter
 ELSE
	r_soil = min(10000._DBL,max(50._DBL,r_d*(1.0_DBL/soil_diff-1.0_DBL)))
	totalR_soil = r_soil + r_d + r_litter
 END IF
	
!=====================================================================================================================
! Evaporation from soil top underneath the canopy
!=====================================================================================================================
 evapo_s  = vege_YES(i,j)*(q_gr-qsfc)*rhow(1)/totalR_soil
	
!=====================================================================================================================
! Evaporation from canopy
!=====================================================================================================================
 wet_canop = 1._DBL 
 qsat_canop = qsatw(real(t_canop(i,j),KIND(1.)),pres0)
 IF(vegetype(i,j).ne.0.and.qsat_canop.gt.qsfc) THEN
	! wet portion of canopy
	wet_canop = min(1.0_DBL,mw(i,j)/mw_mx(i,j))
 END IF

 IF(vegetype(i,j).eq.0) wet_canop = 0._DBL

 ! direct evaporation from the water held on canopy
 evapo_wet = (qsat_canop-qsfc)*rhow(1)*wet_canop/(2.0_DBL*r_b)*vege_YES(i,j)

 ! increment/decrement of the water amount held on leaves following the direct evaporation/dew formation
 mw_inc = mw_inc-dtn_dbl*evapo_wet  ! evapo_wet [kg/m2s=mm/s] 
	
 ! Trainspiration 
 ! For dew formation situation (qsat_canop < qsfc), wet_canop = 1. automatically makes evapo_dry = 0.
 evapo_dry = (qsat_canop-qsfc)*rhow(1)*(1.0_DBL-wet_canop)/(2.0_DBL*r_b+r_c)*vege_YES(i,j)
       
 ! Check soil moisture availability for transpiration
 do k = 1, nsoil
       if(soilw(k,i,j).lt.0.1) then
                evapo_dry = evapo_dry -evapo_dry*rootF(k,i,j)
       end if
 end do! k

!=====================================================================================================================
! Convert evaporation (kg/m2/s) to latent heat flux (W/m2)
!=====================================================================================================================
 lhf_canop(i,j) = lcond*(evapo_wet+evapo_dry)
 lhf_soil(i,j) = lcond*evapo_s	
 lhf_air(i,j) = lcond*(qsfc-DBLE(qr))*rhow(1)/r_a
		
 IF(vegetype(i,j).eq.0) THEN
 ! 	For baresoil case, lhf_air is equal to lhf_soil, lhf_air is updated to include soil diffusion factor 
	lhf_air(i,j) = lhf_air(i,j) * r_a/(r_soil+r_a)
	evapo_s = lhf_air(i,j)/lcond
	lhf_soil(i,j) = lhf_air(i,j)
 else 
	lhf_air(i,j) = lhf_canop(i,j) + lhf_soil(i,j)
 end if

!energy_bal = energy_bal - lhf_air
!energy_bal_c = energy_bal_c - lhf_canop 
!energy_bal_s = energy_bal_s - lhf_soil

END SUBROUTINE vapor_fluxes		
