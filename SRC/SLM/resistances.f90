!=============================================================================
! subroutine resistances
! Note:
! 	compute aerodynamic resistances (r_b, r_d), stomatal resistance (r_c)
!=============================================================================

SUBROUTINE resistances(i,j)
use slm_vars, only : DBL, r_a, r_b, r_c, r_d, r_litter, vegetype, z0_soil, ustar, LAI, ggr, soilt, &
		Rgl, LAI, Rc_min, Rc_max, hs_rc, t_cas, pres0, q_cas, T_opt, &
		rootF, soilw, w_s_FC, s_depth, w_s_WP, theta_FC, theta_WP, poro_soil, &
		ztop, nsoil, dtfactor_dbl
use rad, only: swdsxy
IMPLICIT NONE

INTEGER, INTENT(IN) :: i,j
	
REAL (KIND=DBL) :: Cs, 	    & ! Total Turbulent transfer coef between canopy and underlying soil
	  	   Cs_bare ,& !  turbulent transfer coefficient over baresoil
		   Cs_dense   !  turbulent transfer coefficient under dense canopy : Sakaguchi and Zeng (2009)
REAL (KIND=DBL) :: z_m
REAL (KIND=DBL) :: rb_correc_fac, rd_correc_fac, temp_diff	
REAL (KIND=DBL) :: rc_fac_rad,&  ! stomatal resistance factor from radiation
		   rc_fac_vpd,&  ! stomatal resistance factor from vapor pressure deficit
		   rc_fac_t,&    ! stomatal rersistance factor from temperature
		   rc_fac_sw     ! stomatal resistance factor from soil moisture

REAL (KIND=DBL) :: tmp, tmp2, temp, d_root
INTEGER :: k
real, external :: qsatw


!==============================================================================
! Note:
! Initialize resistances with large number to prevent "divide by zero error" 
! when surface heat and water vapor fluxes are computed
!==============================================================================
r_b = 1.e4_DBL
r_c = 1.e4_DBL
r_d = 1.e4_DBL
r_litter = 1.e8_DBL

SELECT CASE (vegetype(i,j))
CASE(0) ! baresoil case
r_d = r_a ! under canopy resistance

CASE(1:) ! vegetated land surfaces
!=================================================================================
! Aerodynamic resistance for heat and vapor transfer under canopy space : r_d 
!=================================================================================
! temp_diff > 0 :: stable undercanopy
! temp_diff < 0 :: unstable undercanopy
temp_diff = t_cas(i,j) - soilt(1,i,j)
	
! rd_correc_fac : undercanopy stability parameter, in effect only for stable undercanopy
rd_correc_fac = ggr*ztop(i,j)*max(0.0_DBL,temp_diff)/soilt(1,i,j)/(ustar(i,j)**2)		

! turbulent transfer coefficient under dense canopy
IF(temp_diff.lt.0.0_DBL) THEN  
	Cs_dense = 0.004_DBL
ELSE                   
	! stable undercanopy : Cs_dense becomes smaller than 0.004
	Cs_dense = 0.004_DBL/(1.0_DBL+0.5_DBL*MIN(10.0_DBL,rd_correc_fac))
END IF

! turbulent transfer coefficient over the exposed topsoil
! typical value of Cs_bare ~0.2
Cs_bare = 0.4_DBL/0.13_DBL*(z0_soil*ustar(i,j)/(1.5e-5_DBL))**(-0.45_DBL) 

! Turbulence transfer coefficient undercanopy 
! LAI-weighed sum of Cs_bare & Cs_dense
Cs = Cs_bare*exp(-1.0_DBL*LAI(i,j))+Cs_dense*(1.0_DBL-exp(-1._DBL*LAI(i,j)))
		
! Undercanopy aerodynamic resistance depends on the weighed sum of the dense canopy covered soil 
! and baresoil turbulent transfer coefficient and friction velocity 
! Reference :[Oleson et al., 2004] [Zeng et al., 2005]
r_d = 1.0_DBL/ustar(i,j)/Cs

!==================================================================================
! Leaf boundary layer resistance : r_b
!==================================================================================
! turbulent transfer coefficient between canopy surface and canopy air : Cv = 0.01m/s^-0.5
! characteristic dimension of the elaves in the direction of wind flox : d_leaf = 0.04m
r_b = 1.0_DBL/0.01_DBL*((ustar(i,j)/0.04_DBL)**(-0.5_DBL))/max(0.1_DBL,LAI(i,j))

!=================================================================================
! Stomatal resistance : r_c
!=================================================================================
! radiation factor
tmp = 0.55_DBL*swdsxy(i,j)*2.0_DBL/Rgl(i,j)/LAI(i,j)
rc_fac_rad = (Rc_min(i,j)/Rc_max + tmp)/(1.0_DBL+tmp)

! vapor pressure deficit factor
rc_fac_vpd = 1.0_DBL/(1._DBL+hs_rc(i,j)*(qsatw(real(t_cas(i,j),KIND(1.)),pres0)-q_cas(i,j)))

! temperature factor
rc_fac_t = 1._DBL-0.0016_DBL*(T_opt-t_cas(i,j))**2

! rootzone soil moisture factor
rc_fac_sw = 0._DBL
d_root = 0.
DO k = 1,nsoil
	if(rootF(k,i,j).gt.0.0_DBL) then ! soil layer with the root
		if(soilw(k,i,j).gt.w_s_FC(k,i,j)) then
			! no water stree if water level exceed the value at field capacity
			temp = 1.*s_depth(k,i,j) 
                        d_root = d_root + s_depth(k,i,j)
		else if(soilw(k,i,j).lt.w_s_WP(k,i,j)) then
			! below wilting point
			temp = 0.0_DBL
				
		else   ! otherwise
			temp = s_depth(k,i,j) &
				*(soilw(k,i,j)*poro_soil(k,i,j)-theta_WP(k,i,j)) &
				/(theta_FC(k,i,j)-theta_WP(k,i,j))
                        d_root = d_root + s_depth(k,i,j)
		end if
		rc_fac_sw = rc_fac_sw + temp
	end if
END DO
rc_fac_sw = rc_fac_sw/max(1.e-6,d_root)
	
tmp2 = max(1.e-6_DBL,rc_fac_rad*rc_fac_vpd*rc_fac_t*rc_fac_sw)
r_c = min(Rc_max,(Rc_min(i,j)/LAI(i,j)/tmp2))

! collect statistics for 2D output
!s_rc_f1(i,j) = s_rc_f1(i,j) + rc_fac_rad*dtfactor_dbl
!s_rc_f11(i,j) = s_rc_f11(i,j) + swnsxy(i,j)*dtfactor_dbl
!s_rc_f2(i,j) = s_rc_f2(i,j) + rc_fac_vpd*dtfactor_dbl
!s_rc_f21(i,j) = s_rc_f21(i,j) + qsatw(real(t_cas(i,j),4),pres0)*dtfactor_dbl
!s_rc_f3(i,j) = s_rc_f3(i,j) + rc_fac_t*dtfactor_dbl
!s_rc_f4(i,j) = s_rc_f4(i,j) + rc_fac_sw*dtfactor_dbl

!=================================================================================
! r_litter : litter resistance : not in use in this version
!=================================================================================
r_litter = 0._DBL		
! Ref.[Sakaguchi and Zeng 2009]
! set litter LAI as 0.5
!r_litter = 1.0_DBL/0.004_DBL/ustar*(1._DBL-exp(-0.5_DBL))		
END SELECT
END SUBROUTINE resistances

