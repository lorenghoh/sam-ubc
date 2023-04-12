!=================================================================================================
! module soil_proc
! Note: module contains subroutines for soil temperature and moisture calculations
!=================================================================================================

MODULE soil_proc
use slm_vars, only : DBL,soilw, soilw_inc, soilt, precip_in, precip_sfc, ks, run_off_sfc,nsoil,&
		sh_eff_cond, s_depth_mm,Bconst,m_pot_sat, poro_soil, sh_eff_vel, evapo_dry, &
		dtn_dbl,sst_cond, sst_capa,st_cond, st_capa, s_depth, st_eff_cond, alpha, beta, &
		lhf_soil, rootF, lcond, grflux0, w_s_WP, dtfactor_dbl, landtype, tfriz, &
                drainage_flux, mws, mws_mx, mws_inc
IMPLICIT NONE
INTEGER :: k
REAL (KIND=DBL) :: aa, bb, cc, dd
PRIVATE :: k, aa, bb, cc, dd
PUBLIC :: soil_water, soil_temperature

CONTAINS

SUBROUTINE soil_water(i,j)

INTEGER, INTENT(IN) :: i,j
REAL(KIND=DBL) :: sw_wgt(nsoil), drain, puddle

if(landtype(i,j).eq.15) then

  soilw(:,i,j) = 0.
  soilw_inc(:) = 0.
  precip_in = 0.
  run_off_sfc = 0.
  sh_eff_cond(:) = 0.
  sh_eff_vel(:) = 0.

else

  soilw_inc(:) = soilw(:,i,j)

!=================================================================================================
! Calculate precipitation infiltration rate into the first soil layer
! Note:
! Infiltration occurs only the top layer is unsaturated and ground temperature is above freezing point
! Max. infiltration rate is equal to top soil hydraulic conductivity @ sat.(ks).
! Marat: modification Feb 2017:
! before running off, the excess of water accumulates in puddle until maximum amount is reached.
! puddle water continues to infiltrate even  when there is no rain.
!=================================================================================================
  puddle = mws(i,j)/dtn_dbl
  IF((soilw(1,i,j).lt.1.0).and.(soilt(1,i,j).ge.tfriz)) THEN
        precip_in = MIN(precip_sfc+puddle, ks(1,i,j))  ! unit : kg/m2/s
  ELSE
        precip_in = 0.0_DBL
  END IF
  puddle = puddle + precip_sfc - precip_in
  drain = max(0.,puddle - mws_mx(i,j)/dtn_dbl)
  puddle = puddle - drain
  mws_inc = puddle*dtn_dbl - mws(i,j)

  run_off_sfc = drain

!=================================================================================================
! calculate diffusion coefficient and velocity for soil moisture transfer 
! at each interfacial layer(= between adjacent soil layers)
! Note:
! depth-weighed averaging method is used
! sh_eff_cond : diffusion coefficient in diffusion equation ; mm*mm/s
! sh_eff_vel : effective velocity in advection term at interface layer : mm/s
! ks = hydraulic conductivity at saturation
! Bconst  = constant
! m_pot_sat = moisture potential at saturation
! soilw = soil wetness		  
! poro_soil = porosity		
! above parameter values are calculated based on percentage of SAND and CLAY as specified in input file
!=================================================================================================
  DO k = 1, nsoil-1
    if(soilt(k,i,j).ge.tfriz.and.soilt(k+1,i,j).ge.tfriz) then
      sh_eff_cond(k) &
	= (s_depth_mm(k,i,j)*(soilw(k,i,j)**(Bconst(k,i,j)+2._DBL))&
	+s_depth_mm(k+1,i,j)*(soilw(k+1,i,j)**(Bconst(k,i,j)+2._DBL)))&
	/(s_depth_mm(k,i,j)+s_depth_mm(k+1,i,j)) &
        *ks(k,i,j)*Bconst(k,i,j)*ABS(m_pot_sat(k,i,j))/poro_soil(k,i,j)

      sh_eff_vel(k) = &
	(s_depth_mm(k,i,j)*(soilw(k,i,j)**(2._DBL*Bconst(k,i,j)+2._DBL))&
	+s_depth_mm(k+1,i,j)*(soilw(k+1,i,j)**(2._DBL*Bconst(k,i,j)+2._DBL)))&
	/(s_depth_mm(k,i,j)+s_depth_mm(k+1,i,j)) &
        *ks(k,i,j)/poro_soil(k,i,j) 
    else ! no water movment between two layers one of which is frozen
      sh_eff_cond(k) = 0.
      sh_eff_vel(k) = 0.
    end if
  END DO 	

!=================================================================================================
! from FDE by the implicit method, Thomas algorithm is applied
! aa: terms related with soilw(k-1)
! bb: terms related with soilw(k)
! cc: terms related with soilw(k+1)
! dd: current soil wetness - sink + source
!=================================================================================================
  sw_wgt = 1.
  where(soilw(:,i,j).lt.w_s_WP(:,i,j)) sw_wgt = 0.


  aa = 0._DBL
  cc = -1.*2.*sh_eff_cond(1)*dtn_dbl/(s_depth_mm(1,i,j)*(s_depth_mm(1,i,j)+s_depth_mm(2,i,j)))  
  dd = soilw(1,i,j)&  			! soil wetness at time step n
   -(lhf_soil(i,j)/lcond & 		! mm/s ; soil top evaporation
   +rootF(1,i,j)*evapo_dry*sw_wgt(1) & 	! mm/s ; transpiration
   -precip_in) & ! mm/s 		! precipitation infiltration
   *dtn_dbl/poro_soil(1,i,j)/s_depth_mm(1,i,j) ! s/mm

  bb = 1._DBL-cc + sh_eff_vel(1)*dtn_dbl/(s_depth_mm(1,i,j))
  alpha(1) = cc/bb
  beta(1) = dd/bb

! for 2-(nsoil-1) layer
  DO k = 2, nsoil-1
	aa = -1.*dtn_dbl/s_depth_mm(k,i,j) & ! s/mm
		* (sh_eff_vel(k-1) & ! mm/s
		+ sh_eff_cond(k-1)*2._DBL/(s_depth_mm(k-1,i,j)+s_depth_mm(k,i,j))) ! mm/s
	
	cc = -1.*dtn_dbl/s_depth_mm(k,i,j) & ! s/mm
	        * 2._DBL*sh_eff_cond(k)/(s_depth_mm(k,i,j)+s_depth_mm(k+1,i,j)) ! mm/s
	bb = 1._DBL-cc &
		+ dtn_dbl/s_depth_mm(k,i,j)*sh_eff_vel(k) &
		+ 2._DBL*dtn_dbl/s_depth_mm(k,i,j)*sh_eff_cond(k-1)/(s_depth_mm(k-1,i,j)+s_depth_mm(k,i,j))
        
	dd = soilw(k,i,j) & ! current time step
	     - rootF(k,i,j)*evapo_dry*sw_wgt(k) & ! tranpiration mm/s
		*dtn_dbl/poro_soil(k,i,j)/s_depth_mm(k,i,j)  ! s/mm

	alpha(k) = cc/(bb-aa*alpha(k-1))
	beta(k) = (dd-aa*beta(k-1))/(bb-aa*alpha(k-1))
  END DO

  ! for bottom layer
  aa = -1.*dtn_dbl/s_depth_mm(nsoil,i,j) & ! s/mm
		* (sh_eff_vel(nsoil-1) & ! mm/s
		+ sh_eff_cond(nsoil-1)*2._DBL/(s_depth_mm(nsoil-1,i,j)+s_depth_mm(nsoil,i,j))) ! mm/s
  cc = 0._DBL	
  bb = 1._DBL &
	+ (dtn_dbl/s_depth_mm(nsoil,i,j) &  ! s/mm
	* 2._DBL*sh_eff_cond(nsoil-1)/(s_depth_mm(nsoil-1,i,j)+s_depth_mm(nsoil,i,j))) ! mm/s

! drainge when it exceed 1. mm/s
  drainage_flux = max(soilw(nsoil,i,j)-1._DBL, 0._DBL)*poro_soil(nsoil,i,j)*s_depth_mm(nsoil,i,j)/dtn_dbl  
  dd = soilw(nsoil,i,j) & !current time step
     - (rootF(nsoil,i,j)*evapo_dry*sw_wgt(nsoil) & ! tranpiration mm/s
     + drainage_flux) & ! drainge when it exceed 1. mm/s
     *dtn_dbl/poro_soil(nsoil,i,j)/s_depth_mm(nsoil,i,j)  ! s/mm
  alpha(nsoil) = 0._DBL
  beta(nsoil) = (dd-aa*beta(nsoil-1))/(bb-aa*alpha(nsoil-1))

  ! (n+1) time step soil wetness
  soilw(nsoil,i,j) = beta(nsoil)
  DO k = nsoil-1, 1, -1
	soilw(k,i,j) = max(0.,beta(k)-alpha(k)*soilw(k+1,i,j))
  END DO 

  if(any(soilw(:,i,j).gt.1._DBL)) then
  ! fix the levels where wetness exceeds 1 preserving total water:
    dd = 0.
  ! compute the excess:
    do k=1,nsoil
     if(soilw(k,i,j).gt.1._dbl) then
       dd = dd + (soilw(k,i,j)-1._DBL)*poro_soil(k,i,j)*s_depth_mm(k,i,j) 
       soilw(k,i,j) = 1._DBL 
     end if
    end do
  ! distribute  the excess among layer into deepest layers first:
    do k=nsoil,1,-1
     if(soilw(k,i,j).lt.1._DBL) then
       cc = min(dd,(1.-soilw(k,i,j))*poro_soil(k,i,j)*s_depth_mm(k,i,j))
       soilw(k,i,j) = soilw(k,i,j) + cc/(poro_soil(k,i,j)*s_depth_mm(k,i,j))
       dd = dd  - cc
       if(dd.lt.0._DBL) exit
     end if
    end do
  end if

  soilw_inc(:) = soilw(:,i,j)-soilw_inc(:)


end if ! 

END SUBROUTINE soil_water

!=============================================================

SUBROUTINE soil_temperature(i,j)

INTEGER, INTENT(IN) :: i,j
! List of Local variables
REAL (KIND=DBL) :: k_dry, k_sat, Ke
REAL (KIND=DBL) :: temp

if(landtype(i,j).eq.15) then

  DO k = 1,nsoil
     st_cond(k) = 1.6_DBL ! thermal conductivity of ice
     st_capa(k) = 917.0_DBL*2030.0_DBL ! ice heat capacity
  END DO

else

  DO k = 1,nsoil
  ! dry soil density [kg/m3] ! 2700kg/m3 = soilds unit weight 
  ! Empirical method from FIFE site observation
  temp = (1.0_DBL-poro_soil(k,i,j))*2700.0_DBL 
 
  ! dry thermal conductivity [W/mK]
  k_dry = (0.135_DBL*temp+64.7_DBL)/(2700._DBL-0.947_DBL*temp) 

  ! Saturated soil thermal conductivity : unsing percent of sand and clay content
  if(soilt(k,i,j).gt.tfriz) then
      k_sat = ((sst_cond(k,i,j))**(1.0_DBL-poro_soil(k,i,j))) &
             *(0.57_DBL**poro_soil(k,i,j)) ! 0.57W/mK = thermal conductivity of water
  else
      k_sat = ((sst_cond(k,i,j))**(1.0_DBL-poro_soil(k,i,j))) &
             *(1.6_DBL**poro_soil(k,i,j)) ! 1.6W/mK = thermal conductivity of ice
  end if

  ! Weighing factor between dry and saturated soil thermal conductivity : Kersten number
  Ke = log10(max(0.1_DBL,soilw(k,i,j)))+1.0_DBL

  ! Total soil thermal conductivity at each node_z
  ! Johansen (1975) 
  st_cond(k) = Ke*(k_sat - k_dry) + k_dry

  !soil volumetric heat capacity at each node_z depth
  if(soilt(k,i,j).gt.tfriz) then
     st_capa(k) = (1.-poro_soil(k,i,j))*sst_capa(k,i,j) &  
	   + (998.0_DBL*4182.0_DBL)*soilw(k,i,j)*poro_soil(k,i,j) ! water heat capacity
  else
     st_capa(k) = (1.-poro_soil(k,i,j))*sst_capa(k,i,j) &  
	   + (917.0_DBL*2030.0_DBL)*soilw(k,i,j)*poro_soil(k,i,j) ! ice heat capacity
  end if
								  ! 998kg/m3 = water density
								  ! 917kg/m3 = ice density
								  ! 4182 J/kg/K = water specific heat
								  ! 2030 J/kg/K = ice specific heat
  END DO

end if

DO k = 1,nsoil-1
! Calculate effective conductiviy at the adjacent soil layer interface [at interface_z]
! depth weighed average
   st_eff_cond(k) = (st_cond(k+1)*s_depth(k+1,i,j) + st_cond(k)*s_depth(k,i,j))/(s_depth(k+1,i,j)+s_depth(k,i,j))
END DO! k

! from FDE by the implicit method, Thomas algorithm is applied
! aa: terms related with T(k-1)
! bb: terms related with T(k)
! cc: terms related with T(k+1)
! dd: current soil temperature + (additional source/sink on soil top)
! for first layer
aa = 0._DBL
cc = -1.*2.*st_eff_cond(1)*dtn_dbl/st_capa(1)/(s_depth(1,i,j)*(s_depth(1,i,j)+s_depth(2,i,j)))
bb = 1.-cc 
dd = soilt(1,i,j)&  ! soil temperature at time step n
	-grflux0*dtn_dbl/s_depth(1,i,j)/st_capa(1) ! soil top boundary layer contidion
alpha(1) = cc/bb
beta(1) = dd/bb


! for 2-(nsoil-1) layer
DO k = 2, nsoil-1
	aa = -2._DBL*st_eff_cond(k-1)*dtn_dbl/s_depth(k,i,j)/st_capa(k)/(s_depth(k-1,i,j)+s_depth(k,i,j))
	cc = -2._DBL*st_eff_cond(k)*dtn_dbl/s_depth(k,i,j)/st_capa(k)/(s_depth(k,i,j)+s_depth(k+1,i,j))
	bb = 1._DBL-cc-aa 
	dd = soilt(k,i,j) ! current time step

	alpha(k) = cc/(bb-aa*alpha(k-1))
	beta(k) = (dd-aa*beta(k-1))/(bb-aa*alpha(k-1))
END DO
! for bottom layer
aa = -2._DBL*st_eff_cond(nsoil-1)*dtn_dbl/s_depth(nsoil,i,j)/st_capa(nsoil)/(s_depth(nsoil-1,i,j)+s_depth(nsoil,i,j))
cc = 0._DBL	
bb = 1._DBL-aa 
dd = soilt(nsoil,i,j)

alpha(nsoil) = 0._DBL
beta(nsoil) = (dd-aa*beta(nsoil-1))/(bb-aa*alpha(nsoil-1))

! (n+1) time step soil temperature

soilt(nsoil,i,j) = beta(nsoil)
DO k = nsoil-1, 1, -1
	soilt(k,i,j) = beta(k)-alpha(k)*soilt(k+1,i,j)
END DO 

end subroutine soil_temperature
 
END MODULE soil_proc
