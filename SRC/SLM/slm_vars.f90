!==============================================================================
! module slm_vars
! Note:
! contains the declaration of variables and parameters
! contains function:
!	function fh_calc	computes 
! contains subroutines:
!	subroutine init_soil_tw		1.read 'soil' under case directory
!					2.allocate dimension for the soil related 
!					  variables following the number of soil specified in 'soil'
!	subroutine vege_root_init	1. computes soil depths of each layer center and layer's interface
!					2. computes the density fraction of total root present in each soil layer
!	subroutine slm_init		1. read vegetation type map under case directory
!					2. Assign vegetion type specific parameters for each grid following the predefined vegetation parameter file under RUNDATA
!					3. Assign soil parameter values following the soil type read from 'soil' 
!	subroutine init_slm_vars
!	subroutine collect_2D_stat_vars 1. collect 2D statistics
!	subroutine slm_stepout		1. printout temperature and surface heat fluxes computed by SLM, called in SRC/stepout()
!	subroutine slm_setparm		1. read parameters specified in prm
!	subroutine slm_stat_2Dinit    1. initialize 2D statistics
!	subroutine fminmax_print_slm	1. compute min/max across land grids
!	subroutine slm_printout	1. printout vegetation and soil parameters  	
!
! History: March, 2016 ! implmenetation in SAM6.10.9  by Jungmin Lee
!          December 2016 ! Implmenetation by Marat, extending to 16 land types, etc.
!==============================================================================

MODULE slm_vars
use grid, only : save2dbin,case, caseid,rank,restart_sep,nsubdomains,dompi,&
		 nx,ny,nx_gl, ny_gl, dtn,dt,day,pres0,dtfactor, masterproc, dostatis,nstep,nprint,&
		 nsave2D,nstat,nrestart, nstop, nelapse,nrestart_skip,&
                 case_restart, caseid_restart
use params, only : lcond, ggr, rv, cp, rgas, epsv, SLM
use vars, only : rhow, landmask 

IMPLICIT NONE

INTEGER, parameter :: nsoil = 9   ! number of soil levels
!================================================================
! parameter lists
!================================================================
!INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=12) ! double precision
INTEGER, PARAMETER :: DBL = SIZEOF(0.)
REAL (KIND=DBL), PARAMETER :: FC = -3.3_DBL 		! field capacity ~-0.33bar=-3.3m
REAL (KIND=DBL), PARAMETER :: WP = -150._DBL 		!wilting point
REAL (KIND=DBL), PARAMETER :: pii = 3.141592653589793 	! acos(-1)

! parameters for stomatal_resistance
REAL (KIND=DBL), PARAMETER :: Rc_max = 5000.0_DBL ! maximum stomatal resistance
REAL (KIND=DBL), PARAMETER :: T_opt = 298.0_DBL   ! optimum temperature for transpiration

! Heat capacity of water
REAL (KIND=DBL), PARAMETER :: cp_water = 4.1796e6_DBL ! [J/m3K]

REAL (KIND=DBL), PARAMETER :: sigma = 5.67e-8_DBL   	! Stefen- Boltzmann constant [W/m^2/K^4]
REAL (KIND=DBL), PARAMETER :: k0 = 0.41_DBL     	! Von Karman constant	

REAL (KIND=DBL), PARAMETER :: z0_soil = 0.0387_DBL     	! baresoil roughness length (m)
REAL (KIND=DBL), PARAMETER :: seaicedepth = 0.5_DBL     ! prescribed depth of seaice (m)
REAL (KIND=DBL), PARAMETER :: tfriz = 273.15_DBL     ! prescribed depth of seaice (m)
REAL (KIND=DBL), PARAMETER :: mws_mx0 = 50._DBL! maximum puddle water storage (mm)

!========================================================
! NAMELIST elements for LAND, variables that are specified in prm
!========================================================
LOGICAL :: dosoilwnudging = .false.
LOGICAL :: dosoiltnudging = .false.
REAL :: tausoil = 86400.
integer:: landtype0 = 0 ! landtype over whole domain
real:: LAI0 = 0. ! default LAI
real:: clay0 = 0. ! uniform clay content (%)
real:: sand0 = 0. ! uniform sand content (%)
real:: sw0 = 0. ! uniform soil wetness (fraction)
real:: st0 = 0. ! uniform soil temperature (k)
logical:: readlandtype = .false. ! read landtype map from file
logical:: readLAI = .false. ! read LAI file
character(80):: landtypefile = ""
character(80):: LAIfile = ""

!================================================================
! variable lists
!================================================================
REAL (KIND=DBL) :: dtn_dbl, dtfactor_dbl 

LOGICAL :: flag_vege_init = .true. ! flag for SLM variables initialization, controlled in main() of SAM
INTEGER :: nstep_restart = 0 

REAL precip_slm(nx,ny)   ! reference precipitation from SAM, used in SRC/precip_fall()
REAL precip_slm2(nx,ny)

REAL , DIMENSION(nx,ny) :: canopt,casq,cast
REAL, DIMENSION(nx,ny) :: energy_bal, energy_bal_c, energy_bal_s

! prognostics variables
REAL (KIND=DBL), DIMENSION(nx,ny) :: t_canop , & ! canopy leaf temperature
				     mw, &	 ! water storage on leaves, kg/m2 = mm of water
				     mws, &	 ! water storage in puddle, kg/m2 = mm of water
				     t_skin	 ! skin temperature, surface level viewed from atmosphere

REAL (KIND=DBL) soilt(nsoil,nx,ny) ! soil temperature at each soil layer
REAL (KIND=DBL)	soilw(nsoil,nx,ny) ! soil wetness at each soil layer

! diagnostic variables 
REAL (KIND=DBL) t_cas(nx,ny) ! canopy air space temperature
REAL (KIND=DBL) q_cas(nx,ny) ! canopy air space specific humidity

REAL (KIND=DBL) soilw_inc(nsoil), t_canop_inc, mw_inc, mws_inc

integer landtype(nx,ny) ! IGBP landtypeA (= 1 to 16)
REAL (KIND=DBL) LAI(nx,ny) ! Leaf Area Index

! vegetation parameter lists : read from input file, RUNDATA/landtype_param_set1*, * is an integer assigned to a specific vegetation
INTEGER, DIMENSION(nx,ny) :: vegetype ! vegetation type (0-baresoil, non zero for vegetated)

REAL (KIND=DBL), DIMENSION(nx,ny) :: vege_YES, & ! zero- for baresoil, 1-vegetated
				     cp_vege, &  ! heat capacity of vegetation
				     z0_sfc, &   ! surface level roughness : for vegetated z0_sfc=z0_vege, 
                                                 !for non vegetated, z0_sfc=z0_soil
				     Khai_L, &	 ! leaf distribution factor
				     phi_1, &    ! 
				     phi_2, &    !
				     IR_emis_vege, & ! vegetation emissivity for infrared radiation
				     IR_emis_soil, & ! soil emissivity for infrared radiation
                                     ztop, &     ! height of vegetation
                                     disp_hgt, & ! displacement height
				     Rgl, & 	 ! parameter for stomatal resistance
				     Rc_min, &   ! minimum stomatal resistance
				     hs_rc, &    ! parameter for stomatal resistance
				     rootL, &    ! root length
				     root_a , &  ! coef a for root density
				     root_b , &  ! coef b for root density
				     precip_extinc, & ! precipitation extinction coef.
				     mw_mx, &    ! max of mw
				     mws_mx, &   ! max of mws
                                     BAI         ! Basal area index, sq.feet/acre

REAL (KIND=DBL) rootF(nsoil,nx,ny)   ! fraction of total root length of a vegetation in each soil layer:
REAL (KIND=DBL)	w_s_FC(nsoil,nx,ny)   ! wetness at field capacity
REAL (KIND=DBL) w_s_WP(nsoil,nx,ny)  ! wetness at wilting point

REAL (KIND=DBL) alpha(nsoil), beta(nsoil)

! soil parameters 
! read from input file 'soil' under the case directory
REAL (KIND=DBL) SAND(nsoil,nx,ny)           ! sand percentage
REAL (KIND=DBL) CLAY(nsoil,nx,ny)           ! clay percentage
REAL (KIND=DBL) s_depth(nsoil,nx,ny)        ! thickness of each soil layer in meter	
REAL (KIND=DBL) soil_relax_hgt(nsoil,nx,ny) ! depth dependent soil relaxation function

REAL (KIND=DBL) tau_soil(nx,ny)   ! soil nudging time in seconds, for dosoiltnudging.or.dosoilwnudging = .true.

! soil parameters calculated from SAND, CLAY
REAL (KIND=DBL), DIMENSION(nsoil,nx,ny) :: sst_capa, & ! heat capacity for soil solids
					  sst_cond, & ! heat conductivity of soil solids
					  poro_soil, & ! saturation moisture content, porosity
					  theta_FC, &  ! volumetric moisture content at field capacity
					  theta_WP, &   ! Volumetric moisture content at wilting point
					  m_pot_sat, & ! moisture potential at saturation
					  Bconst, &  ! B constant
					  ks  ! Hydraulic Conductivity at saturation [m/s]
REAL (KIND=DBL) soilw_inc_tot(nx,ny)

! radiation
REAL (KIND=DBL), DIMENSION(2) :: net_rad, & ! net radiation
				 net_lw, &  
				 net_lwup,&
				 net_lwdn,&
				 net_sw,&
				 net_swup,&
				 net_swdn
! resistances
REAL (KIND=DBL) :: r_a, & 
		   r_b, &
		   r_c, &
		   r_d, &
		   r_soil, &
		   r_litter

! aerodynamic heat conductance
REAL (KIND=DBL) :: cond_heat, &
		   cond_href, &
		   cond_hcnp, &
		   cond_hundercnp

! aerodynamic vapor conductance
REAL (KIND=DBL) :: cond_vapor, &
		   cond_vref, &
		   cond_vcnp, &
		   cond_vundercnp

REAL (KIND=DBL), DIMENSION(nx,ny) :: ustar, &! 
				     tstar!
REAL (KIND=DBL) :: mom_trans_coef, & ! surface momentum transfer coef.
		   heat_trans_coef, & ! surface heat transfer coef.
		   xsi, & ! M-O stability parameter
		   RiB, & ! bulk Richardson number
		   vel_m ! reference level wind speed 

! surface momentum fluxes
REAL (KIND=DBL) :: taux_sfc, tauy_sfc

! 2D statistics variables
REAL (KIND=DBL), DIMENSION(2,nx,ny) :: s_net_swup,& ! net 
				       s_net_lwup,&
				       s_net_rad
REAL (KIND=DBL), DIMENSION(nx,ny) :: s_taux_sfc, & ! surface momentum flux in x directiron
				     s_tauy_sfc, & ! surface momentum flux in y direction
				     s_precip, & !precipitation intercepted by canopy
				     s_drain, & ! drainage rate from canopy
				     s_precip_sfc, & ! preciptation rate at the actual soil surface
				     s_ustar, & ! 
				     s_r_a, & ! resistance Ra
				     s_r_b, & ! resistance Rb
				     s_r_d, & ! resistance Rd
				     s_r_soil, & ! resistance Rsoil
				     s_r_c, & ! stomatal resistance
				     s_shf_canop, & ! sensible heat flux from leaf surface to air
				     s_shf_soil, & ! sensible heat flux from soil to air
				     s_shf_air, & ! sensible heat flux from air to ref. level
				     s_lhf_canop, & ! canopy latent heat flux 
				     s_lhf_soil, & ! soil latent heat flux
				     s_lhf_air, & ! total latent heat flux
				     s_lhf_wet, & ! direct evaporation from water storage on leaves
				     s_lhf_tr, &  ! transpiration
			 	     s_precip_ref,& ! reference level precipitation
				     s_precip_in,& ! precipitation infiltration rate 
				     s_run_off_sfc, & ! surface run off rate
				     s_drainage, & ! bottom soil layer drainage
				     s_tcnp, & ! canopy T
				     s_tcas,&
				     s_qcas,&
				     s_mw, &
				     s_mws, &
				     s_grflux0, &
				     s_RiB, &
                                     s_LAI, &
                                     s_Cveg ! canopy heat capacity

REAL (KIND=DBL), DIMENSION(nsoil,nx,ny):: s_soilw, & ! soil wetness
					  s_soilw_nudge,& ! soil wetness nudging
					  s_soilt_nudge,& ! soil temperature nudging
					  s_soilt  ! soil temperature

! soil hydrology variables
REAL (KIND=DBL) sh_eff_cond(nsoil-1), sh_eff_vel(nsoil-1)
REAL (KIND=DBL) :: drainage_flux ! drainage from the bottom layer

! variables for soil temperature calculations
REAL (KIND=DBL) st_cond(nsoil) ! soil heat conductivity at each node depth
REAL (KIND=DBL) st_capa(nsoil) ! soil heat capacity at each node depth	
REAL (KIND=DBL) st_eff_cond(nsoil-1) ! effective soil thermal conductivity 
REAL (KIND=DBL) :: grflux0  !heat flux at soil top


REAL (KIND=DBL), DIMENSION(nsoil,nx,ny)::s_depth_mm, & ! thickness of each soil layer in mm
					 node_z, &  ! depth of the center of each soil layer
					 interface_z, & ! depth of soil layer interface
					 soilw_nudge, & ! nudging of soilw
					 soilt_nudge, & ! nudging of soilt
					 sw_over_sat ! soilw over saturation

! Sensible heat fluxes
REAL (KIND=DBL), DIMENSION(nx,ny):: shf_canop, & ! flux transfer between canopy and air space [W/m^2]
			   	     shf_soil , & ! flux transfer between soil surface and air space
			   	     shf_air      ! flux transfer between air space and the reference level

! Latent heat fluxes
REAL (KIND=DBL), DIMENSION(nx,ny) :: lhf_canop , & ! flux transfer between canopy and air space W/m2
 			   	     lhf_soil,&
			   	     lhf_air      ! flux transfer between air space to the reference level
REAL (KIND=DBL)::		      evapo_dry , & ! transpiration rate kg/m2/s = mm/s
				      evapo_wet , & ! direct evaporation rate from 'mw'
				      evapo_s ! direct evaporation from soil top

! Storage for initial soilT and soilW for soil nudging
REAL (KIND=DBL) soilt_obs(nsoil,nx,ny) ! profile of temperature to nudge to
REAL (KIND=DBL) soilw_obs(nsoil,nx,ny) ! profile of wetness to nudge to

REAL (KIND=DBL), DIMENSION(nx,ny) :: albedovis_v, albedonir_v ! vis and nir albedos for vegatation
REAL (KIND=DBL), DIMENSION(nx,ny) :: albedovis_s, albedonir_s ! vis and nir albedos for soil

REAL (KIND=DBL) :: precip, & !precipitation rate intercepted by vegetation layer, mm/s
	           drain, &! drainage rate from 1st vegetation layer
		   precip_sfc,& !precipitation rate reaching to soil surface
		   precip_in,& ! rain infiltration rate [kg/m^2/s]
		   wet_canop,& ! wetted fraction of canopy
		   run_off_sfc,& !surface  run off rate [kg/m^2/s=mm/s]
		   cnp_mw_drip ! drip from leaves
REAL (KIND=DBL) :: q_gr,& ! q at ground
		   qsat_canop,&! sat. q for t_canop and p_ref
		   qsat_skin  ! sat. q for soil surface


REAL (KIND=DBL) :: soilvwc
LOGICAL :: VEGYES

CONTAINS

!================================================================================
!================================================================================
!================================================================================

SUBROUTINE slm_init

! Initialize vegetation and soil parameters  from the value read from prm
! Assign soil depth
! Assign vegetation root fractions in each soil layer

IMPLICIT NONE

REAL (KIND=DBL) :: temp1
INTEGER i,j,k,it,jt,nx1,ny1,nsoil1,iz
real(4) fld(1:nx_gl,1:ny_gl)
real(4), allocatable :: zsoil(:), flds(:,:,:)
real(4) mask(1:nx_gl,1:ny_gl)
real tmp(1),tmp1(1)
integer maski(1:nx_gl,1:ny_gl)

 
if(.not.SLM) return

call slm_setparm()

if(readlandtype.and.landtype0.eq.0) then
  open(11,file=trim(landtypefile),status='old',form='unformatted')
  if(masterproc) print*,'reading land_type data from file:',trim(landtypefile)
  read(11) nx1
  read(11) ny1
  if(nx1.ne.nx_gl.or.ny1.ne.ny_gl) then
      if(masterproc) print*,'dimensions of domain in land_type file are ' // &
                           'different from numerical domain sizes: nx=',nx1, &
                            'ny=',ny1,'  Stop...'
      call task_abort()
  end if
  read(11) maski(1:nx_gl,1:ny_gl)
  call task_rank_to_index(rank,it,jt)
  landtype(1:nx,1:ny) = maski(1+it:nx+it,1+jt:ny+jt)
  close (11)
else
  landtype(1:nx,1:ny) = landtype0
end if
if(minval(landtype).lt.0.or.maxval(landtype).gt.16) then
   print*,'landtype values for process rank',rank,' are outside 0-16 range:', &
           'min:',minval(landtype),'max:',maxval(landtype)
   call task_abort()
end if
if(readLAI.and.LAI0.eq.0.) then
  open(11,file=trim(LAIfile),status='old',form='unformatted')
  if(masterproc) print*,'reading LAI data from file:',trim(LAIfile)
  read(11) nx1
  read(11) ny1
  if(nx1.ne.nx_gl.or.ny1.ne.ny_gl) then
      if(masterproc) print*,'dimensions of domain in LAI file are ' // &
                           'different from numerical domain sizes: nx=',nx1, &
                            'ny=',ny1,'  Stop...'
      call task_abort()
  end if
  read(11) mask(1:nx_gl,1:ny_gl)
  call task_rank_to_index(rank,it,jt)
  LAI(1:nx,1:ny) = mask(1+it:nx+it,1+jt:ny+jt)
  close (11)
else
  LAI(:,:) = LAI0
end if
if(minval(LAI).lt.0.or.maxval(LAI).gt.10.) then
   print*,'LAI values for process rank',rank,' are outside 0-10 range:', &
           'min:',minval(landtype),'max:',maxval(landtype)
   call task_abort()
end if

if(landtype0.ne.0.and.landtype0.ne.15.and.landtype0.ne.16.and.landtype0.ne.13.and.LAI0.eq.0.) then
       if(masterproc) print*,'LAI0 cannot be 0 for landtype0',landtype0,' Stop'
       call task_abort()
end if

dtn_dbl = DBLE(dtn)
	
do j=1,ny 
 do i=1,nx

     select case(landtype(i,j))

! properties according to IGBP class
! albedos from http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.531.9279&rep=rep1&type=pdf
! roughness lengthes from http://journals.ametsoc.org/doi/full/10.1175/JHM536.1
! below from CLM (Zhou et al 2003, JGR)
! soil albedos also from CLM type 16


     case (0)  ! water
       
           landmask(i,j) = 0    
           vegetype(i,j) = 0.
 
     case (1) ! evergreen needleaf forest
                      !   albidovis_v(i,j) = 0.013
                      !   albidonir_v(i,j) = 0.184
                      !   albedovis_v(i,j) = 0.0893
                      !   albedonir_v(i,j) = 0.1729
                         albedovis_v(i,j) = 0.04
                         albedonir_v(i,j) = 0.20
                         albedovis_s(i,j) = 0.19
                         albedonir_s(i,j) = 0.38
                         ztop(i,j) = 20.
                         disp_hgt(i,j) = 0.68*ztop(i,j)
                         z0_sfc(i,j) = 1.09
                         Khai_L(i,j) = 0.25
                         rootL(i,j) = 1.5
                         root_a(i,j) = 5.558
                         root_b(i,j) = 2.614
                         Rc_min(i,j) = 100.
                         Rgl(i,j) = 30.
                         hs_rc(i,j) = 54.53
                         BAI(i,j) = 100.
                         vegetype(i,j)=1.
     case (2) ! evergreen broadleaf forest
                      !   albedovis_v(i,j) = 0.017
                      !   albedonir_v(i,j) = 0.230
                      !   albedovis_v(i,j) = 0.0991
                      !   albedonir_v(i,j) = 0.1728
                         albedovis_v(i,j) = 0.04
                         albedonir_v(i,j) = 0.20
                         albedovis_s(i,j) = 0.19
                         albedonir_s(i,j) = 0.38
                         ztop(i,j) = 20.
                         disp_hgt(i,j) = 0.68*ztop(i,j)
                         z0_sfc(i,j) = 2.65
                         Khai_L(i,j) = 0.25
                         rootL(i,j) = 1.5
                         root_a(i,j) = 5.558
                         root_b(i,j) = 2.614
                         Rc_min(i,j) = 150.
                         Rgl(i,j) = 30.
                         hs_rc(i,j) = 54.53
                         BAI(i,j) = 100.
                         vegetype(i,j)=1.
     case (3) ! deciduous needleaf forest
                      !   albedovis_v(i,j) = 0.026
                      !   albedonir_v(i,j) = 0.216
                      !   albedovis_v(i,j) = 0.0893
                      !   albedonir_v(i,j) = 0.1729
                         albedovis_v(i,j) = 0.05
                         albedonir_v(i,j) = 0.23
                         albedovis_s(i,j) = 0.19
                         albedonir_s(i,j) = 0.38
                         ztop(i,j) = 20.
                         disp_hgt(i,j) = 0.68*ztop(i,j)
                         z0_sfc(i,j) = 0.85
                         Khai_L(i,j) = 0.25
                         rootL(i,j) = 1.5
                         root_a(i,j) = 5.558
                         root_b(i,j) = 2.614
                         Rc_min(i,j) = 100.
                         Rgl(i,j) = 30.
                         hs_rc(i,j) = 54.53
                         BAI(i,j) = 100.
                         vegetype(i,j)=1.
     case (4) ! deciduous broadleaf forest
                      !   albedovis_v(i,j) = 0.029
                      !   albedonir_v(i,j) = 0.302
                      !   albedovis_v(i,j) = 0.0961
                      !   albedonir_v(i,j) = 0.2086
                         albedovis_v(i,j) = 0.07
                         albedonir_v(i,j) = 0.24
                         albedovis_s(i,j) = 0.19
                         albedonir_s(i,j) = 0.38
                         ztop(i,j) = 20.
                         disp_hgt(i,j) = 0.68*ztop(i,j)
                         z0_sfc(i,j) = 0.8
                         Khai_L(i,j) = 0.25
                         rootL(i,j) = 1.5
                         root_a(i,j) = 5.558
                         root_b(i,j) = 2.614
                         Rc_min(i,j) = 100.
                         Rgl(i,j) = 30.
                         hs_rc(i,j) = 54.53
                         BAI(i,j) = 100.
                         vegetype(i,j)=1.
     case (5) ! mixed forest
                      !   albedovis_v(i,j) = 0.020
                      !   albedonir_v(i,j) = 0.268
                      !   albedovis_v(i,j) = 0.0978
                      !   albedonir_v(i,j) = 0.1886
                         albedovis_v(i,j) = 0.06
                         albedonir_v(i,j) = 0.24
                         albedovis_s(i,j) = 0.19
                         albedonir_s(i,j) = 0.38
                         ztop(i,j) = 20.
                         disp_hgt(i,j) = 0.68*ztop(i,j)
                         z0_sfc(i,j) = 0.8
                         Khai_L(i,j) = 0.25
                         rootL(i,j) = 1.5
                         root_a(i,j) = 5.558
                         root_b(i,j) = 2.614
                         Rc_min(i,j) = 100.
                         Rgl(i,j) = 30.
                         hs_rc(i,j) = 54.53
                         BAI(i,j) = 100.
                         vegetype(i,j)=1.
     case (6) ! closed shrublands
                      !   albedovis_v(i,j) = 0.040
                      !   albedonir_v(i,j) = 0.219
                      !   albedovis_v(i,j) = 0.1337
                      !   albedonir_v(i,j) = 0.1883
                         albedovis_v(i,j) = 0.07
                         albedonir_v(i,j) = 0.26
                         albedovis_s(i,j) = 0.19
                         albedonir_s(i,j) = 0.38
                         ztop(i,j) = 1.
                         disp_hgt(i,j) = 0.68*ztop(i,j)
                         z0_sfc(i,j) = 0.03
                         Khai_L(i,j) = 0.25
                         rootL(i,j) = 1.5
                         root_a(i,j) = 5.558
                         root_b(i,j) = 2.614
                         Rc_min(i,j) = 100.
                         Rgl(i,j) = 30.
                         hs_rc(i,j) = 54.53
                         BAI(i,j) = 30.
                         vegetype(i,j)=1.
     case (7) ! open shrublands
                      !   albedovis_v(i,j) = 0.034
                      !   albedonir_v(i,j) = 0.243
                      !   albedovis_v(i,j) = 0.2133
                      !   albedonir_v(i,j) = 0.4527
                         albedovis_v(i,j) = 0.14
                         albedonir_v(i,j) = 0.32
                         albedovis_s(i,j) = 0.19
                         albedonir_s(i,j) = 0.38
                         ztop(i,j) = 1.
                         disp_hgt(i,j) = 0.68*ztop(i,j)
                         z0_sfc(i,j) = 0.03
                         Khai_L(i,j) = 0.25
                         rootL(i,j) = 1.5
                         root_a(i,j) = 5.558
                         root_b(i,j) = 2.614
                         Rc_min(i,j) = 50.
                         Rgl(i,j) = 30.
                         hs_rc(i,j) = 54.53
                         BAI(i,j) = 30.
                         vegetype(i,j)=1.
     case (8) ! woody savannas
                      !   albedovis_v(i,j) = 0.022
                      !   albedonir_v(i,j) = 0.255
                      !   albedovis_v(i,j) = 0.1067
                      !   albedonir_v(i,j) = 0.1284
                         albedovis_v(i,j) = 0.06
                         albedonir_v(i,j) = 0.21
                         albedovis_s(i,j) = 0.19
                         albedonir_s(i,j) = 0.38
                         ztop(i,j) = 5.
                         disp_hgt(i,j) = 0.68*ztop(i,j)
                         z0_sfc(i,j) = 0.86
                         Khai_L(i,j) = 0.25
                         rootL(i,j) = 1.5
                         root_a(i,j) = 5.558
                         root_b(i,j) = 2.614
                         Rc_min(i,j) = 100.
                         Rgl(i,j) = 30.
                         hs_rc(i,j) = 54.53
                         BAI(i,j) = 50.
                         vegetype(i,j)=1.
     case (9) ! savannas
                      !   albedovis_v(i,j) = 0.031
                      !   albedonir_v(i,j) = 0.272
                      !   albedovis_v(i,j) = 0.1260
                      !   albedonir_v(i,j) = 0.3059
                         albedovis_v(i,j) = 0.07
                         albedonir_v(i,j) = 0.26
                         albedovis_s(i,j) = 0.19
                         albedonir_s(i,j) = 0.38
                         ztop(i,j) = 5.
                         disp_hgt(i,j) = 0.68*ztop(i,j)
                         z0_sfc(i,j) = 0.86
                         Khai_L(i,j) = 0.25
                         rootL(i,j) = 1.5
                         root_a(i,j) = 5.558
                         root_b(i,j) = 2.614
                         Rc_min(i,j) = 100.
                         Rgl(i,j) = 30.
                         hs_rc(i,j) = 54.53
                         BAI(i,j) = 50.
                         vegetype(i,j)=1.
     case (10) ! grasslands
                      !   albedovis_v(i,j) = 0.045
                      !   albedonir_v(i,j) = 0.284
                      !   albedovis_v(i,j) = 0.1670
                      !   albedonir_v(i,j) = 0.3743
                         albedovis_v(i,j) = 0.07
                         albedonir_v(i,j) = 0.25
                         albedovis_s(i,j) = 0.19
                         albedonir_s(i,j) = 0.38
                         ztop(i,j) = 0.5
                         disp_hgt(i,j) = 0.68*ztop(i,j)
                         z0_sfc(i,j) = 0.08
                         Khai_L(i,j) = -0.3
                         rootL(i,j) = 1.5
                         root_a(i,j) = 5.558
                         root_b(i,j) = 2.614
                         Rc_min(i,j) = 170.
                         Rgl(i,j) = 100.
                         hs_rc(i,j) = 39.18
                         BAI(i,j) = 30.
                         vegetype(i,j)=1.
     case (11) ! permanent wetlands
                      !   albedovis_v(i,j) = 0.019
                      !   albedonir_v(i,j) = 0.193
                      !   albedovis_v(i,j) = 0.0874
                      !   albedonir_v(i,j) = 0.1943
                         albedovis_v(i,j) = 0.06
                         albedonir_v(i,j) = 0.18
                         albedovis_s(i,j) = 0.19
                         albedonir_s(i,j) = 0.38
                         ztop(i,j) = 0.5
                         disp_hgt(i,j) = 0.68*ztop(i,j)
                         z0_sfc(i,j) = 0.04
                         Khai_L(i,j) = -0.3
                         rootL(i,j) = 1.5
                         root_a(i,j) = 5.558
                         root_b(i,j) = 2.614
                         Rgl(i,j) = 0.
                         hs_rc(i,j) = 0.
                         BAI(i,j) = 0.
                         vegetype(i,j)=0.
     case (12) ! croplands
                      !   albedovis_v(i,j) = 0.035
                      !   albedonir_v(i,j) = 0.347
                      !   albedovis_v(i,j) = 0.1109
                      !   albedonir_v(i,j) = 0.2546
                         albedovis_v(i,j) = 0.06
                         albedonir_v(i,j) = 0.24
                         albedovis_s(i,j) = 0.19
                         albedonir_s(i,j) = 0.38
                         ztop(i,j) = 0.5
                         disp_hgt(i,j) = 0.68*ztop(i,j)
                         z0_sfc(i,j) = 0.075
                         Khai_L(i,j) = -0.3
                         rootL(i,j) = 1.5
                         root_a(i,j) = 5.558
                         root_b(i,j) = 2.614
                         Rc_min(i,j) = 100.
                         Rgl(i,j) = 100.
                         hs_rc(i,j) = 39.18
                         BAI(i,j) = 50.
                         vegetype(i,j)=1.
     case (13) ! urban
                      !   albedovis_v(i,j) = 0.031
                      !   albedonir_v(i,j) = 0.234
                      !   albedovis_v(i,j) = 0.0996
                      !   albedonir_v(i,j) = 0.1598
                         albedovis_v(i,j) = 0.06
                         albedonir_v(i,j) = 0.22
                         albedovis_s(i,j) = 0.19
                         albedonir_s(i,j) = 0.38
                         ztop(i,j) = 10.
                         disp_hgt(i,j) = 0.68*ztop(i,j)
                         z0_sfc(i,j) = 1.0
                         Khai_L(i,j) = -0.3
                         rootL(i,j) = 1.5
                         root_a(i,j) = 5.558
                         root_b(i,j) = 2.614
                         Rc_min(i,j) = 0.
                         Rgl(i,j) = 0.
                         hs_rc(i,j) = 0.
                         BAI(i,j) = 0.
                         vegetype(i,j)=0.
     case (14) ! croplands/natural mozaics
                      !   albedovis_v(i,j) = 0.030
                      !   albedonir_v(i,j) = 0.314
                      !   albedovis_v(i,j) = 0.1023
                      !   albedonir_v(i,j) = 0.2440
                         albedovis_v(i,j) = 0.06
                         albedonir_v(i,j) = 0.22
                         albedovis_s(i,j) = 0.19
                         albedonir_s(i,j) = 0.38
                         ztop(i,j) = 0.5
                         disp_hgt(i,j) = 0.68*ztop(i,j)
                         z0_sfc(i,j) = 0.07
                         Khai_L(i,j) = -0.3
                         rootL(i,j) = 1.5
                         root_a(i,j) = 5.558
                         root_b(i,j) = 2.614
                         Rc_min(i,j) = 100.
                         Rgl(i,j) = 100.
                         hs_rc(i,j) = 39.18
                         BAI(i,j) = 50.
                         vegetype(i,j)=1.
     case (15) ! sbow/ice
                         albedovis_v(i,j) = 0. ! actually depends on moisture content
                         albedonir_v(i,j) = 0.
                         albedovis_s(i,j) = 0.95
                         albedonir_s(i,j) = 0.65
                         ztop(i,j) = 0.
                         z0_sfc(i,j) = 0.01
                         disp_hgt(i,j) = 0.
                         Khai_L(i,j) = 0.
                         rootL(i,j) = 0.
                         root_a(i,j) = 0.
                         root_b(i,j) = 0.
                         Rc_min(i,j) = 0.
                         Rgl(i,j) = 0.
                         hs_rc(i,j) = 0.
                         BAI(i,j) = 0.
                         vegetype(i,j)=0.
     case (16) ! baresoil
                         albedovis_v(i,j) = 0. ! actually depends on moisture content
                         albedonir_v(i,j) = 0.
                         albedovis_s(i,j) = 0.19
                         albedonir_s(i,j) = 0.38
                         ztop(i,j) = 0.
                         disp_hgt(i,j) = 0.
                         z0_sfc(i,j) = z0_soil
                         Khai_L(i,j) = 0.
                         rootL(i,j) = 0.
                         root_a(i,j) = 0.
                         root_b(i,j) = 0.
                         Rc_min(i,j) = 0.
                         Rgl(i,j) = 0.
                         hs_rc(i,j) = 0.
                         BAI(i,j) = 0.
                         vegetype(i,j)=0.
     case default
       print*,'landtype is not defined in slm_init. rank=',rank
       call task_abort()
     end select 
 
     if(vegetype(i,j).eq.1..and.LAI(i,j).eq.0.) then
       print*,'LAI is not set for vegetated land point','rank=',rank,'i,j=',i,j
       call task_abort()
     end if
 end do
end do

vege_YES = 1._DBL ! 1-vegetated, 0-non vegetated
WHERE(vegetype.eq.0) vege_YES = 0._DBL
WHERE(vegetype.eq.1) LAI=max(LAI,0.001_DBL) ! set minimum LAI for vegetated land 

if(dompi) then
  tmp1(1) = real(maxval(vegetype))
  call task_max_real(tmp1,tmp,1)
  VEGYES=tmp(1).gt.0
else
  VEGYES=maxval(vegetype).gt.0.
end if

IR_emis_vege = (1._DBL-exp(-1._DBL*LAI))*0.97_DBL    
phi_1 = 0.5_DBL-0.633_DBL*khai_L-0.33_DBL*(khai_L)**2
phi_2 = 0.877_DBL*(1._DBL-2._DBL*phi_1)
precip_extinc = phi_1 + phi_2 ! precipitation extinction coefficient
mw_mx = 0.1_DBL*LAI ! maximum water storage on leaves ; equivalent with 0.1mm of water
mws_mx = mws_mx0 ! maximum water storage in puddle ; 
!============= Initialize soil parameters ====================================!
call init_soil_tw	

IR_emis_soil = 0.98_DBL ! Soil IR emissivity

DO j = 1,ny
 DO i = 1,nx
  if(landmask(i,j).eq.1) then
   DO k = 1,nsoil

   ! soil solids thermal conductivity(Johansen 1975) ! quartz content = SAND content
   ! quartz(=SAND content) thermal conductivity = 7.7 W/mK

    if(SAND(k,i,j).gt.20._DBL) then
 	temp1 = 2.0_DBL   ! thermal conductivity of other minerals [W/mK]
    else
	temp1 = 3.0_DBL
    end if
    sst_cond(k,i,j) = (7.7_DBL**(SAND(k,i,j)/100.0_DBL))*(temp1**(1.0_DBL-(SAND(k,i,j)/100.0_DBL))) 

   ! Calculated following Cosby et al. 1984 ( hydraulic properties) 
   ! hydraulic conductivity at satuation , mm/s
    ks(k,i,j) = (10.0_DBL**(0.0153_DBL*SAND(k,i,j)-0.884_DBL))*(25.4_DBL/3600.0_DBL) ![mm/s] from [inch/hour]

   ! constant B
    Bconst(k,i,j) = 0.159_DBL*CLAY(k,i,j) + 2.91_DBL

   ! porosity (or saturation volumetric water content)
    poro_soil(k,i,j) = -0.00126_DBL*SAND(k,i,j) + 0.489_DBL ! [volume/volume]

   ! moisture potential at saturation, mm
    m_pot_sat(k,i,j) = min(-150.,-10.0_DBL*(10.**(1.88_DBL-0.0131_DBL*SAND(k,i,j)))) ! [mm] from [cm]

   ! soil heat capacity, J/m3/K
   ! [J/m3/K]Following de Vries(1963) using SAND 34% CLAY 63%
    sst_capa(k,i,j) = (2.128_DBL*SAND(k,i,j)+2.385_DBL*CLAY(k,i,j))/(SAND(k,i,j)+CLAY(k,i,j))*(1.e6_DBL) 

   ! volumetric moisture content at field capacity 
   ! field capacity is assumed to be the occasion when hydraulic conductivity is 0.1mm/d
    theta_FC(k,i,j) = poro_soil(k,i,j)*(0.1_DBL/86400._DBL/ks(k,i,j))**(1._DBL/(2._DBL*Bconst(k,i,j)+3._DBL))

   ! volumetric moisture content at wilting point
    theta_WP(k,i,j) = poro_soil(k,i,j) *(-150000._DBL/m_pot_sat(k,i,j))**(-1._DBL/Bconst(k,i,j)) 

    ! soil wetness at field capacity
    w_s_FC(k,i,j) = theta_FC(k,i,j)/poro_soil(k,i,j)

    ! soil wetness at wilting point
    w_s_WP(k,i,j) = theta_WP(k,i,j)/poro_soil(k,i,j)

   END DO 
  end if ! landmask
 END DO 
END DO  

! Calculate fraction of root in each soil layer
call vege_root_init()

! default nudging profiles are just initial 
soilt_obs(:,:,:) = soilt(:,:,:)
soilw_obs(:,:,:) = soilw(:,:,:)


END SUBROUTINE slm_init

!================================================================================
!================================================================================
!================================================================================

SUBROUTINE slm_setparm
implicit none
 	INTEGER ios
	CHARACTER(LEN=80) :: msg
	
	NAMELIST /SLM/dosoilwnudging,dosoiltnudging,tausoil,landtype0,LAI0,clay0,sand0,sw0,st0, &
                      readlandtype,landtypefile,readLAI,LAIfile

	! READ in from prm
	open(55, file='./'//trim(case)//'/prm',status='old',form='formatted')
	read(55,SLM,IOSTAT=ios,IOMSG=msg)	     
	
	if(ios.gt.0) then
		write(*,*) '***** ERROR : bad specification in SLM namelist for the reason of', TRIM(msg)
		call task_abort()
	else if(ios.lt.0) then
		write(*,*) '****** SLM is not activated in prm : make sure that SLM is not being used for ', trim(case), ' case'
	end if
	close(55)
      ! write namelist values out to file for documentation
        if(masterproc) then
        open(unit=55,file='./OUT_STAT/'//trim(case)//'_'//trim(caseid)//'.nml',&
             form='formatted', position='append')
        write (55,nml=SLM)
        write(55,*)
        close(55)
end if


END SUBROUTINE slm_setparm
		
!================================================================================
!================================================================================
!================================================================================

SUBROUTINE init_soil_tw
use vars, only: sstxy,t00
use params, only: tabs_s
IMPLICIT NONE
integer :: dimsoil, nlsf,i,j,k,it,jt,nx1,ny1,nsoil1
real (KIND=DBL) :: tmp
real(4) fld(1:nx_gl,1:ny_gl),sd(nsoil),st(nsoil),sw(nsoil),snd(nsoil),cly(nsoil),srh(nsoil)
!================================================================================================================
! read 'soil' input file in case directory
! allocate soil variables with nsoil dimension read from soil
! assign initial soil wetness and temperature from soil
! In current version, there is no x-y distribution of the soil type in SAND and CLAY percentages.
!================================================================================================================

open(77,file='./'//trim(case)//'/soil',status='old',form='formatted')
read(77,*)
do k=1,nsoil
   read(77,*) sd(k), st(k), sw(k), snd(k), cly(k), srh(k)
end do
close(77)
DO j = 1,ny
DO i = 1,nx
   if(landmask(i,j).eq.1) then
    tau_soil(i,j) = tausoil
    s_depth(:,i,j) = sd(:)
    if(clay0.eq.0) then
     CLAY(:,i,j) = cly(:)
    else
     CLAY(:,i,j) = clay0
    end if
    if(sand0.eq.0) then
     SAND(:,i,j) = snd(:)
    else
     SAND(:,i,j) = sand0
    end if
    if(st0.eq.0) then
     soilt(:,i,j) = st(:)
    else
     soilt(:,i,j) = st0
    end if
    if(sw0.eq.0) then
     soilw(:,i,j) = sw(:)
    else
     soilw(:,i,j) = sw0
    end if
    soil_relax_hgt(:,i,j) = srh(:)
    if(nrestart.eq.0)sstxy(i,j) = soilt(1,i,j) - t00
   else
    if(nrestart.eq.0)sstxy(i,j) = tabs_s - t00
   end if
END DO ! i
END DO ! j

s_depth_mm = s_depth*1000._DBL ! soil depth in mm unit

! Calculate node_z (depth of the center of each soil lyaer)
 node_z(1,:,:) = 0.5_DBL*s_depth(1,:,:)
 DO k = 2,nsoil
   node_z(k,:,:) = 0.5_DBL*s_depth(k,:,:)
   DO j = 1,k-1
      node_z(k,:,:) = node_z(k,:,:) + s_depth(j,:,:)
   END DO
 END DO
! interface_z : soil depth at the interface level
 DO k = 1,nsoil
   interface_z(k,:,:) = node_z(k,:,:) + 0.5_DBL*s_depth(k,:,:)
 END DO


END SUBROUTINE init_soil_tw

!================================================================================
!================================================================================
!================================================================================

SUBROUTINE vege_root_init
! Fraction of total root in each soil layer is determined based on the soil depth and vegetation root parameters
IMPLICIT NONE

INTEGER :: nrootind, i, j,k
REAL (KIND=DBL) :: tot_root_density, new_root_density

! Assign root fraction in each soil layer 
DO i = 1, nx
DO j = 1, ny 
    if(landmask(i,j).eq.1) then
	rootF(1:nsoil,i,j) = 0.
	nrootind = 1 
	DO k = nsoil,2,-1
		if(interface_z(k,i,j).ge.rootL(i,j)) then
			if(interface_z(k-1,i,j).lt.rootL(i,j)) then
				nrootind =k 
			end if
		end if
	END DO
	rootF(nrootind,i,j) = 1._DBL-0.5_DBL*(exp(-1._DBL*root_a(i,j)*rootL(i,j)) &
			    + exp(-1._DBL*root_b(i,j)*rootL(i,j)))
	
	tot_root_density = rootF(nrootind,i,j)
	DO k = 1, nrootind-1
		rootF(k,i,j) = 1._DBL-0.5_DBL*(exp(-1._DBL*root_a(i,j)*interface_z(k,i,j)) &
					+ exp(-1._DBL*root_b(i,j)*interface_z(k,i,j)))  
	END DO
	
	if(nrootind.gt.1) then
	DO k = nrootind, 2,-1
		rootF(k,i,j) = rootF(k,i,j)-rootF(k-1,i,j)
	END DO

	! to ensure total root density equals to 1 
	new_root_density = 0._DBL
	DO k = 1,nrootind
		rootF(k,i,j) = rootF(k,i,j)/tot_root_density
			new_root_density = new_root_density + rootF(k,i,j)
	END DO
	end if
   end if
END DO ! j
END DO	 ! i
END SUBROUTINE vege_root_init

!================================================================================
!================================================================================
!================================================================================

SUBROUTINE init_slm_vars(tr, ts, qr)
IMPLICIT NONE

REAL, INTENT(IN), DIMENSION(nx,ny) :: tr, ts, qr
INTEGER :: i,j
INTEGER :: it, jt
real, external :: qsatw,qsati

soilw_nudge = 0._DBL
soilt_nudge = 0._DBL
r_a = 0._DBL
r_b = 0._DBL
r_d = 0._DBL
r_c = 0._DBL
DO j = 1, ny
DO i = 1,nx
    if(landmask(i,j).eq.1) then
        ! soil surface specific humidity
        if(soilt(1,i,j).ge.tfriz) then
          q_gr = qsatw(ts(i,j), pres0)*soilw(1,i,j)
        else
          q_gr = qsati(ts(i,j), pres0)*soilw(1,i,j)
        end if
        ! canopy temperature : initialize as (ref. level temperature + soil surface temperature)/2.
        t_canop(i,j) = (DBLE(tr(i,j))+DBLE(ts(i,j)))*0.5_DBL
        ! canopy air space temperature 
        t_cas(i,j)   = (DBLE(tr(i,j))+DBLE(ts(i,j)))*0.5_DBL
        ! amount of water held on vegetation 
        mw(i,j) = 0.0_DBL
        mws(i,j) = 0.0_DBL
        ! specific humidity in canopy air space
        q_cas(i,j)   = (DBLE(qr(i,j))+q_gr)*0.5_DBL
        ustar(i,j) = 0.1_DBL
        tstar(i,j) = 0._DBL
   end if
END DO
END DO
END SUBROUTINE init_slm_vars

!====================================================================================
!====================================================================================
SUBROUTINE collect_2D_stat_vars(i,j)
IMPLICIT NONE
integer,intent(in) :: i,j
integer ::k
dtfactor_dbl = DBLE(dtfactor)

	s_taux_sfc(i,j) = s_taux_sfc(i,j) + taux_sfc*dtfactor_dbl
	s_tauy_sfc(i,j) = s_tauy_sfc(i,j) + tauy_sfc*dtfactor_dbl

	s_precip(i,j) = s_precip(i,j) + precip*dtfactor_dbl
	s_drain(i,j) = s_drain(i,j) + drain*dtfactor_dbl
	s_precip_sfc(i,j) = s_precip_sfc(i,j) + precip_sfc*dtfactor_dbl
	
	s_net_swup(:,i,j) = s_net_swup(:,i,j) + net_swup*dtfactor_dbl
	s_net_lwup(:,i,j) = s_net_lwup(:,i,j) + net_lwup*dtfactor_dbl
	s_net_rad(:,i,j) = s_net_rad(:,i,j) + net_rad*dtfactor_dbl
	
	s_RiB(i,j) = s_RiB(i,j) + RiB*dtfactor_dbl
	s_ustar(i,j) = s_ustar(i,j) + ustar(i,j)*dtfactor_dbl
	
	s_LAI(i,j) = s_LAI(i,j) + LAI(i,j)*dtfactor_dbl
        s_Cveg(i,j) = s_Cveg(i,j) + Cp_vege(i,j)*dtfactor_dbl

	s_r_a(i,j) = s_r_a(i,j) + r_a*dtfactor_dbl
	s_r_b(i,j) = s_r_b(i,j) + r_b*dtfactor_dbl
	s_r_c(i,j) = s_r_c(i,j) + r_c*dtfactor_dbl
	s_r_d(i,j) = s_r_d(i,j) + r_d*dtfactor_dbl
	s_r_soil(i,j) = s_r_soil(i,j) + r_soil*dtfactor_dbl

	s_shf_canop(i,j) = s_shf_canop(i,j) + shf_canop(i,j)*dtfactor_dbl
	s_shf_soil(i,j) = s_shf_soil(i,j) + shf_soil(i,j)*dtfactor_dbl
	s_shf_air(i,j) = s_shf_air(i,j) + shf_air(i,j)*dtfactor_dbl

	s_lhf_canop(i,j) = s_lhf_canop(i,j) + lhf_canop(i,j)*dtfactor_dbl
	s_lhf_soil(i,j) = s_lhf_soil(i,j) + lhf_soil(i,j)*dtfactor_dbl
	s_lhf_air(i,j) = s_lhf_air(i,j) + lhf_air(i,j)*dtfactor_dbl
	s_lhf_wet(i,j) = s_lhf_wet(i,j) + evapo_wet*DBLE(lcond*dtfactor)
	s_lhf_tr(i,j) = s_lhf_tr(i,j) + evapo_dry*DBLE(lcond*dtfactor)

	s_tcnp(i,j)	= s_tcnp(i,j) + t_canop(i,j) * dtfactor_dbl
	s_tcas(i,j)	= s_tcas(i,j) + t_cas(i,j) * dtfactor_dbl
	s_qcas(i,j)	= s_qcas(i,j) + q_cas(i,j) * dtfactor_dbl
	s_mw(i,j)		= s_mw(i,j) + mw(i,j) * dtfactor_dbl
	s_mws(i,j)		= s_mws(i,j) + mws(i,j) * dtfactor_dbl
	s_precip_in(i,j)	= s_precip_in(i,j) + precip_in*dtfactor_dbl
	s_run_off_sfc(i,j)	= s_run_off_sfc(i,j) + run_off_sfc*dtfactor_dbl
	s_drainage(i,j) = s_drainage(i,j) + drainage_flux*dtfactor_dbl
	

	DO k = 1,nsoil
		s_soilw(k,i,j) = s_soilw(k,i,j) + soilw(k,i,j)*dtfactor_dbl
		s_soilt(k,i,j) = s_soilt(k,i,j) + soilt(k,i,j)*dtfactor_dbl
		s_soilw_nudge(k,i,j) = s_soilw_nudge(k,i,j) + soilw_nudge(k,i,j)*dtfactor_dbl
		s_soilt_nudge(k,i,j) = s_soilt_nudge(k,i,j) + soilt_nudge(k,i,j)*dtfactor_dbl
	END DO
	s_grflux0(i,j) = s_grflux0(i,j) + grflux0*dtfactor_dbl
        
END SUBROUTINE collect_2D_stat_vars 

!================================================================================
!================================================================================
!================================================================================


SUBROUTINE slm_printout


implicit none

integer :: nunit, ncnt, i,j,k
if(masterproc)print*, '========== SLM parameters ================'
if(masterproc.and..not.readlandtype)print*,'landtype0=',landtype0
if(masterproc.and..not.readLAI) print*,'LAI0=',LAI0

call fminmax_print_slm('landtype:', real(landtype), 1,nx, 1,ny)
if(VEGYES)call fminmax_print_slm('LAI:', real(LAI), 1,nx, 1,ny)
if(VEGYES)call fminmax_print_slm('albedovis_v:', real(albedovis_v), 1,nx, 1,ny)
if(VEGYES)call fminmax_print_slm('albedonir_v:', real(albedonir_v), 1,nx, 1,ny)
call fminmax_print_slm('albedovis_s:', real(albedovis_s), 1,nx, 1,ny)
call fminmax_print_slm('albedonir_s:', real(albedonir_s), 1,nx, 1,ny)
if(VEGYES)call fminmax_print_slm('ztop:', real(ztop), 1,nx, 1,ny)
if(VEGYES)call fminmax_print_slm('disp_hgt:', real(disp_hgt), 1,nx, 1,ny,.true.)
call fminmax_print_slm('z0_sfc:', real(z0_sfc), 1,nx, 1,ny)
if(VEGYES)call fminmax_print_slm('Khai_L:', real(Khai_L), 1,nx, 1,ny,.true.)
if(VEGYES)call fminmax_print_slm('rootL:', real(rootL), 1,nx, 1,ny,.true.)
if(VEGYES)call fminmax_print_slm('root_a:', real(root_a), 1,nx, 1,ny,.true.)
if(VEGYES)call fminmax_print_slm('root_b:', real(root_b), 1,nx, 1,ny,.true.)
if(VEGYES)call fminmax_print_slm('Rc_min:', real(Rc_min), 1,nx, 1,ny,.true.)
if(VEGYES)call fminmax_print_slm('Rgl:', real(Rgl), 1,nx, 1,ny,.true.)
if(VEGYES)call fminmax_print_slm('hs_rc:', real(hs_rc), 1,nx, 1,ny,.true.)
if(VEGYES)call fminmax_print_slm('BAI:', real(BAI), 1,nx, 1,ny,.true.)
if(masterproc)print*,'number of soil levels:',nsoil
do k=1,nsoil
end do
do k=1,nsoil
  call fminmax_print_slm('layer_depth:', real(node_z(k,:,:)), 1,nx, 1,ny)
end do
do k=1,nsoil
  call fminmax_print_slm('layer_thickness:', real(s_depth(k,:,:)), 1,nx, 1,ny)
end do
do k=1,nsoil
  call fminmax_print_slm('layer_interface:', real(interface_z(k,:,:)), 1,nx, 1,ny)
end do
do k=1,nsoil
  call fminmax_print_slm('soilt:', real(soilt(k,:,:)), 1,nx, 1,ny)
end do
do k=1,nsoil
  call fminmax_print_slm('soilw:', real(soilw(k,:,:)), 1,nx, 1,ny)
end do
do k=1,nsoil
  if(masterproc)print*,'soil level',k
  call fminmax_print_slm('CLAY:', real(CLAY(k,:,:)), 1,nx, 1,ny)
  call fminmax_print_slm('SAND:', real(SAND(k,:,:)), 1,nx, 1,ny)
  call fminmax_print_slm('sst_cond:', real(sst_cond(k,:,:)), 1,nx, 1,ny)
  call fminmax_print_slm('poro_soil:', real(poro_soil(k,:,:)), 1,nx, 1,ny)
  call fminmax_print_slm('m_pot_sat:', real(m_pot_sat(k,:,:)), 1,nx, 1,ny)
  call fminmax_print_slm('sst_capa:', real(sst_capa(k,:,:)), 1,nx, 1,ny)
  call fminmax_print_slm('theta_FC:', real(theta_FC(k,:,:)), 1,nx, 1,ny)
  call fminmax_print_slm('theta_WP:', real(theta_WP(k,:,:)), 1,nx, 1,ny)
end do

if(masterproc)print*, 'soil moisture nudging:', dosoilwnudging
if(masterproc)print*, 'soil temperature nudging:', dosoiltnudging
if(masterproc)print*, '----------------------------------------------------------'

end subroutine slm_printout

!================================================================================
!================================================================================
!================================================================================

SUBROUTINE slm_stat_2Dinit
IMPLICIT NONE

 s_net_swup = 0. 
 s_net_lwup = 0.
 s_net_rad = 0.
 s_taux_sfc = 0. ! surface momentum flux in x directiron
 s_tauy_sfc = 0. ! surface momentum flux in y direction
 s_precip = 0.  !precipitation intercepted by canopy
 s_drain = 0.  ! drainage rate from canopy
 s_precip_sfc = 0.  ! preciptation rate at the actual soil surface
 s_RiB = 0. ! Bulk Richardson parameter
 s_ustar = 0.  ! 
 s_r_a = 0.  ! resistance Ra
 s_r_b = 0. ! resistance Rb
 s_r_d = 0. ! resistance Rd
 s_r_soil = 0. ! resistance Rsoil
 s_r_c = 0. ! stomatal resistance
 s_LAI = 0.  
 s_Cveg = 0. 
 s_shf_canop = 0. ! sensible heat flux from leaf surface to air
 s_shf_soil = 0. ! sensible heat flux from soil to air
 s_shf_air = 0. ! sensible heat flux from air to ref. level
 s_lhf_canop = 0. ! canopy latent heat flux 
 s_lhf_soil = 0. ! soil latent heat flux
 s_lhf_air = 0. ! total latent heat flux
 s_lhf_wet = 0. ! direct evaporation from water storage on leaves
 s_lhf_tr=0.  ! transpiration
 s_precip_in = 0. ! precipitation infiltration rate 
 s_precip_ref = 0. ! precipitation at reference level 
 s_run_off_sfc = 0. ! surface run off rate
 s_drainage = 0. ! bottom soil layer drainage
 s_tcnp = 0. ! canopy T
 s_tcas = 0.
 s_qcas = 0.
 s_mw = 0.
 s_mws = 0.
 
  s_soilw = 0. ! soil wetness
  s_soilt = 0. ! soil temperature
 s_soilw_nudge = 0. ! soil wetness nudging
 s_soilt_nudge = 0. ! soil temperature nudging
 s_soilt =0. ! soil temperature
 s_grflux0 = 0. ! soil heat flux
END SUBROUTINE slm_stat_2Dinit

!================================================================================
!================================================================================
!================================================================================

SUBROUTINE fminmax_print_slm(name,f,dimx1,dimx2,dimy1,dimy2,vege)
implicit none
logical,optional,intent(in):: vege 
integer, intent(in):: dimx1, dimx2, dimy1, dimy2
real, intent(in):: f(dimx1:dimx2, dimy1:dimy2)
real :: fmn, fmx
character(len=*),intent(in) :: name
real ::fmin(1),fmax(1),fff(1)
integer ::i,j
if(dimx2.eq.1.and.dimy2.eq.1) then
   fmn = f(1,1)	
   fmx = f(1,1)	
else
   fmn = 1.e30
   fmx =-1.e30
end if 
do j=1,ny
  do i=1,nx
     if(landmask(i,j).eq.1) then
	if(.not.present(vege)) then ! find min-max over entire land domain
     		fmn = min(fmn,f(i,j))
     		fmx = max(fmx,f(i,j))
	else 
		if(vegetype(i,j).ne.0) then ! fin min-max over the only domain which has a vegetation layer
     			fmn = min(fmn,f(i,j))
     			fmx = max(fmx,f(i,j))
		end if ! vegetype.ne.0
	end if ! .not.present(vege)
     end if
  end do
enddo
fmin(1) = 1.e30
fmax(1) =-1.e30
fmin(1) = min(fmin(1),fmn)
fmax(1) = max(fmax(1),fmx)
	
if(dompi) then
  fff(1)=fmax(1)
  call task_max_real(fff(1),fmax(1),1)
  fff(1)=fmin(1)
  call task_min_real(fff(1),fmin(1),1)
end if
if(masterproc) print *,name,fmin,fmax
END SUBROUTINE fminmax_print_slm

!================================================================================
!================================================================================
!================================================================================

SUBROUTINE slm_stepout
       
        if(masterproc)  print*, '========== SLM : OVER LAND ============'
        if(VEGYES)call fminmax_print_slm('canopT:', canopt, 1,nx, 1,ny,.true.)
        if(VEGYES)call fminmax_print_slm('casT:', cast, 1,nx, 1,ny,.true.)
        if(VEGYES)call fminmax_print_slm('casQ:', casq, 1,nx, 1,ny,.true.)
        call fminmax_print_slm('shf_soil:',real(shf_soil),1,nx,1,ny)
        if(VEGYES)call fminmax_print_slm('shf_canop:',real(shf_canop),1,nx,1,ny,.true.)
        call fminmax_print_slm('total shf:',real(shf_air),1,nx,1,ny)
        call fminmax_print_slm('lhf_soil:',real(lhf_soil),1,nx,1,ny)
        if(VEGYES)call fminmax_print_slm('lhf_canop:',real(lhf_canop),1,nx,1,ny,.true.)
        call fminmax_print_slm('total lhf:',real(lhf_air),1,nx,1,ny)
        call fminmax_print_slm('soil T top:',real(soilt(1,:,:)),1,nx,1,ny)
        call fminmax_print_slm('soil T deep:',real(soilt(nsoil,:,:)),1,nx,1,ny)
        call fminmax_print_slm('soil wetness top:',real(soilw(1,:,:)),1,nx,1,ny)
        call fminmax_print_slm('soil wetness deep:',real(soilw(nsoil,:,:)),1,nx,1,ny)
        if(masterproc) print*, '========== SLM ============'
END SUBROUTINE slm_stepout

!================================================================================
!================================================================================
!================================================================================

REAL (KIND=DBL) FUNCTION fh_calc(t,  mps, sw, B)
implicit none
real (KIND=DBL) :: t,mps, sw, B
real (KIND=DBL) :: moist_pot1
	moist_pot1 = mps/(max(0.0001,sw)**B)/1000._DBL
	moist_pot1 = max(-100000._DBL, moist_pot1)
	fh_calc = min(1._DBL, exp(moist_pot1*ggr/rv/t))
END FUNCTION fh_calc

END MODULE slm_vars
