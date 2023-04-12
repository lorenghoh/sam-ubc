module microphysics
!
!====================================================
!
! Main interface to P3 microphysics
! - original implementation by Guangxing Lin, PNNL
! - process rate outputs added by Blaz Gasparini, UW
! - consistency of rho and thermodynamic parameters with SAM, Peter Blossey, UW
! - control of P3 options through MICRO_P3 namelist, Blaz Gasparini and Peter Blossey, UW
! - update of P3 to WRFV4, Peter Blossey, UW
!=====================================================

!bloss: Gather all use statements at top of module.
use grid, only: nx, ny, nz, nzm, & ! grid dimensions; nzm=nz-1 - # of levels for all scalars
                dimx1_s,dimx2_s,dimy1_s,dimy2_s, & ! actual scalar-array dimensions
                rank, masterproc, dostatis, &  !BG a la Bloss
                dz, adz, nstep, nprint, z, zi, pres, &
                nrestart, case, caseid, &
                compute_reffc, compute_reffi, & !bloss compute_reffl, &
                nstatis, icycle, dt, dtn
use params, only: doprecip, docloud, fac_cond, fac_sub, lcond
use vars, only: t, tabs, w, qv, qcl, qci, qpl, qpi, &
     fluxbq, fluxtq, precsfc, prec_xy, precflux, total_water_prec, &
     rho, rhow, gamaz, q0, tabs0, t0, &
     tlat, tlatqi, qpfall, &
     condavg_mask, ncondavg, condavgname, condavglongname, precinst
use micro_params !BG added
use MODULE_MP_P3, only : p3_init, p3_main, n_qiType

implicit none

logical :: isallocatedMICRO = .false.

!bloss: 9 prognostic variables, 3 diagnostic variables (tabs_old, qv_old, ice_lambda)
integer, parameter :: nmicro_fields = 12   ! total number of prognostic water vars

! Now we can allocate the required memory for all the prognostic microphysics arrays:

real micro_field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, nmicro_fields)

!bloss: Move these to micro_params.f90, so that they can be used both here
!   and in module_mp_p3.f90.
!bloss! indices of water quantities in micro_field, e.g. qv = micro_field(:,:,:,iqv)
!bloss integer :: iqv, iqcl, incl, inci, inr, iqr, iqit, iqir, iqib

! iqv= 1  (vapor + cloud liq) mass mixing ratio [kg H2O / kg dry air]
!  incl = 2 ! cloud water number mixing ratio  (#/kg dry air)
!  iqr = 3   ! rain mass mixing ratio (kg/kg)
!  inr = 4   ! rain number mixing ratio (#/kg)
!  iqit = 5  ! total ice mass mixing raio (kg/kg) =iqci in Bloss (?)
!  inci = 6  ! total ice number mixing ratio (#/kg)
!  iqir = 7  ! rime ice mass mixing raio (kg/kg)
!  iqib = 8  ! rime ice volume mixing ratio (m^-3/kg)
!  iqcl = 9  ! liquid cloud mass mixing ratio (kg H2O/kg dry air)
!  iqv_old = 10 ! water vapor from end of last time step
!  itabs_old = 10 ! temperature from end of last time step

! For many reasons, for example, to be able to compute water budget, we may need to
! know which variables among the prognostic ones represent water mixing ratio, regardless
! of water species. We use a simple array of flags with 1 marking the water mass
! variable:

!bloss: 
!  - Do not include qcl (index 9) in water mass since it is included in index 1 (vapor + cloud).
!  - Do not include qir (index 7) in water mass since it is included in index 5 (total ice mass).
integer, parameter :: flag_wmass(nmicro_fields) = (/1,0,1,0,1,0,0,0,0,0,0,0/) 

! To implement large-scale forcing, surface fluxes, etc, SAM needs to know
! which variable has a water vapor information. In our example, it is variable #3:

integer, parameter :: index_water_vapor = 1 ! index for variable that contains water vapor
logical, parameter :: is_water_vapor(nmicro_fields) = (/.true.,.false.,.false., & !bloss
                                                        .false.,.false.,.false., &
                                                        .false.,.false.,.false., &
                                                        .false.,.false.,.false. /)

! Now, we need to specify which variables describe precipitation. This is needed because
! SAM has two logical flags to deal with microphysics proceses - docloud and doprecip.
! docloud set to true means that condensation/sublimation processes are allowed to
! form clouds. However, the possibility of rain, snow, heil, etc., is controled by
! a second flag: doprecip. If doprecip=.false. than no precipitation is allowed, hence
! no advection, diffusion, and fallout of corresponding variables should be done;
! therefore, SAM needs an array of flags that mark the prognostic variables which
! only make sense when doprecip=.true. :

integer, parameter :: flag_precip(nmicro_fields) = (/0,0,1,1,1,1,1,1,1,0,0,0/)

!bloss(2020-11): If flag_advect(n)==1, then micro_field(:,:,:,n) will be
!   advected and diffused.  P3 requires some non-advected fields (tabs_old, qv_old)
!   to be saved for the next timestep and including them in micro_field makes that work
!   properly across restarts.  I also added lambda_ice (not advected) as an additional output.
!bloss: add diagnostic variables tabs_old, qv_old to micro_field()
!   These are _not_ advected and diffused and so have flag_advect == 0.
integer, parameter :: flag_advect(nmicro_fields) = (/1,1,1,1,1,1,1,1,1,0,0,0/)

! Sometimes, a cloud ice (or even cloud water) is allowed to be a subject of
! gravitational sedimentation, usually quite slow compared to the precipitation
! drops. SAM calls a special routine, ice_fall() that computes sedimentation of cloud ice.
! However, it is a rudiment from SAM's original single-moment microphysics.
! Instead, you may want to handle sedimentation of cloud water/ice yourself similarly
! to precipitation variables. In this case, set the index for falling cloud ice to -1, which
! means that no default ice mixing ration sedimentation is done.

integer, parameter :: index_cloud_ice = -1   ! index for cloud ice (sedimentation)

! The following arrays are needed to hold the turbulent surface and domain-top fluxes
! for the microphysics prognostic variables:

real fluxbmk (nx, ny, 1:nmicro_fields) ! surface fluxes
real fluxtmk (nx, ny, 1:nmicro_fields) ! top boundary fluxes

! these arrays are needed for output statistics from advection and diffusion routines:

real mkwle(nz,1:nmicro_fields)  ! resolved vertical flux
real mkwsb(nz,1:nmicro_fields)  ! SGS vertical flux
real mkadv(nz,1:nmicro_fields)  ! tendency due to vertical advection
real mkdiff(nz,1:nmicro_fields)  ! tendency due to vertical diffusion
real mklsadv(nz,1:nmicro_fields) ! tendency due to large-scale vertical advection

!! statistical arrays
real mtend(nzm,1:nmicro_fields) ! tendency due to microphysical processes
real stend(nzm,1:nmicro_fields) ! tendency due to sedimentation
real mksed(nzm,1:nmicro_fields) ! ! sedimentation vertical flux
real mk0(nzm,1:nmicro_fields) ! horizontal average of each array in micro_fields !bloss
real mk_ref(nzm,1:nmicro_fields) ! bloss: placeholder.  Could hold reference microphysical profiles
!real micro_proc_rates_stat(nzm,1:nmicro_proc) !BG make column out of the 3d proc rates


! number of fields output from micro_write_fields2D, micro_write_fields3D
integer :: nfields2D_micro=0
integer :: nfields3D_micro=0

!------------------------------------------------------------------


! You may also want to have some additional, diagnostic, arrays;
! for example, total nonprecipitating cloud water, etc:

!bloss/duplicates qcl(:,:,:), qci(:,:,:) and micro_field(:,:,:,iqcl), so not needed
!bloss real qn(nx,ny,nzm)  ! cloud condensate (liquid + ice)
!bloss real cloudliq(nx,ny,nzm)  ! cloud condensate (liquid ) !

character(len=1024), parameter :: lookup_file_dir = './RUNDATA'    ! the directory of lookup tables, should be put in the namelist later

!integer :: n_diag_2d ! the number of 2-D diagnositic variables output in P3_MAIN subroutine
!integer :: n_diag_3d ! the number of 3-D diagnositic variables output in P3_MAIN subroutine
!logical :: log_predictNc !logical variable,  .T. (.F.) for prediction (specification) of Nc
!logical :: typeDiags_ON   !logical variable, for diagnostic hydrometeor/precip rate types
!character(len=1024) :: model    ! type of model that is used in P3_MAIN subroutine
logical :: first              !BG ???

real effR_ice(nx,ny,nzm)  ! effective radius for ice cloud in micro meter

!placeholder variables that allow Morrison and Thompson microphysics
!  couple with the radiation.  Not needed here since only a single
!  ice species is radiatively active
logical :: dosnow_radiatively_active = .false. !P3 microphysics does not predict snow seperately
logical :: dorrtm_cloud_optics_from_effrad_LegacyOption = .false.
! If true, use Dge computed within microphysics for computing radiative properties
!   of snow and cloud ice.  Dge=Generalized effective size
!   See Fu (1996, J. Climate, 9, 2058-2082).
logical :: reff_ice_holds_Dge = .false.

real, allocatable, dimension(:,:,:) :: reffc, reffi, reffs, reffr
real, allocatable, dimension(:,:,:) :: &
     CloudLiquidMassMixingRatio, CloudLiquidGammaExponent, CloudLiquidLambda, &
     CloudIceMassMixingRatio, SnowMassMixingRatio
real, allocatable, dimension (:,:,:,:) :: micro_proc_rates !BG added

real, allocatable, dimension (:,:) ::  trtau, qxoeffr !BG approx optical depth

real :: tmpw_3d(nx,ny,nzm) !BG added as used in 2 routines

! arrays with names/units for microphysical outputs in statistics.
character*3, allocatable, dimension(:) :: mkname
character*80, allocatable, dimension(:) :: mklongname
character*10, allocatable, dimension(:) :: mkunits
real, allocatable, dimension(:) :: mkoutputscale

integer, allocatable, dimension(:) :: flag_number
integer, allocatable, dimension(:) :: flag_micro3Dout
integer, allocatable, dimension(:) :: flag_micro2DZout

CONTAINS

! Below are the required subroutines and functions that you need to fill in.

!copied from Bloss
!----------------------------------------------------------------------
function micro_scheme_name()
  character(len=32) :: micro_scheme_name
  ! Return the scheme name, normally the same as the directory name with leading
  ! "MICRO_" removed
  micro_scheme_name = "p3"
end function   micro_scheme_name
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!!! Read microphysics options from prm file
subroutine micro_setparm()
  ! read any options in from prm file here using a namelist named
  !   after your microphysics routine, e.g. MICRO_P3 for the
  !   P3 microphysics.
!bloss  use vars !BG added
  integer ierr, ios, ios_missing_namelist, place_holder
  character(len=32) :: radname
!BG a la' Bloss
   NAMELIST /MICRO_P3/ &
      iWarmRainScheme, & ! choice of warm rain microphysics
      log_predictNc,  &  ! logic variable,  .T. (.F.) for prediction
      Nc0,            &  ! initial/specified cloud droplet number conc (#/cm3)
      dofix_pgam, pgam_fixed, & ! fix value of gamma exponent for cloud DSD (default = 10.3)
      aerosol_mode1_radius, & ! if log_predictNc==true, the cloud droplets will be activated
      aerosol_mode1_sigmag, & !   from one aerosol mode whose properties can be input here.
      aerosol_mode1_number, & !   If the activation code is extended, a second mode will be
      aerosol_mode2_radius, & !   used as well.  radius in m, sigmag dimensionless, 
      aerosol_mode2_sigmag, & !   number in kg^(-1).
      aerosol_mode2_number, &
      nCat,           &  ! number of free ice categories
      MaxTotalIceNumber, &  ! maximum number concentration for all ice categories combined.
      typeDiags_ON,   &  ! logic variable, for diagnostic hydrometeor/precip
      model,          &  !
      n_diag_2d,      &  ! the number of 2-D diagnositic variables output in P3_MAIN subroutine
      n_diag_3d,      &  ! the number of 3-D diagnositic variables output in P3_MAIN subroutine
      douse_reffc, douse_reffi, & ! leave this as true (the default) except for testing
      do_output_micro_process_rates

   !bloss: Create dummy namelist, so that we can figure out error code
   !       for a mising namelist.  This lets us differentiate between
   !       missing namelists and those with an error within the namelist.
   NAMELIST /BNCUIODSBJCB/ place_holder

   !bloss(2015-02): Default values for namelist variables moved to micro_params.f90

  !----------------------------------
  !  Read namelist for microphysics options from prm file:
  !------------
  open(55,file='./'//trim(case)//'/prm', status='old',form='formatted')

  !bloss: get error code for missing namelist (by giving the name for
  !       a namelist that doesn't exist in the prm file).
  read (UNIT=55,NML=BNCUIODSBJCB,IOSTAT=ios_missing_namelist)
  rewind(55) !note that one must rewind before searching for new namelists

  !bloss: read in MICRO_P3 namelist
  read (55,MICRO_P3,IOSTAT=ios)

  if (ios.ne.0) then
     !namelist error checking
     if(ios.ne.ios_missing_namelist) then
        write(*,*) '****** ERROR: bad specification in MICRO_P3 namelist'
        rewind(55)
        read (55,MICRO_P3) ! this should give a useful error message
        call task_abort()
     elseif(masterproc) then
        write(*,*) '****************************************************'
        write(*,*) '****** No MICRO_P3 namelist in prm file *********'
        write(*,*) '****************************************************'
     end if
  end if
  close(55)

  if(nCat.gt.1) then
    write(*,*) '=============================================================='
    write(*,*) 'Error in P3 microphysics:'
    write(*,*) ' SAM is not set up to run with more than one ice category yet.'
    write(*,*) '=============================================================='
    call task_abort()
  end if

   ! write namelist values out to file for documentation
   if(masterproc) then
      open(unit=55,file='./OUT_STAT/'//trim(case)//'_'//trim(caseid)//'.nml', form='formatted', position='append')
      write (unit=55,nml=MICRO_P3,IOSTAT=ios)
      write(55,*) ' '
      close(unit=55)
   end if
!BG end copy Bloss###################################

  iqv = 1 ! water vapor + liquid cloud water mass mixing ratio (kg H2O/kg dry air)
  !iqcl = 2 ! cloud water mass mixing ratio (kg /kg)
  incl = 2 ! cloud water number mixing ratio  (#/kg dry air)
  iqr = 3   ! rain mass mixing ratio (kg/kg)
  inr = 4   ! rain number mixing ratio (#/kg)
  iqit = 5  ! total ice mass mixing raio [kg H2O / kg dry air] =iqci in Bloss
  inci = 6  ! total ice number mixing ratio (#/kg)
  iqir = 7  ! rime ice mass mixing raio (kg/kg)
  iqib = 8  ! rime ice volume mixing ratio (m^-3/kg)
  iqcl = 9   ! liquid cloud water mass mixing ratio (kg/kg)

  ! include these two variables in the micro_field array even though they
  !   should not be advected.  This is done to preserve consistent restarts
  !   as these variables are used to compute the saturation ratio at the beginning 
  !   of the microphysics.
  iqv_old = 10
  itabs_old = 11
  ilambda_ice = 12

  n_diag_2d = 1
  n_diag_3d = 1

  nmicro_proc = nmicro_process_rates !BG no need for renaming if always all micro proc

  if(.not.isallocatedMICRO) then
    ! allocate microphysical variables
    allocate(reffc(nx,ny,nzm), reffi(nx,ny,nzm), reffs(nx,ny,nzm), &
         flag_number(nmicro_fields), flag_micro3Dout(nmicro_fields), &
         flag_micro2DZout(nmicro_fields), &
         mkname(nmicro_fields), mklongname(nmicro_fields), &
         mkunits(nmicro_fields), mkoutputscale(nmicro_fields), &
         micro_proc_rates(nx,ny,nzm,nmicro_proc), &
         trtau(nzm,nmicro_fields),qxoeffr(nzm,nmicro_fields), &    
         STAT=ierr)
    if(ierr.ne.0) then
      write(*,*) 'Failed to allocate microphysical arrays on proc ', rank
      call task_abort()
    end if

   !BG in namelist Nc0 = 100. ! default droplet number concentration

    !bloss/duplicates micro_field(:,:,:,iqcl)
    !bloss  cloudliq =0.0
    micro_field =0.0

    ! Initialize variables
    tlatqi = 0. ! tendency for Marat's ice_fall routine (not used here)

    fluxbmk = 0. ! surface/top-of-model microphysical fluxes
    fluxtmk = 0.
    mkwle = 0. ! tendencies/fluxes
    mkwsb = 0.
    mkadv = 0.
    mkdiff = 0.

    ! No further allocation necessary
    isallocatedMICRO = .true.
  end if

  ! initialize fields useful for radiation
  reffc = 25.
  reffi = 25.
  reffs = 25.

  flag_micro3Dout = 0
  flag_micro2DZout = 0

  flag_number(:) = 0
  flag_number(incl) = 1
  flag_number(inr) = 1
  flag_number(inci) = 1

  if(douse_reffi) then
      ! Since the RRTMG lookup table for ice optical properties is based on generalized
      !   effective size (Dge), we will compute Dge using the microphysical properties
      !   of ice and pass this along to the radiation in the reffi and reffs arrays.
      reff_ice_holds_Dge = .true.
  end if

end subroutine micro_setparm


!----------------------------------------------------------------------
!!! Initialize microphysics:
!
! this one is guaranteed to be called by SAM at the beginning of each run, initial or restart:

subroutine micro_init()

  implicit none

  integer :: k, status
  real :: tmp_pgam, tmp_lambda
  real :: qc0(nzm)

  ! This tells the radiation whether or not to use the cloud effective sizes
  !   computed by the microphysics.  This is true by default.
  compute_reffc = douse_reffc
  compute_reffi = douse_reffi
!bloss  compute_reffl = .false.

  flag_micro3Dout = 1
  flag_micro3Dout(iqv_old) = 0
  flag_micro3Dout(itabs_old) = 0
  flag_micro2DZout = 1
  flag_micro2DZout(iqv_old) = 0
  flag_micro2DZout(itabs_old) = 0
  flag_micro2DZout(ilambda_ice) = 0


  if(nrestart.eq.0) then
     micro_field = 0.
     reffc = 25.
     reffi = 25.
     !bloss/duplicates micro_field(:,:,:,iqcl)
     !bloss     cloudliq =0.0

     ! compute initial profiles of liquid water - M.K.
     call satadj_liquid(nzm,tabs0,q0,qc0,pres*100.)

     do k=1,nzm
        micro_field(:,:,k,iqv) = q0(k) + qc0(k) ! total water (vapor + cloud liquid)
        micro_field(:,:,k,iqv_old) = q0(k)
        tabs(:,:,k) = tabs0(k)
        micro_field(:,:,k,itabs_old) = tabs0(k)

        !bloss: approx initialization of effective radius based on 
        !  Hugh's formula.  Here, I'm taking the ratio of the gamma functions
        !  to be about two when they're all pulled inside the cube root.  
        !  Not perfect, but should be a reasonable approximation, I think.
        if (qc0(k).gt.0.) then
          micro_field(:,:,k,iqcl) = qc0(k) ! cloud liquid mass mixing ratio

          if(log_predictNc) then
            micro_field(:,:,k,incl) = Nc0*1.e6 ! number mixing ratio
          end if

          if(dofix_pgam) then
            tmp_pgam = pgam_fixed
          else
            tmp_pgam=0.0005714*(Nc0*RHO(K))+0.2714
            tmp_pgam = MAX(2.,MIN(10.,1./(tmp_pgam**2)-1.))
          end if
          
          tmp_lambda = ( (3.14159*1000./6.)*1.e6*Nc0 &
               *(tmp_pgam+3.)*(tmp_pgam+2.)*(tmp_pgam+1.) / qc0(k) )**(1./3.)
          tmp_lambda = MAX((tmp_pgam+1.)/60.e-6,MIN((tmp_pgam+1.)/1.e-6, &
               tmp_lambda))

          reffc(:,:,k) = 1.e6 *(tmp_pgam+3.)/tmp_lambda/2.
          if(masterproc) write(*,*) 'Experimental reffc initialization: ', reffc(1,1,k)
        end if
     end do

     !bloss/duplicates qcl(:,:,:), qci(:,:,:), so not needed
!bloss      qn = 0.
     fluxbmk = 0.
     fluxtmk = 0.

! your initialization calls are here. What is known at this point are the following
! vertical profiles:

! temperature tabs0(1:nzm),
! air density rho(1:nzm),
! pressure pres(1:nzm),
! water vapor (possibly supersaturated, so handle with caution) q0(1:nzm).
! Height of scalar levels is given by z(1:nzm),
! height of leyer interfaces zi(1:nz).
! Thickness of each layer k is computed as dz*adz(k)
! So, the mass of each layer k per unit area is dz*adz(k)*rho(k)

! All the arrays above are available here through the module vars (vars.f90).

! Your additional initialization calls are placed below.
! Remember that all your new files that contain actual microphysics subroutines
! and functions should be added only to the microphysics directory.

!  call ... ! your calls

  end if !nrestart .eq. 0


  mkwle = 0.
  mkwsb = 0.
  mkadv = 0.
  mkdiff = 0.

  !qpsrc = 0.
  !qpevp = 0.


  call p3_init(lookup_file_dir,nCat,'WRF',status)

  if(status.ne.0) then
    write(*,*) 'Initialization of P3 microphysics not successful.  Stopping...'
    call task_abort()
  end if

  if(docloud) call micro_diagnose()   ! leave this call here


  ! set up names, units and scales for these microphysical quantities
  mkname(iqv) = 'QTO'
  mklongname(iqv) = 'TOTAL WATER (VAPOR + CLOUD LIQUID)'
  mkunits(iqv) = 'g/kg'
  mkoutputscale(iqv) = 1.e3

  mkname(incl) = 'NC'
  mklongname(incl) = 'CLOUD WATER NUMBER CONCENTRATION'
  mkunits(incl) = '#/cm3' ! bloss: output as number mixing ratio
  mkoutputscale(incl) = 1.e-6

  mkname(iqr) = 'QR'
  mklongname(iqr) = 'RAIN'
  mkunits(iqr) = 'g/kg'
  mkoutputscale(iqr) = 1.e3

  mkname(inr) = 'NR'
  mklongname(inr) = 'RAIN NUMBER CONCENTRATION'
  mkunits(inr) = '#/cm3' ! bloss: output as number mixing ratio
  mkoutputscale(inr) = 1.e-6

  mkname(iqit) = 'QI'
  mklongname(iqit) = 'CLOUD ICE mass mixing ratio (Deposition + Rime Ice)'
  mkunits(iqit) = 'g/kg'
  mkoutputscale(iqit) = 1.e3

  mkname(inci) = 'NI'
  mklongname(inci) = 'CLOUD ICE NUMBER CONCENTRATION'
  mkunits(inci) = '#/cm3' ! bloss: output as number mixing ratio
  mkoutputscale(inci) = 1.e-6

  mkname(iqir) = 'QIR'
  mklongname(iqir) = 'rime ice mass mixing ratio '
  mkunits(iqir) = 'g/kg'
  mkoutputscale(iqir) = 1.e3

  mkname(iqib) = 'QIB'
  mklongname(iqib) = 'rime ice volume mixing ratio'
  mkunits(iqib) = 'cm3/kg' !bloss: Should be volume, not 1/volume
  mkoutputscale(iqib) = 1.e6 !bloss: Scale by 1e6 to give cm3/kg, not m3/kg.

  mkname(iqcl) = 'QC'
  mklongname(iqcl) = 'liquid cloud water mass mixing ratio (diagnosed)'
  mkunits(iqcl) = 'g/kg'
  mkoutputscale(iqcl) = 1.e3

  mkname(itabs_old) = 'TP'
  mklongname(iqv_old) = 'last time step absolute temperature (used by P3 to compute Sw/Si)'
  mkunits(itabs_old) = 'K'
  mkoutputscale(itabs_old) = 1.

  mkname(iqv_old) = 'QVP'
  mklongname(iqv_old) = 'last time step water vapor MMR (used by P3 to compute Sw/Si)'
  mkunits(iqv_old) = 'g/kg'
  mkoutputscale(iqv_old) = 1.e3

  mkname(ilambda_ice) = 'ILA'
  mklongname(ilambda_ice) = 'lambda_ice: slope parameter of ice size distribution'
  mkunits(ilambda_ice) = '1/um'
  mkoutputscale(ilambda_ice) = 1e-6

end subroutine micro_init

!----------------------------------------------------------------------
!!! fill-in the surface and top boundary fluxes here:
!
! Obviously, for liquid/ice water variables those fluxes are zero. They are not zero
! only for water vapor variable and, possibly, for CCN and IN if you have those.

subroutine micro_flux()

  implicit none

  fluxbmk(:,:,:) = 0. ! initialize all fluxes at surface to zero
  fluxtmk(:,:,:) = 0. ! initialize all fluxes at top of domain to zero
  fluxbmk(:,:,index_water_vapor) = fluxbq(:,:)
  fluxtmk(:,:,index_water_vapor) = fluxtq(:,:)

end subroutine micro_flux

!----------------------------------------------------------------------
!!! compute local microphysics processes (beyond advection and SGS diffusion):
!
!  This is the place where the condensation/sublimation, accretion, coagulation, freezing,
!  melting, etc., are done, that is  all the microphysics processes except for the spatial
!  transport and mixing.

! IMPORTANT: For consistancy, you need to use thermodynamic constants like specific heat,
! specific heat of condensation, gas constant, etc, the same as in file params.f90.
! Also, you should assume that the conservative thermodynamic variable during these
! proceses is the liquid/ice water static energy: t = tabs + gz - Lc (qc+qr) - Ls (qi+qs+qg)
! It should not change during all of your point microphysical processes!

subroutine micro_proc()

   implicit none
   real :: tmpqc(nx, nzm), tmpnc(nx, nzm), tmpqr(nx, nzm), tmpnr(nx, nzm), &
           tmpqv(nx, nzm), tmptabs(nx, nzm), tmpth(nx, nzm), tmpw(nx, nzm), &
           tmppres(nx, nzm), tmprho(nx, nzm), tmpdz(nx, nzm), ssat(nx, nzm), &
           diag_zdbz(nx, nzm), diag_effc(nx, nzm), diag_effi(nx, nzm), diag_vmi(nx, nzm), &
           diag_di(nx, nzm), diag_rhopo(nx, nzm),tmpqv_old(nx,nzm),tmpth_old(nx,nzm), &
           tmp_cl_pgam(nx, nzm), tmp_cl_lambda(nx,nzm), tmp_ice_lambda(nx,nzm,nCat)

   !bloss: P3 allows multiple ice species.  We are restricting to one for now, but to make
   !  the interface consistent, declare these as three-dimensional arrays for now.
   real :: tmpqit(nx, nzm, nCat), tmpqir(nx, nzm, nCat), tmpnit(nx, nzm, nCat), tmpbir(nx, nzm, nCat)

   real :: tmp_micro_proc_rates(nx,nzm,nmicro_proc) !!!BG what is passed from p3 main

   real :: diag_2d(nx,n_diag_2d), diag_3d(nx,nzm,n_diag_3d)

   real ::  tmptabs_3d(nx,ny,nzm) !,tmpw_3d(nx,ny,nzm) BG

   real :: pcprt_liq(nx), pcprt_sol(nx)
   integer :: i, j, k, its, ite, kts, kte,itstep,m, n, kk

   ! extra inputs for p3_main in WRFV4
   logical :: debug_on      = .false. !switch for internal debug checking
   real    :: clbfact_dep   = 1.0     !calibration factor for deposition
   real    :: clbfact_sub   = 1.0     !calibration factor for sublimation

   ! Note that SCPF (partial cloud fraction) scheme is not enabled at present
   logical                 :: scpf_on = .false.     ! switch for activation of SCPF scheme
   real                    :: scpf_pfrac            ! precipitation fraction factor (SCPF)
   real                    :: scpf_resfact          ! model resolution factor (SCPF)
   real, dimension(nx,nzm) :: SCF_out        ! cloud fraction computed by SCPF

   real, dimension(nx) :: prt_drzl, prt_rain, prt_crys, prt_snow, prt_grpl, prt_pell, prt_hail, prt_sndp
   real, dimension(nx,nzm,n_qiType  ):: qi_type      ! mass mixing ratio, diag ice type    kg kg-1


   !bloss: get microphysical and sedimentation tendencies of all fields.
   real, dimension(nx, nzm, nmicro_fields) :: tmp_mtend, tmp_stend
   integer :: icol, lchnk
   real :: tmpc, tmpr, tmpi

   logical :: do_accumulate_micro_proc_rates !BG added a la M2005/PB

   call t_startf ('micro_proc')

   if(mod(nstep-1,nstatis).eq.0.and.icycle.eq.1) then
     do j=1,ny
       do i=1,nx
         precsfc(i,j)=0.
       end do
     end do
     do k=1,nzm
       precflux(k) = 0.
     end do
   end if

   if(doprecip) total_water_prec = total_water_prec + total_water()

!BG added a la M2005/bloss
if(dostatis) then ! initialize arrays for statistics
   micro_proc_rates(:,:,:,:) = 0.
   !micro_proc_rates_stat(:,:)= 0.
   do_accumulate_micro_proc_rates = dostatis.AND.do_output_micro_process_rates !do_output_micro_process_rates

   trtau(:,:) = 0.
   qxoeffr(:,:) = 0.
   qpfall(:)=0.
   tlat(:) = 0.
end if
!BG end

   mtend(:,:) = 0. !bloss
   stend(:,:) = 0. !bloss
   mksed(:,:) = 0.
   do j = 1, ny

     tmpqc(1:nx, 1:nzm ) = micro_field(1:nx, j, 1:nzm, iqcl)
     tmpnc(1:nx, 1:nzm) = micro_field(1:nx, j, 1:nzm, incl)
     tmpqr(1:nx, 1:nzm) = micro_field(1:nx, j, 1:nzm, iqr)
     tmpnr(1:nx, 1:nzm) = micro_field(1:nx, j, 1:nzm, inr)
     tmpqv(1:nx, 1:nzm) = micro_field(1:nx, j, 1:nzm, iqv) - tmpqc(1:nx, 1:nzm)  !water vapor
     tmpqv_old(1:nx, 1:nzm) = micro_field(1:nx,j,1:nzm,iqv_old)  !water vapor from last time step

     !bloss: These temporary arrays have an extra argument in case we permit multiple ice
     !  species later
     tmpqit(1:nx, 1:nzm, 1) = micro_field(1:nx, j, 1:nzm, iqit)! total ice mass mixing ratio
     tmpqir(1:nx, 1:nzm, 1) = micro_field(1:nx, j, 1:nzm, iqir)! rime ice mass mixing ratio
     tmpnit(1:nx, 1:nzm, 1) = micro_field(1:nx, j, 1:nzm, inci)! total ice number mixing ratio
     tmpbir(1:nx, 1:nzm, 1) = micro_field(1:nx, j, 1:nzm, iqib)! rime ice volume mixing ratio

    tmp_micro_proc_rates(1:nx,1:nzm,1:nmicro_proc) = micro_proc_rates(1:nx,j,1:nzm,1:nmicro_proc) !BG initialize

    !bloss: add look over k so that memory access is more efficient
    do k = 1,nzm
      do i = 1, nx
        tmptabs(i, k) = t(i,j,k)  &  ! liquid water-ice static energy over Cp
             - gamaz(k) &                             ! potential energy
             + fac_cond * (tmpqr(i,k)+tmpqc(i,k)) &  ! qt: liquid latent energy
             + fac_sub  * (tmpqit(i, k, 1)) ! ice latent energy

        !convert to potential temperature:
        tmpth(i, k) = tmptabs(i,k)*(1.e+5/(pres(k)*100.0))**0.286 !note pres is mb
        tmpth_old(i, k) = micro_field(i,j,k,itabs_old)*(1.e+5/(pres(k)*100.0))**0.286 !note pres is mb

        !copy from the Morrison microphysics
        tmpw(i,k) = ((zi(k+1)-z(k))*w(i,j,k)+ &
             (z(k)-zi(k))*w(i,j,k+1))/(zi(k+1)-zi(k))

        tmpw_3d(i,j,k) = ((zi(k+1)-z(k))*w(i,j,k)+ &
             (z(k)-zi(k))*w(i,j,k+1))/(zi(k+1)-zi(k))

        tmppres(i,k) = 100.0*pres(k)
        tmprho(i,k) = rho(k)

        tmpdz(i,k) = adz(k)*dz
      end do
    end do


     its = 1
     ite = nx
     kts = 1
     kte = nzm

     ! P3 arguments added in WRFV4 !!! TODO: Declare these above !!!
     clbfact_dep = 1. ! calibration factor for deposition
     clbfact_sub = 1. ! calibration factor for sublimation
     debug_on = .false.
     scpf_on = .false. ! enable partial cloudiness for microphysics
     scpf_pfrac = 0.
     scpf_resfact = 0.
     SCF_out(:,:) = 0.
     diag_zdbz(:,:) = 0.
     diag_effc(:,:) = 0.
     diag_effi(:,:) = 0.
     diag_vmi(:,:) = 0.
     diag_di(:,:) = 0.
     diag_rhopo(:,:) = 0.
     diag_2d(:,:) = 0.
     diag_3d(:,:,:) = 0.

     ! note: code for prediction of ssat not currently avaiable, set 2D array to 0
      ssat=0.

      !bloss: These arrays now hold entries for all microphysical fields.
      tmp_stend(:,:,:) = 0.0
      tmp_mtend(:,:,:) = 0.0

      itstep = nstep
      do i=1,nx
        kk = 1
        call p3_main( qc=tmpqc(i:i,kk:), nc=tmpnc(i:i,kk:), qr=tmpqr(i:i,kk:), &
              nr=tmpnr(i:i,kk:), th_old=tmpth_old(i:i,kk:), th=tmpth(i:i,kk:), &
              qv_old=tmpqv_old(i:i,kk:), qv=tmpqv(i:i,kk:), &
              dt=dtn, qitot=tmpqit(i:i,kk:,:), qirim=tmpqir(i:i,kk:,:), nitot=tmpnit(i:i,kk:,:), &
              birim=tmpbir(i:i,kk:,:), ssat=ssat(i:i,kk:), uzpl=tmpw(i:i,kk:), &
              pres=tmppres(i:i,kk:), dzq=tmpdz(i:i,kk:), it=itstep, &
              prt_liq=pcprt_liq(i:i), prt_sol=pcprt_sol(i:i), &
              its=1, ite=1, kts=1, kte=nzm-kk+1, nCat=nCat,            &
              diag_ze=diag_zdbz(i:i,kk:), diag_effc=diag_effc(i:i,kk:), diag_effi=diag_effi(i:i,kk:), &
              diag_vmi=diag_vmi(i:i,kk:), diag_di=diag_di(i:i,kk:), diag_rhoi=diag_rhopo(i:i,kk:), &
              n_diag_2d=n_diag_2d, diag_2d=diag_2d(i:i,:), &
              n_diag_3d=n_diag_3d, diag_3d=diag_3d(i:i,kk:,:), &
              log_predictNc=log_predictNc, typeDiags_ON=typeDiags_ON, model=trim(model), &
              clbfact_dep=clbfact_dep, clbfact_sub=clbfact_sub, debug_on=debug_on, &
              scpf_on=scpf_on, scpf_pfrac=scpf_pfrac, scpf_resfact=scpf_resfact, &
              SCF_out=SCF_out(i:i,kk:), &
              prt_drzl=prt_drzl(i:i), prt_rain=prt_rain(i:i), &
              prt_crys=prt_crys(i:i), prt_snow=prt_snow(i:i), &
              prt_grpl=prt_grpl(i:i), prt_pell=prt_pell(i:i), &
              prt_hail=prt_hail(i:i), prt_sndp=prt_sndp(i:i), &
              qi_type=qi_type(i:i,kk:,:), rho=tmprho(i:i,kk:), &
              pgam_out=tmp_cl_pgam(i:i,kk:), lamc_out=tmp_cl_lambda(i:i,kk:), lami_out=tmp_ice_lambda(i:i,kk:,:), &
              mtend=tmp_mtend(i:i,kk:,:), stend=tmp_stend(i:i,kk:,:), &
              micro_proc_rates=tmp_micro_proc_rates(i:i,kk:,:), &
              do_accumulate_micro_proc_rates=do_accumulate_micro_proc_rates) !BG added
      end do ! i


    ! update microphysical quantities in this grid column
!bloss/duplicates micro_field(:,:,:,iqcl)
!bloss     cloudliq(1:nx,j,1:nzm) = tmpqc(1:nx,1:nzm)
     micro_field(1:nx,j,1:nzm,iqcl) = tmpqc(1:nx,1:nzm)
     micro_field(1:nx, j, 1:nzm, incl) = tmpnc(1:nx, 1:nzm)
     micro_field(1:nx, j, 1:nzm, iqr) = tmpqr(1:nx, 1:nzm)
     micro_field(1:nx, j, 1:nzm, inr) = tmpnr(1:nx, 1:nzm)
     micro_field(1:nx, j, 1:nzm, iqv_old) =tmpqv_old(1:nx,  1:nzm)  !water vapor
     
     ! iqv holds both water vapor and cloud liquid water
     micro_field(1:nx, j, 1:nzm, iqv) = tmpqv(1:nx, 1:nzm) + tmpqc(1:nx,1:nzm)  ! water vapor + liquid cloud water


     micro_field(1:nx, j, 1:nzm, iqit) = tmpqit(1:nx, 1:nzm, 1)  ! total ice mass mixing ratio
     micro_field(1:nx, j, 1:nzm, iqir) = tmpqir(1:nx, 1:nzm, 1)  ! rime ice mass mixing ratio
     micro_field(1:nx, j, 1:nzm, inci) = tmpnit(1:nx, 1:nzm, 1)  ! total ice number mixing ratio
     micro_field(1:nx, j, 1:nzm, iqib) = tmpbir(1:nx, 1:nzm, 1)  ! rime ice volumne mixing ratio

     reffi(1:nx,j,1:nzm) = diag_effi(1:nx,1:nzm)*1.e6 !from meter to micor meter
     reffc(1:nx,j,1:nzm) = diag_effc(1:nx,1:nzm)*1.e6 !from meter to micro meter

     micro_proc_rates(1:nx,j,1:nzm,1:nmicro_proc) = tmp_micro_proc_rates(1:nx,1:nzm,1:nmicro_proc) !BG we are in a j loop!

     !bloss: If using RRTMG, use option to convert effective radius to generalized effective size
     ! Here is the code from create_p3_LookupTable_1.f90
     !bloss ! hm 4/9/09, calculate effective size following Fu (1996)
     !bloss !            eff(i_Qnorm,i_Fr) = sum1/(1.1547*916.7*sum2) <-- This is what we want
     !bloss ! hm, calculate for eff rad for twp ice
     !bloss              eff(i_Qnorm,i_Fr) = 3.*sum1/(4.*sum2*916.7) <-- This is diag_effi and now reffi
     if(reff_ice_holds_Dge) then
       !write(*,*) 'max(reffi) before = ', MAXVAL(reffi(:,:,:))
       reffi(1:nx,j,1:nzm) = reffi(1:nx,j,1:nzm)*4./3./1.1547
       !write(*,*) 'max(reffi) after  = ', MAXVAL(reffi(:,:,:))
     end if


!bloss: We are not using the CAM cloud optics with P3 right now, so we are ignoring these variables.
!bloss     CloudLiquidMassMixingRatio(1:nx,j,1:nzm) = tmpqc(1:nx, 1:nzm)
!bloss     CloudLiquidGammaExponent(1:nx,j,1:nzm) = tmp_cl_pgam(1:nx, 1:nzm)
!bloss     CloudLiquidLambda(1:nx,j,1:nzm) = tmp_cl_lambda(1:nx,1:nzm)
     micro_field(1:nx,j,1:nzm,ilambda_ice) = tmp_ice_lambda(1:nx,1:nzm,1)

     !update surface precipitation
     precinst(1:nx, j) = precinst(1:nx,j) + pcprt_liq(1:nx)*1000.0 &
                                +  pcprt_sol(1:nx)*1000.0 
     precsfc(1:nx, j) = precsfc(1:nx,j) + pcprt_liq(1:nx)*1000.0*dtn/dz & !m/s-->mm/s->mm/dz
                                +  pcprt_sol(1:nx)*1000.0*dtn/dz !m/s-->mm/s->mm/dz
     prec_xy(1:nx, j) = prec_xy(1:nx,j) + pcprt_liq(1:nx)*1000.0*dtn/dz & !m/s-->mm/s->mm/dz
                                +  pcprt_sol(1:nx)*1000.0*dtn/dz !m/s-->mm/s->mm/dz

    ! precssfc(1:nx, j) = precssfc(1:nx, j) + pcprt_sol(1:nx)*1000.0*dtn/dz !m/s-->mm/s->mm/dz
     do k = 1,nzm
       do i = 1, nx
         !=====================================================
         ! update liquid-ice static energy due to precipitation
         t(i,j,k) = t(i,j,k) &
              - dtn*fac_cond*(tmp_stend(i,k,iqv) + tmp_stend(i,k,iqr)) & !bloss(stend)
              - dtn*fac_sub*(tmp_stend(i,k,iqit)) !bloss(stend)
         !=====================================================

         micro_field(i,j,k,itabs_old) =  t(i,j,k) &
              -gamaz(k) + fac_cond * (tmpqc(i,k)+tmpqr(i,k)) &
              + fac_sub *(tmpqit(i,k, 1) )

         stend(k,:) = stend(k,:) + tmp_stend(i,k,:) !bloss(stend)
         mtend(k,:) = mtend(k,:) + tmp_mtend(i,k,:) !bloss(mtend)
       end do !i=1,nx
     end do !k=1,nzm

  !   do i = 1, nx
  !     total_water_prec = total_water_prec +(pcprt_liq(i) + pcprt_sol(i))*1000.0*dtn !mm
  !   end do !i=1,nx

     if(dostatis) then

       !bloss: gather these for the statistics HLLAT and QPFALL
       do k = 1,nzm
         do i = 1,nx
           tlat(k) = tlat(k) &
                - dtn*fac_cond*( tmp_stend(i,k,iqv) + tmp_stend(i,k,iqr) ) &
                - dtn*fac_sub*tmp_stend(i,k,iqit)
           qpfall(k) = qpfall(k) + dtn*tmp_stend(i,k,iqv) &
                + dtn*tmp_stend(i,k,iqr) + dtn*tmp_stend(i,k,iqit)
         end do
       end do

       !BG start here tau and radius
       ! approximate optical depth = 1.5e-3*lwp/effrad !bloss(2018-07): Should be 3/2, not 1.8
       !  integrated up to level at which output
       do i=1,nx
         !bloss: move i to outer loop, so that we are accumulating optical depth in each column.
         !   move initialization of trtau to top of micro_proc subroutine so that it will be
         !   properly accumulated in 3D simulations.
         tmpc = 0.
         tmpi = 0.
         do k = 1,nzm    !BG diag_effc/effi is in meters!
          tmpc = tmpc + 1.5e-3*rho(k)*dz*adz(k)*tmpqc(i,k)/(1.e-20+ diag_effc(i,k))!effc1d(k))
          !no reff rain passed from p3
          !tmpr = tmpr + 1.5e-3*rho(k)*dz*adz(k)*tmpqr(i,k)/(1.e-20+ 1.e-6*diag_effr(i,k ))!effr1d(k))
          !bloss/qt: put cloud liquid optical depth in trtau(:,iqv)
          qxoeffr(k,iqcl) = qxoeffr(k,iqcl) + tmpqc(i,k)/(1.e-20+ diag_effc(i,k))
          trtau(k,iqcl) = trtau(k,iqcl) + tmpc
          !bg if(doprecip) trtau(k,iqr) = trtau(k,iqr) + tmpr

          !if(doicemicro) then
             if((tmpqit(i,k,1).gt.1.e-6).AND.(diag_effi(i,k).gt.1.e-6)) then
               tmpi = tmpi + 1.5e-3*rho(k)*dz*adz(k)*tmpqit(i,k,1)/(1.e-20+ diag_effi(i,k))
               qxoeffr(k,iqit) = qxoeffr(k,iqit) + tmpqit(i,k,1)/(1.e-20+ diag_effi(i,k))
             end if
             trtau(k,iqit) = trtau(k,iqit) + tmpi
          !end if
           end do!k
         end do!i
       !BG end
     end if !dostatis

   end do ! j =1,ny

   ! back sedimentation flux out from sedimentation tendencies
   do n = 1,nmicro_fields
     if(n.ne.iqcl) then
       tmpc = 0.
       do k = 1,nzm
         m = nz-k
         tmpc = tmpc + stend(m,n)*rho(m)*dz*adz(m)
         mksed(m,n) = tmpc
       end do
     end if
   end do
   precflux(1:nzm) = precflux(1:nzm) &
        - dtn/dz*( mksed(:,iqv) + mksed(:,iqr) + mksed(:,iqit) )

   if(doprecip) total_water_prec = total_water_prec - total_water()

   if (docloud)  call micro_diagnose()   ! leave this line here

   call t_stopf ('micro_proc')


end subroutine micro_proc

!----------------------------------------------------------------------
!!! Diagnose arrays nessesary for dynamical core and radiation:
!
! This is the pace where the microphysics field that SAM actually cares about
! are diagnosed. You need to compute all the arrays on the left-hand-side in the loop below
! for SAM dynamical core to see your microphysics (that is to see the cloud and precipitation).

subroutine micro_diagnose()

   real omn, omp
   integer i,j,k

   do k=1,nzm
    do j=1,ny
     do i=1,nx
       qv(i,j,k) = micro_field(i,j,k,iqv)- micro_field(i,j,k,iqcl)
       qcl(i,j,k) = micro_field(i,j,k,iqcl)
       qci(i,j,k) = micro_field(i,j,k,iqit) 
       qpl(i,j,k) = micro_field(i,j,k,iqr)
       qpi(i,j,k) = 0.
     end do
    end do
   end do

end subroutine micro_diagnose

!----------------------------------------------------------------------
!!! functions to compute terminal velocity for precipitating variables:
!
! you need to supply your own functions functions to compute terminal velocity
! for all of your precipitating prognostic variables. Note that all functions should
! compute vertical velocity given two microphysics parameters var1, var2,
! and temperature, and water vapor (single values, not arrays). Also, for
! bin-microphysics implementation, there is a fifth variable with the type of
! integer that can be used for bin index. Var1 and var2
! are some microphysics variables like water content and concentration.
! IMPORTANT: Don't change the number of arguments or their meaning!

!real function term_vel_qr(i,j,k,ind)
! .......
!end function term_vel_qr

!real function term_vel_Nr(i,j,k,ind)
! .......
!end function term_vel_Nr

!real function term_vel_qs(i,j,k,ind)
! .......
!end function term_vel_qs

! etc.

!----------------------------------------------------------------------
!!! compute sedimentation
!
!  The purpose of this subroutine is to prepare variables needed to call
! the precip_all() for each of the falling hydrometeor varibles

subroutine micro_precip_fall()

! before calling precip_fall() for each of falling prognostic variables,
! you need to set hydro_type and omega(:,:,:) variables.
! hydro_type can have four values:
! 0 - variable is liquid water mixing ratio
! 1 - hydrometeor is ice mixing ratio
! 2 - hydrometeor is mixture-of-liquid-and-ice mixing ratio. (As in original SAM microphysics).
! 3 - variable is not mixing ratio, but, for example, rain drop concentration
! OMEGA(:,:,:) is used only for hydro_type=2, and is the fraction of liquid phase (0-1).
! for our hypothetical case, there is no mixed hydrometeor, so omega is not actually used.
! In default SAM microphysics, omega is a mass partition between liquid and ice phases.

  integer hydro_type
  real omega(nx,ny,nzm)
  integer ind ! variable that is reserved for bin-microphysics use (bin index).

  integer i,j,k

  return ! do not need this routine -- sedimentation done in p3_main.
! Initialize arrays that accumulate surface precipitation flux

!# if(mod(nstep-1,nstatis).eq.0.and.icycle.eq.1) then
!#   do j=1,ny
!#    do i=1,nx
!#     precsfc(i,j)=0.
!#    end do
!#   end do
!#   do k=1,nzm
!#    precflux(k) = 0.
!#   end do
!# end if
!#
!# do k = 1,nzm ! Initialize arrays which hold precipitation fluxes for stats.
!#    qpfall(k)=0.
!#    tlat(k) = 0.
!# end do
!#
!#! Compute sedimentation of falling variables:
!#
!# hydro_type=0
!# call precip_fall(qr, term_vel_qr, hydro_type, omega, ind)
!# hydro_type=3
!# call precip_fall(Nr, term_vel_Nr, hydro_type, omega, ind)
!# hydro_type=1
!# call precip_fall(qs, term_vel_qs, hydro_type, omega, ind)
!# hydro_type=3
!# call precip_fall(Ns, term_vel_Ns, hydro_type, omega, ind)
!# hydro_type=1
!# call precip_fall(qg, term_vel_qg, hydro_type, omega, ind)
!# hydro_type=3
!# call precip_fall(Ng, term_vel_Ng, hydro_type, omega, ind)

end subroutine micro_precip_fall

!----------------------------------------------------------------------
!!! Initialize the list of microphysics statistics that will be outputted
!!  to *.stat statistics file

subroutine micro_hbuf_init(namelist,deflist,unitlist,status,average_type,count,trcount)


  character(*) namelist(*), deflist(*), unitlist(*)
  integer status(*),average_type(*),count,trcount
  integer ntr
  integer n, ii, jj, ncond

  character*8 name
  character*80 longname
  character*10 units
  trcount =0

  name = 'QTFLUX'
  longname = 'Total (resolved + subgrid) total water (vapor+cloud) flux'
  units = 'W/m2'
  call add_to_namelist(count,trcount,name,longname,units,0)

  do n = 1,nmicro_fields
    if(flag_advect(n).eq.1) then
      ! add mean value of microphysical field to statistics
      !   if they are advected variables (does not include qv_old, tabs_old)
      name = trim(mkname(n))
      longname = trim(mklongname(n))
      units = trim(mkunits(n))
      call add_to_namelist(count,trcount,name,longname,units,0)
      if(n.eq.iqv) then
        ! add variance of ONLY total water (vapor + cloud liq) field to statistics
        !   cloud water variance and cloud ice variance
        !   already output in statistics.f90
        name = trim(mkname(n))//'2'
        longname = 'Variance of '//trim(mklongname(n))
        units = '('//trim(mkunits(n))//')^2'
        call add_to_namelist(count,trcount,name,longname,units,0)
      end if

      if(n.ne.iqcl) then
        ! output budget tendencies for all advected species (everything but cloud liquid).
        name = trim(mkname(n))//'ADV'
        longname = 'Advective Tendency of '//trim(mklongname(n))
        units = trim(mkunits(n))//'/day'
        call add_to_namelist(count,trcount,name,longname,units,0)

        name = trim(mkname(n))//'DIFF'
        longname = 'Subgrid Mixing Tendency of '//trim(mklongname(n))
        units = trim(mkunits(n))//'/day'
        call add_to_namelist(count,trcount,name,longname,units,0)

        name = trim(mkname(n))//'MPHY'
        longname = 'Microphysical Tendency of '//trim(mklongname(n))
        units = trim(mkunits(n))//'/day'
        call add_to_namelist(count,trcount,name,longname,units,0)

        name = trim(mkname(n))//'SED'
        longname = 'Sedimentation Tendency of '//trim(mklongname(n))
        units = trim(mkunits(n))//'/day'
        call add_to_namelist(count,trcount,name,longname,units,0)

        name = trim(mkname(n))//'SDFLX'
        longname = 'Sedimentation flux of '//trim(mklongname(n))
        units = trim(mkunits(n))//' m/s'
        call add_to_namelist(count,trcount,name,longname,units,0)

        name = trim(mkname(n))//'FLXR'
        longname = 'Flux of '//trim(mklongname(n))//' due to resolved eddies'
        units = trim(mkunits(n))//' m/s'
        call add_to_namelist(count,trcount,name,longname,units,0)

        name = trim(mkname(n))//'FLXS'
        longname = 'Flux of '//trim(mklongname(n))//' due to unresolved (subgrid) eddies'
        units = trim(mkunits(n))//' m/s'
        call add_to_namelist(count,trcount,name,longname,units,0)
      end if

      if ((n.eq.iqcl).OR.(n.eq.iqit)) then
        ! add approximate optical depth of hydrometeor fields
        name = 'TAU'//trim(mkname(n))
        longname = 'Approx optical depth of '//trim(mklongname(n))
        units = '1'
        call add_to_namelist(count,trcount,name,longname,units,0)

        ! add field which can be used to recover mean effective radius.
        name = trim(mkname(n))//'OEFFR'
        longname = 'Mixing ratio of '//trim(mklongname(n)) &
             //' over effective radius, EFFR = ' &
             //trim(mkname(n))//'/'//trim(mkname(n))//'OEFFR'
        units = 'g/kg/microns'
        call add_to_namelist(count,trcount,name,longname,units,0)
      end if

    end if
  end do !n=1, nmicro_fields

  ! output conditional averages of each advected microphysical variable
  do ncond = 1,ncondavg
    do n = 1,nmicro_fields
      if(flag_advect(n).eq.1) then
        call add_to_namelist(count,trcount,trim(mkname(n)) // TRIM(condavgname(ncond)), &
             trim(mklongname(n)) // ' in ' // TRIM(condavglongname(ncond)), &
             trim(mkunits(n)),ncond)
      end if
    end do
  end do

  !BG a la' bloss for process rate - not known whethere mass or number
  if(do_output_micro_process_rates) then
    do n = 1,nmicro_proc
      call add_to_namelist(count,trcount, &
           trim(micro_process_rate_names(n)), &
           trim(micro_process_rate_longnames(n)), &
           '?/kg/day',0)
    end do
  end if
  !BG end----------------------------------------------------

end subroutine micro_hbuf_init

!----------------------------------------------------------------------
!!!! Collect microphysics history statistics (vertical profiles)
!! Note that only the fields declared in micro_hbuf_init() are allowed to
! be collected

subroutine micro_statistics()

  use hbuffer, only: hbuf_put

  real, dimension(nzm) :: tr0, tr2, tr0u,tr0d
  !tr0 = for xy sums of micro fields which are later averaged by factor_xy (number of gridpoints)
  !tr2 = for variance of specific fields
  real, dimension(nx,ny,nzm) :: tr3u, tr3d !BG

  real tmp(2), factor_xy
  real qcz(nzm), qiz(nzm), qrz(nzm), qsz(nzm), qgz(nzm), omg
  integer i,j,k,m,n, ncond


  factor_xy = 1./float(nx*ny)

  do n =1, nmicro_fields
    if(flag_advect(n).eq.1) then
      ! only output statistics for advected variables
      do k=1,nzm
        tmp(1) = dz/rhow(k)
        tmp(2) = tmp(1) / dtn
        tr0(k) = SUM(micro_field(1:nx,1:ny,k,n))
        tr2(k) = SUM(micro_field(1:nx,1:ny,k,n)*micro_field(1:nx,1:ny,k,n))
        mkwsb(k,n) = mkwsb(k,n) * tmp(1) * rhow(k) !subgrid flux
        mkwle(k,n) = mkwle(k,n)*tmp(2)*rhow(k)     !resolved flux
      end do

      !bloss: Convert number concentrations /cm3.
      if(flag_number(n).eq.1) then
        ! remove factor of rho from number concentrations
        tr0(:) = tr0(:)*rho(:)
        tr2(:) = tr2(:)*rho(:)**2
        mkadv(1:nzm,n) = mkadv(1:nzm,n)*rho(:)
        mkdiff(1:nzm,n) = mkdiff(1:nzm,n)*rho(:)
        mtend(1:nzm,n) = mtend(1:nzm,n)*rho(:)
        stend(1:nzm,n) = stend(1:nzm,n)*rho(:)
        mklsadv(1:nzm,n) = mklsadv(1:nzm,n)*rho(:)
      end if

      ! output all microphysical fields
      call hbuf_put(trim(mkname(n)),tr0,mkoutputscale(n)*factor_xy)
      if(n.eq.iqv) then
        ! variance of microphysical field,  only for QTO (qv+qcl)
        call hbuf_put(trim(mkname(n))//'2',tr2,mkoutputscale(n)**2*factor_xy)
      end if

      if(n.ne.iqcl) then
        ! output budget tendencies for all advected species (everything but cloud liquid).
        call hbuf_put(trim(mkname(n))//'ADV',mkadv(1,n),86400.*mkoutputscale(n)*factor_xy)
        call hbuf_put(trim(mkname(n))//'DIFF',mkdiff(1,n),86400.*mkoutputscale(n)*factor_xy)
        call hbuf_put(trim(mkname(n))//'MPHY',mtend(1,n),86400.*mkoutputscale(n)*factor_xy)
        call hbuf_put(trim(mkname(n))//'SED',stend(1,n),86400.*mkoutputscale(n)*factor_xy)
        call hbuf_put(trim(mkname(n))//'SDFLX',mksed(1,n),mkoutputscale(n)*factor_xy)
        call hbuf_put(trim(mkname(n))//'FLXR',mkwle(1,n),mkoutputscale(n)*factor_xy)
        call hbuf_put(trim(mkname(n))//'FLXS',mkwsb(1,n),mkoutputscale(n)*factor_xy)
      end if

      !BG added
      if (n .eq. (iqit .or. iqcl )) then
        ! approx optical depth
        call hbuf_put('TAU'//trim(mkname(n)),trtau(1,n),factor_xy)
        !multiply through by 1e3/1e6 (conversion of qi from kg/kg to g/kg in numerator
        !   and meters to microns in denominator)
        call hbuf_put(trim(mkname(n))//'OEFFR',qxoeffr(1,n),factor_xy*1.e3/1.e6)
      end if
      
      do ncond = 1,ncondavg
        do k = 1,nzm
          tr0(k) = SUM(micro_field(1:nx,1:ny,k,n)*condavg_mask(1:nx,1:ny,k,ncond))
        end do
        call hbuf_put(TRIM(mkname(n)) // TRIM(condavgname(ncond)), &
             tr0,mkoutputscale(n))
      end do

    end if !if(flag_advect(n).eq.1)

  end do !n=1, nmicro_fields

  !bloss: Add output for QTFLUX.  Here, we interpret QT as including all vapor, liquid and ice.
  tr0(:) = 0.
  do n = 1,nmicro_fields
    if (flag_wmass(n).eq.1) then
      tr0(:) = tr0(:) + lcond*(mkwle(1:nzm,n) + mkwsb(1:nzm,n))
    end if
  end do
  call hbuf_put('QTFLUX',tr0,factor_xy)

!!$  ! add separate output for cloud liquid water
!!$  !           and approx cloud liquid optical depth.
!!$  do k = 1,nzm
!!$    tr0(k) = SUM(cloudliq(1:nx,1:ny,k))
!!$  end do
!!$  call hbuf_put('QC',tr0,factor_xy*mkoutputscale(iqv))

  !BG a la' bloss from M2005 for micro process rates output
  if(do_output_micro_process_rates) then
    do n = 1,nmicro_proc
      do k=1,nzm
        ! that's how was before tr0(1:nzm) = micro_proc_rates(1:nzm,n)*86400.*factor_xy
        tr0(k) = SUM(micro_proc_rates(1:nx,1:ny,k,n))*86400*factor_xy !units: per day, field averages
      end do
      call hbuf_put(TRIM(micro_process_rate_names(n)),tr0,1.)
    end do
  end if

end subroutine micro_statistics

!-----------------------------------------
subroutine micro_stat_2Dinit(ResetStorage)
  implicit none
  integer, intent(in) :: ResetStorage

  ! initialize microphysical 2D outputs as necessary

  if(ResetStorage.eq.1) then
    !bloss: If computing storage terms for individual hydrometeors,
    !         store initial profiles for computation of storage terms in budgets
    !bloss mkstor(:,:) = mk0(:,:)
  end if

end subroutine micro_stat_2Dinit

!-----------------------------------------------------------------------
! This one is called when stepout() is called

subroutine micro_print()
  implicit none
  integer :: k

  ! print out min/max values of all microphysical variables
  do k=1,nmicro_fields
     call fminmax_print(trim(mkname(k))//':', &
          micro_field(:,:,:,k),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
  end do

end subroutine micro_print

!-----------------------------------------
subroutine micro_write_fields2D(nfields1)
  implicit none
  integer, intent(inout) :: nfields1

end subroutine micro_write_fields2D

!-----------------------------------------
subroutine micro_write_fields3D(nfields1)
  implicit none
  integer, intent(inout) :: nfields1
!!$  character *80 long_name
!!$  character *8 name
!!$  character *10 units
!!$  integer :: i, j, k, tens, ones
!!$  real(4), dimension(nx,ny,nzm) :: tmp

end subroutine micro_write_fields3D

!-----------------------------------------
subroutine satadj_liquid(nzm,tabs,qt,qc,pres)
  !bloss/qt: Utility routine based on cloud.f90 in 
  !  MICRO_SAM1MOM that was written by Marat Khairoutdinov.
  !  This routine performs a saturation adjustment for
  !  cloud liquid water only using a Newton method.
  !  While 20 iterations are allowed, most often this
  !  routine should exit in five iterations or less.
  !  Only a single calculation of the saturation vapor
  !  pressure is required in subsaturated air.

  implicit none

  integer, intent(in) :: nzm
  real, intent(inout), dimension(nzm) :: tabs ! absolute temperature, K
  real, intent(inout), dimension(nzm) :: qt  ! on input: qt; on output: qv
  real, intent(out), dimension(nzm) :: qc ! cloud liquid water, kg/kg
  real, intent(in), dimension(nzm) :: pres ! pressure, Pa

  real tabs1, dtabs, thresh, esat1, qsat1, fff, dfff
  integer k, niter

  integer, parameter :: maxiter = 20
  real, external :: qsatw, dtqsatw

  !bloss/qt: quick saturation adjustment to compute cloud liquid water content.
  do k = 1,nzm
    if(pres(k).lt.5000.) cycle
    tabs1 = tabs(k) 
    qsat1 = qsatw( tabs(k), 0.01*pres(k) ) ! convert to hPa
    qc(k) = 0. ! no cloud unless qt > qsat
    
    if (qt(k).gt.qsat1) then

      ! if unsaturated, nothing to do (i.e., qv=qt, T=Tl) --> just exit.
      ! if saturated, do saturation adjustment 
      !    (modeled after Marat's cloud.f90).

      ! generate initial guess based on above calculation of qsat
      dtabs = + fac_cond*MAX(0.,qt(k) - qsat1) &
           / ( 1. + fac_cond*dtqsatw( tabs1, 0.01*pres(k) ) ) ! convert to hPa
      tabs1 = tabs1 + dtabs
      niter = 1

      ! convergence threshold: min of 0.01K and latent heating due to
      !    condensation of 1% of saturation mixing ratio.
      thresh = MIN(0.01, 0.01*fac_cond*qsat1)

      ! iterate while temperature increment > thresh and niter < maxiter
      do while((ABS(dtabs).GT.thresh) .AND. (niter.lt.maxiter))

        qsat1 = qsatw(tabs1, 0.01*pres(k) ) ! saturation mixing ratio, convert pressure to hPa

        fff = tabs(k) - tabs1 + fac_cond*MAX(0.,qt(k) - qsat1)
        dfff = 1. + fac_cond*dtqsatw( tabs1, 0.01*pres(k) ) ! convert to hPa
        dtabs = fff/dfff
        tabs1 = tabs1 + dtabs

        niter = niter + 1

      end do

      qc(k) = MAX( 0.,tabs1 - tabs(k) )/fac_cond ! cloud liquid mass mixing ratio
      qt(k) = qt(k) - qc(k) ! This now holds the water vapor mass mixing ratio.
      tabs(k) = tabs1 ! update temperature.
      
      if(niter.gt.maxiter-1) write(*,*) 'Reached iteration limit in satadj_liquid'

    end if ! qt_in > qsat

  end do ! k = 1,nzm

end subroutine satadj_liquid

!-----------------------------------------------------------------------
! Function that computes total water in the domain:
! Don't change this one.

real(8) function total_water()

 use grid, only: nx_gl,ny_gl
  use vars, only : nstep,nprint,adz,dz,rho
  real(8) tmp
  integer i,j,k,m


  total_water = 0.0

  if(mod(nstep,nprint).ne.0) return


  do m=1,nmicro_fields
   if(flag_wmass(m).eq.1) then
    do k=1,nzm
      tmp = 0.
      do j=1,ny
        do i=1,nx
          tmp = tmp + micro_field(i,j,k,m)
        end do
      end do
      total_water = total_water + tmp*adz(k)*dz*rho(k)
    end do
   end if
  end do
  total_water = total_water/dble(nx_gl*ny_gl)


end function total_water

! -------------------------------------------------------------------------------
! If your microphysics allows you to compute drop/ice effective radius,
! insert the code here. If note, leave blank functions for compilation perposes.
logical function micro_provides_reffc()
  micro_provides_reffc = douse_reffc
end function micro_provides_reffc

logical function micro_provides_reffi()
  micro_provides_reffi = douse_reffi
end function micro_provides_reffi

function Get_reffc() ! liquid water
  real, dimension(nx,ny,nzm) :: Get_reffc
  Get_reffc = reffc
end function Get_reffc

function Get_reffi() ! ice
  real, dimension(nx,ny,nzm) :: Get_reffi
  Get_reffi = reffi
end function Get_reffi

function Get_nca() ! aerosol
  real, pointer, dimension(:,:,:) :: Get_nca
  Get_nca = 0.
end function Get_nca

end module microphysics


! You are done! Good luck!
