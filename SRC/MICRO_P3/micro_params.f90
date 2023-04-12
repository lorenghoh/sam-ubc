module micro_params

!bloss: move indices of different microphysical species here, so that they
!  can be accessed in both microphysics.f90 and module_mp_p3.f90.

! indices of water quantities in micro_field, e.g. qv = micro_field(:,:,:,iqv)
integer :: iqv, iqcl, incl, inci, inr, iqr, iqit, iqir, iqib

! these non-advected fields hold the value of qv and tabs from the last time step.
integer :: iqv_old, itabs_old, ilambda_ice

! SETTINGS FOR WARM RAIN MICROPHYSICS (AND DEFAULT VALUES)

! choice of warm rain microphysics scheme
!   = 1 for Seifert-Beheng (2001)
!   = 2 for Beheng (1994)
!   = 3 for Khairoutdinov-Kogan (2000) <-- DEFAULT
integer, public :: iWarmRainScheme = 3 

! control properties of cloud drop size distribution.
logical, public :: log_predictNc = .false. ! if true, predict droplet number
real, public :: Nc0 = 100. ! specified cloud droplet number conc (CDNC, #/cm3)

logical, public :: dofix_pgam = .false. ! option to fix value of exponent in cloud water gamma distn
real, public ::    pgam_fixed = 10.3 ! Geoffroy et al (2010, doi:10.5194/acp-10-4835-2010)

! background aerosol properties for activation when log_predictNc==.true.
! This is mode 1 (nominally, accumulation mode aerosol)
real, public :: aerosol_mode1_radius = 50.e-9 ! aerosol mean size (m)
real, public :: aerosol_mode1_sigmag = 2.0 ! geometric standard deviation
real, public :: aerosol_mode1_number = 300.e6 ! aerosol number mixing ratio (kg-1)

! This is mode 2.  (Default is coarse mode with zero concentration. Could be Aitken as well.)
real, public :: aerosol_mode2_radius = 1.3e-6 ! aerosol mean size (m)
real, public :: aerosol_mode2_sigmag = 2.5 ! geometric standard deviation
real, public :: aerosol_mode2_number = 0.e6 ! aerosol number mixing ratio (kg-1)

! SETTINGS FOR ICE MICROPHYSICS
integer :: nCat =1  ! number of free ice categories

real :: MaxTotalIceNumber = 2.e6 ! max ice concentration in #/m3 (default = 2/cm3).

!BG done a la' Bloss
!more parameters etc to be added, if needed
!(or put here from elswhere) -> set_param.f90 for instance
   logical, public :: typeDiags_ON = .false.! logic variable, for diagnostic hydrometeor/precip rate types
   character(len=10) :: model ="WRF"    ! type of model that is used in P3_MAIN subroutine

   integer :: n_diag_2d = 1       ! the number of 2-D diagnositic variables output in P3_MAIN subroutine
   integer :: n_diag_3d = 1       ! the number of 3-D diagnositic variables output in P3_MAIN subroutine
   !lookup_file_dir = '/global/homes/g/gxlin/MMF_code/wrfv4_p3_release_092217'

   !options for coupling cloud radiative properties to information
   !  from the microphysics
   logical :: douse_reffc = .true., douse_reffi = .true.

   !!!!!!BG do that properly with P3 lookup tables!!!
   real, parameter :: rho_cloud_ice = 500.
   real, parameter :: rho_snow      = 100.

!BG- a la Bloss in M2005 (9Jun2018): Add outputs of microphysical process rates
!   If icemicro==.true., this amount to an extra XXX outputs
!   in the statistics file.  Only 14 for warm cloud microphysics.
   logical, public :: do_output_micro_process_rates = .false.

   integer :: nmicro_proc
   integer, parameter :: nmicro_process_rates = 8 ! out of 43
   !no need for that: integer, parameter :: nmicro_process_rates_warm = 14
   character(len=8), dimension(nmicro_process_rates), parameter, public :: &
        micro_process_rate_names = (/ &
! liquid-phase microphysical process rates:
!  (all Q process rates in kg kg-1 s-1)
!  (all N process rates in # kg-1)
  !    'qrcon    ', & ! rain condensation
  !   'qcacc    ', & ! cloud droplet accretion by rain
  !   'qcaut    ', & ! cloud droplet autoconversion to rain
  !   'ncacc    '/)!, & ! change in cloud droplet number from accretion by rain
   !  'ncautc   ', & ! change in cloud droplet number from autoconversion
   !  'ncslf    ', & ! change in cloud droplet number from self-collection
   !  'nrslf    ', & ! change in rain number from self-collection
   !  'ncnuc    ', & ! change in cloud droplet number from activation of CCN
   !  'qccon    ', & ! cloud droplet condensation
   !  'qcnuc    ', & ! activation of cloud droplets from CCN
   !  'qrevp    ', & ! rain evaporation
   !  'qcevp    ', & ! cloud droplet evaporation
   !  'nrevp    ', & ! change in rain number from evaporation
   !  'ncautr   ', & ! change in rain number from autoconversion of cloud water
! i!ce-phase microphysical process rates:
!  !(all Q process rates in kg kg-1 s-1)
!  !(all N process rates in # kg-1)
!     'qccol     ', & ! collection of cloud water by ice
!     'qwgrth    ', & ! wet growth rate
     'qidep     ', & ! vapor deposition
!     'qrcol     ', & ! collection rain mass by ice
   !  'qinuc     ', & ! deposition/condensation freezing nuc !for T< -15 C and Sice> 5%, not really sure how
   !  'nccol     ', & ! change in cloud droplet number from collection by ice
   !  'nrcol     ', & ! change in rain number from collection by ice
   !  'ninuc     ', & ! change in ice number from deposition/cond-freezing nucleation
     'qisub     ', & ! sublimation of ice
     'qimlt     ', & ! melting of ice
     'qcrfrz    ', & ! freezing of cloud droplets (qinuc + qchetc + qcheti) and rain (qrhti+qrhetc) !BG
     !number rates
     'nimlt     ', & ! melting of ice
     'nisub     ', & ! change in ice number from sublimation
     'nislf     ', & ! change in ice number from collection within a category
     'ncrfrz    '  /) ! number freezing of cloud droplets (qinuc + qchetc + qcheti) and rain (qrhti+qrhetc) !BG
   !  'qchetc    ', & ! contact freezing droplets  !turned off
   !  'qcheti    ', & ! immersion freezing droplets !turned on
   !  'qrhetc    ', & ! contact freezing rain   !turned off
   !  'qrheti    ', & ! immersion freezing rain !turned on
   !  'nchetc    ', & ! contact freezing droplets
   !  'ncheti    ', & ! immersion freezing droplets
   !  'nrhetc    ', & ! contact freezing rain
   !  'nrheti    ', & ! immersion freezing rain
   !  'nrshdr    ', & ! source for rain number from collision of rain/ice above freezing and shedding
   !  'qcshd     ', & ! source for rain mass due to cloud water/ice collision above freezing and shedding or wet growth and shedding
   !  'qcmul     ', & ! change in q, ice multiplication from rime-splitnering of cloud water (not included in the paper)
   !  'qrmul     ', & ! change in q, ice multiplication from rime-splitnering of rain (not included in the paper)
   !  'nimul     ', & ! change in Ni, ice multiplication from rime-splintering (not included in the paper)
   !  'ncshdc    ', & ! source for rain number due to cloud water/ice collision above freezing  and shedding (combined with NRSHD in the paper)
   !  'rhorime_c ', & ! density of rime (from cloud)
   !  'rhorime_r '  /)! density of rime (from rain)
 !only if you have multiple ice categories (i.e. nCat>1)
 !nicol ! change of N due to ice-ice collision between categories
 !qicol ! change of q due to ice-ice collision between categories

   character(len=80), dimension(nmicro_process_rates), parameter, public :: &
        micro_process_rate_longnames = (/ &
! liquid-phase microphysical process rates:
!  (all Q process rates in kg kg-1 s-1)
!  (all N process rates in # kg-1)
!      'qrcon  , rain condensation [kg kg-1 s-1]                                       ', &
!     'qcacc  , cloud droplet accretion by rain                                       ', &
!     'qcaut  , cloud droplet autoconversion to rain                                  ', &
!     'ncacc  , change in cloud droplet number from accretion by rain # kg-1          ' /)!, &
!     'ncautc , change in cloud droplet number from autoconversion                    ', &
!     'ncslf  , change in cloud droplet number from self-collection                   ', &
!     'nrslf  , change in rain number from self-collection                            ', &
!     'ncnuc  , change in cloud droplet number from activation of CCN                 ', &
!     'qccon  , cloud droplet condensation                                            ', &
!     'qcnuc  , activation of cloud droplets from CCN                                 ', &
!     'qrevp  , rain evaporation                                                      ', &
!     'qcevp  , cloud droplet evaporation                                             ', &
!     'nrevp  , change in rain number from evaporation                                ', &
!     'ncautr , change in rain number from autoconversion of cloud water              ', &
!! ice-phase microphysical process rates:
!!  (all Q process rates in kg kg-1 s-1)
!!  (all N process rates in # kg-1)
!     'qccol   , collection of cloud water by ice                                      ', &
!     'qwgrth  , wet growth rate                                                       ', &
     'qidep   , vapor deposition                                                      ', &
!     'qrcol   , collection rain mass by ice                                           ', &
!     'qinuc   , deposition/condensation freezing nuc                                  ', &
!     'nccol   , change in cloud droplet number from collection by ice                 ', &
!     'nrcol   , change in rain number from collection by ice                          ', &
!     'ninuc   , change in ice number from deposition/cond-freezing nucleation         ', &
     'qisub   , sublimation of ice                                                    ', &
     'qimlt   , melting of ice                                                        ', &
     'qcrfrz  , freezing of cloud droplets and rain drops                             ', &
     'nimlt   , melting of ice                                                        ', &
     'nisub   , change in ice number from sublimation                                 ', &
     'nislf   , change in ice number from collection within a category                ', &
     'ncrfrz  , number freezing of cloud droplets and rain drops                      ' /)! &
!     'qchetc  , contact freezing droplets                                             ', &
!     'qcheti  , immersion freezing droplets                                           ', &
!     'qrhetc  , contact freezing rain                                                 ', &
!     'qrheti  , immersion freezing rain                                               ', &
!     'nchetc  , contact freezing droplets                                             ', &
!     'ncheti  , immersion freezing droplets                                           ', &
!     'nrhetc  , contact freezing rain                                                 ', &
!     'nrheti  , immersion freezing rain                                               ', &
!     'nrshdr  , source for rain number from collision of rain/ice above freezing and shedding', &
!     'qcshd   , source for rain mass due to cloud water/ice collision above freezing and shedding or wet growth and shedding', &
!     'qcmul   , change in q, ice multiplication from rime-splitnering of cloud water  ', & ! (not included in the paper)
!     'qrmul   , change in q, ice multiplication from rime-splitnering of rain         ', & ! (not included in the paper)
!     'nimul   , change in Ni, ice multiplication from rime-splintering                ', & !(not included in the paper)
!     'ncshdc  , source for rain number due to cloud water/ice collision above freezing  and shedding', & ! (combined with NRSHD in the paper)
!     'rhorime_c,density of rime (from cloud)                                          ', &
!     'rhorime_r,density of rime (from rain)                                           '/)
 !only if you have multiple ice categories (i.e. nCat>1)
 !nicol ! change of N due to ice-ice collision between categories
 !qicol ! change of q due to ice-ice collision between categories

 end module micro_params
