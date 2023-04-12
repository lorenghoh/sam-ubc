program crm

!       Main module.

use vars
use hbuffer
use microphysics
use sgs
use tracers
use movies, only: init_movies
implicit none

integer k, icyc, nn, nstatsteps
double precision cputime, oldtime, init_time, elapsed_time !bloss wallclocktime
double precision usrtime, systime
!-------------------------------------------------------------------
! determine the rank of the current task and of the neighbour's ranks

call task_init() 
!------------------------------------------------------------------
! print time, version, etc

if(masterproc) call header()	
!------------------------------------------------------------------
! Initialize timing library.  2nd arg 0 means disable, 1 means enable

   call t_setoptionf (1, 0)
   call t_initializef ()

   call t_startf ('total')
   call t_startf ('initialize')
!------------------------------------------------------------------
! Get initial time of job

!------------------------------------------------------------------

call init()     ! initialize some statistics arrays
call setparm()	! set all parameters and constants

!------------------------------------------------------------------
! Initialize or restart from the save-dataset:

if(nrestart.eq.0) then
   day=day0 
   call setgrid() ! initialize vertical grid structure
   call setdata() ! initialize all variables
elseif(nrestart.eq.1) then
   call read_all()
   call setgrid() ! initialize vertical grid structure
   call diagnose()
   call sgs_init()
   call micro_init()  !initialize microphysics
elseif(nrestart.eq.2) then  ! branch run
   call read_all()
   call setgrid() ! initialize vertical grid structure
   call diagnose()
   call setparm() ! overwrite the parameters
   call sgs_init()
   call micro_init()  !initialize microphysics
   nstep = 0
   day0 = day
else
   print *,'Error: confused by value of NRESTART'
   call task_abort() 
endif

call init_movies()
call stat_2Dinit(1) ! argument of 1 means storage terms in stats are reset
call tracers_init() ! initialize tracers
call setforcing()
if(masterproc) call printout()
!------------------------------------------------------------------
!  Initialize statistics buffer:

call hbuf_init()
	
!------------------------------------------------------------------
total_water_before = total_water()
total_water_after = total_water()
call stepout(-1)

nstatis = nstat/nstatfrq
nstat = nstatis * nstatfrq
nstatsteps = 0
call t_stopf ('initialize')
!------------------------------------------------------------------
!   Main time loop    
!------------------------------------------------------------------

do while(nstep.lt.nstop.and.nelapse.gt.0) 
        
  nstep = nstep + 1
  time = time + dt
  day = day0 + nstep*dt/86400.
  nelapse = nelapse - 1
!------------------------------------------------------------------
!  Check if the dynamical time step should be decreased 
!  to handle the cases when the flow being locally linearly unstable
!------------------------------------------------------------------

  ncycle = 1

  call kurant()

  total_water_before = total_water()
  total_water_evap = 0.
  total_water_prec = 0.
  total_water_ls = 0.

  do icyc=1,ncycle

     icycle = icyc
     dtn = dt/ncycle
     dt3(na) = dtn
     dtfactor = dtn/dt

     if(mod(nstep,nstatis).eq.0.and.icycle.eq.ncycle) then
        nstatsteps = nstatsteps + 1
        dostatis = .true.
        if(masterproc) print *,'Collecting statistics...'
     else
        dostatis = .false.
     endif

     !bloss:make special statistics flag for radiation,since it's only updated at icycle==1.
     dostatisrad = .false.
     if(mod(nstep,nstatis).eq.0.and.icycle.eq.1) dostatisrad = .true.

!---------------------------------------------
!  	the Adams-Bashforth scheme in time

     call abcoefs()
 
!---------------------------------------------
!  	initialize stuff: 
	
     call zero()

!-----------------------------------------------------------
!       Buoyancy term:
	     
     call buoyancy()

!------------------------------------------------------------

     total_water_ls =  total_water_ls - total_water()

!------------------------------------------------------------
!       Large-scale and surface forcing:

     call forcing()

!----------------------------------------------------------
!       Nadging to sounding:

     call nudging()

!----------------------------------------------------------
!   	spange-layer damping near the upper boundary:

     if(dodamping) call damping()

!----------------------------------------------------------

     total_water_ls =  total_water_ls + total_water()

!-----------------------------------------------------------
!       Radiation

      if(dolongwave.or.doshortwave) then
        call radiation()
      end if

!----------------------------------------------------------
!     Update scalar boundaries after large-scale processes:

     call boundaries(2)

!-----------------------------------------------
!     surface fluxes:

     if(dosurface) call surface()

!-----------------------------------------------------------
!  SGS physics:

     if (dosgs) call sgs_proc()

!----------------------------------------------------------
!     Fill boundaries for SGS diagnostic fields:

     call boundaries(4)
!-----------------------------------------------
!       advection of momentum:

     call advect_mom()

!----------------------------------------------------------
!	SGS effects on momentum:

     if(dosgs) call sgs_mom()

!-----------------------------------------------------------
!       Coriolis force:
	     
     if(docoriolis) call coriolis()
	 
!---------------------------------------------------------
!       compute rhs of the Poisson equation and solve it for pressure. 

     call pressure()

!---------------------------------------------------------
!       find velocity field at n+1/2 timestep needed for advection of scalars:
!  Note that at the end of the call, the velocities are in nondimensional form.
	 
     call adams()

!---------------------------------------------------------
!      advection of scalars :

     call advect_all_scalars()

!---------------------------------------------------------
!   Ice fall-out
   
      if(docloud) then
          call ice_fall()
      end if

!----------------------------------------------------------
!     Update boundaries for scalars to prepare for SGS effects:

     call boundaries(3)
   
!---------------------------------------------------------
!      SGS effects on scalars :

     if (dosgs) call sgs_scalars()

!-----------------------------------------------------------
!       Handle upper boundary for scalars

     if(doupperbound) call upperbound()

!-----------------------------------------------------------
!       Cloud condensation/evaporation and precipitation processes:

      if(docloud.or.dosmoke) call micro_proc()

!----------------------------------------------------------
!  Tracers' physics:

      call tracers_physics()

!-----------------------------------------------------------
!    Compute diagnostic fields:

      call diagnose()

!----------------------------------------------------------

! Rotate the dynamic tendency arrays for Adams-bashforth scheme:

      nn=na
      na=nc
      nc=nb
      nb=nn

   end do ! icycle	
          
  total_water_after = total_water()
!----------------------------------------------------------
!  collect statistics, write save-file, etc.

   call stepout(nstatsteps)
  
!----------------------------------------------------------


end do ! main loop

!----------------------------------------------------------
!----------------------------------------------------------

   call t_stopf('total')
   if(masterproc) call t_prf(rank)

if(masterproc) write(*,*) 'Finished with SAM, exiting...'
call task_stop()

end program crm
