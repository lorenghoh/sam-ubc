module tracers


! This module serves as a template for adding tracer transport in the model. The tracers can be 
! chemical tracers, or bin microphysics drop/ice categories, etc. 
! The number of tracers is set by the parameter ntracers which is set in domain.f90.
! Also, the logical flag dotracers should be set to .true. in namelist (default is .false.).
! The model will transport the tracers around automatically (advection and SGS diffusion).
! The user must supply the initialization in the subroutine tracers_init() in this module.
! By default, the surface flux of all tracers is zero. Nonzero values can be set in tracers_flux().
! The local sinks/sources of tracers should be supplied in tracers_physics().



 use grid
 implicit none

 real tracer  (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 0:ntracers) 
 real fluxbtr (nx, ny, 0:ntracers) ! surface flux of tracers
 real fluxttr (nx, ny, 0:ntracers) ! top boundary flux of tracers
 real trwle(nz,0:ntracers)  ! resolved vertical flux 
 real trwsb(nz,0:ntracers)  ! SGS vertical flux 
 real tradv(nz,0:ntracers)  ! tendency due to vertical advection 
 real trdiff(nz,0:ntracers)  ! tendency due to vertical diffusion 
 real trphys(nz,0:ntracers)  ! tendency due to physics 
 character *4 tracername(0:ntracers)
 character *10 tracerunits(0:ntracers)

CONTAINS

 subroutine tracers_init()

  integer k,ntr
  character *2 ntrchar
  integer, external :: lenstr

 tracer = 0.
 fluxbtr = 0.
 fluxttr = 0.

! Add your initialization code here. Default is to set to 0 in setdata.f90.

 if(nrestart.eq.0) then

! here ....

 end if

! Specify te tracers' default names:

   ! Default names are TRACER01, TRACER02, etc:

   do ntr = 1,ntracers
     write(ntrchar,'(i2)') ntr
     do k=1,3-lenstr(ntrchar)-1
        ntrchar(k:k)='0'
     end do
     tracername(ntr) = 'TR'//ntrchar(1:2)
     tracerunits(ntr) = '[TR]'
  end do

 end subroutine tracers_init



 subroutine tracers_flux()

! Set surface and top fluxes of tracers. Default is 0 set in setdata.f90

  ! UBC ENT
  fluxbtr(:, :, 1) = 1.
  ! End UBC ENT 

 end subroutine tracers_flux


 
 subroutine tracers_physics()

 ! add here a call to a subroutine that does something to tracers besides advection and diffusion.
 ! The transport is done automatically. 

   trphys = 0. ! Default tendency due to physics. You code should compute this to output statistics.

   ! UBC ENT
   tracer(:, :, :, 1) = tracer(:, :, :, 1) - tracer(:, :, :, 1)/15./60.*dtn
   ! End UBC ENT 

 end subroutine tracers_physics



 subroutine tracers_hbuf_init(namelist,deflist,unitlist,status,average_type,count,trcount)

! Initialize the list of tracers statistics variables written in statistics.f90

   character(*) namelist(*), deflist(*), unitlist(*)
   integer status(*),average_type(*),count,trcount
   integer ntr


   do ntr=1,ntracers

     count = count + 1
     trcount = trcount + 1
     namelist(count) = trim(tracername(ntr))
     deflist(count) = trim(tracername(ntr))
     unitlist(count) = trim(tracerunits(ntr))
     status(count) = 1
     average_type(count) = 0		
     count = count + 1
     trcount = trcount + 1
     namelist(count) = trim(tracername(ntr))//'FLX'
     deflist(count) = 'Total flux of '//trim(tracername(ntr))
     unitlist(count) = trim(tracerunits(ntr))//' kg/m2/s'
     status(count) = 1
     average_type(count) = 0
     count = count + 1
     trcount = trcount + 1
     namelist(count) = trim(tracername(ntr))//'FLXS'
     deflist(count) = 'SGS flux of '//trim(tracername(ntr))
     unitlist(count) = trim(tracerunits(ntr))//' kg/m2/s'
     status(count) = 1
     average_type(count) = 0
     count = count + 1
     trcount = trcount + 1
     namelist(count) = trim(tracername(ntr))//'ADV'
     deflist(count) = 'Tendency of '//trim(tracername(ntr)//'due to vertical advection')
     unitlist(count) = trim(tracerunits(ntr))//'/day'
     status(count) = 1
     average_type(count) = 0
     count = count + 1
     trcount = trcount + 1
     namelist(count) = trim(tracername(ntr))//'DIFF'
     deflist(count) = 'Tendency of '//trim(tracername(ntr)//'due to vertical SGS transport')
     unitlist(count) = trim(tracername(ntr))//'/day'
     status(count) = 1
     average_type(count) = 0
     count = count + 1
     trcount = trcount + 1
     namelist(count) = trim(tracername(ntr))//'PHYS'
     deflist(count) = 'Tendency of '//trim(tracername(ntr)//'due to physics')
     unitlist(count) = trim(tracername(ntr))//'/day'
     status(count) = 1
     average_type(count) = 0
   end do

 end subroutine tracers_hbuf_init

end module tracers
