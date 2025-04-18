subroutine task_init
		
!   Check things, initialize multitasking:

use grid
implicit none

integer itasks,ntasks

if(YES3D .ne. 1 .and. YES3D .ne. 0) then
  print*,'YES3D is not 1 or 0. STOP'
  stop
endif

if(YES3D .eq. 1 .and. ny_gl .lt. 4) then
  print*,'ny_gl is too small for a 3D case.STOP'
  stop
endif

if(YES3D .eq. 0 .and. ny_gl .ne. 1) then
  print*,'ny_gl should be 1 for a 2D case. STOP'
  stop
endif

if(nsubdomains.eq.1) then

  rank =0
  ntasks = 1
  dompi = .false.

else

  call task_start(rank, ntasks)

  dompi = .true.

!  call execute_command_line('hostname')

  if(ntasks.ne.nsubdomains) then
    if(rank.eq.0) print *,'number of processors is not equal to nsubdomains!',&
             '  ntasks=',ntasks,'   nsubdomains=',nsubdomains
    call task_abort() 
  endif
        
  call task_barrier()

  call task_ranks()
        
end if ! nsubdomains.eq.1

do itasks=0,nsubdomains-1
   call task_barrier()
   if(itasks.eq.rank) then
    open(8,file='./CaseName',status='old',form='formatted')
    read(8,'(a)') case
    close (8)
   endif
end do

masterproc = rank.eq.0

if(masterproc) print *,'number of MPI tasks:',ntasks
	

end
