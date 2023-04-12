subroutine t_setoptionf (i, j)
integer i,j
end

subroutine t_initializef()
 use timer
 if(initialized) call t_error("t_initialize should be called only once")
 call cpu_time(time0)
 initialized = .true.
 ntimers = 1
 name(1) = "_total"
 active(1) = .true.
 current_active = 1
 current_parent = 0
end

subroutine t_startf(str)
 use timer
 use grid, only: masterproc
 implicit none
 character(*), intent(in) :: str
 integer id
 logical:: found 
 real ttt
 call task_barrier()
 !if(masterproc) print*,'>>>start  ',trim(str)
 found = .false.
 call cpu_time(ttt)
 if(.not.initialized) call t_error("t_initialize should be called first")
 do id=1,ntimers
  if(trim(str).eq.trim(name(id)).and.parent(id).eq.current_active) then
    found = .true.
    exit
  end if
 end do
 if(.not.found) then
  ntimers = ntimers+1
  id = ntimers
  name(id) = trim(str)
 end if
 if(active(id)) call t_error('timer '//str//' is already active.')
 parent(id) = current_active
 current_parent = parent(id)
 active(id) = .true.
 current_active = id
 time_start(id) = ttt
end 

subroutine t_stopf(str)
 use timer
 implicit none
 character(*), intent(in) ::  str
 integer id
 real ttt
 logical:: found 
 call task_barrier()
 found = .false.
 do id=1,ntimers
  if(trim(name(id)).eq.trim(str).and.parent(id).eq.current_parent) then
   if(.not.active(id)) &
     call t_error('t_stopf for time '//str//' is called without t_startf first')
   found = .true.
   active(id) = .false.
   current_active = parent(id)
   current_parent = parent(current_active)
   call cpu_time(ttt)
   time(id) = time(id) + ttt-time_start(id)
   ncalls(id) = ncalls(id)+1
   exit
  end if
 end do
 if(.not.found) &
  call t_error('timer '//str//' cannot be stopped as it was not started')
end 

subroutine t_prf(rank)
use timer
use grid, only: masterproc
implicit none
integer, intent(in) :: rank
real ttt
integer i, iparents
call cpu_time(ttt)
time(1) = ttt-time0
if(.not.masterproc) return
print*,'------------------------------------------------------'
print*,'ntimers=',ntimers
print*,'------------ timer ----------Ncalls--------time(s)-------- % of total----'
call print_children(1,'a,T30,i9,T40,g15.7,T60,f6.2') 
ttt = 0.
do i=1,ntimers
 if(parent(i).eq.1) ttt = ttt+time(i)
end do
write(*,'(a,T40,g15.7,T60,f6.2)') 'TOTAL:',ttt,ttt/time(1)*100.
end 



