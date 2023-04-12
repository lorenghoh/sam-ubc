! module to collect timing information of a running program
! to replace C-based t_* code
! (C) Marat Khairoutdinov, October 2019

module timer

implicit none

! I decided not to use derived data types (too verbose) =MK
integer, parameter :: NTIMERS_MAX = 100 !maximum allowed number of timers
character(20):: name(NTIMERS_MAX) = ' ' ! unique name of the given timer
logical:: active(NTIMERS_MAX) =.false. ! to tell if it is active now
integer:: parent(NTIMERS_MAX) = 0      ! parrent (called from) timer
real:: time_start(NTIMERS_MAX)  ! time when t_startf is called 
real:: time(NTIMERS_MAX) = 0. ! accumulated time
integer:: ncalls(NTIMERS_MAX) = 0 ! number of calls of t_startf

integer:: ntimers = 0 ! number of created timers
integer:: current_active   ! currently active timer
integer:: current_parent   ! self-explanatory

logical:: initialized = .false.
real:: time0 ! time of initialization of timers

contains

recursive subroutine print_children(iparent,fmt)
implicit none
integer, intent(in) :: iparent
character(*) ::  fmt
integer i
! find children
do i=1,ntimers
  if(parent(i).eq.iparent) then
   write(*,'('//fmt//')') trim(name(i)),ncalls(i),time(i),time(i)/time(1)*100.
   call print_children(i,'TR3,'//trim(fmt))
  end if
end do
end subroutine print_children


subroutine t_error(str)
 implicit none
 character(*) str
  print*,str,current_active,name(current_active)
  stop
end subroutine t_error

end module timer

