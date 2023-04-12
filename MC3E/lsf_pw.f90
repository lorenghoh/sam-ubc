program lsf_pw

implicit none

real pw(368), day,ps
integer k,l, lm
character*46 header
real*8 x(7)

open(unit=2,file='pw.dat',form='formatted')
open(unit=11,file='lsf.cont',form='formatted')
open(unit=12,file='lsf',form='formatted')

read(2,*)pw
close(2)

  read(11,'(a46)')header
  write(12,'(a46)')header
do k = 1,368
  read(11,*) day,lm,ps
  write(12,'(f11.4,i16,2f11.4)') day,lm,ps, pw(k)
  do l = 1,42
    read(11,*)x
    write(12,'(2f10.3,2e14.5,3f10.4)')x
  enddo
enddo

end program lsf_pw
