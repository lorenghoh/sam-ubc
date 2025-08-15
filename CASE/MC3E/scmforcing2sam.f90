program scmforcing2sam

implicit none

! physical constants (values from sam)
real*8, parameter :: g = 9.81, cp=1004., rv=461., rgas=287.,    &
                     kappa = rgas / cp, epsv =0.61

! netcdf variables
integer rc, ncid, varid, dimid
integer, dimension(2) :: start, count

! index variables
integer n, l

! dimensions
integer, parameter :: ntime = 368 ! number of time records
integer, parameter :: nlev =180 ! number of vertical levels

! data variables 
real*8, dimension(1) :: time
real*4, dimension(1) :: ps, tvbar, exners
real*4, allocatable, dimension(:) :: z, t, q, u, v,rtadvt, rqadvt, w, tv, level
real*8, allocatable, dimension(:) :: p
real*4, allocatable, dimension(:) :: exner
integer :: yy, mm, dd, hh, k
real*4 :: xtem
real*4 :: us, vs, ts, qs, rtadvts, rqadvts


!open input files files
  open(unit=41,file='srcdata/fields.v1a.txt',form='formatted')
  open(unit=42,file='srcdata/lsf.v1a.txt',form='formatted')


! allocate variables and initialize as necessary
allocate(p(nlev))
allocate(z(nlev))
allocate(t(nlev))
allocate(tv(nlev))
allocate(q(nlev))
allocate(u(nlev))
allocate(v(nlev))
allocate(rtadvt(nlev))
allocate(rqadvt(nlev))
allocate(exner(nlev))
allocate(w(nlev))
w = 0.0

!open sam files (these need header added to line 1)
open(unit=11,file='./snd',form='formatted')
open(unit=12,file='./lsf_hor',form='formatted')


do n = 1,ntime
 ! time record
 read(41,*) yy, mm, dd, hh
 read(42,*) yy, mm, dd, hh
 if(mm == 4) then
   time = 90. + real(dd,8) + real(hh,8)/8.
 else if (mm==5) then
   time = 120. + real(dd,8) + real(hh,8)/8.
 else if (mm == 6) then
   time = 151. + real(dd,8) + real(hh,8)/8.
 endif
 read(41,*)ps, xtem, us, vs, xtem, xtem, ts, qs, xtem
 read(42,*)ps, rtadvts, xtem, rqadvts
 do k = 1,nlev
   read(41,*)p(k), z(k), u(k), v(k), w(k), xtem, t(k), q(k), xtem
   read(42,*)p(k), rtadvt(k) , xtem, rqadvt(k)
 enddo
  ! convert to sam units 
  w = w * 100./3600.                ! convert omega to Pa/s
  rqadvt(:) = rqadvt(:) / 1000. ! convert to kg/kg/s
  rqadvts = rqadvts / 1000.
  
  
  
  ! write sam files
  write(11,'(f12.4,i16,f11.4,a21)')time,nlev,ps,'     day,levels,pres0'
  write(12,'(f12.4,i16,f11.4,a21)')time,nlev,ps,'     day,levels,pres0'
  z = -999.
  do l = 1,nlev
    if(p(l) < ps(1)) then
       write(11,'(6f10.3)')z(l),p(l),t(l),q(l),u(l),v(l)
       write(12,'(2f10.3,2e14.5,3f10.3)')z(l),p(l),rtadvt(l),rqadvt(l),  &
                                      u(l),v(l),w(l)
    else
       write(11,'(6f10.3)')z(l),p(l),ts,qs,us,vs
       write(12,'(2f10.3,2e14.5,3f10.3)')z(l),p(l),rtadvts,rqadvts,  &
                                      us,vs,0.0
    endif
  enddo
enddo

!close files
close(11)
close(12)
close(41)
close(42)

end program scmforcing2sam
