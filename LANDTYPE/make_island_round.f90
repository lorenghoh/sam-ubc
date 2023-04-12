! make input landtype and LAI files for a round islando

implicit none

integer, parameter :: nx=128
integer, parameter :: ny=128
real, parameter :: radius = 0.2 ! island radius relative to domain size
integer landtype(nx,ny)
real LAI(nx,ny)
integer i,j

do j=1,ny
 do i=1,nx
   if(((i-1.)-nx/2.)**2+((j-1)-ny/2.)**2.lt.(radius*nx)**2) then
      if(i.le.nx/2) then
          landtype(i,j) = 2
          LAI(i,j) = 5.
      else
          landtype(i,j) = 16
          LAI(i,j) = 0.
      end if
   else
      landtype(i,j) = 0
      LAI(i,j) = 0
   end if
 end do
end do
open(1,file='island_landtype.bin',form='unformatted')
write(1) nx
write(1) ny
write(1) landtype
close(1)
open(1,file='island_LAI.bin',form='unformatted')
write(1) nx
write(1) ny
write(1) LAI
close(1)

end

