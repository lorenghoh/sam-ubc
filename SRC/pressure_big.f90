
! A version with non-blocking receive and blocking send. Doesn't need 
! EAGER_LIMIT set for communication.

subroutine pressure_big
	
!   parallel pressure-solver for 3D large domains.
!   the FFT is done on vertical slabs, first in x then in y.
!   Each processor gets its own slab.
!   This routine should be used only when the number of processors is larger
!   then the number of level. Otherwise, use pressure_orig
!   (C) 2002, Marat Khairoutdinov
!   Update 2015:
!     This version is based on the version of  pressure_big.f90 code by Don Dazlich
!     It is designed to minimize MPI communication for very large number of precessors.
!     DD: The transpose form x-z slabs to y-z slabs has each process doing 
!         NSUBDOMAINS send/receive pairs. This is replaced with an inverse 
!         transpose back to grid space and a transpose to y-z slabs, each requiring
!         only SQRT(NSUBDOMAINS) send/receive pairs per process. Additionally, the transpose
!         to y-z slabs switches dimensions so that the fft is passed contiguous memory space.
!   Update: Feb 2022: rm ff variable. All array sizes are nx_s and ny_s. Don's send around the
!           Fourier coeffs that are zero or not needed. Impemented dowallx = .true. option
!   Update: Mar 2022: made it possible for maximum number of preocessors to be NX_GL when NY_gl is
!           smaller but NX_GL should be equal or  divisible by NY_GL. 


use vars
use params, only: dowallx, dowally, docolumn
implicit none

integer, parameter :: nslab = max(1,nsubdomains/ny_gl) !number of slabs in the vertical
integer, parameter :: nx_s=nx_gl/nsubdomains ! width of the x-slabs
integer, parameter :: ny_s=nslab*ny_gl/nsubdomains ! width of the y-slabs
integer, parameter :: nzm2 = nzm/nslab

! Slabs:
real fx(nx_gl, ny_s, nzm2) ! slab for x-pass Fourier coefs
real fy(ny_gl, nx_s, nzm) ! slab for y-pass Fourier coefs
real(8) gx(nx_gl+2, ny_s) ! array to perform FFT in x
real(8) gy(ny_gl+2, nx_s) ! array to perform FFT in y

! Message buffers:
real bufx1(nx, ny_s, nzm2)
real bufx2(nx, ny_s, nzm2, max(1,nsubdomains_x))
real bufy1(nx_s, ny, nzm)
real bufy2(nx_s, ny, nzm, max(1,nsubdomains_y))
! save some memory
equivalence (bufx1(1,1,1),bufy1(1,1,1))
equivalence (bufx2(1,1,1,1),bufy2(1,1,1,1))

! rhs and solution (without extra boundaries)
real ppp(nx,ny,nzm)

! FFT stuff:
real(8) work(max((nx_gl+3)*(ny_s+1),(nx_s+1)*(ny_gl+2)))
real(8) trigxi(3*nx_gl/2+1),trigxj(3*ny_gl/2+1)
integer ifaxj(100),ifaxi(100)

! Tri-diagonal matrix solver coefficients:
real(8) a(nzm),b(nx_s,ny_gl),c(nzm),e	
real(8) xi,xj,xnx,xny,ddx2,ddy2,pii,factx,facty,eign(nx_s,ny_gl)
real(8) alfa(nx_s,ny_gl,nzm),beta(nx_s,ny_gl,nzm)

integer reqs_in(nsubdomains)
integer i, j, k, id, jd, m, n, it, jt, tag
integer irank, rnk
integer n_in, count
logical flag(nsubdomains)

! for wrapping p
real buff_ew1(ny,nzm), buff_ns1(nx,nzm)
real buff_ew2(ny,nzm), buff_ns2(nx,nzm)
integer rf, tagrf
logical waitflag

! Make sure that the grid is suitable for the solver:

if(mod(nx_gl,nsubdomains).ne.0) then
  if(masterproc) print*,'pressure_big: nx_gl/nsubdomains is not round number. STOP'
  call task_abort
endif
if(mod(nx_gl,ny_gl).ne.0) then
  if(masterproc) print*,'pressure_big: nx_gl/ny_gl is not round number. STOP'
  call task_abort
endif
if(nslab.gt.1.and.mod(nzm,nslab).ne.0) then
  if(masterproc) print*,'pressure_big: nzm should be divisible by (nx_gl/ny_gl). STOP'
  call task_abort
end if

!-----------------------------------------------------------------

if(docolumn) return
if(RUN2D) then
  print*,'pressure3D() cannot be called for 2D domains. Quitting...'
  call task_abort()
endif

!==========================================================================
!  Compute the r.h.s. of the Poisson equation for pressure

call press_rhs()

! variable p will also be used as grid decomposition placeholder between transposes
!   for the fourier coefficients

!==========================================================================
!   Form the vertical slabs (x-z) of right-hand-sides of Poisson equation 
!   for the FFT - one slab per a processor.

   ppp(1:nx,1:ny,1:nzm) = p(1:nx,1:ny,1:nzm)
   call transpose_x(ppp)

!==========================================================================
! Perform Fourier transformation n x-direction for a slab:

 call fftfax_crm(nx_gl,ifaxi,trigxi)

if(dowallx) then
 do k=1,nzm2
  gx(1:nx_gl,1:ny_s) = fx(1:nx_gl,1:ny_s,k)
  call cosft_crm(gx,work,trigxi,ifaxi,1,nx_gl+2,nx_gl,ny_s,-1)
  fx(1:nx_gl,1:ny_s,k) = gx(1:nx_gl,1:ny_s)
 end do
else
 do k=1,nzm2
  gx(1:nx_gl,1:ny_s) = fx(1:nx_gl,1:ny_s,k)
  call fft991_crm(gx,work,trigxi,ifaxi,1,nx_gl+2,nx_gl,ny_s,-1)
  fx(1,1:ny_s,k) = gx(1,1:ny_s)
  fx(2:nx_gl,1:ny_s,k) = gx(3:nx_gl+1,1:ny_s)
 end do
end if

call task_barrier()

!==========================================================================
!   Form the vertical slabs (y-z) of Fourier coefs  
!   for the FFT - in y, one slab per a processor.

call transpose_x_inv(ppp)
call transpose_y()

call fftfax_crm(ny_gl,ifaxj,trigxj)

do k=1,nzm
   gy(1:ny_gl,1:nx_s) = fy(1:ny_gl,1:nx_s,k)
   if(dowally) then
    call cosft_crm(gy,work,trigxj,ifaxj,1,ny_gl+2,ny_gl,nx_s,-1)
    fy(1:ny_gl,1:nx_s,k) = gy(1:ny_gl,1:nx_s)
   else
    call fft991_crm(gy,work,trigxj,ifaxj,1,ny_gl+2,ny_gl,nx_s,-1)
    fy(1,1:nx_s,k) = gy(1,1:nx_s)
    fy(2:ny_gl,1:nx_s,k) = gy(3:ny_gl+1,1:nx_s)
   end if
end do 


!==========================================================================
!   Solve the tri-diagonal system for Fourier coeffiecients 
!   in the vertical for each slab:


do k=1,nzm
    a(k)=rhow(k)/rho(k)/(adz(k)*adzw(k)*dz*dz)
    c(k)=rhow(k+1)/rho(k)/(adz(k)*adzw(k+1)*dz*dz)	 
end do 

ddx2=dx*dx
ddy2=dy*dy
pii = dacos(-1.d0)
xny=ny_gl      
xnx=nx_gl
it=rank*nx_s
do j=1,ny_gl
 if(dowally) then
    jd=j-1
    facty = 1.d0
 else
    if(j.eq.1) then
     jd = 0
    else
     jd=(j+1-0.1)/2.
    end if
    facty = 2.d0
 end if
 xj=jd
 do i=1,nx_s
   if(dowallx) then
      id=i+it-1
      factx = 1.d0
   else
      if(i+it.eq.1) then
       id = 0
      else
       id=(i+it+1-0.1)/2.
      end if
      factx = 2.d0
   end if
   xi=id
   eign(i,j)=(2.d0*cos(factx*pii/xnx*xi)-2.d0)/ddx2+ &
            (2.d0*cos(facty*pii/xny*xj)-2.d0)/ddy2
   if(id+jd.eq.0) then
     b(i,j)=eign(i,j)-a(1)-c(1)
     alfa(i,j,1)=-c(1)/b(i,j)
     beta(i,j,1)=fy(j,i,1)/b(i,j)
   else
     b(i,j)=eign(i,j)-c(1)
     alfa(i,j,1)=-c(1)/b(i,j)
     beta(i,j,1)=fy(j,i,1)/b(i,j)
   end if
 end do
end do

do k=2,nzm-1
 do j=1,ny_gl
  do i=1,nx_s
    e=eign(i,j)-a(k)-c(k)+a(k)*alfa(i,j,k-1)
    alfa(i,j,k)=-c(k)/e
    beta(i,j,k)=(fy(j,i,k)-a(k)*beta(i,j,k-1))/e
  end do
 end do
end do

do j=1,ny_gl
  do i=1,nx_s
     fy(j,i,nzm)=(fy(j,i,nzm)-a(nzm)*beta(i,j,nzm-1))/ &
                (eign(i,j)-a(nzm)+a(nzm)*alfa(i,j,nzm-1))
  end do
end do

do k=nzm-1,1,-1
  do j=1,ny_gl
    do i=1,nx_s
       fy(j,i,k)=alfa(i,j,k)*fy(j,i,k+1)+beta(i,j,k)
    end do
  end do
end do


!==========================================================================
! Perform inverse Fourier transf in y-direction for a slab:

call fftfax_crm(ny_gl,ifaxj,trigxj)

 do k=1,nzm
   if(dowally) then
     gy(1:ny_gl,1:nx_s) = fy(1:ny_gl,1:nx_s,k)
     gy(ny_gl+1:ny_gl+2,1:nx_s) = 0.
     call cosft_crm(gy,work,trigxj,ifaxj,1,ny_gl+2,ny_gl,nx_s,1)
   else
     gy(1,1:nx_s) = fy(1,1:nx_s,k)
     gy(2,1:nx_s) = 0.
     gy(3:ny_gl+1,1:nx_s) = fy(2:ny_gl,1:nx_s,k)
     gy(ny_gl+2,1:nx_s) = 0.
     call fft991_crm(gy,work,trigxj,ifaxj,1,ny_gl+2,ny_gl,nx_s,1)
   end if
   fy(1:ny_gl,1:nx_s,k) = gy(1:ny_gl,1:nx_s)
 end do

call task_barrier()

!==========================================================================
!   Form the vertical slabs (x-z) of Fourier coefs
!   for the inverse FFT - in x, one slab per a processor.


call transpose_y_inv()
call transpose_x(ppp)

! Perform inverse Fourier transform n x-direction for a slab:

 call fftfax_crm(nx_gl,ifaxi,trigxi)

if(dowallx) then
 do k=1,nzm2
  gx(1:nx_gl,1:ny_s) = fx(1:nx_gl,1:ny_s,k)
  gx(nx_gl+1:nx_gl+2,:) = 0.
  call cosft_crm(gx,work,trigxi,ifaxi,1,nx_gl+2,nx_gl,ny_s,1)
  fx(1:nx_gl,1:ny_s,k) = gx(1:nx_gl,1:ny_s)
 end do
else
 do k=1,nzm2
  gx(1,1:ny_s) = fx(1,1:ny_s,k)
  gx(2,1:ny_s) = 0.
  gx(3:nx_gl+1,1:ny_s) = fx(2:nx_gl,1:ny_s,k)
  gx(nx_gl+2,1:ny_s) = 0.
  call fft991_crm(gx,work,trigxi,ifaxi,1,nx_gl+2,nx_gl,ny_s,1)
  fx(1:nx_gl,1:ny_s,k) = gx(1:nx_gl,1:ny_s)
 end do
end if

call task_barrier()

call transpose_x_inv(ppp)

p(1:nx,1:ny,1:nzm) = ppp(:,:,:)

!==========================================================================
!  Update the pressure fields in the subdomains
!

! when we cut back on the ffts wrap p here - look to sib for the model
!DD temporary measure for dompi
  if(dompi) then
!DD custom build a wrap that sends from the north and east edges, and receives at the
!DD    south and west edges

      if(rank==rankee) then
        p(0,1:ny,1:nzm) = p(nx,1:ny,1:nzm)
      else
        call task_receive_float(buff_ew2(:,:),ny*nzm,reqs_in(1))
          buff_ew1(1:ny,1:nzm) = p(nx,1:ny,1:nzm)
          call task_bsend_float(rankee,buff_ew1(:,:),ny*nzm,1)
          waitflag = .false.
        do while (.not.waitflag)
            call task_test(reqs_in(1),waitflag,rf,tagrf)
        end do
        call task_barrier()
        p(0,1:ny,1:nzm) = buff_ew2(1:ny,1:nzm)
      endif

      if(rank==ranknn) then
         p(:,1-YES3D,1:nzm) = p(:,ny,1:nzm)
      else
             call task_receive_float(buff_ns2(:,:),nx*nzm,reqs_in(1))
             buff_ns1(1:nx,1:nzm) = p(1:nx,ny,1:nzm)
             call task_bsend_float(ranknn,buff_ns1(:,:),nx*nzm,1)
             waitflag = .false.
         do while (.not.waitflag)
               call task_test(reqs_in(1),waitflag,rf,tagrf)
             end do
         call task_barrier()
         p(1:nx,1-YES3D,1:nzm) = buff_ns2(1:nx,1:nzm)
      endif

  else
    p(0,:,1:nzm) = p(nx,:,1:nzm)
    p(:,1-YES3D,1:nzm) = p(:,ny,1:nzm)
  endif

  ! overwrite for the case of walls in x and y
  call task_rank_to_index(rank,it,jt)
  if(dowallx.and.it.eq.0) then
    p(0,:,1:nzm) = p(1,:,1:nzm)
  end if
  if(dowally.and.jt.eq.0) then
    p(:,1-YES3D,1:nzm) = p(:,1,1:nzm)
  end if

!DD end ugly wrap code.

!==========================================================================
!  Add pressure gradient term to the rhs of the momentum equation:

call press_grad()

!==========================================================================
!==========================================================================
!==========================================================================

contains

   
!==========================================================================
   subroutine transpose_x(pm)

! transpose from blocks to x-z slabs

      REAL, INTENT(in) :: pm(nx, nslab*ny, nzm2)

      irank = rank-mod(rank,nsubdomains_x)

      n_in = 0
      do m = irank, irank+nsubdomains_x-1

        if(m.ne.rank) then

          n_in = n_in + 1
          call task_receive_float(bufx2(1:nx,1:ny_s,1:nzm2,n_in),nx*ny_s*nzm2,reqs_in(n_in))
          flag(n_in) = .false.

        end if

      end do ! m

      do m = irank, irank+nsubdomains_x-1

        if(m.ne.rank) then

          n = m-irank

          bufx1(1:nx,1:ny_s,1:nzm2) = pm(1:nx,n*ny_s+1:n*ny_s+ny_s,1:nzm2)
          call task_bsend_float(m,bufx1(1:nx,1:ny_s,1:nzm2),nx*ny_s*nzm2, 33)

        endif

      end do ! m


! don't sent a buffer to itself, just fill directly.

      n = rank-irank
      call task_rank_to_index(rank,it,jt)
      fx(1+it:nx+it,1:ny_s,1:nzm2) = pm(1:nx,n*ny_s+1:n*ny_s+ny_s,1:nzm2)


      ! Fill slabs when receive buffers are full:

      count = n_in
      do while (count .gt. 0)
        do m = 1,n_in
         if(.not.flag(m)) then
            call task_test(reqs_in(m), flag(m), rnk, tag)
              if(flag(m)) then
                 count=count-1
                 call task_rank_to_index(rnk,it,jt)
                 fx(1+it:nx+it,1:ny_s,1:nzm2) = bufx2(1:nx,1:ny_s,1:nzm2,m)
              endif
          endif
         end do
      end do
      call task_barrier()

   end subroutine transpose_x

!==========================================================================
   subroutine transpose_x_inv(pm)

! transpose from x-z slabs to blocks

      REAL, INTENT(out) :: pm(nx, nslab*ny, nzm2)

      irank = rank-mod(rank,nsubdomains_x)
      n_in = 0
      do m = irank, irank+nsubdomains_x-1

        if(m.ne.rank) then

          n_in = n_in + 1
          call task_receive_float(bufx2(1:nx,1:ny_s,1:nzm2,n_in),nx*ny_s*nzm2,reqs_in(n_in))
          flag(n_in) = .false.

        endif

      end do ! m

      do m = irank, irank+nsubdomains_x-1

        if(m.ne.rank) then

          call task_rank_to_index(m,it,jt)
          bufx1(1:nx,1:ny_s,1:nzm2) = fx(1+it:it+nx,1:ny_s,1:nzm2)
          call task_bsend_float(m,bufx1(1:nx,1:ny_s,1:nzm2),nx*ny_s*nzm2, 33)

        endif

      end do ! m

! don't sent a buffer to itself, just fill directly.

      n = rank-irank
      call task_rank_to_index(rank,it,jt)
      pm(1:nx,n*ny_s+1:n*ny_s+ny_s,1:nzm2) = fx(1+it:nx+it,1:ny_s,1:nzm2)

! Fill slabs when receive buffers are full:

      count = n_in
      do while (count .gt. 0)
        do m = 1,n_in
         if(.not.flag(m)) then
              call task_test(reqs_in(m), flag(m), rnk, tag)
              if(flag(m)) then
                 count=count-1
                 n = rnk-irank
                 pm(1:nx,n*ny_s+1:n*ny_s+ny_s,1:nzm2) = bufx2(1:nx,1:ny_s,1:nzm2,m)
              endif
         endif
        end do
      end do

      call task_barrier()
   end subroutine transpose_x_inv


!==========================================================================
   subroutine transpose_y()

! transpose from blocks to y-z slabs

      irank = rank / nsubdomains_y  

      n_in = 0
      do m = irank, nsubdomains-1, nsubdomains_x

        if(m.ne.rank) then

          n_in = n_in + 1
          call task_receive_float(bufy2(:,:,:,n_in),ny*nx_s*nzm,reqs_in(n_in))
          flag(n_in) = .false.
          
        else
! don't sent a buffer to itself, just fill directly.

          n = mod(rank,nsubdomains_y) 
          call task_rank_to_index(rank,it,jt)
          do i = 1,nx_s	  
            fy(1+jt:ny+jt,i,1:nzm) = ppp(n*nx_s+i,1:ny,1:nzm)
          enddo

        end if

      end do ! m

      irank = nsubdomains_y*mod(rank,nsubdomains_x)
      do m = irank, irank+nsubdomains_y-1

        if(m.ne.rank) then

          n = m-irank

          bufy1(:,:,:) = ppp(n*nx_s+1:n*nx_s+nx_s,1:ny,1:nzm)
          call task_bsend_float(m,bufy1(:,:,:),ny*nx_s*nzm, 33) 

        endif

      end do ! m


      ! Fill slabs when receive buffers are full:

      count = n_in
      do while (count .gt. 0)
        do m = 1,n_in
         if(.not.flag(m)) then
      	    call task_test(reqs_in(m), flag(m), rnk, tag)
              if(flag(m)) then 
            	 count=count-1
                 call task_rank_to_index(rnk,it,jt)
                 do i = 1,nx_s	  
                   fy(1+jt:ny+jt,i,1:nzm) = bufy2(i,1:ny,1:nzm,m)
                 enddo
              endif   
          endif
         end do
      end do

      call task_barrier()
   end subroutine transpose_y
   
!==========================================================================
   subroutine transpose_y_inv()

! transpose from y-z slabs to blocks

      n_in = 0
      irank = nsubdomains_y*mod(rank,nsubdomains_x)
      do m = irank, irank+nsubdomains_y-1

        if(m.ne.rank) then

          n_in = n_in + 1
          call task_receive_float(bufy2(:,:,:,n_in),ny*nx_s*nzm,reqs_in(n_in))
          flag(n_in) = .false.
          
        else

! don't sent a buffer to itself, just fill directly.

          n = rank-irank
          call task_rank_to_index(rank,it,jt)
          do i = 1,nx_s
            ppp(n*nx_s+i,1:ny,1:nzm) = fy(1+jt:ny+jt,i,1:nzm) 
          enddo 

        endif

      end do ! m

      irank = rank / nsubdomains_y  
      do m = irank, nsubdomains-1, nsubdomains_x

        if(m.ne.rank) then

          call task_rank_to_index(m,it,jt)
          do i = 1,nx_s
            bufy1(i,:,:) = fy(1+jt:jt+ny,i,1:nzm)
          enddo
          call task_bsend_float(m,bufy1(:,:,:),ny*nx_s*nzm, 33)

        endif

      end do ! m

! Fill slabs when receive buffers are full:

      irank = nsubdomains_y*mod(rank,nsubdomains_x)
      count = n_in
      do while (count .gt. 0)
        do m = 1,n_in
         if(.not.flag(m)) then
              call task_test(reqs_in(m), flag(m), rnk, tag)
              if(flag(m)) then
                 count=count-1
                 n = rnk-irank
                 ppp(n*nx_s+1:n*nx_s+nx_s,1:ny,1:nzm) = bufy2(1:nx_s,1:ny,1:nzm,m)
              endif
         endif
        end do
      end do

      call task_barrier()
   end subroutine transpose_y_inv

end subroutine pressure_big
