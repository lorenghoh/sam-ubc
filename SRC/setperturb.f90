
subroutine setperturb

!  Random noise

use vars
use params
use microphysics, only: micro_field, index_water_vapor
use sgs, only: setperturb_sgs

implicit none

integer i,j,k,ptype,it,jt,n,ii,jj
real rrr,ranf_
real xxx,yyy,zzz,sss

call ranset_(3+nensemble)

ptype = perturb_type

call setperturb_sgs(ptype)  ! set sgs fields

call task_rank_to_index(rank,it,jt)

select case (ptype)

  case(-1)  ! no perturbation is set
 

  case(0)

     sss = 0.
     do k=1,nzm
      do jj=1,ny_gl
       j=jj-jt
       do ii=1,nx_gl
        i=ii-it 
        rrr=1.-2.*ranf_()
        sss = sss + rrr
        if(i.ge.1.and.i.le.nx.and.j.ge.1.and.j.le.ny) then
         if(k.le.5) then
            t(i,j,k)=t(i,j,k)+0.02*rrr*(6-k)
         endif
        end if
       end do
      end do
     end do
     if(masterproc) print*,'>>>>>> sss=',sss

  case(1)

     do k=1,nzm
      do jj=1,ny_gl
       j=jj-jt
       do ii=1,nx_gl
        i=ii-it 
        rrr=1.-2.*ranf_()
        if(i.ge.1.and.i.le.nx.and.j.ge.1.and.j.le.ny) then
         if(q0(k).gt.6.e-3) then
            t(i,j,k)=t(i,j,k)+0.1*rrr
         endif
        end if
       end do
      end do
     end do

  case(2) ! warm bubble

     if(masterproc) then
       print*, 'initialize with warm bubble:'
       print*, 'bubble_x0=',bubble_x0
       print*, 'bubble_y0=',bubble_y0
       print*, 'bubble_z0=',bubble_z0
       print*, 'bubble_radius_hor=',bubble_radius_hor
       print*, 'bubble_radius_ver=',bubble_radius_ver
       print*, 'bubble_dtemp=',bubble_dtemp
       print*, 'bubble_dq=',bubble_dq
     end if

     do k=1,nzm
       zzz = z(k)
       do j=1,ny
         yyy = dy*(j+jt)
         do i=1,nx
          xxx = dx*(i+it)
           if((xxx-bubble_x0)**2+YES3D*(yyy-bubble_y0)**2.lt.bubble_radius_hor**2 &
            .and.(zzz-bubble_z0)**2.lt.bubble_radius_ver**2) then
              rrr = cos(pi/2.*(xxx-bubble_x0)/bubble_radius_hor)**2 &
               *cos(pi/2.*(yyy-bubble_y0)/bubble_radius_hor)**2 &
               *cos(pi/2.*(zzz-bubble_z0)/bubble_radius_ver)**2
              t(i,j,k) = t(i,j,k) + bubble_dtemp*rrr
              micro_field(i,j,k,index_water_vapor) = &
                  micro_field(i,j,k,index_water_vapor) + bubble_dq*rrr
           end if
         end do
       end do
     end do

  case(3)   ! gcss wg1 smoke-cloud case

     do k=1,nzm
      do jj=1,ny_gl
       j=jj-jt
       do ii=1,nx_gl
        i=ii-it 
        rrr=1.-2.*ranf_()
        if(i.ge.1.and.i.le.nx.and.j.ge.1.and.j.le.ny) then
         if(q0(k).gt.0.5e-3) then
            t(i,j,k)=t(i,j,k)+0.1*rrr
         endif
        end if
       end do
      end do
     end do

  case(4)  ! gcss wg1 arm case

     do k=1,nzm
      do jj=1,ny_gl
       j=jj-jt
       do ii=1,nx_gl
        i=ii-it 
        rrr=1.-2.*ranf_()
        if(i.ge.1.and.i.le.nx.and.j.ge.1.and.j.le.ny) then
         if(z(k).le.200.) then
            t(i,j,k)=t(i,j,k)+0.1*rrr*(1.-z(k)/200.)
         endif
        end if
       end do
      end do
     end do

  case(5)  ! gcss wg1 BOMEX case

     do k=1,nzm
      do jj=1,ny_gl
       j=jj-jt
       do ii=1,nx_gl
        i=ii-it 
        rrr=1.-2.*ranf_()
        if(i.ge.1.and.i.le.nx.and.j.ge.1.and.j.le.ny) then
         if(z(k).le.1600.) then
            t(i,j,k)=t(i,j,k)+0.1*rrr
            micro_field(i,j,k,index_water_vapor)= &
                      micro_field(i,j,k,index_water_vapor)+0.025e-3*rrr
         endif
        end if
       end do
      end do
     end do

  case(6)  ! GCSS Lagragngian ASTEX


     do k=1,nzm
      do jj=1,ny_gl
       j=jj-jt
       do ii=1,nx_gl
        i=ii-it 
        rrr=1.-2.*ranf_()
        if(i.ge.1.and.i.le.nx.and.j.ge.1.and.j.le.ny) then
         if(q0(k).gt.6.e-3) then
            t(i,j,k)=t(i,j,k)+0.1*rrr
            micro_field(i,j,k,index_water_vapor)= &
                      micro_field(i,j,k,index_water_vapor)+2.5e-5*rrr
         endif
        end if
       end do
      end do
     end do


  case default

       if(masterproc) print*,'perturb_type is not defined in setperturb(). Exitting...'
       call task_abort()

end select


end

