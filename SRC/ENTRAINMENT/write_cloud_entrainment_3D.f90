subroutine write_cloud_entrainment_3D

use vars
use entrainment

implicit none

character *120 filename
character *80 long_name
character *8 name
character *10 timechar
character *4 rankchar
character *5 sepchar
character *6 filetype
character *10 units
character *12 c_z(nzm),c_p(nzm),c_dx, c_dy, c_time
integer i,j,k,n,nfields,nfields1
real(4) tmp(nx,ny,nzm)

integer, external :: lenstr

nfields = 10 ! number of 3D fields to save
nfields1 = 0

if(masterproc.or.output_sep) then

  if(output_sep) then
     write(rankchar,'(i4)') rank
     sepchar="_"//rankchar(5-lenstr(rankchar):4)
  else
     sepchar=""
  end if
  write(rankchar,'(i4)') nsubdomains
  write(timechar,'(i10)') nstep
  do k=1,11-lenstr(timechar)-1
    timechar(k:k)='0'
  end do

  if(RUN3D) then
    if(save3Dbin) then
      filetype = '.bin3D'
    else
      filetype = '.com3D'
    end if
    filename='./OUT_3D/'//trim(case)//'_CLOUD_ENTRAIN_'//trim(caseid)//'_'// &
        rankchar(5-lenstr(rankchar):4)//'_'//timechar(1:10)//filetype//sepchar
    open(46,file=filename,status='unknown',form='unformatted', &
            BUFFERED='YES', & ! use for intel compiler
            ACTION='WRITE')
  else
    if(save3Dbin) then
     if(save3Dsep) then
       filetype = '.bin3D'
     else
       filetype = '.bin2D'
     end if
    else
     if(save3Dsep) then
       filetype = '.com3D'
     else
       filetype = '.com2D'
     end if
    end if
    if(save3Dsep) then
      filename='./OUT_3D/'//trim(case)//'_CLOUD_ENTRAIN_'//trim(caseid)//'_'// &
        rankchar(5-lenstr(rankchar):4)//'_'//timechar(1:10)//filetype//sepchar
      open(46,file=filename,status='unknown',form='unformatted')	
    else
      filename='./OUT_3D/'//trim(case)//'_CLOUD_ENTRAIN_'//trim(caseid)//'_'// &
        rankchar(5-lenstr(rankchar):4)//filetype//sepchar
      if(nrestart.eq.0.and.notopened3D) then
         open(46,file=filename,status='unknown',form='unformatted')	
      else
         open(46,file=filename,status='unknown', &
                              form='unformatted', position='append')
      end if
      notopened3D=.false.
    end if  

  end if

  if(masterproc) then

    if(save3Dbin) then

      write(46) nx,ny,nzm,nsubdomains,nsubdomains_x,nsubdomains_y,nfields
      do k=1,nzm
        write(46) real(z(k),4) 
      end do
      do k=1,nzm
        write(46) real(pres(k),4)
      end do
      write(46) real(dx,4)
      write(46) real(dy,4)
      write(46) real(float(nstep)*dt/(3600.*24.)+day0,4)

    else

      write(long_name,'(8i4)') nx,ny,nzm,nsubdomains, &
                                   nsubdomains_x,nsubdomains_y,nfields
      do k=1,nzm
         write(c_z(k),'(f12.3)') z(k)
      end do
      do k=1,nzm
         write(c_p(k),'(f12.3)') pres(k)
      end do
      write(c_dx,'(f12.5)') dx
      write(c_dy,'(f12.5)') dy
      write(c_time,'(f12.5)') nstep*dt/(3600.*24.)+day0
	
      write(46) long_name(1:32)
      write(46) c_time,c_dx,c_dy, (c_z(k),k=1,nzm),(c_p(k),k=1,nzm)

    end if ! save3Dbin

  end if ! masterproc
 
end if ! masterproc.or.output_sep

! JTD
  nfields1=nfields1+1
  tmp = volume_cloud(1:nx, 1:ny, 1:nzm)
  name='VTETCLD'
  long_name='Tetrahedral Cloud Volume'
  units='m3/m3'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  tmp = w_adv_fluxes(1:nx, 1:ny, 1:nzm, 1)*w_area_cloud_frac(1:nx, 1:ny, 1:nzm)
  name='MFTETCLD'
  long_name='Tetrahedral Cloud Vertical Mass Flux'
  units='kg/s/m2'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                  save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  tmp = E_tetra_cloud(1:nx, 1:ny, 1:nzm, 1)
  name='ETETCLD'
  long_name='Tetrahedral Cloud Mass Entrainment'
  units='kg/m3/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  tmp = D_tetra_cloud(1:nx, 1:ny, 1:nzm, 1)
  name='DTETCLD'
  long_name='Tetrahedral Cloud Mass Detrainment'
  units='kg/m3/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                  save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  tmp = E_tetra_cloud(1:nx, 1:ny, 1:nzm, 2)
  name='EQTETCLD'
  long_name='Tetrahedral Cloud Specific Humidity Entrainment'
  units='kg/m3/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                  save3Dbin,dompi,rank,nsubdomains)
 
  nfields1=nfields1+1
  tmp = D_tetra_cloud(1:nx, 1:ny, 1:nzm, 2)
  name='DQTETCLD'
  long_name='Tetrahedral Cloud Specific Humidity Detrainment'
  units='kg/m3/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                  save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  tmp = E_tetra_cloud(1:nx, 1:ny, 1:nzm, 3)
  name='ETTETCLD'
  long_name='Tetrahedral Cloud Liquid Water Moist Static Energy Entrainment'
  units='K kg/m3/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                  save3Dbin,dompi,rank,nsubdomains)
  
  nfields1=nfields1+1
  tmp = D_tetra_cloud(1:nx, 1:ny, 1:nzm, 3)
  name='DTTETCLD'
  long_name='Tetrahedral Cloud Liquid Water Moist Static Energy Detrainment'
  units='K kg/m3/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                  save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  tmp = E_tetra_cloud(1:nx, 1:ny, 1:nzm, 4)
  name='EWTETCLD'
  long_name='Tetrahedral Cloud Vertical Velocity Entrainment'
  units='m/s kg/m3/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                  save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  tmp = D_tetra_cloud(1:nx, 1:ny, 1:nzm, 4)
  name='DWTETCLD'
  long_name='Tetrahedral Cloud Vertical Velocity Detrainment'
  units='kg/m3/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                  save3Dbin,dompi,rank,nsubdomains)

  call task_barrier()

  if(nfields.ne.nfields1) then
    if(masterproc) print*,'write_fields3D error: nfields'
    call task_abort()
  end if
  if(masterproc) then
    close (46)
    if(RUN3D.or.save3Dsep) then
       if(dogzip3D) call system('gzip -f '//filename)
       print*, 'Writting 3D data. file:'//filename
    else
       print*, 'Appending 3D data. file:'//filename
    end if
  endif

end
