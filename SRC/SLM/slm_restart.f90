!==================================================================================
! Based on restart.f90 of SAM
! Contians write and read restart file subroutines for a certain SLM variables
!==================================================================================
	subroutine write_statement_slm()
	! Write slm variables in slm restart file
	use slm_vars
        use rad, only: coszrsxy, swdsvisxy, swdsnirxy, swdsvisdxy, swdsnirdxy, lwdsxy

	implicit none
	character *4 rankchar
        character *256 filename
	integer irank
        integer lenstr
        external lenstr

        if(masterproc) print*,'Writting SLM restart file...'

        if(restart_sep) then

          write(rankchar,'(i4)') rank

          filename='./RESTART/'//case(1:len_trim(case))//'_'// &
               caseid(1:len_trim(caseid))//'_'// &
               rankchar(5-lenstr(rankchar):4)//'_restart_slm.bin'
          open(176,file=trim(filename),status='unknown',form='unformatted', &
!              BUFFERED='YES', & ! use for intel compiler
              ACTION='WRITE')
          write(176) nsubdomains
	  write(176) soilt, soilw, mw, mws, t_canop, t_cas, q_cas, soilt_obs, soilw_obs,&
		dosoilwnudging, dosoiltnudging,tausoil,landtype0,LAI0,clay0,sand0, &
		ustar, tstar, nsoil
          close(176)

        else  ! not restart_sep

          write(rankchar,'(i4)') nsubdomains

          do irank=0,nsubdomains-1

             call task_barrier()

             if(irank.eq.rank) then

               if(masterproc) then

                  open(176,file='./RESTART/'//case(1:len_trim(case))//'_'// &
                      caseid(1:len_trim(caseid))//'_'// &
                      rankchar(5-lenstr(rankchar):4)//'_restart_slm.bin', &
                      status='unknown',form='unformatted', &
!                      BUFFERED='YES',  & ! use for intel compiler
                      ACTION='WRITE')
                  write(176) nsubdomains

               else

                  open(176,file='./RESTART/'//case(1:len_trim(case))//'_'// & 
                      caseid(1:len_trim(caseid))//'_'// &
                      rankchar(5-lenstr(rankchar):4)//'_restart_slm.bin', &
                      status='unknown',form='unformatted', position='append', &
!                      BUFFERED='YES', & ! use for intel compiler
                      ACTION='WRITE')

               end if
              write(176) soilt, soilw, mw, mws, t_canop, t_cas, q_cas, soilt_obs, soilw_obs,&
                dosoilwnudging, dosoiltnudging,tausoil,landtype0,LAI0,clay0,sand0, &
                ustar, tstar, nsoil
               close(176)
             end if
          end do

        end if ! restart_sep

	if(masterproc) then
           print *,'Saved SLM restart file. nstep=',nstep
	endif

        call task_barrier()

        return
        end
 
 
 
 
     
	subroutine read_statement_slm()
	
	use slm_vars
        use rad, only: coszrsxy, swdsvisxy, swdsnirxy, swdsvisdxy, swdsnirdxy, lwdsxy

	implicit none
	character *4 rankchar
        character *256 filename
	integer irank,ii,nsoil1
        integer lenstr
        external lenstr
	
        if(masterproc) print*,'Reading SLM restart file...'

        if(restart_sep) then

          write(rankchar,'(i4)') rank

          if(nrestart.ne.2) then
                filename ='./RESTART/'//case(1:len_trim(case))//'_'// &
                        caseid(1:len_trim(caseid))//'_'// &
                         rankchar(5-lenstr(rankchar):4)//'_restart_slm.bin'
          else
                filename ='./RESTART/'//case_restart(1:len_trim(case_restart))//'_'// &
                        caseid_restart(1:len_trim(caseid_restart))//'_'// &
                         rankchar(5-lenstr(rankchar):4)//'_restart_slm.bin'
           end if

          open(176,file=filename, status='unknown',form='unformatted', &
!                   BUFFERED='YES', & ! use for intel compiler
                   ACTION='READ')
          read (176)
          read (176) soilt, soilw, mw, mws, t_canop, t_cas, q_cas, soilt_obs, soilw_obs,&
                dosoilwnudging, dosoiltnudging,tausoil,landtype0,LAI0,clay0,sand0, &
                ustar, tstar, nsoil1
          close(176)

        else

          write(rankchar,'(i4)') nsubdomains
          if(nrestart.ne.2) then
                filename ='./RESTART/'//case(1:len_trim(case))//'_'// &
                        caseid(1:len_trim(caseid))//'_'// &
                         rankchar(5-lenstr(rankchar):4)//'_restart_slm.bin'
          else
                filename ='./RESTART/'//case_restart(1:len_trim(case_restart))//'_'// &
                        caseid_restart(1:len_trim(caseid_restart))//'_'// &
                         rankchar(5-lenstr(rankchar):4)//'_restart_slm.bin'
           end if

          open(176,file=filename,status='unknown',form='unformatted', &
!                   BUFFERED='YES', & ! use for intel compiler
                   ACTION='READ')

          do irank=0,nsubdomains-1

             call task_barrier()

             if(irank.eq.rank) then

               read (176)

               do ii=0,irank-1 ! skip records
                  read (176)
               end do
               read (176) soilt, soilw, mw, mws, t_canop, t_cas, q_cas, soilt_obs, soilw_obs,&
                dosoilwnudging, dosoiltnudging,tausoil,landtype0,LAI0,clay0,sand0, &
                ustar, tstar, nsoil1
               close(176)
             end if

          end do

        end if ! restart_sep

        if(nsoil1.ne.nsoil) then
           if(masterproc)print*,'number of soil levels is wrong for the restart.'
           if(masterproc)print*,'Required nsoil:',nsoil,' Quitting...'
           call task_abort()
        end if
        if(rank.eq.nsubdomains-1) then
             if(masterproc)print *,'Case:',caseid
             if(masterproc)print *,'Restart SLM at step:',nstep
             if(masterproc)print *,'Time:',nstep*dt
        endif


        call task_barrier()


        return
        end
