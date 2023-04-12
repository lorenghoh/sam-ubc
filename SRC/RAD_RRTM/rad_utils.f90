module rad_utils
  implicit none

 contains
  !----------------------------------------------------------------------
  function rad_scheme_name()
    character(len=32) :: rad_scheme_name
    ! Return the scheme name, normally the same as the directory name with leading "RAD_" removed  
    rad_scheme_name = "rrtm" 
  end function   rad_scheme_name
  ! ----------------------------------------------------------------------------
end module rad_utils
