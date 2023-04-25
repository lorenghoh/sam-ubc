! Wrapper call to avoid circular dependency between hbuffer and entrainment module

! UBC ENT
subroutine hbuf_ent_init(namelist,deflist,unitlist,status,average_type,count,trcount)
  use entrainment
  implicit none
  character(*) namelist(*), deflist(*), unitlist(*)
  integer status(*), average_type(*), count, trcount

  call ent_hbuf_init(namelist,deflist,unitlist,status,average_type,count,trcount)
end
! End UBC ENT
