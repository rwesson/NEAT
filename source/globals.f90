!globals.f90, a few things to be available everywhere
!(C) Roger Wesson
module mod_globals

  implicit none

! precision

  integer,parameter :: sp = kind(1.0)
  integer,parameter :: dp = kind(1.d0)

  integer,parameter :: spc = kind( (1.0_sp,1.0_sp) )
  integer,parameter :: dpc = kind( (1.0_dp,1.0_dp) )

end module mod_globals
