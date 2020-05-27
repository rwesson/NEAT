!globals.f90, a few things to be available everywhere
!(C) Roger Wesson
module mod_globals

  implicit none

! precision

  integer,parameter :: sp = kind(1.0)
  integer,parameter :: dp = kind(1.d0)

  integer,parameter :: spc = kind( (1.0_sp,1.0_sp) )
  integer,parameter :: dpc = kind( (1.0_dp,1.0_dp) )

! command line and file names

  character(len=2048) :: commandline
  character(len=512) :: filename,configfile,defaultconfigfile

end module mod_globals
