!commandline.f90, routine to process command line options
!(C) Roger Wesson
module mod_commandline
use mod_globals
use mod_functions

contains

subroutine readcommandline(runs,switch_ext,switch_he,switch_icf,meanextinction,diagnostics,verbosity,R,identifylines,identifyconfirm,nbins,normalise,norp,calculate_extinction,subtract_recombination,configfile,nperbin)

  implicit none

  character(len=2048), dimension(:), allocatable :: options
  character(len=512) :: configfile
  integer :: i,Narg,runs,verbosity,nbins,nperbin,subtract_recombination
  character :: switch_ext !switch for extinction laws
  character :: switch_he  !switch for helium atomic data
  character :: switch_icf !switch for which ICF scheme to use
  logical :: file_exists,identifylines,identifyconfirm,norp,calculate_extinction
  type(diagnostic_array) :: diagnostics
  real(kind=dp) :: meanextinction, R, normalise

#ifdef CO
        print *,"subroutine: readcommandline"
#endif

  Narg = IARGC() !count input arguments

  if (Narg .eq. 0) then
     print *,"Syntax: neat [option1 value1] [option2 value2] .. [optionx valuex]"
     print *,"type  man neat  for help"
     stop
  endif

  call get_command(commandline)
  print *,gettime(),"command line: ",trim(commandline)

  allocate (options(Narg))

  do i=1,Narg
          call getarg(i,options(i))
  enddo

  do i=1,Narg
          if ((trim(options(i))=="-n" .or. trim(options(i))=="--n-iterations") .and. (i+1) .le. Narg) then
             if (runs .ne. 20000) then
               read (options(i+1),*) runs
             endif
          endif
          if ((trim(options(i))=="-i" .or. trim(options(i))=="--input") .and. (i+1) .le. Narg) then
            filename=trim(options(i+1))
          endif
          if ((trim(options(i))=="-e" .or. trim(options(i))=="--extinction-law") .and. (i+1) .le. Narg) then
            if (trim(options(i+1)) == "LMC")then
              switch_ext = "H"
            elseif (trim(options(i+1)) == "CCM")then
              switch_ext = "C"
            elseif (trim(options(i+1)) == "SMC")then
              switch_ext = "P"
            elseif (trim(options(i+1)) == "Fitz")then
              switch_ext = "F"
            endif
          endif
          if (trim(options(i))=="-c" .and. (i+1) .le. Narg) then
             read (options(i+1),*) meanextinction
             calculate_extinction = .false.
          endif
          if (trim(options(i))=="-nelow" .and. (i+1) .le. Narg) then
             read (options(i+1),*) diagnostics%lowdens
          endif
          if (trim(options(i))=="-nemed" .and. (i+1) .le. Narg) then
             read (options(i+1),*) diagnostics%meddens
          endif
          if (trim(options(i))=="-nehigh" .and. (i+1) .le. Narg) then
             read (options(i+1),*) diagnostics%highdens
          endif
          if (trim(options(i))=="-telow" .and. (i+1) .le. Narg) then
             read (options(i+1),*) diagnostics%lowtemp
          endif
          if (trim(options(i))=="-temed" .and. (i+1) .le. Narg) then
             read (options(i+1),*) diagnostics%medtemp
          endif
          if (trim(options(i))=="-tehigh" .and. (i+1) .le. Narg) then
             read (options(i+1),*) diagnostics%hightemp
          endif
          if ((trim(options(i))=="-he" .or. trim(options(i))=="--helium-data") .and. (i+1) .le. Narg) then
            if (trim(options(i+1))=="S96") then
              switch_he="S"
            endif
          endif
          if (trim(options(i))=="-u" .or. trim(options(i))=="--uncertainties") then
              runs=10000
          endif
          if ((trim(options(i))=="-icf" .or. trim(options(i))=="--ionisation-correction-scheme") .and. (i+1) .le. Narg) then
            if (trim(options(i+1))=="PT92") then
              switch_icf="P"
            elseif (trim(options(i+1))=="KB94") then
              switch_icf="K"
            elseif (trim(options(i+1))=="DI14mod") then
              switch_icf="E"
            endif
          endif
          if ((trim(options(i))=="-v" .or.  trim(options(i))=="--verbosity") .and. (i+1) .le. Narg) then
            read (options(i+1),*) verbosity
            if (verbosity .lt. 1 .or. verbosity .gt. 3) then
              print *,gettime(),"warning: verbosity outside allowed range of 1-3. Set to 3."
              verbosity=3
            endif
          endif
          if (trim(options(i))=="-id" .or. trim(options(i))=="--identify") then
            identifylines=.true.
          endif
          if (trim(options(i))=="-idc" .or. trim(options(i))=="--identify-confirm") then
            identifylines=.true.
            identifyconfirm=.true.
          endif
          if ((trim(options(i))=="-R") .and. (i+1) .le. Narg) then
            read (options(i+1),*) R
          endif
          if ((trim(options(i))=="-nbins") .and. (i+1).le. Narg) then
             read (options(i+1),*) nbins
          endif
          if (trim(options(i))=="-rp") then
             norp=.false.
          endif
          if ((trim(options(i))=="-cf" .or. trim(options(i))=="--configuration-file") .and. (i+1) .le. Narg) then
             configfile=trim(options(i+1))
          endif
          if ((trim(options(i))=="-sr" .or. trim(options(i))=="--subtract-recombination")) then
             subtract_recombination = 2
          endif
          if (trim(options(i))=="--citation") then
             print *
             print *,"NEAT was described in Wesson, Stock and Scicluna, MNRAS, 2012, 422, 3516.  The bibtex data for this paper is:"
             print *
             print *,"@ARTICLE{2012MNRAS.422.3516W,"
             print *,"   author = {{Wesson}, R. and {Stock}, D.~J. and {Scicluna}, P.},"
             print *,"    title = ""{Understanding and reducing statistical uncertainties in nebular abundance determinations}"","
             print *,"  journal = {\mnras},"
             print *,"archivePrefix = ""arXiv"","
             print *,"   eprint = {1203.0567},"
             print *," keywords = {atomic processes, methods: statistical, ISM: abundances},"
             print *,"     year = 2012,"
             print *,"    month = jun,"
             print *,"   volume = 422,"
             print *,"    pages = {3516-3526},"
             print *,"      doi = {10.1111/j.1365-2966.2012.20863.x},"
             print *,"   adsurl = {http://adsabs.harvard.edu/abs/2012MNRAS.422.3516W},"
             print *,"  adsnote = {Provided by the SAO/NASA Astrophysics Data System}"
             print *,"}"
             call exit(0)
          endif
          if ((trim(options(i))=="-o" .or. trim(options(i))=="--output-dir") .and. (i+1) .le. Narg) then
            outputdirectory=trim(options(i+1))
          endif
          if ((trim(options(i))=="-of" .or. trim(options(i))=="--output-format") .and. (i+1) .le. Narg) then
            outputformat=trim(options(i+1))
            if (outputformat.ne."fits" .and. outputformat.ne."text") then
              print *,gettime(),"invalid output format. valid options are 'text' and 'fits'"
              call exit(1)
            endif
          endif
  !  to be fully implemented:
  !  -R                     : R (default 3.1) - only used with CCM at the moment
  !  -b                     : batch mode
   enddo

   if (Narg .eq. 1) then
     filename=trim(options(1))
   endif

   if (filename=="") then
          print *,gettime(),"error: No input file specified"
          call exit(1)
   endif

   inquire(file=filename, exist=file_exists) ! see if the input file is present

   if (.not. file_exists) then
          print *,gettime(),"error: input file ",trim(filename)," does not exist"
          call exit(1)
   endif

   if (configfile .ne. "") then
     inquire(file=configfile, exist=file_exists) ! see if the configuration file is present
   endif

   if (.not. file_exists) then
          print *,gettime(),"error: configuration file ",trim(configfile)," does not exist"
          call exit(1)
   endif

!check number of runs

  if (runs .gt. 1 .and. runs .lt. 5000) print *,gettime(),"warning: number of iterations is low.  At least 5000 is recommended for good sampling of probability distributions"
  if (runs .gt. 1) nperbin=runs/nbins
  if (mod(runs,nbins) .gt. 0 .and. runs .gt. 1) then
     print*,gettime(),': error: number of iterations does not divide exactly by number of bins'
     print*,'            Please set number of iterations to an exact multiple of ',nbins
     print*,'            or modify the number of bins using -nbins'
     call exit(1)
  endif

  deallocate(options)

end subroutine readcommandline

end module mod_commandline
