! NEAT, the nebular abundance analysis tool
! (C) 2006-2012 Roger Wesson, Dave Stock, Peter Scicluna
! NEAT incorporates aspects of several codes developed over decades at UCL:
! equib, by I. Howarth and S. Adams, updated significantly by B. Ercolano
! MIDAS scripts written by X-W. Liu for calculating recombination line abundances
! It also uses the quicksort algorithm as implemented in F90 by Alberto Ramos

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

program neat

        use mod_abundtypes
        use mod_resultarrays
        use mod_extinction
        use mod_quicksort
        use mod_abundIO
        use mod_atomicdata
        use mod_atomic_read
        use mod_recombination_lines
        use mod_linefinder

        IMPLICIT NONE

        CHARACTER :: switch_ext !switch for extinction laws
        CHARACTER :: switch_he  !switch for helium atomic data
        CHARACTER :: switch_icf !switch for which ICF scheme to use
        INTEGER :: I, j, runs, Narg !runs = number of runs for randomiser

!input options

        CHARACTER(len=2048), DIMENSION(:), allocatable :: options
        CHARACTER(len=2048) :: commandline

!file reading variables

        LOGICAL :: file_exists
        LOGICAL :: identifylines
        TYPE(LINE),DIMENSION(:), allocatable :: linelist
        TYPE(LINE),DIMENSION(:), allocatable :: linelist_original
        TYPE(LINE),dimension(:,:), allocatable :: all_linelists
        CHARACTER(len=512) :: filename
        CHARACTER(len=1) :: blank
        INTEGER :: listlength, errstat
        real(kind=dp) :: normalise

!results and result processing

        type(resultarray), dimension(:), allocatable :: all_results
        type(resultarray), dimension(1) :: iteration_result
        real(kind=dp), dimension(:), allocatable :: quantity_result
        real(kind=dp), dimension(:,:), allocatable :: resultprocessingarray
        character(len=40), dimension(:,:), allocatable :: resultprocessingtext

!atomic data

        character(len=10) :: ionlist(40) !list of ion names
        integer :: iion !# of ions in Ilines
        integer :: Iint 
        integer :: maxlevs,maxtemps
        type(atomic_data),allocatable :: atomicdata(:)
        real(kind=dp), dimension(:,:,:), allocatable :: heidata

!extinction

        logical :: calculate_extinction=.true.
        real(kind=dp) :: meanextinction, R

!RP effect switch

        logical :: norp=.false.

!binning and uncertainties

        real(kind=dp), dimension(3) :: uncertainty_array=0d0
        real(kind=dp), dimension (:,:), allocatable :: binned_quantity_result
        logical :: unusual
        integer :: verbosity,nbins,nperbin

!CEL array

        TYPE(line), DIMENSION(:), allocatable :: ILs

!diagnostic array

        real(kind=dp), dimension(6) :: diagnostic_array

!output formats

        character(len=35) :: extinction_format, diagnostic_format, diagnostic_ratio_format, abundances_format, adf_format

!multiple output formats defined as variables so they can be passed to
!the printout subroutine

        extinction_format = "(X,A,F9.3,SP,F9.3,F9.3,S)"
        diagnostic_format = "(X,A,F9.1,SP,F9.1,F9.1,S)"
        diagnostic_ratio_format = "(X,A,F9.3,SP,F9.3,F9.3,S)"
        abundances_format = "(X,A,ES14.3,SP,ES14.3,ES14.3,S)"
        adf_format = "(X,A,F8.2,SP,F8.2,F8.2,S)"

!read command line arguments

        Narg = IARGC() !count input arguments

        if (Narg .eq. 0) then
           print *,"Syntax: ./neat.exe [option1 value1] [option2 value2] .. [optionx valuex]"
           print *,"Available options:"
           print *,"  -i / --input"
           print *,"       Input file"
           print *,"       No default"
           print *,"  -u / --uncertainties"
           print *,"       Calculate uncertainties, using 20,000 iterations"
           print *,"       Default: no calculation of uncertainties"
           print *,"       If this option is specified, the -n option will be ignored"
           print *,"  -n / --n-iterations"
           print *,"       Number of iterations"
           print *,"       Default: 1"
           print *,"  -e / --extinction-law"
           print *,"       Extinction law"
           print *,"       Values:"
           print *,"          How:  Galactic law of Howarth (1983, MNRAS, 203, 301)"
           print *,"          CCM:  Galactic law of Cardelli, Clayton, Mathis (1989, ApJ, 345, 245)"
           print *,"          Fitz: Galactic law of Fitzpatrick & Massa (1990, ApJS, 72, 163)"
           print *,"          LMC:  LMC law of Howarth (1983, MNRAS, 203, 301)"
           print *,"          SMC:  SMC law of Prevot et al. (984, A&A, 132, 389)"
           print *,"       Default: How"
           print *,"  -c"
           print *,"       The logarithmic extinction at H beta"
           print *,"       Default: calculated from Balmer line ratios"
           print *,"  -nelow, -nemed, -nehigh"
           print *,"  -telow, -temed, -tehigh"
           print *,"       The electron densities and temperatures to be used."
           print *,"       Units: cm-3 for densities, K for temperatures"
           print *,"       Default: calculated from available diagnostics."
           print *,"  -he / --helium-data"
           print *,"       The atomic data to use for He I abundances"
           print *,"       Values:"
           print *,"          S96: Smits, 1996, MNRAS, 278, 683"
           print *,"          P12: Porter et al., 2012, MNRAS, 425, 28"
           print *,"       Default: S96"
           print *,"  -icf / --ionisation-correction-scheme"
           print *,"       The ICF scheme to be used to correct for unseen ions"
           print *,"       Values:"
           print *,"          KB94: Kingsburgh & Barlow (1994, MNRAS, 271, 257)"
           print *,"          PT92: Peimbert, Torres-Peimbert & Ruiz (1992, RMxAA, 24, 155)"
           print *,"          DI14: Delgado-Inglada, Morisset & Stasinska (2014, MNRAS, 440, 536)"
           print *,"       Default: KB94"
           print *,"  -v / --verbosity"
           print *,"       Amount of output to write for each derived quantity"
           print *,"       This option has no effect unless -n or -u is specified"
           print *,"       Values:"
           print *,"         1: write out summary files, binned results and complete results"
           print *,"         2: write out summary files and binned results"
           print *,"         3: write out summary files only"
           print *,"       Default: 3"
           print *,"  -id / --identify"
           print *,"       This option triggers an algorithm which attempts to convert observed wavelengths into rest wavelengths"
           print *,"       If the option is not present, NEAT assumes that the input line list already has rest wavelengths in the first column"
           print *,"       No default"
           print *,"  -norp"
           print *,"       When calculating Monte Carlo uncertainties, NEAT's default behaviour is to compensate for the upward bias"
           print *,"       affecting weak lines described by Rola and Pelat (1994).  This option overrides that and no compensation is calculated."
        !  to be fully implemented:
        !  -R                     : R (default 3.1)
        !  --nbins                : user can set number of bins for routine
        !  -b                     : batch mode
           stop
        endif

        call get_command(commandline)

        ALLOCATE (options(Narg))

        do i=1,Narg
                call getarg(i,options(i))
        enddo

! set defaults

        runs=1
        switch_ext="S" !Howarth 1983 Galactic law
        switch_he="S"  !Smits 1996 data
        switch_icf="K" !KB94 ICF
        filename=""
        meanextinction=0.D0
        diagnostic_array=0.D0
        verbosity=3
        R=3.1
        identifylines=.false.
        nbins=25
        normalise = 0.d0

! start the logging output to terminal

        print *,"NEAT, the Nebular Empirical Analysis Tool"
        print *,"-----------------------------------------"

        print *
        print *,gettime(),": starting code"
        print *,gettime(),": command line: ",trim(commandline)

! process command line arguments

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
                   read (options(i+1),*) diagnostic_array(1)
                endif
                if (trim(options(i))=="-nemed" .and. (i+1) .le. Narg) then
                   read (options(i+1),*) diagnostic_array(2)
                endif
                if (trim(options(i))=="-nehigh" .and. (i+1) .le. Narg) then
                   read (options(i+1),*) diagnostic_array(3)
                endif
                if (trim(options(i))=="-telow" .and. (i+1) .le. Narg) then
                   read (options(i+1),*) diagnostic_array(4)
                endif
                if (trim(options(i))=="-temed" .and. (i+1) .le. Narg) then
                   read (options(i+1),*) diagnostic_array(5)
                endif
                if (trim(options(i))=="-tehigh" .and. (i+1) .le. Narg) then
                   read (options(i+1),*) diagnostic_array(6)
                endif
                if ((trim(options(i))=="-he" .or. trim(options(i))=="--helium-data") .and. (i+1) .le. Narg) then
                  if (trim(options(i+1))=="P12") then
                    switch_he="P"
                  endif
                endif
                if (trim(options(i))=="-u" .or. trim(options(i))=="--uncertainties") then
                    runs=10000
                endif
                if ((trim(options(i))=="-icf" .or. trim(options(i))=="--ionisation-correction-scheme") .and. (i+1) .le. Narg) then
                  if (trim(options(i+1))=="PT92") then
                    switch_icf="P"
                  elseif (trim(options(i+1))=="DI14") then
                    switch_icf="D"
                  endif
                endif
                if ((trim(options(i))=="-v" .or.  trim(options(i))=="--verbosity") .and. (i+1) .le. Narg) then
                  read (options(i+1),*) verbosity
                  if (verbosity .lt. 1 .or. verbosity .gt. 3) then
                    print *,gettime(),": warning: verbosity outside allowed range of 1-3. Set to 3."
                    verbosity=3
                  endif
                endif
                if ((trim(options(i))=="-id" .or. trim(options(i))=="--identify") .and. (filename.ne."")) then
                  identifylines=.true.
                endif
                if ((trim(options(i))=="-R") .and. (i+1) .le. Narg) then
                  read (options(i+1),*) R
                endif
                if ((trim(options(i))=="-nbins") .and. (i+1).le. Narg) then
                   read (options(i+1),*) nbins
                endif
                if (trim(options(i))=="-norp") then
                   norp=.true.
                endif
         enddo

         if (Narg .eq. 1) then
           filename=trim(options(1))
         endif

         if (filename=="") then
                print *,gettime(),": error: No input file specified"
                stop
         endif

         inquire(file=filename, exist=file_exists) ! see if the input file is present

         if (.not. file_exists) then
                print *,gettime(),": error: input file ",trim(filename)," does not exist"
                stop
         endif

!check number of runs

        if (runs .gt. 1 .and. runs .lt. 5000) print *,gettime(),": warning: number of iterations is low.  At least 5000 is recommended for good sampling of probability distributions"
        if (runs .gt. 1) nperbin=runs/nbins
        if (mod(runs,nbins) .gt. 0 .and. runs .gt. 1) then
           print*,gettime(),': error: number of iterations does not divide exactly by number of bins'
           print*,'            Please set number of iterations to an exact multiple of ',nbins
           print*,'            or modify the number of bins using -nbins'
           stop
        endif

        deallocate(options)

! read in the line list, allocate original linelist array

        call read_linelist(filename,linelist,listlength, errstat)

        if (errstat .eq. 1) then
                print *,gettime(),": error: line list reading failed"
                print *,"            This can happen if it doesn't have three columns"
                stop
        elseif (errstat .eq. 2) then
                PRINT*, gettime(),": cheese shop error - no inputs"
                STOP
        else
                print "(X,A,A,A,A,I3,A)", gettime(),": line list file ",trim(filename)," read successfully (",listlength," lines)"
        endif

        allocate (linelist_original(listlength))

! run line identifier if required

        if (identifylines) then
                print *,gettime()," : running line finder"
                print *,"---------------------------------"
                call linefinder(linelist, listlength)
                print *
                print *,gettime()," : line finder finished"
                print *,gettime()," : WARNING!!!  The line finding algorithm is intended as an aid only and is not designed to be highly robust"
                print *,gettime()," : check your line list very carefully for potentially wrongly identified lines!!"
                print *,"---------------------------------"
                print *,"Are these line IDs ok? (y/n)"
                read (5,*) blank
                if (blank .eq. "n" .or. blank .eq. "N") then
                  print *,gettime()," : analysis cancelled."
                  stop
                endif
        endif

! normalise to Hbeta

        if (minval(abs(linelist(:)%wavelength - 4861.33)) .lt. 0.001) then
          normalise = 100.d0/linelist(minloc(abs(linelist(:)%wavelength - 4861.33),1))%intensity
        else
          print *,gettime(),": error: no H beta detected"
          print *,"            no further analysis possible"
          stop
        endif

        linelist%intensity = linelist%intensity * normalise
        linelist%int_err = linelist%int_err * normalise

! check for and remove negative line fluxes and uncertainties

        if (minval(linelist%intensity) .lt. 0.) then
          where (linelist%intensity .lt. 0)
            linelist%intensity = 0.D0
          endwhere
        print *,gettime(),": warning: negative line fluxes set to zero"
        endif

        if (minval(linelist%int_err) .lt. 0.) then
          where (linelist%int_err .lt. 0)
            linelist%int_err = abs(linelist%int_err)
          endwhere
        print *,gettime(),": warning: negative flux uncertainties converted to positive"
        endif

        print *

        if (switch_ext == "S") then
                print *,gettime(), ": using Howarth (1983) galactic extinction law"
        elseif (switch_ext == "H") then
                print *,gettime(), ": using Howarth (1983) LMC extinction law"
        elseif (switch_ext == "C") then
                print *,gettime(), ": using CCM (1989) galactic extinction law"
        elseif (switch_ext == "P") then
                print *,gettime(), ": using Prevot et al. (1984) SMC extinction law"
        elseif (switch_ext == "F") then
                print *,gettime(), ": using Fitzpatrick (1990) galactic extinction law"
        endif

        if (switch_he == "S") then
                print *,gettime(), ": using Smits (1996) He I emissivities"
        else
                print *,gettime(), ": using Porter et al. (2012) He I emissivities"
        endif

        if (switch_ICF == "K") then
                print *,gettime(), ": using Kingsburgh & Barlow (1994) ICF"
        elseif (switch_ICF == "D") then
                print *,gettime(), ": using Delgado-Inglada et al. (2014) ICF"
        else
                print *,gettime(), ": using Peimbert et al. (1992) ICF"
        endif

        print *

! read the CEL data

        call read_ilines(ILs, Iint,iion,ionlist)

        print *,gettime(), ": reading atomic data"

!CELs read, can now read in atomic data

        allocate(atomicdata(iion))
        DO I=1,iion
            atomicdata(I)%ion = ionlist(I)
            call read_atomic_data(atomicdata(I))
        ENDDO

!read ORL data
        call read_orl_data

!read helium emissivities

        if (switch_he .eq. "P") then
          allocate(heidata(21,15,44))
          heidata = 0.D0
          call read_porter(heidata)
        elseif (switch_he .eq. "S") then
          allocate(heidata(3,6,44))
          heidata = 0.D0
          call read_smits(heidata)
        endif

!find maximum #levels and temperatures - pass to equib to reduce footprint

        maxlevs = atomicdata(1)%nlevs
        maxtemps = atomicdata(1)%ntemps
        do i=2,iion
            if(atomicdata(i)%nlevs .gt. maxlevs) maxlevs = atomicdata(i)%nlevs
            if(atomicdata(i)%ntemps .gt. maxtemps) maxtemps = atomicdata(i)%ntemps
        enddo
!        print*,maxlevs,maxtemps
!now check number of iterations.  If 1, line list is fine as is.  If more than one, randomize the fluxes

        allocate(all_results(runs))

        if(runs == 1)then !calculates abundances without uncertainties
                call abundances(linelist, switch_ext, listlength, iteration_result, R, meanextinction, calculate_extinction, ILs, Iint, diagnostic_array,iion,atomicdata,maxlevs,maxtemps, heidata, switch_he, switch_icf)
                all_results(1)=iteration_result(1) ! copy the iteration result to all_results to simplify the writing out of results later

        else if(runs > 1)then

!save unrandomised line list

                linelist_original = linelist

                call init_random_seed()!sets seed for randomiser

                !allocate array to store all line info

                allocate(all_linelists(size(linelist),runs))

                !main loop

                print *
                print "(X,A9,X,A)",gettime(), ": starting Monte Carlo calculations"
                print *,gettime(), ": completed ",0,"%"
                DO I=1,runs
                         if ( (10.0*dble(i)/dble(runs)) == int(10*i/runs) ) print *,gettime(),": completed ",100*i/runs,"%"
!                        print*, "iteration ", i, "of", runs

                        call randomizer(linelist, listlength, norp)
                        call abundances(linelist, switch_ext, listlength, iteration_result, R, meanextinction, calculate_extinction, ILs, Iint, diagnostic_array,iion,atomicdata,maxlevs,maxtemps, heidata, switch_he, switch_icf)

                        !store all line and derived quantity in arrays
                        all_linelists(:,i)=linelist 
                        all_results(i)=iteration_result(1)
                        !restore the unrandomized line list ready for the next iteration
                        linelist = linelist_original

                END DO

        else
                print*, gettime(), ": error: I didn't want to be a barber anyway. I wanted to be... a lumberjack!   Also, a positive number of runs helps.."
                stop
        endif

!now write all the lines to line list files, plain text and latex
!linelist first

        print *,gettime(),": Writing line list"

        allocate(quantity_result(runs))
        quantity_result=0d0
        open (650,FILE=trim(filename)//"_linelist", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
        open (651,FILE=trim(filename)//"_linelist.tex", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')

        write (650,*) "Lambda  Ion          F(line)             I(line)               Abundance"
        write (651,*) "\begin{longtable}{lrlrlllllll}"
        write (651,*) "\hline"
        write (651,*) "$ \lambda $ & Ion & $F \left( \lambda \right) $ && $I \left( \lambda \right) $ & Ion & Multiplet & Lower term & Upper term & g$_1$ & g$_2$ \\"
        write (651,*) "\hline"

        if (runs .gt. 1) then
                do j=1, listlength

!rest wavelength, ion name for plain text file

                write (650,"(X,F7.2,X,A11)", advance='no') all_linelists(j,1)%wavelength,all_linelists(j,1)%name
                write (651,"(X,F7.2,' & ',A15,' & ')", advance='no') all_linelists(j,1)%wavelength

!line flux

               if (linelist_original(j)%intensity .eq. 0.d0) then
                 write (650,"(A)", advance='no') " * "
                 write (651,"(A)", advance='no') " *     &             &"
               else
                 write (650,"(F7.3,A,F7.2,3X)", advance='no') linelist_original(j)%intensity," +-",linelist_original(j)%int_err
                 write (651,"(F7.3,'& $\pm$',F7.2, '&')", advance='no') linelist_original(j)%intensity,linelist_original(j)%int_err
               endif

!dereddened flux
                quantity_result = all_linelists(j,:)%int_dered

                call get_uncertainties(quantity_result, binned_quantity_result, uncertainty_array, unusual,nbins,nperbin,runs)

                if (linelist_original(j)%intensity .ne. 0.d0) then
                  if (uncertainty_array(1) .ne. uncertainty_array(3)) then
                    write (650,"(F7.3,SP,F7.2,SP,F7.2)", advance='no') uncertainty_array(2),uncertainty_array(1),-uncertainty_array(3)
                    write (651,"(F7.3,'& $^{',SP,F7.2,'}_{',SP,F7.2,'}$')", advance='no') uncertainty_array(2),uncertainty_array(1),-uncertainty_array(3)
                  else
                    write (650,"(F7.3,A,F7.2,4X)", advance='no') uncertainty_array(2)," +-",uncertainty_array(1)
                    write (651,"(F7.3,'& $\pm$',F7.2)", advance='no') uncertainty_array(2),uncertainty_array(1)
                  endif
                else
                  write (650,"(A)", advance='no') " * "
                  write (651,"(A)", advance='no') " *     &        "
                endif

! transition data

!                write (650,*)
                write (651,*) linelist_original(j)%linedata, "\\"

!abundance - write out if there is an abundance for the line, don't write
!anything except a line break if there is no abundance for the line.
!todo: add an option to choose whether or not to put abundances in the line list table

                quantity_result = all_linelists(j,:)%abundance
                call get_uncertainties(quantity_result, binned_quantity_result, uncertainty_array, unusual,nbins,nperbin,runs)
                if (uncertainty_array(2) .ne. 0.D0) then
                  if (uncertainty_array(1) .ne. uncertainty_array(3)) then
                    write (650,"(ES10.2,SP,ES10.2,SP,ES10.2)") uncertainty_array(2),uncertainty_array(1),-uncertainty_array(3)
!                    write (651,"(' & ${',A,'}$ & $^{+',A,'}_{',A,'}$ \\')") trim(latex_number(uncertainty_array(2))),trim(latex_number(uncertainty_array(1))),trim(latex_number(-uncertainty_array(3)))
                  else
                    write (650,"(ES10.2,A,ES10.2)") uncertainty_array(2)," +-",uncertainty_array(1)
!                    write (651,"(' & $',A,'$ & $\pm',A,'$\\')") trim(latex_number(uncertainty_array(2))),trim(latex_number(uncertainty_array(1)))
                  endif
                else
                  write (650,*)
!                  write (651,*) "\\"
                endif

                end do
        else ! runs == 1, no uncertainties to write out

                do i=1,listlength
                  if (linelist(i)%intensity .ne. 0.d0) then
                    if (linelist(i)%abundance .gt. 0.0) then
                      write (650,"(X,F7.2,X,A11,F7.3,X,F7.3,X,ES14.3)") linelist(i)%wavelength,linelist(i)%name,linelist(i)%intensity,linelist(i)%int_dered, linelist(i)%abundance
                      write (651,"(X,F7.2,X,'&',X,F7.3,X,'&',X,F7.3,X,A,'\\')") linelist(i)%wavelength,linelist(i)%intensity,linelist(i)%int_dered,linelist(i)%linedata
                    else
                      write (650,"(X,F7.2,X,A11,F7.3,X,F7.3)") linelist(i)%wavelength,linelist(i)%name,linelist(i)%intensity,linelist(i)%int_dered
                      write (651,"(X,F7.2,X,'&',X,F7.3,X,'&',X,F7.3,X,A,'\\')") linelist(i)%wavelength,linelist(i)%intensity,linelist(i)%int_dered,linelist(i)%linedata
                    endif
                  else
                    write (650,"(X,F7.2,X,A11,A7,X,A7)") linelist(i)%wavelength,linelist(i)%name,"*      ","*      "
                    write (651,"(X,F7.2,X,'&',X,A7,X,'&',X,A7,X,A,'\\')") linelist(i)%wavelength,"*      ","*      ",linelist(i)%linedata
                  endif
                end do

        endif

        write (651,*) "\hline"
        !write (651,*) "\caption{}"
        write (651,*) "\label{tab:",trim(filename)//"_linelist}"
        write (651,*) "\end{longtable}"

        close(650)
        close(651)

        print *,gettime(),": Linelist written to files:"
        print *,"              ",trim(filename),"_linelist"
        print *,"              ",trim(filename),"_linelist.tex"

!now write out the summary files and all the binned data

        print *
        print *,gettime(),": Writing summary files"

!first, define arrays with links to all the data that needs processing.
!extinction, diagnostics, cel abundances, orl abundances, strong line
!abundances, adfs

        allocate(resultprocessingarray(155,runs))
        resultprocessingarray=0d0
        allocate(resultprocessingtext(155,4))

!extinction

        resultprocessingarray(1,:) = all_results%mean_cHb
        resultprocessingtext(1,:) = (/"c(Hb)                              ","c(H$\beta)$:                       ", extinction_format, "mean_chb                           "/)

!diagnostics
!low density

        resultprocessingarray(2,:) = all_results%oii_density
        resultprocessingtext(2,:) = (/"[OII] density                      ","{}[O~{\sc ii}] density:            ", diagnostic_format, "density_oii                        "/)
        resultprocessingarray(3,:) = all_results%oii_density_ratio
        resultprocessingtext(3,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "density_oii_ratio                  "/)
        resultprocessingarray(4,:) = all_results%SII_density
        resultprocessingtext(4,:) = (/"[SII] density                      ","{}[S~{\sc ii}] density:            ", diagnostic_format, "density_sii                        "/)
        resultprocessingarray(5,:) = all_results%sii_density_ratio
        resultprocessingtext(5,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "density_sii_ratio                  "/)
        resultprocessingarray(6,:) = all_results%low_density
        resultprocessingtext(6,:) = (/"Low ionisation density             ","Low ionisation density:            ", diagnostic_format, "density_low                        "/)

!low temperature

        resultprocessingarray(7,:) = all_results%oii_temp
        resultprocessingtext(7,:) = (/"[OII] temperature                  ","{}[O~{\sc ii}] temperature         ", diagnostic_format, "temp_oii                           "/)
        resultprocessingarray(8,:) = all_results%oii_temp_ratio
        resultprocessingtext(8,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "temp_oii_ratio                     "/)
        resultprocessingarray(9,:) = all_results%SII_temp
        resultprocessingtext(9,:) = (/"[SII] temperature                  ","{}[S~{\sc ii}] temperature         ", diagnostic_format, "temp_sii                           "/)
        resultprocessingarray(10,:) = all_results%sii_temp_ratio
        resultprocessingtext(10,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "temp_sii_ratio                     "/)
        resultprocessingarray(11,:) = all_results%NII_temp
        resultprocessingtext(11,:) = (/"[NII] temperature                  ","{}[N~{\sc ii}] temperature         ", diagnostic_format, "temp_nii                           "/)
        resultprocessingarray(12,:) = all_results%nii_temp_ratio
        resultprocessingtext(12,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "temp_nii_ratio                     "/)
        resultprocessingarray(13,:) = all_results%OI_temp
        resultprocessingtext(13,:) = (/"[OI] temperature                   ","{}[O~{\sc i}] temperature          ", diagnostic_format, "temp_oi                            "/)
        resultprocessingarray(14,:) = all_results%oi_temp_ratio
        resultprocessingtext(14,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "temp_oi_ratio                      "/)
        resultprocessingarray(15,:) = all_results%CI_temp
        resultprocessingtext(15,:) = (/"[CI] temperature                   ","{}[C~{\sc i}] temperature          ", diagnostic_format, "temp_ci                            "/)
        resultprocessingarray(16,:) = all_results%ci_temp_ratio
        resultprocessingtext(16,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "temp_ci_ratio                      "/)
        resultprocessingarray(17,:) = all_results%low_temp
        resultprocessingtext(17,:) = (/"Low ionisation temperature         ","Low ionisation temperature         ", diagnostic_format, "temp_low                           "/)


!medium density

        resultprocessingarray(18,:) = all_results%cliii_density
        resultprocessingtext(18,:) = (/"[ClIII] density                    ","{}[Cl~{\sc iii}] density:          ", diagnostic_format, "density_cliii                      "/)
        resultprocessingarray(19,:) = all_results%cliii_density_ratio
        resultprocessingtext(19,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "density_cliii_ratio                "/)
        resultprocessingarray(20,:) = all_results%ArIV_density
        resultprocessingtext(20,:) = (/"[ArIV] density                     ","{}[Ar~{\sc iv}] density:           ", diagnostic_format, "density_ariv                       "/)
        resultprocessingarray(21,:) = all_results%ariv_density_ratio
        resultprocessingtext(21,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "density_sii_ratio                  "/)
        resultprocessingarray(22,:) = all_results%CIII_density
        resultprocessingtext(22,:) = (/"[CIII] density                     ","{}[C~{\sc iii}] density            ", diagnostic_format, "density_ciii                       "/)
        resultprocessingarray(23,:) = all_results%ciii_density_ratio
        resultprocessingtext(23,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "density_ciii_ratio                 "/)
        resultprocessingarray(24,:) = all_results%OIII_IR_density
        resultprocessingtext(24,:) = (/"[OIII] IR density                  ","{}[O~{\sc iii}] IR density         ", diagnostic_format, "density_oiii_ir                    "/)
        resultprocessingarray(25,:) = all_results%oiii_ir_density_ratio
        resultprocessingtext(25,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "density_oiii_ir_ratio              "/)
        resultprocessingarray(26,:) = all_results%SIII_IR_density
        resultprocessingtext(26,:) = (/"[SIII] IR density                  ","{}[S~{\sc iii}] IR density         ", diagnostic_format, "density_siii_ir                    "/)
        resultprocessingarray(27,:) = all_results%siii_ir_density_ratio
        resultprocessingtext(27,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "density_siii_ir_ratio              "/)
        resultprocessingarray(28,:) = all_results%ArIII_IR_density
        resultprocessingtext(28,:) = (/"[ArIII] IR density                 ","{}[Ar~{\sc iii}] IR density        ", diagnostic_format, "density_ariii_ir                   "/)
        resultprocessingarray(29,:) = all_results%ariii_ir_density_ratio
        resultprocessingtext(29,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "density_ariii_ir_ratio             "/)
        resultprocessingarray(30,:) = all_results%NeIII_IR_density
        resultprocessingtext(30,:) = (/"[NeIII] IR density                 ","{}[Ne~{\sc iii}] IR density        ", diagnostic_format, "density_neiii_ir                   "/)
        resultprocessingarray(31,:) = all_results%neiii_ir_density_ratio
        resultprocessingtext(31,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "density_neiii_ir_ratio             "/)
        resultprocessingarray(32,:) = all_results%med_density
        resultprocessingtext(32,:) = (/"Medium ionisation density          ","Medium ionisation density          ", diagnostic_format, "density_med                        "/)

! medium temperature

        resultprocessingarray(33,:) = all_results%OIII_temp
        resultprocessingtext(33,:) = (/"[OIII] temperature                 ","{}[O~{\sc iii}] temperature        ", diagnostic_format, "temp_oiii                          "/)
        resultprocessingarray(34,:) = all_results%oiii_temp_ratio
        resultprocessingtext(34,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "temp_oiii_ratio                    "/)
        resultprocessingarray(35,:) = all_results%OIII_IR_temp
        resultprocessingtext(35,:) = (/"[OIII] IR temperature              ","{}[O~{\sc iii}] IR temperature     ", diagnostic_format, "temp_oiii_ir                       "/)
        resultprocessingarray(36,:) = all_results%oiii_ir_temp_ratio
        resultprocessingtext(36,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "temp_oiii_ir_ratio                 "/)
        resultprocessingarray(37,:) = all_results%OIII_UV_temp
        resultprocessingtext(37,:) = (/"[OIII] UV temperature              ","{}[O~{\sc iii}] UV temperature     ", diagnostic_format, "temp_oiii_uv                       "/)
        resultprocessingarray(38,:) = all_results%oiii_uv_temp_ratio
        resultprocessingtext(38,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "temp_oiii_uv_ratio                 "/)
        resultprocessingarray(39,:) = all_results%NeIII_temp
        resultprocessingtext(39,:) = (/"[NeIII] temperature                ","{}[Ne~{\sc iii}] temperature       ", diagnostic_format, "temp_neiii                         "/)
        resultprocessingarray(40,:) = all_results%neiii_temp_ratio
        resultprocessingtext(40,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "temp_neiii_ratio                   "/)
        resultprocessingarray(41,:) = all_results%NeIII_IR_temp
        resultprocessingtext(41,:) = (/"[NeIII] IR temperature             ","{}[Ne~{\sc iii}] IR temperature    ", diagnostic_format, "temp_neiii_ir                      "/)
        resultprocessingarray(42,:) = all_results%neiii_ir_temp_ratio
        resultprocessingtext(42,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "temp_neiii_ir_ratio                "/)
        resultprocessingarray(43,:) = all_results%ArIII_temp
        resultprocessingtext(43,:) = (/"[ArIII] temperature                ","{}[Ar~{\sc iii}] temperature       ", diagnostic_format, "temp_ariii                         "/)
        resultprocessingarray(44,:) = all_results%ariii_temp_ratio
        resultprocessingtext(44,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "temp_ariii_ratio                   "/)
        resultprocessingarray(45,:) = all_results%SIII_temp
        resultprocessingtext(45,:) = (/"[SIII] temperature                 ","{}[S~{\sc iii}] temperature        ", diagnostic_format, "temp_siii                          "/)
        resultprocessingarray(46,:) = all_results%siii_temp_ratio
        resultprocessingtext(46,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "temp_siii_ratio                    "/)
        resultprocessingarray(47,:) = all_results%med_temp
        resultprocessingtext(47,:) = (/"Medium ionisation temperature      ","Medium ionisation temperature      ", diagnostic_format, "temp_med                           "/)

!high density

        resultprocessingarray(48,:) = all_results%neiv_density
        resultprocessingtext(48,:) = (/"[NeIV] density                     ","{}[Ne~{\sc iv}] density            ", diagnostic_format, "density_neiv                       "/)
        resultprocessingarray(49,:) = all_results%neiv_density_ratio
        resultprocessingtext(49,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "density_neiv_ratio                 "/)
        resultprocessingarray(50,:) = all_results%high_density
        resultprocessingtext(50,:) = (/"High ionisation density            ","High ionisation density            ", diagnostic_format, "density_high                       "/)

!high temperature

        resultprocessingarray(51,:) = all_results%ArV_temp
        resultprocessingtext(51,:) = (/"[ArV] temperature                  ","{}[Ar~{\sc v}] temperature         ", diagnostic_format, "temp_arv                           "/)
        resultprocessingarray(52,:) = all_results%arv_temp_ratio
        resultprocessingtext(52,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "temp_arv_ratio                     "/)
        resultprocessingarray(53,:) = all_results%NeV_temp
        resultprocessingtext(53,:) = (/"[NeV] temperature                  ","{}[Ne~{\sc v}] temperature         ", diagnostic_format, "temp_nev                           "/)
        resultprocessingarray(54,:) = all_results%nev_temp_ratio
        resultprocessingtext(54,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "temp_nev_ratio                     "/)
        resultprocessingarray(55,:) = all_results%high_temp
        resultprocessingtext(55,:) = (/"High temperature                   ","High temperature                   ", diagnostic_format, "temp_high                          "/)

!Recombination line diagnostics

        resultprocessingarray(56,:) = all_results%Bal_jump_temp
        resultprocessingtext(56,:) = (/"Balmer jump temperature            ","Balmer jump temperature            ", diagnostic_format, "temp_bal_jump                      "/)

        resultprocessingarray(57,:) = all_results%Balmerdec_density
        resultprocessingtext(57,:) = (/"Balmer decrement density           ","Balmer decrement density           ", diagnostic_format, "density_bal_dec                    "/)

        resultprocessingarray(58,:) = all_results%te_5876_4471
        resultprocessingtext(58,:) = (/"He I temperature (5876/4471)       ","He temperature (5876/4471)         ", diagnostic_format, "temp_he_5876_4471                  "/)

        resultprocessingarray(59,:) = all_results%ratio_5876_4471
        resultprocessingtext(59,:) = (/"5876/4471 ratio                    ","5876/4471 ratio                    ", diagnostic_ratio_format, "temp_he_5876_4471_ratio            "/)

        resultprocessingarray(60,:) = all_results%te_6678_4471
        resultprocessingtext(60,:) = (/"He temperature (6678/4471)         ","He temperature (6678/4471)         ", diagnostic_format, "temp_he_6678_4471                  "/)

        resultprocessingarray(61,:) = all_results%ratio_6678_4471
        resultprocessingtext(61,:) = (/"6678/4471 ratio                    ","6678/4471 ratio                    ", diagnostic_ratio_format, "temp_he_6678_4471_ratio            "/)

        resultprocessingarray(62,:) = all_results%oii_te
        resultprocessingtext(62,:) = (/"OII temperature                    ","OII temperature                    ", diagnostic_format, "temp_oiirls                        "/)

        resultprocessingarray(63,:) = all_results%oii_ne
        resultprocessingtext(63,:) = (/"OII density                        ","OII density                        ", diagnostic_format, "density_oiirls                     "/)

! recombination contribution to CELs

        resultprocessingarray(64,:) = all_results%nii5754recCEL
        resultprocessingtext(64,:) = (/"NII 5754 assuming N2+ from CELs    ","N~{\sc ii}5754$_R$ (N$^{2+}$ CELs) ", adf_format, "recombination_5754_CELabund        "/)
        resultprocessingarray(65,:) = all_results%nii5754recRL
        resultprocessingtext(65,:) = (/"NII 5754 assuming N2+ from RLs     ","N~{\sc ii}5754$_R$ (N$^{2+}$ CELs) ", adf_format, "recombination_5754_RLabund         "/)
        resultprocessingarray(66,:) = all_results%oii7325recCEL
        resultprocessingtext(66,:) = (/"OII 7320,30 assuming O2+ from CELs ","O~{\sc ii}7325$_R$ (O$^{2+}$ CELs) ", adf_format, "recombination_7325_CELabund        "/)
        resultprocessingarray(67,:) = all_results%oii7325recRL
        resultprocessingtext(67,:) = (/"OII 7320,30 assuming O2+ from RLs  ","O~{\sc ii}7325$_R$ (O$^{2+}$ CELs) ", adf_format, "recombination_7325_RLabund         "/)
        resultprocessingarray(68,:) = all_results%oiii4363recCEL
        resultprocessingtext(68,:) = (/"OIII 4363 assuming O3+ from CELs   ","O~{\sc ii}4363$_R$ (O$^{3+}$ CELs) ", adf_format, "recombination_4363_CELabund        "/)
        resultprocessingarray(69,:) = all_results%oiii4363recRL
        resultprocessingtext(69,:) = (/"OIII 4363 assuming O3+ from RLs    ","O~{\sc ii}4363$_R$ (O$^{3+}$ CELs) ", adf_format, "recombination_4363_RLabund         "/)

!CEL abundances

        resultprocessingarray(70,:) = all_results%NC_abund_CEL
        resultprocessingtext(70,:) = (/"C0/H                               ","C$^{0}$/H                          ", abundances_format, "abund_cel_nc                       "/)
        resultprocessingarray(71,:) = all_results%cii_abund_CEL
        resultprocessingtext(71,:) = (/"C+/H                               ","C$^{+}$/H                          ", abundances_format, "abund_cel_cii                      "/)
        resultprocessingarray(72,:) = all_results%ciii_abund_CEL
        resultprocessingtext(72,:) = (/"C2+/H                              ","C$^{2+}$/H                         ", abundances_format, "abund_cel_ciii                     "/)
        resultprocessingarray(73,:) = all_results%civ_abund_CEL
        resultprocessingtext(73,:) = (/"C3+/H                              ","C$^{3+}$/H                         ", abundances_format, "abund_cel_civ                      "/)
        resultprocessingarray(74,:) = all_results%c_icf_CEL
        resultprocessingtext(74,:) = (/"icf(C)                             ","icf(C)                             ", abundances_format, "icf_cel_c                          "/)
        resultprocessingarray(75,:) = all_results%C_abund_CEL
        resultprocessingtext(75,:) = (/"C/H                                ","C$^{}$/H                           ", abundances_format, "abund_cel_c                        "/)
        resultprocessingarray(76,:) = all_results%nii_abund_CEL
        resultprocessingtext(76,:) = (/"N+/H                               ","N$^{+}$/H                          ", abundances_format, "abund_cel_nii                      "/)
        resultprocessingarray(77,:) = all_results%niii_abund_CEL
        resultprocessingtext(77,:) = (/"N2+/H                              ","N$^{2+}$/H                         ", abundances_format, "abund_cel_niii                     "/)
        resultprocessingarray(78,:) = all_results%niv_abund_CEL
        resultprocessingtext(78,:) = (/"N3+/H                              ","N$^{3+}$/H                         ", abundances_format, "abund_cel_niv                      "/)
        resultprocessingarray(79,:) = all_results%nv_abund_CEL
        resultprocessingtext(79,:) = (/"N4+/H                              ","N$^{4+}$/H                         ", abundances_format, "abund_cel_nv                       "/)
        resultprocessingarray(80,:) = all_results%n_icf_CEL
        resultprocessingtext(80,:) = (/"icf(N)                             ","icf(N)                             ", abundances_format, "icf_cel_n                          "/)
        resultprocessingarray(81,:) = all_results%N_abund_CEL
        resultprocessingtext(81,:) = (/"N/H                                ","N$^{}$/H                           ", abundances_format, "abund_cel_n                        "/)
        resultprocessingarray(82,:) = all_results%NO_abund_CEL
        resultprocessingtext(82,:) = (/"O0/H                               ","O$^{0}$/H                          ", abundances_format, "abund_cel_no                       "/)
        resultprocessingarray(83,:) = all_results%Oii_abund_CEL
        resultprocessingtext(83,:) = (/"O+/H                               ","O$^{+}$/H                          ", abundances_format, "abund_cel_oii                      "/)
        resultprocessingarray(84,:) = all_results%Oiii_abund_CEL
        resultprocessingtext(84,:) = (/"O2+/H                              ","O$^{2+}$/H                         ", abundances_format, "abund_cel_oiii                     "/)
        resultprocessingarray(85,:) = all_results%Oiv_abund_CEL
        resultprocessingtext(85,:) = (/"O3+/H                              ","O$^{3+}$/H                         ", abundances_format, "abund_cel_oiv                      "/)
        resultprocessingarray(86,:) = all_results%o_icf_CEL
        resultprocessingtext(86,:) = (/"icf(O)                             ","icf(O)                             ", abundances_format, "icf_cel_o                          "/)
        resultprocessingarray(87,:) = all_results%O_abund_CEL
        resultprocessingtext(87,:) = (/"O/H                                ","O$^{}$/H                           ", abundances_format, "abund_cel_o                        "/)
        resultprocessingarray(88,:) = all_results%NeII_abund_CEL
        resultprocessingtext(88,:) = (/"Ne+/H                              ","Ne$^{+}$/H                         ", abundances_format, "abund_cel_neii                     "/)
        resultprocessingarray(89,:) = all_results%NeIII_abund_CEL
        resultprocessingtext(89,:) = (/"Ne2+/H                             ","Ne$^{2+}$/H                        ", abundances_format, "abund_cel_neiii                    "/)
        resultprocessingarray(90,:) = all_results%NeIV_abund_CEL
        resultprocessingtext(90,:) = (/"Ne3+/H                             ","Ne$^{3+}$/H                        ", abundances_format, "abund_cel_neiv                     "/)
        resultprocessingarray(91,:) = all_results%NeV_abund_CEL
        resultprocessingtext(91,:) = (/"Ne4+/H                             ","Ne$^{4+}$/H                        ", abundances_format, "abund_cel_nev                      "/)
        resultprocessingarray(92,:) = all_results%ne_icf_CEL
        resultprocessingtext(92,:) = (/"icf(Ne)                            ","icf(Ne)                            ", abundances_format, "icf_cel_ne                         "/)
        resultprocessingarray(93,:) = all_results%Ne_abund_CEL
        resultprocessingtext(93,:) = (/"Ne/H                               ","Ne$^{}$/H                          ", abundances_format, "abund_cel_ne                       "/)
        resultprocessingarray(94,:) = all_results%ArIII_abund_CEL
        resultprocessingtext(94,:) = (/"Ar2+/H                             ","Ar$^{2+}$/H                        ", abundances_format, "abund_cel_ariii                    "/)
        resultprocessingarray(95,:) = all_results%ArIV_abund_CEL
        resultprocessingtext(95,:) = (/"Ar3+/H                             ","Ar$^{3+}$/H                        ", abundances_format, "abund_cel_ariv                     "/)
        resultprocessingarray(96,:) = all_results%ArV_abund_CEL
        resultprocessingtext(96,:) = (/"Ar4+/H                             ","Ar$^{4+}$/H                        ", abundances_format, "abund_cel_arv                      "/)
        resultprocessingarray(97,:) = all_results%ar_icf_CEL
        resultprocessingtext(97,:) = (/"icf(Ar)                            ","icf(Ar)                            ", abundances_format, "icf_cel_ar                         "/)
        resultprocessingarray(98,:) = all_results%Ar_abund_CEL
        resultprocessingtext(98,:) = (/"Ar/H                               ","Ar$^{}$/H                          ", abundances_format, "abund_cel_ar                       "/)
        resultprocessingarray(99,:) = all_results%SII_abund_CEL
        resultprocessingtext(99,:) = (/"S+/H                               ","S$^{+}$/H                          ", abundances_format, "abund_cel_sii                      "/)
        resultprocessingarray(100,:) = all_results%SIII_abund_CEL
        resultprocessingtext(100,:) = (/"S2+/H                              ","S$^{2+}$/H                         ", abundances_format, "abund_cel_siii                     "/)
        resultprocessingarray(101,:) = all_results%s_icf_CEL
        resultprocessingtext(101,:) = (/"icf(S)                             ","icf(S)                             ", abundances_format, "icf_cel_s                          "/)
        resultprocessingarray(102,:) = all_results%S_abund_CEL
        resultprocessingtext(102,:) = (/"S/H                                ","S$^{}$/H                           ", abundances_format, "abund_cel_s                        "/)
        resultprocessingarray(103,:) = all_results%ClIII_abund_CEL
        resultprocessingtext(103,:) = (/"Cl2+/H                             ","Cl$^{2+}$/H                        ", abundances_format, "abund_cel_cliii                    "/)
        resultprocessingarray(104,:) = all_results%cl_icf_CEL
        resultprocessingtext(104,:) = (/"icf(Cl)                            ","icf(Cl)                            ", abundances_format, "icf_cel_cl                         "/)
        resultprocessingarray(105,:) = all_results%Cl_abund_CEL
        resultprocessingtext(105,:) = (/"Cl/H                               ","Cl$^{}$/H                          ", abundances_format, "abund_cel_cl                       "/)

!ORL abundances

        resultprocessingarray(106,:) = all_results%Hei_abund_ORL
        resultprocessingtext(106,:) = (/"He+/H                              ","He$^{+}$/H                         ", abundances_format, "abund_orl_hei                      "/)
        resultprocessingarray(107,:) = all_results%Heii_abund_ORL
        resultprocessingtext(107,:) = (/"He2+/H                             ","He$^{2+}$/H                        ", abundances_format, "abund_orl_heii                     "/)
        resultprocessingarray(108,:) = all_results%He_abund_ORL
        resultprocessingtext(108,:) = (/"He/H                               ","He/H                               ", abundances_format, "abund_orl_he                       "/)
        resultprocessingarray(109,:) = all_results%Cii_abund_ORL
        resultprocessingtext(109,:) = (/"C2+/H                              ","C$^{2+}$/H                         ", abundances_format, "abund_orl_cii                      "/)
        resultprocessingarray(110,:) = all_results%Ciii_abund_ORL
        resultprocessingtext(110,:) = (/"C3+/H                              ","C$^{3+}$/H                         ", abundances_format, "abund_orl_ciii                     "/)
        resultprocessingarray(111,:) = all_results%c_icf_ORL
        resultprocessingtext(111,:) = (/"icf(C)                             ","icf(C)                             ", abundances_format, "icf_orl_c                          "/)
        resultprocessingarray(112,:) = all_results%C_abund_ORL
        resultprocessingtext(112,:) = (/"C/H                                ","C/H                                ", abundances_format, "abund_orl_c                        "/)
        resultprocessingarray(113,:) = all_results%Nii_v3_abund_ORL
        resultprocessingtext(113,:) = (/"N2+/H (V3)                         ","N$^{2+}$/H (V3)                    ", abundances_format, "abund_orl_nii_v3                   "/)
        resultprocessingarray(114,:) = all_results%Nii_v5_abund_ORL
        resultprocessingtext(114,:) = (/"N2+/H (V5)                         ","N$^{2+}$/H (V5)                    ", abundances_format, "abund_orl_nii_v5                   "/)
        resultprocessingarray(115,:) = all_results%Nii_v8_abund_ORL
        resultprocessingtext(115,:) = (/"N2+/H (V8)                         ","N$^{2+}$/H (V8)                    ", abundances_format, "abund_orl_nii_v8                   "/)
        resultprocessingarray(116,:) = all_results%Nii_v12_abund_ORL
        resultprocessingtext(116,:) = (/"N2+/H (V12)                        ","N$^{2+}$/H (V12)                   ", abundances_format, "abund_orl_nii_v12                  "/)
        resultprocessingarray(117,:) = all_results%Nii_v20_abund_ORL
        resultprocessingtext(117,:) = (/"N2+/H (V20)                        ","N$^{2+}$/H (V20)                   ", abundances_format, "abund_orl_nii_v20                  "/)
        resultprocessingarray(118,:) = all_results%Nii_v28_abund_ORL
        resultprocessingtext(118,:) = (/"N2+/H (V28)                        ","N$^{2+}$/H (V28)                   ", abundances_format, "abund_orl_nii_v28                  "/)
        resultprocessingarray(119,:) = all_results%Nii_3d4f_abund_ORL
        resultprocessingtext(119,:) = (/"N2+/H (3d-4f)                      ","N$^{2+}$/H (3d-4f)                 ", abundances_format, "abund_orl_nii_3d4f                 "/)
        resultprocessingarray(120,:) = all_results%Nii_abund_ORL
        resultprocessingtext(120,:) = (/"N2+/H                              ","N$^{2+}$/H                         ", abundances_format, "abund_orl_nii                      "/)
        resultprocessingarray(121,:) = all_results%Niii_abund_ORL
        resultprocessingtext(121,:) = (/"N3+/H                              ","N$^{3+}$/H                         ", abundances_format, "abund_orl_niii                     "/)
        resultprocessingarray(122,:) = all_results%n_icf_ORL
        resultprocessingtext(122,:) = (/"icf(N)                             ","icf(N)                             ", abundances_format, "icf_orl_n                          "/)
        resultprocessingarray(123,:) = all_results%N_abund_ORL
        resultprocessingtext(123,:) = (/"N/H                                ","N/H                                ", abundances_format, "abund_orl_n                        "/)
        resultprocessingarray(124,:) = all_results%Oii_v1_abund_ORL
        resultprocessingtext(124,:) = (/"O2+/H (V1)                         ","O$^{2+}$/H (V1)                    ", abundances_format, "abund_orl_oii_v1                   "/)
        resultprocessingarray(125,:) = all_results%Oii_v2_abund_ORL
        resultprocessingtext(125,:) = (/"O2+/H (V2)                         ","O$^{2+}$/H (V2)                    ", abundances_format, "abund_orl_oii_v2                   "/)
        resultprocessingarray(126,:) = all_results%Oii_v5_abund_ORL
        resultprocessingtext(126,:) = (/"O2+/H (V5)                         ","O$^{2+}$/H (V5)                    ", abundances_format, "abund_orl_oii_v5                   "/)
        resultprocessingarray(127,:) = all_results%Oii_v10_abund_ORL
        resultprocessingtext(127,:) = (/"O2+/H (V10)                        ","O$^{2+}$/H (V10)                   ", abundances_format, "abund_orl_oii_v10                  "/)
        resultprocessingarray(128,:) = all_results%Oii_v11_abund_ORL
        resultprocessingtext(128,:) = (/"O2+/H (V11)                        ","O$^{2+}$/H (V11)                   ", abundances_format, "abund_orl_oii_v11                  "/)
        resultprocessingarray(129,:) = all_results%Oii_v12_abund_ORL
        resultprocessingtext(129,:) = (/"O2+/H (V12)                        ","O$^{2+}$/H (V12)                   ", abundances_format, "abund_orl_oii_v12                  "/)
        resultprocessingarray(130,:) = all_results%Oii_v19_abund_ORL
        resultprocessingtext(130,:) = (/"O2+/H (V19)                        ","O$^{2+}$/H (V19)                   ", abundances_format, "abund_orl_oii_v19                  "/)
        resultprocessingarray(131,:) = all_results%Oii_v20_abund_ORL
        resultprocessingtext(131,:) = (/"O2+/H (V20)                        ","O$^{2+}$/H (V20)                   ", abundances_format, "abund_orl_oii_v20                  "/)
        resultprocessingarray(132,:) = all_results%Oii_v25_abund_ORL
        resultprocessingtext(132,:) = (/"O2+/H (V25)                        ","O$^{2+}$/H (V25)                   ", abundances_format, "abund_orl_oii_v25                  "/)
        resultprocessingarray(133,:) = all_results%Oii_v28_abund_ORL
        resultprocessingtext(133,:) = (/"O2+/H (V28)                        ","O$^{2+}$/H (V28)                   ", abundances_format, "abund_orl_oii_v28                  "/)
        resultprocessingarray(134,:) = all_results%Oii_v33_abund_ORL
        resultprocessingtext(134,:) = (/"O2+/H (V33)                        ","O$^{2+}$/H (V33)                   ", abundances_format, "abund_orl_oii_v33                  "/)
        resultprocessingarray(135,:) = all_results%Oii_3d4f_abund_ORL
        resultprocessingtext(135,:) = (/"O2+/H (3d-4f)                      ","O$^{2+}$/H (3d-4f)                 ", abundances_format, "abund_orl_oii_3d4f                 "/)
        resultprocessingarray(136,:) = all_results%Oii_abund_ORL
        resultprocessingtext(136,:) = (/"O2+/H                              ","O$^{2+}$/H                         ", abundances_format, "abund_orl_oii                      "/)
        resultprocessingarray(137,:) = all_results%O_icf_ORL
        resultprocessingtext(137,:) = (/"icf(O)                             ","icf(O)                             ", abundances_format, "icf_orl_o                          "/)
        resultprocessingarray(138,:) = all_results%O_abund_ORL
        resultprocessingtext(138,:) = (/"O/H                                ","O/H                                ", abundances_format, "abund_orl_o                        "/)
        resultprocessingarray(139,:) = all_results%Neii_abund_ORL
        resultprocessingtext(139,:) = (/"Ne2+/H                             ","Ne$^{2+}$/H                        ", abundances_format, "abund_orl_neii                     "/)
        resultprocessingarray(140,:) = all_results%ne_icf_ORL
        resultprocessingtext(140,:) = (/"icf(Ne)                            ","icf(Ne)                            ", abundances_format, "icf_orl_ne                         "/)
        resultprocessingarray(141,:) = all_results%Ne_abund_ORL
        resultprocessingtext(141,:) = (/"Ne/H                               ","Ne/H                               ", abundances_format, "abund_orl_ne                       "/)

!strong line abundances

        resultprocessingarray(142,:) = all_results%O_R23_upper
        resultprocessingtext(142,:) = (/"O/H (R23 upper)                    ","O/H (R23 upper)                    ", abundances_format, "abund_o_r23_upper                  "/)
        resultprocessingarray(143,:) = all_results%O_R23_lower
        resultprocessingtext(143,:) = (/"O/H (R23 lower)                    ","O/H (R23 lower)                    ", abundances_format, "abund_o_r23_lower                  "/)
        resultprocessingarray(144,:) = all_results%O_N2
        resultprocessingtext(144,:) = (/"O/H (N2)                           ","O/H (N2)                           ", abundances_format, "abund_o_n2                         "/)
        resultprocessingarray(145,:) = all_results%O_O3N2
        resultprocessingtext(145,:) = (/"O/H (O3N2)                         ","O/H (O3N2)                         ", abundances_format, "abund_o_o3n2                       "/)
        resultprocessingarray(146,:) = all_results%O_Ar3O3
        resultprocessingtext(146,:) = (/"O/H (Ar3O3)                        ","O/H (Ar3O3)                        ", abundances_format, "abund_o_ar3o3                      "/)
        resultprocessingarray(147,:) = all_results%O_S3O3
        resultprocessingtext(147,:) = (/"O/H (S3O3)                         ","O/H (S3O3)                         ", abundances_format, "abund_o_s3o3                       "/)

!adfs

        resultprocessingarray(148,:) = all_results%adf_o2plus
        resultprocessingtext(148,:) = (/"adf (O2+/H)                        ","adf (O$^{2+}$/H)                   ", adf_format, "adf_o2plus                         "/)
        resultprocessingarray(149,:) = all_results%adf_o
        resultprocessingtext(149,:) = (/"adf (O/H)                          ","adf (O/H)                          ", adf_format, "adf_o                              "/)
        resultprocessingarray(150,:) = all_results%adf_n2plus
        resultprocessingtext(150,:) = (/"adf (N2+/H)                        ","adf (N$^{2+}$/H)                   ", adf_format, "adf_n2plus                         "/)
        resultprocessingarray(151,:) = all_results%adf_n
        resultprocessingtext(151,:) = (/"adf (N/H)                          ","adf (N/H)                          ", adf_format, "adf_n                              "/)
        resultprocessingarray(152,:) = all_results%adf_c2plus
        resultprocessingtext(152,:) = (/"adf (C2+/H)                        ","adf (C$^{2+}$/H)                   ", adf_format, "adf_c2plus                         "/)
        resultprocessingarray(153,:) = all_results%adf_c
        resultprocessingtext(153,:) = (/"adf (C/H)                          ","adf (C/H)                          ", adf_format, "adf_c                              "/)
        resultprocessingarray(154,:) = all_results%adf_ne2plus
        resultprocessingtext(154,:) = (/"adf (Ne2+/H)                       ","adf (Ne$^{2+}$/H)                  ", adf_format, "adf_ne2plus                        "/)
        resultprocessingarray(155,:) = all_results%adf_ne
        resultprocessingtext(155,:) = (/"adf (Ne/H)                         ","adf (Ne/H)                         ", adf_format, "adf_ne                             "/)

!open the files and write the headers

        open (650,FILE=trim(filename)//"_results", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
        open (651,FILE=trim(filename)//"_results.tex", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')

        write (650,*) "NEAT (nebular empirical analysis tool)"
        write (650,*) "======================================"
        write (650,*)
        write (650,*) "Analysis of file ",trim(filename)
        write (650,*) "Command line: ",trim(commandline)
        write (650,*)

        write (651,*) "\noindent{\Large {\sc neat} (nebular empirical analysis tool)}"
        write (651,*) "\hrule"
        write (651,*) "\vspace{0.3cm}"
        write (651,*) "\noindent Analysis of file {\tt ",trim(filename),"}\newline"
        write (651,*) "\noindent Command line: {\tt ",trim(commandline),"}\newline"
        write (651,*) "\begin{longtable}[l]{ll}"

!next, loop through the results, processing and printing

        do j=1,155

! here we put some if statements to put things into conveniently separate bits

          if (j .eq. 1) then
            write (650,"(/A,/A/)") "Extinction","=========="
            write (651,*) "\multicolumn{2}{l}{Extinction}\\ \hline"
          elseif (j .eq. 2) then
            write (650,"(/A,/A)") "Diagnostics","==========="
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Diagnostics}\\ \hline"
            write (650,"(/A,/A/)") "Low ionisation densities","-----------"
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Low ionisation densities}\\ \hline"
          elseif (j .eq. 7) then
            write (650,"(/A,/A/)") "Low ionisation temperatures","-----------"
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Low ionisation temperatures}\\ \hline"
          elseif (j .eq. 18) then
            write (650,"(/A,/A/)") "Medium ionisation densities","-----------"
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Medium ionisation densities}\\ \hline"
          elseif (j .eq. 33) then
            write (650,"(/A,/A/)") "Medium ionisation temperatures","-----------"
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Medium ionisation temperatures}\\ \hline"
          elseif (j .eq. 48) then
            write (650,"(/A,/A/)") "High ionisation densities","-----------"
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{High ionisation densities}\\ \hline"
          elseif (j .eq. 51) then
            write (650,"(/A,/A/)") "High ionisation temperatures","-----------"
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{High ionisation temperatures}\\ \hline"
          elseif (j .eq. 56) then
            write (650,"(/A,/A/)") "Recombination line diagnostics","-----------"
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Recombination line diagnostics}\\ \hline"
          elseif (j .eq. 64) then
            write (650,"(/A,/A/)") "Recombination contribution to CELs (%)","-----------"
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Recombination contribution to CELs (\%)}\\ \hline"
          elseif (j .eq. 70) then
            write (650,"(/A,/A/)") "CEL abundances","=============="
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{CEL abundances}\\ \hline"
          elseif (j .eq. 108) then
            write (650,"(/A,/A/)") "ORL abundances","=============="
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{ORL abundances}\\ \hline" 
          elseif (j .eq. 109 .or. j .eq. 113 .or. j .eq. 124 .or. j .eq. 139) then 
            write (650,"(/A,/A/)") 
            write (651,*) "\\"
          elseif (j .eq. 142) then
            write (650,"(/A,/A/)") "Strong line abundances","======================"
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Strong line abundances}\\ \hline" 
          elseif (j .eq. 148) then
            write (650,"(/A,/A/)") "Abundance discrepancy factors","============================="
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Abundance discrepancy factors}\\ \hline"
          endif

! this writes the results to the plain text and latex summary files
          
          quantity_result=resultprocessingarray(j,:)
          call write_uncertainties(quantity_result,uncertainty_array,resultprocessingtext(j,1),resultprocessingtext(j,2),resultprocessingtext(j,3),filename, resultprocessingtext(j,4), verbosity,nbins,nperbin,runs)

        enddo

!write ends of files, close

        write (651,*) "\end{longtable}"

        close(650)
        close(651)

        print *,gettime(),": Results written to files:"
        print *,"              ",trim(filename),"_results"
        print *,"              ",trim(filename),"_results.tex"

        print *
        print *,gettime(),": all done"
        print *

contains

        subroutine randomizer(linelist, listlength, norp)
        ! from http://www.netlib.org/random/random.f90

        TYPE(line), dimension(listlength) :: linelist
        INTEGER :: IO, I, j, listlength
        real(kind=dp) :: temp4
        LOGICAL, intent(in) :: norp

        real(kind=dp) :: fn_val
        real(kind=dp)     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           &
                    r1 = 0.27597, r2 = 0.27846, u, v, x, y, q
        real(kind=dp) :: half
        real(kind=dp) :: newmean, newsnr, snr

        half = 0.5

        I = 1
        IO=0

        do j = 1,listlength
          DO
            CALL RANDOM_NUMBER(u)
            CALL RANDOM_NUMBER(v)
            v = 1.7156 * (v - half)
            x = u - s
            y = ABS(v) - t
            q = x**2 + y*(a*y - b*x)
            IF (q < r1) EXIT
            IF (q > r2) CYCLE
            IF (v**2 < -4.0*LOG(u)*u**2) EXIT
          END DO
          fn_val = v/u

!                if (j==1) R=3.1+(0.15*fn_val)

          if (linelist(j)%intensity/linelist(j)%int_err .gt. 6.0 .or. norp) then !normal distribution

            temp4=linelist(j)%intensity+(fn_val*linelist(j)%int_err)
            if(temp4 < 0) temp4 = 0.D0
            linelist(j)%intensity = temp4

          elseif (linelist(j)%int_err .ge. linelist(j)%intensity) then !it's an upper limit, take number from semi-gaussian distribution with peak at zero and 5 sigma = intensity
            linelist(j)%intensity = abs(fn_val)*0.2*linelist(j)%intensity
          else !if SN<6, then take lognormal distribution, parameters from Rola & Pelat (1994)
                     !for SN<6, the actual mean is derived from the observed mean using
            snr = linelist(j)%intensity/linelist(j)%int_err
            newmean = 0.0765957/(snr**2) + 1.86037/snr - 0.309695
                     !the actual standard deviation is derived from the observed using
            newsnr = -1.11329/(snr**3) + 1.8542/(snr**2) - 0.288222/snr + 0.18018
                     !(fits to the data in Rola & Pelat's table 6)
                     !the distributions in table 6 give the mean and sigma of log-normal distributions of S/N(obs), given S/N(true).  We don't know S/N(true) but using the distributions as that of the factor by which line fluxes are overestimated is equivalent.  So,
            temp4 = exp(fn_val*newsnr + newmean)
            if (temp4 < 0 ) temp4 = 0.D0
            linelist(j)%intensity = linelist(j)%intensity / temp4

          endif
        end do

        end subroutine

        SUBROUTINE init_random_seed()
          INTEGER :: i, n, clock
          INTEGER, DIMENSION(:), ALLOCATABLE :: seed

          n=20
          i=n
          CALL RANDOM_SEED(size = n)
          ALLOCATE(seed(n))

          CALL SYSTEM_CLOCK(COUNT=clock)

          seed = clock + 37 * (/ (i - 1, i = 1, n) /)
          CALL RANDOM_SEED(PUT = seed)

          DEALLOCATE(seed)
        END SUBROUTINE

subroutine write_uncertainties(input_array, uncertainty_array, plaintext, latextext, itemformat, filename, suffix, verbosity,nbins,nperbin,runs)

!wrapper for the get_uncertainties routine, if called it will do the
!uncertainty calculation and also write the binned and unbinned results to files
!as necessary.
implicit none
real(kind=dp) :: input_array(:)
real(kind=dp), intent(out) :: uncertainty_array(3)
real(kind=dp), dimension (:,:), allocatable :: binned_quantity_result
character(len=35), intent(in) :: plaintext, latextext
character(len=35), intent(in) :: itemformat
character(len=512), intent(in) :: filename
character(len=25), intent(in) :: suffix
logical :: unusual
integer :: verbosity !1 = write summaries, binned and full results; 2=write summaries and binned results; 3=write summaries only
integer :: nbins,nperbin,runs

if(size(input_array) .eq. 1) then
!just one iteration, write out the result without uncertainties
  if (input_array(1) .gt. 0.) then
    write (650,itemformat) plaintext,input_array(1)
    write (651,*) latextext," & $", trim(latex_number(input_array(1))),"$\\"
  else
    write (650,*) plaintext,"--"
    write (651,*) latextext," & -- \\"
  endif
else

call get_uncertainties(input_array, binned_quantity_result, uncertainty_array, unusual,nbins,nperbin,runs)

!write out the binned results with the mode+-uncertainties at the top

if (verbosity .lt. 3) then

  if (allocated(binned_quantity_result) .and. maxval(uncertainty_array) .gt. 0.) then
    OPEN(850, FILE=trim(filename)//"_"//trim(suffix)//"_binned", STATUS='REPLACE',ACCESS='SEQUENTIAL', ACTION='WRITE')
    write(unit = 850,FMT=*) uncertainty_array(2),uncertainty_array(2)-uncertainty_array(1),uncertainty_array(3)+uncertainty_array(2)
    write(unit = 850,FMT=*)
  !  do i=1,ii-1
    do i=1,size(binned_quantity_result, 1)
  ! this hacky condition is because for reasons I can't work out right now, the
  ! number of bins allocated to the array is always too large and the last few end
  ! up being full of zeros.  to be fixed soon hopefully.  RW 16/11/2012
      if (binned_quantity_result(i,1) .gt. 0. .and. binned_quantity_result(i,2) .gt. 0) then
        write(unit = 850,FMT=*) binned_quantity_result(i,1),binned_quantity_result(i,2)
      endif
    end do
    close(850)
  endif

endif

!write the unbinned results to file

if (verbosity .lt. 2) then

  OPEN(850, FILE=trim(filename)//"_"//trim(suffix), STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
  do i=1,size(input_array)
    write(unit = 850,FMT=*) input_array(i)
  end do
  close(850)

endif

!write derived value and uncertainties to the summary files

  if (maxval(uncertainty_array) .gt. 0.) then !if this condition is not true, array will be full of zeroes
    if (unusual) write (650,itemformat) "Warning! Unusual probability distribution.  You should inspect this one:"
    write (650,itemformat) plaintext,uncertainty_array(2),uncertainty_array(3),-uncertainty_array(1)
!    write (651,*) latextext," & $",trim(latex_number(uncertainty_array(2))),"^{",trim(latex_number(uncertainty_array(3))),"}_{",trim(latex_number(-uncertainty_array(1))),"}$ \\"
    if (uncertainty_array(1) .eq. uncertainty_array(3)) then
      write (651,"(A,' & ${',A,'}\pm{',A,'}$ \\')") latextext,trim(latex_number(uncertainty_array(2))),trim(latex_number(uncertainty_array(1)))
    else
      write (651,"(A,' & ${',A,'}^{+',A,'}_{',A,'}$ \\')") latextext,trim(latex_number(uncertainty_array(2))),trim(latex_number(uncertainty_array(3))),trim(latex_number(-uncertainty_array(1)))
    endif
  else
    write (650,*) plaintext,"--"
    write (651,*) latextext,"& -- \\"
  endif

endif !end of condition checking whether one or more iterations were done

end subroutine write_uncertainties

subroutine get_uncertainties(input_array, binned_quantity_result, uncertainty_array, unusual,nbins,nperbin,runs)

implicit none
real(kind=dp) :: input_array(:)
real(kind=dp), intent(out) :: uncertainty_array(3)
real(kind=dp), dimension(:), allocatable :: bintemp
real(kind=dp) :: binsize=0d0
real(kind=dp), dimension (:,:), allocatable, intent(out) :: binned_quantity_result
integer :: arraysize, abovepos, belowpos
integer :: nbins, nperbin, runs
integer, dimension(1) :: maxpos
real(kind=dp) :: mean=0d0, sd=0d0
real(kind=dp) :: mean_log=0d0, sd_log=0d0
real(kind=dp) :: mean_exp=0d0, sd_exp=0d0
real(kind=dp), dimension(3,3) :: sds=0d0
real(kind=dp) :: tolerance=0d0
real(kind=dp) :: binstart,binend,binwidth
logical :: unusual !if not true, then the distribution is normal, log normal or exp normal


unusual = .false.

if (allocated(binned_quantity_result)) deallocate(binned_quantity_result)

tolerance=0.02 ! to determine whether a probability distribution is normal, log normal or exp normal, we calculate the mean and standard deviation of the original array, the logarithm of the values, and the exponent of the values.  We then determine the fractions of the distribution lying within 1, 2 and 3 standard deviations of the mean in each case.  If these correspond to the 0.68-0.95-0.99 expected for a normal distribution then we know that the distribution is normal, log normal or exp normal and record the uncertainties accordingly.  The fractions must be (0.6827+-tol), (0.9545-tol/2) and (0.9973+-tol/5) for the subroutine to assign a distribution.
sds=0.D0

!binned_quantity_result = 0.D0
uncertainty_array = (/0.0,0.0,0.0/)
arraysize = size(input_array)

!todo: include a verbosity switch
!verbosity: 0=summary file only, 1=summary file + binned results, 2=summary, binned and full
!set verbosity to 0 automatically for single iteration, then we can use the same routine for any n

!sort the input array into ascending order

call qsort(input_array)

! bin the array

arraysize = size(input_array)
binsize=(input_array(nint(0.841*arraysize)) - input_array(nint(0.159*arraysize)))/5d0

if (binsize .gt. 0d0) then

!quantize the input array by taking the integer of each value divided by the binsize, and multiplying by the binsize.

  allocate(bintemp(arraysize))
!  bintemp = binsize*nint(input_array/binsize)
!  nbins=abs(nint((maxval(bintemp)-minval(bintemp))/binsize))+1
  allocate(binned_quantity_result(nbins,2))

  binned_quantity_result=0.D0

  do i=1,nbins
     if(i.eq.1) then
        binstart=input_array(1)/2d0
        binend=(input_array(1+(i)*nperbin)+input_array((i)*nperbin))/2d0
     elseif(i.eq.nbins) then
        binstart=(input_array(1+(i-1)*nperbin)+input_array((i-1)*nperbin))/2d0
        binend=input_array(arraysize)*2d0
     else
        binstart=(input_array(1+(i-1)*nperbin)+input_array((i-1)*nperbin))/2d0
        binend=(input_array(1+(i)*nperbin)+input_array((i)*nperbin))/2d0
     endif
     binwidth=binend-binstart
     binned_quantity_result(i,1)=(binstart+binend)/2d0
     binned_quantity_result(i,2)=dble(nperbin)/dble(runs)/binwidth !probability density in the current bin
  enddo

!then, go through the quantized array and count the number of occurrences of each value
!!$
!!$  comp=bintemp(1)
!!$  bincount=0
!!$  ii=1
!!$
!!$  do i=1,arraysize
!!$    if (bintemp(i).eq.comp) then
!!$      bincount=bincount+1 
!!$    else 
!!$      binned_quantity_result(ii,:)=(/comp, dble(bincount)/) 
!!$      comp=bintemp(i)
!!$      bincount=0
!!$      ii=ii+1 
!!$    endif
!!$  enddo

!mode of distribution is the value with the highest bin count

  maxpos=maxloc(binned_quantity_result(:,2))
  uncertainty_array(2)=binned_quantity_result(maxpos(1),1)

!now find the value in the unbinned array that's closest to the mode

  bintemp=0.D0
  bintemp = abs(input_array - uncertainty_array(2))
  i = minloc(bintemp,1)
  uncertainty_array(2)=input_array(i)

!now find the values such that 68.2% of all values lie within the uncertainties

  abovepos = i+nint(0.341*arraysize)
  belowpos = i-nint(0.341*arraysize)

  if (abovepos>arraysize) then
    uncertainty_array(3) = input_array(arraysize) - uncertainty_array(2)
  else
    uncertainty_array(3) = input_array(abovepos) - uncertainty_array(2)
  endif

  if (belowpos<1) then
    uncertainty_array(1) = uncertainty_array(2) ! no lower limit so the negative uncertainty is equal to the value
  else
    uncertainty_array(1) = uncertainty_array(2) - input_array(belowpos)
  endif

!simple test for normality - calculate the mean and standard deviation of the array, determine the number of values within 1, 2 and 3 sigma of the mean.
!if the fractions are close to 68.27, 95.45 and 99.73 then it's normal
!also calculate for log and exp

  mean=sum(input_array)/arraysize
  mean_log=sum(log(input_array))/arraysize
  mean_exp=sum(exp(input_array))/arraysize

  do i=1,arraysize
    sd=sd+(input_array(i)-mean)**2
    sd_log=sd_log+(log(input_array(i))-mean_log)**2
    sd_exp=sd_exp+(exp(input_array(i))-mean_exp)**2
  end do

  sd=(sd/arraysize)**0.5
  sd_log=(sd_log/arraysize)**0.5
  sd_exp=(sd_exp/arraysize)**0.5

  do i=1,arraysize
    if (abs(input_array(i)-mean) .lt. sd) sds(1,1)=sds(1,1)+1
    if (abs(input_array(i)-mean) .lt. (2*sd)) sds(1,2)=sds(1,2)+1
    if (abs(input_array(i)-mean) .lt. (3*sd)) sds(1,3)=sds(1,3)+1

    if (abs(log(input_array(i))-mean_log) .lt. sd_log) sds(2,1)=sds(2,1)+1
    if (abs(log(input_array(i))-mean_log) .lt. (2*sd_log)) sds(2,2)=sds(2,2)+1
    if (abs(log(input_array(i))-mean_log) .lt. (3*sd_log)) sds(2,3)=sds(2,3)+1

    if (abs(exp(input_array(i))-mean_exp) .lt. sd_exp) sds(3,1)=sds(3,1)+1
    if (abs(exp(input_array(i))-mean_exp) .lt. (2*sd_exp)) sds(3,2)=sds(3,2)+1
    if (abs(exp(input_array(i))-mean_exp) .lt. (3*sd_exp)) sds(3,3)=sds(3,3)+1
  end do

  sds=sds/arraysize

  if (abs(sds(1,1)-0.6827).lt.tolerance .and. abs(sds(1,2)-0.9545).lt.(tolerance*0.5) .and. abs(sds(1,3)-0.9973).lt. (tolerance*0.1)) then 
    uncertainty_array(:)=(/sd,mean,sd/)
  elseif (abs(sds(2,1)-0.6827).lt.tolerance .and. abs(sds(2,2)-0.9545).lt.(tolerance*0.5) .and. abs(sds(2,3)-0.9973).lt. (tolerance*0.1)) then 
    uncertainty_array(:)=(/exp(mean_log)-exp(mean_log-sd_log),exp(mean_log),exp(mean_log+sd_log)-exp(mean_log)/)
  elseif (abs(sds(3,1)-0.6827).lt.tolerance .and. abs(sds(3,2)-0.9545).lt.(tolerance*0.5) .and. abs(sds(3,3)-0.9973).lt. (tolerance*0.1)) then
    uncertainty_array(:)=(/log(mean_exp+sd_exp)-log(mean_exp),log(mean_exp),log(mean_exp-sd_exp)-log(mean_exp)/)
  else
    unusual = .true.
  endif

  deallocate(bintemp)

else !all results are identical
     !or, almost all of the results are an upper or lower limit with a few
     !actual measurements
     !to account for this second possibility, take the value from the middle of
     !the array as the output quantity.  This will be a useful value in all
     !cases.
  uncertainty_array(1) = 0.D0
  uncertainty_array(2) = input_array(size(input_array)/2)
  uncertainty_array(3) = 0.D0
endif

end subroutine get_uncertainties

character(len=10) function gettime()
implicit none
character(len=10) :: time

  call DATE_AND_TIME(TIME=time)
  gettime = time(1:2)//":"//time(3:4)//":"//time(5:6)
  return

end function gettime

character(len=100) function latex_number(inputnumber)
! for any number in the form aEx, write it in latex format, ie a$\times$10$^{x}$
real(kind=dp) :: inputnumber, mantissa
integer :: exponent, pos

write (latex_number,"(ES14.6)") inputnumber
pos = INDEX (latex_number,'E')
read (latex_number(1:pos-1), '(F6.3)') mantissa
read (latex_number(pos+1:14), '(I3)') exponent

if (exponent .ge. -2 .and. exponent .le. 1) then
  write (latex_number,"(F6.2)") inputnumber !just print out normal number if it's between 0.01 & 10
elseif (exponent .ge. 2 .and. exponent .le. 4) then
  write (latex_number,"(I5)") 10*nint(inputnumber/10) ! write out integer rounded to nearest 10 if it's between 10 and 10,000
else !otherwise, write out a formatted exponent
  write (latex_number,"(F6.2,'\times 10^{',I3,'}')") mantissa,exponent
endif
return

end function latex_number
end program
