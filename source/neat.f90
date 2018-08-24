! NEAT, the nebular abundance analysis tool
! (C) 2006- Roger Wesson, Dave Stock, Peter Scicluna
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
        use mod_hydrogen
        use mod_helium
        use mod_oii_diagnostics
        use mod_weights

        implicit none

        character :: switch_ext !switch for extinction laws
        character :: switch_he  !switch for helium atomic data
        character :: switch_icf !switch for which ICF scheme to use
        integer :: I, j, runs, Narg !runs = number of runs for randomiser

!input options

        character(len=2048), dimension(:), allocatable :: options
        character(len=2048) :: commandline

!file reading variables

        logical :: file_exists
        logical :: identifylines
        type(line),dimension(:), allocatable :: linelist
        type(line),dimension(:), allocatable :: linelist_original
        type(line),dimension(:,:), allocatable :: all_linelists
        character(len=512) :: filename,configfile,defaultconfigfile
        character(len=1) :: blank
        integer :: listlength, ncols, errstat
        real(kind=dp) :: normalise

!results and result processing

        type(resultarray), dimension(:), allocatable :: all_results
        type(resultarray), dimension(1) :: iteration_result
        real(kind=dp), dimension(:), allocatable :: quantity_result
        real(kind=dp), dimension(:,:), allocatable :: resultprocessingarray
        character(len=40), dimension(:,:), allocatable :: resultprocessingtext

!atomic data

        character(len=10), dimension(24) :: ionlist !list of ion names
        integer :: iion !# of ions in Ilines
        integer :: maxlevs,maxtemps
        type(atomic_data),allocatable :: atomicdata(:)
        real(kind=dp), dimension(:,:,:), allocatable :: heidata

!line indexing arrays - these contain the location within the main array of the lines of interest, plus levels and zones for CELs

        integer, dimension(3:40) :: H_Balmer
        integer, dimension(4:39) :: H_Paschen
        integer, dimension(44) :: HeI_lines
        integer, dimension(20,2:6) :: HeII_lines
        type(cel), dimension(88) :: ILs !todo: work out why this becomes undefined on entry if its shape is assumed

!extinction

        logical :: calculate_extinction=.true.
        real(kind=dp) :: meanextinction, R

!RP effect switch

        logical :: norp=.true.

!subtract recombination contribution to auroral lines switch - 1=no subtraction, 2=subtraction

        integer :: subtract_recombination=1

!diagnostics

        type(weightingarray) :: weights

!binning and uncertainties

        real(kind=dp), dimension(3) :: uncertainty_array=0d0
        type(arraycount), dimension (:), allocatable :: binned_quantity_result
        logical :: unusual
        integer :: verbosity,nbins,nperbin

!OpenMP

        integer :: omp_get_num_threads

!diagnostic array

        real(kind=dp), dimension(6) :: diagnostic_array

!output formats

        character(len=35) :: extinction_format, diagnostic_format, diagnostic_ratio_format, abundances_format, adf_format

        configfile=""
        defaultconfigfile=trim(PREFIX)//"/share/neat/default.cfg"

        iion=24 ! number of ions for which we can calculate abundances
                ! todo: calculate this automatically
                ! and replace it in calls to subroutines with getting the array size within the subroutine

!multiple output formats defined as variables so they can be passed to
!the printout subroutine

        extinction_format = "(X,A,F9.3,SP,F9.3,F9.3,S)"
        diagnostic_format = "(X,A,F9.1,SP,F9.1,F9.1,S)"
        diagnostic_ratio_format = "(X,A,F9.3,SP,F9.3,F9.3,S)"
        abundances_format = "(X,A,ES14.3,SP,ES14.3,ES14.3,S)"
        adf_format = "(X,A,F8.2,SP,F8.2,F8.2,S)"

        print *,"NEAT, the Nebular Empirical Analysis Tool"
        print *,"version ",VERSION
        print *

        Narg = IARGC() !count input arguments

        if (Narg .eq. 0) then
           print *,"Syntax: neat [option1 value1] [option2 value2] .. [optionx valuex]"
           print *,"type  man neat  for help"
           stop
        endif

        call get_command(commandline)

        allocate (options(Narg))

        do i=1,Narg
                call getarg(i,options(i))
        enddo

! set defaults

        runs=1
        switch_ext="S" !Howarth 1983 Galactic law
        switch_he="P"  !Porter 2012 data
        switch_icf="D" !DI14 ICF
        filename=""
        meanextinction=0.D0
        diagnostic_array=0.D0
        verbosity=3
        R=3.1
        identifylines=.false.
        nbins=25
        normalise = 0.d0

! start the logging output to terminal

        print *,gettime(),"starting code"
        print *,gettime(),"command line: ",trim(commandline)

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

! read in the line list, allocate original linelist array

        call read_linelist(filename,linelist,listlength,ncols,errstat)

! error status was incremented by 1, 2 or 4 depending on error so btest function determines which errors were encountered
! test the two fatal errors first, then the warning.

        if (btest(errstat,0)) then
                print *,gettime(),"error: line list reading failed"
                print *,"            This can happen if it doesn't have three columns"
                call exit(1)
        elseif (btest(errstat,1)) then
                print*, gettime(),"cheese shop error - no inputs"
                call exit(1)
        elseif (btest(errstat,2) .and. runs .gt. 1) then !only relevant if uncertainties were requested
                print*, gettime(),"warning: no uncertainties given, arbitrarily assuming 10 per cent for everything"
        elseif (btest(errstat,3)) then
                print *,gettime(),"warning: duplicate line measurements detected. Only the second measurement is used in abundance determinations"
        else
                print "(X,A,A,A,A,I3,A)", gettime(),"line list file ",trim(filename)," read successfully (",listlength," lines)"
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
                if (blank .ne. "y" .and. blank .ne. "Y") then
                  print *,gettime()," : analysis cancelled."
                  call exit(1)
                endif
        endif

! normalise to Hbeta

        if (minval(abs(linelist(:)%wavelength - 4861.33)) .lt. 0.001) then
          normalise = 100.d0/linelist(minloc(abs(linelist(:)%wavelength - 4861.33),1))%intensity
        else
          print *,gettime(),"error: no H beta detected. No further analysis possible"
          call exit(1)
        endif

        linelist%intensity = linelist%intensity * normalise
        linelist%int_err = linelist%int_err * normalise
        linelist%blend_intensity = linelist%blend_intensity * normalise
        linelist%blend_int_err = linelist%blend_int_err * normalise

! check for and remove negative line fluxes and uncertainties

        if (minval(linelist%intensity) .lt. 0.) then
          where (linelist%intensity .lt. 0)
            linelist%intensity = 0.D0
          endwhere
        print *,gettime(),"warning: negative line fluxes set to zero"
        endif

        if (minval(linelist%int_err) .lt. 0.) then
          where (linelist%int_err .lt. 0)
            linelist%int_err = abs(linelist%int_err)
          endwhere
        print *,gettime(),"warning: negative flux uncertainties converted to positive"
        endif

        print *

! select reddening law, get flambda from functions if one of the five, otherwise read from file

        if (switch_ext == "S") then
                print *,gettime(),"using Howarth (1983) galactic extinction law"
        elseif (switch_ext == "H") then
                print *,gettime(),"using Howarth (1983) LMC extinction law"
        elseif (switch_ext == "C") then
                print *,gettime(),"using CCM (1989) galactic extinction law"
        elseif (switch_ext == "P") then
                print *,gettime(),"using Prevot et al. (1984) SMC extinction law"
        elseif (switch_ext == "F") then
                print *,gettime(),"using Fitzpatrick (1990) galactic extinction law"
        endif

        call get_flambda(linelist, switch_ext, R)

        if (switch_he == "S") then
                print *,gettime(),"using Smits (1996) He I emissivities"
        else
                print *,gettime(),"using Porter et al. (2012) He I emissivities"
        endif

        if (switch_ICF == "K") then
                print *,gettime(),"using Kingsburgh & Barlow (1994) ICF"
        elseif (switch_ICF == "D") then
                print *,gettime(),"using Delgado-Inglada et al. (2014) ICF"
        else
                print *,gettime(),"using Peimbert et al. (1992) ICF"
        endif

        print *

! read the CEL data

        call read_celdata(ILs,ionlist)

! do line identifying here

        !first find H & He lines, and CELs
        call get_H(H_Balmer, H_Paschen, linelist)
        call get_Hei(Hei_lines, linelist)
        call get_Heii(Heii_lines, linelist)
        call get_cels(ILs,linelist)

! read in the weights to be used in the analysis.

        print *,gettime(),"reading in abundance analysis weights from ",trim(defaultconfigfile)
        call setweights(defaultconfigfile,weights,linelist,ILs,H_Balmer,H_Paschen,HeII_lines)
        if (configfile .ne. "") then
          print *,gettime(),"reading in abundance analysis weights from ",trim(configfile)
          print *,gettime(),"(defaults still apply for any weights not specified in this file)"
          call setweights(configfile,weights,linelist,ILs,H_Balmer,H_Paschen,HeII_lines)
        endif

!CELs read, can now read in atomic data

        print *,gettime(),"reading atomic data from ",trim(PREFIX),"/share/neat"
        allocate(atomicdata(iion))
        do I=1,iion
            atomicdata(I)%ion = ionlist(I)
            call read_atomic_data(atomicdata(I))
        enddo

!read ORL data
        call read_orl_data

!read hydrogen emissivities

        print *,gettime(),"reading H emissivities"
        call read_hydrogen

!read helium emissivities

        print *,gettime(),"reading He I emissivities"
        if (switch_he .eq. "P") then
          allocate(heidata(21,14,44))
          heidata = 0.D0
          call read_porter(heidata)
        elseif (switch_he .eq. "S") then
          allocate(heidata(3,6,44))
          heidata = 0.D0
          call read_smits(heidata)
        endif

        print *,gettime(),"reading He II emissivities"
        call read_heii

!read oii data for diagnostics

        print *,gettime(),"reading OII diagnostics data"
        call read_oii_diagnosticdata

!find maximum #levels and temperatures - pass to equib to reduce footprint

        maxlevs = maxval(atomicdata%nlevs)
        maxtemps = maxval(atomicdata%ntemps)

!now check number of iterations.  If 1, line list is fine as is.  If more than one, randomize the fluxes

        allocate(all_results(runs))

        if(runs == 1)then !calculates abundances without uncertainties
                print *
                print *,gettime(),"doing abundance calculations"
                call abundances(linelist, listlength, iteration_result, meanextinction, calculate_extinction, ILs, diagnostic_array,iion,atomicdata,maxlevs,maxtemps, heidata, switch_he, switch_icf, H_Balmer, H_Paschen, HeI_lines, HeII_lines, weights,subtract_recombination)
                all_results(1)=iteration_result(1) ! copy the iteration result to all_results to simplify the writing out of results later
                print *,gettime(),"finished abundance calculations"
        else if(runs .gt. 1)then

!save unrandomised line list

                linelist_original = linelist

                call init_random_seed()!sets seed for randomiser

                !allocate array to store all line info

                allocate(all_linelists(size(linelist),runs))

                !main loop
!$OMP PARALLEL default(firstprivate) shared(all_linelists,listlength,norp,calculate_extinction,ILs,diagnostic_array,iion,atomicdata,maxlevs,maxtemps,switch_he,switch_icf,all_results,subtract_recombination)
!$OMP MASTER

                print *
                if (omp_get_num_threads().gt.1) then
                  print "(X,A10,X,A,I2,A)",gettime(),"starting Monte Carlo calculations, using ",omp_get_num_threads()," processors"
                else
                  print "(X,A10,X,A)",gettime(),"starting Monte Carlo calculations"
                endif
                print *,gettime(),"completed ",0,"%"
!$OMP END MASTER
!$OMP DO schedule(dynamic)
                do I=1,runs

                         if ( (10.0*dble(i)/dble(runs)) == int(10*i/runs) ) print *,gettime(),"completed ",100*i/runs,"%"
!                        print*, "iteration ", i, "of", runs

                        call randomizer(linelist, listlength, norp)
                        call abundances(linelist, listlength, iteration_result, meanextinction, calculate_extinction, ILs, diagnostic_array,iion,atomicdata,maxlevs,maxtemps, heidata, switch_he, switch_icf, H_Balmer, H_Paschen, HeI_lines, HeII_lines, weights, subtract_recombination)

                        !store all line and derived quantity in arrays
                        all_linelists(:,i)=linelist
                        all_results(i)=iteration_result(1)
                        !restore the unrandomized line list ready for the next iteration
                        linelist = linelist_original
                enddo
!$OMP END DO
!$OMP END PARALLEL
        else
                print*, gettime(),"error: I didn't want to be a barber anyway. I wanted to be... a lumberjack!   Also, a positive number of runs helps.."
                call exit(1)
        endif

!now write all the lines to line list files, plain text and latex
!linelist first

        print *
        print *,gettime(),"Writing line list"

        allocate(quantity_result(runs))
        quantity_result=0d0
        open (650,file=trim(filename)//"_linelist", status='replace', access='sequential', action='write')
        open (651,file=trim(filename)//"_linelist.tex", status='replace', access='sequential', action='write')

        write (650,*) "Lambda  Ion        F(line) I(line)      Abundance"
        write (651,*) "\begin{longtable}{lrlrlllllll}"
        write (651,*) "\hline"
        write (651,*) "$ \lambda $ & Ion & $F \left( \lambda \right) $ && $I \left( \lambda \right) $ & Ion & Multiplet & Lower term & Upper term & g$_1$ & g$_2$ \\"
        write (651,*) "\hline"

        if (runs .gt. 1) then
                do j=1, listlength

!observed wavelength if known

                if (ncols .ge. 4) then
                  if (linelist(j)%wavelength_observed .gt. 0.) then
                    write (650,"(X,F8.2,X)", advance='no') all_linelists(j,1)%wavelength_observed
                    write (651,"(X,F8.2,' & ')", advance='no') all_linelists(j,1)%wavelength_observed
                  else
                    write (650,"(X,A,X)", advance='no') "       *"
                    write (651,"(X,A,' & ')", advance='no') "       *"
                  endif
                endif

!rest wavelength, ion name for plain text file

                write (650,"(X,F8.2,X,A11)", advance='no') all_linelists(j,1)%wavelength,all_linelists(j,1)%name
                write (651,"(X,F8.2,' & ',A15,' & ')", advance='no') all_linelists(j,1)%wavelength

!line flux

               if (all_linelists(j,1)%intensity .eq. 0.d0 .and. all_linelists(j,1)%blend_intensity .eq. 0.d0) then
                 write (650,"(A)", advance='no') " * "
                 write (651,"(A)", advance='no') " *     &             &"
               else
                 write (650,"(F8.3,A,F8.3,3X)", advance='no') all_linelists(j,1)%intensity," +-",all_linelists(j,1)%int_err
                 write (651,"(F8.3,'& $\pm$',F8.3, '&')", advance='no') all_linelists(j,1)%intensity,all_linelists(j,1)%int_err
               endif

!dereddened flux
                quantity_result = all_linelists(j,:)%int_dered

                call get_uncertainties(quantity_result, binned_quantity_result, uncertainty_array, unusual,nbins)

                if (all_linelists(j,1)%intensity .ne. 0.d0 .or. all_linelists(j,1)%blend_intensity .ne. 0.d0) then
                  if (uncertainty_array(1) .ne. uncertainty_array(3)) then
                    write (650,"(F8.3,SP,F8.3,SP,F8.3)", advance='no') uncertainty_array(2),uncertainty_array(1),-uncertainty_array(3)
                    write (651,"(F8.3,'& $^{',SP,F8.3,'}_{',SP,F8.3,'}$')", advance='no') uncertainty_array(2),uncertainty_array(1),-uncertainty_array(3)
                  else
                    write (650,"(F8.3,A,F8.3,4X)", advance='no') uncertainty_array(2)," +-",uncertainty_array(1)
                    write (651,"(F8.3,'& $\pm$',F8.3)", advance='no') uncertainty_array(2),uncertainty_array(1)
                  endif
                else
                  write (650,"(A)", advance='no') " * "
                  write (651,"(A)", advance='no') " *     &        "
                endif

! transition data

!                write (650,*)
                write (651,*) all_linelists(j,1)%linedata, "\\"

!abundance - write out if there is an abundance for the line, don't write
!anything except a line break if there is no abundance for the line.
!todo: add an option to choose whether or not to put abundances in the line list table

                quantity_result = all_linelists(j,:)%abundance
                call get_uncertainties(quantity_result, binned_quantity_result, uncertainty_array, unusual,nbins)
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

                enddo
        else ! runs == 1, no uncertainties to write out

                do i=1,listlength
                  if (ncols .ge. 4) then
                    if (linelist(i)%wavelength_observed .gt. 0.) then
                      write (650,"(X,F8.2,X)", advance='no') linelist(i)%wavelength_observed
                      write (651,"(X,F8.2,' & ')", advance='no') linelist(i)%wavelength_observed
                    else
                      write (650,"(X,A,X)", advance='no') "       *"
                      write (651,"(X,A,' & ')", advance='no') "       *"
                    endif
                  endif
                  if (linelist(i)%intensity .ne. 0.d0) then
                    if (linelist(i)%abundance .gt. 0.0) then
                      write (650,"(X,F8.2,X,A11,F8.3,X,F8.3,X,ES14.3)") linelist(i)%wavelength,linelist(i)%name,linelist(i)%intensity,linelist(i)%int_dered, linelist(i)%abundance
                      write (651,"(X,F8.2,X,'&',X,F8.3,X,'&',X,F8.3,X,A,'\\')") linelist(i)%wavelength,linelist(i)%intensity,linelist(i)%int_dered,linelist(i)%linedata
                    else
                      write (650,"(X,F8.2,X,A11,F8.3,X,F8.3)") linelist(i)%wavelength,linelist(i)%name,linelist(i)%intensity,linelist(i)%int_dered
                      write (651,"(X,F8.2,X,'&',X,F8.3,X,'&',X,F8.3,X,A,'\\')") linelist(i)%wavelength,linelist(i)%intensity,linelist(i)%int_dered,linelist(i)%linedata
                    endif
                  else
                    write (650,"(X,F8.2,X,A11,A7,X,A7)") linelist(i)%wavelength,linelist(i)%name,"*      ","*      "
                    write (651,"(X,F8.2,X,'&',X,A7,X,'&',X,A7,X,A,'\\')") linelist(i)%wavelength,"*      ","*      ",linelist(i)%linedata
                  endif
                enddo

        endif

        write (651,*) "\hline"
        !write (651,*) "\caption{}"
        write (651,*) "\label{tab:",trim(filename)//"_linelist}"
        write (651,*) "\end{longtable}"

        close(650)
        close(651)

        print *,gettime(),"Linelist written to files:"
        print *,"              ",trim(filename),"_linelist"
        print *,"              ",trim(filename),"_linelist.tex"

!now write out the summary files and all the binned data

        if (verbosity .eq. 1) then
          print *,gettime(),"Writing summary files, binned results and complete results"
        elseif (verbosity .eq. 2) then
         print *,gettime(),"Writing summary files and binned results"
        else
          print *,gettime(),"Writing summary files"
        endif

!first, define arrays with links to all the data that needs processing.
!extinction, diagnostics, cel abundances, orl abundances, strong line
!abundances, adfs

        allocate(resultprocessingarray(164,runs))
        resultprocessingarray=0d0
        allocate(resultprocessingtext(164,4))

!extinction

        resultprocessingarray(1,:) = all_results%cHb_ha
        resultprocessingtext(1,:) = (/"c(Hb) (Ha/Hb)                      ","c(H$\beta)$ (H$\alpha/H$\beta)     ", extinction_format, "chb_ha_hb                          "/)
        resultprocessingarray(2,:) = all_results%cHb_hg
        resultprocessingtext(2,:) = (/"c(Hb) (Hg/Hb)                      ","c(H$\beta)$ (H$\gamma/H$\beta)     ", extinction_format, "chb_hg_hb                          "/)
        resultprocessingarray(3,:) = all_results%cHb_hd
        resultprocessingtext(3,:) = (/"c(Hb) (Hd/Hb)                      ","c(H$\beta)$ (H$\delta/H$\beta)     ", extinction_format, "chb_hd_hb                          "/)
        resultprocessingarray(4,:) = all_results%mean_cHb
        resultprocessingtext(4,:) = (/"mean c(Hb)                         ","c(H$\beta)$                        ", extinction_format, "mean_chb                           "/)

!diagnostics
!low density

        resultprocessingarray(5,:) = all_results%oii_density
        resultprocessingtext(5,:) = (/"[OII] density                      ","{}[O~{\sc ii}] density             ", diagnostic_format, "density_oii                        "/)
        resultprocessingarray(6,:) = all_results%oii_density_ratio
        resultprocessingtext(6,:) = (/"[OII] 3729/3726 ratio              ","{}[O~{\sc ii}] 3729/3726 ratio     ", diagnostic_ratio_format, "density_oii_ratio                  "/)
        resultprocessingarray(7,:) = all_results%SII_density
        resultprocessingtext(7,:) = (/"[SII] density                      ","{}[S~{\sc ii}] density             ", diagnostic_format, "density_sii                        "/)
        resultprocessingarray(8,:) = all_results%sii_density_ratio
        resultprocessingtext(8,:) = (/"[SII] 6731/6717 ratio              ","{}[S~{\sc ii}] 6731/6717 ratio     ", diagnostic_ratio_format, "density_sii_ratio                  "/)
        resultprocessingarray(9,:) = all_results%low_density
        resultprocessingtext(9,:) = (/"Low ionisation density             ","Low ionisation density             ", diagnostic_format, "density_low                        "/)

!low temperature

        resultprocessingarray(10,:) = all_results%oii_temp
        resultprocessingtext(10,:) = (/"[OII] temperature                  ","{}[O~{\sc ii}] temperature         ", diagnostic_format, "temp_oii                           "/)
        resultprocessingarray(11,:) = all_results%oii_temp_ratio
        resultprocessingtext(11,:) = (/"[OII] (7320+7330)/(3726+3729) ratio","{}[O~{\sc ii}] 7320+7330/3726+3729 ", diagnostic_ratio_format, "temp_oii_ratio                     "/)
        resultprocessingarray(12,:) = all_results%SII_temp
        resultprocessingtext(12,:) = (/"[SII] temperature                  ","{}[S~{\sc ii}] temperature         ", diagnostic_format, "temp_sii                           "/)
        resultprocessingarray(13,:) = all_results%sii_temp_ratio
        resultprocessingtext(13,:) = (/"[SII] (6717+6731)/(4068+4076) ratio","{}[S~{\sc ii}] 6717+6731/4068+4076 ", diagnostic_ratio_format, "temp_sii_ratio                     "/)
        resultprocessingarray(14,:) = all_results%NII_temp
        resultprocessingtext(14,:) = (/"[NII] temperature                  ","{}[N~{\sc ii}] temperature         ", diagnostic_format, "temp_nii                           "/)
        resultprocessingarray(15,:) = all_results%nii_temp_ratio
        resultprocessingtext(15,:) = (/"[NII] (6548+6584)/5754 ratio       ","{}[N~{\sc ii}] 6548+6584/5754 ratio", diagnostic_ratio_format, "temp_nii_ratio                     "/)
        resultprocessingarray(16,:) = all_results%OI_temp
        resultprocessingtext(16,:) = (/"[OI] temperature                   ","{}[O~{\sc i}] temperature          ", diagnostic_format, "temp_oi                            "/)
        resultprocessingarray(17,:) = all_results%oi_temp_ratio
        resultprocessingtext(17,:) = (/"[OI] (6300+6363)/5577 ratio        ","{}[O~{\sc i}] 6300+6363/5577 ratio ", diagnostic_ratio_format, "temp_oi_ratio                      "/)
        resultprocessingarray(18,:) = all_results%CI_temp
        resultprocessingtext(18,:) = (/"[CI] temperature                   ","{}[C~{\sc i}] temperature          ", diagnostic_format, "temp_ci                            "/)
        resultprocessingarray(19,:) = all_results%ci_temp_ratio
        resultprocessingtext(19,:) = (/"[CI] (9824+9850)/8727 ratio        ","{}[C~{\sc i}] 9824+9850/8727 ratio ", diagnostic_ratio_format, "temp_ci_ratio                      "/)
        resultprocessingarray(20,:) = all_results%low_temp
        resultprocessingtext(20,:) = (/"Low ionisation temperature         ","Low ionisation temperature         ", diagnostic_format, "temp_low                           "/)


!medium density

        resultprocessingarray(21,:) = all_results%cliii_density
        resultprocessingtext(21,:) = (/"[ClIII] density                    ","{}[Cl~{\sc iii}] density           ", diagnostic_format, "density_cliii                      "/)
        resultprocessingarray(22,:) = all_results%cliii_density_ratio
        resultprocessingtext(22,:) = (/"[ClIII] 5537/5517 ratio            ","{}[Cl~{\sc iii}] 5537/5517 ratio   ", diagnostic_ratio_format, "density_cliii_ratio                "/)
        resultprocessingarray(23,:) = all_results%ArIV_density
        resultprocessingtext(23,:) = (/"[ArIV] density                     ","{}[Ar~{\sc iv}] density            ", diagnostic_format, "density_ariv                       "/)
        resultprocessingarray(24,:) = all_results%ariv_density_ratio
        resultprocessingtext(24,:) = (/"[ArIV] 4740/4711 ratio             ","{}[Ar~{\sc iv}] 4740/4711 ratio    ", diagnostic_ratio_format, "density_sii_ratio                  "/)
        resultprocessingarray(25,:) = all_results%CIII_density
        resultprocessingtext(25,:) = (/"[CIII] density                     ","{}[C~{\sc iii}] density            ", diagnostic_format, "density_ciii                       "/)
        resultprocessingarray(26,:) = all_results%ciii_density_ratio
        resultprocessingtext(26,:) = (/"[CIII] 1909/1907 ratio             ","{}[C~{\sc iii}] 1909/1907 ratio    ", diagnostic_ratio_format, "density_ciii_ratio                 "/)
        resultprocessingarray(27,:) = all_results%OIII_IR_density
        resultprocessingtext(27,:) = (/"[OIII] IR density                  ","{}[O~{\sc iii}] IR density         ", diagnostic_format, "density_oiii_ir                    "/)
        resultprocessingarray(28,:) = all_results%oiii_ir_density_ratio
        resultprocessingtext(28,:) = (/"[OIII] IR 88um/52um ratio          ","{}[O~{\sc iii}] IR 88um/52um ratio ", diagnostic_ratio_format, "density_oiii_ir_ratio              "/)
        resultprocessingarray(29,:) = all_results%SIII_IR_density
        resultprocessingtext(29,:) = (/"[SIII] IR density                  ","{}[S~{\sc iii}] IR density         ", diagnostic_format, "density_siii_ir                    "/)
        resultprocessingarray(30,:) = all_results%siii_ir_density_ratio
        resultprocessingtext(30,:) = (/"[SIII] IR 33um/18um ratio          ","{}[S~{\sc iii}] IR 33um/18um ratio ", diagnostic_ratio_format, "density_siii_ir_ratio              "/)
        resultprocessingarray(31,:) = all_results%ArIII_IR_density
        resultprocessingtext(31,:) = (/"[ArIII] IR density                 ","{}[Ar~{\sc iii}] IR density        ", diagnostic_format, "density_ariii_ir                   "/)
        resultprocessingarray(32,:) = all_results%ariii_ir_density_ratio
        resultprocessingtext(32,:) = (/"[ArIII] 9um/22um ratio             ","{}[Ar~{\sc iii}] IR 9um/22um ratio ", diagnostic_ratio_format, "density_ariii_ir_ratio             "/)
        resultprocessingarray(33,:) = all_results%NeIII_IR_density
        resultprocessingtext(33,:) = (/"[NeIII] IR density                 ","{}[Ne~{\sc iii}] IR density        ", diagnostic_format, "density_neiii_ir                   "/)
        resultprocessingarray(34,:) = all_results%neiii_ir_density_ratio
        resultprocessingtext(34,:) = (/"[Ne III] 15um/36um ratio           ","{}[Ne~{\sc iii}] IR 15um/36um ratio", diagnostic_ratio_format, "density_neiii_ir_ratio             "/)
        resultprocessingarray(35,:) = all_results%med_density
        resultprocessingtext(35,:) = (/"Medium ionisation density          ","Medium ionisation density          ", diagnostic_format, "density_med                        "/)

! medium temperature

        resultprocessingarray(36,:) = all_results%OIII_temp
        resultprocessingtext(36,:) = (/"[OIII] temperature                 ","{}[O~{\sc iii}] temperature        ", diagnostic_format, "temp_oiii                          "/)
        resultprocessingarray(37,:) = all_results%oiii_temp_ratio
        resultprocessingtext(37,:) = (/"[OIII] (4959+5007)/4363 ratio      ","{}[O~{\sc iii}] 4959+5007/4363     ", diagnostic_ratio_format, "temp_oiii_ratio                    "/)
        resultprocessingarray(38,:) = all_results%OIII_IR_temp
        resultprocessingtext(38,:) = (/"[OIII] IR temperature              ","{}[O~{\sc iii}] IR temperature     ", diagnostic_format, "temp_oiii_ir                       "/)
        resultprocessingarray(39,:) = all_results%oiii_ir_temp_ratio
        resultprocessingtext(39,:) = (/"[OIII] (4959+5007)/52um ratio      ","{}[O~{\sc iii}] 4959+5007/52um     ", diagnostic_ratio_format, "temp_oiii_ir_ratio                 "/)
        resultprocessingarray(40,:) = all_results%OIII_UV_temp
        resultprocessingtext(40,:) = (/"[OIII] UV temperature              ","{}[O~{\sc iii}] UV temperature     ", diagnostic_format, "temp_oiii_uv                       "/)
        resultprocessingarray(41,:) = all_results%oiii_uv_temp_ratio
        resultprocessingtext(41,:) = (/"[OIII] (4959+5007)/1666 ratio      ","{}[O~{\sc iii}] 4959+5007/1666     ", diagnostic_ratio_format, "temp_oiii_uv_ratio                 "/)
        resultprocessingarray(42,:) = all_results%NeIII_temp
        resultprocessingtext(42,:) = (/"[NeIII] temperature                ","{}[Ne~{\sc iii}] temperature       ", diagnostic_format, "temp_neiii                         "/)
        resultprocessingarray(43,:) = all_results%neiii_temp_ratio
        resultprocessingtext(43,:) = (/"[NeIII] (3868+3967)/3342 ratio     ","{}[Ne~{\sc iii}] 3868+3967/3342    ", diagnostic_ratio_format, "temp_neiii_ratio                   "/)
        resultprocessingarray(44,:) = all_results%NeIII_IR_temp
        resultprocessingtext(44,:) = (/"[NeIII] IR temperature             ","{}[Ne~{\sc iii}] IR temperature    ", diagnostic_format, "temp_neiii_ir                      "/)
        resultprocessingarray(45,:) = all_results%neiii_ir_temp_ratio
        resultprocessingtext(45,:) = (/"[NeIII] (3868+3967)/15um ratio     ","{}[Ne~{\sc iii}] 3868+3967/15um    ", diagnostic_ratio_format, "temp_neiii_ir_ratio                "/)
        resultprocessingarray(46,:) = all_results%ArIII_temp
        resultprocessingtext(46,:) = (/"[ArIII] temperature                ","{}[Ar~{\sc iii}] temperature       ", diagnostic_format, "temp_ariii                         "/)
        resultprocessingarray(47,:) = all_results%ariii_temp_ratio
        resultprocessingtext(47,:) = (/"[ArIII] (7135+7751)/5192 ratio     ","{}[Ar~{\sc iii}] 7135+7751/5192    ", diagnostic_ratio_format, "temp_ariii_ratio                   "/)
        resultprocessingarray(48,:) = all_results%SIII_temp
        resultprocessingtext(48,:) = (/"[SIII] temperature                 ","{}[S~{\sc iii}] temperature        ", diagnostic_format, "temp_siii                          "/)
        resultprocessingarray(49,:) = all_results%siii_temp_ratio
        resultprocessingtext(49,:) = (/"[SIII] (9069+9531)/6312 ratio      ","{}[S~{\sc iii}] 9069+9531/6312     ", diagnostic_ratio_format, "temp_siii_ratio                    "/)
        resultprocessingarray(50,:) = all_results%med_temp
        resultprocessingtext(50,:) = (/"Medium ionisation temperature      ","Medium ionisation temperature      ", diagnostic_format, "temp_med                           "/)

!high density

        resultprocessingarray(51,:) = all_results%neiv_density
        resultprocessingtext(51,:) = (/"[NeIV] density                     ","{}[Ne~{\sc iv}] density            ", diagnostic_format, "density_neiv                       "/)
        resultprocessingarray(52,:) = all_results%neiv_density_ratio
        resultprocessingtext(52,:) = (/"[NeIV] 2425/2423 ratio             ","{}[Ne~{\sc iv}] 2425/2423 ratio    ", diagnostic_ratio_format, "density_neiv_ratio                 "/)
        resultprocessingarray(53,:) = all_results%high_density
        resultprocessingtext(53,:) = (/"High ionisation density            ","High ionisation density            ", diagnostic_format, "density_high                       "/)

!high temperature

        resultprocessingarray(54,:) = all_results%ArV_temp
        resultprocessingtext(54,:) = (/"[ArV] temperature                  ","{}[Ar~{\sc v}] temperature         ", diagnostic_format, "temp_arv                           "/)
        resultprocessingarray(55,:) = all_results%arv_temp_ratio
        resultprocessingtext(55,:) = (/"[ArV] (6435+7005)/4625 ratio       ","{}[Ar~{\sc v}] 6435+7005/4625 ratio", diagnostic_ratio_format, "temp_arv_ratio                     "/)
        resultprocessingarray(56,:) = all_results%NeV_temp
        resultprocessingtext(56,:) = (/"[NeV] temperature                  ","{}[Ne~{\sc v}] temperature         ", diagnostic_format, "temp_nev                           "/)
        resultprocessingarray(57,:) = all_results%nev_temp_ratio
        resultprocessingtext(57,:) = (/"[NeV] (3426+3345)/2975 ratio       ","{}[Ne~{\sc v}] 3426+3345/2975 ratio", diagnostic_ratio_format, "temp_nev_ratio                     "/)
        resultprocessingarray(58,:) = all_results%high_temp
        resultprocessingtext(58,:) = (/"High temperature                   ","High temperature                   ", diagnostic_format, "temp_high                          "/)

!Recombination line diagnostics

        resultprocessingarray(59,:) = all_results%Balmer_jump_temp
        resultprocessingtext(59,:) = (/"Balmer jump temperature            ","Balmer jump temperature            ", diagnostic_format, "temp_balmerjump                    "/)

        resultprocessingarray(60,:) = all_results%Paschen_jump_temp
        resultprocessingtext(60,:) = (/"Paschen jump temperature           ","Paschen jump temperature           ", diagnostic_format, "temp_paschenjump                   "/)

        resultprocessingarray(61,:) = all_results%Balmerdec_density
        resultprocessingtext(61,:) = (/"Balmer decrement density           ","Balmer decrement density           ", diagnostic_format, "density_balmerdec                  "/)

        resultprocessingarray(62,:) = all_results%Paschendec_density
        resultprocessingtext(62,:) = (/"Paschen decrement density          ","Paschen decrement density          ", diagnostic_format, "density_paschendec                 "/)

        resultprocessingarray(63,:) = all_results%te_5876_4471
        resultprocessingtext(63,:) = (/"He I temperature (5876/4471)       ","He I temperature (5876/4471)       ", diagnostic_format, "temp_he_5876_4471                  "/)

        resultprocessingarray(64,:) = all_results%ratio_5876_4471
        resultprocessingtext(64,:) = (/"He I 5876/4471 ratio               ","He I 5876/4471 ratio               ", diagnostic_ratio_format, "temp_he_5876_4471_ratio            "/)

        resultprocessingarray(65,:) = all_results%te_6678_4471
        resultprocessingtext(65,:) = (/"He I temperature (6678/4471)       ","He I temperature (6678/4471)       ", diagnostic_format, "temp_he_6678_4471                  "/)

        resultprocessingarray(66,:) = all_results%ratio_6678_4471
        resultprocessingtext(66,:) = (/"He I 6678/4471 ratio               ","He I 6678/4471 ratio               ", diagnostic_ratio_format, "temp_he_6678_4471_ratio            "/)

        resultprocessingarray(67,:) = all_results%oii_rl_R1
        resultprocessingtext(67,:) = (/"OII 4649/4089 ratio                ","OII 4649/4089 ratio                ", diagnostic_ratio_format, "temp_oiirls_R1                     "/)

        resultprocessingarray(68,:) = all_results%oii_rl_R2
        resultprocessingtext(68,:) = (/"OII 4649/4662 ratio                ","OII 4649/4662 ratio                ", diagnostic_ratio_format, "temp_oiirls_R2                     "/)

        resultprocessingarray(69,:) = all_results%oii_te
        resultprocessingtext(69,:) = (/"OII temperature                    ","OII temperature                    ", diagnostic_format, "temp_oiirls                        "/)

        resultprocessingarray(70,:) = all_results%oii_ne
        resultprocessingtext(70,:) = (/"OII density                        ","OII density                        ", diagnostic_format, "density_oiirls                     "/)

! recombination contribution to CELs

        resultprocessingarray(71,:) = all_results%nii5754recCEL
        resultprocessingtext(71,:) = (/"NII 5754 assuming N2+ from CELs    ","N~{\sc ii}5754$_R$ (N$^{2+}$ CELs) ", adf_format, "recombination_5754_CELabund        "/)
        resultprocessingarray(72,:) = all_results%nii5754recRL
        resultprocessingtext(72,:) = (/"NII 5754 assuming N2+ from RLs     ","N~{\sc ii}5754$_R$ (N$^{2+}$ CELs) ", adf_format, "recombination_5754_RLabund         "/)
        resultprocessingarray(73,:) = all_results%oii7325recCEL
        resultprocessingtext(73,:) = (/"OII 7320,30 assuming O2+ from CELs ","O~{\sc ii}7325$_R$ (O$^{2+}$ CELs) ", adf_format, "recombination_7325_CELabund        "/)
        resultprocessingarray(74,:) = all_results%oii7325recRL
        resultprocessingtext(74,:) = (/"OII 7320,30 assuming O2+ from RLs  ","O~{\sc ii}7325$_R$ (O$^{2+}$ CELs) ", adf_format, "recombination_7325_RLabund         "/)
        resultprocessingarray(75,:) = all_results%oiii4363recCEL
        resultprocessingtext(75,:) = (/"OIII 4363 assuming O3+ from CELs   ","O~{\sc ii}4363$_R$ (O$^{3+}$ CELs) ", adf_format, "recombination_4363_CELabund        "/)
        resultprocessingarray(76,:) = all_results%oiii4363recRL
        resultprocessingtext(76,:) = (/"OIII 4363 assuming O3+ from RLs    ","O~{\sc ii}4363$_R$ (O$^{3+}$ CELs) ", adf_format, "recombination_4363_RLabund         "/)

!CEL abundances

        resultprocessingarray(77,:) = all_results%NC_abund_CEL
        resultprocessingtext(77,:) = (/"C0/H                               ","C$^{0}$/H                          ", abundances_format, "abund_cel_nc                       "/)
        resultprocessingarray(78,:) = all_results%cii_abund_CEL
        resultprocessingtext(78,:) = (/"C+/H                               ","C$^{+}$/H                          ", abundances_format, "abund_cel_cii                      "/)
        resultprocessingarray(79,:) = all_results%ciii_abund_CEL
        resultprocessingtext(79,:) = (/"C2+/H                              ","C$^{2+}$/H                         ", abundances_format, "abund_cel_ciii                     "/)
        resultprocessingarray(80,:) = all_results%civ_abund_CEL
        resultprocessingtext(80,:) = (/"C3+/H                              ","C$^{3+}$/H                         ", abundances_format, "abund_cel_civ                      "/)
        resultprocessingarray(81,:) = all_results%c_icf_CEL
        resultprocessingtext(81,:) = (/"icf(C)                             ","icf(C)                             ", abundances_format, "icf_cel_c                          "/)
        resultprocessingarray(82,:) = all_results%C_abund_CEL
        resultprocessingtext(82,:) = (/"C/H                                ","C$^{}$/H                           ", abundances_format, "abund_cel_c                        "/)
        resultprocessingarray(83,:) = all_results%nii_abund_CEL
        resultprocessingtext(83,:) = (/"N+/H                               ","N$^{+}$/H                          ", abundances_format, "abund_cel_nii                      "/)
        resultprocessingarray(84,:) = all_results%niii_abund_CEL
        resultprocessingtext(84,:) = (/"N2+/H                              ","N$^{2+}$/H                         ", abundances_format, "abund_cel_niii                     "/)
        resultprocessingarray(85,:) = all_results%niv_abund_CEL
        resultprocessingtext(85,:) = (/"N3+/H                              ","N$^{3+}$/H                         ", abundances_format, "abund_cel_niv                      "/)
        resultprocessingarray(86,:) = all_results%nv_abund_CEL
        resultprocessingtext(86,:) = (/"N4+/H                              ","N$^{4+}$/H                         ", abundances_format, "abund_cel_nv                       "/)
        resultprocessingarray(87,:) = all_results%n_icf_CEL
        resultprocessingtext(87,:) = (/"icf(N)                             ","icf(N)                             ", abundances_format, "icf_cel_n                          "/)
        resultprocessingarray(88,:) = all_results%N_abund_CEL
        resultprocessingtext(88,:) = (/"N/H                                ","N$^{}$/H                           ", abundances_format, "abund_cel_n                        "/)
        resultprocessingarray(89,:) = all_results%NO_abund_CEL
        resultprocessingtext(89,:) = (/"O0/H                               ","O$^{0}$/H                          ", abundances_format, "abund_cel_no                       "/)
        resultprocessingarray(90,:) = all_results%Oii_abund_CEL
        resultprocessingtext(90,:) = (/"O+/H                               ","O$^{+}$/H                          ", abundances_format, "abund_cel_oii                      "/)
        resultprocessingarray(91,:) = all_results%Oiii_abund_CEL
        resultprocessingtext(91,:) = (/"O2+/H                              ","O$^{2+}$/H                         ", abundances_format, "abund_cel_oiii                     "/)
        resultprocessingarray(92,:) = all_results%Oiv_abund_CEL
        resultprocessingtext(92,:) = (/"O3+/H                              ","O$^{3+}$/H                         ", abundances_format, "abund_cel_oiv                      "/)
        resultprocessingarray(93,:) = all_results%o_icf_CEL
        resultprocessingtext(93,:) = (/"icf(O)                             ","icf(O)                             ", abundances_format, "icf_cel_o                          "/)
        resultprocessingarray(94,:) = all_results%O_abund_CEL
        resultprocessingtext(94,:) = (/"O/H                                ","O$^{}$/H                           ", abundances_format, "abund_cel_o                        "/)
        resultprocessingarray(95,:) = all_results%NeII_abund_CEL
        resultprocessingtext(95,:) = (/"Ne+/H                              ","Ne$^{+}$/H                         ", abundances_format, "abund_cel_neii                     "/)
        resultprocessingarray(96,:) = all_results%NeIII_abund_CEL
        resultprocessingtext(96,:) = (/"Ne2+/H                             ","Ne$^{2+}$/H                        ", abundances_format, "abund_cel_neiii                    "/)
        resultprocessingarray(97,:) = all_results%NeIV_abund_CEL
        resultprocessingtext(97,:) = (/"Ne3+/H                             ","Ne$^{3+}$/H                        ", abundances_format, "abund_cel_neiv                     "/)
        resultprocessingarray(98,:) = all_results%NeV_abund_CEL
        resultprocessingtext(98,:) = (/"Ne4+/H                             ","Ne$^{4+}$/H                        ", abundances_format, "abund_cel_nev                      "/)
        resultprocessingarray(99,:) = all_results%ne_icf_CEL
        resultprocessingtext(99,:) = (/"icf(Ne)                            ","icf(Ne)                            ", abundances_format, "icf_cel_ne                         "/)
        resultprocessingarray(100,:) = all_results%Ne_abund_CEL
        resultprocessingtext(100,:) = (/"Ne/H                               ","Ne$^{}$/H                          ", abundances_format, "abund_cel_ne                       "/)
        resultprocessingarray(101,:) = all_results%ArIII_abund_CEL
        resultprocessingtext(101,:) = (/"Ar2+/H                             ","Ar$^{2+}$/H                        ", abundances_format, "abund_cel_ariii                    "/)
        resultprocessingarray(102,:) = all_results%ArIV_abund_CEL
        resultprocessingtext(102,:) = (/"Ar3+/H                             ","Ar$^{3+}$/H                        ", abundances_format, "abund_cel_ariv                     "/)
        resultprocessingarray(103,:) = all_results%ArV_abund_CEL
        resultprocessingtext(103,:) = (/"Ar4+/H                             ","Ar$^{4+}$/H                        ", abundances_format, "abund_cel_arv                      "/)
        resultprocessingarray(104,:) = all_results%ar_icf_CEL
        resultprocessingtext(104,:) = (/"icf(Ar)                            ","icf(Ar)                            ", abundances_format, "icf_cel_ar                         "/)
        resultprocessingarray(105,:) = all_results%Ar_abund_CEL
        resultprocessingtext(105,:) = (/"Ar/H                               ","Ar$^{}$/H                          ", abundances_format, "abund_cel_ar                       "/)
        resultprocessingarray(106,:) = all_results%SII_abund_CEL
        resultprocessingtext(106,:) = (/"S+/H                               ","S$^{+}$/H                          ", abundances_format, "abund_cel_sii                      "/)
        resultprocessingarray(107,:) = all_results%SIII_abund_CEL
        resultprocessingtext(107,:) = (/"S2+/H                              ","S$^{2+}$/H                         ", abundances_format, "abund_cel_siii                     "/)
        resultprocessingarray(108,:) = all_results%s_icf_CEL
        resultprocessingtext(108,:) = (/"icf(S)                             ","icf(S)                             ", abundances_format, "icf_cel_s                          "/)
        resultprocessingarray(109,:) = all_results%S_abund_CEL
        resultprocessingtext(109,:) = (/"S/H                                ","S$^{}$/H                           ", abundances_format, "abund_cel_s                        "/)

        resultprocessingarray(110,:) = all_results%ClII_abund_CEL
        resultprocessingtext(110,:) = (/"Cl+/H                              ","Cl$^{+}$/H                         ", abundances_format, "abund_cel_clii                     "/)
        resultprocessingarray(111,:) = all_results%ClIII_abund_CEL
        resultprocessingtext(111,:) = (/"Cl2+/H                             ","Cl$^{2+}$/H                        ", abundances_format, "abund_cel_cliii                    "/)
        resultprocessingarray(112,:) = all_results%ClIV_abund_CEL
        resultprocessingtext(112,:) = (/"Cl3+/H                             ","Cl$^{3+}$/H                        ", abundances_format, "abund_cel_cliv                     "/)

        resultprocessingarray(113,:) = all_results%cl_icf_CEL
        resultprocessingtext(113,:) = (/"icf(Cl)                            ","icf(Cl)                            ", abundances_format, "icf_cel_cl                         "/)
        resultprocessingarray(114,:) = all_results%Cl_abund_CEL
        resultprocessingtext(114,:) = (/"Cl/H                               ","Cl$^{}$/H                          ", abundances_format, "abund_cel_cl                       "/)

!ORL abundances

        resultprocessingarray(115,:) = all_results%Hei_abund_ORL
        resultprocessingtext(115,:) = (/"He+/H                              ","He$^{+}$/H                         ", abundances_format, "abund_orl_hei                      "/)
        resultprocessingarray(116,:) = all_results%Heii_abund_ORL
        resultprocessingtext(116,:) = (/"He2+/H                             ","He$^{2+}$/H                        ", abundances_format, "abund_orl_heii                     "/)
        resultprocessingarray(117,:) = all_results%He_abund_ORL
        resultprocessingtext(117,:) = (/"He/H                               ","He/H                               ", abundances_format, "abund_orl_he                       "/)
        resultprocessingarray(118,:) = all_results%Cii_abund_ORL
        resultprocessingtext(118,:) = (/"C2+/H                              ","C$^{2+}$/H                         ", abundances_format, "abund_orl_cii                      "/)
        resultprocessingarray(119,:) = all_results%Ciii_abund_ORL
        resultprocessingtext(119,:) = (/"C3+/H                              ","C$^{3+}$/H                         ", abundances_format, "abund_orl_ciii                     "/)
        resultprocessingarray(120,:) = all_results%c_icf_ORL
        resultprocessingtext(120,:) = (/"icf(C)                             ","icf(C)                             ", abundances_format, "icf_orl_c                          "/)
        resultprocessingarray(121,:) = all_results%C_abund_ORL
        resultprocessingtext(121,:) = (/"C/H                                ","C/H                                ", abundances_format, "abund_orl_c                        "/)
        resultprocessingarray(122,:) = all_results%Nii_v3_abund_ORL
        resultprocessingtext(122,:) = (/"N2+/H (V3)                         ","N$^{2+}$/H (V3)                    ", abundances_format, "abund_orl_nii_v3                   "/)
        resultprocessingarray(123,:) = all_results%Nii_v5_abund_ORL
        resultprocessingtext(123,:) = (/"N2+/H (V5)                         ","N$^{2+}$/H (V5)                    ", abundances_format, "abund_orl_nii_v5                   "/)
        resultprocessingarray(124,:) = all_results%Nii_v8_abund_ORL
        resultprocessingtext(124,:) = (/"N2+/H (V8)                         ","N$^{2+}$/H (V8)                    ", abundances_format, "abund_orl_nii_v8                   "/)
        resultprocessingarray(125,:) = all_results%Nii_v12_abund_ORL
        resultprocessingtext(125,:) = (/"N2+/H (V12)                        ","N$^{2+}$/H (V12)                   ", abundances_format, "abund_orl_nii_v12                  "/)
        resultprocessingarray(126,:) = all_results%Nii_v20_abund_ORL
        resultprocessingtext(126,:) = (/"N2+/H (V20)                        ","N$^{2+}$/H (V20)                   ", abundances_format, "abund_orl_nii_v20                  "/)
        resultprocessingarray(127,:) = all_results%Nii_v28_abund_ORL
        resultprocessingtext(127,:) = (/"N2+/H (V28)                        ","N$^{2+}$/H (V28)                   ", abundances_format, "abund_orl_nii_v28                  "/)
        resultprocessingarray(128,:) = all_results%Nii_3d4f_abund_ORL
        resultprocessingtext(128,:) = (/"N2+/H (3d-4f)                      ","N$^{2+}$/H (3d-4f)                 ", abundances_format, "abund_orl_nii_3d4f                 "/)
        resultprocessingarray(129,:) = all_results%Nii_abund_ORL
        resultprocessingtext(129,:) = (/"N2+/H                              ","N$^{2+}$/H                         ", abundances_format, "abund_orl_nii                      "/)
        resultprocessingarray(130,:) = all_results%Niii_abund_ORL
        resultprocessingtext(130,:) = (/"N3+/H                              ","N$^{3+}$/H                         ", abundances_format, "abund_orl_niii                     "/)
        resultprocessingarray(131,:) = all_results%n_icf_ORL
        resultprocessingtext(131,:) = (/"icf(N)                             ","icf(N)                             ", abundances_format, "icf_orl_n                          "/)
        resultprocessingarray(132,:) = all_results%N_abund_ORL
        resultprocessingtext(132,:) = (/"N/H                                ","N/H                                ", abundances_format, "abund_orl_n                        "/)
        resultprocessingarray(133,:) = all_results%Oii_v1_abund_ORL
        resultprocessingtext(133,:) = (/"O2+/H (V1)                         ","O$^{2+}$/H (V1)                    ", abundances_format, "abund_orl_oii_v1                   "/)
        resultprocessingarray(134,:) = all_results%Oii_v2_abund_ORL
        resultprocessingtext(134,:) = (/"O2+/H (V2)                         ","O$^{2+}$/H (V2)                    ", abundances_format, "abund_orl_oii_v2                   "/)
        resultprocessingarray(135,:) = all_results%Oii_v5_abund_ORL
        resultprocessingtext(135,:) = (/"O2+/H (V5)                         ","O$^{2+}$/H (V5)                    ", abundances_format, "abund_orl_oii_v5                   "/)
        resultprocessingarray(136,:) = all_results%Oii_v10_abund_ORL
        resultprocessingtext(136,:) = (/"O2+/H (V10)                        ","O$^{2+}$/H (V10)                   ", abundances_format, "abund_orl_oii_v10                  "/)
        resultprocessingarray(137,:) = all_results%Oii_v11_abund_ORL
        resultprocessingtext(137,:) = (/"O2+/H (V11)                        ","O$^{2+}$/H (V11)                   ", abundances_format, "abund_orl_oii_v11                  "/)
        resultprocessingarray(138,:) = all_results%Oii_v12_abund_ORL
        resultprocessingtext(138,:) = (/"O2+/H (V12)                        ","O$^{2+}$/H (V12)                   ", abundances_format, "abund_orl_oii_v12                  "/)
        resultprocessingarray(139,:) = all_results%Oii_v19_abund_ORL
        resultprocessingtext(139,:) = (/"O2+/H (V19)                        ","O$^{2+}$/H (V19)                   ", abundances_format, "abund_orl_oii_v19                  "/)
        resultprocessingarray(140,:) = all_results%Oii_v20_abund_ORL
        resultprocessingtext(140,:) = (/"O2+/H (V20)                        ","O$^{2+}$/H (V20)                   ", abundances_format, "abund_orl_oii_v20                  "/)
        resultprocessingarray(141,:) = all_results%Oii_v25_abund_ORL
        resultprocessingtext(141,:) = (/"O2+/H (V25)                        ","O$^{2+}$/H (V25)                   ", abundances_format, "abund_orl_oii_v25                  "/)
        resultprocessingarray(142,:) = all_results%Oii_v28_abund_ORL
        resultprocessingtext(142,:) = (/"O2+/H (V28)                        ","O$^{2+}$/H (V28)                   ", abundances_format, "abund_orl_oii_v28                  "/)
        resultprocessingarray(143,:) = all_results%Oii_v33_abund_ORL
        resultprocessingtext(143,:) = (/"O2+/H (V33)                        ","O$^{2+}$/H (V33)                   ", abundances_format, "abund_orl_oii_v33                  "/)
        resultprocessingarray(144,:) = all_results%Oii_3d4f_abund_ORL
        resultprocessingtext(144,:) = (/"O2+/H (3d-4f)                      ","O$^{2+}$/H (3d-4f)                 ", abundances_format, "abund_orl_oii_3d4f                 "/)
        resultprocessingarray(145,:) = all_results%Oii_abund_ORL
        resultprocessingtext(145,:) = (/"O2+/H                              ","O$^{2+}$/H                         ", abundances_format, "abund_orl_oii                      "/)
        resultprocessingarray(146,:) = all_results%O_icf_ORL
        resultprocessingtext(146,:) = (/"icf(O)                             ","icf(O)                             ", abundances_format, "icf_orl_o                          "/)
        resultprocessingarray(147,:) = all_results%O_abund_ORL
        resultprocessingtext(147,:) = (/"O/H                                ","O/H                                ", abundances_format, "abund_orl_o                        "/)
        resultprocessingarray(148,:) = all_results%Neii_abund_ORL
        resultprocessingtext(148,:) = (/"Ne2+/H                             ","Ne$^{2+}$/H                        ", abundances_format, "abund_orl_neii                     "/)
        resultprocessingarray(149,:) = all_results%ne_icf_ORL
        resultprocessingtext(149,:) = (/"icf(Ne)                            ","icf(Ne)                            ", abundances_format, "icf_orl_ne                         "/)
        resultprocessingarray(150,:) = all_results%Ne_abund_ORL
        resultprocessingtext(150,:) = (/"Ne/H                               ","Ne/H                               ", abundances_format, "abund_orl_ne                       "/)

!strong line abundances

        resultprocessingarray(151,:) = all_results%O_R23_upper
        resultprocessingtext(151,:) = (/"O/H (R23 upper)                    ","O/H (R23 upper)                    ", abundances_format, "abund_o_r23_upper                  "/)
        resultprocessingarray(152,:) = all_results%O_R23_lower
        resultprocessingtext(152,:) = (/"O/H (R23 lower)                    ","O/H (R23 lower)                    ", abundances_format, "abund_o_r23_lower                  "/)
        resultprocessingarray(153,:) = all_results%O_N2
        resultprocessingtext(153,:) = (/"O/H (N2)                           ","O/H (N2)                           ", abundances_format, "abund_o_n2                         "/)
        resultprocessingarray(154,:) = all_results%O_O3N2
        resultprocessingtext(154,:) = (/"O/H (O3N2)                         ","O/H (O3N2)                         ", abundances_format, "abund_o_o3n2                       "/)
        resultprocessingarray(155,:) = all_results%O_Ar3O3
        resultprocessingtext(155,:) = (/"O/H (Ar3O3)                        ","O/H (Ar3O3)                        ", abundances_format, "abund_o_ar3o3                      "/)
        resultprocessingarray(156,:) = all_results%O_S3O3
        resultprocessingtext(156,:) = (/"O/H (S3O3)                         ","O/H (S3O3)                         ", abundances_format, "abund_o_s3o3                       "/)

!adfs

        resultprocessingarray(157,:) = all_results%adf_o2plus
        resultprocessingtext(157,:) = (/"adf (O2+/H)                        ","adf (O$^{2+}$/H)                   ", adf_format, "adf_o2plus                         "/)
        resultprocessingarray(158,:) = all_results%adf_o
        resultprocessingtext(158,:) = (/"adf (O/H)                          ","adf (O/H)                          ", adf_format, "adf_o                              "/)
        resultprocessingarray(159,:) = all_results%adf_n2plus
        resultprocessingtext(159,:) = (/"adf (N2+/H)                        ","adf (N$^{2+}$/H)                   ", adf_format, "adf_n2plus                         "/)
        resultprocessingarray(160,:) = all_results%adf_n
        resultprocessingtext(160,:) = (/"adf (N/H)                          ","adf (N/H)                          ", adf_format, "adf_n                              "/)
        resultprocessingarray(161,:) = all_results%adf_c2plus
        resultprocessingtext(161,:) = (/"adf (C2+/H)                        ","adf (C$^{2+}$/H)                   ", adf_format, "adf_c2plus                         "/)
        resultprocessingarray(162,:) = all_results%adf_c
        resultprocessingtext(162,:) = (/"adf (C/H)                          ","adf (C/H)                          ", adf_format, "adf_c                              "/)
        resultprocessingarray(163,:) = all_results%adf_ne2plus
        resultprocessingtext(163,:) = (/"adf (Ne2+/H)                       ","adf (Ne$^{2+}$/H)                  ", adf_format, "adf_ne2plus                        "/)
        resultprocessingarray(164,:) = all_results%adf_ne
        resultprocessingtext(164,:) = (/"adf (Ne/H)                         ","adf (Ne/H)                         ", adf_format, "adf_ne                             "/)

!open the files and write the headers

        open (650,file=trim(filename)//"_results", status='replace', access='sequential', action='write')
        open (651,file=trim(filename)//"_results.tex", status='replace', access='sequential', action='write')

        write (650,*) "NEAT (nebular empirical analysis tool)"
        write (650,*) "======================================"
        write (650,*) "Version ",VERSION
        write (650,*)
        write (650,*) "Analysis of file ",trim(filename)
        write (650,*) "Command line: ",trim(commandline)
        write (650,*)

        write (651,*) "\noindent{\Large {\sc neat} (nebular empirical analysis tool)}"
        write (651,*) "\noindent Version ",VERSION
        write (651,*) "\hrule"
        write (651,*) "\vspace{0.3cm}"
        write (651,*) "\noindent Analysis of file {\tt ",trim(filename),"}\newline"
        write (651,*) "\noindent Command line: {\tt ",trim(commandline),"}\newline"
        write (651,*) "\begin{longtable}[l]{ll}"

!next, loop through the results, processing and printing

        do j=1,164

! here we put some if statements to put things into conveniently separate bits

          if (j .eq. 1) then
            write (650,"(/A,/A/)") "Extinction","=========="
            write (651,*) "\multicolumn{2}{l}{Extinction}\\ \hline"
          elseif (j .eq. 5) then
            write (650,"(/A,/A)") "Diagnostics","==========="
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Diagnostics}\\ \hline"
            write (650,"(/A,/A/)") "Low ionisation densities","-----------"
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Low ionisation densities}\\ \hline"
          elseif (j .eq. 10) then
            write (650,"(/A,/A/)") "Low ionisation temperatures","-----------"
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Low ionisation temperatures}\\ \hline"
          elseif (j .eq. 21) then
            write (650,"(/A,/A/)") "Medium ionisation densities","-----------"
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Medium ionisation densities}\\ \hline"
          elseif (j .eq. 36) then
            write (650,"(/A,/A/)") "Medium ionisation temperatures","-----------"
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Medium ionisation temperatures}\\ \hline"
          elseif (j .eq. 51) then
            write (650,"(/A,/A/)") "High ionisation densities","-----------"
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{High ionisation densities}\\ \hline"
          elseif (j .eq. 54) then
            write (650,"(/A,/A/)") "High ionisation temperatures","-----------"
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{High ionisation temperatures}\\ \hline"
          elseif (j .eq. 59) then
            write (650,"(/A,/A/)") "Recombination line diagnostics","-----------"
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Recombination line diagnostics}\\ \hline"
          elseif (j .eq. 71) then
            if (subtract_recombination .eq. 2) then
              write (650,"(/A,/A/)") "Recombination contribution to CELs (%) (RL value has been subtracted from line fluxes)","-----------"
              write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Recombination contribution to CELs (\%)}\\ \hline"
            else
              write (650,"(/A,/A/)") "Recombination contribution to CELs (%) (for information only, not subtracted from line fluxes)","-----------"
              write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Recombination contribution to CELs (\%)}\\ \hline"
            endif
          elseif (j .eq. 77) then
            write (650,"(/A,/A/)") "CEL abundances","=============="
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{CEL abundances}\\ \hline"
          elseif (j .eq. 115) then
            write (650,"(/A,/A/)") "ORL abundances","=============="
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{ORL abundances}\\ \hline"
          elseif (j .eq. 118 .or. j .eq. 122 .or. j .eq. 133 .or. j .eq. 148) then
            write (650,"(/A,/A/)")
            write (651,*) "\\"
          elseif (j .eq. 151) then
            write (650,"(/A,/A/)") "Strong line abundances","======================"
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Strong line abundances}\\ \hline"
          elseif (j .eq. 157) then
            write (650,"(/A,/A/)") "Abundance discrepancy factors","============================="
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Abundance discrepancy factors}\\ \hline"
          endif

! this writes the results to the plain text and latex summary files

          quantity_result=resultprocessingarray(j,:)
          call write_uncertainties(quantity_result,uncertainty_array,resultprocessingtext(j,1),resultprocessingtext(j,2),resultprocessingtext(j,3),filename, resultprocessingtext(j,4), verbosity,nbins)

        enddo

!write ends of files, close

        write (651,*) "\end{longtable}"

        close(650)
        close(651)

        print *,gettime(),"Results written to files:"
        print *,"              ",trim(filename),"_results"
        print *,"              ",trim(filename),"_results.tex"

        print *
        print *,gettime(),"all done"
        print *

contains

subroutine randomizer(linelist, listlength, norp)
        ! from http://www.netlib.org/random/random.f90

        type(line), dimension(listlength) :: linelist
        integer :: IO, I, j, listlength
        real(kind=dp) :: temp4
        logical, intent(in) :: norp

        real(kind=dp) :: fn_val
        real(kind=dp)     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           &
                    r1 = 0.27597, r2 = 0.27846, u, v, x, y, q
        real(kind=dp) :: half
        real(kind=dp) :: newmean, newsnr, snr

!debugging
#ifdef CO
        print *,"subroutine: randomizer"
#endif

        half = 0.5

        I = 1
        IO=0

        do j = 1,listlength
          do
            call random_number(u)
            call random_number(v)
            v = 1.7156 * (v - half)
            x = u - s
            y = abs(v) - t
            q = x**2 + y*(a*y - b*x)
            if (q .lt. r1) exit
            if (q .gt. r2) cycle
            if (v**2 .lt. -4.0*log(u)*u**2) exit
          enddo
          fn_val = v/u

!                if (j==1) R=3.1+(0.15*fn_val)

          if (linelist(j)%int_err .eq. 0.d0 .and. .not. norp) then ! can't calculate signal to noise for RP calculations
            linelist(j)%intensity = 0.d0
          elseif (linelist(j)%intensity/linelist(j)%int_err .gt. 6.0 .or. norp) then !normal distribution
            temp4=linelist(j)%intensity+(fn_val*linelist(j)%int_err)
            if(temp4 .lt. 0) temp4 = 0.D0
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
            if (temp4 .lt. 0 ) temp4 = 0.D0
            linelist(j)%intensity = linelist(j)%intensity / temp4

          endif
        enddo

end subroutine randomizer

subroutine init_random_seed()
          integer :: i, n, clock
          integer, dimension(:), allocatable :: seed

!debugging
#ifdef CO
        print *,"subroutine: init_random_seed"
#endif

          n=20
          i=n
          call random_seed(size = n)
          allocate(seed(n))

          call system_clock(count=clock)

          seed = clock + 37 * (/ (i - 1, i = 1, n) /)
          call random_seed(put = seed)

          deallocate(seed)

end subroutine init_random_seed

subroutine write_uncertainties(input_array, uncertainty_array, plaintext, latextext, itemformat, filename, suffix, verbosity,nbins)

!wrapper for the get_uncertainties routine, if called it will do the
!uncertainty calculation and also write the binned and unbinned results to files
!as necessary.
implicit none
real(kind=dp) :: input_array(:)
real(kind=dp), intent(out) :: uncertainty_array(3)
type(arraycount), dimension (:), allocatable :: binned_quantity_result
character(len=35), intent(in) :: plaintext, latextext
character(len=35), intent(in) :: itemformat
character(len=512), intent(in) :: filename
character(len=25), intent(in) :: suffix
logical :: unusual
integer :: verbosity !1 = write summaries, binned and full results; 2=write summaries and binned results; 3=write summaries only
integer :: nbins

!debugging
#ifdef CO
        print *,"subroutine: write_uncertainties, quantity=",suffix
#endif

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

call get_uncertainties(input_array, binned_quantity_result, uncertainty_array, unusual,nbins)

!write out the binned results with the mode+-uncertainties at the top

if (verbosity .lt. 3) then

  if (allocated(binned_quantity_result) .and. maxval(uncertainty_array) .ne. minval(uncertainty_array)) then
    open(850, file=trim(filename)//"_"//trim(suffix)//"_binned", status='replace',access='sequential', action='write')
    write(unit = 850,fmt=*) uncertainty_array(2),uncertainty_array(2)-uncertainty_array(1),uncertainty_array(3)+uncertainty_array(2)
    write(unit = 850,fmt=*)
  !  do i=1,ii-1
    do i=1,size(binned_quantity_result, 1)
  ! this hacky condition is because for reasons I can't work out right now, the
  ! number of bins allocated to the array is always too large and the last few end
  ! up being full of zeros.  to be fixed soon hopefully.  RW 16/11/2012
!fixed 12/06/2017 - I guess 4.5 years wasn't quite within the hoped for "soon"
!      if (binned_quantity_result(i)%value .gt. 0. .and. binned_quantity_result(i)%counts .gt. 0) then
        write(unit = 850,fmt=*) binned_quantity_result(i)%value,binned_quantity_result(i)%counts
!      endif
    enddo
    close(850)
  endif

endif

!write the unbinned results to file

if (verbosity .lt. 2) then

  open(850, file=trim(filename)//"_"//trim(suffix), status='replace', access='sequential', action='write')
  do i=1,size(input_array)
    write(unit = 850,fmt=*) input_array(i)
  enddo
  close(850)

endif

!write derived value and uncertainties to the summary files

  if (maxval(uncertainty_array) .gt. 0.) then !if this condition is not true, array will be full of zeroes
!    if (unusual) write (650,itemformat) "Warning! Unusual probability distribution.  You should inspect this one:"
    write (650,itemformat) plaintext,uncertainty_array(2),uncertainty_array(3),-uncertainty_array(1)
!    write (651,*) latextext," & $",trim(latex_number(uncertainty_array(2))),"^{",trim(latex_number(uncertainty_array(3))),"}_{",trim(latex_number(-uncertainty_array(1))),"}$ \\"
    if (uncertainty_array(1) .eq. uncertainty_array(3)) then
      write (651,"(A,' & ${',A,'}\pm{',A,'}$ \\')") latextext,trim(latex_number(uncertainty_array(2))),trim(latex_number(uncertainty_array(1)))
    else
      write (651,"(A,' & ${',A,'}^{+',A,'}_{',A,'}$ \\')") latextext,trim(latex_number(uncertainty_array(2))),trim(latex_number(uncertainty_array(3))),trim(latex_number(-uncertainty_array(1)))
    endif
  else
    write (650,*) plaintext,"--  --  --"
    write (651,*) latextext,"& -- & -- & -- \\"
  endif

endif !end of condition checking whether one or more iterations were done

end subroutine write_uncertainties

subroutine get_uncertainties(input_array, binned_quantity_result, uncertainty_array, unusual,nbins)

implicit none
real(kind=dp) :: input_array(:)
real(kind=dp),dimension(:),allocatable :: logarray,exparray
real(kind=dp), intent(out) :: uncertainty_array(3)
real(kind=dp), dimension(:), allocatable :: bintemp
real(kind=dp) :: binsize=0d0
type(arraycount), dimension (:), allocatable, intent(out) :: binned_quantity_result
integer :: arraysize
integer :: nbins
integer :: bincount, ii
real(kind=dp) :: comp
real(kind=dp) :: mean=0d0, sd=0d0
real(kind=dp) :: mean_log=0d0, sd_log=0d0
real(kind=dp) :: mean_exp=0d0, sd_exp=0d0
real(kind=dp), dimension(3,3) :: sds=0d0
real(kind=dp) :: tolerance=0d0
real(kind=dp) :: rounding ! value to round to nearest to
logical :: unusual !if not true, then the distribution is normal, log normal or exp normal

!debugging
#ifdef CO
        print *,"subroutine: get_uncertainties. "
#endif

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

! bin the array. use 1/25 x 1 sigma above median as bin size

arraysize = size(input_array)
binsize=input_array(nint(0.841*arraysize))/25d0
nbins=0

if (binsize .gt. 0d0) then

!quantize the input array by taking the integer of each value divided by the binsize, and multiplying by the binsize.

  allocate(bintemp(arraysize))
  bintemp = binsize*nint(input_array/binsize)

  ! count the unique values
  comp=bintemp(1)
  ii=0
  do i=1,arraysize
    if (comp.ne.bintemp(i)) then
      comp=bintemp(i)
      ii=ii+1
    endif
  enddo

  nbins=ii

endif

if (nbins.gt.0) then
  allocate(binned_quantity_result(nbins))
  binned_quantity_result%value=0.D0
  binned_quantity_result%counts=0

!then, go through the quantized array and count the number of occurrences of each value

  comp=bintemp(1)
  bincount=0
  ii=1

  do i=1,arraysize
    if (bintemp(i).ne.comp) then
      binned_quantity_result(ii)%value=comp
      binned_quantity_result(ii)%counts=count(bintemp.eq.comp)
      comp=bintemp(i)
      ii=ii+1
    endif
  enddo

!get median and 1 sigma limits

  uncertainty_array = (/ input_array(nint(0.5d0*arraysize))-input_array(nint(0.159*arraysize)), input_array(nint(0.5d0*arraysize)), input_array(nint(0.841*arraysize))-input_array(nint(0.5d0*arraysize)) /)

!simple test for normality - calculate the mean and standard deviation of the array, determine the number of values within 1, 2 and 3 sigma of the mean.
!if the fractions are close to 68.27, 95.45 and 99.73 then it's normal
!also calculate for log and exp. Temperatures often have an exp-normal distribution but 10**(a typical temperature) exceeds HUGE(0.d0).  Scale to avoid this.

  allocate(logarray(size(input_array)))
  allocate(exparray(size(input_array)))

  where (input_array.gt.0)
    logarray=log10(input_array)
  elsewhere
    logarray=0.d0
  endwhere

  exparray=input_array/maxval(input_array)
  exparray=10**exparray

  mean=sum(input_array)/arraysize
  mean_log=sum(logarray)/arraysize
  mean_exp=sum(exparray)/arraysize

  sd=(sum((input_array-mean)**2)/arraysize)**0.5
  sd_log=(sum((logarray-mean_log)**2,input_array.gt.0.d0)/arraysize)**0.5
  sd_exp=(sum((exparray-mean_exp)**2)/arraysize)**0.5

  sds(1,1)=count(abs(input_array-mean).lt.sd)
  sds(1,2)=count(abs(input_array-mean).lt.2*sd)
  sds(1,3)=count(abs(input_array-mean).lt.3*sd)

  sds(2,1)=count(abs(logarray-mean_log).lt.sd_log)
  sds(2,2)=count(abs(logarray-mean_log).lt.2*sd_log)
  sds(2,3)=count(abs(logarray-mean_log).lt.3*sd_log)

  sds(1,1)=count(abs(exparray-mean_exp).lt.sd_exp)
  sds(1,2)=count(abs(exparray-mean_exp).lt.2*sd_exp)
  sds(1,3)=count(abs(exparray-mean_exp).lt.3*sd_exp)

  sds=sds/arraysize

!scale exp values back up

  mean_exp=mean_exp*maxval(input_array)
  sd_exp=sd_exp*maxval(input_array)

  if (abs(sds(1,1)-0.6827).lt.tolerance .and. abs(sds(1,2)-0.9545).lt.(tolerance*0.5) .and. abs(sds(1,3)-0.9973).lt. (tolerance*0.1)) then
    uncertainty_array(:)=(/sd,mean,sd/)
  elseif (abs(sds(2,1)-0.6827).lt.tolerance .and. abs(sds(2,2)-0.9545).lt.(tolerance*0.5) .and. abs(sds(2,3)-0.9973).lt. (tolerance*0.1)) then
    uncertainty_array(:)=(/10**(mean_log)-10**(mean_log-sd_log),10**(mean_log),10**(mean_log+sd_log)-10**(mean_log)/)
  elseif (abs(sds(3,1)-0.6827).lt.tolerance .and. abs(sds(3,2)-0.9545).lt.(tolerance*0.5) .and. abs(sds(3,3)-0.9973).lt. (tolerance*0.1)) then
    uncertainty_array(:)=(/log10(mean_exp+sd_exp)-log10(mean_exp),log10(mean_exp),log10(mean_exp-sd_exp)-log10(mean_exp)/)
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

!round values to 3 significant figures
!if uncertainties get rounded to zero, increment them

if (uncertainty_array(2) .ne. 0.d0 .and. uncertainty_array(2) .ne. 1.d0) then
  if (uncertainty_array(1).ge.0.001 .and.uncertainty_array(1).lt.1.0) then
    rounding=0.001d0
  else
    rounding=10**dble(floor(log10(uncertainty_array(2))))/100.d0
  endif

  uncertainty_array = rounding * nint(uncertainty_array/rounding)
  uncertainty_array(1) = max(rounding,uncertainty_array(1))
  uncertainty_array(3) = max(rounding,uncertainty_array(3))

endif

end subroutine get_uncertainties

character(len=11) function gettime()
implicit none
character(len=10) :: time
character(len=11), save :: oldtime

!debugging
#ifdef CO
        !print *,"function: gettime"
#endif

  call date_and_time(TIME=time)
  gettime = time(1:2)//":"//time(3:4)//":"//time(5:6)//" : "
  if (gettime .eq. oldtime) then
    gettime = "          "
  else
    oldtime = gettime
  endif
  return

end function gettime

character(len=100) function latex_number(inputnumber)
! for any number in the form aEx, write it in latex format, ie a$\times$10$^{x}$
real(kind=dp) :: inputnumber, mantissa
integer :: exponent, pos

!debugging
#ifdef CO
        print *,"function: latex_number, ",inputnumber
#endif

write (latex_number,"(ES14.6)") inputnumber
pos = index (latex_number,'E')
read (latex_number(1:pos-1), '(F6.3)') mantissa
read (latex_number(pos+1:14), '(I3)') exponent

if (exponent .ge. -2 .and. exponent .le. 1) then
  write (latex_number,"(F6.2)") inputnumber !just print out normal number if it's between 0.01 & 10
elseif (exponent .ge. 2 .and. exponent .le. 4) then
  write (latex_number,"(I6)") 10*nint(inputnumber/10) ! write out integer rounded to nearest 10 if it's between 10 and 10,000
else !otherwise, write out a formatted exponent
  write (latex_number,"(F6.2,'\times 10^{',I3,'}')") mantissa,exponent
endif
return

end function latex_number
end program
