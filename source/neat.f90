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
        !use mod_common_data

        IMPLICIT NONE

        CHARACTER :: switch_ext !switch for extinction laws
        CHARACTER :: switch_he  !switch for helium atomic data
        CHARACTER :: switch_icf !switch for which ICF scheme to use
        INTEGER :: I, j, runs, Narg !runs = number of runs for randomiser

        !input options

        CHARACTER*2048, DIMENSION(:), allocatable :: options
        CHARACTER*2048 :: commandline

        !file reading variables

        LOGICAL :: file_exists
        TYPE(LINE),DIMENSION(:), allocatable :: linelist
        TYPE(LINE),DIMENSION(:), allocatable :: linelist_original
        TYPE(LINE),dimension(:,:), allocatable :: all_linelists
        CHARACTER*80 :: filename
        CHARACTER*1 :: null
        INTEGER :: IO, listlength

        !results and result processing

        type(resultarray), dimension(:), allocatable :: all_results
        type(resultarray), dimension(1) :: iteration_result
        double precision, dimension(:), allocatable :: quantity_result
        double precision, dimension(:,:), allocatable :: resultprocessingarray
        character*40, dimension(:,:), allocatable :: resultprocessingtext

        !atomic data

        character*10 :: ionlist(40) !list of ion names
        integer :: iion !# of ions in Ilines
        integer :: Iint 
        integer :: maxlevs,maxtemps
        type(atomic_data),allocatable :: atomicdata(:)
        double precision, dimension(:,:,:), allocatable :: heidata

        !extinction

        logical :: calculate_extinction=.true.
        DOUBLE PRECISION :: temp1,temp2,temp3, meanextinction, R

        !binning and uncertainties

        double precision, dimension(3) :: uncertainty_array
        double precision, dimension (:,:), allocatable :: binned_quantity_result
        logical :: unusual
        integer :: verbosity

        !CEL array

        TYPE(line), DIMENSION(:), allocatable :: ILs

        !diagnostic array

        double precision, dimension(6) :: diagnostic_array

        !output formats

        character*35 :: extinction_format, diagnostic_format, diagnostic_ratio_format, abundances_format, adf_format

        !multiple output formats defined as variables so they can be passed to
        !the printout subroutine

        extinction_format = "(X,A,F9.3,SP,F9.3,F9.3,S)"
        diagnostic_format = "(X,A,F9.1,SP,F9.1,F9.1,S)"
        diagnostic_ratio_format = "(X,A,F9.3,SP,F9.3,F9.3,S)"
        abundances_format = "(X,A,ES14.3,SP,ES14.3,ES14.3,S)"
        adf_format = "(X,A,F8.2,SP,F8.2,F8.2,S)"

        !single iteration formats defined as normal format statements

        700 FORMAT(X,A20,F5.3) !extinction format
        701 FORMAT(X,A20,I5) !diagnostics format
        705 FORMAT(X,A20,I5,5X,F8.3) !diagnostics and ratio format
        702 FORMAT(X,A20,ES14.3) !abundances format
        703 FORMAT(X,A20,F5.2) !strong line format
        704 FORMAT(X,A20,F5.2) !adf format
        710 FORMAT(X,A20,X,"&",X,F5.3,X,"\\") !latex table version of extinction format
        711 FORMAT(X,A20,X,"&",X,I5,X,"\\") !as above, diagnostics format
        715 FORMAT(X,A20,X,"&",X,I5,5X,"&",X,F8.3,X,"\\") !as above, diagnostics and ratio format
        712 FORMAT(X,A20,X,"&",X,ES14.3,X,"\\") !as above, abundances format
        713 FORMAT(X,A20,X,"&",X,F5.2,X,"\\") !as above, strong line format
        714 FORMAT(X,A20,X,"&",X,F5.2,X,"\\") !as above, adf format
        R=3.1

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
           print *,"       Default: KB94"
           print *,"  -v / --verbosity"
           print *,"       Amount of output to write for each derived quantity"
           print *,"       This option has no effect unless -n or -u is specified"
           print *,"       Values:"
           print *,"         1: write out summary files, binned results and complete results"
           print *,"         2: write out summary files and binned results"
           print *,"         3: write out summary files only"
           print *,"       Default: 1"
        !  to be fully implemented:
        !  -R                     : R (default 3.1)
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
                  elseif (trim(options(i+1)) == "Fit")then
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
                    runs=20000
                endif
                if ((trim(options(i))=="-icf" .or. trim(options(i))=="--ionisation-correction-scheme") .and. (i+1) .le. Narg) then
                  if (trim(options(i+1))=="PT92") then
                    switch_icf="P"
                  endif
                endif
                if ((trim(options(i))=="-v" .or.  trim(options(i))=="--verbosity") .and. (i+1) .le. Narg) then
                  read (options(i+1),*) verbosity
                  if (verbosity .lt. 1 .or. verbosity .gt. 3) then
                    print *,"Error: verbosity outside allowed range of 1-3"
                    stop
                  endif
                endif
         enddo

         if (Narg .eq. 1) then
           filename=trim(options(1))
         endif

         if (filename=="") then
                print *,"Error: No input file specified"
                stop
         endif

         inquire(file=filename, exist=file_exists) ! see if the input file is present

         if (.not. file_exists) then
                print *,"Error: input file ",trim(filename)," does not exist"
                stop
         endif

        !first, read in the line list

        print *,"NEAT, the Nebular Empirical Analysis Tool"
        print *,"-----------------------------------------"

        print *
        print *,gettime(),": starting code"
        print *,gettime(),": command line: ",trim(commandline)
        if (runs .gt. 1 .and. runs .lt. 5000) print *,gettime(),": number of iterations is low.  At least 5000 is recommended for good sampling of probability distributions"

        deallocate(options)

        I = 0
        OPEN(199, file=filename, iostat=IO, status='old')
                DO WHILE (IO >= 0)
                        READ(199,*,end=111) null
                        I = I + 1
                END DO
        111 print *
        listlength=I

!then allocate and read
        allocate (linelist(listlength))
        allocate (linelist_original(listlength))

        REWIND (199)
        DO I=1,listlength
                READ(199,*,end=110) temp1, temp2, temp3
                linelist(i)%wavelength = temp1
                linelist(i)%intensity = temp2
                linelist(i)%int_err = temp3
                linelist(i)%latextext = ""
        END DO
        CLOSE(199)

        110 PRINT "(X,A9,A11,I4,A15,I4,A9)", gettime(),": read in ", I - 1," lines (out of ",listlength," in file)"

        if (I - 1 .ne. listlength) then
                print *,gettime(),": line list reading failed"
                print *,"This can happen if it doesn't have three columns"
                stop
        endif

        if(linelist(1)%wavelength == 0)then
                PRINT*, gettime(),": cheese shop error - no inputs"
                STOP
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

                        call randomizer(linelist, listlength, R)
                        R=3.1 ! no randomisation
                        call abundances(linelist, switch_ext, listlength, iteration_result, R, meanextinction, calculate_extinction, ILs, Iint, diagnostic_array,iion,atomicdata,maxlevs,maxtemps, heidata, switch_he, switch_icf)

                        !store all line and derived quantity in arrays
                        all_linelists(:,i)=linelist 
                        all_results(i)=iteration_result(1)
                        !restore the unrandomized line list ready for the next iteration
                        linelist = linelist_original

                END DO

        else
                print*, gettime(), ": I didn't want to be a barber anyway. I wanted to be... a lumberjack!   Also, a positive number of runs helps.."
        endif

!now write all the lines to line list files, plain text and latex
!linelist first

        print *,gettime(),": Writing line list"

        allocate(quantity_result(runs))

        open (650,FILE=trim(filename)//"_linelist", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
        open (651,FILE=trim(filename)//"_linelist.tex", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')

        write (650,*) "Lambda  Ion           F(line)  I(line) X(line)/H"
        write (651,*) "\centering"
        write (651,*) "\small "
        write (651,*) "\begin{longtable}{lllll}"
        write (651,*) "\hline"
        write (651,*) " $ \lambda $ & Ion & $F \left( \lambda \right) $ & $I \left( \lambda \right) $ & $\frac{X(line)}{H}$ \\ \hline \hline "

        if (runs .gt. 1) then
                do j=1, listlength
!line flux - calculate the uncertainties analytically, direct from input if
!SNR>6, from the equations if SNR<6.

                write (650,"(X,F7.2,X,A11)", advance='no') all_linelists(j,1)%wavelength,all_linelists(j,1)%name
                write (651,"(X,F7.2,' & ',A15,' & ')", advance='no') all_linelists(j,1)%wavelength,all_linelists(j,1)%latextext

!dereddened flux
                quantity_result = all_linelists(j,:)%int_dered
                call get_uncertainties(quantity_result, binned_quantity_result, uncertainty_array, unusual)
                if (uncertainty_array(1) .ne. uncertainty_array(3)) then
                  write (650,"(F8.3,SP,F8.3,SP,F8.3)", advance='no') uncertainty_array(2),uncertainty_array(1),-uncertainty_array(3)
                  write (651,"(F8.3,'$^{',SP,F8.3,'}_{',SP,F8.3,'}$')", advance='no') uncertainty_array(2),uncertainty_array(1),-uncertainty_array(3)
                else
                  write (650,"(F8.3,A,F8.3,5X)", advance='no') uncertainty_array(2)," +-",uncertainty_array(1)
                  write (651,"(F8.3,'$\pm$',F8.3)", advance='no') uncertainty_array(2),uncertainty_array(1)
                endif

!abundance - write out if there is an abundance for the line, don't write
!anything except a line break if there is no abundance for the line.

                quantity_result = all_linelists(j,:)%abundance
                call get_uncertainties(quantity_result, binned_quantity_result, uncertainty_array, unusual)
                if (uncertainty_array(2) .ne. 0.D0) then
                  if (uncertainty_array(1) .ne. uncertainty_array(3)) then
                    write (650,"(ES10.2,SP,ES10.2,SP,ES10.2)") uncertainty_array(2),uncertainty_array(1),-uncertainty_array(3)
                    write (651,"(' & ${',A,'}^{+',A,'}_{',A,'}$ \\')") trim(latex_number(uncertainty_array(2))),trim(latex_number(uncertainty_array(1))),trim(latex_number(-uncertainty_array(3)))
                  else
                    write (650,"(ES10.2,A,ES10.2)") uncertainty_array(2)," +-",uncertainty_array(1)
                    write (651,"(' & $',A,'\pm',A,'$\\')") trim(latex_number(uncertainty_array(2))),trim(latex_number(uncertainty_array(1)))
                  endif
                else
                  write (650,*)
                  write (651,*) "\\"
                endif

                end do
        else ! runs == 1, no uncertainties to write out

                do i=1,listlength
                  if (linelist(i)%abundance .gt. 0.0) then
                     write (650,"(X,F7.2,X,A11,F8.3,X,F8.3,X,ES14.3)") linelist(i)%wavelength,linelist(i)%name,linelist(i)%intensity,linelist(i)%int_dered, linelist(i)%abundance
                     write (651,"(X,F7.2,X,'&',A15,'&',X,F8.3,X,'&',X,F8.3,X,'&',X,'$',A,'$',X,'\\')") linelist(i)%wavelength,linelist(i)%latextext,linelist(i)%intensity,linelist(i)%int_dered, trim(latex_number(linelist(i)%abundance))
                  else
                     write (650,"(X,F7.2,X,A11,F8.3,X,F8.3)") linelist(i)%wavelength,linelist(i)%name,linelist(i)%intensity,linelist(i)%int_dered
                     write (651,"(X,F7.2,X,'&',A15,'&',X,F8.3,X,'&',X,F8.3,X,'&',X,'\\')") linelist(i)%wavelength,linelist(i)%latextext,linelist(i)%intensity,linelist(i)%int_dered
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

        allocate(resultprocessingarray(140,runs))
        allocate(resultprocessingtext(140,4))

!extinction

        resultprocessingarray(1,:) = all_results%mean_cHb
        resultprocessingtext(1,:) = (/"c(Hb)                              ","c(H$\beta)$:                       ", extinction_format, "mean_chb                           "/)

!diagnostics
!low density

        resultprocessingarray(2,:) = all_results%oii_density
        resultprocessingtext(2,:) = (/"[OII] density                      ","{}[O~{\sc ii}] density:            ", diagnostic_format, "oii_density                        "/)
        resultprocessingarray(3,:) = all_results%oii_density_ratio
        resultprocessingtext(3,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "oii_density_ratio                  "/)
        resultprocessingarray(4,:) = all_results%SII_density
        resultprocessingtext(4,:) = (/"[SII] density                      ","{}[S~{\sc ii}] density:            ", diagnostic_format, "sii_density                        "/)
        resultprocessingarray(5,:) = all_results%sii_density_ratio
        resultprocessingtext(5,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "sii_density_ratio                  "/)
        resultprocessingarray(6,:) = all_results%low_density
        resultprocessingtext(6,:) = (/"Low ionisation density             ","Low ionisation density:            ", diagnostic_format, "low_density                        "/)

!low temperature

        resultprocessingarray(7,:) = all_results%oii_temp
        resultprocessingtext(7,:) = (/"[OII] temperature                  ","{}[O~{\sc ii}] temperature         ", diagnostic_format, "oii_temp                           "/)
        resultprocessingarray(8,:) = all_results%oii_temp_ratio
        resultprocessingtext(8,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "oii_temp_ratio                     "/)
        resultprocessingarray(9,:) = all_results%SII_temp
        resultprocessingtext(9,:) = (/"[SII] temperature                  ","{}[S~{\sc ii}] temperature         ", diagnostic_format, "SII_temp                           "/)
        resultprocessingarray(10,:) = all_results%sii_temp_ratio
        resultprocessingtext(10,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "sii_temp_ratio                     "/)
        resultprocessingarray(11,:) = all_results%NII_temp
        resultprocessingtext(11,:) = (/"[NII] temperature                  ","{}[N~{\sc ii}] temperature         ", diagnostic_format, "NII_temp                           "/)
        resultprocessingarray(12,:) = all_results%nii_temp_ratio
        resultprocessingtext(12,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "nii_temp_ratio                     "/)
        resultprocessingarray(13,:) = all_results%OI_temp
        resultprocessingtext(13,:) = (/"[OI] temperature                   ","{}[O~{\sc i}] temperature          ", diagnostic_format, "OI_temp                            "/)
        resultprocessingarray(14,:) = all_results%oi_temp_ratio
        resultprocessingtext(14,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "oi_temp_ratio                      "/)
        resultprocessingarray(15,:) = all_results%CI_temp
        resultprocessingtext(15,:) = (/"[CI] temperature                   ","{}[C~{\sc i}] temperature          ", diagnostic_format, "CI_temp                            "/)
        resultprocessingarray(16,:) = all_results%ci_temp_ratio
        resultprocessingtext(16,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "ci_temp_ratio                      "/)
        resultprocessingarray(17,:) = all_results%low_temp
        resultprocessingtext(17,:) = (/"Low ionisation temperature         ","Low ionisation temperature         ", diagnostic_format, "low_temp                           "/)


!medium density

        resultprocessingarray(18,:) = all_results%cliii_density
        resultprocessingtext(18,:) = (/"[ClIII] density                    ","{}[Cl~{\sc iii}] density:          ", diagnostic_format, "cliii_density                      "/)
        resultprocessingarray(19,:) = all_results%cliii_density_ratio
        resultprocessingtext(19,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "cliii_density_ratio                "/)
        resultprocessingarray(20,:) = all_results%ArIV_density
        resultprocessingtext(20,:) = (/"[ArIV] density                     ","{}[Ar~{\sc iv}] density:           ", diagnostic_format, "ariv_density                       "/)
        resultprocessingarray(21,:) = all_results%ariv_density_ratio
        resultprocessingtext(21,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "sii_density_ratio                  "/)
        resultprocessingarray(22,:) = all_results%CIII_density
        resultprocessingtext(22,:) = (/"[CIII] density                     ","{}[C~{\sc iii}] density            ", diagnostic_format, "CIII_density                       "/)
        resultprocessingarray(23,:) = all_results%ciii_density_ratio
        resultprocessingtext(23,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "ciii_density_ratio                 "/)
        resultprocessingarray(24,:) = all_results%OIII_IR_density
        resultprocessingtext(24,:) = (/"[OIII] IR density                  ","{}[O~{\sc iii}] IR density         ", diagnostic_format, "OIII_IR_density                    "/)
        resultprocessingarray(25,:) = all_results%oiii_ir_density_ratio
        resultprocessingtext(25,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "oiii_ir_density_ratio              "/)
        resultprocessingarray(26,:) = all_results%SIII_IR_density
        resultprocessingtext(26,:) = (/"[SIII] IR density                  ","{}[S~{\sc iii}] IR density         ", diagnostic_format, "SIII_IR_density                    "/)
        resultprocessingarray(27,:) = all_results%siii_ir_density_ratio
        resultprocessingtext(27,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "siii_ir_density_ratio              "/)
        resultprocessingarray(28,:) = all_results%ArIII_IR_density
        resultprocessingtext(28,:) = (/"[ArIII] IR density                 ","{}[Ar~{\sc iii}] IR density        ", diagnostic_format, "ArIII_IR_density                   "/)
        resultprocessingarray(29,:) = all_results%ariii_ir_density_ratio
        resultprocessingtext(29,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "ariii_ir_density_ratio             "/)
        resultprocessingarray(30,:) = all_results%NeIII_IR_density
        resultprocessingtext(30,:) = (/"[NeIII] IR density                 ","{}[Ne~{\sc iii}] IR density        ", diagnostic_format, "NeIII_IR_density                   "/)
        resultprocessingarray(31,:) = all_results%neiii_ir_density_ratio
        resultprocessingtext(31,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "neiii_ir_density_ratio             "/)
        resultprocessingarray(32,:) = all_results%med_density
        resultprocessingtext(32,:) = (/"Medium ionisation density          ","Medium ionisation density          ", diagnostic_format, "med_density                        "/)

! medium temperature

        resultprocessingarray(33,:) = all_results%OIII_temp
        resultprocessingtext(33,:) = (/"[OIII] temperature                 ","{}[O~{\sc iii}] temperature        ", diagnostic_format, "OIII_temp                          "/)
        resultprocessingarray(34,:) = all_results%oiii_temp_ratio
        resultprocessingtext(34,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "oiii_temp_ratio                    "/)
        resultprocessingarray(35,:) = all_results%OIII_IR_temp
        resultprocessingtext(35,:) = (/"[OIII] IR temperature              ","{}[O~{\sc iii}] IR temperature     ", diagnostic_format, "OIII_IR_temp                       "/)
        resultprocessingarray(36,:) = all_results%oiii_ir_temp_ratio
        resultprocessingtext(36,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "oiii_ir_temp_ratio                 "/)
        resultprocessingarray(37,:) = all_results%NeIII_temp
        resultprocessingtext(37,:) = (/"[NeIII] temperature                ","{}[Ne~{\sc iii}] temperature       ", diagnostic_format, "NeIII_temp                         "/)
        resultprocessingarray(38,:) = all_results%neiii_temp_ratio
        resultprocessingtext(38,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "neiii_temp_ratio                   "/)
        resultprocessingarray(39,:) = all_results%NeIII_IR_temp
        resultprocessingtext(39,:) = (/"[NeIII] IR temperature             ","{}[Ne~{\sc iii}] IR temperature    ", diagnostic_format, "NeIII_IR_temp                      "/)
        resultprocessingarray(40,:) = all_results%neiii_ir_temp_ratio
        resultprocessingtext(40,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "neiii_ir_temp_ratio                "/)
        resultprocessingarray(41,:) = all_results%ArIII_temp
        resultprocessingtext(41,:) = (/"[ArIII] temperature                ","{}[Ar~{\sc iii}] temperature       ", diagnostic_format, "ArIII_temp                         "/)
        resultprocessingarray(42,:) = all_results%ariii_temp_ratio
        resultprocessingtext(42,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "ariii_temp_ratio                   "/)
        resultprocessingarray(43,:) = all_results%SIII_temp
        resultprocessingtext(43,:) = (/"[SIII] temperature                 ","{}[S~{\sc iii}] temperature        ", diagnostic_format, "SIII_temp                          "/)
        resultprocessingarray(44,:) = all_results%siii_temp_ratio
        resultprocessingtext(44,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "siii_temp_ratio                    "/)
        resultprocessingarray(45,:) = all_results%med_temp
        resultprocessingtext(45,:) = (/"Medium ionisation temperature      ","Medium ionisation temperature      ", diagnostic_format, "med_temp                           "/)

!high density

        resultprocessingarray(46,:) = all_results%neiv_density
        resultprocessingtext(46,:) = (/"[NeIV] density                     ","{}[Ne~{\sc iv}] density            ", diagnostic_format, "neiv_density                       "/)
        resultprocessingarray(47,:) = all_results%neiv_density_ratio
        resultprocessingtext(47,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "neiv_density_ratio                 "/)
        resultprocessingarray(48,:) = all_results%high_density
        resultprocessingtext(48,:) = (/"High ionisation density            ","High ionisation density            ", diagnostic_format, "high_density                       "/)

!high temperature

        resultprocessingarray(49,:) = all_results%ArV_temp
        resultprocessingtext(49,:) = (/"[ArV] temperature                  ","{}[Ar~{\sc v}] temperature         ", diagnostic_format, "ArV_temp                           "/)
        resultprocessingarray(50,:) = all_results%arv_temp_ratio
        resultprocessingtext(50,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "arv_temp_ratio                     "/)
        resultprocessingarray(51,:) = all_results%NeV_temp
        resultprocessingtext(51,:) = (/"[NeV] temperature                  ","{}[Ne~{\sc v}] temperature         ", diagnostic_format, "NeV_temp                           "/)
        resultprocessingarray(52,:) = all_results%nev_temp_ratio
        resultprocessingtext(52,:) = (/"Ratio                              ","Ratio                              ", diagnostic_ratio_format, "nev_temp_ratio                     "/)
        resultprocessingarray(53,:) = all_results%high_temp
        resultprocessingtext(53,:) = (/"High temperature                   ","High temperature                   ", diagnostic_format, "high_temp                          "/)

!balmer jump temperature

        resultprocessingarray(54,:) = all_results%Bal_jump_temp
        resultprocessingtext(54,:) = (/"BJ temperature                     ","BJ temperature                     ", diagnostic_format, "Bal_jump_temp                      "/)

!CEL abundances

        resultprocessingarray(55,:) = all_results%NC_abund_CEL
        resultprocessingtext(55,:) = (/"C0/H                               ","C$^{0}$/H                          ", abundances_format, "NC_abund_CEL                       "/)
        resultprocessingarray(56,:) = all_results%cii_abund_CEL
        resultprocessingtext(56,:) = (/"C+/H                               ","C$^{+}$/H                          ", abundances_format, "cii_abund_CEL                      "/)
        resultprocessingarray(57,:) = all_results%ciii_abund_CEL
        resultprocessingtext(57,:) = (/"C2+/H                              ","C$^{2+}$/H                         ", abundances_format, "ciii_abund_CEL                     "/)
        resultprocessingarray(58,:) = all_results%civ_abund_CEL
        resultprocessingtext(58,:) = (/"C3+/H                              ","C$^{3+}$/H                         ", abundances_format, "civ_abund_CEL                      "/)
        resultprocessingarray(59,:) = all_results%c_icf_CEL
        resultprocessingtext(59,:) = (/"icf(C)                             ","icf(C)                             ", abundances_format, "c_icf_CEL                          "/)
        resultprocessingarray(60,:) = all_results%C_abund_CEL
        resultprocessingtext(60,:) = (/"C/H                                ","C$^{}$/H                           ", abundances_format, "C_abund_CEL                        "/)
        resultprocessingarray(61,:) = all_results%nii_abund_CEL
        resultprocessingtext(61,:) = (/"N+/H                               ","N$^{+}$/H                          ", abundances_format, "nii_abund_CEL                      "/)
        resultprocessingarray(62,:) = all_results%niii_abund_CEL
        resultprocessingtext(62,:) = (/"N2+/H                              ","N$^{2+}$/H                         ", abundances_format, "niii_abund_CEL                     "/)
        resultprocessingarray(63,:) = all_results%niv_abund_CEL
        resultprocessingtext(63,:) = (/"N3+/H                              ","N$^{3+}$/H                         ", abundances_format, "niv_abund_CEL                      "/)
        resultprocessingarray(64,:) = all_results%nv_abund_CEL
        resultprocessingtext(64,:) = (/"N4+/H                              ","N$^{4+}$/H                         ", abundances_format, "nv_abund_CEL                       "/)
        resultprocessingarray(65,:) = all_results%n_icf_CEL
        resultprocessingtext(65,:) = (/"icf(N)                             ","icf(N)                             ", abundances_format, "n_icf_CEL                          "/)
        resultprocessingarray(66,:) = all_results%N_abund_CEL
        resultprocessingtext(66,:) = (/"N/H                                ","N$^{}$/H                           ", abundances_format, "N_abund_CEL                        "/)
        resultprocessingarray(67,:) = all_results%NO_abund_CEL
        resultprocessingtext(67,:) = (/"O0/H                               ","O$^{0}$/H                          ", abundances_format, "NO_abund_CEL                       "/)
        resultprocessingarray(68,:) = all_results%Oii_abund_CEL
        resultprocessingtext(68,:) = (/"O+/H                               ","O$^{+}$/H                          ", abundances_format, "Oii_abund_CEL                      "/)
        resultprocessingarray(69,:) = all_results%Oiii_abund_CEL
        resultprocessingtext(69,:) = (/"O2+/H                              ","O$^{2+}$/H                         ", abundances_format, "Oiii_abund_CEL                     "/)
        resultprocessingarray(70,:) = all_results%Oiv_abund_CEL
        resultprocessingtext(70,:) = (/"O3+/H                              ","O$^{3+}$/H                         ", abundances_format, "Oiv_abund_CEL                      "/)
        resultprocessingarray(71,:) = all_results%o_icf_CEL
        resultprocessingtext(71,:) = (/"icf(O)                             ","icf(O)                             ", abundances_format, "o_icf_CEL                          "/)
        resultprocessingarray(72,:) = all_results%O_abund_CEL
        resultprocessingtext(72,:) = (/"O/H                                ","O$^{}$/H                           ", abundances_format, "O_abund_CEL                        "/)
        resultprocessingarray(73,:) = all_results%NeII_abund_CEL
        resultprocessingtext(73,:) = (/"Ne+/H                              ","Ne$^{+}$/H                         ", abundances_format, "NeII_abund_CEL                     "/)
        resultprocessingarray(74,:) = all_results%NeIII_abund_CEL
        resultprocessingtext(74,:) = (/"Ne2+/H                             ","Ne$^{2+}$/H                        ", abundances_format, "NeIII_abund_CEL                    "/)
        resultprocessingarray(75,:) = all_results%NeIV_abund_CEL
        resultprocessingtext(75,:) = (/"Ne3+/H                             ","Ne$^{3+}$/H                        ", abundances_format, "NeIV_abund_CEL                     "/)
        resultprocessingarray(76,:) = all_results%NeV_abund_CEL
        resultprocessingtext(76,:) = (/"Ne4+/H                             ","Ne$^{4+}$/H                        ", abundances_format, "NeV_abund_CEL                      "/)
        resultprocessingarray(77,:) = all_results%ne_icf_CEL
        resultprocessingtext(77,:) = (/"icf(Ne)                            ","icf(Ne)                            ", abundances_format, "ne_icf_CEL                         "/)
        resultprocessingarray(78,:) = all_results%Ne_abund_CEL
        resultprocessingtext(78,:) = (/"Ne/H                               ","Ne$^{}$/H                          ", abundances_format, "Ne_abund_CEL                       "/)
        resultprocessingarray(79,:) = all_results%ArIII_abund_CEL
        resultprocessingtext(79,:) = (/"Ar2+/H                             ","Ar$^{2+}$/H                        ", abundances_format, "ArIII_abund_CEL                    "/)
        resultprocessingarray(80,:) = all_results%ArIV_abund_CEL
        resultprocessingtext(80,:) = (/"Ar3+/H                             ","Ar$^{3+}$/H                        ", abundances_format, "ArIV_abund_CEL                     "/)
        resultprocessingarray(81,:) = all_results%ArV_abund_CEL
        resultprocessingtext(81,:) = (/"Ar4+/H                             ","Ar$^{4+}$/H                        ", abundances_format, "ArV_abund_CEL                      "/)
        resultprocessingarray(82,:) = all_results%ar_icf_CEL
        resultprocessingtext(82,:) = (/"icf(Ar)                            ","icf(Ar)                            ", abundances_format, "ar_icf_CEL                         "/)
        resultprocessingarray(83,:) = all_results%Ar_abund_CEL
        resultprocessingtext(83,:) = (/"Ar/H                               ","Ar$^{}$/H                          ", abundances_format, "Ar_abund_CEL                       "/)
        resultprocessingarray(84,:) = all_results%SII_abund_CEL
        resultprocessingtext(84,:) = (/"S+/H                               ","S$^{+}$/H                          ", abundances_format, "SII_abund_CEL                      "/)
        resultprocessingarray(85,:) = all_results%SIII_abund_CEL
        resultprocessingtext(85,:) = (/"S2+/H                              ","S$^{2+}$/H                         ", abundances_format, "SIII_abund_CEL                     "/)
        resultprocessingarray(86,:) = all_results%s_icf_CEL
        resultprocessingtext(86,:) = (/"icf(S)                             ","icf(S)                             ", abundances_format, "s_icf_CEL                          "/)
        resultprocessingarray(87,:) = all_results%S_abund_CEL
        resultprocessingtext(87,:) = (/"S/H                                ","S$^{}$/H                           ", abundances_format, "S_abund_CEL                        "/)
        resultprocessingarray(88,:) = all_results%ClIII_abund_CEL
        resultprocessingtext(88,:) = (/"Cl2+/H                             ","Cl$^{2+}$/H                        ", abundances_format, "ClIII_abund_CEL                    "/)
        resultprocessingarray(89,:) = all_results%cl_icf_CEL
        resultprocessingtext(89,:) = (/"icf(Cl)                            ","icf(Cl)                            ", abundances_format, "cl_icf_CEL                         "/)
        resultprocessingarray(90,:) = all_results%Cl_abund_CEL
        resultprocessingtext(90,:) = (/"Cl/H                               ","Cl$^{}$/H                          ", abundances_format, "Cl_abund_CEL                       "/)

!ORL abundances

        resultprocessingarray(91,:) = all_results%Hei_abund_ORL
        resultprocessingtext(91,:) = (/"He+/H                              ","He$^{+}$/H                         ", abundances_format, "Hei_abund_ORL                      "/)
        resultprocessingarray(92,:) = all_results%Heii_abund_ORL
        resultprocessingtext(92,:) = (/"He2+/H                             ","He$^{2+}$/H                        ", abundances_format, "Heii_abund_ORL                     "/)
        resultprocessingarray(93,:) = all_results%He_abund_ORL
        resultprocessingtext(93,:) = (/"He/H                               ","He/H                               ", abundances_format, "He_abund_ORL                       "/)
        resultprocessingarray(94,:) = all_results%Cii_abund_ORL
        resultprocessingtext(94,:) = (/"C2+/H                              ","C$^{2+}$/H                         ", abundances_format, "Cii_abund_ORL                      "/)
        resultprocessingarray(95,:) = all_results%Ciii_abund_ORL
        resultprocessingtext(95,:) = (/"C3+/H                              ","C$^{3+}$/H                         ", abundances_format, "Ciii_abund_ORL                     "/)
        resultprocessingarray(96,:) = all_results%c_icf_ORL
        resultprocessingtext(96,:) = (/"icf(C)                             ","icf(C)                             ", abundances_format, "c_icf_ORL                          "/)
        resultprocessingarray(97,:) = all_results%C_abund_ORL
        resultprocessingtext(97,:) = (/"C/H                                ","C/H                                ", abundances_format, "C_abund_ORL                        "/)
        resultprocessingarray(98,:) = all_results%Nii_v3_abund_ORL
        resultprocessingtext(98,:) = (/"N2+/H (V3)                         ","N$^{2+}$/H (V3)                    ", abundances_format, "Nii_v3_abund_ORL                   "/)
        resultprocessingarray(99,:) = all_results%Nii_v5_abund_ORL
        resultprocessingtext(99,:) = (/"N2+/H (V5)                         ","N$^{2+}$/H (V5)                    ", abundances_format, "Nii_v5_abund_ORL                   "/)
        resultprocessingarray(100,:) = all_results%Nii_v8_abund_ORL
        resultprocessingtext(100,:) = (/"N2+/H (V8)                         ","N$^{2+}$/H (V8)                    ", abundances_format, "Nii_v8_abund_ORL                   "/)
        resultprocessingarray(101,:) = all_results%Nii_v12_abund_ORL
        resultprocessingtext(101,:) = (/"N2+/H (V12)                        ","N$^{2+}$/H (V12)                   ", abundances_format, "Nii_v12_abund_ORL                  "/)
        resultprocessingarray(102,:) = all_results%Nii_v20_abund_ORL
        resultprocessingtext(102,:) = (/"N2+/H (V20)                        ","N$^{2+}$/H (V20)                   ", abundances_format, "Nii_v20_abund_ORL                  "/)
        resultprocessingarray(103,:) = all_results%Nii_v28_abund_ORL
        resultprocessingtext(103,:) = (/"N2+/H (V28)                        ","N$^{2+}$/H (V28)                   ", abundances_format, "Nii_v28_abund_ORL                  "/)
        resultprocessingarray(104,:) = all_results%Nii_3d4f_abund_ORL
        resultprocessingtext(104,:) = (/"N2+/H (3d-4f)                      ","N$^{2+}$/H (3d-4f)                 ", abundances_format, "Nii_3d4f_abund_ORL                 "/)
        resultprocessingarray(105,:) = all_results%Nii_abund_ORL
        resultprocessingtext(105,:) = (/"N2+/H                              ","N$^{2+}$/H                         ", abundances_format, "Nii_abund_ORL                      "/)
        resultprocessingarray(106,:) = all_results%Niii_abund_ORL
        resultprocessingtext(106,:) = (/"N3+/H                              ","N$^{3+}$/H                         ", abundances_format, "Niii_abund_ORL                     "/)
        resultprocessingarray(107,:) = all_results%n_icf_ORL
        resultprocessingtext(107,:) = (/"icf(N)                             ","icf(N)                             ", abundances_format, "n_icf_ORL                          "/)
        resultprocessingarray(108,:) = all_results%N_abund_ORL
        resultprocessingtext(108,:) = (/"N/H                                ","N/H                                ", abundances_format, "N_abund_ORL                        "/)
        resultprocessingarray(109,:) = all_results%Oii_v1_abund_ORL
        resultprocessingtext(109,:) = (/"O2+/H (V1)                         ","O$^{2+}$/H (V1)                    ", abundances_format, "Oii_v1_abund_ORL                   "/)
        resultprocessingarray(110,:) = all_results%Oii_v2_abund_ORL
        resultprocessingtext(110,:) = (/"O2+/H (V2)                         ","O$^{2+}$/H (V2)                    ", abundances_format, "Oii_v2_abund_ORL                   "/)
        resultprocessingarray(111,:) = all_results%Oii_v5_abund_ORL
        resultprocessingtext(111,:) = (/"O2+/H (V5)                         ","O$^{2+}$/H (V5)                    ", abundances_format, "Oii_v5_abund_ORL                   "/)
        resultprocessingarray(112,:) = all_results%Oii_v10_abund_ORL
        resultprocessingtext(112,:) = (/"O2+/H (V10)                        ","O$^{2+}$/H (V10)                   ", abundances_format, "Oii_v10_abund_ORL                  "/)
        resultprocessingarray(113,:) = all_results%Oii_v11_abund_ORL
        resultprocessingtext(113,:) = (/"O2+/H (V11)                        ","O$^{2+}$/H (V11)                   ", abundances_format, "Oii_v11_abund_ORL                  "/)
        resultprocessingarray(114,:) = all_results%Oii_v12_abund_ORL
        resultprocessingtext(114,:) = (/"O2+/H (V12)                        ","O$^{2+}$/H (V12)                   ", abundances_format, "Oii_v12_abund_ORL                  "/)
        resultprocessingarray(115,:) = all_results%Oii_v19_abund_ORL
        resultprocessingtext(115,:) = (/"O2+/H (V19)                        ","O$^{2+}$/H (V19)                   ", abundances_format, "Oii_v19_abund_ORL                  "/)
        resultprocessingarray(116,:) = all_results%Oii_v20_abund_ORL        
        resultprocessingtext(116,:) = (/"O2+/H (V20)                        ","O$^{2+}$/H (V20)                   ", abundances_format, "Oii_v20_abund_ORL                  "/)
        resultprocessingarray(117,:) = all_results%Oii_v25_abund_ORL
        resultprocessingtext(117,:) = (/"O2+/H (V25)                        ","O$^{2+}$/H (V25)                   ", abundances_format, "Oii_v25_abund_ORL                  "/)
        resultprocessingarray(118,:) = all_results%Oii_v28_abund_ORL
        resultprocessingtext(118,:) = (/"O2+/H (V28)                        ","O$^{2+}$/H (V28)                   ", abundances_format, "Oii_v28_abund_ORL                  "/)
        resultprocessingarray(119,:) = all_results%Oii_v33_abund_ORL
        resultprocessingtext(119,:) = (/"O2+/H (V33)                        ","O$^{2+}$/H (V33)                   ", abundances_format, "Oii_v33_abund_ORL                  "/)
        resultprocessingarray(120,:) = all_results%Oii_3d4f_abund_ORL
        resultprocessingtext(120,:) = (/"O2+/H (3d-4f)                      ","O$^{2+}$/H (3d-4f)                 ", abundances_format, "Oii_3d4f_abund_ORL                 "/)
        resultprocessingarray(121,:) = all_results%Oii_abund_ORL
        resultprocessingtext(121,:) = (/"O2+/H                              ","O$^{2+}$/H                         ", abundances_format, "Oii_abund_ORL                      "/)
        resultprocessingarray(122,:) = all_results%O_icf_ORL
        resultprocessingtext(122,:) = (/"icf(O)                             ","icf(O)                             ", abundances_format, "o_icf_ORL                          "/)
        resultprocessingarray(123,:) = all_results%O_abund_ORL
        resultprocessingtext(123,:) = (/"O/H                                ","O/H                                ", abundances_format, "O_abund_ORL                        "/)
        resultprocessingarray(124,:) = all_results%Neii_abund_ORL
        resultprocessingtext(124,:) = (/"Ne2+/H                             ","Ne$^{2+}$/H                        ", abundances_format, "Neii_abund_ORL                     "/)
        resultprocessingarray(125,:) = all_results%ne_icf_ORL
        resultprocessingtext(125,:) = (/"icf(Ne)                            ","icf(Ne)                            ", abundances_format, "ne_icf_ORL                         "/)
        resultprocessingarray(126,:) = all_results%Ne_abund_ORL
        resultprocessingtext(126,:) = (/"Ne/H                               ","Ne/H                               ", abundances_format, "Ne_abund_ORL                       "/)

!strong line abundances

        resultprocessingarray(127,:) = all_results%O_R23_upper
        resultprocessingtext(127,:) = (/"O/H (R23 upper)                    ","O/H (R23 upper)                    ", abundances_format, "O_R23_upper                        "/)
        resultprocessingarray(128,:) = all_results%O_R23_lower
        resultprocessingtext(128,:) = (/"O/H (R23 lower)                    ","O/H (R23 lower)                    ", abundances_format, "O_R23_lower                        "/)
        resultprocessingarray(129,:) = all_results%O_N2
        resultprocessingtext(129,:) = (/"O/H (N2)                           ","O/H (N2)                           ", abundances_format, "O_N2                               "/)
        resultprocessingarray(130,:) = all_results%O_O3N2
        resultprocessingtext(130,:) = (/"O/H (O3N2)                         ","O/H (O3N2)                         ", abundances_format, "O_O3N2                             "/)
        resultprocessingarray(131,:) = all_results%O_Ar3O3
        resultprocessingtext(131,:) = (/"O/H (Ar3O3)                        ","O/H (Ar3O3)                        ", abundances_format, "O_Ar3O3                            "/)
        resultprocessingarray(132,:) = all_results%O_S3O3
        resultprocessingtext(132,:) = (/"O/H (S3O3)                         ","O/H (S3O3)                         ", abundances_format, "O_S3O3                             "/)

!adfs

        resultprocessingarray(133,:) = all_results%adf_o2plus
        resultprocessingtext(133,:) = (/"adf (O2+/H)                        ","adf (O$^{2+}$/H)                   ", adf_format, "adf_o2plus                         "/)
        resultprocessingarray(134,:) = all_results%adf_o
        resultprocessingtext(134,:) = (/"adf (O/H)                          ","adf (O/H)                          ", adf_format, "adf_o                              "/)
        resultprocessingarray(135,:) = all_results%adf_n2plus
        resultprocessingtext(135,:) = (/"adf (N2+/H)                        ","adf (N$^{2+}$/H)                   ", adf_format, "adf_n2plus                         "/)
        resultprocessingarray(136,:) = all_results%adf_n
        resultprocessingtext(136,:) = (/"adf (N/H)                          ","adf (N/H)                          ", adf_format, "adf_n                              "/)
        resultprocessingarray(137,:) = all_results%adf_c2plus
        resultprocessingtext(137,:) = (/"adf (C2+/H)                        ","adf (C$^{2+}$/H)                   ", adf_format, "adf_c2plus                         "/)
        resultprocessingarray(138,:) = all_results%adf_c
        resultprocessingtext(138,:) = (/"adf (C/H)                          ","adf (C/H)                          ", adf_format, "adf_c                              "/)
        resultprocessingarray(139,:) = all_results%adf_ne2plus
        resultprocessingtext(139,:) = (/"adf (Ne2+/H)                       ","adf (Ne$^{2+}$/H)                  ", adf_format, "adf_ne2plus                        "/)
        resultprocessingarray(140,:) = all_results%adf_ne
        resultprocessingtext(140,:) = (/"adf (Ne/H)                         ","adf (Ne/H)                         ", adf_format, "adf_ne                             "/)

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

        do j=1,140

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
          elseif (j .eq. 46) then
            write (650,"(/A,/A/)") "High ionisation densities","-----------"
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{High ionisation densities}\\ \hline"
          elseif (j .eq. 49) then
            write (650,"(/A,/A/)") "High ionisation temperatures","-----------"
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{High ionisation temperatures}\\ \hline"
          elseif (j .eq. 54) then
            write (650,"(/A,/A/)") "Balmer jump temperature","-----------"
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Balmer jump temperature}\\ \hline"
          elseif (j .eq. 55) then
            write (650,"(/A,/A/)") "CEL abundances","=============="
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{CEL abundances}\\ \hline"
          elseif (j .eq. 91) then
            write (650,"(/A,/A/)") "ORL abundances","=============="
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{ORL abundances}\\ \hline" 
          elseif (j .eq. 94 .or. j .eq. 98 .or. j .eq. 109 .or. j .eq. 124) then 
            write (650,"(/A,/A/)") 
            write (651,*) "\\"
          elseif (j .eq. 127) then
            write (650,"(/A,/A/)") "Strong line abundances","======================"
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Strong line abundances}\\ \hline" 
          elseif (j .eq. 133) then
            write (650,"(/A,/A/)") "Abundance discrepancy factors","============================="
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Abundance discrepancy factors}\\ \hline"
          endif

! this writes the results to the plain text and latex summary files

          quantity_result=resultprocessingarray(j,:)
          call write_uncertainties(quantity_result,uncertainty_array,resultprocessingtext(j,1),resultprocessingtext(j,2),resultprocessingtext(j,3),filename, resultprocessingtext(j,4), verbosity)

        enddo

!write ends of files, close

        write (651,*) "\end{longtable}"

        close(650)
        close(651)

        print *,gettime(),": Results written to files:"
        print *,"              ",trim(filename),"_results"
        print *,"              ",trim(filename),"_results.tex"

        print *
        print *,gettime(),": Finished."

contains

        subroutine randomizer(linelist, listlength, R)

                TYPE(line), dimension(listlength) :: linelist
                INTEGER :: IO, I, j, listlength
                DOUBLE PRECISION :: temp4, R

                REAL :: fn_val

                !     Local variables
                REAL     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           &
                            r1 = 0.27597, r2 = 0.27846, u, v, x, y, q
                REAL :: half
                REAL :: newmean, newsnr, snr

                half = 0.5

                I = 1
                IO=0

                        do j = 1,listlength

                                ! from http://www.netlib.org/random/random.f90

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

!                                if (j==1) R=3.1+(0.15*fn_val)

                                if (linelist(j)%intensity/linelist(j)%int_err .gt. 6.0) then !normal distribution

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
                !end do
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

subroutine write_uncertainties(input_array, uncertainty_array, plaintext, latextext, itemformat, filename, suffix, verbosity)

!wrapper for the get_uncertainties routine, if called it will do the
!uncertainty calculation and also write the binned and unbinned results to files
!as necessary.
implicit none
double precision :: input_array(:)
double precision, intent(out) :: uncertainty_array(3)
double precision, dimension (:,:), allocatable :: binned_quantity_result
character*35, intent(in) :: plaintext, latextext
character*35, intent(in) :: itemformat
character*80, intent(in) :: filename
character*25, intent(in) :: suffix
logical :: unusual
integer :: verbosity !1 = write summaries, binned and full results; 2=write summaries and binned results; 3=write summaries only

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

call get_uncertainties(input_array, binned_quantity_result, uncertainty_array, unusual)

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
        write(unit = 850,FMT=*) binned_quantity_result(i,1),int(binned_quantity_result(i,2))
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
      write (651,"(A,' & ${',A,'}^{+',A,'}_{',A,'}$ \\')") latextext,trim(latex_number(uncertainty_array(2))),trim(latex_number(uncertainty_array(1))),trim(latex_number(-uncertainty_array(3)))
    endif
  else
    write (650,*) plaintext,"--"
    write (651,*) latextext,"& -- \\"
  endif

endif !end of condition checking whether one or more iterations were done

end subroutine write_uncertainties

subroutine get_uncertainties(input_array, binned_quantity_result, uncertainty_array, unusual)

implicit none
double precision :: input_array(:)
double precision, intent(out) :: uncertainty_array(3)
double precision, dimension(:), allocatable :: bintemp
double precision :: binsize, comp
double precision, dimension (:,:), allocatable, intent(out) :: binned_quantity_result
integer :: ii, bincount, arraysize, abovepos, belowpos, nbins
integer, dimension(1) :: maxpos
double precision :: mean, sd
double precision :: mean_log, sd_log
double precision :: mean_exp, sd_exp
double precision, dimension(3,3) :: sds
double precision :: tolerance
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
binsize=(input_array(nint(0.841*arraysize)) - input_array(nint(0.159*arraysize)))/5

if (binsize .gt. 0) then

!quantize the input array by taking the integer of each value divided by the binsize, and multiplying by the binsize.

  allocate(bintemp(arraysize))
  bintemp = binsize*nint(input_array/binsize)
  nbins=nint((maxval(bintemp)-minval(bintemp))/binsize)+1
  allocate(binned_quantity_result(nbins,2))
  binned_quantity_result=0.D0

!then, go through the quantized array and count the number of occurrences of each value

  comp=bintemp(1)
  bincount=0
  ii=1

  do i=1,arraysize
    if (bintemp(i).eq.comp) then
      bincount=bincount+1 
    else 
      binned_quantity_result(ii,:)=(/comp, dble(bincount)/) 
      comp=bintemp(i)
      bincount=0
      ii=ii+1 
    endif
  enddo

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
  uncertainty_array(1) = 0.D0
  uncertainty_array(2) = input_array(1)
  uncertainty_array(3) = 0.D0
endif

end subroutine get_uncertainties

character*10 function gettime()
implicit none
character*10 :: time

  call DATE_AND_TIME(TIME=time)
  gettime = time(1:2)//":"//time(3:4)//":"//time(5:6)
  return

end function gettime

character*100 function latex_number(inputnumber)
! for any number in the form aEx, write it in latex format, ie a$\times$10$^{x}$
double precision :: inputnumber, mantissa
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
