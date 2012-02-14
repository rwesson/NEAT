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

        CHARACTER :: switch_ext !switch for extinction laws
        INTEGER :: I, runs, Narg !runs = number of runs for randomiser

        !time variables

        character*8 :: date
        character*10 :: time

        !input options

        CHARACTER*2048, DIMENSION(:), allocatable :: options

        !file reading variables

        TYPE(LINE),DIMENSION(:), allocatable :: linelist
        TYPE(LINE),DIMENSION(:), allocatable :: linelist_original
        CHARACTER*80 :: filename
        CHARACTER*1 :: null
        INTEGER :: IO, listlength
        type(resultarray), dimension(:), allocatable :: all_results
        type(resultarray), dimension(1) :: iteration_result
        double precision, dimension(:), allocatable :: quantity_result 
        double precision :: binvalue

        !extinction

        logical :: calculate_extinction=.true.
        DOUBLE PRECISION :: temp1,temp2,temp3, meanextinction, R

        !binning

        double precision :: mode
        double precision, dimension(3) :: median_array

        !CEL array

        TYPE(line), DIMENSION(:), allocatable :: ILs

        R=3.1

        !read command line arguments

        Narg = IARGC() !count input arguments

        if (Narg .eq. 0) then
           print *,"Syntax: ./abundances.exe [option1 value1] [option2 value2] .. [optionx valuex]"
           print *,"Available options:" 
           print *,"  -i / --input"
           print *,"       Input file"
           print *,"       No default"
           print *,"  -n / --n-iterations"
           print *,"       Number of iterations"
           print *,"       Default: 1"
           print *,"  -e / --extinction-law"
           print *,"       Extinction law"
           print *,"       Default: Howarth (1983, MNRAS, 203, 301)"
           print *,"       Values:"
           print *,"          How:  Galactic law of Howarth (1983, MNRAS, 203, 301)"
           print *,"          CCM:  Galactic law of Cardelli, Clayton, Mathis (1989, ApJ, 345, 245)"
           print *,"          Fitz: Galactic law of Fitzpatrick & Massa (1990, ApJS, 72, 163)"
           print *,"          LMC:  LMC law of Howarth (1983, MNRAS, 203, 301)"
           print *,"          SMC:  SMC law of Prevot et al. (984, A&A, 132, 389)"
           print *,"  -c"
           print *,"       The logarithmic extinction at H beta"
           print *,"       Default: calculated from Balmer line ratios"
        !  to be implemented:
        !  -R                     : R (default 3.1)
        !  -nelow / --density-low : low ionisation zone density (default - calculate from line list)
        !  -telow / --temperature-low : low i. zone temperature
        !  -nemed / --density-med : medium i. density
        !  -temed / --temperature-med : medium i. temperature
        !  -nehigh / --density-high : high i. density
        !  -tehigh / --temperature-high : high i. temperature
        !  -c                     : value of c(Hb), the logarithmic extinction at H beta
           stop
        endif

        ALLOCATE (options(Narg))

        do i=1,Narg
                call getarg(i,options(i))
        enddo

        ! set defaults

        runs=1
        switch_ext="S"
        filename=""
        meanextinction=0.

        ! process command line arguments

        do i=1,Narg 
                if ((trim(options(i))=="-n" .or. trim(options(i))=="--n-iterations") .and. (i+1) .le. Narg) then
                   read (options(i+1),*) runs
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
         enddo

         if (Narg .eq. 1) then
           filename=trim(options(1))
         endif

         deallocate(options)

         if (filename=="") then
                print *,"Error: No input file specified"
                stop
         endif

        !first, read in the line list

        print *,"Initialising"
        print *,"------------"

        call DATE_AND_TIME(date,time)
        print *
        print *,"Start time: ",time(1:2),":",time(3:4),":",time(5:10)," on ",date(7:8),"/",date(5:6),"/",date(1:4)
        print *,"Input file: ",filename
        print *,"Iterations: ",runs

        I = 1
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
        END DO
        CLOSE(199)

        110 PRINT "(A9,I4,A15,I4,A9)", "Read in ", I," lines (out of ",listlength," in file)"

        if (I .ne. listlength) then
                print *,"Line list reading failed"
                print *,"This can happen if it doesn't have three columns"
                stop
        endif

        if(linelist(1)%wavelength == 0)then
                PRINT*, "Cheese shop error: no inputs"
                STOP
        endif

        if (switch_ext == "S") then
                print *,"Using Howarth (1983) galactic law"
        elseif (switch_ext == "H") then
                print *,"Using Howarth (1983) LMC law"
        elseif (switch_ext == "C") then
                print *,"Using CCM (1989) galactic law"
        elseif (switch_ext == "P") then
                print *,"Using Prevot et al. (1984) SMC law"
        elseif (switch_ext == "F") then
                print *,"Using Fitzpatrick (1990) galactic law"
        endif

! read the CEL data

        call read_ilines(ILs, Iint)

        !now check number of iterations.  If 1, line list is fine as is.  If more than one, randomize the fluxes

        if(runs == 1)then !calculates abundances without uncertainties
                call abundances(linelist, switch_ext, listlength, filename, iteration_result, R, meanextinction, calculate_extinction, ILs, Iint)

                !generate outputs

                700 FORMAT(X,A20,F5.3) !extinction format
                701 FORMAT(X,A20,I5) !diagnostics format
                702 FORMAT(X,A20,ES14.2) !abundances format
                703 FORMAT(X,A20,F5.2) !strong line format
                704 FORMAT(X,A20,F5.2) !adf format

                print *
                print *,"Extinction"
                print *,"=========="
                print *
                write (*,700) "mean_cHb: ",iteration_result(1)%mean_cHb
                print *
                print *,"Diagnostics"
                print *,"==========="
                print *
                write (*,701) "OII_density: ",INT(iteration_result(1)%OII_density)
                write (*,701) "SII_density: ",INT(iteration_result(1)%SII_density)
                write (*,701) "low_density: ",INT(iteration_result(1)%low_density)
                print *
                write (*,701) "OII_temp: ",INT(iteration_result(1)%OII_temp)
                write (*,701) "NII_temp: ",INT(iteration_result(1)%NII_temp)
                write (*,701) "SII_temp: ",INT(iteration_result(1)%SII_temp)
                write (*,701) "OI_temp: ",INT(iteration_result(1)%OI_temp)
                write (*,701) "CI_temp: ",INT(iteration_result(1)%CI_temp)
                write (*,701) "low_temp: ",INT(iteration_result(1)%low_temp)
                print *
                write (*,701) "ClIII_density: ",INT(iteration_result(1)%ClIII_density)
                write (*,701) "ArIV_density: ",INT(iteration_result(1)%ArIV_density)
                write (*,701) "CIII_density: ",INT(iteration_result(1)%CIII_density)
                write (*,701) "OIII_IR_density: ",INT(iteration_result(1)%OIII_IR_density)
                write (*,701) "SIII_IR_density: ",INT(iteration_result(1)%SIII_IR_density)
                write (*,701) "ArIII_IR_density: ",INT(iteration_result(1)%ArIII_IR_density)
                write (*,701) "NeIII_IR_density: ",INT(iteration_result(1)%NeIII_IR_density)
                write (*,701) "med_density: ",INT(iteration_result(1)%med_density)
                print *
                write (*,701) "OIII_temp: ",INT(iteration_result(1)%OIII_temp)
                write (*,701) "OIII_IR_temp: ",INT(iteration_result(1)%OIII_IR_temp)
                write (*,701) "NeIII_temp: ",INT(iteration_result(1)%NeIII_temp)
                write (*,701) "NeIII_IR_temp: ",INT(iteration_result(1)%NeIII_IR_temp)
                write (*,701) "ArIII_temp: ",INT(iteration_result(1)%ArIII_temp)
                write (*,701) "SIII_temp: ",INT(iteration_result(1)%SIII_temp)
                write (*,701) "med_temp: ",INT(iteration_result(1)%med_temp)
                print *
                write (*,701) "NeIV_density: ",INT(iteration_result(1)%NeIV_density)
                write (*,701) "high_density: ",INT(iteration_result(1)%high_density)
                print *
                write (*,701) "ArV_temp: ",INT(iteration_result(1)%ArV_temp)
                write (*,701) "NeV_temp: ",INT(iteration_result(1)%NeV_temp)
                write (*,701) "high_temp: ",INT(iteration_result(1)%high_temp)
                print *
                print *,"Abundances"
                print *,"=========="
                print *
                print *,"Collisionally excited lines"
                write (*,702) "NC_abund_CEL: ",iteration_result(1)%NC_abund_CEL
                write (*,702) "cii_abund_CEL: ",iteration_result(1)%cii_abund_CEL
                write (*,702) "ciii_abund_CEL: ",iteration_result(1)%ciii_abund_CEL
                write (*,702) "civ_abund_CEL: ",iteration_result(1)%civ_abund_CEL
                write (*,702) "C_abund_CEL: ",iteration_result(1)%C_abund_CEL
                write (*,702) "Nii_abund_CEL: ",iteration_result(1)%Nii_abund_CEL
                write (*,702) "Niii_abund_CEL: ",iteration_result(1)%Niii_abund_CEL
                write (*,702) "Niv_abund_CEL: ",iteration_result(1)%Niv_abund_CEL
                write (*,702) "Nv_abund_CEL: ",iteration_result(1)%Nv_abund_CEL
                write (*,702) "N_abund_CEL: ",iteration_result(1)%N_abund_CEL
                write (*,702) "NO_abund_CEL: ",iteration_result(1)%NO_abund_CEL
                write (*,702) "Oii_abund_CEL: ",iteration_result(1)%Oii_abund_CEL
                write (*,702) "Oiii_abund_CEL: ",iteration_result(1)%Oiii_abund_CEL
                write (*,702) "Oiv_abund_CEL: ",iteration_result(1)%Oiv_abund_CEL
                write (*,702) "O_abund_CEL: ",iteration_result(1)%O_abund_CEL
                write (*,702) "Neii_abund_CEL: ",iteration_result(1)%Neii_abund_CEL
                write (*,702) "Neiii_abund_CEL: ",iteration_result(1)%Neiii_abund_CEL
                write (*,702) "Neiv_abund_CEL: ",iteration_result(1)%Neiv_abund_CEL
                write (*,702) "Nev_abund_CEL: ",iteration_result(1)%Nev_abund_CEL
                write (*,702) "Ne_abund_CEL: ",iteration_result(1)%Ne_abund_CEL
                write (*,702) "Ariii_abund_CEL: ",iteration_result(1)%Ariii_abund_CEL
                write (*,702) "Ariv_abund_CEL: ",iteration_result(1)%Ariv_abund_CEL
                write (*,702) "Arv_abund_CEL: ",iteration_result(1)%Arv_abund_CEL
                write (*,702) "Ar_abund_CEL: ",iteration_result(1)%Ar_abund_CEL
                write (*,702) "Sii_abund_CEL: ",iteration_result(1)%Sii_abund_CEL
                write (*,702) "Siii_abund_CEL: ",iteration_result(1)%Siii_abund_CEL
                write (*,702) "Cliii_abund_CEL: ",iteration_result(1)%Cliii_abund_CEL
                write (*,702) "Cl_abund_CEL: ",iteration_result(1)%Cl_abund_CEL
                write (*,702) "S_abund_CEL: ",iteration_result(1)%S_abund_CEL
                print *
                print *,"Recombination lines"
                print *,"-------------------"
                print *
                write (*,702) "He_abund_ORL: ",iteration_result(1)%He_abund_ORL
                write (*,702) "C_abund_ORL: ",iteration_result(1)%C_abund_ORL
                write (*,702) "N_abund_ORL: ",iteration_result(1)%N_abund_ORL
                write (*,702) "O_abund_ORL: ",iteration_result(1)%O_abund_ORL
                write (*,702) "Ne_abund_ORL: ",iteration_result(1)%Ne_abund_ORL
                print *
                print *,"Strong line methods"
                print *,"-------------------"
                print *
                write (*,703) "O_R23_upper: ",iteration_result(1)%O_R23_upper
                write (*,703) "O_R23_lower: ",iteration_result(1)%O_R23_lower
                write (*,703) "O_N2: ",iteration_result(1)%O_N2
                write (*,703) "O_O3N2: ",iteration_result(1)%O_O3N2
                write (*,703) "O_Ar3O3: ",iteration_result(1)%O_Ar3O3
                write (*,703) "O_S3O3: ",iteration_result(1)%O_S3O3
                print *
                print *,"Abundance discrepancy factors"
                print *,"-----------------------------"
                print *
                write (*,704) "adf_O: ",iteration_result(1)%adf_O
                write (*,704) "adf_O2plus: ",iteration_result(1)%adf_O2plus
                write (*,704) "adf_N: ",iteration_result(1)%adf_N
                write (*,704) "adf_N2plus: ",iteration_result(1)%adf_N2plus
                write (*,704) "adf_C: ",iteration_result(1)%adf_C
                write (*,704) "adf_C2plus: ",iteration_result(1)%adf_C2plus
                write (*,704) "adf_Ne: ",iteration_result(1)%adf_Ne
                write (*,704) "adf_Ne2plus : ",iteration_result(1)%adf_Ne2plus


        else if(runs > 1)then

                !save unrandomised line list

                linelist_original = linelist

                call init_random_seed()!sets seed for randomiser
                allocate(all_results(runs))

                !main loop

                DO I=1,runs 
                        print*, "iteration ", i, "of", runs 

                        call randomizer(linelist, listlength, R)
                        R=3.1 ! no randomisation
                        call abundances(linelist, switch_ext, listlength, filename, iteration_result, R, meanextinction, calculate_extinction, ILs, Iint)
                        linelist = linelist_original
                        all_results(i)=iteration_result(1)
                END DO

                ! now process outputs
                print *, "Processing results..."

                OPEN(841, FILE=trim(filename)//"_NC_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(842, FILE=trim(filename)//"_C_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(843, FILE=trim(filename)//"_Nii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(844, FILE=trim(filename)//"_N_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(845, FILE=trim(filename)//"_NO_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(846, FILE=trim(filename)//"_Oii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(847, FILE=trim(filename)//"_Oiii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(848, FILE=trim(filename)//"_O_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(849, FILE=trim(filename)//"_Neii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(850, FILE=trim(filename)//"_Neiii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(851, FILE=trim(filename)//"_Neiv_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(852, FILE=trim(filename)//"_Nev_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(853, FILE=trim(filename)//"_Ne_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(854, FILE=trim(filename)//"_Ariii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(855, FILE=trim(filename)//"_Ariv_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(856, FILE=trim(filename)//"_Arv_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(857, FILE=trim(filename)//"_Ar_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(858, FILE=trim(filename)//"_Sii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(859, FILE=trim(filename)//"_Siii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(860, FILE=trim(filename)//"_S_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(861, FILE=trim(filename)//"_He_abund_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(862, FILE=trim(filename)//"_C_abund_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(863, FILE=trim(filename)//"_N_abund_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(864, FILE=trim(filename)//"_O_abund_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(865, FILE=trim(filename)//"_Ne_abund_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(866, FILE=trim(filename)//"_[OII]_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(867, FILE=trim(filename)//"_[SII]_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(868, FILE=trim(filename)//"_low_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(869, FILE=trim(filename)//"_[OII]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(870, FILE=trim(filename)//"_[NII]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(871, FILE=trim(filename)//"_[SII]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(872, FILE=trim(filename)//"_[OI]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(873, FILE=trim(filename)//"_[CI]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(874, FILE=trim(filename)//"_low_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(875, FILE=trim(filename)//"_[ClIII]_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(876, FILE=trim(filename)//"_[ArIV]_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(877, FILE=trim(filename)//"_CIII]_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(878, FILE=trim(filename)//"_med_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(879, FILE=trim(filename)//"_[OIII]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(880, FILE=trim(filename)//"_[NeIII]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(881, FILE=trim(filename)//"_[ArIII]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(882, FILE=trim(filename)//"_[SIII]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(883, FILE=trim(filename)//"_med_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(884, FILE=trim(filename)//"_[NeIV]_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(885, FILE=trim(filename)//"_high_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(886, FILE=trim(filename)//"_[ArV]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(887, FILE=trim(filename)//"_[NeV]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(888, FILE=trim(filename)//"_high_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(889, FILE=trim(filename)//"_mean_cHb", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(890, FILE=trim(filename)//"_adf_O", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(891, FILE=trim(filename)//"_adf_O2plus", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(892, FILE=trim(filename)//"_adf_C", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(893, FILE=trim(filename)//"_adf_C2plus", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(894, FILE=trim(filename)//"_adf_N", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(895, FILE=trim(filename)//"_adf_N2plus", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(896, FILE=trim(filename)//"_adf_Ne", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(897, FILE=trim(filename)//"_adf_Ne2plus", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
!XXXX add Cl III, Cl/H, Niii, cii, ciii, ArIII IR dens, NeIII IR dens, strong line, ICF files

! sort all the arrays into ascending order

                call qsort(all_results%NC_abund_CEL)
                call qsort(all_results%Cii_abund_CEL)
                call qsort(all_results%Ciii_abund_CEL)
                call qsort(all_results%C_abund_CEL)
                call qsort(all_results%nii_abund_CEL)
                call qsort(all_results%niii_abund_CEL)
                call qsort(all_results%N_abund_CEL)
                call qsort(all_results%NO_abund_CEL)
                call qsort(all_results%oii_abund_CEL)
                call qsort(all_results%oiii_abund_CEL)
                call qsort(all_results%O_abund_CEL)
                call qsort(all_results%Neii_abund_CEL)
                call qsort(all_results%neiii_abund_CEL)
                call qsort(all_results%neiv_abund_CEL)
                call qsort(all_results%nev_abund_CEL)
                call qsort(all_results%Ne_abund_CEL)
                call qsort(all_results%ariii_abund_CEL)
                call qsort(all_results%ariv_abund_CEL)
                call qsort(all_results%arv_abund_CEL)
                call qsort(all_results%Ar_abund_CEL)
                call qsort(all_results%sii_abund_CEL)
                call qsort(all_results%siii_abund_CEL)
                call qsort(all_results%S_abund_CEL)
                call qsort(all_results%cliii_abund_CEL)
                call qsort(all_results%cl_abund_CEL)
                call qsort(all_results%He_abund_ORL)
                call qsort(all_results%C_abund_ORL)
                call qsort(all_results%N_abund_ORL)
                call qsort(all_results%O_abund_ORL)
                call qsort(all_results%Ne_abund_ORL)
                call qsort(all_results%oii_density)
                call qsort(all_results%sii_density)
                call qsort(all_results%low_density)
                call qsort(all_results%nii_temp)
                call qsort(all_results%oii_temp)
                call qsort(all_results%sii_temp)
                call qsort(all_results%oi_temp)
                call qsort(all_results%ci_temp)
                call qsort(all_results%low_temp)
                call qsort(all_results%cliii_density)
                call qsort(all_results%ariv_density)
                call qsort(all_results%ciii_density)
                call qsort(all_results%ariii_ir_density)
                call qsort(all_results%Neiii_IR_density)
                call qsort(all_results%med_density)
                call qsort(all_results%oiii_temp)
                call qsort(all_results%neiii_temp)
                call qsort(all_results%ariii_temp)
                call qsort(all_results%siii_temp)
                call qsort(all_results%med_temp)
                call qsort(all_results%neiv_density)
                call qsort(all_results%high_density)
                call qsort(all_results%arv_temp)
                call qsort(all_results%nev_temp)
                call qsort(all_results%high_temp)
                call qsort(all_results%mean_cHb)
                call qsort(all_results%adf_C2plus)
                call qsort(all_results%adf_C)
                call qsort(all_results%adf_O2plus)
                call qsort(all_results%adf_O)
                call qsort(all_results%adf_N2plus)
                call qsort(all_results%adf_N)
                call qsort(all_results%adf_Ne2plus)
                call qsort(all_results%adf_Ne)
                call qsort(all_results%O_R23_upper)
                call qsort(all_results%O_R23_lower)
                call qsort(all_results%O_N2)
                call qsort(all_results%O_o3n2)
                call qsort(all_results%O_Ar3O3)
                call qsort(all_results%O_S3O3)


                do i=1,runs
                        write(unit = 841,FMT=*) all_results(i)%NC_abund_CEL
                        write(unit = 842,FMT=*) all_results(i)%C_abund_CEL
                        write(unit = 843,FMT=*) all_results(i)%nii_abund_CEL
                        write(unit = 844,FMT=*) all_results(i)%N_abund_CEL
                        write(unit = 845,FMT=*) all_results(i)%NO_abund_CEL
                        write(unit = 846,FMT=*) all_results(i)%oii_abund_CEL
                        write(unit = 847,FMT=*) all_results(i)%oiii_abund_CEL
                        write(unit = 848,FMT=*) all_results(i)%O_abund_CEL
                        write(unit = 849,FMT=*) all_results(i)%Neii_abund_CEL
                        write(unit = 850,FMT=*) all_results(i)%neiii_abund_CEL
                        write(unit = 851,FMT=*) all_results(i)%neiv_abund_CEL
                        write(unit = 852,FMT=*) all_results(i)%nev_abund_CEL
                        write(unit = 853,FMT=*) all_results(i)%Ne_abund_CEL
                        write(unit = 854,FMT=*) all_results(i)%ariii_abund_CEL
                        write(unit = 855,FMT=*) all_results(i)%ariv_abund_CEL
                        write(unit = 856,FMT=*) all_results(i)%arv_abund_CEL
                        write(unit = 857,FMT=*) all_results(i)%Ar_abund_CEL
                        write(unit = 858,FMT=*) all_results(i)%sii_abund_CEL
                        write(unit = 859,FMT=*) all_results(i)%siii_abund_CEL
                        write(unit = 860,FMT=*) all_results(i)%S_abund_CEL
                        write(unit = 861,FMT=*) all_results(i)%He_abund_ORL
                        write(unit = 862,FMT=*) all_results(i)%C_abund_ORL
                        write(unit = 863,FMT=*) all_results(i)%N_abund_ORL
                        write(unit = 864,FMT=*) all_results(i)%O_abund_ORL
                        write(unit = 865,FMT=*) all_results(i)%Ne_abund_ORL
                        write(unit = 866,FMT=*) all_results(i)%oii_density
                        write(unit = 867,FMT=*) all_results(i)%sii_density
                        write(unit = 868,FMT=*) all_results(i)%low_density
                        write(unit = 869,FMT=*) all_results(i)%nii_temp
                        write(unit = 870,FMT=*) all_results(i)%oii_temp
                        write(unit = 871,FMT=*) all_results(i)%sii_temp
                        write(unit = 872,FMT=*) all_results(i)%oi_temp
                        write(unit = 873,FMT=*) all_results(i)%ci_temp
                        write(unit = 874,FMT=*) all_results(i)%low_temp
                        write(unit = 875,FMT=*) all_results(i)%cliii_density
                        write(unit = 876,FMT=*) all_results(i)%ariv_density
                        write(unit = 877,FMT=*) all_results(i)%ciii_density
                        write(unit = 878,FMT=*) all_results(i)%med_density
                        write(unit = 879,FMT=*) all_results(i)%oiii_temp
                        write(unit = 880,FMT=*) all_results(i)%neiii_temp
                        write(unit = 881,FMT=*) all_results(i)%ariii_temp
                        write(unit = 882,FMT=*) all_results(i)%siii_temp
                        write(unit = 883,FMT=*) all_results(i)%med_temp
                        write(unit = 884,FMT=*) all_results(i)%neiv_density
                        write(unit = 885,FMT=*) all_results(i)%high_density
                        write(unit = 886,FMT=*) all_results(i)%arv_temp
                        write(unit = 887,FMT=*) all_results(i)%nev_temp
                        write(unit = 888,FMT=*) all_results(i)%high_temp
                        write(unit = 889,FMT=*) all_results(i)%mean_cHb
                        write(unit = 890,FMT=*) all_results(i)%adf_O2plus
                        write(unit = 891,FMT=*) all_results(i)%adf_O
                        write(unit = 892,FMT=*) all_results(i)%adf_N2plus
                        write(unit = 893,FMT=*) all_results(i)%adf_N
                        write(unit = 894,FMT=*) all_results(i)%adf_C2plus
                        write(unit = 895,FMT=*) all_results(i)%adf_C
                        write(unit = 896,FMT=*) all_results(i)%adf_Ne2plus
                        write(unit = 897,FMT=*) all_results(i)%adf_Ne
                end do

                DO I=841,897
                        CLOSE(unit=I)
                END DO
! get median +- pseudo gaussian 34.1% and mode
                allocate (quantity_result(runs))
!cHb
                print *,"Extinction"
                print *,"=========="
                print *
                quantity_result = all_results%mean_cHb
                call get_median_mode(quantity_result, median_array, mode)
                print "(X,A,F5.3,A,F5.3,A,F5.3)","c(Hb): median = ",median_array(2)," +",median_array(3),"-",median_array(1) 
                print "(X,A,F5.3)","       mode   = ",mode
print *
print *,"Diagnostics"
print *,"==========="
print *
711 format (X,3(A,I5)) ! diagnostic format
712 format (X,A,I5)! diagnostic format 2 
!low densities

                quantity_result = all_results%oii_density
                call get_median_mode(quantity_result, median_array, mode)
                print 711,"[OII] density: median =",int(median_array(2))," +",int(median_array(3)),"-",int(median_array(1)) 
                print 712,"               mode   = ",int(mode)

                quantity_result = all_results%SII_density
                call get_median_mode(quantity_result, median_array, mode)
                print 711,"[SII] density: median = ",int(median_array(2))," +",int(median_array(3)),"-",int(median_array(1)) 
                print 712,"               mode   = ",int(mode)

                quantity_result = all_results%low_density
                call get_median_mode(quantity_result, median_array, mode)
                print 711,"low density  : median =",int(median_array(2))," +",int(median_array(3)),"-",int(median_array(1)) 
                print 712,"               mode   = ",int(mode)

!low temperatures
print *

                quantity_result = all_results%oii_temp
                call get_median_mode(quantity_result, median_array, mode)
                print 711,"[OII] temperature: median =",int(median_array(2))," +",int(median_array(3)),"-",int(median_array(1)) 
                print 712,"                   mode   = ",int(mode)

                quantity_result = all_results%SII_temp
                call get_median_mode(quantity_result, median_array, mode)
                print 711,"[SII] temperature: median = ",int(median_array(2))," +",int(median_array(3)),"-",int(median_array(1)) 
                print 712,"                   mode   = ",int(mode)

                quantity_result = all_results%NII_temp
                call get_median_mode(quantity_result, median_array, mode)
                print 711,"[NII] temperature: median =",int(median_array(2))," +",int(median_array(3)),"-",int(median_array(1)) 
                print 712,"                   mode   = ",int(mode)

                quantity_result = all_results%OI_temp
                call get_median_mode(quantity_result, median_array, mode)
                print 711,"[OI] temperature:  median =",int(median_array(2))," +",int(median_array(3)),"-",int(median_array(1)) 
                print 712,"                   mode   = ",int(mode)

                quantity_result = all_results%CI_temp
                call get_median_mode(quantity_result, median_array, mode)
                print 711,"[CI] temperature:  median = ",int(median_array(2))," +",int(median_array(3)),"-",int(median_array(1)) 
                print 712,"                   mode   = ",int(mode)

                quantity_result = all_results%low_temp
                call get_median_mode(quantity_result, median_array, mode)
                print 711,"low temperature  : median =",int(median_array(2))," +",int(median_array(3)),"-",int(median_array(1)) 
                print 712,"                   mode   = ",int(mode)

!medium density
print *

                quantity_result = all_results%cliii_density
                call get_median_mode(quantity_result, median_array, mode)
                print 711,"[ClIII] density:    median =",int(median_array(2))," +",int(median_array(3)),"-",int(median_array(1)) 
                print 712,"                    mode   = ",int(mode)

                quantity_result = all_results%ArIV_density
                call get_median_mode(quantity_result, median_array, mode)
                print 711,"[ArIV] density:     median = ",int(median_array(2))," +",int(median_array(3)),"-",int(median_array(1)) 
                print 712,"                    mode   = ",int(mode)

                quantity_result = all_results%CIII_density
                call get_median_mode(quantity_result, median_array, mode)
                print 711,"[CIII] density:     median =",int(median_array(2))," +",int(median_array(3)),"-",int(median_array(1)) 
                print 712,"                    mode   = ",int(mode)

                quantity_result = all_results%OIII_IR_density
                call get_median_mode(quantity_result, median_array, mode)
                print 711,"[OIII] IR density:  median =",int(median_array(2))," +",int(median_array(3)),"-",int(median_array(1)) 
                print 712,"                    mode   = ",int(mode)

                quantity_result = all_results%SIII_IR_density
                call get_median_mode(quantity_result, median_array, mode)
                print 711,"[SIII] IRdensity:   median = ",int(median_array(2))," +",int(median_array(3)),"-",int(median_array(1)) 
                print 712,"                    mode   = ",int(mode)

                quantity_result = all_results%ArIII_IR_density
                call get_median_mode(quantity_result, median_array, mode)
                print 711,"[ArIII] IR density: median =",int(median_array(2))," +",int(median_array(3)),"-",int(median_array(1)) 
                print 712,"                    mode   = ",int(mode)

                quantity_result = all_results%NeIII_IR_density
                call get_median_mode(quantity_result, median_array, mode)
                print 711,"[NeIII] IR density: median =",int(median_array(2))," +",int(median_array(3)),"-",int(median_array(1)) 
                print 712,"                    mode   = ",int(mode)

                quantity_result = all_results%med_density
                call get_median_mode(quantity_result, median_array, mode)
                print 711,"medium density:     median = ",int(median_array(2))," +",int(median_array(3)),"-",int(median_array(1)) 
                print 712,"                    mode   = ",int(mode)

!medium temperature
print *

                quantity_result = all_results%OIII_temp
                call get_median_mode(quantity_result, median_array, mode)
                print 711,"[OIII] temperature:      median =",int(median_array(2))," +",int(median_array(3)),"-",int(median_array(1)) 
                print 712,"                         mode   = ",int(mode)

                quantity_result = all_results%OIII_IR_temp
                call get_median_mode(quantity_result, median_array, mode)
                print 711,"[OIII] IR temperature:   median =",int(median_array(2))," +",int(median_array(3)),"-",int(median_array(1)) 
                print 712,"                         mode   = ",int(mode)

                quantity_result = all_results%NeIII_temp
                call get_median_mode(quantity_result, median_array, mode)
                print 711,"[NeIII] temperature:    median = ",int(median_array(2))," +",int(median_array(3)),"-",int(median_array(1)) 
                print 712,"                        mode   = ",int(mode)

                quantity_result = all_results%NeIII_IR_temp
                call get_median_mode(quantity_result, median_array, mode)
                print 711,"[NeIII] IR temperature: median =",int(median_array(2))," +",int(median_array(3)),"-",int(median_array(1)) 
                print 712,"                        mode   = ",int(mode)

                quantity_result = all_results%ArIII_temp
                call get_median_mode(quantity_result, median_array, mode)
                print 711,"[ArIII] temperature:    median =",int(median_array(2))," +",int(median_array(3)),"-",int(median_array(1)) 
                print 712,"                        mode   = ",int(mode)

                quantity_result = all_results%SIII_temp
                call get_median_mode(quantity_result, median_array, mode)
                print 711,"[SIII] temperature:     median = ",int(median_array(2))," +",int(median_array(3)),"-",int(median_array(1)) 
                print 712,"                        mode   = ",int(mode)

                quantity_result = all_results%med_temp
                call get_median_mode(quantity_result, median_array, mode)
                print 711,"medium temperature:     median =",int(median_array(2))," +",int(median_array(3)),"-",int(median_array(1)) 
                print 712,"                        mode   = ",int(mode)

!high density
print *

                quantity_result = all_results%neiv_density
                call get_median_mode(quantity_result, median_array, mode)
                print 711,"[NeIV] density:  median =",int(median_array(2))," +",int(median_array(3)),"-",int(median_array(1)) 
                print 712,"                 mode   = ",int(mode)

                quantity_result = all_results%high_density
                call get_median_mode(quantity_result, median_array, mode)
                print 711,"high density:    median = ",int(median_array(2))," +",int(median_array(3)),"-",int(median_array(1)) 
                print 712,"                 mode   = ",int(mode)

!high temperature

                quantity_result = all_results%ArV_temp
                call get_median_mode(quantity_result, median_array, mode)
                print 711,"[ArV] temperature: median =",int(median_array(2))," +",int(median_array(3)),"-",int(median_array(1)) 
                print 712,"                   mode   = ",int(mode)

                quantity_result = all_results%NeV_temp
                call get_median_mode(quantity_result, median_array, mode)
                print 711,"[NeV] temperature: median =",int(median_array(2))," +",int(median_array(3)),"-",int(median_array(1)) 
                print 712,"                   mode   = ",int(mode)

                quantity_result = all_results%high_temp
                call get_median_mode(quantity_result, median_array, mode)
                print 711,"High temperature:  median = ",int(median_array(2))," +",int(median_array(3)),"-",int(median_array(1)) 
                print 712,"                   mode   = ",int(mode)

!CEL abundances
print *
print *,"Abundances (collisionally excited lines)"
print *,"========================================"
print *

713 format (X,3(A,ES10.2)) ! abundances formats
714 format (X,A,ES10.2)      ! abundances formats

                quantity_result = all_results%NC_abund_CEL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"[CI] abundance: median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                   mode   = ",mode

                quantity_result = all_results%cii_abund_CEL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"[CII] abundance:  median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                   mode   = ",mode

                quantity_result = all_results%ciii_abund_CEL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"[CIII] abundance:  median = ",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                   mode   = ",mode

                quantity_result = all_results%civ_abund_CEL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"[CIV] abundance  : median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                   mode   = ",mode

                quantity_result = all_results%C_abund_CEL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"C/H abundance: median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                   mode   = ",mode

                quantity_result = all_results%nii_abund_CEL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"[NII] abundance:  median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                   mode   = ",mode

                quantity_result = all_results%niii_abund_CEL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"[NIII] abundance:  median = ",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                   mode   = ",mode

                quantity_result = all_results%niv_abund_CEL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"[NIV] abundance  : median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                   mode   = ",mode

                quantity_result = all_results%nv_abund_CEL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"[NV] abundance: median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                   mode   = ",mode

                quantity_result = all_results%N_abund_CEL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"N/H abundance:  median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                   mode   = ",mode

                quantity_result = all_results%NO_abund_CEL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"[OI] abundance:  median = ",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                   mode   = ",mode

                quantity_result = all_results%Oii_abund_CEL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"[OII] abundance  : median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                   mode   = ",mode

                quantity_result = all_results%Oiii_abund_CEL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"[OIII] abundance: median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                   mode   = ",mode

                quantity_result = all_results%Oiv_abund_CEL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"[OIV] abundance:  median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                   mode   = ",mode

                quantity_result = all_results%O_abund_CEL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"O/H abundance:  median = ",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                   mode   = ",mode

!                quantity_result = all_results%NeII_abund_CEL
!                call get_median_mode(quantity_result, median_array, mode)
!                print 713,"[NeII] abundance  : median =",median_array(2)," +",median_array(3),"-",median_array(1) 
!                print 714,"                   mode   = ",mode

                quantity_result = all_results%NeIII_abund_CEL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"[NeIII] abundance: median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                   mode   = ",mode

                quantity_result = all_results%NeIV_abund_CEL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"[NeIV] abundance:  median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                   mode   = ",mode

                quantity_result = all_results%NeV_abund_CEL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"[NeV] abundance:  median = ",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                   mode   = ",mode

                quantity_result = all_results%Ne_abund_CEL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"Ne/H abundance  : median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                   mode   = ",mode

                quantity_result = all_results%ArIII_abund_CEL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"[ArIII] abundance: median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                   mode   = ",mode

                quantity_result = all_results%ArIV_abund_CEL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"[ArIV] abundance:  median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                   mode   = ",mode

                quantity_result = all_results%ArV_abund_CEL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"[ArV] abundance:  median = ",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                   mode   = ",mode

                quantity_result = all_results%Ar_abund_CEL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"Ar/H abundance  : median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                   mode   = ",mode

                quantity_result = all_results%SII_abund_CEL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"[SII] abundance: median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                   mode   = ",mode

                quantity_result = all_results%SIII_abund_CEL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"[SIII] abundance:  median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                   mode   = ",mode

                quantity_result = all_results%S_abund_CEL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"S/H abundance:  median = ",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                   mode   = ",mode

                quantity_result = all_results%ClIII_abund_CEL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"[ClIII] abundance  : median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                   mode   = ",mode

                quantity_result = all_results%Cl_abund_CEL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"Cl/H abundance: median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                   mode   = ",mode

!RL abundances
print *
print *,"Abundances (recombination lines)"
print *,"================================"
print *

                quantity_result = all_results%He_abund_ORL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"He/H abundance:  median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                 mode   = ",mode

                quantity_result = all_results%C_abund_ORL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"C/H abundance:   median = ",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                 mode   = ",mode

                quantity_result = all_results%N_abund_ORL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"N/H abundance:  median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                mode   = ",mode

                quantity_result = all_results%O_abund_ORL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"O/H abundance:  median = ",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                mode   = ",mode

                quantity_result = all_results%Ne_abund_ORL
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"Ne/H abundance: median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                mode   = ",mode

!strong line abundances
print *
print *,"Abundances (strong line methods)"
print *,"================================"
print *

                quantity_result = all_results%O_R23_upper
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"O/H (R23 upper): median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                 mode   = ",mode

                quantity_result = all_results%O_R23_lower
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"O/H (R23 lower): median = ",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                 mode   = ",mode

                quantity_result = all_results%O_N2
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"O/H (N2):       median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                mode   = ",mode

                quantity_result = all_results%O_O3N2
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"O/H (O3N2):     median = ",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                mode   = ",mode

                quantity_result = all_results%O_Ar3O3
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"O/H (Ar3O3):    median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                mode   = ",mode

                quantity_result = all_results%O_S3O3
                call get_median_mode(quantity_result, median_array, mode)
                print 713,"O/H (S3O3):     median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 714,"                mode   = ",mode

!adfs
print *
print *,"Abundance discrepancy factors"
print *,"============================="
print *

715 format (X,3(A,F5.2))
716 format (X,A,F5.2)

                quantity_result = all_results%adf_o2plus
                call get_median_mode(quantity_result, median_array, mode)
                print 715,"adf(O2+/H):   median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 716,"              mode   = ",mode

                quantity_result = all_results%adf_o
                call get_median_mode(quantity_result, median_array, mode)
                print 715,"adf(O/H+):    median = ",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 716,"              mode   = ",mode

                quantity_result = all_results%adf_n2plus
                call get_median_mode(quantity_result, median_array, mode)
                print 715,"adf(N2+/H):   median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 716,"              mode   = ",mode

                quantity_result = all_results%adf_n
                call get_median_mode(quantity_result, median_array, mode)
                print 715,"adf(N/H):     median = ",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 716,"              mode   = ",mode

                quantity_result = all_results%adf_c2plus
                call get_median_mode(quantity_result, median_array, mode)
                print 715,"adf(C2+/H+):  median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 716,"              mode   = ",mode

                quantity_result = all_results%adf_c
                call get_median_mode(quantity_result, median_array, mode)
                print 715,"adf(C/H):     median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 716,"              mode   = ",mode

                quantity_result = all_results%adf_ne2plus
                call get_median_mode(quantity_result, median_array, mode)
                print 715,"adf(Ne2+/H+): median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 716,"            : mode   = ",mode

                quantity_result = all_results%adf_ne
                call get_median_mode(quantity_result, median_array, mode)
                print 715,"adf(Ne/H):    median =",median_array(2)," +",median_array(3),"-",median_array(1) 
                print 716,"              mode   = ",mode

        else
                print*, "I didn't want to be a barber anyway. I wanted to be... a lumberjack!   Also, a positive number of runs helps.."
        endif

        call DATE_AND_TIME(date,time)
        print *
        print *,"End time:   ",time(1:2),":",time(3:4),":",time(5:10)," on ",date(7:8),"/",date(5:6),"/",date(1:4)

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

SUBROUTINE deredden_ll(switch_ext, linelist, listlength, meanextinction )
        INTEGER :: iii, listlength
        CHARACTER :: switch_ext !switch for extinction laws
        TYPE(LINE),DIMENSION(:), allocatable :: linelist
        double precision :: meanextinction

        if (switch_ext == "S") then
                CALL deredden(linelist, listlength, meanextinction)
        elseif (switch_ext == "H") then
                CALL deredden_LMC(linelist, listlength, meanextinction)
        elseif (switch_ext == "C") then
                CALL deredden_CCM(linelist, listlength, meanextinction, R)
        elseif (switch_ext == "P") then
                CALL deredden_SMC(linelist, listlength, meanextinction)
        elseif (switch_ext == "F") then
                CALL deredden_Fitz(linelist, listlength, meanextinction)
        endif


        500 FORMAT (5(f10.4))

        OPEN(801, FILE=trim(filename)//"_dered", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
        do iii=1, listlength
                if(linelist(iii)%int_dered .ne. 0)then
                        write(801,500) linelist(iii)%wavelength, linelist(iii)%intensity, linelist(iii)%int_err, linelist(iii)%int_dered
                endif
        end do
               CLOSE(801)
               call system("sort "//trim(filename)//"_dered > "//trim(filename)//"_dered_sort")
               call system("rm "//trim(filename)//"_dered")


END SUBROUTINE

subroutine get_median_mode(input_array, median_array, mode)
double precision, intent(in) :: input_array(:)
double precision, intent(out) :: median_array(3)
double precision :: binsize
double precision, intent(out) :: mode
double precision, dimension (:,:), allocatable :: binned_quantity_result
integer :: ii, bincount, bincountmax, arraysize

mode=0.D0
median_array = (/0.0,0.0,0.0/)

arraysize = size(input_array)

median_array(1) = input_array(int(0.500*size(input_array))) - input_array(int(0.159*size(input_array)))
median_array(2) = input_array(int(0.500*size(input_array)))
median_array(3) = input_array(int(0.841*size(input_array))) - input_array(int(0.500*size(input_array)))

binsize = (median_array(3)+median_array(1))/20

if (binsize .gt. 0) then
!  allocate(binned_quantity_result(1+(int(quantity_result(arraysize)/binsize)-int(quantity_result(1)/binsize)),2)) 
  allocate(binned_quantity_result(arraysize,2))
  binned_quantity_result = 0.D0

  ii=1
  bincount=1 !(why does this need to be one and not zero??)
  binvalue = int(quantity_result(1)/binsize)
  mode=0.D0
  bincountmax=0

  do i=1,runs 
    if (int(quantity_result(i)/binsize) == binvalue) then 
      bincount = bincount + 1
    else 
      binned_quantity_result(ii,1) = binvalue*binsize
      binned_quantity_result(ii,2) = bincount 

      if (bincount>bincountmax) then
        bincountmax=bincount
        mode = binned_quantity_result(ii,1)
      endif

      ii=ii+1
      bincount = 1
      binvalue = binvalue + 1
    endif
  enddo

  deallocate(binned_quantity_result)
endif

end subroutine get_median_mode

end program
