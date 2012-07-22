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

        CHARACTER :: switch_ext !switch for extinction laws
        INTEGER :: I, runs, Narg !runs = number of runs for randomiser

        !input options

        CHARACTER*2048, DIMENSION(:), allocatable :: options
        CHARACTER*2048 :: commandline

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

        !atomic data

        character*10 :: ionlist(40) !list of ion names
        integer :: iion !# of ions in Ilines
        integer :: maxlevs,maxtemps
        type(atomic_data),allocatable :: atomicdata(:)

        !extinction

        logical :: calculate_extinction=.true.
        DOUBLE PRECISION :: temp1,temp2,temp3, meanextinction, R

        !binning

        double precision, dimension(3) :: uncertainty_array

        !CEL array

        TYPE(line), DIMENSION(:), allocatable :: ILs

        !diagnostic array

        double precision, dimension(6) :: diagnostic_array

        R=3.1

        !read command line arguments

        Narg = IARGC() !count input arguments

        if (Narg .eq. 0) then
           print *,"Syntax: ./neat.exe [option1 value1] [option2 value2] .. [optionx valuex]"
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
           print *,"  -nelow, -nemed, -nehigh"
           print *,"  -telow, -temed, -tehigh"
           print *,"       The electron densities and temperatures to be used."
           print *,"       Units: cm-3 for densities, K for temperatures"
           print *,"       Default: calculated from available diagnostics."
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
        switch_ext="S"
        filename=""
        meanextinction=0.D0
        diagnostic_array=0.D0

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
         enddo

         if (Narg .eq. 1) then
           filename=trim(options(1))
         endif

         if (filename=="") then
                print *,"Error: No input file specified"
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

        if (switch_ext == "S") then
                print *,gettime(), ": using Howarth (1983) galactic law"
        elseif (switch_ext == "H") then
                print *,gettime(), ": using Howarth (1983) LMC law"
        elseif (switch_ext == "C") then
                print *,gettime(), ": using CCM (1989) galactic law"
        elseif (switch_ext == "P") then
                print *,gettime(), ": using Prevot et al. (1984) SMC law"
        elseif (switch_ext == "F") then
                print *,gettime(), ": using Fitzpatrick (1990) galactic law"
        endif

! read the CEL data

        call read_ilines(ILs, Iint,iion,ionlist)

        print *,gettime(), ": Reading atomic data"

        !CELs read, can now read in atomic data

        allocate(atomicdata(iion))
        DO I=1,iion
            atomicdata(I)%ion = ionlist(I)
            call read_atomic_data(atomicdata(I))
        ENDDO

        !read ORL data
        call read_orl_data

        !find maximum #levels and temperatures - pass to equib to reduce footprint

        maxlevs = atomicdata(1)%nlevs
        maxtemps = atomicdata(1)%ntemps
        do i=2,iion
            if(atomicdata(i)%nlevs .gt. maxlevs) maxlevs = atomicdata(i)%nlevs
            if(atomicdata(i)%ntemps .gt. maxtemps) maxtemps = atomicdata(i)%ntemps
        enddo
!        print*,maxlevs,maxtemps
        !now check number of iterations.  If 1, line list is fine as is.  If more than one, randomize the fluxes

        if(runs == 1)then !calculates abundances without uncertainties
                call abundances(linelist, switch_ext, listlength, iteration_result, R, meanextinction, calculate_extinction, ILs, Iint, diagnostic_array,iion,atomicdata,maxlevs,maxtemps)

                !generate outputs

                print *,gettime(),": writing results to file ",trim(filename),"_results"
                OPEN(650, FILE=trim(filename)//"_results", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')

                700 FORMAT(X,A20,F5.3) !extinction format
                701 FORMAT(X,A20,I5) !diagnostics format
                705 FORMAT(X,A20,I5,5X,F8.3) !diagnostics and ratio format
                702 FORMAT(X,A20,ES14.3) !abundances format
                703 FORMAT(X,A20,F5.2) !strong line format
                704 FORMAT(X,A20,F5.2) !adf format

                write (650,*) "NEAT (nebular empirical analysis tool)"
                write (650,*) "======================================"
                write (650,*)
                write (650,*) "Analysis of file ",trim(filename)
                write (650,*)

                write (650,*),"Extinction"
                write (650,*),"=========="
                write (650,*)
                write (650,700) "mean_cHb :         ",iteration_result(1)%mean_cHb
                write (650,*)
                write (650,*),"Diagnostics"
                write (650,*),"==========="
                write (650,*)
                write (650,705) "OII_density :      ",INT(iteration_result(1)%OII_density), iteration_result(1)%OII_density_ratio
                write (650,705) "SII_density :      ",INT(iteration_result(1)%SII_density), iteration_result(1)%SII_density_ratio
                write (650,701) "low_density :      ",INT(iteration_result(1)%low_density)
                write (650,*)
                write (650,705) "OII_temp :         ",INT(iteration_result(1)%OII_temp), iteration_result(1)%OII_temp_ratio
                write (650,705) "NII_temp :         ",INT(iteration_result(1)%NII_temp), iteration_result(1)%NII_temp_ratio
                write (650,705) "SII_temp :         ",INT(iteration_result(1)%SII_temp), iteration_result(1)%SII_temp_ratio
                write (650,705) "OI_temp :          ",INT(iteration_result(1)%OI_temp), iteration_result(1)%OI_temp_ratio
                write (650,705) "CI_temp :          ",INT(iteration_result(1)%CI_temp), iteration_result(1)%CI_temp_ratio
                write (650,701) "low_temp :         ",INT(iteration_result(1)%low_temp)
                write (650,*)
                write (650,705) "ClIII_density :    ",INT(iteration_result(1)%ClIII_density), iteration_result(1)%ClIII_density_ratio
                write (650,705) "ArIV_density :     ",INT(iteration_result(1)%ArIV_density), iteration_result(1)%ArIV_density_ratio
                write (650,705) "CIII_density :     ",INT(iteration_result(1)%CIII_density), iteration_result(1)%CIII_density_ratio
                write (650,705) "OIII_IR_density :  ",INT(iteration_result(1)%OIII_IR_density), iteration_result(1)%OIII_IR_density_ratio
                write (650,705) "SIII_IR_density :  ",INT(iteration_result(1)%SIII_IR_density), iteration_result(1)%SIII_IR_density_ratio
                write (650,705) "ArIII_IR_density : ",INT(iteration_result(1)%ArIII_IR_density), iteration_result(1)%ArIII_IR_density_ratio
                write (650,705) "NeIII_IR_density : ",INT(iteration_result(1)%NeIII_IR_density), iteration_result(1)%NeIII_IR_density_ratio
                write (650,701) "med_density :      ",INT(iteration_result(1)%med_density)
                write (650,*)
                write (650,705) "OIII_temp :        ",INT(iteration_result(1)%OIII_temp), iteration_result(1)%OIII_temp_ratio
                write (650,705) "OIII_IR_temp :     ",INT(iteration_result(1)%OIII_IR_temp), iteration_result(1)%OIII_IR_temp_ratio
                write (650,705) "NeIII_temp :       ",INT(iteration_result(1)%NeIII_temp), iteration_result(1)%NeIII_temp_ratio
                write (650,705) "NeIII_IR_temp :    ",INT(iteration_result(1)%NeIII_IR_temp), iteration_result(1)%NeIII_IR_temp_ratio
                write (650,705) "ArIII_temp :       ",INT(iteration_result(1)%ArIII_temp), iteration_result(1)%ArIII_temp_ratio
                write (650,705) "SIII_temp :        ",INT(iteration_result(1)%SIII_temp), iteration_result(1)%SIII_temp_ratio
                write (650,701) "med_temp :         ",INT(iteration_result(1)%med_temp)
                write (650,*)
                write (650,705) "NeIV_density :     ",INT(iteration_result(1)%NeIV_density), iteration_result(1)%NeIV_density_ratio
                write (650,701) "high_density :     ",INT(iteration_result(1)%high_density)
                write (650,*)
                write (650,705) "ArV_temp :         ",INT(iteration_result(1)%ArV_temp), iteration_result(1)%ArV_temp_ratio
                write (650,705) "NeV_temp :         ",INT(iteration_result(1)%NeV_temp), iteration_result(1)%NeV_temp_ratio
                write (650,701) "high_temp :        ",INT(iteration_result(1)%high_temp)
                write (650,701) "Balmer_jump_temp : ",INT(iteration_result(1)%Bal_jump_temp)
                write (650,*)
                write (650,*),"Abundances"
                write (650,*),"=========="
                write (650,*)
                write (650,*),"Collisionally excited lines"
                write (650,*),"-------------------"
                write (650,*)
                write (650,702) "NC_abund_CEL :     ",iteration_result(1)%NC_abund_CEL
                write (650,702) "cii_abund_CEL :    ",iteration_result(1)%cii_abund_CEL
                write (650,702) "ciii_abund_CEL :   ",iteration_result(1)%ciii_abund_CEL
                write (650,702) "civ_abund_CEL :    ",iteration_result(1)%civ_abund_CEL
                write (650,704) "C_icf_CEL :        ",iteration_result(1)%C_icf_CEL
                write (650,702) "C_abund_CEL :      ",iteration_result(1)%C_abund_CEL
                write (650,702) "Nii_abund_CEL :    ",iteration_result(1)%Nii_abund_CEL
                write (650,702) "Niii_abund_CEL :   ",iteration_result(1)%Niii_abund_CEL
                write (650,702) "Niv_abund_CEL :    ",iteration_result(1)%Niv_abund_CEL
                write (650,702) "Nv_abund_CEL :     ",iteration_result(1)%Nv_abund_CEL
                write (650,704) "N_icf_CEL :        ",iteration_result(1)%N_icf_CEL
                write (650,702) "N_abund_CEL :      ",iteration_result(1)%N_abund_CEL
                write (650,702) "NO_abund_CEL :     ",iteration_result(1)%NO_abund_CEL
                write (650,702) "Oii_abund_CEL :    ",iteration_result(1)%Oii_abund_CEL
                write (650,702) "Oiii_abund_CEL :   ",iteration_result(1)%Oiii_abund_CEL
                write (650,702) "Oiv_abund_CEL :    ",iteration_result(1)%Oiv_abund_CEL
                write (650,704) "O_icf_CEL :        ",iteration_result(1)%O_icf_CEL
                write (650,702) "O_abund_CEL :      ",iteration_result(1)%O_abund_CEL
                write (650,702) "Neii_abund_CEL :   ",iteration_result(1)%Neii_abund_CEL
                write (650,702) "Neiii_abund_CEL :  ",iteration_result(1)%Neiii_abund_CEL
                write (650,702) "Neiv_abund_CEL :   ",iteration_result(1)%Neiv_abund_CEL
                write (650,702) "Nev_abund_CEL :    ",iteration_result(1)%Nev_abund_CEL
                write (650,704) "Ne_icf_CEL :       ",iteration_result(1)%Ne_icf_CEL
                write (650,702) "Ne_abund_CEL :     ",iteration_result(1)%Ne_abund_CEL
                write (650,702) "Ariii_abund_CEL :  ",iteration_result(1)%Ariii_abund_CEL
                write (650,702) "Ariv_abund_CEL :   ",iteration_result(1)%Ariv_abund_CEL
                write (650,702) "Arv_abund_CEL :    ",iteration_result(1)%Arv_abund_CEL
                write (650,704) "Ar_icf_CEL :       ",iteration_result(1)%Ar_icf_CEL
                write (650,702) "Ar_abund_CEL :     ",iteration_result(1)%Ar_abund_CEL
                write (650,702) "Sii_abund_CEL :    ",iteration_result(1)%Sii_abund_CEL
                write (650,702) "Siii_abund_CEL :   ",iteration_result(1)%Siii_abund_CEL
                write (650,704) "S_icf_CEL :        ",iteration_result(1)%S_icf_CEL
                write (650,702) "S_abund_CEL :      ",iteration_result(1)%S_abund_CEL
                write (650,702) "Cliii_abund_CEL :  ",iteration_result(1)%Cliii_abund_CEL
                write (650,704) "Cl_icf_CEL :       ",iteration_result(1)%Cl_icf_CEL
                write (650,702) "Cl_abund_CEL :     ",iteration_result(1)%Cl_abund_CEL
                write (650,*)
                write (650,*),"Recombination lines"
                write (650,*),"-------------------"
                write (650,*)
                write (650,702) "Hei_abund_ORL :    ",iteration_result(1)%Hei_abund_ORL
                write (650,702) "Heii_abund_ORL :   ",iteration_result(1)%Heii_abund_ORL
                write (650,702) "He_abund_ORL :     ",iteration_result(1)%He_abund_ORL
                write (650,702) "Cii_abund_ORL :    ",iteration_result(1)%Cii_abund_ORL
                write (650,702) "Ciii_abund_ORL :   ",iteration_result(1)%Ciii_abund_ORL
                write (650,704) "C_icf_ORL :        ",iteration_result(1)%C_icf_ORL
                write (650,702) "C_abund_ORL :      ",iteration_result(1)%C_abund_ORL
                write (650,702) "Nii_abund_ORL :    ",iteration_result(1)%Nii_abund_ORL
                write (650,702) "Niii_abund_ORL :   ",iteration_result(1)%Niii_abund_ORL
                write (650,704) "N_icf_ORL :        ",iteration_result(1)%N_icf_ORL
                write (650,702) "N_abund_ORL :      ",iteration_result(1)%N_abund_ORL
                write (650,702) "Oii_abund_ORL :    ",iteration_result(1)%Oii_abund_ORL
                write (650,704) "O_icf_ORL :        ",iteration_result(1)%O_icf_ORL
                write (650,702) "O_abund_ORL :      ",iteration_result(1)%O_abund_ORL
                write (650,702) "Neii_abund_ORL :   ",iteration_result(1)%Neii_abund_ORL
                write (650,704) "Ne_icf_ORL :       ",iteration_result(1)%Ne_icf_ORL
                write (650,702) "Ne_abund_ORL :     ",iteration_result(1)%Ne_abund_ORL
                write (650,*)
                write (650,*),"Strong line methods"
                write (650,*),"-------------------"
                write (650,*)
                write (650,703) "O_R23_upper :      ",iteration_result(1)%O_R23_upper
                write (650,703) "O_R23_lower :      ",iteration_result(1)%O_R23_lower
                write (650,703) "O_N2 :             ",iteration_result(1)%O_N2
                write (650,703) "O_O3N2 :           ",iteration_result(1)%O_O3N2
                write (650,703) "O_Ar3O3 :          ",iteration_result(1)%O_Ar3O3
                write (650,703) "O_S3O3 :           ",iteration_result(1)%O_S3O3
                write (650,*)
                write (650,*),"Abundance discrepancy factors"
                write (650,*),"-----------------------------"
                write (650,*)
                write (650,704) "adf_O :            ",iteration_result(1)%adf_O
                write (650,704) "adf_O2plus :       ",iteration_result(1)%adf_O2plus
                write (650,704) "adf_N :            ",iteration_result(1)%adf_N
                write (650,704) "adf_N2plus :       ",iteration_result(1)%adf_N2plus
                write (650,704) "adf_C :            ",iteration_result(1)%adf_C
                write (650,704) "adf_C2plus :       ",iteration_result(1)%adf_C2plus
                write (650,704) "adf_Ne :           ",iteration_result(1)%adf_Ne
                write (650,704) "adf_Ne2plus :      ",iteration_result(1)%adf_Ne2plus

                close (650)

        else if(runs > 1)then

                !save unrandomised line list

                linelist_original = linelist

                call init_random_seed()!sets seed for randomiser
                allocate(all_results(runs))

                !main loop

                print *
                print "(X,A9,X,A)",gettime(), ": starting Monte Carlo calculations"
                print *,gettime(), ": completed ",0,"%"
                DO I=1,runs
                         if ( (10.0*dble(i)/dble(runs)) == int(10*i/runs) ) print *,gettime(),": completed ",100*i/runs,"%"
!                        print*, "iteration ", i, "of", runs

                        call randomizer(linelist, listlength, R)
                        R=3.1 ! no randomisation
                        call abundances(linelist, switch_ext, listlength, iteration_result, R, meanextinction, calculate_extinction, ILs, Iint, diagnostic_array,iion,atomicdata,maxlevs,maxtemps)
                        linelist = linelist_original
                        all_results(i)=iteration_result(1)
                END DO

                ! now process outputs
                print *, gettime(), ": processing results"

                OPEN(841, FILE=trim(filename)//"_NC_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(842, FILE=trim(filename)//"_Cii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(843, FILE=trim(filename)//"_Ciii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(844, FILE=trim(filename)//"_Civ_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(845, FILE=trim(filename)//"_C_icf_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(846, FILE=trim(filename)//"_C_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(847, FILE=trim(filename)//"_Nii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(848, FILE=trim(filename)//"_N_icf_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(849, FILE=trim(filename)//"_N_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(850, FILE=trim(filename)//"_NO_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(851, FILE=trim(filename)//"_Oii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(852, FILE=trim(filename)//"_Oiii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(853, FILE=trim(filename)//"_O_icf_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(854, FILE=trim(filename)//"_O_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(855, FILE=trim(filename)//"_Neii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(856, FILE=trim(filename)//"_Neiii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(857, FILE=trim(filename)//"_Neiv_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(858, FILE=trim(filename)//"_Nev_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(859, FILE=trim(filename)//"_Ne_icf_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(860, FILE=trim(filename)//"_Ne_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(861, FILE=trim(filename)//"_Ariii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(862, FILE=trim(filename)//"_Ariv_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(863, FILE=trim(filename)//"_Arv_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(864, FILE=trim(filename)//"_Ar_icf_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(865, FILE=trim(filename)//"_Ar_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(866, FILE=trim(filename)//"_Sii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(867, FILE=trim(filename)//"_Siii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(868, FILE=trim(filename)//"_S_icf_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(869, FILE=trim(filename)//"_S_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(870, FILE=trim(filename)//"_Cliii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(871, FILE=trim(filename)//"_Cl_icf_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(872, FILE=trim(filename)//"_Cl_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(873, FILE=trim(filename)//"_Hei_abund_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(874, FILE=trim(filename)//"_Heii_abund_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(875, FILE=trim(filename)//"_He_abund_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(876, FILE=trim(filename)//"_Cii_abund_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(877, FILE=trim(filename)//"_Ciii_abund_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(878, FILE=trim(filename)//"_C_icf_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(879, FILE=trim(filename)//"_C_abund_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(880, FILE=trim(filename)//"_Nii_abund_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(881, FILE=trim(filename)//"_Niii_abund_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(882, FILE=trim(filename)//"_N_icf_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(883, FILE=trim(filename)//"_N_abund_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(884, FILE=trim(filename)//"_Oii_abund_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(885, FILE=trim(filename)//"_O_icf_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(886, FILE=trim(filename)//"_O_abund_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(887, FILE=trim(filename)//"_Neii_abund_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(888, FILE=trim(filename)//"_Ne_icf_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(889, FILE=trim(filename)//"_Ne_abund_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(890, FILE=trim(filename)//"_[OII]_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(891, FILE=trim(filename)//"_[SII]_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(892, FILE=trim(filename)//"_low_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(893, FILE=trim(filename)//"_[OII]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(894, FILE=trim(filename)//"_[NII]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(895, FILE=trim(filename)//"_[SII]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(896, FILE=trim(filename)//"_[OI]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(897, FILE=trim(filename)//"_[CI]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(898, FILE=trim(filename)//"_low_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(899, FILE=trim(filename)//"_[ClIII]_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(900, FILE=trim(filename)//"_[ArIV]_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(901, FILE=trim(filename)//"_CIII]_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(902, FILE=trim(filename)//"_[ArIII]_IR_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(903, FILE=trim(filename)//"_[NeIII]_IR_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(904, FILE=trim(filename)//"_med_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(905, FILE=trim(filename)//"_[OIII]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(906, FILE=trim(filename)//"_[NeIII]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(907, FILE=trim(filename)//"_[ArIII]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(908, FILE=trim(filename)//"_[SIII]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(909, FILE=trim(filename)//"_med_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(910, FILE=trim(filename)//"_[NeIV]_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(911, FILE=trim(filename)//"_high_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(912, FILE=trim(filename)//"_[ArV]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(913, FILE=trim(filename)//"_[NeV]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(914, FILE=trim(filename)//"_high_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(915, FILE=trim(filename)//"_mean_cHb", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(916, FILE=trim(filename)//"_adf_O", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(917, FILE=trim(filename)//"_adf_O2plus", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(918, FILE=trim(filename)//"_adf_C", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(919, FILE=trim(filename)//"_adf_C2plus", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(920, FILE=trim(filename)//"_adf_N", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(921, FILE=trim(filename)//"_adf_N2plus", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(922, FILE=trim(filename)//"_adf_Ne", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(923, FILE=trim(filename)//"_adf_Ne2plus", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(924, FILE=trim(filename)//"_O_R23_upper", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(925, FILE=trim(filename)//"_O_R23_lower", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(926, FILE=trim(filename)//"_O_N2", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(927, FILE=trim(filename)//"_O_O3N2", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(928, FILE=trim(filename)//"_O_Ar3O3", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(929, FILE=trim(filename)//"_O_S3O3", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(930, FILE=trim(filename)//"_OII_density_ratio", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(931, FILE=trim(filename)//"_SII_density_ratio", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(932, FILE=trim(filename)//"_NII_temp_ratio", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(933, FILE=trim(filename)//"_OII_temp_ratio", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(934, FILE=trim(filename)//"_SII_temp_ratio", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(935, FILE=trim(filename)//"_CI_temp_ratio", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(936, FILE=trim(filename)//"_OI_temp_ratio", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(937, FILE=trim(filename)//"_ClIII_density_ratio", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(938, FILE=trim(filename)//"_ArIV_density_ratio", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(939, FILE=trim(filename)//"_CIII_density_ratio", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(940, FILE=trim(filename)//"_OIII_IR_density_ratio", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(941, FILE=trim(filename)//"_ArIII_IR_density_ratio", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(942, FILE=trim(filename)//"_SIII_IR_density_ratio", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(943, FILE=trim(filename)//"_NeIII_IR_density_ratio", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(944, FILE=trim(filename)//"_OIII_temp_ratio", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(945, FILE=trim(filename)//"_NeIII_temp_ratio", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(946, FILE=trim(filename)//"_ArIII_temp_ratio", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(947, FILE=trim(filename)//"_SIII_temp_ratio", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(948, FILE=trim(filename)//"_OIII_IR_temp_ratio", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(949, FILE=trim(filename)//"_NeIII_IR_temp_ratio", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(950, FILE=trim(filename)//"_NeIV_density_ratio", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(951, FILE=trim(filename)//"_ArV_temp_ratio", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(952, FILE=trim(filename)//"_NeV_temp_ratio", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(953, FILE=trim(filename)//"_Bal_jump_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')

!XXXX add Cl/H, Niii, cii, ciii, ArIII IR dens, NeIII IR dens, strong line, ICF files

! sort all the arrays into ascending order

                call qsort(all_results%NC_abund_CEL)
                call qsort(all_results%Cii_abund_CEL)
                call qsort(all_results%Ciii_abund_CEL)
                call qsort(all_results%C_icf_CEL)
                call qsort(all_results%C_abund_CEL)
                call qsort(all_results%nii_abund_CEL)
                call qsort(all_results%niii_abund_CEL)
                call qsort(all_results%N_icf_CEL)
                call qsort(all_results%N_abund_CEL)
                call qsort(all_results%NO_abund_CEL)
                call qsort(all_results%oii_abund_CEL)
                call qsort(all_results%oiii_abund_CEL)
                call qsort(all_results%O_icf_CEL)
                call qsort(all_results%O_abund_CEL)
                call qsort(all_results%Neii_abund_CEL)
                call qsort(all_results%neiii_abund_CEL)
                call qsort(all_results%neiv_abund_CEL)
                call qsort(all_results%nev_abund_CEL)
                call qsort(all_results%Ne_icf_CEL)
                call qsort(all_results%Ne_abund_CEL)
                call qsort(all_results%ariii_abund_CEL)
                call qsort(all_results%ariv_abund_CEL)
                call qsort(all_results%arv_abund_CEL)
                call qsort(all_results%Ar_icf_CEL)
                call qsort(all_results%Ar_abund_CEL)
                call qsort(all_results%sii_abund_CEL)
                call qsort(all_results%siii_abund_CEL)
                call qsort(all_results%S_icf_CEL)
                call qsort(all_results%S_abund_CEL)
                call qsort(all_results%cliii_abund_CEL)
                call qsort(all_results%Cl_icf_CEL)
                call qsort(all_results%cl_abund_CEL)
                call qsort(all_results%Hei_abund_ORL)
                call qsort(all_results%Heii_abund_ORL)
                call qsort(all_results%He_abund_ORL)
                call qsort(all_results%Cii_abund_ORL)
                call qsort(all_results%Ciii_abund_ORL)
                call qsort(all_results%C_icf_ORL)
                call qsort(all_results%C_abund_ORL)
                call qsort(all_results%Nii_abund_ORL)
                call qsort(all_results%Niii_abund_ORL)
                call qsort(all_results%N_icf_ORL)
                call qsort(all_results%N_abund_ORL)
                call qsort(all_results%Oii_abund_ORL)
                call qsort(all_results%O_icf_ORL)
                call qsort(all_results%O_abund_ORL)
                call qsort(all_results%Neii_abund_ORL)
                call qsort(all_results%Ne_icf_ORL)
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
                call qsort(all_results%OII_density_ratio)
                call qsort(all_results%SII_density_ratio)
                call qsort(all_results%NII_temp_ratio)
                call qsort(all_results%OII_temp_ratio)
                call qsort(all_results%SII_temp_ratio)
                call qsort(all_results%CI_temp_ratio)
                call qsort(all_results%OI_temp_ratio)
                call qsort(all_results%ClIII_density_ratio)
                call qsort(all_results%ArIV_density_ratio)
                call qsort(all_results%CIII_density_ratio)
                call qsort(all_results%OIII_IR_density_ratio)
                call qsort(all_results%ArIII_IR_density_ratio)
                call qsort(all_results%SIII_IR_density_ratio)
                call qsort(all_results%NeIII_IR_density_ratio)
                call qsort(all_results%OIII_temp_ratio)
                call qsort(all_results%NeIII_temp_ratio)
                call qsort(all_results%ArIII_temp_ratio)
                call qsort(all_results%SIII_temp_ratio)
                call qsort(all_results%OIII_IR_temp_ratio)
                call qsort(all_results%NeIII_IR_temp_ratio)
                call qsort(all_results%NeIV_density_ratio)
                call qsort(all_results%ArV_temp_ratio)
                call qsort(all_results%NeV_temp_ratio)
                call qsort(all_results%Bal_jump_temp)

print *, gettime(), ": results processed.  Now writing to files"

                do i=1,runs
                    write(unit = 841,FMT=*) all_results(i)%NC_abund_CEL
                    write(unit = 842,FMT=*) all_results(i)%Cii_abund_CEL
                    write(unit = 843,FMT=*) all_results(i)%Ciii_abund_CEL
                    write(unit = 844,FMT=*) all_results(i)%Civ_abund_CEL
                    write(unit = 845,FMT=*) all_results(i)%C_icf_CEL
                    write(unit = 846,FMT=*) all_results(i)%C_abund_CEL
                    write(unit = 847,FMT=*) all_results(i)%Nii_abund_CEL
                    write(unit = 848,FMT=*) all_results(i)%N_icf_CEL
                    write(unit = 849,FMT=*) all_results(i)%N_abund_CEL
                    write(unit = 850,FMT=*) all_results(i)%NO_abund_CEL
                    write(unit = 851,FMT=*) all_results(i)%Oii_abund_CEL
                    write(unit = 852,FMT=*) all_results(i)%Oiii_abund_CEL
                    write(unit = 853,FMT=*) all_results(i)%O_icf_CEL
                    write(unit = 854,FMT=*) all_results(i)%O_abund_CEL
                    write(unit = 855,FMT=*) all_results(i)%Neii_abund_CEL
                    write(unit = 856,FMT=*) all_results(i)%Neiii_abund_CEL
                    write(unit = 857,FMT=*) all_results(i)%Neiv_abund_CEL
                    write(unit = 858,FMT=*) all_results(i)%Nev_abund_CEL
                    write(unit = 859,FMT=*) all_results(i)%Ne_icf_CEL
                    write(unit = 860,FMT=*) all_results(i)%Ne_abund_CEL
                    write(unit = 861,FMT=*) all_results(i)%Ariii_abund_CEL
                    write(unit = 862,FMT=*) all_results(i)%Ariv_abund_CEL
                    write(unit = 863,FMT=*) all_results(i)%Arv_abund_CEL
                    write(unit = 864,FMT=*) all_results(i)%Ar_icf_CEL
                    write(unit = 865,FMT=*) all_results(i)%Ar_abund_CEL
                    write(unit = 866,FMT=*) all_results(i)%Sii_abund_CEL
                    write(unit = 867,FMT=*) all_results(i)%Siii_abund_CEL
                    write(unit = 868,FMT=*) all_results(i)%S_icf_CEL
                    write(unit = 869,FMT=*) all_results(i)%S_abund_CEL
                    write(unit = 870,FMT=*) all_results(i)%Cliii_abund_CEL
                    write(unit = 871,FMT=*) all_results(i)%Cl_icf_CEL
                    write(unit = 872,FMT=*) all_results(i)%Cl_abund_CEL
                    write(unit = 873,FMT=*) all_results(i)%Hei_abund_ORL
                    write(unit = 874,FMT=*) all_results(i)%Heii_abund_ORL
                    write(unit = 875,FMT=*) all_results(i)%He_abund_ORL
                    write(unit = 876,FMT=*) all_results(i)%Cii_abund_ORL
                    write(unit = 877,FMT=*) all_results(i)%Ciii_abund_ORL
                    write(unit = 878,FMT=*) all_results(i)%C_icf_ORL
                    write(unit = 879,FMT=*) all_results(i)%C_abund_ORL
                    write(unit = 880,FMT=*) all_results(i)%Nii_abund_ORL
                    write(unit = 881,FMT=*) all_results(i)%Niii_abund_ORL
                    write(unit = 882,FMT=*) all_results(i)%N_icf_ORL
                    write(unit = 883,FMT=*) all_results(i)%N_abund_ORL
                    write(unit = 884,FMT=*) all_results(i)%Oii_abund_ORL
                    write(unit = 885,FMT=*) all_results(i)%O_icf_ORL
                    write(unit = 886,FMT=*) all_results(i)%O_abund_ORL
                    write(unit = 887,FMT=*) all_results(i)%Neii_abund_ORL
                    write(unit = 888,FMT=*) all_results(i)%Ne_icf_ORL
                    write(unit = 889,FMT=*) all_results(i)%Ne_abund_ORL
                    write(unit = 890,FMT=*) all_results(i)%OII_density
                    write(unit = 891,FMT=*) all_results(i)%SII_density
                    write(unit = 892,FMT=*) all_results(i)%low_density
                    write(unit = 893,FMT=*) all_results(i)%OII_temp
                    write(unit = 894,FMT=*) all_results(i)%NII_temp
                    write(unit = 895,FMT=*) all_results(i)%SII_temp
                    write(unit = 896,FMT=*) all_results(i)%OI_temp
                    write(unit = 897,FMT=*) all_results(i)%CI_temp
                    write(unit = 898,FMT=*) all_results(i)%low_temp
                    write(unit = 899,FMT=*) all_results(i)%ClIII_density
                    write(unit = 900,FMT=*) all_results(i)%ArIV_density
                    write(unit = 901,FMT=*) all_results(i)%CIII_density
                    write(unit = 902,FMT=*) all_results(i)%ArIII_IR_density
                    write(unit = 903,FMT=*) all_results(i)%NeIII_IR_density
                    write(unit = 904,FMT=*) all_results(i)%med_density
                    write(unit = 905,FMT=*) all_results(i)%OIII_temp
                    write(unit = 906,FMT=*) all_results(i)%NeIII_temp
                    write(unit = 907,FMT=*) all_results(i)%ArIII_temp
                    write(unit = 908,FMT=*) all_results(i)%SIII_temp
                    write(unit = 909,FMT=*) all_results(i)%med_temp
                    write(unit = 910,FMT=*) all_results(i)%NeIV_density
                    write(unit = 911,FMT=*) all_results(i)%high_density
                    write(unit = 912,FMT=*) all_results(i)%ArV_temp
                    write(unit = 913,FMT=*) all_results(i)%NeV_temp
                    write(unit = 914,FMT=*) all_results(i)%high_temp
                    write(unit = 915,FMT=*) all_results(i)%mean_cHb
                    write(unit = 916,FMT=*) all_results(i)%adf_O
                    write(unit = 917,FMT=*) all_results(i)%adf_O2plus
                    write(unit = 918,FMT=*) all_results(i)%adf_C
                    write(unit = 919,FMT=*) all_results(i)%adf_C2plus
                    write(unit = 920,FMT=*) all_results(i)%adf_N
                    write(unit = 921,FMT=*) all_results(i)%adf_N2plus
                    write(unit = 922,FMT=*) all_results(i)%adf_Ne
                    write(unit = 923,FMT=*) all_results(i)%adf_Ne2plus
                    write(unit = 924,FMT=*) all_results(i)%O_R23_upper
                    write(unit = 925,FMT=*) all_results(i)%O_R23_lower
                    write(unit = 926,FMT=*) all_results(i)%O_N2
                    write(unit = 927,FMT=*) all_results(i)%O_O3N2
                    write(unit = 928,FMT=*) all_results(i)%O_Ar3O3
                    write(unit = 929,FMT=*) all_results(i)%O_S3O3
                    write(unit = 930,FMT=*) all_results(i)%OII_density_ratio
                    write(unit = 931,FMT=*) all_results(i)%SII_density_ratio
                    write(unit = 932,FMT=*) all_results(i)%NII_temp_ratio
                    write(unit = 933,FMT=*) all_results(i)%OII_temp_ratio
                    write(unit = 934,FMT=*) all_results(i)%SII_temp_ratio
                    write(unit = 935,FMT=*) all_results(i)%CI_temp_ratio
                    write(unit = 936,FMT=*) all_results(i)%OI_temp_ratio
                    write(unit = 937,FMT=*) all_results(i)%ClIII_density_ratio
                    write(unit = 938,FMT=*) all_results(i)%ArIV_density_ratio
                    write(unit = 939,FMT=*) all_results(i)%CIII_density_ratio
                    write(unit = 940,FMT=*) all_results(i)%OIII_IR_density_ratio
                    write(unit = 941,FMT=*) all_results(i)%ArIII_IR_density_ratio
                    write(unit = 942,FMT=*) all_results(i)%SIII_IR_density_ratio
                    write(unit = 943,FMT=*) all_results(i)%NeIII_IR_density_ratio
                    write(unit = 944,FMT=*) all_results(i)%OIII_temp_ratio
                    write(unit = 945,FMT=*) all_results(i)%NeIII_temp_ratio
                    write(unit = 946,FMT=*) all_results(i)%ArIII_temp_ratio
                    write(unit = 947,FMT=*) all_results(i)%SIII_temp_ratio
                    write(unit = 948,FMT=*) all_results(i)%OIII_IR_temp_ratio
                    write(unit = 949,FMT=*) all_results(i)%NeIII_IR_temp_ratio
                    write(unit = 950,FMT=*) all_results(i)%NeIV_density_ratio
                    write(unit = 951,FMT=*) all_results(i)%ArV_temp_ratio
                    write(unit = 952,FMT=*) all_results(i)%NeV_temp_ratio
                    write(unit = 953,FMT=*) all_results(i)%Bal_jump_temp
                end do

                DO I=841,953
                        CLOSE(unit=I)
                END DO
! get median +- pseudo gaussian 34.1%
                allocate (quantity_result(runs))

                open (650,FILE=trim(filename)//"_results", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')

                write (650,*) "NEAT (nebular empirical analysis tool)"
                write (650,*) "======================================"
                write (650,*)
                write (650,*) "Analysis of file ",trim(filename)
                write (650,*)

!cHb
                write (650,*) "Extinction"
                write (650,*) "=========="
                write (650,*)
                quantity_result = all_results%mean_cHb
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,"(X,A,F5.3,A,F5.3,A,F5.3)") "c(Hb) :                  ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

write (650,*)
write (650,*) "Diagnostics"
write (650,*) "==========="
write (650,*)
711 format (X,3(A,I8)) ! diagnostic format
715 format (X,3(A,G8.3)) ! diagnostic ratio format
!low densities

                quantity_result = all_results%oii_density
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,711) "[OII] density :         ",int(uncertainty_array(2))," +",int(uncertainty_array(3)),"-",int(uncertainty_array(1))
                quantity_result = all_results%oii_density_ratio
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,715) "Ratio :                  ",uncertainty_array(2),"+",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%SII_density
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,711) "[SII] density :         ",int(uncertainty_array(2))," +",int(uncertainty_array(3)),"-",int(uncertainty_array(1))
                quantity_result = all_results%sii_density_ratio
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,715) "Ratio :                  ",uncertainty_array(2),"+",uncertainty_array(3),"-",uncertainty_array(1)

write (650,*)

                quantity_result = all_results%low_density
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,711) "low density :           ",int(uncertainty_array(2))," +",int(uncertainty_array(3)),"-",int(uncertainty_array(1))

!low temperatures
write (650,*)

                quantity_result = all_results%oii_temp
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,711) "[OII] temperature :     ",int(uncertainty_array(2))," +",int(uncertainty_array(3)),"-",int(uncertainty_array(1))
                quantity_result = all_results%oii_temp_ratio
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,715) "Ratio :                  ",uncertainty_array(2),"+",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%SII_temp
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,711) "[SII] temperature :     ",int(uncertainty_array(2))," +",int(uncertainty_array(3)),"-",int(uncertainty_array(1))
                quantity_result = all_results%sii_temp_ratio
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,715) "Ratio :                  ",uncertainty_array(2),"+",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%NII_temp
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,711) "[NII] temperature :     ",int(uncertainty_array(2))," +",int(uncertainty_array(3)),"-",int(uncertainty_array(1))
                quantity_result = all_results%nii_temp_ratio
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,715) "Ratio :                  ",uncertainty_array(2),"+",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%OI_temp
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,711) "[OI] temperature :      ",int(uncertainty_array(2))," +",int(uncertainty_array(3)),"-",int(uncertainty_array(1))
                quantity_result = all_results%oi_temp_ratio
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,715) "Ratio :                  ",uncertainty_array(2),"+",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%CI_temp
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,711) "[CI] temperature :      ",int(uncertainty_array(2))," +",int(uncertainty_array(3)),"-",int(uncertainty_array(1))
                quantity_result = all_results%ci_temp_ratio
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,715) "Ratio :                  ",uncertainty_array(2),"+",uncertainty_array(3),"-",uncertainty_array(1)

write (650,*)

                quantity_result = all_results%low_temp
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,711) "low temperature :       ",int(uncertainty_array(2))," +",int(uncertainty_array(3)),"-",int(uncertainty_array(1))

!medium density
write (650,*)

                quantity_result = all_results%cliii_density
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,711) "[ClIII] density :       ",int(uncertainty_array(2))," +",int(uncertainty_array(3)),"-",int(uncertainty_array(1))
                quantity_result = all_results%cliii_density_ratio
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,715) "Ratio :                  ",uncertainty_array(2),"+",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%ArIV_density
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,711) "[ArIV] density :        ",int(uncertainty_array(2))," +",int(uncertainty_array(3)),"-",int(uncertainty_array(1))
                quantity_result = all_results%Ariv_density_ratio
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,715) "Ratio :                  ",uncertainty_array(2),"+",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%CIII_density
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,711) "[CIII] density :        ",int(uncertainty_array(2))," +",int(uncertainty_array(3)),"-",int(uncertainty_array(1))
                quantity_result = all_results%ciii_density_ratio
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,715) "Ratio :                  ",uncertainty_array(2),"+",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%OIII_IR_density
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,711) "[OIII] IR density :     ",int(uncertainty_array(2))," +",int(uncertainty_array(3)),"-",int(uncertainty_array(1))
                quantity_result = all_results%oiii_ir_density_ratio
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,715) "Ratio :                  ",uncertainty_array(2),"+",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%SIII_IR_density
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,711) "[SIII] IRdensity :      ",int(uncertainty_array(2))," +",int(uncertainty_array(3)),"-",int(uncertainty_array(1))
                quantity_result = all_results%siii_ir_density_ratio
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,715) "Ratio :                  ",uncertainty_array(2),"+",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%ArIII_IR_density
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,711) "[ArIII] IR density :    ",int(uncertainty_array(2))," +",int(uncertainty_array(3)),"-",int(uncertainty_array(1))
                quantity_result = all_results%ariii_ir_density_ratio
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,715) "Ratio :                  ",uncertainty_array(2),"+",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%NeIII_IR_density
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,711) "[NeIII] IR density :    ",int(uncertainty_array(2))," +",int(uncertainty_array(3)),"-",int(uncertainty_array(1))
                quantity_result = all_results%neiii_ir_density_ratio
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,715) "Ratio :                  ",uncertainty_array(2),"+",uncertainty_array(3),"-",uncertainty_array(1)

write (650,*)

                quantity_result = all_results%med_density
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,711) "medium density :        ",int(uncertainty_array(2))," +",int(uncertainty_array(3)),"-",int(uncertainty_array(1))

!medium temperature
write (650,*)

                quantity_result = all_results%OIII_temp
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,711) "[OIII] temperature :    ",int(uncertainty_array(2))," +",int(uncertainty_array(3)),"-",int(uncertainty_array(1))
                quantity_result = all_results%oiii_temp_ratio
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,715) "Ratio :                  ",uncertainty_array(2),"+",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%OIII_IR_temp
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,711) "[OIII] IR temperature : ",int(uncertainty_array(2))," +",int(uncertainty_array(3)),"-",int(uncertainty_array(1))
                quantity_result = all_results%oiii_ir_temp_ratio
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,715) "Ratio :                  ",uncertainty_array(2),"+",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%NeIII_temp
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,711) "[NeIII] temperature :   ",int(uncertainty_array(2))," +",int(uncertainty_array(3)),"-",int(uncertainty_array(1))
                quantity_result = all_results%neiii_temp_ratio
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,715) "Ratio :                  ",uncertainty_array(2),"+",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%NeIII_IR_temp
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,711) "[NeIII] IR temperature :",int(uncertainty_array(2))," +",int(uncertainty_array(3)),"-",int(uncertainty_array(1))
                quantity_result = all_results%neiii_ir_temp_ratio
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,715) "Ratio :                  ",uncertainty_array(2),"+",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%ArIII_temp
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,711) "[ArIII] temperature :   ",int(uncertainty_array(2))," +",int(uncertainty_array(3)),"-",int(uncertainty_array(1))
                quantity_result = all_results%ariii_temp_ratio
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,715) "Ratio :                  ",uncertainty_array(2),"+",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%SIII_temp
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,711) "[SIII] temperature :    ",int(uncertainty_array(2))," +",int(uncertainty_array(3)),"-",int(uncertainty_array(1))
                quantity_result = all_results%siii_temp_ratio
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,715) "Ratio :                  ",uncertainty_array(2),"+",uncertainty_array(3),"-",uncertainty_array(1)

write (650,*)

                quantity_result = all_results%med_temp
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,711) "medium temperature :    ",int(uncertainty_array(2))," +",int(uncertainty_array(3)),"-",int(uncertainty_array(1))

!high density
write (650,*)

                quantity_result = all_results%neiv_density
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,711) "[NeIV] density :        ",int(uncertainty_array(2))," +",int(uncertainty_array(3)),"-",int(uncertainty_array(1))
                quantity_result = all_results%neiv_density_ratio
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,715) "Ratio :                  ",uncertainty_array(2),"+",uncertainty_array(3),"-",uncertainty_array(1)

write (650,*)

                quantity_result = all_results%high_density
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,711) "high density :          ",int(uncertainty_array(2))," +",int(uncertainty_array(3)),"-",int(uncertainty_array(1))

!high temperature

                quantity_result = all_results%ArV_temp
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,711) "[ArV] temperature :     ",int(uncertainty_array(2))," +",int(uncertainty_array(3)),"-",int(uncertainty_array(1))
                quantity_result = all_results%arv_temp_ratio
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,715) "Ratio :                  ",uncertainty_array(2),"+",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%NeV_temp
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,711) "[NeV] temperature :     ",int(uncertainty_array(2))," +",int(uncertainty_array(3)),"-",int(uncertainty_array(1))
                quantity_result = all_results%nev_temp_ratio
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,715) "Ratio :                  ",uncertainty_array(2),"+",uncertainty_array(3),"-",uncertainty_array(1)

write (650,*)

                quantity_result = all_results%high_temp
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,711) "High temperature :      ",int(uncertainty_array(2))," +",int(uncertainty_array(3)),"-",int(uncertainty_array(1))
                
                quantity_result = all_results%Bal_jump_temp
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,711) "Balmer jump temperature :      ",int(uncertainty_array(2))," +",int(uncertainty_array(3)),"-",int(uncertainty_array(1))

!CEL abundances
write (650,*)
write (650,*) "Abundances (collisionally excited lines)"
write (650,*) "========================================"
write (650,*)

713 format (X,3(A,ES10.2)) ! abundances formats


                quantity_result = all_results%NC_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "[CI] abundance :        ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%cii_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "[CII] abundance :       ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%ciii_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "[CIII] abundance :      ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%civ_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "[CIV] abundance :       ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%C_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "C/H abundance :         ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%nii_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "[NII] abundance :       ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%niii_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "[NIII] abundance :      ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%niv_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "[NIV] abundance :       ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%nv_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "[NV] abundance :        ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%N_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "N/H abundance :         ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%NO_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "[OI] abundance :        ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%Oii_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "[OII] abundance :       ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%Oiii_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "[OIII] abundance :      ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%Oiv_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "[OIV] abundance :       ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%O_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "O/H abundance :         ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

!                quantity_result = all_results%NeII_abund_CEL
!                call get_uncertainties(quantity_result, uncertainty_array)
!                write (650,713) "[NeII] abundance :  ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)
!
                quantity_result = all_results%NeIII_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "[NeIII] abundance :     ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%NeIV_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "[NeIV] abundance :      ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%NeV_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "[NeV] abundance :       ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%Ne_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "Ne/H abundance :        ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%ArIII_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "[ArIII] abundance :     ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%ArIV_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "[ArIV] abundance :      ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%ArV_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "[ArV] abundance :       ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%Ar_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "Ar/H abundance :        ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%SII_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "[SII] abundance :       ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%SIII_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "[SIII] abundance :      ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%S_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "S/H abundance :         ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%ClIII_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "[ClIII] abundance :     ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%Cl_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "Cl/H abundance :        ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

!RL abundances
write (650,*)
write (650,*) "Abundances (recombination lines)"
write (650,*) "================================"
write (650,*)

                quantity_result = all_results%He_abund_ORL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "He/H abundance :        ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%C_abund_ORL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "C/H abundance :         ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%N_abund_ORL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "N/H abundance :         ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%O_abund_ORL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "O/H abundance :         ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%Ne_abund_ORL
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "Ne/H abundance :        ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

!strong line abundances
write (650,*)
write (650,*) "Abundances (strong line methods)"
write (650,*) "================================"
write (650,*)

                quantity_result = all_results%O_R23_upper
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "O/H (R23 upper) :       ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%O_R23_lower
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "O/H (R23 lower) :       ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%O_N2
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "O/H (N2) :              ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%O_O3N2
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "O/H (O3N2) :            ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%O_Ar3O3
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "O/H (Ar3O3) :           ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%O_S3O3
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,713) "O/H (S3O3) :            ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

!adfs
write (650,*)
write (650,*) "Abundance discrepancy factors"
write (650,*) "============================="
write (650,*)

714 format (X,3(A,F5.2))


                quantity_result = all_results%adf_o2plus
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,714) "adf(O2+/H) :            ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%adf_o
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,714) "adf(O/H+) :             ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%adf_n2plus
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,714) "adf(N2+/H) :            ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%adf_n
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,714) "adf(N/H) :              ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%adf_c2plus
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,714) "adf(C2+/H+) :           ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%adf_c
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,714) "adf(C/H) :              ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%adf_ne2plus
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,714) "adf(Ne2+/H+) :          ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                quantity_result = all_results%adf_ne
                call get_uncertainties(quantity_result, uncertainty_array)
                write (650,714) "adf(Ne/H) :             ",uncertainty_array(2)," +",uncertainty_array(3),"-",uncertainty_array(1)

                close (650)

        else
                print*, gettime(), ": I didn't want to be a barber anyway. I wanted to be... a lumberjack!   Also, a positive number of runs helps.."
        endif

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

subroutine get_uncertainties(input_array, uncertainty_array)

double precision, intent(in) :: input_array(:)
double precision, intent(out) :: uncertainty_array(3)
double precision :: binsize, comp
double precision, dimension (:,:), allocatable :: binned_quantity_result
integer :: ii, bincount, bincountmax, arraysize, abovepos, belowpos

uncertainty_array = (/0.0,0.0,0.0/)

arraysize = size(input_array)
binsize=(input_array(int(0.841*size(input_array))) - input_array(int(0.159*size(input_array))))/20

if (binsize .gt. 0) then
  allocate(binned_quantity_result(arraysize,2))
  binned_quantity_result = 0.D0

  ii=1
  bincount=1 !(why does this need to be one and not zero??)
  binvalue = int(quantity_result(1)/binsize)
  bincountmax=0

  do i=1,runs
    if (int(quantity_result(i)/binsize) == binvalue) then
      bincount = bincount + 1
    else
      binned_quantity_result(ii,1) = binvalue*binsize
      binned_quantity_result(ii,2) = bincount

      if (bincount>bincountmax) then
        bincountmax=bincount
        uncertainty_array(2) = binned_quantity_result(ii,1) + 0.5*binsize ! otherwise it is the value of the edge of the bin and not the centre
      endif

      ii=ii+1
      bincount = 1
      binvalue = binvalue + 1
    endif
  enddo

  deallocate(binned_quantity_result)

!find value in array closest to mode

  comp = 1.e10
  do i=1,arraysize
    if (abs(input_array(i)-uncertainty_array(2))>comp) exit
    comp = abs(input_array(i)-uncertainty_array(2))
  enddo

  abovepos = i+int(0.341*arraysize)
  belowpos = i-int(0.341*arraysize)

  if (abovepos>arraysize) then
    uncertainty_array(3) = 99999999
  else
    uncertainty_array(3) = input_array(abovepos) - uncertainty_array(2)
  endif

  if (belowpos<1) then
    uncertainty_array(1) = 0.D0
  else
    uncertainty_array(1) = uncertainty_array(2) - input_array(belowpos)
  endif

else !all results are identical
  uncertainty_array(1) = 0.D0
  uncertainty_array(2) = input_array(1)
  uncertainty_array(3) = 0.D0
endif

end subroutine get_uncertainties

character*10 function gettime()

character*10 :: time

  call DATE_AND_TIME(TIME=time)
  gettime = time(1:2)//":"//time(3:4)//":"//time(5:6)
  return

end function gettime

end program
