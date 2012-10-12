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
        CHARACTeR :: switch_he  !switch for helium atomic data
        INTEGER :: I, runs, Narg !runs = number of runs for randomiser

        !input options

        CHARACTER*2048, DIMENSION(:), allocatable :: options
        CHARACTER*2048 :: commandline

        !file reading variables

        LOGICAL :: file_exists
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
        double precision, dimension(21,15,44) :: heidata

        !extinction

        logical :: calculate_extinction=.true.
        DOUBLE PRECISION :: temp1,temp2,temp3, meanextinction, R

        !binning

        double precision, dimension(3) :: uncertainty_array

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
           print *,"  -he / --helium-data"
           print *,"       The atomic data to use for He I abundances"
           print *,"       Default: Smits, 1996, MNRAS, 278, 683"
           print *,"       Values:"
           print *,"          S96: Smits, 1996, MNRAS, 278, 683"
           print *,"          P12: Porter et al., 2012, MNRAS, 425, 28"
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
        switch_he="S"
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
                if ((trim(options(i))=="-he" .or. trim(options(i))=="--helium-data") .and. (i+1) .le. Narg) then
                  if (trim(options(i+1))=="P12") then
                    switch_he="P"
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

        !read Porter et al helium emissivities
        call read_porter(heidata)

        !find maximum #levels and temperatures - pass to equib to reduce footprint

        maxlevs = atomicdata(1)%nlevs
        maxtemps = atomicdata(1)%ntemps
        do i=2,iion
            if(atomicdata(i)%nlevs .gt. maxlevs) maxlevs = atomicdata(i)%nlevs
            if(atomicdata(i)%ntemps .gt. maxtemps) maxtemps = atomicdata(i)%ntemps
        enddo
!        print*,maxlevs,maxtemps
        !now check number of iterations.  If 1, line list is fine as is.  If more than one, randomize the fluxes


!formats for writing output
!single iteration formats

        if(runs == 1)then !calculates abundances without uncertainties
                call abundances(linelist, switch_ext, listlength, iteration_result, R, meanextinction, calculate_extinction, ILs, Iint, diagnostic_array,iion,atomicdata,maxlevs,maxtemps, heidata, switch_he)

                !generate outputs

                print *,gettime(),": writing results to file ",trim(filename),"_results"
                OPEN(650, FILE=trim(filename)//"_results", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')

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
                write (650,*)
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
                        call abundances(linelist, switch_ext, listlength, iteration_result, R, meanextinction, calculate_extinction, ILs, Iint, diagnostic_array,iion,atomicdata,maxlevs,maxtemps, heidata, switch_he)
                        linelist = linelist_original
                        all_results(i)=iteration_result(1)
                END DO

                ! now process outputs
                print *, gettime(), ": processing results"

!XXXX add Cl/H, Niii, cii, ciii, ArIII IR dens, NeIII IR dens, strong line, ICF files to output

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
                call get_uncertainties(quantity_result, uncertainty_array,"c(Hb) :                 ", extinction_format, filename, "mean_chb                 ")

write (650,*)
write (650,*) "Diagnostics"
write (650,*) "==========="
write (650,*)

!low densities

                quantity_result = all_results%oii_density
                call get_uncertainties(quantity_result, uncertainty_array, "[OII] density :         ", diagnostic_format, filename, "oii_density              ")
                quantity_result = all_results%oii_density_ratio
                call get_uncertainties(quantity_result, uncertainty_array, "Ratio :                 ", diagnostic_ratio_format, filename, "oii_density_ratio        ")

                quantity_result = all_results%SII_density
                call get_uncertainties(quantity_result, uncertainty_array, "[SII] density :         ", diagnostic_format, filename, "sii_density              ")
                quantity_result = all_results%sii_density_ratio
                call get_uncertainties(quantity_result, uncertainty_array, "Ratio :                 ", diagnostic_ratio_format, filename, "sii_density_ratio        ")

write (650,*)

                quantity_result = all_results%low_density
                call get_uncertainties(quantity_result, uncertainty_array, "low density :           ", diagnostic_format, filename, "low_density              ")

!low temperatures
write (650,*)

                quantity_result = all_results%oii_temp
                call get_uncertainties(quantity_result, uncertainty_array, "[OII] temperature :     ", diagnostic_format, filename, "oii_temp                 ")
                quantity_result = all_results%oii_temp_ratio
                call get_uncertainties(quantity_result, uncertainty_array, "Ratio :                 ", diagnostic_ratio_format, filename, "oii_temp_ratio           ")

                quantity_result = all_results%SII_temp
                call get_uncertainties(quantity_result, uncertainty_array, "[SII] temperature :     ", diagnostic_format, filename, "sii_temp                 ")
                quantity_result = all_results%sii_temp_ratio
                call get_uncertainties(quantity_result, uncertainty_array, "Ratio :                 ", diagnostic_ratio_format, filename, "sii_temp_ratio           ")

                quantity_result = all_results%NII_temp
                call get_uncertainties(quantity_result, uncertainty_array, "[NII] temperature :     ", diagnostic_format, filename, "nii_temp                 ")
                quantity_result = all_results%nii_temp_ratio
                call get_uncertainties(quantity_result, uncertainty_array, "Ratio :                 ", diagnostic_ratio_format, filename, "nii_temp_ratio           ")

                quantity_result = all_results%OI_temp
                call get_uncertainties(quantity_result, uncertainty_array, "[OI] temperature :      ", diagnostic_format, filename, "oi_temp                  ")
                quantity_result = all_results%oi_temp_ratio
                call get_uncertainties(quantity_result, uncertainty_array, "Ratio :                 ", diagnostic_ratio_format, filename, "oi_temp_ratio            ")

                quantity_result = all_results%CI_temp
                call get_uncertainties(quantity_result, uncertainty_array, "[CI] temperature :      ", diagnostic_format, filename, "ci_temp                  ")
                quantity_result = all_results%ci_temp_ratio
                call get_uncertainties(quantity_result, uncertainty_array, "Ratio :                 ", diagnostic_ratio_format, filename, "ci_temp_ratio            ")

write (650,*)

                quantity_result = all_results%low_temp
                call get_uncertainties(quantity_result, uncertainty_array, "low temperature :       ", diagnostic_format, filename, "low_temp                 ")

!medium density
write (650,*)

                quantity_result = all_results%cliii_density
                call get_uncertainties(quantity_result, uncertainty_array, "[ClIII] density :       ", diagnostic_format, filename, "cliii_density            ")
                quantity_result = all_results%cliii_density_ratio
                call get_uncertainties(quantity_result, uncertainty_array, "Ratio :                 ", diagnostic_ratio_format, filename, "cliii_density_ratio      ")

                quantity_result = all_results%ArIV_density
                call get_uncertainties(quantity_result, uncertainty_array, "[ArIV] density :        ", diagnostic_format, filename, "ariv_density             ")
                quantity_result = all_results%Ariv_density_ratio
                call get_uncertainties(quantity_result, uncertainty_array, "Ratio :                 ", diagnostic_ratio_format, filename, "ariv_density_ratio       ")

                quantity_result = all_results%CIII_density
                call get_uncertainties(quantity_result, uncertainty_array, "[CIII] density :        ", diagnostic_format, filename, "ciii_density             ")
                quantity_result = all_results%ciii_density_ratio
                call get_uncertainties(quantity_result, uncertainty_array, "Ratio :                 ", diagnostic_ratio_format, filename, "ciii_density_ratio       ")

                quantity_result = all_results%OIII_IR_density
                call get_uncertainties(quantity_result, uncertainty_array, "[OIII] IR density :     ", diagnostic_format, filename, "oiii_ir_density          ")
                quantity_result = all_results%oiii_ir_density_ratio
                call get_uncertainties(quantity_result, uncertainty_array, "Ratio :                 ", diagnostic_ratio_format, filename, "oiii_ir_density_ratio    ")

                quantity_result = all_results%SIII_IR_density
                call get_uncertainties(quantity_result, uncertainty_array, "[SIII] IR density :     ", diagnostic_format, filename, "siii_ir_density          ")
                quantity_result = all_results%siii_ir_density_ratio
                call get_uncertainties(quantity_result, uncertainty_array, "Ratio :                 ", diagnostic_ratio_format, filename, "siii_ir_density_ratio    ")

                quantity_result = all_results%ArIII_IR_density
                call get_uncertainties(quantity_result, uncertainty_array, "[ArIII] IR density :    ", diagnostic_format, filename, "ariii_ir_density         ")
                quantity_result = all_results%ariii_ir_density_ratio
                call get_uncertainties(quantity_result, uncertainty_array, "Ratio :                 ", diagnostic_ratio_format, filename, "ariii_ir_density_ratio   ")

                quantity_result = all_results%NeIII_IR_density
                call get_uncertainties(quantity_result, uncertainty_array, "[NeIII] IR density :    ", diagnostic_format, filename, "neiii_ir_density         ")
                quantity_result = all_results%neiii_ir_density_ratio
                call get_uncertainties(quantity_result, uncertainty_array, "Ratio :                 ", diagnostic_ratio_format, filename, "neiii_ir_density_ratio   ")

write (650,*)

                quantity_result = all_results%med_density
                call get_uncertainties(quantity_result, uncertainty_array, "medium density :        ", diagnostic_format, filename, "med_density              ")

!medium temperature
write (650,*)

                quantity_result = all_results%OIII_temp
                call get_uncertainties(quantity_result, uncertainty_array, "[OIII] temperature :    ", diagnostic_format, filename, "oiii_temp                ")
                quantity_result = all_results%oiii_temp_ratio
                call get_uncertainties(quantity_result, uncertainty_array, "Ratio :                 ", diagnostic_ratio_format, filename, "oiii_temp_ratio          ")

                quantity_result = all_results%OIII_IR_temp
                call get_uncertainties(quantity_result, uncertainty_array, "[OIII] IR temperature : ", diagnostic_format, filename, "oiii_ir_temp             ")
                quantity_result = all_results%oiii_ir_temp_ratio
                call get_uncertainties(quantity_result, uncertainty_array, "Ratio :                 ", diagnostic_ratio_format, filename, "oiii_ir_temp_ratio       ")

                quantity_result = all_results%NeIII_temp
                call get_uncertainties(quantity_result, uncertainty_array, "[NeIII] temperature :   ", diagnostic_format, filename, "neiii_temp               ")
                quantity_result = all_results%neiii_temp_ratio
                call get_uncertainties(quantity_result, uncertainty_array, "Ratio :                 ", diagnostic_ratio_format, filename, "neiii_temp_ratio         ")

                quantity_result = all_results%NeIII_IR_temp
                call get_uncertainties(quantity_result, uncertainty_array, "[NeIII] IR temperature :", diagnostic_format, filename, "neiii_ir_temp            ")
                quantity_result = all_results%neiii_ir_temp_ratio
                call get_uncertainties(quantity_result, uncertainty_array, "Ratio :                 ", diagnostic_ratio_format, filename, "neiii_ir_temp_ratio      ")

                quantity_result = all_results%ArIII_temp
                call get_uncertainties(quantity_result, uncertainty_array, "[ArIII] temperature :   ", diagnostic_format, filename, "ariii_temp               ")
                quantity_result = all_results%ariii_temp_ratio
                call get_uncertainties(quantity_result, uncertainty_array, "Ratio :                 ", diagnostic_ratio_format, filename, "ariii_temp_ratio         ")

                quantity_result = all_results%SIII_temp
                call get_uncertainties(quantity_result, uncertainty_array, "[SIII] temperature :    ", diagnostic_format, filename, "siii_temp                ")
                quantity_result = all_results%siii_temp_ratio
                call get_uncertainties(quantity_result, uncertainty_array, "Ratio :                 ", diagnostic_ratio_format, filename, "siii_temp_ratio          ")

write (650,*)

                quantity_result = all_results%med_temp
                call get_uncertainties(quantity_result, uncertainty_array, "medium temperature :    ", diagnostic_format, filename, "med_temp                 ")

!high density
write (650,*)

                quantity_result = all_results%neiv_density
                call get_uncertainties(quantity_result, uncertainty_array, "[NeIV] density :        ", diagnostic_format, filename, "neiv_density             ")
                quantity_result = all_results%neiv_density_ratio
                call get_uncertainties(quantity_result, uncertainty_array, "Ratio :                 ", diagnostic_ratio_format, filename, "neiv_density_ratio       ")

write (650,*)

                quantity_result = all_results%high_density
                call get_uncertainties(quantity_result, uncertainty_array, "high density :          ", diagnostic_format, filename, "high_density             ")

!high temperature

                quantity_result = all_results%ArV_temp
                call get_uncertainties(quantity_result, uncertainty_array, "[ArV] temperature :     ", diagnostic_format, filename, "arv_temp                 ")
                quantity_result = all_results%arv_temp_ratio
                call get_uncertainties(quantity_result, uncertainty_array, "Ratio :                 ", diagnostic_ratio_format, filename, "arv_temp_ratio           ")

                quantity_result = all_results%NeV_temp
                call get_uncertainties(quantity_result, uncertainty_array, "[NeV] temperature :     ", diagnostic_format, filename, "nev_temp                 ")
                quantity_result = all_results%nev_temp_ratio
                call get_uncertainties(quantity_result, uncertainty_array, "Ratio :                 ", diagnostic_ratio_format, filename, "nev_temp_ratio           ")

write (650,*)

                quantity_result = all_results%high_temp
                call get_uncertainties(quantity_result, uncertainty_array, "High temperature :      ", diagnostic_format, filename, "high_temp                ")

write (650,*)
                
                quantity_result = all_results%Bal_jump_temp
                call get_uncertainties(quantity_result, uncertainty_array, "Balmer jump temp :      ", diagnostic_format, filename, "bal_jump_temp            ")

!CEL abundances
write (650,*)
write (650,*) "Abundances (collisionally excited lines)"
write (650,*) "========================================"
write (650,*)

                quantity_result = all_results%NC_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array, "[CI] abundance :        ", abundances_format, filename, "nc_abund_cel             ")

                quantity_result = all_results%cii_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array, "[CII] abundance :       ", abundances_format, filename, "cii_abund_cel            ")

                quantity_result = all_results%ciii_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array, "[CIII] abundance :      ", abundances_format, filename, "ciii_abund_cel           ")

                quantity_result = all_results%civ_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array, "[CIV] abundance :       ", abundances_format, filename, "civ_abund_cel            ")

                quantity_result = all_results%C_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array, "C/H abundance :         ", abundances_format, filename, "c_abund_cel              ")

                quantity_result = all_results%nii_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array, "[NII] abundance :       ", abundances_format, filename, "nii_abund_cel            ")

                quantity_result = all_results%niii_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array, "[NIII] abundance :      ", abundances_format, filename, "niii_abund_cel           ")

                quantity_result = all_results%niv_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array, "[NIV] abundance :       ", abundances_format, filename, "niv_abund_cel            ")

                quantity_result = all_results%nv_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array, "[NV] abundance :        ", abundances_format, filename, "nv_abund_cel             ")

                quantity_result = all_results%N_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array, "N/H abundance :         ", abundances_format, filename, "n_abund_cel              ")

                quantity_result = all_results%NO_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array, "[OI] abundance :        ", abundances_format, filename, "no_abund_cel             ")

                quantity_result = all_results%Oii_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array, "[OII] abundance :       ", abundances_format, filename, "oii_abund_cel            ")

                quantity_result = all_results%Oiii_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array, "[OIII] abundance :      ", abundances_format, filename, "oiii_abund_cel           ")

                quantity_result = all_results%Oiv_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array, "[OIV] abundance :       ", abundances_format, filename, "oiv_abund_cel            ")

                quantity_result = all_results%O_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array, "O/H abundance :         ", abundances_format, filename, "o_abund_cel              ")

!                quantity_result = all_results%NeII_abund_CEL
!                call get_uncertainties(quantity_result, uncertainty_array, filename, "neii_abund_cel           ")
!                write (650,713) "[NeII] abundance :  ",uncertainty_array(2),uncertainty_array(3),-uncertainty_array(1)
!
                quantity_result = all_results%NeIII_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array, "[NeIII] abundance :     ", abundances_format, filename, "neiii_abund_cel          ")

                quantity_result = all_results%NeIV_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array, "[NeIV] abundance :      ", abundances_format, filename, "neiv_abund_cel           ")

                quantity_result = all_results%NeV_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array, "[NeV] abundance :       ", abundances_format, filename, "nev_abund_cel            ")

                quantity_result = all_results%Ne_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array, "Ne/H abundance :        ", abundances_format, filename, "ne_abund_cel             ")

                quantity_result = all_results%ArIII_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array, "[ArIII] abundance :     ", abundances_format, filename, "ariii_abund_cel          ")

                quantity_result = all_results%ArIV_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array, "[ArIV] abundance :      ", abundances_format, filename, "ariv_abund_cel           ")

                quantity_result = all_results%ArV_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array, "[ArV] abundance :       ", abundances_format, filename, "arv_abund_cel            ")

                quantity_result = all_results%Ar_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array, "Ar/H abundance :        ", abundances_format, filename, "ar_abund_cel             ")

                quantity_result = all_results%SII_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array, "[SII] abundance :       ", abundances_format, filename, "sii_abund_cel            ")

                quantity_result = all_results%SIII_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array, "[SIII] abundance :      ", abundances_format, filename, "siii_abund_cel           ")

                quantity_result = all_results%S_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array, "S/H abundance :         ", abundances_format, filename, "s_abund_cel              ")

                quantity_result = all_results%ClIII_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array, "[ClIII] abundance :     ", abundances_format, filename, "cliii_abund_cel          ")

                quantity_result = all_results%Cl_abund_CEL
                call get_uncertainties(quantity_result, uncertainty_array, "Cl/H abundance :        ", abundances_format, filename, "cl_abund_cel             ")

!RL abundances
write (650,*)
write (650,*) "Abundances (recombination lines)"
write (650,*) "================================"
write (650,*)

                quantity_result = all_results%He_abund_ORL
                call get_uncertainties(quantity_result, uncertainty_array, "He/H abundance :        ", abundances_format, filename, "he_abund_orl             ")

                quantity_result = all_results%C_abund_ORL
                call get_uncertainties(quantity_result, uncertainty_array, "C/H abundance :         ", abundances_format, filename, "c_abund_orl              ")

                quantity_result = all_results%N_abund_ORL
                call get_uncertainties(quantity_result, uncertainty_array, "N/H abundance :         ", abundances_format, filename, "n_abund_orl              ")

                quantity_result = all_results%O_abund_ORL
                call get_uncertainties(quantity_result, uncertainty_array, "O/H abundance :         ", abundances_format, filename, "o_abund_orl              ")

                quantity_result = all_results%Ne_abund_ORL
                call get_uncertainties(quantity_result, uncertainty_array, "Ne/H abundance :        ", abundances_format, filename, "ne_abund_orl             ")

!strong line abundances
write (650,*)
write (650,*) "Abundances (strong line methods)"
write (650,*) "================================"
write (650,*)

                quantity_result = all_results%O_R23_upper
                call get_uncertainties(quantity_result, uncertainty_array, "O/H (R23 upper) :       ", abundances_format, filename, "o_r23_upper              ")

                quantity_result = all_results%O_R23_lower
                call get_uncertainties(quantity_result, uncertainty_array, "O/H (R23 lower) :       ", abundances_format, filename, "o_r23_lower              ")

                quantity_result = all_results%O_N2
                call get_uncertainties(quantity_result, uncertainty_array, "O/H (N2) :              ", abundances_format, filename, "o_n2                     ")

                quantity_result = all_results%O_O3N2
                call get_uncertainties(quantity_result, uncertainty_array, "O/H (O3N2) :            ", abundances_format, filename, "o_o3n2                   ")

                quantity_result = all_results%O_Ar3O3
                call get_uncertainties(quantity_result, uncertainty_array, "O/H (Ar3O3) :           ", abundances_format, filename, "o_ar3o3                  ")

                quantity_result = all_results%O_S3O3
                call get_uncertainties(quantity_result, uncertainty_array, "O/H (S3O3) :            ", abundances_format, filename, "o_s3o3                   ")

!adfs
write (650,*)
write (650,*) "Abundance discrepancy factors"
write (650,*) "============================="
write (650,*)

                quantity_result = all_results%adf_o2plus
                call get_uncertainties(quantity_result, uncertainty_array, "adf(O2+/H) :            ", adf_format, filename, "adf_o2plus               ")

                quantity_result = all_results%adf_o
                call get_uncertainties(quantity_result, uncertainty_array, "adf(O/H+) :             ", adf_format, filename, "adf_o                    ")

                quantity_result = all_results%adf_n2plus
                call get_uncertainties(quantity_result, uncertainty_array, "adf(N2+/H) :            ", adf_format, filename, "adf_n2plus               ")

                quantity_result = all_results%adf_n
                call get_uncertainties(quantity_result, uncertainty_array, "adf(N/H) :              ", adf_format, filename, "adf_n                    ")

                quantity_result = all_results%adf_c2plus
                call get_uncertainties(quantity_result, uncertainty_array, "adf(C2+/H+) :           ", adf_format, filename, "adf_c2plus               ")

                quantity_result = all_results%adf_c
                call get_uncertainties(quantity_result, uncertainty_array, "adf(C/H) :              ", adf_format, filename, "adf_c                    ")

                quantity_result = all_results%adf_ne2plus
                call get_uncertainties(quantity_result, uncertainty_array, "adf(Ne2+/H+) :          ", adf_format, filename, "adf_ne2plus              ")

                quantity_result = all_results%adf_ne
                call get_uncertainties(quantity_result, uncertainty_array, "adf(Ne/H) :             ", adf_format, filename, "adf_ne                   ")

                close (650)

                print *, gettime(), ": results summary file ",trim(filename)//"_results written"

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

subroutine get_uncertainties(input_array, uncertainty_array, itemtext,itemformat, filename, suffix)

implicit none
double precision :: input_array(:)
double precision, intent(out) :: uncertainty_array(3)
double precision, dimension(:), allocatable :: bintemp
double precision :: binsize, comp
double precision, dimension (:,:), allocatable :: binned_quantity_result
integer :: ii, bincount, bincountmax, arraysize, abovepos, belowpos, nbins
integer, dimension(1) :: maxpos
character*24, intent(in) :: itemtext
character*35, intent(in) :: itemformat
character*80, intent(in) :: filename
character*25, intent(in) :: suffix

!binned_quantity_result = 0.D0
uncertainty_array = (/0.0,0.0,0.0/)
arraysize = size(input_array)

!todo: include a verbosity switch
!verbosity: 0=summary file only, 1=summary file + binned results, 2=summary, binned and full
!set verbosity to 0 automatically for single iteration, then we can use the same routine for any n

!sort the input array into ascending order

call qsort(input_array)

!write this to file

OPEN(850, FILE=trim(filename)//"_"//trim(suffix), STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
do i=1,arraysize
  write(unit = 850,FMT=*) input_array(i)
end do
close(850)

! bin the array

arraysize = size(input_array)
binsize=(input_array(nint(0.841*size(input_array))) - input_array(nint(0.159*size(input_array))))/20

if (binsize .gt. 0) then

!quantize the input array by taking the integer of each value divided by the binsize, and multiplying by the binsize.
allocate(bintemp(arraysize))
bintemp = binsize*nint(input_array/binsize)
nbins=nint((maxval(bintemp)-minval(bintemp))/binsize)+1
allocate(binned_quantity_result(nbins,2))

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

  comp = 1.e10
  do i=1,arraysize
    if (abs(input_array(i)-uncertainty_array(2))>comp) exit
    comp = abs(input_array(i)-uncertainty_array(2))
  enddo

  abovepos = i+nint(0.341*arraysize)
  belowpos = i-nint(0.341*arraysize)

  if (abovepos>arraysize) then
    uncertainty_array(3) = 99999999
  else
    uncertainty_array(3) = input_array(abovepos) - uncertainty_array(2)
  endif

  if (belowpos<1) then
    uncertainty_array(1) = uncertainty_array(2) ! no lower limit so the negative uncertainty is equal to the value
  else
    uncertainty_array(1) = uncertainty_array(2) - input_array(belowpos)
  endif

!write out the binned results with the mode+-uncertainties at the top

OPEN(850, FILE=trim(filename)//"_"//trim(suffix)//"_binned", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
write(unit = 850,FMT=*) uncertainty_array(2),uncertainty_array(2)-uncertainty_array(1), uncertainty_array(3)+uncertainty_array(2)
write(unit = 850,FMT=*)
do i=1,ii-1
  write(unit = 850,FMT=*) binned_quantity_result(i,1), int(binned_quantity_result(i,2))
end do
close(850)

deallocate(binned_quantity_result)

else !all results are identical
  uncertainty_array(1) = 0.D0
  uncertainty_array(2) = input_array(1)
  uncertainty_array(3) = 0.D0
endif

if (maxval(uncertainty_array) .gt. 0) then !if this condition is not true, array will be full of zeroes
  write (650,itemformat) itemtext,uncertainty_array(2),uncertainty_array(3),-uncertainty_array(1)
else
  write (650,*) itemtext,"--"
endif

end subroutine get_uncertainties

character*10 function gettime()

character*10 :: time

  call DATE_AND_TIME(TIME=time)
  gettime = time(1:2)//":"//time(3:4)//":"//time(5:6)
  return

end function gettime

end program
