
program wrapper

        use mod_abundtypes
        use mod_resultarrays
        use mod_extinction

        CHARACTER*10 :: temp, tempa
        CHARACTER :: switch_ext !switch for extinction laws
        INTEGER :: I, runs, Narg !runs = number of runs for randomiser
        character*6 :: no

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
        DOUBLE PRECISION :: temp1,temp2,temp3, mean_ext, R
        type(resultarray), dimension(:), allocatable :: all_results
        type(resultarray), dimension(1) :: iteration_result

        make_dered_ll=0
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
         enddo

         if (Narg .eq. 1) then
           filename=trim(options(1))
         endif

         deallocate(options)

         if (filename=="") then
                print *,"Error: No input file specified"
                stop
         endif

        ! read all arguments into array
        ! loop through and read the relevant variables
        ! print warning for unrecognised options
        ! options are:
        !  -n / --n-iterations    : number of iterations (default 1)
        !  -i / --input           : input file (no default)
        !  -e / --extinction-law  : extinction law (default Howarth 1983)
        !  to be implemented:
        !  -R                     : R (default 3.1)
        !  -nelow / --density-low : low ionisation zone density (default - calculate from line list)
        !  -telow / --temperature-low : low i. zone temperature
        !  -nemed / --density-med : medium i. density
        !  -temed / --temperature-med : medium i. temperature
        !  -nehigh / --density-high : high i. density
        !  -tehigh / --temperature-high : high i. temperature
        !  -chb                   : value of c(Hb), the logarithmic extinction at H beta


        !first, read in the line list

        print *,"Initialising"
        print *,"------------"

        call DATE_AND_TIME(date,time)
        print *
        print *,"Start time: ",time(1:2),":",time(3:4),":",time(5:10)," on ",date(7:8),"/",date(5:6),"/",date(1:4)
        print *,"Input file: ",filename

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

        !now check number of iterations.  If 1, line list is fine as is.  If more than one, randomize the fluxes

        if(runs == 1)then !calculates abundances without uncertainties
                call abundances(linelist, 1, switch_ext, listlength, filename, iteration_result, R)

                if(make_dered_ll ==  1)then
                        mean_ext=DBLE(0)
                        CALL deredden_ll(switch_ext, linelist, listlength, mean_ext )
                endif


        else if(runs > 1)then

                linelist_original = linelist

                call init_random_seed()!sets seed for randomiser
                allocate(all_results(runs))

                !open/create files here for adundances

                DO I=1,runs
                        print*, "-=-=-=-=-=-=-=-"
                        print*, "iteration ", i
                        print*, "=-=-=-=-=-=-=-="
                        print*, " "
                        write (no, '(I6)') i

                        call randomizer(linelist, listlength, R)
                        R=3.1 ! no randomisation
                        call abundances(linelist, 0, switch_ext, listlength, filename, iteration_result, R)
                        linelist = linelist_original
                        all_results(i)=iteration_result(1)
                END DO

                OPEN(841, FILE=trim(filename)//"_NC_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(842, FILE=trim(filename)//"_C_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(843, FILE=trim(filename)//"_Nii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(844, FILE=trim(filename)//"_N_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(845, FILE=trim(filename)//"_NO_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(846, FILE=trim(filename)//"_Oii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(847, FILE=trim(filename)//"_Oiii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(848, FILE=trim(filename)//"_O_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(849, FILE=trim(filename)//"_Neiii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(850, FILE=trim(filename)//"_Neiv_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(851, FILE=trim(filename)//"_Nev_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(852, FILE=trim(filename)//"_Ne_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(853, FILE=trim(filename)//"_Ariii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(854, FILE=trim(filename)//"_Ariv_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(855, FILE=trim(filename)//"_Arv_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(856, FILE=trim(filename)//"_Ar_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(857, FILE=trim(filename)//"_Sii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(858, FILE=trim(filename)//"_Siii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(859, FILE=trim(filename)//"_S_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(860, FILE=trim(filename)//"_He_abund_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(861, FILE=trim(filename)//"_C_abund_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(862, FILE=trim(filename)//"_N_abund_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(863, FILE=trim(filename)//"_O_abund_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(864, FILE=trim(filename)//"_Ne_abund_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(865, FILE=trim(filename)//"_[OII]_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(866, FILE=trim(filename)//"_[SII]_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(867, FILE=trim(filename)//"_low_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(868, FILE=trim(filename)//"_[OII]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(869, FILE=trim(filename)//"_[NII]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(870, FILE=trim(filename)//"_[SII]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(871, FILE=trim(filename)//"_[OI]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(872, FILE=trim(filename)//"_[CI]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(873, FILE=trim(filename)//"_low_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(874, FILE=trim(filename)//"_[ClIII]_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(875, FILE=trim(filename)//"_[ArIV]_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(876, FILE=trim(filename)//"_CIII]_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(877, FILE=trim(filename)//"_med_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(878, FILE=trim(filename)//"_[OIII]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(879, FILE=trim(filename)//"_[NeIII]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(880, FILE=trim(filename)//"_[ArIII]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(881, FILE=trim(filename)//"_[SIII]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(882, FILE=trim(filename)//"_med_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(883, FILE=trim(filename)//"_[NeIV]_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(884, FILE=trim(filename)//"_high_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(885, FILE=trim(filename)//"_[ArV]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(886, FILE=trim(filename)//"_[NeV]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(887, FILE=trim(filename)//"_high_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(888, FILE=trim(filename)//"_mean_cHb", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')

                do i=1,runs
                        write(unit = 841,FMT=*) all_results(i)%NC_abund_CEL
                        write(unit = 842,FMT=*) all_results(i)%C_abund_CEL
                        write(unit = 843,FMT=*) all_results(i)%nii_abund_CEL
                        write(unit = 844,FMT=*) all_results(i)%N_abund_CEL
                        write(unit = 845,FMT=*) all_results(i)%NO_abund_CEL
                        write(unit = 846,FMT=*) all_results(i)%oii_abund_CEL
                        write(unit = 847,FMT=*) all_results(i)%oiii_abund_CEL
                        write(unit = 848,FMT=*) all_results(i)%O_abund_CEL
                        write(unit = 849,FMT=*) all_results(i)%neiii_abund_CEL
                        write(unit = 850,FMT=*) all_results(i)%neiv_abund_CEL
                        write(unit = 851,FMT=*) all_results(i)%nev_abund_CEL
                        write(unit = 852,FMT=*) all_results(i)%Ne_abund_CEL
                        write(unit = 853,FMT=*) all_results(i)%ariii_abund_CEL
                        write(unit = 854,FMT=*) all_results(i)%ariv_abund_CEL
                        write(unit = 855,FMT=*) all_results(i)%arv_abund_CEL
                        write(unit = 856,FMT=*) all_results(i)%Ar_abund_CEL
                        write(unit = 857,FMT=*) all_results(i)%sii_abund_CEL
                        write(unit = 858,FMT=*) all_results(i)%siii_abund_CEL
                        write(unit = 859,FMT=*) all_results(i)%S_abund_CEL
                        write(unit = 860,FMT=*) all_results(i)%He_abund_ORL
                        write(unit = 861,FMT=*) all_results(i)%C_abund_ORL
                        write(unit = 862,FMT=*) all_results(i)%N_abund_ORL
                        write(unit = 863,FMT=*) all_results(i)%O_abund_ORL
                        write(unit = 864,FMT=*) all_results(i)%Ne_abund_ORL
                        write(unit = 865,FMT=*) all_results(i)%oii_density
                        write(unit = 866,FMT=*) all_results(i)%sii_density
                        write(unit = 867,FMT=*) all_results(i)%low_density
                        write(unit = 868,FMT=*) all_results(i)%nii_temp
                        write(unit = 869,FMT=*) all_results(i)%oii_temp
                        write(unit = 870,FMT=*) all_results(i)%sii_temp
                        write(unit = 871,FMT=*) all_results(i)%oi_temp
                        write(unit = 872,FMT=*) all_results(i)%ci_temp
                        write(unit = 873,FMT=*) all_results(i)%low_temp
                        write(unit = 874,FMT=*) all_results(i)%cliii_density
                        write(unit = 875,FMT=*) all_results(i)%ariv_density
                        write(unit = 876,FMT=*) all_results(i)%ciii_density
                        write(unit = 877,FMT=*) all_results(i)%med_density
                        write(unit = 878,FMT=*) all_results(i)%oiii_temp
                        write(unit = 879,FMT=*) all_results(i)%neiii_temp
                        write(unit = 880,FMT=*) all_results(i)%ariii_temp
                        write(unit = 881,FMT=*) all_results(i)%siii_temp
                        write(unit = 882,FMT=*) all_results(i)%med_temp
                        write(unit = 883,FMT=*) all_results(i)%neiv_density
                        write(unit = 884,FMT=*) all_results(i)%high_density
                        write(unit = 885,FMT=*) all_results(i)%arv_temp
                        write(unit = 886,FMT=*) all_results(i)%nev_temp
                        write(unit = 887,FMT=*) all_results(i)%high_temp
                        write(unit = 888,FMT=*) all_results(i)%mean_cHb
                end do

                DO I=841,888
                        CLOSE(unit=I)
                END DO

                if(make_dered_ll ==  1)then
                        mean_ext=DBLE(0)! fix this.
                        CALL deredden_ll(switch_ext, linelist, listlength, mean_ext )
                endif


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
                DOUBLE PRECISION :: temp1,temp2,temp3,temp4, R

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

                                if (j==1) R=3.1+(0.15*fn_val)

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
                PRINT*, "Randomizing complete"
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


end program
