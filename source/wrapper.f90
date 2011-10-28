module mod_abundtypes

TYPE line
        CHARACTER*11 :: name
        CHARACTER*20 :: ion
        DOUBLE PRECISION :: wavelength
        DOUBLE PRECISION :: intensity
        DOUBLE PRECISION :: int_err
        DOUBLE PRECISION :: freq
        DOUBLE PRECISION :: int_dered
        CHARACTER*20 :: transition
        DOUBLE PRECISION :: abundance
        CHARACTER*4 :: zone
END TYPE

end module mod_abundtypes

module mod_resultarrays

TYPE resultarray
        double precision :: NC_abund_CEL
        double precision :: cii_abund_CEL
        double precision :: ciii_abund_CEL
        double precision :: civ_abund_CEL
        double precision :: C_abund_CEL
        double precision :: Nii_abund_CEL
        double precision :: Niii_abund_CEL
        double precision :: Niv_abund_CEL
        double precision :: Nv_abund_CEL
        double precision :: N_abund_CEL
        double precision :: NO_abund_CEL
        double precision :: Oii_abund_CEL
        double precision :: Oiii_abund_CEL
        double precision :: Oiv_abund_CEL
        double precision :: O_abund_CEL
        double precision :: Neii_abund_CEL
        double precision :: Neiii_abund_CEL
        double precision :: Neiv_abund_CEL
        double precision :: Nev_abund_CEL
        double precision :: Ne_abund_CEL
        double precision :: Ariii_abund_CEL
        double precision :: Ariv_abund_CEL
        double precision :: Arv_abund_CEL
        double precision :: Ar_abund_CEL
        double precision :: Sii_abund_CEL
        double precision :: Siii_abund_CEL
        double precision :: S_abund_CEL
        double precision :: He_abund_ORL
        double precision :: C_abund_ORL
        double precision :: N_abund_ORL
        double precision :: O_abund_ORL
        double precision :: Ne_abund_ORL
        double precision :: OII_density
        double precision :: SII_density
        double precision :: low_density
        double precision :: OII_temp
        double precision :: NII_temp
        double precision :: SII_temp
        double precision :: OI_temp
        double precision :: CI_temp
        double precision :: low_temp
        double precision :: ClIII_density
        double precision :: ArIV_density
        double precision :: CIII_density
        double precision :: OIII_IR_density
        double precision :: SIII_IR_density
        double precision :: ArIII_IR_density
        double precision :: NeIII_IR_density
        double precision :: med_density
        double precision :: OIII_temp
        double precision :: OIII_IR_temp
        double precision :: NeIII_temp
        double precision :: NeIII_IR_temp
        double precision :: ArIII_temp
        double precision :: SIII_temp
        double precision :: med_temp
        double precision :: NeIV_density
        double precision :: high_density
        double precision :: ArV_temp
        double precision :: NeV_temp
        double precision :: high_temp
        double precision :: mean_cHb
end type

end module mod_resultarrays

program wrapper

        use mod_abundtypes
        use mod_resultarrays

        CHARACTER*10 :: temp
        CHARACTER :: switch_ext !switch for extinction laws
        INTEGER :: I, runs, doublext, Narg !runs = number of runs for randomiser
        character*6 :: no

        !file reading variables

        TYPE(LINE),DIMENSION(:), allocatable :: linelist 
        TYPE(LINE),DIMENSION(:), allocatable :: linelist_original
        CHARACTER*80 :: filename 
        CHARACTER*1 :: null
        INTEGER :: IO, listlength
        DOUBLE PRECISION :: temp1,temp2,temp3
        type(resultarray), dimension(:), allocatable :: all_results
        type(resultarray), dimension(1) :: iteration_result

        !read command line arguments

        Narg = IARGC() !count input arguments

        if(Narg < 2) then
                PRINT*, "Incorrect arguments"
                PRINT*, "syntax: ./abundances.exe #iterations linelist_file"
                PRINT*, " "
                STOP
        endif

        CALL getarg(1,temp) !get info from input arguments
        read (temp,*) runs 
        CALL getarg(2,filename) 
        if(Narg > 2) then !check additional input arguments for switches (currently extinction laws only)
               CALL getarg (3, temp) !get argument for extinction laws, allowed values -How for Howarth LMC, -CCM for CCM galactic, -Pre for Prevot SMC, default is Seaton/Howarth galactic.
               if (temp == "-How") then
                       switch_ext = "H"
               elseif (temp == "-CCM") then
                       switch_ext = "C"
               elseif (temp == "-Pre") then
                       switch_ext = "P"
               elseif (temp == "-Fit") then
                       switch_ext = "F"
               else
                       switch_ext = "S"
               endif
        else
               switch_ext = "S"
        endif

        !first, read in the line list 

        print *,"Initialising"
        print *,"------------"

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
                call abundances(linelist, 1, switch_ext, listlength, filename, iteration_result)

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

                        call randomizer(linelist, listlength)
                        call abundances(linelist, 0, switch_ext, listlength, filename, iteration_result)
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


        else
                print*, "I didn't want to be a barber anyway. I wanted to be... a lumberjack!   Also, a positive number of runs helps.."
        endif


contains

        subroutine randomizer(linelist, listlength)

                TYPE(line), dimension(listlength) :: linelist 
                INTEGER :: IO, I, j, listlength
                DOUBLE PRECISION :: temp1,temp2,temp3,temp4

                REAL :: fn_val

                !     Local variables
                REAL     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           &
                            r1 = 0.27597, r2 = 0.27846, u, v, x, y, q
                REAL :: half

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
                                temp4=linelist(j)%intensity+(fn_val*linelist(j)%int_err)
                                if(temp4 < 0) temp4 = 0
                                linelist(j)%intensity = temp4 
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

end program
