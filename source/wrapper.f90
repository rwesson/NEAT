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

program wrapper

        use mod_abundtypes

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

        !read command line arguments

        Narg = IARGC() !count input arguments
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
        111 print *,"File contains ",I," lines"
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

        110 PRINT*, "Successfully read ", I," lines"

        if(linelist(1)%wavelength == 0)then
                PRINT*, "Cheese shop error: no inputs"
                STOP
        endif

        linelist_original = linelist

        !now check number of iterations.  If 1, line list is fine as is.  If more than one, randomize the fluxes

        if(runs == 1)then !calculates abundances without uncertainties
                call abundances(linelist, 1, switch_ext, listlength, filename)

        else if(runs > 1)then

                call init_random_seed()!sets seed for randomiser 

                !open/create files here for adundances
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
                DO I=1,runs
                        print*, "-=-=-=-=-=-=-=-"
                        print*, "iteration ", i
                        print*, "=-=-=-=-=-=-=-="
                        print*, " "
                        write (no, '(I6)') i

                        call randomizer(linelist, listlength)
                        call abundances(linelist, 0, switch_ext, listlength, filename)
                        linelist = linelist_original

                END DO
                DO I=841,864
                        CLOSE(unit=I)
                END DO 
        else
                print*, "I didn't want to be a barber anyway. I wanted to be... a lumberjack!   Also, a positive number of runs helps.."
        endif


contains

        subroutine randomizer(linelist, listlength)

                TYPE(line), dimension(listlength) :: linelist 
                INTEGER :: IO, I, j, listlength!, mcnumber
                DOUBLE PRECISION :: temp1,temp2,temp3,temp4

                REAL :: fn_val

                !     Local variables
                REAL     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           &
                            r1 = 0.27597, r2 = 0.27846, u, v, x, y, q
                REAL :: half
                !call init_random_seed()
                !read in the file

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
