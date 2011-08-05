program wrapper

        CHARACTER*80 :: fname1, ofname !filenames - ofname = output from randomiser, fname1 = input line list
        CHARACTER*10 :: temp
        CHARACTER :: switch_ext !switch for extinction laws
        INTEGER :: I, runs, doublext, Narg !runs = number of runs for randomiser
        character*6 :: no
        Narg = IARGC() !count input arguments
        CALL getarg(1,temp) !get info from input arguments
        read (temp,*) runs 
        CALL getarg(2,fname1) 
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

        !print*, runs

        call init_random_seed()!sets seed for randomiser
        ofname = "r_out"
                !"ball1                                                                           " !length of 80char
        if(runs > 1)then

                !open/create files here for adundances
                OPEN(841, FILE=trim(fname1)//"_NC_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(842, FILE=trim(fname1)//"_C_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(843, FILE=trim(fname1)//"_Nii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(844, FILE=trim(fname1)//"_N_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(845, FILE=trim(fname1)//"_NO_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(846, FILE=trim(fname1)//"_Oii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(847, FILE=trim(fname1)//"_Oiii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(848, FILE=trim(fname1)//"_O_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(849, FILE=trim(fname1)//"_Neiii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(850, FILE=trim(fname1)//"_Neiv_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(851, FILE=trim(fname1)//"_Nev_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(852, FILE=trim(fname1)//"_Ne_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(853, FILE=trim(fname1)//"_Ariii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(854, FILE=trim(fname1)//"_Ariv_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(855, FILE=trim(fname1)//"_Arv_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(856, FILE=trim(fname1)//"_Ar_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(857, FILE=trim(fname1)//"_Sii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(858, FILE=trim(fname1)//"_Siii_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(859, FILE=trim(fname1)//"_S_abund_CEL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(860, FILE=trim(fname1)//"_He_abund_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(861, FILE=trim(fname1)//"_C_abund_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(862, FILE=trim(fname1)//"_N_abund_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(863, FILE=trim(fname1)//"_O_abund_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(864, FILE=trim(fname1)//"_Ne_abund_ORL", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(865, FILE=trim(fname1)//"_[OII]_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(866, FILE=trim(fname1)//"_[SII]_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(867, FILE=trim(fname1)//"_low_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(868, FILE=trim(fname1)//"_[OII]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(869, FILE=trim(fname1)//"_[NII]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(870, FILE=trim(fname1)//"_[SII]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(871, FILE=trim(fname1)//"_[OI]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(872, FILE=trim(fname1)//"_[CI]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(873, FILE=trim(fname1)//"_low_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(874, FILE=trim(fname1)//"_[ClIII]_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(875, FILE=trim(fname1)//"_[ArIV]_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(876, FILE=trim(fname1)//"_CIII]_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(877, FILE=trim(fname1)//"_med_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(878, FILE=trim(fname1)//"_[OIII]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(879, FILE=trim(fname1)//"_[NeIII]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(880, FILE=trim(fname1)//"_[ArIII]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(881, FILE=trim(fname1)//"_[SIII]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(882, FILE=trim(fname1)//"_med_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(883, FILE=trim(fname1)//"_[NeIV]_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(884, FILE=trim(fname1)//"_high_density", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(885, FILE=trim(fname1)//"_[ArV]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(886, FILE=trim(fname1)//"_[NeV]_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(887, FILE=trim(fname1)//"_high_temp", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                OPEN(888, FILE=trim(fname1)//"_mean_cHb", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
                DO I=1,runs
                        print*, "-=-=-=-=-=-=-=-"
                        print*, "iteration ", i
                        print*, "=-=-=-=-=-=-=-="
                        print*, " "
                        write (no, '(I6)') i
                        no = adjustl(no)
                        ofname = "r_out"//no
                        print*, ofname

                        call randomizer(fname1, ofname)
                        call abundances(ofname, 0, switch_ext)!, doublext)
                        call system("rm "//ofname)
                END DO
                DO I=841,864
                        CLOSE(unit=I)
                END DO
        else if(runs == 1)then !calculates abundances without uncertainties
                call abundances(fname1, 1, switch_ext)!, doublext)
        else
                print*, "I didn't want to be a barber anyway. I wanted to be... a lumberjack!   Also, a positive number of runs helps.."
        endif


contains

        subroutine randomizer(filename, outfilename)

                TYPE line
                        CHARACTER*11 :: name
                        CHARACTER*20 :: ion
                        DOUBLE PRECISION :: wavelength
                        DOUBLE PRECISION :: intensity
                        DOUBLE PRECISION :: int_err
                        DOUBLE PRECISION :: freq
                        DOUBLE PRECISION :: int_dered
                        DOUBLE PRECISION :: int_dered_err
                        CHARACTER*20 :: transition
                        DOUBLE PRECISION :: abundance
                        DOUBLE PRECISION :: abund_err
                        CHARACTER*4 :: zone
                END TYPE


                TYPE(line), dimension(:), allocatable :: linelist
                CHARACTER*80, INTENT(IN) :: filename, outfilename
                character*1 :: null
                !character*10 :: Ich, mcnumberch
                INTEGER :: IO, I, j, listlength!, mcnumber
                DOUBLE PRECISION :: temp,temp2,temp3,temp4

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
                !PRINT*, filename

        ! first get the number of lines

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

                REWIND (199)
                DO I=1,listlength
                        READ(199,*,end=110) temp, temp2, temp3
                        linelist(i)%wavelength = temp
                        linelist(i)%intensity = temp2
                        linelist(i)%int_err = temp3
                END DO
                CLOSE(199)

                110 PRINT*, "Successfully read ", I," lines"

                if(linelist(1)%wavelength == 0)then
                        PRINT*, "Cheese shop error: no inputs"
                        STOP
                endif


                !create [mcnumber] of files, each with a randomly varied flux

                !do i = 1,mcnumber
                        !write (ich,'(I3)') i
                        !open(unit=101, name='randomized_'//Ich//'.dat',type='unknown')
                        print*, outfilename

                        open(unit=401, file=outfilename,action="WRITE")!,type='unknown')
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
                                write (401,"(F7.2,1X,F10.5,1X,F7.2)") linelist(j)%wavelength, temp4, linelist(j)%int_err
                        end do
                        close(unit=401)
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
