program wrapper

        CHARACTER*80 :: fname1, fname2, fname3, ofname !filenames - ofname = output from randomiser, fname1 2 3 = inputs
	CHARACTER*10 :: temp
	INTEGER :: I, runs, doublext !runs = number of runs for randomiser
	character*6 :: no
	CALL getarg(1,temp) !get info from input arguments
	read (temp,*) runs
	CALL getarg(2,temp)
	read (temp,*) doublext	
	CALL getarg(3,fname1) 
        CALL getarg(4,fname2)
        CALL getarg(5,fname3)
	
	!print*, runs
	
	call init_random_seed()!sets seed for randomiser
	ofname = "r_out"
                !"ball1                                                                           " !length of 80char
	if(runs > 1)then 	
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
			call abundances(ofname, fname2, fname3, 0, doublext)
			call system("rm "//ofname)
		END DO
	else if(runs == 1)then !calculates abundances without uncertainties
		call abundances(fname1, fname2, fname3, 1, doublext)
	else
		print*, "I didn't want to be a barber anyway. I wanted to be... a lumberjack!   Also, a positive number of runs helps.."	
	endif		
	
	
contains

	subroutine randomizer(filename, outfilename)

		TYPE line
			DOUBLE PRECISION :: wavelength
			DOUBLE PRECISION :: flux
			DOUBLE PRECISION :: uncertainty
		end TYPE
		double precision, dimension(3,335) :: linelist
		CHARACTER*80, INTENT(IN) :: filename, outfilename
		!character*10 :: Ich, mcnumberch
		INTEGER :: IO, I, j, nlines!, mcnumber
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
		
		REWIND 902
		OPEN(902, file=filename, iostat=IO, action="READ")
	        DO WHILE (IO>=0)
	                READ(902,*,end=101) temp, temp2, temp3
	                linelist(:,I) = (/ temp, temp2, temp3 /)
	                I = I + 1
	        END DO
		CLOSE(902)

		101 PRINT*, "done reading lines, no. = ", I
	
		nlines = I-1

		!create [mcnumber] of files, each with a randomly varied flux
	
		!do i = 1,mcnumber
			!write (ich,'(I3)') i
			!open(unit=101, name='randomized_'//Ich//'.dat',type='unknown')
			print*, outfilename

			open(unit=401, file=outfilename,action="WRITE")!,type='unknown')
			do j = 1,nlines
		
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
				temp4=linelist(2,j)+(fn_val*linelist(3,j))
				if(temp4 < 0) temp4 = 0
				write (401,"(F7.2,1X,F10.5,1X,F7.2)") linelist(1,j), temp4, linelist(3,j)
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
