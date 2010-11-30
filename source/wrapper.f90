program wrapper

        CHARACTER*80 :: fname1, fname2, fname3, ofname
	CHARACTER*10 :: temp
	INTEGER :: I, runs
	character*4 :: no
	CALL getarg(1,temp)
	read (temp,*) runs
	CALL getarg(2,fname1) 
        CALL getarg(3,fname2)
        CALL getarg(4,fname3)
	
	!print*, runs
	
	call init_random_seed()
	ofname = "r_out"
                !"ball1                                                                           " !length of 80char
		
	DO I=1,runs
		print*, "-=-=-=-=-=-=-=-"
		print*, "iteration ", i 
		print*, "=-=-=-=-=-=-=-="
		print*, " "
		write (no, '(I4)') i
		no = adjustl(no)
		ofname = "r_out"//no
		print*, ofname

		call randomizer(fname1, ofname)
		call abundances(ofname, fname2, fname3, 0)
		call system("rm "//ofname)
	END DO
	
	
	
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
		DOUBLE PRECISION :: temp,temp2,temp3

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

			open(unit=401, name=outfilename,action="WRITE")!,type='unknown')
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
				!print*, fn_val
				write (401,"(F7.2,1X,F10.5,1X,F7.2)") linelist(1,j), linelist(2,j)+(fn_val*linelist(3,j)), 0.0
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
