module mod_abundtypes

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

end module mod_abundtypes

module mod_abundIO
use mod_abundtypes
implicit none!

contains

subroutine read_ilines(ILs, Iint)        
        TYPE(line), DIMENSION(62) :: ILs
        INTEGER :: Iint 

        print *, "Initialisation"
        print *, "=============="

        Iint = 1

        301 FORMAT(A11, 1X, A6, 1X, F7.2, 1X, A20,1X,A4)
        OPEN(201, file="source/Ilines_levs", status='old')
                DO WHILE (Iint < 62)!(.true.)
                        READ(201,301,end=401) ILs(Iint)%name, ILs(Iint)%ion, ILs(Iint)%wavelength, ILs(Iint)%transition ,ILs(Iint)%zone!end condition breaks loop.  
			Iint = Iint + 1
                END DO
                401 PRINT*, "done reading important lines, Iint = ", Iint 
        CLOSE(201)
end subroutine        

subroutine fileread(optical, infrared, ultraviolet, filename1, filename2, filename3)
        REAL*8, DIMENSION(3,335), intent(OUT) :: optical, infrared, ultraviolet
        REAL*8, DIMENSION(3,335) ::  temp1, temp2, temp3
        CHARACTER*80 :: filename1, filename2, filename3
        CHARACTER*20 :: message
        INTEGER :: I
                
        if(filename1 .ne. "          ")then
		temp1 = readfile(filename1)
                call w_assign(optical, infrared, ultraviolet, temp1, message)
                PRINT*, message, filename1
        endif
                
        if(filename2 .ne. "          ")then
                temp2 = readfile(filename2)
                call w_assign(optical, infrared, ultraviolet, temp2, message)
                PRINT*, message, filename2
        endif
                
        if(filename3 .ne. "          ")then
                temp3 = readfile(filename3)
                call w_assign(optical, infrared, ultraviolet, temp3, message)
                PRINT*, message, filename3
        endif
        
        if(optical(1,1) == 0)then
                PRINT*, "Cheese shop error: no inputs"
                STOP
        endif
        
        !do I = 1, 329
        !PRINT*, optical(:,I), temp1(:,I)
        !END DO
        
end subroutine

function readfile(filename)
        REAL*8, DIMENSION(3,335) :: readfile
        CHARACTER*80 :: filename
        INTEGER :: IO, I
        DOUBLE PRECISION :: temp,temp2,temp3

        I = 1
        OPEN(199, file=filename, iostat=IO, status='old')
                DO WHILE (IO >= 0)
                        READ(199,*,end=110) temp, temp2, temp3
                        readfile(:,I) =  (/ temp, temp2, temp3 /) 
                        I = I + 1
                END DO
        CLOSE(199)

        110 PRINT*, "Reached end of file, I = ", I
        
	readfile(:,(I+1):335) = 0
	        
	
end function

subroutine w_assign(optical, infrared, ultraviolet, temp, message)
        REAL*8, DIMENSION(3,335) :: optical, infrared, ultraviolet, temp
        CHARACTER*20 :: message 
                
        if(temp(1,1) < 100)then
                infrared = temp
                infrared(1,:) = infrared(1,:) * 10000 !converting microns to angstroms
                message = "file = infrared"
        
        elseif(temp(1,1) > 3000)then
                
                optical = temp
                message = "file = optical"
        
        else
                
                ultraviolet = temp
                message = "file = ultraviolet"
        
        endif        

end subroutine

end module


module mod_abundmaths
use mod_abundtypes
implicit none!

contains

!this fantastically ugly function gets the location of certain ions in the important ions array using their name as a key.

integer function get_ion(ionname, iontable, Iint)
        CHARACTER*11 :: ionname
        TYPE(line), DIMENSION(62) :: iontable 
        INTEGER :: i
	INTEGER, INTENT(IN) :: Iint
        
        do i = 1, Iint
                
                !PRINT*, trim(iontable(i)%name), trim(ionname)
                
                if(trim(iontable(i)%name) == trim(ionname))then
                        get_ion = i
                        return
                endif
        end do

        PRINT*, "Nudge Nudge, wink, wink error. Ion not found, say no more.", ionname
        
end function        
        

subroutine element_assign(ILs, O, IR, UV, Iint)
        TYPE(line), DIMENSION(62), INTENT(OUT) :: ILs
        REAL*8, DIMENSION(3,335), INTENT(IN) :: O, IR, UV
        INTEGER, INTENT(IN) :: Iint
        INTEGER :: i, j

        do i = 1, Iint
                do j = 1, 335 
			if(O(1,j) == ILs(i)%wavelength)then 
                                ILs(i)%intensity = O(2,j)
                                ILs(i)%int_err   = O(3,j) 
                                cycle
                        endif        
                end do
                
                do j = 1, 335
                        if(IR(1,j) == ILs(i)%wavelength)then
                                ILs(i)%intensity = IR(2,j)
                                ILs(i)%int_err   = UV(3,j)
                                cycle
                        endif        
                end do
                
                do j = 1, 335
                        if(UV(1,j) > ILs(i)%wavelength-101 .AND.UV(1,j) < ILs(i)%wavelength+101 )then
                                ILs(i)%intensity = UV(2,j)
                                ILs(i)%int_err   = UV(3,j)
                                cycle
                        endif        
                end do
        end do
        
end subroutine

subroutine get_H(H_BS, O)
        TYPE(line), DIMENSION(4), INTENT(OUT) :: H_BS
        REAL*8, DIMENSION(3,335), INTENT(IN) :: O
        INTEGER :: i, j
        REAL*8 :: HW = 0.00000000
        CHARACTER*10 :: blank 
        !another ugly kludge, but it works.

        do i = 1, 4
                if(i == 1)then
                        blank = "Halpha     "
                        HW = 6562.77D0 
                elseif(i == 2)then
                        blank = "Hbeta      "
                        HW = 4861.33D0
                elseif(i == 3)then
                        blank = "Hgamma     "
                        HW = 4340.47D0
                elseif(i == 4)then
                        blank = "Hdelta     "
                        HW = 4101.74D0
                else
                        PRINT*, "This is an EX-PARROT!!"
                endif        

                do j = 1, 335 
                         if (int(O(1,j)-HW)==0) then
                                H_BS(i)%name = blank
                                H_BS(i)%wavelength = O(1,j)
                                H_BS(i)%intensity = O(2,j)
                                H_BS(i)%int_err = O(3,j) 
                        endif
                end do
        end do
        
end subroutine

subroutine get_He(He_lines, O)
        TYPE(line), DIMENSION(4), INTENT(OUT) :: He_lines
        REAL*8, DIMENSION(3,335), INTENT(IN) :: O
        INTEGER :: i, j
        REAL*8 :: HW
        CHARACTER*10 :: blank
        !another ugly kludge, but it works.  
        do i = 1, 4
                if(i == 1)then
                        blank = "HeII4686   " 
                        HW = 4685.68D0 
                elseif(i == 2)then
                        blank = "HeI4471    "
                        HW = 4471.50D0
                elseif(i == 3)then
                        blank = "HeI5876    "
                        HW = 5875.66D0
                elseif(i == 4)then
                        blank = "HeI6678    "
                        HW = 6678.16D0
                else
                        PRINT*, "This is an EX-PARROT!!"
                endif
                He_lines(i)%name = blank
                He_lines(i)%wavelength = HW
                He_lines(i)%intensity = 0.0
                He_lines(i)%int_err = 0.0
                do j = 1, 335
                        if(O(1,j) == HW)then 
                                He_lines(i)%intensity = O(2,j)
                                He_lines(i)%int_err = O(3,j)
                        endif
                end do
        end do

end subroutine

subroutine calc_extinction_coeffs(H_BS, c1, c2, c3, c1_err, c2_err, c3_err, cerror, meanextinction)
        TYPE(line), DIMENSION(4) :: H_BS
        DOUBLE PRECISION :: c1, c2, c3, meanextinction, cerror
        DOUBLE PRECISION :: c1_err, c2_err, c3_err

if (H_BS(1)%intensity .gt. 0 .and. H_BS(2)%intensity .gt. 0) then
        c1 = log10( ( DBLE(H_BS(1)%intensity) / DBLE(H_BS(2)%intensity) )/2.850 )/0.32
        c1_err = 0.434 * ((H_BS(1)%int_err/H_BS(1)%intensity)**2 + (H_BS(2)%int_err/H_BS(2)%intensity)**2)**0.5
else
        c1 = 0.0
        c1_err = 0.0
endif
if (H_BS(3)%intensity .gt. 0 .and. H_BS(2)%intensity .gt. 0) then
        c2 = log10( ( DBLE(H_BS(3)%intensity) / DBLE(H_BS(2)%intensity) )/0.469 )/(-0.127)
        c2_err = 0.434 * ((H_BS(3)%int_err/H_BS(3)%intensity)**2 + (H_BS(2)%int_err/H_BS(2)%intensity)**2)**0.5
else
        c2=0.0
        c2_err = 0.0
endif
if (H_BS(4)%intensity .gt. 0 .and. H_BS(2)%intensity .gt. 0) then
        c3 = log10( ( DBLE(H_BS(4)%intensity) / DBLE(H_BS(2)%intensity) )/0.259 )/(-0.180)
        c3_err = 0.434 * ((H_BS(4)%int_err/H_BS(4)%intensity)**2 + (H_BS(2)%int_err/H_BS(2)%intensity)**2)**0.5
else
        c3 = 0.0
        c3_err = 0.0
endif

        meanextinction = (c1*H_BS(1)%intensity + c2*H_BS(3)%intensity + c3*H_BS(4)%intensity) / (H_BS(1)%intensity + H_BS(3)%intensity + H_BS(4)%intensity)   
        
        cerror = (H_BS(1)%intensity*(c1_err/c1)**2 + H_BS(3)%intensity*(c2_err/c2)**2 + H_BS(4)%intensity*(c3_err/c3)**2)**0.5 / (H_BS(1)%intensity + H_BS(3)%intensity + H_BS(4)%intensity)
!        cerror = (( (meanextinction-c1)**2 + (meanextinction-c2)**2 + (meanextinction-c3)**2 )/3) ** 0.5

end subroutine

double precision function flambda(X,switch)
        DOUBLE PRECISION :: X
        INTEGER :: switch
        !Howarth 1983
        if(switch == 1) flambda = ((16.07 - (3.20 * X) + (0.2975*X**2))) !far UV
        if(switch == 2) flambda = ((2.19 + (0.848*X) + (1.01/(((X-4.60)**2) + 0.280)))) !mid UV
        if(switch == 3) flambda = ((1.46 + (1.048*X) + (1.01/(((X-4.60)**2) + 0.280)))) !near UV
        if(switch == 4) flambda = ((3.1 + (2.56*(X-1.83)) - (0.993*(X-1.83)**2) )) ! optical
        if(switch == 5) flambda = (( (1.86*X**2) - (0.48*X**3) - (0.1*X))) ! IR
        
        flambda = (flambda / 3.63) - 1

end function        

subroutine deredden(lines, number, m_ext, cerror)
        TYPE(line), DIMENSION(:) :: lines
        INTEGER :: number
        DOUBLE PRECISION :: m_ext, fl, cerror
        INTEGER :: i
        
        do i = 1,number

                lines(i)%freq = DBLE(10000) / DBLE(lines(i)%wavelength)

                if( (lines(i)%freq .gt. 7.14) .AND. (lines(i)%freq .lt. 10.0))then !far UV
                        fl = flambda(lines(i)%freq, 1)
                        lines(i)%int_dered = lines(i)%intensity * 10**(m_ext*fl) 
                        lines(i)%int_dered_err = ((lines(i)%int_err / lines(i)%intensity)**2 + (log(10.0)*cerror*fl)**2)**0.5
		
		elseif((lines(i)%freq .gt. 3.65) .AND. (lines(i)%freq .lt. 7.14))then ! mid UV
                        fl = flambda(lines(i)%freq, 2)
                        lines(i)%int_dered = lines(i)%intensity * 10**(m_ext*fl) 
                        lines(i)%int_dered_err = ((lines(i)%int_err / lines(i)%intensity)**2 + (log(10.0)*cerror*fl)**2)**0.5
		
		elseif((lines(i)%freq .gt. 2.75) .AND. (lines(i)%freq .lt. 3.65))then ! near UV
                        fl = flambda(lines(i)%freq, 3)
                        lines(i)%int_dered = lines(i)%intensity * 10**(m_ext*fl) 
                        lines(i)%int_dered_err = ((lines(i)%int_err / lines(i)%intensity)**2 + (log(10.0)*cerror*fl)**2)**0.5
		
		elseif((lines(i)%freq .gt. 1.83) .AND. (lines(i)%freq .lt. 2.75))then ! optical
                        fl = flambda(lines(i)%freq, 4)
                        
			lines(i)%int_dered = lines(i)%intensity * 10**(m_ext*fl) 
                        lines(i)%int_dered_err = ((lines(i)%int_err / lines(i)%intensity)**2 + (log(10.0)*cerror*fl)**2)**0.5
		
		elseif(lines(i)%freq .lt. 1.83)then !IR
                        fl = flambda(lines(i)%freq, 5)
                        lines(i)%int_dered = lines(i)%intensity * 10**(m_ext*fl) 
                        lines(i)%int_dered_err = ((lines(i)%int_err / lines(i)%intensity)**2 + (log(10.0)*cerror*fl)**2)**0.5
		
		endif
                
        end do

end subroutine

subroutine deredden_O(O_red, O_dered, m_ext, cerror)
        REAL*8, DIMENSION(3,335) :: O_red, O_dered
        DOUBLE PRECISION :: m_ext, fl, cerror
        INTEGER :: I
        REAL*8 :: X

        do I = 1, 335
            X =  DBLE(DBLE(10000)/DBLE(O_red(1,I)))
            if (X .lt. 1.83) then 
                fl = flambda(X, 5) 
                O_dered(:,I) = (/ O_red(1,I), O_red(2,I)*(10**(m_ext*fl)), ((O_red(3,I)/O_red(2,I))**2 + (log(10.0)**(cerror*fl) )**2)**0.5 /)
            elseif (X .gt. 1.83 .and. X .lt. 2.75) then
                fl = flambda(X, 4)
                O_dered(:,I) = (/ O_red(1,I), O_red(2,I)*(10**(m_ext*fl)), O_red(3,I)*(10**(m_ext*fl) ) /)
	    elseif (X .gt. 2.75 .and. X .lt. 3.63) then
                fl = flambda(X, 3)
                O_dered(:,I) = (/ O_red(1,I), O_red(2,I)*(10**(m_ext*fl)), O_red(3,I)*(10**(m_ext*fl) ) /)
            endif 
        end do

end subroutine  

subroutine deredden_2(lines, number, b_ext, r_ext) !this is specific to Dave's NTT/EFOSC data which is poorly calibrated between red/blue. Remove
        TYPE(line), DIMENSION(:) :: lines
        INTEGER :: number
        DOUBLE PRECISION :: b_ext, r_ext, fl
        INTEGER :: i
        
        do i = 1,number

                lines(i)%freq = DBLE(10000) / DBLE(lines(i)%wavelength)

  
                if((lines(i)%freq .gt. 1.83) .AND. (lines(i)%freq .lt. 2.75))then ! optical 3636A - 5464
                        fl = flambda(lines(i)%freq, 4)
                        lines(i)%int_dered = lines(i)%intensity * 10**(b_ext*fl) 
			lines(i)%int_dered_err = lines(i)%int_err / lines(i)%intensity
                elseif(lines(i)%freq .lt. 1.83)then !IR   5464 - end
                        fl = flambda(lines(i)%freq, 5)
                        lines(i)%int_dered = lines(i)%intensity * 10**(r_ext*fl) 
			lines(i)%int_dered_err = lines(i)%int_err / lines(i)%intensity
                endif
                
        end do

end subroutine


      
        
recursive function getTD(temp, dens, dt, dd, ratio1, ratio2, ion_name1, ion_name2, run) result(TD)
        DOUBLE PRECISION, DIMENSION(2) :: TD
        DOUBLE PRECISION, INTENT(IN) :: temp, dens, dt, dd, ratio1, ratio2
        DOUBLE PRECISION :: T, D
        CHARACTER*6 :: ion_name1, ion_name2
        INTEGER :: run
        
        T = temp
        D = dens
        
        run = run * -1
        
        !if(run > 0)call ratioT(ratio1, dens, ion_name1, T) ! roger must create these
        !if(run < 0)call ratioD(ratio2, temp, ion_name2, D)
        
        if( ((T - temp) < dt) .AND. (D-dens < dd))then
                TD = (/T, D/)
        else                
                TD = getTD(TD(1), TD(2), dt, dd, ratio1, ratio2, ion_name1, ion_name2, run)
        endif
        
end function

end module mod_abundmaths 
