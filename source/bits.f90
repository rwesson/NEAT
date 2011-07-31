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

!extinction laws now in mod_extinction
        
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
