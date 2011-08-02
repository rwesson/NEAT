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

subroutine fileread(linelist, filename)
        REAL*8, DIMENSION(3,335), intent(OUT) :: linelist 
        CHARACTER*80 :: filename 
        INTEGER :: IO, I
        DOUBLE PRECISION :: temp,temp2,temp3

        I = 1
        OPEN(199, file=filename, iostat=IO, status='old')
                DO WHILE (IO >= 0)
                        READ(199,*,end=110) temp, temp2, temp3
                        linelist(:,I) =  (/ temp, temp2, temp3 /)
                        I = I + 1
                END DO
        CLOSE(199)

        110 PRINT*, "Reached end of file, I = ", I

        if(linelist(1,1) == 0)then
                PRINT*, "Cheese shop error: no inputs"
                STOP
        endif

        linelist(:,(I+1):335) = 0

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


subroutine element_assign(ILs, linelist, Iint)
        TYPE(line), DIMENSION(62), INTENT(OUT) :: ILs
        REAL*8, DIMENSION(3,335), INTENT(IN) :: linelist
        INTEGER, INTENT(IN) :: Iint
        INTEGER :: i, j

        do i = 1, Iint
                do j = 1, 335 
                        if(linelist(1,j) == ILs(i)%wavelength)then 
                                ILs(i)%intensity = linelist(2,j)
                                ILs(i)%int_err   = linelist(3,j) 
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

end module mod_abundmaths 
