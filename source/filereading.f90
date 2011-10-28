module mod_abundIO
use mod_abundtypes
implicit none!

contains

subroutine read_ilines(ILs, Iint)        
        TYPE(line), DIMENSION(:), allocatable :: ILs
        INTEGER :: Iint, Iread

        Iint = 1

        301 FORMAT(A11, 1X, A6, 1X, F7.2, 1X, A20,1X,A4)
        OPEN(201, file="source/Ilines_levs", status='old')
                READ (201,*) Iread
                ALLOCATE (ILs(Iread))
                ILs%intensity=0.D0 !otherwise it seems you can get random very small numbers in the array.
                DO WHILE (Iint .le. Iread)!(.true.)
                        READ(201,301,end=401) ILs(Iint)%name, ILs(Iint)%ion, ILs(Iint)%wavelength, ILs(Iint)%transition ,ILs(Iint)%zone!end condition breaks loop.  
                        Iint = Iint + 1
                END DO
                Iint = Iint - 1 !count ends up one too high
                401 PRINT "(A19,I3,A5)", " Read in CEL list, ",Iint," lines"
        CLOSE(201)
end subroutine        

end module

module mod_abundmaths
use mod_abundtypes
implicit none!

contains

!this fantastically ugly function gets the location of certain ions in the important ions array using their name as a key.

integer function get_ion(ionname, iontable, Iint)
        CHARACTER*11 :: ionname
        TYPE(line), DIMENSION(:) :: iontable 
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


subroutine element_assign(ILs, linelist, Iint, listlength)
        TYPE(line), DIMENSION(:), INTENT(OUT) :: ILs
        TYPE(line), DIMENSION(:) :: linelist 
        INTEGER, INTENT(IN) :: Iint, listlength
        INTEGER :: i, j

        do i = 1, Iint
                do j = 1, listlength
                        if(linelist(j)%wavelength == ILs(i)%wavelength)then 
                                ILs(i)%intensity = linelist(j)%intensity
                                ILs(i)%int_err   = linelist(j)%int_err
                                cycle
                        endif        
                end do 
        end do

end subroutine

subroutine get_H(H_BS, linelist, listlength)
        TYPE(line), DIMENSION(4), INTENT(OUT) :: H_BS
        TYPE(line), DIMENSION(:) :: linelist 
        INTEGER :: i, j, listlength
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

                do j = 1, listlength
                         if (linelist(j)%wavelength-HW==0) then
                                H_BS(i)%name = blank
                                H_BS(i)%wavelength = linelist(j)%wavelength
                                H_BS(i)%intensity = linelist(j)%intensity
                                H_BS(i)%int_err = linelist(j)%int_err
                        endif
                end do
        end do

end subroutine

subroutine get_He(He_lines, linelist,listlength)
        TYPE(line), DIMENSION(4), INTENT(OUT) :: He_lines
        TYPE(line), DIMENSION(:), INTENT(IN) :: linelist
        INTEGER :: i, j, listlength
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
                do j = 1, listlength
                        if(linelist(j)%wavelength == HW) then 
                                He_lines(i)%intensity = linelist(j)%intensity
                                He_lines(i)%int_err = linelist(j)%int_err
                        endif
                end do
        end do

end subroutine

!extinction laws now in mod_extinction

end module mod_abundmaths 
