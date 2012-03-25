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
                        READ(201,301) ILs(Iint)%name, ILs(Iint)%ion, ILs(Iint)%wavelength, ILs(Iint)%transition ,ILs(Iint)%zone!end condition breaks loop.  
                        Iint = Iint + 1
                END DO
                Iint = Iint - 1 !count ends up one too high 
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
        TYPE(line), DIMENSION(38), INTENT(OUT) :: H_BS
        TYPE(line), DIMENSION(:) :: linelist 
        double precision, dimension(38) :: balmerlines
        INTEGER :: i, j, listlength
        REAL*8 :: HW = 0.00000000 
        !another ugly kludge, but it works.

        balmerlines = (/ 6562.77D0, 4861.33D0, 4340.47D0, 4101.74D0, 3970.07D0, 3889.05D0, 3835.38D0, 3797.90D0, 3770.63D0, 3750.15D0, 3734.37D0, 3721.94D0, 3711.97D0, 3703.85D0, 3697.15D0, 3691.55D0, 3686.83D0, 3682.81D0, 3679.35D0, 3676.36D0, 3673.76D0, 3671.48D0, 3669.46D0, 3667.68D0, 3666.10D0, 3664.68D0, 3663.40D0, 3662.26D0, 3661.22D0, 3660.28D0, 3659.42D0, 3658.64D0, 3657.92D0, 3657.27D0, 3656.66D0, 3656.11D0, 3655.59D0, 3655.12D0 /)

        do i = 1, 38
                HW = balmerlines(i)

                do j = 1, listlength
                         if (linelist(j)%wavelength-HW==0) then
                                H_BS(i)%name = "Hbalmer    "
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

module mod_atomic_read


contains
subroutine read_atomic_data(ion)
use mod_atomicdata
	type(atomic_data) :: ion
	integer :: I,J,K,NCOMS,ID(2),JD(2),KP1,NLEV1
	character*1 :: comments(78)
	real*8 :: GX,WN,AX,QX
	
	OPEN(unit=1, status = 'OLD', file = 'Atomic-data/'//ion%ion(1(INDEX( &
& 		ION%ion,' ') - 1))//'.dat',ACTION='READ')

!read # of comment lines and skip them
	READ(1,*)NCOMS
	do I = 1,NCOMS
		read(1,1003) comments
	end do
	
!read # levels and temps, then allocate arrays
	read(1,*) ion%NLEVS,ion%NTEMPS 
	
	allocate(ion%labels(ion%nlevs))
	allocate(ion%temps(ion%ntemps))
	allocate(ion%roott(ion%ntemps))
	allocate(ion%G(ion%nlevs))
	allocate(ion%waveno(ion%nlevs))
	allocate(ion%col_str(ion%ntemps,ion%nlevs,ion%nlevs))
	allocate(ion%A_coeffs(ion%nlevs,ion%nlevs))
	
	ion%col_str = 0d0
	ion%A_coeffs = 0d0
	ion%G = 0d0
	ion%waveno= 0d0
	
	!read levels and temperatures
	do I = 1,ion%NLEVS
	read(1,1002) ion%labels(I)
	enddo
	
	do I = 1,ion%NTEMPS
	read(1,*) ion%temps(I)
	enddo
	
	!read collision strengths
        QX=1
        K = 1
        DO WHILE (QX .gt. 0)
                READ(1,*) ID(2), JD(2), QX
                IF (QX.eq.0.D0) exit
                if (ID(2) .eq. 0) then
                   ID(2) = ID(1)
                   K = K + 1
                else
                   ID(1) = ID(2)
                   K = 1
                endif
                if (JD(2) .eq. 0) then
                   JD(2) = JD(1)
                else
                   JD(1) = JD(2)
                endif
                if (QX .ne. 0.D0) then
                I = ID(2) 
                J = JD(2) 
                ion%col_str(K,I,J) = QX
                endif
        enddo

    NLEV1 = ion%NLEVS-1
      DO K = 1,NLEV1
        KP1 = K + 1 
          DO L = KP1, ion%NLEVS
            READ (1,*) I, J, AX  !read transition probabilities
            ion%A_coeffs(J,I) = AX 
          ENDDO 

    ENDDO 

	DO I=1,NLEVS
 	 	READ(1,*) N, GX, WN !read wavenumbers
        ion%G(N) = GX
        ion%waveno(N) = WN
    enddo
	
    CLOSE(UNIT=1)
	
1002 FORMAT(A20)
1003 FORMAT(78A1)
end subroutine read_atomic_data
	
end module mod_atomic_read
