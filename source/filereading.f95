module mod_abundIO
use mod_abundtypes
implicit none!

contains

subroutine read_ilines(ILs, Iint,iion,ionlist)
        IMPLICIT NONE
        TYPE(line), DIMENSION(:), allocatable :: ILs
        INTEGER :: Iint, Iread,iion
        character(len=10) :: ionlist(40)

        Iint = 1

        301 FORMAT(A11, 1X, A6, 1X, F7.2, 1X, A20,1X,A4,1X,A15)
        OPEN(201, file="source/Ilines_levs", status='old')
                READ (201,*) Iread
                ALLOCATE (ILs(Iread))
                ILs%intensity=0.D0 !otherwise it seems you can get random very small numbers in the array.
                Ils%wavelength=0d0
                ils%int_err=0d0
                ils%abundance=0d0
                ils%freq=0d0
                ils%int_dered=0d0
                ils%zone='    '
                ils%name='           '
                ils%transition='                    '
                ils%location=0
                ils%ion='                   '
                ils%latextext='               '
                DO WHILE (Iint .le. Iread)!(.true.)
                        READ(201,301) ILs(Iint)%name, ILs(Iint)%ion, ILs(Iint)%wavelength, ILs(Iint)%transition ,ILs(Iint)%zone, ILs(Iint)%latextext!end condition breaks loop.
                        if(Iint .eq. 1) then
                            Iion = 1
                            Ionlist(iion) = ILs(Iint)%ion(1:10)
                        elseif(ILs(Iint)%ion .ne. ILs(Iint - 1)%ion) then
                            Iion = iion + 1
                            Ionlist(iion) = ILs(Iint)%ion(1:10)
                        endif
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
        IMPLICIT NONE
        CHARACTER(len=11) :: ionname
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

        get_ion = 0
        PRINT*, "Nudge Nudge, wink, wink error. Ion not found, say no more.", ionname

end function


subroutine element_assign(ILs, linelist, Iint, listlength)
        IMPLICIT NONE
        TYPE(line), DIMENSION(:), INTENT(OUT) :: ILs
        TYPE(line), DIMENSION(:) :: linelist
        INTEGER, INTENT(IN) :: Iint, listlength
        INTEGER :: i, j

        do i = 1, Iint
                ILs(i)%location = 0 !initialise first to prevent random integers appearing and breaking things
                do j = 1, listlength
                        if(linelist(j)%wavelength == ILs(i)%wavelength)then
                                ILs(i)%intensity = linelist(j)%intensity
                                ILs(i)%int_err   = linelist(j)%int_err
                                ILs(i)%location = j !recording where the line is in linelist array so that its abundance can be copied back in
                                cycle
                        endif
                end do
        end do

end subroutine

subroutine get_H(H_BS, linelist, listlength)
        IMPLICIT NONE
        TYPE(line), DIMENSION(38), INTENT(OUT) :: H_BS
        TYPE(line), DIMENSION(:) :: linelist
        double precision, dimension(38) :: balmerlines
        INTEGER :: i, j, listlength
        REAL*8 :: HW = 0.00000000
        !another ugly kludge, but it works.

        balmerlines = (/ 6562.77D0, 4861.33D0, 4340.47D0, 4101.74D0, 3970.07D0, 3889.05D0, 3835.38D0, 3797.90D0, 3770.63D0, 3750.15D0, 3734.37D0, 3721.94D0, 3711.97D0, 3703.85D0, 3697.15D0, 3691.55D0, 3686.83D0, 3682.81D0, 3679.35D0, 3676.36D0, 3673.76D0, 3671.48D0, 3669.46D0, 3667.68D0, 3666.10D0, 3664.68D0, 3663.40D0, 3662.26D0, 3661.22D0, 3660.28D0, 3659.42D0, 3658.64D0, 3657.92D0, 3657.27D0, 3656.66D0, 3656.11D0, 3655.59D0, 3655.12D0 /)
        H_BS%location = 0.0
        H_Bs%intensity=0.D0 !otherwise it seems you can get random very small numbers in the array.
        H_Bs%wavelength=0d0
        H_Bs%int_err=0d0
        H_Bs%abundance=0d0
        H_Bs%freq=0d0
        H_Bs%int_dered=0d0
        H_Bs%zone='    '
        H_Bs%name='           '
        H_Bs%transition='                    '
        H_Bs%ion='                   '
        H_Bs%latextext='               '

        do i = 1, 38
                HW = balmerlines(i)

                do j = 1, listlength
                         if (linelist(j)%wavelength-HW==0) then
                                H_BS(i)%wavelength = linelist(j)%wavelength
                                H_BS(i)%intensity = linelist(j)%intensity
                                H_BS(i)%int_err = linelist(j)%int_err
                                H_BS(i)%location = j
                        endif
                end do
        end do

end subroutine

subroutine get_HeI(HeI_lines, linelist,listlength)
        IMPLICIT NONE
        TYPE(line), DIMENSION(44), INTENT(OUT) :: HeI_lines
        TYPE(line), DIMENSION(:), INTENT(IN) :: linelist
        INTEGER :: i, j, listlength
        double precision, dimension(44) :: hei_wavelengths

        hei_wavelengths = (/ 2945.10D0,3188.74D0,3613.64D0,3888.65D0,3964.73D0,4026.21D0,4120.82D0,4387.93D0,4437.55D0,4471.50D0,4713.17D0,4921.93D0,5015.68D0,5047.74D0,5875.66D0,6678.16D0,7065.25D0,7281.35D0,9463.58D0,10830.25D0,11013.07D0,11969.06D0,12527.49D0,12755.69D0,12784.92D0,12790.50D0,12845.98D0,12968.43D0,12984.88D0,13411.69D0,15083.65D0,17002.40D0,18555.57D0,18685.33D0,18697.21D0,19089.36D0,19543.19D0,20424.97D0,20581.28D0,20601.76D0,21120.12D0,21132.03D0,21607.80D0,21617.01D0 /)
        HeI_lines%location = 0 !initialise to prevent bugs!
        HeI_lines%intensity=0.D0 !otherwise it seems you can get random very small numbers in the array.
        HeI_lines%wavelength=0d0
        HeI_lines%int_err=0d0
        HeI_lines%abundance=0d0
        HeI_lines%freq=0d0
        HeI_lines%int_dered=0d0
        HeI_lines%zone='    '
        HeI_lines%name='           '
        HeI_lines%transition='                    '
        HeI_lines%ion='                   '
        HeI_lines%latextext='               '

        do i = 1, 44
                HeI_lines(i)%wavelength = hei_wavelengths(i)
                do j = 1, listlength
                        if(linelist(j)%wavelength == HeI_lines(i)%wavelength) then
                                Hei_lines(i)%wavelength = linelist(j)%wavelength
                                Hei_lines(i)%intensity = linelist(j)%intensity
                                Hei_lines(i)%int_err = linelist(j)%int_err
                                Hei_lines(i)%location = j
                        endif
                end do
        end do

end subroutine

subroutine get_HeII(HeII_lines, linelist,listlength)
        IMPLICIT NONE
        TYPE(line), DIMENSION(1), INTENT(OUT) :: HeII_lines
        TYPE(line), DIMENSION(:), INTENT(IN) :: linelist
        INTEGER :: i, j, listlength
        REAL*8 :: HW
        double precision, dimension(1) :: heii_wavelengths

        heii_wavelengths = (/ 4685.68D0 /)
        heii_lines%location = 0
        HeIi_lines%intensity=0.D0 !otherwise it seems you can get random very small numbers in the array.
        HeIi_lines%wavelength=0d0
        HeIi_lines%int_err=0d0
        HeIi_lines%abundance=0d0
        HeIi_lines%freq=0d0
        HeIi_lines%int_dered=0d0
        HeIi_lines%zone='    '
        HeIi_lines%name='           '
        HeIi_lines%transition='                    '
        HeIi_lines%ion='                   '
        HeIi_lines%latextext='               '

        do i = 1, 1
                HW = heii_wavelengths(i)
                do j = 1, listlength
                        if(linelist(j)%wavelength == HW) then
                                Heii_lines(i)%wavelength = linelist(j)%wavelength
                                Heii_lines(i)%intensity = linelist(j)%intensity
                                Heii_lines(i)%int_err = linelist(j)%int_err
                                Heii_lines(i)%location = j
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
    IMPLICIT NONE
    type(atomic_data) :: ion
    integer :: I,J,K,L,N,NCOMS,ID(2),JD(2),KP1,NLEV1,ionl,dummy
    character(len=1) :: comments(78)
    character(len=10) :: ionname
    character(len=25) :: filename
    real*8 :: GX,WN,AX,QX

    id = 0
    jd = 0
    ionname = ion%ion
!    print*,'Reading atomic data ion',ionname
    ionl = index(ionname,' ') - 1
    filename = 'Atomic-data/'//ionname(1:IONL)//'.dat'
    OPEN(unit=1, status = 'OLD', file=filename,ACTION='READ')

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
        ion%temps=0d0
        ion%roott=0d0

        !read levels and temperatures
        do I = 1,ion%NLEVS
        read(1,1002) ion%labels(I)
        enddo

        do I = 1,ion%NTEMPS
        read(1,*) ion%temps(I)
        enddo

        read(1,*) dummy

        !read collision strengths
        QX=1
        K = 1
!        print*,'Reading collision strengths'
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
!                print*,k,i,j
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

    DO I=1,ion%NLEVS
          READ(1,*) N, GX, WN !read wavenumbers
        ion%G(N) = GX
        ion%waveno(N) = WN
    enddo

    CLOSE(UNIT=1)

1002 FORMAT(A20)
1003 FORMAT(78A1)
end subroutine read_atomic_data

!read in tables of helium emissivities from Porter et al.
!http://cdsads.u-strasbg.fr/abs/2012MNRAS.425L..28P

subroutine read_porter(heidata)

implicit none
double precision, dimension(21,15,44) :: heidata
integer :: i,j,tpos,npos,io
double precision, dimension(46) :: temp

!read data

OPEN(100, file='Atomic-data/RHei_porter2012.dat', iostat=IO, status='old')

! read in the data

do i=1,294
  read (100,*) temp
  tpos=nint((temp(1)/1000)-4)
  npos=nint(temp(2))
  do j=1,44
    heidata(tpos,npos,j)=temp(j+2)
  end do
end do

close(100)

end subroutine

subroutine read_smits(heidata)

implicit none
double precision, dimension(3,6,44) :: heidata
integer :: i,j,k,io
double precision, dimension(18) :: temp

!read data
!fitted fourth order polynomials to the Smits 1996 emissivities
!the data file contains the coefficients for log(ne)=2,4,6

OPEN(100, file='Atomic-data/RHei_smits1996_coeffs.dat', iostat=IO, status='old')

! read in the data

do i=1,44
  read (100,*) temp
  do j=1,3
    do k=1,6
      heidata(j,k,i)=temp(k+((j-1)*6))
    end do
  end do
end do

close(100)

end subroutine


end module mod_atomic_read
