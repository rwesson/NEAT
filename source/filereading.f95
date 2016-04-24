module mod_abundIO
use mod_abundtypes
implicit none
private :: dp
integer, parameter :: dp = kind(1.d0)

contains

subroutine read_linelist(filename,linelist,listlength,ncols,errstat)
        implicit none
        integer :: i, j, io, nlines, listlength, errstat, ncols
        character(len=1) :: blank
        type(line), dimension(:), allocatable :: linelist
        character(len=512) :: filename, rowdata
        CHARACTER(len=15),dimension(4) :: invar !line fluxes and uncertainties are read as strings into these variables
        real(kind=dp),dimension(5) :: rowdata2 !to check number of columns, first row of file is read as character, then read as real into this array

        TYPE neat_line
          real(kind=dp) :: wavelength
          CHARACTER(len=85) :: linedata
        END TYPE

        type(neat_line), dimension(:), allocatable :: neatlines

! first get number of rows

        errstat=0
        I = 0
        OPEN(199, file=filename, iostat=IO, status='old')
          DO WHILE (IO >= 0)
            READ(199,*,end=111) blank
            if (blank.ne."#") then
              I = I + 1
            endif
          END DO
        111 continue
        listlength=I

!then allocate and read
        allocate (linelist(listlength))

        linelist%intensity = 0.D0
        linelist%abundance = 0.D0
        linelist%freq=0d0
        linelist%wavelength=0d0
        linelist%wavelength_observed=0d0
        linelist%int_dered=0d0
        linelist%int_err=0d0
        linelist%blend_intensity=0.d0
        linelist%blend_int_err=0d0
        linelist%zone='    '
        linelist%name='           '
        linelist%transition='                    '
        linelist%location=0
        linelist%ion='                   '
        linelist%latextext='               '
        linelist%linedata='                                                                           '

! now count the columns
! if 2 - assume rest wavelength and intensity, read in and restrict to single iteration
! if 3 - assume rest wavelength, intensity, uncertainty, read in as we do currently
! if 4 - assume observed wavelength, rest wavelength, intensity, uncertainty, read in and add obs wlen to output table

        rewind(199)
        rowdata2=-1.0d-27

        do while (io>=0)
          read (199,"(A)") rowdata
          if (index(rowdata,"#") .ne. 1) then
            read (rowdata,*,end=113) rowdata2(:)
          endif
        enddo

113     ncols=count(rowdata2 .ne. -1.0d-27)
        invar="               "

        REWIND (199)
        I=1
        DO while (I.le.listlength)
          read(199,"(A)",end=110) rowdata

          if (index(rowdata,"#") .ne. 1) then !not a comment, read in the columns
            READ(rowdata,*) invar(1:ncols)
          else !do nothing with comment lines
            cycle
          endif

          if (ncols .eq. 4) then
            read (invar(2),*) linelist(i)%wavelength
          else
            read (invar(1),*) linelist(i)%wavelength
          endif
          if (index(invar(2),"*") .gt. 0 .or. index(invar(3),"*") .gt. 0) then
!line is blended, its intensity will be removed from abundance and diagnostic calculations but retained in linelist
            linelist(i-1)%blend_intensity=linelist(i-1)%intensity
            linelist(i-1)%blend_int_err=linelist(i-1)%int_err
            linelist(i-1:i)%intensity = 0.d0
            linelist(i-1:i)%int_err = 0.d0
          else
            if (ncols .eq. 2) then
              read (invar(2),*) linelist(i)%intensity
            elseif (ncols .eq. 3) then
              read (invar(2),*) linelist(i)%intensity
              read (invar(3),*) linelist(i)%int_err
            elseif (ncols .eq. 4) then
              read (invar(1),*) linelist(i)%wavelength_observed
              read (invar(3),*) linelist(i)%intensity
              read (invar(4),*) linelist(i)%int_err
            endif
          endif
          linelist(i)%latextext = ""
          i=i+1
        END DO

        110 continue
        CLOSE(199)

! check for errors

        if (I - 1 .ne. listlength) then
          errstat=errstat+1
          return
        endif

        if (linelist(1)%wavelength == 0) then
          errstat=errstat+2
          return
        endif

        if (ncols .eq. 2) then !set uncertainties to 10 per cent, warn
                linelist%int_err=linelist%intensity*0.1
                errstat=errstat+4
        endif

!if no fatal errors, proceed to copying line data into the array

        I = 1
        OPEN(100, file=trim(PREFIX)//'/share/neat/complete_line_list', iostat=IO, status='old')
          DO WHILE (IO >= 0)
            READ(100,"(A1)",end=101) blank
            I = I + 1
          END DO
        101 nlines=I-1

!then allocate and read
        allocate (neatlines(nlines))

        REWIND (100)
        DO I=1,nlines
          READ(100,"(F8.2,A85)",end=102) neatlines(i)%wavelength,neatlines(i)%linedata
          do j=1,listlength
            if (abs(linelist(j)%wavelength - neatlines(i)%wavelength) .lt. 0.011) then
              linelist(j)%linedata = neatlines(i)%linedata
            endif
          enddo
        END DO
        102 print *
        CLOSE(100)

end subroutine read_linelist

subroutine read_celdata(ILs, ionlist)
!this subroutine reads in the parameters of CEL diagnostic lines - the line id, ion name, wavelength, energy levels, and zone.  the index locating the line in the main line list is filled in later.
!the routine also populates the list of ions which is then used to read in the relevant atomic data

        implicit none
        type(cel), dimension(82) :: ILs !todo: restore this to allocatable once I've worked out why assumed shape caused it to become undefined on entry to abundances
        character(len=20), dimension(22) :: ionlist !todo: find a clever way of counting this instead of hard coding it if possible.
        integer :: i, iint, iion, numberoflines

        Iint = 1

        301 format(A11, 1X, A6, 1X, F7.2, 1X, A20,1X,A4,1X,A15)
        open(201, file=trim(PREFIX)//"/share/neat/Ilines_levs", status='old')
                read (201,*) numberoflines !number of lines to read in is given at the top of the file
!                allocate (ILs(numberoflines)) todo: restore
                Ils%name='           '
                Ils%ion='                   '
                Ils%wavelength=0.d0
                Ils%transition='                    '
                Ils%zone='    '
                Ils%latextext='               '
                Ils%location=0
                do while (Iint .le. numberoflines)
                        read(201,301) ILs(Iint)%name, ILs(Iint)%ion, ILs(Iint)%wavelength, ILs(Iint)%transition ,ILs(Iint)%zone, ILs(Iint)%latextext!end condition breaks loop.
                        if(Iint .eq. 1) then
                            Iion = 1
                            Ionlist(iion) = ILs(Iint)%ion
                        elseif(ILs(Iint)%ion .ne. ILs(Iint - 1)%ion) then
                            Iion = iion + 1
                            Ionlist(iion) = ILs(Iint)%ion
                        endif
                        Iint = Iint + 1
                END DO
                Iint = Iint - 1 !count ends up one too high
        CLOSE(201)

end subroutine read_celdata

end module

module mod_abundmaths
use mod_abundtypes
use mod_atomicdata
implicit none
integer, parameter :: dp=kind(1.d0)

contains

!this fantastically ugly function gets the location of certain ions in the important ions array using their name as a key.
!April 2016: made it yet uglier so that it:
! - gets the ion's location in the CELs array
! - if the line is in the CEL array, then look up its dereddened flux in the main linelist array
!todo: check if we ever need the undereddened flux

real(kind=dp) function get_cel_flux(ionname, linelist, ILs)
        implicit none
        character(len=11) :: ionname
        type(cel), dimension(:) :: ILs
        type(line), dimension(:) :: linelist
        integer :: i

        get_cel_flux=0.d0

        do i = 1, size(ILs)
          if(trim(ILs(i)%name) == trim(ionname))then
            if (ILs(i)%location .gt. 0) then
              get_cel_flux=linelist(ILs(i)%location)%int_dered
            else
              get_cel_flux=0.d0
            endif
            return
          endif
        end do

end function get_cel_flux

real(kind=dp) function get_cel_abundance(ionname, linelist, ILs)
!another ugly function to get an abundance for a cel
        implicit none
        character(len=11) :: ionname
        type(cel), dimension(:) :: ILs
        type(line), dimension(:) :: linelist
        integer :: i

        get_cel_abundance=0.d0

        do i = 1, size(ILs)
          if(trim(ILs(i)%name) == trim(ionname))then
            if (ILs(i)%location .gt. 0) then
              get_cel_abundance=linelist(ILs(i)%location)%abundance
            else
              get_cel_abundance=0.d0
            endif
            return
          endif
        end do

end function get_cel_abundance

!this fantastically ugly function gets the location of certain ions in the important ions array using their name as a key.

integer function get_ion(ionname, ILs)
        implicit none
        character(len=11) :: ionname
        type(cel), dimension(:) :: ILs
        integer :: i

        do i = 1, size(ILs)
          if(trim(ILs(i)%name) == trim(ionname))then
            get_ion = i
            return
          endif
        end do

        get_ion = 0
        print*, "Nudge Nudge, wink, wink error. Ion not found, say no more.", ionname

end function

!same as above for getting the location of ion within atomic data array. equally ugly.

integer function get_atomicdata(ionname, atomicdatatable)
        implicit none
        character(len=20) :: ionname
        type(atomic_data), dimension(:) :: atomicdatatable
        integer :: i

        do i = 1, size(atomicdatatable)
          if(trim(atomicdatatable(i)%ion) == trim(ionname))then
            get_atomicdata = i
            return
          endif
        end do

        get_atomicdata = 0
        print*, "My hovercraft is full of eels.  Atomic data not found. ", ionname

end function

subroutine get_cels(ILs, linelist)
!index the locations of CELs within the main line list
        implicit none
        type(cel), dimension(:) :: ILs
        type(line), dimension(:) :: linelist
        integer :: i, j
        character(len=11) :: temp

        ILs%location=0

        do i = 1, size(ILs)
          do j = 1, size(linelist)
            if(linelist(j)%wavelength .eq. ILs(i)%wavelength)then
              ILs(i)%location = j
              cycle
            endif
          end do
        end do

end subroutine

subroutine get_H(H_Balmer, H_paschen, linelist)
!index the location of Balmer and Paschen lines in the main array
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer, dimension(3:40), intent(out) :: H_balmer ! indexing starts at 3 so that it represents the upper level of the line
  integer, dimension(4:39), intent(out) :: H_paschen ! indexing starts at 4 for the same reason
  real(kind=dp), dimension(3:40) :: balmerwavelengths
  real(kind=dp), dimension(4:39) :: paschenwavelengths
  type(line), dimension(:) :: linelist
  integer :: i,j

  balmerwavelengths = (/ 6562.77D0, 4861.33D0, 4340.47D0, 4101.74D0, 3970.07D0, 3889.05D0, 3835.38D0, 3797.90D0, 3770.63D0, 3750.15D0, 3734.37D0, 3721.94D0, 3711.97D0, 3703.85D0, 3697.15D0, 3691.55D0, 3686.83D0, 3682.81D0, 3679.35D0, 3676.36D0, 3673.76D0, 3671.48D0, 3669.46D0, 3667.68D0, 3666.10D0, 3664.68D0, 3663.40D0, 3662.26D0, 3661.22D0, 3660.28D0, 3659.42D0, 3658.64D0, 3657.92D0, 3657.27D0, 3656.66D0, 3656.11D0, 3655.59D0, 3655.12D0 /)

  paschenwavelengths = (/ 18751.01d0, 12818.08d0, 10938.10d0, 10049.37d0, 9545.97d0, 9229.01d0, 9014.91d0, 8862.78d0, 8750.47d0, 8665.02d0, 8598.39d0, 8545.38d0, 8502.48d0, 8467.25d0, 8437.95d0, 8413.32d0, 8392.40d0, 8374.48d0, 8359.00d0, 8345.47d0, 8333.78d0, 8323.42d0, 8314.26d0, 8306.11d0, 8298.83d0, 8292.31d0, 8286.43d0, 8281.12d0, 8276.31d0, 8271.93d0, 8267.94d0, 8264.28d0, 8260.93d0, 8255.02d0, 8252.40d0, 8249.97d0 /)

  H_balmer = 0
  H_paschen = 0

  do i=1,size(linelist)
    do j=3,40
      if (abs(linelist(i)%wavelength - balmerwavelengths(j)) .lt. 0.005) then
        H_balmer(j)=i
        cycle
      endif
    enddo
  enddo

  do i=1,size(linelist)
    do j=4,39
      if (abs(linelist(i)%wavelength - paschenwavelengths(j)) .lt. 0.005) then
        H_paschen(j)=i
        cycle
      endif
    enddo
  enddo

end subroutine

subroutine get_HeI(HeI_lines, linelist)
!index the locations of He I lines in the main linelist array
        implicit none
        integer, parameter :: dp = kind(1.d0)
        integer, dimension(44), intent(out) :: HeI_lines
        type(line), dimension(:), intent(in) :: linelist
        real(kind=dp), dimension(44) :: wavelengths
        integer :: i, j


        wavelengths = (/ 2945.10D0,3188.74D0,3613.64D0,3888.65D0,3964.73D0,4026.21D0,4120.82D0,4387.93D0,4437.55D0,4471.50D0,4713.17D0,4921.93D0,5015.68D0,5047.74D0,5875.66D0,6678.16D0,7065.25D0,7281.35D0,9463.58D0,10830.25D0,11013.07D0,11969.06D0,12527.49D0,12755.69D0,12784.92D0,12790.50D0,12845.98D0,12968.43D0,12984.88D0,13411.69D0,15083.65D0,17002.40D0,18555.57D0,18685.33D0,18697.21D0,19089.36D0,19543.19D0,20424.97D0,20581.28D0,20601.76D0,21120.12D0,21132.03D0,21607.80D0,21617.01D0 /)

        HeI_lines = 0

        do i = 1, 44
          do j = 1, size(linelist)
            if(abs(linelist(j)%wavelength - wavelengths(i)) .lt.  0.005) then
              Hei_lines(i) = j
              cycle
            endif
          end do
        end do

end subroutine

subroutine get_HeII(HeII_lines, linelist)
!index the location of He II lines (only 4686 at the moment) in the linelist array
        implicit none
        integer, parameter :: dp = kind(1.d0)
        integer, dimension(1), intent(out) :: HeII_lines
        type(line), dimension(:), intent(in) :: linelist
        real(kind=dp), dimension(1) :: wavelengths
        integer :: i, j

        wavelengths = (/ 4685.68D0 /)

        heii_lines = 0

        do i = 1, 1
          do j = 1, size(linelist)
            if(abs(linelist(j)%wavelength-wavelengths(i)).lt.0.005) then
              Heii_lines(i) = j
              cycle
            endif
          end do
        end do

end subroutine

!extinction laws now in mod_extinction

end module mod_abundmaths

module mod_atomic_read
use mod_atomicdata

private :: dp
integer, parameter :: dp = kind(1.d0)

contains
subroutine read_atomic_data(ion)
use mod_atomicdata
    implicit none
    type(atomic_data) :: ion
    integer :: I,J,K,L,N,NCOMS,ID(2),JD(2),KP1,NLEV1,GX,ionl,dummy
    character(len=1) :: comments(78)
    character(len=10) :: ionname
    character(len=128) :: filename
    real(kind=dp) :: WN,AX,QX

    id = 0
    jd = 0
    ionname = ion%ion
!    print*,'Reading atomic data ion',ionname
    ionl = index(ionname,' ') - 1
    filename = trim(PREFIX)//'/share/neat/'//ionname(1:IONL)//'.dat'
    open(unit=1, status = 'old', file=filename,action='read')

!read # of comment lines and skip them
        read(1,*)NCOMS
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
        ion%G = 0
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
        DO while (QX .gt. 0)
                read(1,*) ID(2), JD(2), QX
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
            read (1,*) I, J, AX  !read transition probabilities
            ion%A_coeffs(J,I) = AX
          enddo

    enddo

    DO I=1,ion%NLEVS
          read(1,*) N, GX, WN !read wavenumbers
        ion%G(N) = GX
        ion%waveno(N) = WN
    enddo

    close(unit=1)

1002 format(A20)
1003 format(78A1)
end subroutine read_atomic_data

!read in tables of helium emissivities from Porter et al.
!http://cdsads.u-strasbg.fr/abs/2012MNRAS.425L..28P

subroutine read_porter(heidata)

implicit none
real(kind=dp), dimension(21,15,44) :: heidata
integer :: i,j,tpos,npos,io
real(kind=dp), dimension(46) :: temp

!read data

open(100, file=trim(PREFIX)//'/share/neat/RHei_porter2012.dat', iostat=IO, status='old')

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
real(kind=dp), dimension(3,6,44) :: heidata
integer :: i,j,k,io
real(kind=dp), dimension(18) :: temp

!read data
!fitted fourth order polynomials to the Smits 1996 emissivities
!the data file contains the coefficients for log(ne)=2,4,6

open(100, file=trim(PREFIX)//'/share/neat/RHei_smits1996_coeffs.dat', iostat=IO, status='old')

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
