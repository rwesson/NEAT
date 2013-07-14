module mod_linefinder
contains

subroutine linefinder(linelist, listlength)

use mod_abundtypes
use mod_quicksort

implicit none

TYPE neat_line
        double precision :: wavelength 
        CHARACTER*20 :: ion
END TYPE

TYPE assigned_line
        double precision :: wavelength
        double precision :: flux
        double precision :: uncertainty
        double precision :: diff
end type assigned_line

type spectrum
        double precision :: wavelength
        double precision :: flux
end type spectrum

type(neat_line), dimension(:), allocatable :: neatlines
type(assigned_line), dimension(:), allocatable :: assignedlines
double precision, dimension(10) :: xcorr_reference
double precision, dimension(2001) :: xcorr
integer, intent(in) :: listlength
TYPE(line), dimension(listlength) :: linelist
TYPE(line), dimension(listlength) :: linelist_copy
double precision, dimension(20) :: linelist_compare

integer :: I, J, n_neatlines, IO, assign_1, assign_2, count
double precision :: temp_wave, temp_flux, temp_unc, diff, shift, rms
character*1 :: null
character*5 :: temp_ion1, temp_ion2

xcorr = 0.D0

!algorithm:

!1. read NEAT line list in

        I = 1
        OPEN(100, file='utilities/complete_line_list', iostat=IO, status='old')
                DO WHILE (IO >= 0)
                        READ(100,"(A200)",end=101) null
                        I = I + 1
                END DO
        101 n_neatlines=I-1

!then allocate and read
        allocate (neatlines(n_neatlines))

        REWIND (100)
        DO I=1,n_neatlines
                READ(100,*,end=102) temp_wave, temp_ion1, temp_ion2
                neatlines(i)%wavelength = temp_wave
                neatlines(i)%ion = temp_ion1//temp_ion2 
        END DO
        102 print *,n_neatlines," read from NEAT database"
        CLOSE(100)

!first, find shift

!use these lines for cross correlation - H I, He I, OIII and ArIII as they will almost always be present in spectra covering their wavelengths

xcorr_reference = (/4101.74, 4340.47, 4471.50, 4861.33, 5006.84, 5875.66, 6562.77, 6678.16, 7065.25, 7751.43/)

! find the 20 strongest lines in the observed line list which have wavelengths within the range of reference lines

linelist_copy = linelist

! first, remove lines outside the range of reference lines

do i=1,listlength
  if (linelist_copy(i)%wavelength .lt. minval(xcorr_reference)-10 .or. linelist_copy(i)%wavelength .gt. maxval(xcorr_reference)+10) then
    linelist_copy(i)%intensity = 0.D0
  end if
end do

! then, copy across the 20 strongest lines

do i=1,20
  !copy wavelength of line with highest intensity
  linelist_compare(i)=linelist_copy(maxloc(linelist_copy%intensity,1))%wavelength 
  !replace that line with an intensity of zero so that we can repeat and get the next strongest
  linelist_copy(maxloc(linelist_copy%intensity))%intensity = 0.D0
end do

!sort into wavelength order

call qsort(linelist_compare)

!calculate cross correlation function
!find closest observed line to each rest line for a given shift, calculate sum
!of differences as means of quantifying the match over all lines
xcorr = 0.D0

do i=-1000,1000
  do j=1,10
    xcorr(i+1001)=xcorr(i+1001)+minval(abs(linelist_compare - xcorr_reference(j) - 0.01*dble(i)))
  end do
end do

shift=(minloc(xcorr,1)-1001)*0.01
linelist_compare = linelist_compare + shift

print "(A28,F7.3)", "estimated wavelength shift: ",shift

!now calculate a tolerance from the scatter among closely identified lines
!if difference to a recognised line after shift applied is less than 2, then
!include the difference in calculation of RMS difference
!this value can then be used as a tolerance

rms= 0.D0
count=0

do j=1,10
  if (minval(abs(linelist_compare + shift - xcorr_reference(j))) .lt. 2.0) then
    rms = rms + (minval(abs(linelist%wavelength + shift - xcorr_reference(j))))**2
    count = count + 1
  endif
end do

if (count .gt. 0) then
  rms = (rms/count)**0.5
endif

if (rms .lt. 0.05) then
  rms = 0.05
endif

print *,count," reference lines detected"
print "(A54,F6.3)","RMS difference between observed and rest wavelengths: ",rms

count=0

105 format (F8.2,1X,F15.3,1X,F15.3)

do I=1,listlength
  diff = 1.e10
  do j=1,n_neatlines
    if (abs(linelist(i)%wavelength+shift-neatlines(j)%wavelength).lt.diff) then
      diff = abs(linelist(i)%wavelength+shift-neatlines(j)%wavelength)
      assign_1 = i
      assign_2 = j
    endif 
  enddo
  if (diff .le. rms) then
    print "(F8.2,A4,F8.2,1X,A10,F6.3)",linelist(assign_1)%wavelength," ->  ",neatlines(assign_2), diff
!            write (100,105) neatlines(assign_2)%wavelength, linelist(assign_1)%flux, linelist(assign_1)%uncertainty
    linelist(assign_1)%wavelength = neatlines(assign_2)%wavelength
    count=count+1
  elseif (diff .gt. rms .and. diff .le. 2*rms) then
    print "(F8.2,A4,F8.2,1X,A10,F6.3)",linelist(assign_1)%wavelength," -> ",neatlines(assign_2), diff
!            write (100,105) neatlines(assign_2)%wavelength, linelist(assign_1)%flux, linelist(assign_1)%uncertainty
    linelist(assign_1)%wavelength = neatlines(assign_2)%wavelength
    count=count+1
  elseif (diff .gt. 2*rms .and. diff .le. 3*rms) then
    print "(F8.2,A4,F8.2,1X,A10,F6.3)",linelist(assign_1)%wavelength," -> ",neatlines(assign_2), diff
!            write (100,105) neatlines(assign_2)%wavelength, linelist(assign_1)%flux, linelist(assign_1)%uncertainty
    linelist(assign_1)%wavelength = neatlines(assign_2)%wavelength
    count=count+1
  else
!            write (100,105) linelist(assign_1)%wavelength, linelist(assign_1)%flux, linelist(assign_1)%uncertainty
    print "(F8.2,A30,F8.2,1X,A10,A21)",linelist(assign_1)%wavelength," unrecognised. nearest known: ",neatlines(assign_2)
  endif
enddo

print *,count," lines identified; ",listlength-count," unidentified"
print *,"Wavelengths of identified lines changed as necessary"

end subroutine linefinder
end module mod_linefinder
