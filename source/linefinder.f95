module mod_linefinder
contains

subroutine linefinder(linelist, listlength)

use mod_abundtypes
use mod_quicksort

implicit none

TYPE xcorrarray
        double precision :: restwavelength
        double precision :: observedwavelength
        integer :: match ! 0 if the reference line is not observed, 1 if it is
end TYPE

TYPE neat_line
        double precision :: wavelength 
        CHARACTER(len=20) :: ion
END TYPE

type(neat_line), dimension(:), allocatable :: neatlines
type(xcorrarray), dimension(10) :: xcorr_array
double precision, dimension(2001) :: xcorr
integer, intent(in) :: listlength
TYPE(line), dimension(listlength) :: linelist
TYPE(line), dimension(listlength) :: linelist_copy
double precision, dimension(20) :: linelist_compare

integer :: I, J, n_neatlines, IO, assign_1, assign_2, count
double precision :: temp_wave, diff, shift, rms
character(len=1) :: null
character(len=5) :: temp_ion1, temp_ion2

xcorr = 0.D0
xcorr_array%restwavelength = 0.D0
xcorr_array%observedwavelength = 0.D0
xcorr_array%match = 0.D0

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
        102 print *
        CLOSE(100)

!first, find shift

!use these lines for cross correlation - H I, He I, OIII and ArIII as they will almost always be present in spectra covering their wavelengths

xcorr_array%restwavelength = (/4101.74, 4340.47, 4471.50, 4861.33, 5006.84, 5875.66, 6562.77, 6678.16, 7065.25, 7751.43/)

! find the 20 strongest lines in the observed line list which have wavelengths within the range of reference lines

linelist_copy = linelist

! first, remove lines outside the range of reference lines

do i=1,listlength
  if (linelist_copy(i)%wavelength .lt. minval(xcorr_array%restwavelength)-10 .or. linelist_copy(i)%wavelength .gt. maxval(xcorr_array%restwavelength)+10) then
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

!find nearest observed line within maximum shift range of the reference lines

do i=1,10
  if (minval(abs(linelist_compare - xcorr_array(i)%restwavelength)) .lt. 10) then
    xcorr_array(i)%observedwavelength = linelist_compare(minloc(abs(linelist_compare - xcorr_array(i)%restwavelength),1))
    xcorr_array(i)%match = 1.D0
  endif
enddo

!calculate cross correlation function
!find closest observed line to each rest line for a given shift, calculate sum
!of differences as means of quantifying the match over all lines
xcorr = 0.D0

do i=-1000,1000
  do j=1,10
    xcorr(i+1001)=xcorr(i+1001)+(abs(xcorr_array(j)%observedwavelength - xcorr_array(j)%restwavelength - 0.01*dble(i)) * xcorr_array(j)%match)
  end do
end do

shift=(minloc(xcorr,1)-1001)*0.01
linelist_compare = linelist_compare - shift

!now calculate the rms scatter between rest wavelengths and observed after shift
!this value can then be used as a tolerance when assigning line IDs

rms= 0.D0
count=0

do j=1,10
    rms = rms + (xcorr_array(j)%match*(xcorr_array(j)%observedwavelength - shift - xcorr_array(j)%restwavelength)**2)
    count = count + xcorr_array(j)%match
end do

if (count .gt. 0) then
  rms = (rms/count)**0.5
endif

if (rms .lt. 0.05) then
  rms = 0.05
endif

print "(I2,A25)",count," reference lines detected"
print "(X,A54,F5.3)","Average offset between observed and rest wavelengths: ",shift
print "(X,A74,F6.3)","RMS difference between rest wavelengths and shifted observed wavelengths: ",rms
print *
print *,"The following line IDs are suggested:"
print *
print *,"Obs         rest    ID         offset"
print *,"-------------------------------------"

count=0

do I=1,listlength
  diff = 1.e10
  do j=1,n_neatlines
    if (abs(linelist(i)%wavelength-shift-neatlines(j)%wavelength).lt.diff) then
      diff = abs(linelist(i)%wavelength-shift-neatlines(j)%wavelength)
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
print *,"-------------------------------------"
print *
print *,count," lines identified; ",listlength-count," unidentified"
print *,"Wavelengths of identified lines changed as necessary"

end subroutine linefinder
end module mod_linefinder
