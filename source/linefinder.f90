!linefinder.f90, a rough and definitely non-rigorous way to suggest line identifications
!(C) Roger Wesson, Dave Stock, Peter Scicluna
module mod_linefinder
use mod_globals
use mod_functions

contains

subroutine linefinder(linelist, listlength, identifyconfirm)

use mod_types
use mod_quicksort

implicit none

type xcorrarray
        real(kind=dp) :: restwavelength
        real(kind=dp) :: observedwavelength
        integer :: match ! 0 if the reference line is not observed, 1 if it is
end type xcorrarray

type(neat_line), dimension(:), allocatable :: neatlines
type(xcorrarray), dimension(10) :: xcorr_array
real(kind=dp), dimension(2001) :: xcorr
integer, intent(in) :: listlength
type(line), dimension(listlength) :: linelist
type(line), dimension(listlength) :: linelist_copy
real(kind=dp), dimension(20) :: linelist_compare

integer :: I, J, n_neatlines, IO, assign_1, assign_2, count
real(kind=dp) :: diff, shift, rms
character(len=1) :: readchar
logical :: identifyconfirm

!debugging
#ifdef CO
        print *,"subroutine: linefinder"
#endif

print *,gettime()," : running line finder"
print *,gettime(),"---------------------------------"

xcorr = 0.D0
xcorr_array%restwavelength = 0.D0
xcorr_array%observedwavelength = 0.D0
xcorr_array%match = 0

!algorithm:

!1. read NEAT line list in

        I = 1
        open(100, file=trim(PREFIX)//'/share/neat/complete_line_list', iostat=IO, status='old')
                do while (IO .ge. 0)
                        read(100,"(A200)",end=101) readchar
                        I = I + 1
                enddo
        101 n_neatlines=I-1

!then allocate and read
        allocate (neatlines(n_neatlines))

        rewind (100)
        do I=1,n_neatlines
                read(100,"(F8.2,2X,A11,X,A12,X,A12,X,A12,X,I12,X,I9)",end=102) neatlines(i)%wavelength,neatlines(i)%ion,neatlines(i)%multiplet,neatlines(i)%lowerterm,neatlines(i)%upperterm,neatlines(i)%g1,neatlines(i)%g2
        enddo
        102 print *
        close(100)

!first, find shift

!use these lines for cross correlation - H I, He I, OIII and ArIII as they will almost always be present in spectra covering their wavelengths

xcorr_array%restwavelength = (/4101.74, 4340.47, 4471.50, 4861.33, 5006.84, 5875.66, 6562.77, 6678.16, 7065.25, 7751.43/)

! find the 20 strongest lines in the observed line list which have wavelengths within the range of reference lines

linelist_copy = linelist

! first, remove lines outside the range of reference lines

do i=1,listlength
  if (linelist_copy(i)%wavelength .lt. minval(xcorr_array%restwavelength)-10 .or. linelist_copy(i)%wavelength .gt. maxval(xcorr_array%restwavelength)+10) then
    linelist_copy(i)%intensity = 0.D0
  endif
enddo

! then, copy across the 20 strongest lines

do i=1,20
  !copy wavelength of line with highest intensity
  linelist_compare(i)=linelist_copy(maxloc(linelist_copy%intensity,1))%wavelength
  !replace that line with an intensity of zero so that we can repeat and get the next strongest
  linelist_copy(maxloc(linelist_copy%intensity))%intensity = 0.D0
enddo

!sort into wavelength order

call qsort(linelist_compare)

!find nearest observed line within maximum shift range of the reference lines

do i=1,10
  if (minval(abs(linelist_compare - xcorr_array(i)%restwavelength)) .lt. 10) then
    xcorr_array(i)%observedwavelength = linelist_compare(minloc(abs(linelist_compare - xcorr_array(i)%restwavelength),1))
    xcorr_array(i)%match = 1
  endif
enddo

!calculate cross correlation function
!find closest observed line to each rest line for a given shift, calculate sum
!of differences as means of quantifying the match over all lines
xcorr = 0.D0

do i=-1000,1000
  do j=1,10
    xcorr(i+1001)=xcorr(i+1001)+(abs(xcorr_array(j)%observedwavelength - xcorr_array(j)%restwavelength - 0.01*dble(i)) * dble(xcorr_array(j)%match))
  enddo
enddo

shift=(minloc(xcorr,1)-1001)*0.01
linelist_compare = linelist_compare - shift

!now calculate the rms scatter between rest wavelengths and observed after shift
!this value can then be used as a tolerance when assigning line IDs

rms= 0.D0
count=0

do j=1,10
    rms = rms + (xcorr_array(j)%match*(xcorr_array(j)%observedwavelength - shift - xcorr_array(j)%restwavelength)**2)
    count = count + xcorr_array(j)%match
enddo

if (count .gt. 0) then
  rms = (rms/count)**0.5
endif

if (rms .lt. 0.05) then
  rms = 0.05
endif

print "(I2,A25)",count," reference lines detected"
print "(A55,F5.3)"," Average offset between observed and rest wavelengths: ",shift
print "(A75,F6.3)"," RMS difference between rest wavelengths and shifted observed wavelengths: ",rms
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
    print "(F8.2,A4,F8.2,1X,A10,F6.3)",linelist(assign_1)%wavelength," ->  ",neatlines(assign_2)%wavelength,neatlines(assign_2)%ion,diff
!            write (100,105) neatlines(assign_2)%wavelength, linelist(assign_1)%flux, linelist(assign_1)%uncertainty
    linelist(assign_1)%wavelength = neatlines(assign_2)%wavelength
    count=count+1
  elseif (diff .gt. rms .and. diff .le. 2*rms) then
    print "(F8.2,A4,F8.2,1X,A10,F6.3)",linelist(assign_1)%wavelength," ->  ",neatlines(assign_2)%wavelength,neatlines(assign_2)%ion,diff
!            write (100,105) neatlines(assign_2)%wavelength, linelist(assign_1)%flux, linelist(assign_1)%uncertainty
    linelist(assign_1)%wavelength = neatlines(assign_2)%wavelength
    count=count+1
  elseif (diff .gt. 2*rms .and. diff .le. 3*rms) then
    print "(F8.2,A4,F8.2,1X,A10,F6.3)",linelist(assign_1)%wavelength," ->  ",neatlines(assign_2)%wavelength,neatlines(assign_2)%ion,diff
!            write (100,105) neatlines(assign_2)%wavelength, linelist(assign_1)%flux, linelist(assign_1)%uncertainty
    linelist(assign_1)%wavelength = neatlines(assign_2)%wavelength
    count=count+1
  else
!            write (100,105) linelist(assign_1)%wavelength, linelist(assign_1)%flux, linelist(assign_1)%uncertainty
    print "(F8.2,A30,F8.2,1X,A10,A11)",linelist(assign_1)%wavelength," unrecognised. nearest known: ",neatlines(assign_2)%wavelength,neatlines(assign_2)%ion
  endif
enddo
print *,"-------------------------------------"
print *
print *,gettime(),count," lines identified; ",listlength-count," unidentified"
print *,gettime()," : Wavelengths of identified lines changed as necessary"

print *
print *,gettime()," : line finder finished"
print *,gettime()," : WARNING!!!  The line finding algorithm is intended as an aid only and is not designed to be highly robust"
print *,gettime()," : check your line list very carefully for potentially wrongly identified lines!!"
print *,gettime(),"---------------------------------"
if (.not. identifyconfirm) then
  print *,gettime()," : Are these line IDs ok? (y/n)"
  read (5,*) readchar
  if (readchar .ne. "y" .and. readchar .ne. "Y") then
    print *,gettime()," : analysis cancelled."
    call exit(1)
  endif
else
  print *,gettime()," : line finder finished, assignments automatically accepted"
endif

end subroutine linefinder
end module mod_linefinder
