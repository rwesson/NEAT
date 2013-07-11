program linefinder

implicit none

TYPE neat_line
        double precision :: wavelength 
        CHARACTER*20 :: ion
END TYPE

TYPE user_line
        double precision :: wavelength
        double precision :: flux
        double precision :: uncertainty
end TYPE

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
type(user_line), dimension(:), allocatable :: userlines
type(assigned_line), dimension(:), allocatable :: assignedlines
type(spectrum), dimension(10000) :: refspec, obsspec
double precision, dimension(10) :: xcorr_lines
double precision, dimension(2001) :: xcorr

integer :: I, J, n_neatlines, n_userlines, IO, assign_1, assign_2, count
double precision :: temp_wave, temp_flux, temp_unc, diff, shift, rms
character*1 :: null
character*5 :: temp_ion1, temp_ion2
character*1024 :: filename

obsspec%flux = 0.D0
refspec%flux = 0.D0
xcorr = 0.D0

I = IARGC()
if (I .ne. 1) then
  print *,"No line list specified"
  stop
endif

CALL getarg(1,filename)

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


!2. read observed line list in

        I = 1
        OPEN(200, file=filename, iostat=IO, status='old')
                DO WHILE (IO >= 0)
                        READ(200,"(A200)",end=201) null 
                        I = I + 1
                END DO 
        201 n_userlines=I-1

!then allocate and read
        allocate (userlines(n_userlines))
        allocate (assignedlines(n_userlines))

        REWIND (200)
        DO I=1,n_userlines
                READ(200,*,end=202) temp_wave, temp_flux, temp_unc
                userlines(i)%wavelength = temp_wave
                userlines(i)%flux = temp_flux
                userlines(i)%uncertainty = temp_unc
        END DO
        202 print *,n_userlines," lines read from user linelist"
        CLOSE(200)

!first, find shift

!use these lines for cross correlation - H I, He I, OIII and ArIII as they will almost always be present in spectra covering their wavelengths

xcorr_lines = (/4101.74, 4340.47, 4471.50, 4861.33, 5006.84, 5875.66, 6562.77, 6678.16, 7065.25, 7751.43/)

!calculate cross correlation function
!find closest observed line to each rest line for a given shift, calculate sum
!of differences as means of quantifying the match over all lines
xcorr = 0.D0

do i=-1000,1000
  do j=1,10
    !only consider lines where there is a chance of an observed line lying
    !within the maximum shift:
    if (xcorr_lines(j) .gt. userlines(1)%wavelength - 10 .and. xcorr_lines(j) .lt. userlines(n_userlines)%wavelength + 10) then
      xcorr(i+1001)=xcorr(i+1001)+minval(abs(userlines%wavelength - xcorr_lines(j) - 0.01*dble(i)))
    endif
  end do
end do

shift=(minloc(xcorr,1)-1001)*0.01
userlines%wavelength = userlines%wavelength - shift

print *, "estimated wavelength shift: ",shift

!now calculate a tolerance from the scatter among closely identified lines
!if difference to a recognised line after shift applied is less than 1, then
!include the difference in calculation of RMS difference
!this value can then be used as a tolerance

rms= 0.D0
count=0

do j=1,10
  if (abs(minval(userlines%wavelength - xcorr_lines(j))) .lt. 10.0) then
    rms = rms + (minval(userlines%wavelength - xcorr_lines(j)))**2
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
print *,"RMS difference between observed and rest wavelengths: ",rms

!find the three closest matches in wavelength to given
!if diff = 0 then assign
!if diff < 0.1 then assign and flag
!if diff > 0.1 then don't assign

        count=0
        open (100, file=trim(filename)//".proc", status='replace', action='write')
        105 format (F8.2,1X,F15.3,1X,F15.3)

        do I=1,n_userlines
          diff = 1.e10
          do j=1,n_neatlines
            if (abs(userlines(i)%wavelength-neatlines(j)%wavelength).lt.diff) then
              diff = abs(userlines(i)%wavelength-neatlines(j)%wavelength)
              assign_1 = i
              assign_2 = j
            endif 
          enddo
          if (diff .le. rms) then
            print "(F8.2,A30,F8.2,1X,A10,A21)",userlines(assign_1)%wavelength+shift," =  ",neatlines(assign_2), "(ID = certain)"
            write (100,105) neatlines(assign_2)%wavelength, userlines(assign_1)%flux, userlines(assign_1)%uncertainty
            count=count+1
          elseif (diff .gt. rms .and. diff .le. 2*rms) then
            print "(F8.2,A30,F8.2,1X,A10,A21)",userlines(assign_1)%wavelength+shift," = ",neatlines(assign_2), "(ID = very likely)"
            write (100,105) neatlines(assign_2)%wavelength, userlines(assign_1)%flux, userlines(assign_1)%uncertainty
            count=count+1
          elseif (diff .gt. 2*rms .and. diff .le. 3*rms) then
            print "(F8.2,A30,F8.2,1X,A10,A21)",userlines(assign_1)%wavelength+shift," = ",neatlines(assign_2), "(ID = probable)"
            write (100,105) neatlines(assign_2)%wavelength, userlines(assign_1)%flux, userlines(assign_1)%uncertainty
            count=count+1
          else
            write (100,105) userlines(assign_1)%wavelength+shift, userlines(assign_1)%flux, userlines(assign_1)%uncertainty
            print "(F8.2,A30,F8.2,1X,A10,A21)",userlines(assign_1)%wavelength," unrecognised. nearest known: ",neatlines(assign_2)
          endif
        enddo

print *,count," lines identified; ",n_userlines-count," unidentified"
print *,"All lines written to file ",trim(filename),".proc"
print *,"Wavelengths of identified lines changed as necessary"

end program linefinder
