program linefinder

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

type(neat_line), dimension(:), allocatable :: neatlines
type(user_line), dimension(:), allocatable :: userlines
type(assigned_line), dimension(:), allocatable :: assignedlines
integer :: I, J, n_neatlines, n_userlines, listlength, Iint, Iread, IO, assign_1, assign_2, count
double precision :: temp_wave, temp_flux, temp_unc, diff
character*5 :: temp_ion1, temp_ion2
character*1024 :: filename

I = IARGC()
if (I .ne. 1) then
  print *,"No line list specified"
  stop
endif

CALL getarg(1,filename)

!algorithm:

!1. read NEAT line list in

        I = 1
        OPEN(100, file='manual/lines', iostat=IO, status='old')
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

!find the three closest matches in wavelength to given
!if diff = 0 then assign
!if diff < 0.1 then assign and flag
!if diff > 0.1 then don't assign

        count=0
        open (100, file=trim(filename)//".proc", status='replace', action='write')
        105 format (F8.2,1X,F10.5,1X,F10.5)

        do I=1,n_userlines
          diff = 1.e10
          do j=1,n_neatlines
            if (abs(userlines(i)%wavelength-neatlines(j)%wavelength).lt.diff) then
              diff = abs(userlines(i)%wavelength-neatlines(j)%wavelength)
              assign_1 = i
              assign_2 = j
            endif 
          enddo
          if (diff .le. 1.e-3) then
            print "(F8.2,A29,F8.2,1X,A10)",userlines(assign_1)%wavelength," identified as ",neatlines(assign_2)
            write (100,105) neatlines(assign_2)%wavelength, userlines(assign_1)%flux, userlines(assign_1)%uncertainty
            count=count+1
          elseif (diff .gt. 1.e-3 .and. diff .le. 0.02) then
            print "(F8.2,A29,F8.2,1X,A10)",userlines(assign_1)%wavelength," very likely ",neatlines(assign_2)
            write (100,105) neatlines(assign_2)%wavelength, userlines(assign_1)%flux, userlines(assign_1)%uncertainty
            count=count+1
          elseif (diff .gt. 0.02 .and. diff .le. 0.1) then
            print "(F8.2,A29,F8.2,1X,A10)",userlines(assign_1)%wavelength," could be ",neatlines(assign_2)
            write (100,105) neatlines(assign_2)%wavelength, userlines(assign_1)%flux, userlines(assign_1)%uncertainty
            count=count+1
          else
            write (100,105) userlines(assign_1)%wavelength, userlines(assign_1)%flux, userlines(assign_1)%uncertainty
            print "(F8.2,A29,F8.2,1X,A10)",userlines(assign_1)%wavelength," unidentified. nearest line: ",neatlines(assign_2)
          endif
        enddo

print *,count," lines identified; ",n_userlines-count," unidentified"
print *,"All lines written to file ",trim(filename),".proc"
print *,"Wavelengths of identified lines changed as necessary"

end program linefinder
