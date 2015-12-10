module mod_hydrogen
use mod_abundtypes

implicit none
private :: dp
integer, parameter :: dp = kind(1.d0)
real(kind=dp), dimension(:,:,:,:), allocatable :: hidata
real(kind=dp), dimension(:), allocatable :: temperatures
integer :: ntemps, ndens, nlevs

contains

subroutine read_hydrogen
!get Case B emissivities from file, values are from Storey and Hummer 1995.
implicit none
character(len=1) :: junk
character(len=20), dimension(8) :: invar
integer :: i,j,k,l !counters

open (unit=357, file="Atomic-data/RHi.dat")
read (357,*) junk !first line is a comment
read (357,"(I3,I3)") ntemps, ndens !second line has number of temperatures and densities
nlevs=25 ! maximum number of levels shown in intrat data file

!allocate the emissivities array, dimensions are temperature, density, level 1, level 2

allocate(hidata(ntemps, ndens, nlevs, nlevs))
allocate(temperatures(ntemps))

hidata(:,:,:,:)=0.d0

!now, loop through the temperatures and densities and read in the emissivities

do i=1,ntemps
  do j=1,ndens
    read (357,*) invar(1:6) !line at top of each block has density, charge, temperature, case, maximum level calculated and maximum level displayed
    if (j.eq.1) then
      read(invar(3),"(E9.2)") temperatures(i) ! get the temperatures for later use in interpolating
    endif
    do k=nlevs,1,-1
      do l=1,k-1
        read (357,"(E10.3)", advance="no") hidata(i,j,k,l)
      enddo
    enddo
  enddo
enddo
close (357)

end subroutine read_hydrogen

subroutine balmer_densities(H_Balmer,medtemp,density)

implicit none
integer, parameter :: dp = kind(1.d0)

type(line), dimension(:) :: H_Balmer
real(kind=dp), dimension(10:25) :: ratios,densities
real(kind=dp), dimension(:,:,:), allocatable :: searcharray
real(kind=dp) :: density, interp_t, weight
real(kind=dp), intent(in) :: medtemp
integer :: i,j !counters
real(kind=dp) :: r1, r2
integer :: H ! temporary, which line to look at

density=0.d0
weight=0.d0

!nb H_Balmer indexing is such that entry 1 contains H3, etc
!Thus H10 is in entry 8
!allocate arrays to store values, high order lines (n>10) only
ratios = H_Balmer(8:23)%int_dered/H_Balmer(2)%int_dered
densities=0.d0

!allocate the search array.  emissivities array has dimensions of temperature, density, level 1, level 2
!search array just has dimensions of density, level 1, level 2

allocate(searcharray(size(hidata,2),size(hidata,3),size(hidata,4)))

!now find the temperature

do i=1,ntemps
  if (medtemp .lt. temperatures(i)) then
    exit
  endif
end do

!interpolate linearly between temperatures to get line ratios with density at specified temperature
!if temperature is outside data limits, it is just set to the limit
!lower limit is 500K, upper is 30000K

if (medtemp .lt. temperatures(1) .or. medtemp .gt. temperatures(size(temperatures))) then
  searcharray(:,:,:)=hidata(i,:,:,:)
else
  interp_t=(medtemp-temperatures(i-1))/(temperatures(i)-temperatures(i-1))
  searcharray(:,:,:)=hidata(i-1,:,:,:)+interp_t*(hidata(i,:,:,:)-hidata(i,:,:,:))
endif

!now go through H_Balmer from 10 to 25, calculate the density for each line ratio, store it
!intrat data only goes up to n=25

do H=10,25
  if (H_Balmer(H-2)%int_dered .gt. 0.d0) then !line is present, calculate a density
    if (ratios(H) .lt. (searcharray(1,H,2)/searcharray(1,4,2))) then
      densities(H) = 2.
    elseif (ratios(H) .gt. (searcharray(ndens,H,2)/searcharray(ndens,4,2))) then
      densities(H) = 8.
    else !search for the value
      do j=1,ndens-1
        r1=searcharray(j,H,2)/searcharray(j,4,2)
        r2=searcharray(j+1,H,2)/searcharray(j+1,4,2)
        if (ratios(H) .gt. r1 .and. ratios(H) .lt. r2) then
          densities(H)=j+1+(ratios(H)-r1)/(r2-r1)
        endif
      enddo
    endif
  endif
end do

!now we have an array with all the densities implied by the ratios.  Derive a flux weighted value

do i=10,25
  if (densities(i).gt.0.d0) then
    density=density+H_Balmer(i-2)%int_dered*densities(i)
    weight=weight+H_Balmer(i-2)%int_dered
  endif
enddo

if (weight .gt. 0.d0) then
  density=10**(density/weight)
else
  density=0.d0
endif

!check if lower limit, Storey and Hummer go down to 100cm-3

if (density .lt. 100. .and. density .ne. 0.d0) then
  density = 100.d0
endif

!deallocate search array

deallocate(searcharray)

end subroutine balmer_densities

subroutine paschen_densities(H_Paschen,medtemp,density)
! computed using ratios of Paschen lines to Hbeta
! will change to P9, checking if it is observed
implicit none
integer, parameter :: dp = kind(1.d0)

type(line), dimension(:) :: H_Paschen
real(kind=dp), dimension(10:25) :: ratios,densities
real(kind=dp), dimension(:,:,:), allocatable :: searcharray
real(kind=dp) :: density, interp_t, weight
real(kind=dp), intent(in) :: medtemp
integer :: i,j !counters
real(kind=dp) :: r1, r2
integer :: H ! temporary, which line to look at

density=0.d0
weight=0.d0

!nb H_Paschen indexing is such that entry 1 contains H4, etc
!Thus H10 is in entry 7
!allocate arrays to store values, high order lines (n>10) only
ratios = H_Paschen(7:22)%int_dered/100.d0
densities=0.d0

!allocate the search array.  emissivities array has dimensions of temperature, density, level 1, level 2
!search array just has dimensions of density, level 1, level 2

allocate(searcharray(size(hidata,2),size(hidata,3),size(hidata,4)))

!now find the temperature

do i=1,ntemps
  if (medtemp .lt. temperatures(i)) then
    exit
  endif
end do

!interpolate linearly between temperatures to get line ratios with density at specified temperature
!if temperature is outside data limits, it is just set to the limit
!lower limit is 500K, upper is 30000K

if (medtemp .lt. temperatures(1) .or. medtemp .gt. temperatures(size(temperatures))) then
  searcharray(:,:,:)=hidata(i,:,:,:)
else
  interp_t=(medtemp-temperatures(i-1))/(temperatures(i)-temperatures(i-1))
  searcharray(:,:,:)=hidata(i-1,:,:,:)+interp_t*(hidata(i,:,:,:)-hidata(i,:,:,:))
endif

!now go through H_Paschen from 10 to 25, calculate the density for each line ratio, store it
!intrat data only goes up to n=25

do H=10,25
  if (H_Paschen(H-3)%int_dered .gt. 0.d0) then !line is present, calculate a density
    if (ratios(H) .lt. (searcharray(1,H,3)/searcharray(1,4,2))) then
      densities(H) = 2.
    elseif (ratios(H) .gt. (searcharray(ndens,H,3)/searcharray(ndens,4,2))) then
      densities(H) = 8.
    else !search for the value
      do j=1,ndens-1
        r1=searcharray(j,H,3)/searcharray(j,4,2)
        r2=searcharray(j+1,H,3)/searcharray(j+1,4,2)
        if (ratios(H) .gt. r1 .and. ratios(H) .lt. r2) then
          densities(H)=j+1+(ratios(H)-r1)/(r2-r1)
        endif
      enddo
    endif
  endif
end do

!now we have an array with all the densities implied by the ratios.  Derive a flux weighted value

do i=10,25
  if (densities(i).gt.0.d0) then
    density=density+H_Paschen(i-3)%int_dered*densities(i)
    weight=weight+H_Paschen(i-3)%int_dered
  endif
enddo

if (weight .gt. 0.d0) then
  density=10**(density/weight)
else
  density=0.d0
endif

!check if lower limit, Storey and Hummer go down to 100cm-3

if (density .lt. 100. .and. density .ne. 0.d0) then
  density = 100.d0
endif

!deallocate search array

deallocate(searcharray)

end subroutine paschen_densities

end module mod_hydrogen
