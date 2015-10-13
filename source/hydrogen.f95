module mod_hydrogen
use mod_abundtypes
contains

subroutine balmer_densities(H_BS,medtemp,density)

implicit none
integer, parameter :: dp = kind(1.d0)

type(line), dimension(:) :: H_BS
real(kind=dp), dimension(10:25) :: ratios,densities
real(kind=dp), dimension(:,:,:,:), allocatable :: emissivities
real(kind=dp), dimension(:,:,:), allocatable :: searcharray
integer :: ntemps, ndens, nlevs
real(kind=dp) :: density, interp_t, weight
real(kind=dp), intent(in) :: medtemp
character(len=1) :: junk
character(len=20), dimension(8) :: invar
integer :: i,j,k,l !counters
real(kind=dp) :: r1, r2
integer :: H ! temporary, which line to look at
real(kind=dp), dimension(:), allocatable :: temperatures

density=0.d0
weight=0.d0

!nb H_BS indexing is such that entry 1 contains H3, etc
!Thus H10 is in entry 8
!allocate arrays to store values, high order lines (n>10) only
ratios = H_BS(8:23)%int_dered/H_BS(2)%int_dered
densities=0.d0

!first read in the data

open (unit=357, file="Atomic-data/RHi.dat")
read (357,*) junk !first line is a comment
read (357,"(I3,I3)") ntemps, ndens !second line has number of temperatures and densities
nlevs=25 ! maximum number of levels shown in intrat data file

!allocate the emissivities array, dimensions are temperature, density, level 1, level 2

allocate(emissivities(ntemps, ndens, nlevs, nlevs))
allocate(searcharray(ndens,nlevs,nlevs))
allocate(temperatures(ntemps))

emissivities(:,:,:,:)=0.d0

!now, loop through the temperatures and densities and read in the emissivities

do i=1,ntemps
  do j=1,ndens
    read (357,*) invar(1:6) !line at top of each block has density, charge, temperature, case, maximum level calculated and maximum level displayed
    if (j.eq.1) then
      read(invar(3),"(E9.2)") temperatures(i) ! get the temperatures for later use in interpolating
    endif
    do k=nlevs,1,-1
      do l=1,k-1
        read (357,"(E10.3)", advance="no") emissivities(i,j,k,l)
      enddo
    enddo
  enddo
enddo
close (357)

do i=1,ntemps
  if (medtemp .lt. temperatures(i)) then
    exit
  endif
end do

!interpolate linearly between temperatures to get line ratios with density at specified temperature
!if temperature is outside data limits, it is just set to the limit
!lower limit is 500K, upper is 30000K

if (medtemp .lt. temperatures(1) .or. medtemp .gt. temperatures(size(temperatures))) then
  searcharray(:,:,:)=emissivities(i,:,:,:) 
else
  interp_t=(medtemp-temperatures(i-1))/(temperatures(i)-temperatures(i-1))
  searcharray(:,:,:)=emissivities(i,:,:,:)+interp_t*(emissivities(i+1,:,:,:)-emissivities(i,:,:,:))
endif

!now go through H_BS from 10 to 25, calculate the density for each line ratio, store it
!intrat data only goes up to n=25

do H=10,25
  if (H_BS(H-2)%int_dered .gt. 0.d0) then !line is present, calculate a density
    if (ratios(H) .lt. (searcharray(1,H,2)/searcharray(1,4,2))) then
      densities(H) = 1.
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
    density=density+H_BS(i-2)%int_dered*densities(i)
    weight=weight+H_BS(i-2)%int_dered
  endif
enddo

density=10**(density/weight)

!check if lower limit, Storey and Hummer go down to 100cm-3

if (density .lt. 100.) then
  density = 100.d0
endif

!deallocate arrays

deallocate(emissivities)
deallocate(searcharray)
deallocate(temperatures)

end subroutine balmer_densities

end module mod_hydrogen
