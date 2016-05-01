module mod_weights
use mod_abundtypes
use mod_abundIO
contains

subroutine setweights(weights,linelist,ILs)
!will read configuration file in containing all the weights. currently just sets them.

        implicit none
        integer, parameter :: dp=kind(1.d0)
        type(weightingarray) :: weights
        type(line),dimension(:) :: linelist
        type(cel),dimension(82) :: ILs
        real(kind=dp), dimension(82) :: celweights
        integer :: i

! define weighting array. to be replaced with reading it in from a configuration file
!low ionisation densities
        weights%oiiDens = 1.d0
        weights%siiDens = 1.d0
!low ionisation temperatures
        weights%oiiTemp = 1.d0
        weights%siiTemp = 1.d0
        weights%niiTemp = 5.d0
        weights%ciTemp = 1.d0
        weights%oiTemp = 1.d0
!medium ionisation densities
        weights%cliiiDens = 1.d0
        weights%ciiiDens = 1.d0
        weights%arivDens = 1.d0
        weights%oiiiIRDens = 0.d0
        weights%ariiiIRDens = 0.d0
        weights%siiiIRDens = 0.d0
        weights%neiiiIRDens = 0.d0
!medium ionisation temperatures
        weights%oiiiTemp = 4.d0
        weights%siiiTemp = 1.d0
        weights%ariiiTemp = 2.d0
        weights%neiiiTemp = 2.d0
        weights%neiiiIRTemp = 0.d0
        weights%oiiiIRTemp = 0.d0
        weights%oiiiUVTemp = 0.d0
!high ionisation densities
        weights%neivDens = 1.d0
!high ionisation temperatures
        weights%arvTemp = 1.d0
        weights%nevTemp = 1.d0
!balmer weights for extinction calculation
        weights%ha = -1.d0
        weights%hg = -1.d0
        weights%hd = -1.d0
!helium abundance weights
        weights%he4471 = 1.d0
        weights%he5876 = 3.d0
        weights%he6678 = 1.d0
!CEL abundance weights
!this horrifically impenetrable and undebuggable array will be replaced with sensible and intelligible file reading
        celweights = (/ 0.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, 0.d0, -1.d0, -1.d0, 0.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, 0.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, 0.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, 0.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, 0.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, 0.d0, -1.d0, -1.d0, -1.d0 /)

!for each ion, get_ion(name,ILs)%location tells us if it was detected
!if it was, then set its weight. if it wasn't, do nothing.

        do i=1,size(ILs)
          if (ILs(i)%location .gt. 0) then
            if (celweights(i) .eq. -1.d0) then !weight by observed intensity
              linelist(ILs(i)%location)%weight = linelist(ILs(i)%location)%intensity
            else !weight by specified weighting
              linelist(ILs(i)%location)%weight = celweights(i)
            endif
          endif
        enddo

end subroutine setweights

end module
