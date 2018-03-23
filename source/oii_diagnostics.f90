!oii_diagnostics, routines to interpolate data from Storey et al 2017 to estimate the temperature and density from OII recombination line ratios
!(C) Roger Wesson
module mod_oii_diagnostics

implicit none
private :: dp
integer, parameter :: dp = kind(1.d0)
real(kind=dp), dimension(9,8,4) :: oiidiagnostics

contains

subroutine read_oii_diagnosticdata

implicit none
integer, parameter :: dp = kind(1.d0)
integer :: i,j
character(len=7) :: x
character(len=5) :: in1, in2

!debugging
#ifdef CO
        print *,"subroutine: read_oii_diagnosticdata"
#endif

!read data files

open (999, file=trim(PREFIX)//"/share/neat/te_oii.dat")
do i=1,4
  read (999,*) x !4 lines of comments ignored
enddo
do i=1,9
  do j=0,8
    if (j.eq.0) then
      read (999,*) x
    else
      oiidiagnostics(i,j,1)=2.4+0.2*i !temperature
      oiidiagnostics(i,j,2)=1.6+0.4*j !density
      read (999,*) in1, in2
      read (in1,"(F5.3)") oiidiagnostics(i,j,3)
      read (in2,"(F5.3)") oiidiagnostics(i,j,4)
    endif
  enddo
enddo
close(999)

oiidiagnostics(:,8,2)=5.0 ! all lines separated by 0.4dex except top density, 0.6 dex higher

end subroutine read_oii_diagnosticdata

subroutine oii_rl_diagnostics(obs4649_4089, obs4649_4662,oii_te,oii_ne)

implicit none
real(kind=dp), dimension(4) :: a,b,c, x0, y0, d !equations of sides, locations of closest points on lines to value, distances to closest points
real(kind=dp) :: obs4649_4089, obs4649_4662, oii_te, oii_ne
integer :: i,j

!debugging
#ifdef CO
        print *,"subroutine: oii_rl_diagnostics"
#endif

oii_te=0.d0
oii_ne=0.d0
x0=0.d0
y0=0.d0

!find where our value is by first finding which box it's in
!get equation for line joining points on each side of the box
!get point on each of those lines closest to our value
!use this to determine if point is between two lines in both directions
!if it is, then we're in the box and the magnitude of the vectors gives us an interpolated value

do i=1,8
  do j=1,7

    !left hand side is (i,j,3),(i,j,4) to (i,j+1,3),(i,j+1,4)
    !right hand side is (i+1,j,3),(i+1,j,4) to (i+1,j+1,3),(i+1,j+1,4)
    !bottom is (i,j,3),(i,j,4) to (i+1,j,3),(i+1,j,4)
    !top is (i,j+1,3),(i,j+1,4) to (i+1,j+1,3),(i+1,j+1,4)

    !equation of left hand side, ax+by+c=0, a, b and c are
    !a=(i,j+1,4)-(i,j,4)/(i,j+1,3)-(i,j,3), b=-1, c=-a*(i,j,3)-b*(i,j,4)

    if ((oiidiagnostics(i,j+1,3)-oiidiagnostics(i,j,3)).gt.0) a(1)=(oiidiagnostics(i,j+1,4)-oiidiagnostics(i,j,4))/(oiidiagnostics(i,j+1,3)-oiidiagnostics(i,j,3))
    if ((oiidiagnostics(i+1,j+1,3)-oiidiagnostics(i+1,j,3)).gt.0) a(2)=(oiidiagnostics(i+1,j+1,4)-oiidiagnostics(i+1,j,4))/(oiidiagnostics(i+1,j+1,3)-oiidiagnostics(i+1,j,3))
    if ((oiidiagnostics(i+1,j+1,3)-oiidiagnostics(i,j,3)).gt.0) a(3)=(oiidiagnostics(i+1,j,4)-oiidiagnostics(i,j,4))/(oiidiagnostics(i+1,j+1,3)-oiidiagnostics(i,j,3))
    if ((oiidiagnostics(i+1,j+1,3)-oiidiagnostics(i,j+1,3)).gt.0) a(4)=(oiidiagnostics(i+1,j+1,4)-oiidiagnostics(i,j+1,4))/(oiidiagnostics(i+1,j+1,3)-oiidiagnostics(i,j+1,3))

    b=-1

    c(1)=-a(1)*oiidiagnostics(i,j,3)-b(1)*oiidiagnostics(i,j,4)
    c(2)=-a(2)*oiidiagnostics(i+1,j,3)-b(2)*oiidiagnostics(i+1,j,4)
    c(3)=-a(3)*oiidiagnostics(i,j,3)-b(3)*oiidiagnostics(i,j,4)
    c(4)=-a(4)*oiidiagnostics(i,j+1,3)-b(4)*oiidiagnostics(i,j+1,4)

    where(a**2+b**2.gt.0)
      x0=(b*(b*obs4649_4089-a*obs4649_4662)-(a*c))/(a**2+b**2)
      y0=(a*(-b*obs4649_4089+a*obs4649_4662)-(b*c))/(a**2+b**2)
    endwhere

    if (x0(1).eq.0.d0) then !condition to check for vertical
      x0(1)=oiidiagnostics(i,j,3)
      y0(1)=1.
    endif
    if (x0(2).eq.0.d0) then
      x0(2)=oiidiagnostics(i+1,j,3)
      y0(2)=1.
    endif

    d=((obs4649_4089-x0)**2+(obs4649_4662-y0)**2)**0.5

    if (obs4649_4089 - x0(1) .ge. 0 .and. obs4649_4089 - x0(2) .le. 0 .and. obs4649_4662 - y0(3) .ge. 0 .and. obs4649_4662 - y0(4) .le. 0) then
      oii_te=10**(oiidiagnostics(i,j,1)+(oiidiagnostics(i+1,j,1)-oiidiagnostics(i,j,1))*(d(1)/(d(1)+d(2))))
      oii_ne=10**(oiidiagnostics(i,j,2)+(oiidiagnostics(i,j+1,2)-oiidiagnostics(i,j,2))*(d(3)/(d(3)+d(4))))
    endif

  enddo
enddo

end subroutine oii_rl_diagnostics
end module mod_oii_diagnostics
