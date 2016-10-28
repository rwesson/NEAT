module mod_helium
use mod_abundtypes

      implicit none
      private :: dp
      integer, parameter :: dp = kind(1.d0)
      real(kind=dp), dimension(:,:,:,:), allocatable :: heiidata
      real(kind=dp), dimension(:), allocatable :: temperatures
      integer :: ntemps, ndens, nlevs

      contains

      subroutine get_heii_abund(te, ne, IHeII4686, heiiabund)

      implicit none
      integer, parameter :: dp = kind(1.d0)

      real(kind=dp) :: A4686,TE,NE,IHeII4686,heiiabund

!debugging
#ifdef CO
        print *,"subroutine: get_heii_abund"
#endif

      A4686=10.**(GAMM4861(TE,NE)-GAMM4686(TE,NE))
      heiiabund=IHeII4686/100.*A4686

end subroutine get_heii_abund

subroutine get_hei_porter(linelist,Te, ne, he_lines, heidata, Heiabund, weight4471, weight5876, weight6678)

      use mod_abundtypes

      IMPLICIT NONE
      real(kind=dp) :: te, ne
      real(kind=dp) :: AB4471,AB5876,AB6678,heiabund
      real(kind=dp) :: weight4471, weight5876, weight6678 ! weights for abundance determination

      type(line), dimension(:) :: linelist
      integer, dimension(44) :: he_lines
      real(kind=dp), dimension(21,14,44), intent(in) :: heidata
      real(kind=dp), dimension(44) :: emissivities
      integer :: i

!debugging
#ifdef CO
        print *,"subroutine: get_hei_porter. Te=",Te,", ne=",ne,",weights=",weight4471, weight5876, weight6678
#endif

!      data is for the following lines, in this order: 2945.10,3188.74,3613.64,3888.65,3964.73,4026.21,4120.82,4387.93,4437.55,4471.50,4713.17,4921.93,5015.68,5047.74,5875.66,6678.16,7065.25,7281.35,9463.58,10830.25,11013.07,11969.06,12527.49,12755.69,12784.92,12790.50,12845.98,12968.43,12984.88,13411.69,15083.65,17002.40,18555.57,18685.33,18697.21,19089.36,19543.19,20424.97,20581.28,20601.76,21120.12,21132.03,21607.80,21617.01 /)

      do i = 1,44
        call get_emissivity_porter(te,ne, i, emissivities(i), heidata)
        if (He_lines(i).gt.0) then
          linelist(He_lines(i))%abundance = linelist(He_lines(i))%int_dered/100. * 10.**(GAMM4861(TE,NE)-emissivities(i))
        endif
      enddo

      if (He_lines(10).gt.0) then
        AB4471 = linelist(He_lines(10))%abundance
      else
        AB4471 = 0.d0
        weight4471 = 0.d0
      endif

      if (He_lines(15).gt.0) then
        AB5876 = linelist(He_lines(15))%abundance
      else
        AB5876 = 0.d0
        weight5876 = 0.d0
      endif

      if (He_lines(16).gt.0) then
        AB6678 = linelist(He_lines(16))%abundance
      else
        AB6678 = 0.d0
        weight6678 = 0.d0
      endif

      if (weight4471+weight5876+weight6678.gt.0.0) then
         heiabund = ((weight4471*AB4471) + (weight5876*AB5876) + (weight6678*AB6678)) / (weight4471 + weight5876 + weight6678)
      else
         heiabund = 0.D0
      endif

end subroutine get_hei_porter

subroutine get_emissivity_porter(te, ne, line, emissivity, heidata)

        implicit none
        integer, parameter :: dp = kind(1.d0)

        real(kind=dp), intent(in) :: te, ne
        real(kind=dp) :: testart, nestart, logne, te_copied
        real(kind=dp) :: interp_factor_te, interp_factor_ne, interp_t1, interp_t2, emissivity
        real(kind=dp), dimension(21,14,44), intent(in) :: heidata
        integer :: i,j,line

!debugging
#ifdef CO
        print "(A,F8.0,A,F8.0)"," subroutine: get_emissivity_porter. te=",te,", ne=",ne
#endif

        ! local variables. ne needs to be logarithmic, input is linear

        te_copied=te
        logne=log10(ne)

        !do bilinear interpolation of the log values
        !find starting values

        testart=25000.
        nestart=15.

        !find temperature box
        i=floor((te-4000)/1000.)
        testart=1000.*floor(te/1000.)

        !find density box
        j=floor(logne)
        nestart=floor(logne)

        !now we have the array positions in i and j and the values in nestart and testart
        !check if we are within the Porter limits of 5000<Te<25000, 1<log(ne)<14

        if (testart .le. 5000.) then !set variables to return emissivity for 5000K. todo: print warning
          testart=5000.
          te_copied=5000.
          i=1
        elseif (testart .ge. 25000.) then !set variables to return emissivitity for 25000K. todo: print warning
          testart=24000.
          te_copied=25000.
          i=20
        endif
        if (nestart .le. 1) then !return emissivity for log(ne)=1
          nestart=1.
          logne=1.
          j=1
        elseif (nestart .ge. 14) then !return emissivity for log(ne)=14
          nestart=13
          logne=14
          j=13
        endif

        !now interpolate first in T and then in ne

        interp_factor_te = (te_copied-testart)/1000.
        interp_factor_ne = (logne-nestart)

        !get two interpolated values for each ne, then interpolate between those

        interp_t1 = heidata(i,j,line)+interp_factor_te*(heidata(i+1,j,line)-heidata(i,j,line))
        interp_t2 = heidata(i,j+1,line)+interp_factor_te*(heidata(i+1,j+1,line)-heidata(i,j+1,line))

        emissivity = interp_t1+(interp_factor_ne*(interp_t2-interp_t1))

end subroutine get_emissivity_porter

subroutine get_hei_smits_new(linelist,Te, ne, he_lines, heidata, Heiabund, weight4471, weight5876, weight6678)

      use mod_abundtypes

      IMPLICIT NONE
      real(kind=dp) :: te, ne, tereduced
      real(kind=dp) :: AB4471,AB5876,AB6678,heiabund
      real(kind=dp) :: c4471, c5876, c6678, d ! corrections for collisional excitation
      real(kind=dp) :: weight4471, weight5876, weight6678 ! weights for abundance determination

      type(line), dimension(:) :: linelist
      integer, dimension(44), intent(in) :: he_lines
      real(kind=dp), dimension(3,6,44), intent(in) :: heidata
      real(kind=dp), dimension(44,3) :: emissivities
      real(kind=dp) :: interpolatedemissivity
      integer :: i,j

!debugging
#ifdef CO
        print *,"subroutine: get_hei_smits_new. Te=",Te,", ne=",ne,",weights=",weight4471, weight5876, weight6678
#endif

!      data is for the following lines, in this order: 2945.10,3188.74,3613.64,3888.65,3964.73,4026.21,4120.82,4387.93,4437.55,4471.50,4713.17,4921.93,5015.68,5047.74,5875.66,6678.16,7065.25,7281.35,9463.58,10830.25,11013.07,11969.06,12527.49,12755.69,12784.92,12790.50,12845.98,12968.43,12984.88,13411.69,15083.65,17002.40,18555.57,18685.33,18697.21,19089.36,19543.19,20424.97,20581.28,20601.76,21120.12,21132.03,21607.80,21617.01 /)

      do i = 1,44
        do j = 1,3
          emissivities(i,j)=heidata(j,1,i) + heidata(j,2,i)*log10(te) + heidata(j,3,i)*log10(te)**2 + heidata(j,4,i)*log10(te)**3 + heidata(j,5,i)*log10(te)**4 + heidata(j,6,i)
        enddo

        if (log10(ne).lt.2) then
          interpolatedemissivity=emissivities(i,1)
        elseif (log10(ne).ge.2..and.log10(ne).lt.4) then
          interpolatedemissivity=emissivities(i,1) + (emissivities(i,2) - emissivities(i,1)) * (log10(ne) - 2)
        elseif (log10(ne).ge.4..and.log10(ne).lt.6) then
          interpolatedemissivity=emissivities(i,2) + (emissivities(i,3) - emissivities(i,2)) * (log10(ne) - 4)
        else
          interpolatedemissivity=emissivities(i,3)
        endif

        if (interpolatedemissivity.ne. 0.D0 .and. He_lines(i) .gt. 0) then
          linelist(He_lines(i))%abundance = linelist(He_lines(i))%int_dered/100. * 10.**(GAMM4861(TE,NE)-interpolatedemissivity)
        else
          linelist(He_lines(i))%abundance = 0.D0
        endif

      enddo

      tereduced = te/10000.

!correct 4471, 5876 and 6678 for collisional contributions

      D=1.+3130.*tereduced**(-0.50)/ne
      C4471=(  6.95*tereduced**(0.15)*exp(-4.545/tereduced)  + 0.22*tereduced**(-0.55)*exp(-4.884/tereduced)  )/D
      C5876=(  6.78*tereduced**(0.07)*exp(-3.776/tereduced)  + 1.67*tereduced**(-0.15)*exp(-4.545/tereduced) + 0.60*tereduced**(-0.34)*exp(-4.901/tereduced)  )/D
      C6678=(  3.15*tereduced**(-0.54)*exp(-3.776/tereduced) + 0.51*tereduced**(-0.51)*exp(-4.545/tereduced) + 0.20*tereduced**(-0.66)*exp(-4.901/tereduced)  )/D

      linelist(He_lines(10))%abundance=linelist(He_lines(10))%abundance/(1.+C4471)
      linelist(He_lines(15))%abundance=linelist(He_lines(15))%abundance/(1.+C5876)
      linelist(He_lines(16))%abundance=linelist(He_lines(16))%abundance/(1.+C6678)

      AB4471 = linelist(He_lines(10))%abundance
      AB5876 = linelist(He_lines(15))%abundance
      AB6678 = linelist(He_lines(16))%abundance

           if (AB4471 .eq. 0) weight4471 = 0.d0
           if (AB5876 .eq. 0) weight5876 = 0.d0
           if (AB6678 .eq. 0) weight6678 = 0.d0

           if (weight4471+weight5876+weight6678.gt.0.0) then
               heiabund = ((weight4471*AB4471) + (weight5876*AB5876) + (weight6678*AB6678)) / (weight4471 + weight5876 + weight6678)
           else
               heiabund = 0.D0
           endif

end subroutine get_hei_smits_new

real(kind=dp) function GAMM4861(TE,NE)
!     This function determines the value of Log10 (gamm(H Beta))
!     = Log10( 4*Pai*j(HBeta)/NpNe) at temperature Te and density Ne
!     Storey P. J., Hummer D. G., 1995, MNRAS, 272, 41

      implicit none
      integer, parameter :: dp = kind(1.d0)

      real(kind=dp) TE,NE,AE2,AE3,AE4,AE5,AE6,AE7,AE8,AEFF,HCLL
      real(kind=dp) LNE, LTE

!debugging
#ifdef CO
        !print *,"function: GAMM4861"
#endif

      HCLL=-11.38871 ! = Log10 ( h * c / lambda(H-beta) ) - [ cgs units
      LNE = log10(NE)
      LTE = log10(TE)

      AE2 =(-9.06524E+00)+                                              &
     &     (-2.69954E+00)* LTE +                                        &
     &       8.80123E-01 * LTE ** 2 +                                   &
     &     (-1.57946E-01)* LTE ** 3 +                                   &
     &       9.25920E-03 * LTE ** 4

      AE3 =(-8.13757E+00)+                                              &
     &     (-3.57392E+00)* LTE +                                        &
     &       1.19331E+00 * LTE ** 2 +                                   &
     &     (-2.08362E-01)* LTE ** 3 +                                   &
     &       1.23303E-02 * LTE ** 4

      AE4 =(-6.87230E+00)+                                              &
     &     (-4.72312E+00)* LTE +                                        &
     &       1.58890E+00 * LTE ** 2 +                                   &
     &     (-2.69447E-01)* LTE ** 3 +                                   &
     &       1.58955E-02 * LTE ** 4

      AE5 =(-5.15059E+00)+                                              &
     &     (-6.24549E+00)* LTE +                                        &
     &       2.09801E+00 * LTE ** 2 +                                   &
     &     (-3.45649E-01)* LTE ** 3 +                                   &
     &       2.01962E-02 * LTE ** 4

      AE6 =(-2.35923E+00)+                                              &
     &     (-8.75565E+00)* LTE +                                        &
     &       2.95600E+00 * LTE ** 2 +                                   &
     &     (-4.77584E-01)* LTE ** 3 +                                   &
     &       2.78852E-02 * LTE ** 4

      AE7 =  1.55373E+00 +                                              &
     &     (-1.21894E+01)* LTE +                                        &
     &       4.10096E+00 * LTE ** 2 +                                   &
     &     (-6.49318E-01)* LTE ** 3 +                                   &
     &       3.76487E-02 * LTE ** 4

      AE8 =  6.59883E+00 +                                              &
     &     (-1.64030E+01)* LTE +                                        &
     &       5.43844E+00 * LTE ** 2 +                                   &
     &     (-8.40253E-01)* LTE ** 3 +                                   &
     &       4.79786E-02 * LTE ** 4


      if (LNE.lt.2) then
        AEFF=AE2
      elseif (LNE.ge.2..AND.LNE.lt.3) then
        AEFF=AE2 + (AE3 - AE2) * (LNE - 2)
      elseif (LNE.ge.3..AND.LNE.lt.4) then
        AEFF=AE3 + (AE4 - AE3) * (LNE - 3)
      elseif (LNE.ge.4..AND.LNE.lt.5) then
        AEFF=AE4 + (AE5 - AE4) * (LNE - 4)
      elseif (LNE.ge.5..AND.LNE.lt.6) then
        AEFF=AE5 + (AE6 - AE5) * (LNE - 5)
      elseif (LNE.ge.6..AND.LNE.lt.7) then
        AEFF=AE6 + (AE7 - AE6) * (LNE - 6)
      elseif (LNE.ge.7..AND.LNE.lt.8) then
        AEFF=AE7 + (AE8 - AE7) * (LNE - 7)
      else
        AEFF=AE8
      endif

      GAMM4861=AEFF + HCLL

end function gamm4861

real(kind=dp) function GAMM6563(TE,NE)
!     This function determines the value of Log10 (gamm(H Alpha))
!     = Log10( 4*Pai*j(HBeta)/NpNe) at temperature Te and density Ne
!     Storey P. J., Hummer D. G., 1995, MNRAS, 272, 41

      implicit none
      integer, parameter :: dp = kind(1.d0)

      real(kind=dp) TE,NE,AE2,AE3,AE4,AE5,AE6,AE7,AE8,AEFF,HCLL
      real(kind=dp) LNE, LTE

!debugging
#ifdef CO
        !print *,"function: GAMM6563"
#endif

      HCLL=-11.51901 ! = Log10 ( h * c / lambda(H-alpha) ) - [ cgs units
      LNE = log10(NE)
      LTE = log10(TE)

      AE2 =(-9.78342E+00)+                                              &
     &     (-9.79048E-01)* LTE +                                        &
     &     ( 1.28883E-01)* LTE ** 2 +                                   &
     &     (-2.68907E-02)* LTE ** 3 +                                   &
     &     ( 1.15588E-03)* LTE ** 4

      AE3 =(-8.96834E+00)+                                              &
     &     (-1.81758E+00)* LTE +                                        &
     &     ( 4.52249E-01)* LTE ** 2 +                                   &
     &     (-8.22940E-02)* LTE ** 3 +                                   &
     &     ( 4.71377E-03)* LTE ** 4

      AE4 =(-7.54040E+00)+                                              &
     &     (-3.26881E+00)* LTE +                                        &
     &     ( 1.00538E+00)* LTE ** 2 +                                   &
     &     (-1.75962E-01)* LTE ** 3 +                                   &
     &     ( 1.06554E-02)* LTE ** 4

      AE5 =(-5.54310E+00)+                                              &
     &     (-5.24291E+00)* LTE +                                        &
     &     ( 1.73956E+00)* LTE ** 2 +                                   &
     &     (-2.97779E-01)* LTE ** 3 +                                   &
     &     ( 1.82635E-02)* LTE ** 4

      AE6 =(-1.66024E+00)+                                              &
     &     (-9.07781E+00)* LTE +                                        &
     &     ( 3.16498E+00)* LTE ** 2 +                                   &
     &     (-5.33982E-01)* LTE ** 3 +                                   &
     &     ( 3.29781E-02)* LTE ** 4

      AE7 =( 2.44816E+00)+                                              &
     &     (-1.27787E+01)* LTE +                                        &
     &     ( 4.41371E+00)* LTE ** 2 +                                   &
     &     (-7.20967E-01)* LTE ** 3 +                                   &
     &     ( 4.34582E-02)* LTE ** 4

      AE8 =( 8.09651E+00)+                                              &
     &     (-1.76766E+01)* LTE +                                        &
     &     ( 6.00896E+00)* LTE ** 2 +                                   &
     &     (-9.52069E-01)* LTE ** 3 +                                   &
     &     ( 5.60108E-02)* LTE ** 4


      if (LNE.lt.2) then
        AEFF=AE2
      elseif (LNE.ge.2..AND.LNE.lt.3) then
        AEFF=AE2 + (AE3 - AE2) * (LNE - 2)
      elseif (LNE.ge.3..AND.LNE.lt.4) then
        AEFF=AE3 + (AE4 - AE3) * (LNE - 3)
      elseif (LNE.ge.4..AND.LNE.lt.5) then
        AEFF=AE4 + (AE5 - AE4) * (LNE - 4)
      elseif (LNE.ge.5..AND.LNE.lt.6) then
        AEFF=AE5 + (AE6 - AE5) * (LNE - 5)
      elseif (LNE.ge.6..AND.LNE.lt.7) then
        AEFF=AE6 + (AE7 - AE6) * (LNE - 6)
      elseif (LNE.ge.7..AND.LNE.lt.8) then
        AEFF=AE7 + (AE8 - AE7) * (LNE - 7)
      else
        AEFF=AE8
      endif

      GAMM6563=AEFF + HCLL

end function gamm6563

real(kind=dp) function GAMM4686(TE,NE)
!     This function determines the value of Log10 (gamm(HeII4686))
!     = Log10( 4*Pai*j(HeII 4686)/N(He++)Ne) at temperature Te and densi
!     Storey P. J., Hummer D. G., 1995, MNRAS, 272, 41

      implicit none
      integer, parameter :: dp = kind(1.d0)

      real(kind=dp) TE,NE,AE2,AE3,AE4,AE5,AE6,AE7,AE8,AEFF,HCLL
      real(kind=dp) LNE, LTE

!debugging
#ifdef CO
        !print *,"function: GAMM4686"
#endif

      HCLL=-11.37272 ! = Log10 ( h * c / lambda(4686) ) - [ cgs units ]
      LNE = log10(NE)
      LTE = log10(TE)

      AE2 =(-8.14491E+00)+                                              &
     &     (-2.10904E+00)* LTE +                                        &
     &       6.31303E-01 * LTE ** 2 +                                   &
     &     (-1.25035E-01)* LTE ** 3 +                                   &
     &       7.95257E-03 * LTE ** 4

      AE3 =(-7.41098E+00)+                                              &
     &     (-2.82798E+00)* LTE +                                        &
     &       8.90834E-01 * LTE ** 2 +                                   &
     &     (-1.66196E-01)* LTE ** 3 +                                   &
     &       1.03803E-02 * LTE ** 4

      AE4 =(-6.20402E+00)+                                              &
     &     (-3.99470E+00)* LTE +                                        &
     &       1.30756E+00 * LTE ** 2 +                                   &
     &     (-2.31744E-01)* LTE ** 3 +                                   &
     &       1.42229E-02 * LTE ** 4

      AE5 =(-4.11426E+00)+                                              &
     &     (-5.97955E+00)* LTE +                                        &
     &       2.00560E+00 * LTE ** 2 +                                   &
     &     (-3.39997E-01)* LTE ** 3 +                                   &
     &       2.04854E-02 * LTE ** 4

      AE6 =(-7.88408E-01)+                                              &
     &     (-9.06439E+00)* LTE +                                        &
     &       3.06790E+00 * LTE ** 2 +                                   &
     &     (-5.01667E-01)* LTE ** 3 +                                   &
     &       2.96826E-02 * LTE ** 4

      AE7 =  4.22286E+00 +                                              &
     &     (-1.35870E+01)* LTE +                                        &
     &       4.58664E+00 * LTE ** 2 +                                   &
     &     (-7.27438E-01)* LTE ** 3 +                                   &
     &       4.22467E-02 * LTE ** 4

      AE8 =  1.10404E+01 +                                              &
     &     (-1.94998E+01)* LTE +                                        &
     &       6.49776E+00 * LTE ** 2 +                                   &
     &     (-1.00125E+00)* LTE ** 3 +                                   &
     &       5.69519E-02 * LTE ** 4


      if (LNE.lt.2) then
        AEFF=AE2
      elseif (LNE.ge.2..AND.LNE.lt.3) then
        AEFF=AE2 + (AE3 - AE2) * (LNE - 2)
      elseif (LNE.ge.3..AND.LNE.lt.4) then
        AEFF=AE3 + (AE4 - AE3) * (LNE - 3)
      elseif (LNE.ge.4..AND.LNE.lt.5) then
        AEFF=AE4 + (AE5 - AE4) * (LNE - 4)
      elseif (LNE.ge.5..AND.LNE.lt.6) then
        AEFF=AE5 + (AE6 - AE5) * (LNE - 5)
      elseif (LNE.ge.6..AND.LNE.lt.7) then
        AEFF=AE6 + (AE7 - AE6) * (LNE - 6)
      elseif (LNE.ge.7..AND.LNE.lt.8) then
        AEFF=AE7 + (AE8 - AE7) * (LNE - 7)
      else
        AEFF=AE8
      endif

      GAMM4686=AEFF + HCLL

end function gamm4686

real(kind=dp) function GAMM6683(TE,NE)
!     This function determines the value of Log10 (gamm(HeII6683))
!     = Log10( 4*Pai*j(HeII 6683)/N(He++)Ne) at temperature Te and densi
!     Storey P. J., Hummer D. G., 1995, MNRAS, 272, 41

      implicit none
      integer, parameter :: dp = kind(1.d0)

      real(kind=dp) TE,NE,AE2,AE3,AE4,AE5,AE6,AE7,AE8,AEFF,HCLL
      real(kind=dp) LNE, LTE

!debugging
#ifdef CO
        !print *,"function: GAMM6683"
#endif

      HCLL=-11.52688452 ! = Log10 ( h * c / lambda(6683) ) - [ cgs units
      LNE = log10(NE)
      LTE = log10(TE)

      AE2 =(-8.19880E+00)+                                              &
     &     (-4.85499E+00)* LTE +                                        &
     &       1.80862E+00 * LTE ** 2 +                                   &
     &     (-3.27542E-01)* LTE ** 3 +                                   &
     &       2.01978E-02 * LTE ** 4

      AE3 =(-7.36105E+00)+                                              &
     &     (-5.57358E+00)* LTE +                                        &
     &       2.04065E+00 * LTE ** 2 +                                   &
     &     (-3.60960E-01)* LTE ** 3 +                                   &
     &       2.20083E-02 * LTE ** 4

      AE4 =(-5.99396E+00)+                                              &
     &     (-6.75645E+00)* LTE +                                        &
     &       2.42722E+00 * LTE ** 2 +                                   &
     &     (-4.17511E-01)* LTE ** 3 +                                   &
     &       2.51313E-02 * LTE ** 4

      AE5 =(-4.04880E+00)+                                              &
     &     (-8.39805E+00)* LTE +                                        &
     &       2.94942E+00 * LTE ** 2 +                                   &
     &     (-4.91759E-01)* LTE ** 3 +                                   &
     &       2.91127E-02 * LTE ** 4

      AE6 =(-1.41066E+00)+                                              &
     &     (-1.05026E+01)* LTE +                                        &
     &       3.58120E+00 * LTE ** 2 +                                   &
     &     (-5.76551E-01)* LTE ** 3 +                                   &
     &       3.34122E-02 * LTE ** 4

      AE7 =  1.98527E+00 +                                              &
     &     (-1.30886E+01)* LTE +                                        &
     &       4.33145E+00 * LTE ** 2 +                                   &
     &     (-6.75637E-01)* LTE ** 3 +                                   &
     &       3.84651E-02 * LTE ** 4

      AE8 =  7.37380E+00 +                                              &
     &     (-1.74412E+01)* LTE +                                        &
     &       5.70033E+00 * LTE ** 2 +                                   &
     &     (-8.74093E-01)* LTE ** 3 +                                   &
     &       4.96046E-02 * LTE ** 4


      if (LNE.lt.2) then
        AEFF=AE2
      elseif (LNE.ge.2..AND.LNE.lt.3) then
        AEFF=AE2 + (AE3 - AE2) * (LNE - 2)
      elseif (LNE.ge.3..AND.LNE.lt.4) then
        AEFF=AE3 + (AE4 - AE3) * (LNE - 3)
      elseif (LNE.ge.4..AND.LNE.lt.5) then
        AEFF=AE4 + (AE5 - AE4) * (LNE - 4)
      elseif (LNE.ge.5..AND.LNE.lt.6) then
        AEFF=AE5 + (AE6 - AE5) * (LNE - 5)
      elseif (LNE.ge.6..AND.LNE.lt.7) then
        AEFF=AE6 + (AE7 - AE6) * (LNE - 6)
      elseif (LNE.ge.7..AND.LNE.lt.8) then
        AEFF=AE7 + (AE8 - AE7) * (LNE - 7)
      else
        AEFF=AE8
      endif

      GAMM6683=AEFF + HCLL

end function gamm6683

subroutine read_heii
!get Case B emissivities from file, values are from Storey and Hummer 1995.
implicit none
character(len=1) :: junk
character(len=20), dimension(8) :: invar
integer :: i,j,k,l !counters

!debugging
#ifdef CO
        print *,"subroutine: read_heii"
#endif

open (unit=357, file=trim(PREFIX)//"/share/neat/RHeii.dat")
read (357,*) junk !first line is a comment
read (357,"(I3,I3)") ntemps, ndens !second line has number of temperatures and densities
nlevs=25 ! maximum number of levels shown in intrat data file

!allocate the emissivities array, dimensions are temperature, density, level 1, level 2

allocate(heiidata(ntemps, ndens, nlevs, nlevs))
allocate(temperatures(ntemps))

heiidata(:,:,:,:)=0.d0

!now, loop through the temperatures and densities and read in the emissivities

do i=1,ntemps
  do j=1,ndens
    read (357,*) invar(1:6) !line at top of each block has density, charge, temperature, case, maximum level calculated and maximum level displayed
    if (j.eq.1) then
      read(invar(3),"(E9.2)") temperatures(i) ! get the temperatures for later use in interpolating
    endif
    do k=nlevs,1,-1
      do l=1,k-1
        read (357,"(E10.3)", advance="no") heiidata(i,j,k,l)
      enddo
    enddo
  enddo
enddo
close (357)

!convert to logarithms, makes it easier for the abundance calculation
where (heiidata .ne. 0.d0)
  heiidata=log10(heiidata)
endwhere

end subroutine read_heii

subroutine get_heii_abund_new(linelist,HeII_lines,medtemp,density,heiiabund)!,weights)
!get an array of He II emissivities interpolated to the correct temperature and density
!then use that to calculate the abundance from all the detected HeII lines
!then get the final abundance from a weighted average of individual line abundances (default is 4686 has weight of 1, all others are zero)
implicit none

real(kind=dp), dimension(:,:), allocatable :: emissivityarray
real(kind=dp) :: medtemp, density
real(kind=dp) :: x1,x2,y1,y2,y,factor
real(kind=dp) :: abundsum, weightsum, heiiabund
type(line), dimension(:) :: linelist
integer, dimension(20,2:6) :: HeII_lines
integer :: i,j,d
integer :: t1, t2, d1, d2 !locations on temperature and density axes of heiidata between which interpolation will take place

!debugging
#ifdef CO
        print *,"subroutine: get_heii_emissivities"
#endif

!initialise

heiiabund = 0.d0

!allocate the search array.  emissivities array has dimensions of temperature, density, level 1, level 2
!search array just has dimensions of level 1, level 2

allocate(emissivityarray(size(heiidata,3),size(heiidata,4)))

!for each transition, first interpolate in temperature then interpolate in density to get the emissivity and put it in the result array.
!get the box it's in in te vs. ne grid:

do i=1,ntemps-1
  if (medtemp .ge. temperatures(i) .and. medtemp .le. temperatures(i+1)) then
    exit
  endif
enddo

d=floor(log10(density))

! bilinear interpolation of the logarithm of the emissivity.
! degrades to linear or no interpolation if the temperature and density are outside the limits of the data

  if (i .lt. 1) then
    t1=1
    t2=1
  elseif (i .eq. size(temperatures)) then
    t1=size(temperatures)
    t2=size(temperatures)
  else
    t1=i
    t2=i+1
  endif
! for density, array index is from 1 to 7, containing densities from 10**2 to 10**8. todo: change the array indexing to 2:8
  if (d .lt. 2) then
    d1=1
    d2=1
  elseif (d .ge. 8) then
    d1=7
    d2=7
  else
    d1=d-1
    d2=d
  endif

  x1=temperatures(t1)
  x2=temperatures(t2)
  y1=real(d1+1)
  y2=real(d2+1)
  y=log10(density)

  if (t1.ne.t2 .and. d1.ne.d2) then ! interpolate bilinearly in temperature and density
    factor = 1.d0/((x2-x1)*(y2-y1))
    emissivityarray = (((x2-medtemp)*(y2-y)) * heiidata(t1,d1,:,:) + &
                   & ((medtemp-x1)*(y2-y)) * heiidata(t2,d1,:,:) + &
                   & ((x2-medtemp)*(y-y1)) * heiidata(t1,d2,:,:) + &
                   & ((medtemp-x1)*(y-y1)) * heiidata(t2,d2,:,:)) * factor
  elseif (t1.eq.t2 .and. d1.ne.d2) then ! interpolate in density only
    emissivityarray = heiidata(t1,d1,:,:) + (heiidata(t1,d2,:,:) - heiidata(t1,d1,:,:))*(y-y1)/(y2-y1)
  elseif (t1.ne.t2 .and. d1.eq.d2) then ! interpolate in temperature only
    emissivityarray = heiidata(t1,d1,:,:) + (heiidata(t2,d1,:,:) - heiidata(t1,d1,:,:))*(medtemp-x1)/(x2-x1)
  else !no interpolation. todo: warn about being outside the limits.
    emissivityarray = heiidata(t1,d1,:,:)
  endif

!now to get the abundances, go through the helium line location array and for all lines present, get their abundance.
!then get weighted overall abundance
abundsum = 0.d0
weightsum = 0.d0

do i=2,6
  do j=1,20
    if (HeII_lines(j,i).gt.0) then !the line is detected
      if (linelist(HeII_lines(j,i))%int_dered .gt. 0.d0) then !and has non zero flux so calculate abundance
        linelist(HeII_lines(j,i))%abundance = linelist(HeII_lines(j,i))%int_dered/100.0 * 10.**(gamm4861(medtemp,density)-emissivityarray(j,i))
        abundsum = abundsum + linelist(HeII_lines(j,i))%weight*linelist(HeII_lines(j,i))%abundance
        weightsum = weightsum + linelist(HeII_lines(j,i))%weight
!!print *,linelist(HeII_lines(j,i))%wavelength,linelist(HeII_lines(j,i))%abundance
      endif
    endif
  enddo
enddo

if (weightsum .gt. 0.d0) then
  heiiabund = abundsum / weightsum
endif

end subroutine get_heii_abund_new

end module mod_helium
