      module mod_helium
      contains

      subroutine get_heii_abund(te, ne, IHeII4686, heiiabund)

      IMPLICIT NONE
      REAL A4686,TE,NE,IHeII4686,heiiabund

      A4686=10.**(GAMM4861(TE,NE)-GAMM4686(TE,NE))
      heiiabund=IHeII4686/100.*A4686

      end subroutine

      subroutine get_hei_porter(Te, ne, he_lines, heidata, Heiabund)
      use mod_abundtypes

      IMPLICIT NONE
      real :: te, ne
      real :: AB4471,AB5876,AB6678,heiabund 

      type(line), dimension(44) :: he_lines 
      double precision, dimension(3) :: weights
      double precision, dimension(21,15,44), intent(in) :: heidata
      real, dimension(44) :: emissivities
      integer :: i

!      data is for the following lines, in this order: 2945.10,3188.74,3613.64,3888.65,3964.73,4026.21,4120.82,4387.93,4437.55,4471.50,4713.17,4921.93,5015.68,5047.74,5875.66,6678.16,7065.25,7281.35,9463.58,10830.25,11013.07,11969.06,12527.49,12755.69,12784.92,12790.50,12845.98,12968.43,12984.88,13411.69,15083.65,17002.40,18555.57,18685.33,18697.21,19089.36,19543.19,20424.97,20581.28,20601.76,21120.12,21132.03,21607.80,21617.01 /)


      do i = 1,44
        call get_emissivity_porter(te,ne, i, emissivities(i), heidata)
        He_lines(i)%abundance = He_lines(i)%int_dered/100. * 10.**(GAMM4861(TE,NE)-emissivities(i))
      end do

      AB4471 = He_lines(10)%abundance
      AB5876 = He_lines(15)%abundance
      AB6678 = He_lines(16)%abundance

           weights=0.D0

           if (AB4471 .gt. 0) weights(1) = 1.
           if (AB5876 .gt. 0) weights(2) = 3.
           if (AB6678 .gt. 0) weights(3) = 1.

           if (weights(1)+weights(2)+weights(3).gt.0.0) then
               heiabund = ((weights(1)*AB4471) + (weights(2)*AB5876) + (weights(3)*AB6678)) / (weights(1) + weights(2) + weights(3))
           else
               heiabund = 0.D0
           endif

      end subroutine get_hei_porter

        subroutine get_emissivity_porter(te, ne, line, emissivity, heidata)

        implicit none
        real :: te, ne, testart, nestart, logne
        real :: interp_factor_te, interp_factor_ne, interp_t1, interp_t2, emissivity
        double precision, dimension(21,15,44), intent(in) :: heidata
        integer :: i,j, line

        ! ne needs to be logarithmic, input is linear

        logne=log10(ne)

        ! check that data is within the ranges calculated by Porter

        if (te .lt. 5000) te=5000.
        if (te .gt. 25000) te=25000.
        if (logne .lt. 1.) logne=1.
        if (logne .gt. 14.) logne=14.

        !do bilinear interpolation of the log values
        !find starting values

        testart=25000.
        nestart=15.

        !find temperature box
        do i=1,21
          if (te .lt. dble(i*1000+4000)) then
            testart=dble((i-1)*1000+4000)
            exit
          endif
        end do

        !find density box
        do j=1,14
          if (logne .lt. dble(j)) then
            nestart=dble(j-1)
            exit
          endif
        end do

        i=i-1
        j=j-1

        !now we have the array positions in i and j and the values in nestart and
        !testart

        !now interpolate first in T and then in ne

        interp_factor_te = (te-testart)/1000.
        interp_factor_ne = (logne-nestart)

        !get two interpolated values for each ne, then interpolate between those

        interp_t1 = heidata(i,j,line)+interp_factor_te*(heidata(i+1,j,line)-heidata(i,j,line))
        interp_t2 = heidata(i,j+1,line)+interp_factor_te*(heidata(i+1,j+1,line)-heidata(i,j+1,line))

        emissivity = interp_t1+(interp_factor_ne*(interp_t2-interp_t1))

      end subroutine get_emissivity_porter

      subroutine get_hei_smits_new(Te, ne, he_lines, heidata, Heiabund)
      use mod_abundtypes

      IMPLICIT NONE
      real :: te, ne, tereduced
      real :: AB4471,AB5876,AB6678,heiabund 
      real :: c4471, c5876, c6678, d ! corrections for collisional excitation

      type(line), dimension(44) :: he_lines 
      double precision, dimension(3) :: weights
      double precision, dimension(3,6,44), intent(in) :: heidata
      real, dimension(44,3) :: emissivities
      real :: interpolatedemissivity
      integer :: i,j

!      data is for the following lines, in this order: 2945.10,3188.74,3613.64,3888.65,3964.73,4026.21,4120.82,4387.93,4437.55,4471.50,4713.17,4921.93,5015.68,5047.74,5875.66,6678.16,7065.25,7281.35,9463.58,10830.25,11013.07,11969.06,12527.49,12755.69,12784.92,12790.50,12845.98,12968.43,12984.88,13411.69,15083.65,17002.40,18555.57,18685.33,18697.21,19089.36,19543.19,20424.97,20581.28,20601.76,21120.12,21132.03,21607.80,21617.01 /)

      do i = 1,44
        do j = 1,3
          emissivities(i,j)=heidata(j,1,i) + heidata(j,2,i)*log10(te) + heidata(j,3,i)*log10(te)**2 + heidata(j,4,i)*log10(te)**3 + heidata(j,5,i)*log10(te)**4 + heidata(j,6,i)
        end do

        if (log10(ne).lt.2) then
          interpolatedemissivity=emissivities(i,1)
        elseif (log10(ne).ge.2..and.log10(ne).lt.4) then
          interpolatedemissivity=emissivities(i,1) + (emissivities(i,2) - emissivities(i,1)) * (log10(ne) - 2)
        elseif (log10(ne).ge.4..and.log10(ne).lt.6) then
          interpolatedemissivity=emissivities(i,2) + (emissivities(i,3) - emissivities(i,2)) * (log10(ne) - 4)
        else
          interpolatedemissivity=emissivities(i,3)
        endif

        if (interpolatedemissivity.ne. 0.D0) then
          He_lines(i)%abundance = He_lines(i)%int_dered/100. * 10.**(GAMM4861(TE,NE)-interpolatedemissivity)
        else
          He_lines(i)%abundance = 0.D0
        endif

      end do

      tereduced = te/10000.

!correct 4471, 5876 and 6678 for collisional contributions

      D=1.+3130.*tereduced**(-0.50)/ne
      C4471=(  6.95*tereduced**(0.15)*exp(-4.545/tereduced)  + 0.22*tereduced**(-0.55)*exp(-4.884/tereduced)  )/D
      C5876=(  6.78*tereduced**(0.07)*exp(-3.776/tereduced)  + 1.67*tereduced**(-0.15)*exp(-4.545/tereduced) + 0.60*tereduced**(-0.34)*exp(-4.901/tereduced)  )/D
      C6678=(  3.15*tereduced**(-0.54)*exp(-3.776/tereduced) + 0.51*tereduced**(-0.51)*exp(-4.545/tereduced) + 0.20*tereduced**(-0.66)*exp(-4.901/tereduced)  )/D

      He_lines(10)%abundance=He_lines(10)%abundance/(1.+C4471)
      He_lines(15)%abundance=He_lines(15)%abundance/(1.+C5876)
      He_lines(16)%abundance=He_lines(16)%abundance/(1.+C6678)

      AB4471 = He_lines(10)%abundance
      AB5876 = He_lines(15)%abundance
      AB6678 = He_lines(16)%abundance

           weights=0.D0

           if (AB4471 .gt. 0) weights(1) = 1.
           if (AB5876 .gt. 0) weights(2) = 3.
           if (AB6678 .gt. 0) weights(3) = 1.

           if (weights(1)+weights(2)+weights(3).gt.0.0) then
               heiabund = ((weights(1)*AB4471) + (weights(2)*AB5876) + (weights(3)*AB6678)) / (weights(1) + weights(2) + weights(3))
           else
               heiabund = 0.D0
           endif

      end subroutine get_hei_smits_new

      REAL FUNCTION GAMM4861(TE,NE)
!     This function determines the value of Log10 (gamm(H Beta))
!     = Log10( 4*Pai*j(HBeta)/NpNe) at temperature Te and density Ne
!     Storey P. J., Hummer D. G., 1995, MNRAS, 272, 41
!
      IMPLICIT NONE
      REAL TE,NE,AE2,AE3,AE4,AE5,AE6,AE7,AE8,AEFF,HCLL
      REAL LNE, LTE

                     ! = Log10 ( h * c / lambda(H-beta) ) - [ cgs units
      HCLL=-11.38871
      LNE = ALOG10(NE)
      LTE = ALOG10(TE)
!
      AE2 =(-9.06524E+00)+                                              &
     &     (-2.69954E+00)* LTE +                                        &
     &       8.80123E-01 * LTE ** 2 +                                   &
     &     (-1.57946E-01)* LTE ** 3 +                                   &
     &       9.25920E-03 * LTE ** 4
!
      AE3 =(-8.13757E+00)+                                              &
     &     (-3.57392E+00)* LTE +                                        &
     &       1.19331E+00 * LTE ** 2 +                                   &
     &     (-2.08362E-01)* LTE ** 3 +                                   &
     &       1.23303E-02 * LTE ** 4
!
      AE4 =(-6.87230E+00)+                                              &
     &     (-4.72312E+00)* LTE +                                        &
     &       1.58890E+00 * LTE ** 2 +                                   &
     &     (-2.69447E-01)* LTE ** 3 +                                   &
     &       1.58955E-02 * LTE ** 4
!
      AE5 =(-5.15059E+00)+                                              &
     &     (-6.24549E+00)* LTE +                                        &
     &       2.09801E+00 * LTE ** 2 +                                   &
     &     (-3.45649E-01)* LTE ** 3 +                                   &
     &       2.01962E-02 * LTE ** 4
!
      AE6 =(-2.35923E+00)+                                              &
     &     (-8.75565E+00)* LTE +                                        &
     &       2.95600E+00 * LTE ** 2 +                                   &
     &     (-4.77584E-01)* LTE ** 3 +                                   &
     &       2.78852E-02 * LTE ** 4
!
      AE7 =  1.55373E+00 +                                              &
     &     (-1.21894E+01)* LTE +                                        &
     &       4.10096E+00 * LTE ** 2 +                                   &
     &     (-6.49318E-01)* LTE ** 3 +                                   &
     &       3.76487E-02 * LTE ** 4
!
      AE8 =  6.59883E+00 +                                              &
     &     (-1.64030E+01)* LTE +                                        &
     &       5.43844E+00 * LTE ** 2 +                                   &
     &     (-8.40253E-01)* LTE ** 3 +                                   &
     &       4.79786E-02 * LTE ** 4
!
!
      IF (LNE.LT.2) THEN
        AEFF=AE2
      ELSEIF (LNE.GE.2..AND.LNE.LT.3) THEN
        AEFF=AE2 + (AE3 - AE2) * (LNE - 2)
      ELSEIF (LNE.GE.3..AND.LNE.LT.4) THEN
        AEFF=AE3 + (AE4 - AE3) * (LNE - 3)
      ELSEIF (LNE.GE.4..AND.LNE.LT.5) THEN
        AEFF=AE4 + (AE5 - AE4) * (LNE - 4)
      ELSEIF (LNE.GE.5..AND.LNE.LT.6) THEN
        AEFF=AE5 + (AE6 - AE5) * (LNE - 5)
      ELSEIF (LNE.GE.6..AND.LNE.LT.7) THEN
        AEFF=AE6 + (AE7 - AE6) * (LNE - 6)
      ELSEIF (LNE.GE.7..AND.LNE.LT.8) THEN
        AEFF=AE7 + (AE8 - AE7) * (LNE - 7)
      ELSE
        AEFF=AE8
      ENDIF
!
      GAMM4861=AEFF + HCLL
      END function gamm4861
!
!
      REAL FUNCTION GAMM4686(TE,NE)
!     This function determines the value of Log10 (gamm(HeII4686))
!     = Log10( 4*Pai*j(HeII 4686)/N(He++)Ne) at temperature Te and densi
!     Storey P. J., Hummer D. G., 1995, MNRAS, 272, 41
!
      IMPLICIT NONE
      REAL TE,NE,AE2,AE3,AE4,AE5,AE6,AE7,AE8,AEFF,HCLL
      REAL LNE, LTE

                     ! = Log10 ( h * c / lambda(4686) ) - [ cgs units ]
      HCLL=-11.37272
      LNE = ALOG10(NE)
      LTE = ALOG10(TE)
!
      AE2 =(-8.14491E+00)+                                              &
     &     (-2.10904E+00)* LTE +                                        &
     &       6.31303E-01 * LTE ** 2 +                                   &
     &     (-1.25035E-01)* LTE ** 3 +                                   &
     &       7.95257E-03 * LTE ** 4
!
      AE3 =(-7.41098E+00)+                                              &
     &     (-2.82798E+00)* LTE +                                        &
     &       8.90834E-01 * LTE ** 2 +                                   &
     &     (-1.66196E-01)* LTE ** 3 +                                   &
     &       1.03803E-02 * LTE ** 4
!
      AE4 =(-6.20402E+00)+                                              &
     &     (-3.99470E+00)* LTE +                                        &
     &       1.30756E+00 * LTE ** 2 +                                   &
     &     (-2.31744E-01)* LTE ** 3 +                                   &
     &       1.42229E-02 * LTE ** 4
!
      AE5 =(-4.11426E+00)+                                              &
     &     (-5.97955E+00)* LTE +                                        &
     &       2.00560E+00 * LTE ** 2 +                                   &
     &     (-3.39997E-01)* LTE ** 3 +                                   &
     &       2.04854E-02 * LTE ** 4
!
      AE6 =(-7.88408E-01)+                                              &
     &     (-9.06439E+00)* LTE +                                        &
     &       3.06790E+00 * LTE ** 2 +                                   &
     &     (-5.01667E-01)* LTE ** 3 +                                   &
     &       2.96826E-02 * LTE ** 4
!
      AE7 =  4.22286E+00 +                                              &
     &     (-1.35870E+01)* LTE +                                        &
     &       4.58664E+00 * LTE ** 2 +                                   &
     &     (-7.27438E-01)* LTE ** 3 +                                   &
     &       4.22467E-02 * LTE ** 4
!
      AE8 =  1.10404E+01 +                                              &
     &     (-1.94998E+01)* LTE +                                        &
     &       6.49776E+00 * LTE ** 2 +                                   &
     &     (-1.00125E+00)* LTE ** 3 +                                   &
     &       5.69519E-02 * LTE ** 4
!
!
      IF (LNE.LT.2) THEN
        AEFF=AE2
      ELSEIF (LNE.GE.2..AND.LNE.LT.3) THEN
        AEFF=AE2 + (AE3 - AE2) * (LNE - 2)
      ELSEIF (LNE.GE.3..AND.LNE.LT.4) THEN
        AEFF=AE3 + (AE4 - AE3) * (LNE - 3)
      ELSEIF (LNE.GE.4..AND.LNE.LT.5) THEN
        AEFF=AE4 + (AE5 - AE4) * (LNE - 4)
      ELSEIF (LNE.GE.5..AND.LNE.LT.6) THEN
        AEFF=AE5 + (AE6 - AE5) * (LNE - 5)
      ELSEIF (LNE.GE.6..AND.LNE.LT.7) THEN
        AEFF=AE6 + (AE7 - AE6) * (LNE - 6)
      ELSEIF (LNE.GE.7..AND.LNE.LT.8) THEN
        AEFF=AE7 + (AE8 - AE7) * (LNE - 7)
      ELSE
        AEFF=AE8
      ENDIF
!
      GAMM4686=AEFF + HCLL
      END function gamm4686
!
!
      REAL FUNCTION GAMM6683(TE,NE)
!     This function determines the value of Log10 (gamm(HeII6683))
!     = Log10( 4*Pai*j(HeII 6683)/N(He++)Ne) at temperature Te and densi
!     Storey P. J., Hummer D. G., 1995, MNRAS, 272, 41
!
      IMPLICIT NONE
      REAL TE,NE,AE2,AE3,AE4,AE5,AE6,AE7,AE8,AEFF,HCLL
      REAL LNE, LTE

                        ! = Log10 ( h * c / lambda(6683) ) - [ cgs units
      HCLL=-11.52688452
      LNE = ALOG10(NE)
      LTE = ALOG10(TE)
!
      AE2 =(-8.19880E+00)+                                              &
     &     (-4.85499E+00)* LTE +                                        &
     &       1.80862E+00 * LTE ** 2 +                                   &
     &     (-3.27542E-01)* LTE ** 3 +                                   &
     &       2.01978E-02 * LTE ** 4
!
      AE3 =(-7.36105E+00)+                                              &
     &     (-5.57358E+00)* LTE +                                        &
     &       2.04065E+00 * LTE ** 2 +                                   &
     &     (-3.60960E-01)* LTE ** 3 +                                   &
     &       2.20083E-02 * LTE ** 4
!
      AE4 =(-5.99396E+00)+                                              &
     &     (-6.75645E+00)* LTE +                                        &
     &       2.42722E+00 * LTE ** 2 +                                   &
     &     (-4.17511E-01)* LTE ** 3 +                                   &
     &       2.51313E-02 * LTE ** 4
!
      AE5 =(-4.04880E+00)+                                              &
     &     (-8.39805E+00)* LTE +                                        &
     &       2.94942E+00 * LTE ** 2 +                                   &
     &     (-4.91759E-01)* LTE ** 3 +                                   &
     &       2.91127E-02 * LTE ** 4
!
      AE6 =(-1.41066E+00)+                                              &
     &     (-1.05026E+01)* LTE +                                        &
     &       3.58120E+00 * LTE ** 2 +                                   &
     &     (-5.76551E-01)* LTE ** 3 +                                   &
     &       3.34122E-02 * LTE ** 4
!
      AE7 =  1.98527E+00 +                                              &
     &     (-1.30886E+01)* LTE +                                        &
     &       4.33145E+00 * LTE ** 2 +                                   &
     &     (-6.75637E-01)* LTE ** 3 +                                   &
     &       3.84651E-02 * LTE ** 4
!
      AE8 =  7.37380E+00 +                                              &
     &     (-1.74412E+01)* LTE +                                        &
     &       5.70033E+00 * LTE ** 2 +                                   &
     &     (-8.74093E-01)* LTE ** 3 +                                   &
     &       4.96046E-02 * LTE ** 4
!
!
      IF (LNE.LT.2) THEN
        AEFF=AE2
      ELSEIF (LNE.GE.2..AND.LNE.LT.3) THEN
        AEFF=AE2 + (AE3 - AE2) * (LNE - 2)
      ELSEIF (LNE.GE.3..AND.LNE.LT.4) THEN
        AEFF=AE3 + (AE4 - AE3) * (LNE - 3)
      ELSEIF (LNE.GE.4..AND.LNE.LT.5) THEN
        AEFF=AE4 + (AE5 - AE4) * (LNE - 4)
      ELSEIF (LNE.GE.5..AND.LNE.LT.6) THEN
        AEFF=AE5 + (AE6 - AE5) * (LNE - 5)
      ELSEIF (LNE.GE.6..AND.LNE.LT.7) THEN
        AEFF=AE6 + (AE7 - AE6) * (LNE - 6)
      ELSEIF (LNE.GE.7..AND.LNE.LT.8) THEN
        AEFF=AE7 + (AE8 - AE7) * (LNE - 7)
      ELSE
        AEFF=AE8
      ENDIF
!
      GAMM6683=AEFF + HCLL
      END function gamm6683

      end module mod_helium
