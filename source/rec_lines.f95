! conversion of MIDAS script Roii.prg, written by XWL, to F90
! RW May 2009

      module mod_recombination_lines

      TYPE oiiRL
            CHARACTER(len=1) :: Hyb
            CHARACTER(len=1) :: n_E1
            CHARACTER(len=1) :: n_E1GA
            CHARACTER(len=1) :: n_E2
            CHARACTER(len=1) :: n_E2GA
            CHARACTER(len=1) :: n_g1
            CHARACTER(len=1) :: n_g2
            CHARACTER(len=1) :: Rem1
            CHARACTER(len=1) :: Rem2
            CHARACTER(len=1) :: Rem3
            CHARACTER(len=1) :: Rem4
            CHARACTER(len=3) :: q_gf1
            CHARACTER(len=3) :: q_gf2
            CHARACTER(len=7) :: Mult
            CHARACTER(len=9) :: Term1
            CHARACTER(len=9) :: Term2
            INTEGER :: g1
            INTEGER :: g2
            INTEGER :: ION
            DOUBLE PRECISION :: Wave
            DOUBLE PRECISION :: E1
            DOUBLE PRECISION :: E2
            DOUBLE PRECISION :: Em
            DOUBLE PRECISION :: Int
            DOUBLE PRECISION :: Br_A
            DOUBLE PRECISION :: Br_B
            DOUBLE PRECISION :: Br_C
            DOUBLE PRECISION :: gf1
            DOUBLE PRECISION :: gf2
            DOUBLE PRECISION :: Obs
            DOUBLE PRECISION :: abundance
      END TYPE

      TYPE(oiiRL), DIMENSION(415) :: oiiRLs

      TYPE niiRL
            CHARACTER(len=1) :: Hyb
            CHARACTER(len=1) :: n_E1
            CHARACTER(len=1) :: n_E1GA
            CHARACTER(len=1) :: n_E2
            CHARACTER(len=1) :: n_E2GA
            CHARACTER(len=1) :: n_g1
            CHARACTER(len=1) :: n_g2
            CHARACTER(len=1) :: Rem1
            CHARACTER(len=1) :: Rem2
            CHARACTER(len=1) :: Rem3
            CHARACTER(len=1) :: Rem4
            CHARACTER(len=3) :: q_gf1
            CHARACTER(len=3) :: q_gf2
            CHARACTER(len=7) :: Mult
            CHARACTER(len=9) :: Term1
            CHARACTER(len=9) :: Term2
            INTEGER :: g1
            INTEGER :: g2
            INTEGER :: ION
            DOUBLE PRECISION :: Wave
            DOUBLE PRECISION :: E1
            DOUBLE PRECISION :: E2
            DOUBLE PRECISION :: Em
            DOUBLE PRECISION :: Int
            DOUBLE PRECISION :: Br_LS
            DOUBLE PRECISION :: gf1
            DOUBLE PRECISION :: gf2
            DOUBLE PRECISION :: Obs
            DOUBLE PRECISION :: abundance
      END TYPE

      TYPE(niiRL), DIMENSION(99) :: niiRLs

      TYPE ciiRL
            DOUBLE PRECISION :: Wave
            DOUBLE PRECISION :: a
            DOUBLE PRECISION :: b
            DOUBLE PRECISION :: c
            DOUBLE PRECISION :: d
            DOUBLE PRECISION :: f
            DOUBLE PRECISION :: aeff
            DOUBLE PRECISION :: Int
            DOUBLE PRECISION :: Obs
            DOUBLE PRECISION :: abundance
      END TYPE

      TYPE(ciiRL), DIMENSION(57) :: ciiRLs

      TYPE neiiRL
            DOUBLE PRECISION :: Wave
            DOUBLE PRECISION :: a
            DOUBLE PRECISION :: b
            DOUBLE PRECISION :: c
            DOUBLE PRECISION :: d
            DOUBLE PRECISION :: f
            DOUBLE PRECISION :: Br
            DOUBLE PRECISION :: aeff
            DOUBLE PRECISION :: Int
            DOUBLE PRECISION :: Obs
            DOUBLE PRECISION :: abundance
      END TYPE

      TYPE(neiiRL), DIMENSION(38) :: neiiRLs

      TYPE xiiiRL
            CHARACTER(len=3) :: Ion
            DOUBLE PRECISION :: Wave
            DOUBLE PRECISION :: a
            DOUBLE PRECISION :: b
            DOUBLE PRECISION :: c
            DOUBLE PRECISION :: d
            DOUBLE PRECISION :: Br
            DOUBLE PRECISION :: aeff
            DOUBLE PRECISION :: Int
            DOUBLE PRECISION :: Obs
            DOUBLE PRECISION :: abundance
      END TYPE

      TYPE(xiiiRL), DIMENSION(6) :: xiiiRLs

      contains

      subroutine read_orl_data

        IMPLICIT NONE
        integer :: i
        ! read in OII data

            301 FORMAT (I5, 1X, F9.4, 1X, A1, A1, A1, A1, A1, F7.4,     &
     & 1X, A3, 1X, F7.4, 1X, A3, 1X, A7, 3X, F11.4, A1, A1, 1X, I2, &
     &1X, A1, 1X, A9, 1X, F13.4, 1X, A1, A1, 1X, I2, 1X, A1, 1X, A9, 1X,&
     & F7.4, 1X, F7.4, 1X, F7.4)!, 1X, E10.4, 1X, E10.4, 1X)
            OPEN(201, file="Atomic-data/Roii.dat", status='old')
            DO i = 1,415
            READ(201,301) oiiRLs(i)%ION, oiiRLs(i)%Wave, oiiRLs(i)%Hyb, &
     &oiiRLs(i)%Rem1, oiiRLs(i)%Rem2, oiiRLs(i)%Rem3, oiiRLs(i)%Rem4,   &
     &oiiRLs(i)%gf1, oiiRLs(i)%q_gf1, oiiRLs(i)%gf2, oiiRLs(i)%q_gf2,   &
     &oiiRLs(i)%Mult, oiiRLs(i)%E1, oiiRLs(i)%n_E1, oiiRLs(i)%n_E1GA,   &
     &oiiRLs(i)%g1, oiiRLs(i)%n_g1, oiiRLs(i)%Term1, oiiRLs(i)%E2,      &
     &oiiRLs(i)%n_E2, oiiRLs(i)%n_E2GA, oiiRLs(i)%g2, oiiRLs(i)%n_g2,   &
     &oiiRLs(i)%Term2, oiiRLs(i)%Br_A, oiiRLs(i)%Br_B, oiiRLs(i)%Br_C
            END DO
      CLOSE(201)

                ! read in NII data

            302 FORMAT (I5, 1X, F9.4, 1X, A1, A1, A1, A1, A1, 1X, F7.4, &
     & 1X, A3, 1X, F7.4, 1X, A3, 1X, A7, 3X, F11.4, A1, A1, 1X, I2, &
     &1X, A1, 1X, A9, 1X, F13.4, 1X, A1, A1, 1X, I2, 1X, A1, 1X, A9, 1X,&
     & F7.4, 1X, F7.4, 1X, F7.4)!, 1X, E10.4, 1X, E10.4, 1X)
            OPEN(201, file="Atomic-data/Rnii.dat", status='old')
            DO i = 1,99
            READ(201,302) niiRLs(i)%ION, niiRLs(i)%Wave, niiRLs(i)%Hyb, &
     &niiRLs(i)%Rem1, niiRLs(i)%Rem2, niiRLs(i)%Rem3, niiRLs(i)%Rem4,   &
     &niiRLs(i)%gf1, niiRLs(i)%q_gf1, niiRLs(i)%gf2, niiRLs(i)%q_gf2,   &
     &niiRLs(i)%Mult, niiRLs(i)%E1, niiRLs(i)%n_E1, niiRLs(i)%n_E1GA,   &
     &niiRLs(i)%g1, niiRLs(i)%n_g1, niiRLs(i)%Term1, niiRLs(i)%E2,      &
     &niiRLs(i)%n_E2, niiRLs(i)%n_E2GA, niiRLs(i)%g2, niiRLs(i)%n_g2,   &
     &niiRLs(i)%Term2, niiRLs(i)%Br_LS
            END DO
      CLOSE(201)

! read in CII data

       303 FORMAT (F7.2, 1X, F6.4, 1X, F7.4, 1X, F7.4, 1X, F7.4, 1X, F7.4)
       OPEN(201, file="Atomic-data/Rcii.dat", status='old')
       DO i = 1,57
         READ(201,303) ciiRLs(i)%Wave, ciiRLs(i)%a, ciiRLs(i)%b, &
         & ciiRLs(i)%c, ciiRLs(i)%d, ciiRLs(i)%f
       END DO
       CLOSE(201)

           ! read in NeII data

       304 FORMAT (F7.2, 1X, F6.3, 1X, F6.3, 1X, F6.3, 1X, F6.3, 1X, F7.4, 1X, F6.3)
       OPEN(201, file="Atomic-data/Rneii.dat", status='old')
       DO i = 1,38
         READ(201,304) neiiRLs(i)%Wave, neiiRLs(i)%a, neiiRLs(i)%b, &
         & neiiRLs(i)%c, neiiRLs(i)%d, neiiRLs(i)%f, neiiRLs(i)%Br
       END DO
       CLOSE(201)

        ! read in XIII data

       305 FORMAT (A3,1X,F7.2, 1X, F5.3, 1X, F6.3, 1X, F5.3, 1X, F5.3, 1X, F5.4)
       OPEN(201, file="Atomic-data/Rxiii.dat", status='old')
       DO i = 1,6
         READ(201,305) xiiiRLs(i)%ion, xiiiRLs(i)%Wave, xiiiRLs(i)%a, &
         & xiiiRLs(i)%b, xiiiRLs(i)%c, xiiiRLs(i)%d, xiiiRLs(i)%Br
       END DO
       CLOSE(201)


      end subroutine

      subroutine oii_rec_lines(te,ne,abund,oiiRLs)

      IMPLICIT NONE
      DOUBLE PRECISION :: aeff, aeff_hb, Em_Hb, &
     & Te, Ne, abund
      DOUBLE PRECISION :: a, b, c, d, an(4)

      INTEGER :: i


      TYPE(oiiRL), DIMENSION(415) :: oiiRLs

      call get_aeff_hb(te,ne, aeff_hb, em_hb)

!! read in OII data
!
!            301 FORMAT (I5, 1X, F9.4, 1X, A1, A1, A1, A1, A1, F7.4,     &
!     & 1X, A3, 1X, F7.4, 1X, A3, 1X, A7, 3X, F11.4, A1, A1, 1X, I2, &
!     &1X, A1, 1X, A9, 1X, F13.4, 1X, A1, A1, 1X, I2, 1X, A1, 1X, A9, 1X,&
!     & F7.4, 1X, F7.4, 1X, F7.4)!, 1X, E10.4, 1X, E10.4, 1X)
!            OPEN(201, file="Atomic-data/Roii.dat", status='old')
!            DO i = 1,415
!            READ(201,301) oiiRLs(i)%ION, oiiRLs(i)%Wave, oiiRLs(i)%Hyb, &
!     &oiiRLs(i)%Rem1, oiiRLs(i)%Rem2, oiiRLs(i)%Rem3, oiiRLs(i)%Rem4,   &
!     &oiiRLs(i)%gf1, oiiRLs(i)%q_gf1, oiiRLs(i)%gf2, oiiRLs(i)%q_gf2,   &
!     &oiiRLs(i)%Mult, oiiRLs(i)%E1, oiiRLs(i)%n_E1, oiiRLs(i)%n_E1GA,   &
!     &oiiRLs(i)%g1, oiiRLs(i)%n_g1, oiiRLs(i)%Term1, oiiRLs(i)%E2,      &
!     &oiiRLs(i)%n_E2, oiiRLs(i)%n_E2GA, oiiRLs(i)%g2, oiiRLs(i)%n_g2,   &
!     &oiiRLs(i)%Term2, oiiRLs(i)%Br_A, oiiRLs(i)%Br_B, oiiRLs(i)%Br_C
!            END DO
!      CLOSE(201)

! 4f-3d transitions

      a = 0.232
      b = -0.92009
      c = 0.15526
      d = 0.03442
      an = (/0.236,0.232,0.228,0.222/)

      if (log10(ne) .le. 2) then
        a = an(1)
      elseif (log10(ne) .gt. 2 .and. log10(ne) .le. 4) then
        a = an(1) + (an(2) - an(1)) / 2. * (log10(ne) - 2.)
      elseif (log10(ne) .gt. 4 .and. log10(ne) .le. 5) then
        a = an(2) + (an(3) - an(2)) * (log10(ne) - 2.)
      elseif (log10(ne) .gt. 6 .and. log10(ne) .le. 6) then
        a = an(3) + (an(4) - an(3)) * (log10(ne) - 2.)
      else
        a = an(4)
      endif

      te = te/10000
      aeff = 1.e-14 * a * te ** (b)
      aeff = aeff * (1. + c * (1. - te) + d * (1. - te) ** 2)

! New eff. recomb. data for 1000 < T < 5000 K, Ne = 100/cm3

      if (te .le. 0.5) then
        a = 0.236
        b = -1.07552
        c = -0.04843
        aeff = 1.e-14 * a * te ** (b + c * log(te))
      endif

      do i = 1,183
        oiiRLs(i)%Em = aeff * 1.98648E-08 /oiiRLs(i)%Wave * &
        & oiiRLs(i)%g2 * oiiRLs(i)%Br_B
        oiiRLs(i)%Int = 100. * oiiRLs(i)%Em / Em_hb * abund
      enddo

! 3d-3p ^4F transitions (Case A=B=C for a,b,c,d; Br diff. slightly, adopt Case B)
      a =  0.876
      b =  -0.73465
      c =  0.13689
      d =  0.06220
      an = (/0.876,0.876,0.877,0.880/) !a for logNe = 2,4,5,6 (LSBC95, Tab.5a)
      if ( log(ne) .le. 2) then
            a = an(1)
      elseif ( log(ne) .gt. 2 .and. log(ne) .le. 4) then
            a = an(1) + (an(2) - an(1)) / 2.*(log(ne) - 2.)
      elseif ( log(ne) .gt. 4 .and. log(ne) .le. 5) then
            a = an(2) + (an(3) - an(2))*(log(ne) - 2.)
      elseif ( log(ne) .gt. 6 .and. log(ne) .le. 6) then
            a = an(3) + (an(4) - an(3))*(log(ne) - 2.)
      else
              a = an(4)
      endif
!
      aeff = 1.e-14*a*te ** (b)
      aeff = aeff*(1. + c*(1. - te) + d*(1. - te) ** 2)
!
! New eff. recomb. data for 1000 < T < 5000 K, Ne = 100/cm3
! c.f. ~/source/OII_a_eff_at_low_Te/00readme/
      if ( te .le. 0.5) then
            a = 0.878
            b = -0.86175
            c = -0.02470
            aeff = 1.e-14*a*te ** (b + c*log(te))
      endif
!
     do i = 184,219
        oiiRLs(i)%Em = aeff*1.98648E-08 / oiiRLs(i)%Wave*&
      & oiiRLs(i)%g2*oiiRLs(i)%Br_B
        oiiRLs(i)%Int = 100.*oiiRLs(i)%Em / Em_hb*abund
      enddo
!
! 3d-3p ^4D, ^4P transitions
      a =  0.745
      b =  -0.74621
      c =  0.15710
      d =  0.07059
!      an = (/0.727,0.726,0.725,0.726/) !a for logNe = 2,4,5,6 (LSBC95, Tab.5a) Case A
      an = (/0.747,0.745,0.744,0.745/) !a for logNe = 2,4,5,6 (LSBC95, Tab.5a) Case B
!      an = (/0.769,0.767,0.766,0.766/) !a for logNe = 2,4,5,6 (LSBC95, Tab.5a) Case C
      if ( log(ne) .le. 2) then
            a = an(1)
      elseif ( log(ne) .gt. 2 .and. log(ne) .le. 4) then
            a = an(1) + (an(2) - an(1)) / 2.*(log(ne) - 2.)
      elseif ( log(ne) .gt. 4 .and. log(ne) .le. 5) then
            a = an(2) + (an(3) - an(2))*(log(ne) - 2.)
      elseif ( log(ne) .gt. 6 .and. log(ne) .le. 6) then
            a = an(3) + (an(4) - an(3))*(log(ne) - 2.)
      else
              a = an(4)
      endif
!
      aeff = 1.e-14*a*te ** (b)
      aeff = aeff*(1. + c*(1. - te) + d*(1. - te) ** 2)
!
! New eff. recomb. data for 1000 < T < 5000 K, Ne = 100/cm3
! c.f. ~/source/OII_a_eff_at_low_Te/00readme/
      if ( te .le. 0.5) then
            a = 0.747
            b = -0.89382
            c = -0.02906
            aeff = 1.e-14*a*te ** (b + c*log(te))
      endif
!
     do i = 220,310
        oiiRLs(i)%Em = aeff*1.98648E-08 / oiiRLs(i)%Wave*&
      & oiiRLs(i)%g2*oiiRLs(i)%Br_B
        oiiRLs(i)%Int = 100.*oiiRLs(i)%Em / Em_hb*abund
      enddo
!
! 3d-3p ^2F transitions
      a =  0.745
      b =  -0.74621
      c =  0.15710
      d =  0.07059
      an = (/0.727,0.726,0.725,0.726/) !a for logNe = 2,4,5,6 (LSBC95, Tab.5a) Case A
!      an = (/0.747,0.745,0.744,0.745/) !a for logNe = 2,4,5,6 (LSBC95, Tab.5a) Case B
!      an = (/0.769,0.767,0.766,0.766/) !a for logNe = 2,4,5,6 (LSBC95, Tab.5a) Case C
      if ( log(ne) .le. 2) then
            a = an(1)
      elseif ( log(ne) .gt. 2 .and. log(ne) .le. 4) then
            a = an(1) + (an(2) - an(1)) / 2.*(log(ne) - 2.)
      elseif ( log(ne) .gt. 4 .and. log(ne) .le. 5) then
            a = an(2) + (an(3) - an(2))*(log(ne) - 2.)
      elseif ( log(ne) .gt. 6 .and. log(ne) .le. 6) then
            a = an(3) + (an(4) - an(3))*(log(ne) - 2.)
      else
              a = an(4)
      endif
!
      aeff = 1.e-14*a*te ** (b)
      aeff = aeff*(1. + c*(1. - te) + d*(1. - te) ** 2)
!
! New eff. recomb. data for 1000 < T < 5000 K, Ne = 100/cm3
! c.f. ~/source/OII_a_eff_at_low_Te/00readme/
      if ( te .le. 0.5) then
            a = 0.747
            b = -0.89382
            c = -0.02906
            aeff = 1.e-14*a*te ** (b + c*log(te))
      endif
!
     do i = 311,328
        oiiRLs(i)%Em = aeff*1.98648E-08 / oiiRLs(i)%Wave*&
      & oiiRLs(i)%g2*oiiRLs(i)%Br_A
        oiiRLs(i)%Int = 100.*oiiRLs(i)%Em / Em_hb*abund
      enddo
!
! 3d-3p ^2D transitions
      a =  0.601
      b =  -0.79533
      c =  0.15314
      d =  0.05322
      an = (/0.603,0.601,0.600,0.599/) !a for logNe = 2,4,5,6 (LSBC95, Tab.5a) Case A
!      an = (/0.620,0.618,0.616,0.615/) !a for logNe = 2,4,5,6 (LSBC95, Tab.5a) Case C
      if ( log(ne) .le. 2) then
            a = an(1)
      elseif ( log(ne) .gt. 2 .and. log(ne) .le. 4) then
            a = an(1) + (an(2) - an(1)) / 2.*(log(ne) - 2.)
      elseif ( log(ne) .gt. 4 .and. log(ne) .le. 5) then
            a = an(2) + (an(3) - an(2))*(log(ne) - 2.)
      elseif ( log(ne) .gt. 6 .and. log(ne) .le. 6) then
            a = an(3) + (an(4) - an(3))*(log(ne) - 2.)
      else
              a = an(4)
      endif
!
      aeff = 1.e-14*a*te ** (b)
      aeff = aeff*(1. + c*(1. - te) + d*(1. - te) ** 2)
!
! New eff. recomb. data for 1000 < T < 5000 K, Ne = 100/cm3
! c.f. ~/source/OII_a_eff_at_low_Te/00readme/
      if ( te .le. 0.5) then
            a = 0.603
            b = -0.94025
            c = -0.03467
            aeff = 1.e-14*a*te ** (b + c*log(te))
      endif
!
     do i = 329,358
        oiiRLs(i)%Em = aeff*1.98648E-08 / oiiRLs(i)%Wave*&
      & oiiRLs(i)%g2*oiiRLs(i)%Br_A
        oiiRLs(i)%Int = 100.*oiiRLs(i)%Em / Em_hb*abund
      enddo
!
! 3d-3p ^2P transitions
      a =  0.524
      b =  -0.78448
      c =  0.13681
      d =  0.05608
      an = (/0.526,0.524,0.523,0.524/) !a for logNe = 2,4,5,6 (LSBC95, Tab.5a) Case A
!      an = (/0.538,0.536,0.535,0.536/) !a for logNe = 2,4,5,6 (LSBC95, Tab.5a) Case C
      if ( log(ne) .le. 2) then
            a = an(1)
      elseif ( log(ne) .gt. 2 .and. log(ne) .le. 4) then
            a = an(1) + (an(2) - an(1)) / 2.*(log(ne) - 2.)
      elseif ( log(ne) .gt. 4 .and. log(ne) .le. 5) then
            a = an(2) + (an(3) - an(2))*(log(ne) - 2.)
      elseif ( log(ne) .gt. 6 .and. log(ne) .le. 6) then
            a = an(3) + (an(4) - an(3))*(log(ne) - 2.)
      else
              a = an(4)
      endif
!
      aeff = 1.e-14*a*te ** (b)
      aeff = aeff*(1. + c*(1. - te) + d*(1. - te) ** 2)
!
! New eff. recomb. data for 1000 < T < 5000 K, Ne = 100/cm3
! c.f. ~/source/OII_a_eff_at_low_Te/00readme/
      if ( te .le. 0.5) then
            a = 0.526
            b = -0.91758
            c = -0.03120
            aeff = 1.e-14*a*te ** (b + c*log(te))
      endif
!
     do i = 359,388
        oiiRLs(i)%Em = aeff*1.98648E-08 / oiiRLs(i)%Wave*&
      & oiiRLs(i)%g2*oiiRLs(i)%Br_A
        oiiRLs(i)%Int = 100.*oiiRLs(i)%Em / Em_hb*abund
      enddo
!
! 3p-3s ^4D - ^4P transitions
!      an = (/34.7,34.9,35.1,35.0/) !a for logNe = 2,4,5,6 Case A
!      a =  36.2
!      b =  -0.749
!      c =  0.023
!      d =  0.074
      an = (/36.0,36.2,36.4,36.3/) !a for logNe = 2,4,5,6 Case B
      a =  36.2
      b =  -0.736
      c =  0.033
      d =  0.077
      if ( log(ne) .le. 2) then
            a = an(1)
      elseif ( log(ne) .gt. 2 .and. log(ne) .le. 4) then
            a = an(1) + (an(2) - an(1)) / 2.*(log(ne) - 2.)
      elseif ( log(ne) .gt. 4 .and. log(ne) .le. 5) then
            a = an(2) + (an(3) - an(2))*(log(ne) - 2.)
      elseif ( log(ne) .gt. 6 .and. log(ne) .le. 6) then
            a = an(3) + (an(4) - an(3))*(log(ne) - 2.)
      else
              a = an(4)
      endif
!
      aeff = 1.e-14*a*te ** (b)
      aeff = aeff*(1. + c*(1. - te) + d*(1. - te) ** 2)
!
! New eff. recomb. data for 1000 < T < 5000 K, Ne = 100/cm3
! c.f. ~/source/OII_a_eff_at_low_Te/00readme/
      if ( te .le. 0.5) then
            a = 36.288
            b = -0.75421
            c = 0.02883
            d = 0.01213
            aeff = 1.e-14*a*te ** (b + c*log(te) + d* log(te) ** 2)
      endif
!
     do i = 389,396
        oiiRLs(i)%Em = aeff*1.98648E-08 / oiiRLs(i)%Wave*&
      & oiiRLs(i)%g2*oiiRLs(i)%Br_B
        oiiRLs(i)%Int = 100.*oiiRLs(i)%Em / Em_hb*abund
      enddo
!
! 3p-3s ^4P - ^4P transitions
!      an = (/10.4,10.4,10.5,10.4/) !a for logNe = 2,4,5,6 Case A
!      a =  10.4
!      b =  -0.721
!      c =  0.073
!      d =  0.072
      an = (/14.6,14.6,14.7,14.6/) !a for logNe = 2,4,5,6 Case B
      a =  14.6
      b =  -0.732
      c =  0.081
      d =  0.066
      if ( log(ne) .le. 2) then
            a = an(1)
      elseif ( log(ne) .gt. 2 .and. log(ne) .le. 4) then
            a = an(1) + (an(2) - an(1)) / 2.*(log(ne) - 2.)
      elseif ( log(ne) .gt. 4 .and. log(ne) .le. 5) then
            a = an(2) + (an(3) - an(2))*(log(ne) - 2.)
      elseif ( log(ne) .gt. 6 .and. log(ne) .le. 6) then
            a = an(3) + (an(4) - an(3))*(log(ne) - 2.)
      else
              a = an(4)
      endif
!
      aeff = 1.e-14*a*te ** (b)
      aeff = aeff*(1. + c*(1. - te) + d*(1. - te) ** 2)
!
! New eff. recomb. data for 1000 < T < 5000 K, Ne = 100/cm3
! c.f. ~/source/OII_a_eff_at_low_Te/00readme/
      if ( te .le. 0.5) then
            a = 14.656
            b = -0.80449
            c = 0.00018
            d = 0.00517
            aeff = 1.e-14*a*te ** (b + c*log(te) + d* log(te) ** 2)
      endif
!
     do i = 397,403
        oiiRLs(i)%Em = aeff*1.98648E-08 / oiiRLs(i)%Wave*&
      & oiiRLs(i)%g2*oiiRLs(i)%Br_B
        oiiRLs(i)%Int = 100.*oiiRLs(i)%Em / Em_hb*abund
      enddo
!
! 3p-3s ^4S - ^4P transitions
!      an = (/0.90,0.90,0.90,1.00/) !a for logNe = 2,4,5,6 Case A
!      a =  0.90
!      b =  -0.485
!      c =  -0.047
!      d =  0.140
      an = (/4.80,4.90,4.90,4.90/) !a for logNe = 2,4,5,6 Case B
      a =  4.90
      b =  -0.730
      c =  -0.003
      d =  0.057
      if ( log(ne) .le. 2) then
            a = an(1)
      elseif ( log(ne) .gt. 2 .and. log(ne) .le. 4) then
            a = an(1) + (an(2) - an(1)) / 2.*(log(ne) - 2.)
      elseif ( log(ne) .gt. 4 .and. log(ne) .le. 5) then
            a = an(2) + (an(3) - an(2))*(log(ne) - 2.)
      elseif ( log(ne) .gt. 6 .and. log(ne) .le. 6) then
            a = an(3) + (an(4) - an(3))*(log(ne) - 2.)
      else
              a = an(4)
      endif
!
      aeff = 1.e-14*a*te ** (b)
      aeff = aeff*(1. + c*(1. - te) + d*(1. - te) ** 2)
!
! New eff. recomb. data for 1000 < T < 5000 K, Ne = 100/cm3
! c.f. ~/source/OII_a_eff_at_low_Te/00readme/
      if ( te .le. 0.5) then
            a = 4.8340
            b = -0.71947
            c = 0.02544
            d = 0.00936
            aeff = 1.e-14*a*te ** (b + c*log(te) + d* log(te) ** 2)
      endif
!
     do i = 404,406
        oiiRLs(i)%Em = aeff*1.98648E-08 / oiiRLs(i)%Wave*&
      & oiiRLs(i)%g2*oiiRLs(i)%Br_B
        oiiRLs(i)%Int = 100.*oiiRLs(i)%Em / Em_hb*abund
      enddo
!
! 3p-3s ^2D - ^2P transitions
      an = (/2.40,2.40,2.50,2.60/) !a for logNe = 2,4,5,6 Case A
      a =  2.40
      b =  -0.550
      c =  -0.051
      d =  0.178
!      an = (/14.5,14.6,14.5,14.3/) !a for logNe = 2,4,5,6 Case C
!      a =  14.6
!      b =  -0.736
!      c =  0.068
!      d =  0.066
      if ( log(ne) .le. 2) then
            a = an(1)
      elseif ( log(ne) .gt. 2 .and. log(ne) .le. 4) then
            a = an(1) + (an(2) - an(1)) / 2.*(log(ne) - 2.)
      elseif ( log(ne) .gt. 4 .and. log(ne) .le. 5) then
            a = an(2) + (an(3) - an(2))*(log(ne) - 2.)
      elseif ( log(ne) .gt. 6 .and. log(ne) .le. 6) then
            a = an(3) + (an(4) - an(3))*(log(ne) - 2.)
      else
              a = an(4)
      endif
!
      aeff = 1.e-14*a*te ** (b)
      aeff = aeff*(1. + c*(1. - te) + d*(1. - te) ** 2)
!
! New eff. recomb. data for 1000 < T < 5000 K, Ne = 100/cm3
! c.f. ~/source/OII_a_eff_at_low_Te/00readme/
      if ( te .le. 0.5) then
            a = 2.3616
            b = -0.46263
            c = 0.14697
            d = 0.03856
            aeff = 1.e-14*a*te ** (b + c*log(te) + d* log(te) ** 2)
      endif
!
     do i = 407,409
        oiiRLs(i)%Em = aeff*1.98648E-08 / oiiRLs(i)%Wave*&
      & oiiRLs(i)%g2*oiiRLs(i)%Br_A
        oiiRLs(i)%Int = 100.*oiiRLs(i)%Em / Em_hb*abund
      enddo
!
! 3p-3s ^2P - ^2P transitions
      an = (/1.10,1.20,1.20,1.20/) !a for logNe = 2,4,5,6 Case A
      a =  1.20
      b =  -0.523
      c =  -0.044
      d =  0.173
!      an = (/1.30,1.40,1.40,1.40/) !a for logNe = 2,4,5,6 Case C
!      a =  1.40
!      b =  -0.565
!      c =  -0.042
!      d =  0.158
      if ( log(ne) .le. 2) then
            a = an(1)
      elseif ( log(ne) .gt. 2 .and. log(ne) .le. 4) then
            a = an(1) + (an(2) - an(1)) / 2.*(log(ne) - 2.)
      elseif ( log(ne) .gt. 4 .and. log(ne) .le. 5) then
            a = an(2) + (an(3) - an(2))*(log(ne) - 2.)
      elseif ( log(ne) .gt. 6 .and. log(ne) .le. 6) then
            a = an(3) + (an(4) - an(3))*(log(ne) - 2.)
      else
              a = an(4)
      endif
!
      aeff = 1.e-14*a*te ** (b)
      aeff = aeff*(1. + c*(1. - te) + d*(1. - te) ** 2)
!
! New eff. recomb. data for 1000 < T < 5000 K, Ne = 100/cm3
! c.f. ~/source/OII_a_eff_at_low_Te/00readme/
      if ( te .le. 0.5) then
            a = 1.1198
            b = -0.44147
            c = 0.13837
            d = 0.03191
            aeff = 1.e-14*a*te ** (b + c*log(te) + d* log(te) ** 2)
      endif
!
     do i = 410,413
        oiiRLs(i)%Em = aeff*1.98648E-08 / oiiRLs(i)%Wave*&
      & oiiRLs(i)%g2*oiiRLs(i)%Br_A
        oiiRLs(i)%Int = 100.*oiiRLs(i)%Em / Em_hb*abund
      enddo
!
! 3p-3s ^2S - ^2P transitions
      an = (/0.40,0.40,0.40,0.40/) !a for logNe = 2,4,5,6 Case A
      a =  0.40
      b =  -0.461
      c =  -0.083
      d =  0.287
!      an = (/0.50,0.50,0.50,0.60/) !a for logNe = 2,4,5,6 Case C
!      a =  0.50
!      b =  -0.547
!      c =  -0.074
!      d =  0.244
      if ( log(ne) .le. 2) then
            a = an(1)
      elseif ( log(ne) .gt. 2 .and. log(ne) .le. 4) then
            a = an(1) + (an(2) - an(1)) / 2.*(log(ne) - 2.)
      elseif ( log(ne) .gt. 4 .and. log(ne) .le. 5) then
            a = an(2) + (an(3) - an(2))*(log(ne) - 2.)
      elseif ( log(ne) .gt. 6 .and. log(ne) .le. 6) then
            a = an(3) + (an(4) - an(3))*(log(ne) - 2.)
      else
              a = an(4)
      endif
!
      aeff = 1.e-14*a*te ** (b)
      aeff = aeff*(1. + c*(1. - te) + d*(1. - te) ** 2)
!
! New eff. recomb. data for 1000 < T < 5000 K, Ne = 100/cm3
! c.f. ~/source/OII_a_eff_at_low_Te/00readme/
      if ( te .le. 0.5) then
            a = 0.3922
            b = -0.35043
            c = 0.26366
            d = 0.06666
            aeff = 1.e-14*a*te ** (b + c*log(te) + d* log(te) ** 2)
      endif
!
      do i = 414,415
        oiiRLs(i)%Em = aeff*1.98648E-08 / oiiRLs(i)%Wave*&
      & oiiRLs(i)%g2*oiiRLs(i)%Br_A
        oiiRLs(i)%Int = 100.*oiiRLs(i)%Em / Em_hb*abund
      enddo

      te = te * 10000

      end subroutine oii_rec_lines

      subroutine nii_rec_lines(te, ne, abund, niiRLs)

      IMPLICIT NONE
      DOUBLE PRECISION :: aeff, aeff_hb, Em_Hb, &
     & Te, Ne, abund, &
     & Br_term, z
      DOUBLE PRECISION :: a, b, c, d

      INTEGER :: ii, i


      TYPE(niiRL), DIMENSION(99) :: niiRLs

      call get_aeff_hb(te,ne, aeff_hb, em_hb)

!! read in NII data
!
!            301 FORMAT (I5, 1X, F9.4, 1X, A1, A1, A1, A1, A1, 1X, F7.4, &
!     & 1X, A3, 1X, F7.4, 1X, A3, 1X, A7, 3X, F11.4, A1, A1, 1X, I2, &
!     &1X, A1, 1X, A9, 1X, F13.4, 1X, A1, A1, 1X, I2, 1X, A1, 1X, A9, 1X,&
!     & F7.4, 1X, F7.4, 1X, F7.4)!, 1X, E10.4, 1X, E10.4, 1X)
!            OPEN(201, file="Atomic-data/Rnii.dat", status='old')
!            DO i = 1,99
!            READ(201,301) niiRLs(i)%ION, niiRLs(i)%Wave, niiRLs(i)%Hyb, &
!     &niiRLs(i)%Rem1, niiRLs(i)%Rem2, niiRLs(i)%Rem3, niiRLs(i)%Rem4,   &
!     &niiRLs(i)%gf1, niiRLs(i)%q_gf1, niiRLs(i)%gf2, niiRLs(i)%q_gf2,   &
!     &niiRLs(i)%Mult, niiRLs(i)%E1, niiRLs(i)%n_E1, niiRLs(i)%n_E1GA,   &
!     &niiRLs(i)%g1, niiRLs(i)%n_g1, niiRLs(i)%Term1, niiRLs(i)%E2,      &
!     &niiRLs(i)%n_E2, niiRLs(i)%n_E2GA, niiRLs(i)%g2, niiRLs(i)%n_g2,   &
!     &niiRLs(i)%Term2, niiRLs(i)%Br_LS
!            END DO
!      CLOSE(201)

      te = te/10000

!       2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p E1 3P* - 3D  M03  transitions
!      i = 0002     !case A
      i = 0003     !case B
      a = -12.7289
      b = -0.689816
      c = 0.022005
!
      aeff = 10. ** (a + b * log10(te) + c * log10(te) ** 2)
      do ii = 1,6
        niiRLs(ii)%Em = aeff * 1.98648E-08 / niiRLs(ii)%Wave * niiRLs(ii)%Br_LS
        niiRLs(ii)%Int = 100 * niiRLs(ii)%Em / Em_Hb * abund
      enddo
!
!       2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p 3P* - 3S     M04 transitions
!      i = 0004     !case A
      i = 0005     !case B
      a = -13.8161
      b = -0.778606
      c = -0.028944
!
      aeff = 10. ** (a + b * log10(te) + c * log10(te) ** 2)
      do ii = 7,9
        niiRLs(ii)%Em = aeff * 1.98648E-08 / niiRLs(ii)%Wave * niiRLs(ii)%Br_LS
        niiRLs(ii)%Int = 100 * niiRLs(ii)%Em / Em_Hb * abund
      enddo
!
!       2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p 3P* - 3P     M05 transitions
!      i = 0006     !case A
      i = 0007     !case B
      a = -13.0765
      b = -0.734594
      c = -0.0251909
!
      aeff = 10. ** (a + b * log10(te) + c * log10(te) ** 2)
      do ii = 10,15
        niiRLs(ii)%Em = aeff * 1.98648E-08 / niiRLs(ii)%Wave * niiRLs(ii)%Br_LS
        niiRLs(ii)%Int = 100 * niiRLs(ii)%Em / Em_Hb * abund
      enddo
!
!       2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p 1P* - 1P     M08 transitions
      i = 0008     !case A
!      i = 0009     !case B
      a = -14.1211
      b = -0.608107
      c = 0.0362301
!
      aeff = 10. ** (a + b * log10(te) + c * log10(te) ** 2)
      do ii = 16,16
        niiRLs(ii)%Em = aeff * 1.98648E-08 / niiRLs(ii)%Wave * niiRLs(ii)%Br_LS
        niiRLs(ii)%Int = 100 * niiRLs(ii)%Em / Em_Hb * abund
      enddo
!
!       2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p 1P* - 1D     M12 transitions
      i = 0010     !case A
!      i = 0011     !case B
      a = -13.7473
      b = -0.509595
      c = 0.0255685
!
      aeff = 10. ** (a + b * log10(te) + c * log10(te) ** 2)
      do ii = 17,17
        niiRLs(ii)%Em = aeff * 1.98648E-08 / niiRLs(ii)%Wave * niiRLs(ii)%Br_LS
        niiRLs(ii)%Int = 100 * niiRLs(ii)%Em / Em_Hb * abund
      enddo
!
!       2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p 1P* - 1S     M13 transitions
      i = 0012     !case A
!      i = 0013     !case B
      a = -14.3753
      b = -0.515547
      c = 0.0100966
!
      aeff = 10. ** (a + b * log10(te) + c * log10(te) ** 2)
      do ii = 18,18
        niiRLs(ii)%Em = aeff * 1.98648E-08 / niiRLs(ii)%Wave * niiRLs(ii)%Br_LS
        niiRLs(ii)%Int = 100 * niiRLs(ii)%Em / Em_Hb * abund
      enddo
!
!       2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 1P - 1D*     M15 transitions
      i = 0014     !case A
!      i = 0015     !case B
      a = -14.3932
      b = -0.887946
      c = -0.0525855
!
      aeff = 10. ** (a + b * log10(te) + c * log10(te) ** 2)
      do ii = 19,19
        niiRLs(ii)%Em = aeff * 1.98648E-08 / niiRLs(ii)%Wave * niiRLs(ii)%Br_LS
        niiRLs(ii)%Int = 100 * niiRLs(ii)%Em / Em_Hb * abund
      enddo
!
!       2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 1P - 1P*     M17 transitions
      i = 0016     !case A
!      i = 0017     !case B
      a = -15.0052
      b = -0.89811
      c = -0.0581789
!
      aeff = 10. ** (a + b * log10(te) + c * log10(te) ** 2)
      do ii = 20,20
        niiRLs(ii)%Em = aeff * 1.98648E-08 / niiRLs(ii)%Wave * niiRLs(ii)%Br_LS
        niiRLs(ii)%Int = 100 * niiRLs(ii)%Em / Em_Hb * abund
      enddo
!
!       2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3D - 3F*     M19 transitions
!      i = 0018     !case A
      i = 0019     !case B
      a = -12.6183
      b = -0.840727
      c = -0.0229685
!
      aeff = 10. ** (a + b * log10(te) + c * log10(te) ** 2)
      do ii = 21,26
        niiRLs(ii)%Em = aeff * 1.98648E-08 / niiRLs(ii)%Wave * niiRLs(ii)%Br_LS
        niiRLs(ii)%Int = 100 * niiRLs(ii)%Em / Em_Hb * abund
      enddo
!
!       2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3D - 3D*     M20 transitions
!      i = 0020     !case A
      i = 0021     !case B
      a = -13.3184
      b = -0.884034
      c = -0.0512093
!
      aeff = 10. ** (a + b * log10(te) + c * log10(te) ** 2)
      do ii = 27,33
        niiRLs(ii)%Em = aeff * 1.98648E-08 / niiRLs(ii)%Wave * niiRLs(ii)%Br_LS
        niiRLs(ii)%Int = 100 * niiRLs(ii)%Em / Em_Hb * abund
      enddo
!
!       2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3D - 3P*     M21 transitions
!      i = 0022     !case A
      i = 0023     !case B
      a = -14.5113
      b = -0.87792
      c = -0.0552785
!
      aeff = 10. ** (a + b * log10(te) + c * log10(te) ** 2)
      do ii = 34,39
        niiRLs(ii)%Em = aeff * 1.98648E-08 / niiRLs(ii)%Wave * niiRLs(ii)%Br_LS
        niiRLs(ii)%Int = 100 * niiRLs(ii)%Em / Em_Hb * abund
      enddo
!
!       2s2.2p.(2P*).3p - 2s2.2p.(2P*).4s 3D - 3P*     M22 transitions
!      i = 0024     !case A
      i = 0025     !case B
      a = -14.1305
      b = -0.487037
      c = 0.0354135
!
      aeff = 10. ** (a + b * log10(te) + c * log10(te) ** 2)
      do ii = 40,45
        niiRLs(ii)%Em = aeff * 1.98648E-08 / niiRLs(ii)%Wave * niiRLs(ii)%Br_LS
        niiRLs(ii)%Int = 100 * niiRLs(ii)%Em / Em_Hb * abund
      enddo
!
!       2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3S - 3P*     M24 transitions
!      i = 0026     !case A
      i = 0027     !case B
      a = -13.3527
      b = -0.878224
      c = -0.0557112
!
      aeff = 10. ** (a + b * log10(te) + c * log10(te) ** 2)
      do ii = 46,48
        niiRLs(ii)%Em = aeff * 1.98648E-08 / niiRLs(ii)%Wave * niiRLs(ii)%Br_LS
        niiRLs(ii)%Int = 100 * niiRLs(ii)%Em / Em_Hb * abund
      enddo
!
!       2s2.2p.(2P*).3p - 2s2.2p.(2P*).4s 3S - 3P*     M26 transitions
!      i = 0028     !case A
      i = 0029     !case B
      a = -14.9628
      b = -0.486746
      c = 0.0358261
!
      aeff = 10. ** (a + b * log10(te) + c * log10(te) ** 2)
      do ii = 49,51
        niiRLs(ii)%Em = aeff * 1.98648E-08 / niiRLs(ii)%Wave * niiRLs(ii)%Br_LS
        niiRLs(ii)%Int = 100 * niiRLs(ii)%Em / Em_Hb * abund
      enddo
!
!       2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3P - 3D*     M28 transitions
!      i = 0030     !case A
      i = 0031     !case B
      a = -13.0871
      b = -0.883624
      c = -0.0506882
!
      aeff = 10. ** (a + b * log10(te) + c * log10(te) ** 2)
      do ii = 52,57
        niiRLs(ii)%Em = aeff * 1.98648E-08 / niiRLs(ii)%Wave * niiRLs(ii)%Br_LS
        niiRLs(ii)%Int = 100 * niiRLs(ii)%Em / Em_Hb * abund
      enddo
!
!       2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3P - 3P*     M29 transitions
!      i = 0032     !case A
      i = 0033     !case B
      a = -13.5581
      b = -0.878488
      c = -0.0557583
!
      aeff = 10. ** (a + b * log10(te) + c * log10(te) ** 2)
      do ii = 58,63
        niiRLs(ii)%Em = aeff * 1.98648E-08 / niiRLs(ii)%Wave * niiRLs(ii)%Br_LS
        niiRLs(ii)%Int = 100 * niiRLs(ii)%Em / Em_Hb * abund
      enddo
!
!       2s2.2p.(2P*).3p - 2s2.2p.(2P*).4s 3P - 3P*     M30 transitions
!      i = 0034     !case A
      i = 0035     !case B
      a = -14.3521
      b = -0.487527
      c = 0.0355516
!
      aeff = 10. ** (a + b * log10(te) + c * log10(te) ** 2)
      do ii = 64,69
        niiRLs(ii)%Em = aeff * 1.98648E-08 / niiRLs(ii)%Wave * niiRLs(ii)%Br_LS
        niiRLs(ii)%Int = 100 * niiRLs(ii)%Em / Em_Hb * abund
      enddo
!
!       2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 1D - 1F*     M31 transitions
      i = 0036     !case A
!      i = 0037     !case B
      a = -15.0026
      b = -0.923093
      c = -0.0588371
!
      aeff = 10. ** (a + b * log10(te) + c * log10(te) ** 2)
      do ii = 70,70
        niiRLs(ii)%Em = aeff * 1.98648E-08 / niiRLs(ii)%Wave * niiRLs(ii)%Br_LS
        niiRLs(ii)%Int = 100 * niiRLs(ii)%Em / Em_Hb * abund
      enddo
!
!       2s2.2p.(2P*).3d - 2s2.2p.(2P*).4p 3F* - 3D     M36 transitions
!      i = 0038     !case A
      i = 0039     !case B
      a = -13.8636
      b = -0.569144
      c = 0.0068655
!
      aeff = 10. ** (a + b * log10(te) + c * log10(te) ** 2)
      do ii = 71,76
        niiRLs(ii)%Em = aeff * 1.98648E-08 / niiRLs(ii)%Wave * niiRLs(ii)%Br_LS
        niiRLs(ii)%Int = 100 * niiRLs(ii)%Em / Em_Hb * abund
      enddo
!
!       2s2.2p.(2P*).3d - 2s2.2p.(2P*<3/2>).4f 3F* - 3G M39 transitions
!      i = 0040     !case A
      i = 0041     !case B
      a = -13.035
      b = -1.12035
      c = -0.10642
!
      aeff = 10. ** (a + b * log10(te) + c * log10(te) ** 2)
      do ii = 77,82
        niiRLs(ii)%Em = aeff * 1.98648E-08 / niiRLs(ii)%Wave * niiRLs(ii)%Br_LS
        niiRLs(ii)%Int = 100 * niiRLs(ii)%Em / Em_Hb * abund
      enddo
!
!       2s2.2p.(2P*).3d - 2s2.2p.(2P*<3/2>).4f 1F* - 1G M58 transitions
      i = 0042     !case A
!      i = 0043     !case B
      a = -13.5484
      b = -1.11909
      c = -0.105123
!
      aeff = 10. ** (a + b * log10(te) + c * log10(te) ** 2)
      do ii = 83,83
        niiRLs(ii)%Em = aeff * 1.98648E-08 / niiRLs(ii)%Wave * niiRLs(ii)%Br_LS
        niiRLs(ii)%Int = 100 * niiRLs(ii)%Em / Em_Hb * abund
      enddo
!
!       3d 3D* - 4f 3F 4242 M48 transitions
!      i = 0044     !case A
      i = 0045     !case B
      a = -13.2548
      b = -1.12902
      c = -0.110368
!
      aeff = 10. ** (a + b * log10(te) + c * log10(te) ** 2)
      do ii = 84,89
        niiRLs(ii)%Em = aeff * 1.98648E-08 / niiRLs(ii)%Wave * niiRLs(ii)%Br_LS
        niiRLs(ii)%Int = 100 * niiRLs(ii)%Em / Em_Hb * abund
      enddo
!
!       3d 3P* - 4f 3D 4435 M55 transitions
!      i = 0046     !case A
      i = 0047     !case B
      a = -13.5656
      b = -1.11989
      c = -0.105818
!
      aeff = 10. ** (a + b * log10(te) + c * log10(te) ** 2)
      do ii = 90,95
        niiRLs(ii)%Em = aeff * 1.98648E-08 / niiRLs(ii)%Wave * niiRLs(ii)%Br_LS
        niiRLs(ii)%Int = 100 * niiRLs(ii)%Em / Em_Hb * abund
      enddo
!
!       3d 1D* - 4f 1F 4176 M43 (RMT M42) transitions
      i = 0048     !case A
!      i = 0049     !case B
      a = -13.7426
      b = -1.13351
      c = -0.111146
!
      aeff = 10. ** (a + b * log10(te) + c * log10(te) ** 2)
      do ii = 96,96
        niiRLs(ii)%Em = aeff * 1.98648E-08 / niiRLs(ii)%Wave * niiRLs(ii)%Br_LS
        niiRLs(ii)%Int = 100 * niiRLs(ii)%Em / Em_Hb * abund
      enddo
!
!       3d 1P* - 4f 1D 4677 M61 (RMT M62) transitions
      i = 0049     !case A
!      i = 0051     !case B
      a = -13.7373
      b = -1.12695
      c = -0.108158
!
      aeff = 10. ** (a + b * log10(te) + c * log10(te) ** 2)
      do ii = 97,97
        niiRLs(ii)%Em = aeff * 1.98648E-08 / niiRLs(ii)%Wave * niiRLs(ii)%Br_LS
        niiRLs(ii)%Int = 100 * niiRLs(ii)%Em / Em_Hb * abund
      enddo
!
!       3d 3F* - 4f 1G 4026 M39b transitions
!      case A (PPB):
      a = 0.108
      b = -0.754
      c = 2.587
      d = 0.719
      z = 2.
      Br_term = 0.350
!
      aeff = 1.e-13 * z * a  * (te/z**2) ** (b)
      aeff = aeff / (1. + c * (te/z**2) ** (d)) * Br_term
      do ii = 98,98
        niiRLs(ii)%Em = aeff * 1.98648E-08 / niiRLs(ii)%Wave * niiRLs(ii)%Br_LS
        niiRLs(ii)%Int = 100 * niiRLs(ii)%Em / Em_Hb * abund
      enddo
!
!       3d 1F* - 4f 3G 4552 M58a transitions
!      case A (PPB):
      a = 0.326
      b = -0.754
      c = 2.587
      d = 0.719
      z = 2.
      Br_term = 0.074
!
      aeff = 1.e-13 * z * a  * (te/z**2) ** (b)
      aeff = aeff / (1. + c * (te/z**2) ** (d)) * Br_term
      do ii = 99,99
        niiRLs(ii)%Em = aeff * 1.98648E-08 / niiRLs(ii)%Wave * niiRLs(ii)%Br_LS
        niiRLs(ii)%Int = 100 * niiRLs(ii)%Em / Em_Hb * abund
      enddo

      te = te * 10000

      end subroutine nii_rec_lines

      subroutine cii_rec_lines(te, ne, abund, ciiRLs)

      IMPLICIT NONE
      DOUBLE PRECISION :: aeff_Hb, Em_Hb, &
     & Te, Ne, abund

      INTEGER :: i


      TYPE(ciiRL), DIMENSION(57) :: ciiRLs

      call get_aeff_hb(te,ne, aeff_hb, em_hb)

!! read in CII data
!
!       301 FORMAT (F7.2, 1X, F6.4, 1X, F7.4, 1X, F7.4, 1X, F7.4, 1X, F7.4)
!       OPEN(201, file="Atomic-data/Rcii.dat", status='old')
!       DO i = 1,57
!         READ(201,301) ciiRLs(i)%Wave, ciiRLs(i)%a, ciiRLs(i)%b, &
!         & ciiRLs(i)%c, ciiRLs(i)%d, ciiRLs(i)%f
!       END DO
!       CLOSE(201)

      te = te/10000

      do i = 1,57
        ciiRLs(i)%aeff = 1e-14 * (ciiRLs(i)%a*(te**ciiRLs(i)%f)) * (1 &
        &+ (ciiRLs(i)%b*(1-te)) &
        &+ (ciiRLs(i)%c * ((1-te)**2) ) &
        &+ (ciiRLs(i)%d * ((1-te)**3) ) &
        &)
        ciiRLs(i)%Int = 100 * (ciiRLs(i)%aeff/aeff_hb) * (4861.33/ciiRLs(i)%Wave) * abund
      enddo

      te = te * 10000

      end subroutine cii_rec_lines

      subroutine neii_rec_lines(te, ne, abund, neiiRLs)

      IMPLICIT NONE
      DOUBLE PRECISION :: aeff_Hb, Em_Hb,                         &
     &  Te, Ne, abund

      INTEGER :: i


      TYPE(neiiRL), DIMENSION(38) :: neiiRLs

      call get_aeff_hb(te,ne, aeff_hb, em_hb)

!! read in NeII data
!
!       301 FORMAT (F7.2, 1X, F6.3, 1X, F6.3, 1X, F6.3, 1X, F6.3, 1X, F7.4, 1X, F6.3)
!       OPEN(201, file="Atomic-data/Rneii.dat", status='old')
!       DO i = 1,38
!         READ(201,301) neiiRLs(i)%Wave, neiiRLs(i)%a, neiiRLs(i)%b, &
!         & neiiRLs(i)%c, neiiRLs(i)%d, neiiRLs(i)%f, neiiRLs(i)%Br
!       END DO
!       CLOSE(201)

      te = te/10000

      do i = 1,38
        neiiRLs(i)%aeff = neiiRLs(i)%Br * 1e-14 * &
        &(neiiRLs(i)%a*(te**neiiRLs(i)%f)) * (1 &
        &+ (neiiRLs(i)%b*(1-te)) &
        &+ (neiiRLs(i)%c * ((1-te)**2) ) &
        &+ (neiiRLs(i)%d * ((1-te)**3) ) &
        &)
        neiiRLs(i)%Int = 100 * (neiiRLs(i)%aeff/aeff_hb) * (4861.33/neiiRLs(i)%Wave) * abund
      enddo

      te = te * 10000

      end subroutine neii_rec_lines

      subroutine xiii_rec_lines(te, ne, abund, xiiiRLs)

      IMPLICIT NONE
      DOUBLE PRECISION :: aeff_Hb, Em_Hb, &
     & Te, Ne, abund

      INTEGER :: i

      TYPE(xiiiRL), DIMENSION(6) :: xiiiRLs

      call get_aeff_hb(te,ne, aeff_hb, em_hb)

!! read in XIII data
!
!       301 FORMAT (A3,1X,F7.2, 1X, F5.3, 1X, F6.3, 1X, F5.3, 1X, F5.3, 1X, F5.4)
!       OPEN(201, file="Atomic-data/Rxiii.dat", status='old')
!       DO i = 1,6
!         READ(201,301) xiiiRLs(i)%ion, xiiiRLs(i)%Wave, xiiiRLs(i)%a, &
!         & xiiiRLs(i)%b, xiiiRLs(i)%c, xiiiRLs(i)%d, xiiiRLs(i)%Br
!       END DO
!       CLOSE(201)

      te = te/90000 !ionic charge=3 so divide by 9

      do i = 1,4
        xiiiRLs(i)%aeff = xiiiRLs(i)%Br * 1e-13 * 3 * &
      & (xiiiRLs(i)%a*(te**xiiiRLs(i)%b)) / &
      & (1 + (xiiiRLs(i)%c * (te**xiiiRLs(i)%d)))
        xiiiRLs(i)%Int = 100 * (xiiiRLs(i)%aeff/aeff_hb) * (4861.33/xiiiRLs(i)%Wave) * abund
      enddo

      do i = 5,6
        xiiiRLs(i)%aeff = xiiiRLs(i)%Br * 1e-13 * 3 * &
      & (xiiiRLs(i)%a*(te**xiiiRLs(i)%b)) / &
      & (1 + (xiiiRLs(i)%c * (te**xiiiRLs(i)%d)))
        xiiiRLs(i)%Int = 100 * (xiiiRLs(i)%aeff/aeff_hb) * (4861.33/xiiiRLs(i)%Wave) * abund
      enddo

      te = te * 90000

      end subroutine xiii_rec_lines

      subroutine get_aeff_hb(te, ne, aeff_hb, em_hb)
      IMPLICIT NONE
      double precision :: Te, Ne, AE2, AE3, AE4, AE5, AE6, AE7, AE8, aeff_hb, Em_Hb, logem

      AE2 = -9.06524E+00 -2.69954E+00 * log10(te) + 8.80123E-01 * &
      &log10(te) ** 2 -1.57946E-01 * log10(te) ** 3 + &
      &9.25920E-03 * log10(te) ** 4
      AE3 = -8.13757E+00 -3.57392E+00 * log10(te) + 1.19331E+00 * &
      &log10(te) ** 2 -2.08362E-01 * log10(te) ** 3 + &
      &1.23303E-02 * log10(te) ** 4
      AE4 = -6.87230E+00 -4.72312E+00 * log10(te) + 1.58890E+00 * &
      &log10(te) ** 2 -2.69447E-01 * log10(te) ** 3 + &
      &1.58955E-02 * log10(te) ** 4
      AE5 = -5.15059E+00 -6.24549E+00 * log10(te) + 2.09801E+00 * &
      &log10(te) ** 2 -3.45649E-01 * log10(te) ** 3 + &
      &2.01962E-02 * log10(te) ** 4
      AE6 = -2.35923E+00 -8.75565E+00 * log10(te) + 2.95600E+00 * &
      &log10(te) ** 2 -4.77584E-01 * log10(te) ** 3 + &
      &2.78852E-02 * log10(te) ** 4
      AE7 =  1.55373E+00 -1.21894E+01 * log10(te) + 4.10096E+00 * &
      &log10(te) ** 2 -6.49318E-01 * log10(te) ** 3 + &
      &3.76487E-02 * log10(te) ** 4
      AE8 =  6.59883E+00 -1.64030E+01 * log10(te) + 5.43844E+00 * &
      &log10(te) ** 2 -8.40253E-01 * log10(te) ** 3 + &
      &4.79786E-02 * log10(te) ** 4

      if (log10(ne) .lt. 2) then
            aeff_hb = ae2
      elseif (log10(ne) .GE. 2 .AND. log10(ne) .LT. 3) then
            aeff_hb = AE2 + (AE3 - AE2) * (log10(ne) - 2)
      elseif (log10(ne) .GE. 3 .AND. log10(ne) .LT. 4) then
            aeff_hb = AE3 + (AE4 - AE3) * (log10(ne) - 3)
      elseif (log10(ne) .GE. 4 .AND. log10(ne) .LT. 5) then
            aeff_hb = AE4 + (AE5 - AE4) * (log10(ne) - 4)
      elseif (log10(ne) .GE. 5 .AND. log10(ne) .LT. 6) then
            aeff_hb = AE5 + (AE6 - AE5) * (log10(ne) - 5)
      elseif (log10(ne) .GE. 6 .AND. log10(ne) .LT. 7) then
            aeff_hb = AE6 + (AE7 - AE6) * (log10(ne) - 6)
      elseif (log10(ne) .GE. 7 .AND. log10(ne) .LT. 8) then
            aeff_hb = AE7 + (AE8 - AE7) * (log10(ne) - 7)
      else
            aeff_hb = AE8
      endif

      LogEm = aeff_hb - 11.38871 ! = log10(hc/lambda in cgs)
      aeff_hb = 10**aeff_hb
      em_hb = 10**logem

      end subroutine get_aeff_hb

      end module mod_recombination_lines
