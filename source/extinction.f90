!extinction.f90, calculation of extinction coefficients for various interstellar reddening laws
!(C) Roger Wesson, Dave Stock, Peter Scicluna
module mod_extinction
use mod_types
use mod_globals
use mod_hydrogen

implicit none

contains

subroutine calc_extinction_coeffs(linelist,H_Balmer, c1, c2, c3, meanextinction, temp, dens, weightha, weighthg, weighthd,hidata)
        IMPLICIT NONE
        type(line), dimension(:), intent(in) :: linelist
        integer, dimension(3:40) :: H_Balmer
        real(kind=dp) :: c1, c2, c3, weightha, weighthg, weighthd, meanextinction
        real(kind=dp) :: temp, dens
        real(kind=dp), dimension(:,:,:,:), allocatable :: hidata
        real(kind=dp), dimension(3:40) :: balmerratios,extinctioncoefficients,r1,r2,r3,r4,r5,r6
        integer :: temperaturestart,densitystart,temperatureend,densityend,i
        real(kind=dp) :: temperatureinterpolation,densityinterpolation

!debugging
#ifdef CO
        print *,"subroutine: calc_extinction_coeffs. Te=",temp,", ne=",dens,", weights=",weightha, weighthg, weighthd
#endif

!interpolate HI emissivities to temp and dens
!first temperature

        do temperatureend=1,ntemps
          if (temp .lt. temperatures(temperatureend)) then
            exit
          endif
        enddo

        if (temperatureend.eq.1) then
          temperatureend=1
          temperatureinterpolation=0.
        elseif (temperatureend.ge.ntemps) then ! because in do i=1,10, final value of i is 11
          temperaturestart=ntemps
          temperatureend=ntemps
          temperatureinterpolation=0.
        else
          temperaturestart=temperatureend-1
          temperatureinterpolation=(temp-temperatures(temperaturestart))/(temperatures(temperatureend)-temperatures(temperaturestart))
        endif

!then density which runs from 2 to 8 in the data file
!so that log10(density)=2 corresponds to array index 1
!note: the following line if uncommented causes a segmentation fault. why?
!dens=1200.

        densitystart=floor(log10(dens))-1
        if (densitystart.lt.1) then
          densitystart=1
          densityend=1
          densityinterpolation=0.
        elseif (densitystart.gt.7) then
          densitystart=7
          densityend=7
          densityinterpolation=0.
        else
          densityend=densitystart+1
          densityinterpolation=(log10(dens)-densitystart)/(densityend-densitystart)
        endif

!now interpolate

        r1=hidata(temperaturestart,densitystart,3:40,2) / hidata(temperaturestart,densitystart,4,2)
        r2=hidata(temperaturestart,densityend,3:40,2) / hidata(temperaturestart,densityend,4,2)
        r3=hidata(temperatureend,densitystart,3:40,2) / hidata(temperatureend,densitystart,4,2)
        r4=hidata(temperatureend,densityend,3:40,2) / hidata(temperatureend,densityend,4,2)

        r5=densityinterpolation*(r2-r1)+r1
        r6=densityinterpolation*(r4-r3)+r3

        balmerratios=temperatureinterpolation*(r6-r5)+r5

! calculate extinction coefficients

        extinctioncoefficients=0
        where (linelist(H_Balmer(:))%intensity .gt. 0)
          extinctioncoefficients=log10( ( linelist(H_Balmer(:))%intensity / linelist(H_Balmer(4))%intensity )/ balmerratios(:) )/(-linelist(H_Balmer(:))%flambda)
        endwhere

! set negative values to zero

        if (any(extinctioncoefficients.lt.0)) then
          where(extinctioncoefficients.lt.0)
            extinctioncoefficients=0
          endwhere
!          print *,"some extinction coefficients were negative, set to zero"
        endif

! calculate weighted mean, Ha to Hd for now

!        meanextinction=sum(extinctioncoefficients*linelist(H_Balmer(:))%weight)/sum(linelist(H_Balmer(:))%weight)

        if (sum(linelist(H_Balmer(3:6))%weight).gt.0) then
          meanextinction=sum(extinctioncoefficients(3:6)*linelist(H_Balmer(3:6))%weight)/(sum(linelist(H_Balmer(3:6))%weight)-linelist(H_Balmer(4))%weight)
        else
          meanextinction=0
        endif

        c1=extinctioncoefficients(3)
        c2=extinctioncoefficients(5)
        c3=extinctioncoefficients(6)

end subroutine calc_extinction_coeffs

subroutine get_flambda(linelist, switch_ext, R)
!calculat flambda for the chosen law, or if it was not recognised, attempt to read in from file
!todo: change length of switch_ext variable to allow for a filename to be given
        implicit none
        type(line), dimension(:) :: linelist
        character(len=1) :: switch_ext
        real(kind=dp) :: R
        real(kind=dp) :: x_0, c1, c2, c3, c4, gamma, F ! Fitzpatrick constants

!debugging
#ifdef CO
        print *,"subroutine: get_flambda, switch=",switch_ext
#endif

        c1 = -0.38
        c2 =  0.74
        c3 =  3.96
        c4 = 0.26
        x_0 = 4.596
        gamma = 1.051
        F = 0

        where (linelist%wavelength.lt.500) ! allows line list to contain IR lines with wavelengths in microns.
                linelist%freq = DBLE(1) / DBLE(linelist%wavelength)
        elsewhere
                linelist%freq = DBLE(10000) / DBLE(linelist%wavelength)
        endwhere

        if (switch_ext == "S") then !Howarth 1983 galactic law

                where( (linelist%freq .gt. 7.14) .AND. (linelist%freq .le. 10.0)) !far UV
                        linelist%flambda = ((16.07 - (3.20 * linelist%freq) + (0.2975*linelist%freq**2))) !far UV
                elsewhere((linelist%freq .gt. 3.65) .AND. (linelist%freq .le. 7.14)) ! mid UV
                        linelist%flambda = ((2.19 + (0.848*linelist%freq) + (1.01/(((linelist%freq-4.60)**2) + 0.280)))) !mid UV
                elsewhere((linelist%freq .gt. 2.75) .AND. (linelist%freq .le. 3.65)) ! near UV
                        linelist%flambda = ((1.46 + (1.048*linelist%freq) + (1.01/(((linelist%freq-4.60)**2) + 0.280)))) !near UV
                elsewhere((linelist%freq .gt. 1.83) .AND. (linelist%freq .le. 2.75)) ! optical
                        linelist%flambda = ((3.1 + (2.56*(linelist%freq-1.83)) - (0.993*(linelist%freq-1.83)**2) )) ! optical (3)
                elsewhere(linelist%freq .le. 1.83) !IR
                        linelist%flambda = (( (1.86*linelist%freq**2) - (0.48*linelist%freq**3) - (0.1*linelist%freq))) ! IR (4)
                elsewhere
!                       print *, "He's not the messiah!"
                        linelist%flambda = 0.d0
                endwhere

                linelist%flambda = (linelist%flambda / 3.63) - 1 ! R+0.53

        elseif (switch_ext == "H") then !Howarth 1983 LMC law

                where (linelist%freq .gt. 2.75) !UV
                    linelist%flambda = ( (3.1 - 0.236 + (0.462*linelist%freq) + (0.105*(linelist%freq**2)) + (0.454/((linelist%freq-4.557)**2 + 0.293))) )
                elsewhere ((linelist%freq .gt. 1.83) .AND. (linelist%freq .lt. 2.75)) !Optical
                    linelist%flambda = ((3.1 + (2.04*(linelist%freq - 1.83)) + 0.094*(linelist%freq - 1.83)**2))
                elsewhere(linelist%freq .lt. 1.83) !IR
                    linelist%flambda = (((1.86*(linelist%freq**2)) - (0.48*(linelist%freq**3)) - (0.1*linelist%freq)))
                elsewhere
!                   print *,"Your mother was a hamster and your father smelt of elderberries."
                    linelist%flambda = 0.d0
                endwhere

                linelist%flambda = (linelist%flambda / 3.57) - 1 ! R+0.47

        elseif (switch_ext == "C") then !CCM 1989 galactic law.
                where(linelist%freq .gt. 8) !far UV
                        linelist%flambda = R* (0.137*((linelist%freq - 8)**2) - 1.073 - 0.628*(linelist%freq - 8) - 0.070*((linelist%freq - 8)**3)) &
        &                                   + 13.670 + 4.257*(linelist%freq - 8) - 0.420*((linelist%freq - 8)**2) + 0.374*((linelist%freq - 8)**3)
                elsewhere((linelist%freq .gt. 5.9) .AND. (linelist%freq .lt. 8)) !UV
                        linelist%flambda = R* &
                                           &(1.752 - 0.316*linelist%freq - (0.104/(((linelist%freq - 4.67)**2) + 0.341)) + ((-1)*0.04473*((linelist%freq - 5.9)**2) - 0.009779*((linelist%freq - 5.9)**3))) + &
                                           &(1.825*linelist%freq + (1.206/(((linelist%freq - 4.62)**2) + 0.263)) + 0.2130*((linelist%freq - 5.9)**2) + 0.1207*((linelist%freq - 5.9)**3) - 3.090)
                elsewhere((linelist%freq .gt. 3.3) .AND. (linelist%freq .lt. 5.9))
                        linelist%flambda = R* (1.752 - 0.316*linelist%freq - (0.104/(((linelist%freq - 4.67)**2) + 0.341))) &
                                           & + (1.825*linelist%freq + (1.206/(((linelist%freq - 4.62)**2) + 0.263)) - 3.090)
                elsewhere((linelist%freq .gt. 1.1) .AND. (linelist%freq .lt. 3.3))  !Optical & NIR
                        linelist%flambda = R* &
                                          &(1 + (0.17699*(linelist%freq-1.82)) - (0.50447*((linelist%freq-1.82)**2)) - (0.02427*((linelist%freq-1.82)**3)) + (0.72085*((linelist%freq-1.82)**4)) + (0.01979*((linelist%freq-1.82)**5)) - (0.77530*((linelist%freq-1.82)**6)) + (0.32999*((linelist%freq-1.82)**7))) + &
                                          &(1.41338*(linelist%freq-1.82) + 2.28305*((linelist%freq-1.82)**2) + 1.07233*((linelist%freq-1.82)**3) - 5.38434*((linelist%freq-1.82)**4) - 0.62251*((linelist%freq-1.82)**5) + 5.30260*((linelist%freq-1.82)**6) - 2.09002*((linelist%freq-1.82)**7))
                elsewhere(linelist%freq .lt. 1.1) !IR
                        linelist%flambda = R* (0.574 * (linelist%freq**1.61)) + &
                                           & ((-1)*(0.527)*(linelist%freq**1.61))
                elsewhere
!                       print *, "What... is the air-speed velocity of an unladen swallow?"
                        linelist%flambda = 0.d0
                endwhere

                linelist%flambda = (linelist%flambda / ((1.015452*R) + 0.461000)) - 1 ! 1.015452R + 0.461000

        elseif (switch_ext == "P") then !Prevot 1985 SMC law

                where(linelist%freq .gt. 6.72)
                    linelist%flambda = 3.1 + ((3.184*linelist%freq) - 11.45) !Far UV
                elsewhere((linelist%freq .gt. 1.83) .AND. (linelist%freq .lt. 6.72))
                    linelist%flambda = 3.1 + ((2.067*linelist%freq) - 4.110) !Optical/UV
                elsewhere(linelist%freq .lt. 1.83)
                    linelist%flambda = (((1.86*(linelist%freq**2)) - (0.48*(linelist%freq**3)) - (0.1*linelist%freq))) !IR
                elsewhere
!                   print *,"We interrupt this program to annoy you and make things generally irritating."
                    linelist%flambda = 0.d0
                endwhere

                linelist%flambda = (linelist%flambda / 3.242) - 1 ! R+0.142

        elseif (switch_ext == "F") then !Fitzpatrick 1990 galactic law

                where(linelist%freq .lt. 1.83)
                    linelist%flambda = (( (1.86*linelist%freq**2) - (0.48*linelist%freq**3) - (0.1*linelist%freq))) ! IR (4)
                elsewhere((linelist%freq .gt. 1.83) .AND. (linelist%freq .lt. 2.75))
                    linelist%flambda = ((3.1 + (2.56*(linelist%freq-1.83)) - (0.993*(linelist%freq-1.83)**2) )) ! optical (3)
                elsewhere((linelist%freq .gt. 2.75) .AND. (linelist%freq .lt. 3.65))
                    linelist%flambda = ((1.46 + (1.048*linelist%freq) + (1.01/(((linelist%freq-4.60)**2) + 0.280)))) !near UV
                elsewhere((linelist%freq .gt. 3.65) .AND. (linelist%freq .lt. 3.70))
                    linelist%flambda = ((2.19 + (0.848*linelist%freq) + (1.01/(((linelist%freq-4.60)**2) + 0.280)))) !mid UV
                elsewhere((linelist%freq .gt. 3.70) .AND. (linelist%freq .lt. 5.90))
                    linelist%flambda = 3.1 + c1 + c2*linelist%freq + c3*&
                    & (linelist%freq**2)/((((linelist%freq**2) + (x_0**2))**2) +((linelist%freq**2)*(gamma**2)))
                elsewhere(linelist%freq .gt. 5.90)
                    linelist%flambda = 3.1 + c1 + c2*linelist%freq + &
                                       & c3*(linelist%freq**2)/((((linelist%freq**2) + (x_0**2))**2) +((linelist%freq**2)*(gamma**2))) + &
                                       & c4*(0.539*((linelist%freq-5.9)**2.0)) + (0.0564*((linelist%freq-5.9)**3.0))
                elsewhere
!                   print *,"No-one expects the spanish inquisition."
                    linelist%flambda = 0.d0
                endwhere

                linelist%flambda = (linelist%flambda / 3.63) - 1 ! should be R+1.20, why is it same as Howarth galactic? RW

        endif

end subroutine get_flambda

subroutine deredden(lines, extinction)
!this subroutine calls the appropriate dereddening subroutine so that in abundances.f90, we can just call a single routine
        IMPLICIT NONE
        TYPE(line), DIMENSION(:) :: lines
        real(kind=dp) :: extinction

!debugging
#ifdef CO
        print *,"subroutine: deredden"
#endif

        where (lines%intensity .gt. 0.d0 .or. lines%blend_intensity .gt. 0.d0)
          lines%int_dered = lines%intensity * 10**(extinction*lines%flambda)
          lines%blend_int_dered = lines%blend_intensity * 10**(extinction*lines%flambda)
        elsewhere
          lines%int_dered = 0.d0
          lines%blend_int_dered = 0.d0
        endwhere

end subroutine deredden

end module mod_extinction
