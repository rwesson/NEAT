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

!debugging
#ifdef CO
        print *,"subroutine: calc_extinction_coeffs. Te=",temp,", ne=",dens,", weights=",weightha, weighthg, weighthd
#endif

!set up weights

        if (weightha .eq. -1 .and. H_Balmer(3) .gt. 0) then
          weightha = linelist(H_Balmer(3))%intensity
        endif
        if (weighthg .eq. -1 .and. H_Balmer(5) .gt. 0) then
          weighthg = linelist(H_Balmer(5))%intensity
        endif
        if (weighthd .eq. -1 .and. H_Balmer(6) .gt. 0) then
          weighthd = linelist(H_Balmer(6))%intensity
        endif

!Section with interpolations for Balmer ratios for T_e (temp) and N_e (dens)
!interpolating over 6 Te points, 5,7.5,10,12.5,15,20kK and 3 Ne points 10^2, 10^3 10^4 cm^(-3)

        if (H_Balmer(3) .gt. 0 .and. H_Balmer(4) .gt. 0) then
          c1 = log10( ( linelist(H_Balmer(3))%intensity / linelist(H_Balmer(4))%intensity )/ calc_balmer_ratios(temp, dens, 1) )/(-linelist(H_Balmer(3))%flambda)
        else
          c1 = 0.d0
          weightha = 0.d0
        endif

        if (H_Balmer(5) .gt. 0 .and. H_Balmer(4) .gt. 0) then
          c2 = log10( ( linelist(H_Balmer(5))%intensity / linelist(H_Balmer(4))%intensity )/ calc_balmer_ratios(temp, dens, 2) )/(-linelist(H_Balmer(5))%flambda)
        else
          c2=0.d0
          weighthg=0.d0
        endif

        if (H_Balmer(6) .gt. 0 .and. H_Balmer(4) .gt. 0) then
          c3 = log10( ( linelist(H_Balmer(6))%intensity / linelist(H_Balmer(4))%intensity )/ calc_balmer_ratios(temp, dens, 3) )/(-linelist(H_Balmer(6))%flambda)
        else
          c3 = 0.d0
          weighthd = 0.d0
        endif

        if (c1<0) then
!        print *,"Warning: c(Ha) less than zero"
          c1=0
        endif
        if (c2<0) then
!        print *,"Warning: c(Hg) less than zero"
          c2=0
        endif
        if (c3<0) then
!        print *,"Warning: c(Hd) less than zero"
          c3=0
        endif

        if ((weightha + weighthg + weighthd) .gt. 0) then
          meanextinction = (c1*weightha + c2*weighthg + c3*weighthd) / (weightha + weighthg + weighthd)
        else
          meanextinction = 0.d0 ! todo: print a warning here
        endif

end subroutine calc_extinction_coeffs

real(kind=dp) function calc_balmer_ratios(temp, dens, line)
        implicit none
        integer :: i, j, line
        real(kind=dp) :: d1, d2, temp, dens
        real(kind=dp), dimension(6,3,4) :: HS

!debugging
#ifdef CO
        !print *,"function: calc_balmer_ratios"
#endif

        !Ha/Hb
        HS(1,1,:) = (/ 5000., 3.04, 3.02, 3.00/)
        HS(2,1,:) = (/ 7500., 2.93, 2.92, 2.90/)
        HS(3,1,:) = (/10000., 2.86, 2.86, 2.85/)
        HS(4,1,:) = (/12500., 2.82, 2.82, 2.81/)
        HS(5,1,:) = (/15000., 2.79, 2.79, 2.78/)
        HS(6,1,:) = (/20000., 2.75, 2.75, 2.74/)
        !Hg/Hb
        HS(1,2,:) = (/ 5000., 0.458, 0.459, 0.460/)
        HS(2,2,:) = (/ 7500., 0.465, 0.465, 0.466/)
        HS(3,2,:) = (/10000., 0.468, 0.469, 0.469/)
        HS(4,2,:) = (/12500., 0.471, 0.471, 0.472/)
        HS(5,2,:) = (/15000., 0.473, 0.473, 0.473/)
        HS(6,2,:) = (/20000., 0.475, 0.475, 0.476/)
        !Hd/Hb
        HS(1,3,:) = (/ 5000., 0.251, 0.252, 0.253/)
        HS(2,3,:) = (/ 7500., 0.256, 0.256, 0.257/)
        HS(3,3,:) = (/10000., 0.259, 0.259, 0.260/)
        HS(4,3,:) = (/12500., 0.261, 0.261, 0.261/)
        HS(5,3,:) = (/15000., 0.262, 0.263, 0.263/)
        HS(6,3,:) = (/20000., 0.264, 0.264, 0.264/)


        do i = 1,5
                if(temp .ge. HS(i,line,1) .and. temp .lt. HS(i+1,line, 1) )then
                        exit
                endif
        enddo

        do j = 2,3
                if(dens .ge. 10**j .and. dens .le. 10**(j+1) )then
                        exit
                endif
        enddo

        if(temp .lt. HS(1,line,1))then
                if(dens .lt. 10**2)then
                        calc_balmer_ratios = HS(1,line,2)
                        return
                elseif(dens .ge. 10**4)then
                        calc_balmer_ratios = HS(1,line,4)
                        return
                endif
                !print*, j
                calc_balmer_ratios=( (HS(1,line,j+1)-HS(1,line,j)) / (10**(j+1)-10**j) )*(dens-(10**j))+HS(1,line,j)
                return
        elseif(temp .ge. HS(6,line,1))then
                if(dens .lt. 10**2)then
                        calc_balmer_ratios = HS(6,line,2)
                        return
                elseif(dens .gt. 10**4)then
                        calc_balmer_ratios = HS(6,line,4)
                        return
                endif
                calc_balmer_ratios =( ( HS(6,line,j+1)-HS(6,line,j) )/ (10**(j+1)-10**j) )*(dens-(10**j) ) + HS(6,line,j)
                return
        endif

        if(dens .lt. 10**2)then
                calc_balmer_ratios = ((HS(i+1,line,2)-HS(i,line,2))/(HS(i+1,line,1)-HS(i,line,1)))*(temp-HS(i,line,1))+HS(i,line,2)
                return
        elseif(dens .gt. 10**4)then
                calc_balmer_ratios = ((HS(i+1,line,4)-HS(i,line,4))/(HS(i+1,line,1)-HS(i,line,1)))*(temp-HS(i,line,1))+HS(i,line,4)
                return
        endif

        !        print*, j
        !PRINT*, i, j

        d1=((HS(i+1,line,j)-HS(i,line,j))/(HS(i+1,line,1)-HS(i,line,1)))*(temp-HS(i,line,1))+HS(i,line,j)
        d2=((HS(i+1,line,j+1)-HS(i,line,j+1))/(HS(i+1,line,1)-HS(i,line,1)))*(temp-HS(i,line,1))+HS(i,line,j+1)

        !print*, d1, d2

        calc_balmer_ratios =  ( (d2-d1) /(10**(j+1) - 10**j))*(dens - (10**j) ) + d1  !(d1+d2)/2

end function

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
