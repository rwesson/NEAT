module mod_extinction
use mod_abundtypes
implicit none

contains

subroutine calc_extinction_coeffs(H_BS, c1, c2, c3, meanextinction, switch_ext, temp, dens, R)
        IMPLICIT NONE
        TYPE(line), DIMENSION(38) :: H_BS
        DOUBLE PRECISION :: c1, c2, c3, meanextinction, R
        double precision :: fl_ha, fl_hg, fl_hd
        character :: switch_ext
        double precision :: temp, dens

!determine f(lambda) for the balmer lines, depending on the law used
!galactic, howarth

if (switch_ext=="S") then ! howarth galactic
        fl_ha = flambda(dble(10000./6562.77),5)
        fl_hg = flambda(dble(10000./4340.47),4)
        fl_hd = flambda(dble(10000./4101.74),4)
elseif (switch_ext=="H") then ! howarth lmc
        fl_ha = flambdaLMC(dble(10000./6562.77),3)
        fl_hg = flambdaLMC(dble(10000./4340.47),2)
        fl_hd = flambdaLMC(dble(10000./4101.74),2)
elseif (switch_ext=="C") then ! CCM galactic
        fl_ha = flambdaCCM(dble(10000./6562.77),3, R)
        fl_hg = flambdaCCM(dble(10000./4340.47),3, R)
        fl_hd = flambdaCCM(dble(10000./4101.74),3, R)
elseif (switch_ext=="P") then ! Prevot SMC
        fl_ha = flambdaSMC(dble(10000./6562.77),3)
        fl_hg = flambdaSMC(dble(10000./4340.47),2)
        fl_hd = flambdaSMC(dble(10000./4101.74),2)
elseif (switch_ext=="F") then ! Fitzpatrick galactic
        fl_ha = flambdaFitz(dble(10000./6562.77),1)
        fl_hg = flambdaFitz(dble(10000./4340.47),2)
        fl_hd = flambdaFitz(dble(10000./4101.74),2)
else
        print *, "He's a very naughty boy!"
        fl_ha = 0.0
        fl_hg = 0.0
        fl_hd = 0.0
endif

!Section with interpolations for Balmer ratios for T_e (temp) and N_e (dens)
!interpolating over 6 Te points, 5,7.5,10,12.5,15,20kK and 3 Ne points 10^2, 10^3 10^4 cm^(-3)

if (H_BS(1)%intensity .gt. 0 .and. H_BS(2)%intensity .gt. 0) then
        c1 = log10( ( DBLE(H_BS(1)%intensity) / DBLE(H_BS(2)%intensity) )/ calc_balmer_ratios(temp, dens, 1)    )/(-fl_ha)
else
        c1 = 0.0
endif

if (H_BS(3)%intensity .gt. 0 .and. H_BS(2)%intensity .gt. 0) then
        c2 = log10( ( DBLE(H_BS(3)%intensity) / DBLE(H_BS(2)%intensity) )/ calc_balmer_ratios(temp, dens, 2) )/(-fl_hg)
else
        c2=0.0
endif

if (H_BS(4)%intensity .gt. 0 .and. H_BS(2)%intensity .gt. 0) then
        c3 = log10( ( DBLE(H_BS(4)%intensity) / DBLE(H_BS(2)%intensity) )/ calc_balmer_ratios(temp, dens, 3) )/(-fl_hd)
else
        c3 = 0.0
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

        meanextinction = (c1*H_BS(1)%intensity + c2*H_BS(3)%intensity + c3*H_BS(4)%intensity) / (H_BS(1)%intensity + H_BS(3)%intensity + H_BS(4)%intensity)

end subroutine calc_extinction_coeffs

double precision function calc_balmer_ratios(temp, dens, line)
        implicit none
        integer :: i, j, line
        double precision :: d1, d2, temp, dens
        double precision, dimension(6,3,4) :: HS

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
        end do

        do j = 2,3
                if(dens .ge. 10**j .and. dens .le. 10**(j+1) )then
                        exit
                endif
        end do

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

!-------SEATON GALACTIC LAW-------------------------------!

double precision function flambda(X,switch)
        IMPLICIT NONE
        DOUBLE PRECISION :: X
        INTEGER :: switch
        !Howarth 1983 Galactic + Seaton 1979
        if(switch == 1) then
            flambda = ((16.07 - (3.20 * X) + (0.2975*X**2))) !far UV
        elseif(switch == 2) then
            flambda = ((2.19 + (0.848*X) + (1.01/(((X-4.60)**2) + 0.280)))) !mid UV
        elseif(switch == 3) then
            flambda = ((1.46 + (1.048*X) + (1.01/(((X-4.60)**2) + 0.280)))) !near UV
        elseif(switch == 4) then
            flambda = ((3.1 + (2.56*(X-1.83)) - (0.993*(X-1.83)**2) )) ! optical (3)
        elseif(switch == 5) then
            flambda = (( (1.86*X**2) - (0.48*X**3) - (0.1*X))) ! IR (4)
        else
            print *, "He's not the messiah!"
            flambda = 0.0
        endif

        flambda = (flambda / 3.63) - 1 ! R+0.53

end function

subroutine deredden(lines, number, m_ext)
        IMPLICIT NONE
        TYPE(line), DIMENSION(:) :: lines
        INTEGER :: number
        DOUBLE PRECISION :: m_ext, fl
        INTEGER :: i

        do i = 1,number

                if (lines(i)%wavelength<500) then ! allows line list to contain IR lines with wavelengths in microns.
                        lines(i)%freq = DBLE(1) / DBLE(lines(i)%wavelength)
                else
                        lines(i)%freq = DBLE(10000) / DBLE(lines(i)%wavelength)
                endif

                if( (lines(i)%freq .gt. 7.14) .AND. (lines(i)%freq .lt. 10.0))then !far UV
                        fl = flambda(lines(i)%freq, 1)
                        lines(i)%int_dered = lines(i)%intensity * 10**(m_ext*fl)

                elseif((lines(i)%freq .gt. 3.65) .AND. (lines(i)%freq .lt. 7.14))then ! mid UV
                        fl = flambda(lines(i)%freq, 2)
                        lines(i)%int_dered = lines(i)%intensity * 10**(m_ext*fl)

                elseif((lines(i)%freq .gt. 2.75) .AND. (lines(i)%freq .lt. 3.65))then ! near UV
                        fl = flambda(lines(i)%freq, 3)
                        lines(i)%int_dered = lines(i)%intensity * 10**(m_ext*fl)

                elseif((lines(i)%freq .gt. 1.83) .AND. (lines(i)%freq .lt. 2.75))then ! optical
                        fl = flambda(lines(i)%freq, 4)
                        lines(i)%int_dered = lines(i)%intensity * 10**(m_ext*fl)

                elseif(lines(i)%freq .lt. 1.83)then !IR
                        fl = flambda(lines(i)%freq, 5)
                        lines(i)%int_dered = lines(i)%intensity * 10**(m_ext*fl)

                endif

        end do

end subroutine

!-------HOWARTH LMC LAW-----------------------------------!

double precision function flambdaLMC(X,switch)
        IMPLICIT NONE
        DOUBLE PRECISION :: X
        INTEGER :: switch
        !Howarth 1983 LMC
        if(switch == 1) then
            flambdaLMC = ( (3.1 - 0.236 + (0.462*X) + (0.105*(X**2)) + (0.454/((X-4.557)**2 + 0.293))) ) !UV (1)
        elseif(switch == 2) then
            flambdaLMC = ((3.1 + (2.04*(X - 1.83)) + 0.094*(X - 1.83)**2)) !Optical (2)
        elseif(switch == 3) then
            flambdaLMC = (((1.86*(X**2)) - (0.48*(X**3)) - (0.1*X))) !IR (4)
        else
            print *,"Your mother was a hamster and your father smelt of elderberries."
            flambdaLMC = 0.0
        endif

        flambdaLMC = (flambdaLMC / 3.57) - 1 ! R+0.47

end function

subroutine deredden_LMC(lines, number, m_ext)
        IMPLICIT NONE
        TYPE(line), DIMENSION(:) :: lines
        INTEGER :: number
        DOUBLE PRECISION :: m_ext, fl
        INTEGER :: i

        do i = 1,number

                if (lines(i)%wavelength<500) then ! allows line list to contain IR lines with wavelengths in microns.
                        lines(i)%freq = DBLE(1) / DBLE(lines(i)%wavelength)
                else
                        lines(i)%freq = DBLE(10000) / DBLE(lines(i)%wavelength)
                endif

                if(lines(i)%freq .gt. 2.75)then ! UV
                        fl = flambdaLMC(lines(i)%freq, 1)
                        lines(i)%int_dered = lines(i)%intensity * 10**(m_ext*fl)

                elseif((lines(i)%freq .gt. 1.83) .AND. (lines(i)%freq .lt. 2.75))then ! optical
                        fl = flambdaLMC(lines(i)%freq, 2)
                        lines(i)%int_dered = lines(i)%intensity * 10**(m_ext*fl)

                elseif(lines(i)%freq .lt. 1.83)then !IR
                        fl = flambdaLMC(lines(i)%freq, 3)
                        lines(i)%int_dered = lines(i)%intensity * 10**(m_ext*fl)

                endif

        end do

end subroutine

!-------CCM GALACTIC LAW----------------------------------!

double precision function flambdaCCM(X,switch, R)
        IMPLICIT NONE
        DOUBLE PRECISION :: X, a, b, y, Fa, Fb, R
        INTEGER :: switch
        !CCM 1989 Galactic
        if(switch == 1) then !far UV
                a = 0.137*((X - 8)**2) - 1.073 - 0.628*(X - 8) - 0.070*((X - 8)**3)
                b = 13.670 + 4.257*(X - 8) - 0.420*((X - 8)**2) + 0.374*((X - 8)**3)
        elseif(switch == 2) then !UV
                Fa = 0
                Fb = 0
                if (X .gt. 5.9) then
                Fa = (-1)*0.04473*((X - 5.9)**2) - 0.009779*((X - 5.9)**3)
                Fb = 0.2130*((X - 5.9)**2) + 0.1207*((X - 5.9)**3)
                endif
                a = 1.752 - 0.316*X - (0.104/(((X - 4.67)**2) + 0.341)) +Fa
                b = 1.825*X + (1.206/(((X - 4.62)**2) + 0.263)) + Fb - 3.090
        elseif(switch == 3) then  !Optical & NIR
                y = X-1.82
                a = 1 + (0.17699*y) - (0.50447*(y**2)) - (0.02427*(y**3)) + (0.72085*(y**4)) + (0.01979*(y**5)) - (0.77530*(y**6)) + (0.32999*(y**7))
                b = 1.41338*y + 2.28305*(y**2) + 1.07233*(y**3) - 5.38434*(y**4) - 0.62251*(y**5) + 5.30260*(y**6) - 2.09002*(y**7)
        elseif(switch == 4) then !IR
                a = 0.574 * (X**1.61)
                b = (-1)*(0.527)*(X**1.61)
        else
                print *, "What... is the air-speed velocity of an unladen swallow?"
                a = 0.0
                b = 0.0
        endif

        flambdaCCM = R*a + b
        flambdaCCM = (flambdaCCM / ((1.015452*R) + 0.461000)) - 1 ! 1.015452R + 0.461000

end function

subroutine deredden_CCM(lines, number, m_ext, R)
        IMPLICIT NONE
        TYPE(line), DIMENSION(:) :: lines
        INTEGER :: number
        DOUBLE PRECISION :: m_ext, fl, R
        INTEGER :: i

        do i = 1,number

                if (lines(i)%wavelength<500) then ! allows line list to contain IR lines with wavelengths in microns.
                        lines(i)%freq = DBLE(1) / DBLE(lines(i)%wavelength)
                else
                        lines(i)%freq = DBLE(10000) / DBLE(lines(i)%wavelength)
                endif

                if(lines(i)%freq .gt. 8)then ! Far UV
                        fl = flambdaCCM(lines(i)%freq, 1, R)
                        lines(i)%int_dered = lines(i)%intensity * 10**(m_ext*fl)

                elseif((lines(i)%freq .gt. 3.3) .AND. (lines(i)%freq .lt. 8))then ! UV
                        fl = flambdaCCM(lines(i)%freq, 2, R)
                        lines(i)%int_dered = lines(i)%intensity * 10**(m_ext*fl)

                elseif((lines(i)%freq .gt. 1.1) .AND. (lines(i)%freq .lt. 3.3))then !optical & NIR
                        fl = flambdaCCM(lines(i)%freq, 3, R)
                        lines(i)%int_dered = lines(i)%intensity * 10**(m_ext*fl)

                elseif(lines(i)%freq .lt. 1.1)then !IR
                        fl = flambdaCCM(lines(i)%freq, 4, R)
                        lines(i)%int_dered = lines(i)%intensity * 10**(m_ext*fl)

                endif

        end do

end subroutine

!-------PREVOT SMC LAW------------------------------------!

double precision function flambdaSMC(X,switch)
        IMPLICIT NONE
        DOUBLE PRECISION :: X
        INTEGER :: switch
        !Prevot 1984 SMC
        if(switch == 1) then
            flambdaSMC = 3.1 + ((3.184*X) - 11.45) !Far UV
        elseif(switch == 2) then
            flambdaSMC = 3.1 + ((2.067*X) - 4.110) !Optical/UV
        elseif(switch == 3) then
            flambdaSMC = (((1.86*(X**2)) - (0.48*(X**3)) - (0.1*X))) !IR
        else
            print *,"We interrupt this program to annoy you and make things generally irritating."
            flambdaSMC = 0.0
        endif

        flambdaSMC = (flambdaSMC / 3.242) - 1 ! R+0.142

end function

subroutine deredden_SMC(lines, number, m_ext)
        IMPLICIT NONE
        TYPE(line), DIMENSION(:) :: lines
        INTEGER :: number
        DOUBLE PRECISION :: m_ext, fl
        INTEGER :: i

        do i = 1,number

                if (lines(i)%wavelength<500) then ! allows line list to contain IR lines with wavelengths in microns.
                        lines(i)%freq = DBLE(1) / DBLE(lines(i)%wavelength)
                else
                        lines(i)%freq = DBLE(10000) / DBLE(lines(i)%wavelength)
                endif

                if(lines(i)%freq .gt. 6.72)then ! Far UV
                        fl = flambdaSMC(lines(i)%freq, 1)
                        lines(i)%int_dered = lines(i)%intensity * 10**(m_ext*fl)

                elseif((lines(i)%freq .gt. 1.83) .AND. (lines(i)%freq .lt. 6.72))then ! optical/UV
                        fl = flambdaSMC(lines(i)%freq, 2)
                        lines(i)%int_dered = lines(i)%intensity * 10**(m_ext*fl)

                elseif(lines(i)%freq .lt. 1.83)then !IR
                        fl = flambdaSMC(lines(i)%freq, 3)
                        lines(i)%int_dered = lines(i)%intensity * 10**(m_ext*fl)

                endif

        end do

end subroutine

!-------Fitzpatrick Galactic law--------------------------!


double precision function flambdaFitz(X,switch) !based on FM90 with values taken from Fitxpatrick 1992
        IMPLICIT NONE
        DOUBLE PRECISION :: X, x_0, c1, c2, c3, c4, gamma, D, F
        INTEGER :: switch
       c1 = -0.38
       c2 =  0.74
       c3 =  3.96
       c4 = 0
       x_0 = 4.596
       gamma = 1.051
       F = 0
       D = (X**2)/((((X**2) + (x_0**2))**2) +((X**2)*(gamma**2)))

        if(switch == 1) then
            flambdaFitz = (( (1.86*X**2) - (0.48*X**3) - (0.1*X))) ! IR (4)
        elseif(switch == 2) then
            flambdaFitz = ((3.1 + (2.56*(X-1.83)) - (0.993*(X-1.83)**2) )) ! optical (3)
        elseif(switch == 3) then
            flambdaFitz = ((1.46 + (1.048*X) + (1.01/(((X-4.60)**2) + 0.280)))) !near UV
        elseif(switch == 4) then
            flambdaFitz = ((2.19 + (0.848*X) + (1.01/(((X-4.60)**2) + 0.280)))) !mid UV
        elseif(switch == 5) then
                        flambdaFitz = 3.1 + c1 + c2*X + c3*D
        elseif(switch == 6) then
                c4 = 0.26
                F=(0.539*((X-5.9)**2.0)) + (0.0564*((X-5.9)**3.0))
                flambdaFitz = 3.1 + c1 + c2*X + c3*D + c4*F
        else
                print *,"No-one expects the spanish inquisition."
                flambdaFitz = 0.0
        endif

        flambdaFitz = (flambdaFitz / 3.63) - 1 ! should be R+1.20, why is it same as Howarth galactic? RW

end function

subroutine deredden_Fitz(lines, number, m_ext) !Uses Seaton/Howarth for IR, optical and near UV, as fit for Fitzpatrick law is poor at wavelengths longer than 2700 Angstroms
        IMPLICIT NONE
        TYPE(line), DIMENSION(:) :: lines
        INTEGER :: number
        DOUBLE PRECISION :: m_ext, fl
        INTEGER :: i

        do I = 1, number

                if (lines(i)%wavelength<500) then ! allows line list to contain IR lines with wavelengths in microns.
                        lines(i)%freq = DBLE(1) / DBLE(lines(i)%wavelength)
                else
                        lines(i)%freq = DBLE(10000) / DBLE(lines(i)%wavelength)
                endif

                if (lines(i)%freq .lt. 1.83) then
                        fl = flambdaFitz(lines(i)%freq, 1)
                        lines(i)%int_dered = lines(i)%intensity * 10**(m_ext*fl)
                elseif((lines(i)%freq .gt. 1.83) .AND. (lines(i)%freq .lt. 2.75))then ! optical
                        fl = flambdaFitz(lines(i)%freq, 2)
                        lines(i)%int_dered = lines(i)%intensity * 10**(m_ext*fl)
                elseif((lines(i)%freq .gt. 2.75) .AND. (lines(i)%freq .lt. 3.65))then ! near UV
                        fl = flambdaFitz(lines(i)%freq, 3)
                        lines(i)%int_dered = lines(i)%intensity * 10**(m_ext*fl)
                elseif((lines(i)%freq .gt. 3.65) .AND. (lines(i)%freq .lt. 3.70))then ! mid UV
                        fl = flambdaFitz(lines(i)%freq, 4)
                        lines(i)%int_dered = lines(i)%intensity * 10**(m_ext*fl)
                elseif((lines(i)%freq .gt. 3.70) .AND. (lines(i)%freq .lt. 5.90))then ! mid UV
                        fl = flambdaFitz(lines(i)%freq, 5)
                        lines(i)%int_dered = lines(i)%intensity * 10**(m_ext*fl)
                elseif(lines(i)%freq .gt. 5.90)then ! far UV
                        fl = flambdaFitz(lines(i)%freq, 6)
                        lines(i)%int_dered = lines(i)%intensity * 10**(m_ext*fl)
                endif

        end do
end subroutine

end module mod_extinction
