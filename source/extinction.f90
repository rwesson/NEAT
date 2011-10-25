module mod_extinction
use mod_abundtypes
implicit none

contains

subroutine calc_extinction_coeffs(H_BS, c1, c2, c3, meanextinction, switch_ext)
        TYPE(line), DIMENSION(4) :: H_BS
        DOUBLE PRECISION :: c1, c2, c3, meanextinction
        double precision :: fl_ha, fl_hg, fl_hd
        character :: switch_ext

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
        fl_ha = flambdaCCM(dble(10000./6562.77),3)
        fl_hg = flambdaCCM(dble(10000./4340.47),3)
        fl_hd = flambdaCCM(dble(10000./4101.74),3)
elseif (switch_ext=="P") then ! Prevot SMC
        fl_ha = flambdaSMC(dble(10000./6562.77),3)
        fl_hg = flambdaSMC(dble(10000./4340.47),2)
        fl_hd = flambdaSMC(dble(10000./4101.74),2)
elseif (switch_ext=="F") then ! Fitzpatrick galactic
        fl_ha = flambdaFitz(dble(10000./6562.77),1)
        fl_hg = flambdaFitz(dble(10000./4340.47),2)
        fl_hd = flambdaFitz(dble(10000./4101.74),2)
endif

if (H_BS(1)%intensity .gt. 0 .and. H_BS(2)%intensity .gt. 0) then
        c1 = log10( ( DBLE(H_BS(1)%intensity) / DBLE(H_BS(2)%intensity) )/2.850 )/(-fl_ha)
else
        c1 = 0.0
endif
if (H_BS(3)%intensity .gt. 0 .and. H_BS(2)%intensity .gt. 0) then
        c2 = log10( ( DBLE(H_BS(3)%intensity) / DBLE(H_BS(2)%intensity) )/0.469 )/(-fl_hg)
else
        c2=0.0
endif
if (H_BS(4)%intensity .gt. 0 .and. H_BS(2)%intensity .gt. 0) then
        c3 = log10( ( DBLE(H_BS(4)%intensity) / DBLE(H_BS(2)%intensity) )/0.259 )/(-fl_hd)
else
        c3 = 0.0
endif

        meanextinction = (c1*H_BS(1)%intensity + c2*H_BS(3)%intensity + c3*H_BS(4)%intensity) / (H_BS(1)%intensity + H_BS(3)%intensity + H_BS(4)%intensity)   

end subroutine

!-------SEATON GALACTIC LAW-------------------------------!

double precision function flambda(X,switch)
        DOUBLE PRECISION :: X
        INTEGER :: switch
        !Howarth 1983 Galactic + Seaton 1979
        if(switch == 1) flambda = ((16.07 - (3.20 * X) + (0.2975*X**2))) !far UV
        if(switch == 2) flambda = ((2.19 + (0.848*X) + (1.01/(((X-4.60)**2) + 0.280)))) !mid UV
        if(switch == 3) flambda = ((1.46 + (1.048*X) + (1.01/(((X-4.60)**2) + 0.280)))) !near UV
        if(switch == 4) flambda = ((3.1 + (2.56*(X-1.83)) - (0.993*(X-1.83)**2) )) ! optical (3)
        if(switch == 5) flambda = (( (1.86*X**2) - (0.48*X**3) - (0.1*X))) ! IR (4)

        flambda = (flambda / 3.63) - 1 ! R+0.53

end function        

subroutine deredden(lines, number, m_ext)
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

subroutine deredden_O(O_red, O_dered, m_ext)
        REAL*8, DIMENSION(3,335) :: O_red
        REAL*8, DIMENSION(2,335) :: O_dered
        DOUBLE PRECISION :: m_ext, fl
        INTEGER :: I
        REAL*8 :: X

        do I = 1, 335
            X =  DBLE(DBLE(10000)/DBLE(O_red(1,I)))
            if (X .lt. 1.83) then 
                fl = flambda(X, 5) 
                O_dered(:,I) = (/ O_red(1,I), O_red(2,I)*(10**(m_ext*fl)) /)
            elseif (X .gt. 1.83 .and. X .lt. 2.75) then
                fl = flambda(X, 4)
                O_dered(:,I) = (/ O_red(1,I), O_red(2,I)*(10**(m_ext*fl)) /)
            elseif (X .gt. 2.75 .and. X .lt. 3.63) then
                fl = flambda(X, 3)
                O_dered(:,I) = (/ O_red(1,I), O_red(2,I)*(10**(m_ext*fl)) /)
            endif 
        end do

end subroutine  

!-------HOWARTH LMC LAW-----------------------------------!

double precision function flambdaLMC(X,switch)
        DOUBLE PRECISION :: X
        INTEGER :: switch
        !Howarth 1983 LMC
        if(switch == 1) flambdaLMC = ( (3.1 - 0.236 + (0.462*X) + (0.105*(X**2)) + (0.454/((X-4.557)**2 + 0.293))) ) !UV (1)
        if(switch == 2) flambdaLMC = ((3.1 + (2.04*(X - 1.83)) + 0.094*(X - 1.83)**2)) !Optical (2)
        if(switch == 3) flambdaLMC = (((1.86*(X**2)) - (0.48*(X**3)) - (0.1*X))) !IR (4)

       flambdaLMC = (flambdaLMC / 3.57) - 1 ! R+0.47

end function

subroutine deredden_LMC(lines, number, m_ext)
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

subroutine deredden_O_LMC(O_red, O_dered, m_ext)
        REAL*8, DIMENSION(3,335) :: O_red
        REAL*8, DIMENSION(2,335) :: O_dered
        DOUBLE PRECISION :: m_ext, fl
        INTEGER :: I
        REAL*8 :: X

        do I = 1, 335
            X =  DBLE(DBLE(10000)/DBLE(O_red(1,I)))
            if (X .lt. 1.83) then 
                fl = flambdaLMC(X, 3) 
                O_dered(:,I) = (/ O_red(1,I), O_red(2,I)*(10**(m_ext*fl)) /)
            elseif (X .gt. 1.83 .and. X .lt. 2.75) then
                fl = flambdaLMC(X, 2)
                O_dered(:,I) = (/ O_red(1,I), O_red(2,I)*(10**(m_ext*fl)) /)
            elseif (X .gt. 2.75) then
                fl = flambdaLMC(X, 1)
                O_dered(:,I) = (/ O_red(1,I), O_red(2,I)*(10**(m_ext*fl)) /)
            endif 
        end do

end subroutine  

!-------CCM GALACTIC LAW----------------------------------!

double precision function flambdaCCM(X,switch)
        DOUBLE PRECISION :: X, a, b, y, Fa, Fb
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
        endif

        flambdaCCM = 3.1*a + b
        flambdaCCM = (flambdaCCM / 3.61) - 1 ! 1.015452R + 0.461000

end function

subroutine deredden_CCM(lines, number, m_ext)
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

                if(lines(i)%freq .gt. 8)then ! Far UV
                        fl = flambdaCCM(lines(i)%freq, 1)
                        lines(i)%int_dered = lines(i)%intensity * 10**(m_ext*fl) 

                elseif((lines(i)%freq .gt. 3.3) .AND. (lines(i)%freq .lt. 8))then ! UV
                        fl = flambdaCCM(lines(i)%freq, 2)
                        lines(i)%int_dered = lines(i)%intensity * 10**(m_ext*fl) 

                elseif((lines(i)%freq .gt. 1.1) .AND. (lines(i)%freq .lt. 3.3))then !optical & NIR
                        fl = flambdaCCM(lines(i)%freq, 3)
                        lines(i)%int_dered = lines(i)%intensity * 10**(m_ext*fl) 

                elseif(lines(i)%freq .lt. 1.1)then !IR
                        fl = flambdaCCM(lines(i)%freq, 4)
                        lines(i)%int_dered = lines(i)%intensity * 10**(m_ext*fl) 

                endif

        end do

end subroutine

subroutine deredden_O_CCM(O_red, O_dered, m_ext)
        REAL*8, DIMENSION(3,335) :: O_red
        REAL*8, DIMENSION(2,335) :: O_dered
        DOUBLE PRECISION :: m_ext, fl
        INTEGER :: I
        REAL*8 :: X

        do I = 1, 335
            X =  DBLE(DBLE(10000)/DBLE(O_red(1,I)))
            if (X .lt. 1.1) then 
                fl = flambdaCCM(X, 4) 
                O_dered(:,I) = (/ O_red(1,I), O_red(2,I)*(10**(m_ext*fl)) /)
            elseif (X .gt. 1.1 .and. X .lt. 3.3) then
                fl = flambda(X, 3)
                O_dered(:,I) = (/ O_red(1,I), O_red(2,I)*(10**(m_ext*fl)) /)
            elseif (X .gt. 3.3 .and. X .lt. 8) then
                fl = flambda(X, 2)
                O_dered(:,I) = (/ O_red(1,I), O_red(2,I)*(10**(m_ext*fl)) /)
            endif 
        end do

end subroutine  

!-------PREVOT SMC LAW------------------------------------!

double precision function flambdaSMC(X,switch)
        DOUBLE PRECISION :: X
        INTEGER :: switch
        !Prevot 1984 SMC
        if(switch == 1) flambdaSMC = 3.1 + ((3.184*X) - 11.45) !Far UV
        if(switch == 2) flambdaSMC = 3.1 + ((2.067*X) - 4.110) !Optical/UV
        if(switch == 3) flambdaSMC = (((1.86*(X**2)) - (0.48*(X**3)) - (0.1*X))) !IR

       flambdaSMC = (flambdaSMC / 3.242) - 1 ! R+0.142

end function

subroutine deredden_SMC(lines, number, m_ext)
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

subroutine deredden_O_SMC(O_red, O_dered, m_ext)
        REAL*8, DIMENSION(3,335) :: O_red
        REAL*8, DIMENSION(2,335) :: O_dered
        DOUBLE PRECISION :: m_ext, fl
        INTEGER :: I
        REAL*8 :: X

        do I = 1, 335
            X =  DBLE(DBLE(10000)/DBLE(O_red(1,I)))
            if (X .lt. 1.83) then 
                fl = flambdaSMC(X, 3) 
                O_dered(:,I) = (/ O_red(1,I), O_red(2,I)*(10**(m_ext*fl)) /)
            elseif (X .gt. 1.83 .and. X .lt. 6.72) then
                fl = flambdaSMC(X, 2)
                O_dered(:,I) = (/ O_red(1,I), O_red(2,I)*(10**(m_ext*fl)) /)
            elseif (X .gt. 6.72) then
                fl = flambdaSMC(X, 1)
                O_dered(:,I) = (/ O_red(1,I), O_red(2,I)*(10**(m_ext*fl)) /)
            endif 
        end do

end subroutine  


!-------Fitzpatrick Galactic law--------------------------!


double precision function flambdaFitz(X,switch) !based on FM90 with values taken from Fitxpatrick 1992
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

        if(switch == 1) flambdaFitz = (( (1.86*X**2) - (0.48*X**3) - (0.1*X))) ! IR (4)
        if(switch == 2) flambdaFitz = ((3.1 + (2.56*(X-1.83)) - (0.993*(X-1.83)**2) )) ! optical (3)
        if(switch == 3) flambdaFitz = ((1.46 + (1.048*X) + (1.01/(((X-4.60)**2) + 0.280)))) !near UV
        if(switch == 4) flambdaFitz = ((2.19 + (0.848*X) + (1.01/(((X-4.60)**2) + 0.280)))) !mid UV
        if(switch == 5) then
                        flambdaFitz = 3.1 + c1 + c2*X + c3*D
        elseif(switch == 6) then
                c4 = 0.26
                F=(0.539*((X-5.9)**2.0)) + (0.0564*((X-5.9)**3.0))
                flambdaFitz = 3.1 + c1 + c2*X + c3*D + c4*F
        endif

        flambdaFitz = (flambdaFitz / 3.63) - 1 ! should be R+1.20, why is it same as Howarth galactic? RW

end function

subroutine deredden_Fitz(lines, number, m_ext) !Uses Seaton/Howarth for IR, optical and near UV, as fit for Fitzpatrick law is poor at wavelengths longer than 2700 Angstroms
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
