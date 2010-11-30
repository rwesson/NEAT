!INPUT: CALL ./PROGRAM OPTICAL.DAT IR.DAT UV.DAT
program abundances
use mod_abundmaths
use mod_abundtypes
use mod_diagnostics
use mod_getabunds
use mod_abundIO
use mod_helium
use mod_recombination_lines

implicit none

        INTEGER :: count, Iint, fIL, i, j, ion_no1, ion_no2, ion_no3, ion_no4
        INTEGER :: opt
        CHARACTER*20 :: fname1, fname2, fname3, fnameIL, blank
        CHARACTER*8 :: lion
        DOUBLE PRECISION :: normalise, oiiNratio, oiiDens, oiiiTratio, oiiiTemp, oiiiIRNratio, oiiiIRTratio, niiTratio, niiTemp, arivNratio, arivDens, cliiiNratio, cliiiDens, siiNratio, siiDens, siiTratio, siiTemp, oiiTratio, oiiTemp, neiiiTratio, neiiiIRTratio, neiiiTemp, neiiiIRTemp, abund, meandensity, meantemp
        DOUBLE PRECISION :: ciiiNratio,neivNratio,nevTratio,siiiTratio,ariiiTratio,arvTratio,lowtemp,lowdens,medtemp,ciiidens,meddens,siiitemp,ariiitemp,hightemp,neivdens,highdens,arvtemp,nevtemp
        DOUBLE PRECISION :: oiiRLabund, niiRLabund, ciiRLabund, neiiRLabund, ciiiRLabund, niiiRLabund, RLabundtemp, weight
        DOUBLE PRECISION :: ciiiCELabund, niiCELabund, niiiIRCELabund, niiiUVCELabund, oiiCELabund, oiiiCELabund, oiiiIRCELabund, neiiIRCELabund, neiiiIRCELabund, neiiiCELabund, neivCELabund, siiCELabund, siiiCELabund, siiiIRCELabund, sivIRCELabund, cliiiCELabund, ariiiCELabund, arivCELabund, ariiiIRCELabund, nivCELabund, niCELabund, niiiCELabund, ciiCELabund, civCELabund, nvCELabund, nevCELabund, arvCELabund, CELabundtemp
        DOUBLE PRECISION :: CELicfO, CELicfC, CELicfN, CELicfNe, CELicfAr, CELicfS
        DOUBLE PRECISION :: RLicfO, RLicfC, RLicfN, RLicfNe
        DOUBLE PRECISION :: CabundRL, CabundCEL, NabundRL, NabundCEL, OabundRL, OabundCEL, NeabundRL, NeabundCEL, SabundCEL, ArabundCEL
        DOUBLE PRECISION :: adfC, adfN, adfO, adfNe
        DOUBLE PRECISION :: adfC2plus, adfN2plus, adfO2plus, adfNe2plus
        DOUBLE PRECISION :: c1, c1_err, c2, c2_err, c3, c3_err, meanextinction, cerror, fl, ratob, tempi, temp, temp2
        REAL :: heiabund,heiiabund,Hetotabund
        REAL*8 :: HW
                
        DOUBLE PRECISION, DIMENSION(2) :: conditions 
        REAL*8, DIMENSION(3,335) :: O, IR, UV, O_dered
        REAL*8 :: result       
 
        TYPE(line), DIMENSION(51) :: ILs
        TYPE(line), DIMENSION(4) :: H_BS
        TYPE(line), DIMENSION(4) :: He_lines

! recombination line variables

        TYPE RLabund
           CHARACTER*7 :: Multiplet
           REAL*8 :: Abundance
        END TYPE

        TYPE (RLabund), DIMENSION(20) :: oiimultiplets
        TYPE (RLabund), DIMENSION(20) :: niimultiplets

        O = 0

        !file reading stuff
        
        CALL getarg(1,fname1) 
        CALL getarg(2,fname2)
        CALL getarg(3,fname3)
        
        !reading in Rogers "important" lines list 
        
        CALL read_ilines(ILs, Iint) 
        CALL fileread(O, IR, UV, fname1, fname2, fname3) ! see above 
        CALL element_assign(ILs, O, IR, UV, Iint)

        !dereddening

        !first lets find some hydrogen lines 
        CALL get_H(H_BS, O) 
        call get_He(He_lines, O)

        !aside: normalisation check, correction
        
        if(H_BS(2)%intensity .ne. 100)then
                normalise =  DBLE(100) / DBLE(H_BS(2)%intensity)
                do i = 1, Iint !normalising important ions
                        ILs(i)%intensity = ILs(i)%intensity * normalise 
                end do
                do i = 1, 4 !normalising balmer series
                        H_BS(i)%intensity = H_BS(i)%intensity*normalise 
                end do
                do i = 1,4 !normalise helium
                        He_lines(i)%intensity = He_lines(i)%intensity * normalise
                end do
        endif

        print *, ""
        print *, "Extinction"
        print *, "----------"
        
        CALL calc_extinction_coeffs(H_BS, c1, c2, c3, c1_err, c2_err, c3_err, cerror, meanextinction)

        !need to write output/ input stuff so user can insert own c(Hb)
        !assume we go on with calculated extinctions
                
        !actual dereddening
        
        CALL deredden(ILs, Iint, meanextinction)        
        CALL deredden(H_BS, 4, meanextinction)
        call deredden(He_lines, 4, meanextinction)
        CALL deredden_O(O, O_dered, meanextinction)

        print "(1X,A17,F4.2,A4,F4.2)","Ha/Hb => c(Hb) = ",c1," +- ",c1_err
        print "(1X,A17,F4.2,A4,F4.2)","Hg/Hb => c(Hb) = ",c2," +- ",c2_err
        print "(1X,A17,F4.2,A4,F4.2)","Hd/Hb => c(Hb) = ",c3," +- ",c3_err
   
        PRINT "(1X,A13,F4.2,A4,F4.2)", "Mean c(Hb) = ",meanextinction," +- ",cerror

!diagnostics
!ciii ratio

      ion_no1 = get_ion("ciii1909   ", ILs)
      ion_no2 = get_ion("ciii1907   ", ILs)

      if((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0))then
              ciiiNratio = DBLE(ILs(ion_no1)%int_dered) / DBLE(ILs(ion_no2)%int_dered)
      else
              ciiiNratio = 0.0
      endif

!oii ratio

      ion_no1 = get_ion("oii3729    ", ILs)
      ion_no2 = get_ion("oii3726    ", ILs)

      if((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0))then
              oiiNratio = DBLE(ILs(ion_no1)%int_dered) / DBLE(ILs(ion_no2)%int_dered)
      else
              oiiNratio = 0.0
      endif

! neiv ratio

      ion_no1 = get_ion("neiv2425   ", ILs)
      ion_no2 = get_ion("neiv2423   ", ILs)

      if((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0))then
              neivNratio = DBLE(ILs(ion_no1)%int_dered) / DBLE(ILs(ion_no2)%int_dered)
      else
              neivNratio = 0.0
      endif

!s ii ratio

      ion_no1 = get_ion("sii6731    ", ILs)
      ion_no2 = get_ion("sii6716    ", ILs)

      if((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0))then
              siiNratio = DBLE(ILs(ion_no1)%int_dered) / DBLE(ILs(ion_no2)%int_dered)
      else
              siiNratio = 0.0
      endif

! Cl iii ratio

      ion_no1 = get_ion("cliii5537  ", ILs)
      ion_no2 = get_ion("cliii5517  ", ILs)

      if((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0))then
              cliiiNratio = DBLE(ILs(ion_no1)%int_dered) / DBLE(ILs(ion_no2)%int_dered)
      else
              cliiiNratio = 0.0
      endif

! Ar iv ratio

      ion_no1 = get_ion("ariv4740   ", ILs)
      ion_no2 = get_ion("ariv4711   ", ILs)

      if((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0))then
              arivNratio = DBLE(ILs(ion_no1)%int_dered) / DBLE(ILs(ion_no2)%int_dered)
      else
              arivNratio = 0.0
      endif

! temperature ratios:

! N II

      ion_no1 = get_ion("nii6548    ", ILs)
      ion_no2 = get_ion("nii6584    ", ILs)
      ion_no3 = get_ion("nii5754    ", ILs)

      if(((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0)) .and. (ILs(ion_no3)%int_dered .gt. 0 ))then
              niiTratio = (ILs(ion_no1)%int_dered + ILs(ion_no2)%int_dered) / ILs(ion_no3)%int_dered
      elseif(((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .eq. 0)) .and. (ILs(ion_no3)%int_dered .gt. 0 ))then
              niiTratio = (ILs(ion_no1)%int_dered * 4.054) / ILs(ion_no3)%int_dered
      elseif(((ILs(ion_no1)%int_dered .eq. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0)) .and. (ILs(ion_no3)%int_dered .gt. 0 ))then
              niiTratio = (ILs(ion_no2)%int_dered * 1.3274) / ILs(ion_no3)%int_dered
      else
              niiTratio = 0.0
      endif

! O III

      ion_no1 = get_ion("oiii5007   ", ILs)
      ion_no2 = get_ion("oiii4959   ", ILs)
      ion_no3 = get_ion("oiii4363   ", ILs)

      if(((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0)) .and. (ILs(ion_no3)%int_dered .gt. 0 ))then
              oiiiTratio = (ILs(ion_no1)%int_dered + ILs(ion_no2)%int_dered) / ILs(ion_no3)%int_dered
      elseif(((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .eq. 0)) .and. (ILs(ion_no3)%int_dered .gt. 0 ))then 
              oiiiTratio = (ILs(ion_no1)%int_dered * 1.3356) / ILs(ion_no3)%int_dered
      elseif(((ILs(ion_no1)%int_dered .eq. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0)) .and. (ILs(ion_no3)%int_dered .gt. 0 ))then 
              oiiiTratio = (ILs(ion_no2)%int_dered * 3.98) / ILs(ion_no3)%int_dered
      else
              oiiiTratio = 0.0
      endif

! Ne III

      ion_no1 = get_ion("neiii3868  ", ILs)
      ion_no2 = get_ion("neiii3967  ", ILs)
      ion_no3 = get_ion("neiii3342  ", ILs)

      if(((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0)) .and. (ILs(ion_no3)%int_dered .gt. 0 ))then
              neiiiTratio = (ILs(ion_no1)%int_dered + ILs(ion_no2)%int_dered) / ILs(ion_no3)%int_dered
      elseif(((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .eq. 0)) .and. (ILs(ion_no3)%int_dered .gt. 0 ))then
              neiiiTratio = (ILs(ion_no1)%int_dered * 1.3013) / ILs(ion_no3)%int_dered
      elseif(((ILs(ion_no1)%int_dered .eq. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0)) .and. (ILs(ion_no3)%int_dered .gt. 0 ))then
              neiiiTratio = (ILs(ion_no2)%int_dered * 4.319) / ILs(ion_no3)%int_dered
      else
              neiiiTratio = 0.0
      endif

! Ne III ir

      ion_no1 = get_ion("neiii3868  ", ILs)
      ion_no2 = get_ion("neiii3967  ", ILs)
      ion_no3 = get_ion("neiii15p5um", ILs)

      if(((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0)) .and. (ILs(ion_no3)%int_dered .gt. 0 ))then
               neiiiIRTratio = (ILs(ion_no1)%int_dered + ILs(ion_no2)%int_dered) / ILs(ion_no3)%int_dered
      elseif(((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .eq. 0)) .and. (ILs(ion_no3)%int_dered .gt. 0 ))then
               neiiiIRTratio = (ILs(ion_no1)%int_dered * 1.3013) / ILs(ion_no3)%int_dered
      elseif(((ILs(ion_no1)%int_dered .eq. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0)) .and. (ILs(ion_no3)%int_dered .gt. 0 ))then
               neiiiIRTratio = (ILs(ion_no2)%int_dered * 4.319) / ILs(ion_no3)%int_dered
      else
              neiiiIRTratio = 0.0
      endif

!!ne v

      ion_no1 = get_ion("nev3426    ", ILs)
      ion_no2 = get_ion("nev3345    ", ILs)
      ion_no3 = get_ion("nev2975    ", ILs)

      if(((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0)) .and. (ILs(ion_no3)%int_dered .gt. 0 ))then
              nevTratio = (ILs(ion_no1)%int_dered + ILs(ion_no2)%int_dered) / ILs(ion_no3)%int_dered
      elseif(((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .eq. 0)) .and. (ILs(ion_no3)%int_dered .gt. 0 ))then
              nevTratio = (ILs(ion_no1)%int_dered * 1.3571) / ILs(ion_no3)%int_dered
      elseif(((ILs(ion_no1)%int_dered .eq. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0)) .and. (ILs(ion_no3)%int_dered .gt. 0 ))then
              nevTratio = (ILs(ion_no2)%int_dered * 3.800) / ILs(ion_no3)%int_dered
      else
              nevTratio = 0.0
      endif

!s iii

      ion_no1 = get_ion("siii9069   ", ILs)
      ion_no2 = get_ion("siii9531   ", ILs)
      ion_no3 = get_ion("siii6312   ", ILs)

      if(((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0)) .and. (ILs(ion_no3)%int_dered .gt. 0 ))then
              siiiTratio = (ILs(ion_no1)%int_dered + ILs(ion_no2)%int_dered) / ILs(ion_no3)%int_dered
      elseif(((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .eq. 0)) .and. (ILs(ion_no3)%int_dered .gt. 0 ))then
              siiiTratio = (ILs(ion_no1)%int_dered * 6.685) / ILs(ion_no3)%int_dered
      elseif(((ILs(ion_no1)%int_dered .eq. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0)) .and. (ILs(ion_no3)%int_dered .gt. 0 ))then
              siiiTratio = (ILs(ion_no2)%int_dered * 1.176) / ILs(ion_no3)%int_dered
      else
              siiiTratio = 0.0
      endif

!ar iii

      ion_no1 = get_ion("ariii7135  ", ILs)
      ion_no2 = get_ion("ariii7751  ", ILs)
      ion_no3 = get_ion("ariii5192  ", ILs)

      if(((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0)) .and. (ILs(ion_no3)%int_dered .gt. 0 ))then
              ariiiTratio = (ILs(ion_no1)%int_dered + ILs(ion_no2)%int_dered) / ILs(ion_no3)%int_dered
      elseif(((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .eq. 0)) .and. (ILs(ion_no3)%int_dered .gt. 0 ))then
              ariiiTratio = (ILs(ion_no1)%int_dered * 1.240) / ILs(ion_no3)%int_dered
      elseif(((ILs(ion_no1)%int_dered .eq. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0)) .and. (ILs(ion_no3)%int_dered .gt. 0 ))then
              ariiiTratio = (ILs(ion_no2)%int_dered * 5.174) / ILs(ion_no3)%int_dered
      else
              ariiiTratio = 0.0
      endif

! ar v

      ion_no1 = get_ion("arv6435    ", ILs)
      ion_no2 = get_ion("arv7005    ", ILs)
      ion_no3 = get_ion("arv4625    ", ILs)

      if(((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0)) .and. (ILs(ion_no3)%int_dered .gt. 0 ))then
              arvTratio = (ILs(ion_no1)%int_dered + ILs(ion_no2)%int_dered) / ILs(ion_no3)%int_dered
      elseif(((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .eq. 0)) .and. (ILs(ion_no3)%int_dered .gt. 0 ))then
              arvTratio = (ILs(ion_no1)%int_dered * 3.125) / ILs(ion_no3)%int_dered
      elseif(((ILs(ion_no1)%int_dered .eq. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0)) .and. (ILs(ion_no3)%int_dered .gt. 0 ))then
              arvTratio = (ILs(ion_no2)%int_dered * 1.471) / ILs(ion_no3)%int_dered
      else
              arvTratio = 0.0
      endif

! O II

      ion_no1 = get_ion("oii7319    ",ILs)
      ion_no2 = get_ion("oii7330    ",ILs)
      ion_no3 = get_ion("oii3726    ",ILs)
      ion_no4 = get_ion("oii3729    ",ILs)

      if (ion_no1 .gt. 0 .and. ion_no2 .gt. 0 .and. ion_no3 .gt. 0 .and. ion_no4 .gt. 0) then
           oiiTratio = (ILs(ion_no1)%int_dered+ILs(ion_no2)%int_dered)/(ILs(ion_no3)%int_dered+ILs(ion_no4)%int_dered) 
       else
           oiiTratio = 0.0
       endif

! S II

      ion_no1 = get_ion("sii6716    ",ILs)
      ion_no2 = get_ion("sii6731    ",ILs)
      ion_no3 = get_ion("sii4068    ",ILs)
      ion_no4 = get_ion("sii4076    ",ILs)

      if (ion_no1 .gt. 0 .and. ion_no2 .gt. 0 .and. ion_no3 .gt. 0 .and. ion_no4 .gt. 0) then
           siiTratio = (ILs(ion_no1)%int_dered+ILs(ion_no2)%int_dered)/(ILs(ion_no3)%int_dered+ILs(ion_no4)%int_dered) 
      else
           siiTratio = 0.0
      endif

! now get diagnostics zone by zone.

! low ionisation

      lowtemp = 10000.0

      do i = 1,2

         count = 0
         if (oiiNratio .gt. 0) then 
           call get_diagnostic("oii       ","1,2/                ","1,3/                ",oiiNratio,"D",lowtemp, oiiDens) 
           count = count + 1
         endif
         if (siiNratio .gt. 0) then
           call get_diagnostic("sii       ","1,2/                ","1,3/                ",siiNratio,"D",lowtemp, siiDens)
           count = count + 1
         endif

         if (count .eq. 0) then
           lowdens = 1000.0
         else
           lowdens = (oiiDens + siiDens) / count
         endif

         count = 0

         if (oiiTratio .gt. 0) then
           call get_diagnostic("oii       ","2,4,2,5,3,4,3,5/    ","1,2,1,3/            ",oiiTratio,"T",lowdens,oiiTemp)
           count = count + 1
         else
           oiiTemp = 0.0
         endif
         if (siiTratio .gt. 0) then
           call get_diagnostic("sii       ","1,2,1,3/            ","1,4,1,5/            ",siiTratio,"T",lowdens,siiTemp)
           count = count + 1
         else
           siiTemp = 0.0
         endif
         if (niiTratio .gt. 0) then
           call get_diagnostic("nii       ","2,4,3,4/            ","4,5/                ",niiTratio,"T",lowdens,niitemp)
           count = count + 5
         else
           niitemp = 0.0
         endif

         if (count .gt. 0) then 
           lowtemp = ((5*niitemp) + siitemp + oiitemp) / count
         else
           lowtemp = 10000.0
         endif

      enddo

! medium ionisation

      medtemp = 10000.0

      do i = 1,2

         count = 0
         if (cliiiNratio .gt. 0) then
           call get_diagnostic("cliii     ","1,2/                ","1,3/                ",cliiiNratio,"D",medtemp, cliiiDens)
           count = count + 1
         endif
         if (ciiiNratio .gt. 0) then
           call get_diagnostic("ciii      ","1,2/                ","1,3/                ",ciiiNratio,"D",medtemp, ciiiDens)
           count = count + 1
         endif
         if (arivNratio .gt. 0) then
           call get_diagnostic("ariv      ","1,2/                ","1,3/                ",arivNratio,"D",medtemp, arivDens)
           count = count + 1
         endif

         if (count .eq. 0) then
           meddens = 1000.0
         else
           meddens = (ciiiDens + cliiiDens + arivDens) / count
         endif

         count = 0

         if (oiiiTratio .gt. 0) then
           call get_diagnostic("oiii      ","2,4,3,4/            ","4,5/                ",oiiiTratio,"T",meddens,oiiiTemp)
           count = count + 4
         else
           oiiiTemp = 0.0
         endif
         if (siiiTratio .gt. 0) then
           call get_diagnostic("siii      ","2,4,3,4/            ","4,5/                ",siiiTratio,"T",meddens,siiiTemp)
           count = count + 1
         else
           siiiTemp = 0.0
         endif
         if (ariiiTratio .gt. 0) then
           call get_diagnostic("ariii     ","1,4,2,4/            ","4,5/                ",ariiiTratio,"T",meddens,ariiitemp)
           count = count + 2
         else
           ariiitemp = 0.0
         endif
         if (neiiiTratio .gt. 0) then
           call get_diagnostic("neiii     ","1,4,2,4/            ","4,5/                ",neiiiTratio,"T",meddens,neiiitemp)
           count = count + 2
         else
           neiiitemp = 0.0
         endif

         if (count .gt. 0) then
           medtemp = (4*oiiitemp + siiitemp + 2*ariiitemp + 2*neiiitemp) / count
         else
           medtemp = 10000.0
         endif

      enddo

! high ionisation

      hightemp = medtemp

      do i = 1,2

         if (neivNratio .gt. 0) then
           call get_diagnostic("neiv      ","1,2/                ","1,3/                ",neivNratio,"D",hightemp, neivDens)
           highdens = neivdens
         else
           neivDens = 0.0
           highdens = 1000.0
         endif

         count = 0

         if (arvTratio .gt. 0) then
           call get_diagnostic("arv       ","2,4,3,4/            ","4,5/                ",arvTratio,"T",highdens,arvTemp)
           count = count + 1
         else
           arvTemp = 0.0
         endif
         if (nevTratio .gt. 0) then
           call get_diagnostic("nev       ","2,4,3,4/            ","4,5/                ",nevTratio,"T",highdens,nevtemp)
           count = count + 1
         else
           nevtemp = 0.0
         endif

         if (count .gt. 0) then
           hightemp = (arvtemp + nevtemp) / count
         else
           hightemp = medtemp
         endif

      enddo

!done calculating, now write out.

      print *,""
      print *,"Diagnostics"
      print *,"==========="
      print *,""

      print *,"Diagnostic       Zone      Value"
      print *,""
      print "(A28,F8.0)","[O II] density   Low       ",oiidens
      print "(A28,F8.0)","[S II] density   Low       ",siidens
      print "(A28,F8.0)","   adopted       Low       ",lowdens
      print *,""

      print "(A28,F8.0)","[N II] temp      Low       ",niitemp
      print "(A28,F8.0)","[O II] temp      Low       ",oiitemp
      print "(A28,F8.0)","[S II] temp      Low       ",siitemp
      print "(A28,F8.0)","   adopted       Low       ",lowtemp
      print *,""

      print "(A28,F8.0)","[Cl III] density Medium    ",cliiidens
      print "(A28,F8.0)","[Ar IV] density  Medium    ",arivdens
      print "(A28,F8.0)","C III] density   Medium    ",ciiidens
      print "(A28,F8.0)","   adopted       Medium    ",meddens
      print *,""

      print "(A28,F8.0)","[O III] temp     Medium    ",oiiitemp
      print "(A28,F8.0)","[Ne III] temp    Medium    ",neiiitemp
      print "(A28,F8.0)","[Ar III] temp    Medium    ",ariiitemp
      print "(A28,F8.0)","[S III] temp     Medium    ",siiitemp
      print "(A28,F8.0)","   adopted       Medium    ",medtemp
      print *,""

      print "(A28,F8.0)","[Ne IV] density  High      ",neivdens
      print "(A28,F8.0)","   adopted       High      ",highdens
      print *,""
      print "(A28,F8.0)","[Ar V] temp      High      ",arvtemp
      print "(A28,F8.0)","[Ne V] temp      High      ",nevtemp
      print "(A28,F8.0)","   adopted       High      ",hightemp

! later, make this check with user whether to adopt these values

!      print *,"Enter an option:"
!      print *,"1. Use these diagnostics"
!      print *,"2. Input your own"
!
!      read (5,*) opt
!      if (opt .ne. 1) then
!        print *,"Input low-ionisation zone density: "
!        read (5,*) lowdens
!        print *,"Input low-ionisation zone temperature: "
!        read (5,*) lowtemp
!        print *,"Input medium-ionisation zone density: "
!        read (5,*) meddens
!        print *,"Input medium-ionisation zone temperature: "
!        read (5,*) medtemp
!        print *,"Input high-ionisation zone density: "
!        read (5,*) highdens
!        print *,"Input high-ionisation zone temperature: "
!        read (5,*) hightemp
!      endif

! Helium abundances

        print *,""
        print *,"Abundances"
        print *,"----------"

        print *,"Helium"
        print *,"------"
        call get_helium(REAL(oiiiTemp),REAL(oiiDens),REAL(He_lines(1)%int_dered),REAL(He_lines(2)%int_dered),REAL(He_lines(3)%int_dered),REAL(He_lines(4)%int_dered),heiabund,heiiabund,Hetotabund)
        print "(1X,A10,F5.3)", "He+/H+ =  ",heiabund
        print "(1X,A10,F5.3)", "He++/H+ = ",heiiabund
        print "(1X,A10,F5.3)", "He/H =    ",Hetotabund

! get abundances for all CELs

        print *,""
        print *,"CELs"
        print *,"----"
        print *,"Ion         I(lambda)   Abundance"

        do i = 1,Iint-1 !XXXX why does Int end up 1 too high?
!                 print *,ILs(i)%ion,ILs(i)%transition,ILs(i)%int_dered
           if (ILs(i)%zone .eq. "low ") then
                 call get_abundance(ILs(i)%ion, ILs(i)%transition, lowtemp, lowdens,ILs(i)%int_dered, ILs(i)%abundance)
           elseif (ILs(i)%zone .eq. "med ") then
                 call get_abundance(ILs(i)%ion, ILs(i)%transition, medtemp, meddens,ILs(i)%int_dered, ILs(i)%abundance)
           elseif (ILs(i)%zone .eq. "high") then
                 call get_abundance(ILs(i)%ion, ILs(i)%transition, hightemp, highdens,ILs(i)%int_dered, ILs(i)%abundance)
           endif
                 if (ILs(i)%abundance .gt. 0) then
                       PRINT "(1X, A11, 1X, F7.3, 5X, ES9.3)",ILs(i)%name,ILs(i)%int_dered,ILs(i)%abundance
                 endif
        enddo

! calculate averages

        ciiiCELabund = ILs(1)%abundance

        celabundtemp = 0.
        weight = 0.
        do i=2,3
          niiCELabund = niiCELabund + ILs(i)%abundance
          if (ILs(i)%abundance .gt. 0) then
            weight = weight + 1
          endif
        enddo
        if (weight .gt. 0) then
          niiCELabund = niiCELabund / weight
        else
          niiCELabund = 0.0
        endif

        niiiIRCELabund = ILs(5)%abundance
        niiiUVCELabund = ILs(6)%abundance

        if (niiiIRCELabund .gt. 0 .and. niiiUVCELabund .gt. 0) then
          niiiCELabund = (niiiIRCELabund + niiiUVCELabund)/2
        elseif (niiiIRCELabund .gt. 0) then
          niiiCELabund = niiiIRCELabund
        elseif (niiiUVCELabund .gt. 0) then
          niiiCELabund = niiiUVCELabund
        else
          niiiCELabund = 0
        endif

        celabundtemp = 0.
        weight = 0.
        do i=7,10
          oiiCELabund = oiiCELabund + ILs(i)%abundance
          if (ILs(i)%abundance .gt. 0) then
            weight = weight + 1
          endif
        enddo
        if (weight .gt. 0) then
          oiiCELabund = oiiCELabund / weight
        else
          oiiCELabund = 0.0
        endif

        celabundtemp = 0.
        weight = 0.
        do i=12,13
          oiiiCELabund = oiiiCELabund + ILs(i)%abundance
          if (ILs(i)%abundance .gt. 0) then
            weight = weight + 1
          endif
        enddo
        if (weight .gt. 0) then
          oiiiCELabund = oiiiCELabund / weight
        else
          oiiiCELabund = 0.0
        endif

        celabundtemp = 0.
        weight = 0.
        do i=14,15
          oiiiIRCELabund = oiiiIRCELabund + ILs(i)%abundance
          if (ILs(i)%abundance .gt. 0) then
            weight = weight + 1
          endif
        enddo
        if (weight .gt. 0) then
          oiiiIRCELabund = oiiiIRCELabund / weight
        else
          oiiiIRCELabund = 0.0
        endif

        neiiIRCELabund = ILs(16)%abundance
        neiiiIRCELabund = ILs(17)%abundance

        celabundtemp = 0.
        weight = 0.
        do i=18,19
          neiiiCELabund = neiiiCELabund + ILs(i)%abundance
          if (ILs(i)%abundance .gt. 0) then
            weight = weight + 1
          endif
        enddo
        if (weight .gt. 0) then
          neiiiCELabund = neiiiCELabund / weight
        else
          neiiiCELabund = 0.0
        endif


        celabundtemp = 0.
        weight = 0.
        do i=20,23
          neivCELabund = neivCELabund + ILs(i)%abundance
          if (ILs(i)%abundance .gt. 0) then
            weight = weight + 1
          endif
        enddo
        if (weight .gt. 0) then
          neivCELabund = neivCELabund / weight
        else
          neivCELabund = 0.0
        endif

        celabundtemp = 0.
        weight = 0.
        do i=24,27
          siiCELabund = siiCELabund + ILs(i)%abundance
          if (ILs(i)%abundance .gt. 0) then
            weight = weight + 1
          endif
        enddo
        if (weight .gt. 0) then
          siiCELabund = siiCELabund / weight
        else
          siiCELabund = 0.0
        endif

        siiiCELabund = ILs(28)%abundance
        siiiIRCELabund = ILs(29)%abundance
        sivIRCELabund = ILs(30)%abundance

        celabundtemp = 0.
        weight = 0.
        do i=31,32
          cliiiCELabund = cliiiCELabund + ILs(i)%abundance
          if (ILs(i)%abundance .gt. 0) then
            weight = weight + 1
          endif
        enddo
        if (weight .gt. 0) then
          cliiiCELabund = cliiiCELabund / weight
        else
          cliiiCELabund = 0.0
        endif

        celabundtemp = 0.
        weight = 0.
        do i=33,34
          ariiiCELabund = ariiiCELabund + ILs(i)%abundance
          if (ILs(i)%abundance .gt. 0) then
            weight = weight + 1
          endif
        enddo
        if (weight .gt. 0) then
          ariiiCELabund = ariiiCELabund / weight
        else
          ariiiCELabund = 0.0
        endif

        celabundtemp = 0.
        weight = 0.
        do i=35,36
          arivCELabund = arivCELabund + ILs(i)%abundance
          if (ILs(i)%abundance .gt. 0) then
            weight = weight + 1
          endif
        enddo
        if (weight .gt. 0) then
          arivCELabund = arivCELabund / weight
        else
          arivCELabund = 0.0
        endif

        ariiiIRCELabund = ILs(37)%abundance

        celabundtemp = 0.
        weight = 0.
        do i=38,39 
          if (ILs(i)%abundance .gt. 0) then
            ciiiCELabund = ciiiCELabund + ILs(i)%abundance
            weight = weight + 1
          endif
        enddo
        if (weight .gt. 0) then
          ciiiCELabund = ciiiCELabund / weight
        else
          ciiiCELabund = 0.0
        endif

        celabundtemp = 0.
        weight = 0.
        do i=40,41
          neivCELabund = neivCELabund + ILs(i)%abundance
          if (ILs(i)%abundance .gt. 0) then
            weight = weight + 1
          endif
        enddo
        if (weight .gt. 0) then
          neivCELabund = neivCELabund / weight
        else
          neivCELabund = 0.0
        endif

        celabundtemp = 0.
        weight = 0.
        do i=44,45
          nevCELabund = nevCELabund + ILs(i)%abundance
          if (ILs(i)%abundance .gt. 0) then
            weight = weight + 1
          endif
        enddo
        if (weight .gt. 0) then
          nevCELabund = nevCELabund / weight
        else
          nevCELabund = 0.0
        endif

        celabundtemp = 0.
        weight = 0.
        do i=46,47
          siiiCELabund = siiiCELabund + ILs(i)%abundance
          if (ILs(i)%abundance .gt. 0) then
            weight = weight + 1
          endif
        enddo
        if (weight .gt. 0) then
          siiiCELabund = siiiCELabund / weight
        else
          siiiCELabund = 0.0
        endif

        celabundtemp = 0.
        weight = 0.
        do i=50,51
          arvCELabund = arvCELabund + ILs(i)%abundance
          if (ILs(i)%abundance .gt. 0) then
            weight = weight + 1
          endif
        enddo
        if (weight .gt. 0) then
          arvCELabund = arvCELabund / weight
        else
          arvCELabund = 0.0
        endif

! now get abundances for ORLs
! o2+
       call oii_rec_lines(oiiiTemp,oiiDens,DBLE(1),oiiRLs)

       do i = 1,329 
         do j = 1,415 
          if (abs(O_dered(1,i)-oiiRLs(j)%Wave) .le. 0.005) then
            oiiRLs(j)%Obs = O_dered(2,i) 
            oiiRLs(j)%abundance = oiiRLs(j)%obs/oiiRLs(j)%Int 
          endif
         enddo
       enddo

!N2+

       call nii_rec_lines(oiiiTemp,oiiDens,DBLE(1),niiRLs)

       do i = 1,329
         do j = 1,99
          if (abs(O_dered(1,i)-niiRLs(j)%Wave) .le. 0.005) then
            niiRLs(j)%Obs = O_dered(2,i)
            niiRLs(j)%abundance = niiRLs(j)%obs/niiRLs(j)%Int
          endif
         enddo
       enddo

!C2+ 
       call cii_rec_lines(oiiiTemp,oiiDens,DBLE(1),ciiRLs)

       do i = 1,329
         do j = 1,57
          if (abs(O_dered(1,i)-ciiRLs(j)%Wave) .le. 0.005) then
            ciiRLs(j)%Obs = O_dered(2,i)
            ciiRLs(j)%abundance = ciiRLs(j)%obs/ciiRLs(j)%Int
          endif
         enddo
       enddo

!Ne2+
       call neii_rec_lines(oiiiTemp,oiiDens,DBLE(1),neiiRLs)

       do i = 1,329
         do j = 1,38
          if (abs(O_dered(1,i)-neiiRLs(j)%Wave) .le. 0.005) then
            neiiRLs(j)%Obs = O_dered(2,i)
            neiiRLs(j)%abundance = neiiRLs(j)%obs/neiiRLs(j)%Int
          endif
         enddo
       enddo

!C3+, N3+
       call xiii_rec_lines(oiiiTemp,oiiDens,DBLE(1),xiiiRLs)

       do i = 1,329
         do j = 1,6
          if (abs(O_dered(1,i)-xiiiRLs(j)%Wave) .le. 0.005) then
            xiiiRLs(j)%Obs = O_dered(2,i)
            xiiiRLs(j)%abundance = xiiiRLs(j)%obs/xiiiRLs(j)%Int
          endif
         enddo
       enddo

      print *,""
      print *,"Recombination lines"
      print *,"-------------------"

      rlabundtemp = 0.0
      weight = 0.0

!cii recombination lines

      print *,""
      print *,"CII"
      print *,"lambda   Int   Abund"
      do i = 1,57
        if (ciiRLs(i)%abundance .gt. 0) then 
          print "(1X,F7.2,1X,F6.3,1X,ES9.3)",ciiRLs(i)%wave,ciiRLs(i)%obs,ciiRLs(i)%abundance
          rlabundtemp = rlabundtemp + ciiRLs(i)%obs
          weight = weight + ciiRLs(i)%Int
        endif 
        if (ciiRLs(i)%wave .eq. 4267.15D0) then
          ciirlabund = ciiRLs(i)%abundance
        endif
      enddo

      print "(A34,ES9.3)","Abundance (all lines co-added): ",rlabundtemp/weight
      if (ciirlabund .gt. 0) then
        print "(A34,ES9.3)","Abundance (4267 line only):     ",ciirlabund
      else
        ciirlabund = rlabundtemp/weight
      endif

!nii recombination lines

      print *,""
      print *,"NII"
!      print *,"lambda   Mult   Int   Abund"
      do i = 1,99
        if (niiRLs(i)%abundance .gt. 0) then
!          print "(F7.2,1X,A7,1X,F6.3,1X,ES9.3)",niiRLs(i)%wave,niiRLs(i)%Mult,niiRLs(i)%obs,niiRLs(i)%abundance
          rlabundtemp = rlabundtemp + niiRLs(i)%obs
          weight = weight + niiRLs(i)%Int
        endif
      enddo

      print *,"Abundance from co-added intensity: "
      print "(ES9.2)",rlabundtemp / weight

      niimultiplets%Multiplet = (/"V3     ","V5     ","V8     " ,"V12    ","V20    ","V28    ","3d-4f  ","","","","","","","","","","","","",""/)

! get multiplet abundances from coadded intensity

      print *,"Mult    Intensity   N2+/H+"

      do j = 1,6
        rlabundtemp = 0.
        weight = 1.
        do i = 1,99
          if (niiRLs(i)%Mult .eq. niimultiplets(j)%Multiplet .and. niiRLs(i)%obs .gt. 0) then
!            rlabundtemp = rlabundtemp + (niiRLs(i)%obs * niiRLs(i)%abundance)
!            weight = weight + niiRLs(i)%obs
             rlabundtemp = rlabundtemp + niiRLs(i)%obs
             weight = weight + niiRLs(i)%Int 
          endif
        enddo 
        print "(1X,A7,F6.3,7X,ES9.3)",niimultiplets(j)%Multiplet,rlabundtemp, rlabundtemp/weight
        niimultiplets(j)%Abundance = rlabundtemp/weight
      enddo

      rlabundtemp = 0.
      weight = 0.
      do i = 77,99
        if (niiRLs(i)%obs .gt. 0) then
!          rlabundtemp = rlabundtemp + (niiRLs(i)%obs * niiRLs(i)%abundance)
!          weight = weight + niiRLs(i)%obs
           rlabundtemp = rlabundtemp + niiRLs(i)%obs
           weight = weight + niiRLs(i)%Int
        endif
      enddo

      niimultiplets(7)%abundance = rlabundtemp/weight

      print "(1X,A7,F6.3,7X,ES9.3)",niimultiplets(7)%Multiplet,rlabundtemp, niimultiplets(7)%abundance
!      print "(F6.3,16X,ES9.3)",rlabundtemp, rlabundtemp/weight

      rlabundtemp = 0.0
      weight = 0
      do i = 1,7
        rlabundtemp = rlabundtemp + niimultiplets(i)%abundance
        if (niimultiplets(i)%abundance .gt. 0) then
          weight = weight + 1
        endif
      enddo

      print *,"Abundance - mean of each multiplet's abundance:"
      niiRLabund = rlabundtemp/weight
      print "(ES9.3)",niiRLabund

!oii recombination lines

      rlabundtemp = 0.00
      weight = 0.00

      print *,""
      print *,"OII"
!      print *,"lambda   Mult   Int   Abund"
      do i = 1,415
        if (oiiRLs(i)%abundance .gt. 0) then
!          print "(F7.2,1X,A7,1X,F6.3,1X,ES9.3)",oiiRLs(i)%wave,oiiRLs(i)%Mult,oiiRLs(i)%obs,oiiRLs(i)%abundance
          rlabundtemp = rlabundtemp + oiiRLs(i)%obs
          weight = weight + oiiRLs(i)%Int
        endif
      enddo

      print *,"Abundance from co-added intensity: "
      print "(ES9.2)",rlabundtemp/weight

      oiimultiplets%Multiplet = (/" V1    "," V2    "," V5    " ," V10   "," V11   "," V12   "," V19   "," V20   "," V25   "," V28   "," V33   "," 3d-4f ","","","","","","","",""/)

! get multiplet abundances from coadded intensity

      print *,"Co-added intensity   O2+/H+"

      do j = 1,11
        rlabundtemp = 0.
        weight = 0.
        do i = 1,415
          if (oiiRLs(i)%Mult .eq. oiimultiplets(j)%Multiplet .and. oiiRLs(i)%obs .gt. 0) then
!            rlabundtemp = rlabundtemp + (oiiRLs(i)%obs * oiiRLs(i)%abundance)
!            weight = weight + oiiRLs(i)%obs
             rlabundtemp = rlabundtemp + oiiRLs(i)%obs
             weight = weight + oiiRLs(i)%Int
          endif
        enddo 
        if (weight .gt. 0) then
          oiimultiplets(j)%Abundance = rlabundtemp/weight 
        else
          oiimultiplets(j)%Abundance = 0.0
        endif
        print "(1X,A7,F6.3,7X,ES9.3)",oiimultiplets(j)%Multiplet,rlabundtemp,oiimultiplets(j)%abundance
      enddo

      rlabundtemp = 0.
      weight = 0.
      do i = 1,182
        if (oiiRLs(i)%Mult .ne. "       " .and. oiiRLs(i)%obs .gt. 0) then
!          rlabundtemp = rlabundtemp + (oiiRLs(i)%obs * oiiRLs(i)%abundance)
!          weight = weight + oiiRLs(i)%obs
           rlabundtemp = rlabundtemp + oiiRLs(i)%obs
           weight = weight + oiiRLs(i)%Int
        endif
      enddo

      oiimultiplets(12)%abundance = rlabundtemp/weight

      print "(1X,A7,F6.3,7X,ES9.3)",oiimultiplets(j)%Multiplet,rlabundtemp, rlabundtemp/weight
!      print *,"3d-4f :"
!      print *,"Co-added intensity   O2+/H+"
!      print "(F6.3,16X,ES9.3)",rlabundtemp, rlabundtemp/weight

      rlabundtemp = 0.0
      weight = 0
      do i = 1,7
        rlabundtemp = rlabundtemp + oiimultiplets(i)%abundance
        if (oiimultiplets(i)%abundance .gt. 0) then
          weight = weight + 1
        endif
      enddo

      print *,"Abundance - mean of each multiplet's abundance:"
      oiiRLabund = rlabundtemp/weight
      print "(ES9.3)",oiiRLabund

!neii recombination lines

      rlabundtemp = 0.0
      weight = 0.0

      print *,""
      print *,"NeII"
!      print *,"lambda   Mult   Int   Abund"
      do i = 1,38
        if (neiiRLs(i)%abundance .gt. 0) then
!          print "(F7.2,1X,F6.3,1X,ES9.3)",neiiRLs(i)%wave,neiiRLs(i)%obs,neiiRLs(i)%abundance 
           rlabundtemp = rlabundtemp + neiiRLs(i)%obs
           weight = weight + neiiRLs(i)%Int 
        endif
      enddo
      neiiRLabund = rlabundtemp/weight

      print *,"Abundance from co-added intensities: "
      print "(ES9.3)",neiiRLabund

      rlabundtemp = 0.0
      weight = 0.0
      print *,""
      print *,"CIII"
      print *,"lambda   Mult   Int   Abund"
      do i = 1,4
        if (xiiiRLs(i)%abundance .gt. 0) then
          print "(F7.2,1X,F6.3,1X,ES9.3)",xiiiRLs(i)%wave,xiiiRLs(i)%obs,xiiiRLs(i)%abundance
          rlabundtemp = rlabundtemp + xiiiRLs(i)%obs
          weight = weight + xiiiRLs(i)%Int
        endif
      enddo
      if (weight .gt. 0) then
        ciiiRLabund = rlabundtemp / weight
      else
        ciiiRLabund = 0.0
      endif

      print *,""
      print *,"NIII"
      print *,"lambda   Mult   Int   Abund"
      do i = 5,6
        if (xiiiRLs(i)%abundance .gt. 0) then
           print "(F7.2,1X,F6.3,1X,ES9.3)",xiiiRLs(i)%wave,xiiiRLs(i)%obs,xiiiRLs(i)%abundance
        endif
      enddo

     if (xiiiRLs(6)%abundance .gt. 0) then
        niiiRLabund = xiiiRLs(6)%abundance
     else
        niiiRLabund = 0.0
     endif

! ICFs (Kingsburgh + Barlow 1994)

! oxygen, just the simple case at the moment

     CELicfO = ((heiabund + heiiabund)/heiabund)**(2/3)
     OabundCEL = CELicfO * (oiiCELabund + oiiiCELabund)

! nitrogen

     if (niiCELabund .gt. 0 .and. niiiUVCELabund .gt. 0 .and. nivCELabund .gt. 0) then
       CELicfN = 1.
       NabundCEL = niCELabund + niiCELabund + niiiCELabund
     elseif (niiCELabund .gt. 0 .and. niiiUVCELabund .eq. 0 .and. nivCELabund .gt. 0) then
       CELicfN = 1.5
       NabundCEL = 1.5*(niCELabund + niiiUVCELabund)
     elseif (niiCELabund .gt. 0 .and. niiiUVCELabund .eq. 0 .and. nivCELabund .eq.0) then
       CELicfN = OabundCEL/oiiCELabund
       NabundCEL = niiCELabund * CELicfN
     endif

! carbon - couple of cases still to add.

     if (ciiCELabund .gt. 0 .and. ciiiCELabund .gt. 0 .and. civCELabund .gt. 0 .and. heiiabund .eq. 0) then 
       CELicfC = 1.
       CabundCEL = ciiCELabund + ciiiCELabund + civCELabund
     elseif (ciiCELabund .eq. 0 .and. ciiiCELabund .gt. 0 .and. civCELabund .gt. 0 .and. heiiabund .eq. 0) then
       CELicfC = (oiiCELabund + oiiiCELabund) / oiiiCELabund
       CabundCEL = CELicfC * (ciiCELabund + ciiiCELabund)
     elseif (ciiCELabund .eq. 0 .and. ciiiCELabund .gt. 0 .and. civCELabund .eq. 0 .and. oiiiCELabund .gt. 0) then
       CELicfC = OabundCEL/oiiiCELabund
       CabundCEL = CELicfC * ciiiCELabund
     elseif (nvCELabund .gt. 0 .and. heiiabund .gt. 0) then
       CELicfC = 1/(1-(nvCELabund/(niiCELabund + niiiCELabund + nivCELabund + nvCELabund)))
       if (CELicfC .gt. 5) then
         CELicfC = (niiCELabund + niiiCELabund + nivCELabund + nvCELabund) / (niiCELabund + niiiCELabund + nivCELabund)
       endif
       CabundCEL = CELicfC * (ciiCELabund + ciiiCELabund + civCELabund)
     endif

! Neon

     if (neiiiCELabund .gt. 0 .and. neivCELabund .gt. 0 .and. nevCELabund .gt. 0) then
       CELicfNe = 1.
       NeabundCEL = neiiiCELabund + neivCELabund + nevCELabund
     elseif (neiiiCELabund .gt. 0 .and. neivCELabund .eq. 0 .and. nevCELabund .gt. 0) then
       CELicfNe = 1.5
       NeabundCEL = CELicfNe * (neiiiCELabund + nevCELabund)
     elseif (neiiiCELabund .gt. 0 .and. neivCELabund .eq. 0 .and. nevCELabund .eq. 0) then
       CELicfNe = OabundCEL / oiiiCELabund
       NeabundCEL = CELicfNe * neiiiCELabund
     endif

! Argon

     if (ariiiCELabund .gt. 0 .and. arivCELabund .eq. 0 .and. arvCELabund .eq. 0) then
       CELicfAr = 1.87
       ArabundCEL = CELicfAr * ariiiCELabund
     elseif (ariiiCELabund .eq. 0 .and. arivCELabund .gt. 0 .and. arvCELabund .eq. 0) then
       CELicfAr = NeabundCEL / neiiiCELabund
       ArabundCEL = CELicfAr * arivCELabund
     elseif (ariiiCELabund .gt. 0 .and. arivCELabund .gt. 0) then
       CELicfAr = 1.0
       ArabundCEL = ariiiCELabund + arivCELabund + arvCELabund
     endif

! Sulphur

     if (siiCELabund .gt. 0 .and. siiiCELabund .gt. 0 .and. sivIRCELabund .eq. 0) then
       CELicfS = (1-(1-((oiiCELabund/OabundCEL)**3)))**(-1/3)
       SabundCEL = CELicfS * (siiCELabund + siiiCELabund)
     elseif (siiCELabund .gt. 0 .and. siiiCELabund .gt. 0 .and. sivIRCELabund .gt. 0) then
       CELicfS = 1.
       SabundCEL = siiCELabund + siiiCELabund + sivIRCELabund
     elseif (siiCELabund .gt. 0 .and. siiiCELabund .eq. 0 .and. sivIRCELabund .eq. 0) then
       CELicfS = ((oiiiCELabund/oiiCELabund)**0.433 + 4.677) * (1-(1-((oiiCELabund/OabundCEL)**3)))**(-1/3)
       SabundCEL = CELicfS * siiCELabund
     endif

!XXXX finish very high excitation case.

!ORLs
!Oxygen

     if (oiiRLabund .gt. 0) then
       RLicfO = ((heiabund + heiiabund)/heiabund)**(2/3) * (1+(oiiCELabund/oiiiCELabund))
       OabundRL = RLicfO * oiiRLabund
     else
       RLicfO = 0
       OabundRL = 0
     endif

!Nitrogen

     RLicfN = 1.0
     NabundRL = niiRLabund + niiiRLabund

!Carbon

     RLicfC = 1.0
     CabundRL = ciiRLabund + ciiiRLabund

!Neon

     RLicfNe = OabundRL / oiiRLabund
     NeabundRL = RLicfNe * neiiRLabund

!finish these XXXX

print *,""
print *,"Total abundances"
print *,"================"
print *,""

print *,"CELs"
print *,""
print *,"Element     ICF     X/H"
print *,"-------     ---     ---"
print "(A12,F5.2,4X,ES8.2,2X,F5.2)"," Carbon      ",CELicfC,CabundCEL,12+log10(CabundCEL)
print "(A12,F5.2,4X,ES8.2,2X,F5.2)"," Nitrogen    ",CELicfN,NabundCEL,12+log10(NabundCEL)
print "(A12,F5.2,4X,ES8.2,2X,F5.2)"," Oxygen      ",CELicfO,OabundCEL,12+log10(OabundCEL)
print "(A12,F5.2,4X,ES8.2,2X,F5.2)"," Neon        ",CELicfNe,NeabundCEL,12+log10(NeabundCEL)
print "(A12,F5.2,4X,ES8.2,2X,F5.2)"," Argon       ",CELicfAr,ArabundCEL,12+log10(ArabundCEL)
print "(A12,F5.2,4X,ES8.2,2X,F5.2)"," Sulphur     ",CELicfS,SabundCEL,12+log10(SabundCEL)

print *,""
print *,"ORLs"
print *,""
print *,"Element     ICF     X/H"
print *,"-------     ---     ---"
print "(A12,F5.2,4X,ES8.2,2X,F5.2)"," Helium      ",1.0,Hetotabund,12+log10(Hetotabund)
print "(A12,F5.2,4X,ES8.2,2X,F5.2)"," Carbon      ",RLicfC,CabundRL,12+log10(CabundRL)
print "(A12,F5.2,4X,ES8.2,2X,F5.2)"," Nitrogen    ",RLicfN,NabundRL,12+log10(NabundRL)
print "(A12,F5.2,4X,ES8.2,2X,F5.2)"," Oxygen      ",RLicfO,OabundRL,12+log10(OabundRL)
print "(A12,F5.2,4X,ES8.2,2X,F5.2)"," Neon        ",RLicfNe,NeabundRL,12+log10(NeabundRL)


!abundance discrepancy factors

print *,""
print *,"Abundance Discrepancy Factors"
print *,"============================="
print *,""

  if (oiiiCELabund .gt. 0) then
    adfO2plus = oiiRLabund/oiiiCELabund
  else
    adfO2plus = 0.0
  endif

  if (oabundCEL .gt. 0) then
    adfO = OabundRL/OabundCEL
  else
    adfO = 0.0
  endif


  if (ciiiCELabund .gt. 0) then
    adfC2plus = ciiRLabund/ciiiCELabund
  else
    adfC2plus = 0.0
  endif

  if (CabundCEL .gt. 0) then
    adfC = CabundRL/CabundCEL
  else
    adfC = 0.0
  endif


  if (NiiiCELabund .gt. 0) then
    adfN2plus = NiiRLabund/NiiiCELabund
  else
    adfN2plus = 0.0
  endif

  if (NabundCEL .gt. 0) then
    adfN = NabundRL/NabundCEL
  else
    adfN = 0.0
  endif


  if (NeiiiCELabund .gt. 0) then
    adfNe2plus = NeiiRLabund/NeiiiCELabund
  else
    adfNe2plus = 0.0
  endif

  if (NeabundCEL .gt. 0) then
    adfNe = NeabundRL/NeabundCEL
  else
    adfNe = 0.0
  endif

print "(A12,F5.2)","adf (O2+) = ", adfo2plus
print "(A12,F5.2)","adf (O)   = ", adfO
print *,""
print "(A12,F5.2)","adf (N2+) = ", adfn2plus
print "(A12,F5.2)","adf (N)   = ", adfn
print *,""
print "(A12,F5.2)","adf (C2+) = ", adfc2plus
print "(A12,F5.2)","adf (C)   = ", adfc
print *,""
print "(A12,F5.2)","adf (Ne2+)= ", adfne2plus
print "(A12,F5.2)","adf (Ne)  = ", adfne

end program abundances 
