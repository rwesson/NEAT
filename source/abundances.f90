
subroutine abundances(linelist, switch_ext, listlength, iteration_result, R, meanextinction, calculate_extinction, ILs, Iint, diagnostic_array,iion,atomicdata,maxlevs,maxtemps)
use mod_abundmaths
use mod_abundtypes
use mod_equib
use mod_abundIO
use mod_helium
use mod_recombination_lines
use mod_extinction
use mod_resultarrays
use mod_atomicdata

implicit none

        INTEGER :: count, Iint, i, j, ion_no1, ion_no2, ion_no3, ion_no4, ion_no5, ion_no6
        integer, intent(in) :: listlength
        TYPE(line), dimension(listlength) :: linelist, linelist_orig
        CHARACTER :: switch_ext !switch for extinction laws
        type(resultarray), dimension(1) :: iteration_result
        double precision, dimension(6), intent(in) :: diagnostic_array

        DOUBLE PRECISION :: normalise, oiiNratio, oiiDens, oiiiTratio, oiiiTemp, oiiiIRNratio, oiiiIRTratio, oiiiIRtemp, oiiiIRdens, niiTratio, niiTemp, ariiiIRNratio, ariiiIRdens, arivNratio, arivDens, cliiiNratio, cliiiDens, siiNratio, siiDens, siiTratio, siiTemp, siiiIRNratio, siiiIRdens, oiiTratio, oiiTemp, neiiiTratio, neiiiIRTratio, neiiiIRNratio, neiiiIRdens, neiiiTemp, neiiiIRTemp, oitemp, citemp
        DOUBLE PRECISION :: ciiiNratio,neivNratio,nevTratio,siiiTratio,ariiiTratio,arvTratio,lowtemp,lowdens,medtemp,ciiidens,meddens,siiitemp,ariiitemp,hightemp,neivdens,highdens,arvtemp,nevtemp,oiTratio,ciTratio
        DOUBLE PRECISION :: oiiRLabund, niiRLabund, ciiRLabund, cii4267rlabund, neiiRLabund, ciiiRLabund, niiiRLabund, RLabundtemp, weight
        DOUBLE PRECISION :: ciiiCELabund, niiCELabund, niiiIRCELabund, niiiUVCELabund, oiiCELabund, oiiiCELabund, oiiiIRCELabund, oivCELabund, neiiIRCELabund, neiiiIRCELabund, neiiiCELabund, neivCELabund, siiCELabund, siiiCELabund, siiiIRCELabund, sivIRCELabund, cliiiCELabund, ariiiCELabund, arivCELabund, ariiiIRCELabund, nivCELabund, niiiCELabund, ciiCELabund, civCELabund, nvCELabund, nevCELabund, arvCELabund, CELabundtemp, ciCELabund, oiCELabund
        DOUBLE PRECISION :: fn4
        DOUBLE PRECISION :: CELicfO, CELicfC, CELicfN, CELicfNe, CELicfAr, CELicfS, CELicfCl
        DOUBLE PRECISION :: RLicfO, RLicfC, RLicfN, RLicfNe
        DOUBLE PRECISION :: CabundRL, CabundCEL, NabundRL, NabundCEL, OabundRL, OabundCEL, NeabundRL, NeabundCEL, SabundCEL, ArabundCEL, NOabundCEL, NCabundCEL, ClabundCEL
        DOUBLE PRECISION :: adfC, adfN, adfO, adfNe, w1, w2, w3, w4
        DOUBLE PRECISION :: adfC2plus, adfN2plus, adfO2plus, adfNe2plus
        DOUBLE PRECISION :: c1, c2, c3, meanextinction, A4471, A4686, A6678, A5876, R
        REAL :: heiabund,heiiabund,Hetotabund

        logical :: calculate_extinction

        TYPE(line), DIMENSION(Iint) :: ILs
        TYPE(line), DIMENSION(38) :: H_BS
        TYPE(line), DIMENSION(4) :: He_lines

        !atomic data
        character*10 :: ionlist(40) !list of ion names
        integer :: iion !# of ions in Ilines
        integer :: maxlevs,maxtemps
        type(atomic_data) :: atomicdata(iion)

! recombination line variables

        TYPE RLabund
           CHARACTER*7 :: Multiplet
           REAL*8 :: Abundance
        END TYPE

        TYPE (RLabund), DIMENSION(12) :: oiimultiplets
        TYPE (RLabund), DIMENSION(7) :: niimultiplets

! strong line variables
        DOUBLE PRECISION :: X23,O_R23upper, O_R23lower, N2,O_N2, O3N2, O_O3N2, Ar3O3, O_Ar3O3, S3O3, O_S3O3, x23temp1, x23temp2, x23temp3, x23temp4

        ! initialise some variables to zero

        oiiRLabund = 0.d0
        niiRLabund = 0.d0
        CELicfCl = 0.d0
        CELicfS = 0.d0
        CELicfAr = 0.d0
        CELicfNe = 0.d0
        CELicfN = 0.d0
        CELicfO = 0.d0

        linelist_orig = linelist

        H_BS%intensity = 0
        He_lines%intensity = 0

        !file reading stuff

        !assign IDs to CELs

        CALL element_assign(ILs, linelist, Iint, listlength)

        !dereddening

        ILs%abundance = 0
        ILs%int_dered = 0
        H_BS%abundance = 0
        H_BS%int_dered = 0
        He_lines%abundance = 0
        He_lines%int_dered = 0

        !first lets find some hydrogen lines
        CALL get_H(H_BS, linelist, listlength)
        call get_He(He_lines, linelist, listlength)

        !aside: normalisation check, correction

        if(H_BS(2)%intensity .ne. 100)then
                normalise =  DBLE(100) / DBLE(H_BS(2)%intensity)
                do i = 1, Iint !normalising important ions
                        ILs(i)%intensity = ILs(i)%intensity * normalise
                end do
                do i = 1, 38 !normalising balmer series
                        H_BS(i)%intensity = H_BS(i)%intensity*normalise
                end do
                do i = 1,4 !normalise helium
                        He_lines(i)%intensity = He_lines(i)%intensity * normalise
                end do
        endif

        if (calculate_extinction) then
                CALL calc_extinction_coeffs(H_BS, c1, c2, c3, meanextinction, switch_ext, DBLE(10000.),DBLE(1000.), R)

                if (meanextinction .lt. 0.0) then
                   meanextinction = 0.0
                endif

        endif
        !actual dereddening

        if (switch_ext == "S") then
                CALL deredden(ILs, Iint, meanextinction)
                CALL deredden(H_BS, 4, meanextinction)
                call deredden(He_lines, 4, meanextinction)
                CALL deredden(linelist, listlength, meanextinction)
        elseif (switch_ext == "H") then
                CALL deredden_LMC(ILs, Iint, meanextinction)
                CALL deredden_LMC(H_BS, 4, meanextinction)
                call deredden_LMC(He_lines, 4, meanextinction)
                CALL deredden_LMC(linelist, listlength, meanextinction)
        elseif (switch_ext == "C") then
                CALL deredden_CCM(ILs, Iint, meanextinction, R)
                CALL deredden_CCM(H_BS, 4, meanextinction, R)
                call deredden_CCM(He_lines, 4, meanextinction, R)
                CALL deredden_CCM(linelist, listlength, meanextinction, R)
        elseif (switch_ext == "P") then
                CALL deredden_SMC(ILs, Iint, meanextinction)
                CALL deredden_SMC(H_BS, 4, meanextinction)
                call deredden_SMC(He_lines, 4, meanextinction)
                CALL deredden_SMC(linelist, listlength, meanextinction)
        elseif (switch_ext == "F") then
                CALL deredden_Fitz(ILs, Iint, meanextinction)
                CALL deredden_Fitz(H_BS, 4, meanextinction)
                call deredden_Fitz(He_lines, 4, meanextinction)
                CALL deredden_Fitz(linelist, listlength, meanextinction)
        endif


!diagnostics
        call get_diag("ciii1909   ","ciii1907   ", ILs, ciiiNratio)        ! ciii ratio
        call get_diag("oii3729    ","oii3726    ", ILs, oiiNratio )        ! oii ratio
        call get_diag("neiv2425   ","neiv2423   ", ILs, neivNratio )       ! neiv ratio
        call get_diag("sii6731    ","sii6716    ", ILs, siiNratio )        ! s ii ratio
        call get_diag("cliii5537  ","cliii5517  ", ILs, cliiiNratio )      ! Cl iii ratio
        call get_diag("ariv4740   ","ariv4711   ", ILs, arivNratio )       ! Ar iv ratio
        call get_diag("oiii88um   ","oiii52um   ", ILs, oiiiIRNratio )     ! oiii ir ratio
        call get_diag("siii33p5um ","siii18p7um ", ILs, siiiIRNratio )       ! siii ir ratio
        call get_diag("neiii15p5um","neiii36p0um", ILs, neiiiIRNratio )       ! neiii ir ratio
        call get_diag("ariii9um   ","ariii21p8um", ILs, ariiiIRNratio )       ! ariii ir ratio

! temperature ratios:
        !TODO: try to calculate from atomic data at start
        CALL get_Tdiag("nii6548    ","nii6584    ","nii5754    ", ILs, DBLE(4.054), DBLE(1.3274), niiTratio)        ! N II
        CALL get_Tdiag("oiii5007   ","oiii4959   ","oiii4363   ", ILs, DBLE(1.3356), DBLE(3.98), oiiiTratio)        ! O III
        CALL get_Tdiag("neiii3868  ","neiii3967  ","neiii3342  ", ILs, DBLE(1.3013), DBLE(4.319), neiiiTratio)        ! Ne III
        CALL get_Tdiag("neiii3868  ","neiii3967  ","neiii15p5um", ILs, DBLE(1.3013), DBLE(4.319), neiiiIRTratio)! Ne III ir
        CALL get_Tdiag("nev3426    ","nev3345    ","nev2975    ", ILs, DBLE(1.3571), DBLE(3.800), nevTratio)        !!ne v
        CALL get_Tdiag("siii9069   ","siii9531   ","siii6312   ", ILs, DBLE(3.47), DBLE(1.403), siiiTratio)        !s iii
        CALL get_Tdiag("ariii7135  ","ariii7751  ","ariii5192  ",ILs, DBLE(1.24), DBLE(5.174), ariiiTratio)        !ar iii
        CALL get_Tdiag("arv6435    ","arv7005    ","arv4625    ",ILs, DBLE(3.125), DBLE(1.471), arvTratio)        !ar v
        CALL get_Tdiag("ci9850     ","ci9824     ","ci8727     ",ILs, DBLE(1.337), DBLE(3.965), ciTratio)      !C I
        CALL get_Tdiag("oi6364     ","oi6300     ","oi5577     ",ILs, DBLE(4.127), DBLE(1.320), oiTratio)      !O I
        CALL get_Tdiag("oiii4959   ","oiii5007   ","oiii52um   ", ILS, DBLE(1.3356), DBLE(3.98), oiiiIRTratio) ! OIII ir
        !Fixed, DJS

! O II


        if(ILs(get_ion("oii7319b    ",ILs, Iint))%int_dered > 0) then

                ion_no1 = get_ion("oii7319b    ",ILs, Iint)
                ion_no2 = get_ion("oii7330b    ",ILs, Iint)
                ion_no3 = get_ion("oii3726    ",ILs, Iint)
                ion_no4 = get_ion("oii3729    ",ILs, Iint)



                if (ion_no1 .gt. 0 .and. ion_no2 .gt. 0 .and. ion_no3 .gt. 0 .and. ion_no4 .gt. 0) then
                        oiiTratio = (ILs(ion_no1)%int_dered+ILs(ion_no2)%int_dered)/(ILs(ion_no3)%int_dered+ILs(ion_no4)%int_dered)
                else
                        oiiTratio = 0.0
                endif

       elseif(ILs(get_ion("oii7319     ",ILs, Iint))%int_dered > 0)then

                ion_no1 = get_ion("oii7319    ",ILs, Iint)
                ion_no2 = get_ion("oii7320    ",ILs, Iint)
                ion_no3 = get_ion("oii7330    ",ILs, Iint)
                ion_no4 = get_ion("oii7331    ",ILs, Iint)
                ion_no5 = get_ion("oii3726    ",ILs, Iint)
                ion_no6 = get_ion("oii3729    ",ILs, Iint)

                if (ion_no1 .gt. 0 .and. ion_no2 .gt. 0 .and. ion_no3 .gt. 0 .and. ion_no4 .gt. 0 .and. ion_no5 .gt. 0 .and. ion_no6 .gt. 0) then
                        oiiTratio = ((ILs(ion_no1)%int_dered+ILs(ion_no2)%int_dered)+(ILs(ion_no3)%int_dered+ILs(ion_no4)%int_dered))/(ILs(ion_no5)%int_dered+ILs(ion_no6)%int_dered)
                else
                        oiiTratio = 0.0
                endif
        !add condition for 3727 blend


       else
                       oiiTratio=0.0
       endif

! S II
      ion_no1 = get_ion("sii6716    ",ILs, Iint)
      ion_no2 = get_ion("sii6731    ",ILs, Iint)
      ion_no3 = get_ion("sii4068    ",ILs, Iint)
      ion_no4 = get_ion("sii4076    ",ILs, Iint)

      if (ion_no1 .gt. 0 .and. ion_no2 .gt. 0 .and. ion_no3 .gt. 0 .and. ion_no4 .gt. 0) then
           siiTratio = (ILs(ion_no1)%int_dered+ILs(ion_no2)%int_dered)/(ILs(ion_no3)%int_dered+ILs(ion_no4)%int_dered)
      else
           siiTratio = 0.0
      endif

! now get diagnostics zone by zone.

! low ionisation
        ! Edited to stop high limits being included in diagnostic averages. DJS
      lowtemp = 10000.0

      do i = 1,2
        oiiDens=0
        siiDens=0
         count = 0
         if (oiiNratio .gt. 0 .and. oiiNratio .lt. 1e10) then
           call get_diagnostic("oii       ","1,2/                ","1,3/                ",oiiNratio,"D",lowtemp, oiiDens,maxlevs,maxtemps,atomicdata,iion)
           count = count + 1
         endif

         if (siiNratio .gt. 0 .and. siiNratio .lt. 1e10) then
           call get_diagnostic("sii       ","1,2/                ","1,3/                ",siiNratio,"D",lowtemp, siiDens,maxlevs,maxtemps,atomicdata,iion)
           count = count + 1
         endif

         if (count .eq. 0 .and. diagnostic_array(1) .eq. 0) then
           lowdens = 1000.0
         elseif (diagnostic_array(1) .gt. 0.0) then
           lowdens = diagnostic_array(1)
         else
           lowdens = (oiiDens + siiDens) / count
         endif

         count = 0

         if (oiiTratio .gt. 0 .and. oiiTratio .lt. 1e10) then
           call get_diagnostic("oii       ","2,4,2,5,3,4,3,5/    ","1,2,1,3/            ",oiiTratio,"T",lowdens,oiiTemp,maxlevs,maxtemps,atomicdata,iion)
           count = count + 1

                 if(oiitemp .ge. 19999.5)then
                        count=count-1
                        oiitemp=-1
                 endif

         else
           oiiTemp = 0.0
         endif

         if (siiTratio .gt. 0 .and. siiTratio .lt. 1e10) then
           call get_diagnostic("sii       ","1,2,1,3/            ","1,4,1,5/            ",siiTratio,"T",lowdens,siiTemp,maxlevs,maxtemps,atomicdata,iion)
           count = count + 1

                if(siitemp .ge. 19999.5)then
                        count=count-1
                        siitemp=-1
                 endif

         else
           siiTemp = 0.0
         endif

         if (niiTratio .gt. 0 .and. niiTratio .lt. 1e10) then
           call get_diagnostic("nii       ","2,4,3,4/            ","4,5/                ",niiTratio,"T",lowdens,niitemp,maxlevs,maxtemps,atomicdata,iion)
           count = count + 5

                 if(niitemp .ge. 19999.5)then
                        count=count-5
                        niitemp=-1
                 endif

         else
           niitemp = 0.0
         endif

         if (ciTratio .gt. 0 .and. ciTratio .lt. 1e10) then
           call get_diagnostic("ci        ","2,4,3,4/            ","4,5/                ",ciTratio,"T",lowdens,citemp,maxlevs,maxtemps,atomicdata,iion)
           count = count + 1

                if(citemp .ge. 19999.5)then
                        count=count-1
                        citemp=-1
                 endif

          else
           citemp = 0.0
         endif

         if (oiTratio .gt. 0 .and. oiTratio .lt. 1e10) then
           call get_diagnostic("oi        ","1,4,2,4/            ","4,5/                ",oiTratio,"T",lowdens,oitemp,maxlevs,maxtemps,atomicdata,iion)
           count = count + 1

         if(oitemp .ge. 19999.5)then
                 count=count-1
                 oitemp=-1
         endif
         else
           oitemp = 0.0
         endif

         if (count .gt. 0 .and. diagnostic_array(4) .eq. 0) then
           lowtemp = ((5*niitemp) + siitemp + oiitemp + oitemp + citemp) / count
         elseif (diagnostic_array(4) .gt. 0.0) then
           lowtemp = diagnostic_array(4)
         else
           lowtemp = 10000.0
         endif

      enddo

! medium ionisation
        cliiiDens = 0
        ciiiDens = 0
        arivDens = 0


      medtemp = lowtemp

      do i = 1,2

         count = 0
         if (cliiiNratio .gt. 0 .and. cliiiNratio .lt. 1e10) then
           call get_diagnostic("cliii     ","1,2/                ","1,3/                ",cliiiNratio,"D",medtemp, cliiiDens,maxlevs,maxtemps,atomicdata,iion)
           count = count + 1
         else
           cliiiDens=0
         endif
         if (ciiiNratio .gt. 0 .and. ciiiNratio .lt. 1e10) then
           call get_diagnostic("ciii      ","1,2/                ","1,3/                ",ciiiNratio,"D",medtemp, ciiiDens,maxlevs,maxtemps,atomicdata,iion)
           count = count + 1
         else
           ciiiDens=0
         endif
         if (arivNratio .gt. 0 .and. arivNratio .lt. 1e10) then
           call get_diagnostic("ariv      ","1,2/                ","1,3/                ",arivNratio,"D",medtemp, arivDens,maxlevs,maxtemps,atomicdata,iion)
           count = count + 1
         else
           arivDens=0
         endif

! IR densities, not included in average
!Ar, S, Ne, O

         if (oiiiIRNratio .gt. 0 .and. oiiiIRNratio .lt. 1e10) then
           call get_diagnostic("oiii      ","1,2/                ","2,3/                ",oiiiIRNratio,"D",medtemp, oiiiIRDens,maxlevs,maxtemps,atomicdata,iion)
         endif

         if (ariiiIRNratio .gt. 0 .and. ariiiIRNratio .lt. 1e10) then
           call get_diagnostic("ariii     ","1,2/                ","2,3/                ",ariiiIRNratio,"D",medtemp, ariiiIRDens,maxlevs,maxtemps,atomicdata,iion)
         endif

         if (siiiIRNratio .gt. 0 .and. siiiIRNratio .lt. 1e10) then
           call get_diagnostic("siii      ","1,2/                ","2,3/                ",siiiIRNratio,"D",medtemp, siiiIRDens,maxlevs,maxtemps,atomicdata,iion)
         endif

         if (neiiiIRNratio .gt. 0 .and. neiiiIRNratio .lt. 1e10) then
           call get_diagnostic("neiii     ","1,2/                ","2,3/                ",neiiiIRNratio,"D",medtemp, neiiiIRDens,maxlevs,maxtemps,atomicdata,iion)
         endif

         if (count .eq. 0 .and. diagnostic_array(2) .eq. 0) then
           meddens = 1000.0
         elseif (diagnostic_array(2) .gt. 0.0) then
           meddens = diagnostic_array(2)
         else
           meddens = (ciiiDens + cliiiDens + arivDens) / count
         endif

         count = 0

         if (oiiiTratio .gt. 0 .and. oiiiTratio .lt. 1e10) then
           call get_diagnostic("oiii      ","2,4,3,4/            ","4,5/                ",oiiiTratio,"T",meddens,oiiiTemp,maxlevs,maxtemps,atomicdata,iion)
           count = count + 4
                 if(oiiitemp .ge. 19999.5)then
                         count=count-4
                         oiiitemp=-1
                 endif
         else
           oiiiTemp = 0.0
         endif

         if (siiiTratio .gt. 0 .and. siiiTratio .lt. 1e10) then
           call get_diagnostic("siii      ","2,4,3,4/            ","4,5/                ",siiiTratio,"T",meddens,siiiTemp,maxlevs,maxtemps,atomicdata,iion)
           count = count + 1

                if(siiitemp .ge. 19999.5)then
                         count=count-1
                         siiitemp=-1
                 endif

         else
           siiiTemp = 0.0
         endif

         if (ariiiTratio .gt. 0 .and. ariiiTratio .lt. 1e10) then
           call get_diagnostic("ariii     ","1,4,2,4/            ","4,5/                ",ariiiTratio,"T",meddens,ariiitemp,maxlevs,maxtemps,atomicdata,iion)
           count = count + 2
                 if(ariiitemp .ge. 19999.5)then
                         count=count-2
                         ariiitemp=-1
                 endif
         else
           ariiitemp = 0.0
         endif

         if (neiiiTratio .gt. 0 .and. neiiiTratio .lt. 1e10) then
           call get_diagnostic("neiii     ","1,4,2,4/            ","4,5/                ",neiiiTratio,"T",meddens,neiiitemp,maxlevs,maxtemps,atomicdata,iion)
           count = count + 2
                if(neiiitemp .ge. 19999.5)then
                         count=count-2
                         neiiitemp=-1
                 endif

         else
           neiiitemp = 0.0
         endif

!extra IR temperatures, not included in the average at the moment but could be if we decide that would be good

         if (neiiiIRTratio .gt. 0 .and. neiiiIRTratio .lt. 1.e10) then
           call get_diagnostic("neiii     ","1,4,2,4/            ","1,2/                ",neiiiIRTratio,"T",meddens,neiiiIRtemp,maxlevs,maxtemps,atomicdata,iion)
         else
           neiiiIRtemp = 0.0
         endif

         if (oiiiIRTratio .gt. 0 .and. oiiiIRTratio .lt. 1.e10) then
           call get_diagnostic("oiii      ","2,4,3,4/            ","2,3/                ",oiiiIRTratio,"T",oiiiIRdens,oiiiIRtemp,maxlevs,maxtemps,atomicdata,iion)
         else
           oiiiIRtemp = 0.0
         endif


!averaging

         if (count .gt. 0 .and. diagnostic_array(5) .eq. 0.0) then
           medtemp = (4*oiiitemp + siiitemp + 2*ariiitemp + 2*neiiitemp) / count
         elseif (diagnostic_array(5) .gt. 0.0) then
           medtemp = diagnostic_array(5)
         else
           medtemp = lowtemp
         endif

        !dereddening again

        if (calculate_extinction) then

          ILs%int_dered = 0
          H_BS%int_dered = 0
          He_lines%int_dered = 0


        !update extinction. DS 22/10/11
          meanextinction=0
          CALL calc_extinction_coeffs(H_BS, c1, c2, c3, meanextinction, switch_ext, medtemp, lowdens, R)

          if (meanextinction .lt. 0.0) then
             meanextinction = 0.0
          endif

          linelist = linelist_orig
          if (switch_ext == "S") then
                  CALL deredden(ILs, Iint, meanextinction)
                  CALL deredden(H_BS, 4, meanextinction)
                  call deredden(He_lines, 4, meanextinction)
                  CALL deredden(linelist, listlength, meanextinction)
          elseif (switch_ext == "H") then
                  CALL deredden_LMC(ILs, Iint, meanextinction)
                  CALL deredden_LMC(H_BS, 4, meanextinction)
                  call deredden_LMC(He_lines, 4, meanextinction)
                  CALL deredden_LMC(linelist, listlength, meanextinction)
          elseif (switch_ext == "C") then
                  CALL deredden_CCM(ILs, Iint, meanextinction, R)
                  CALL deredden_CCM(H_BS, 4, meanextinction, R)
                  call deredden_CCM(He_lines, 4, meanextinction, R)
                  CALL deredden_CCM(linelist, listlength, meanextinction, R)
          elseif (switch_ext == "P") then
                  CALL deredden_SMC(ILs, Iint, meanextinction)
                  CALL deredden_SMC(H_BS, 4, meanextinction)
                  call deredden_SMC(He_lines, 4, meanextinction)
                  CALL deredden_SMC(linelist, listlength, meanextinction)
          elseif (switch_ext == "F") then
                  CALL deredden_Fitz(ILs, Iint, meanextinction)
                  CALL deredden_Fitz(H_BS, 4, meanextinction)
                  call deredden_Fitz(He_lines, 4, meanextinction)
                  CALL deredden_Fitz(linelist, listlength, meanextinction)
          endif
        endif ! end of exinction calculating

      enddo ! end of diagnostic iteration

        iteration_result(1)%mean_cHb = meanextinction


! high ionisation
      hightemp = medtemp

      do i = 1,2

         if (neivNratio .gt. 0 .and. neivNratio .lt. 1e10) then
           call get_diagnostic("neiv      ","1,2/                ","1,3/                ",neivNratio,"D",hightemp, neivDens,maxlevs,maxtemps,atomicdata,iion)
           highdens = neivdens
         else
           neivDens = 0.0
           highdens = meddens
         endif

         if (diagnostic_array(3) .gt. 0.0) then
           highdens = diagnostic_array(3)
         endif

         count = 0

         if (arvTratio .gt. 0 .and. arvTratio .lt. 1e10) then
           call get_diagnostic("arv       ","2,4,3,4/            ","4,5/                ",arvTratio,"T",highdens,arvTemp,maxlevs,maxtemps,atomicdata,iion)
           count = count + 1
         else
           arvTemp = 0.0
         endif

         if (nevTratio .gt. 0 .and. nevTratio .lt. 1e10) then
           call get_diagnostic("nev       ","2,4,3,4/            ","4,5/                ",nevTratio,"T",highdens,nevtemp,maxlevs,maxtemps,atomicdata,iion)
           count = count + 1
         else
           nevtemp = 0.0
         endif

         if (count .gt. 0 .and. diagnostic_array(6) .eq. 0) then
           hightemp = (arvtemp + nevtemp) / count
         elseif (diagnostic_array(6) .gt. 0.0) then
           hightemp = diagnostic_array(6)
         else
           hightemp = medtemp
         endif

      enddo
!done calculating, now write results and ratios to arrays.

!edited diagnostic section to reflect above changes which stopped high/low limit densities/temperatures being included in averages. High limit cases set to 0.1 so that we know that it was in the high limit and not the low limit. DJS

iteration_result(1)%OII_density_ratio = oiiNratio
iteration_result(1)%OII_density = oiidens
iteration_result(1)%SII_density_ratio = siiNratio
iteration_result(1)%SII_density = siidens

iteration_result(1)%low_density = lowdens
if(niitemp > 0.2)then
        iteration_result(1)%NII_temp = niitemp
else if(INT(niitemp) == -1)then
        iteration_result(1)%NII_temp = 20000
endif
iteration_result(1)%NII_temp_ratio = niiTratio

if(oiitemp >0.2)then
        iteration_result(1)%OII_temp = oiitemp
else if(INT(oiitemp) == -1)then
        iteration_result(1)%OII_temp = 20000
endif
iteration_result(1)%OII_temp_ratio = oiiTratio

if(siitemp >0.2 )then
        iteration_result(1)%SII_temp = siitemp
else if(INT(siitemp) == -1)then
        iteration_result(1)%SII_temp = 20000
endif
iteration_result(1)%SII_temp_ratio = siiTratio

if(oitemp >0.2 )then
        iteration_result(1)%OI_temp = oitemp
else if(INT(oitemp) == -1)then
        iteration_result(1)%OI_temp = 20000
endif
iteration_result(1)%OI_temp_ratio = oiTratio

if(citemp >0.2 )then
        iteration_result(1)%CI_temp = citemp
else if(INT(citemp) == -1)then
        iteration_result(1)%CI_temp = 20000
endif
iteration_result(1)%CI_temp_ratio = ciTratio

iteration_result(1)%low_temp = lowtemp

iteration_result(1)%ClIII_density = cliiidens
iteration_result(1)%ArIV_density = arivdens
iteration_result(1)%CIII_density = ciiidens
iteration_result(1)%OIII_IR_density = oiiiIRdens
iteration_result(1)%ArIII_IR_density = ariiiIRdens
iteration_result(1)%SIII_IR_density = siiiIRdens
iteration_result(1)%NeIII_IR_density = neiiiIRdens

iteration_result(1)%ClIII_density_ratio = cliiiNratio
iteration_result(1)%ArIV_density_ratio = arivNratio
iteration_result(1)%CIII_density_ratio = ciiiNratio
iteration_result(1)%OIII_IR_density_ratio = oiiiIRNratio
iteration_result(1)%ArIII_IR_density_ratio = ariiiIRNratio
iteration_result(1)%SIII_IR_density_ratio = siiiIRNratio
iteration_result(1)%NeIII_IR_density_ratio = neiiiIRNratio


iteration_result(1)%med_density = meddens
if(oiiitemp >0.2)then
        iteration_result(1)%OIII_temp = oiiitemp
else if(INT(oiiitemp) == -1)then
        iteration_result(1)%OIII_temp = 20000
endif
iteration_result(1)%OIII_temp_ratio = oiiiTratio

if(neiiitemp>0.2)then
        iteration_result(1)%NeIII_temp = neiiitemp
else if(INT(neiiitemp) == -1)then
        iteration_result(1)%NeIII_temp = 20000
endif
iteration_result(1)%NeIII_temp_ratio = neiiiTratio

if(ariiitemp>0.2)then
        iteration_result(1)%ArIII_temp = ariiitemp
else if(INT(ariiitemp) == -1)then
        iteration_result(1)%ArIII_temp = 20000
endif
iteration_result(1)%ArIII_temp_ratio = AriiiTratio

if(siiitemp > 0.2)then
        iteration_result(1)%SIII_temp = siiitemp
else if(int(siiitemp) == -1)then
        iteration_result(1)%SIII_temp = 20000
endif
iteration_result(1)%SIII_temp_ratio = siiiTratio

if(oiiiIRtemp > 0.2)then
        iteration_result(1)%OIII_IR_temp = oiiiIRtemp
else if(int(oiiiIRtemp) == -1)then
        iteration_result(1)%OIII_IR_temp = 20000
endif
iteration_result(1)%OIII_IR_temp_ratio = oiiiIRTratio

if(neiiiIRtemp > 0.2)then
        iteration_result(1)%NeIII_IR_temp = neiiiIRtemp
else if(int(neiiiIRtemp) == -1)then
        iteration_result(1)%NeIII_IR_temp = 20000
endif
iteration_result(1)%NeIII_IR_temp_ratio = NeiiiTratio

iteration_result(1)%med_temp = medtemp
iteration_result(1)%NeIV_density = neivdens
iteration_result(1)%high_density = highdens
iteration_result(1)%ArV_temp = arvtemp
iteration_result(1)%NeV_temp = nevtemp
iteration_result(1)%high_temp = hightemp

iteration_result(1)%NeIV_density_ratio = neivNratio
iteration_result(1)%ArV_temp_ratio = arvTratio
iteration_result(1)%NeV_temp_ratio = nevTratio

! Helium abundances

        call get_helium(REAL(medtemp),REAL(meddens),REAL(He_lines(1)%int_dered),REAL(He_lines(2)%int_dered),REAL(He_lines(3)%int_dered),REAL(He_lines(4)%int_dered),heiabund,heiiabund,Hetotabund, A4471, A4686, A6678, A5876)

        w1=0
        w2=0
        w3=0
        w4=0
        if( (A4471 > 0 .or. A5876 > 0 ) .or. A6678 > 0)then

                if(He_lines(2)%intensity > 0) w1 = 1/((He_lines(2)%int_err / He_lines(2)%intensity)**2)
                if(He_lines(3)%intensity > 0) w2 = 1/((He_lines(3)%int_err / He_lines(3)%intensity)**2)
                if(He_lines(4)%intensity > 0) w3 = 1/((He_lines(4)%int_err / He_lines(4)%intensity)**2)

                !PRINT*, w1, " ", w2, " ", w3

                heiabund = (w1*A4471 + w2*A5876 + w3*A6678)/(w1+w2+w3)

        else

                heiabund = 0.0
        endif


! get abundances for all CELs

        !This routine is too simple. I have been changing the temperatures /densities which are input to each zone to disable the zone schtick.
        !Make a better routine that allows using or not using the zone thing.. Its not always appropriate.


        !if(oiiitemp > 0 )then
        !        medtemp = oiiitemp
        !        siiitemp = oiiitemp
        !endif

        !if(int(siiitemp) == -1) siiitemp = oiiitemp


        do i = 1,Iint !This used to be Iint-1 but I think that's corrected in the file reading routine now (RW 25/10/2011)
!                 print *,ILs(i)%ion,ILs(i)%transition,ILs(i)%int_dered
           if (ILs(i)%zone .eq. "low ") then
                !PRINT*, siiitemp, lowdens
                 call get_abundance(ILs(i)%ion, ILs(i)%transition, lowtemp, lowdens,ILs(i)%int_dered, ILs(i)%abundance,maxlevs,maxtemps,atomicdata,iion)
                 ! elseif ( ( i== 47 .or. (i == 46 .or. i == 28 ) ) .and. siiitemp > 1.0 ) then
          !       call get_abundance(ILs(i)%ion, ILs(i)%transition, siiitemp, meddens,ILs(i)%int_dered, ILs(i)%abundance)
                !this makes the code use siii temperatures for siii
                ! print*, "using this bit"
           elseif (ILs(i)%zone .eq. "med ") then
                 call get_abundance(ILs(i)%ion, ILs(i)%transition, medtemp, meddens,ILs(i)%int_dered, ILs(i)%abundance,maxlevs,maxtemps,atomicdata,iion)
           elseif (ILs(i)%zone .eq. "high") then
                 call get_abundance(ILs(i)%ion, ILs(i)%transition, hightemp, highdens,ILs(i)%int_dered, ILs(i)%abundance,maxlevs,maxtemps,atomicdata,iion)
           endif
                 if ((ILs(i)%abundance .ge. 1e-10) .and. (ILs(i)%abundance .lt. 10 ) ) then
!                       PRINT "(1X, A11, 1X, F8.3, 5X, ES10.4)",ILs(i)%name,ILs(i)%int_dered,ILs(i)%abundance
                 elseif( (ILs(i)%abundance .ge. 10 ) .or. (ILs(i)%abundance .lt. 1E-10 ) )then
                         ILs(i)%abundance = 0
                 endif
        enddo

! calculate averages

        celabundtemp = 0.
                niiCELabund = 0.
        weight = 0.
        do i= get_ion("nii5754    ", ILs, Iint), get_ion("nii6584    ", ILs, Iint)
          if (ILs(i)%abundance .ge. 1e-10) niiCELabund = niiCELabund + ILs(i)%abundance/ ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          if (ILs(i)%abundance .ge. 1e-10) weight = weight + 1/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
        enddo
        if (weight .ge. 1e-20) then
          niiCELabund = niiCELabund / weight
        else
          niiCELabund = 0.0
        endif

        niiiIRCELabund = ILs( get_ion("niii57um   ", ILs, Iint)  )%abundance
        niiiUVCELabund = ILs( get_ion("niii1751   ", ILs, Iint)   )%abundance

        if (niiiIRCELabund .ge. 1e-20 .and. niiiUVCELabund .ge. 1e-20) then
          niiiCELabund = (niiiIRCELabund + niiiUVCELabund)/2
        elseif (niiiIRCELabund .ge. 1e-20) then
          niiiCELabund = niiiIRCELabund
        elseif (niiiUVCELabund .ge. 1e-20) then
          niiiCELabund = niiiUVCELabund
        else
          niiiCELabund = 0
        endif

                nivCELabund = 0.
        do i= get_ion("niv1483    ", ILs, Iint), get_ion("niv1485b   ", ILs, Iint) ! would screw up if blend and non blends were both given
          if (ILs(i)%abundance .ge. 1e-20) nivCELabund = nivCELabund + ILs(i)%abundance/ ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          if (ILs(i)%abundance .ge. 1e-20) weight = weight + 1/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
        enddo
        if (weight .ge. 1e-20) then
          nivCELabund = nivCELabund / weight
        else
          nivCELabund = 0.0
        endif

        nvCELabund = ILs(get_ion("nv1240     ",ILs,Iint))%abundance

        !OII CEL routine fixed. DJS 23/11/10

        if(ILs(get_ion("oii3728b   ", ILs, Iint))%int_dered > 0.0)then
                !calc abundance from doublet blend
                oiiCELabund = ILs(get_ion("oii3728b   ", ILs, Iint))%abundance

        else if(ILs(get_ion("oii3729    ", ILs, Iint))%int_dered > 0.0 .and. ILs(get_ion("oii3726    ", ILs, Iint))%int_dered > 0.0 )then
                !calc abundance from doublet
                w1 = 1/(( ILs(get_ion("oii3729    ", ILs, Iint))%int_err / ILs(get_ion("oii3729    ", ILs, Iint))%intensity   )**2)
                w2 = 1/(( ILs(get_ion("oii3726    ", ILs, Iint))%int_err / ILs(get_ion("oii3726    ", ILs, Iint))%intensity   )**2)

                oiiCELabund = (w1*ILs(get_ion("oii3729    ", ILs, Iint))%abundance + w2*ILs(get_ion("oii3726    ", ILs, Iint))%abundance)/(w1+w2)


        else if((ILs(get_ion("oii3728b   ", ILs, Iint))%int_dered == 0.0 .and. (ILs(get_ion("oii3729    ", ILs, Iint))%int_dered ==0.0 .and. ILs(get_ion("oii3726    ", ILs, Iint))%int_dered == 0.0 )) .and.  (ILs(get_ion("oii7330b   ", ILs, Iint))%abundance > 0.0 .or. ILs(get_ion("oii7319b   ", ILs, Iint))%abundance > 0.0)  )then
                !calc abundance based on far red blends
                w1 = 1/(( ILs(get_ion("oii7330b   ", ILs, Iint))%int_err / ILs(get_ion("oii7330b   ", ILs, Iint))%intensity   )**2)
                w2 = 1/(( ILs(get_ion("oii7319b   ", ILs, Iint))%int_err / ILs(get_ion("oii7319b   ", ILs, Iint))%intensity   )**2)

                oiiCELabund = (w1*ILs(get_ion("oii7330b   ", ILs, Iint))%abundance + w2*ILs(get_ion("oii7319b   ", ILs, Iint))%abundance)/(w1+w2)


        else if        ((ILs(get_ion("oii3728b   ", ILs, Iint))%int_dered == 0.0 .and. (ILs(get_ion("oii3729    ", ILs, Iint))%int_dered ==0.0 .and. ILs(get_ion("oii3726    ", ILs, Iint))%int_dered == 0.0 )) .and. ( (ILs(get_ion("oii7320    ", ILs, Iint))%abundance > 0.0 .or. ILs(get_ion("oii7319    ", ILs, Iint))%abundance > 0.0) .or.  (ILs(get_ion("oii7330    ", ILs, Iint))%abundance > 0.0 .or. ILs(get_ion("oii7331    ", ILs, Iint))%abundance > 0.0)) )then
                !calc abundance based on far red quadruplet
                if(ILs(get_ion("oii7319    ", ILs, Iint))%int_err > 0) w1 = 1/(( ILs(get_ion("oii7319    ", ILs, Iint))%int_err / ILs(get_ion("oii7319    ", ILs, Iint))%intensity   )**2)
                if(ILs(get_ion("oii7320    ", ILs, Iint))%int_err > 0) w2 = 1/(( ILs(get_ion("oii7320    ", ILs, Iint))%int_err / ILs(get_ion("oii7320    ", ILs, Iint))%intensity   )**2)
                if(ILs(get_ion("oii7330    ", ILs, Iint))%int_err > 0) w3 = 1/(( ILs(get_ion("oii7330    ", ILs, Iint))%int_err / ILs(get_ion("oii7330    ", ILs, Iint))%intensity   )**2)
                if(ILs(get_ion("oii7331    ", ILs, Iint))%int_err > 0) w4 = 1/(( ILs(get_ion("oii7331    ", ILs, Iint))%int_err / ILs(get_ion("oii7320    ", ILs, Iint))%intensity   )**2)

                !if statements stop non existent lines being granted infinite weight 1/(0/0)^2 = infinity, defaults to zero if no line detected which keeps the following calculation honest

                oiiCELabund = (w1*ILs(get_ion("oii7319    ", ILs, Iint))%abundance + w2*ILs(get_ion("oii7320    ", ILs, Iint))%abundance + w3*ILs(get_ion("oii7330    ", ILs, Iint))%abundance + w4*ILs(get_ion("oii7331    ", ILs, Iint))%abundance)/(w1+w2+w3+w4)

        else
                oiiCELabund = 0.0
        endif

        !The above routine for oii replaces the following as it was decided that the weighting scheme only works if the lines are originating from the same energy levels.

!       celabundtemp = 0.
!        weight = 0.
!        if( ILs(get_ion("oii7330b   ", ILs, Iint))%int_dered > 0.0  )then
!                 do i=get_ion("oii3729    ", ILs, Iint), get_ion("oii7330b   ", ILs, Iint)
!                  if (ILs(i)%abundance .gt. 0) oiiCELabund = oiiCELabund + ILs(i)%abundance/ ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
!                  if (ILs(i)%abundance .gt. 0) weight = weight + 1/ ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
!                enddo
!                if (weight .gt. 0) then
!                  oiiCELabund = oiiCELabund / weight
!                else
!                  oiiCELabund = 0.0
!                endif
!        elseif( ILs(get_ion("oii7330    ", ILs, Iint))%int_dered >0   ) then
!                do i=get_ion("oii3729    ", ILs, Iint), get_ion("oii3726    ", ILs, Iint)
!                  if (ILs(i)%abundance .gt. 0) oiiCELabund = oiiCELabund + ILs(i)%abundance/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
!                  if (ILs(i)%abundance .gt. 0) weight = weight + 1/ ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
!                enddo
!                do i=get_ion("oii7320    ", ILs, Iint), get_ion("oii7330    ", ILs, Iint)
!                  if (ILs(i)%abundance .gt. 0) oiiCELabund = oiiCELabund + ILs(i)%abundance/ ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
!                  if (ILs(i)%abundance .gt. 0) weight = weight + 1/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
!                enddo
!
!
!                if (weight .gt. 0) then
!                  oiiCELabund = oiiCELabund / weight
!                else
!                  oiiCELabund = 0.0
!                endif
!        else
!                oiiCELabund = 0.0
!        endif

        celabundtemp = 0.
                oiiiCELabund = 0.0
        weight = 0.
        do i=get_ion("oiii4959   ", ILs, Iint), get_ion("oiii5007   ", ILs, Iint)
          if (ILs(i)%abundance .ge. 1e-20) oiiiCELabund = oiiiCELabund + ILs(i)%abundance /((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          if (ILs(i)%abundance .ge. 1e-20) weight = weight + 1 / ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
        enddo
        if (weight .ge. 1e-20) then
          oiiiCELabund = oiiiCELabund / weight
        else
          oiiiCELabund = 0.0
        endif

        celabundtemp = 0.
                oiiiIRCELabund = 0.0
        weight = 0.
        do i=get_ion("oiii52um   ", ILs, Iint), get_ion("oiii88um   ", ILs, Iint)
          if (ILs(i)%abundance .ge. 1e-20) oiiiIRCELabund = oiiiIRCELabund + ILs(i)%abundance / ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          if (ILs(i)%abundance .ge. 1e-20) weight = weight + 1/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
        enddo
        if (weight .ge. 1e-20) then
          oiiiIRCELabund = oiiiIRCELabund / weight
        else
          oiiiIRCELabund = 0.0
        endif

        oivCELabund = ILS( get_ion("oiv25p9um  ", ILS, Iint))%abundance

        neiiIRCELabund = ILs(  get_ion("neii12p8um ", ILs, Iint)  )%abundance
        neiiiIRCELabund = ILs(  get_ion("neiii15p5um ", ILs, Iint)  )%abundance

        celabundtemp = 0.
                neiiiCELabund = 0.
        weight = 0.
        do i=get_ion("neiii3868  ", ILs, Iint), get_ion("neiii3967  ", ILs, Iint)
          if (ILs(i)%abundance .ge. 1e-20) neiiiCELabund = neiiiCELabund + ILs(i)%abundance/ ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          if (ILs(i)%abundance .ge. 1e-20) weight = weight + 1/ ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
        enddo
        if (weight .ge. 1e-20) then
          neiiiCELabund = neiiiCELabund / weight
        else
          neiiiCELabund = 0.0
        endif


        celabundtemp = 0.
                neivCELabund = 0.
        weight = 0.
        do i=get_ion("neiv2423   ", ILs, Iint), get_ion("neiv4725b  ", ILs, Iint) ! would screw up if blends and non blends were given
          if (ILs(i)%abundance .ge. 1e-20) neivCELabund = neivCELabund + ILs(i)%abundance/ ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          if (ILs(i)%abundance .ge. 1e-20) weight = weight + 1/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
        enddo
        if (weight .ge. 1e-20) then
          neivCELabund = neivCELabund / weight
        else
          neivCELabund = 0.0
        endif

        celabundtemp = 0.
                siiCELabund = 0.
        weight = 0.
        do i=get_ion("sii4068    ", ILs, Iint), get_ion("sii6731    ", ILs, Iint)
          if (ILs(i)%abundance .ge. 1e-20) siiCELabund = siiCELabund + ILs(i)%abundance/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          if (ILs(i)%abundance .ge. 1e-20) weight = weight + 1/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
        enddo
        if (weight .ge. 1e-20) then
          siiCELabund = siiCELabund / weight
        else
          siiCELabund = 0.0
        endif

        !siiiCELabund = 0 ! ILs(28)%abundance
        !celabundtemp = 0.
        !weight = 0.
        !do i=get_ion("siii9069   ", ILs, Iint), get_ion("siii9531   ", ILs, Iint)
        !  if (ILs(i)%abundance .gt. 0) siiiCELabund = siiiCELabund + ILs(i)%abundance/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
        !  if (ILs(i)%abundance .gt. 0) weight = weight + 1/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
        !enddo
             !if (weight .gt. 0) then
        !  siiiCELabund = siiiCELabund / weight
        !else
        !  siiiCELabund = 0.0
        !endif

        ! SIII abundance section previously buggered. SIII 9069 and 9531 doublet special due to absorption lines. See Liu, Barlow etc 1995 (Far Red IR Lines in a PN), one of 9069/9531 is ALWAYS absorbed however the other is then always fine.  We can tell which is is by looking at the ratio 9069/9531.
        if (ILs(get_ion("siii9069   ", ILs, Iint))%abundance .eq. 0 .and. ILs(get_ion("siii9531   ", ILs, Iint))%abundance .eq. 0) then
                siiiCELabund = ILs(get_ion("siii6312   ", ILs, Iint))%abundance
        else !only calculate all the IR SIII telluric absorption if the lines are present.
                siiiCELabund = ILs(get_ion("siii9069   ", ILs, Iint))%int_dered / ILs(get_ion("siii9531   ", ILs, Iint))%int_dered
        !I am using siiiCELabund as a switch for the following if statement, replace this with another variable if you like but I don't think it matters.. (DJS)

                if(siiiCELabund < (1.05 * 0.403) .and. siiiCELabund > (0.90 * 0.403))then !this case should never occur

                        if(ILs(get_ion("siii9069   ", ILs, Iint))%intensity > 0) w1 = 1/(ILs(get_ion("siii9069   ", ILs, Iint))%int_err / ILs(get_ion("siii9069   ", ILs, Iint))%intensity)**2
                        if(ILs(get_ion("siii9531   ", ILs, Iint))%intensity > 0) w2 = 1/(ILs(get_ion("siii9531   ", ILs, Iint))%int_err / ILs(get_ion("siii9531   ", ILs, Iint))%intensity)**2

                        siiiCELabund= ( w1*ILs(get_ion("siii9069   ", ILs, Iint))%abundance + w2*ILs(get_ion("siii9531   ", ILs, Iint))%abundance )/(w1+w2)

                elseif(siiiCELabund > (1.1 * 0.403) )then
                        !9531 absorbed
                        siiiCELabund = ILs( get_ion("siii9069   ", ILs, Iint) )%abundance

                elseif(siiiCELabund < (0.9 * 0.403) )then
                        !9069 absorbed
                        siiiCELabund = ILs( get_ion("siii9531   ", ILs, Iint) )%abundance
                else
                        siiiCELabund=0.0
                endif
        endif

        siiiIRCELabund = ILs( get_ion("siii18p7um ", ILs, Iint)   )%abundance
        sivIRCELabund = ILs(  get_ion("siv10p5um  ", ILs, Iint) )%abundance

        celabundtemp = 0.
                cliiiCELabund = 0.
        weight = 0.
        do i=get_ion("cliii5517  ", ILs, Iint), get_ion("cliii5537  ", ILs, Iint)
          if (ILs(i)%abundance .ge. 1e-20) cliiiCELabund = cliiiCELabund + ILs(i)%abundance/ ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          if (ILs(i)%abundance .ge. 1e-20) then
            weight = weight + 1/ ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          endif
        enddo
        if (weight .ge. 1e-20) then
          cliiiCELabund = cliiiCELabund / weight
        else
          cliiiCELabund = 0.0
        endif

        celabundtemp = 0.
                ariiiCELabund = 0.
        weight = 0.
        do i=get_ion("ariii7135  ", ILs, Iint), get_ion("ariii7751  ", ILs, Iint)
          if (ILs(i)%abundance .ge. 1e-20) ariiiCELabund = ariiiCELabund + ILs(i)%abundance/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          if (ILs(i)%abundance .ge. 1e-20) then
            weight = weight + 1/ ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          endif
        enddo
        if (weight .ge. 1e-20) then
          ariiiCELabund = ariiiCELabund / weight
        else
          ariiiCELabund = 0.0
        endif

        celabundtemp = 0.
                arivCELabund = 0.
        weight = 0.
        do i=get_ion("ariv4711   ", ILs, Iint), get_ion("ariv4740   ", ILs, Iint)
        arivCELabund = arivCELabund + ILs(i)%abundance/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          if (ILs(i)%abundance .ge. 1e-20) then
            weight = weight + 1/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          endif
        enddo
        if (weight .ge. 1e-20) then
          arivCELabund = arivCELabund / weight
        else
          arivCELabund = 0.0
        endif

        ariiiIRCELabund = ILs(get_ion("ariii9um   ", ILs, Iint))%abundance

        ciiCELabund = ILs(get_ion("cii2325    ", ILs, Iint))%abundance
        civCELabund = ILs(get_ion("civ1548    ", ILs, Iint))%abundance

        celabundtemp = 0.
                ciiiCELabund = 0.
        weight = 0.
        do i=get_ion("ciii1907   ", ILs, Iint), get_ion("ciii1909b  ", ILs, Iint) ! would screw up if blend and non-blend were given.
          if (ILs(i)%abundance .ge. 1e-20) then
            ciiiCELabund = ciiiCELabund + ILs(i)%abundance/ ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
            weight = weight + 1/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          endif
        enddo
        if (weight .ge. 1e-20) then
          ciiiCELabund = ciiiCELabund / weight
        else
          ciiiCELabund = 0.0
        endif

        celabundtemp = 0.
                neivCELabund = 0.
        weight = 0.
        do i=get_ion("neiv2423   ", ILs, Iint), get_ion("neiv2425   ", ILs, Iint)
          if (ILs(i)%abundance .ge. 1e-20) neivCELabund = neivCELabund + ILs(i)%abundance/ ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          if (ILs(i)%abundance .ge. 1e-20) then
            weight = weight + 1/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          endif
        enddo
        if (weight .ge. 1e-20) then
          neivCELabund = neivCELabund / weight
        else
          neivCELabund = 0.0
        endif

        celabundtemp = 0.
                nevCELabund = 0.
        weight = 0.
        do i=get_ion("nev3345    ", ILs, Iint), get_ion("nev3426    ", ILs, Iint)
          if (ILs(i)%abundance .ge. 1e-20) nevCELabund = nevCELabund + ILs(i)%abundance/ ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          if (ILs(i)%abundance .ge. 1e-20) then
            weight = weight + 1/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          endif
        enddo
        if (weight .ge. 1e-20) then
          nevCELabund = nevCELabund / weight
        else
          nevCELabund = 0.0
        endif

        celabundtemp = 0.
                arvCELabund = 0.
        weight = 0.
        do i=get_ion("arv6435    ", ILs, Iint), get_ion("arv7005    ", ILs, Iint)
          if (ILs(i)%abundance .ge. 1e-20) arvCELabund = arvCELabund + ILs(i)%abundance/ ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          if (ILs(i)%abundance .ge. 1e-20) then
            weight = weight + 1/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          endif
        enddo
        if (weight .ge. 1e-20) then
          arvCELabund = arvCELabund / weight
        else
          arvCELabund = 0.0
        endif

         celabundtemp = 0.
                 ciCELabund = 0.
        weight = 0.
        do i=get_ion("ci9850     ", ILs, Iint), get_ion("ci8727     ", ILs, Iint)
          if (ILs(i)%abundance .ge. 1e-20) ciCELabund = ciCELabund + ILs(i)%abundance/ ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          if (ILs(i)%abundance .ge. 1e-20) then
            weight = weight + 1/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          endif
        enddo
        if (weight .ge. 1e-20) then
          ciCELabund = ciCELabund / weight
        else
          ciCELabund = 0.0
        endif

        NCabundCEL = ciCELabund

         celabundtemp = 0.
                 oiCELabund = 0.
        weight = 0.
        do i=get_ion("oi6300     ", ILs, Iint), get_ion("oi5577     ", ILs, Iint)
          if (ILs(i)%abundance .ge. 1e-20) oiCELabund = oiCELabund + ILs(i)%abundance/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          if (ILs(i)%abundance .ge. 1e-20) then
            weight = weight + 1/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          endif
        enddo
        if (weight .ge. 1e-20) then
          oiCELabund = oiCELabund / weight
        else
          oiCELabund = 0.0
        endif

        NOabundCEL = oiCELabund

! now get abundances for ORLs
! o2+
       call oii_rec_lines(oiiiTemp,oiiDens,DBLE(1),oiiRLs)

       do i = 1,listlength
         do j = 1,415
          if (abs(linelist(i)%wavelength-oiiRLs(j)%Wave) .le. 0.005) then
            oiiRLs(j)%Obs = linelist(i)%int_dered
            oiiRLs(j)%abundance = oiiRLs(j)%obs/oiiRLs(j)%Int
          endif
         enddo
       enddo

!N2+

       call nii_rec_lines(oiiiTemp,oiiDens,DBLE(1),niiRLs)

       do i = 1,listlength
         do j = 1,99
          if (abs(linelist(i)%wavelength-niiRLs(j)%Wave) .le. 0.005) then
            niiRLs(j)%Obs = linelist(i)%int_dered
            niiRLs(j)%abundance = niiRLs(j)%obs/niiRLs(j)%Int
          endif
         enddo
       enddo

!C2+
       call cii_rec_lines(oiiiTemp,oiiDens,DBLE(1),ciiRLs)

       do i = 1,listlength
         do j = 1,57
          if (abs(linelist(i)%wavelength-ciiRLs(j)%Wave) .le. 0.005) then
            ciiRLs(j)%Obs = linelist(i)%int_dered
            ciiRLs(j)%abundance = ciiRLs(j)%obs/ciiRLs(j)%Int
          endif
         enddo
       enddo

!Ne2+
       call neii_rec_lines(oiiiTemp,oiiDens,DBLE(1),neiiRLs)

       do i = 1,listlength
         do j = 1,38
          if (abs(linelist(i)%wavelength-neiiRLs(j)%Wave) .le. 0.005) then
            neiiRLs(j)%Obs = linelist(i)%int_dered
            neiiRLs(j)%abundance = neiiRLs(j)%obs/neiiRLs(j)%Int
          endif
         enddo
       enddo

!C3+, N3+
       call xiii_rec_lines(oiiiTemp,oiiDens,DBLE(1),xiiiRLs)

       do i = 1,listlength
         do j = 1,6
          if (abs(linelist(i)%wavelength-xiiiRLs(j)%Wave) .le. 0.005) then
            xiiiRLs(j)%Obs = linelist(i)%int_dered
            xiiiRLs(j)%abundance = xiiiRLs(j)%obs/xiiiRLs(j)%Int
          endif
         enddo
       enddo

      rlabundtemp = 0.0
      weight = 0.0

!cii recombination lines

      do i = 1,57
        if (ciiRLs(i)%abundance .ge. 1e-20) then
!          print "(1X,F7.2,1X,F6.3,1X,ES9.3)",ciiRLs(i)%wave,ciiRLs(i)%obs,ciiRLs(i)%abundance
          rlabundtemp = rlabundtemp + ciiRLs(i)%obs
          weight = weight + ciiRLs(i)%Int
        endif
        if (ciiRLs(i)%wave .eq. 4267.15D0) then
          cii4267rlabund = ciiRLs(i)%abundance
        endif
      enddo

      if (weight .gt. 0) then
        ciirlabund = rlabundtemp/weight
      else
        ciirlabund = 0.
      endif

!nii recombination lines

!      print *,"lambda   Mult   Int   Abund"
      do i = 1,99
        if (niiRLs(i)%abundance .ge. 1e-20) then
!          print "(F7.2,1X,A7,1X,F6.3,1X,ES9.3)",niiRLs(i)%wave,niiRLs(i)%Mult,niiRLs(i)%obs,niiRLs(i)%abundance
          rlabundtemp = rlabundtemp + niiRLs(i)%obs
          weight = weight + niiRLs(i)%Int
        endif
      enddo

  if (weight .gt. 0) then
!      print "(ES9.2)",rlabundtemp / weight

      niimultiplets%Multiplet = (/"V3     ","V5     ","V8     " ,"V12    ","V20    ","V28    ","3d-4f  "/)

! get multiplet abundances from coadded intensity

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
      if (isnan((rlabundtemp/weight))) then
        niimultiplets(j)%Abundance = 0.
      else
        niimultiplets(j)%Abundance = rlabundtemp/weight
      endif
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

      if (isnan((rlabundtemp/weight))) then
      niimultiplets(7)%abundance = 0
      else
      niimultiplets(7)%abundance = rlabundtemp/weight
      endif

!      print "(F6.3,16X,ES9.3)",rlabundtemp, rlabundtemp/weight

      rlabundtemp = 0.0
      weight = 0
      do i = 1,7
        rlabundtemp = rlabundtemp + niimultiplets(i)%abundance
        if (niimultiplets(i)%abundance .ge. 1e-20) then
          weight = weight + 1
        endif
      enddo

      if (weight .gt. 0) then
        niiRLabund = rlabundtemp/weight
      else
        niiRLabund = 0.D0
      endif
  endif

!oii recombination lines

      rlabundtemp = 0.00
      weight = 0.00

!      print *,"lambda   Mult   Int   Abund"
      do i = 1,415
        if (oiiRLs(i)%abundance .ge. 1e-20) then
!          print "(F7.2,1X,A7,1X,F6.3,1X,ES9.3)",oiiRLs(i)%wave,oiiRLs(i)%Mult,oiiRLs(i)%obs,oiiRLs(i)%abundance
          rlabundtemp = rlabundtemp + oiiRLs(i)%obs
          weight = weight + oiiRLs(i)%Int
        endif
      enddo

  if (weight .gt. 0) then

!     print "(ES9.2)",rlabundtemp/weight

      oiimultiplets%Multiplet = (/" V1    "," V2    "," V5    " ," V10   "," V11   "," V12   "," V19   "," V20   "," V25   "," V28   "," V33   "," 3d-4f "/)

! get multiplet abundances from coadded intensity

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

      if ( isnan( (rlabundtemp/weight) ) ) then
      oiimultiplets(12)%abundance = 0
      else
      oiimultiplets(12)%abundance = rlabundtemp/weight
      endif

!      print *,"3d-4f :"
!      print *,"Co-added intensity   O2+/H+"
!      print "(F6.3,16X,ES9.3)",rlabundtemp, rlabundtemp/weight

      rlabundtemp = 0.0
      weight = 0
      do i = 1,7
        rlabundtemp = rlabundtemp + oiimultiplets(i)%abundance
        if (oiimultiplets(i)%abundance .ge. 1e-20) then
          weight = weight + 1
        endif
      enddo

      if (weight .gt. 0) then
        oiiRLabund = rlabundtemp/weight
      else
        oiiRLabund = 0.D0
      endif

   endif

!neii recombination lines

      rlabundtemp = 0.0
      weight = 0.0

!      print *,"lambda   Mult   Int   Abund"
      do i = 1,38
        if (neiiRLs(i)%abundance .ge. 1e-20) then
!          print "(F7.2,1X,F6.3,1X,ES9.3)",neiiRLs(i)%wave,neiiRLs(i)%obs,neiiRLs(i)%abundance
           rlabundtemp = rlabundtemp + neiiRLs(i)%obs
           weight = weight + neiiRLs(i)%Int
        endif
      enddo

   if (weight .gt. 0) then
      neiiRLabund = rlabundtemp/weight
   else
      neiiRLabund = 0.D0
   endif


      rlabundtemp = 0.0
      weight = 0.0

      do i = 1,4
        if (xiiiRLs(i)%abundance .ge. 1e-20) then
          rlabundtemp = rlabundtemp + xiiiRLs(i)%obs
          weight = weight + xiiiRLs(i)%Int
        endif
      enddo
      if (weight .gt. 0) then
        ciiiRLabund = rlabundtemp / weight
      else
        ciiiRLabund = 0.0
      endif

!      do i = 5,6
!        if (xiiiRLs(i)%abundance .ge. 1e-20) then
!           print "(F7.2,1X,F6.3,1X,ES9.3)",xiiiRLs(i)%wave,xiiiRLs(i)%obs,xiiiRLs(i)%abundance
!        endif
!      enddo

     if (xiiiRLs(6)%abundance .ge. 1e-20) then
        niiiRLabund = xiiiRLs(6)%abundance
     else
        niiiRLabund = 0.0
     endif

! ICFs (Kingsburgh + Barlow 1994)

! oxygen - complete
     OabundCEL = 0.
        if (oiiCELabund .ge. 1e-20 .and. oiiiCELabund .ge. 1e-20 .and. oivCELabund .ge. 1e-20 .and. nvCELabund .ge. 1e-20)then ! O3+ and N4+
                fn4 = (nvCELabund)/(niiCELabund + niiiCELabund + nivCELabund + nvCELabund) !A4
                CELicfO = 1./(1.-0.95*fn4)                                                 !A5
                OabundCEL = CELicfO * (oiiCELabund + oiiiCELabund + oivCELabund)           !A6
        elseif (oiiCELabund .ge. 1e-20 .and. oiiiCELabund .ge. 1e-20 .and. oivCELabund .ge. 1e-20 .and. nvCELabund .lt. 1e-20) then ! O3+ but no N4+
                CELicfO = 1
                OabundCEL = oiiCELabund + oiiiCELabund + oivCELabund !sum of visible ionisation stages
        elseif (oiiCELabund .ge. 1e-20 .and. oiiiCELabund .ge. 1e-20 .and. oivCELabund .lt. 1e-20 .and. nvCELabund .ge. 1e-20) then !no O3+ but N4+
                CELicfO = (niiCELabund + niiiCELabund + nivCELabund + nvCELabund) / (niiCELabund + niiiCELabund) ! A7
                OabundCEL = CELicfO * (oiiCELabund + oiiiCELabund)                                              ! A8
        elseif (oiiCELabund .ge. 1e-20 .and. oiiiCELabund .ge. 1e-20 .and. oivCELabund .lt. 1e-20 .and. nvCELabund .lt. 1e-20)then !no O3+ or N4+ seen
                CELicfO = ((heiabund + heiiabund)/heiabund)**(2./3.) !KB94 A9
                OabundCEL = CELicfO * (oiiCELabund + oiiiCELabund) !A10
        endif

! nitrogen - complete
     NabundCEL = 0.
     if (niiCELabund .ge. 1e-20 .and. niiiUVCELabund .ge. 1e-20 .and. nivCELabund .ge. 1e-20) then !all ionisation stages seen
       CELicfN = 1.
       NabundCEL = niiCELabund + niiiCELabund + nivCELabund
     elseif (niiCELabund .ge. 1e-20 .and. niiiUVCELabund .lt. 1e-20 .and. nivCELabund .ge. 1e-20) then !no N2+ seen
       CELicfN = 1.5
       NabundCEL = 1.5*(niiCELabund + nivCELabund)
     elseif (niiCELabund .ge. 1e-20 .and. niiiUVCELabund .lt. 1e-20 .and. nivCELabund .lt. 1e-20) then !Only N+ seen
       CELicfN = OabundCEL/oiiCELabund
       NabundCEL = niiCELabund * CELicfN
     endif

! carbon - complete
     CabundCEL = 0.
     if (ciiCELabund .ge. 1e-20 .and. ciiiCELabund .ge. 1e-20 .and. civCELabund .ge. 1e-20 .and. heiiabund .lt. 1e-20) then !No C4+ but all other seen
       CELicfC = 1.
       CabundCEL = ciiCELabund + ciiiCELabund + civCELabund
     elseif (ciiCELabund .lt. 1e-20 .and. ciiiCELabund .ge. 1e-20 .and. civCELabund .ge. 1e-20 .and. heiiabund .lt. 1e-20) then !No C4+, and no CII lines seen
       CELicfC = (oiiCELabund + oiiiCELabund) / oiiiCELabund !A11
       CabundCEL = CELicfC * (ciiCELabund + ciiiCELabund) !A12
     elseif (ciiCELabund .lt. 1e-20 .and. ciiiCELabund .ge. 1e-20 .and. civCELabund .lt. 1e-20 .and. oiiiCELabund .ge. 1e-20) then !Only C2+ seen, but O2+ also seen
       CELicfC = OabundCEL/oiiiCELabund !A13
       CabundCEL = CELicfC * ciiiCELabund !A14
     elseif (nvCELabund .ge. 1e-20 .and. heiiabund .ge. 1e-20) then !N4+ and He2+
       fn4 = (nvCELabund)/(niiCELabund + niiiCELabund + nivCELabund + nvCELabund) !A4, A15
       if (fn4 .lt. 0.29629) then !condition in KB94 is if icf(C)>5, but eqn A16 can go negative at high fn4 so this is a better check for the high-excitation case
         CELicfC = 1/(1-(2.7*fn4))!A16
       else
         CELicfC = (niiCELabund + niiiCELabund + nivCELabund + nvCELabund) / (niiCELabund + niiiCELabund + nivCELabund) !A18
       endif
       CabundCEL = CELicfC * (ciiCELabund + ciiiCELabund + civCELabund) !A19
     elseif (ciiCELabund .ge. 1e-20 .and. ciiiCELabund .ge. 1e-20 .and. civCELabund .ge. 1e-20 .and. heiiabund .ge. 1e-20 .and. nvCELabund .lt. 1e-20) then !PN is hot enough for He2+ but not for N4+
        CELicfC = ((heiiabund + heiabund)*heiabund)**(1./3.) !A20
        CabundCEL = CELicfC * (ciiCELabund + ciiiCELabund + civCELabund) !A21
     elseif (ciiCELabund .lt. 1e-20 .and. ciiiCELabund .ge. 1e-20 .and. civCELabund .ge. 1e-20) then !C+ not seen
        CELicfC = 1/(1 - (niiCELabund/NabundCEL) - 2.7*(nvCELabund/NabundCEL)) !A22
        if (CELicfC .gt. 5) then !if ICF is greater than 5
           CELicfC = (niiCELabund + niiiCELabund + nivCELabund + nvCELabund)/(niiiCELabund + nivCELabund) !A23
        endif
        CabundCEL = CELicfC * (ciiiCELabund+ civCELabund) !A24
     elseif (ciiCELabund .ge. 1e-20 .and. ciiiCELabund .ge. 1e-20 .and. civCELabund .ge. 1e-20 .and. heiiabund .ge. 1e-20 .and. (nivCELabund .lt. 1e-20 .or. nvCELabund .lt. 1e-20)) then !final case i think.
        CELicfC = ((oiiCELabund + oiiiCELabund)/oiiiCELabund)*((heiiabund + heiabund)*heiabund)**(1./3.) !A25
        CabundCEL = CELicfC * (ciiCELabund + ciiiCELabund + civCELabund) !A26
     else
        CabundCEL = 0.D0 ! nothing at all seen
        CELicfC = 1.0
     endif

! Neon - complete
     NeabundCEL = 0.
     if (neiiiCELabund .ge. 1e-20 .and. neivCELabund .ge. 1e-20 .and. nevCELabund .ge. 1e-20) then !all stages seen
       CELicfNe = 1.
       NeabundCEL = neiiiCELabund + neivCELabund + nevCELabund
     elseif (neiiiCELabund .ge. 1e-20 .and. neivCELabund .lt. 1e-20 .and. nevCELabund .ge. 1e-20) then !no Ne IV seen
       CELicfNe = 1.5
       NeabundCEL = CELicfNe * (neiiiCELabund + nevCELabund) !KB94 A27
     elseif (neiiiCELabund .ge. 1e-20 .and. neivCELabund .lt. 1e-20 .and. nevCELabund .lt. 1e-20) then !Only Ne2+ seen
       CELicfNe = OabundCEL / oiiiCELabund !KB94 A28
       NeabundCEL = CELicfNe * neiiiCELabund
     endif

! Argon - complete
     ArabundCEL = 0.
     if (ariiiCELabund .ge. 1e-20 .and. arivCELabund .lt. 1e-20 .and. arvCELabund .lt. 1e-20) then !only Ar2+ seen
       CELicfAr = 1.87 !KB94 A32
       ArabundCEL = CELicfAr * ariiiCELabund !KB94 A33
     elseif (ariiiCELabund .lt. 1e-20 .and. arivCELabund .ge. 1e-20 .and. arvCELabund .lt. 1e-20) then !Only Ar3+ seen
       CELicfAr = NeabundCEL / neiiiCELabund !KB94 A34
       ArabundCEL = CELicfAr * arivCELabund !KB94 A35
     elseif (ariiiCELabund .ge. 1e-20 .and. arivCELabund .ge. 1e-20) then !Ar 2+ and 3+ seen
       CELicfAr = 1./(1.-(niiCELabund/NabundCEL))
       ArabundCEL = ariiiCELabund + arivCELabund + arvCELabund !KB94 A31
     endif

! Sulphur
     SabundCEL = 0.
     if (siiCELabund .ge. 1e-20 .and. siiiCELabund .ge. 1e-20 .and. sivIRCELabund .lt. 1e-20) then !both S+ and S2+
       CELicfS = (1 - (  (1-(oiiCELabund/OabundCEL))**3.0  )   )**(-1.0/3.0) !KB94 A36
       SabundCEL = CELicfS * (siiCELabund + siiiCELabund) !KB94 A37
     elseif (siiCELabund .ge. 1e-20 .and. siiiCELabund .ge. 1e-20 .and. sivIRCELabund .ge. 1e-20) then !all states observed
       CELicfS = 1.
       SabundCEL = siiCELabund + siiiCELabund + sivIRCELabund
     elseif (siiCELabund .ge. 1e-20 .and. siiiCELabund .lt. 1e-20 .and. sivIRCELabund .lt. 1e-20) then !Only S+ observed
       CELicfS = (((oiiiCELabund/oiiCELabund)**0.433) * 4.677) * (1-(1-((oiiCELabund/OabundCEL)**3)))**(-1./3.) ! KB94 A37 with S2+/S+ from A38
       SabundCEL = CELicfS * siiCELabund
     endif

!very high excitation cases, all He is doubly ionised

     if (heiabund .lt. 1e-20 .and. heiiabund .ge. 1e-20) then
       CELicfO = NeabundCEL / neiiiCELabund                              !A39
       CELicfC = CELicfO                                                 !A39
       OabundCEL = CELicfO * (oiiCELabund + oiiiCELabund + oivCELabund)  !A40
       CabundCEL = CELicfC * (ciiiCELabund + civCELabund)                !A41
     endif

!Chlorine - not included in KB94, this prescription is from Liu et al. (2000)
     ClabundCEL = 0.
    if (cliiiCELabund .ge. 1e-20 .and. siiiCELabund .ge. 1e-20) then
       CELicfCl = SabundCEL/siiiCELabund
       ClabundCEL = CELicfCl * cliiiCELabund
    endif

!ORLs
!Oxygen

     if (oiiRLabund .ge. 1e-20) then
       RLicfO = ((heiabund + heiiabund)/heiabund)**(2./3.) * (1+(oiiCELabund/oiiiCELabund))
       OabundRL = RLicfO * oiiRLabund
     else
       RLicfO = 1.0
       OabundRL = 0
     endif

!Nitrogen

     RLicfN = 1.0
     NabundRL = niiRLabund + niiiRLabund

!Carbon

     RLicfC = 1.0
     CabundRL = ciiRLabund + ciiiRLabund

!Neon

     if (oiiRLabund .ge. 1e-20) then
       RLicfNe = OabundRL / oiiRLabund
       NeabundRL = RLicfNe * neiiRLabund
     else
       RLicfNe = 1.0
       NeabundRL = 0.0
     endif

!finish these TODO

!Printout edited to include weighted averages of Ionic species as these are neccessary for paper tables. DJS

!carbon

        iteration_result(1)%NC_abund_CEL = NCabundCEL
        iteration_result(1)%cii_abund_CEL = ciiCELabund
        iteration_result(1)%ciii_abund_CEL = ciiiCELabund
        iteration_result(1)%civ_abund_CEL = civCELabund
        iteration_result(1)%c_icf_CEL = CELicfC
        iteration_result(1)%C_abund_CEL = CabundCEL
!nitrogen
        iteration_result(1)%Nii_abund_CEL = niiCELabund
        iteration_result(1)%Niii_abund_CEL = niiiCELabund
        iteration_result(1)%Niv_abund_CEL = nivCELabund
        iteration_result(1)%Nv_abund_CEL = nvCELabund
        iteration_result(1)%N_icf_CEL = CELicfN
        iteration_result(1)%N_abund_CEL = NabundCEL
!oxygen
        iteration_result(1)%NO_abund_CEL = NOabundCEL
        iteration_result(1)%Oii_abund_CEL = oiiCELabund
        iteration_result(1)%Oiii_abund_CEL = oiiiCELabund
        iteration_result(1)%Oiv_abund_CEL = oivCELabund
        iteration_result(1)%O_icf_CEL = CELicfO
        iteration_result(1)%O_abund_CEL = OabundCEL
!neon
        iteration_result(1)%Neii_abund_CEL = neiiIRCELabund
        iteration_result(1)%Neiii_abund_CEL = neiiiCELabund
        iteration_result(1)%Neiv_abund_CEL = neivCELabund
        iteration_result(1)%Nev_abund_CEL = nevCELabund
        iteration_result(1)%Ne_icf_CEL = CELicfNe
        iteration_result(1)%Ne_abund_CEL = NeabundCEL
!argon
        iteration_result(1)%Ariii_abund_CEL = ariiiCELabund
        iteration_result(1)%Ariv_abund_CEL = arivCELabund
        iteration_result(1)%Arv_abund_CEL = arvCELabund
        iteration_result(1)%Ar_icf_CEL = CELicfAr
        iteration_result(1)%Ar_abund_CEL = ArabundCEL
!sulphur
        iteration_result(1)%Sii_abund_CEL = siiCELabund
        iteration_result(1)%Siii_abund_CEL = siiiCELabund
        iteration_result(1)%S_icf_CEL = CELicfS
        iteration_result(1)%S_abund_CEL = SabundCEL
!chlorine
        iteration_result(1)%Cliii_abund_CEL = cliiiCELabund
        iteration_result(1)%Cl_icf_CEL = CELicfCl
        iteration_result(1)%Cl_abund_CEL = ClabundCEL

!ORLs
        iteration_result(1)%Hei_abund_ORL = Heiabund
        iteration_result(1)%Heii_abund_ORL = Heiiabund
        iteration_result(1)%He_abund_ORL = Hetotabund
        iteration_result(1)%Cii_abund_ORL = ciiRLabund
        iteration_result(1)%Ciii_abund_ORL = ciiiRLabund
        iteration_result(1)%C_icf_ORL = RLicfC
        iteration_result(1)%C_abund_ORL = CabundRL
        iteration_result(1)%Nii_abund_ORL = niiRLabund
        iteration_result(1)%Niii_abund_ORL = niiiRLabund
        iteration_result(1)%N_icf_ORL = RLicfN
        iteration_result(1)%N_abund_ORL = NabundRL
        iteration_result(1)%Oii_abund_ORL = oiiRLabund
        iteration_result(1)%O_icf_ORL = RLicfO
        iteration_result(1)%O_abund_ORL = OabundRL
        iteration_result(1)%Neii_abund_ORL = NeiiRLabund
        iteration_result(1)%Ne_icf_ORL = RLicfNe
        iteration_result(1)%Ne_abund_ORL = NeabundRL

!Strong line methods

!O R23 Pilyugin 2000
x23temp1 = ILs(get_ion("oii3726    ",ILs, Iint))%int_dered
x23temp2 = ILs(get_ion("oii3729    ",ILs, Iint))%int_dered
x23temp3 = ILs(get_ion("oiii4959   ",ILs, Iint))%int_dered
x23temp4 = ILs(get_ion("oiii5007   ",ILs, Iint))%int_dered

if (ILs(get_ion("oii3728b   ",ILs, Iint))%int_dered .gt. 0 .and. x23temp3 .gt. 0 .and. x23temp4 .gt. 0) then ! OII blended
        X23 = log10(ILs(get_ion("oii3728b   ",ILs, Iint))%int_dered/(x23temp3 + x23temp4))
elseif (x23temp1 .gt. 0 .and. x23temp2 .gt. 0 .and. x23temp3 .gt. 0 .and. x23temp4 .gt. 0) then
        X23 = log10((x23temp1+x23temp2)/(x23temp3+x23temp4))
else
        X23 = 0.
endif

if (X23 .gt. 0) then
  O_R23upper = 9.50 - (1.4 * X23)
  O_R23lower = 6.53 + (1.45 * X23)
  iteration_result(1)%O_R23_upper = O_R23upper
  iteration_result(1)%O_R23_lower = O_R23lower
endif

!O N2 Pettini + Pagel 2004

ion_no1 = get_ion("nii6584    ",ILs, Iint)
if (ILs(ion_no1)%int_dered .gt. 0 .and. H_BS(1)%int_dered .gt. 0) then
  N2 = ILs(ion_no1)%int_dered / H_BS(1)%int_dered
  O_N2 = 8.90 + (0.57 * N2)
  iteration_result(1)%O_N2 = O_N2
endif

!O O3N2 Pettini + Pagel 2004

ion_no2 = get_ion("oiii5007   ",ILs, Iint)
if (ILS(ion_no1)%int_dered .gt. 0 .and. ILs(ion_no2)%int_dered .gt. 0) then
  O3N2 = log10((ILS(ion_no2)%int_dered*H_BS(1)%int_dered)/(ILS(ion_no1)%int_dered * H_BS(2)%int_dered))
  O_O3N2 = 8.73 - (0.32*O3N2)
  iteration_result(1)%O_O3N2 = O_O3N2
endif

!O Ar3O3 Stasinska 2006

ion_no1 = get_ion("ariii7135  ",ILs, Iint)
ion_no2 = get_ion("oiii5007   ",ILs, Iint)
if (ILs(ion_no1)%int_dered .gt. 0 .and. ILs(ion_no2)%int_dered .gt. 0) then
  Ar3O3 = ILs(ion_no1)%int_dered / ILs(ion_no2)%int_dered
  O_Ar3O3 = 8.91 + (0.34*Ar3O3) + (0.27*Ar3O3**2) + (0.2*Ar3O3**3)
  iteration_result(1)%O_Ar3O3 = O_Ar3O3
endif

!O S3O3 Stasinska 2006

ion_no1 = get_ion("siii9069   ",ILs, Iint)
if (ILs(ion_no1)%int_dered .gt. 0 .and. ILs(ion_no2)%int_dered .gt. 0) then
  S3O3 = ILs(ion_no1)%int_dered / ILs(ion_no2)%int_dered
  O_S3O3 = 9.37 + (2.03*S3O3) + (1.26*S3O3**2) + (0.32*S3O3**3)
  iteration_result(1)%O_S3O3 = O_S3O3
endif

!abundance discrepancy factors

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

iteration_result(1)%adf_o2plus = adfo2plus
iteration_result(1)%adf_O = adfO

iteration_result(1)%adf_n2plus = adfn2plus
iteration_result(1)%adf_n = adfn

iteration_result(1)%adf_c2plus = adfc2plus
iteration_result(1)%adf_c = adfc

iteration_result(1)%adf_ne2plus = adfne2plus
iteration_result(1)%adf_ne = adfne


contains

        SUBROUTINE get_diag(name1, name2, lines, diag)
                TYPE(line), DIMENSION(:), INTENT(IN) :: lines
                CHARACTER*11 :: name1, name2
                INTEGER :: ion_no1, ion_no2
                DOUBLE PRECISION :: diag

                ion_no1 = get_ion(name1, ILs, Iint)
                ion_no2 = get_ion(name2, ILs, Iint)

                if((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0))then
                        diag = DBLE(ILs(ion_no1)%int_dered) / DBLE(ILs(ion_no2)%int_dered)
                else
                        diag = 0.0
                endif

        END SUBROUTINE

        SUBROUTINE get_Tdiag(name1, name2, name3, lines, factor1, factor2, ratio)
                                                        !3.47,   1.403   SIII
                TYPE(line), DIMENSION(:), INTENT(IN) :: lines
                CHARACTER*11 :: name1, name2, name3
                INTEGER :: ion_no1, ion_no2, ion_no3
                DOUBLE PRECISION :: factor1, factor2, ratio, ratio2


                ion_no1 = get_ion(name1, ILs, Iint)
                ion_no2 = get_ion(name2, ILs, Iint)
                ion_no3 = get_ion(name3, ILs, Iint)

                if(((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0)) .and. (ILs(ion_no3)%int_dered .gt. 0 ))then

                        ratio = (ILs(ion_no1)%int_dered + ILs(ion_no2)%int_dered) / ILs(ion_no3)%int_dered

                        if(name1 == "siii9069   ")then
                                ratio2 = (ILs(ion_no1)%int_dered / ILs(ion_no2)%int_dered)
                                !print*, ratio2
                                if(ratio2 > (factor2-1)*0.95 .and. ratio2 < (factor2-1)*1.05)then
                                        !PRINT*, ratio2, " ", factor2-1, " 1 ", name1
                                        ratio = (ILs(ion_no1)%int_dered + ILs(ion_no2)%int_dered) / ILs(ion_no3)%int_dered
                                else if(ratio2 < (factor2-1)*0.95)then
                                        !PRINT*, ratio2, " ", factor2-1, " 2 ", name1
                                        ratio = (factor2 * ILs(ion_no2)%int_dered) / ILs(ion_no3)%int_dered
                                        !PRINT*, (factor2 * ILs(ion_no2)%int_dered), " ", (ILs(ion_no1)%int_dered + ILs(ion_no2)%int_dered)
                                else if(ratio2 > (factor2-1)*1.05)then
                                        !PRINT*, ratio2, " ", factor2-1, " 3 ", name1
                                        ratio = (factor1 * ILs(ion_no1)%int_dered) / ILs(ion_no3)%int_dered
                                else
                                        ratio = 0.0
                                end if
                        end if

                elseif(((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .eq. 0)) .and. (ILs(ion_no3)%int_dered .gt. 0 ))then
                        ratio = (ILs(ion_no1)%int_dered * factor1) / ILs(ion_no3)%int_dered
                elseif(((ILs(ion_no1)%int_dered .eq. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0)) .and. (ILs(ion_no3)%int_dered .gt. 0 ))then
                        ratio = (ILs(ion_no2)%int_dered * factor2) / ILs(ion_no3)%int_dered
                else
                        ratio = 0.0
                endif


        END SUBROUTINE



end subroutine abundances
