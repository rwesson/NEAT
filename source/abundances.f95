subroutine abundances(linelist, listlength, iteration_result, meanextinction, calculate_extinction, ILs, diagnostic_array,iion,atomicdata,maxlevs,maxtemps, heidata, switch_he, switch_icf, H_Balmer, H_Paschen, HeI_lines, HeII_lines)
use mod_abundtypes
use mod_equib
use mod_abundIO
use mod_helium
use mod_recombination_lines
use mod_extinction
use mod_resultarrays
use mod_atomicdata
use mod_oii_diagnostics
use mod_hydrogen

        implicit none
        integer, parameter :: dp=kind(1.d0)
        integer :: counter, i, j
        integer, intent(in) :: listlength
        type(line), dimension(listlength) :: linelist, linelist_orig
        character :: switch_he !switch for helium atomic data
        character :: switch_icf !switch for which ICF to use
        type(resultarray), dimension(1) :: iteration_result
        real(kind=dp), dimension(6), intent(in) :: diagnostic_array
        real(kind=dp) :: flux_no1, flux_no2, flux_no3, flux_no4, flux_no5, flux_no6

        real(kind=dp) :: oiiNratio, oiiDens, oiiiTratio, oiiiTemp, oiiiIRNratio, oiiiIRTratio, oiiiIRtemp, oiiiUVTratio, oiiiUVtemp, oiiiIRdens, niiTratio, niiTemp, ariiiIRNratio, ariiiIRdens, arivNratio, arivDens, cliiiNratio, cliiiDens, siiNratio, siiDens, siiTratio, siiTemp, siiiIRNratio, siiiIRdens, oiiTratio, oiiTemp, neiiiTratio, neiiiIRTratio, neiiiIRNratio, neiiiIRdens, neiiiTemp, neiiiIRTemp, oitemp, citemp
        real(kind=dp) :: ciiiNratio,neivNratio,nevTratio,siiiTratio,ariiiTratio,arvTratio,lowtemp,lowdens,medtemp,ciiidens,meddens,siiitemp,ariiitemp,hightemp,neivdens,highdens,arvtemp,nevtemp,oiTratio,ciTratio
        real(kind=dp) :: oiiRLabund, niiRLabund, ciiRLabund, cii4267rlabund, neiiRLabund, ciiiRLabund, niiiRLabund, RLabundtemp, weight
        real(kind=dp) :: ciiiCELabund, niiCELabund, niiiIRCELabund, niiiUVCELabund, oiiCELabund, oiiiCELabund, oiiiIRCELabund, oivCELabund, neiiIRCELabund, neiiiIRCELabund, neiiiCELabund, neivCELabund, siiCELabund, siiiCELabund, siiiIRCELabund, sivIRCELabund, cliiCELabund, cliiiCELabund, clivCELabund, ariiiCELabund, arivCELabund, ariiiIRCELabund, nivCELabund, niiiCELabund, ciiCELabund, civCELabund, nvCELabund, nevCELabund, arvCELabund, CELabundtemp, ciCELabund, oiCELabund
        real(kind=dp) :: fn4
        real(kind=dp) :: CELicfO, CELicfC, CELicfN, CELicfNe, CELicfAr, CELicfS, CELicfCl
        real(kind=dp) :: RLicfO, RLicfC, RLicfN, RLicfNe, RLicfHe
        real(kind=dp) :: CabundRL, CabundCEL, NabundRL, NabundCEL, OabundRL, OabundCEL, NeabundRL, NeabundCEL, SabundCEL, ArabundCEL, NOabundCEL, NCabundCEL, ClabundCEL
        real(kind=dp) :: adfC, adfN, adfO, adfNe, w1, w2, w3, w4
        real(kind=dp) :: adfC2plus, adfN2plus, adfO2plus, adfNe2plus
        real(kind=dp) :: c1, c2, c3, meanextinction
        real(kind=dp) :: heiabund,heiiabund,Hetotabund, heICFfactor, OICFfactor, upsilon, upsilonprime
        real(kind=dp) :: Te_balmer, Te_paschen
        real(kind=dp) :: oii4649, oii4089,oii4662,oii_te,oii_ne
        real(kind=dp) :: ratio_5876_4471, ratio_6678_4471, te_5876_4471, te_6678_4471

        logical :: calculate_extinction

        type(line), dimension(2) :: Balmer_jump, Paschen_jump

! index arrays to locate lines of interest in the main line list

        integer, dimension(3:40) :: H_Balmer
        integer, dimension(4:39) :: H_Paschen
        integer, dimension(44) :: HeI_lines
        integer, dimension(1) :: HeII_lines
        type(cel), dimension(82) :: ILs !todo: work out why this becomes undefined on entry if its shape is assumed

!atomic data

        integer :: iion !# of ions in Ilines
        integer :: maxlevs,maxtemps
        type(atomic_data),dimension(22) :: atomicdata
        real(kind=dp), dimension(21,14,44) :: heidata

! recombination line variables

        type RLabund
           character(len=7) :: Multiplet
           real(kind=dp) :: Abundance
        end type RLabund

        type (RLabund), dimension(12) :: oiimultiplets
        type (RLabund), dimension(7) :: niimultiplets

        real(kind=dp) :: balmerdec_density, paschendec_density

! strong line variables
        real(kind=dp) :: X23,O_R23upper, O_R23lower, N2,O_N2, O3N2, O_O3N2, Ar3O3, O_Ar3O3, S3O3, O_S3O3, x23temp1, x23temp2, x23temp3, x23temp4

! recombination contribution to CEL line fluxes for both RL and CEL abundances

        real(kind=dp) :: nii5754recCEL=0.d0, oii7325recCEL=0.d0, oiii4363recCEL=0.d0, nii5754recRL=0.d0, oii7325recRL=0.d0, oiii4363recRL=0.d0

! debugging

#ifdef CO
        print *,"subroutine: abundances"
#endif

! initialise some variables

        oiiRLabund = 0.d0
        niiRLabund = 0.d0
        CELicfCl = 0.d0
        CELicfS = 0.d0
        CELicfAr = 0.d0
        CELicfNe = 0.d0
        CELicfN = 0.d0
        CELicfO = 0.d0
        Ar3O3 = 0.D0
        O_Ar3O3 = 0.D0
        oiimultiplets(:)%abundance = 0.d0
        niimultiplets(:)%abundance = 0.d0

        linelist%name="           "

        linelist_orig = linelist

        !store fluxes of blends for later retrieval

        where (linelist%blend_intensity .gt. 0.d0)
          linelist%intensity=0.d0
          linelist%int_err=0.d0
        endwhere

        !is H beta detected? if not, we can't do anything.
        !todo: the calculation of expected Hbeta from Halpha and c(Hb) needs fixing as Hb=0 now triggers an error outside the loop - 15/05/2015

        if (linelist(H_Balmer(4))%intensity .eq. 0.D0) then
          print *,"   No H beta found"
          if (linelist(H_Balmer(3))%intensity .gt. 0.D0 .and. calculate_extinction .eqv. .false.) then
            print *, "   Estimating from observed H alpha and specified c(Hb)" !do we really want this to be printed out 10000 times?
            print "(4X,F6.2,A5,F6.3,A13)",linelist(H_Balmer(3))%intensity," and ",meanextinction," respectively"
            linelist(H_Balmer(4))%intensity=(linelist(H_Balmer(3))%intensity/2.85)*10.**(meanextinction*linelist(H_Balmer(3))%flambda)
            print "(4X,A9,F6.2)","H beta = ",(linelist(H_Balmer(4))%intensity)
          else
            print *,"   Specify the extinction from the command line with the -c option to proceed"
            stop
          endif
        endif

! note that extinction is recalculated later on after diagnostics are known, so
! changes here may also need to be made in the subsequent section too

        if (calculate_extinction) then
                call calc_extinction_coeffs(linelist,H_Balmer, c1, c2, c3, meanextinction, dble(10000.),dble(1000.))

                if (meanextinction .lt. 0.0 .or. isnan(meanextinction)) then
                   meanextinction = 0.d0
                endif
! NaN can happen if the code tries to calculate the extinction but there is no
! H alpha line
! ideally an error message should be written here, but as it's inside the loop
! this would result in 10,000 error messages when uncertainties are being
! calculated
! to be improved...

        endif
        !actual dereddening

        call deredden(linelist, meanextinction)

!diagnostics
        call get_diag("ciii1909   ","ciii1907   ", ciiiNratio)        ! ciii ratio
        call get_diag("oii3729    ","oii3726    ", oiiNratio )        ! oii ratio
        call get_diag("neiv2425   ","neiv2423   ", neivNratio )       ! neiv ratio
        call get_diag("sii6731    ","sii6716    ", siiNratio )        ! s ii ratio
        call get_diag("cliii5537  ","cliii5517  ", cliiiNratio )      ! Cl iii ratio
        call get_diag("ariv4740   ","ariv4711   ", arivNratio )       ! Ar iv ratio
        call get_diag("oiii88um   ","oiii52um   ", oiiiIRNratio )     ! oiii ir ratio
        call get_diag("siii33p5um ","siii18p7um ", siiiIRNratio )       ! siii ir ratio
        call get_diag("neiii15p5um","neiii36p0um", neiiiIRNratio )       ! neiii ir ratio
        call get_diag("ariii9um   ","ariii21p8um", ariiiIRNratio )       ! ariii ir ratio

! temperature ratios:
        !For diagnostic ratios with the sum of two lines on top, the get_Tdiag
        !subroutine will properly calculate the ratio if one of the two lines is
        !missing, from the theoretical expected line strengths.

        call get_Tdiag("nii6548    ","nii6584    ","nii5754    ", "nii                 ", niiTratio)        ! N II
        call get_Tdiag("oiii5007   ","oiii4959   ","oiii4363   ", "oiii                ", oiiiTratio)        ! O III
        call get_Tdiag("neiii3868  ","neiii3967  ","neiii3342  ", "neiii               ", neiiiTratio)        ! Ne III
        call get_Tdiag("neiii3868  ","neiii3967  ","neiii15p5um", "neiii               ", neiiiIRTratio)! Ne III ir
        call get_Tdiag("nev3426    ","nev3345    ","nev2975    ", "nev                 ", nevTratio)        !!ne v
        call get_Tdiag("siii9069   ","siii9531   ","siii6312   ", "siii                ", siiiTratio)        !s iii
        call get_Tdiag("ariii7135  ","ariii7751  ","ariii5192  ", "ariii               ", ariiiTratio)        !ar iii
        call get_Tdiag("arv6435    ","arv7005    ","arv4625    ", "arv                 ", arvTratio)        !ar v
        call get_Tdiag("ci9850     ","ci9824     ","ci8727     ", "ci                  ", ciTratio)      !C I
        call get_Tdiag("oi6364     ","oi6300     ","oi5577     ", "oi                  ", oiTratio)      !O I
        call get_Tdiag("oiii5007   ","oiii4959   ","oiii52um   ", "oiii                ", oiiiIRTratio) ! OIII ir
        call get_Tdiag("oiii5007   ","oiii4959   ","oiii1666   ", "oiii                ", oiiiUVTratio) ! OIII UV

        !Fixed, DJS

! O II


        if(get_cel_flux("oii7319b   ",linelist,ILs) .gt. 0.d0) then

                flux_no1 = get_cel_flux("oii7319b    ",linelist,ILs)
                flux_no2 = get_cel_flux("oii7330b    ",linelist,ILs)
                flux_no3 = get_cel_flux("oii3726    ",linelist,ILs)
                flux_no4 = get_cel_flux("oii3729    ",linelist,ILs)

                if (flux_no1 .gt. 0 .and. flux_no2 .gt. 0 .and. flux_no3 .gt. 0 .and. flux_no4 .gt. 0) then
                        oiiTratio = (flux_no1+flux_no2)/(flux_no3+flux_no4)
                else
                        oiiTratio = 0.d0
                endif

       elseif(get_cel_flux("oii7319    ",linelist,ILs) .gt. 0)then

                flux_no1 = get_cel_flux("oii7319    ",linelist,ILs)
                flux_no2 = get_cel_flux("oii7320    ",linelist,ILs)
                flux_no3 = get_cel_flux("oii7330    ",linelist,ILs)
                flux_no4 = get_cel_flux("oii7331    ",linelist,ILs)
                flux_no5 = get_cel_flux("oii3726    ",linelist,ILs)
                flux_no6 = get_cel_flux("oii3729    ",linelist,ILs)

                if (flux_no1 .gt. 0 .and. flux_no2 .gt. 0 .and. flux_no3 .gt. 0 .and. flux_no4 .gt. 0 .and. flux_no5 .gt. 0 .and. flux_no6 .gt. 0) then
                        oiiTratio = ((flux_no1+flux_no2)+(flux_no3+flux_no4))/(flux_no5+flux_no6)
                else
                        oiiTratio = 0.d0
                endif
        !add condition for 3727 blend


       else
                       oiiTratio=0.d0
       endif

! S II
      flux_no1 = get_cel_flux("sii6716    ",linelist,ILs)
      flux_no2 = get_cel_flux("sii6731    ",linelist,ILs)
      flux_no3 = get_cel_flux("sii4068    ",linelist,ILs)
      flux_no4 = get_cel_flux("sii4076    ",linelist,ILs)

      if (flux_no1 .gt. 0 .and. flux_no2 .gt. 0 .and. flux_no3 .gt. 0 .and. flux_no4 .gt. 0) then
           siiTratio = (flux_no1+flux_no2)/(flux_no3+flux_no4)
      elseif (flux_no1 .gt. 0 .and. flux_no2 .gt. 0 .and. flux_no3 .gt. 0 .and. flux_no4 .eq. 0.d0) then
! if 4076 is not seen, assume that 4076/4068 = 0.338.  This is not an exact value but for 5000<Te<15000 and 0<ne<1e5 it varies by less than 5%
           siiTratio = (flux_no1+flux_no2)/(1.338*flux_no3)
      else
           siiTratio = 0.d0
      endif

! now get diagnostics zone by zone.

! low ionisation
        ! Edited to stop high limits being included in diagnostic averages. DJS
      lowtemp = 10000.d0

      do i = 1,2

        oiiDens=0.d0
        siiDens=0.d0
        counter=0

        call get_diagnostic("oii       ","1,2/                ","1,3/                ",oiiNratio,"D",lowtemp, oiiDens,maxlevs,maxtemps,atomicdata,iion)
        call get_diagnostic("sii       ","1,2/                ","1,3/                ",siiNratio,"D",lowtemp, siiDens,maxlevs,maxtemps,atomicdata,iion)

        if (oiidens .gt. 0.d0) counter=counter+1
        if (siidens .gt. 0.d0) counter=counter+1

        if (counter .eq. 0 .and. diagnostic_array(1) .eq. 0) then
          lowdens = 1000.d0
        elseif (diagnostic_array(1) .gt. 0.0) then
          lowdens = diagnostic_array(1)
        else
          lowdens = (oiiDens + siiDens) / counter
        endif

        counter=0

        call get_diagnostic("oii       ","2,4,2,5,3,4,3,5/    ","1,2,1,3/            ",oiiTratio,"T",lowdens,oiiTemp,maxlevs,maxtemps,atomicdata,iion)
        call get_diagnostic("sii       ","1,2,1,3/            ","1,4,1,5/            ",siiTratio,"T",lowdens,siiTemp,maxlevs,maxtemps,atomicdata,iion)
        call get_diagnostic("nii       ","2,4,3,4/            ","4,5/                ",niiTratio,"T",lowdens,niitemp,maxlevs,maxtemps,atomicdata,iion)
        call get_diagnostic("ci        ","2,4,3,4/            ","4,5/                ",ciTratio,"T",lowdens,citemp,maxlevs,maxtemps,atomicdata,iion)
        call get_diagnostic("oi        ","1,4,2,4/            ","4,5/                ",oiTratio,"T",lowdens,oitemp,maxlevs,maxtemps,atomicdata,iion)

        if (oiiTemp .gt. 0.d0) counter=counter+1
        if (siiTemp .gt. 0.d0) counter=counter+1
        if (niiTemp .gt. 0.d0) counter=counter+5
        if (ciTemp .gt. 0.d0) counter=counter+1
        if (oiTemp .gt. 0.d0) counter=counter+1

        if (counter .gt. 0 .and. diagnostic_array(4) .eq. 0) then
          lowtemp = ((5*niitemp) + siitemp + oiitemp + oitemp + citemp) / counter
        elseif (diagnostic_array(4) .gt. 0.0) then
          lowtemp = diagnostic_array(4)
        else
          lowtemp = 10000.d0
        endif

      enddo

! medium ionisation
      cliiiDens = 0.D0
      ciiiDens = 0.D0
      arivDens = 0.D0
      oiiiIRdens = 0.D0
      ariiiIRdens = 0.D0
      siiiIRdens = 0.D0
      neiiiIRdens = 0.D0

      medtemp = lowtemp

      do i = 1,2

        counter = 0
        call get_diagnostic("cliii     ","1,2/                ","1,3/                ",cliiiNratio,"D",medtemp, cliiiDens,maxlevs,maxtemps,atomicdata,iion)
        call get_diagnostic("ciii      ","1,3/                ","1,4/                ",ciiiNratio,"D",medtemp, ciiiDens,maxlevs,maxtemps,atomicdata,iion)
        call get_diagnostic("ariv      ","1,2/                ","1,3/                ",arivNratio,"D",medtemp, arivDens,maxlevs,maxtemps,atomicdata,iion)
! IR densities, not included in average
!Ar, S, Ne, O
        call get_diagnostic("oiii      ","1,2/                ","2,3/                ",oiiiIRNratio,"D",medtemp, oiiiIRDens,maxlevs,maxtemps,atomicdata,iion)
        call get_diagnostic("ariii     ","1,2/                ","2,3/                ",ariiiIRNratio,"D",medtemp, ariiiIRDens,maxlevs,maxtemps,atomicdata,iion)
        call get_diagnostic("siii      ","1,2/                ","2,3/                ",siiiIRNratio,"D",medtemp, siiiIRDens,maxlevs,maxtemps,atomicdata,iion)
        call get_diagnostic("neiii     ","1,2/                ","2,3/                ",neiiiIRNratio,"D",medtemp, neiiiIRDens,maxlevs,maxtemps,atomicdata,iion)

        if (cliiiDens .gt. 0) counter=counter+1
        if (ciiiDens .gt. 0) counter=counter+1
        if (arivDens .gt. 0) counter=counter+1

        if (counter .eq. 0 .and. diagnostic_array(2) .eq. 0) then
          meddens = lowdens
        elseif (diagnostic_array(2) .gt. 0.0) then
          meddens = diagnostic_array(2)
        else
          meddens = (ciiiDens + cliiiDens + arivDens) / counter
        endif

        if (meddens .le. 1.d0) meddens = 1.d0 !todo: output warning

        counter = 0

        call get_diagnostic("oiii      ","2,4,3,4/            ","4,5/                ",oiiiTratio,"T",meddens,oiiiTemp,maxlevs,maxtemps,atomicdata,iion)
        call get_diagnostic("siii      ","2,4,3,4/            ","4,5/                ",siiiTratio,"T",meddens,siiiTemp,maxlevs,maxtemps,atomicdata,iion)
        call get_diagnostic("ariii     ","1,4,2,4/            ","4,5/                ",ariiiTratio,"T",meddens,ariiitemp,maxlevs,maxtemps,atomicdata,iion)
        call get_diagnostic("neiii     ","1,4,2,4/            ","4,5/                ",neiiiTratio,"T",meddens,neiiitemp,maxlevs,maxtemps,atomicdata,iion)
!extra IR temperatures, not included in the average at the moment but could be if we decide that would be good
        call get_diagnostic("neiii     ","1,4,2,4/            ","1,2/                ",neiiiIRTratio,"T",meddens,neiiiIRtemp,maxlevs,maxtemps,atomicdata,iion)
        call get_diagnostic("oiii      ","2,4,3,4/            ","2,3/                ",oiiiIRTratio,"T",oiiiIRdens,oiiiIRtemp,maxlevs,maxtemps,atomicdata,iion)
!oiii UV temperature, also not included in the average at the moment
        call get_diagnostic("oiii      ","2,4,3,4/            ","3,6/                ",oiiiUVTratio,"T",oiidens,oiiiUVtemp,maxlevs,maxtemps,atomicdata,iion)

!averaging

       if (oiiiTemp .gt. 0.d0) counter=counter+4
       if (siiiTemp .gt. 0.d0) counter=counter+1
       if (ariiiTemp .gt. 0.d0) counter=counter+2
       if (neiiiTemp .gt. 0.d0) counter=counter+2

       if (counter .gt. 0 .and. diagnostic_array(5) .eq. 0.0) then
         medtemp = (4*oiiitemp + siiitemp + 2*ariiitemp + 2*neiiitemp) / counter
       elseif (diagnostic_array(5) .gt. 0.0) then
         medtemp = diagnostic_array(5)
       else
         medtemp = lowtemp
       endif

!if there is no lowtemp, set it to medtemp
!this would break if lowtemp was actually calculated to be 10000.d0 but I think the probability of that is low enough that we can live without checking for it

        if (lowtemp .eq. 10000.d0 .and. medtemp .ne. 10000.d0) lowtemp = medtemp

        !dereddening again

        if (calculate_extinction) then

        !update extinction. DS 22/10/11
          meanextinction=0
          call calc_extinction_coeffs(linelist,H_Balmer, c1, c2, c3, meanextinction, medtemp, lowdens)

          if (meanextinction .lt. 0.0 .or. isnan(meanextinction)) then
             meanextinction = 0.d0
          endif

          linelist = linelist_orig
          call deredden(linelist, meanextinction)

        endif ! end of extinction calculating

      enddo ! end of diagnostic iteration

        iteration_result(1)%mean_cHb = meanextinction


! high ionisation
      hightemp = medtemp

      do i = 1,2

        call get_diagnostic("neiv      ","1,2/                ","1,3/                ",neivNratio,"D",hightemp, neivDens,maxlevs,maxtemps,atomicdata,iion)

        if (diagnostic_array(3) .gt. 0.0) then
          highdens = diagnostic_array(3)
        elseif (neivdens .le. 0d0) then
          neivdens = 0d0
          highdens = meddens
        else
          highdens = neivDens
        endif

        if (highdens .lt. 1.d0) highdens = 1.d0 !todo: output warning

        counter = 0

        call get_diagnostic("arv       ","2,4,3,4/            ","4,5/                ",arvTratio,"T",highdens,arvTemp,maxlevs,maxtemps,atomicdata,iion)
        call get_diagnostic("nev       ","2,4,3,4/            ","4,5/                ",nevTratio,"T",highdens,nevtemp,maxlevs,maxtemps,atomicdata,iion)

        if (arvTemp .gt. 0.d0) counter=counter+1
        if (nevTemp .gt. 0.d0) counter=counter+1

        if (counter .gt. 0 .and. diagnostic_array(6) .eq. 0) then
          hightemp = (arvtemp + nevtemp) / counter
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
if(niitemp .gt. 0.2)then
        iteration_result(1)%NII_temp = niitemp
else if(INT(niitemp) .eq. -1)then
        iteration_result(1)%NII_temp = 35000.
endif
iteration_result(1)%NII_temp_ratio = niiTratio

if(oiitemp .gt. 0.2)then
        iteration_result(1)%OII_temp = oiitemp
else if(INT(oiitemp) .eq. -1)then
        iteration_result(1)%OII_temp = 20000
endif
iteration_result(1)%OII_temp_ratio = oiiTratio

if(siitemp .gt. 0.2 )then
        iteration_result(1)%SII_temp = siitemp
else if(INT(siitemp) .eq. -1)then
        iteration_result(1)%SII_temp = 35000.
endif
iteration_result(1)%SII_temp_ratio = siiTratio

if(oitemp .gt. 0.2 )then
        iteration_result(1)%OI_temp = oitemp
else if(INT(oitemp) .eq. -1)then
        iteration_result(1)%OI_temp = 35000.
endif
iteration_result(1)%OI_temp_ratio = oiTratio

if(citemp .gt. 0.2 )then
        iteration_result(1)%CI_temp = citemp
else if(INT(citemp) .eq. -1)then
        iteration_result(1)%CI_temp = 35000.
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
if(oiiitemp .gt. 0.2)then
        iteration_result(1)%OIII_temp = oiiitemp
else if(INT(oiiitemp) .eq. -1)then
        iteration_result(1)%OIII_temp = 35000.
endif
iteration_result(1)%OIII_temp_ratio = oiiiTratio

if(neiiitemp .gt. 0.2)then
        iteration_result(1)%NeIII_temp = neiiitemp
else if(INT(neiiitemp) .eq. -1)then
        iteration_result(1)%NeIII_temp = 35000.
endif
iteration_result(1)%NeIII_temp_ratio = neiiiTratio

if(ariiitemp .gt. 0.2)then
        iteration_result(1)%ArIII_temp = ariiitemp
else if(INT(ariiitemp) .eq. -1)then
        iteration_result(1)%ArIII_temp = 35000.
endif
iteration_result(1)%ArIII_temp_ratio = AriiiTratio

if(siiitemp .gt. 0.2)then
        iteration_result(1)%SIII_temp = siiitemp
else if(int(siiitemp) .eq. -1)then
        iteration_result(1)%SIII_temp = 20000
endif
iteration_result(1)%SIII_temp_ratio = siiiTratio

if(oiiiIRtemp .gt. 0.2)then
        iteration_result(1)%OIII_IR_temp = oiiiIRtemp
else if(int(oiiiIRtemp) .eq. -1)then
        iteration_result(1)%OIII_IR_temp = 35000.
endif
iteration_result(1)%OIII_IR_temp_ratio = oiiiIRTratio

if(oiiiUVtemp .gt. 0.2)then
        iteration_result(1)%OIII_UV_temp = oiiiUVtemp
else if(int(oiiiUVtemp) .eq. -1)then
        iteration_result(1)%OIII_UV_temp = 35000.
endif
iteration_result(1)%OIII_UV_temp_ratio = oiiiUVTratio

if(neiiiIRtemp .gt. 0.2)then
        iteration_result(1)%NeIII_IR_temp = neiiiIRtemp
else if(int(neiiiIRtemp) .eq. -1)then
        iteration_result(1)%NeIII_IR_temp = 35000.
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
! He II

        if (Heii_lines(1) .gt. 0) then
          call get_heii_abund(medtemp,meddens,linelist(Heii_lines(1))%int_dered,heiiabund)
          linelist(Heii_lines(1))%abundance = heiiabund
        else
          heiiabund = 0.d0
        endif

! He I

        if (switch_he .eq. "S") then
          call get_hei_smits_new(linelist,medtemp,meddens,HeI_lines,heidata, heiabund)
        else
          call get_hei_porter(linelist,medtemp,meddens,HeI_lines,heidata, heiabund)
        endif

        hetotabund = heiabund + heiiabund

!  Calculate T_e from Balmer Jump using equation 3 of Liu et al. (2001)

        Te_balmer = 0.d0

        Balmer_jump(1)%int_dered = 0.d0
        Balmer_jump(2)%int_dered = 0.d0

        do i=1,listlength
            if(linelist(i)%wavelength .eq. 3645.50) Balmer_jump(1) = linelist(i)
            if(linelist(i)%wavelength .eq. 3646.50) Balmer_jump(2) = linelist(i)
        enddo

        if(Balmer_jump(1)%int_dered .gt. 0 .and. Balmer_jump(2)%int_dered .gt. 0 .and. H_Balmer(11) .gt. 0) then
        Te_balmer = (Balmer_jump(1)%int_dered - Balmer_jump(2)%int_dered)/linelist(H_Balmer(11))%int_dered
        Te_balmer = Te_balmer**(-1.5)
        Te_balmer = Te_balmer*368
        Te_balmer = Te_balmer*(1+0.256*heiabund+3.409*heiiabund)
        endif

        iteration_result%Balmer_jump_temp = Te_balmer

!and from Paschen jump, using equation 7 of Fang et al. (2011)

        Te_paschen = 0.d0
        Paschen_jump(1)%int_dered = 0.d0
        Paschen_jump(2)%int_dered = 0.d0

        do i=1,listlength
            if(linelist(i)%wavelength .eq. 8100.d0) Paschen_jump(1) = linelist(i)
            if(linelist(i)%wavelength .eq. 8400.d0) Paschen_jump(2) = linelist(i)
        enddo

        if(Paschen_jump(1)%int_dered .gt. 0 .and. Paschen_jump(2)%int_dered .gt. 0 .and. H_Paschen(11) .gt. 0) then
        Te_paschen = (Paschen_jump(1)%int_dered - Paschen_jump(2)%int_dered)/linelist(H_Paschen(11))%int_dered
        Te_paschen = Te_paschen**(-1.77)
        Te_paschen = Te_paschen*8.72
        Te_paschen = Te_paschen*(1+0.52*heiabund+4.40*heiiabund)
        endif

        iteration_result%Paschen_jump_temp = Te_paschen

! calculate n_e from Balmer decrement

        call balmer_densities(linelist,H_Balmer,medtemp,balmerdec_density)
        iteration_result%balmerdec_density=balmerdec_density

! and from Paschen decrement. todo: add to iteration result

        call paschen_densities(linelist,H_Paschen,medtemp,paschendec_density)
        iteration_result%paschendec_density=paschendec_density

! get Te from He line ratios
! equations derived from Smits 1996 data, at ne=5000.
! to be replaced with proper interpolation in density and using whichever dataset was selected
! Te(5876/4471) = a*r**4 + b*r**3 + c*r**2 + d*r + e + f/(r-g)
! a=6858; b=-99452.0; c=541944.8; d=-1317918.7; e=1210016.1; f=23.21; g=2.833
! accurate to +-50K
! Te(6678/4471) = a*r**2+b*r+c (r<0.775)
!                 a*r**4+b**r*3+c**r*2+d**r+e (r>0.775)

        ratio_5876_4471=0.d0
        ratio_6678_4471=0.d0
        Te_5876_4471=0.d0
        Te_6678_4471=0.d0

        if (linelist(HeI_lines(10))%int_dered .gt. 0.d0 .and. linelist(HeI_lines(15))%int_dered .gt.0.d0) then
          ratio_5876_4471=linelist(HeI_lines(15))%int_dered/linelist(HeI_lines(10))%int_dered
          if (ratio_5876_4471.lt. 4.2 .and. ratio_5876_4471.gt. 2.84) then
            Te_5876_4471 = 6858*ratio_5876_4471**4 - 99452.0*ratio_5876_4471**3 + 541944.8*ratio_5876_4471**2 - 1317918.7*ratio_5876_4471 + 1210016.1 + 23.21/(ratio_5876_4471-2.833)
          endif
        endif

        if (linelist(HeI_lines(10))%int_dered .gt. 0.d0 .and. linelist(HeI_lines(16))%int_dered .gt.0.d0) then
          ratio_6678_4471=linelist(HeI_lines(16))%int_dered/linelist(HeI_lines(10))%int_dered
          if (ratio_6678_4471 .lt. 0.775 .and. ratio_6678_4471 .gt. 0.59) then
            Te_6678_4471 = -64152.4*ratio_6678_4471**2 + 33203.3*ratio_6678_4471 + 22553.1
          elseif (ratio_6678_4471 .ge. 0.775 .and. ratio_6678_4471 .lt. 1.2) then
            Te_6678_4471 = 1.06363e+06*ratio_6678_4471**4 - 4.4416e+06*ratio_6678_4471**3 + 6.96343e+06*ratio_6678_4471**2 - 4.86749e+06*ratio_6678_4471 + 1.28347e+06
          endif
        endif

        iteration_result%ratio_5876_4471=ratio_5876_4471
        iteration_result%ratio_6678_4471=ratio_6678_4471
        iteration_result%Te_5876_4471=Te_5876_4471
        iteration_result%Te_6678_4471=Te_6678_4471

! get te and ne from OII RL diagnostics if we have them

        oii4089=0.d0
        oii4649=0.d0
        oii4662=0.d0

        do i=1,listlength
          if (abs(linelist(i)%wavelength-4089.29) .lt. 0.005) oii4089=linelist(i)%int_dered
          if (abs(linelist(i)%wavelength-4649.13) .lt. 0.005) oii4649=linelist(i)%int_dered
          if (abs(linelist(i)%wavelength-4661.63) .lt. 0.005) oii4662=linelist(i)%int_dered
        enddo

        if (oii4089.gt.0 .and. oii4649.gt.0 .and. oii4662.gt.0) then
          call oii_rl_diagnostics(oii4649/oii4089,oii4649/oii4662,oii_te,oii_ne)
          iteration_result%oii_te=oii_te
          iteration_result%oii_ne=oii_ne
        else
          iteration_result%oii_te=0.d0
          iteration_result%oii_ne=0.d0
        endif

! get abundances for all CELs

        do i = 1,size(ILs)
          if (ILs(i)%location .gt. 0) then !the line is in the spectrum
             linelist(ILs(i)%location)%name=ILs(i)%ion
             if (ILs(i)%zone .eq. "low ") then
               call get_abundance(ILs(i)%ion, ILs(i)%transition, lowtemp, lowdens,linelist(ILs(i)%location)%int_dered, linelist(ILs(i)%location)%abundance,maxlevs,maxtemps,atomicdata,iion)
             elseif (ILs(i)%zone .eq. "med ") then
               call get_abundance(ILs(i)%ion, ILs(i)%transition, medtemp, meddens,linelist(ILs(i)%location)%int_dered, linelist(ILs(i)%location)%abundance,maxlevs,maxtemps,atomicdata,iion)
             elseif (ILs(i)%zone .eq. "high") then
               call get_abundance(ILs(i)%ion, ILs(i)%transition, hightemp, highdens,linelist(ILs(i)%location)%int_dered, linelist(ILs(i)%location)%abundance,maxlevs,maxtemps,atomicdata,iion)
             endif
             if ((linelist(ILs(i)%location)%abundance .ge. 1e-10) .and. (linelist(ILs(i)%location)%abundance .lt. 10 ) ) then
             elseif ((linelist(ILs(i)%location)%abundance .ge. 10 ) .or. (linelist(ILs(i)%location)%abundance .lt. 1E-10 ) )then
               linelist(ILs(i)%location)%abundance = 0
             endif
!print *,ILs(i)%location,get_cel_flux(ILs(i)%name,linelist,ILs),linelist(ILs(i)%location)%wavelength,linelist(ILs(i)%location)%abundance
          endif
        enddo

! calculate averages

        celabundtemp = 0.d0
                niiCELabund = 0.d0
        weight = 0.d0
        do i=get_ion("nii6548    ", ILs), get_ion("nii6584    ", ILs)
          if (linelist(ILs(i)%location)%abundance .ge. 1e-10) niiCELabund = niiCELabund + linelist(ILs(i)%location)%abundance*linelist(ILs(i)%location)%intensity
          if (linelist(ILs(i)%location)%abundance .ge. 1e-10) weight = weight + linelist(ILs(i)%location)%intensity
        enddo
        if (weight .ge. 1e-20) then
          niiCELabund = niiCELabund / weight
        else
          niiCELabund = 0.d0
        endif

        niiiIRCELabund = get_cel_abundance("niii57um   ",linelist,ILs)
        niiiUVCELabund = get_cel_abundance("niii1751   ",linelist,ILs)

        if (niiiIRCELabund .ge. 1e-20 .and. niiiUVCELabund .ge. 1e-20) then
          niiiCELabund = (niiiIRCELabund + niiiUVCELabund)/2
        elseif (niiiIRCELabund .ge. 1e-20) then
          niiiCELabund = niiiIRCELabund
        elseif (niiiUVCELabund .ge. 1e-20) then
          niiiCELabund = niiiUVCELabund
        else
          niiiCELabund = 0.d0
        endif

        nivCELabund = 0.d0

        do i= get_ion("niv1483    ", ILs), get_ion("niv1485b   ", ILs) ! would screw up if blend and non blends were both given
          if (ILs(i)%location .gt. 0) nivCELabund = nivCELabund + linelist(ILs(i)%location)%abundance*linelist(ILs(i)%location)%intensity
          if (ILs(i)%location .gt. 0) weight = weight + linelist(ILs(i)%location)%intensity
        enddo
        if (weight .ge. 1e-20) then
          nivCELabund = nivCELabund / weight
        else
          nivCELabund = 0.d0
        endif

        nvCELabund = linelist(get_ion("nv1240     ",ILs))%abundance

        !OII CEL routine fixed. DJS 23/11/10

        if(get_cel_flux("oii3728b   ",linelist,ILs) .gt. 0.0)then
                !calc abundance from doublet blend
                oiiCELabund = get_cel_abundance("oii3728b   ",linelist,ILs)

        else if(get_cel_flux("oii3729    ",linelist,ILs) .gt. 0.0 .and. get_cel_flux("oii3726    ",linelist,ILs) .gt. 0.0 )then
                !calc abundance from doublet
                w1 = linelist(get_ion("oii3729    ", ILs))%intensity
                w2 = linelist(get_ion("oii3726    ", ILs))%intensity

                oiiCELabund = (w1*get_cel_abundance("oii3729    ",linelist,ILs) + w2*get_cel_abundance("oii3726    ",linelist,ILs))/(w1+w2)

        else if((get_cel_flux("oii3728b   ",linelist,ILs) .eq. 0.0 .and. (get_cel_flux("oii3729    ",linelist,ILs) .eq. 0.0 .and. get_cel_flux("oii3726    ",linelist,ILs) .eq. 0.0 )) .and.  (get_cel_abundance("oii7330b   ",linelist,ILs) .gt. 0.0 .or. get_cel_abundance("oii7319b   ",linelist,ILs) .gt. 0.0)  )then
                !calc abundance based on far red blends
                w1 = linelist(get_ion("oii7330b   ", ILs))%intensity
                w2 = linelist(get_ion("oii7319b   ", ILs))%intensity

                oiiCELabund = (w1*get_cel_abundance("oii7330b   ",linelist,ILs) + w2*get_cel_abundance("oii7319b   ",linelist,ILs))/(w1+w2)

        else if ((get_cel_flux("oii3728b   ",linelist,ILs) .eq. 0.0 .and. (get_cel_flux("oii3729    ",linelist,ILs) .eq. 0.0 .and. get_cel_flux("oii3726    ",linelist,ILs) .eq. 0.0 )) .and. ( (get_cel_abundance("oii7320    ",linelist,ILs) .gt. 0.0 .or. get_cel_abundance("oii7319    ",linelist,ILs) .gt. 0.0) .or.  (get_cel_abundance("oii7330    ",linelist,ILs) .gt. 0.0 .or. get_cel_abundance("oii7331    ",linelist,ILs) .gt. 0.0)) )then
                !calc abundance based on far red quadruplet
                if(linelist(get_ion("oii7319    ", ILs))%intensity .gt. 0) w1 = linelist(get_ion("oii7319    ", ILs))%intensity
                if(linelist(get_ion("oii7320    ", ILs))%intensity .gt. 0) w2 = linelist(get_ion("oii7320    ", ILs))%intensity
                if(linelist(get_ion("oii7330    ", ILs))%intensity .gt. 0) w3 = linelist(get_ion("oii7330    ", ILs))%intensity
                if(linelist(get_ion("oii7331    ", ILs))%intensity .gt. 0) w4 = linelist(get_ion("oii7320    ", ILs))%intensity

                !if statements stop non existent lines being granted infinite weight 1/(0/0)^2 = infinity, defaults to zero if no line detected which keeps the following calculation honest

                oiiCELabund = (w1*get_cel_abundance("oii7319    ",linelist,ILs) + w2*get_cel_abundance("oii7320    ",linelist,ILs) + w3*get_cel_abundance("oii7330    ",linelist,ILs) + w4*get_cel_abundance("oii7331    ",linelist,ILs))/(w1+w2+w3+w4)

        else
                oiiCELabund = 0.d0
        endif

        celabundtemp = 0.d0
                oiiiCELabund = 0.d0
        weight = 0.d0
        do i=get_ion("oiii4959   ", ILs), get_ion("oiii5007   ", ILs)
          if (ILs(i)%location .gt. 0) oiiiCELabund = oiiiCELabund + linelist(ILs(i)%location)%abundance *linelist(ILs(i)%location)%intensity
          if (ILs(i)%location .gt. 0) weight = weight + linelist(ILs(i)%location)%intensity
        enddo
        if (weight .ge. 1e-20) then
          oiiiCELabund = oiiiCELabund / weight
        else
          oiiiCELabund = 0.d0
        endif

        celabundtemp = 0.d0
                oiiiIRCELabund = 0.d0
        weight = 0.d0
        do i=get_ion("oiii52um   ", ILs), get_ion("oiii88um   ", ILs)
          if (ILs(i)%location .gt. 0) oiiiIRCELabund = oiiiIRCELabund + linelist(ILs(i)%location)%abundance *linelist(ILs(i)%location)%intensity
          if (ILs(i)%location .gt. 0) weight = weight + linelist(ILs(i)%location)%intensity
        enddo
        if (weight .ge. 1e-20) then
          oiiiIRCELabund = oiiiIRCELabund / weight
        else
          oiiiIRCELabund = 0.d0
        endif

        oivCELabund = get_cel_abundance("oiv25p9um  ",linelist,ILs)

        neiiIRCELabund = linelist(get_ion("neii12p8um ", ILs)  )%abundance
        neiiiIRCELabund = linelist(get_ion("neiii15p5um ", ILs)  )%abundance

        celabundtemp = 0.d0
                neiiiCELabund = 0.d0
        weight = 0.d0
        do i=get_ion("neiii3868  ", ILs), get_ion("neiii3967  ", ILs)
          if (ILs(i)%location .gt. 0) neiiiCELabund = neiiiCELabund + linelist(ILs(i)%location)%abundance*linelist(ILs(i)%location)%intensity
          if (ILs(i)%location .gt. 0) weight = weight + linelist(ILs(i)%location)%intensity
        enddo
        if (weight .ge. 1e-20) then
          neiiiCELabund = neiiiCELabund / weight
        else
          neiiiCELabund = 0.d0
        endif


        celabundtemp = 0.d0
                neivCELabund = 0.d0
        weight = 0.d0
        do i=get_ion("neiv2423   ", ILs), get_ion("neiv4725b  ", ILs) ! would screw up if blends and non blends were given
          if (ILs(i)%location .gt. 0) neivCELabund = neivCELabund + linelist(ILs(i)%location)%abundance*linelist(ILs(i)%location)%intensity
          if (ILs(i)%location .gt. 0) weight = weight + linelist(ILs(i)%location)%intensity
        enddo
        if (weight .ge. 1e-20) then
          neivCELabund = neivCELabund / weight
        else
          neivCELabund = 0.d0
        endif

        celabundtemp = 0.d0
                siiCELabund = 0.d0
        weight = 0.d0
        do i=get_ion("sii4068    ", ILs), get_ion("sii6731    ", ILs)
          if (ILs(i)%location .gt. 0) siiCELabund = siiCELabund + linelist(ILs(i)%location)%abundance*linelist(ILs(i)%location)%intensity
          if (ILs(i)%location .gt. 0) weight = weight + linelist(ILs(i)%location)%intensity
        enddo
        if (weight .ge. 1e-20) then
          siiCELabund = siiCELabund / weight
        else
          siiCELabund = 0.d0
        endif

        !siiiCELabund = 0 ! ILs(28)%abundance
        !celabundtemp = 0.
        !weight = 0.
        !do i=get_ion("siii9069   ", ILs), get_ion("siii9531   ", ILs)
        !  if (linelist(ILs(i)%location)%abundance .gt. 0) siiiCELabund = siiiCELabund + linelist(ILs(i)%location)%abundance/  ((ILs(i)%int_err/ linelist(ILs(i)%location)%intensity)  **2)
        !  if (linelist(ILs(i)%location)%abundance .gt. 0) weight = weight + 1/  ((ILs(i)%int_err/ linelist(ILs(i)%location)%intensity)  **2)
        !enddo
             !if (weight .gt. 0) then
        !  siiiCELabund = siiiCELabund / weight
        !else
        !  siiiCELabund = 0.d0
        !endif

        ! SIII abundance section previously buggered. SIII 9069 and 9531 doublet special due to absorption lines. See Liu, Barlow etc 1995 (Far Red IR Lines in a PN), one of 9069/9531 is ALWAYS absorbed however the other is then always fine.  We can tell which is is by looking at the ratio 9069/9531.
        if (get_cel_abundance("siii9069   ",linelist,ILs) .eq. 0 .and. get_cel_abundance("siii9531   ",linelist,ILs) .eq. 0) then
                siiiCELabund = get_cel_abundance("siii6312   ",linelist,ILs)
        else !only calculate all the IR SIII telluric absorption if the lines are present.
                siiiCELabund = get_cel_flux("siii9069   ",linelist,ILs) / get_cel_flux("siii9531   ",linelist,ILs)
        !I am using siiiCELabund as a switch for the following if statement, replace this with another variable if you like but I don't think it matters.. (DJS)

                if(siiiCELabund .lt. (1.05 * 0.403) .and. siiiCELabund .gt. (0.90 * 0.403))then !this case should never occur

                        if(linelist(get_ion("siii9069   ", ILs))%intensity .gt. 0) w1 = linelist(get_ion("siii9069   ", ILs))%intensity
                        if(linelist(get_ion("siii9531   ", ILs))%intensity .gt. 0) w2 = linelist(get_ion("siii9531   ", ILs))%intensity

                        siiiCELabund= ( w1*get_cel_abundance("siii9069   ",linelist,ILs) + w2*get_cel_abundance("siii9531   ",linelist,ILs) )/(w1+w2)

                elseif(siiiCELabund .gt. (1.1 * 0.403) )then
                        !9531 absorbed
                        siiiCELabund = linelist(get_ion("siii9069   ", ILs) )%abundance

                elseif(siiiCELabund .lt. (0.9 * 0.403) )then
                        !9069 absorbed
                        siiiCELabund = linelist(get_ion("siii9531   ", ILs) )%abundance
                else
                        siiiCELabund=0.d0
                endif
        endif

        siiiIRCELabund = linelist(get_ion("siii18p7um ", ILs)   )%abundance
        sivIRCELabund = linelist(get_ion("siv10p5um  ", ILs) )%abundance

        celabundtemp = 0.d0
                cliiiCELabund = 0.d0
        weight = 0.d0
        do i=get_ion("cliii5517  ", ILs), get_ion("cliii5537  ", ILs)
          if (ILs(i)%location .gt. 0) cliiiCELabund = cliiiCELabund + linelist(ILs(i)%location)%abundance*linelist(ILs(i)%location)%intensity
          if (ILs(i)%location .gt. 0) then
            weight = weight + linelist(ILs(i)%location)%intensity
          endif
        enddo
        if (weight .ge. 1e-20) then
          cliiiCELabund = cliiiCELabund / weight
        else
          cliiiCELabund = 0.d0
        endif

        celabundtemp = 0.d0
                ariiiCELabund = 0.d0
        weight = 0.d0
        do i=get_ion("ariii7135  ", ILs), get_ion("ariii7751  ", ILs)
          if (ILs(i)%location .gt. 0) ariiiCELabund = ariiiCELabund + linelist(ILs(i)%location)%abundance*linelist(ILs(i)%location)%intensity
          if (ILs(i)%location .gt. 0) then
            weight = weight + linelist(ILs(i)%location)%intensity
          endif
        enddo
        if (weight .ge. 1e-20) then
          ariiiCELabund = ariiiCELabund / weight
        else
          ariiiCELabund = 0.d0
        endif

        celabundtemp = 0.d0
                arivCELabund = 0.d0
        weight = 0.d0
        do i=get_ion("ariv4711   ", ILs), get_ion("ariv4740   ", ILs)
        arivCELabund = arivCELabund + linelist(ILs(i)%location)%abundance*linelist(ILs(i)%location)%intensity
          if (ILs(i)%location .gt. 0) then
            weight = weight + linelist(ILs(i)%location)%intensity
          endif
        enddo
        if (weight .ge. 1e-20) then
          arivCELabund = arivCELabund / weight
        else
          arivCELabund = 0.d0
        endif

        ariiiIRCELabund = get_cel_abundance("ariii9um   ",linelist,ILs)

        ciiCELabund = get_cel_abundance("cii2325    ",linelist,ILs)
        civCELabund = get_cel_abundance("civ1548    ",linelist,ILs)

        celabundtemp = 0.d0
                ciiiCELabund = 0.d0
        weight = 0.d0
        do i=get_ion("ciii1907   ", ILs), get_ion("ciii1909b  ", ILs) ! would screw up if blend and non-blend were given.
          if (ILs(i)%location .gt. 0) then
            ciiiCELabund = ciiiCELabund + linelist(ILs(i)%location)%abundance*linelist(ILs(i)%location)%intensity
            weight = weight + linelist(ILs(i)%location)%intensity
          endif
        enddo
        if (weight .ge. 1e-20) then
          ciiiCELabund = ciiiCELabund / weight
        else
          ciiiCELabund = 0.d0
        endif

        celabundtemp = 0.d0
                neivCELabund = 0.d0
        weight = 0.d0
        do i=get_ion("neiv2423   ", ILs), get_ion("neiv2424b  ", ILs)! would break if blend and individual lines were both specified
          if (ILs(i)%location .gt. 0) neivCELabund = neivCELabund + linelist(ILs(i)%location)%abundance*linelist(ILs(i)%location)%intensity
          if (ILs(i)%location .gt. 0) then
            weight = weight + linelist(ILs(i)%location)%intensity
          endif
        enddo
        if (weight .ge. 1e-20) then
          neivCELabund = neivCELabund / weight
        else
          neivCELabund = 0.d0
        endif

        celabundtemp = 0.d0
                nevCELabund = 0.d0
        weight = 0.d0
        do i=get_ion("nev3345    ", ILs), get_ion("nev3426    ", ILs)
          if (ILs(i)%location .gt. 0) nevCELabund = nevCELabund + linelist(ILs(i)%location)%abundance*linelist(ILs(i)%location)%intensity
          if (ILs(i)%location .gt. 0) then
            weight = weight + linelist(ILs(i)%location)%intensity
          endif
        enddo
        if (weight .ge. 1e-20) then
          nevCELabund = nevCELabund / weight
        else
          nevCELabund = 0.d0
        endif

        celabundtemp = 0.d0
                arvCELabund = 0.d0
        weight = 0.d0
        do i=get_ion("arv6435    ", ILs), get_ion("arv7005    ", ILs)
          if (ILs(i)%location .gt. 0) arvCELabund = arvCELabund + linelist(ILs(i)%location)%abundance*linelist(ILs(i)%location)%intensity
          if (ILs(i)%location .gt. 0) then
            weight = weight + linelist(ILs(i)%location)%intensity
          endif
        enddo
        if (weight .ge. 1e-20) then
          arvCELabund = arvCELabund / weight
        else
          arvCELabund = 0.d0
        endif

         celabundtemp = 0.d0
                 ciCELabund = 0.d0
        weight = 0.d0
        do i=get_ion("ci9824     ", ILs), get_ion("ci9850     ", ILs)
          if (ILs(i)%location .gt. 0) ciCELabund = ciCELabund + linelist(ILs(i)%location)%abundance*linelist(ILs(i)%location)%intensity
          if (ILs(i)%location .gt. 0) then
            weight = weight + linelist(ILs(i)%location)%intensity
          endif
        enddo
        if (weight .ge. 1e-20) then
          ciCELabund = ciCELabund / weight
        else
          ciCELabund = 0.d0
        endif

        NCabundCEL = ciCELabund

         celabundtemp = 0.d0
                 oiCELabund = 0.d0
        weight = 0.d0
        do i=get_ion("oi6300     ", ILs), get_ion("oi6364     ", ILs)
          if (ILs(i)%location .gt. 0) then
            oiCELabund = oiCELabund + linelist(ILs(i)%location)%abundance*linelist(ILs(i)%location)%intensity
            weight = weight + linelist(ILs(i)%location)%intensity
          endif
        enddo
        if (weight .ge. 1e-20) then
          oiCELabund = oiCELabund / weight
        else
          oiCELabund = 0.d0
        endif

        NOabundCEL = oiCELabund

! now get abundances for ORLs
! o2+
       call oii_rec_lines(medtemp,meddens,1.d0,oiiRLs)

       do i = 1,listlength
         do j = 1,size(oiiRLs)
          if (abs(linelist(i)%wavelength-oiiRLs(j)%Wave) .le. 0.005) then
            oiiRLs(j)%Obs = linelist(i)%int_dered
            if (oiiRLs(j)%Int .eq. 0.0) then
              oiiRLs(j)%abundance = 0.D0
            else
              oiiRLs(j)%abundance = oiiRLs(j)%obs/oiiRLs(j)%Int
            endif
            linelist(i)%abundance = oiiRLs(j)%abundance
            linelist(i)%name="O II       "
            linelist(i)%latextext="O~{\sc ii}     "
          endif
         enddo
       enddo

!N2+

       call nii_rec_lines(medtemp,meddens,1.d0,niiRLs)

       do i = 1,listlength
         do j = 1,size(niiRLs)
          if (abs(linelist(i)%wavelength-niiRLs(j)%Wave) .le. 0.005) then
            niiRLs(j)%Obs = linelist(i)%int_dered
            niiRLs(j)%abundance = niiRLs(j)%obs/niiRLs(j)%Int
            linelist(i)%abundance = niiRLs(j)%abundance
            linelist(i)%name="N II       "
            linelist(i)%latextext="N~{\sc ii}     "
          endif
         enddo
       enddo

!C2+
       call cii_rec_lines(medtemp,meddens,1.d0,ciiRLs)

       do i = 1,listlength
         do j = 1,size(ciiRLs)
          if (abs(linelist(i)%wavelength-ciiRLs(j)%Wave) .le. 0.005) then
            ciiRLs(j)%Obs = linelist(i)%int_dered
            ciiRLs(j)%abundance = ciiRLs(j)%obs/ciiRLs(j)%Int
            linelist(i)%abundance = ciiRLs(j)%abundance
            linelist(i)%name="C II       "
            linelist(i)%latextext="C~{\sc ii}     "
          endif
         enddo
       enddo

!Ne2+
       call neii_rec_lines(medtemp,meddens,1.d0,neiiRLs)

       do i = 1,listlength
         do j = 1,size(neiiRLs)
          if (abs(linelist(i)%wavelength-neiiRLs(j)%Wave) .le. 0.005) then
            neiiRLs(j)%Obs = linelist(i)%int_dered
            neiiRLs(j)%abundance = neiiRLs(j)%obs/neiiRLs(j)%Int
            linelist(i)%abundance = neiiRLs(j)%abundance
            linelist(i)%name="Ne II      "
            linelist(i)%latextext="Ne~{\sc ii}    "
          endif
         enddo
       enddo

!C3+, N3+
       call xiii_rec_lines(medtemp,meddens,1.d0,xiiiRLs)

       do i = 1,listlength
         do j = 1,size(xiiiRLs)
          if (abs(linelist(i)%wavelength-xiiiRLs(j)%Wave) .le. 0.005) then
            xiiiRLs(j)%Obs = linelist(i)%int_dered
            xiiiRLs(j)%abundance = xiiiRLs(j)%obs/xiiiRLs(j)%Int
            linelist(i)%abundance = xiiiRLs(j)%abundance
            if (j .le. 4) then
              linelist(i)%name="C III      "
              linelist(i)%latextext="C~{\sc iii}    "
            else
              linelist(i)%name="N III      "
              linelist(i)%latextext="N~{\sc iii}    "
            endif
          endif
         enddo
       enddo

      rlabundtemp = 0.d0
      weight = 0.d0

!cii recombination lines

      do i = 1,size(ciiRLs)
        if (ciiRLs(i)%wave .eq. 4267.15D0) then
          cii4267rlabund = ciiRLs(i)%abundance
        endif
      enddo

      ciirlabund = sum(ciiRLs%obs)/sum(ciiRLs%int, mask=ciiRLs%obs.gt.0.d0)

!nii recombination lines

  if (sum(niiRLs%int) .gt. 0) then
!      print "(ES9.2)",rlabundtemp / weight

      niimultiplets%Multiplet = (/"V3     ","V5     ","V8     " ,"V12    ","V20    ","V28    ","3d-4f  "/)

! get multiplet abundances from coadded intensity

      do j = 1,6
        rlabundtemp = 0.d0
        weight = 0.d0
        do i = 1,size(niiRLs)
          if (niiRLs(i)%Mult .eq. niimultiplets(j)%Multiplet .and. niiRLs(i)%obs .gt. 0) then
!            rlabundtemp = rlabundtemp + (niiRLs(i)%obs * niiRLs(i)%abundance)
!            weight = weight + niiRLs(i)%obs
             rlabundtemp = rlabundtemp + niiRLs(i)%obs
             weight = weight + niiRLs(i)%Int
          endif
        enddo
      if (isnan((rlabundtemp/weight))) then
        niimultiplets(j)%Abundance = 0.d0
      else
        niimultiplets(j)%Abundance = rlabundtemp/weight
      endif
      enddo

!3d-4f transitions

      rlabundtemp = 0.d0
      weight = 0.d0
      do i = 1,size(niiRLs)
        if (niiRLs(i)%obs .gt. 0 .and. (niiRLs(i)%Mult(4:4) .eq. "a" .or. niiRLs(i)%Mult(4:4).eq."b")) then
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

      rlabundtemp = 0.d0
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

  if (sum(oiiRLs%Int).gt. 0) then

!     print "(ES9.2)",rlabundtemp/weight

      oiimultiplets%Multiplet = (/" V1    "," V2    "," V5    " ," V10   "," V11   "," V12   "," V19   "," V20   "," V25   "," V28   "," V33   "," 3d-4f "/)

! get multiplet abundances from coadded intensity

      do j = 1,11 !10 XXX
        rlabundtemp = 0.d0
        weight = 0.d0
        do i = 1,size(oiiRLs)
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
          oiimultiplets(j)%Abundance = 0.d0
        endif
      enddo

      rlabundtemp = 0.d0
      weight = 0.d0

      do i=1,size(oiiRLs)
        if (oiiRLs(i)%Term1(4:5) .eq. "3d" .and. oiiRLs(i)%Term2(3:4) .eq. "4f" .and. oiiRLs(i)%Mult .ne. "       " .and. oiiRLs(i)%obs .gt. 0) then
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

      rlabundtemp = 0.d0
      weight = 0
      do i = 1,12
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

      rlabundtemp = 0.d0
      weight = 0.d0

      if (sum(neiiRLs%Int) .gt. 0) then
        neiiRLabund = sum(neiiRLs%obs)/sum(neiiRLs%Int, mask=neiiRLs%obs.gt.0.d0)
      else
        neiiRLabund = 0.D0
      endif

!N3+ and C3+

      rlabundtemp=0.d0
      weight=0.d0

      do i=1,4
        rlabundtemp = rlabundtemp + xiiiRLs(i)%obs
        weight = weight + xiiiRLs(i)%Int
      enddo
      if (weight .gt. 0) then
        ciiiRLabund = rlabundtemp/weight
      else
        ciiiRLabund = 0.d0
      endif

      do i=5,6
        rlabundtemp = rlabundtemp + xiiiRLs(i)%obs
        weight = weight + xiiiRLs(i)%Int
      enddo
      if (weight .gt. 0) then
        niiiRLabund = rlabundtemp/weight
      else
        niiiRLabund = 0.d0
      endif

! calculate recombination contributions to CELs using equations 1-3 from Liu et al. 2000
! equations give recombination flux relative to Hb=1

      if (linelist(get_ion("nii5754    ", ILs))%Int_dered .gt. 0.d0) then
        if (niiicelabund .gt. 0.d0) then
          nii5754recCEL=10000.*(3.19*(medtemp/1.e4)**0.30*niiicelabund)/linelist(get_ion("nii5754    ", ILs))%Int_dered
        endif
        if (niirlabund .gt. 0.d0) then
          nii5754recRL=10000.*(3.19*(medtemp/1.e4)**0.30*niirlabund)/linelist(get_ion("nii5754    ", ILs))%Int_dered
        endif
      endif

      if (linelist(get_ion("oii7320    ", ILs))%Int_dered .gt. 0 .and. linelist(get_ion("oii7330    ", ILs))%Int_dered .gt. 0) then
        if (oiiicelabund .gt. 0) then
          oii7325recCEL=10000.*(9.36*(medtemp/1.e4)**0.44*oiiicelabund)/(linelist(get_ion("oii7320    ", ILs))%Int_dered+linelist(get_ion("oii7330    ", ILs))%Int_dered)
        endif
        if (oiiRLabund .gt. 0) then
          oii7325recRL=10000.*(9.36*(medtemp/1.e4)**0.44*oiirlabund)/(linelist(get_ion("oii7320    ", ILs))%Int_dered+linelist(get_ion("oii7330    ", ILs))%Int_dered)
        endif
      endif

      if (linelist(get_ion("oiii4363   ", ILs))%Int_dered .gt. 0 .and. hetotabund .gt. 0) then
        oiii4363recCEL=10000.*12.4*(hightemp/1.e4)**0.59*(((hetotabund/heiabund)**0.66666)-1)*(oiicelabund+oiiicelabund)/linelist(get_ion("oiii4363   ", ILs))%Int_dered
        oiii4363recRL=10000.*12.4*(hightemp/1.e4)**0.59*(((hetotabund/heiabund)**0.66666)-1)*(oiiRLabund)/linelist(get_ion("oiii4363   ", ILs))%Int_dered
      endif

! ICFs

! first calculate the helium factor that appears in the ORL ICFs and several KB94 ICFS
! if no he lines are seen set this factor to zero, so that abundances relying on
! this ratio are not calculated.

     if (heiabund .gt. 0.) then
        heICFfactor = (heiabund + heiiabund)/heiabund
     elseif (heiabund .eq.0. .and. heiiabund .eq.0.0) then
        heICFfactor = 0.d0
     endif

! now apply the scheme chosen by the user

if (switch_icf .eq. "K") then

! ICFs (Kingsburgh + Barlow 1994)

! oxygen - complete
     OabundCEL = 0.d0
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
                CELicfO = heICFfactor**(2./3.) !KB94 A9
                OabundCEL = CELicfO * (oiiCELabund + oiiiCELabund) !A10
        endif

! nitrogen - complete
     NabundCEL = 0.d0
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
     CabundCEL = 0.d0
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
        CELicfC = heICFfactor**(1./3.) !A20
        CabundCEL = CELicfC * (ciiCELabund + ciiiCELabund + civCELabund) !A21
     elseif (ciiCELabund .lt. 1e-20 .and. ciiiCELabund .ge. 1e-20 .and. civCELabund .ge. 1e-20) then !C+ not seen
        CELicfC = 1/(1 - (niiCELabund/NabundCEL) - 2.7*(nvCELabund/NabundCEL)) !A22
        if (CELicfC .gt. 5) then !if ICF is greater than 5
           CELicfC = (niiCELabund + niiiCELabund + nivCELabund + nvCELabund)/(niiiCELabund + nivCELabund) !A23
        endif
        CabundCEL = CELicfC * (ciiiCELabund+ civCELabund) !A24
     elseif (ciiCELabund .ge. 1e-20 .and. ciiiCELabund .ge. 1e-20 .and. civCELabund .ge. 1e-20 .and. heiiabund .ge. 1e-20 .and. (nivCELabund .lt. 1e-20 .or. nvCELabund .lt. 1e-20)) then !final case i think.
        CELicfC = ((oiiCELabund + oiiiCELabund)/oiiiCELabund)*heICFfactor**(1./3.) !A25
        CabundCEL = CELicfC * (ciiCELabund + ciiiCELabund + civCELabund) !A26
     else
        CabundCEL = 0.D0 ! nothing at all seen
        CELicfC = 1.0
     endif

! Neon - complete
     NeabundCEL = 0.d0
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
     ArabundCEL = 0.d0
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
     SabundCEL = 0.d0
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
     ClabundCEL = 0.d0
    if (cliiiCELabund .ge. 1e-20 .and. siiiCELabund .ge. 1e-20) then
       CELicfCl = SabundCEL/siiiCELabund
       ClabundCEL = CELicfCl * cliiiCELabund
    endif

elseif (switch_icf .eq. "P") then
 !PTPR92 ICF - equation numbers in the paper are given in brackets

!CELs
!Carbon - no ICF specified
     CabundCEL = 0.d0
     CELicfC = 0.d0
!Chlorine - no ICF specified
     ClabundCEL = 0.d0
     CELicfCl = 0.d0
!Oxygen
     OabundCEL = 0.d0
     OabundCEL = oiiCELabund + oiiiCELabund ! (13)

!Nitrogen
     NabundCEL = 0.d0
     if (oiiCELabund .gt. 1e-20) then
       CELicfN = OabundCEL/oiiCELabund
       NabundCEL = niiCELabund * CELicfN ! (14)
     endif

!Neon
     NeabundCEL = 0.d0
     if (oiiiCELabund .gt. 1e-20) then
       CELicfNe = (oiiCELabund + oiiiCELabund)/oiiiCELabund
       NeabundCEL = CELicfNe * NeiiiCELabund ! (15)
     endif

!Sulphur - no ICF given in paper, see discussion on P168
     SabundCEL = SiiCELabund + SiiiCELabund

!Argon - no ICF given in paper, see discussion on P168
     ArabundCEL = AriiiCELabund + ArivCELabund

!Iron - TODO - to be added once we've added iron ionic abundance calculations
!     if (niiCELabund .gt. 1.e-20) then
!       CELicfFe = NabundCEL/niiCELabund
!       FeabundCEL = CELicfFe * FeiiCELabund ! (16)
!     endif

!Helium
     RLicfHe = 1.0
     if (SabundCEL - SiiCELabund .gt. 1.e-20) then
       RLicfHe = 1+(SiiCELabund/(SabundCEL-SiiCELabund))
       Hetotabund = RLicfHe * Heiabund ! (21) - note that this icf is not meaningful if only S+ is seen
     endif

elseif (switch_icf .eq. "D") then
!Delgado-Inglada 2014

!starting definitions

!Equation 4:
  if (heiabund+heiiabund .gt. 0.D0) then
    upsilon = heiiabund/(heiabund+heiiabund)
  else
    upsilon=0.d0
  endif
!Equation 5:
  if (oiiCELabund + oiiiCELabund .gt. 0.D0) then
    OICFfactor = oiiiCELabund/(oiiCELabund + oiiiCELabund)
  else
    oICFfactor = 0.d0
  endif

!helium - no ICF for neutral so it's already been calculated earlier in this
!routine

!oxygen
!equation 12:

  if (oiiCELabund .gt. 0.D0 .and. oiiiCELabund .gt. 0.D0) then
    celICFO = 10.**((0.08*upsilon + 0.006*upsilon**2)/(0.34-0.27*upsilon)) !should check if upsilon is greater than 0.95 before using
    oabundCEL = celICFO * (oiiCELabund + oiiiCELabund)
  else
    celICFO = 1.0
    oabundCEL = 0.D0
  endif

!nitrogen
!equation 14:

  if (niiCELabund .gt. 0.D0 .and. oiiCELabund .gt. 0.D0 .and. upsilon .gt. 0.D0) then
    CELicfN = 10.**(-0.16*oICFfactor*(1+log10(upsilon)))
    nabundCEL = CELicfN * niiCELabund * OabundCEL / oiiCELabund
  elseif (niiCELabund .gt. 0.D0 .and. oiiCELabund .gt. 0.D0 .and. upsilon .eq. 0.D0) then
    CELicfN = 0.64*Oicffactor
    nabundCEL = CELicfN * niiCELabund * OabundCEL / oiiCELabund
  else
    CELicfN = 1.0
    nabundCEL = 0.D0
  endif

!neon

  if (neiiiCELabund .gt. 0.D0 .and. oiiiCELabund .gt. 0.D0 .and. nevCELabund .eq. 0.D0) then
    if (heiiabund .eq. 0.D0) then
      upsilonprime = 0.01
    elseif (heiiabund .gt. 0.D0 .and. upsilon .lt. 0.015) then
      upsilonprime = 0.015
    else
      upsilonprime = upsilon
    endif
!equation 17:
    CELicfNe = oICFfactor + ((0.014/upsilonprime + 2*upsilonprime**2.7)**3 * (0.7 + 0.2*oICFfactor-0.8*oICFfactor**2))
    NeabundCEL = CELicfne * OabundCEL * (neiiiCELabund / oiiiCELabund)
!equation 20:
  elseif (neiiiCELabund .gt. 0.D0 .and. nevCELabund .gt. 0.D0) then
    CELicfNe = (1.31+12.68*upsilon**2.57)**0.27
    NeabundCEL = CELicfNe * OabundCEL * (neiiiCELabund + nevCELabund)
  else
    NeabundCEL = 0.D0
    CELicfNe = 1.0
  endif

!sulphur

  if (siiiCELabund .eq. 0.D0 .and. siiCELabund .gt. 0.D0) then !equation 23:
    CELicfS = 10.**(0.31-0.51*upsilon)
    SabundCEL = CELicfS * OabundCEL * (siiCELabund/oiiCELabund)
  elseif (siiiCELabund .gt. 0.D0 .and. siiCELabund .gt. 0.D0) then !equation 26:
    CELicfS = 10.**((-0.02 - 0.03*OICFfactor - 2.31*OICFfactor**2 + 2.19*OICFfactor**3) / (90.69 + 2.09*OICFfactor - 2.69*OICFfactor**2))
    SabundCEL = CELicfS * OabundCEL * (siiCELabund+siiiCELabund)/oiiCELabund
  else
    CELicfS = 1.0
    SabundCEL = 0.D0
  endif

!chlorine
!we currently don't calculate clii or cliv so they are set to zero
cliiCELabund = 0.d0
clivCELabund = 0.d0

  if (cliiCELabund .gt. 0.D0 .and. cliiiCELabund .gt. 0.D0 .and. clivCELabund .gt. 0.D0) then !equation 32:
    CELicfCl = 0.98+(0.56-0.57*upsilon)**7.64
    ClabundCEL = CELicfCl * (cliiCELabund + cliiiCELabund + clivCELabund)
  elseif (OICFfactor .le. 0.02 .and. cliiCELabund .gt. 0.D0 .and. cliiiCELabund .gt. 0.D0) then
    CELicfCl = 1.0
    ClabundCEL = cliiCELabund + cliiiCELabund
  elseif (OICFfactor .gt. 0.02 .and. OICFfactor .lt. 0.95 .and. cliiiCELabund .gt. 0.D0 .and. oiiCELabund .gt. 0.D0) then !equation 29:
    CELicfCL = (4.1620 - (4.1622*OICFfactor**0.21))**0.75
    ClabundCEL = CELicfCL * (cliiiCELabund / oiiCELabund)
  else
    CELicfCl = 1.0
    CLabundCEL = 0.D0
  endif

!argon

  if ((oiiCELabund+oiiiCELabund) .gt. 0.D0 .and. ariiiCELabund .gt. 0.D0 .and. OICFfactor .le. 0.5) then !equation 35
    CELicfAr = 10.**(0.05/(0.06+OICFfactor) - 0.07)
    ArabundCEL = CELicfAr * OabundCEL * ariiiCELabund/(oiiCELabund+oiiiCELabund)
  elseif ((oiiCELabund+oiiiCELabund) .gt. 0.D0 .and. ariiiCELabund .gt. 0.D0 .and. OICFfactor .gt. 0.5) then !equation 36
    CELicfAr = 10.**(0.03/(0.4-0.3*OICFfactor) - 0.05)
    ArabundCEL = CELicfAr * OabundCEL * ariiiCELabund/(oiiCELabund+oiiiCELabund)
  else
    CELicfAr = 1.0
    ArabundCEL = 0.D0
  endif

!carbon
!none given for CELs, D-I 2014 is for optical lines only

    celICFC = 1.0
    cabundCEL = 0.D0

endif

!ORLs - not from KB94 but based on Liu et al. 2000
!redefine He ICF factor, we use the same variable name for different things in DI14 and KB94
!todo: fix this ugly workaround
     if (heiabund .gt. 0.) then
        heICFfactor = (heiabund + heiiabund)/heiabund
     elseif (heiabund .eq.0. .and. heiiabund .eq.0.0) then
        heICFfactor = 0.d0
     endif
!Oxygen

     if (oiiRLabund .ge. 1e-20) then
       RLicfO = heICFfactor**(2./3.) * (1+(oiiCELabund/oiiiCELabund))
       OabundRL = RLicfO * oiiRLabund
     else
       RLicfO = 1.0
       OabundRL = 0
     endif

!Nitrogen

     RLicfN = 1.0
     NabundRL = niiRLabund + niiiRLabund

!Carbon

     if (switch_icf .eq. "D") then !equation 39
       RLicfC = 0.05 + 2.21*OICFfactor - 2.77*OICFfactor**2 + 1.74*OICFfactor**3
       CabundRL = RLicfC * ciiRLabund / oiiRLabund
     else
       RLicfC = 1.0
       CabundRL = ciiRLabund + ciiiRLabund
     endif

!Neon

     if (oiiRLabund .ge. 1e-20) then
       RLicfNe = OabundRL / oiiRLabund
       NeabundRL = RLicfNe * neiiRLabund
     else
       RLicfNe = 1.0
       NeabundRL = 0.d0
     endif

!rewrite all the ICFs so that they are simply the factor by which the summed ionic abundances are multipled to get the total elemental abundance

if (ciiCELabund + ciiiCELabund + civCELabund .gt. 0.d0) CELicfC = CabundCEL/(ciiCELabund + ciiiCELabund + civCELabund)
if (niiCELabund + niiiCELabund + nivCELabund + nvCELabund .gt. 0.d0) CELicfN = NabundCEL/(niiCELabund + niiiCELabund + nivCELabund + nvCELabund)
if (oiiCELabund + oiiiCELabund + oiiiIRCELabund + oivCELabund .gt. 0.d0) CELicfO = OabundCEL/(oiiCELabund + oiiiCELabund + oiiiIRCELabund + oivCELabund)
if (neiiiCELabund + neivCELabund + nevCELabund .gt. 0.d0) CELicfNe = NeabundCEL/(neiiiCELabund + neivCELabund + nevCELabund)
if (siiCELabund + siiiCELabund + siiiIRCELabund + sivIRCELabund .gt. 0.d0) CELicfS = SabundCEL/(siiCELabund + siiiCELabund + siiiIRCELabund + sivIRCELabund)
if (cliiCELabund + cliiiCELabund + clivCELabund .gt. 0.d0) CELicfCl = ClabundCEL/(cliiCELabund + cliiiCELabund + clivCELabund)
if (ariiiCELabund + arivCELabund + arvCELabund .gt. 0.d0) CELicfAr = ArabundCEL/(ariiiCELabund + arivCELabund + arvCELabund)

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

!recombination contributions

        iteration_result(1)%nii5754recCEL = nii5754recCEL
        iteration_result(1)%oii7325recCEL = oii7325recCEL
        iteration_result(1)%oiii4363recCEL = oiii4363recCEL
        iteration_result(1)%nii5754recRL = nii5754recRL
        iteration_result(1)%oii7325recRL = oii7325recRL
        iteration_result(1)%oiii4363recRL = oiii4363recRL

!multiplet abundances

        iteration_result(1)%oii_v1_abund_orl = oiimultiplets(1)%abundance
        iteration_result(1)%oii_v2_abund_orl = oiimultiplets(2)%abundance
        iteration_result(1)%oii_v5_abund_orl = oiimultiplets(3)%abundance
        iteration_result(1)%oii_v10_abund_orl = oiimultiplets(4)%abundance
        iteration_result(1)%oii_v11_abund_orl = oiimultiplets(5)%abundance
        iteration_result(1)%oii_v12_abund_orl = oiimultiplets(6)%abundance
        iteration_result(1)%oii_v19_abund_orl = oiimultiplets(7)%abundance
        iteration_result(1)%oii_v20_abund_orl = oiimultiplets(8)%abundance
        iteration_result(1)%oii_v25_abund_orl = oiimultiplets(9)%abundance
        iteration_result(1)%oii_v28_abund_orl = oiimultiplets(10)%abundance
        iteration_result(1)%oii_v33_abund_orl = oiimultiplets(11)%abundance
        iteration_result(1)%oii_3d4f_abund_orl = oiimultiplets(12)%abundance
        iteration_result(1)%nii_v3_abund_orl = niimultiplets(1)%abundance
        iteration_result(1)%nii_v5_abund_orl = niimultiplets(2)%abundance
        iteration_result(1)%nii_v8_abund_orl = niimultiplets(3)%abundance
        iteration_result(1)%nii_v12_abund_orl = niimultiplets(4)%abundance
        iteration_result(1)%nii_v20_abund_orl = niimultiplets(5)%abundance
        iteration_result(1)%nii_v28_abund_orl = niimultiplets(6)%abundance
        iteration_result(1)%nii_3d4f_abund_orl = niimultiplets(7)%abundance

!Strong line methods

!O R23 Pilyugin 2000
x23temp1 = get_cel_flux("oii3726    ",linelist,ILs)
x23temp2 = get_cel_flux("oii3729    ",linelist,ILs)
x23temp3 = get_cel_flux("oiii4959   ",linelist,ILs)
x23temp4 = get_cel_flux("oiii5007   ",linelist,ILs)

if (get_cel_flux("oii3728b   ",linelist,ILs) .gt. 0 .and. x23temp3 .gt. 0 .and. x23temp4 .gt. 0) then ! OII blended
        X23 = log10(get_cel_flux("oii3728b   ",linelist,ILs)/(x23temp3 + x23temp4))
elseif (x23temp1 .gt. 0 .and. x23temp2 .gt. 0 .and. x23temp3 .gt. 0 .and. x23temp4 .gt. 0) then
        X23 = log10((x23temp1+x23temp2)/(x23temp3+x23temp4))
else
        X23 = 0.d0
endif

if (X23 .gt. 0) then
  O_R23upper = 9.50 - (1.4 * X23)
  O_R23lower = 6.53 + (1.45 * X23)
  iteration_result(1)%O_R23_upper = O_R23upper
  iteration_result(1)%O_R23_lower = O_R23lower
endif

!O N2 Pettini + Pagel 2004

flux_no1 = get_cel_flux("nii6584    ",linelist,ILs)
if (flux_no1 .gt. 0 .and. linelist(H_Balmer(3))%int_dered .gt. 0) then
  N2 = flux_no1 / linelist(H_Balmer(3))%int_dered
  O_N2 = 8.90 + (0.57 * N2)
  iteration_result(1)%O_N2 = O_N2
endif

!O O3N2 Pettini + Pagel 2004

flux_no2 = get_cel_flux("oiii5007   ",linelist,ILs)
if (flux_no1 .gt. 0 .and. flux_no2 .gt. 0 .and. linelist(H_Balmer(3))%int_dered .gt. 0) then
  O3N2 = log10((flux_no2*linelist(H_Balmer(3))%int_dered)/(flux_no1 * linelist(H_Balmer(3))%int_dered))
  O_O3N2 = 8.73 - (0.32*O3N2)
  iteration_result(1)%O_O3N2 = O_O3N2
endif

!O Ar3O3 Stasinska 2006

flux_no1 = get_cel_flux("ariii7135  ",linelist,ILs)
flux_no2 = get_cel_flux("oiii5007   ",linelist,ILs)
if (flux_no1 .gt. 0 .and. flux_no2 .gt. 0) then
  Ar3O3 = flux_no1 / flux_no2
  O_Ar3O3 = 8.91 + (0.34*Ar3O3) + (0.27*Ar3O3**2) + (0.2*Ar3O3**3)
  iteration_result(1)%O_Ar3O3 = O_Ar3O3
endif

!O S3O3 Stasinska 2006

flux_no1 = get_cel_flux("siii9069   ",linelist,ILs)
if (flux_no1 .gt. 0 .and. flux_no2 .gt. 0) then
  S3O3 = flux_no1 / flux_no2
  O_S3O3 = 9.37 + (2.03*S3O3) + (1.26*S3O3**2) + (0.32*S3O3**3)
  iteration_result(1)%O_S3O3 = O_S3O3
endif

!abundance discrepancy factors

  if (oiiiCELabund .gt. 0) then
    adfO2plus = oiiRLabund/oiiiCELabund
  else
    adfO2plus = 0.d0
  endif

  if (oabundCEL .gt. 0) then
    adfO = OabundRL/OabundCEL
  else
    adfO = 0.d0
  endif


  if (ciiiCELabund .gt. 0) then
    adfC2plus = ciiRLabund/ciiiCELabund
  else
    adfC2plus = 0.d0
  endif

  if (CabundCEL .gt. 0) then
    adfC = CabundRL/CabundCEL
  else
    adfC = 0.d0
  endif


  if (NiiiCELabund .gt. 0) then
    adfN2plus = NiiRLabund/NiiiCELabund
  else
    adfN2plus = 0.d0
  endif

  if (NabundCEL .gt. 0) then
    adfN = NabundRL/NabundCEL
  else
    adfN = 0.d0
  endif


  if (NeiiiCELabund .gt. 0) then
    adfNe2plus = NeiiRLabund/NeiiiCELabund
  else
    adfNe2plus = 0.d0
  endif

  if (NeabundCEL .gt. 0) then
    adfNe = NeabundRL/NeabundCEL
  else
    adfNe = 0.d0
  endif

iteration_result(1)%adf_o2plus = adfo2plus
iteration_result(1)%adf_O = adfO

iteration_result(1)%adf_n2plus = adfn2plus
iteration_result(1)%adf_n = adfn

iteration_result(1)%adf_c2plus = adfc2plus
iteration_result(1)%adf_c = adfc

iteration_result(1)%adf_ne2plus = adfne2plus
iteration_result(1)%adf_ne = adfne

!copy CEL abundances back into main linelist array

do i=1,size(ILs) !todo: remove this if necessary
  if (ILs(i)%location .gt. 0) then
    linelist(ILs(i)%location)%abundance = linelist(ILs(i)%location)%abundance
    linelist(ILs(i)%location)%latextext = ILs(i)%latextext
  endif
enddo

!copy blend fluxes back into their original location

where (linelist%blend_intensity .gt. 0.d0)
  linelist%intensity=linelist%blend_intensity
  linelist%int_dered=linelist%blend_int_dered
  linelist%int_err=linelist%blend_int_err
endwhere

!write ion names to array

do i=3,40
  if (H_Balmer(i) .gt. 0) then
    linelist(H_Balmer(i))%name = "H I"
    linelist(H_Balmer(i))%latextext = "H~{\sc i}"
  endif
enddo
do i=4,39
  if (H_Paschen(i) .gt. 0) then
    linelist(H_Paschen(i))%name = "H I"
    linelist(H_Paschen(i))%latextext = "H~{\sc i}"
  endif
enddo

do i=1,44
  if (HeI_lines(i) .gt. 0) then
    linelist(HeI_lines(i))%name = "He I"
    linelist(HeI_lines(i))%latextext = "He~{\sc i}"
  endif
enddo

if (HeII_lines(1) .gt. 0) then
    linelist(HeII_lines(1))%name = "He II"
    linelist(HeII_lines(1))%latextext = "He~{\sc ii}"
endif

contains

subroutine get_diag(name1, name2, diag)
        implicit none
        integer, parameter :: dp = kind(1.d0)

        character(len=11) :: name1, name2
        real(kind=dp) :: flux_no1, flux_no2
        real(kind=dp) :: diag

!debugging
#ifdef CO
        print *,"subroutine: get_diag. lines=",name1,name2
#endif

         flux_no1 = get_cel_flux(name1, linelist, ILs)
         flux_no2 = get_cel_flux(name2, linelist, ILs)

         if (flux_no1 .gt. 0 .and. flux_no2 .gt. 0) then
           diag = flux_no1 / flux_no2
         else
           diag = 0.d0
         endif

end subroutine get_diag

subroutine get_Tdiag(name1, name2, name3, ion, ratio)
        !this routine gets the ratio for nebular to auroral diagnostics.  In case one of nebular pair is not observed, it assumes the intensity of the other is given by the theoretical ratio
        implicit none
        integer, parameter :: dp = kind(1.d0)

        character(len=11) :: name1, name2, name3
        character(len=20) :: ion
        integer :: ion_no1, ion_no2, ion_no3
        real(kind=dp) :: flux1, flux2, flux3
        real(kind=dp) :: factor1, factor2, ratio, ratio2
        integer :: level1, level2, level3 !level2 is the upper level, level1 the lower of the two lower levels

!debugging
#ifdef CO
        print *,"subroutine: get_Tdiag. lines=",name1, name2, name3
#endif

        ion_no1 = get_ion(name1, ILs)
        ion_no2 = get_ion(name2, ILs)
        ion_no3 = get_ion(name3, ILs)

        flux1 = get_cel_flux(name1, linelist, ILs)
        flux2 = get_cel_flux(name2, linelist, ILs)
        flux3 = get_cel_flux(name3, linelist, ILs)

        !get the levels of the transitions. upper level is common to both transitions
        read (ILs(ion_no1)%transition(1:index(ILs(ion_no1)%transition,',')-1),"(i2)") level1
        read (ILs(ion_no1)%transition(index(ILs(ion_no1)%transition,',')+1:index(ILs(ion_no1)%transition,'/')-1),"(i2)") level2
        read (ILs(ion_no2)%transition(1:index(ILs(ion_no2)%transition,',')-1),"(i2)") level3

        !now, in possibly the ugliest lines of code in all of NEAT, get the theoretical ratio, which is ratio of A values multiplied by ratio of wavelengths

        factor1=1.+(atomicdata(get_atomicdata(ion,atomicdata))%A_coeffs(level2,level3)/atomicdata(get_atomicdata(ion,atomicdata))%A_coeffs(level2,level1) * (atomicdata(get_atomicdata(ion,atomicdata))%waveno(level2)-atomicdata(get_atomicdata(ion,atomicdata))%waveno(level3))/(atomicdata(get_atomicdata(ion,atomicdata))%waveno(level2)-atomicdata(get_atomicdata(ion,atomicdata))%waveno(level1)))
        factor2=1.+(1./(factor1-1.))

        !now calculate the diagnostic ratio

        if(((flux1 .gt. 0) .and. (flux2 .gt. 0)) .and. (flux3 .gt. 0 ))then

          ratio = (flux1 + flux2) / flux3

          if(name1 .eq. "siii9069   ")then
            ratio2 = (flux1 / flux2)
            !print*, ratio2
            if(ratio2 .gt. (factor2-1)*0.95 .and. ratio2 .lt. (factor2-1)*1.05)then
              !PRINT*, ratio2, " ", factor2-1, " 1 ", name1
              ratio = (flux1 + flux2) / flux3
            else if(ratio2 .lt. (factor2-1)*0.95)then
              !PRINT*, ratio2, " ", factor2-1, " 2 ", name1
              ratio = (factor2 * flux2) / flux3
              !PRINT*, (factor2 * flux2), " ", (flux1 + flux2)
            else if(ratio2 .gt. (factor2-1)*1.05)then
              !PRINT*, ratio2, " ", factor2-1, " 3 ", name1
              ratio = (factor1 * flux1) / flux3
            else
              ratio = 0.d0
            endif
          endif

        elseif(((flux1 .gt. 0) .and. (flux2 .eq. 0)) .and. (flux3 .gt. 0 ))then
          ratio = (flux1 * factor1) / flux3
        elseif(((flux1 .eq. 0) .and. (flux2 .gt. 0)) .and. (flux3 .gt. 0 ))then
          ratio = (flux2 * factor2) / flux3
        else
          ratio = 0.d0
        endif

end subroutine get_Tdiag

end subroutine abundances
