subroutine abundances(fname1, run, switch_ext)
use mod_abundmaths
use mod_abundtypes
use mod_diagnostics
use mod_getabunds
use mod_abundIO
use mod_helium
use mod_recombination_lines
use mod_extinction

implicit none

        INTEGER :: count, Iint, fIL, i, j, ion_no1, ion_no2, ion_no3, ion_no4, iii, ion_no5, ion_no6
        INTEGER :: opt, runonce, run
        CHARACTER*80, INTENT(IN) :: fname1 !input filenames
        CHARACTER*80 :: fnameIL, blank
        CHARACTER*8 :: lion
        CHARACTER :: switch_ext !switch for extinction laws
        DOUBLE PRECISION :: normalise, oiiNratio, oiiDens, oiiiTratio, oiiiTemp, oiiiIRNratio, oiiiIRTratio, niiTratio, niiTemp, arivNratio, arivDens, cliiiNratio, cliiiDens, siiNratio, siiDens, siiTratio, siiTemp, oiiTratio, oiiTemp, neiiiTratio, neiiiIRTratio, neiiiTemp, neiiiIRTemp, abund, meandensity, meantemp, oitemp, citemp
        DOUBLE PRECISION :: ciiiNratio,neivNratio,nevTratio,siiiTratio,ariiiTratio,arvTratio,lowtemp,lowdens,medtemp,ciiidens,meddens,siiitemp,ariiitemp,hightemp,neivdens,highdens,arvtemp,nevtemp,oiTratio,ciTratio
        DOUBLE PRECISION :: oiiRLabund, niiRLabund, ciiRLabund, cii4267rlabund, neiiRLabund, ciiiRLabund, niiiRLabund, RLabundtemp, weight
        DOUBLE PRECISION :: ciiiCELabund, niiCELabund, niiiIRCELabund, niiiUVCELabund, oiiCELabund, oiiiCELabund, oiiiIRCELabund, neiiIRCELabund, neiiiIRCELabund, neiiiCELabund, neivCELabund, siiCELabund, siiiCELabund, siiiIRCELabund, sivIRCELabund, cliiiCELabund, ariiiCELabund, arivCELabund, ariiiIRCELabund, nivCELabund, niCELabund, niiiCELabund, ciiCELabund, civCELabund, nvCELabund, nevCELabund, arvCELabund, CELabundtemp, ciCELabund, oiCELabund 
        DOUBLE PRECISION :: CELicfO, CELicfC, CELicfN, CELicfNe, CELicfAr, CELicfS
        DOUBLE PRECISION :: RLicfO, RLicfC, RLicfN, RLicfNe
        DOUBLE PRECISION :: CabundRL, CabundCEL, NabundRL, NabundCEL, OabundRL, OabundCEL, NeabundRL, NeabundCEL, SabundCEL, ArabundCEL, NOabundCEL, NCabundCEL 
        DOUBLE PRECISION :: adfC, adfN, adfO, adfNe, w1, w2, w3, w4
        DOUBLE PRECISION :: adfC2plus, adfN2plus, adfO2plus, adfNe2plus
        DOUBLE PRECISION :: c1, c2, c3, meanextinction, fl, ratob, tempi, temp, temp2, A4471, A4686, A6678, A5876
        REAL :: heiabund,heiiabund,Hetotabund
        REAL*8 :: HW

        DOUBLE PRECISION, DIMENSION(2) :: conditions 
        REAL*8, DIMENSION(3,335) :: linelist
        REAL*8, DIMENSION(2,335) :: linelist_dered
        REAL*8 :: result       

        TYPE(line), DIMENSION(62) :: ILs
        TYPE(line), DIMENSION(4) :: H_BS
        TYPE(line), DIMENSION(4) :: He_lines

! recombination line variables

        TYPE RLabund
           CHARACTER*7 :: Multiplet
           REAL*8 :: Abundance
        END TYPE

        TYPE (RLabund), DIMENSION(12) :: oiimultiplets
        TYPE (RLabund), DIMENSION(7) :: niimultiplets

! strong line variables
        DOUBLE PRECISION :: X23,O_R23upper, O_R23lower, N2,O_N2, O3N2, O_O3N2, Ar3O3, O_Ar3O3, S3O3, O_S3O3

        linelist = 0
        ILs%intensity = 0
        H_BS%intensity = 0
        He_lines%intensity = 0
        !runonce = 1 !allows printing of supplementary files
        runonce = run !suppresses supplementary files and enables monte-carlo error estimation
        nivCELabund = 0.0
        nvCELabund = 0.0
        civCELabund = 0.0

        !file reading stuff

        !reading in Rogers "important" lines list 

        CALL read_ilines(ILs, Iint) 
        CALL fileread(linelist, fname1) ! see above 
        CALL element_assign(ILs, linelist, Iint)

        !dereddening

        ILs%abundance = 0
        ILs%int_dered = 0
        ILs%int_dered_err = 0
        H_BS%abundance = 0
        H_BS%int_dered = 0
        H_BS%int_dered_err = 0
        He_lines%abundance = 0
        He_lines%int_dered = 0
        He_lines%int_dered_err = 0

        !first lets find some hydrogen lines 
        CALL get_H(H_BS, linelist) 
        call get_He(He_lines, linelist)

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
        print *, "=========="

        CALL calc_extinction_coeffs(H_BS, c1, c2, c3, meanextinction)

        !need to write output/ input stuff so user can insert own c(Hb)
        !assume we go on with calculated extinctions

        print "(1X,A17,F4.2,A4,F4.2)","Ha/Hb => c(Hb) = ",c1
        print "(1X,A17,F4.2,A4,F4.2)","Hg/Hb => c(Hb) = ",c2
        print "(1X,A17,F4.2,A4,F4.2)","Hd/Hb => c(Hb) = ",c3

        PRINT "(1X,A13,F4.2,A4,F4.2)", "Mean c(Hb) = ",meanextinction
        if (runonce == 0 .and. meanextinction > 0) write(UNIT=888, FMT=*)meanextinction

        if (meanextinction .lt. 0.0) then
           print *,"Derived extinction <0 ; assuming 0"
           meanextinction = 0.0
        endif

        !actual dereddening
        if (switch_ext == "S") then
                CALL deredden(ILs, Iint, meanextinction)
                CALL deredden(H_BS, 4, meanextinction)
                call deredden(He_lines, 4, meanextinction) 
                CALL deredden_O(linelist, linelist_dered, meanextinction)
        elseif (switch_ext == "H") then
                CALL deredden_LMC(ILs, Iint, meanextinction)
                CALL deredden_LMC(H_BS, 4, meanextinction)
                call deredden_LMC(He_lines, 4, meanextinction) 
                CALL deredden_O_LMC(linelist, linelist_dered, meanextinction)
        elseif (switch_ext == "C") then
                CALL deredden_CCM(ILs, Iint, meanextinction)
                CALL deredden_CCM(H_BS, 4, meanextinction)
                call deredden_CCM(He_lines, 4, meanextinction) 
                CALL deredden_O_CCM(linelist, linelist_dered, meanextinction)
        elseif (switch_ext == "P") then
                CALL deredden_SMC(ILs, Iint, meanextinction)
                CALL deredden_SMC(H_BS, 4, meanextinction)
                call deredden_SMC(He_lines, 4, meanextinction) 
                CALL deredden_O_SMC(linelist, linelist_dered, meanextinction)
        endif

        500 FORMAT (5(f10.4))
        if(runonce == 1) OPEN(801, FILE=trim(fname1)//"_dered", STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
        !WRITE(801,'(A40)') "sjfnb;ashfb"
        if(runonce == 1)then
                do iii=1, Iint

                        if(ILs(iii)%int_dered .ne. 0) write(801,500) ILs(iii)%wavelength, ILs(iii)%intensity, ILs(iii)%int_err, ILs(iii)%int_dered, ILs(iii)%int_dered_err*ILs(iii)%int_dered
                end do
                do iii=1,4
                        if( He_lines(iii)%int_dered .ne. 0 ) write(801,500) He_lines(iii)%wavelength, He_lines(iii)%intensity, He_lines(iii)%int_err, He_lines(iii)%int_dered, He_lines(iii)%int_dered_err*He_lines(iii)%int_dered

                        if( H_BS(iii)%int_dered .ne. 0 ) write(801,500) H_BS(iii)%wavelength, H_BS(iii)%intensity, H_BS(iii)%int_err, H_BS(iii)%int_dered, H_BS(iii)%int_dered_err*H_BS(iii)%int_dered
                end do
        endif        
        if(runonce == 1) CLOSE(801)
        if(runonce == 1) call system("sort "//trim(fname1)//"_dered > "//trim(fname1)//"_dered_sort")
        if(runonce == 1) call system("rm "//trim(fname1)//"_dered")

!diagnostics
        call get_diag("ciii1909   ","ciii1907   ", ILs, ciiiNratio)        ! ciii ratio   
        call get_diag("oii3729    ","oii3726    ", ILs, oiiNratio )        ! oii ratio
        call get_diag("neiv2425   ","neiv2423   ", ILs, neivNratio )        ! neiv ratio
        call get_diag("sii6731    ","sii6716    ", ILs, siiNratio )        ! s ii ratio
        call get_diag("cliii5537  ","cliii5517  ", ILs, cliiiNratio )        ! Cl iii ratio
        call get_diag("ariv4740   ","ariv4711   ", ILs, arivNratio )        ! Ar iv ratio

! temperature ratios:
        !Try to calculate from atomic data at start
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
                PRINT*, "OII missing, Wi n�t trei a h�liday in Sweden this yer? "
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

         count = 0
         if (oiiNratio .gt. 0 .and. oiiNratio .lt. 1e10) then 
           call get_diagnostic("oii       ","1,2/                ","1,3/                ",oiiNratio,"D",lowtemp, oiiDens) 
           count = count + 1
         endif
         if (siiNratio .gt. 0 .and. siiNratio .lt. 1e10) then
           call get_diagnostic("sii       ","1,2/                ","1,3/                ",siiNratio,"D",lowtemp, siiDens)
           count = count + 1
         endif

         if (count .eq. 0) then
           lowdens = 1000.0
         else
           lowdens = (oiiDens + siiDens) / count
         endif

         count = 0

         if (oiiTratio .gt. 0 .and. oiiTratio .lt. 1e10) then
           call get_diagnostic("oii       ","2,4,2,5,3,4,3,5/    ","1,2,1,3/            ",oiiTratio,"T",lowdens,oiiTemp)
           count = count + 1
                 if(oiitemp == 20000)then
                        count=count-1
                        oiitemp=-1
                 endif 

         else
           oiiTemp = 0.0
         endif
         if (siiTratio .gt. 0 .and. siiTratio .lt. 1e10) then
           call get_diagnostic("sii       ","1,2,1,3/            ","1,4,1,5/            ",siiTratio,"T",lowdens,siiTemp)
           count = count + 1

                if(siitemp == 20000)then
                        count=count-1
                        siitemp=-1
                 endif 

         else
           siiTemp = 0.0
         endif
         if (niiTratio .gt. 0 .and. niiTratio .lt. 1e10) then
           call get_diagnostic("nii       ","2,4,3,4/            ","4,5/                ",niiTratio,"T",lowdens,niitemp)
           count = count + 5

                 if(niitemp == 20000)then
                        count=count-1
                        niitemp=-1
                 endif 

         else
           niitemp = 0.0
         endif

         if (ciTratio .gt. 0 .and. ciTratio .lt. 1e10) then
           call get_diagnostic("ci        ","2,4,3,4/            ","4,5/                ",ciTratio,"T",lowdens,citemp)
           count = count + 1

                if(citemp == 20000)then
                        count=count-1
                        citemp=-1
                 endif 

          else
           citemp = 0.0
         endif

          if (oiTratio .gt. 0 .and. oiTratio .lt. 1e10) then
           call get_diagnostic("oi        ","1,4,2,4/            ","4,5/                ",oiTratio,"T",lowdens,oitemp)
           count = count + 1

         if(oitemp == 20000)then
                 count=count-1
                 oitemp=-1
         endif
         else
           oitemp = 0.0
         endif         

         if (count .gt. 0) then 
           lowtemp = ((5*niitemp) + siitemp + oiitemp + oitemp + citemp) / count
         else
           lowtemp = 10000.0
         endif

      enddo



! medium ionisation
        cliiiDens = 0
        ciiiDens = 0
        arivDens = 0


      medtemp = 10000.0

      do i = 1,2

         count = 0
         if (cliiiNratio .gt. 0 .and. cliiiNratio .lt. 1e10) then
           call get_diagnostic("cliii     ","1,2/                ","1,3/                ",cliiiNratio,"D",medtemp, cliiiDens)

           count = count + 1
         endif
         if (ciiiNratio .gt. 0 .and. ciiiNratio .lt. 1e10) then
           call get_diagnostic("ciii      ","1,2/                ","1,3/                ",ciiiNratio,"D",medtemp, ciiiDens)
           count = count + 1
         endif
         if (arivNratio .gt. 0 .and. arivNratio .lt. 1e10) then
           call get_diagnostic("ariv      ","1,2/                ","1,3/                ",arivNratio,"D",medtemp, arivDens)
           count = count + 1
         endif

         if (count .eq. 0) then
           meddens = 1000.0
         else
           meddens = (ciiiDens + cliiiDens + arivDens) / count
         endif

         count = 0

         if (oiiiTratio .gt. 0 .and. oiiiTratio .lt. 1e10) then
           call get_diagnostic("oiii      ","2,4,3,4/            ","4,5/                ",oiiiTratio,"T",meddens,oiiiTemp)
           count = count + 4
                 if(oiiitemp == 20000)then
                         count=count-4
                         oiiitemp=-1
                 endif

         else
           oiiiTemp = 0.0
         endif
         if (siiiTratio .gt. 0 .and. siiiTratio .lt. 1e10) then
           call get_diagnostic("siii      ","2,4,3,4/            ","4,5/                ",siiiTratio,"T",meddens,siiiTemp)
           count = count + 1

                if(siiitemp == 20000)then
                         count=count-1
                         siiitemp=-1
                 endif

         else
           siiiTemp = 0.0
         endif
         if (ariiiTratio .gt. 0 .and. ariiiTratio .lt. 1e10) then
           call get_diagnostic("ariii     ","1,4,2,4/            ","4,5/                ",ariiiTratio,"T",meddens,ariiitemp)
           count = count + 2
                 if(ariiitemp == 20000)then
                         count=count-2
                         ariiitemp=-1
                 endif
         else
           ariiitemp = 0.0
         endif
         if (neiiiTratio .gt. 0 .and. neiiiTratio .lt. 1e10) then
           call get_diagnostic("neiii     ","1,4,2,4/            ","4,5/                ",neiiiTratio,"T",meddens,neiiitemp)
           count = count + 2
                if(neiiitemp == 20000)then
                         count=count-2
                         neiiitemp=-1
                 endif

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

         if (neivNratio .gt. 0 .and. neivNratio .lt. 1e10) then
           call get_diagnostic("neiv      ","1,2/                ","1,3/                ",neivNratio,"D",hightemp, neivDens)
           highdens = neivdens
         else
           neivDens = 0.0
           highdens = 1000.0
         endif

         count = 0

         if (arvTratio .gt. 0 .and. arvTratio .lt. 1e10) then
           call get_diagnostic("arv       ","2,4,3,4/            ","4,5/                ",arvTratio,"T",highdens,arvTemp)
           count = count + 1
         else
           arvTemp = 0.0
         endif
         if (nevTratio .gt. 0 .and. nevTratio .lt. 1e10) then
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


!edited diagnostic section to reflect above changes which stopped high/low limit densities/temperatures being included in averages. High limit cases set to 0.1 so that we know that it was in the high limit and not the low limit. DJS


      print *,""
      print *,"Diagnostics"
      print *,"==========="
      print *,""

      print *,"Diagnostic       Zone      Value    Diagnostic ratio"
      print *,""
if(oiidens >0)     print "(A28,F8.0,A1,F8.3)","[O II] density   Low       ",oiidens, " ", oiiNratio
        if(runonce == 0 .and. oiidens > 0) WRITE(unit=865,FMT=*)oiidens
if(siidens >0)     print "(A28,F8.0,A1,F8.3)","[S II] density   Low       ",siidens, " ", REAL(1/siiNratio)
        if(runonce == 0 .and. siidens > 0) WRITE(unit=866,FMT=*)siidens
if(lowdens >0)     print "(A28,F8.0)"," density adopted Low       ",lowdens
        if(runonce == 0 .and. lowdens > 0) WRITE(unit=867,FMT=*)lowdens
if(lowdens >0)     print *,""



if(niitemp > 0.2)then
                   print "(A28,F8.0,A1,F8.3)","[N II] temp      Low       ",niitemp, " ", niitratio 
        if(runonce == 0)WRITE(UNIT=868,FMT=*)niitemp
else if(INT(niitemp) == -1)then
                   print "(A28,F8.0,A1,F8.3)","[N II] temp      Low       ",20000.0, " ", niitratio 
        if(runonce == 0)WRITE(UNIT=868,FMT=*)"20000"
else

endif

if(oiitemp >0.2)then
                   print "(A28,F8.0,A1,F8.3)","[O II] temp      Low       ",oiitemp, " ", REAL(1/oiitratio) 
        if(runonce == 0)WRITE(UNIT=869,FMT=*)oiitemp
else if(INT(oiitemp) == -1)then
                   print "(A28,F8.0,A1,F8.3)","[O II] temp      Low       ",20000.0, " ", REAL(1/oiitratio) 
        if(runonce == 0)WRITE(UNIT=869,FMT=*)"20000"
else
endif 

if(siitemp >0.2 )then
                   print "(A28,F8.0,A1,F8.3)","[S II] temp      Low       ",siitemp, " ", siitratio 
        if(runonce == 0)WRITE(UNIT=870,FMT=*)siitemp
else if(INT(siitemp) == -1)then
                   print "(A28,F8.0,A1,F8.3)","[S II] temp      Low       ",20000.0, " ", siitratio 
        if(runonce == 0)WRITE(UNIT=870,FMT=*)"20000"
else
endif 

if(oitemp >0.2 )then
                   print "(A28,F8.0,A1,F8.3)","[O I]  temp      Low       ",oitemp, " ", oitratio  
        if(runonce == 0)WRITE(UNIT=871,FMT=*)oitemp
else if(INT(oitemp) == -1)then
                   print "(A28,F8.0,A1,F8.3)","[O I]  temp      Low       ",20000.0, " ", oitratio 
        if(runonce == 0)WRITE(UNIT=871,FMT=*)"20000"
else
endif     

if(citemp >0.2 )then
                   print "(A28,F8.0,A1,F8.3)","[C I]  temp      Low       ",citemp, " ", citratio  
        if(runonce == 0)WRITE(UNIT=872,FMT=*)citemp
else if(INT(citemp) == -1)then
                   print "(A28,F8.0,A1,F8.3)","[C I]  temp      Low       ",20000.0, " ", citratio 
        if(runonce == 0)WRITE(UNIT=872,FMT=*)"20000"
else
endif                    


if(lowtemp >0)     print "(A28,F8.0)"," temp adopted    Low       ",lowtemp
        if(runonce == 0 .and. lowtemp > 0)WRITE(UNIT=873,FMT=*)lowtemp
if(lowtemp >0)     print *,""

if(cliiidens > 0 )    print "(A28,F8.0,A1,F8.3)","[Cl III] density Medium    ",cliiidens," ", cliiinratio
        if(runonce == 0 .and. cliiidens > 0 )WRITE(UNIT=874,FMT=*)cliiidens
if(arivdens  > 0 )    print "(A28,F8.0,A1,F8.3)","[Ar IV] density  Medium    ",arivdens," ", arivnratio
        if(runonce == 0 .and. arivdens > 0 )WRITE(UNIT=875,FMT=*)arivdens
if(ciiidens  > 0 )    print "(A28,F8.0,A1,F8.3)","C III] density   Medium    ",ciiidens," ", ciiinratio
        if(runonce == 0 .and. ciiidens > 0 )WRITE(UNIT=876,FMT=*)ciiidens
if(meddens   > 0 )    print "(A28,F8.0)"," density adopted Medium    ",meddens
        if(runonce == 0 .and. meddens > 0 )WRITE(UNIT=877,FMT=*)meddens
if(meddens   > 0 )    print *,""

if(oiiitemp >0.2)then
                   print "(A28,F8.0,A1,F8.3)","[O III] temp     Medium    ",oiiitemp, " ",oiiitratio
        if(runonce == 0)WRITE(UNIT=878,FMT=*)oiiitemp
else if(INT(oiiitemp) == -1)then 
                   print "(A28,F8.0,A1,F8.3)","[O III] temp     Medium    ",20000.0, " ",oiiitratio
        if(runonce == 0)WRITE(UNIT=878,FMT=*)"20000"
else
endif

if(neiiitemp>0.2)then
                   print "(A28,F8.0,A1,F8.3)","[Ne III] temp    Medium    ",neiiitemp, " ",neiiitratio
        if(runonce == 0)WRITE(UNIT=879,FMT=*)neiiitemp
else if(INT(neiiitemp) == -1)then
                   print "(A28,F8.0,A1,F8.3)","[Ne III] temp    Medium    ",20000.0, " ",neiiitratio
        if(runonce == 0)WRITE(UNIT=879,FMT=*)"20000"
else
endif                   
if(ariiitemp>0.2)then
                   print "(A28,F8.0,A1,F8.3)","[Ar III] temp    Medium    ",ariiitemp, " ",ariiitratio
        if(runonce == 0)WRITE(UNIT=880,FMT=*)ariiitemp
else if(INT(ariiitemp) == -1)then
                   print "(A28,F8.0,A1,F8.3)","[Ar III] temp    Medium    ",20000.0, " ",ariiitratio
        if(runonce == 0)WRITE(UNIT=880,FMT=*)"20000"
else
endif
if(siiitemp > 0.2)then
                   print "(A28,F8.0,A1,F8.3)","[S III] temp     Medium    ",siiitemp, " ",siiitratio
        if(runonce == 0)WRITE(UNIT=881,FMT=*)siiitemp
else if(int(siiitemp) == -1)then
                   print "(A28,F8.0,A1,F8.3)","[S III] temp     Medium    ",20000.0, " ",siiitratio
        if(runonce == 0)WRITE(UNIT=881,FMT=*)"20000"
else
endif                   

if(medtemp  >0)    print "(A28,F8.0)"," temp adopted    Medium    ",medtemp
        if(runonce == 0 .and. medtemp > 0)WRITE(UNIT=882,FMT=*)medtemp
if(medtemp  >0)    print *,""

if(neivdens >0)    print "(A28,F8.0,A1,F8.3)","[Ne IV] density  High      ",neivdens, " ",neivnratio
        if(runonce == 0 .and. neivdens > 0)WRITE(UNIT=883,FMT=*)neivdens
if(highdens >0)    print "(A28,F8.0)"," density adopted High      ",highdens
        if(runonce == 0 .and. highdens > 0)WRITE(UNIT=884,FMT=*)highdens
if(highdens >0)    print *,""
if(arvtemp  >0)    print "(A28,F8.0,A1,F8.3)","[Ar V] temp      High      ",arvtemp, " ",arvtratio
        if(runonce == 0 .and. arvtemp > 0)WRITE(UNIT=885,FMT=*)arvtemp
if(nevtemp  >0)    print "(A28,F8.0,A1,F8.3)","[Ne V] temp      High      ",nevtemp, " ",nevtratio
        if(runonce == 0 .and. nevtemp > 0)WRITE(UNIT=886,FMT=*)nevtemp
if(hightemp >0)    print "(A28,F8.0)"," temp adopted    High      ",hightemp
        if(runonce == 0 .and. hightemp > 0)WRITE(UNIT=882,FMT=*)hightemp

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
        print *,"Ionic abundances"
        print *,"=========="

        print *,"Helium"
        print *,"------"

        call get_helium(REAL(medtemp),REAL(meddens),REAL(He_lines(1)%int_dered),REAL(He_lines(2)%int_dered),REAL(He_lines(3)%int_dered),REAL(He_lines(4)%int_dered),heiabund,heiiabund,Hetotabund, A4471, A4686, A6678, A5876)

if(A4471 > 0)   print "(1x,A17,F6.4)", " He+ (4471)/H+ = ", A4471
if(A5876 > 0)        print "(1x,A17,F6.4)", " He+ (5876)/H+ = ", A5876        
if(A6678 > 0)        print "(1x,A17,F6.4)", " He+ (6678)/H+ = ", A6678
if(A4686 > 0)        print "(1x,A17,F6.4)", "He++ (4686)/H+ = ", A4686 

        if( (A4471 > 0 .or. A5876 > 0 ) .or. A6678 > 0)then

                if(He_lines(2)%intensity > 0) w1 = 1/((He_lines(2)%int_err / He_lines(2)%intensity)**2)
                if(He_lines(3)%intensity > 0) w2 = 1/((He_lines(3)%int_err / He_lines(3)%intensity)**2)
                if(He_lines(4)%intensity > 0) w3 = 1/((He_lines(4)%int_err / He_lines(4)%intensity)**2)        

                !PRINT*, w1, " ", w2, " ", w3

                heiabund = (w1*A4471 + w2*A5876 + w3*A6678)/(w1+w2+w3)

        else

                heiabund = 0.0
        endif                                

        print "(1X,A17,F6.4)", "        He+/H+ = ",heiabund
        print "(1X,A17,F6.4)", "       He++/H+ = ",heiiabund
        print "(1X,A17,F6.4)", "          He/H = ",heiabund + heiiabund

        w1=0
        w2=0
        w3=0
        w4=0


! get abundances for all CELs

        print *,""
        print *,"CELs"
        print *,"----"
        print *,"Ion         I(lambda)   Abundance"


        !This routine is too simple. I have been changing the temperatures /densities which are input to each zone to disable the zone schtick.
        !Make a better routine that allows using or not using the zone thing.. Its not always appropriate.


        !if(oiiitemp > 0 )then
        !        medtemp = oiiitemp
        !        siiitemp = oiiitemp
        !endif

        !if(int(siiitemp) == -1) siiitemp = oiiitemp


        do i = 1,Iint-1 !XXXX why does Int end up 1 too high?
!                 print *,ILs(i)%ion,ILs(i)%transition,ILs(i)%int_dered
           if (ILs(i)%zone .eq. "low ") then
                !PRINT*, siiitemp, lowdens
                 call get_abundance(ILs(i)%ion, ILs(i)%transition, lowtemp, lowdens,ILs(i)%int_dered, ILs(i)%abundance)
                 ! elseif ( ( i== 47 .or. (i == 46 .or. i == 28 ) ) .and. siiitemp > 1.0 ) then
          !       call get_abundance(ILs(i)%ion, ILs(i)%transition, siiitemp, meddens,ILs(i)%int_dered, ILs(i)%abundance)
                !this makes the code use siii temperatures for siii
                ! print*, "using this bit" 
           elseif (ILs(i)%zone .eq. "med ") then
                 call get_abundance(ILs(i)%ion, ILs(i)%transition, medtemp, meddens,ILs(i)%int_dered, ILs(i)%abundance)
           elseif (ILs(i)%zone .eq. "high") then
                 call get_abundance(ILs(i)%ion, ILs(i)%transition, hightemp, highdens,ILs(i)%int_dered, ILs(i)%abundance)
           endif
                 if ((ILs(i)%abundance .gt. 0) .and. (ILs(i)%abundance < 1 ) ) then
                       PRINT "(1X, A11, 1X, F7.3, 5X, ES10.4)",ILs(i)%name,ILs(i)%int_dered,ILs(i)%abundance
                 elseif( (ILs(i)%abundance > 1 ) .or. (ILs(i)%abundance < 1E-10 ) )then
                         ILs(i)%abundance = 0        
                 endif
        enddo

! calculate averages

        ciiiCELabund = ILs( get_ion("ciii1909   ", ILs, Iint)  )%abundance


        celabundtemp = 0.
        weight = 0.
        do i= get_ion("nii5754    ", ILs, Iint), get_ion("nii6584    ", ILs, Iint)
          if (ILs(i)%abundance .gt. 0) niiCELabund = niiCELabund + ILs(i)%abundance/ ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          if (ILs(i)%abundance .gt. 0) weight = weight + 1/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
        enddo
        if (weight .gt. 0) then
          niiCELabund = niiCELabund / weight
        else
          niiCELabund = 0.0
        endif

        niiiIRCELabund = ILs( get_ion("niii57um   ", ILs, Iint)  )%abundance
        niiiUVCELabund = ILs( get_ion("niii1751   ", ILs, Iint)   )%abundance

        if (niiiIRCELabund .gt. 0 .and. niiiUVCELabund .gt. 0) then
          niiiCELabund = (niiiIRCELabund + niiiUVCELabund)/2
        elseif (niiiIRCELabund .gt. 0) then
          niiiCELabund = niiiIRCELabund
        elseif (niiiUVCELabund .gt. 0) then
          niiiCELabund = niiiUVCELabund
        else
          niiiCELabund = 0
        endif


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
        weight = 0.
        do i=get_ion("oiii4959   ", ILs, Iint), get_ion("oiii5007   ", ILs, Iint)
          if (ILs(i)%abundance .gt. 0) oiiiCELabund = oiiiCELabund + ILs(i)%abundance /((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          if (ILs(i)%abundance .gt. 0) weight = weight + 1 / ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
        enddo
        if (weight .gt. 0) then
          oiiiCELabund = oiiiCELabund / weight
        else
          oiiiCELabund = 0.0
        endif

        celabundtemp = 0.
        weight = 0.
        do i=get_ion("oiii52um   ", ILs, Iint), get_ion("oiii88um   ", ILs, Iint)
          if (ILs(i)%abundance .gt. 0) oiiiIRCELabund = oiiiIRCELabund + ILs(i)%abundance / ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          if (ILs(i)%abundance .gt. 0) weight = weight + 1/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
        enddo
        if (weight .gt. 0) then
          oiiiIRCELabund = oiiiIRCELabund / weight
        else
          oiiiIRCELabund = 0.0
        endif

        neiiIRCELabund = ILs(  get_ion("neii12p8um ", ILs, Iint)  )%abundance
        neiiiIRCELabund = ILs(  get_ion("neiii15p5um ", ILs, Iint)  )%abundance

        celabundtemp = 0.
        weight = 0.
        do i=get_ion("neiii3868  ", ILs, Iint), get_ion("neiii3967  ", ILs, Iint)
          if (ILs(i)%abundance .gt. 0) neiiiCELabund = neiiiCELabund + ILs(i)%abundance/ ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          if (ILs(i)%abundance .gt. 0) weight = weight + 1/ ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
        enddo
        if (weight .gt. 0) then
          neiiiCELabund = neiiiCELabund / weight
        else
          neiiiCELabund = 0.0
        endif


        celabundtemp = 0.
        weight = 0.
        do i=get_ion("neiv4724   ", ILs, Iint), get_ion("neiv4715   ", ILs, Iint)
          if (ILs(i)%abundance .gt. 0) neivCELabund = neivCELabund + ILs(i)%abundance/ ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          if (ILs(i)%abundance .gt. 0) weight = weight + 1/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
        enddo
        if (weight .gt. 0) then
          neivCELabund = neivCELabund / weight
        else
          neivCELabund = 0.0
        endif

        celabundtemp = 0.
        weight = 0.
        do i=get_ion("sii4068    ", ILs, Iint), get_ion("sii6731    ", ILs, Iint)
          if (ILs(i)%abundance .gt. 0) siiCELabund = siiCELabund + ILs(i)%abundance/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          if (ILs(i)%abundance .gt. 0) weight = weight + 1/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
        enddo
        if (weight .gt. 0) then
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


        siiiIRCELabund = ILs( get_ion("siii18p7um ", ILs, Iint)   )%abundance
        sivIRCELabund = ILs(  get_ion("siv10p5um  ", ILs, Iint) )%abundance

        celabundtemp = 0.
        weight = 0.
        do i=get_ion("cliii5517  ", ILs, Iint), get_ion("cliii5537  ", ILs, Iint)
          if (ILs(i)%abundance .gt. 0) cliiiCELabund = cliiiCELabund + ILs(i)%abundance/ ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          if (ILs(i)%abundance .gt. 0) then
            weight = weight + 1/ ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          endif
        enddo
        if (weight .gt. 0) then
          cliiiCELabund = cliiiCELabund / weight
        else
          cliiiCELabund = 0.0
        endif

        celabundtemp = 0.
        weight = 0.
        do i=get_ion("ariii7135  ", ILs, Iint), get_ion("ariii7751  ", ILs, Iint)
          if (ILs(i)%abundance .gt. 0) ariiiCELabund = ariiiCELabund + ILs(i)%abundance/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          if (ILs(i)%abundance .gt. 0) then
            weight = weight + 1/ ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          endif
        enddo
        if (weight .gt. 0) then
          ariiiCELabund = ariiiCELabund / weight
        else
          ariiiCELabund = 0.0
        endif

        celabundtemp = 0.
        weight = 0.
         do i=get_ion("ariv4711   ", ILs, Iint), get_ion("ariv4740   ", ILs, Iint)
          arivCELabund = arivCELabund + ILs(i)%abundance/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          if (ILs(i)%abundance .gt. 0) then
            weight = weight + 1/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          endif
        enddo
        if (weight .gt. 0) then
          arivCELabund = arivCELabund / weight
        else
          arivCELabund = 0.0
        endif

        ariiiIRCELabund = ILs(get_ion("ariii9um   ", ILs, Iint))%abundance

        celabundtemp = 0.
        weight = 0.
        do i=get_ion("ciii1907   ", ILs, Iint), get_ion("ciii1909   ", ILs, Iint)
          if (ILs(i)%abundance .gt. 0) then
            ciiiCELabund = ciiiCELabund + ILs(i)%abundance/ ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
            weight = weight + 1/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          endif
        enddo
        if (weight .gt. 0) then
          ciiiCELabund = ciiiCELabund / weight
        else
          ciiiCELabund = 0.0
        endif

        celabundtemp = 0.
        weight = 0.
        do i=get_ion("neiv2423   ", ILs, Iint), get_ion("neiv2425   ", ILs, Iint)
          if (ILs(i)%abundance .gt. 0) neivCELabund = neivCELabund + ILs(i)%abundance/ ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          if (ILs(i)%abundance .gt. 0) then
            weight = weight + 1/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          endif
        enddo
        if (weight .gt. 0) then
          neivCELabund = neivCELabund / weight
        else
          neivCELabund = 0.0
        endif

        celabundtemp = 0.
        weight = 0.
        do i=get_ion("nev3345    ", ILs, Iint), get_ion("nev3426    ", ILs, Iint)
          if (ILs(i)%abundance .gt. 0) nevCELabund = nevCELabund + ILs(i)%abundance/ ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          if (ILs(i)%abundance .gt. 0) then
            weight = weight + 1/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          endif
        enddo
        if (weight .gt. 0) then
          nevCELabund = nevCELabund / weight
        else
          nevCELabund = 0.0
        endif

        celabundtemp = 0.
        weight = 0.
        do i=get_ion("arv6435    ", ILs, Iint), get_ion("arv7005    ", ILs, Iint)
          if (ILs(i)%abundance .gt. 0) arvCELabund = arvCELabund + ILs(i)%abundance/ ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          if (ILs(i)%abundance .gt. 0) then
            weight = weight + 1/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          endif
        enddo
        if (weight .gt. 0) then
          arvCELabund = arvCELabund / weight
        else
          arvCELabund = 0.0
        endif

         celabundtemp = 0.
        weight = 0.
        do i=get_ion("ci9850     ", ILs, Iint), get_ion("ci8727     ", ILs, Iint)
          if (ILs(i)%abundance .gt. 0) ciCELabund = ciCELabund + ILs(i)%abundance/ ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          if (ILs(i)%abundance .gt. 0) then
            weight = weight + 1/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          endif
        enddo
        if (weight .gt. 0) then
          ciCELabund = ciCELabund / weight
        else
          ciCELabund = 0.0
        endif

        NCabundCEL = ciCELabund

         celabundtemp = 0.
        weight = 0.
        do i=get_ion("oi6300     ", ILs, Iint), get_ion("oi5577     ", ILs, Iint)
          if (ILs(i)%abundance .gt. 0) oiCELabund = oiCELabund + ILs(i)%abundance/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          if (ILs(i)%abundance .gt. 0) then
            weight = weight + 1/  ((ILs(i)%int_err/ ILs(i)%intensity)  **2)
          endif
        enddo
        if (weight .gt. 0) then
          oiCELabund = oiCELabund / weight
        else
          oiCELabund = 0.0
        endif        

        NOabundCEL = oiCELabund

! now get abundances for ORLs
! o2+
       call oii_rec_lines(oiiiTemp,oiiDens,DBLE(1),oiiRLs)

       do i = 1,329 
         do j = 1,415 
          if (abs(linelist_dered(1,i)-oiiRLs(j)%Wave) .le. 0.005) then
            oiiRLs(j)%Obs = linelist_dered(2,i) 
            oiiRLs(j)%abundance = oiiRLs(j)%obs/oiiRLs(j)%Int 
          endif
         enddo
       enddo

!N2+

       call nii_rec_lines(oiiiTemp,oiiDens,DBLE(1),niiRLs)

       do i = 1,329
         do j = 1,99
          if (abs(linelist_dered(1,i)-niiRLs(j)%Wave) .le. 0.005) then
            niiRLs(j)%Obs = linelist_dered(2,i)
            niiRLs(j)%abundance = niiRLs(j)%obs/niiRLs(j)%Int
          endif
         enddo
       enddo

!C2+ 
       call cii_rec_lines(oiiiTemp,oiiDens,DBLE(1),ciiRLs)

       do i = 1,329
         do j = 1,57
          if (abs(linelist_dered(1,i)-ciiRLs(j)%Wave) .le. 0.005) then
            ciiRLs(j)%Obs = linelist_dered(2,i)
            ciiRLs(j)%abundance = ciiRLs(j)%obs/ciiRLs(j)%Int
          endif
         enddo
       enddo

!Ne2+
       call neii_rec_lines(oiiiTemp,oiiDens,DBLE(1),neiiRLs)

       do i = 1,329 
         do j = 1,38
          if (abs(linelist_dered(1,i)-neiiRLs(j)%Wave) .le. 0.005) then
           neiiRLs(j)%Obs = linelist_dered(2,i)
            neiiRLs(j)%abundance = neiiRLs(j)%obs/neiiRLs(j)%Int
          endif
         enddo
       enddo

!C3+, N3+
       call xiii_rec_lines(oiiiTemp,oiiDens,DBLE(1),xiiiRLs)

       do i = 1,329
         do j = 1,6
          if (abs(linelist_dered(1,i)-xiiiRLs(j)%Wave) .le. 0.005) then
            xiiiRLs(j)%Obs = linelist_dered(2,i)
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

      do i = 1,57
        if (ciiRLs(i)%abundance .gt. 0) then 
          print "(1X,F7.2,1X,F6.3,1X,ES9.3)",ciiRLs(i)%wave,ciiRLs(i)%obs,ciiRLs(i)%abundance
          rlabundtemp = rlabundtemp + ciiRLs(i)%obs
          weight = weight + ciiRLs(i)%Int
        endif 
        if (ciiRLs(i)%wave .eq. 4267.15D0) then
          cii4267rlabund = ciiRLs(i)%abundance
        endif
      enddo

      if (weight .gt. 0) then
        ciirlabund = rlabundtemp/weight
        print *,""
        print *,"CII"
        print *,"lambda   Int   Abund"
        print "(A34,ES9.3)","Abundance (all lines co-added): ",ciirlabund 
        print "(A34,ES9.3)","Abundance (4267 line only):     ",cii4267rlabund
      else
        print *,"No CII recombination lines"
      endif

!nii recombination lines

!      print *,"lambda   Mult   Int   Abund"
      do i = 1,99
        if (niiRLs(i)%abundance .gt. 0) then
!          print "(F7.2,1X,A7,1X,F6.3,1X,ES9.3)",niiRLs(i)%wave,niiRLs(i)%Mult,niiRLs(i)%obs,niiRLs(i)%abundance
          rlabundtemp = rlabundtemp + niiRLs(i)%obs
          weight = weight + niiRLs(i)%Int
        endif
      enddo

  if (weight .gt. 0) then
      print *,""
      print *,"NII"

      print *,"Abundance from co-added intensity: "
      print "(ES9.2)",rlabundtemp / weight

      niimultiplets%Multiplet = (/"V3     ","V5     ","V8     " ,"V12    ","V20    ","V28    ","3d-4f  "/)

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

      if (isnan((rlabundtemp/weight))) then
      niimultiplets(7)%abundance = 0
      else
      niimultiplets(7)%abundance = rlabundtemp/weight
      endif

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
  else
    print *,"No NII recombination lines"
  endif

!oii recombination lines

      rlabundtemp = 0.00
      weight = 0.00

!      print *,"lambda   Mult   Int   Abund"
      do i = 1,415
        if (oiiRLs(i)%abundance .gt. 0) then
!          print "(F7.2,1X,A7,1X,F6.3,1X,ES9.3)",oiiRLs(i)%wave,oiiRLs(i)%Mult,oiiRLs(i)%obs,oiiRLs(i)%abundance
          rlabundtemp = rlabundtemp + oiiRLs(i)%obs
          weight = weight + oiiRLs(i)%Int
        endif
      enddo

  if (weight .gt. 0) then

      print *,""
      print *,"OII"


      print *,"Abundance from co-added intensity: "
      print "(ES9.2)",rlabundtemp/weight

      oiimultiplets%Multiplet = (/" V1    "," V2    "," V5    " ," V10   "," V11   "," V12   "," V19   "," V20   "," V25   "," V28   "," V33   "," 3d-4f "/)

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

      if ( isnan( (rlabundtemp/weight) ) ) then
      oiimultiplets(12)%abundance = 0
      else
      oiimultiplets(12)%abundance = rlabundtemp/weight
      endif

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
   else
      print *,"No OII recombination lines"
   endif

!neii recombination lines

      rlabundtemp = 0.0
      weight = 0.0

!      print *,"lambda   Mult   Int   Abund"
      do i = 1,38
        if (neiiRLs(i)%abundance .gt. 0) then
!          print "(F7.2,1X,F6.3,1X,ES9.3)",neiiRLs(i)%wave,neiiRLs(i)%obs,neiiRLs(i)%abundance 
           rlabundtemp = rlabundtemp + neiiRLs(i)%obs
           weight = weight + neiiRLs(i)%Int 
        endif
      enddo

   if (weight .gt. 0) then
      print *,""
      print *,"NeII"

      neiiRLabund = rlabundtemp/weight

      print *,"Abundance from co-added intensities: "
      print "(ES9.3)",neiiRLabund
   else
      print *,"No NeII recombination lines"
   endif

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
        print *,"No CIII recombination lines"
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
        print *,"No NIII recombination lines"
     endif

! ICFs (Kingsburgh + Barlow 1994)

! oxygen, just the simple case at the moment
        !if (oiiCELabund .gt. 0 .and. oiiiCELabund .gt. 0 .and. oivCELabund .eq. 0 .and. nvCELabund .eq. 0)then !no O3+ or N4+ seen
                     CELicfO = ((heiabund + heiiabund)/heiabund)**(2./3.) !KB94 A9
                    OabundCEL = CELicfO * (oiiCELabund + oiiiCELabund) !A10
        !elseif (oiiCELabund .gt. 0 .and. oiiiCELabund .gt. 0 .and. oivCELabund .gt. 0 .and. nvCELabund .eq. 0) then !no N4+ seen
                !CELicfO = 1
                !OabundCEL = oiiCELabund + oiiiCELabund + oivCELabund !sum of visible ionisation stages
        !endif



 !three other cases to add for oxygen, A6, A8, and sum of observed abundances, commented case causes errors due to undefined variables

! nitrogen - complete

     if (niiCELabund .gt. 0 .and. niiiUVCELabund .gt. 0 .and. nivCELabund .gt. 0) then !all ionisation stages seen
       CELicfN = 1.
       NabundCEL = niiCELabund + niiiCELabund + nivCELabund
     elseif (niiCELabund .gt. 0 .and. niiiUVCELabund .eq. 0 .and. nivCELabund .gt. 0) then !no N2+ seen
       CELicfN = 1.5
       NabundCEL = 1.5*(niiCELabund + nivCELabund)
     elseif (niiCELabund .gt. 0 .and. niiiUVCELabund .eq. 0 .and. nivCELabund .eq.0) then !Only N+ seen
       CELicfN = OabundCEL/oiiCELabund
       NabundCEL = niiCELabund * CELicfN 
     endif

! carbon - couple of cases still to add.

     if (ciiCELabund .gt. 0 .and. ciiiCELabund .gt. 0 .and. civCELabund .gt. 0 .and. heiiabund .eq. 0) then !No C4+ but all other seen
       CELicfC = 1.
       CabundCEL = ciiCELabund + ciiiCELabund + civCELabund
     elseif (ciiCELabund .eq. 0 .and. ciiiCELabund .gt. 0 .and. civCELabund .gt. 0 .and. heiiabund .eq. 0) then !No C4+, and no CII lines seen
       CELicfC = (oiiCELabund + oiiiCELabund) / oiiiCELabund !A11
       CabundCEL = CELicfC * (ciiCELabund + ciiiCELabund) !A12
     elseif (ciiCELabund .eq. 0 .and. ciiiCELabund .gt. 0 .and. civCELabund .eq. 0 .and. oiiiCELabund .gt. 0) then !Only C2+ seen, but O2+ also seen
       CELicfC = OabundCEL/oiiiCELabund !A13
       CabundCEL = CELicfC * ciiiCELabund !A14
     elseif (nvCELabund .gt. 0 .and. heiiabund .gt. 0) then !N4+ and He2+
       CELicfC = 1/(1-(nvCELabund/(niiCELabund + niiiCELabund + nivCELabund + nvCELabund)))!A16
       if (CELicfC .gt. 5) then !High-excitation case
         CELicfC = (niiCELabund + niiiCELabund + nivCELabund + nvCELabund) / (niiCELabund + niiiCELabund + nivCELabund) !A18
       endif
       CabundCEL = CELicfC * (ciiCELabund + ciiiCELabund + civCELabund) !A19
     elseif (ciiCELabund .gt. 0 .and. ciiiCELabund .gt. 0 .and. civCELabund .gt. 0 .and. heiiabund .gt. 0 .and. nvCELabund .eq. 0) then !PN is hot enough for He2+ but not for N4+
        CELicfC = ((heiiabund + heiabund)*heiabund)**(1./3.) !A20
        CabundCEL = CELicfC * (ciiCELabund + ciiiCELabund + civCELabund) !A21
     elseif (ciiCELabund .eq. 0 .and. ciiiCELabund .gt. 0 .and. civCELabund .gt. 0) then !C+ not seen
        CELicfC = 1/(1 - (niiCELabund/NabundCEL) - 2.7*(nvCELabund/NabundCEL)) !A22
        if (CELicfC .gt. 5) then !if ICF is greater than 5
           CELicfC = (niiCELabund + niiiCELabund + nivCELabund + nvCELabund)/(niiiCELabund + nivCELabund) !A23
        endif
        CabundCEL = CELicfC * (ciiiCELabund+ civCELabund) !A24
     elseif (ciiCELabund .gt. 0 .and. ciiiCELabund .gt. 0 .and. civCELabund .gt. 0 .and. heiiabund .gt. 0 .and. (nivCELabund .eq. 0 .or. nvCELabund .eq. 0)) then !final case i think.
        CELicfC = ((oiiCELabund + oiiiCELabund)/oiiiCELabund)*((heiiabund + heiabund)*heiabund)**(1./3.) !A25
        CabundCEL = CELicfC * (ciiCELabund + ciiiCELabund + civCELabund) !A26
     endif 


! Neon - complete

     if (neiiiCELabund .gt. 0 .and. neivCELabund .gt. 0 .and. nevCELabund .gt. 0) then !all stages seen
       CELicfNe = 1.
       NeabundCEL = neiiiCELabund + neivCELabund + nevCELabund
     elseif (neiiiCELabund .gt. 0 .and. neivCELabund .eq. 0 .and. nevCELabund .gt. 0) then !no Ne IV seen
       CELicfNe = 1.5
       NeabundCEL = CELicfNe * (neiiiCELabund + nevCELabund) !KB94 A27
     elseif (neiiiCELabund .gt. 0 .and. neivCELabund .eq. 0 .and. nevCELabund .eq. 0) then !Only Ne2+ seen
       CELicfNe = OabundCEL / oiiiCELabund !KB94 A28
       NeabundCEL = CELicfNe * neiiiCELabund
     endif

! Argon - complete

     if (ariiiCELabund .gt. 0 .and. arivCELabund .eq. 0 .and. arvCELabund .eq. 0) then !only Ar2+ seen
       CELicfAr = 1.87 !KB94 A32
       ArabundCEL = CELicfAr * ariiiCELabund !KB94 A33
     elseif (ariiiCELabund .eq. 0 .and. arivCELabund .gt. 0 .and. arvCELabund .eq. 0) then !Only Ar3+ seen
       CELicfAr = NeabundCEL / neiiiCELabund !KB94 A34
       ArabundCEL = CELicfAr * arivCELabund !KB94 A35
     elseif (ariiiCELabund .gt. 0 .and. arivCELabund .gt. 0) then !Ar 2+ and 3+ seen
       CELicfAr = 1.0 !??
       ArabundCEL = ariiiCELabund + arivCELabund + arvCELabund !KB94 A31
     endif

! Sulphur

     if (siiCELabund .gt. 0 .and. siiiCELabund .gt. 0 .and. sivIRCELabund .eq. 0) then !both S+ and S2+
       CELicfS = (1 - (  (1-(oiiCELabund/OabundCEL))**3.0  )   )**(-1.0/3.0) !KB94 A36

       SabundCEL = CELicfS * (siiCELabund + siiiCELabund) !KB94 A37
     elseif (siiCELabund .gt. 0 .and. siiiCELabund .gt. 0 .and. sivIRCELabund .gt. 0) then !all states observed
       CELicfS = 1.
       SabundCEL = siiCELabund + siiiCELabund + sivIRCELabund
     elseif (siiCELabund .gt. 0 .and. siiiCELabund .eq. 0 .and. sivIRCELabund .eq. 0) then !Only S+ observed
       CELicfS = (((oiiiCELabund/oiiCELabund)**0.433) * 4.677) * (1-(1-((oiiCELabund/OabundCEL)**3)))**(-1./3.) ! KB94 A37 with S2+/S+ from A38
       SabundCEL = CELicfS * siiCELabund
     endif

!XXXX finish very high excitation case. i.e A39-A41 for C and O

!ORLs
!Oxygen

     if (oiiRLabund .gt. 0) then
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

     if (oiiRLabund .gt. 0) then
       RLicfNe = OabundRL / oiiRLabund
       NeabundRL = RLicfNe * neiiRLabund
     else
       RLicfNe = 1.0
       NeabundRL = 0.0
     endif

!finish these XXXX


!Printout edited to include weighted averages of Ionic species as these are neccessary for paper tables. DJS

print *,""
print *,"Total abundances"
print *,"================"
print *,""

print *,"CELs"
print *,""
print *,"Element           ICF     X/H"
print *,"-------           ---     ---"

if(NCabundCEL > 0) print "(A15,F5.2,4X,ES8.2,2X,F5.2)"," Neutral Carbon ",0.0,NCabundCEL,12+log10(NCabundCEL)
        if(runonce == 0 .and. NCabundCEL > 0) WRITE(unit=841,FMT=*)NCabundCEL
if(CabundCEL > 0)  print "(A15,F5.2,4X,ES8.2,2X,F5.2)"," Carbon         ",CELicfC,CabundCEL,12+log10(CabundCEL)
        if(runonce == 0 .and. CabundCEL > 0) WRITE(unit = 842,FMT=*)CabundCEL

if(niiCELabund > 0) print "(A24,ES8.2)"                ,"  N+/H+                 ",niiCELabund
        if(runonce == 0 .and. niiCELabund > 0) WRITE(unit = 843,FMT=*)niiCELabund
if(NabundCEL > 0)  print "(A15,F5.2,4X,ES8.2,2X,F5.2)"," Nitrogen       ",CELicfN,NabundCEL,12+log10(NabundCEL)
        if(runonce == 0 .and. NabundCEL > 0) WRITE(unit = 844,FMT=*)NabundCEL

if(NOabundCEL > 0) print "(A15,F5.2,4X,ES8.2,2X,F5.2)"," Neutral Oxygen ",0.0,NOabundCEL,12+log10(NOabundCEL)
        if(runonce == 0 .and. NOabundCEL > 0) WRITE(unit = 845,FMT=*)NOabundCEL
if(oiiCELabund >0) print "(A24,ES8.2)"                , "  O+/H+                 ",oiiCELabund
        if(runonce == 0 .and. oiiCELabund >0) WRITE(unit = 846,FMT=*)oiiCELabund
if(oiiiCELabund >0) print "(A24,ES8.2)"               , "  O++/H+                ",oiiiCELabund
        if(runonce == 0 .and. oiiiCELabund >0) WRITE(unit = 847,FMT=*)oiiiCELabund
if(OabundCEL > 0)  print "(A15,F5.2,4X,ES8.2,2X,F5.2)"," Oxygen         ",CELicfO,OabundCEL,12+log10(OabundCEL)
        if(runonce == 0 .and. OabundCEL > 0) WRITE(unit = 848,FMT=*)OabundCEL

if(neiiiCELabund >0) print "(A24,ES8.2)"               , "  Ne++/H+               ",neiiiCELabund
        if(runonce == 0 .and. neiiiCELabund >0) WRITE(unit = 849,FMT=*)neiiiCELabund
if(neivCELabund >0)  print "(A24,ES8.2)"               , "  Ne+++/H+              ",neivCELabund
        if(runonce == 0 .and. neivCELabund >0) WRITE(unit = 850,FMT=*)neivCELabund
if(nevCELabund >0)   print "(A24,ES8.2)"               , "  Ne++++/H+             ",nevCELabund
        if(runonce == 0 .and. nevCELabund >0) WRITE(unit = 851,FMT=*)nevCELabund
if(NeabundCEL > 0) print "(A15,F5.2,4X,ES8.2,2X,F5.2)"," Neon           ",CELicfNe,NeabundCEL,12+log10(NeabundCEL)
        if(runonce == 0 .and. NeabundCEL > 0) WRITE(unit = 852,FMT=*)NeabundCEL

if(ariiiCELabund >0) print "(A24,ES8.2)"               , "  Ar++/H+               ",ariiiCELabund
        if(runonce == 0 .and. ariiiCELabund >0) WRITE(unit = 853,FMT=*)ariiiCELabund
if(arivCELabund >0)  print "(A24,ES8.2)"               , "  Ar+++/H+              ",arivCELabund
        if(runonce == 0 .and. arivCELabund >0) WRITE(unit = 854,FMT=*)arivCELabund
if(arvCELabund >0)   print "(A24,ES8.2)"               , "  Ar++++/H+             ",arvCELabund
        if(runonce == 0 .and. arvCELabund >0) WRITE(unit = 855,FMT=*)arvCELabund
if(ArabundCEL > 0) print "(A15,F5.2,4X,ES8.2,2X,F5.2)"," Argon          ",CELicfAr,ArabundCEL,12+log10(ArabundCEL)
        if(runonce == 0 .and. ArabundCEL > 0) WRITE(unit = 856,FMT=*)ArabundCEL

if(siiCELabund >0) print "(A24,ES8.2)"                , "  S+/H+                 ",siiCELabund
        if(runonce == 0 .AND. siiCELabund >0) WRITE(unit = 857,FMT=*)siiCELabund
if(siiiCELabund >0) print "(A24,ES8.2)"               , "  S++/H+                ",siiiCELabund
        if(runonce == 0 .AND. siiiCELabund >0) WRITE(unit = 858,FMT=*)siiiCELabund
if(SabundCEL > 0)  print "(A15,F5.2,4X,ES8.2,2X,F5.2)"," Sulphur        ",CELicfS,SabundCEL,12+log10(SabundCEL)
        if(runonce == 0 .AND. SabundCEL > 0) WRITE(unit = 859,FMT=*)SabundCEL

if(OabundCEL > 0 .and. NabundCEL > 0)  print*, " "
if(OabundCEL > 0 .and. NabundCEL > 0)  print*, " "
if(OabundCEL > 0 .and. NabundCEL > 0)  print "(A15,F5.2,4X,ES8.2,2X,F5.2)"," N/O            ", log10(NabundCEL/OabundCEL)
if(NCabundCEL > 0 .and. NOabundCEL > 0) print*, " "
if(NCabundCEL > 0 .and. NOabundCEL > 0) print "(A15,F5.2,4X,ES8.2,2X,F5.2)"," nC/nO            ", log10(NCabundCEL/NOabundCEL)


print *,""
print *,"ORLs"
print *,""
print *,"Element     ICF     X/H"
print *,"-------     ---     ---"
if(Hetotabund > 0) print "(A12,F5.2,4X,ES8.2,2X,F5.2)"," Helium      ",1.0,Hetotabund,12+log10(Hetotabund)
        if(runonce == 0 .AND. Hetotabund > 0) WRITE(unit = 860,FMT=*)Hetotabund
if(CabundRL > 0) print "(A12,F5.2,4X,ES8.2,2X,F5.2)"," Carbon      ",RLicfC,CabundRL,12+log10(CabundRL)
        if(runonce == 0 .AND. CabundRL > 0) WRITE(unit = 861,FMT=*)CabundRL
if(NabundRL > 0) print "(A12,F5.2,4X,ES8.2,2X,F5.2)"," Nitrogen    ",RLicfN,NabundRL,12+log10(NabundRL)
        if(runonce == 0 .AND. NabundRL > 0) WRITE(unit = 862,FMT=*)NabundRL
if(OabundRL > 0) print "(A12,F5.2,4X,ES8.2,2X,F5.2)"," Oxygen      ",RLicfO,OabundRL,12+log10(OabundRL)
        if(runonce == 0 .AND. OabundRL > 0) WRITE(unit = 863,FMT=*)OabundRL
if(NeabundRL > 0) print "(A12,F5.2,4X,ES8.2,2X,F5.2)"," Neon        ",RLicfNe,NeabundRL,12+log10(NeabundRL)
        if(runonce == 0 .AND. NeabundRL > 0) WRITE(unit = 864,FMT=*)NeabundRL


!Strong line methods

ion_no1 = get_ion("oiii4959   ",ILs, Iint)
ion_no2 = get_ion("oiii5007   ",ILs, Iint)
ion_no3 = get_ion("oii3726    ",ILs, Iint)
ion_no4 = get_ion("oii3729    ",ILs, Iint)

        print *,""
        print *,"O/H (strong line methods)"
        print *,"-------------------"
        print *,"Calibration         Reference          O/H"

if (ion_no1 .gt. 0 .and. ion_no2 .gt. 0) then
  X23 = log10((ILs(ion_no3)%int_dered + ILs(ion_no4)%int_dered)/(ILs(ion_no1)%int_dered + ILs(ion_no2)%int_dered))
  O_R23upper = 9.50 - (1.4 * X23)
  O_R23lower = 6.53 + (1.45 * X23)
  print "(1X,A37,2X,F5.2)","R23 (upper branch)  Pilyugin 2000    ",O_R23upper
  print "(1X,A37,2X,F5.2)","R23 (lower branch)  Pilyugin 2000    ",O_R23lower
endif

ion_no1 = get_ion("nii6584    ",ILs, Iint)
if (ILs(ion_no1)%int_dered .gt. 0 .and. H_BS(1)%int_dered .gt. 0) then
  N2 = ILs(ion_no1)%int_dered / H_BS(1)%int_dered
  O_N2 = 8.90 + (0.57 * N2)
  print "(1X,A37,2X,F5.2)","N2                  Pet. + Pag. 2004 ",O_N2
endif

ion_no2 = get_ion("oiii5007   ",ILs, Iint)
if (ILS(ion_no1)%int_dered .gt. 0 .and. ILs(ion_no2)%int_dered .gt. 0) then
  O3N2 = log10((ILS(ion_no2)%int_dered*H_BS(1)%int_dered)/(ILS(ion_no1)%int_dered * H_BS(2)%int_dered))
  O_O3N2 = 8.73 - (0.32*O3N2)
  print "(1X,A37,2X,F5.2)","O3N2                Pet. + Pag. 2004 ",O_O3N2
endif

ion_no1 = get_ion("ariii7135  ",ILs, Iint)
ion_no2 = get_ion("oiii5007   ",ILs, Iint)
if (ILs(ion_no1)%int_dered .gt. 0 .and. ILs(ion_no2)%int_dered .gt. 0) then
  Ar3O3 = ILs(ion_no1)%int_dered / ILs(ion_no2)%int_dered
  O_Ar3O3 = 8.91 + (0.34*Ar3O3) + (0.27*Ar3O3**2) + (0.2*Ar3O3**3)
  print "(1X,A37,2X,F5.2)","Ar3O3               Stasinska 2006   ",O_Ar3O3
endif

ion_no1 = get_ion("siii9069   ",ILs, Iint)
if (ILs(ion_no1)%int_dered .gt. 0 .and. ILs(ion_no2)%int_dered .gt. 0) then
  S3O3 = ILs(ion_no1)%int_dered / ILs(ion_no2)%int_dered
  O_S3O3 = 9.37 + (2.03*S3O3) + (1.26*S3O3**2) + (0.32*S3O3**3)
  print "(1X,A37,2X,F5.2)","S3O3                Stasinska 2006   ",O_S3O3
endif

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

if(adfo2plus >0) print "(A12,F5.2)","adf (O2+) = ", adfo2plus
if(adfO>0) print "(A12,F5.2)","adf (O)   = ", adfO
if(adfn2plus>0) print *,""
if(adfn2plus>0) print "(A12,F5.2)","adf (N2+) = ", adfn2plus
if(adfn>0) print "(A12,F5.2)","adf (N)   = ", adfn
if(adfc2plus>0) print *,""
if(adfc2plus>0) print "(A12,F5.2)","adf (C2+) = ", adfc2plus
if(adfc>0) print "(A12,F5.2)","adf (C)   = ", adfc
if(adfne2plus>0) print *,""
if(adfne2plus>0) print "(A12,F5.2)","adf (Ne2+)= ", adfne2plus
if(adfne>0) print "(A12,F5.2)","adf (Ne)  = ", adfne




contains

        SUBROUTINE get_diag(name1, name2, lines, diag)
                TYPE(line), DIMENSION(51), INTENT(IN) :: lines 
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
                TYPE(line), DIMENSION(51), INTENT(IN) :: lines 
                CHARACTER*11 :: name1, name2, name3
                INTEGER :: ion_no1, ion_no2, ion_no3
                DOUBLE PRECISION :: diag, factor1, factor2, ratio, ratio2


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
