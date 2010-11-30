        print *,""
        print *,"Diagnostics"
        print *,"-----------"

        ! start with Oii density, assuming Te=10,000K
        ! note for roger, this now works by having a master structure with all of the lines you said were interesting, whenever you want a line strength from it you can call this function with the line name and it will return the line number in the structure which you can then use..

        ion_no1 = get_ion("oii3729    ", ILs)
        ion_no2 = get_ion("oii3726    ", ILs)

        if((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0))then
                oiiNratio = DBLE(ILs(ion_no1)%int_dered) / DBLE(ILs(ion_no2)%int_dered)
                call get_diagnostic("oii       ","1,2/                ","1,3/                ",oiiNratio,"D",DBLE(10000.0),result)
                oiiDens = result
        else
                oiiNratio = 0
                print *,"No [OII] density"
                print *,"Enter a density value: "
                read (5,"(F7.0)") oiiDens
                stop
        endif

        ![oiii] temperature, using either oii density or user input
                
        ion_no1 = get_ion("oiii5007   ", ILs)
        ion_no2 = get_ion("oiii4959   ", ILs)
        ion_no3 = get_ion("oiii4363   ", ILs)

        if(((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0)) .and. (ILs(ion_no3)%int_dered .gt. 0 ))then
                
                oiiiTratio = (ILs(ion_no1)%int_dered + ILs(ion_no2)%int_dered) / ILs(ion_no3)%int_dered
        elseif(((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .eq. 0)) .and. (ILs(ion_no3)%int_dered .gt. 0 ))then
                ILs(ion_no2)%int_dered = ILs(ion_no1)%int_dered * 1.3356
                oiiiTratio = (ILs(ion_no1)%int_dered * 1.3356) / ILs(ion_no3)%int_dered
        elseif(((ILs(ion_no1)%int_dered .eq. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0)) .and. (ILs(ion_no3)%int_dered .gt. 0 ))then
                ILs(ion_no1)%int_dered = ILs(ion_no2)%int_dered * 3.98
                oiiiTratio = (ILs(ion_no2)%int_dered * 3.98) / ILs(ion_no3)%int_dered
        else
                print *,"Can't calculate [O III] temperature"
                stop
        endif

! quick iteration to get oii density and oiii temperature
        
        call get_diagnostic("oiii      ","2,4,3,4/            ","4,5/                ",oiiiTratio,"T",oiiDens,oiiiTemp)

!        print *, "[O III] temperature (1st iteration) = ", oiiiTemp

      if (oiiNratio .gt. 0) then
        call get_diagnostic("oii       ","1,2/                ","1,3/                ",oiiNratio,"D",oiiiTemp,oiiDens)
        print "(A20,F7.0)","[O II] density = ",oiiDens
        count = 1

        call get_diagnostic("oiii      ","2,4,3,4/            ","4,5/                ",oiiiTratio,"T",oiiDens,oiiiTemp)
        print "(A20,F7.0)","[O III] temperature = ",oiiiTemp
      endif
        !dunno what this commented out bit does..
        !conditions = getTD(DBLE(8000), DBLE(0), DBLE(10), DBLE(10), oiiNratio, oiiiTratio, "oii", "oiii", -1) result(TD)

        !oiii IR density

        ion_no1 = get_ion("oiii52um   ", ILs)
        ion_no2 = get_ion("oiii88um   ", ILs)
        
        if (((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0)))then
                oiiiIRNratio = (ILs(ion_no1)%int_dered) / (ILs(ion_no2)%int_dered)
                call get_diagnostic("oiii       ","2,3/                ","1,2/                ",oiiiIRNratio,"D",oiiiTemp,result)
                print *,"[O III infrared density: ",result
        else
                print *,"[O III] IR density = n/a"
        endif

! [N II] temperature

        ion_no1 = get_ion("nii6548    ", ILs)
        ion_no2 = get_ion("nii6584    ", ILs)
        ion_no3 = get_ion("nii5754    ", ILs)

        if(((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0)) .and. (ILs(ion_no3)%int_dered .gt. 0 ))then
                niiTratio = (ILs(ion_no1)%int_dered + ILs(ion_no2)%int_dered) / ILs(ion_no3)%int_dered
                call get_diagnostic("nii       ","2,4,3,4/            ","4,5/                ",niiTratio,"T",oiiDens,result)
                print "(A20,F7.0)","[N II] temperature = ",result
        endif

! [O II] temperature

        ion_no1 = get_ion("oii7319    ",ILs)
        ion_no2 = get_ion("oii7330    ",ILs)
        ion_no3 = get_ion("oii3726    ",ILs)
        ion_no4 = get_ion("oii3729    ",ILs)

        if (ion_no1 .gt. 0 .and. ion_no2 .gt. 0 .and. ion_no3 .gt. 0 .and. ion_no4 .gt. 0) then
             oiiTratio = (ILs(ion_no1)%int_dered+ILs(ion_no2)%int_dered)/(ILs(ion_no3)%int_dered+ILs(ion_no4)%int_dered)
             call get_diagnostic("oii       ","2,4,2,5,3,4,3,5/    ","1,2,1,3/            ",oiiTratio,"T",oiiDens,oiiTemp)
             print "(A20,F7.0)","[O II] temperature = ",oiiTemp
         else
            print *,"[O II] temperature = n/a"
         endif

! Various densities
! Argon
        ion_no1 = get_ion("ariv4740   ", ILs)
        ion_no2 = get_ion("ariv4711   ", ILs)

        if(((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0))) then
                arivNratio = ILs(ion_no1)%int_dered / ILs(ion_no2)%int_dered
                call get_diagnostic("ariv      ","1,2/                ","1,3/                ",arivNratio,"D",oiiiTemp, arivDens)
                print "(A20,F7.0)","[Ar IV] density = ",arivDens
                count = count + 1
        endif
! Chlorine
        ion_no1 = get_ion("cliii5537  ", ILs)
        ion_no2 = get_ion("cliii5517  ", ILs)

        if(((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0))) then
                cliiiNratio = ILs(ion_no1)%int_dered / ILs(ion_no2)%int_dered
                call get_diagnostic("cliii     ","1,2/                ","1,3/                ",cliiiNratio,"D",oiiiTemp, cliiiDens)
                print "(A20,F7.0)","[Cl III] density = ",cliiiDens
                count = count + 1
        endif

!Sulphur
        ion_no1 = get_ion("sii6731    ", ILs)
        ion_no2 = get_ion("sii6716    ", ILs)

        if(((ILs(ion_no1)%int_dered .gt. 0) .AND. (ILs(ion_no2)%int_dered .gt. 0))) then
                siiNratio = ILs(ion_no1)%int_dered / ILs(ion_no2)%int_dered 
                call get_diagnostic("sii       ","1,2/                ","1,3/                ",siiNratio,"D",oiiiTemp, siiDens)
                print "(A20,F7.0)","[S II] density = ",siiDens
                count = count + 1
        endif

! [SII] temperature

        ion_no3 = get_ion("sii4068    ",ILs)
        ion_no4 = get_ion("sii4076    ",ILs)

        if (siiNratio .gt. 0 .and. ILs(ion_no3)%int_dered .gt. 0 .and. ILs(ion_no4)%int_dered .gt. 0) then
            siiTratio = (ILs(ion_no1)%int_dered+ILs(ion_no2)%int_dered)/(ILs(ion_no3)%int_dered+ILs(ion_no4)%int_dered)
            call get_diagnostic("sii       ","1,2,1,3/            ","1,4,1,5/            ",siiTratio, "T",oiiDens,siiTemp)
            print "(A20,F7.0)","[S II] temperature = ",siiTemp
        else
            print *,"[S II] temperature = n/a"
        endif

![Ne III] temperature

        ion_no1 = get_ion("neiii3868  ",ILs)
        ion_no2 = get_ion("neiii3967  ",ILs)
        ion_no3 = get_ion("neiii15p5um",ILs)

        if (ILs(ion_no1)%int_dered .gt. 0 .and. ILs(ion_no2)%int_dered .gt. 0 .and. ILs(ion_no3)%int_dered .gt. 0) then
           neiiiTratio = (ILs(ion_no1)%int_dered+ILs(ion_no2)%int_dered)/ILs(ion_no3)%int_dered
           call get_diagnostic("neiii     ","1,4,2,4/            ","1,2/                ",neiiiTratio,"T",oiiiTemp,neiiiTemp)
           print "(A20,F7.0)","[Ne III] temperature = ",neiiiTemp
        else
           print *,"[Ne III] temperature = n/a"
        endif

! mean density

        meandensity = (oiiDens + arivDens + siiDens + cliiiDens) / count
        print "(A14,F7.0)","Mean density: ",meandensity 
