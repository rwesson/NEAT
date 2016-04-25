! For input atomic data file formats see Readme file attached at the end
! **********************************************************************
!          Program EQUIB  (FORTRAN 77)
!
!    Programming history:
!
!    1981 May 3    IDH    Version 1
!    1981 May 5    IDH    Minibug fixed!
!    1981 May 7    IDH    Now takes collision rates or strengths
!    1981 Aug 3    SA     Interpolates collision strengths
!    1981 Aug 7    SA     Input method changed
!    1984 Nov 19   RESC   SA files entombed in scratch disk. Logical
!                         filenames given to SA's data files.
!    1995 Aug      DPR    Changed input file format. Increased matrices.
!    1996 Feb      XWL    Tidy up. SUBROUTINES SPLMAT, HGEN, CFY and CFD
!                         modified such that matrix sizes (i.e. maximum
!                         of Te and maximum no of levels) can now be cha
!                         by modifying the parameters NDIM1, NDIM2 and N
!                         in the Main program. EASY!
!                         Now takes collision rates as well.
!                         All variables are declared explicitly
!                         Generate two extra files (ionpop.lis and ionra
!                         of plain stream format for plotting
!    1996 June     CJP    Changed input data format for cases IBIG=1,2.
!                         Fixed readin bug for IBIG=2 case.
!                         Now reads reformatted upsilons (easier to see
!                         and the 0 0 0 data end is excluded for these c
!                         The A values have a different format for IBIG=
!    2009 April    RW     Converted to F90, inputs from cmd line, version
!                         written purely to do diagnostics.
!    2012 February RW     Tidying up, combining into simpler routines for NEAT
! ***** N.B!!  NO TRAPS FOR BAD DATA!!  TAKE CARE!! ****C
!
      module mod_equib
      contains

      subroutine get_diagnostic(ion,levu,levl,inratio,diagtype,fixedq,result,ndim2,ndim1,atomicdata,iion)
      use mod_atomicdata
      implicit none
      integer, parameter :: dp = kind(1.d0)

      integer :: NDIM1, NDIM2, NDIM1T3, MAXND
                                                      !Maximum no of Te & levels
                                             !NDIM1T3 should be at least 3*NDIM1
                                                   !Maximum no. of Ne increments
      parameter (MAXND=100)
      integer :: G(NDIM2),                                                 &
     &  ITRANA(2,NDIM2),ITRANB(2,NDIM2),ITRANC(2,NDIM2),LOOP
      type(atomic_data),dimension(:),intent(in) :: atomicdata
      integer :: iion,nion
      real(kind=dp) :: N(NDIM2)
      real(kind=dp) :: & !TDRAT(2,MAXND)
     & TNIJ(NDIM2,NDIM2), FINTIJ(NDIM2,NDIM2),                          &
     & WAVA(NDIM2), WAVB(NDIM2), WAVC(NDIM2), CS(NDIM2,NDIM2),          &
     & QEFF(NDIM2,NDIM2), QQ(NDIM1),                                    &
     & QOM(NDIM1,NDIM2,NDIM2), A(NDIM2,NDIM2), E(NDIM2), T(NDIM1),      &
     & ROOTT(NDIM1), X(NDIM2,NDIM2), Y(NDIM2),                          &
     & X2(NDIM2,NDIM2), XKEEP(NDIM2,NDIM2), Y2(NDIM2), YKEEP(NDIM2),    &
     & HMH(NDIM1,NDIM1), D(NDIM1)
      character(len=20) :: LABEL(NDIM2)
      character(len=10) :: ION
      integer :: I, I1, I2, J, KK, LL, JT, JJD,                            &
     & NLEV, NTEMP, IBIG, IRATS,                                        &
     & NLEV1, INT, IND, IOPT, IT, IM1, JM1, IP1,                        &
     & IAPR, IBPR, ICPR, IKT, IA, IB, IC, IA1, IA2, IB1, IB2, IC1, IC2
      real(kind=dp) :: TEMPI, TINC, DENSI, DINC, DENS, DLOGD, TEMP, TLOGT,        &
     & TEMP2, DD, DELTEK, EXPE, VALUE, SUMN, TTT, TTP, AHB, EJI, WAV,   &
     & RLINT, FINT, SUMA, SUMB, FRAT, DEE

      real(kind=dp) :: fixedq
      real(kind=dp) :: inratio,result
      character(len=20) :: levu,levl
      character(len=1) :: diagtype
      real(kind=dp), DIMENSION(:,:), ALLOCATABLE :: RESULTS
      real(kind=dp) :: valtest(3)
      integer :: test

!debugging
#ifdef CO
        print *,"subroutine: get_diagnostic. ion=",ion
#endif

! check if ratio is meaningful

      if (inratio .le. 0.d0 .or. inratio .gt. 1.e10) then
        result=0.d0
        return
      endif

      ndim1t3=3*ndim1
      g=0
      itrana=0
      itranb=0
      itranc=0
      valtest=0
      test=0

      read(levu,*) ((ITRANA(LL,KK),LL=1,2),KK=1,ndim2)!150)
      read(levl,*) ((ITRANB(LL,KK),LL=1,2),KK=1,ndim2)!150)


!Transfer atomic data to local variables
      nion = 0
      do i = 1,iion
         if(atomicdata(i)%ion .eq. ion) nion=i
      enddo
      if (nion .eq. 0) then
          print *,"I'm afraid. I'm afraid, Dave."
          nion = 1
      endif
!          print*,nion,atomicdata(nion)%ion,ion
      nlev=atomicdata(nion)%nlevs
      ntemp=atomicdata(nion)%ntemps
          T(1:ntemp)= log10(atomicdata(nion)%Temps(1:ntemp))
      ROOTT(1:ntemp)=atomicdata(nion)%rootT(1:ntemp)
      Label(1:nlev)=atomicdata(nion)%labels(1:nlev)
      QOM(1:ntemp,1:nlev,1:nlev)=atomicdata(nion)%col_str(1:ntemp,1:nlev,1:nlev)
      A(1:nlev,1:nlev)=atomicdata(nion)%A_coeffs(1:nlev,1:nlev)
      E(1:nlev)=atomicdata(nion)%Waveno(1:nlev)
      G(1:nlev)=atomicdata(nion)%G(1:nlev)

      irats=0
      ibig=0

      NLEV1 = NLEV - 1
      ITRANC = 0

      !*****LOOP STARTS HERE*************************
      do LOOP = 1, 9
      if (diagtype .eq. "t" .or. diagtype .eq. "T") then

        if (LOOP .eq. 1) then
                TEMPI=5000.
        else
                TEMPI= valtest(1)
        endif

        INT=4
        TINC=(30000.)/((INT-1)**(LOOP))

        densi=fixedq
        dinc=0
        ind=1

        allocate(RESULTS(3,INT))

      else

        if (LOOP .eq. 1) then
          densi=0
        else
          densi=valtest(2)
        endif

        IND=4
        DINC=(100000.)/((IND-1)**(LOOP))

        TempI=fixedq
        TINC=0
        INT=1

        allocate(results(3,IND))
      endif

      if (densi .le. 0) densi=1
          if (tempi .lt. 5000) tempi=5000

                                                               !Start of Te loop
      do JT = 1, INT
        TEMP=TEMPI+(JT-1)*TINC
                                                               !Start of Ne loop

        do JJD = 1, IND
          DENS=DENSI+(JJD-1)*DINC
          if (TEMP.LE.0.D0.OR.DENS.LE.0.D0) then
            write (6,6100)
                print *,"Temp = ", TEMP, ", Dens = ", DENS, ", Ion = ",ion,diagtype
            stop
          endif
          DLOGD = LOG10 (DENS)
          TLOGT = LOG10 (TEMP)
          TEMP2= SQRT (TEMP)
                                                                 !Form matrices

          X=0.0
          CS=0.0
          QEFF=0.0
          TNIJ=0.0
          Y=0.0

          IOPT=0
          if (NTEMP.EQ.1) then
            write (6,*)
            write (6,*)                                                 &
     &      'Coll. strengths available for 1 Te only - assuming const'
          elseif (NTEMP.EQ.2) then
            write (6,*)
            write (6,*)                                                 &
     &      'Coll. strengths available for 2 Te only - linear interp'
          else
            call SPLMAT(T, NTEMP, IOPT, NDIM1, NDIM1T3, HMH)
            call CFD(TLOGT,T,NTEMP,NDIM1,HMH,D)
          endif
          do I = 2, NLEV
            do J = I, NLEV
                                                           !Negative!
              DELTEK = (E(I-1)-E(J))*1.4388463D0
              EXPE = EXP(DELTEK/TEMP)
              do IT = 1, NTEMP
                if (IRATS.EQ.0.D+00) then
                  QQ(IT) = QOM(IT,I-1,J)
                else
                                                      !Take out the exp. depend.
                  QQ(IT) = QOM(IT,I-1,J) / EXPE
                                                      !before interpolation
                endif
              enddo
              if (NTEMP.EQ.1) then
                 DD = QQ(1)
              elseif (NTEMP.EQ.2) then
                 DD = QQ(1) +                                           &
     &            (QQ(2) - QQ(1))/(T(2) - T(1)) * (TLOGT - T(1))
              else
                call CFY(TLOGT, DD, T, QQ, NTEMP, NDIM1, D)
              endif
              if (IRATS.EQ.0.D+00) then
                CS(I-1,J) = DD
              else
                CS(I-1,J) = DD * EXPE
              endif
              if (IRATS .EQ. 0.D+00) then
                QEFF(I-1,J) = 8.63D-06*CS(I-1,J) * EXPE / (G(I-1)*TEMP2)
                QEFF(J,I-1) = 8.63D-06 * CS(I-1,J) / (G(J)*TEMP2)
              else
                QEFF(I-1,J) = CS(I-1,J) * 10. ** IRATS
                                                                     !Be careful
                QEFF(J,I-1) = G(I-1) * QEFF(I-1,J) / (EXPE * G(J))
                                                                     !G integer!
              endif
            enddo
          enddo
          do I = 2, NLEV
            do J = 1, NLEV
              if (J.NE.I) then
                X(I,J) = X(I,J) + DENS * QEFF(J,I)
                X(I,I) = X(I,I) - DENS * QEFF(I,J)
                if (J.GT.I) then
                  X(I,J) = X(I,J) + A(J,I)
                else
                  X(I,I) = X(I,I) - A(I,J)
                endif
              endif
            enddo
          enddo
          do I = 2, NLEV
            IM1 = I - 1
            VALUE = 0.D0 - X(I,1)
            Y(IM1) = VALUE
            Y2(IM1) = VALUE
            YKEEP(IM1) = VALUE
            do J = 2, NLEV
              JM1 = J - 1
              VALUE = X(I,J)
              X(IM1,JM1) = VALUE
              X2(IM1,JM1) = VALUE
              XKEEP(IM1,JM1) = VALUE
            enddo
          enddo
                                              !Solve matrices for populations
          call LUSLV(X,Y,NLEV1,NDIM2)
          do I = NLEV, 2, -1
            N(I) = Y(I-1)
          enddo
          SUMN = 1.D0
          do I = 2, NLEV
            SUMN = SUMN + N(I)
          enddo
          do I = 2, NLEV
            N(I) = N(I) / SUMN
          enddo
          N(1) = 1.D0 / SUMN
                                                                  !Output data
          TTT=TEMP*1.0D-4
          TTP=TTT**(-0.87D0)
                                       !Eff. recombination coef. of Hb
          AHB=3.036D-14*TTP
          do I = 1, NLEV1
            IP1 = I + 1
            do J = IP1, NLEV
               if (A(J,I).NE.0.D0) then
                 EJI = E(J) - E(I)
                 WAV = 1.D8 / EJI
                 RLINT = A(J,I) * EJI
                 RLINT = RLINT *N(J)
                 TNIJ(I,J)=RLINT
                 FINT=N(J)*A(J,I)*4861.D0/(DENS*AHB*WAV)
                 FINTIJ(I,J)=FINT
               endif
            enddo
          enddo
                     !Search ITRANA, ITRANB & ITRANC for transitions & sum up
          SUMA=0.D0
          SUMB=0.D0

          IAPR=0
          IBPR=0
          ICPR=0
          do IKT = 1, NDIM2
            IA1=ITRANA(1,IKT)
            IA2=ITRANA(2,IKT)
            if (IA1.NE.0.AND.IA2.NE.0) then
              SUMA=SUMA+TNIJ(IA1,IA2)
              IAPR=IAPR+1
            endif
            IB1=ITRANB(1,IKT)
            IB2=ITRANB(2,IKT)
            if (IB1.NE.0.AND.IB2.NE.0) then
             IBPR=IBPR+1
             SUMB=SUMB+TNIJ(IB1,IB2)
            endif

            IC1=ITRANC(1,IKT)
            IC2=ITRANC(2,IKT)
            if (IC1.NE.0.AND.IC2.NE.0) then
             ICPR=ICPR+1
            endif
          enddo
          FRAT=SUMA/SUMB

          if (diagtype .eq. "t" .or. diagtype .eq. "T") then
            RESULTS(1, JT) = TEMP
            RESULTS(2, JT) = DENS
            RESULTS(3, JT) = FRAT-inratio
          else
            RESULTS(1, JJD) = TEMP
            RESULTS(2, JJD) = DENS
            RESULTS(3, JJD) = FRAT-inratio
          endif                                   !End of the Ne loop
        enddo

        do IA = 1, IAPR
          I1=ITRANA(1,IA)
          I2=ITRANA(2,IA)
          DEE=E(I2)-E(I1)
          WAVA(IA)=1.D8/DEE
        enddo
        do IB = 1, IBPR
          I1=ITRANB(1,IB)
          I2=ITRANB(2,IB)
          DEE=E(I2)-E(I1)
          WAVB(IB)=1.D8/DEE
        enddo
        do IC = 1, ICPR
          I1=ITRANC(1,IC)
          I2=ITRANC(2,IC)
          DEE=E(I2)-E(I1)
          WAVC(IC)=1.D8/DEE
        enddo

      enddo                                              !End of the Te loop

! here, find the value in RESULTS which is closest to zero
! sort values in results to find two lowest values
!print*,results
!print*," "
      if (diagtype .eq. "D" .or. diagtype .eq. "d") then
        INT = ind
      endif

        ! loop through array and find out where the sign changes.

        do I=2,INT
            test=0
                if (sign(results(3,I),results(3,1)) .ne. results(3,I)) then !when this condition is fulfilled, the values in the array are now a different sign to the first value in the array
                        valtest(:) = (results(:,I-1)) ! return the value before the sign change so that the next loop starts at a sensible value
                        test=1
                        exit
                endif
        enddo

        if(test .eq. 0 .and. loop .lt. 9) then !test fails if no change of sign
                             !this kicks in then, and checks if it should be upper or lower limit
            if(abs(results(3,1)) .lt. abs(results(3,INT))) then
                valtest(:)=results(:,1)
            elseif(abs(results(3,INT)) .lt. abs(results(3,1))) then
                valtest(:)=results(:,INT-1)
            else                   !A simplistic work-around for this problem:
                result = -1d0      !it flags the value as ill-defined so that it
                deallocate(results)!can be dealt with without breaking the code.
                return             !It's set to zero later and excluded from
                                   !the averaging. Turns out long train rides
                                   !are good for sorting such problems.
                !print*,"Valtest failed"
                !print*,ion,levu,levl,loop,inratio,diagtype
                !print*,results
                !STOP
            endif
        elseif(test .eq. 0 .and. loop .eq. 9) then !test fails if no change of sign
                             !this kicks in then, and checks if it should be upper or lower limit
            if(abs(results(3,1)) .lt. abs(results(3,INT))) then
                valtest(:)=results(:,1)
            elseif(abs(results(3,INT)) .lt. abs(results(3,1))) then
                valtest(:)=results(:,INT)
            else
                result = -1d0
                deallocate(results)
                return
                !print*,"Valtest failed"
                !print*,ion,levu,levl,loop,inratio,diagtype
                !print*,results
                !STOP
            endif
        endif

        !LOOP = LOOP + 1
        DEALLOCATE(RESULTS) !thanks Bruce!
      enddo
      !********LOOP WOULD END HERE**********************

      if (diagtype .eq. "D" .or. diagtype .eq. "d") then
        result = valtest(2)
      else
        result = valtest(1)
      endif

      return

 6100 FORMAT (' PROCESSING COMPLETED'/                                  &
     & ' GOODBYE!!'///)


      END subroutine get_diagnostic

      subroutine get_abundance(ion,levels,tempi,densi,iobs,abund,ndim2,ndim1,atomicdata,iion)
      use mod_atomicdata
      implicit none
      integer, parameter :: dp = kind(1.d0)

          !INTEGER maxlevs,maxtemps
      integer :: NDIM1, NDIM2, NDIM1T3, MAXND
                                                      !Maximum no of Te & levels
      !PARAMETER (NDIM1=maxtemps, NDIM2=maxlevs)!35, NDIM2=150)
                                             !NDIM1T3 should be at least 3*NDIM1
      !PARAMETER (NDIM1T3 = 3*NDIM1)!105)
                                                   !Maximum no. of Ne increments
      PARAMETER (MAXND=100)
      integer :: G(NDIM2),                                                 &
     &  ITRANA(2,NDIM2),ITRANB(2,NDIM2),ITRANC(2,NDIM2)
      type(atomic_data),dimension(:),intent(in) :: atomicdata
      integer :: nion,iion
      real(kind=dp) :: N(NDIM2)
      real(kind=dp) :: TDRAT(2,MAXND), TNIJ(NDIM2,NDIM2), FINTIJ(NDIM2,NDIM2),    &
     & WAVA(NDIM2), WAVB(NDIM2), WAVC(NDIM2), CS(NDIM2,NDIM2),          &
     & QEFF(NDIM2,NDIM2), QQ(NDIM1),                                    &
     & QOM(NDIM1,NDIM2,NDIM2), A(NDIM2,NDIM2), E(NDIM2), T(NDIM1),      &
     & ROOTT(NDIM1), X(NDIM2,NDIM2), Y(NDIM2),                          &
     & X2(NDIM2,NDIM2), XKEEP(NDIM2,NDIM2), Y2(NDIM2), YKEEP(NDIM2),    &
     & HMH(NDIM1,NDIM1), D(NDIM1)
      character(len=20) :: LABEL(NDIM2)
      character(len=10) :: ION
      integer :: I, I1, I2, J, KK, LL, JT, JJD,                            &
     & NLEV, NTEMP, IBIG, IRATS,                                        &
     & NLEV1, INT, IND, IOPT, IT, IM1, JM1, IP1,                        &
     & IAPR, IBPR, ICPR, IKT, IA, IB, IC, IA1, IA2, IB1, IB2, IC1, IC2
      real(kind=dp) :: TEMPI, TINC, DENSI, DINC, DENS, DLOGD, TEMP, TLOGT,        &
     & TEMP2, DD, DELTEK, EXPE, VALUE, SUMN, TTT, TTP, AHB, EJI, WAV,   &
     & RLINT, FINT, SUMA, SUMB, SUMC, FRAT, DEE

      character(len=20) :: levels
      real(kind=dp) :: iobs, abund

!debugging
#ifdef CO
        print *,"subroutine: get_abundance. ion=",ion
#endif

      ndim1t3=3*ndim1
      g=0
      itrana=0
      itranb=0
      itranc=0

      read(levels,*) ((ITRANC(LL,KK),LL=1,2),KK=1,ndim2)!150)

      tinc=0
      dinc=0
      int=1
      ind=1

      nion = 0
      do i = 1,iion
         if(atomicdata(i)%ion .eq. ion(1:10)) nion=i
      enddo
      if (nion .eq. 0) then
          print *, "Dave, my mind is going. I can feel it."
          nion = 1
      endif
!          print*,nion,atomicdata(nion)%ion,ion
      nlev=atomicdata(nion)%nlevs
      ntemp=atomicdata(nion)%ntemps
          T(1:ntemp)=log10(atomicdata(nion)%Temps(1:ntemp))
      ROOTT(1:ntemp)=atomicdata(nion)%rootT(1:ntemp)
      Label(1:nlev)=atomicdata(nion)%labels(1:nlev)
      QOM(1:ntemp,1:nlev,1:nlev)=atomicdata(nion)%col_str(1:ntemp,1:nlev,1:nlev)
      A(1:nlev,1:nlev)=atomicdata(nion)%A_coeffs(1:nlev,1:nlev)
      E(1:nlev)=atomicdata(nion)%Waveno(1:nlev)
      G(1:nlev)=atomicdata(nion)%G(1:nlev)

      irats=0
      ibig=0

      NLEV1 = NLEV - 1
                                                               !Start of Te loop
      do JT = 1, INT
        TEMP=TEMPI+(JT-1)*TINC

                                                               !Start of Ne loop
        do JJD = 1, IND
          DENS=DENSI+(JJD-1)*DINC
          if (TEMP.LE.0.D0.OR.DENS.LE.0.D0) then
            WRITE (6,6100)
            STOP
          endif
          DLOGD = LOG10 (DENS)
          TLOGT = LOG10 (TEMP)
          TEMP2= SQRT (TEMP)
                                                                  !Form matrices
          X = 0.D0
          CS=0.D0
          QEFF=0.D0
          TNIJ = 0.D0
          Y=0.D0

          IOPT=0
          if (NTEMP.EQ.1) then
            WRITE (6,*)
            WRITE (6,*)                                                 &
     &      'Coll. strengths available for 1 Te only - assuming const'
          elseif (NTEMP.EQ.2) then
            WRITE (6,*)
            WRITE (6,*)                                                 &
     &      'Coll. strengths available for 2 Te only - linear interp'
          else
            call SPLMAT(T, NTEMP, IOPT, NDIM1, NDIM1T3, HMH)
            call CFD(TLOGT,T,NTEMP,NDIM1,HMH,D)
          endif
          do I = 2, NLEV
            do J = I, NLEV
                                                           !Negative!
              DELTEK = (E(I-1)-E(J))*1.4388463D0
              EXPE = EXP(DELTEK/TEMP)
              do IT = 1, NTEMP
                if (IRATS.EQ.0.D+00) then
                  QQ(IT) = QOM(IT,I-1,J)
                else
                                                      !Take out the exp. depend.
                  QQ(IT) = QOM(IT,I-1,J) / EXPE
                                                      !before interpolation
                endif
              enddo
              if (NTEMP.EQ.1) then
                 DD = QQ(1)
              elseif (NTEMP.EQ.2) then
                 DD = QQ(1) +                                           &
     &            (QQ(2) - QQ(1))/(T(2) - T(1)) * (TLOGT - T(1))
              else
                call CFY(TLOGT, DD, T, QQ, NTEMP, NDIM1, D)
              endif
              if (IRATS.EQ.0.D+00) then
                CS(I-1,J) = DD
              else
                CS(I-1,J) = DD * EXPE
              endif
              if (IRATS .EQ. 0.D+00) then
                QEFF(I-1,J) = 8.63D-06*CS(I-1,J) * EXPE / (G(I-1)*TEMP2)
                QEFF(J,I-1) = 8.63D-06 * CS(I-1,J) / (G(J)*TEMP2)
              else
                QEFF(I-1,J) = CS(I-1,J) * 10. ** IRATS
                                                                     !Be careful
                QEFF(J,I-1) = G(I-1) * QEFF(I-1,J) / (EXPE * G(J))
                                                                     !G integer!
              endif
            enddo
          enddo
          do I = 2, NLEV
            do J = 1, NLEV
              if (J.NE.I) then
                X(I,J) = X(I,J) + DENS * QEFF(J,I)
                X(I,I) = X(I,I) - DENS * QEFF(I,J)
                if (J.GT.I) then
                  X(I,J) = X(I,J) + A(J,I)
                else
                  X(I,I) = X(I,I) - A(I,J)
                endif
              endif
            enddo
          enddo
          do I = 2, NLEV
            IM1 = I - 1
            VALUE = 0.D0 - X(I,1)
            Y(IM1) = VALUE
            Y2(IM1) = VALUE
            YKEEP(IM1) = VALUE
            do J = 2, NLEV
              JM1 = J - 1
              VALUE = X(I,J)
              X(IM1,JM1) = VALUE
              X2(IM1,JM1) = VALUE
              XKEEP(IM1,JM1) = VALUE
            enddo
          enddo
                                              !Solve matrices for populations
          call LUSLV(X,Y,NLEV1,NDIM2)
          do I = NLEV, 2, -1
            N(I) = Y(I-1)
          enddo
          SUMN = 1.D0
          do I = 2, NLEV
            SUMN = SUMN + N(I)
          enddo
          do I = 2, NLEV
            N(I) = N(I) / SUMN
          enddo
          N(1) = 1.D0 / SUMN
                                                                  !Output data
          TTT=TEMP*1.0D-4
          TTP=TTT**(-0.87D0)
                                       !Eff. recombination coef. of Hb
          AHB=3.036D-14*TTP
          do I = 1, NLEV1
            IP1 = I + 1
            do J = IP1, NLEV
               if (A(J,I).NE.0.D0) then
                 EJI = E(J) - E(I)
                 WAV = 1.D8 / EJI
                 RLINT = A(J,I) * EJI
                 RLINT = RLINT *N(J)
                 TNIJ(I,J)=RLINT
                 FINT=N(J)*A(J,I)*4861.D0/(DENS*AHB*WAV)
                 FINTIJ(I,J)=FINT
               endif
            enddo
          enddo
                     !Search ITRANA, ITRANB & ITRANC for transitions & sum up
          SUMA=0.D0
          SUMB=0.D0
          SUMC=0.D0
          IAPR=0
          IBPR=0
          ICPR=0
          do IKT = 1, NDIM2
            IA1=ITRANA(1,IKT)
            IA2=ITRANA(2,IKT)
            if (IA1.NE.0.AND.IA2.NE.0) then
              SUMA=SUMA+TNIJ(IA1,IA2)
              IAPR=IAPR+1
            endif
            IB1=ITRANB(1,IKT)
            IB2=ITRANB(2,IKT)
            if (IB1.NE.0.AND.IB2.NE.0) then
             IBPR=IBPR+1
             SUMB=SUMB+TNIJ(IB1,IB2)
            endif

            IC1=ITRANC(1,IKT)
            IC2=ITRANC(2,IKT)
            if (IC1.NE.0.AND.IC2.NE.0) then
             ICPR=ICPR+1
             SUMC=SUMC+FINTIJ(IC1,IC2)
            endif
          enddo
          FRAT=SUMA/SUMB
          SUMC = 1./SUMC
          TDRAT(1,JJD)=DENS
          TDRAT(2,JJD)=FRAT
          abund = sumc*iobs/100
                                                       !End of the Ne loop
        enddo
        do IA = 1, IAPR
          I1=ITRANA(1,IA)
          I2=ITRANA(2,IA)
          DEE=E(I2)-E(I1)
          WAVA(IA)=1.D8/DEE
        enddo
        do IB = 1, IBPR
          I1=ITRANB(1,IB)
          I2=ITRANB(2,IB)
          DEE=E(I2)-E(I1)
          WAVB(IB)=1.D8/DEE
        enddo
        do IC = 1, ICPR
          I1=ITRANC(1,IC)
          I2=ITRANC(2,IC)
          DEE=E(I2)-E(I1)
          WAVC(IC)=1.D8/DEE
        enddo
                                                         !End of the Te loop
      enddo

      return

 6100 FORMAT (' PROCESSING COMPLETED'/                                  &
     & ' GOODBYE!!'///)

      END subroutine get_abundance

!---- PROC LUSLV
                                                       !Solving linear equations
      subroutine LUSLV(A,B,N,M)
      implicit none
      integer, parameter :: dp = kind(1.d0)

      integer :: M, N
      real(kind=dp) :: A(M,M),B(M)
#ifdef CO
        print *,"subroutine: luslv"
#endif
      call LURED(A,N,M)
      call RESLV(A,B,N,M)
      return
      END subroutine luslv
!
!---- PROC LURED
      subroutine LURED(A,N,NR)
      implicit none
      integer, parameter :: dp = kind(1.d0)

      integer :: N, NR, NM1, I, J, K, IP1
      real(kind=dp) :: A(NR,NR), FACT

#ifdef CO
        print *,"subroutine: lured"
#endif

      if (N.EQ.1) return
      NM1=N-1
      do I=1,NM1
        IP1=I+1
        do K=IP1,N
          FACT=A(K,I)/A(I,I)
          do J=IP1,N
            A(K,J)=A(K,J)-A(I,J)*FACT
          enddo
        enddo
      enddo
      return
      END subroutine lured
!
!---- PROC RESLV
                                                               !Resolve A with B
      subroutine RESLV(A,B,N,NR)
      implicit none
      integer, parameter :: dp = kind(1.d0)

      integer :: N, NR, NM1, I, J, K, L, IP1
      real(kind=dp) :: A(NR,NR),B(NR)

#ifdef CO
        print *,"subroutine: reslv"
#endif

      if (N.EQ.1) GOTO 1
      NM1=N-1
      do I=1,NM1
        IP1=I+1
        do J=IP1,N
          B(J)=B(J)-B(I)*A(J,I)/A(I,I)
        enddo
      enddo
      B(N)=B(N)/A(N,N)
      do I=1,NM1
        K=N-I
        L=K+1
        do J=L,N
          B(K)=B(K)-B(J)*A(K,J)
        enddo
        B(K)=B(K)/A(K,K)
      enddo
      return
    1 B(N)=B(N)/A(N,N)
      return
      END subroutine reslv
!
!---- PROC SPLMAT
      subroutine SPLMAT(XX,NPT,IOPT,NDIM, NDIMT3, HMH)
      implicit none
      integer, parameter :: dp = kind(1.d0)

      integer :: NDIM, NDIMT3, NPT, IOPT, NPM, NELEM
      real(kind=dp) :: XX(NDIM),GH(NDIMT3),Y(NDIM), HMH(NDIM,NDIM)

#ifdef CO
        print *,"subroutine: splmat"
#endif

      NPM=NPT-2
      call GHGEN(GH,XX,NPT,IOPT,NDIM,NDIMT3)
      NELEM=3*NPM-2
      call ELU(GH,NPM,NDIMT3)
      call HGEN(XX,GH,Y,NPT,IOPT,NDIM,NDIMT3,HMH)
      return
      END subroutine splmat
!
!---- PROC DERIV
!     Calculate the first derivative of the lagrangian interpolator
!     of a function F, tabulated at the N points XY(I), I=1 to N.
!     The derivative is given as the coefficients of F(I), I=1 to N,
!     in the array D(I), I=1 to N.
      subroutine DERIV(XY,D,X,N,NDIM)
      implicit none
      integer, parameter :: dp = kind(1.d0)

      integer :: N ,NDIM, I, J, K
      real(kind=dp) :: XY(NDIM),D(NDIM), X, P1, P2, S

#ifdef CO
        print *,"subroutine: deriv"
#endif

      do I=1,N
        P1=1.
        S=0.
        do J=1,N
          if (J.NE.I) then
            P1=P1*(XY(I)-XY(J))
            P2=1.
            do K=1,N
              if (K.NE.I.AND.K.NE.J) P2=P2*(X-XY(K))
            enddo
            S=S+P2
          endif
        enddo
        D(I)=S/P1
      enddo
      return
      END subroutine deriv
!
!---- PROC HGEN
!     Cubic spline interpolation
!     The equation for the second derivatives at internal points
!     is of the form G*YPP=B, where G has been evaluated and LU
!     decomposed.
!     this routine writes B=HMH*Y and then solves YPP=G**(-1)*HMH*Y,
!     =HMH*Y.
!     Three options are provided for boundary conditions-
!     IOPT = 0  YPP=0 at end points
!     IOPT = 1  YP=0  at end points
!     IOPT = 2  YP at end points from lagarnge interpolant of a set of
!     internal points.
      subroutine HGEN(XX,GH,Y,NPT,IOPT,NDIM,NDIMT3,HMH)
      implicit none
      integer, parameter :: dp = kind(1.d0)

      integer :: NPT, IOPT, NDIM, NDIMT3, NDIM3, NIP, I, J, K, NPM,        &
     & INDX
      real(kind=dp) :: XX(NDIM), GH(NDIMT3), Y(NDIM), HMH(NDIM,NDIM),             &
     & XY(5),D(5),C(2,5), A0, AN1, H1, H2

#ifdef CO
        print *,"subroutine: hgen"
#endif

      if (IOPT.EQ.2) then          !Case of derivative boundary condition, with
        NDIM3=5                    !derivatives from NIP-point Lagrange at
        NIP=3                      !internal points
        do J=1,2
          do I=1,NIP
            K=(NPT-NIP)*(J-1)
            XY(I)=XX(K+I)
          enddo
          K=1+(NPT-1)*(J-1)
          call DERIV(XY,D,XX(K),NIP,NDIM3)
          do I=1,NIP
            C(J,I)=D(I)
          enddo
        enddo
      endif
                                             !Set up matrix equation G*YPP=HMH*Y
      A0=XX(2)-XX(1)
      AN1=XX(NPT)-XX(NPT-1)
      NPM=NPT-2
      hmh=0d0
      do I=1,NPM
        H1=6./(XX(I+1)-XX(I))
        H2=6./(XX(I+2)-XX(I+1))
        HMH(I,I)=H1
        HMH(I,I+1)=-H1-H2
        HMH(I,I+2)=H2
!        do J=1,NPT
!          HMH(I,J)=0.
!          if (J.EQ.I) HMH(I,J)=H1
!          if (J.EQ.I+2) HMH(I,J)=H2
!          if (J.EQ.I+1) HMH(I,J)=-H1-H2
!        enddo
      enddo
                                                 !Correct matrix for case of
      if (IOPT.EQ.1.OR.IOPT.EQ.2) then
                                                 !derivative boundary conditions
        HMH(1,1)=HMH(1,1)+3/A0
        HMH(1,2)=HMH(1,2)-3/A0
        HMH(NPM,NPT-1)=HMH(NPM,NPT-1)-3/AN1
        HMH(NPM,NPT)=HMH(NPM,NPT)+3/AN1
      endif
      if (IOPT.EQ.2) then
        do J=1,NIP
          HMH(1,J)=HMH(1,J)+3*C(1,J)
          K=NPT+J-NIP
          HMH(NPM,K)=HMH(NPM,K)-3*C(2,J)
        enddo
      endif
                                 !Solve matrix equation with results in the form
      do I=1,NPT
                                 !YPP=HMH*Y. matrix g has been LU decomposed
        Y(1)=HMH(1,I)
        INDX=0
        do J=2,NPM
          INDX=INDX+3
          Y(J)=HMH(J,I)-GH(INDX)*Y(J-1)
        enddo
        INDX=INDX+1
        Y(NPM)=Y(NPM)/GH(INDX)
        do J=2,NPM
          K=NPM-J+1
          INDX=INDX-3
          Y(K)=(Y(K)-GH(INDX+1)*Y(K+1))/GH(INDX)
        enddo

        HMH(2:npm+1,I)=Y(1:npm)
                                    !Insert values for second derivative at end
        HMH(1,I)=0.
                                    !points: first and last rows of the matrix
        HMH(NPT,I)=0.
      enddo
                                         !Case of derivative boundary conditions
      if (IOPT.GT.0) then
        do J=1,NPT
          HMH(1,J)=-0.5*HMH(2,J)
          HMH(NPT,J)=-0.5*HMH(NPT-1,J)
        enddo
        HMH(1,1)=HMH(1,1)-3/(A0*A0)
        HMH(1,2)=HMH(1,2)+3/(A0*A0)
        HMH(NPT,NPT-1)=HMH(NPT,NPT-1)+3/(AN1*AN1)
        HMH(NPT,NPT)=HMH(NPT,NPT)-3/(AN1*AN1)
      endif
      if (IOPT.EQ.2) then
        do J=1,NIP
          HMH(1,J)=HMH(1,J)-3*C(1,J)/A0
          K=NPT+J-NIP
          HMH(NPT,K)=HMH(NPT,K)+3*C(2,J)/AN1
        enddo
      endif
      return
      END subroutine hgen
!
!---- PROC GHGEN
      subroutine GHGEN(GH,XX,NPT,IOPT,NDIM,NDIMT3)
      implicit none
      integer, parameter :: dp = kind(1.d0)

      integer :: NPT, IOPT, NDIM, NDIMT3, INDX, NPTM, I, J, IP, JP, IK
      real(kind=dp) :: XX(NDIM),GH(NDIMT3)

#ifdef CO
        print *,"subroutine: ghgen"
#endif

      INDX=0
      NPTM=NPT-1
      do I=2,NPTM
        IP=I-1
        do J=1,3
          JP=IP+J-2
          if (JP.GE.1.AND.JP.LE.NPTM-1) then
            INDX=INDX+1
            if (J.EQ.2) then
              GH(INDX)=2*(XX(I+1)-XX(I-1))
            else
              IK=I+(J-1)/2
              GH(INDX)=XX(IK)-XX(IK-1)
            endif
          endif
        enddo
      enddo
      if (IOPT.GE.1) then
        GH(1)=GH(1)-(XX(2)-XX(1))/2.
        GH(INDX)=GH(INDX)-(XX(NPT)-XX(NPT-1))/2.
      endif
      return
      END subroutine ghgen
!
!---- PROC ELU
      subroutine ELU(GH,N,NDIM)
      implicit none
      integer, parameter :: dp = kind(1.d0)

      integer :: N, NDIM, INDX, I, J, JP
      real(kind=dp) :: GH(NDIM)

#ifdef CO
        print *,"subroutine: elu"
#endif

      INDX=0
      do I=1,N
        do J=1,3
          JP=I+J-2
          if (JP.GE.1.AND.JP.LE.N) then
            INDX=INDX+1
            if (I.GT.1) then
              if (J.EQ.1) then
                GH(INDX)=GH(INDX)/GH(INDX-2)
              endif
              if (J.EQ.2) then
                GH(INDX)=GH(INDX)-GH(INDX-1)*GH(INDX-2)
              endif
            endif
          endif
        enddo
      enddo
      return
      END subroutine elu
!
!---- PROC CFY
      subroutine CFY(X,Y,XX,YY,NPT,NDIM,D)
      implicit none
      integer, parameter :: dp = kind(1.d0)

      integer :: NPT, NDIM, J
      real(kind=dp) :: XX(NDIM),YY(NDIM), D(NDIM), X, Y, TT

#ifdef CO
        print *,"subroutine: cfy"
#endif

      if (X.LT.XX(1)) then
        Y=YY(1)
      endif
      if (X.GT.XX(NPT)) then
        Y=YY(NPT)
      endif
      TT=0.
      do J=1,NPT
        TT=TT+D(J)*YY(J)
      enddo
      Y=TT
      return

      END subroutine cfy
!
!---- PROC CFD
      subroutine CFD(X,XX,NPT,NDIM, HMH, D)
      implicit none
      integer, parameter :: dp = kind(1.d0)

      integer :: NPT, NDIM, NPTM, I, J
      real(kind=dp) :: X, XX(NDIM), HMH(NDIM,NDIM), D(NDIM), X1, X2, A1, A2, HI

#ifdef CO
        print *,"subroutine: cfd"
#endif

      if (X.LT.XX(1)) then
        !WRITE(6,400) XX(1)
        return
      endif
      if (X.GT.XX(NPT)) then
        !WRITE(6,401) XX(NPT)
        return
      endif
      NPTM=NPT-1
      do I=1,NPTM
        if (X.LT.XX(I+1)) then
          X1=XX(I+1)-X
          X2=X-XX(I)
          HI=XX(I+1)-XX(I)
          A1=X1*(X1*X1/(6*HI)-HI/6)
          A2=X2*(X2*X2/(6*HI)-HI/6)
          do J=1,NPT
            D(J)=(A1*HMH(I,J)+A2*HMH(I+1,J))
          enddo
          D(I)=D(I)+X1/HI
          D(I+1)=D(I+1)+X2/HI
          return
        endif
      enddo

      END subroutine cfd

      end module mod_equib
