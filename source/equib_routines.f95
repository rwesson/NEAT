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
      IMPLICIT NONE

      INTEGER NDIM1, NDIM2, NDIM1T3, MAXND
                                                      !Maximum no of Te & levels
      !PARAMETER (NDIM1=maxtemps, NDIM2=maxlevs)!35, NDIM2=150)
                                             !NDIM1T3 should be at least 3*NDIM1
      !PARAMETER (NDIM1T3 = 3*NDIM1)!105)
                                                   !Maximum no. of Ne increments
      PARAMETER (MAXND=100)
      INTEGER G(NDIM2),                                                 &
     &  ITRANA(2,NDIM2),ITRANB(2,NDIM2),ITRANC(2,NDIM2),LOOP
      type(atomic_data),dimension(:),intent(in) :: atomicdata
      integer iion,nion
      REAL*8 N(NDIM2)
      REAL*8 & !TDRAT(2,MAXND)
     & TNIJ(NDIM2,NDIM2), FINTIJ(NDIM2,NDIM2),                          &
     & WAVA(NDIM2), WAVB(NDIM2), WAVC(NDIM2), CS(NDIM2,NDIM2),          &
     & QEFF(NDIM2,NDIM2), QQ(NDIM1),                                    &
     & QOM(NDIM1,NDIM2,NDIM2), A(NDIM2,NDIM2), E(NDIM2), T(NDIM1),      &
     & ROOTT(NDIM1), X(NDIM2,NDIM2), Y(NDIM2),                          &
     & X2(NDIM2,NDIM2), XKEEP(NDIM2,NDIM2), Y2(NDIM2), YKEEP(NDIM2),    &
     & HMH(NDIM1,NDIM1), D(NDIM1)
      CHARACTER(len=20) LABEL(NDIM2)
      CHARACTER(len=10) ION 
      INTEGER I, I1, I2, J, KK, LL, JT, JJD,                            &
     & NLEV, NTEMP, IBIG, IRATS,                                        &
     & NLEV1, INT, IND, IOPT, IT, IM1, JM1, IP1,                        &
     & IAPR, IBPR, ICPR, IKT, IA, IB, IC, IA1, IA2, IB1, IB2, IC1, IC2
      REAL*8 TEMPI, TINC, DENSI, DINC, DENS, DLOGD, TEMP, TLOGT,        &
     & TEMP2, DD, DELTEK, EXPE, VALUE, SUMN, TTT, TTP, AHB, EJI, WAV,   &
     & RLINT, FINT, SUMA, SUMB, SUMC, FRAT, DEE

      DOUBLE PRECISION :: fixedq
      REAL*8 inratio,result
      CHARACTER(len=20) levu,levl
      CHARACTER(len=1) diagtype
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: RESULTS
      REAL*8 valtest(3)
      integer test

          ndim1t3=3*ndim1
      g=0
      itrana=0
      itranb=0
      itranc=0
      valtest=0
      test=0

      READ(levu,*) ((ITRANA(LL,KK),LL=1,2),KK=1,ndim2)!150)
      READ(levl,*) ((ITRANB(LL,KK),LL=1,2),KK=1,ndim2)!150)


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



!      IONL = INDEX(ION,' ') - 1
!      OPEN(UNIT=1,STATUS='OLD',                                         &
!     & FILE='Atomic-data/'//ION(1:IONL)//'.dat')
                                                      !Read in no. comment lines
!      READ(1,*) NLINES
!      DO I = 1, NLINES
                                                                       !Comments
!        READ(1,1003) LTEXT
!      ENDDO
                                          !Read no. of levels (max=NDIM2) NLEV,
!      READ (1,*) NLEV, NTEMP
                                          !no. of Te (max=NDIM1) NTEMP and the
!      DO I = 1, NLEV
                                          !input format (cf Readme)
!         READ (1,1002) LABEL(I)
!      ENDDO
!     be
      ibig=0
                                            !Read in Te's where coll. strengths are tabulated
!      DO I = 1, NTEMP
!           READ (1,*) T(I)
!           T(I) = LOG10 (T(I))
!           ROOTT(I) = SQRT(T(I))
!      ENDDO



                             !If IRATS=0, what tabulated are collision strengths
!      READ(1,*) IRATS
!                            !Else Coll. rates = tabulated values * 10 ** IRATS
!      IF(IBIG.EQ.0) THEN
!   10   READ (1,*) ID(2), JD(2), QX
!        IF (QX.EQ.0.D0) GOTO 20
!        IF (ID(2).EQ.0) THEN
!          ID(2) = ID(1)
!          K = K + 1
!        ELSE
!          ID(1) = ID(2)
!          K = 1
!        ENDIF
!        IF (JD(2).EQ.0) THEN
!          JD(2) = JD(1)
!        ELSE
!          JD(1) = JD(2)
!        ENDIF
!        I = ID(2)
!        J = JD(2)
!        QOM(K,I,J) = QX
!        GO TO 10
!      ENDIF
!   20 IF(IBIG.EQ.1.OR.IBIG.EQ.2) THEN
!        READ(1,*) NTRA
!        DO IN = 1, NTRA
!          READ(1,*) I,J,(QOM(ITEMP,I,J),ITEMP=1,NTEMP)
!        ENDDO
!      ENDIF
                                                  !Read transition probabilities
      NLEV1 = NLEV - 1
!      IF (IBIG.EQ.1) THEN
!       READ(1,7000) ((I,J,A(J,I),L=K+1,NLEV),K=1,NLEV1)
!      ELSE
!      DO K = 1, NLEV1
!        KP1 = K + 1
!          DO L = KP1, NLEV
!            READ (1,*) I, J, AX
!            A(J,I) = AX
!          ENDDO
!      ENDDO
!      ENDIF
                                 !Read statistical weights, energy levels (cm-1)
!      DO J = 1, NLEV
!        READ (1,*) I, GX, EX
!        G(I) = GX
!        E(I) = EX
!      ENDDO
!      CLOSE (UNIT=1)
                                !Get levels for ratio
                                                           !150 large enough
      ITRANC = 0

!newbit

!print*,ion,diagtype,fixedq,int
! set up T and D loops depending on input.
                                               !Read in Te and Ne where the line
                                               !ratio is to be calculated

      !*****LOOP STARTS HERE*************************
      DO LOOP = 1, 9
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

        ALLOCATE(RESULTS(3,INT))

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
      DO JT = 1, INT
        TEMP=TEMPI+(JT-1)*TINC
                                                               !Start of Ne loop

        DO JJD = 1, IND
          DENS=DENSI+(JJD-1)*DINC
!          IF(DENSI.LT.30.D0) THEN
!            DENS=10.D0**DENS
!          ENDIF
          IF (TEMP.LE.0.D0.OR.DENS.LE.0.D0) THEN
            WRITE (6,6100)
                print *,"Temp = ", TEMP, ", Dens = ", DENS, ", Ion = ",ion,diagtype
            STOP
          ENDIF
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
          IF (NTEMP.EQ.1) THEN
            WRITE (6,*)
            WRITE (6,*)                                                 &
     &      'Coll. strengths available for 1 Te only - assuming const'
          ELSEIF (NTEMP.EQ.2) THEN
            WRITE (6,*)
            WRITE (6,*)                                                 &
     &      'Coll. strengths available for 2 Te only - linear interp'
          ELSE
            CALL SPLMAT(T, NTEMP, IOPT, NDIM1, NDIM1T3, HMH)
            CALL CFD(TLOGT,T,NTEMP,NDIM1,HMH,D)
          ENDIF
          DO I = 2, NLEV
            DO J = I, NLEV
                                                           !Negative!
              DELTEK = (E(I-1)-E(J))*1.4388463D0
              EXPE = EXP(DELTEK/TEMP)
              DO IT = 1, NTEMP
                IF (IRATS.EQ.0.D+00) THEN
                  QQ(IT) = QOM(IT,I-1,J)
                ELSE
                                                      !Take out the exp. depend.
                  QQ(IT) = QOM(IT,I-1,J) / EXPE
                                                      !before interpolation
                ENDIF
              ENDDO
              IF (NTEMP.EQ.1) THEN
                 DD = QQ(1)
              ELSEIF (NTEMP.EQ.2) THEN
                 DD = QQ(1) +                                           &
     &            (QQ(2) - QQ(1))/(T(2) - T(1)) * (TLOGT - T(1))
              ELSE
                CALL CFY(TLOGT, DD, T, QQ, NTEMP, NDIM1,HMH,D)
              ENDIF
              IF (IRATS.EQ.0.D+00) THEN
                CS(I-1,J) = DD
              ELSE
                CS(I-1,J) = DD * EXPE
              ENDIF
              IF (IRATS .EQ. 0.D+00) THEN
                QEFF(I-1,J) = 8.63D-06*CS(I-1,J) * EXPE / (G(I-1)*TEMP2)
                QEFF(J,I-1) = 8.63D-06 * CS(I-1,J) / (G(J)*TEMP2)
              ELSE
                QEFF(I-1,J) = CS(I-1,J) * 10. ** IRATS
                                                                     !Be careful
                QEFF(J,I-1) = G(I-1) * QEFF(I-1,J) / (EXPE * G(J))
                                                                     !G integer!
              ENDIF
            ENDDO
          ENDDO
          DO I = 2, NLEV
            DO J = 1, NLEV
              IF (J.NE.I) THEN
                X(I,J) = X(I,J) + DENS * QEFF(J,I)
                X(I,I) = X(I,I) - DENS * QEFF(I,J)
                IF (J.GT.I) THEN
                  X(I,J) = X(I,J) + A(J,I)
                ELSE
                  X(I,I) = X(I,I) - A(I,J)
                ENDIF
              ENDIF
            ENDDO
          ENDDO
          DO I = 2, NLEV
            IM1 = I - 1
            VALUE = 0.D0 - X(I,1)
            Y(IM1) = VALUE
            Y2(IM1) = VALUE
            YKEEP(IM1) = VALUE
            DO J = 2, NLEV
              JM1 = J - 1
              VALUE = X(I,J)
              X(IM1,JM1) = VALUE
              X2(IM1,JM1) = VALUE
              XKEEP(IM1,JM1) = VALUE
            ENDDO
          ENDDO
                                              !Solve matrices for populations
          CALL LUSLV(X,Y,NLEV1,NDIM2)
          DO I = NLEV, 2, -1
            N(I) = Y(I-1)
          ENDDO
          SUMN = 1.D0
          DO I = 2, NLEV
            SUMN = SUMN + N(I)
          ENDDO
          DO I = 2, NLEV
            N(I) = N(I) / SUMN
          ENDDO
          N(1) = 1.D0 / SUMN
                                                                  !Output data
          TTT=TEMP*1.0D-4
          TTP=TTT**(-0.87D0)
                                       !Eff. recombination coef. of Hb
          AHB=3.036D-14*TTP
          DO I = 1, NLEV1
            IP1 = I + 1
            DO J = IP1, NLEV
               IF (A(J,I).NE.0.D0) THEN
                 EJI = E(J) - E(I)
                 WAV = 1.D8 / EJI
                 RLINT = A(J,I) * EJI
                 RLINT = RLINT *N(J)
                 TNIJ(I,J)=RLINT
                 FINT=N(J)*A(J,I)*4861.D0/(DENS*AHB*WAV)
                 FINTIJ(I,J)=FINT
               ENDIF
            ENDDO
          ENDDO
                     !Search ITRANA, ITRANB & ITRANC for transitions & sum up
          SUMA=0.D0
          SUMB=0.D0
          SUMC=0.D0
          IAPR=0
          IBPR=0
          ICPR=0
          DO IKT = 1, NDIM2
            IA1=ITRANA(1,IKT)
            IA2=ITRANA(2,IKT)
            IF(IA1.NE.0.AND.IA2.NE.0) THEN
              SUMA=SUMA+TNIJ(IA1,IA2)
              IAPR=IAPR+1
            ENDIF
            IB1=ITRANB(1,IKT)
            IB2=ITRANB(2,IKT)
            IF(IB1.NE.0.AND.IB2.NE.0) THEN
             IBPR=IBPR+1
             SUMB=SUMB+TNIJ(IB1,IB2)
            ENDIf

            IC1=ITRANC(1,IKT)
            IC2=ITRANC(2,IKT)
            IF(IC1.NE.0.AND.IC2.NE.0) THEN
             ICPR=ICPR+1
             SUMC=SUMC+FINTIJ(IC1,IC2)
            ENDIf
          ENDDO
          FRAT=SUMA/SUMB
          SUMC = 1./SUMC
!          TDRAT(1,JJD)=DENS  !are these lines necessary,
!          TDRAT(2,JJD)=FRAT  !TDRAT is now never used again?
!          write(6,*),jd,suma,sumb,sumc,dens,frat
!          WRITE(7,1017) TEMP, DENS, SUMC
!          WRITE(8,1017) TEMP, DENS, FRAT, FRAT-inratio
          if (diagtype .eq. "t" .or. diagtype .eq. "T") then
            RESULTS(1, JT) = TEMP
            RESULTS(2, JT) = DENS
            RESultS(3, JT) = FRAT-inratio
          else
            RESULTS(1, JJD) = TEMP
            RESULTS(2, JJD) = DENS
            RESultS(3, JJD) = FRAT-inratio
          endif                                   !End of the Ne loop
        ENDDO

        DO IA = 1, IAPR
          I1=ITRANA(1,IA)
          I2=ITRANA(2,IA)
          DEE=E(I2)-E(I1)
          WAVA(IA)=1.D8/DEE
        ENDDO
        DO IB = 1, IBPR
          I1=ITRANB(1,IB)
          I2=ITRANB(2,IB)
          DEE=E(I2)-E(I1)
          WAVB(IB)=1.D8/DEE
        ENDDO
        DO IC = 1, ICPR
          I1=ITRANC(1,IC)
          I2=ITRANC(2,IC)
          DEE=E(I2)-E(I1)
          WAVC(IC)=1.D8/DEE
        ENDDO

                                                         !End of the Te loop
      ENDDO
! here, find the value in RESULTS which is closest to zero
! sort values in results to find two lowest values
!print*,results
!print*," "
      if (diagtype .eq. "D" .or. diagtype .eq. "d") then
        INT = ind
      endif

        ! loop through array and find out where the sign changes.

        DO I=2,INT
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
      END DO
      !********LOOP WOULD END HERE**********************

      if (diagtype .eq. "D" .or. diagtype .eq. "d") then
        result = valtest(2)
      else
        result = valtest(1)
      endif

      RETURN

 6100 FORMAT (' PROCESSING COMPLETED'/                                  &
     & ' GOODBYE!!'///)


      END subroutine get_diagnostic

      subroutine get_abundance(ion,levels,tempi,densi,iobs,abund,ndim2,ndim1,atomicdata,iion)
      use mod_atomicdata
      IMPLICIT NONE

          !INTEGER maxlevs,maxtemps
      INTEGER NDIM1, NDIM2, NDIM1T3, MAXND
                                                      !Maximum no of Te & levels
      !PARAMETER (NDIM1=maxtemps, NDIM2=maxlevs)!35, NDIM2=150)
                                             !NDIM1T3 should be at least 3*NDIM1
      !PARAMETER (NDIM1T3 = 3*NDIM1)!105)
                                                   !Maximum no. of Ne increments
      PARAMETER (MAXND=100)
      INTEGER G(NDIM2),                                                 &
     &  ITRANA(2,NDIM2),ITRANB(2,NDIM2),ITRANC(2,NDIM2)
      type(atomic_data),dimension(:),intent(in) :: atomicdata
      integer :: nion,iion
      REAL*8 N(NDIM2)
      REAL*8 TDRAT(2,MAXND), TNIJ(NDIM2,NDIM2), FINTIJ(NDIM2,NDIM2),    &
     & WAVA(NDIM2), WAVB(NDIM2), WAVC(NDIM2), CS(NDIM2,NDIM2),          &
     & QEFF(NDIM2,NDIM2), QQ(NDIM1),                                    &
     & QOM(NDIM1,NDIM2,NDIM2), A(NDIM2,NDIM2), E(NDIM2), T(NDIM1),      &
     & ROOTT(NDIM1), X(NDIM2,NDIM2), Y(NDIM2),                          &
     & X2(NDIM2,NDIM2), XKEEP(NDIM2,NDIM2), Y2(NDIM2), YKEEP(NDIM2),    &
     & HMH(NDIM1,NDIM1), D(NDIM1)
      CHARACTER(len=20) LABEL(NDIM2), ION 
      INTEGER I, I1, I2, J, KK, LL, JT, JJD,                            &
     & NLEV, NTEMP, IBIG, IRATS,                                        &
     & NLEV1, INT, IND, IOPT, IT, IM1, JM1, IP1,                        &
     & IAPR, IBPR, ICPR, IKT, IA, IB, IC, IA1, IA2, IB1, IB2, IC1, IC2
      REAL*8 TEMPI, TINC, DENSI, DINC, DENS, DLOGD, TEMP, TLOGT,        &
     & TEMP2, DD, DELTEK, EXPE, VALUE, SUMN, TTT, TTP, AHB, EJI, WAV,   &
     & RLINT, FINT, SUMA, SUMB, SUMC, FRAT, DEE

      CHARACTER(len=20) levels
      real*8 iobs, abund
!
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

                                                        !Interrogate for input
!      IONL = INDEX(ION,' ') - 1
!      OPEN(UNIT=1,STATUS='OLD',                                         &
!     & FILE='Atomic-data/'//ION(1:IONL)//'.dat')
!!     + NAME='atomic_data/'//ION(1:IONL)//'.dat')
                                                      !Read in no. comment lines
!      READ(1,*) NLINES
!      DO I = 1, NLINES
                                                                       !comments
!        READ(1,1003) LTEXT
!      ENDDO
                                          !Read no. of levels (max=NDIM2) NLEV,
!      READ (1,*) NLEV, NTEMP
                                          !no. of Te (max=NDIM1) NTEMP and the
!      DO I = 1, NLEV
                                          !input format (cf Readme)
!         READ (1,1002) LABEL(I)
!      ENDDO
!     be
      ibig=0
!Read in Te's where coll. strengths are tabulated
!      DO I = 1, NTEMP
!           READ (1,*) T(I)
!           T(I) = LOG10 (T(I))
!           ROOTT(I) = SQRT(T(I))
!      ENDDO
                             !If IRATS=0, what tabulated are collision strengths
!      READ(1,*) IRATS
!                            !Else Coll. rates = tabulated values * 10 ** IRATS
!      IF(IBIG.EQ.0) THEN
!   10   READ (1,*) ID(2), JD(2), QX
!        IF (QX.EQ.0.D0) GOTO 20
!        IF (ID(2).EQ.0) THEN
!          ID(2) = ID(1)
!          K = K + 1
!        ELSE
!          ID(1) = ID(2)
!          K = 1
!        ENDIF
!        IF (JD(2).EQ.0) THEN
!          JD(2) = JD(1)
!        ELSE
!          JD(1) = JD(2)
!        ENDIF
!        I = ID(2)
!        J = JD(2)
!        QOM(K,I,J) = QX
!        GO TO 10
!      ENDIF
!   20 IF(IBIG.EQ.1.OR.IBIG.EQ.2) THEN
!        READ(1,*) NTRA
!        DO IN = 1, NTRA
!          READ(1,*) I,J,(QOM(ITEMP,I,J),ITEMP=1,NTEMP)
!        ENDDO
!      ENDIF
                                                  !Read transition probabilities
      NLEV1 = NLEV - 1
!      IF (IBIG.EQ.1) THEN
!       READ(1,7000) ((I,J,A(J,I),L=K+1,NLEV),K=1,NLEV1)
!      ELSE
!      DO K = 1, NLEV1
!        KP1 = K + 1
!          DO L = KP1, NLEV
!            READ (1,*) I, J, AX
!            A(J,I) = AX
!          ENDDO
!      ENDDO
!      ENDIF
                                 !Read statistical weights, energy levels (cm-1)
!      DO J = 1, NLEV
!        READ (1,*) I, GX, EX
!        G(I) = GX
!        E(I) = EX
!      ENDDO
!      CLOSE (UNIT=1)
                                !Get levels for ratio
                                                           !150 large enough

                                               !Read in Te and Ne where the line
                                               !ratio is to be calculated

                                                               !Start of Te loop
      DO JT = 1, INT
        TEMP=TEMPI+(JT-1)*TINC

                                                               !Start of Ne loop
        DO JJD = 1, IND
          DENS=DENSI+(JJD-1)*DINC
!          IF(DENSI.LT.30.D0) THEN !commented out as using linear
!          densities only now
!            DENS=10.D0**DENS
!          ENDIF
          IF (TEMP.LE.0.D0.OR.DENS.LE.0.D0) THEN
            WRITE (6,6100)
            STOP
          ENDIF
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
          IF (NTEMP.EQ.1) THEN
            WRITE (6,*)
            WRITE (6,*)                                                 &
     &      'Coll. strengths available for 1 Te only - assuming const'
          ELSEIF (NTEMP.EQ.2) THEN
            WRITE (6,*)
            WRITE (6,*)                                                 &
     &      'Coll. strengths available for 2 Te only - linear interp'
          ELSE
            CALL SPLMAT(T, NTEMP, IOPT, NDIM1, NDIM1T3, HMH)
            CALL CFD(TLOGT,T,NTEMP,NDIM1,HMH,D)
          ENDIF
          DO I = 2, NLEV
            DO J = I, NLEV
                                                           !Negative!
              DELTEK = (E(I-1)-E(J))*1.4388463D0
              EXPE = EXP(DELTEK/TEMP)
              DO IT = 1, NTEMP
                IF (IRATS.EQ.0.D+00) THEN
                  QQ(IT) = QOM(IT,I-1,J)
                ELSE
                                                      !Take out the exp. depend.
                  QQ(IT) = QOM(IT,I-1,J) / EXPE
                                                      !before interpolation
                ENDIF
              ENDDO
              IF (NTEMP.EQ.1) THEN
                 DD = QQ(1)
              ELSEIF (NTEMP.EQ.2) THEN
                 DD = QQ(1) +                                           &
     &            (QQ(2) - QQ(1))/(T(2) - T(1)) * (TLOGT - T(1))
              ELSE
                CALL CFY(TLOGT, DD, T, QQ, NTEMP, NDIM1,HMH,D)
              ENDIF
              IF (IRATS.EQ.0.D+00) THEN
                CS(I-1,J) = DD
              ELSE
                CS(I-1,J) = DD * EXPE
              ENDIF
              IF (IRATS .EQ. 0.D+00) THEN
                QEFF(I-1,J) = 8.63D-06*CS(I-1,J) * EXPE / (G(I-1)*TEMP2)
                QEFF(J,I-1) = 8.63D-06 * CS(I-1,J) / (G(J)*TEMP2)
              ELSE
                QEFF(I-1,J) = CS(I-1,J) * 10. ** IRATS
                                                                     !Be careful
                QEFF(J,I-1) = G(I-1) * QEFF(I-1,J) / (EXPE * G(J))
                                                                     !G integer!
              ENDIF
            ENDDO
          ENDDO
          DO I = 2, NLEV
            DO J = 1, NLEV
              IF (J.NE.I) THEN
                X(I,J) = X(I,J) + DENS * QEFF(J,I)
                X(I,I) = X(I,I) - DENS * QEFF(I,J)
                IF (J.GT.I) THEN
                  X(I,J) = X(I,J) + A(J,I)
                ELSE
                  X(I,I) = X(I,I) - A(I,J)
                ENDIF
              ENDIF
            ENDDO
          ENDDO
          DO I = 2, NLEV
            IM1 = I - 1
            VALUE = 0.D0 - X(I,1)
            Y(IM1) = VALUE
            Y2(IM1) = VALUE
            YKEEP(IM1) = VALUE
            DO J = 2, NLEV
              JM1 = J - 1
              VALUE = X(I,J)
              X(IM1,JM1) = VALUE
              X2(IM1,JM1) = VALUE
              XKEEP(IM1,JM1) = VALUE
            ENDDO
          ENDDO
                                              !Solve matrices for populations
          CALL LUSLV(X,Y,NLEV1,NDIM2)
          DO I = NLEV, 2, -1
            N(I) = Y(I-1)
          ENDDO
          SUMN = 1.D0
          DO I = 2, NLEV
            SUMN = SUMN + N(I)
          ENDDO
          DO I = 2, NLEV
            N(I) = N(I) / SUMN
          ENDDO
          N(1) = 1.D0 / SUMN
                                                                  !Output data
          TTT=TEMP*1.0D-4
          TTP=TTT**(-0.87D0)
                                       !Eff. recombination coef. of Hb
          AHB=3.036D-14*TTP
          DO I = 1, NLEV1
            IP1 = I + 1
            DO J = IP1, NLEV
               IF (A(J,I).NE.0.D0) THEN
                 EJI = E(J) - E(I)
                 WAV = 1.D8 / EJI
                 RLINT = A(J,I) * EJI
                 RLINT = RLINT *N(J)
                 TNIJ(I,J)=RLINT
                 FINT=N(J)*A(J,I)*4861.D0/(DENS*AHB*WAV)
                 FINTIJ(I,J)=FINT
               ENDIF
            ENDDO
          ENDDO
                     !Search ITRANA, ITRANB & ITRANC for transitions & sum up
          SUMA=0.D0
          SUMB=0.D0
          SUMC=0.D0
          IAPR=0
          IBPR=0
          ICPR=0
          DO IKT = 1, NDIM2
            IA1=ITRANA(1,IKT)
            IA2=ITRANA(2,IKT)
            IF(IA1.NE.0.AND.IA2.NE.0) THEN
              SUMA=SUMA+TNIJ(IA1,IA2)
              IAPR=IAPR+1
            ENDIF
            IB1=ITRANB(1,IKT)
            IB2=ITRANB(2,IKT)
            IF(IB1.NE.0.AND.IB2.NE.0) THEN
             IBPR=IBPR+1
             SUMB=SUMB+TNIJ(IB1,IB2)
            ENDIf

            IC1=ITRANC(1,IKT)
            IC2=ITRANC(2,IKT)
            IF(IC1.NE.0.AND.IC2.NE.0) THEN
             ICPR=ICPR+1
             SUMC=SUMC+FINTIJ(IC1,IC2)
            ENDIf
          ENDDO
          FRAT=SUMA/SUMB
          SUMC = 1./SUMC
          TDRAT(1,JJD)=DENS
          TDRAT(2,JJD)=FRAT
          abund = sumc*iobs/100
                                                       !End of the Ne loop
        ENDDO
        DO IA = 1, IAPR
          I1=ITRANA(1,IA)
          I2=ITRANA(2,IA)
          DEE=E(I2)-E(I1)
          WAVA(IA)=1.D8/DEE
        ENDDO
        DO IB = 1, IBPR
          I1=ITRANB(1,IB)
          I2=ITRANB(2,IB)
          DEE=E(I2)-E(I1)
          WAVB(IB)=1.D8/DEE
        ENDDO
        DO IC = 1, ICPR
          I1=ITRANC(1,IC)
          I2=ITRANC(2,IC)
          DEE=E(I2)-E(I1)
          WAVC(IC)=1.D8/DEE
        ENDDO
                                                         !End of the Te loop
      ENDDO

      RETURN

 6100 FORMAT (' PROCESSING COMPLETED'/                                  &
     & ' GOODBYE!!'///)

      END subroutine get_abundance

!---- PROC LUSLV
                                                       !Solving linear equations
      SUBROUTINE LUSLV(A,B,N,M)
      IMPLICIT NONE
      INTEGER M, N
      REAL*8 A(M,M),B(M)
      CALL LURED(A,N,M)
      CALL RESLV(A,B,N,M)
      RETURN
      END subroutine luslv
!
!---- PROC LURED
      SUBROUTINE LURED(A,N,NR)
      IMPLICIT NONE
      INTEGER N, NR, NM1, I, J, K, IP1
      REAL*8 A(NR,NR), FACT
      IF(N.EQ.1) RETURN
      NM1=N-1
      DO I=1,NM1
        IP1=I+1
        DO K=IP1,N
          FACT=A(K,I)/A(I,I)
          DO J=IP1,N
            A(K,J)=A(K,J)-A(I,J)*FACT
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END subroutine lured
!
!---- PROC RESLV
                                                               !Resolve A with B
      SUBROUTINE RESLV(A,B,N,NR)
      IMPLICIT NONE
      INTEGER N, NR, NM1, I, J, K, L, IP1
      REAL*8 A(NR,NR),B(NR)
      IF(N.EQ.1) GOTO 1
      NM1=N-1
      DO I=1,NM1
        IP1=I+1
        DO J=IP1,N
          B(J)=B(J)-B(I)*A(J,I)/A(I,I)
        ENDDO
      ENDDO
      B(N)=B(N)/A(N,N)
      DO I=1,NM1
        K=N-I
        L=K+1
        DO J=L,N
          B(K)=B(K)-B(J)*A(K,J)
        ENDDO
        B(K)=B(K)/A(K,K)
      ENDDO
      RETURN
    1 B(N)=B(N)/A(N,N)
      RETURN
      END subroutine reslv
!
!---- PROC SPLMAT
      SUBROUTINE SPLMAT(XX,NPT,IOPT,NDIM, NDIMT3, HMH)
      IMPLICIT NONE
      INTEGER NDIM, NDIMT3, NPT, IOPT, NPM, NELEM
      REAL*8 XX(NDIM),GH(NDIMT3),Y(NDIM), HMH(NDIM,NDIM)
      NPM=NPT-2
      CALL GHGEN(GH,XX,NPT,IOPT,NDIM,NDIMT3)
      NELEM=3*NPM-2
      CALL ELU(GH,NPM,NDIMT3)
      CALL HGEN(XX,GH,Y,NPT,IOPT,NDIM,NDIMT3,HMH)
      RETURN
      END subroutine splmat
!
!---- PROC DERIV
!     Calculate the first derivative of the lagrangian interpolator
!     of a function F, tabulated at the N points XY(I), I=1 to N.
!     The derivative is given as the coefficients of F(I), I=1 to N,
!     in the array D(I), I=1 to N.
      SUBROUTINE DERIV(XY,D,X,N,NDIM)
      IMPLICIT NONE
      INTEGER N ,NDIM, I, J, K
      REAL*8 XY(NDIM),D(NDIM), X, P1, P2, S
      DO I=1,N
        P1=1.
        S=0.
        DO J=1,N
          IF(J.NE.I) THEN
            P1=P1*(XY(I)-XY(J))
            P2=1.
            DO K=1,N
              IF(K.NE.I.AND.K.NE.J) P2=P2*(X-XY(K))
            ENDDO
            S=S+P2
          ENDIF
        ENDDO
        D(I)=S/P1
      ENDDO
      RETURN
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
      SUBROUTINE HGEN(XX,GH,Y,NPT,IOPT,NDIM,NDIMT3,HMH)
      IMPLICIT NONE
      INTEGER NPT, IOPT, NDIM, NDIMT3, NDIM3, NIP, I, J, K, NPM,        &
     & INDX
      REAL*8 XX(NDIM), GH(NDIMT3), Y(NDIM), HMH(NDIM,NDIM),             &
     & XY(5),D(5),C(2,5), A0, AN1, H1, H2
                                   !Case of derivative boundary condition, with
      IF(IOPT.EQ.2) THEN
                                   !derivatives from NIP-point Lagrange at
        NDIM3=5
                                   !internal points
        NIP=3
        DO J=1,2
          DO I=1,NIP
            K=(NPT-NIP)*(J-1)
            XY(I)=XX(K+I)
          ENDDO
          K=1+(NPT-1)*(J-1)
          CALL DERIV(XY,D,XX(K),NIP,NDIM3)
          DO I=1,NIP
            C(J,I)=D(I)
          ENDDO
        ENDDO
      ENDIF
                                             !Set up matrix equation G*YPP=HMH*Y
      A0=XX(2)-XX(1)
      AN1=XX(NPT)-XX(NPT-1)
      NPM=NPT-2
      hmh=0d0
      DO I=1,NPM
        H1=6./(XX(I+1)-XX(I))
        H2=6./(XX(I+2)-XX(I+1))
        HMH(I,I)=H1
        HMH(I,I+1)=-H1-H2
        HMH(I,I+2)=H2
!        DO J=1,NPT
!          HMH(I,J)=0.
!          IF(J.EQ.I) HMH(I,J)=H1
!          IF(J.EQ.I+2) HMH(I,J)=H2
!          IF(J.EQ.I+1) HMH(I,J)=-H1-H2
!        ENDDO
      ENDDO
                                                 !Correct matrix for case of
      IF(IOPT.EQ.1.OR.IOPT.EQ.2) THEN
                                                 !derivative boundary conditions
        HMH(1,1)=HMH(1,1)+3/A0
        HMH(1,2)=HMH(1,2)-3/A0
        HMH(NPM,NPT-1)=HMH(NPM,NPT-1)-3/AN1
        HMH(NPM,NPT)=HMH(NPM,NPT)+3/AN1
      ENDIF
      IF(IOPT.EQ.2) THEN
        DO J=1,NIP
          HMH(1,J)=HMH(1,J)+3*C(1,J)
          K=NPT+J-NIP
          HMH(NPM,K)=HMH(NPM,K)-3*C(2,J)
        ENDDO
      ENDIF
                                 !Solve matrix equation with results in the form
      DO I=1,NPT
                                 !YPP=HMH*Y. matrix g has been LU decomposed
        Y(1)=HMH(1,I)
        INDX=0
        DO J=2,NPM
          INDX=INDX+3
          Y(J)=HMH(J,I)-GH(INDX)*Y(J-1)
        ENDDO
        INDX=INDX+1
        Y(NPM)=Y(NPM)/GH(INDX)
        DO J=2,NPM
          K=NPM-J+1
          INDX=INDX-3
          Y(K)=(Y(K)-GH(INDX+1)*Y(K+1))/GH(INDX)
        ENDDO

        HMH(2:npm+1,I)=Y(1:npm)
!        DO J=1,NPM
!        HMH(J+1,I)=Y(J)
!        ENDDO
                                    !Insert values for second derivative at end
        HMH(1,I)=0.
                                    !points: first and last rows of the matrix
        HMH(NPT,I)=0.
      ENDDO
                                         !Case of derivative boundary conditions
      IF(IOPT.GT.0) THEN
        DO J=1,NPT
          HMH(1,J)=-0.5*HMH(2,J)
          HMH(NPT,J)=-0.5*HMH(NPT-1,J)
        ENDDO
        HMH(1,1)=HMH(1,1)-3/(A0*A0)
        HMH(1,2)=HMH(1,2)+3/(A0*A0)
        HMH(NPT,NPT-1)=HMH(NPT,NPT-1)+3/(AN1*AN1)
        HMH(NPT,NPT)=HMH(NPT,NPT)-3/(AN1*AN1)
      ENDIF
      IF(IOPT.EQ.2) THEN
        DO J=1,NIP
          HMH(1,J)=HMH(1,J)-3*C(1,J)/A0
          K=NPT+J-NIP
          HMH(NPT,K)=HMH(NPT,K)+3*C(2,J)/AN1
        ENDDO
      ENDIF
      RETURN
      END subroutine hgen
!
!---- PROC GHGEN
      SUBROUTINE GHGEN(GH,XX,NPT,IOPT,NDIM,NDIMT3)
      IMPLICIT NONE
      INTEGER NPT, IOPT, NDIM, NDIMT3, INDX, NPTM, I, J, IP, JP, IK
      REAL*8 XX(NDIM),GH(NDIMT3)
      INDX=0
      NPTM=NPT-1
      DO I=2,NPTM
        IP=I-1
        DO J=1,3
          JP=IP+J-2
          IF(JP.GE.1.AND.JP.LE.NPTM-1) THEN
            INDX=INDX+1
            IF(J.EQ.2) THEN
              GH(INDX)=2*(XX(I+1)-XX(I-1))
            ELSE
              IK=I+(J-1)/2
              GH(INDX)=XX(IK)-XX(IK-1)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      IF(IOPT.GE.1) THEN
        GH(1)=GH(1)-(XX(2)-XX(1))/2.
        GH(INDX)=GH(INDX)-(XX(NPT)-XX(NPT-1))/2.
      ENDIF
      RETURN
      END subroutine ghgen
!
!---- PROC ELU
      SUBROUTINE ELU(GH,N,NDIM)
      IMPLICIT NONE
      INTEGER N, NDIM, INDX, I, J, JP
      REAL*8 GH(NDIM)
      INDX=0
      DO I=1,N
        DO J=1,3
          JP=I+J-2
          IF(JP.GE.1.AND.JP.LE.N) THEN
            INDX=INDX+1
            IF(I.GT.1) THEN
              IF(J.EQ.1) THEN
                GH(INDX)=GH(INDX)/GH(INDX-2)
              ENDIF
              IF(J.EQ.2) THEN
                GH(INDX)=GH(INDX)-GH(INDX-1)*GH(INDX-2)
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      RETURN
      END subroutine elu
!
!---- PROC CFY
      SUBROUTINE CFY(X,Y,XX,YY,NPT,NDIM, HMH, D)
      IMPLICIT NONE
      INTEGER NPT, NDIM, J
      REAL*8 XX(NDIM),YY(NDIM), HMH(NDIM,NDIM), D(NDIM),                &
     & X, Y, TT
      IF(X.LT.XX(1)) THEN
        Y=YY(1)
      ENDIF
      IF(X.GT.XX(NPT)) THEN
        Y=YY(NPT)
      ENDIF
      TT=0.
      DO J=1,NPT
        TT=TT+D(J)*YY(J)
      ENDDO
      Y=TT
      RETURN

      END subroutine cfy
!
!---- PROC CFD
      SUBROUTINE CFD(X,XX,NPT,NDIM, HMH, D)
      IMPLICIT NONE
      INTEGER NPT, NDIM, NPTM, I, J
      REAL*8 X, XX(NDIM), HMH(NDIM,NDIM), D(NDIM), X1, X2, A1, A2, HI
      IF(X.LT.XX(1)) THEN
        !WRITE(6,400) XX(1)
        RETURN
      ENDIF
      IF(X.GT.XX(NPT)) THEN
        !WRITE(6,401) XX(NPT)
        RETURN
      ENDIF
      NPTM=NPT-1
      DO I=1,NPTM
        IF(X.LT.XX(I+1)) THEN
          X1=XX(I+1)-X
          X2=X-XX(I)
          HI=XX(I+1)-XX(I)
          A1=X1*(X1*X1/(6*HI)-HI/6)
          A2=X2*(X2*X2/(6*HI)-HI/6)
          DO J=1,NPT
            D(J)=(A1*HMH(I,J)+A2*HMH(I+1,J))
          ENDDO
          D(I)=D(I)+X1/HI
          D(I+1)=D(I+1)+X2/HI
          RETURN
        ENDIF
      ENDDO

      END subroutine cfd

      end module mod_equib
