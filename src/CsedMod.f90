    MODULE CSedMod
    use RivVarMod
    use CalcCond
    IMPLICIT NONE
    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:) :: QTOT
    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:) :: RATIO, d, DUM1
    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:) :: DHDN, DHDS, USTRESS, VSTRESS
    !	REAL :: rsin, rcos, ux, uy, uxnew, uynew
    !	REAL, ALLOCATABLE, DIMENSION(:,:) :: CON
    CONTAINS
    SUBROUTINE alloc_csed()
    !		INTEGER, INTENT(IN) :: ns, nn, nz
    INTEGER :: status
    !Allocate REAL types
    ALLOCATE(QTOT(ns), STAT = status)
    ALLOCATE(RATIO(ns, nn), STAT = status)
    ALLOCATE(d(ns, nn), STAT = status)
    !			ALLOCATE(QS(ns, nn), STAT = status)
    !			ALLOCATE(QN(ns, nn), STAT = status)
    ALLOCATE(DUM1(ns, nn), STAT = status)
    ALLOCATE(DHDN(ns, nn), STAT = status)
    ALLOCATE(DHDS(ns, nn), STAT = status)
    ALLOCATE(USTRESS(ns, nn), STAT = status)
    ALLOCATE(VSTRESS(ns, nn), STAT = status)
    !			ALLOCATE(CON(ns, nn), STAT = status)
    END SUBROUTINE alloc_csed

    SUBROUTINE dealloc_csed()
    INTEGER :: status
    !Allocate REAL types
    DEALLOCATE(QTOT, STAT = status)
    DEALLOCATE(RATIO, STAT = status)
    DEALLOCATE(d, STAT = status)
    !			DEALLOCATE(QS, STAT = status)
    !			DEALLOCATE(QN, STAT = status)
    DEALLOCATE(DUM1, STAT = status)
    DEALLOCATE(DHDN, STAT = status)
    DEALLOCATE(DHDS, STAT = status)
    DEALLOCATE(USTRESS, STAT = status)
    DEALLOCATE(VSTRESS, STAT = status)
    END SUBROUTINE dealloc_csed
    !


    SUBROUTINE csed(DT,nct,nsteps,newdt)
    REAL(KIND=mp), INTENT(IN) :: DT
    INTEGER, INTENT(IN) :: nct, nsteps
    REAL(KIND=mp), INTENT(INOUT) :: newdt
    !		parameter(ns=41,nn=25,nz=11)
    !	   REAL mo,hav
    !	   common taus(ns,nn),taun(ns,nn),hl(ns,nn),eta(ns,nn),rn(ns,nn),
    !	  & ibc(ns,nn),r(ns),w(ns),cd,ds,dn,mo,q,hwt,urelax,erelax,
    !	  & arelax,itm,hav,iplinc,xo(ns),yo(ns),uz(ns,nn,nz),vz(ns,nn,nz)
    !	dimension RATIO(ns,nn),d(ns,nn),
    !	 &	  QS(ns,nn),QN(ns,nn),DUM1(ns,nn),
    !	 &	  DHDN(ns,nn),DHDS(ns,nn),USTRESS(ns,nn),
    !	 &	  VSTRESS(ns,nn),QTOT(ns),CON(ns,nn)
    CHARACTER(LEN=20) RUNID
    REAL(KIND=mp) :: CDune, VKC, GSW, PI, GAMC, TB, TC, G, RHO
    REAL(KIND=mp) :: RAT, ZOSF, RHOFAC, DUM, GAMG, TAUG
    REAL(KIND=mp) :: A2, DTAU, TAN1, TAN2, QTOTAL, SC, weight
    REAL(KIND=mp) :: fac
    REAL(KIND=mp) :: PLINC
    INTEGER :: NM, I, J, ityp, ISMOO

    REAL(KIND=mp) :: mintdt, maxtdt, tdt
    REAL(KIND=mp) :: tmpdepth
    REAL(KIND=mp) :: tdd
    REAL(KIND =mp) :: t, csf, p, s, ws, kv

    !			CALL alloc_csed()
    runid='test         '
    !OPEN(7,FILE='seddat')
    !OPEN(10,FILE='top2')
    !			OPEN(11,FILE='topser')
    !			DIN=0.02

    CDune=.212
    VKC=0.4
    GSW=1.
    NM=(NN+1)/2
    PI=ACOS(-1.)
    GAMC=SUBANGLEOFREPOSE*PI/180.
    TC=(200.*(DIN**2))+(26.*DIN)+1.
    G=980.
    RHO=1.
    T = 20.
    kv = 1.79D-6/(1.0D0+3.37D-2*t+2.21D-5*t*t) !MKS  !Kinematic Viscosity

    !	 WRITE(6,*) 'ENTER THE DUNE HEIGHT IN CM (0 FOR PLANE BED):'
    !	 READ(5,*) HD
    !			hd=0.
    IF(HD.LE.0) THEN
        RAT=1.
    ELSE
        !			WRITE(6,*) 'ENTER THE AVERAGE DUNE WAVELENGTH IN CM:'
        !			READ(5,*) WD
        IF(WD.LE.0) THEN
            WD=1.
        ENDIF
        ZOSF=.2*DIN
        RAT=1.+(CDune/(2.*VKC**2.))*(HD/WD)*((LOG(HD/(ZOSF)))-1.)**2.
    ENDIF
    DO I=1,ns
        DO J=1,nn
            D(I,J)=DIN
            !	 hd=0.2*h(i,j)
            !	 wd=20*hd
            !	 RAT=1.+(CDune/(2.*VKC**2.))*(HD/WD)*((ALOG(HD/(ZOSF)))-1.)**2.
            RATIO(I,J)=RAT
            USTRESS(I,J)=taus(I,J)/RATIO(I,J)
            VSTRESS(I,J)=taun(I,J)/RATIO(I,J)
        ENDDO
    ENDDO
    !	 write(6,*) 'Yalin (0) or Engelund-Hansen (1):'
    !	 read(5,*) ityp
    write(6,*) 'itype', TRANSEQTYPE, 'sedsmoo',SEDSMOO,'din',din
    !			ityp=0
    RHOFAC=1.65
    DO I=1,ns
        DO J=1,nn
            IF(I.EQ.1) THEN
                DHDS(I,J)=(hl(2,J)-hl(1,J))/(ds*(rn(I,J)))
            ELSE
                IF(I.EQ.ns) THEN
                    DHDS(I,J)=(hl(ns,J)-hl(ns-1,J))/(ds*(RN(I,J)))
                ELSE
                    DHDS(I,J)=(hl(I+1,J)-hl(I,J))/(ds*(RN(I,J)))
                ENDIF
            ENDIF
            IF(J.EQ.1) THEN
                !       DHDN(I,J)=(hl(I,2)-hl(I,1))/dn
                DHDN(I,J) = 0.
            ELSE
                IF(J.EQ.nn) THEN
                    !        DHDN(I,J)=(hl(I,nn)-hl(I,nn-1))/dn
                    DHDN(I,J) = 0.
                ELSE
                    DHDN(I,J)=(hl(I,J+1)-hl(I,J))/dn
                ENDIF
            ENDIF
            ! now fix things for parts sticking out of channel
            !     do i=1,ns
            !        do j=2,nn-1

            IF(j.ne.1.and.j.ne.nn) THEN
                if(ibc(i,j+1).eq.0.and.ibc(i,j-1).eq.0) then
                    dhdn(i,j) = 0.
                else
                    if(ibc(i,j-1)==0) then !I think this statement below is wrong.  Need to check.
                        !!                    dhdn(i,j) = abs(hl(i,j+1)-hl(i,j))/dn
                        !
                        dhdn(i,j) = -1*(eta(i,j+1)-eta(i,j))/dn
                    endif
                    if(ibc(i,j+1)==0) then
                        !                    dhdn(i,j) = (hl(i,j)-hl(i,j-1))/dn
                        dhdn(i,j) = -1*(eta(i,j)-eta(i,j-1))/dn
                    endif
                endif
            ENDIF
            !        enddo
            !     enddo
            DUM=((DHDS(I,J)**2)+(DHDN(I,J)**2))**.5
            GAMG=ATAN(DUM)
            Taug = 0
            if(CALCGRAVCORR == 1) then
                TAUG=1.0*TC*sin(GAMG)/sin(GAMC)
                TAN2= 0.
            else
                TAUG = 0.
                TAN2 = 0.
            endif

            ! rmcd hard wire taug=0 at j=1 and j=nn for flume walls
            if(j == 1 .or. j == nn .or. ibc(i,j) == 0) then
                taug = 0.
            endif
            ! end PAN change

            IF(DUM.NE.0) THEN
                USTRESS(I,J)=USTRESS(I,J)+(TAUG*(DHDS(I,J)/DUM))
                VSTRESS(I,J)=VSTRESS(I,J)+(TAUG*(DHDN(I,J)/DUM))
            ENDIF
        ENDDO
    ENDDO
    !612 CONTINUE

    !!!!!  YALIN EQUATION   !!!!!
    if(TRANSEQTYPE == 0) then !Yalin
        DO  I=1,ns
            DO  J=1,nn
                TC=(200.*(D(I,J)**2.))+(26.*D(I,J))+1.
                tdd = D(I,J)/100. !convert to MKS
                CALL tcrit2(tdd, kv, TC) !MKS
                TC = TC*10 !Convert to CGS
                A2=1.66*((TC/(1617.*D(I,J)))**.5)
                TB=((USTRESS(I,J)**2)+(VSTRESS(I,J)**2))**.5
                IF (TB.LE.TC.or.ibc(i,j).eq.0) THEN
                    QS(I,J)=0.
                    QN(I,J)=0.
                ELSE
                    DTAU=(TB-TC)/TC
                    DUM=1.-(LOG(1.+(A2*DTAU)))/(A2*DTAU)
                    DUM=.635*D(I,J)*DUM*DTAU
                    QS(I,J)=(DUM*((TB)**.5))
                    QS(I,J)=QS(I,J)*USTRESS(I,J)/(TB)
                    TAN1=VSTRESS(I,J)/USTRESS(I,J)
                    if(CALCGRAVCORR == 2) then
                        TAN2=GRAVFLATBEDCORRCOEF*((TC/TB)**.5)*DHDN(I,J)
                        QN(I,J)=QS(I,J)*(TAN1+TAN2)

                    else
                        QN(I,J)=QS(I,J)*TAN1

                    ENDIF

                    !			write(6,*) i, j, tan1, tan2, qn(i,j)
                    IF(J.EQ.1.OR.J.EQ.nn) THEN
                        QN(I,J)=0.
                    ENDIF
                ENDIF
            ENDDO
        ENDDO
    ELSEIF (transeqtype == 1) then !Engelund and Hanssen
        do i=1,ns
            do j=1,nn
                TC=(200.*(D(I,J)**2.))+(26.*D(I,J))+1.
                tdd = D(I,J)/100. !convert to MKS
                CALL tcrit2(tdd, kv, TC) !MKS
                TC = TC*10 !Convert to CGS
                TB=((USTRESS(I,J)**2)+(VSTRESS(I,J)**2))**.5
                IF (TB.LE.TC.or.ibc(i,j).eq.0) THEN
                    QS(I,J)=0.
                    QN(I,J)=0.
                else
                    fac = ((1.65*980.0*d(i,j)**3.)**0.5)*cd(i,j)
                    qs(i,j)=fac*(tb/tc)**2.5
                    qs(i,j)=qs(i,j)*ustress(i,j)/TB
                    TAN1 = VSTRESS(I,J)/USTRESS(I,J)
                    ! if(CalcGravCorr == 1 .AND. GRAVCORRTYPE == 2) then !pan - gravcorrtype doesn't exist anymore.
                    if(CalcGravCorr == 2) then
                        TAN2 = gravflatbedcorrcoef*((tc/tb)**.5)*dhdn(i,j)
                        QN(I,J) = QS(I,J)*(TAN1+TAN2)
                    else
                        QN(I,J) = QS(I,J)*TAN1
                    endif
                    if(j.eq.1.or.j.eq.nn) then
                        qn(i,j)=0.
                    endif
                endif
            enddo
        enddo
    ENDIF
    !	set upstream sed flux (can change later to allow recirculation)
    DO J = 1,nn
        QS(SEDBCNODE, J) = BCFRACTION*QS(SEDBCNODE, J)
        QN(SEDBCNODE, J) = BCFRACTION*QN(SEDBCNODE, J)
    ENDDO

710 QTOTAL=0.
    DO I=1,ns
        QTOT(I)=0.
        DO J=2,nn
            QTOT(I)=QTOT(I)+(.5*dn)*(QS(I,J)+QS(I,J-1))
        ENDDO
        QTOTAL=QTOTAL+QTOT(I)/ns
    ENDDO


    DO  I=1,ns
        DO  J=1,nn
            SC=rn(I,J)
            IF(I.eq.1) THEN
                CON(I,J)=0.
                !			        CON(I,J)=(QS(2,J)-QS(1,J))/(ds*SC)
                !	              GO TO 643
            ELSEIF(I.EQ.ns) THEN
                !	              CON(I,J)=(QS(I,J)-QS(I-1,J))/(ds*SC)
                !	              CON(I,J)=(QS(2,J)-QS(1,J))/(ds*SC)
                CON(I,J)=CON(I-1,J)
                !			         CON(I,J)=CON(I,J)-(QN(I,J)/(R(I)*SC))
            ELSE
                CON(I,J)=(QS(I+1,J)-QS(I-1,J))/(2.*ds*SC)
            ENDIF

            CON(I,J)=CON(I,J)-(QN(I,J)/(R(I)*SC))

            IF(J.EQ.1) THEN
                !	        CON(I,1)=CON(I,1)+(QN(I,2)-QN(I,1))/dn
                con(i,1)=0.
            ELSEIF(J.EQ.nn) THEN
                !	         CON(I,nn)=CON(I,nn)+(QN(I,nn)-QN(I,nn-1))/dn
                con(i,nn)=0.
            ELSE
                !first two if statements below worked well with parabolic bed
                IF(IBC(I,J-1).ne.-1) THEN
                    CON(I,J)=CON(I,J)+(QN(I,J)-QN(I,J-1))/(dn)
                    !	        ELSEIF(IBC(I,J+1).ne.-1)THEN
                    !	            CON(I,J)=CON(I,J)+(QN(I,J+1)-QN(I,J))/(dn)
                ELSE
                    CON(I,J)=CON(I,J)+(QN(I,J+1)-QN(I,J-1))/(2.*dn)
                ENDIF
            ENDIF
            if(i.lt.sedbcnode) then
                con(i,j) = 0.
            endif

        ENDDO
    ENDDO


    !			do 652 j=1,nn
    !652			con(1,j)=con(2,j)
    DO I=1,ns
        DO J=1,nn
            CON(I,J)=1.54*CON(I,J)
        ENDDO
    ENDDO

    !           write(6,*) 'con',(con(10,j),j=1,nn)
    !             write(6,*) 'rn',(rn(10,j),j=1,nn)
    IF (SEDSMOO.ge.1) then
        DO ISMOO=1,SEDSMOO
            DO I=2,ns-1
                DO J=1,nn
                    if(ibc(i,j).eq.0) then
                        dum1(i,j)=0.
                    elseif(j.eq.1) then
                        DUM1(I,1)=(CON(I-1,1)+CON(I,1)+CON(I+1,1)+con(I,j+1))/4.
                        !             DUM1(I,1)=(.5*CON(I-1,1)+CON(I,1)+.5*CON(I+1,1))/2.
                    elseif(j.eq.nn) then
                        DUM1(I,nn)=(CON(I-1,nn)+CON(I,nn)+CON(I+1,nn)+con(i,nn-1))/4.
                        !             DUM1(I,nn)=(.5*CON(I-1,nn)+CON(I,nn)+.5*CON(I+1,nn))/2.
                    else
                        weight=1
                        if(con(i,j-1).ne.0) then
                            weight=weight+SEDSMOOWGHT
                        endif
                        if(con(i,j+1).ne.0) then
                            weight=weight+SEDSMOOWGHT
                        endif
                        if(con(i-1,j).ne.0) then
                            weight=weight+SEDSMOOWGHT
                        endif
                        if(con(i+1,j).ne.0) then
                            weight=weight+SEDSMOOWGHT
                        endif
                        DUM1(I,J)=1.*CON(I,J)+SEDSMOOWGHT*CON(I+1,J)+SEDSMOOWGHT*CON(I-1,J)
                        DUM1(I,J)=DUM1(I,J)+SEDSMOOWGHT*CON(I,J-1)+SEDSMOOWGHT*CON(I,J+1)
                        DUM1(I,J)=DUM1(I,J)/weight
                    endif
                ENDDO
            ENDDO
            !            write(6,*) 'smooth pass', ismoo, (dum1(10,j),j=1,nn)
            DO I=2,ns-1
                DO J=1,nn
                    CON(I,J)=DUM1(I,J)
                ENDDO
            ENDDO

        ENDDO
    endif
    !            write(6,*) 'con',(con(10,j),j=1,nn)
    !            write(6,*) 'rn',(rn(10,j),j=1,nn)
    !            pause
    !	 DO 685 ISMOO=1,2
    !	 DO 680 I=2,ns-1
    !	 DUM1(I,1)=(CON(I-1,1)+CON(I,1)+CON(I+1,1)+con(I,j+1))/4.
    !680	 DUM1(I,25)=(CON(I-1,nn)+CON(I,nn)+CON(I+1,nn)+con(i,nn-1))/4.
    !	 DO 685 I=2,ns-1
    !	 CON(I,1)=DUM1(I,1)
    !685	  CON(I,nn)=DUM1(I,nn)
    DO  J=1,nn
        !			CON(ns,J)=(.2*CON(ns-1,J)+CON(ns,J))/1.2
        con(ns,j)= con(ns-1,j)
        CON(1,J)=0.
    ENDDO

    !     FIND MINIMUM TIME STEP TO SATISFY CRITERIA FOR MAXIMUM ELEVATION CHANGE
    !     AS A FUNTION OF DEPTH !9/24/06 - rmcd
    mintdt = 1e12
    maxtdt = -1e12
    DO I = 1,ns
        DO J = 1, nn
            IF(ibc(i,j).ne.0.and.con(i,j).ne.0) THEN
                tdt = abs((TSFracDepth*hl(i,j))/con(i,j))
                IF(tdt < mintdt) THEN
                    mintdt = tdt
                ENDIF
                IF(tdt > maxtdt) THEN
                    maxtdt = tdt
                ENDIF
            ENDIF
        ENDDO
    ENDDO
    !            WRITE(6,*) 'Maxdt', maxtdt, 'Mindt', mintdt
    newdt = 1.0*mintdt

    DO I=1,ns
        DO J=1,nn
            if(ibc(i,j).eq.0) then
                con(i,j)=0.
            endif
            if(dt*con(i,j) > hl(i,j)) then
                write(6,*)'Error in Sediment Transport', i, j, hl(i,j)/con(i,j)
            ENDIF

            if(i.lt.sedbcnode) then
                con(i,j) = 0
            endif
            tmpdepth = hl(i,j)+(dt*con(i,j))
            !			if(tmpdepth < hmin) then
            !			    eta(i,j) = e(i,j)-hmin
            !			    hl(i,j) = hmin
            !			else
            eta(i,j)=eta(i,j)-(dt*con(i,j))
            hl(I,J)=hl(I,J)+(dt*con(I,J))
            !			endif
            !			eta(i,j)=eta(i,j)-(dt*con(i,j))
            !			hl(I,J)=hl(I,J)+(dt*con(I,J))
            if(hl(i,j).le.hmin) then
                hl(i,j)=hmin
            endif
        ENDDO
    ENDDO

    !          if(roughnesstype == 1) then
    !		    CALL ZOTOCDTWO()
    !          endif

    !			    Write(11,2010) RUNID
    !			    write(11,*) Q,MO
    !			    write(11,*) R,W,eta
    !			    write(11,*) xo,yo
2000 if(nct.eq.nsteps) then
        !WRITE(7,2010) RUNID
        !WRITE(7,*) QS,QN,CON,hl,QTOT
        !write(7,*) ((i,j,con(i,j),j=1,nn),i=10,15)
        !write(10,*) hwt,erelax,urelax,arelax,itm,cd,q,mo
        !write(10,2010) runid
        !write(10,*) q,mo,hav
        !write(10,*) r,w,eta,ibc
        !write(10,*) xo,yo
    endif
2010 format(A20)
    return
    end subroutine csed

    SUBROUTINE TCRIT2(D, KV, TC)
    !! Based on Shields 1984 all units in MKS.
    IMPLICIT NONE
    REAL(KIND=mp), INTENT(IN) :: D, KV
    REAL(KIND=mp), INTENT(OUT) :: TC

    REAL(KIND=mp) :: DSTAR, TAUCRIT
    REAL(KIND=mp) :: RELDEN, G

    RELDEN = 1.650
    G = 9.81

    DSTAR = D*((G*RELDEN)/(KV*KV))**.333

    IF(DSTAR.LE.4) THEN
        TAUCRIT = 0.24/DSTAR
    ELSE IF (DSTAR.GT.4.AND.DSTAR.LE.10) THEN
        TAUCRIT = 0.15/(DSTAR**0.64)
    ELSE IF (DSTAR.GT.10.AND.DSTAR.LE.20) THEN
        TAUCRIT = 0.04/(DSTAR**0.10)
    ELSE IF (DSTAR.GT.20.AND.DSTAR.LE.150) THEN
        TAUCRIT = 0.013*(DSTAR**0.29)
    ELSE
        TAUCRIT = 0.055
    ENDIF
    TC = TAUCRIT*1650*G*D
    END SUBROUTINE
    END MODULE
