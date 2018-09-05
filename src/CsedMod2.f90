    MODULE CSedMod2
    use RivVarMod2
    use CalcCond2
    IMPLICIT NONE
    type csed
        REAL(KIND=mp), ALLOCATABLE, DIMENSION(:) :: QTOT
        REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:) :: RATIO, d, DUM1
        REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:) :: DHDN, DHDS, USTRESS, VSTRESS
    end type
    !	REAL :: rsin, rcos, ux, uy, uxnew, uynew
    !	REAL, ALLOCATABLE, DIMENSION(:,:) :: CON
    CONTAINS
    SUBROUTINE alloc_csed(this, ns, nn)
    implicit none
    type(csed), intent(inout) :: this
    INTEGER, INTENT(IN) :: ns, nn
    INTEGER :: status
    !Allocate REAL types
    ALLOCATE(this%QTOT(ns), STAT = status)
    ALLOCATE(this%RATIO(ns, nn), STAT = status)
    ALLOCATE(this%d(ns, nn), STAT = status)
    !			ALLOCATE(QS(ns, nn), STAT = status)
    !			ALLOCATE(QN(ns, nn), STAT = status)
    ALLOCATE(this%DUM1(ns, nn), STAT = status)
    ALLOCATE(this%DHDN(ns, nn), STAT = status)
    ALLOCATE(this%DHDS(ns, nn), STAT = status)
    ALLOCATE(this%USTRESS(ns, nn), STAT = status)
    ALLOCATE(this%VSTRESS(ns, nn), STAT = status)
    !			ALLOCATE(CON(ns, nn), STAT = status)
    END SUBROUTINE alloc_csed

    SUBROUTINE dealloc_csed(this)
    implicit none
    type(csed), intent(inout) :: this
    INTEGER :: status
    !Allocate REAL types
    DEALLOCATE(this%QTOT, STAT = status)
    DEALLOCATE(this%RATIO, STAT = status)
    DEALLOCATE(this%d, STAT = status)
    !			DEALLOCATE(QS, STAT = status)
    !			DEALLOCATE(QN, STAT = status)
    DEALLOCATE(this%DUM1, STAT = status)
    DEALLOCATE(this%DHDN, STAT = status)
    DEALLOCATE(this%DHDS, STAT = status)
    DEALLOCATE(this%USTRESS, STAT = status)
    DEALLOCATE(this%VSTRESS, STAT = status)
    END SUBROUTINE dealloc_csed
    !


    SUBROUTINE csed2(this, rvo, cco, DT,nct,nsteps,newdt)
    implicit none
    type(csed), intent(inout) :: this
    type(rivvar), intent(inout) :: rvo !rivar mod variables
    type(calccond), intent(inout) :: cco !calc condition mod variables
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
    REAL(KIND=mp) :: CDune, GSW, PI, GAMC, TB, TC, G, RHO
    REAL(KIND=mp) :: RAT, ZOSF, RHOFAC, DUM, GAMG, TAUG
    REAL(KIND=mp) :: A2, DTAU, TAN1, TAN2, QTOTAL, SC, weight
    REAL(KIND=mp) :: fac
    REAL(KIND=mp) :: PLINC
    INTEGER :: NM, I, J, ityp, ISMOO

    REAL(KIND=mp) :: mintdt, maxtdt, tdt
    REAL(KIND=mp) :: tmpdepth
    REAL(KIND=mp) :: tdd
    REAL(KIND =mp) :: t, csf, p, s, ws, kv
    REAL(kind = mp) :: vkc = 0.4
    integer :: tns, tnn
    !			CALL alloc_csed()
    runid='test         '
    tns = rvo%ns
    tnn = rvo%nn
    !OPEN(7,FILE='seddat')
    !OPEN(10,FILE='top2')
    !			OPEN(11,FILE='topser')
    !			DIN=0.02

    CDune=.212
    !VKC=0.4
    GSW=1.
    NM=(tnn+1)/2
    PI=ACOS(-1.)
    GAMC=cco%SUBANGLEOFREPOSE*PI/180.
    TC=(200.*(cco%DIN**2))+(26.*cco%DIN)+1.
    G=980.
    RHO=1.
    T = 20.
    kv = 1.79D-6/(1.0D0+3.37D-2*t+2.21D-5*t*t) !MKS  !Kinematic Viscosity

    !	 WRITE(6,*) 'ENTER THE DUNE HEIGHT IN CM (0 FOR PLANE BED):'
    !	 READ(5,*) HD
    !			hd=0.
    IF(cco%HD.LE.0) THEN
        RAT=1.
    ELSE
        !			WRITE(6,*) 'ENTER THE AVERAGE DUNE WAVELENGTH IN CM:'
        !			READ(5,*) WD
        IF(cco%WD.LE.0) THEN
            cco%WD=1.
        ENDIF
        ZOSF=.2*cco%DIN
        RAT=1.+(CDune/(2.*VKC**2.))*(cco%HD/cco%WD)*((LOG(cco%HD/(ZOSF)))-1.)**2.
    ENDIF
    DO I=1,tns
        DO J=1,tnn
            this%d(i,j)=cco%DIN
            !	 hd=0.2*h(i,j)
            !	 wd=20*hd
            !	 RAT=1.+(CDune/(2.*VKC**2.))*(HD/WD)*((ALOG(HD/(ZOSF)))-1.)**2.
            this%RATIO(I,J)=RAT
            this%USTRESS(I,J)=rvo%taus(I,J)/this%RATIO(I,J)
            this%VSTRESS(I,J)=rvo%taun(I,J)/this%RATIO(I,J)
        ENDDO
    ENDDO
    !	 write(6,*) 'Yalin (0) or Engelund-Hansen (1):'
    !	 read(5,*) ityp
    write(6,*) 'itype', cco%TRANSEQTYPE, 'cco%sedsmoo',cco%sedsmoo,'din',cco%din
    !			ityp=0
    RHOFAC=1.65
    DO I=1,tns
        DO J=1,tnn
            IF(I.EQ.1) THEN
                this%dhds(I,J)=(rvo%hl(2,J)-rvo%hl(1,J))/(rvo%ds*(rvo%rn(I,J)))
            ELSE
                IF(I.EQ.tns) THEN
                    this%dhds(I,J)=(rvo%hl(tns,J)-rvo%hl(tns-1,J))/(rvo%ds*(rvo%RN(I,J)))
                ELSE
                    this%dhds(I,J)=(rvo%hl(I+1,J)-rvo%hl(I,J))/(rvo%ds*(rvo%RN(I,J)))
                ENDIF
            ENDIF
            IF(J.EQ.1) THEN
                !       DHDN(I,J)=(hl(I,2)-hl(I,1))/rvo%dn
                this%dhdn(I,J) = 0.
            ELSE
                IF(J.EQ.tnn) THEN
                    !        DHDN(I,J)=(hl(I,nn)-hl(I,nn-1))/rvo%dn
                    this%dhdn(I,J) = 0.
                ELSE
                    this%DHDN(I,J)=(rvo%hl(I,J+1)-rvo%hl(I,J))/rvo%dn
                ENDIF
            ENDIF
            ! now fix things for parts sticking out of channel
            !     do i=1,ns
            !        do j=2,nn-1

            IF(j.ne.1.and.j.ne.tnn) THEN
                if(rvo%ibc(i,j+1).eq.0.and.rvo%ibc(i,j-1).eq.0) then
                    this%dhdn(i,j) = 0.
                else
                    if(rvo%ibc(i,j-1)==0) then !I think this statement below is wrong.  Need to check.
                        !!                    dhdn(i,j) = abs(hl(i,j+1)-hl(i,j))/rvo%dn
                        !
                        this%dhdn(i,j) = -1*(rvo%eta(i,j+1)-rvo%eta(i,j))/rvo%dn
                    endif
                    if(rvo%ibc(i,j+1)==0) then
                        !                    dhdn(i,j) = (hl(i,j)-hl(i,j-1))/rvo%dn
                        this%dhdn(i,j) = -1*(rvo%eta(i,j)-rvo%eta(i,j-1))/rvo%dn
                    endif
                endif
            ENDIF
            !        enddo
            !     enddo
            DUM=((this%dhds(I,J)**2)+(this%dhdn(I,J)**2))**.5
            GAMG=ATAN(DUM)
            Taug = 0
            if(cco%CALCGRAVCORR == 1) then
                TAUG=1.0*TC*sin(GAMG)/sin(GAMC)
                TAN2= 0.
            else
                TAUG = 0.
                TAN2 = 0.
            endif

            ! rmcd hard wire taug=0 at j=1 and j=nn for flume walls
            if(j == 1 .or. j == tnn .or. rvo%ibc(i,j) == 0) then
                taug = 0.
            endif
            ! end PAN change

            IF(DUM.NE.0) THEN
                this%ustress(I,J)=this%ustress(I,J)+(TAUG*(this%dhds(I,J)/DUM))
                this%vstress(I,J)=this%vstress(I,J)+(TAUG*(this%dhdn(I,J)/DUM))
            ENDIF
        ENDDO
    ENDDO
    !612 CONTINUE

    !!!!!  YALIN EQUATION   !!!!!
    if(cco%TRANSEQTYPE == 0) then !Yalin
        DO  I=1,tns
            DO  J=1,tnn
                TC=(200.*(this%d(i,j)**2.))+(26.*this%d(i,j))+1.
                tdd = this%d(i,j)/100. !convert to MKS
                CALL tcrit2(tdd, kv, TC) !MKS
                TC = TC*10 !Convert to CGS
                A2=1.66*((TC/(1617.*this%d(i,j)))**.5)
                TB=((this%ustress(I,J)**2)+(this%vstress(I,J)**2))**.5
                IF (TB.LE.TC.or.rvo%ibc(i,j).eq.0) THEN
                    rvo%qs(I,J)=0.
                    rvo%qn(I,J)=0.
                ELSE
                    DTAU=(TB-TC)/TC
                    DUM=1.-(LOG(1.+(A2*DTAU)))/(A2*DTAU)
                    DUM=.635*this%d(i,j)*DUM*DTAU
                    rvo%qs(I,J)=(DUM*((TB)**.5))
                    rvo%qs(I,J)=rvo%qs(I,J)*this%ustress(I,J)/(TB)
                    TAN1=this%vstress(I,J)/this%ustress(I,J)
                    if(cco%calcgravcorr == 2) then
                        TAN2=cco%GRAVFLATBEDCORRCOEF*((TC/TB)**.5)*this%dhdn(I,J)
                        rvo%qn(I,J)=rvo%qs(I,J)*(TAN1+TAN2)

                    else
                        rvo%qn(I,J)=rvo%qs(I,J)*TAN1

                    ENDIF

                    !			write(6,*) i, j, tan1, tan2, rvo%qn(i,j)
                    IF(J.EQ.1.OR.J.EQ.tnn) THEN
                        rvo%qn(I,J)=0.
                    ENDIF
                ENDIF
            ENDDO
        ENDDO
    ELSEIF (cco%transeqtype == 1) then !Engelund and Hanssen
        do i=1,tns
            do j=1,tnn
                TC=(200.*(this%d(i,j)**2.))+(26.*this%d(i,j))+1.
                tdd = this%d(i,j)/100. !convert to MKS
                CALL tcrit2(tdd, kv, TC) !MKS
                TC = TC*10 !Convert to CGS
                TB=((this%ustress(I,J)**2)+(this%vstress(I,J)**2))**.5
                IF (TB.LE.TC.or.rvo%ibc(i,j).eq.0) THEN
                    rvo%qs(I,J)=0.
                    rvo%qn(I,J)=0.
                else
                    fac = ((1.65*980.0*this%d(i,j)**3.)**0.5)*rvo%cd(i,j)
                    rvo%qs(i,j)=fac*(tb/tc)**2.5
                    rvo%qs(i,j)=rvo%qs(i,j)*this%ustress(i,j)/TB
                    TAN1 = this%vstress(I,J)/this%ustress(I,J)
                    ! if(cco%calcgravcorr == 1 .AND. GRAVCORRTYPE == 2) then !pan - gravcorrtype doesn't exist anymore.
                    if(cco%calcgravcorr == 2) then
                        TAN2 = cco%gravflatbedcorrcoef*((tc/tb)**.5)*this%dhdn(i,j)
                        rvo%qn(I,J) = rvo%qs(I,J)*(TAN1+TAN2)
                    else
                        rvo%qn(I,J) = rvo%qs(I,J)*TAN1
                    endif
                    if(j.eq.1.or.j.eq.tnn) then
                        rvo%qn(i,j)=0.
                    endif
                endif
            enddo
        enddo
    ENDIF
    !	set upstream sed flux (can change later to allow recirculation)
    DO J = 1,tnn
        rvo%qs(cco%sedbcnode, J) = cco%BCFRACTION*rvo%qs(cco%sedbcnode, J)
        rvo%qn(cco%sedbcnode, J) = cco%BCFRACTION*rvo%qn(cco%sedbcnode, J)
    ENDDO

710 QTOTAL=0.
    DO I=1,tns
        this%qtot(I)=0.
        DO J=2,tnn
            this%qtot(I)=this%qtot(I)+(.5*rvo%dn)*(rvo%qs(I,J)+rvo%qs(I,J-1))
        ENDDO
        QTOTAL=QTOTAL+this%qtot(I)/tns
    ENDDO


    DO  I=1,tns
        DO  J=1,tnn
            SC=rvo%rn(I,J)
            IF(I.eq.1) THEN
                rvo%con(I,J)=0.
                !			        rvo%con(I,J)=(rvo%qs(2,J)-rvo%qs(1,J))/(rvo%ds*SC)
                !	              GO TO 643
            ELSEIF(I.EQ.tns) THEN
                !	              rvo%con(I,J)=(rvo%qs(I,J)-rvo%qs(I-1,J))/(rvo%ds*SC)
                !	              rvo%con(I,J)=(rvo%qs(2,J)-rvo%qs(1,J))/(rvo%ds*SC)
                rvo%con(I,J)=rvo%con(I-1,J)
                !			         rvo%con(I,J)=rvo%con(I,J)-(rvo%qn(I,J)/(R(I)*SC))
            ELSE
                rvo%con(I,J)=(rvo%qs(I+1,J)-rvo%qs(I-1,J))/(2.*rvo%ds*SC)
            ENDIF

            rvo%con(I,J)=rvo%con(I,J)-(rvo%qn(I,J)/(rvo%R(I)*SC))

            IF(J.EQ.1) THEN
                !	        rvo%con(I,1)=rvo%con(I,1)+(rvo%qn(I,2)-rvo%qn(I,1))/rvo%dn
                rvo%con(i,1)=0.
            ELSEIF(J.EQ.tnn) THEN
                !	         rvo%con(I,nn)=rvo%con(I,nn)+(rvo%qn(I,nn)-rvo%qn(I,nn-1))/rvo%dn
                rvo%con(i,tnn)=0.
            ELSE
                !first two if statements below worked well with parabolic bed
                IF(rvo%ibc(I,J-1).ne.-1) THEN
                    rvo%con(I,J)=rvo%con(I,J)+(rvo%qn(I,J)-rvo%qn(I,J-1))/(rvo%dn)
                    !	        ELSEIF(rvo%ibc(I,J+1).ne.-1)THEN
                    !	            rvo%con(I,J)=rvo%con(I,J)+(rvo%qn(I,J+1)-rvo%qn(I,J))/(rvo%dn)
                ELSE
                    rvo%con(I,J)=rvo%con(I,J)+(rvo%qn(I,J+1)-rvo%qn(I,J-1))/(2.*rvo%dn)
                ENDIF
            ENDIF
            if(i.lt.cco%sedbcnode) then
                rvo%con(i,j) = 0.
            endif

        ENDDO
    ENDDO


    !			do 652 j=1,nn
    !652			rvo%con(1,j)=rvo%con(2,j)
    DO I=1,tns
        DO J=1,tnn
            rvo%con(I,J)=1.54*rvo%con(I,J)
        ENDDO
    ENDDO

    !           write(6,*) 'rvo%con',(rvo%con(10,j),j=1,nn)
    !             write(6,*) 'rn',(rn(10,j),j=1,nn)
    IF (cco%sedsmoo.ge.1) then
        DO ISMOO=1,cco%sedsmoo
            DO I=2,tns-1
                DO J=1,tnn
                    if(rvo%ibc(i,j).eq.0) then
                        this%dum1(i,j)=0.
                    elseif(j.eq.1) then
                        this%dum1(I,1)=(rvo%con(I-1,1)+rvo%con(I,1)+rvo%con(I+1,1)+rvo%con(I,j+1))/4.
                        !             this%dum1(I,1)=(.5*rvo%con(I-1,1)+rvo%con(I,1)+.5*rvo%con(I+1,1))/2.
                    elseif(j.eq.tnn) then
                        this%dum1(I,tnn)=(rvo%con(I-1,tnn)+rvo%con(I,tnn)+rvo%con(I+1,tnn)+rvo%con(i,tnn-1))/4.
                        !             this%dum1(I,nn)=(.5*rvo%con(I-1,nn)+rvo%con(I,nn)+.5*rvo%con(I+1,nn))/2.
                    else
                        weight=1
                        if(rvo%con(i,j-1).ne.0) then
                            weight=weight+cco%sedsmoowght
                        endif
                        if(rvo%con(i,j+1).ne.0) then
                            weight=weight+cco%sedsmoowght
                        endif
                        if(rvo%con(i-1,j).ne.0) then
                            weight=weight+cco%sedsmoowght
                        endif
                        if(rvo%con(i+1,j).ne.0) then
                            weight=weight+cco%sedsmoowght
                        endif
                        this%dum1(I,J)=1.*rvo%con(I,J)+cco%sedsmoowght*rvo%con(I+1,J)+cco%sedsmoowght*rvo%con(I-1,J)
                        this%dum1(I,J)=this%dum1(I,J)+cco%sedsmoowght*rvo%con(I,J-1)+cco%sedsmoowght*rvo%con(I,J+1)
                        this%dum1(I,J)=this%dum1(I,J)/weight
                    endif
                ENDDO
            ENDDO
            !            write(6,*) 'smooth pass', ismoo, (this%dum1(10,j),j=1,nn)
            DO I=2,tns-1
                DO J=1,tnn
                    rvo%con(I,J)=this%dum1(I,J)
                ENDDO
            ENDDO

        ENDDO
    endif
    !            write(6,*) 'rvo%con',(rvo%con(10,j),j=1,nn)
    !            write(6,*) 'rn',(rn(10,j),j=1,nn)
    !            pause
    !	 DO 685 ISMOO=1,2
    !	 DO 680 I=2,ns-1
    !	 this%dum1(I,1)=(rvo%con(I-1,1)+rvo%con(I,1)+rvo%con(I+1,1)+rvo%con(I,j+1))/4.
    !680	 this%dum1(I,25)=(rvo%con(I-1,nn)+rvo%con(I,nn)+rvo%con(I+1,nn)+rvo%con(i,nn-1))/4.
    !	 DO 685 I=2,ns-1
    !	 rvo%con(I,1)=this%dum1(I,1)
    !685	  rvo%con(I,nn)=this%dum1(I,nn)
    DO  J=1,tnn
        !			rvo%con(ns,J)=(.2*rvo%con(ns-1,J)+rvo%con(ns,J))/1.2
        rvo%con(tns,j)= rvo%con(tns-1,j)
        rvo%con(1,J)=0.
    ENDDO

    !     FIND MINIMUM TIME STEP TO SATISFY CRITERIA FOR MAXIMUM ELEVATION CHANGE
    !     AS A FUNTION OF DEPTH !9/24/06 - rmcd
    mintdt = 1e12
    maxtdt = -1e12
    DO I = 1,tns
        DO J = 1, tnn
            IF(rvo%ibc(i,j).ne.0.and.rvo%con(i,j).ne.0) THEN
                tdt = abs((cco%tsfracdepth*rvo%hl(i,j))/rvo%con(i,j))
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

    DO I=1,tns
        DO J=1,tnn
            if(rvo%ibc(i,j).eq.0) then
                rvo%con(i,j)=0.
            endif
            if(dt*rvo%con(i,j) > rvo%hl(i,j)) then
                write(6,*)'Error in Sediment Transport', i, j, rvo%hl(i,j)/rvo%con(i,j)
            ENDIF

            if(i.lt.cco%sedbcnode) then
                rvo%con(i,j) = 0
            endif
            tmpdepth = rvo%hl(i,j)+(dt*rvo%con(i,j))
            !			if(tmpdepth < hmin) then
            !			    rvo%eta(i,j) = e(i,j)-hmin
            !			    rvo%hl(i,j) = hmin
            !			else
            rvo%eta(i,j)=rvo%eta(i,j)-(dt*rvo%con(i,j))
            rvo%hl(I,J)=rvo%hl(I,J)+(dt*rvo%con(I,J))
            !			endif
            !			rvo%eta(i,j)=rvo%eta(i,j)-(dt*rvo%con(i,j))
            !			rvo%hl(I,J)=rvo%hl(I,J)+(dt*rvo%con(I,J))
            if(rvo%hl(i,j).le.cco%hmin) then
                rvo%hl(i,j)=cco%hmin
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
!2000 if(nct.eq.nsteps) then
!        WRITE(7,2010) RUNID
!        WRITE(7,*) QS,QN,CON,hl,QTOT
!        write(7,*) ((i,j,con(i,j),j=1,nn),i=10,15)
!        write(10,*) hwt,erelax,urelax,arelax,itm,cd,q,mo
!        write(10,2010) runid
!        write(10,*) q,mo,hav
!        write(10,*) r,w,eta,ibc
!        write(10,*) xo,yo
!    endif
!2010 format(A20)
    return
    end subroutine csed2

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
    END MODULE CSedMod2
