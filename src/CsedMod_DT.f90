    MODULE CSedMod_DT_SUSP
    use RivVarMod
    use RivVertMod
    IMPLICIT NONE
    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:) :: QTOT
    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:) :: RATIO, d, DUM1
    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:) :: DHDN, DHDS, USTRESS, VSTRESS

    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:) :: hmax, lsub, lsubactdepth !DT
    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:) :: av ,bv, abh, bbh  !DT
    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:) :: imeta, dro, dref  !DT
    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:) :: Vsandsub, Vsandsurf   !DT
    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:) :: Vsandsubmax, Vsandsurfmax !DT
    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:) :: FRACS0, HFINE0 !DT

    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:,:,:) :: TCS
    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:,:) :: QSUSPN, QSUSPS, TOTCSZ
    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:) :: QSUSP, TOTCS
    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:) :: QSUSPXS, AVGCS

    !! For Suspended Sediment Calculation
    REAL(KIND=mp), DIMENSION(8) :: sedsize, sedsizefrac

    !! For Non-equilibrium sediment supply
    REAL(KIND=mp) :: NONEQSMOOTH

    CONTAINS
    SUBROUTINE alloc_csed3d_DT()
    INTEGER :: status
    ALLOCATE(tcs(ns, nn, nz,8), STAT = status)
    ALLOCATE(qsusps(ns, nn, nz), STAT = status)
    ALLOCATE(qsuspn(ns, nn, nz), STAT = status)
    ALLOCATE(totcsz(ns, nn, nz), STAT = status)
    ALLOCATE(qsuspxs(ns), STAT = status)
    ALLOCATE(avgcs(ns), STAT = status)
    ALLOCATE(qsusp(ns, nn), STAT = status)
    ALLOCATE(totcs(ns, nn), STAT = status)
    END SUBROUTINE

    SUBROUTINE dealloc_csed3d_DT()
    INTEGER :: status
    DEALLOCATE(tcs, STAT = status)
    DEALLOCATE(qsusps, STAT = status)
    DEALLOCATE(qsuspn, STAT = status)
    DEALLOCATE(totcsz, STAT = status)
    DEALLOCATE(qsuspxs, STAT = status)
    DEALLOCATE(avgcs, STAT = status)
    DEALLOCATE(qsusp, STAT = status)
    DEALLOCATE(totcs, STAT = status)
    END SUBROUTINE

    SUBROUTINE alloc_csed_DT()
    !		INTEGER, INTENT(IN) :: ns, nn, nz
    INTEGER :: status
    !Allocate REAL types
    ALLOCATE(QTOT(ns), STAT = status)
    ALLOCATE(RATIO(ns, nn), STAT = status)
    ALLOCATE(d(ns, nn), STAT = status)
    ALLOCATE(QS(ns, nn), STAT = status)
    !			ALLOCATE(QN(ns, nn), STAT = status)
    ALLOCATE(DUM1(ns, nn), STAT = status)
    ALLOCATE(DHDN(ns, nn), STAT = status)
    ALLOCATE(DHDS(ns, nn), STAT = status)
    ALLOCATE(USTRESS(ns, nn), STAT = status)
    ALLOCATE(VSTRESS(ns, nn), STAT = status)
    !			ALLOCATE(CON(ns, nn), STAT = status)
    !					!DT
    !			ALLOCATE(Fracs(ns, nn), STAT = status)
    ALLOCATE(hmax(ns, nn) , STAT = status)
    ALLOCATE(lsub(ns, nn) , STAT = status)
    ALLOCATE(lsubactdepth(ns, nn) , STAT = status)
    ALLOCATE(av(ns, nn) , STAT = status)
    ALLOCATE(bv(ns, nn) , STAT = status)
    ALLOCATE(abh(ns, nn) , STAT = status)
    ALLOCATE(bbh(ns, nn) , STAT = status)
    !			ALLOCATE(h(ns, nn) , STAT = status)
    ALLOCATE(imeta(ns, nn) , STAT = status)
    ALLOCATE(dro(ns, nn) , STAT = status)
    ALLOCATE(dref(ns,nn) , STAT = status)
    ALLOCATE(Vsandsub(ns, nn), STAT = status)
    ALLOCATE(Vsandsurf(ns, nn) , STAT = status)
    ALLOCATE(Vsandsubmax(ns, nn) , STAT = status)
    ALLOCATE(Vsandsurfmax(ns, nn) , STAT = status)
    !                       SUSPENDED SEDIMENT
    !			ALLOCATE(tcs(ns, nn, nz) , STAT = status)

    !NEW VARIABLE FOR THE UPSTREAM SAND INPUT
    ALLOCATE (fracs0(NN), STAT=status)
    ALLOCATE (Hfine0(NN), STAT=status)
    !END NEW VARIABLES

    !            ALLOCATE(qsint(ns), STAT=status)

    END SUBROUTINE alloc_csed_DT

    SUBROUTINE dealloc_csed_DT()
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
    !					!DT
    !			DEALLOCATE(Fracs, STAT = status)
    DEALLOCATE(hmax , STAT = status)
    DEALLOCATE(lsub , STAT = status)
    DEALLOCATE(lsubactdepth , STAT = status)
    DEALLOCATE(av , STAT = status)
    DEALLOCATE(bv , STAT = status)
    DEALLOCATE(abh , STAT = status)
    DEALLOCATE(bbh , STAT = status)
    !			DEALLOCATE(h , STAT = status)
    DEALLOCATE(imeta , STAT = status)
    DEALLOCATE(dro , STAT = status)
    DEALLOCATE(dref , STAT = status)
    DEALLOCATE(Vsandsub, STAT = status)
    DEALLOCATE(Vsandsurf , STAT = status)
    DEALLOCATE(Vsandsubmax , STAT = status)
    DEALLOCATE(Vsandsurfmax , STAT = status)
    !
    !            DEALLOCATE(tcs, STAT = status)

    !NEW VARIABLE FOR THE UPSTREAM SAND INPUT
    DEALLOCATE (fracs0, STAT=status)
    DEALLOCATE (Hfine0, STAT=status)
    !END NEW VARIABLES

    !            DEALLOCATE(qsint, STAT=status)

    END SUBROUTINE dealloc_csed_DT
    !


    SUBROUTINE csed_DT(DT,nct,nsteps,newdt)
    REAL(KIND=mp), INTENT(IN) :: DT
    INTEGER, INTENT(IN) :: nct, nsteps
    REAL(KIND=mp), INTENT(INOUT) :: newdt
    CHARACTER(LEN=20) RUNID
    REAL(KIND=mp) :: CDune, GSW, PI, GAMC, TB, TC, G, RHO
    REAL(KIND=mp) :: RAT, ZOSF, RHOFAC, DUM, GAMG, TAUG
    REAL(KIND=mp) :: A2, DTAU, TAN1, TAN2, QTOTAL, SC, weight
    REAL(KIND=mp) :: fac
    REAL(KIND=mp) :: PLINC
    INTEGER :: NM, I, J, K,ngc, ityp, ISMOO

    REAL(KIND=mp) :: mintdt, maxtdt, tdt
    REAL(KIND=mp) :: tmpdepth

    REAL(KIND=mp) :: t, csf, p, s, ws, kv
    REAL(KIND=mp) :: tdd

    !!             Daniele variables for sedmodel

    REAL(KIND=mp) :: bedporosity, he
    REAL(KIND=mp) :: taustar_rs0
    REAL(KIND=mp) :: Wdensity, Sdensity, deltaD
    REAL(KIND=mp) :: taustar_rs, tau_rs
    REAL(KIND=mp) :: phi_sand, Wprime_sand, Qsand
    REAL(KIND=mp) :: check, drough, ustar, vs

    REAL(KIND=mp), DIMENSION(8) :: cb

    REAL(KIND=mp) :: porosity, bedcoverage, gamma
    REAL(KIND=mp) ::  sstc, ss, ca, ustart, rexp, za
    REAL(KIND=mp) :: dz, tqsusp, csdepth, zcoord, h1, h2, tsusp
    REAL(KIND=mp) :: cs_washload, tdz, tdn
    REAL(kind = mp) :: vkc = 0.4
    INTEGER :: ierr
    bedporosity = 0.3

    runid='test     '
    !OPEN(7,FILE='seddat')
    !OPEN(10,FILE='top2')
    !OPEN(12,FILE='SedFlux')
    !			OPEN(11,FILE='topser')
    !			DIN=0.02

    !            count = 0
    CDune=.212
    !VKC=0.4
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
        RAT=1.+(CDune/(2.*VKC**2.))*(HD/WD)*((DLOG(HD/(ZOSF)))-1.)**2.
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
    !	 write(6,*) 'Yalin (0), Engelund-Hansen (1), Wilcock-Kenworthy just sand (2):'   !Daniele
    !	 read(5,*) ityp
    !            write(6,*) 'itype', TRANSEQTYPE, 'sedsmoo',SEDSMOO,'din',din
    !			ityp=0
    RHOFAC=1.65
    !			DO 612 I=1,ns
    !			DO 612 J=1,nn
    !			IF(I.EQ.1) THEN
    !			 DHDS(I,J)=(hl(2,J)-hl(1,J))/(ds*(rn(I,J)))
    !			ELSE IF(I.EQ.ns) THEN
    !			  DHDS(I,J)=(hl(ns,J)-hl(ns-1,J))/(ds*(RN(I,J)))
    !			ELSE
    !			  DHDS(I,J)=(hl(I+1,J)-hl(I,J))/(ds*(RN(I,J)))
    !			ENDIF
    !			IF(J.EQ.1) THEN
    !			 DHDN(I,J)=(hl(I,2)-hl(I,1))/dn
    !			ELSE IF(J.EQ.nn) THEN
    !			  DHDN(I,J)=(hl(I,nn)-hl(I,nn-1))/dn
    !			ELSE
    !			  DHDN(I,J)=(hl(I,J+1)-hl(I,J))/dn
    !!			  DHDN(I,J)=((hl(I,J+1)-hl(I,J)) + (hl(I,J)-hl(I,J-1)))/2*dn
    !			ENDIF
    !			DUM=((DHDS(I,J)**2)+(DHDN(I,J)**2))**.5
    !			GAMG=ATAN(DUM)
    !			if(CalcGravCorr) then
    !			    if(GRAVCORRTYPE == 0) then
    !	              TAUG=1.0*TC*sin(GAMG)/sin(GAMC)
    !	              TAN2= 0.
    !	          else
    !	              TAUG=0.
    !	          endif
    !	      else
    !	          TAUG = 0.
    !	          TAN2 = 0.
    !	      endif
    !
    !			IF(DUM.NE.0) THEN
    !			USTRESS(I,J)=USTRESS(I,J)+(TAUG*(DHDS(I,J)/DUM))
    !			VSTRESS(I,J)=VSTRESS(I,J)+(TAUG*(DHDN(I,J)/DUM))
    !			ENDIF
    !612			CONTINUE
    DO I=1,ns
        DO J=1,nn
            IF(I.EQ.1) THEN
                DHDS(I,J)=((-1.)*eta(ns,J)-eta(ns-1,J))/(ds*(RN(I,J)))
            ELSE IF(I.EQ.ns) THEN
                DHDS(I,J)=(-1.)*(eta(I+1,J)-eta(I,J))/(ds*(RN(I,J)))
            ELSE
                !			  DHDS(I,J)= (-1.)*(eta(I+1,J)-eta(I,J))/(ds*(RN(I,J)))
                DHDS(I,J)= (-1.)*(eta(I+1,J)-eta(I-1,J))/(ds*(RN(I,J)))
            ENDIF
            IF(J.EQ.1) THEN
                DHDN(I,J)=(-1.)*(eta(I,2)-eta(I,1))/dn
            ELSE IF(J.EQ.nn) THEN
                DHDN(I,J)=(-1.)*(eta(I,nn)-eta(I,nn-1))/dn
            ELSE
                !			  DHDN(I,J)=(-1.)*(eta(I,J+1)-eta(I,J))/dn
                DHDN(I,J)=((hl(I,J+1)-hl(I,J)) + (hl(I,J)-hl(I,J-1)))/2*dn
            ENDIF
            DUM=((DHDS(I,J)**2)+(DHDN(I,J)**2))**.5
            GAMG=ATAN(DUM)
            if(CalcGravCorr == 1) then
                if(GRAVCORRTYPE == 0) then
                    TAUG=1.0*TC*sin(GAMG)/sin(GAMC)
                    TAN2= 0.
                else
                    TAUG=0.
                    !			        TC=(200.*(D(I,J)**2.))+(26.*D(I,J))+1.
                    !	              TB=((USTRESS(I,J)**2)+(VSTRESS(I,J)**2))**.5
                    !	              TAN2=GRAVFLATBEDCORRCOEF*((TC/TB)**.5)*DHDN(I,J)
                    !	              write(6,*) i,j,tan2, tc, tb, dhdn(i,j)
                endif
            else
                TAUG = 0.
                TAN2 = 0.
            endif

            IF(DUM.NE.0) THEN
                USTRESS(I,J)=USTRESS(I,J)+(TAUG*(DHDS(I,J)/DUM))
                VSTRESS(I,J)=VSTRESS(I,J)+(TAUG*(DHDN(I,J)/DUM))

            ENDIF
        ENDDO
    ENDDO



777 IF (nct.EQ.0) THEN

        hmax = tmphmax
        abh = tmpabh
        bbh = tmpbbh
        lsub = tmplsub
        lsubactdepth = tmplsubactdepth

        DO  i=1,ns
            DO  j=1,nn



                !      INITIALIZATION HEIGHT OF SAND FROM FS AND FROM INPUT FILE

                CALL INITHDT( Fracs(I,J) , abh(i,j) , bbh(i,j) &
                    , Fsmin , hmax(I,J) , hfine(I,J) )




                !      CALCULATES Vsandsurfmax, AV, BV

                CALL BETAVOL( hmax(I,J) &
                    , abh(I,J), bbh(I,J), Vsandsurfmax(I,J) &
                    , AV(I,J) , BV(I,J),harea(I,J))   !changed on 10/16/09 for new betavol



                !      INITIALIZATION OF THE VARIABLES NEED FOR WILCKOCK AND KENTWORTHY
                !      SEDIMENT TRANSPORT ANALYSIS
                !	 CALCULATES Vsandsubmax, Vsandsub, Vsandsurf, dref


                CALL INITDT( Fracs(I,J) , Fsmin , hmax(I,J) , lsub(I,J) &
                    , bedporosity , Lsubactdepth(I,J)  &
                    , Vsandsurfmax(I,J), hfine(I,J) , harea(i,j) &
                    , AV(I,J) , BV(I,J), znaught(i,j), Z0Min &
                    , Z0Max , WK_RoughShapeParam , Vsandsubmax(I,J)  &
                    ,  Vsandsub(I,J) , Vsandsurf(I,J), dref(i,j))


                !      CALCULATES THE ELEVATION OF THE IMMOBILE GRAVEL BED FROM THE
                !      GIVEN HEIGHT OF SAND FOR WHICH THE HYDRAULIC MODEL WAS GENERATED

                CALL  INITIMETA( Fracs(i,j) ,harea(i,j),Fsmin ,hmax(i,j) &
                    , hfine(i,j), vsandsurf(i,j), vsandsurfmax(i,j) &
                    , he)

                imeta(i,j)=eta(i,j)-he



            ENDDO
        ENDDO

        !     SET THE INITIAL SAND FRACTION AND DEPTH AT THE UPSTREAM SEDIMENT BOUNDARY
        !     FOR WILCOCK MODEL FOR UPSTREAM FINE SEDIMENT INPUTS
        !     THIS ENSURSS THAT THE SAND IS ALWAYS AVAILABLE FROM UPSTREAM
        !     we probably need a switch here to enter the loop only if upstream
        !     input is checked

        DO J = 1,nn
            Fracs0(J) = Fracs(SEDBCNODE, J)
            Hfine0(J) = Hfine(SEDBCNODE, J)
        ENDDO

    ENDIF


    !	----WILCOCK AND KENWORTHY TRANSPORT MODEL WRR [2002]
    !
    !	only tracking movement of sand not gravel!!! Gravel always immobile
    !
    !	VARIABLES ARREYS
    !
    !	Fs( ns , nn )  = fraction of fine material (sand) variece with time
    !     Vsandsub( ns , nn ) = volume of sand in the subsurface layer
    !     Vsandsurf( ns , nn ) = volume of sand in the surface
    !
    !     VARIABLES SCALARS
    !	Dsand  = diamater of the fine sediment class
    !	Dg  = diameter of the coarse class
    !	taustar_rs0 = dimensionless shear stress for sand Fs=0
    !	taustar_rs  = dimensionless sand reference shear stress at Fs fraction of sand
    !     tau_rs      =  dimensional sand reference shear stress at Fs frac of sand
    !	phi_sand    = ratio applied reference shear stresses for sand
    !     Wprime_sand = dimensionless sand trasnport
    !     Qsand       = sand trasnport
    !	TB = applied shear stress as module of taus and taun

    !     PARAMETERS

    !      Can be modified by user
    !	taustar_rg0 = dimensionless shear stress for gravel Fs=0
    !	taustar_rg1 = dimensionless shear stress for gravel Fs=1
    !	taustar_rs1 = dimensionless shear stress for sand Fs=1
    !	alpha       = a parameter
    !	chi         = a parameter
    !	AK          = a parameter
    !
    !	Hard coded
    !	MATERIAL PROPERTIES
    !
    !	Wdensity = density of water
    !	Sdensity = density of sediment
    !	deltaD = submerged specific gravity of sediment (Sdensity-Wdensity)/Wdensity

    !
    !
    !
    !
    !	Sediment discharge is per unit width

    !     -------constant for sediment and water

    Wdensity=1.0         !density of water
    Sdensity=2.650      ! density of sediment
    deltaD=(Sdensity-Wdensity)/Wdensity !sediment submerged specific gravity

    !       ---Parameters in Wilcock & Kenworthy model

    !			taustar_rg0=0.035
    !			taustar_rg1=0.011
    !			taustar_rs1=0.065
    !			alpha=1


    !			AK=115
    !			chi=0.923
    !			phi_prime=1.27



    DO i=1,ns
        DO j=1,nn
            IF(I.lt.SEDBCNODE) THEN
                QS(I,J) = 0.
                QN(I,J)=0.

                !				No Qsand if there is not sand
            ELSEIF (ibc(i,j).EQ.0.OR.(Vsandsub(i,j)+Vsandsurf(i,j)).EQ.0.0) THEN
                QS(I,J)=0.
                QN(I,J)=0.
            ELSE






                TB=((USTRESS(I,J)**2)+(VSTRESS(I,J)**2))**.5 !shear stress module
                ustar=SQRT(TB/Wdensity)

                IF (USTRESS(I,J).eq.0) THEN !vector direction in radiant
                    TAN1=0.5*PI
                ELSE
                    TAN1=ATAN(VSTRESS(I,J)/USTRESS(I,J))
                ENDIF

                taustar_rs0=alpha*taustar_rg0*(Dg/Dsand)
                taustar_rs=taustar_rs1 &
                    + (taustar_rs0-taustar_rs1)*exp(-14*Fracs(i,j))
                tau_rs=taustar_rs*deltaD*Wdensity*981*Dsand
                phi_sand=TB/tau_rs


                IF (phi_sand.lt.phi_prime) THEN
                    Wprime_sand=0.002*(phi_sand)**(7.5)
                ELSE
                    Wprime_sand=AK &
                        * (1-chi/(phi_sand)**(0.25))**(4.5)
                ENDIF

                Qsand=Wprime_sand*Fracs(i,j) &
                    *(ustar)**(3.0)/(deltaD*981)


                !!!!!!!!!!!!!!!!!!!!!!gravel section
                !
                !				taustar_rg=taustar_rg1+(taustar_rg0-taustar_rg1)*exp(-14*Fracs(i,j))
                !				tau_rg=taustar_rg*deltaD*Wdensity*981*Dg
                !				phi_gravel=TB/tau_rg
                !				IF (phi_gravel.lt.phi_prime) THEN
                !						Wprime_gravel=0.002*(phi_gravel)**(7.5)
                !				ELSE
                !						Wprime_gravel=A*(1-chi/(phi_gravel)**(0.25))**(4.5)
                !				ENDIF
                !				Qgravel=Wprime_gravel(i,j)
                !    &				*(1-Fracs(i,j))*(ustar)**(3.0)*Sdensity/(deltaD*981)
                !				QTotal=Qsand+Qgravel
                !
                !!!!!!!!!!!!!!!!!!!!!!END

                if(dt == 0) then
                    tdt = 1.
                else
                    tdt = dt
                endif

                !                check=Qsand*1.54*tdt*dn    ! check if volume stransported is higher than available 10/09/09
                !  check on the natural coordinate system
                QN(i,j)=Qsand*sin(TAN1) !just for sand
                QS(i,j)=Qsand*cos(TAN1) !just for sand
                !			    TB=((USTRESS(I,J)**2)+(VSTRESS(I,J)**2))**.5
                !				qs(i,j)=Qsand * ustress(i,j)/TB !just for sand
                !				qn(i,j)=qs(i,j) * (vstress(i,j)/ustress(i,j)) !just for sand

                check=(abs(QS(i,j)/(rn(I,J)*ds) &
                    -(QN(I,J)/(R(I)*rn(I,J))))+abs(QN(I,J)/dn)) &
                    *1.54*harea(i,j)*tdt


                IF (check.GT.(Vsandsub(i,j)).AND.(Vsandsurf(i,j)).EQ.0.0)THEN

                    Qsand=Vsandsub(i,j)/(1.54*harea(i,j)*tdt) !(tdt*dn*1.54)!modified 10/09/09

                ELSEIF (check.GT.(Vsandsub(i,j)+Vsandsurf(i,J)) )then

                    Qsand=(Vsandsub(i,j)+Vsandsurf(i,j))/(1.54*harea(i,j)*tdt) !(tdt*dn*1.54)!modified 10/09/09
                    !      write(6,*)'high transport in cell', i,j,'Fine fraction',fracs(i,j)
                    !      pause


                ENDIF

                qn(i,j)=Qsand*sin(TAN1) !just for sand
                qs(i,j)=Qsand*cos(TAN1) !just for sand
                !			    TB=((USTRESS(I,J)**2)+(VSTRESS(I,J)**2))**.5
                !				qs(i,j)=Qsand * ustress(i,j)/TB !just for sand
                !				qn(i,j)=qs(i,j) * (vstress(i,j)/ustress(i,j)) !just for sand

                ! 25 June 2010-
                !!commented out 24 OCTOBER 2012 This is moved in the what is exiting the cell
                !           NO SAND TRANSPORT FROM A WET TO A DRY CELL
                !          IF (I.LT.NS) THEN
                !				IF (ibc(i+1,j).EQ.0.AND.QS(I,J).GT.0.0) THEN
                !					 QS(I,J)=0.0
                !					 GO TO 778
                !			    ENDIF
                !			ENDIF
                !			IF (I.gT.1) THEN
                !			 	IF (ibc(i-1,j).EQ.0.AND.QS(I,J).LT.0.0) THEN
                !				     QS(I,J)=0.0
                !					 GO TO 778
                !				ENDIF
                !			ENDIF
                !			IF (J.LT.NN) THEN
                !			  	IF (ibc(i,j+1).EQ.0.AND.QN(I,J).GT.0.0) THEN
                !					 QN(I,J)=0.0
                !					 GO TO 778
                !				ENDIF
                !		    ENDIF
                !			IF (J.GT.1) THEN
                !				IF (ibc(i,j-1).EQ.0.AND.QN(I,J).LT.0.0) THEN
                !				     QN(I,J)=0.0
                !					 GO TO 778
                !				ENDIF
                !			ENDIF
                ! 25 JUNE 2010 END
                !commented out 24 OCTOBER 2012 This is moved in the what is exiting the cell



                IF (j.eq.1.or.j.eq.nn) THEN
                    qn(i,j)=0.
                    qs(i,j)=0.
                ENDIF
            ENDIF
        ENDDO
    ENDDO


    !        END WILCOCK AND KENWORTHY TRANSPORT MODEL


    !         Rate of change of bed or div of Qs
    if(TRANSEQTYPE.eq.5) THEN !12/04/09 Exner equation.
        !                                For WK model the second order is not used but we used a mass balance on the fluxes
        !                                this is for accounting for the limited material in a cell

        DO I=1,ns
            DO J=1,nn
                ! moved from bottom of Do loop on 24 October 2012
                ! if dry nothing change
                IF(ibc(i,j).EQ.0) THEN
                    CON(i,j) = 0.0
                    cycle
                ENDIF
                ! end moved from bottom of Do loop on 24 October 2012

                SC=rn(I,J)
                ! what is exiting the cell is positive=erosion
                !           NO SAND TRANSPORT FROM A WET TO A DRY CELL  added 24 october 2012
                check=QS(I,J)/(ds*rn(I,J))-QN(I,J)/(R(I)*rn(I,J))
                IF (I.LT.NS) THEN
                    IF (ibc(i+1,j).EQ.0.AND.check.GT.0.0) THEN
                        check=0.0
                    ENDIF
                ENDIF
                IF (I.gT.1) THEN
                    IF (ibc(i-1,j).EQ.0.AND.check.LT.0.0) THEN
                        check=0.0
                    ENDIF
                ENDIF
                IF (J.LT.NN) THEN
                    IF (ibc(i,j+1).EQ.0.AND.QN(I,J).GT.0.0) THEN
                        QN(I,J)=0.0
                    ENDIF
                ENDIF
                IF (J.GT.1) THEN
                    IF (ibc(i,j-1).EQ.0.AND.QN(I,J).LT.0.0) THEN
                        QN(I,J)=0.0
                    ENDIF
                ENDIF

                con(i,j)=(abs(check)+abs(QN(I,J)/dn))*harea(i,j)

                IF (j.eq.1.or.j.eq.nn) THEN
                    con(i,j)=0.0
                ENDIF
                ! NEW 24 OCTOBER 2012 END


                !              con(i,j)=(abs(QS(I,J)/(ds*rn(I,J))-QN(I,J)/(R(I)*rn(I,J)))
                !     &                  +abs(QN(I,J)/dn))*harea(i,j)       modified 24 october 2012


                ! end what is exiting the cell
                IF(I.lt.SEDBCNODE) THEN

                    CON(I,J)=0.

                    !                   CON(I,J)=(QS(2,J)-QS(1,J))/(ds*SC)
                    !                    GO TO 643
                    cycle

                ELSEIF(I.EQ.ns) THEN
                    !                 CON(I,J)=(QS(I,J)-QS(I-1,J))/(ds*SC)
                    !                 CON(I,J)=(QS(2,J)-QS(1,J))/(ds*SC)
                    CON(I,J)=CON(I-1,J)
                    !GO TO 1643
                ELSE

                    ! sediment entering along s negative=deposition
                    check=QS(I-1,J)/(ds*rn(I-1,J))-QN(I-1,J)/(R(I-1)*rn(I-1,J))
                    if (check.GT.0.0) then
                        CON(I,J)=CON(I,J)-check*harea(i-1,j)
                    endif
                    check=QS(I+1,J)/(ds*rn(I+1,J))-QN(I+1,J)/(R(I+1)*rn(I+1,J))
                    if (check.LT.0.0) then
                        CON(I,J)=CON(I,J)+check*harea(i+1,j)
                    endif
                    ! end sediment entering along s


                ENDIF


                ! sediment entering along n
            IF(J.EQ.1) THEN
                    check=QN(I,2)/dn
                    if (check.LT.0.0) then
                        CON(I,1)=CON(I,1)+check*harea(I,2)
                    endif
                    CON(I,J)=CON(I,J)/(harea(i,j))  ! added 24 January 2012 correct for the lateral nodes
                    cycle
                ELSE
                    IF(J.EQ.nn) THEN
                        check=QN(I,nn-1)/dn
                        if (check.GT.0.0) then
                            CON(I,nn)=CON(I,nn)-check*harea(I,nn-1)
                        endif
                        CON(I,J)=CON(I,J)/(harea(i,j))  ! added 24 January 2012 correct for the lateral nodes
                        cycle
                    ENDIF
                    check=QN(I,J+1)/dn
                    if (check.LT.0.0) then
                        CON(I,J)=CON(I,J)+check*harea(I,J+1)
                    endif
                    check=QN(I,J-1)/dn
                    if (check.GT.0.0) then
                        CON(I,J)=CON(I,J)-check*harea(i,J-1)
                    endif

                    CON(I,J)=CON(I,J)/(harea(i,j))
                ENDIF
                IF(ibc(i,j).EQ.0) THEN
                    CON(i,j) = 0.0
                ENDIF
1645        enddo
1650    enddo

        !     End the Exner for WK

    else

        !         Rate of change of bed or div of Qs

        !For non-equilibrium case smooth change in sediment supply at boundary
710     IF(SedEqNode > SedBCNode) THEN
            DO I = SEDBCNODE,SEDEQNODE

                NONEQSMOOTH=(((1.0-BCFRACTION)/(SEDEQNODE-SEDBCNODE)) &
                    *(I-SEDBCNODE)) &
                    +BCFRACTION
                DO J = 1,nn
                    QS(I, J) = NONEQSMOOTH*QS(I, J)
                    QN(I, J) = NONEQSMOOTH*QN(I, J)
                ENDDO
            ENDDO
        ENDIF

        !            ENDDO
        !
        !
        !            CLOSE(11)

        DO I=1,ns
            DO J=1,nn
                SC=rn(I,J)
                IF(I.lt.SEDBCNODE) THEN

                    CON(I,J)=0.

                    !                   CON(I,J)=(QS(2,J)-QS(1,J))/(ds*SC)
                    !                    GO TO 643
                    cycle

                ELSEIF(I.EQ.ns) THEN
                    !                 CON(I,J)=(QS(I,J)-QS(I-1,J))/(ds*SC)
                    !                 CON(I,J)=(QS(2,J)-QS(1,J))/(ds*SC)
                    CON(I,J)=CON(I-1,J)
                    !GO TO 643
                ELSE
                    CON(I,J)=(QS(I+1,J)-QS(I-1,J))/(2.*ds*SC)
                ENDIF

643             CON(I,J)=CON(I,J)-(QN(I,J)/(R(I)*SC))

                IF(J.EQ.1) THEN
                    CON(I,1)=CON(I,1)+(QN(I,2)-QN(I,1))/dn
                    cycle
                ELSE

                    IF(J.EQ.nn) THEN
                        CON(I,nn)=CON(I,nn)+(QN(I,nn)-QN(I,nn-1))/dn
                        cycle
                    ENDIF

                    CON(I,J)=CON(I,J)+(QN(I,J+1)-QN(I,J-1))/(2.*dn)  !end equation 26 in Nelson&Smith AGU12

                ENDIF
645         enddo
650     enddo
        !            DO 630 I=1,ns
        !            DO 630 J=1,nn
        !               SC=rn(I,J)
        !               IF(I.lt.SEDBCNODE) THEN
        !
        !                 CON(I,J)=0.
        !
        !c                   CON(I,J)=(QS(2,J)-QS(1,J))/(ds*SC)
        !c                    GO TO 643
        !                   GO TO 625
        !
        !               ELSEIF(I.EQ.ns) THEN
        !C                 CON(I,J)=(QS(I,J)-QS(I-1,J))/(ds*SC)
        !c                 CON(I,J)=(QS(2,J)-QS(1,J))/(ds*SC)
        !                    CON(I,J)=CON(I-1,J)
        !                    GO TO 623
        !               ELSE
        !               IF(QS(i,j) < 0) THEN
        !                   CON(I,J)=(QS(I+1,J)-QS(I,J))/(1.*ds*SC)
        !               ELSE
        !                   CON(I,J)=(QS(I,J)-QS(I-1,J))/(1.*ds*SC)
        !               ENDIF
        !               ENDIF
        !
        !623                   CON(I,J)=CON(I,J)-(QN(I,J)/(R(I)*SC))
        !
        !               IF(J.EQ.1) THEN
        !                   CON(I,1)=CON(I,1)+(QN(I,2)-QN(I,1))/dn
        !                   GO TO 625
        !               ELSE
        !
        !                IF(J.EQ.nn) THEN
        !                    CON(I,nn)=CON(I,nn)+(QN(I,nn)-QN(I,nn-1))/dn
        !                    GO TO 625
        !                ENDIF
        !
        !                IF(Qn(i,j) < 0) THEN
        !               CON(I,J)=CON(I,J)+(QN(I,J+1)-QN(I,J))/(1.*dn)  !end equation 26 in Nelson&Smith AGU12
        !                ELSE
        !               CON(I,J)=CON(I,J)+(QN(I,J)-QN(I,J-1))/(1.*dn)  !end equation 26 in Nelson&Smith AGU12
        !                ENDIF
        !                ENDIF
        !625            CONTINUE
        !630            CONTINUE
    endif   !Daniele 12/02/09 endif on the Exner equation

    !
    QTOTAL=0.
    !
    QTOTAL=0.
    DO  I=1,ns
        QTOT(I)=0.
        DO  J=2,nn
            QTOT(I)=QTOT(I)+(.5*dn)*(QS(I,J)+QS(I,J-1))
        ENDDO
        !WRITE(12,*) I, QTOT(i)
        QTOTAL=QTOTAL+QTOT(I)/ns
    ENDDO
    !WRITE(12,*) 'AVERAGE', QTOTAL
    !CLOSE(12)

    !			do 652 j=1,nn
    !652			con(1,j)=con(2,j)
    DO I=1,ns
        DO J=1,nn
            CON(I,J)=1.54*CON(I,J)
        ENDDO
    ENDDO

    !         CON is div of Qs/porosity; 1/porosity=1.54;

    ! smoothing the DIV Qs over the grid


    IF (SEDSMOO.ge.1) then
        DO ISMOO=1,SEDSMOO
            DO I=2,ns-1
                DO J=1,nn
                    IF(ibc(i,j).ne.-1) THEN
                        con(i,j)=0.

                    ELSEIF(I.LT.SEDBCNODE) THEN
                        CON(i,j) = 0.
                    ELSE
                        if(j.eq.1) then
                            DUM1(I,1)=(CON(I-1,1)+CON(I,1)+CON(I+1,1)+con(I,j+1))/4.
                        elseif(j.eq.nn) then
                            DUM1(I,NN)=(CON(I-1,nn)+CON(I,nn)+CON(I+1,nn)+con(i,nn-1))/4.
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
                    ENDIF
660             ENDDO
            ENDDO
            write(6,*) 'smooth pass', ismoo, (dum1(10,j),j=1,nn)
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





    !      Daniele set of subroutine for updating Fs in the surface and subsurface
    !      for updating bed elevation and updating bed roughness



    DO I = 1,ns
        DO J = 1, nn

            IF (con(i,j).EQ.0.0) GOTO 779

            !			UPDATES THE SAND FRACTION ON THE SURFACE
            !
            !	Vsandsub(ns,nn) = Volume of sand in the subsurface active layer
            !     Vsandsurf(ns,nn)= Volume of sand in the surface
            !     Vsandsurfmax = max volume of sand in the surface !!it may be a vector
            !
            !     Vsandsubmax = max volume of sand in the subsurface active layer
            !     Fsmin = minimum fraction of sand
            !	VS = volume of sand moved each time step
            !	porosity = porosity of subsurface gravel matrix = 0.30
            !	Lsubactdepth = depth of the active subsurface layer
            !
            !	av = alpha parameter for the Beta model of Volume vs. Surface fraction of Sand
            !     bv = beta parameter for the Beta model of Volume vs. Surface fraction of Sand




            VS = CON(I,J)*dt*harea(i,j)
            CALL VOLFRAC(av(i,j),bv(i,j), VS, Vsandsub(i,j) &
                , Vsandsurf(i,j),Vsandsurfmax(i,j), Vsandsubmax(I,J) &
                , Fsmin, Fracs(i,j))

            !             UPDATES THE ELEVATION OF SAND OVER THE GRAVEL
            !             hfine(ns,nn) is the relative change in sand depth from lowest point into the gravel

            !
            !	abh = alpha parameter for the Beta model of Sand height vs. Surface fraction of Sand
            !	bbh = beta parameter for the Beta model of Sand height vs. Surface fraction of Sand
            !	hmax = maximum height of the gravel surface  !!
            !
            !


            CALL ELEVDT(abh(i,j), bbh(i,j), hmax(i,j) &
                , Vsandsurf(I,J), Vsandsurfmax(i,j) &
                , Fsmin, Fracs(I,J), harea(i,j), hfine(I,J), he)

            if(he.GT.hl(i,j).and.ibc(i,j).eq.-1) then
                write(6,*)'Error in Sediment Transport', &
                    i, j, VS/(con(i,j)*harea(i,j))
            ENDIF

            tmpdepth = hl(i,j)+he
            !			if(ibc(i,j).ne.0) then
            !			    if(tmpdepth < hmin) then
            !			        eta(i,j) = e(i,j)-hmin
            !			        hl(i,j) = hmin
            !			    else
            !			        eta(i,j)=imeta(i,j)+he
            !			        hl(I,J)=hl(I,J)-he
            !			    endif
            !			endif
            if(ibc(i,j).eq.-1) then
                if(he.ge.hl(i,j)) then
                    eta(i,j) = imeta(i,j)+(0.5*hl(i,j))
                    hl(i,j) = hl(i,j)-(0.5*hl(i,j))
                else
                    eta(i,j)=imeta(i,j)+he
                    hl(I,J)=hl(I,J)-he
                endif
            endif

            !
            !		UPDATES BED ROUGHNESS COEFFICIENT
            !
            !	Dro(ns,nn) = roughness height for each cell
            !
            !	Dref = reference diameter for hydraulic roughness
            !	Drough = roughness height dimension as a function of Dref and anout of sand
            !
            !
            !
            !
            IF(WK_RoughType == 1) THEN

                CALL roughdt( Fracs(I,J) , Fsmin , Z0Min &
                    , Dref(i,j) , hmax(i,j) &
                    , hfine(i,j) , WK_RoughShapeParam, Drough )

                znaught(i,j) = Drough
            ENDIF

779     ENDDO
    ENDDO

    !     SET CONSTANT SAND FRACTION AND DEPTH FOR WILCOCK MODEL
    !     THIS ENSURSS THAT THE SAND IS ALWAYS AVAILABLE FROM UPSTREAM
    !     we probably need a switch here to enter the loop only if upstream
    !     input is checked

    IF(WKConstBndry) THEN
        DO J = 1,nn
            Fracs(SEDBCNODE, J) = Fracs0(J)
            HFINE(SEDBCNODE, J) = HFINE0(J)
        ENDDO
    ENDIF
    if(WK_RoughType == 1) then
        CALL ZOTOCDTWO()
    endif

!    if (TRANSEQTYPE.eq.2) go to 2000         !Daniele
!
!2000 if(nct.eq.nsteps) then
!        !WRITE(7,2010) RUNID
!        !WRITE(7,*) QS,QN,CON,hl,QTOT
!        !write(7,*) ((i,j,con(i,j),j=1,nn),i=10,15)
!        !			    write(10,*) hwt,erelax,urelax,arelax,itm,cd,q,mo
!        !			    write(10,2010) runid
!        !			    write(10,*) q,mo,hav
!        !			    write(10,*) r,w,eta,ibc
!        !			    write(10,*) xo,yo
!    endif
2010 format(A20)
    return
    end subroutine csed_DT

    REAL(KIND=8) FUNCTION falveld(d,nu,s,g,p,csf)
    IMPLICIT NONE
    ! Computes particle fall velocity according to W. Dietrich, WRR 18(6),
    ! 1615-1626, 1982. falveld in m/s.
    ! Input parameters:
    !   d    particle diameter (m)
    !   nu   kinematic viscosity (m^2/s)
    !   s    specific gravity
    !   g    acceleration due to gravity (m/s^2)
    !   p    Powers roundness factor
    !   csf  Corey shape factor (must be >= 0.2)
    !
    ! NOTE: for a sphere, p=6.0, csf=1.0; for natural sediments csf=0.7.
    ! Function by Francisco Simoes, March 1999.

    REAL(KIND=mp), INTENT(IN) :: d,nu,s,g,p,csf
    REAL(KIND=mp) :: dstar,logdstar,r1,r2,r3,tanhdstar,wstar
    REAL(KIND=mp), PARAMETER :: one=1.0

    dstar = (s-one)*g*d**3/(nu*nu)
    logdstar = LOG10(dstar)
    tanhdstar = DTANH(logdstar-4.6D0)
    r1 = -3.76715D0+1.92944*logdstar-9.815D-2*logdstar*logdstar &
        -5.75D-3*logdstar**3+5.6D-4*logdstar**4
    r2 = LOG10(one-(one-csf)/8.5D-1)
    r2 = r2-tanhdstar*(one-csf)**2.3
    r2 = r2+0.3D0*(0.5D0-csf)*(one-csf)*(one-csf)*(logdstar-4.6D0)
    r3 = (6.5D-1-csf*tanhdstar/2.83D0)**(one+(3.5D0-p)/2.5D0)
    wstar = r3*10.0D0**(r1+r2)
    falveld = (wstar*(s-one)*g*nu)**(one/3.0D0)

    END FUNCTION falveld

    SUBROUTINE TCRIT(D, KV, TC)
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
