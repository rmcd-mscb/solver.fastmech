      MODULE CSedMod_DT_SUSP
	use RivVarMod
	use RivVertMod
	IMPLICIT NONE
	REAL, ALLOCATABLE, DIMENSION(:) :: QTOT
	REAL, ALLOCATABLE, DIMENSION(:,:) :: RATIO, d, DUM1
	REAL, ALLOCATABLE, DIMENSION(:,:) :: DHDN, DHDS, USTRESS, VSTRESS
	
	REAL, ALLOCATABLE, DIMENSION(:,:) :: hmax, lsub, lsubactdepth !DT
	REAL, ALLOCATABLE, DIMENSION(:,:) :: av ,bv, abh, bbh			  !DT
	REAL, ALLOCATABLE, DIMENSION(:,:) :: imeta, dro, dref	  !DT
	REAL, ALLOCATABLE, DIMENSION(:,:) :: Vsandsub, Vsandsurf		  !DT
	REAL, ALLOCATABLE, DIMENSION(:,:) :: Vsandsubmax, Vsandsurfmax	  !DT
	REAL, ALLOCATABLE, DIMENSION(:) :: FRACS0, HFINE0 !DT
	
	REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: TCS
	REAL, ALLOCATABLE, DIMENSION(:,:,:) :: QSUSPN, QSUSPS, TOTCSZ
	REAL, ALLOCATABLE, DIMENSION(:,:) :: QSUSP, TOTCS
	REAL, ALLOCATABLE, DIMENSION(:) :: QSUSPXS, AVGCS
	
!	REAL, ALLOCATABLE, DIMENSION(:) :: qsint

	
!	REAL :: AK, ALPHA, CHI, PHI_PRIME
!	REAL :: TAUSTAR_RG0, TAUSTAR_RG1, TAUSTAR_RS1
!	REAL :: Dsand, Dg, Fsmin
!	REAL :: Z0Min, Z0Max
!	LOGICAL :: WKConstBndry
!    
!      REAL :: WK_RoughShapeParam
!      INTEGER :: WK_RoughType
      !! For Suspended Sediment Calculation
      REAL, DIMENSION(8) :: sedsize, sedsizefrac
      
      !! For Non-equilibrium sediment supply
      REAL :: NONEQSMOOTH

		!! For DT Wilcock-Kenworthy Transport !! rmcd 5/2/07
!		REAL :: tmphmax, tmpabh, tmpbbh, tmplsub, tmplsubactdepth

!		!!For Rouse susp load Nelson
!		REAL :: cs,ca,gamma,za
!		
!		REAL :: cb, porosity, bedcoverage, gamma
!		REAL :: ws, sstc, s, ca, ustart, rexp, za
		!!Need to add z values, uz and vz back in for conc and load calc



!	REAL :: rsin, rcos, ux, uy, uxnew, uynew
!	REAL, ALLOCATABLE, DIMENSION(:,:) :: CON
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
			REAL, INTENT(IN) :: DT
			INTEGER, INTENT(IN) :: nct, nsteps
			REAL, INTENT(INOUT) :: newdt
!		parameter(ns=41,nn=25,nz=11)
	!	   REAL mo,hav
	!	   common taus(ns,nn),taun(ns,nn),hl(ns,nn),eta(ns,nn),rn(ns,nn),
	!	  & ibc(ns,nn),r(ns),w(ns),cd,ds,dn,mo,q,hwt,urelax,erelax,
	!	  & arelax,itm,hav,iplinc,xo(ns),yo(ns),uz(ns,nn,nz),vz(ns,nn,nz)
	!	dimension RATIO(ns,nn),d(ns,nn),
	!	 &	  QS(ns,nn),QN(ns,nn),DUM1(ns,nn),
	!	 &	  DHDN(ns,nn),DHDS(ns,nn),USTRESS(ns,nn),
	!	 &	  VSTRESS(ns,nn),QTOT(ns),CON(ns,nn)

			CHARACTER*20 RUNID
			REAL :: CDune, VKC, GSW, PI, GAMC, TB, TC, G, RHO
			REAL :: RAT, ZOSF, RHOFAC, DUM, GAMG, TAUG
			REAL :: A2, DTAU, TAN1, TAN2, QTOTAL, SC, weight
			REAL :: fac
			REAL :: PLINC
			INTEGER :: NM, I, J, K,ngc, ityp, ISMOO
			
	      REAL :: mintdt, maxtdt, tdt
	      REAL :: tmpdepth
	      
		REAL(KIND =8) :: t, csf, p, s, ws, kv
		REAL(KIND=8) :: tdd

!!             Daniele variables for sedmodel
			
			REAL :: bedporosity, he
			REAL :: taustar_rs0 
			REAL :: Wdensity, Sdensity, deltaD 
			REAL :: taustar_rs, tau_rs         
			REAL :: phi_sand, Wprime_sand, Qsand
			REAL :: check, drough, ustar, vs

		!!For Rouse susp load Nelson
!		REAL :: cs,ca,gamma,za
!		
        REAL, DIMENSION(8) :: cb
        
		REAL :: porosity, bedcoverage, gamma
		REAL ::  sstc, ss, ca, ustart, rexp, za
		REAL :: dz, tqsusp, csdepth, zcoord, h1, h2, tsusp
        REAL :: cs_washload, tdz, tdn

        INTEGER :: ierr
        
!        INTEGER :: count
!        REAL :: tqsint
        
!			REAL :: tau_rg, taustar_rg  !variable for gravel class
!			REAL :: phi_gravel, Wprime_gravel, Qgravel, QTotal  !variables for gravel class
!
!			End daniele variables for WK
!            CANCELLALE DOPO 

!			INTEGER CALCGRAVCORR
!			REAL qn(10,10), QS(10,10)
!			real con(10,10), eta(10,10), hl(10,10)
            bedporosity = 0.3
	
!			CALL alloc_csed()
			runid='test 			  '
			OPEN(7,FILE='seddat')
			OPEN(10,FILE='top2')
			OPEN(12,FILE='SedFlux')
!			OPEN(11,FILE='topser')
!			DIN=0.02
            
!            count = 0
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
c	 WRITE(6,*) 'ENTER THE DUNE HEIGHT IN CM (0 FOR PLANE BED):'
c	 READ(5,*) HD
!			hd=0.
			IF(HD.LE.0) THEN
			 RAT=1.
			 GO TO 600
			ENDIF
!			WRITE(6,*) 'ENTER THE AVERAGE DUNE WAVELENGTH IN CM:'
!			READ(5,*) WD
			IF(WD.LE.0) THEN
			 WD=1.
			ENDIF
			ZOSF=.2*DIN
			RAT=1.+(CDune/(2.*VKC**2.))*(HD/WD)*((ALOG(HD/(ZOSF)))-1.)**2.
600			DO 610 I=1,ns
			DO 610 J=1,nn
			D(I,J)=DIN
c	 hd=0.2*h(i,j)
c	 wd=20*hd
c	 RAT=1.+(CDune/(2.*VKC**2.))*(HD/WD)*((ALOG(HD/(ZOSF)))-1.)**2.
			RATIO(I,J)=RAT
			USTRESS(I,J)=taus(I,J)/RATIO(I,J)
610			VSTRESS(I,J)=taun(I,J)/RATIO(I,J)
c	 write(6,*) 'Yalin (0), Engelund-Hansen (1), Wilcock-Kenworthy just sand (2):'   !Daniele
c	 read(5,*) ityp
!            write(6,*) 'itype', TRANSEQTYPE, 'sedsmoo',SEDSMOO,'din',din
!			ityp=0
			RHOFAC=1.65
			DO 612 I=1,ns
			DO 612 J=1,nn
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
			 DHDN(I,J)=(hl(I,2)-hl(I,1))/dn
			ELSE
			 IF(J.EQ.nn) THEN
			  DHDN(I,J)=(hl(I,nn)-hl(I,nn-1))/dn
			 ELSE
			  DHDN(I,J)=(hl(I,J+1)-hl(I,J))/dn
			 ENDIF
			ENDIF
			DUM=((DHDS(I,J)**2)+(DHDN(I,J)**2))**.5
			GAMG=ATAN(DUM)
			if(CalcGravCorr) then
			    if(GRAVCORRTYPE == 0) then
	              TAUG=1.0*TC*SIN(GAMG)/SIN(GAMC)
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
612			CONTINUE
			
		if (TRANSEQTYPE.eq.1) THEN
		    go to 700
		else if(TRANSEQTYPE.eq.2) THEN
		    go to 777         !Daniele 
        ENDIF
!		!!!!  YALIN EQUATION   !!!!!

			DO 615 I=1,ns
			DO 615 J=1,nn
			
			TC=(200.*(D(I,J)**2.))+(26.*D(I,J))+1.
			tdd = D(I,J)/100. !convert to MKS
			CALL tcrit(tdd, kv, TC) !MKS 
			TC = TC*10 !Convert to CGS
			A2=1.66*((TC/(1617.*D(I,J)))**.5)
			TB=((USTRESS(I,J)**2)+(VSTRESS(I,J)**2))**.5
			IF (TB.LE.TC.or.ibc(i,j).eq.0) THEN
			 QS(I,J)=0.
			 QN(I,J)=0.
			 GO TO 615
			ENDIF
			DTAU=(TB-TC)/TC
			DUM=1.-(ALOG(1.+(A2*DTAU)))/(A2*DTAU)
			DUM=.635*D(I,J)*DUM*DTAU
		    QS(I,J)=(DUM*((TB)**.5))
			QS(I,J)=QS(I,J)*USTRESS(I,J)/(TB)
			TAN1=VSTRESS(I,J)/USTRESS(I,J)
			if(CalcGravCorr .AND. GRAVCORRTYPE == 1) then
				TAN2=GRAVFLATBEDCORRCOEF*((TC/TB)**.5)*DHDN(I,J)
				QN(I,J)=QS(I,J)*(TAN1+TAN2)
		    else
                QN(I,J)=QS(I,J)*TAN1
		    ENDIF
            QS(I,J) = QS(I,J) * SEDEQMULT
            QN(I,J) = QN(I,J) * SEDEQMULT
!			TAN2=GSW*1.50*((TC/TB)**.5)*DHDN(I,J)
!       	TAN2=0.
!			write(6,*) i, j, tan1, tan2, qn(i,j)
			IF(J.EQ.1.OR.J.EQ.nn) THEN
			 QN(I,J)=0.
			ENDIF
615		 CONTINUE
			go to 710

700		do 705 i=1,ns
			do 705 j=1,nn 



			!!!!! Engelund and Hansen !!!!
			
			TC=(200.*(D(I,J)**2.))+(26.*D(I,J))+1.
			tdd = D(I,J)/100. !convert to MKS
			CALL tcrit(tdd, kv, TC) !MKS 
			TC = TC*10 !Convert to CGS
			TB=((USTRESS(I,J)**2)+(VSTRESS(I,J)**2))**.5
			IF (TB.LE.TC.or.ibc(i,j).eq.0) THEN
			 QS(I,J)=0.
			 QN(I,J)=0.
			 GO TO 705
			ENDIF
			
			fac = ((1.65*980.0*d(i,j)**3.)**0.5)*cd(i,j)
			qs(i,j)=fac*(tb/tc)**2.5
			qs(i,j)=qs(i,j)*ustress(i,j)/TB
			TAN1 = VSTRESS(I,J)/USTRESS(I,J)
			if(CalcGravCorr .AND. GRAVCORRTYPE == 1) then
			    TAN2 = gravflatbedcorrcoef*((tc/tb)**.5)*dhdn(i,j)
			    QN(I,J) = QS(I,J)*(TAN1+TAN2)
			else
			    QN(I,J) = QS(I,J)*TAN1
			endif
						
            QS(I,J) = QS(I,J) * SEDEQMULT
            QN(I,J) = QN(I,J) * SEDEQMULT
			
			if(j.eq.1.or.j.eq.nn) then
			 qn(i,j)=0.
			endif
			
705			continue

!! Suspended load routines (probably need to add switch for Rouse vs ad/diff)
!! Code uses reduced stress for ref conc and total stress for ev (Rouse exp)
!!      ! mean grainsize fraction for upper monument creek
!!      sedsize = (/.063, .125, .25, .5, 1., 2., 4., 8./)  !In millimeters
!!      sedsizefrac=(/.0159,.0359,.0815,.1558,.2257,.2317,.1374,.1084/) 
!!      !mean greainsize of coarsest fraction upper monument creek
!!       sedsize = (/.063, .125, .25, .5, 1., 2., 4., 8./)  !In millimeters
!!      sedsizefrac=(/.00149,.00746,.0593,.1492,.2502,.3088,.1443,.0772/) 
!!      !mean greainsize of finest fraction upper monument creek
!!       sedsize = (/.063, .125, .25, .5, 1., 2., 4., 8./)  !In millimeters
!!      sedsizefrac=(/.2414,.4139,.1889,.05157,.0131,.0047,.0023,.0/) 
!      !median greainsize without finest fraction bartop sample upper monument creek
!       sedsize = (/.063, .125, .25, .5, 1., 2., 4., 8./)  !In millimeters
!      sedsizefrac=(/.00226,.0133,.0645,.1561,.2277,.2428,.1491,.1059/) 
!
!      !These three should be set in a dialogue
!	cb=1.-bedporosity  !Set by size fraction below
!	bedcoverage=1.
!	gamma=0.0027
!	t = 20.   !Temperature    
!	csf = 0.7 !Cory Shape Factor
!	p = 6.    !Powers Roundness Factor
!	s = 2.65  !Specific Gravidy
!	kv = 1.79D-6/(1.0D0+3.37D-2*t+2.21D-5*t*t) !MKS  !Kinematic Viscosity
!	dg = 9.81 !Gravity MKS units
!!	cs_washload = 2.8e-4 !unitless Here based on median silt/clay fraction
!!	                     !of both upper and lower monument creek samples
!	cs_washload = 0. !unitless Here based on average silt/clay fraction
!	                     !of both upper and lower monument creek samples
!!	call wset(d,ws)  !need to add correct call for settling vel
!!	call tcrit(d,tc)  !need to add correct call of tau critical
!!      cs = 0.
!
!      !initialize arrays
!      tcs = 0.0
!      qsusps = 0.0
!      qsuspn = 0.0
!      qsusp = 0.0   !streamwise flux at each node
!      totcs = 0.0
!      totcsz = 0.0
!      qsuspxs = 0.0 !Total flux at each cross-section
!      avgcs = 0.0 !Average concentration at each cross-section
!      !
!	do i=1,ns
!	    do j=1,nn
!	      do k = 1,nz
!	        do ngc = 1,8
!	          if(ibc(i,j).ne.0) THEN
!	          cb(ngc) = (1-bedporosity)*sedsizefrac(ngc)
!	          tdd = sedsize(ngc)/1000. !convert mm to m
!	          ws = falveld(tdd,kv,s,dg,p,csf) !MKS
!	          ws = ws*100 !CGS
!	          CALL tcrit(tdd, kv, TC) !MKS
!	          TC = TC*10 !CGS
!	          ss=(((ustress(i,j)**2.+vstress(i,j)**2.)**.5)/tc)-1  !assumes dune correction has been made above
!	          ca=bedcoverage*(cb(ngc)*gamma*ss)/(1+gamma*ss)
! 	          ustart=(((taus(i,j)**2.+taun(i,j)**2.)**.5)/rho)**.5
! 	          rexp=ws/(0.4*ustart)
!	          za=znaught(i,j)
!    	          IF(ss.gt.0) THEN
!!	            csdepth = zz(i,j,nz)-eta(i,j)
!	            zcoord = zz(i,j,k) - eta(i,j)
!	            h1 = hl(i,j)-zcoord
!	            h2 = hl(i,j)-za
!	            if(h2.lt.0.or.h1.lt.0) then
!	              tcs(i,j,k,ngc) = 0.0
!	            else
!                    tcs(i,j,k,ngc)=ca*((za/zcoord)*(h1/h2))**rexp
!                  endif
!!	            write(*,*) za/zcoord, h1/h2, rexp, tcs(i,j,k,ngc)
!!    	            tcs(i,j,k,ngc)=ca*(za/(zz(i,j,k)-eta(i,j)))**rexp
!                ELSE
!                  tcs(i,j,k,ngc) = 0.0
!                ENDIF
!	          ELSE 
!	            tcs(i,j,k,ngc) = 0.0
!	          ENDIF
!	          ENDDO !END NGC LOOP
!	          DO NGC = 1,8
!        	      qsusps(i,j,k)=qsusps(i,j,k)+(tcs(i,j,k,ngc))
!!        	      qsuspn(i,j,k)=qsuspn(i,j,k)+tcs(i,j,k,ngc)
!	          ENDDO
!	          qsusps(i,j,k) = qsusps(i,j,k)+cs_washload
!	          totcsz(i,j,k) = qsusps(i,j,k)!*(uz(i,j,k)*u(i,j)) !to simulate velocity weighted concentrations of measuremnts
!        	    qsusps(i,j,k)=qsusps(i,j,k)*uz(i,j,k)
!!        	    qsuspn(i,j,k)=qsuspn(i,j,k)*vz(i,j,k)
!	             
!	       ENDDO ! END K LOOP
!   	    ENDDO ! END J LOOP
!	ENDDO ! END K LOOP
!   	!Calc total load in vertical then integrate across each row
!   	DO i = 1,ns
!       DO j = 1,nn
!        tqsusp = 0.
!        tsusp = 0.
!        tdz = 0.
!        IF(ibc(i,j).NE.0) THEN
!            DO k = 2,nz
!            dz = zz(i,j,k)-zz(i,j,k-1)
!            tdz = tdz+dz
!!            write(*,*) i,j,k,dz
!            tqsusp = tqsusp +((qsusps(i,j,k)+qsusps(i,j,k-1))*dz/2.0)
!            tsusp = tsusp+((totcsz(i,j,k)+totcsz(i,j,k-1))*dz/2.0)
!            ENDDO ! END K LOOP - VERTICAL INTEGRATION
!            qsusp(i,j) = tqsusp
!            if(tdz.eq.0.0) then
!                totcs(i,j) = 0.0
!            else
!                totcs(i,j) = tsusp/tdz
!            endif
!            ELSE
!         qsusp(i,j) = 0.0
!         totcs(i,j) = 0.0
!        ENDIF
!       ENDDO
!       
!        !INTEGRATE ACROSS EACH ROW OF NODES
!        tdn = 0.
!       DO j = 2,nn
!       IF(ibc(i,j).ne.0) THEN
!        tdn = tdn+dn
!        qsuspxs(i) = qsuspxs(i)+((qsusp(i,j)+qsusp(i,j-1))*dn/2.0)
!        avgcs(i) = avgcs(i) + ((totcs(i,j)+totcs(i,j-1))*dn/2.0)
!       ENDIF
!       ENDDO
!       avgcs(i) = avgcs(i)/tdn;
!      ENDDO
!	OPEN(23,file = "XSSuspLoad.dat", status = 'unknown', iostat = ierr)
!	WRITE(23,*)'VARIABLES = "XS","Suspended Sediment Load"'
!	WRITE(23,*)'ZONE I=', ns, ',DATAPACKING=POINT'
!	  DO i = 1, ns
!	          WRITE(23,*) i, qsuspxs(i)/q, avgcs(i)
!        ENDDO
!      CLOSE(23)
!	OPEN(22,file = "SuspLoad.dat", status = 'unknown')
!	WRITE(22,*)'VARIABLES = "X", "Y", "Suspended Sediment Concentration", 
!     &                         "Suspended Sediment Load"'
!	WRITE(22,*)'ZONE I=', ns, ',J=', nn, ',DATAPACKING=POINT'
!	  DO i = 1, ns
!	      DO j = 1, nn
!	          WRITE(22,*) x(i,j), y(i,j), totcs(i,j), qsusp(i,j)
!	      ENDDO
!        ENDDO
!      CLOSE(22)
!! END OF SUSPENDED SEDIMENT CALCULATION
            go to 710

!!	INITIALIZE BED TOPOGRAPHY FOR SEDIMENT TRASNPORT
!!    DANIELE TONINA      



777	IF (nct.EQ.0) THEN
			
        hmax = tmphmax
        abh = tmpabh
        bbh = tmpbbh
        lsub = tmplsub
        lsubactdepth = tmplsubactdepth

		DO  i=1,ns
			DO  j=1,nn
			
			

!      INITIALIZATION HEIGHT OF SAND FROM FS AND FROM INPUT FILE

				CALL INITHDT( Fracs(I,J) , abh(i,j) , bbh(i,j)
     &                        , Fsmin , hmax(I,J) , hfine(I,J) )   




!      CALCULATES Vsandsurfmax, AV, BV
	
                CALL BETAVOL( hmax(I,J) 
     &                 , abh(I,J), bbh(I,J), Vsandsurfmax(I,J)
     &                 , AV(I,J) , BV(I,J),harea(I,J))   !changed on 10/16/09 for new betavol



!      INITIALIZATION OF THE VARIABLES NEED FOR WILCKOCK AND KENTWORTHY 
!      SEDIMENT TRANSPORT ANALYSIS
!	 CALCULATES Vsandsubmax, Vsandsub, Vsandsurf, dref


				CALL INITDT( Fracs(I,J) , Fsmin , hmax(I,J) , lsub(I,J)
     &			  , bedporosity , Lsubactdepth(I,J)   
     &			  , Vsandsurfmax(I,J), hfine(I,J) , harea(i,j) 
     &             , AV(I,J) , BV(I,J), znaught(i,j), Z0Min
     &			  , Z0Max , WK_RoughShapeParam , Vsandsubmax(I,J)  
     &            ,  Vsandsub(I,J) , Vsandsurf(I,J), dref(i,j))


!      CALCULATES THE ELEVATION OF THE IMMOBILE GRAVEL BED FROM THE 
!      GIVEN HEIGHT OF SAND FOR WHICH THE HYDRAULIC MODEL WAS GENERATED

             		CALL  INITIMETA( Fracs(i,j) ,harea(i,j),Fsmin ,hmax(i,j) 
     &                   , hfine(i,j), vsandsurf(i,j), vsandsurfmax(i,j)
     &					 , he)  

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
			Sdensity=2.650		  ! density of sediment
	      deltaD=(Sdensity-Wdensity)/Wdensity !sediment submerged specific gravity 

!       ---Parameters in Wilcock & Kenworthy model

!			taustar_rg0=0.035
!			taustar_rg1=0.011
!			taustar_rs1=0.065
!			alpha=1


!			AK=115
!			chi=0.923
!			phi_prime=1.27



		DO 778 i=1,ns
			DO 778 j=1,nn
			    IF(I.lt.SEDBCNODE) THEN
			        QS(I,J) = 0.
			        QN(I,J)=0.
					 GO TO 778
			    ENDIF

!				No Qsand if there is not sand
				IF (ibc(i,j).EQ.0.OR.(Vsandsub(i,j)+Vsandsurf(i,j)).EQ.0.0) THEN   
					 QS(I,J)=0.
					 QN(I,J)=0.
					 GO TO 778
				ENDIF

			
				
				
				
				
				TB=((USTRESS(I,J)**2)+(VSTRESS(I,J)**2))**.5 !shear stress module
				ustar=SQRT(TB/Wdensity)
				
				IF (USTRESS(I,J).eq.0) THEN !vector direction in radiant
					TAN1=0.5*PI                          
				ELSE
					TAN1=ATAN(VSTRESS(I,J)/USTRESS(I,J)) 
				ENDIF
				
				taustar_rs0=alpha*taustar_rg0*(Dg/Dsand)
				taustar_rs=taustar_rs1+
     &                     (taustar_rs0-taustar_rs1)*exp(-14*Fracs(i,j))
				tau_rs=taustar_rs*deltaD*Wdensity*981*Dsand
				phi_sand=TB/tau_rs


				IF (phi_sand.lt.phi_prime) THEN
						Wprime_sand=0.002*(phi_sand)**(7.5)
				ELSE
						Wprime_sand=AK*
     &                               (1-chi/(phi_sand)**(0.25))**(4.5)
				ENDIF
			
				Qsand=Wprime_sand*Fracs(i,j)
     &			           *(ustar)**(3.0)/(deltaD*981)	
 

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
                 QN(i,j)=Qsand*SIN(TAN1) !just for sand
                QS(i,j)=Qsand*COS(TAN1) !just for sand
!			    TB=((USTRESS(I,J)**2)+(VSTRESS(I,J)**2))**.5
!				qs(i,j)=Qsand * ustress(i,j)/TB !just for sand
!				qn(i,j)=qs(i,j) * (vstress(i,j)/ustress(i,j)) !just for sand
                
                check=(abs(QS(i,j)/(rn(I,J)*ds)
     &          -(QN(I,J)/(R(I)*rn(I,J))))+abs(QN(I,J)/dn
     &                ))*1.54*harea(i,j)*tdt


      IF (check.GT.(Vsandsub(i,j)).AND.(Vsandsurf(i,j)).EQ.0.0)THEN
                    
        Qsand=Vsandsub(i,j)/(1.54*harea(i,j)*tdt) !(tdt*dn*1.54)!modified 10/09/09
                
      ELSEIF (check.GT.(Vsandsub(i,j)+Vsandsurf(i,J)) )then

        Qsand=(Vsandsub(i,j)+Vsandsurf(i,j))/(1.54*harea(i,j)*tdt) !(tdt*dn*1.54)!modified 10/09/09
!      write(6,*)'high transport in cell', i,j,'Fine fraction',fracs(i,j)
!      pause
        

      ENDIF

				qn(i,j)=Qsand*SIN(TAN1) !just for sand
				qs(i,j)=Qsand*COS(TAN1) !just for sand
!			    TB=((USTRESS(I,J)**2)+(VSTRESS(I,J)**2))**.5
!				qs(i,j)=Qsand * ustress(i,j)/TB !just for sand
!				qn(i,j)=qs(i,j) * (vstress(i,j)/ustress(i,j)) !just for sand
	
! 25 June 2010- 
!           NO SAND TRANSPORT FROM A WET TO A DRY CELL
            IF (I.LT.NS) THEN
				IF (ibc(i+1,j).EQ.0.AND.QS(I,J).GT.0.0) THEN   
					 QS(I,J)=0.0
					 GO TO 778
			    ENDIF
			ENDIF
			IF (I.gT.1) THEN
			 	IF (ibc(i-1,j).EQ.0.AND.QS(I,J).LT.0.0) THEN 	
				     QS(I,J)=0.0
					 GO TO 778 
				ENDIF
			ENDIF
			IF (J.LT.NN) THEN	  	 
			  	IF (ibc(i,j+1).EQ.0.AND.QN(I,J).GT.0.0) THEN   
					 QN(I,J)=0.0
					 GO TO 778
				ENDIF
		    ENDIF
			IF (J.GT.1) THEN
				IF (ibc(i,j-1).EQ.0.AND.QN(I,J).LT.0.0) THEN 	
				     QN(I,J)=0.0
					 GO TO 778
				ENDIF
			ENDIF
! 25 JUNE 2010 END
	
	
	
	
				IF (j.eq.1.or.j.eq.nn) THEN
					qn(i,j)=0.
					qs(i,j)=0.   
				ENDIF
778			CONTINUE


!        END WILCOCK AND KENWORTHY TRANSPORT MODEL


!         Rate of change of bed or div of Qs
      if(TRANSEQTYPE.eq.2) THEN !12/04/09 Exner equation. 
!                                For WK model the second order is not used but we used a mass balance on the fluxes
!                                this is for accounting for the limited material in a cell

            DO 1650 I=1,ns
            DO 1645 J=1,nn

               SC=rn(I,J)
                ! what is exiting the cell is positive=erosion
                  
              con(i,j)=(abs(QS(I,J)/(ds*rn(I,J))-QN(I,J)/(R(I)*rn(I,J)))
     &                  +abs(QN(I,J)/dn))*harea(i,j)                  
  
                ! end what is exiting the cell
               IF(I.lt.SEDBCNODE) THEN
                       
                 CON(I,J)=0.
                 
c                   CON(I,J)=(QS(2,J)-QS(1,J))/(ds*SC)
c                    GO TO 643
                   GO TO 1645
                    
               ELSEIF(I.EQ.ns) THEN
C                 CON(I,J)=(QS(I,J)-QS(I-1,J))/(ds*SC)
c                 CON(I,J)=(QS(2,J)-QS(1,J))/(ds*SC)
                    CON(I,J)=CON(I-1,J)
                    GO TO 1643
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
1643           IF(J.EQ.1) THEN
                check=QN(I,2)/dn
                if (check.LT.0.0) then  
                    CON(I,1)=CON(I,1)+check*harea(I,2)
                endif
				CON(I,J)=CON(I,J)/(harea(i,j))  ! added 24 January 2012 correct for the lateral nodes
                GO TO 1645
               ELSE
                IF(J.EQ.nn) THEN
                    check=QN(I,nn-1)/dn
                    if (check.GT.0.0) then
                        CON(I,nn)=CON(I,nn)-check*harea(I,nn-1)
                    endif
				CON(I,J)=CON(I,J)/(harea(i,j))  ! added 24 January 2012 correct for the lateral nodes
                    GO TO 1645
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
1645            CONTINUE
1650            CONTINUE           

!     End the Exner for WK

        else

!         Rate of change of bed or div of Qs

         !For non-equilibrium case smooth change in sediment supply at boundary
710   IF(SedEqNode > SedBCNode) THEN
      DO I = SEDBCNODE,SEDEQNODE
       
      NONEQSMOOTH=(((1.0-BCFRACTION)/(SEDEQNODE-SEDBCNODE))
     &*(I-SEDBCNODE))
     & +BCFRACTION
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

            DO 650 I=1,ns
            DO 645 J=1,nn
               SC=rn(I,J)
               IF(I.lt.SEDBCNODE) THEN
                       
                 CON(I,J)=0.
                 
c                   CON(I,J)=(QS(2,J)-QS(1,J))/(ds*SC)
c                    GO TO 643
                   GO TO 645
                    
               ELSEIF(I.EQ.ns) THEN
C                 CON(I,J)=(QS(I,J)-QS(I-1,J))/(ds*SC)
c                 CON(I,J)=(QS(2,J)-QS(1,J))/(ds*SC)
                    CON(I,J)=CON(I-1,J)
                    GO TO 643
               ELSE
                   CON(I,J)=(QS(I+1,J)-QS(I-1,J))/(2.*ds*SC)
               ENDIF
               
643                   CON(I,J)=CON(I,J)-(QN(I,J)/(R(I)*SC))
            
               IF(J.EQ.1) THEN
                   CON(I,1)=CON(I,1)+(QN(I,2)-QN(I,1))/dn
                   GO TO 645
               ELSE
               
                IF(J.EQ.nn) THEN
                    CON(I,nn)=CON(I,nn)+(QN(I,nn)-QN(I,nn-1))/dn
                    GO TO 645
                ENDIF
               
               CON(I,J)=CON(I,J)+(QN(I,J+1)-QN(I,J-1))/(2.*dn)  !end equation 26 in Nelson&Smith AGU12
               
               ENDIF
645            CONTINUE
650            CONTINUE

      endif   !Daniele 12/02/09 endif on the Exner equation

!710			QTOTAL=0.
!
			DO 617 I=1,ns
			QTOT(I)=0.
			DO 616 J=2,nn
616			QTOT(I)=QTOT(I)+(.5*dn)*(QS(I,J)+QS(I,J-1))
            WRITE(12,*) I, QTOT(i)
			QTOTAL=QTOTAL+QTOT(I)/ns
617			CONTINUE
            WRITE(12,*) 'AVERAGE', QTOTAL
            CLOSE(12)
!			DO 650 I=1,ns
!			DO 645 J=1,nn
!			    SC=rn(I,J)
!			    IF(I.lt.SEDBCNODE) THEN
!			            
!	              CON(I,J)=0.
!	              
!c			        CON(I,J)=(QS(2,J)-QS(1,J))/(ds*SC)
!c					GO TO 643
!			        GO TO 645
!!			    ELSEIF(I.EQ.SEDBCNODE) THEN
!!			        CON(I,J)=(QS(i+1,J)-QS(i,J))/(ds*SC)
!!			        GO TO 645
!			    ELSEIF(I.EQ.ns) THEN
!C	              CON(I,J)=(QS(I,J)-QS(I-1,J))/(ds*SC)
!c	              CON(I,J)=(QS(2,J)-QS(1,J))/(ds*SC)
!			         CON(I,J)=CON(I-1,J)
!			         GO TO 643
!			    ELSE
!			        CON(I,J)=(QS(I+1,J)-QS(I-1,J))/(2.*ds*SC)
!			    ENDIF
!			    
!643		            CON(I,J)=CON(I,J)-(QN(I,J)/(R(I)*SC))
!			
!			    IF(J.EQ.1) THEN
!			            CON(I,1)=CON(I,1)+(QN(I,2)-QN(I,1))/dn
!				    GO TO 645
!			    ELSE
!			        IF(J.EQ.nn) THEN
!      			            CON(I,nn)=CON(I,nn)+(QN(I,nn)-QN(I,nn-1))/dn
!			            GO TO 645
!			        ENDIF
!			    
!			    CON(I,J)=CON(I,J)+(QN(I,J+1)-QN(I,J-1))/(2.*dn)  !end equation 26 in Nelson&Smith AGU12
!			    
!			    ENDIF
!645			CONTINUE
!650			CONTINUE

!      DO J = 1,nn
!        CON(SEDBCNODE, J) = BCFRACTION*CON(SEDBCNODE, J)
!      ENDDO

!         CON is div of Qs


c			do 652 j=1,nn
c652			con(1,j)=con(2,j)
			DO 655 I=1,ns
			DO 655 J=1,nn
			CON(I,J)=1.54*CON(I,J)
655			CONTINUE

!         CON is div of Qs/porosity; 1/porosity=1.54;



! smoothing the DIV Qs over the grid

 !           write(6,*) 'con',(con(10,j),j=1,nn)
!             write(6,*) 'rn',(rn(10,j),j=1,nn)
            IF (SEDSMOO.ge.1) then
			DO 670 ISMOO=1,SEDSMOO
			DO 660 I=2,ns-1
			DO 660 J=1,nn
			if(ibc(i,j).eq.0) then
			    con(i,j)=0.
			    go to 660
			endif
			IF(I.LT.SEDBCNODE) THEN
			    CON(i,j) = 0.
			    goto 660
			ENDIF
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
660			continue
!            write(6,*) 'smooth pass', ismoo, (dum1(10,j),j=1,nn)
			DO 670 I=2,ns-1
			DO 670 J=1,nn
670			CON(I,J)=DUM1(I,J)
            endif
!            write(6,*) 'con',(con(10,j),j=1,nn)
!            write(6,*) 'rn',(rn(10,j),j=1,nn)
!            pause
c	 DO 685 ISMOO=1,2
c	 DO 680 I=2,ns-1
c	 DUM1(I,1)=(CON(I-1,1)+CON(I,1)+CON(I+1,1)+con(I,j+1))/4.
c680	 DUM1(I,25)=(CON(I-1,nn)+CON(I,nn)+CON(I+1,nn)+con(i,nn-1))/4.
c	 DO 685 I=2,ns-1
c	 CON(I,1)=DUM1(I,1)
c685	  CON(I,nn)=DUM1(I,nn)
			DO 687 J=1,nn
!			CON(ns,J)=(.2*CON(ns-1,J)+CON(ns,J))/1.2
            con(ns,j)= con(ns-1,j)
687			CON(1,J)=0.

	
	





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
 


	if (TRANSEQTYPE.NE.2) go to   780     !Daniele 	

 
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
				CALL VOLFRAC(av(i,j),bv(i,j), VS, Vsandsub(i,j) 
     &              , Vsandsurf(i,j),Vsandsurfmax(i,j), Vsandsubmax(I,J)
     &              , Fsmin, Fracs(i,j))
				
!             UPDATES THE ELEVATION OF SAND OVER THE GRAVEL 
!             hfine(ns,nn) is the relative change in sand depth from lowest point into the gravel

!		
!	abh = alpha parameter for the Beta model of Sand height vs. Surface fraction of Sand
!	bbh = beta parameter for the Beta model of Sand height vs. Surface fraction of Sand
!	hmax = maximum height of the gravel surface  !!
!	
!	
	

				CALL ELEVDT(abh(i,j), bbh(i,j), hmax(i,j)
     &                      , Vsandsurf(I,J), Vsandsurfmax(i,j) 
     &                  , Fsmin, Fracs(I,J), harea(i,j), hfine(I,J), he)

			if(((imeta(i,j)+he)-eta(i,j)).GT.hl(i,j)) then
			    write(6,*)'Error in Sediment Transport', 
     &                        i, j, hl(i,j)/con(i,j)
			ENDIF
			
			tmpdepth = hl(i,j)+con(i,j)*dt
			if(tmpdepth < hmin) then
			    eta(i,j) = e(i,j)-hmin
			    hl(i,j) = hmin
			else
			    eta(i,j)=imeta(i,j)+he
			    hl(I,J)=hl(I,J)+he
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

				    CALL roughdt( Fracs(I,J) , Fsmin , Z0Min 
     &				        , Dref(i,j) , hmax(i,j) 
     &                      , hfine(i,j) , WK_RoughShapeParam, Drough )
    			    		
				        znaught(i,j) = Drough
				ENDIF

779       ENDDO		
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

	if (TRANSEQTYPE.eq.2) go to 2000         !Daniele 	

780			DO 720 I=1,ns
			DO 720 J=1,nn
			if(ibc(i,j).eq.0) then
			 con(i,j)=0.
			endif
			if(dt*con(i,j) > hl(i,j)) then
			    write(6,*)'Error in Sediment Transport deposition'
     &                       , i, j, hl(i,j)/con(i,j)
			ENDIF
			if(dt*abs(con(i,j)) > hl(i,j)) then
			    write(6,*) 'Error in Sediment Transport erosion'
     &    , i,j,ibc(i,j), ibc(i,j), cd(i,j), cd(i-1,j), hl(i,j),con(i,j)
               
			endif
			tmpdepth = hl(i,j)+(dt*con(i,j))
			if(ibc(i,j).ne.0.and.tmpdepth.le.hmin) then
			    eta(i,j) = e(i,j)-hmin
			    hl(i,j) = hmin
			else if (dt*abs(con(i,j)) > hl(i,j)) then
			    if(con(i,j).lt.0) then
			        eta(i,j) = eta(i,j)+(TSFracDepth*hl(i,j))
			        hl(i,j) = hl(i,j)-(TSFracDepth*hl(i,j))
			    else
			        eta(i,j) = eta(i,j)-(TSFracDepth*hl(i,j))
			        hl(i,j) = hl(i,j)+(TSFracDepth*hl(i,j))			    
			    endif
			else
			    eta(i,j)=eta(i,j)-(dt*con(i,j))
			    hl(I,J)=hl(I,J)+(dt*con(I,J))
			endif
!			eta(i,j)=eta(i,j)-(dt*con(i,j))
!			hl(I,J)=hl(I,J)+(dt*con(I,J))
!			if(hl(i,j).le.2) then
!			 hl(i,j)=2.
!			endif
720			CONTINUE

!			    Write(11,2010) RUNID
!			    write(11,*) Q,MO
!			    write(11,*) R,W,eta
!			    write(11,*) xo,yo
	if(roughnesstype == 1) then
		CALL ZOTOCDTWO()
	endif
2000			if(nct.eq.nsteps) then
			    WRITE(7,2010) RUNID
			    WRITE(7,*) QS,QN,CON,hl,QTOT
			    write(7,*) ((i,j,con(i,j),j=1,nn),i=10,15)
!			    write(10,*) hwt,erelax,urelax,arelax,itm,cd,q,mo
!			    write(10,2010) runid
!			    write(10,*) q,mo,hav
!			    write(10,*) r,w,eta,ibc
!			    write(10,*) xo,yo
			endif
2010			format(A20)
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

        REAL(KIND=8), INTENT(IN) :: d,nu,s,g,p,csf
        REAL(KIND=8) :: dstar,logdstar,r1,r2,r3,tanhdstar,wstar
        REAL(KIND=8), PARAMETER :: one=1.0

        dstar = (s-one)*g*d**3/(nu*nu)
        logdstar = LOG10(dstar)
        tanhdstar = DTANH(logdstar-4.6D0)
        r1 = -3.76715D0+1.92944*logdstar-9.815D-2*logdstar*logdstar- 
     &       5.75D-3*logdstar**3+5.6D-4*logdstar**4
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
            REAL(KIND=8), INTENT(IN) :: D, KV
            REAL, INTENT(OUT) :: TC
            
            REAL :: DSTAR, TAUCRIT
            REAL :: RELDEN, G
            
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
