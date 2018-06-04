MODULE RivStagr4Mod_jmn

	USE EInitMod
	USE RivVarMod
	USE CalcCond
	USE CSedMod_DT_SUSP
	USE CSedMod
	USE ReadCGNS
	USE RivVarWMod
	USE RivVarVertMod
	USE WriteCGNS
	USE RivVertMod
	USE StressDivMod
!	USE CSedMod
    USE UInitMod
	USE TRIDIAG
	
	USE RivRoughnessMod
	USE RivVarTimeMod
	USE RivConnectivityMod
	IMPLICIT NONE

	CONTAINS 
	
	SUBROUTINE STAGR4(STR_IN )
		CHARACTER(*), INTENT(IN) :: STR_IN
	character*20 runid
	INTEGER :: solIndex
	INTEGER :: n, n2
	REAL ::	pi, g, rho, vkc, ustar
	REAL :: dinc
	REAL :: area
	REAL :: ba, bb, bc, be, bf, bg, bot, bot2
	REAL :: be1, be2, bf1, bf2
	REAL :: va, ta, tb, tc, td, te, tf 
	REAL :: ua, uu, vv
	REAL :: dube
	REAL :: sumqpdiff
	INTEGER, PARAMETER :: num_string = 4
	INTEGER :: i, j, nct, jswi, j1, j2, j3
	REAL :: plinc
	LOGICAL :: lwet, doagain

!	output vars
	REAL :: rcos, rsin, ux, uy, uxnew, uynew, qx, qy
	REAL :: dhav2

	REAL :: rnr1, rnrm1, rnrm1j, rnrm1ij !Changed as per CD
	REAL :: dsvpimj, dsvmimj, ekapimj, ekamimj
	INTEGER :: ip, im, jp, jm
	
	REAL :: hlrewe
    REAL :: dr
    
	INTEGER :: imod, iwetdry
	INTEGER :: nclusters, nwetnodes, noldwetnodes
	REAL :: nodechange
	REAL :: dDisch, dStage, newStage, oldStage, newDisch, oldDisch, dsstage
!   tmpVars
	REAL :: tmpLEV, tmpvardt
	INTEGER :: ier

	CHARACTER(20) text
	INTEGER(4) iret
	INTEGER :: ibccount
	
 !    OPEN(11,FILE='topser')
    nwetnodes = 0
	errorcode = 0
    uvinterp = 0
	pi=acos(-1.)
	g=980.
	rho=1.
	vkc=0.40
    solIndex = 0
	CALL WELCOME
	CALL read_iRIC_CGNS(STR_IN)
	   
	   imod = hiterinterval
	   if(hcalcwetting == .TRUE.) then
		iwetdry = 1
	   else
		iwetdry = 0
	   endif


3	format(A20)
	CALL Calc_Area(x, y, xo, yo, nm, dn, harea)

!pause
65	format(25f6.2)
    if(varDischType == 1) then !Discharge Time Series
        CALL getInterpTimeSeriesValue(1, VarDischStartTime, newDisch)
        q = newDisch*1e6
    else
        newDisch = q/1e6
    endif
    
    if(varStageType == 1) then !Stage Time Series
    !in this case the times series is stored in the rating curve
        CALL getInterpRatingCurveValue(1, VarDischStartTime, newStage)
!        newStage = (newStage-elevoffset) * 100.
    else if (varStageType == 2) then !Stage Rating Curve
        Call getInterpRatingCurveValue(1, newDisch, newStage)
!        newStage = (newStage-elevoffset) * 100.
    endif
    
!    if(cdtype == 1) then !variable roughness
!        CALL setRoughness(ns, nn, roughnessType, newDisch, cd, znaught)
!    endif

!    CALL Label_clusters(-1, ibc, nclusters, icon)
!    CALL Delete_clusters(ibc, icon, nclusters, nwetnodes)
    noldwetnodes = nwetnodes

	if(HotStart.eq..FALSE.) then 
		CALL einit(e, hl, eta, ibc, w, hav)
		call uinit(u,v,hav,w,eta,q,ibc)
		if(errorcode < 0) then
			!call write_error_CGNS(STR_IN)
            CALL dealloc_common2d()
            if(vbc) then
                CALL DEALLOC_VELBC()
            endif
            if(CALCQUASI3D) then
                CALL dealloc_common3d()
                if(TRANSEQTYPE) THEN
                    CAll dealloc_csed3d_dt()
                ENDIF
            endif
            CALL dealloc_working()
            !pause
 			return
		endif
	endif
		
	if(roughnesstype == 1) then !Roughness as Z0
!		CALL ZOTOCDTWO() !Calculate cd from log profile
		CALL Z0TOCD() !Calculate cd from 2-part profile
	endif
	
    nsteps = (VarDischEndTime-VarDischStartTime)/dt

	if(IterationOut) Then
	    if(itm > nsteps) then
	        CALL alloc_TSNames(itm)
	    else    
	        CALL alloc_TSNames(nsteps+10)
	    endif
	else
        CALL alloc_TSNames(nsteps+10)
    endif
    
    
!    IF(HOTSTART) then
!        oldStage = newstage
!        dsstage = newstage
!    ELSE
        oldStage = hav(ns)
        dsstage = hav(ns)
!    ENDIF
    newStage = oldStage
    oldDisch = q
    nct = -1
    totTime = VarDischStartTime
    ptime = totTime+iplinc*dt
    
 4   do while(totTime <= VarDischEndTime)
        nct = nct+1
        
        if(nct > 0.and.solType == 0) then
            Exit
        endif
        if(nct > 0.and.VarDischEndTime.eq.VarDischStartTime) Then
            EXIT
        ENDIF
        
        if(nct == 0) then
            vardt = dt
        else
            if(vardt > dt) then
                vardt = dt
             endif
        endif
        
        if(CALCCSED.eq..FALSE.) then
            vardt = dt
        endif

        oldDisch = q
        oldStage = newStage
        
        if(nct > 0) then
            totTime = totTime+vardt
            if(totTime > VarDischEndTime) then
                RETURN
            endif
        endif
        
       if(varDischType == 1) then !Update Discharge
            CALL getInterpTimeSeriesValue(1, totTime, newDisch)
            q = newDisch*1e6
        endif
        dDisch = q-oldDisch
        
        if(varStageType == 1) then !Stage Time Series
            CALL getInterpRatingCurveValue(1, totTime, newStage)
!            newStage = newStage * 100.
        else if (varStageType == 2) then !Stage Rating Curve
            Call getInterpRatingCurveValue(1, newDisch, newStage)
!            newStage = newStage * 100.
        endif
        
        dStage = newStage-oldStage
        
!        if(cdtype == 1) then !variable roughness
!            CALL setRoughness(ns, nn, roughnessType, newDisch, cd, znaught)
!        endif
        
        if(varDischType == 1) then
        !!  Adjust water surface elevation on all wet nodes by dStage !!
        !!  Check for minimum depth and turn nodes on or off.         !!
            CALL VarDischEInit(e, hl, u, v, eta, ibc, dStage, dsstage)
        endif
 
!        if(vbcds == 1) then
!            CALL VarDischUBInit(u, v, w, eta, e, ibc, q)
!        endif
        
        !! TURN NODES ON - CHECK RE-WETTING !!
        Call UpdateWetting(e,hl,u,v,eta, ibc, dsstage)
        
        !! SET UPSTREAM BOUNDARY VELOCITY  !!
        if(varDischType == 1) then
            CALL VarDischUInit(u, v, w, eta, e, ibc, q, dStage)
        endif
                
        !! FINDS ISOLATED POCKETS OF WET NODES AND DELETES THEM
        CALL Label_clusters(-1, ibc, nclusters, icon)
        CALL Delete_clusters(ibc, icon, nclusters, nwetnodes)
        
        !!  CHECK FOR IBC == 6 OR IBC == 4  !!
        CALL updateIBC()

        !!  IF GRID EXTENSION THEN RESET EXTENSION TOPOTRAPHY !!
!        if(nct > 1) then
!            CALL RESETGRIDEXTENSION()
!        endif
	if(roughnesstype == 1) then !Roughness as Z0
!		CALL ZOTOCDTWO() !Calculate cd from log profile
		CALL Z0TOCD() !Calculate cd from 2-part profile
	endif
			       

 	    DO i = 1,10
 	        write(*,*)
 	    ENDDO
 	    
	    write(*,*)
        Write(*,*) 'Variabledt', tmpvardt, 'DT', dt
        write(*,*) 'TIME STEP', nct, 'TIME', tottime, 'PrintTime', ptime
        write(*,*)
        IF(CALCQUASI3D) THEN
            write(*,*) 'Quasi-3D Calculation is on'
            if(CALCQUASI3DRS) then
                write(*,*) 'Streamline Curvature is on'
            ENDIF
        ELSE
            write(*,*) 'Quasi-3D Calculation is off'
        ENDIF
        write(*,*)
        IF(CALCCSED) THEN
            write(*,*) 'Sediment-transport is on'
            IF(CALCSEDAUTO) THEN
                write(*,*) 'Time-stepping method is automatic'
            ELSE
                write(*,*) 'Time-stepping method is fixed'
            ENDIF
            write(*,*) 'Old wet nodes: ', noldwetnodes, ' New wet nodes: ', nwetnodes

        ELSE 
            write(*,*) 'Sediment-transport is off'
        ENDIF
            
	    write(*,*)        
        write(*,*) 'Discharge', q/1e6, 'Change Discharge', dDisch/1e6
        write(*,*) 'Stage', newStage/100.+elevoffset, 'Change Stage', dStage/100.
	    write(*,*)        
        
        if(nct.ge.1) then
            itm=interitm
            nodechange = abs((real(noldwetnodes)-real(nwetnodes)))
            itm = real(nodechange*interitm)
            IF(itm < interitm)THEN
                itm = interitm
            else if (itm >= 5*interitm) then
                itm = maxInterIterMult*interitm
            ENDIF
            noldwetnodes = nwetnodes
        endif 

	    if(itm == 0) goto 601
	    ITER_LOOP: DO iter = 1,itm  
    
        IF(iwetdry.eq.1.and.MOD(iter,imod).eq.0.and.varDischType == 0.and.iter.lt.hiterstop.and.nct.lt.1)then
!        IF(iwetdry.eq.1.and.MOD(iter,imod).eq.0.and.varDischType == 0)then
            CALL    UpdateWETTING(e, hl, u, v, eta, ibc, dsstage)
        ENDIF
       
        CALL UPDATEIBC()

        DO I = 1,NS
            DO J = 1,NN
	            if(ibc(i,j).eq.1) then
		            ustar=(totcd(i,j)*(u(i,j)**2.))**.5
	            else
		            ustar=(totcd(i,j)*(u(i,j)**2.+v(i,j)**2.))**.5
	            endif

                !eka(i,j)=.005*200.*200.+vkc*hl(i,j)*ustar/6.
		        IF(LEVType == 0) THEN
			        eka(i,j) = evc+vkc*hl(i,j)*ustar/6.
		        ELSE IF(LEVType.ne.0.and.NCT.eq.0) THEN 
			        IF(iter < LEVBegIter) THEN
				        tmpLEV = StartLEV
			        ELSEIF (iter >= LEVBegIter .and. iter < LEVEndIter) THEN
				        tmpLEV = &
				        (((StartLEV-EndLEV)/(LEVBegIter-LEVEndIter))*&
					        (iter-LEVBegIter))+StartLEV
			        ELSEIF (iter > LEVEndIter) THEN
				        tmpLEV = EndLEV
			        ENDIF
			        eka(i,j) = tmpLEV+vkc*hl(i,j)*ustar/6.
		        
		        ELSE
		            eka(i,j) = EndLEV+vkc*hl(i,j)*ustar/6.
		        ENDIF
            ENDDO
        ENDDO

    
	do 100 i=2,ns
        do 100 jswi=1,2
        if(jswi.eq.1) then
         j1=1
         j2=nn
         j3=1
        else
         j1=nn
         j2=1
         j3=-1
        endif
	do 100 j=j1,j2,j3
        if(ibc(i,j).eq.2) then
         va=(v(i-1,j)+v(i,j))/4.
        else if(j.lt.nn) then
         va=(v(i-1,j)+v(i,j)+v(i,j+1)+v(i-1,j+1))/4.
        endif
	if(ibc(i,j).eq.4.or.ibc(i,j).eq.0) then
	 dude(i,j)=0.
	 u(i,j)=0.
	 go to 100
	endif
	
    if(u(i,j).ge.0.or.i.eq.ns) then
	    ba=u(i,j)/dsu(i,j)
	    ta=ba*u(i-1,j)
    else
        ba=-1.*u(i,j)/dsu(i+1,j)
        ta=ba*u(i+1,j)
    endif
    
    if((va.ge.0.and.ibc(i,j).ne.1).or.ibc(i,j).eq.2) then
        bb=va/dnu(i,j)
        tb=bb*u(i,j-1)
    else
        bb=-1.*va/dnu(i,j+1)
        tb=bb*u(i,j+1)
    endif
    
    bc=-1.*va/(rn(i,j)*r(i))
    if(ibc(i-1,j).eq.0) then
        td = 0.
    else
	    td=-1.*g*(e(i,j)-e(i-1,j))/dse(i,j)
	endif
    
    if(i.eq.ns) then
        be=(2./(dsu(i,j)))*eka(i,j)/dsu(i,j)
        te=be*u(i-1,j)
    else
        be=(4./(dsu(i,j)+dsu(i+1,j)))*((eka(i,j)/dsu(i,j))+&
            (eka(i+1,j)/dsu(i+1,j)))
        te=(4./(dsu(i,j)+dsu(i+1,j)))*((eka(i,j)*u(i-1,j)/dsu(i,j))+&
            (eka(i+1,j)*u(i+1,j)/dsu(i+1,j)))
    endif
    
    if(ibc(i,j).eq.2) then
         bf=(2./dnu(i,j))*(totcd(i,j)*abs(u(i,j))+eka(i,j)/dnu(i,j))
         tf=(2./dnu(i,j))*(eka(i,j)*u(i,j-1)/dnu(i,j))   
    elseif(ibc(i,j).eq.1) then
         bf=(4./(dnu(i,j)+dnu(i,j+1)))*((eka(i,j+1)/dnu(i,j+1))+&
            totcd(i,j)*abs(u(i,j)))
         tf=(4./(dnu(i,j)+dnu(i,j+1)))*(eka(i,j+1)*u(i,j+1)/dnu(i,j+1)) 
    else
         bf=-1.*eka(i,j)/(dnu(i,j)*rn(i,j)*r(i))+(2./(dnu(i,j)+&
            dnu(i,j+1)))*((eka(i,j+1)/dnu(i,j+1))+(eka(i,j)/dnu(i,j)))
         tf=-1.*eka(i,j)*u(i,j-1)/(dnu(i,j-1)*rn(i,j-1)*r(i))+&
            (2./(dnu(i,j)+dnu(i,j+1)))*((eka(i,j+1)*u(i,j+1)/dnu(i,j+1))+&
            (eka(i,j)*u(i,j-1)/dnu(i,j)))
    endif
	
	bot2=rho*(hl(i,j)+hl(i-1,j))
	
	if(bot2.eq.0.) then
        bg=10.**12.
	else
	    bg=2.*(totcd(i,j)*(u(i,j)**2.+va**2.)**.5)/bot2
	endif

	bot=ba+bb+bc+be+bf+bg

    if(abs(bot).gt..000001) then
	    u(i,j)=u(i,j)*(1.-urelax)+urelax*(ta+tb+td+te+tf)/bot
	    if(i.eq.ns.and.u(i,j).lt.0) then
	        if(vbcds .ne. 0) then
	            u(i,j) = 0.
	        endif
	    endif
	    dude(i,j)=-1.*g/(dse(i,j)*bot)   
    endif

100     continue

        do 120 i=2,ns
        do 120 jswi=1,2
        if(jswi.eq.1) then
         j1=1
         j2=nn
         j3=1
        else
         j1=nn
         j2=1
         j3=-1
        endif
        do 120 j=j1,j2,j3
        if(ibc(i,j).eq.0.or.ibc(i,j).eq.1) then
         v(i,j)=0.
         dvde(i,j)=0.
         go to 120
        endif
        if(i.eq.ns) then
         ua=(u(i,j)+u(i,j-1))/2.
        else
         ua=(u(i,j)+u(i,j-1)+u(i+1,j)+u(i+1,j-1))/4.
        endif
        if(ua.ge.0.or.i.eq.ns) then
         ba=ua/dsv(i,j)
         ta=ba*v(i-1,j)
        else
         ba=-1.*ua/dsv(i+1,j)
         ta=ba*v(i+1,j)
        endif
        if(v(i,j).ge.0.or.ibc(i,j).eq.2) then
         bb=v(i,j)/dnv(i,j)
         tb=bb*v(i,j-1)
        else
         bb=-1.*v(i,j)/dnv(i,j+1)
         tb=bb*v(i,j+1)
        endif
       
        tc=-1.*(ua**2.)/(rn(i,j)*r(i))
        if(ibc(i,j-1).eq.0) then
            td = 0.
        else
            td=-1.*g*(e(i,j)-e(i,j-1))/dne(i,j)
        endif
        if(i.eq.ns) then
         be=(1./dsv(i,j))*eka(i,j)/dsv(i,j)
         te=be*v(i-1,j)
        else
         be1=(2./(dsv(i,j)+dsv(i+1,j)))*eka(i+1,j)/dsv(i+1,j) 
         be2=(2./(dsv(i,j)+dsv(i+1,j)))*eka(i,j)/dsv(i,j)
         be=be1+be2
         te=be1*v(i+1,j)+be2*v(i-1,j) 
        endif
        if(ibc(i,j).eq.2) then
         bf1=(1./dnv(i,j))*eka(i,j)/dnv(i,j)
         bf2=bf1
         bf=bf1+bf2
         tf=bf1*v(i,j-1)
        else
         bf1=(2./(dnv(i,j)+dnv(i,j+1)))*eka(i,j+1)/dnv(i,j+1)
         bf2=(2./(dnv(i,j)+dnv(i,j+1)))*eka(i,j)/dnv(i,j)
         bf=bf1+bf2 
         tf=bf1*v(i,j+1)+bf2*v(i,j-1)
        endif
	 bot2=rho*(hl(i,j)+hl(i,j-1))
	 if(bot2.eq.0.) then
	 bg=10.**12.
	 else
		  bg=2.*(totcd(i,j)*(ua**2.+v(i,j)**2.)**.5)/bot2
	 endif
	 bot=ba+bb+be+bf+bg
		  if(abs(bot).gt..000001) then
		  if(vbcds .eq. 0) then
	            v(i,j)=v(i,j)*(1.-urelax)+urelax*(ta+tb+tc+td+te+tf)/bot
	      else
	            if(i.ne.ns) then
	                v(i,j)=v(i,j)*(1.-urelax)+urelax*(ta+tb+tc+td+te+tf)/bot
	            endif
	      endif
	            dvde(i,j)=-1.*g/(dn*bot)
		  endif
120	continue


        do 150 i=1,ns
        do 150 j=1,nn
        if(i.eq.1) then
        area=dnq(i,j)*hl(i,j)
        else 
        area=dnq(i,j)*(hl(i,j)+hl(i-1,j))/2.
        endif
        if(j.eq.1.or.j.eq.nn) area=area/2.
        qu(i,j)=u(i,j)*area
        dqsde(i,j)=dude(i,j)*area
        if(j.eq.1) then
        area=dsq(i,j)*hl(i,j)
        else
!        area=.5*(dsq(i,j-1)+dsq(i,j))*(hl(i,j)+hl(i,j-1))/2.
        area=dsq(i,j)*(hl(i,j)+hl(i,j-1))/2.

        endif
        qv(i,j)=v(i,j)*area
        dqnde(i,j)=dvde(i,j)*area
150     continue
        do 200 i=1,ns-1
        do 200 j=1,nn
        if(j.eq.nn) then
         dq(i,j)=qu(i+1,j)-qu(i,j)-qv(i,j)
        else
         dq(i,j)=qu(i+1,j)-qu(i,j)+qv(i,j+1)-qv(i,j)
        endif
200     continue
        do 250 i=1,ns-1
        do 250 j=1,nn
        if(j.eq.nn) then
         dqe(i,j)=dqsde(i,j)+dqsde(i+1,j)+dqnde(i,j)
        else
         dqe(i,j)=dqsde(i,j)+dqsde(i+1,j)+dqnde(i,j)+dqnde(i,j+1)
        endif
        if(abs(dqe(i,j)).lt..000001) then
         de(i,j)=0.
        else
         de(i,j)=dq(i,j)/dqe(i,j)
        endif
250     continue
        do 300 i=ns-1,2,-1
        n=0.
        do 280 j=1,nn
!c       if(abs(dqe(i,j)).gt..000001) then
        if(abs(dqe(i,j)).gt..000001.and.hl(i,j).gt.hmin) then
        n=n+1
        if(j.eq.1) go to 260
        am(n)=-1.*dqnde(i,j)
260     bm(n)=dqe(i,j)
        if(j.eq.nn) then
         ccm(n)=0.
        else
         ccm(n)=-1.*dqnde(i,j+1)
        endif
        dm(n)=dq(i,j)+dqsde(i+1,j)*de(i+1,j)+dqsde(i,j)*de(i-1,j)
        isw(j)=1
        else
        isw(j)=0
        endif
280     continue
        if(n.ge.1) then
        call tridag(am,bm,ccm,dm,em,n)
		if(errorcode.eq.-1) then
			!call write_error_CGNS(STR_IN)
			CALL dealloc_common2d()
            if(vbc) then
                CALL DEALLOC_VELBC()
            endif

            if(CALCQUASI3D) then
                CALL dealloc_common3d()
                if(TRANSEQTYPE) THEN
                    CAll dealloc_csed3d_dt()
                ENDIF
            endif
			CALL dealloc_working()
			WRITE(6,*) 'Error at dqs', 'i = ', i

			!pause
			return
		endif
        endif
        n2=0
        do 300 j=1,nn
        if(isw(j).eq.1) then
         n2=n2+1
         de(i,j)=em(n2)
        else
         de(i,j)=0.
        endif
300     continue 
        do 400 j=1,nn
        n=0
				do 380 i=ns-1,2,-1
					if(abs(dqe(i,j)).gt..000001.and.hl(i,j).gt.hmin) then
					n=n+1
					am(n)=-1.*dqsde(i,j)
					bm(n)=dqe(i,j)
					ccm(n)=-1.*dqsde(i+1,j)
					if(j.eq.1) then
					 dm(n)=dq(i,j)+dqnde(i,j+1)*de(i,j+1)
					elseif(j.eq.nn) then
					 dm(n)=dq(i,j)+dqnde(i,j)*de(i,j-1)
					else
					dm(n)=dq(i,j)+dqnde(i,j)*de(i,j-1)+dqnde(i,j+1)*de(i,j+1)
					endif
					isw(i)=1        
					else
					isw(i)=0
					endif
380				continue
			if(n.ge.1) then
				call tridag(am,bm,ccm,dm,em,n)
		        if(errorcode.eq.-1) then
			        !call write_error_CGNS(STR_IN)
			        CALL dealloc_common2d()
                    if(vbc) then
                        CALL DEALLOC_VELBC()
                    endif

                    if(CALCQUASI3D) then
                        CALL dealloc_common3d()
                        if(TRANSEQTYPE) THEN
                            CAll dealloc_csed3d_dt()
                        ENDIF
                    endif
			        CALL dealloc_working()
					        WRITE(6,*) 'Error at dqn'

			        !pause
			        return
		        endif
            endif
	    n2=0
		do 400 i=ns-1,2,-1
			if(isw(i).eq.1) then
				n2=n2+1
				de(i,j)=em(n2)
			else
				de(i,j)=0.
			endif
400     continue
        do i=1,ns
        qp(i)=0.
        eav(i)=0.
        ibccount = 0
        do j=1,nn
!        qp(i)=qp(i)+qu(i,j)
!        e(i,j)=e(i,j)+erelax*de(i,j)
        if(ibc(i,j).ne.0) then
            e(i,j)=e(i,j)+erelax*de(i,j)
            ibccount = ibccount+1
            eav(i)=eav(i)+e(i,j)
        endif
        qp(i)=qp(i)+qu(i,j)
        hl(i,j)=e(i,j)-eta(i,j)
        if(ibc(i,j).eq.0.or.hl(i,j).le.hmin) then
            if(dryType.eq.0) then
            if(ibc(i,j).ne.0) then
                u(i,j) = 0
                v(i,j) = 0
            endif
               ibc(i,j)=0
            endif
           hl(i,j)=hmin
        endif
        !UNSURE OF THE EFFECT OF THE CODE BLOCK BELOW rmcd 12/23/05
        !Change rmcd made sure i>1 4/30/10
        if(i.gt.1) THEN
        if(abs(u(i,j)).lt..00001) then
         u(i,j)=0.
        endif
        if(abs(v(i,j)).lt..00001) then
         v(i,j)=0.
        endif
        ENDIF
        !!!!!!!!!!!!
        
        enddo
            eav(i)=eav(i)/ibccount
            do j=1,nn
                if(ibc(i,j).eq.0)then
                    e(i,j) = eav(i)
                endif
            enddo

        enddo
500     continue

        do 550 j=1,nn
550         u(1,j)=u(1,j)*q/qp(1)

        dea(ns)=0.
        do  i=ns-1,1,-1
            dinc=arelax*erelax*abs(eav(i)-eav(i+1))*(1.-qp(i+1)/q)
            dea(i)=dea(i+1)+dinc
            do  j=1,nn
            !! code change by rmcd !!
                if(ibc(i,j).ne.0)then
600                 e(i,j)=e(i,j)+dea(i)
!                else
!                    e(i,j) = eav(i)+dea(i) !! added so that e(i,j) in dry nodes follows 
!                                           !! average eav(i) during time-dependent runs
                 endif
            !! end code change by rmcd !!
            enddo
        enddo
		sumqpdiff = 0.
        do 70 i=1,ns
			sumqpdiff = sumqpdiff + (((qp(i)-q)/q*100.))**2.
70      continue
	
		qprms(iter) = sqrt(sumqpdiff/ns)
		if(isnan(qprms(iter))) then
			errorcode=-1
			!call write_error_CGNS(STR_IN)
!			call dealloc_all()
           CALL dealloc_common2d()
            if(vbc) then
                CALL DEALLOC_VELBC()
            endif

            if(CALCQUASI3D) then
                CALL dealloc_common3d()
                if(TRANSEQTYPE) THEN
                    CAll dealloc_csed3d_dt()
                ENDIF
            endif
              CALL dealloc_working()
		WRITE(6,*) 'Error at qprms'
 			!pause

			return

		endif
        if(nct == DbgTimeStep) then
                if(DEBUGSTOP.eq.1.and.iter == DbgIterNum) then
                    solIndex = solIndex+1
                    CALL CG_IRIC_WRITE_SOL_TIME_F(tottime, ier)

                    CALL Write_CGNS2(tottime, q)
!			        CALL write_TimeStep_CGNS(STR_IN, nct, tottime) 
!			        CALL write_timeiter_cgns(STR_IN)
                    CALL dealloc_all()
 	                !pause

	                return
			    endif
		endif
!        IF(LEVType == 0) THEN
!!		    WRITE(*,65, ADVANCE='NO') 'Iteration:', iter, ' Mean Error on Discharge:', qprms(iter)
!            write(6,68)'+ '
!		    WRITE(6,66,ADVANCE='YES')'+ ', 'Iteration: ', iter,' Mean Error on Discharge:', qprms(iter)
!!		    pause
!		ELSE
!		    write(6,68)'+ '
!		    WRITE(*,67, ADVANCE='NO')'+ ', 'Iteration:', iter, ' LEV:', tmpLEV,' Mean Error on Discharge: ', qprms(iter)
!	    ENDIF
        IF(LEVType == 0) THEN
!		    WRITE(*,65, ADVANCE='NO') 'Iteration:', iter, ' Mean Error on Discharge:', qprms(iter)
!            write(6,68)'+ '
		    WRITE(*,*), 'Iteration: ', iter,' Mean Error on Discharge:', qprms(iter)
!		    pause
		ELSE
!		    write(6,68)'+ '
		    WRITE(*,*)'Iteration:', iter, ' LEV:', tmpLEV,' Mean Error on Discharge: ', qprms(iter)
	    ENDIF
	    
	    IF(IterationOut.and.mod(iter,IterPlOut).eq.0) THEN

        do i = 1,ns
            do j = 1,nn
                if(ibc(i,j).eq.0) then
                    v(i,j) = 0.
                    u(i,j) = 0.
                endif
            enddo
        enddo
        if (uvinterp == 0) then

        do  i=1,ns
            do  j=1,nn
!                if(ibc(i,j).eq.0.or.hl(i,j).le.hmin) then
                if(ibc(i,j).ne.-1.or.hl(i,j).le.hmin) then
                    u(i,j)=0.
                    vout(i,j)=0.
                    go to 888
                endif
                if(i.eq.1) then 
                    vout(i,j)=v(i,j)
                    go to 888
                endif

                if(j.eq.nn) then
                vout(i,j)=(v(i-1,j)+v(i,j))/4.
                elseif(j.eq.1) then
                vout(i,j) = 0.
                elseif(ibc(i,j).ne.0.and.ibc(i,j+1).eq.0) then
                vout(i,j)=v(i,j)
                elseif(ibc(i,j).ne.0.and.ibc(i,j-1).eq.0) then
                vout(i,j)=v(i,j+1)
                else
                vout(i,j)=(v(i-1,j)+v(i,j)+v(i,j+1)+v(i-1,j+1))/4. !jmn
                endif
                
        888     if(i.eq.1) then
                    dube=(cd(i,j))*sqrt(u(i,j)**2+vout(i,j)**2)        
                else
                    dube=(cd(i,j)+cd(i-1,j))*.5*sqrt(u(i,j)**2+vout(i,j)**2)
                endif
                taus(i,j)=dube*u(i,j)
                taun(i,j)=dube*vout(i,j)
            enddo
        enddo
    else
        do  i=1,ns
            do  j=1,nn
!                if(ibc(i,j).eq.0.or.hl(i,j).le.hmin) then
                if(ibc(i,j).ne.-1.or.hl(i,j).le.hmin) then
                    uout(i,j)=0.
                    vout(i,j)=0.
                    go to 889
                endif
                if(i.eq.1) then 
                    vout(i,j)=v(i,j)
                    uout(i,j) = u(i,j)
                    go to 889
                endif

                if(j.eq.nn) then
                vout(i,j)=0.0
                elseif(j.eq.1) then
                vout(i,j) = 0.0
!                elseif(ibc(i,j).ne.0.and.ibc(i,j+1).eq.0) then
!                vout(i,j)=v(i,j)
!                elseif(ibc(i,j).ne.0.and.ibc(i,j-1).eq.0) then
!                vout(i,j)=v(i,j+1)
                else
                vout(i,j) = (v(i,j)+v(i,j-1))/2.
                uout(i,j) = (u(i,j)+u(i-1,j))/2. 
                endif
                
        889     if(i.eq.1) then
                    dube=(cd(i,j))*sqrt(uout(i,j)**2+vout(i,j)**2)        
                else
                    dube=(cd(i,j)+cd(i-1,j))*.5*sqrt(uout(i,j)**2+vout(i,j)**2)
                endif
                taus(i,j)=dube*uout(i,j)
                taun(i,j)=dube*vout(i,j)
            enddo
        enddo
    endif
            solIndex = solIndex+1
            CALL CG_IRIC_WRITE_SOL_TIME_F(tottime, ier)
            CALL Write_CGNS2(tottime, q)

   ENDIF

        END DO ITER_LOOP
  
64      format(25f7.3)
!66      format(A15,I6,A25,F15.5)
66      format(A, A, I5,A, F15.5)
67      format(A, A, I5,A, F15.5,A, F15.5)
68      format(A)
75      format(25f7.0)
85      format(25f6.4)
86      format(25f8.3)
84      FORMAT(3i5,f6.2,2f12.1,2g12.4)

        do i = 1,ns
            do j = 1,nn
                if(ibc(i,j).eq.0) then
                    v(i,j) = 0.
                    u(i,j) = 0.
                endif
            enddo
        enddo
601 continue
    if (uvinterp == 0) then
      do 87 i=1,ns
        do 87 j=1,nn
!          if(ibc(i,j).eq.0.or.hl(i,j).le.hmin) then
          if(ibc(i,j).ne.-1) then
           u(i,j)=0.
           vout(i,j)=0.
           go to 88
        endif
        if(i.eq.1) then 
         vout(i,j)=v(i,j)
         go to 88
        endif

        if(j.eq.nn) then
         vout(i,j)=(v(i-1,j)+v(i,j))/4.
        elseif(j.eq.1) then
         vout(i,j) = 0.
        elseif(ibc(i,j).ne.0.and.ibc(i,j+1).eq.0) then
         vout(i,j)=v(i,j)
        elseif(ibc(i,j).ne.0.and.ibc(i,j-1).eq.0) then
         vout(i,j)=v(i,j+1)
        else
         vout(i,j)=(v(i-1,j)+v(i,j)+v(i,j+1)+v(i-1,j+1))/4. !jmn
        endif
88      if(i.eq.1) then
        dube=(cd(i,j))*sqrt(u(i,j)**2+vout(i,j)**2)        
        else
        dube=(cd(i,j)+cd(i-1,j))*.5*sqrt(u(i,j)**2+vout(i,j)**2)
        endif
        taus(i,j)=dube*u(i,j)
        taun(i,j)=dube*vout(i,j)
87      continue
    else
         do  i=1,ns
            do  j=1,nn
!                if(ibc(i,j).eq.0.or.hl(i,j).le.hmin) then
                if(ibc(i,j).ne.-1.or.hl(i,j).le.hmin) then
                    uout(i,j)=0.
                    vout(i,j)=0.
                    go to 8888
                endif
                if(i.eq.1) then 
                    vout(i,j)=v(i,j)
                    uout(i,j) = u(i,j)
                    go to 8888
                endif

                if(j.eq.nn) then
                vout(i,j)=0.0
                elseif(j.eq.1) then
                vout(i,j) = 0.0
!                elseif(ibc(i,j).ne.0.and.ibc(i,j+1).eq.0) then
!                vout(i,j)=v(i,j)
!                elseif(ibc(i,j).ne.0.and.ibc(i,j-1).eq.0) then
!                vout(i,j)=v(i,j+1)
                else
                vout(i,j) = (v(i,j)+v(i,j-1))/2.
                uout(i,j) = (u(i,j)+u(i-1,j))/2. 
                endif
                
        8888     if(i.eq.1) then
                    dube=(cd(i,j))*sqrt(uout(i,j)**2+vout(i,j)**2)        
                else
                    dube=(cd(i,j)+cd(i-1,j))*.5*sqrt(uout(i,j)**2+vout(i,j)**2)
                endif
                taus(i,j)=dube*uout(i,j)
                taun(i,j)=dube*vout(i,j)
            enddo
        enddo
    endif
!c        call vert
	  if(CALCQUASI3D) then
		call vert
		
	  endif
	  if(CALCCSED) then
	    if(TRANSEQTYPE == 2) then	   
	     call csed_DT(vardt,nct,nsteps, tmpvardt)
	    else
	     call csed(vardt,nct,nsteps, tmpvardt)
	    endif
	  else
	    call StressDiv()
	  endif
	  IF(CALCSEDAUTO) THEN
	    vardt = tmpvardt
	  ENDIF

4000 format(6f10.2)
            if(nct == 0) then
                solIndex = solIndex+1
            	CALL Calc_Area(x, y, xo, yo, nm, dn, harea)
            	CALL CG_IRIC_WRITE_SOL_TIME_F(tottime, ier)
                CALL Write_CGNS2(tottime, q)
!                IF(CALCQUASI3D.and.IO_3DOUTPUT) THEN
!                    CAll Write_CGNS3D_Grid()
!                ENDIF
                IF(CALCQUASI3D.and.IO_3DOUTPUT) THEN
                    CAll Write_CGNS3D_Grid()
!                    CAll Write_CGNS3D_SolGrid()
                    CALL Write_CGNS3D_FixedBed(solIndex, tottime, q)
                ENDIF

			else
			    if(tottime >= ptime)then
			        solIndex = solIndex+1
			    	CALL Calc_Area(x, y, xo, yo, nm, dn, harea)
                    CALL CG_IRIC_WRITE_SOL_TIME_F(tottime, ier)			    	
                    CALL Write_CGNS2(tottime, q)
                    IF(CALCQUASI3D.and.IO_3DOUTPUT) THEN
!                        CAll Write_CGNS3D_SolGrid()
                        CAll Write_CGNS3D_FixedBed(solIndex, tottime, q)
!                    ELSEIF(CALCQUASI3D.and.IO_3DOUTPUT.and.CALCCSED) THEN
!!                        CALL Write_CGNS3D_MoveableBed(tottime, q)
                    ENDIF
			        ptime = ptime+(iplinc*dt)
			    endif
			endif  

5000    ENDDO
 !       close(11)
805     format(2f11.2,6f8.2)		
!        write(3,3) runid 
!        write(3,*) taus,taun
!        write(3,*) u,vout
!        write(3,*) e
    CALL dealloc_all()
    IF(FID > 0) THEN
        CALL cg_close_f(FID, IER)
    ENDIF
 
 !	close(2)
!	close(4)
	close(9)
!        write(3,*) uz,vz

	END SUBROUTINE STAGR4
	
	SUBROUTINE dealloc_all()
	IMPLICIT NONE
		CALL dealloc_working()
	    CALL dealloc_common2d()
        if(vbc) then
            CALL DEALLOC_VELBC()
        endif
	   
	    CALL dealloc_init2d()
	    CALL dealloc_roughness()
	    CALL dealloc_TimeSeries()
	    CALL dealloc_RatingCurves()
    	CALL dealloc_TSNames()
		
	
            if(CALCQUASI3D) then
                CALL dealloc_common3d()
                if(TRANSEQTYPE) THEN
                    CAll dealloc_csed3d_dt()
                ENDIF
            endif
!	    if(CALCCSED == 1) then
!	        call dealloc_csed()
!	        call dealloc_csed_DT()
!	    endif 
        if(TRANSEQTYPE == 2) then	   
         call dealloc_csed_DT()
        else
         call dealloc_csed()
        endif

	END SUBROUTINE
	
	SUBROUTINE updateIBC()
	IMPLICIT NONE
	INTEGER :: i,j
	tibc = ibc
	    do i = 1,ns
	        do j = 1,nn
	        if(ibc(i,j).gt.0) then
	            tibc(i,j) = -1
!	        else if(j == 1) then
!!	            ibc(i,j) = 1
!	        else if(j == nn) then
!!	            ibc(i,j) = 2
	        endif
	        
	        enddo
	    enddo
	    DO i = 2,ns-1
	        Do j = 2, nn-1
	        if(ibc(i,j).ne.0.and.ibc(i-1,j).eq.0) then
	            tibc(i,j)=4 
!	        elseif(ibc(i,j).eq.-1.and.ibc(i,j-1).eq.0)then
!	            tibc(i,j)=1
!	        elseif(ibc(i,j).eq.-1.and.ibc(i,j+1).eq.0)then
!	            tibc(i,j)=2
	        endif
!	        if(ibc(i,j).eq.-1.and.ibc(i+1,j).eq.0) then
!	            tibc(i+1,j)=6
!	        endif
	        ENDDO
	    ENDDO
	    !Enforce Lateral boundary nodes
	    if(FLUMEBNDRY) then
	        Do i = 1, ns
	            tibc(i,1) = 1;
	            tibc(i,nn) = 2;
	        ENDDO
	    else
	        Do i = 1, ns
                tibc(i,1) = 0
	            u(i,1) = 0;
	            v(i,1) = 0;
    	        
                tibc(i,nn) = 0
	            u(i,nn) = 0;
	            v(i,nn) = 0;
	        ENDDO
	   endif
	    
	    ibc = tibc
	    	    
	END SUBROUTINE updateIBC
	
    SUBROUTINE RESETGRIDEXTENSION()
    IMPLICIT NONE

    INTEGER :: i,j
    REAL :: mag
    	do i=ns2+1,ns2+nsext
    	    do j = 1,nn
	            eta(i,j)=eta2(ns2,j) - ((i-ns2)*ds*nsextslope)
	            !mineta(i,j)=mineta2(ns2,j) - ((i-ns2)*ds*nsextslope)
	            if(ibc(i,j).ne.0) then
	                mag = sqrt(u(ns2,j)**2.+v(ns2,j)**2.)
                    v(i,j) = v(i,j) - (i*(-v(ns2,j)/(ns2+nsext)))
                    u(i,j) = sqrt(mag**2. - v(i,j)**2.)
                endif
            enddo
        enddo
    END SUBROUTINE RESETGRIDEXTENSION

!    SUBROUTINE INITHOTSTART()
!    IMPLICIT NONE
!     !DEC$ ATTRIBUTES REFERENCE, C, VARYING :: CG_GOTO_F
!!DEC$ ATTRIBUTES REFERENCE, C, VARYING :: CG_ARRAY_READ_F
!
!INCLUDE "cgnslib_f.h"
!   INTEGER :: i,j,k
!	CHARACTER(250) :: FILENAME, ZONENAME, BASENAME, USERNAME, FIELDNAME
!    INTEGER :: datatype
!	INTEGER :: FID, BID, ZID, IER, iret
!	INTEGER :: NUSER_DATA, USER_ITER, NZONES, ZONES_ITER, SOLN_ITER
!	INTEGER :: FIELD_ITER, nfields
!	INTEGER :: NBASES, BASES_ITER
!	INTEGER :: CELLDIM, PHYSDIM
!	INTEGER :: tnn, tns
!	INTEGER, DIMENSION(2) :: irmin2, irmax2
!	INTEGER, DIMENSION(3,3) :: isize
!	INTEGER, DIMENSION(3) :: irmin, irmax
!	INTEGER, DIMENSION(6) :: irinddata
!	CHARACTER(LEN = 250) :: name
!	INTEGER :: namelen
!	INTEGER :: count
!	REAL, ALLOCATABLE, DIMENSION(:) :: rtemp
!	INTEGER, ALLOCATABLE, DIMENSION(:) :: itemp
!	
!	!HOTSTART VARS
!	INTEGER :: nsols
!!		namelen = len(InputFile)
!!		name = TRIM(InputFile(1:namelen-3)//'cgns')
!	CALL cg_open_f(InputFile, MODE_READ, FID, IER)
!		IF(IER .NE. 0) THEN
!            call cg_error_print_f()			
!        ENDIF
!	BID = 1
!	CALL cg_nbases_f(FID, NBASES, IER)
!	DO BASES_ITER = 1, NBASES
!		CALL cg_base_read_f(FID, BASES_ITER, BASENAME, CELLDIM, PHYSDIM, IER)
!		
!		SELECT CASE(TRIM(BASENAME))
!		CASE('MD_SWMS_1D')
!
!		CASE('MD_SWMS_2D')
!		CALL cg_nzones_f(FID, BASES_ITER, NZONES, IER)
!		DO ZONES_ITER = 1, NZONES
!			CALL cg_zone_read_f(FID, BASES_ITER, ZONES_ITER, ZONENAME, isize, IER)
!                tns = isize(1,1)
!                tnn = isize(2,1)
!                irmin2(1) = 1
!                irmin2(2) = 1
!                irmax2(1)=isize(1,1)
!                irmax2(2)=isize(2,1)
!                ALLOCATE(rtemp(tns*tnn), STAT = IER)
!                ALLOCATE(itemp(tns*tnn), STAT = IER)
!			SELECT CASE(TRIM(ZONENAME))
!			CASE('CurvilinearOrthogonal')
!			    CALL cg_nsols_f(FID,BASES_ITER,ZONES_ITER,nsols,IER)
!			    DO SOLN_ITER = 1,nsols
!                IF(SOLN_ITER == SolnIndex) THEN
!                    Call cg_nfields_f(FID, BASES_ITER, ZONES_ITER, SOLN_ITER, nfields, IER)
!                    DO FIELD_ITER = 1,nfields
!                        CALL cg_field_info_f(FID, BASES_ITER, ZONES_ITER, SOLN_ITER, FIELD_ITER, datatype, FIELDNAME , IER)
!                        SELECT CASE(TRIM(FIELDNAME))
!                            CASE('Sand_Fraction')
!                                CALL cg_field_read_f(FID, BASES_ITER, ZONES_ITER, SOLN_ITER,'Sand_Fraction',RealSingle, &
!                                                    irmin2,irmax2,rtemp,IER)
!                                 DO j=1,nn
!				                    DO i= 1,ns2
!				                    	count = i + (j-1)*ns2 
!								        Fracs(i,j) = rtemp(count) 
!								    ENDDO
!								 ENDDO
!
!                            CASE('Sand_Depth')
!                                CALL cg_field_read_f(FID, BASES_ITER, ZONES_ITER, SOLN_ITER,'Sand_Depth',RealSingle, &
!                                                    irmin2,irmax2,rtemp,IER)
!                                 DO j=1,nn
!				                    DO i= 1,ns2
!				                    	count = i + (j-1)*ns2 
!                                        hfine(i,j) = rtemp(count)*100.
!   								    ENDDO
!								 ENDDO
!
!                            CASE('VelocityInitS')
!                                CALL cg_field_read_f(FID, BASES_ITER, ZONES_ITER, SOLN_ITER,'VelocityInitS',RealSingle, &
!                                                    irmin2,irmax2,rtemp,IER)
!                                 DO j=1,nn
!				                    DO i= 1,ns2
!				                    	count = i + (j-1)*ns2 
!				                        iu(i,j) = rtemp(count)*100.
!				                    ENDDO
!				                 ENDDO
!                            CASE('VelocityInitN')
!                                CALL cg_field_read_f(FID, BASES_ITER, ZONES_ITER, SOLN_ITER,'VelocityInitN',RealSingle, &
!                                                        irmin2,irmax2,rtemp,IER)
!                                 DO j=1,nn
!				                    DO i= 1,ns2
!				                    	count = i + (j-1)*ns2 
!				                        iv(i,j) = rtemp(count)*100.
!				                    ENDDO
!				                 ENDDO
!                            CASE('WaterSurfaceElevation')
!                                CALL cg_field_read_f(FID, BASES_ITER, ZONES_ITER, SOLN_ITER,'WaterSurfaceElevation',RealSingle, &
!                                                    irmin2,irmax2,rtemp,IER)
!                                 DO j=1,nn
!				                    DO i= 1,ns2
!				                    	count = i + (j-1)*ns2 
!				                        iwse(i,j) = rtemp(count)*100.
!				                    ENDDO
!				                 ENDDO
!                            CASE('Elevation')
!                                CALL cg_field_read_f(FID, BASES_ITER, ZONES_ITER, SOLN_ITER,'Elevation',RealSingle, &
!                                                    irmin2,irmax2,rtemp,IER)
!                                 DO j=1,nn
!				                    DO i= 1,ns2
!				                    	count = i + (j-1)*ns2 
!				                        ie(i,j) = rtemp(count)*100.
!				                    ENDDO
!				                 ENDDO
!                            CASE('IBC')
!                                CALL cg_field_read_f(FID, BASES_ITER, ZONES_ITER, SOLN_ITER,'IBC',INTEGER, &
!                                                    irmin2,irmax2,itemp,IER)
!                                 DO j=1,nn
!				                    DO i= 1,ns2
!				                    	count = i + (j-1)*ns2 
!				                        iibc(i,j) = itemp(count)
!				                    ENDDO
!				                 ENDDO
!                        END SELECT
!                    ENDDO !field_iter
!                ENDIF
!			    
!			    
!			    ENDDO
!			
!			END SELECT
!		ENDDO	!ZONES
!		END SELECT
!	ENDDO !BASES
!	CALL cg_close_f(FID, IER)
!
!    
!    
!    END SUBROUTINE INITHOTSTART
    

       
    SUBROUTINE Calc_Area(x, y, xo, yo, nm, dn, area)
         IMPLICIT NONE
       REAL, INTENT (OUT), DIMENSION (:,:) :: area
       REAL, INTENT(INOUT), DIMENSION(:,:):: x,y
       REAL, INTENT (IN), DIMENSION (:) :: xo, yo
       INTEGER, INTENT(IN) :: nm
       REAL, INTENT(IN) :: dn
       
        INTEGER :: I, J
        REAL :: x1, x2, x3, x4
        REAL :: y1, y2, y3, y4
        REAL :: area1, area2
        
        REAL :: rsin, rcos, ux, uy, uxnew, uynew
        
        INTEGER :: icount
         do i=1,ns2
			rcos=cos(phirotation(i))
            rsin=sin(phirotation(i))

			do j=1,nn
!			    ux=(u(i,j)*rcos-vout(i,j)*rsin)
!                uy=(u(i,j)*rsin+vout(i,j)*rcos)
!                uxnew=ux*fcos-uy*fsin
!                uynew=ux*fsin+uy*fcos
!!                        write(15,*) i,j,xo(i), nm, dn, phirotation(i), u(i,j), vout(i,j), fcos, fsin, x(i,j), y(i,j), xo(i), yo(i)
			    x(i,j)=xo(i)+(nm-j)*dn*sin(phirotation(i))
                y(i,j)=yo(i)+(j-nm)*dn*cos(phirotation(i))
!!     write(11,4000)x(i,j),y(i,j),eta(i,j),hl(i,j),uxnew,uynew
!4000 format(6f10.2)
			        enddo
			    enddo
! Loop through the mesh and calculate areas for each point except the edges of the mesh

        do i=2,ns-1
        do j=2,nn-1
           x1=(x(i-1,j-1)+x(i-1,j)+x(i,j)+x(i,j-1))/4.
           x2=(x(i-1,j)+x(i-1,j+1)+x(i,j+1)+x(i,j))/4.
           x3=(x(i,j)+x(i,j+1)+x(i+1,j+1)+x(i+1,j))/4.
           x4=(x(i,j-1)+x(i,j)+x(i+1,j)+x(i+1,j-1))/4.
           y1=(y(i-1,j-1)+y(i-1,j)+y(i,j)+y(i,j-1))/4.
           y2=(y(i-1,j)+y(i-1,j+1)+y(i,j+1)+y(i,j))/4.
           y3=(y(i,j)+y(i,j+1)+y(i+1,j+1)+y(i+1,j))/4.
           y4=(y(i,j-1)+y(i,j)+y(i+1,j)+y(i+1,j-1))/4.
!c          using area of a triangle algorithm in standard mathmatical tables book page 362
!c          0.5*(x1*y2+x2*y3+x3*y1-y1*x2-y2*x3-y3*x1)
!c          Where in the first case x1=>x1,x2=>x4,x3=>x2
!c          Where in the second case x1=>x3,x2=>x2,x3=>x4
           area1=0.5*(x1*y4+x4*y2+x2*Y1-y1*x4-y4*x2-y2*x1)
           area2=0.5*(x3*y2+x2*y4+x4*y3-y3*x2-y2*x4-y4*x3)
           area(i,j)=area1+area2
!c           WRITE(9,*) i,j,x1,x2,x3,x4,y1,y2,y3,y4,area1,area2  ! debug
        end do
        end do

!C Loop through the mesh and estimate the cell areas for the edge of the mesh
!C and write the node number,x coord, y coord, and area out to a file

        icount=0
! c       WRITE(9,*)'node,x,y,area'

        do i=1,ns
        do j=1,nn
        icount=icount+1
!c       deal with left edge of the mesh (take area from interior)
        IF(j.eq.1)then
          IF(i.eq.1)then
           area(i,j)=area(2,2)
          ELSEIF(i.eq.ns)then
           area(i,j)=area(ns-1,2)
          else
           area(i,j)=area(i,2)
          endif
!c       deal with right edge of the mesh (take area from interior)
        ELSEIF(j.eq.nn)then
          IF(i.eq.1)then
           area(i,j)=area(2,nn-1)
          ELSEIF(i.eq.ns)then
           area(i,j)=area(ns-1,nn-1)
          else
           area(i,j)=area(i,nn-1)
          endif
        ENDIF
!c       deal with the interior ends of the mesh (top and bottom of the mesh) areas
        IF(i.eq.1)then
         IF(j.gt.1.and.j.lt.nn)then
           area(i,j)=area(2,j)
         endif
        ELSEIF(i.eq.ns)then
         IF(j.gt.1.and.j.lt.nn)then
           area(i,j)=area(ns-1,j)
         endif
        ENDIF

!!c       write out the x coordinates ,y coordinates, and area values
!        x(i,j)=x(i,j)+offsetx
!        y(i,j)=y(i,j)+offsety
!        WRITE(9,210) icount,x(i,j),y(i,j),area(i,j)

        end do
        end do
    END SUBROUTINE Calc_Area
    
    


	END MODULE RivStagr4Mod_jmn
