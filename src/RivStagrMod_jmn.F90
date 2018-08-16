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
    USE fm_helpers
    IMPLICIT NONE

    CONTAINS

    SUBROUTINE STAGR4(STR_IN )
    IMPLICIT NONE
    INCLUDE "iriclib_f.h"
    CHARACTER(LEN=*), INTENT(IN) :: STR_IN

    !character*20 runid
    INTEGER :: solIndex
    INTEGER :: n, n2
    REAL(KIND=mp) :: ustar
    REAL(KIND=mp) :: dinc
    REAL(KIND=mp) :: area
    REAL(KIND=mp) :: ba, bb, bc, be, bf, bg, bot, bot2
    REAL(KIND=mp) :: be1, be2, bf1, bf2
    REAL(KIND=mp) :: va, ta, tb, tc, td, te, tf
    REAL(KIND=mp) :: ua, uu, vv
    REAL(KIND=mp) :: dube
    REAL(KIND=mp) :: sumqpdiff
    INTEGER, PARAMETER :: num_string = 4
    INTEGER :: i, j, jswi, j1, j2, j3
    REAL(KIND=mp) :: plinc
    LOGICAL :: lwet, doagain

    !	output vars
    REAL(KIND=mp) :: rcos, rsin, ux, uy, uxnew, uynew, qx, qy
    REAL(KIND=mp) :: dhav2

    REAL(KIND=mp) :: rnr1, rnrm1, rnrm1j, rnrm1ij !Changed as per CD
    REAL(KIND=mp) :: dsvpimj, dsvmimj, ekapimj, ekamimj
    INTEGER :: ip, im, jp, jm

    REAL(KIND=mp) :: hlrewe
    REAL(KIND=mp) :: dr



    REAL(KIND=mp) :: nodechange
    !   tmpVars
    REAL(KIND=mp) :: tmpLEV, tmpvardt
    INTEGER :: ier

    CHARACTER(20) text
    INTEGER(4) iret
    INTEGER :: ibccount

    !    OPEN(11,FILE='topser')
    nwetnodes = 0
    errorcode = 0
    !pi=acos(-1.)
    g=980.
    rho=1.
    vkc=0.40
    solIndex = 0
    CALL WELCOME
    CALL read_iRIC_CGNS(STR_IN)




3   format(A20)
    CALL Calc_Area(x, y, xo, yo, nm, dn, harea)

    !pause
65  format(25f6.2)
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

    if(.not.HotStart) then
        CALL einit(e, hl, eta, ibc, w, hav)
        call uinit(u,v,hav,w,eta,q,ibc)
        if(errorcode < 0) then
            !call write_error_CGNS(STR_IN)
            CALL dealloc_common2d()
            if(vbc == 1) then
                CALL DEALLOC_VELBC()
            endif
            if(CALCQUASI3D) then
                CALL dealloc_common3d()
                if(TRANSEQTYPE == 1) THEN
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

    nsteps = (VarDischEndTime-VarDischStartTime)/fmdt

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
    ptime = totTime+iplinc*fmdt

4   do while(totTime <= VarDischEndTime)
        nct = nct+1

        if(nct > 0.and.solType == 0) then
            Exit
        endif
        if(nct > 0.and.VarDischEndTime.eq.VarDischStartTime) Then
            EXIT
        ENDIF

        if(nct == 0) then
            vardt = fmdt
        else
            if(vardt > fmdt) then
                vardt = fmdt
            endif
        endif

        if(CALCCSED.eqv..FALSE.) then
            vardt = fmdt
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

        !!  IF GRID EXTENSION THEN RESET EXTENSION TOPOGRAPHY !!
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
        Write(*,*) 'Variabledt', tmpvardt, 'DT', fmdt
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

        if(itm /= 0) then
            ITER_LOOP: DO iter = 1,itm
                !       Check if user cancelled simulation and if so exit
                CALL iric_check_cancel_f(errorcode)
                IF(errorcode.eq.1) THEN
                    CALL dealloc_all()
                    return
                ENDIF

                IF(iwetdry.and.MOD(iter,imod).eq.0.and.varDischType == 0.and.iter.lt.hiterstop.and.nct.lt.1)then
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


                do i=2,ns
                    do jswi=1,2
                        if(jswi.eq.1) then
                            j1=1
                            j2=nn
                            j3=1
                        else
                            j1=nn
                            j2=1
                            j3=-1
                        endif
                        do j=j1,j2,j3
                            if(ibc(i,j).eq.2) then
                                va=(v(i-1,j)+v(i,j))/4.
                            else if(j.lt.nn) then
                                va=(v(i-1,j)+v(i,j)+v(i,j+1)+v(i-1,j+1))/4.
                            endif
                            if(ibc(i,j).eq.4.or.ibc(i,j).eq.0) then
                                dude(i,j)=0.
                                u(i,j)=0.
                                CYCLE
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
                        enddo
                    enddo
                enddo

                do i=2,ns
                    do jswi=1,2
                        if(jswi.eq.1) then
                            j1=1
                            j2=nn
                            j3=1
                        else
                            j1=nn
                            j2=1
                            j3=-1
                        endif
                        do j=j1,j2,j3
                            if(ibc(i,j).eq.0.or.ibc(i,j).eq.1) then
                                v(i,j)=0.
                                dvde(i,j)=0.
                                CYCLE ! will need to make sure this change is correct
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
                        enddo
                    enddo
                enddo

                do i=1,ns
                    do j=1,nn
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
                    enddo
                enddo

                do i=1,ns-1
                    do j=1,nn
                        if(j.eq.nn) then
                            dq(i,j)=qu(i+1,j)-qu(i,j)-qv(i,j)
                        else
                            dq(i,j)=qu(i+1,j)-qu(i,j)+qv(i,j+1)-qv(i,j)
                        endif
                    enddo
                enddo

                do i=1,ns-1
                    do j=1,nn
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
                    enddo
                enddo

                do i=ns-1,2,-1
                    n=0.
                    do j=1,nn
                        !c       if(abs(dqe(i,j)).gt..000001) then
                        if(abs(dqe(i,j)).gt..000001.and.hl(i,j).gt.hmin) then
                            n=n+1
                            if(j.ne.1) then
                                am(n)=-1.*dqnde(i,j)
                            endif
                            bm(n)=dqe(i,j)

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
                    enddo
                    if(n.ge.1) then
                        call tridag(am,bm,ccm,dm,em,n,errorcode)
                        if(errorcode.eq.-1) then
                            !call write_error_CGNS(STR_IN)
                            CALL dealloc_common2d()
                            if(vbc == 1) then
                                CALL DEALLOC_VELBC()
                            endif

                            if(CALCQUASI3D) then
                                CALL dealloc_common3d()
                                if(TRANSEQTYPE == 1) THEN
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
                    do j=1,nn
                        if(isw(j).eq.1) then
                            n2=n2+1
                            de(i,j)=em(n2)
                        else
                            de(i,j)=0.
                        endif
                    enddo
                enddo

                do j=1,nn
                    n=0
                    do i=ns-1,2,-1
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
                    enddo
                    if(n.ge.1) then
                        call tridag(am,bm,ccm,dm,em,n,errorcode)
                        if(errorcode.eq.-1) then
                            !call write_error_CGNS(STR_IN)
                            CALL dealloc_common2d()
                            if(vbc == 1) then
                                CALL DEALLOC_VELBC()
                            endif

                            if(CALCQUASI3D) then
                                CALL dealloc_common3d()
                                if(TRANSEQTYPE == 1) THEN
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
                    do i=ns-1,2,-1
                        if(isw(i).eq.1) then
                            n2=n2+1
                            de(i,j)=em(n2)
                        else
                            de(i,j)=0.
                        endif
                    enddo
                enddo

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

                do  j=1,nn
                    u(1,j)=u(1,j)*q/qp(1)
                enddo

                dea(ns)=0.
                do  i=ns-1,1,-1
                    dinc=arelax*erelax*abs(eav(i)-eav(i+1))*(1.-qp(i+1)/q)
                    dea(i)=dea(i+1)+dinc
                    do  j=1,nn
                        !! code change by rmcd !!
                        if(ibc(i,j).ne.0)then
600                         e(i,j)=e(i,j)+dea(i)
                            !                else
                            !                    e(i,j) = eav(i)+dea(i) !! added so that e(i,j) in dry nodes follows
                            !                                           !! average eav(i) during time-dependent runs
                        endif
                        !! end code change by rmcd !!
                    enddo
                enddo

                sumqpdiff = 0.
                do i=1,ns
                    sumqpdiff = sumqpdiff + (((qp(i)-q)/q*100.))**2.
                enddo

                qprms(iter) = sqrt(sumqpdiff/ns)
                !write(*,*) 'after qprms'
                if(isnan(qprms(iter))) then
                    errorcode=-1
                    !call write_error_CGNS(STR_IN)
                    !			call dealloc_all()
                    CALL dealloc_common2d()
                    if(vbc == 1) then
                        CALL DEALLOC_VELBC()
                    endif

                    if(CALCQUASI3D) then
                        CALL dealloc_common3d()
                        if(TRANSEQTYPE == 1) THEN
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
                IF(LEVType == 0) THEN
                    WRITE(*,*) 'Iteration: ', iter,' Mean Error on Discharge:', qprms(iter)
                ELSE
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

                    do  i=1,ns
                        do  j=1,nn
                            !                if(ibc(i,j).eq.0.or.hl(i,j).le.hmin) then
                            if(ibc(i,j).ne.-1.or.hl(i,j).le.hmin) then
                                u(i,j)=0.
                                vout(i,j)=0.
                                !    go to 888
                                !endif
                            elseif(i.eq.1) then
                                vout(i,j)=v(i,j)
                                !    go to 888
                                !endif

                            elseif(j.eq.nn) then
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

888                         if(i.eq.1) then
                                dube=(cd(i,j))*sqrt(u(i,j)**2+vout(i,j)**2)
                            else
                                dube=(cd(i,j)+cd(i-1,j))*.5*sqrt(u(i,j)**2+vout(i,j)**2)
                            endif
                            taus(i,j)=dube*u(i,j)
                            taun(i,j)=dube*vout(i,j)
                        enddo
                    enddo

                    solIndex = solIndex+1
                    CALL CG_IRIC_WRITE_SOL_TIME_F(tottime, ier)
                    CALL Write_CGNS2(tottime, q)

                ENDIF

            END DO ITER_LOOP

            IF(i_re_flag_o.eq.1.and.nct.eq.0) THEN
                IF(i_tmp_count <= n_rest) THEN
                    IF(iter.eq.opt_tmp(i_tmp_count)) THEN
                        tmp_file_o(i_tmp_count)=trim(tmp_pass)//tmp_file_o(i_tmp_count)  !i110419
                        open(502,file=tmp_file_o(i_tmp_count) &
                            ,status='unknown',form='unformatted')
                        !
                        write(502) tottime,solIndex,fmdt
                        write(502) ns,nn
                        !
                        write(502) ((eta(i,j),i=1,ns),j=1,nn)
                        write(502) ((u(i,j),i=1,ns),j=1,nn)
                        write(502) ((v(i,j),i=1,ns),j=1,nn)
                        write(502) ((ibc(i,j),i=1,ns),j=1,nn)
                        write(502) ((e(i,j),i=1,ns),j=1,nn)
                        write(502) ((hl(i,j),i=1,ns),j=1,nn)
                        write(502) (hav(i),i=1,ns)

                        close(502)

                        i_tmp_count = i_tmp_count +1
                    ENDIF
                ENDIF
            ENDIF

64          format(25f7.3)
            !66      format(A15,I6,A25,F15.5)
66          format(A, A, I5,A, F15.5)
67          format(A, A, I5,A, F15.5,A, F15.5)
68          format(A)
75          format(25f7.0)
85          format(25f6.4)
86          format(25f8.3)
84          FORMAT(3i5,f6.2,2f12.1,2g12.4)

            do i = 1,ns
                do j = 1,nn
                    if(ibc(i,j).eq.0) then
                        v(i,j) = 0.
                        u(i,j) = 0.
                    endif
                enddo
            enddo
        else

            !601 continue
            do i=1,ns
                do j=1,nn
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
88                  if(i.eq.1) then
                        dube=(cd(i,j))*sqrt(u(i,j)**2+vout(i,j)**2)
                    else
                        dube=(cd(i,j)+cd(i-1,j))*.5*sqrt(u(i,j)**2+vout(i,j)**2)
                    endif
                    taus(i,j)=dube*u(i,j)
                    taun(i,j)=dube*vout(i,j)
                enddo
            enddo

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

4000        format(6f10.2)
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
                    ptime = ptime+(iplinc*fmdt)
                endif
                IF(i_re_flag_o.eq.1.and.nct.ne.0) THEN
                    IF(i_tmp_count <= n_rest) THEN
                        IF(tottime.ge.opt_tmp(i_tmp_count)) THEN
                            tmp_file_o(i_tmp_count)=trim(tmp_pass)//tmp_file_o(i_tmp_count)  !i110419
                            open(502,file=tmp_file_o(i_tmp_count) &
                                ,status='unknown',form='unformatted')
                            !
                            write(502) tottime,solIndex,fmdt
                            write(502) ns,nn
                            !
                            write(502) ((eta(i,j),i=1,ns),j=1,nn)
                            write(502) ((u(i,j),i=1,ns),j=1,nn)
                            write(502) ((v(i,j),i=1,ns),j=1,nn)
                            write(502) ((ibc(i,j),i=1,ns),j=1,nn)
                            write(502) ((e(i,j),i=1,ns),j=1,nn)
                            write(502) ((hl(i,j),i=1,ns),j=1,nn)
                            write(502) (hav(i),i=1,ns)

                            close(502)

                            i_tmp_count = i_tmp_count +1
                        ENDIF
                    ENDIF
                ENDIF
            endif
        endif
    ENDDO
    !       close(11)
805 format(2f11.2,6f8.2)
    !        write(3,3) runid
    !        write(3,*) taus,taun
    !        write(3,*) u,vout
    !        write(3,*) e
    CALL dealloc_all()
    IF(CGNSFILEID > 0) THEN
        CALL cg_close_f(CGNSFILEID, IER)
    ENDIF

    !	close(2)
    !	close(4)
    close(9)
    !        write(3,*) uz,vz

    END SUBROUTINE STAGR4










    END MODULE RivStagr4Mod_jmn
