MODULE RivCalcInitCond
USE RivVarMod
USE RivVarWMod
USE GridCoord
USE CalcCond
USE NetPreisMain

IMPLICIT NONE

CONTAINS

SUBROUTINE calcWSInitCond()

REAL(KIND=mp) :: tslope
INTEGER :: i,j, ier, count
REAL(KIND=mp), DIMENSION(:), ALLOCATABLE :: tx, ty, ttopo
REAL(KIND=mp) :: rcos, rsin, ux, uy, uxnew, uynew
REAL(KIND=mp) :: tmpwselev, tmpq

!Calc Init Water Surface Elevation

IF(wsType.eq.0) THEN !Constant slope using upper water surface elevation
    
    tslope = (wsupelev-wselev)/(ns*scals)
    DO i = ns,1,-1
        hav(i) = wselev +((scals*(ns-i))*tslope)
    ENDDO
    DO I = 1,NS
        DO J = 1,NN
            e(i,j) = hav(i)
        ENDDO
    ENDDO

ELSE IF(wsType.eq.1) THEN !Constant slope using given water surface slope
    
    DO i = ns,1,-1
        hav(i) = wselev +((scals*(ns-i))*wsslope)
    ENDDO
    DO I = 1,NS
        DO J = 1,NN
            e(i,j) = hav(i)
        ENDDO
    ENDDO

ELSE IF(wsType.eq.2) THEN

    ALLOCATE(tx(ns*nn), ty(ns*nn), ttopo(ns*nn), STAT = ier)
		do j=1,nn
			do i=1,ns2+nsext
			    rcos=cos(phirotation(i))
                rsin=sin(phirotation(i))
!                count = i + (j-1)*(ns2+nsext)
				count = ((i-1)*nn)+j
				ux = xo(i)+(dn)*(nm-j)*sin(phirotation(i))
				uy = yo(i)-(dn)*(nm-j)*cos(phirotation(i))
                uxnew=ux*fcos-uy*fsin
                uynew=ux*fsin+uy*fcos
                tx(count) = (uxnew+xshift)/100.
                ty(count) = (uynew+yshift)/100.
                If(i.le.ns2) THEN
                    ttopo(count) = eta(i,j)/100.
                ELSE
                    ttopo(count) = eta(ns2,j)/100. - ((i-ns2)*ds*nsextslope)
                ENDIF
			enddo
		enddo

!    DO i = 1,ns
!        DO j = 1,nn
!            count = ((i-1))*nn +j
!            IF(i.le.ns2) THEN
!                ttopo(count) = eta(i,j)/100.
!                tx(count) = x(i,j)/100.
!                ty(count) = y(i,j)/100.
!            ELSE
!                ttopo(count) = eta(i,j)/100.
!                tx(count) = x(i,j)/100.
!                ty(count) = y(i,j)/100.
!              
!        ENDDO
!    ENDDO
    if(varDischType == 1) then !Discharge Time Series
        CALL getInterpTimeSeriesValue(1, VarDischStartTime, tmpq)
    else
        tmpq = q/1e6
    endif
    
    if(varStageType == 1) then !Stage Time Series
        CALL getInterpRatingCurveValue(1, VarDischStartTime, tmpwselev)
        tmpwselev = tmpwselev/100.
!        newStage = (newStage-elevoffset) * 100.
    else if (varStageType == 2) then !Stage Rating Curve
        Call getInterpRatingCurveValue(1, tmpq, tmpwselev)
        tmpwselev = tmpwselev/100.
!        newStage = (newStage-elevoffset) * 100.
    else
        tmpwselev = wselev/100 
    endif
!    CALL NETPREIS2(ns, nn, ttopo, tx, ty, ONEDCD, tmpq/1e6, tmpwselev/100., hav )
    CALL NETPREIS2(ns, nn, ttopo, tx, ty, ONEDCD, tmpq, tmpwselev, hav )
    hav = hav*100
    DO I = 1,NS
        DO J = 1,NN
            e(i,j) = hav(i)
        ENDDO
    ENDDO
    
ELSE IF(wstype.eq.3) THEN
!    CALL INITHOTSTART2(CGNSHSFile)
    DO I = 1,NS
        DO J = 1,NN
!		    if(i.le.ns2) then
!		        ibc2(i,j) = iibc(i,j)
!		        eta2(i,j) = ie(i,j)
!		        e2(i,j) = iwse(i,j)
!                u2(i,j) = iu(i,j)
!                v2(i,j) = iv(i,j)
!		    endif

		    if(i.le.ns2) then
			    u(i,j) = iu(i,j)
			    v(i,j) = iv(i,j)
			    e(i,j) = iwse(i,j)
			    ibc(i,j) = iibc(i,j)
			    eta(i,j) = ie(i,j)
		    else 
			    u(i,j) = u(ns2,j)
			    v(i,j) = v(ns2,j)
			    e(i,j) = e(ns2,j)
			    ibc(i,j) = ibc(ns2,j)
	            eta(i,j) = eta(ns2,j) - ((i-ns2)*ds*nsextslope)
	            fracs(i,j) = fracs(ns2,j)
	            hfine(i,j) = hfine(ns2,j)
	            
			    
		    endif
		    hl(i,j) = e(i,j)-eta(i,j)
		    IF(j == nn/2+1) THEN
		        hav(i) = e(i,j)
		    ENDIF
        if(ibc(i,j).eq.0.or.hl(i,j).le.hmin) then
            if(dryType.eq.0) then !turn node off
            if(ibc(i,j).ne.0) then
                u(i,j) = 0
                v(i,j) = 0
            endif
               ibc(i,j)=0
            endif
           hl(i,j)=hmin
        endif
        ENDDO
    ENDDO       
ENDIF

DO i = 1,ns
    DO j = 1,nn
        IF(j == 1) THEN
            ibc(i,j) = 0
        ELSE IF(j == nn) THEN
            ibc(i,j) = 0
        
        ELSE IF(e(i,j).gt.eta(i,j)) THEN
            ibc(i,j) = -1
        ELSE
            ibc(i,j) = 0
        ENDIF
    ENDDO
ENDDO

END SUBROUTINE

!SUBROUTINE INITHOTSTART2(InputFile)
!    USE iriclibf
!    IMPLICIT NONE
!!    INCLUDE "cgnswin_f.h"
!    INCLUDE "cgnslib_f.h"
!
!    CHARACTER(*), INTENT(IN) ::InputFile
!    
!    INTEGER :: FID
!    INTEGER :: i,j,k, IER
!    INTEGER :: NX, NY
!	REAL, ALLOCATABLE, DIMENSION(:) :: tmpreal4
!	INTEGER, ALLOCATABLE, DIMENSION(:) :: itemp
!	INTEGER :: count, countji
!	
!	!HOTSTART VARS
!	INTEGER :: nsols
!	
!    CALL cg_open_f(InputFile, MODE_MODIFY, FID, IER)
!    IF(IER .NE. 0) THEN
!        call cg_error_print_f()			
!    ENDIF
!
!    CALL CG_IRIC_GOTOSOL_F(FID, SolnIndex, NX, NY, IER)
!    ALLOCATE(tmpreal4(nx*ny), STAT = ier)
!    ALLOCATE(itemp(nx*ny), STAT = ier)
!    
!    
!    CALL cg_iRIC_Read_SolRealNode(SolnIndex, 'Elevation', tmpreal4, IER)
!!    DO I = 1, nx*ny
!!        IF( tmpreal4(i) < elevoffset) THEN
!!            elevoffset = tmpreal4(i)
!!        ENDIF
!!    ENDDO
!    DO I= 1,NX
!        DO J=1,NY
!            COUNT = ((I-1)*NY)+J
!            countji = ((j-1)*NX)+i
!            ie(I,J) = (tmpreal4(countji) - elevoffset)*100.
!        ENDDO
!    ENDDO
!    
!    CALL cg_iRIC_Read_SolRealNode(SolnIndex, 'roughness', tmpreal4, IER)
!    IF(ier.eq.0) THEN
!        DO I= 1,NX
!            DO J=1,NY
!                COUNT = ((I-1)*NY)+J
!                countji = ((j-1)*NX)+i
!			    cd2(i,j) = tmpreal4(countji)
!			    znaught2(i,j) = tmpreal4(countji)*100.
!            ENDDO
!        ENDDO
!    ENDIF
!    cdv2 = 0.0
!    CALL cg_iRIC_Read_SolRealNode(SolnIndex, 'vegroughness', tmpreal4, IER)
!    IF(ier.eq.0) THEN
!        DO I= 1,NX
!            DO J=1,NY
!                COUNT = ((I-1)*NY)+J
!                countji = ((j-1)*NX)+i
!			    cdv2(i,j) = tmpreal4(countji)
!            ENDDO
!        ENDDO
!    ENDIF
!
!    CALL cg_iRIC_Read_SolRealNode(SolnIndex, 'sandfraction', tmpreal4, IER)
!    IF(ier.eq.0) THEN
!        DO I= 1,NX
!            DO J=1,NY
!                COUNT = ((I-1)*NY)+J
!                countji = ((j-1)*NX)+i
!			    Fracs(i,j) = tmpreal4(countji)
!            ENDDO
!        ENDDO
!    ENDIF
!    CALL cg_iRIC_Read_SolRealNode(SolnIndex, 'sanddepth', tmpreal4, IER)
!    IF(ier.eq.0) THEN
!        DO I= 1,NX
!            DO J=1,NY
!                COUNT = ((I-1)*NY)+J
!                countji = ((j-1)*NX)+i
!			    hfine(i,j) = tmpreal4(countji)*100. !Convert to cm
!            ENDDO
!        ENDDO
!    ENDIF
!    
!    CALL cg_iRIC_Read_SolRealNode(SolnIndex, 'VelocityInitS', tmpreal4, IER)
!    IF(IER.ne.0) THEN
!        STOP 'No Initial Velocity in hotstart solution'
!        pause
!    ENDIF
!    DO I= 1,NX
!        DO J=1,NY
!            COUNT = ((I-1)*NY)+J
!            countji = ((j-1)*NX)+i
!			iu(i,j) = tmpreal4(countji)*100. !Convert to cm
!        ENDDO
!    ENDDO
!
!    CALL cg_iRIC_Read_SolRealNode(SolnIndex, 'VelocityInitN', tmpreal4, IER)
!    DO I= 1,NX
!        DO J=1,NY
!            COUNT = ((I-1)*NY)+J
!            countji = ((j-1)*NX)+i
!			iv(i,j) = tmpreal4(countji)*100. !Convert to cm
!        ENDDO
!    ENDDO
!
!    CALL cg_iRIC_Read_SolRealNode(SolnIndex, 'WaterSurfaceElevation', tmpreal4, IER)
!    DO I= 1,NX
!        DO J=1,NY
!            COUNT = ((I-1)*NY)+J
!            countji = ((j-1)*NX)+i
!			iwse(i,j) = (tmpreal4(countji) - elevoffset)*100. !Convert to cm
!			if(j == ny/2+1) then
!			    hav(i) = (tmpreal4(countji) - elevoffset)*100 !Convert to cm
!			ENDIF
!        ENDDO
!    ENDDO
!    IF(nx < ns) THEN
!        DO I = NX+1,NS
!            hav(i) = hav(nx)
!        ENDDO
!    ENDIF
!    
!    CALL cg_iRIC_Read_SolIntegerNode(SolnIndex, 'IBC', itemp, IER)
!    DO I= 1,NX
!        DO J=1,NY
!            COUNT = ((I-1)*NY)+J
!            countji = ((j-1)*NX)+i
!			iibc(i,j) = itemp(countji) 
!        ENDDO
!    ENDDO
!	CALL cg_close_f(FID, IER)
!
!    END SUBROUTINE INITHOTSTART2
    
    SUBROUTINE INITARRAYS()
    IMPLICIT NONE
    INTEGER :: i,j
    REAL :: rcos, rsin, uu, vv
!    open(4,FILE='RadiusOfCurvature.txt')
    	IF (cdtype == 0) THEN	!Constant CD	
			DO i = 1, ns2
				DO j = 1, nn
					znaught2(i,j) = constcd*100.
					cd2(i,j) = constcd 
					cdv2(i,j) = 0
				END DO
			END DO
		ENDIF

    do 25 i=1,ns
		rcos = cos(phirotation(i))
		rsin = sin(phirotation(i))

!    lengthen grid to eliminate eddies at downstream end
        if(i.le.ns2) then
	        r(i)=r2(i)
	        w(i)=w2(i)
!	        hav(i)=hav2(i)
!	        xo(i)=xo2(i)
!	        yo(i)=yo2(i)
	    else
		    r(i) = r2(i)
!		   		r(i) = r2(ns2)+ (i-ns2)*dr

	        w(i)=w2(ns2)
!	        hav(i)=hav2(ns2)
	    endif

!
        do 25 j=1,nn
!   lengthen grid to elimate eddies at downstream end
            if(i.le.ns2) then
	            eta(i,j)=eta2(i,j)
	            mineta(i,j) = mineta2(i,j)
!	            ibc(i,j)=ibc2(i,j)
!	            ribc(i,j) = ribc2(i,j)
	            cd(i,j)=cd2(i,j)
	            cdv(i,j)=cdv2(i,j)
	            totcd(i,j) = cd(i,j)+cdv(i,j)
		        znaught(i,j)=znaught2(i,j)
	        else
	            eta(i,j)=eta2(ns2,j) - ((i-ns2)*ds*nsextslope)
	            mineta(i,j)=mineta2(ns2,j) - ((i-ns2)*ds*nsextslope)
!	            ibc(i,j)=ibc2(ns2,j)
!	            ribc(i,j) = ribc2(ns2, j)
	            cd(i,j)=cd2(ns2,j)
	            cdv(i,j)=cdv2(ns2,j)
	            totcd(i,j) = cd(i,j)+cdv(i,j)
		        znaught(i,j)=znaught2(ns2,j)
		        fracs(i,j) = fracs(ns2,j)
		        hfine(i,j) = hfine(ns2,j)
	        endif
		    rn(i,j)=1.-(j-nm)*dn/r(i)
!		    if(j.eq.1)then
!    	        write(4,*) i,r(i)/100.
!    	    endif
		    dude(i,j)=0.
		    dvde(i,j)=0.
		    dsu(i,j)=ds*rn(i,j)
		    dsv(i,j)=ds*rn(i,j)
		    dse(i,j)=ds*rn(i,j)
		    dsq(i,j)=ds*rn(i,j)
		    dnu(i,j)=dn
		    dnv(i,j)=dn
		    dne(i,j)=dn
		    dnq(i,j)=dn
25      continue    
!    close(4)
    END SUBROUTINE INITARRAYS


END MODULE