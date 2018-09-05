    MODULE RivCalcInitCond2
    use fm_global
    USE RivVarMod2
    USE RivVarWMod2
    USE RivVarTimeMod2
    !USE GridCoord
    USE CalcCond2
    USE NetPreisMain

    IMPLICIT NONE

    CONTAINS

    SUBROUTINE calcWSInitCond(rvo, rwvo, cco, rvto)
    implicit none
    type(rivvar), intent(inout) :: rvo
    type(riv_w_var), intent(inout) :: rwvo
    type(calccond), intent(inout) :: cco
    type(rivvartime), intent(inout) :: rvto
    REAL(KIND=mp) :: tslope
    INTEGER :: i,j, ier, count
    REAL(KIND=mp), DIMENSION(:), ALLOCATABLE :: tx, ty, ttopo
    REAL(KIND=mp) :: rcos, rsin, ux, uy, uxnew, uynew
    REAL(KIND=mp) :: tmpwselev, tmpq
    REAL(KIND=mp) :: ttime, dt2
    INTEGER :: nx2, ny2,icount
    integer :: tmpns, tmpnn
    !Calc Init Water Surface Elevation
    tmpns = rvo%ns
    tmpnn = rvo%nn
    IF(cco%wstype.eq.0) THEN !Constant slope using upper water surface elevation

        tslope = (cco%wsupelev-cco%wselev)/(tmpns*rvo%scals)
        DO i = tmpns,1,-1
            rvo%hav(i) = cco%wselev +((rvo%scals*(tmpns-i))*tslope)
        ENDDO
        DO I = 1,tmpns
            DO J = 1,tmpnn
                rvo%e(i,j) = rvo%hav(i)
            ENDDO
        ENDDO

    ELSE IF(cco%wstype.eq.1) THEN !Constant slope using given water surface slope

        DO i = tmpns,1,-1
            rvo%hav(i) = cco%wselev +((rvo%scals*(tmpns-i))*cco%wsslope)
        ENDDO
        DO I = 1,tmpns
            DO J = 1,tmpnn
                rvo%e(i,j) = rvo%hav(i)
            ENDDO
        ENDDO

    ELSE IF(cco%wstype.eq.2) THEN

        ALLOCATE(tx(tmpns*tmpnn), ty(tmpns*tmpnn), ttopo(tmpns*tmpnn), STAT = ier)
        do j=1,tmpnn
            do i=1,rvo%ns2+rvo%nsext
                rcos=cos(rvo%phirotation(i))
                rsin=sin(rvo%phirotation(i))
                !                count = i + (j-1)*(rvo%ns2+rvo%nsext)
                count = ((i-1)*tmpnn)+j
                ux = rvo%xo(i)+(rvo%dn)*(rvo%nm-j)*sin(rvo%phirotation(i))
                uy = rvo%yo(i)-(rvo%dn)*(rvo%nm-j)*cos(rvo%phirotation(i))
                uxnew=ux*rvo%fcos-uy*rvo%fsin
                uynew=ux*rvo%fsin+uy*rvo%fcos
                tx(count) = (uxnew+rvo%xshift)/100.
                ty(count) = (uynew+rvo%yshift)/100.
                If(i.le.rvo%ns2) THEN
                    ttopo(count) = rvo%eta(i,j)/100.
                ELSE
                    ttopo(count) = rvo%eta(rvo%ns2,j)/100. - ((i-rvo%ns2)*rvo%ds*cco%nsextslope)
                ENDIF
            enddo
        enddo

        if(rvto%varDischType == 1) then !Discharge Time Series
            CALL getInterpTimeSeriesValue(rvto, 1, rvto%VarDischStartTime, tmpq)
        else
            tmpq = cco%q/1e6
        endif

        if(rvto%varStageType == 1) then !Stage Time Series
            CALL getInterpRatingCurveValue(rvto, 1, rvto%VarDischStartTime, tmpwselev)
            tmpwselev = tmpwselev/100.
            !        newStage = (newStage-elevoffset) * 100.
        else if (rvto%varStageType == 2) then !Stage Rating Curve
            Call getInterpRatingCurveValue(rvto, 1, tmpq, tmpwselev)
            tmpwselev = tmpwselev/100.
            !        newStage = (newStage-elevoffset) * 100.
        else
            tmpwselev = cco%wselev/100
        endif
        !    CALL NETPREIS2(tmpns, tmpnn, ttopo, tx, ty, ONEDCD, tmpq/1e6, tmpwselev/100., rvo%hav )
        CALL NETPREIS2(tmpns, tmpnn, ttopo, tx, ty, cco%ONEDCD, tmpq, tmpwselev, rvo%hav )
        rvo%hav = rvo%hav*100
        DO I = 1,tmpns
            DO J = 1,tmpnn
                rvo%e(i,j) = rvo%hav(i)
            ENDDO
        ENDDO

    ELSE IF(cco%wstype.eq.3) THEN
        CAll alloc_init2D(rvo, tmpns,tmpnn)
        open(501,file=cco%tmp_file_i,status='old',iostat = ier,form='unformatted')
        if(ier /= 0) then
            write(6,*) 'Input file error!'
            write(6,*) 'Temporary file for Hot Start does not exist'
            write(6,*) 'Press Enter to continue'
            read(6,*)
            stop
        end if
        !
        read(501) ttime,icount,dt2
        !
        !       time=(icount-1)*dt2
        !       icount=time/fmdt
        !
        read(501) nx2, ny2
        if(tmpns /= nx2.or.tmpnn /= ny2) then
            write(6,*) 'Number of grid is different between grid file and temporary file!'
            write(6,*) 'Press Enter to continue'
            read(6,*)
            stop
        end if
        !
        !       read(501) ((x(i,j),i=0,nx2),j=0,ny2)
        read(501) ((rvo%ie(i,j),i=1,nx2),j=1,ny2)
        read(501) ((rvo%iu(i,j),i=1,nx2),j=1,ny2)
        read(501) ((rvo%iv(i,j),i=1,nx2),j=1,ny2)
        read(501) ((rvo%iibc(i,j),i=1,nx2),j=1,ny2)
        read(501) ((rvo%iwse(i,j),i=1,nx2),j=1,ny2)
        read(501) ((rvo%ihl(i,j),i=1,nx2),j=1,ny2)
        read(501) (rvo%hav(i),i=1,nx2)

        Close(501)
        DO I = 1,tmpns
            DO J = 1,tmpnn
                rvo%u(i,j) = rvo%iu(i,j)
                rvo%v(i,j) = rvo%iv(i,j)
                rvo%e(i,j) = rvo%iwse(i,j)
                rvo%ibc(i,j) = rvo%iibc(i,j)
                rvo%eta(i,j) = rvo%ie(i,j)
                rvo%hl(i,j) = rvo%ihl(i,j)
                IF(j == tmpnn/2+1) THEN
                    rvo%hav(i) = rvo%e(i,j)
                ENDIF
            ENDDO
        ENDDO

        CALL dealloc_init2D(rvo)
    ENDIF
    DO i = 1,tmpns
        DO j = 1,tmpnn
            IF(j == 1) THEN
                rvo%ibc(i,j) = 0
            ELSE IF(j == tmpnn) THEN
                rvo%ibc(i,j) = 0

            ELSE IF(rvo%e(i,j).gt.rvo%eta(i,j)) THEN
                rvo%ibc(i,j) = -1
            ELSE
                rvo%ibc(i,j) = 0
            ENDIF
        ENDDO
    ENDDO
    END SUBROUTINE calcWSInitCond


    SUBROUTINE INITARRAYS(rvo, rwvo, cco)
    IMPLICIT NONE
    type(rivvar), intent(inout) :: rvo
    type(riv_w_var), intent(inout) :: rwvo
    type(calccond), intent(inout) :: cco
    INTEGER :: i,j
    REAL(KIND=mp) :: rcos, rsin, uu, vv
    integer :: tmpnn, tmpns
    tmpnn = rvo%nn
    tmpns = rvo%ns
    !    open(4,FILE='RadiusOfCurvature.txt')
    IF (cco%cdtype == 0) THEN !Constant CD
        DO i = 1, rvo%ns2
            DO j = 1, tmpnn
                rvo%znaught2(i,j) = cco%constcd*100.
                rvo%cd2(i,j) = cco%constcd
                rvo%cdv2(i,j) = 0
            END DO
        END DO
    ENDIF

    do i=1,tmpns
        rcos = cos(rvo%phirotation(i))
        rsin = sin(rvo%phirotation(i))

        !    lengthen grid to eliminate eddies at downstream end
        if(i.le.rvo%ns2) then
            rvo%r(i)=rvo%r2(i)
            rvo%w(i)=rvo%w2(i)
            !	        rvo%hav(i)=hav2(i)
            !	        xo(i)=xo2(i)
            !	        yo(i)=yo2(i)
        else
            rvo%r(i) = rvo%r2(i)
            !		   		r(i) = r2(rvo%ns2)+ (i-rvo%ns2)*dr

            rvo%w(i)=rvo%w2(rvo%ns2)
            !	        rvo%hav(i)=hav2(rvo%ns2)
        endif

        !
        do  j=1,tmpnn
            !   lengthen grid to elimate eddies at downstream end
            if(i.le.rvo%ns2) then
                rvo%eta(i,j)=rvo%eta2(i,j)
                !mineta(i,j) = mineta2(i,j)
                !	            ibc(i,j)=ibc2(i,j)
                !	            ribc(i,j) = ribc2(i,j)
                rvo%cd(i,j)=rvo%cd2(i,j)
                rvo%cdv(i,j)=rvo%cdv2(i,j)
                rvo%totcd(i,j) = rvo%cd(i,j)+rvo%cdv(i,j)
                rvo%znaught(i,j)=rvo%znaught2(i,j)
            else
                rvo%eta(i,j)=rvo%eta2(rvo%ns2,j) - ((i-rvo%ns2)*rvo%ds*cco%nsextslope)
                !mineta(i,j)=mineta2(rvo%ns2,j) - ((i-rvo%ns2)*ds*nsextslope)
                !	            ibc(i,j)=ibc2(rvo%ns2,j)
                !	            ribc(i,j) = ribc2(rvo%ns2, j)
                rvo%cd(i,j)=rvo%cd2(rvo%ns2,j)
                rvo%cdv(i,j)=rvo%cdv2(rvo%ns2,j)
                rvo%totcd(i,j) = rvo%cd(i,j)+rvo%cdv(i,j)
                rvo%znaught(i,j)=rvo%znaught2(rvo%ns2,j)
                rvo%fracs(i,j) = rvo%fracs(rvo%ns2,j)
                rvo%hfine(i,j) = rvo%hfine(rvo%ns2,j)
            endif
            rvo%rn(i,j)=1.-(j-rvo%nm)*rvo%dn/rvo%r(i)
            !		    if(j.eq.1)then
            !    	        write(4,*) i,r(i)/100.
            !    	    endif
            rwvo%dude(i,j)=0.
            rwvo%dvde(i,j)=0.
            rwvo%dsu(i,j)=rvo%ds*rvo%rn(i,j)
            rwvo%dsv(i,j)=rvo%ds*rvo%rn(i,j)
            rwvo%dse(i,j)=rvo%ds*rvo%rn(i,j)
            rwvo%dsq(i,j)=rvo%ds*rvo%rn(i,j)
            rwvo%dnu(i,j)=rvo%dn
            rwvo%dnv(i,j)=rvo%dn
            rwvo%dne(i,j)=rvo%dn
            rwvo%dnq(i,j)=rvo%dn
        ENDDO
        ENDDO
        !    close(4)
    END SUBROUTINE INITARRAYS

    END MODULE
