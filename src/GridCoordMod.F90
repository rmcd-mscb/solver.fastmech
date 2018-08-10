    Module GridCoord
    USE RivVarMod
    USE RivVarWMod
    USE CSedMod_DT_SUSP
    USE CalcCond
    IMPLICIT NONE

    REAL(KIND=mp) :: stot, scals

    CONTAINS
    SUBROUTINE CGNS2_READ_GRIDCOORD(IER)
    IMPLICIT NONE
    !	INCLUDE "cgnslib_f.h"
    !    INCLUDE "cgnswin_f.h"

    INTEGER, INTENT(OUT) :: IER
    INTEGER :: status, i, j, count, countji, ierror
    INTEGER :: tmpnn, tmpns
    REAL(kind = mp), ALLOCATABLE, DIMENSION(:,:) :: tmpx1, tmpy1
    !	REAL*8, ALLOCATABLE, DIMENSION(:) :: tmpx, tmpy

    CALL CG_IRIC_GOTOGRIDCOORD2D_F(tmpns, tmpnn, IER)
    ns2 = tmpns
    ns = tmpns+nsext
    nn = tmpnn

    write(*,*) tmpns, nsext, tmpnn
    CALL alloc_common2D(tmpns, nsext, tmpnn)
    ! CALL ALLOC_INIT2D(ns2, nn)
    !	call alloc_csed_DT()
    CALL alloc_working()

    !	IF(CALCQUASI3D) THEN
    !	    CALL ALLOC_COMMON3D(ns2, nn, nz)
    !	    IF(TRANSEQTYPE == 2) THEN
    !		    CALL alloc_csed3d_dt() !array for Daniele Tonina Wilcock-Kenworthy
    !		ENDIF
    !    ENDIF

    !    write(*,*) 'before allocation \n'
    ALLOCATE(tmpx1(tmpns,tmpnn), STAT=IER)
    !    write(*,*) 'after allocation1 \n %d', ier
    ALLOCATE(tmpy1(tmpns,tmpnn), STAT=IER)
    !	write(*,*) 'after allocation2\n %d', ier
    CALL CG_IRIC_GETGRIDCOORD2D_F(tmpx1,tmpy1,IER)
    !	write(*,*) 'after iRIC CALL \n %d', ier
    !	CALL CG_IRIC_GETGRIDCOORD2D_F(x,y,IER)
    ns2 = tmpns
    nn = tmpnn
    DO I=1,ns2
        DO J=1,nn
            !	        COUNT = ((I-1)*nn)+J
            !	        countji = ((j-1)*ns2)+i
            !	        X(I,J) = tmpx(COUNTji)*100.
            !	        Y(I,J) = tmpy(COUNTji)*100.
            X(I,J) = tmpx1(i,j)*100.
            Y(I,J) = tmpy1(i,j)*100.
        ENDDO
    ENDDO
    !    x = x*100.
    !    y = y*100.

    CALL CalcGridMetrics()

    DEALLOCATE(tmpx1, STAT = IER)
    DEALLOCATE(tmpy1, STAT = IER)
    END SUBROUTINE CGNS2_READ_GRIDCOORD




    SUBROUTINE CalcGridMetrics()
    IMPLICIT NONE

    REAL(kind=mp) :: tTheta, tdy, tdx, xext, yext, tds, tdr, dtphi
    REAL(kind=mp) :: ds,slin, dx, dy, dphi, tdphi, tmpwidth
    DOUBLE PRECISION :: xc1, yc1

    INTEGER :: i

    INTEGER :: jmid

    jmid = (nn+1)/2
    r2 = 0
    ! Find length along centerline
    stot = 0.
    DO i = 2,ns2
        ds = sqrt((x(i,jmid)-x(i-1,jmid))**2 + (y(i,jmid)-y(i-1,jmid))**2)
        stot = stot+ds
    ENDDO


    scals = stot/(ns2-1)
    xshift = x(1,jmid)
    yshift = y(1,jmid)

    tmpwidth = sqrt((x(1,1)-x(1,nn))**2 + (y(1,1)-y(1,nn))**2)
    DO i = 1,ns2
        xo(i) = x(i,jmid)-xshift
        yo(i) = y(i,jmid)-yshift
        w2(i) = tmpwidth
    ENDDO

    slin = sqrt(xo(ns2)**2+yo(ns2)**2)
    fcos = xo(ns2)/slin
    fsin = yo(ns2)/slin

    DO i = 1,ns2
        xc1 = xo(i)*fcos+yo(i)*fsin
        yc1 = yo(i)*fcos-xo(i)*fsin
        xo(i) = xc1
        yo(i) = yc1
    ENDDO

    DO i=2,ns2
        dx = xo(i)-xo(i-1)
        dy = yo(i)-yo(i-1)
        IF(dx.eq.0) THEN
            IF(DY.gt.0) THEN
                phirotation(i) = acos(-1.)/2.
            ELSE
                phirotation(i) = -1. * acos(-1.)/2.
            ENDIF
        ELSE
            phirotation(i) = atan2(dy,dx)
        ENDIF
        phi2(i) = phirotation(i)
    ENDDO

    phirotation(1) = (2.*phirotation(2))-phirotation(3)
    phi2(1) = phirotation(1)

    DO i = 2,ns2
        dphi = phirotation(i)-phirotation(i-1)
        IF(dphi == 0.) THEN
            if(r2(i-1) < 0) THEN
                r2(i) = -100000000.
            ELSE
                r2(i) = 100000000.
            ENDIF
        ELSE
            r2(i) = scals/dphi
        ENDIF
    ENDDO

    r2(1) = (2.*r2(2))-r2(3)

    tdphi = phirotation(ns2)-phirotation(ns2-1)
    tdr = (tdphi)/nsext
    DO i=ns2+1,ns2+nsext
        dphi = ((ns2+nsext+1)-i)*tdr
        phirotation(i) = phirotation(i-1)+dphi

        IF(dphi == 0) THEN
            IF(r2(i-1) < 0) THEN
                r2(i) = -100000000.
            ELSE
                r2(i) = 100000000.
            ENDIF
        ELSE
            r2(i) = scals/dphi
        ENDIF
        dx = scals*cos(phirotation(i))
        dy = scals*sin(phirotation(i))
        xo(i) = xo(i-1) + dx
        yo(i) = yo(i-1) + dy
        dx = xo(i)-xo(i-1)
        dy = yo(i)-yo(i-1)
        IF(dx.eq.0) THEN
            IF(dy.gt.0) THEN
                phirotation(i)=acos(-1.)/2.
            ELSEIF(dy.le.0) THEN
                phirotation(i)=-1.*acos(-1.)/2.
            ENDIF
        ELSE
            phirotation(i)=atan2(dy,dx)
        ENDIF

    ENDDO

    END SUBROUTINE

    END MODULE GridCoord