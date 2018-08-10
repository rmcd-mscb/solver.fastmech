    MODULE RivVarMod
    IMPLICIT NONE
    !	SAVE
    !	INTEGER, DIMENSION(ns, nn) :: ibc
    !Allocate REAL common types
    ! Machine precision:
    INTEGER, PARAMETER :: mp = KIND(1.0D0) ! = KIND(1.0D0) for double precision.
    integer,parameter :: strMax = 250
    INTEGER :: nclusters, nwetnodes, noldwetnodes

    REAL(kind=mp) :: ElevOffset=1e12
    REAL(kind=mp) :: ds,dn,mo,hwt
    REAL(kind=mp) :: fcos, fsin
    REAL(kind=mp) :: lbwse, rbwse, vdcds, vacds
    !			REAL :: hrwdepth
    !			REAL ::	startLEV, endLEV
    REAL(kind=mp) :: totTime, vardt, ptime
    REAL(kind=mp) :: wmax
    !Variables for CSed
    !			    REAL :: HD, WD, DIN
    !			    INTEGER :: SEDSMOO, SEDSMOOWGHT, TRANSEQTYPE, CALCCSED, SEDBCNODE
    !			    LOGICAL :: CALCSEDAUTO
    !			    INTEGER :: GRAVCORRTYPE, CALCGRAVCORR
    !			    REAL :: GRAVFLATBEDCORRCOEF, SUBANGLEOFREPOSE
    !			    REAL :: TSFracDepth, BCFRACTION
    !Variables for Vert
    !			    INTEGER ::  CALCQUASI3D
    !			    REAL :: MinRS
    !Allocate INTEGER common types
    INTEGER :: iter
    INTEGER :: nm
    INTEGER :: nsteps, nct
    INTEGER :: numVBCPts
    INTEGER :: errorcode
    LOGICAL :: FLUMEBNDRY
    REAL(KIND=mp) :: dDisch, dStage, newStage, oldStage, newDisch, oldDisch, dsstage


    !			INTEGER :: LEVEndIter, LEVBegIter, LEVChangeIter
    !!			INTEGER :: LEVType
    !!			INTEGER :: dryType
    !!			LOGICAL :: hcalcwetting



    DOUBLE PRECISION:: xshift, yshift

    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ibc, iibc, ribc, icon, tibc
    REAL(kind=mp), ALLOCATABLE, DIMENSION(:) :: r, w, xo, yo, hav
    REAL(kind=mp), ALLOCATABLE, DIMENSION(:,:) :: havn
    REAL(kind=mp), ALLOCATABLE, DIMENSION(:,:) :: taus, taun, hl, eta, rn, cd, cdv
    REAL(kind=mp), ALLOCATABLE, DIMENSION(:,:) :: totcd, mineta
    REAL(kind=mp), ALLOCATABLE, DIMENSION(:,:) :: znaught
    REAL(kind=mp), ALLOCATABLE, DIMENSION(:,:) :: u, v, e, iu, iv, iwse, ie, ihl
    REAL(kind=mp), ALLOCATABLE, DIMENSION(:,:) :: x, y
    REAL(kind=mp), ALLOCATABLE, DIMENSION(:,:) :: harea

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: uz, vz, zz
    REAL(kind=mp), ALLOCATABLE, DIMENSION(:,:) :: con, qs, qn

    REAL(kind=mp), ALLOCATABLE, DIMENSION(:,:) :: Fracs, hfine

    !Allocate Seperate arrays to deal with extend lower boundary
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ibc2
    REAL(kind=mp), ALLOCATABLE, DIMENSION(:) :: r2, w2, hav2
    !	REAL, ALLOCATABLE, DIMENSION(:,:) :: hav2n
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: xo2, yo2, phi2
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: phirotation
    REAL(kind=mp), ALLOCATABLE, DIMENSION(:,:) :: u2, v2, e2
    REAL(kind=mp), ALLOCATABLE, DIMENSION(:,:) :: eta2, cd2, cdv2, znaught2, mineta2

    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: iedger, iedgel
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: icountl, icountr
    INTEGER, ALLOCATABLE, DIMENSION(:) :: iswitch

    REAL(kind=mp) mo2



    INTEGER, PUBLIC :: ns2
    INTEGER, PUBLIC :: ns, nn, nz
    INTEGER, PUBLIC :: ns2a, nsa, nna !Added for CD Code (rmcd 12/25/02)

    INTEGER :: CGNSFILEID

    CONTAINS
    SUBROUTINE alloc_common3D(ns2, nn, nz)
    INTEGER, INTENT(IN) :: ns2, nn, nz
    INTEGER :: status, tmpns, tmp
    !			tmpns = ns2+nsext
    !			tmp = ns
    !Allocate REAL 3D common types
    ALLOCATE(uz(ns, nn, nz), STAT = status)
    ALLOCATE(vz(ns, nn, nz), STAT = status)
    ALLOCATE(zz(ns, nn, nz), STAT = status)
    uz = 0.
    vz = 0.
    zz = 0.
    END SUBROUTINE

    SUBROUTINE alloc_init2D(ns2, nn)
    INTEGER, INTENT(IN) :: ns2, nn
    INTEGER :: status
    ALLOCATE(iv(ns2, nn), STAT = status)
    ALLOCATE(iu(ns2, nn), STAT = status)
    ALLOCATE(ie(ns2, nn), STAT = status)
    ALLOCATE(iwse(ns2, nn), STAT = status)
    Allocate(ihl(ns2, nn), STAT = status)
    Allocate(iibc(ns2, nn), STAT = status)
    END SUBROUTINE

    SUBROUTINE dealloc_init2D()
    INTEGER :: status
    DEALLOCATE(iv, STAT = status)
    DEALLOCATE(iu, STAT = status)
    DEALLOCATE(iwse, STAT = status)
    DEALLOCATE(ie, STAT = status)
    DEALLOCATE(iibc, STAT = status)
    DEALLOCATE(ihl, STAT = status)
    END SUBROUTINE

    SUBROUTINE alloc_common2D(ns2, nsext, nn)
    INTEGER, INTENT(IN) :: ns2, nsext, nn
    INTEGER :: status
    ns = ns2 + nsext
    nsa = ns    !Added for CD Code (rmcd 12/25/02)
    nna = nn    !Added for CD Code (rmcd 12/25/02)
    ns2a = ns2  !Added for CD Code (rmcd 12/25/02)
    !Allocate INTEGER common types
    ALLOCATE(ibc(ns, nn), STAT = status)
    ALLOCATE(tibc(ns, nn), STAT = status)
    ALLOCATE(ribc(ns, nn), STAT = status)
    ALLOCATE(icon(ns, nn), STAT = status)
    !Allocate REAL 1D common types
    ALLOCATE(r(ns), STAT = status)
    ALLOCATE(w(ns), STAT = status)
    ALLOCATE(xo(ns), STAT = status)
    ALLOCATE(yo(ns), STAT = status)
    ALLOCATE(hav(ns), STAT = status)
    ALLOCATE(phirotation(ns), STAT = status)
    ALLOCATE(havn(ns, nn), STAT = status)
    !Allocate REAL 2D common types
    ALLOCATE(u(ns, nn), STAT = status)
    ALLOCATE(v(ns, nn), STAT = status)
    ALLOCATE(e(ns, nn), STAT = status)
    ALLOCATE(taus(ns, nn), STAT = status)
    ALLOCATE(taun(ns, nn), STAT = status)
    ALLOCATE(hl(ns, nn), STAT = status)
    ALLOCATE(eta(ns, nn), STAT = status)
    !ALLOCATE(mineta(ns, nn), STAT = status)
    ALLOCATE(rn(ns, nn), STAT = status)
    ALLOCATE(cd(ns, nn), STAT = status)
    ALLOCATE(cdv(ns, nn), STAT = status)
    ALLOCATE(totcd(ns, nn), STAT = status)
    ALLOCATE(znaught(ns, nn), STAT = status)
    ALLOCATE(con(ns, nn), STAT = status)
    ALLOCATE(qs(ns, nn), STAT = status)
    ALLOCATE(qn(ns, nn), STAT = status)
    ALLOCATE(x(ns, nn), STAT = status)
    ALLOCATE(y(ns, nn), STAT = status)
    ALLOCATE(harea(ns, nn), STAT = status)
    ALLOCATE(Fracs(ns, nn), STAT = status)
    ALLOCATE(hfine(ns, nn), STAT = status)

    !Allocate for CRAIG DIXON CHANGES
    ALLOCATE(iedger(ns, 15), STAT = status)
    ALLOCATE(iedgel(ns, 15), STAT = status)
    ALLOCATE(icountr(ns, nn), STAT = status)
    ALLOCATE(icountl(ns, nn), STAT = status)
    ALLOCATE(iswitch(0:nn), STAT = status)

    !Allocate variables for extended bounds
    ALLOCATE(r2(ns2+nsext), STAT = status)
    ALLOCATE(w2(ns2), STAT = status)
    ALLOCATE(xo2(ns2), STAT = status)
    ALLOCATE(yo2(ns2), STAT = status)
    ALLOCATE(hav2(ns2), STAT = status)
    !			ALLOCATE(hav2n(ns2, nn), STAT = status)
    ALLOCATE(phi2(ns2), STAT = status)
    ALLOCATE(cd2(ns2, nn), STAT = status)
    ALLOCATE(cdv2(ns2, nn), STAT = status)
    ALLOCATE(znaught2(ns2, nn), STAT = status)

    ALLOCATE(eta2(ns2, nn), STAT = status)
    !ALLOCATE(mineta2(ns2, nn), STAT = status)
    ALLOCATE(ibc2(ns2, nn), STAT = status)
    ALLOCATE(u2(ns2, nn), STAT = status)
    ALLOCATE(v2(ns2, nn), STAT = status)
    ALLOCATE(e2(ns2, nn), STAT = status)
    END SUBROUTINE

    SUBROUTINE dealloc_common3D()
    INTEGER :: status
    !Allocate REAL 3D common types
    DEALLOCATE(uz, STAT = status)
    DEALLOCATE(vz, STAT = status)
    DEALLOCATE(zz, STAT = status)
    END SUBROUTINE

    SUBROUTINE dealloc_common2D()
    INTEGER :: status
    !Allocate INTEGER common types
    DEALLOCATE(ibc, STAT = status)
    DEALLOCATE(tibc, STAT = status)
    DEALLOCATE(ribc, STAT = status)
    DEALLOCATE(icon, STAT = status)
    !Allocate REAL 1D common types
    DEALLOCATE(r, STAT = status)
    DEALLOCATE(w, STAT = status)
    DEALLOCATE(xo, STAT = status)
    DEALLOCATE(yo, STAT = status)
    DEALLOCATE(hav, STAT = status)
    DEALLOCATE(havn, STAT = status)
    !Allocate REAL 2D common types
    DEALLOCATE(u, STAT = status)
    DEALLOCATE(v, STAT = status)
    DEALLOCATE(e, STAT = status)
    DEALLOCATE(taus, STAT = status)
    DEALLOCATE(taun, STAT = status)
    DEALLOCATE(hl, STAT = status)
    DEALLOCATE(eta, STAT = status)
    !DEALLOCATE(mineta, STAT = status)
    DEALLOCATE(rn, STAT = status)
    DEALLOCATE(cd, STAT = status)
    DEALLOCATE(cdv, STAT = status)
    DEALLOCATE(totcd, STAT = status)
    DEALLOCATE(znaught, STAT = status)
    DEALLOCATE(con, STAT = status)
    DEALLOCATE(qs, STAT = status)
    DEALLOCATE(qn, STAT = status)
    DEALLOCATE(x, STAT = status)
    DEALLOCATE(y, STAT = status)
    DEALLOCATE(harea, STAT = status)
    DEALLOCATE(Fracs, stat = status)
    DEALLOCATE(hfine, stat = status)
    !Allocate for CRAIG DIXON CHANGES
    DEALLOCATE(iedger, STAT = status)
    DEALLOCATE(iedgel, STAT = status)
    DEALLOCATE(icountr, STAT = status)
    DEALLOCATE(icountl, STAT = status)
    DEALLOCATE(iswitch, STAT = status)

    !Allocate variables for extended bounds
    DEALLOCATE(r2, STAT = status)
    DEALLOCATE(w2, STAT = status)
    DEALLOCATE(xo2, STAT = status)
    DEALLOCATE(yo2, STAT = status)
    DEALLOCATE(hav2, STAT = status)
    !			DEALLOCATE(hav2n, STAT = status)
    DEALLOCATE(phi2, STAT = status)
    DEALLOCATE(phirotation, STAT = status)
    DEALLOCATE(cd2, STAT = status)
    DEALLOCATE(cdv2, STAT = status)
    DEALLOCATE(znaught2, STAT = status)

    DEALLOCATE(eta2, STAT = status)
    !DEALLOCATE(mineta2, STAT = status)
    DEALLOCATE(ibc2, STAT = status)
    DEALLOCATE(u2, STAT = status)
    DEALLOCATE(v2, STAT = status)
    DEALLOCATE(e2, STAT = status)
    !TO DO:  THESE are Temporarily Commented out
    !			IF(vbc == 1) THEN
    !				DEALLOCATE(vbcdist, STAT = status)
    !				DEALLOCATE(vbcvel, STAT = status)
    !				DEALLOCATE(vbcang, STAT = status)
    !			END IF
    END SUBROUTINE

    SUBROUTINE Calc_Area(x, y, xo, yo, nm, dn, area)
    IMPLICIT NONE
    REAL(KIND=mp), INTENT (OUT), DIMENSION (:,:) :: area
    REAL(KIND=mp), INTENT(INOUT), DIMENSION(:,:):: x,y
    REAL(KIND=mp), INTENT (IN), DIMENSION (:) :: xo, yo
    INTEGER, INTENT(IN) :: nm
    REAL(KIND=mp), INTENT(IN) :: dn

    INTEGER :: I, J
    REAL(KIND=mp) :: x1, x2, x3, x4
    REAL(KIND=mp) :: y1, y2, y3, y4
    REAL(KIND=mp) :: area1, area2

    REAL(KIND=mp) :: rsin, rcos, ux, uy, uxnew, uynew

    INTEGER :: icount
    do i=1,ns2
        rcos=cos(phirotation(i))
        rsin=sin(phirotation(i))

        do j=1,nn
            x(i,j)=xo(i)+(nm-j)*dn*sin(phirotation(i))
            y(i,j)=yo(i)+(j-nm)*dn*cos(phirotation(i))
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

    !SUBROUTINE RESETGRIDEXTENSION()
    !IMPLICIT NONE
    !
    !INTEGER :: i,j
    !REAL(KIND=mp) :: mag
    !do i=ns2+1,ns2+nsext
    !    do j = 1,nn
    !        eta(i,j)=eta2(ns2,j) - ((i-ns2)*ds*nsextslope)
    !        !mineta(i,j)=mineta2(ns2,j) - ((i-ns2)*ds*nsextslope)
    !        if(ibc(i,j).ne.0) then
    !            mag = sqrt(u(ns2,j)**2.+v(ns2,j)**2.)
    !            v(i,j) = v(i,j) - (i*(-v(ns2,j)/(ns2+nsext)))
    !            u(i,j) = sqrt(mag**2. - v(i,j)**2.)
    !        endif
    !    enddo
    !enddo
    !END SUBROUTINE RESETGRIDEXTENSION
    END MODULE RivVarMod
