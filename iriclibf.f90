    !  iricdllf.f90
    !
    !  FUNCTIONS/SUBROUTINES exported from iricdllf.dll:
    !  iricdllf      - subroutine
    !

    MODULE iriclibf
    IMPLICIT NONE
    INTEGER, PARAMETER :: mp = KIND(1.0D0)
    CHARACTER(LEN=21), PARAMETER :: CCNODE = 'CalculationConditions'
    CHARACTER(LEN=14), PARAMETER :: GCNODE = 'GridConditions'
    CHARACTER(LEN=12), PARAMETER :: SOLNODE = 'FlowSolution'
    INTEGER, PARAMETER :: NAME_MAXLENGTH = 250
    INTEGER :: FILEID, CELLDIM, PHYSDIM
    INTEGER :: BASEID, BASE3DID, ZONEID, ZONE3DID
    INTEGER :: FLOWSOL2DID, FLOWSOL3DID
    INTEGER :: ZONESIZE(9)
    INTEGER, DIMENSION(3,3) :: isize
    INTEGER, DIMENSION(3) :: irmin, irmax

    CONTAINS

    SUBROUTINE CG_IRIC_CREATE_3DBASEZONE_F(FID, NS, NN, NZ, X, Y, Z, IER)
    IMPLICIT NONE
    !! INCLUDE "cgnswin_f.h"
    INCLUDE "cgnslib_f.h"
    INTEGER, INTENT(IN) :: FID, NS, NN, NZ
    REAL(kind = mp), DIMENSION(:,:), INTENT(IN) :: X,Y
    REAL(kind = mp), DIMENSION(:,:,:), INTENT(IN) :: Z
    INTEGER, INTENT(OUT) :: IER
    INTEGER :: I,J,K, countkji, CID
    REAL(kind = mp), ALLOCATABLE, DIMENSION(:) :: tmp1, tmp2, tmp3

    celldim = 3
    physdim = 3

    isize(1,1) = ns
    isize(2,1) = nn
    isize(3,1) = nz

    isize(1,2) =  isize(1,1)-1
    isize(2,2) = isize(2,1)-1
    isize(3,2) = isize(3,1)-1

    isize(1,3) = 0
    isize(2,3) = 0
    isize(3,3) = 0

    CALL cg_base_write_f(FID, "iRIC3D", celldim, physdim, BASE3DID, ier )
    CALL cg_zone_write_f(FID, BASE3DID, "iRICZone", isize, Structured, ZONE3DID , ier )

    ALLOCATE(tmp1(ns*nn*nz), STAT = ier)
    ALLOCATE(tmp2(ns*nn*nz), STAT = ier)
    ALLOCATE(tmp3(ns*nn*nz), STAT = ier)


    DO K=1,NZ
        DO J=1,NN
            DO I=1,NS
                countkji = ((k-1)*ns*nn)+((j-1)*ns)+i
                tmp1(countkji) = x(i,j)/100.
                tmp2(countkji) = y(i,j)/100.
                tmp3(countkji) = z(i,j,k)/100.
            ENDDO
        ENDDO
    ENDDO

    call cg_coord_write_f(FID,BASE3DID,ZONE3DID,RealSingle,'CoordinateX',tmp1, CID, IER)
    call cg_coord_write_f(FID,BASE3DID,ZONE3DID,RealSingle,'CoordinateY',tmp2, CID, IER)
    call cg_coord_write_f(FID,BASE3DID,ZONE3DID,RealSingle,'CoordinateZ',tmp3, CID, IER)


    END SUBROUTINE CG_IRIC_CREATE_3DBASEZONE_F

    SUBROUTINE CG_IRIC_CREATE_3DFLOWSOL_F(FID, INDEX, NX, NY, NZ, IER)
    IMPLICIT NONE
    ! INCLUDE "cgnswin_f.h"
    INCLUDE "cgnslib_f.h"

    INTEGER, INTENT(IN) :: FID
    INTEGER, INTENT(IN) :: INDEX
    INTEGER, INTENT(OUT) :: NX, NY, NZ
    INTEGER, INTENT(OUT) :: IER
    ! Variables
    !INTEGER, PARAMETER :: NAME_MAXLENGTH = 250
    CHARACTER (LEN = 4), PARAMETER :: BASENAME = 'iRIC3D'
    CHARACTER (LEN = 8), PARAMETER :: ZONENAME = 'iRICZone'
    CHARACTER (LEN = NAME_MAXLENGTH) :: TMPBASENAME, TMPZONENAME
    INTEGER :: I, J, K, IERROR, NBASES, NZONES
    INTEGER :: BASES_ITER, ZONES_ITER
    CHARACTER(LEN = 14) :: flowsol3D
    CHARACTER(LEN = 5) :: tsindex
    CHARACTER(LEN = 32) :: tmp_flowsol3D
    flowsol3D = '3DFlowSolution'

    WRITE(tsindex, '(I5)') INDEX
    tsindex = ADJUSTL(tsindex)

    tmp_flowsol3D = flowsol3D//TRIM(tsindex)
    tmp_flowsol3D = TRIM(tmp_flowsol3D)

    ! Body of iricdllf
    FILEID = FID

    CALL cg_nbases_f(FID,NBASES, IER)
    IF(IER.ne.0) RETURN

    DO BASES_ITER = 1, NBASES
        CALL cg_base_read_f(FID, BASES_ITER, TMPBASENAME, CELLDIM, PHYSDIM, IER)
        IF(IER.ne.0) RETURN
        SELECT CASE(TRIM(TMPBASENAME))
        CASE('iRIC3D')
            BASEID = BASES_ITER
            CALL cg_nzones_f(FID, BASES_ITER, NZONES, IER)
            IF(IER.ne.0) RETURN
            DO ZONES_ITER = 1, NZONES
                CALL cg_zone_read_f(FID, BASES_ITER, ZONES_ITER, TMPZONENAME, isize, IER)
                irmin(1)=1
                irmin(2)=1
                irmin(3)=1
                irmax(1)=isize(1,1)
                irmax(2)=isize(2,1)
                irmax(3)=isize(3,1)

                IF(IER.ne.0) RETURN
                SELECT CASE(TRIM(TMPZONENAME))
                CASE('iRICZone')
                    ZONEID = ZONES_ITER
                    NX = isize(1,1)
                    NY = isize(2,1)
                    NZ = isize(3,1)
                    CALL cg_sol_write_f(FID, BASES_ITER, ZONES_ITER, tmp_flowsol3D, Vertex, FLOWSOL3DID, IER)

                END SELECT
            ENDDO
        END SELECT
    ENDDO
    END SUBROUTINE CG_IRIC_CREATE_3DFLOWSOL_F

    SUBROUTINE CG_IRIC_CREATEFLOWSOL_F(FID, INDEX, NX, NY, IER)
    IMPLICIT NONE
    ! INCLUDE "cgnswin_f.h"
    INCLUDE "cgnslib_f.h"

    INTEGER, INTENT(IN) :: FID
    INTEGER, INTENT(IN) :: INDEX
    INTEGER, INTENT(OUT) :: NX, NY
    INTEGER, INTENT(OUT) :: IER
    ! Variables
    !INTEGER, PARAMETER :: NAME_MAXLENGTH = 250
    CHARACTER (LEN = 4), PARAMETER :: BASENAME = 'iRIC'
    CHARACTER (LEN = 8), PARAMETER :: ZONENAME = 'iRICZone'
    CHARACTER (LEN = NAME_MAXLENGTH) :: TMPBASENAME, TMPZONENAME
    INTEGER :: I, J, IERROR, NBASES, NZONES
    INTEGER :: BASES_ITER, ZONES_ITER
    CHARACTER(LEN = 14) :: flowsol2D
    CHARACTER(LEN = 5) :: tsindex
    CHARACTER(LEN = 32) :: tmp_flowsol2D
    flowsol2D = '2DFlowSolution'

    WRITE(tsindex, '(I5)') INDEX
    tsindex = ADJUSTL(tsindex)

    tmp_flowsol2D = flowsol2D//TRIM(tsindex)
    tmp_flowsol2D = TRIM(tmp_flowsol2D)

    ! Body of iricdllf
    FILEID = FID

    CALL cg_nbases_f(FID,NBASES, IER)
    IF(IER.ne.0) RETURN

    DO BASES_ITER = 1, NBASES
        CALL cg_base_read_f(FID, BASES_ITER, TMPBASENAME, CELLDIM, PHYSDIM, IER)
        IF(IER.ne.0) RETURN
        SELECT CASE(TRIM(TMPBASENAME))
        CASE('iRIC')
            BASEID = BASES_ITER
            CALL cg_nzones_f(FID, BASES_ITER, NZONES, IER)
            IF(IER.ne.0) RETURN
            DO ZONES_ITER = 1, NZONES
                CALL cg_zone_read_f(FID, BASES_ITER, ZONES_ITER, TMPZONENAME, isize, IER)
                irmin(1)=1
                irmin(2)=1
                irmin(3)=1
                irmax(1)=isize(1,1)
                irmax(2)=isize(2,1)
                irmax(3)=isize(3,1)

                IF(IER.ne.0) RETURN
                SELECT CASE(TRIM(TMPZONENAME))
                CASE('iRICZone')
                    ZONEID = ZONES_ITER
                    NX = isize(1,1)
                    NY = isize(2,1)
                    CALL cg_sol_write_f(FID, BASES_ITER, ZONES_ITER, tmp_flowsol2D, Vertex, FLOWSOL2DID, IER)

                END SELECT
            ENDDO
        END SELECT
    ENDDO
    END SUBROUTINE CG_IRIC_CREATEFLOWSOL_F

    SUBROUTINE CG_IRIC_GOTOSOL_F(FID, INDEX, NX, NY, IER)
    IMPLICIT NONE
    ! INCLUDE "cgnswin_f.h"
    INCLUDE "cgnslib_f.h"

    INTEGER, INTENT(IN) :: FID
    INTEGER, INTENT(IN) :: INDEX
    INTEGER, INTENT(OUT) :: NX, NY
    INTEGER, INTENT(OUT) :: IER
    ! Variables
    !INTEGER, PARAMETER :: NAME_MAXLENGTH = 250
    CHARACTER (LEN = 4), PARAMETER :: BASENAME = 'iRIC'
    CHARACTER (LEN = 8), PARAMETER :: ZONENAME = 'iRICZone'
    CHARACTER (LEN = NAME_MAXLENGTH) :: TMPBASENAME, TMPZONENAME
    INTEGER :: I, J, IERROR, NBASES, NZONES
    INTEGER :: BASES_ITER, ZONES_ITER
    CHARACTER(LEN = 14) :: flowsol2D
    CHARACTER(LEN = 5) :: tsindex
    CHARACTER(LEN = 32) :: tmp_flowsol2D
    flowsol2D = 'FlowSolution'

    WRITE(tsindex, '(I5)') INDEX
    tsindex = ADJUSTL(tsindex)

    tmp_flowsol2D = flowsol2D//TRIM(tsindex)
    tmp_flowsol2D = TRIM(tmp_flowsol2D)

    ! Body of iricdllf
    FILEID = FID

    CALL cg_nbases_f(FID,NBASES, IER)
    IF(IER.ne.0) RETURN

    DO BASES_ITER = 1, NBASES
        CALL cg_base_read_f(FID, BASES_ITER, TMPBASENAME, CELLDIM, PHYSDIM, IER)
        IF(IER.ne.0) RETURN
        SELECT CASE(TRIM(TMPBASENAME))
        CASE('iRIC')
            BASEID = BASES_ITER
            CALL cg_nzones_f(FID, BASES_ITER, NZONES, IER)
            IF(IER.ne.0) RETURN
            DO ZONES_ITER = 1, NZONES
                CALL cg_zone_read_f(FID, BASES_ITER, ZONES_ITER, TMPZONENAME, isize, IER)
                irmin(1)=1
                irmin(2)=1
                irmin(3)=1
                irmax(1)=isize(1,1)
                irmax(2)=isize(2,1)
                irmax(3)=isize(3,1)

                IF(IER.ne.0) RETURN
                SELECT CASE(TRIM(TMPZONENAME))
                CASE('iRICZone')
                    ZONEID = ZONES_ITER
                    NX = isize(1,1)
                    NY = isize(2,1)

                    CALL LOCAL_GOTOFS(INDEX, IER)
                END SELECT
            ENDDO
        END SELECT
    ENDDO
    END SUBROUTINE CG_IRIC_GOTOSOL_F

    SUBROUTINE CG_IRIC_GOTOGC_F(FID, NX, NY, IER)
    IMPLICIT NONE
    ! INCLUDE "cgnswin_f.h"
    INCLUDE "cgnslib_f.h"

    INTEGER, INTENT(IN) :: FID
    INTEGER, INTENT(OUT) :: NX, NY
    INTEGER, INTENT(OUT) :: IER
    ! Variables
    !INTEGER, PARAMETER :: NAME_MAXLENGTH = 250
    CHARACTER (LEN = 4), PARAMETER :: BASENAME = 'iRIC'
    CHARACTER (LEN = 8), PARAMETER :: ZONENAME = 'iRICZone'
    CHARACTER (LEN = NAME_MAXLENGTH) :: TMPBASENAME, TMPZONENAME
    INTEGER :: I, J, IERROR, NBASES, NZONES
    INTEGER :: BASES_ITER, ZONES_ITER

    ! Body of iricdllf
    FILEID = FID

    CALL cg_nbases_f(FID,NBASES, IER)
    IF(IER.ne.0) RETURN

    DO BASES_ITER = 1, NBASES
        CALL cg_base_read_f(FID, BASES_ITER, TMPBASENAME, CELLDIM, PHYSDIM, IER)
        IF(IER.ne.0) RETURN
        SELECT CASE(TRIM(TMPBASENAME))
        CASE('iRIC')
            BASEID = BASES_ITER
            CALL cg_nzones_f(FID, BASES_ITER, NZONES, IER)
            IF(IER.ne.0) RETURN
            DO ZONES_ITER = 1, NZONES
                CALL cg_zone_read_f(FID, BASES_ITER, ZONES_ITER, TMPZONENAME, isize, IER)
                irmin(1)=1
                irmin(2)=1
                irmin(3)=1
                irmax(1)=isize(1,1)
                irmax(2)=isize(2,1)
                irmax(3)=isize(3,1)

                IF(IER.ne.0) RETURN
                SELECT CASE(TRIM(TMPZONENAME))
                CASE('iRICZone')
                    ZONEID = ZONES_ITER
                    NX = isize(1,1)
                    NY = isize(2,1)

                    CALL LOCAL_GOTOGC(IER)
                END SELECT
            ENDDO
        END SELECT
    ENDDO
    END SUBROUTINE CG_IRIC_GOTOGC_F

    SUBROUTINE CG_IRIC_GOTOCC_F(FID, IER)
    IMPLICIT NONE

    ! INCLUDE "cgnswin_f.h"
    INCLUDE "cgnslib_f.h"
    ! Expose subroutine iricdllf to users of this DLL
    !
    INTEGER, INTENT(IN) :: FID
    INTEGER, INTENT(OUT) :: IER

    ! Variables
    !INTEGER, PARAMETER :: NAME_MAXLENGTH = 250
    CHARACTER (LEN = 4), PARAMETER :: BASENAME = 'iRIC'
    CHARACTER (LEN = 8), PARAMETER :: ZONENAME = 'iRICZone'
    CHARACTER (LEN = NAME_MAXLENGTH) :: TMPBASENAME, TMPZONENAME
    INTEGER :: I, J, IERROR, NBASES, NZONES
    INTEGER :: BASES_ITER, ZONES_ITER

    ! Body of iricdllf
    FILEID = FID

    CALL cg_nbases_f(FID,NBASES, IER)
    IF(IER.ne.0) RETURN

    DO BASES_ITER = 1, NBASES
        CALL cg_base_read_f(FID, BASES_ITER, TMPBASENAME, CELLDIM, PHYSDIM, IER)
        IF(IER.ne.0) RETURN
        SELECT CASE(TRIM(TMPBASENAME))
        CASE('iRIC')
            BASEID = BASES_ITER
            CALL LOCAL_GOTOCC(IER)
            !			    CALL cg_nzones_f(FID, BASES_ITER, NZONES, IER)
            !                IF(IER.ne.0) RETURN
            !			    DO ZONES_ITER = 1, NZONES
            !				    CALL cg_zone_read_f(FID, BASES_ITER, ZONES_ITER, TMPZONENAME, isize, IER)
            !				    irmin(1)=1
            !				    irmin(2)=1
            !				    irmin(3)=1
            !				    irmax(1)=isize(1,1)
            !				    irmax(2)=isize(2,1)
            !				    irmax(3)=isize(3,1)
            !                    IF(IER.ne.0) RETURN
            !				    SELECT CASE(TRIM(TMPZONENAME))
            !				    CASE('iRICZone')
            !                        ZONEID = ZONES_ITER
            !                        CALL LOCAL_GOTOCC(IER)
            !                    END SELECT
            !               ENDDO
        END SELECT
    ENDDO

    end subroutine CG_IRIC_GOTOCC_F

    SUBROUTINE CG_IRIC_GOTOGRIDCOORD_F(FID, NX, NY, IER)
    IMPLICIT NONE

    ! INCLUDE "cgnswin_f.h"
    INCLUDE "cgnslib_f.h"
    ! Expose subroutine iricdllf to users of this DLL
    !
    INTEGER, INTENT(IN) :: FID
    INTEGER, INTENT(OUT) :: NX, NY, IER

    ! Variables
    !INTEGER, PARAMETER :: NAME_MAXLENGTH = 250
    CHARACTER (LEN = 4), PARAMETER :: BASENAME = 'iRIC'
    CHARACTER (LEN = 8), PARAMETER :: ZONENAME = 'iRICZone'
    CHARACTER (LEN = NAME_MAXLENGTH) :: TMPBASENAME, TMPZONENAME
    INTEGER :: I, J, IERROR, NBASES, NZONES
    INTEGER :: BASES_ITER, ZONES_ITER

    ! Body of iricdllf
    FILEID = FID

    CALL cg_nbases_f(FID,NBASES, IER)
    IF(IER.ne.0) RETURN

    DO BASES_ITER = 1, NBASES
        CALL cg_base_read_f(FID, BASES_ITER, TMPBASENAME, CELLDIM, PHYSDIM, IER)
        IF(IER.ne.0) RETURN
        SELECT CASE(TRIM(TMPBASENAME))
        CASE('iRIC')
            BASEID = BASES_ITER
            CALL cg_nzones_f(FID, BASES_ITER, NZONES, IER)
            IF(IER.ne.0) RETURN
            DO ZONES_ITER = 1, NZONES
                CALL cg_zone_read_f(FID, BASES_ITER, ZONES_ITER, TMPZONENAME, isize, IER)
                irmin(1)=1
                irmin(2)=1
                irmin(3)=1
                irmax(1)=isize(1,1)
                irmax(2)=isize(2,1)
                irmax(3)=isize(3,1)

                IF(IER.ne.0) RETURN
                SELECT CASE(TRIM(TMPZONENAME))
                CASE('iRICZone')
                    ZONEID = ZONES_ITER
                    NX = isize(1,1)
                    NY = isize(2,1)

                    !                        CALL LOCAL_GOTOCC(IER)
                END SELECT
            ENDDO
        END SELECT
    ENDDO

    end subroutine CG_IRIC_GOTOGRIDCOORD_F

    SUBROUTINE CG_IRIC_GETGRIDCOORD_F(NX,NY,X,Y,IER)
    IMPLICIT NONE

    ! INCLUDE "cgnswin_f.h"
    INCLUDE "cgnslib_f.h"

    INTEGER, INTENT(IN) :: NX, NY
    REAL*4, DIMENSION(NX*NY), INTENT(OUT) :: X,Y
    INTEGER, INTENT(OUT) :: IER
    CHARACTER (LEN = NAME_MAXLENGTH) :: coordname
    INTEGER :: ncoord
    INTEGER :: datatype, n_pts
    INTEGER :: k

    n_pts = nx*ny
    ! Read nodal coordinates.
    CALL cg_ncoords_f(FILEID,BASEID,ZONEID,ncoord,IER)
    DO k = 1,ncoord
        CALL cg_coord_info_f(FILEID,BASEID,ZONEID,k,datatype,coordname,IER)

        SELECT CASE(TRIM(coordname))
        CASE('CoordinateX')
            CALL cg_coord_read_f(FILEID,BASEID,ZONEID,coordname,RealSingle,irmin, &
                irmax,X,IER)
        CASE('CoordinateY')
            CALL cg_coord_read_f(FILEID,BASEID,ZONEID,coordname,RealSingle,irmin, &
                irmax,Y,IER)
        END SELECT
    END DO

    END SUBROUTINE CG_IRIC_GETGRIDCOORD_F

    SUBROUTINE LOCAL_GOTOBASE(IER)
    IMPLICIT NONE
    ! INCLUDE "cgnswin_f.h"
    INCLUDE "cgnslib_f.h"
    INTEGER, INTENT(OUT) :: IER
    CALL CG_GOTO_F(FILEID, BASEID, IER, 'end')

    ENDSUBROUTINE LOCAL_GOTOBASE

    SUBROUTINE LOCAL_GOTOZONE(IER)
    IMPLICIT NONE
    ! INCLUDE "cgnswin_f.h"
    INCLUDE "cgnslib_f.h"
    INTEGER, INTENT(OUT) :: IER
    CALL CG_GOTO_F(FILEID, BASEID, IER, 'Zone_t', ZONEID, 'end')

    ENDSUBROUTINE

    SUBROUTINE LOCAL_GOTOS(PATH, IER)
    IMPLICIT NONE
    ! INCLUDE "cgnswin_f.h"
    INCLUDE "cgnslib_f.h"
    CHARACTER(*), INTENT(IN) :: PATH
    INTEGER, INTENT(OUT) :: IER
    CALL LOCAL_GOTOZONE(IER)
    IF(IER.NE.0) RETURN
    CALL cg_gopath_f(FILEID, PATH, IER)
    IF(IER.NE.0) RETURN
    END SUBROUTINE

    SUBROUTINE LOCAL_GOTOZ(PATH, IER)
    IMPLICIT NONE
    ! INCLUDE "cgnswin_f.h"
    INCLUDE "cgnslib_f.h"
    CHARACTER(*), INTENT(IN) :: PATH
    INTEGER, INTENT(OUT) :: IER
    CALL LOCAL_GOTOZONE(IER)
    IF(IER.NE.0) RETURN
    CALL cg_gopath_f(FILEID, PATH, IER)
    IF(IER.NE.0) RETURN
    END SUBROUTINE

    SUBROUTINE LOCAL_GOTOB(PATH, IER)
    IMPLICIT NONE
    ! INCLUDE "cgnswin_f.h"
    INCLUDE "cgnslib_f.h"
    CHARACTER(*), INTENT(IN) :: PATH
    INTEGER, INTENT(OUT) :: IER

    CHARACTER(LEN=250) :: TmpStr

    CALL LOCAL_GOTOBASE(IER)
    IF(IER.NE.0) RETURN
    Tmpstr = REPEAT(path, 1)
    CALL cg_gopath_f(FILEID,Path, ier)
    IF(IER.NE.0) RETURN
    END SUBROUTINE LOCAL_GOTOB

    SUBROUTINE LOCAL_GOTOCC(IER)
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: IER
    CALL LOCAL_GOTOB(CCNODE, IER)

    END SUBROUTINE LOCAL_GOTOCC

    SUBROUTINE LOCAL_GOTOGC(IER)
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: IER
    CALL LOCAL_GOTOZ(GCNODE, IER)

    END SUBROUTINE LOCAL_GOTOGC

    SUBROUTINE LOCAL_GOTOFS(INDEX, IER)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: INDEX
    INTEGER, INTENT(OUT) :: IER
    CHARACTER(LEN = 12) :: flowsol2D
    CHARACTER(LEN = 5) :: tsindex
    CHARACTER(LEN = 30) :: tmp_flowsol2D
    flowsol2D = 'FlowSolution'

    WRITE(tsindex, '(I5)') INDEX
    tsindex = ADJUSTL(tsindex)

    tmp_flowsol2D = flowsol2D//TRIM(tsindex)
    tmp_flowsol2D = TRIM(tmp_flowsol2D)

    CALL LOCAL_GOTOS(tmp_flowsol2D, IER)


    END SUBROUTINE LOCAL_GOTOFS

    !SUBROUTINE LOCAL_CREATEFLOWSOL(INDEX, IER)
    !IMPLICIT NONE
    !INTEGER, INTENT(IN) :: INDEX
    !INTEGER, INTENT(OUT) :: IER
    !CHARACTER(LEN = 14) :: flowsol2D
    !CHARACTER(LEN = 5) :: tsindex
    !CHARACTER(LEN = 32) :: tmp_flowsol2D,
    !flowsol2D = '2DFlowSolution'
    !
    !WRITE(tsindex, '(I5)') INDEX
    !tsindex = ADJUSTL(tsindex)
    !
    !tmp_flowsol2D = flowsol2D//TRIM(tsindex)
    !tmp_flowsol2D = TRIM(tmp_flowsol2D)
    !
    !
    !    CALL LOCAL_GOTO(tmp_flowsol2D, IER)
    !
    !END SUBROUTINE LOCAL_GOTOFLOWSOL

    SUBROUTINE LOCAL_GOTOCCCHILD(PATH, IER)
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: PATH
    INTEGER, INTENT(OUT) :: IER
    !	char tmpstr[NAME_MAXLENGTH];
    !	int slen;
    !	strcpy(tmpstr, CCNODE);
    !	slen = strlen(CCNODE);
    !	strcpy(tmpstr + slen, "/");
    !	strcpy(tmpstr + slen + 1, path);
    !	return local_goto(tmpstr);

    CHARACTER(LEN=NAME_MAXLENGTH) :: tmpstr
    INTEGER :: slen

    tmpstr = REPEAT(CCNODE,1)
    slen = LEN(TRIM(tmpstr))

    tmpstr = TRIM(tmpstr) // '/' // PATH

    CALL LOCAL_GOTOB(tmpstr, IER)


    ENDSUBROUTINE LOCAL_GOTOCCCHILD

    SUBROUTINE LOCAL_GOTOGCCHILD(PATH, IER)
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: PATH
    INTEGER, INTENT(OUT) :: IER
    !	char tmpstr[NAME_MAXLENGTH];
    !	int slen;
    !	strcpy(tmpstr, CCNODE);
    !	slen = strlen(CCNODE);
    !	strcpy(tmpstr + slen, "/");
    !	strcpy(tmpstr + slen + 1, path);
    !	return local_goto(tmpstr);

    CHARACTER(LEN=NAME_MAXLENGTH) :: tmpstr
    INTEGER :: slen

    tmpstr = REPEAT(GCNODE,1)
    slen = LEN(TRIM(tmpstr))

    tmpstr = TRIM(tmpstr) // '/' // PATH

    CALL LOCAL_GOTOZ(tmpstr, IER)


    ENDSUBROUTINE LOCAL_GOTOGCCHILD

    SUBROUTINE LOCAL_GOTOFSCHILD(INDEX, PATH, IER)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: INDEX
    CHARACTER(*), INTENT(IN) :: PATH
    INTEGER, INTENT(OUT) :: IER
    CHARACTER(LEN=NAME_MAXLENGTH) :: tmpstr
    INTEGER :: slen
    CHARACTER(LEN = 14) :: flowsol2D
    CHARACTER(LEN = 5) :: tsindex
    CHARACTER(LEN = 32) :: tmp_flowsol2D
    flowsol2D = '2DFlowSolution'

    WRITE(tsindex, '(I5)') INDEX
    tsindex = ADJUSTL(tsindex)

    tmp_flowsol2D = flowsol2D//TRIM(tsindex)
    tmpstr = TRIM(tmp_flowsol2D) // '/' // PATH

    CALL LOCAL_GOTOZ(tmpstr, IER)


    ENDSUBROUTINE LOCAL_GOTOFSCHILD


    !SUBROUTINE LOCAL__GETGRIDSIZE(GSIZE, IER)
    !IMPLICIT NONE
    !INTEGER, INTENT(OUT) :: GSIZE
    !INTEGER, INTENT(OUT) :: IER
    !
    !END SUBROUTINE LOCAL__GETGRIDSIZE
    SUBROUTINE cg_iRIC_Read_SolIntegerNode(index, name, intarraypointer, IER)
    IMPLICIT NONE
    ! INCLUDE "cgnswin_f.h"
    INCLUDE "cgnslib_f.h"
    INTEGER, INTENT(IN) :: index
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, DIMENSION(:), INTENT(OUT) :: intarraypointer
    INTEGER, INTENT(out) :: ier

    INTEGER :: narrays
    INTEGER :: datatype, dim, dimvec
    INTEGER :: I
    CHARACTER(LEN=NAME_MAXLENGTH) :: arrayname

    !    CALL local_gotofschild(index, name, ier)
    !    IF(ier.ne.0) RETURN
    !
    !	// under this UserDefinedData_t, array with name "Value" will exists.
    CALL cg_narrays_f(narrays, ier)
    IF(ier.ne.0) RETURN

    DO I=1,narrays
        CALL cg_array_info_f(i, arrayname, datatype, dim, dimvec, IER)
        IF(ier.ne.0) RETURN
        IF(TRIM(arrayname) == TRIM(name)) THEN
            CALL cg_array_read_f(i, intarraypointer, ier)
        ENDIF
    ENDDO


    END SUBROUTINE cg_iRIC_Read_SolIntegerNode

    SUBROUTINE cg_iRIC_Read_SolRealNode(index, name, realarraypointer, IER)
    IMPLICIT NONE
    ! INCLUDE "cgnswin_f.h"
    INCLUDE "cgnslib_f.h"
    INTEGER, INTENT(IN) :: index
    CHARACTER(LEN=*), INTENT(IN) :: name
    Real*8, DIMENSION(:), INTENT(OUT) :: realarraypointer
    INTEGER, INTENT(out) :: ier

    INTEGER :: narrays, found
    INTEGER :: datatype, dim, dimvec
    INTEGER :: I
    CHARACTER(LEN=NAME_MAXLENGTH) :: arrayname

    !CALL local_gotofschild(index, name, ier)
    !IF(ier.ne.0) RETURN

    !	// under this UserDefinedData_t, array with name "Value" will exists.
    CALL cg_narrays_f(narrays, ier)
    IF(ier.ne.0) RETURN

    DO I=1,narrays
        CALL cg_array_info_f(i, arrayname, datatype, dim, dimvec, IER)
        IF(ier.ne.0) RETURN
        ier = -1
        IF(TRIM(arrayname) == TRIM(name)) THEN
            CALL cg_array_read_f(i, realarraypointer, ier)
            EXIT
        ENDIF
    ENDDO


    END SUBROUTINE cg_iRIC_Read_SolRealNode

    SUBROUTINE CG_IRIC_READ_GRIDREALNODE(name, realarraypointer, IER)
    IMPLICIT NONE
    ! INCLUDE "cgnswin_f.h"
    INCLUDE "cgnslib_f.h"
    CHARACTER(LEN=*), INTENT(IN) :: name
    Real*4, DIMENSION(:), INTENT(OUT) :: realarraypointer
    INTEGER, INTENT(out) :: ier

    INTEGER :: narrays, found
    INTEGER :: datatype, dim, dimvec
    INTEGER :: I
    CHARACTER(LEN=NAME_MAXLENGTH) :: arrayname

    CALL local_gotogcchild(name, ier)
    IF(ier.ne.0) RETURN

    !	// under this UserDefinedData_t, array with name "Value" will exists.
    CALL cg_narrays_f(narrays, ier)
    IF(ier.ne.0) RETURN

    found = 0;
    DO I=1,narrays
        CALL cg_array_info_f(i, arrayname, datatype, dim, dimvec, IER)
        IF(ier.ne.0) RETURN
        SELECT CASE(TRIM(arrayname))
        CASE('Value')
            found = 1
            !			// the Value is found. dim must be 1, dimvec must be 1, and dataype must be real.
            !			// TODO: we have to check!
            CALL cg_array_read_f(i, realarraypointer, ier)
            IF(ier.ne.0) THEN
                found = 0;
                RETURN
            ENDIF
        END SELECT
    ENDDO


    END SUBROUTINE cg_iRIC_Read_GridRealNode

    SUBROUTINE CG_IRIC_READ_INTEGER(name, intvalue, ier)
    IMPLICIT NONE
    ! INCLUDE "cgnswin_f.h"
    INCLUDE "cgnslib_f.h"
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(OUT) :: intvalue
    INTEGER, INTENT(out) :: ier
    CHARACTER(LEN=NAME_MAXLENGTH) :: arrayname
    INTEGER :: datatype, arrays, dim, dimvec, ierror, i, found

    CALL local_gotoccchild(name, ier)
    IF(ier.ne.0) RETURN
    !	// under this UserDefinedData_t, array with name "Value" will exists.
    CALL cg_narrays_f(arrays, ier)
    IF(ier.ne.0) RETURN
    found = 0
    DO i = 1,arrays
        CALL cg_array_info_f(i, arrayname, datatype, dim, dimvec, ier)
        IF(ier.ne.0) RETURN
        SELECT CASE(TRIM(arrayname))
        CASE('Value')
            found = 1
            !			// the Value is found. dim must be 1, dimvec must be 1, and dataype must be real.
            !			// TODO: we have to check!
            CALL cg_array_read_f(i, intvalue, ier)
            IF(ier.ne.0) THEN
                found = 0;
                RETURN
            ENDIF
        END SELECT
    ENDDO

    END SUBROUTINE

    SUBROUTINE CG_IRIC_READ_SINGLE(name, realvalue, ier)
    IMPLICIT NONE
    ! INCLUDE "cgnswin_f.h"
    INCLUDE "cgnslib_f.h"
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL*4, INTENT(OUT) :: realvalue
    INTEGER, INTENT(OUT) :: ier
    CHARACTER(LEN=NAME_MAXLENGTH) :: arrayname
    INTEGER :: datatype, arrays, dim, dimvec, i, found
    CALL local_gotoccchild(name, ier)
    IF(ier.ne.0) RETURN
    !	// under this UserDefinedData_t, array with name "Value" will exists.
    CALL cg_narrays_f(arrays, ier)
    IF(ier.ne.0) RETURN
    found = 0
    DO i = 1,arrays
        CALL cg_array_info_f(i, arrayname, datatype, dim, dimvec, ier)
        IF(ier.ne.0) RETURN
        SELECT CASE(TRIM(arrayname))
        CASE('Value')
            found = 1
            !			// the Value is found. dim must be 1, dimvec must be 1, and dataype must be real.
            !			// TODO: we have to check!
            CALL cg_array_read_f(i, realvalue, ier)
            IF(ier.ne.0) THEN
                found = 0;
                RETURN
            ENDIF
        END SELECT
    ENDDO
    END SUBROUTINE

    SUBROUTINE CG_IRIC_READ_DOUBLE(name, dvalue, ier)
    IMPLICIT NONE
    ! INCLUDE "cgnswin_f.h"
    INCLUDE "cgnslib_f.h"
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL*8, INTENT(OUT) :: dvalue
    INTEGER, INTENT(OUT) :: ier
    CHARACTER(LEN=NAME_MAXLENGTH) :: arrayname
    INTEGER :: datatype, arrays, dim, dimvec, i, found
    CALL local_gotoccchild(name, ier)
    IF(ier.ne.0) RETURN
    !	// under this UserDefinedData_t, array with name "Value" will exists.
    CALL cg_narrays_f(arrays, ier)
    IF(ier.ne.0) RETURN
    found = 0
    DO i = 1,arrays
        CALL cg_array_info_f(i, arrayname, datatype, dim, dimvec, ier)
        IF(ier.ne.0) RETURN
        SELECT CASE(TRIM(arrayname))
        CASE('Value')
            found = 1
            !			// the Value is found. dim must be 1, dimvec must be 1, and dataype must be real.
            !			// TODO: we have to check!
            CALL cg_array_read_f(i, dvalue, ier)
            IF(ier.ne.0) THEN
                found = 0;
                RETURN
            ENDIF
        END SELECT
    ENDDO
    END SUBROUTINE

    SUBROUTINE cg_iRIC_Write_3DSolRealNode_f(Index, NX, NY, NZ, name, realarray, IER)
    IMPLICIT NONE
    ! INCLUDE "cgnswin_f.h"
    INCLUDE "cgnslib_f.h"
    INTEGER, INTENT(IN) :: Index
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: NX, NY, NZ
    Real*8, DIMENSION(:,:,:), INTENT(IN) :: realarray
    INTEGER, INTENT(out) :: ier

    INTEGER :: I,J,K, FSOL_ID, count
    CHARACTER(LEN=NAME_MAXLENGTH) :: arrayname
    REAL*4, ALLOCATABLE, DIMENSION(:) :: tmpvar1


    ALLOCATE(tmpvar1(NX*NY*NZ), STAT = IER)
    DO K = 1,NZ
        DO j= 1,NY
            DO i=1,NX
                count = i + ((j-1)*NX) + ((k-1)*NX*NY)
                tmpvar1(count) = realarray(i,j,k)
            END DO
        END DO
    ENDDO

    CALL cg_field_write_f(FILEID, BASE3DID, ZONE3DID, FLOWSOL2DID, RealSingle, name, tmpvar1, FSOL_ID, IER)
    DEALLOCATE(tmpvar1, STAT = IER)

    END SUBROUTINE cg_iRIC_Write_3DSolRealNode_f

    SUBROUTINE cg_iRIC_Write_SolRealNode_f(Index, NX, NY, name, realarray, IER)
    IMPLICIT NONE
    ! INCLUDE "cgnswin_f.h"
    INCLUDE "cgnslib_f.h"
    INTEGER, INTENT(IN) :: Index
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: NX, NY
    Real*8, DIMENSION(:,:), INTENT(IN) :: realarray
    INTEGER, INTENT(out) :: ier

    INTEGER :: I,J,FSOL_ID, count
    CHARACTER(LEN=NAME_MAXLENGTH) :: arrayname
    REAL*4, ALLOCATABLE, DIMENSION(:) :: tmpvar1


    ALLOCATE(tmpvar1(NX*NY), STAT = IER)

    DO j= 1,NY
        DO i=1,NX
            count = i + (j-1)*NX
            tmpvar1(count) = realarray(i,j)
        END DO
    END DO

    CALL cg_field_write_f(FILEID, BASEID, ZONEID, FLOWSOL2DID, RealSingle, name, tmpvar1, FSOL_ID, IER)
    DEALLOCATE(tmpvar1, STAT = IER)

    END SUBROUTINE cg_iRIC_Write_SolRealNode_f

    SUBROUTINE cg_iRIC_Write_3DSolIntegerNode_f(Index, NX, NY, NZ, name, intarray, IER)
    IMPLICIT NONE
    ! INCLUDE "cgnswin_f.h"
    INCLUDE "cgnslib_f.h"
    INTEGER, INTENT(IN) :: Index
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: NX, NY, NZ
    INTEGER, DIMENSION(:,:,:), INTENT(IN) :: intarray
    INTEGER, INTENT(out) :: ier

    INTEGER :: I,J,K, FSOL_ID, count
    CHARACTER(LEN=NAME_MAXLENGTH) :: arrayname
    INTEGER, ALLOCATABLE, DIMENSION(:) :: tmpvar1


    ALLOCATE(tmpvar1(NX*NY*NZ), STAT = IER)
    DO k = 1, NZ
        DO j= 1,NY
            DO i=1,NX
                count = i + ((j-1)*NX) + ((k-1)*NX*NY)
                tmpvar1(count) = intarray(i,j,K)
            END DO
        END DO
    END DO

    CALL cg_field_write_f(FILEID, BASEID, ZONEID, FLOWSOL2DID, Integer, name, tmpvar1, FSOL_ID, IER)
    DEALLOCATE(tmpvar1, STAT = IER)

    END SUBROUTINE cg_iRIC_Write_3DSolIntegerNode_f

    SUBROUTINE cg_iRIC_Write_SolIntegerNode_f(Index, NX, NY, name, intarray, IER)
    IMPLICIT NONE
    ! INCLUDE "cgnswin_f.h"
    INCLUDE "cgnslib_f.h"
    INTEGER, INTENT(IN) :: Index
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: NX, NY
    INTEGER, DIMENSION(:,:), INTENT(IN) :: intarray
    INTEGER, INTENT(out) :: ier

    INTEGER :: I,J,FSOL_ID, count
    CHARACTER(LEN=NAME_MAXLENGTH) :: arrayname
    INTEGER, ALLOCATABLE, DIMENSION(:) :: tmpvar1


    ALLOCATE(tmpvar1(NX*NY), STAT = IER)

    DO j= 1,NY
        DO i=1,NX
            count = i + (j-1)*NX
            tmpvar1(count) = intarray(i,j)
        END DO
    END DO

    CALL cg_field_write_f(FILEID, BASEID, ZONEID, FLOWSOL2DID, Integer, name, tmpvar1, FSOL_ID, IER)
    DEALLOCATE(tmpvar1, STAT = IER)

    END SUBROUTINE cg_iRIC_Write_SolIntegerNode_f

    SUBROUTINE cg_iRIC_Write_BaseIter(FID, sindex, Time, Disch)
    IMPLICIT NONE
    ! INCLUDE "cgnswin_f.h"
    INCLUDE "cgnslib_f.h"
    INTEGER, INTENT(IN) :: FID, sindex
    REAL*8, Dimension(:), INTENT(IN) :: Time, Disch
    INTEGER :: ier
    CALL cg_biter_write_f(FID, BASEID,'TimeIterValues',sindex,ier)
    call cg_goto_f(FID, BASEID,ier,'BaseIterativeData_t',1,'end')
    call cg_array_write_f('TimeValues',RealDouble,1,sindex,Time,ier)
    call cg_array_write_f('DischValues', RealDouble, 1, sindex, Disch, ier)
    END SUBROUTINE

    SUBROUTINE cg_iRIC_Write_ZoneIter(FID, sindex, SolNames)
    IMPLICIT NONE
    ! INCLUDE "cgnswin_f.h"
    INCLUDE "cgnslib_f.h"
    INTEGER, INTENT(IN) :: FID, sindex
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: SolNames
    INTEGER, DIMENSION(2) :: idata
    INTEGER :: ier, ier2, status, i
    CHARACTER*32, Allocatable, Dimension(:) :: tmpsolname
    CHARACTER(LEN = 32) :: tmp_flowsol1D, tmp_flowsol2D, tmp_flowsol3D
    CHARACTER(LEN = 14) :: flowsol1D, flowsol2D, flowsol3D
    CHARACTER(LEN = 5) :: tsindex

    CALL cg_ziter_write_f(FID, BASEID, ZONEID,'ZoneIterativeData',ier)
    call cg_goto_f(FID, BASEID,ier,'Zone_t',ZONEID,'ZoneIterativeData_t',1,'end')
    idata(1)=32
    idata(2)=sindex
    call cg_array_write_f('FlowSolutionPointers',Character,2,idata,solnames,ier2)
    END SUBROUTINE

    SUBROUTINE cg_iRIC_Write_3DBaseIter(FID, sindex, Time, Disch)
    IMPLICIT NONE
    ! INCLUDE "cgnswin_f.h"
    INCLUDE "cgnslib_f.h"
    INTEGER, INTENT(IN) :: FID, sindex
    REAL*8, Dimension(:), INTENT(IN) :: Time, Disch
    INTEGER :: ier
    CALL cg_biter_write_f(FID, BASE3DID,'TimeIterValues',sindex,ier)
    call cg_goto_f(FID, BASE3DID,ier,'BaseIterativeData_t',1,'end')
    call cg_array_write_f('TimeValues',RealDouble,1,sindex,Time,ier)
    call cg_array_write_f('DischValues', RealDouble, 1, sindex, Disch, ier)
    END SUBROUTINE

    SUBROUTINE cg_iRIC_Write_3DZoneIter(FID, sindex, SolNames)
    IMPLICIT NONE
    ! INCLUDE "cgnswin_f.h"
    INCLUDE "cgnslib_f.h"
    INTEGER, INTENT(IN) :: FID, sindex
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: SolNames
    INTEGER, DIMENSION(2) :: idata
    INTEGER :: ier, ier2, status, i
    CHARACTER*32, Allocatable, Dimension(:) :: tmpsolname
    CHARACTER(LEN = 32) :: tmp_flowsol1D, tmp_flowsol2D, tmp_flowsol3D
    CHARACTER(LEN = 14) :: flowsol1D, flowsol2D, flowsol3D
    CHARACTER(LEN = 5) :: tsindex

    CALL cg_ziter_write_f(FID, BASE3DID, ZONE3DID,'ZoneIterativeData',ier)
    call cg_goto_f(FID, BASE3DID,ier,'Zone_t',ZONE3DID,'ZoneIterativeData_t',1,'end')
    idata(1)=32
    idata(2)=sindex
    call cg_array_write_f('FlowSolutionPointers',Character,2,idata,solnames,ier2)
    END SUBROUTINE

    SUBROUTINE cg_iRIC_Read_Functional(name, length, x, y, IER)
    IMPLICIT NONE
    ! INCLUDE "cgnswin_f.h"
    INCLUDE "cgnslib_f.h"
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(OUT) :: length
    REAL*4, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: x, y
    INTEGER, INTENT(OUT) :: IER

    INTEGER :: narrays, paramfound, valuefound
    INTEGER :: datatype, dim
    INTEGER :: I
    CHARACTER(LEN = NAME_MAXLENGTH) :: arrayname
    INTEGER, DIMENSION(3) :: dimvec

    CALL local_gotoccchild(name, ier)
    if(ier.ne.0) RETURN
    ! under this UserDefinedData_t, array with name "Param", Value" will exist
    CALL cg_narrays_f(narrays, ier)
    IF(ier.ne.0) RETURN
    paramfound = 0
    valuefound = 0
    DO I = 1,narrays
        CALL cg_array_info_f(i, arrayname, datatype, dim, dimvec, IER)
        IF(IER.ne.0) RETURN
        SELECT CASE(TRIM(arrayname))
        CASE('Param')
            !the Value is found. dim must be 1, dimvec is the length of array, and datatype must be real.
            ! TODO: we have to check!
            paramfound = 1
            ALLOCATE(x(dimvec(1)), STAT = IER)
            if(IER.ne.0) RETURN
            length = dimvec(1)
            CALL cg_array_read_f(i, x, IER)
            IF(IER.ne.0) THEN
                length = 0
                RETURN
            ENDIF
        CASE('Value')
            !the Value is found. dim must be 1, dimvec is the length of array, and datatype must be real.
            ! TODO: we have to check!
            valuefound = 1
            ALLOCATE(y(dimvec(1)), STAT = IER)
            if(IER.ne.0) RETURN
            CALL cg_array_read_f(i, y, IER)
            IF(IER.ne.0) THEN
                length = 0
                RETURN
            ENDIF
        END SELECT
    ENDDO
    IF(paramfound.ne.1.and.valuefound.ne.1) THEN
        length = 0
        IER = 1
        RETURN
    ENDIF
    END SUBROUTINE

    SUBROUTINE CG_IRIC_READ_STRING(name, strvalue, IER)
    IMPLICIT NONE
    ! INCLUDE "cgnswin_f.h"
    INCLUDE "cgnslib_f.h"
    CHARACTER(LEN=*), INTENT(IN) :: name
    CHARACTER(LEN=*), INTENT(OUT) :: strvalue
    INTEGER, INTENT(OUT) :: IER
    CHARACTER(LEN=NAME_MAXLENGTH) :: arrayname
    INTEGER :: arrays, datatype, dim, dimvec, i, found
    CALL local_gotoccchild(name, ier)
    if(ier.ne.0) RETURN
    CALL cg_narrays_f(arrays, ier)
    IF(ier.ne.0) RETURN
    found = 0
    DO i = 1, arrays
        CALL cg_array_info_f(i, arrayname, datatype, dim, dimvec, ier)
        SELECT CASE(TRIM(arrayname))
        CASE('Value')
            found = 1
            CALL cg_array_read_f(i, strvalue, ier)
            IF(ier.ne.0) THEN
                found = 0
                RETURN
            ENDIF
        END SELECT
    ENDDO
    END SUBROUTINE

    !int cg_iRIC_Read_String(char* name, char* strvalue){
    !	char arrayname[NAME_MAXLENGTH];
    !	DataType_t datatype;
    !	int arrays;
    !	int dim;
    !	int dimvec;
    !	int ier;
    !	int i;
    !	int found;
    !
    !	ier = local_gotoccchild(name);
    !	if (ier != 0){return ier;}
    !	// under this UserDefinedData_t, array with name "Value" will exists.
    !	ier = cg_narrays(&arrays);
    !	if (ier != 0){return ier;}
    !	found = 0;
    !	for (i = 1; i <= arrays; ++i){
    !		ier = cg_array_info(i, arrayname, &datatype, &dim, &dimvec);
    !		if (strcmp(arrayname, "Value") == 0){
    !			// the Value is found. dim must be 1, dimvec must be 1, and dataype must be real.
    !			// TODO: we have to check!
    !			found = 1;
    !			cg_array_read(i, strvalue);
    !			strvalue[dimvec] = '\0';
    !		}
    !	}
    !	if (! found){return 1;}
    !	return 0;
    !}


    end module iriclibf