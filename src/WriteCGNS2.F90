    MODULE WriteCGNS2
    USE RivVarMod2
    USE RivVarWMod2
    USE RivVarVertMod2
    USE CalcCond2
    IMPLICIT NONE

    INTEGER :: tmpnx, tmpny


    CONTAINS

    SUBROUTINE Write_CGNS3D_Grid(rvo)
    IMPLICIT NONE
    INCLUDE "cgnslib_f.h"
    !        INCLUDE "cgnswin_f.h"
    type(rivvar), intent(in) :: rvo
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: tmpx3d, tmpy3d, tmpz3d
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: tmpint
    INTEGER :: ier
    CHARACTER(LEN=250) :: ZONENAME, BASENAME, USERNAME, SOLNAME
    INTEGER :: NUSER_DATA, USER_ITER, NZONES, ZONES_ITER
    INTEGER :: NBASES, BASES_ITER
    INTEGER :: BASE3DID, ZONE3DID
    INTEGER :: NSOLs, SOL_ITER, sindex
    INTEGER :: NGRIDS, GRIDS_ITER
    INTEGER :: F2DSOL_ITER, F3DSOL_ITER, F1DSOL_ITER, F2DSOL_EXT_ITER, FSOL_ITER
    INTEGER :: CELLDIM, PHYSDIM
    INTEGER :: GIndex, CIndex
    INTEGER, DIMENSION(3,3) :: isize
    INTEGER, DIMENSION(3) :: irmin, irmax
    CHARACTER(LEN = 250) :: name
    INTEGER, DIMENSION(2) :: idata
    integer :: ns2, nn, ns, nz
    integer :: fid
    ns2 = rvo%ns2
    ns = rvo%ns
    nn = rvo%nn
    nz = rvo%nz
    fid = rvo%cgnsfileid
    

    ALLOCATE(tmpx3d(ns2, nn, nz), STAT = ier)
    ALLOCATE(tmpy3d(ns2, nn, nz), STAT = ier)
    ALLOCATE(tmpz3d(ns2, nn, nz), STAT = ier)

    CALL GETXYZ3DOUT(rvo, tmpx3d, tmpy3d, tmpz3d)

    CALL cg_nbases_f(fid, NBASES, IER)
    IF(NBASES.EQ.1) THEN
        celldim = 3
        physdim = 3
        isize(1,1) = ns2
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

        CALL cg_grid_write_f(FID, BASE3DID, ZONE3DID, "GridCoordinates", GIndex, ier);
        CALL cg_coord_write_f(FID, BASE3DID, ZONE3DID, RealDouble, "CoordinateX", tmpx3d, CIndex, ier);
        CALL cg_coord_write_f(FID, BASE3DID, ZONE3DID, RealDouble, "CoordinateY", tmpy3d, CIndex, ier);
        CALL cg_coord_write_f(FID, BASE3DID, ZONE3DID, RealDouble, "CoordinateZ", tmpz3d, CIndex, ier);

        CAll cg_biter_write_f(FID, BASE3DID, 'BaseIterativeData', 1, ier);

    ELSE !iRIC3D already exists so delete BaseIterative, ZoneIterative and Solutions
        IF(NBASES.ne.2) THEN
            STOP 'Error in 3D CGNS FILE - Bases'
            !PAUSE
        ENDIF

        BASES_ITER = 2 !Hardwired should be iRIC3D
        ZONES_ITER = 1 !Hardwired should be iRICZone

        CALL cg_goto_f(FID, 2, ier, 'END')
        CALL cg_delete_node_f('BaseIterativeData', ier)
        CALL cg_nzones_f(FID, BASES_ITER, NZONES, IER)
        IF(NZONES.ne.1) THEN
            STOP 'Error in 3D CGNS File - Zones'
            !PAUSE
        ENDIF
        CALL cg_goto_f(FID, BASES_ITER, ier, 'Zone_t', ZONES_ITER, 'end')
        !Delete existing solutions
        CALL cg_nsols_F(FID, BASES_ITER, ZONES_ITER, nsols, ier);
        DO SOL_ITER = 1,nsols
            CALL cg_sol_info_f(FID, BASES_ITER, ZONES_ITER, 1, solname, sindex, ier)
            CALL cg_delete_node_f(TRIM(solname), ier);
        ENDDO
        !Delete existing solution grids
        CALL cg_ngrids_F(FID, BASES_ITER, ZONES_ITER, ngrids, ier);
        DO GRIDS_ITER = 1,ngrids-1
            CALL cg_grid_read_f(FID, BASES_ITER, ZONES_ITER, 2, solname, ier)
            CALL cg_delete_node_F(TRIM(solname), ier);
        ENDDO

        !Delete ZoneIterative Node
        CALL cg_goto_f(FID, BASES_ITER, ier, "Zone_t", ZONES_ITER, "end")
        CALL cg_delete_node_f('ZoneIterativeData', ier)

    ENDIF

    DEALLOCATE(tmpx3d, tmpy3d, tmpz3d, STAT = ier)

    END SUBROUTINE Write_CGNS3D_Grid

    SUBROUTINE Write_CGNS3D_FixedBed(rvo, rvto, solIndex, time, disch)
    IMPLICIT NONE
    INCLUDE "cgnslib_f.h"
    !        INCLUDE "cgnswin_f.h"
    type(rivvar), intent(in) :: rvo
    type(rivvartime), intent(inout) :: rvto
    INTEGER, INTENT(IN) :: solIndex
    REAL(KIND=mp), INTENT(IN) :: time, disch
    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:,:) :: tx,  ty, tz, tmpreal4, tmpreal5
    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:,:) :: tmpx3d, tmpy3d, tmpz3d
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: tmpint
    INTEGER :: ier, tmpibc, i,j,k
    CHARACTER(LEN = 14) :: flowsol3D
    CHARACTER(LEN = 26) :: gridsol3D
    CHARACTER(LEN = 5) :: tsindex
    CHARACTER(LEN = 32) :: tmp_flowsol3D, tmp_gridsol3D
    CHARACTER(LEN=250) :: ZONENAME, BASENAME, USERNAME, SOLNAME
    INTEGER :: NUSER_DATA, USER_ITER, NZONES, ZONES_ITER
    INTEGER :: NBASES, BASES_ITER
    INTEGER :: NSOL, SOL_ITER
    INTEGER :: F2DSOL_ITER, F3DSOL_ITER, F1DSOL_ITER, F2DSOL_EXT_ITER, FSOL_ITER
    INTEGER :: CELLDIM, PHYSDIM
    INTEGER :: GIndex, CIndex
    INTEGER, DIMENSION(3,3) :: isize
    INTEGER, DIMENSION(3) :: irmin, irmax
    CHARACTER(LEN = 250) :: name
    INTEGER, DIMENSION(2) :: idata
    INTEGER, DIMENSION(3) :: dimvec
    integer :: ns, ns2, nn, fid, nz
    ns = rvo%ns
    ns2 = rvo%ns2
    nn = rvo%nn
    fid = rvo%cgnsfileid
    nz = rvo%nz


    flowsol3D = '3DFlowSolution'
    gridsol3D = 'GridCoordinatesForSolution'
    WRITE(tsindex, '(I5)') solIndex
    tsindex = ADJUSTL(tsindex)
    tmp_flowsol3D = flowsol3D//TRIM(tsindex)
    tmp_flowsol3D = TRIM(tmp_flowsol3D)
    tmp_gridsol3D = gridsol3D//TRIM(tsindex)
    tmp_gridsol3D = TRIM(tmp_gridsol3D)

    rvto%SolNames3D(solIndex) = tmp_flowsol3D
    rvto%GridNames3D(solIndex) = tmp_gridsol3D
    rvto%TimeIncrements(solIndex) = time
    rvto%DischIncrements(solIndex) = disch

    ALLOCATE(tx(ns2, nn, nz), STAT = ier)
    ALLOCATE(ty(ns2, nn, nz), STAT = ier)
    ALLOCATE(tz(ns2, nn, nz), STAT = ier)
    ALLOCATE(tmpint(ns2, nn, nz), STAT = ier)

    !        // write base iterative data, with new solid to the current base.
    !	    ier = cg_biter_write(fileid, baseid, BINAME, solid);
    !	    if (ier != 0){return ier;}
    !	    ier = local_gotobaseiterative();
    !	    if (ier != 0){return ier;}
    !	    // write time.
    !	    dimVec[0] = solid;
    !	    ier = cg_array_write("TimeValues", RealDouble, 1, dimVec, timearray);



    CALL cg_nbases_f(FID, NBASES, IER)
    IF(NBASES.ne.2) THEN
        STOP 'ERROR in 3D base node'
        !PAUSE
    ENDIF
    ! Building iRIC3D node
    BASES_ITER = 2
    ZONES_ITER = 1

    !write base iter time data
    CALL cg_biter_write_f(FID, BASES_ITER, 'BaseIterativeData', solIndex, ier)
    CALL cg_goto_f(FID, BASES_ITER, ier, 'BaseIterativeData_t', 1, "end")
    CALL cg_array_write_f('TimeValues', RealDouble, 1, solindex, rvto%timeincrements, ier)

    CALL GETXYZ3DOUT(rvo, tx, ty, tz)
    CALL cg_grid_write_f(FID, BASES_ITER, ZONES_ITER, tmp_gridsol3D, GIndex, ier)
    CALL cg_goto_f(FID, BASES_ITER, ier, 'Zone_t', ZONES_ITER, 'GridCoordinates_t', solIndex+1, 'end')
    dimvec(1) = ns2
    dimvec(2) = nn
    dimvec(3) = nz
    CALL cg_array_write_f('CoordinateX', RealDouble, 3, dimVec, tx, ier)
    CALL cg_array_write_f('CoordinateY', RealDouble, 3, dimVec, ty, ier)
    CALL cg_array_write_f('CoordinateZ', RealDouble, 3, dimVec, tz, ier)

    Call GET3DVELOCITYOUT(rvo, tx, ty, tz)
    CALL cg_sol_write_f(FID, BASES_ITER, ZONES_ITER, tmp_flowsol3D, Vertex, F3DSOL_ITER, IER)
    CALL cg_field_write_f(FID, BASES_ITER, ZONES_ITER, F3DSOL_ITER, RealDouble, 'VelocityX', tx, FSOL_ITER, IER)
    CALL cg_field_write_f(FID, BASES_ITER, ZONES_ITER, F3DSOL_ITER, RealDouble, 'VelocityY', ty, FSOL_ITER, IER)
    CALL cg_field_write_f(FID, BASES_ITER, ZONES_ITER, F3DSOL_ITER, RealDouble, 'VelocityZ', tz, FSOL_ITER, IER)

    CALL GET3DVELOCITYSMAGOUT(rvo, tx)
    CALL cg_field_write_f(FID, BASES_ITER, ZONES_ITER, F3DSOL_ITER, RealDouble, 'SVelocity', tx, FSOL_ITER, IER)

    CALL GET3DVELOCITYNMAGOUT(rvo, tx)
    CALL cg_field_write_f(FID, BASES_ITER, ZONES_ITER, F3DSOL_ITER, RealDouble, 'NVelocity', tx, FSOL_ITER, IER)

    DO i = 1,ns2
        DO j = 1,nn
            DO k = 1,nz
                tmpibc = rvo%ibc(i,j)
                if(tmpibc > 0) then
                    tmpibc = -1
                endif
                tmpint(i,j,k) = tmpibc
            ENDDO
        ENDDO
    ENDDO
    CALL cg_field_write_f(FID, BASES_ITER, ZONES_ITER, F3DSOL_ITER, INTEGER, 'IBC', tmpint, FSOL_ITER, IER)

    CALL cg_ziter_write_f(FID,BASES_ITER,ZONES_ITER,'ZoneIterativeData',ier)
    call cg_goto_f(FID,BASES_ITER,ier,'Zone_t',ZONES_ITER,'ZoneIterativeData_t',1,'end')
    idata(1)=32
    idata(2)=solIndex
    call cg_array_write_f('FlowSolutionPointers',Character,2,idata,rvto%solnames3D,ier)
    call cg_array_write_f('GridCoordinatesPointers', Character, 2, idata, rvto%gridnames3D, ier)

    DEALLOCATE(tx, ty, tz, tmpint, STAT = ier)
    END SUBROUTINE Write_CGNS3D_FixedBed

    !SUBROUTINE Write_CGNS3D_MoveableBed(time, disch)
    !IMPLICIT NONE
    !REAL(KIND=mp), INTENT(IN) :: time, disch
    !REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:,:) :: tx,  ty, tz, tmpreal4, tmpreal5
    !REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:,:) :: tmpx3d, tmpy3d, tmpz3d
    !INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: tmpint
    !INTEGER :: ier
    !
    !!		CALL CG_IRIC_WRITE_SOL_TIME_F(time, ier)
    !ALLOCATE(tx(ns2, nn, nz), STAT = ier)
    !ALLOCATE(ty(ns2, nn, nz), STAT = ier)
    !ALLOCATE(tz(ns2, nn, nz), STAT = ier)
    !
    !call cg_iric_write_sol_time_f(time, ier);
    !
    !CALL GETXYZ3DOUT(tx, ty, tz)
    !call cg_iric_writegridcoord3d_f(tx, ty, tz, ier)
    !
    !Call GET3DVELOCITYOUT(tx, ty, tz)
    !call cg_iric_write_grid_real_node_f('VelocityX', tx, ier)
    !call cg_iric_write_grid_real_node_f('VelocityY', ty, ier)
    !call cg_iric_write_grid_real_node_f('VelocityZ', tz, ier)
    !
    !CALL GET3DVELOCITYSMAGOUT(tx)
    !call cg_iric_write_grid_real_node_f('SVelocity', tx, ier)
    !
    !CALL GET3DVELOCITYNMAGOUT(tx)
    !call cg_iric_write_grid_real_node_f('NVelocity', tx, ier)
    !
    !DEALLOCATE(tx, ty, tz, STAT = ier)
    !
    !END SUBROUTINE Write_CGNS3D_MoveableBed

    SUBROUTINE Write_CGNS2(rvo, cco, rvwo, time, disch)
    IMPLICIT NONE
    type(rivvar), intent(in) :: rvo
    type(calccond), intent(in) :: cco
    type(riv_w_var), intent(in) :: rvwo
    REAL(KIND=mp), INTENT(IN) :: time, disch
    INTEGER :: IER, iret
    INTEGER :: I,J, ns, nn, fid, ns2
    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:) :: tmpreal1,  tmpreal2
    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:,:) :: tmpreal1a,  tmpreal2a, tmpreal3, tmpreal4, tmpreal5
    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:,:) :: tmpx3d, tmpy3d, tmpz3d
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: tmpint
    INTEGER, ALLOCATABLE, DIMEnSION(:,:) :: tmp2dint
    !		CALL CG_IRIC_WRITE_SOL_TIME_F(time, ier)
    ns = rvo%ns
    ns2 = rvo%ns2
    nn = rvo%nn
    fid = rvo%cgnsfileid

    ALLOCATE(tmpreal1(ns2, nn), STAT = ier)
    ALLOCATE(tmpreal2(ns2, nn), STAT = ier)
    ALLOCATE(tmp2dint(ns2, nn), STAT = ier)
    
    IF(cco%CALCQUASI3D.and.cco%IO_3DOUTPUT) THEN
        CALL getXYOut(rvo, tmpreal1, tmpreal2)
        CALL CG_IRIC_WRITE_SOL_GRIDCOORD2D_F(tmpreal1,tmpreal2,IER)
    ENDIF

    CALL GETIBCOUT(rvo, tmp2dint)
    CALL cg_iRIC_Write_Sol_Integer_f("IBC", tmp2dint, IER)
    CALL GETFMIBCOUT(rvo, tmp2dint)
    CALL cg_iRIC_Write_Sol_Integer_f("FMIBC", tmp2dint, IER)

    CALL getVelocityOut(rvo, rvwo, tmpreal1, tmpreal2)
    CALL cg_iRIC_Write_Sol_Real_f("VelocityX", tmpreal1, IER)
    CALL cg_iRIC_Write_Sol_Real_f("VelocityY", tmpreal2, IER)

    CALL getDepthOut(rvo, tmpreal1)
    CALL cg_iRIC_Write_Sol_Real_f("Depth", tmpreal1, IER)

    CALL getWSEOut(rvo, tmpreal1)
    CALL cg_iRIC_Write_Sol_Real_f("WaterSurfaceElevation", tmpreal1, IER)

    CALL getElevationOut(rvo, tmpreal1)
    CALL cg_iRIC_Write_Sol_Real_f("Elevation", tmpreal1, IER)

    IF(cco%IO_VelSN) THEN
        CALL GETVELOCITYSNOUT(rvo, rvwo, tmpreal1, tmpreal2)
        CALL cg_iRIC_Write_Sol_Real_f("VelocityS", tmpreal1, IER)
        CALL cg_iRIC_Write_Sol_Real_f("VelocityN", tmpreal2, IER)
    ENDIF

    IF(cco%IO_UnitDisch) THEN
        CALL getUnitDischOut(rvo, rvwo, tmpreal1, tmpreal2)
        CALL cg_iRIC_Write_Sol_Real_f("UnitDischargeX", tmpreal1, IER)
        CALL cg_iRIC_Write_Sol_Real_f("UnitDischargeY", tmpreal2, IER)
    ENDIF

    IF(cco%IO_HArea) THEN
        CALL getHAreaOut(rvo, tmpreal1)
        CALL cg_iRIC_Write_Sol_Real_f("HabitatArea", tmpreal1, IER)
    ENDIF

    IF(cco%IO_InitVel) THEN
        CALL getInitVelOut(rvo, tmpreal1, tmpreal2)
        CALL cg_iRIC_Write_Sol_Real_f("VelocityInitS", tmpreal1, IER)
        CALL cg_iRIC_Write_Sol_Real_f("VelocityInitN", tmpreal2, IER)
    ENDIF

    IF(cco%IO_ShearXY) THEN
        CALL GETSHEARSTRESSOUT(rvo, tmpreal1, tmpreal2)
        CALL cg_iRIC_Write_Sol_Real_f("ShearStressX", tmpreal1, IER)
        CALL cg_iRIC_Write_Sol_Real_f("ShearStressY", tmpreal2, IER)
    ENDIF

    IF(cco%IO_ShearSN) THEN
        CALL GETSHEARSTRESSSNOUT(rvo, tmpreal1, tmpreal2)
        CALL cg_iRIC_Write_Sol_Real_f("ShearStressS", tmpreal1, IER)
        CALL cg_iRIC_Write_Sol_Real_f("ShearStressN", tmpreal2, IER)
    ENDIF

    IF(cco%IO_CD) THEN
        CALL getCDOut(rvo, tmpreal1)
        CALL cg_iRIC_Write_Sol_Real_f("Drag_Coefficient", tmpreal1, IER)
    ENDIF

    IF(cco%CALCCSED) THEN
        IF(cco%TRANSEQTYPE == 2) THEN
            !CALL GETSANDDEPTHOUT(rvo, tmpreal1)
            !CALL cg_iRIC_Write_Sol_Real_f("Sand_Depth", tmpreal1, IER)
            !CALL GETSANDFRACOUT(rvo, tmpreal1)
            !CALL cg_iRIC_Write_Sol_Real_f("Sand_Fraction", tmpreal1, IER)
            !CALL GETLSUBHOUT(rvo, tmpreal1)
            !CALL cg_iRIC_Write_Sol_Real_f("LSub", tmpreal1, IER)

        ENDIF
        CALL GetTransportRateOut(rvo, tmpreal1, tmpreal2)
        CALL cg_iRIC_Write_Sol_Real_f("SedFluxX", tmpreal1, IER)
        CALL cg_iRIC_Write_Sol_Real_f("SedFluxY", tmpreal2, IER)
    ENDIF

    IF(cco%IO_StressDiv) THEN
        CALL GetStressDivOut(rvo, tmpreal1)
        IF(cco%CALCCSED) THEN
            CALL cg_iRIC_Write_Sol_Real_f("Erosion Rate", tmpreal1, IER)
        ELSE
            CALL cg_iRIC_Write_Sol_Real_f("Shear Stress Divergence", tmpreal1, IER)
        ENDIF
    ENDIF

    IF(cco%CALCQUASI3D.and.cco%IO_HELIX) THEN
        CALL getHelixOut(rvo, tmpreal1)
        CALL cg_iRIC_Write_Sol_Real_f("Helix Strength", tmpreal1, IER)
        !            CALL getRSOut(tmpreal1)
        !            CALL cg_iRIC_Write_Sol_Real_f("RS", tmpreal1, IER)
        !            CALL getThetaOut(tmpreal1)
        !            CALL cg_iRIC_Write_Sol_Real_f("Theta", tmpreal1, IER)
    ENDIF


    DEALLOCATE(tmpreal1, STAT=ier)
    DEALLOCATE(tmpreal2, STAT=ier)
    DEALLOCATE(tmp2dint, STAT=ier)

    END SUBROUTINE




    SUBROUTINE GET3DGRIDZOUT(rvo, gridz)
    IMPLICIT NONE
    type(rivvar), intent(in) :: rvo
    REAL(KIND=mp), DIMENSION(:,:,:), INTENT(OUT) :: gridz
    INTEGER :: I,J,K
    DO k = 1,rvo%nz
        DO j= 1,rvo%nn
            DO i=1,rvo%ns2
                gridz(i,j,k) = (rvo%zz(i,j,k)/100.) + rvo%elevoffset
            END DO
        END DO
    END DO


    END SUBROUTINE

    SUBROUTINE GETIBCOUT(rvo, iibc)
    IMPLICIT NONE
    type(rivvar), intent(in) :: rvo
    INTEGER, DIMENSION(:,:), INTENT(OUT) :: iibc
    INTEGER :: I,J,K
    DO j= 1,rvo%nn
        DO i=1,rvo%ns2
            if(rvo%ibc(i,j).gt.0) then
                iibc(i,j) = -1
            else
                iibc(i,j) = rvo%ibc(i,j)
            endif
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETFMIBCOUT(rvo, iibc)
    IMPLICIT NONE
    type(rivvar), intent(in) :: rvo
    INTEGER, DIMENSION(:,:), INTENT(OUT) :: iibc
    INTEGER :: I,J,K
    DO j= 1,rvo%nn
        DO i=1,rvo%ns2
            !		        if(ibc(i,j).gt.0) then
            !		            iibc(i,j) = -1
            !		        else
            iibc(i,j) = rvo%ibc(i,j)
            !		        endif
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GET3DIBCOUT(rvo, iibc)
    IMPLICIT NONE
    type(rivvar), intent(in) :: rvo
    INTEGER, DIMENSION(:,:,:), INTENT(OUT) :: iibc
    INTEGER :: I,J,K
    DO k = 1,rvo%nz
        DO j= 1,rvo%nn
            DO i=1,rvo%ns2
                iibc(i,j,k) = rvo%ibc(i,j)
            END DO
        END DO
    END DO

    END SUBROUTINE


    SUBROUTINE GET3DVELOCITYOUT(rvo, ivelx, ively, ivelz)
    IMPLICIT NONE
    type(rivvar), intent(in) :: rvo
    REAL(KIND=mp), DIMENSION(:,:,:), INTENT(OUT) :: ivelx, ively, ivelz
    INTEGER :: I,J,K
    REAL(KIND=mp) :: rcos, rsin, xx, yy, ux, uy

    DO i=1,rvo%ns2
        rcos = cos(rvo%phirotation(i))
        rsin = sin(rvo%phirotation(i))
        DO j= 1,rvo%nn
            DO k = 1,rvo%nz

                ux = (rvo%uz(i,j,k)*rcos - rvo%vz(i,j,k)*rsin)/100.
                uy = (rvo%uz(i,j,k)*rsin + rvo%vz(i,j,k)*rcos)/100.
                xx = ux*rvo%fcos - uy*rvo%fsin
                yy = ux*rvo%fsin + uy*rvo%fcos
                ivelx(i,j,k) = xx
                ively(i,j,k) = yy
                ivelz(i,j,k) = 0.0

            END DO
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GET3DVELOCITYSOUT(rvo, ivelx, ively, ivelz)
    IMPLICIT NONE
    type(rivvar), intent(in) :: rvo
    REAL(KIND=mp), DIMENSION(:,:,:), INTENT(OUT) :: ivelx, ively, ivelz
    INTEGER :: i, j, k

    DO j= 1,rvo%nn
        DO i=1,rvo%ns2
            DO k = 1,rvo%nz
                ivelx(i,j,k) = rvo%uz(i,j,k)/100.
                ively(i,j,k) = 0.0
                ivelz(i,j,k) = 0.0
            ENDDO
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GET3DVELOCITYNOUT(rvo, ivelx, ively, ivelz)
    IMPLICIT NONE
    type(rivvar), intent(in) :: rvo
    REAL(KIND=mp), DIMENSION(:,:,:), INTENT(OUT) :: ivelx, ively, ivelz
    INTEGER :: i, j, k

    DO j= 1,rvo%nn
        DO i=1,rvo%ns2
            DO k=1,rvo%nz
                ivelx(i,j,k) = 0.0
                ively(i,j,k) = rvo%vz(i,j,k)/100.
                ivelz(i,j,k) = 0.0
            ENDDO
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GET3DVELOCITYNMAGOUT(rvo, ively)
    IMPLICIT NONE
    type(rivvar), intent(in) :: rvo
    REAL(KIND=mp), DIMENSION(:,:,:), INTENT(OUT) :: ively
    INTEGER :: i, j, k

    DO j= 1,rvo%nn
        DO i=1,rvo%ns2
            DO k=1,rvo%nz
                ively(i,j,k) = rvo%vz(i,j,k)/100.
            ENDDO
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GET3DVELOCITYSMAGOUT(rvo, ivelx)
    IMPLICIT NONE
    type(rivvar), intent(in) :: rvo
    REAL(KIND=mp), DIMENSION(:,:,:), INTENT(OUT) :: ivelx
    INTEGER :: i, j, k

    DO j= 1,rvo%nn
        DO i=1,rvo%ns2
            DO k=1,rvo%nz
                ivelx(i,j,k) = rvo%uz(i,j,k)/100.
            ENDDO
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETVELOCITYSNOUT(rvo, rvwo, ivelx, ively)
    IMPLICIT NONE
    type(rivvar), intent(in) :: rvo
    type(riv_w_var), intent(in) :: rvwo
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: ivelx, ively
    INTEGER :: i, j

    DO j= 1,rvo%nn
        DO i=1,rvo%ns2
            ivelx(i,j) = rvo%u(i,j)/100.
            ively(i,j) = rvwo%vout(i,j)/100.
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETSHEARSTRESSSNOUT(rvo, issx, issy)
    IMPLICIT NONE
    type(rivvar), intent(in) :: rvo
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: issx, issy
    INTEGER :: i, j

    DO j= 1,rvo%nn
        DO i=1,rvo%ns2
            issx(i,j) = rvo%taus(i,j)/10.
            issy(i,j) = rvo%taun(i,j)/10.
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETSHEARSTRESSOUT(rvo, SSXOut, SSYOut)
    IMPLICIT NONE
    type(rivvar), intent(in) :: rvo
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: SSXOut, SSYOut
    INTEGER :: i,j,count
    REAL(KIND=mp) :: rcos, rsin, xx, yy, ux, uy


    DO j= 1,rvo%nn
        DO i=1,rvo%ns2
            rcos = cos(rvo%phirotation(i))
            rsin = sin(rvo%phirotation(i))
            !					count = ((i-1)*nn)+j
            count = i + (j-1)*rvo%ns

            ux = (rvo%taus(i,j)*rcos - rvo%taun(i,j)*rsin)/10.
            uy = (rvo%taus(i,j)*rsin + rvo%taun(i,j)*rcos)/10.
            xx = ux*rvo%fcos - uy*rvo%fsin
            yy = ux*rvo%fsin + uy*rvo%fcos
            SSXOut(i,j) = xx
            SSYOut(i,j) = yy
            !				tmpvar3(count) = u(i,j)/100.
            !				tmpvar4(count) = vout(i,j)/100.
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETVELOCITYOUT(rvo, rvwo, velXOut, velYOut)
    IMPLICIT NONE
    type(rivvar), intent(in) :: rvo
    type(riv_w_var), intent(in) :: rvwo
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: velXOut, velYOut
    INTEGER :: i,j,count
    REAL(KIND=mp) :: rcos, rsin, xx, yy, ux, uy


    DO j= 1,rvo%nn
        DO i=1,rvo%ns2
            rcos = cos(rvo%phirotation(i))
            rsin = sin(rvo%phirotation(i))
            !					count = ((i-1)*nn)+j
            count = i + (j-1)*rvo%ns

            !				    ux = (u(i,j)*rcos - vout(i,j)*rsin)/100.
            !				    uy = (u(i,j)*rsin + vout(i,j)*rcos)/100.
            !if(uvinterp == 0) then
            ux = (rvo%u(i,j)*rcos - rvwo%vout(i,j)*rsin)/100.
            uy = (rvo%u(i,j)*rsin + rvwo%vout(i,j)*rcos)/100.
            !else
            !  ux = (uout(i,j)*rcos - vout(i,j)*rsin)/100.
            !  uy = (uout(i,j)*rsin + vout(i,j)*rcos)/100.
            !endif
            xx = ux*rvo%fcos - uy*rvo%fsin
            yy = ux*rvo%fsin + uy*rvo%fcos
            velXOut(i,j) = xx
            velYOut(i,j) = yy
            !				tmpvar3(count) = u(i,j)/100.
            !				tmpvar4(count) = vout(i,j)/100.
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETTransportRateOUT(rvo, TRXOut, TRYOut)
    IMPLICIT NONE
    type(rivvar), intent(in) :: rvo
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: TRXOut, TRYOut
    INTEGER :: i,j,count
    REAL(KIND=mp) :: rcos, rsin, xx, yy, ux, uy


    DO j= 1,rvo%nn
        DO i=1,rvo%ns2
            rcos = cos(rvo%phirotation(i))
            rsin = sin(rvo%phirotation(i))
            !					count = ((i-1)*nn)+j
            count = i + (j-1)*rvo%ns

            ux = (rvo%qs(i,j)*rcos - rvo%qn(i,j)*rsin)
            uy = (rvo%qs(i,j)*rsin + rvo%qn(i,j)*rcos)
            xx = ux*rvo%fcos - uy*rvo%fsin
            yy = ux*rvo%fsin + uy*rvo%fcos
            TRXOut(i,j) = xx/10000.
            TRYOut(i,j) = yy/10000.
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETINITVELOUT(rvo, ivelx, ively)
    IMPLICIT NONE
    type(rivvar), intent(in) :: rvo
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: ivelx, ively
    INTEGER :: i, j

    DO j= 1,rvo%nn
        DO i=1,rvo%ns2
            ivelx(i,j) = rvo%u(i,j)/100.
            ively(i,j) = rvo%v(i,j)/100.
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETUNITDISCHOUT(rvo, rvwo, ivelx, ively)
    IMPLICIT NONE
    type(rivvar), intent(in) :: rvo
    type(riv_w_var), intent(in) :: rvwo
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: ivelx, ively
    INTEGER :: i, j, count
    REAL(KIND=mp) :: rcos, rsin, xx, yy, ux, uy
    DO j= 1,rvo%nn
        DO i=1,rvo%ns2
            rcos = cos(rvo%phirotation(i))
            rsin = sin(rvo%phirotation(i))
            !					count = ((i-1)*nn)+j
            count = i + (j-1)*rvo%ns

            ux = (rvo%u(i,j)*rcos - rvwo%vout(i,j)*rsin)/100.
            uy = (rvo%u(i,j)*rsin + rvwo%vout(i,j)*rcos)/100.
            xx = ux*rvo%fcos - uy*rvo%fsin
            yy = ux*rvo%fsin + uy*rvo%fcos
            ivelx(i,j) = xx
            ively(i,j) = yy
            !				tmpvar3(count) = u(i,j)/100.
            !				tmpvar4(count) = vout(i,j)/100.
        END DO
    END DO
    DO j= 1,rvo%nn
        DO i=1,rvo%ns2
            ivelx(i,j) = ivelx(i,j) * rvo%hl(i,j)/100.
            ively(i,j) = ively(i,j) * rvo%hl(i,j)/100.
        END DO
    END DO

    END SUBROUTINE

    !SUBROUTINE GETSANDDEPTHOUT(sd)
    !IMPLICIT NONE
    !type(rivvar), intent(in) :: rvo
    !REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: sd
    !INTEGER :: i, j
    !
    !DO j= 1,nn
    !    DO i=1,ns2
    !        sd(i,j) = hfine(i,j)/100.
    !    END DO
    !END DO
    !
    !END SUBROUTINE

    !SUBROUTINE GETLSUBHOUT(tlsub)
    !IMPLICIT NONE
    !type(rivvar), intent(in) :: rvo
    !REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: tlsub
    !INTEGER :: i, j
    !
    !DO j= 1,nn
    !    DO i=1,ns2
    !        tlsub(i,j) = (vsandsub(i,j)/(0.3*harea(i,j)))/100.
    !    END DO
    !END DO
    !
    !END SUBROUTINE


    SUBROUTINE getRSOut(rvo, rvvo, irs)
    IMPLICIT NONE
    type(rivvar), intent(in) :: rvo
    type(rivvarvert), intent(in) :: rvvo
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: irs
    INTEGER :: i, j

    DO j= 1,rvo%nn
        DO i=1,rvo%ns2
            irs(i,j) = rvvo%rs(i,j)/100.
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE getThetaOut(rvo, rvvo, itheta)
    IMPLICIT NONE
    type(rivvar), intent(in) :: rvo
    type(rivvarvert), intent(in) :: rvvo
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: itheta
    INTEGER :: i, j

    DO j= 1,rvo%nn
        DO i=1,rvo%ns2
            itheta(i,j) = rvvo%theta(i,j)/100.
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETSANDFRACOUT(rvo, sf)
    IMPLICIT NONE
    type(rivvar), intent(in) :: rvo
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: sf
    INTEGER :: i, j

    DO j= 1,rvo%nn
        DO i=1,rvo%ns2
            sf(i,j) = rvo%Fracs(i,j)
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETCDOUT(rvo, tcd)
    IMPLICIT NONE
    type(rivvar), intent(in) :: rvo
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: tcd
    INTEGER :: i, j

    DO j= 1,rvo%nn
        DO i=1,rvo%ns2
            tcd(i,j) = rvo%totcd(i,j)
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETDEPTHOUT(rvo, depth)
    IMPLICIT NONE
    type(rivvar), intent(in) :: rvo
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: depth
    INTEGER :: i, j

    DO j= 1,rvo%nn
        DO i=1,rvo%ns2
            depth(i,j) = rvo%hl(i,j)/100.
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETHAREAOUT(rvo, htArea)
    IMPLICIT NONE
    type(rivvar), intent(in) :: rvo
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: htArea
    INTEGER :: i, j

    DO j= 1,rvo%nn
        DO i=1,rvo%ns2
            htArea(i,j) = rvo%harea(i,j)/10000.
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETWSEOUT(rvo, wse)
    IMPLICIT NONE
    type(rivvar), intent(in) :: rvo
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: wse
    INTEGER :: i, j

    DO j= 1,rvo%nn
        DO i=1,rvo%ns2
            wse(i,j) = (rvo%e(i,j)/100.) + rvo%elevoffset
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETELEVATIONOUT(rvo, elev)
    IMPLICIT NONE
    type(rivvar), intent(in) :: rvo
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: elev
    INTEGER :: i, j

    DO j= 1,rvo%nn
        DO i=1,rvo%ns2
            elev(i,j) = (rvo%eta(i,j)/100.) + rvo%elevoffset
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETXYOUT(rvo, tx, ty)
    IMPLICIT NONE
    type(rivvar), intent(in) :: rvo
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: tx, ty
    REAL(KIND=mp) :: xx, yy, ux, uy
    double precision :: rcos, rsin
    INTEGER :: i, j
    DO j= 1,rvo%nn
        DO i=1,rvo%ns2
            rcos = cos(rvo%phirotation(i))
            rsin = sin(rvo%phirotation(i))
            !ux = rvo%x(i,j)*rcos - rvo%y(i,j)*rsin
            !uy = rvo%x(i,j)*rsin + rvo%y(i,j)*rcos

            xx = rvo%x(i,j)*rvo%fcos - rvo%y(i,j)*rvo%fsin
            yy = rvo%x(i,j)*rvo%fsin + rvo%y(i,j)*rvo%fcos
            tx(i,j) = (xx + rvo%xshift)/100.
            ty(i,j) = (yy + rvo%yshift)/100.

        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETXYZ3DOUT(rvo, tx, ty, tz)
    IMPLICIT NONE
    type(rivvar), intent(in) :: rvo
    Double Precision, DIMENSION(:,:,:), INTENT(OUT) :: tx, ty, tz
    INTEGER :: i, j, k
    Double Precision :: xx, yy, ux, uy
    double precision :: rcos, rsin

    DO k = 1,rvo%nz
        DO j=1,rvo%nn
            DO i= 1,rvo%ns2
                rcos = cos(rvo%phirotation(i))
                rsin = sin(rvo%phirotation(i))

                !ux = rvo%x(i,j)*rcos - rvo%y(i,j)*rsin
                !uy = rvo%x(i,j)*rsin + rvo%y(i,j)*rcos

                xx = rvo%x(i,j)*rvo%fcos - rvo%y(i,j)*rvo%fsin
                yy = rvo%x(i,j)*rvo%fsin + rvo%y(i,j)*rvo%fcos

                tx(i,j,k) = (xx + rvo%xshift)/100.
                ty(i,j,k) = (yy + rvo%yshift)/100.
                tz(i,j,k) = (rvo%zz(i,j,k)/100.) + rvo%elevoffset
            END DO
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE getHelixOut(rvo, helix)
    IMPLICIT NONE
    type(rivvar), intent(in) :: rvo
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: helix
    REAL(KIND=mp) :: xx, yy
    INTEGER :: i, j
    DO j= 1,rvo%nn
        DO i=1,rvo%ns2
            xx = atan2(rvo%taun(i,j), rvo%taus(i,j))*(360./(2.*3.14159))
            yy = atan2(rvo%vz(i,j,rvo%nz), rvo%uz(i,j, rvo%nz))*(360./(2.*3.14159))
            !			    if(vz(i,j,2).eq.0.0.or.uz(i,j,2).eq.0) then
            !			        helix(i,j) = 0
            !			    else
            helix(i,j) = yy-xx
            if(helix(i,j) > 180) then
                helix(i,j) = helix(i,j)-360.
            elseif(helix(i,j) < -180) then
                helix(i,j) = helix(i,j)+360.
            endif
            !		        Endif
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GetStressDivOut(rvo, stressdiv)
    IMPLICIT NONE
    type(rivvar), intent(in) :: rvo
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: stressdiv
    INTEGER :: i, j
    DO j= 1,rvo%nn
        DO i=1,rvo%ns2
            stressdiv(i,j) = rvo%con(i,j)
        END DO
    END DO

    END SUBROUTINE

    END MODULE
