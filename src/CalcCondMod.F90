    MODULE CalcCond
    USE RivVarMod
    USE RivVarTimeMod
    IMPLICIT NONE

    REAL(KIND=mp) :: erelax, urelax, arelax
    REAL(KIND=mp) :: dt, q
    REAL(KIND=mp) :: vdc, vac
    REAL(KIND=mp) :: wselev, wsupelev, wsslope
    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:) :: vbcdist, vbcvel, vbcang

    LOGICAL :: HotStart
    LOGICAL :: FLUMEBNDRY


    INTEGER :: itm, interitm, iplinc
    INTEGER :: debugstop, DbgTimeStep, DbgIterNum
    INTEGER :: vbc
    INTEGER :: solType

    !Variables for iteration output
    LOGICAL :: IterationOut
    INTEGER :: IterPlOut

    INTEGER :: maxInterIterMult

    LOGICAL :: IO_IBC, IO_VelXY, IO_VelSN, IO_InitVel, IO_UnitDisch
    LOGICAL :: IO_ShearXY, IO_ShearSN, IO_WSE, IO_Helix, IO_3DOutput
    LOGICAL :: IO_Depth, IO_StressDiv, IO_CD, IO_Elev
    LOGICAL :: IO_Vort, IO_KEG, IO_VelStrain
    LOGICAL :: IO_HArea

    !Variables for HOTSTART
    CHARACTER(LEN=250) :: CGNSHSFile
    INTEGER :: SolnIndex

    INTEGER :: cdtype, roughnesstype, wstype
    REAL(KIND=mp) :: constcd, ONEDCD
    REAL(KIND=mp) :: cdmin, cdmax

    INTEGER, PUBLIC :: nsext
    INTEGER :: ShowGridExt, vbcds

    REAL :: nsextslope

    !Wetting and Drying
    INTEGER :: dryType, hiterInterval, hiterstop
    LOGICAL :: hcalcwetting
    REAL(KIND=mp) :: hmin, hwmin

    !LEV
    INTEGER :: LEVType, LEVChangeIter, LEVBegIter, LEVEndIter
    REAL(KIND=mp) :: evc, startLEV, endLEV

    !Quasi3D
    LOGICAL :: CALCQUASI3D, CALCQUASI3DRS
    REAL(KIND=mp) :: MinRS

    !Sediment Transport SEDSMOOWGHT
    INTEGER :: SEDSMOO, RSSMOO, SEDSMOOWGHT, TRANSEQTYPE
    INTEGER :: SEDBCNODE, CALCGRAVCORR, GRAVCORRTYPE, SEDEQNODE
    REAL(KIND=mp) :: HD, WD, DIN, BCFRACTION, SUBANGLEOFREPOSE
    REAL(KIND=mp) :: TSFracDepth, GRAVFLATBEDCORRCOEF
    !REAL(KIND=mp) :: SEDEQMULT
    LOGICAL :: CALCCSED, CALCSEDAUTO

    !Daniele Tonina Wilcock-Kenworthy transport function
    REAL(KIND=mp) :: taustar_rg0, taustar_rg1, taustar_rs1, alpha, AK, Chi
    REAL(KIND=mp) :: phi_prime, tmphmax, tmpabh, tmpbbh, tmplsub, tmplsubactdepth
    REAL(KIND=mp) :: Dsand, Dg, Z0Min, Z0Max, WK_RoughShapeParam, Fsmin
    INTEGER :: WK_RoughType
    LOGICAL :: WKConstBndry

    INTEGER :: uvinterp

    ! For Hot Start
    integer :: i_re_flag_i, i_re_flag_o, n_rest, i_tmp_count
    real(8) :: opt_tmp(0:9)

    character(len = strMax) :: tmp_file_o(0:9), tmp_caption(0:9) &
        ,tmp_file_i, tmp_dummy, tmp_pass

    CONTAINS

    SUBROUTINE ALLOC_VELBC(size)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: size
    INTEGER :: status
    ALLOCATE(vbcdist(size), STAT = status)
    ALLOCATE(vbcvel(size), STAT = status)
    ALLOCATE(vbcang(size), STAT = status)
    END SUBROUTINE

    SUBROUTINE DEALLOC_VELBC()
    IMPLICIT NONE
    INTEGER :: status
    DEALLOCATE(vbcdist, STAT = status)
    DEALLOCATE(vbcvel, STAT = status)
    DEALLOCATE(vbcang, STAT = status)
    END SUBROUTINE


    SUBROUTINE CGNS2_Read_CC_ForAlloc(IER)
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: IER
    INTEGER :: tmpint

    CALL CG_IRIC_READ_INTEGER(fid, 'FM_SolAttNSExt', nsext, ier)
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_SolAttItm', itm, ier)
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_Quasi3D_NZ', nz, ier)
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_Quasi3D', tmpint, ier)
    IF(tmpint == 1) THEN
        CALCQUASI3D = .TRUE.
    ELSE
        CALCQUASI3D = .FALSE.
    ENDIF


    ENDSUBROUTINE

    SUBROUTINE CGNS2_Read_CalcCondition(IER)
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: IER

    INTEGER :: status, i, j, count, ierror
    REAL(kind = mp) :: rwork(1), ttt
    REAL(kind = mp) :: tmpreal
    INTEGER, DIMENSION(1) :: tint
    REAL(kind = mp), DIMENSION(:), ALLOCATABLE :: treal
    INTEGER :: iwork(1), tmpint
    REAL(KIND=mp), DIMENSION(:), ALLOCATABLE :: xtmp, ytmp, ztmp
    REAL(KIND=mp), DIMENSION(:), ALLOCATABLE :: dischT, dischQ, stageT, StageQ, stageH
    INTEGER :: sizeDischTQ, sizeStageTQ, sizeStageHQ

    INTEGER :: ii, iii, count2
    !For hotstart
    i_tmp_count = 0
    do ii=0,9
        tmp_dummy = "tmp0.d"
        write(tmp_dummy(4:4),'(i1)') ii
        tmp_file_o(ii) = tmp_dummy
        tmp_dummy = "outtime_tmp0"
        write(tmp_dummy(12:12),'(i1)') ii
        tmp_caption(ii) = tmp_dummy
    end do

    !These are hard wired
    IO_IBC= .TRUE.
    IO_VELXY= .TRUE.
    IO_WSE= .TRUE.
    IO_DEPTH= .TRUE.
    IO_ELEV= .TRUE.

    CALL CG_IRIC_READ_INTEGER(fid, 'FM_SolWrtVelSN', tmpint, ier)
    if(tmpint == 0) then
        IO_VELSN = .FALSE.
    else
        IO_VELSN = .TRUE.
    endif
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_SolWrtShearXY', tmpint, ier)
    if(tmpint == 0) then
        IO_SHEARXY = .FALSE.
    else
        IO_SHEARXY= .TRUE.
    endif
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_SolWrtShearSN', tmpint, ier)
    if(tmpint == 0) then
        IO_SHEARSN = .FALSE.
    else
        IO_SHEARSN= .TRUE.
    endif
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_SolWrtShearDiv', tmpint, ier)
    if(tmpint == 0) then
        IO_STRESSDIV = .FALSE.
    else
        IO_STRESSDIV= .TRUE.
    endif
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_SolWrtCD', tmpint, ier)
    if(tmpint == 0) then
        IO_CD = .FALSE.
    else
        IO_CD= .TRUE.
    endif
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_SolWrtUnitDischarge', tmpint, ier)
    if(tmpint == 0) then
        IO_UnitDisch = .FALSE.
    else
        IO_UnitDisch = .TRUE.
    endif
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_SolWrtHArea', tmpint, ier)
    if(tmpint == 0) then
        IO_HArea = .FALSE.
    else
        IO_HArea = .TRUE.
    endif
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_SolWrtInitVel', tmpint, ier)
    if(tmpint == 0) then
        IO_InitVel = .FALSE.
    else
        IO_InitVel = .TRUE.
    endif
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_SolWrt3DOutput', tmpint, ier)
    if(tmpint == 0) then
        IO_3DOUTPUT = .FALSE.
    else
        IO_3DOUTPUT = .TRUE.
    endif
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_SolWrtHelix', tmpint, ier)
    if(tmpint == 0) then
        IO_HELIX = .FALSE.
    else
        IO_HELIX = .TRUE.
    endif




    IO_VORT= .FALSE.
    IO_KEG= .FALSE.
    IO_VELSTRAIN= .FALSE.

    !CALL CG_IRIC_GOTOCC_F(FID, IER)

    CALL cg_iRIC_Read_Integer(fid, 'FM_SolAttType', solType, ier)

    CALL cg_iRIC_Read_Real('FM_SolAttERlx', erelax, ier)
    CALL cg_iRIC_Read_Real('FM_SolAttURlx', urelax, ier)
    CALL cg_iRIC_Read_Real('FM_SolAttARlx', arelax, ier)

    CALL CG_IRIC_READ_INTEGER(fid, 'FM_SolAttItm', itm, ier)

    CALL cg_iRIC_Read_Real('FM_SolAttDT', dt, ier)

    CALL CG_IRIC_READ_INTEGER(fid, 'FM_SolAttInterItm', interitm, ier)

    CALL CG_IRIC_READ_INTEGER(fid, 'FM_SolAttPlinc', iplinc, ier)

    !NEXT TWO DEFINED IN RivVarTimeMod
    CALL cg_iRIC_Read_Real('FM_SolAttVarDischSTime', vardischstarttime, ier)
    CALL cg_iRIC_Read_Real('FM_SolAttVarDischETime', vardischendtime, ier)

    CALL CG_IRIC_READ_INTEGER(fid, 'FM_SolAttFlumeBndry', tmpint, ier)
    IF(tmpint == 1) THEN
        FLUMEBNDRY = .TRUE.
    ELSE
        FLUMEBNDRY = .FALSE.
    ENDIF

    CALL CG_IRIC_READ_INTEGER(fid, 'FM_SolAttDbgStop', debugstop, ier)
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_SolAttDbgTimeStep', DbgTimeStep, ier)
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_SolAttDbgIterNum', DbgIterNum, ier)

    CALL CG_IRIC_READ_INTEGER(fid, 'FM_SolAttIterOut', tmpint, ier)
    IF(tmpint == 1) THEN
        IterationOut = .TRUE.
    ELSE
        IterationOut = .FALSE.
    ENDIF
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_SolAttIncPlOut', IterPlOut, ier)

    CALL CG_IRIC_READ_INTEGER(fid, 'FM_SolAttMaxInterIter', maxInterIterMult, ier)

    CALL cg_iRIC_Read_Real('FM_HydAttQ', q, ier)
    q = q*1.e6
    !NEXT TWO DEFINED IN RivVarTimeMod
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_HydAttVarDischType', VarDischType, ier)
    IF(VarDischType.eq.1)THEN
        !dischT and dischQ are allocated in _READ_FUNCTIONAL
        CALL CG_IRIC_READ_FUNCTIONALSIZE('FM_HydAttVarDischarge', sizeDischTQ, IER)
        ALLOCATE(dischT(sizeDischTQ), dischQ(sizeDischTQ), STAT=IER)
        CALL CG_IRIC_READ_FUNCTIONAL('FM_HydAttVarDischarge', dischT, dischQ, ier)
        tint(1) = sizeDischTQ
        CALL SetTimeSeriesNumPts(1, tint)
        ALLOCATE(treal(2*sizeDischTQ), STAT = ier)
        count2=0
        DO ii = 1, sizeDischTQ
            count2 = (ii-1)*2
            treal(count2+1) = dischT(ii)
            treal(count2+2) = dischQ(ii)
        ENDDO
        Call SetTimeSeries((2*sizeDischTQ), treal)
        tint(1) = 0
        CALL SetTimeSeriesType(1,tint)
        DEALLOCATE(treal, STAT = ier)
        DEALLOCATE(dischT, STAT = ier)
        DEALLOCATE(dischQ, STAT = ier)
    ENDIF

    CALL CG_IRIC_READ_INTEGER(fid, 'FM_HydAttVarStgType', VarStageType, ier)
    IF(VarStageType.eq.1)THEN
        CALL CG_IRIC_READ_FUNCTIONALSIZE('FM_HydAttVarStageTS', sizeStageTQ, IER)
        ALLOCATE(stageT(sizeStageTQ), stageH(sizeStageTQ), STAT=IER)
        CALL CG_IRIC_READ_FUNCTIONAL('FM_HydAttVarStageTS', stageT, stageH, ier)
        tint(1) = sizeStageTQ
        CALL SetRatingCurveNumPts(1, tint)
        ALLOCATE(treal(2*sizeStageTQ), STAT = ier)
        count2=0
        DO ii = 1, sizeStageTQ
            count2 = (ii-1)*2
            treal(count2+1) = stageT(ii)
            treal(count2+2) = (stageH(ii)-elevOffset)*100.
        ENDDO
        Call SetRatingCurves((2*sizeStageTQ), treal)
        tint(1) = 0
        CALL SetRatingCurveType(1,tint)
        DEALLOCATE(treal, STAT = ier)
        DEALLOCATE(stageT, STAT = ier)
        DEALLOCATE(stageH, STAT = ier)

    ELSE IF(VarStageType.eq.2) THEN
        CALL CG_IRIC_READ_FUNCTIONALSIZE('FM_HydAttVarStageRC', sizeStageHQ, IER)
        ALLOCATE(stageQ(sizeStageHQ), stageH(sizeStageHQ), STAT=IER)
        CALL CG_IRIC_READ_FUNCTIONAL('FM_HydAttVarStageRC', stageQ, stageH, ier)
        tint(1) = sizeStageHQ
        CALL SetRatingCurveNumPts(1, tint)
        ALLOCATE(treal(2*sizeStageHQ), STAT = ier)
        count2=0
        DO ii = 1, sizeStageHQ
            count2 = (ii-1)*2
            treal(count2+1) = stageQ(ii)
            treal(count2+2) = (stageH(ii) - elevoffset)*100
        ENDDO
        Call SetRatingCurves((2*sizeStageHQ), treal)
        tint(1) = 1
        CALL SetRatingCurveType(1,tint)
        DEALLOCATE(treal, STAT = ier)
        DEALLOCATE(stageH, STAT = ier)
        DEALLOCATE(stageQ, STAT = ier)

    ENDIF

    CALL cg_iRIC_Read_Real('FM_HydAttVelDepthCoef', vdc, ier)
    CALL cg_iRIC_Read_Real('FM_HydAttVelAngleCoef', vac, ier)
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_HydAttVelBC', vbc, ier)
    IF(vbc == 1) THEN
        CALL CG_IRIC_READ_FUNCTIONALSIZE('FM_HydAttVarVelBC', tmpint, IER)
        Allocate(xtmp(tmpint), ytmp(tmpint), ztmp(tmpint), STAT=IER)
        CALL ALLOC_VELBC(tmpint)
        !   CALL CG_IRIC_READ_FUNCTIONAL_F('FM_HydAttVarVelBC', xtmp, ytmp, ztmp, ier)
        !   CALL CG_IRIC_READ_FUNCTIONAL_F('FM_HydAttVarVelBC', xtmp, ytmp, ier)
        CALL CG_IRIC_READ_FUNCTIONALWITHNAME("FM_HydAttVarVelBC","Normalized_Distance",xtmp,ier)
        CALL CG_IRIC_READ_FUNCTIONALWITHNAME("FM_HydAttVarVelBC","Normalized_Velocity",ytmp,ier)
        CALL CG_IRIC_READ_FUNCTIONALWITHNAME("FM_HydAttVarVelBC","Angle",ztmp,ier)
        IF(ier.eq.0) THEN
            DO i = 1, tmpint
                vbcdist(i) = xtmp(i)
                vbcvel(i) = ytmp(i)
                vbcang(i) = ztmp(i)
            ENDDO
        ENDIF
        DEALLOCATE(xtmp, STAT=ier)
        DEALLOCATE(ytmp, STAT=ier)
        DEALLOCATE(ztmp, STAT=ier)
    ENDIF

    CALL cg_iRIC_Read_Real('FM_HydAttWS', wselev, ier)
    wselev=(wselev - elevoffset)*100.
    CALL cg_iRIC_Read_Real('FM_HydAttWS2', wsupelev, ier)
    wsupelev=(wsupelev - elevoffset)*100.
    CALL cg_iRIC_Read_Real('FM_HydAttWSSlope', wsslope, ier)
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_HydAttWSType', wstype, ier)
    IF(wstype == 3) THEN
        HotStart = .TRUE.
    ELSE
        HotStart = .FALSE.
    ENDIF

    !CALL CG_IRIC_READ_INTEGER(fid, 'FM_HydAttHotStart', tmpint, ier)
    CALL cg_iRIC_Read_Real('FM_HydAttWS1DCD', ONEDCD, ier)

    CALL CG_IRIC_READ_STRING('FM_HydAttHSFile', CGNSHSFile, ier)
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_HydAttSolIndex', SolnIndex, ier)

    CALL CG_IRIC_READ_INTEGER(fid, 'FM_HydAttRoughnessType', roughnesstype, ier)
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_HydAttCDType', cdtype, ier)
    CALL cg_iRIC_Read_Real('FM_HydAttCD', constcd, ier)
    CALL cg_iRIC_Read_Real('FM_HydAttCDMin', cdmin, ier)
    CALL cg_iRIC_Read_Real('FM_HydAttCDMax', cdmax, ier)

    CALL cg_iRIC_Read_Real('FM_SolAttNSExtSlope', nsextslope, ier)
    !CALL CG_IRIC_READ_INTEGER_F('FM_SolAttNSExtShow', ShowGridExt, ier)
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_HydAttVelBCDS', vbcds, ier)
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_HydAttDryType', dryType, ier)
    CALL cg_iRIC_Read_Real('FM_HydAttDryMinDepth', hmin, ier)
    hmin = hmin*100.
    CALL cg_iRIC_Read_Real('FM_HydAttWetMinDepth', hwmin, ier)
    hwmin = hwmin*100.
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_HydAttWetIterInterval', hiterInterval, ier)
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_HydAttWetIterStop', hiterstop, ier)
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_HydAttCalcWetting', tmpint, ier)
    IF(tmpint == 1) THEN
        hcalcwetting = .TRUE.
    ELSE
        hcalcwetting = .FALSE.
    ENDIF

    CALL CG_IRIC_READ_INTEGER(fid, 'FM_HydAttLEVType', LEVType, ier)
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_HydAttLEVChangeIter', LEVChangeIter, ier)
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_LEVStartIter', LEVBegIter, ier)
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_LEVEndIter', LEVEndIter, ier)
    CALL cg_iRIC_Read_Real('FM_StartLEV', startLEV, ier)
    startLEV = startLEV*100*100
    CALL cg_iRIC_Read_Real('FM_EndLEV', endLEV, ier)
    endLEV = endLEV*100*100
    CALL cg_iRIC_Read_Real('FM_HydAttEVC', evc, ier)
    evc = evc*100.*100.

    CALL CG_IRIC_READ_INTEGER(fid, 'FM_Quasi3D', tmpint, ier)
    IF(tmpint == 1) THEN
        CALCQUASI3D = .TRUE.
    ELSE
        CALCQUASI3D = .FALSE.
    ENDIF

    CALL CG_IRIC_READ_INTEGER(fid, 'FM_Quasi3D_NZ', nz, ier)
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_Quasi3D_CalcRS', tmpint, ier)
    IF(tmpint == 1) THEN
        CALCQUASI3DRS = .TRUE.
    ELSE
        CALCQUASI3DRS = .FALSE.
    ENDIF

    CALL cg_iRIC_Read_Real('FM_Quasi3D_MinRS', MinRS, ier)
    MinRs = MinRs * 100.

    CALL CG_IRIC_READ_INTEGER(fid, 'FM_SedTrans', tmpint, ier)
    IF(tmpint == 1) THEN
        CALCCSED = .TRUE.
    ELSE
        CALCCSED = .FALSE.
    ENDIF
    CALL cg_iRIC_Read_Real('FM_SedDuneH', HD, ier)
    HD = HD * 100.
    CALL cg_iRIC_Read_Real('FM_SedDuneWL', WD, ier)
    WD = WD * 100.

    !CALL CG_IRIC_READ_INTEGER_F('FM_RsSmoothLvl', RSSMOO, ier)
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_SedSmoothLvl', SEDSMOO, ier)
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_SedSmoothWeight', SEDSMOOWGHT, ier)

    CALL CG_IRIC_READ_INTEGER(fid, 'FM_SedTrans_Type', TRANSEQTYPE, ier)
    CALL cg_iRIC_Read_Real('FM_SedGrainSize', din, ier)
    din = din*100.
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_SedBCNode', SEDBCNODE, ier)
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_SedBCFractionTaper', SEDEQNODE, ier)
    SEDEQNODE = SEDBCNODE+SEDEQNODE
    CALL cg_iRIC_Read_Real('FM_SedBCFraction', BCFRACTION, ier)
    !CALL cg_iRIC_Read_Real_F('FM_SedBCEquiMult', SEDEQMULT, ier)

    CALL CG_IRIC_READ_INTEGER(fid, 'FM_SedCalcGravCorr', CALCGRAVCORR, ier)
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_SedGravCorrType', GRAVCORRTYPE, ier)
    CALL cg_iRIC_Read_Real('FM_SedFlatBedCorrCoef', GRAVFLATBEDCORRCOEF, ier)
    CALL cg_iRIC_Read_Real('FM_SedSAngleOfRepose', SUBANGLEOFREPOSE, ier)
    CALL cg_iRIC_Read_Real('FM_SedTSFracDepth', TSFracDepth, ier)

    CALL CG_IRIC_READ_INTEGER(fid, 'FM_SedTransAuto', tmpint, ier)
    IF(tmpint == 1) THEN
        CALCSEDAUTO = .TRUE.
    ELSE
        CALCSEDAUTO = .FALSE.
    ENDIF

    !! variables for DT Wilcock-Kenworthy transport function !! rmcd 5/2/07

    CALL cg_iRIC_Read_Real('FM_Sed_WK_RG0', taustar_rg0, ier)
    CALL cg_iRIC_Read_Real('FM_Sed_WK_RG1', taustar_rg1, ier)
    CALL cg_iRIC_Read_Real('FM_Sed_WK_RS1', taustar_rs1, ier)
    CALL cg_iRIC_Read_Real('FM_Sed_WK_Alpha', alpha, ier)
    CALL cg_iRIC_Read_Real('FM_Sed_WK_AK', AK, ier)
    CALL cg_iRIC_Read_Real('FM_Sed_WK_Chi', Chi, ier)
    CALL cg_iRIC_Read_Real('FM_Sed_WK_PhiPrime', phi_prime, ier)
    CALL cg_iRIC_Read_Real('FM_Sed_WK_HMax', tmphmax, ier)
    tmphmax = tmphmax * 100.
    CALL cg_iRIC_Read_Real('FM_Sed_WK_ABH', tmpabh, ier)
    CALL cg_iRIC_Read_Real('FM_Sed_WK_BBH', tmpbbh, ier)
    CALL cg_iRIC_Read_Real('FM_Sed_WK_LSub', tmplsub, ier)
    tmplsub = tmplsub*100.
    CALL cg_iRIC_Read_Real('FM_Sed_WK_LSubActDepth', tmplsubactdepth, ier)
    tmplsubactdepth = tmplsubactdepth*100.
    CALL cg_iRIC_Read_Real('FM_Sed_WK_DSand', Dsand, ier)
    Dsand = Dsand *100.
    CALL cg_iRIC_Read_Real('FM_Sed_WK_DGravel', Dg, ier)
    Dg = Dg*100.
    CALL cg_iRIC_Read_Real('FM_Sed_WK_Z0Min', Z0Min, ier)
    Z0Min = Z0Min*100.
    CALL cg_iRIC_Read_Real('FM_Sed_WK_Z0Max', Z0Max, ier)
    Z0Max = Z0Max*100.

    CALL CG_IRIC_READ_INTEGER(fid, 'FM_Sed_WK_RoughType', WK_RoughType, ier)
    CALL cg_iRIC_Read_Real('FM_Sed_WK_RoughShapeParam', WK_RoughShapeParam, ier)
    CALL cg_iRIC_Read_Real('FM_Sed_WK_FsMin', Fsmin, ier)
    CALL CG_IRIC_READ_INTEGER(fid, 'FM_Sed_WK_ConstBndry', tmpint, ier)
    IF(tmpint == 1) THEN
        WKConstBndry = .TRUE.
    ELSE
        WKConstBndry = .FALSE.
    ENDIF
    !
    ! --- Parameters for Hot Start ---
    !
    CALL CG_IRIC_READ_INTEGER(fid, 'write_flag', i_re_flag_o, ier)
    !     CALL CG_IRIC_READ_INTEGER_F('read_flag', i_re_flag_i, ier)
    CALL CG_IRIC_READ_INTEGER(fid, 'n_tempfile', n_rest, ier)
    CALL CG_IRIC_READ_STRING('tmp_readfile', tmp_file_i, ier)
    CALL CG_IRIC_READ_STRING('tmp_pass', tmp_pass, ier)

    do ii=0,9
        CALL CG_IRIC_READ_REAL(tmp_caption(ii), opt_tmp(ii), ier)
    end do
    !
    do iii=1,n_rest
        do ii=0,n_rest-1
            if(opt_tmp(ii) /= opt_tmp(ii+1) &
                .or.opt_tmp(ii+1) < opt_tmp(ii)+dt) then
            if(opt_tmp(ii) > opt_tmp(ii+1)) then
                ttt = opt_tmp(ii)
                opt_tmp(ii) = opt_tmp(ii+1)
                opt_tmp(ii+1) = ttt
            end if
            else
                opt_tmp(ii+1) = opt_tmp(ii+1)+dt
            end if
        end do
    end do

    END SUBROUTINE




    END MODULE
