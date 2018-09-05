    MODULE CalcCond2
    !use fm_global
    USE rivvarmod2
    USE RivVarTimeMod2
    IMPLICIT NONE

    type calccond
        REAL(KIND=mp) :: erelax, urelax, arelax
        REAL(KIND=mp) :: fmdt, q
        REAL(KIND=mp) :: vdc, vac
        REAL(KIND=mp) :: wselev, wsupelev, wsslope
        REAL(KIND=mp), ALLOCATABLE, DIMENSION(:) :: vbcdist, vbcvel, vbcang

        LOGICAL :: HotStart


        INTEGER, pointer :: itm, interitm, iplinc
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
        INTEGER :: dryType, hiterInterval, hiterstop, imod
        LOGICAL :: iwetdry
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
    end type

    CONTAINS

    SUBROUTINE ALLOC_VELBC(object, size)
    IMPLICIT NONE
    type(calccond), intent(inout) :: object
    INTEGER, INTENT(IN) :: size
    INTEGER :: status
    ALLOCATE(object%vbcdist(size), STAT = status)
    ALLOCATE(object%vbcvel(size), STAT = status)
    ALLOCATE(object%vbcang(size), STAT = status)
    END SUBROUTINE

    SUBROUTINE DEALLOC_VELBC(object)
    IMPLICIT NONE
    type(calccond), intent(inout) :: object
    INTEGER :: status
    DEALLOCATE(object%vbcdist, STAT = status)
    DEALLOCATE(object%vbcvel, STAT = status)
    DEALLOCATE(object%vbcang, STAT = status)
    END SUBROUTINE


    SUBROUTINE CGNS2_Read_CC_ForAlloc(cc_object, rv_object, IER)
    IMPLICIT NONE
    type(calccond), intent(inout) :: cc_object
    type(rivvar), intent(inout) :: rv_object
    INTEGER, INTENT(OUT) :: IER
    INTEGER :: tmpint

    CALL CG_IRIC_READ_INTEGER_f('FM_SolAttNSExt', cc_object%nsext, ier)
    CALL CG_IRIC_READ_INTEGER_f('FM_SolAttItm', cc_object%itm, ier)
    CALL CG_IRIC_READ_INTEGER_F('FM_Quasi3D_NZ', tmpint, ier)
    rv_object%nz = tmpint
    CALL CG_IRIC_READ_INTEGER_F('FM_Quasi3D', tmpint, ier)
    IF(tmpint == 1) THEN
        cc_object%CALCQUASI3D = .TRUE.
    ELSE
        cc_object%CALCQUASI3D = .FALSE.
    ENDIF

    ENDSUBROUTINE

    SUBROUTINE CGNS2_Read_CalcCondition(cc_object, rv_object, rvt_object, IER)
    IMPLICIT NONE
    type(calccond), intent(inout) :: cc_object
    type(rivvar), intent(inout) :: rv_object
    type(rivvartime), intent(inout) :: rvt_object
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
    cc_object%i_tmp_count = 0
    do ii=0,9
        cc_object%tmp_dummy = "tmp0.d"
        write(cc_object%tmp_dummy(4:4),'(i1)') ii
        cc_object%tmp_file_o(ii) = cc_object%tmp_dummy
        cc_object%tmp_dummy = "outtime_tmp0"
        write(cc_object%tmp_dummy(12:12),'(i1)') ii
        cc_object%tmp_caption(ii) = cc_object%tmp_dummy
    end do

    !These are hard wired
    cc_object%IO_IBC= .TRUE.
    cc_object%IO_VELXY= .TRUE.
    cc_object%IO_WSE= .TRUE.
    cc_object%IO_DEPTH= .TRUE.
    cc_object%IO_ELEV= .TRUE.

    CALL CG_IRIC_READ_INTEGER_F('FM_SolWrtVelSN', tmpint, ier)
    if(tmpint == 0) then
        cc_object%IO_VELSN = .FALSE.
    else
        cc_object%IO_VELSN = .TRUE.
    endif
    CALL CG_IRIC_READ_INTEGER_F('FM_SolWrtShearXY', tmpint, ier)
    if(tmpint == 0) then
        cc_object%IO_SHEARXY = .FALSE.
    else
        cc_object%IO_SHEARXY= .TRUE.
    endif
    CALL CG_IRIC_READ_INTEGER_F('FM_SolWrtShearSN', tmpint, ier)
    if(tmpint == 0) then
        cc_object%IO_SHEARSN = .FALSE.
    else
        cc_object%IO_SHEARSN= .TRUE.
    endif
    CALL CG_IRIC_READ_INTEGER_F('FM_SolWrtShearDiv', tmpint, ier)
    if(tmpint == 0) then
        cc_object%IO_STRESSDIV = .FALSE.
    else
        cc_object%IO_STRESSDIV= .TRUE.
    endif
    CALL CG_IRIC_READ_INTEGER_F('FM_SolWrtCD', tmpint, ier)
    if(tmpint == 0) then
        cc_object%IO_CD = .FALSE.
    else
        cc_object%IO_CD= .TRUE.
    endif
    CALL CG_IRIC_READ_INTEGER_F('FM_SolWrtUnitDischarge', tmpint, ier)
    if(tmpint == 0) then
        cc_object%IO_UnitDisch = .FALSE.
    else
        cc_object%IO_UnitDisch = .TRUE.
    endif
    CALL CG_IRIC_READ_INTEGER_F('FM_SolWrtHArea', tmpint, ier)
    if(tmpint == 0) then
        cc_object%IO_HArea = .FALSE.
    else
        cc_object%IO_HArea = .TRUE.
    endif
    CALL CG_IRIC_READ_INTEGER_F('FM_SolWrtInitVel', tmpint, ier)
    if(tmpint == 0) then
        cc_object%IO_InitVel = .FALSE.
    else
        cc_object%IO_InitVel = .TRUE.
    endif
    CALL CG_IRIC_READ_INTEGER_F('FM_SolWrt3DOutput', tmpint, ier)
    if(tmpint == 0) then
        cc_object%IO_3DOUTPUT = .FALSE.
    else
        cc_object%IO_3DOUTPUT = .TRUE.
    endif
    CALL CG_IRIC_READ_INTEGER_F('FM_SolWrtHelix', tmpint, ier)
    if(tmpint == 0) then
        cc_object%IO_HELIX = .FALSE.
    else
        cc_object%IO_HELIX = .TRUE.
    endif




    cc_object%IO_VORT= .FALSE.
    cc_object%IO_KEG= .FALSE.
    cc_object%IO_VELSTRAIN= .FALSE.

    !CALL CG_IRIC_GOTOCC_F(FID, IER)

    CALL cg_iRIC_Read_Integer_f('FM_SolAttType', cc_object%solType, ier)

    CALL cg_iRIC_Read_Real_F('FM_SolAttERlx', cc_object%erelax, ier)
    CALL cg_iRIC_Read_Real_F('FM_SolAttURlx', cc_object%urelax, ier)
    CALL cg_iRIC_Read_Real_F('FM_SolAttARlx', cc_object%arelax, ier)

    CALL CG_IRIC_READ_INTEGER_F('FM_SolAttItm', cc_object%itm, ier)

    CALL cg_iRIC_Read_Real_F('FM_SolAttDT', cc_object%fmdt, ier)

    CALL CG_IRIC_READ_INTEGER_F('FM_SolAttInterItm', cc_object%interitm, ier)

    CALL CG_IRIC_READ_INTEGER_F('FM_SolAttPlinc', cc_object%iplinc, ier)

    !NEXT TWO DEFINED IN RivVarTimeMod
    CALL cg_iRIC_Read_Real_F('FM_SolAttVarDischSTime', rvt_object%vardischstarttime, ier)
    CALL cg_iRIC_Read_Real_F('FM_SolAttVarDischETime', rvt_object%vardischendtime, ier)

    CALL CG_IRIC_READ_INTEGER_F('FM_SolAttFlumeBndry', tmpint, ier)
    IF(tmpint == 1) THEN
        rv_object%FLUMEBNDRY = .TRUE.
    ELSE
        rv_object%FLUMEBNDRY = .FALSE.
    ENDIF

    CALL CG_IRIC_READ_INTEGER_F('FM_SolAttDbgStop', cc_object%debugstop, ier)
    CALL CG_IRIC_READ_INTEGER_F('FM_SolAttDbgTimeStep', cc_object%DbgTimeStep, ier)
    CALL CG_IRIC_READ_INTEGER_F('FM_SolAttDbgIterNum', cc_object%DbgIterNum, ier)

    CALL CG_IRIC_READ_INTEGER_F('FM_SolAttIterOut', tmpint, ier)
    IF(tmpint == 1) THEN
        cc_object%IterationOut = .TRUE.
    ELSE
        cc_object%IterationOut = .FALSE.
    ENDIF
    CALL CG_IRIC_READ_INTEGER_F('FM_SolAttIncPlOut', cc_object%IterPlOut, ier)

    CALL CG_IRIC_READ_INTEGER_F('FM_SolAttMaxInterIter', cc_object%maxInterIterMult, ier)

    CALL cg_iRIC_Read_Real_F('FM_HydAttQ', cc_object%q, ier)
    cc_object%q = cc_object%q*1.e6
    !NEXT TWO DEFINED IN RivVarTimeMod
    CALL CG_IRIC_READ_INTEGER_F('FM_HydAttVarDischType', rvt_object%VarDischType, ier)
    IF(rvt_object%VarDischType.eq.1)THEN
        !dischT and dischQ are allocated in _READ_FUNCTIONAL
        CALL CG_IRIC_READ_FUNCTIONALSIZE_F('FM_HydAttVarDischarge', sizeDischTQ, IER)
        ALLOCATE(dischT(sizeDischTQ), dischQ(sizeDischTQ), STAT=IER)
        CALL CG_IRIC_READ_FUNCTIONAL_F('FM_HydAttVarDischarge', dischT, dischQ, ier)
        tint(1) = sizeDischTQ
        CALL SetTimeSeriesNumPts(rvt_object, 1, tint)
        ALLOCATE(treal(2*sizeDischTQ), STAT = ier)
        count2=0
        DO ii = 1, sizeDischTQ
            count2 = (ii-1)*2
            treal(count2+1) = dischT(ii)
            treal(count2+2) = dischQ(ii)
        ENDDO
        Call SetTimeSeries(rvt_object, (2*sizeDischTQ), treal)
        tint(1) = 0
        CALL SetTimeSeriesType(rvt_object, 1,tint)
        DEALLOCATE(treal, STAT = ier)
        DEALLOCATE(dischT, STAT = ier)
        DEALLOCATE(dischQ, STAT = ier)
    ENDIF

    CALL CG_IRIC_READ_INTEGER_F('FM_HydAttVarStgType', rvt_object%VarStageType, ier)
    IF(rvt_object%VarStageType.eq.1)THEN
        CALL CG_IRIC_READ_FUNCTIONALSIZE_F('FM_HydAttVarStageTS', sizeStageTQ, IER)
        ALLOCATE(stageT(sizeStageTQ), stageH(sizeStageTQ), STAT=IER)
        CALL CG_IRIC_READ_FUNCTIONAL_F('FM_HydAttVarStageTS', stageT, stageH, ier)
        tint(1) = sizeStageTQ
        CALL SetRatingCurveNumPts(rvt_object, 1, tint)
        ALLOCATE(treal(2*sizeStageTQ), STAT = ier)
        count2=0
        DO ii = 1, sizeStageTQ
            count2 = (ii-1)*2
            treal(count2+1) = stageT(ii)
            treal(count2+2) = (stageH(ii)-rv_object%elevOffset)*100.
        ENDDO
        Call SetRatingCurves(rvt_object, (2*sizeStageTQ), treal)
        tint(1) = 0
        CALL SetRatingCurveType(rvt_object, 1,tint)
        DEALLOCATE(treal, STAT = ier)
        DEALLOCATE(stageT, STAT = ier)
        DEALLOCATE(stageH, STAT = ier)

    ELSE IF(rvt_object%VarStageType.eq.2) THEN
        CALL CG_IRIC_READ_FUNCTIONALSIZE_F('FM_HydAttVarStageRC', sizeStageHQ, IER)
        ALLOCATE(stageQ(sizeStageHQ), stageH(sizeStageHQ), STAT=IER)
        CALL CG_IRIC_READ_FUNCTIONAL_F('FM_HydAttVarStageRC', stageQ, stageH, ier)
        tint(1) = sizeStageHQ
        CALL SetRatingCurveNumPts(rvt_object, 1, tint)
        ALLOCATE(treal(2*sizeStageHQ), STAT = ier)
        count2=0
        DO ii = 1, sizeStageHQ
            count2 = (ii-1)*2
            treal(count2+1) = stageQ(ii)
            treal(count2+2) = (stageH(ii) - rv_object%elevoffset)*100
        ENDDO
        Call SetRatingCurves(rvt_object, (2*sizeStageHQ), treal)
        tint(1) = 1
        CALL SetRatingCurveType(rvt_object, 1,tint)
        DEALLOCATE(treal, STAT = ier)
        DEALLOCATE(stageH, STAT = ier)
        DEALLOCATE(stageQ, STAT = ier)

    ENDIF

    CALL cg_iRIC_Read_Real_F('FM_HydAttVelDepthCoef', cc_object%vdc, ier)
    CALL cg_iRIC_Read_Real_F('FM_HydAttVelAngleCoef', cc_object%vac, ier)
    CALL CG_IRIC_READ_INTEGER_F('FM_HydAttVelBC', cc_object%vbc, ier)
    IF(cc_object%vbc == 1) THEN
        CALL CG_IRIC_READ_FUNCTIONALSIZE_F('FM_HydAttVarVelBC', tmpint, IER)
        Allocate(xtmp(tmpint), ytmp(tmpint), ztmp(tmpint), STAT=IER)
        CALL ALLOC_VELBC(cc_object, tmpint)
        !   CALL CG_IRIC_READ_FUNCTIONAL_F('FM_HydAttVarVelBC', xtmp, ytmp, ztmp, ier)
        !   CALL CG_IRIC_READ_FUNCTIONAL_F('FM_HydAttVarVelBC', xtmp, ytmp, ier)
        CALL CG_IRIC_READ_FUNCTIONALWITHNAME_F("FM_HydAttVarVelBC","Normalized_Distance",xtmp,ier)
        CALL CG_IRIC_READ_FUNCTIONALWITHNAME_F("FM_HydAttVarVelBC","Normalized_Velocity",ytmp,ier)
        CALL CG_IRIC_READ_FUNCTIONALWITHNAME_F("FM_HydAttVarVelBC","Angle",ztmp,ier)
        IF(ier.eq.0) THEN
            DO i = 1, tmpint
                cc_object%vbcdist(i) = xtmp(i)
                cc_object%vbcvel(i) = ytmp(i)
                cc_object%vbcang(i) = ztmp(i)
            ENDDO
        ENDIF
        DEALLOCATE(xtmp, STAT=ier)
        DEALLOCATE(ytmp, STAT=ier)
        DEALLOCATE(ztmp, STAT=ier)
    ENDIF

    CALL cg_iRIC_Read_Real_F('FM_HydAttWS', cc_object%wselev, ier)
    cc_object%wselev=(cc_object%wselev - rv_object%elevoffset)*100.
    CALL cg_iRIC_Read_Real_F('FM_HydAttWS2', cc_object%wsupelev, ier)
    cc_object%wsupelev=(cc_object%wsupelev - rv_object%elevoffset)*100.
    CALL cg_iRIC_Read_Real_F('FM_HydAttWSSlope', cc_object%wsslope, ier)
    CALL CG_IRIC_READ_INTEGER_F('FM_HydAttWSType', cc_object%wstype, ier)
    IF(cc_object%wstype == 3) THEN
        cc_object%HotStart = .TRUE.
    ELSE
        cc_object%HotStart = .FALSE.
    ENDIF

    !CALL CG_IRIC_READ_INTEGER_F('FM_HydAttHotStart', tmpint, ier)
    CALL cg_iRIC_Read_Real_F('FM_HydAttWS1DCD', cc_object%ONEDCD, ier)

    CALL CG_IRIC_READ_STRING_F('FM_HydAttHSFile', cc_object%CGNSHSFile, ier)
    CALL CG_IRIC_READ_INTEGER_F('FM_HydAttSolIndex', cc_object%SolnIndex, ier)

    CALL CG_IRIC_READ_INTEGER_F('FM_HydAttRoughnessType', cc_object%roughnesstype, ier)
    CALL CG_IRIC_READ_INTEGER_F('FM_HydAttCDType', cc_object%cdtype, ier)
    CALL cg_iRIC_Read_Real_F('FM_HydAttCD', cc_object%constcd, ier)
    CALL cg_iRIC_Read_Real_F('FM_HydAttCDMin', cc_object%cdmin, ier)
    CALL cg_iRIC_Read_Real_F('FM_HydAttCDMax', cc_object%cdmax, ier)

    CALL cg_iRIC_Read_Real_F('FM_SolAttNSExtSlope', cc_object%nsextslope, ier)
    !CALL CG_IRIC_READ_INTEGER_F('FM_SolAttNSExtShow', ShowGridExt, ier)
    CALL CG_IRIC_READ_INTEGER_F('FM_HydAttVelBCDS', cc_object%vbcds, ier)
    CALL CG_IRIC_READ_INTEGER_F('FM_HydAttDryType', cc_object%dryType, ier)
    CALL cg_iRIC_Read_Real_F('FM_HydAttDryMinDepth', cc_object%hmin, ier)
    cc_object%hmin = cc_object%hmin*100.
    CALL cg_iRIC_Read_Real_F('FM_HydAttWetMinDepth', cc_object%hwmin, ier)
    cc_object%hwmin = cc_object%hwmin*100.
    CALL CG_IRIC_READ_INTEGER_F('FM_HydAttWetIterInterval', cc_object%hiterInterval, ier)
    cc_object%imod = cc_object%hiterinterval
    CALL CG_IRIC_READ_INTEGER_F('FM_HydAttWetIterStop', cc_object%hiterstop, ier)
    CALL CG_IRIC_READ_INTEGER_F('FM_HydAttCalcWetting', tmpint, ier)
    IF(tmpint == 1) THEN
        cc_object%iwetdry = .TRUE.
    ELSE
        cc_object%iwetdry = .FALSE.
    ENDIF

    CALL CG_IRIC_READ_INTEGER_F('FM_HydAttLEVType', cc_object%LEVType, ier)
    CALL CG_IRIC_READ_INTEGER_F('FM_HydAttLEVChangeIter', cc_object%LEVChangeIter, ier)
    CALL CG_IRIC_READ_INTEGER_F('FM_LEVStartIter', cc_object%LEVBegIter, ier)
    CALL CG_IRIC_READ_INTEGER_F('FM_LEVEndIter', cc_object%LEVEndIter, ier)
    CALL cg_iRIC_Read_Real_F('FM_StartLEV', cc_object%startLEV, ier)
    cc_object%startLEV = cc_object%startLEV*100*100
    CALL cg_iRIC_Read_Real_F('FM_EndLEV', cc_object%endLEV, ier)
    cc_object%endLEV = cc_object%endLEV*100*100
    CALL cg_iRIC_Read_Real_F('FM_HydAttEVC', cc_object%evc, ier)
    cc_object%evc = cc_object%evc*100.*100.

    CALL CG_IRIC_READ_INTEGER_F('FM_Quasi3D', tmpint, ier)
    IF(tmpint == 1) THEN
        cc_object%CALCQUASI3D = .TRUE.
    ELSE
        cc_object%CALCQUASI3D = .FALSE.
    ENDIF

    CALL CG_IRIC_READ_INTEGER_F('FM_Quasi3D_NZ', rv_object%nz, ier)
    CALL CG_IRIC_READ_INTEGER_F('FM_Quasi3D_CalcRS', tmpint, ier)
    IF(tmpint == 1) THEN
        cc_object%CALCQUASI3DRS = .TRUE.
    ELSE
        cc_object%CALCQUASI3DRS = .FALSE.
    ENDIF

    CALL cg_iRIC_Read_Real_F('FM_Quasi3D_MinRS', cc_object%MinRS, ier)
    cc_object%MinRs = cc_object%MinRs * 100.

    CALL CG_IRIC_READ_INTEGER_F('FM_SedTrans', tmpint, ier)
    IF(tmpint == 1) THEN
        cc_object%CALCCSED = .TRUE.
    ELSE
        cc_object%CALCCSED = .FALSE.
    ENDIF
    CALL cg_iRIC_Read_Real_F('FM_SedDuneH', cc_object%HD, ier)
    cc_object%HD = cc_object%HD * 100.
    CALL cg_iRIC_Read_Real_F('FM_SedDuneWL', cc_object%WD, ier)
    cc_object%WD = cc_object%WD * 100.

    CALL CG_IRIC_READ_INTEGER_F('FM_RsSmoothLvl', cc_object%RSSMOO, ier)
    CALL CG_IRIC_READ_INTEGER_F('FM_SedSmoothLvl', cc_object%SEDSMOO, ier)
    CALL CG_IRIC_READ_INTEGER_F('FM_SedSmoothWeight', cc_object%SEDSMOOWGHT, ier)

    CALL CG_IRIC_READ_INTEGER_F('FM_SedTrans_Type', cc_object%TRANSEQTYPE, ier)
    CALL cg_iRIC_Read_Real_F('FM_SedGrainSize', cc_object%din, ier)
    cc_object%din = cc_object%din*100.
    CALL CG_IRIC_READ_INTEGER_F('FM_SedBCNode', cc_object%SEDBCNODE, ier)
    CALL CG_IRIC_READ_INTEGER_F('FM_SedBCFractionTaper', cc_object%SEDEQNODE, ier)
    cc_object%SEDEQNODE = cc_object%SEDBCNODE+cc_object%SEDEQNODE
    CALL cg_iRIC_Read_Real_F('FM_SedBCFraction', cc_object%BCFRACTION, ier)
    !CALL cg_iRIC_Read_Real_F('FM_SedBCEquiMult', SEDEQMULT, ier)

    CALL CG_IRIC_READ_INTEGER_F('FM_SedCalcGravCorr', cc_object%CALCGRAVCORR, ier)
    CALL CG_IRIC_READ_INTEGER_F('FM_SedGravCorrType', cc_object%GRAVCORRTYPE, ier)
    CALL cg_iRIC_Read_Real_F('FM_SedFlatBedCorrCoef', cc_object%GRAVFLATBEDCORRCOEF, ier)
    CALL cg_iRIC_Read_Real_F('FM_SedSAngleOfRepose', cc_object%SUBANGLEOFREPOSE, ier)
    CALL cg_iRIC_Read_Real_F('FM_SedTSFracDepth', cc_object%TSFracDepth, ier)

    CALL CG_IRIC_READ_INTEGER_F('FM_SedTransAuto', tmpint, ier)
    IF(tmpint == 1) THEN
        cc_object%CALCSEDAUTO = .TRUE.
    ELSE
        cc_object%CALCSEDAUTO = .FALSE.
    ENDIF

    !! variables for DT Wilcock-Kenworthy transport function !! rmcd 5/2/07

    CALL cg_iRIC_Read_Real_F('FM_Sed_WK_RG0', cc_object%taustar_rg0, ier)
    CALL cg_iRIC_Read_Real_F('FM_Sed_WK_RG1', cc_object%taustar_rg1, ier)
    CALL cg_iRIC_Read_Real_F('FM_Sed_WK_RS1', cc_object%taustar_rs1, ier)
    CALL cg_iRIC_Read_Real_F('FM_Sed_WK_Alpha', cc_object%alpha, ier)
    CALL cg_iRIC_Read_Real_F('FM_Sed_WK_AK', cc_object%AK, ier)
    CALL cg_iRIC_Read_Real_F('FM_Sed_WK_Chi', cc_object%Chi, ier)
    CALL cg_iRIC_Read_Real_F('FM_Sed_WK_PhiPrime', cc_object%phi_prime, ier)
    CALL cg_iRIC_Read_Real_F('FM_Sed_WK_HMax', cc_object%tmphmax, ier)
    cc_object%tmphmax = cc_object%tmphmax * 100.
    CALL cg_iRIC_Read_Real_F('FM_Sed_WK_ABH', cc_object%tmpabh, ier)
    CALL cg_iRIC_Read_Real_F('FM_Sed_WK_BBH', cc_object%tmpbbh, ier)
    CALL cg_iRIC_Read_Real_F('FM_Sed_WK_LSub', cc_object%tmplsub, ier)
    cc_object%tmplsub = cc_object%tmplsub*100.
    CALL cg_iRIC_Read_Real_F('FM_Sed_WK_LSubActDepth', cc_object%tmplsubactdepth, ier)
    cc_object%tmplsubactdepth = cc_object%tmplsubactdepth*100.
    CALL cg_iRIC_Read_Real_F('FM_Sed_WK_DSand', cc_object%Dsand, ier)
    cc_object%Dsand = cc_object%Dsand *100.
    CALL cg_iRIC_Read_Real_F('FM_Sed_WK_DGravel', cc_object%Dg, ier)
    cc_object%Dg = cc_object%Dg*100.
    CALL cg_iRIC_Read_Real_F('FM_Sed_WK_Z0Min', cc_object%Z0Min, ier)
    cc_object%Z0Min = cc_object%Z0Min*100.
    CALL cg_iRIC_Read_Real_F('FM_Sed_WK_Z0Max', cc_object%Z0Max, ier)
    cc_object%Z0Max = cc_object%Z0Max*100.

    CALL CG_IRIC_READ_INTEGER_F('FM_Sed_WK_RoughType', cc_object%WK_RoughType, ier)
    CALL cg_iRIC_Read_Real_F('FM_Sed_WK_RoughShapeParam', cc_object%WK_RoughShapeParam, ier)
    CALL cg_iRIC_Read_Real_F('FM_Sed_WK_FsMin', cc_object%Fsmin, ier)
    CALL CG_IRIC_READ_INTEGER_F('FM_Sed_WK_ConstBndry', tmpint, ier)
    IF(tmpint == 1) THEN
        cc_object%WKConstBndry = .TRUE.
    ELSE
        cc_object%WKConstBndry = .FALSE.
    ENDIF
    !
    ! --- Parameters for Hot Start ---
    !
    CALL CG_IRIC_READ_INTEGER_F('write_flag', cc_object%i_re_flag_o, ier)
    !     CALL CG_IRIC_READ_INTEGER_F('read_flag', i_re_flag_i, ier)
    CALL CG_IRIC_READ_INTEGER_F('n_tempfile', cc_object%n_rest, ier)
    CALL CG_IRIC_READ_STRING_F('tmp_readfile', cc_object%tmp_file_i, ier)
    CALL CG_IRIC_READ_STRING_F('tmp_pass', cc_object%tmp_pass, ier)

    do ii=0,9
        CALL CG_IRIC_READ_REAL_F(cc_object%tmp_caption(ii), cc_object%opt_tmp(ii), ier)
    end do
    !
    do iii=1,cc_object%n_rest
        do ii=0,cc_object%n_rest-1
            if(cc_object%opt_tmp(ii) /= cc_object%opt_tmp(ii+1) &
                .or.cc_object%opt_tmp(ii+1) < cc_object%opt_tmp(ii)+cc_object%fmdt) then
            if(cc_object%opt_tmp(ii) > cc_object%opt_tmp(ii+1)) then
                ttt = cc_object%opt_tmp(ii)
                cc_object%opt_tmp(ii) = cc_object%opt_tmp(ii+1)
                cc_object%opt_tmp(ii+1) = ttt
            end if
            else
                cc_object%opt_tmp(ii+1) = cc_object%opt_tmp(ii+1)+cc_object%fmdt
            end if
        end do
    end do

    END SUBROUTINE




    END MODULE
