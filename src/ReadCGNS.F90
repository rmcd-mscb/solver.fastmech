    MODULE ReadCGNS
    USE RivVarMod
    USE CalcCond
    USE GridCoord
    USE GridCond
    USE RivCalcInitCond
    USE CSedMod_DT_SUSP
    USE CSedMod
    USE iric
    IMPLICIT NONE
    INTEGER :: FID, BID, ZID

    CONTAINS

    SUBROUTINE read_iRIC_CGNS(InputFile)
    USE unlink_hdf5
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) ::InputFile
    CHARACTER(250) :: ZONENAME, BASENAME, USERNAME

    INTEGER :: IER, IRET

    CALL cg_iric_open(InputFile, IRIC_MODE_MODIFY, FID, IER)
    IF(IER .NE. 0) THEN
        !SRC call cg_error_print()
        !pause
        return
    ENDIF
    ! uncomment lines below for split solution
    !call iric_initoption_f(IRIC_OPTION_DIVIDESOLUTIONS, ier)
    !    if (ier /=0) STOP "*** Initialize option error***"

    call unlink_solutions(fid, ier)
    call iric_initoption(IRIC_OPTION_CANCEL, ier)
    if (ier /= 0) STOP "*** Initialize option error***"

    !    CALL CGNS_Read_CalcCondtion_ForAlloc(FID, IER)
    CALL CGNS2_Read_CC_ForAlloc(fid, IER)
    CALL CGNS2_READ_GRIDCOORD(fid, IER)
    CALL CGNS2_Read_GridCondition(fid, IER)
    CALL CGNS2_Read_CalcCondition(fid, IER)

    if(TRANSEQTYPE == 2) then
        call alloc_csed_DT()
    else
        call alloc_csed()
    endif

    IF(CALCQUASI3D) THEN
        CALL ALLOC_COMMON3D(ns2, nn, nz)
        IF(TRANSEQTYPE == 2) THEN
            CALL alloc_csed3d_dt() !array for Daniele Tonina Wilcock-Kenworthy
        ENDIF
    ENDIF


    mo=stot*(ns-1)/(ns2-1)
    wmax=w2(1)
    dn=wmax/(nn-1)
    ds=mo/(ns-1)
    nm=(nn+1)/2.

    CALL initRatingCurves()
    CALL initTimeSeries()

    Call initArrays()
    CALL calcWSInitCond()

    !    CALL cg_close_f(FID, IER)
    END SUBROUTINE

    END MODULE