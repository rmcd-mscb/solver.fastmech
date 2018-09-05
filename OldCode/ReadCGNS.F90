    MODULE ReadCGNS
    USE RivVarMod
    USE CalcCond
    USE GridCoord
    USE GridCond
    USE RivCalcInitCond
    USE CSedMod_DT_SUSP
    USE CSedMod
    IMPLICIT NONE
    INTEGER :: FID, BID, ZID

    CONTAINS

    SUBROUTINE read_iRIC_CGNS(InputFile)
    IMPLICIT NONE
    !INCLUDE "cgnswin_f.h"
    INCLUDE "cgnslib_f.h"
    INCLUDE "iriclib_f.h"
    CHARACTER(*), INTENT(IN) ::InputFile
    CHARACTER(250) :: ZONENAME, BASENAME, USERNAME

    INTEGER :: IER, IRET

    CALL cg_open_f(InputFile, CG_MODE_MODIFY, FID, IER)
    IF(IER .NE. 0) THEN
        call cg_error_print_f()
        !pause
        return
    ENDIF
    CGNSFILEID = FID
    ! uncomment lines below for split solution
    !call iric_initoption_f(IRIC_OPTION_DIVIDESOLUTIONS, ier)
    !    if (ier /=0) STOP "*** Initialize option error***"

    CALL cg_iric_init_f(fid, ier)
    IF(IER .NE. 0) THEN
        call cg_error_print_f()
        !pause
        return
    ENDIF
    call iric_initoption_f(IRIC_OPTION_CANCEL, ier)
    if (ier /=0) STOP "*** Initialize option error***"

    !    CALL CGNS_Read_CalcCondtion_ForAlloc(FID, IER)
    CALL CGNS2_Read_CC_ForAlloc(IER)
    CALL CGNS2_READ_GRIDCOORD(IER)
    CALL CGNS2_Read_GridCondition(IER)
    CALL CGNS2_Read_CalcCondition(IER)

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