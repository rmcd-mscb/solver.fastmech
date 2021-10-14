    MODULE ReadCGNS
    USE RivVarMod
    USE CalcCond
    USE GridCoord
    USE GridCond
    USE RivCalcInitCond
    USE CSedMod_DT_SUSP
    USE CSedMod
    USE IRICMI
    IMPLICIT NONE
    INTEGER :: FID, BID, ZID

    CONTAINS

    SUBROUTINE read_iRIC_CGNS(InputFile, MIFLAG)
    IMPLICIT NONE
    !INCLUDE "cgnswin_f.h"
    INCLUDE "cgnslib_f.h"
    INCLUDE "iriclib_f.h"
    CHARACTER(*), INTENT(IN) ::InputFile
    INTEGER, INTENT(IN) :: MIFLAG
    CHARACTER(250) :: ZONENAME, BASENAME, USERNAME

    INTEGER :: IER, IRET
    write(0, *) 'Entered ric_iRIC_CGNS'
    write(0, *) 'MIFLAG', MIFLAG
    ier=0
    IF(MIFLAG.eq.0) then
        CALL cg_open_f(InputFile, CG_MODE_MODIFY, FID, IER)
        IF(IER .NE. 0) THEN
            call cg_error_print_f()
            !pause
            return
        ENDIF
!    ELSE
!        CALL iricmi_model_init(IER)
    ENDIF
    ! uncomment lines below for split solution
    !call iric_initoption_f(IRIC_OPTION_DIVIDESOLUTIONS, ier)
    !    if (ier /=0) STOP "*** Initialize option error***"
    IF(MIFLAG.eq.0) then
        CALL cg_iric_init_f(fid, ier)
        IF(IER .NE. 0) THEN
            call cg_error_print_f()
            !pause
            return
        ENDIF
        call iric_initoption_f(IRIC_OPTION_CANCEL, ier)
        if (ier /=0) STOP "*** Initialize option error***"
!    ELSE
!        CALL iricmi_model_init(ier)
    ENDIF
    !    CALL CGNS_Read_CalcCondtion_ForAlloc(FID, IER)
    write(0, *) 'Entering CGNS2_Read_CC_ForAlloc with MIFLAG', MIFLAG, IER
    CALL CGNS2_Read_CC_ForAlloc(MIFLAG,IER)
    write(0, *) 'Leaving CGNS2_Read_CC_ForAlloc with MIFLAG', MIFLAG, IER
    CALL CGNS2_READ_GRIDCOORD(MIFLAG,IER)
    write(0,*) 'leaving read grid coords'
    CALL CGNS2_Read_GridCondition(MIFLAG,IER)
    write(0,*) 'leaving read grid conds'
    CALL CGNS2_Read_CalcCondition(MIFLAG,IER)
    write(0,*) 'leaving read calc conds'

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
    write(0,*) 'Leaving initrating curves'
    CALL initTimeSeries()
    write(0,*) 'Leaving inittime series'

    Call initArrays()
    write(0,*) 'Leaving initarrays'
    CALL calcWSInitCond()
    write(0,*) 'Leaving calc swinitconds'

    !    CALL cg_close_f(FID, IER)
    END SUBROUTINE

    END MODULE