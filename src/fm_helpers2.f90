module fm_helpers2
use RivVarWMod2
use RivVarMod2
use CalcCond2
!use RivCalcInitCond2
use RivVarTimeMod2
!use CSedMod_DT_SUSP
use CSedMod2
contains
    SUBROUTINE dealloc_all(rvo,rwvo,cco,rvto,csedo)
    IMPLICIT NONE
    type(rivvar), intent(inout) :: rvo
    type(riv_w_var), intent(inout) :: rwvo
    type(calccond), intent(inout) :: cco
    type(rivvartime), intent(inout) :: rvto
    type(csed), intent(inout) :: csedo
    CALL dealloc_working(rwvo)
    CALL dealloc_common2d(rvo)
    if(cco%vbc == 1) then
        CALL DEALLOC_VELBC(cco)
    endif
    call dealloc_hotstart(cco)
    CALL dealloc_init2d(rvo)
    !CALL dealloc_roughness()
    CALL dealloc_TimeSeries(rvto)
    CALL dealloc_RatingCurves(rvto)
    CALL dealloc_TSNames(rvto)


    if(cco%CALCQUASI3D) then
        CALL dealloc_common3d(rvo)
        if(cco%TRANSEQTYPE == 2) THEN
            !CAll dealloc_csed3d_dt()
        ENDIF
    endif
    if(cco%CALCCSED == 1) then

        if(cco%TRANSEQTYPE == 2) then
            !call dealloc_csed_DT()
        else
            call dealloc_csed(csedo)
        endif
    endif
    IF(rvo%CGNSFILEID > 0) THEN
        CALL cg_close_f(rvo%CGNSFILEID, rvo%errorcode)
        rvo%CGNSFILEID = 0
    ENDIF
    END SUBROUTINE dealloc_all
    
    
end module fm_helpers2
