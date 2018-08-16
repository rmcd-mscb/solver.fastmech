module fm_helpers
use RivVarWMod
use RivVarMod
use CalcCond
use RivCalcInitCond
use RivRoughnessMod
use RivVarTimeMod
use CSedMod_DT_SUSP
use CSedMod
contains
    SUBROUTINE dealloc_all()
    IMPLICIT NONE
    CALL dealloc_working()
    CALL dealloc_common2d()
    if(vbc == 1) then
        CALL DEALLOC_VELBC()
    endif

    CALL dealloc_init2d()
    CALL dealloc_roughness()
    CALL dealloc_TimeSeries()
    CALL dealloc_RatingCurves()
    CALL dealloc_TSNames()


    if(CALCQUASI3D) then
        CALL dealloc_common3d()
        if(TRANSEQTYPE == 1) THEN
            CAll dealloc_csed3d_dt()
        ENDIF
    endif
    !	    if(CALCCSED == 1) then
    !	        call dealloc_csed()
    !	        call dealloc_csed_DT()
    !	    endif
    if(TRANSEQTYPE == 2) then
        call dealloc_csed_DT()
    else
        call dealloc_csed()
    endif
    IF(CGNSFILEID > 0) THEN
        CALL cg_close_f(CGNSFILEID, errorcode)
        CGNSFILEID = 0
    ENDIF
    END SUBROUTINE dealloc_all
    
    SUBROUTINE Get_GRID_2D_COORD(ctype, tmp)
    IMPLICIT NONE
    CHARACTER(LEN=*) :: ctype
    REAL, DIMENSION(:), INTENT(OUT) :: tmp !changed this to just real but could create problems
    REAL(KIND=mp) :: xx, yy, ux, uy
    double precision :: rcos, rsin
    INTEGER :: i, j, count, ier
    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:) :: tmpval
    ALLOCATE(tmpval(ns2,nn), stat=ier)
    DO i=1,ns2
        DO j= 1,nn
            rcos = cos(phirotation(i))
            rsin = sin(phirotation(i))
            ux = x(i,j)*rcos - y(i,j)*rsin
            uy = x(i,j)*rsin + y(i,j)*rcos

            xx = x(i,j)*fcos - y(i,j)*fsin
            yy = x(i,j)*fsin + y(i,j)*fcos
            SELECT CASE (ctype)
                CASE('x')
                    tmpval(i,j) = (xx + xshift)/100.0D0
                CASE('y')
                    tmpval(i,j) = (yy + yshift)/100.0D0
                CASE DEFAULT
                    tmpval(i,j) = -9999
                    return
            END SELECT
        END DO
    END DO
    tmp = reshape(tmpval, (/ns2*nn/))
    deallocate(tmpval, stat=ier)
    END SUBROUTINE Get_GRID_2D_COORD

    SUBROUTINE Get_GRID_3D_COORD(ctype, tmpx)
    IMPLICIT NONE
    CHARACTER(LEN=*) :: ctype
    real, DIMENSION(:), INTENT(OUT) :: tmpx
    INTEGER :: i, j, k, ier
    Double Precision :: xx, yy, ux, uy
    double precision :: rcos, rsin
    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:,:) :: tmpval
    ALLOCATE(tmpval(ns2,nn,nz), stat=ier)
    DO k = 1,nz
        DO j=1,nn
            DO i= 1,ns2
                rcos = cos(phirotation(i))
                rsin = sin(phirotation(i))

                ux = x(i,j)*rcos - y(i,j)*rsin
                uy = x(i,j)*rsin + y(i,j)*rcos

                xx = x(i,j)*fcos - y(i,j)*fsin
                yy = x(i,j)*fsin + y(i,j)*fcos
                SELECT CASE (ctype)
                CASE('x')
                        tmpval(i,j,k) = (xx + xshift)/100.
                CASE('y')
                    tmpval(i,j,k) = (yy + yshift)/100.
                CASE('z')
                    tmpval(i,j,k) = (zz(i,j,k)/100.) + elevoffset
                END SELECT
                !ty(i,j,k) = (yy + yshift)/100.
                !tz(i,j,k) = (zz(i,j,k)/100.) + elevoffset
            END DO
        END DO
    END DO
    tmpx = reshape(tmpval, (/ns2*nn*nz/))
    DEALLOCATE(tmpval, stat=ier)
    END SUBROUTINE Get_GRID_3D_COORD
    
end module fm_helpers
