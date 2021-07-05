    MODULE GridCond
    USE RivVarMod
    IMPLICIT NONE


    CONTAINS
    SUBROUTINE CGNS2_Read_GridCondition(IER)
    IMPLICIT NONE

    INTEGER, INTENT(OUT) :: IER
    INTEGER :: status, i, j, count, countji, ierror
    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:) :: tmpreal4


    !    CALL CG_IRIC_GOTOGC_F(FID, NX, NY, IER)

    !    CALL ALLOCATE_GRIDCONDITION(NX+1, NY) !Need to increment NX by 1


    ALLOCATE(tmpreal4(ns2*nn), STAT=IER)

    CALL cg_iRIC_Read_Grid_Real_Node_F('Elevation', tmpreal4, IER)
    DO I = 1, ns2*nn
        IF( tmpreal4(i) < elevoffset) THEN
            elevoffset = tmpreal4(i)
        ENDIF
        IF(elevoffset < 0) THEN
            elevoffset = 0
        ENDIF
    ENDDO
    DO I= 1,ns2
        DO J=1,nn
            COUNT = ((I-1)*nn)+J
            countji = ((j-1)*ns2)+i
            eta2(I,J) = (tmpreal4(countji) - elevoffset)*100.
        ENDDO
    ENDDO

    CALL cg_iRIC_Read_Grid_Real_Node_F('roughness', tmpreal4, IER)
    DO I= 1,ns2
        DO J=1,nn
            COUNT = ((I-1)*nn)+J
            countji = ((j-1)*ns2)+i
            cd2(i,j) = tmpreal4(countji)
            znaught2(i,j) = tmpreal4(countji)*100.
        ENDDO
    ENDDO
    CALL cg_iRIC_Read_Grid_Real_Node_F('vegroughness', tmpreal4, IER)
    DO I= 1,ns2
        DO J=1,nn
            COUNT = ((I-1)*nn)+J
            countji = ((j-1)*ns2)+i
            cdv2(i,j) = tmpreal4(countji)
            !cdv2(i,j) = 0
        ENDDO
    ENDDO

    !CALL cg_iRIC_Read_Grid_Real_Node_F('minelevation', tmpreal4, IER)
    !DO I= 1,ns2
    !    DO J=1,nn
    !        COUNT = ((I-1)*nn)+J
    !        countji = ((j-1)*ns2)+i
    !        mineta2(I,J) = (tmpreal4(countji) - elevoffset)*100.
    !    ENDDO
    !ENDDO

    CALL cg_iRIC_Read_Grid_Real_Node_F('sandfraction', tmpreal4, IER)
    DO I= 1,ns2
        DO J=1,nn
            COUNT = ((I-1)*nn)+J
            countji = ((j-1)*ns2)+i
            Fracs(i,j) = tmpreal4(countji)
        ENDDO
    ENDDO

    CALL cg_iRIC_Read_Grid_Real_Node_F('sanddepth', tmpreal4, IER)
    DO I= 1,ns2
        DO J=1,nn
            COUNT = ((I-1)*nn)+J
            countji = ((j-1)*ns2)+i
            hfine(i,j) = tmpreal4(countji)*100. !Convert to cm
        ENDDO
    ENDDO

    !    CALL CG_IRIC_READ_GRIDREALNODE('FixedBedElevation', tmpreal4, IER)
    !    DO I= 1,NX
    !        DO J=1,NY
    !            COUNT = ((I-1)*NY)+J
    !            countji = ((j-1)*NX)+i
    !            ZbRock(I,J) = tmpreal4(countji)
    !        ENDDO
    !    ENDDO
    !
    !    CALL CG_IRIC_READ_GRIDREALNODE('VegetationDensity', tmpreal4, IER)
    !    DO I= 1,NX
    !        DO J=1,NY
    !            COUNT = ((I-1)*NY)+J
    !            countji = ((j-1)*NX)+i
    !            Vege(I,J) = tmpreal4(countji)
    !        ENDDO
    !    ENDDO
    !
    !    CALL CG_IRIC_READ_GRIDREALNODE('VegetationHeight', tmpreal4, IER)
    !    DO I= 1,NX
    !        DO J=1,NY
    !            COUNT = ((I-1)*NY)+J
    !            countji = ((j-1)*NX)+i
    !            Vegeh(I,J) = tmpreal4(countji)
    !        ENDDO
    !    ENDDO

    END SUBROUTINE CGNS2_Read_GridCondition



    END MODULE GridCond