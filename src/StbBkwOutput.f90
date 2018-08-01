    MODULE StpBckOutputMod

    USE Support
    IMPLICIT NONE
    !INCLUDE "netcdf.inc"
    CONTAINS

    SUBROUTINE writeStpBwOutput(wse)
    INTEGER, PARAMETER :: mp2 = KIND(1.0D0)
    REAL(kind = mp2), DIMENSION(:), INTENT(OUT) :: wse
    INTEGER :: i
    DO i = 1, nsc
        wse(i) = curxsc(i,1)
    ENDDO
    END SUBROUTINE


    !SUBROUTINE write_output(OutputFile)
    !	CHARACTER(*), INTENT(IN) :: OutputFile
    !	REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar1, tmpvar2
    !	INTEGER :: mdl_StpBkwWSElevID, NCID, status, i
    !		status = NF_OPEN(OutputFile, NF_WRITE, NCID)
    !		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
    !		status = NF_INQ_VARID(NCID, 'mdl_StpBkwWSElev', mdl_StpBkwWSElevID)
    !			IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
    !		ALLOCATE(tmpvar1(nsc), STAT = status)
    !		DO i = 1, nsc
    !			tmpvar1(i) = curxsc(i,1)
    !		ENDDO
    !		status = NF_PUT_VAR_REAL(NCID, mdl_StpBkwWSElevID, tmpvar1)
    !		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
    !
    !		DEALLOCATE(tmpvar1, STAT=status)
    !
    !		status = NF_CLOSE(NCID)
    !		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
    !
    !END SUBROUTINE write_output
    !SUBROUTINE HANDLE_ERR(STATUS)
    !	INCLUDE "netcdf.inc"
    !	INTEGER STATUS
    !	IF (STATUS .NE. NF_NOERR) THEN
    !		PRINT *, NF_STRERROR(STATUS)
    !		STOP 'Stopped'
    !	ENDIF
    !END SUBROUTINE HANDLE_ERR

    END