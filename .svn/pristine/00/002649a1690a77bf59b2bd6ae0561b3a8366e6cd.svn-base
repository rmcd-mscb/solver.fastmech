MODULE writeNetCDFMod
USE Support
IMPLICIT NONE
INCLUDE "netcdf.inc"
CONTAINS
	SUBROUTINE write_NetCDFSBWOutput(OutputFile)
		CHARACTER(*), INTENT(IN) :: OutputFile
		INTEGER :: status, NCID, WSEID, i, QID
		INTEGER :: start1, end1
		status = NF_OPEN(OutputFile, NF_WRITE, NCID)
			IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
		status = NF_INQ_VARID(NCID,'1DSBW_Elevation', WSEID)
			IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
		end1 = 1
		DO i = 1, nsc
			start1 = i
			status = NF_PUT_VARA_REAL(NCID, WSEID, start1, end1, curxsc(i,1))
		ENDDO
		status = NF_CLOSE(NCID)

	END SUBROUTINE

	SUBROUTINE write_NetCDFOutput(OutputFile)

	CHARACTER(*), INTENT(IN) :: OutputFile
	INTEGER :: status, NCID, WSEID, i, QID
	INTEGER :: RTWSXID, RTWSYID, LTWSXID, LTWSYID
	INTEGER :: TimeID
	INTEGER, DIMENSION(2) :: istart, iend
	INTEGER :: start1, end1
	REAL:: minv, maxv
	tstep = tstep+1
	end1 = 1
	iend(1) = 1
	iend(2) = 1
	status = NF_OPEN(OutputFile, NF_WRITE, NCID)
		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
	status = NF_INQ_VARID(NCID,'1DWater_Surface_Elevation', WSEID)
		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
	status = NF_INQ_VARID(NCID,'1DDischarge', QID)
		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
	status = NF_INQ_VARID(NCID,'1DTime', TimeID)
		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
	status = NF_INQ_VARID(NCID,'1DRTWSX', RTWSXID)
		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
	status = NF_INQ_VARID(NCID,'1DRTWSY', RTWSYID)
		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
	status = NF_INQ_VARID(NCID,'1DLTWSX', LTWSXID)
		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
	status = NF_INQ_VARID(NCID,'1DLTWSY', LTWSYID)
		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
	DO i = 1, nsc
		istart(2) = tstep
		istart(1) = i
		status = NF_PUT_VARA_REAL(NCID, WSEID, istart, iend, curxsc(i, 1))
			IF (status .NE. NF_NOERR)	CALL HANDLE_ERR(status)
		status = NF_PUT_VARA_REAL(NCID, QID, istart, iend, curxsc(i, 7))
			IF (status .NE. NF_NOERR)	CALL HANDLE_ERR(status)
		status = NF_PUT_VARA_REAL(NCID, RTWSXID, istart, iend, plnxsc(i,3))
			IF (status .NE. NF_NOERR)	CALL HANDLE_ERR(status)
		status = NF_PUT_VARA_REAL(NCID, RTWSYID, istart, iend, plnxsc(i,4))
			IF (status .NE. NF_NOERR)	CALL HANDLE_ERR(status)
		status = NF_PUT_VARA_REAL(NCID, LTWSXID, istart, iend, plnxsc(i,5))
			IF (status .NE. NF_NOERR)	CALL HANDLE_ERR(status)
		status = NF_PUT_VARA_REAL(NCID, LTWSYID, istart, iend, plnxsc(i,6))
			IF (status .NE. NF_NOERR)	CALL HANDLE_ERR(status)
	ENDDO
	start1 = tstep
	status = NF_PUT_VARA_REAL(NCID, TimeID, start1, end1, ttime)
		IF (status .NE. NF_NOERR)	CALL HANDLE_ERR(status)
		
	status = NF_CLOSE(NCID)
!	status = NF_PUT_VAR_FLOAT(

	END SUBROUTINE write_NetCDFOutput

	SUBROUTINE HANDLE_ERR(STATUS)
		INCLUDE "netcdf.inc"
		INTEGER STATUS
		IF (STATUS .NE. NF_NOERR) THEN
!			write(2,*) NF_STRERROR(STATUS)
			STOP 'Stopped'
		ENDIF
	END SUBROUTINE HANDLE_ERR


END MODULE writeNetCDFMod




