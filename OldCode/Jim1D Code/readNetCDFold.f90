MODULE readNetCDFMod
USE Support
IMPLICIT NONE
INCLUDE "netcdf.inc"
CONTAINS


	SUBROUTINE read_NetCDFInput(InputFile)

	CHARACTER(*), INTENT(IN) :: InputFile
	INTEGER :: status, NCID, tmpint
	REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar, tmpvar2
	INTEGER, ALLOCATABLE, DIMENSION(:) :: tmpvari
!	INTEGER, Dimension(1)::numCh
	INTEGER::numCh
	INTEGER::netID, topoDimID, topoVarID
	INTEGER::BCDimID, BCVarID, jctID, jctDimID, jctVarID
	INTEGER::numTopoPts, numXSects, maxelem
	INTEGER::chindex, bcindex, jctindex
	INTEGER::maxsectchan, totnumxsect, totnumgeopts 
	INTEGER::numjct, totnumjct, numchjct, maxchanjct, numbcts, numbc, totnumbc, tmp
	character(len = 14) strTopoDim
	character(len = 14) strTopo
	character(len = 12) strBCDim
	character(len = 6) strBC
	character(len = 18) strJctDim
	character(len = 12) strJct
	character(len = 3) indstr 
	character(len = 3) trmstr
	character(40) tmpDimstr
	character(40) tmpVarstr

	strTopoDim = '1D_NumTopoPts_'
	strTopo = '1D_Topography_'
	strBCDim = '1D_NumBCPts_'
	strBC = '1D_BC_'
	strJct = '1D_Junction_'
	strJctDim = '1D_NumJctChannels_'
	

	maxsectchan = -1e6
	maxchanjct = -1e6
	status = NF_OPEN(InputFile, NF_NOWRITE, NCID)
		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
!Get Channel params
	status = NF_INQ_VARID(NCID,'1DNetwork', netID)
		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
	status = NF_GET_VAR_INT(NCID, netID, numCH)
	DO chindex = 1, numCh
		WRITE(indstr,'(I3)') chindex-1
		indstr = ADJUSTL(indstr)
		tmpDimstr = strTopoDim//TRIM(indstr)
		tmpVarstr = strTopo//TRIM(indstr)
		
		status = NF_INQ_DIMID(NCID,TRIM(tmpDimstr), topoDimID)
		status = NF_INQ_DIMLEN(NCID, topoDimID, numTopoPts)
		totnumgeopts = totnumgeopts + (numTopoPts/3)
		status = NF_INQ_VARID(NCID,TRIM(tmpVarstr), topoVarID) 
		status = NF_GET_ATT_INT(NCID, topoVarID, 'numXSections', numXSects)
		totnumxsect = totnumxsect + numXSects
		ALLOCATE(tmpvari(numXSects), STAT = status)
		status = NF_GET_ATT_INT(NCID, topoVarID, 'numPtsXSections', tmpvari)
		maxsectchan = max(maxsectchan, maxval(tmpvari))
		DEALLOCATE(tmpvari)
		tmpDimstr = strBCDim//TRIM(indstr)
		tmpVarstr = strBC//TRIM(indstr)
		status = NF_INQ_DIMID(NCID, TRIM(tmpDimstr), BCDimID)
		status = NF_INQ_DIMLEN(NCID, BCDimID, numbc)
		totnumbc = totnumbc+numbc
		status = NF_INQ_VARID(NCID, TRIM(tmpVarstr), BCVarID)
		ALLOCATE(tmpvari(numbc), STAT = status)
		status = NF_GET_VAR_INT(NCID, BCVARID, tmpvari)
		DO bcindex = 1,numbc
			if(tmpvari(bcindex) < 10) then
				numbcts = numbcts+1
			endif
		ENDDO
		DEALLOCATE(tmpvari)

	ENDDO

	status = NF_INQ_VARID(NCID, '1DJunction', jctID)
	status = NF_GET_VAR_INT(NCID, jctID, numjct)
	DO jctindex = 1,numjct
		WRITE(indstr,'(I3)') jctindex-1
		indstr = ADJUSTL(indstr)
		tmpDimstr = strJctDim//TRIM(indstr)
		tmpVarstr = strJct//TRIM(indstr)
			
		status = NF_INQ_DIMID(NCID, TRIM(tmpDimstr), jctDimID)
		status = NF_INQ_DIMLEN(NCID, jctDimID, numchjct)
		status = NF_INQ_VARID(NCID, TRIM(tmpVarstr), jctVarID)
		ALLOCATE(tmpvari(numchjct), STAT = status)
		status = NF_GET_VAR_INT(NCID, jctVarID, tmpvari)
		maxchanjct = max(maxchanjct, numchjct)
		DEALLOCATE(tmpvari)
	ENDDO

	CALL AllocateSecs(numCh, maxsectchan, totnumxsect, totnumgeopts, numjct, maxchanjct, numbcts, totnumbc)
	tmp = 0




	END SUBROUTINE read_NetCDFInput

SUBROUTINE HANDLE_ERR(STATUS)
	INCLUDE "netcdf.inc"
	INTEGER STATUS
	IF (STATUS .NE. NF_NOERR) THEN
		PRINT *, NF_STRERROR(STATUS)
		STOP 'Stopped'
	ENDIF
END SUBROUTINE HANDLE_ERR


END MODULE readNetCDFMod




