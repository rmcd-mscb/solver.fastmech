MODULE StpBckInputMod
USE Support
IMPLICIT NONE
!INCLUDE "netcdf.inc"
CONTAINS

SUBROUTINE setStpBwInput(nsec, nsecpts, tmptopo, tmpx, tmpy, tmpcd1, q, wselev)
    INTEGER, INTENT(IN) :: nsec, nsecpts
    REAL, DIMENSION(:), INTENT(IN)  :: tmptopo, tmpx, tmpy
    REAL, INTENT(IN) :: tmpcd1, q, wselev
    
    INTEGER :: i,j, jp, k
    INTEGER :: npts, ns, nss, nps, nsp, nsf, counti, countij 
    INTEGER :: ncd, nnpt, npf, nonds, maxbw, nods, ntbc
    INTEGER :: nchx
    
    REAL :: emin, eav
    npts = nsec*nsecpts

	CALL AllocateSecs(1, nsec, nsec, npts, 0, 0, 2)
	CALL AllocateTBC(2)
	ellim = -1.e9
	nchx = 1
	ns = nsec
	chbc(1,1) = 2
	chbc(1,2) = 1
	nss = 1
	nps = 1
	nsf=nss+ns-1
	nsce(1)=nsf

	do i = 1,nsec
		counti = ((i-1)*nsecpts)
		xsloc(i,1) = tmpx(counti+1)
		xsloc(i,2) = tmpy(counti+1)
		xsloc(i,3) = tmpx(counti+nsecpts)
		xsloc(i,4) = tmpy(counti+nsecpts)
		xsloc(i,5) = 0.
		xsloc(i,6) = sqrt((xsloc(i,1)-xsloc(i,3))**2+(xsloc(i,2)-xsloc(i,4))**2) ! length

			npt = nsecpts
			ncd = 1
			nnpt = nsecpts
		npf=nps+npt-1
		nscpte(i)=npf
		emin=1.e9
		if(i==1)eav = 0
		do j = 1,npt
			countij = ((i-1)*nsecpts) + j
			jp = nps+j-1

			if(j == 1) then
				xspts(jp,1) = 0.
			else
				xspts(jp,1) = sqrt((xsloc(i,1)-tmpx(countij))**2 + &
							(xsloc(i,2)-tmpy(countij))**2)
			endif
			xspts(jp,2) = tmptopo(countij)
			emin = min(emin, xspts(jp,2))
			if(i==1)eav = eav+(tmptopo(countij)-wselev)
			
!			if(tmpcdtype.ne.0) then
				xspts(jp,3) = tmpcd1
!			else
!				xspts(jp,3) = tmpcd(countij)
!			endif
		
		enddo
		if(i==1)eav = eav/npt
		do j = 1,npt
			countij = ((i-1)*nsecpts) + j
			jp = nps+j-1
			
!			if(tmpcdtype.ne.0) then
!				xspts(jp,3) = sqrt(tmpcd1*(eav)**.333/9.81);
				xspts(jp,3) = tmpcd1;
!			else
!!				xspts(jp,3) = sqrt(tmpcd(countij)*(eav)**.333/9.81);
!				xspts(jp,3) = tmpcd(countij);
!			endif
		
		enddo

		xsloc(i,5) = emin
		if(emin>ellim(1)) then
			ellim(1) = emin
			limsec(1) = i
		endif
		nps = npf+1
		nss = nsf+1

	enddo
	
	nonds = 0
	jctsecs = 0
	maxbw = 1
	CALL SetBW(2*(maxbw+1))
	nods = nsce(nchx)
	DO i = 1, nods
		CALL Farea(i, xsloc(i,5)+etol)
	ENDDO
	simstrt = 0
	simend = 0

	ntbc = 0
	tsbcloc = 0

	DO i = 1, nchx
		DO j = 1,2
		if(chbc(i,j)==1.or.chbc(i,j)==2) then
			tsbcloc(i,2*j-1) = ntbc+1
			DO k = 1,1
				jp = ntbc + k
				tmsrbc(jp, 1) = 0
				if(j==1) then
					tmsrbc(jp, 2) = q
				else
					tmsrbc(jp, 2) = wselev
				endif
			ENDDO
			ntbc = ntbc+1
			tsbcloc(i,2*k) = ntbc
		endif
		ENDDO
	ENDDO


END SUBROUTINE

!SUBROUTINE read_input(InputFile)
!	CHARACTER(*), INTENT(IN) :: InputFile
!	INTEGER :: status, NCID, tmpint
!	REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar, tmpvar2
!	INTEGER, ALLOCATABLE, DIMENSION(:) :: tmpvari
!
!	INTEGER :: i, j, counti, countij, k, jp
!	INTEGER :: nsecID, nsecptsID
!	INTEGER :: nsec, nsecpts, npts
!	INTEGER :: nchx, ns, nss, nps, nsf, npf
!	REAL, ALLOCATABLE, DIMENSION(:) :: tmptopo, tmpx, tmpy, tmpcd
!	INTEGER :: mdl_GridXID, mdl_GridYID, mdl_GridTopoID, mdl_GridCDID
!	INTEGER :: tmpcdtype, ncd, nnpt, maxbw, nonds, jctsecs, nods, ntbc
!	REAL :: tmpcd1, emin, eav
!	REAL :: q, wselev
!
!
!
!	status = NF_OPEN(InputFile, NF_NOWRITE, NCID)
!		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
!!Get grid params
!	status = NF_INQ_DIMID(NCID, 'mdl_GridNSPts', nsecID)
!		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
!	status = NF_INQ_DIMLEN(NCID, nsecID, nsec)
!		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
!
!	status = NF_INQ_DIMID(NCID, 'mdl_GridNNPts', nsecptsID)
!		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
!	status = NF_INQ_DIMLEN(NCID, nsecptsID, nsecpts)
!		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
!
!	npts = nsec*nsecpts
!
!!Get Topography
!	status = NF_INQ_VARID(NCID, 'mdl_GridTopo', mdl_GridTopoID)
!		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
!	ALLOCATE(tmptopo(npts), STAT = status)
!	status = NF_GET_VAR_REAL(NCID, mdl_GridTopoID, tmptopo)
!		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
!!Get Grid X and Y
!	status = NF_INQ_VARID(NCID, 'mdl_GridX', mdl_GridXID)
!		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
!	ALLOCATE(tmpx(npts), STAT = status)
!	status = NF_GET_VAR_REAL(NCID, mdl_GridXID, tmpx)
!		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
!
!	status = NF_INQ_VARID(NCID, 'mdl_GridY', mdl_GridYID)
!		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
!	ALLOCATE(tmpy(npts), STAT = status)
!	status = NF_GET_VAR_REAL(NCID, mdl_GridYID, tmpy)
!		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
!! temporary solution to cd
!	status = NF_GET_ATT_REAL(NCID, NF_GLOBAL, 'r2d_HydAttCD', tmpcd1)
!		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
!
!	status = NF_GET_ATT_INT(NCID, NF_GLOBAL, 'r2d_HydAttCDType', tmpcdtype)
!		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
!	
!	status = NF_GET_ATT_REAL(NCID, NF_GLOBAL, 'r1d_Disch', q)
!		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
!		
!	status = NF_GET_ATT_REAL(NCID, NF_GLOBAL, 'r1d_WSElev', wselev)
!		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
!
!
!	IF (tmpcdtype.ne.0) THEN		
!	ELSE
!	
!
!		status = NF_INQ_VARID(NCID, 'mdl_GridCD', mdl_GridCDID)
!			IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
!		ALLOCATE(tmpcd(npts), STAT = status)
!		status = NF_GET_VAR_REAL(NCID, mdl_GridCDID, tmpcd)
!			IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
!		
!	ENDIF
!
!
!	CALL AllocateSecs(1, nsec, nsec, npts, 0, 0, 2)
!	CALL AllocateTBC(2)
!	ellim = -1.e9
!	nchx = 1
!	ns = nsec
!	chbc(1,1) = 2
!	chbc(1,2) = 1
!	nss = 1
!	nps = 1
!	nsf=nss+ns-1
!	nsce(1)=nsf
!
!	do i = 1,nsec
!		counti = ((i-1)*nsecpts)
!		xsloc(i,1) = tmpx(counti+1)
!		xsloc(i,2) = tmpy(counti+1)
!		xsloc(i,3) = tmpx(counti+nsecpts)
!		xsloc(i,4) = tmpy(counti+nsecpts)
!		xsloc(i,5) = 0.
!		xsloc(i,6) = sqrt((xsloc(i,1)-xsloc(i,3))**2+(xsloc(i,2)-xsloc(i,4))**2) ! length
!
!			npt = nsecpts
!			ncd = 1
!			nnpt = nsecpts
!		npf=nps+npt-1
!		nscpte(i)=npf
!		emin=1.e9
!		if(i==1)eav = 0
!		do j = 1,npt
!			countij = ((i-1)*nsecpts) + j
!			jp = nps+j-1
!
!			if(j == 1) then
!				xspts(jp,1) = 0.
!			else
!				xspts(jp,1) = sqrt((xsloc(i,1)-tmpx(countij))**2 + &
!							(xsloc(i,2)-tmpy(countij))**2)
!			endif
!			xspts(jp,2) = tmptopo(countij)
!			emin = min(emin, xspts(jp,2))
!			if(i==1)eav = eav+(tmptopo(countij)-wselev)
!			
!			if(tmpcdtype.ne.0) then
!				xspts(jp,3) = tmpcd1
!			else
!				xspts(jp,3) = tmpcd(countij)
!			endif
!		
!		enddo
!		if(i==1)eav = eav/npt
!		do j = 1,npt
!			countij = ((i-1)*nsecpts) + j
!			jp = nps+j-1
!			
!			if(tmpcdtype.ne.0) then
!!				xspts(jp,3) = sqrt(tmpcd1*(eav)**.333/9.81);
!				xspts(jp,3) = tmpcd1;
!			else
!!				xspts(jp,3) = sqrt(tmpcd(countij)*(eav)**.333/9.81);
!				xspts(jp,3) = tmpcd(countij);
!			endif
!		
!		enddo
!
!		xsloc(i,5) = emin
!		if(emin>ellim(1)) then
!			ellim(1) = emin
!			limsec(1) = i
!		endif
!		nps = npf+1
!		nss = nsf+1
!
!	enddo
!	
!	nonds = 0
!	jctsecs = 0
!	maxbw = 1
!!	CALL SetBW(2*(maxbw+1))
!	nods = nsce(nchx)
!	DO i = 1, nods
!		CALL Farea(i, xsloc(i,5)+etol)
!	ENDDO
!	simstrt = 0
!	simend = 0
!
!	ntbc = 0
!	tsbcloc = 0
!
!	DO i = 1, nchx
!		DO j = 1,2
!		if(chbc(i,j)==1.or.chbc(i,j)==2) then
!			tsbcloc(i,2*j-1) = ntbc+1
!			DO k = 1,1
!				jp = ntbc + k
!				tmsrbc(jp, 1) = 0
!				if(j==1) then
!					tmsrbc(jp, 2) = q
!				else
!					tmsrbc(jp, 2) = wselev
!				endif
!			ENDDO
!			ntbc = ntbc+1
!			tsbcloc(i,2*k) = ntbc
!		endif
!		ENDDO
!	ENDDO
!	
!	DEALLOCATE(tmptopo, STAT = status)
!	DEALLOCATE(tmpx, STAT = status)
!	DEALLOCATE(tmpy, STAT = status)
!	DEALLOCATE(tmpcd, STAT = status)
!
!	status = NF_CLOSE(NCID)
!		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
!END SUBROUTINE read_input
	
!SUBROUTINE HANDLE_ERR(STATUS)
!	INCLUDE "netcdf.inc"
!	INTEGER STATUS
!	IF (STATUS .NE. NF_NOERR) THEN
!		PRINT *, NF_STRERROR(STATUS)
!		STOP 'Stopped'
!	ENDIF
!END SUBROUTINE HANDLE_ERR

END
