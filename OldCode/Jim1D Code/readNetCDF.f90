MODULE readNetCDFMod
USE Support
IMPLICIT NONE
INCLUDE "netcdf.inc"
CONTAINS


	SUBROUTINE read_NetCDFInput(InputFile)

	CHARACTER(*), INTENT(IN) :: InputFile
	INTEGER :: status, NCID, tmpint
	REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar, tmpvar2
	INTEGER, ALLOCATABLE, DIMENSION(:) :: tmpvari,numptsinxsect
	REAL, ALLOCATABLE, DIMENSION(:) ::  geopts, xnpts
	
!	INTEGER, Dimension(1)::numCh
	INTEGER::numCh
	INTEGER::netID, topoDimID, topoVarID
	INTEGER::BCDimID, BCVarID, jctID, jctDimID, jctVarID
	INTEGER::RghDimID, rghID
	INTEGER::numTopoPts, numXSects, maxelem
	INTEGER::chindex, bcindex, jctindex, tmpindex
	INTEGER::maxsectchan, totnumxsect, totnumgeopts 
	INTEGER::numjct, totnumjct, numchjct, maxchanjct, numbcts, numbc, totnumbc, tmp
	INTEGER::nnptsch
	!for topography
	INTEGER::nss, nps, npsch, npstot, n, k, ns, nsf, ict, i, numsecpts, ngptend, ngptbgn, ncd
	INTEGER::nnpt, npf, j, jp
	REAL::x, y, z, xmnn, zz, xx, yy, min, emin
	!for junction
	INTEGER::nonds, maxbw, ic, nus, nds, nsec, ibc, nin, nc, jc, ncc, jcc, nloc, indx, nu
	INTEGER::nd, jmax, jmin, nos, nods
	!for bcs
	INTEGER::ntbc, nbcs
	INTEGER::textlen, tnbcs
	!for step-backwater calculations
	REAL:: tDisch, tStage, tbc1, tbc2
	REAL:: tcd
	INTEGER:: tcdtype
	INTEGER :: slen

	character(len = 14) strTopoDim
	character(len = 14) strTopo
	character(len = 12) strBCDim
	character(len = 6) strBC
	character(len = 15) strUSTS
	character(len = 15) strDSTS
	character(len = 18) strJctDim
	character(len = 12) strJct
	character(len = 15) strRoughDim
	character(len = 13) strRough
	character(len = 3) indstr 
	character(len = 3) trmstr
	character(40) tmpDimstr
	character(40) tmpVarstr
	character(225) tmpTSstr

	strTopoDim = '1D_NumTopoPts_'
	strTopo = '1D_Topography_'
	strBCDim = '1D_NumBCPts_'
	strBC = '1D_BC_'
	strJct = '1D_Junction_'
	strJctDim = '1D_NumJctChannels_'
	strRoughDim = '1D_NumRoughPts_'
	strRough = '1D_Roughness_'
	strUSTS = '1D_BC_US_TSFile'
	strDSTS = '1D_BC_DS_TSFile'
	maxsectchan = -1e6
	maxchanjct = -1e6
	totnumxsect = 0
	totnumgeopts = 0
	totnumbc = 0
	numbcts = 0
	nsf = 0
	jp = 0
	npsch = 0
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
			  
!!!!!!!!!!!!!! should be <3
!!!!!!!!!!!!!! should not use Rich's boundary code 4=slope given as of this version

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
			DO tmpindex = 1,numchjct
				maxchanjct = max(maxchanjct, abs(tmpvari(tmpindex)))
			ENDDO
			IF (maxchanjct < 0) maxchanjct = 1
!			maxchanjct = max(maxchanjct, maxval(tmpvari))
			DEALLOCATE(tmpvari)
		ENDDO
  
  
  	CALL AllocateSecs(numCh, maxsectchan, totnumxsect, totnumgeopts, numjct, maxchanjct, numbcts)
	tmp = 0


	ellim=-1.e9
	limsec=0
	nss=1
	nps=1
	xsloc=0.

! Rich, the geopts array has to be opened here--I've called it geopts() below
  ! I assume it has 3*totnumgeopts entries in x-y-z order

! Rich, the n-value array has to be opened here--I've called it xnpts() below
  ! I assume it has totnumgeopts entries
  
  ! it would be a hell of a lot easier for me if geopts and xnpts  
  ! could be opened and allocated or at least included in Support 
  !  --I've assumed this to be the case below  
	npstot = 1
	do n=1,numCh  !  Channel loop
		npstot = npstot + npsch
		npsch = 1
		WRITE(indstr,'(I3)') n-1
!
!Get Boundary Conditions
!
		indstr = ADJUSTL(indstr)
		tmpDimstr = strBCDim//TRIM(indstr)
		tmpVarstr = strBC//TRIM(indstr)
		status = NF_INQ_DIMID(NCID, TRIM(tmpDimstr), BCDimID)
		status = NF_INQ_DIMLEN(NCID, BCDimID, numbc)  !numbc should always be 2 is this necessary
		status = NF_INQ_VARID(NCID, TRIM(tmpVarstr), BCVarID)
		ALLOCATE(tmpvari(numbc), STAT = status)
		status = NF_GET_VAR_INT(NCID, BCVARID, tmpvari)
			do k=1,2
				chbc(n,k)=tmpvari(k)
				if(chbc(n,k)>3)chbc(n,k)=0 !MD_SWMS writes codes 0-5 4 is internal stage
										   !                         5 is internal discharge
			enddo	 
		DEALLOCATE(tmpvari)
!
! chbc=0 node stage
! chbc=1 external stage
! chbc=2 external discharge
! chbc=3 critical depth
!		
!
!Get Channel Topography
!
		tmpDimstr = strTopoDim//TRIM(indstr)
		tmpVarstr = strTopo//TRIM(indstr)
		status = NF_INQ_DIMID(NCID, TRIM(tmpDimstr), topoDimID)
		status = NF_INQ_DIMLEN(NCID, topoDimID, numTopoPts)
!		read(1,*)ns,chbc(n,1),chbc(n,2) ! number sections, us and ds bc types
		status = NF_INQ_VARID(NCID,TRIM(tmpVarstr), topoVarID) 
		status = NF_GET_ATT_INT(NCID, topoVarID, 'numXSections', ns)
		ALLOCATE(numptsinxsect(ns), STAT = status)
		status = NF_GET_ATT_INT(NCID, topoVarID, 'numPtsXSections', numptsinxsect)
		ALLOCATE(geopts(numTopoPts), STAT = status)
		status = NF_GET_VAR_REAL(NCID, topoVarID, geopts)
!
!GET ROUGHNESS VALUES
!
		tmpDimstr = strRoughDim//TRIM(indstr)
		tmpVarstr = strRough//TRIM(indstr)
		status = NF_INQ_DIMID(NCID, TRIM(tmpDimstr), RghDimID)
		status = NF_INQ_DIMLEN(NCID, RghDimID, nnptsch)
		status = NF_INQ_VARID(NCID,TRIM(tmpVarstr), RghID)
		ALLOCATE(xnpts(nnptsch), STAT = status)
		status = NF_GET_VAR_REAL(NCID, RghID, xnpts)

		if(simend <= 0.) then
			status = nf_get_att_int(NCID, NF_GLOBAL, "1DSBW_CDTYPE", tcdtype);
			status = nf_get_att_real(NCID, NF_GLOBAL, "1DSBW_CD", tcd);
		endif


		nsf=nss+ns-1
		nsce(n)=nsf
		ict=0
		do i=nss,nsf  ! Sections in channel loop
			ict=ict+1

			numsecpts=numptsinxsect(ict)  !Rich, is this the number of points in the section?
			
			ngptend=npsch+numsecpts-1
			ngptbgn= npsch
			
	!		read(1,*) (xsloc(i,j),j=1,5) !xr,yr,xl,yl,elev adj

			call ReadGpts(ngptbgn,xsloc(i,1),xsloc(i,2),zz,xmnn)
			call ReadGpts(ngptend,xsloc(i,3),xsloc(i,4),zz,xmnn)
			! elev adj is ==0.
			xsloc(i,6)=sqrt((xsloc(i,1)-xsloc(i,3))**2+(xsloc(i,2)-xsloc(i,4))**2) ! length
				! of section base line

	!		read(1,*) npt,ncd,nnpt
			npt=numsecpts
			ncd=1 ! pts always entered as elevs
			nnpt=npt ! number on manning n entries

			npf=nps+npt-1
			nscpte(i)=npf
			emin=1.e9

	!		else  !elevs are entered directly
	!			read(1,*)(dum(j),j=1,2*npt)
				do j=1,npt
	!	jp=nps+j-1
					jp = npstot+npsch+j-1-n
					call ReadGpts(npsch+j-1,xx,yy,zz,xmnn)
					xspts(jp,1)=sqrt((xx-xsloc(i,1))**2 + (yy-xsloc(i,2))**2)
					xspts(jp,2)=zz
					if(simend <= 0.) then
						xspts(jp,3)=tcd
					else
						xspts(jp,3)=xmnn
					endif
					emin=min(emin,xspts(jp,2))
					
				enddo
			xsloc(i,5)=emin  ! xsloc( ,5) set to min ch elev
			if(emin>ellim(n))then
				ellim(n)=emin
				limsec(n)=i
			endif

	!		read(1,*)(dum(j),j=1,nnpt)  ! Manning's n's set in distance-elev loop above
			nps=npf+1
			npsch = npsch+npt
			enddo  ! Sections in channel loop
			nss=nsf+1
			DEALLOCATE(geopts)
			DEALLOCATE(xnpts)
			DEALLOCATE(numptsinxsect)
	enddo  !  Channel loop
!
! START JUNTION INFORMATION
!
!Get Juntion params
 	!	read(1,*)nonds
	nonds=nnd

	jctsecs=0
	maxbw=1

	if(nonds>0)then

	nochjcts=0
	qsign=1.

	do n=1,nonds

!		read(1,*)nprjct(n)
!		read(1,*)(jctentrs(n,j),j=1,nprjct(n))

		jctindex=n
			WRITE(indstr,'(I3)') jctindex-1
			indstr = ADJUSTL(indstr)
			tmpDimstr = strJctDim//TRIM(indstr)
			tmpVarstr = strJct//TRIM(indstr)
			
			status = NF_INQ_DIMID(NCID, TRIM(tmpDimstr), jctDimID)
			status = NF_INQ_DIMLEN(NCID, jctDimID, numchjct)
		nprjct(n)=numchjct	
			status = NF_INQ_VARID(NCID, TRIM(tmpVarstr), jctVarID)
			status = NF_GET_ATT_REAL(NCID, jctVarID, "initVal", jctelev(n))
			ALLOCATE(tmpvari(numchjct), STAT = status)
			status = NF_GET_VAR_INT(NCID, jctVarID, tmpvari)
		do j=1,nprjct(n)
			jctentrs(n,j)=tmpvari(j)
		enddo
			DEALLOCATE(tmpvari)

		do j=1,nprjct(n)
			if(jctentrs(n,j)<0)then
				qsign(n,j)=-1.
				jctentrs(n,j)=-jctentrs(n,j)
				nochjcts(jctentrs(n,j),1)=n
			else
				nochjcts(jctentrs(n,j),2)=n
			endif
		enddo
	enddo
		do ic=1,nch  !find cross-sections that connect at nodes and place in jctsecs(ic,2,ndsmx) 
			call ChLims(ic,nus,nds)
			nsec=nus
			do ibc=1,2  ! do for both ends of channel
				if(ibc==2)nsec=nds
				if(chbc(ic,ibc)==0)then ! do only for channel ends at junctions
					do n=1,nonds  ! check all nodes for this channel (end)
						nin=nprjct(n)
						do nc=1,nin  ! look at each channel entering node n
							jc=jctentrs(n,nc)
!							print *, n,nc,jc,ic
							if(jc==ic.and.((ibc==1.and.qsign(n,nc)<0.).or.(ibc==2.and.qsign(n,nc)>0.)))then  ! the channel at node n matches ic
!							print *, n,nc,jc
								jctsecs(ic,ibc,1)=nin  ! no at junction
								jctsecs(ic,ibc,2)=nc  ! location ic channel
									do ncc=1,nin  ! find section numbers for other channels at jct
										jcc=jctentrs(n,ncc)
										nloc=ncc+2
										indx=int(qsign(n,ncc))
										call ChLims(jcc,nu,nd)
										nss=nu
										if(indx>0)nss=nd
										jctsecs(ic,ibc,nloc)=indx*nss
										jctsecs(ic,ibc,nloc+nin)=jcc
									enddo !ncc=1,nin  ! find section numbers for other channels at jct
							endif  !(jc==ic)then  ! the channel at node n matches ic
						enddo !nc=1,nin  ! look at each channel entering node n
					enddo !n=1,nonds  ! check all nodes for this channel (end)
				endif  !(chbc(ic,ibc)==0)then ! do only for channel ends at junctions
!		print *,ic,nsec
! 		write(2,*)ic,(jctsecs(ic,ibc,jj),jj=1,maxentrs+2)
			enddo !ibc=1,2  ! do for both ends of channel
		enddo !ic=1,nch
		do ic=1,nch  !find maxbandwidth
			do ibc=1,2
!		write(2,*)ic,(jctsecs(ic,ibc,jj),jj=1,2*maxentrs+2)
				jmax=0
				jmin=100000
				do nos=1,jctsecs(ic,ibc,1)
					jmax=max(jmax,abs(jctsecs(ic,ibc,nos+2)))
					jmin=min(jmin,abs(jctsecs(ic,ibc,nos+2)))
				enddo
			maxbw=max(maxbw,jmax-jmin)
!			write(2,*)jmax,jmin,maxbw
			enddo 
		enddo !ic=1,nch  !find maxbandwidth
	endif ! (nonds>0)then

!		call SetBW(2*(maxbw+1))

! set initial channel locations and areas
	nods=nsce(nch)
	do i=1,nods
		call Farea(i,xsloc(i,5)+etol)
	enddo
! Rich, start of boundary run time and boundary condition input info
! Pls supply Net CDF Patches to get this stuff

!	rewind(3)
!		read (3,*)simstrt,simend  ! simulation starting and ending times in days

	ntbc=0
	tsbcloc=0
	tnbcs=0
	do n=1,nch
		do k=1,2
!
!
! chbc=0 node stage
! chbc=1 external stage
! chbc=2 external discharge
! 
!
		WRITE(indstr,'(I3)') n-1
		indstr = ADJUSTL(indstr)
		tmpVarstr = strBC//TRIM(indstr)
		status = NF_INQ_VARID(NCID, TRIM(tmpVarstr), BCVarID)
		
			if(chbc(n,k)==1.or.chbc(n,k)==2)then
!				tsbcloc(n,2*k-1)=ntbc+1
					if(k == 1) then
						status = nf_inq_attlen(NCID, BCVarID, strUSTS, slen);
						if(slen > 0) then
							status = nf_get_att_text(NCID, BCVarID, strUSTS, tmpTSstr);
 							open (3,file=ADJUSTL(tmpTSstr))
							read(3,*)nbcs  
							close(3)
						else
							nbcs = 2
						endif
					else
						status = nf_inq_attlen(NCID, BCVarID, strDSTS, slen);
						if(slen > 0) then
							status = nf_get_att_text(NCID, BCVarID, strDSTS, tmpTSstr);
 							open (3,file=ADJUSTL(tmpTSstr))
							read(3,*)nbcs  
							close(3)
						else
							nbcs = 2
						endif
					endif
					tnbcs = tnbcs + nbcs
			endif
		enddo
	enddo
	CALL AllocateTBC(tnbcs)
	ntbc = 0
	do n = 1, nch
		do k = 1,2
		WRITE(indstr,'(I3)') n-1
		indstr = ADJUSTL(indstr)
		tmpVarstr = strBC//TRIM(indstr)
		status = NF_INQ_VARID(NCID, TRIM(tmpVarstr), BCVarID)
		
			if(chbc(n,k)==1.or.chbc(n,k)==2)then
				tsbcloc(n,2*k-1)=ntbc+1
				if(simend > 0.) then
					if(k == 1) then
						status = nf_get_att_text(NCID, BCVarID, strUSTS, tmpTSstr);
 						open (3,file=ADJUSTL(tmpTSstr))

					else
						status = nf_inq_attlen(NCID, BCVarID, strDSTS, slen);
						status = nf_get_att_text(NCID, BCVarID, strDSTS, tmpTSstr);
 						open (3,file=ADJUSTL(tmpTSstr))

					endif
					read(3,*)nbcs  
					do j=1,nbcs
						jp=ntbc+j
						read(3,*)tmsrbc(jp,1),tmsrbc(jp,2)
					enddo
					close(3)
					ntbc=ntbc+nbcs
					tsbcloc(n,2*k)=ntbc
				else	! will be step back water solution only
!					status = nf_get_att_real(NCID, NF_GLOBAL, "1DSBW_DISCH", tDisch);
!					status = nf_get_att_real(NCID, NF_GLOBAL, "1DSBW_STAGE", tStage);
					status = nf_get_att_real(NCID, BCVarID, "USSS_Value", tbc1)
					status = nf_get_att_real(NCID, BCVarID, "DSSS_Value", tbc2)
					if(k == 1) then
						nbcs = 1
						jp = ntbc+1
						tmsrbc(jp, 1) = 0
						tmsrbc(jp, 2) = tbc1
						ntbc = ntbc+nbcs
						tsbcloc(n,2*k) = ntbc
					else
						nbcs = 1
						jp = ntbc+1
						tmsrbc(jp, 1) = 0
						tmsrbc(jp, 2) = tbc2
						ntbc = ntbc+nbcs
						tsbcloc(n,2*k) = ntbc
					endif
				endif

						
			endif
		enddo
	enddo

! Locate and enter the initail ws elevations for the junctions --from a dialog probably
!	read(3,*)(jctelev(k),k=1,nnd)

! Ditto for parameters, remember that default values are set in the main program and 
! we may need to provide a link here or in main program to Rich's dialog
	status = NF_CLOSE(NCID)

	CONTAINS
		Subroutine ReadGpts(n,x,y,z,xmnn)
		IMPLICIT NONE
		INTEGER:: n
		REAL:: x, y, z, xmnn
			xmnn=xnpts(n)
			x=geopts(3*n-2)
			y=geopts(3*n-1)
			z=geopts(3*n)
		endSubroutine ReadGpts

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




