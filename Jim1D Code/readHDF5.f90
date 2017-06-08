MODULE readHDF5Mod
USE HDF5
USE HDF5UtilsMod
USE Support
IMPLICIT NONE

CHARACTER(LEN=3) :: g_Index, g_Index2
INTEGER :: numCh, numJct
!ALLOCATION DIMENSIONS
INTEGER :: totNumXSects, TotNumXSectPts, MaxNumXSects

CONTAINS

SUBROUTINE readHDF5(InputFile)
	IMPLICIT NONE
	CHARACTER(*), INTENT(IN) :: InputFile
	INTEGER(HID_T) :: file_id
	INTEGER :: error
		CALL h5open_f(error)
			CALL h5fopen_f(InputFile, H5F_ACC_RDONLY_F, file_id, error)

				CALL initializeParameters(file_id)
				CALL initializeArrays(file_id)
				CALL initializeChannelData(file_id)
				CALL initializeJunctionData(file_id)

			CALL h5fclose_f(file_id, error)
		CALL h5close_f(error)

END SUBROUTINE readHDF5

SUBROUTINE initializeParameters(file_id)
	IMPLICIT NONE
	INTEGER(HID_T), INTENT(IN) :: file_id
	INTEGER(HID_T) :: n_id, p_id, et_id, bt_id, mp_id
	INTEGER(HID_T) :: attid, dtype
	INTEGER, DIMENSION(7) :: data_dims
	REAL, DIMENSION(9) :: params
	INTEGER :: error, nprm, k
	CALL getGroup(file_id, g_1DNetwork, n_id);
		CALL getGroup(n_id, g_Params, p_id);

			CALL h5aopen_name_f(p_id, "Sim Start", attid, error)
			CALL h5aget_type_f(attid, dtype, error)
			data_dims(1) = 1
			CALL h5aread_f(attid, dtype, simstrt, data_dims, error)
			CALL h5aclose_f(dtype, error)
			CALL h5aclose_f(attid, error)

			CALL h5aopen_name_f(p_id, "Sim End", attid, error)
			CALL h5aget_type_f(attid, dtype, error)
			data_dims(1) = 1
			CALL h5aread_f(attid, dtype, simend, data_dims, error)
			CALL h5aclose_f(dtype, error)
			CALL h5aclose_f(attid, error)

			CALL h5aopen_name_f(p_id, "Model Parameters", attid, error)
			CALL h5aget_type_f(attid, dtype, error)
			data_dims(1) = 1
			CALL h5aread_f(attid, dtype, params, data_dims, error)
			CALL h5aclose_f(dtype, error)
			CALL h5aclose_f(attid, error)

	  nprm=8
	  prmt=0.
	  prmt(1,1)=params(1) !courant time step multiplier
	  prmt(1,2)=params(2) !dry-wet channel tolerance, m
	  prmt(1,3)=params(3)  !theta, advanced time step weight in Preissman scheme
	  prmt(1,4)=params(4)  !Number of points in time series plots, max=200
	  prmt(1,5)=params(5) !Mass concentration tolerance level
	  prmt(1,6)=params(6)  !Maximum number iterations, preissman scheme
	  prmt(1,7)=params(7)  !Relaxation coefficient--discharge, preissman scheme
	  prmt(1,8)=params(8)  !Relaxation coefficient--stage, preissman scheme
	  dragtype = int(params(9)) !Drag type variable
		corstpmlt=prmt(1,1)  ! min courant time step multiplier for setting dt used in time stepping
		etol=prmt(1,2)
		thet=prmt(1,3)	! time step weighting factor applied to future step in iterative solution
		itersmx=int(prmt(1,6))
		voltol=prmt(1,5)
		rlxdschg=prmt(1,7)
		rlxstg=prmt(1,8)
!	  nodebalindx=1   ! run nodebalance during simulation if index==1
!	  indxresist=0    ! use user specified alluvail n values rather than program calculated   
		do k=1,nprm
			prmt(3,k)=prmt(1,k)
			prmt(2,k)=prmt(1,k)
		end do

END SUBROUTINE

SUBROUTINE initializeJunctionData(file_id)
	IMPLICIT NONE
	INTEGER(HID_T), INTENT(IN) :: file_id
	INTEGER(HID_T) :: gj_id, j_id, n_id
	INTEGER :: numJcts, n
	INTEGER :: maxbw
		maxbw = 1
		nochjcts=0
		qsign=1.
		CALL getGroup(file_id, g_1DNetwork, n_id)
			CALL getGroup(n_id, g_Juncts, gj_id)

				CALL getNumJunctions(gj_id, numJcts)
				DO n = 1, numJcts
					CALL setJunctionValues(gj_id, n)
				ENDDO
				IF(numJcts > 0) THEN
					CALL setJunctionSections(numJcts)
					CALL findBandWidth(maxbw)
				ENDIF
!				CALL SetBW(2*(maxbw+1))
				Do n = 1, numJcts
					CALL Farea(n, xsloc(n,5)+etol)
				END DO

			CALL closeGroup(gj_id)
		CALL closeGroup(n_id)

END SUBROUTINE

SUBROUTINE setJunctionValues(gj_id, n)
	IMPLICIT NONE
	INTEGER(HID_T), INTENT(IN) :: gj_id
	INTEGER, INTENT(IN) :: n
	INTEGER(HID_T) :: j_id
	INTEGER :: index, numJCh
	INTEGER :: status, j
	INTEGER, ALLOCATABLE, DIMENSION(:) :: jdata
	REAL :: initval
		CALL getIndexGroup(n, gj_id, g_Junct, j_id)
			CALL getNumJChannels(j_id, numJCh)
			ALLOCATE(jdata(numJCh), STAT = status)
			CALL getJunctionChannels(j_id, jdata)

				nprjct(n) = numJCh
				DO j = 1, numJCh
					jctentrs(n,j) = jdata(j)
				ENDDO
				do j=1,nprjct(n)
					if(jctentrs(n,j)<0)then
						qsign(n,j)=-1.
						jctentrs(n,j)=-jctentrs(n,j)
						nochjcts(jctentrs(n,j),1)=n
					else
						nochjcts(jctentrs(n,j),2)=n
					endif
				enddo
			CALL getJInitValue(j_id, initval)
			jctelev(n) = initval
			DEALLOCATE(jdata, STAT=status)
		CALL closeGroup(j_id)
END SUBROUTINE

SUBROUTINE getJunctionChannels(j_id, jdata)
	IMPLICIT NONE
	INTEGER(HID_T), INTENT(IN) :: j_id
	INTEGER, DIMENSION(:) :: jdata
	INTEGER(HID_T) :: attid, dtype
	INTEGER :: error
	INTEGER, DIMENSION(7) :: data_dims
	CHARACTER(LEN=11) :: jstring = 'jctChannels'

		CALL h5aopen_name_f(j_id, jstring, attid, error)
		CALL h5aget_type_f(attid, dtype, error)
		data_dims(1) = 1
		CALL h5aread_f(attid, dtype, jdata, data_dims, error)
		CALL h5aclose_f(dtype, error)
		CALL h5aclose_f(attid, error)

END SUBROUTINE

SUBROUTINE setJunctionSections(nonds)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: nonds
	INTEGER :: ic, nus, nds, nsec, ibc, n, nin, nc, jc, ncc
	INTEGER :: jcc, nloc, indx, nu, nd, nss
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
END SUBROUTINE

SUBROUTINE findBandWidth(maxbw)
	IMPLICIT NONE
	INTEGER, INTENT(INOUT) :: maxbw
	INTEGER :: ic, ibc, jmax, jmin, nos

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

END SUBROUTINE

SUBROUTINE initializeChannelData(file_id)
	IMPLICIT NONE
	INTEGER(HID_T), INTENT(IN) :: file_id
	INTEGER(HID_T) :: gn_id, gch_id, ch_id
	INTEGER :: n, error
	INTEGER :: nss, nps
	INTEGER :: nch
	INTEGER :: ntbc
	ntbc = 0
	nss = 1; nps = 1; 
	ellim=-1.e9
	limsec=0
	xsloc=0.

	CALL getGroup(file_id, g_1DNetwork, gn_id)
		CALL getGroup(gn_id, g_Channels,  gch_id)
		CALL getNumChannels(gch_id, nch)
		DO n = 1, nch
			CALL getIndexGroup(n, gch_id, g_Channel, ch_id)
				CALL setChannelBCValues(n, ch_id, ntbc)
				CALL setChannelValues(n, ch_id, nss, nps)
			CALL closeGroup(ch_id)
		END DO
		CALL closeGroup(gch_id)
	CALL closeGroup(gn_id)

END SUBROUTINE

SUBROUTINE setChannelValues(chan_num, g_id, nss, nps)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: chan_num
	INTEGER(HID_T), INTENT(IN) :: g_id
	INTEGER, INTENT(INOUT) :: nss, nps
	INTEGER(HID_T) :: ch_id
	INTEGER :: ns, i
	INTEGER :: nsf
	INTEGER :: count = 0
		count = 0
		CALL getNumXSections(g_id, ns)
		nsf=nss+ns-1
		nsce(chan_num)=nsf
		DO i = nss,nsf
			count = count+1
			CALL setXSectionValues(g_id, chan_num, i, count, nps)

		END DO
		nss = nsf+1
		
END SUBROUTINE

SUBROUTINE setXSectionValues(ch_id, nch, chindex, sindex, nps)
	IMPLICIT NONE
	INTEGER(HID_T), INTENT(IN) :: ch_id
	INTEGER, INTENT(IN) :: nch, chindex, sindex
	INTEGER, INTENT(INOUT) :: nps
	INTEGER(HID_T) :: xs_id
	INTEGER :: status, j, jp
	INTEGER :: numXSPts, numRPts, ncd, nnpt
	INTEGER :: npf, jj
	REAL :: emin
	REAL, ALLOCATABLE, DIMENSION(:) :: x, y, dist, rough, topo
		CALL getIndexGroup(sindex, ch_id, g_XSection, xs_id)
			CALL getNumXSectPts(xs_id, numXSPts)
			CALL getNumRoughPts(xs_id, numRPts)
				npt = numXSPts
				ncd = 1
				nnpt = numRPts

				npf=nps+npt-1
				nscpte(chindex)=npf
				emin=1.e9
				ALLOCATE(x(numXSPts), STAT = status)
				ALLOCATE(y(numXSPts), STAT = status)
				ALLOCATE(dist(numXSPts), STAT = status)
				ALLOCATE(topo(numXSPts), STAT = status)
				ALLOCATE(rough(numRPts), STAT = status)
				CALL getXSectionData(xs_id, 'Roughness', rough)
				CALL getXSectionData(xs_id, 'Topography', topo)
				CALL getXSectionData(xs_id, 'SectX', x)
				CALL getXSectionData(xs_id, 'SectY', y)
				CALL getXSectionData(xs_id, 'Distance', dist)
		CALL closeGroup(xs_id)
		xsloc(chindex,1) = x(1); xsloc(chindex,2) = y(1)
		xsloc(chindex,3) = x(numXSPts); xsloc(chindex,4) = y(numXSPts)
		xsloc(chindex,6) = dist(numXSpts)
		DO j=1,npt
			jp=nps+j-1
			xspts(jp,1)=dist(j)
			xspts(jp,2)=topo(j)
			emin=min(emin,xspts(jp,2))
		END DO
		xsloc(chindex,5)=emin  ! xsloc( ,5) set to min ch elev
		if(emin>ellim(nch))then
			ellim(nch)=emin
			limsec(nch)=chindex
		endif
		jj=0
		do j=nps,npf
			jj=jj+1
			if(jj<=nnpt)then
			xspts(j,3)=rough(jj)
			else
			xspts(j,3)=rough(nnpt)
			endif
		enddo
		nps = npf+1
		DEALLOCATE(x, STAT = status)
		DEALLOCATE(y, STAT = status)
		DEALLOCATE(dist, STAT = status)
		DEALLOCATE(topo, STAT = status)
		DEALLOCATE(rough, STAT = status)
END SUBROUTINE

SUBROUTINE getXSectionData(xs_id, string, xsdata)
	IMPLICIT NONE
	INTEGER(HID_T), INTENT(IN) :: xs_id
	CHARACTER(LEN=*), INTENT(IN) :: string
	REAL, DIMENSION(:) :: xsdata
	INTEGER(HID_T) :: dset_id, dspace_id
	INTEGER :: error
	INTEGER, DIMENSION(7) :: data_dims

		CALL h5dopen_f(xs_id, string, dset_id, error)
		CALL h5dget_type_f(dset_id, dspace_id, error)
		CALL h5dextend_f(dset_id, data_dims, error)
		CALL h5dread_f(dset_id, dspace_id, xsdata, data_dims, error)
		CALL h5dclose_f(dspace_id, error)
		CALL h5dclose_f(dset_id, error)

END SUBROUTINE

SUBROUTINE setChannelBCValues(chindex, ch_id, ntbc)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: chindex
	INTEGER(HID_T), INTENT(IN) :: ch_id
	INTEGER, INTENT(INOUT) :: ntbc
	INTEGER(HID_T) :: bc_id, usbc_id, dsbc_id
	INTEGER :: bccode, nbcs, j, jp
	REAL :: ssval, tsimend
	CHARACTER(LEN=120) :: ts
		CALL getGroup(ch_id, g_BoundCond, bc_id)
			CALL getGroup(bc_id, g_USBC, usbc_id)
			CALL getBCValues(usbc_id, bccode, ssval, ts, nbcs)
				if(bccode <= 3) then
					chbc(chindex,1) = bccode
				else
					chbc(chindex,1) = 0
				end if
				if(bccode == 1 .or. bccode == 2) then
					tsbcloc(chindex, 1) = ntbc+1
					if(simend > 0.) then
						open(3, file=ADJUSTL(ts))
						read(3,*) nbcs
						DO j = 1,nbcs
							jp = ntbc+j
							read(3,*) tmsrbc(jp,1), tmsrbc(jp,2)
						END DO
						ntbc = ntbc + nbcs
						tsbcloc(chindex, 2) = ntbc
					else
						nbcs = 1
						jp = ntbc+1
						tmsrbc(jp, 1) = 0
						tmsrbc(jp, 2) = ssval
						ntbc = ntbc + nbcs
						tsbcloc(chindex, 2) = ntbc
					endif
				END IF
			CALL closeGroup(usbc_id)

			CALL getGroup(bc_id, g_DSBC, dsbc_id)
			CALL getBCValues(dsbc_id, bccode, ssval, ts, nbcs)
				IF(bccode <= 3) THEN
					chbc(chindex,2) = bccode
				ELSE
					chbc(chindex,2) = 0
				END IF
				IF(bccode == 1 .or. bccode == 2) THEN
					tsbcloc(chindex, 3) = ntbc+1
					IF(simend > 0.) THEN
						OPEN(3, file=ADJUSTL(ts))
						READ(3,*) nbcs
						DO j = 1,nbcs
							jp = ntbc+j
							READ(3,*) tmsrbc(jp,1), tmsrbc(jp,2)
						END DO
						ntbc = ntbc + nbcs
						tsbcloc(chindex, 4) = ntbc
					ELSE
						nbcs = 1
						jp = ntbc+1
						tmsrbc(jp, 1) = 0
						tmsrbc(jp, 2) = ssval
						ntbc = ntbc + nbcs
						tsbcloc(chindex, 4) = ntbc
					END IF
				END IF
			CALL closeGroup(dsbc_id)
		CALL closeGroup(bc_id)

END SUBROUTINE

SUBROUTINE initializeArrays(file_id)
	IMPLICIT NONE
	INTEGER(HID_T), INTENT(IN) :: file_id
	INTEGER(HID_T) :: gn_id, gch_id, ch_id, gj_id
	INTEGER :: i
	INTEGER :: tmpnxsect, tmpmnxsect, tmpnumpts
	INTEGER :: tmpnumbcts, tmpnumbctspts
	INTEGER :: tmpmnchjct

		tmpnxsect = 0; tmpnumpts = 0
		tmpmnxsect = -1e9; tmpmnchjct = -1e9
		tmpnumbcts = 0; tmpnumbctspts = 0

		CALL getGroup(file_id, g_1DNetwork, gn_id)
			CALL getGroup(gn_id, g_Channels, gch_id)
				CALL getNumChannels(gch_id, numCh)
				DO i = 1,numCH
					CALL getChannelDimensions(gch_id, i, tmpnxsect, tmpmnxsect, tmpnumpts, &
												tmpnumbcts, tmpnumbctspts)
				END DO
			CALL closeGroup(gch_id)

			CALL getGroup(gn_id, g_Juncts, gj_id)
				CALL getNumJunctions(gj_id, numJct)
				DO i = 1,numJct
					CALL getJunctionDimensions(gj_id, i, tmpmnchjct)
				END DO
			CALL closeGroup(gj_id)

			CALL AllocateSecs(numCh, tmpmnxsect, tmpnxsect, tmpnumpts, numJct, tmpmnchjct, tmpnumbcts)
			CALL AllocateTBC(tmpnumbctspts)

		CALL closeGroup(gn_id)

END SUBROUTINE initializeArrays

SUBROUTINE getJunctionDimensions(gn_id, index, mnchjct)
	IMPLICIT NONE
	INTEGER(HID_T), INTENT(IN) :: gn_id
	INTEGER, INTENT(IN) :: index
	INTEGER, INTENT(INOUT) :: mnchjct
	INTEGER(HID_T) :: j_id
	INTEGER :: error
	INTEGER :: i
	INTEGER :: tmpnjct=0
		tmpnjct=0
		CALL getIndexGroup(index, gn_id, g_Junct, j_id)
			CALL getNumJChannels(j_id, tmpnjct)
			mnchjct = MAX(mnchjct, tmpnjct)
		CALL closeGroup(j_id)
END SUBROUTINE

SUBROUTINE getNumJunctions(gj_id, numJnct)
	IMPLICIT NONE
	INTEGER(HID_T), INTENT(IN) :: gj_id
	INTEGER, INTENT(OUT) :: numJnct
	INTEGER(HID_T) :: attid, dtype
	INTEGER, DIMENSION(7) :: data_dims
	INTEGER :: error
		CALL h5aopen_name_f(gj_id, "numJunctions", attid, error)
		CALL h5aget_type_f(attid, dtype, error)
		data_dims(1) = 1
		CALL h5aread_f(attid, dtype, numJnct, data_dims, error)
		CALL h5aclose_f(dtype, error)
		CALL h5aclose_f(attid, error)
END SUBROUTINE

SUBROUTINE getNumJChannels(j_id, numJChn)
	IMPLICIT NONE
	INTEGER(HID_T), INTENT(IN) :: j_id
	INTEGER, INTENT(OUT) :: numJChn
	INTEGER(HID_T) :: attid, dtype
	INTEGER, DIMENSION(7) :: data_dims
	INTEGER :: error
		CALL h5aopen_name_f(j_id, "numChannels", attid, error)
		CALL h5aget_type_f(attid, dtype, error)
		data_dims(1) = 1
		CALL h5aread_f(attid, dtype, numJChn, data_dims, error)
		CALL h5aclose_f(dtype, error)
		CALL h5aclose_f(attid, error)
END SUBROUTINE

SUBROUTINE getJInitValue(j_id, initval)
	IMPLICIT NONE
	INTEGER(HID_T), INTENT(IN) :: j_id
	REAL, INTENT(OUT) :: initval
	INTEGER(HID_T) :: attid, dtype
	INTEGER, DIMENSION(7) :: data_dims
	INTEGER :: error
		CALL h5aopen_name_f(j_id, "initVal", attid, error)
		CALL h5aget_type_f(attid, dtype, error)
		data_dims(1) = 1
		CALL h5aread_f(attid, dtype, initval, data_dims, error)
		CALL h5aclose_f(dtype, error)
		CALL h5aclose_f(attid, error)
END SUBROUTINE


SUBROUTINE getBCValues(bc_id, code, ssVal, TS, nbcs)
	IMPLICIT NONE
	INTEGER(HID_T), INTENT(IN) :: bc_id
	INTEGER, INTENT(INOUT) :: code, nbcs
	REAL, INTENT(INOUT) :: ssVal
	CHARACTER(LEN=*), INTENT(INOUT) ::TS
	INTEGER(HID_T) :: attid, atype 
	INTEGER :: error
	INTEGER, DIMENSION(7) :: data_dims
		CALL h5aopen_name_f(bc_id, "Code", attid, error)
		CALL h5aget_type_f(attid, atype, error)
		data_dims(1) = 1
		CALL h5aread_f(attid, atype, code, data_dims, error)
		CALL h5aclose_f(atype, error)
		CALL h5aclose_f(attid, error)

		CALL h5aopen_name_f(bc_id, "Steady State Value", attid, error)
		CALL h5aget_type_f(attid, atype, error)
		data_dims(1) = 1
		CALL h5aread_f(attid, atype, ssVal, data_dims, error)
		CALL h5aclose_f(atype, error)
		CALL h5aclose_f(attid, error)

		CALL h5aopen_name_f(bc_id, "Time Series", attid, error)
		CALL h5aget_type_f(attid, atype, error)
		data_dims(1) = 1
		CALL h5aread_f(attid, atype, TS, data_dims, error)
		CALL h5aclose_f(atype, error)
		CALL h5aclose_f(attid, error)
		IF(simend > 0 .and. code < 10) then
			OPEN(3, file = ADJUSTL(ts))
			READ(3,*) nbcs
			CLOSE(3)
		ELSE
			nbcs = 1
		END IF


END SUBROUTINE




SUBROUTINE getNumXSectPts(xs_id, numXSPts)
	IMPLICIT NONE
	INTEGER(HID_T), INTENT(IN) :: xs_id
	INTEGER, INTENT(OUT) :: numXSPts
	INTEGER(HID_T) :: attid, dtype
	INTEGER, DIMENSION(7) :: data_dims
	INTEGER :: error
		CALL h5aopen_name_f(xs_id, "numSectPts", attid, error)
		CALL h5aget_type_f(attid, dtype, error)
		data_dims(1) = 1
		CALL h5aread_f(attid, dtype, numXSPts, data_dims, error)
		CALL h5aclose_f(dtype, error)
		CALL h5aclose_f(attid, error)
END SUBROUTINE

SUBROUTINE getNumRoughPts(xs_id, numRPts)
	IMPLICIT NONE
	INTEGER(HID_T), INTENT(IN) :: xs_id
	INTEGER, INTENT(OUT) :: numRPts
	INTEGER(HID_T) :: attid, dtype
	INTEGER, DIMENSION(7) :: data_dims
	INTEGER :: error
		CALL h5aopen_name_f(xs_id, "numRoughPts", attid, error)
		CALL h5aget_type_f(attid, dtype, error)
		data_dims(1) = 1
		CALL h5aread_f(attid, dtype, numRPts, data_dims, error)
		CALL h5aclose_f(dtype, error)
		CALL h5aclose_f(attid, error)
END SUBROUTINE

SUBROUTINE getChannelDimensions(g_id, index, numXSects, maxNumSect, numXSectPts, numBCTS, numBCTSPts)
	IMPLICIT NONE
	INTEGER(HID_T), INTENT(IN) :: g_id	!This is the Channel_INDEX group id
	INTEGER, INTENT(IN) :: index	!This the Channel number
	INTEGER, INTENT(INOUT) :: numXSects, maxNumSect, numXSectPts, numbcts, numbctspts
	INTEGER(HID_T) :: ch_id, xs_id
	INTEGER :: error
	INTEGER :: i
	INTEGER :: tmpnxsectpts=0, tmpnxsects=0, tmpnbcpts=0, tmpnbc=0
	tmpnxsectpts=0; tmpnxsects=0; tmpnbcpts=0; tmpnbc=0;
		CALL getIndexGroup(index, g_id, g_Channel, ch_id)
			CALL getBCDimensions(ch_id, tmpnbc, tmpnbcpts)
				numBCTS = numBCTS + tmpnbc
				numBCTSPts = numBCTSPts + tmpnbcpts
			CALL getNumXSections(ch_id, tmpnxsects)
			numXSects = numXSects+tmpnxsects
			DO i = 1,tmpnxsects
				CALL getIndexGroup(i, ch_id, g_XSection, xs_id)
				CALL getNumXSectPts(xs_id, tmpnxsectpts)
					numXSectPts = numXSectPts+tmpnxsectpts
					maxNumSect = MAX(maxNumSect, tmpnxsectpts)
				CALL closeGroup(xs_id)
			END DO
		CALL closeGroup(ch_id)
END SUBROUTINE

SUBROUTINE getBCDimensions(ch_id, numbctimeser, numbcpts)
	IMPLICIT NONE
	INTEGER(HID_T), INTENT(IN) :: ch_id
	INTEGER, INTENT(INOUT) :: numbctimeser, numbcpts
	INTEGER(HID_T) :: bc_id, usbc_id, dsbc_id
	INTEGER :: bccode=-1, nbcs=0
	REAL :: ssval
	CHARACTER(LEN=120) :: ts
		CALL getGroup(ch_id, g_BoundCond, bc_id)

			CALL getGroup(bc_id, g_USBC, usbc_id)
			nbcs = 0
			CALL getBCValues(usbc_id, bccode, ssval, ts, nbcs)
				IF(bccode < 10) then
					numbctimeser = numbctimeser+1
				END IF
				numbcpts = numbcpts + nbcs
			CALL closeGroup(usbc_id)

			CALL getGroup(bc_id, g_DSBC, dsbc_id)
			nbcs = 0
			CALL getBCValues(dsbc_id, bccode, ssval, ts, nbcs)
				IF(bccode < 10) then
					numbctimeser = numbctimeser+1
				END IF
				numbcpts = numbcpts + nbcs
			CALL closeGroup(dsbc_id)

		CALL closeGroup(bc_id)

END SUBROUTINE


END MODULE readHDF5Mod
