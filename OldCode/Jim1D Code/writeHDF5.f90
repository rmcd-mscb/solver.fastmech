MODULE writeHDF5Mod
USE HDF5
USE HDF5UtilsMod
USE Support
IMPLICIT NONE

CHARACTER(LEN=3) :: g_Index, g_Index2

CONTAINS

SUBROUTINE createHDF5Solutions(InputFile)
	IMPLICIT NONE
	CHARACTER(LEN=*), INTENT(IN) :: InputFile
	INTEGER(HID_T) :: file_id, dataspace, dataset, cparms
	INTEGER(HID_T) :: n_id, ch_id, ich_id
	INTEGER :: RANK = 2
	INTEGER :: error, n, numCh
		CALL h5open_f(error)
			CALL h5fopen_f(InputFile, H5F_ACC_RDWR_F, file_id, error)
				CALL getGroup(file_id, g_1DNetwork, n_id)
					CALL getGroup(n_id, g_Channels, ch_id)
						CALL createTimeSolution(ch_id)
						CALL getNumChannels(ch_id, numCh)
						DO n = 1, numCh
							CALL getIndexGroup(n, ch_id, g_Channel, ich_id)
								CALL createChannelSolutions(ich_id)

							CALL closeGroup(ich_id)
						END DO
					CALL closeGroup(ch_id)
				CALL closeGroup(n_id)
			CALL h5fclose_f(file_id, error)
		CALL h5close_f(error)

END SUBROUTINE

SUBROUTINE createTimeSolution(chs_id)
	IMPLICIT NONE
	INTEGER(HID_T), INTENT(IN) :: chs_id
	INTEGER(HID_T) :: chsg_id
	CALL createGroup(chs_id, g_Solutions, chsg_id)
			CALL createTimeSolDataSet(chsg_id, gs_Time)
	CALL closeGroup(chsg_id)
END SUBROUTINE


SUBROUTINE createChannelSolutions(ich_id)
	IMPLICIT NONE
	INTEGER(HID_T), INTENT(IN) :: ich_id
	INTEGER(HID_T) :: chsg_id, xssg_id, xs_id
	INTEGER :: error, numXS, i

		CALL getNumXSections(ich_id, numXS)
		CALL createGroup(ich_id, g_Solutions, chsg_id)

			CALL createCHSolDataSets(chsg_id, numXS)
			DO i = 1,numXS
				CALL getIndexGroup(i, ich_id, g_XSection, xs_id)
				CALL createXSSolDataSets(xs_id, numXS)
				CALL closeGroup(xs_id)

!! add xsect solution here???
			ENDDO
			CALL createXSSolDataSets(chsg_id, numXS)

		CALL closeGroup(chsg_id)

END SUBROUTINE

SUBROUTINE createXSSolDataSets(s_id, numXS)
	IMPLICIT NONE
	INTEGER(HID_T), INTENT(IN) :: s_id
	INTEGER, INTENT(IN) :: numXS
	INTEGER(HID_T) :: xss_id
	INTEGER :: RANK = 2
	INTEGER(HID_T), DIMENSION(2) :: maxdims, dims, chunkdims
	INTEGER :: error
		CALL createGroup(s_id, g_Solutions, xss_id)
		
			CALL createSolTimeDataSet_REAL(xss_id, 1, gs_RTWSX)
			CALL createSolTimeDataSet_REAL(xss_id, 1, gs_LTWSX)
			CALL createSolTimeDataSet_REAL(xss_id, 1, gs_RTWSY)
			CALL createSolTimeDataSet_REAL(xss_id, 1, gs_LTWSY)

		CALL closeGroup(xss_id)
		
END SUBROUTINE

SUBROUTINE createCHSolDataSets(s_id, numXS)
	IMPLICIT NONE
	INTEGER(HID_T), INTENT(IN) :: s_id
	INTEGER, INTENT(IN) :: numXS
	INTEGER(HID_T) :: chsol_id
	INTEGER :: RANK = 2
	INTEGER(HID_T), DIMENSION(2) :: maxdims, dims, chunkdims
	INTEGER :: error
		CALL createGroup(s_id, g_ChanSol, chsol_id)
		
			CALL createSolTimeDataSet_REAL(chsol_id, numXS, gs_WSE)
			CALL createSolTimeDataSet_REAL(chsol_id, numXS, gs_XSArea)
			CALL createSolTimeDataSet_REAL(chsol_id, numXS, gs_HydRadius)
			CALL createSolTimeDataSet_REAL(chsol_id, numXS, gs_HydDepth)
			CALL createSolTimeDataSet_REAL(chsol_id, numXS, gs_TopWidth)
			CALL createSolTimeDataSet_REAL(chsol_id, numXS, gs_Conveyance)
			CALL createSolTimeDataSet_REAL(chsol_id, numXS, gs_Discharge)
			CALL createSolTimeDataSet_REAL(chsol_id, numXS, gs_ShearStress)

		CALL closeGroup(chsol_id)
		
END SUBROUTINE

SUBROUTINE createTimeSolDataSet(chs_id, name)
	IMPLICIT NONE
	INTEGER(HID_T), INTENT(IN) :: chs_id
	CHARACTER(LEN=*) :: name

	INTEGER(HID_T) :: cparms, memspace, dataspace, dset_id
	INTEGER :: RANK = 1
	INTEGER(HID_T), DIMENSION(1) :: maxdims, dims, chunkdims
	INTEGER :: error

!		dims = (/numXS, 1/)
!		maxdims = (/numXS,H5S_UNLIMITED_F/)
		dims = (/1/)
		maxdims = (/H5S_UNLIMITED_F/)
		CALL h5screate_simple_f(RANK, dims, dataspace, error, maxdims)
		CALL h5pcreate_f(H5P_DATASET_CREATE_F, cparms, error)
		CALL h5pset_chunk_f(cparms, RANK, dims, error)
		CALL h5dcreate_f(chs_id, name, H5T_NATIVE_REAL, dataspace, dset_id, error, cparms)

		CALL h5sclose_f(dataspace, error)
		CALL h5dclose_f(dset_id, error)
		CALL h5pclose_f(cparms, error)


END SUBROUTINE

SUBROUTINE createSolTimeDataSet_REAL(chs_id, numXS, name)
	IMPLICIT NONE
	INTEGER(HID_T), INTENT(IN) :: chs_id
	INTEGER, INTENT(IN) :: numXS
	CHARACTER(LEN=*) :: name

	INTEGER(HID_T) :: cparms, memspace, dataspace, dset_id
	INTEGER :: RANK = 2
	INTEGER(HID_T), DIMENSION(2) :: maxdims, dims, chunkdims
	INTEGER :: error

!		dims = (/numXS, 1/)
!		maxdims = (/numXS,H5S_UNLIMITED_F/)
		if (numXS == 1) then 
			RANK = 1
		else 
			RANK =2 
		endif
		dims = (/1, numXS/)
		maxdims = (/H5S_UNLIMITED_F, numXS/)
		CALL h5screate_simple_f(RANK, dims, dataspace, error, maxdims)
		CALL h5pcreate_f(H5P_DATASET_CREATE_F, cparms, error)
		CALL h5pset_chunk_f(cparms, RANK, dims, error)
		CALL h5dcreate_f(chs_id, name, H5T_NATIVE_REAL, dataspace, dset_id, error, cparms)

		CALL h5sclose_f(dataspace, error)
		CALL h5dclose_f(dset_id, error)
		CALL h5pclose_f(cparms, error)


END SUBROUTINE
	

SUBROUTINE writeHDF5Output(InputFile)
	IMPLICIT NONE
	CHARACTER(*), INTENT(IN) :: InputFile
	INTEGER(HID_T) :: file_id, n_id, ch_id, ich_id, chs_id
	INTEGER :: numCH, n
	INTEGER :: error
	tstep = tstep+1
		CALL h5open_f(error)
			CALL h5fopen_f(InputFile, H5F_ACC_RDWR_F, file_id, error)
				CALL getGroup(file_id, g_1DNetwork, n_id)
					CALL getGroup(n_id, g_Channels, ch_id)
						CALL getGroup(ch_id, g_Solutions, chs_id)
						CALL writeSolTime(chs_id, gs_Time)
						CALL getNumChannels(ch_id, numCh)
						DO n = 1, numCh
							CALL getIndexGroup(n, ch_id, g_Channel, ich_id)

								CALL writeChannelSolutions(ich_id, n)

							CALL closeGroup(ich_id)
						END DO
						CALL closeGroup(chs_id)
					CALL closeGroup(ch_id)
				CALL closeGroup(n_id)
			CALL h5fclose_f(file_id, error)
		CALL h5close_f(error)

END SUBROUTINE writeHDF5Output

SUBROUTINE writeChannelSolutions(ich_id, chan_index)
	IMPLICIT NONE
	INTEGER(HID_T), INTENT(IN) :: ich_id
	INTEGER, INTENT(IN) :: chan_index
	INTEGER(HID_T) :: chsg_id, xssg_id, xs_id
	INTEGER :: error, numXS, i

		CALL getNumXSections(ich_id, numXS)
			CALL getGroup(ich_id, g_Solutions, chsg_id)

				CALL writeCHSolDataSets(chsg_id, chan_index, numXS)
!				CALL writeXSSolDataSets(chsg_id, chan_index, numXS)
			DO i = 1,numXS
				CALL getIndexGroup(i, ich_id, g_XSection, xs_id)
					CALL writeXSSolDataSets(xs_id, chan_index, i, numXS)
				CALL closeGroup(xs_id)

!! add xsect solution here???
			ENDDO

			CALL closeGroup(chsg_id)

END SUBROUTINE

SUBROUTINE writeXSSolDataSets(chsg_id, chan_index, xs_index, numXS)
	IMPLICIT NONE
	INTEGER(HID_T), INTENT(IN) :: chsg_id
	INTEGER, INTENT(IN) :: chan_index, xs_index, numXS
	INTEGER :: xss_id
		CALL getGroup(chsg_id, g_Solutions, xss_id)
		
			CALL writeXSSolTimeDataSet_REAL2(xss_id, chan_index, xs_index, numXS, gs_RTWSX, 3)
			CALL writeXSSolTimeDataSet_REAL2(xss_id, chan_index, xs_index, numXS, gs_RTWSY, 4)
			CALL writeXSSolTimeDataSet_REAL2(xss_id, chan_index, xs_index, numXS, gs_LTWSX, 5)
			CALL writeXSSolTimeDataSet_REAL2(xss_id, chan_index, xs_index, numXS, gs_LTWSY, 6)

		CALL closeGroup(xss_id)

END SUBROUTINE

SUBROUTINE writeCHSolDataSets(chsg_id, chan_index, numXS)
	IMPLICIT NONE
	INTEGER(HID_T), INTENT(IN) :: chsg_id
	INTEGER, INTENT(IN) :: chan_index, numXS
	INTEGER :: chsol_id
		CALL getGroup(chsg_id, g_ChanSol, chsol_id)
		
			CALL writeSolTimeDataSet_REAL(chsol_id, chan_index, numXS, gs_WSE, 1)
			CALL writeSolTimeDataSet_REAL(chsol_id, chan_index, numXS, gs_XSArea, 2)
			CALL writeSolTimeDataSet_REAL(chsol_id, chan_index, numXS, gs_HydRadius, 3)
			CALL writeSolTimeDataSet_REAL(chsol_id, chan_index, numXS, gs_HydDepth, 4)
			CALL writeSolTimeDataSet_REAL(chsol_id, chan_index, numXS, gs_TopWidth, 5)
			CALL writeSolTimeDataSet_REAL(chsol_id, chan_index, numXS, gs_Conveyance, 6)
			CALL writeSolTimeDataSet_REAL(chsol_id, chan_index, numXS, gs_Discharge, 7)
			CALL writeSolTimeDataSet_ShearStess_REAL(chsol_id, chan_index, numXS, gs_ShearStress)

		CALL closeGroup(chsol_id)

END SUBROUTINE

SUBROUTINE writeXSSolTimeDataSet_REAL2(chs_id, chan_index, xs_index, numXS, name, solIndex)
	IMPLICIT NONE
	INTEGER(HID_T), INTENT(IN) :: chs_id
	INTEGER, INTENT(IN) :: chan_index, xs_index, numXS, solIndex
	CHARACTER(LEN=*) :: name

	INTEGER(HID_T) :: cparms, memspace, dataspace, dset_id, filespace
	INTEGER :: RANK = 1
	INTEGER(HID_T), DIMENSION(2) :: chunkdims, offset, dimsw, hdims
	INTEGER(HID_T), DIMENSION(2) :: extdims, start, dcount
	INTEGER(HID_T), DIMENSION(1) :: dims, maxdims
	INTEGER :: error, i, icount
	INTEGER, DIMENSION(7) :: data_dims
	REAL :: value

!		dims = (/numXS, 1/)
!		maxdims = (/numXS,H5S_UNLIMITED_F/)
		dims = (/1/)
		maxdims = (/H5S_UNLIMITED_F/)
		dcount(1) = 1
		dcount(2) = 1
		CALL h5dopen_f(chs_id, name, dset_id, error)
		extdims(2) = numXS
		extdims(1) = tstep
		CALL h5dextend_f(dset_id, extdims, error)
		CALL h5dget_space_f(dset_id, filespace, error)
		icount = 0
!		DO i = nsce(chan_index), nsce(chan_index)-(numXS-1), -1
!		DO i = nsce(chan_index)-(numXS-1),nsce(chan_index)
			i = (nsce(chan_index) -numXS)+ xs_index
			CALL h5dget_space_f(dset_id, filespace, error)
			icount = icount+1
			start(2) = icount-1
			start(1) = tstep-1
			CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, start, dcount, error)
			CALL h5screate_simple_f(RANK, dcount, memspace, error)
			
			data_dims(1) = dcount(1)
			data_dims(2) = dcount(2)
			value = plnxsc(i, solIndex)
			CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, value, data_dims, error, &
							mem_space_id=memspace, file_space_id = filespace)

			CALL h5sclose_f(filespace, error)
			CALL h5sclose_f(memspace, error)

!		END DO
		
		CALL h5dclose_f(dset_id, error)

END SUBROUTINE

SUBROUTINE writeXSSolTimeDataSet_REAL(chs_id, chan_index, numXS, name, solIndex)
	IMPLICIT NONE
	INTEGER(HID_T), INTENT(IN) :: chs_id
	INTEGER, INTENT(IN) :: chan_index, numXS, solIndex
	CHARACTER(LEN=*) :: name

	INTEGER(HID_T) :: cparms, memspace, dataspace, dset_id, filespace
	INTEGER :: RANK = 2
	INTEGER(HID_T), DIMENSION(2) :: maxdims, dims, chunkdims, offset, dimsw, hdims
	INTEGER(HID_T), DIMENSION(2) :: extdims, start, dcount
	INTEGER :: error, i, icount
	INTEGER, DIMENSION(7) :: data_dims
	REAL :: value

!		dims = (/numXS, 1/)
!		maxdims = (/numXS,H5S_UNLIMITED_F/)
		dims = (/1, numXS/)
		maxdims = (/H5S_UNLIMITED_F, numXS/)
		dcount(1) = 1
		dcount(2) = 1
		CALL h5dopen_f(chs_id, name, dset_id, error)
		extdims(2) = numXS
		extdims(1) = tstep
		CALL h5dextend_f(dset_id, extdims, error)
		CALL h5dget_space_f(dset_id, filespace, error)
		icount = 0
!		DO i = nsce(chan_index), nsce(chan_index)-(numXS-1), -1
		DO i = nsce(chan_index)-(numXS-1),nsce(chan_index)
			CALL h5dget_space_f(dset_id, filespace, error)
			icount = icount+1
			start(2) = icount-1
			start(1) = tstep-1
			CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, start, dcount, error)
			CALL h5screate_simple_f(RANK, dcount, memspace, error)
			
			data_dims(1) = dcount(1)
			data_dims(2) = dcount(2)
			value = plnxsc(i, solIndex)
			CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, value, data_dims, error, &
							mem_space_id=memspace, file_space_id = filespace)

			CALL h5sclose_f(filespace, error)
			CALL h5sclose_f(memspace, error)

		END DO
		
		CALL h5dclose_f(dset_id, error)

END SUBROUTINE

SUBROUTINE writeSolTime(chs_id, name)
	IMPLICIT NONE
	INTEGER(HID_T), INTENT(IN) :: chs_id
	CHARACTER(LEN=*) :: name

	INTEGER(HID_T) :: cparms, memspace, dataspace, dset_id, filespace
	INTEGER :: RANK = 1
	INTEGER(HID_T), DIMENSION(1) :: maxdims, dims, chunkdims, offset, dimsw, hdims
	INTEGER(HID_T), DIMENSION(1) :: extdims, start, dcount
	INTEGER :: error, i, icount
	INTEGER, DIMENSION(7) :: data_dims
	REAL :: value

!		dims = (/numXS, 1/)
!		maxdims = (/numXS,H5S_UNLIMITED_F/)
		dims = (/1/)
		maxdims = (/H5S_UNLIMITED_F/)
		dcount(1) = 1
		CALL h5dopen_f(chs_id, name, dset_id, error)
		extdims(1) = tstep
		CALL h5dextend_f(dset_id, extdims, error)
		CALL h5dget_space_f(dset_id, filespace, error)
			CALL h5dget_space_f(dset_id, filespace, error)
			start(1) = tstep-1
			CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, start, dcount, error)
			CALL h5screate_simple_f(RANK, dcount, memspace, error)
			
			data_dims(1) = dcount(1)
			
			CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, ttime, data_dims, error, &
							mem_space_id=memspace, file_space_id = filespace)

			CALL h5sclose_f(filespace, error)
			CALL h5sclose_f(memspace, error)

		CALL h5dclose_f(dset_id, error)

END SUBROUTINE

SUBROUTINE writeSolTimeDataSet_REAL(chs_id, chan_index, numXS, name, solIndex)
	IMPLICIT NONE
	INTEGER(HID_T), INTENT(IN) :: chs_id
	INTEGER, INTENT(IN) :: chan_index, numXS, solIndex
	CHARACTER(LEN=*) :: name

	INTEGER(HID_T) :: cparms, memspace, dataspace, dset_id, filespace
	INTEGER :: RANK = 2
	INTEGER(HID_T), DIMENSION(2) :: maxdims, dims, chunkdims, offset, dimsw, hdims
	INTEGER(HID_T), DIMENSION(2) :: extdims, start, dcount
	INTEGER :: error, i, icount
	INTEGER, DIMENSION(7) :: data_dims
	REAL :: value

!		dims = (/numXS, 1/)
!		maxdims = (/numXS,H5S_UNLIMITED_F/)
		dims = (/1, numXS/)
		maxdims = (/H5S_UNLIMITED_F, numXS/)
		dcount(1) = 1
		dcount(2) = 1
		CALL h5dopen_f(chs_id, name, dset_id, error)
		extdims(2) = numXS
		extdims(1) = tstep
		CALL h5dextend_f(dset_id, extdims, error)
		CALL h5dget_space_f(dset_id, filespace, error)
		icount = 0
!		DO i = nsce(chan_index), nsce(chan_index)-(numXS-1), -1
		DO i = nsce(chan_index)-(numXS-1),nsce(chan_index)
			CALL h5dget_space_f(dset_id, filespace, error)
			icount = icount+1
			start(2) = icount-1
			start(1) = tstep-1
			CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, start, dcount, error)
			CALL h5screate_simple_f(RANK, dcount, memspace, error)
			
			data_dims(1) = dcount(1)
			data_dims(2) = dcount(2)
			value = curxsc(i, solIndex)
			CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, value, data_dims, error, &
							mem_space_id=memspace, file_space_id = filespace)

			CALL h5sclose_f(filespace, error)
			CALL h5sclose_f(memspace, error)

		END DO
		
		CALL h5dclose_f(dset_id, error)

END SUBROUTINE
	
SUBROUTINE writeSolTimeDataSet_ShearStess_REAL(chs_id, chan_index, numXS, name)
	IMPLICIT NONE
	INTEGER(HID_T), INTENT(IN) :: chs_id
	INTEGER, INTENT(IN) :: chan_index, numXS
	CHARACTER(LEN=*) :: name

	INTEGER(HID_T) :: cparms, memspace, dataspace, dset_id, filespace
	INTEGER :: RANK = 2
	INTEGER(HID_T), DIMENSION(2) :: maxdims, dims, chunkdims, offset, dimsw, hdims
	INTEGER(HID_T), DIMENSION(2) :: extdims, start, dcount
	INTEGER :: error, i, icount
	INTEGER, DIMENSION(7) :: data_dims
	REAL :: value, rise, run

!		dims = (/numXS, 1/)
!		maxdims = (/numXS,H5S_UNLIMITED_F/)
		dims = (/1, numXS/)
		maxdims = (/H5S_UNLIMITED_F, numXS/)
		dcount(1) = 1
		dcount(2) = 1
		CALL h5dopen_f(chs_id, name, dset_id, error)
		extdims(2) = numXS
		extdims(1) = tstep
		CALL h5dextend_f(dset_id, extdims, error)
		CALL h5dget_space_f(dset_id, filespace, error)
!		icount = nsce(chan_index)
		icount = numXS
		DO i = nsce(chan_index), nsce(chan_index)-(numXS-1), -1
			CALL h5dget_space_f(dset_id, filespace, error)
			start(2) = icount-1
			start(1) = tstep-1
			CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, start, dcount, error)
			CALL h5screate_simple_f(RANK, dcount, memspace, error)
			
			data_dims(1) = dcount(1)
			data_dims(2) = dcount(2)
			IF(icount > 1) THEN
				rise = (curxsc(i-1, 1) - curxsc(i, 1))
				run = sqrt(((plnxsc(i-1, 1)-plnxsc(i,1))**2. + (plnxsc(i-1,2)-plnxsc(i,2))**2.))
			ELSE
				rise = (curxsc(i, 1) - curxsc(i+1, 1))
				run = sqrt(((plnxsc(i, 1)-plnxsc(i+1,1))**2. + (plnxsc(i,2)-plnxsc(i+1,2))**2.))
			END IF
				value = 1000.*gg*curxsc(i, 3)*sin(atan(rise/run))

			CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, value, data_dims, error, &
							mem_space_id=memspace, file_space_id = filespace)

			CALL h5sclose_f(filespace, error)
			CALL h5sclose_f(memspace, error)
			icount = icount-1

		END DO
		
		CALL h5dclose_f(dset_id, error)

END SUBROUTINE

END MODULE
