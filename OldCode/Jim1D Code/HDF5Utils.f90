MODULE HDF5UtilsMod
USE HDF5
USE Support
IMPLICIT NONE
	CHARACTER(LEN=9), PARAMETER :: g_1DNetwork = "1DNetwork"
	CHARACTER(LEN=8), PARAMETER :: g_Channels = "Channels"
	CHARACTER(LEN=8), PARAMETER :: g_Channel = "Channel_"
	CHARACTER(LEN=9), PARAMETER :: g_XSection = "XSection_"
	CHARACTER(LEN=19), PARAMETER :: g_BoundCond = 'Boundary Conditions'
	CHARACTER(LEN=11), PARAMETER :: g_USBC = 'Upstream'
	CHARACTER(LEN=13), PARAMETER :: g_DSBC = 'Downstream'
	CHARACTER(LEN=9), PARAMETER :: g_Juncts = 'Junctions'
	CHARACTER(LEN=9), PARAMETER :: g_Junct = 'Junction_'
	CHARACTER(LEN=10), PARAMETER :: g_Params = 'Parameters'
	CHARACTER(LEN=9), PARAMETER :: g_Solutions = 'Solutions'
	CHARACTER(LEN=17), PARAMETER :: g_ChanSol = 'Channel Solutions'
	CHARACTER(LEN=23), PARAMETER :: g_XSSol = 'Cross-Section Solutions'
	!    curxsc(ns,1)=ws elevation
	!    curxsc(ns,2)=xsec area
	!    curxsc(ns,3)=hydraulic radius =area/wet perim
	!    curxsc(ns,4)=hydraulic depth =area/top width
	!    curxsc(ns,5)=top width
	!    curxsc(ns,6)=conveyance
!	!    curxsc(ns,7)=discharge
!  xsloc(i,5)=emin  ! xsloc( ,5) set to min ch elev

	CHARACTER(LEN=23), PARAMETER :: gs_WSE = 'Water Surface Elevation'
	CHARACTER(LEN=19), PARAMETER :: gs_WSS = 'Water Surface Slope'
	CHARACTER(LEN=18), PARAMETER :: gs_XSArea = 'Cross-section Area'
	CHARACTER(LEN=16), PARAMETER :: gs_HydRadius = 'Hydraulic Radius'
	CHARACTER(LEN=15), PARAMETER :: gs_HydDepth = 'Hydraulic Depth'
	CHARACTER(LEN=9), PARAMETER :: gs_TopWidth = 'Top Width'
	CHARACTER(LEN=10), PARAMETER :: gs_Conveyance = 'Conveyance'
	CHARACTER(LEN=9), PARAMETER :: gs_Discharge = 'Discharge'
	 CHARACTER(LEN=16), PARAMETER :: gs_ShearStress = 'Bed Shear Stress'
	CHARACTER(LEN=5), PARAMETER :: gs_RTWSX = 'RTWSX'
	CHARACTER(LEN=5), PARAMETER :: gs_LTWSX = 'LTWSX'
	CHARACTER(LEN=5), PARAMETER :: gs_RTWSY = 'RTWSY'
	CHARACTER(LEN=5), PARAMETER :: gs_LTWSY = 'LTWSY'
	CHARACTER(LEN=4), PARAMETER :: gs_Time = 'Time'

CONTAINS
	SUBROUTINE createGroup(gi_id, name, go_id)
		IMPLICIT NONE
		INTEGER(HID_T), INTENT(IN) :: gi_id
		INTEGER(HID_T), INTENT(OUT) :: go_id
		CHARACTER(LEN=*), INTENT(IN) :: name
		INTEGER :: error
			CALL h5gcreate_f(gi_id, name, go_id, error)
	END SUBROUTINE

	SUBROUTINE getGroup(gi_id, name, go_id)
		IMPLICIT NONE
		INTEGER(HID_T), INTENT(IN) :: gi_id
		CHARACTER(LEN=*), INTENT(IN) :: name
		INTEGER(HID_T), INTENT(OUT) :: go_id
		INTEGER :: error

			CALL h5gopen_f(gi_id, name, go_id, error)
	END SUBROUTINE

	SUBROUTINE closeGroup(g_id)
		IMPLICIT NONE
		INTEGER(HID_T), INTENT(IN) :: g_id
		INTEGER :: error
			CALL h5gclose_f(g_id, error)
	END SUBROUTINE

	SUBROUTINE getIndexGroup(index, gi_id, name, go_id)
		IMPLICIT NONE
		INTEGER :: index
		INTEGER(HID_T), INTENT(IN) :: gi_id
		CHARACTER(LEN=*), INTENT(IN) :: name
		INTEGER(HID_T), INTENT(OUT) :: go_id
		INTEGER :: error
		CHARACTER(LEN=7) :: g_Index
		CHARACTER(LEN = 120) :: tmpstr

			WRITE(g_Index, '(I7)') index-1
			g_Index = ADJUSTL(g_Index)
			tmpstr = name//TRIM(g_Index)
			CALL h5gopen_f(gi_id, TRIM(tmpstr), go_id, error)
	END SUBROUTINE

	SUBROUTINE getNumChannels(g_ch_id, numCh)
		IMPLICIT NONE
		INTEGER(HID_T), INTENT(IN) :: g_ch_id
		INTEGER, INTENT(OUT) :: numCH
		INTEGER(HID_T) :: attid, dtype
		INTEGER, DIMENSION(7) :: data_dims
		INTEGER :: error
			CALL h5aopen_name_f(g_ch_id, "numChannels", attid, error)
			CALL h5aget_type_f(attid, dtype, error)
			data_dims(1) = 1
			CALL h5aread_f(attid, dtype, numCH, data_dims, error)
			CALL h5aclose_f(dtype, error)
			CALL h5aclose_f(attid, error)
	END SUBROUTINE

	SUBROUTINE getNumXSections(ch_id, numXS)
		IMPLICIT NONE
		INTEGER(HID_T), INTENT(IN) :: ch_id
		INTEGER, INTENT(OUT) :: numXS
		INTEGER(HID_T) :: attid, dtype
		INTEGER, DIMENSION(7) :: data_dims
		INTEGER :: error
			CALL h5aopen_name_f(ch_id, "numXSections", attid, error)
			CALL h5aget_type_f(attid, dtype, error)
			data_dims(1) = 1
			CALL h5aread_f(attid, dtype, numXS, data_dims, error)
			CALL h5aclose_f(dtype, error)
			CALL h5aclose_f(attid, error)
	END SUBROUTINE

END MODULE