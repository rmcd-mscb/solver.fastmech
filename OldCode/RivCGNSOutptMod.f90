MODULE RivCGNSOutptMod
USE RivVarMod
USE TSplineMod
USE RivCGNSOutptSubRtMod
USE RivVarTimeMod
IMPLICIT NONE
!INCLUDE "cgnslib_f.h"

	CONTAINS
	SUBROUTINE write_TimeIter_CGNS(InputFile)
!		USE USER32
!!DEC$ ATTRIBUTES REFERENCE, C, VARYING :: cg_goto_f
!!DEC$ ATTRIBUTES REFERENCE, C, VARYING :: cg_array_read_f
		IMPLICIT NONE
		CHARACTER(*), INTENT(IN) :: InputFile
		CHARACTER(LEN=250) :: ZONENAME, BASENAME, USERNAME, SOLNAME

		INTEGER :: FID, BID, ZID, IER, iret, STATUS
		INTEGER :: NUSER_DATA, USER_ITER, NZONES, ZONES_ITER
		INTEGER :: NBASES, BASES_ITER
		INTEGER :: NSOL, SOL_ITER
		INTEGER :: F2DSOL_ITER, F3DSOL_ITER, F1DSOL_ITER, F2DSOL_EXT_ITER
		INTEGER :: CELLDIM, PHYSDIM
		INTEGER, DIMENSION(3,3) :: isize
		INTEGER, DIMENSION(3) :: irmin, irmax
		CHARACTER(LEN = 250) :: name
		INTEGER ::namelen

		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar1, tmpvar2
		INTEGER, ALLOCATABLE, DIMENSION(:) :: tmpIntVar1
		INTEGER :: i, j, k, count
		REAL :: xx, yy
		REAL :: ux, uy, rcos, rsin

		INTEGER :: defOutput, newDefOutput

		INTEGER, DIMENSION(2) :: Dims

		INTEGER :: location
        FID = CGNSFILEID
! 		CALL cg_open_f(InputFile, MODE_MODIFY, FID, IER)
!			IF(IER .NE. 0) THEN
!!				iret = MESSAGEBOX(0, "CGNS ERROR"C, "Error"C, MB_OK)
!			ENDIF
		BID = 1
		CALL cg_nbases_f(FID, NBASES, IER)
		DO BASES_ITER = 1, NBASES
			CALL cg_base_read_f(FID, BASES_ITER, BASENAME, CELLDIM, PHYSDIM, IER)
			SELECT CASE(TRIM(BASENAME))
			CASE('MD_SWMS_1D')
    			CALL CGNS_BaseIter(FID, BASES_ITER)
				CALL cg_nzones_f(FID, BASES_ITER, NZONES, IER)
				DO ZONES_ITER = 1, NZONES
					CALL cg_zone_read_f(FID, BASES_ITER, ZONES_ITER, ZONENAME, isize, IER)
					SELECT CASE(TRIM(ZONENAME))
					    CASE('CurvilinearOrthogonal')
					        CALL CGNS_ZoneIter(FID, BASES_ITER, ZONES_ITER, 1)
					END SELECT
				END DO

			CASE('MD_SWMS_2D')
    			CALL CGNS_BaseIter(FID, BASES_ITER)
 				CALL cg_nzones_f(FID, BASES_ITER, NZONES, IER)
				DO ZONES_ITER = 1, NZONES
					CALL cg_zone_read_f(FID, BASES_ITER, ZONES_ITER, ZONENAME, isize, IER)
					SELECT CASE(TRIM(ZONENAME))
					CASE('CurvilinearOrthogonal_Ext')
					    if(ShowGridExt == 1) then
					        CALL CGNS_ZoneIter(FID, BASES_ITER, ZONES_ITER, 2)
					    endif
					CASE('CurvilinearOrthogonal')
                        if(ShowGridExt == 0) then
					        CALL CGNS_ZoneIter(FID, BASES_ITER, ZONES_ITER, 2)
					    endif
					END SELECT
 				END DO

			CASE('MD_SWMS_3D')
    			CALL CGNS_BaseIter(FID, BASES_ITER)
				CALL cg_nzones_f(FID, BASES_ITER, NZONES, IER)
				DO ZONES_ITER = 1, NZONES
					CALL cg_zone_read_f(FID, BASES_ITER, ZONES_ITER, ZONENAME, isize, IER)
					SELECT CASE(TRIM(ZONENAME))
					    CASE('CurvilinearOrthogonal')
					        CALL CGNS_ZoneIter(FID, BASES_ITER, ZONES_ITER, 3)
					END SELECT
				END DO
			END SELECT
		ENDDO
!		CALL cg_close_f(FID, IER)
			
	END SUBROUTINE
	
	SUBROUTINE write_TimeStep_CGNS(InputFile, timeStep, time)
!		USE USER32
!!DEC$ ATTRIBUTES REFERENCE, C, VARYING :: cg_goto_f
!!DEC$ ATTRIBUTES REFERENCE, C, VARYING :: cg_array_read_f
		IMPLICIT NONE
		CHARACTER(*), INTENT(IN) :: InputFile
		INTEGER, INTENT(IN) :: timeStep
		REAL, INTENT(IN) :: time
		CHARACTER(LEN=250) :: ZONENAME, BASENAME, USERNAME, SOLNAME

		INTEGER :: FID, BID, ZID, IER, iret, STATUS
		INTEGER :: NUSER_DATA, USER_ITER, NZONES, ZONES_ITER
		INTEGER :: NBASES, BASES_ITER
		INTEGER :: NSOL, SOL_ITER
		INTEGER :: F2DSOL_ITER, F3DSOL_ITER, F1DSOL_ITER, F2DSOL_EXT_ITER
		INTEGER :: CELLDIM, PHYSDIM
		INTEGER, DIMENSION(3,3) :: isize
		INTEGER, DIMENSION(3) :: irmin, irmax
		CHARACTER(LEN = 250) :: name
		INTEGER ::namelen

		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar1, tmpvar2
		INTEGER, ALLOCATABLE, DIMENSION(:) :: tmpIntVar1
		INTEGER :: i, j, k, count
		REAL :: xx, yy
		REAL :: ux, uy, rcos, rsin

		INTEGER :: defOutput, newDefOutput

		INTEGER, DIMENSION(2) :: Dims

		INTEGER :: location

        CHARACTER(LEN = 14) :: flowsol1D, flowsol2D, flowsol3D
        CHARACTER(LEN = 5) :: tsindex
        CHARACTER(LEN = 32) :: tmp_flowsol1D, tmp_flowsol2D, tmp_flowsol3D
        
        FID = CGNSFILEID
        
        TSPrintCount = TSPrintCount+1
        flowsol1D = '1DFlowSolution'
        flowsol2D = '2DFlowSolution'
        flowsol3D = '3DFlowSolution'
        
        WRITE(tsindex, '(I5)') TSPrintCount
        tsindex = ADJUSTL(tsindex)
        tmp_flowsol1D = flowsol1D//TRIM(tsindex)
        tmp_flowsol1D = TRIM(tmp_flowsol1D)

        tmp_flowsol2D = flowsol2D//TRIM(tsindex)
        tmp_flowsol2D = TRIM(tmp_flowsol2D)
        
        tmp_flowsol3D = flowsol3D//TRIM(tsindex)
        tmp_flowsol3D = TRIM(tmp_flowsol3D)

        SolNames(TSPrintCount) = tmp_flowsol2D
        SolNames1D(TSPrintCount) = tmp_flowsol1D
        SolNames3D(TSPrintCount) = tmp_flowsol3D
       
        TimeIncrements(TSPrintCount) = time
        DischIncrements(TSPrintCount) = q/1e6

!		namelen = len(InputFile)
!		name = TRIM(InputFile(1:namelen-3)//'cgns')
        if(timeStep.eq.0) then
		call write_error_CGNS(InputFile)
		endif

!		CALL cg_open_f(InputFile, MODE_MODIFY, FID, IER)
!			IF(IER .NE. 0) THEN
!!				iret = MESSAGEBOX(0, "CGNS ERROR"C, "Error"C, MB_OK)
!			ENDIF
		BID = 1
		CALL cg_nbases_f(FID, NBASES, IER)
		DO BASES_ITER = 1, NBASES
			CALL cg_base_read_f(FID, BASES_ITER, BASENAME, CELLDIM, PHYSDIM, IER)
			SELECT CASE(TRIM(BASENAME))
			CASE('MD_SWMS_1D')
				CALL cg_nzones_f(FID, BASES_ITER, NZONES, IER)
				DO ZONES_ITER = 1, NZONES
					CALL cg_zone_read_f(FID, BASES_ITER, ZONES_ITER, ZONENAME, isize, IER)
					SELECT CASE(TRIM(ZONENAME))
					CASE('CurvilinearOrthogonal')
!						CALL cg_nsols_f(FID,  BASES_ITER, ZONES_ITER, NSOL, IER)
!						DO SOL_ITER = 1, NSOL
!							CALL cg_sol_info_f(FID, BASES_ITER, ZONES_ITER, SOL_ITER, SOLNAME, location, IER)
!							SELECT CASE(SOLNAME)
!							CASE('1DFlowSolution')
								CALL cg_sol_write_f(FID, BASES_ITER, ZONES_ITER, tmp_flowsol1D, Vertex, F1DSOL_ITER, IER)

								CALL CGNS_Write_XSConvergence(FID, BASES_ITER, ZONES_ITER, F1DSOL_ITER)
								CALL CGNS_Write_AveWSElev(FID, BASES_ITER, ZONES_ITER, F1DSOL_ITER)
                                Call CGNS_Write_RadiusCurv(FID, BASES_ITER, ZONES_ITER, F1DSOL_ITER)

!							END SELECT
!						ENDDO

					END SELECT
				END DO

			CASE('MD_SWMS_2D')
!			CALL CGNS_BaseIter(FID, BASES_ITER)
				CALL CGNS_Write_Convergence(FID, BASES_ITER)
				IF(ShowGridExt == 1) THEN
!				    CALL CGNS_WRITE_2DGRIDEXTENSION(FID, BASES_ITER, ZONES_ITER)
				END IF
				CALL cg_nzones_f(FID, BASES_ITER, NZONES, IER)
				DO ZONES_ITER = 1, NZONES
					CALL cg_zone_read_f(FID, BASES_ITER, ZONES_ITER, ZONENAME, isize, IER)
					SELECT CASE(TRIM(ZONENAME))
					CASE('CurvilinearOrthogonal_Ext')
					    if(ShowGridExt == 1) then
					        if(timestep.le.0) then
            				    CALL CGNS_WRITE_2DGRIDEXTENSION(FID, BASES_ITER, ZONES_ITER)
            				    write(6,*)'writing grid extension', FID, BASES_ITER, ZONES_ITER
					        endif
					        CALL cg_sol_write_f(FID, BASES_ITER, ZONES_ITER, tmp_flowsol2D, Vertex, F2DSOL_EXT_ITER, IER)
            				    write(6,*)'writing solution extension', FID, BASES_ITER, ZONES_ITER, F2DSOL_EXT_ITER, IER
					            
!            					CALL CGNS_ZoneIter(FID, BASES_ITER, ZONES_ITER)
            					write(6,*) 'write ibc in grid extension', FID, BASES_ITER, ZONES_ITER, F2DSOL_EXT_ITER
								CALL CGNS_Write_GridIBC(FID, BASES_ITER, ZONES_ITER, F2DSOL_EXT_ITER)
								CALL CGNS_Write_WSElev(FID, BASES_ITER, ZONES_ITER, F2DSOL_EXT_ITER)
								CALL CGNS_Write_Velocity(FID, BASES_ITER, ZONES_ITER, F2DSOL_EXT_ITER)
								CALL CGNS_Write_Elevation(FID, BASES_ITER, ZONES_ITER, F2DSOL_EXT_ITER)
					    endif
					CASE('CurvilinearOrthogonal')
!						CALL cg_nsols_f(FID,  BASES_ITER, ZONES_ITER, NSOL, IER)
!						DO SOL_ITER = 1, NSOL
!							CALL cg_sol_info_f(FID, BASES_ITER, ZONES_ITER, SOL_ITER, SOLNAME, location, IER)
!							SELECT CASE(SOLNAME)
!							CASE('2DFlowSolution')
                                if(ShowGridExt == 0) then
								CALL cg_sol_write_f(FID, BASES_ITER, ZONES_ITER, tmp_flowsol2D, Vertex, F2DSOL_ITER, IER)

!					CALL CGNS_ZoneIter(FID, BASES_ITER, ZONES_ITER)
                                IF(IO_IBC) THEN
 								    CALL CGNS_Write_GridIBC(FID, BASES_ITER, ZONES_ITER, F2DSOL_ITER)
 								ENDIF
! 								IF(IO_VelXY) THEN
								    CALL CGNS_Write_Velocity(FID, BASES_ITER, ZONES_ITER, F2DSOL_ITER)
!								ENDIF
								IF(IO_InitVel) THEN
								    CALL CGNS_Write_InitVelocity(FID, BASES_ITER, ZONES_ITER, F2DSOL_ITER)
								ENDIF
!								IF(IO_ShearXY) THEN
								    CALL CGNS_Write_ShearStress(FID, BASES_ITER, ZONES_ITER, F2DSOL_ITER)
!								ENDIF
								IF(IO_WSE) THEN
								    CALL CGNS_Write_WSElev(FID, BASES_ITER, ZONES_ITER, F2DSOL_ITER)
								ENDIF
								IF(IO_Depth) THEN
								    CALL CGNS_Write_Depth(FID, BASES_ITER, ZONES_ITER, F2DSOL_ITER)
								ENDIF
								
								CALL CGNS_Write_FroudeNumber(FID, BASES_ITER, ZONES_ITER, F2DSOL_ITER)
!                                Call CGNS_Write_SuspSed(FID, BASES_ITER, ZONES_ITER, F2DSOL_ITER)
!                                
                                IF(IO_StressDiv) THEN
								    CALL CGNS_Write_StressDiv(FID, BASES_ITER, ZONES_ITER, F2DSOL_ITER)
								ENDIF
								IF(IO_CD) THEN
								    CALL CGNS_Write_CD(FID, BASES_ITER, ZONES_ITER, F2DSOL_ITER)
								ENDIF
								IF(IO_Elev) THEN
								    CALL CGNS_Write_Elevation(FID, BASES_ITER, ZONES_ITER, F2DSOL_ITER)
								ENDIF
!								    CALL CGNS_Write_ICON(FID, BASES_ITER, ZONES_ITER, F2DSOL_ITER)
								    CALL CGNS_Write_HAREA(FID, BASES_ITER, ZONES_ITER, F2DSOL_ITER)
				!				CALL CGNS_Write_XSConvergence(FID, BASES_ITER, ZONES_ITER, F2DSOL_ITER)
				!				CALL CGNS_Write_AveWSElev(FID, BASES_ITER, ZONES_ITER, F2DSOL_ITER)
				                endif

								IF(CALCQUASI3D == 1) THEN
								    IF(IO_Helix) THEN
									CALL CGNS_Write_HelixStrength(FID, BASES_ITER, ZONES_ITER, F2DSOL_ITER)
									ENDIF
									if(CALCCSED == 1) THEN
									    CALL CGNS_Write_SedFluxComp(FID, BASES_ITER, ZONES_ITER, F2DSOL_ITER)
									    IF(TRANSEQTYPE == 2) THEN
									        CALL CGNS_Write_SandDepth(FID, BASES_ITER, ZONES_ITER, F2DSOL_ITER)
									        CALL CGNS_Write_SandFraction(FID, BASES_ITER, ZONES_ITER, F2DSOL_ITER)
!									        CALL CGNS_Write_LSub(FID, BASES_ITER, ZONES_ITER, F2DSOL_ITER)
									    ENDIF
									ENDIF
								ENDIF
!							END SELECT
!						ENDDO

!						CALL cg_goto_f(FID, BASES_ITER, IER, 'Zone_t', ZONES_ITER, 'end')
!						CALL cg_nuser_data_f(NUSER_DATA, IER)
!						DO USER_ITER = 1, NUSER_DATA
!							CALL cg_user_data_read_f(USER_ITER, USERNAME, IER)
!							SELECT CASE(TRIM(USERNAME))
!								CASE('GridInputValues')
!								CALL CGNS_Write_GridIBC(FID, BASES_ITER, ZONES_ITER, USER_ITER)
!							END SELECT
!						ENDDO
					END SELECT
				END DO

			CASE('MD_SWMS_3D')
				CALL cg_nzones_f(FID, BASES_ITER, NZONES, IER)
				DO ZONES_ITER = 1, NZONES
					CALL cg_zone_read_f(FID, BASES_ITER, ZONES_ITER, ZONENAME, isize, IER)
					SELECT CASE(TRIM(ZONENAME))
					CASE('CurvilinearOrthogonal')
					    IF (IO_3DOutput) THEN
!					    if(timestep.le.0) then
						    CALL CGNS_Write_CoordinateZ(FID, BASES_ITER, ZONES_ITER)
!            				write(6,*)'writing grid CoordinateZ', FID, BASES_ITER, ZONES_ITER
!					    endif

!						CALL cg_nsols_f(FID,  BASES_ITER, ZONES_ITER, NSOL, IER)
!						IF(NSOL == 0) THEN
							CALL cg_sol_write_f(FID, BASES_ITER, ZONES_ITER, tmp_flowsol3D, Vertex, F3DSOL_ITER, IER)
							CALL CGNS_Write_Velocity3D(FID, BASES_ITER, ZONES_ITER, F3DSOL_ITER)
							CALL CGNS_Write_GridIBC3D(FID, BASES_ITER, ZONES_ITER, F3DSOL_ITER)
							CALL CGNS_Write_GridZ(FID, BASES_ITER, ZONES_ITER, F3DSOL_ITER)
!							CALL cg_sol_write_f(FID, BASES_ITER, ZONES_ITER, 'IBC', Vertex, F3DSOL_ITER, IER)
!							CALL CGNS_Write_GridIBC(FID, BASES_ITER, ZONES_ITER, F3DSOL_ITER)

!						ELSE
!							DO SOL_ITER = 1, NSOL
!								CALL cg_sol_info_f(FID, BASES_ITER, ZONES_ITER, SOL_ITER, SOLNAME, location, IER)
!								SELECT CASE(SOLNAME)
!								CASE('3DFlowSolution')
!									CALL CGNS_Write_Velocity3D(FID, BASES_ITER, ZONES_ITER, F3DSOL_ITER)
!!								CASE('IBC')
!!									CALL CGNS_Write_GridIBC(FID, BASES_ITER, ZONES_ITER, F3DSOL_ITER)
!
!								END SELECT
!							ENDDO
!						ENDIF
                        ENDIF
					END SELECT
				END DO


			END SELECT
		ENDDO
!		CALL cg_close_f(FID, IER)
			
	END SUBROUTINE


	SUBROUTINE write_error_CGNS(OutputFile)
!		USE USER32
!DEC$ ATTRIBUTES REFERENCE, C, VARYING :: cg_goto_f
!DEC$ ATTRIBUTES REFERENCE, C, VARYING :: cg_array_write_f
		IMPLICIT NONE
		CHARACTER(*), INTENT(IN) :: OutputFile
		CHARACTER(250) :: ZONENAME, BASENAME, USERNAME

		INTEGER :: FID, BID, ZID, IER, iret
		INTEGER :: NUSER_DATA, USER_ITER, NZONES, ZONES_ITER
		INTEGER :: NBASES, BASES_ITER
		INTEGER :: CELLDIM, PHYSDIM
		INTEGER, DIMENSION(3,3) :: isize
		INTEGER, DIMENSION(3) :: irmin, irmax
		INTEGER, DIMENSION(6) :: irinddata
		CHARACTER(LEN = 250) :: name, errorstring
		INTEGER ::namelen
		INTEGER, DIMENSION(1) :: dimvals
!		namelen = len(InputFile)
!		name = TRIM(InputFile(1:namelen-3)//'cgns')
        FID = CGNSFILEID
!		CALL cg_open_f(OutputFile, MODE_MODIFY, FID, IER)
!			IF(IER .NE. 0) THEN
!!				iret = MESSAGEBOX(0, "CGNS ERROR"C, "Error"C, MB_OK)
!			ENDIF
		BID = 1
		CALL cg_nbases_f(FID, NBASES, IER)
		DO BASES_ITER = 1, NBASES
			CALL cg_base_read_f(FID, BASES_ITER, BASENAME, CELLDIM, PHYSDIM, IER)
			
			SELECT CASE(TRIM(BASENAME))
			CASE('MD_SWMS_1D')

			CASE('MD_SWMS_2D')
				CALL cg_nzones_f(FID, BASES_ITER, NZONES, IER)
				DO ZONES_ITER = 1, NZONES
					CALL cg_zone_read_f(FID, BASES_ITER, ZONES_ITER, ZONENAME, isize, IER)

					SELECT CASE(TRIM(ZONENAME))
					CASE('CurvilinearOrthogonal')
						CALL cg_goto_f(FID, BASES_ITER, IER, 'Zone_t', ZONES_ITER, 'end')
						CALL cg_nuser_data_f(NUSER_DATA, IER)
						DO USER_ITER = 1, NUSER_DATA
							CALL cg_user_data_read_f(USER_ITER, USERNAME, IER)

							SELECT CASE(TRIM(USERNAME))
								CASE('Error')
				call cg_goto_f(FID, BASES_ITER, IER, 'Zone_t', ZONES_ITER,'UserDefinedData_t', USER_ITER, 'end')
				dimvals(1) = 1
				CALL cg_array_write_f('ErrorCode', Integer, 1, 1, errorcode, ier)
									SELECT CASE (errorcode)
									CASE (-1) 
				Write(errorstring, '(A70,I5)') 'Error: Parameters chosen do not result in convergence at iteration:', iter	
				dimvals(1) = len(errorstring)
				CALL cg_array_write_f('ErrorString', Character, 1, dimvals, errorstring, ier)
									CASE (-10)
		!								status = NF_PUT_ATT_TEXT(NCID, NF_GLOBAL, "r2d_ErrorText",85, 
		!						 +"At least one row of nodes are completely dry - check water surface 
		!						 + boundary condition")
				Write(errorstring, '(A100)') 'At least one row of nodes are completely dry - check water surface boundary condition'	
				dimvals(1) = len(errorstring)
				CALL cg_array_write_f('ErrorString', Character, 1, dimvals, errorstring, ier)

									CASE DEFAULT
		!								status = NF_PUT_ATT_TEXT(NCID, NF_GLOBAL, "r2d_ErrorText",
		!						 +		23, "An unkown error occured")
				Write(errorstring, '(A70)') 'An unkown error occured'	
				dimvals(1) = len(errorstring)
				CALL cg_array_write_f('ErrorString', Character, 1, dimvals, errorstring, ier)
									END SELECT
!                                CALL cg_close_f(FID, IER)
                                return
							END SELECT

						ENDDO
					END SELECT
				ENDDO
			END SELECT
		ENDDO

!		CALL cg_close_f(FID, IER)
			
	END SUBROUTINE
	
	SUBROUTINE OPEN_CGNS(OutputFile)
		IMPLICIT NONE
		CHARACTER(*), INTENT(IN) :: OutputFile

		INTEGER :: IER

		CALL cg_open_f(OutputFile, MODE_MODIFY, CGNSFILEID, IER)
			IF(IER .NE. 0) THEN
!				iret = MESSAGEBOX(0, "CGNS ERROR"C, "Error"C, MB_OK)
			ENDIF
    END SUBROUTINE
    
    SUBROUTINE CLOSE_CGNS()
		IMPLICIT NONE
		INTEGER :: IER
		
		CALL cg_close_f(CGNSFILEID, IER)

    END SUBROUTINE
END MODULE