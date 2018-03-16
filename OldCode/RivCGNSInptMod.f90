MODULE RivCGNSInptMod
USE RivVarMod
!USE CSedMod
USE CSedMod_DT_SUSP
USE TSplineMod
USE RivCGNSInptSubRtMod
USE RivRoughnessMod
USE NetPreisMain
IMPLICIT NONE
 	CONTAINS

 	SUBROUTINE read_CGNS(InputFile)
!		USE USER32
		IMPLICIT NONE
		CHARACTER(*), INTENT(IN) :: InputFile
		CHARACTER(250) :: ZONENAME, BASENAME, USERNAME

		INTEGER :: FID, BID, ZID, IER, iret
		INTEGER :: NUSER_DATA, USER_ITER, NZONES, ZONES_ITER
		INTEGER :: NBASES, BASES_ITER
		INTEGER :: CELLDIM, PHYSDIM
		INTEGER, DIMENSION(3,3) :: isize
		INTEGER, DIMENSION(3) :: irmin, irmax
		INTEGER, DIMENSION(6) :: irinddata
		CHARACTER(LEN = 250) :: name
		INTEGER ::namelen
!		namelen = len(InputFile)
!		name = TRIM(InputFile(1:namelen-3)//'cgns')
		CALL cg_open_f(InputFile, MODE_READ, FID, IER)
			IF(IER .NE. 0) THEN
                call cg_error_print_f()			
            ENDIF
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
						ns2 = isize(1, 1)
						nn = isize(2, 1)
					IF(PHYSDIM == 3) THEN
						nz = isize(3,1)
					ENDIF
!					CALL cg_goto_f(FID, BASES_ITER, IER, 'Zone_t', ZONES_ITER, 'FlowSolution_t',1, 'end')
!					CALL cg_rind_read_f(irinddata, IER)
!					nsext = irinddata(2)
				
!					CALL ALLOC_COMMON2D(ns2, nn)
!					CALL ALLOC_INIT2D(ns2, nn)
!					CALL ALLOC_ROUGHNESS(ns2, nn)

 					CALL cg_goto_f(FID, BASES_ITER, IER, 'Zone_t', ZONES_ITER, 'end')
					CALL cg_nuser_data_f(NUSER_DATA, IER)
					DO USER_ITER = 1, NUSER_DATA
						CALL cg_user_data_read_f(USER_ITER, USERNAME, IER)
						SELECT CASE(TRIM(USERNAME))
						CASE('FASTMECHParams')
!							Call CGNS_READ_FASTMECHPARAMS(FID, BASES_ITER, IER, ZONES_ITER, USER_ITER)
					    END SELECT
					ENDDO
					
					CALL ALLOC_COMMON2D(ns2, nsext, nn)
					CALL ALLOC_INIT2D(ns2, nn)
!					if (CALCCSED == 1) then
!	                    call alloc_csed()
	                    call alloc_csed_DT()
!	                endif

					DO USER_ITER = 1, NUSER_DATA
						CALL cg_user_data_read_f(USER_ITER, USERNAME, IER)
						SELECT CASE(TRIM(USERNAME))
!						CASE('FASTMECHParams')
!							Call CGNS_READ_FASTMECHPARAMS(FID, BASES_ITER, IER, ZONES_ITER, USER_ITER)
							!Called here to ensure nsext is defined
						CASE('GridInitialValues')
							CALL CGNS_READ_GRIDINITIALVALUES(FID, BASES_ITER, IER, ZONES_ITER, USER_ITER)
						CASE('GridInputValues')
							CALL CGNS_READ_GRIDINPUTVALUES(FID, BASES_ITER, IER, ZONES_ITER, USER_ITER)
						CASE('GridCenterlineValues')
							CALL CGNS_Read_CentLine(FID, BASES_ITER, IER, ZONES_ITER, USER_ITER)
						CASE('TimeSeriesandRatingCurve')
						    CALL CGNS_Read_TimeAndRatingValues(FID, BASES_ITER, IER, ZONES_ITER, USER_ITER)
!CALL NETPREIS2(ns, nn, topo, x, y, cd, q, stage, wse )
CALL NETPREIS2B(FID)
						END SELECT
					END DO
				END SELECT
			ENDDO	!ZONES
			CASE('MD_SWMS_3D')
				CALL cg_nzones_f(FID, BASES_ITER, NZONES, IER)
				DO ZONES_ITER = 1, NZONES
					CALL cg_zone_read_f(FID, BASES_ITER, ZONES_ITER, ZONENAME, isize, IER)

					SELECT CASE(TRIM(ZONENAME))
					CASE('CurvilinearOrthogonal')
						ns2 = isize(1, 1)
						nn = isize(2, 1)
						nz = isize(3,1)
						 CALL ALLOC_COMMON3D(ns2, nn, nz)
						 CALL alloc_csed3d_dt()
					END SELECT
				ENDDO

			END SELECT
		ENDDO !BASES
		CALL cg_close_f(FID, IER)
	END SUBROUTINE


END MODULE
