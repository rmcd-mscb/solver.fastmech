MODULE RivCGNSInptSubRtMod
USE RivVarMod
USE CalcCond
USE RivVarVertMod
USE RivVarTimeMod
USE RivRoughnessMod
!USE CSedMod
USE CSedMod_DT_SUSP
USE TSplineMod
IMPLICIT NONE
!DEC$ ATTRIBUTES REFERENCE, C, VARYING :: CG_GOTO_F
!DEC$ ATTRIBUTES REFERENCE, C, VARYING :: CG_ARRAY_READ_F

INCLUDE "cgnslib_f.h"
!INCLUDE "cgnswin_f.h"
	CONTAINS

!		SUBROUTINE CGNS_Read_FASTMECHParams(FID, BID, IER, ZONES_ITER, USER_ITER)
!		IMPLICIT NONE
!		INTEGER, INTENT(IN) :: FID, BID, IER
!		INTEGER, INTENT(IN) :: USER_ITER, ZONES_ITER
!		INTEGER :: NARRAYS, IARRAY
!		CHARACTER(LEN = 250) :: name
!		INTEGER :: datatype, nndim
!		INTEGER, DIMENSION(3) :: dim_vals
!		INTEGER :: status, i, j, count
!		INTEGER :: tmpcdtype, tmpnsext, tcalcwet
!		REAL :: tmpcd, tmpevc, tmpnsextslope
!		
!		!! For DT Wilcock-Kenworthy Transport !! rmcd 5/2/07
!!		REAL :: tmphmax, tmpabh, tmpbbh, tmplsub, tmplsubactdepth
!
!		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar, tmpvar2
!		INTEGER, ALLOCATABLE, DIMENSION(:) :: tmpvari
!		REAL, ALLOCATABLE, DIMENSION(:) :: tmpcdregion
!
!
!			call cg_goto_f(FID, BID, IER, 'Zone_t', ZONES_ITER,'UserDefinedData_t', USER_ITER, 'end')
!			call cg_narrays_f(NARRAYS, IER)
!			DO IARRAY = 1, NARRAYS
!				call cg_array_info_f(iarray, name, datatype, nndim, dim_vals, ier)
!					SELECT CASE(TRIM(name))
!!					CASE('FM_SolAttRlx')
!!						ALLOCATE(tmpvar(dim_vals(1)), STAT = status)
!!						call cg_array_read_f(IARRAY, tmpvar, IER)
!!						erelax = tmpvar(1);
!!						urelax = tmpvar(2);
!!						arelax = tmpvar(3);
!!						DEALLOCATE(tmpvar, STAT = status)
!!					CASE('FM_SolAttItm')
!!						call cg_array_read_f(IARRAY, itm, IER)
!!	                CASE('FM_SolAttNTimeStep')
!!						call cg_array_read_f(IARRAY, nsteps, IER)
!!	                CASE('FM_SolAttDT')
!!						call cg_array_read_f(IARRAY, dt, IER)
!!	                CASE('FM_SolAttInterItm')
!!						call cg_array_read_f(IARRAY, interitm, IER)
!!	                CASE('FM_SolAttPlinc')
!!						call cg_array_read_f(IARRAY, iplinc, IER)
!!					CASE('FM_SolAttVarDischSTime')
!!					    call cg_array_read_f(IARRAY, VarDischStartTime, IER)
!!					CASE('FM_SolAttVarDischETime')
!!					    call cg_array_read_f(IARRAY, VarDischEndTime, IER)
!!					CASE('FM_SolAttDbgStop')
!!					    call cg_array_read_f(IARRAY, DEBUGSTOP, IER)
!!					CASE('FM_SolAttIterOut')
!!					    call cg_array_read_f(IARRAY, IterationOut, IER)
!!					CASE('FM_SolAttIncPlOut')
!!					    call cg_array_read_f(IARRAY, IterPlOut, IER)
!!					CASE('FM_SolAttDbgTimeStep')
!!					    call cg_array_read_f(IARRAY, DbgTimeStep, IER)
!!					CASE('FM_SolAttDbgIterNum')
!!					    call cg_array_read_f(IARRAY, DbgIterNum, IER)
!!					CASE('FM_SolAttMaxInterIter')
!!					    call cg_array_read_f(IARRAY, maxInterIterMult, IER)
!					    
!!					CASE('FM_SolWrtElev')
!!					    call cg_array_read_f(IARRAY, IO_Elev, IER)
!!					    
!!					CASE('FM_SolWrtIBC')
!!					    call cg_array_read_f(IARRAY, IO_IBC, IER)
!!					
!!					CASE('FM_SolWrtWSE')
!!					    call cg_array_read_f(IARRAY, IO_WSE, IER)
!!					
!!					CASE('FM_SolWrtVelXY')
!!					    call cg_array_read_f(IARRAY, IO_VelXY, IER)
!!					
!!					CASE('FM_SolWrtVelSN')
!!					    call cg_array_read_f(IARRAY, IO_VelSN, IER)
!!					
!!					CASE('FM_SolWrtShearXY')
!!					    call cg_array_read_f(IARRAY, IO_ShearXY, IER)
!!					
!!					CASE('FM_SolWrtShearSN')
!!					    call cg_array_read_f(IARRAY, IO_ShearSN, IER)
!!					
!!					CASE('FM_SolWrtShearDiv')
!!					    call cg_array_read_f(IARRAY, IO_StressDiv, IER)
!!					
!!					CASE('FM_SolWrtInitVel')
!!					    call cg_array_read_f(IARRAY, IO_InitVel, IER)
!!					
!!					CASE('FM_SolWrtDepth')
!!					    call cg_array_read_f(IARRAY, IO_Depth, IER)
!!					
!!					CASE('FM_SolWrtCD')
!!					    call cg_array_read_f(IARRAY, IO_CD, IER)
!!					    
!!					CASE('FM_SolWrtHelix')
!!					    call cg_array_read_f(IARRAY, IO_Helix, IER)
!!					
!!					CASE('FM_SolWrt3DOutput')
!!					    call cg_array_read_f(IARRAY, IO_3DOutput, IER)
!!					
!!					CASE('FM_SolWrtKEG')
!!					    call cg_array_read_f(IARRAY, IO_KEG, IER)
!!					
!!					CASE('FM_SolWrtVort')
!!					    call cg_array_read_f(IARRAY, IO_Vort, IER)
!!					
!!					CASE('FM_SolWrtVelStrain')
!!					    call cg_array_read_f(IARRAY, IO_VelStrain, IER)
!
!!					CASE('FM_HydAttQ')
!!						call cg_array_read_f(IARRAY, q, IER)
!!						q = q*1.e6 
!!					CASE('FM_HydAttVarDischType')
!!					    call cg_array_read_f(IARRAY, VarDischType, IER)
!					CASE('FM_HydAttVarDischTSIndex')
!					    call cg_array_read_f(IARRAY, DischTSNum, IER)
!!					CASE('FM_HydAttVelDepthCoef')
!!						call cg_array_read_f(IARRAY, vdc, IER)
!!					CASE('FM_HydAttVelAngleCoef')
!!						call cg_array_read_f(IARRAY, vac, IER)
!!					CASE('FM_HydAttVelBC')
!!						call cg_array_read_f(IARRAY, vbc, IER)
!!					CASE('FM_VelBCDistance')
!!						ALLOCATE(vbcdist(dim_vals(1)), STAT = status)
!!						call cg_array_read_f(IARRAY, vbcdist, IER)
!!					CASE('FM_VelBCVelocity')
!!						ALLOCATE(vbcvel(dim_vals(1)), STAT = status)
!!						call cg_array_read_f(IARRAY, vbcvel, IER)
!!					CASE('FM_VelBCAngle')
!!						ALLOCATE(vbcang(dim_vals(1)), STAT = status)
!!						call cg_array_read_f(IARRAY, vbcang, IER)
!!					CASE('FM_HydAttSType')
!!						call cg_array_read_f(IARRAY, stagetype, IER)
!!					CASE('FM_HydAttRBSWE')
!!						call cg_array_read_f(IARRAY, rbwse, IER)
!!						rbwse = rbwse*100.
!!					CASE('FM_HydAttLBWSE')
!!						call cg_array_read_f(IARRAY, lbwse, IER)
!!						lbwse = lbwse*100.
!!					CASE('FM_HydAttVarStgType')
!!					    call cg_array_read_f(IARRAY, VarStageType, IER)
!!					CASE('FM_HydAttVarStgTSIndex')
!!					    call cg_array_read_f(IARRAY, StageTSNum, IER)
!!					CASE('FM_HydAttVarStgRCIndex')
!!					    call cg_array_read_f(IARRAY, StageRCNum, IER)
!!					CASE('FM_HydAttWS')
!!						call cg_array_read_f(IARRAY, wselev, IER)
!!						wselev = wselev*100.
!!					CASE('FM_HydAttWSType')
!!						call cg_array_read_f(IARRAY, wstype, IER)
!!            CASE('FM_HydAttHotStart')
!!                 CALL cg_array_read_f(iarray,HotStart,IER)
!!            
!!            CASE('FM_HydAttHSFile')
!!                 CALL cg_array_read_f(iarray,CGNSHSFile,IER)
!!            
!!            CASE('FM_HydAttSolIndex')
!!                 CALL cg_array_read_f(iarray,SolnIndex,IER)
!
!!					CASE('FM_HydAttCD')
!!						call cg_array_read_f(IARRAY, constcd, IER)
!!					CASE('FM_HydAttRoughnessType')
!!						call cg_array_read_f(IARRAY, roughnesstype, IER)
!!					CASE('FM_HydAttCDType')
!!						call cg_array_read_f(IARRAY, cdtype, IER)
!!							IF (cdtype.ne.0) THEN		
!!								DO i = 1, ns2
!!									DO j = 1, nn
!!										znaught2(i,j) = tmpcd*100.
!!										cd2(i,j) = tmpcd 
!!									END DO
!!								END DO
!!							ENDIF
!!                    CASE('FM_HydAttCDMin')
!!                        call cg_array_read_f(IARRAY, cdmin, IER)
!!                    CASE('FM_HydAttCDMax')
!!                        call cg_array_read_f(IARRAY, cdmax, IER)
!!					CASE('FM_HydAttEVC')
!!						call cg_array_read_f(IARRAY, tmpevc, IER)
!!						evc = tmpevc*100.*100.
!!					CASE('FM_SolAttNSExt')
!!						call cg_array_read_f(IARRAY, tmpnsext, IER)
!!						nsext = tmpnsext
!!					CASE('FM_SolAttNSExtSlope')
!!						call cg_array_read_f(IARRAY, tmpnsextslope, IER)
!!						nsextslope = tmpnsextslope
!!					CASE('FM_SolAttNSExtShow')
!!					    call cg_array_read_f(IARRAY, ShowGridExt, IER)
!!					CASE('FM_HydAttVelBCDS')
!!					    call cg_array_read_f(IARRAY, vbcds, IER)
!!					CASE('FM_HydAttVelDepthCoefDS')
!!					    call cg_array_read_f(IARRAY, vdcds, IER)
!!					CASE('FM_HydAttVelAngleCoefDS')
!!					    call cg_array_read_f(IARRAY, vacds, IER)
!!					CASE('FM_HydAttDryMinDepth')
!!						call cg_array_read_f(IARRAY, hmin, IER)
!!						hmin = hmin*100.
!!				    CASE('FM_HydAttDryType')
!!				        call cg_array_read_f(IARRAY, dryType, IER)
!!					CASE('FM_HydAttWetMinDepth')
!!						call cg_array_read_f(IARRAY, hwmin, IER)
!!						hwmin = hwmin*100.
!!					CASE('FM_HydAttWetIterInterval')
!!						call cg_array_read_f(IARRAY, hiterInterval, IER)
!!					CASE('FM_HydAttWetIterStop')
!!						call cg_array_read_f(IARRAY, hiterstop, IER)
!!					CASE('FM_HydAttCalcWetting')
!!						call cg_array_read_f(IARRAY, tcalcwet, IER)
!!						IF(tcalcwet == 1) THEN
!!							hcalcwetting = .TRUE.
!!						ELSE
!!							hcalcwetting = .FALSE.
!!						ENDIF
!!					CASE('FM_HydAttReWetDepth')
!!						call cg_array_read_f(IARRAY, hrwdepth, IER)
!!						hrwdepth = hrwdepth*100.
!
!!					CASE('FM_HydAttLEVType')
!!						call cg_array_read_f(IARRAY, LEVType, IER)
!!					CASE('FM_HydAttLEVChangeIter')
!!						call cg_array_read_f(IARRAY, LEVChangeIter, IER)
!!
!!					CASE('FM_LEVStartIter')
!!						call cg_array_read_f(IARRAY, LEVBegIter, IER)
!!
!!					CASE('FM_LEVEndIter')
!!						call cg_array_read_f(IARRAY, LEVEndIter, IER)
!!
!!					CASE('FM_StartLEV')
!!						call cg_array_read_f(IARRAY, startLEV, IER)
!!						startLEV = startLEV*100.*100.
!!					CASE('FM_EndLEV')
!!						call cg_array_read_f(IARRAY, endLEV, IER)
!!						endLEV = endLEV*100.*100.
!
!					CASE('FM_Z0')
!						call cg_array_read_f(IARRAY, zo, IER)
!!					CASE('FM_Quasi3D')
!!					    call cg_array_read_f(IARRAY, CALCQUASI3D, IER)
!!					CASE('FM_Quasi3D_CalcRS')
!!					    call cg_array_read_f(IARRAY, CALCQUASI3DRS, IER)
!!
!!					CASE('FM_Quasi3D_MinRS')
!!					    call cg_array_read_f(IARRAY, MinRS, IER)
!!					    MinRs = MinRs * 100.
!					    
!!					CASE('FM_SedTrans')
!!						call cg_array_read_f(IARRAY, CALCCSED, IER)
!!					CASE('FM_SedDuneH')
!!						call cg_array_read_f(IARRAY, HD, IER)
!!					CASE('FM_SedDuneWL')
!!						call cg_array_read_f(IARRAY, WD, IER)
!! 					CASE('FM_SedSmoothLvl')
!!						call cg_array_read_f(IARRAY, SEDSMOO, IER)
!!				    CASE('FM_SedSmoothWeight')
!!				        call cg_array_read_f(IARRAY, SEDSMOOWGHT, IER)
!!
!!					CASE('FM_SedTransEq')
!!						call cg_array_read_f(IARRAY, TRANSEQTYPE, IER)
!!					CASE('FM_SedGrainSize')
!!						call cg_array_read_f(IARRAY, DIN, IER)
!!						din = din*100.
!!					CASE('FM_SedBCNode')
!!					    call cg_array_read_f(IARRAY,SEDBCNODE , IER)
!!					CASE('FM_SedCalcGravCorr')
!!					    call cg_array_read_f(IARRAY, CALCGRAVCORR, IER)
!!					CASE('FM_SedGravCorrType')
!!					    call cg_array_read_f(IARRAY, GRAVCORRTYPE, IER)
!!					CASE('FM_SedFlatBedCorrCoef')
!!					    call cg_array_read_f(IARRAY, GRAVFLATBEDCORRCOEF, IER)
!!					CASE('FM_SedSAngleOfRepose')
!!					    call cg_array_read_f(IARRAY, SUBANGLEOFREPOSE, IER)
!!					CASE('FM_SedTSFracDepth')
!!					    call cg_array_read_f(IARRAY, TSFracDepth, IER)
!!					CASE('FM_SedTransAuto')
!!					    call cg_array_read_f(IARRAY, CALCSEDAUTO, IER)
!!					CASE('FM_SedBCFraction')
!!					    call cg_array_read_f(IARRAY, BCFRACTION, IER)
!
!!					!! variables for DT Wilcock-Kenworthy transport function !! rmcd 5/2/07
!!					CASE('FM_Sed_WK_RG0')
!!					    call cg_array_read_F(IARRAY, taustar_rg0, IER)
!!					CASE('FM_Sed_WK_RG1')
!!					    call cg_array_read_F(IARRAY, taustar_rg1, IER)
!!					CASE('FM_Sed_WK_RS1')
!!					    call cg_array_read_F(IARRAY, taustar_rs1, IER)
!!					CASE('FM_Sed_WK_Alpha')
!!					    call cg_array_read_F(IARRAY, alpha, IER)
!!					CASE('FM_Sed_WK_AK')
!!					    call cg_array_read_F(IARRAY, AK, IER)
!!					CASE('FM_Sed_WK_Chi')
!!					    call cg_array_read_F(IARRAY, Chi, IER)
!!					CASE('FM_Sed_WK_PhiPrime')
!!					    call cg_array_read_F(IARRAY, phi_prime, IER)
!!					CASE('FM_Sed_WK_HMax')
!!					    call cg_array_read_F(IARRAY, tmphmax, IER)
!!					    tmphmax = tmphmax * 100.
!!					CASE('FM_Sed_WK_ABH')
!!					    call cg_array_read_F(IARRAY, tmpabh, IER)
!!					CASE('FM_Sed_WK_BBH')
!!					    call cg_array_read_F(IARRAY, tmpbbh, IER)
!!					CASE('FM_Sed_WK_LSub')
!!					    call cg_array_read_F(IARRAY, tmplsub, IER)
!!					    tmplsub = tmplsub*100.
!!					CASE('FM_Sed_WK_LSubActDepth')
!!					    call cg_array_read_F(IARRAY, tmplsubactdepth, IER)
!!					    tmplsubactdepth = tmplsubactdepth*100.
!!					    
!!					CASE('FM_Sed_WK_DSand')
!!					    call cg_array_read_F(IARRAY, Dsand, IER)
!!					    Dsand = Dsand * 100.
!!					CASE('FM_Sed_WK_DGravel')
!!					    call cg_array_read_F(IARRAY, Dg, IER)
!!					    Dg = Dg*100.
!!					CASE('FM_Sed_WK_Z0Min')
!!					    call cg_array_read_F(IARRAY, Z0Min, IER)
!!					        Z0Min = Z0Min*100.
!!					CASE('FM_Sed_WK_Z0Max')
!!					    call cg_array_read_F(IARRAY, Z0Max, IER)
!!					    Z0Max = Z0Max*100.		
!!					CASE('FM_Sed_WK_RoughType')
!!					    call cg_array_read_F(IARRAY, WK_RoughType, IER)			    
!!					CASE('FM_Sed_WK_RoughShapeParam')
!!					    call cg_array_read_F(IARRAY, WK_RoughShapeParam, IER)		
!!							    
!!							    
!!					CASE('FM_Sed_WK_FsMin')
!!					    call cg_array_read_F(IARRAY, Fsmin, IER)
!!					    
!!					CASE('FM_Sed_WK_ConstBndry')
!!					    call cg_array_read_F(IARRAY, WKConstBndry, IER)
!					    
!					
!
!				END SELECT
!			END DO
!			CALL cg_goto_f(FID, BID, IER, 'Zone_t', ZONES_ITER, 'end')
!
!		END SUBROUTINE

		SUBROUTINE CGNS_Read_GridInitialValues(FID, BID, IER, ZONES_ITER, USER_ITER)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BID, IER
		INTEGER, INTENT(IN) :: USER_ITER, ZONES_ITER
		INTEGER :: NARRAYS, IARRAY
		CHARACTER(LEN = 250) :: name
		INTEGER :: datatype, nndim
		INTEGER, DIMENSION(3) :: dim_vals
		INTEGER :: status, i, j, count
		INTEGER :: countij, countji

		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar, tmpvar2, vx, vy
		INTEGER, ALLOCATABLE, DIMENSION(:) :: tmpvari
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpcdregion


			call cg_goto_f(FID, BID, IER, 'Zone_t', ZONES_ITER,'UserDefinedData_t', USER_ITER, 'end')
			call cg_narrays_f(NARRAYS, IER)
			DO IARRAY = 1, NARRAYS
				call cg_array_info_f(iarray, name, datatype, nndim, dim_vals, ier)
					SELECT CASE(TRIM(name))
					CASE('WaterSurfaceElevation')
						ALLOCATE(tmpvar(dim_vals(1)), STAT = status)
						call cg_array_read_f(IARRAY, tmpvar, IER)
						DO i = 1, ns2
							DO j = 1, nn
								countij = ((i-1)*nn)+j
								countji = ((j-1)*ns2)+i
								IF(j == (nn/2)+1) THEN
									hav2(i) = tmpvar(countji)*100.
									!convert to cm from m
								ENDIF
!								hav2n(i,j) = tmpvar(countji)*100.
							ENDDO
						ENDDO
						DEALLOCATE(tmpvar, STAT = status)

					CASE('IBC')
						ALLOCATE(tmpvari(dim_vals(1)), STAT = status)
						call cg_array_read_f(IARRAY, tmpvari, IER)
						DO i = 1, ns2
							DO j = 1, nn
								countji = ((j-1)*ns2)+i
								ibc2(i,j) = tmpvari(countji)
							ENDDO
						ENDDO
						DEALLOCATE(tmpvari, STAT = status)
					CASE('Init_WaterSurfaceElevation')
						ALLOCATE(tmpvar(dim_vals(1)), STAT = status)
						call cg_array_read_f(IARRAY, tmpvar, IER)
						DO i = 1, ns2
							DO j = 1, nn
								countji = ((j-1)*ns2)+i
								e2(i,j) = tmpvar(countji)*100.
								IF(tmpvar(countji) > -999999.) THEN
									e2(i,j) = tmpvar(countji)*100.
								ELSE 
									e2(i,j) = 0.
								ENDIF
							ENDDO
						ENDDO
						DEALLOCATE(tmpvar, STAT = status)

					CASE('Init_VelocityX')
						ALLOCATE(tmpvar(dim_vals(1)), STAT = status)
						call cg_array_read_f(IARRAY, tmpvar, IER)
						DO i = 1, ns2
							DO j = 1, nn
								countji = ((j-1)*ns2)+i
								IF(tmpvar(countji) > -999999.) THEN
									iu(i,j) = tmpvar(countji)*100.
								ELSE 
									iu(i,j) = 0.
								ENDIF
							ENDDO
						ENDDO
						DEALLOCATE(tmpvar, STAT = status)

					CASE('Init_VelocityY')
						ALLOCATE(tmpvar(dim_vals(1)), STAT = status)
						call cg_array_read_f(IARRAY, tmpvar, IER)
						DO i = 1, ns2
							DO j = 1, nn
								countji = ((j-1)*ns2)+i
								IF(tmpvar(countji) > -999999.) THEN
									iv(i,j) = tmpvar(countji)*100.
								ELSE 
									iv(i,j) = 0.
								ENDIF
							ENDDO
						ENDDO
						DEALLOCATE(tmpvar, STAT = status)

				END SELECT
			END DO
			CALL cg_goto_f(FID, BID, IER, 'Zone_t', ZONES_ITER, 'end')

		END SUBROUTINE

		SUBROUTINE CGNS_Read_GridInputValues(FID, BID, IER, ZONES_ITER, USER_ITER)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BID, IER
		INTEGER, INTENT(IN) :: USER_ITER, ZONES_ITER
		INTEGER :: NARRAYS, IARRAY
		CHARACTER(LEN = 250) :: name
		INTEGER :: datatype, nndim
		INTEGER, DIMENSION(3) :: dim_vals
		INTEGER :: status, i, j, count
		INTEGER :: countij, countji

		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar, tmpvar2
		INTEGER, ALLOCATABLE, DIMENSION(:) :: tmpvari
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpcdregion


			call cg_goto_f(FID, BID, IER, 'Zone_t', ZONES_ITER,'UserDefinedData_t', USER_ITER, 'end')
			call cg_narrays_f(NARRAYS, IER)
			DO IARRAY = 1, NARRAYS
				call cg_array_info_f(iarray, name, datatype, nndim, dim_vals, ier)
					SELECT CASE(TRIM(name))
					CASE('Elevation')
						ALLOCATE(tmpvar(dim_vals(1)), STAT = status)
						call cg_array_read_f(IARRAY, tmpvar, IER)
						!First determin the elevOffset as the min elevation value
						DO i = 1, ns2
							DO j = 1, nn
								countji = ((j-1)*ns2)+i
								if(tmpvar(countji) < elevOffset) then
								    elevOffset = tmpvar(countji)
								endif
							ENDDO
						ENDDO
						
						DO i = 1, ns2
							DO j = 1, nn
								countji = ((j-1)*ns2)+i
								eta2(i,j) = (tmpvar(countji)-elevOffset)*100. !Convert to cm
							ENDDO
						ENDDO
						DEALLOCATE(tmpvar, STAT = status)

					CASE('Roughness')
						IF(cdtype == 0) THEN
							ALLOCATE(tmpvar(dim_vals(1)), STAT = status)
							call cg_array_read_f(IARRAY, tmpvar, IER)
							DO i = 1, ns2
								DO j = 1, nn
								    countji = ((j-1)*ns2)+i
									cd2(i,j) = tmpvar(countji)
									znaught2(i,j) = tmpvar(countji)*100.
								ENDDO
							ENDDO
							call setRoughnessIBC(tmpvar)
							DEALLOCATE(tmpvar, STAT = status)
						ELSE IF (cdtype == 1) THEN
						    ALLOCATE(tmpvar(dim_vals(1)), STAT = status)
							call cg_array_read_f(IARRAY, tmpvar, IER)
							call setRoughnessIBC(tmpvar)
							DO i = 1, ns2
								DO j = 1, nn
								    countji = ((j-1)*ns2)+i
									cd2(i,j) = tmpvar(countji)
									znaught2(i,j) = tmpvar(countji)*100.
								ENDDO
							ENDDO
							DEALLOCATE(tmpvar, STAT = status)
						ENDIF
					CASE('RoughnessRegionVals')
					    ALLOCATE(tmpvar(dim_vals(1)), STAT = status)
					    call alloc_roughnessvals(dim_vals(1))
					    call cg_array_read_f(IARRAY, tmpvar, IER)
						call setRoughnessIBCVals(tmpvar)
						DEALLOCATE(tmpvar, STAT = status)
					CASE('Sand_Fraction')
					    ALLOCATE(tmpvar(dim_vals(1)), STAT = status)
					    call cg_array_read_f(IARRAY, tmpvar, IER)
						DO i = 1, ns2
							DO j = 1, nn
								countji = ((j-1)*ns2)+i
								Fracs(i,j) = tmpvar(countji) 
							ENDDO
						ENDDO
						DEALLOCATE(tmpvar, STAT = status)
					CASE('Sand_Depth')
					    ALLOCATE(tmpvar(dim_vals(1)), STAT = status)
					    call cg_array_read_f(IARRAY, tmpvar, IER)
						DO i = 1, ns2
							DO j = 1, nn
								countji = ((j-1)*ns2)+i
								hfine(i,j) = tmpvar(countji)*100. !Convert to cm
							ENDDO
						ENDDO
						DEALLOCATE(tmpvar, STAT = status)
					CASE('Grainsize_d50')
					    ALLOCATE(tmpvar(dim_vals(1)), STAT = status)
					    call cg_array_read_f(IARRAY, tmpvar, IER)
						DO i = 1, ns2
							DO j = 1, nn
								countji = ((j-1)*ns2)+i
								d(i,j) = tmpvar(countji)*100. !Convert to cm
							ENDDO
						ENDDO
						DEALLOCATE(tmpvar, STAT = status)

    
				END SELECT
			END DO
			CALL cg_goto_f(FID, BID, IER, 'Zone_t', ZONES_ITER, 'end')
		END SUBROUTINE

		SUBROUTINE CGNS_Read_CentLine(FID, BID, IER, ZONES_ITER, USER_ITER)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BID, IER
		INTEGER, INTENT(IN) :: USER_ITER, ZONES_ITER
		INTEGER :: NARRAYS, IARRAY
		CHARACTER(LEN = 250) :: name
		INTEGER :: datatype, nndim
		INTEGER, DIMENSION(3) :: dim_vals
		DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: xc, yc
		REAL, ALLOCATABLE, DIMENSION(:) :: xct, yct, xot, yot
		INTEGER :: status

		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar, tmpvar2
		INTEGER, ALLOCATABLE, DIMENSION(:) :: tmpvari
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpcdregion
		REAL, ALLOCATABLE, DIMENSION(:) :: si, yp, temp, sout
		REAL, ALLOCATABLE, DIMENSION(:) :: xotmp, yotmp, potmp
		DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: phitmp
		REAL :: tens
		REAL :: tTheta, tdy, tdx, xext, yext, tds, tdr, tdphi

		INTEGER :: nsid, nnid, nzid
		INTEGER :: i, j, count, tmpnsext, tmpcdtype
		INTEGER :: numVBCPts, tcalcwet
		INTEGER :: numRegions
		INTEGER :: nclpts, nc

		REAL :: tmpwidth, tmpcd, slin, tmpevc
		REAL :: uu, vv, thmin
		REAL :: xc1, yc1, scals, stot
		REAL :: rsin, rcos
		DOUBLE PRECISION :: dphi, dx, dy

			call cg_goto_f(FID, BID, IER, 'Zone_t', ZONES_ITER,'UserDefinedData_t', USER_ITER, 'end')
			call cg_narrays_f(NARRAYS, IER)
			DO IARRAY = 1, NARRAYS
				call cg_array_info_f(iarray, name, datatype, nndim, dim_vals, ier)
!				nc = dim_vals(1)
					SELECT CASE(TRIM(name))
					CASE('CenterlineX')
!					    if(nsext > 0) then
!					        nc = dim_vals(1) + 1
!					        ALLOCATE(xc(dim_vals(1)+1), STAT = status)
!					    else
						    nc = dim_vals(1)
    						ALLOCATE(xc(dim_vals(1)), STAT = status)
 !   					endif
					    
						ALLOCATE(xct(dim_vals(1)), STAT = status)
						call cg_array_read_f(IARRAY, xct, IER)
					CASE('CenterlineY')
!					    if(nsext > 0) then
!						    ALLOCATE(yc(dim_vals(1)+1), STAT = status)
!						else
						    ALLOCATE(yc(dim_vals(1)), STAT = status)
!						endif
					        
						ALLOCATE(yct(dim_vals(1)), STAT = status)
						call cg_array_read_f(IARRAY, yct, IER)
					CASE('GridLength')
						call cg_array_read_f(IARRAY, mo2, IER)
						mo2 = mo2*100.
					CASE('GridWidth')
						call cg_array_read_f(IARRAY, tmpwidth, IER)
						tmpwidth = tmpwidth * 100.
					CASE('GridTension')
						call cg_array_read_f(IARRAY, tens, IER)
					CASE('GridCenterlineX')
						ALLOCATE(xot(dim_vals(1)), STAT = status)
						call cg_array_read_f(IARRAY, xot, IER)
					CASE('GridCenterlineY')
						ALLOCATE(yot(dim_vals(1)), STAT = status)
						call cg_array_read_f(IARRAY, yot, IER)
					END SELECT
					
			END DO
			CALL cg_goto_f(FID, BID, IER, 'Zone_t', ZONES_ITER, 'end')
! Create Grid from Centerline
!			nc = dim_vals(1)
            if(nsext > 0) then
				tdy = yct(nc-1)-yct(nc-2)
				tdx = xct(nc-1)-xct(nc-2)
				tTheta = atan2(tdy,tdx)
				tds = sqrt((xot(ns2)-xot(ns2-1))**2.+(yot(ns2)-yot(ns2-1))**2.)
				xext = xct(nc-1)+(tds*nsext*cos(tTheta))
		   		yext = yct(nc-1)+(tds*nsext*sin(tTheta))
			endif
            IF(nc == 2) THEN
                DEALLOCATE(xc)
                DEALLOCATE(yc)
                ALLOCATE(xc(3), STAT = status)
                ALLOCATE(yc(3), STAT = status)
                DO I = 1, 3
                    IF(I == 2) THEN
                        xc(I) = (xct(1)+((xct(2)-xct(1))/2))*100.
                        yc(I) = (yct(1)+((yct(2)-yct(1))/2))*100.
                    ELSE IF (i == 3) THEN
                       xc(I) = xct(2)*100.
                        yc(I) = yct(2)*100.
!                        xc(I) = xext*100.
!                        yc(I) = yext*100.
                   else 
                        xc(I) = xct(I)*100.
                        yc(I) = yct(I)*100.
                    END IF
                END DO
                nc = 3
            ELSE
			    DO I = 1, nc
				    xc(I) = xct(I)*100.
				    yc(I) = yct(I)*100.
!				    if(i == nc) then
!				        xc(i) = xext*100.
!				        yc(i) = yext*100.
!				    endif
			    ENDDO
			END IF

			DEALLOCATE(xct)
			DEALLOCATE(yct)

			ALLOCATE(si(nc))
			DO i = 1, nc
				if(i == 1) then
					si(i) = 0
				else 
						ds=((xc(i)-xc(i-1))**2.+(yc(i)-yc(i-1))**2.)**.5
						si(i) = si(i-1)+ds
				endif
			END DO

			stot = si(nc)
			xshift = xc(1)
			yshift = yc(1)
			DO j = 1, ns2
!				IF(j == 1) THEN
!					xshift = xot(j)*100.
!					yshift = yot(j)*100.
!				ENDIF
				xo2(j) = (xot(j)*100.)-xshift
				yo2(j) = (yot(j)*100.)-yshift
				w2(j) = tmpwidth
			ENDDO

			DEALLOCATE(xot)
			DEALLOCATE(yot)

			DO j = 1,nc
				xc(j) = xc(j)-xshift
				yc(j) = yc(j)-yshift
			ENDDO


		open(10, file='debugk.dat')
			slin=((xc(nc)**2)+yc(nc)**2)**.5
			fcos=xc(nc)/slin
			fsin=yc(nc)/slin
			ALLOCATE(yp(ns2))
			ALLOCATE(temp(ns2))
			ALLOCATE(sout(ns2))
			ALLOCATE(phitmp(nc))
			DO j = 1,ns2
				xc1 = xo2(j)*fcos + yo2(j)*fsin
				yc1 = yo2(j)*fcos - xo2(j)*fsin
				xo2(j) = xc1
				yo2(j) = yc1
			END DO
			DO j = 1,nc
				xc1 = xc(j)*fcos + yc(j)*fsin
				yc1 = yc(j)*fcos - xc(j)*fsin
				xc(j) = xc1
				yc(j) = yc1
			ENDDO
							
			DO j=1,ns2
				sout(j)=(j-1)*stot/(ns2-1)
			END DO

			DO j=2,nc
				dx = xc(j) - xc(j-1)
				dy = yc(j) - yc(j-1)
				if(dx == 0) then
					if(dy > 0) then
						phitmp(j) = acos(-1.)/2.
					else
						phitmp(j) = -1. * acos(-1.)/2.
					endif
				else
					phitmp(j) = atan2(dy,dx)
				endif
			ENDDO
				phitmp(1) = (2.*phitmp(2))-phitmp(3)
			ALLOCATE(xotmp(ns2))
			ALLOCATE(yotmp(ns2))
			ALLOCATE(potmp(ns2))

			call tspline(si,phitmp,nc,sout,potmp,ns2,tens,yp,temp)
			call tspline(si,xc,nc,sout,xotmp,ns2,tens,yp,temp)
			call tspline(si,yc,nc,sout,yotmp,ns2,tens,yp,temp)
			
			DO i = 1, ns2
			    xo(i) = xotmp(i)
			    yo(i) = yotmp(i)
			ENDDO
 			scals=stot/(ns2-1)
            
			DO 100 i=1,ns2
	!        w(i)=width*100.
				if(i.eq.1) go to 100
    !				dx=xo2(i)-xo2(i-1)
    !				dy=yo2(i)-yo2(i-1)
				    dx = xotmp(i)-xotmp(i-1)
				    dy = yotmp(i) - yotmp(i-1)
				    if(dx.eq.0) then
					    if(dy.gt.0) then
						    phirotation(i)=acos(-1.)/2.
					    elseif(dy.le.0) then
						    phirotation(i)=-1.*acos(-1.)/2.
					    endif
				    else
					    phirotation(i)=atan2(dy,dx)
				    endif
				    
	100		continue
			phirotation(1)=(2.*phirotation(2))-phirotation(3)
			do i = 1,ns2
				phi2(i) = phirotation(i)
			enddo
			do 140 i=2,ns2
				dphi=phirotation(i)-phirotation(i-1)
				if(dphi.eq.0) then
				    if(r2(i-1) < 0) then
				        r2(i) = -100000000.
				    else
					    r2(i)=100000000.
					endif
				else
					r2(i)=scals/dphi
				endif
	140		continue
			r2(1)=2.*r2(2)-r2(3)

		
		tdphi = phirotation(ns2)-phirotation(ns2-1)
		tdr = (tdphi)/nsext
		do i=ns2+1,ns2+nsext
		    dphi = ((ns2+nsext+1)-i)*tdr
		    phirotation(i) = phirotation(i-1)+dphi
		    
		    if(dphi == 0) then
		        if(r2(i-1) < 0) then
		            r2(i) = -10000000.
		        else
		            r2(i) = 10000000.
		         endif
		    else
		        r2(i) = scals/dphi
		    endif
		    dx = scals*cos(phirotation(i))
		    dy = scals*sin(phirotation(i))
		    xo(i) = xo(i-1) + dx
		    yo(i) = yo(i-1) + dy
		    dx = xo(i)-xo(i-1)
			dy = yo(i)-yo(i-1)
			if(dx.eq.0) then
				if(dy.gt.0) then
					phirotation(i)=acos(-1.)/2.
				elseif(dy.le.0) then
					phirotation(i)=-1.*acos(-1.)/2.
				endif
			else
				phirotation(i)=atan2(dy,dx)
			endif

		enddo


		do i=1,ns2+nsext
		write(10,*) xo(i), yo(i), phirotation(i), r2(i)
!		write(10,*) r2(i), xo2(i), yo2(i)
		enddo
		
		
1000    format(4f11.2)
		close(10)

		DEALLOCATE(xc)
		DEALLOCATE(yc)
		DEALLOCATE(si)
		DEALLOCATE(yp)
		DEALLOCATE(temp)
		DEALLOCATE(sout)
		DEALLOCATE(phitmp)
		DEALLOCATE(yotmp)
		DEALLOCATE(xotmp)
		DEALLOCATE(potmp)



		END SUBROUTINE
		
		SUBROUTINE CGNS_Read_TimeAndRatingValues(FID, BID, IER, ZONES_ITER, USER_ITER)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BID, IER
		INTEGER, INTENT(IN) :: USER_ITER, ZONES_ITER
		
		INTEGER :: NARRAYS, IARRAY
		CHARACTER(LEN = 250) :: name
		INTEGER :: datatype, nndim
		INTEGER, DIMENSION(3) :: dim_vals
		INTEGER :: status, i, j, count

		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar, tmpvar2, vx, vy
		INTEGER, ALLOCATABLE, DIMENSION(:) :: tmpvari
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpcdregion


			call cg_goto_f(FID, BID, IER, 'Zone_t', ZONES_ITER,'UserDefinedData_t', USER_ITER, 'end')
			call cg_narrays_f(NARRAYS, IER)
			
			
			DO IARRAY = 1, NARRAYS
				call cg_array_info_f(iarray, name, datatype, nndim, dim_vals, ier)
					SELECT CASE(TRIM(name))
					CASE('RatingCurveNumVals')
						ALLOCATE(tmpvari(dim_vals(1)), STAT = status)
						call cg_array_read_f(IARRAY, tmpvari, IER)
                        call setRatingCurveNumPts(dim_vals(1), tmpvari)
   						DEALLOCATE(tmpvari, STAT = status)

					CASE('TimeSeriesNumVals')
						ALLOCATE(tmpvari(dim_vals(1)), STAT = status)
						call cg_array_read_f(IARRAY, tmpvari, IER)
                        call setTimeSeriesNumPts(dim_vals(1), tmpvari)
   						DEALLOCATE(tmpvari, STAT = status)
					
					CASE('RatingCurveVals')
						ALLOCATE(tmpvar(dim_vals(1)), STAT = status)
						call cg_array_read_f(IARRAY, tmpvar, IER)
                        call setRatingCurves(dim_vals(1), tmpvar)
   						DEALLOCATE(tmpvar, STAT = status)
					
					CASE('TimeSeriesVals')
						ALLOCATE(tmpvar(dim_vals(1)), STAT = status)
						call cg_array_read_f(IARRAY, tmpvar, IER)
                        call setTimeSeries(dim_vals(1), tmpvar)
   						DEALLOCATE(tmpvar, STAT = status)
					
					CASE('TimeSeriesType')
						ALLOCATE(tmpvari(dim_vals(1)), STAT = status)
						call cg_array_read_f(IARRAY, tmpvari, IER)
                        call setTimeSeriesType(dim_vals(1), tmpvari)
   						DEALLOCATE(tmpvari, STAT = status)
					
					CASE('RatingCurveType')
						ALLOCATE(tmpvari(dim_vals(1)), STAT = status)
						call cg_array_read_f(IARRAY, tmpvari, IER)
                        call setRatingCurveType(dim_vals(1), tmpvari)
   						DEALLOCATE(tmpvari, STAT = status)

				END SELECT
			END DO
			CALL initRatingCurves()
			CALL initTimeSeries()
			CALL cg_goto_f(FID, BID, IER, 'Zone_t', ZONES_ITER, 'end')
						
		END SUBROUTINE


END MODULE