MODULE RivCGNSOutptSubRtMod
USE RivVarMod
USE RivVarWMod
USE CSedMod_DT_SUSP
!USE CSedMod
USE TSplineMod
USE RivVarTimeMod
IMPLICIT NONE
!!DEC$ ATTRIBUTES REFERENCE, C, VARYING :: cg_goto_f
!!DEC$ ATTRIBUTES REFERENCE, C, VARYING :: cg_array_read_f
!!DEC$ ATTRIBUTES REFERENCE, C, VARYING :: cg_array_write_f

INCLUDE "cgnslib_f.h"
	CONTAINS

    SUBROUTINE CGNS_BaseIter(FID, BASES_ITER)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: FID, BASES_ITER
        INTEGER :: ier
    		CALL cg_biter_write_f(FID,BASES_ITER,'TimeIterValues',TSPrintCount,ier)
            call cg_goto_f(FID,BASES_ITER,ier,'BaseIterativeData_t',1,'end')
            call cg_array_write_f('TimeValues',RealSingle,1,TSPrintCount,TimeIncrements,ier)  
            call cg_array_write_f('DischValues', RealSingle, 1, TSPrintCount, DischIncrements, ier)      
    END SUBROUTINE
    
    SUBROUTINE CGNS_ZoneIter(FID, BASES_ITER, ZONES_ITER, INDEX)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: FID, BASES_ITER, ZONES_ITER, INDEX
        INTEGER, DIMENSION(2) :: idata
        INTEGER :: ier, ier2, status, i
        CHARACTER*32, Allocatable, Dimension(:) :: tmpsolname
        CHARACTER(LEN = 32) :: tmp_flowsol1D, tmp_flowsol2D, tmp_flowsol3D
        CHARACTER(LEN = 14) :: flowsol1D, flowsol2D, flowsol3D
         CHARACTER(LEN = 5) :: tsindex
       write(6,*)'in zoneiter', FID
!        flowsol1D = '1DFlowSolution'
!        flowsol2D = '2DFlowSolution'
!        flowsol3D = '3DFlowSolution'
!        
!        WRITE(tsindex, '(I5)') TSPrintCount
!        tsindex = ADJUSTL(tsindex)
!        tmp_flowsol1D = flowsol1D//TRIM(tsindex)
!        tmp_flowsol1D = TRIM(tmp_flowsol1D)
!
!        tmp_flowsol2D = flowsol2D//TRIM(tsindex)
!        tmp_flowsol2D = TRIM(tmp_flowsol2D)
! 
!        tmp_flowsol3D = flowsol3D//TRIM(tsindex)
!        tmp_flowsol3D = TRIM(tmp_flowsol3D)
!       
!            Allocate(tmpsolname(TSPrintCount+1), STAT = status)
!            Do i = 1, tsprintcount
!                SELECT CASE(INDEX)
!                CASE(1)
!                    WRITE(tsindex, '(I5)') i
!                    tsindex = ADJUSTL(tsindex)
!                    tmp_flowsol1D = flowsol1D//TRIM(tsindex)
!                    tmp_flowsol1D = TRIM(tmp_flowsol1D)                
!                    tmpsolname(i) = tmp_flowsol1D
!                CASE(2)
!                    WRITE(tsindex, '(I5)') i
!                    tsindex = ADJUSTL(tsindex)
!                    tmp_flowsol2D = flowsol2D//TRIM(tsindex)
!                    tmp_flowsol2D = TRIM(tmp_flowsol2D)                
!                    tmpsolname(i) = tmp_flowsol2D
!                
!                CASE(3)
!                    WRITE(tsindex, '(I5)') i
!                    tsindex = ADJUSTL(tsindex)
!                    tmp_flowsol3D = flowsol3D//TRIM(tsindex)
!                    tmp_flowsol3D = TRIM(tmp_flowsol3D)                
!                    tmpsolname(i) = tmp_flowsol3D
!                END SELECT
!            enddo
        ier2 = 0
            Allocate(tmpsolname(TSPrintCount+1), STAT = status)
!            Do i = 1, tsprintcount
!                WRITE(tsindex, '(I5)') i
!                tsindex = ADJUSTL(tsindex)
!                tmp_flowsol2D = flowsol2D//TRIM(tsindex)
!                tmp_flowsol2D = TRIM(tmp_flowsol2D)                
!                tmpsolname(i) = tmp_flowsol2D
!            enddo
    		CALL cg_ziter_write_f(FID,BASES_ITER,ZONES_ITER,'ZoneIterativeData',ier)
            call cg_goto_f(FID,BASES_ITER,ier,'Zone_t',ZONES_ITER,'ZoneIterativeData_t',1,'end')
            idata(1)=32
            idata(2)=TSPrintCount
!            call cg_array_write_f('FlowSolutionPointers',Character,2,idata,solnames,ier2)
!            deallocate(tmpsolname, STAT = status)
!             ier2 = 0
!    		CALL cg_ziter_write_f(FID,BASES_ITER,ZONES_ITER,'ZoneIterativeData',ier)
!            call cg_goto_f(FID,BASES_ITER,ier,'Zone_t',ZONES_ITER,'ZoneIterativeData_t',1,'end')
!            idata(1)=32
!            idata(2)=TSPrintCount
!            call cg_array_write_f('FlowSolutionPointers',Character,2,idata,solnames,ier2)
            select case(index)
            case(1)
                call cg_array_write_f('FlowSolutionPointers',Character,2,idata,solnames1D, ier)
            case(2)
                call cg_array_write_f('FlowSolutionPointers',Character,2,idata,solnames,ier2)
            case(3)
               call cg_array_write_f('FlowSolutionPointers',Character,2,idata,solnames3D,ier2)
           END SELECT
            deallocate(tmpsolname, STAT = status)
           ! add SimulationType
!            call cg_simulation_type_write_f(FID,BASES_ITER,TimeAccurate,ier)            
END SUBROUTINE
   
	SUBROUTINE CGNS_Write_HelixStrength(FID, BID, ZID, F2DSOL_ITER)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BID, ZID
		INTEGER, INTENT(IN) :: F2DSOL_ITER
		INTEGER :: FSOL_ITER, IER, STATUS, I, J, K, COUNT
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar1

		REAL :: dx, dy, rcos, rsin, xx, yy
		REAL :: mag1, mag2, tmp1

			ALLOCATE(tmpvar1(nn*ns2*nz), STAT = status)
!			DO i = 1, ns2
!				rcos = cos(phirotation(i))
!				rsin = sin(phirotation(i))
!				DO j = 1,nn
!					count = ((i-1)*nn)+j
!!					if(ibc(i,j) == 0) then
!!						tmpvar1(count) = -999999.
!!					else
!						xx = atan2(vz(i,j,1), uz(i,j,1))*(360./(2.*3.14159))
!						yy = atan2(vz(i,j,nz), uz(i,j, nz))*(360./(2.*3.14159))	
!						tmpvar1(count) = xx-yy
!!					endif
!				END DO
!			ENDDO
			DO j = 1, nn
				DO i = 1,ns2
				    rcos = cos(phirotation(i))
				    rsin = sin(phirotation(i))
!					count = ((i-1)*nn)+j
					count = i+(j-1)*ns2
!					if(ibc(i,j) == 0) then
!						tmpvar1(count) = -999999.
!					else
!                        mag1 = sqrt((uz(i,j,1)*uz(i,j,1))+(vz(i,j,1)*vz(i,j,1)))
!                        mag2 = sqrt((uz(i,j,nz)*uz(i,j,nz))+(vz(i,j,nz)*vz(i,j,nz))) 
!                        tmp1 = (uz(i,j,1)*uz(i,j,nz))+(vz(i,j,1)*vz(i,j,nz))
                        
                        
						xx = atan2(vz(i,j,2), uz(i,j,2))*(360./(2.*3.14159))
						yy = atan2(vz(i,j,nz), uz(i,j, nz))*(360./(2.*3.14159))	
						tmpvar1(count) = xx-yy
!						tmpvar1(count) = acos(tmp1/(mag1*mag2))*(360./(2.*3.14159))
!					endif
				END DO
			ENDDO
			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'HelixStrength', tmpvar1, FSOL_ITER, IER)
			DEALLOCATE(tmpvar1, STAT = status)

	END SUBROUTINE

	SUBROUTINE CGNS_Write_Velocity3D(FID, BID, ZID, F3DSOL_ITER)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BID, ZID
		INTEGER, INTENT(IN) :: F3DSOL_ITER
		INTEGER :: FSOL_ITER, IER, STATUS, I, J, K, COUNT
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar1, tmpvar2
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar3, tmpvar4
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar5

		REAL :: ux, uy, rcos, rsin, xx, yy

			ALLOCATE(tmpvar1(nn*ns2*nz), STAT = status)
			ALLOCATE(tmpvar2(nn*ns2*nz), STAT = status)
			ALLOCATE(tmpvar3(nn*ns2*nz), STAT = status)
			ALLOCATE(tmpvar4(nn*ns2*nz), STAT = status)
			ALLOCATE(tmpvar5(nn*ns2*nz), STAT = status)
			
 !          DO i = 1, ns2
 !                  rcos = cos(phirotation(i))
 !                  rsin = sin(phirotation(i))
 !                  DO j = 1,nn
 !                      DO k = 1,nz
 !                          count = ((i-1)*nn*nz)+((j-1)*nz)+k
 !  !                       if(ibc(i,j) == 0) then
 !  !                           tmpvar1(count) = -999999.
 !  !                           tmpvar2(count) = -999999.
 !  !                           tmpvar3(count) = -999999.
 !  !                           tmpvar4(count) = -999999.
 !  !                       else
 !  
 !                              ux = (uz(i,j,k)*rcos - vz(i,j,k)*rsin)/100.
 !                              uy = (uz(i,j,k)*rsin + vz(i,j,k)*rcos)/100.
 !                              xx = ux*fcos - uy*fsin
 !                              yy = ux*fsin + uy*fcos
 !                              tmpvar1(count) = xx
 !                              tmpvar2(count) = yy
 !                              tmpvar3(count) = uz(i,j,k)/100.
 !                              tmpvar4(count) = vz(i,j,k)/100.
 !  !                       endif
 !                      END DO
 !                  END DO
 !              ENDDO
 			DO k = 1, nz
				DO j = 1,nn
					DO i = 1,ns2
				        rcos = cos(phirotation(i))
				        rsin = sin(phirotation(i))
!						count = ((i-1)*nn*nz)+((j-1)*nz)+k
						count = i + (j-1)*ns2 + (k-1)*ns2*nn
!						if(ibc(i,j) == 0) then
!							tmpvar1(count) = -999999.
!							tmpvar2(count) = -999999.
!							tmpvar3(count) = -999999.
!							tmpvar4(count) = -999999.
!						else

							ux = (uz(i,j,k)*rcos - vz(i,j,k)*rsin)/100.
							uy = (uz(i,j,k)*rsin + vz(i,j,k)*rcos)/100.
							xx = ux*fcos - uy*fsin
							yy = ux*fsin + uy*fcos
							tmpvar1(count) = xx
							tmpvar2(count) = yy
							tmpvar3(count) = uz(i,j,k)/100.
							tmpvar4(count) = vz(i,j,k)/100.
							tmpvar5(count) = 0.0
!						endif
					END DO
				END DO
			ENDDO
			IF(IO_VelXY) THEN
			CALL cg_field_write_f(FID, BID, ZID, F3DSOL_ITER, RealSingle, 'VelocityX', tmpvar1, FSOL_ITER, IER)
			CALL cg_field_write_f(FID, BID, ZID, F3DSOL_ITER, RealSingle, 'VelocityY', tmpvar2, FSOL_ITER, IER)
			CALL cg_field_write_f(FID, BID, ZID, F3DSOL_ITER, RealSingle, 'VelocityZ', tmpvar5, FSOL_ITER, IER)
			ENDIF
			IF(IO_VelSN) THEN
			CALL cg_field_write_f(FID, BID, ZID, F3DSOL_ITER, RealSingle, 'VelocityS', tmpvar3, FSOL_ITER, IER)
			CALL cg_field_write_f(FID, BID, ZID, F3DSOL_ITER, RealSingle, 'VelocityN', tmpvar4, FSOL_ITER, IER)
			ENDIF
			DEALLOCATE(tmpvar1, STAT = status)
			DEALLOCATE(tmpvar2, STAT = status)
			DEALLOCATE(tmpvar3, STAT = status)
			DEALLOCATE(tmpvar4, STAT = status)
			DEALLOCATE(tmpvar5, STAT = status)

	END SUBROUTINE

	SUBROUTINE CGNS_Write_CoordinateZ(FID, BID, ZID)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BID, ZID
		INTEGER :: CID, IER, STATUS, I, J, K, COUNT
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar1, tmpvar2, tmpvar3

		ALLOCATE (tmpvar1(ns2*nn*nz), STAT = STATUS)
 !      DO I= 1, ns2
 !              DO J = 1, nn
 !                  DO K = 1, nz
 !                      count = ((i-1)*nn*nz)+((j-1)*nz)+k
 !                      tmpvar1(count) = zz(i,j,k)/100.
 !  !                   tmpvar2(count) = xo2(i)*((nn/2.)-j)*sin(phirotation(i))
 !  !                   tmpvar3(count) = yo2(i)*((nn/2.)-j)*cos(phirotation(i))
 !                  ENDDO
 !              ENDDO
 !          ENDDO
        DO k= 1, nz
 			DO J = 1, nn
				DO i = 1, ns2
!					count = ((i-1)*nn*nz)+((j-1)*nz)+k
						count = i + (j-1)*ns2 + (k-1)*ns2*nn
					tmpvar1(count) = zz(i,j,k)/100.
				ENDDO
			ENDDO
		ENDDO		
		call cg_coord_write_f(FID,BID,ZID,RealSingle,'CoordinateZ',tmpvar1, CID, IER)

		DEALLOCATE(tmpvar1, STAT = status)
	
	END SUBROUTINE
	SUBROUTINE CGNS_Write_GRIDZ(FID, BID, ZID, F3DSOL_ITER)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BID, ZID
		INTEGER, INTENT(IN) :: F3DSOL_ITER
		INTEGER :: CID, IER, FSOL_ITER, STATUS, I, J, K, COUNT
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar1, tmpvar2, tmpvar3

		ALLOCATE (tmpvar1(ns2*nn*nz), STAT = STATUS)
 !      DO I= 1, ns2
 !              DO J = 1, nn
 !                  DO K = 1, nz
 !                      count = ((i-1)*nn*nz)+((j-1)*nz)+k
 !                      tmpvar1(count) = zz(i,j,k)/100.
 !  !                   tmpvar2(count) = xo2(i)*((nn/2.)-j)*sin(phirotation(i))
 !  !                   tmpvar3(count) = yo2(i)*((nn/2.)-j)*cos(phirotation(i))
 !                  ENDDO
 !              ENDDO
 !          ENDDO
        DO k= 1, nz
 			DO J = 1, nn
				DO i = 1, ns2
!					count = ((i-1)*nn*nz)+((j-1)*nz)+k
						count = i + (j-1)*ns2 + (k-1)*ns2*nn
					tmpvar1(count) = (zz(i,j,k)/100.) + elevOffset
				ENDDO
			ENDDO
		ENDDO		
!		call cg_coord_write_f(FID,BID,ZID,RealSingle,'GRIDZ',tmpvar1, CID, IER)
		CALL cg_field_write_f(FID, BID, ZID, F3DSOL_ITER, RealSingle, 'GridZ', tmpvar1, FSOL_ITER, IER)

		DEALLOCATE(tmpvar1, STAT = status)
	
	END SUBROUTINE

	SUBROUTINE CGNS_Write_ShearStress(FID, BID, ZID, F2DSOL_ITER)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BID, ZID
		INTEGER, INTENT(IN) :: F2DSOL_ITER
		INTEGER :: FSOL_ITER, IER, STATUS, I, J, COUNT
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar1, tmpvar2
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar3, tmpvar4
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar5

		REAL :: ux, uy, rcos, rsin, xx, yy
		
		REAL (KIND = 8):: ws, nu, t
		REAL (KIND = 8) :: g, csf, p, s, dd
		REAl :: sws, ustar, tb
		t = 20.
		g = 9.81
		csf = 0.7
		p = 6.
		s = 2.65

			ALLOCATE(tmpvar1(nn*ns2), STAT = status)
			ALLOCATE(tmpvar2(nn*ns2), STAT = status)
			ALLOCATE(tmpvar3(nn*ns2), STAT = status)
			ALLOCATE(tmpvar4(nn*ns2), STAT = status)
			ALLOCATE(tmpvar5(nn*ns2), STAT = status)
 !          DO i=1,ns2
 !                  rcos = cos(phirotation(i))
 !                  rsin = sin(phirotation(i))
 !                  DO j= 1,nn
 !                      count = ((i-1)*nn)+j
 !  !                   if(ibc(i,j) == 0) then
 !  !                       tmpvar1(count) = -999999.
 !  !                       tmpvar2(count) = -999999.
 !  !                       tmpvar3(count) = -999999.
 !  !                       tmpvar4(count) = -999999.
 !  !                   else
 !  
 !                          ux = (taus(i,j)*rcos - taun(i,j)*rsin)/10.
 !                          uy = (taus(i,j)*rsin + taun(i,j)*rcos)/10.
 !                          xx = ux*fcos - uy*fsin
 !                          yy = ux*fsin + uy*fcos
 !                          tmpvar1(count) = xx
 !                          tmpvar2(count) = yy
 !                          tmpvar3(count) = taus(i,j)/10.
 !                          tmpvar4(count) = taun(i,j)/10.
 !  !                   endif
 !                  END DO
 !              END DO
 
 			DO j=1,nn
				DO i= 1,ns2
				    IF(CALCCSED) THEN
                        nu = 1.792D-6/(1.0D0+3.37D-2*t+2.21D-4*t*t) ! Kinematic viscosity.
                        dd = d(i,j)/100.
				        ws = falveld(dd, nu, s, g, p, csf)
				        sws = ws
				    ENDIF
				    rcos = cos(phirotation(i))
				    rsin = sin(phirotation(i))
!					count = ((i-1)*nn)+j
					count = i + (j-1)*ns2 
!					if(ibc(i,j) == 0) then
!						tmpvar1(count) = -999999.
!						tmpvar2(count) = -999999.
!						tmpvar3(count) = -999999.
!						tmpvar4(count) = -999999.
!					else

						ux = (taus(i,j)*rcos - taun(i,j)*rsin)/10.
						uy = (taus(i,j)*rsin + taun(i,j)*rcos)/10.
						xx = ux*fcos - uy*fsin
						yy = ux*fsin + uy*fcos
						tmpvar1(count) = xx
						tmpvar2(count) = yy
						tmpvar3(count) = taus(i,j)/10.
						tmpvar4(count) = taun(i,j)/10.
						IF(CALCCSED) THEN
						    tb = sqrt((tmpvar3(count)*tmpvar3(count)) + (tmpvar4(count)*tmpvar4(count)))
						    ustar = sqrt(tb/1000.)
!						    denom = sqrt((tmpvar3(count)*tmpvar3(count) + tmpvar4(count)*tmpvar4(count))/1000.)
						    if(ibc(i,j) == 0.or.tb == 0) then
						        tmpvar5(count) = 0
						    else
						    tmpvar5(count) = sws/(0.408*ustar)
						    endif
						ENDIF
						
!					endif
				END DO
			END DO
			Write(6,*) 'Grainsize', dd, 'Settling Velocity', sws
			
			IF (CALCCSED) THEN
			    CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'Rouse_Number', tmpvar5, FSOL_ITER, IER)
			ENDIF
			IF(IO_ShearXY) THEN
			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'ShearStressX', tmpvar1, FSOL_ITER, IER)
			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'ShearStressY', tmpvar2, FSOL_ITER, IER)
			ENDIF
			IF(IO_ShearSN) THEN
			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'ShearStressS', tmpvar3, FSOL_ITER, IER)
			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'ShearStressN', tmpvar4, FSOL_ITER, IER)
			ENDIF
			DEALLOCATE(tmpvar1, STAT = status)
			DEALLOCATE(tmpvar2, STAT = status)
			DEALLOCATE(tmpvar3, STAT = status)
			DEALLOCATE(tmpvar4, STAT = status)
			DEALLOCATE(tmpvar5, STAT = status)

	END SUBROUTINE
	SUBROUTINE CGNS_Write_RouseNumber(FID, BID, ZID, F2DSOL_ITER)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BID, ZID
		INTEGER, INTENT(IN) :: F2DSOL_ITER
		INTEGER :: FSOL_ITER, IER, STATUS, I, J, COUNT
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar1, tmpvar2
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar3, tmpvar4
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar5

		REAL :: ux, uy, rcos, rsin, xx, yy
		
		REAL (KIND = 8):: ws, nu, t
		REAL (KIND = 8) :: g, csf, p, s, dd
		REAl :: sws, denom
		t = 20.
		g = 9.81
		csf = 0.7
		p = 6.
		s = 2.65

			ALLOCATE(tmpvar1(nn*ns2), STAT = status)
			ALLOCATE(tmpvar2(nn*ns2), STAT = status)
			ALLOCATE(tmpvar3(nn*ns2), STAT = status)
			ALLOCATE(tmpvar4(nn*ns2), STAT = status)
			ALLOCATE(tmpvar5(nn*ns2), STAT = status)
 
 			DO j=1,nn
				DO i= 1,ns2
                    nu = 1.792D-6/(1.0D0+3.37D-2*t+2.21D-4*t*t) ! Kinematic viscosity.
                    dd = d(i,j)/100.
				    ws = falveld(dd, nu, s, g, p, csf)
				    sws = ws
				    rcos = cos(phirotation(i))
				    rsin = sin(phirotation(i))
					count = i + (j-1)*ns2 

						ux = (taus(i,j)*rcos - taun(i,j)*rsin)/10.
						uy = (taus(i,j)*rsin + taun(i,j)*rcos)/10.
						xx = ux*fcos - uy*fsin
						yy = ux*fsin + uy*fcos
						tmpvar1(count) = xx
						tmpvar2(count) = yy
						tmpvar3(count) = taus(i,j)/10.
						tmpvar4(count) = taun(i,j)/10.
						denom = sqrt((tmpvar3(count)*tmpvar3(count) + tmpvar4(count)*tmpvar4(count))/1000.)
						if(ibc(i,j) == 0.or.denom == 0) then
						    tmpvar5(count) = 0
						else
						tmpvar5(count) = sws/(0.408*denom)
						endif
!					endif
				END DO
			END DO
			Write(6,*) 'Grainsize', dd, 'Settling Velocity', sws
			
			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'Rouse_Number', tmpvar5, FSOL_ITER, IER)
			
			IF(IO_ShearXY) THEN
			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'ShearStressX', tmpvar1, FSOL_ITER, IER)
			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'ShearStressY', tmpvar2, FSOL_ITER, IER)
			ENDIF
			IF(IO_ShearSN) THEN
			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'ShearStressS', tmpvar3, FSOL_ITER, IER)
			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'ShearStressN', tmpvar4, FSOL_ITER, IER)
			ENDIF
			DEALLOCATE(tmpvar1, STAT = status)
			DEALLOCATE(tmpvar2, STAT = status)
			DEALLOCATE(tmpvar3, STAT = status)
			DEALLOCATE(tmpvar4, STAT = status)
			DEALLOCATE(tmpvar5, STAT = status)

	END SUBROUTINE
	
!REAL FUNCTION falveld(d,nu,s,g,p,csf)
!  IMPLICIT NONE
!! Computes particle fall velocity according to W. Dietrich, WRR 18(6),
!! 1615-1626, 1982. falveld in m/s.
!! Input parameters:
!!   d    particle diameter (m)
!!   nu   kinematic viscosity (m^2/s)
!!   s    specific gravity
!!   g    acceleration due to gravity (m/s^2)
!!   p    Powers roundness factor
!!   csf  Corey shape factor (must be >= 0.2)
!!
!! NOTE: for a sphere, p=6.0, csf=1.0; for natural sediments csf=0.7.
!! Function by Francisco Simoes, March 1999.
!
!  REAL, INTENT(IN) :: d,nu,s,g,p,csf
!  REAL :: dstar,logdstar,r1,r2,r3,tanhdstar,wstar
!  REAL, PARAMETER :: one=1.0
!
!  dstar = ((s-1.)*g*(d**3.))/(nu*nu)
!  logdstar = LOG10(dstar)
!  tanhdstar = TANH(logdstar-4.6)
!  r1 = -3.76715+1.92944*logdstar-9.815-2*logdstar*logdstar- &
!       5.75-3*logdstar**3+5.6-4*logdstar**4
!  r2 = ALOG10(1.-(1.-csf)/8.5-1.)
!  r2 = r2-tanhdstar*(1.-csf)**2.3
!  r2 = r2+0.30*(0.50-csf)*(1.-csf)*(1.-csf)*(logdstar-4.60)
!  r3 = (6.5-1.-csf*tanhdstar/2.830)**(1.+(3.50-p)/2.50)
!  wstar = r3*10.00**(r1+r2)
!  falveld = (wstar*(s-1.)*g*nu)**(1./3.00)
!
!END FUNCTION falveld
!REAL(KIND=8) FUNCTION falveld(d,nu,s,g,p,csf)
!  IMPLICIT NONE
!! Computes particle fall velocity according to W. Dietrich, WRR 18(6),
!! 1615-1626, 1982. falveld in m/s.
!! Input parameters:
!!   d    particle diameter (m)
!!   nu   kinematic viscosity (m^2/s)
!!   s    specific gravity
!!   g    acceleration due to gravity (m/s^2)
!!   p    Powers roundness factor
!!   csf  Corey shape factor (must be >= 0.2)
!!
!! NOTE: for a sphere, p=6.0, csf=1.0; for natural sediments csf=0.7.
!! Function by Francisco Simoes, March 1999.
!
!  REAL(KIND=8), INTENT(IN) :: d,nu,s,g,p,csf
!  REAL(KIND=8) :: dstar,logdstar,r1,r2,r3,tanhdstar,wstar
!  REAL(KIND=8), PARAMETER :: one=1.0
!
!  dstar = (s-one)*g*d**3/(nu*nu)
!  logdstar = LOG10(dstar)
!  tanhdstar = DTANH(logdstar-4.6D0)
!  r1 = -3.76715D0+1.92944*logdstar-9.815D-2*logdstar*logdstar- &
!       5.75D-3*logdstar**3+5.6D-4*logdstar**4
!  r2 = LOG10(one-(one-csf)/8.5D-1)
!  r2 = r2-tanhdstar*(one-csf)**2.3
!  r2 = r2+0.3D0*(0.5D0-csf)*(one-csf)*(one-csf)*(logdstar-4.6D0)
!  r3 = (6.5D-1-csf*tanhdstar/2.83D0)**(one+(3.5D0-p)/2.5D0)
!  wstar = r3*10.0D0**(r1+r2)
!  falveld = (wstar*(s-one)*g*nu)**(one/3.0D0)
!
!END FUNCTION falveld



	SUBROUTINE CGNS_Write_Velocity(FID, BID, ZID, F2DSOL_ITER)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BID, ZID
		INTEGER, INTENT(IN) :: F2DSOL_ITER
		INTEGER :: FSOL_ITER, IER, STATUS, I, J, COUNT, next
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar1, tmpvar2
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar3, tmpvar4

		REAL :: ux, uy, rcos, rsin, xx, yy
		
     REAL, ALLOCATABLE, DIMENSION(:) :: KEG,VORT,STRAIN
     REAL, ALLOCATABLE, DIMENSION(:,:) :: UU,DUU2DS,DVDS,DVDN,DUDS,DUDN
     REAL:: SC
    
    if(ShowGridExt == 1) then
        next = ns2+nsext
    else
        next = ns2
    endif
    ALLOCATE(KEG(nn*next), STAT = status)
    ALLOCATE(VORT(nn*next), STAT = status)
    ALLOCATE(STRAIN(nn*next), STAT = status)
    ALLOCATE(UU(next, nn), STAT = status)
    ALLOCATE(DUU2DS(next, nn), STAT = status)
    ALLOCATE(DVDS(next, nn), STAT = status)
    ALLOCATE(DVDN(next, nn), STAT = status)
    ALLOCATE(DUDS(next, nn), STAT = status)
    ALLOCATE(DUDN(next, nn), STAT = status)
    

    DO j = 1,nn
        DO i = 1,next
          xx = u(i,j)
          yy = vout(i,j)
          uu(i,j) = sqrt((xx*xx)+(yy*yy))
        ENDDO
    ENDDO

     DO 1778 j=1,nn
        DO 1778 i=1,next
          count = i + (j-1)*next 

          SC=rn(I,J)
!          xx = u(i,j)
!          yy = v(i,j)
!          uu(i,j) = sqrt((xx*xx)+(yy*yy))
          IF (ibc(i,j).EQ.0) THEN   
              KEG(count)    = 0.0
              Vort(count)   = 0.0
              Strain(count) = 0.0
              uu(i,j) = 0.
              GO TO 1778
          ENDIF
                
    
          IF (I.eq.1) THEN   !at upper boundary
                DUU2DS(I,J)=(UU(I+1,J)*UU(I+1,J)-UU(I,J)*UU(I,J))/(2.0*ds*SC)
               dVds(I,J)= (V(I+1,J)-V(I,J))/(ds*SC)  
                dUds(I,J)= (U(I+1,J)-U(I,J))/(ds*SC)
                IF(J.EQ.1) THEN  
                   dUdn(I,J)= (U(I,J+1)-U(I,J))/(dn)
                   dVdn(I,J)= (V(I,J+1)-V(I,J))/(dn)
                ELSEIF (J.EQ.nn) THEN       
                  dUdn(I,J)= (0.0-U(I,J))/(dn)
                   dVdn(I,J)= (0.0-V(I,J))/(dn)
                ELSE
                   dUdn(I,J)= (U(I,J+1)-U(I,J-1))/(2.0*dn)
                   dVdn(I,J)= (V(I,J+1)-V(I,J-1))/(2.0*dn)
                ENDIF
               GO TO 1643
          ELSEIF (I.EQ.ns) THEN !at lower boundary
               DUU2ds(I,J)=DUU2ds(I-1,J)
               dVds(I,J)= dVds(I-1,J)  
                dUds(I,J)= dUds(I-1,J)
                IF(J.EQ.1) THEN  
                   dUdn(I,J)= (U(I,J+1)-U(I,J))/(dn)
                   dVdn(I,J)= (V(I,J+1)-V(I,J))/(dn)
                ELSEIF (J.EQ.nn) THEN       
                  dUdn(I,J)= (0.0-U(I,J))/(dn)
                   dVdn(I,J)= (0.0-V(I,J))/(dn)
                ELSE
                   dUdn(I,J)= (U(I,J+1)-U(I,J-1))/(2.0*dn)
                   dVdn(I,J)= (V(I,J+1)-V(I,J-1))/(2.0*dn)
                ENDIF
                GO TO 1643
           ELSE
                DUU2ds(I,J)=(UU(I+1,J)*UU(I+1,J)-UU(I-1,J)*UU(I-1,J))/(2.0*ds*SC)
                dVds(I,J)= (V(I+1,J)-V(I-1,J))/(2.0*ds*SC)  
                dUds(I,J)= (U(I+1,J)-U(I-1,J))/(2.0*ds*SC)
                IF(J.EQ.1) THEN  
                   dUdn(I,J)= (U(I,J+1)-U(I,J))/(dn)
                   dVdn(I,J)= (V(I,J+1)-V(I,J))/(dn)
                ELSEIF (J.EQ.nn) THEN       
                  dUdn(I,J)= (0.0-U(I,J))/(dn)
                   dVdn(I,J)= (0.0-V(I,J))/(dn)
                ELSE
                   dUdn(I,J)= (U(I,J+1)-U(I,J-1))/(2.0*dn)
                   dVdn(I,J)= (V(I,J+1)-V(I,J-1))/(2.0*dn)
                ENDIF  
           ENDIF          
               
1643      KEG(count)    =DUU2ds(I,J)/2.0
          Vort(count)   = DVDS(I,J)-DUDN(I,J)
          Strain(count) = ABS(DUDS(I,J))+ABS(DUDN(I,J))+ABS(DVDS(I,J))+ABS(DVDN(I,J))
                              
1778  CONTINUE
      IF(IO_KEG) THEN
       CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'KineticEnergyGradient', KEG, FSOL_ITER, IER)
      ENDIF
	  IF(IO_Vort) THEN
	    CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'Vorticity', Vort, FSOL_ITER, IER)
	  ENDIF
	  IF(IO_VelStrain) THEN
	    CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'VelocityStrain', Strain, FSOL_ITER, IER)
	  ENDIF

    DEALLOCATE(KEG, STAT = status)
    DEALLOCATE(VORT, STAT = status)
    DEALLOCATE(STRAIN, STAT = status)
    DEALLOCATE(UU, STAT = status)
    DEALLOCATE(DUU2DS, STAT = status)
    DEALLOCATE(DVDS, STAT = status)
    DEALLOCATE(DVDN, STAT = status)
    DEALLOCATE(DUDS, STAT = status)
    DEALLOCATE(DUDN, STAT = status)

        if(ShowGridExt == 1) then
            next = ns2+nsext
        else
            next = ns2
        endif
			ALLOCATE(tmpvar1(nn*next), STAT = status)
			ALLOCATE(tmpvar2(nn*next), STAT = status)
			ALLOCATE(tmpvar3(nn*next), STAT = status)
			ALLOCATE(tmpvar4(nn*next), STAT = status)
 
              DO j= 1,nn
				DO i=1,next
 				    rcos = cos(phirotation(i))
				    rsin = sin(phirotation(i))
!					count = ((i-1)*nn)+j
					count = i + (j-1)*next 

						ux = (u(i,j)*rcos - vout(i,j)*rsin)/100.
						uy = (u(i,j)*rsin + vout(i,j)*rcos)/100.
						xx = ux*fcos - uy*fsin
						yy = ux*fsin + uy*fcos
						tmpvar1(count) = xx
						tmpvar2(count) = yy
						tmpvar3(count) = u(i,j)/100.
						tmpvar4(count) = vout(i,j)/100.
				END DO
			END DO
			
			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'VelocityX', tmpvar1, FSOL_ITER, IER)
			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'VelocityY', tmpvar2, FSOL_ITER, IER)
			IF(IO_VelSN) THEN
			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'VelocityS', tmpvar3, FSOL_ITER, IER)
			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'VelocityN', tmpvar4, FSOL_ITER, IER)
			ENDIF
			DEALLOCATE(tmpvar1, STAT = status)
			DEALLOCATE(tmpvar2, STAT = status)
			DEALLOCATE(tmpvar3, STAT = status)
			DEALLOCATE(tmpvar4, STAT = status)

	END SUBROUTINE
	
	SUBROUTINE CGNS_Write_InitVelocity(FID, BID, ZID, F2DSOL_ITER)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BID, ZID
		INTEGER, INTENT(IN) :: F2DSOL_ITER
		INTEGER :: FSOL_ITER, IER, STATUS, I, J, COUNT, next
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar1, tmpvar2
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar3, tmpvar4

		REAL :: ux, uy, rcos, rsin, xx, yy

        if(ShowGridExt == 1) then
            next = ns2+nsext
        else
            next = ns2
        endif
			ALLOCATE(tmpvar1(nn*next), STAT = status)
			ALLOCATE(tmpvar2(nn*next), STAT = status)
			ALLOCATE(tmpvar3(nn*next), STAT = status)
			ALLOCATE(tmpvar4(nn*next), STAT = status)
 
              DO j= 1,nn
				DO i=1,next
!					count = ((i-1)*nn)+j
					count = i + (j-1)*next 
						tmpvar3(count) = u(i,j)/100.
						tmpvar4(count) = v(i,j)/100.
				END DO
			END DO
			
			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'VelocityInitS', tmpvar3, FSOL_ITER, IER)
			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'VelocityInitN', tmpvar4, FSOL_ITER, IER)
			DEALLOCATE(tmpvar3, STAT = status)
			DEALLOCATE(tmpvar4, STAT = status)

	END SUBROUTINE
	
	SUBROUTINE CGNS_Write_AveWSElev(FID, BID, ZID, F2DSOL_ITER)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BID, ZID
		INTEGER, INTENT(IN) :: F2DSOL_ITER
		INTEGER :: FSOL_ITER, IER, STATUS, I, J, COUNT
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar1

		REAL :: U, D

			ALLOCATE(tmpvar1(ns2), STAT = status)
			DO i=1,ns2
				tmpvar1(i) = (eav(i)/100.) + elevOffset
			END DO
			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'WaterSurfaceElevation', tmpvar1, FSOL_ITER, IER)
			DEALLOCATE(tmpvar1, STAT = status)

		END SUBROUTINE
	
	SUBROUTINE CGNS_Write_RadiusCurv(FID, BID, ZID, F2DSOL_ITER)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BID, ZID
		INTEGER, INTENT(IN) :: F2DSOL_ITER
		INTEGER :: FSOL_ITER, IER, STATUS, I, J, COUNT
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar1

		REAL :: U, D

			ALLOCATE(tmpvar1(ns2), STAT = status)
			DO i=1,ns2
				tmpvar1(i) = r(i)/100.
			END DO
			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'Radius Of Curvature', tmpvar1, FSOL_ITER, IER)
			DEALLOCATE(tmpvar1, STAT = status)

		END SUBROUTINE

	SUBROUTINE CGNS_Write_XSConvergence(FID, BID, ZID, F2DSOL_ITER)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BID, ZID
		INTEGER, INTENT(IN) :: F2DSOL_ITER
		INTEGER :: FSOL_ITER, IER, STATUS, I, J, COUNT
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar1

		REAL :: U, D

			ALLOCATE(tmpvar1(ns2), STAT = status)
			DO i=1,ns2
				tmpvar1(i) = ((qp(i)/q)-1.)*100.
			END DO
			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'XSConvergence', tmpvar1, FSOL_ITER, IER)
			DEALLOCATE(tmpvar1, STAT = status)

		END SUBROUTINE

	SUBROUTINE CGNS_Write_Convergence(FID, BASES_ITER)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BASES_ITER
		INTEGER :: IER, STATUS

			CALL cg_goto_f(FID, BASES_ITER, IER, 'end')
			CALL cg_convergence_write_f(itm, 'NormalizedXSMassConv', IER)
			CALL cg_goto_f(FID, BASES_ITER, IER, 'ConvergenceHistory_t', 1, 'end')
			CALL cg_array_write_f('NormalizedXSMass', RealSingle, 1, itm, qprms, IER)

	END SUBROUTINE

	SUBROUTINE CGNS_Write_FroudeNumber(FID, BID, ZID, F2DSOL_ITER)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BID, ZID
		INTEGER, INTENT(IN) :: F2DSOL_ITER
		INTEGER :: FSOL_ITER, IER, STATUS, I, J, COUNT
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar1

		REAL :: UU, D

			ALLOCATE(tmpvar1(nn*ns2), STAT = status)
 !          DO i=1,ns2
 !                  DO j= 1,nn
 !                  count = ((i-1)*nn)+j
 !                  if(ibc(i,j) == 0) then
 !                      tmpvar1(count) = 0.
 !                  else
 !                      uu = sqrt(u(i,j)**2 +vout(i,j)**2)/100.
 !                      d = hl(i,j)/100.
 !                      tmpvar1(count) = uu/sqrt(9.8*d)
 !                  endif
 !                  END DO
 !              END DO
 			
			DO j= 1,nn
				DO i=1,ns2
!				count = ((i-1)*nn)+j
				count = i + (j-1)*ns2 
				if(ibc(i,j) == 0) then
					tmpvar1(count) = 0.
				else
					uu = sqrt(u(i,j)**2 +vout(i,j)**2)/100.
					d = hl(i,j)/100.
					tmpvar1(count) = uu/sqrt(9.8*d)
				endif
				END DO
			END DO
			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'FroudeNumber', tmpvar1, FSOL_ITER, IER)
			DEALLOCATE(tmpvar1, STAT = status)

		END SUBROUTINE

	SUBROUTINE CGNS_Write_Depth(FID, BID, ZID, F2DSOL_ITER)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BID, ZID
		INTEGER, INTENT(IN) :: F2DSOL_ITER
		INTEGER :: FSOL_ITER, IER, STATUS, I, J, COUNT
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar1

			ALLOCATE(tmpvar1(nn*ns2), STAT = status)
 !          DO i=1,ns2
 !                  DO j= 1,nn
 !                  count = ((i-1)*nn)+j
 !  !               if(ibc(i,j) == 0) then
 !  !                   tmpvar1(count) = -999999.
 !  !               else
 !                      tmpvar1(count) = hl(i,j)/100.
 !  !               endif
 !                  END DO
 !              END DO
 
			DO j= 1,nn
				DO i=1,ns2
				count = i + (j-1)*ns2 
!				count = ((i-1)*nn)+j
!				if(ibc(i,j) == 0) then
!					tmpvar1(count) = -999999.
!				else
					tmpvar1(count) = hl(i,j)/100.
!				endif
				END DO
			END DO
			
			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'Depth', tmpvar1, FSOL_ITER, IER)
			DEALLOCATE(tmpvar1, STAT = status)

		END SUBROUTINE

	SUBROUTINE CGNS_Write_StressDiv(FID, BID, ZID, F2DSOL_ITER)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BID, ZID
		INTEGER, INTENT(IN) :: F2DSOL_ITER
		INTEGER :: FSOL_ITER, IER, STATUS, I, J, COUNT
		REAL :: min, max
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar1
            
            min = 1e12
            max = -1e12
			ALLOCATE(tmpvar1(nn*ns2), STAT = status)
 !          DO i=1,ns2
 !                  DO j= 1,nn
 !                  count = ((i-1)*nn)+j
 !  !               if(ibc(i,j) == 0 ) then
 !  !                   tmpvar1(count) = -999999.
 !  !               else
 !                      tmpvar1(count) = con(i,j)
 !  !                   if(tmpvar1(count) < min) min = tmpvar1(count)
 !  !                   if(tmpvar1(count) > max) max = tmpvar1(count)
 !  !               endif
 !                  END DO
 !              END DO
 			
			DO j= 1,nn
				DO i=1,ns2
				count = i + (j-1)*ns2 
!				count = ((i-1)*nn)+j
!				if(ibc(i,j) == 0 ) then
!					tmpvar1(count) = -999999.
!				else
					tmpvar1(count) = con(i,j)
!					if(tmpvar1(count) < min) min = tmpvar1(count)
!					if(tmpvar1(count) > max) max = tmpvar1(count)
!				endif
				END DO
			END DO
			
			if(CALCCSED == 1) then
			    CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'Erosion Rate', tmpvar1, FSOL_ITER, IER)
			else
			    CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'Shear Stress Divergence', tmpvar1, FSOL_ITER, IER)
			endif
			DEALLOCATE(tmpvar1, STAT = status)

		END SUBROUTINE
	
	SUBROUTINE CGNS_Write_SuspSed(FID, BID, ZID, F2DSOL_ITER)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BID, ZID
		INTEGER, INTENT(IN) :: F2DSOL_ITER
		INTEGER :: FSOL_ITER, IER, STATUS, I, J, COUNT
		REAL :: min, max
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar1
            
            min = 1e12
            max = -1e12
			ALLOCATE(tmpvar1(nn*ns2), STAT = status)
 			
			DO j= 1,nn
				DO i=1,ns2
				count = i + (j-1)*ns2 
					tmpvar1(count) = totcs(i,j)
				END DO
			END DO
			    
			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'Susp. Sed. Conc.', tmpvar1, FSOL_ITER, IER)
			
			DEALLOCATE(tmpvar1, STAT = status)

		END SUBROUTINE

	SUBROUTINE CGNS_Write_WSElev(FID, BID, ZID, F2DSOL_ITER)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BID, ZID
		INTEGER, INTENT(IN) :: F2DSOL_ITER
		INTEGER :: FSOL_ITER, IER, STATUS, I, J, COUNT, next
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar1
        
        if(ShowGridExt == 1) then
            next = ns2+nsext
        else
            next = ns2
        endif
			ALLOCATE(tmpvar1(nn*next), STAT = status)
 !          DO i=1,next
 !                  DO j= 1,nn
 !                  count = ((i-1)*nn)+j
 !  !               if(ibc(i,j) == 0) then
 !  !                   tmpvar1(count) = -999999.
 !  !               else
 !                      tmpvar1(count) = e(i,j)/100.
 !  !               endif
 !                  END DO
 !              END DO
 
			DO j= 1,nn
				DO i=1,next
				count = i + (j-1)*next 
!				count = ((i-1)*nn)+j
!				if(ibc(i,j) == 0) then
!					tmpvar1(count) = -999999.
!				else
					tmpvar1(count) = (e(i,j)/100.) + elevOffset
!				endif
				END DO
			END DO
			
			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'WaterSurfaceElevation', tmpvar1, FSOL_ITER, IER)
			DEALLOCATE(tmpvar1, STAT = status)

		END SUBROUTINE

	SUBROUTINE CGNS_Write_CD(FID, BID, ZID, F2DSOL_ITER)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BID, ZID
		INTEGER, INTENT(IN) :: F2DSOL_ITER
		INTEGER :: FSOL_ITER, IER, STATUS, I, J, COUNT
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar1

			ALLOCATE(tmpvar1(nn*ns2), STAT = status)
 !          DO i=1,ns2
 !                  DO j= 1,nn
 !                  count = ((i-1)*nn)+j
 !  !               if(ibc(i,j) == 0) then
 !  !                   tmpvar1(count) = -999999.
 !  !               else
 !                      tmpvar1(count) = cd(i,j)
 !  !               endif
 !                  END DO
 !              END DO
 !  
 			DO j= 1,nn
				DO i=1,ns2
				count = i + (j-1)*ns2 
!				count = ((i-1)*nn)+j
!				if(ibc(i,j) == 0) then
!					tmpvar1(count) = -999999.
!				else
					tmpvar1(count) = cd(i,j)
!				endif
				END DO
			END DO
			
			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'Drag_Coefficient', tmpvar1, FSOL_ITER, IER)
			DEALLOCATE(tmpvar1, STAT = status)

		END SUBROUTINE
		
	SUBROUTINE CGNS_Write_GridIBC3D(FID, BID, ZID, F3DSOL_ITER)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BID, ZID
		INTEGER, INTENT(IN) :: F3DSOL_ITER
		INTEGER :: FSOL_ITER, IER, STATUS, I, J, K, COUNT, next
		INTEGER, ALLOCATABLE, DIMENSION(:) :: tmpIntVar1
			ALLOCATE(tmpIntVar1(nn*ns2*nz), STAT = status)
 !          DO i=1,next
 !                  DO j= 1,nn
 !                  count = ((i-1)*nn)+j
 !  !               if(ibc(i,j) == 0) then
 !  !                   tmpIntVar1(count) = -999999 
 !  !               else
 !                      tmpIntVar1(count) = ibc(i,j)
 !  !               endif
 !                  END DO
 !              END DO
            DO k = 1,nz
			    DO j= 1,nn
				    DO i=1,ns2
				    	count = i + (j-1)*ns2 + (k-1)*ns2*nn

!				        count = i + (j-1)*next 
!				count = ((i-1)*nn)+j
!				if(ibc(i,j) == 0) then
!					tmpIntVar1(count) = -999999	
!				else
					tmpIntVar1(count) = ibc(i,j)
!				endif
                    END DO
				END DO
			END DO
			
			CALL cg_field_write_f(FID, BID, ZID, F3DSOL_ITER, Integer, 'IBC', tmpIntVar1, FSOL_ITER, IER)
			write(6,*) 'Write IBC Solution', FID, BID, ZID, F3DSOL_ITER, FSOL_ITER, IER
			DEALLOCATE(tmpIntVar1, STAT = status)

		END SUBROUTINE
				
	SUBROUTINE CGNS_Write_GridIBC(FID, BID, ZID, F2DSOL_ITER)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BID, ZID
		INTEGER, INTENT(IN) :: F2DSOL_ITER
		INTEGER :: FSOL_ITER, IER, STATUS, I, J, COUNT, next
		INTEGER, ALLOCATABLE, DIMENSION(:) :: tmpIntVar1
            IF(ShowGridExt == 1) then
                next = ns2+nsext
            else
                next = ns2
            endif
			ALLOCATE(tmpIntVar1(nn*next), STAT = status)
 !          DO i=1,next
 !                  DO j= 1,nn
 !                  count = ((i-1)*nn)+j
 !  !               if(ibc(i,j) == 0) then
 !  !                   tmpIntVar1(count) = -999999 
 !  !               else
 !                      tmpIntVar1(count) = ibc(i,j)
 !  !               endif
 !                  END DO
 !              END DO
 
			DO j= 1,nn
				DO i=1,next
				count = i + (j-1)*next 
!				count = ((i-1)*nn)+j
!				if(ibc(i,j) == 0) then
!					tmpIntVar1(count) = -999999	
!				else
					tmpIntVar1(count) = ibc(i,j)
!				endif
				END DO
			END DO
			
			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, Integer, 'IBC', tmpIntVar1, FSOL_ITER, IER)
			write(6,*) 'Write IBC Solution', FID, BID, ZID, F2DSOL_ITER, FSOL_ITER, IER
			DEALLOCATE(tmpIntVar1, STAT = status)

		END SUBROUTINE
		
	SUBROUTINE CGNS_Write_HAREA(FID, BID, ZID, F2DSOL_ITER)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BID, ZID
		INTEGER, INTENT(IN) :: F2DSOL_ITER
		INTEGER :: FSOL_ITER, IER, STATUS, I, J, COUNT, next
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpVar
            IF(ShowGridExt == 1) then
                next = ns2+nsext
            else
                next = ns2
            endif
			ALLOCATE(tmpvar(nn*next), STAT = status)
			DO j= 1,nn
				DO i=1,next
				count = i + (j-1)*next 
					tmpvar(count) = harea(i,j)/10000.
				END DO
			END DO
			
			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'Habitat Area', tmpvar, FSOL_ITER, IER)
			DEALLOCATE(tmpvar, STAT = status)
			
	END SUBROUTINE CGNS_Write_HAREA

	SUBROUTINE CGNS_Write_ICON(FID, BID, ZID, F2DSOL_ITER)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BID, ZID
		INTEGER, INTENT(IN) :: F2DSOL_ITER
		INTEGER :: FSOL_ITER, IER, STATUS, I, J, COUNT, next
		INTEGER, ALLOCATABLE, DIMENSION(:) :: tmpIntVar1
            IF(ShowGridExt == 1) then
                next = ns2+nsext
            else
                next = ns2
            endif
			ALLOCATE(tmpIntVar1(nn*next), STAT = status)
			DO j= 1,nn
				DO i=1,next
				count = i + (j-1)*next 
					tmpIntVar1(count) = icon(i,j)
				END DO
			END DO
			
			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, Integer, 'Connectivity', tmpIntVar1, FSOL_ITER, IER)
			DEALLOCATE(tmpIntVar1, STAT = status)
			
	END SUBROUTINE CGNS_Write_ICON
	
	SUBROUTINE CGNS_Write_SedFluxComp(FID, BID, ZID, F2DSOL_ITER)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BID, ZID
		INTEGER, INTENT(IN) :: F2DSOL_ITER
		INTEGER :: FSOL_ITER, IER, STATUS, I, J, COUNT
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar1, tmpvar2
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar3, tmpvar4

		REAL :: ux, uy, rcos, rsin, xx, yy

			ALLOCATE(tmpvar1(nn*ns2), STAT = status)
			ALLOCATE(tmpvar2(nn*ns2), STAT = status)
 !          DO i=1,ns2
 !                  rcos = cos(phirotation(i))
 !                  rsin = sin(phirotation(i))
 !                  DO j= 1,nn
 !                      count = ((i-1)*nn)+j
 !  !                   if(ibc(i,j) == 0) then
 !  !                       tmpvar1(count) = -999999.
 !  !                       tmpvar2(count) = -999999.
 !  !                       tmpvar3(count) = -999999.
 !  !                       tmpvar4(count) = -999999.
 !  !                   else
 !  
 !                          ux = (qs(i,j)*rcos - qn(i,j)*rsin)
 !                          uy = (qs(i,j)*rsin + qn(i,j)*rcos)
 !                          xx = ux*fcos - uy*fsin
 !                          yy = ux*fsin + uy*fcos
 !                          tmpvar1(count) = xx
 !                          tmpvar2(count) = yy
 !  !                   endif
 !                  END DO
 !              END DO
 
			DO j= 1,nn
				DO i=1,ns2
				    rcos = cos(phirotation(i))
				    rsin = sin(phirotation(i))
    				count = i + (j-1)*ns2 
!					count = ((i-1)*nn)+j
!					if(ibc(i,j) == 0) then
!						tmpvar1(count) = -999999.
!						tmpvar2(count) = -999999.
!						tmpvar3(count) = -999999.
!						tmpvar4(count) = -999999.
!					else

						ux = (qs(i,j)*rcos - qn(i,j)*rsin)
						uy = (qs(i,j)*rsin + qn(i,j)*rcos)
						xx = ux*fcos - uy*fsin
						yy = ux*fsin + uy*fcos
						tmpvar1(count) = xx/10000.
						tmpvar2(count) = yy/10000.
!					endif
				END DO
			END DO
			
			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'SedFluxX', tmpvar1, FSOL_ITER, IER)
			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'SedFluxY', tmpvar2, FSOL_ITER, IER)
			DEALLOCATE(tmpvar1, STAT = status)
			DEALLOCATE(tmpvar2, STAT = status)
			DEALLOCATE(tmpvar3, STAT = status)
			DEALLOCATE(tmpvar4, STAT = status)

	END SUBROUTINE
	
	SUBROUTINE CGNS_Write_Elevation(FID, BID, ZID, F2DSOL_ITER)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BID, ZID
		INTEGER, INTENT(IN) :: F2DSOL_ITER
		INTEGER :: FSOL_ITER, IER, STATUS, I, J, COUNT, next
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar1

        if(ShowGridExt == 1) then
            next = ns2+nsext
        else
            next = ns2
        endif
			ALLOCATE(tmpvar1(nn*next), STAT = status)
 !          DO i=1,next
 !                  DO j= 1,nn
 !                  count = ((i-1)*nn)+j
 !  !               if(ibc(i,j) == 0) then
 !  !                   tmpvar1(count) = -999999.
 !  !               else
 !                      tmpvar1(count) = eta(i,j)/100.
 !  !               endif
 !                  END DO
 !              END DO
 
			DO j= 1,nn
				DO i=1,next
    			    count = i + (j-1)*next 
!				count = ((i-1)*nn)+j
!				if(ibc(i,j) == 0) then
!					tmpvar1(count) = -999999.
!				else
					tmpvar1(count) = (eta(i,j)/100.) + elevOffset
!				endif
				END DO
			END DO


			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'Elevation', tmpvar1, FSOL_ITER, IER)
			DEALLOCATE(tmpvar1, STAT = status)

	END SUBROUTINE
	
	SUBROUTINE CGNS_Write_SandDepth(FID, BID, ZID, F2DSOL_ITER)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BID, ZID
		INTEGER, INTENT(IN) :: F2DSOL_ITER
		INTEGER :: FSOL_ITER, IER, STATUS, I, J, COUNT, next
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar1

        if(ShowGridExt == 1) then
            next = ns2+nsext
        else
            next = ns2
        endif
			ALLOCATE(tmpvar1(nn*next), STAT = status)
 
			DO j= 1,nn
				DO i=1,next
    			    count = i + (j-1)*next 
					tmpvar1(count) = hfine(i,j)/100.
				END DO
			END DO


			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'Sand_Depth', tmpvar1, FSOL_ITER, IER)
			DEALLOCATE(tmpvar1, STAT = status)

	END SUBROUTINE	
	
	SUBROUTINE CGNS_Write_SandFraction(FID, BID, ZID, F2DSOL_ITER)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BID, ZID
		INTEGER, INTENT(IN) :: F2DSOL_ITER
		INTEGER :: FSOL_ITER, IER, STATUS, I, J, COUNT, next
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar1

        if(ShowGridExt == 1) then
            next = ns2+nsext
        else
            next = ns2
        endif
			ALLOCATE(tmpvar1(nn*next), STAT = status)
 
			DO j= 1,nn
				DO i=1,next
    			    count = i + (j-1)*next 
					tmpvar1(count) = Fracs(i,j)
				END DO
			END DO


			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'Sand_Fraction', tmpvar1, FSOL_ITER, IER)
			DEALLOCATE(tmpvar1, STAT = status)

	END SUBROUTINE
!		SUBROUTINE CGNS_Write_LSub(FID, BID, ZID, F2DSOL_ITER)
!		IMPLICIT NONE
!		INTEGER, INTENT(IN) :: FID, BID, ZID
!		INTEGER, INTENT(IN) :: F2DSOL_ITER
!		INTEGER :: FSOL_ITER, IER, STATUS, I, J, COUNT, next
!		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar1
!
!        if(ShowGridExt == 1) then
!            next = ns2+nsext
!        else
!            next = ns2
!        endif
!			ALLOCATE(tmpvar1(nn*next), STAT = status)
! 
!			DO j= 1,nn
!				DO i=1,next
!    			    count = i + (j-1)*next 
!					tmpvar1(count) = lsub(i,j)
!				END DO
!			END DO
!
!
!			CALL cg_field_write_f(FID, BID, ZID, F2DSOL_ITER, RealSingle, 'L_Sub', tmpvar1, FSOL_ITER, IER)
!			DEALLOCATE(tmpvar1, STAT = status)
!
!	END SUBROUTINE
	
	SUBROUTINE CGNS_WRITE_2DGRIDEXTENSION(FID, BID, ZID)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BID
		INTEGER, INTENT(INOUT) :: ZID
		INTEGER :: IER, STATUS, I, J, COUNT, GID, SID
		INTEGER :: CELLDIM, PHYSDIM
		INTEGER, DIMENSION(2,3) :: isize
		REAL, ALLOCATABLE, DIMENSION(:) :: tmpvar1, tmpvar2
        REAL :: rcos, rsin, ux, uy, uxnew, uynew
		
		ALLOCATE(tmpvar1(nn*(ns2+nsext)), STAT = status)
		ALLOCATE(tmpvar2(nn*(ns2+nsext)), STAT = status)
        ZID = 0
		
		isize(1,1) = ns2+nsext
		isize(2,1) = nn
		
		isize(1,2) =  isize(1,1)-1
		isize(2,2) = isize(2,1)-1
		
		isize(1,3) = 0
		isize(2,3) = 0
		
		physdim = 2
		celldim = 2
		
		CALL cg_zone_write_f(FID, BID, "CurvilinearOrthogonal_Ext", isize, Structured, ZID, IER) 
		
 !      do i=1,ns2+nsext
 !              rcos=cos(phirotation(i))
 !              rsin=sin(phirotation(i))
 !  
 !              do j=1,nn
 !                  count = ((i-1)*nn)+j
 !                  ux = xo(i)+(dn)*(nm-j)*sin(phirotation(i))
 !                  uy = yo(i)-(dn)*(nm-j)*cos(phirotation(i))
 !  !                ux = (xo(i)-xshift)/100.
 !  !                uy = (yo(i)-yshift)/100.
 !                  uxnew=ux*fcos-uy*fsin
 !                  uynew=ux*fsin+uy*fcos
 !  !                tmpvar1(count) = (uxnew)+(dn/100.)*(nm-j)*sin(phirotation(i))
 !  !                tmpvar2(count) = (uynew)-(dn/100.)*(nm-j)*cos(phirotation(i))
 !                  tmpvar1(count) = (uxnew+xshift)/100.
 !                  tmpvar2(count) = (uynew+yshift)/100.
 !  
 !              enddo
 !          enddo
 
		do j=1,nn
			do i=1,ns2+nsext
			    rcos=cos(phirotation(i))
                rsin=sin(phirotation(i))
                count = i + (j-1)*(ns2+nsext)
!				count = ((i-1)*nn)+j
				ux = xo(i)+(dn)*(nm-j)*sin(phirotation(i))
				uy = yo(i)-(dn)*(nm-j)*cos(phirotation(i))
!                ux = (xo(i)-xshift)/100.
!                uy = (yo(i)-yshift)/100.
                uxnew=ux*fcos-uy*fsin
                uynew=ux*fsin+uy*fcos
!                tmpvar1(count) = (uxnew)+(dn/100.)*(nm-j)*sin(phirotation(i))
!                tmpvar2(count) = (uynew)-(dn/100.)*(nm-j)*cos(phirotation(i))
                tmpvar1(count) = (uxnew+xshift)/100.
                tmpvar2(count) = (uynew+yshift)/100.

			enddo
		enddo


        CALL cg_coord_write_f(FID, BID, ZID, RealSingle, "CoordinateX", tmpvar1, GID, IER)
        CALL cg_coord_write_f(FID, BID, ZID, RealSingle, "CoordinateY", tmpvar2, GID, IER)
        
 !       Call cg_sol_write_f(FID,  BID, ZID, "2DFlowSolution_Ext", Vertex, SID, IER)
        
		
	END SUBROUTINE
	
	SUBROUTINE CGNS_WRITE_3DGRIDEXTENSION(FID, BID, ZID)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: FID, BID, ZID
		INTEGER :: IER, STATUS, I, J, COUNT
		INTEGER :: CELLDIM, PHYSDIM
		INTEGER, DIMENSION(3,3) :: isize
		
		isize(1,1) = ns2+nsext
		isize(2,1) = nn
		isize(3,1) = nz
		
		isize(1,2) =  isize(1,1)-1
		isize(2,2) = isize(2,1)-1
		isize(3,2) = isize(3,1)-1
		
		isize(1,3) = 0
		isize(2,3) = 0
		isize(3,3) = 0
		
		physdim = 3
		celldim = 3
		
	END SUBROUTINE
	
	END MODULE
