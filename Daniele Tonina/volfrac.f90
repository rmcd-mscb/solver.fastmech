	SUBROUTINE VOLFRAC(av,bv, VS, Vsandsub, &
                        Vsandsurf,Vsandsurfmax, Vsandsubmax, Fsmin, Fs)
!
!
!
!	CALCULATE THE FRACTION OF SURFACE AREA COVERED BY SAND
!
!	INPUT: av, bv, VS, Vsandsub, Vsandsurf,Vsandsurfmax, Vsandsubmax, Fsmin
!	OUTPUT: Fs
!
!     Fs > Fsmin
!	Fs < 1
!
!	
!	av =alpha parameter for beta model for Volume to % Area
!	bv =beta  parameter for beta model for Volume to % Area
!	VS = CON*dt*area with area=ds*dn (dimension of the cell)	
!	Vsandsub = the volume of sand in the subsurface actve layer
!	Vsandsurf = the volume of sand above the gravel
!	Vsandsurfmax = the min volume of sand that cover all gravel
!	Vsandsubmax = the max volume of sand within the sublayer
!	Fsmin = the minimum fraction of sand 
!	Fs = fraction fo sand
!

	IMPLICIT NONE
	INTEGER, PARAMETER :: mp = KIND(1.0D0)
	INTEGER iflag
	REAL(kind = mp)	VS, Fsmin, eps
	REAL(kind = mp)	Fs, CDF, av, bv, X 
	REAL(kind = mp)	Vsandsub, Vsandsurf
	REAL(kind = mp)	Vsandsubmax, Vsandsurfmax
	REAL(kind = mp)	Vsandtotal, DV, VOID
	
	Vsandtotal=Vsandsub + Vsandsurf
	DV=Vsandtotal-VS
	
	IF (DV.LE.0.0) THEN          !situation non possible
		
		Fs=Fsmin
		Vsandsub=0.0
		Vsandsurf=0.0
		
!       sediment transport time is too larger! This cell cannot transport! 
	
	
	ELSE                      ! Erosion
	
		DV=Vsandsurf-VS
		
		IF (DV.EQ.0.0) THEN   ! just cleaned the surface

			Fs=Fsmin
			Vsandsurf=0.0
		
		ELSEIF (DV.lt.0.0) THEN   ! eroded surface and start winnowing the sublayer
			
			Fs=Fsmin
			Vsandsurf=0.0
			Vsandsub=Vsandsub+DV
		
		ELSE							!Deposition
			
			IF (Vsandsub.lt.Vsandsubmax) THEN !filling up the sublayer
				
				VOID=Vsandsubmax-Vsandsub

				IF (VOID.gt.DV) THEN          !partial filling
					
					Vsandsub=Vsandsub+DV 
					DV=0.0
					GOTO 333

				ELSEIF (VOID.eq.DV) THEN      !total filling
					
					Vsandsub=Vsandsubmax
					DV=0.0
					GOTO 333

				ELSE						!total filling with deposit on top

					Vsandsub=Vsandsubmax
					DV=DV-VOID
			
				ENDIF
			
			ENDIF
		    
			Vsandsurf=DV					  
				
			IF (Vsandsurf.ge.Vsandsurfmax) THEN !total deposition on top
				
				Fs=1
				Vsandsurf = vsandsurf - VS
			
			ELSE
				
				X=Vsandsurf/Vsandsurfmax      !Partial deposition = partial cover of sand
				eps=0.00001 !accuracy of the fraction area
        			 
				CALL CDFBET ( X, av, bv, EPS, IFLAG, CDF )
					
					IF (IFLAG.GE.1) THEN
						
						write(*,*) 'Error in Beta evaluation of area'
						write(*,*) 'Error:', IFLAG
						write(*,*) 'in subroutine CDFBET'					 

					ENDIF


				Fs=cdf

				IF (Fs.LE.Fsmin) Fs=Fsmin
				
			ENDIF
		ENDIF
	ENDIF

333	RETURN

	END
			

			
			
