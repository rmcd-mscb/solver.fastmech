	SUBROUTINE  INITIMETA( Fs , area, Fsmin , hmax &
                          , h, vsandsurf, vsandsurfmax, he)   
!
!	RUNS ONLY ONCE BEFORE SEDIMENT TRANSPORT CALCULATION BEGAN
!	CALCULATES THE ELEVATION OF TEH IMMOBILE GRAVEL BED FROM FS 
!	IF FS=1 HE=(VSANDSURF-VSANDSURFMAX)/(DN*DS)
!     IF FS<FSMIN THEN HE==0
!     IF FSMIN<FS<1 THEN HE==0
!
!	INPUTS: Fs, ds, dn, Fsmin, hmax, h
!
!	OUTPUTS: HE
!
!
!               VARIABLES
!
!
!	Fs = sand fraction
!	Fsmin = minimum sand fraction
!	area = area of cell
!	hmax = maximum height of gravel
!	h = height of the sand
!	Vsandsurf = volume of sand
!     Vsandsurfmax = minimum volume of sand required to cover gravel
!	he = elevation change due to sand deposit over immobile gravel bed
!	
!
!
!
!
	IMPLICIT NONE

	REAL	Fs , area ,  Fsmin , hmax , abh , bbh , H
	REAL    vsandsurf, vsandsurfmax
	REAL	HE

 

	IF (Fs.LE.Fsmin) THEN
		
		he=0.0


	ELSEIF (Fs.GE.1.0) THEN
			
			he=(VSANDSURF-VSANDSURFMAX)/(area)

			
	ELSE
		
		he=0.0
		
	
	ENDIF

	RETURN

	END
