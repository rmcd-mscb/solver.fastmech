	SUBROUTINE INITHDT( Fs , abh, bbh, Fsmin , hmax , h)   
!
!	RUNS ONLY ONCE ABEFORE SEDIMENT TRANSPORT CALCULATION BEGAN
!	CALCULATES THE SAND HEIGHT 
!	IF FS=1 THEN KEEPS THE HEIGHT READ FROM INPUT
!     IF FS<FSMIN THEN H==0
!     IF FSMIN<FS<1 THEN CLAUCLATES THE HEIGHT FROM BETA MODEL
!
!	INPUTS: Fs, abh, bbh, Fsmin, hmax
!
!	OUTPUTS: H
!
!
!               VARIABLES
!
!
!	Fs = sand fraction
!	Fsmin = minimum sand fraction
!	abh = alpha parameter for Beta model of Sand height VS. Sand area
!	bbh = beta parameter for Beta model of Sand height VS. Sand area
!	hmax = maximum height of gravel
!
!
!	h = height of the sand
!	
!	tol = accuracy for the beta model
!	IFLAG = error number
!
!
!
!
	IMPLICIT NONE
	INTEGER IFLAG
	REAL	Fs , Fsmin , hmax, abh, bbh 
	REAL	TOL, H, XX

 

	IF (Fs.LE.Fsmin) THEN
		
		h=0.0
		Fs = Fsmin


	ELSEIF (Fs.GE.1.0) THEN
        IF(h.lt.(hmax-.001)) then
            h = hmax
!      WRITE(6,*) 'Sand_Depth needs to be larger than hmax' 
!      WRITE(6,*) ' if Sand_Fractions = 1'
!!      PAUSE
!!      EXIT
        ELSE
			h=h
	  ENDIF

			
	ELSE
		


		tol=0.00001
		CALL PPFBET (Fs, abh, bbh, tol, IFLAG, XX )

				IF (IFLAG.GE.1) THEN
						
					write(*,*) 'Error in Beta evaluation of area'
					write(*,*) 'Error:', IFLAG
					write(*,*) 'in subroutine CDFBET'					 

				ENDIF

		h=XX
		h=h*hmax
	
	ENDIF

	RETURN

	END
