	SUBROUTINE ELEVDT(abh, bbh, hmax, Vsandsurf, Vsandsurfmax, &
                        Fsmin, Fs, area, h, he)
!

!	CALCULATE THE BED ELEVATION RELATIVE TO THE ORIGINAL	
!	POSITION
!
!
!
!	INPUT: abh, bbh, hmax, Vsandsurf, Vsandsurfmax, Fsmin Fs, dn, ds
!	OUTPUT: h, he
!
!	h>=0.0
!
!	examples if h = 0 then bed is at the original elevation
!			 if h = hmax  then bed elevation is original + hmax elevation
!
!	we suppose that the change in elevation for the entire cell eccurs when h>hmax
!
!      VARIABLES
!	
!	ah =alpha parameter for beta model for elevation to Sand fraction
!	bh =beta  parameter for beta model for elevation to Sand fraction
!	hmax = maximum elevation of largest gravel given	
!	Vsandsurf = total volume of sand deposited over gravel
!	Vsandsurfmax = minimum volume of sand to cover all gravel
!	Fsmin = the minimum fraction of sand 
!	Fs = fraction fo sand
!	h = sand deposit elevation
!	he = change in elevation of the entire bed cell
!	dn = dimension cell along n
!     ds=  dimension cell along s
!     area = area of cell
!
!       END

	IMPLICIT NONE


 	REAL		VS, Fsmin, tol
	REAL		Fs, CDF, av, bv, XX, abh, bbh
	REAL		Vsandsub, Vsandsurf
	REAL		Vsandsubmax, Vsandsurfmax
	REAL		Vsandtotal, DV
	REAL		area, hin, h, hmax, he
	INTEGER		Flag
	

	tol=0.00001



	IF	(Fs.LE.Fsmin) THEN

		h=0.0
		he=0.0         ! change in elevation

	ELSEIF (Fs.GE.0.99) THEN

		DV=Vsandsurf-Vsandsurfmax
		h=DV/(area)
		he=h			! change in elevation
		h=hmax+h	
	ELSE
		
		he=0.0           ! change in elevation
	
		CALL PPFBET (Fs, abh, bbh,tol,Flag, XX )
		h=XX
		h=h*hmax

		IF (FLAG.GE.1) THEN
						
			write(*,*) 'Error in Beta evaluation of sand depth'
			write(*,*) 'Error:', FLAG
			write(*,*) 'in subroutine PPFBET'					 

		ENDIF




	ENDIF






	RETURN
	END