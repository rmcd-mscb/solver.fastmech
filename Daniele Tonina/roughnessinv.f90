	real function roughnessinv(Droughmin, Dro,  P, Fs)
!-----------------------------------------
!	we used a logistic curve of shape y=(A1-A2)/(1+(x/x0)^P)+A2
!	where A1 and A2 are the roughnesses for Fs=0 and Fs=1 respectively
!	x0 is the location parameter which is set to 0.5,
!    		where a rapid change in roughness
!		behavior happens
!	P is the shape paramenter 
!	Fs is the areal fraction of sand
!	Dro present Z0 from calibration
!	Droughmin lower Z0 due to sand
!	roughnessinv roughness of the bed without sand (Z0max)


	IMPLICIT NONE

	REAL Droughmin, Dro, P
	REAL FS, x0
	
	x0=0.5


	roughnessinv=(Dro-Droughmin)*(1+(Fs/x0)**P)+Droughmin


				

	RETURN

	END
