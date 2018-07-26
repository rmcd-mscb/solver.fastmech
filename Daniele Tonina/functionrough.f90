	real function roughness(A1, A2, P, FS)
!-----------------------------------------
! we used a logistic curve of shape y=(A1-A2)/(1+(x/x0)^P)+A2
! where A1 and A2 are the roughnesses for Fs=0 and Fs=1 respectively
! x0 is the location parameter which is set to 0.5,
!    where a rapid change in roughness
!    behavior happens
! p is the shape paramenter which we set to 7
! Fs is the areal fraction of sand	


    INTEGER, PARAMETER :: mp = KIND(1.0D0)
	REAL(kind = mp) A1, A2, P
	REAL(kind = mp) FS
	
	x0=0.5


	roughness=(A1-A2)/(1+(FS/x0)**P)+A2

	RETURN

	END
