	SUBROUTINE BETAVOL( hmax , abh, bbh, Vsandsurfmax, AV , BV, area) 
               
!
!	RUNS ONLY ONCE BEFORE SEDIMENT TRANSPORT CALCULATION BEGAN
!	IT CALCULATES THE ALPHA AND BETA PARAMETERS FOR THE MODEL VOLUME VS AREA
!	CALCULATIONS
!
!	INPUTS: Fs, Fsmin, hmax, abh, bbh
!			
!
!	OUTPUTS: Vsandsurfmax, av,  bv
!
!
!               VARIABLES
!
!
!	Fs = sand fraction
!	Fsmin = minimum sand fraction
!	hmax = maximum height of gravel
!	abh = alpha parameter for Beta model Sand Height  VS Sand area
!	bbh = beta parameter for Beta model Sand Height  VS Sand area
!	area= area of the cell
!	Vsandsurfmax = minimum sand volume that covers the gravel
!	av = alpha parameter for Beta model Sand Volume  VS Sand area
!	bv = beta parameter for Beta model Sand Volume  VS Sand area!	
!	
!	
!	
!	EPS = accuracy for the beta model
!	IFLAG = error number
!
!
!
!
	IMPLICIT NONE

	INTEGER i , IFLAG
	REAL	hmax, abh, bbh
	REAL    Vsandsurfmax, av, bv
	REAL	hh1, hh2, a1, a2, EPS, CDF
 	REAL	VVT,  VV,  std1, m1
	REAL	VTS(100),  VTSbeta(100)
	REAL	meandt, stdevdt
	REAL, INTENT(IN) :: area


	
	
	VVT=0.0
	hh1=0.0
	hh2=0.01
	EPS=0.000001
	a1=0.0
	DO 100 i=1,100
		
!		CALL CDFBET ( hh1, abh, bbh, EPS, IFLAG, CDF )
!		a1 = CDF
		CALL CDFBET ( hh2, abh, bbh, EPS, IFLAG, CDF )
		a2 = CDF

		VV=0.5*(hh1+hh2)*(a2-a1)*hmax
		VVT=VVT+VV
	    VTS(i)=hh2*a2*hmax-VVT
		
		hh1=hh1+0.01
		hh2=hh2+0.01
		a1=a2

100	CONTINUE
	
	Vsandsurfmax = VTS(100)

	DO 101 i = 1 , 100
	
		VTSbeta(i)=VTS(i)/Vsandsurfmax
	
101	CONTINUE
	
		

	m1=meandt(VTSbeta,100)
	std1=stdevdt(VTSbeta,100)


	
	av=m1*((m1*(1-m1)/(std1*std1))-1)
	bv=(1-m1)*(av/m1)
	
	Vsandsurfmax=Vsandsurfmax*area !added 10/16/2009
	
	IF (av.LE.0.0.OR.bv.LE.0.0) THEN

		write (*,*) 'error to find Beta parameter for Volume vs Area'

	ENDIF

	RETURN

	END


	
	real function  meandt(X,rows)
	
	IMPLICIT NONE

	INTEGER i, rows
	REAL  X(rows)
	REAL sum
	

	sum=0.0
	DO i=1,rows
		sum=sum+X(i)
	ENDDO
	
	meandt=sum/real(rows)
	
	RETURN

	END

	real function  stdevdt(X,rows)
	
	IMPLICIT NONE

	INTEGER i, rows
	REAL  X(rows), mean
	REAL sum, RN, meandt
	
	mean=meandt(X,rows)

	sum=0.0
	DO i=1,rows
		
		sum=sum+(X(i)-mean)*(X(i)-mean)

	ENDDO
	
	RN=real((rows-1))
	stdevdt=SQRT(sum/RN)
	
	
	RETURN

	END