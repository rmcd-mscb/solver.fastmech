! SUBROUTINE TRIDAG
	MODULE TRIDIAG
	USE RIVVARMOD
	IMPLICIT NONE
	CONTAINS
 	subroutine tridag(a,b,c,d,e,n)
	INTEGER, INTENT(IN) :: n
 	REAL, DIMENSION(n) :: a,b,c,d,e,wk2
	REAL :: bet
	INTEGER :: j
	INTEGER(4) :: iret
	errorcode = 0
30	bet=b(1)
	e(1)=d(1)/bet
	do 35 j=2,n
	 wk2(j)=c(j-1)/bet
	 bet=b(j)-a(j)*wk2(j)
	 if(bet.eq.0 .or. isnan(bet)) then
	  write(6,*) 'oops',a(j),b(j),wk2(j),j,n
		errorcode = -1
	  return
	 endif
	 e(j)=(d(j)-a(j)*e(j-1))/bet
35	continue
	do 40 j=n-1,1,-1
40	e(j)=e(j)-wk2(j+1)*e(j+1)
50	return
	end SUBROUTINE
        
	END MODULE
