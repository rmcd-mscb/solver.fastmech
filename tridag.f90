! SUBROUTINE TRIDAG
MODULE TRIDIAG
USE RIVVARMOD
IMPLICIT NONE
CONTAINS
subroutine tridag(a,b,c,d,e,n,code)
INTEGER, INTENT(IN) :: n
INTEGER, INTENT(INOUT) :: code
REAL(kind = mp), DIMENSION(n) :: a,b,c,d,e,wk2
REAL(kind = mp) :: bet
INTEGER :: j
INTEGER(4) :: iret
code = 0
IF(b(1).EQ.0.D0) THEN
    CODE=-1
    RETURN
END IF

bet=b(1)
e(1)=d(1)/bet
do j=2,n
    wk2(j)=c(j-1)/bet
    bet=b(j)-a(j)*wk2(j)
    if(bet.eq.0.D0) then
        write(6,*) 'oops',a(j),b(j),wk2(j),j,n
        code = -1
        return
    endif
    e(j)=(d(j)-a(j)*e(j-1))/bet
enddo
do j=n-1,1,-1
    e(j)=e(j)-wk2(j+1)*e(j+1)
enddo
return
end SUBROUTINE

END MODULE
