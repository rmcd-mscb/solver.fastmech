Module Median
! http://jean-pierre.moreau.pagesperso-orange.fr/Fortran/tmdian_f90.txt
CONTAINS
!*******************************************************
!* Given an array X of N numbers, returns their median *
!* value XMED. The array X is not modified, and is     *
!* accessed sequentially in each consecutive pass.     *
!*******************************************************
SUBROUTINE MDIAN1(X,N,XMED)
    IMPLICIT NONE                             ! SRC
    INTEGER, PARAMETER :: mp = KIND(1.0D0)    ! SRC
    REAL(KIND=mp)  XMED                       ! SRC  Median.f90(9)
    INTEGER        N                          ! SRC  Median.f90(10)
    REAL(KIND=mp)  A                          ! SRC  Median.f90(16)
    REAL(KIND=mp)  EPS                        ! SRC  Median.f90(17)
    REAL(KIND=mp)  AP                         ! SRC  Median.f90(18)
    REAL(KIND=mp)  AM                         ! SRC  Median.f90(19)
    REAL(KIND=mp)  SUM                        ! SRC  Median.f90(20)
    REAL(KIND=mp)  SUMX                       ! SRC  Median.f90(21)
    INTEGER        NP                         ! SRC  Median.f90(22)
    INTEGER        NM                         ! SRC  Median.f90(23)
    REAL(KIND=mp)  XP                         ! SRC  Median.f90(24)
    REAL(KIND=mp)  XM                         ! SRC  Median.f90(25)
    INTEGER        J                          ! SRC  Median.f90(26)
    REAL(KIND=mp)  XX                         ! SRC  Median.f90(27)
    REAL(KIND=mp)  DUM                        ! SRC  Median.f90(36)
    REAL(KIND=mp)  AA                         ! SRC  Median.f90(44)
  real X(N)
  real, parameter :: BIG = 1.e30, AFAC=1.5, AMP=1.5
! Here, AMP is an overconvergence factor: on each iteration,
! we move the guess by this factor. AFAC is a factor used to
! optimize the size of the "smoothing constant" EPS at each
! iteration.
  A=0.5*(X(1)+X(N))
  EPS=ABS(X(N)-X(1))
  AP=BIG
  AM=-BIG
1 SUM=0.
  SUMX=0.
  NP=0
  NM=0
  XP=BIG
  XM=-BIG
  do J=1, N
    XX=X(J)
    if(XX.ne.A)then
      if(XX.gt.A)then
        NP=NP+1
	if(XX.lt.XP) XP=XX
      else if(XX.lt.A)then
        NM=NM+1
	if(XX.gt.XM) XM=XX
      endif
      DUM=1./(EPS+ABS(XX-A))
      SUM=SUM+DUM
      SUMX=SUMX+XX*DUM
    endif
  end do
  print *,' NP=',NP,'  NM=',NM
  if(NP-NM.ge.2)then      !guess is too low, make another pass
    AM=A
    AA=XP+MAX(0.,SUMX/SUM-A)*AMP
    if(AA.gt.AP) AA=0.5*(A+AP)
    EPS=AFAC*ABS(AA-A)
    A=AA
    goto 1
  else if(NM-NP.ge.2)then !guess is too high
    AM=A
    AA=XM+MIN(0.,SUMX/SUM-A)*AMP
    if(AA.lt.AM) AA=0.5*(A+AM)
    EPS=AFAC*ABS(AA-A)
    A=AA
    goto 1
  else                    !guess is ok
    if(MOD(N,2).eq.0)then !for even N median is an average
      if(NP.eq.NM)then
        XMED=0.5*(XP+XM)
      else if(NP.gt.NM)then
        XMED=0.5*(A+XP)
      else
        XMED=0.5*(XM+A)
      endif
    else
      if(NP.eq.NM)then
        XMED=A
      else if(NP.gt.NM)then
        XMED=XP
      else
        XMED=XM
      endif
    endif
  endif
  return
END SUBROUTINE

!*******************************************************
!* Given an array X of N numbers, returns their median *
!* value XMED. The array X is modified and returned    *
!* sorted in ascending order.                          *
!*******************************************************
SUBROUTINE MDIAN(X,N,XMED)
    IMPLICIT NONE                             ! SRC
    INTEGER, PARAMETER :: mp = KIND(1.0D0)    ! SRC
    REAL(KIND=mp)  XMED                       ! SRC  Median.f90(83)
    INTEGER        N                          ! SRC  Median.f90(84)
    INTEGER        N2                         ! SRC  Median.f90(86)
  real X(N)
  call hpsort(N,X)
  N2=N/2
  if (2*N2.eq.N) then
    XMED = 0.5*(X(N2)+X(N2+1))
  else
    XMED = X(N2+1)
  endif
  return
END SUBROUTINE


!*****************************************************
!*  Sorts an array RA of length N in ascending order *
!*                by the Heapsort method             *
!* ------------------------------------------------- *
!* INPUTS:                                           *
!*	    N	  size of table RA                       *
!*      RA	  table to be sorted                     *
!* OUTPUT:                                           *
!*	    RA    table sorted in ascending order        *
!*                                                   *
!* NOTE: The Heapsort method is a N Log2 N routine,  *
!*       and can be used for very large arrays.      *
!*****************************************************         
SUBROUTINE HPSORT(N,RA)
    IMPLICIT NONE                             ! SRC
    INTEGER, PARAMETER :: mp = KIND(1.0D0)    ! SRC
    INTEGER        N                          ! SRC  Median.f90(110)
    INTEGER        L                          ! SRC  Median.f90(111)
    INTEGER        IR                         ! SRC  Median.f90(112)
    REAL(KIND=mp)  RRA                        ! SRC  Median.f90(120)
    INTEGER        I                          ! SRC  Median.f90(130)
    INTEGER        J                          ! SRC  Median.f90(131)
  real RA(N)
  L=N/2+1
  IR=N
  !The index L will be decremented from its initial value during the
  !"hiring" (heap creation) phase. Once it reaches 1, the index IR 
  !will be decremented from its initial value down to 1 during the
  !"retirement-and-promotion" (heap selection) phase.
10 continue
  if(L > 1)then
    L=L-1
    RRA=RA(L)
  else
    RRA=RA(IR)
    RA(IR)=RA(1)
    IR=IR-1
    if(IR.eq.1)then
      RA(1)=RRA
      return
    end if
  end if
  I=L
  J=L+L
20 if(J.le.IR)then
  if(J < IR)then
    if(RA(J) < RA(J+1))  J=J+1
  end if
  if(RRA < RA(J))then
    RA(I)=RA(J)
    I=J; J=J+J
  else
    J=IR+1
  end if
  goto 20
  end if
  RA(I)=RRA
  goto 10
END SUBROUTINE

END MODULE