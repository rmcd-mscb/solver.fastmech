    MODULE StressDivMod2
    use RivVarMod2
    USE CalcCond2
    IMPLICIT NONE
    !REAL, ALLOCATABLE, DIMENSION(:,:) :: DUM1
    CONTAINS
    SUBROUTINE StressDiv(rvo, cco, ibc, qs, qn, taus, taun, con, rn, r)
    implicit none
    type(rivvar), intent(in) :: rvo
    type(calccond), intent(in) :: cco
    REAL(kind=mp), ALLOCATABLE, DIMENSION(:,:) :: DUM1
    real(kind=mp), dimension(:,:), intent(in) :: taus, taun, rn
    real(kind=mp), dimension(:), intent(in) :: r
    integer, dimension(:,:), intent(in) :: ibc
    real(kind=mp), dimension(:,:), intent(inout) :: con, qs, qn
    REAL(kind=mp) ::  SC, weight, ds, dn
    INTEGER :: I, J, ISMOO
    integer :: ns, nn, ier
    ns = rvo%ns
    nn = rvo%nn
    ds = rvo%ds
    dn = rvo%dn
    allocate(dum1(ns, nn), stat=ier)
    DO I=1,ns
        DO J=1,nn

            IF (ibc(i,j).eq.0) THEN
                QS(I,J)=0.
                QN(I,J)=0.
            ELSE
                QS(I,J)=taus(I,J)
                QN(I,J)=taun(I,J)
                IF(J.EQ.1.OR.J.EQ.nn) THEN
                    QN(I,J)=0.
                ENDIF
            ENDIF
        ENDDO
    ENDDO
    DO I=1,ns
        DO J=1,nn
            SC=rn(I,J)
            IF(I.EQ.1.) THEN
                CON(I,J)=(QS(2,J)-QS(1,J))/(ds*SC)
                CYCLE
            ELSEIF(I.EQ.ns) THEN
                CON(I,J)=CON(I-1,J)
            ELSE
                CON(I,J)=(QS(I+1,J)-QS(I-1,J))/(2.*ds*SC)
            ENDIF

            CON(I,J)=CON(I,J)-(QN(I,J)/(R(I)*SC))

            IF(J.EQ.1) THEN
                CON(I,1)=CON(I,1)+(QN(I,2)-QN(I,1))/dn
                CYCLE
            ELSE
                IF(J.EQ.nn) THEN
                    CON(I,nn)=CON(I,nn)+(QN(I,nn)-QN(I,nn-1))/dn
                ELSE
                    CON(I,J)=CON(I,J)+(QN(I,J+1)-QN(I,J-1))/(2.*dn)
                ENDIF
            ENDIF
        ENDDO
    ENDDO

    DO ISMOO=1,cco%sedsmoo
        DO I=2,ns-1
            DO  J=1,nn
                if(ibc(i,j).eq.0) then
                    con(i,j)=0.

                elseif(j.eq.1) then
                    DUM1(I,1)=(CON(I-1,1)+CON(I,1)+CON(I+1,1)+con(I,j+1))/4.
                elseif(j.eq.nn) then
                    DUM1(I,nn)=(CON(I-1,nn)+CON(I,nn)+CON(I+1,nn)+con(i,nn-1))/4.
                else
                    weight=1
                    if(con(i,j-1).ne.0) then
                        weight=weight+0.5
                    endif
                    if(con(i,j+1).ne.0) then
                        weight=weight+0.5
                    endif
                    if(con(i-1,j).ne.0) then
                        weight=weight+0.5
                    endif
                    if(con(i+1,j).ne.0) then
                        weight=weight+0.5
                    endif
                    DUM1(I,J)=1.*CON(I,J)+.5*CON(I+1,J)+.5*CON(I-1,J)
                    DUM1(I,J)=DUM1(I,J)+.5*CON(I,J-1)+.5*CON(I,J+1)
                    DUM1(I,J)=DUM1(I,J)/weight
                endif
            ENDDO
        ENDDO


        DO I=2,ns-1
            DO J=1,nn
                CON(I,J)=DUM1(I,J)
            ENDDO
        ENDDO
    ENDDO
    deallocate(dum1, stat=ier)
    end subroutine

    END MODULE
