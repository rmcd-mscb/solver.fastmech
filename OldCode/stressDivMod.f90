    MODULE StressDivMod
    use RivVarMod2
    USE CalcCond2
    IMPLICIT NONE
    REAL, ALLOCATABLE, DIMENSION(:,:) :: DUM1
    CONTAINS
    SUBROUTINE StressDiv(ibc, qs, qn, taus, taun, con, rn, r)
    REAL ::  SC, weight
    INTEGER :: I, J, ISMOO
    CALL alloc_StressDiv()
    allocate(dum1(
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

    DO ISMOO=1,sedsmoo
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
    CALL dealloc_StressDiv
    end subroutine

    SUBROUTINE alloc_StressDiv()
    INTEGER :: status
    !Allocate REAL types
    !			ALLOCATE(QS(ns, nn), STAT = status)
    !			ALLOCATE(QN(ns, nn), STAT = status)
    ALLOCATE(DUM1(ns, nn), STAT = status)
    END SUBROUTINE alloc_StressDiv

    SUBROUTINE dealloc_StressDiv()
    INTEGER :: status
    !			DEALLOCATE(QS, STAT = status)
    !			DEALLOCATE(QN, STAT = status)
    DEALLOCATE(DUM1, STAT = status)
    END SUBROUTINE dealloc_StressDiv


    END MODULE
