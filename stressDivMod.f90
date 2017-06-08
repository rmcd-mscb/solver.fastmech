    MODULE StressDivMod
	use RivVarMod
	USE CalcCond
	IMPLICIT NONE
	REAL, ALLOCATABLE, DIMENSION(:,:) :: DUM1
	CONTAINS
		SUBROUTINE StressDiv()
			REAL ::  SC, weight
			INTEGER :: I, J, ISMOO
			CALL alloc_StressDiv()

			DO 615 I=1,ns
			DO 615 J=1,nn
			
			IF (ibc(i,j).eq.0) THEN
			 QS(I,J)=0.
			 QN(I,J)=0.
			 GO TO 615
			ENDIF
			QS(I,J)=taus(I,J)
			QN(I,J)=taun(I,J)
			IF(J.EQ.1.OR.J.EQ.nn) THEN
			 QN(I,J)=0.
			ENDIF
615		 CONTINUE
			DO 650 I=1,ns
			DO 650 J=1,nn
			SC=rn(I,J)
			IF(I.EQ.1.) THEN
			    CON(I,J)=(QS(2,J)-QS(1,J))/(ds*SC)
			    GO TO 645
			ELSEIF(I.EQ.ns) THEN
			    CON(I,J)=CON(I-1,J)
			    GO TO 643
			ELSE
			    CON(I,J)=(QS(I+1,J)-QS(I-1,J))/(2.*ds*SC)
			ENDIF
			
643		    CON(I,J)=CON(I,J)-(QN(I,J)/(R(I)*SC))

			IF(J.EQ.1) THEN
			   CON(I,1)=CON(I,1)+(QN(I,2)-QN(I,1))/dn
				GO TO 645
			ELSE
			    IF(J.EQ.nn) THEN
			        CON(I,nn)=CON(I,nn)+(QN(I,nn)-QN(I,nn-1))/dn
			        GO TO 645
			    ENDIF
			    
			    CON(I,J)=CON(I,J)+(QN(I,J+1)-QN(I,J-1))/(2.*dn)
			    
			ENDIF
645			CONTINUE
650			CONTINUE

			DO 670 ISMOO=1,sedsmoo
			DO 660 I=2,ns-1
			DO 660 J=1,nn
			if(ibc(i,j).eq.0) then
			 con(i,j)=0.
			 go to 660
			endif
			if(j.eq.1) then
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
660			continue
			DO 670 I=2,ns-1
			DO 670 J=1,nn
670			CON(I,J)=DUM1(I,J)
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
