 	MODULE EInitMod
	USE CalcCond
	USE RivVarMod
	IMPLICIT NONE

CONTAINS 
	SUBROUTINE EInit(e, hl, eta, ibc, w, hav)
	IMPLICIT NONE
	REAL, DIMENSION(ns, nn), INTENT(INOUT) :: e, hl
	REAL, DIMENSION(ns, nn), INTENT(IN) :: eta
	INTEGER, DIMENSION(ns, nn), INTENT(IN) :: ibc
	REAL, DIMENSION(ns), INTENT(IN) :: hav, w
	REAL :: deln, twidth, dz
	INTEGER :: i, j, ibtot, iinc
		write(6,*) "WSE", lbwse, rbwse
	DO 30 i = 1,ns
		deln = w(i)/nn-1
		ibtot = 0
		DO 25 j = 1,nn	
			IF(ibc(i,j).ne.0) THEN
				ibtot = ibtot+1
			ENDIF
25		CONTINUE
		twidth = deln*ibtot
		dz = (lbwse - rbwse)/twidth
		iinc = 0
		DO 27 j = 1,nn
!			IF(i.eq.ns.and.stagetype.eq.1) THEN
!				IF(ibc(i,j).ne.0) then
!					e(i,j) = rbwse + (iinc*deln)*dz
!					iinc = iinc + 1
!				ELSE
!					IF(j.lt.nn/2) THEN
!						e(i,j) = rbwse
!					ELSE
!						e(i,j) = lbwse
!					ENDIF
!				ENDIF
!			ELSE
				e(i,j) = hav(i)
!			ENDIF
			hl(i,j) = e(i,j) - eta(i,j)
			IF(hl(i,j).lt.hmin) THEN 
				hl(i,j) = hmin
			ENDIF
27		CONTINUE
30	CONTINUE
	RETURN
	END SUBROUTINE EInit
		
	
	SUBROUTINE VarDischEInit(e, hl, u, v, eta, ibc, stageChange, dsstage)
	IMPLICIT NONE
	    REAL, DIMENSION(ns, nn), INTENT(INOUT) :: e, hl, u, v
	    INTEGER, DIMENSION(ns, nn), INTENT(INOUT) :: ibc
	    REAL, DIMENSION(ns, nn), INTENT(INOUT) :: eta
	    REAL, INTENT(IN) :: stageChange
	    REAL, INTENT(INOUT) :: dsstage
	    REAL :: depth
	    INTEGER :: i, j

		!First adjust stage for all previously wet nodes
		DO i = 1,ns
		    DO j = 1, nn
		        if(ibc(i,j).ne.0) then
		            e(i,j) = e(i,j) + stageChange
		            if(i.eq.ns) then
		                dsstage = e(i,j)
		            endif
		            depth = e(i,j) - eta(i,j)
		            if(depth > hmin) then
		                hl(i,j) = e(i,j) - eta(i,j)
!		                if(dryType.eq.1.and.hl(i,j).lt.hmin) then
!		                    hl(i,j) = hmin
!		                endif
		                ibc(i,j) = -1
		                
		            else
!		                e(i,j) = eta(i,j)
		                hl(i,j) = hmin
		                ibc(i,j) = 0
		                u(i,j) = 0.
		                v(i,j) = 0.
		            endif
		        else
		            e(i,j) = eta(i,j)
		        endif
		      
		    enddo
		enddo
!		DO j = 1,nn
!            e(1,j) = e(1,j) + stageChange
!            depth = e(1,j) - eta(1,j)
!            if(depth > hmin) then
!                hl(1,j) = e(1,j) - eta(1,j)
!!		                if(dryType.eq.1.and.hl(i,j).lt.hmin) then
!!		                    hl(i,j) = hmin
!!		                endif
!                ibc(1,j) = -1
!                
!            else
!!		                e(i,j) = eta(i,j)
!                hl(1,j) = hmin
!                ibc(1,j) = 0
!                u(1,j) = 0.
!                v(1,j) = 0.
!            endif		
!        END DO
	    RETURN
	END SUBROUTINE VarDischEInit
	
	SUBROUTINE UpdateWETTING(e, hl, u, v, eta, ibc, dsstage)
	    IMPLICIT NONE
	    REAL, DIMENSION(ns, nn), INTENT(INOUT) :: e, hl,u,v
	    INTEGER, DIMENSION(ns, nn), INTENT(INOUT) :: ibc
	    REAL, DIMENSION(ns, nn), INTENT(IN) :: eta
	    REAL, INTENT(in) :: dsstage
	    INTEGER :: I, J, IM, JP, JM, IP 
	    INTEGER :: change = -1
	    change = -1
	    do while (change == -1)  !CURRENTLY FORCED TO PASS ONLY ONCE THROUGH LOOP
	        change = 1
	        do  i=ns,1,-1
                do  j=2,nn-1
                    IF(i.EQ.1) THEN
                        im=i
                    ELSE
                        im=i-1
                    END IF
                    IF(i.EQ.ns) THEN
                        ip=ns
                    ELSE
                        ip=i+1
                    ENDIF
                    jp=j+1
                    jm=j-1

		            if (ibc(i,j).ne.0) then !added by rmcd to account for ibc == 4 or 6
                        if(ibc(im,j).eq.0.and.e(i,j).gt.eta(im,j)+hwmin.and.i.gt.1)then              
                            hl(im,j)=e(i,j)-eta(im,j)
                            e(im,j)=e(i,j)
                            ibc(im,j)=-1
                            u(im,j) = 0.
                            v(im,j) = 0.
                            change = -1
!                        endif
                        elseif(ibc(im,jp).eq.0.and.e(i,j).gt.eta(im,jp)+hwmin.and.i.gt.1)then      
 		                    if(jp < nn) then
                           hl(im,jp)=e(i,j)-eta(im,jp)
                            e(im,jp)=e(i,j)
		                    	ibc(im,jp)=-1
		                    	u(im,jp) = 0.
		                    	v(im,jp) = 0.
		                    	change = -1
		                    endif
!                        endif
                        elseif(ibc(im,jm).eq.0.and.e(i,j).gt.eta(im,jm)+hwmin.and.i.gt.1)then      
 		                    if(jm > 1)	then
                           hl(im,jm)=e(i,j)-eta(im,jm)
                            e(im,jm)=e(i,j)
		                        ibc(im,jm)=-1
		                        u(im,jm) = 0.
		                        v(im,jm) = 0.
		                        change = -1
		                    endif
!                        endif
                        elseif(ibc(i,jp).eq.0.and.e(i,j).gt.eta(i,jp)+hwmin.and.i.le.ns.and.i.gt.1)then  
		                    if(jp < nn)	then
                            hl(i,jp)=e(i,j)-eta(i,jp)                                
                            e(i,jp)=e(i,j)
		                        ibc(i,jp)=-1
		                        u(i,jp) = 0.
		                        v(i,jp) = 0.
		                        change = -1
		                    endif
!                        endif
                        elseif(ibc(i,jm).eq.0.and.e(i,j).gt.eta(i,jm)+hwmin.and.i.le.ns.and.i.gt.1)then  
		                    if(jm > 1)	then
                            hl(i,jm)=e(i,j)-eta(i,jm)
                            e(i,jm)=e(i,j)
		                        ibc(i,jm)=-1
		                        u(i,jm) = 0.
		                        v(i,jm) = 0.
		                        change = -1
		                    endif
!                        endif
                        elseif(ibc(ip,jp).eq.0.and.e(i,j).gt.eta(ip,jp)+hwmin.AND.i .le. ns-2)then     
		                    if(jp < nn)	then
                            hl(ip,jp)=(e(i,j)-eta(ip,jp))                         
!                            e(ip,jp)=hl(ip,jp)+eta(ip,jp)
!                            if(i.eq.ns-1) then
!                                e(ip,jp) = dsstage
!                            else
                                e(ip,jp)=e(i,j)
!                            endif
		                        ibc(ip,jp)=-1
		                        u(ip,jp) = 0.
		                        v(ip,jp) = 0.
		                        change = -1
                            endif		                        
!                        endif
                        elseif(ibc(ip,j).eq.0.and.e(i,j).gt.eta(ip,j)+hwmin .AND.i .le. ns-2)then                
                            hl(ip,j)=(e(i,j)-eta(ip,j))
!                            e(ip,j)=hl(ip,j)+eta(ip,j)
!                            e(ip,j)=e(i,j)
!                            if(i.eq.ns-1) then
!                                e(ip,j) = dsstage
!                            else
                                e(ip,j)=e(i,j)
!                            endif
                            ibc(ip,j)=-1
		                        u(ip,j) = 0.
		                        v(ip,j) = 0.
		                        change = -1
                            
!                        endif
                        elseif(ibc(ip,jm).eq.0.and.e(i,j).gt.eta(ip,jm)+hwmin .and.i .le. ns-2)then 
                           if(jm > 1)	then
                            hl(ip,jm)=(e(i,j)-eta(ip,jm))                         
 !                           e(ip,jm)=hl(ip,jm)+eta(ip,jm)
!                            e(ip,jm)= e(i,j)
!                            if(i.eq.ns-1) then
!                                e(ip,jm) = dsstage
!                            else
                                e(ip,jm)=e(i,j)
!                            endif
                                ibc(ip,jm)=-1
		                        u(ip,jm) = 0.
		                        v(ip,jm) = 0.
		                        change = -1
                            endif
                                                          
                        endif 
                  endif  
                ENDDO
            ENDDO
            change = 1
        ENDDO
    END SUBROUTINE UpdateWETTING
    
    SUBROUTINE CheckNodeContinuity(e, hl, u, v, eta, ibc)
	    IMPLICIT NONE
	    REAL, DIMENSION(ns, nn), INTENT(INOUT) :: e, hl,u,v
	    INTEGER, DIMENSION(ns, nn), INTENT(INOUT) :: ibc
	    REAL, DIMENSION(ns, nn), INTENT(IN) :: eta
	    INTEGER :: I, J, IM, JP, JM, IP 
	    INTEGER :: change = -1
	    INTEGER :: Count = 0
	    change = -1
	    do while (change == -1)
	        change = 1
!	        do  i=ns,1,-1
 	        do  i=1,ns
               do  j=2,nn-1
                    IF(i.EQ.1) THEN
                        im=i
                    ELSE
                        im=i-1
                    END IF
                    IF(i.EQ.ns) THEN
                        ip=ns
                    ELSE
                        ip=i+1
                    ENDIF
                    jp=j+1
                    jm=j-1
		            if (ibc(i,j).ne.0) then !added by rmcd to account for ibc == 4 or 6
		                Count = 0
		                if(i.eq.ns) then
		                     if(ibc(im,j).ne.0) then
		                        Count = Count+3
		                     endif
		                     if(ibc(i,jp).ne.0) then 
		                        Count = Count+1
		                     endif
		                     if(ibc(i,jm).ne.0) then
		                        Count = Count+1
		                     endif
		                     if(Count < 4) then
		                        u(i,j) = 0.
		                        v(i,j) = 0.
		                        ibc(i,j) = 0
!		                        change = -1
                             endif
                        else if(i.eq.1) then
                             if(ibc(ip,j).ne.0) then
		                        Count = Count+3
		                     endif
		                     if(ibc(i,jp).ne.0) then 
		                        Count = Count+1
		                     endif
		                     if(ibc(i,jm).ne.0) then
		                        Count = Count+1
		                     endif
		                     if(Count < 4) then
		                        u(i,j) = 0.
		                        v(i,j) = 0.
		                        ibc(i,j) = 0
!		                        change = -1
                             endif		                
		                else 
		                    if(ibc(im,j).ne.0)then
		                        Count = Count+1
		                    endif
    		                if(ibc(ip,j).ne.0) then
		                        Count = Count+1
		                    endif
     		                if(ibc(i,jm).ne.0) then
		                        Count = Count+1
		                    endif
		                    if(ibc(i,jp).ne.0) then
		                        Count = Count+1
		                    endif
		                    if(count <= 1) then
		                        u(i,j) = 0.
		                        v(i,j) = 0.
		                        ibc(i,j) = 0
!		                        change = -1
		                    endif
		                endif
		            endif
		        enddo
		    enddo
        enddo
    END SUBROUTINE CheckNodeContinuity

END MODULE 
