      MODULE UInitMod
	USE RivVarMod
	USE CalcCond
	IMPLICIT NONE

	CONTAINS
	  subroutine uinit(u,v,hav,w,eta,q,ibc)
	  IMPLICIT NONE
        REAL, DIMENSION(ns, nn), INTENT(INOUT) :: u, v
	  REAL, DIMENSION(ns, nn), INTENT(IN) :: eta
        INTEGER, DIMENSION(ns, nn), INTENT(IN) :: ibc
	  REAL, DIMENSION(ns), INTENT(IN) ::w, hav	
	  REAL, INTENT(IN) :: q
	  INTEGER :: i, j, ibtot, count
	  REAL :: deln, depave, twidth, ubar, qpred
        do 30 i=1,ns
        deln=w(i)/(nn-1)
        depave=0.
        ibtot=0
		count = 0
        do 25 j=1,nn
        if(ibc(i,j).ne.0) then
         ibtot=ibtot+1
         depave=depave+GetDepth(i,j) 
        endif
25      continue
	  if(ibtot.eq.0) then
		errorcode = -10
          WRITE(6,*) 'At least one row of nodes are completely dry - '
          WRITE(6,*) 'check water surface boundary condition'
		return
	  end if
        depave=depave/(ibtot)
        twidth=deln*(ibtot)
        ubar=q/(twidth*depave)
        qpred=0.
	count = 0
        do 27 j=1,nn
        if(ibc(i,j).ne.0) then

c		  u(i,j)=ubar*(hav(i)-eta(i,j))/((j**5.)*depave)
		  if(vbc == 1 .and. i == 1) then
			u(i,j) = ubar*GetVelocityWieght(count,ibtot)
		  else
			u(i,j)=ubar*GetDepthWieght(i,j,depave)
		  endif
c		  u(i,j)=ubar
		  qpred=qpred+u(i,j)*GetDepth(i,j)*deln 
		  count = count+1
        else
		  u(i,j)=0.
        endif
c        if(i.eq.1) then
c        v(1,j)=-1.*u(1,j)*tand(30.)
c        else
        v(i,j)=0.
c        endif
27      continue
        count = 0
        do 28 j=1,nn
        u(i,j)=u(i,j)*q/qpred
	  if(i.eq.1) then
	      if(vbc == 1 .and. ibc(i,j).ne.0) then
	          v(1,j) = u(1,j)*tand(GetAngleWieght(count, ibtot))
	          count = count+1
	      else
		        v(1,j)=u(1,j)*tand(vac)
		    endif
        else
        v(i,j)=0.
         endif
28      continue
30 	continue
        return
        END SUBROUTINE uinit
        
	  subroutine VarDischUInit(u,v,w,eta,e,ibc, q, stageChange)
	  IMPLICIT NONE
        REAL, DIMENSION(ns, nn), INTENT(INOUT) :: u, v, e
	  REAL, DIMENSION(ns), INTENT(IN) ::w
	  REAL, DIMENSION(ns, nn), INTENT(IN) :: eta
        INTEGER, DIMENSION(ns, nn), INTENT(INOUT) :: ibc
	  REAL, INTENT(IN) :: stageChange	
	  REAL, INTENT(IN) :: q
	  
	  INTEGER :: i, j, ibtot, count
	  REAL :: deln, depave, twidth, ubar, qpred, depth
	  
        do 30 i=1,1
            deln=w(i)/(nn-1)
            depave=0.
            ibtot=0
		    count = 0
            do 25 j=1,nn
!                e(i,j) = e(i,j)
!                e(i,j) = e(i,j) + stageChange
                depth = e(i,j) - eta(i,j)
                if(depth > hmin) then
                    ibc(i,j) = -1
                else
                    ibc(i,j) = 0
                endif
                if(ibc(i,j).ne.0) then
                    ibtot=ibtot+1
                    depave=depave+depth
                endif
25          continue
	      if(ibtot.eq.0) then
              errorcode = -10
              WRITE(6,*) 'At least one row of nodes are completely dry-'
              WRITE(6,*) 'check water surface boundary condition'  
              return
	      end if
            depave=depave/(ibtot)
            twidth=deln*(ibtot)
            ubar=q/(twidth*depave)
            qpred=0.
	      count = 0
            do 27 j=1,nn
                if(ibc(i,j).ne.0) then
			        u(i,j)=ubar*((e(i,j)-eta(i,j))**vdc)/depave
!		            qpred=qpred+u(i,j)*GetDepth(i,j)*deln 
		            qpred=qpred+u(i,j)*(e(i,j) - eta(i,j))*deln 
		            count = count+1
                else
		            u(i,j)=0.
                endif
c        if(i.eq.1) then
c           v(1,j)=-1.*u(1,j)*tand(30.)
c        else
                 v(i,j)=0.
c        endif
27          continue
            do 28 j=1,nn
                u(i,j)=u(i,j)*q/qpred
	          if(i.eq.1) then
		            v(1,j)=-1.*u(1,j)*tand(vac)
                else
                    v(i,j)=0.
                endif
28          continue
30 	  continue
        return
        END SUBROUTINE VarDischUInit
        
	  subroutine VarDischUBInit(u,v,w,eta,e,ibc, q)
	  IMPLICIT NONE
        REAL, DIMENSION(ns, nn), INTENT(INOUT) :: u, v, e
	  REAL, DIMENSION(ns), INTENT(IN) ::w
	  REAL, DIMENSION(ns, nn), INTENT(IN) :: eta
        INTEGER, DIMENSION(ns, nn), INTENT(INOUT) :: ibc
	  REAL, INTENT(IN) :: q
	  
	  INTEGER :: i, j, ibtot, count
	  REAL :: deln, depave, twidth, ubar, qpred, depth
	  
        do 30 i=ns,ns
            deln=w(i)/(nn-1)
            depave=0.
            ibtot=0
		    count = 0
            do 25 j=1,nn
                depth = e(i,j) - eta(i,j)
                if(depth > 0.0) then
                    ibc(i,j) = -1
                else
                    ibc(i,j) = 0
                    u(i,j) = 0
                    v(i,j) = 0
                    hl(i,j) = 0
                endif
                if(ibc(i,j).ne.0) then
                    ibtot=ibtot+1
                    depave=depave+depth
                endif
25          continue
	      if(ibtot.eq.0) then
              errorcode = -10
              WRITE(6,*) 'At least one row of nodes are completely dry-'
              WRITE(6,*) 'check water surface boundary condition'  
              return
	      end if
            depave=depave/(ibtot)
            twidth=deln*(ibtot)
            ubar=q/(twidth*depave)
            qpred=0.
	      count = 0
            do 27 j=1,nn
                if(ibc(i,j).ne.0) then
			        u(i,j)=ubar*((e(i,j)-eta(i,j))**vdcds)/depave
!		            qpred=qpred+u(i,j)*GetDepth(i,j)*deln 
		            qpred=qpred+u(i,j)*(e(i,j) - eta(i,j))*deln 
		            count = count+1
                else
		            u(i,j)=0.
                endif
c        if(i.eq.1) then
c           v(1,j)=-1.*u(1,j)*tand(30.)
c        else
                 v(i,j)=0.
c        endif
27          continue
            do 28 j=1,nn
                u(i,j)=u(i,j)*q/qpred
	          if(i.eq.1) then
		            v(1,j)=-1.*u(1,j)*tand(vacds)
                else
                    v(i,j)=0.
                endif
28          continue
30 	  continue
        return
        END SUBROUTINE VarDischUBInit


	REAL FUNCTION GetDepth(i,j)
		INTEGER, INTENT(IN)::i, j
		GetDepth = hav(i)-eta(i,j)
	END FUNCTION

	REAL FUNCTION GetDepthWieght(i, j, depave)
		INTEGER, INTENT(IN):: i, j
		REAL, INTENT(IN)::depave
		GetDepthWieght = ((hav(i)-eta(i,j))**vdc)/depave
	END FUNCTION

	REAL FUNCTION GetVelocityWieght(i,nn)
		INTEGER, INTENT(IN) :: i,nn
		INTEGER :: count
		REAL :: val
		val = real(i)/real(nn)
		count = 1
		DO WHILE(vbcdist(count) <= val)
			count=count+1
		END DO
		count = count -1
		GetVelocityWieght = vbcvel(count) +			!y
     &	 (val-vbcdist(count)) *		!x-xo
     &	 (((vbcvel(count+1)-vbcvel(count))/			!slope
     &      (vbcdist(count+1)-vbcdist(count))))
	END FUNCTION
	
	REAL FUNCTION GetAngleWieght(i,nn)
	INTEGER, INTENT(IN) :: i, nn
	INTEGER :: count
	REAL :: val
	val = real(i)/real(nn)
	count = count - 1
	DO WHILE (vbcdist(count) <= val)
			count=count+1
	END DO
		count = count -1
		GetAngleWieght = vbcang(count) +			!y
     &	 (val-vbcdist(count)) *		!x-xo
     &	 (((vbcang(count+1)-vbcang(count))/			!slope
     &      (vbcdist(count+1)-vbcdist(count))))
      END FUNCTION

	END MODULE UInitMod
