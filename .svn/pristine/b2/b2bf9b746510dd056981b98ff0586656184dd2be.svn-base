	Subroutine tsbc(nc,nloc,t,q)
	use Support
! times t are in days
	ns=tsbcloc(nc,2*nloc-1)
	ne=tsbcloc(nc,2*nloc)
! 	print *,ns,ne,t
	if(t<=tmsrbc(ns,1))then
		q=tmsrbc(ns,2)
		return
	elseif(t>=tmsrbc(ne,1))then
		q=tmsrbc(ne,2)
		return
	else
		do n=ns+1,ne
			if(t<=tmsrbc(n,1)) then
				r=(t-tmsrbc(n-1,1))/(tmsrbc(n,1)-tmsrbc(n-1,1))
				q=tmsrbc(n-1,2)*(1.-r)+tmsrbc(n,2)*r
				return
			endif
		enddo
	endif
	end Subroutine tsbc
