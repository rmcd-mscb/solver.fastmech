	Module TSplineMod
	CONTAINS
	  subroutine tspline(x,y,n,xout,yout,iout,sigma,yp,temp)
		  dimension x(n),yp(n),temp(n),xout(iout),yout(iout)
		  double precision y(n)
		  nm1=n-1
		  np1=n+1
		  delx1=x(2)-x(1)
		  dx1=(y(2)-y(1))/delx1
		  delx2=x(3)-x(2)
		  delx12=x(3)-x(1)
		  c1=-(delx12+delx1)/delx12/delx1
		  c2=delx12/delx1/delx2
		  c3=-delx1/delx12/delx2
		  slpp1=c1*y(1)+c2*y(2)+c3*y(3)
		  deln=x(n)-x(nm1)
		  delnm1=x(nm1)-x(n-2)
		  delnn=x(n)-x(n-2)
		  c1=(delnn+deln)/delnn/deln
		  c2=-delnn/deln/delnm1
		  c3=deln/delnn/delnm1
		  slppn=c3*y(n-2)+c2*y(nm1)+c1*y(n)
10		  sigmap=abs(sigma)*float(n-1)/(x(n)-x(1))
		  dels=sigmap*delx1
		  exps=exp(dels)
		  sinhs=.5*(exps-1./exps)
		  sinhin=1./(delx1*sinhs)
		  diag1=sinhin*(dels*.5*(exps+1./exps)-sinhs)
		  diagin=1./diag1
		  yp(1)=diagin*(dx1-slpp1)
		  spdiag=sinhin*(sinhs-dels)
		  temp(1)=diagin*spdiag
		  do 20 i=2,nm1
			   delx2=x(i+1)-x(i)
			   dx2=(y(i+1)-y(i))/delx2
			   dels=sigmap*delx2
			   exps=exp(dels)
			   sinhs=.5*(exps-1./exps)
			   sinhin=1./(delx2*sinhs)
			   diag2=sinhin*(dels*(.5*(exps+1./exps))-sinhs)
			   diagin=1./(diag1+diag2-spdiag*temp(i-1))
			   yp(i)=diagin*(dx2-dx1-spdiag*yp(i-1))
			   spdiag=sinhin*(sinhs-dels)
			   temp(i)=diagin*spdiag
			   dx1=dx2
			   diag1=diag2
20		  continue
30		  diagin=1./(diag1-spdiag*temp(nm1))
		  yp(n)=diagin*(slppn-dx2-spdiag*yp(nm1))
		  do 40 i=2,n
				ibak=np1-i
				yp(ibak)=yp(ibak)-temp(ibak)*yp(ibak+1)
40		  continue
		  a=x(1)
		  b=x(2)
		  nj=2
		  do 150 i=1,iout
145			  if(xout(i).gt.b) then
				   a=b
				   nj=nj+1
				   b=x(nj)
				   go to 145
			  endif
			  del1=xout(i)-a
			  del2=b-xout(i)
			  dels=b-a
			  exps1=exp(sigmap*del1)
			  sinhd1=.5*(exps1-1./exps1)
			  exps=exp(sigmap*del2)
			  sinhd2=.5*(exps-1./exps)
			  exps=exps*exps1
			  sinhs=.5*(exps-1./exps)
150			  yout(i)=(yp(nj)*sinhd1+yp(nj-1)*sinhd2)/sinhs+
     &			  ((y(nj)-yp(nj))*del1+(y(nj-1)-yp(nj-1))*del2)/dels
		return
        end subroutine
	end module
