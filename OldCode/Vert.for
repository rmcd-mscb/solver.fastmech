        subroutine vert
	  use RivVarMod
!        parameter (ns=41,nn=25,nz=11)
!        real mo,hav
!        common taus(ns,nn),taun(ns,nn),hl(ns,nn),eta(ns,nn),rn(ns,nn),
!     &  ibc(ns,nn),r(ns),w(ns),cd,ds,dn,mo,q,hwt,urelax,erelax,
!     &  arelax,itm,hav,iplinc,xo(ns),yo(ns),uz(ns,nn,nz),vz(ns,nn,nz)
        dimension ecoef(ns,nn,nz),zeta0(ns,nn),f1(ns,nn,nz),zeta(nz),
     &  fone(ns,nn),f2(ns,nn,nz),ftwo(ns,nn),f3(ns,nn,nz),
     &  dzeta(nz),theta(ns,nn),ustr2(ns,nn),rs(ns,nn)
        pi=acos(-1.)
        vkc=0.4
        beta=6.24
        do 10 i=1,ns
        do 10 j=1,nn
        zeta0(i,j)=.007
        sczeta=(alog(1./zeta0(i,j)))/(nz-1)
        do 4 k=1,nz
        zeta(k)=zeta0(i,j)*exp((k-1)*sczeta)
        if(zeta(k).lt.0.2) then
         ecoef(i,j,k)=vkc*zeta(k)*(1.-zeta(k))
        else
         ecoef(i,j,k)=vkc/beta
        endif
        if(k.eq.1) then
         f1(i,j,k)=0.
         fone(i,j)=0.
         f3(i,j,k)=0.
        else
         dzeta(k)=zeta(k)-zeta(k-1)
         finc1=(1.-zeta(k))/ecoef(i,j,k)
         finc2=(1.-zeta(k-1))/ecoef(i,j,k-1)
         finc=(finc1+finc2)
         dzp=0.5*dzeta(k)
         f1(i,j,k)=f1(i,j,k-1)+finc*dzp
         fone(i,j)=fone(i,j)+(f1(i,j,k)+f1(i,j,k-1))*dzp
         f3(i,j,k)=f3(i,j,k-1)+(f1(i,j,k)**2.+f1(i,j,k-1)**2.)*dzp
        endif
4       continue
        do 6 k=1,nz
        f3(i,j,k)=f3(i,j,nz)-f3(i,j,k)
        if(k.eq.1) then
         f2(i,j,k)=0.
         ftwo(i,j)=0.
        else
         finc1=f3(i,j,k)/ecoef(i,j,k)
         finc2=f3(i,j,k-1)/ecoef(i,j,k-1)
         f2(i,j,k)=f2(i,j,k-1)+(finc1+finc2)*0.5*dzeta(k)
         ftwo(i,j)=ftwo(i,j)+(f2(i,j,k)+f2(i,j,k-1))*0.5*dzeta(k)
        endif 
6       continue       
        if(taus(i,j).eq.0.and.taun(i,j).eq.0) then
         theta(i,j)=1000.
        else
         theta(i,j)=atan2(taun(i,j),taus(i,j))
        endif
        ustr2(i,j)=(taus(i,j)**2.+taun(i,j)**2.)**0.5
10      continue
        do 50 i=1,ns
        do 50 j=1,nn
        if(i.eq.1.or.theta(i-1,j).eq.1000.) then
         dth=theta(i+1,j)-theta(i,j)
         scn=1.
        elseif(i.eq.ns.or.theta(i+1,j).eq.1000.) then
         dth=theta(i,j)-theta(i-1,j)
         scn=1.
        else
         dth=(theta(i+1,j)-theta(i-1,j))
         scn=2.
        endif
        if(dth.gt.pi) then
         dth=dth-2.*pi
        elseif(dth.lt.(-1.*pi)) then
         dth=2.*pi+dth
        endif
        dthds=dth/(scn*ds*rn(i,j))
        if(j.eq.1.or.theta(i,j-1).eq.1000.) then
         dth=theta(i,j+1)-theta(i,j)
         scn=1.
        elseif(j.eq.nn.or.theta(i,j+1).eq.1000.) then
         dth=theta(i,j)-theta(i,j-1)
         scn=1.
        else
         dth=(theta(i,j+1)-theta(i,j-1))
         scn=2.
        endif
        if(dth.gt.pi) then
         dth=dth-2.*pi
        elseif(dth.lt.(-1.*pi)) then
         dth=2.*pi+dth
        endif
        dthdn=dth/(scn*dn)
        if(ustr2(i,j).eq.0.or.ibc(i,j).ne.-1) then 
         rs(i,j)=10**8.
        else
         curvs=(dthds+(1./r(i)))*taus(i,j)/ustr2(i,j)
         curvn=dthdn*taun(i,j)/ustr2(i,j)
         curv=curvs+curvn
         if(curv.eq.0) then
          rs(i,j)=10**8.
         else 
          rs(i,j)=1./curv
         endif
        endif
        if(abs(rs(i,j)).lt.5000.) then
         rs(i,j)=5000.*sign(1.,rs(i,j))
        endif
        taux=ustr2(i,j)*hl(i,j)/(rn(i,j)*rs(i,j)) 
        vs=((ftwo(i,j)/fone(i,j))-f3(i,j,1))
        taux=taux*vs
        do 40 k=1,nz
        uu1=(ustr2(i,j)**.5)*cos(theta(i,j))*f1(i,j,k)
        vv1=(ustr2(i,j)**.5)*sin(theta(i,j))*f1(i,j,k)
        g1=(ustr2(i,j)**.5)*hl(i,j)/(rn(i,j)*rs(i,j))
        g2=((ftwo(i,j)/fone(i,j))*f1(i,j,k))-f2(i,j,k)
        uu2=-1.*g1*g2*sin(theta(i,j))
        vv2=g1*g2*cos(theta(i,j)) 
        uz(i,j,k)=uu1+uu2
40      vz(i,j,k)=vv1+vv2
        taus(i,j)=taus(i,j)-taux*sin(theta(i,j))
        taun(i,j)=taun(i,j)+taux*cos(theta(i,j))
50      continue
c      DO 670 ISMOO=1,2
c      DO 660 I=2,ns-1
c      DO 660 J=1,nn
c      if(j.eq.1) then
c       DUM1(I,1)=(rs(I-1,1)+rs(I,1)+rs(I+1,1)+rs(I,j+1))/4.
c      elseif(j.eq.nn) then
c       DUM1(I,25)=(rs(I-1,nn)+rs(I,nn)+rs(I+1,nn)+rs(i,nn-1))/4.
c      else
c       DUM1(I,J)=rs(I,J)+rs(I+1,J)+rs(I-1,J)
c       DUM1(I,J)=DUM1(I,J)+rs(I,J-1)+rs(I,J+1)
c       DUM1(I,J)=DUM1(I,J)/5.
c      endif
c660   continue
c      DO 670 I=2,ns-1
c      DO 670 J=1,nn
c670    rs(I,J)=DUM1(I,J)
c        do 750 i=1,ns
c        do 750 j=1,nn
c        if(abs(rs(i,j)).lt.3000.) then
c         rs(i,j)=3000.*sign(1.,rs(i,j))
c        endif
c        taux=ustr2(i,j)*hl(i,j)/(rn(i,j)*rs(i,j)) 
c        vs=((ftwo(i,j)/fone(i,j))-f3(i,j,1))
c        taux=taux*vs
c        taus(i,j)=taus(i,j)-taux*sin(theta(i,j))
c750     taun(i,j)=taun(i,j)+taux*cos(theta(i,j))
      return
      end
