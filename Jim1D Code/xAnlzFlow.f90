    Subroutine AnlzFlow(ic)
    !    curxsc(ns,1)=ws elevation
    !    curxsc(ns,2)=xsec area
    !    curxsc(ns,3)=hydraulic radius =area/wet perim
    !    curxsc(ns,4)=hydraulic depth =area/top width
    !    curxsc(ns,5)=top width
    !    curxsc(ns,6)=conveyance
    !	!    curxsc(ns,7)=discharge
    !  xsloc(i,5)=emin  ! xsloc( ,5) set to min ch elev
    use Support
    real(kind = mp2) kds
    call ChLims(ic,nus,nds)
    !	print *,nus,nds
    kds=curxsc(nds,6)**2
    ads=curxsc(nds,2)
    indxcrt=0
    do n=nds,nus+1,-1
        ns=n-1
        if(curxsc(n,2)<=0..or.curxsc(n,5)<=0.)then
            print *,'area error in AnlzFlow',n,curxsc(n,1),curxsc(n,2),curxsc(n,5)
            !				write(2,*)(curxsc(1,j),j=1,nsc)
            !				write(2,*)(curxsc(2,j),j=1,nsc)
            !				write(2,*)(curxsc(7,j),j=1,nsc)
            stop
        endif !(curxsc(n,2)<=0..or.curxsc(n,5)<=0.)then
        frd=abs(curxsc(n,7))/sqrt(gg*curxsc(n,2)**3/curxsc(n,5))
        q=curxsc(ns,7)
        call dxc(n,dx)
        sfdnom=q*abs(q)*2.
        atrm=2.*q*q/gg
        zdum=curxsc(ns,1)
        call zCrit(ns,abs(q),zc)
        fcrit=curxsc(n,1)-zc+dx*sfdnom/(kds+curxsc(ns,6)**2)+ &
            atrm*(1./ads-1./curxsc(ns,2))/(ads+curxsc(ns,2))
        call farea(ns,zdum)
        !			write(2,*)n,frd,fcrit,ads
        if(fcrit<0.)indxcrt(ns)=1
        kds=curxsc(ns,6)**2
        ads=curxsc(ns,2)
    enddo
    n=nds
    if(curxsc(n,2)<=0..or.curxsc(n,5)<=0.)then
        print *,'area error in AnlzFlow',n,curxsc(n,1),curxsc(n,2),curxsc(n,5)
        stop
    endif !(curxsc(n,2)<=0..or.curxsc(n,5)<=0.)then
    frd=abs(curxsc(nds,7))/sqrt(gg*curxsc(nds,2)**3/curxsc(nds,5))
    if(frd>=1.)indxcrt(nds)=1
    end Subroutine AnlzFlow !(ic)