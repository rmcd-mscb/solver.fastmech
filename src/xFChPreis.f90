!    Subroutine CompOnly(ti,dt)
!    use support
!    common/RstDbg/rstrtset
!    !    curxsc(ns,1)=ws elevation
!    !    curxsc(ns,2)=xsec area
!    !    curxsc(ns,3)=hydraulic radius =area/wet perim
!    !    curxsc(ns,4)=hydraulic depth =area/top width
!    !    curxsc(ns,5)=top width
!    !    curxsc(ns,6)=conveyance
!    !	!    curxsc(ns,7)=discharge
!    !  xsloc(i,5)=emin  ! xsloc( ,5) set to min ch elev
!    !dimension rdum(6),cmdum(6,6)
!    real kb
!    ! Continuity functions
!    bbar(j)=.5*(thet*(curxsc(j,5)+curxsc(j-1,5))+(1.-thet)*(oldxsc(j,5)+oldxsc(j-1,5)))
!    pqpx(j)=(thet*(curxsc(j,7)-curxsc(j-1,7))+(1.-thet)*(oldxsc(j,7)-oldxsc(j-1,7)))/dx
!    pypt(j)=.5*(curxsc(j-1,1)-oldxsc(j-1,1)+curxsc(j,1)-oldxsc(j,1))/dt
!    ! Momentum functions
!    abar(j)=.5*(thet*(curxsc(j,2)+curxsc(j-1,2))+(1.-thet)*(oldxsc(j,2)+oldxsc(j-1,2)))
!    qabsqbar(j)=.5*(thet*(curxsc(j,7)*abs(curxsc(j,7))+curxsc(j-1,7)*abs(curxsc(j-1,7)))&
!        +(1.-thet)*(oldxsc(j,7)*abs(oldxsc(j,7))+oldxsc(j-1,7)*abs(oldxsc(j-1,7))))
!    pypx(j)=(thet*(curxsc(j,1)-curxsc(j-1,1))+(1.-thet)*(oldxsc(j,1)-oldxsc(j-1,1)))/dx
!    pqsovapx(j)=(thet*(curxsc(j,7)**2/curxsc(j,2)-curxsc(j-1,7)**2/curxsc(j-1,2))+&
!        (1.-thet)*(oldxsc(j,7)**2/oldxsc(j,2)-oldxsc(j-1,7)**2/oldxsc(j-1,2)))/dx
!    pqpt(j)=.5*(curxsc(j-1,7)-oldxsc(j-1,7)+curxsc(j,7)-oldxsc(j,7))/dt
!    call ChLims(nch,nus,nds)
!    do ns=1,nds
!        do j=1,7
!            curxsc(ns,j)=oldxsc(ns,j)
!            if(j/=7)plnxsc(ns,j)=oplnxsc(ns,j)
!        enddo
!    enddo
!    dtinv=.5/dt
!    do iter=1,itersmx
!        cm=0.
!        do ic=1,nch
!            call ChLims(ic,nus,nds)
!            ivar=2*nus-1  ! row location first variable (del q)
!            ncus=0
!            ncds=0
!            if(flwexst(ic)>0)then
!                call AnlzFlow(ic)
!                do n=nus,nds
!                    if(ncus==0.and.indxcrt(n)==1)ncus=n
!                    if(indxcrt(n)==1)ncds=n
!                enddo
!            endif !(flwexst(ic)>0)then
!            !write(2,*)ic,ncus,ncds
!            if(flwexst(ic)<0)then ! this channel is not active
!                iloc(ivar)=ivar  ! set first del q this chnl==0.
!                cm(ivar,1)=1.
!                rhs(ivar)=0.
!            else ! ! this channel is  active, do bc
!                select case(chbc(ic,1))
!                case(1) ! usbc is stage
!                    iloc(ivar)=ivar+1  ! column location of first non-zero coef
!                    cm(ivar,1)=1.
!                    call tsbc(ic,1,ti,zbc)
!                    rhs(ivar)=zbc-curxsc(nus,1)
!                case(2) ! usbc is discharge
!                    iloc(ivar)=ivar
!                    cm(ivar,1)=1.
!                    call tsbc(ic,1,ti,qbc)
!                    rhs(ivar)=qbc-curxsc(nus,7)
!                case(0) ! usbc is junction
!                    call ActvLocs(ic,1,indx,nfrst)  ! find sequence of active channels
!                    if(nfrst==ic)then ! this is first secn for jct, assign q continuity
!                        !  000731  this version assumes that first ch into jct always has lowest ch number of any there
!                        !  000731  this version assumes that order of chs into jct is monotonic increasing
!                        iloc(ivar)=ivar
!                        rhs(ivar)=0.
!                        isb=abs(jctsecs(ic,1,iactv(1)))
!                        do k=1,indx
!                            isec=abs(jctsecs(ic,1,iactv(k)))
!                            icol=2*(isec-isb)+1
!                            coef=sign(1.,real(jctsecs(ic,1,iactv(k))))
!                            cm(ivar,icol)=coef
!                            rhs(ivar)=rhs(ivar)-curxsc(isec,7)*coef
!                        enddo
!                    else ! this is subsequent secn for jct, assign z continuity
!                        iprv=abs(jctsecs(ic,1,iactv(1)))
!                        ishft=2*(nus-iprv)
!                        iloc(ivar)=ivar-ishft+1
!                        cm(ivar,1)=1.  ! previous z
!                        cm(ivar,ishft+1)=-1.  ! this z
!                        rhs(ivar)=curxsc(nus,1)-curxsc(iprv,1) ! -(previous z- this z)
!                    end if !(jctsecs(ic,1,2)==1)then ! this is first secn for jct, assign q continuity
!                endselect
!                ! write(2,*)ivar,iloc(ivar),rhs(ivar)
!                ! write(2,*)(cm(ivar,jj),jj=1,ntbw)
!            endif  !(flwexst(ic)<0)then ! this channel is not active
!            if(flwexst(ic)<0)then  ! this is a brute force way to handle flow in channels==0.
!                do nn=nus+1,nds
!                    ivar=2*(nn-1)
!                    iloc(ivar)=ivar
!                    cm(ivar,1)=1.
!                    rhs(ivar)=0.  ! del z=0., cause z already set==zbottom +etol
!                    ivar=ivar+1
!                    iloc(ivar)=ivar
!                    cm(ivar,1)=1.
!                    rhs(ivar)=0.  ! del q=0., cause flow already set==0.
!                enddo
!            else !if(flwexst(ic)<0)then
!                do ns=nus+1,nds
!                    ivar=2*(ns-1)
!                    nsm=ns-1
!                    if(indxcrt(ns)==1)then  !there was a tendency
!                        ! to produce critical flow at this section last iteration
!                        ncmns=ncus-1
!                        iloc(ivar)=ivar-1
!                        qq=curxsc(ncmns,7)
!                        call zCrit(ns,qq,zc)
!                        rhs(ivar)=zc-curxsc(ns,1)
!                        cm(ivar,2)=1.
!                        !				rhs(ivar)=1.-curxsc(ncmns,7)/sqrt(gg*curxsc(ns,2)**3/curxsc(ns,5))
!                        !				cm(ivar,2)=-3.*curxsc(ncmns,7)*sqrt(curxsc(ns,5)**3/curxsc(ns,2)**5/gg)/2.
!                        ivar=ivar+1
!                        iloc(ivar)=ivar-2
!                        rhs(ivar)=curxsc(ncmns,7)-curxsc(ns,7)
!                        cm(ivar,3)=1.
!                    else !if(indxcrt(ns)==1)then
!                        call dxc(ns,dx)
!                        ! continuity, even rows of coef mtrx, solve for delz
!                        bb=bbar(ns)
!                        pq=pqpx(ns)
!                        qcoef=thet/bb/dx
!                        ycoef=.5*thet*pq/bb/bb
!                        fc=pypt(ns)+pq/bb
!                        iloc(ivar)=ivar-1
!                        rhs(ivar)=-fc
!                        cm(ivar,1)=-qcoef
!                        call pbpy(nsm,py)
!                        !				print*,py
!                        cm(ivar,2)=dtinv-ycoef*py  ! diagonal term
!                        cm(ivar,3)=qcoef
!                        call pbpy(ns,py)
!                        !				print*,py
!                        cm(ivar,4)=dtinv-ycoef*py
!                        ! momentum, odd rows of coef mtrx, solve for delq
!                        ivar=ivar+1
!                        iloc(ivar)=ivar-2
!                        if(indxcrt(nsm)==1)then
!                            rhs(ivar)=curxsc(ncmns,7)-curxsc(ns,7)
!                            cm(ivar,3)=1.
!                        else !if(indxcrt(nsm)==1)then
!                            call kbar(ns,kb)
!                            ab=abar(ns)
!                            qabsb=qabsqbar(ns)
!                            xmc=pypx(ns)+qabsb/kb
!                            aj=curxsc(nsm,2)
!                            twj=curxsc(nsm,5)
!                            qj=curxsc(nsm,7)
!                            ajp=curxsc(ns,2)
!                            twjp=curxsc(ns,5)
!                            qjp=curxsc(ns,7)
!                            fm=pqpt(ns)+pqsovapx(ns)+gg*ab*xmc
!                            rhs(ivar)=-fm
!                            cm(ivar,1)=dtinv+qj*thet*(gg*ab/kb-2./dx/aj)
!                            call pkpy(nsm,pkp)
!                            cm(ivar,2)=thet*(twj*(qj/aj)**2/dx+gg*(twj*xmc/2.-ab*&
!                                (1./dx+qabsb*pkp/kb/kb/2.)))
!                            cm(ivar,3)=dtinv+qjp*thet*(gg*ab/kb+2./dx/ajp)  ! diagonal term
!                            call pkpy(ns,pkp)
!                            cm(ivar,4)=thet*(-twjp*(qjp/ajp)**2/dx+gg*(twjp*xmc/2.+ab*&
!                                (1./dx-qabsb*pkp/kb/kb/2.)))
!                        end if  !(indxcrt(nsm)==1)then
!                    endif !(indxcrt(ns)==1)then
!                enddo ! ns, section loop
!            endif !(flwexst(ic)<0)then
!            ivar=ivar+1
!            if(flwexst(ic)<0)then ! this ch is not active
!                iloc(ivar)=ivar ! set last del z==0
!                cm(ivar,1)=1.
!                rhs(ivar)=0.
!            else !if(flwexst(ic)>0)then ! this ch is  active
!                select case(chbc(ic,2))
!                case(1) ! dsbc is stage
!                    iloc(ivar)=ivar  ! column location of first non-zero coef
!                    cm(ivar,1)=1.
!                    call tsbc(ic,2,ti,zbc)
!                    rhs(ivar)=zbc-curxsc(nds,1)
!                case(3) !critical depth dsbc
!                    if(curxsc(nds,2)<=0..or.curxsc(nds,5)<=0.)then
!                        print *,nds,'area<0, crit depth ChPreis',curxsc(nds,1),curxsc(nds,2),curxsc(nds,5)
!                        stop
!                    end if  !(curxsc(nds,2)<=0..or.curxsc(nds,5<=0.)then
!                    rhs(ivar)=sqrt(gg*curxsc(nds,2)**3/curxsc(nds,5))-curxsc(nds,7)
!                    iloc(ivar)=ivar-1
!                    cm(ivar,1)= 1.
!                    cm(ivar,2)=-3.*sqrt(gg*curxsc(nds,2)*curxsc(nds,5))/2.
!                case(2) ! dsbc is discharge
!                    iloc(ivar)=ivar-1
!                    cm(ivar,1)=1.
!                    call tsbc(ic,2,ti,qbc)
!                    rhs(ivar)=qbc-curxsc(nds,7)
!                case(0) ! dsbc is junction
!                    call ActvLocs(ic,2,indx,nfrst)  ! find sequence of active channels
!                    if(nfrst==ic)then ! this is first active secn for jct, assign q continuity
!                        !  000731  this version assumes that first ch into jct always has lowest ch number of any there
!                        !  000731  this version assumes that order of chs into jct is monotonic increasing
!                        iloc(ivar)=ivar-1
!                        rhs(ivar)=0.
!                        isb=abs(jctsecs(ic,2,iactv(1)))
!                        do k=1,indx
!                            isec=abs(jctsecs(ic,2,iactv(k)))
!                            icol=2*(isec-isb)+1
!                            coef=sign(1.,real(jctsecs(ic,2,iactv(k))))
!                            cm(ivar,icol)=coef
!                            rhs(ivar)=rhs(ivar)-curxsc(isec,7)*coef
!                        enddo
!                    else ! this is subsequent secn for jct, assign z continuity
!                        iprv=abs(jctsecs(ic,2,iactv(1)))
!                        zdum=curxsc(nds,1)
!                        call zcrit(nds,curxsc(nds,7),zc)
!                        if(zc>curxsc(iprv,1))then
!                            cm(ivar,1)=1.
!                            rhs(ivar)=zc-curxsc(nds,1)
!                        else !if(zc>curxsc(iprv,1))then
!                            call Farea(nds,zdum)
!                            ishft=2*(nds-iprv)
!                            iloc(ivar)=ivar-ishft
!                            cm(ivar,1)=1.
!                            cm(ivar,ishft+1)=-1.
!                            rhs(ivar)=curxsc(nds,1)-curxsc(iprv,1)
!                        end if !(zc>curxsc(iprv,1))then
!                    end if !(jctsecs(ic,1,2)==1)then ! this is first secn for jct, assign q continuity
!                endselect
!            endif !(flwexst(ic)<0)then ! this ch is not active
!            ! write(2,*)ivar,iloc(ivar),rhs(ivar)
!            ! write(2,*)(cm(ivar,jj),jj=1,ntbw)
!        enddo ! ic=1,nch  channel loop
!
!        !	do n=1,2*nsc
!        !		write(2,*)n,iloc(n)
!        !		write(2,'(5es14.5)')(cm(n,j),j=1,4),rhs(n)
!        !	enddo
!        call matinv()
!        do n=1,nsc
!            nindx=2*n
!            curxsc(n,7)=curxsc(n,7)+rhs(nindx-1)*rlxdschg
!            curxsc(n,1)=curxsc(n,1)+rhs(nindx)*rlxstg
!            call Farea(n,curxsc(n,1))
!            ! 		write(2,*)n,curxsc(n,1),curxsc(n,7)
!        enddo
!        xmsct=0.  !total mass conservation error
!        xint=0.  !total abs mass inflow
!        do ic=1,nch
!            call ChLims(ic,nus,nds)
!            xmscon=dt*(curxsc(nus,7)-curxsc(nds,7)) ! postv if inflow exceeds outflow-volume change
!            xinflw=abs(dt*curxsc(nus,7))
!            xint=xint+xinflw
!            do n=nus+1,nsc
!                nm=n-1
!                call dxc(n,dx)
!                a1=curxsc(nm,5)*(curxsc(nm,1)-oldxsc(nm,1))
!                a2=curxsc(n,5)*(curxsc(n,1)-oldxsc(n,1))
!                xmscon=xmscon-dx*(a1+a2)/2.
!            enddo
!        end	do !ic=1,nch
!        xmsct=xmsct+xmscon
!        if(abs(xmsct)<voltol*xint.and.iter>=5)goto 1000
!        !			if(abs(xmsct)<.00001*xint.and.iter>=5)goto 1000
!    enddo ! iter loop
!1000 continue !oldxsc=curxsc
!    !			write(2,*)'iter rslts',iter,xmsct/xinflw
!    ! oplnxsc=plnxsc
!    end Subroutine CompOnly
    ! continuity subroutines
    Subroutine pbpy(j,py)
    use support
    IMPLICIT NONE                             ! SRC
    INTEGER, PARAMETER :: mp = KIND(1.0D0)    ! SRC
    REAL(KIND=mp)  DZ                         ! SRC  xFChPreis.f90(271)
    INTEGER        J                          ! SRC  xFChPreis.f90(271)
    REAL(KIND=mp)  TW                         ! SRC  xFChPreis.f90(272)
    REAL(KIND=mp)  PY                         ! SRC  xFChPreis.f90(274)
    dz=.05*curxsc(j,4)
    tw=curxsc(j,5)
    call Farea(j,curxsc(j,1)+dz)
    py=(curxsc(j,5)-tw)/dz
    call Farea(j,curxsc(j,1)-dz)
    end subroutine pbpy
    ! momentum subroutines
    Subroutine kbar(j,py)
    use support
    IMPLICIT NONE                             ! SRC
    INTEGER, PARAMETER :: mp = KIND(1.0D0)    ! SRC
    REAL(KIND=mp)  PY                         ! SRC  xFChPreis.f90(280)
    INTEGER        J                          ! SRC  xFChPreis.f90(280)
    py=.5*(thet*(curxsc(j,6)**2+curxsc(j-1,6)**2)+(1.-thet)*(oldxsc(j,6)**2+oldxsc(j-1,6)**2))
    end subroutine kbar
    Subroutine pkpy(j,py)
    use support
    IMPLICIT NONE                             ! SRC
    INTEGER, PARAMETER :: mp = KIND(1.0D0)    ! SRC
    REAL(KIND=mp)  DZ                         ! SRC  xFChPreis.f90(284)
    INTEGER        J                          ! SRC  xFChPreis.f90(284)
    REAL(KIND=mp)  XK                         ! SRC  xFChPreis.f90(285)
    REAL(KIND=mp)  PY                         ! SRC  xFChPreis.f90(287)
    dz=.05*curxsc(j,4)
    xk=curxsc(j,6)
    call Farea(j,curxsc(j,1)+dz)
    py=2.*xk*(curxsc(j,6)-xk)/dz
    call Farea(j,curxsc(j,1)-dz)
    end subroutine pkpy
    Subroutine ActvLocs(ic,iend,indx,nfrst)  ! find sequence of active channels
    use Support
    IMPLICIT NONE                             ! SRC
    INTEGER        NJ                         ! SRC  xFChPreis.f90(292)
    INTEGER        IC                         ! SRC  xFChPreis.f90(292)
    INTEGER        IEND                       ! SRC  xFChPreis.f90(292)
    INTEGER        INDX                       ! SRC  xFChPreis.f90(294)
    INTEGER        NN                         ! SRC  xFChPreis.f90(295)
    INTEGER        ICN                        ! SRC  xFChPreis.f90(296)
    INTEGER        NFRST                      ! SRC  xFChPreis.f90(300)
    nj=jctsecs(ic,iend,1)
    iactv=0
    indx=0
    do nn=3,nj+2
        icn=jctsecs(ic,iend,nn+nj)  ! entering or leaving channel nos
        if(flwexst(icn)>0)then
            indx=indx+1
            iactv(indx)=nn  ! location in junction sequence
            if(indx==1)nfrst=icn ! # of first active channel
        endif
    enddo
    !	print *,ic,iend,iactv
    end  Subroutine ActvLocs   ! find sequence of active channels
    Subroutine CheckCrit(nus,nds,icrit)
    use Support
    IMPLICIT NONE                             ! SRC
    INTEGER, PARAMETER :: mp = KIND(1.0D0)    ! SRC
    INTEGER        ICRIT                      ! SRC  xFChPreis.f90(307)
    REAL(KIND=mp)  FRDMAX                     ! SRC  xFChPreis.f90(308)
    INTEGER        IC                         ! SRC  xFChPreis.f90(309)
    INTEGER        NUS                        ! SRC  xFChPreis.f90(309)
    INTEGER        NDS                        ! SRC  xFChPreis.f90(309)
    icrit=0
    frdmax=0.
    do ic=nus,nds
        if(curxsc(ic,2)<=0.or.curxsc(ic,5)<=0.)goto 1000
        frdmax=max(abs(curxsc(ic,7))/sqrt(gg*curxsc(ic,2)**3/curxsc(ic,5)),frdmax)
    enddo !ic=nus+1,nds
    if(frdmax>.8)icrit=1
    return
1000 continue
    do ic=nus,nds
        !		write(2,'(i3,7f10.3)')ic,(curxsc(ic,j),j=1,7)
    enddo !ic=nus+1,nds
    !	write(2,*) 'stop in CheckCrit'
    stop
    end Subroutine CheckCrit
