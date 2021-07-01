    Subroutine StdyStt(ich,itp,qus,zds)
    !    curxsc(ns,1)=ws elevation
    !    curxsc(ns,2)=xsec area
    !    curxsc(ns,3)=hydraulic radius =area/wet perim
    !    curxsc(ns,4)=hydraulic depth =area/top width
    !    curxsc(ns,5)=top width
    !    curxsc(ns,6)=conveyance
    !	!    curxsc(ns,7)=discharge
    !  xsloc(i,5)=emin  ! xsloc( ,5) set to min ch elev
    use support
    IMPLICIT NONE                             ! SRC
    INTEGER, PARAMETER :: mp = KIND(1.0D0)    ! SRC
    INTEGER        NUS                        ! SRC  xFlowSupport.F90(11)
    INTEGER        ICH                        ! SRC  xFlowSupport.F90(12)
    INTEGER        NDS                        ! SRC  xFlowSupport.F90(13)
    INTEGER        ITP                        ! SRC  xFlowSupport.F90(14)
    REAL(KIND=mp)  QUS                        ! SRC  xFlowSupport.F90(17)
    REAL(KIND=mp)  ZDS                        ! SRC  xFlowSupport.F90(17)
    REAL(KIND=mp)  USSTG                      ! SRC  xFlowSupport.F90(20)
    INTEGER        I                          ! SRC  xFlowSupport.F90(22)
    REAL(KIND=mp)  QCRIT                      ! SRC  xFlowSupport.F90(29)
    INTEGER        INDRCN                     ! SRC  xFlowSupport.F90(30)
    REAL(KIND=mp)  BCSTG                      ! SRC  xFlowSupport.F90(31)
    REAL(KIND=mp)  STGFCT                     ! SRC  xFlowSupport.F90(32)
    INTEGER        ISSTG                      ! SRC  xFlowSupport.F90(33)
    REAL(KIND=mp)  QFCT                       ! SRC  xFlowSupport.F90(42)
    REAL(KIND=mp)  TOL                        ! SRC  xFlowSupport.F90(43)
    INTEGER        ITER                       ! SRC  xFlowSupport.F90(44)
    REAL(KIND=mp)  FO                         ! SRC  xFlowSupport.F90(46)
    REAL(KIND=mp)  DQ                         ! SRC  xFlowSupport.F90(47)
    REAL(KIND=mp)  F                          ! SRC  xFlowSupport.F90(48)
    nus=1
    if(ich>1)nus=nsce(ich-1)+1
    nds=nsce(ich)
    select case(itp)
    case(1)
        ! itp==1, discharge-stage stepbw
        call StepBW(1,ich,qus,zds)
    case(2)
        ! itp==2, stage-stage stepbw, solve for flow in channel, qus==us stage
        usstg=qus
        if(usstg==zds)then ! for usstg>zds
            do i=nus,nds
                call Farea(i,usstg)
                curxsc(i,7)=0.
            enddo
        else ! ! for usstg/=zds
            if(usstg>zds)then
                call Farea(nus,usstg)
                qcrit=sqrt(gg*curxsc(nus,2)**3/curxsc(nus,5))
                indrcn=1
                bcstg=zds
                stgfct=usstg
                isstg=nus
            else !if(usstg>zds)then
                call Farea(nds,zds)
                qcrit=sqrt(gg*curxsc(nds,2)**3/curxsc(nds,5))
                indrcn=-1
                bcstg=usstg
                stgfct=zds
                isstg=nds
            endif !(usstg>zds)then
            qfct=qcrit*.95
            tol=.001*curxsc(nus,4)
            iter=0
            call StepBw(indrcn,ich,qfct,bcstg)
            fo= curxsc(isstg,1)-stgfct
            dq=.01*qfct
            f=1.d9
            ! 	  write(2,*)indrcn,isstg,usstg,zds,stgfct,fo,qfct,dq
            do while(abs(f)>tol.and.iter<100)
                iter=iter+1
                if(qfct+dq>qcrit)then
                    dq=qcrit-qfct
                elseif(qfct+dq<0.)then
                    dq=-qfct+tol
                endif
                qfct=qfct+dq
                call StepBw(indrcn,ich,qfct,bcstg)
                !      		write(2,*)iter,qfct,f,curxsc(nus,1),dq
                f=curxsc(isstg,1)-stgfct
                if(f==fo)return
                dq=-.1*f*dq/(f-fo) !
                fo=f
            enddo
        endif !(usstg==zds)then ! for usstg>zds
    end select
    end Subroutine StdyStt
    subroutine StepBW(inddrcn,nc,q,zds)
    !  if inddrcn==1, usual step bw, compute ds>>us
    !  if inddrcn==-1, ds wselev>us, compute us>>ds
    use support
    IMPLICIT NONE                             ! SRC
    INTEGER, PARAMETER :: mp = KIND(1.0D0)    ! SRC
    INTEGER        NC                         ! SRC  xFlowSupport.F90(74)
    INTEGER        NFRST                      ! SRC  xFlowSupport.F90(75)
    INTEGER        NLST                       ! SRC  xFlowSupport.F90(77)
    INTEGER        NUS                        ! SRC  xFlowSupport.F90(78)
    INTEGER        NDS                        ! SRC  xFlowSupport.F90(79)
    INTEGER        INSTP                      ! SRC  xFlowSupport.F90(80)
    INTEGER        INDDRCN                    ! SRC  xFlowSupport.F90(81)
    REAL(KIND=mp)  SFDNOM                     ! SRC  xFlowSupport.F90(86)
    REAL(KIND=mp)  Q                          ! SRC  xFlowSupport.F90(86)
    REAL(KIND=mp)  ATRM                       ! SRC  xFlowSupport.F90(87)
    REAL(KIND=mp)  ZC                         ! SRC  xFlowSupport.F90(88)
    REAL(KIND=mp)  ZDS                        ! SRC  xFlowSupport.F90(89)
    REAL(KIND=mp)  ZS                         ! SRC  xFlowSupport.F90(96)
    REAL(KIND=mp)  ADS                        ! SRC  xFlowSupport.F90(98)
    INTEGER        NS                         ! SRC  xFlowSupport.F90(101)
    INTEGER        NSP                        ! SRC  xFlowSupport.F90(102)
    REAL(KIND=mp)  DZ                         ! SRC  xFlowSupport.F90(104)
    REAL(KIND=mp)  DX                         ! SRC  xFlowSupport.F90(105)
    REAL(KIND=mp)  FCRIT                      ! SRC  xFlowSupport.F90(108)
    REAL(KIND=mp)  Z                          ! SRC  xFlowSupport.F90(111)
    REAL(KIND=mp)  ZPREC                      ! SRC  xFlowSupport.F90(123)
    REAL(KIND=mp)  FO                         ! SRC  xFlowSupport.F90(124)
    INTEGER        ITER                       ! SRC  xFlowSupport.F90(125)
    REAL(KIND=mp)  F                          ! SRC  xFlowSupport.F90(126)
    real(kind = mp2) kds,kus

    oflwcrtcl(nc)=1
    nfrst=1
    if(nc>1) nfrst=nsce(nc-1)+1
    nlst=nsce(nc)
    nus=nfrst
    nds=nlst
    instp=-1
    if(inddrcn<0)then
        nus=nlst
        nds=nfrst
        instp=1
    endif !(inddrcn<0)then
    sfdnom=q*abs(q)*2.
    atrm=2.*q*q/gg
    call zCrit(nds,abs(q),zc)
    if(zc>zds)then
        ! 		print *,'using critical downstream elev',nc,nds,zc,zds,q
        ! 		write(2,1001)nc,nds,zc,zds,q
        oflwcrtcl(nc)=-1
1001    format('Using critical downstream elev ',2i6,3f10.3)
    endif
    curxsc(nds,1)=max(zds,zc+.001)
    zs=curxsc(nds,1)
    call Farea(nds,zs)
    ads=curxsc(nds,2)
    kds=curxsc(nds,6)**2
    curxsc(nds,7)=q*real(inddrcn)
    do ns=nds+instp,nus,instp
        nsp=ns+1
        if(inddrcn<0)nsp=ns
        dz=.1
        call dxc(nsp,dx)
        call zCrit(ns,abs(q),zc)
        call Farea(ns,zc)
        fcrit=zs-zc+dx*sfdnom/(kds+curxsc(ns,6)**2)+atrm*(1./ads-1./curxsc(ns,2))/(ads+curxsc(ns,2))
        ! 	print *,ns,fcrit,zs,dx*q*q/kds,q/sqrt(GG*curxsc(ns,2)**3/curxsc(ns,5))
        if(fcrit<=0)then  ! there is no subcritical upstream solution
            z=zc+.001
            call Farea(ns,z)
            !		print *,'Supercritical flow solution'
            ! 		write(2,1002)nc,ns,zc,zs,q
            oflwcrtcl(nc)=-1
1002        format('Supercritical flow solution',2i6,3f10.3)
            goto 1000
        endif !(fcrit<0)then  ! there is no subcritical upstream solution
        z=max(zc+dz,curxsc(ns-instp,1))
        call Farea(ns,z)
        call dxc(nsp,dx)
        !  	print *,ns,q,ffroud,z
        zprec=curxsc(ns,4)  !hd
        fo=zs-z+dx*sfdnom/(kds+curxsc(ns,6)**2)+atrm*(1./ads-1./curxsc(ns,2))/(ads+curxsc(ns,2))
        iter=0
        f=1.d9
        do while(abs(f)>.001*zprec.and.iter<100)
            iter=iter+1
            if(z+dz<zc)then
                dz=zc+.001-z
            endif
            z=z+dz
            call Farea(ns,z)
            call dxc(nsp,dx)
            f=zs-z+dx*sfdnom/(kds+curxsc(ns,6)**2)+atrm*(1./ads-1./curxsc(ns,2))/(ads+curxsc(ns,2))
            if(f==fo)goto 1000
            dz=-f*dz/(f-fo)
            fo=f
        enddo
1000    continue
        !	print *,ns,q/sqrt(GG*curxsc(ns,2)**3/curxsc(ns,5))
        zs=z
        curxsc(ns,7)=q*real(inddrcn)
        ads=curxsc(ns,2)
        kds=curxsc(ns,6)**2
    enddo  !section loop, ns
    end subroutine StepBW
    subroutine zCrit(is,q,zc)
    use Support
    IMPLICIT NONE                             ! SRC
    INTEGER, PARAMETER :: mp = KIND(1.0D0)    ! SRC
    REAL(KIND=mp)  ZSV                        ! SRC  xFlowSupport.F90(150)
    INTEGER        IS                         ! SRC  xFlowSupport.F90(150)
    REAL(KIND=mp)  QTOL                       ! SRC  xFlowSupport.F90(151)
    REAL(KIND=mp)  Q                          ! SRC  xFlowSupport.F90(151)
    REAL(KIND=mp)  ZC                         ! SRC  xFlowSupport.F90(153)
    REAL(KIND=mp)  F                          ! SRC  xFlowSupport.F90(154)
    REAL(KIND=mp)  DZ                         ! SRC  xFlowSupport.F90(155)
    REAL(KIND=mp)  FO                         ! SRC  xFlowSupport.F90(166)
    INTEGER        ITER                       ! SRC  xFlowSupport.F90(168)
    zsv=curxsc(is,1)
    qtol=.001*q
    qtol=min(qtol,.001)
    zc=xsloc(is,5)+.01 !etol
    f=-1.d9
    dz=0.1
    do while (f<0.)
        zc=zc+dz
        call Farea(is,zc)
        if(curxsc(is,2)<=0..or.curxsc(is,5)<=0.)then
            print *,' 1 area negative in zCrit',is,curxsc(is,1),curxsc(is,2),curxsc(is,5)
            stop
        endif
        f=rtgg*sqrt(curxsc(is,2)**3/curxsc(is,5))-q
    enddo
    call Farea(is,zc)
    fo=rtgg*sqrt(curxsc(is,2)**3/curxsc(is,5))-q
    f=1.d9
    iter=0
    do while (abs(f)>qtol.and.iter<100)
        iter=iter+1
        if(zc+dz<xsloc(is,5)+.01)then  !etol
            dz=xsloc(is,5)-zc+.01 !etol
        endif
        zc=zc+dz
        call Farea(is,zc)
        if(curxsc(is,2)<=0..or.curxsc(is,5)<=0.)then
            print *,'f area negative in zCrit',is,curxsc(is,1),curxsc(is,2),curxsc(is,5)
            stop
        endif
        ! 		print *,is,iter,zc,f,xsloc(is,5),curxsc(is,2)
        f=rtgg*sqrt(curxsc(is,2)**3/curxsc(is,5))-q
        if(f==fo)return
        dz=-f*dz/(f-fo)
        fo=f
    enddo
    end subroutine zCrit
    subroutine MinCourTS(ich,tsmin)
    use support
    IMPLICIT NONE                             ! SRC
    INTEGER, PARAMETER :: mp = KIND(1.0D0)    ! SRC
    INTEGER        ICH                        ! SRC  xFlowSupport.F90(189)
    INTEGER        NUS                        ! SRC  xFlowSupport.F90(189)
    INTEGER        NDS                        ! SRC  xFlowSupport.F90(189)
    INTEGER        IFUCK                      ! SRC  xFlowSupport.F90(190)
    INTEGER        I                          ! SRC  xFlowSupport.F90(191)
    REAL(KIND=mp)  TSMIN                      ! SRC  xFlowSupport.F90(202)
    REAL(KIND=mp)  DX                         ! SRC  xFlowSupport.F90(204)
    REAL(KIND=mp)  VEL                        ! SRC  xFlowSupport.F90(205)
    REAL(KIND=mp)  C                          ! SRC  xFlowSupport.F90(206)
    Call ChLims(ich,nus,nds)
    ifuck=0
    do i=nus,nds
        if(curxsc(i,2)<=0..or.curxsc(i,5)<=0.)then
            ifuck=1
        endif
    enddo
    if(ifuck==1)then
        do i=nus,nds
            !	write(2,'(i4,7f10.3)')i,(curxsc(i,j),j=1,7)
        enddo
        stop
    endif
    tsmin=1.d9
    do i=nus+1,nds
        call dxc(i,dx)
        vel=abs(curxsc(i,7)/curxsc(i,2))
        c=sqrt(gg*curxsc(i,2)/curxsc(i,5))
        tsmin=min(tsmin,dx/(vel+c))
        !		print *,i,dx,vel,c
    enddo
    end subroutine MinCourTS
    subroutine CheckFlow(ich,elus,elds)
    use support
    IMPLICIT NONE                             ! SRC
    INTEGER, PARAMETER :: mp = KIND(1.0D0)    ! SRC
    REAL(KIND=mp)  ELTRY                      ! SRC  xFlowSupport.F90(213)
    INTEGER        ICH                        ! SRC  xFlowSupport.F90(213)
    REAL(KIND=mp)  ELUS                       ! SRC  xFlowSupport.F90(214)
    REAL(KIND=mp)  ELDS                       ! SRC  xFlowSupport.F90(214)
    INTEGER        NUS                        ! SRC  xFlowSupport.F90(215)
    INTEGER        NDS                        ! SRC  xFlowSupport.F90(215)
    INTEGER        I                          ! SRC  xFlowSupport.F90(217)
    eltry=ellim(ich)+etol
    if( eltry<elus.or.eltry<elds)return
    call ChLims(ich,nus,nds)
    flwexst(ich)=-1
    do i=nus,nds
        call Farea(i,xsloc(i,5)+.01) !+etol)
        curxsc(i,7)=0.
    enddo
    end subroutine CheckFlow
