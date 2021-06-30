    !    curxsc(ns,1)=ws elevation
    !    curxsc(ns,2)=xsec area
    !    curxsc(ns,3)=hydraulic radius =area/wet perim
    !    curxsc(ns,4)=hydraulic depth =area/top width
    !    curxsc(ns,5)=top width
    !    curxsc(ns,6)=conveyance
    !	!    curxsc(ns,7)=discharge
    !  xsloc(i,5)=emin  ! xsloc( ,5) set to min ch elev
    MODULE NetPreisMain
    use Support
    use ParamSetMod
    use StpBckInputMod
    use StpBckOutputMod
    !use readNetCDFMod
    !use writeNetCDFMod
    !USE DFLOGM
    !use msflib
    !use dflib

    CONTAINS

    SUBROUTINE NETPREIS2B(tmp)
    INTEGER, INTENT(IN) :: tmp

    ENDSUBROUTINE

    SUBROUTINE NETPREIS2(ns, nn, topo, tx, ty, cd, q, stage, wse )
    IMPLICIT NONE                             ! SRC
    INTEGER, PARAMETER :: mp = KIND(1.0D0)    ! SRC
    REAL(KIND=mp)  RSTRTSET                   ! SRC  xFNetPreisMain.f90(47)
    INTEGER, INTENT(IN) :: ns, nn
    REAL(kind = mp2), DIMENSION(:), INTENT(IN)  :: topo, tx, ty
    REAL(kind = mp2), INTENT(IN) :: cd, q, stage
    REAL(kind = mp2), DIMENSION(:), INTENT(OUT) :: wse


    common/RstDbg/rstrtset
    !common/prmtrs/nprm,prmt(4,8),ndiams(4),pdiams(4,20),&
    !  gomfl(3),bcfl(3),cgomfl,cbcfl
    character *1 char
    char='X'
    call paramset()
    call setStpBwInput(ns, nn, topo, tx, ty, cd, q, stage)

    !call read_input(STR_IN)
    tstep = 0


    !Call ParamSetNetCDF(STR_IN)
    rstrtset=0.
    !	open(2,file='stuff.txt')
    !	open(8,file='wsmap.txt')
    !	call input()
    !	zus=96.388
    !	zds=95.842
    !	call StdyStt(3,2,zus,zds)
    !	print *,'restart',ic,zus,zds,curxsc(21,7)
    !	stop
    call NetInitCond(simstrt)
    call writeStpBwOutput(wse)
    !	do i = 1,nsc
    !!		write(2,'(i4,7f10.3)')i,curxsc(i,1),curxsc(i,1)+(curxsc(i,7)/curxsc(i,2))**2/2./gg &
    !!		,curxsc(i,2),curxsc(i,7),curxsc(i,5),curxsc(i,4),curxsc(i,7)/curxsc(i,2)/sqrt(gg*curxsc(i,4))
    !	enddo
    !	call ChLims(nch,ics,ice)
    !	oldxsc=curxsc
    !	oplnxsc=plnxsc
    !	hydint=1./(prmt(3,4)+1.)
    !	hydstart=	86400.*simstrt
    !	hydtime=86400.*simend
    !	dtsprint=hydint*(hydtime-hydstart)
    !	tprint=hydstart+dtsprint
    !	t=hydstart
    !!		write(2,*)t/86400.,dt
    !!		write(2,'(i4,7f10.3)')(i,curxsc(i,1),curxsc(i,1)+(curxsc(i,7)/curxsc(i,2))**2/2./gg &
    !!		,curxsc(i,2),curxsc(i,7),curxsc(i,5),curxsc(i,4),curxsc(i,7)/curxsc(i,2)/sqrt(gg*curxsc(i,4)),i=1,ice)
    !!		write(21,*)t/86400.
    !		do i=1,ice
    !!			write(21,'(i6,4f13.3,3f10.3)')i,(plnxsc(i,jj),jj=3,6),curxsc(i,1),curxsc(i,2),curxsc(i,7)
    !		enddo
    !	do while (t<hydtime)
    !		tsmin=1.e9
    !		do nc=1,nch
    !			if(flwexst(nc)>0)call MinCourTS(nc,tsm)
    !			tsmin=min(tsm,tsmin)
    !		enddo
    !		dt=tsmin
    ! 		dt=corstpmlt*tsmin
    !		t=t+dt
    !		tdys=t/86400.
    !		ttime = t/86400.
    !! 		write(2,*)t/86400.,dt
    !!		write(2,*)t/86400.,dt ,(jctelev(jj),jj=1,nnd)
    !		call SetFlowStruc(tdys)
    !		call CompOnly(tdys,dt)
    !		if(t>=tprint)then  !.or.rstrtset>0.
    !			call tsbc(1,1,tdys,elevm)
    !			nzprnt=int(25.*(tdys- simstrt)/(simend- simstrt))
    !!			write(2,*)tdys,dt
    !!			write(2,'(i4,7f10.3)')(i,curxsc(i,1),curxsc(i,1)+(curxsc(i,7)/curxsc(i,2))**2/2./gg &
    !!				,curxsc(i,2),curxsc(i,7),curxsc(i,5),curxsc(i,4),curxsc(i,7)/sqrt(gg*curxsc(i,2)**3/curxsc(i,5)),i=1,ice)
    !!			write(21,*)tdys
    !			do i=1,ice
    !!				write(21,'(i6,4f13.3,3f10.3)')i,(plnxsc(i,jj),jj=3,6),curxsc(i,1),curxsc(i,2),curxsc(i,7)
    !			enddo
    !			tprint=t+dtsprint
    !		endif
    !		if(nnd>0)call SetJctElevs()
    !! check if any new wse's have dropped below etol above bottom and set flow in
    !!   those branches to zero
    !!		call CheckFlowExist()
    !! check if any new wse's have dropped below etol above bottom and set flow in
    !!   those branches to zero
    !! advance 'old' variables to current time step
    !! advance 'old' variables to current time step
    !		oldxsc=curxsc
    !		oplnxsc=plnxsc
    !! advance 'old' variables to current time step
    !! advance 'old' variables to current time step
    !	enddo  !while (t<100000.)
    !!	call write_output(STR_IN)
    CALL FreeUp()
    END SUBROUTINE NetPreis2
    Subroutine SetJctElevs()
    use Support
    IMPLICIT NONE                             ! SRC
    INTEGER        NN                         ! SRC  xFNetPreisMain.f90(123)
    INTEGER        NJ                         ! SRC  xFNetPreisMain.f90(124)
    INTEGER        J                          ! SRC  xFNetPreisMain.f90(125)
    INTEGER        IC                         ! SRC  xFNetPreisMain.f90(126)
    INTEGER        NUS                        ! SRC  xFNetPreisMain.f90(128)
    INTEGER        NDS                        ! SRC  xFNetPreisMain.f90(128)
    do nn=1,nnd
        nj=nprjct(nn)
        do j=1,nj
            ic=jctentrs(nn,j)
            if(flwexst(ic)>0)then
                call ChLims(ic,nus,nds)
                jctelev(nn)=curxsc(nds,1)
                if(qsign(nn,j)<0.)jctelev(nn)=curxsc(nus,1)
            endif
        enddo
    enddo
    end Subroutine SetJctElevs
    Subroutine CheckFlowExist()
    use Support
    IMPLICIT NONE                             ! SRC
    INTEGER, PARAMETER :: mp = KIND(1.0D0)    ! SRC
    INTEGER        IC                         ! SRC  xFNetPreisMain.f90(137)
    INTEGER        NUS                        ! SRC  xFNetPreisMain.f90(138)
    INTEGER        NDS                        ! SRC  xFNetPreisMain.f90(138)
    INTEGER        N                          ! SRC  xFNetPreisMain.f90(140)
    INTEGER        NN                         ! SRC  xFNetPreisMain.f90(147)
    REAL(KIND=mp)  XT                         ! SRC  xFNetPreisMain.f90(157)
    INTEGER        I                          ! SRC  xFNetPreisMain.f90(158)
    REAL(KIND=mp)  DX                         ! SRC  xFNetPreisMain.f90(159)
    REAL(KIND=mp)  X                          ! SRC  xFNetPreisMain.f90(162)
    REAL(KIND=mp)  ELUS                       ! SRC  xFNetPreisMain.f90(163)
    REAL(KIND=mp)  ELDS                       ! SRC  xFNetPreisMain.f90(164)
    REAL(KIND=mp)  R                          ! SRC  xFNetPreisMain.f90(169)
    REAL(KIND=mp)  EL                         ! SRC  xFNetPreisMain.f90(170)
    INTEGER        IMIN                       ! SRC  xFNetPreisMain.f90(177)
    REAL(KIND=mp)  AMIN                       ! SRC  xFNetPreisMain.f90(178)
    REAL(KIND=mp)  SLP                        ! SRC  xFNetPreisMain.f90(200)
    REAL(KIND=mp)  QCHNL                      ! SRC  xFNetPreisMain.f90(201)
    do ic=1,nch
        call CHLims(ic,nus,nds)
        if(flwexst(ic)>0)then  ! check to see if flow should be shut off
            do n=nus,nds
                if(curxsc(n,1)<xsloc(n,5)+etol)then
                    print *,'Terminate Flow Channel',ic,n,curxsc(n,1),xsloc(n,5)
                    !				write(2,*)'Terminate Flow Channel',ic,n,curxsc(n,1),xsloc(n,5)
                    !				write(2,*)(curxsc(nn,1),nn=nus,nds)
                    !				write(2,*)(curxsc(nn,7),nn=nus,nds)
                    flwexst(ic)=-1
                    do nn=nus,nds
                        curxsc(nn,7)=0.
                        curxsc(nn,1)=xsloc(nn,5)+etol
                        call Farea(nn,curxsc(nn,1))
                    enddo
                    goto 1000
                endif  !(curxsc(n,1)<xsloc(n,5)+etol)then
            enddo
1000        continue
        else !if(flwexst(ic)<0)then  ! check to see if flow should be turned on
            xt=0.
            do i=nus+1,nds
                call dxc(i,dx)
                xt=xt+dx
            enddo
            x=0.
            elus=jctelev(nochjcts(ic,1))
            elds=jctelev(nochjcts(ic,2))
            do i=nus,nds
                dx=0.
                if(i>nus)call dxc(i,dx)
                x=x+dx
                r=x/xt
                el=elus*(1.-r)+elds*r
                if(xsloc( i,5)+etol>el)goto 1010  ! do not restart channel
            enddo
            ! restart channel
            !			write(2,*)elus,elds
            flwexst(ic)=1
            x=0.
            imin=nus
            amin=1.e9
            do i=nus,nds
                dx=0.
                if(i>nus)call dxc(i,dx)
                x=x+dx
                r=x/xt
                el=elus*(1.-r)+elds*r
                curxsc(i,1)=el
                call Farea(i,el)
                if(curxsc(i,2)<amin)then
                    amin=curxsc(i,2)
                    imin=i
                endif
            enddo
            !			call StdyStt(ic,2,elus,elds)
            print *,'Restart Flow Channel',ic,elus,elds,curxsc(nus,7)
            !				write(2,*)'Restart Flow Channel',ic,elus,elds,curxsc(nus,7)
            ! 			write(2,*)(curxsc(i,1),i=nus,nds)
            !			write(2,*)(curxsc(i,2),i=nus,nds)
            !			write(2,*)(curxsc(i,7),i=nus,nds)
            !			goto 1010
            ! for now, leave q at 0.
            slp=(elus-elds)/x
            call normq(imin,slp,qchnl)
            x=0.
            do i=nus,nds
                dx=0.
                if(i>nus)call dxc(i,dx)
                x=x+dx
                r=x/xt
                curxsc(i,1)=elus*(1.-r)+elds*r
                call Farea(i,curxsc(i,1))
                curxsc(i,7)=qchnl
            enddo
1010        continue
        endif  !(flwexst(ic)>0)thne
    enddo !ic=1,nch
    end Subroutine CheckFlowExist
    Subroutine normq(i,s,q) ! determines normal dischrge for existing wselev and given s
    use Support
    IMPLICIT NONE                             ! SRC
    INTEGER, PARAMETER :: mp = KIND(1.0D0)    ! SRC
    INTEGER        I                          ! SRC  xFNetPreisMain.f90(216)
    REAL(KIND=mp)  S                          ! SRC  xFNetPreisMain.f90(216)
    REAL(KIND=mp)  Q                          ! SRC  xFNetPreisMain.f90(216)
    REAL(KIND=mp)  QS                         ! SRC  xFNetPreisMain.f90(219)
    REAL(KIND=mp)  FO                         ! SRC  xFNetPreisMain.f90(220)
    REAL(KIND=mp)  DQ                         ! SRC  xFNetPreisMain.f90(221)
    INTEGER        ITER                       ! SRC  xFNetPreisMain.f90(222)
    REAL(KIND=mp)  F                          ! SRC  xFNetPreisMain.f90(225)
    q=sign(1.,s)*curxsc(i,2)
    qs=q
    fo=q*abs(q)/curxsc(i,6)**2-s
    dq=.1*q
    iter=0
    do while (iter<100.and.abs(dq)>.001*qs)
        q=q+dq
        f=q*abs(q)/curxsc(i,6)**2-s
        dq=-f*dq/(f-fo)
        fo=f
    enddo !while (iter<100.and.abs(dq>.001*qs)
    !	print *,i,s,q,curxsc(i,2)
    end Subroutine normq !(i,q)
    Subroutine SetFlowStruc(tdys)
    use Support
    IMPLICIT NONE                             ! SRC
    INTEGER, PARAMETER :: mp = KIND(1.0D0)    ! SRC
    REAL(KIND=mp)  TDYS                       ! SRC  xFNetPreisMain.f90(231)
    INTEGER        IC                         ! SRC  xFNetPreisMain.f90(234)
    INTEGER        NUS                        ! SRC  xFNetPreisMain.f90(236)
    INTEGER        NDS                        ! SRC  xFNetPreisMain.f90(236)
    REAL(KIND=mp)  ZUS                        ! SRC  xFNetPreisMain.f90(237)
    REAL(KIND=mp)  ZDS                        ! SRC  xFNetPreisMain.f90(237)
    REAL(KIND=mp)  ELTRY                      ! SRC  xFNetPreisMain.f90(238)
    INTEGER        N                          ! SRC  xFNetPreisMain.f90(245)
    INTEGER        J                          ! SRC  xFNetPreisMain.f90(246)
    INTEGER        NNNNN                      ! SRC  xFNetPreisMain.f90(251)
    INTEGER        NFNL                       ! SRC  xFNetPreisMain.f90(251)
    INTEGER        NS                         ! SRC  xFNetPreisMain.f90(260)
    !	write(2,*)flwexst
    do ic=1,nch
        if(chbc(ic,1)/=2.and.chbc(ic,2)/=2)then
            call ChLims(ic,nus,nds)
            call FindChElBC(ic,tdys,zus,zds)
            eltry=ellim(ic)+etol
            if(flwexst(ic)<1)then  !flow not in this channel last time step
                ! channels w q bc not involved here
                if(eltry<zus.or.eltry<zds)then
                    flwexst(ic)=1
                    !			write(2,*)ic,zus,zds
                    call StdyStt(ic,2,zus,zds)
                    do n=nus,nds
                        do j=1,7
                            oldxsc(n,j)=curxsc(n,j)
                        enddo
                    enddo
                    !			write(2,*)ic,zus,zds,curxsc(nds,7)
                    call ChLims(nch,nnnnn,nfnl)
                    !			do n=1,nfnl
                    !			write(2,*)n,curxsc(n,1),curxsc(n,2),curxsc(n,7),curxsc(n,7)/sqrt(gg*curxsc(n,2)**3/curxsc(n,5))
                    !			enddo
                endif  !(eltry<zus.or.eltry<zds)then
            else !ie flwexst(ic)>1)!flow was in this channel last time step
                !	write(2,*)eltry,zus,zds
                if(eltry>zus.and.eltry>zds)then
                    flwexst(ic)=-1
                    do ns=nus,nds
                        curxsc(ns,1)=xsloc(ns,5) +.01 !+etol
                        call Farea(ns,curxsc(ns,1))
                        curxsc(ns,7)=0.
                        do j=1,7
                            oldxsc(ns,j)=curxsc(ns,j)
                        enddo
                    enddo
                endif  !(eltry<zus.and.eltry<zds)then
            endif !(flwexst(ic)<1)then
        endif  !(chbc(ic,1)/=2.and.chbc(ic,2)/=2)then
    enddo
    !	write(2,*)flwexst
    end Subroutine SetFlowStruc
    Subroutine FindChElBC(ic,tdys,zus,zds)
    use Support
    IMPLICIT NONE                             ! SRC
    INTEGER, PARAMETER :: mp = KIND(1.0D0)    ! SRC
    INTEGER        IC                         ! SRC  xFNetPreisMain.f90(274)
    REAL(KIND=mp)  TDYS                       ! SRC  xFNetPreisMain.f90(274)
    REAL(KIND=mp)  ZUS                        ! SRC  xFNetPreisMain.f90(274)
    REAL(KIND=mp)  ZDS                        ! SRC  xFNetPreisMain.f90(274)
    if(chbc(ic,1)==0)then
        zus=jctelev(nochjcts(ic,1))
    else  !chbc==1, q bc not allowed here
        call tsbc(ic,1,tdys,zus)
    endif
    if(chbc(ic,2)==0)then
        zds=jctelev(nochjcts(ic,2))
    elseif(chbc(ic,2)==1)then  !chbc==1, q bc not allowed here
        call tsbc(ic,2,tdys,zds)
    else  ! chbc==3
        zds=-1.e9  !critical ds bc set in CompOnly
    endif
    end Subroutine FindChElBC
    Subroutine ReStart(ic,zus,zds,sq)
    use Support
    IMPLICIT NONE                             ! SRC
    INTEGER, PARAMETER :: mp = KIND(1.0D0)    ! SRC
    INTEGER        IC                         ! SRC  xFNetPreisMain.f90(289)
    REAL(KIND=mp)  ZUS                        ! SRC  xFNetPreisMain.f90(289)
    REAL(KIND=mp)  ZDS                        ! SRC  xFNetPreisMain.f90(289)
    REAL(KIND=mp)  SQ                         ! SRC  xFNetPreisMain.f90(289)
    INTEGER        NUS                        ! SRC  xFNetPreisMain.f90(293)
    INTEGER        NDS                        ! SRC  xFNetPreisMain.f90(293)
    INTEGER        I                          ! SRC  xFNetPreisMain.f90(300)
    INTEGER        J                          ! SRC  xFNetPreisMain.f90(301)
    ! 	print *,'restart',ic,zus,zds
    ! 	write(2,*)'restart',ic,zus,zds
    call ChLims(ic,nus,nds)
    call StdyStt(ic,2,zus,zds)
    !  	write(2,'(a10,2i4,5f10.3)')'restart',ic,oflwcrtcl(ic),zus,curxsc(nus,1),zds,curxsc(nds,1),curxsc(nds,7)
    if(abs(curxsc(nds,7))>.01*sq)then
        ! 	oflwexst(ic)=1
        !  	write(2,'(a10,2i4,7f10.3)')'restarthit',ic,oflwcrtcl(ic),zus,curxsc(nus,1),zds,curxsc(nds,1),curxsc(nds,7)
        !	write(2,'(7f10.3)')(curxsc(jj,7)/sqrt(gg*curxsc(jj,2)**3/curxsc(jj,5)),jj=nus,nds)
        do i=nus,nds
            do j=1,7
                oldxsc(i,j)=curxsc(i,j)
                if(j/=7)oplnxsc(i,j)=plnxsc(i,j)
            enddo
            !  		write(2,*)i,(oldxsc(i,j),j=1,7)
        enddo
    else
        flwexst(ic)=-1
        do i=nus,nds
            call Farea(i,xsloc(i,5)+etol)
            curxsc(i,7)=0.
        enddo
    endif
    end Subroutine ReStart
    end

