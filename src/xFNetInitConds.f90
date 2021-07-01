    Subroutine ChLims(nc,ics,ice)
    use Support
    IMPLICIT NONE                             ! SRC
    INTEGER, PARAMETER :: mp = KIND(1.0D0)    ! SRC
    INTEGER        ICS                        ! SRC  xFNetInitConds.f90(3)
    INTEGER        NC                         ! SRC  xFNetInitConds.f90(4)
    INTEGER        ICE                        ! SRC  xFNetInitConds.f90(5)
    ics=1
    if(nc>1)ics=nsce(nc-1)+1
    ice=nsce(nc)
    end Subroutine ChLims
    Subroutine NetInitCond(ti)
    use Support
    IMPLICIT NONE                             ! SRC
    INTEGER, PARAMETER :: mp = KIND(1.0D0)    ! SRC
    INTEGER        NC                         ! SRC  xFNetInitConds.f90(14)
    REAL(KIND=mp)  TI                         ! SRC  xFNetInitConds.f90(16)
    REAL(KIND=mp)  QQ                         ! SRC  xFNetInitConds.f90(16)
    REAL(KIND=mp)  ZZ                         ! SRC  xFNetInitConds.f90(17)
    INTEGER        ICS                        ! SRC  xFNetInitConds.f90(18)
    INTEGER        ICE                        ! SRC  xFNetInitConds.f90(18)
    REAL(KIND=mp)  ZUS                        ! SRC  xFNetInitConds.f90(27)
    REAL(KIND=mp)  ZDS                        ! SRC  xFNetInitConds.f90(28)
    REAL(KIND=mp)  AVABSELDF                  ! SRC  xFNetInitConds.f90(37)
    INTEGER        ITEROUT                    ! SRC  xFNetInitConds.f90(38)
    INTEGER        NN                         ! SRC  xFNetInitConds.f90(42)
    REAL(KIND=mp)  SQ                         ! SRC  xFNetInitConds.f90(45)
    REAL(KIND=mp)  FO                         ! SRC  xFNetInitConds.f90(48)
    REAL(KIND=mp)  DZ                         ! SRC  xFNetInitConds.f90(49)
    INTEGER        ITER                       ! SRC  xFNetInitConds.f90(50)
    REAL(KIND=mp)  F                          ! SRC  xFNetInitConds.f90(54)
    !	print *,nch,nnd,nsc
    flwexst=1  ! master index for flow in channels, set to -1 if depth at any section
    ! drops below etol from bottom
    oflwcrtcl=1
    if(nch==1)then
        nc=1
        if(chbc(nc,1)==2.and.chbc(nc,2)==1)then  !discharge-stage
            call tsbc(nc,1,ti,qq)
            call tsbc(nc,2,ti,zz)
            call ChLims(nc,ics,ice)
            call StdyStt(nc,1,qq,zz)
        elseif(chbc(nc,1)==2.and.chbc(nc,2)==3)then  !discharge-stage
            call tsbc(nc,1,ti,qq)
            call tsbc(nc,2,ti,zz)
            call ChLims(nc,ics,ice)
            call zCrit(ice,qq,zz)
            call StdyStt(nc,1,qq,zz)
        elseif(chbc(nc,1)==1.and.chbc(nc,2)==1)then   !stage-stage
            call tsbc(nc,1,ti,zus)
            call tsbc(nc,2,ti,zds)
            call StdyStt(nc,2,zus,zds)
        else
            print *,'An inappropriate boundary condition combination has been chosen'
            print *,'for a single-channel simulation.  Stop'
            stop
        endif
        return
    endif
    AvAbsElDf=1.d9
    iterout=0
    do while (AvAbsElDf>.001.and.iterout<100)
        iterout=iterout+1
        ojctelev=jctelev
        do nn=nnd,1,-1
            zz=jctelev(nn)
            !     print *,'Node',nn,zz
            call StgFct(nn,ti,zz,sq)  ! node balance q based on any jct elev and bc or
            ! jct elevs at other end of ch
            !	print *,sq
            fo=-sq
            dz=.1
            iter=0
            do while (iter<100.and.abs(dz)>.001)
                zz=zz+dz
                call StgFct(nn,ti,zz,sq)
                f=-sq
                dz=-f*dz/(f-fo)
                iter=iter+1
                fo=f
                !		  print *,iter,zz,dz,f
            enddo
            jctelev(nn)=zz
        enddo !nn=nnd,1,-1, node sums
        AvAbsElDf=0.
        do nn=1,nnd
            AvAbsElDf=AvAbsElDf+abs(jctelev(nn)-ojctelev(nn))
        enddo
        AvAbsElDf=AvAbsElDf/real(nnd)
        !   	   print *,' iterout',iterout,AvAbsElDf
    enddo !while (AvAbsElDf>.001.and.iter<100)
    !   				print *,jctelev
    !  				print *,(curxsc(i,7),i=1,nsce(nch))
    end Subroutine NetInitCond
    Subroutine StgFct(nn,ti,e,sq)
    ! calculates node balance q based on any jct elev and bc or
    ! jct elevs at other end of ch
    use Support
    IMPLICIT NONE                             ! SRC
    INTEGER, PARAMETER :: mp = KIND(1.0D0)    ! SRC
    INTEGER        NJ                         ! SRC  xFNetInitConds.f90(76)
    INTEGER        NN                         ! SRC  xFNetInitConds.f90(76)
    REAL(KIND=mp)  SQ                         ! SRC  xFNetInitConds.f90(77)
    INTEGER        I                          ! SRC  xFNetInitConds.f90(78)
    INTEGER        IC                         ! SRC  xFNetInitConds.f90(79)
    INTEGER        ICS                        ! SRC  xFNetInitConds.f90(81)
    INTEGER        ICE                        ! SRC  xFNetInitConds.f90(81)
    REAL(KIND=mp)  TI                         ! SRC  xFNetInitConds.f90(84)
    REAL(KIND=mp)  QQ                         ! SRC  xFNetInitConds.f90(84)
    REAL(KIND=mp)  E                          ! SRC  xFNetInitConds.f90(85)
    INTEGER        ICRIT                      ! SRC  xFNetInitConds.f90(86)
    REAL(KIND=mp)  ZUS                        ! SRC  xFNetInitConds.f90(92)
    REAL(KIND=mp)  ZDS                        ! SRC  xFNetInitConds.f90(108)
    nj=nprjct(nn)
    sq=0.
    do i=1,nj  !do for all channels entering jct
        ic=jctentrs(nn,i)
        if(flwexst(ic)>0)then
            call ChLims(ic,ics,ice)
            if(qsign(nn,i)>0.)then !ds end of ch is at jct
                if(chbc(ic,1)==2)then
                    call tsbc(ic,1,ti,qq)
                    call StdyStt(ic,1,qq,e)
                    call CheckCrit(ics,ice,icrit)
                    if(icrit>0)oflwcrtcl(ic)=-1
                    sq=sq+qq
                    !				print *,'1 ',ics,ice,sq
                else
                    if(chbc(ic,1)==1)then
                        call tsbc(ic,1,ti,zus)
                    else
                        zus=jctelev(nochjcts(ic,1))
                    endif
                    call CheckFlow(ic,zus,e)
                    if(flwexst(ic)>0)then
                        call StdyStt(ic,2,zus,e)
                        call CheckCrit(ics,ice,icrit)
                        if(icrit>0)oflwcrtcl(ic)=-1
                        sq=sq+curxsc(ice,7)
                    endif
                    !				print *,'2 ',ics,ice,sq,zus,e
                endif
            else !(qsign(nn,i)<0.)us end of ch is at jct
                ! 000616 current programming here does not allow for downstream q specified
                if(chbc(ic,2)==0)then
                    zds=jctelev(nochjcts(ic,2))
                elseif(chbc(ic,1)==1)then
                    call tsbc(ic,2,ti,zds)
                else !chbc(ic,1)==3
                    zds=-1.d9  ! Step BW will handle setting ds bc at crit depth
                endif
                !				print *,'3a ',ics,ice,sq,e,zds
                call CheckFlow(ic,e,zds)
                if(flwexst(ic)>0)then
                    call StdyStt(ic,2,e,zds)
                    call CheckCrit(ics,ice,icrit)
                    if(icrit>0)oflwcrtcl(ic)=-1
                    sq=sq-curxsc(ics,7)
                endif
                !				print *,'3 ',ics,ice,sq,e,zds
            endif !(qsign(nn,i)>0.)which? end of ch is at jct
        endif !(flwexst(ic)>0)then
    enddo !i=1,nj  !do for all channels entering jct
    end Subroutine StgFct

