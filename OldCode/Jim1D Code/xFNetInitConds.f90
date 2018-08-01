    Subroutine ChLims(nc,ics,ice)
    use Support
    ics=1
    if(nc>1)ics=nsce(nc-1)+1
    ice=nsce(nc)
    end Subroutine ChLims
    Subroutine NetInitCond(ti)
    use Support
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
    AvAbsElDf=1.e9
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
                    zds=-1.e9  ! Step BW will handle setting ds bc at crit depth
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

