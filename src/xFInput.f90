    SUBROUTINE INPUT()
    use Support
    IMPLICIT NONE                             ! SRC
    INTEGER, PARAMETER :: mp = KIND(1.0D0)    ! SRC
    INTEGER        NCHX                       ! SRC  xFInput.f90(7)
    INTEGER        NSEC                       ! SRC  xFInput.f90(9)
    INTEGER        NSPT                       ! SRC  xFInput.f90(10)
    INTEGER        NSCCHMAX                   ! SRC  xFInput.f90(11)
    INTEGER        N                          ! SRC  xFInput.f90(12)
    INTEGER        NS                         ! SRC  xFInput.f90(13)
    INTEGER        I                          ! SRC  xFInput.f90(17)
    REAL(KIND=mp)  GAR                        ! SRC  xFInput.f90(18)
    REAL(KIND=mp)  GAR1                       ! SRC  xFInput.f90(18)
    REAL(KIND=mp)  GAR2                       ! SRC  xFInput.f90(18)
    REAL(KIND=mp)  GAR3                       ! SRC  xFInput.f90(18)
    INTEGER        NCD                        ! SRC  xFInput.f90(20)
    INTEGER        NNPT                       ! SRC  xFInput.f90(20)
    INTEGER        J                          ! SRC  xFInput.f90(23)
    INTEGER        NONDS                      ! SRC  xFInput.f90(29)
    INTEGER        MAXENTRS                   ! SRC  xFInput.f90(31)
    INTEGER        NOENTRS                    ! SRC  xFInput.f90(34)
    INTEGER        NGAR                       ! SRC  xFInput.f90(36)
    INTEGER        NTBC                       ! SRC  xFInput.f90(40)
    INTEGER        NBCLOCS                    ! SRC  xFInput.f90(41)
    INTEGER        NBCENT                     ! SRC  xFInput.f90(43)
    INTEGER        NSS                        ! SRC  xFInput.f90(58)
    INTEGER        NPS                        ! SRC  xFInput.f90(59)
    INTEGER        NSF                        ! SRC  xFInput.f90(70)
    INTEGER        NPF                        ! SRC  xFInput.f90(77)
    REAL(KIND=mp)  EMIN                       ! SRC  xFInput.f90(79)
    REAL(KIND=mp)  DUM                        ! SRC  xFInput.f90(81)
    INTEGER        JP                         ! SRC  xFInput.f90(84)
    INTEGER        JJ                         ! SRC  xFInput.f90(105)
    INTEGER        MAXBW                      ! SRC  xFInput.f90(124)
    INTEGER        IC                         ! SRC  xFInput.f90(142)
    INTEGER        NUS                        ! SRC  xFInput.f90(143)
    INTEGER        NDS                        ! SRC  xFInput.f90(143)
    INTEGER        IBC                        ! SRC  xFInput.f90(145)
    INTEGER        NIN                        ! SRC  xFInput.f90(149)
    INTEGER        NC                         ! SRC  xFInput.f90(150)
    INTEGER        JC                         ! SRC  xFInput.f90(151)
    INTEGER        NCC                        ! SRC  xFInput.f90(158)
    INTEGER        JCC                        ! SRC  xFInput.f90(159)
    INTEGER        NLOC                       ! SRC  xFInput.f90(160)
    INTEGER        INDX                       ! SRC  xFInput.f90(161)
    INTEGER        NU                         ! SRC  xFInput.f90(162)
    INTEGER        ND                         ! SRC  xFInput.f90(162)
    INTEGER        JMAX                       ! SRC  xFInput.f90(179)
    INTEGER        JMIN                       ! SRC  xFInput.f90(180)
    INTEGER        NOS                        ! SRC  xFInput.f90(181)
    INTEGER        NODS                       ! SRC  xFInput.f90(193)
    INTEGER        K                          ! SRC  xFInput.f90(203)
    INTEGER        NBCS                       ! SRC  xFInput.f90(213)
    dimension dum(10000)
    character(LEN=1) :: FileType
    open (1,file='..\..\BoundaryConditions\LongRch010125.txt')
    open (3,file='..\..\BoundaryConditions\BCLongRch.txt')
    read(1,*)nchx
    write(10,*)'No Channels',nchx
    nsec=0
    nspt=0
    nscchmax=0
    do n=1,nchx
        read(1,*)ns
        write(10,*)'CH ',n,' No Secns ',ns
        nsec=nsec+ns
        nscchmax=max(nscchmax,ns)
        do i=1,ns
            read(1,*) gar,gar1,gar2,gar3
            write(10,*)i,gar,gar1,gar2,gar3
            read(1,*) npt,ncd,nnpt
            nspt=nspt+npt
            !			print *,n, npt,ncd,nnpt,nspt
            if(ncd<1)read(1,*)gar,(gar1,gar2,j=1,npt)
            if(ncd>=1)read(1,*)(gar1,gar2,j=1,npt)
            read(1,*)(gar,j=1,nnpt)
        enddo
    enddo
    !	print *,nspt
    read(1,*)nonds
    write(10,*)'No nods',nonds
    maxentrs=1
    if(nonds>0)then
        do n=1,nonds
            read(1,*)noentrs
            maxentrs=max(maxentrs,noentrs)
            read(1,*)ngar
        enddo
    endif !(nonds>0)then
    read (3,*)simstrt,simend
    ntbc=0
    nbclocs=0
    do
        read(3,*,end=100)nbcent
        nbclocs=nbclocs+1
        ntbc=ntbc+nbcent
        !		print *,nbcent,ntbc
        do j=1,nbcent
            read(3,*,end=100)gar
        enddo
    enddo
100 continue
    call AllocateSecs(nchx,nscchmax,nsec,nspt,nonds,maxentrs,nbclocs)
    call AllocateTBC(ntbc)
    rewind(1)
    ellim=-1.d9
    limsec=0
    read(1,*)nchx
    nss=1
    nps=1
    do n=1,nchx
        read(1,*)ns,chbc(n,1),chbc(n,2) ! number sections, us and ds bc types
        !
        !
        ! chbc=0 node stage
        ! chbc=1 external stage
        ! chbc=2 external discharge
        ! chbc=3 critical depth
        !
        !
        nsf=nss+ns-1
        nsce(n)=nsf
        do i=nss,nsf
            read(1,*) (xsloc(i,j),j=1,5) !xr,yr,xl,yl,elev adj
            xsloc(i,6)=sqrt((xsloc(i,1)-xsloc(i,3))**2+(xsloc(i,2)-xsloc(i,4))**2) ! length
            ! of section base line
            read(1,*) npt,ncd,nnpt
            npf=nps+npt-1
            nscpte(i)=npf
            emin=1.d9
            if(ncd<1)then  !elevations are entered as offsets
                read(1,*)gar,(dum(j),j=1,2*npt)
                gar=gar+xsloc(i,5)
                do j=1,npt
                    jp=nps+j-1
                    xspts(jp,1)=dum(2*j-1)
                    xspts(jp,2)=gar-dum(2*j)
                    emin=min(emin,xspts(jp,2))
                enddo
            else  !elevs are entered directly
                read(1,*)(dum(j),j=1,2*npt)
                do j=1,npt
                    jp=nps+j-1
                    xspts(jp,1)=dum(2*j-1)
                    xspts(jp,2)=dum(2*j)+xsloc(i,5)
                    emin=min(emin,xspts(jp,2))
                enddo
            endif
            xsloc(i,5)=emin  ! xsloc( ,5) set to min ch elev
            if(emin>ellim(n))then
                ellim(n)=emin
                limsec(n)=i
            endif
            read(1,*)(dum(j),j=1,nnpt)  ! number of n's may not correspond to
            !number of pts, but there must be at least 1
            jj=0
            do j=nps,npf
                jj=jj+1
                if(jj<=nnpt)then
                    xspts(j,3)=dum(jj)
                else
                    xspts(j,3)=dum(nnpt)
                endif
            enddo
            ! 	do j=1,3
            !	print *,(xspts(jj,j),jj=nps,npf)
            !	enddo
            nps=npf+1
        enddo !end section loop
        nss=nsf+1
    enddo !end channel loop
    read(1,*)nonds
    !	print *,nonds
    jctsecs=0
    maxbw=1
    if(nonds>0)then
        nochjcts=0
        qsign=1.
        do n=1,nonds
            read(1,*)nprjct(n)
            read(1,*)(jctentrs(n,j),j=1,nprjct(n))
            do j=1,nprjct(n)
                if(jctentrs(n,j)<0)then
                    qsign(n,j)=-1.
                    jctentrs(n,j)=-jctentrs(n,j)
                    nochjcts(jctentrs(n,j),1)=n
                else
                    nochjcts(jctentrs(n,j),2)=n
                endif
            enddo
            !	print *,nprjct(n),(jctentrs(n,j),qsign(n,j),j=1,nprjct(n))
        enddo
        do ic=1,nch  !find cross-sections that connect at nodes and place in jctsecs(ic,2,ndsmx)
            call ChLims(ic,nus,nds)
            nsec=nus
            do ibc=1,2  ! do for both ends of channel
                if(ibc==2)nsec=nds
                if(chbc(ic,ibc)==0)then ! do only for channel ends at junctions
                    do n=1,nonds  ! check all nodes for this channel (end)
                        nin=nprjct(n)
                        do nc=1,nin  ! look at each channel entering node n
                            jc=jctentrs(n,nc)
                            !							print *, n,nc,jc,ic
                            if(jc==ic.and.((ibc==1.and.qsign(n,nc)<0.).or.(ibc==2.and.qsign(n,nc)>0.)))then 
                                !the channel at node n matches ic
                                !							print *, n,nc,jc
                                jctsecs(ic,ibc,1)=nin  ! no at junction
                                jctsecs(ic,ibc,2)=nc  ! location ic channel
                                do ncc=1,nin  ! find section numbers for other channels at jct
                                    jcc=jctentrs(n,ncc)
                                    nloc=ncc+2
                                    indx=int(qsign(n,ncc))
                                    call ChLims(jcc,nu,nd)
                                    nss=nu
                                    if(indx>0)nss=nd
                                    jctsecs(ic,ibc,nloc)=indx*nss
                                    jctsecs(ic,ibc,nloc+nin)=jcc
                                enddo !ncc=1,nin  ! find section numbers for other channels at jct
                            endif  !(jc==ic)then  ! the channel at node n matches ic
                        enddo !nc=1,nin  ! look at each channel entering node n
                    enddo !n=1,nonds  ! check all nodes for this channel (end)
                endif  !(chbc(ic,ibc)==0)then ! do only for channel ends at junctions
                !		print *,ic,nsec
                ! 		write(2,*)ic,(jctsecs(ic,ibc,jj),jj=1,maxentrs+2)
            enddo !ibc=1,2  ! do for both ends of channel
        enddo !ic=1,nch
        do ic=1,nch  !find maxbandwidth
            do ibc=1,2
                !		write(2,*)ic,(jctsecs(ic,ibc,jj),jj=1,2*maxentrs+2)
                jmax=0
                jmin=100000
                do nos=1,jctsecs(ic,ibc,1)
                    jmax=max(jmax,abs(jctsecs(ic,ibc,nos+2)))
                    jmin=min(jmin,abs(jctsecs(ic,ibc,nos+2)))
                enddo
                maxbw=max(maxbw,jmax-jmin)
                !			write(2,*)jmax,jmin,maxbw
            enddo
        enddo !ic=1,nch  !find maxbandwidth
    endif ! (nonds>0)then
    !		call SetBW(2*(maxbw+1))
    !	print *,(nm,(nochjcts(nm,j),j=1,2),nm=1,nch)
    ! set initial channel locations and areas
    nods=nsce(nch)
    do i=1,nods
        call Farea(i,xsloc(i,5)+etol)
    enddo
    rewind(3)
    read (3,*)simstrt,simend
    ntbc=0
    tsbcloc=0
    do n=1,nch
        !	print*,chbc(n,1),chbc(n,2)
        do k=1,2
            !
            !
            ! chbc=0 node stage
            ! chbc=1 external stage
            ! chbc=2 external discharge
            !
            !
            if(chbc(n,k)==1.or.chbc(n,k)==2)then
                tsbcloc(n,2*k-1)=ntbc+1
                read(3,*)nbcs
                do j=1,nbcs
                    jp=ntbc+j
                    read(3,*)tmsrbc(jp,1),tmsrbc(jp,2)
                enddo
                ntbc=ntbc+nbcs
                tsbcloc(n,2*k)=ntbc
            endif
        enddo
    enddo
    read(3,*)(jctelev(k),k=1,nnd)
    END SUBROUTINE INPUT
