    Module Support
    INTEGER, PARAMETER :: mp2 = KIND(1.0D0)
    real(kind = mp2),allocatable,dimension(:)::jctelev,ojctelev,ellim
    real(kind = mp2),allocatable,dimension(:)::rhs
    integer,allocatable,dimension(:)::indxcrt,nsce,nscpte,nprjct,flwexst,oflwexst,oflwcrtcl,iloc &
        ,iactv,limsec
    integer,allocatable,dimension(:,:)::chbc,nochjcts,jctentrs,tsbcloc
    real(kind = mp2),allocatable,dimension(:,:)::xsloc,xspts,xsecd,curxsc,oldxsc,plnxsc,oplnxsc,qsign,&
        tmsrbc,wtbc
    real(kind = mp2),allocatable,dimension (:,:):: cm !coefmatrix
    real(kind = mp2),allocatable,dimension (:,:):: xend !locations of x-sec endpoints
    real(kind = mp2) :: ttime
    integer,allocatable,dimension(:,:,:)::jctsecs
    integer(2) vwx,vwy
    integer nch,nsc,npt,nnd,nbcloc,ntsbc,ntbw,itersmx,dragtype
    real(kind = mp2) gg,thet,etol,rtgg,simstrt,simend,corstpmlt,voltol,rlxdschg,rlxstg
    REAL(kind = mp2)::prmt(4,8)


    INTEGER:: tstep

    contains
    subroutine AllocateTBC(ntbc)
    ntsbc = ntbc
    allocate(tmsrbc(ntsbc,2))
    allocate(wtbc(ntsbc,nnd))
    end subroutine AllocateTBC

    subroutine AllocateSecs(nchs,nsmax,nsecs,npts,nonds,ndsmx,nbclocs)
    gg=9.8
    rtgg=sqrt(gg)
    etol=.2
    !	thet=0.7	! time step weighting factor applied to future step in iterative solution
    !	corstpmlt=5.  ! min courant time step multiplier for setting dt used in time stepping
    nch=nchs
    nsc=nsecs
    npt=npts
    nnd=nonds
    nbcloc=nbclocs
    ! section variables
    allocate(xsloc(nsecs,6))
    allocate(xend(nsecs,4))
    allocate(xspts(npts,3))
    allocate(nscpte(nsc))
    ! channel variables
    allocate(nsce(nch))
    allocate(ellim(nch))
    allocate(limsec(nch))
    allocate(chbc(nch,2))
    allocate(flwexst(nch))
    allocate(oflwexst(nch))
    allocate(oflwcrtcl(nch))
    allocate(nochjcts(nch,2))
    allocate(tsbcloc(nch,4))
    allocate(indxcrt(nsecs))
    allocate(curxsc(nsecs,7))
    allocate(oldxsc(nsecs,7))
    allocate(plnxsc(nsecs,6))
    allocate(oplnxsc(nsecs,6))
    allocate(jctsecs(nch,2,2*ndsmx+2))
    ! junction variables
    allocate(nprjct(nnd))
    allocate(jctentrs(nnd,ndsmx))
    allocate(qsign(nnd,ndsmx))
    allocate(iactv(ndsmx))
    allocate(jctelev(nnd))
    allocate(ojctelev(nnd))
    return
    entry SetBW(mbw)
    ntbw=mbw
     allocate(rhs(2*nsc))
     allocate(cm(2*nsc,ntbw))
   
    allocate(iloc(2*nsc))
    end subroutine AllocateSecs
    end Module Support
    Subroutine FreeUp()
    use Support
    !	deallocate(jctelev,ojctelev,rhs)
    deallocate(jctelev,ojctelev)
    !	deallocate(nsce,nscpte,nprjct,flwexst,oflwexst,oflwcrtcl,iloc ,iactv)
    deallocate(nsce,nscpte,nprjct,flwexst,oflwexst,oflwcrtcl,iactv)
    deallocate(chbc,nochjcts,jctentrs,tsbcloc,indxcrt)
    deallocate(xsloc,xspts,curxsc,oldxsc,plnxsc,oplnxsc,qsign,tmsrbc,wtbc)
    !	deallocate(cm,xend, jctsecs,limsec,ellim)
    deallocate(xend, jctsecs,limsec,ellim)
    endSubroutine FreeUp