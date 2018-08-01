    MODULE RivVarWMod
    USE RivVarMod
    USE CalcCond
    IMPLICIT NONE

    INTEGER, ALLOCATABLE, DIMENSION(:) :: isw
    REAL(kind = mp), ALLOCATABLE, DIMENSION(:) :: dea, eav, am, bm, ccm, dm, em, qp
    REAL(kind = mp), ALLOCATABLE, DIMENSION(:) :: qprms
    !	REAL(kind = mp), ALLOCATABLE, DIMENSION(:,:) :: u, v, e
    REAL(kind = mp), ALLOCATABLE, DIMENSION(:,:) :: dqsde, dqnde, qv, qu, dq, de
    REAL(kind = mp), ALLOCATABLE, DIMENSION(:,:) :: dude, dvde, dqe, vout, uout, eka
    REAL(kind = mp), ALLOCATABLE, DIMENSION(:,:) :: dsu, dne, dnq, dsv, dse, dsq
    REAL(kind = mp), ALLOCATABLE, DIMENSION(:,:) :: dnu, dnv
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: npos
    CONTAINS
    SUBROUTINE alloc_working()
    INTEGER :: status, tmpval
    IF (ns.gt.nn) THEN
        ALLOCATE(isw(ns), STAT = status)
        tmpval = ns
    ELSE
        ALLOCATE(isw(nn), STAT = status)
        tmpval = nn
    END IF
    ALLOCATE(qprms(itm*10), STAT = status)
    ALLOCATE(dea(ns), STAT = status)
    ALLOCATE(eav(ns), STAT = status)
    ALLOCATE(am(tmpval), STAT = status)
    ALLOCATE(bm(tmpval), STAT = status)
    ALLOCATE(ccm(tmpval), STAT = status)
    ALLOCATE(dm(tmpval), STAT = status)
    ALLOCATE(em(tmpval), STAT = status)
    ALLOCATE(qp(ns), STAT = status)
    !		ALLOCATE(u(ns, nn), STAT = status)
    !		ALLOCATE(v(ns, nn), STAT = status)
    !		ALLOCATE(e(ns, nn), STAT = status)
    ALLOCATE(dqsde(ns, nn), STAT = status)
    ALLOCATE(dqnde(ns, nn), STAT = status)
    ALLOCATE(qv(ns, nn), STAT = status)
    ALLOCATE(qu(ns, nn), STAT = status)
    ALLOCATE(dq(ns, nn), STAT = status)
    ALLOCATE(de(ns, nn), STAT = status)
    ALLOCATE(dude(ns, nn), STAT = status)
    ALLOCATE(dvde(ns, nn), STAT = status)
    ALLOCATE(dqe(ns, nn), STAT = status)
    ALLOCATE(vout(ns, nn), STAT = status)
    ALLOCATE(uout(ns, nn), STAT = STATUS)
    ALLOCATE(eka(ns, nn), STAT = status)
    ALLOCATE(dsu(ns, nn), STAT = status)
    ALLOCATE(dne(ns, nn), STAT = status)
    ALLOCATE(dnq(ns, nn), STAT = status)
    ALLOCATE(dsv(ns, nn), STAT = status)
    ALLOCATE(dse(ns, nn), STAT = status)
    ALLOCATE(dsq(ns, nn), STAT = status)
    ALLOCATE(dnu(ns, nn), STAT = status)
    ALLOCATE(dnv(ns, nn), STAT = status)

    ALLOCATE(npos(ns, 2), STAT = status)

    END SUBROUTINE alloc_working
    SUBROUTINE dealloc_working()
    INTEGER :: status
    DEALLOCATE(isw, STAT = status)
    DEALLOCATE(qprms, STAT = status)
    DEALLOCATE(dea, STAT = status)
    DEALLOCATE(eav, STAT = status)
    DEALLOCATE(am, STAT = status)
    DEALLOCATE(bm, STAT = status)
    DEALLOCATE(ccm, STAT = status)
    DEALLOCATE(dm, STAT = status)
    DEALLOCATE(em, STAT = status)
    DEALLOCATE(qp, STAT = status)
    !		DEALLOCATE(u, STAT = status)
    !		DEALLOCATE(v, STAT = status)
    !		DEALLOCATE(e, STAT = status)
    DEALLOCATE(dqsde, STAT = status)
    DEALLOCATE(dqnde, STAT = status)
    DEALLOCATE(qv, STAT = status)
    DEALLOCATE(qu, STAT = status)
    DEALLOCATE(dq, STAT = status)
    DEALLOCATE(de, STAT = status)
    DEALLOCATE(dude, STAT = status)
    DEALLOCATE(dvde, STAT = status)
    DEALLOCATE(dqe, STAT = status)
    DEALLOCATE(vout, STAT = status)
    DEALLOCATE(uout, STAT = status)
    DEALLOCATE(eka, STAT = status)
    DEALLOCATE(dsu, STAT = status)
    DEALLOCATE(dne, STAT = status)
    DEALLOCATE(dnq, STAT = status)
    DEALLOCATE(dsv, STAT = status)
    DEALLOCATE(dse, STAT = status)
    DEALLOCATE(dsq, STAT = status)
    DEALLOCATE(dnu, STAT = status)
    DEALLOCATE(dnv, STAT = status)

    DEALLOCATE(npos, STAT = status)

    END SUBROUTINE dealloc_working

    END MODULE RivVarWMod
