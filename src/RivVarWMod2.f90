    MODULE RivVarWMod2
    !USE RivVarMod
    !USE CalcCond
    use fm_global
    IMPLICIT NONE
    type riv_w_var
        INTEGER, ALLOCATABLE, DIMENSION(:) :: isw
        REAL(kind = mp), ALLOCATABLE, DIMENSION(:) :: dea, eav, am, bm, ccm, dm, em, qp
        REAL(kind = mp), ALLOCATABLE, DIMENSION(:) :: qprms
        !	REAL(kind = mp), ALLOCATABLE, DIMENSION(:,:) :: u, v, e
        REAL(kind = mp), ALLOCATABLE, DIMENSION(:,:) :: dqsde, dqnde, qv, qu, dq, de
        REAL(kind = mp), ALLOCATABLE, DIMENSION(:,:) :: dude, dvde, dqe, vout, uout, eka
        REAL(kind = mp), ALLOCATABLE, DIMENSION(:,:) :: dsu, dne, dnq, dsv, dse, dsq
        REAL(kind = mp), ALLOCATABLE, DIMENSION(:,:) :: dnu, dnv
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: npos
    end type
    CONTAINS
    SUBROUTINE alloc_working(object, ns, nn, itm)
    implicit none
    type(riv_w_var), intent(inout) :: object
    integer, intent(in) :: ns, nn, itm
    INTEGER :: status, tmpval
    IF (ns.gt.nn) THEN
        ALLOCATE(object%isw(ns), STAT = status)
        tmpval = ns
    ELSE
        ALLOCATE(object%isw(nn), STAT = status)
        tmpval = nn
    END IF
    ALLOCATE(object%qprms(itm*10), STAT = status)
    ALLOCATE(object%dea(ns), STAT = status)
    ALLOCATE(object%eav(ns), STAT = status)
    ALLOCATE(object%am(tmpval), STAT = status)
    ALLOCATE(object%bm(tmpval), STAT = status)
    ALLOCATE(object%ccm(tmpval), STAT = status)
    ALLOCATE(object%dm(tmpval), STAT = status)
    ALLOCATE(object%em(tmpval), STAT = status)
    ALLOCATE(object%qp(ns), STAT = status)
    !		ALLOCATE(u(ns, nn), STAT = status)
    !		ALLOCATE(v(ns, nn), STAT = status)
    !		ALLOCATE(e(ns, nn), STAT = status)
    ALLOCATE(object%dqsde(ns, nn), STAT = status)
    ALLOCATE(object%dqnde(ns, nn), STAT = status)
    ALLOCATE(object%qv(ns, nn), STAT = status)
    ALLOCATE(object%qu(ns, nn), STAT = status)
    ALLOCATE(object%dq(ns, nn), STAT = status)
    ALLOCATE(object%de(ns, nn), STAT = status)
    ALLOCATE(object%dude(ns, nn), STAT = status)
    ALLOCATE(object%dvde(ns, nn), STAT = status)
    ALLOCATE(object%dqe(ns, nn), STAT = status)
    ALLOCATE(object%vout(ns, nn), STAT = status)
    ALLOCATE(object%uout(ns, nn), STAT = STATUS)
    ALLOCATE(object%eka(ns, nn), STAT = status)
    ALLOCATE(object%dsu(ns, nn), STAT = status)
    ALLOCATE(object%dne(ns, nn), STAT = status)
    ALLOCATE(object%dnq(ns, nn), STAT = status)
    ALLOCATE(object%dsv(ns, nn), STAT = status)
    ALLOCATE(object%dse(ns, nn), STAT = status)
    ALLOCATE(object%dsq(ns, nn), STAT = status)
    ALLOCATE(object%dnu(ns, nn), STAT = status)
    ALLOCATE(object%dnv(ns, nn), STAT = status)

    ALLOCATE(object%npos(ns, 2), STAT = status)

    END SUBROUTINE alloc_working
    
    SUBROUTINE dealloc_working(object)
    implicit none
    type(riv_w_var), intent(inout) :: object 
    INTEGER :: status
    DEALLOCATE(object%isw, STAT = status)
    DEALLOCATE(object%qprms, STAT = status)
    DEALLOCATE(object%dea, STAT = status)
    DEALLOCATE(object%eav, STAT = status)
    DEALLOCATE(object%am, STAT = status)
    DEALLOCATE(object%bm, STAT = status)
    DEALLOCATE(object%ccm, STAT = status)
    DEALLOCATE(object%dm, STAT = status)
    DEALLOCATE(object%em, STAT = status)
    DEALLOCATE(object%qp, STAT = status)
    !		DEALLOCATE(u, STAT = status)
    !		DEALLOCATE(v, STAT = status)
    !		DEALLOCATE(e, STAT = status)
    DEALLOCATE(object%dqsde, STAT = status)
    DEALLOCATE(object%dqnde, STAT = status)
    DEALLOCATE(object%qv, STAT = status)
    DEALLOCATE(object%qu, STAT = status)
    DEALLOCATE(object%dq, STAT = status)
    DEALLOCATE(object%de, STAT = status)
    DEALLOCATE(object%dude, STAT = status)
    DEALLOCATE(object%dvde, STAT = status)
    DEALLOCATE(object%dqe, STAT = status)
    DEALLOCATE(object%vout, STAT = status)
    DEALLOCATE(object%uout, STAT = status)
    DEALLOCATE(object%eka, STAT = status)
    DEALLOCATE(object%dsu, STAT = status)
    DEALLOCATE(object%dne, STAT = status)
    DEALLOCATE(object%dnq, STAT = status)
    DEALLOCATE(object%dsv, STAT = status)
    DEALLOCATE(object%dse, STAT = status)
    DEALLOCATE(object%dsq, STAT = status)
    DEALLOCATE(object%dnu, STAT = status)
    DEALLOCATE(object%dnv, STAT = status)

    DEALLOCATE(object%npos, STAT = status)

    END SUBROUTINE dealloc_working

    END MODULE RivVarWMod2
