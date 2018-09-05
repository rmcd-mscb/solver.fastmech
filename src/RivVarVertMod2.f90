    MODULE RivVarVertMod2
    USE RivVarMod2
    USE CalcCond2
    IMPLICIT NONE
    type rivvarvert
        REAL(kind = mp), ALLOCATABLE, DIMENSION(:,:,:) :: ecoef, f1, f2, f3
        REAL(kind = mp), ALLOCATABLE, DIMENSION(:,:) :: zeta0, fone, ftwo, theta, ustr2, rs, dum1a
        REAL(kind = mp), ALLOCATABLE, DIMENSION(:) :: zeta, dzeta
        !	REAL, ALLOCATABLE, DIMENSION(:,:) :: E_Corr, Strmln_Corr

        REAL(kind = mp) :: PI, BETA, sczeta
        REAL(kind = mp) :: finc, finc1, finc2
        REAL(kind = mp) :: dzp, dth, scn, dthds, dthdn
        REAL(kind = mp) :: curv, curvn, curvs
        REAL(kind = mp) :: taux, vs, uu1, vv1, vv2, uu2, g1, g2
        REAL(kind = mp) :: zo
    end type
    CONTAINS
    SUBROUTINE alloc_vert(this, rvo)
    implicit none
    type(rivvarvert),  intent(inout) :: this
    type(rivvar), intent(inout) :: rvo
    INTEGER :: status, ns, nn, nz
    if(rvo%nz <= 1) then
        rvo%nz = 11
        nz = 11
    endif
    ns = rvo%ns
    nn = rvo%nn
    ALLOCATE(this%ecoef(ns, nn, nz), STAT = status)
    ALLOCATE(this%f1(ns, nn, nz), STAT = status)
    ALLOCATE(this%f2(ns, nn, nz), STAT = status)
    ALLOCATE(this%f3(ns, nn, nz), STAT = status)

    ALLOCATE(this%zeta0(ns, nn), STAT = status)
    ALLOCATE(this%fone(ns, nn), STAT = status)
    ALLOCATE(this%ftwo(ns, nn), STAT = status)
    ALLOCATE(this%theta(ns, nn), STAT = status)
    ALLOCATE(this%ustr2(ns, nn), STAT = status)
    ALLOCATE(this%rs(ns, nn), STAT = status)
    Allocate(this%dum1a(ns, nn), Stat = status)
    !		ALLOCATE(E_Corr(ns, nn), STAT = status)
    !		ALLOCATE(Strmln_Corr(ns, nn), STAT = status)

    ALLOCATE(this%zeta(nz), STAT = status)
    ALLOCATE(this%dzeta(nz), STAT = status)
    END SUBROUTINE alloc_vert

    SUBROUTINE dealloc_vert(this)
    implicit none
    type(rivvarvert), intent(inout) :: this
    INTEGER :: status
    DEALLOCATE(this%ecoef, STAT = status)
    DEALLOCATE(this%f1, STAT = status)
    DEALLOCATE(this%f2, STAT = status)
    DEALLOCATE(this%f3, STAT = status)

    DEALLOCATE(this%zeta0, STAT = status)
    DEALLOCATE(this%fone, STAT = status)
    DEALLOCATE(this%ftwo, STAT = status)
    DEALLOCATE(this%theta, STAT = status)
    DEALLOCATE(this%ustr2, STAT = status)
    DEALLOCATE(this%rs, STAT = status)
    DEALLOCATE(this%dum1a, STAT = status)
    !		DEALLOCATE(E_Corr, STAT = status)
    !		DEALLOCATE(Strmln_Corr, STAT = status)

    DEALLOCATE(this%zeta, STAT = status)
    DEALLOCATE(this%dzeta, STAT = status)

    END SUBROUTINE dealloc_vert
    END MODULE