	MODULE RivVarVertMod
	USE RivVarMod
	USE CalcCond
	IMPLICIT NONE
	REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ecoef, f1, f2, f3
	REAL, ALLOCATABLE, DIMENSION(:,:) :: zeta0, fone, ftwo, theta, ustr2, rs, dum1a
	REAL, ALLOCATABLE, DIMENSION(:) :: zeta, dzeta
!	REAL, ALLOCATABLE, DIMENSION(:,:) :: E_Corr, Strmln_Corr
	
	REAL :: PI, VKC, BETA, sczeta
	REAL :: finc, finc1, finc2
	REAL :: dzp, dth, scn, dthds, dthdn
	REAL :: curv, curvn, curvs
	REAL :: taux, vs, uu1, vv1, vv2, uu2, g1, g2
	REAL :: zo
	CONTAINS
	SUBROUTINE alloc_vert()
		INTEGER :: status
		if(nz <= 1) then
		    nz = 11
		endif
		ALLOCATE(ecoef(ns, nn, nz), STAT = status)
		ALLOCATE(f1(ns, nn, nz), STAT = status)
		ALLOCATE(f2(ns, nn, nz), STAT = status)
		ALLOCATE(f3(ns, nn, nz), STAT = status)

		ALLOCATE(zeta0(ns, nn), STAT = status)
		ALLOCATE(fone(ns, nn), STAT = status)
		ALLOCATE(ftwo(ns, nn), STAT = status)
		ALLOCATE(theta(ns, nn), STAT = status)
		ALLOCATE(ustr2(ns, nn), STAT = status)
		ALLOCATE(rs(ns, nn), STAT = status)
		Allocate(dum1a(ns, nn), Stat = status)
!		ALLOCATE(E_Corr(ns, nn), STAT = status)
!		ALLOCATE(Strmln_Corr(ns, nn), STAT = status)

		ALLOCATE(zeta(nz), STAT = status)
		ALLOCATE(dzeta(nz), STAT = status)
	END SUBROUTINE alloc_vert

	SUBROUTINE dealloc_vert()
		INTEGER :: status
		DEALLOCATE(ecoef, STAT = status)
		DEALLOCATE(f1, STAT = status)
		DEALLOCATE(f2, STAT = status)
		DEALLOCATE(f3, STAT = status)

		DEALLOCATE(zeta0, STAT = status)
		DEALLOCATE(fone, STAT = status)
		DEALLOCATE(ftwo, STAT = status)
		DEALLOCATE(theta, STAT = status)
		DEALLOCATE(ustr2, STAT = status)
		DEALLOCATE(rs, STAT = status)
        DEALLOCATE(dum1a, STAT = status)
!		DEALLOCATE(E_Corr, STAT = status)
!		DEALLOCATE(Strmln_Corr, STAT = status)

		DEALLOCATE(zeta, STAT = status)
		DEALLOCATE(dzeta, STAT = status)

	END SUBROUTINE dealloc_vert
END MODULE