MODULE ParamSetMod
USE Support
IMPLICIT NONE
!INCLUDE "netcdf.inc"

INTEGER::nprm, ndiams(4)
!REAL::prmt(4,8), pdiams(4,20), gomfl(3), bcfl(3), cbfl
REAL::tmpprmt(9), tmptime(2)
INTEGER::nodebalindx, indxresist, k
CONTAINS

Subroutine ParamSet()
 !   common/prmtrs/nprm,prmt(4,8),ndiams(4),pdiams(4,20),&
 ! gomfl(3),bcfl(3),cgomfl,cbcfl
	nprm=8
  prmt=0.
  prmt(1,1)=6. !courant time step multiplier
  prmt(1,2)=.2 !dry-wet channel tolerance, m
  prmt(1,3)=.7  !theta, advanced time step weight in Preissman scheme
  prmt(1,4)=100.  !Number of points in time series plots, max=200
  prmt(1,5)=1.e-5  !Mass concentration tolerance level
  prmt(1,6)=50.  !Maximum number iterations, preissman scheme
  prmt(1,7)=.4  !Relaxation coefficient--discharge, preissman scheme
  prmt(1,8)=.6  !Relaxation coefficient--stage, preissman scheme
	dragtype = 1;
	corstpmlt=prmt(1,1)  ! min courant time step multiplier for setting dt used in time stepping
	etol=prmt(1,2)
	thet=prmt(1,3)	! time step weighting factor applied to future step in iterative solution
	itersmx=int(prmt(1,6))
	voltol=prmt(1,5)
	rlxdschg=prmt(1,7)
	rlxstg=prmt(1,8)
  nodebalindx=1   ! run nodebalance during simulation if index==1
  indxresist=0    ! use user specified alluvail n values rather than program calculated   
	do k=1,nprm
	prmt(3,k)=prmt(1,k)
	prmt(2,k)=prmt(1,k)
	end do


end Subroutine ParamSet

!subroutine ParamSetNetCDF(InputFile)
!
!	CHARACTER(*), INTENT(IN) :: InputFile
!	INTEGER :: status, NCID, tmpint
!	status = NF_OPEN(InputFile, NF_NOWRITE, NCID)
!		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
!
!	status = nf_get_att_real(NCID, NF_GLOBAL, "1DPreiss_ModelParams", tmpprmt)
!		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
!	nprm=8
!	  prmt=0.
!	  prmt(1,1)=tmpprmt(1) !courant time step multiplier
!	  prmt(1,2)=tmpprmt(2) !dry-wet channel tolerance, m
!	  prmt(1,3)=tmpprmt(3)  !theta, advanced time step weight in Preissman scheme
!	  prmt(1,4)=tmpprmt(4)  !Number of points in time series plots, max=200
!	  prmt(1,5)=tmpprmt(5) !Mass concentration tolerance level
!	  prmt(1,6)=tmpprmt(6)  !Maximum number iterations, preissman scheme
!	  prmt(1,7)=tmpprmt(7)  !Relaxation coefficient--discharge, preissman scheme
!	  prmt(1,8)=tmpprmt(8)  !Relaxation coefficient--stage, preissman scheme
!	  dragtype = int(tmpprmt(9)) !Drag type variable
!		corstpmlt=prmt(1,1)  ! min courant time step multiplier for setting dt used in time stepping
!		etol=prmt(1,2)
!		thet=prmt(1,3)	! time step weighting factor applied to future step in iterative solution
!		itersmx=int(prmt(1,6))
!		voltol=prmt(1,5)
!		rlxdschg=prmt(1,7)
!		rlxstg=prmt(1,8)
!	  nodebalindx=1   ! run nodebalance during simulation if index==1
!	  indxresist=0    ! use user specified alluvail n values rather than program calculated   
!		do k=1,nprm
!			prmt(3,k)=prmt(1,k)
!			prmt(2,k)=prmt(1,k)
!		end do
!
!	status = nf_get_att_real(NCID, NF_GLOBAL, "1DPreiss_TimeParams", tmptime)
!		IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
!	simstrt = tmptime(1)
!	simend = tmptime(2)
!	status = nf_close(NCID)
!
!end Subroutine ParamSetNetCDF
!
!SUBROUTINE HANDLE_ERR(STATUS)
!	INCLUDE "netcdf.inc"
!	INTEGER STATUS
!	IF (STATUS .NE. NF_NOERR) THEN
!		PRINT *, NF_STRERROR(STATUS)
!		STOP 'Stopped'
!	ENDIF
!END SUBROUTINE HANDLE_ERR

end Module ParamSetMod