    Subroutine tsbc(nc,nloc,t,q)
    use Support
    IMPLICIT NONE                             ! SRC
    INTEGER, PARAMETER :: mp = KIND(1.0D0)    ! SRC
    INTEGER        NS                         ! SRC  xFTSBoundCond.f90(4)
    INTEGER        NC                         ! SRC  xFTSBoundCond.f90(4)
    INTEGER        NLOC                       ! SRC  xFTSBoundCond.f90(4)
    INTEGER        NE                         ! SRC  xFTSBoundCond.f90(5)
    REAL(KIND=mp)  T                          ! SRC  xFTSBoundCond.f90(7)
    REAL(KIND=mp)  Q                          ! SRC  xFTSBoundCond.f90(8)
    INTEGER        N                          ! SRC  xFTSBoundCond.f90(14)
    REAL(KIND=mp)  R                          ! SRC  xFTSBoundCond.f90(16)
    ! times t are in days
    ns=tsbcloc(nc,2*nloc-1)
    ne=tsbcloc(nc,2*nloc)
    ! 	print *,ns,ne,t
    if(t<=tmsrbc(ns,1))then
        q=tmsrbc(ns,2)
        return
    elseif(t>=tmsrbc(ne,1))then
        q=tmsrbc(ne,2)
        return
    else
        do n=ns+1,ne
            if(t<=tmsrbc(n,1)) then
                r=(t-tmsrbc(n-1,1))/(tmsrbc(n,1)-tmsrbc(n-1,1))
                q=tmsrbc(n-1,2)*(1.-r)+tmsrbc(n,2)*r
                return
            endif
        enddo
    endif
    end Subroutine tsbc
