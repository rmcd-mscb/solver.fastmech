    MODULE UInitMod2
    USE RivVarMod2
    USE CalcCond2
    IMPLICIT NONE

    CONTAINS
    subroutine uinit(rvo, cco, u,v,hav,w,eta,q,ibc)
    IMPLICIT NONE
    type(rivvar), intent(inout) :: rvo
    type(calccond), intent(inout) :: cco
    REAL(kind = mp), DIMENSION(:, :), INTENT(INOUT) :: u, v
    REAL(kind = mp), DIMENSION(:, :), INTENT(IN) :: eta
    INTEGER, DIMENSION(:, :), INTENT(IN) :: ibc
    REAL(kind = mp), DIMENSION(:), INTENT(IN) ::w, hav
    REAL(kind = mp), INTENT(IN) :: q
    INTEGER :: i, j, ibtot, count
    REAL(kind = mp) :: deln, depave, twidth, ubar, qpred
    integer :: ns, nn
    ns = rvo%ns
    nn = rvo%nn
    do 30 i=1,ns
        deln=w(i)/(nn-1)
        depave=0.
        ibtot=0
        count = 0
        do 25 j=1,nn
            if(ibc(i,j).ne.0) then
                ibtot=ibtot+1
                depave=depave+GetDepth(i,j,hav,eta)
            endif
25      continue
        if(ibtot.eq.0) then
            rvo%errorcode = -10
            WRITE(6,*) 'At least one row of nodes are completely dry - '
            WRITE(6,*) 'check water surface boundary condition'
            return
        end if
        depave=depave/(ibtot)
        twidth=deln*(ibtot)
        ubar=q/(twidth*depave)
        qpred=0.
        count = 0
        do 27 j=1,nn
            if(ibc(i,j).ne.0) then

                !		  u(i,j)=ubar*(hav(i)-eta(i,j))/((j**5.)*depave)
                if(cco%vbc == 1 .and. i == 1) then
                    u(i,j) = ubar*GetVelocityWieght(cco, count,ibtot)
                else
                    u(i,j)=ubar*GetDepthWieght(cco, i,j,hav,eta,depave)
                endif
                !		  u(i,j)=ubar
                qpred=qpred+u(i,j)*GetDepth(i,j,hav,eta)*deln
                count = count+1
            else
                u(i,j)=0.
            endif
            !        if(i.eq.1) then
            !        v(1,j)=-1.*u(1,j)*tand(30.)
            !        else
            v(i,j)=0.
            !        endif
27      continue
        count = 0
        do 28 j=1,nn
            u(i,j)=u(i,j)*q/qpred
            if(i.eq.1) then
                if(cco%vbc == 1 .and. ibc(i,j).ne.0) then
                    v(1,j) = u(1,j)*tan(GetAngleWieght(cco, count, ibtot))
                    count = count+1
                else
                    v(1,j)=u(1,j)*dtan(cco%vac)
                endif
            else
                v(i,j)=0.
            endif
28      continue
30  continue
    return
    END SUBROUTINE uinit

    subroutine VarDischUInit(rvo, cco, u,v,w,eta,e,ibc, q, stageChange)
    IMPLICIT NONE
    type(rivvar), intent(inout) :: rvo
    type(calccond), intent(inout) :: cco
    REAL(kind = mp), DIMENSION(:, :), INTENT(INOUT) :: u, v, e
    REAL(kind = mp), DIMENSION(:), INTENT(IN) ::w
    REAL(kind = mp), DIMENSION(:, :), INTENT(IN) :: eta
    INTEGER, DIMENSION(:, :), INTENT(INOUT) :: ibc
    REAL(kind = mp), INTENT(IN) :: stageChange
    REAL(kind = mp), INTENT(IN) :: q

    INTEGER :: i, j, ibtot, count
    REAL(kind = mp) :: deln, depave, twidth, ubar, qpred, depth
    integer :: ns, nn
    ns = rvo%ns
    nn = rvo%nn

    do 30 i=1,1
        deln=w(i)/(nn-1)
        depave=0.
        ibtot=0
        count = 0
        do 25 j=1,nn
            !                e(i,j) = e(i,j)
            !                e(i,j) = e(i,j) + stageChange
            depth = e(i,j) - eta(i,j)
            if(depth > cco%hmin) then
                ibc(i,j) = -1
            else
                ibc(i,j) = 0
            endif
            if(ibc(i,j).ne.0) then
                ibtot=ibtot+1
                depave=depave+depth
            endif
25      continue
        if(ibtot.eq.0) then
            rvo%errorcode = -10
            WRITE(6,*) 'At least one row of nodes are completely dry-'
            WRITE(6,*) 'check water surface boundary condition'
            return
        end if
        depave=depave/(ibtot)
        twidth=deln*(ibtot)
        ubar=q/(twidth*depave)
        qpred=0.
        count = 0
        do 27 j=1,nn
            if(ibc(i,j).ne.0) then
                u(i,j)=ubar*((e(i,j)-eta(i,j))**cco%vdc)/depave
                !		            qpred=qpred+u(i,j)*GetDepth(i,j)*deln
                qpred=qpred+u(i,j)*(e(i,j) - eta(i,j))*deln
                count = count+1
            else
                u(i,j)=0.
            endif
            !        if(i.eq.1) then
            !           v(1,j)=-1.*u(1,j)*tand(30.)
            !        else
            v(i,j)=0.
            !        endif
27      continue
        do 28 j=1,nn
            u(i,j)=u(i,j)*q/qpred
            if(i.eq.1) then
                v(1,j)=-1.*u(1,j)*dtan(cco%vac)
            else
                v(i,j)=0.
            endif
28      continue
30  continue
    return
    END SUBROUTINE VarDischUInit

    FUNCTION GetDepth(i,j,hav,eta) result(r)
    INTEGER, INTENT(IN)::i, j
    real(kind=mp), dimension(:), intent(in) :: hav
    real(kind=mp), dimension(:,:), intent(in) :: eta
    real(kind=mp) :: r
    r = hav(i)-eta(i,j)
    END FUNCTION

    FUNCTION GetDepthWieght(cco, i, j, hav,eta,depave) result(r)
    type(calccond), intent(in) :: cco
    INTEGER, INTENT(IN):: i, j
    real(kind=mp), dimension(:), intent(in) :: hav
    real(kind=mp), dimension(:,:), intent(in) :: eta
    real(kind=mp) :: r
    REAL(kind = mp), INTENT(IN)::depave
    r = ((hav(i)-eta(i,j))**cco%vdc)/depave
    END FUNCTION

    FUNCTION GetVelocityWieght(cco, i,nn) result(r)
    type(calccond), intent(in) :: cco
    INTEGER, INTENT(IN) :: i,nn
    INTEGER :: count
    REAL(kind = mp) :: val
    real(kind=mp) :: r
    val = real(i)/real(nn)
    count = 1
    DO WHILE(cco%vbcdist(count) <= val)
        count=count+1
    END DO
    count = count -1
    r = cco%vbcvel(count) & !y
    + (val-cco%vbcdist(count)) & !x-xo
    * (((cco%vbcvel(count+1)-cco%vbcvel(count)) & !slope
    / (cco%vbcdist(count+1)-cco%vbcdist(count))))
    END FUNCTION

    FUNCTION GetAngleWieght(cco,i,nn) result(r)
    type(calccond), intent(in) :: cco
    INTEGER, INTENT(IN) :: i, nn
    INTEGER :: count
    REAL(kind = mp) :: val
    real(kind = mp) :: r
    val = real(i)/real(nn)
    count = count - 1
    DO WHILE (cco%vbcdist(count) <= val)
        count=count+1
    END DO
    count = count -1
    r = cco%vbcang(count) &    !y
    + (val-cco%vbcdist(count)) & !x-xo
    * (((cco%vbcang(count+1)-cco%vbcang(count)) &   !slope
    / (cco%vbcdist(count+1)-cco%vbcdist(count))))
    END FUNCTION

    END MODULE UInitMod2
