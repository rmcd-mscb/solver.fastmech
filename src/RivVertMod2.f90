    MODULE RivVertMod2
    USE RivVarVertMod2
    use rivvarmod2
    USE CalcCond2
    IMPLICIT NONE
    CONTAINS

    SUBROUTINE ZOTOCDTWO(rvo, cco)
    IMPLICIT NONE
    type(rivvar), intent(inout) :: rvo
    type(calccond), intent(in) :: cco
    INTEGER :: i,j, ns, nn
    REAL(kind = mp) :: tmp1, tmp2
    REAL(kind = mp) :: vkc = 0.4
    ns = rvo%ns
    nn = rvo%nn
    do i = 1,ns
        do j = 1,nn
            tmp1 = log(rvo%hl(i,j)/rvo%znaught(i,j)) - 1
            !	            tmp2 =
            rvo%cd(i,j) = ((1/rvo%vkc)*tmp1)**(-2)
            if(rvo%cd(i,j) < cco%cdmin)then
                rvo%cd(i,j) = cco%cdmin
            endif
            if(rvo%cd(i,j) > cco%cdmax) then
                rvo%cd(i,j) = cco%cdmax
            endif
            rvo%totcd(i,j) = rvo%cd(i,j)

        enddo
    enddo
    END SUBROUTINE



    SUBROUTINE Z0TOCD(rvo, rvvo, cco)
    implicit none
    type(rivvar), intent(inout) :: rvo
    type(rivvarvert), intent(inout) :: rvvo
    type(calccond), intent(in) :: cco
    INTEGER :: i, j, k, ns, nn, nz
    REAL(kind = mp) :: vkc = 0.4
    ns = rvo%ns
    nn = rvo%nn
    nz = rvo%nz
    !	cdmin = 0.001
    !	cdmax = 0.1
    call alloc_vert(rvvo, rvo)
    rvvo%pi=acos(-1.)
    rvvo%beta=6.24
    DO i=1,ns
        DO j=1,nn
            if(rvo%hl(i,j) > 0) then
                rvvo%zeta0(i,j)= rvo%znaught(i,j)/rvo%hl(i,j)
            else
                rvvo%zeta0(i,j)= rvo%znaught(i,j)/.01
            endif
            rvvo%sczeta=(log(1./rvvo%zeta0(i,j)))/(nz-1)
            DO k=1,nz
                rvvo%zeta(k)=rvvo%zeta0(i,j)*exp((k-1)*rvvo%sczeta)
                if(rvvo%zeta(k).lt.0.2) then
                    rvvo%ecoef(i,j,k)=rvo%vkc*rvvo%zeta(k)*(1.-rvvo%zeta(k))
                else
                    rvvo%ecoef(i,j,k)=rvo%vkc/rvvo%beta
                endif
                if(k.eq.1) then
                    rvvo%f1(i,j,k)=0.
                    rvvo%fone(i,j)=0.
                    rvvo%f3(i,j,k)=0.
                else
                    rvvo%dzeta(k)=rvvo%zeta(k)-rvvo%zeta(k-1)
                    rvvo%finc1=(1.-rvvo%zeta(k))/rvvo%ecoef(i,j,k)
                    rvvo%finc2=(1.-rvvo%zeta(k-1))/rvvo%ecoef(i,j,k-1)
                    rvvo%finc=(rvvo%finc1+rvvo%finc2)
                    rvvo%dzp=0.5*rvvo%dzeta(k)
                    rvvo%f1(i,j,k)=rvvo%f1(i,j,k-1)+rvvo%finc*rvvo%dzp
                    rvvo%fone(i,j)=rvvo%fone(i,j)+(rvvo%f1(i,j,k)+rvvo%f1(i,j,k-1))*rvvo%dzp
                    rvvo%f3(i,j,k)=rvvo%f3(i,j,k-1)+(rvvo%f1(i,j,k)**2.+rvvo%f1(i,j,k-1)**2.)*rvvo%dzp
                endif
            ENDDO

            rvo%cd(i,j) = 1/(rvvo%fone(i,j)**2.)
            if(rvo%cd(i,j) > cco%cdmax) then
                rvo%cd(i,j) = cco%cdmax
            endif
            if(rvo%cd(i,j) < cco%cdmin) then
                rvo%cd(i,j) = cco%cdmin
            endif
            if(rvo%ibc(i,j).eq.0) then
                rvo%cd(i,j) = 0
            endif
            rvo%totcd(i,j) = rvo%cd(i,j)

            !		 write(6,*) i, j, fone(i,j), cd(i,j)
        ENDDO
    ENDDO

    call dealloc_vert(rvvo)


    END SUBROUTINE

    SUBROUTINE vert(rvo, rvvo, cco)
    implicit none
    type(rivvar), intent(inout) :: rvo
    type(rivvarvert), intent(inout) :: rvvo
    type(calccond), intent(in) :: cco
    INTEGER :: i, j, k, ismoo, ns, nn, nz
    REAL(kind = mp) :: tmp, tmp1, tmp2
    REAL(kind = mp) :: vkc = 0.4
    
    call alloc_vert(rvvo, rvo)
    rvvo%pi=acos(-1.)
    rvvo%beta=6.24
    DO i=1,ns
        DO j=1,nn
            if(cco%transeqtype == 4) then !pan
                tmp2 = (-1)*((rvo%vkc/sqrt(rvo%totcd(i,j)))+1)
                rvvo%zeta0(i,j) = exp(tmp2)
            else
                if(cco%roughnesstype == 1) then
                    if(rvo%hl(i,j) > 0) then
                        rvvo%zeta0(i,j)= rvo%znaught(i,j)/rvo%hl(i,j)
                    else
                        rvvo%zeta0(i,j)= .001
                    endif

                else
                    !				if(hl(i,j) > 0) then
                    tmp2 = (-1)*((rvo%vkc/sqrt(rvo%totcd(i,j)))+1)
                    rvvo%zeta0(i,j) = exp(tmp2)
                    !				else
                    !					zeta0(i,j) = 0.001
                    !				endif
                endif
            endif !pan
            !		zeta0 = 0.0146
            rvvo%sczeta=(log(1./rvvo%zeta0(i,j)))/(nz-1)

            DO k=1,nz
                rvvo%zeta(k)=rvvo%zeta0(i,j)*exp((k-1)*rvvo%sczeta)
                rvo%zz(i, j, k) = rvo%eta(i,j) + (rvo%hl(i,j)*rvvo%zeta(k))
                if(rvvo%zeta(k).lt.0.2) then
                    rvvo%ecoef(i,j,k)=rvo%vkc*rvvo%zeta(k)*(1.-rvvo%zeta(k))
                else
                    rvvo%ecoef(i,j,k)=rvo%vkc/rvvo%beta
                endif
                if(k.eq.1) then
                    rvvo%f1(i,j,k)=0.
                    rvvo%fone(i,j)=0.
                    rvvo%f3(i,j,k)=0.
                else
                    rvvo%dzeta(k)=rvvo%zeta(k)-rvvo%zeta(k-1)
                    rvvo%finc1=(1.-rvvo%zeta(k))/rvvo%ecoef(i,j,k)
                    rvvo%finc2=(1.-rvvo%zeta(k-1))/rvvo%ecoef(i,j,k-1)
                    rvvo%finc=(rvvo%finc1+rvvo%finc2)
                    rvvo%dzp=0.5*rvvo%dzeta(k)
                    rvvo%f1(i,j,k)=rvvo%f1(i,j,k-1)+rvvo%finc*rvvo%dzp
                    rvvo%fone(i,j)=rvvo%fone(i,j)+(rvvo%f1(i,j,k)+rvvo%f1(i,j,k-1))*rvvo%dzp
                    rvvo%f3(i,j,k)=rvvo%f3(i,j,k-1)+(rvvo%f1(i,j,k)**2.+rvvo%f1(i,j,k-1)**2.)*rvvo%dzp
                endif
            ENDDO
            DO k=1,nz
                rvvo%f3(i,j,k)=rvvo%f3(i,j,nz)-rvvo%f3(i,j,k)
                if(k.eq.1) then
                    rvvo%f2(i,j,k)=0.
                    rvvo%ftwo(i,j)=0.
                else
                    rvvo%finc1=rvvo%f3(i,j,k)/rvvo%ecoef(i,j,k)
                    rvvo%finc2=rvvo%f3(i,j,k-1)/rvvo%ecoef(i,j,k-1)
                    rvvo%f2(i,j,k)=rvvo%f2(i,j,k-1)+(rvvo%finc1+rvvo%finc2)*0.5*rvvo%dzeta(k)
                    rvvo%ftwo(i,j)=rvvo%ftwo(i,j)+(rvvo%f2(i,j,k)+rvvo%f2(i,j,k-1))*0.5*rvvo%dzeta(k)
                endif
            ENDDO
            if(rvo%taus(i,j).eq.0.and.rvo%taun(i,j).eq.0) then
                rvvo%theta(i,j)=1000.
            else
                rvvo%theta(i,j)=atan2(rvo%taun(i,j),rvo%taus(i,j))
            endif
            rvvo%ustr2(i,j)=(rvo%taus(i,j)**2.+rvo%taun(i,j)**2.)**0.5
        ENDDO
    ENDDO

    DO i=1,ns
        DO j=1,nn
            if(i.eq.1) then
                rvvo%dth=rvvo%theta(i+1,j)-rvvo%theta(i,j)
                rvvo%scn=1.
            else if(i.eq.ns) then
                rvvo%dth=rvvo%theta(i,j)-rvvo%theta(i-1,j)
                rvvo%scn=1.
            else if(i.ne.1.and.i.ne.ns) then
                if(rvvo%theta(i-1,j).eq.1000.) then
                    rvvo%dth=rvvo%theta(i+1,j)-rvvo%theta(i,j)
                    rvvo%scn=1.
                else if(rvvo%theta(i+1,j).eq.1000.) then
                    rvvo%dth=rvvo%theta(i,j)-rvvo%theta(i-1,j)
                    rvvo%scn=1.
                else
                    rvvo%dth=(rvvo%theta(i+1,j)-rvvo%theta(i-1,j))
                    rvvo%scn=2.
                endif
            endif
            if(rvvo%dth.gt.rvvo%pi) then
                rvvo%dth=rvvo%dth-2.*rvvo%pi
            else if(rvvo%dth.lt.(-1.*rvvo%pi)) then
                rvvo%dth=2.*rvvo%pi+rvvo%dth
            endif
            rvvo%dthds=rvvo%dth/(rvvo%scn*rvo%ds*rvo%rn(i,j))
            if(j.eq.1) then
                rvvo%dth=rvvo%theta(i,j+1)-rvvo%theta(i,j)
                rvvo%scn=1.
            else if(j.eq.nn) then
                rvvo%dth=rvvo%theta(i,j)-rvvo%theta(i,j-1)
                rvvo%scn=1.
            else if(j.ne.1.or.j.ne.nn) then
                if(rvvo%theta(i,j-1).eq.1000.) then
                    rvvo%dth=rvvo%theta(i,j+1)-rvvo%theta(i,j)
                    rvvo%scn=1.
                elseif(rvvo%theta(i,j+1).eq.1000.) then
                    rvvo%dth=rvvo%theta(i,j)-rvvo%theta(i,j-1)
                    rvvo%scn=1.
                else
                    rvvo%dth=(rvvo%theta(i,j+1)-rvvo%theta(i,j-1))
                    rvvo%scn=2.
                endif
            endif
            if(rvvo%dth.gt.rvvo%pi) then
                rvvo%dth=rvvo%dth-2.*rvvo%pi
            elseif(rvvo%dth.lt.(-1.*rvvo%pi)) then
                rvvo%dth=2.*rvvo%pi+rvvo%dth
            endif
            rvvo%dthdn=rvvo%dth/(rvvo%scn*rvo%dn)
            if(rvvo%ustr2(i,j).eq.0.or.rvo%ibc(i,j).ne.-1) then
                rvvo%rs(i,j)=10**8.
            else
                !		Adjusted for Kootenai 7/22/03 jmn & rmcd
                !		Center line and stream line curvature accounted for
                If(cco%CALCQUASI3DRS) Then
                    rvvo%curvs=(rvvo%dthds+(1./rvo%r(i)))*rvo%taus(i,j)/rvvo%ustr2(i,j)
                    rvvo%curvn=rvvo%dthdn*rvo%taun(i,j)/rvvo%ustr2(i,j)
                ELSE
                    rvvo%curvs=((1./rvo%r(i)))*rvo%taus(i,j)/rvvo%ustr2(i,j)
                    rvvo%curvn = 0
                ENDIF
                !         rvvo%curvs=(dthds+(1./r(i)))*taus(i,j)/ustr2(i,j)
                !         rvvo%curvn=dthdn*taun(i,j)/ustr2(i,j)
                rvvo%curv=rvvo%curvs+rvvo%curvn
                if(rvvo%curv.eq.0.) then
                    rvvo%rs(i,j)=10**8.
                else
                    rvvo%rs(i,j)=1./rvvo%curv
                endif
            endif
            if(abs(rvvo%rs(i,j)).lt.cco%minrs) then
                rvvo%rs(i,j)=cco%minrs*dsign(1.0D0,rvvo%rs(i,j))
            endif
            rvvo%taux=rvvo%ustr2(i,j)*rvo%hl(i,j)/(rvo%rn(i,j)*rvvo%rs(i,j))
            rvvo%vs=((rvvo%ftwo(i,j)/rvvo%fone(i,j))-rvvo%f3(i,j,1))
            rvvo%taux=rvvo%taux*rvvo%vs
            DO k=1,nz
                rvvo%uu1=(rvvo%ustr2(i,j)**.5)*cos(rvvo%theta(i,j))*rvvo%f1(i,j,k)
                rvvo%vv1=(rvvo%ustr2(i,j)**.5)*sin(rvvo%theta(i,j))*rvvo%f1(i,j,k)
                rvvo%g1=(rvvo%ustr2(i,j)**.5)*rvo%hl(i,j)/(rvo%rn(i,j)*rvvo%rs(i,j))
                rvvo%g2=((rvvo%ftwo(i,j)/rvvo%fone(i,j))*rvvo%f1(i,j,k))-rvvo%f2(i,j,k)
                rvvo%uu2=-1.*rvvo%g1*rvvo%g2*sin(rvvo%theta(i,j))
                rvvo%vv2=rvvo%g1*rvvo%g2*cos(rvvo%theta(i,j))
                rvo%uz(i,j,k)=rvvo%uu1+rvvo%uu2
                rvo%vz(i,j,k)=rvvo%vv1+rvvo%vv2
            ENDDO
            rvo%taus(i,j)=rvo%taus(i,j)-rvvo%taux*sin(rvvo%theta(i,j))
            rvo%taun(i,j)=rvo%taun(i,j)+rvvo%taux*cos(rvvo%theta(i,j))
            !!        taus(i,j)=taus(i,j)-taux*sin(theta(i,j))
            !        taun(i,j)=taus(i,j)*7.*(hl(i,j)/rs(i,j))
            !        E_Corr(i,j) = atan(7*hl(i,j)/rs(i,j)) *(360./(2.*3.14159))
            !        Strmln_Corr(i,j) = atan2(taun(i,j),taus(i,j)) *(360./(2.*3.14159))
        ENDDO
    ENDDO
    !        return
    IF(cco%RSSMOO.gt.0) THEN
        DO  ISMOO=1,cco%RSSMOO
            DO  I=2,ns-1
                DO  J=1,nn
                    if(j.eq.1) then
                        rvvo%dum1a(I,1)=(rvvo%rs(I-1,1)+rvvo%rs(I,1)+rvvo%rs(I+1,1)+rvvo%rs(I,j+1))/4.
                    elseif(j.eq.nn) then
                        rvvo%dum1a(I,nn)=(rvvo%rs(I-1,nn)+rvvo%rs(I,nn)+rvvo%rs(I+1,nn)+rvvo%rs(i,nn-1))/4.
                    else
                        rvvo%dum1a(I,J)=rvvo%rs(I,J)+cco%sedsmoowght*rvvo%rs(I+1,J)+cco%sedsmoowght*rvvo%rs(I-1,J)
                        rvvo%dum1a(I,J)=rvvo%dum1a(I,J)+cco%sedsmoowght*rvvo%rs(I,J-1)+cco%sedsmoowght*rvvo%rs(I,J+1)
                        rvvo%dum1a(I,J)=rvvo%dum1a(I,J)/5.
                    endif
                ENDDO
            ENDDO
            DO I=2,ns-1
                DO J=1,nn
                    rvvo%rs(I,J)=rvvo%dum1a(I,J)
                ENDDO
            ENDDO
        ENDDO


        do i=1,ns
            do j=1,nn
                if(abs(rvvo%rs(i,j)).lt.cco%minrs) then
                   rvvo%rs(i,j)=cco%minrs*dsign(1.0D0,rvvo%rs(i,j))
                endif
                rvvo%taux=rvvo%ustr2(i,j)*rvo%hl(i,j)/(rvo%rn(i,j)*rvvo%rs(i,j))
                rvvo%vs=((rvvo%ftwo(i,j)/rvvo%fone(i,j))-rvvo%f3(i,j,1))
                rvvo%taux=rvvo%taux*rvvo%vs
                rvo%taus(i,j)=rvo%taus(i,j)-rvvo%taux*sin(rvvo%theta(i,j))
                rvo%taun(i,j)=rvo%taun(i,j)+rvvo%taux*cos(rvvo%theta(i,j))
            ENDDO
        ENDDO

    ENDIF

    END SUBROUTINE vert

    END MODULE RivVertMod2
