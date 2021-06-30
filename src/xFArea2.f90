    SUBROUTINE FAREA(ns,e)
    ! xk == conveyance, must be squared to => Sf
    use Support
    IMPLICIT NONE                             ! SRC
    INTEGER, PARAMETER :: mp = KIND(1.0D0)    ! SRC
    INTEGER        NPS                        ! SRC  xFArea2.f90(4)
    INTEGER        NS                         ! SRC  xFArea2.f90(5)
    INTEGER        NPF                        ! SRC  xFArea2.f90(6)
    REAL(KIND=mp)  A                          ! SRC  xFArea2.f90(7)
    REAL(KIND=mp)  AM                         ! SRC  xFArea2.f90(8)
    REAL(KIND=mp)  XK                         ! SRC  xFArea2.f90(9)
    REAL(KIND=mp)  TW                         ! SRC  xFArea2.f90(10)
    REAL(KIND=mp)  WP                         ! SRC  xFArea2.f90(11)
    INTEGER        NFRST                      ! SRC  xFArea2.f90(12)
    REAL(KIND=mp)  Y1                         ! SRC  xFArea2.f90(13)
    REAL(KIND=mp)  YFRST                      ! SRC  xFArea2.f90(14)
    REAL(KIND=mp)  Z1                         ! SRC  xFArea2.f90(15)
    REAL(KIND=mp)  XN1                        ! SRC  xFArea2.f90(17)
    REAL(KIND=mp)  E                          ! SRC  xFArea2.f90(19)
    REAL(KIND=mp)  YLST                       ! SRC  xFArea2.f90(23)
    INTEGER        I                          ! SRC  xFArea2.f90(24)
    REAL(KIND=mp)  Y2                         ! SRC  xFArea2.f90(25)
    REAL(KIND=mp)  Z2                         ! SRC  xFArea2.f90(26)
    REAL(KIND=mp)  XN2                        ! SRC  xFArea2.f90(28)
    REAL(KIND=mp)  TWI                        ! SRC  xFArea2.f90(38)
    REAL(KIND=mp)  DA                         ! SRC  xFArea2.f90(39)
    REAL(KIND=mp)  DZ                         ! SRC  xFArea2.f90(47)
    REAL(KIND=mp)  DAT                        ! SRC  xFArea2.f90(52)
    REAL(KIND=mp)  HR                         ! SRC  xFArea2.f90(108)
    REAL(KIND=mp)  HD                         ! SRC  xFArea2.f90(109)
    REAL(KIND=mp)  CTRD                       ! SRC  xFArea2.f90(119)
    REAL(KIND=mp)  R                          ! SRC  xFArea2.f90(120)
    REAL(KIND=mp)  RM                         ! SRC  xFArea2.f90(121)
    nps=1
    if(ns>1)nps=nscpte(ns-1)+1
    npf=nscpte(ns)
    a=0.
    am=0.
    xk=0.
    tw=0.
    wp=0.
    nfrst=0
    y1=xspts(nps,1)
    yfrst=y1
    z1=xspts(nps,2)
    if(dragtype == 0) then
        xn1=xspts(nps,3)
    else
        CALL CnvrtCD(nps,e,xn1)
    endif

    !	xn1=xspts(nps,3)
    ylst=xspts(npf,1)
    do i=nps+1,npf
        y2=xspts(i,1)
        z2=xspts(i,2)
        if(dragtype == 0) then
            xn2=xspts(i,3)
        else
            CALL CnvrtCD(i, e, xn2)
        endif

        !		xn2=xspts(i,3)
        !		print *,i,y1,z1,xn1,y2,z2,xn2
        if((z1<e.or.z2<e).and.y2>y1)then
            !		print *,i,y1,z1,xn1,y2,z2,xn2
            if(z1==z2)then  ! segment of channel bottom is flat
                twi=y2-y1
                da=(e-z1)*twi
                a=a+da
                am=am+da*(y1+y2)/2.
                xk=xk+da**1.66667/xn2/twi**.66667
                tw=tw+twi
                wp=wp+twi
            elseif(z1<z2)then
                twi=y2-y1
                dz=z2-z1
                if(e<z2)then
                    twi=twi*(e-z1)/dz
                    dz=e-z1
                endif
                dat=twi*dz/2.
                a=a+dat
                am=am+dat*(y1+twi/3.)
                tw=tw+twi
                wp=wp+sqrt(dz*dz+twi*twi)
                if(e<=z2)then
                    xk=xk+dat**1.66667/xn2/twi**.66667
                else !e>z2
                    da=twi*(e-z2)
                    a=a+da
                    am=am+da*(y1+y2)/2.
                    xk=xk+(dat+da)**1.66667/xn2/twi**.66667
                endif
            else !z2<z1
                twi=y2-y1
                dz=z1-z2
                if(e<z1)then
                    twi=twi*(e-z2)/dz
                    dz=e-z2
                endif
                dat=twi*dz/2.
                a=a+dat
                am=am+dat*(y2-twi/3.)
                tw=tw+twi
                wp=wp+sqrt(dz*dz+twi*twi)
                if(e<=z1)then
                    xk=xk+dat**1.66667/xn2/twi**.66667
                else !e>z1
                    da=twi*(e-z1)
                    a=a+da
                    am=am+da*(y1+y2)/2.
                    xk=xk+(dat+da)**1.66667/xn2/twi**.66667
                endif
            endif
        endif  ! this i contributes to area
        ! locate the rt and lt ends of ws
        if(nfrst==0.and.z2<=z1.and.(e>=z2.and.e<=z1))then
            nfrst=1
            if(z1==z2)then
                yfrst=y2
            else
                yfrst=y1+(z1-e)*(y2-y1)/(z1-z2)
            endif
        endif
        if(z1<=z2.and.(e>=z1.and.e<=z2))then
            if(z1==z2)then
                ylst=y1
            else
                ylst=y1+(e-z1)*(y2-y1)/(z2-z1)
            endif
        endif
        y1=y2
        z1=z2
        xn1=xn2
    enddo
    !	print *,yfrst,ylst
    hr=0.
    hd=0.
    curxsc(ns,1)=e
    curxsc(ns,2)=a
    curxsc(ns,6)=xk
    if(tw>0.)hd=a/tw
    if(wp>0.)hr=a/wp
    curxsc(ns,3)=hr
    curxsc(ns,4)=hd
    curxsc(ns,5)=tw
    !	curxsc(ns,7)=q
    if(a>0.)ctrd=am/a
    r=ctrd/xsloc(ns,6)
    rm=1.-r
    plnxsc(ns,1)=r*xsloc(ns,3)+rm*xsloc(ns,1) !ctrx
    plnxsc(ns,2)=r*xsloc(ns,4)+rm*xsloc(ns,2) !ctry
    r=yfrst/xsloc(ns,6)
    rm=1.-r
    plnxsc(ns,3)=r*xsloc(ns,3)+rm*xsloc(ns,1) !rtwsx
    plnxsc(ns,4)=r*xsloc(ns,4)+rm*xsloc(ns,2) !rtwsy
    r=ylst/xsloc(ns,6)
    rm=1.-r
    plnxsc(ns,5)=r*xsloc(ns,3)+rm*xsloc(ns,1) !ltwsx
    plnxsc(ns,6)=r*xsloc(ns,4)+rm*xsloc(ns,2) !ltwsy
    end SUBROUTINE FAREA
    SUBROUTINE AREAOnly(ns,e,a)
    ! xk == conveyance, must be squared to => Sf
    use Support
    IMPLICIT NONE                             ! SRC
    INTEGER, PARAMETER :: mp = KIND(1.0D0)    ! SRC
    INTEGER        NPS                        ! SRC  xFArea2.f90(136)
    INTEGER        NS                         ! SRC  xFArea2.f90(137)
    INTEGER        NPF                        ! SRC  xFArea2.f90(138)
    REAL(KIND=mp)  A                          ! SRC  xFArea2.f90(140)
    REAL(KIND=mp)  Y1                         ! SRC  xFArea2.f90(141)
    REAL(KIND=mp)  Z1                         ! SRC  xFArea2.f90(142)
    REAL(KIND=mp)  XN1                        ! SRC  xFArea2.f90(144)
    REAL(KIND=mp)  E                          ! SRC  xFArea2.f90(146)
    INTEGER        I                          ! SRC  xFArea2.f90(150)
    REAL(KIND=mp)  Y2                         ! SRC  xFArea2.f90(151)
    REAL(KIND=mp)  Z2                         ! SRC  xFArea2.f90(152)
    REAL(KIND=mp)  XN2                        ! SRC  xFArea2.f90(154)
    REAL(KIND=mp)  TWI                        ! SRC  xFArea2.f90(164)
    REAL(KIND=mp)  DA                         ! SRC  xFArea2.f90(165)
    REAL(KIND=mp)  DZ                         ! SRC  xFArea2.f90(169)
    REAL(KIND=mp)  DAT                        ! SRC  xFArea2.f90(174)
    REAL(KIND=mp)  AM                         ! SRC  xFArea2.f90(176)
    REAL(KIND=mp)  XK                         ! SRC  xFArea2.f90(178)
    nps=1
    if(ns>1)nps=nscpte(ns-1)+1
    npf=nscpte(ns)
    !	print *,ns,nps,npf
    a=0.
    y1=xspts(nps,1)
    z1=xspts(nps,2)
    if(dragtype == 0) then
        xn1=xspts(nps,3)
    else
        CALL CnvrtCD(nps,e,xn1)
    endif

    !	xn1=xspts(nps,3)
    do i=nps+1,npf
        y2=xspts(i,1)
        z2=xspts(i,2)
        if(dragtype == 0) then
            xn2=xspts(i,3)
        else
            CALL CnvrtCD(i, e, xn2)
        endif

        !		xn2=xspts(i,3)
        !		print *,i,y1,z1,xn1,y2,z2,xn2
        if((z1<e.or.z2<e).and.y2>y1)then
            !		print *,i,y1,z1,xn1,y2,z2,xn2
            if(z1==z2)then  ! segment of channel bottom is flat
                twi=y2-y1
                da=(e-z1)*twi
                a=a+da
            elseif(z1<z2)then
                twi=y2-y1
                dz=z2-z1
                if(e<z2)then
                    twi=twi*(e-z1)/dz
                    dz=e-z1
                endif
                dat=twi*dz/2.
                a=a+dat
                am=am+dat*(y1+twi/3.)
                if(e<=z2)then
                    xk=1.
                else !e>z2
                    da=twi*(e-z2)
                    a=a+da
                endif
            else !z2<z1
                twi=y2-y1
                dz=z1-z2
                if(e<z1)then
                    twi=twi*(e-z2)/dz
                    dz=e-z2
                endif
                dat=twi*dz/2.
                a=a+dat
                if(e<=z1)then
                    xk=1.
                else !e>z1
                    da=twi*(e-z1)
                    a=a+da
                endif
            endif
        endif  ! this i contributes to area
        y1=y2
        z1=z2
    enddo
    end SUBROUTINE AREAOnly
    SUBROUTINE FFindElev(ns,a,e)
    use Support
    IMPLICIT NONE                             ! SRC
    INTEGER, PARAMETER :: mp = KIND(1.0D0)    ! SRC
    REAL(KIND=mp)  DZ                         ! SRC  xFArea2.f90(207)
    REAL(KIND=mp)  E                          ! SRC  xFArea2.f90(208)
    INTEGER        NS                         ! SRC  xFArea2.f90(208)
    REAL(KIND=mp)  AS                         ! SRC  xFArea2.f90(209)
    REAL(KIND=mp)  A                          ! SRC  xFArea2.f90(211)
    REAL(KIND=mp)  FO                         ! SRC  xFArea2.f90(217)
    REAL(KIND=mp)  F                          ! SRC  xFArea2.f90(218)
    INTEGER        ITER                       ! SRC  xFArea2.f90(219)
    !	print *,xsloc(ns,5)
    dz=.1
    if(e<xsloc(ns,5)+.1.or.e>xsloc(ns,5)+5.)then
        as=-1.
        e=xsloc(ns,5)+.1
        do while (as<a)
            e=e+dz
            call areaonly(ns,e,as)
        enddo
    endif
    call areaonly(ns,e,as)
    fo=a-as
    f=1.e9
    iter=0
    do while(abs(f)>.001*a.and.iter<100)
        iter=iter+1
        e=e+dz
        call areaonly(ns,e,as)
        f=a-as
        if(f==0.)return
        dz=-f*dz/(f-fo)
        fo=f
        !		print *,iter,f,dz,as
    enddo
    end SUBROUTINE FFindElev

    SUBROUTINE dxc(ns,dx)
    use Support
    IMPLICIT NONE                             ! SRC
    INTEGER, PARAMETER :: mp = KIND(1.0D0)    ! SRC
    REAL(KIND=mp)  DX                         ! SRC  xFArea2.f90(234)
    INTEGER        NS                         ! SRC  xFArea2.f90(234)
    dx=sqrt((plnxsc(ns,1)-plnxsc(ns-1,1))**2+(plnxsc(ns,2)-plnxsc(ns-1,2))**2)
    end SUBROUTINE dxc

    SUBROUTINE CnvrtCD(n, e, xn1)
    use Support
    IMPLICIT NONE                             ! SRC
    INTEGER, PARAMETER :: mp = KIND(1.0D0)    ! SRC
    REAL(KIND=mp)  XN1                        ! SRC  xFArea2.f90(239)
    REAL(KIND=mp)  H                          ! SRC  xFArea2.f90(240)
    REAL(KIND=mp)  E                          ! SRC  xFArea2.f90(240)
    INTEGER        N                          ! SRC  xFArea2.f90(240)
    xn1 = 1.e6
    h = e - xspts(n,2)
    if(h > 0.001) then
        !			xn1 = (xspts(n,3)*xspts(n,3)*gg)/h**.33333
        xn1 = sqrt(xspts(n,3)*h**.33333/gg)
    endif
    return
    END SUBROUTINE CnvrtCD
