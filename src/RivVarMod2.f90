    module rivvarmod2
    use fm_global
    implicit none
    !	save
    !	integer, dimension(ns, nn) :: ibc
    !allocate real common types
    ! machine precision:
    type rivvar
        real(kind = mp) :: g, rho, vkc
        integer :: nclusters, nwetnodes, noldwetnodes

        !real(kind=mp) :: elevoffset=1e12
        real(kind=mp) :: ds,dn,mo!,hwt
        real(kind=mp) :: fcos, fsin
        !real(kind=mp) :: lbwse, rbwse, vdcds, vacds
        real(kind=mp) :: tottime, vardt, ptime
        real(kind=mp) :: wmax
        !allocate integer common types
        integer :: iter
        integer :: nm
        integer :: nsteps, nct
        integer :: numvbcpts
        integer :: errorcode
        logical :: flumebndry
        real(kind=mp) :: ddisch, dstage, newstage, oldstage, newdisch, olddisch, dsstage

        double precision:: xshift, yshift

        integer, allocatable, dimension(:,:) :: ibc, iibc, ribc, icon, tibc
        real(kind=mp), allocatable, dimension(:) :: r, w, xo, yo, hav
        real(kind=mp), dimension(:,:), pointer :: taus, taun, hl, eta, rn, cd, cdv, totcd
        real(kind=mp), allocatable, dimension(:,:) :: mineta
        real(kind=mp), allocatable, dimension(:,:) :: znaught
        real(kind=mp), allocatable, dimension(:,:) :: u, v, e
        real(kind=mp), allocatable, dimension(:,:) :: iu, iv, iwse, ie, ihl
        real(kind=mp), allocatable, dimension(:,:) :: x, y
        real(kind=mp), allocatable, dimension(:,:) :: harea

        real(kind=mp), allocatable, dimension(:,:,:) :: uz, vz, zz
        real(kind=mp), allocatable, dimension(:,:) :: con, qs, qn

        real(kind=mp), allocatable, dimension(:,:) :: fracs, hfine

        !allocate seperate arrays to deal with extend lower boundary
        integer, allocatable, dimension(:,:) :: ibc2
        real(kind=mp), allocatable, dimension(:) :: r2, w2, hav2
        real(kind=mp), allocatable, dimension(:) :: xo2, yo2, phi2
        real(kind=mp), allocatable, dimension(:) :: phirotation
        real(kind=mp), allocatable, dimension(:,:) :: u2, v2, e2
        real(kind=mp), allocatable, dimension(:,:) :: eta2, cd2, cdv2, znaught2, mineta2

        !integer, allocatable, dimension(:,:) :: iedger, iedgel
        !integer, allocatable, dimension(:,:) :: icountl, icountr
        !integer, allocatable, dimension(:) :: iswitch

        real(kind=mp) mo2

        integer :: ns2
        integer :: ns, nn, nz
        INTEGER, PUBLIC :: nsext
        integer :: cgnsfileid

        REAL(KIND=mp) :: stot, scals

    end type

    contains
    subroutine alloc_common3d(this)
    type(rivvar), intent(inout) :: this
    integer :: status, tmpns, tmpnn, tmpnz
    !			tmpns = ns2+nsext
    !			tmp = ns
    !allocate real 3d common types
    tmpns = this%ns
    tmpnn = this%nn
    tmpnz = this%nz
    allocate(this%uz(tmpns, tmpnn, tmpnz), stat = status)
    allocate(this%vz(tmpns, tmpnn, tmpnz), stat = status)
    allocate(this%zz(tmpns, tmpnn, tmpnz), stat = status)
    this%uz = 0.
    this%vz = 0.
    this%zz = 0.
    end subroutine

    subroutine alloc_init2d(object, ns2, nn)
    type(rivvar), intent(inout) :: object
    integer, intent(in) :: ns2, nn
    integer :: status
    allocate(object%iv(ns2, nn), stat = status)
    allocate(object%iu(ns2, nn), stat = status)
    allocate(object%ie(ns2, nn), stat = status)
    allocate(object%iwse(ns2, nn), stat = status)
    allocate(object%ihl(ns2, nn), stat = status)
    allocate(object%iibc(ns2, nn), stat = status)
    end subroutine

    subroutine dealloc_init2d(object)
    type(rivvar), intent(inout) :: object
    integer :: status
    deallocate(object%iv, stat = status)
    deallocate(object%iu, stat = status)
    deallocate(object%iwse, stat = status)
    deallocate(object%ie, stat = status)
    deallocate(object%iibc, stat = status)
    deallocate(object%ihl, stat = status)
    end subroutine

    subroutine alloc_common2d(object, ns2, nsext, nn)
    type(rivvar), intent(inout) :: object
    integer, intent(in) :: ns2, nsext, nn
    integer :: status, ns
    object%ns = ns2 + nsext
    ns = ns2+nsext !added to avoid using object% for each instance below
    !allocate integer common types
    allocate(object%ibc(ns, nn), stat = status)
    allocate(object%tibc(ns, nn), stat = status)
    allocate(object%ribc(ns, nn), stat = status)
    allocate(object%icon(ns, nn), stat = status)
    !allocate real 1d common types
    allocate(object%r(ns), stat = status)
    allocate(object%w(ns), stat = status)
    allocate(object%xo(ns), stat = status)
    allocate(object%yo(ns), stat = status)
    allocate(object%hav(ns), stat = status)
    allocate(object%phirotation(ns), stat = status)
    !allocate real 2d common types
    allocate(object%u(ns, nn), stat = status)
    allocate(object%v(ns, nn), stat = status)
    allocate(object%e(ns, nn), stat = status)
    allocate(object%taus(ns, nn), stat = status)
    allocate(object%taun(ns, nn), stat = status)
    allocate(object%hl(ns, nn), stat = status)
    allocate(object%eta(ns, nn), stat = status)
    !allocate(mineta(ns, nn), stat = status)
    allocate(object%rn(ns, nn), stat = status)
    allocate(object%cd(ns, nn), stat = status)
    allocate(object%cdv(ns, nn), stat = status)
    allocate(object%totcd(ns, nn), stat = status)
    allocate(object%znaught(ns, nn), stat = status)
    allocate(object%con(ns, nn), stat = status)
    allocate(object%qs(ns, nn), stat = status)
    allocate(object%qn(ns, nn), stat = status)
    allocate(object%x(ns, nn), stat = status)
    allocate(object%y(ns, nn), stat = status)
    allocate(object%harea(ns, nn), stat = status)
    allocate(object%fracs(ns, nn), stat = status)
    allocate(object%hfine(ns, nn), stat = status)

    !!allocate for craig dixon changes
    !allocate(object%iedger(ns, 15), stat = status)
    !allocate(object%iedgel(ns, 15), stat = status)
    !allocate(object%icountr(ns, nn), stat = status)
    !allocate(object%icountl(ns, nn), stat = status)
    !allocate(object%iswitch(0:nn), stat = status)

    !allocate variables for extended bounds
    allocate(object%r2(ns2+nsext), stat = status)
    allocate(object%w2(ns2), stat = status)
    allocate(object%xo2(ns2), stat = status)
    allocate(object%yo2(ns2), stat = status)
    allocate(object%hav2(ns2), stat = status)
    allocate(object%phi2(ns2), stat = status)
    allocate(object%cd2(ns2, nn), stat = status)
    allocate(object%cdv2(ns2, nn), stat = status)
    allocate(object%znaught2(ns2, nn), stat = status)

    allocate(object%eta2(ns2, nn), stat = status)
    !allocate(mineta2(ns2, nn), stat = status)
    allocate(object%ibc2(ns2, nn), stat = status)
    allocate(object%u2(ns2, nn), stat = status)
    allocate(object%v2(ns2, nn), stat = status)
    allocate(object%e2(ns2, nn), stat = status)
    end subroutine

    subroutine dealloc_common3d(object)
    type(rivvar), intent(inout) :: object
    integer :: status
    !allocate real 3d common types
    deallocate(object%uz, stat = status)
    deallocate(object%vz, stat = status)
    deallocate(object%zz, stat = status)
    end subroutine

    subroutine dealloc_common2d(object)
    type(rivvar), intent(inout) :: object
    integer :: status
    !allocate integer common types
    deallocate(object%ibc, stat = status)
    deallocate(object%tibc, stat = status)
    deallocate(object%ribc, stat = status)
    deallocate(object%icon, stat = status)
    !allocate real 1d common types
    deallocate(object%r, stat = status)
    deallocate(object%w, stat = status)
    deallocate(object%xo, stat = status)
    deallocate(object%yo, stat = status)
    deallocate(object%hav, stat = status)
    !allocate real 2d common types
    deallocate(object%u, stat = status)
    deallocate(object%v, stat = status)
    deallocate(object%e, stat = status)
    deallocate(object%taus, stat = status)
    deallocate(object%taun, stat = status)
    deallocate(object%hl, stat = status)
    deallocate(object%eta, stat = status)
    !deallocate(mineta, stat = status)
    deallocate(object%rn, stat = status)
    deallocate(object%cd, stat = status)
    deallocate(object%cdv, stat = status)
    deallocate(object%totcd, stat = status)
    deallocate(object%znaught, stat = status)
    deallocate(object%con, stat = status)
    deallocate(object%qs, stat = status)
    deallocate(object%qn, stat = status)
    deallocate(object%x, stat = status)
    deallocate(object%y, stat = status)
    deallocate(object%harea, stat = status)
    deallocate(object%fracs, stat = status)
    deallocate(object%hfine, stat = status)
    !!allocate for craig dixon changes
    !deallocate(iedger, stat = status)
    !deallocate(iedgel, stat = status)
    !deallocate(icountr, stat = status)
    !deallocate(icountl, stat = status)
    !deallocate(iswitch, stat = status)

    !allocate variables for extended bounds
    deallocate(object%r2, stat = status)
    deallocate(object%w2, stat = status)
    deallocate(object%xo2, stat = status)
    deallocate(object%yo2, stat = status)
    deallocate(object%hav2, stat = status)
    deallocate(object%phi2, stat = status)
    deallocate(object%phirotation, stat = status)
    deallocate(object%cd2, stat = status)
    deallocate(object%cdv2, stat = status)
    deallocate(object%znaught2, stat = status)

    deallocate(object%eta2, stat = status)
    !deallocate(mineta2, stat = status)
    deallocate(object%ibc2, stat = status)
    deallocate(object%u2, stat = status)
    deallocate(object%v2, stat = status)
    deallocate(object%e2, stat = status)
    end subroutine

    subroutine calc_area(tns, tnn, phirotation, x, y, xo, yo, nm, dn, area)
    implicit none
    integer, intent (in) :: tns, tnn
    real(kind=mp), intent(in), dimension(:) :: phirotation
    real(kind=mp), intent (out), dimension (:,:) :: area
    real(kind=mp), intent(inout), dimension(:,:):: x,y
    real(kind=mp), intent (in), dimension (:) :: xo, yo
    integer, intent(in) :: nm
    real(kind=mp), intent(in) :: dn

    integer :: i, j
    real(kind=mp) :: x1, x2, x3, x4
    real(kind=mp) :: y1, y2, y3, y4
    real(kind=mp) :: area1, area2

    real(kind=mp) :: rsin, rcos, ux, uy, uxnew, uynew

    integer :: icount
    do i=1,tns
        rcos=COS(phirotation(i))
        rsin=sin(phirotation(i))

        do j=1,tnn
            x(i,j)=xo(i)+(nm-j)*dn*sin(phirotation(i))
            y(i,j)=yo(i)+(j-nm)*dn*cos(phirotation(i))
        enddo
    enddo
    ! loop through the mesh and calculate areas for each point except the edges of the mesh

    do i=2,tns-1
        do j=2,tnn-1
            x1=(x(i-1,j-1)+x(i-1,j)+x(i,j)+x(i,j-1))/4.
            x2=(x(i-1,j)+x(i-1,j+1)+x(i,j+1)+x(i,j))/4.
            x3=(x(i,j)+x(i,j+1)+x(i+1,j+1)+x(i+1,j))/4.
            x4=(x(i,j-1)+x(i,j)+x(i+1,j)+x(i+1,j-1))/4.
            y1=(y(i-1,j-1)+y(i-1,j)+y(i,j)+y(i,j-1))/4.
            y2=(y(i-1,j)+y(i-1,j+1)+y(i,j+1)+y(i,j))/4.
            y3=(y(i,j)+y(i,j+1)+y(i+1,j+1)+y(i+1,j))/4.
            y4=(y(i,j-1)+y(i,j)+y(i+1,j)+y(i+1,j-1))/4.
            !c          using area of a triangle algorithm in standard mathmatical tables book page 362
            !c          0.5*(x1*y2+x2*y3+x3*y1-y1*x2-y2*x3-y3*x1)
            !c          where in the first case x1=>x1,x2=>x4,x3=>x2
            !c          where in the second case x1=>x3,x2=>x2,x3=>x4
            area1=0.5*(x1*y4+x4*y2+x2*y1-y1*x4-y4*x2-y2*x1)
            area2=0.5*(x3*y2+x2*y4+x4*y3-y3*x2-y2*x4-y4*x3)
            area(i,j)=area1+area2
            !c           write(9,*) i,j,x1,x2,x3,x4,y1,y2,y3,y4,area1,area2  ! debug
        end do
    end do

    !c loop through the mesh and estimate the cell areas for the edge of the mesh
    !c and write the node number,x coord, y coord, and area out to a file

    icount=0
    ! c       write(9,*)'node,x,y,area'

    do i=1,tns
        do j=1,tnn
            icount=icount+1
            !c       deal with left edge of the mesh (take area from interior)
            if(j.eq.1)then
                if(i.eq.1)then
                    area(i,j)=area(2,2)
                elseif(i.eq.tns)then
                    area(i,j)=area(tns-1,2)
                else
                    area(i,j)=area(i,2)
                endif
                !c       deal with right edge of the mesh (take area from interior)
            elseif(j.eq.tnn)then
                if(i.eq.1)then
                    area(i,j)=area(2,tnn-1)
                elseif(i.eq.tns)then
                    area(i,j)=area(tns-1,tnn-1)
                else
                    area(i,j)=area(i,tnn-1)
                endif
            endif
            !c       deal with the interior ends of the mesh (top and bottom of the mesh) areas
            if(i.eq.1)then
                if(j.gt.1.and.j.lt.tnn)then
                    area(i,j)=area(2,j)
                endif
            elseif(i.eq.tns)then
                if(j.gt.1.and.j.lt.tnn)then
                    area(i,j)=area(tns-1,j)
                endif
            endif

            !!c       write out the x coordinates ,y coordinates, and area values
            !        x(i,j)=x(i,j)+offsetx
            !        y(i,j)=y(i,j)+offsety
            !        write(9,210) icount,x(i,j),y(i,j),area(i,j)

        end do
    end do
    end subroutine calc_area

    subroutine updateibc(object)
    implicit none
    type(rivvar), intent(inout) :: object
    integer :: i,j
    object%tibc = object%ibc
    do i = 1,object%ns
        do j = 1,object%nn
            if(object%ibc(i,j).gt.0) then
                object%tibc(i,j) = -1
                !	        else if(j == 1) then
                !!	            ibc(i,j) = 1
                !	        else if(j == nn) then
                !!	            ibc(i,j) = 2
            endif

        enddo
    enddo
    do i = 2,object%ns-1
        do j = 2, object%nn-1
            if(object%ibc(i,j).ne.0.and.object%ibc(i-1,j).eq.0) then
                object%tibc(i,j)=4
                !	        elseif(ibc(i,j).eq.-1.and.ibc(i,j-1).eq.0)then
                !	            tibc(i,j)=1
                !	        elseif(ibc(i,j).eq.-1.and.ibc(i,j+1).eq.0)then
                !	            tibc(i,j)=2
            endif
            !	        if(ibc(i,j).eq.-1.and.ibc(i+1,j).eq.0) then
            !	            tibc(i+1,j)=6
            !	        endif
        enddo
    enddo
    !enforce lateral boundary nodes
    if(object%flumebndry) then
        do i = 1, object%ns
            object%tibc(i,1) = 1;
            object%tibc(i,object%nn) = 2;
        enddo
    else
        do i = 1, object%ns
            object%tibc(i,1) = 0
            object%u(i,1) = 0;
            object%v(i,1) = 0;

            object%tibc(i,object%nn) = 0
            object%u(i,object%nn) = 0;
            object%v(i,object%nn) = 0;
        enddo
    endif

    object%ibc = object%tibc

    end subroutine updateibc
    
    SUBROUTINE CalcGridMetrics(object)
    IMPLICIT NONE
    type(rivvar), intent(inout) :: object
    REAL(kind=mp) :: tTheta, tdy, tdx, xext, yext, tds, tdr, dtphi
    REAL(kind=mp) :: ds,slin, dx, dy, dphi, tdphi, tmpwidth
    DOUBLE PRECISION :: xc1, yc1

    INTEGER :: i, tnn, tns

    INTEGER :: jmid
    tnn = object%nn
    tns = object%ns
    jmid = (object%nn+1)/2
    object%r2 = 0
    ! Find length along centerline
    object%stot = 0.
    DO i = 2,object%ns2
        ds = sqrt((object%x(i,jmid)-object%x(i-1,jmid))**2 + (object%y(i,jmid)-object%y(i-1,jmid))**2)
        object%stot = object%stot+ds
    ENDDO


    object%scals = object%stot/(object%ns2-1)
    object%xshift = object%x(1,jmid)
    object%yshift = object%y(1,jmid)

    tmpwidth = sqrt((object%x(1,1)-object%x(1,tnn))**2 + (object%y(1,1)-object%y(1,tnn))**2)
    DO i = 1,object%ns2
        object%xo(i) = object%x(i,jmid)-object%xshift
        object%yo(i) = object%y(i,jmid)-object%yshift
        object%w2(i) = tmpwidth
    ENDDO

    slin = sqrt(object%xo(object%ns2)**2+object%yo(object%ns2)**2)
    object%fcos = object%xo(object%ns2)/slin
    object%fsin = object%yo(object%ns2)/slin

    DO i = 1,object%ns2
        xc1 = object%xo(i)*object%fcos+object%yo(i)*object%fsin
        yc1 = object%yo(i)*object%fcos-object%xo(i)*object%fsin
        object%xo(i) = xc1
        object%yo(i) = yc1
    ENDDO

    DO i=2,object%ns2
        dx = object%xo(i)-object%xo(i-1)
        dy = object%yo(i)-object%yo(i-1)
        IF(dx.eq.0) THEN
            IF(DY.gt.0) THEN
                object%phirotation(i) = acos(-1.)/2.
            ELSE
                object%phirotation(i) = -1. * acos(-1.)/2.
            ENDIF
        ELSE
            object%phirotation(i) = atan2(dy,dx)
        ENDIF
        object%phi2(i) = object%phirotation(i)
    ENDDO

    object%phirotation(1) = (2.*object%phirotation(2))-object%phirotation(3)
    object%phi2(1) = object%phirotation(1)

    DO i = 2,object%ns2
        dphi = object%phirotation(i)-object%phirotation(i-1)
        IF(dphi == 0.) THEN
            if(object%r2(i-1) < 0) THEN
                object%r2(i) = -100000000.
            ELSE
                object%r2(i) = 100000000.
            ENDIF
        ELSE
            object%r2(i) = object%scals/dphi
        ENDIF
    ENDDO

    object%r2(1) = (2.*object%r2(2))-object%r2(3)

    tdphi = object%phirotation(object%ns2)-object%phirotation(object%ns2-1)
    tdr = (tdphi)/object%nsext
    DO i=object%ns2+1,object%ns2+object%nsext
        dphi = ((object%ns2+object%nsext+1)-i)*tdr
        object%phirotation(i) = object%phirotation(i-1)+dphi

        IF(dphi == 0) THEN
            IF(object%r2(i-1) < 0) THEN
                object%r2(i) = -100000000.
            ELSE
                object%r2(i) = 100000000.
            ENDIF
        ELSE
            object%r2(i) = object%scals/dphi
        ENDIF
        dx = object%scals*cos(object%phirotation(i))
        dy = object%scals*sin(object%phirotation(i))
        object%xo(i) = object%xo(i-1) + dx
        object%yo(i) = object%yo(i-1) + dy
        dx = object%xo(i)-object%xo(i-1)
        dy = object%yo(i)-object%yo(i-1)
        IF(dx.eq.0) THEN
            IF(dy.gt.0) THEN
                object%phirotation(i)=acos(-1.)/2.
            ELSEIF(dy.le.0) THEN
                object%phirotation(i)=-1.*acos(-1.)/2.
            ENDIF
        ELSE
            object%phirotation(i)=atan2(dy,dx)
        ENDIF

    ENDDO

    END SUBROUTINE CalcGridMetrics
    !subroutine resetgridextension()
    !implicit none
    !
    !integer :: i,j
    !real(kind=mp) :: mag
    !do i=ns2+1,ns2+nsext
    !    do j = 1,nn
    !        eta(i,j)=eta2(ns2,j) - ((i-ns2)*ds*nsextslope)
    !        !mineta(i,j)=mineta2(ns2,j) - ((i-ns2)*ds*nsextslope)
    !        if(ibc(i,j).ne.0) then
    !            mag = sqrt(u(ns2,j)**2.+v(ns2,j)**2.)
    !            v(i,j) = v(i,j) - (i*(-v(ns2,j)/(ns2+nsext)))
    !            u(i,j) = sqrt(mag**2. - v(i,j)**2.)
    !        endif
    !    enddo
    !enddo
    !end subroutine resetgridextension
    end module rivvarmod2
