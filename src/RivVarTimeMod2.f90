    MODULE RivVarTimeMod2
    !USE RivVarMod2
    use fm_global
    USE GeometryMod
    IMPLICIT NONE

    TYPE point
        REAL(KIND=mp) :: x, y
    END TYPE point

    TYPE Segment
        TYPE(POINT) p1, p2
    END TYPE Segment

    TYPE BoundingBox
        REAL(KIND=mp) :: xmin, ymin, xmax, ymax
    END TYPE BoundingBox

    type rivvartime
            REAL(KIND=mp), ALLOCATABLE, DIMENSION(:) :: stage, discharge, time_stage, time_discharge
            INTEGER :: nstage, ndischarge

            !!Variables for running FaSTMECH in time varying mode
            REAL(KIND=mp), ALLOCATABLE, DIMENSION(:) :: RatingCurves, TimeSeries
            INTEGER, ALLOCATABLE, DIMENSION(:) ::  RatingCurveType, TimeSeriesType
            INTEGER, ALLOCATABLE, DIMENSION(:) :: TimeSeriesPts, RatingCurvePts
            Type(point), ALLOCATABLE, DIMENSION(:) :: TSPts, RCPts
            Type(BoundingBox),  ALLOCATABLE, DIMENSION(:) :: RCBBx, TSBBx
            CHARACTER(LEN=32), ALLOCATABLE, DIMENSION(:) :: SolNames, SolNames1D, SolNames3D, GridNames3D
            REAL(KIND=mp), ALLOCATABLE, DIMENSION(:) :: TimeIncrements, DischIncrements
            INTEGER :: TSPrintCount = 0;
            INTEGER :: NumTimeSeries, NumRatingCurves

            INTEGER, ALLOCATABLE, DIMENSION(:) :: RCStartPts, TSStartPts

            INTEGER :: VarDischType, DischTSNum, VarStageType, StageTSNum, StageRCNum

            REAL(KIND=mp) :: VarDischStartTime, VarDischEndTime
    end type rivvartime

    !!  VarDischType:
    !!  -0  Constant Disch
    !!  -1  Variable Disch
    !!
    !!  DischTSNum:
    !!      The time series used (0-indexed!) for upstream boundary
    !!
    !!  VarStageType:
    !!  -0  Constant
    !!  -1  Time Series
    !!  -2  Rating Curve
    !!
    !!  StageTSNum:
    !!      The time series number (0-indexed!) to use for downstream boundary
    !!
    !!  StageRCNum
    !!      The Rating Curve number (0-indexed!) to use for downstream boundary

    CONTAINS
    !	            SUBROUTINE  initRatingCurves()
    !	INTEGER     FUNCTION    getNumRC()
    !	            SUBROUTINE  inSegment(P, S, ins)
    !	            SUBROUTINE  getRCSegment(index, segNum, s)
    !	            SUBROUTINE  getRCPt(index, ptNum, p)
    !   INTEGER     FUNCTION    getRCNumSegments(index)
    !	INTEGER     FUNCTION    getRCNumPts(index)
    !	            SUBROUTINE  getInterpRatingCurveValue(index, position, value)
    !
    SUBROUTINE alloc_TSNames(object, numSol)
    IMPLICIT NONE
    type(rivvartime), intent(inout) :: object
    INTEGER, INTENT(IN) :: numSol
    INTEGER :: I, STATUS

    ALLOCATE(object%SolNames(numSol+1), STAT = status)
    ALLOCATE(object%SolNames1D(numSol+1), STAT = status)
    ALLOCATE(object%SolNames3D(numSol+1), STAT = status)
    ALLOCATE(object%GridNames3D(numSol+1), STAT = status)
    ALLOCATE(object%TimeIncrements(numSol+1), STAT = status)
    ALLOCATE(object%DischIncrements(numSol+1), STAT = status)

    END SUBROUTINE
    
    SUBROUTINE dealloc_TSNames(object)
    IMPLICIT NONE
    type(rivvartime), intent(inout) :: object
    INTEGER :: I, STATUS

    DEALLOCATE(object%SolNames, STAT = status)
    DEALLOCATE(object%SolNames1D, STAT = status)
    DEALLOCATE(object%SolNames3D, STAT = status)
    DEALLOCATE(object%GridNames3D, STAT = status)
    DEALLOCATE(object%TimeIncrements, STAT = status)
    DEALLOCATE(object%DischIncrements, STAT = status)

    END SUBROUTINE
    
    SUBROUTINE initRatingCurves(object)
    IMPLICIT NONE
    type(rivvartime), intent(inout) :: object
    REAL(KIND=mp) :: minx, maxx, miny, maxy
    INTEGER :: i,j
    Type(point) :: p1
    minx = 1e12
    miny = 1e12
    maxx = -1e12
    maxy = -1e12
    DO i = 1, getnumrc(object)
        minx = 1e12
        miny = 1e12
        maxx = -1e12
        maxy = -1e12
        Do j = 1, getRCNumPts(object, i)
            call getRCPt(object, i, j, p1)
            if(p1%x < minx) minx = p1%x
            if(p1%x > maxx) maxx = p1%x
            if(p1%y < miny) miny = p1%y
            if(p1%y > maxy) maxy = p1%y
        enddo
        object%RCBBx(i) = BoundingBox(minx, miny, maxx, maxy)
    enddo
    END SUBROUTINE

    SUBROUTINE initTimeSeries(object)
    IMPLICIT NONE
    type(rivvartime), intent(inout) :: object
    REAL(KIND=mp) :: minx, maxx, miny, maxy
    INTEGER :: i,j
    Type(point) :: p1
    minx = 1e12
    miny = 1e12
    maxx = -1e12
    maxy = -1e12
    DO i = 1, getNumTS(object)
        minx = 1e12
        miny = 1e12
        maxx = -1e12
        maxy = -1e12
        Do j = 1, getTSNumPts(object, i)
            call getTSPt(object, i, j, p1)
            if(p1%x < minx) minx = p1%x
            if(p1%x > maxx) maxx = p1%x
            if(p1%y < miny) miny = p1%y
            if(p1%y > maxy) maxy = p1%y
        enddo
        object%TSBBx(i) = BoundingBox(minx, miny, maxx, maxy)
    enddo
    END SUBROUTINE

    INTEGER FUNCTION getNumRC(object)
    type(rivvartime), intent(inout) :: object
    getNumRC = object%NumRatingCurves
    END FUNCTION

    INTEGER FUNCTION getNumTS(object)
    type(rivvartime), intent(inout) :: object
    getNumTS = object%NumTimeSeries
    END FUNCTION

    SUBROUTINE inSegment(P, S, ins)
    IMPLICIT NONE
    TYPE(POINT), INTENT(IN) :: P
    TYPE(SEGMENT), INTENT(IN) :: S
    INTEGER, INTENT(OUT) :: ins

    IF (s%p1%x .ne. s%p2%x) then
        if(s%p1%x<=p%x.and.p%x<=s%p2%x) then
            ins = 1
            return
        endif
        if(s%p1%x >= p%x.and.p%x >= s%p2%x) then
            ins = 1
            return
        endif
    ELSE
        if(s%p1%y <= p%y.and.p%x <= s%p2%y) then
            ins = 1
            return
        endif
        if(s%p1%y >= p%y.and.p%y >= s%p2%y) then
            ins = 1
            return
        endif
    ENDIF
    ins = 0
    END SUBROUTINE

    SUBROUTINE getRCSegment(object, index, segNum, s)
    IMPLICIT NONE
    type(rivvartime), intent(inout) :: object
    INTEGER, INTENT(IN) :: index, segNum
    TYPE(Segment), INTENT(OUT) :: s

    s%p1 = object%RCPts(object%RCStartPts(index)+(segNum-1))
    s%p2 = object%RCPts(object%RCStartPts(index)+(segNum))
    END SUBROUTINE

    SUBROUTINE getTSSegment(object, index, segNum, s)
    IMPLICIT NONE
    type(rivvartime), intent(inout) :: object
    INTEGER, INTENT(IN) :: index, segNum
    TYPE(Segment), INTENT(OUT) :: s

    s%p1 = object%TSPts(object%TSStartPts(index)+(segNum-1))
    s%p2 = object%TSPts(object%TSStartPts(index)+(segNum))
    END SUBROUTINE

    SUBROUTINE getRCPt(object, index, ptNum, p)
    IMPLICIT NONE
    type(rivvartime), intent(inout) :: object
    INTEGER, INTENT(IN) :: index, ptNum
    Type(POINT), INTENT(OUT) :: p

    p%x = object%RCPts(object%RCStartPts(index)+(ptNum-1))%x
    p%y = object%RCPts(object%RCStartPts(index)+(ptNum-1))%y
    END SUBROUTINE

    SUBROUTINE getTSPt(object, index, ptNum, p)
    IMPLICIT NONE
    type(rivvartime), intent(inout) :: object
    INTEGER, INTENT(IN) :: index, ptNum
    Type(POINT), INTENT(INOUT) :: p

    p%x = object%TSPts(object%TSStartPts(index)+(ptNum-1))%x
    p%y = object%TSPts(object%TSStartPts(index)+(ptNum-1))%y
    END SUBROUTINE

    INTEGER FUNCTION getRCNumSegments(object, index)
    IMPLICIT NONE
    type(rivvartime), intent(inout) :: object
    INTEGER, INTENT(IN) :: index
    getRCNumSegments = object%RatingCurvePts(index)-1
    END FUNCTION

    INTEGER FUNCTION getTSNumSegments(object, index)
    IMPLICIT NONE
    type(rivvartime), intent(inout) :: object
    INTEGER, INTENT(IN) :: index
    getTSNumSegments = object%TimeSeriesPts(index)-1
    END FUNCTION

    INTEGER FUNCTION getRCNumPts(object, index)
    IMPLICIT NONE
    type(rivvartime), intent(inout) :: object
    INTEGER, INTENT(IN) :: index
    getRCNumPTs =  object%RatingCurvePts(index)
    END FUNCTION

    INTEGER FUNCTION getTSNumPts(object, index)
    IMPLICIT NONE
    type(rivvartime), intent(inout) :: object
    INTEGER, INTENT(IN) :: index
    getTSNumPTs =  object%TimeSeriesPts(index)
    END FUNCTION

    SUBROUTINE getInterpRatingCurveValue(object, index, position, value)
    IMPLICIT NONE
    type(rivvartime), intent(inout) :: object
    INTEGER, INTENT(IN) :: index
    REAL(KIND=mp), INTENT(IN) :: position
    REAL(KIND=mp), INTENT(OUT) :: value
    INTEGER :: numSeg, i
    integer, parameter :: dim_num = 2
    REAL(KIND=mp), DIMENSION(2) :: p1, p2, q1, q2, v
    TYPE(Segment) :: seg
    TYPE(Point) :: pnt
    INTEGER :: retval
    seg = segment(point(0,0), point(0,0))
    p1(1) = position
    p1(2) = object%RCBBx(index)%ymin
    pnt%x = position
    pnt%y = object%RCBBx(index)%ymin

    p2(1) = position
    p2(2) = object%RCBBx(index)%ymax

    DO i = 1,getRCNumSegments(object, index)
        CALL getRCSegment(object, index, i, seg)
        CALL inSegment(pnt, seg, retval)
        if(retval == 1) then
            q1(1) = seg%p1%x
            q1(2) = seg%p1%y
            q2(1) = seg%p2%x
            q2(2) = seg%p2%y
            call lines_exp_int_2d ( p1, p2, q1, q2, retval, v )
            if(retval == 1) then
                value = v(2)
                exit
            else if(retval == 2) then
                value = (seg%p2%y - seg%p1%y)/2.
                exit
            else
                value = -1.0
                exit
            endif
        endif
    ENDDO
    END SUBROUTINE

    SUBROUTINE getInterpTimeSeriesValue(object, index, position, value)
    IMPLICIT NONE
    type(rivvartime), intent(inout) :: object
    INTEGER, INTENT(IN) :: index
    REAL(KIND=mp), INTENT(IN) :: position
    REAL(KIND=mp), INTENT(OUT) :: value
    INTEGER :: numSeg, i
    REAL(KIND=mp), DIMENSION(2) :: p1, p2, q1, q2, v
    TYPE(Segment) :: seg
    TYPE(Point) :: pnt
    INTEGER :: retval
    seg = segment(point(0,0), point(0,0))
    p1(1) = position
    p1(2) = object%TSBBx(index)%ymin
    pnt%x = position
    pnt%y = object%TSBBx(index)%ymin

    p2(1) = position
    p2(2) = object%TSBBx(index)%ymax

    DO i = 1,getTSNumSegments(object, index)
        CALL getTSSegment(object, index, i, seg)
        CALL inSegment(pnt, seg, retval)
        if(retval == 1) then
            q1(1) = seg%p1%x
            q1(2) = seg%p1%y
            q2(1) = seg%p2%x
            q2(2) = seg%p2%y
            call lines_exp_int_2d ( p1, p2, q1, q2, retval, v )
            if(retval == 1) then
                value = v(2)
                exit
            else if(retval == 2) then
                value = (seg%p2%y - seg%p1%y)/2.
                exit
            else
                value = -1.0
                exit
            endif
        endif
    ENDDO
    END SUBROUTINE

    SUBROUTINE setRatingCurveNumPts(object, dim, val)
    IMPLICIT NONE
    type(rivvartime), intent(inout) :: object
    INTEGER, INTENT(IN) :: dim
    INTEGER, INTENT(IN) :: val(:)
    INTEGER :: I, STATUS

    ALLOCATE(object%RatingCurvePts(dim), STAT = status)
    ALLOCATE(object%RCStartPts(dim), STAT = status)
    ALLOCATE(object%RCBBx(dim), STAT = status)

    DO I = 1,dim
        object%RatingCurvePts(I)= val(i)
        if(i == 1) then
            object%RCStartPts(I) = 1
        ELSE
            object%RCStartPts(I) = object%RCStartPts(i-1)+(val(i-1))
        ENDIF
    ENDDO

    object%NumRatingCurves = dim
    END SUBROUTINE

    SUBROUTINE setTimeSeriesNumPts(object, dim, val)
    IMPLICIT NONE
    type(rivvartime), intent(inout) :: object
    INTEGER, INTENT(IN) :: dim
    INTEGER, INTENT(IN) :: val(:)
    INTEGER :: I, STATUS

    ALLOCATE(object%TimeSeriesPts(dim), STAT = status)
    ALLOCATE(object%TSStartPts(dim), STAT = status)
    ALLOCATE(object%TSBBx(dim), STAT = status)

    DO I = 1,dim
        object%TimeSeriesPts(I)= val(i)
        if(i==1) then
            object%TSStartPts(i) = 1
        else
            object%TSStartPts(i) = object%TSStartPts(i-1)+(val(i-1))
        endif
    ENDDO

    object%NumTimeSeries = dim
    END SUBROUTINE

    SUBROUTINE setTimeSeries(object, dim, val)
    IMPLICIT NONE
    type(rivvartime), intent(inout) :: object
    INTEGER, INTENT(IN) :: dim
    REAL(KIND=mp), INTENT(IN) :: val(:)
    INTEGER :: I, STATUS, count

    ALLOCATE(object%TimeSeries(dim), STAT = status)
    ALLOCATE(object%TSPts((dim/2)), STAT = status)

    DO I = 1,dim/2
        count = (i-1)*2
        object%TSPts(i)%x = val(count+1)
        object%TSPts(i)%y = val(count+2)
    ENDDO

    DO I = 1,dim
        object%TimeSeries(I)= val(i)
    ENDDO

    END SUBROUTINE

    SUBROUTINE setRatingCurves(object, dim, val)
    IMPLICIT NONE
    type(rivvartime), intent(inout) :: object
    INTEGER, INTENT(IN) :: dim
    REAL(KIND=mp), INTENT(IN) :: val(:)
    INTEGER :: I, STATUS, count

    ALLOCATE(object%RatingCurves(dim), STAT = status)
    ALLOCATE(object%RCPts((dim/2)), STAT = status)

    DO I = 1,dim/2
        count = (i-1)*2
        object%RCPts(i)%x = val(count+1)
        object%RCPts(i)%y = val(count+2)
    ENDDO

    DO I = 1,dim
        object%RatingCurves(I)= val(i)
    ENDDO

    END SUBROUTINE

    SUBROUTINE setTimeSeriesType(object, dim, val)
    IMPLICIT NONE
    type(rivvartime), intent(inout) :: object
    INTEGER, INTENT(IN) :: dim
    INTEGER, INTENT(IN) :: val(:)
    INTEGER :: I, STATUS

    ALLOCATE(object%TimeSeriesType(dim), STAT = status)

    DO I = 1,dim
        object%TimeSeriesType(I)= val(i)
    ENDDO

    END SUBROUTINE

    SUBROUTINE setRatingCurveType(object, dim, val)
    IMPLICIT NONE
    type(rivvartime), intent(inout) :: object
    INTEGER, INTENT(IN) :: dim
    INTEGER, INTENT(IN) :: val(:)
    INTEGER :: I, STATUS

    ALLOCATE(object%RatingCurveType(dim), STAT = status)

    DO I = 1,dim
        object%RatingCurveType(I)= val(i)
    ENDDO

    END SUBROUTINE

    SUBROUTINE dealloc_TimeSeries(object)
    IMPLICIT NONE
    type(rivvartime), intent(inout) :: object
    INTEGER :: status
    DEALLOCATE(object%TimeSeriesType, STAT = status)
    DEALLOCATE(object%TimeSeries, STAT = status)
    DEALLOCATE(object%TSPts, STAT = status)
    DEALLOCATE(object%TimeSeriesPts, STAT = status)
    DEALLOCATE(object%TSStartPts, STAT = status)
    DEALLOCATE(object%TSBBx, STAT = status)

    END SUBROUTINE

    SUBROUTINE dealloc_RatingCurves(object)
    IMPLICIT NONE
    type(rivvartime), intent(inout) :: object
    INTEGER :: status
    DEALLOCATE(object%RatingCurveType, STAT = status)
    DEALLOCATE(object%RatingCurves, STAT = status)
    DEALLOCATE(object%RCPts, STAT = status)
    DEALLOCATE(object%RatingCurvePts, STAT = status)
    DEALLOCATE(object%RCStartPts, STAT = status)
    DEALLOCATE(object%RCBBx, STAT = status)

    END SUBROUTINE

    END MODULE
