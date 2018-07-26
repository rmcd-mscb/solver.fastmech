MODULE RivVarTimeMod
	USE RivVarMod
	USE GeometryMod
	IMPLICIT NONE
	
	TYPE point
        REAL :: x, y
    END TYPE point
    
    TYPE Segment
        TYPE(POINT) p1, p2
    END TYPE Segment
    
    TYPE BoundingBox
        REAL :: xmin, ymin, xmax, ymax
    END TYPE BoundingBox
	

	REAL, ALLOCATABLE, DIMENSION(:) :: stage, discharge, time_stage, time_discharge
	INTEGER :: nstage, ndischarge
	
!!Variables for running FaSTMECH in time varying mode
    REAL, ALLOCATABLE, DIMENSION(:) :: RatingCurves, TimeSeries
    INTEGER, ALLOCATABLE, DIMENSION(:) ::  RatingCurveType, TimeSeriesType
    INTEGER, ALLOCATABLE, DIMENSION(:) :: TimeSeriesPts, RatingCurvePts
    Type(point), ALLOCATABLE, DIMENSION(:) :: TSPts, RCPts
    Type(BoundingBox),  ALLOCATABLE, DIMENSION(:) :: RCBBx, TSBBx
    CHARACTER*32, ALLOCATABLE, DIMENSION(:) :: SolNames, SolNames1D, SolNames3D, GridNames3D
    REAL, ALLOCATABLE, DIMENSION(:) :: TimeIncrements, DischIncrements
    INTEGER :: TSPrintCount = 0;
    INTEGER :: NumTimeSeries, NumRatingCurves
    
    INTEGER, ALLOCATABLE, DIMENSION(:) :: RCStartPts, TSStartPts

    INTEGER :: VarDischType, DischTSNum, VarStageType, StageTSNum, StageRCNum
    
    REAL :: VarDischStartTime, VarDischEndTime
    
    
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
	SUBROUTINE alloc_TSNames(numSol)
	     IMPLICIT NONE
	    INTEGER, INTENT(IN) :: numSol
	    INTEGER :: I, STATUS
	    
		ALLOCATE(SolNames(numSol+1), STAT = status)
		ALLOCATE(SolNames1D(numSol+1), STAT = status)
		ALLOCATE(SolNames3D(numSol+1), STAT = status)
		ALLOCATE(GridNames3D(numSol+1), STAT = status)
		ALLOCATE(TimeIncrements(numSol+1), STAT = status)
		ALLOCATE(DischIncrements(numSol+1), STAT = status)

    END SUBROUTINE
	SUBROUTINE dealloc_TSNames()
	    IMPLICIT NONE
	    INTEGER :: I, STATUS
	    
		DEALLOCATE(SolNames, STAT = status)
		DEALLOCATE(SolNames1D, STAT = status)
		DEALLOCATE(SolNames3D, STAT = status)
		DEALLOCATE(GridNames3D, STAT = status)
		DEALLOCATE(TimeIncrements, STAT = status)
		DEALLOCATE(DischIncrements, STAT = status)

    END SUBROUTINE
    SUBROUTINE initRatingCurves()
	    IMPLICIT NONE
	    REAL :: minx, maxx, miny, maxy
	    INTEGER :: i,j
	    Type(point) :: p1
	    minx = 1e12
	    miny = 1e12
	    maxx = -1e12
	    maxy = -1e12
	    DO i = 1, getnumrc()
	        minx = 1e12
	        miny = 1e12
	        maxx = -1e12
	        maxy = -1e12
	        Do j = 1, getRCNumPts(i)
	            call getRCPt(i, j, p1)
	                if(p1%x < minx) minx = p1%x
	                if(p1%x > maxx) maxx = p1%x
	                if(p1%y < miny) miny = p1%y
	                if(p1%y > maxy) maxy = p1%y
	        enddo
	        RCBBx(i) = BoundingBox(minx, miny, maxx, maxy)
	    enddo
    END SUBROUTINE

	SUBROUTINE initTimeSeries()
	    IMPLICIT NONE
	    REAL :: minx, maxx, miny, maxy
	    INTEGER :: i,j
	    Type(point) :: p1
	    minx = 1e12
	    miny = 1e12
	    maxx = -1e12
	    maxy = -1e12
	    DO i = 1, getNumTS()
	        minx = 1e12
	        miny = 1e12
	        maxx = -1e12
	        maxy = -1e12
	        Do j = 1, getTSNumPts(i)
	            call getTSPt(i, j, p1)
	                if(p1%x < minx) minx = p1%x
	                if(p1%x > maxx) maxx = p1%x
	                if(p1%y < miny) miny = p1%y
	                if(p1%y > maxy) maxy = p1%y
	        enddo
	        TSBBx(i) = BoundingBox(minx, miny, maxx, maxy)
	    enddo
    END SUBROUTINE
    
	INTEGER FUNCTION getNumRC()
	    getNumRC = NumRatingCurves
	END FUNCTION
	
	INTEGER FUNCTION getNumTS()
	    getNumTS = NumTimeSeries
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
	
	SUBROUTINE getRCSegment(index, segNum, s)
	    IMPLICIT NONE
	    INTEGER, INTENT(IN) :: index, segNum
	    TYPE(Segment), INTENT(OUT) :: s
	    
	    s%p1 = RCPts(RCStartPts(index)+(segNum-1))
        s%p2 = RCPts(RCStartPts(index)+(segNum))
    END SUBROUTINE
	
	SUBROUTINE getTSSegment(index, segNum, s)
	    IMPLICIT NONE
	    INTEGER, INTENT(IN) :: index, segNum
	    TYPE(Segment), INTENT(OUT) :: s
	    
	    s%p1 = TSPts(TSStartPts(index)+(segNum-1))
        s%p2 = TSPts(TSStartPts(index)+(segNum))
    END SUBROUTINE
    
   	SUBROUTINE getRCPt(index, ptNum, p)
	    IMPLICIT NONE
	    INTEGER, INTENT(IN) :: index, ptNum
	    Type(POINT), INTENT(OUT) :: p
	    
	    p%x = RCPts(RCStartPts(index)+(ptNum-1))%x
	    p%y = RCPts(RCStartPts(index)+(ptNum-1))%y
    END SUBROUTINE

    SUBROUTINE getTSPt(index, ptNum, p)
	    IMPLICIT NONE
	    INTEGER, INTENT(IN) :: index, ptNum
	    Type(POINT), INTENT(INOUT) :: p
	    
	    p%x = TSPts(TSStartPts(index)+(ptNum-1))%x
	    p%y = TSPts(TSStartPts(index)+(ptNum-1))%y
    END SUBROUTINE
    
    INTEGER FUNCTION getRCNumSegments(index)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: index
        getRCNumSegments = RatingCurvePts(index)-1
    END FUNCTION
	    
    INTEGER FUNCTION getTSNumSegments(index)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: index
        getTSNumSegments = TimeSeriesPts(index)-1
    END FUNCTION
    
   	INTEGER FUNCTION getRCNumPts(index)
	    IMPLICIT NONE
	    INTEGER, INTENT(IN) :: index
	    getRCNumPTs =  RatingCurvePts(index)
	END FUNCTION
	
	INTEGER FUNCTION getTSNumPts(index)
	    IMPLICIT NONE
	    INTEGER, INTENT(IN) :: index
	    getTSNumPTs =  TimeSeriesPts(index)
	END FUNCTION

	SUBROUTINE getInterpRatingCurveValue(index, position, value)
	    IMPLICIT NONE
	    INTEGER, INTENT(IN) :: index
	    REAL, INTENT(IN) :: position
	    REAL, INTENT(OUT) :: value
	    INTEGER :: numSeg, i
        integer, parameter :: dim_num = 2
        REAL(kind = 16), DIMENSION(2) :: p1, p2, q1, q2, v
        TYPE(Segment) :: seg
        TYPE(Point) :: pnt
        INTEGER :: retval
	    seg = segment(point(0,0), point(0,0))
	    p1(1) = position
	    p1(2) = RCBBx(index)%ymin
	    pnt%x = position
	    pnt%y = RCBBx(index)%ymin
	    
	    p2(1) = position
	    p2(2) = RCBBx(index)%ymax
	    
	    DO i = 1,getRCNumSegments(index)
	        CALL getRCSegment(index, i, seg)
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

	SUBROUTINE getInterpTimeSeriesValue(index, position, value)
	    IMPLICIT NONE
	    INTEGER, INTENT(IN) :: index
	    REAL, INTENT(IN) :: position
	    REAL, INTENT(OUT) :: value
	    INTEGER :: numSeg, i
        REAL(kind = 16), DIMENSION(2) :: p1, p2, q1, q2, v
        TYPE(Segment) :: seg
        TYPE(Point) :: pnt
        INTEGER :: retval
	    seg = segment(point(0,0), point(0,0))
	    p1(1) = position
	    p1(2) = TSBBx(index)%ymin
	    pnt%x = position
	    pnt%y = TSBBx(index)%ymin
	    
	    p2(1) = position
	    p2(2) = TSBBx(index)%ymax
	    
	    DO i = 1,getTSNumSegments(index)
	        CALL getTSSegment(index, i, seg)
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
		
	SUBROUTINE setRatingCurveNumPts(dim, val)
	    IMPLICIT NONE
	    INTEGER, INTENT(IN) :: dim
	    INTEGER, INTENT(IN) :: val(:)
	    INTEGER :: I, STATUS
	    
		ALLOCATE(RatingCurvePts(dim), STAT = status)
		ALLOCATE(RCStartPts(dim), STAT = status)
		ALLOCATE(RCBBx(dim), STAT = status)
		
		DO I = 1,dim
		    RatingCurvePts(I)= val(i)
		    if(i == 1) then
		        RCStartPts(I) = 1
		    ELSE
		        RCStartPts(I) = RCStartPts(i-1)+(val(i-1))
		    ENDIF
		ENDDO
		
		NumRatingCurves = dim
	END SUBROUTINE
	
	SUBROUTINE setTimeSeriesNumPts(dim, val)
	    IMPLICIT NONE
	    INTEGER, INTENT(IN) :: dim
	    INTEGER, INTENT(IN) :: val(:)
	    INTEGER :: I, STATUS
	    
		ALLOCATE(TimeSeriesPts(dim), STAT = status)
		ALLOCATE(TSStartPts(dim), STAT = status)
		ALLOCATE(TSBBx(dim), STAT = status)
	
		DO I = 1,dim
		    TimeSeriesPts(I)= val(i)
		    if(i==1) then
		        TSStartPts(i) = 1
		    else
		        TSStartPts(i) = TSStartPts(i-1)+(val(i-1))
		    endif
		ENDDO
		
		NumTimeSeries = dim
	END SUBROUTINE
	
	SUBROUTINE setTimeSeries(dim, val)
	    IMPLICIT NONE
	    INTEGER, INTENT(IN) :: dim
	    REAL, INTENT(IN) :: val(:)
	    INTEGER :: I, STATUS, count
	    
		ALLOCATE(TimeSeries(dim), STAT = status)
		ALLOCATE(TSPts((dim/2)), STAT = status)
		
		DO I = 1,dim/2
		    count = (i-1)*2
		    TSPts(i)%x = val(count+1)
		    TSPts(i)%y = val(count+2)
		ENDDO
		    
		DO I = 1,dim
		    TimeSeries(I)= val(i)
		ENDDO

	END SUBROUTINE
	
	SUBROUTINE setRatingCurves(dim, val)
	    IMPLICIT NONE
	    INTEGER, INTENT(IN) :: dim
	    REAL, INTENT(IN) :: val(:)
	    INTEGER :: I, STATUS, count
	    
		ALLOCATE(RatingCurves(dim), STAT = status)
		ALLOCATE(RCPts((dim/2)), STAT = status)
		
		DO I = 1,dim/2
		    count = (i-1)*2
		    RCPts(i)%x = val(count+1)
		    RCPts(i)%y = val(count+2)
		ENDDO
		
		DO I = 1,dim
		    RatingCurves(I)= val(i)
		ENDDO
		
	END SUBROUTINE
	
	SUBROUTINE setTimeSeriesType(dim, val)
	    IMPLICIT NONE
	    INTEGER, INTENT(IN) :: dim
	    INTEGER, INTENT(IN) :: val(:)
	    INTEGER :: I, STATUS
	    
		ALLOCATE(TimeSeriesType(dim), STAT = status)
		
		DO I = 1,dim
		    TimeSeriesType(I)= val(i)
		ENDDO
		
	END SUBROUTINE
	
	SUBROUTINE setRatingCurveType(dim, val)
	    IMPLICIT NONE
	    INTEGER, INTENT(IN) :: dim
	    INTEGER, INTENT(IN) :: val(:)
	    INTEGER :: I, STATUS
	    
		ALLOCATE(RatingCurveType(dim), STAT = status)
		
		DO I = 1,dim
		    RatingCurveType(I)= val(i)
		ENDDO
		
	END SUBROUTINE
	
	SUBROUTINE dealloc_TimeSeries()
        IMPLICIT NONE
        INTEGER :: status
		DEALLOCATE(TimeSeriesType, STAT = status)
		DEALLOCATE(TimeSeries, STAT = status)
		DEALLOCATE(TSPts, STAT = status)
		DEALLOCATE(TimeSeriesPts, STAT = status)
		DEALLOCATE(TSStartPts, STAT = status)
		DEALLOCATE(TSBBx, STAT = status)

    END SUBROUTINE
    
    SUBROUTINE dealloc_RatingCurves()
        IMPLICIT NONE
        INTEGER :: status
		DEALLOCATE(RatingCurveType, STAT = status)
		DEALLOCATE(RatingCurves, STAT = status)
		DEALLOCATE(RCPts, STAT = status)
		DEALLOCATE(RatingCurvePts, STAT = status)
		DEALLOCATE(RCStartPts, STAT = status)
		DEALLOCATE(RCBBx, STAT = status)
        
    END SUBROUTINE
	
END MODULE
