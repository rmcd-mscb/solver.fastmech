    MODULE RivRoughnessMod2
    USE RivVarTimeMod2
    USE RivVarMod2
    USE CalcCond2
    IMPLICIT NONE
    !	REAL, ALLOCATABLE, DIMENSION(:,:) :: ribc
    type rivrough
        REAL(KIND=mp), ALLOCATABLE, DIMENSION(:) :: ribcvals
        INTEGER :: numRoughnessRegions
    end type
    !! ribc:         Is an array which specifies the region of the node
    !! ribcvals:     Is the value of the region
    !!  Can either be constant or a rating curve.  If rating curve
    !!  then value is negative and its integer value represents the
    !!  rating curve to use.

    CONTAINS
    SUBROUTINE setRoughness(ns, nn, roughnessType, q, cd, znaught)
    INTEGER, INTENT(IN) :: ns, nn, roughnessType
    REAL(KIND=mp), INTENT(IN) :: q
    REAL(KIND=mp), INTENT(INOUT) :: cd(:,:), znaught(:,:)
    INTEGER :: i,j, tmpcd
    REAL(KIND=mp) :: ratingval

    DO i = 1,ns
        DO j = 1, nn
            if(ribc(i,j) >= 0) then
                !	                write(6,*) ribc(i,j), ribcvals(ribc(i,j))
                if(DMod(cd(i,j), 1.0D0) == 0) then
                    ratingval = ribcvals(ribc(i,j)+1)
                else if (roughnessType == 0) then
                    ratingval = cd(i,j)
                else if (roughnessType == 1) then
                    ratingval = znaught(i,j)
                endif

            else
                tmpcd = ABS(ribcvals(ribc(i,j)+1))
                CALL getInterpRatingCurveValue(tmpcd, q, ratingval)
            endif

            SELECT CASE(roughnessType)
            CASE(0)
                cd(i,j) = ratingval
            CASE(1)
                znaught(i,j) = ratingval
            END SELECT

        ENDDO
    ENDDO

    END SUBROUTINE

    !	    SUBROUTINE setRoughness(ns, nn, roughnessType, q, cd, znaught)
    !	        INTEGER, INTENT(IN) :: ns, nn, roughnessType
    !	        REAL, INTENT(IN) :: q
    !	        REAL, INTENT(INOUT) :: cd(:,:), znaught(:,:)
    !	        INTEGER :: i,j, tmpcd
    !	        REAL :: ratingval
    !
    !	        DO i = 1,ns
    !	            DO j = 1, nn
    !	                if(ribc(i,j) >= 0) then
    !!	                write(6,*) ribc(i,j), ribcvals(ribc(i,j))
    !	                    ratingval = ribcvals(ribc(i,j)+1)
    !	                else
    !	                    tmpcd = ABS(ribcvals(ribc(i,j)))
    !	                    CALL getInterpRatingCurveValue(tmpcd, q, ratingval)
    !	                endif
    !
    !	                SELECT CASE(roughnessType)
    !	                CASE(0)
    !	                    cd(i,j) = ratingval
    !	                CASE(1)
    !	                    znaught(i,j) = ratingval*100. !convert to cm from m
    !	                END SELECT
    !
    !	            ENDDO
    !	        ENDDO
    !
    !	    END SUBROUTINE

    !      SUBROUTINE alloc_roughness(ns2, nn)
    !              INTEGER, INTENT(IN) :: ns2, nn
    !              INTEGER :: status
    !
    !              ALLOCATE(ribc(ns, nn), STAT = status)
    !
    ! 		END SUBROUTINE

    SUBROUTINE alloc_roughnessvals(this, dim)
    implicit none
    type(rivrough), intent(inout) :: this
    INTEGER, INTENT(IN) :: dim
    INTEGER :: status
    this%numRoughnessRegions = dim
    ALLOCATE(this%ribcvals(dim), STAT = status)

    END SUBROUTINE

    SUBROUTINE dealloc_roughness(this)
    implicit none
    type(rivrough), intent(inout) :: this
    INTEGER :: status
    !			DEALLOCATE(ribc, STAT = status)
    DEALLOCATE(this%ribcvals, STAT = status)
    END SUBROUTINE

    SUBROUTINE setRoughnessIBC(val)
    IMPLICIT NONE
    REAL(KIND=mp), INTENT(IN) :: val(:)
    INTEGER :: i, j, count, countji
    DO i = 1, ns2
        DO j = 1, nn
            countji = ((j-1)*ns2)+i
            count = ((i-1)*nn)+j
            !					ribc2(i,j) = val(countji)
        END DO
    END DO
    END SUBROUTINE

    SUBROUTINE setRoughnessIBCVals(val)
    IMPLICIT NONE
    REAL(KIND=mp), INTENT(IN) :: val(:)
    INTEGER :: i, j, count, status
    DO i = 1, numRoughnessRegions
        SELECT CASE(roughnessType)
        CASE(0)
            !	                cd(i,j) = ratingval
            ribcvals(i) = val(i)
        CASE(1)
            !	                znaught(i,j) = ratingval
            ribcvals(i) = val(i)*100.
        END SELECT
    END DO
    END SUBROUTINE

    END MODULE
