    SUBROUTINE ROUGHDT(Fs, Fsmin, droughmin, dref, hmax, h, P, drough)

    !	CHANGED 7/19/2007
    !
    !	CALCULATE THE REDUCTION OF BED ROUGHNESS
    !	WITH AS A sigmoidal REDUCTION FROM SAND DEPOSITION
    !
    !	INPUT : Fs, Fsmin, droughmin, Dref, hmax, h, P
    !	OUTPUT: Drough
    !
    !	drough>=droughmin
    !	drough<=dref
    !
    !	VARIABLES
    !
    !	Fs = fraction of sand
    !	Fsmin = minimal fraction of sand
    !	Dref = referece diameter of roughness
    !     droughmin= roughness due to the fine sediment class
    !	hmax = maximum protrusion of gravel
    !	h = depth of sand relative to original topography
    !	Drough = diameter for the roughness of the bed
    !
    IMPLICIT NONE
    INTEGER, PARAMETER :: mp = KIND(1.0D0)
    REAL(kind = mp) Fs, Fsmin, droughmin, dref
    REAL(kind = mp) hmax, h
    REAL(kind = mp) drough, m, roughness
    REAL(kind = mp) P

    IF (Fs.LE.Fsmin) THEN

        drough=dref
        GOTO 333

    ELSEIF (Fs.GE.0.99) THEN

        drough=droughmin
        GOTO 333

    ELSE

        IF (h.LE.0.0) THEN

            drough=dref
            GOTO 333

        ENDIF


        drough=roughness(dref, droughmin ,P ,Fs)

        IF (drough.LE.droughmin) THEN

            drough=droughmin

        ENDIF



    ENDIF

333 RETURN

    END