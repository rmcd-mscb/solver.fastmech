    SUBROUTINE INITDT( Fs , Fsmin , hmax , lsub &
        , porosity , Lsubactdepth , Vsandsurfmax &
        , h , area , AV , BV , dro , Droughmin &
        , Droughmax , P , Vsandsubmax , Vsandsub &
        , Vsandsurf , dref)
    !
    !	CHANGED 7/20/2007
    !
    !	RUNS ONLY ONCE BEFORE SEDIMENT TRANSPORT CALCULATION BEGAN
    !	IT INIZIALIZES THE VARIABLES NEEDED IN THE SEDIMENT TRASNPORT
    !	CALCULATIONS
    !
    !     IF FS>FSMIN 	Vsandsub =	Vsandsubmax
    !
    !	INPUTS: Fs, Fsmin, hmax, lsub, porosity, Lsubactdepth, Vsandsurfmax
    !			h, ds, dn, av, bv, dro, Droughmin, Droughmax, P
    !
    !	OUTPUTS: Vsandsubmax, Vsandsub, Vsandsurf, dref
    !
    !
    !               VARIABLES
    !
    !
    !	Fs = sand fraction
    !	Fsmin = minimum sand fraction
    !	hmax = maximum height of gravel
    !	lsub = sand depth within the subsurface active layer
    !	porosity = porosity of the subsurface active layer
    !	Vsandsurfmax = minimum sand volume that covers the gravel
    !	h = height of the sand
    !	area = area of cell
    !     Droughmin = minimum roughness when sand cover all
    !
    !	tol = accuracy for the beta model
    !	IFLAG = error number
    !
    !	Vsandsubmax = maximum volume of sand in the subsurface active layer
    !	Vsandsub = volume of sand in the subsurface
    !	Vsandsurf = volume of sand over the sand
    !
    !!-----------------------------------------
    ! we used a logistic curve of shape y=(A1-A2)/(1+(Fs/x0)^P)+A2
    ! where A1 and A2 are the roughnesses for Fs=0 and Fs=1 respectively
    ! x0 is the location parameter which is set to 0.5,
    !    where a rapid change in roughness
    !    behavior happens
    ! p is the shape paramenter which we set to 7
    ! Fs is the areal fraction of sand
    !
    !
    IMPLICIT NONE
    INTEGER, PARAMETER :: mp = KIND(1.0D0)
    INTEGER i , IFLAG
    REAL(kind = mp) Fs, Fsmin, hmax, lsub, porosity, Lsubactdepth
    REAL(kind = mp) Vsandsubmax , Vsandsurfmax , Vsandsub , Vsandsurf
    REAL(kind = mp) area, DV, TOL, AV, BV, H, XX, dro, dref, Droughmin
    REAL(kind = mp) Droughmax, x0, P
    REAL(kind = mp) roughnessinv



    Vsandsubmax = porosity*Lsubactdepth*area



    IF (Fs.LE.Fsmin) THEN

        Vsandsub = porosity*Lsub*area
        Vsandsurf=0.0
        !		dref=Dro

    ELSEIF (Fs.GE.1.0) THEN

        Vsandsub = Vsandsubmax
        DV = h-hmax
        DV = DV*area
        Vsandsurf = Vsandsurfmax+DV
        !	    dref=Dsand

    ELSEIF (Fs.GE.0.99.AND.Fs.LT.1.0) THEN

        Vsandsurf = Vsandsurfmax
        !		dref=Dsand
        Vsandsub = Vsandsubmax

    ELSE

        tol=0.00001
        CALL PPFBET (Fs, av, bv, tol, IFLAG, XX )

        IF (IFLAG.GE.1) THEN

            write(*,*) 'Error in Beta evaulation of area'
            write(*,*) 'Error:', IFLAG
            write(*,*) 'in subroutine CDFBET'

        ENDIF

        DV=XX
        Vsandsurf = DV*Vsandsurfmax
        !		dref=hmax/(hmax-h)*dro
        Vsandsub = Vsandsubmax

    ENDIF


    IF (h.GE.hmax) THEN

        dref=Droughmax

    ELSEIF (Fs.LE.Fsmin) THEN

        dref=Dro

    ELSE

        x0=0.5
        !		P=7.0
        dref = roughnessinv(Droughmin, dro, P, Fs)
        !		dref=(Dro-Droughmin)*(1+(Fs/x0)**P)+Droughmin

    ENDIF



    RETURN

    END
