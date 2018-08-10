    module fastmech
    USE EInitMod
    USE RivVarMod
    USE CalcCond
    USE CSedMod_DT_SUSP
    USE CSedMod
    USE ReadCGNS
    USE RivVarWMod
    USE RivVarVertMod
    USE WriteCGNS
    USE RivVertMod
    USE StressDivMod
    !	USE CSedMod
    USE UInitMod
    USE TRIDIAG

    USE RivRoughnessMod
    USE RivVarTimeMod
    USE RivConnectivityMod
    implicit none

    type :: fastmech_model
        real :: itermax
    end type fastmech_model

    private :: initialize

    contains

    ! Initializes the model with values read from a file.
    subroutine initialize_from_file(model, config_file)
    character (len=*), intent (in) :: config_file
    type(fastmech_model), intent(out) :: model
    integer :: tmp
    real(kind=mp) :: newDisch, newStage
    noldwetnodes = 0
    nwetnodes = 0
    errorcode = 0
    CALL WELCOME
    CALL read_iRIC_CGNS(config_file)
    CALL Calc_Area(x, y, xo, yo, nm, dn, harea)
    if(varDischType == 1) then !Discharge Time Series
        CALL getInterpTimeSeriesValue(1, VarDischStartTime, newDisch)
        q = newDisch*1e6
    else
        newDisch = q/1e6
    endif
    if(varStageType == 1) then !Stage Time Series
        !in this case the times series is stored in the rating curve
        CALL getInterpRatingCurveValue(1, VarDischStartTime, newStage)
    else if (varStageType == 2) then !Stage Rating Curve
        Call getInterpRatingCurveValue(1, newDisch, newStage)
    endif
    if(.not.HotStart) then
        CALL einit(e, hl, eta, ibc, w, hav)
        call uinit(u,v,hav,w,eta,q,ibc)
        if(errorcode < 0) then
            !call write_error_CGNS(STR_IN)
            CALL dealloc_common2d()
            if(vbc == 1) then
                CALL DEALLOC_VELBC()
            endif
            if(CALCQUASI3D) then
                CALL dealloc_common3d()
                if(TRANSEQTYPE == 1) THEN
                    CAll dealloc_csed3d_dt()
                ENDIF
            endif
            CALL dealloc_working()
            !pause
            return
        endif
    endif
    if(roughnesstype == 1) then !Roughness as Z0
        !		CALL ZOTOCDTWO() !Calculate cd from log profile
        CALL Z0TOCD() !Calculate cd from 2-part profile
    endif
    nsteps = (VarDischEndTime-VarDischStartTime)/dt

    if(IterationOut) Then
        if(itm > nsteps) then
            CALL alloc_TSNames(itm)
        else
            CALL alloc_TSNames(nsteps+10)
        endif
    else
        CALL alloc_TSNames(nsteps+10)
    endif
    oldStage = hav(ns)
    dsstage = hav(ns)
    !    ENDIF
    newStage = oldStage
    oldDisch = q
    
    nct = -1
    totTime = VarDischStartTime
    ptime = totTime+iplinc*dt
    
    end subroutine initialize_from_file

    ! Initializes the model with default hardcoded values.
    subroutine initialize_from_defaults(model)
    type (fastmech_model), intent (out) :: model
    integer :: tmp
    tmp = 1
    !model%alpha = 0.75
    !model%t_end = 20.
    !model%n_x = 10
    !model%n_y = 20
    !call initialize(model)
    end subroutine initialize_from_defaults

    ! Allocates memory and sets values for either initialization technique.
    subroutine initialize(model)
    type (fastmech_model), intent (inout) :: model
    integer :: tmp
    tmp = 1

    !model%t = 0.
    !model%dt = 1.
    !model%dx = 1.
    !model%dy = 1.
    !
    !allocate(model%temperature(model%n_y, model%n_x))
    !allocate(model%temperature_tmp(model%n_y, model%n_x))
    !
    !model%temperature = 0.
    !model%temperature_tmp = 0.
    !
    !call set_boundary_conditions(model%temperature)
    !call set_boundary_conditions(model%temperature_tmp)
    end subroutine initialize

    end module fastmech