    !  Stagr_Console.f90
    !
    !  FUNCTIONS:
    !  Stagr_Console      - Entry point of console application.
    !

    !****************************************************************************
    !
    !  PROGRAM: Stagr_Console
    !
    !  PURPOSE:  Entry point for the console application.
    !
    !****************************************************************************

    program Stagr_Console
    USE RivStagr4Mod_jmn
    use RivStagrMod_bmi
    implicit none
    INTEGER(4) count, num, i, status, cptArg
    logical :: lookForBMI=.FALSE.
    logical :: fileExist
    CHARACTER(LEN=250) buf
    
    !Check if arguments are found
    count = COMMAND_ARGUMENT_COUNT()

    if(count > 0) then
        !      CALL GETARG(1, buf, status)
        do cptArg=1,count
        CALL get_command_argument(1, buf) !gnu fortran only take 2 args
        select case(adjustl(buf))
        case("--BMI")
            lookForBMI = .TRUE.
        case default
            inquire(file=buf, exist=fileExist)
            if(.not.fileExist)then
                write(*,*)'file ', buf, ' not found'
                stop
            endif
            call stagr4(buf)
        end select
        if(lookForBMI) then
            inquire(file=buf, exist=fileExist)
            if(.not.fileExist)then
                write(*,*)'file ', buf, ' not found'
                stop
            endif
            call STAGRBMI(buf)
        endif
        end do
    endif

    ! Variables


    ! Body of Stagr_Console

    end program Stagr_Console

