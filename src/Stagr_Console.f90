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
    USE IRICMI
    implicit none
    INTEGER(4) count, num, i, status
    CHARACTER(LEN=250) buf
    count = COMMAND_ARGUMENT_COUNT()
    IF(count.lt.1) then
        buf='No CGNS file'
        write(0,*) buf
    else
        CALL GET_COMMAND_ARGUMENT(1, buf)
    END IF
    CALL STAGR4(buf)

    ! Variables


    ! Body of Stagr_Console

    end program Stagr_Console

