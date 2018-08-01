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
    implicit none
    INTEGER(4) count, num, i, status
    CHARACTER(LEN=250) buf
    !      count = NARGS()

    !      CALL GETARG(1, buf, status)
    CALL GETARG(1, buf) !gnu fortran only take 2 args
    IF (status .lt. 0) THEN
        WRITE (*,*) 'GETARG error - exiting'
        !            EXIT
    END IF
    CALL STAGR4(buf)

    ! Variables


    ! Body of Stagr_Console

    end program Stagr_Console

