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
    !USE RivStagr4Mod_jmn
    use bmifastmech
    implicit none
    type (bmi_fastmech) :: fm
    INTEGER(4) count, num, i, status, cptArg, s
    logical :: lookForBMI=.FALSE.
    logical :: fileExist
    CHARACTER(LEN=250) buf
    CHARACTER(LEN=250) file
    double precision :: time0,time1,time
    
    !Check if arguments are found
    count = COMMAND_ARGUMENT_COUNT()

    if(count > 0) then
        !      CALL GETARG(1, buf, status)
        do cptArg=1,count
        CALL get_command_argument(cptArg, buf) !gnu fortran only take 2 args
        select case(buf)
        case('--BMI')
            lookForBMI = .TRUE.
            write(*,*)'BMI = ', lookForBMI
        case default
            inquire(file=adjustl(buf), exist=fileExist)
            if(.not.fileExist)then
                write(*,*)'file ', buf, ' not found'
                stop
            else
                file=adjustl(buf)
                write(*,*)'file name: ', file
                write (*,"(a)",advance="no") "Initializing..."
                s = fm%initialize(file)
                write (*,*) "Done."
                s = fm%get_start_time(time0)
                write (*,"(a30, f8.2)") "Start time:", time0
                s = fm%get_end_time(time1)
                write (*,"(a30, f8.2)") "End time:", time1
                s = fm%get_current_time(time)
                write (*,"(a30, f8.2)") "Current time:", time
                write (*,"(a)") "Running steady version"
                s = fm%update()
                s = fm%get_current_time(time)
                write (*,*) "Done with steady solution"
                if(time < time1) then
                    write (*,"(a)") "Running time solution"
                endif
                do while (time < time1)
                    s = fm%update()
                    s = fm%get_current_time(time)
                end do
                write (*,"(a)", advance="no") "Finalizing..."
                s = fm%finalize()
                write (*,*) "Done"

                endif
            !call stagr4(buf)
            write(*,*)'file name: ', file
        end select
        enddo
        if(lookForBMI) then
            !inquire(file=file, exist=fileExist)
            !if(.not.fileExist)then
            !    write(*,*)'file ', buf, ' not found'
            !    stop
            !endif
            write(*,*)'--BMI Arguement not needed'
        endif
    endif

    end program Stagr_Console

