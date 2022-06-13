!  test1.f90 
!
!  FUNCTIONS:
!  test1 - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: test1
!
!  PURPOSE:  Test unlinking solutions from hdf5 file
!
!****************************************************************************

    program test1
    
    use hdf5
    use unlink_hdf5
    use iso_c_binding
    implicit none

    ! Variables
    integer        :: status
    integer(hid_t) :: file_id

    ! Body of test1
    call execute_command_line("copy /Y Case1.cgn.orig Case1.cgn")
    call h5open_f(status)
    call h5fopen_f("Case1.cgn", H5F_ACC_RDWR_F, file_id, status)
    call unlink_solutions(file_id, status)
    call h5fclose_f(file_id, status)
    call h5close_f(status)

    !call exit(status)   ! THIS CAUSES THE HDF5 FILE TO NOT BE SAVED

    end program test1
