module unlink_hdf5

  use hdf5
  use iso_c_binding
  implicit none
  
  integer(hid_t), private :: g_file_id
  integer(hsize_t), dimension(1:2) :: g_maxdims

  
  ! Not in /iRIC/iRICZone/ZoneIterativeData
  ! ParticleGroupSolutionPointers
  ! PolydataSolutionPointers
  
  logical, private :: g_particle_group_solution
  logical, private :: g_polydata_solution

contains

!!************************************************************
!!
!!  Operator function.  unlinks each given object by name
!!
!! ************************************************************
!
!  integer function op_cb(loc_id, name, info, op_data) bind(C)
!
!    use hdf5
!    use iso_c_binding
!    implicit none
!
!    integer(hid_t), VALUE :: loc_id
!    character(len=1), dimension(1:255) :: name ! must have LEN=1 for bind(C) strings
!    type(c_ptr) :: info
!    type(c_ptr) :: op_data
!
!    integer :: status, i, size
!    character(len=255) :: name_string
!    logical :: link_exists
!
!    size = 0
!    do i = 1, 255
!        if(name(i)(1:1).EQ.C_NULL_CHAR) exit
!        name_string(i:i) = name(i)(1:1)
!        size = i
!    enddo
!
!    call h5lexists_f(loc_id, name_string(1:size), link_exists, status)
!    if (link_exists) then
!        call h5ldelete_f(loc_id, name_string(1:size), status)
!    endif
!
!    op_cb = 0
!
!  end function op_cb
!
!!************************************************************
!!
!!  Strips output groups -- unlinks groups, file size doesn't
!!  change
!!
!! ************************************************************
!
!  subroutine strip_hdf5(file_id, ier)
!    use hdf5
!    use iso_c_binding
!    implicit none
!
!    integer(hid_t) :: file_id
!    integer :: ier
!    integer(hid_t) :: grp_id
!    logical :: link_exists
!
!    integer :: status
!
!    type(c_funptr) :: funptr
!    type(c_ptr) :: ptr
!    integer(hsize_t) :: idx
!    integer :: ret_value
!
!    ! /iRIC3D
!    call h5lexists_f(file_id, "iRIC3D", link_exists, status)
!    if (link_exists) then
!        call h5ldelete_f(file_id, "iRIC3D", status)
!    endif
!
!    ! /iRIC
!    call h5lexists_f(file_id, "/iRIC", link_exists, status)
!    if (link_exists) then
!        ! /iRIC/BaseIterativeData
!        call h5lexists_f(file_id, "/iRIC/BaseIterativeData", link_exists, status)
!        if (link_exists) then
!            CALL h5ldelete_f(file_id, "/iRIC/BaseIterativeData", status)
!        endif
!        ! /iRIC/iRICZone
!        call h5lexists_f(file_id, "/iRIC/iRICZone", link_exists, status)
!        if (link_exists) then
!            ! /iRIC/iRICZone/FlowCellSolution1
!            call h5lexists_f(file_id, "/iRIC/iRICZone/FlowCellSolution1", link_exists, status)
!            if (link_exists) then
!                call h5ldelete_f(file_id, "/iRIC/iRICZone/FlowCellSolution1", status)
!            endif
!            ! /iRIC/iRICZone/FlowIFaceSolution1
!            call h5lexists_f(file_id, "/iRIC/iRICZone/FlowIFaceSolution1", link_exists, status)
!            if (link_exists) then
!                CALL h5ldelete_f(file_id, "/iRIC/iRICZone/FlowIFaceSolution1", status)
!            endif
!            ! /iRIC/iRICZone/FlowJFaceSolution1
!            call h5lexists_f(file_id, "/iRIC/iRICZone/FlowJFaceSolution1", link_exists, status)
!            if (link_exists) then
!                call h5ldelete_f(file_id, "/iRIC/iRICZone/FlowJFaceSolution1", status)
!            endif
!            ! /iRIC/iRICZone/FlowSolution1
!            call h5lexists_f(file_id, "/iRIC/iRICZone/FlowSolution1", link_exists, status)
!            if (link_exists) then
!                call h5ldelete_f(file_id, "/iRIC/iRICZone/FlowSolution1", status)
!            endif
!            ! /iRIC/iRICZone/GridCoordinatesForSolution1
!            call h5lexists_f(file_id, "/iRIC/iRICZone/GridCoordinatesForSolution1", link_exists, status)
!            if (link_exists) then
!                CALL h5ldelete_f(file_id, "/iRIC/iRICZone/GridCoordinatesForSolution1", status)
!            endif
!            ! /iRIC/iRICZone/ParticleGroupSolution1
!            call h5lexists_f(file_id, "/iRIC/iRICZone/ParticleGroupSolution1", link_exists, status)
!            iF (link_exists) then
!                call h5ldelete_f(file_id, "/iRIC/iRICZone/ParticleGroupSolution1", status)
!            endif
!            ! /iRIC/iRICZone/PolydataSolution1
!            call h5lexists_f(file_id, "/iRIC/iRICZone/PolydataSolution1", link_exists, status)
!            iF (link_exists) then
!                call h5ldelete_f(file_id, "/iRIC/iRICZone/PolydataSolution1", status)
!            endif
!            ! /iRIC/iRICZone/ZoneIterativeData/*
!            call h5lexists_f(file_id, "/iRIC/iRICZone/ZoneIterativeData", link_exists, status)
!            if (link_exists) then
!                ! CALL h5ldelete_f(file_id, "/iRIC/iRICZone/ZoneIterativeData", status)
!                call h5gopen_f(file_id, "/iRIC/iRICZone/ZoneIterativeData", grp_id, status)
!                idx    = 0
!                funptr = C_FUNLOC(op_cb)
!                ptr    = C_NULL_PTR
!                call H5Literate_f(grp_id, H5_INDEX_NAME_F, H5_ITER_NATIVE_F, idx, funptr, ptr, ret_value, status)
!           endif
!        endif
!    endif
!  end subroutine strip_hdf5
  
!************************************************************
!
!  Operator function.  unlinks each given object by name
!
! ************************************************************

  !!INTEGER FUNCTION pointers_cb(loc_id, name, info, op_data) bind(C)
  integer function pointers_cb(loc_id, name, info, file_id) bind(C)

    use hdf5
    use iso_c_binding
    implicit none

    integer(hid_t), value :: loc_id
    integer(hid_t)        :: file_id
    integer(hid_t)        :: dset_id
    integer(hid_t)        :: type_id
    integer(hid_t)        :: space_id
    character(len=1), dimension(1:255) :: name ! must have LEN=1 for bind(C) strings
    type(c_ptr) :: info
    
    integer(hsize_t), dimension(1:2) :: dims
    integer :: rankr
    
    integer :: status, i, size
    character(len=255) :: name_string
    logical :: link_exists

    character(len=20) :: num_str

    ! works w/ intel fortran on windows
    if (file_id .ne. g_file_id) then
        print "(A)", "Warning: file_id doesn't match"
    endif

    size = 0
    do i = 1, 255
        if(name(i)(1:1) .EQ. c_null_char) exit
        name_string(i:i) = name(i)(1:1)
        size = i
    enddo

    call h5lexists_f(loc_id, name_string(1:size), link_exists, status)
    if (link_exists) then
        call h5dopen_f(loc_id, name_string(1:size) // "/ data" , dset_id, status)

        call h5dget_space_f(dset_id, space_id, status)
        call h5sget_simple_extent_dims_f(space_id, dims, g_maxdims, status)
        call h5sclose_f(space_id, status)

        do i = 1, g_maxdims(2)
            write(num_str, "(I)"), i
            if (name_string(1:size) .EQ. "GridCoordinatesPointers") then
                ! Note /iRIC/iRICZone/ZoneIterativeData/GridCoordinatesPointers -> GridCoordinatesForSolution1 .. GridCoordinatesForSolutionN
                call h5ldelete_f(g_file_id, "/iRIC/iRICZone/" // name_string(1:size-8) // "ForSolution" // trim(adjustl(num_str)), status)
            else
                call h5ldelete_f(g_file_id, "/iRIC/iRICZone/" // name_string(1:size-8) // trim(adjustl(num_str)), status)
            endif
            
            ! remove solutions not stored in /iRIC/iRICZone/ZoneIterativeData
            if (g_particle_group_solution) then
                call h5ldelete_f(g_file_id, "/iRIC/iRICZone/ParticleGroupSolution" // trim(adjustl(num_str)), status)
            endif
            if (g_polydata_solution) then
                call h5ldelete_f(g_file_id, "/iRIC/iRICZone/PolydataSolution" // trim(adjustl(num_str)), status)
            endif
        enddo

        g_particle_group_solution = .FALSE.
        g_polydata_solution = .FALSE.

        call h5dclose_f(dset_id, status)
        call h5ldelete_f(loc_id, name_string(1:size), status)
    endif

    pointers_cb = 0

  end function pointers_cb


!************************************************************
!
!  Strips output groups -- unlinks groups, file size doesn't
!  change
!
! ************************************************************

  subroutine unlink_solutions(file_id, ier)

    use hdf5
    use iso_c_binding
    implicit none

    integer(hid_t) :: file_id     
    !! integer(hid_t) :: fid

    integer :: ier
    integer(hid_t) :: grp_id
    logical :: link_exists

    integer :: status

    type(c_funptr) :: funptr
    type(c_ptr) :: ptr
    integer(hsize_t) :: idx
    integer :: ret_value

    g_file_id = file_id
    g_particle_group_solution = .TRUE.
    g_polydata_solution = .TRUE.    

    ! /iRIC
    call h5lexists_f(file_id, "/iRIC", link_exists, status)
    if (link_exists) then

        ! /iRIC/BaseIterativeData
        call h5lexists_f(file_id, "/iRIC/BaseIterativeData", link_exists, status)
        if (link_exists) then
            CALL h5ldelete_f(file_id, "/iRIC/BaseIterativeData", status)
        endif

        ! /iRIC/iRICZone
        call h5lexists_f(file_id, "/iRIC/iRICZone", link_exists, status)
        if (link_exists) then
            
            ! /iRIC/iRICZone/ZoneIterativeData/*
            call h5lexists_f(file_id, "/iRIC/iRICZone/ZoneIterativeData", link_exists, status)
            if (link_exists) then
                call h5gopen_f(file_id, "/iRIC/iRICZone/ZoneIterativeData", grp_id, status)
                idx    = 0
                funptr = C_FUNLOC(pointers_cb)
                !!ptr    = C_NULL_PTR
                ptr    = C_LOC(file_id)
                call h5literate_f(grp_id, H5_INDEX_NAME_F, H5_ITER_NATIVE_F, idx, funptr, ptr, ret_value, status)
           endif
        endif
    endif

    ! /iRIC3D
    call h5lexists_f(file_id, "iRIC3D", link_exists, status)
    if (link_exists) then
        call h5ldelete_f(file_id, "iRIC3D", status)
    endif

  end subroutine unlink_solutions
    
  !subroutine TestStripHDF5()
  !  use hdf5
  !  use iso_c_binding
  !  implicit none
  !
  !  integer(hid_t) :: file_id
  !  integer :: ier
  !  integer(hid_t) :: grp_id
  !  logical :: link_exists
  !
  !  integer :: status
  !
  !  type(c_funptr) :: funptr
  !  type(c_ptr) :: ptr
  !  integer(hsize_t) :: idx
  !  integer :: ret_value
  !
  !  call execute_command_line("copy /Y Case1.cgn.orig Case1.cgn")
  !  call h5open_f(status)
  !  call h5fopen_f("Case1.cgn", H5F_ACC_RDWR_F, file_id, status)
  !  call strip_hdf5(file_id, status)
  !  call h5fclose_f(file_id, status)
  !  call h5close_f(status)
  !
  !end subroutine TestStripHDF5

end module unlink_hdf5
