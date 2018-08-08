    module fastmech

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
    tmp = 1
    
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