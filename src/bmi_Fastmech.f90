    module bmifastmech

    use fastmech
    use bmif
    use, intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_f_pointer
    implicit none

    type, extends (bmi) :: bmi_fastmech
        private
        type (fastmech_model) :: model
    contains
    procedure :: get_component_name => fm_component_name
    procedure :: get_input_var_names => fm_input_var_names
    procedure :: get_output_var_names => fm_output_var_names
    procedure :: initialize => fm_initialize
    procedure :: finalize => fm_finalize
    procedure :: get_start_time => fm_start_time
    procedure :: get_end_time => fm_end_time
    procedure :: get_current_time => fm_current_time
    procedure :: get_time_step => fm_time_step
    procedure :: get_time_units => fm_time_units
    procedure :: update => fm_update
    procedure :: update_frac => fm_update_frac
    procedure :: update_until => fm_update_until
    procedure :: get_var_grid => fm_var_grid
    procedure :: get_grid_type => fm_grid_type
    procedure :: get_grid_rank => fm_grid_rank
    procedure :: get_grid_shape => fm_grid_shape
    procedure :: get_grid_size => fm_grid_size
    !procedure :: get_grid_spacing => fm_grid_spacing
    !procedure :: get_grid_origin => fm_grid_origin
    procedure :: get_grid_x => fm_grid_x
    procedure :: get_grid_y => fm_grid_y
    procedure :: get_grid_z => fm_grid_z
    !procedure :: get_grid_connectivity => fm_grid_connectivity
    !procedure :: get_grid_offset => fm_grid_offset
    procedure :: get_var_type => fm_var_type
    procedure :: get_var_units => fm_var_units
    procedure :: get_var_itemsize => fm_var_itemsize
    procedure :: get_var_nbytes => fm_var_nbytes
    procedure :: get_value_int => fm_get_int
    !procedure :: get_value_float => fm_get_float
    procedure :: get_value_double => fm_get_double
    generic :: get_value => &
        get_value_int, &
        !get_value_float, &
        get_value_double
    procedure :: get_value_ref_int => fm_get_ref_int
    !procedure :: get_value_ref_float => fm_get_ref_float
    procedure :: get_value_ref_double => fm_get_ref_double
    generic :: get_value_ref => &
        get_value_ref_int, &
        !get_value_ref_float, &
        get_value_ref_double
    procedure :: get_value_at_indices_int => fm_get_at_indices_int
    !procedure :: get_value_at_indices_float => fm_get_at_indices_float
    procedure :: get_value_at_indices_double => fm_get_at_indices_double
    generic :: get_value_at_indices => &
         get_value_at_indices_int, &
    !     get_value_at_indices_float, &
         get_value_at_indices_double
    procedure :: set_value_int => fm_set_int
    !procedure :: set_value_float => fm_set_float
    procedure :: set_value_double => fm_set_double
    generic :: set_value => &
         set_value_int, &
    !     set_value_float, &
         set_value_double
    procedure :: set_value_at_indices_int => fm_set_at_indices_int
    !procedure :: set_value_at_indices_float => fm_set_at_indices_float
    procedure :: set_value_at_indices_double => fm_set_at_indices_double
    generic :: set_value_at_indices => &
         set_value_at_indices_int, &
    !     set_value_at_indices_float, &
         set_value_at_indices_double
    procedure :: print_model_info
    end type bmi_fastmech

    private
    public :: bmi_fastmech

    character (len=BMI_MAX_COMPONENT_NAME), target :: &
        component_name = "Fastmech"

    ! Exchange items
    integer, parameter :: input_item_count = 3
    integer, parameter :: output_item_count = 10
    character (len=BMI_MAX_VAR_NAME), target, &
        dimension(input_item_count) :: input_items = (/ &
        'Elevation   ', &
        'roughness   ', &
        'vegroughness'/)
    character (len=BMI_MAX_VAR_NAME), target, &
        dimension(output_item_count) :: &
        output_items = (/ &
        'Depth                ', &
        'Drag_Coefficient     ', &
        'Elevation            ', &
        'FMIBC                ', &
        'IBC                  ', &
        'ShearStressX         ', &
        'ShearStressY         ', &
        'VelocityX            ', &
        'VelocityY            ', &
        'WaterSurfaceElevation'/)

    contains

    ! Get the name of the model.
    function fm_component_name(self, name) result (bmi_status)
    class (bmi_fastmech), intent(in) :: self
    character (len=*), pointer, intent(out) :: name
    integer :: bmi_status

    name => component_name
    bmi_status = BMI_SUCCESS
    end function fm_component_name

    ! List input variables.
    function fm_input_var_names(self, names) result (bmi_status)
    class (bmi_fastmech), intent(in) :: self
    character (*), pointer, intent(out) :: names(:)
    integer :: bmi_status

    names => input_items
    bmi_status = BMI_SUCCESS
    end function fm_input_var_names

    ! List output variables.
    function fm_output_var_names(self, names) result (bmi_status)
    class (bmi_fastmech), intent(in) :: self
    character (*), pointer, intent(out) :: names(:)
    integer :: bmi_status

    names => output_items
    bmi_status = BMI_SUCCESS
    end function fm_output_var_names

    !BMI initializer.
    function fm_initialize(self, config_file) result (bmi_status)
    class (bmi_fastmech), intent(out) :: self
    character (len=*), intent(in) :: config_file
    integer :: bmi_status

    if (len(config_file) > 0) then
        call initialize_from_file(self%model, config_file)
    else
        call initialize_from_defaults(self%model)
    end if
    bmi_status = BMI_SUCCESS
    end function fm_initialize

    ! BMI finalizer.
    function fm_finalize(self) result (bmi_status)
    class (bmi_fastmech), intent(inout) :: self
    integer :: bmi_status

    call cleanup(self%model)
    bmi_status = BMI_SUCCESS
    end function fm_finalize

    ! Model start time.
    function fm_start_time(self, time) result (bmi_status)
    class (bmi_fastmech), intent(in) :: self
    double precision, intent(out) :: time
    integer :: bmi_status

    if(self%model%t_calccond%soltype == 0) then
        time = 0.d0
    else
        time = dble(self%model%t_rivvartime%vardischstarttime)
    endif
    bmi_status = BMI_SUCCESS
    end function fm_start_time

    ! Model end time.
    function fm_end_time(self, time) result (bmi_status)
    class (bmi_fastmech), intent(in) :: self
    double precision, intent(out) :: time
    integer :: bmi_status

    time = dble(self%model%t_end)
    bmi_status = BMI_SUCCESS
    end function fm_end_time

    ! Model current time.
    function fm_current_time(self, time) result (bmi_status)
    class (bmi_fastmech), intent(in) :: self
    double precision, intent(out) :: time
    integer :: bmi_status

    time = dble(self%model%t)
    bmi_status = BMI_SUCCESS
    end function fm_current_time

    ! Model time step.
    function fm_time_step(self, time_step) result (bmi_status)
    class (bmi_fastmech), intent(in) :: self
    double precision, intent(out) :: time_step
    integer :: bmi_status

    time_step = dble(self%model%dt)
    bmi_status = BMI_SUCCESS
    end function fm_time_step

    ! Model time units.
    function fm_time_units(self, time_units) result (bmi_status)
    class (bmi_fastmech), intent(in) :: self
    character (len=*), intent(out) :: time_units
    integer :: bmi_status

    time_units = "s"
    bmi_status = BMI_SUCCESS
    end function fm_time_units

    ! Advance model by one time step.
    function fm_update(self) result (bmi_status)
    class (bmi_fastmech), intent(inout) :: self
    integer :: bmi_status

    call advance_in_time(self%model)
    bmi_status = BMI_SUCCESS
    end function fm_update

    ! Advance the model by a fraction of a time step.
    function fm_update_frac(self, time_frac) result (bmi_status)
    class (bmi_fastmech), intent(inout) :: self
    double precision, intent(in) :: time_frac
    integer :: bmi_status
    real :: time_step

    if (time_frac > 0.0) then
        time_step = self%model%dt
        self%model%dt = time_step*real(time_frac)
        call advance_in_time(self%model)
        self%model%dt = time_step
    endif
    bmi_status = BMI_SUCCESS
    end function fm_update_frac

    ! Advance the model until the given time.
    function fm_update_until(self, time) result (bmi_status)
    class (bmi_fastmech), intent(inout) :: self
    double precision, intent(in) :: time
    integer :: bmi_status
    double precision :: n_steps_real
    integer :: n_steps, i, s

    if (time > self%model%t) then
        n_steps_real = (time - self%model%t) / self%model%dt
        n_steps = floor(n_steps_real)
        do i = 1, n_steps
            s = self%update()
        end do
        s = self%update_frac(n_steps_real - dble(n_steps))
    end if
    bmi_status = BMI_SUCCESS
    end function fm_update_until

    ! Get the grid id for a particular variable.
    function fm_var_grid(self, var_name, grid_id) result (bmi_status)
    class (bmi_fastmech), intent(in) :: self
    character (len=*), intent(in) :: var_name
    integer, intent(out) :: grid_id
    integer :: bmi_status

    select case (var_name)
    case ('Elevation', 'roughness', 'vegroughness', &
        'Depth', 'Drag_Coefficient', 'FMIBC', &
        'IBC', 'ShearStressX', 'ShearStressY', &
        'WaterSurfaceElevation')
        grid_id = 0
        bmi_status = BMI_SUCCESS
    case('VelocityX', 'VelocityY')
        if(self%model%t_calccond%CALCQUASI3D) then
            grid_id = 1
        else
            grid_id = 0
        endif
        bmi_status = BMI_SUCCESS
        case default
        grid_id = -1
        bmi_status = BMI_FAILURE
    end select
    end function fm_var_grid

    ! The type of a variable's grid.
    function fm_grid_type(self, grid_id, grid_type) result (bmi_status)
    class (bmi_fastmech), intent(in) :: self
    integer, intent(in) :: grid_id
    character (len=*), intent(out) :: grid_type
    integer :: bmi_status

    select case (grid_id)
    case (0)
        grid_type = "structured_quadrilateral"
        bmi_status = BMI_SUCCESS
    case (1)
        grid_type = "structured_quadrilateral"
        bmi_status = BMI_SUCCESS
        case default
        grid_type = "-"
        bmi_status = BMI_FAILURE
    end select
    end function fm_grid_type

    ! The number of dimensions of a grid.
    function fm_grid_rank(self, grid_id, grid_rank) result (bmi_status)
    class (bmi_fastmech), intent(in) :: self
    integer, intent(in) :: grid_id
    integer, intent(out) :: grid_rank
    integer :: bmi_status

    select case (grid_id)
    case (0)
        grid_rank = 2
        bmi_status = BMI_SUCCESS
    case (1)
        grid_rank = 3
        bmi_status = BMI_SUCCESS
        case default
        grid_rank = -1
        bmi_status = BMI_FAILURE
    end select
    end function fm_grid_rank

    ! The dimensions of a grid.
    function fm_grid_shape(self, grid_id, grid_shape) result (bmi_status)
    class (bmi_fastmech), intent(in) :: self
    integer, intent(in) :: grid_id
    integer, dimension(:), intent(out) :: grid_shape
    integer :: bmi_status

    select case (grid_id)
    case (0)
        grid_shape = [self%model%n_x, self%model%n_y]
        bmi_status = BMI_SUCCESS
    case (1)
        grid_shape = [self%model%n_x, self%model%n_y, self%model%n_z]
        bmi_status = BMI_SUCCESS
        case default
        grid_shape = [-1, -1]
        bmi_status = BMI_FAILURE
    end select
    end function fm_grid_shape

    ! The total number of elements in a grid.
    function fm_grid_size(self, grid_id, grid_size) result (bmi_status)
    class (bmi_fastmech), intent(in) :: self
    integer, intent(in) :: grid_id
    integer, intent(out) :: grid_size
    integer :: bmi_status

    select case (grid_id)
    case (0)
        grid_size = self%model%n_y * self%model%n_x
        bmi_status = BMI_SUCCESS
    case (1)
        grid_size = self%model%n_y * self%model%n_x * self%model%n_z
        bmi_status = BMI_SUCCESS
        case default
        grid_size = -1
        bmi_status = BMI_FAILURE
    end select
    end function fm_grid_size

    !! The distance between nodes of a grid.
    !function fm_grid_spacing(self, grid_id, grid_spacing) result (bmi_status)
    !  class (bmi_fastmech), intent(in) :: self
    !  integer, intent(in) :: grid_id
    !  real, dimension(:), intent(out) :: grid_spacing
    !  integer :: bmi_status
    !
    !  select case(grid_id)
    !  case(0)
    !     grid_spacing = [self%model%dy, self%model%dx]
    !     bmi_status = BMI_SUCCESS
    !  case default
    !     grid_spacing = -1
    !     bmi_status = BMI_FAILURE
    !  end select
    !end function fm_grid_spacing
    !
    !! Coordinates of grid origin.
    !function fm_grid_origin(self, grid_id, grid_origin) result (bmi_status)
    !  class (bmi_fastmech), intent(in) :: self
    !  integer, intent(in) :: grid_id
    !  real, dimension(:), intent(out) :: grid_origin
    !  integer :: bmi_status
    !
    !  select case(grid_id)
    !  case(0)
    !     grid_origin = [0.0, 0.0]
    !     bmi_status = BMI_SUCCESS
    !  case default
    !     grid_origin = [-1.0]
    !     bmi_status = BMI_FAILURE
    !  end select
    !end function fm_grid_origin
    !
    ! X-coordinates of grid nodes.
    function fm_grid_x(self, grid_id, grid_x) result (bmi_status)
    class (bmi_fastmech), intent(in) :: self
    integer, intent(in) :: grid_id
    real, dimension(:), intent(out) :: grid_x
    integer :: bmi_status

    select case (grid_id)
    case (0)
        CALL Get_GRID_2D_COORD(self%model, 'x', grid_x)
        bmi_status = BMI_SUCCESS

    case(1)
        CALL Get_GRID_3D_COORD(self%model, 'x', grid_x)
        bmi_status = BMI_SUCCESS
        case default
        grid_x = [-1.0]
        bmi_status = BMI_FAILURE
    end select

    end function fm_grid_x

    ! Y-coordinates of grid nodes.
    function fm_grid_y(self, grid_id, grid_y) result (bmi_status)
    class (bmi_fastmech), intent(in) :: self
    integer, intent(in) :: grid_id
    real, dimension(:), intent(out) :: grid_y
    integer :: bmi_status

    select case (grid_id)
    case (0)
        CALL Get_GRID_2D_COORD(self%model,'y', grid_y)
        bmi_status = BMI_SUCCESS

    case(1)
        CALL Get_GRID_3D_COORD(self%model,'y', grid_y)
        bmi_status = BMI_SUCCESS
        case default
        grid_y = -[1.0]
        bmi_status = BMI_FAILURE
    end select
    end function fm_grid_y

    ! Z-coordinates of grid nodes.
    function fm_grid_z(self, grid_id, grid_z) result (bmi_status)
    class (bmi_fastmech), intent(in) :: self
    integer, intent(in) :: grid_id
    real, dimension(:), intent(out) :: grid_z
    integer :: bmi_status

    select case (grid_id)
    case (0)
        CALL Get_GRID_2D_COORD(self%model,'z', grid_z)
        bmi_status = BMI_SUCCESS

    case(1)
        CALL Get_GRID_3D_COORD(self%model,'z', grid_z)
        bmi_status = BMI_SUCCESS
        case default
        grid_z = [-1.0]
        bmi_status = BMI_FAILURE
    end select
    end function fm_grid_z

    !! Connectivity array of unstructured grid nodes.
    !function fm_grid_connectivity(self, grid_id, grid_conn) &
    !     result (bmi_status)
    !  class (bmi_fastmech), intent(in) :: self
    !  integer, intent(in) :: grid_id
    !  integer, dimension(:), intent(out) :: grid_conn
    !  integer :: bmi_status
    !
    !  select case(grid_id)
    !  case(1)
    !     grid_conn = [0]
    !     bmi_status = BMI_SUCCESS
    !  case default
    !     grid_conn = [-1]
    !     bmi_status = BMI_FAILURE
    !  end select
    !end function fm_grid_connectivity
    !
    !! Offsets of unstructured grid nodes.
    !function fm_grid_offset(self, grid_id, grid_offset) &
    !     result (bmi_status)
    !  class (bmi_fastmech), intent(in) :: self
    !  integer, intent(in) :: grid_id
    !  integer, dimension(:), intent(out) :: grid_offset
    !  integer :: bmi_status
    !
    !  select case(grid_id)
    !  case(1)
    !     grid_offset = [0]
    !     bmi_status = BMI_SUCCESS
    !  case default
    !     grid_offset = [-1]
    !     bmi_status = BMI_FAILURE
    !  end select
    !end function fm_grid_offset

    ! The data type of the variable, as a string.
    function fm_var_type(self, var_name, var_type) result (bmi_status)
    class (bmi_fastmech), intent(in) :: self
    character (len=*), intent(in) :: var_name
    character (len=*), intent(out) :: var_type
    integer :: bmi_status
    select case (var_name)
    case ('Elevation', 'roughness', 'vegroughness', &
        'Depth', 'Drag_Coefficient','ShearStressX', &
        'ShearStressY', 'WaterSurfaceElevation', &
        'VelocityX', 'VelocityY')
        var_type = "DOUBLE PRECISION"
        bmi_status = BMI_SUCCESS
    case ('IBC', 'FMIBC')
        var_type = "integer"
        bmi_status = BMI_SUCCESS
        case default
        var_type = "-"
        bmi_status = BMI_FAILURE
    end select

    end function fm_var_type

    ! The units of the given variable.
    function fm_var_units(self, var_name, var_units) result (bmi_status)
    class (bmi_fastmech), intent(in) :: self
    character (len=*), intent(in) :: var_name
    character (len=*), intent(out) :: var_units
    integer :: bmi_status

    select case (var_name)
    case ("Elevation", "Depth", "WaterSurfaceElevation")
        var_units = "cm"
        bmi_status = BMI_SUCCESS
    case ("roughness", "vegroughness", "Drag_Coefficient", &
        "IBC", "FMIBC")
        var_units = "unitless"
        bmi_status = BMI_SUCCESS
    case ("VelocityX", "VelocityY")
        var_units = "cm s-1"
        bmi_status = BMI_SUCCESS
    case ("ShearStressX", "ShearStressY")
        var_units = "D cm-2"
        case default
        var_units = "-"
        bmi_status = BMI_FAILURE
    end select
    end function fm_var_units

    ! Memory use per array element.
    function fm_var_itemsize(self, var_name, var_size) result (bmi_status)
    class (bmi_fastmech), intent(in) :: self
    character (len=*), intent(in) :: var_name
    integer, intent(out) :: var_size
    integer :: bmi_status
    select case (var_name)

    case('Elevation')
        !var_size = sizeof(self%model%t_rivvar%eta(1,1))
        var_size = sizeof(self%model%fm_elevation(1,1))
        bmi_status = BMI_SUCCESS
    case('roughness')
        var_size = sizeof(self%model%fm_roughness(1,1))
        bmi_status = BMI_SUCCESS
    case('vegroughness')
        var_size = sizeof(self%model%fm_vegroughness(1,1))
        bmi_status = BMI_SUCCESS
    case('Depth')
        var_size = sizeof(self%model%fm_depth(1,1))
        bmi_status = BMI_SUCCESS
    case('Drag_Coefficient')
        var_size = sizeof(self%model%fm_dragcoefficient(1,1))
        bmi_status = BMI_SUCCESS
    case('ShearStressX')
        var_size = sizeof(self%model%fm_shearstress_x(1,1))
        bmi_status = BMI_SUCCESS
    case('ShearStressY')
        var_size = sizeof(self%model%fm_shearstress_y(1,1))
        bmi_status = BMI_SUCCESS
    case('WaterSurfaceElevation')
        var_size = sizeof(self%model%fm_wse(1,1))
        bmi_status = BMI_SUCCESS
    case('VelocityX')
        var_size = sizeof(self%model%fm_velocity_x(1,1))
        bmi_status = BMI_SUCCESS
    case('VelocityY')
        var_size = sizeof(self%model%fm_velocity_y(1,1))
        bmi_status = BMI_SUCCESS
    case('IBC')
        var_size = sizeof(self%model%fm_ibc(1,1))
        bmi_status = BMI_SUCCESS
        case default
        var_size = -1
        bmi_status = BMI_FAILURE
    end select

    end function fm_var_itemsize

    ! The size of the given variable.
    function fm_var_nbytes(self, var_name, var_nbytes) result (bmi_status)
    class (bmi_fastmech), intent(in) :: self
    character (len=*), intent(in) :: var_name
    integer, intent(out) :: var_nbytes
    integer :: bmi_status
    integer :: s1, s2, s3, grid_id, grid_size, item_size

    s1 = self%get_var_grid(var_name, grid_id)
    s2 = self%get_grid_size(grid_id, grid_size)
    s3 = self%get_var_itemsize(var_name, item_size)

    if ((s1 == BMI_SUCCESS).and.(s2 == BMI_SUCCESS).and.(s3 == BMI_SUCCESS)) then
        var_nbytes = item_size * grid_size
        bmi_status = BMI_SUCCESS
    else
        var_nbytes = -1
        bmi_status = BMI_FAILURE
    end if
    end function fm_var_nbytes

    ! Get a copy of a integer variable's values, flattened.
    function fm_get_int(self, var_name, dest) result (bmi_status)
    class (bmi_fastmech), intent(in) :: self
    character (len=*), intent(in) :: var_name
    integer, intent(inout) :: dest(:)
    integer :: bmi_status

    select case(var_name)
    case("IBC")
        dest = reshape(self%model%fm_ibc, [self%model%n_x*self%model%n_y])
        bmi_status = BMI_SUCCESS
        case default
        dest = [-1]
        bmi_status = BMI_FAILURE
    end select
    end function fm_get_int

    ! Get a copy of a real variable's values, flattened.
    !function fm_get_float(self, var_name, dest) result (bmi_status)
    !  class (bmi_fastmech), intent(in) :: self
    !  character (len=*), intent(in) :: var_name
    !  real, intent(inout) :: dest(:)
    !  integer :: bmi_status
    !
    !  select case(var_name)
    !  !case("plate_surface__temperature")
    !  !   ! This would be safe, but subject to indexing errors.
    !  !   ! do j = 1, self%model%n_y
    !  !   !    do i = 1, self%model%n_x
    !  !   !       k = j + self%model%n_y*(i-1)
    !  !   !       dest(k) = self%model%temperature(j,i)
    !  !   !    end do
    !  !   ! end do
    !  !
    !  !   ! This is an equivalent, elementwise copy into `dest`.
    !  !   ! See https://stackoverflow.com/a/11800068/1563298
    !  !   dest = reshape(self%model%temperature, [self%model%n_x*self%model%n_y])
    !  !   bmi_status = BMI_SUCCESS
    !  !case("plate_surface__thermal_diffusivity")
    !  !   dest = [self%model%alpha]
    !  !   bmi_status = BMI_SUCCESS
    !  case default
    !     dest = [-1.0]
    !     bmi_status = BMI_FAILURE
    !  end select
    !end function fm_get_float

    ! Get a copy of a double variable's values, flattened.
    function fm_get_double(self, var_name, dest) result (bmi_status)
    class (bmi_fastmech), intent(in) :: self
    character (len=*), intent(in) :: var_name
    double precision, intent(inout) :: dest(:)
    integer :: bmi_status

    select case(var_name)
    case('Elevation')
        dest = reshape(self%model%t_rivvar%eta, [self%model%n_x*self%model%n_y])
        bmi_status = BMI_SUCCESS
    case('roughness')
        dest = reshape(self%model%fm_roughness, [self%model%n_x*self%model%n_y])
        bmi_status = BMI_SUCCESS
    case('vegroughness')
        dest = reshape(self%model%fm_vegroughness, [self%model%n_x*self%model%n_y])
        bmi_status = BMI_SUCCESS
    case('Depth')
        dest = reshape(self%model%fm_depth, [self%model%n_x*self%model%n_y])
        bmi_status = BMI_SUCCESS
    case('Drag_Coefficient')
        dest = reshape(self%model%fm_dragcoefficient, [self%model%n_x*self%model%n_y])
        bmi_status = BMI_SUCCESS
    case('ShearStressX')
        !call get_grid_2d_ssvec(self%model, 'ShearStressX', dest)
        dest = reshape(self%model%fm_shearstress_x, [self%model%n_x*self%model%n_y])
        bmi_status = BMI_SUCCESS
    case('ShearStressY')
        !call get_grid_2d_ssvec(self%model, 'ShearStressY', dest)
        dest = reshape(self%model%fm_shearstress_y, [self%model%n_x*self%model%n_y])
        bmi_status = BMI_SUCCESS
    case('WaterSurfaceElevation')
        dest = reshape(self%model%fm_wse, [self%model%n_x*self%model%n_y])
        bmi_status = BMI_SUCCESS
    case('VelocityX')
        !call get_grid_2d_velvec(self%model, 'VelocityX', dest)
        dest = reshape(self%model%fm_velocity_x, [self%model%n_x*self%model%n_y])
        bmi_status = BMI_SUCCESS
    case('VelocityY')
        !call get_grid_2d_velvec(self%model, 'VelocityY', dest)
        dest = reshape(self%model%fm_velocity_y, [self%model%n_x*self%model%n_y])
        bmi_status = BMI_SUCCESS
    case('IBC')
        dest = reshape(self%model%fm_ibc, [self%model%n_x*self%model%n_y])
        bmi_status = BMI_SUCCESS

        case default
        dest = [-1.d0]
        bmi_status = BMI_FAILURE
    end select
    end function fm_get_double

    ! Get a reference to an integer-valued variable, flattened.
    function fm_get_ref_int(self, var_name, dest) result (bmi_status)
    class (bmi_fastmech), intent(in) :: self
    character (len=*), intent(in) :: var_name
    integer, pointer, intent(inout) :: dest(:)
    integer :: bmi_status
    type (c_ptr) :: src
    integer :: n_elements

    select case(var_name)
    case('IBC')
        src = c_loc(self%model%fm_ibc(1,1))
        n_elements = self%model%n_y * self%model%n_x
        call c_f_pointer(src, dest, [n_elements])
        bmi_status = BMI_SUCCESS
        case default
        bmi_status = BMI_FAILURE
    end select
    end function fm_get_ref_int

    !! Get a reference to a real-valued variable, flattened.
    !function fm_get_ref_float(self, var_name, dest) result (bmi_status)
    !  class (bmi_fastmech), intent(in) :: self
    !  character (len=*), intent(in) :: var_name
    !  real, pointer, intent(inout) :: dest(:)
    !  integer :: bmi_status
    !  type (c_ptr) :: src
    !  integer :: n_elements
    !
    !  select case(var_name)
    !  case("plate_surface__temperature")
    !     src = c_loc(self%model%temperature(1,1))
    !     n_elements = self%model%n_y * self%model%n_x
    !     call c_f_pointer(src, dest, [n_elements])
    !     bmi_status = BMI_SUCCESS
    !  case default
    !     bmi_status = BMI_FAILURE
    !  end select
    !end function fm_get_ref_float

    ! Get a reference to an double-valued variable, flattened.
    function fm_get_ref_double(self, var_name, dest) result (bmi_status)
    class (bmi_fastmech), intent(in) :: self
    character (len=*), intent(in) :: var_name
    double precision, pointer, intent(inout) :: dest(:)
    integer :: bmi_status
    type (c_ptr) :: src
    integer :: n_elements

    select case(var_name)
    case('Elevation')
        src = c_loc(self%model%t_rivvar%eta(1,1))
        n_elements = self%model%n_y * self%model%n_x
        call c_f_pointer(src, dest, [n_elements])
        bmi_status = BMI_SUCCESS
    case('roughness')
        src = c_loc(self%model%fm_roughness(1,1))
        n_elements = self%model%n_y * self%model%n_x
        call c_f_pointer(src, dest, [n_elements])
        bmi_status = BMI_SUCCESS
    case('vegroughness')
        src = c_loc(self%model%fm_vegroughness(1,1))
        n_elements = self%model%n_y * self%model%n_x
        call c_f_pointer(src, dest, [n_elements])
        bmi_status = BMI_SUCCESS
    case('Depth')
        src = c_loc(self%model%fm_depth(1,1))
        n_elements = self%model%n_y * self%model%n_x
        call c_f_pointer(src, dest, [n_elements])
        bmi_status = BMI_SUCCESS
    case('Drag_Coefficient')
        src = c_loc(self%model%fm_dragcoefficient(1,1))
        n_elements = self%model%n_y * self%model%n_x
        call c_f_pointer(src, dest, [n_elements])
        bmi_status = BMI_SUCCESS
    case('ShearStressX')
        src = c_loc(self%model%fm_shearstress_x(1,1))
        n_elements = self%model%n_y * self%model%n_x
        call c_f_pointer(src, dest, [n_elements])
        bmi_status = BMI_SUCCESS
    case('ShearStressY')
        src = c_loc(self%model%fm_shearstress_y(1,1))
        n_elements = self%model%n_y * self%model%n_x
        call c_f_pointer(src, dest, [n_elements])
        bmi_status = BMI_SUCCESS
    case('WaterSurfaceElevation')
        src = c_loc(self%model%fm_wse(1,1))
        n_elements = self%model%n_y * self%model%n_x
        call c_f_pointer(src, dest, [n_elements])
        bmi_status = BMI_SUCCESS
    case('VelocityX')
        src = c_loc(self%model%fm_velocity_x(1,1))
        n_elements = self%model%n_y * self%model%n_x
        call c_f_pointer(src, dest, [n_elements])
        bmi_status = BMI_SUCCESS
    case('VelocityY')
        src = c_loc(self%model%fm_velocity_y(1,1))
        n_elements = self%model%n_y * self%model%n_x
        call c_f_pointer(src, dest, [n_elements])
        bmi_status = BMI_SUCCESS
        case default
        bmi_status = BMI_FAILURE
    end select
    end function fm_get_ref_double

    ! Get values of an integer variable at the given locations.
    function fm_get_at_indices_int(self, var_name, dest, indices) &
         result (bmi_status)
      class (bmi_fastmech), intent(in) :: self
      character (len=*), intent(in) :: var_name
      integer, intent(inout) :: dest(:)
      integer, intent(in) :: indices(:)
      integer :: bmi_status
      type (c_ptr) src
      integer, pointer :: src_flattened(:)
      integer :: i, n_elements
    
      select case(var_name)
      case('IBC')
         src = c_loc(self%model%fm_ibc(1,1))
         call c_f_pointer(src, src_flattened, [self%model%n_y * self%model%n_x])
         n_elements = size(indices)
         do i = 1, n_elements
            dest(i) = src_flattened(indices(i))
         end do
         bmi_status = BMI_SUCCESS
      case default
         bmi_status = BMI_FAILURE
      end select
    end function fm_get_at_indices_int
    
    !! Get values of a real variable at the given locations.
    !function fm_get_at_indices_float(self, var_name, dest, indices) &
    !     result (bmi_status)
    !  class (bmi_fastmech), intent(in) :: self
    !  character (len=*), intent(in) :: var_name
    !  real, intent(inout) :: dest(:)
    !  integer, intent(in) :: indices(:)
    !  integer :: bmi_status
    !  type (c_ptr) src
    !  real, pointer :: src_flattened(:)
    !  integer :: i, n_elements
    !
    !  select case(var_name)
    !  case("plate_surface__temperature")
    !     src = c_loc(self%model%temperature(1,1))
    !     call c_f_pointer(src, src_flattened, [self%model%n_y * self%model%n_x])
    !     n_elements = size(indices)
    !     do i = 1, n_elements
    !        dest(i) = src_flattened(indices(i))
    !     end do
    !     bmi_status = BMI_SUCCESS
    !  case default
    !     bmi_status = BMI_FAILURE
    !  end select
    !end function fm_get_at_indices_float
    !
    ! Get values of a double variable at the given locations.
    function fm_get_at_indices_double(self, var_name, dest, indices) &
         result (bmi_status)
      class (bmi_fastmech), intent(in) :: self
      character (len=*), intent(in) :: var_name
      double precision, intent(inout) :: dest(:)
      integer, intent(in) :: indices(:)
      integer :: bmi_status
      type (c_ptr) src
      double precision, pointer :: src_flattened(:)
      integer :: i, n_elements
    
      select case(var_name)
      case('Elevation')
         src = c_loc(self%model%t_rivvar%eta(1,1))
         call c_f_pointer(src, src_flattened, [self%model%n_y * self%model%n_x])
         n_elements = size(indices)
         do i = 1, n_elements
            dest(i) = src_flattened(indices(i))
         end do
         bmi_status = BMI_SUCCESS          
      case('roughness')
         src = c_loc(self%model%fm_roughness(1,1))
         call c_f_pointer(src, src_flattened, [self%model%n_y * self%model%n_x])
         n_elements = size(indices)
         do i = 1, n_elements
            dest(i) = src_flattened(indices(i))
         end do
         bmi_status = BMI_SUCCESS          
      case('vegroughness')
         src = c_loc(self%model%fm_vegroughness(1,1))
         call c_f_pointer(src, src_flattened, [self%model%n_y * self%model%n_x])
         n_elements = size(indices)
         do i = 1, n_elements
            dest(i) = src_flattened(indices(i))
         end do
         bmi_status = BMI_SUCCESS          
      case('Depth')
         src = c_loc(self%model%fm_depth(1,1))
         call c_f_pointer(src, src_flattened, [self%model%n_y * self%model%n_x])
         n_elements = size(indices)
         do i = 1, n_elements
            dest(i) = src_flattened(indices(i))
         end do
         bmi_status = BMI_SUCCESS          
      case('Drag_Coefficient')
         src = c_loc(self%model%fm_dragcoefficient(1,1))
         call c_f_pointer(src, src_flattened, [self%model%n_y * self%model%n_x])
         n_elements = size(indices)
         do i = 1, n_elements
            dest(i) = src_flattened(indices(i))
         end do
         bmi_status = BMI_SUCCESS          
      case('ShearStressX')
         src = c_loc(self%model%fm_shearstress_x(1,1))
         call c_f_pointer(src, src_flattened, [self%model%n_y * self%model%n_x])
         n_elements = size(indices)
         do i = 1, n_elements
            dest(i) = src_flattened(indices(i))
         end do
         bmi_status = BMI_SUCCESS          
      case('ShearStressY')
         src = c_loc(self%model%fm_shearstress_y(1,1))
         call c_f_pointer(src, src_flattened, [self%model%n_y * self%model%n_x])
         n_elements = size(indices)
         do i = 1, n_elements
            dest(i) = src_flattened(indices(i))
         end do
         bmi_status = BMI_SUCCESS          
      case('WaterSurfaceElevation')
         src = c_loc(self%model%fm_wse(1,1))
         call c_f_pointer(src, src_flattened, [self%model%n_y * self%model%n_x])
         n_elements = size(indices)
         do i = 1, n_elements
            dest(i) = src_flattened(indices(i))
         end do
         bmi_status = BMI_SUCCESS          
      case('VelocityX')
         src = c_loc(self%model%fm_velocity_x(1,1))
         call c_f_pointer(src, src_flattened, [self%model%n_y * self%model%n_x])
         n_elements = size(indices)
         do i = 1, n_elements
            dest(i) = src_flattened(indices(i))
         end do
         bmi_status = BMI_SUCCESS          
      case('VelocityY')
         src = c_loc(self%model%fm_velocity_y(1,1))
         call c_f_pointer(src, src_flattened, [self%model%n_y * self%model%n_x])
         n_elements = size(indices)
         do i = 1, n_elements
            dest(i) = src_flattened(indices(i))
         end do
         bmi_status = BMI_SUCCESS          
         
      case default
         bmi_status = BMI_FAILURE
      end select
    end function fm_get_at_indices_double
    
    ! Set new integer values.
    function fm_set_int(self, var_name, src) result (bmi_status)
      class (bmi_fastmech), intent(inout) :: self
      character (len=*), intent(in) :: var_name
      integer, intent(in) :: src(:)
      integer :: bmi_status
    
      select case(var_name)
      case("IBC")
         self%model%fm_ibc = reshape(src, [self%model%n_y, self%model%n_x])
         bmi_status = BMI_SUCCESS
      case default
         bmi_status = BMI_FAILURE
      end select
    end function fm_set_int
    
    !! Set new real values.
    !function fm_set_float(self, var_name, src) result (bmi_status)
    !  class (bmi_fastmech), intent(inout) :: self
    !  character (len=*), intent(in) :: var_name
    !  real, intent(in) :: src(:)
    !  integer :: bmi_status
    !
    !  select case(var_name)
    !  case("plate_surface__temperature")
    !     self%model%temperature = reshape(src, [self%model%n_y, self%model%n_x])
    !     bmi_status = BMI_SUCCESS
    !  case("plate_surface__thermal_diffusivity")
    !     self%model%alpha = src(1)
    !     bmi_status = BMI_SUCCESS
    !  case default
    !     bmi_status = BMI_FAILURE
    !  end select
    !end function fm_set_float
    !
    ! Set new double values.
    function fm_set_double(self, var_name, src) result (bmi_status)
      class (bmi_fastmech), intent(inout) :: self
      character (len=*), intent(in) :: var_name
      double precision, intent(in) :: src(:)
      integer :: bmi_status
    
      select case(var_name)
      case('Elevation')
          self%model%t_rivvar%eta = reshape(src, [self%model%n_y, self%model%n_x])
          bmi_status = BMI_SUCCESS
      case('roughness')
          self%model%fm_roughness = reshape(src, [self%model%n_y, self%model%n_x])
          bmi_status = BMI_SUCCESS
      case('vegroughness')
          self%model%fm_vegroughness = reshape(src, [self%model%n_y, self%model%n_x])
          bmi_status = BMI_SUCCESS
      case('Depth')
          self%model%fm_depth = reshape(src, [self%model%n_y, self%model%n_x])
          bmi_status = BMI_SUCCESS
      case('Drag_Coefficient')
          self%model%fm_dragcoefficient = reshape(src, [self%model%n_y, self%model%n_x])
          bmi_status = BMI_SUCCESS
      case('ShearStressX')
          self%model%fm_shearstress_x = reshape(src, [self%model%n_y, self%model%n_x])
          bmi_status = BMI_SUCCESS
      case('ShearStressY')
          self%model%fm_shearstress_y = reshape(src, [self%model%n_y, self%model%n_x])
          bmi_status = BMI_SUCCESS
      case('WaterSurfaceElevation')
          self%model%fm_wse = reshape(src, [self%model%n_y, self%model%n_x])
          bmi_status = BMI_SUCCESS
      case('VelocityX')
          self%model%fm_velocity_x = reshape(src, [self%model%n_y, self%model%n_x])
          bmi_status = BMI_SUCCESS
      case('VelocityY')
          self%model%fm_velocity_y = reshape(src, [self%model%n_y, self%model%n_x])
          bmi_status = BMI_SUCCESS

      case default
         bmi_status = BMI_FAILURE
      end select
    end function fm_set_double
    
    ! Set integer values at particular locations.
    function fm_set_at_indices_int(self, var_name, indices, src) &
         result (bmi_status)
      class (bmi_fastmech), intent(inout) :: self
      character (len=*), intent(in) :: var_name
      integer, intent(in) :: indices(:)
      integer, intent(in) :: src(:)
      integer :: bmi_status
      type (c_ptr) dest
      integer, pointer :: dest_flattened(:)
      integer :: i
    
      select case(var_name)
      case('IBC')
         dest = c_loc(self%model%fm_ibc(1,1))
         call c_f_pointer(dest, dest_flattened, [self%model%n_y * self%model%n_x])
         do i = 1, size(indices)
            dest_flattened(indices(i)) = src(i)
         end do
         bmi_status = BMI_SUCCESS
      case default
         bmi_status = BMI_FAILURE
      end select
    end function fm_set_at_indices_int
    
    !! Set real values at particular locations.
    !function fm_set_at_indices_float(self, var_name, indices, src) &
    !     result (bmi_status)
    !  class (bmi_fastmech), intent(inout) :: self
    !  character (len=*), intent(in) :: var_name
    !  integer, intent(in) :: indices(:)
    !  real, intent(in) :: src(:)
    !  integer :: bmi_status
    !  type (c_ptr) dest
    !  real, pointer :: dest_flattened(:)
    !  integer :: i
    !
    !  select case(var_name)
    !  case("plate_surface__temperature")
    !     dest = c_loc(self%model%temperature(1,1))
    !     call c_f_pointer(dest, dest_flattened, [self%model%n_y * self%model%n_x])
    !     do i = 1, size(indices)
    !        dest_flattened(indices(i)) = src(i)
    !     end do
    !     bmi_status = BMI_SUCCESS
    !  case default
    !     bmi_status = BMI_FAILURE
    !  end select
    !end function fm_set_at_indices_float
    !
    ! Set double values at particular locations.
    function fm_set_at_indices_double(self, var_name, indices, src) &
         result (bmi_status)
      class (bmi_fastmech), intent(inout) :: self
      character (len=*), intent(in) :: var_name
      integer, intent(in) :: indices(:)
      double precision, intent(in) :: src(:)
      integer :: bmi_status
      type (c_ptr) dest
      double precision, pointer :: dest_flattened(:)
      integer :: i
    
      select case(var_name)
      case('IBC')
         dest = c_loc(self%model%fm_ibc(1,1))
         call c_f_pointer(dest, dest_flattened, [self%model%n_y * self%model%n_x])
         do i = 1, size(indices)
            dest_flattened(indices(i)) = src(i)
         end do
         bmi_status = BMI_SUCCESS
      case('roughness')
         dest = c_loc(self%model%fm_roughness(1,1))
         call c_f_pointer(dest, dest_flattened, [self%model%n_y * self%model%n_x])
         do i = 1, size(indices)
            dest_flattened(indices(i)) = src(i)
         end do
         bmi_status = BMI_SUCCESS
      case('vegroughness')
         dest = c_loc(self%model%fm_vegroughness(1,1))
         call c_f_pointer(dest, dest_flattened, [self%model%n_y * self%model%n_x])
         do i = 1, size(indices)
            dest_flattened(indices(i)) = src(i)
         end do
         bmi_status = BMI_SUCCESS
      case('Depth')
         dest = c_loc(self%model%fm_depth(1,1))
         call c_f_pointer(dest, dest_flattened, [self%model%n_y * self%model%n_x])
         do i = 1, size(indices)
            dest_flattened(indices(i)) = src(i)
         end do
         bmi_status = BMI_SUCCESS
      case('Drag_Coefficient')
         dest = c_loc(self%model%fm_dragcoefficient(1,1))
         call c_f_pointer(dest, dest_flattened, [self%model%n_y * self%model%n_x])
         do i = 1, size(indices)
            dest_flattened(indices(i)) = src(i)
         end do
         bmi_status = BMI_SUCCESS
      case('ShearStressX')
         dest = c_loc(self%model%fm_shearstress_x(1,1))
         call c_f_pointer(dest, dest_flattened, [self%model%n_y * self%model%n_x])
         do i = 1, size(indices)
            dest_flattened(indices(i)) = src(i)
         end do
         bmi_status = BMI_SUCCESS
      case('ShearStressY')
         dest = c_loc(self%model%fm_shearstress_y(1,1))
         call c_f_pointer(dest, dest_flattened, [self%model%n_y * self%model%n_x])
         do i = 1, size(indices)
            dest_flattened(indices(i)) = src(i)
         end do
         bmi_status = BMI_SUCCESS
      case('WaterSurfaceElevation')
         dest = c_loc(self%model%fm_wse(1,1))
         call c_f_pointer(dest, dest_flattened, [self%model%n_y * self%model%n_x])
         do i = 1, size(indices)
            dest_flattened(indices(i)) = src(i)
         end do
         bmi_status = BMI_SUCCESS
      case('VelocityX')
         dest = c_loc(self%model%fm_velocity_x(1,1))
         call c_f_pointer(dest, dest_flattened, [self%model%n_y * self%model%n_x])
         do i = 1, size(indices)
            dest_flattened(indices(i)) = src(i)
         end do
         bmi_status = BMI_SUCCESS
      case('VelocityY')
         dest = c_loc(self%model%fm_velocity_y(1,1))
         call c_f_pointer(dest, dest_flattened, [self%model%n_y * self%model%n_x])
         do i = 1, size(indices)
            dest_flattened(indices(i)) = src(i)
         end do
         bmi_status = BMI_SUCCESS
          
      case default
         bmi_status = BMI_FAILURE
      end select
    end function fm_set_at_indices_double
    
    ! A non-BMI procedure for model introspection.
    subroutine print_model_info(self)
      class (bmi_fastmech), intent(in) :: self
    
      call print_info(self%model)
    end subroutine print_model_info

    end module bmifastmech
