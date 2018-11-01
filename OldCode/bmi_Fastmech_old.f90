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
    procedure :: get_grid_x => fm_grid_x
    procedure :: get_grid_y => fm_grid_y
    procedure :: get_grid_z => fm_grid_z
    procedure :: get_var_type => fm_var_type
    procedure :: get_var_units => fm_var_units
    procedure :: get_var_itemsize => fm_var_itemsize
    procedure :: get_var_nbytes => fm_var_nbytes
    procedure :: get_value => fm_get
    procedure :: get_value_ref => fm_get_ref
    procedure :: get_value_at_indices => fm_get_at_indices
    procedure :: set_value => fm_set
    procedure :: set_value_at_indices => fm_set_at_indices
    end type bmi_fastmech

    private :: fm_component_name, fm_input_var_names, fm_output_var_names
    private :: fm_initialize, fm_finalize
    private :: fm_start_time, fm_end_time, fm_current_time
    private :: fm_time_step, fm_time_units
    private :: fm_update, fm_update_frac, fm_update_until
    private :: fm_var_grid
    private :: fm_grid_type, fm_grid_rank, fm_grid_shape
    private :: fm_grid_size
    private :: fm_grid_x, fm_grid_y, fm_grid_z
    private :: fm_var_type, fm_var_units, fm_var_itemsize, fm_var_nbytes
    private :: fm_get, fm_get_ref, fm_get_at_indices
    private :: fm_set, fm_set_at_indices

    character (len=BMI_MAXCOMPNAMESTR), target :: &
        component_name = "Fastmech"

    ! Exchange items
    integer, parameter :: input_item_count = 3
    integer, parameter :: output_item_count = 10
    character (len=BMI_MAXVARNAMESTR), target, &
        dimension (input_item_count) :: &
        input_items = (/ &
        'Elevation   ', &
        'roughness   ', &
        'vegroughness'/)
    character (len=BMI_MAXVARNAMESTR), target, &
        dimension (output_item_count) :: &
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
    class (bmi_fastmech), intent (in) :: self
    character (len=BMI_MAXCOMPNAMESTR), pointer, intent (out) :: name
    integer :: bmi_status

    name => component_name
    bmi_status = BMI_SUCCESS
    end function fm_component_name

    ! List input variables.
    function fm_input_var_names(self, names) result (bmi_status)
    class (bmi_fastmech), intent (in) :: self
    character (len=BMI_MAXVARNAMESTR), pointer, intent (out) :: names(:)
    integer :: bmi_status

    names => input_items
    bmi_status = BMI_SUCCESS
    end function fm_input_var_names

    ! List output variables.
    function fm_output_var_names(self, names) result (bmi_status)
    class (bmi_fastmech), intent (in) :: self
    character (len=BMI_MAXVARNAMESTR), pointer, intent (out) :: names(:)
    integer :: bmi_status

    names => output_items
    bmi_status = BMI_SUCCESS
    end function fm_output_var_names

    ! BMI initializer.
    function fm_initialize(self, config_file) result (bmi_status)
    class (bmi_fastmech), intent (out) :: self
    character (len=*), intent (in) :: config_file
    integer :: bmi_status

    if (len (config_file) > 0) then
        call initialize_from_file(self%model, config_file)
    else
        call initialize_from_defaults(self%model)
    end if
    bmi_status = BMI_SUCCESS
    end function fm_initialize

    ! BMI finalizer.
    function fm_finalize(self) result (bmi_status)
    class (bmi_fastmech), intent (inout) :: self
    integer :: bmi_status

    call cleanup(self%model)
    bmi_status = BMI_SUCCESS
    end function fm_finalize
    !
    ! Model start time.
    function fm_start_time(self, time) result (bmi_status)
    class (bmi_fastmech), intent (in) :: self
    real, intent (out) :: time
    integer :: bmi_status
    if(self%model%t_calccond%soltype == 0) then
        time = 0
    else
        time = self%model%t_rivvartime%vardischstarttime
    endif
    bmi_status = BMI_SUCCESS
    end function fm_start_time

    ! Model end time.
    function fm_end_time(self, time) result (bmi_status)
    class (bmi_fastmech), intent (in) :: self
    real, intent (out) :: time
    integer :: bmi_status

    time = self%model%t_end
    bmi_status = BMI_SUCCESS
    end function fm_end_time

    ! Model current time.
    function fm_current_time(self, time) result (bmi_status)
    class (bmi_fastmech), intent (in) :: self
    real, intent (out) :: time
    integer :: bmi_status

    time = self%model%t
    bmi_status = BMI_SUCCESS
    end function fm_current_time

    ! Model time step.
    function fm_time_step(self, time_step) result (bmi_status)
    class (bmi_fastmech), intent (in) :: self
    real, intent (out) :: time_step
    integer :: bmi_status

    time_step = self%model%dt
    bmi_status = BMI_SUCCESS
    end function fm_time_step

    ! Model time units.
    function fm_time_units(self, time_units) result (bmi_status)
    class (bmi_fastmech), intent (in) :: self
    character (len=BMI_MAXUNITSSTR), intent (out) :: time_units
    integer :: bmi_status

    time_units = "s"
    bmi_status = BMI_SUCCESS
    end function fm_time_units

    ! Advance model by one time step.
    function fm_update(self) result (bmi_status)
    class (bmi_fastmech), intent (inout) :: self
    integer :: bmi_status

    call advance_in_time(self%model)
    bmi_status = BMI_SUCCESS
    end function fm_update

    ! Advance the model by a fraction of a time step.
    function fm_update_frac(self, time_frac) result (bmi_status)
    class (bmi_fastmech), intent (inout) :: self
    real, intent (in) :: time_frac
    integer :: bmi_status
    real :: time_step

    if (time_frac > 0.0) then
        time_step = self%model%dt
        self%model%dt = time_step*time_frac
        call advance_in_time(self%model)
        self%model%dt = time_step
    end if
    bmi_status = BMI_SUCCESS
    end function fm_update_frac

    ! Advance the model until the given time.
    function fm_update_until(self, time) result (bmi_status)
    class (bmi_fastmech), intent (inout) :: self
    real, intent (in) :: time
    integer :: bmi_status
    real :: n_steps_real
    integer :: n_steps, i, s

    if (time > self%model%t) then
        n_steps_real = (time - self%model%t) / self%model%dt
        n_steps = floor(n_steps_real)
        do i = 1, n_steps
            s = self%update()
        end do
        s = self%update_frac(n_steps_real - real(n_steps))
    end if
    bmi_status = BMI_SUCCESS
    end function fm_update_until

    ! Get the grid id for a particular variable.
    function fm_var_grid(self, var_name, grid_id) result (bmi_status)
    class (bmi_fastmech), intent (in) :: self
    character (len=BMI_MAXVARNAMESTR), intent (in) :: var_name
    integer, intent (out) :: grid_id
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
    class (bmi_fastmech), intent (in) :: self
    integer, intent (in) :: grid_id
    character (len=*), intent (out) :: grid_type
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
    function fm_grid_rank(self, grid_id, rank) result (bmi_status)
    class (bmi_fastmech), intent (in) :: self
    integer, intent (in) :: grid_id
    integer, intent (out) :: rank
    integer :: bmi_status

    select case (grid_id)
    case (0)
        rank = 2
        bmi_status = BMI_SUCCESS
    case (1)
        rank = 3
        bmi_status = BMI_SUCCESS
        case default
        rank = -1
        bmi_status = BMI_FAILURE
    end select
    end function fm_grid_rank

    ! The dimensions of a grid.
    function fm_grid_shape(self, grid_id, shape) result (bmi_status)
    class (bmi_fastmech), intent (in) :: self
    integer, intent (in) :: grid_id
    integer, dimension(:), intent (out) :: shape
    integer :: bmi_status

    select case (grid_id)
    case (0)
        shape = [self%model%n_x, self%model%n_y]
        bmi_status = BMI_SUCCESS
    case (1)
        shape = [self%model%n_x, self%model%n_y, self%model%n_z]
        bmi_status = BMI_SUCCESS
        case default
        shape = [-1, -1]
        bmi_status = BMI_FAILURE
    end select
    end function fm_grid_shape

    ! The total number of elements in a grid.
    function fm_grid_size(self, grid_id, size) result (bmi_status)
    class (bmi_fastmech), intent (in) :: self
    integer, intent (in) :: grid_id
    integer, intent (out) :: size
    integer :: bmi_status

    select case (grid_id)
    case (0)
        size = self%model%n_y * self%model%n_x
        bmi_status = BMI_SUCCESS
    case (1)
        size = self%model%n_y * self%model%n_x * self%model%n_z
        bmi_status = BMI_SUCCESS
        case default
        size = -1
        bmi_status = BMI_FAILURE
    end select
    end function fm_grid_size

    ! The x-coordinate nodes of a grid.
    function fm_grid_x(self, grid_id, grid_x) result (bmi_status)
    class (bmi_fastmech), intent (in) :: self
    integer, intent (in) :: grid_id
    real, pointer, intent (inout) :: grid_x(:)
    integer :: bmi_status, ier
    integer :: size

    select case (grid_id)
    case (0)
        bmi_status = self%get_grid_size(grid_id, size)
        CALL Get_GRID_2D_COORD(self%model, 'x', grid_x)
        bmi_status = BMI_SUCCESS

    case(1)
        bmi_status = self%get_grid_size(grid_id, size)
        CALL Get_GRID_3D_COORD(self%model, 'x', grid_x)
        bmi_status = BMI_SUCCESS
        case default
        grid_x = -1.0
        bmi_status = BMI_FAILURE
    end select
    end function fm_grid_x

! The y-coordinate nodes of a grid.
    function fm_grid_y(self, grid_id, grid_y) result (bmi_status)
    class (bmi_fastmech), intent (in) :: self
    integer, intent (in) :: grid_id
    real, pointer, intent (inout) :: grid_y(:)
    integer :: bmi_status, ier
    integer :: size

    select case (grid_id)
    case (0)
        bmi_status = self%get_grid_size(grid_id, size)
        CALL Get_GRID_2D_COORD(self%model,'y', grid_y)
        bmi_status = BMI_SUCCESS

    case(1)
        bmi_status = self%get_grid_size(grid_id, size)
        CALL Get_GRID_3D_COORD(self%model,'y', grid_y)
        bmi_status = BMI_SUCCESS
        case default
        grid_y = -1.0
        bmi_status = BMI_FAILURE
    end select
    end function fm_grid_y
    
! The y-coordinate nodes of a grid.
    function fm_grid_z(self, grid_id, grid_z) result (bmi_status)
    class (bmi_fastmech), intent (in) :: self
    integer, intent (in) :: grid_id
    real, pointer, intent (inout) :: grid_z(:)
    integer :: bmi_status, ier
    integer :: size

    select case (grid_id)
    case (0)
        bmi_status = self%get_grid_size(grid_id, size)
        !ALLOCATE(grid_z(size), stat=ier)
        CALL Get_GRID_2D_COORD(self%model,'z', grid_z)
        bmi_status = BMI_SUCCESS

    case(1)
        bmi_status = self%get_grid_size(grid_id, size)
        CALL Get_GRID_3D_COORD(self%model,'z', grid_z)
        bmi_status = BMI_SUCCESS
        case default
        grid_z = -1.0
        bmi_status = BMI_FAILURE
    end select
    end function fm_grid_z
    
	! The data type of the variable, as a string.
    function fm_var_type(self, var_name, type) result (bmi_status)
    class (bmi_fastmech), intent (in) :: self
    character (len=BMI_MAXVARNAMESTR), intent (in) :: var_name
    character (len=BMI_MAXUNITSSTR), intent (out) :: type
    integer :: bmi_status

    select case (var_name)
    case ('Elevation', 'roughness', 'vegroughness', &
        'Depth', 'Drag_Coefficient','ShearStressX', &
        'ShearStressY', 'WaterSurfaceElevation', &
        'VelocityX', 'VelocityY')
        type = "float64"
        bmi_status = BMI_SUCCESS
    case ('IBC', 'FMIBC')
        type = "int32"
        bmi_status = BMI_SUCCESS
    case default
        type = "-"
        bmi_status = BMI_FAILURE
    end select
    end function fm_var_type

    ! The units of the given variable.
    function fm_var_units(self, var_name, units) result (bmi_status)
      class (bmi_fastmech), intent (in) :: self
      character (len=BMI_MAXVARNAMESTR), intent (in) :: var_name
      character (len=BMI_MAXUNITSSTR), intent (out) :: units
      integer :: bmi_status
    
      select case (var_name)
      case ("Elevation", "Depth", "WaterSurfaceElevation")
         units = "cm"
         bmi_status = BMI_SUCCESS
      case ("roughness", "vegroughness", "Drag_Coefficient", &
            "IBC", "FMIBC")
          units = "unitless"
          bmi_status = BMI_SUCCESS
      case ("VelocityX", "VelocityY")
          units = "cm s-1"
          bmi_status = BMI_SUCCESS
      case ("ShearStressX", "ShearStressY")
          units = "D cm-2"
      case default
         units = "-"
         bmi_status = BMI_FAILURE
      end select
    end function fm_var_units
    
    ! Memory use per array element.
    function fm_var_itemsize(self, var_name, size) result (bmi_status)
      class (bmi_fastmech), intent (in) :: self
      character (len=BMI_MAXVARNAMESTR), intent (in) :: var_name
      integer, intent (out) :: size
      integer :: bmi_status
    
      select case (var_name)
      case ('Elevation', 'roughness', 'vegroughness', &
        'Depth', 'Drag_Coefficient','ShearStressX', &
        'ShearStressY', 'WaterSurfaceElevation', &
        'VelocityX', 'VelocityY')
         size = BMI_DOUBLE
         bmi_status = BMI_SUCCESS
    case ('IBC', 'FMIBC')
        size = BMI_LONG
        bmi_status = BMI_SUCCESS
      case default
         size = -1
         bmi_status = BMI_FAILURE
      end select
    end function fm_var_itemsize
    
    ! The size of the given variable.
    function fm_var_nbytes(self, var_name, size) result (bmi_status)
      class (bmi_fastmech), intent (in) :: self
      character (len=BMI_MAXVARNAMESTR), intent (in) :: var_name
      integer, intent (out) :: size
      integer :: bmi_status
      integer :: s1, s2, s3, grid_id, grid_size, item_size
    
      s1 = self%get_var_grid(var_name, grid_id)
      s2 = self%get_grid_size(grid_id, grid_size)
      s3 = self%get_var_itemsize(var_name, item_size)
    
      if ((s1 == BMI_SUCCESS).and.(s2 == BMI_SUCCESS).and.(s3 == BMI_SUCCESS)) then
         size = item_size * grid_size
         bmi_status = BMI_SUCCESS
      else
         size = -1
         bmi_status = BMI_FAILURE
      end if
    end function fm_var_nbytes
    
    ! Get a copy of a variable's values, flattened.
    function fm_get(self, var_name, dest) result (bmi_status)
    class (bmi_fastmech), intent (in) :: self
    !type(rivvar), pointer :: rvo
    character (len=*), intent (in) :: var_name
    character (len=BMI_MAXVARNAMESTR) :: tmpname
    real, pointer, intent (inout) :: dest(:)
    integer :: bmi_status
    integer :: n_elements, s1, s2, grid_id
    !rvo => self%model%t_rivvar
    select case (var_name)
    case ('Elevation')
        tmpname = adjustl(var_name) !did it htis way because get_var_grid uses string the BMI_MAVARNAMESTR
        s1 = self%get_var_grid(tmpname, grid_id)
        s2 = self%get_grid_size(grid_id, n_elements)
        !n_elements = self%model%n_y * self%model%n_x
        allocate(dest(n_elements))
        dest = reshape(self%model%t_rivvar%eta, [n_elements])
        bmi_status = BMI_SUCCESS
    case ('roughness')
        tmpname = adjustl(var_name) !did it htis way because get_var_grid uses string the BMI_MAVARNAMESTR
        s1 = self%get_var_grid(tmpname, grid_id)
        s2 = self%get_grid_size(grid_id, n_elements)
        !n_elements = self%model%n_y * self%model%n_x
        allocate(dest(n_elements))
        dest = reshape(self%model%t_rivvar%cd, [n_elements])
        bmi_status = BMI_SUCCESS
    case ('vegroughness')
        tmpname = adjustl(var_name) !did it htis way because get_var_grid uses string the BMI_MAVARNAMESTR
        s1 = self%get_var_grid(tmpname, grid_id)
        s2 = self%get_grid_size(grid_id, n_elements)
        !n_elements = self%model%n_y * self%model%n_x
        allocate(dest(n_elements))
        dest = reshape(self%model%t_rivvar%cdv, [n_elements])
        bmi_status = BMI_SUCCESS
    case ('Depth')
        tmpname = adjustl(var_name) !did it htis way because get_var_grid uses string the BMI_MAVARNAMESTR
        s1 = self%get_var_grid(tmpname, grid_id)
        s2 = self%get_grid_size(grid_id, n_elements)
        !n_elements = self%model%n_y * self%model%n_x
        allocate(dest(n_elements))
        dest = reshape(self%model%t_rivvar%hl, [n_elements])
        bmi_status = BMI_SUCCESS
    case ('Drag_Coefficient')
        tmpname = adjustl(var_name) !did it htis way because get_var_grid uses string the BMI_MAVARNAMESTR
        s1 = self%get_var_grid(tmpname, grid_id)
        s2 = self%get_grid_size(grid_id, n_elements)
        !n_elements = self%model%n_y * self%model%n_x
        allocate(dest(n_elements))
        dest = reshape(self%model%t_rivvar%totcd, [n_elements])
        bmi_status = BMI_SUCCESS
    case ('ShearStressX')
        tmpname = adjustl(var_name) !did it htis way because get_var_grid uses string the BMI_MAVARNAMESTR
        s1 = self%get_var_grid(tmpname, grid_id)
        s2 = self%get_grid_size(grid_id, n_elements)
        !n_elements = self%model%n_y * self%model%n_x
        allocate(dest(n_elements))
        dest = reshape(self%model%t_rivvar%taus, [n_elements])
        bmi_status = BMI_SUCCESS
    case ('ShearStressY')
        tmpname = adjustl(var_name) !did it htis way because get_var_grid uses string the BMI_MAVARNAMESTR
        s1 = self%get_var_grid(tmpname, grid_id)
        s2 = self%get_grid_size(grid_id, n_elements)
        !n_elements = self%model%n_y * self%model%n_x
        allocate(dest(n_elements))
        dest = reshape(self%model%t_rivvar%taun, [n_elements])
        bmi_status = BMI_SUCCESS
    case ('WaterSurfaceElevation')
        tmpname = adjustl(var_name) !did it htis way because get_var_grid uses string the BMI_MAVARNAMESTR
        s1 = self%get_var_grid(tmpname, grid_id)
        s2 = self%get_grid_size(grid_id, n_elements)
        !n_elements = self%model%n_y * self%model%n_x
        allocate(dest(n_elements))
        dest = reshape(self%model%t_rivvar%e, [n_elements])
        bmi_status = BMI_SUCCESS
    case ('VelocityX')
        tmpname = adjustl(var_name) !did it htis way because get_var_grid uses string the BMI_MAVARNAMESTR
        s1 = self%get_var_grid(tmpname, grid_id)
        s2 = self%get_grid_size(grid_id, n_elements)
        !n_elements = self%model%n_y * self%model%n_x
        allocate(dest(n_elements))
        dest = reshape(self%model%t_rivvar%u, [n_elements])
        bmi_status = BMI_SUCCESS
    case ('VelocityY')
        tmpname = adjustl(var_name) !did it htis way because get_var_grid uses string the BMI_MAVARNAMESTR
        s1 = self%get_var_grid(tmpname, grid_id)
        s2 = self%get_grid_size(grid_id, n_elements)
        !n_elements = self%model%n_y * self%model%n_x
        allocate(dest(n_elements))
        dest = reshape(self%model%t_rivvar%v, [n_elements])
        bmi_status = BMI_SUCCESS
    case default
        n_elements = 1
        allocate(dest(n_elements))
        dest = -1.0
        bmi_status = BMI_FAILURE
    end select
    end function fm_get
    
    ! Get a reference to a variable's values, flattened.
    function fm_get_ref(self, var_name, dest) result (bmi_status)
    class (bmi_fastmech), intent (in) :: self
    character (len=*), intent (in) :: var_name
    character (len=BMI_MAXVARNAMESTR) :: tmpname
    real, pointer, intent (inout) :: dest(:)
    integer :: bmi_status
    type (c_ptr) :: src
    real, pointer :: src_flattened(:)
    integer :: n_elements, s1, s2, grid_id
    tmpname = adjustl(var_name) !did it htis way because get_var_grid uses string the BMI_MAVARNAMESTR
    s1 = self%get_var_grid(tmpname, grid_id)
    s2 = self%get_grid_size(grid_id, n_elements)

    select case (var_name)
    case ("Elevation")
        src = c_loc (self%model%t_rivvar%eta(1,1))
        call c_f_pointer(src, dest, [n_elements])
        bmi_status = BMI_SUCCESS
    case ("Depth")
        src = c_loc (self%model%t_rivvar%hl(1,1))
        call c_f_pointer(src, dest, [n_elements])
        bmi_status = BMI_SUCCESS
    case ("WaterSurfaceElevation")
        src = c_loc (self%model%t_rivvar%e(1,1))
        call c_f_pointer(src, dest, [n_elements])
        bmi_status = BMI_SUCCESS
    case ("DragCoefficient")
        src = c_loc (self%model%t_rivvar%totcd(1,1))
        call c_f_pointer(src, dest, [n_elements])
        bmi_status = BMI_SUCCESS
    case ("roughness")
        src = c_loc (self%model%t_rivvar%cd(1,1))
        call c_f_pointer(src, dest, [n_elements])
        bmi_status = BMI_SUCCESS
    case ("vegroughness")
        src = c_loc (self%model%t_rivvar%cdv(1,1))
        call c_f_pointer(src, dest, [n_elements])
        bmi_status = BMI_SUCCESS
    case ("VelocityX")
        src = c_loc (self%model%t_rivvar%u(1,1))
        call c_f_pointer(src, dest, [n_elements])
        bmi_status = BMI_SUCCESS
    case ("VelocityY")
        src = c_loc (self%model%t_rivvar%v(1,1))
        call c_f_pointer(src, dest, [n_elements])
        bmi_status = BMI_SUCCESS
    case ("ShearStressX")
        src = c_loc (self%model%t_rivvar%taus(1,1))
        call c_f_pointer(src, dest, [n_elements])
        bmi_status = BMI_SUCCESS
    case ("ShearStressY")
        src = c_loc (self%model%t_rivvar%taun(1,1))
        call c_f_pointer(src, dest, [n_elements])
        bmi_status = BMI_SUCCESS

    case default
        bmi_status = BMI_FAILURE
    end select
    end function fm_get_ref
    
    ! Get values of a variable at the given locations.
    function fm_get_at_indices(self, var_name, dest, indices) result (bmi_status)
    class (bmi_fastmech), intent (in) :: self
    character (len=*), intent (in) :: var_name
    character (len=BMI_MAXVARNAMESTR) :: tmpname
    real, pointer, intent (inout) :: dest(:)
    integer, intent (in) :: indices(:)
    integer :: bmi_status
    type (c_ptr) src
    real, pointer :: src_flattened(:)
    integer :: i, n_elements, s1, s2, grid_id
    tmpname = adjustl(var_name) !did it htis way because get_var_grid uses string with length BMI_MAVARNAMESTR
    s1 = self%get_var_grid(tmpname, grid_id)
    s2 = self%get_grid_size(grid_id, n_elements)

    select case (var_name)
    case ("Elevation")
        src = c_loc (self%model%t_rivvar%eta(1,1))
        bmi_status = BMI_SUCCESS
    case ("Depth")
        src = c_loc (self%model%t_rivvar%hl(1,1))
        bmi_status = BMI_SUCCESS
    case ("WaterSurfaceElevation")
        src = c_loc (self%model%t_rivvar%e(1,1))
        bmi_status = BMI_SUCCESS
    case ("DragCoefficient")
        src = c_loc (self%model%t_rivvar%totcd(1,1))
        bmi_status = BMI_SUCCESS
    case ("roughness")
        src = c_loc (self%model%t_rivvar%cd(1,1))
        bmi_status = BMI_SUCCESS
    case ("vegroughness")
        src = c_loc (self%model%t_rivvar%cdv(1,1))
        bmi_status = BMI_SUCCESS
    case ("VelocityX")
        src = c_loc (self%model%t_rivvar%u(1,1))
        bmi_status = BMI_SUCCESS
    case ("VelocityY")
        src = c_loc (self%model%t_rivvar%v(1,1))
        bmi_status = BMI_SUCCESS
    case ("ShearStressX")
        src = c_loc (self%model%t_rivvar%taus(1,1))
        bmi_status = BMI_SUCCESS
    case ("ShearStressY")
        src = c_loc (self%model%t_rivvar%taun(1,1))
        bmi_status = BMI_SUCCESS
    case default
        bmi_status = BMI_FAILURE
    end select
    if(bmi_status == BMI_SUCCESS) then
        call c_f_pointer(src, src_flattened, [n_elements])
        n_elements = size (indices)
        allocate(dest(n_elements))
        !select case (var_name)
        !case ('Elevation', 'Depth','ShearStressX', &
        !'ShearStressY', 'WaterSurfaceElevation', &
        !'VelocityX', 'VelocityY')
        !    ! convert cgs to mks units
        !    do i = 1, n_elements
        !        dest(i) = src_flattened(indices(i))/100.0D0
        !    end do
        !case ('ShearStressX', 'ShearStressY')
        !    ! convert cgs to mks units
        !    do i = 1, n_elements
        !        dest(i) = src_flattened(indices(i))/10.0D0
        !    end do
        !
        !case default
            do i = 1, n_elements
                dest(i) = src_flattened(indices(i))
            end do
        !end select
    
    endif

    end function fm_get_at_indices
    
    ! Set new values.
    function fm_set(self, var_name, src) result (bmi_status)
    class (bmi_fastmech), intent (inout) :: self
    character (len=*), intent (in) :: var_name
    character (len=BMI_MAXVARNAMESTR) :: tmpname
    real, intent (in) :: src(:)
    integer :: dims(2)
    integer :: bmi_status
    integer :: n_elements, s1, s2, grid_id
    integer :: i
    tmpname = adjustl(var_name) !did it htis way because get_var_grid uses string with length BMI_MAVARNAMESTR
    s1 = self%get_var_grid(tmpname, grid_id)
    s2 = self%get_grid_shape(grid_id, dims)

    select case (var_name)
    case ("Elevation")
        self%model%t_rivvar%eta = reshape(src, dims)
        bmi_status = BMI_SUCCESS
    case ("Depth")
        self%model%t_rivvar%hl = reshape(src, dims)
        bmi_status = BMI_SUCCESS
    case ("WaterSurfaceElevation")
        self%model%t_rivvar%e = reshape(src, dims)
        bmi_status = BMI_SUCCESS
    case ("DragCoefficient")
        self%model%t_rivvar%totcd = reshape(src, dims)
        bmi_status = BMI_SUCCESS
    case ("roughness")
        self%model%t_rivvar%cd = reshape(src, dims)
        bmi_status = BMI_SUCCESS
    case ("vegroughness")
        self%model%t_rivvar%cdv = reshape(src, dims)
        bmi_status = BMI_SUCCESS
    case ("VelocityX")
        self%model%t_rivvar%u = reshape(src, dims)
        bmi_status = BMI_SUCCESS
    case ("VelocityY")
        self%model%t_rivvar%v = reshape(src, dims)
        bmi_status = BMI_SUCCESS
    case ("ShearStressX")
        self%model%t_rivvar%taus = reshape(src, dims)
        bmi_status = BMI_SUCCESS
    case ("ShearStressy")
        self%model%t_rivvar%taun = reshape(src, dims)
        bmi_status = BMI_SUCCESS
    case default
        bmi_status = BMI_FAILURE
    end select
    end function fm_set

    ! Set new values at particular locations.
    function fm_set_at_indices(self, var_name, indices, src) result (bmi_status)
    class (bmi_fastmech), intent (inout) :: self
    character (len=*), intent (in) :: var_name
    character (len=BMI_MAXVARNAMESTR) :: tmpname
    integer, intent (in) :: indices(:)
    real, intent (in) :: src(:)
    integer :: bmi_status
    type (c_ptr) dest
    real, pointer :: dest_flattened(:)
    integer :: i, n_elements, s1, s2, grid_id
    tmpname = adjustl(var_name) !did it htis way because get_var_grid uses string with length BMI_MAVARNAMESTR
    s1 = self%get_var_grid(tmpname, grid_id)
    s2 = self%get_grid_size(grid_id, n_elements)

    select case (var_name)
    case ("Elevation")
        dest = c_loc (self%model%t_rivvar%eta(1,1))
        bmi_status = BMI_SUCCESS
    case ("Depth")
        dest = c_loc (self%model%t_rivvar%hl(1,1))
        bmi_status = BMI_SUCCESS
    case ("WaterSurfaceElevation")
        dest = c_loc (self%model%t_rivvar%e(1,1))
        bmi_status = BMI_SUCCESS
    case ("DragCoefficient")
        dest = c_loc (self%model%t_rivvar%totcd(1,1))
        bmi_status = BMI_SUCCESS
    case ("roughness")
        dest = c_loc (self%model%t_rivvar%cd(1,1))
        bmi_status = BMI_SUCCESS
    case ("vegroughness")
        dest = c_loc (self%model%t_rivvar%cdv(1,1))
        bmi_status = BMI_SUCCESS
    case ("VelocityX")
        dest = c_loc (self%model%t_rivvar%u(1,1))
        bmi_status = BMI_SUCCESS
    case ("VelocityY")
        dest = c_loc (self%model%t_rivvar%v(1,1))
        bmi_status = BMI_SUCCESS
    case ("ShearStressX")
        dest = c_loc (self%model%t_rivvar%taus(1,1))
        bmi_status = BMI_SUCCESS
    case ("ShearStressY")
        dest = c_loc (self%model%t_rivvar%taun(1,1))
        bmi_status = BMI_SUCCESS
    case default
        bmi_status = BMI_FAILURE
    end select
    if(bmi_status == BMI_SUCCESS) then
        call c_f_pointer(dest, dest_flattened, [n_elements])
        do i = 1, size (indices)
            dest_flattened(indices(i)) = src(i)
        end do
    endif
    end function fm_set_at_indices

    end module bmifastmech
