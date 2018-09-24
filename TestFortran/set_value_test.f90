! Test the set_value and set_value_at_indices functions.
program set_value_test

  use bmifastmech
  use testing_helpers
  implicit none

  type (bmi_fastmech) :: m
  integer :: s, i, j, grid_id
  character (len=BMI_MAXVARNAMESTR), pointer :: names(:)
  integer :: dims(2), locations(11)
  real :: values(11)
  real, pointer :: z(:), y(:)
  character(len=30) :: rowfmt

  write (*,"(a)",advance="no") "Initializing..."
  s = m%initialize(".\..\TestBMI\Test1.cgn")
  write (*,*) "Done."

  s = m%get_output_var_names(names)
  write (*,"(a, a)") "Output variables: ", names

  s = m%get_var_grid(names(1), grid_id)
  s = m%get_grid_shape(grid_id, dims)
  write(rowfmt,'(a,i4,a)') '(', dims(2), '(1x,f6.1))'

  write (*, "(a)") "Initial values:"
  !s = m%get_value("Elevation", z)
  s = m%get_value("WaterSurfaceElevation", z)
  call print_array(z, dims)

  write (*,"(a)",advance="no") "Setting new values..."
  z = 10.2d0*100.d0
  s = m%set_value("WaterSurfaceElevation", z)
  write (*,*) "Done."
  write (*, "(a)") "New values:"
  s = m%get_value("WaterSurfaceElevation", y)
  call print_array(y, dims)

  write (*, "(a)") "Set values at three locations:"
  locations = [2201, 2202, 2203, 2204, 2205, 2206, 2207, 2208, 2209, 2210, 2211]
  values = [10.1, 10.1, 10.1, 10.1, 10.1, 10.1, 10.1, 10.1, 10.1, 10.1, 10.1]*100.d0
  write (*,*) "Locations: ", locations
  write (*,*) "Values: ", values
  s = m%set_value_at_indices("WaterSurfaceElevation", locations, values)
  write (*, "(a)") "New values:"
  s = m%get_value("WaterSurfaceElevation", y)
  call print_array(y, dims)

  s = m%update()
  
  write (*,"(a)", advance="no") "Finalizing..."
  s = m%finalize()
  write (*,*) "Done"

end program set_value_test
