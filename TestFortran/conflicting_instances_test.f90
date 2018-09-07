! Do two instances of bmi_heat conflict?
program conflicting_instances_test

  use bmifastmech
  use testing_helpers
  implicit none

  type (bmi_fastmech) :: m1
  type (bmi_fastmech) :: m2
  character (len=BMI_MAXVARNAMESTR) :: &
       cfg_file1 = ".\..\TestBMI\Case2.cgn", cfg_file2 = ".\..\TestBMI\Case2_1.cgn"
  integer :: s
  integer :: grid_id1, grid_id2
  character (len=BMI_MAXVARNAMESTR), pointer :: names1(:), names2(:)
  integer :: dims1(2), dims2(2)
  real, pointer :: z1(:), z2(:)
  character(len=30) :: rowfmt1, rowfmt2

  write(*, "(a, a10, a10)") "Configuration files: ", cfg_file1, cfg_file2

  write (*,"(a)",advance="no") "Initializing..."
  s = m1%initialize(cfg_file1)
  s = m2%initialize(cfg_file2)
  write (*,*) "Done."

  s = m1%get_output_var_names(names1)
  s = m1%get_var_grid(names1(1), grid_id1)
  s = m1%get_grid_shape(grid_id1, dims1)
  write(rowfmt1,'(a,i4,a)') '(', dims1(2), '(1x,f6.1))'

  write (*, "(a)") "Initial values, model 1:"
  s = m1%get_value("WaterSurfaceElevation", z1)
  call print_array(z1, dims1)

  s = m2%get_output_var_names(names2)
  s = m2%get_var_grid(names2(1), grid_id2)
  s = m2%get_grid_shape(grid_id2, dims2)
  write(rowfmt2,'(a,i4,a)') '(', dims2(2), '(1x,f6.1))'

  write (*, "(a)") "Initial values, model 2:"
  s = m2%get_value("WaterSurfaceElevation", z2)
  call print_array(z2, dims2)

  write (*, "(a)") "Set new value in model 2; does it affect model 1?"
  s = m2%set_value_at_indices("WaterSurfaceElevation", & ! array ordered ji, changing downstream boundary by 0.1
      ![201, 402, 603, 804, 1005, 1206, 1407, 1608, 1809, 2010, 2211], &
      [2201, 2202, 2203, 2204, 2205, 2206, 2207, 2208, 2209, 2210, 2211], & 
      [10.1, 10.1, 10.1, 10.1, 10.1, 10.1, 10.1, 10.1, 10.1, 10.1, 10.1])
  write (*, "(a)") "New values, model 2:"
  s = m2%get_value("WaterSurfaceElevation", z2)
  call print_array(z2, dims2)
  write (*, "(a)") "New values, model 1:"
  s = m1%get_value("WaterSurfaceElevation", z1)
  call print_array(z1, dims1)

  write (*, "(a)") "Update both models by one time step..."
  s = m1%update()
  s = m2%update()
  write (*, "(a)") "Updated values, model 1:"
  s = m1%get_value("WaterSurfaceElevation", z1)
  call print_array(z1, dims1)
  write (*, "(a)") "Updated values, model 2:"
  s = m2%get_value("WaterSurfaceElevation", z2)
  call print_array(z2, dims2)

  write (*,"(a)", advance="no") "Finalizing..."
  s = m1%finalize()
  s = m2%finalize()
  write (*,*) "Done"

end program conflicting_instances_test
