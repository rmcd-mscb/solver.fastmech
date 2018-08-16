!  GetValueTest.f90 
!
!  FUNCTIONS:
!  GetValueTest - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: GetValueTest
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program GetValueTest
    use bmifastmech
    use testing_helpers

    implicit none

    ! Variables
type (bmi_fastmech) :: m
  integer :: s, i, j, grid_id, countji
  character (len=BMI_MAXVARNAMESTR), pointer :: names(:)
  integer :: dims(2), locations(11)
  real, pointer :: z(:), twse(:), telev(:), tdepth(:)
  real :: time
    ! Body of GetValueTest


  write (6,"(a)",advance="no") "Initializing..."
  
  s = m%initialize(".\..\TestBMI\Test1.cgn")
  write (*,*) "Done."

  s = m%get_output_var_names(names)
  write (*,"(a, a)") "Output variables: ", names

  s = m%get_var_grid(names(3), grid_id)
  s = m%get_grid_shape(grid_id, dims)
  write(*,'(a,2i4)') 'Grid shape (ny,nx): ', dims
  
  write (*, "(a)") "Initial values:"
  s = m%get_value("Elevation", z)
  call print_array(z, dims)
  write (*, "(a, i5)") "Shape of returned values:", shape(z)
  
  write (*,"(a)") "Running (using get_value)..."
  do j = 1, 4
     s = m%update()
     s = m%get_value("Elevation", z)
     s = m%get_current_time(time)
     write (*,"(a, f6.1)") "Current time:", time
     call print_array(z, dims)
  end do
  write (*,"(a)") "Done."
  
  write (*, "(a)") "Values at three locations:"
  i = 1
  do j = 1, 11
      countji = ((j-1)*dims(2))+i
      locations(j) = countji
  enddo
  !locations = [1,2,3,4,5,6,7,8,9,10,11] !upstream boundary
  write (*,*) "Locations: ", locations
  s = m%get_value_at_indices("WaterSurfaceElevation", twse, locations)
  s = m%get_value_at_indices("Elevation", telev, locations)
  s = m%get_value_at_indices("Depth", tdepth, locations)
  do i = 1,size(locations)
    write (*,*) "Values: ", twse(i), telev(i), tdepth(i), twse(i)-telev(i)
  enddo 
  write (*,"(a)") "Running (using get_value_ref)..."
  s = m%get_value_ref("Depth", tdepth)
  do j = 1, 4
     s = m%update()
     s = m%get_current_time(time)
     write (*,"(a, f6.1)") "Current time:", time
     call print_array(tdepth, dims)
  end do
  write (*,"(a)") "Done."
  
  write (*, "(a)") "Values at three locations:"
  !locations = [21, 41, 62]
  write (*,*) "Locations: ", locations
  s = m%get_value_at_indices("WaterSurfaceElevation", twse, locations)
  s = m%get_value_at_indices("Elevation", telev, locations)
  s = m%get_value_at_indices("Depth", tdepth, locations)
  do i = 1,size(locations)
    write (*,*) "Values: ", twse(i), telev(i), tdepth(i), twse(i)-telev(i)
  enddo 


  write (*,"(a)", advance="no") "Finalizing..."
  s = m%finalize()
  write (*,*) "Done"


    end program GetValueTest

