Module GeometryMod

Contains
logical function line_exp_is_degenerate_nd ( dim_num, p1, p2 )

!*******************************************************************************
!
!! LINE_EXP_IS_DEGENERATE_ND finds if an explicit line is degenerate in ND.
!
!  Discussion:
!
!    The explicit form of a line in ND is:
!
!      the line through the points P1 and P2.
!
!    An explicit line is degenerate if the two defining points are equal.
!
!  Modified:
!
!    06 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) P1(DIM_NUM), P2(DIM_NUM), two points on the line.
!
!    Output, logical LINE_EXP_IS_DEGENERATE_ND, is TRUE if the line
!    is degenerate.
!
  implicit none
  INTEGER, PARAMETER :: mp = KIND(1.0D0)
  integer dim_num

!  logical line_exp_is_degenerate_nd
  real ( kind = mp ) p1(dim_num)
  real ( kind = mp ) p2(dim_num)

  line_exp_is_degenerate_nd = ( all ( p1(1:dim_num) == p2(1:dim_num) ) )

  return
end function

logical function line_imp_is_degenerate_2d ( a, b, c )

!*******************************************************************************
!
!! LINE_IMP_IS_DEGENERATE_2D finds if an implicit point is degenerate in 2D.
!
!  Discussion:
!
!    The implicit form of a line in 2D is:
!
!      A * X + B * Y + C = 0
!
!  Modified:
!
!    06 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the implicit line parameters.
!
!    Output, logical LINE_IMP_IS_DEGENERATE_2D, is true if the
!    line is degenerate.
!
  implicit none
  INTEGER, PARAMETER :: mp = KIND(1.0D0)
  integer, parameter :: dim_num = 2

  real ( kind = mp ) a
  real ( kind = mp ) b
  real ( kind = mp ) c
!  logical line_imp_is_degenerate_2d

  line_imp_is_degenerate_2d = ( a * a + b * b == 0.0D+00 )

  return
end function

subroutine lines_exp_int_2d ( p1, p2, q1, q2, ival, p )

!*******************************************************************************
!
!! LINES_EXP_INT_2D determines where two explicit lines intersect in 2D.
!
!  Discussion:
!
!    The explicit form of a line in 2D is:
!
!      the line through the points P1 and P2.
!
!  Modified:
!
!    02 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P1(2), P2(2), two points on the first line.
!
!    Input, real ( kind = 8 ) Q1(2), Q2(2), two points on the second line.
!
!    Output, integer IVAL, reports on the intersection:
!    0, no intersection, the lines may be parallel or degenerate.
!    1, one intersection point, returned in P.
!    2, infinitely many intersections, the lines are identical.
!
!    Output, real ( kind = 8 ) P(2), if IVAl = 1, P is
!    the intersection point.  Otherwise, P = 0.
!
  implicit none
  INTEGER, PARAMETER :: mp = KIND(1.0D0)
  integer, parameter :: dim_num = 2

  real ( kind = mp ) a1
  real ( kind = mp ) a2
  real ( kind = mp ) b1
  real ( kind = mp ) b2
  real ( kind = mp ) c1
  real ( kind = mp ) c2
  integer ival
  logical point_1
  logical point_2
  real ( kind = mp ) p(dim_num)
  real ( kind = mp ) p1(dim_num)
  real ( kind = mp ) p2(dim_num)
  real ( kind = mp ) q1(dim_num)
  real ( kind = mp ) q2(dim_num)

  ival = 0
  p(1:dim_num) = 0.0D+00
!
!  Check whether either line is a point.
!
  if ( all ( p1(1:dim_num) == p2(1:dim_num) ) ) then
    point_1 = .true.
  else
    point_1 = .false.
  end if

  if ( all ( q1(1:dim_num) == q2(1:dim_num) ) ) then
    point_2 = .true.
  else
    point_2 = .false.
  end if
!
!  Convert the lines to ABC format.
!
  if ( .not. point_1 ) then
    call line_exp2imp_2d ( p1, p2, a1, b1, c1 )
  end if

  if ( .not. point_2 ) then
    call line_exp2imp_2d ( q1, q2, a2, b2, c2 )
  end if
!
!  Search for intersection of the lines.
!
  if ( point_1 .and. point_2 ) then
    if ( all ( p1(1:dim_num) == q1(1:dim_num) ) ) then
      ival = 1
      p(1:dim_num) = p1(1:dim_num)
    end if
  else if ( point_1 ) then
    if ( a2 * p1(1) + b2 * p1(2) == c2 ) then
      ival = 1
      p(1:dim_num) = p1(1:dim_num)
    end if
  else if ( point_2 ) then
    if ( a1 * q1(1) + b1 * q1(2) == c1 ) then
      ival = 1
      p(1:dim_num) = q1(1:dim_num)
    end if
  else
    call lines_imp_int_2d ( a1, b1, c1, a2, b2, c2, ival, p )
  end if

  return
end subroutine

subroutine line_exp2imp_2d ( p1, p2, a, b, c )

!*******************************************************************************
!
!! LINE_EXP2IMP_2D converts an explicit line to implicit form in 2D.
!
!  Discussion:
!
!    The explicit form of a line in 2D is:
!
!      the line through the points P1 and P2.
!
!    The implicit form of a line in 2D is:
!
!      A * X + B * Y + C = 0
!
!  Modified:
!
!    06 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P1(2), P2(2), two points on the line.
!
!    Output, real ( kind = 8 ) A, B, C, the implicit form of the line.
!
  implicit none
  INTEGER, PARAMETER :: mp = KIND(1.0D0)
  integer, parameter :: dim_num = 2

  real ( kind = mp ) a
  real ( kind = mp ) b
  real ( kind = mp ) c
!  logical line_exp_is_degenerate_nd
  real ( kind = mp ) norm
  real ( kind = mp ) p1(dim_num)
  real ( kind = mp ) p2(dim_num)
!
!  Take care of degenerate cases.
!
  if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINE_EXP2IMP_2D - Warning!'
    write ( *, '(a)' ) '  The line is degenerate.'
  end if

  a = p2(2) - p1(2)
  b = p1(1) - p2(1)
  c = p2(1) * p1(2) - p1(1) * p2(2)

  norm = a * a + b * b + c * c

  if ( 0.0D+00 < norm ) then
    a = a / norm
    b = b / norm
    c = c / norm
  end if

  if ( a < 0.0D+00 ) then
    a = -a
    b = -b
    c = -c
  end if

  return
end Subroutine
subroutine lines_imp_int_2d ( a1, b1, c1, a2, b2, c2, ival, p )

!*******************************************************************************
!
!! LINES_IMP_INT_2D determines where two implicit lines intersect in 2D.
!
!  Discussion:
!
!    The implicit form of a line in 2D is:
!
!      A * X + B * Y + C = 0
!
!  Modified:
!
!    25 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A1, B1, C1, define the first line.
!    At least one of A1 and B1 must be nonzero.
!
!    Input, real ( kind = 8 ) A2, B2, C2, define the second line.
!    At least one of A2 and B2 must be nonzero.
!
!    Output, integer IVAL, reports on the intersection.
!
!    -1, both A1 and B1 were zero.
!    -2, both A2 and B2 were zero.
!     0, no intersection, the lines are parallel.
!     1, one intersection point, returned in P.
!     2, infinitely many intersections, the lines are identical.
!
!    Output, real ( kind = 8 ) P(2), if IVAL = 1, then P is
!    the intersection point.  Otherwise, P = 0.
!
  implicit none
  INTEGER, PARAMETER :: mp = KIND(1.0D0)
  integer, parameter :: dim_num = 2

  real ( kind = mp ) a(dim_num,dim_num+1)
  real ( kind = mp ) a1
  real ( kind = mp ) a2
  real ( kind = mp ) b1
  real ( kind = mp ) b2
  real ( kind = mp ) c1
  real ( kind = mp ) c2
  integer info
  integer ival
!  logical line_imp_is_degenerate_2d
  real ( kind = mp ) p(dim_num)

  p(1:dim_num) = 0.0D+00
!
!  Refuse to handle degenerate lines.
!
  if ( line_imp_is_degenerate_2d ( a1, b1, c1 ) ) then
    ival = -1
    return
  end if

  if ( line_imp_is_degenerate_2d ( a2, b2, c2 ) ) then
    ival = -2
    return
  end if
!
!  Set up and solve a linear system.
!
  a(1,1) = a1
  a(1,2) = b1
  a(1,3) = -c1

  a(2,1) = a2
  a(2,2) = b2
  a(2,3) = -c2

  call dmat_solve ( 2, 1, a, info )
!
!  If the inverse exists, then the lines intersect at the solution point.
!
  if ( info == 0 ) then

    ival = 1
    p(1:dim_num) = a(1:dim_num,3)
!
!  If the inverse does not exist, then the lines are parallel
!  or coincident.  Check for parallelism by seeing if the
!  C entries are in the same ratio as the A or B entries.
!
  else

    ival = 0

    if ( a1 == 0.0D+00 ) then
      if ( b2 * c1 == c2 * b1 ) then
        ival = 2
      end if
    else
      if ( a2 * c1 == c2 * a1 ) then
        ival = 2
      end if
    end if

  end if

  return
end subroutine

subroutine dmat_solve ( n, rhs_num, a, info )

!*******************************************************************************
!
!! DMAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.
!
!  Modified:
!
!    29 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer rhs_num, the number of right hand sides.  rhs_num
!    must be at least 0.
!
!    Input/output, real ( kind = 8 ) A(N,N+rhs_num), contains in rows and columns 1
!    to N the coefficient matrix, and in columns N+1 through
!    N+rhs_num, the right hand sides.  On output, the coefficient matrix
!    area has been destroyed, while the right hand sides have
!    been overwritten with the corresponding solutions.
!
!    Output, integer INFO, singularity flag.
!    0, the matrix was not singular, the solutions were computed;
!    J, factorization failed on step J, and the solutions could not
!    be computed.
!
  implicit none
  INTEGER, PARAMETER :: mp = KIND(1.0D0)
  integer n
  integer rhs_num

  real ( kind = mp ) a(n,n+rhs_num)
  real ( kind = mp ) apivot
  real ( kind = mp ) factor
  integer i
  integer info
  integer ipivot
  integer j

  info = 0

  do j = 1, n
!
!  Choose a pivot row.
!
    ipivot = j
    apivot = a(j,j)

    do i = j+1, n
      if ( abs ( apivot ) < abs ( a(i,j) ) ) then
        apivot = a(i,j)
        ipivot = i
      end if
    end do

    if ( apivot == 0.0D+00 ) then
      info = j
      return
    end if
!
!  Interchange.
!
    do i = 1, n + rhs_num
      call d_swap ( a(ipivot,i), a(j,i) )
    end do
!
!  A(J,J) becomes 1.
!
    a(j,j) = 1.0D+00
    a(j,j+1:n+rhs_num) = a(j,j+1:n+rhs_num) / apivot
!
!  A(I,J) becomes 0.
!
    do i = 1, n

      if ( i /= j ) then

        factor = a(i,j)
        a(i,j) = 0.0D+00
        a(i,j+1:n+rhs_num) = a(i,j+1:n+rhs_num) - factor * a(j,j+1:n+rhs_num)

      end if

    end do

  end do

  return
end subroutine

subroutine d_swap ( x, y )

!*******************************************************************************
!
!! D_SWAP switches two real values.
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none
  INTEGER, PARAMETER :: mp = KIND(1.0D0)
  real ( kind = mp ) x
  real ( kind = mp ) y
  real ( kind = mp ) z

  z = x
  x = y
  y = z

  return
end subroutine

End Module
