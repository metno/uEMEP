program test_area_interpolation_functions

    use area_interpolation_functions

    implicit none

    real :: xgrid(3,3), ygrid(3,3), zgrid(3,3), delta(2), delta_val(2)
    real :: xgrid_vec(3), ygrid_vec(3)
    real :: xval, yval, result, expected
    logical :: ok
    integer :: xdim, ydim

    ok = .true.
    xdim = 3
    ydim = 3

    ! Define a simple grid for testing
    xgrid = reshape([1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0], shape(xgrid))
    ygrid = reshape([1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0], shape(ygrid))
    zgrid = reshape([10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0], shape(zgrid))
    delta = [1.0, 1.0]
    delta_val = [1.0, 1.0]

    ! Define vector grids for testing
    xgrid_vec = [1.0, 2.0, 3.0]
    ygrid_vec = [1.0, 2.0, 3.0]

    ! Test case 1: Point-based area interpolation in the center of the grid
    xval = 2.0
    yval = 2.0
    expected = 50.0
    result = area_weighted_interpolation_function(xgrid, ygrid, zgrid, delta, xval, yval)
    if (abs(result - expected) > 1.0e-5) then
        ok = .false.
        print "(a)", "Test case 1 failed! Routine: area_weighted_interpolation_function"
    end if

    ! Test case 2: Point-based area interpolation at the upper-left edge of the grid
    xval = 0.50001
    yval = 0.50001
    expected = 10.0
    result = area_weighted_interpolation_function(xgrid, ygrid, zgrid, delta, xval, yval)
    if (abs(result - expected) > 1.0e-5) then
        ok = .false.
        print "(a)", "Test case 2 failed! Routine: area_weighted_interpolation_function"
    end if

    ! Test case 3: Point-based area interpolation at the lower-right edge of the grid
    xval = 3.49999
    yval = 3.49999
    expected = 90.0
    result = area_weighted_interpolation_function(xgrid, ygrid, zgrid, delta, xval, yval)
    if (abs(result - expected) > 1.0e-5) then
        ok = .false.
        print "(a)", "Test case 3 failed! Routine: area_weighted_interpolation_function"
    end if

    ! Test case 4: Rectangle-based area interpolation in the center of the grid
    xval = 2.0
    yval = 2.0
    delta_val = [2.0, 2.0]
    expected = 50.0
    result = area_weighted_extended_interpolation_function(xgrid, ygrid, zgrid, delta, xval, yval, delta_val)
    if (abs(result - expected) > 1.0e-5) then
        ok = .false.
        print "(a)", "Test case 4 failed! Routine: area_weighted_extended_interpolation_function"
    end if

    ! Test case 5: Rectangle-based area interpolation at the upper-left edge of the grid
    xval = 0.50001
    yval = 0.50001
    expected = 10.0002
    result = area_weighted_extended_interpolation_function(xgrid, ygrid, zgrid, delta, xval, yval, delta_val)
    if (abs(result - expected) > 1.0e-5) then
        ok = .false.
        print "(a)", "Test case 5 failed! Routine: area_weighted_extended_interpolation_function"
    end if

    ! Test case 6: Rectangle-based area interpolation at the lower-right edge of the grid
    xval = 3.49999
    yval = 3.49999
    expected = 89.99980
    result = area_weighted_extended_interpolation_function(xgrid, ygrid, zgrid, delta, xval, yval, delta_val)
    if (abs(result - expected) > 1.0e-5) then
        ok = .false.
        print "(a)", "Test case 6 failed! Routine: area_weighted_extended_interpolation_function"
    end if

    ! Test case 7: Rectangle-based area interpolation with smaller delta_val
    xval = 2.15
    yval = 1.8
    delta_val = [1.0, 1.0]
    expected = 45.5
    result = area_weighted_extended_interpolation_function(xgrid, ygrid, zgrid, delta, xval, yval, delta_val)
    if (abs(result - expected) > 1.0e-5) then
        ok = .false.
        print "(a)", "Test case 7 failed! Routine: area_weighted_extended_interpolation_function"
    end if

    ! Test case 8: Vector-based area interpolation in the center of the grid
    xval = 2.0
    yval = 2.0
    delta_val = [2.0, 2.0]
    expected = 50.0
    result = area_weighted_extended_vectorgrid_interpolation_function(xgrid_vec, ygrid_vec, zgrid, delta, xval, yval, delta_val)
    if (abs(result - expected) > 1.0e-5) then
        ok = .false.
        print "(a)", "Test case 8 failed! Routine: area_weighted_extended_vectorgrid_interpolation_function"
    end if

    ! Test case 9: Vector-based area interpolation at the upper-left edge of the grid
    xval = 0.50001
    yval = 0.50001
    expected = 10.0002
    result = area_weighted_extended_vectorgrid_interpolation_function(xgrid_vec, ygrid_vec, zgrid, delta, xval, yval, delta_val)
    if (abs(result - expected) > 1.0e-5) then
        ok = .false.
        print "(a)", "Test case 9 failed! Routine: area_weighted_extended_vectorgrid_interpolation_function"
    end if

    ! Test case 10: Vector-based area interpolation at the lower-right edge of the grid
    xval = 3.49999
    yval = 3.49999
    expected = 89.99980
    result = area_weighted_extended_vectorgrid_interpolation_function(xgrid_vec, ygrid_vec, zgrid, delta, xval, yval, delta_val)
    if (abs(result - expected) > 1.0e-5) then
        ok = .false.
        print "(a)", "Test case 10 failed! Routine: area_weighted_extended_vectorgrid_interpolation_function"
    end if

    ! Test case 11: Vector-based area interpolation using smaller delta_val
    xval = 1.88
    yval = 1.93
    delta_val = [1.0, 1.0]
    expected = 46.7
    result = area_weighted_extended_vectorgrid_interpolation_function(xgrid_vec, ygrid_vec, zgrid, delta, xval, yval, delta_val)
    if (abs(result - expected) > 1.0e-5) then
        ok = .false.
        print "(a)", "Test case 11 failed! Routine: area_weighted_extended_vectorgrid_interpolation_function"
    end if

    ! Return the test results summary
    if (ok) then
        print "(a)", "test_area_interpolation_functions: All tests passed."
    else
        print "(a)", "test_area_interpolation_functions: One or more tests failed."
        stop 1
    end if

end program test_area_interpolation_functions