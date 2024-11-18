module area_interpolation_functions

    implicit none
    private

    public :: area_weighted_interpolation_function, &
        area_weighted_extended_interpolation_function, &
        area_weighted_extended_vectorgrid_interpolation_function

contains

    function area_weighted_interpolation_function(xgrid, ygrid, zgrid, delta, xval, yval) result(res)
        !! Returns the area weighted value for a point with fixed delta
        !!
        !! This function performs an area-weighted interpolation to estimate the value of a field 
        !! at a specific point `(xval, yval)` based on a grid of values (`zgrid`) and positions
        !! (`xgrid` and `ygrid`)
        !!
        !! If no valid weights are found (`sum_weight == 0.0`), the function returns `0.0`
        real, intent(in) :: xgrid(:,:) !! Array of x-coordinates of the grid points
        real, intent(in) :: ygrid(:,:) !! Array of y-coordinates of the grid points
        real, intent(in) :: zgrid(:,:) !! Array of values at each grid point to be interpolated
        real, intent(in) :: delta(2) !! Array specifying the grid resolution
        real, intent(in) :: xval !! The x-coordinate of the target point for interpolation
        real, intent(in) :: yval !! The y-coordinate of the target point for interpolation
        real :: res !! Interpolated value of `zgrid` at point `(xval, yval)`

        ! Local variables
        integer :: xdim, ydim
        integer :: i, j, ii, jj, iii, jjj
        real :: zval, sum_weight, weighting
        real :: xpos_area_max, xpos_area_min, ypos_area_max, ypos_area_min
        real :: xpos_max, xpos_min, ypos_max, ypos_min

        ! Find the size of the arrays
        xdim = size(xgrid, 1)
        ydim = size(xgrid, 2)

        ! If only on grid available then return the value of that grid
        if (xdim == 1 .and. ydim == 1) then
            res = zgrid(xdim,ydim)
            return
        end if

        ! Find grid index for position val
        i = 1 + floor((xval - xgrid(1,1))/delta(1) + 0.5)
        j = 1 + floor((yval - ygrid(1,1))/delta(2) + 0.5)

        if (i < 1 .or. j < 1 .or. i > xdim .or. j > ydim) then
            write(*,'(A,4i6)') 'Interpolation out of range in area_weighted_interpolation_function. Stopping. (i,j,xdim,ydim)', i, j, xdim, ydim
            write(*,'(4f12.2)') xval, yval, xgrid(1,1), ygrid(1,1)
            stop 1
        else
            xpos_area_max = xval + delta(1)/2.0
            xpos_area_min = xval - delta(1)/2.0
            ypos_area_max = yval + delta(2)/2.0
            ypos_area_min = yval - delta(2)/2.0

            zval = 0.0
            sum_weight = 0.0
            do jjj = j-1, j+1
                do iii = i-1, i+1
                    jj = max(jjj,1); jj = min(jj,ydim)
                    ii = max(iii,1); ii = min(ii,xdim)
                    xpos_min = max(xpos_area_min, xgrid(ii,jj) - delta(1)/2.0)
                    xpos_max = min(xpos_area_max, xgrid(ii,jj) + delta(1)/2.0)
                    ypos_min = max(ypos_area_min, ygrid(ii,jj) - delta(2)/2.0)
                    ypos_max = min(ypos_area_max, ygrid(ii,jj) + delta(2)/2.0)

                    if (xpos_max > xpos_min .and. ypos_max > ypos_min) then
                        weighting = (ypos_max - ypos_min)*(xpos_max - xpos_min)/delta(1)/delta(2)
                    else
                        weighting = 0.0
                    end if

                    zval = zval + zgrid(ii,jj)*weighting
                    sum_weight = sum_weight + weighting
                end do
            end do
        end if
        if (sum_weight > 0.0) then
            res = zval/sum_weight
        else
            res = 0.0
        end if
    end function area_weighted_interpolation_function

    function area_weighted_extended_interpolation_function(xgrid, ygrid, zgrid, delta, xval, yval, delta_val) result(res)
        !! Returns the area weighted value for a point with variable delta
        !!
        !! This function performs an area-weighted interpolation to estimate the value of a field
        !! at a specific point `(xval,yval)` within an rectangle `delta_val(2)` based on grid of
        !! values (`zgrid`) and positions (`xgrid` and `ygrid`)
        !!
        !! If no valid weights are found (`sum_weight == 0.0`), the function returns `0.0`
        real, intent(in) :: xgrid(:,:) !! Array of x-coordinates of the grid points
        real, intent(in) :: ygrid(:,:) !! Array of y-coordinates of the grid points
        real, intent(in) :: zgrid(:,:) !! Array of values at each grid point to be interpolated
        real, intent(in) :: delta(2) !! Array specifying the grid resolution
        real, intent(in) :: delta_val(2) !! Array specifying the size of the grid used for interpolation
        real, intent(in) :: xval !! The x-coordinate of the target point for interpolation
        real, intent(in) :: yval !! The y-coordinate of the target point for interpolation
        real :: res !! Interpolated value of `zgrid` at point `(xval, yval)`

        ! Local variables
        real :: zval, sum_weight, weighting
        real :: xpos_max, xpos_min, ypos_max, ypos_min
        real :: xpos_area_max, xpos_area_min, ypos_area_max, ypos_area_min
        integer :: xdim, ydim
        integer :: i, j, ii, jj, iii, jjj
        integer :: ii_delta, jj_delta

        ! Find the size of the arrays
        xdim = size(xgrid, 1)
        ydim = size(ygrid, 2)

        ! If only on grid available then return the value of that grid
        if (xdim == 1 .and. ydim == 1) then
            res = zgrid(xdim,ydim)
            return
        end if

        ! Find grid index for position val
        i = 1 + floor((xval - xgrid(1,1))/delta(1) + 0.5)
        j = 1 + floor((yval - ygrid(1,1))/delta(2) + 0.5)

        if (i < 1 .or. j < 1 .or. i > xdim .or. j > ydim) then
            write(*,'(A,4i6)') 'Interpolation out of range in area_weighted_extended_interpolation_function. Stopping. (i,j,xdim,ydim)', i, j, xdim, ydim
            write(*,'(4f12.2)') xval, yval, xgrid(1,1), ygrid(1,1)
            stop 1
        else
            xpos_area_max = xval + delta_val(1)/2.0
            xpos_area_min = xval - delta_val(1)/2.0
            ypos_area_max = yval + delta_val(2)/2.0
            ypos_area_min = yval - delta_val(2)/2.0

            jj_delta = 1 + floor(0.5*(delta_val(2)/delta(2) - 1.0))
            ii_delta = 1 + floor(0.5*(delta_val(1)/delta(1) - 1.0))

            zval = 0.0
            sum_weight = 0.0
            do jjj = j-jj_delta, j+jj_delta
                do iii = i-ii_delta, i+ii_delta
                    jj = max(jjj,1); jj = min(jj,ydim)
                    ii = max(iii,1); ii = min(ii,xdim)
                    xpos_min = max(xpos_area_min, xgrid(ii,jj) - delta(1)/2.0)
                    xpos_max = min(xpos_area_max, xgrid(ii,jj) + delta(1)/2.0)
                    ypos_min = max(ypos_area_min, ygrid(ii,jj) - delta(2)/2.0)
                    ypos_max = min(ypos_area_max, ygrid(ii,jj) + delta(2)/2.0)

                    if (xpos_max > xpos_min .and. ypos_max > ypos_min) then
                        weighting = (ypos_max - ypos_min)*(xpos_max - xpos_min)/delta_val(1)/delta_val(2)
                    else
                        weighting = 0.0
                    endif
                    zval = zval + zgrid(ii,jj)*weighting
                    sum_weight = sum_weight + weighting
                end do
            end do
        end if
        if (sum_weight > 0.0) then
            res = zval/sum_weight
        else
            res = 0.0
        end if
    end function area_weighted_extended_interpolation_function

    function area_weighted_extended_vectorgrid_interpolation_function(xgrid, ygrid, zgrid, delta, xval, yval, delta_val) result(res)
        !! Returns the area weighted value for a point with variable delta (vector version)
        !!
        !! This function performs an area-weighted interpolation to estimate the value of a field
        !! at a specific point `(xval,yval)` within an rectangle `delta_val(2)` based on grid of
        !! values (`zgrid`) and vectors of positions (`xgrid` and `ygrid`)
        !!
        !! If no valid weights are found (`sum_weight == 0.0`), the function returns `0.0`
        real, intent(in) :: xgrid(:) !! Vector of x-coordinates of the grid points
        real, intent(in) :: ygrid(:) !! Vector of y-coordinates of the grid points
        real, intent(in) :: zgrid(:,:) !! Array of values at each grid point to be interpolated
        real, intent(in) :: delta(2) !! Array specifying the grid resolution
        real, intent(in) :: delta_val(2) !! Array specifying the size of the grid used for interpolation
        real, intent(in) :: xval !! The x-coordinate of the target point for interpolation
        real, intent(in) :: yval !! The y-coordinate of the target point for interpolation
        real :: res !! Interpolated value of `zgrid` at point `(xval, yval)`

        ! Local variables
        real :: zval, sum_weight, weighting
        real :: xpos_area_max, xpos_area_min, ypos_area_max, ypos_area_min
        real :: xpos_max, xpos_min, ypos_max, ypos_min
        integer :: xdim, ydim
        integer :: i, j, ii, jj, iii, jjj
        integer :: ii_delta, jj_delta

        ! Find the size of the arrays
        xdim = size(xgrid)
        ydim = size(ygrid)

        ! If only on grid available then return the value of that grid
        if (xdim == 1 .and. ydim == 1) then
            res = zgrid(xdim,ydim)
            return
        end if

        ! Find grid index for position val
        i = 1 + floor((xval-xgrid(1))/delta(1) + 0.5)
        j = 1 + floor((yval-ygrid(1))/delta(2) + 0.5)

        if (i < 1 .or. j < 1 .or. i > xdim .or. j > ydim) then
            write(*,'(A,4i6)') 'Interpolation out of range in area_weighted_extended_interpolation_function. Stopping. (i,j,xdim,ydim)', i, j, xdim, ydim
            write(*,'(4f12.2)') xval, yval, xgrid(1), ygrid(1)
            stop
        else
            xpos_area_max = xval + delta_val(1)/2.0
            xpos_area_min = xval - delta_val(1)/2.0
            ypos_area_max = yval + delta_val(2)/2.0
            ypos_area_min = yval - delta_val(2)/2.0

            jj_delta = 1 + floor(0.5*(delta_val(2)/delta(2) - 1.0))
            ii_delta = 1 + floor(0.5*(delta_val(1)/delta(1) - 1.0))

            zval = 0.0
            sum_weight = 0.0
            do jjj = j-jj_delta, j+jj_delta
                do iii = i-ii_delta, i+ii_delta
                    jj = max(jjj,1); jj = min(jj,ydim)
                    ii = max(iii,1); ii = min(ii,xdim)
                    xpos_min = max(xpos_area_min, xgrid(ii) - delta(1)/2.0)
                    xpos_max = min(xpos_area_max, xgrid(ii) + delta(1)/2.0)
                    ypos_min = max(ypos_area_min, ygrid(jj) - delta(2)/2.0)
                    ypos_max = min(ypos_area_max, ygrid(jj) + delta(2)/2.0)

                    if (xpos_max > xpos_min .and. ypos_max > ypos_min) then
                        weighting = (ypos_max - ypos_min)*(xpos_max - xpos_min)/delta_val(1)/delta_val(2)
                    else
                        weighting = 0.0
                    end if

                    zval = zval + zgrid(ii,jj)*weighting
                    sum_weight = sum_weight + weighting
                end do
            end do
        end if
        if (sum_weight > 0.0) then
            res = zval/sum_weight
        else
            res = 0.0
        end if
    end function area_weighted_extended_vectorgrid_interpolation_function

end module area_interpolation_functions

