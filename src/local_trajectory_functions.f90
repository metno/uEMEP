module local_trajectory_functions
    !! This module contains subroutines for calculating local emission trajectories

    use uEMEP_definitions
    use uemep_constants, only: NODATA_value, epsilon0
    use utility_functions, only: distrl_sqr

    implicit none
    private

    public :: calculate_all_trajectory, calculate_minimum_distance_trajectory

contains

    subroutine calculate_all_trajectory(x_emis, y_emis, t, traj_max_index, dr_traj, x_traj, y_traj)
        !! Calculates the trajectory of an emission from a source point
        !!
        !! The subroutine iteratively identifies the meteorological grid the trajectory
        !! point is in and update the trajectory coordinates (`x_traj` and `y_traj`)
        !!
        !! The subroutine will return if the trajectory point moves outside the grid, or
        !! if the maximum number of trajectory points (`traj_max_index`) is reached
        real, intent(in) :: x_emis !! Emission source x-coordinate
        real, intent(in) :: y_emis !! Emission source y-coordinate
        integer, intent(in) :: t !! Time index
        integer, intent(in) :: traj_max_index !! Maximum number of trajectory points
        real, intent(in) :: dr_traj !! Trajectory step size
        real, intent(out) :: x_traj(traj_max_index) !! Trajectory point x-coordinates
        real, intent(out) :: y_traj(traj_max_index) !! Trajectory point y-coordinates

        ! Local variables
        integer :: k, i_integral, j_integral
        logical :: exit_traj

        ! Initialize trajectory position to the emission source
        k = 1
        x_traj(k) = x_emis
        y_traj(k) = y_emis

        ! Calculate trajectory for each point in the grid
        exit_traj = .false.
        do while (.not. exit_traj .and. k < traj_max_index)
            ! Identify the meteorological grid the current trajectory point is in
            i_integral = 1 + floor((x_traj(k) - integral_subgrid_min(x_dim_index))/integral_subgrid_delta(x_dim_index))
            j_integral = 1 + floor((y_traj(k) - integral_subgrid_min(y_dim_index))/integral_subgrid_delta(y_dim_index))

            ! Check if the trajectory point is still within the grid boundaries
            if (i_integral >= 1 .and. i_integral <= integral_subgrid_dim(x_dim_index) &
                .and. j_integral >= 1 .and. j_integral <= integral_subgrid_dim(y_dim_index)) then

                k = k + 1
                ! Update trajectory point coordinates based on current wind direction
                x_traj(k) = x_traj(k-1) + dr_traj*meteo_subgrid(i_integral,j_integral,t,cos_subgrid_index)
                y_traj(k) = y_traj(k-1) + dr_traj*meteo_subgrid(i_integral,j_integral,t,sin_subgrid_index)
            else
                ! Exit the loop if the trajectory point has moved outside the grid
                exit_traj = .true.
            end if
        end do
    end subroutine calculate_all_trajectory

    subroutine calculate_minimum_distance_trajectory(x_r, y_r, traj_max_index, dr_traj, x_traj, y_traj, x_loc, y_loc, valid_traj)
        !! Calculates the minimum distance from a given point to a set of trajectory points
        !!
        !! This subroutine calculates the minimum squared distance from a given point (`x_r`, `y_r`)
        !! to a set of trajectory points (`x_traj`, `y_traj`), and determines the trajectory point 
        !! that's closest to the given point
        !!
        !! If the given point coincides with an emission point, the subroutine returns
        !! If the minimum distance is upwind, the `valid_traj` flag is set to false
        real, intent(in) :: x_r !! Point x-coordinate
        real, intent(in) :: y_r !! Point y-coordinate
        integer, intent(in) :: traj_max_index !! Maximum number of trajectory points
        real, intent(in) :: dr_traj !! Trajectory step size
        real, intent(in) :: x_traj(:) !! Trajectory point x-coordinates
        real, intent(in) :: y_traj(:) !! Trajectory point y-coordinates
        real, intent(out) :: x_loc !! Cumulative distance from the emission source to the closest trajectory point
        real, intent(out) :: y_loc !! Minimum distance from the given point to the closest trajectory point
        logical, intent(out) :: valid_traj !! A logical indicating whether the closest trajectory point is valid
        
        ! Local variables
        integer :: k
        real :: distance_traj, distance_intercept_traj
        real :: x_intercept_traj, y_intercept_traj, frac_length_traj
        real :: distance_intercept_min

        ! Set the distance for the initial emission grid
        k=1
        distance_traj = (x_traj(k) - x_r)*(x_traj(k) - x_r) + (y_traj(k) - y_r)*(y_traj(k) - y_r)

        ! Return if the receptor is the same as the emission grid
        if (abs(distance_traj) < epsilon0) then
            y_loc = 0.0
            x_loc = 0.0
            valid_traj = .true.
            return
        end if

        ! Reset return values
        y_loc = 0.0
        x_loc = 0.0
        valid_traj = .false.
        distance_intercept_min = distance_traj

        ! Calculate distance for each trajectory point
        do k = 2, traj_max_index
            if (x_traj(k) /= NODATA_value) then
                call distrl_sqr(x_r, y_r, x_traj(k-1), y_traj(k-1), x_traj(k), y_traj(k), &
                    x_intercept_traj, y_intercept_traj, distance_intercept_traj, frac_length_traj)

                if (distance_intercept_traj < distance_intercept_min) then
                    distance_intercept_min = distance_intercept_traj
                    y_loc = distance_intercept_traj
                    x_loc = dr_traj*(k-2) + frac_length_traj*dr_traj
                    valid_traj = .true.
                end if
            end if
        end do

        ! Invalidate results if they are upwind
        if (abs(x_loc) < epsilon0 .and. y_loc > dr_traj) then
            valid_traj = .false.
        end if

        ! Convert to actual distance
        y_loc = sqrt(y_loc)
    end subroutine calculate_minimum_distance_trajectory

end module local_trajectory_functions

