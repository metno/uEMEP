module error_handling

    !! This module provides general error handling, checking, and assertion procedures
    !!
    !! Copyright (C) 2007 Free Software Foundation.
    !! License GNU LGPL-3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>.
    !! This is free software: you are free to change and redistribute it.
    !!
    !! Developed and maintained at the Norwegian Meteorological Institute.
    !! Contribute at: <https://github.com/metno/uEMEP>

    use uEMEP_definitions, only: unit_logfile
    use uemep_constants, only: dp

    implicit none

    private

    ! Error codes
    integer, parameter, public :: no_error = 0
    integer, parameter, public :: default_error = -1
    integer, parameter, public :: file_not_found = -2
    integer, parameter, public :: read_error = -3

    ! Precision tolerances
    real, parameter, public :: tol_real = 1.0e-5
    real, parameter, public :: tol_dp = 1.0e-12

    ! Public interfaces
    public :: assert, check_equality

    interface assert
        !! Asserts a condition, and raises an error if violated
        module procedure assert_true
        module procedure assert_code
    end interface

    interface check_equality
        !! Checks if two values are equal within precision tolerances
        module procedure check_equality_integer
        module procedure check_equality_real
        module procedure check_equality_dp
    end interface check_equality

    interface is_between
        !! Checks if a value is between a minimum and maximum threshold (min <= value <= max)
        module procedure is_between_integer
        module procedure is_between_real
        module procedure is_between_dp
    end interface is_between

contains

    subroutine print_error(message)
        !! Writes an error message to the log
        character(len=*), intent(in) :: message

        write(unit_logfile, "(3a)") "ERROR: ", message
    end subroutine print_error

    logical function check_equality_integer(first_value, second_value) result(are_equal)
        !! Compares two integers for equality
        integer, intent(in) :: first_value
        integer, intent(in) :: second_value

        are_equal = (first_value == second_value)
    end function check_equality_integer

    logical function check_equality_real(first_value, second_value) result(are_equal)
        !! Compares two single precision reals for equality within a set tolerance
        real, intent(in) :: first_value
        real, intent(in) :: second_value

        are_equal = (abs(second_value - first_value) <= tol_real)
    end function check_equality_real

    logical function check_equality_dp(first_value, second_value) result(are_equal)
        !! Compares two double precision reals for equality within a set tolerance
        real(dp), intent(in) :: first_value
        real(dp), intent(in) :: second_value

        are_equal = (abs(second_value - first_value) <= tol_dp)
    end function check_equality_dp

    logical function is_between_integer(value, min_threshold, max_threshold) result(is_within_range)
        !! Checks if an integer value is between a minimum and maximum threshold (min <= value <= max)
        integer, intent(in) :: value
        integer, intent(in) :: min_threshold
        integer, intent(in) :: max_threshold

        call assert((min_threshold > max_threshold), "Minimum threshold cannot exceed maximum threshold")
        is_within_range = (value >= min_threshold .and. value <= max_threshold)
    end function is_between_integer

    logical function is_between_real(value, min_threshold, max_threshold) result(is_within_range)
        !! Checks if a real value is between a minimum and maximum threshold (min <= value <= max)
        real, intent(in) :: value
        real, intent(in) :: min_threshold
        real, intent(in) :: max_threshold

        call assert((min_threshold > max_threshold), "Minimum threshold cannot exceed maximum threshold")
        is_within_range = (value >= min_threshold .and. value <= max_threshold)
    end function is_between_real

    logical function is_between_dp(value, min_threshold, max_threshold) result(is_within_range)
        !! Checks if a double precision real value is between a minimum and maximum threshold (min <= value <= max)
        real(dp), intent(in) :: value
        real(dp), intent(in) :: min_threshold
        real(dp), intent(in) :: max_threshold

        call assert((min_threshold > max_threshold), "Minimum threshold cannot exceed maximum threshold")
        is_within_range = (value >= min_threshold .and. value <= max_threshold)
    end function is_between_dp

    subroutine assert_true(condition, message, code)
        !! Asserts a condition and terminates with an error message and code if violated
        logical, intent(in) :: condition
        character(len=*), intent(in), optional :: message
        integer, intent(in), optional :: code

        integer :: error_code
        if (.not. condition) then
            if (present(code)) then
                error_code = code
            else
                error_code = default_error
            end if

            if (present(message)) then
                call print_error(message)
            else
                call print_error("Assertion failed!")
            end if
            stop error_code
        end if
    end subroutine assert_true

    subroutine assert_code(code, message, no_error_code)
        !! Asserts the error code, and terminates with an error message if different from no_error
        integer, intent(in) :: code
        character(len=*), intent(in), optional :: message
        integer, intent(in), optional :: no_error_code

        logical :: assert_condition

        integer :: local_error_code
        if (present(no_error_code)) then
            local_error_code = no_error_code
        else
            local_error_code = no_error
        end if

        assert_condition = check_equality(code, local_error_code)

        if (present(message)) then
            call assert(assert_condition, message, code)
        else
            call assert(assert_condition, code=code)
        end if
    end subroutine assert_code

end module error_handling