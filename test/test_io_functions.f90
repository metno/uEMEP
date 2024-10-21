program test_io_functions

    use io_functions

    implicit none

    ! Local variables
    character(len=256) :: pwd
    logical :: expected, result, ok
    ok = .true.

    ! Get current work directory
    call get_command_argument(0, pwd)
    
    ! Test case 1 (check_dir_exist): Directory does not exist
    expected = .false.
    result = check_dir_exist(path=pwd//"nothing_here/")
    if (result /= expected) then
        ok = .false.
        print "(a)", "Test case 1 failed! Routine: check_dir_exist"
    end if

    ! Test case 2 (check_dir_exist): Directory exists
    expected = .true.
    result = check_dir_exist(path=pwd)
    if (result /= expected) then
        ok = .false.
        print "(a)", "Test case 2 failed! Routine: check_dir_exist"
    end if

    ! Test case 3 (check_dir_exist): Directory exist, but no write access
    ! This will fail is tests are run as root!
    expected = .false.
    result = check_dir_exist(path="/")
    if (result /= expected) then
        ok = .false.
        print "(a)", "Test case 3 failed! Routine: check_dir_exist"
    end if

    ! Return test results summary
    if (ok) then
        print "(a)", "test_io_functions: All tests passed."
    else
        print "(a)", "test_io_functions: One ore more tests failed."
        stop 1
    end if
    
end program test_io_functions