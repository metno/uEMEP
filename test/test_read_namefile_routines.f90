program test_namefile_routines

    use read_namefile_routines

    implicit none
    
    logical :: ok = .true.

    integer, parameter :: dp = selected_real_kind(15, 307)

    integer :: unit_in, unit_out
    real :: real_default, real_out, real_test
    real(dp) :: dp_default, dp_out, dp_test
    integer :: int_default, int_out, int_test
    logical :: bool_default, bool_out, bool_test
    character(len=:), allocatable :: char_default, char_out, char_test
    character(len=:), allocatable :: str_in
    character(len=2048) :: path, relpath

    ! Test case 1 (read_name_real): Check that indented variable names are found
    str_in = "param_real_2"
    real_out = -999.0
    real_test = 2.0
    real_default = 42.0
    unit_in = 10
    unit_out = 20
    call get_environment_variable("PWD", path)
    relpath = trim(path)//"/../test/test_data/"
    open(unit=unit_in, file=trim(relpath)//"test_read_namelist.txt", status="old")
    real_out = read_name_real(str_in, real_default, unit_in, unit_out)
    close(unit_in)
    deallocate(str_in)
    if (abs(real_out - real_test) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 1 failed! Routine: read_name_real"
    end if

    ! Test case 2 (read_name_real): Return default value if variable name is not found
    str_in = "non-existent-variable"
    real_out = -999.0
    real_test = 42.0
    real_default = 42.0
    unit_in = 10
    unit_out = 20
    call get_environment_variable("PWD", path)
    relpath = trim(path)//"/../test/test_data/"
    open(unit=unit_in, file=trim(relpath)//"test_read_namelist.txt", status="old")
    real_out = read_name_real(str_in, real_default, unit_in, unit_out)
    close(unit_in)
    deallocate(str_in)
    if (abs(real_out - real_test) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 2 failed! Routine: read_name_real"
    end if

    ! Test case 3 (read_name_real): Check that commented lines are skipped
    str_in = "param_real_1"
    real_out = -999.0
    real_test = 11.0
    real_default = 42.0
    unit_in = 10
    unit_out = 20
    call get_environment_variable("PWD", path)
    relpath = trim(path)//"/../test/test_data/"
    open(unit=unit_in, file=trim(relpath)//"test_read_namelist.txt", status="old")
    real_out = read_name_real(str_in, real_default, unit_in, unit_out)
    close(unit_in)
    deallocate(str_in)
    if (abs(real_out - real_test) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 3 failed! Routine: read_name_real"
    end if

    ! Test case 4 (read_name_real): Check that right value of named variable is returned
    str_in = "param_real_3"
    real_out = -999.0
    real_test = 3.0
    real_default = 42.0
    unit_in = 10
    unit_out = 20
    call get_environment_variable("PWD", path)
    relpath = trim(path)//"/../test/test_data/"
    open(unit=unit_in, file=trim(relpath)//"test_read_namelist.txt", status="old")
    real_out = read_name_real(str_in, real_default, unit_in, unit_out)
    close(unit_in)
    deallocate(str_in)
    if (abs(real_out - real_test) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 4 failed! Routine: read_name_real"
    end if

    ! Test case 5 (read_name_double): Check that indented variable names are found
    str_in = "param_real_2"
    dp_out = -999.0
    dp_test = 2.0
    dp_default = 42.0
    unit_in = 10
    unit_out = 20
    call get_environment_variable("PWD", path)
    relpath = trim(path)//"/../test/test_data/"
    open(unit=unit_in, file=trim(relpath)//"test_read_namelist.txt", status="old")
    dp_out = read_name_double(str_in, dp_default, unit_in, unit_out)
    close(unit_in)
    deallocate(str_in)
    if (abs(dp_out - dp_test) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 5 failed! Routine: read_name_double"
    end if

    ! Test case 6 (read_name_double): Return default value if variable name is not found
    str_in = "non-existent-variable"
    dp_out = -999.0
    dp_test = 42.0
    dp_default = 42.0
    unit_in = 10
    unit_out = 20
    call get_environment_variable("PWD", path)
    relpath = trim(path)//"/../test/test_data/"
    open(unit=unit_in, file=trim(relpath)//"test_read_namelist.txt", status="old")
    dp_out = read_name_double(str_in, dp_default, unit_in, unit_out)
    close(unit_in)
    deallocate(str_in)
    if (abs(dp_out - dp_test) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 6 failed! Routine: read_name_double"
    end if

    ! Test case 7 (read_name_double): Check that commented lines are skipped
    str_in = "param_real_1"
    dp_out = -999.0
    dp_test = 11.0
    dp_default = 42.0
    unit_in = 10
    unit_out = 20
    call get_environment_variable("PWD", path)
    relpath = trim(path)//"/../test/test_data/"
    open(unit=unit_in, file=trim(relpath)//"test_read_namelist.txt", status="old")
    dp_out = read_name_double(str_in, dp_default, unit_in, unit_out)
    close(unit_in)
    deallocate(str_in)
    if (abs(dp_out - dp_test) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 7 failed! Routine: read_name_double"
    end if

    ! Test case 8 (read_name_double): Check that right value of named variable is returned
    str_in = "param_real_3"
    dp_out = -999.0
    dp_test = 3.0
    dp_default = 42.0
    unit_in = 10
    unit_out = 20
    call get_environment_variable("PWD", path)
    relpath = trim(path)//"/../test/test_data/"
    open(unit=unit_in, file=trim(relpath)//"test_read_namelist.txt", status="old")
    dp_out = read_name_double(str_in, dp_default, unit_in, unit_out)
    close(unit_in)
    deallocate(str_in)
    if (abs(dp_out - dp_test) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 8 failed! Routine: read_name_double"
    end if

    ! Test case 9 (read_name_integer): Check that indented variable names are found
    str_in = "param_int_2"
    int_out = -999
    int_test = 2
    int_default = 42
    unit_in = 10
    unit_out = 20
    call get_environment_variable("PWD", path)
    relpath = trim(path)//"/../test/test_data/"
    open(unit=unit_in, file=trim(relpath)//"test_read_namelist.txt", status="old")
    int_out = read_name_integer(str_in, int_default, unit_in, unit_out)
    close(unit_in)
    deallocate(str_in)
    if (abs(int_out - int_test) /= 0) then
        ok = .false.
        print "(a)", "Test case 9 failed! Routine: read_name_integer"
    end if

    ! Test case 10 (read_name_integer): Return default value if variable name is not found
    str_in = "non-existent-variable"
    int_out = -999
    int_test = 42
    int_default = 42
    unit_in = 10
    unit_out = 20
    call get_environment_variable("PWD", path)
    relpath = trim(path)//"/../test/test_data/"
    open(unit=unit_in, file=trim(relpath)//"test_read_namelist.txt", status="old")
    int_out = read_name_integer(str_in, int_default, unit_in, unit_out)
    close(unit_in)
    deallocate(str_in)
    if (abs(int_out - int_test) /= 0) then
        ok = .false.
        print "(a)", "Test case 10 failed! Routine: read_name_integer"
    end if

    ! Test case 11 (read_name_integer): Check that commented lines are skipped
    str_in = "param_int_1"
    int_out = -999
    int_test = 11
    int_default = 42
    unit_in = 10
    unit_out = 20
    call get_environment_variable("PWD", path)
    relpath = trim(path)//"/../test/test_data/"
    open(unit=unit_in, file=trim(relpath)//"test_read_namelist.txt", status="old")
    int_out = read_name_integer(str_in, int_default, unit_in, unit_out)
    close(unit_in)
    deallocate(str_in)
    if (abs(int_out - int_test) /= 0) then
        ok = .false.
        print "(a)", "Test case 11 failed! Routine: read_name_integer"
    end if

    ! Test case 12 (read_name_integer): Check that right value of named variable is returned
    str_in = "param_int_3"
    int_out = -999
    int_test = 3
    int_default = 42
    unit_in = 10
    unit_out = 20
    call get_environment_variable("PWD", path)
    relpath = trim(path)//"/../test/test_data/"
    open(unit=unit_in, file=trim(relpath)//"test_read_namelist.txt", status="old")
    int_out = read_name_integer(str_in, int_default, unit_in, unit_out)
    close(unit_in)
    deallocate(str_in)
    if (abs(int_out - int_test) /= 0) then
        ok = .false.
        print "(a)", "Test case 12 failed! Routine: read_name_integer"
    end if

    ! Test case 13 (read_name_char): Check that indented variable names are found
    str_in = "param_char_2"
    char_out = "-999"
    char_test = "2"
    char_default = "42"
    unit_in = 10
    unit_out = 20
    call get_environment_variable("PWD", path)
    relpath = trim(path)//"/../test/test_data/"
    open(unit=unit_in, file=trim(relpath)//"test_read_namelist.txt", status="old")
    char_out = read_name_char(str_in, char_default, unit_in, unit_out)
    close(unit_in)
    if (char_out /= char_test) then
        ok = .false.
        print "(a)", "Test case 13 failed! Routine: read_name_char"
    end if
    deallocate(str_in, char_default, char_out, char_test)

    ! Test case 14 (read_name_char): Return default value if variable name is not found
    str_in = "non-existent-variable"
    char_out = "-999"
    char_test = "42"
    char_default = "42"
    unit_in = 10
    unit_out = 20
    call get_environment_variable("PWD", path)
    relpath = trim(path)//"/../test/test_data/"
    open(unit=unit_in, file=trim(relpath)//"test_read_namelist.txt", status="old")
    char_out = read_name_char(str_in, char_default, unit_in, unit_out)
    close(unit_in)
    if (char_out /= char_test) then
        ok = .false.
        print "(a)", "Test case 14 failed! Routine: read_name_char"
    end if
    deallocate(str_in, char_default, char_out, char_test)

    ! Test case 15 (read_name_char): Check that commented lines are skipped
    str_in = "param_char_1"
    char_out = "-999"
    char_test = "11"
    char_default = "42"
    unit_in = 10
    unit_out = 20
    call get_environment_variable("PWD", path)
    relpath = trim(path)//"/../test/test_data/"
    open(unit=unit_in, file=trim(relpath)//"test_read_namelist.txt", status="old")
    char_out = read_name_char(str_in, char_default, unit_in, unit_out)
    close(unit_in)
    if (char_out /= char_test) then
        ok = .false.
        print "(a)", "Test case 15 failed! Routine: read_name_char"
    end if
    deallocate(str_in, char_default, char_out, char_test)

    ! Test case 16 (read_name_char): Check that right value of named variable is returned
    str_in = "param_char_3"
    char_out = "-999"
    char_test = "3"
    char_default = "42"
    unit_in = 10
    unit_out = 20
    call get_environment_variable("PWD", path)
    relpath = trim(path)//"/../test/test_data/"
    open(unit=unit_in, file=trim(relpath)//"test_read_namelist.txt", status="old")
    char_out = read_name_char(str_in, char_default, unit_in, unit_out)
    close(unit_in)
    if (char_out /= char_test) then
        ok = .false.
        print "(a)", "Test case 16 failed! Routine: read_name_char"
    end if
    deallocate(str_in, char_default, char_out, char_test)

    ! Test case 17 (read_name_logical): Check that indented variable names are found
    str_in = "param_logical_2"
    bool_out = .false.
    bool_test = .true.
    bool_default = .false.
    unit_in = 10
    unit_out = 20
    call get_environment_variable("PWD", path)
    relpath = trim(path)//"/../test/test_data/"
    open(unit=unit_in, file=trim(relpath)//"test_read_namelist.txt", status="old")
    bool_out = read_name_logical(str_in, bool_default, unit_in, unit_out)
    close(unit_in)
    if (bool_out .neqv. bool_test) then
        print *, bool_out, bool_test
        ok = .false.
        print "(a)", "Test case 17 failed! Routine: read_name_logical"
    end if
    deallocate(str_in)

    ! Test case 18 (read_name_logical): Return default value if variable name is not found
    str_in = "non-existent-variable"
    bool_out = .false.
    bool_test = .true.
    bool_default = .true.
    unit_in = 10
    unit_out = 20
    call get_environment_variable("PWD", path)
    relpath = trim(path)//"/../test/test_data/"
    open(unit=unit_in, file=trim(relpath)//"test_read_namelist.txt", status="old")
    bool_out = read_name_logical(str_in, bool_default, unit_in, unit_out)
    close(unit_in)
    if (bool_out .neqv. bool_test) then
        ok = .false.
        print "(a)", "Test case 18 failed! Routine: read_name_logical"
    end if
    deallocate(str_in)

    ! Test case 19 (read_name_logical): Check that commented lines are skipped
    str_in = "param_logical_1"
    bool_out = .false.
    bool_test = .true.
    bool_default = .false.
    unit_in = 10
    unit_out = 20
    call get_environment_variable("PWD", path)
    relpath = trim(path)//"/../test/test_data/"
    open(unit=unit_in, file=trim(relpath)//"test_read_namelist.txt", status="old")
    bool_out = read_name_logical(str_in, bool_default, unit_in, unit_out)
    close(unit_in)
    if (bool_out .neqv. bool_test) then
        ok = .false.
        print "(a)", "Test case 19 failed! Routine: read_name_logical"
    end if
    deallocate(str_in)

    ! Test case 20 (read_name_logical): Check that right value of named variable is returned
    str_in = "param_logical_3"
    bool_out = .false.
    bool_test = .true.
    bool_default = .false.
    unit_in = 10
    unit_out = 20
    call get_environment_variable("PWD", path)
    relpath = trim(path)//"/../test/test_data/"
    open(unit=unit_in, file=trim(relpath)//"test_read_namelist.txt", status="old")
    bool_out = read_name_logical(str_in, bool_default, unit_in, unit_out)
    close(unit_in)
    if (bool_out .neqv. bool_test) then
        ok = .false.
        print "(a)", "Test case 20 failed! Routine: read_name_logical"
    end if
    deallocate(str_in)

    ! Return test results summary
    if (ok) then
        print "(a)", "test_utility_functions: All tests passed."
    else
        print "(a)", "test_utility_functions: One ore more tests failed."
        stop 1
    end if

end program test_namefile_routines