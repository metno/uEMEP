program test_uemep_logger

    use uemep_logger

    implicit none

    logical :: ok = .true.

    character(len=32) :: logfile_name = "test_logfile.txt"
    logical :: exists, isopen
    integer :: io, unit_in
    character(len=1028) :: msg_in, msg_log, msg_test

    ! Test case 1 (open_log_file): Test that a log file is opened
    call open_log_file(trim(logfile_name), io_err=io)
    inquire(file=trim(logfile_name), exist=exists, opened=isopen)
    if (.not. exists .or. isopen .eqv. .false. .or. io /= 0) then
        ok = .false.
        print "(a)", "Test case 1 failed! Routine: open_log_file"
    end if
    call close_log_file()
    call execute_command_line("rm " // logfile_name)

    ! Test case 2 (close_log_file): Check that the log file is closed properly
    call open_log_file(trim(logfile_name))
    call close_log_file(io_err=io)
    inquire(file=trim(logfile_name), opened=isopen)
    if (io /= 0 .or. isopen .eqv. .true.) then
        ok = .false.
        print "(a)", "Test case 2 failed! Routine: close_log_file"
    end if
    call execute_command_line("rm " // logfile_name)

    ! Test case 3 (log_message): Check output against expected values using INFO
    msg_in = "This is an INFO message"
    msg_log = "-999"
    msg_test = " [INFO]  This is an INFO message"
    call open_log_file(trim(logfile_name))
    call set_log_level(INFO)
    call log_message(msg_in, INFO)
    call close_log_file()
    open(newunit=unit_in, file=logfile_name, status="old")
    read(unit_in,"(a)") msg_log
    close(unit_in)
    if (msg_log(18:) /= trim(msg_test)) then ! have to remove the timestamp
        ok = .false.
        print "(a)", "Test case 3 failed! Routine: log_message"
    end if
    call execute_command_line("rm " // logfile_name)

    ! Test case 4 (log_message): Check output against expected values using WARNING
    msg_in = "This is a WARNING message"
    msg_log = "-999"
    msg_test = " [WARN]  This is a WARNING message"
    call open_log_file(trim(logfile_name))
    call set_log_level(INFO)
    call log_message(msg_in, WARNING)
    call close_log_file()
    open(newunit=unit_in, file=logfile_name, status="old")
    read(unit_in,"(a)") msg_log
    close(unit_in)
    if (msg_log(18:) /= trim(msg_test)) then ! have to remove the timestamp
        ok = .false.
        print "(a)", "Test case 4 failed! Routine: log_message"
    end if
    call execute_command_line("rm " // logfile_name)

    ! Test case 5 (log_message): Check output against expected values using ERROR
    msg_in = "This is an ERROR message"
    msg_log = "-999"
    msg_test = " [ERROR] This is an ERROR message"
    call open_log_file(trim(logfile_name))
    call set_log_level(INFO)
    call log_message(msg_in, ERROR)
    call close_log_file()
    open(newunit=unit_in, file=logfile_name, status="old")
    read(unit_in,"(a)") msg_log
    close(unit_in)
    if (msg_log(18:) /= trim(msg_test)) then ! have to remove the timestamp
        ok = .false.
        print "(a)", "Test case 5 failed! Routine: log_message"
    end if
    call execute_command_line("rm " // logfile_name)

    ! Test case 6 (log_message): Check that DEBUG is NOT written with log level INFO (i.e., > DEBUG)
    msg_in = "This is a DEBUG message"
    msg_log = "-999"
    msg_test = " [DEBUG] This is a DEBUG message"
    call open_log_file(trim(logfile_name))
    call set_log_level(INFO)
    call log_message(msg_in, DEBUG)
    call close_log_file()
    open(newunit=unit_in, file=logfile_name, status="old")
    read(unit_in,"(a)", iostat=io) msg_log
    if (io == 0 .or. msg_log(18:) == trim(msg_test)) then
        ok = .false.
        print "(a)", "Test case 6 failed! Routine: log_message"
    end if
    call execute_command_line("rm " // logfile_name)

    ! Test case 7 (log_message): Check that DEBUG is written with log level DEBUG
    msg_in = "This is a DEBUG message"
    msg_log = "-999"
    msg_test = " [DEBUG] This is a DEBUG message"
    call open_log_file(trim(logfile_name))
    call set_log_level(DEBUG)
    call log_message(msg_in, DEBUG)
    call close_log_file()
    open(newunit=unit_in, file=logfile_name, status="old")
    read(unit_in,"(a)", iostat=io) msg_log
    if (io /= 0 .or. msg_log(18:) /= trim(msg_test)) then
        ok = .false.
        print "(a)", "Test case 7 failed! Routine: log_message"
    end if
    call execute_command_line("rm " // logfile_name)

    ! Return test results summary
    if (ok) then
        print "(a)", "test_utility_functions: All tests passed."
    else
        print "(a)", "test_utility_functions: One ore more tests failed."
        stop 1
    end if

end program test_uemep_logger