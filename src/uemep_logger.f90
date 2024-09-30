module uemep_logger
    ! A simple logger module for uEMEP

    implicit none
    private

    public :: DEBUG, INFO, WARNING, ERROR
    public :: open_log_file, close_log_file, set_log_level, log_message, log_header
    public :: log_msg

    integer, parameter :: DEBUG = 1
    integer, parameter :: INFO = 2
    integer, parameter :: WARNING = 3
    integer, parameter :: ERROR = 4
    
    integer, save :: log_level = INFO
    character(len=256), save :: log_name
    logical, save :: file_opened = .false.
    integer :: log_unit
    character(len=1028) :: log_msg

    integer :: unit

contains

    subroutine open_log_file(logfile_name, io_err)
        !! Opens a new log file for writing
        character(len=*), intent(in) :: logfile_name
        integer, intent(inout), optional :: io_err

        ! Store unit and file name information in module
        log_name = trim(logfile_name)

        ! Open log file
        if (present(io_err)) then
            open(newunit=log_unit, file=log_name, status="replace", iostat=io_err)
        else
            open(newunit=log_unit, file=log_name, status="replace")
        end if
        file_opened = .true.
    end subroutine open_log_file

    subroutine close_log_file(io_err)
        integer, intent(inout), optional :: io_err
        !! Closes the log file
        if (file_opened) then
            if (present(io_err)) then
                close(log_unit, iostat=io_err)
            else
                close(log_unit)
            end if
            log_unit = -1
        end if
    end subroutine close_log_file

    subroutine set_log_level(level)
        !! Sets the log level
        !!
        !! The logger will write all message at or above this level
        integer, intent(in) :: level
        if (level >= DEBUG .and. level <= ERROR) then
            log_level = level
        else
            print *, "ERROR! Invalid log level: ", level
            stop 1
        end if
    end subroutine set_log_level

    subroutine log_header(message, level, upper_space, lower_space)
        !! Send a log header message to the log file
        !!
        !! Log levels: DEBUG, INFO, WARNING and ERROR        
        character(len=*), intent(in) :: message !! Message to log
        integer, intent(in) :: level !! Log level
        logical, intent(in), optional :: upper_space, lower_space

        ! Local variables
        logical :: u_space, l_space

        if (present(upper_space)) then
            u_space = upper_space
        else
            u_space = .true.
        end if

        if (present(lower_space)) then
            l_space = lower_space
        else
            l_space = .true.
        end if
        
        if (u_space) call log_message("", level)
        call log_message("================================================================", level)
        call log_message(message, level)
        call log_message("================================================================", level)
        if (l_space) call log_message("", level)
    end subroutine log_header

    subroutine log_message(message, level)
        !! Send a log message to the log file
        !!
        !! Log levels: DEBUG, INFO, WARNING and ERROR
        character(len=*), intent(in) :: message !! Message to log
        integer, intent(in) :: level !! Log level
        
        ! Local variables
        character(len=18) :: timestamp
        character(len=8) :: datestr
        character(len=10) :: timestr

        ! Get current date and time
        call date_and_time(date=datestr, time=timestr)
        timestamp = '[' // datestr // ' ' // timestr(1:6) // '] '

        if (level >= log_level) then
            select case(level)
            case(DEBUG)
                call write_log(timestamp // '[DEBUG] ' // message)
            case(INFO)
                call write_log(timestamp // '[INFO]  ' // message)
            case(WARNING)
                call write_log(timestamp // '[WARN]  ' // message)
            case(ERROR)
                call write_log(timestamp // '[ERROR] ' // message)
            case default
                print "(a)", "ERROR: Invalid log level: ", level
                stop 1
            end select
        end if
    end subroutine log_message

    subroutine write_log(message)
        !! Writes the message to the log file
        character(len=*), intent(in) :: message !! Message to log

        if (.not. file_opened) then
            print "(a)", "ERROR: No file has been opened for writing"
            stop 1
        end if

        write(log_unit, "(a)") trim(adjustl(message))
    end subroutine write_log

end module uemep_logger