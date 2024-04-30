module uemep_configuration

    use uemep_logger

    implicit none
    private

    public :: initialize_uemep

    ! Model version
    character(len=256) :: model_version_str

    ! Command line config
    integer, parameter :: n_max_config_files = 10
    integer, protected :: n_config_files
    character(len=256) :: name_config_file(n_max_config_files) = ""
    character(len=256), parameter :: filename_log_file = "uEMEP_log.txt"
    character(len=256) :: pathname_log_file = ""
    character(len=256) :: config_date_str = ""
    character(len=256) :: emission_date_str = ""

    ! Local module variables
    character(len=1028) :: log_msg

contains

    subroutine initialize_uemep()

        call read_command_line()

    end subroutine initialize_uemep

    subroutine read_command_line()
        !! Reads the configuration file name and substitution date_str from the command line
        !!
        !! Can have up to 10 config files, limitted by the array definition
        !! Last command line string is always the date string

        ! Local variables
        integer :: nb, i

        ! Get number of command line arguments
        nb = command_argument_count()

        if (nb > n_max_config_files + 1) then
            write(log_msg,"(a,i0,a)") "Too many command line arguments. Maximum is ", n_max_config_files, &
                " configuration files plus one date_str. Stopping uEMEP"
            call log_message(log_msg, ERROR)
            call close_log_file()
            stop 1
        end if

        if (nb >= 2) then
            n_config_files = nb + 1
            ! Get config file names
            do i = 1, n_config_files
                call get_command_argument(i, name_config_file(i))
                write(log_msg,"(a,i0,2a)") "name_config_file(", i, ") = ", trim(name_config_file(i))
                call log_message(log_msg, INFO)
            end do
            ! Get date string
            call get_command_argument(nb, config_date_str)
            write(log_msg,"(2a)") "config_date_str = ", trim(config_date_str)
            call log_message(log_msg, INFO)
        else
            write(log_msg,"(a)") "Insufficient number of command line arguments. Stopping uEMEP"
            call log_message(log_msg, ERROR)
            call close_log_file()
            stop 1
        end if
    end subroutine read_command_line

end module uemep_configuration