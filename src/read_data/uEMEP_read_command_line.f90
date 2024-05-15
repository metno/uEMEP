module read_command_line
    !! Routines for interacting with the command line

    use uemep_configuration, only: name_config_file, config_date_str, n_config_files, n_max_config_files
    use uEMEP_definitions, only: model_version_str
    use uemep_logger

    implicit none
    private

    public :: uEMEP_read_command_line, check_command_line

contains

    subroutine check_command_line()
        !! Checks that a suitable number of command line arguments has been supplied
        !! and handles some special cases of command line inputs
        !!
        !! 'check_command_line() will write to stdout instead of the log file
        !! to give direct feedback to the user

        ! Local variables
        integer :: n_args, i
        character(len=256) :: arg

        ! Check if no arguments are found
        n_args = command_argument_count()
        if (n_args <= 0) then
            write(*,"(a)") " Insufficient number of command line arguments"
            write(*,"(a)") " Try 'uemep --help' for more information"
            stop
        end if

        ! Loop over input arguments to handle special cases
        do i = 1, n_args
            call get_command_argument(i, arg)
            select case(adjustl(arg))
            case("--help")
                call print_help_page()
                stop
            case("--version")
                call print_version()
                stop
            end select
        end do

        ! After checking that no special cases are found, check if the number of 
        ! arguments are within acceptable bounds (2:n_max_config_files+1)
        if (n_args < 2) then
            write(*,"(a)") " Insufficient number of command line arguments"
            write(*,"(a)") " Try 'uemep --help' for more information"
            stop
        else if (n_args > n_max_config_files + 1) then
            write(*,"(a)") " Too many command line arguments"
            write(*,"(a)") " Try 'uemep --help' for more information"
            stop
        end if
    end subroutine check_command_line

    subroutine uEMEP_read_command_line()
        !! Assigns the configuration file name(s) and substitution date_str from the command line

        ! Local variables
        integer :: n_args, i

        ! Get number of command line arguments
        n_args = command_argument_count()
        n_config_files = n_args - 1

        ! Read file names
        do i = 1, n_config_files
            call get_command_argument(i, name_config_file(i))
            write(log_msg,"(a,i0,2a)"), "name_config_file(", i, ") = ", trim(name_config_file(i))
            call log_message(log_msg, INFO)
        end do

        ! Read date string
        call get_command_argument(n_args, config_date_str)
        write(log_msg,"(2a)") "config_date_str = ", trim(config_date_str)
        call log_message(log_msg, INFO)
    end subroutine uEMEP_read_command_line

    subroutine print_help_page()
        !! Prints a help message to the terminal console
        write(*,"(a)") " uEMEP: Air quality dispersion model for high resolution downscaling of EMEP MSC-W"
        write(*,"(a)") " "
        write(*,"(a)") " To run the uEMEP model:"
        write(*,"(a)") " "
        write(*,"(a)") " ./uemep config_file_1 config_file_2 ... config_file_10 date_string"
        write(*,"(a)") " "
        write(*,"(a)") " Where config_file_n is the name of the configuration file(s) which specify"
        write(*,"(a)") " the model calculations and with a maximum number of 10 config files, and"
        write(*,"(a)") " the date_string, required, takes the form 'yyyymmdd'."
        write(*,"(a)") " "
        write(*,"(a)") " Please read the manual on how to configure uEMEP at: <URL>."
        write(*,"(a)") " "
        write(*,"(a)") " --help                     displays this help message and exits"
        write(*,"(a)") " --version                  outputs the uEMEP version and exits"
    end subroutine print_help_page

    subroutine print_version()
        !! Prints the program version to the terminal console
        write(*,"(a)") " uEMEP: Air quality dispersion model for high resolution downscaling of EMEP MSC-W"
        write(*,"(a)") " "
        write(*,"(2a)") " Version: ", trim(model_version_str)
        write(*,"(a)") " Copyright (C) 2007 Free Software Foundation."
        write(*,"(a)") " License GNU LGPL-3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>."
        write(*,"(a)") " This is free software: you are free to change and redistribute it."
        write(*,"(a)") " "
        write(*,"(a)") " Developed and maintained at the Norwegian Meteorological Institute."
        write(*,"(a)") " Contribute at: <https://github.com/metno/uEMEP>"
    end subroutine print_version

end module read_command_line

