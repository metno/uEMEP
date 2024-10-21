module io_functions
    !! This module contains procedures that support IO operations

    use uEMEP_definitions, only: unit_logfile

    implicit none
    private

    public :: check_dir_exist

contains

    function check_dir_exist(path) result(exists)
        !! Checks if a directory exists and/or if we have write permission
        !! in that directory.
        !! 
        !! The specific 'path' must end with '/'
        !!
        !! Returns .true. if directory exists and is accessible
        character(len=*), intent(in) :: path !! Directory path
        logical:: exists 

        ! Local variables
        integer :: u, ios

        ! Attempt to open a temporary file in the specified directory
        open(newunit=u, file=path//"tmp", status="new", iostat=ios)
        if (ios /= 0) then
            exists = .false.
        else
            ! If succesfull, delete the temporary file
            close(unit=u, status="delete")
            exists = .true.
        end if
    end function check_dir_exist

end module io_functions