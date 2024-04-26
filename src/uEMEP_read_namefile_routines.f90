module read_namefile_routines
    !! Subroutines and functions for reading in name files

    use uemep_constants, only: dp

    implicit none
    private

    public :: read_name_char, read_name_logical, read_name_integer, &
        read_name_real, read_name_double

    interface read_name
        module procedure read_name_real
        module procedure read_name_double
        module procedure read_name_integer
        module procedure read_name_char
        module procedure read_name_logical
    end interface read_name

contains

    function read_name_real(name_str, default_val, unit_in, unit_out) result(res)
        !! Reads single precision value from name file
        !!
        !! Skips comments (!) and returns default value if name is not present in name file
        character(len=*), intent(in) :: name_str !! Name of value in name file
        real, intent(in) :: default_val !! Default value
        integer, intent(in) :: unit_in !! Name file unit
        integer, intent(in) :: unit_out !! Log file unit
        real :: res !! Value found in name file (or default value)

        ! Local variables
        integer :: i, index_val, io
        character(len=256) :: temp_str, temp_str1, temp_str2

        ! Initially, set to default value
        res = default_val

        rewind(unit_in)
        do
            read(unit_in, "(a)", iostat=io) temp_str
            if (io /= 0) then
                exit
            end if

            ! Remove tabs
            index_val = 0
            do i = 1, len(temp_str)
                if (ichar(temp_str(i:i)) .ne. 9) then
                    index_val = index_val + 1
                    temp_str1(index_val:index_val) = temp_str(i:i)
                end if
            end do
            temp_str = ADJUSTL(temp_str1)
            temp_str1 = ''

            ! If not a comment
            if (trim(temp_str(1:1)) .ne. '!') then
                ! Find the position of the equals sign if there is one
                index_val = index(temp_str, '=', back=.false.)
                if (index_val .gt. 1) then
                    ! Create the string before the equals sign
                    temp_str1 = trim(temp_str(1:index_val-1))

                    ! Check to see if it is a matching string
                    if (trim(temp_str1) .eq. trim(name_str)) then
                        ! Create the string after the equals sign
                        temp_str2=temp_str(index_val+1:)
                        temp_str2=adjustl(temp_str2)
                        if (len(trim(temp_str2)) .ge. 1) then
                            read(temp_str2,*, iostat=io) res
                            if (io /= 0) then
                                cycle
                            end if
                            write(unit_out,'(A,es12.4)') 'Setting: '//trim(name_str)//' = ', res
                        end if
                    end if
                end if
            end if
        end do
    end function read_name_real

    function read_name_double(name_str, default_val, unit_in, unit_out) result(res)
        !! Reads double precision value from name file
        !!
        !! Skips comments (!) and returns default value if name is not present in name file
        character(len=*), intent(in) :: name_str !! Name of value in name file 
        real(dp), intent(in) :: default_val !! Default value
        integer, intent(in) :: unit_in !! Name file unit
        integer, intent(in) :: unit_out !! Log file unit
        real(dp) :: res !! Value found in name file (or default value)

        ! Local variables
        integer :: i, index_val, io
        character(len=256) :: temp_str, temp_str1, temp_str2

        ! Initiallly set to default value
        res=default_val

        rewind(unit_in)
        !do while (.not.eof(unit_in))
        do
            read(unit_in,'(A)', iostat=io) temp_str
            if (io /= 0) then
                exit
            end if

            ! Remove tabs
            index_val = 0
            do i = 1, len(temp_str)
                if (ichar(temp_str(i:i)) .ne. 9) then
                    index_val = index_val + 1
                    temp_str1(index_val:index_val) = temp_str(i:i)
                end if
            end do
            temp_str = ADJUSTL(temp_str1)
            temp_str1 = ''

            ! If not a comment
            if (trim(temp_str(1:1)) .ne. '!') then
                ! Find the position of the equals sign if there is one
                index_val = index(temp_str, '=', back=.false.)
                if (index_val.gt.1) then
                    ! Create the string before the equals sign
                    temp_str1 = trim(temp_str(1:index_val-1))

                    ! Check to see if it is a matching string
                    if (trim(temp_str1) .eq. trim(name_str)) then
                        ! Create the string after the equals sign
                        temp_str2 = temp_str(index_val+1:)
                        temp_str2 = adjustl(temp_str2)
                        if (len(trim(temp_str2)) .ge. 1) then
                            read(temp_str2,*, iostat=io) res
                            if (io /= 0) then
                                cycle
                            end if
                            write(unit_out,'(A,es12.4)') 'Setting: '//trim(name_str)//' = ', res
                        end if
                    end if
                end if
            end if
        end do
    end function read_name_double

    function read_name_integer(name_str, default_val, unit_in, unit_out) result(res)
        !! Reads integer value from name file
        !!
        !! Skips comments (!) and returns default value if name is not present in name file
        character(len=*), intent(in) :: name_str !! Name of value in name file
        integer, intent(in) :: default_val !! Default value
        integer, intent(in) :: unit_in !! Name file unit
        integer, intent(in) :: unit_out !! Log file unit
        integer :: res !! Value found in name file (or default value)

        ! Local variables
        integer :: i, index_val, io
        character(256) temp_str,temp_str1,temp_str2

        ! Initially set return value to default value
        res = default_val

        rewind(unit_in)
        do
            read(unit_in,'(A)', iostat=io) temp_str
            if (io /= 0) then
                exit
            end if

            ! Remove tabs
            index_val = 0
            do i = 1, len(temp_str)
                if (ichar(temp_str(i:i)) .ne. 9) then
                    index_val = index_val + 1
                    temp_str1(index_val:index_val) = temp_str(i:i)
                end if
            end do
            temp_str = ADJUSTL(temp_str1)
            temp_str1 = ''

            !If not a comment
            if (trim(temp_str(1:1)) .ne. '!') then
                ! Find the position of the equals sign if there is one
                index_val = index(temp_str, '=', back=.false.)
                if (index_val.gt.1) then
                    ! Create the string before the equals sign
                    temp_str1 = trim(temp_str(1:index_val-1))

                    ! Check to see if it is a matching string
                    if (trim(temp_str1) .eq. trim(name_str)) then
                        ! Create the string after the equals sign
                        temp_str2=temp_str(index_val+1:)
                        temp_str2=adjustl(temp_str2)
                        if (len(trim(temp_str2)) .ge. 1) then
                            read(temp_str2,*, iostat=io) res
                            if (io /= 0) then
                                cycle
                            end if
                            write(unit_out,'(A,i12)') 'Setting: '//trim(name_str)//' = ', res
                        end if
                    end if
                end if
            end if
        end do
    end function read_name_integer

    function read_name_char(name_str, default_val, unit_in, unit_out) result(res)
        !! Reads string from name file
        !!
        !! Skips comments (!) and returns default string if name is not present in name file
        character(len=*), intent(in) :: name_str !! Name of string in name file
        character(len=*), intent(in) :: default_val !! Default string
        integer, intent(in) :: unit_in !! Name file unit
        integer, intent(in) :: unit_out !! Log file unit
        character(len=:), allocatable :: res !! String found in name file (or default string)

        ! Local variables
        integer :: i, index_val, io
        character(256) :: temp_str, temp_str1, temp_str2
        character(256) :: call_str = 'read_name_char'

        ! Initially set default string as return string
        res = default_val

        rewind(unit_in)
        do
            read(unit_in,'(A)',iostat=io) temp_str
            if (io /= 0) then
                exit
            end if

            ! Remove tabs
            index_val = 0
            do i = 1, len(temp_str)
                if (ichar(temp_str(i:i)) .ne. 9) then
                    index_val = index_val + 1
                    temp_str1(index_val:index_val) = temp_str(i:i)
                end if
            end do
            temp_str = ADJUSTL(temp_str1)
            temp_str1 = ''

            !If not a comment
            if (trim(temp_str(1:1)) .ne. '!') then
                !Find the position of the equals sign if there is one
                index_val = index(temp_str, '=', back=.false.)
                if (index_val .gt. 1) then
                    ! Create the string before the equals sign
                    temp_str1 = trim(temp_str(1:index_val-1))

                    ! Check to see if it is a matching string
                    if (trim(temp_str1) .eq. trim(name_str)) then
                        !Create the string after the equals sign
                        temp_str2 = temp_str(index_val+1:)
                        temp_str2 = adjustl(temp_str2)
                        if (len(trim(temp_str2)) .ge. 1) then
                            ! Special for characters so it doesn't read it's own call
                            if (trim(temp_str2(1:min(len(trim(temp_str2)),len(trim(call_str))))) .ne. trim(call_str)) then
                                read(temp_str2,*, iostat=io) res
                                if (io /= 0) then
                                    cycle
                                end if
                                res = adjustl(res)
                                write(unit_out,'(A,A)') 'Setting: '//trim(name_str)//' = ', trim(res)
                            end if
                        end if
                    end if
                end if
            end if
        end do
    end function read_name_char

    function read_name_logical(name_str, default_val, unit_in, unit_out) result(res)
        !! Reads boolean value from name file
        !!
        !! Skips comments (!) and returns default value if name is not present in name file
        character(len=*), intent(in) :: name_str !! Name of value in name file
        logical, intent(in) :: default_val !! Default value
        integer, intent(in) :: unit_in !! Name file unit
        integer, intent(in) :: unit_out !! Log file unit
        logical :: res

        ! Local variables
        integer :: i, index_val, io
        character(len=256) :: temp_str,temp_str1,temp_str2

        ! Initially set return value as default
        res=default_val

        rewind(unit_in)
        do
            read(unit_in,'(A)', iostat=io) temp_str
            if (io /= 0) then
                exit
            end if

            ! Remove tabs
            index_val = 0
            do i = 1, len(temp_str)
                if (ichar(temp_str(i:i)) .ne. 9) then
                    index_val = index_val + 1
                    temp_str1(index_val:index_val) = temp_str(i:i)
                end if
            end do
            temp_str = ADJUSTL(temp_str1)
            temp_str1 = ''

            ! If not a comment
            if (trim(temp_str(1:1)) .ne. '!') then
                ! Find the position of the equals sign if there is one
                index_val = index(temp_str, '=', back=.false.)
                if (index_val .gt. 1) then
                    ! Create the string before the equals sign
                    temp_str1 = trim(temp_str(1:index_val-1))

                    ! Check to see if it is a matching string
                    if (trim(temp_str1) .eq. trim(name_str)) then
                        ! Create the string after the equals sign
                        temp_str2 = temp_str(index_val+1:)
                        temp_str2 = adjustl(temp_str2)
                        if (len(trim(temp_str2)) .ge. 1) then
                            read(temp_str2,*, iostat=io) res
                            if (io /= 0) then
                                cycle
                            end if
                            write(unit_out,'(A,L)') 'Setting: '//trim(name_str)//' = ', res
                        end if
                    end if
                end if
            end if
        end do
    end function read_name_logical

end module read_namefile_routines

