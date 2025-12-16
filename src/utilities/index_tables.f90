module index_tables

    implicit none
    private
    
    public :: index_table_t

    integer, parameter :: default_key_len = 128
    integer, parameter :: default_key_value_pairs = 128

    type :: index_table_t
        integer :: n = 0
        integer :: cap = 0
        character(len=default_key_len), allocatable :: keys(:)
        integer, allocatable :: vals(:)
    contains
        procedure :: init => index_table_init
        procedure :: add => index_table_add
        procedure :: get => index_table_get
    end type index_table_t

contains

    subroutine index_table_init(self, hint)
        class(index_table_t), intent(inout) :: self
        integer, intent(in), optional :: hint
        integer :: c

        if (present(hint)) then
            c = hint
        else
            c = default_key_value_pairs
        end if

        if (allocated(self%keys)) deallocate(self%keys)
        if (allocated(Self%vals)) deallocate(self%vals)
        allocate(self%keys(c))
        allocate(self%vals(c))
        self%keys = ""
        self%vals = 0
        self%n = 0
        self%cap = c
    end subroutine index_table_init

    subroutine index_table_add(self, key)
        class(index_table_t), intent(inout) :: self
        character(len=*), intent(in) :: key
        integer :: idx
        
        if (.not. allocated(self%keys)) call self%init()
        
        idx = find_index(self, key)
        if (idx /= 0) return

        if (self%n + 1 <= self%cap) then
            self%n = self%n + 1
            self%keys(self%n) = key
            self%vals(self%n) = self%n
        else
            print *, "ERROR: hash table full"
            stop 1
        end if
    end subroutine index_table_add

    function index_table_get(self, key, found) result(val)
        class(index_table_t), intent(in) :: self
        character(len=*), intent(in) :: key
        logical, intent(out), optional :: found
        integer :: val
        integer :: idx

        if (.not. allocated(self%keys)) then
            if (present(found)) found = .false.
            val = -1
            return
        end if

        idx = find_index(self, key)
        if (idx == 0) then
            if (present(found)) found = .false.
            val = -1
        else
            if (present(found)) found = .true.
            val = self%vals(idx)
        end if
    end function index_table_get

    function find_index(self, key) result(idx)
        class(index_table_t), intent(in) :: self
        character(len=*), intent(in) :: key
        integer :: idx, i

        idx = 0
        if (.not. allocated(self%keys)) return
        do i = 1, self%n
            if (trim(self%keys(i)) == trim(key)) then
                idx = i
                return
            end if
        end do
    end function find_index

end module index_tables