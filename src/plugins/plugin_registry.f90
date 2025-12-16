module plugin_registry

    use plugin_base, only: plugin_t

    implicit none
    private

    public :: init_registry, register_plugin_instance, run_plugins, finalize_plugins
    public :: nreg, registry

    integer, parameter :: max_plugins = 12

    type :: registered_plugin_t
        class(plugin_t), pointer :: inst => null()
    end type registered_plugin_t

    type(registered_plugin_t), allocatable :: registry(:)
    integer :: nreg = 0

contains

    subroutine init_registry()
        if (.not. allocated(registry)) then
            allocate(registry(max_plugins))
            nreg = 0
        end if
    end subroutine init_registry

    subroutine register_plugin_instance(instance)
        class(plugin_t), pointer :: instance
        call init_registry()
        if (nreg >= size(registry)) then
            write(*,"(a)") "ERROR: Too many plugins registrered. Increase 'nreg'."
            stop 1
        end if
        nreg = nreg + 1
        registry(nreg)%inst => instance
        call instance%init()
    end subroutine register_plugin_instance

    subroutine run_plugins()
        integer :: i
        do i = 1, nreg
            if (associated(registry(i)%inst)) call registry(i)%inst%update()
        end do
    end subroutine run_plugins

    subroutine finalize_plugins()
        integer :: i
        do i = 1, nreg
            if (associated(registry(i)%inst)) call registry(i)%inst%finalize()
        end do
    end subroutine finalize_plugins

end module plugin_registry