module plugin_base

    implicit none
    private

    public :: plugin_t

    type, abstract :: plugin_t
        character(len=32) :: source_name
        character(len=32) :: compound_names
    contains
        procedure(init_if), deferred :: init
        procedure(update_if), deferred :: update
        procedure(finalize_if), deferred :: finalize
    end type plugin_t

    abstract interface
        subroutine init_if(self)
            import :: plugin_t
            class(plugin_t), intent(inout) :: self
        end subroutine init_if
        subroutine update_if(self)
            import :: plugin_t
            class(plugin_t), intent(inout) :: self
        end subroutine update_if
        subroutine finalize_if(self)
            import :: plugin_t
            class(plugin_t), intent(inout) :: self
        end subroutine finalize_if
    end interface

end module plugin_base