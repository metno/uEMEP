module plugin_birch

    use plugin_base, only: plugin_t
    use plugin_registry, only: register_plugin_instance

    implicit none
    private

    public :: register_plugin_birch

    type, extends(plugin_t) :: plugin_birch_t
    contains
        procedure :: init => plugin_init
        procedure :: update => plugin_update
        procedure :: finalize => plugin_finalize
    end type plugin_birch_t

    class(plugin_birch_t), pointer :: my_plugin => null()

contains

    subroutine register_plugin_birch()
        class(plugin_t), pointer :: base
        if (.not. associated(my_plugin)) allocate(plugin_birch_t :: my_plugin)
        base => my_plugin
        call register_plugin_instance(base)
    end subroutine register_plugin_birch

    subroutine plugin_init(self)
        class(plugin_birch_t), intent(inout) :: self

        self%source_name = "birch"
    end subroutine plugin_init

    subroutine plugin_update(self)
        class(plugin_birch_t), intent(inout) :: self
    
    end subroutine plugin_update

    subroutine plugin_finalize(self)
        class(plugin_birch_t), intent(inout) :: self

    end subroutine plugin_finalize

end module plugin_birch