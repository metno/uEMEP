module enable_plugins

    implicit none
    private

    public :: init_plugins

contains

    subroutine init_plugins()
        use plugin_birch, only: register_plugin_birch

        call register_plugin_birch()
    end subroutine init_plugins

end module enable_plugins