module uemep_constants

    implicit none
    private

    integer, parameter, public :: dp = selected_real_kind(15, 307)

    real, parameter, public :: secphour = 3600.0
    real, parameter, public :: secpday = 86400.0

end module uemep_constants
