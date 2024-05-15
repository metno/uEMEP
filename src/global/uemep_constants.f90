module uemep_constants

    implicit none

    ! Type constants
    integer, parameter, public :: dp = selected_real_kind(15, 307)

    ! Mathematical constants
    real, parameter :: pi = 3.14159265358979323

    ! Time constants
    real, parameter, public :: secphour = 3600.0
    real, parameter, public :: secpday = 86400.0

    ! Other constants
    real, parameter :: NODATA_value = -128.0
    real, parameter :: epsilon0 = 1.0e-6


end module uemep_constants
