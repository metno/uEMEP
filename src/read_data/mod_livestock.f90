module mod_livestock

    use uEMEP_definitions, only: unit_logfile, x_dim_index, y_dim_index, emission_max_subgrid_dim, livestock_index, &
        proxy_emission_subgrid, emission_subgrid_dim
    use uemep_configuration, only: filename_livestock, pathname_livestock, limit_livestock_delta, livestock_var_name
    use mod_emission_utils

    implicit none
    private

    public :: initialize_livestock

    character(len=9), parameter :: sector_name = "livestock"
    real, allocatable :: livestock_subgrid(:,:) ! Array to read in livestock (proxy) data
    integer, allocatable :: crossref_emission_to_livestock_subgrid(:,:,:)
    integer :: livestock_subgrid_dim(2)
    real :: livestock_subgrid_min(2)
    real :: livestock_subgrid_max(2)
    real :: livestock_subgrid_delta(2)
    real, allocatable :: x_livestock_subgrid(:,:)
    real, allocatable :: y_livestock_subgrid(:,:)
    integer :: livestock_buffer_index(2)
    real :: livestock_buffer_size(2)

contains

    subroutine initialize_livestock()
        write(unit_logfile, "(a)") ""
        write(unit_logfile, "(a)") "================================================================"
        write(unit_logfile, "(a)") "Initializing livestock emissions (GNFR sector 11 (K))"
        write(unit_logfile, "(a)") "================================================================"

        call setup_livestock_arrays_and_variables()
        call read_livestock_data()
        call redistribute_livestock_emissions()
    end subroutine initialize_livestock

    subroutine setup_livestock_arrays_and_variables()

        call set_subgrid_dimensions(sector_name, livestock_subgrid_delta, livestock_subgrid_min, &
            livestock_subgrid_max, livestock_subgrid_dim, limit_livestock_delta)
        call set_buffer_zone(sector_name, livestock_buffer_index, livestock_buffer_size, &
            livestock_subgrid_delta, livestock_subgrid_min, livestock_subgrid_max, livestock_subgrid_dim)
        
        if (allocated(livestock_subgrid)) deallocate(livestock_subgrid)
        allocate(livestock_subgrid(livestock_subgrid_dim(x_dim_index),livestock_subgrid_dim(y_dim_index)))
        livestock_subgrid = 0.0
        if (allocated(x_livestock_subgrid)) deallocate(x_livestock_subgrid)
        allocate(x_livestock_subgrid(livestock_subgrid_dim(x_dim_index),livestock_subgrid_dim(y_dim_index)))
        x_livestock_subgrid = 0.0
        if (allocated(y_livestock_subgrid)) deallocate(y_livestock_subgrid)
        allocate(y_livestock_subgrid(livestock_subgrid_dim(x_dim_index),livestock_subgrid_dim(y_dim_index)))
        y_livestock_subgrid = 0.0
        if (allocated(crossref_emission_to_livestock_subgrid)) deallocate(crossref_emission_to_livestock_subgrid)
        allocate(crossref_emission_to_livestock_subgrid(emission_max_subgrid_dim(x_dim_index), &
            emission_max_subgrid_dim(y_dim_index), 2))

        call set_subgrid_xy(sector_name, livestock_subgrid_dim, livestock_subgrid_min, &
            livestock_subgrid_delta, x_livestock_subgrid, y_livestock_subgrid)
        call set_crossref_grid(sector_name, livestock_index, crossref_emission_to_livestock_subgrid, &
            livestock_subgrid_delta, livestock_subgrid_min, livestock_subgrid_dim)
    end subroutine setup_livestock_arrays_and_variables

    subroutine read_livestock_data()
        call read_netcdf_data(sector_name, pathname_livestock, filename_livestock, livestock_subgrid, &
            livestock_subgrid_min, livestock_subgrid_max, livestock_subgrid_delta, livestock_subgrid_dim, &
            x_livestock_subgrid, y_livestock_subgrid, livestock_var_name)
    end subroutine read_livestock_data

    subroutine redistribute_livestock_emissions()
        integer :: i, j, i_livestock, j_livestock
        proxy_emission_subgrid(:,:,livestock_index,:) = 0.0
        do j = 1, emission_subgrid_dim(y_dim_index, livestock_index)
            do i = 1, emission_subgrid_dim(x_dim_index, livestock_index)
                i_livestock = crossref_emission_to_livestock_subgrid(i,j,x_dim_index)
                j_livestock = crossref_emission_to_livestock_subgrid(i,j,y_dim_index)
                proxy_emission_subgrid(i,j,livestock_index,:) = livestock_subgrid(i_livestock,j_livestock)
            end do
        end do
    end subroutine redistribute_livestock_emissions

end module mod_livestock