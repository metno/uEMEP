module mod_agriculture

    use uemep_definitions, only: unit_logfile, x_dim_index, y_dim_index, emission_max_subgrid_dim, agriculture_index, &
        proxy_emission_subgrid, emission_subgrid_dim
    use uemep_configuration, only: filename_agriculture, pathname_agriculture, limit_agriculture_delta, agriculture_var_name
    use mod_emission_utils

    implicit none
    private

    public :: initialize_agriculture

    character(len=11), parameter :: sector_name = "agriculture"
    real, allocatable :: agriculture_subgrid(:,:)
    integer, allocatable :: crossref_emission_to_agriculture_subgrid(:,:,:)
    integer :: agriculture_subgrid_dim(2)
    real :: agriculture_subgrid_min(2)
    real :: agriculture_subgrid_max(2)
    real :: agriculture_subgrid_delta(2)
    real, allocatable :: x_agriculture_subgrid(:,:)
    real, allocatable :: y_agriculture_subgrid(:,:)
    integer :: agriculture_buffer_index(2)
    real :: agriculture_buffer_size(2)

contains

    subroutine initialize_agriculture()
        write(unit_logfile, "(a)") ""
        write(unit_logfile, "(a)") "================================================================"
        write(unit_logfile, "(a)") "Initializing agriculture emissions (GNFR sector 12 (K))"
        write(unit_logfile, "(a)") "================================================================"

        call setup_agriculture_arrays_and_variables()
        call read_agriculture_data()
        call redistribute_agriculture_emissions()
    end subroutine initialize_agriculture

    subroutine setup_agriculture_arrays_and_variables()

        call setup_subgrid_dimensions(sector_name, agriculture_subgrid_delta, agriculture_subgrid_min, &
            agriculture_subgrid_max, agriculture_subgrid_dim, limit_agriculture_delta)
        call setup_buffer_zone(sector_name, agriculture_buffer_index, agriculture_buffer_size, &
            agriculture_subgrid_delta, agriculture_subgrid_min, agriculture_subgrid_max, agriculture_subgrid_dim)
        
        if (allocated(agriculture_subgrid)) deallocate(agriculture_subgrid)
        allocate(agriculture_subgrid(agriculture_subgrid_dim(x_dim_index), agriculture_subgrid_dim(y_dim_index)))
        agriculture_subgrid = 0.0
        if (allocated(x_agriculture_subgrid)) deallocate(x_agriculture_subgrid)
        allocate(x_agriculture_subgrid(agriculture_subgrid_dim(x_dim_index), agriculture_subgrid_dim(y_dim_index)))
        x_agriculture_subgrid = 0.0
        if (allocated(y_agriculture_subgrid)) deallocate(y_agriculture_subgrid)
        allocate(y_agriculture_subgrid(agriculture_subgrid_dim(x_dim_index), agriculture_subgrid_dim(y_dim_index)))
        y_agriculture_subgrid = 0.0
        if (allocated(crossref_emission_to_agriculture_subgrid)) deallocate(crossref_emission_to_agriculture_subgrid)
        allocate(crossref_emission_to_agriculture_subgrid(emission_max_subgrid_dim(x_dim_index), &
            emission_max_subgrid_dim(y_dim_index), 2))
        
        call set_subgrid_xy(sector_name, agriculture_subgrid_dim, agriculture_subgrid_min, &
            agriculture_subgrid_delta, x_agriculture_subgrid, y_agriculture_subgrid)
        call setup_crossref_grid(sector_name, crossref_emission_to_agriculture_subgrid, &
            agriculture_subgrid_delta, agriculture_subgrid_min, agriculture_subgrid_dim)
    end subroutine setup_agriculture_arrays_and_variables

    subroutine read_agriculture_data()
        call read_netcdf_data(sector_name, pathname_agriculture, filename_agriculture, agriculture_subgrid, &
            agriculture_subgrid_min, agriculture_subgrid_max, agriculture_subgrid_delta, agriculture_subgrid_dim, &
            x_agriculture_subgrid, y_agriculture_subgrid, agriculture_var_name)
    end subroutine read_agriculture_data

    subroutine redistribute_agriculture_emissions()
        integer :: i, j, i_agriculture, j_agriculture
        proxy_emission_subgrid(:,:,agriculture_index,:) = 0.0
        do j = 1, emission_subgrid_dim(y_dim_index, agriculture_index)
            do i = 1, emission_subgrid_dim(x_dim_index, agriculture_index)
                i_agriculture = crossref_emission_to_agriculture_subgrid(i,j,x_dim_index)
                j_agriculture = crossref_emission_to_agriculture_subgrid(i,j,y_dim_index)
                proxy_emission_subgrid(i,j,agriculture_index,:) = agriculture_subgrid(i_agriculture,j_agriculture)
            end do
        end do
    end subroutine redistribute_agriculture_emissions

end module mod_agriculture