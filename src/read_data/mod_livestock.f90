module mod_livestock

    use mod_emission_utils, only: open_netcdf_file, read_netcdf_data, setup_subgrid_dimensions, setup_buffer_zone, &
        setup_crossref_grid, set_subgrid_xy
    use uemep_configuration, only: filename_livestock, pathname_livestock, local_subgrid_method_flag, &
        subgrid_min, subgrid_max, subgrid_delta, limit_livestock_delta, livestock_var_name
    use uEMEP_definitions, only: x_dim_index, y_dim_index, livestock_index, &
        emission_subgrid_dim, x_emission_subgrid, y_emission_subgrid, proxy_emission_subgrid, &
        emission_max_subgrid_dim, buffer_index_scale, subgrid_dim, use_buffer_zone
    use define_subgrid, only: dx_temp, dy_temp

    implicit none
    private

    public :: initialize_livestock

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
        use uemep_definitions, only: unit_logfile
        write(unit_logfile, "(a)") ""
        write(unit_logfile, "(a)") "================================================================"
        write(unit_logfile, "(a)") "Initializing livestock emissions (GNFR sector 11 (K))"
        write(unit_logfile, "(a)") "================================================================"

        call setup_livestock_arrays_and_variables()
        call read_livestock_data()
        call redistribute_livestock_emissions()
    end subroutine initialize_livestock

    subroutine setup_livestock_arrays_and_variables()
        ! Setup subgrid
        call setup_subgrid_dimensions("livestock", livestock_subgrid_delta, livestock_subgrid_min, &
            livestock_subgrid_max, livestock_subgrid_dim, limit_livestock_delta)
        call setup_buffer_zone("livestock", livestock_buffer_index, livestock_buffer_size, &
            livestock_subgrid_delta, livestock_subgrid_min, livestock_subgrid_max, livestock_subgrid_dim)
        
        ! Allocate arrays 
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

        ! Setup crossrefence subgrids
        call set_subgrid_xy("livestock", livestock_subgrid_dim, livestock_subgrid_min, livestock_subgrid_delta, x_livestock_subgrid, y_livestock_subgrid)
        call setup_crossref_grid("livestock", crossref_emission_to_livestock_subgrid, livestock_subgrid_delta, livestock_subgrid_min, livestock_subgrid_dim)
    end subroutine setup_livestock_arrays_and_variables

    subroutine read_livestock_data()
        call read_netcdf_data("livestock", pathname_livestock, filename_livestock, livestock_subgrid, &
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