module pollen_proxy_data

    !! This module provides procedures for reading pollen proxy data
    !!
    !! Copyright (C) 2007 Free Software Foundation.
    !! License GNU LGPL-3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>.
    !! This is free software: you are free to change and redistribute it.
    !!
    !! Developed and maintained at the Norwegian Meteorological Institute.
    !! Contribute at: <https://github.com/metno/uEMEP>

    use error_handling
    use uemep_constants, only: dp, pi
    use uemep_configuration, only: pathfilename_pollen, pathname_pollen, filename_pollen, var_name_pollen_nc, &
        projection_attributes, projection_type, calculate_source
    use uEMEP_definitions, only: unit_logfile, x_dim_nc_index, y_dim_nc_index, x_dim_index, y_dim_index, &
        num_dims_pollen_nc, dim_name_pollen_nc, dim_length_pollen_nc, dim_start_pollen_nc, &
        pollen_subgrid, pollen_subgrid_dim, x_pollen_subgrid, y_pollen_subgrid, &
        pollen_subgrid_delta, pollen_subgrid_min, pollen_subgrid_max, birch_proxy_index, &
        n_source_index, emission_subgrid_dim, crossreference_emission_to_pollen_subgrid, proxy_emission_subgrid
    use netcdf
    use mod_lambert_projection, only: proj2ll

    implicit none

    private

    public :: read_pollen_proxy, redistribute_pollen_emissions

    interface read_pollen_proxy
        !! Read pollen proxy data from netcdf
        module procedure read_pollen_proxy_data_latlon
    end interface read_pollen_proxy

contains

    subroutine read_pollen_proxy_data_latlon(i_pollen)
        !! Read pollen proxy data from netcdf and puts it in the pollen proxy subgrid
        integer, intent(in) :: i_pollen !! Pollen proxy index

        integer :: ncstat, ncid, dimid, varid
        integer :: i_dim, i, j, i_nearest, j_nearest
        logical :: file_exists, cond
        character(len=256) :: tmp_string, varname
        real :: tmp_lon(3), tmp_lat(3), x, y
        real :: tmp_delta(num_dims_pollen_nc), delta_nc(num_dims_pollen_nc)
        real(dp), allocatable :: lonlat_nc(:,:)
        real(dp), allocatable :: pollen_nc(:,:)

        write(unit_logfile, "(a)") ""
        write(unit_logfile, "(a)") "================================================================"
        write(unit_logfile, "(a,i0)") "Reading pollen proxy data for source number: ", i_pollen
        write(unit_logfile, "(a)") "================================================================"

        pathfilename_pollen(i_pollen) = trim(pathname_pollen(i_pollen))//trim(filename_pollen(i_pollen))
        inquire(file=trim(pathfilename_pollen(i_pollen)), exist=file_exists)
        call assert(file_exists, " NetCDF file does not exist"//trim(pathfilename_pollen(i_pollen)), code=file_not_found)

        write(unit_logfile, "(2a)") " Opening NetCDF file: ", trim(pathfilename_pollen(i_pollen))
        ncstat = nf90_open(pathfilename_pollen(i_pollen), nf90_nowrite, ncid)
        call assert((ncstat == nf90_noerr), " Could not open NetCDF file: "//trim(pathfilename_pollen(i_pollen)), code=read_error)

        ! Get dimensions for the entire dataset from the NetCDF
        do i_dim = 1, num_dims_pollen_nc
            ncstat = nf90_inq_dimid(ncid, dim_name_pollen_nc(i_dim), dimid)
            call assert((ncstat == nf90_noerr), " No dimension with name: " &
                // trim(dim_name_pollen_nc(i_dim)), code=read_error)
            ncstat = nf90_inquire_dimension(ncid, dimid, tmp_string, dim_length_pollen_nc(i_dim))
            call assert((ncstat == nf90_noerr), " No dimension information available for: " &
                // trim(dim_name_pollen_nc(i_dim)), code=read_error)
        end do

        write(unit_logfile,"(a,2i8)") " Size of pollen proxy dimensions (lon,lat): ", dim_length_pollen_nc

        ! We reduce the domain to read according to the target subgrid to reduce read time and memory use
        call reduce_pollen_proxy_region(ncid, i_pollen, buffer_delta=10.0, padding=10)

        if (.not. allocated(lonlat_nc)) then
            allocate(lonlat_nc(max(dim_length_pollen_nc(x_dim_nc_index), dim_length_pollen_nc(y_dim_nc_index)), num_dims_pollen_nc))
        end if
        if (.not. allocated(pollen_nc)) then
            allocate(pollen_nc(dim_length_pollen_nc(x_dim_nc_index), dim_length_pollen_nc(y_dim_nc_index)))
        end if

        ! Read from NetCDF

        ! First we read lon/lat (again) to calculate the delta
        do i = 1, num_dims_pollen_nc
            varname = dim_name_pollen_nc(i)
            ncstat = nf90_inq_varid(ncid, trim(varname), varid)
            call assert((ncstat == nf90_noerr), " No dimension with name: "//trim(varname), code=read_error)
            ncstat = nf90_get_var(ncid, varid, lonlat_nc(1:dim_length_pollen_nc(i),i), &
                start=[dim_start_pollen_nc(i)], count=[dim_length_pollen_nc(i)])
            call assert((ncstat == nf90_noerr), " Could not read dimension with name: "//trim(varname), code=read_error)
        end do
        delta_nc = lonlat_nc(2,:) - lonlat_nc(1,:)

        ! Then we read the pollen proxy data
        varname = var_name_pollen_nc(i_pollen)
        ncstat = nf90_inq_varid(ncid, trim(varname), varid)
        call assert((ncstat == nf90_noerr), " No variable with name: "//trim(varname), code=read_error)
        ncstat = nf90_get_var(ncid, varid, pollen_nc, start=[dim_start_pollen_nc], count=[dim_length_pollen_nc])
        call assert((ncstat == nf90_noerr), " Cannot read variable: "//trim(varname), code=read_error)
        call assert((minval(pollen_nc) < 0.0), " Pollen proxy cannot contain negative values: "//trim(varname), code=invalid_value)

        write(unit_logfile,"(3a,2f12.2)") " Pollen proxy min and max: ", trim(varname), " ", minval(pollen_nc), maxval(pollen_nc)

        ! Loop through the pollen proxy grid and put them in the pollen proxy subgrid grid
        ! Lat/lon is converted to subgrid coordinates and the inserted value is found using nearest neighbour
        pollen_subgrid(:,:,i_pollen) = 0.0
        do j = 1, pollen_subgrid_dim(y_dim_nc_index)
            do i = 1, pollen_subgrid_dim(x_dim_nc_index)
                ! Project the center position to lon/lat
                x = x_pollen_subgrid(i,j)
                y = y_pollen_subgrid(i,j)
                call proj2ll(x, y, tmp_lon(1), tmp_lat(1), projection_attributes, projection_type)
                
                ! Project both sides to get delta x
                x = x_pollen_subgrid(i,j) - 0.5*pollen_subgrid_delta(x_dim_index)
                y = y_pollen_subgrid(i,j)
                call proj2ll(x, y, tmp_lon(2), tmp_lat(2), projection_attributes, projection_type)
                x = x_pollen_subgrid(i,j) + 0.5*pollen_subgrid_delta(x_dim_index)
                y = y_pollen_subgrid(i,j)
                call PROJ2LL(x, y, tmp_lon(3), tmp_lat(3), projection_attributes, projection_type)
                tmp_delta(x_dim_index) = tmp_lon(3) - tmp_lon(2)

                ! Project both sides to get delta y
                x = x_pollen_subgrid(i,j)
                y = y_pollen_subgrid(i,j) - 0.5*pollen_subgrid_delta(y_dim_index)
                call proj2ll(x, y, tmp_lon(2), tmp_lat(2), projection_attributes, projection_type)
                x = x_pollen_subgrid(i,j)
                y = y_pollen_subgrid(i,j) + 0.5*pollen_subgrid_delta(y_dim_index)
                call proj2ll(x, y, tmp_lon(3), tmp_lat(3), projection_attributes, projection_type)
                tmp_delta(y_dim_index) = tmp_lat(3) - tmp_lat(2)

                ! Find nearest neighbour and insert value in subgrid
                i_nearest = 1 + floor((tmp_lon(1) - lonlat_nc(1,x_dim_nc_index))/delta_nc(1) + 0.5)
                j_nearest = 1 + floor((tmp_lat(1) - lonlat_nc(1,y_dim_nc_index))/delta_nc(2) + 0.5)
                pollen_subgrid(i,j,i_pollen) = pollen_nc(i_nearest,j_nearest)!/100.0

                if (pollen_subgrid(i,j,i_pollen) > 0.5) print *, pollen_subgrid(i,j,i_pollen)

                if (pollen_subgrid(i,j,i_pollen) < 0.0) pollen_subgrid(i,j,i_pollen) = 0.0

                cond = .not. isnan(pollen_subgrid(i,j,i_pollen))
                call assert(cond, " NaN in pollen proxy subgrid", code=invalid_value)

                ! cond = (pollen_subgrid(i,j,i_pollen) > 0.0)
                ! print *, cond
                ! call assert(cond, " Negative number of pollen proxy subgrid", code=invalid_value)
            end do
        end do

        if (allocated(lonlat_nc)) deallocate(lonlat_nc)
        if (allocated(pollen_nc)) deallocate(pollen_nc)
    end subroutine read_pollen_proxy_data_latlon


    subroutine reduce_pollen_proxy_region(ncid, i_pollen, buffer_delta, padding)
        !! Reduce the area that should be read to the size of the target subgrid
        integer, intent(in) :: ncid !! ID of open NetCDF file
        integer, intent(in) :: i_pollen !! Index of the pollen species
        real, intent(in) :: buffer_delta !! Relative buffer delta to add when reading the domain corner coordinates
        integer, intent(in) :: padding !! Number of grid cells to add as padding around the read area

        integer :: i_dim, ncstat, varid
        integer :: i_min, i_max, j_min, j_max
        real :: x_min, x_max, y_min, y_max
        real :: tmp_lon(4), tmp_lat(4), tmp_x(4), tmp_y(4)
        real :: delta(num_dims_pollen_nc)
        character(len=256) :: varname

        real(dp), allocatable :: lonlat(:,:)

        call assert((ncid > 0), " NetCDF ID has to be positive", code=read_error)
        call assert((i_pollen > 0), " Invalid pollen index", code=index_error)
        call assert((buffer_delta >= 0.0), " Buffer delta cannot be a negative number", code=invalid_value)
        call assert(.not. allocated(lonlat), " Temporary reading array is already allocated", code=allocation_error)

        write(unit_logfile, "(a)") " Reducing pollen domain for reading"

        ! Retrieve the four corners of the target grid in lon/lat
        x_min = pollen_subgrid_min(x_dim_index) - buffer_delta*pollen_subgrid_delta(x_dim_index)
        x_max = pollen_subgrid_max(x_dim_index) + buffer_delta*pollen_subgrid_delta(x_dim_index)
        y_min = pollen_subgrid_min(y_dim_index) - buffer_delta*pollen_subgrid_delta(y_dim_index)
        y_max = pollen_subgrid_max(y_dim_index) + buffer_delta*pollen_subgrid_delta(y_dim_index)
        
        call proj2ll(x_min, y_min, tmp_lon(1), tmp_lat(1), projection_attributes, projection_type)
        call proj2ll(x_max, y_max, tmp_lon(2), tmp_lat(2), projection_attributes, projection_type)
        call proj2ll(x_min, y_max, tmp_lon(3), tmp_lat(3), projection_attributes, projection_type)
        call proj2ll(x_max, y_min, tmp_lon(4), tmp_lat(4), projection_attributes, projection_type)

        tmp_x = tmp_lon
        tmp_y = tmp_lat
        x_min = minval(tmp_x)
        x_max = maxval(tmp_y)
        y_min = minval(tmp_y)
        y_max = maxval(tmp_y)

        write(unit_logfile, "(a,2f12.2)") " Min (lon,lat): ", x_min, y_min
        write(unit_logfile, "(a,2f12.2)") " Max (lon,lat): ", x_max, y_max

        ! Allocate field for storing lon/lat
        if (.not. allocated(lonlat)) then
            allocate(lonlat(max(dim_length_pollen_nc(x_dim_index), dim_length_pollen_nc(y_dim_index)), num_dims_pollen_nc))
        end if

        ! Read dimensions from NetCDF
        dim_start_pollen_nc = 1
        do i_dim = 1, num_dims_pollen_nc
            varname = dim_name_pollen_nc(i_dim)
            ncstat = nf90_inq_varid(ncid, trim(varname), varid)
            call assert((ncstat == nf90_noerr), " No dimension with name: "//trim(varname), code=read_error)
            ncstat = nf90_get_var(ncid, varid, lonlat(1:dim_length_pollen_nc(i_dim), i_dim), &
                start=[dim_start_pollen_nc(i_dim)], count=[dim_length_pollen_nc(i_dim)])
            call assert((ncstat == nf90_noerr), " Could not read dimension: "//trim(varname), code=read_error)
        end do

        delta = lonlat(2,:) - lonlat(1,:)
        write(unit_logfile,"(a,2f12.6)") " Pollen proxy delta (degrees): ", delta

        ! We re-use the variables, and convert lon/lat into grid indices
        i_min = 1 + floor((x_min - lonlat(1,1))/delta(1) + 0.5)
        i_max = 1 + floor((x_max - lonlat(1,1))/delta(1) + 0.5)
        j_min = 1 + floor((y_min - lonlat(1,2))/delta(2) + 0.5)
        j_max = 1 + floor((y_max - lonlat(1,2))/delta(2) + 0.5)

        ! Constrain the new domain to be within the available data region
        i_min = max(1, i_min - padding)
        i_max = min(dim_length_pollen_nc(x_dim_nc_index), i_max + padding)
        j_min = max(1, j_min - padding)
        j_max = min(dim_length_pollen_nc(y_dim_nc_index), j_max + padding)

        call assert((i_max >= i_min), " Minimum domain extend for i cannot exceed the maximum", code=index_error)
        call assert((j_max >= j_min), " Minimum domain extend for j cannot exceed the maximum", code=index_error)

        ! Set the new, reduced domain size to be read
        dim_length_pollen_nc(x_dim_nc_index) = i_max - i_min + 1
        dim_length_pollen_nc(y_dim_nc_index) = j_max - j_min + 1
        dim_start_pollen_nc(x_dim_nc_index) = i_min
        dim_start_pollen_nc(y_dim_nc_index) = j_min

        write(unit_logfile, "(a,3i)") " Reading pollen proxy i grids: ", i_min, i_max, dim_length_pollen_nc(x_dim_nc_index)
        write(unit_logfile, "(a,3i)") " Reading pollen proxy j grids: ", j_min, j_max, dim_length_pollen_nc(y_dim_nc_index)
        write(unit_logfile, "(a,2f12.2)") " Reading pollen proxy lon grids: ", &
            lonlat(i_min,x_dim_nc_index), lonlat(i_max,x_dim_nc_index)
        write(unit_logfile, "(a,2f12.2)") " Reading pollen proxy lat grids: ", &
            lonlat(j_min,y_dim_nc_index), lonlat(j_max,y_dim_nc_index)

        if (allocated(lonlat)) deallocate(lonlat)
    end subroutine reduce_pollen_proxy_region

    subroutine redistribute_pollen_emissions()
        integer :: i_source, i, j, i_pollen_index, j_pollen_index

        do i_source = 1, n_source_index
            if (calculate_source(i_source)) then
                proxy_emission_subgrid(:,:,i_source,:) = 0.0
                do j = 1, emission_subgrid_dim(y_dim_nc_index,i_source)
                    do i = 1, emission_subgrid_dim(x_dim_nc_index,i_source)
                        
                        ! Note that this only works with birch as a source for now - extending later
                        i_pollen_index = crossreference_emission_to_pollen_subgrid(i,j,x_dim_index)
                        j_pollen_index = crossreference_emission_to_pollen_subgrid(i,j,y_dim_index)

                        proxy_emission_subgrid(i,j,i_source,:) = pollen_subgrid(i_pollen_index,j_pollen_index,birch_proxy_index)
                    end do
                end do
            end if
        end do

        write(unit_logfile,"(a,2f12.2)") " Pollen emission proxy min and max: ", &
            minval(proxy_emission_subgrid(:,:,17,:)), maxval(proxy_emission_subgrid(:,:,17,:))
    end subroutine redistribute_pollen_emissions

end module pollen_proxy_data