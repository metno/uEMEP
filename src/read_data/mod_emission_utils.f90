module mod_emission_utils
    !! The module contains procedures for setting up and reading data for emission proxy subgrids

    use uemep_constants, only: dp
    use uEMEP_definitions, only: unit_logfile, x_dim_index, y_dim_index, x_dim_nc_index, y_dim_nc_index
    use uemep_configuration, only: projection_type, projection_attributes
    use mod_lambert_projection, only: proj2ll
    use netcdf

    implicit none
    private

    public :: read_netcdf_data, set_subgrid_dimensions, set_buffer_zone, set_crossref_grid, set_subgrid_xy

contains

    function open_netcdf_file(dirname, filename) result(ncid)
        !! Opens a netcdf file
        character(len=*), intent(in) :: dirname
        character(len=*), intent(in) :: filename
        integer :: ncid, ncstat

        character(len=512) :: filepath
        logical :: exist

        filepath = trim(dirname)//trim(filename)
        inquire(file=trim(filepath), exist=exist)
        if (.not. exist) then
            write(unit_logfile, "(2a)") "ERROR: NetCDF file not found: ", trim(filepath)
            stop 1
        end if

        write(unit_logfile, "(2a)") "Opening NetCDF file: ", trim(filepath)
        ncstat = nf90_open(filepath, nf90_nowrite, ncid)
        if (ncstat /= nf90_noerr) then
            write(unit_logfile, "(2a)") "ERROR: Could not open NetCDF file: ", trim(filepath)
            stop 1
        end if
    end function open_netcdf_file

    subroutine read_netcdf_data(data_name, dirname, filename, subgrid, subgrid_min, subgrid_max, subgrid_delta, subgrid_dim, &
            x_subgrid, y_subgrid, var_name)
        !! Reads in emission data from a netcdf file
        character(len=*), intent(in) :: data_name
        character(len=*), intent(in) :: dirname
        character(len=*), intent(in) :: filename
        real, intent(out) :: subgrid(:,:)
        real, intent(inout) :: subgrid_min(2), subgrid_max(2), subgrid_delta(2)
        integer, intent(in) :: subgrid_dim(2)
        real, intent(in) :: x_subgrid(:,:), y_subgrid(:,:)
        character(len=*), intent(in) :: var_name

        ! Local variables
        integer :: i, j
        integer :: ncstat, ncid, dimid, varid
        character(len=256) :: dummy
        character(len=3) :: dim_names(2) = ["lon", "lat"]
        integer :: dim_length(2), dim_start(2)
        real(dp), allocatable :: lonlat(:,:), ncdata(:,:)
        real :: delta(2), tmp_delta(2)
        real :: x, y, lon(3), lat(3)
        integer :: i_nearest, j_nearest

        write(unit_logfile, "(3a)") "Reading ", trim(data_name), " data"

        ! Open NetCDF file
        ncid = open_netcdf_file(dirname, filename)

        ! Get lon and lat dimension sizes for the entire domain
        do i = 1, 2
            ncstat = nf90_inq_dimid(ncid, dim_names(i), dimid)
            if (ncstat /= nf90_noerr) then
                write(unit_logfile, "(2a)") "ERROR: Dimension name not found: ", trim(dim_names(i))
                stop 1
            end if
            ncstat = nf90_inquire_dimension(ncid, dimid, dummy, dim_length(i))
            if (ncstat /= nf90_noerr) then
                write(unit_logfile, "(2a)") "ERROR: Could not read dimension for: ", trim(dim_names(i))
                stop 1
            endif
        end do
        write(unit_logfile, "(3a,2i8)") "Size of ", trim(data_name), " data (lon,lat): ", dim_length(1), dim_length(2)

        ! Reduce domain size
        call reduce_data_to_target_subgrid(data_name, ncid, 10.0, 10, dim_names, dim_start, dim_length, subgrid_min, subgrid_max, subgrid_delta)

        ! Allocate working arrays
        if (allocated(lonlat)) deallocate(lonlat)
        allocate(lonlat(max(dim_length(x_dim_nc_index), dim_length(y_dim_nc_index)), 2))
        if (allocated(ncdata)) deallocate(ncdata)
        allocate(ncdata(dim_length(x_dim_nc_index), dim_length(y_dim_nc_index)))

        ! Calculate delta
        do i = 1, 2
            ncstat = nf90_inq_varid(ncid, trim(dim_names(i)), varid)
            if (ncstat /= nf90_noerr) then
                write(unit_logfile, "(2a)") "ERROR: Variable for dimension name not found: ", trim(dim_names(i))
                stop 1
            end if
            ncstat = nf90_get_var(ncid, varid, lonlat(1:dim_length(i),i), start=[dim_start(i)], count=[dim_length(i)])
            if (ncstat /= nf90_noerr) then
                write(unit_logfile, "(2a)") "ERROR: Could not read variable for dimension name: ", trim(dim_names(i))
                stop 1
            end if
        end do
        delta = lonlat(2,:) - lonlat(1,:)

        ! Read data
        ncstat = nf90_inq_varid(ncid, trim(var_name), varid)
        if (ncstat /= nf90_noerr) then
            write(unit_logfile, "(2a)") "ERROR: No variable with name: ", trim(var_name)
            stop 1
        end if
        ncstat = nf90_get_var(ncid, varid, ncdata, start=[dim_start], count=[dim_length])
        if (ncstat /= nf90_noerr) then
            write(unit_logfile, "(2a)") "ERROR: Could not read variable with name: ", trim(var_name)
            stop 1
        end if

        ! Contrain data to be larger than zero if necessary
        if (minval(ncdata) < 0.0) then
            write(unit_logfile, "(3a)") "WARNING: Data for variable ", trim(var_name), " contains negative values"
            write(unit_logfile, "(a)") "WARNING: Setting negative values to zero"
        end if
        write(unit_logfile, "(3a,2f12.2)") "Data range for ", trim(var_name), " (min/max): ", &
            minval(ncdata), maxval(ncdata)

        subgrid(:,:) = 0.0
        do j = 1, subgrid_dim(y_dim_nc_index)
            do i = 1, subgrid_dim(x_dim_nc_index)
                ! Project the center position to lon/lat
                x = x_subgrid(i,j)
                y = y_subgrid(i,j)
                call proj2ll(x, y, lon(1), lat(1), projection_attributes, projection_type)

                ! Project both sides to get delta x
                x = x_subgrid(i,j) - 0.5*subgrid_delta(x_dim_index)
                y = y_subgrid(i,j)
                call proj2ll(x, y, lon(2), lat(2), projection_attributes, projection_type)
                x = x_subgrid(i,j) + 0.5*subgrid_delta(x_dim_index)
                y = y_subgrid(i,j)
                call proj2ll(x, y, lon(3), lat(3), projection_attributes, projection_type)
                tmp_delta(x_dim_index) = lon(3) - lon(2)

                ! Again for delta y
                x = x_subgrid(i,j)
                y = y_subgrid(i,j) - 0.5*subgrid_delta(y_dim_index)
                call proj2ll(x, y, lon(2), lat(2), projection_attributes, projection_type)
                x = x_subgrid(i,j)
                y = y_subgrid(i,j) + 0.5*subgrid_delta(y_dim_index)
                call proj2ll(x, y, lon(3), lat(3), projection_attributes, projection_type)
                tmp_delta(y_dim_index) = lat(3) - lat(2)

                ! Find nearest neighbour and insert value in subgrid
                i_nearest = 1 + floor((lon(1) - lonlat(1,x_dim_nc_index))/delta(1) + 0.5)
                j_nearest = 1 + floor((lat(1) - lonlat(1,y_dim_nc_index))/delta(2) + 0.5)

                subgrid(i,j) = ncdata(i_nearest,j_nearest)

                ! Constrain values
                if (subgrid(i,j) < 0.0) then
                    write(unit_logfile, "(a,2i4,2a)") "WARNING: Negative value at (i,j): ", i, j, " for variable: ", trim(var_name)
                    write(unit_logfile, "(a)") "WARNING: Setting negative value to zero"
                    subgrid(i,j) = 0.0
                end if

                ! Check for NaN
                if (isnan(subgrid(i,j))) then
                    write(unit_logfile, "(a,2i4,2a)") "ERROR: NaN value at (i,j): ", i, j, " for variable: ", trim(var_name)
                    stop 1
                end if
            end do
        end do
        write(unit_logfile, "(5a,2f12.2)") "Data range for variable: ", trim(var_name), " in ", trim(data_name), " subgrid (min/max): ", minval(subgrid), maxval(subgrid)

        if (allocated(lonlat)) deallocate(lonlat)
        if (allocated(ncdata)) deallocate(ncdata)
    end subroutine read_netcdf_data

    subroutine reduce_data_to_target_subgrid(data_name, ncid, buffer_delta, padding, dim_names, dim_start, dim_length, &
            subgrid_min, subgrid_max, subgrid_delta)
        !! Reduces the size of the domain to read so it corresponds to the target subgrid
        character(len=*), intent(in) :: data_name
        integer, intent(in) :: ncid
        real, intent(in) :: buffer_delta
        integer, intent(in) :: padding
        character(len=3), intent(in) :: dim_names(2)
        integer, intent(inout) :: dim_start(2), dim_length(2)
        real, intent(in) :: subgrid_min(2), subgrid_max(2), subgrid_delta(2)

        ! Local variables
        integer :: i
        integer :: ncstat, varid
        real :: x_min, x_max, y_min, y_max
        integer :: i_min, i_max, j_min, j_max
        real :: lon(4), lat(4)
        real(dp), allocatable :: lonlat(:,:)
        real :: delta(2)

        ! Sanity checks
        if (ncid <= 0) then
            write(unit_logfile, "(a)") "ERROR: NetCDF id has to be a positive number"
            stop 1
        end if
        if (buffer_delta < 0) then
            write(unit_logfile, "(a)") "ERROR: buffer delta cannot be a negative number"
            stop 1
        end if
        if (padding < 0) then
            write(unit_logfile, "(a)") "ERROR: padding cannot be a negative number"
            stop 1
        end if

        write(unit_logfile, "(3a)") "Reducing domain size for ", trim(data_name), " data before reading"

        ! Retrieve the four corners of the target grid in lon/lat
        x_min = subgrid_min(x_dim_index) - buffer_delta*subgrid_delta(x_dim_index)
        x_max = subgrid_max(x_dim_index) + buffer_delta*subgrid_delta(x_dim_index)
        y_min = subgrid_min(y_dim_index) - buffer_delta*subgrid_delta(y_dim_index)
        y_max = subgrid_max(y_dim_index) + buffer_delta*subgrid_delta(y_dim_index)
        call proj2ll(x_min, y_min, lon(1), lat(1), projection_attributes, projection_type)
        call proj2ll(x_max, y_max, lon(2), lat(2), projection_attributes, projection_type)
        call proj2ll(x_min, y_max, lon(3), lat(3), projection_attributes, projection_type)
        call proj2ll(x_max, y_min, lon(4), lat(4), projection_attributes, projection_type)
        x_min = minval(lon)
        x_max = maxval(lon)
        y_min = minval(lat)
        y_max = maxval(lat)
        write(unit_logfile, "(a,2f12.2)") "Longitude (min/max): ", x_min, x_max
        write(unit_logfile, "(a,2f12.2)") "Latitude (min/max): ", y_min, y_max

        ! Allocate working array
        if (allocated(lonlat)) deallocate(lonlat)
        allocate(lonlat(max(dim_length(x_dim_nc_index), dim_length(y_dim_nc_index)), 2))

        ! Read dimensions from NetCDF
        dim_start = 1
        do i = 1, 2
            ncstat = nf90_inq_varid(ncid, trim(dim_names(i)), varid)
            if (ncstat /= nf90_noerr) then
                write(unit_logfile, "(2a)") "ERROR: No dimension varible with name: ", trim(dim_names(i))
                stop 1
            end if
            ncstat = nf90_get_var(ncid, varid, lonlat(1:dim_length(i), i), start=[dim_start(i)], count=[dim_length(i)])
            if (ncstat /= nf90_noerr) then
                write(unit_logfile, "(2a)") "ERROR: Could not read dimension variable: ", trim(dim_names(i))
                stop 1
            end if
        end do
        delta = lonlat(2,:) - lonlat(1,:)
        write(unit_logfile, "(a,2f12.6)") "Delta degrees: ", delta

        ! Convert from lon/lat to grid indices
        i_min = 1 + floor((x_min - lonlat(1,1))/delta(1) + 0.5)
        i_max = 1 + floor((x_max - lonlat(1,1))/delta(1) + 0.5)
        j_min = 1 + floor((y_min - lonlat(1,2))/delta(2) + 0.5)
        j_max = 1 + floor((y_max - lonlat(1,2))/delta(2) + 0.5)

        ! Constrain new domain to be within the available data region
        i_min = max(1, i_min - padding)
        i_max = min(dim_length(x_dim_nc_index), i_max + padding)
        j_min = max(1, j_min - padding)
        j_max = min(dim_length(y_dim_nc_index), j_max + padding)

        ! Sanity checks
        if (i_min >= i_max) then
            write(unit_logfile, "(a)") "ERROR: Minimum exceeds maximum for index i"
            stop 1
        end if
        if (j_min >= j_max) then
            write(unit_logfile, "(a)") "ERROR: Mimimum exceeds maximum for index j"
            stop 1
        end if

        ! Set new reduced domain size
        dim_length(x_dim_nc_index) = i_max - i_min + 1
        dim_length(y_dim_nc_index) = j_max - j_min + 1
        dim_start(x_dim_nc_index) = i_min
        dim_start(y_dim_nc_index) = j_min

        write(unit_logfile, "(3a,3i)") "Reading ", trim(data_name), " data i grids (min/max/length): ", &
            i_min, i_max, dim_length(x_dim_nc_index)
        write(unit_logfile, "(3a,3i)") "Reading ", trim(data_name), " data j grids (min/max/length): ", &
            j_min, j_max, dim_length(y_dim_nc_index)
        write(unit_logfile, "(3a,2f12.2)") "Reading ", trim(data_name), " data longitudes (min/max): ", &
            lonlat(i_min,x_dim_nc_index), lonlat(i_max,x_dim_nc_index)
        write(unit_logfile, "(3a,2f12.2)") "Reading ", trim(data_name), " data latitudes (min/max): ", &
            lonlat(j_min,y_dim_nc_index), lonlat(j_max,y_dim_nc_index)

        if (allocated(lonlat)) deallocate(lonlat)
    end subroutine reduce_data_to_target_subgrid

    subroutine set_subgrid_dimensions(name, sg_delta, sg_min, sg_max, sg_dim, delta_limit)
        !! Sets up the dimensions of the proxy subgrid
        use uemep_configuration, only: subgrid_delta, subgrid_min, subgrid_max
        use uemep_definitions, only: subgrid_dim
        character(len=*), intent(in) :: name
        real, intent(out) :: sg_delta(:)
        real, intent(out) :: sg_min(:)
        real, intent(out) :: sg_max(:)
        integer, intent(out) :: sg_dim(:)
        real, intent(in) :: delta_limit

        write(unit_logfile,"(3a)") "Setting up dimensions for the ", trim(name), " subgrid"
        sg_delta(x_dim_index) = max(subgrid_delta(x_dim_index), delta_limit)
        sg_delta(y_dim_index) = max(subgrid_delta(y_dim_index), delta_limit)
        sg_min = subgrid_min
        sg_max = subgrid_max
        sg_dim(x_dim_index) = floor((sg_max(x_dim_index) - sg_min(x_dim_index))/sg_delta(x_dim_index))
        sg_dim(y_dim_index) = floor((sg_max(y_dim_index) - sg_min(y_dim_index))/sg_delta(y_dim_index))
        sg_dim(x_dim_index) = max(min(sg_dim(x_dim_index), subgrid_dim(x_dim_index)), 1)
        sg_dim(y_dim_index) = max(min(sg_dim(y_dim_index), subgrid_dim(y_dim_index)), 1)
    end subroutine set_subgrid_dimensions

    subroutine set_buffer_zone(name, b_index, b_size, sg_delta, sg_min, sg_max, sg_dim)
        !! Sets up a buffer zone around the proxy subgrid
        use uemep_configuration, only: local_subgrid_method_flag
        use uemep_definitions, only: use_buffer_zone, buffer_index_scale
        use define_subgrid, only: dx => dx_temp, dy => dy_temp
        character(len=*), intent(in) :: name
        integer, intent(out) :: b_index(:)
        real, intent(out) :: b_size(:)
        real, intent(in) :: sg_delta(:)
        real, intent(inout) :: sg_min(:)
        real, intent(inout) :: sg_max(:)
        integer, intent(inout) :: sg_dim(:)

        ! Local variables
        real :: offset
        
        if (use_buffer_zone) then
            write(unit_logfile,"(3a)") "Setting up buffer zone for the ", trim(name), " subgrid"
            if (local_subgrid_method_flag == 3) then
                offset = 1.0
            else
                offset = 0.5
            end if
            b_index(x_dim_index) = floor(dx/sg_delta(x_dim_index)*(buffer_index_scale + offset))
            b_index(y_dim_index) = floor(dy/sg_delta(y_dim_index)*(buffer_index_scale + offset))
        else
            b_index = 0.0
        end if

        b_size = b_index * sg_delta
        sg_dim(1:2) = sg_dim(1:2) + b_index(1:2)*2
        sg_min(1:2) = sg_min(1:2) - b_size(1:2)
        sg_max(1:2) = sg_max(1:2) + b_size(1:2)
    end subroutine set_buffer_zone

    subroutine set_crossref_grid(name, crossref_sg, sg_delta, sg_min, sg_dim)
        !! Sets up a crossreference grid to map between the emissions and proxy subgrids
        use uemep_definitions, only: livestock_index, emission_subgrid_dim, x_emission_subgrid, y_emission_subgrid
        character(len=*), intent(in) :: name
        integer, allocatable, intent(inout) :: crossref_sg(:,:,:)
        real, intent(in) :: sg_delta(:)
        real, intent(in) :: sg_min(:)
        integer, intent(in) :: sg_dim(:)

        ! Local variables
        integer :: i, j

        if (.not. allocated(crossref_sg)) then
            write(unit_logfile,"(3a)") "ERROR: Crossreference grid for ", trim(name), " is not allocated"
            stop 1
        end if

        write(unit_logfile, "(3a)") "Crossreferencing emissions to ", trim(name), " subgrid"

        do j = 1, emission_subgrid_dim(y_dim_index,livestock_index)
            do i = 1, emission_subgrid_dim(x_dim_index,livestock_index)
                crossref_sg(i,j,x_dim_index) = &
                    1 + floor((x_emission_subgrid(i,j,livestock_index) &
                    - sg_min(x_dim_index))/sg_delta(x_dim_index))
                crossref_sg(i,j,y_dim_index) = &
                    1 + floor((y_emission_subgrid(i,j,livestock_index) &
                    - sg_min(y_dim_index))/sg_delta(y_dim_index))

                ! Avoid invalid values at the grid edge
                crossref_sg(i,j,x_dim_index) = &
                    max(min(crossref_sg(i,j,x_dim_index), &
                    sg_dim(x_dim_index)), 1)
                crossref_sg(i,j,y_dim_index) = &
                    max(min(crossref_sg(i,j,y_dim_index), &
                    sg_dim(y_dim_index)), 1)
            end do
        end do
    end subroutine set_crossref_grid

    subroutine set_subgrid_xy(name, sg_dim, sg_min, sg_delta, x_sg, y_sg)
        !! Sets up the proxy subgrid x and y values
        character(len=*), intent(in) :: name
        integer, intent(in) :: sg_dim(:)
        real, intent(in) :: sg_min(:)
        real, intent(in) :: sg_delta(:)
        real, allocatable, intent(inout) :: x_sg(:,:)
        real, allocatable, intent(inout) :: y_sg(:,:)
        
        ! Local variables
        integer :: i, j

        if (.not. allocated(x_sg)) then
            write(unit_logfile, "(3a)") "ERROR: x_", trim(name), "_subgrid is not allocated"
            stop 1
        end if
        if (.not. allocated(y_sg)) then
            write(unit_logfile, "(3a)") "ERROR: y_", trim(name), "_subgrid is not allocated"
            stop 1
        end if
        
        do j = 1, sg_dim(y_dim_index)
            do i = 1, sg_dim(x_dim_index)
                x_sg(i,j) = sg_min(x_dim_index) + sg_delta(x_dim_index)*(i - 0.5)
                y_sg(i,j) = sg_min(y_dim_index) + sg_delta(y_dim_index)*(j - 0.5)
            end do
        end do
    end subroutine set_subgrid_xy

end module mod_emission_utils