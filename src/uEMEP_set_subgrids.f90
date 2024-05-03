module set_subgrids

    use uemep_configuration
    use uEMEP_definitions
    use utility_functions, only: ll2utm, ll2ltm
    use mod_lambert_projection, only: LL2LAEA

    implicit none
    private

    public :: uEMEP_set_subgrids, uEMEP_set_subgrid_select_latlon_centre

contains

    subroutine uEMEP_set_subgrids()
        
        ! Local variables
        integer :: i
        real :: dim_check(2)

        ! In the case of interpolation and auto subgridding we need to extend the domain by dx and dy to get the different grids to fit
        ! Assume the maximum is 1 km
        if (use_emission_positions_for_auto_subgrid_flag(allsource_index)) then
            ! Check that it is divisable by the largest grid size
            dim_check(x_dim_index) = mod((subgrid_max(x_dim_index) - subgrid_min(x_dim_index)),max_interpolation_subgrid_size)
            dim_check(y_dim_index) = mod((subgrid_max(y_dim_index) - subgrid_min(y_dim_index)),max_interpolation_subgrid_size)
            if (dim_check(x_dim_index) .eq. 0) then
                subgrid_max(x_dim_index) = subgrid_max(x_dim_index)+subgrid_delta(x_dim_index)
                write(unit_logfile,'(A,f12.2)') 'Setting subgrids for auto emission gridding. Adding to x grid max: ', subgrid_delta(x_dim_index)
            else
                subgrid_max(x_dim_index) = subgrid_max(x_dim_index) + (max_interpolation_subgrid_size - dim_check(x_dim_index) + subgrid_delta(x_dim_index))
                write(unit_logfile,'(A,f12.2)') 'Setting subgrids for auto emission gridding. Adding to x grid max: ', (max_interpolation_subgrid_size - dim_check(x_dim_index) + subgrid_delta(x_dim_index))
            end if
            if (dim_check(x_dim_index) .eq. 0) then
                subgrid_max(y_dim_index) = subgrid_max(y_dim_index) + subgrid_delta(y_dim_index)
                write(unit_logfile,'(A,f12.2)') 'Setting subgrids for auto emission gridding. Adding to y grid max: ', subgrid_delta(y_dim_index)
            else
                subgrid_max(y_dim_index) = subgrid_max(y_dim_index) + (max_interpolation_subgrid_size - dim_check(y_dim_index) + subgrid_delta(y_dim_index))
                write(unit_logfile,'(A,f12.2)') 'Setting subgrids for auto emission gridding. Adding to y grid max: ', (max_interpolation_subgrid_size - dim_check(y_dim_index) + subgrid_delta(y_dim_index))
            end if
        end if

        ! Reset min and max with the buffer and calculate dimensions
        subgrid_dim(x_dim_index) = floor((subgrid_max(x_dim_index) - subgrid_min(x_dim_index))/subgrid_delta(x_dim_index))
        subgrid_dim(y_dim_index) = floor((subgrid_max(y_dim_index) - subgrid_min(y_dim_index))/subgrid_delta(y_dim_index))

        ! Set all integral subgrids relative to the target subgrid
        if (integral_subgrid_delta_ref .eq. 0.) then
            integral_subgrid_delta = subgrid_delta*integral_subgrid_step
        else
            integral_subgrid_delta(x_dim_index) = max(integral_subgrid_delta_ref, subgrid_delta(x_dim_index))
            integral_subgrid_delta(y_dim_index) = max(integral_subgrid_delta_ref, subgrid_delta(y_dim_index))
            integral_subgrid_step = floor(integral_subgrid_delta_ref/subgrid_delta(x_dim_index) + 0.5)
        end if

        integral_subgrid_min = subgrid_min
        integral_subgrid_max = subgrid_max
        integral_subgrid_dim(x_dim_index) = floor((subgrid_max(x_dim_index) - subgrid_min(x_dim_index))/integral_subgrid_delta(x_dim_index))
        integral_subgrid_dim(y_dim_index) = floor((subgrid_max(y_dim_index) - subgrid_min(y_dim_index))/integral_subgrid_delta(y_dim_index))
        integral_subgrid_dim(t_dim_index) = subgrid_dim(t_dim_index)
        
        ! Set the integral subgrid dimensions so they cannot be larger than the target subgrid
        integral_subgrid_dim(x_dim_index) = min(integral_subgrid_dim(x_dim_index), subgrid_dim(x_dim_index))
        integral_subgrid_dim(y_dim_index) = min(integral_subgrid_dim(y_dim_index), subgrid_dim(y_dim_index))

        ! Set all population subgrids relative to the target subgrid
        if (population_data_type .eq. population_index .and. .not. use_region_select_and_mask_flag) then
            ! When not using regional mask then 250 m population data is used then set this as a limit
            population_subgrid_delta(x_dim_index) = max(subgrid_delta(x_dim_index), limit_population_delta)
            population_subgrid_delta(y_dim_index) = max(subgrid_delta(y_dim_index), limit_population_delta)
        else
            ! Allow the population data to have the same grid as the target grid
            population_subgrid_delta(x_dim_index) = subgrid_delta(x_dim_index)
            population_subgrid_delta(y_dim_index) = subgrid_delta(y_dim_index)
        end if

        population_subgrid_min = subgrid_min
        population_subgrid_max = subgrid_max
        population_subgrid_dim(x_dim_index) = floor((population_subgrid_max(x_dim_index) - population_subgrid_min(x_dim_index))/population_subgrid_delta(x_dim_index))
        population_subgrid_dim(y_dim_index) = floor((population_subgrid_max(y_dim_index) - population_subgrid_min(y_dim_index))/population_subgrid_delta(y_dim_index))
        
        ! Set the population subgrid dimensions so they cannot be larger than the target subgrid. Not certain why I do this.
        population_subgrid_dim(x_dim_index) = min(population_subgrid_dim(x_dim_index), subgrid_dim(x_dim_index))
        population_subgrid_dim(y_dim_index) = min(population_subgrid_dim(y_dim_index), subgrid_dim(y_dim_index))
        
        ! Set population subgrid so it has a minimum of 1 dimensions, to avoid problems when running receptor calculations
        population_subgrid_dim(x_dim_index) = max(population_subgrid_dim(x_dim_index), 1)
        population_subgrid_dim(y_dim_index) = max(population_subgrid_dim(y_dim_index), 1)

        ! Set all emission subgrids to be the same as the target subgrid
        emission_max_subgrid_dim = subgrid_dim
        do i = 1, n_source_index
            emission_subgrid_delta(:,i) = subgrid_delta
            emission_subgrid_min(:,i) = subgrid_min
            emission_subgrid_max(:,i) = subgrid_max
            emission_subgrid_dim(:,i) = subgrid_dim
        end do

        ! Set shipping data to a minimum value for all sources (Cannot be smaller than the target subgrid)
        emission_subgrid_delta(x_dim_index,shipping_index) = max(subgrid_delta(x_dim_index), limit_shipping_delta)
        emission_subgrid_delta(y_dim_index,shipping_index) = max(subgrid_delta(y_dim_index), limit_shipping_delta)
        emission_subgrid_delta(x_dim_index,heating_index) = max(subgrid_delta(x_dim_index), limit_heating_delta)
        emission_subgrid_delta(y_dim_index,heating_index) = max(subgrid_delta(y_dim_index), limit_heating_delta)
        emission_subgrid_delta(x_dim_index,industry_index) = max(subgrid_delta(x_dim_index), limit_industry_delta)
        emission_subgrid_delta(y_dim_index,industry_index) = max(subgrid_delta(y_dim_index), limit_industry_delta)

        ! Set all the emission subgrid dimensions after changes
        do i = 1, n_source_index
            emission_subgrid_dim(x_dim_index,i) = floor((emission_subgrid_max(x_dim_index,i) - emission_subgrid_min(x_dim_index,i))/emission_subgrid_delta(x_dim_index,i))
            emission_subgrid_dim(y_dim_index,i) = floor((emission_subgrid_max(y_dim_index,i) - emission_subgrid_min(y_dim_index,i))/emission_subgrid_delta(y_dim_index,i))
            write(unit_logfile,'(A,I6,A5,2I6)') 'Emission grid dimensions for source ',i,': ', emission_subgrid_dim(1:2,i)
        end do

        ! Set the landuse and deposition grids to have the same extent as the target grid
        if (calculate_deposition_flag) then
            if (deposition_subgrid_delta(x_dim_index) .eq. 0) deposition_subgrid_delta(x_dim_index) = subgrid_delta(x_dim_index)
            if (deposition_subgrid_delta(y_dim_index) .eq. 0) deposition_subgrid_delta(y_dim_index) = subgrid_delta(y_dim_index)
            deposition_subgrid_min(:) = subgrid_min
            deposition_subgrid_max(:) = subgrid_max
            deposition_subgrid_dim(x_dim_index) = floor((subgrid_max(x_dim_index) - subgrid_min(x_dim_index))/deposition_subgrid_delta(x_dim_index))
            deposition_subgrid_dim(y_dim_index) = floor((subgrid_max(y_dim_index) - subgrid_min(y_dim_index))/deposition_subgrid_delta(y_dim_index))
            deposition_subgrid_dim(t_dim_index) = subgrid_dim(t_dim_index)
        end if

        if (read_landuse_flag) then
            if (landuse_subgrid_delta(x_dim_index) .eq. 0) landuse_subgrid_delta(x_dim_index) = subgrid_delta(x_dim_index)
            if (landuse_subgrid_delta(y_dim_index) .eq. 0) landuse_subgrid_delta(y_dim_index) = subgrid_delta(y_dim_index)
            landuse_subgrid_min(:) = subgrid_min
            landuse_subgrid_max(:) = subgrid_max
            landuse_subgrid_dim(x_dim_index) = floor((subgrid_max(x_dim_index) - subgrid_min(x_dim_index))/landuse_subgrid_delta(x_dim_index))
            landuse_subgrid_dim(y_dim_index) = floor((subgrid_max(y_dim_index) - subgrid_min(y_dim_index))/landuse_subgrid_delta(y_dim_index))
            landuse_subgrid_dim(t_dim_index) = subgrid_dim(t_dim_index)
        end if

        write(unit_logfile,'(A,2I6)') 'Concentration grid dimensions: ', subgrid_dim(1:2)
        write(unit_logfile,'(A,2I6)') 'Integral grid dimensions: ', integral_subgrid_dim(1:2)
        write(unit_logfile,'(A,2f10.1)') 'Concentration subgrid grid sizes: ', subgrid_delta
        write(unit_logfile,'(A,I6,2f10.1)') 'Integral subgrid step and grid sizes: ', integral_subgrid_step, integral_subgrid_delta

        if (calculate_deposition_flag) then
            write(unit_logfile,'(A,2I6)') 'Deposition grid dimensions: ', deposition_subgrid_dim(1:2)
            write(unit_logfile,'(A,2f10.1)') 'Deposition subgrid grid sizes: ', deposition_subgrid_delta
        end if

        if (read_landuse_flag) then
            write(unit_logfile,'(A,2I6)') 'Landuse grid dimensions: ', landuse_subgrid_dim(1:2)
            write(unit_logfile,'(A,2f10.1)') 'Landuse subgrid grid sizes: ', landuse_subgrid_delta
        end if

    end subroutine uEMEP_set_subgrids

    subroutine uEMEP_set_subgrid_select_latlon_centre()
        ! If specified using select_latlon_centre_domain_position_flag then this routines specifies
        ! the grid according to the lat lon position and the width and height
        ! This is intended to make life easier for users and to implement uEMEP in a more global context

        ! Local variables
        real :: x_out, y_out

        ! Find centre position in specified coordinates
        if (projection_type .eq. UTM_projection_index) then
            call ll2utm(1, utm_zone, select_lat_centre_position, select_lon_centre_position, y_out, x_out)
        else if (projection_type .eq. LTM_projection_index) then
            call ll2ltm(1, ltm_lon0, select_lat_centre_position, select_lon_centre_position, y_out, x_out)
        else if (projection_type .eq. LAEA_projection_index) then
            call LL2LAEA(x_out, y_out, select_lon_centre_position, select_lat_centre_position, projection_attributes)
        end if

        ! Snap to nearest 1 km
        x_out = floor(x_out/1000.0 + 0.5)*1000.0
        y_out = floor(y_out/1000.0 + 0.5)*1000.0

        ! Set max and min values
        subgrid_min(x_dim_index) = x_out - select_domain_width_EW_km*1000.0/2.0
        subgrid_min(y_dim_index) = y_out - select_domain_height_NS_km*1000.0/2.0
        subgrid_max(x_dim_index) = x_out + select_domain_width_EW_km*1000.0/2.0
        subgrid_max(y_dim_index) = y_out + select_domain_height_NS_km*1000.0/2.0

        write(unit_logfile,'(A,2f12.1)') 'Setting domain centre (lon,lat) to: ', select_lon_centre_position, select_lat_centre_position
        write(unit_logfile,'(A,2f12.1)') 'Setting min domain (x,y) to: ', subgrid_min(1:2)
        write(unit_logfile,'(A,2f12.1)') 'Setting max domain (x,y) to: ', subgrid_max(1:2)

        ! Reset the initial subgrid as well, needed for EMEP and receptor selection
        init_subgrid_min(x_dim_index) = subgrid_min(x_dim_index)
        init_subgrid_min(y_dim_index) = subgrid_min(y_dim_index)
        init_subgrid_max(x_dim_index) = subgrid_max(x_dim_index)
        init_subgrid_max(y_dim_index) = subgrid_max(y_dim_index)

    end subroutine uEMEP_set_subgrid_select_latlon_centre

end module set_subgrids
