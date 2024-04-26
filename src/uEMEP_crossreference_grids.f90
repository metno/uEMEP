module crossreference_grids

    use uEMEP_definitions
    use mod_lambert_projection, only: lb2lambert2_uEMEP, LL2PS_spherical

    implicit none
    private

    public :: uEMEP_crossreference_grids, uEMEP_assign_region_coverage_to_EMEP

contains

    subroutine uEMEP_crossreference_grids()
        integer :: i, j, k
        integer :: ii, jj
        integer :: i_source
        real :: x_temp, y_temp

        ! Cross referencing must be done for each new grid when using multiple grids
        if (allocated(crossreference_target_to_emep_subgrid)) then
            deallocate (crossreference_target_to_emep_subgrid)
        end if
        if (allocated(crossreference_target_to_localfraction_subgrid)) then
            deallocate (crossreference_target_to_localfraction_subgrid)
        end if
        if (allocated(crossreference_integral_to_emep_subgrid)) then
            deallocate (crossreference_integral_to_emep_subgrid)
        end if
        if (allocated(crossreference_target_to_integral_subgrid)) then
            deallocate (crossreference_target_to_integral_subgrid)
        end if
        if (allocated(crossreference_target_to_emission_subgrid)) then
            deallocate (crossreference_target_to_emission_subgrid)
        end if
        if (allocated(crossreference_emission_to_EMEP_subgrid)) then
            deallocate (crossreference_emission_to_EMEP_subgrid)
        end if
        if (allocated(crossreference_integral_to_emission_subgrid)) then
            deallocate (crossreference_integral_to_emission_subgrid)
        end if
        if (allocated(crossreference_emission_to_integral_subgrid)) then
            deallocate (crossreference_emission_to_integral_subgrid)
        end if
        if (allocated(crossreference_target_to_population_subgrid)) then
            deallocate (crossreference_target_to_population_subgrid)
        end if
        if (use_alternative_meteorology_flag) then
            if (allocated(crossreference_integral_to_meteo_nc_subgrid)) then
                deallocate (crossreference_integral_to_meteo_nc_subgrid)
            end if
        end if
        if (calculate_deposition_flag) then
            if (allocated(crossreference_emission_to_deposition_subgrid)) then
                deallocate (crossreference_emission_to_deposition_subgrid)
            end if
            if (allocated(crossreference_target_to_deposition_subgrid)) then
                deallocate (crossreference_target_to_deposition_subgrid)
            end if
            if (allocated(crossreference_deposition_to_emep_subgrid)) then
                deallocate (crossreference_deposition_to_emep_subgrid)
            end if
        end if
        if (read_landuse_flag) then
            if (allocated(crossreference_emission_to_landuse_subgrid)) then
                deallocate (crossreference_emission_to_landuse_subgrid)
            end if
        end if

        ! Allocate arrays
        allocate (crossreference_target_to_emep_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),2))
        allocate (crossreference_target_to_localfraction_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),2,n_local_fraction_grids))
        allocate (crossreference_integral_to_emep_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),2))
        allocate (crossreference_target_to_integral_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),2))
        allocate (crossreference_target_to_emission_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),2,n_source_index))
        allocate (crossreference_emission_to_EMEP_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),2,n_source_index))
        allocate (crossreference_integral_to_emission_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),2,n_source_index))
        allocate (crossreference_emission_to_integral_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),2,n_source_index))
        allocate (crossreference_target_to_population_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),2))

        if (use_alternative_meteorology_flag) then
            allocate (crossreference_integral_to_meteo_nc_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),2))
        end if
        if (calculate_deposition_flag) then
            allocate (crossreference_emission_to_deposition_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),2,n_source_index))
            allocate (crossreference_target_to_deposition_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),2))
            allocate (crossreference_deposition_to_emep_subgrid(deposition_subgrid_dim(x_dim_index),deposition_subgrid_dim(y_dim_index),2))
        end if
        if (read_landuse_flag) then
            allocate (crossreference_emission_to_landuse_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),2,n_source_index))
        end if

        write(unit_logfile,'(A)') 'Allocating EMEP grid index to subgrid index'
        
        ! Loop through subgrid and find those subgrids within EMEP grids and allocate concentrations directly from EMEP grids.
        do j = 1, subgrid_dim(y_dim_index)
            do i = 1, subgrid_dim(x_dim_index)
                if (EMEP_projection_type .eq. LL_projection_index) then
                    ii = 1 + floor((lon_subgrid(i,j) - var1d_nc(1,lon_nc_index))/dgrid_nc(lon_nc_index) + 0.5)
                    jj = 1 + floor((lat_subgrid(i,j) - var1d_nc(1,lat_nc_index))/dgrid_nc(lat_nc_index) + 0.5)
                else if (EMEP_projection_type .eq. LCC_projection_index) then
                    ! When EMEP is read as x,y projection then var1d_nc(:,lon/lat_nc_index) are the x, y projection indexes, actually
                    call lb2lambert2_uEMEP(x_temp, y_temp, lon_subgrid(i,j), lat_subgrid(i,j), EMEP_projection_attributes)
                    ii = 1 + floor((x_temp - var1d_nc(1,lon_nc_index))/dgrid_nc(lon_nc_index) + 0.5)
                    jj = 1 + floor((y_temp - var1d_nc(1,lat_nc_index))/dgrid_nc(lat_nc_index) + 0.5)
                else if (EMEP_projection_type .eq. PS_projection_index) then
                    call LL2PS_spherical(x_temp,y_temp,lon_subgrid(i,j),lat_subgrid(i,j),EMEP_projection_attributes)
                    ii = 1 + floor((x_temp - var1d_nc(1,lon_nc_index))/dgrid_nc(lon_nc_index) + 0.5)
                    jj = 1 + floor((y_temp - var1d_nc(1,lat_nc_index))/dgrid_nc(lat_nc_index) + 0.5)
                else
                    write(unit_logfile,'(A)') 'No valid projection in use. Stopping'
                    stop 1
                end if
                crossreference_target_to_emep_subgrid(i,j,x_dim_index) = ii
                crossreference_target_to_emep_subgrid(i,j,y_dim_index) = jj
            end do
        end do

        write(unit_logfile,'(A)') 'Allocating EMEP local fraction grid index to subgrid index'
        ! Loop through subgrid and find those subgrids within EMEP grids and allocate concentrations directly from EMEP grids.
        do k = 1, n_local_fraction_grids
            do j = 1, subgrid_dim(y_dim_index)
                do i = 1, subgrid_dim(x_dim_index)
                    if (EMEP_projection_type .eq. LL_projection_index) then
                        ii = 1 + floor((lon_subgrid(i,j) - var1d_nc(1,lon_nc_index))/dgrid_nc(lon_nc_index)/local_fraction_grid_size(k) + 0.5)
                        jj = 1 + floor((lat_subgrid(i,j) - var1d_nc(1,lat_nc_index))/dgrid_nc(lat_nc_index)/local_fraction_grid_size(k) + 0.5)
                    else if (EMEP_projection_type .eq. LCC_projection_index) then
                        !When EMEP is read as x,y projection then var1d_nc(:,lon/lat_nc_index) are the x, y projection indexes, actually
                        call lb2lambert2_uEMEP(x_temp, y_temp, lon_subgrid(i,j), lat_subgrid(i,j), EMEP_projection_attributes)
                        ii = 1 + floor((x_temp - var1d_nc(1,lon_nc_index))/dgrid_nc(lon_nc_index)/local_fraction_grid_size(k) + 0.5)
                        jj = 1 + floor((y_temp - var1d_nc(1,lat_nc_index))/dgrid_nc(lat_nc_index)/local_fraction_grid_size(k) + 0.5)
                    else if (EMEP_projection_type .eq. PS_projection_index) then
                        call LL2PS_spherical(x_temp, y_temp, lon_subgrid(i,j), lat_subgrid(i,j),EMEP_projection_attributes)
                        ii = 1 + floor((x_temp - var1d_nc(1,lon_nc_index))/dgrid_nc(lon_nc_index)/local_fraction_grid_size(k) + 0.5)
                        jj = 1 + floor((y_temp - var1d_nc(1,lat_nc_index))/dgrid_nc(lat_nc_index)/local_fraction_grid_size(k) + 0.5)
                    else
                        write(unit_logfile,'(A)') 'No valid projection in use. Stopping'
                        stop 1
                    end if
                    crossreference_target_to_localfraction_subgrid(i,j,x_dim_index,k) = ii
                    crossreference_target_to_localfraction_subgrid(i,j,y_dim_index,k) = jj
                end do
            end do
        end do

        write(unit_logfile,'(A)') 'Allocating integral grid index to subgrid index'
        do j = 1, subgrid_dim(y_dim_index)
            do i = 1, subgrid_dim(x_dim_index)
                crossreference_target_to_integral_subgrid(i,j,x_dim_index) = 1 + floor((x_subgrid(i,j) - integral_subgrid_min(x_dim_index))/integral_subgrid_delta(x_dim_index))
                crossreference_target_to_integral_subgrid(i,j,y_dim_index) = 1 + floor((y_subgrid(i,j) - integral_subgrid_min(y_dim_index))/integral_subgrid_delta(y_dim_index))
            end do
        end do
        write(unit_logfile,'(A)') 'Allocating population grid index to subgrid index'
        do j = 1, subgrid_dim(y_dim_index)
            do i = 1, subgrid_dim(x_dim_index)
                crossreference_target_to_population_subgrid(i,j,x_dim_index) = 1 + floor((x_subgrid(i,j) - population_subgrid_min(x_dim_index))/population_subgrid_delta(x_dim_index))
                crossreference_target_to_population_subgrid(i,j,y_dim_index) = 1 + floor((y_subgrid(i,j) - population_subgrid_min(y_dim_index))/population_subgrid_delta(y_dim_index))
            end do
        end do

        write(unit_logfile,'(A)') 'Allocating EMEP grid index to integral subgrid index'
        do j = 1, integral_subgrid_dim(y_dim_index)
            do i = 1, integral_subgrid_dim(x_dim_index)
                if (EMEP_projection_type .eq. LL_projection_index) then
                    ii = 1 + floor((lon_integral_subgrid(i,j) - var1d_nc(1,lon_nc_index))/dgrid_nc(lon_nc_index) + 0.5)
                    jj = 1 + floor((lat_integral_subgrid(i,j) - var1d_nc(1,lat_nc_index))/dgrid_nc(lat_nc_index) + 0.5)
                else if (EMEP_projection_type .eq. LCC_projection_index) then
                    call lb2lambert2_uEMEP(x_temp, y_temp, lon_integral_subgrid(i,j), lat_integral_subgrid(i,j), EMEP_projection_attributes)
                    ii = 1 + floor((x_temp - var1d_nc(1,lon_nc_index))/dgrid_nc(lon_nc_index) + 0.5)
                    jj = 1 + floor((y_temp - var1d_nc(1,lat_nc_index))/dgrid_nc(lat_nc_index) + 0.5)
                else if (EMEP_projection_type .eq. PS_projection_index) then
                    call LL2PS_spherical(x_temp, y_temp, lon_integral_subgrid(i,j), lat_integral_subgrid(i,j), EMEP_projection_attributes)
                    ii = 1 + floor((x_temp - var1d_nc(1,lon_nc_index))/dgrid_nc(lon_nc_index) + 0.5)
                    jj = 1 + floor((y_temp - var1d_nc(1,lat_nc_index))/dgrid_nc(lat_nc_index) + 0.5)
                else
                    write(unit_logfile,'(A)') 'No valid projection in use. Stopping'
                    stop 1
                end if
                crossreference_integral_to_emep_subgrid(i,j,x_dim_index) = ii
                crossreference_integral_to_emep_subgrid(i,j,y_dim_index) = jj
            end do
        end do

        if (use_alternative_meteorology_flag) then
            write(unit_logfile,'(A)') 'Allocating alternative meteo nc grid index to integral subgrid index'
            do j = 1, integral_subgrid_dim(y_dim_index)
                do i = 1, integral_subgrid_dim(x_dim_index)
                    if (meteo_nc_projection_type .eq. LL_projection_index) then
                        ii = 1 + floor((lon_integral_subgrid(i,j) - meteo_var1d_nc(1,lon_nc_index))/meteo_dgrid_nc(lon_nc_index) + 0.5)
                        jj = 1 + floor((lat_integral_subgrid(i,j) - meteo_var1d_nc(1,lat_nc_index))/meteo_dgrid_nc(lat_nc_index) + 0.5)
                    else if (meteo_nc_projection_type .eq. LCC_projection_index) then
                        call lb2lambert2_uEMEP(x_temp, y_temp, lon_integral_subgrid(i,j), lat_integral_subgrid(i,j), meteo_nc_projection_attributes)
                        ii = 1 + floor((x_temp - meteo_var1d_nc(1,lon_nc_index))/meteo_dgrid_nc(lon_nc_index) + 0.5)
                        jj = 1 + floor((y_temp - meteo_var1d_nc(1,lat_nc_index))/meteo_dgrid_nc(lat_nc_index) + 0.5)
                    else if (meteo_nc_projection_type .eq. PS_projection_index) then
                        call LL2PS_spherical(x_temp, y_temp, lon_integral_subgrid(i,j), lat_integral_subgrid(i,j), meteo_nc_projection_attributes)
                        ii = 1 + floor((x_temp - meteo_var1d_nc(1,lon_nc_index))/meteo_dgrid_nc(lon_nc_index) + 0.5)
                        jj = 1 + floor((y_temp - meteo_var1d_nc(1,lat_nc_index))/meteo_dgrid_nc(lat_nc_index) + 0.5)
                    else
                        write(unit_logfile,'(A)')'No valid projection in use. Stopping'
                        stop 1
                    end if
                    crossreference_integral_to_meteo_nc_subgrid(i,j,x_dim_index) = ii
                    crossreference_integral_to_meteo_nc_subgrid(i,j,y_dim_index) = jj
                end do
            end do
        end if

        do i_source = 1, n_source_index
            if (calculate_source(i_source)) then
                do j = 1, emission_subgrid_dim(y_dim_index,i_source)
                    do i = 1, emission_subgrid_dim(x_dim_index,i_source)
                        if (EMEP_projection_type .eq. LL_projection_index) then
                            ii = 1 + floor((lon_emission_subgrid(i,j,i_source) - var1d_nc(1,lon_nc_index))/dgrid_nc(lon_nc_index) + 0.5)
                            jj = 1 + floor((lat_emission_subgrid(i,j,i_source) - var1d_nc(1,lat_nc_index))/dgrid_nc(lat_nc_index) + 0.5)
                        else if (EMEP_projection_type .eq. LCC_projection_index) then
                            call lb2lambert2_uEMEP(x_temp, y_temp, lon_emission_subgrid(i,j,i_source), lat_emission_subgrid(i,j,i_source), EMEP_projection_attributes)
                            ii = 1 + floor((x_temp - var1d_nc(1,lon_nc_index))/dgrid_nc(lon_nc_index) + 0.5)
                            jj = 1 + floor((y_temp-var1d_nc(1,lat_nc_index))/dgrid_nc(lat_nc_index) + 0.5)
                        else if (EMEP_projection_type .eq. PS_projection_index) then
                            call LL2PS_spherical(x_temp, y_temp, lon_emission_subgrid(i,j,i_source), lat_emission_subgrid(i,j,i_source), EMEP_projection_attributes)
                            ii = 1 + floor((x_temp - var1d_nc(1,lon_nc_index))/dgrid_nc(lon_nc_index) + 0.5)
                            jj = 1 + floor((y_temp - var1d_nc(1,lat_nc_index))/dgrid_nc(lat_nc_index) + 0.5)
                        else
                            write(unit_logfile,'(A)') 'No valid projection in use. Stopping'
                            stop 1
                        end if
                        crossreference_emission_to_emep_subgrid(i,j,x_dim_index,i_source) = ii
                        crossreference_emission_to_emep_subgrid(i,j,y_dim_index,i_source) = jj
                    end do
                end do
                do j = 1, subgrid_dim(y_dim_index)
                    do i = 1, subgrid_dim(x_dim_index)
                        crossreference_target_to_emission_subgrid(i,j,x_dim_index,i_source) = 1 + floor((x_subgrid(i,j) - emission_subgrid_min(x_dim_index,i_source))/emission_subgrid_delta(x_dim_index,i_source))
                        crossreference_target_to_emission_subgrid(i,j,y_dim_index,i_source) = 1 + floor((y_subgrid(i,j) - emission_subgrid_min(y_dim_index,i_source))/emission_subgrid_delta(y_dim_index,i_source))
                    end do
                end do
                do j = 1, emission_subgrid_dim(y_dim_index,i_source)
                    do i = 1, emission_subgrid_dim(x_dim_index,i_source)
                        crossreference_emission_to_integral_subgrid(i,j,x_dim_index,i_source) = 1 + floor((x_emission_subgrid(i,j,i_source) - integral_subgrid_min(x_dim_index))/integral_subgrid_delta(x_dim_index))
                        crossreference_emission_to_integral_subgrid(i,j,y_dim_index,i_source) = 1 + floor((y_emission_subgrid(i,j,i_source) - integral_subgrid_min(y_dim_index))/integral_subgrid_delta(y_dim_index))

                        ! At edge this can return negative distances due to the different sizes of emission and integral grids and buffer zones. Set the limits here. Should not be a problem
                        crossreference_emission_to_integral_subgrid(i,j,x_dim_index,i_source) = max(min(crossreference_emission_to_integral_subgrid(i,j,x_dim_index,i_source), integral_subgrid_dim(x_dim_index)), 1)
                        crossreference_emission_to_integral_subgrid(i,j,y_dim_index,i_source) = max(min(crossreference_emission_to_integral_subgrid(i,j,y_dim_index,i_source), integral_subgrid_dim(y_dim_index)), 1)

                        if (crossreference_emission_to_integral_subgrid(i,j,x_dim_index,i_source) .lt. 1 .or. crossreference_emission_to_integral_subgrid(i,j,x_dim_index,i_source) .gt. integral_subgrid_dim(x_dim_index) &
                            .or. crossreference_emission_to_integral_subgrid(i,j,y_dim_index,i_source) .lt. 1 .or. crossreference_emission_to_integral_subgrid(i,j,y_dim_index,i_source) .gt. integral_subgrid_dim(y_dim_index)) then

                            write(unit_logfile,'(A,4i,4f)') 'WARNING: crossreference_emission_to_integral_subgrid is out of bounds (i_emis,j_emis,i_integral,j_integral,x_emis,y_emis)', i, j, &
                                crossreference_emission_to_integral_subgrid(i,j,x_dim_index,i_source), crossreference_emission_to_integral_subgrid(i,j,y_dim_index,i_source), &
                                x_emission_subgrid(i,j,i_source)/1000, y_emission_subgrid(i,j,i_source)/1000, (x_emission_subgrid(i,j,i_source) - integral_subgrid_min(x_dim_index))/integral_subgrid_delta(x_dim_index) + 0.5, &
                                (y_emission_subgrid(i,j,i_source) - integral_subgrid_min(y_dim_index))/integral_subgrid_delta(y_dim_index) + 0.5
                        end if
                    end do
                end do
                do j = 1, integral_subgrid_dim(y_dim_index)
                    do i = 1, integral_subgrid_dim(x_dim_index)
                        crossreference_integral_to_emission_subgrid(i,j,x_dim_index,i_source) = 1 + floor((x_integral_subgrid(i,j) - emission_subgrid_min(x_dim_index,i_source))/emission_subgrid_delta(x_dim_index,i_source))
                        crossreference_integral_to_emission_subgrid(i,j,y_dim_index,i_source) = 1 + floor((y_integral_subgrid(i,j) - emission_subgrid_min(y_dim_index,i_source))/emission_subgrid_delta(y_dim_index,i_source))
                    end do
                end do

                if (calculate_deposition_flag) then
                    do j = 1, emission_subgrid_dim(y_dim_index,i_source)
                        do i = 1, emission_subgrid_dim(x_dim_index,i_source)
                            crossreference_emission_to_deposition_subgrid(i,j,x_dim_index,i_source) = 1 + floor((x_emission_subgrid(i,j,i_source)-deposition_subgrid_min(x_dim_index))/deposition_subgrid_delta(x_dim_index))
                            crossreference_emission_to_deposition_subgrid(i,j,y_dim_index,i_source) = 1 + floor((y_emission_subgrid(i,j,i_source)-deposition_subgrid_min(y_dim_index))/deposition_subgrid_delta(y_dim_index))
                            
                            ! At edge this can return negative distances due to the different sizes of emission and integral grids and buffer zones. Set the limits here. Should not be a problem
                            crossreference_emission_to_deposition_subgrid(i,j,x_dim_index,i_source) = max(min(crossreference_emission_to_deposition_subgrid(i,j,x_dim_index,i_source), deposition_subgrid_dim(x_dim_index)), 1)
                            crossreference_emission_to_deposition_subgrid(i,j,y_dim_index,i_source) = max(min(crossreference_emission_to_deposition_subgrid(i,j,y_dim_index,i_source), deposition_subgrid_dim(y_dim_index)), 1)
                        end do
                    end do
                end if

                if (read_landuse_flag) then
                    do j = 1, emission_subgrid_dim(y_dim_index,i_source)
                        do i = 1, emission_subgrid_dim(x_dim_index,i_source)
                            crossreference_emission_to_landuse_subgrid(i,j,x_dim_index,i_source) = 1 + floor((x_emission_subgrid(i,j,i_source)-landuse_subgrid_min(x_dim_index))/landuse_subgrid_delta(x_dim_index))
                            crossreference_emission_to_landuse_subgrid(i,j,y_dim_index,i_source) = 1 + floor((y_emission_subgrid(i,j,i_source)-landuse_subgrid_min(y_dim_index))/landuse_subgrid_delta(y_dim_index))
                            
                            ! At edge this can return negative distances due to the different sizes of emission and integral grids and buffer zones. Set the limits here. Should not be a problem
                            crossreference_emission_to_landuse_subgrid(i,j,x_dim_index,i_source) = max(min(crossreference_emission_to_landuse_subgrid(i,j,x_dim_index,i_source), landuse_subgrid_dim(x_dim_index)), 1)
                            crossreference_emission_to_landuse_subgrid(i,j,y_dim_index,i_source) = max(min(crossreference_emission_to_landuse_subgrid(i,j,y_dim_index,i_source), landuse_subgrid_dim(y_dim_index)), 1)
                        end do
                    end do
                end if
            end if
        end do

        if (calculate_deposition_flag) then
            do j = 1, subgrid_dim(y_dim_index)
                do i = 1, subgrid_dim(x_dim_index)
                    crossreference_target_to_deposition_subgrid(i,j,x_dim_index) = 1 + floor((x_subgrid(i,j)-deposition_subgrid_min(x_dim_index))/deposition_subgrid_delta(x_dim_index))
                    crossreference_target_to_deposition_subgrid(i,j,y_dim_index) = 1 + floor((y_subgrid(i,j)-deposition_subgrid_min(y_dim_index))/deposition_subgrid_delta(y_dim_index))
                    
                    ! At edge this can return negative distances due to the different sizes of emission and integral grids and buffer zones. Set the limits here. Should not be a problem
                    crossreference_target_to_deposition_subgrid(i,j,x_dim_index) = max(min(crossreference_target_to_deposition_subgrid(i,j,x_dim_index), deposition_subgrid_dim(x_dim_index)), 1)
                    crossreference_target_to_deposition_subgrid(i,j,y_dim_index) = max(min(crossreference_target_to_deposition_subgrid(i,j,y_dim_index), deposition_subgrid_dim(y_dim_index)), 1)
                end do
            end do

            write(unit_logfile,'(A)') 'Allocating EMEP grid index to deposition subgrid index'
            do j = 1, deposition_subgrid_dim(y_dim_index)
                do i = 1, deposition_subgrid_dim(x_dim_index)
                    if (EMEP_projection_type .eq. LL_projection_index) then
                        ii = 1 + floor((lon_deposition_subgrid(i,j) - var1d_nc(1,lon_nc_index))/dgrid_nc(lon_nc_index) + 0.5)
                        jj = 1 + floor((lat_deposition_subgrid(i,j) - var1d_nc(1,lat_nc_index))/dgrid_nc(lat_nc_index) + 0.5)
                    else if (EMEP_projection_type .eq. LCC_projection_index) then
                        call lb2lambert2_uEMEP(x_temp, y_temp, lon_deposition_subgrid(i,j), lat_deposition_subgrid(i,j), EMEP_projection_attributes)
                        ii = 1 + floor((x_temp - var1d_nc(1,lon_nc_index))/dgrid_nc(lon_nc_index) + 0.5)
                        jj = 1 + floor((y_temp - var1d_nc(1,lat_nc_index))/dgrid_nc(lat_nc_index) + 0.5)
                    else if (EMEP_projection_type .eq. PS_projection_index) then
                        call LL2PS_spherical(x_temp, y_temp, lon_deposition_subgrid(i,j), lat_deposition_subgrid(i,j), EMEP_projection_attributes)
                        ii = 1 + floor((x_temp - var1d_nc(1,lon_nc_index))/dgrid_nc(lon_nc_index) + 0.5)
                        jj = 1 + floor((y_temp - var1d_nc(1,lat_nc_index))/dgrid_nc(lat_nc_index) + 0.5)
                    else
                        write(unit_logfile,'(A)') 'No valid projection in use. Stopping'
                        stop 1
                    end if
                    crossreference_deposition_to_emep_subgrid(i,j,x_dim_index) = ii
                    crossreference_deposition_to_emep_subgrid(i,j,y_dim_index) = jj
                end do
            end do
        end if

    end subroutine uEMEP_crossreference_grids

    subroutine uEMEP_assign_region_coverage_to_EMEP()
        integer :: i, j
        integer :: ii, jj, iii, jjj, iiii, jjjj, ii_nc, jj_nc
        real, allocatable :: count_EMEP_grid_fraction_in_region(:,:,:)
        integer :: i_source, chosen_source
        integer :: ii_start, jj_start

        write(unit_logfile,'(A)') ''
        write(unit_logfile,'(A)') '================================================================'
        write(unit_logfile,'(A)') 'Assigning regional coverage to EMEP grids (uEMEP_assign_region_coverage_to_EMEP)'
        write(unit_logfile,'(A)') '================================================================'

        ! Allocate the fraction of each EMEP grid that is within a region, if required
        if (trace_emissions_from_in_region) then
            if (allocated(EMEP_grid_fraction_in_region)) then
                deallocate (EMEP_grid_fraction_in_region)
            end if
            if ( .not. allocated(EMEP_grid_fraction_in_region)) then
                allocate (EMEP_grid_fraction_in_region(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),n_source_nc_index,2))
            end if
            if (allocated(lf_EMEP_grid_fraction_in_region)) then
                deallocate (lf_EMEP_grid_fraction_in_region)
            end if
            if ( .not. allocated(lf_EMEP_grid_fraction_in_region)) then
                allocate (lf_EMEP_grid_fraction_in_region(dim_length_nc(xdist_dim_nc_index),dim_length_nc(ydist_dim_nc_index),dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),n_source_nc_index,2))
            end if
            if (allocated(count_EMEP_grid_fraction_in_region)) then
                deallocate (count_EMEP_grid_fraction_in_region)
            end if
            if ( .not. allocated(count_EMEP_grid_fraction_in_region)) then
                allocate (count_EMEP_grid_fraction_in_region(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),n_source_nc_index))
            end if
            
            EMEP_grid_fraction_in_region = 0.0
            lf_EMEP_grid_fraction_in_region = 0.0
            count_EMEP_grid_fraction_in_region = 0.0

            do i_source = 1, n_source_index
                if (calculate_source(i_source) .or. calculate_EMEP_source(i_source)) then
                    ! Loop through the subgrid and find the EMEP grid it is in
                    do j = 1, emission_subgrid_dim(y_dim_index,i_source)
                        do i = 1, emission_subgrid_dim(x_dim_index,i_source)
                            ii = crossreference_emission_to_emep_subgrid(i,j,x_dim_index,i_source)
                            jj = crossreference_emission_to_emep_subgrid(i,j,y_dim_index,i_source)
                            if (ii .ge. 1 .and. ii .le. dim_length_nc(x_dim_nc_index) .and. jj .ge. 1 .and. jj .le. dim_length_nc(y_dim_nc_index)) then
                                if (use_subgrid_region(i,j,i_source)) then
                                    EMEP_grid_fraction_in_region(ii,jj,i_source,1) = EMEP_grid_fraction_in_region(ii,jj,i_source,1) + 1.0
                                end if
                                count_EMEP_grid_fraction_in_region(ii,jj,i_source) = count_EMEP_grid_fraction_in_region(ii,jj,i_source) + 1.0
                            end if
                        end do
                    end do

                    where (count_EMEP_grid_fraction_in_region(:,:,i_source) .gt. 0) EMEP_grid_fraction_in_region(:,:,i_source,1) = EMEP_grid_fraction_in_region(:,:,i_source,1)/count_EMEP_grid_fraction_in_region(:,:,i_source)

                    ! Loop through all the EMEP grids and fill in the local fraction arrays
                    do ii = 1, dim_length_nc(x_dim_nc_index)
                        do jj = 1, dim_length_nc(y_dim_nc_index)
                            do iii = 1, dim_length_nc(xdist_dim_nc_index)
                                do jjj = 1, dim_length_nc(ydist_dim_nc_index)
                                    iiii = iii - xdist_centre_nc
                                    jjjj = jjj - ydist_centre_nc
                                    ii_nc = ii + iiii
                                    jj_nc = jj + jjjj
                                    if (ii_nc .ge. 1 .and. ii_nc .le. dim_length_nc(x_dim_nc_index) .and. jj_nc .ge. 1 .and. jj_nc .le. dim_length_nc(y_dim_nc_index)) then
                                        lf_EMEP_grid_fraction_in_region(iii,jjj,ii,jj,:,1) = EMEP_grid_fraction_in_region(ii_nc,jj_nc,:,1)
                                    end if
                                end do
                            end do
                        end do
                    end do

                    ! Allocate the additional subgrid fraction based on the EMEP fraction
                    if (EMEP_additional_grid_interpolation_size .gt. 0) then
                        ! Use the starting position of the read in EMEP file to initialise the starting point
                        ii_start = mod(dim_start_EMEP_nc(x_dim_nc_index)-1,local_fraction_grid_size(2))
                        jj_start = mod(dim_start_EMEP_nc(y_dim_nc_index)-1,local_fraction_grid_size(2))
                        do ii = 1, dim_length_nc(x_dim_nc_index)
                            do jj = 1, dim_length_nc(y_dim_nc_index)
                                ! Bottom left corner of the additional grid associatedd with thatEMEP grid
                                iii = int((ii-1+ii_start)/local_fraction_grid_size(2))*local_fraction_grid_size(2) + 1 - ii_start
                                jjj = int((jj-1+jj_start)/local_fraction_grid_size(2))*local_fraction_grid_size(2) + 1 - jj_start
                                if (iii .ge. 1 .and. iii .le. dim_length_nc(x_dim_nc_index) .and. jjj .ge. 1 .and. jjj .le. dim_length_nc(y_dim_nc_index)) then
                                    ! Add the region fractions from EMEP to the additional
                                    do iiii = iii, iii-1+local_fraction_grid_size(2)
                                        do jjjj = jjj, jjj-1+local_fraction_grid_size(2)
                                            if (iiii .ge. 1 .and. iiii .le. dim_length_nc(x_dim_nc_index) .and. jjjj .ge. 1 .and. jjjj .le. dim_length_nc(y_dim_nc_index)) then
                                                EMEP_grid_fraction_in_region(ii,jj,i_source,2) = EMEP_grid_fraction_in_region(ii,jj,i_source,2) + EMEP_grid_fraction_in_region(iiii,jjjj,i_source,1)
                                            end if
                                        end do
                                    end do
                                end if
                                !Normalize. Does not account for edges but should not be a problem. For safety limit it to 1.
                                EMEP_grid_fraction_in_region(ii,jj,i_source,2) = min(1.0, EMEP_grid_fraction_in_region(ii,jj,i_source,2)*local_fraction_grid_size(1)**2/local_fraction_grid_size(2)**2)
                            end do
                        end do

                        ! When using additional grid
                        ! Loop through all the additional EMEP grids and fill in the local fraction arrays
                        do ii = 1, dim_length_nc(x_dim_nc_index)
                            do jj = 1, dim_length_nc(y_dim_nc_index)
                                do iii = 1, dim_length_nc(xdist_dim_nc_index)
                                    do jjj = 1, dim_length_nc(ydist_dim_nc_index)
                                        iiii = iii - xdist_centre_nc
                                        jjjj = jjj - ydist_centre_nc
                                        ii_nc = ii + iiii*local_fraction_grid_size(2)
                                        jj_nc = jj + jjjj*local_fraction_grid_size(2)
                                        if (ii_nc .ge. 1 .and. ii_nc .le. dim_length_nc(x_dim_nc_index) .and. jj_nc .ge. 1 .and. jj_nc .le. dim_length_nc(y_dim_nc_index)) then
                                            lf_EMEP_grid_fraction_in_region(iii,jjj,ii,jj,:,2) = EMEP_grid_fraction_in_region(ii_nc,jj_nc,:,2)
                                        end if
                                    end do
                                end do
                            end do
                        end do
                    end if
                end if
            end do

            if (allocated(count_EMEP_grid_fraction_in_region)) deallocate (count_EMEP_grid_fraction_in_region)

            ! Fill in missing sources for EMEP
            ! Choose a source with highest resolution, the last one in this case as they all should be the same
            chosen_source = 0
            do i_source = 1, n_source_index
                if (calculate_source(i_source)) then
                    chosen_source = i_source
                end if
            end do
            do i_source = 1, n_source_nc_index
                if (chosen_source .ne. i_source) then
                    EMEP_grid_fraction_in_region(:,:,i_source,:) = EMEP_grid_fraction_in_region(:,:,chosen_source,:)
                    lf_EMEP_grid_fraction_in_region(:,:,:,:,i_source,:) = lf_EMEP_grid_fraction_in_region(:,:,:,:,chosen_source,:)
                end if
            end do
        end if

    end subroutine uEMEP_assign_region_coverage_to_EMEP

end module crossreference_grids

