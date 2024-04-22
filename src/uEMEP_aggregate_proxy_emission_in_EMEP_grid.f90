module aggregate_proxy_emission_in_emep_grid
    ! NOTE This module is not currently in use. Should it be deleted?

    use uEMEP_definitions
    use mod_read_esri_ascii_file, only: write_esri_ascii_file

    implicit none
    private

contains

    subroutine uEMEP_aggregate_proxy_emission_in_EMEP_grid()
        ! This routine takes subgrid emissions and aggregates them in the EMEP grid
        ! This is used for cross checking emissions
        integer :: i,j
        integer :: i_source
        real, allocatable :: EMEP_aggregated_subgid_emission(:,:)
        real, allocatable :: EMEP_aggregated_emission(:,:)
        integer, allocatable :: EMEP_aggregated_subgid_emission_count(:,:)
        real, allocatable :: lon_array(:,:), lat_array(:,:)
        integer :: iii, jjj
        character(256) :: temp_name
        real :: var3d_nc_local_temp
        integer :: t
        integer :: i_nc_source
        integer :: i_subsource=1
        integer :: i_pollutant, p_loop

        allocate(EMEP_aggregated_subgid_emission_count(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index)))
        allocate(EMEP_aggregated_subgid_emission(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index)))
        allocate(EMEP_aggregated_emission(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index)))
        allocate(lon_array(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index)))
        allocate(lat_array(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index)))

        EMEP_aggregated_subgid_emission = 0.0
        EMEP_aggregated_subgid_emission_count = 0

        t = 1
        unit_conversion = 1.0
        ! Units should be specified in the reading emission files

        ! Loop through the emission grid and aggregate data
        ! Units are the same as EMEP mg/m2
        do i_pollutant = 1, n_pollutant_loop
            p_loop = pollutant_loop_index(i_pollutant)
            do i_source = 1, n_source_index
                if (calculate_source(i_source)) then
                    if (make_EMEP_grid_emission_data(i_source)) then

                        write(unit_logfile,'(A)') ''
                        write(unit_logfile,'(A)') '================================================================'
                        write(unit_logfile,'(A)') 'Aggregating proxy emissions in EMEP grids for diagnostics (uEMEP_aggregate_proxy_emission_in_EMEP_grid)'
                        write(unit_logfile,'(A)') '================================================================'

                        write(unit_logfile,'(A)') 'Saving proxy subgrid emissions in EMEP grid ('//trim(source_file_str(i_source))//'_'//trim(pollutant_file_str(p_loop))//')'
                        if (i_source .eq. agriculture_index) then
                            i_nc_source = agriculture_nc_index
                        elseif (i_source .eq. traffic_index) then
                            i_nc_source = traffic_nc_index
                        elseif (i_source .eq. shipping_index) then
                            i_nc_source = shipping_nc_index
                        elseif (i_source .eq. heating_index) then
                            i_nc_source = heating_nc_index
                        else
                            write(unit_logfile,'(A)') 'Undefined source in routine uEMEP_aggregate_proxy_emission_in_EMEP_grid. Stopping'
                            stop 1
                        end if


                        EMEP_aggregated_subgid_emission = 0.0
                        EMEP_aggregated_emission = 0.0
                        EMEP_aggregated_subgid_emission_count = 0
                        do j = 1, emission_subgrid_dim(y_dim_index,i_source)
                            do i = 1, emission_subgrid_dim(x_dim_index,i_source)

                                ! Get indexes
                                iii = crossreference_emission_to_emep_subgrid(i,j,x_dim_index,i_source)
                                jjj = crossreference_emission_to_emep_subgrid(i,j,y_dim_index,i_source)
                                EMEP_aggregated_subgid_emission_count(iii,jjj) = EMEP_aggregated_subgid_emission_count(iii,jjj) + 1

                                ! Calculate as tonnes per year
                                EMEP_aggregated_subgid_emission(iii,jjj) = EMEP_aggregated_subgid_emission(iii,jjj) &
                                    + proxy_emission_subgrid(i,j,i_source,i_pollutant)*emission_factor_conversion(compound_index,i_source,i_subsource) &
                                    *1.0e-12*3600.0*24.0*365.0 ! Conversion from ug/sec/subgrid to ton/year
                                if (hourly_calculations) then
                                    EMEP_aggregated_emission(iii,jjj) = EMEP_aggregated_emission(iii,jjj) &
                                        + sum(var3d_nc(iii,jjj,:,emis_nc_index,allsource_index,i_pollutant))/dim_length_nc(time_dim_nc_index) &
                                        *(emission_subgrid_delta(y_dim_index,i_source)*emission_subgrid_delta(x_dim_index,i_source))*1.0e-9*24.0*365.0 ! Conversion from mg/m2/hr to ton/year (EMEP)
                                    ! NOTE: Have put in allsource_index because on this occasion there were no sectors in the EMEP file for this. Beware for later comparisons!!!
                                end if
                                if (annual_calculations) then
                                    EMEP_aggregated_emission(iii,jjj) = EMEP_aggregated_emission(iii,jjj)+sum(var3d_nc(iii,jjj,:,emis_nc_index,i_nc_source,i_pollutant))/dim_length_nc(time_dim_nc_index) &
                                        *(emission_subgrid_delta(y_dim_index,i_source)*emission_subgrid_delta(x_dim_index,i_source))*1.0e-9 ! Conversion from mg/m2/yr to ton/year (EMEP)
                                end if
                            end do
                        end do

                        do j = 1, dim_length_nc(y_dim_nc_index)
                            do i = 1, dim_length_nc(x_dim_nc_index)
                                if (EMEP_aggregated_subgid_emission_count(i,j) .gt. 0) then
                                    EMEP_aggregated_subgid_emission(i,j) = EMEP_aggregated_subgid_emission(i,j) !/EMEP_aggregated_subgid_emission_count(i,j)
                                    EMEP_aggregated_emission(i,j) = EMEP_aggregated_emission(i,j) !/EMEP_aggregated_subgid_emission_count(i,j)
                                else
                                    EMEP_aggregated_subgid_emission(i,j) = 0.0
                                    EMEP_aggregated_emission(i,j) = 0.0
                                end if
                                lon_array(i,j) = var1d_nc(i,x_dim_nc_index)
                                lat_array(i,j) = var1d_nc(j,y_dim_nc_index)
                            end do
                        end do

                        ! Check this. This is no longer valid and does not take into account neighbouring grids. Do not use
                        if (replace_EMEP_local_with_subgrid_local(i_source)) then
                            do j = 1, dim_length_nc(y_dim_nc_index)
                                do i = 1, dim_length_nc(x_dim_nc_index)
                                    if (EMEP_aggregated_subgid_emission(i,j) .ge. 0 .and. var3d_nc(i,j,t,emis_nc_index,i_nc_source,i_pollutant) .gt. 0) then
                                        var3d_nc_local_temp = var3d_nc(i,j,t,conc_nc_index,i_nc_source,i_pollutant)*(1.0 - var3d_nc(i,j,t,frac_nc_index,i_nc_source,i_pollutant)) ! nonlocal contribution
                                        var3d_nc(i,j,t,frac_nc_index,i_nc_source,i_pollutant) = EMEP_aggregated_subgid_emission(i,j)/var3d_nc(i,j,t,emis_nc_index,i_nc_source,i_pollutant)*var3d_nc(i,j,t,frac_nc_index,i_nc_source,i_pollutant) ! New local fraction
                                        var3d_nc(i,j,t,conc_nc_index,i_nc_source,i_pollutant) = var3d_nc(i,j,t,conc_nc_index,i_nc_source,i_pollutant)*var3d_nc(i,j,t,frac_nc_index,i_nc_source,i_pollutant) + var3d_nc_local_temp ! new total contribution
                                        var3d_nc(i,j,t,conc_nc_index,allsource_index,i_pollutant) = var3d_nc(i,j,t,conc_nc_index,i_nc_source,i_pollutant)
                                    end if
                                end do
                            end do
                        end if

                        temp_name = trim(pathname_grid(emission_file_index(i_source)))//trim(filename_grid(emission_file_index(i_source)))//'_'//trim(pollutant_file_str(p_loop))//'_aggregated_proxy_EMEP_'//trim(file_tag)//'.asc'
                        write(unit_logfile,'(a)') 'Writing to: '//trim(temp_name)
                        write(unit_logfile,'(a,f12.2)') 'Total local emissions (ton/year): ', sum(EMEP_aggregated_subgid_emission)
                        call write_esri_ascii_file(temp_name, dim_length_nc(x_dim_nc_index), dim_length_nc(y_dim_nc_index), dgrid_nc(lat_nc_index), EMEP_aggregated_subgid_emission(:,:), lon_array, lat_array)
                        
                        temp_name = trim(pathname_grid(emission_file_index(i_source)))//trim(filename_grid(emission_file_index(i_source)))//'_'//trim(pollutant_file_str(p_loop))//'_EMEP_'//trim(file_tag)//'.asc'
                        write(unit_logfile,'(a)') 'Writing to: '//trim(temp_name)
                        write(unit_logfile,'(a,f12.2)') 'Total EMEP emissions (ton/year): ', sum(EMEP_aggregated_emission)
                        call write_esri_ascii_file(temp_name,dim_length_nc(x_dim_nc_index), dim_length_nc(y_dim_nc_index), dgrid_nc(lat_nc_index), EMEP_aggregated_emission(:,:),lon_array,lat_array)

                        write(unit_logfile,'(a,f12.4)') 'Ratio of local to EMEP emissions: ', sum(EMEP_aggregated_subgid_emission)/sum(EMEP_aggregated_emission)
                    end if
                end if
            end do
        end do

        if (allocated(EMEP_aggregated_subgid_emission_count)) deallocate(EMEP_aggregated_subgid_emission_count)
        if (allocated(EMEP_aggregated_subgid_emission)) deallocate(EMEP_aggregated_subgid_emission)
        if (allocated(EMEP_aggregated_emission)) deallocate(EMEP_aggregated_emission)
        if (allocated(lon_array)) deallocate(lon_array)
        if (allocated(lat_array)) deallocate(lat_array)

    end subroutine uEMEP_aggregate_proxy_emission_in_EMEP_grid

end module aggregate_proxy_emission_in_emep_grid

