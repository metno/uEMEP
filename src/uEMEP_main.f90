program uEMEP
    !! ****************************************************************************
    !!   uEMEP
    !!
    !!   Bruce rolstad Denby (brucerd@met.no)
    !!   MET Norway
    !! ****************************************************************************
    !! Reminder note for compilation on Intel
    !! To add this library to your linker input in the IDE, open the context menu for the project node, choose Properties, then in the Project Properties dialog box, choose Linker
    !! , and edit the Linker Input to add legacy_stdio_definitions.lib to the semi-colon-separated list
    !! Tools/options/intel compilers and tools/visual fortran/compilers and add bin, include and lib, e.g. C:\Program Files (x86)\netcdf 4.3.3.1\bin;
    !!   Control programme for running the downscaling routine uEMEP
    !! ****************************************************************************
    !!
    !! ****************************************************************************
    !! To link to netcdf in visual studio
    !! Tools - options - Intel compilers - VisusalFortran - Compilers - Libraries/includes/executables
    !! C:\Program Files (x86)\netcdf 4.3.3.1\include
    !! C:\Program Files (x86)\netcdf 4.3.3.1\bin
    !! C:\Program Files (x86)\netcdf 4.3.3.1\lib
    
    use uemep_configuration
    use uEMEP_definitions
    use read_command_line, only: uEMEP_read_command_line, check_command_line
    use set_constants, only: uEMEP_set_constants, uEMEP_set_pollutant_loop, &
        uEMEP_reset_constants, uEMEP_set_species_loop
    use read_config, only: uEMEP_read_config
    use read_emep, only: uEMEP_read_EMEP
    use read_meteo_nc, only: uEMEP_read_meteo_nc
    use read_rwc_heating_data, only: uEMEP_read_RWC_heating_data
    use save_emission_netcdf, only: uEMEP_calculate_emissions_for_EMEP
    use set_subgrids, only: uEMEP_set_subgrids, uEMEP_set_subgrid_select_latlon_centre
    use read_landuse_rivm_data, only: uEMEP_read_landuse_rivm_data, &
        uEMEP_set_landuse_classes, uEMEP_read_netcdf_landuse_latlon
    use read_roadlink_data_ascii, only: read_country_bounding_box_data, &
        uEMEP_read_roadlink_data_ascii, uEMEP_change_road_data, uEMEP_read_roadlink_emission_data
    use set_filenames, only: uEMEP_set_filenames
    use read_receptor_data, only: uEMEP_read_receptor_data, uEMEP_set_loop_receptor_grid, &
        uEMEP_grid_receptor_data
    use read_ssb_data, only: uEMEP_read_netcdf_population, uEMEP_read_SSB_data, &
        uEMEP_read_netcdf_population_latlon
    use read_agriculture_asi_data, only: uEMEP_read_agriculture_rivm_data, &
        uEMEP_read_emission_rivm_data
    use read_industry_data, only: uEMEP_read_industry_data
    use read_shipping_asi_data, only: uEMEP_preaggregate_shipping_asi_data, &
        uEMEP_read_netcdf_shipping_latlon, uEMEP_read_weekly_shipping_asi_data, &
        uEMEP_read_monthly_and_daily_shipping_asi_data, uEMEP_read_shipping_asi_data
    use read_time_profiles, only: uEMEP_read_time_profiles
    use redistribute_data, only: uEMEP_redistribute_local_source, uEMEP_disperse_local_source, &
        uEMEP_combine_local_source
    use save_netcdf_file, only: uEMEP_save_netcdf_control
    use subgrid_deposition, only: uEMEP_subgrid_deposition
    use subgrid_dispersion, only: uEMEP_subgrid_dispersion
    use set_emission_factors, only: uEMEP_set_emission_factors, uEMEP_convert_proxy_to_emissions, &
        uEMEP_nox_emission_temperature
    use subgrid_emep, only: uEMEP_subgrid_EMEP, uEMEP_subgrid_EMEP_from_in_region
    use subgrid_deposition_emep, only: uEMEP_set_deposition_velocities, &
        uEMEP_subgrid_deposition_EMEP, uEMEP_calculate_deposition
    use subgrid_emission_emep, only: uEMEP_subgrid_emission_EMEP
    use subgrid_meteo_emep, only: uEMEP_subgrid_meteo_EMEP
    use tiling_routines, only: uEMEP_set_tile_grids, uEMEP_set_region_tile_grids
    use chemistry_no2, only: uEMEP_chemistry, correct_annual_mean_chemistry
    use crossreference_grids, only: uEMEP_crossreference_grids
    use grid_roads, only: uEMEP_grid_roads
    use define_subgrid, only: uEMEP_define_subgrid_extent, uEMEP_define_subgrid
    use calculate_exposure, only: uEMEP_calculate_exposure
    use auto_subgrid, only: uEMEP_region_mask_new

    use uemep_logger

    implicit none

    integer :: source_index
    real :: start_time_cpu, end_time_cpu
    logical :: have_read_emep = .false.
    logical :: use_default_config = .false.

    ! Start timer
    call cpu_Time(start_time_cpu)

        ! Set model version
    model_version_str='7.0.7'

    ! Check command line arguments and handle special cases that have to be printed to stdout
    call check_command_line(use_default_config)

    write(*,*) ''
    write(*,*) '------------------------------------------------------------------------'
    write(*,*) 'Starting program uEMEP v'//trim(model_version_str)
    write(*,*) '------------------------------------------------------------------------'

    ! Read the command line, assigning the configuration file names and the substitution date_str
    call uEMEP_read_command_line(use_default_config)

    ! Set constants and variable names to be read from EMEP and meteo files
    call uEMEP_set_constants()

    ! Read the configuration files. Hard coded to be up to 5 files. Log file opened in this routine
    call uEMEP_read_config()

    ! If selected then specify subgrid using the lat and lon coordinates
    if (select_latlon_centre_domain_position_flag) then
        call uEMEP_set_subgrid_select_latlon_centre()
    end if

    ! Set the landuse if required
    if (use_landuse_as_proxy .or. read_landuse_flag) then
        call uEMEP_set_landuse_classes()
    end if

    ! Set the pollutant and compound loop definitions
    call uEMEP_set_pollutant_loop()

    !Reset any constants needed based on the configuration input
    call uEMEP_reset_constants()

    ! Autoselect files and countries if required. Place here because it changes config data
    if (auto_select_OSM_country_flag .or. trim(select_country_by_name) .ne. '') then
        call read_country_bounding_box_data()
    end if

    ! Set the EMEP species definitions if they are to be read
    call uEMEP_set_species_loop()

    ! Set the names of files to be written to when saving intermediate files
    call uEMEP_set_filenames()

    ! Read positions of receptor points (usually observations) for specifying multiple receptor grids or calculation points within a single grid
    call uEMEP_read_receptor_data()

    ! Enter the routine for saving emissions used in uEMEP for EMEP in netcdf files defined for the Norwegian domain in Lambert coordinates. Will stop after this
    if (save_emissions_for_EMEP(allsource_index)) then
        call uEMEP_calculate_emissions_for_EMEP()
    end if

    ! We set up an initial emission grid parameter set that can be used to first select the outester region
    ! This has been done to enable reading of multiple road link files but only keeping those in the initial defined emission area
    call uEMEP_set_subgrids()
    init_emission_subgrid_min = emission_subgrid_min
    init_emission_subgrid_max = emission_subgrid_max
    init_emission_subgrid_dim = emission_subgrid_dim
    init_emission_subgrid_delta = emission_subgrid_delta

    ! Set the grid loop (g_loop) extent based on use_multiple_receptor_grids_flag or not
    if (use_multiple_receptor_grids_flag) then
        start_grid_loop_index = 1
        end_grid_loop_index = n_receptor_in
        n_receptor = 1
        n_valid_receptor = 1
        valid_receptor_index(1) = 1
    else
        start_grid_loop_index = 1
        end_grid_loop_index = 1
        n_receptor = n_receptor_in
        use_receptor(start_grid_loop_index) = .true.
    end if

    first_g_loop = .true.

    ! If the use_single_time_loop_flag is true (Reads and calculates one time step at a time to save memory) then set these parameters
    if (use_single_time_loop_flag) then
        start_time_loop_index = 1
        end_time_loop_index = end_time_nc_index - start_time_nc_index + 1
        subgrid_dim(t_dim_index) = 1
        dim_length_nc(time_dim_nc_index) = 1
    else
        start_time_loop_index = 1
        end_time_loop_index = 1
        subgrid_dim(t_dim_index) = end_time_nc_index - start_time_nc_index + 1
        dim_length_nc(time_dim_nc_index) = subgrid_dim(t_dim_index)
    end if

    ! Start internal grid receptor loop using only those receptor grids specified in uEMEP_read_receptor_data
    do g_loop = start_grid_loop_index, end_grid_loop_index
        if (use_receptor(g_loop)) then
            ! Set the grid definitions according to the receptor/observation positions
            call uEMEP_set_loop_receptor_grid()

            ! Create the subgrid
            call uEMEP_set_subgrids()

            ! Set emission factors for the current subgrid
            call uEMEP_set_emission_factors()

            ! Start the internal time loop
            do t_loop = start_time_loop_index, end_time_loop_index
                ! Write progress in time and receptor grid loop to screen
                write(*,*) 'REC LOOP= ', g_loop, ' OF ', end_grid_loop_index
                if (unit_logfile .ne. 0) then
                    write(unit_logfile,*) 'REC LOOP= ', g_loop,' OF ', end_grid_loop_index
                end if
                write(*,*) 'TIME LOOP=', t_loop,' OF ', end_time_loop_index
                if (unit_logfile .ne. 0) then
                    write(unit_logfile,*) 'TIME LOOP=', t_loop,' OF ', end_time_loop_index
                end if

                ! For the first time loop set the initial subgrid range values used in reading EMEP and meteo data
                if (t_loop .ge. start_time_loop_index) then
                    init_subgrid_min = subgrid_min
                    init_subgrid_max = subgrid_max
                end if

                ! Read EMEP data from netcdf files. Time stamps based on this
                if ( .not. have_read_emep) then
                    call uEMEP_read_EMEP()
                end if

                ! If read EMEP only once flag is on then turn off the EMEP reading
                ! This is intended for use with multiple receptor files and requires alot of memory so is permanently turned off
                if (read_EMEP_only_once_flag) have_read_emep = .true.

                ! Read meteo grid from netcdf files if required
                if (use_alternative_meteorology_flag .or. use_alternative_z0_flag) then
                    call uEMEP_read_meteo_nc()
                end if

                ! Set the following for the first internal time step only
                if (t_loop .eq. start_time_loop_index) then
                    ! Define subgrid positions and buffer zones. Must be done after reading EMEP data as is based on EMEP grid sizes
                    call uEMEP_define_subgrid_extent()
                    call uEMEP_define_subgrid()

                    ! Define and allocate cross reference subgrids used to transfer data between different subgrids
                    call uEMEP_crossreference_grids()

                    ! Read all road link data from ascii files
                    if (calculate_source(traffic_index) .and. .not. read_subgrid_emission_data) then
                        ! Do this only for the first receptor grid loop
                        if (first_g_loop) then
                            call uEMEP_read_roadlink_data_ascii()
                            call uEMEP_change_road_data()
                            
                            ! Read in the NORTRIP emission data for traffic in the first g_loop if required
                            if (use_NORTRIP_emission_data) then
                                call uEMEP_read_roadlink_emission_data()
                            end if
                        end if
                    end if

                    ! Read in and grid industry data
                    if (calculate_source(industry_index) .and. .not. read_subgrid_emission_data) then
                        call uEMEP_read_industry_data()
                    end if

                    ! Read and subgrid shipping data
                    if (calculate_source(shipping_index) .and. .not. read_subgrid_emission_data) then
                        ! If necessary aggregate shipping data first
                        call uEMEP_preaggregate_shipping_asi_data()
                        
                        ! Read in shipping data
                        if (read_shipping_from_netcdf_flag) then
                            call uEMEP_read_netcdf_shipping_latlon()
                        else
                            if (read_weekly_shipping_data_flag) then
                                call uEMEP_read_weekly_shipping_asi_data()
                            else if (read_monthly_and_daily_shipping_data_flag) then
                                call uEMEP_read_monthly_and_daily_shipping_asi_data()
                            else
                                call uEMEP_read_shipping_asi_data()
                            end if
                        end if
                    end if

                    ! Read in proxy data for home heating. Currently dwelling density
                    if (calculate_source(heating_index) .and. .not. read_subgrid_emission_data) then
                        ! If calculating tiles then read only the dwelling data
                        if (calculate_tiling_flag .or. calculate_region_tiling_flag) then
                            use_RWC_emission_data = .false.
                        end if
                        
                        ! Read the Residential Wood Combustion data from MetVed
                        if (use_RWC_emission_data) then
                            call uEMEP_read_RWC_heating_data()
                        else
                            ! Read and subgrid SSB dwelling data
                            SSB_data_type = dwelling_index
                            if (read_population_from_netcdf_flag) then
                                call uEMEP_read_netcdf_population_latlon()
                            else if (read_population_from_netcdf_local_flag) then
                                call uEMEP_read_netcdf_population()
                            else
                                call uEMEP_read_SSB_data()
                            end if
                        end if
                    end if

                    ! Read and subgrid agriculture data
                    if (calculate_source(agriculture_index) .and. use_rivm_agricuture_emission_data .and. .not. read_subgrid_emission_data) then
                        ! Currently only data from RIVM here
                        call uEMEP_read_agriculture_rivm_data()
                    end if
                    if (read_rivm_landuse_flag) then
                        call uEMEP_read_landuse_rivm_data()
                    end if
                    if (read_subgrid_emission_data) then
                        ! Special routine for reading in RIVM point source emission data
                        if (use_rivm_subgrid_emission_format) then
                            call uEMEP_read_emission_rivm_data()
                        else
                            ! Nothing else available yet
                        end if
                    end if

                    ! Read in population data
                    if (calculate_population_exposure_flag .or. use_population_positions_for_auto_subgrid_flag .or. save_population) then
                        ! Read and subgrid SSB population data
                        SSB_data_type = population_data_type
                        if (read_population_from_netcdf_flag) then
                            call uEMEP_read_netcdf_population_latlon()
                        else if (read_population_from_netcdf_local_flag) then
                            call uEMEP_read_netcdf_population()
                        else
                            call uEMEP_read_SSB_data()
                        end if
                    end if

                    if (use_landuse_as_proxy .or. read_landuse_flag) then
                        call uEMEP_read_netcdf_landuse_latlon()
                    end if

                    ! Autogrid setting for selecting which subgrids to calculate
                    if (use_emission_positions_for_auto_subgrid_flag(allsource_index)) then
                        call uEMEP_grid_roads()
                        write(unit_logfile,*) "'uEMEP_auto_subgrid' has been disabled, because array 'use_subgrid_val' is disabled"
                        stop
                        !call uEMEP_auto_subgrid()
                    end if

                    ! No longer call uEMEP_region_mask, as use_subgrid_val is deactivated and use_subgrid is set elsewhere
                    !if (use_region_select_and_mask_flag) then
                    !    call uEMEP_region_mask()
                    !end if

                    ! New subroutine for reading region mask and region fraction
                    if (trace_emissions_from_in_region .or. use_region_select_and_mask_flag) then
                        call uEMEP_region_mask_new()
                    endif

                    ! Specify the subgrids sizes to be calculated using use_receptor_region
                    call uEMEP_grid_receptor_data

                    ! Carry out tiling. Programme will stop here
                    if (calculate_tiling_flag) then
                        call uEMEP_grid_roads()
                        call uEMEP_set_tile_grids()
                    end if

                    ! Carry out regional tiling. Programme will stop here
                    if (calculate_region_tiling_flag) then
                        call uEMEP_set_region_tile_grids()
                    end if
                end if

                ! Read time profiles for emissions
                call uEMEP_read_time_profiles()

                ! Call grid_roads again to include the time variation from NORTRIP
                if ( .not. read_subgrid_emission_data) then
                    call uEMEP_grid_roads()
                end if

                ! Interpolate meteo data to subgrid. Placed on the integral subgrid
                call uEMEP_subgrid_meteo_EMEP()

                ! Replaces proxy emissions with distributed EMEP emissions
                call uEMEP_subgrid_emission_EMEP()

                ! Convert proxies to emissions including time profiles
                call uEMEP_convert_proxy_to_emissions()

                ! Adjust traffic emissions of NOx based on temperature
                if (use_traffic_nox_emission_temperature_dependency) then
                    call uEMEP_nox_emission_temperature()
                end if

                ! Places EMEP deposition velocities into the deposition_subgrid
                if (calculate_deposition_flag) then
                    call uEMEP_set_deposition_velocities()
                end if

                ! Set travel_time values to 0 outside of the source loop as these are aggregated over all sources
                traveltime_subgrid = 0.0
                
                ! Subgrid dispersion calculation
                do source_index = 1, n_source_index
                    if (calculate_source(source_index) .and. .not. use_plume_dispersion_deposition_flag) then
                        call uEMEP_subgrid_dispersion(source_index)
                    end if
                end do

                do source_index = 1, n_source_index
                    if (calculate_source(source_index) .and. use_plume_dispersion_deposition_flag) then
                        call uEMEP_subgrid_deposition(source_index)
                    end if
                end do

                ! Interpolate local_subgrid if necessary
                if (interpolate_subgrids_flag) then
                    write(unit_logfile,*) "'uEMEP_interpolate_auto_subgrid' has been disabled, because array 'use_subgrid_val' is disabled"
                    stop
                    !call uEMEP_interpolate_auto_subgrid()
                end if

                ! Old diagnostic for comparing EMEP and proxy data emissions. Working only on lat lon EMEP grids. Do not use
                if (make_EMEP_grid_emission_data(allsource_index)) then
                    !call uEMEP_aggregate_proxy_emission_in_EMEP_grid
                end if

                ! Put EMEP data into the additional subgrids for all sources.
                ! Must be run first
                if (EMEP_additional_grid_interpolation_size .gt. 0.0) then
                    calculate_EMEP_additional_grid_flag = .true.
                    call uEMEP_subgrid_EMEP()
                    calculate_EMEP_additional_grid_flag = .false.
                end if

                ! Put EMEP data into subgrids for all sources
                call uEMEP_subgrid_EMEP()

                ! Call the new subroutine for calculating more precise estimates of the contributions from outside moving window but within region
                if (trace_emissions_from_in_region) then
                    if (EMEP_grid_interpolation_flag == 0 .or. EMEP_grid_interpolation_flag == 6) then
                        call uEMEP_subgrid_EMEP_from_in_region()
                    end if
                    ! NB: Only implemented to be consistent with interpolation flag 0 and 6
                end if

                if (calculate_deposition_flag) then
                    call uEMEP_subgrid_deposition_EMEP()
                end if

                ! Interpolate EMEP to sub-grid
                do source_index = 1, n_source_index
                    if (calculate_source(source_index)) then
                        ! Redistributes proxy subgrid data into the EMEP grid concentrations only when local_subgrid_method_flag=1 (based on EMEP concentration redistribution scaling factor)
                        call uEMEP_redistribute_local_source(source_index)
                        
                        ! Places the proxy_subgrid data into the local_subgrid when local_subgrid_method_flag<>1
                        call uEMEP_disperse_local_source(source_index)
                    end if
                end do

                ! Combine and save sources in local and total values
                call uEMEP_combine_local_source()

                ! Calculate the nonlocal depositions
                if (calculate_deposition_flag) then
                    call uEMEP_calculate_deposition()
                end if

                ! Calculate chemistry for NO2 and O3
                call uEMEP_chemistry()

                ! Correct annual mean chemistry for pdf
                if (use_annual_mean_pdf_chemistry_correction) then
                    call correct_annual_mean_chemistry()
                end if

                ! Calculate exposure
                if (calculate_population_exposure_flag) then
                    call uEMEP_calculate_exposure()
                end if

                ! Save results to netcdf
                if (save_netcdf_file_flag .or. save_netcdf_receptor_flag) then
                    call uEMEP_save_netcdf_control()
                end if

            end do ! t_loop

            ! Update first_g_loop flag
            if (first_g_loop) first_g_loop = .false.

        end if ! use_receptor
    end do ! g_loop

    call cpu_time(end_time_cpu)

    if (unit_logfile .ne. 0) then
        write(unit_logfile,*) ''
        write(unit_logfile,*) '------------------------------------------------------------------------'
        write(unit_logfile,*) 'Ending program '//trim(model_version_str)
        write(unit_logfile,'(a,i5,a,i2)') ' CPU time taken (MM:SS): ', floor((end_time_cpu - start_time_cpu)/60.0), ':', floor(mod(end_time_cpu - start_time_cpu, 60.0))
        write(unit_logfile,*) '------------------------------------------------------------------------'
    end if

    if (unit_logfile .gt. 0) then
        close(unit_logfile, status='keep')
    end if

    ! Save finished file
    if (trim(finished_filename) .ne. '') then
        if (save_netcdf_receptor_flag .and. n_valid_receptor .ne. 0) then
            write(*,'(2A)') 'Writing finished file for uEMEP output: ', trim(finished_file_rec)
            open(unit_finishedfile, file=finished_file_rec, status='replace')
            close(unit_finishedfile)
        end if
        if (save_netcdf_file_flag) then
            write(*,'(2A)') 'Writing finished file for uEMEP output: ', trim(finished_file)
            open(unit_finishedfile, file=finished_file, status='replace')
            close(unit_finishedfile)
        end if
    end if

    write(*,*) ''
    write(*,*) '------------------------------------------------------------------------'
    write(*,*) 'Ending program '//trim(model_version_str)
    write(*,'(a,i5,a,i2)') ' CPU time taken (MM:SS): ', floor((end_time_cpu - start_time_cpu)/60.0),':', floor(mod(end_time_cpu - start_time_cpu, 60.0))
    write(*,*) '------------------------------------------------------------------------'

end program uEMEP

