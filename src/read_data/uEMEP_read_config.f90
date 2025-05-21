module read_config

    use uemep_configuration
    use uEMEP_definitions
    use read_namefile_routines, only: read_name_char, read_name_real, &
        read_name_logical, read_name_integer, read_name_double
    use time_functions, only: date_to_datestr_bracket, date_to_number, &
        number_to_date, date_to_datestr_squarebracket, date_to_datestr, &
        datestr_to_date
    use io_functions, only: check_dir_exist

    implicit none
    private

    public :: uEMEP_read_config, replace_string_char

contains

!   Set up of uEMEP before starting calculations
    subroutine uEMEP_read_config

        implicit none

        integer i
        character(256) pathfilename_log_file
        logical :: exists
        integer a(6)
        character(256) format_temp

        integer :: unit_in=30
        integer i_config,i_source,i_landuse
        double precision datenum_temp
        character(256) yesterday_date_str
        character(256) a_str,b_str
        character(256) temp_str

!==========================================================================
!   uEMEP model setup
!==========================================================================

        write(*,'(A)') ''
        write(*,'(A)') '================================================================'
        write(*,'(A)') 'Reading model run configuration (uEMEP_read_config)'
        write(*,'(A)') '================================================================'


        do i_config=1,n_config_files

            !Temporary hardcoding of config file name. Will be read in as command line string
            if (len(name_config_file(i_config)).eq.0) then
                !name_config_file='C:\uEMEP\Fortran\application\config\uEMEP_config_test.txt'
                write (*,'(a)') 'ERROR: No configuration file available. Stopping.'
                stop
            endif


            write (*,'(a)') 'Reading configuration file: '//trim(name_config_file(i_config))

            !Check existence of file
            inquire(file=trim(name_config_file(i_config)),exist=exists)
            if (.not.exists) then
                write(*,'(A)')'ERROR: Configuration file '//trim(name_config_file(i_config))//' does not exist.'
                stop
            endif

            !Open the config file for reading
            open(unit_in,file=name_config_file(i_config),access='sequential',status='old',readonly)

            !First read log file name and open it
            filename_log_file=read_name_char('filename_log_file',filename_log_file,unit_in,unit_logfile)
            pathname_log_file=read_name_char('pathname_log_file',pathname_log_file,unit_in,unit_logfile)

            !Open log file when reading the first configuration file
            if (i_config.eq.1) then
                if (len(trim(filename_log_file)).gt.0) then
                    unit_logfile=10
                    !Check existence of path
                    exists = check_dir_exist(path=trim(pathname_log_file))
                    if (.not.exists) then
                        write(unit_logfile,'(A)')'ERROR: Log file directory path '//trim(pathname_log_file)//' does not exist.'
                        stop
                    endif
                    !Write to screen if writing to log file
                    pathfilename_log_file=trim(pathname_log_file)//trim(filename_log_file)
                    write(*,'(A,A)') 'Writing to log file: ', trim(pathfilename_log_file)
                    write(*,'(A)') '================================================================'
                    open(unit_logfile,file=trim(pathfilename_log_file),access='sequential',form='formatted',status='unknown')
                    if (unit_logfile.ne.0) then
                        write(unit_logfile,*) '------------------------------------------------------------------------'
                        write(unit_logfile,*) 'Starting programm ',trim(model_version_str)
                        write(unit_logfile,*) '------------------------------------------------------------------------'
                    endif
                else
                    unit_logfile=0
                endif

            endif

            file_tag=read_name_char('file_tag',file_tag,unit_in,unit_logfile)

            replacement_date_str=read_name_char('replacement_date_str',replacement_date_str,unit_in,unit_logfile)
            replacement_yesterday_date_str=read_name_char('replacement_yesterday_date_str',replacement_yesterday_date_str,unit_in,unit_logfile)
            replacement_hour_str=read_name_char('replacement_hour_str',replacement_hour_str,unit_in,unit_logfile)
            NORTRIP_replacement_hour_str=read_name_char('NORTRIP_replacement_hour_str',NORTRIP_replacement_hour_str,unit_in,unit_logfile)

            input_comp_name=read_name_char('input_comp_name',input_comp_name,unit_in,unit_logfile)

            hourly_calculations=read_name_logical('hourly_calculations',hourly_calculations,unit_in,unit_logfile)
            annual_calculations=read_name_logical('annual_calculations',annual_calculations,unit_in,unit_logfile)
            !Not implemented
            !start_month_in_annual_calculations=read_name_integer('start_month_in_annual_calculations',start_month_in_annual_calculations,unit_in,unit_logfile)
            !end_month_in_annual_calculations=read_name_integer('end_month_in_annual_calculations',end_month_in_annual_calculations,unit_in,unit_logfile)

            start_time_nc_index=read_name_integer('start_time_nc_index',start_time_nc_index,unit_in,unit_logfile)
            end_time_nc_index=read_name_integer('end_time_nc_index',end_time_nc_index,unit_in,unit_logfile)
            start_time_meteo_nc_index=read_name_integer('start_time_meteo_nc_index',start_time_meteo_nc_index,unit_in,unit_logfile)
            end_time_meteo_nc_index=read_name_integer('end_time_meteo_nc_index',end_time_meteo_nc_index,unit_in,unit_logfile)

            use_single_time_loop_flag=read_name_logical('use_single_time_loop_flag',use_single_time_loop_flag,unit_in,unit_logfile)
            reduce_EMEP_region_flag=read_name_logical('reduce_EMEP_region_flag',reduce_EMEP_region_flag,unit_in,unit_logfile)
            use_multiple_receptor_grids_flag=read_name_logical('use_multiple_receptor_grids_flag',use_multiple_receptor_grids_flag,unit_in,unit_logfile)
            use_receptor_region=read_name_integer('use_receptor_region',use_receptor_region,unit_in,unit_logfile)
            reduce_roadlink_region_flag=read_name_logical('reduce_roadlink_region_flag',reduce_roadlink_region_flag,unit_in,unit_logfile)

            !Choose which sources to calculate
            !calculate_source(:)=read_name_logical('calculate_source(:)',calculate_source(allsource_index),unit_in,unit_logfile)
            calculate_source(traffic_index)=read_name_logical('calculate_source(traffic_index)',calculate_source(traffic_index),unit_in,unit_logfile)
            calculate_source(shipping_index)=read_name_logical('calculate_source(shipping_index)',calculate_source(shipping_index),unit_in,unit_logfile)
            calculate_source(heating_index)=read_name_logical('calculate_source(heating_index)',calculate_source(heating_index),unit_in,unit_logfile)
            calculate_source(agriculture_index)=read_name_logical('calculate_source(agriculture_index)',calculate_source(agriculture_index),unit_in,unit_logfile)
            calculate_source(industry_index)=read_name_logical('calculate_source(industry_index)',calculate_source(industry_index),unit_in,unit_logfile)
            !Additional GNFR sources
            calculate_source(publicpower_index)=read_name_logical('calculate_source(publicpower_index)',calculate_source(publicpower_index),unit_in,unit_logfile)
            calculate_source(fugitive_index)=read_name_logical('calculate_source(fugitive_index)',calculate_source(fugitive_index),unit_in,unit_logfile)
            calculate_source(solvents_index)=read_name_logical('calculate_source(solvents_index)',calculate_source(solvents_index),unit_in,unit_logfile)
            calculate_source(aviation_index)=read_name_logical('calculate_source(aviation_index)',calculate_source(aviation_index),unit_in,unit_logfile)
            calculate_source(offroad_index)=read_name_logical('calculate_source(offroad_index)',calculate_source(offroad_index),unit_in,unit_logfile)
            calculate_source(waste_index)=read_name_logical('calculate_source(waste_index)',calculate_source(waste_index),unit_in,unit_logfile)
            calculate_source(livestock_index)=read_name_logical('calculate_source(livestock_index)',calculate_source(livestock_index),unit_in,unit_logfile)
            calculate_source(other_index)=read_name_logical('calculate_source(other_index)',calculate_source(other_index),unit_in,unit_logfile)

            ! Pollen "sources"
            calculate_source(birch_source_index) = read_name_logical("calculate_source(birch_source_index)", calculate_source(birch_source_index), unit_in, unit_logfile)

            !Choose which EMEP sources to include/calculate. Will not be downscaled but will included as gridded source contributions
            !calculate_source(:)=read_name_logical('calculate_source(:)',calculate_source(allsource_index),unit_in,unit_logfile)
            calculate_EMEP_source(traffic_index)=read_name_logical('calculate_EMEP_source(traffic_index)',calculate_EMEP_source(traffic_index),unit_in,unit_logfile)
            calculate_EMEP_source(shipping_index)=read_name_logical('calculate_EMEP_source(shipping_index)',calculate_EMEP_source(shipping_index),unit_in,unit_logfile)
            calculate_EMEP_source(heating_index)=read_name_logical('calculate_EMEP_source(heating_index)',calculate_EMEP_source(heating_index),unit_in,unit_logfile)
            calculate_EMEP_source(agriculture_index)=read_name_logical('calculate_EMEP_source(agriculture_index)',calculate_EMEP_source(agriculture_index),unit_in,unit_logfile)
            calculate_EMEP_source(industry_index)=read_name_logical('calculate_EMEP_source(industry_index)',calculate_EMEP_source(industry_index),unit_in,unit_logfile)
            !Additional GNFR13 sources
            calculate_EMEP_source(publicpower_index)=read_name_logical('calculate_EMEP_source(publicpower_index)',calculate_EMEP_source(publicpower_index),unit_in,unit_logfile)
            calculate_EMEP_source(fugitive_index)=read_name_logical('calculate_EMEP_source(fugitive_index)',calculate_EMEP_source(fugitive_index),unit_in,unit_logfile)
            calculate_EMEP_source(solvents_index)=read_name_logical('calculate_EMEP_source(solvents_index)',calculate_EMEP_source(solvents_index),unit_in,unit_logfile)
            calculate_EMEP_source(aviation_index)=read_name_logical('calculate_EMEP_source(aviation_index)',calculate_EMEP_source(aviation_index),unit_in,unit_logfile)
            calculate_EMEP_source(offroad_index)=read_name_logical('calculate_EMEP_source(offroad_index)',calculate_EMEP_source(offroad_index),unit_in,unit_logfile)
            calculate_EMEP_source(waste_index)=read_name_logical('calculate_EMEP_source(waste_index)',calculate_EMEP_source(waste_index),unit_in,unit_logfile)
            calculate_EMEP_source(livestock_index)=read_name_logical('calculate_EMEP_source(livestock_index)',calculate_EMEP_source(livestock_index),unit_in,unit_logfile)
            calculate_EMEP_source(other_index)=read_name_logical('calculate_EMEP_source(other_index)',calculate_EMEP_source(other_index),unit_in,unit_logfile)
            !GNFR19 sources. Note nc_index not in the input
            calculate_EMEP_source(traffic_gasoline_nc_index)=read_name_logical('calculate_EMEP_source(traffic_gasoline_index)',calculate_EMEP_source(traffic_gasoline_nc_index),unit_in,unit_logfile)
            calculate_EMEP_source(traffic_diesel_nc_index)=read_name_logical('calculate_EMEP_source(traffic_diesel_index)',calculate_EMEP_source(traffic_diesel_nc_index),unit_in,unit_logfile)
            calculate_EMEP_source(traffic_gas_nc_index)=read_name_logical('calculate_EMEP_source(traffic_gas_index)',calculate_EMEP_source(traffic_gas_nc_index),unit_in,unit_logfile)
            calculate_EMEP_source(traffic_nonexhaust_nc_index)=read_name_logical('calculate_EMEP_source(traffic_nonexhaust_index)',calculate_EMEP_source(traffic_nonexhaust_nc_index),unit_in,unit_logfile)
            calculate_EMEP_source(publicpower_point_nc_index)=read_name_logical('calculate_EMEP_source(publicpower_point_index)',calculate_EMEP_source(publicpower_point_nc_index),unit_in,unit_logfile)
            calculate_EMEP_source(publicpower_area_nc_index)=read_name_logical('calculate_EMEP_source(publicpower_area_index)',calculate_EMEP_source(publicpower_area_nc_index),unit_in,unit_logfile)


            !do i_source=1,n_source_index
            !    if (calculate_source(i_source)) calculate_source(allsource_index)=.true.
            !enddo

            !For aggregating proxy emission data to EMEP grids in ascii test routines that no longer exist. Default is false
            make_EMEP_grid_emission_data(traffic_index)=read_name_logical('make_EMEP_grid_emission_data(traffic_index)',make_EMEP_grid_emission_data(traffic_index),unit_in,unit_logfile)
            make_EMEP_grid_emission_data(shipping_index)=read_name_logical('make_EMEP_grid_emission_data(shipping_index)',make_EMEP_grid_emission_data(shipping_index),unit_in,unit_logfile)
            make_EMEP_grid_emission_data(heating_index)=read_name_logical('make_EMEP_grid_emission_data(heating_index)',make_EMEP_grid_emission_data(heating_index),unit_in,unit_logfile)
            make_EMEP_grid_emission_data(agriculture_index)=read_name_logical('make_EMEP_grid_emission_data(agriculture_index)',make_EMEP_grid_emission_data(agriculture_index),unit_in,unit_logfile)
            make_EMEP_grid_emission_data(industry_index)=read_name_logical('make_EMEP_grid_emission_data(industry_index)',make_EMEP_grid_emission_data(industry_index),unit_in,unit_logfile)
            !Additional GNFR sources
            make_EMEP_grid_emission_data(publicpower_index)=read_name_logical('make_EMEP_grid_emission_data(publicpower_index)',make_EMEP_grid_emission_data(publicpower_index),unit_in,unit_logfile)
            make_EMEP_grid_emission_data(fugitive_index)=read_name_logical('make_EMEP_grid_emission_data(fugitive_index)',make_EMEP_grid_emission_data(fugitive_index),unit_in,unit_logfile)
            make_EMEP_grid_emission_data(solvents_index)=read_name_logical('make_EMEP_grid_emission_data(solvents_index)',make_EMEP_grid_emission_data(solvents_index),unit_in,unit_logfile)
            make_EMEP_grid_emission_data(aviation_index)=read_name_logical('make_EMEP_grid_emission_data(aviation_index)',make_EMEP_grid_emission_data(aviation_index),unit_in,unit_logfile)
            make_EMEP_grid_emission_data(offroad_index)=read_name_logical('make_EMEP_grid_emission_data(offroad_index)',make_EMEP_grid_emission_data(offroad_index),unit_in,unit_logfile)
            make_EMEP_grid_emission_data(waste_index)=read_name_logical('make_EMEP_grid_emission_data(waste_index)',make_EMEP_grid_emission_data(waste_index),unit_in,unit_logfile)
            make_EMEP_grid_emission_data(livestock_index)=read_name_logical('make_EMEP_grid_emission_data(livestock_index)',make_EMEP_grid_emission_data(livestock_index),unit_in,unit_logfile)
            make_EMEP_grid_emission_data(other_index)=read_name_logical('make_EMEP_grid_emission_data(other_index)',make_EMEP_grid_emission_data(other_index),unit_in,unit_logfile)
            do i_source=1,n_source_index
                if (make_EMEP_grid_emission_data(i_source)) make_EMEP_grid_emission_data(allsource_index)=.true.
            enddo

            !For scaling EMEP concentrations with EMEP/local emission ratio. for testing purposes only. Default is false
            replace_EMEP_local_with_subgrid_local(traffic_index)=read_name_logical('replace_EMEP_local_with_subgrid_local(traffic_index)',replace_EMEP_local_with_subgrid_local(traffic_index),unit_in,unit_logfile)
            replace_EMEP_local_with_subgrid_local(shipping_index)=read_name_logical('replace_EMEP_local_with_subgrid_local(shipping_index)',replace_EMEP_local_with_subgrid_local(shipping_index),unit_in,unit_logfile)
            replace_EMEP_local_with_subgrid_local(heating_index)=read_name_logical('replace_EMEP_local_with_subgrid_local(heating_index)',replace_EMEP_local_with_subgrid_local(heating_index),unit_in,unit_logfile)
            replace_EMEP_local_with_subgrid_local(agriculture_index)=read_name_logical('replace_EMEP_local_with_subgrid_local(agriculture_index)',replace_EMEP_local_with_subgrid_local(agriculture_index),unit_in,unit_logfile)
            replace_EMEP_local_with_subgrid_local(industry_index)=read_name_logical('replace_EMEP_local_with_subgrid_local(industry_index)',replace_EMEP_local_with_subgrid_local(industry_index),unit_in,unit_logfile)

            projection_type=read_name_integer('projection_type',projection_type,unit_in,unit_logfile)
            EMEP_projection_type=read_name_integer('EMEP_projection_type',EMEP_projection_type,unit_in,unit_logfile)
            utm_zone=read_name_integer('utm_zone',utm_zone,unit_in,unit_logfile)
            !Present UTM central lon position if not overridden by input
            utm_lon0=abs(utm_zone)*6-180-3
            utm_lon0=read_name_real('utm_lon0',utm_lon0,unit_in,unit_logfile)
            ltm_lon0=read_name_real('ltm_lon0',ltm_lon0,unit_in,unit_logfile)

            !Read the projection attributes for uEMEP if they are available
            projection_attributes(1)=read_name_double('projection_attributes(1)',projection_attributes(1),unit_in,unit_logfile)
            projection_attributes(2)=read_name_double('projection_attributes(2)',projection_attributes(2),unit_in,unit_logfile)
            projection_attributes(3)=read_name_double('projection_attributes(3)',projection_attributes(3),unit_in,unit_logfile)
            projection_attributes(4)=read_name_double('projection_attributes(4)',projection_attributes(4),unit_in,unit_logfile)
            projection_attributes(5)=read_name_double('projection_attributes(5)',projection_attributes(5),unit_in,unit_logfile)
            projection_attributes(6)=read_name_double('projection_attributes(6)',projection_attributes(6),unit_in,unit_logfile)
            projection_attributes(7)=read_name_double('projection_attributes(7)',projection_attributes(7),unit_in,unit_logfile)

            !Read the projection attributes for EMEP if they are available
            !THese will be reset if EMEP data is read in
            EMEP_projection_attributes(1)=read_name_double('EMEP_projection_attributes(1)',EMEP_projection_attributes(1),unit_in,unit_logfile)
            EMEP_projection_attributes(2)=read_name_double('EMEP_projection_attributes(2)',EMEP_projection_attributes(2),unit_in,unit_logfile)
            EMEP_projection_attributes(3)=read_name_double('EMEP_projection_attributes(3)',EMEP_projection_attributes(3),unit_in,unit_logfile)
            EMEP_projection_attributes(4)=read_name_double('EMEP_projection_attributes(4)',EMEP_projection_attributes(4),unit_in,unit_logfile)
            EMEP_projection_attributes(5)=read_name_double('EMEP_projection_attributes(5)',EMEP_projection_attributes(5),unit_in,unit_logfile)
            EMEP_projection_attributes(6)=read_name_double('EMEP_projection_attributes(6)',EMEP_projection_attributes(6),unit_in,unit_logfile)
            EMEP_projection_attributes(7)=read_name_double('EMEP_projection_attributes(7)',EMEP_projection_attributes(7),unit_in,unit_logfile)

            EMEP_grid_interpolation_flag=read_name_integer('EMEP_grid_interpolation_flag',EMEP_grid_interpolation_flag,unit_in,unit_logfile)
            EMEP_meteo_grid_interpolation_flag=read_name_integer('EMEP_meteo_grid_interpolation_flag',EMEP_meteo_grid_interpolation_flag,unit_in,unit_logfile)
            EMEP_emission_grid_interpolation_flag=read_name_integer('EMEP_emission_grid_interpolation_flag',EMEP_emission_grid_interpolation_flag,unit_in,unit_logfile)
            subgrid_emission_distribution_flag=read_name_logical('subgrid_emission_distribution_flag',subgrid_emission_distribution_flag,unit_in,unit_logfile)
            !Not used
            EMEP_grid_interpolation_simple_flag=read_name_logical('EMEP_grid_interpolation_simple_flag',EMEP_grid_interpolation_simple_flag,unit_in,unit_logfile)

            EMEP_grid_interpolation_size=read_name_real('EMEP_grid_interpolation_size',EMEP_grid_interpolation_size,unit_in,unit_logfile)
            EMEP_additional_grid_interpolation_size=read_name_real('EMEP_additional_grid_interpolation_size',EMEP_additional_grid_interpolation_size,unit_in,unit_logfile)
            use_downwind_position_flag=read_name_logical('use_downwind_position_flag',use_downwind_position_flag,unit_in,unit_logfile)
            local_subgrid_method_flag=read_name_integer('local_subgrid_method_flag',local_subgrid_method_flag,unit_in,unit_logfile)

            stability_scheme_flag=read_name_integer('stability_scheme_flag',stability_scheme_flag,unit_in,unit_logfile)
            average_zc_h_in_Kz_flag=read_name_logical('average_zc_h_in_Kz_flag',average_zc_h_in_Kz_flag,unit_in,unit_logfile)

            wind_level_flag=read_name_integer('wind_level_flag',wind_level_flag,unit_in,unit_logfile)
            wind_level_integral_flag=read_name_integer('wind_level_integral_flag',wind_level_flag,unit_in,unit_logfile) !Default is wind_level_flag
            no2_chemistry_scheme_flag=read_name_integer('no2_chemistry_scheme_flag',no2_chemistry_scheme_flag,unit_in,unit_logfile)
            no2_background_chemistry_scheme_flag=read_name_integer('no2_background_chemistry_scheme_flag',no2_background_chemistry_scheme_flag,unit_in,unit_logfile)

            use_emission_positions_for_auto_subgrid_flag(traffic_index)=read_name_logical('use_emission_positions_for_auto_subgrid_flag(traffic_index)',use_emission_positions_for_auto_subgrid_flag(traffic_index),unit_in,unit_logfile)
            use_emission_positions_for_auto_subgrid_flag(shipping_index)=read_name_logical('use_emission_positions_for_auto_subgrid_flag(shipping_index)',use_emission_positions_for_auto_subgrid_flag(shipping_index),unit_in,unit_logfile)
            use_emission_positions_for_auto_subgrid_flag(heating_index)=read_name_logical('use_emission_positions_for_auto_subgrid_flag(heating_index)',use_emission_positions_for_auto_subgrid_flag(heating_index),unit_in,unit_logfile)
            use_emission_positions_for_auto_subgrid_flag(agriculture_index)=read_name_logical('use_emission_positions_for_auto_subgrid_flag(agriculture_index)',use_emission_positions_for_auto_subgrid_flag(agriculture_index),unit_in,unit_logfile)
            use_emission_positions_for_auto_subgrid_flag(industry_index)=read_name_logical('use_emission_positions_for_auto_subgrid_flag(industry_index)',use_emission_positions_for_auto_subgrid_flag(industry_index),unit_in,unit_logfile)
            !Set all source index to true if any of the sources are to be auto gridded. allsource_index defines if the routine is called or not
            do i_source=1,n_source_index
                if (use_emission_positions_for_auto_subgrid_flag(i_source).and.i_source.ne.allsource_index) use_emission_positions_for_auto_subgrid_flag(allsource_index)=.true.
            enddo
            !do i_source=1,n_source_index
            !write(*,*) 'Using auto subgrid for source ',trim(source_file_str(i_source)),use_emission_positions_for_auto_subgrid_flag(i_source)
            !enddo

            use_receptor_positions_for_auto_subgrid_flag=read_name_logical('use_receptor_positions_for_auto_subgrid_flag',use_receptor_positions_for_auto_subgrid_flag,unit_in,unit_logfile)
            use_population_positions_for_auto_subgrid_flag=read_name_logical('use_population_positions_for_auto_subgrid_flag',use_population_positions_for_auto_subgrid_flag,unit_in,unit_logfile)
            interpolate_subgrids_flag=read_name_logical('interpolate_subgrids_flag',interpolate_subgrids_flag,unit_in,unit_logfile)

            !use_trajectory_flag(:)=read_name_logical('use_trajectory_flag',use_trajectory_flag(allsource_index),unit_in,unit_logfile)
            use_trajectory_flag(shipping_index)=read_name_logical('use_trajectory_flag(shipping_index)',use_trajectory_flag(shipping_index),unit_in,unit_logfile)
            use_trajectory_flag(traffic_index)=read_name_logical('use_trajectory_flag(traffic_index)',use_trajectory_flag(traffic_index),unit_in,unit_logfile)
            use_trajectory_flag(heating_index)=read_name_logical('use_trajectory_flag(heating_index)',use_trajectory_flag(heating_index),unit_in,unit_logfile)
            use_trajectory_flag(agriculture_index)=read_name_logical('use_trajectory_flag(agriculture_index)',use_trajectory_flag(agriculture_index),unit_in,unit_logfile)
            use_trajectory_flag(industry_index)=read_name_logical('use_trajectory_flag(industry_index)',use_trajectory_flag(industry_index),unit_in,unit_logfile)
            !Additional GNFR sources
            use_trajectory_flag(publicpower_index)=read_name_logical('use_trajectory_flag(publicpower_index)',use_trajectory_flag(publicpower_index),unit_in,unit_logfile)
            use_trajectory_flag(fugitive_index)=read_name_logical('use_trajectory_flag(fugitive_index)',use_trajectory_flag(fugitive_index),unit_in,unit_logfile)
            use_trajectory_flag(solvents_index)=read_name_logical('use_trajectory_flag(solvents_index)',use_trajectory_flag(solvents_index),unit_in,unit_logfile)
            use_trajectory_flag(aviation_index)=read_name_logical('use_trajectory_flag(aviation_index)',use_trajectory_flag(aviation_index),unit_in,unit_logfile)
            use_trajectory_flag(offroad_index)=read_name_logical('use_trajectory_flag(offroad_index)',use_trajectory_flag(offroad_index),unit_in,unit_logfile)
            use_trajectory_flag(waste_index)=read_name_logical('use_trajectory_flag(waste_index)',use_trajectory_flag(waste_index),unit_in,unit_logfile)
            use_trajectory_flag(livestock_index)=read_name_logical('use_trajectory_flag(livestock_index)',use_trajectory_flag(livestock_index),unit_in,unit_logfile)
            use_trajectory_flag(other_index)=read_name_logical('use_trajectory_flag(other_index)',use_trajectory_flag(other_index),unit_in,unit_logfile)

            traj_step_scale=read_name_real('traj_step_scale',traj_step_scale,unit_in,unit_logfile)

            calculate_aggregated_shipping_emissions_flag=read_name_logical('calculate_aggregated_shipping_emissions_flag',calculate_aggregated_shipping_emissions_flag,unit_in,unit_logfile)
            use_aggregated_shipping_emissions_flag=read_name_logical('use_aggregated_shipping_emissions_flag',use_aggregated_shipping_emissions_flag,unit_in,unit_logfile)

            calculate_population_exposure_flag=read_name_logical('calculate_population_exposure_flag',calculate_population_exposure_flag,unit_in,unit_logfile)


            !Grid specific parameters
            integral_subgrid_step=read_name_integer('integral_subgrid_step',integral_subgrid_step,unit_in,unit_logfile)

            !Specify the site grid information
            subgrid_delta(x_dim_index)=read_name_real('subgrid_delta(x_dim_index)',subgrid_delta(x_dim_index),unit_in,unit_logfile)
            subgrid_delta(y_dim_index)=read_name_real('subgrid_delta(y_dim_index)',subgrid_delta(y_dim_index),unit_in,unit_logfile)
            subgrid_min(x_dim_index)=read_name_real('subgrid_min(x_dim_index)',subgrid_min(x_dim_index),unit_in,unit_logfile)
            subgrid_min(y_dim_index)=read_name_real('subgrid_min(y_dim_index)',subgrid_min(y_dim_index),unit_in,unit_logfile)
            subgrid_max(x_dim_index)=read_name_real('subgrid_max(x_dim_index)',subgrid_max(x_dim_index),unit_in,unit_logfile)
            subgrid_max(y_dim_index)=read_name_real('subgrid_max(y_dim_index)',subgrid_max(y_dim_index),unit_in,unit_logfile)

            !Save the read in grid data. This will be used to select EMEP region and receptor points
            init_subgrid_delta(x_dim_index)=subgrid_delta(x_dim_index)
            init_subgrid_delta(y_dim_index)=subgrid_delta(y_dim_index)
            init_subgrid_min(x_dim_index)=subgrid_min(x_dim_index)
            init_subgrid_min(y_dim_index)=subgrid_min(y_dim_index)
            init_subgrid_max(x_dim_index)=subgrid_max(x_dim_index)
            init_subgrid_max(y_dim_index)=subgrid_max(y_dim_index)

            !These will probably be set by the input data
            deposition_subgrid_delta(x_dim_index)=read_name_real('deposition_subgrid_delta(x_dim_index)',deposition_subgrid_delta(x_dim_index),unit_in,unit_logfile)
            deposition_subgrid_delta(y_dim_index)=read_name_real('deposition_subgrid_delta(y_dim_index)',deposition_subgrid_delta(y_dim_index),unit_in,unit_logfile)
            landuse_subgrid_delta(x_dim_index)=read_name_real('landuse_subgrid_delta(x_dim_index)',landuse_subgrid_delta(x_dim_index),unit_in,unit_logfile)
            landuse_subgrid_delta(y_dim_index)=read_name_real('landuse_subgrid_delta(y_dim_index)',landuse_subgrid_delta(y_dim_index),unit_in,unit_logfile)

            !Specifies the number of subsources for each source. Not used as default is 1
            n_subsource(traffic_index)=read_name_integer('n_subsource(traffic_index)',n_subsource(traffic_index),unit_in,unit_logfile)
            n_subsource(shipping_index)=read_name_integer('n_subsource(shipping_index)',n_subsource(shipping_index),unit_in,unit_logfile)
            n_subsource(heating_index)=read_name_integer('n_subsource(heating_index)',n_subsource(heating_index),unit_in,unit_logfile)
            n_subsource(agriculture_index)=read_name_integer('n_subsource(agriculture_index)',n_subsource(agriculture_index),unit_in,unit_logfile)
            n_subsource(industry_index)=read_name_integer('n_subsource(industry_index)',n_subsource(industry_index),unit_in,unit_logfile)

            h_emis(traffic_index,1)=read_name_real('h_emis(traffic_index,1)',h_emis(traffic_index,1),unit_in,unit_logfile)
            h_emis(shipping_index,1)=read_name_real('h_emis(shipping_index,1)',h_emis(shipping_index,1),unit_in,unit_logfile)
            h_emis(heating_index,1)=read_name_real('h_emis(heating_index,1)',h_emis(heating_index,1),unit_in,unit_logfile)
            h_emis(agriculture_index,1)=read_name_real('h_emis(agriculture_index,1)',h_emis(agriculture_index,1),unit_in,unit_logfile)
            h_emis(industry_index,1)=read_name_real('h_emis(industry_index,1)',h_emis(industry_index,1),unit_in,unit_logfile)
            !Additional GNFR sources
            h_emis(publicpower_index,1)=read_name_real('h_emis(publicpower_index,1)',h_emis(publicpower_index,1),unit_in,unit_logfile)
            h_emis(fugitive_index,1)=read_name_real('h_emis(fugitive_index,1)',h_emis(fugitive_index,1),unit_in,unit_logfile)
            h_emis(solvents_index,1)=read_name_real('h_emis(solvents_index,1)',h_emis(solvents_index,1),unit_in,unit_logfile)
            h_emis(aviation_index,1)=read_name_real('h_emis(aviation_index,1)',h_emis(aviation_index,1),unit_in,unit_logfile)
            h_emis(offroad_index,1)=read_name_real('h_emis(offroad_index,1)',h_emis(offroad_index,1),unit_in,unit_logfile)
            h_emis(waste_index,1)=read_name_real('h_emis(waste_index,1)',h_emis(waste_index,1),unit_in,unit_logfile)
            h_emis(livestock_index,1)=read_name_real('h_emis(livestock_index,1)',h_emis(livestock_index,1),unit_in,unit_logfile)
            h_emis(other_index,1)=read_name_real('h_emis(other_index,1)',h_emis(other_index,1),unit_in,unit_logfile)

            !These second subsources do not exist but still possible to implement
            h_emis(traffic_index,2)=read_name_real('h_emis(traffic_index,2)',h_emis(traffic_index,2),unit_in,unit_logfile)
            h_emis(shipping_index,2)=read_name_real('h_emis(shipping_index,2)',h_emis(shipping_index,2),unit_in,unit_logfile)
            h_emis(heating_index,2)=read_name_real('h_emis(heating_index,2)',h_emis(heating_index,2),unit_in,unit_logfile)
            h_emis(agriculture_index,2)=read_name_real('h_emis(agriculture_index,2)',h_emis(agriculture_index,2),unit_in,unit_logfile)
            h_emis(industry_index,2)=read_name_real('h_emis(industry_index,2)',h_emis(industry_index,2),unit_in,unit_logfile)

            sig_y_00(traffic_index,1)=read_name_real('sig_y_00(traffic_index,1)',sig_y_00(traffic_index,1),unit_in,unit_logfile)
            sig_y_00(shipping_index,1)=read_name_real('sig_y_00(shipping_index,1)',sig_y_00(shipping_index,1),unit_in,unit_logfile)
            sig_y_00(heating_index,1)=read_name_real('sig_y_00(heating_index,1)',sig_y_00(heating_index,1),unit_in,unit_logfile)
            sig_y_00(agriculture_index,1)=read_name_real('sig_y_00(agriculture_index,1)',sig_y_00(agriculture_index,1),unit_in,unit_logfile)
            sig_y_00(industry_index,1)=read_name_real('sig_y_00(industry_index,1)',sig_y_00(industry_index,1),unit_in,unit_logfile)
            !Additional GNFR sources
            sig_y_00(publicpower_index,1)=read_name_real('sig_y_00(publicpower_index,1)',sig_y_00(publicpower_index,1),unit_in,unit_logfile)
            sig_y_00(fugitive_index,1)=read_name_real('sig_y_00(fugitive_index,1)',sig_y_00(fugitive_index,1),unit_in,unit_logfile)
            sig_y_00(solvents_index,1)=read_name_real('sig_y_00(solvents_index,1)',sig_y_00(solvents_index,1),unit_in,unit_logfile)
            sig_y_00(aviation_index,1)=read_name_real('sig_y_00(aviation_index,1)',sig_y_00(aviation_index,1),unit_in,unit_logfile)
            sig_y_00(offroad_index,1)=read_name_real('sig_y_00(offroad_index,1)',sig_y_00(offroad_index,1),unit_in,unit_logfile)
            sig_y_00(waste_index,1)=read_name_real('sig_y_00(waste_index,1)',sig_y_00(waste_index,1),unit_in,unit_logfile)
            sig_y_00(livestock_index,1)=read_name_real('sig_y_00(livestock_index,1)',sig_y_00(livestock_index,1),unit_in,unit_logfile)
            sig_y_00(other_index,1)=read_name_real('sig_y_00(other_index,1)',sig_y_00(other_index,1),unit_in,unit_logfile)

            !These second subsources do not exist but still possible to implement
            sig_y_00(traffic_index,2)=read_name_real('sig_y_00(traffic_index,2)',sig_y_00(traffic_index,2),unit_in,unit_logfile)
            sig_y_00(shipping_index,2)=read_name_real('sig_y_00(shipping_index,2)',sig_y_00(shipping_index,2),unit_in,unit_logfile)
            sig_y_00(heating_index,2)=read_name_real('sig_y_00(heating_index,2)',sig_y_00(heating_index,2),unit_in,unit_logfile)
            sig_y_00(agriculture_index,2)=read_name_real('sig_y_00(agriculture_index,2)',sig_y_00(agriculture_index,2),unit_in,unit_logfile)
            sig_y_00(industry_index,2)=read_name_real('sig_y_00(industry_index,2)',sig_y_00(industry_index,2),unit_in,unit_logfile)

            !This scales the subgrid size to give to give sig_y_0. Optimal value is 0.8 so sig_y_0=0.8/2*delta_y
            sigy_0_subgid_width_scale=read_name_real('sigy_0_subgid_width_scale',sigy_0_subgid_width_scale,unit_in,unit_logfile)

            sig_z_00(traffic_index,1)=read_name_real('sig_z_00(traffic_index,1)',sig_z_00(traffic_index,1),unit_in,unit_logfile)
            sig_z_00(shipping_index,1)=read_name_real('sig_z_00(shipping_index,1)',sig_z_00(shipping_index,1),unit_in,unit_logfile)
            sig_z_00(heating_index,1)=read_name_real('sig_z_00(heating_index,1)',sig_z_00(heating_index,1),unit_in,unit_logfile)
            sig_z_00(agriculture_index,1)=read_name_real('sig_z_00(agriculture_index,1)',sig_z_00(agriculture_index,1),unit_in,unit_logfile)
            sig_z_00(industry_index,1)=read_name_real('sig_z_00(industry_index,1)',sig_z_00(industry_index,1),unit_in,unit_logfile)
            !Additional GNFR sources
            sig_z_00(publicpower_index,1)=read_name_real('sig_z_00(publicpower_index,1)',sig_z_00(publicpower_index,1),unit_in,unit_logfile)
            sig_z_00(fugitive_index,1)=read_name_real('sig_z_00(fugitive_index,1)',sig_z_00(fugitive_index,1),unit_in,unit_logfile)
            sig_z_00(solvents_index,1)=read_name_real('sig_z_00(solvents_index,1)',sig_z_00(solvents_index,1),unit_in,unit_logfile)
            sig_z_00(aviation_index,1)=read_name_real('sig_z_00(aviation_index,1)',sig_z_00(aviation_index,1),unit_in,unit_logfile)
            sig_z_00(offroad_index,1)=read_name_real('sig_z_00(offroad_index,1)',sig_z_00(offroad_index,1),unit_in,unit_logfile)
            sig_z_00(waste_index,1)=read_name_real('sig_z_00(waste_index,1)',sig_z_00(waste_index,1),unit_in,unit_logfile)
            sig_z_00(livestock_index,1)=read_name_real('sig_z_00(livestock_index,1)',sig_z_00(livestock_index,1),unit_in,unit_logfile)
            sig_z_00(other_index,1)=read_name_real('sig_z_00(other_index,1)',sig_z_00(other_index,1),unit_in,unit_logfile)

            !These second subsources do not exist but still possible to implement
            sig_z_00(traffic_index,2)=read_name_real('sig_z_00(traffic_index,2)',sig_z_00(traffic_index,2),unit_in,unit_logfile)
            sig_z_00(shipping_index,2)=read_name_real('sig_z_00(shipping_index,2)',sig_z_00(shipping_index,2),unit_in,unit_logfile)
            sig_z_00(heating_index,2)=read_name_real('sig_z_00(heating_index,2)',sig_z_00(heating_index,2),unit_in,unit_logfile)
            sig_z_00(agriculture_index,2)=read_name_real('sig_z_00(agriculture_index,2)',sig_z_00(agriculture_index,2),unit_in,unit_logfile)
            sig_z_00(industry_index,2)=read_name_real('sig_z_00(industry_index,2)',sig_z_00(industry_index,2),unit_in,unit_logfile)

            !Read output grid path for all data
            pathname_output_grid=read_name_char('pathname_output_grid',pathname_output_grid,unit_in,unit_logfile)
            filename_date_output_grid=read_name_char('filename_date_output_grid',filename_date_output_grid,unit_in,unit_logfile)

            !Read in file names. Only 2 choices for most file types
            pathname_rl(1)=read_name_char('pathname_rl(1)',pathname_rl(1),unit_in,unit_logfile)
            pathname_rl(2)=read_name_char('pathname_rl(2)',pathname_rl(2),unit_in,unit_logfile)
            filename_rl(1)=read_name_char('filename_rl(1)',filename_rl(1),unit_in,unit_logfile)
            filename_rl(2)=read_name_char('filename_rl(2)',filename_rl(2),unit_in,unit_logfile)

            pathname_mrl(1)=read_name_char('pathname_mrl(1)',pathname_mrl(1),unit_in,unit_logfile)
            pathname_mrl(2)=read_name_char('pathname_mrl(2)',pathname_mrl(2),unit_in,unit_logfile)
            pathname_mrl(3)=read_name_char('pathname_mrl(3)',pathname_mrl(3),unit_in,unit_logfile)
            pathname_mrl(4)=read_name_char('pathname_mrl(4)',pathname_mrl(4),unit_in,unit_logfile)
            pathname_mrl(5)=read_name_char('pathname_mrl(5)',pathname_mrl(5),unit_in,unit_logfile)
            pathname_mrl(6)=read_name_char('pathname_mrl(6)',pathname_mrl(6),unit_in,unit_logfile)
            pathname_mrl(7)=read_name_char('pathname_mrl(7)',pathname_mrl(7),unit_in,unit_logfile)
            pathname_mrl(8)=read_name_char('pathname_mrl(8)',pathname_mrl(8),unit_in,unit_logfile)
            pathname_mrl(9)=read_name_char('pathname_mrl(9)',pathname_mrl(9),unit_in,unit_logfile)
            pathname_mrl(10)=read_name_char('pathname_mrl(10)',pathname_mrl(10),unit_in,unit_logfile)
            filename_mrl(1)=read_name_char('filename_mrl(1)',filename_mrl(1),unit_in,unit_logfile)
            filename_mrl(2)=read_name_char('filename_mrl(2)',filename_mrl(2),unit_in,unit_logfile)
            filename_mrl(3)=read_name_char('filename_mrl(3)',filename_mrl(3),unit_in,unit_logfile)
            filename_mrl(4)=read_name_char('filename_mrl(4)',filename_mrl(4),unit_in,unit_logfile)
            filename_mrl(5)=read_name_char('filename_mrl(5)',filename_mrl(5),unit_in,unit_logfile)
            filename_mrl(6)=read_name_char('filename_mrl(6)',filename_mrl(6),unit_in,unit_logfile)
            filename_mrl(7)=read_name_char('filename_mrl(7)',filename_mrl(7),unit_in,unit_logfile)
            filename_mrl(8)=read_name_char('filename_mrl(8)',filename_mrl(8),unit_in,unit_logfile)
            filename_mrl(9)=read_name_char('filename_mrl(9)',filename_mrl(9),unit_in,unit_logfile)
            filename_mrl(10)=read_name_char('filename_mrl(10)',filename_mrl(10),unit_in,unit_logfile)
            num_multiple_roadlink_files=read_name_integer('num_multiple_roadlink_files',num_multiple_roadlink_files,unit_in,unit_logfile)

            pathname_EMEP(1)=read_name_char('pathname_EMEP(1)',pathname_EMEP(1),unit_in,unit_logfile)
            pathname_EMEP(2)=read_name_char('pathname_EMEP(2)',pathname_EMEP(2),unit_in,unit_logfile)
            pathname_EMEP(3)=read_name_char('pathname_EMEP(3)',pathname_EMEP(3),unit_in,unit_logfile)
            pathname_EMEP(4)=read_name_char('pathname_EMEP(4)',pathname_EMEP(4),unit_in,unit_logfile)
            filename_EMEP(1)=read_name_char('filename_EMEP(1)',filename_EMEP(1),unit_in,unit_logfile)
            filename_EMEP(2)=read_name_char('filename_EMEP(2)',filename_EMEP(2),unit_in,unit_logfile)
            filename_EMEP(3)=read_name_char('filename_EMEP(3)',filename_EMEP(3),unit_in,unit_logfile)
            filename_EMEP(4)=read_name_char('filename_EMEP(4)',filename_EMEP(4),unit_in,unit_logfile)
            original_pathname_EMEP=pathname_EMEP
            original_filename_EMEP=filename_EMEP


            pathname_ship(1)=read_name_char('pathname_ship(1)',pathname_ship(1),unit_in,unit_logfile)
            pathname_ship(2)=read_name_char('pathname_ship(2)',pathname_ship(2),unit_in,unit_logfile)
            filename_ship(1)=read_name_char('filename_ship(1)',filename_ship(1),unit_in,unit_logfile)
            filename_ship(2)=read_name_char('filename_ship(2)',filename_ship(2),unit_in,unit_logfile)

            pathname_agriculture(1)=read_name_char('pathname_agriculture(1)',pathname_agriculture(1),unit_in,unit_logfile)
            pathname_agriculture(2)=read_name_char('pathname_agriculture(2)',pathname_agriculture(2),unit_in,unit_logfile)
            filename_agriculture(1)=read_name_char('filename_agriculture(1)',filename_agriculture(1),unit_in,unit_logfile)
            filename_agriculture(2)=read_name_char('filename_agriculture(2)',filename_agriculture(2),unit_in,unit_logfile)

            pathname_emission_rivm(1)=read_name_char('pathname_emission_rivm(1)',pathname_emission_rivm(1),unit_in,unit_logfile)
            pathname_emission_rivm(2)=read_name_char('pathname_emission_rivm(2)',pathname_emission_rivm(2),unit_in,unit_logfile)
            filename_emission_rivm(1)=read_name_char('filename_emission_rivm(1)',filename_emission_rivm(1),unit_in,unit_logfile)
            filename_emission_rivm(2)=read_name_char('filename_emission_rivm(2)',filename_emission_rivm(2),unit_in,unit_logfile)

            pathname_industry(1)=read_name_char('pathname_industry(1)',pathname_industry(1),unit_in,unit_logfile)
            pathname_industry(2)=read_name_char('pathname_industry(2)',pathname_industry(2),unit_in,unit_logfile)
            filename_industry(1)=read_name_char('filename_industry(1)',filename_industry(1),unit_in,unit_logfile)
            filename_industry(2)=read_name_char('filename_industry(2)',filename_industry(2),unit_in,unit_logfile)

            pathname_heating(dwelling_index)=read_name_char('pathname_heating(dwelling_index)',pathname_heating(dwelling_index),unit_in,unit_logfile)
            pathname_heating(population_index)=read_name_char('pathname_heating(population_index)',pathname_heating(population_index),unit_in,unit_logfile)
            filename_heating(dwelling_index)=read_name_char('filename_heating(dwelling_index)',filename_heating(dwelling_index),unit_in,unit_logfile)
            filename_heating(population_index)=read_name_char('filename_heating(population_index)',filename_heating(population_index),unit_in,unit_logfile)
            pathname_heating(RWC_heating_index)=read_name_char('pathname_heating(RWC_heating_index)',pathname_heating(RWC_heating_index),unit_in,unit_logfile)
            filename_heating(RWC_heating_index)=read_name_char('filename_heating(RWC_heating_index)',filename_heating(RWC_heating_index),unit_in,unit_logfile)

            pathname_population(dwelling_index)=read_name_char('pathname_population(dwelling_index)',pathname_population(dwelling_index),unit_in,unit_logfile)
            pathname_population(population_index)=read_name_char('pathname_population(population_index)',pathname_population(population_index),unit_in,unit_logfile)
            pathname_population(establishment_index)=read_name_char('pathname_population(establishment_index)',pathname_population(establishment_index),unit_in,unit_logfile)
            pathname_population(school_index)=read_name_char('pathname_population(school_index)',pathname_population(school_index),unit_in,unit_logfile)
            pathname_population(kindergaten_index)=read_name_char('pathname_population(kindergaten_index)',pathname_population(kindergaten_index),unit_in,unit_logfile)
            pathname_population(home_index)=read_name_char('pathname_population(home_index)',pathname_population(home_index),unit_in,unit_logfile)
            pathname_population(municipality_index)=read_name_char('pathname_population(municipality_index)',pathname_population(municipality_index),unit_in,unit_logfile)
            filename_population(dwelling_index)=read_name_char('filename_population(dwelling_index)',filename_population(dwelling_index),unit_in,unit_logfile)
            filename_population(population_index)=read_name_char('filename_population(population_index)',filename_population(population_index),unit_in,unit_logfile)
            filename_population(establishment_index)=read_name_char('filename_population(establishment_index)',filename_population(establishment_index),unit_in,unit_logfile)
            filename_population(school_index)=read_name_char('filename_population(school_index)',filename_population(school_index),unit_in,unit_logfile)
            filename_population(kindergaten_index)=read_name_char('filename_population(kindergaten_index)',filename_population(kindergaten_index),unit_in,unit_logfile)
            filename_population(home_index)=read_name_char('filename_population(home_index)',filename_population(home_index),unit_in,unit_logfile)
            filename_population(municipality_index)=read_name_char('filename_population(municipality_index)',filename_population(municipality_index),unit_in,unit_logfile)

            ! Pollen
            pathname_pollen(birch_proxy_index) = read_name_char( &
                "pathname_pollen(birch_proxy_index)", pathname_pollen(birch_proxy_index), unit_in, unit_logfile)
            filename_pollen(birch_proxy_index) = read_name_char( &
                "filename_pollen(birch_proxy_index)", filename_pollen(birch_proxy_index), unit_in, unit_logfile)
            

            pathname_receptor=read_name_char('pathname_receptor',pathname_receptor,unit_in,unit_logfile)
            filename_receptor=read_name_char('filename_receptor',filename_receptor,unit_in,unit_logfile)

            pathname_timeprofile=read_name_char('pathname_timeprofile',pathname_timeprofile,unit_in,unit_logfile)
            filename_timeprofile=read_name_char('filename_timeprofile',filename_timeprofile,unit_in,unit_logfile)

            population_data_type=read_name_integer('population_data_type',population_data_type,unit_in,unit_logfile)

            FF_min_dispersion=read_name_real('FF_min_dispersion',FF_min_dispersion,unit_in,unit_logfile)
            emission_timeprofile_hour_shift=read_name_integer('emission_timeprofile_hour_shift',emission_timeprofile_hour_shift,unit_in,unit_logfile)

            use_last_meteo_in_dispersion=read_name_logical('use_last_meteo_in_dispersion',use_last_meteo_in_dispersion,unit_in,unit_logfile)
            use_meandering_in_dispersion=read_name_logical('use_meandering_in_dispersion',use_meandering_in_dispersion,unit_in,unit_logfile)

            use_traffic_for_sigma0_flag=read_name_logical('use_traffic_for_sigma0_flag',use_traffic_for_sigma0_flag,unit_in,unit_logfile)
            !use_traffic_for_minFF_flag=read_name_logical('use_traffic_for_minFF_flag',use_traffic_for_minFF_flag,unit_in,unit_logfile)

            use_alternative_meteorology_flag=read_name_logical('use_alternative_meteorology_flag',use_alternative_meteorology_flag,unit_in,unit_logfile)
            ustar_min=read_name_real('ustar_min',ustar_min,unit_in,unit_logfile)
            hmix_min=read_name_real('hmix_min',hmix_min,unit_in,unit_logfile)
            hmix_max=read_name_real('hmix_max',hmix_max,unit_in,unit_logfile)
            use_alternative_z0_flag=read_name_logical('use_alternative_z0_flag',use_alternative_z0_flag,unit_in,unit_logfile)
            alternative_meteorology_type=read_name_char('alternative_meteorology_type',alternative_meteorology_type,unit_in,unit_logfile)


            !Read emission factors for traffic
            emission_factor(nox_index,traffic_index,:)=read_name_real('emission_factor(nox_index,traffic_index,:)',emission_factor(nox_index,traffic_index,1),unit_in,unit_logfile)
            emission_factor(nox_index,traffic_index,1)=read_name_real('emission_factor(nox_index,traffic_index,1)',emission_factor(nox_index,traffic_index,1),unit_in,unit_logfile)
            emission_factor(nox_index,traffic_index,2)=read_name_real('emission_factor(nox_index,traffic_index,2)',emission_factor(nox_index,traffic_index,2),unit_in,unit_logfile)
            emission_factor(no2_index,traffic_index,:)=read_name_real('emission_factor(no2_index,traffic_index,:)',emission_factor(no2_index,traffic_index,1),unit_in,unit_logfile)
            emission_factor(no2_index,traffic_index,1)=read_name_real('emission_factor(no2_index,traffic_index,1)',emission_factor(no2_index,traffic_index,1),unit_in,unit_logfile)
            emission_factor(no2_index,traffic_index,2)=read_name_real('emission_factor(no2_index,traffic_index,2)',emission_factor(no2_index,traffic_index,2),unit_in,unit_logfile)
            emission_factor(pm25_index,traffic_index,:)=read_name_real('emission_factor(pm25_index,traffic_index,:)',emission_factor(pm25_index,traffic_index,1),unit_in,unit_logfile)
            emission_factor(pm25_index,traffic_index,1)=read_name_real('emission_factor(pm25_index,traffic_index,1)',emission_factor(pm25_index,traffic_index,1),unit_in,unit_logfile)
            emission_factor(pm25_index,traffic_index,2)=read_name_real('emission_factor(pm25_index,traffic_index,2)',emission_factor(pm25_index,traffic_index,2),unit_in,unit_logfile)
            emission_factor(pm10_index,traffic_index,:)=read_name_real('emission_factor(pm10_index,traffic_index,:)',emission_factor(pm10_index,traffic_index,1),unit_in,unit_logfile)
            emission_factor(pm10_index,traffic_index,1)=read_name_real('emission_factor(pm10_index,traffic_index,1)',emission_factor(pm10_index,traffic_index,1),unit_in,unit_logfile)
            emission_factor(pm10_index,traffic_index,2)=read_name_real('emission_factor(pm10_index,traffic_index,2)',emission_factor(pm10_index,traffic_index,2),unit_in,unit_logfile)
            emission_factor(pmex_index,traffic_index,:)=read_name_real('emission_factor(pmex_index,traffic_index,:)',emission_factor(pmex_index,traffic_index,1),unit_in,unit_logfile)
            emission_factor(pmex_index,traffic_index,1)=read_name_real('emission_factor(pmex_index,traffic_index,1)',emission_factor(pmex_index,traffic_index,1),unit_in,unit_logfile)
            emission_factor(pmex_index,traffic_index,2)=read_name_real('emission_factor(pmex_index,traffic_index,2)',emission_factor(pmex_index,traffic_index,2),unit_in,unit_logfile)

            ratio_truck_car_emission(nox_index)=read_name_real('ratio_truck_car_emission(nox_index)',ratio_truck_car_emission(nox_index),unit_in,unit_logfile)
            ratio_truck_car_emission(pm25_index)=read_name_real('ratio_truck_car_emission(pm25_index)',ratio_truck_car_emission(pm25_index),unit_in,unit_logfile)
            ratio_truck_car_emission(pm10_index)=read_name_real('ratio_truck_car_emission(pm10_index)',ratio_truck_car_emission(pm10_index),unit_in,unit_logfile)
            !NO2 ratio does not do anything but is included for possible future changes
            ratio_truck_car_emission(no2_index)=read_name_real('ratio_truck_car_emission(no2_index)',ratio_truck_car_emission(no2_index),unit_in,unit_logfile)

            z_rec=read_name_real('z_rec',z_rec(allsource_index,1),unit_in,unit_logfile)

            replace_z0=read_name_real('replace_z0',replace_z0,unit_in,unit_logfile)
            replace_invL=read_name_real('replace_invL',replace_invL,unit_in,unit_logfile)
            replace_hmix=read_name_real('replace_hmix',replace_hmix,unit_in,unit_logfile)
            FF_scale=read_name_real('FF_scale',FF_scale,unit_in,unit_logfile)
            FF10_offset=read_name_real('FF10_offset',FF10_offset,unit_in,unit_logfile)
            DD_offset=read_name_real('DD_offset',DD_offset,unit_in,unit_logfile)
            J_scale=read_name_real('J_scale',J_scale,unit_in,unit_logfile)

            save_netcdf_file_flag=read_name_logical('save_netcdf_file_flag',save_netcdf_file_flag,unit_in,unit_logfile)
            save_netcdf_receptor_flag=read_name_logical('save_netcdf_receptor_flag',save_netcdf_receptor_flag,unit_in,unit_logfile)
            save_netcdf_fraction_as_contribution_flag=read_name_logical('save_netcdf_fraction_as_contribution_flag',save_netcdf_fraction_as_contribution_flag,unit_in,unit_logfile)

            calculate_tiling_flag=read_name_logical('calculate_tiling_flag',calculate_tiling_flag,unit_in,unit_logfile)
            calculate_region_tiling_flag=read_name_logical('calculate_region_tiling_flag',calculate_region_tiling_flag,unit_in,unit_logfile)

            pathname_region_id=read_name_char('pathname_region_id',pathname_region_id,unit_in,unit_logfile)
            filename_region_id=read_name_char('filename_region_id',filename_region_id,unit_in,unit_logfile)
            region_name=read_name_char('region_name',region_name,unit_in,unit_logfile)
            region_id=read_name_integer('region_id',region_id,unit_in,unit_logfile)
            region_index=read_name_integer('region_index',region_index,unit_in,unit_logfile)
            region_subgrid_delta=read_name_real('region_subgrid_delta',region_subgrid_delta,unit_in,unit_logfile)
            use_region_select_and_mask_flag=read_name_logical('use_region_select_and_mask_flag',use_region_select_and_mask_flag,unit_in,unit_logfile)

            max_interpolation_subgrid_size=read_name_real('max_interpolation_subgrid_size',max_interpolation_subgrid_size,unit_in,unit_logfile)

            pathname_tiles=read_name_char('pathname_tiles',pathname_tiles,unit_in,unit_logfile)
            filename_tiles=read_name_char('filename_tiles',filename_tiles,unit_in,unit_logfile)
            tile_tag=read_name_char('tile_tag',tile_tag,unit_in,unit_logfile)
            save_tile_tag=read_name_char('save_tile_tag',save_tile_tag,unit_in,unit_logfile)

            use_NORTRIP_emission_data=read_name_logical('use_NORTRIP_emission_data',use_NORTRIP_emission_data,unit_in,unit_logfile)
            use_NORTRIP_emission_pollutant(nox_index)=read_name_logical('use_NORTRIP_emission_pollutant(nox_index)',use_NORTRIP_emission_pollutant(nox_index),unit_in,unit_logfile)
            use_NORTRIP_emission_pollutant(pm10_index)=read_name_logical('use_NORTRIP_emission_pollutant(pm10_index)',use_NORTRIP_emission_pollutant(pm10_index),unit_in,unit_logfile)
            use_NORTRIP_emission_pollutant(pm25_index)=read_name_logical('use_NORTRIP_emission_pollutant(pm25_index)',use_NORTRIP_emission_pollutant(pm25_index),unit_in,unit_logfile)
            use_NORTRIP_emission_pollutant(pmex_index)=read_name_logical('use_NORTRIP_emission_pollutant(pmex_index)',use_NORTRIP_emission_pollutant(pmex_index),unit_in,unit_logfile)

            use_RWC_emission_data=read_name_logical('use_RWC_emission_data',use_RWC_emission_data,unit_in,unit_logfile)
            HDD_threshold_value=read_name_integer('HDD_threshold_value',HDD_threshold_value,unit_in,unit_logfile)
            DMT_min_value=read_name_real('DMT_min_value',DMT_min_value,unit_in,unit_logfile)
            inpath_region_heating_scaling=read_name_char('inpath_region_heating_scaling',inpath_region_heating_scaling,unit_in,unit_logfile)
            infile_region_heating_scaling=read_name_char('infile_region_heating_scaling',infile_region_heating_scaling,unit_in,unit_logfile)

            integral_subgrid_delta_ref=read_name_real('integral_subgrid_delta_ref',integral_subgrid_delta_ref,unit_in,unit_logfile)

            pathname_rl_change=read_name_char('pathname_rl_change',pathname_rl_change,unit_in,unit_logfile)
            filename_rl_change=read_name_char('filename_rl_change',filename_rl_change,unit_in,unit_logfile)

            forecast_hour_str=read_name_char('forecast_hour_str',forecast_hour_str,unit_in,unit_logfile)
            NORTRIP_hour_str=read_name_char('NORTRIP_hour_str',NORTRIP_hour_str,unit_in,unit_logfile)

            include_o3_in_aqi_index=read_name_logical('include_o3_in_aqi_index',include_o3_in_aqi_index,unit_in,unit_logfile)

            n_kz_iterations=read_name_integer('n_kz_iterations',n_kz_iterations,unit_in,unit_logfile)

            read_weekly_shipping_data_flag=read_name_logical('read_weekly_shipping_data_flag',read_weekly_shipping_data_flag,unit_in,unit_logfile)
            read_monthly_and_daily_shipping_data_flag=read_name_logical('read_monthly_and_daily_shipping_data_flag',read_monthly_and_daily_shipping_data_flag,unit_in,unit_logfile)

            use_tunnel_emissions_flag=read_name_logical('use_tunnel_emissions_flag',use_tunnel_emissions_flag,unit_in,unit_logfile)
            use_tunnel_deposition_flag=read_name_logical('use_tunnel_deposition_flag',use_tunnel_deposition_flag,unit_in,unit_logfile)
            ventilation_factor=read_name_real('ventilation_factor',ventilation_factor,unit_in,unit_logfile)
            min_ADT_ventilation_factor=read_name_real('min_ADT_ventilation_factor',min_ADT_ventilation_factor,unit_in,unit_logfile)
            min_length_ventilation_factor=read_name_real('min_length_ventilation_factor',min_length_ventilation_factor,unit_in,unit_logfile)
            windspeed_tunnel=read_name_real('windspeed_tunnel',windspeed_tunnel,unit_in,unit_logfile)


            save_emissions_for_EMEP(traffic_index)=read_name_logical('save_emissions_for_EMEP(traffic_index)',save_emissions_for_EMEP(traffic_index),unit_in,unit_logfile)
            save_emissions_for_EMEP(shipping_index)=read_name_logical('save_emissions_for_EMEP(shipping_index)',save_emissions_for_EMEP(shipping_index),unit_in,unit_logfile)
            save_emissions_for_EMEP(industry_index)=read_name_logical('save_emissions_for_EMEP(industry_index)',save_emissions_for_EMEP(industry_index),unit_in,unit_logfile)
            save_emissions_for_EMEP(heating_index)=read_name_logical('save_emissions_for_EMEP(heating_index)',save_emissions_for_EMEP(heating_index),unit_in,unit_logfile)
            save_emissions_for_EMEP(agriculture_index)=read_name_logical('save_emissions_for_EMEP(agriculture_index)',save_emissions_for_EMEP(agriculture_index),unit_in,unit_logfile)
            !Additional GNFR sources
            save_emissions_for_EMEP(publicpower_index)=read_name_logical('save_emissions_for_EMEP(publicpower_index)',save_emissions_for_EMEP(publicpower_index),unit_in,unit_logfile)
            save_emissions_for_EMEP(fugitive_index)=read_name_logical('save_emissions_for_EMEP(fugitive_index)',save_emissions_for_EMEP(fugitive_index),unit_in,unit_logfile)
            save_emissions_for_EMEP(solvents_index)=read_name_logical('save_emissions_for_EMEP(solvents_index)',save_emissions_for_EMEP(solvents_index),unit_in,unit_logfile)
            save_emissions_for_EMEP(aviation_index)=read_name_logical('save_emissions_for_EMEP(aviation_index)',save_emissions_for_EMEP(aviation_index),unit_in,unit_logfile)
            save_emissions_for_EMEP(offroad_index)=read_name_logical('save_emissions_for_EMEP(offroad_index)',save_emissions_for_EMEP(offroad_index),unit_in,unit_logfile)
            save_emissions_for_EMEP(waste_index)=read_name_logical('save_emissions_for_EMEP(waste_index)',save_emissions_for_EMEP(waste_index),unit_in,unit_logfile)
            save_emissions_for_EMEP(livestock_index)=read_name_logical('save_emissions_for_EMEP(livestock_index)',save_emissions_for_EMEP(livestock_index),unit_in,unit_logfile)
            save_emissions_for_EMEP(other_index)=read_name_logical('save_emissions_for_EMEP(other_index)',save_emissions_for_EMEP(other_index),unit_in,unit_logfile)
            !Set all source index to true if any of the sources are to be saved. allsource_index defines if the routine is called or not
            do i_source=1,n_source_index
                if (save_emissions_for_EMEP(i_source)) save_emissions_for_EMEP(allsource_index)=.true.
            enddo

            pathname_emissions_for_EMEP=read_name_char('pathname_emissions_for_EMEP',pathname_emissions_for_EMEP,unit_in,unit_logfile)

            save_emissions_start_index=read_name_integer('save_emissions_start_index',save_emissions_start_index,unit_in,unit_logfile)
            save_emissions_end_index=read_name_integer('save_emissions_end_index',save_emissions_end_index,unit_in,unit_logfile)

            save_emissions_for_EMEP_projection=read_name_char('save_emissions_for_EMEP_projection',save_emissions_for_EMEP_projection,unit_in,unit_logfile)
            save_emissions_for_EMEP_region=read_name_char('save_emissions_for_EMEP_region',save_emissions_for_EMEP_region,unit_in,unit_logfile)

            save_compounds=read_name_logical('save_compounds',save_compounds,unit_in,unit_logfile)
            save_source_contributions=read_name_logical('save_source_contributions',save_source_contributions,unit_in,unit_logfile)
            save_emep_source_contributions=read_name_logical('save_emep_source_contributions',save_emep_source_contributions,unit_in,unit_logfile)
            save_emep_additional_source_contributions=read_name_logical('save_emep_additional_source_contributions',save_emep_additional_source_contributions,unit_in,unit_logfile)
            save_total_source_contributions=read_name_logical('save_total_source_contributions',save_total_source_contributions,unit_in,unit_logfile)
            save_local_source_contributions_from_in_region=read_name_logical('save_local_source_contributions_from_in_region',save_local_source_contributions_from_in_region,unit_in,unit_logfile)
            save_semilocal_source_contributions_from_in_region=read_name_logical('save_semilocal_source_contributions_from_in_region',save_semilocal_source_contributions_from_in_region,unit_in,unit_logfile)
            save_total_source_contributions_from_in_region=read_name_logical('save_total_source_contributions_from_in_region',save_total_source_contributions_from_in_region,unit_in,unit_logfile)
            save_no2_source_contributions=read_name_logical('save_no2_source_contributions',save_no2_source_contributions,unit_in,unit_logfile)
            save_o3_source_contributions=read_name_logical('save_o3_source_contributions',save_o3_source_contributions,unit_in,unit_logfile)
            save_wind_vectors=read_name_logical('save_wind_vectors',save_wind_vectors,unit_in,unit_logfile)
            save_other_meteo=read_name_logical('save_other_meteo',save_other_meteo,unit_in,unit_logfile)
            save_emep_original=read_name_logical('save_emep_original',save_emep_original,unit_in,unit_logfile)
            save_emissions=read_name_logical('save_emissions',save_emissions,unit_in,unit_logfile)
            save_for_chemistry=read_name_logical('save_for_chemistry',save_for_chemistry,unit_in,unit_logfile)
            save_population=read_name_logical('save_population',save_population,unit_in,unit_logfile)
            save_aqi=read_name_logical('save_aqi',save_aqi,unit_in,unit_logfile)
            save_emep_species=read_name_logical('save_emep_species',save_emep_species,unit_in,unit_logfile)
            save_deposition=read_name_logical('save_deposition',save_deposition,unit_in,unit_logfile)
            save_seasalt=read_name_logical('save_seasalt',save_seasalt,unit_in,unit_logfile)
            save_pollen=read_name_logical("save_pollen", save_pollen, unit_in, unit_logfile)

            lowest_stable_L=read_name_real('lowest_stable_L',lowest_stable_L,unit_in,unit_logfile)
            lowest_unstable_L=read_name_real('lowest_unstable_L',lowest_unstable_L,unit_in,unit_logfile)

            tunnel_sig_z_00=read_name_real('tunnel_sig_z_00',tunnel_sig_z_00,unit_in,unit_logfile)
            bridge_h_emis=read_name_real('bridge_h_emis',bridge_h_emis,unit_in,unit_logfile)

            !Input variable names for meteo data read from EMEP
            var_name_nc(hmix_nc_index,all_nc_index,allsource_nc_index)=read_name_char('var_name_nc(hmix_nc_index)',var_name_nc(hmix_nc_index,all_nc_index,allsource_nc_index),unit_in,unit_logfile)
            var_name_nc(u10_nc_index,all_nc_index,allsource_nc_index)=read_name_char('var_name_nc(u10_nc_index)',var_name_nc(u10_nc_index,all_nc_index,allsource_nc_index),unit_in,unit_logfile)
            var_name_nc(v10_nc_index,all_nc_index,allsource_nc_index)=read_name_char('var_name_nc(v10_nc_index)',var_name_nc(v10_nc_index,all_nc_index,allsource_nc_index),unit_in,unit_logfile)
            var_name_nc(FF10_nc_index,all_nc_index,allsource_nc_index)=read_name_char('var_name_nc(FF10_nc_index)',var_name_nc(FF10_nc_index,all_nc_index,allsource_nc_index),unit_in,unit_logfile)
            var_name_nc(ugrid_nc_index,all_nc_index,allsource_nc_index)=read_name_char('var_name_nc(ugrid_nc_index)',var_name_nc(ugrid_nc_index,all_nc_index,allsource_nc_index),unit_in,unit_logfile)
            var_name_nc(vgrid_nc_index,all_nc_index,allsource_nc_index)=read_name_char('var_name_nc(vgrid_nc_index)',var_name_nc(vgrid_nc_index,all_nc_index,allsource_nc_index),unit_in,unit_logfile)
            var_name_nc(FFgrid_nc_index,all_nc_index,allsource_nc_index)=read_name_char('var_name_nc(FFgrid_nc_index)',var_name_nc(FFgrid_nc_index,all_nc_index,allsource_nc_index),unit_in,unit_logfile)
            var_name_nc(inv_FFgrid_nc_index,all_nc_index,allsource_nc_index)=read_name_char('var_name_nc(inv_FFgrid_nc_index)',var_name_nc(inv_FFgrid_nc_index,all_nc_index,allsource_nc_index),unit_in,unit_logfile)
            var_name_nc(inv_FF10_nc_index,all_nc_index,allsource_nc_index)=read_name_char('var_name_nc(inv_FF10_nc_index)',var_name_nc(inv_FF10_nc_index,all_nc_index,allsource_nc_index),unit_in,unit_logfile)
            var_name_nc(kz_nc_index,all_nc_index,allsource_nc_index)=read_name_char('var_name_nc(kz_nc_index)',var_name_nc(kz_nc_index,all_nc_index,allsource_nc_index),unit_in,unit_logfile)
            var_name_nc(ustar_nc_index,all_nc_index,allsource_nc_index)=read_name_char('var_name_nc(ustar_nc_index)',var_name_nc(ustar_nc_index,all_nc_index,allsource_nc_index),unit_in,unit_logfile)
            var_name_nc(logz0_nc_index,all_nc_index,allsource_nc_index)=read_name_char('var_name_nc(logz0_nc_index)',var_name_nc(logz0_nc_index,all_nc_index,allsource_nc_index),unit_in,unit_logfile)
            var_name_nc(invL_nc_index,all_nc_index,allsource_nc_index)=read_name_char('var_name_nc(invL_nc_index)',var_name_nc(invL_nc_index,all_nc_index,allsource_nc_index),unit_in,unit_logfile)
            var_name_nc(ZTOP_nc_index,all_nc_index,allsource_nc_index)=read_name_char('var_name_nc(ZTOP_nc_index)',var_name_nc(ZTOP_nc_index,all_nc_index,allsource_nc_index),unit_in,unit_logfile)
            var_name_nc(t2m_nc_index,all_nc_index,allsource_nc_index)=read_name_char('var_name_nc(t2m_nc_index)',var_name_nc(t2m_nc_index,all_nc_index,allsource_nc_index),unit_in,unit_logfile)
            var_name_nc(precip_nc_index,all_nc_index,allsource_nc_index)=read_name_char('var_name_nc(precip_nc_index)',var_name_nc(precip_nc_index,all_nc_index,allsource_nc_index),unit_in,unit_logfile)
            var_name_nc(J_nc_index,all_nc_index,allsource_nc_index)=read_name_char('var_name_nc(J_nc_index)',var_name_nc(J_nc_index,all_nc_index,allsource_nc_index),unit_in,unit_logfile)
            var_name_nc(phi_nc_index,all_nc_index,allsource_nc_index)=read_name_char('var_name_nc(phi_nc_index)',var_name_nc(phi_nc_index,all_nc_index,allsource_nc_index),unit_in,unit_logfile)

            save_netcdf_average_flag=read_name_logical('save_netcdf_average_flag',save_netcdf_average_flag,unit_in,unit_logfile)

            use_traffic_nox_emission_temperature_dependency=read_name_logical('use_traffic_nox_emission_temperature_dependency',use_traffic_nox_emission_temperature_dependency,unit_in,unit_logfile)
            traffic_nox_emission_temperature_ref_temperature(1)=read_name_real('traffic_nox_emission_temperature_ref_temperature(1)',traffic_nox_emission_temperature_ref_temperature(1),unit_in,unit_logfile)
            traffic_nox_emission_temperature_ref_temperature(2)=read_name_real('traffic_nox_emission_temperature_ref_temperature(2)',traffic_nox_emission_temperature_ref_temperature(2),unit_in,unit_logfile)
            traffic_nox_emission_temperature_ref_scaling(1)=read_name_real('traffic_nox_emission_temperature_ref_scaling(1)',traffic_nox_emission_temperature_ref_scaling(1),unit_in,unit_logfile)
            traffic_nox_emission_temperature_ref_scaling(2)=read_name_real('traffic_nox_emission_temperature_ref_scaling(2)',traffic_nox_emission_temperature_ref_scaling(2),unit_in,unit_logfile)

            calculate_deposition_flag=read_name_logical('calculate_deposition_flag',calculate_deposition_flag,unit_in,unit_logfile)
            calculate_source_depletion_flag=read_name_logical('calculate_source_depletion_flag',calculate_source_depletion_flag,unit_in,unit_logfile)
            read_landuse_flag=read_name_logical('read_landuse_flag',read_landuse_flag,unit_in,unit_logfile)
            use_plume_dispersion_deposition_flag=read_name_logical('use_plume_dispersion_deposition_flag',use_plume_dispersion_deposition_flag,unit_in,unit_logfile)

            pathname_landuse=read_name_char('pathname_landuse',pathname_landuse,unit_in,unit_logfile)
            filename_landuse=read_name_char('filename_landuse',filename_landuse,unit_in,unit_logfile)

            adjust_wetdepo_integral_to_lowest_layer_flag=read_name_logical('adjust_wetdepo_integral_to_lowest_layer_flag',adjust_wetdepo_integral_to_lowest_layer_flag,unit_in,unit_logfile)

            auto_adjustment_for_summertime=read_name_logical('auto_adjustment_for_summertime',auto_adjustment_for_summertime,unit_in,unit_logfile)

            use_EMEP_surface_ozone_flag=read_name_logical('use_EMEP_surface_ozone_flag',use_EMEP_surface_ozone_flag,unit_in,unit_logfile)
            use_EMEP_surface_compounds_flag=read_name_logical('use_EMEP_surface_compounds_flag',use_EMEP_surface_compounds_flag,unit_in,unit_logfile)
            use_water_in_EMEP_surface_pm_flag=read_name_logical('use_water_in_EMEP_surface_pm_flag',use_water_in_EMEP_surface_pm_flag,unit_in,unit_logfile)

            save_compounds_as_ascii=read_name_logical('save_compounds_as_ascii',save_compounds_as_ascii,unit_in,unit_logfile)

            use_GNFR_emissions_from_EMEP_flag=read_name_logical('use_GNFR_emissions_from_EMEP_flag',use_GNFR_emissions_from_EMEP_flag,unit_in,unit_logfile)
            use_GNFR19_emissions_from_EMEP_flag=read_name_logical('use_GNFR19_emissions_from_EMEP_flag',use_GNFR19_emissions_from_EMEP_flag,unit_in,unit_logfile)
            use_alphabetic_GNFR_emissions_from_EMEP_flag=read_name_logical('use_alphabetic_GNFR_emissions_from_EMEP_flag',use_alphabetic_GNFR_emissions_from_EMEP_flag,unit_in,unit_logfile)

            use_emission_naming_template_flag=read_name_logical('use_emission_naming_template_flag',use_emission_naming_template_flag,unit_in,unit_logfile)
            emission_naming_template_str=read_name_char('emission_naming_template_str',emission_naming_template_str,unit_in,unit_logfile)

            read_OSM_roadlink_data_flag=read_name_logical('read_OSM_roadlink_data_flag',read_OSM_roadlink_data_flag,unit_in,unit_logfile)
            no_header_roadlink_data_flag=read_name_logical('no_header_roadlink_data_flag',no_header_roadlink_data_flag,unit_in,unit_logfile)

            EMEP_surface_level_nc=read_name_integer('EMEP_surface_level_nc',EMEP_surface_level_nc,unit_in,unit_logfile)
            EMEP_surface_level_nc_2=read_name_integer('EMEP_surface_level_nc_2',EMEP_surface_level_nc_2,unit_in,unit_logfile)

            limit_industry_delta=read_name_real('limit_industry_delta',limit_industry_delta,unit_in,unit_logfile)
            limit_shipping_delta=read_name_real('limit_shipping_delta',limit_shipping_delta,unit_in,unit_logfile)
            limit_heating_delta=read_name_real('limit_heating_delta',limit_heating_delta,unit_in,unit_logfile)
            limit_population_delta=read_name_real('limit_population_delta',limit_population_delta,unit_in,unit_logfile)
            limit_pollen_delta=read_name_real("limit_pollen_delta", limit_pollen_delta, unit_in, unit_logfile)

            use_user_specified_sectors_flag=read_name_logical('use_user_specified_sectors_flag',use_user_specified_sectors_flag,unit_in,unit_logfile)
            if (use_user_specified_sectors_flag) then
                uEMEP_to_EMEP_replace_sector(traffic_index)=read_name_integer('uEMEP_to_EMEP_replace_sector(traffic_index)',uEMEP_to_EMEP_replace_sector(traffic_index),unit_in,unit_logfile)
                uEMEP_to_EMEP_replace_sector(shipping_index)=read_name_integer('uEMEP_to_EMEP_replace_sector(shipping_index)',uEMEP_to_EMEP_replace_sector(shipping_index),unit_in,unit_logfile)
                uEMEP_to_EMEP_replace_sector(agriculture_index)=read_name_integer('uEMEP_to_EMEP_replace_sector(agriculture_index)',uEMEP_to_EMEP_replace_sector(agriculture_index),unit_in,unit_logfile)
                uEMEP_to_EMEP_replace_sector(heating_index)=read_name_integer('uEMEP_to_EMEP_replace_sector(heating_index)',uEMEP_to_EMEP_replace_sector(heating_index),unit_in,unit_logfile)
                uEMEP_to_EMEP_replace_sector(industry_index)=read_name_integer('uEMEP_to_EMEP_replace_sector(industry_index)',uEMEP_to_EMEP_replace_sector(industry_index),unit_in,unit_logfile)
                !Additional GNFR sources
                uEMEP_to_EMEP_replace_sector(publicpower_index)=read_name_integer('uEMEP_to_EMEP_replace_sector(publicpower_index)',uEMEP_to_EMEP_replace_sector(publicpower_index),unit_in,unit_logfile)
                uEMEP_to_EMEP_replace_sector(fugitive_index)=read_name_integer('uEMEP_to_EMEP_replace_sector(fugitive_index)',uEMEP_to_EMEP_replace_sector(fugitive_index),unit_in,unit_logfile)
                uEMEP_to_EMEP_replace_sector(solvents_index)=read_name_integer('uEMEP_to_EMEP_replace_sector(solvents_index)',uEMEP_to_EMEP_replace_sector(solvents_index),unit_in,unit_logfile)
                uEMEP_to_EMEP_replace_sector(aviation_index)=read_name_integer('uEMEP_to_EMEP_replace_sector(aviation_index)',uEMEP_to_EMEP_replace_sector(aviation_index),unit_in,unit_logfile)
                uEMEP_to_EMEP_replace_sector(offroad_index)=read_name_integer('uEMEP_to_EMEP_replace_sector(offroad_index)',uEMEP_to_EMEP_replace_sector(offroad_index),unit_in,unit_logfile)
                uEMEP_to_EMEP_replace_sector(waste_index)=read_name_integer('uEMEP_to_EMEP_replace_sector(waste_index)',uEMEP_to_EMEP_replace_sector(waste_index),unit_in,unit_logfile)
                uEMEP_to_EMEP_replace_sector(livestock_index)=read_name_integer('uEMEP_to_EMEP_replace_sector(livestock_index)',uEMEP_to_EMEP_replace_sector(livestock_index),unit_in,unit_logfile)
                uEMEP_to_EMEP_replace_sector(other_index)=read_name_integer('uEMEP_to_EMEP_replace_sector(other_index)',uEMEP_to_EMEP_replace_sector(other_index),unit_in,unit_logfile)
            endif

            EMEP_emission_aggregation_period=read_name_real('EMEP_emission_aggregation_period',EMEP_emission_aggregation_period,unit_in,unit_logfile)
            read_population_from_netcdf_flag=read_name_logical('read_population_from_netcdf_flag',read_population_from_netcdf_flag,unit_in,unit_logfile)
            read_population_from_netcdf_local_flag=read_name_logical('read_population_from_netcdf_local_flag',read_population_from_netcdf_flag,unit_in,unit_logfile)

            auto_select_OSM_country_flag=read_name_logical('auto_select_OSM_country_flag',auto_select_OSM_country_flag,unit_in,unit_logfile)
            pathname_boundingbox=read_name_char('pathname_boundingbox',pathname_boundingbox,unit_in,unit_logfile)
            filename_boundingbox=read_name_char('filename_boundingbox',filename_boundingbox,unit_in,unit_logfile)
            select_country_by_name=read_name_char('select_country_by_name',select_country_by_name,unit_in,unit_logfile)

            select_latlon_centre_domain_position_flag=read_name_logical('select_latlon_centre_domain_position_flag',select_latlon_centre_domain_position_flag,unit_in,unit_logfile)
            select_lat_centre_position=read_name_real('select_lat_centre_position',select_lat_centre_position,unit_in,unit_logfile)
            select_lon_centre_position=read_name_real('select_lon_centre_position',select_lon_centre_position,unit_in,unit_logfile)
            select_domain_width_EW_km=read_name_real('select_domain_width_EW_km',select_domain_width_EW_km,unit_in,unit_logfile)
            select_domain_height_NS_km=read_name_real('select_domain_height_NS_km',select_domain_height_NS_km,unit_in,unit_logfile)

            osm_adt_power_scale=read_name_real('osm_adt_power_scale',osm_adt_power_scale,unit_in,unit_logfile)

            romberg_parameters(1)=read_name_real('romberg_parameters(1)',romberg_parameters(1),unit_in,unit_logfile)
            romberg_parameters(2)=read_name_real('romberg_parameters(2)',romberg_parameters(2),unit_in,unit_logfile)
            romberg_parameters(3)=read_name_real('romberg_parameters(3)',romberg_parameters(3),unit_in,unit_logfile)
            SRM_parameters(1)=read_name_real('SRM_parameters(1)',SRM_parameters(1),unit_in,unit_logfile)!beta
            SRM_parameters(2)=read_name_real('SRM_parameters(2)',SRM_parameters(2),unit_in,unit_logfile)!K
            SRM_parameters(3)=read_name_real('SRM_parameters(3)',SRM_parameters(3),unit_in,unit_logfile)!F

            sig_y_scaling_factor=read_name_real('sig_y_scaling_factor',sig_y_scaling_factor,unit_in,unit_logfile)

            read_shipping_from_netcdf_flag=read_name_logical('read_shipping_from_netcdf_flag',read_shipping_from_netcdf_flag,unit_in,unit_logfile)
            min_proxy_emission_shipping_value=read_name_real('min_proxy_emission_shipping_value',min_proxy_emission_shipping_value,unit_in,unit_logfile)

            population_power_scale=read_name_real('population_power_scale',population_power_scale,unit_in,unit_logfile)
            H_emep=read_name_real('H_emep',H_emep,unit_in,unit_logfile)

            !Allow the user to change the EMEP PM used for nonlocal contribution.
            !Will be overridden by use_EMEP_surface_compounds_flag and use_water_in_EMEP_surface_pm_flag if they are set to true
            comp_name_nc(pm10_nc_index)=read_name_char('comp_name_nc(pm10_nc_index)',comp_name_nc(pm10_nc_index),unit_in,unit_logfile)
            comp_name_nc(pm25_nc_index)=read_name_char('comp_name_nc(pm25_nc_index)',comp_name_nc(pm25_nc_index),unit_in,unit_logfile)
            comp_name_nc(o3_nc_index)=read_name_char('comp_name_nc(o3_nc_index)',comp_name_nc(o3_nc_index),unit_in,unit_logfile)
            comp_name_nc(no2_nc_index)=read_name_char('comp_name_nc(no2_nc_index)',comp_name_nc(no2_nc_index),unit_in,unit_logfile)
            comp_name_nc(nox_nc_index)=read_name_char('comp_name_nc(nox_nc_index)',comp_name_nc(nox_nc_index),unit_in,unit_logfile)
            comp_name_nc(nh3_nc_index)=read_name_char('comp_name_nc(nh3_nc_index)',comp_name_nc(nh3_nc_index),unit_in,unit_logfile)
            comp_name_nc(nh4_nc_index)=read_name_char('comp_name_nc(nh4_nc_index)',comp_name_nc(nh4_nc_index),unit_in,unit_logfile)
            comp_name_nc(pmex_nc_index)=read_name_char('comp_name_nc(pmex_nc_index)',comp_name_nc(pmex_nc_index),unit_in,unit_logfile)
            comp_name_nc(co_nc_index)=read_name_char('comp_name_nc(co_nc_index)',comp_name_nc(co_nc_index),unit_in,unit_logfile)
            comp_name_nc(bap_nc_index)=read_name_char('comp_name_nc(bap_nc_index)',comp_name_nc(bap_nc_index),unit_in,unit_logfile)
            comp_name_nc(c6h6_nc_index)=read_name_char('comp_name_nc(c6h6_nc_index)',comp_name_nc(c6h6_nc_index),unit_in,unit_logfile)
            comp_name_nc(so2_nc_index)=read_name_char('comp_name_nc(so2_nc_index)',comp_name_nc(so2_nc_index),unit_in,unit_logfile)

            comp_name_nc(somo35_nc_index)=read_name_char('comp_name_nc(somo35_nc_index)',comp_name_nc(somo35_nc_index),unit_in,unit_logfile)
            comp_name_nc(comax_nc_index)=read_name_char('comp_name_nc(comax_nc_index)',comp_name_nc(comax_nc_index),unit_in,unit_logfile)
            comp_name_nc(o3max_nc_index)=read_name_char('comp_name_nc(o3max_nc_index)',comp_name_nc(o3max_nc_index),unit_in,unit_logfile)
            comp_name_nc(o3_26th_nc_index)=read_name_char('comp_name_nc(o3_26th_nc_index)',comp_name_nc(o3_26th_nc_index),unit_in,unit_logfile)

            read_RWC_file_with_extra_HDD=read_name_logical('read_RWC_file_with_extra_HDD',read_RWC_file_with_extra_HDD,unit_in,unit_logfile)
            read_RWC_file_with_extra_HDD_and_height=read_name_logical('read_RWC_file_with_extra_HDD_and_height',read_RWC_file_with_extra_HDD_and_height,unit_in,unit_logfile)

            !Allows a scaling of EMEP input ozone. For testing.
            comp_scale_nc(o3_nc_index)=read_name_real('comp_scale_nc(o3_nc_index)',comp_scale_nc(o3_nc_index),unit_in,unit_logfile)
            comp_scale_nc(nox_nc_index)=read_name_real('comp_scale_nc(nox_nc_index)',comp_scale_nc(nox_nc_index),unit_in,unit_logfile)
            comp_scale_nc(no2_nc_index)=read_name_real('comp_scale_nc(no2_nc_index)',comp_scale_nc(no2_nc_index),unit_in,unit_logfile)

            use_alternative_traveltime_weighting=read_name_logical('use_alternative_traveltime_weighting',use_alternative_traveltime_weighting,unit_in,unit_logfile)
            traveltime_power=read_name_real('traveltime_power',traveltime_power,unit_in,unit_logfile)
            traveltime_scaling=read_name_real('traveltime_scaling',traveltime_scaling,unit_in,unit_logfile)
            use_straightline_traveltime_distance=read_name_logical('use_straightline_traveltime_distance',use_straightline_traveltime_distance,unit_in,unit_logfile)

            !Name of the netcdf variable read for population or dwelling proxy
            var_name_population_nc(population_nc_index)=read_name_char('var_name_population_nc(population_nc_index)',var_name_population_nc(population_nc_index),unit_in,unit_logfile)
            var_name_population_nc(dwelling_nc_index)=read_name_char('var_name_population_nc(dwelling_nc_index)',var_name_population_nc(dwelling_nc_index),unit_in,unit_logfile)

            ! Pollen
            var_name_pollen_nc(birch_proxy_index) = read_name_char( &
                "var_name_pollen_nc(birch_proxy_index)", var_name_pollen_nc(birch_proxy_index), unit_in, unit_logfile)

            f_no2_emep=read_name_real('f_no2_emep',f_no2_emep,unit_in,unit_logfile)

            limit_emep_grid_interpolation_region_to_calculation_region=read_name_logical('limit_emep_grid_interpolation_region_to_calculation_region',limit_emep_grid_interpolation_region_to_calculation_region,unit_in,unit_logfile)

            use_local_fraction_naming_template_flag=read_name_logical('use_local_fraction_naming_template_flag',use_local_fraction_naming_template_flag,unit_in,unit_logfile)
            use_local_fraction_grid_size_in_template_flag=read_name_logical('use_local_fraction_grid_size_in_template_flag',use_local_fraction_grid_size_in_template_flag,unit_in,unit_logfile)

            local_fraction_grid_size(1)=read_name_integer('local_fraction_grid_size(1)',local_fraction_grid_size(1),unit_in,unit_logfile)
            local_fraction_grid_size(2)=read_name_integer('local_fraction_grid_size(2)',local_fraction_grid_size(2),unit_in,unit_logfile)
            local_fraction_grid_size(3)=read_name_integer('local_fraction_grid_size(3)',local_fraction_grid_size(3),unit_in,unit_logfile)
            n_local_fraction_grids=read_name_integer('n_local_fraction_grids',n_local_fraction_grids,unit_in,unit_logfile)

            local_fraction_naming_template_str=read_name_char('local_fraction_naming_template_str',local_fraction_naming_template_str,unit_in,unit_logfile)

            local_fraction_grid_for_EMEP_grid_interpolation=read_name_integer('local_fraction_grid_for_EMEP_grid_interpolation',local_fraction_grid_for_EMEP_grid_interpolation,unit_in,unit_logfile)
            local_fraction_grid_for_EMEP_additional_grid_interpolation=read_name_integer('local_fraction_grid_for_EMEP_additional_grid_interpolation',local_fraction_grid_for_EMEP_additional_grid_interpolation,unit_in,unit_logfile)

            save_traffic_emissions_for_EMEP_as_exhaust_nonexhaust_flag=read_name_logical('save_traffic_emissions_for_EMEP_as_exhaust_nonexhaust_flag',save_traffic_emissions_for_EMEP_as_exhaust_nonexhaust_flag,unit_in,unit_logfile)

            n_var_av=read_name_integer('n_var_av',n_var_av,unit_in,unit_logfile)

            finished_filename=read_name_char('finished_filename',finished_filename,unit_in,unit_logfile)
            finished_subpath=read_name_char('finished_subpath',finished_subpath,unit_in,unit_logfile)

            use_annual_mean_pdf_chemistry_correction=read_name_logical('use_annual_mean_pdf_chemistry_correction',use_annual_mean_pdf_chemistry_correction,unit_in,unit_logfile)
            quick_annual_mean_pdf_chemistry_correction=read_name_logical('quick_annual_mean_pdf_chemistry_correction',quick_annual_mean_pdf_chemistry_correction,unit_in,unit_logfile)
            ox_sigma_ratio_pdf=read_name_real('ox_sigma_ratio_pdf',ox_sigma_ratio_pdf,unit_in,unit_logfile)
            nox_sigma_ratio_pdf=read_name_real('nox_sigma_ratio_pdf',nox_sigma_ratio_pdf,unit_in,unit_logfile)
            min_bin_pdf = read_name_real('min_bin_pdf', min_bin_pdf, unit_in, unit_logfile)
            max_bin_pdf=read_name_real('max_bin_pdf',max_bin_pdf,unit_in,unit_logfile)
            log10_step_bin_pdf=read_name_real('log10_step_bin_pdf',log10_step_bin_pdf,unit_in,unit_logfile)

            use_landuse_as_proxy=read_name_logical('use_landuse_as_proxy',use_landuse_as_proxy,unit_in,unit_logfile)
            read_rivm_landuse_flag=read_name_logical('read_rivm_landuse_flag',read_rivm_landuse_flag,unit_in,unit_logfile)
            var_name_landuse_nc(num_var_landuse_nc)=read_name_char('var_name_landuse_nc',var_name_landuse_nc(num_var_landuse_nc),unit_in,unit_logfile)
            use_rivm_agricuture_emission_data=read_name_logical('use_rivm_agricuture_emission_data',use_rivm_agricuture_emission_data,unit_in,unit_logfile)
            read_subgrid_emission_data=read_name_logical('read_subgrid_emission_data',read_subgrid_emission_data,unit_in,unit_logfile)
            use_rivm_subgrid_emission_format=read_name_logical('use_rivm_subgrid_emission_format',use_rivm_subgrid_emission_format,unit_in,unit_logfile)

            !Read landuse weighting this may take some time
            !Source input is numbered as GNFR13 in input but is placed in the uEMEP source sectors
            do i_source=1,n_source_index
                do i_landuse=1,n_clc_landuse_index

                    write(UNIT=a_str,FMT=*) convert_uEMEP_to_GNFR_sector_index(i_source)
                    write(UNIT=b_str,FMT=*) i_landuse
                    temp_str='landuse_proxy_weighting('//trim(adjustl(a_str))//','//trim(adjustl(b_str))//')'
                    !write(*,*) i_source,trim(temp_str)
                    landuse_proxy_weighting(i_source,i_landuse)=read_name_real(trim(temp_str),landuse_proxy_weighting(i_source,i_landuse),unit_in,unit_logfile)

                enddo
            enddo

            !Scale EMEP emission by sector
            !Source input is numbered as GNFR13 in input but is placed in the uEMEP source sectors
            do i_source=1,n_source_index

                write(UNIT=a_str,FMT=*) convert_uEMEP_to_GNFR_sector_index(i_source)
                temp_str='scale_GNFR_emission_source('//trim(adjustl(a_str))//')'
                !write(*,*) i_source,trim(temp_str)
                scale_GNFR_emission_source(i_source)=read_name_real(trim(temp_str),scale_GNFR_emission_source(i_source),unit_in,unit_logfile)

            enddo

            !Save original EMEP values for
            save_EMEP_somo35=read_name_logical('save_EMEP_somo35',save_EMEP_somo35,unit_in,unit_logfile)
            save_EMEP_comax=read_name_logical('save_EMEP_comax',save_EMEP_comax,unit_in,unit_logfile)
            save_EMEP_o3max=read_name_logical('save_EMEP_o3max',save_EMEP_o3max,unit_in,unit_logfile)
            save_EMEP_o3_26th=read_name_logical('save_EMEP_o3_26th',save_EMEP_o3_26th,unit_in,unit_logfile)
            save_EMEP_so2=read_name_logical('save_EMEP_so2',save_EMEP_so2,unit_in,unit_logfile)

            !Read subgrid receptor offset. This is for testing purposes only and applies only to the receptor subgrids
            subgrid_receptor_offset(x_dim_index)=read_name_real('subgrid_receptor_offset(x_dim_index)',subgrid_receptor_offset(x_dim_index),unit_in,unit_logfile)
            subgrid_receptor_offset(y_dim_index)=read_name_real('subgrid_receptor_offset(y_dim_index)',subgrid_receptor_offset(y_dim_index),unit_in,unit_logfile)

            derive_SOA_from_other_species=read_name_logical('derive_SOA_from_other_species',derive_SOA_from_other_species,unit_in,unit_logfile)

            Kz_scheme=read_name_integer('Kz_scheme',Kz_scheme,unit_in,unit_logfile)
            use_phi_for_invL=read_name_logical('use_phi_for_invL',use_phi_for_invL,unit_in,unit_logfile)
            z_invL=read_name_real('z_invL',z_invL,unit_in,unit_logfile)

            save_emission_subgrid_min(x_dim_index)=read_name_real('save_emission_subgrid_min(x_dim_index)',save_emission_subgrid_min(x_dim_index),unit_in,unit_logfile)
            save_emission_subgrid_min(y_dim_index)=read_name_real('save_emission_subgrid_min(y_dim_index)',save_emission_subgrid_min(y_dim_index),unit_in,unit_logfile)
            save_emission_subgrid_delta(x_dim_index)=read_name_real('save_emission_subgrid_delta(x_dim_index)',save_emission_subgrid_delta(x_dim_index),unit_in,unit_logfile)
            save_emission_subgrid_delta(y_dim_index)=read_name_real('save_emission_subgrid_delta(y_dim_index)',save_emission_subgrid_delta(y_dim_index),unit_in,unit_logfile)
            save_emission_subgrid_dim(x_dim_index)=read_name_integer('save_emission_subgrid_dim(x_dim_index)',save_emission_subgrid_dim(x_dim_index),unit_in,unit_logfile)
            save_emission_subgrid_dim(y_dim_index)=read_name_integer('save_emission_subgrid_dim(y_dim_index)',save_emission_subgrid_dim(y_dim_index),unit_in,unit_logfile)

            trace_emissions_from_in_region=read_name_logical('trace_emissions_from_in_region',trace_emissions_from_in_region,unit_in,unit_logfile)

            calc_grid_vertical_average_concentration_annual_flag=read_name_logical('calc_grid_vertical_average_concentration_annual_flag',calc_grid_vertical_average_concentration_annual_flag,unit_in,unit_logfile)

            wind_level_zc_flag=read_name_logical('wind_level_zc_flag',wind_level_zc_flag,unit_in,unit_logfile)

            use_alternative_ppm_variable_for_lf=read_name_logical('use_alternative_ppm_variable_for_lf',use_alternative_ppm_variable_for_lf,unit_in,unit_logfile)
            alternative_ppm_variable_for_lf_dim=read_name_integer('alternative_ppm_variable_for_lf_dim',alternative_ppm_variable_for_lf_dim,unit_in,unit_logfile)

            var_name_nc(conc_nc_index,pm25_nc_index,extrasource_nc_index)=read_name_char('var_name_nc(conc_nc_index,pm25_nc_index,extrasource_nc_index)',var_name_nc(conc_nc_index,pm25_nc_index,extrasource_nc_index),unit_in,unit_logfile)
            var_name_nc(conc_nc_index,pmco_nc_index,extrasource_nc_index)=read_name_char('var_name_nc(conc_nc_index,pmco_nc_index,extrasource_nc_index)',var_name_nc(conc_nc_index,pmco_nc_index,extrasource_nc_index),unit_in,unit_logfile)

            save_emep_OP_species=read_name_logical('save_emep_OP_species',save_emep_OP_species,unit_in,unit_logfile)

            ! Read configs added for the nonlocal from-in-region method
            pathname_region_mask=read_name_char('pathname_region_mask',pathname_region_mask,unit_in,unit_logfile)
            filename_region_mask=read_name_char('filename_region_mask',filename_region_mask,unit_in,unit_logfile)
            varname_region_mask=read_name_char('varname_region_mask',varname_region_mask,unit_in,unit_logfile)

            close (unit_in)

        enddo !End configuration file number loop

        !The rest below was inside the file loop before. Have moved to outside now. Hope that works!

        !Call some error traps
        if (len(trim(pathname_output_grid)).eq.0) then
            write (unit_logfile,'(A)') 'WARNING: No output path given in configuration file. Stopping'
            stop
        endif
        if (local_subgrid_method_flag == 1 .and. trace_emissions_from_in_region) then
            write(unit_logfile,'(A)') 'trace_emissions_from_in_region must be false when local_subgrid_method_flag is 1, since it is not possible to trace emission origin when EMEP concentrations are directly redistributed. Stopping'
            stop
        end if

        !Find the correct compound index based on the compound string
        do i=1,n_pollutant_nc_index
            if (trim(var_name_nc(conc_nc_index,i,allsource_nc_index)).eq.trim(input_comp_name)) then
                compound_index=i
                pollutant_index=i
                write(unit_logfile,*) 'Selected pollutant: ',trim(input_comp_name),i
            endif
        enddo


        !Replace some of the strings with the date_str. Do this twice in case there are two occurences of it in a string
        do i=1,2
            pathname_EMEP(1)=replace_string_char(config_date_str,replacement_date_str,pathname_EMEP(1))
            pathname_EMEP(2)=replace_string_char(config_date_str,replacement_date_str,pathname_EMEP(2))
            pathname_EMEP(3)=replace_string_char(config_date_str,replacement_date_str,pathname_EMEP(3))
            pathname_EMEP(4)=replace_string_char(config_date_str,replacement_date_str,pathname_EMEP(4))
            filename_EMEP(1)=replace_string_char(config_date_str,replacement_date_str,filename_EMEP(1))
            filename_EMEP(2)=replace_string_char(config_date_str,replacement_date_str,filename_EMEP(2))
            filename_EMEP(3)=replace_string_char(config_date_str,replacement_date_str,filename_EMEP(3))
            filename_EMEP(4)=replace_string_char(config_date_str,replacement_date_str,filename_EMEP(4))
            pathname_output_grid=replace_string_char(config_date_str,replacement_date_str,pathname_output_grid)
            filename_date_output_grid=replace_string_char(config_date_str,replacement_date_str,filename_date_output_grid)
            !NORTRIP file and path name
            pathname_rl(2)=replace_string_char(config_date_str,replacement_date_str,pathname_rl(2))
            filename_rl(2)=replace_string_char(config_date_str,replacement_date_str,filename_rl(2))

            pathname_EMEP(1)=replace_string_char(forecast_hour_str,replacement_hour_str,pathname_EMEP(1))
            pathname_EMEP(2)=replace_string_char(forecast_hour_str,replacement_hour_str,pathname_EMEP(2))
            pathname_EMEP(3)=replace_string_char(forecast_hour_str,replacement_hour_str,pathname_EMEP(3))
            pathname_EMEP(4)=replace_string_char(forecast_hour_str,replacement_hour_str,pathname_EMEP(4))
            filename_EMEP(1)=replace_string_char(forecast_hour_str,replacement_hour_str,filename_EMEP(1))
            filename_EMEP(2)=replace_string_char(forecast_hour_str,replacement_hour_str,filename_EMEP(2))
            filename_EMEP(3)=replace_string_char(forecast_hour_str,replacement_hour_str,filename_EMEP(3))
            filename_EMEP(4)=replace_string_char(forecast_hour_str,replacement_hour_str,filename_EMEP(4))
            pathname_output_grid=replace_string_char(forecast_hour_str,replacement_hour_str,pathname_output_grid)
            filename_date_output_grid=replace_string_char(forecast_hour_str,replacement_hour_str,filename_date_output_grid)
            !NORTRIP file and path name
            pathname_rl(2)=replace_string_char(NORTRIP_hour_str,NORTRIP_replacement_hour_str,pathname_rl(2))
            filename_rl(2)=replace_string_char(NORTRIP_hour_str,NORTRIP_replacement_hour_str,filename_rl(2))

        enddo

        !Replace date in the output file if required, 3 times for yyyy mm dd
        format_temp='yyyymmdd'
        call datestr_to_date(config_date_str,format_temp,a)
        !write (unit_logfile,'(2A)') ' Updating output path from: ',trim(pathname_output_grid)

        do i=1,3
            call date_to_datestr_bracket(a,pathname_output_grid,pathname_output_grid)
            call date_to_datestr_bracket(a,filename_date_output_grid,filename_date_output_grid)
            call date_to_datestr_bracket(a,pathname_EMEP(1),pathname_EMEP(1))
            call date_to_datestr_bracket(a,pathname_EMEP(2),pathname_EMEP(2))
            call date_to_datestr_bracket(a,pathname_EMEP(3),pathname_EMEP(3))
            call date_to_datestr_bracket(a,pathname_EMEP(4),pathname_EMEP(4))
            call date_to_datestr_bracket(a,pathname_rl(1),pathname_rl(1))
            call date_to_datestr_bracket(a,pathname_rl(2),pathname_rl(2))

            call date_to_datestr_bracket(a,filename_EMEP(1),filename_EMEP(1))
            call date_to_datestr_bracket(a,filename_EMEP(2),filename_EMEP(2))
            call date_to_datestr_bracket(a,filename_EMEP(3),filename_EMEP(3))
            call date_to_datestr_bracket(a,filename_EMEP(4),filename_EMEP(4))
            call date_to_datestr_bracket(a,filename_rl(1),filename_rl(1))
            call date_to_datestr_bracket(a,filename_rl(2),filename_rl(2))

            call date_to_datestr_bracket(a,pathname_emissions_for_EMEP,pathname_emissions_for_EMEP)

        enddo
        !write (unit_logfile,'(2A)') ' Updating output path to:   ',trim(pathname_output_grid)
        !write (unit_logfile,'(2A)') ' Updating output file to:   ',trim(pathname_EMEP(1))

        !Specify the yesterday date string and replace it if found. Identified wiht a square bracket
        !call datestr_to_date(config_date_str,format_temp,a)
        datenum_temp=date_to_number(a,ref_year_EMEP)
        datenum_temp=datenum_temp-1.
        call number_to_date(datenum_temp,a,ref_year_EMEP)
        call date_to_datestr(a,format_temp,yesterday_date_str)
        write(unit_logfile,'(2a)')'Todays date string: ',trim(config_date_str)
        write(unit_logfile,'(2a)')'Yesterdays date string: ',trim(yesterday_date_str)

        do i=1,3
            call date_to_datestr_squarebracket(a,pathname_output_grid,pathname_output_grid)
            call date_to_datestr_squarebracket(a,filename_date_output_grid,filename_date_output_grid)
            call date_to_datestr_squarebracket(a,pathname_EMEP(1),pathname_EMEP(1))
            call date_to_datestr_squarebracket(a,pathname_EMEP(2),pathname_EMEP(2))
            call date_to_datestr_squarebracket(a,pathname_EMEP(3),pathname_EMEP(3))
            call date_to_datestr_squarebracket(a,pathname_EMEP(4),pathname_EMEP(4))
            call date_to_datestr_squarebracket(a,pathname_rl(1),pathname_rl(1))
            call date_to_datestr_squarebracket(a,pathname_rl(2),pathname_rl(2))
        enddo

        do i=1,2
            pathname_EMEP(1)=replace_string_char(yesterday_date_str,replacement_yesterday_date_str,pathname_EMEP(1))
            pathname_EMEP(2)=replace_string_char(yesterday_date_str,replacement_yesterday_date_str,pathname_EMEP(2))
            pathname_EMEP(3)=replace_string_char(yesterday_date_str,replacement_yesterday_date_str,pathname_EMEP(3))
            pathname_EMEP(4)=replace_string_char(yesterday_date_str,replacement_yesterday_date_str,pathname_EMEP(4))
            filename_EMEP(1)=replace_string_char(yesterday_date_str,replacement_yesterday_date_str,filename_EMEP(1))
            filename_EMEP(2)=replace_string_char(yesterday_date_str,replacement_yesterday_date_str,filename_EMEP(2))
            filename_EMEP(3)=replace_string_char(yesterday_date_str,replacement_yesterday_date_str,filename_EMEP(3))
            filename_EMEP(4)=replace_string_char(yesterday_date_str,replacement_yesterday_date_str,filename_EMEP(4))
            pathname_output_grid=replace_string_char(yesterday_date_str,replacement_yesterday_date_str,pathname_output_grid)
            filename_date_output_grid=replace_string_char(yesterday_date_str,replacement_yesterday_date_str,filename_date_output_grid)
            !NORTRIP file and path name
            pathname_rl(2)=replace_string_char(yesterday_date_str,replacement_yesterday_date_str,pathname_rl(2))
            filename_rl(2)=replace_string_char(yesterday_date_str,replacement_yesterday_date_str,filename_rl(2))
        enddo

        !write(*,*) trim(filename_EMEP(3))
        !write(*,*) trim(pathname_EMEP(3))


        !Place tile_tag in front of file_tag if it has been read
        if (tile_tag.ne.'') then
            file_tag=trim(file_tag)//'_'//trim(tile_tag)
        endif

    end subroutine uEMEP_read_config


!----------------------------------------------------------------------
    function replace_string_char(replace_str,match_str,read_str)
        !Finds a match_str in read_str and replaces it with replace_str to give a new version of read_str
        implicit none

        character(256) replace_string_char
        character (*) match_str,replace_str,read_str
        character(256) temp_str1,temp_str2
        integer index_start,index_stop

        replace_string_char=read_str

        index_start=index(read_str,trim(match_str))
        if (index_start.ne.0) then
            index_stop=index_start+len(trim(match_str))
            temp_str1=read_str(1:index_start-1)
            temp_str2=read_str(index_stop:len(read_str))
            replace_string_char=trim(temp_str1)//trim(replace_str)//trim(temp_str2)
        endif
        !write(*,'(A)') trim(replace_string_char)

    end function replace_string_char
!----------------------------------------------------------------------

end module read_config
