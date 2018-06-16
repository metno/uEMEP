!   Set up of uEMEP before starting calculations
    subroutine uEMEP_read_config
    
    use uEMEP_definitions

    implicit none
    
    integer i,j,k
    real read_name_real
    logical read_name_logical
    integer read_name_integer
    character(256) read_name_char,pathfilename_log_file
    integer exists
    
    integer :: unit_in=30
    integer i_config
    
    !Functions
    character(256) replace_string_char
    
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
        if (len(trim(filename_log_file)).gt.0.and.i_config.eq.1) then
            unit_logfile=10 
            !Check existence of path
            inquire(directory=trim(pathname_log_file),exist=exists)
            if (.not.exists) then
                write(unit_logfile,'(A)')'ERROR: Log file directory path '//trim(pathname_log_file)//' does not exist.'
                stop
            endif
            !Write to screen if writing to log file
            pathfilename_log_file=trim(pathname_log_file)//trim(filename_log_file)
            write(*,'(A,A)') 'Writing to log file: ', trim(pathfilename_log_file)
            write(*,'(A)') '================================================================'
            open(unit_logfile,file=trim(pathfilename_log_file),access='sequential',form='formatted',status='unknown')
            write(unit_logfile,'(A)') ''
            write(unit_logfile,'(A)') '================================================================'
            write(unit_logfile,'(A)') 'Starting program uEMEP v2.1' 
  	        write(unit_logfile,'(A)') '================================================================'
        else
            unit_logfile=0
        endif
              
        file_tag=read_name_char('file_tag',file_tag,unit_in,unit_logfile)
        
        replacement_date_str=read_name_char('replacement_date_str',replacement_date_str,unit_in,unit_logfile)
        
        input_comp_name=read_name_char('input_comp_name',input_comp_name,unit_in,unit_logfile)
        
        hourly_calculations=read_name_logical('hourly_calculations',hourly_calculations,unit_in,unit_logfile)
        annual_calculations=read_name_logical('annual_calculations',annual_calculations,unit_in,unit_logfile)
        
        start_time_nc_index=read_name_integer('start_time_nc_index',start_time_nc_index,unit_in,unit_logfile)
        end_time_nc_index=read_name_integer('end_time_nc_index',end_time_nc_index,unit_in,unit_logfile)
        start_time_meteo_nc_index=read_name_integer('start_time_meteo_nc_index',start_time_meteo_nc_index,unit_in,unit_logfile)
        end_time_meteo_nc_index=read_name_integer('end_time_meteo_nc_index',end_time_meteo_nc_index,unit_in,unit_logfile)

        use_single_time_loop_flag=read_name_logical('use_single_time_loop_flag',use_single_time_loop_flag,unit_in,unit_logfile)
        reduce_EMEP_region_flag=read_name_logical('reduce_EMEP_region_flag',reduce_EMEP_region_flag,unit_in,unit_logfile)
        save_intermediate_files=read_name_logical('save_intermediate_files',save_intermediate_files,unit_in,unit_logfile)
        use_multiple_receptor_grids_flag=read_name_logical('use_multiple_receptor_grids_flag',use_multiple_receptor_grids_flag,unit_in,unit_logfile)
        use_receptor_region=read_name_integer('use_receptor_region',use_receptor_region,unit_in,unit_logfile)

        
        !Read in choice of reading existing proxy emission, proxy dispersion, meteorology and use_subgid file data. These will rarely, if ever, be used
        read_existing_grid_data(proxy_emission_file_index(:))=read_name_logical('read_existing_grid_data(proxy_emission_file_index(:))',read_existing_grid_data(proxy_emission_file_index(allsource_index)),unit_in,unit_logfile)
        read_existing_grid_data(proxy_emission_file_index(traffic_index))=read_name_logical('read_existing_grid_data(proxy_emission_file_index(traffic_index))',read_existing_grid_data(proxy_emission_file_index(traffic_index)),unit_in,unit_logfile)
        read_existing_grid_data(proxy_emission_file_index(shipping_index))=read_name_logical('read_existing_grid_data(proxy_emission_file_index(shipping_index))',read_existing_grid_data(proxy_emission_file_index(shipping_index)),unit_in,unit_logfile)
        read_existing_grid_data(proxy_emission_file_index(heating_index))=read_name_logical('read_existing_grid_data(proxy_emission_file_index(heating_index))',read_existing_grid_data(proxy_emission_file_index(heating_index)),unit_in,unit_logfile)
        read_existing_grid_data(proxy_emission_file_index(agriculture_index))=read_name_logical('read_existing_grid_data(proxy_emission_file_index(agriculture_index))',read_existing_grid_data(proxy_emission_file_index(agriculture_index)),unit_in,unit_logfile)
        read_existing_grid_data(proxy_emission_file_index(industry_index))=read_name_logical('read_existing_grid_data(proxy_emission_file_index(industry_index))',read_existing_grid_data(proxy_emission_file_index(industry_index)),unit_in,unit_logfile)
 
        read_existing_grid_data(proxy_file_index(:))=read_name_logical('read_existing_grid_data(proxy_file_index(:))',read_existing_grid_data(proxy_file_index(allsource_index)),unit_in,unit_logfile)
        read_existing_grid_data(proxy_file_index(traffic_index))=read_name_logical('read_existing_grid_data(proxy_file_index(traffic_index))',read_existing_grid_data(proxy_file_index(traffic_index)),unit_in,unit_logfile)
        read_existing_grid_data(proxy_file_index(shipping_index))=read_name_logical('read_existing_grid_data(proxy_file_index(shipping_index))',read_existing_grid_data(proxy_file_index(shipping_index)),unit_in,unit_logfile)
        read_existing_grid_data(proxy_file_index(heating_index))=read_name_logical('read_existing_grid_data(proxy_file_index(heating_index))',read_existing_grid_data(proxy_file_index(heating_index)),unit_in,unit_logfile)
        read_existing_grid_data(proxy_file_index(agriculture_index))=read_name_logical('read_existing_grid_data(proxy_file_index(agriculture_index))',read_existing_grid_data(proxy_file_index(agriculture_index)),unit_in,unit_logfile)
        read_existing_grid_data(proxy_file_index(industry_index))=read_name_logical('read_existing_grid_data(proxy_file_index(industry_index))',read_existing_grid_data(proxy_file_index(industry_index)),unit_in,unit_logfile)

        read_existing_grid_data(emep_subgrid_file_index(:))=read_name_logical('read_existing_grid_data(emep_subgrid_file_index(:))',read_existing_grid_data(emep_subgrid_file_index(allsource_index)),unit_in,unit_logfile)
        read_existing_grid_data(emep_subgrid_file_index(traffic_index))=read_name_logical('read_existing_grid_data(emep_subgrid_file_index(traffic_index))',read_existing_grid_data(emep_subgrid_file_index(traffic_index)),unit_in,unit_logfile)
        read_existing_grid_data(emep_subgrid_file_index(shipping_index))=read_name_logical('read_existing_grid_data(emep_subgrid_file_index(shipping_index))',read_existing_grid_data(emep_subgrid_file_index(shipping_index)),unit_in,unit_logfile)
        read_existing_grid_data(emep_subgrid_file_index(heating_index))=read_name_logical('read_existing_grid_data(emep_subgrid_file_index(heating_index))',read_existing_grid_data(emep_subgrid_file_index(heating_index)),unit_in,unit_logfile)
        read_existing_grid_data(emep_subgrid_file_index(agriculture_index))=read_name_logical('read_existing_grid_data(emep_subgrid_file_index(agriculture_index))',read_existing_grid_data(emep_subgrid_file_index(agriculture_index)),unit_in,unit_logfile)
        read_existing_grid_data(emep_subgrid_file_index(industry_index))=read_name_logical('read_existing_grid_data(emep_subgrid_file_index(industry_index))',read_existing_grid_data(emep_subgrid_file_index(industry_index)),unit_in,unit_logfile)
        
        read_existing_grid_data(subgrid_meteo_file_index)=read_name_logical('read_existing_grid_data(subgrid_meteo_file_index)',read_existing_grid_data(subgrid_meteo_file_index),unit_in,unit_logfile)

        read_existing_grid_data(use_subgrid_file_index(:))=read_name_logical('read_existing_grid_data(use_subgrid_file_index(:))',read_existing_grid_data(use_subgrid_file_index(allsource_index)),unit_in,unit_logfile)
        read_existing_grid_data(use_subgrid_file_index(traffic_index))=read_name_logical('read_existing_grid_data(use_subgrid_file_index(traffic_index))',read_existing_grid_data(use_subgrid_file_index(traffic_index)),unit_in,unit_logfile)
        read_existing_grid_data(use_subgrid_file_index(shipping_index))=read_name_logical('read_existing_grid_data(use_subgrid_file_index(shipping_index))',read_existing_grid_data(use_subgrid_file_index(shipping_index)),unit_in,unit_logfile)
        read_existing_grid_data(use_subgrid_file_index(heating_index))=read_name_logical('read_existing_grid_data(use_subgrid_file_index(heating_index))',read_existing_grid_data(use_subgrid_file_index(heating_index)),unit_in,unit_logfile)
        read_existing_grid_data(use_subgrid_file_index(agriculture_index))=read_name_logical('read_existing_grid_data(use_subgrid_file_index(agriculture_index))',read_existing_grid_data(use_subgrid_file_index(agriculture_index)),unit_in,unit_logfile)
        read_existing_grid_data(use_subgrid_file_index(industry_index))=read_name_logical('read_existing_grid_data(use_subgrid_file_index(industry_index))',read_existing_grid_data(use_subgrid_file_index(industry_index)),unit_in,unit_logfile)

        !Choose which sources to calculate
        !calculate_source(:)=read_name_logical('calculate_source(:)',calculate_source(allsource_index),unit_in,unit_logfile)
        calculate_source(traffic_index)=read_name_logical('calculate_source(traffic_index)',calculate_source(traffic_index),unit_in,unit_logfile)
        calculate_source(shipping_index)=read_name_logical('calculate_source(shipping_index)',calculate_source(shipping_index),unit_in,unit_logfile)
        calculate_source(heating_index)=read_name_logical('calculate_source(heating_index)',calculate_source(heating_index),unit_in,unit_logfile)
        calculate_source(agriculture_index)=read_name_logical('calculate_source(agriculture_index)',calculate_source(agriculture_index),unit_in,unit_logfile)
        calculate_source(industry_index)=read_name_logical('calculate_source(industry_index)',calculate_source(industry_index),unit_in,unit_logfile)

        !For aggregating proxy emission data to EMEP grids
        !make_EMEP_grid_emission_data(:)=read_name_logical('make_EMEP_grid_emission_data(:)',make_EMEP_grid_emission_data(allsource_index),unit_in,unit_logfile)
        make_EMEP_grid_emission_data(traffic_index)=read_name_logical('make_EMEP_grid_emission_data(traffic_index)',make_EMEP_grid_emission_data(traffic_index),unit_in,unit_logfile)
        make_EMEP_grid_emission_data(shipping_index)=read_name_logical('make_EMEP_grid_emission_data(shipping_index)',make_EMEP_grid_emission_data(shipping_index),unit_in,unit_logfile)
        make_EMEP_grid_emission_data(heating_index)=read_name_logical('make_EMEP_grid_emission_data(heating_index)',make_EMEP_grid_emission_data(heating_index),unit_in,unit_logfile)
        make_EMEP_grid_emission_data(agriculture_index)=read_name_logical('make_EMEP_grid_emission_data(agriculture_index)',make_EMEP_grid_emission_data(agriculture_index),unit_in,unit_logfile)
        make_EMEP_grid_emission_data(industry_index)=read_name_logical('make_EMEP_grid_emission_data(industry_index)',make_EMEP_grid_emission_data(industry_index),unit_in,unit_logfile)
        
        replace_EMEP_local_with_subgrid_local(:)=read_name_logical('replace_EMEP_local_with_subgrid_local',replace_EMEP_local_with_subgrid_local(allsource_index),unit_in,unit_logfile)
        replace_EMEP_local_with_subgrid_local(traffic_index)=read_name_logical('replace_EMEP_local_with_subgrid_local(traffic_index)',replace_EMEP_local_with_subgrid_local(traffic_index),unit_in,unit_logfile)
        replace_EMEP_local_with_subgrid_local(shipping_index)=read_name_logical('replace_EMEP_local_with_subgrid_local(shipping_index)',replace_EMEP_local_with_subgrid_local(shipping_index),unit_in,unit_logfile)
        replace_EMEP_local_with_subgrid_local(heating_index)=read_name_logical('replace_EMEP_local_with_subgrid_local(heating_index)',replace_EMEP_local_with_subgrid_local(heating_index),unit_in,unit_logfile)
        replace_EMEP_local_with_subgrid_local(agriculture_index)=read_name_logical('replace_EMEP_local_with_subgrid_local(agriculture_index)',replace_EMEP_local_with_subgrid_local(agriculture_index),unit_in,unit_logfile)
        replace_EMEP_local_with_subgrid_local(industry_index)=read_name_logical('replace_EMEP_local_with_subgrid_local(industry_index)',replace_EMEP_local_with_subgrid_local(industry_index),unit_in,unit_logfile)

        projection_type=read_name_integer('projection_type',projection_type,unit_in,unit_logfile)
        utm_zone=read_name_integer('utm_zone',utm_zone,unit_in,unit_logfile)
        !Present UTM central lon position if not overridden by input
        utm_lon0=abs(utm_zone)*6-180-3
        utm_lon0=read_name_real('utm_lon0',utm_lon0,unit_in,unit_logfile)

        EMEP_grid_interpolation_flag=read_name_integer('EMEP_grid_interpolation_flag',EMEP_grid_interpolation_flag,unit_in,unit_logfile)
        EMEP_meteo_grid_interpolation_flag=read_name_integer('EMEP_meteo_grid_interpolation_flag',EMEP_meteo_grid_interpolation_flag,unit_in,unit_logfile)
        EMEP_emission_grid_interpolation_flag=read_name_integer('EMEP_emission_grid_interpolation_flag',EMEP_emission_grid_interpolation_flag,unit_in,unit_logfile)
        subgrid_emission_distribution_flag=read_name_logical('subgrid_emission_distribution_flag',subgrid_emission_distribution_flag,unit_in,unit_logfile)
   
        EMEP_grid_interpolation_size=read_name_real('EMEP_grid_interpolation_size',EMEP_grid_interpolation_size,unit_in,unit_logfile)
        use_downwind_position_flag=read_name_logical('use_downwind_position_flag',use_downwind_position_flag,unit_in,unit_logfile)
        local_subgrid_method_flag=read_name_integer('local_subgrid_method_flag',local_subgrid_method_flag,unit_in,unit_logfile)
        
        stability_scheme_flag=read_name_integer('stability_scheme_flag',stability_scheme_flag,unit_in,unit_logfile)
        average_zc_h_in_Kz_flag=read_name_logical('average_zc_h_in_Kz_flag',average_zc_h_in_Kz_flag,unit_in,unit_logfile)
        
        wind_level_flag=read_name_integer('wind_level_flag',wind_level_flag,unit_in,unit_logfile)
        wind_level_integral_flag=read_name_integer('wind_level_integral_flag',wind_level_flag,unit_in,unit_logfile) !Default is wind_level_flag
        no2_chemistry_scheme_flag=read_name_integer('no2_chemistry_scheme_flag',no2_chemistry_scheme_flag,unit_in,unit_logfile)
        
        use_emission_positions_for_auto_subgrid_flag=read_name_logical('use_emission_positions_for_auto_subgrid_flag',use_emission_positions_for_auto_subgrid_flag,unit_in,unit_logfile)
        use_receptor_positions_for_auto_subgrid_flag=read_name_logical('use_receptor_positions_for_auto_subgrid_flag',use_receptor_positions_for_auto_subgrid_flag,unit_in,unit_logfile)
        use_population_positions_for_auto_subgrid_flag=read_name_logical('use_population_positions_for_auto_subgrid_flag',use_population_positions_for_auto_subgrid_flag,unit_in,unit_logfile)
        interpolate_subgrids_flag=read_name_logical('interpolate_subgrids_flag',interpolate_subgrids_flag,unit_in,unit_logfile)

        use_trajectory_flag=read_name_logical('use_trajectory_flag',use_trajectory_flag,unit_in,unit_logfile)
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
        
        !Specifies the number of subsources for each source. Usually not used as default is 1
        n_subsource(:)=read_name_integer('n_subsource(:)',n_subsource(allsource_index),unit_in,unit_logfile)
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
        sig_y_00(traffic_index,2)=read_name_real('sig_y_00(traffic_index,2)',sig_y_00(traffic_index,2),unit_in,unit_logfile)
        sig_y_00(shipping_index,2)=read_name_real('sig_y_00(shipping_index,2)',sig_y_00(shipping_index,2),unit_in,unit_logfile)
        sig_y_00(heating_index,2)=read_name_real('sig_y_00(heating_index,2)',sig_y_00(heating_index,2),unit_in,unit_logfile)
        sig_y_00(agriculture_index,2)=read_name_real('sig_y_00(agriculture_index,2)',sig_y_00(agriculture_index,2),unit_in,unit_logfile)
        sig_y_00(industry_index,2)=read_name_real('sig_y_00(industry_index,2)',sig_y_00(industry_index,2),unit_in,unit_logfile)

        sig_z_00(traffic_index,1)=read_name_real('sig_z_00(traffic_index,1)',sig_z_00(traffic_index,1),unit_in,unit_logfile)
        sig_z_00(shipping_index,1)=read_name_real('sig_z_00(shipping_index,1)',sig_z_00(shipping_index,1),unit_in,unit_logfile)
        sig_z_00(heating_index,1)=read_name_real('sig_z_00(heating_index,1)',sig_z_00(heating_index,1),unit_in,unit_logfile)
        sig_z_00(agriculture_index,1)=read_name_real('sig_z_00(agriculture_index,1)',sig_z_00(agriculture_index,1),unit_in,unit_logfile)
        sig_z_00(industry_index,1)=read_name_real('sig_z_00(industry_index,1)',sig_z_00(industry_index,1),unit_in,unit_logfile)
        sig_z_00(traffic_index,2)=read_name_real('sig_z_00(traffic_index,2)',sig_z_00(traffic_index,2),unit_in,unit_logfile)
        sig_z_00(shipping_index,2)=read_name_real('sig_z_00(shipping_index,2)',sig_z_00(shipping_index,2),unit_in,unit_logfile)
        sig_z_00(heating_index,2)=read_name_real('sig_z_00(heating_index,2)',sig_z_00(heating_index,2),unit_in,unit_logfile)
        sig_z_00(agriculture_index,2)=read_name_real('sig_z_00(agriculture_index,2)',sig_z_00(agriculture_index,2),unit_in,unit_logfile)
        sig_z_00(industry_index,2)=read_name_real('sig_z_00(industry_index,2)',sig_z_00(industry_index,2),unit_in,unit_logfile)
        
        !Read output grid path for all data
        pathname_output_grid=read_name_char('pathname_output_grid',pathname_output_grid,unit_in,unit_logfile)
        
        !Read in file names. Only 2 choices for most file types
        pathname_rl(1)=read_name_char('pathname_rl(1)',pathname_rl(1),unit_in,unit_logfile)
        pathname_rl(2)=read_name_char('pathname_rl(2)',pathname_rl(2),unit_in,unit_logfile)
        filename_rl(1)=read_name_char('filename_rl(1)',filename_rl(1),unit_in,unit_logfile)
        filename_rl(2)=read_name_char('filename_rl(2)',filename_rl(2),unit_in,unit_logfile)

        pathname_EMEP(1)=read_name_char('pathname_EMEP(1)',pathname_EMEP(1),unit_in,unit_logfile)
        pathname_EMEP(2)=read_name_char('pathname_EMEP(2)',pathname_EMEP(2),unit_in,unit_logfile)
        pathname_EMEP(3)=read_name_char('pathname_EMEP(3)',pathname_EMEP(3),unit_in,unit_logfile)
        pathname_EMEP(4)=read_name_char('pathname_EMEP(4)',pathname_EMEP(4),unit_in,unit_logfile)
        filename_EMEP(1)=read_name_char('filename_EMEP(1)',filename_EMEP(1),unit_in,unit_logfile)
        filename_EMEP(2)=read_name_char('filename_EMEP(2)',filename_EMEP(2),unit_in,unit_logfile)
        filename_EMEP(3)=read_name_char('filename_EMEP(3)',filename_EMEP(3),unit_in,unit_logfile)
        filename_EMEP(4)=read_name_char('filename_EMEP(4)',filename_EMEP(4),unit_in,unit_logfile)

        pathname_ship(1)=read_name_char('pathname_ship(1)',pathname_ship(1),unit_in,unit_logfile)
        pathname_ship(2)=read_name_char('pathname_ship(2)',pathname_ship(2),unit_in,unit_logfile)
        filename_ship(1)=read_name_char('filename_ship(1)',filename_ship(1),unit_in,unit_logfile)
        filename_ship(2)=read_name_char('filename_ship(2)',filename_ship(2),unit_in,unit_logfile)

        pathname_agriculture(1)=read_name_char('pathname_agriculture(1)',pathname_agriculture(1),unit_in,unit_logfile)
        pathname_agriculture(2)=read_name_char('pathname_agriculture(2)',pathname_agriculture(2),unit_in,unit_logfile)
        filename_agriculture(1)=read_name_char('filename_agriculture(1)',filename_agriculture(1),unit_in,unit_logfile)
        filename_agriculture(2)=read_name_char('filename_agriculture(2)',filename_agriculture(2),unit_in,unit_logfile)
        
        pathname_heating(dwelling_index)=read_name_char('pathname_heating(dwelling_index)',pathname_heating(dwelling_index),unit_in,unit_logfile)
        pathname_heating(population_index)=read_name_char('pathname_heating(population_index)',pathname_heating(population_index),unit_in,unit_logfile)
        filename_heating(dwelling_index)=read_name_char('filename_heating(dwelling_index)',filename_heating(dwelling_index),unit_in,unit_logfile)
        filename_heating(population_index)=read_name_char('filename_heating(population_index)',filename_heating(population_index),unit_in,unit_logfile)

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

        pathname_receptor=read_name_char('pathname_receptor',pathname_receptor,unit_in,unit_logfile)
        filename_receptor=read_name_char('filename_receptor',filename_receptor,unit_in,unit_logfile)

        pathname_timeprofile=read_name_char('pathname_timeprofile',pathname_timeprofile,unit_in,unit_logfile)
        filename_timeprofile=read_name_char('filename_timeprofile',filename_timeprofile,unit_in,unit_logfile)

        population_data_type=read_name_integer('population_data_type',population_data_type,unit_in,unit_logfile)
        
        combine_emission_subsources_during_dispersion(:)=read_name_logical('combine_emission_subsources_during_dispersion(:)',combine_emission_subsources_during_dispersion(allsource_index),unit_in,unit_logfile)
        combine_emission_subsources_during_dispersion(traffic_index)=read_name_logical('combine_emission_subsources_during_dispersion(traffic_index)',combine_emission_subsources_during_dispersion(traffic_index),unit_in,unit_logfile)
        combine_emission_subsources_during_dispersion(shipping_index)=read_name_logical('combine_emission_subsources_during_dispersion(shipping_index)',combine_emission_subsources_during_dispersion(shipping_index),unit_in,unit_logfile)
        combine_emission_subsources_during_dispersion(heating_index)=read_name_logical('combine_emission_subsources_during_dispersion(heating_index)',combine_emission_subsources_during_dispersion(heating_index),unit_in,unit_logfile)
        combine_emission_subsources_during_dispersion(agriculture_index)=read_name_logical('combine_emission_subsources_during_dispersion(agriculture_index)',combine_emission_subsources_during_dispersion(agriculture_index),unit_in,unit_logfile)
        combine_emission_subsources_during_dispersion(industry_index)=read_name_logical('combine_emission_subsources_during_dispersion(industry_index)',combine_emission_subsources_during_dispersion(industry_index),unit_in,unit_logfile)
        
        FF_min_dispersion=read_name_real('FF_min_dispersion',FF_min_dispersion,unit_in,unit_logfile)
        emission_timeprofile_hour_shift=read_name_real('emission_timeprofile_hour_shift',emission_timeprofile_hour_shift,unit_in,unit_logfile)
        
        use_last_meteo_in_dispersion=read_name_logical('use_last_meteo_in_dispersion',use_last_meteo_in_dispersion,unit_in,unit_logfile)
        use_meandering_in_dispersion=read_name_logical('use_meandering_in_dispersion',use_meandering_in_dispersion,unit_in,unit_logfile)
        
        use_traffic_for_sigma0_flag=read_name_logical('use_traffic_for_sigma0_flag',use_traffic_for_sigma0_flag,unit_in,unit_logfile)
!        use_traffic_for_minFF_flag=read_name_logical('use_traffic_for_minFF_flag',use_traffic_for_minFF_flag,unit_in,unit_logfile)
        use_emission_grid_gradient_flag=read_name_logical('use_emission_grid_gradient_flag',use_emission_grid_gradient_flag,unit_in,unit_logfile)

        use_alternative_meteorology_flag=read_name_logical('use_alternative_meteorology_flag',use_alternative_meteorology_flag,unit_in,unit_logfile)
        ustar_min=read_name_real('ustar_min',ustar_min,unit_in,unit_logfile)
        hmix_min=read_name_real('hmix_min',hmix_min,unit_in,unit_logfile)
        hmix_max=read_name_real('hmix_max',hmix_max,unit_in,unit_logfile)
        !use_alternative_z0_flag=read_name_logical('use_alternative_z0_flag',use_alternative_z0_flag,unit_in,unit_logfile)
       
        
        !Read emission factors
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
        
        save_netcdf_file_flag=read_name_logical('save_netcdf_file_flag',save_netcdf_file_flag,unit_in,unit_logfile)
        save_netcdf_receptor_flag=read_name_logical('save_netcdf_receptor_flag',save_netcdf_receptor_flag,unit_in,unit_logfile)
        
        calculate_tiling_flag=read_name_logical('calculate_tiling_flag',calculate_tiling_flag,unit_in,unit_logfile)        
        pathname_tiles=read_name_char('pathname_tiles','',unit_in,unit_logfile)
        filename_tiles=read_name_char('filename_tiles','',unit_in,unit_logfile)
        tile_tag=read_name_char('tile_tag','',unit_in,unit_logfile)

    close (unit_in)
    
    !Call some error traps
    if (len(trim(pathname_output_grid)).eq.0) then
        write (unit_logfile,'(A)') 'WARNING: No output path given in configuration file. Stopping'
    endif

    !Find the correct compound index based on the compound string
    do i=1,n_compound_nc_index
        if (trim(var_name_nc(conc_nc_index,i,allsource_index)).eq.trim(input_comp_name)) then
            compound_index=i
        endif
        if (trim(var_name_nc(frac_nc_index,i,allsource_index)).eq.trim(input_comp_name)) then
            compound_frac_index=i
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
    enddo

    !Place tile_tag in front of file_tag if it has been read
    if (tile_tag.ne.'') then
        file_tag=trim(tile_tag)//'_'//trim(file_tag)
    endif
    
    enddo !End configuration file number loop
    

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
