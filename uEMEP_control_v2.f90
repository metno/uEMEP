!  uEMEP.f90 
!
!
!****************************************************************************
!
!  PROGRAM: uEMEP
!
!   Bruce rolstad Denby (brucerd@met.no)
!****************************************************************************

    program uEMEP_v2

    use uEMEP_definitions
   
    implicit none
    
    integer source_index

    
    write(*,*) ''
    write(*,*) '------------------------------------------------------------------------'
    write(*,*) 'Starting programm uEMEP_v2.1'
    write(*,*) '------------------------------------------------------------------------'
    
    !Assign the configuration file name and substitution date_str from the command line
    call get_command_argument(1,name_config_file)
    call get_command_argument(2,config_date_str)
    write(*,*) 'name_config_file= ',trim(name_config_file)
    write(*,*) 'config_date_str= ',trim(config_date_str)
    
    !Initialise parameters and read name file
    !call uEMEP_setup

    call uEMEP_set_constants
    call uEMEP_read_config
    call uEMEP_set_filenames
    call uEMEP_set_subgrids    
    call uEMEP_set_emission_factors
        
    if (use_single_time_loop_flag) then
        start_time_loop_index=1
        end_time_loop_index=end_time_nc_index-start_time_nc_index+1
        subgrid_dim(t_dim_index)=1     
        dim_length_nc(time_dim_nc_index)=1
    else
        start_time_loop_index=1
        end_time_loop_index=1
        subgrid_dim(t_dim_index)=end_time_nc_index-start_time_nc_index+1
        dim_length_nc(time_dim_nc_index)=subgrid_dim(t_dim_index)
    endif

    
    do t_loop=start_time_loop_index,end_time_loop_index
   
        write(*,*) 'TIME=',t_loop,' OF ',end_time_loop_index
        if (unit_logfile.ne.0) then 
            write(unit_logfile,*) 'TIME=',t_loop,' OF ',end_time_loop_index
        endif
        
        !Read EMEP data and grid
        call uEMEP_read_EMEP
 
        if (t_loop.eq.start_time_loop_index) then
        
            !Define subgrid positions and buffer zones
            call uEMEP_define_subgrid
        
            !Define and allocate cross reference subgrids
            call uEMEP_crossreference_grids         
        
            if (calculate_source(traffic_index)) then
                !Read road link data from ascii files
                call uEMEP_read_roadlink_data_ascii
                !Sub-grid traffic data    
                call uEMEP_grid_roads
            endif
 
            if (calculate_source(shipping_index)) then
                !Read and subgrid shipping data
                call uEMEP_preaggregate_shipping_asi_data
                call uEMEP_read_shipping_asi_data
            endif

            if (calculate_source(heating_index)) then
                !Read and subgrid SSB dwelling data
                SSB_data_type=dwelling_index
                call uEMEP_read_SSB_data
            endif

            if (calculate_source(agriculture_index)) then
                !Read and subgrid agriculture data
                call uEMEP_read_agriculture_rivm_data
            endif

            if (calculate_population_exposure_flag) then
                !Read and subgrid SSB population data
                SSB_data_type=population_data_type
                call uEMEP_read_SSB_data
            endif

            !Autogrid setting
            call uEMEP_auto_subgrid
    
            !Read receptor data and distribute to subgrids
            call uEMEP_read_receptor_data

        endif
    
        !Interpolate meteo data to subgrid
        call uEMEP_subgrid_meteo_EMEP
        
        !Convert proxies to emissions
        call uEMEP_convert_proxy_to_emissions
    
        !Replaces proxy emissions with distributed EMEP emissions
        call uEMEP_subgrid_emission_EMEP
        
        !Grid pseudo dispersion
        traveltime_subgrid=0.
        do source_index=1,n_source_index
        if (calculate_source(source_index)) then
            call uEMEP_grid_proxy(source_index)
        endif
        enddo
    
        !Diagnostic for comparing EMEP and proxy data emissions
        call uEMEP_aggregate_proxy_emission_in_EMEP_grid
    
        !Put EMEP data into subgrids for all sources.
        call uEMEP_subgrid_EMEP
        
        !Interpolate EMEP to sub-grid
        do source_index=1,n_source_index
        if (calculate_source(source_index)) then
            call uEMEP_redistribute_local_source(source_index)
            !Do not call this routine. It scales the proxy emission dispersion with 1/FF. THis is now done in uEMEP_grid_proxy
            call uEMEP_disperse_local_source(source_index)
        endif
        enddo
    
        !Combine and save sources in local and total values
        call uEMEP_combine_local_source
    
        if (compound_index.eq.nox_index) then
            call uEMEP_chemistry
        endif

        if (calculate_population_exposure_flag) then
            call uEMEP_calculate_exposure
        endif
    
        call uEMEP_save_netcdf_control
    
    enddo !t_loop
    
    if (unit_logfile.gt.0) then
         close(unit_logfile,status='keep')
    endif

    write(*,*) ''
    write(*,*) '------------------------------------------------------------------------'
    write(*,*) 'Ending programm uEMEP'
    write(*,*) '------------------------------------------------------------------------'

    end program uEMEP_v2

