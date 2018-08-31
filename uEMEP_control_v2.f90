!****************************************************************************
!   uEMEP control v2
!
!   Bruce rolstad Denby (brucerd@met.no)
!   MET Norway
!
!   Control programme for running the downscaling routine uEMEP
!****************************************************************************

    program uEMEP_v3

    use uEMEP_definitions
   
    implicit none
    
    integer source_index
    logical first_g_loop
    real start_time_cpu,end_time_cpu
    !real temp_val,area_weighted_interpolation_function
    
    call CPU_TIME(start_time_cpu)
    
    write(*,*) ''
    write(*,*) '------------------------------------------------------------------------'
    write(*,*) 'Starting programm uEMEP_v3.1'
    write(*,*) '------------------------------------------------------------------------'
    
    !Read the command line, assigning the configuration file names and the substitution date_str
    call uEMEP_read_command_line
    
    !Set constants and variable names to be read from EMEP and meteo files
    call uEMEP_set_constants
    
    !Read the configuration files. Hard coded to be up to 5 files. Log file opened in this routine
    call uEMEP_read_config
    
    !Set the pollutant and compound loop definitions
    call uEMEP_set_pollutant_loop
    
    !Set the names of files to be written to when saving intermediate files
    call uEMEP_set_filenames
    
    !Read positions of receptor points (usually observations) for specifying multiple receptor grids or calculation points within a single grid
    call uEMEP_read_receptor_data
    
    !Set the grid loop (g_loop) extent based on use_multiple_receptor_grids_flag or not
    if (use_multiple_receptor_grids_flag) then
        start_grid_loop_index=1
        end_grid_loop_index=n_receptor_in
        n_receptor=1
        n_valid_receptor=1
        valid_receptor_index(1)=1
        reduce_roadlink_region_flag=.false. !Multiple receptor flags reads in all road links and allocates to the different receptor grids. Only reads once
    else
        start_grid_loop_index=1
        end_grid_loop_index=1
        n_receptor=n_receptor_in
        use_receptor(start_grid_loop_index)=.true.
    endif
    
    first_g_loop=.true.
    
    !Start internal grid receptor loop using only those receptor grids specified in uEMEP_read_receptor_data
    do g_loop=start_grid_loop_index,end_grid_loop_index
    if (use_receptor(g_loop)) then
    
        
        !Set the grid definitions according to the receptor/observation positions
        call uEMEP_set_loop_receptor_grid
    
        !Create the subgrid
        call uEMEP_set_subgrids
        
        !Set emission factors for the current subgrid
        call uEMEP_set_emission_factors
        
        !Set the internal time loop (t_loop). If use_single_time_loop_flag=T then time array dimmensions are set to 1 and each time set of data
        !is read individually. This mostly to save memory
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
    
        !Start the internal time loop
        do t_loop=start_time_loop_index,end_time_loop_index
   
            !Write progress in time and receptor grid loop to screen
            write(*,*) 'REC LOOP= ',g_loop,' OF ',end_grid_loop_index
            if (unit_logfile.ne.0) then 
                write(unit_logfile,*) 'REC LOOP= ',g_loop,' OF ',end_grid_loop_index
            endif
            write(*,*) 'TIME LOOP=',t_loop,' OF ',end_time_loop_index
            if (unit_logfile.ne.0) then 
                write(unit_logfile,*) 'TIME LOOP=',t_loop,' OF ',end_time_loop_index
            endif
        
            
            !For the first time loop set the initial subgrid range values used in reading EMEP and meteo data
            if (t_loop.ge.start_time_loop_index) then
                init_subgrid_min=subgrid_min
                init_subgrid_max=subgrid_max
            endif
                
            !Read EMEP data and meteo grid from netcdf files
            call uEMEP_read_EMEP
            if (use_alternative_meteorology_flag.or.use_alternative_z0_flag) call uEMEP_read_meteo_nc
        
            !Set the following for the first internal time step only
            if (t_loop.eq.start_time_loop_index) then
        
                !Define subgrid positions and buffer zones. Must be done after reading EMEP data as is based on EMEP grid sizes
                call uEMEP_define_subgrid
        
                !Define and allocate cross reference subgrids used to transfer data between different subgrids
                call uEMEP_crossreference_grids
        
                !Read all road link data from ascii files
                if (calculate_source(traffic_index)) then
                    !Do this only for the first receptor grid loop
                    if (first_g_loop) then
                        call uEMEP_read_roadlink_data_ascii
                        call uEMEP_change_road_data
                        !Read in the emission data for traffic in the first g_loop if required
                        if (use_NORTRIP_emission_data) then
                            call uEMEP_read_roadlink_emission_data
                        endif       
                    endif
                    !Sub-grid traffic data for every receptor grid loop   
                    !call uEMEP_grid_roads
                endif
 
                !Read and subgrid shipping data
                if (calculate_source(shipping_index)) then
                    !If necessary aggregate shipping ddata first
                    call uEMEP_preaggregate_shipping_asi_data
                    !Read in shipping data
                    call uEMEP_read_shipping_asi_data
                endif

                !Read in proxy data fro home heating. Currently dwelling density
                if (calculate_source(heating_index)) then
                    !Read and subgrid SSB dwelling data
                    !SSB_data_type=dwelling_index
                    !call uEMEP_read_SSB_data
                    if (use_RWC_emission_data) then
                        call uEMEP_read_RWC_heating_data
                    else
                        SSB_data_type=dwelling_index
                        call uEMEP_read_SSB_data
                    endif
                    
                endif

                !Read and subgrid agriculture data
                if (calculate_source(agriculture_index)) then
                    !Currently only data from RIVM here
                    call uEMEP_read_agriculture_rivm_data
                endif

                !Read in population data
                if (calculate_population_exposure_flag) then
                    !Read and subgrid SSB population data
                    SSB_data_type=population_data_type
                    call uEMEP_read_SSB_data
                endif

                !Autogrid setting for selecting which subgrids to calculate
                call uEMEP_auto_subgrid
    
                !Specify the subgrids sizes to be calculated using use_receptor_region
                call uEMEP_grid_receptor_data
                
                !Carry out tiling. Programme will stop here
                if (calculate_tiling_flag) then
                    call uEMEP_grid_roads
                    call uEMEP_set_tile_grids
                endif

            endif
    
            !Read time profiles for emissions
            call uEMEP_read_time_profiles

            !Call grid_roads again to include the time variation from NORTRIP
            call uEMEP_grid_roads

            !Interpolate meteo data to subgrid. Placed on the integral subgrid
            call uEMEP_subgrid_meteo_EMEP
        
            !Replaces proxy emissions with distributed EMEP emissions
            call uEMEP_subgrid_emission_EMEP  
               
            !Convert proxies to emissions including time profiles
            call uEMEP_convert_proxy_to_emissions
    
            !Grid dispersion
            traveltime_subgrid=0.
            do source_index=1,n_source_index
            if (calculate_source(source_index)) then
                call uEMEP_subgrid_dispersion(source_index)
            endif
            enddo
    
            !Diagnostic for comparing EMEP and proxy data emissions
            call uEMEP_aggregate_proxy_emission_in_EMEP_grid
    
            !Put EMEP data into subgrids for all sources
            call uEMEP_subgrid_EMEP
        
            !Interpolate EMEP to sub-grid
            do source_index=1,n_source_index
            if (calculate_source(source_index)) then
                !Redistributes dispersed proxies into the EMEP grid concentrations only when local_subgrid_method_flag=1
                call uEMEP_redistribute_local_source(source_index)
                !Places dispersed proxies into the local subgrid concentrations when local_subgrid_method_flag<>1
                call uEMEP_disperse_local_source(source_index)
            endif
            enddo
    
            !Combine and save sources in local and total values
            call uEMEP_combine_local_source
    
            !Calculate chemistry for NO2
            call uEMEP_chemistry

            !Calculate exposure
            if (calculate_population_exposure_flag) then
                call uEMEP_calculate_exposure
            endif
    
            !if (t_loop.eq.end_time_loop_index.and.g_loop.eq.end_grid_loop_index) then
            !    call CPU_TIME(end_time_cpu)
            !    write(*,*) ''
            !    write(*,'(a,i3,a,i2)') 'CPU time taken until saving (MM:SS): ',floor((end_time_cpu-start_time_cpu)/60.),':',floor(mod(end_time_cpu-start_time_cpu,60.))
            !endif
            
            if (save_netcdf_file_flag.or.save_netcdf_receptor_flag) then
                call uEMEP_save_netcdf_control
            endif
            
    
        enddo !t_loop
    
        if (first_g_loop) first_g_loop=.false.
    
    endif !use_receptor
    enddo !g_loop

    call CPU_TIME(end_time_cpu)

    if (unit_logfile.ne.0) then 
    write(unit_logfile,*) ''
    write(unit_logfile,*) '------------------------------------------------------------------------'
    write(unit_logfile,*) 'Ending programm uEMEP'
    write(unit_logfile,'(a,i3,a,i2)') ' CPU time taken (MM:SS): ',floor((end_time_cpu-start_time_cpu)/60.),':',floor(mod(end_time_cpu-start_time_cpu,60.))
    write(unit_logfile,*) '------------------------------------------------------------------------'
    endif
    
    if (unit_logfile.gt.0) then
         close(unit_logfile,status='keep')
    endif


    write(*,*) ''
    write(*,*) '------------------------------------------------------------------------'
    write(*,*) 'Ending programm uEMEP'
    write(*,'(a,i3,a,i2)') ' CPU time taken (MM:SS): ',floor((end_time_cpu-start_time_cpu)/60.),':',floor(mod(end_time_cpu-start_time_cpu,60.))
    write(*,*) '------------------------------------------------------------------------'

    end program uEMEP_v3

