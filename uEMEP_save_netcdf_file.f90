!Saves data in netcdf format
    
    subroutine uEMEP_save_netcdf_control
    
    use uEMEP_definitions

    implicit none
    
    integer i,j,k,t,l
    integer i_comp,i_file,i_meteo
    character(256) temp_name,unit_str,title_str,title_str_rec,var_name_temp, temp_date_str,station_name_str,temp_name_rec,temp_compound_str
    logical create_file,create_file_rec
    integer i_source
    real :: valid_min=0.
    real, allocatable :: temp_subgrid(:,:,:)
    real, allocatable :: temp_integral_subgrid(:,:,:)
    real, allocatable :: exhaust_subgrid(:,:,:)
    real, allocatable :: aqi_subgrid(:,:,:,:)
    integer, allocatable :: aqi_responsible_pollutant_index(:,:,:)
    real, allocatable :: temp_subgrid_ascii(:,:)
    
    integer ii,jj,tt
    
    real aqi_limits_temp(n_compound_index,1:5)
    real max_aqi
    integer n_aqi_pollutant_index
    parameter (n_aqi_pollutant_index=4)
    integer aqi_pollutant_index(n_aqi_pollutant_index)
    integer i_pollutant,i_loop
    character(256) variable_type
    logical :: receptor_available=.true.
    real scale_factor
    integer n_save_aqi_pollutant_index
    real temp_sum_comp
    integer count
    character(256) filename_ascii
    
    integer i_sp,ii_sp
    
    !Functions
    real DIRECTION
    
    if (include_o3_in_aqi_index) then
        n_save_aqi_pollutant_index=n_aqi_pollutant_index
    else
        n_save_aqi_pollutant_index=n_aqi_pollutant_index-1
    endif   
    
    if (.not.allocated(temp_subgrid)) allocate(temp_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index)))
    if (.not.allocated(temp_integral_subgrid)) allocate(temp_integral_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),subgrid_dim(t_dim_index)))
    if (.not.allocated(exhaust_subgrid)) allocate(exhaust_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index)))
    if (.not.allocated(aqi_subgrid).and.save_aqi) allocate(aqi_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index))
    if (.not.allocated(aqi_responsible_pollutant_index).and.save_aqi) allocate(aqi_responsible_pollutant_index(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index)))
    if (.not.allocated(temp_subgrid_ascii).and.save_compounds_as_ascii) allocate(temp_subgrid_ascii(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index)))
    
    !Save subgrid calculations
    valid_min=0.   
    unit_str="ug/m3"
    variable_type='byte'
    !variable_type='double'
    variable_type='float'
    scale_factor=1.
    
    !Save the data
    i_file=subgrid_total_file_index(allsource_index)
    !temp_name=trim(pathname_grid(i_file))//trim(filename_grid(i_file))//'_'//trim(file_tag)//'.nc'
    !if (len(config_date_str).gt.0) then
    !    temp_date_str='_'//trim(config_date_str)
    !else
    !    temp_date_str=''
    !endif
    
    if (len(filename_date_output_grid).gt.0) then
        temp_date_str='_'//filename_date_output_grid
    else
        temp_date_str=''
    endif
    
    if (use_multiple_receptor_grids_flag) then
        station_name_str=trim(name_receptor_in(g_loop,1))//'_';
    else
        station_name_str=''
    endif
    
    !Do not write 'all' in file name if all compounds are selected
    if (pollutant_index.eq.all_nc_index) then
        temp_compound_str=''
    else
        temp_compound_str='_'//trim(var_name_nc(conc_nc_index,compound_index,allsource_index))
    endif
    !Do not write anything about the compounds
    temp_compound_str=''
    
    !temp_name=trim(pathname_grid(i_file))//trim(station_name_str)//'uEMEP_'//trim(file_tag)//trim(temp_compound_str)//trim(temp_date_str)//'_'//trim(forecast_hour_str)//'.nc'
    !temp_name_rec=trim(pathname_grid(i_file))//'uEMEP_'//trim(file_tag)//'_station'//trim(temp_compound_str)//trim(temp_date_str)//'_'//trim(forecast_hour_str)//'.nc'
    !if (save_netcdf_average_flag) then
    !temp_name=trim(pathname_grid(i_file))//trim(station_name_str)//'uEMEP_'//trim(file_tag)//trim(temp_compound_str)//'_mean'//trim(temp_date_str)//'_'//trim(forecast_hour_str)//'.nc'
    !temp_name_rec=trim(pathname_grid(i_file))//'uEMEP_'//trim(file_tag)//'_station'//trim(temp_compound_str)//'_mean'//trim(temp_date_str)//'_'//trim(forecast_hour_str)//'.nc'
    !endif
    temp_name=trim(pathname_grid(i_file))//trim(station_name_str)//'uEMEP_'//trim(file_tag)//trim(temp_compound_str)//trim(temp_date_str)//'.nc'
    temp_name_rec=trim(pathname_grid(i_file))//'uEMEP_'//trim(file_tag)//'_station'//trim(temp_compound_str)//trim(temp_date_str)//'.nc'
    if (save_netcdf_average_flag) then
    temp_name=trim(pathname_grid(i_file))//trim(station_name_str)//'uEMEP_'//trim(file_tag)//trim(temp_compound_str)//'_mean'//trim(temp_date_str)//'.nc'
    temp_name_rec=trim(pathname_grid(i_file))//'uEMEP_'//trim(file_tag)//'_station'//trim(temp_compound_str)//'_mean'//trim(temp_date_str)//'.nc'
    endif
       
    write(unit_logfile,'(A)') ''
    write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Saving netcdf data (uEMEP_save_netcdf_control)'
	write(unit_logfile,'(A)') '================================================================'
    
    i_comp=1
    if (save_netcdf_receptor_flag.and.n_valid_receptor.eq.0) then
        if (i_comp.eq.1.and.t_loop.eq.start_time_loop_index) then
            write(unit_logfile,'(a)')'No receptor positions available. Will not save receptor data.'
            receptor_available=.false.
        endif
    endif
    
    if (save_netcdf_average_flag) then
        if (.not.allocated(val_array_av)) then
            allocate(val_array_av(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),n_var_av))
            val_array_av=0.
        endif
        if (.not.allocated(time_seconds_output_av)) then
            allocate(time_seconds_output_av(n_var_av))
            time_seconds_output_av=0
        endif
        !Reset average file counter for every time step saving occurs
        counter_av=0
        write(unit_logfile,*) 'Saving as mean data'
    endif
    
    !Save the final result of the subgrid calculation in ascii format
    !This should only be used for annual mean calculations since it cannot have a time dimmension
    !Intended for FAIRMODE output
    !Makes a new file for each compound and averages over time, if there is a time dimension
    if (save_compounds_as_ascii) then
    variable_type='float'
    unit_str="ug/m3"
    do i_pollutant=1,n_pollutant_loop
        if (pollutant_loop_index(i_pollutant).ne.pmex_index &
            .and.pollutant_loop_index(i_pollutant).ne.pm25_sand_index.and.pollutant_loop_index(i_pollutant).ne.pm10_sand_index &
            .and.pollutant_loop_index(i_pollutant).ne.pm25_salt_index.and.pollutant_loop_index(i_pollutant).ne.pm10_salt_index) then
            do i_loop=1,n_pollutant_compound_loop(i_pollutant)
            
            i_comp=pollutant_compound_loop_index(i_pollutant,i_loop)
            !write(*,*) i_comp
        
            var_name_temp=trim(var_name_nc(conc_nc_index,i_comp,allsource_index))
            
            title_str=trim(var_name_temp)//'_'//trim(file_tag)//trim(temp_date_str)
            write(unit_logfile,'(a)')'Writing ascii data to: '//trim(title_str)
        
            filename_ascii=trim(pathname_output_grid)//trim(title_str)//'.asc'
        
            temp_subgrid_ascii(:,:)=sum(comp_subgrid(:,:,:,i_comp),3)/subgrid_dim(t_dim_index)
        

            write(unit_logfile,'(a,f12.2)')'Writing ascii array variable: '//trim(var_name_temp),sum(comp_subgrid(:,:,:,i_comp))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index)

            call write_esri_ascii_file(unit_logfile,filename_ascii,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_delta(x_dim_index),temp_subgrid_ascii,x_subgrid,y_subgrid)
            
            enddo
        endif
    enddo
    endif

    !Save the final result of the subgrid calculation
    if (save_compounds) then
    variable_type='float'
    unit_str="ug/m3"
    do i_pollutant=1,n_pollutant_loop
        if (pollutant_loop_index(i_pollutant).ne.pmex_index &
            .and.pollutant_loop_index(i_pollutant).ne.pm25_sand_index.and.pollutant_loop_index(i_pollutant).ne.pm10_sand_index &
            .and.pollutant_loop_index(i_pollutant).ne.pm25_salt_index.and.pollutant_loop_index(i_pollutant).ne.pm10_salt_index) then
        do i_loop=1,n_pollutant_compound_loop(i_pollutant)
            
        i_comp=pollutant_compound_loop_index(i_pollutant,i_loop)
        !write(*,*) i_comp
        
        if (i_pollutant.eq.1.and.i_loop.eq.1.and.t_loop.eq.start_time_loop_index.and.save_netcdf_file_flag) then
            create_file=.true.
            title_str='uEMEP_concentration_'//trim(file_tag)//temp_date_str
            write(unit_logfile,'(a)')'Writing to: '//trim(temp_name)
        else
            create_file=.false.
        endif

        if (i_pollutant.eq.1.and.i_loop.eq.1.and.t_loop.eq.start_time_loop_index.and.first_g_loop.and.save_netcdf_receptor_flag) then
            create_file_rec=.true.
            title_str_rec='uEMEP_receptor_'//trim(file_tag)//temp_date_str
            if (receptor_available) write(unit_logfile,'(a)')'Writing to: '//trim(temp_name_rec)
        else
            create_file_rec=.false.
        endif
        
        var_name_temp=trim(var_name_nc(conc_nc_index,i_comp,allsource_index))//'_concentration'
        if (save_netcdf_file_flag) then
            temp_sum_comp=0.
            count=0
            do j=1,subgrid_dim(y_dim_index)
            do i=1,subgrid_dim(x_dim_index)
            if (use_subgrid(i,j,allsource_index)) then
                temp_sum_comp=temp_sum_comp+sum(comp_subgrid(i,j,:,i_comp))/subgrid_dim(t_dim_index)
                count=count+1
            endif
            enddo
            enddo
            if (count.gt.0) then
                temp_sum_comp=temp_sum_comp/count
            else
                temp_sum_comp=0
            endif            
            write(unit_logfile,'(a,f12.2)')'Writing netcdf array variable:    '//trim(var_name_temp),temp_sum_comp
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
            ,comp_subgrid(:,:,:,i_comp),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
            ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
            write(unit_logfile,'(a,f12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(comp_subgrid(:,:,:,i_comp))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
            ,comp_subgrid(:,:,:,i_comp),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
            ,unit_str,title_str_rec,create_file_rec,valid_min &
            ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
            ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
            ,z_rec(allsource_index,1) &
            ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
        endif
        enddo
        endif
    enddo
    endif
    
    create_file=.false.
    create_file_rec=.false.
    
    
    
    !Save the different local source contributions, not the total though
    if (save_source_contributions) then
        if (save_netcdf_fraction_as_contribution_flag) then
            variable_type='float'
            unit_str="ug/m3"
        else           
            variable_type='byte'
            unit_str="%"
        endif
        
    do i_pollutant=1,n_pollutant_loop
    do i_source=1,n_source_index
        !if (calculate_source(i_source).or.i_source.eq.allsource_index) then
        !Don't save any exhaust, sand or salt pollutant sources. Dealt with later
        if (calculate_source(i_source).and.i_source.ne.allsource_index.and.pollutant_loop_index(i_pollutant).ne.pmex_index &
            .and.pollutant_loop_index(i_pollutant).ne.pm25_sand_index.and.pollutant_loop_index(i_pollutant).ne.pm10_sand_index &
            .and.pollutant_loop_index(i_pollutant).ne.pm25_salt_index.and.pollutant_loop_index(i_pollutant).ne.pm10_salt_index) then
        
            i_file=subgrid_local_file_index(i_source)
            
            !Only save nonexhaust pm and all the other sources
            if ((pollutant_loop_index(i_pollutant).eq.nox_index.and.i_source.eq.traffic_index)) then
            !if ((pollutant_loop_index(i_pollutant).ne.nox_index.or.i_source.ne.traffic_index)) then
            else
                if (i_source.eq.traffic_index.and.(pollutant_loop_index(i_pollutant).eq.pm10_index.or.pollutant_loop_index(i_pollutant).eq.pm25_index)) then
                    var_name_temp=trim(var_name_nc(conc_nc_index,pollutant_loop_index(i_pollutant),allsource_index))//'_'//trim(filename_grid(i_file))//'_nonexhaust'
                    if (save_netcdf_fraction_as_contribution_flag) then
                        temp_subgrid=(subgrid(:,:,:,local_subgrid_index,i_source,i_pollutant)-subgrid(:,:,:,local_subgrid_index,i_source,pollutant_loop_back_index(pmex_index)))
                    else
                        temp_subgrid=(subgrid(:,:,:,local_subgrid_index,i_source,i_pollutant)-subgrid(:,:,:,local_subgrid_index,i_source,pollutant_loop_back_index(pmex_index)))/subgrid(:,:,:,total_subgrid_index,allsource_index,i_pollutant)*100.
                    endif
                else
                    var_name_temp=trim(var_name_nc(conc_nc_index,pollutant_loop_index(i_pollutant),allsource_index))//'_'//trim(filename_grid(i_file))
                    if (save_netcdf_fraction_as_contribution_flag) then
                        temp_subgrid=subgrid(:,:,:,local_subgrid_index,i_source,i_pollutant)
                    else
                        temp_subgrid=subgrid(:,:,:,local_subgrid_index,i_source,i_pollutant)/subgrid(:,:,:,total_subgrid_index,allsource_index,i_pollutant)*100.
                    endif
                endif
                !In case of any 0 concentrations
                where (subgrid(:,:,:,total_subgrid_index,allsource_index,i_pollutant).eq.0) temp_subgrid=0
            
                if (save_netcdf_file_flag) then
                    write(unit_logfile,'(a,f12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                endif
                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                    write(unit_logfile,'(a,f12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str_rec,create_file_rec,valid_min &
                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,z_rec(allsource_index,1) &
                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
                endif
            
            endif
            
            !Special case for exhaust as this must be given as a fraction for both PM2.5 and PM10.
            !if ((pollutant_loop_index(i_pollutant).eq.pm10_index.or.pollutant_loop_index(i_pollutant).eq.pm25_index).and.i_source.eq.traffic_index) then
            !Write all exhaust pollutants except the pmex
            if (i_source.eq.traffic_index) then
            
                !If PM then exhaust output is exhaust
                if (pollutant_loop_index(i_pollutant).eq.pm10_index.or.pollutant_loop_index(i_pollutant).eq.pm25_index) then
                    exhaust_subgrid=subgrid(:,:,:,local_subgrid_index,i_source,pollutant_loop_back_index(pmex_index))
                else
                    exhaust_subgrid=subgrid(:,:,:,local_subgrid_index,i_source,i_pollutant)            
                endif
            
                if (save_netcdf_fraction_as_contribution_flag) then
                    var_name_temp=trim(var_name_nc(conc_nc_index,pollutant_loop_index(i_pollutant),allsource_index))//'_'//'local_contribution_traffic_exhaust'
                    temp_subgrid=exhaust_subgrid
                else
                    var_name_temp=trim(var_name_nc(conc_nc_index,pollutant_loop_index(i_pollutant),allsource_index))//'_'//'local_fraction_traffic_exhaust'
                    temp_subgrid=exhaust_subgrid/subgrid(:,:,:,total_subgrid_index,allsource_index,i_pollutant)*100.
                endif
                where (subgrid(:,:,:,total_subgrid_index,allsource_index,i_pollutant).eq.0) temp_subgrid=0
              
                if (save_netcdf_file_flag) then
                    write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                endif
                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                    write(unit_logfile,'(a,f12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str_rec,create_file_rec,valid_min &
                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,z_rec(allsource_index,1) &
                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
                endif
            
            endif
        endif
            
            !Special case for salt and sand
            if (pollutant_loop_index(i_pollutant).eq.pm25_sand_index.or.pollutant_loop_index(i_pollutant).eq.pm10_sand_index &
            .or.pollutant_loop_index(i_pollutant).eq.pm25_salt_index.or.pollutant_loop_index(i_pollutant).eq.pm10_salt_index) then
                !Save the total nonlocal part from EMEP
                if (i_source.eq.traffic_index) then
               
                    i_file=subgrid_local_file_index(i_source)
                    var_name_temp=trim(var_name_nc(conc_nc_index,pollutant_loop_index(i_pollutant),allsource_index))//'_'//trim(filename_grid(i_file))
                    if (save_netcdf_fraction_as_contribution_flag) then
                        temp_subgrid=subgrid(:,:,:,local_subgrid_index,i_source,i_pollutant)
                    else
                        if (pollutant_loop_index(i_pollutant).eq.pm10_sand_index.or.pollutant_loop_index(i_pollutant).eq.pm10_salt_index) then
                            temp_subgrid=subgrid(:,:,:,local_subgrid_index,i_source,i_pollutant)/subgrid(:,:,:,total_subgrid_index,allsource_index,pollutant_loop_back_index(pm10_index))*100.
                        else
                            temp_subgrid=subgrid(:,:,:,local_subgrid_index,i_source,i_pollutant)/subgrid(:,:,:,total_subgrid_index,allsource_index,pollutant_loop_back_index(pm25_index))*100.
                        endif
                    
                    endif
                    where (subgrid(:,:,:,total_subgrid_index,allsource_index,i_pollutant).eq.0) temp_subgrid=0

                    if (save_netcdf_file_flag) then
                        write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
                        call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                            ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                            ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                    endif
                    if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                    write(unit_logfile,'(a,f12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                        call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                            ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                            ,unit_str,title_str_rec,create_file_rec,valid_min &
                            ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                            ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                            ,z_rec(allsource_index,1) &
                            ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
                    endif
                
                endif
            endif
            
            if (pollutant_loop_index(i_pollutant).ne.pmex_index &
            .and.pollutant_loop_index(i_pollutant).ne.pm25_sand_index.and.pollutant_loop_index(i_pollutant).ne.pm10_sand_index &
            .and.pollutant_loop_index(i_pollutant).ne.pm25_salt_index.and.pollutant_loop_index(i_pollutant).ne.pm10_salt_index) then
            !Save the total nonlocal part from EMEP
            if (i_source.eq.allsource_index) then
               
                i_file=emep_subgrid_nonlocal_file_index(i_source)
                var_name_temp=trim(var_name_nc(conc_nc_index,pollutant_loop_index(i_pollutant),allsource_index))//'_'//trim(filename_grid(i_file))
                if (save_netcdf_fraction_as_contribution_flag) then
                    temp_subgrid=subgrid(:,:,:,emep_nonlocal_subgrid_index,i_source,i_pollutant)
                else
                    temp_subgrid=subgrid(:,:,:,emep_nonlocal_subgrid_index,i_source,i_pollutant)/subgrid(:,:,:,total_subgrid_index,allsource_index,i_pollutant)*100.
                endif
                where (subgrid(:,:,:,total_subgrid_index,allsource_index,i_pollutant).eq.0) temp_subgrid=0

                if (save_netcdf_file_flag) then
                    write(unit_logfile,'(a,f12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                endif
                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str_rec,create_file_rec,valid_min &
                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,z_rec(allsource_index,1) &
                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
                endif
                
            endif
        endif
    enddo
    enddo
    endif

    if (save_no2_source_contributions.or.save_o3_source_contributions) then       
        call uEMEP_source_fraction_chemistry
    endif
    
    if (save_no2_source_contributions) then
    
        if (save_netcdf_fraction_as_contribution_flag) then
            variable_type='float'
            unit_str="ug/m3"
        else           
            variable_type='byte'
            unit_str="%"
        endif

    do i_source=1,n_source_index
        !if (calculate_source(i_source).or.i_source.eq.allsource_index) then
        if (calculate_source(i_source).or.i_source.eq.allsource_index) then
        
            if (i_source.eq.allsource_index) then
                i_file=emep_subgrid_nonlocal_file_index(i_source)
            else
                i_file=subgrid_local_file_index(i_source)
            endif
            
            if (i_source.eq.traffic_index) then
                var_name_temp=trim(var_name_nc(conc_nc_index,no2_nc_index,allsource_nc_index))//'_'//trim(filename_grid(i_file))//'_exhaust'
            else
                var_name_temp=trim(var_name_nc(conc_nc_index,no2_nc_index,allsource_nc_index))//'_'//trim(filename_grid(i_file))
            endif
            if (save_netcdf_fraction_as_contribution_flag) then
                temp_subgrid=comp_source_fraction_subgrid(:,:,:,no2_index,i_source)*comp_subgrid(:,:,:,no2_index)
            else
                temp_subgrid=comp_source_fraction_subgrid(:,:,:,no2_index,i_source)*100.
            endif
            
            
            if (save_netcdf_file_flag) then
                write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
                call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                    ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
            endif
            if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                    ,unit_str,title_str_rec,create_file_rec,valid_min &
                    ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,z_rec(allsource_index,1) &
                    ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
            endif
        endif
    enddo
    
    endif
    
    if (save_o3_source_contributions) then
           
        if (save_netcdf_fraction_as_contribution_flag) then
            variable_type='float'
            unit_str="ug/m3"
        else           
            variable_type='byte'
            unit_str="%"
        endif
    
    valid_min=0.

    do i_source=1,n_source_index
        !if (calculate_source(i_source).or.i_source.eq.allsource_index) then
        !Only save the nonlocal part as 100%
        if (i_source.eq.allsource_index) then
        
            if (i_source.eq.allsource_index) then
                i_file=emep_subgrid_nonlocal_file_index(i_source)
            else
                i_file=subgrid_local_file_index(i_source)
            endif
            
            !Temporary measure is to set nonlocal to 100% and the sources to 0.
            !This correctly says that all ozone comes from outside, it just does not say how much is removed by the local sources
            if (i_source.eq.allsource_index) then
                comp_source_fraction_subgrid(:,:,:,o3_index,i_source)=1.
            else
                comp_source_fraction_subgrid(:,:,:,o3_index,i_source)=0.
            endif
             
            var_name_temp=trim(var_name_nc(conc_nc_index,o3_nc_index,allsource_nc_index))//'_'//trim(filename_grid(i_file))
            if (save_netcdf_fraction_as_contribution_flag) then
                temp_subgrid=comp_subgrid(:,:,:,o3_index)
            else
                temp_subgrid=comp_source_fraction_subgrid(:,:,:,o3_index,i_source)*100.
            endif
            
            if (save_netcdf_file_flag) then
                write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
                call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                    ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
            endif
            if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                    ,unit_str,title_str_rec,create_file_rec,valid_min &
                    ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,z_rec(allsource_index,1) &
                    ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
            endif
        endif
    enddo
    
    valid_min=0.   
    
    endif
    
    !Save the emissions interpolated to the target grid
    if (save_emissions) then
    variable_type='float'
    unit_str="g/s"   
    do i_pollutant=1,n_pollutant_loop
    do i_source=1,n_source_index
        if (calculate_source(i_source)) then
    
            !Do not save emissions if the source is not traffic and the pollutant is road sand and salt or exhaust as these do not exist
            if (i_source.ne.traffic_index.and.(pollutant_loop_index(i_pollutant).eq.pm25_sand_index.or.pollutant_loop_index(i_pollutant).eq.pm10_sand_index &
                .or.pollutant_loop_index(i_pollutant).eq.pm25_salt_index.or.pollutant_loop_index(i_pollutant).eq.pm10_salt_index &
                .or.pollutant_loop_index(i_pollutant).eq.pmex_index)) then
                !Do nothing because these pollutants only exist for traffic
            else
                
                i_file=emission_file_index(i_source)
                var_name_temp=trim(var_name_nc(conc_nc_index,pollutant_loop_index(i_pollutant),allsource_index))//'_'//trim(filename_grid(i_file))
                
                !Calculate the emissions in the target grid
                temp_subgrid=0.
                do j=1,subgrid_dim(y_dim_index)
                do i=1,subgrid_dim(x_dim_index)
                
                    
                    ii=crossreference_target_to_emission_subgrid(i,j,x_dim_index,i_source)
                    jj=crossreference_target_to_emission_subgrid(i,j,y_dim_index,i_source)
            
                    temp_subgrid(i,j,:)=emission_subgrid(ii,jj,:,i_source,i_pollutant)
                    !Subgrid emissions, if relevant, are in units ug/sec/subgrid. Convert to g/s. Acount for the difference in subgrid sizes here
                    temp_subgrid(i,j,:)=1.0e-6*temp_subgrid(i,j,:)*(subgrid_delta(y_dim_index)*subgrid_delta(x_dim_index)) &
                        /(emission_subgrid_delta(y_dim_index,i_source)*emission_subgrid_delta(x_dim_index,i_source))
 
                enddo
                enddo
                
                if (save_netcdf_file_flag) then
                    write(unit_logfile,'(a,es12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                endif
                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,es12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str_rec,create_file_rec,valid_min &
                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,z_rec(allsource_index,1) &
                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
                endif
            
            endif
            
        endif
    enddo
    enddo
    endif

    !Save population interpolated to the target grid
    if (save_population) then
    variable_type='float'
    unit_str="inhabitants/grid"
    
                i_file=population_file_index(population_index)
                var_name_temp=trim(filename_grid(i_file))
                
                !Calculate the emissions in the target grid
                temp_subgrid=0.
                do j=1,subgrid_dim(y_dim_index)
                do i=1,subgrid_dim(x_dim_index)
                                   
                    ii=crossreference_target_to_population_subgrid(i,j,x_dim_index)
                    jj=crossreference_target_to_population_subgrid(i,j,y_dim_index)
            
                    temp_subgrid(i,j,:)=population_subgrid(ii,jj,population_index)
                    !Acount for the difference in subgrid sizes here
                    temp_subgrid(i,j,:)=temp_subgrid(i,j,:)*(subgrid_delta(y_dim_index)*subgrid_delta(x_dim_index)) &
                        /(population_subgrid_delta(y_dim_index)*population_subgrid_delta(x_dim_index))
 
                enddo
                enddo
                
                unit_str='inhabitants/grid'
                if (save_netcdf_file_flag) then
                    write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                endif
                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str_rec,create_file_rec,valid_min &
                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,z_rec(allsource_index,1) &
                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
                endif
    
    endif

    !Save the EMEP data interpolated to the subgrid
    if (save_emep_source_contributions) then
    do i_pollutant=1,n_pollutant_loop
    do i_source=1,n_source_index
        if (calculate_source(i_source).or.i_source.eq.allsource_index) then
            
            !Only save the allsource value here
            if (i_source.eq.allsource_index) then
                variable_type='float'
                i_file=emep_subgrid_file_index(i_source)
                var_name_temp=trim(var_name_nc(conc_nc_index,pollutant_loop_index(i_pollutant),allsource_index))//'_'//trim(filename_grid(i_file))
                unit_str='ug/m3'
                if (save_netcdf_file_flag) then
                    write(unit_logfile,'(a,f12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(subgrid(:,:,:,emep_subgrid_index,i_source,i_pollutant))/size(subgrid,1)/size(subgrid,2)/size(subgrid,3)
                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,subgrid(:,:,:,emep_subgrid_index,i_source,i_pollutant),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                endif
                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(subgrid(:,:,:,emep_subgrid_index,i_source,i_pollutant))/size(subgrid,1)/size(subgrid,2)/size(subgrid,3)
                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,subgrid(:,:,:,emep_subgrid_index,i_source,i_pollutant),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str_rec,create_file_rec,valid_min &
                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,z_rec(allsource_index,1) &
                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
                endif
                
            endif
    
                if (save_netcdf_fraction_as_contribution_flag) then
                    variable_type='float'
                    unit_str="ug/m3"
                else           
                    variable_type='byte'
                    unit_str="%"
                endif
                i_file=emep_subgrid_local_file_index(i_source)
                var_name_temp=trim(var_name_nc(conc_nc_index,pollutant_loop_index(i_pollutant),allsource_index))//'_'//trim(filename_grid(i_file))
                if (save_netcdf_fraction_as_contribution_flag) then
                    temp_subgrid=subgrid(:,:,:,emep_local_subgrid_index,i_source,i_pollutant)
                else
                    temp_subgrid=subgrid(:,:,:,emep_local_subgrid_index,i_source,i_pollutant)/subgrid(:,:,:,emep_subgrid_index,allsource_index,i_pollutant)*100.
                endif
                if (save_netcdf_file_flag) then
                    write(unit_logfile,'(a,f12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                endif
                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str_rec,create_file_rec,valid_min &
                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,z_rec(allsource_index,1) &
                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
                endif
                
            
        endif
    enddo
    enddo
    endif
    
    !Save the other interpolated EMEP compounds used for nox chemistry as well
    if (save_for_chemistry) then
    variable_type='float'
    do i_pollutant=1,n_emep_pollutant_loop
    if (pollutant_loop_index(i_pollutant).ne.pmex_index) then
        do i_loop=1,n_pollutant_compound_loop(i_pollutant)
            
        i_comp=pollutant_compound_loop_index(i_pollutant,i_loop)
        
        i_file=emep_subgrid_file_index(allsource_index)
        var_name_temp=trim(var_name_nc(conc_nc_index,i_comp,allsource_index))//'_interpolated'//'_'//trim(filename_grid(i_file))
        unit_str="ug/m3"
       if (save_netcdf_file_flag) then
            write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,comp_EMEP_subgrid(:,:,:,i_comp),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(comp_EMEP_subgrid(:,:,:,i_comp))/size(subgrid,1)/size(subgrid,2)/size(subgrid,3)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,comp_EMEP_subgrid(:,:,:,i_comp),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,z_rec(allsource_index,1) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
        endif
    enddo
    endif
    enddo
    endif
    
    if (save_emep_species) then
        variable_type='float'
        do i_sp=1,n_species_loop_index 
            ii_sp=species_loop_index(i_sp)
            if (i_sp.ne.sp_pm_index) then
                !Do not include the total
                do i_loop=1,n_pmxx_sp_index
                    if (i_loop.eq.pm25_sp_index.or.i_loop.eq.pm10_sp_index) then            

                if (i_loop.eq.pm25_sp_index) i_pollutant=pollutant_loop_back_index(pm25_index)
                if (i_loop.eq.pm10_sp_index) i_pollutant=pollutant_loop_back_index(pm10_index)
                
                if (save_netcdf_fraction_as_contribution_flag) then
                    variable_type='float'
                    var_name_temp=trim(species_name_nc(i_loop,ii_sp))//'_nonlocal_contribution'
                    unit_str="ug/m3"
                    temp_subgrid=species_EMEP_subgrid(:,:,:,i_loop,i_sp)
                else
                    variable_type='byte'
                    var_name_temp=trim(species_name_nc(i_loop,i_sp))//'_nonlocal_fraction'
                    unit_str="%"
                    temp_subgrid=species_EMEP_subgrid(:,:,:,i_loop,i_sp)/subgrid(:,:,:,total_subgrid_index,allsource_index,i_pollutant)*100.
                endif

                    if (save_netcdf_file_flag) then
                        write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid(:,:,:))/size(subgrid,1)/size(subgrid,2)/size(subgrid,3)
                        call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                            ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                            ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                    endif
                    if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                        write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid(:,:,:))/size(subgrid,1)/size(subgrid,2)/size(subgrid,3)
                        call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                            ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                            ,unit_str,title_str_rec,create_file_rec,valid_min &
                            ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                            ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                            ,z_rec(allsource_index,1) &
                            ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
                    endif
                    endif
                enddo
            endif
        enddo
    endif

    !Only save sea salt in this way when emep species are not saved in the general way
    if (save_seasalt.and..not.save_emep_species) then
        variable_type='float'
        i_sp=n_species_loop_index
        ii_sp=species_loop_index(i_sp)
        do i_loop=1,n_pmxx_sp_index
        if (i_loop.eq.pm25_sp_index.or.i_loop.eq.pm10_sp_index) then            

            if (i_loop.eq.pm25_sp_index) i_pollutant=pollutant_loop_back_index(pm25_index)
            if (i_loop.eq.pm10_sp_index) i_pollutant=pollutant_loop_back_index(pm10_index)
            if (save_netcdf_fraction_as_contribution_flag) then
                variable_type='float'
                var_name_temp=trim(species_name_nc(i_loop,ii_sp))//'_nonlocal_contribution'
                unit_str="ug/m3"
                temp_subgrid=species_EMEP_subgrid(:,:,:,i_loop,i_sp)
            else
                variable_type='byte'
                var_name_temp=trim(species_name_nc(i_loop,i_sp))//'_nonlocal_fraction'
                unit_str="%"
                temp_subgrid=species_EMEP_subgrid(:,:,:,i_loop,i_sp)/subgrid(:,:,:,total_subgrid_index,allsource_index,i_pollutant)*100.
            endif

            if (save_netcdf_file_flag) then
            write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid(:,:,:))/size(subgrid,1)/size(subgrid,2)/size(subgrid,3)
                call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                    ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
            endif
            if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid(:,:,:))/size(subgrid,1)/size(subgrid,2)/size(subgrid,3)
                call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                    ,unit_str,title_str_rec,create_file_rec,valid_min &
                    ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,z_rec(allsource_index,1) &
                    ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
            endif
        endif
        enddo
        
    endif

    !Save the original EMEP compounds
    if (save_emep_original) then
    variable_type='float'
    do i_pollutant=1,n_emep_pollutant_loop
    if (pollutant_loop_index(i_pollutant).ne.pmex_index) then
    do i_loop=1,n_pollutant_compound_loop(i_pollutant)

        i_comp=pollutant_compound_loop_index(i_pollutant,i_loop)

        i_file=emep_subgrid_file_index(allsource_index)
        var_name_temp=trim(var_name_nc(conc_nc_index,i_comp,allsource_index))//'_original_EMEP_concentration'
        unit_str="ug/m3"
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a,f12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(orig_EMEP_subgrid(:,:,:,i_comp))/size(subgrid,1)/size(subgrid,2)/size(subgrid,3)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,orig_EMEP_subgrid(:,:,:,i_comp),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(orig_EMEP_subgrid(:,:,:,i_comp))/size(subgrid,1)/size(subgrid,2)/size(subgrid,3)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,orig_EMEP_subgrid(:,:,:,i_comp),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,z_rec(allsource_index,1) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
        endif
    enddo
    endif
    enddo
    endif
    
    !Save AQI
    if (save_aqi) then
    variable_type='byte'
    scale_factor=0.1
        !Hard coded AQI limits
        aqi_pollutant_index(1)=no2_index;aqi_pollutant_index(2)=pm10_index;aqi_pollutant_index(3)=pm25_index;aqi_pollutant_index(4)=o3_index
        aqi_limits_temp(:,2:4)=aqi_hourly_limits(:,1:3)
        aqi_limits_temp(:,1)=0.
        aqi_limits_temp(:,5)=2.*aqi_hourly_limits(:,3)
        aqi_subgrid=0.
        
        do t=1,subgrid_dim(t_dim_index)
        do j=1,subgrid_dim(y_dim_index)
        do i=1,subgrid_dim(x_dim_index)
            max_aqi=0.
            do l=1,n_save_aqi_pollutant_index !pollutant (no2,pm10,pm2.5,o3)
                do k=1,4 !level 
                    if (comp_subgrid(i,j,t,aqi_pollutant_index(l)).ge.aqi_limits_temp(aqi_pollutant_index(l),k).and.comp_subgrid(i,j,t,aqi_pollutant_index(l)).lt.aqi_limits_temp(aqi_pollutant_index(l),k+1)) then
                        aqi_subgrid(i,j,t,aqi_pollutant_index(l))=k+(comp_subgrid(i,j,t,aqi_pollutant_index(l))-aqi_limits_temp(aqi_pollutant_index(l),k))/(aqi_limits_temp(aqi_pollutant_index(l),k+1)-aqi_limits_temp(aqi_pollutant_index(l),k))
                    endif
                enddo
                aqi_subgrid(i,j,t,aqi_pollutant_index(l))=min(aqi_subgrid(i,j,t,aqi_pollutant_index(l)),4.99)
                if (aqi_subgrid(i,j,t,aqi_pollutant_index(l)).gt.max_aqi) then
                    max_aqi=aqi_subgrid(i,j,t,aqi_pollutant_index(l))
                    aqi_responsible_pollutant_index(i,j,t)=l
                endif        
            !write(*,*)  i,j,t,aqi_responsible_pollutant_index(i,j,t),aqi_subgrid(i,j,t,aqi_pollutant_index(l)),max_aqi
            enddo
        enddo
        enddo        
        enddo
        
        do l=1,n_save_aqi_pollutant_index
            !write(unit_logfile,*)  'MAX AQI in time and space from '//trim(pollutant_file_str(aqi_pollutant_index(l)))//' = ',maxval(aqi_subgrid(:,:,:,aqi_pollutant_index(l)))
        enddo
        
        var_name_temp='AQI'
        unit_str='1'
        !Take the maximum of the pollutants
        temp_subgrid=maxval(aqi_subgrid(:,:,:,aqi_pollutant_index(1:n_save_aqi_pollutant_index)),4)/scale_factor
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid*scale_factor)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,z_rec(allsource_index,1) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
        endif
        
        
    do l=1,n_save_aqi_pollutant_index

        i_comp=aqi_pollutant_index(l)
        var_name_temp='AQI_'//trim(var_name_nc(conc_nc_index,i_comp,allsource_index))
        unit_str='1'
        !Take the maximum of the pollutants
        temp_subgrid=aqi_subgrid(:,:,:,i_comp)/scale_factor
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid*scale_factor)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,z_rec(allsource_index,1) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
        endif
        
    enddo
        
    !Reset scale_factor
    scale_factor=1.
    
    endif

    !Save deposition
    if (save_deposition.and.calculate_deposition_flag) then
    variable_type='float'
    do i_pollutant=1,n_pollutant_loop
    do i_source=1,n_source_index
    if (calculate_source(i_source)) then

        i_comp=pollutant_loop_index(i_pollutant)

        var_name_temp=trim(var_name_nc(conc_nc_index,i_comp,allsource_index))//'_local_dry_deposition_'//source_file_str(i_source)
        unit_str="mg/m2/hr"
        !Save in mg per hour, conversition from ug/m2/s
        temp_subgrid=subgrid(:,:,:,drydepo_local_subgrid_index,i_source,i_pollutant)/1000.*3600.
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a,es12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
            write(unit_logfile,'(a,es12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,z_rec(allsource_index,1) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
        endif

        var_name_temp=trim(var_name_nc(conc_nc_index,i_comp,allsource_index))//'_local_wet_deposition_'//source_file_str(i_source)
        unit_str="mg/m2/hr"
        !Save in mg per hour, conversition from ug/m2/s
        temp_subgrid=subgrid(:,:,:,wetdepo_local_subgrid_index,i_source,i_pollutant)/1000.*3600.
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a,es12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
            write(unit_logfile,'(a,es12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,z_rec(allsource_index,1) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
        endif

    endif
    enddo 
    
        var_name_temp=trim(var_name_nc(conc_nc_index,i_comp,allsource_index))//'_nonlocal_dry_deposition_'//source_file_str(allsource_index)
        unit_str="mg/m2/hr"
        !Save in mg per hour, conversition from ug/m2/s
        temp_subgrid=subgrid(:,:,:,drydepo_nonlocal_subgrid_index,allsource_index,i_pollutant)/1000.*3600.
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a,es12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
            write(unit_logfile,'(a,es12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,z_rec(allsource_index,1) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
        endif

        var_name_temp=trim(var_name_nc(conc_nc_index,i_comp,allsource_index))//'_total_dry_deposition_'//source_file_str(allsource_index)
        unit_str="mg/m2/hr"
        !Save in mg per hour, conversition from ug/m2/s
        temp_subgrid=(subgrid(:,:,:,drydepo_nonlocal_subgrid_index,allsource_index,i_pollutant)+subgrid(:,:,:,drydepo_local_subgrid_index,allsource_index,i_pollutant))/1000.*3600.
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a,es12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
            write(unit_logfile,'(a,es12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,z_rec(allsource_index,1) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
        endif

        var_name_temp=trim(var_name_nc(conc_nc_index,i_comp,allsource_index))//'_nonlocal_wet_deposition_'//source_file_str(allsource_index)
        unit_str="mg/m2/hr"
        !Save in mg per hour, conversition from ug/m2/s
        temp_subgrid=subgrid(:,:,:,wetdepo_nonlocal_subgrid_index,allsource_index,i_pollutant)/1000.*3600.
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a,es12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
            write(unit_logfile,'(a,es12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,z_rec(allsource_index,1) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
        endif

        var_name_temp=trim(var_name_nc(conc_nc_index,i_comp,allsource_index))//'_total_wet_deposition_'//source_file_str(allsource_index)
        unit_str="mg/m2/hr"
        !Save in mg per hour, conversition from ug/m2/s
        temp_subgrid=(subgrid(:,:,:,wetdepo_nonlocal_subgrid_index,allsource_index,i_pollutant)+subgrid(:,:,:,wetdepo_local_subgrid_index,allsource_index,i_pollutant))/1000.*3600.
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a,es12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
            write(unit_logfile,'(a,es12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,z_rec(allsource_index,1) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
        endif

        !Deposition velocity
        temp_subgrid=0.
        do j=1,subgrid_dim(y_dim_index)
        do i=1,subgrid_dim(x_dim_index)
            ii=crossreference_target_to_deposition_subgrid(i,j,x_dim_index);jj=crossreference_target_to_deposition_subgrid(i,j,y_dim_index)
            temp_subgrid(i,j,:)=deposition_subgrid(ii,jj,:,vd_index,i_pollutant)
        enddo
        enddo
        var_name_temp=trim(var_name_nc(conc_nc_index,i_comp,allsource_index))//'_deposition_velocity_'//source_file_str(allsource_index)
        unit_str="cm/s"
        !Save in cm/s, conversition from ug/m2/s
        temp_subgrid=temp_subgrid*100.
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a,es12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
            write(unit_logfile,'(a,es12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,z_rec(allsource_index,1) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
        endif

        !Landuse category
        if (read_landuse_flag) then
        temp_subgrid=0.
        do j=1,subgrid_dim(y_dim_index)
        do i=1,subgrid_dim(x_dim_index)
            ii=crossreference_target_to_deposition_subgrid(i,j,x_dim_index);jj=crossreference_target_to_deposition_subgrid(i,j,y_dim_index)
            temp_subgrid(i,j,:)=landuse_subgrid(ii,jj,grid_index)
        enddo
        enddo
        var_name_temp=trim(var_name_nc(conc_nc_index,i_comp,allsource_index))//'_EMEP_landuse_category'
        unit_str="1"
        variable_type='byte'
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a,es12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
            write(unit_logfile,'(a,es12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,z_rec(allsource_index,1) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
        endif
        endif
    
        variable_type='float'

        !Save the original emep wet and dry deposition
        var_name_temp=trim(var_name_nc(conc_nc_index,i_comp,allsource_index))//'_original_EMEP_wet_deposition_'//source_file_str(allsource_index)
        unit_str="mg/m2/hr"
        !Save in mg per hour, conversition from ug/m2/s
        temp_subgrid=orig_emep_deposition_subgrid(:,:,:,wetdepo_index,i_pollutant)
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a,es12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
            write(unit_logfile,'(a,es12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,z_rec(allsource_index,1) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
        endif
        var_name_temp=trim(var_name_nc(conc_nc_index,i_comp,allsource_index))//'_original_EMEP_dry_deposition_'//source_file_str(allsource_index)
        unit_str="mg/m2/hr"
        !Save in mg per hour, conversition from ug/m2/s
        temp_subgrid=orig_emep_deposition_subgrid(:,:,:,drydepo_index,i_pollutant)
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a,es12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
            write(unit_logfile,'(a,es12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,z_rec(allsource_index,1) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
        endif

    enddo

    endif

    !Save weighted travel time
    if (save_for_chemistry) then
        variable_type='float'
        var_name_temp='Weighted_travel_time'
        unit_str='seconds'
        i_pollutant=1 !Only save the first pollutant index for travel time, which is nox in all.
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,traveltime_subgrid(:,:,:,1,i_pollutant),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(traveltime_subgrid(:,:,:,1,i_pollutant))/size(subgrid,1)/size(subgrid,2)/size(subgrid,3)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,traveltime_subgrid(:,:,:,1,i_pollutant),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,z_rec(allsource_index,1) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
        endif
    endif
 
    !Save the meteo interpolated to the target grid
    valid_min=-1.e24   
    
    if (save_wind_vectors) then
        variable_type='float'
        i_file=subgrid_ugrid_file_index
        i_meteo=ugrid_subgrid_index
        unit_str="m/s"
        !The same for all--------------------
        var_name_temp=trim(filename_grid(i_file))
        temp_subgrid=0.
        do j=1,subgrid_dim(y_dim_index)
        do i=1,subgrid_dim(x_dim_index)
            ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
        enddo
        enddo
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,z_rec(allsource_index,1) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
        endif
        !The same for all--------------------

        i_file=subgrid_vgrid_file_index
        i_meteo=vgrid_subgrid_index
        unit_str="m/s"
        !The same for all--------------------
        var_name_temp=trim(filename_grid(i_file))
        temp_subgrid=0.
        do j=1,subgrid_dim(y_dim_index)
        do i=1,subgrid_dim(x_dim_index)
            ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
        enddo
        enddo
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,z_rec(allsource_index,1) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
        endif
        !The same for all--------------------
    endif
    
    if (save_for_chemistry) then
        variable_type='float'
        i_file=subgrid_J_file_index
        i_meteo=J_subgrid_index
        unit_str="1/s"
        !The same for all--------------------
        var_name_temp=trim(filename_grid(i_file))
        temp_subgrid=0.
        do j=1,subgrid_dim(y_dim_index)
        do i=1,subgrid_dim(x_dim_index)
            ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
        enddo
        enddo
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,z_rec(allsource_index,1) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
        endif
        !The same for all--------------------
    endif

    if (save_other_meteo) then
        variable_type='float'
        i_file=subgrid_hmix_file_index
        i_meteo=hmix_subgrid_index
        unit_str="m"
        !The same for all--------------------
        var_name_temp=trim(filename_grid(i_file))
        temp_subgrid=0.
        do j=1,subgrid_dim(y_dim_index)
        do i=1,subgrid_dim(x_dim_index)
            ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
        enddo
        enddo
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,z_rec(allsource_index,1) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
        endif
        !The same for all--------------------

        variable_type='float'
        i_file=subgrid_hmix_file_index
        i_meteo=precip_subgrid_index
        unit_str="mm/hr"
        !The same for all--------------------
        var_name_temp='precipitation'!trim(filename_grid(i_file))
        temp_subgrid=0.
        do j=1,subgrid_dim(y_dim_index)
        do i=1,subgrid_dim(x_dim_index)
            ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
        enddo
        enddo
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,z_rec(allsource_index,1) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
        endif
        !The same for all--------------------

        
        i_file=subgrid_kz_file_index
        i_meteo=kz_subgrid_index
        unit_str="m2/s"
        !The same for all--------------------
        var_name_temp=trim(filename_grid(i_file))
        temp_subgrid=0.
        do j=1,subgrid_dim(y_dim_index)
        do i=1,subgrid_dim(x_dim_index)
            ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
        enddo
        enddo
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,z_rec(allsource_index,1) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
        endif
        !The same for all--------------------
 
        if (hourly_calculations) then
        
        i_file=subgrid_FFgrid_file_index
        i_meteo=FFgrid_subgrid_index
        unit_str="m/s"
        !The same for all--------------------
        var_name_temp=trim(filename_grid(i_file))
        temp_subgrid=0.
        do j=1,subgrid_dim(y_dim_index)
        do i=1,subgrid_dim(x_dim_index)
            ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
        enddo
        enddo
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,z_rec(allsource_index,1) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
        endif
        !The same for all--------------------

        i_file=subgrid_t2m_file_index
        i_meteo=t2m_subgrid_index
        unit_str="Celcius"
        !The same for all--------------------
        var_name_temp=trim(filename_grid(i_file))
        temp_subgrid=0.
        do j=1,subgrid_dim(y_dim_index)
        do i=1,subgrid_dim(x_dim_index)
            ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)-273.13
        enddo
        enddo
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,z_rec(allsource_index,1) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
        endif
        !The same for all--------------------

        i_file=subgrid_DDgrid_file_index
        !i_meteo=DDgrid_subgrid_index !i_meteo not used in this conversion to wind direction
        unit_str="degrees"
        do tt=1,subgrid_dim(t_dim_index)
        do jj=1,integral_subgrid_dim(y_dim_index)
        do ii=1,integral_subgrid_dim(x_dim_index)
            temp_integral_subgrid(ii,jj,tt)=DIRECTION(meteo_subgrid(ii,jj,tt,ugrid_subgrid_index),meteo_subgrid(ii,jj,tt,vgrid_subgrid_index))        
        enddo
        enddo
        enddo
        !This is not converted to real degrees, but subgrid degrees, so is wrong but close
        
        !The same for all--------------------
        var_name_temp=trim(filename_grid(i_file))
        temp_subgrid=0.
        do j=1,subgrid_dim(y_dim_index)
        do i=1,subgrid_dim(x_dim_index)
            ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);
            !In this case this is not the same for all
            temp_subgrid(i,j,:)=temp_integral_subgrid(ii,jj,:)
        enddo
        enddo
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,z_rec(allsource_index,1) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
        endif
        !The same for all--------------------

        i_file=subgrid_FF10_file_index
        i_meteo=FF10_subgrid_index
        unit_str="m/s"
        !The same for all--------------------
        var_name_temp=trim(filename_grid(i_file))
        temp_subgrid=0.
        do j=1,subgrid_dim(y_dim_index)
        do i=1,subgrid_dim(x_dim_index)
            ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
        enddo
        enddo
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,z_rec(allsource_index,1) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
        endif
        !The same for all--------------------

        endif
    
        i_file=subgrid_invL_file_index
        i_meteo=invL_subgrid_index
        unit_str="1/m"
        !The same for all--------------------
        var_name_temp=trim(filename_grid(i_file))
        temp_subgrid=0.
        do j=1,subgrid_dim(y_dim_index)
        do i=1,subgrid_dim(x_dim_index)
            ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
        enddo
        enddo
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.4)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,z_rec(allsource_index,1) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
        endif
        !The same for all--------------------

        i_file=subgrid_logz0_file_index
        i_meteo=logz0_subgrid_index
        unit_str="log(m)"
        !The same for all--------------------
        var_name_temp=trim(filename_grid(i_file))
        temp_subgrid=0.
        do j=1,subgrid_dim(y_dim_index)
        do i=1,subgrid_dim(x_dim_index)
            ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
        enddo
        enddo
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,z_rec(allsource_index,1) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
        endif
        !The same for all--------------------

        i_file=subgrid_ustar_file_index
        i_meteo=ustar_subgrid_index
        unit_str="m/s"
        !The same for all--------------------
        var_name_temp=trim(filename_grid(i_file))
        temp_subgrid=0.
        do j=1,subgrid_dim(y_dim_index)
        do i=1,subgrid_dim(x_dim_index)
            ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
        enddo
        enddo
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,z_rec(allsource_index,1) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
        endif
        !The same for all--------------------

        if (annual_calculations) then
        
        i_file=subgrid_invFFgrid_file_index
        i_meteo=inv_FFgrid_subgrid_index
        unit_str="s/m"
        !The same for all--------------------
        var_name_temp=trim(filename_grid(i_file))
        temp_subgrid=0.
        do j=1,subgrid_dim(y_dim_index)
        do i=1,subgrid_dim(x_dim_index)
            ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
        enddo
        enddo
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,z_rec(allsource_index,1) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
        endif
        !The same for all--------------------
            
        i_file=subgrid_invFF10_file_index
        i_meteo=inv_FF10_subgrid_index
        unit_str="s/m"
        !The same for all--------------------
        var_name_temp=trim(filename_grid(i_file))
        temp_subgrid=0.
        do j=1,subgrid_dim(y_dim_index)
        do i=1,subgrid_dim(x_dim_index)
            ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
        enddo
        enddo
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,z_rec(allsource_index,1) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)          
        endif
        !The same for all--------------------

        endif
    
    endif
    
    !Deallocate the averaging arrays with each receptor grid
    if (use_multiple_receptor_grids_flag.and..not.use_single_time_loop_flag.and.save_netcdf_average_flag) then
        if (allocated(val_array_av)) deallocate(val_array_av)
        if (allocated(time_seconds_output_av)) deallocate(time_seconds_output_av)
        counter_av=0
    endif
   

    end subroutine uEMEP_save_netcdf_control
    
    
    subroutine uEMEP_save_netcdf_file(unit_logfile_in,filename_netcdf,nx,ny,nt_in,val_array_in,x_array,y_array,lon_array,lat_array,name_array,unit_array,title_str,create_file,valid_min,variable_type,scale_factor)
    
    use uEMEP_definitions
    use netcdf
    
    implicit none
    
    character(256) filename_netcdf,name_array,unit_array,title_str,temp_name
    integer unit_logfile_in
    integer nx,ny,nt_in
    real val_array(nx,ny,nt_in),val_array_in(nx,ny,nt_in)!,val_array_temp(nx,ny,nt)
    real x_array(nx,ny)
    real y_array(nx,ny)
    real lon_array(nx,ny)
    real lat_array(nx,ny) !,lat_array_temp(nx,ny)
    !real time_array(nt)
    real x_vector(nx)
    real y_vector(ny)
    logical create_file
    real valid_min
    character(256) variable_type
    real scale_factor
    
    integer ncid
    integer y_dimid,x_dimid,time_dimid
    integer y_varid,x_varid,lat_varid,lon_varid,val_varid,time_varid,proj_varid
    integer dimids3(3),dimids2(2),chunks3(3)
    integer n_dims(3)
    integer status
    integer nf90_type
    integer t
    integer n_time_total
    integer nt
    integer(8) time_seconds_output_nc(nt_in) !Need integer 8 for the averaging

    if (trim(variable_type).eq.'byte') nf90_type=NF90_BYTE
    if (trim(variable_type).eq.'short') nf90_type=NF90_SHORT
    if (trim(variable_type).eq.'float') nf90_type=NF90_FLOAT
    if (trim(variable_type).eq.'double') nf90_type=NF90_DOUBLE
    
    !Assumes x and y are the dimensions   
    x_vector=x_array(:,1)
    y_vector=y_array(1,:)
    
    n_time_total=end_time_nc_index-start_time_nc_index+1
    nt=nt_in
    val_array=val_array_in
    time_seconds_output_nc=time_seconds_output
    
    !Save averages only
    if (save_netcdf_average_flag) then
        counter_av=counter_av+1
        if (counter_av.gt.n_var_av) then
            write(unit_logfile_in,*) 'ERROR: Array size for saving averages (n_var_av) not large enough. Stopping'
            stop
        endif
        if (use_single_time_loop_flag) then
            val_array_av(:,:,counter_av)=val_array_av(:,:,counter_av)+val_array(:,:,nt) !nt=1 in this case
            time_seconds_output_av(counter_av)=time_seconds_output_av(counter_av)+time_seconds_output_nc(nt)
            if (t_loop.eq.end_time_loop_index) then
                val_array(:,:,nt)=val_array_av(:,:,counter_av)/n_time_total
                time_seconds_output_nc(nt)=time_seconds_output_av(counter_av)/n_time_total
            endif
            !write(unit_logfile_in,'(a,3i)') 'Saving as average single time loop (nt,counter_av):',nt,counter_av,time_seconds_output_nc(nt)
        else
            val_array_av(:,:,counter_av)=sum(val_array,3)/size(val_array,3)
            time_seconds_output_av(counter_av)=sum(time_seconds_output_nc)/size(time_seconds_output_nc,1)
            nt=1
            val_array(:,:,nt)=val_array_av(:,:,counter_av)
            time_seconds_output_nc(nt)=time_seconds_output_av(counter_av)
            !write(unit_logfile_in,'(a,3i)') 'Saving as average multiple time loop (nt,counter_av):',nt,counter_av,time_seconds_output_nc(nt)
            
        endif
        
    endif
    
    !Mask the regions if required
    if (use_region_select_and_mask_flag) then
        do t=1,nt
            where (use_subgrid_val(:,:,allsource_index).eq.outside_region_index) val_array(:,:,t)=NODATA_value
            where (.not.use_subgrid(:,:,allsource_index)) val_array(:,:,t)=NODATA_value
        enddo
    endif
    
    if (create_file) then
        !Create a netcdf file
        !call check(  nf90_create(filename_netcdf, nf90_clobber, ncid) )
        !call check(  nf90_create(filename_netcdf, NF90_HDF5, ncid) )
        call check(  nf90_create(filename_netcdf, IOR(NF90_HDF5, NF90_CLASSIC_MODEL), ncid) ) !New

        !Specify global attributes
        call check(  nf90_put_att(ncid, nf90_global, "Conventions", "CF-1.4" ) )
        call check(  nf90_put_att(ncid, nf90_global, "title", trim(title_str)) )
        call check(  nf90_put_att(ncid, nf90_global, "Model", "uEMEP" ) )        
    
        !Projection data
        if (projection_type.eq.UTM_projection_index) then
        call check(  nf90_def_var(ncid, "projection_utm", NF90_int, proj_varid) )
        call check(  nf90_put_att(ncid, proj_varid, "semi_major_axis", 6378137.0 ) )
        call check(  nf90_put_att(ncid, proj_varid, "inverse_flattening", 298.257222101 ) )

        call check(  nf90_put_att(ncid, proj_varid, "grid_mapping_name", "transverse_mercator" ) )
        call check(  nf90_put_att(ncid, proj_varid, "scale_factor_at_central_meridian", 0.9996 ) )
        call check(  nf90_put_att(ncid, proj_varid, "latitude_of_projection_origin", 0 ) )
        call check(  nf90_put_att(ncid, proj_varid, "false_easting", 500000. ) )
        call check(  nf90_put_att(ncid, proj_varid, "false_northing", 0. ) )
        call check(  nf90_put_att(ncid, proj_varid, "longitude_of_central_meridian", utm_lon0 ) )
        endif
        
        !Define the dimensions
        call check(  nf90_def_dim(ncid,"time",NF90_UNLIMITED, time_dimid) )
        call check(  nf90_def_dim(ncid, "y", ny, y_dimid) )
        call check(  nf90_def_dim(ncid, "x", nx, x_dimid) )
     
        !Define the dimension variables
        !call check(  nf90_def_var(ncid, "time", NF90_DOUBLE, time_dimid, time_varid) )
        call check(  nf90_def_var(ncid, "time", NF90_INT, time_dimid, time_varid) )
        call check(  nf90_def_var(ncid, "y", NF90_REAL, y_dimid, y_varid) )
        call check(  nf90_def_var(ncid, "x", NF90_REAL, x_dimid, x_varid) )
    
        !Define the values
        dimids3 = (/ x_dimid, y_dimid, time_dimid /)
        dimids2 = (/ x_dimid, y_dimid /)
        call check(  nf90_def_var(ncid, "lat", NF90_REAL, dimids2, lat_varid) )
        call check(  nf90_def_var(ncid, "lon", NF90_REAL, dimids2, lon_varid) )

        !Specify the units
        call check(  nf90_put_att(ncid, lat_varid, "units", "degrees_north") )
        call check(  nf90_put_att(ncid, lon_varid, "units", "degrees_east") )
        call check(  nf90_put_att(ncid, y_varid, "units", "m") )
        call check(  nf90_put_att(ncid, x_varid, "units", "m") )
        call check(  nf90_put_att(ncid, time_varid, "units", trim(unit_dim_nc(time_dim_nc_index))) )

        !Specify other dimension attributes
        call check(  nf90_put_att(ncid, y_varid, "standard_name", "projection_y_axis") )
        call check(  nf90_put_att(ncid, x_varid, "standard_name", "projection_x_axis") )
        call check(  nf90_put_att(ncid, y_varid, "axis", "Y") )
        call check(  nf90_put_att(ncid, x_varid, "axis", "X") )
  
        !Close the definitions
        call check( nf90_enddef(ncid) )
        
        !call check( nf90_put_var(ncid, time_varid, val_dim_nc(1:dim_length_nc(time_dim_nc_index),time_dim_nc_index)) )
        !call check( nf90_put_var(ncid, time_varid, time_seconds_output(1:dim_length_nc(time_dim_nc_index))) )
        call check( nf90_put_var(ncid, time_varid, time_seconds_output_nc(1:nt)) )
        call check( nf90_put_var(ncid, y_varid, y_vector) )
        call check( nf90_put_var(ncid, x_varid, x_vector) )
        call check( nf90_put_var(ncid, lat_varid, lat_array) )
        call check( nf90_put_var(ncid, lon_varid, lon_array) )

        call check( nf90_close(ncid) )
    
    endif
        
    !Do not save any average data for single time loops until the last time step where it will save the average
    if (save_netcdf_average_flag.and.use_single_time_loop_flag.and.t_loop.ne.end_time_loop_index) then
        !call check( nf90_close(ncid) )
        return
    endif
        
    !Add to the existing file
    call check( nf90_open(filename_netcdf, NF90_WRITE, ncid) )
    
    !Get the dimensions id from the existing file
    call check( nf90_inq_dimid(ncid,"time",time_dimid) )
    call check( nf90_inq_dimid(ncid, "y", y_dimid) )
    call check( nf90_inq_dimid(ncid, "x", x_dimid) )
    dimids3 = (/ x_dimid, y_dimid, time_dimid /)
    chunks3 = (/ nx, ny, 1 /) !New
    call check( nf90_inquire_dimension(ncid, dimids3(1), temp_name, n_dims(1)) )
    call check( nf90_inquire_dimension(ncid, dimids3(2), temp_name, n_dims(2)) )
    call check( nf90_inquire_dimension(ncid, dimids3(3), temp_name, n_dims(3)) )


    status=nf90_inq_varid(ncid, trim(name_array), val_varid)
    if (status.ne.nf90_NoErr) then
        call check( nf90_redef(ncid) )
        !if the variable does not exist then create a new one
        !write(*,*) 'Creating new: ',trim(name_array)
        call check( nf90_def_var(ncid, trim(name_array), nf90_type, dimids3, val_varid) )
        ! gzip level 3 compression and shuffling
        ! optional _FillValue for values which never have been written, unpacked value
        call check( nf90_def_var_chunking(ncid, val_varid, NF90_CHUNKED, chunks3) ) !New
        call check( nf90_def_var_deflate(ncid, val_varid, 1, 1, 3) ) !New
        call check( nf90_put_att(ncid, val_varid, "units", trim(unit_array)) )
    
        !Specify other variable attributes
        if (nf90_type.eq.NF90_byte) then
            call check(  nf90_put_att(ncid, val_varid, "missing_value", int1(NODATA_value) ) ) !New
            call check(  nf90_put_att(ncid, val_varid, "valid_min", int1(valid_min)) )
        elseif (nf90_type.eq.NF90_short) then
            call check(  nf90_put_att(ncid, val_varid, "missing_value", int2(NODATA_value) ) ) !New
            call check(  nf90_put_att(ncid, val_varid, "valid_min", int2(valid_min)) )
        else
            call check(  nf90_put_att(ncid, val_varid, "missing_value", NODATA_value ) ) !New
            call check(  nf90_put_att(ncid, val_varid, "valid_min", valid_min) )
        endif
        call check(  nf90_put_att(ncid, val_varid, "grid_mapping", "projection_utm") )
        call check(  nf90_put_att(ncid, val_varid, "coordinates", "lon lat") )
        if (scale_factor.ne.1.) call check(  nf90_put_att(ncid, val_varid, "scale_factor", scale_factor) )
        
        !Close the definitions
        call check( nf90_enddef(ncid) )

    endif

    
    if (use_single_time_loop_flag) then
        !Add time to the time dimension       
        call check( nf90_inq_varid(ncid, "time", time_varid) )
        !call check( nf90_inquire_dimension(ncid, time_dimid, temp_name, n_dims(3)) )
        !n_dims(3)=n_dims(3)+1
        if (save_netcdf_average_flag) then
            n_dims(3)=nt
        else
            n_dims(3)=t_loop
        endif
        call check( nf90_put_var(ncid, time_varid, time_seconds_output_nc(1), start = (/n_dims(3)/) ) )
        !write(*,*) n_dims(3),val_dim_nc(1,time_dim_nc_index)
        !write(*,*) n_dims
        
        !Add dimension and array to existing
        call check( nf90_inq_varid(ncid, trim(name_array), val_varid) )
        if (nf90_type.eq.NF90_byte) then
            call check( nf90_put_var(ncid, val_varid, int1(val_array(:,:,1:nt)), start=(/1,1,n_dims(3)/), count=(/n_dims(1),n_dims(2),1/)) )
        elseif (nf90_type.eq.NF90_short) then
            call check( nf90_put_var(ncid, val_varid, int2(val_array(:,:,1:nt)), start=(/1,1,n_dims(3)/), count=(/n_dims(1),n_dims(2),1/)) )
        else
            call check( nf90_put_var(ncid, val_varid, val_array, start=(/1,1,n_dims(3)/), count=(/n_dims(1),n_dims(2),1/)) )
        endif
        
    else
        
        !Write the whole variable to file. Default is float
        if (nf90_type.eq.NF90_byte) then
            call check( nf90_put_var(ncid, val_varid, int1(val_array(:,:,1:nt))) )
        elseif (nf90_type.eq.NF90_short) then
            call check( nf90_put_var(ncid, val_varid, int2(val_array(:,:,1:nt))) )
        else
            call check( nf90_put_var(ncid, val_varid, val_array(:,:,1:nt)) )
        endif

    endif
    
    call check( nf90_close(ncid) )
    

    
    end subroutine uEMEP_save_netcdf_file
    
    
    subroutine check(status)
    
    use netcdf
    implicit none
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      write(*,*) 'Stopping due to netcdf error: '//trim(nf90_strerror(status))
      stop 
    end if
    
    end subroutine check
    
    
! ######################################################################	
	FUNCTION DIRECTION(UD,VD)
!	CALCULATES THE WIND DIRECTION
    !Taken from NBLM1
	IMPLICIT NONE
	REAL DIRECTION,UD,VD,PI
	PI=180./3.14159
	DIRECTION=0.
	IF (UD.GT.0.AND.VD.GE.0) DIRECTION=270.-ATAN(ABS(VD/UD))*PI
	IF (UD.LE.0.AND.VD.GT.0) DIRECTION=180.-ATAN(ABS(UD/VD))*PI
	IF (UD.LT.0.AND.VD.LE.0) DIRECTION=90.-ATAN(ABS(VD/UD))*PI
	IF (UD.GE.0.AND.VD.LT.0) DIRECTION=360.-ATAN(ABS(UD/VD))*PI
	END FUNCTION DIRECTION
! ######################################################################