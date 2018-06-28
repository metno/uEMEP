!Saves data in netcdf format
    
    subroutine uEMEP_save_netcdf_control
    
    use uEMEP_definitions

    implicit none
    
    integer i,j,k,t
    integer i_comp,i_file,i_meteo
    character(256) temp_name,unit_str,title_str,title_str_rec,var_name_temp, temp_date_str,station_name_str,temp_name_rec
    logical create_file,create_file_rec
    integer i_source,i_subsource
    integer :: EMEP_subsource=1
    real :: valid_min=0.
    real, allocatable :: temp_subgrid(:,:,:)
    real, allocatable :: aqi_subgrid(:,:,:)
    integer ii,jj
    logical :: save_compounds=.true.,save_source_contributions=.true.,save_wind_vectors=.true.,save_other_meteo=.true.
    logical :: save_emep_source_contributions=.false.,save_emep_original=.true.,save_emissions=.false.,save_for_chemistry=.false.
    logical :: save_aqi=.true.
    real aqi_limits(5),max_aqi
    
    if (.not.allocated(temp_subgrid)) allocate(temp_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index)))
    if (.not.allocated(aqi_subgrid).and.save_aqi) allocate(aqi_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index)))
    
    !Save subgrid calculations
    valid_min=0.   
    unit_str="ug/m3"
    
    !Save the data
    i_file=subgrid_total_file_index(allsource_index)
    !temp_name=trim(pathname_grid(i_file))//trim(filename_grid(i_file))//'_'//trim(file_tag)//'.nc'
    if (len(config_date_str).gt.0) then
        temp_date_str='_'//trim(config_date_str)
    else
        temp_date_str=''
    endif
    
    if (use_multiple_receptor_grids_flag) then
        station_name_str=trim(name_receptor_in(g_loop,1))//'_';
    else
        station_name_str=''
    endif
    
    
    temp_name=trim(pathname_grid(i_file))//trim(station_name_str)//'uEMEP_output'//'_'//trim(var_name_nc(conc_nc_index,compound_index,allsource_index))//'_'//trim(file_tag)//trim(temp_date_str)//'.nc'
    temp_name_rec=trim(pathname_grid(i_file))//'uEMEP_receptor'//'_'//trim(var_name_nc(conc_nc_index,compound_index,allsource_index))//'_'//trim(file_tag)//trim(temp_date_str)//'.nc'

    write(unit_logfile,'(A)') ''
    write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Saving netcdf data (uEMEP_save_netcdf_control)'
	write(unit_logfile,'(A)') '================================================================'
    
    if (save_netcdf_receptor_flag.and.n_valid_receptor.eq.0) then
        if (i_comp.eq.1.and.t_loop.eq.start_time_loop_index) then
            write(unit_logfile,'(a)')'No receptor positions available. Will not save receptor data.'
        endif
    endif
        
    !Save the final result of the subgrid calculation
    if (save_compounds) then
    do i_comp=1,n_compound_loop
        
        if (i_comp.eq.1.and.t_loop.eq.start_time_loop_index) then
            create_file=.true.
            title_str='uEMEP_concentration_'//trim(file_tag)//temp_date_str
            write(unit_logfile,'(a)')'Writing to: '//trim(temp_name)
        else
            create_file=.false.
        endif

        if (i_comp.eq.1.and.t_loop.eq.start_time_loop_index.and.g_loop.eq.start_grid_loop_index) then
            create_file_rec=.true.
            title_str_rec='uEMEP_receptor_'//trim(file_tag)//temp_date_str
            write(unit_logfile,'(a)')'Writing to: '//trim(temp_name_rec)
        else
            create_file_rec=.false.
        endif
        
        var_name_temp=trim(var_name_nc(conc_nc_index,compound_loop_index(i_comp),allsource_index))//'_concentration'
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a)')'Writing netcdf array variable:    '//trim(var_name_temp)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
            ,comp_subgrid(:,:,:,compound_loop_index(i_comp)),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
            ,unit_str,title_str,create_file,valid_min)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
            write(unit_logfile,'(a)')'Writing netcdf receptor variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
            ,comp_subgrid(:,:,:,compound_loop_index(i_comp)),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
            ,unit_str,title_str_rec,create_file_rec,valid_min &
            ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
            ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
            ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor)          
        endif
        
    enddo
    endif
    
    create_file=.false.
    create_file_rec=.false.
    i_comp=1
    
    !Save the different local source contributions, not the total though
    if (save_source_contributions) then
    do i_source=1,n_source_index
        !if (calculate_source(i_source).or.i_source.eq.allsource_index) then
        if (calculate_source(i_source).and.i_source.ne.allsource_index) then
    
            i_file=subgrid_local_file_index(i_source)
            var_name_temp=trim(var_name_nc(conc_nc_index,compound_loop_index(i_comp),allsource_index))//'_'//trim(filename_grid(i_file))
            if (save_netcdf_file_flag) then
                write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
                call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,subgrid(:,:,:,local_subgrid_index,i_source,emep_subsource),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                    ,unit_str,title_str,create_file,valid_min)
            endif
            if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a)')'Writing netcdf receptor variable: '//trim(var_name_temp)
                call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,subgrid(:,:,:,local_subgrid_index,i_source,emep_subsource),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                    ,unit_str,title_str_rec,create_file_rec,valid_min &
                    ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor)          
            endif
            
        endif
            
            !Save the total nonlocal part from EMEP
            if (i_source.eq.allsource_index) then
               
                i_file=emep_subgrid_nonlocal_file_index(i_source)
                var_name_temp=trim(var_name_nc(conc_nc_index,compound_loop_index(i_comp),allsource_index))//'_'//trim(filename_grid(i_file))
                if (save_netcdf_file_flag) then
                    write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,subgrid(:,:,:,emep_nonlocal_subgrid_index,i_source,emep_subsource),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str,create_file,valid_min)
                endif
                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                    write(unit_logfile,'(a)')'Writing netcdf receptor variable: '//trim(var_name_temp)
                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,subgrid(:,:,:,emep_nonlocal_subgrid_index,i_source,emep_subsource),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str_rec,create_file_rec,valid_min &
                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor)          
                endif
                
            endif
    enddo
    endif

    !Save the emissions interpolated to the target grid
    if (save_emissions) then
    do i_source=1,n_source_index
        if (calculate_source(i_source)) then
    
                i_file=emission_file_index(i_source)
                var_name_temp=trim(var_name_nc(conc_nc_index,compound_loop_index(i_comp),allsource_index))//'_'//trim(filename_grid(i_file))
                
                !Calculate the emissions in the target grid
                temp_subgrid=0.
                do j=1,subgrid_dim(y_dim_index)
                do i=1,subgrid_dim(x_dim_index)
                
                    
                    ii=crossreference_target_to_emission_subgrid(i,j,x_dim_index,i_source)
                    jj=crossreference_target_to_emission_subgrid(i,j,y_dim_index,i_source)
            
                    temp_subgrid(i,j,:)=sum(emission_subgrid(ii,jj,:,i_source,:),2)
                    !Subgrid emissions, if relevant, are in units ug/sec/subgrid. Convert to g/s. Acount for the difference in subgrid sizes here
                    temp_subgrid(i,j,:)=1.0e-6*temp_subgrid(i,j,:)*(subgrid_delta(y_dim_index)*subgrid_delta(x_dim_index)) &
                        /(emission_subgrid_delta(y_dim_index,i_source)*emission_subgrid_delta(x_dim_index,i_source))
 
                enddo
                enddo
                
                unit_str='g/s'
                if (save_netcdf_file_flag) then
                    write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str,create_file,valid_min)
                endif
                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                    write(unit_logfile,'(a)')'Writing netcdf receptor variable: '//trim(var_name_temp)
                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str_rec,create_file_rec,valid_min &
                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor)          
                endif
    
        endif
    enddo
    endif
    
    !Save the EMEP data interpolated to the subgrid
    if (save_EMEP_source_contributions) then
    do i_source=1,n_source_index
        if (calculate_source(i_source).or.i_source.eq.allsource_index) then
            !do subsource_index=1,n_subsource(source_index)
            !temp_name=trim(pathname_grid(i_file))//trim(filename_grid(i_file))//'_'//trim(file_tag)//'.asc'
            !write(unit_logfile,'(a)')'Writing to: '//trim(temp_name)
            !call write_esri_ascii_3d_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),subgrid_delta(1),subgrid(:,:,:,emep_subgrid_index,source_index,emep_subsource),x_subgrid,y_subgrid)
            
            !Only save the allsource value here
            if (i_source.eq.allsource_index) then
                i_file=emep_subgrid_file_index(i_source)
                var_name_temp=trim(var_name_nc(conc_nc_index,compound_loop_index(i_comp),allsource_index))//'_'//trim(filename_grid(i_file))
                if (save_netcdf_file_flag) then
                    write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,subgrid(:,:,:,emep_subgrid_index,i_source,emep_subsource),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str,create_file,valid_min)
                endif
                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                    write(unit_logfile,'(a)')'Writing netcdf receptor variable: '//trim(var_name_temp)
                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,subgrid(:,:,:,emep_subgrid_index,i_source,emep_subsource),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str_rec,create_file_rec,valid_min &
                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor)          
                endif
                
            endif

                i_file=emep_subgrid_local_file_index(i_source)
                var_name_temp=trim(var_name_nc(conc_nc_index,compound_loop_index(i_comp),allsource_index))//'_'//trim(filename_grid(i_file))
                if (save_netcdf_file_flag) then
                    write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,subgrid(:,:,:,emep_local_subgrid_index,i_source,emep_subsource),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str,create_file,valid_min)
                endif
                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                    write(unit_logfile,'(a)')'Writing netcdf receptor variable: '//trim(var_name_temp)
                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,subgrid(:,:,:,emep_local_subgrid_index,i_source,emep_subsource),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str_rec,create_file_rec,valid_min &
                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor)          
                endif
                
            
        endif
    enddo
    endif
    
    !Save the other interpolated EMEP compounds used for nox chemistry as well
    if (save_for_chemistry) then
    do i_comp=1,n_compound_loop
        i_file=emep_subgrid_file_index(allsource_index)
        var_name_temp=trim(var_name_nc(conc_nc_index,compound_loop_index(i_comp),allsource_index))//'_interpolated'//'_'//trim(filename_grid(i_file))
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,comp_EMEP_subgrid(:,:,:,compound_loop_index(i_comp)),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str,create_file,valid_min)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
            write(unit_logfile,'(a)')'Writing netcdf receptor variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,comp_EMEP_subgrid(:,:,:,compound_loop_index(i_comp)),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor)          
        endif
    enddo
    endif
    
    !Save the original EMEP compounds
    if (save_emep_original) then
    do i_comp=1,n_compound_loop
        i_file=emep_subgrid_file_index(allsource_index)
        var_name_temp=trim(var_name_nc(conc_nc_index,compound_loop_index(i_comp),allsource_index))//'_original_EMEP_concentration'
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,orig_EMEP_subgrid(:,:,:,compound_loop_index(i_comp)),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str,create_file,valid_min)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
            write(unit_logfile,'(a)')'Writing netcdf receptor variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,orig_EMEP_subgrid(:,:,:,compound_loop_index(i_comp)),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor)          
        endif
    enddo
    endif
    
    !Save AQI
    if (save_aqi) then
        !Temporary creation of AQI for testing
        aqi_limits(1)=0.;aqi_limits(2)=100.;aqi_limits(3)=200.;aqi_limits(4)=400.;aqi_limits(5)=800.
        aqi_subgrid=5.
        max_aqi=0.
        do t=1,subgrid_dim(t_dim_index)
        do j=1,subgrid_dim(y_dim_index)
        do i=1,subgrid_dim(x_dim_index)
            do k=1,4
                if (comp_subgrid(i,j,t,no2_index).ge.aqi_limits(k).and.comp_subgrid(i,j,t,no2_index).lt.aqi_limits(k+1)) then
                    aqi_subgrid(i,j,t)=k+(comp_subgrid(i,j,t,no2_index)-aqi_limits(k))/(aqi_limits(k+1)-aqi_limits(k))
                endif
            enddo
            aqi_subgrid(i,j,t)=min(aqi_subgrid(i,j,t),4.99)
            if (aqi_subgrid(i,j,t).gt.max_aqi) max_aqi=aqi_subgrid(i,j,t)
            !write(*,*)  aqi_subgrid(i,j,t),comp_subgrid(i,j,t,no2_index),comp_subgrid(i,j,t,nox_index)
        enddo
        enddo
        enddo
        write(unit_logfile,*)  'MAX AQI: ',max_aqi
       
        var_name_temp='AQI'
        unit_str='1'
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,aqi_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str,create_file,valid_min)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
            write(unit_logfile,'(a)')'Writing netcdf receptor variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,aqi_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor)          
        endif
    endif

    !Save weighted travel time
    if (save_for_chemistry) then
        var_name_temp='Weighted_travel_time'
        unit_str='seconds'
        if (save_netcdf_file_flag) then
            write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,traveltime_subgrid(:,:,:,1),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str,create_file,valid_min)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
            write(unit_logfile,'(a)')'Writing netcdf receptor variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,traveltime_subgrid(:,:,:,1),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                ,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor)          
        endif
    endif
 
    !Save the meteo interpolated to the target grid
    valid_min=-1.e24   
    
    if (save_wind_vectors) then
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
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
            write(unit_logfile,'(a)')'Writing netcdf receptor variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor)          
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
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
            write(unit_logfile,'(a)')'Writing netcdf receptor variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor)          
        endif
        !The same for all--------------------
    endif
    
    if (save_for_chemistry) then
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
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
            write(unit_logfile,'(a)')'Writing netcdf receptor variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor)          
        endif
        !The same for all--------------------
    endif

    if (save_other_meteo) then
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
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
            write(unit_logfile,'(a)')'Writing netcdf receptor variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor)          
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
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
            write(unit_logfile,'(a)')'Writing netcdf receptor variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor)          
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
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
            write(unit_logfile,'(a)')'Writing netcdf receptor variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor)          
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
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
            write(unit_logfile,'(a)')'Writing netcdf receptor variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor)          
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
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
            write(unit_logfile,'(a)')'Writing netcdf receptor variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor)          
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
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
            write(unit_logfile,'(a)')'Writing netcdf receptor variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor)          
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
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
            write(unit_logfile,'(a)')'Writing netcdf receptor variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor)          
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
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
            write(unit_logfile,'(a)')'Writing netcdf receptor variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor)          
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
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min)
        endif
        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
            write(unit_logfile,'(a)')'Writing netcdf receptor variable: '//trim(var_name_temp)
            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor)          
        endif
        !The same for all--------------------

        endif
    
    endif
    
    end subroutine uEMEP_save_netcdf_control
    
    
    subroutine uEMEP_save_netcdf_file(unit_logfile_in,filename_netcdf,nx,ny,nt,val_array,x_array,y_array,lon_array,lat_array,name_array,unit_array,title_str,create_file,valid_min)
    
    use uEMEP_definitions
    use netcdf
    
    implicit none
    
    character(256) filename_netcdf,name_array,unit_array,title_str,temp_name,temp_name3(3)
    integer unit_logfile_in
    integer nx,ny,nt
    real val_array(nx,ny,nt)!,val_array_temp(nx,ny,nt)
    real x_array(nx,ny)
    real y_array(nx,ny)
    real lon_array(nx,ny)
    real lat_array(nx,ny) !,lat_array_temp(nx,ny)
    !real time_array(nt)
    real x_vector(nx)
    real y_vector(ny)
    logical create_file
    real valid_min
    
    integer ncid
    integer y_dimid,x_dimid,lat_dimid,lon_dimid,val_dimid,time_dimid
    integer y_varid,x_varid,lat_varid,lon_varid,val_varid,time_varid,proj_varid
    integer dimids3(3),dimids2(2),chunks3(3)
    integer n_dims(3)
    integer status
    

    !integer j
    
    !Assumes x and y are the dimensions   
    x_vector=x_array(:,1)
    y_vector=y_array(1,:)
    
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
        call check(  nf90_def_var(ncid, "projection_utm", NF90_int, proj_varid) )
        call check(  nf90_put_att(ncid, proj_varid, "semi_major_axis", 6378137.0 ) )
        call check(  nf90_put_att(ncid, proj_varid, "inverse_flattening", 298.257222101 ) )

        call check(  nf90_put_att(ncid, proj_varid, "grid_mapping_name", "transverse_mercator" ) )
        call check(  nf90_put_att(ncid, proj_varid, "scale_factor_at_central_meridian", 0.9996 ) )
        call check(  nf90_put_att(ncid, proj_varid, "latitude_of_projection_origin", 0 ) )
        call check(  nf90_put_att(ncid, proj_varid, "false_easting", 500000. ) )
        call check(  nf90_put_att(ncid, proj_varid, "longitude_of_central_meridian", utm_lon0 ) )
        !call check(  nf90_put_att(ncid, proj_varid, "semi_major_axis", 6378140.0 ) )
        !call check(  nf90_put_att(ncid, proj_varid, "semi_minor_axis", 6356750.0 ) )
  
        !Define the dimmensions
        call check(  nf90_def_dim(ncid,"time",NF90_UNLIMITED, time_dimid) )
        call check(  nf90_def_dim(ncid, "y", ny, y_dimid) )
        call check(  nf90_def_dim(ncid, "x", nx, x_dimid) )
     
        !Define the dimmension variables
        call check(  nf90_def_var(ncid, "time", NF90_DOUBLE, time_dimid, time_varid) )
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

        !Specify other dimmension attributes
        call check(  nf90_put_att(ncid, y_varid, "standard_name", "projection_y_axis") )
        call check(  nf90_put_att(ncid, x_varid, "standard_name", "projection_x_axis") )
        call check(  nf90_put_att(ncid, y_varid, "axis", "Y") )
        call check(  nf90_put_att(ncid, x_varid, "axis", "X") )
  
        !Close the definitions
        call check( nf90_enddef(ncid) )

        call check( nf90_put_var(ncid, time_varid, val_dim_nc(1:dim_length_nc(time_dim_nc_index),time_dim_nc_index)) )
        call check( nf90_put_var(ncid, y_varid, y_vector) )
        call check( nf90_put_var(ncid, x_varid, x_vector) )
        call check( nf90_put_var(ncid, lat_varid, lat_array) )
        call check( nf90_put_var(ncid, lon_varid, lon_array) )

        call check( nf90_close(ncid) )
    
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
        call check( nf90_def_var(ncid, trim(name_array), NF90_REAL, dimids3, val_varid) )
        ! gzip level 3 compression and shuffling
        ! optional _FillValue for values which never have been written, unpacked value
        call check( nf90_def_var_chunking(ncid, val_varid, NF90_CHUNKED, chunks3) ) !New
        call check( nf90_def_var_deflate(ncid, val_varid, 1, 1, 3) ) !New
        call check( nf90_put_att(ncid, val_varid, "units", trim(unit_array)) )
    
        !Specify other variable attributes
        call check(  nf90_put_att(ncid, val_varid, "missing_value", NODATA_value ) ) !New
        call check(  nf90_put_att(ncid, val_varid, "grid_mapping", "projection_utm") )
        call check(  nf90_put_att(ncid, val_varid, "coordinates", "lon lat") )
        call check(  nf90_put_att(ncid, val_varid, "valid_min", valid_min) )
        
        !Close the definitions
        call check( nf90_enddef(ncid) )

    endif

    
    if (use_single_time_loop_flag) then
        !Add time to the time dimension       
        call check( nf90_inq_varid(ncid, "time", time_varid) )
        !call check( nf90_inquire_dimension(ncid, time_dimid, temp_name, n_dims(3)) )
        !n_dims(3)=n_dims(3)+1
        n_dims(3)=t_loop
        call check( nf90_put_var(ncid, time_varid, val_dim_nc(1,time_dim_nc_index), start = (/n_dims(3)/) ) )
        !write(*,*) n_dims(3),val_dim_nc(1,time_dim_nc_index)
        !write(*,*) n_dims
        
        !Add dimension and array to existing
        call check( nf90_inq_varid(ncid, trim(name_array), val_varid) )       
        call check( nf90_put_var(ncid, val_varid, val_array, start=(/1,1,n_dims(3)/), count=(/n_dims(1),n_dims(2),1/)) )
    else
        
        !Write the variable to file
        call check( nf90_put_var(ncid, val_varid, val_array) )
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