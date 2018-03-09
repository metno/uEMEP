!Saves data in netcdf format
    
    subroutine uEMEP_save_netcdf_control
    
    use uEMEP_definitions

    implicit none
    
    integer i_comp,i_file,i_meteo
    character(256) temp_name,unit_str,title_str,var_name_temp, temp_date_str,station_name_str
    logical create_file
    integer i_source,i_subsource
    integer :: EMEP_subsource=1
    real :: valid_min=0.
    real, allocatable :: temp_subgrid(:,:,:)
    integer ii,jj
    
    if (.not.allocated(temp_subgrid)) allocate(temp_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index)))
    
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

    write(unit_logfile,'(A)') ''
    write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Saving netcdf data (uEMEP_save_netcdf_control)'
	write(unit_logfile,'(A)') '================================================================'
    !Save the final result of the subgrid calculation
    do i_comp=1,n_compound_loop
        
        if (i_comp.eq.1.and.t_loop.eq.start_time_loop_index) then
            create_file=.true.
            title_str='uEMEP_concentration_'//trim(file_tag)//temp_date_str
            write(unit_logfile,'(a)')'Writing to: '//trim(temp_name)
        else
            create_file=.false.
        endif
        
        var_name_temp=trim(var_name_nc(conc_nc_index,compound_loop_index(i_comp),allsource_index))//'_concentration'
        write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
        call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
            ,comp_subgrid(:,:,:,compound_loop_index(i_comp)),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
            ,unit_str,title_str,create_file,valid_min)
    enddo

    create_file=.false.
    i_comp=1
    
    !Save the different local source contributions
    do i_source=1,n_source_index
        if (calculate_source(i_source).or.i_source.eq.allsource_index) then
    
                i_file=subgrid_local_file_index(i_source)
                var_name_temp=trim(var_name_nc(conc_nc_index,compound_loop_index(i_comp),i_source))//'_'//trim(filename_grid(i_file))
                write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
                call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,subgrid(:,:,:,local_subgrid_index,i_source,emep_subsource),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                    ,unit_str,title_str,create_file,valid_min)
    
        endif
    enddo

    !Save the emissions interpolated to the target grid
    do i_source=1,n_source_index
        if (calculate_source(i_source)) then
    
                i_file=emission_file_index(i_source)
                var_name_temp=trim(var_name_nc(conc_nc_index,compound_loop_index(i_comp),i_source))//'_'//trim(filename_grid(i_file))
                write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
                
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
                call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                    ,unit_str,title_str,create_file,valid_min)
    
        endif
    enddo

    !Save the EMEP data interpolated to the subgrid
    do i_source=1,n_source_index
        if (calculate_source(i_source).or.i_source.eq.allsource_index) then
            !do subsource_index=1,n_subsource(source_index)
            !temp_name=trim(pathname_grid(i_file))//trim(filename_grid(i_file))//'_'//trim(file_tag)//'.asc'
            !write(unit_logfile,'(a)')'Writing to: '//trim(temp_name)
            !call write_esri_ascii_3d_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),subgrid_delta(1),subgrid(:,:,:,emep_subgrid_index,source_index,emep_subsource),x_subgrid,y_subgrid)
            
            !Only save the allsource value here
            if (i_source.eq.allsource_index) then
                i_file=emep_subgrid_file_index(i_source)
                var_name_temp=trim(var_name_nc(conc_nc_index,compound_loop_index(i_comp),i_source))//'_'//trim(filename_grid(i_file))
                write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
                call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,subgrid(:,:,:,emep_subgrid_index,i_source,emep_subsource),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                    ,unit_str,title_str,create_file,valid_min)
            endif
                
                i_file=emep_subgrid_nonlocal_file_index(i_source)
                var_name_temp=trim(var_name_nc(conc_nc_index,compound_loop_index(i_comp),i_source))//'_'//trim(filename_grid(i_file))
                write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
                call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,subgrid(:,:,:,emep_nonlocal_subgrid_index,i_source,emep_subsource),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                    ,unit_str,title_str,create_file,valid_min)
                
                i_file=emep_subgrid_local_file_index(i_source)
                var_name_temp=trim(var_name_nc(conc_nc_index,compound_loop_index(i_comp),i_source))//'_'//trim(filename_grid(i_file))
                write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
                call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,subgrid(:,:,:,emep_local_subgrid_index,i_source,emep_subsource),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                    ,unit_str,title_str,create_file,valid_min)
                
                
            
        endif
    enddo

    
    !Save the other interpolated EMEP compounds used for nox chemistry as well
    do i_comp=1,n_compound_loop
        i_file=emep_subgrid_file_index(allsource_index)
        var_name_temp=trim(var_name_nc(conc_nc_index,compound_loop_index(i_comp),allsource_index))//'_interpolated'//'_'//trim(filename_grid(i_file))
        write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
        call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
            ,comp_EMEP_subgrid(:,:,:,compound_loop_index(i_comp)),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
            ,unit_str,title_str,create_file,valid_min)
    enddo

    !Save the original EMEP compounds
    do i_comp=1,n_compound_loop
        i_file=emep_subgrid_file_index(allsource_index)
        var_name_temp=trim(var_name_nc(conc_nc_index,compound_loop_index(i_comp),allsource_index))//'_original_EMEP_concentration'
        write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
        call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
            ,orig_EMEP_subgrid(:,:,:,compound_loop_index(i_comp)),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
            ,unit_str,title_str,create_file,valid_min)
    enddo
    
    !Save weighted travel time
        var_name_temp='Weighted_travel_time'
        unit_str='seconds'
        write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
        call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
            ,traveltime_subgrid(:,:,:,1),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
            ,unit_str,title_str,create_file,valid_min)

 
    !Save the meteo interpolated to the target grid
    valid_min=-1.e24   
       
    i_file=subgrid_ugrid_file_index
    i_meteo=ugrid_subgrid_index
    unit_str="m/s"
    !The same for all--------------------
    var_name_temp=trim(filename_grid(i_file))
    write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
    temp_subgrid=0.
    do j=1,subgrid_dim(y_dim_index)
    do i=1,subgrid_dim(x_dim_index)
        ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
    enddo
    enddo
    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
            ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min)
    !The same for all--------------------

    i_file=subgrid_vgrid_file_index
    i_meteo=vgrid_subgrid_index
    unit_str="m/s"
    !The same for all--------------------
    var_name_temp=trim(filename_grid(i_file))
    write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
    temp_subgrid=0.
    do j=1,subgrid_dim(y_dim_index)
    do i=1,subgrid_dim(x_dim_index)
        ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
    enddo
    enddo
    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
            ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min)
    !The same for all--------------------

    i_file=subgrid_hmix_file_index
    i_meteo=hmix_subgrid_index
    unit_str="m"
    !The same for all--------------------
    var_name_temp=trim(filename_grid(i_file))
    write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
    temp_subgrid=0.
    do j=1,subgrid_dim(y_dim_index)
    do i=1,subgrid_dim(x_dim_index)
        ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
    enddo
    enddo
    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
            ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min)
    !The same for all--------------------
    
    i_file=subgrid_kz_file_index
    i_meteo=kz_subgrid_index
    unit_str="m2/s"
    !The same for all--------------------
    var_name_temp=trim(filename_grid(i_file))
    write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
    temp_subgrid=0.
    do j=1,subgrid_dim(y_dim_index)
    do i=1,subgrid_dim(x_dim_index)
        ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
    enddo
    enddo
    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
            ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min)
    !The same for all--------------------
 
    i_file=subgrid_FFgrid_file_index
    i_meteo=FFgrid_subgrid_index
    unit_str="m/s"
    !The same for all--------------------
    var_name_temp=trim(filename_grid(i_file))
    write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
    temp_subgrid=0.
    do j=1,subgrid_dim(y_dim_index)
    do i=1,subgrid_dim(x_dim_index)
        ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
    enddo
    enddo
    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
            ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min)
    !The same for all--------------------

    i_file=subgrid_FF10_file_index
    i_meteo=FF10_subgrid_index
    unit_str="m/s"
    !The same for all--------------------
    var_name_temp=trim(filename_grid(i_file))
    write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
    temp_subgrid=0.
    do j=1,subgrid_dim(y_dim_index)
    do i=1,subgrid_dim(x_dim_index)
        ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
    enddo
    enddo
    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
            ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min)
    !The same for all--------------------

    i_file=subgrid_J_file_index
    i_meteo=J_subgrid_index
    unit_str="1/s"
    !The same for all--------------------
    var_name_temp=trim(filename_grid(i_file))
    write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
    temp_subgrid=0.
    do j=1,subgrid_dim(y_dim_index)
    do i=1,subgrid_dim(x_dim_index)
        ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
    enddo
    enddo
    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
            ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min)
    !The same for all--------------------

    i_file=subgrid_invL_file_index
    i_meteo=invL_subgrid_index
    unit_str="1/m"
    !The same for all--------------------
    var_name_temp=trim(filename_grid(i_file))
    write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
    temp_subgrid=0.
    do j=1,subgrid_dim(y_dim_index)
    do i=1,subgrid_dim(x_dim_index)
        ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
    enddo
    enddo
    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
            ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min)
    !The same for all--------------------

    i_file=subgrid_logz0_file_index
    i_meteo=logz0_subgrid_index
    unit_str="log(m)"
    !The same for all--------------------
    var_name_temp=trim(filename_grid(i_file))
    write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
    temp_subgrid=0.
    do j=1,subgrid_dim(y_dim_index)
    do i=1,subgrid_dim(x_dim_index)
        ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
    enddo
    enddo
    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
            ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min)
    !The same for all--------------------

    i_file=subgrid_ustar_file_index
    i_meteo=ustar_subgrid_index
    unit_str="m/s"
    !The same for all--------------------
    var_name_temp=trim(filename_grid(i_file))
    write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
    temp_subgrid=0.
    do j=1,subgrid_dim(y_dim_index)
    do i=1,subgrid_dim(x_dim_index)
        ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
    enddo
    enddo
    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
            ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min)
    !The same for all--------------------
            
    i_file=subgrid_invFFgrid_file_index
    i_meteo=inv_FFgrid_subgrid_index
    unit_str="s/m"
    !The same for all--------------------
    var_name_temp=trim(filename_grid(i_file))
    write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
    temp_subgrid=0.
    do j=1,subgrid_dim(y_dim_index)
    do i=1,subgrid_dim(x_dim_index)
        ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
    enddo
    enddo
    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
            ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min)
    !The same for all--------------------
            
    i_file=subgrid_invFF10_file_index
    i_meteo=inv_FF10_subgrid_index
    unit_str="s/m"
    !The same for all--------------------
    var_name_temp=trim(filename_grid(i_file))
    write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
    temp_subgrid=0.
    do j=1,subgrid_dim(y_dim_index)
    do i=1,subgrid_dim(x_dim_index)
        ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
    enddo
    enddo
    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
            ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min)
    !The same for all--------------------

    end subroutine uEMEP_save_netcdf_control
    
    
    subroutine uEMEP_save_netcdf_file(unit_logfile_in,filename_netcdf,nx,ny,nt,val_array,x_array,y_array,lon_array,lat_array,name_array,unit_array,title_str,create_file,valid_min)
    
    use uEMEP_definitions
    use netcdf
    
    implicit none
    
    character(256) filename_netcdf,name_array,unit_array,title_str,temp_name,temp_name3(3)
    integer unit_logfile_in
    integer nx,ny,nt
    real val_array(nx,ny,nt),val_array_temp(nx,ny,nt)
    real x_array(nx,ny)
    real y_array(nx,ny)
    real lon_array(nx,ny)
    real lat_array(nx,ny),lat_array_temp(nx,ny)
    real time_array(nt)
    real x_vector(nx)
    real y_vector(ny)
    logical create_file
    real valid_min
    
    integer ncid
    integer y_dimid,x_dimid,lat_dimid,lon_dimid,val_dimid,time_dimid
    integer y_varid,x_varid,lat_varid,lon_varid,val_varid,time_varid,proj_varid
    integer dimids3(3),dimids2(2)
    integer n_dims(3)
    integer status
    

    !integer j
    
    !Assumes x and y are the dimensions   
    x_vector=x_array(:,1)
    y_vector=y_array(1,:)
    
    if (create_file) then
        !Create a netcdf file
        call check(  nf90_create(filename_netcdf, nf90_clobber, ncid) )

        !Specify global attributes
        call check(  nf90_put_att(ncid, nf90_global, "Conventions", "CF-1.4" ) )
        call check(  nf90_put_att(ncid, nf90_global, "title", trim(title_str)) )
        call check(  nf90_put_att(ncid, nf90_global, "Model", "uEMEP" ) )
        call check(  nf90_put_att(ncid, nf90_global, "missing_value", NODATA_value ) )
        
    
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
    call check( nf90_redef(ncid) )
    
    !Get the dimensions id from the existing file
    call check( nf90_inq_dimid(ncid,"time",time_dimid) )
    call check( nf90_inq_dimid(ncid, "y", y_dimid) )
    call check( nf90_inq_dimid(ncid, "x", x_dimid) )
    dimids3 = (/ x_dimid, y_dimid, time_dimid /)
    call check( nf90_inquire_dimension(ncid, dimids3(1), temp_name, n_dims(1)) )
    call check( nf90_inquire_dimension(ncid, dimids3(2), temp_name, n_dims(2)) )
    call check( nf90_inquire_dimension(ncid, dimids3(3), temp_name, n_dims(3)) )


    status=nf90_inq_varid(ncid, trim(name_array), val_varid)
    if (status.ne.nf90_NoErr) then
        !if the variable does not exist then create a new one
        !write(*,*) 'Creating new: ',trim(name_array)
        call check( nf90_def_var(ncid, trim(name_array), NF90_REAL, dimids3, val_varid) )
        call check( nf90_put_att(ncid, val_varid, "units", trim(unit_array)) )
    
        !Specify other variable attributes
        call check(  nf90_put_att(ncid, val_varid, "grid_mapping", "projection_utm") )
        call check(  nf90_put_att(ncid, val_varid, "coordinates", "lon lat") )
        call check(  nf90_put_att(ncid, val_varid, "valid_min", valid_min) )
    endif

    !Close the definitions
    call check( nf90_enddef(ncid) )
    
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