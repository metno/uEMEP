!uEMEP_save_emission_netcdf.f90
    !This routine saves the various emission sources in the EMEP grid
    !It first reads in an example EMEP file (z0 file) to get projections and x,y grid dimensions
    !Then extracts the emission data from its original form
    !It should run independently and stop when it is finished
    
    subroutine uEMEP_calculate_emissions_for_EMEP
    
    use uEMEP_definitions
    
    implicit none
    
    !double precision :: EMEP_projection_attributes(10)
    !real, allocatable :: EMEP_emissions_grid(:,:,:,:,:) !x,y,t,source,pollutant
    integer a_start(6),date_array(6),a_start_emission(6)
    character(256) format_temp
    double precision date_num_temp,date_num_start,date_num_start_emission
    integer t,i_source,i,j
    
    !Functions
    double precision date_to_number

    write(unit_logfile,'(A)') ''
    write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Saving emission data for EMEP (uEMEP_calculate_emissions_for_EMEP)'
	write(unit_logfile,'(A)') '================================================================'

 
    !Read in example EMEP file for dimensions and projection information
    !Or just defne them here without reading

    if (trim(save_emissions_for_EMEP_projection).eq.'lambert') then
        projection_type=LCC_projection_index
        !grid_mapping_name = "lambert_conformal_conic";
        EMEP_projection_attributes(1:2) = 63.0
        EMEP_projection_attributes(3) = 15.0
        EMEP_projection_attributes(4) = 63.0
        EMEP_projection_attributes(5) = 6370000.0
        !float i(i=531);
        !standard_name = "projection_x_coordinate";
        !long_name = "x-coordinate in Cartesian system";
        !units = "m";
        !axis = "X";
        if (trim(save_emissions_for_EMEP_region).eq.'NO') then
            emission_subgrid_min(x_dim_index,:)=-6.498834E+05 ! 6.751166E+05
            emission_subgrid_delta(x_dim_index,:)=2500.
            emission_subgrid_dim(x_dim_index,:)=531
            emission_max_subgrid_dim(x_dim_index)=531
            !Move minium to edge, for consistency with normal subgrid definition
            emission_subgrid_min(x_dim_index,:)=emission_subgrid_min(x_dim_index,:)-emission_subgrid_delta(x_dim_index,:)/2.

            !float j(j=671);
            !standard_name = "projection_y_coordinate";
            !long_name = "y-coordinate in Cartesian system";
            !units = "m";
            !axis = "Y";
            emission_subgrid_min(y_dim_index,:)=-6.567275E+05 ! 1.018272E+06
            emission_subgrid_delta(y_dim_index,:)=2500.
            emission_subgrid_dim(y_dim_index,:)=671
            emission_max_subgrid_dim(y_dim_index)=671
            !Move minium to edge, for consistency with normal subgrid definition
            emission_subgrid_min(y_dim_index,:)=emission_subgrid_min(y_dim_index,:)-emission_subgrid_delta(y_dim_index,:)/2.
        else
            write(unit_logfile,'(A)') 'ERROR: Emission region '//trim(save_emissions_for_EMEP_region)//' not currently defined for lambert coordinates'
            stop      
        endif
    
    elseif (trim(save_emissions_for_EMEP_projection).eq.'latlon') then
        projection_type=LL_projection_index
        write(unit_logfile,'(A,i)') 'Projection of emission grid set to '//trim(save_emissions_for_EMEP_projection),projection_type
        if (trim(save_emissions_for_EMEP_region).eq.'NL') then
            emission_subgrid_min(x_dim_index,:)=3.0 
            emission_subgrid_delta(x_dim_index,:)=0.25
            emission_subgrid_dim(x_dim_index,:)=floor((8.-3.)/0.25)+1
            emission_max_subgrid_dim(x_dim_index)=floor((8.-3.)/0.25)+1
            !Full EMEP European domain
            emission_subgrid_min(x_dim_index,:)=-30. 
            emission_subgrid_delta(x_dim_index,:)=0.25
            emission_subgrid_dim(x_dim_index,:)=floor((45.+30.)/0.25)+1
            emission_max_subgrid_dim(x_dim_index)=floor((45.+30.)/0.25)+1
            !Move minium to edge, for consistency with normal subgrid definition
            emission_subgrid_min(x_dim_index,:)=emission_subgrid_min(x_dim_index,:)-emission_subgrid_delta(x_dim_index,:)/2.

            !float j(j=671);
            !standard_name = "projection_y_coordinate";
            !long_name = "y-coordinate in Cartesian system";
            !units = "m";
            !axis = "Y";
            emission_subgrid_min(y_dim_index,:)=50. 
            emission_subgrid_delta(y_dim_index,:)=0.125
            emission_subgrid_dim(y_dim_index,:)=floor((54.-50.)/0.125)+1
            emission_max_subgrid_dim(y_dim_index)=floor((54.-50.)/0.125)+1
            !Full EMEP European domain
            emission_subgrid_min(y_dim_index,:)=30. 
            emission_subgrid_delta(y_dim_index,:)=0.125
            emission_subgrid_dim(y_dim_index,:)=floor((76.-30.)/0.125)+1
            emission_max_subgrid_dim(y_dim_index)=floor((76.-30.)/0.125)+1
            !Move minium to edge, for consistency with normal subgrid definition
            emission_subgrid_min(y_dim_index,:)=emission_subgrid_min(y_dim_index,:)-emission_subgrid_delta(y_dim_index,:)/2.
        else
	            write(unit_logfile,'(A)') 'ERROR: Emission region '//trim(save_emissions_for_EMEP_region)//' not currently defined for lat lon coordinates'
                stop      
        endif   
    
    else
        write(unit_logfile,'(A)') 'ERROR: Emission projection '//trim(save_emissions_for_EMEP_projection)//' not currently defined'
        stop        
    endif
    
    !subgrid_dim(t_dim_index)=save_emissions_end_index-save_emissions_start_index+1
    subgrid_dim(t_dim_index)=save_emissions_end_index

    dim_length_nc(x_dim_nc_index)=emission_subgrid_dim(x_dim_index,allsource_index)
    dim_length_nc(y_dim_nc_index)=emission_subgrid_dim(y_dim_index,allsource_index)
    dim_length_nc(time_dim_nc_index)=subgrid_dim(t_dim_index)
    
    write(unit_logfile,'(a,3i6)') 'EMEP emission grid dimensions: ', emission_subgrid_dim(x_dim_index,allsource_index),emission_subgrid_dim(y_dim_index,allsource_index),dim_length_nc(time_dim_nc_index)
    
    if (.not.allocated(val_dim_nc)) allocate (val_dim_nc(maxval(dim_length_nc),num_dims_nc)) !x, y, z and time dimension values
    if (.not.allocated(unit_dim_nc)) allocate (unit_dim_nc(num_dims_nc)) !x, y, z and time dimension values
    !Define the emission subgrid to correspond to the EMEP grid
    if (.not.allocated(emission_subgrid)) allocate (emission_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_source_index,n_pollutant_loop)) 
    if (.not.allocated(proxy_emission_subgrid)) allocate (proxy_emission_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),n_source_index,n_pollutant_loop)) 
    if (.not.allocated(x_emission_subgrid)) allocate (x_emission_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),n_source_index))
    if (.not.allocated(y_emission_subgrid)) allocate (y_emission_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),n_source_index))
    if (.not.allocated(lon_emission_subgrid)) allocate (lon_emission_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),n_source_index))
    if (.not.allocated(lat_emission_subgrid)) allocate (lat_emission_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),n_source_index))
    if (.not.allocated(emission_time_profile_subgrid)) allocate (emission_time_profile_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_source_index,n_pollutant_loop)) 
    if (.not.allocated(emission_properties_subgrid)) allocate (emission_properties_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),n_emission_index,n_source_index)) 
    
    !Set the timing data. Assumes no single loops with the time based on the input date and array length defined by the length of the EMEP input data
    format_temp='yyyymmddHH'
    call datestr_to_date(config_date_str,format_temp,a_start)
    !Assumes it starts at hour 1
    !a_start(4)=1
    !Do not assume it starts at hour 1 any more
    a_start(5:6)=0
    a_start_emission=a_start
    
    date_num_start=date_to_number(a_start,ref_year_EMEP)
    !Move the starting time according to the index value given (index is the number of hours)
    date_num_start_emission=date_num_start+dble(save_emissions_start_index-1)/24.
    !Set the emission_date_str to be used to name the output file as this may not be the same as the input file
    call number_to_date(date_num_start,a_start,ref_year_EMEP)
    call number_to_date(date_num_start_emission,a_start_emission,ref_year_EMEP)
    call date_to_datestr(a_start_emission,format_temp,emission_date_str)
    !Do not do this now
    emission_date_str=config_date_str
    !write(*,*) config_date_str,emission_date_str
    
    !long_name = "time at middle of period";
    unit_dim_nc(time_dim_nc_index)="days since 1900-1-1 0:0:0";
    do t=1,subgrid_dim(t_dim_index)
            date_num_temp=date_num_start+dble(t-1)/24.+dble(0.0001)/dble(24.)/dble(3600.) !Add 0.0001 of a second to avoid any rounding off errors
            call number_to_date(date_num_temp,date_array,ref_year_EMEP)
            write(*,'(6i6)') date_array(1:6)
            val_dim_nc(t,time_dim_nc_index)=date_num_temp
    enddo
    
    unit_dim_nc(x_dim_nc_index)='m'
    unit_dim_nc(y_dim_nc_index)='m'
    
    !Define emission grids
    do i_source=1,n_source_index
    if (save_emissions_for_EMEP(i_source)) then

        do j=1,emission_subgrid_dim(y_dim_index,i_source)
        do i=1,emission_subgrid_dim(x_dim_index,i_source)
        
            x_emission_subgrid(i,j,i_source)=emission_subgrid_min(x_dim_index,i_source)+emission_subgrid_delta(x_dim_index,i_source)*(i-0.5)
            y_emission_subgrid(i,j,i_source)=emission_subgrid_min(y_dim_index,i_source)+emission_subgrid_delta(y_dim_index,i_source)*(j-0.5)
        
            if (trim(save_emissions_for_EMEP_projection).eq.'lambert') then
                call lambert2lb2_uEMEP(x_emission_subgrid(i,j,i_source),y_emission_subgrid(i,j,i_source) &
                ,lon_emission_subgrid(i,j,i_source),lat_emission_subgrid(i,j,i_source),EMEP_projection_attributes)
            elseif (trim(save_emissions_for_EMEP_projection).eq.'latlon') then
                lon_emission_subgrid(i,j,i_source)=x_emission_subgrid(i,j,i_source)
                lat_emission_subgrid(i,j,i_source)=y_emission_subgrid(i,j,i_source)               
            else
                write(unit_logfile,'(A)') 'ERROR: Emission projection '//trim(save_emissions_for_EMEP_projection)//' not currently defined'
                stop
            endif
            
            !write(*,*) x_emission_subgrid(i,j,i_source),y_emission_subgrid(i,j,i_source),lon_emission_subgrid(i,j,i_source),lat_emission_subgrid(i,j,i_source)

        enddo
        enddo
    
        if (i_source.eq.industry_index) then
            !Read in industry data
            call uEMEP_read_industry_data
        
            call uEMEP_read_time_profiles
            call uEMEP_set_emission_factors
            call uEMEP_convert_proxy_to_emissions
            
        endif
        if (i_source.eq.agriculture_index) then
            !Read agriculture data
            call uEMEP_read_agriculture_rivm_data
 
            call uEMEP_read_time_profiles
            call uEMEP_set_emission_factors
            call uEMEP_convert_proxy_to_emissions         
        endif
        if (i_source.eq.traffic_index) then
            g_loop=1
            !Read inthe road data
            call uEMEP_read_roadlink_data_ascii
            if (use_NORTRIP_emission_data) then
                call uEMEP_read_roadlink_emission_data
            endif
            !Redefine the road link positions to correspond to the EMEP coordinates
        
            !Grid the data. Road link coordinates will be redefined within this routine
            call uEMEP_grid_roads
            !call uEMEP_read_time_profiles
            !call uEMEP_set_emission_factors
            !call uEMEP_convert_proxy_to_emissions
            
            !Adjust traffic emissions of NOx based on temperature
            if (use_traffic_nox_emission_temperature_dependency) then
                use_alternative_meteorology_flag=.true.
                !Set the maximum dimension to that which is necessary. Minimum is not changed as it is selected in uEMEP_save_emission_netcdf
                end_time_meteo_nc_index=start_time_meteo_nc_index+(save_emissions_end_index-1)
                call uEMEP_read_meteo_nc
                !call uEMEP_subgrid_meteo_EMEP
                !call uEMEP_crossreference_grids
                
                call uEMEP_nox_emission_temperature
            endif
            
        endif
        if (i_source.eq.heating_index) then
            !meteo_var3d_nc(i_nc,j_nc,:,t2m_nc_index)
            !Read the heating data. Emission grid coordinates will be redefined within this routine
            !Must set g_loop=1 for it to read
            g_loop=1
            call uEMEP_read_RWC_heating_data
            !Read in the temperature fields from the alternative meteorology always, since EMEP data should not exist yet
            use_alternative_meteorology_flag=.true.
            !Set the maximum dimension to that which is necessary. Minimum is not changed as it is selected in uEMEP_save_emission_netcdf
            !start_time_meteo_nc_index=start_time_meteo_nc_index+(save_emissions_start_index-1)
            end_time_meteo_nc_index=start_time_meteo_nc_index+(save_emissions_end_index-1)

            call uEMEP_read_meteo_nc
            !Need to make a cross reference here or simply skip the two based on an if statement
            call uEMEP_read_time_profiles
            call uEMEP_set_emission_factors
            call uEMEP_convert_proxy_to_emissions
            
        endif
        if (i_source.eq.shipping_index) then
            !meteo_var3d_nc(i_nc,j_nc,:,t2m_nc_index)
            !Read the heating data. Emission grid coordinates will be redefined within this routine
            !Must set g_loop=1 for it to read
            g_loop=1
            if (read_weekly_shipping_data_flag) then
                call uEMEP_read_weekly_shipping_asi_data
            elseif (read_monthly_and_daily_shipping_data_flag) then
                call uEMEP_read_monthly_and_daily_shipping_asi_data
            else
                call uEMEP_read_shipping_asi_data
            endif
            !Read in the temperature fields from the alternative meteorology always, since EMEP data should not exist yet
            !use_alternative_meteorology_flag=.true.
            !call uEMEP_read_meteo_nc
            !Need to make a cross reference here or simply skip the two based on an if statement
            call uEMEP_read_time_profiles
            call uEMEP_set_emission_factors
            call uEMEP_convert_proxy_to_emissions
            
        endif
    
    endif
    enddo

    call uEMEP_save_emission_netcdf
    
	write(unit_logfile,'(A)') '================================================================'
    write(unit_logfile,'(a)') 'Finished saving emission data for EMEP.'
	write(unit_logfile,'(A)') '================================================================'
    
    stop
    
    end subroutine uEMEP_calculate_emissions_for_EMEP
    
    
    
    subroutine uEMEP_save_emission_netcdf
            
    use uEMEP_definitions
    
    implicit none

    character(256) temp_date_str,temp_compound_str,variable_type,unit_str,temp_name,var_name_temp,title_str
    logical create_file
    real scale_factor,valid_min
    integer i_file,i_pollutant,i_source
    real, allocatable :: temp_subgrid(:,:,:)
    integer exists
    integer temp_time_dim
    integer i,j
    
    temp_time_dim=save_emissions_end_index-save_emissions_start_index+1
    write(unit_logfile,'(A,3i6)') 'Time dimensions to be saved: ',save_emissions_start_index,save_emissions_end_index,temp_time_dim
    if (.not.allocated(temp_subgrid)) allocate(temp_subgrid(emission_subgrid_dim(x_dim_index,allsource_index),emission_subgrid_dim(y_dim_index,allsource_index),temp_time_dim))

    valid_min=0.   
    unit_str="ug/m3"
    !variable_type='byte'
    !variable_type='double'
    variable_type='float'
    scale_factor=1.
    
    !Save the data
    !i_file=subgrid_total_file_index(allsource_index)
    !temp_name=trim(pathname_grid(i_file))//trim(filename_grid(i_file))//'_'//trim(file_tag)//'.nc'
    if (len(emission_date_str).gt.0) then
        temp_date_str='_'//trim(emission_date_str)
    else
        temp_date_str=''
    endif
    
    !Do not use the actual emission start date for the file name but use the format specified by filename_date_output_grid
    if (len(filename_date_output_grid).gt.0) then
        temp_date_str='_'//trim(filename_date_output_grid)
    else
        temp_date_str=''
    endif
   
   
    !Do not write 'all' in file name if all compounds are selected
    if (pollutant_index.eq.all_nc_index) then
        temp_compound_str=''
    else
        temp_compound_str='_'//trim(var_name_nc(conc_nc_index,compound_index,allsource_index))
    endif
    
    title_str='uEMEP_emission_'//trim(file_tag)//temp_date_str
    i_file=subgrid_total_file_index(allsource_index)

    write(unit_logfile,'(A)') ''
    write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Saving emission netcdf data (uEMEP_save_emission_netcdf)'
	write(unit_logfile,'(A)') '================================================================'
    
    inquire(directory=trim(pathname_emissions_for_EMEP),exist=exists)
    if (.not.exists) then
        write(unit_logfile,'(A)')'ERROR: Path to EMEP emission output '//trim(pathname_emissions_for_EMEP)//' does not exist.'
        stop
    endif

   !Save the emissions interpolated to the target grid
    variable_type='float'
    unit_str="mg/m2"
    do i_source=1,n_source_index
    if (save_emissions_for_EMEP(i_source).and.i_source.ne.allsource_index) then
        !Create a new file for each source
        create_file=.true.
        !temp_name=trim(pathname_emissions_for_EMEP)//'uEMEP_emission_for_EMEP_'//trim(file_tag)//'_'//trim(source_file_str(i_source))//trim(temp_compound_str)//trim(temp_date_str)//'_'//trim(forecast_hour_str)//'.nc'
        temp_name=trim(pathname_emissions_for_EMEP)//'uEMEP_emission_for_EMEP_'//trim(file_tag)//'_'//trim(source_file_str(i_source))//trim(temp_compound_str)//trim(temp_date_str)//'.nc'
        
        write(unit_logfile,'(a,a)') 'Saving netcdf file: ',trim(temp_name)
        
        do i_pollutant=1,n_pollutant_loop
        !if (pollutant_loop_index(i_pollutant).ne.pmex_nc_index.and.pollutant_loop_index(i_pollutant).ne.pm10_sand_nc_index.and.pollutant_loop_index(i_pollutant).ne.pm10_salt_nc_index &
        !    .and.pollutant_loop_index(i_pollutant).ne.pm25_sand_nc_index.and.pollutant_loop_index(i_pollutant).ne.pm25_salt_nc_index) then
        if (pollutant_loop_index(i_pollutant).eq.pm10_nc_index.or.pollutant_loop_index(i_pollutant).eq.pm25_nc_index.or.pollutant_loop_index(i_pollutant).eq.nox_nc_index.or.pollutant_loop_index(i_pollutant).eq.nh3_nc_index) then
                i_file=emission_file_index(i_source)
                var_name_temp=trim(var_name_nc(emis_nc_index,pollutant_loop_index(i_pollutant),allsource_index)) !//'_'//trim(filename_grid(i_file))
                
                !Calculate the emissions in the target grid            
                temp_subgrid(:,:,:)=emission_subgrid(:,:,save_emissions_start_index:save_emissions_end_index,i_source,i_pollutant)

                !Convert the PM10 to PMco, special case
                if (pollutant_loop_index(i_pollutant).eq.pm10_nc_index) then
                    var_name_temp=trim(var_name_nc(emis_nc_index,pmco_nc_index,allsource_index)) !//'_'//trim(filename_grid(i_file))
                    temp_subgrid(:,:,:)=emission_subgrid(:,:,save_emissions_start_index:save_emissions_end_index,i_source,pollutant_loop_back_index(pm10_nc_index))-emission_subgrid(:,:,save_emissions_start_index:save_emissions_end_index,i_source,pollutant_loop_back_index(pm25_nc_index))
                endif

                !Subgrid emissions are in units ug/sec/subgrid. Convert to mg/m2/hour. Acount for the difference in subgrid sizes here             
                if (trim(save_emissions_for_EMEP_projection).eq.'lambert') then
                    temp_subgrid(:,:,:)=1.0e-3*3600.*temp_subgrid(:,:,:)/(emission_subgrid_delta(y_dim_index,i_source)*emission_subgrid_delta(x_dim_index,i_source))
                else
                    !Temporary estimate of area of lat lon. Needs to be fixed
                    do j=1,emission_subgrid_dim(y_dim_index,i_source)
                    do i=1,emission_subgrid_dim(x_dim_index,i_source)
                    
                        temp_subgrid(i,j,:)=1.0e-3*3600.*temp_subgrid(i,j,:)/(emission_subgrid_delta(y_dim_index,i_source)*emission_subgrid_delta(x_dim_index,i_source) &
                            *110570.*110570.*cos(3.14159/180.*y_emission_subgrid(i,j,i_source)))
                    enddo
                    enddo
                endif
                
                !write(*,'(4i,a)') i_pollutant,i_file,i_source,pollutant_loop_index(i_pollutant),trim(var_name_temp)
                
                if (save_netcdf_file_flag.or.save_netcdf_receptor_flag) then
                    write(unit_logfile,'(a,es12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid(:,:,:))/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_for_EMEP_netcdf_file(unit_logfile,temp_name,emission_subgrid_dim(x_dim_index,i_source),emission_subgrid_dim(y_dim_index,i_source),temp_time_dim &
                        ,temp_subgrid(:,:,:),x_emission_subgrid(:,:,i_source),y_emission_subgrid(:,:,i_source),lon_emission_subgrid(:,:,i_source),lat_emission_subgrid(:,:,i_source),var_name_temp &
                        ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                endif
                !Do not create file after first loop
                create_file=.false.
        endif
        enddo
    endif
    enddo
        
    end subroutine uEMEP_save_emission_netcdf
    
    subroutine uEMEP_save_for_EMEP_netcdf_file(unit_logfile_in,filename_netcdf,nx,ny,nt,val_array,x_array,y_array,lon_array,lat_array,name_array,unit_array,title_str,create_file,valid_min,variable_type,scale_factor)
    
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
    character(256) variable_type
    real scale_factor
    
    integer ncid
    integer y_dimid,x_dimid,lat_dimid,lon_dimid,val_dimid,time_dimid
    integer y_varid,x_varid,lat_varid,lon_varid,val_varid,time_varid,proj_varid
    integer dimids3(3),dimids2(2),chunks3(3)
    integer n_dims(3)
    integer status
    integer nf90_type
    

    if (trim(variable_type).eq.'byte') nf90_type=NF90_BYTE
    if (trim(variable_type).eq.'short') nf90_type=NF90_SHORT
    if (trim(variable_type).eq.'float') nf90_type=NF90_FLOAT
    if (trim(variable_type).eq.'double') nf90_type=NF90_DOUBLE
    
    !Assumes x and y are the dimensions   
    x_vector=x_array(:,1)
    y_vector=y_array(1,:)
    !write(*,*) x_vector
    !write(*,*) y_vector
    
    if (create_file) then
        !Create a netcdf file
        !call check(  nf90_create(filename_netcdf, nf90_clobber, ncid) )
        !call check(  nf90_create(filename_netcdf, NF90_HDF5, ncid) )
        call check(  nf90_create(filename_netcdf, IOR(NF90_HDF5, NF90_CLASSIC_MODEL), ncid) ) !New

        !Specify global attributes
        call check(  nf90_put_att(ncid, nf90_global, "Conventions", "CF-1.4" ) )
        call check(  nf90_put_att(ncid, nf90_global, "title", trim(title_str)) )
        call check(  nf90_put_att(ncid, nf90_global, "Model", "uEMEP emissions for EMEP" ) )        
    
        !Projection data
        if (trim(save_emissions_for_EMEP_projection).eq.'lambert') then
        call check(  nf90_def_var(ncid, "projection_lambert", NF90_int, proj_varid) )
        call check(  nf90_put_att(ncid, proj_varid, "standard_parallel", EMEP_projection_attributes(1:2) ) )
        call check(  nf90_put_att(ncid, proj_varid, "longitude_of_central_meridian", EMEP_projection_attributes(3) ) )

        call check(  nf90_put_att(ncid, proj_varid, "grid_mapping_name", "lambert_conformal_conic" ) )
        call check(  nf90_put_att(ncid, proj_varid, "latitude_of_projection_origin", EMEP_projection_attributes(4) ) )
        call check(  nf90_put_att(ncid, proj_varid, "earth_radius", EMEP_projection_attributes(5) ) )
        endif
        
        !Define the dimensions
        !write(*,*) 'here1',nt
        call check(  nf90_def_dim(ncid,"time",NF90_UNLIMITED, time_dimid) )
        !write(*,*) 'herea'
        !if (trim(save_emissions_for_EMEP_projection).eq.'lambert') then
        call check(  nf90_def_dim(ncid, "y", ny, y_dimid) )
        call check(  nf90_def_dim(ncid, "x", nx, x_dimid) )
        !else
        !call check(  nf90_def_dim(ncid, "lat", ny, y_dimid) )
        !call check(  nf90_def_dim(ncid, "lon", nx, x_dimid) )
        !endif
        
        !write(*,*) 'here2'
        !Define the dimension variables
        call check(  nf90_def_var(ncid, "time", NF90_DOUBLE, time_dimid, time_varid) )
        !call check(  nf90_def_var(ncid, "time", NF90_INT, time_dimid, time_varid) )
        !if (trim(save_emissions_for_EMEP_projection).eq.'lambert') then
        call check(  nf90_def_var(ncid, "y", NF90_REAL, y_dimid, y_varid) )
        call check(  nf90_def_var(ncid, "x", NF90_REAL, x_dimid, x_varid) )
        !else
        !call check(  nf90_def_var(ncid, "lat", NF90_REAL, y_dimid, y_varid) )
        !call check(  nf90_def_var(ncid, "lon", NF90_REAL, x_dimid, x_varid) )
        !endif
        
        !write(*,*) 'here3'
        !Define the values
        dimids3 = (/ x_dimid, y_dimid, time_dimid /)
        dimids2 = (/ x_dimid, y_dimid /)
        !if (trim(save_emissions_for_EMEP_projection).eq.'lambert') then
        call check(  nf90_def_var(ncid, "lat", NF90_REAL, dimids2, lat_varid) )
        call check(  nf90_def_var(ncid, "lon", NF90_REAL, dimids2, lon_varid) )
        !endif
    

        !write(*,*) 'here4'
        !Specify the units
        call check(  nf90_put_att(ncid, lat_varid, "units", "degrees_north") )
        call check(  nf90_put_att(ncid, lon_varid, "units", "degrees_east") )
        if (trim(save_emissions_for_EMEP_projection).eq.'lambert') then
        call check(  nf90_put_att(ncid, y_varid, "units", "m") )
        call check(  nf90_put_att(ncid, x_varid, "units", "m") )
        else
        call check(  nf90_put_att(ncid, y_varid, "units", "degrees_north") )
        call check(  nf90_put_att(ncid, x_varid, "units", "degrees_east") )
        endif
        call check(  nf90_put_att(ncid, time_varid, "units", trim(unit_dim_nc(time_dim_nc_index))) )

        !write(*,*) 'here5'
        !Specify other dimension attributes
        !if (trim(save_emissions_for_EMEP_projection).eq.'lambert') then
        call check(  nf90_put_att(ncid, y_varid, "standard_name", "projection_y_coordinate") )
        call check(  nf90_put_att(ncid, x_varid, "standard_name", "projection_x_coordinate") )
        call check(  nf90_put_att(ncid, y_varid, "axis", "Y") )
        call check(  nf90_put_att(ncid, x_varid, "axis", "X") )
        !else
        !call check(  nf90_put_att(ncid, y_varid, "standard_name", "latitude") )
        !call check(  nf90_put_att(ncid, x_varid, "standard_name", "longitude") )
        !endif
        
        !Close the definitions
        call check( nf90_enddef(ncid) )
        
        !write(*,*) 'here6',shape(val_dim_nc(1:dim_length_nc(time_dim_nc_index),time_dim_nc_index))
        !write(*,*) val_dim_nc(1:dim_length_nc(time_dim_nc_index),time_dim_nc_index)

        call check( nf90_put_var(ncid, time_varid, val_dim_nc(save_emissions_start_index:save_emissions_end_index,time_dim_nc_index)) )
        !call check( nf90_put_var(ncid, time_varid, time_seconds_output(1:dim_length_nc(time_dim_nc_index))) )
        !if (trim(save_emissions_for_EMEP_projection).eq.'lambert') then
        call check( nf90_put_var(ncid, y_varid, y_vector) )
        call check( nf90_put_var(ncid, x_varid, x_vector) )
        call check( nf90_put_var(ncid, lat_varid, lat_array) )
        call check( nf90_put_var(ncid, lon_varid, lon_array) )
        !else
        !call check( nf90_put_var(ncid, y_varid, y_vector) )
        !call check( nf90_put_var(ncid, x_varid, x_vector) )
        !endif

        call check( nf90_close(ncid) )
    
    endif
        
    !Add to the existing file
    call check( nf90_open(filename_netcdf, NF90_WRITE, ncid) )
    
    !Get the dimensions id from the existing file
    call check( nf90_inq_dimid(ncid,"time",time_dimid) )
    !if (trim(save_emissions_for_EMEP_projection).eq.'lambert') then
    call check( nf90_inq_dimid(ncid, "y", y_dimid) )
    call check( nf90_inq_dimid(ncid, "x", x_dimid) )
    !else
    !call check( nf90_inq_dimid(ncid, "lat", y_dimid) )
    !call check( nf90_inq_dimid(ncid, "lon", x_dimid) )
    !endif
    dimids3 = (/ x_dimid, y_dimid, time_dimid /)
    chunks3 = (/ nx, ny, 1 /) !New
    call check( nf90_inquire_dimension(ncid, dimids3(1), temp_name, n_dims(1)) )
    call check( nf90_inquire_dimension(ncid, dimids3(2), temp_name, n_dims(2)) )
    call check( nf90_inquire_dimension(ncid, dimids3(3), temp_name, n_dims(3)) )

    !write(*,*) 'here7'
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
    !write(*,*) 'here8'
        if (trim(save_emissions_for_EMEP_projection).eq.'lambert') then
        call check(  nf90_put_att(ncid, val_varid, "grid_mapping", "projection_lambert") )
        else
        call check(  nf90_put_att(ncid, val_varid, "grid_mapping", "latitude_longitude") )
        endif
    !write(*,*) 'here9'
        call check(  nf90_put_att(ncid, val_varid, "coordinates", "lon lat") )
        if (scale_factor.ne.1.) call check(  nf90_put_att(ncid, val_varid, "scale_factor", scale_factor) )
        
        !Close the definitions
        call check( nf90_enddef(ncid) )
    !write(*,*) 'here10'

    endif

    
    if (use_single_time_loop_flag) then
        write(unit_logfile,'(A)') 'ERROR: Saving emissions for EMEP will not work when use_single_time_loop_flag=.true. Set to false'
        stop      
        !Add time to the time dimension       
        call check( nf90_inq_varid(ncid, "time", time_varid) )
            !write(*,*) 'here11a'
        !call check( nf90_inquire_dimension(ncid, time_dimid, temp_name, n_dims(3)) )
        !n_dims(3)=n_dims(3)+1
        n_dims(3)=t_loop
        !write(*,*) n_dims(3),val_dim_nc(1,time_dim_nc_index)
        call check( nf90_put_var(ncid, time_varid, val_dim_nc(1,time_dim_nc_index), start = (/n_dims(3)/) ) )
        !write(*,*) n_dims(3),val_dim_nc(1,time_dim_nc_index)
        !write(*,*) n_dims
            !write(*,*) 'here11b'
        
        !Add dimension and array to existing
        call check( nf90_inq_varid(ncid, trim(name_array), val_varid) )
        if (nf90_type.eq.NF90_byte) then
            call check( nf90_put_var(ncid, val_varid, int1(val_array), start=(/1,1,n_dims(3)/), count=(/n_dims(1),n_dims(2),1/)) )
        elseif (nf90_type.eq.NF90_short) then
            call check( nf90_put_var(ncid, val_varid, int2(val_array), start=(/1,1,n_dims(3)/), count=(/n_dims(1),n_dims(2),1/)) )
        else
            call check( nf90_put_var(ncid, val_varid, val_array, start=(/1,1,n_dims(3)/), count=(/n_dims(1),n_dims(2),1/)) )
        endif
            !write(*,*) 'here11'

    else
        
        !Write the variable to file. Default is float
        if (nf90_type.eq.NF90_byte) then
            call check( nf90_put_var(ncid, val_varid, int1(val_array)) )
        elseif (nf90_type.eq.NF90_short) then
            call check( nf90_put_var(ncid, val_varid, int2(val_array)) )
        else
            !write(*,*) ncid, val_varid, shape(val_array)
            call check( nf90_put_var(ncid, val_varid, val_array) )
        endif

    endif
    
    call check( nf90_close(ncid) )
    

    
    end subroutine uEMEP_save_for_EMEP_netcdf_file
    
