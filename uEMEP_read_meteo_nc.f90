!uEMEP_read_meteo_nc
    !Reads in AROME data in 2 files, 3 for meteo and 4 for z0
    !These two files must have the same x,y dimensions as they are placed in the same grid
    
    subroutine uEMEP_read_meteo_nc
    
    use uEMEP_definitions
    use netcdf
    
    implicit none
    
    integer i,j,k,t
    integer ii,jj
    logical exists
    character(256) pathfilename_nc
    integer status_nc     !Error message
    integer id_nc
    integer dim_id_nc(num_dims_meteo_nc)
    character(256) dimname_temp,var_name_nc_temp,unit_name_nc_temp
    integer var_id_nc
    real :: local_fraction_scaling=1.0
    integer i_file,i_source,i_conc,i_dim
    integer temp_frac_index,temp_file_index,temp_compound_index,temp_source_index
    integer temp_num_dims
    integer temp_start_time_nc_index,temp_end_time_nc_index
    integer temp_start_time_meteo_nc_index,temp_end_time_meteo_nc_index
    integer i_loop
    integer valid_dim_length_meteo_nc(num_dims_meteo_nc) !dimensions of file 3
    integer numAtts_projection
    logical :: invert_levels_flag=.false.
    integer surface_level_nc_2
    
    real temp_lat(4),temp_lon(4)
    real temp_y(4),temp_x(4)
    real temp_x_min,temp_x_max,temp_y_min,temp_y_max
    integer i_temp_min,i_temp_max,j_temp_min,j_temp_max
    double precision temp_var1d_nc_dp(2,2)
    real temp_delta(2)
    real H_emep_temp
    real scale_grid_interpolation_size(2)
    real EMEP_temp_delta(2)
    real temp_lat_mean
    
    integer n_file,n_file_start
    double precision date_num_temp
    integer date_array(6)
    double precision scale_factor_nc
    integer temp_nc_projection_type
    double precision :: temp_nc_projection_attributes(10)

    logical found_file
    integer :: search_hour_step=6
    integer new_start_date_input(6)
    character(256) format_temp
    double precision date_to_number
    character(256) replace_string_char
    
    real EMEP_grid_interpolation_size_temp
    
    !Temporary reading rvariables
    double precision, allocatable :: var1d_nc_dp(:)
    double precision, allocatable :: var2d_nc_dp(:,:)
    
    !Temporary files for roatating wind field
    real, allocatable :: temp_meteo_var3d_nc(:,:,:,:)
    
    !Daily mean temperature variables
    integer DMT_start_time_nc_index,DMT_end_time_nc_index,DMT_dim_length_nc

    
    write(unit_logfile,'(A)') ''
    write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Reading additional meteo data (uEMEP_read_meteo_nc)'
	write(unit_logfile,'(A)') '================================================================'

    
    !This if statement is already specified in uEMEP_define_subgrid and is not necessary here
    if (hourly_calculations) then
        temp_start_time_meteo_nc_index=start_time_meteo_nc_index
        temp_end_time_meteo_nc_index=end_time_meteo_nc_index
    else
        temp_start_time_meteo_nc_index=1
        temp_end_time_meteo_nc_index=1
    endif
     
    if (use_single_time_loop_flag) then
        temp_start_time_meteo_nc_index=start_time_meteo_nc_index+t_loop-1
        temp_end_time_meteo_nc_index=temp_start_time_meteo_nc_index       
    endif

    !Presettng the surface level to 1. Valid when there is no inverting of layers
    surface_level_nc=1
    write(unit_logfile,'(A,I)') ' Surface level base set to: ',surface_level_nc

        if (allocated(val_dim_meteo_nc)) deallocate (val_dim_meteo_nc)
        if (allocated(var1d_nc_dp)) deallocate (var1d_nc_dp) 
        if (allocated(var2d_nc_dp)) deallocate (var2d_nc_dp)

        if (allocated(meteo_var1d_nc)) deallocate (meteo_var1d_nc)
        if (allocated(meteo_var2d_nc)) deallocate (meteo_var2d_nc)
        if (allocated(meteo_var3d_nc)) deallocate (meteo_var3d_nc)
        if (allocated(meteo_var4d_nc)) deallocate (meteo_var4d_nc)

    !Loop through the meteorological files containing the data
    if (use_alternative_meteorology_flag) then
        n_file_start=3
        n_file=3
    endif
    if (use_alternative_z0_flag) then
        n_file_start=4
        n_file=4
    endif
    if (use_alternative_meteorology_flag.and.use_alternative_z0_flag) then
        n_file_start=3
        n_file=4
    endif
    
    do i_file=n_file_start,n_file
        
        !Set the filename
        pathfilename_EMEP(i_file)=trim(pathname_EMEP(i_file))//trim(filename_EMEP(i_file))
     
        !Test existence of the filename 4. If does not exist then stop
        if (i_file.eq.4) then
        inquire(file=trim(pathfilename_EMEP(i_file)),exist=exists)
        if (.not.exists) then
            write(unit_logfile,'(A,A)') ' ERROR: Netcdf file does not exist: ', trim(pathfilename_EMEP(i_file))
            write(unit_logfile,'(A)') '  STOPPING'
            stop
        endif
        endif

        !Test existence of the filename. If does not exist then try 6 hours before
        if (i_file.eq.3) then
        
        inquire(file=trim(pathfilename_EMEP(i_file)),exist=exists)
        if (.not.exists) then
            write(unit_logfile,'(A,A)') ' WARNING: Meteo netcdf file does not exist: ',  trim(pathfilename_EMEP(i_file))
            write(unit_logfile,'(A)') ' Will try 6 hours before 4 times'
        
            !Start search back 6 hours
            found_file=.false.
            do i=1,4
                if (hourly_calculations) then
                    temp_start_time_meteo_nc_index=start_time_meteo_nc_index+search_hour_step*(i)
                    temp_end_time_meteo_nc_index=end_time_meteo_nc_index+search_hour_step*(i)
                endif     
                if (use_single_time_loop_flag) then
                    temp_start_time_meteo_nc_index=start_time_meteo_nc_index+t_loop-1+search_hour_step*(i)
                    temp_end_time_meteo_nc_index=temp_start_time_meteo_nc_index+search_hour_step*(i)    
                endif
                !Create new date_str
                format_temp='yyyymmddHH'
                call datestr_to_date(config_date_str,format_temp,new_start_date_input)
                date_num_temp=date_to_number(new_start_date_input,ref_year_meteo)
                call number_to_date(date_num_temp-dble(search_hour_step*i+0.5)/dble(24.),new_start_date_input,ref_year_meteo)
                !Replace replacement_date_str with <yyyyhhmm> so the new_start_date_input can be inserted
                format_temp='<yyyymmdd>'
                filename_EMEP(i_file)=replace_string_char(format_temp,replacement_date_str,original_filename_EMEP(i_file))
                pathname_EMEP(i_file)=replace_string_char(format_temp,replacement_date_str,original_pathname_EMEP(i_file))
                !write(*,*) trim(filename_EMEP(i_file)),'  ',trim(replacement_date_str)
                !Replace replacement_hour_str with <HH> so the forecast hour can be inserted
                format_temp='<HH>'
                filename_EMEP(i_file)=replace_string_char(format_temp,replacement_hour_str,filename_EMEP(i_file))
                pathname_EMEP(i_file)=replace_string_char(format_temp,replacement_hour_str,pathname_EMEP(i_file))
                !write(*,*) trim(filename_EMEP(i_file)),'  ',trim(forecast_hour_str)
                !Replace datestr twice for both forecast_hour and config_date
                call date_to_datestr_bracket(new_start_date_input,filename_EMEP(i_file),filename_EMEP(i_file))
                call date_to_datestr_bracket(new_start_date_input,pathname_EMEP(i_file),pathname_EMEP(i_file))
                call date_to_datestr_bracket(new_start_date_input,filename_EMEP(i_file),filename_EMEP(i_file))
                call date_to_datestr_bracket(new_start_date_input,pathname_EMEP(i_file),pathname_EMEP(i_file))
                pathfilename_EMEP(i_file)=trim(pathname_EMEP(i_file))//trim(filename_EMEP(i_file))
                write(unit_logfile,'(A,A)') ' Trying: ', trim(pathfilename_EMEP(i_file))
                inquire(file=trim(pathfilename_EMEP(i_file)),exist=exists)
                if (exists) then
                    found_file=.true.
                    exit
                else 
                    found_file=.false.
                endif
            enddo
        
            if (.not.found_file) then
                write(unit_logfile,'(A,A)') ' ERROR: Meteo netcdf file still does not exist: ', trim(pathfilename_EMEP(i_file))
                write(unit_logfile,'(A)') ' STOPPING'
                stop
            else
                write(unit_logfile,'(A,A)') ' Found earlier meteo netcdf file: ', trim(pathfilename_EMEP(i_file))
                write(unit_logfile,'(A,2i6)') ' New start and end index: ', temp_start_time_meteo_nc_index,temp_end_time_meteo_nc_index
            endif
        
        endif
    
        endif
    
        !Open the netcdf file for reading
        write(unit_logfile,'(2A)') ' Opening netcdf file: ',trim(pathfilename_EMEP(i_file))
        status_nc = NF90_OPEN (pathfilename_EMEP(i_file), nf90_nowrite, id_nc)
        if (status_nc .NE. NF90_NOERR) then
            write(unit_logfile,'(A,I)') 'ERROR opening netcdf file. Stopping: ',status_nc
            stop
        endif
        
        meteo_nc_projection_type=LL_projection_index
        
        !Find the projection. If no projection then in lat lon coordinates
        status_nc = NF90_INQ_VARID (id_nc,'projection_lambert',var_id_nc)
        
        if (status_nc.eq.NF90_NOERR) then
            !If there is a projection then read in the attributes. All these are doubles
            !status_nc = nf90_inquire_variable(id_nc, var_id_nc, natts = numAtts_projection)
                status_nc = nf90_get_att(id_nc, var_id_nc, 'standard_parallel', meteo_nc_projection_attributes(1:2))
                status_nc = nf90_get_att(id_nc, var_id_nc, 'longitude_of_central_meridian', meteo_nc_projection_attributes(3))
                status_nc = nf90_get_att(id_nc, var_id_nc, 'latitude_of_projection_origin', meteo_nc_projection_attributes(4))
                status_nc = nf90_get_att(id_nc, var_id_nc, 'earth_radius', meteo_nc_projection_attributes(5))
                meteo_nc_projection_type=LCC_projection_index
                        
            write(unit_logfile,'(A,5f12.2)') 'Reading lambert_conformal_conic projection. ',meteo_nc_projection_attributes(1:5)
            if (meteo_nc_projection_attributes(1).ne.meteo_nc_projection_attributes(4).or.meteo_nc_projection_attributes(2).ne.meteo_nc_projection_attributes(4)) then
                use_alternative_LCC_projection_flag=.true.
                write(unit_logfile,'(A,l)') 'Using alternative lambert_conformal_conic projection: ',use_alternative_LCC_projection_flag
            else
                use_alternative_LCC_projection_flag=.false.                
            endif
            !Always set to true. i.e. not use anymore
            use_alternative_LCC_projection_flag=.true.
        endif
        
        !Find the projection. If no projection then in lat lon coordinates
        status_nc = NF90_INQ_VARID (id_nc,'Polar_Stereographic',var_id_nc)

        if (status_nc.eq.NF90_NOERR) then

            EMEP_projection_attributes=0.
            EMEP_projection_attributes(5)=6.370e6
            status_nc = nf90_get_att(id_nc, var_id_nc, 'straight_vertical_longitude_from_pole', meteo_nc_projection_attributes(1))
            status_nc = nf90_get_att(id_nc, var_id_nc, 'latitude_of_projection_origin', meteo_nc_projection_attributes(2))
            status_nc = nf90_get_att(id_nc, var_id_nc, 'false_easting', meteo_nc_projection_attributes(3))
            status_nc = nf90_get_att(id_nc, var_id_nc, 'false_northing', meteo_nc_projection_attributes(4))
            status_nc = nf90_get_att(id_nc, var_id_nc, 'earth_radius', meteo_nc_projection_attributes(5))
            status_nc = nf90_get_att(id_nc, var_id_nc, 'scale_factor_at_projection_origin', meteo_nc_projection_attributes(6))
                
            meteo_nc_projection_type=PS_projection_index
                        
            write(unit_logfile,'(A,5f12.2)') 'Reading Polar_Stereographic projection. ',meteo_nc_projection_attributes(1:5)
        
        endif

        !Find the (x,y,z,time) dimensions of the file
        do i_dim=1,num_dims_meteo_nc
            status_nc = NF90_INQ_DIMID (id_nc,dim_name_meteo_nc(i_dim),dim_id_nc(i_dim))
            status_nc = NF90_INQUIRE_DIMENSION (id_nc,dim_id_nc(i_dim),dimname_temp,dim_length_meteo_nc(i_dim))
            if (status_nc .NE. NF90_NOERR) then
                write(unit_logfile,'(A,A,A,I)') 'No dimension information available for ',trim(dim_name_meteo_nc(i_dim)),' Setting to 1 with status: ',status_nc
                dim_length_meteo_nc(i_dim)=1
            endif
        enddo
        
        if (i_file.eq.3) then
            valid_dim_length_meteo_nc=dim_length_meteo_nc
        endif
        
        if (i_file.eq.3) then
            if (subgrid_dim(t_dim_index).gt.dim_length_meteo_nc(time_dim_nc_index)) then
                write(unit_logfile,'(A,2I)') 'ERROR: Specified time dimensions are greater than meteo netcdf dimensions. Stopping ',subgrid_dim(t_dim_index),dim_length_meteo_nc(time_dim_nc_index)
                stop
            endif
            if (temp_end_time_meteo_nc_index.gt.dim_length_meteo_nc(time_dim_nc_index)) then
                write(unit_logfile,'(A,2I)') 'ERROR: Required meteo time dimension larger than available meteo time dimension. Stopping ',temp_end_time_meteo_nc_index,dim_length_meteo_nc(time_dim_nc_index)
                stop
            endif
        endif

               
        write(unit_logfile,'(A,6I)') ' Size of meteo dimensions (x,y,z,t): ',dim_length_meteo_nc
        
        if (i_file.eq.3) then
            dim_start_meteo_nc(time_dim_nc_index)=temp_start_time_meteo_nc_index
            dim_length_meteo_nc(time_dim_nc_index)=min(dim_length_meteo_nc(time_dim_nc_index),subgrid_dim(t_dim_index))
        elseif (i_file.eq.4) then
            dim_start_meteo_nc(time_dim_nc_index)=1
            dim_length_meteo_nc(time_dim_nc_index)=1
        endif

        write(unit_logfile,'(A,6I)') ' New size of meteo dimensions (x,y,z,t): ',dim_length_meteo_nc
        
        
        !Calculate the necessary extent of the meteo_nc grid region and only read these grids
        if (reduce_EMEP_region_flag) then
            !Determine the LL cordinates of the target grid
            !EMEP_grid_interpolation_size_temp=max(EMEP_grid_interpolation_size*local_fraction_grid_size_scaling,EMEP_additional_grid_interpolation_size_original*local_fraction_grid_size_scaling)
            EMEP_grid_interpolation_size_temp=EMEP_grid_interpolation_size*local_fraction_grid_size_scaling

            !Retrieve the four corners of the target grid in lat and lon
            call PROJ2LL(init_subgrid_min(x_dim_index),init_subgrid_min(y_dim_index),temp_lon(1),temp_lat(1),projection_attributes,projection_type)
            call PROJ2LL(init_subgrid_max(x_dim_index),init_subgrid_max(y_dim_index),temp_lon(2),temp_lat(2),projection_attributes,projection_type)
            call PROJ2LL(init_subgrid_min(x_dim_index),init_subgrid_max(y_dim_index),temp_lon(3),temp_lat(3),projection_attributes,projection_type)
            call PROJ2LL(init_subgrid_max(x_dim_index),init_subgrid_min(y_dim_index),temp_lon(4),temp_lat(4),projection_attributes,projection_type)
            !call UTM2LL(utm_zone,init_subgrid_min(y_dim_index),init_subgrid_min(x_dim_index),temp_lat(1),temp_lon(1))
            !call UTM2LL(utm_zone,init_subgrid_max(y_dim_index),init_subgrid_max(x_dim_index),temp_lat(2),temp_lon(2))
            !call UTM2LL(utm_zone,init_subgrid_max(y_dim_index),init_subgrid_min(x_dim_index),temp_lat(3),temp_lon(3))
            !call UTM2LL(utm_zone,init_subgrid_min(y_dim_index),init_subgrid_max(x_dim_index),temp_lat(4),temp_lon(4))
        
            !Find the average for use later
            temp_lat_mean=sum(temp_lat)/4.
            
            temp_x_min=1.e32;temp_y_min=1.e32
            temp_x_max=-1.e32;temp_y_max=-1.e32

                    if (meteo_nc_projection_type.eq.LCC_projection_index) then
                        !Convert lat lon corners to lambert
                        do i=1,4
                            call lb2lambert2_uEMEP(temp_x(i),temp_y(i),temp_lon(i),temp_lat(i),meteo_nc_projection_attributes)
                        enddo            
                    elseif (meteo_nc_projection_type.eq.PS_projection_index) then
                        !Convert lat lon corners to lambert
                        do i=1,4
                            call LL2PS_spherical(temp_x(i),temp_y(i),temp_lon(i),temp_lat(i),meteo_nc_projection_attributes)
                        enddo            
                    elseif (meteo_nc_projection_type.eq.LL_projection_index) then
                        !Set lat lon corners if EMEP is in lat lon
                        temp_x=temp_lon;temp_y=temp_lat
                    else
                        !Otherwise assume the same coordinate system
                        temp_x(1)=init_subgrid_min(x_dim_index);temp_y(1)=init_subgrid_min(y_dim_index)
                        temp_x(2)=init_subgrid_max(x_dim_index);temp_y(2)=init_subgrid_min(y_dim_index)
                        temp_x(3)=init_subgrid_min(x_dim_index);temp_y(3)=init_subgrid_max(y_dim_index)
                        temp_x(4)=init_subgrid_max(x_dim_index);temp_y(4)=init_subgrid_max(y_dim_index)
                    endif
                    
                do i=1,4
                    if (temp_x(i).lt.temp_x_min) temp_x_min=temp_x(i)
                    if (temp_y(i).lt.temp_y_min) temp_y_min=temp_y(i)
                    if (temp_x(i).gt.temp_x_max) temp_x_max=temp_x(i)
                    if (temp_y(i).gt.temp_y_max) temp_y_max=temp_y(i)
                enddo
            
                !Read in the first 2 x and y position values from the nc file to get min values and delta values
                !write(*,*) temp_x_min,temp_x_max,temp_y_min,temp_y_max
                
                status_nc = NF90_INQ_VARID (id_nc, trim(dim_name_meteo_nc(x_dim_nc_index)), var_id_nc)
                status_nc = NF90_GET_VAR (id_nc, var_id_nc,temp_var1d_nc_dp(1,1:2),start=(/1/),count=(/2/))
                status_nc = NF90_INQ_VARID (id_nc, trim(dim_name_meteo_nc(y_dim_nc_index)), var_id_nc)
                status_nc = NF90_GET_VAR (id_nc, var_id_nc,temp_var1d_nc_dp(2,1:2),start=(/1/),count=(/2/))
                status_nc = nf90_get_att(id_nc, var_id_nc, "units", unit_name_nc_temp)
                if (trim(unit_name_nc_temp).eq.'km') then
                    write(unit_logfile,'(A)') 'Units of x y data are in kilometres. Converting to metres'
                    temp_var1d_nc_dp=temp_var1d_nc_dp*1000.
                endif
                
                !HERE FIX. The EMEP_grid_interpolation_size is too small when using EMEP is a different grid to the meteo grid. Need to rescale this somehow
                !By using meteo_dgrid_nc(lon_nc_index), not defined yet, and dgrid_nc(lon_nc_index)
                !For example  dx_temp=111000.*dgrid_nc(lon_nc_index)*cos(lat_temp*pi/180.) and dy_temp=111000.*dgrid_nc(lat_nc_index)
  
                !write(*,*) temp_var1d_nc_dp
                temp_delta(1)=temp_var1d_nc_dp(1,2)-temp_var1d_nc_dp(1,1)
                temp_delta(2)=temp_var1d_nc_dp(2,2)-temp_var1d_nc_dp(2,1)
                !write(*,*) temp_delta

                if ((meteo_nc_projection_type.eq.LCC_projection_index.and.EMEP_projection_type.eq.LCC_projection_index) &
                .or.(meteo_nc_projection_type.eq.LL_projection_index.and.EMEP_projection_type.eq.LL_projection_index) &
                    .or.(meteo_nc_projection_type.eq.PS_projection_index.and.EMEP_projection_type.eq.PS_projection_index)) then
                    !If both EMEP and meteo are in the same coordinates then set the EMEP size to be the same as the meteo siz
                    EMEP_temp_delta(1)=dgrid_nc(lon_nc_index)
                    EMEP_temp_delta(2)=dgrid_nc(lat_nc_index)
                    elseif ((meteo_nc_projection_type.eq.LCC_projection_index.and.EMEP_projection_type.eq.LL_projection_index) &
                    .or.(meteo_nc_projection_type.eq.PS_projection_index.and.EMEP_projection_type.eq.LL_projection_index)) then
                     !EMEP is in latlon, convert to local coordinates
                    EMEP_temp_delta(1)=dgrid_nc(lon_nc_index)*111000.*cos(temp_lat_mean*3.14159/180.) 
                    EMEP_temp_delta(2)=dgrid_nc(lat_nc_index)*111000.
                    elseif ((meteo_nc_projection_type.eq.LL_projection_index.and.EMEP_projection_type.eq.LCC_projection_index)&
                    .or.(meteo_nc_projection_type.eq.LL_projection_index.and.EMEP_projection_type.eq.PS_projection_index)) then
                    !This conversion not available
                    write(unit_logfile,'(A,3I)') 'Use of lat lon projection in meteo data, together with Lambert or PS in EMEP not available. Stopping'
                    stop
                else
                    write(unit_logfile,'(A,3I)') 'Use of current projections in meteo and EMEP data not available. Stopping'
                    stop
                endif
                
                scale_grid_interpolation_size=EMEP_temp_delta/temp_delta
                !write(*,*) dgrid_nc(lon_nc_index),dgrid_nc(lat_nc_index)
                !write(*,*) EMEP_temp_delta
                !write(*,*) scale_grid_interpolation_size
                !stop
                
                !Find grid position of the max and min coordinates and add2 grids*EMEP_grid_interpolation_size
                i_temp_min=1+floor((temp_x_min-temp_var1d_nc_dp(1,1))/temp_delta(1)+0.5)
                i_temp_max=1+floor((temp_x_max-temp_var1d_nc_dp(1,1))/temp_delta(1)+0.5)
                j_temp_min=1+floor((temp_y_min-temp_var1d_nc_dp(2,1))/temp_delta(2)+0.5)
                j_temp_max=1+floor((temp_y_max-temp_var1d_nc_dp(2,1))/temp_delta(2)+0.5)
                !write(unit_logfile,'(A,2I)') ' Reading EMEP i grids: ',i_temp_min,i_temp_max
                !write(unit_logfile,'(A,2I)') ' Reading EMEP j grids: ',j_temp_min,j_temp_max
                i_temp_min=max(1,i_temp_min-1-ceiling(scale_grid_interpolation_size(1)*EMEP_grid_interpolation_size_temp))
                i_temp_max=min(dim_length_meteo_nc(x_dim_nc_index),i_temp_max+1+ceiling(scale_grid_interpolation_size(1)*EMEP_grid_interpolation_size_temp))
                j_temp_min=max(1,j_temp_min-1-ceiling(scale_grid_interpolation_size(2)*EMEP_grid_interpolation_size_temp))
                j_temp_max=min(dim_length_meteo_nc(y_dim_nc_index),j_temp_max+1+ceiling(scale_grid_interpolation_size(2)*EMEP_grid_interpolation_size_temp))
                dim_length_meteo_nc(x_dim_nc_index)=i_temp_max-i_temp_min+1
                dim_length_meteo_nc(y_dim_nc_index)=j_temp_max-j_temp_min+1
                dim_start_meteo_nc(x_dim_nc_index)=i_temp_min
                dim_start_meteo_nc(y_dim_nc_index)=j_temp_min
                write(unit_logfile,'(A,3I)') ' Reading meteo x grids: ',i_temp_min,i_temp_max,dim_length_meteo_nc(x_dim_nc_index)
                write(unit_logfile,'(A,3I)') ' Reading meteo y grids: ',j_temp_min,j_temp_max,dim_length_meteo_nc(y_dim_nc_index)

                !Set the new valid meteo dimensions
                if (i_file.ge.3) then
                    valid_dim_length_meteo_nc(x_dim_nc_index)=dim_length_meteo_nc(x_dim_nc_index)
                    valid_dim_length_meteo_nc(y_dim_nc_index)=dim_length_meteo_nc(y_dim_nc_index)
                endif

        endif
        
        !Allocate the nc arrays for reading
        if (.not.allocated(val_dim_meteo_nc)) allocate (val_dim_meteo_nc(maxval(dim_length_meteo_nc),num_dims_meteo_nc)) !x, y, z and time dimension values
        if (.not.allocated(unit_dim_meteo_nc)) allocate (unit_dim_meteo_nc(num_dims_meteo_nc)) !x, y, z and time dimension values
        if (.not.allocated(var1d_nc_dp)) allocate (var1d_nc_dp(maxval(dim_length_meteo_nc))) 
        if (.not.allocated(var2d_nc_dp)) allocate (var2d_nc_dp(dim_length_meteo_nc(x_dim_nc_index),dim_length_meteo_nc(y_dim_nc_index))) !Lat and lon

        !Allocate array for the alternative meteo files
        if (i_file.ge.3) then
            if (.not.allocated(meteo_var1d_nc)) allocate (meteo_var1d_nc(maxval(dim_length_meteo_nc),num_dims_meteo_nc)) !x, y, z and time maximum dimensions
            if (.not.allocated(meteo_var2d_nc)) allocate (meteo_var2d_nc(dim_length_meteo_nc(x_dim_nc_index),dim_length_meteo_nc(y_dim_nc_index),2)) !Lat and lon
            if (.not.allocated(meteo_var3d_nc)) allocate (meteo_var3d_nc(dim_length_meteo_nc(x_dim_nc_index),dim_length_meteo_nc(y_dim_nc_index),0:dim_length_meteo_nc(time_dim_nc_index),num_var_meteo_nc))
            if (.not.allocated(meteo_var4d_nc)) allocate (meteo_var4d_nc(dim_length_meteo_nc(x_dim_nc_index),dim_length_meteo_nc(y_dim_nc_index),dim_length_meteo_nc(z_dim_nc_index),0:dim_length_meteo_nc(time_dim_nc_index),num_var_meteo_nc))
            
        endif

        !Read in the dimensions and check values of the dimensions.
        do i=1,num_dims_meteo_nc
            status_nc = NF90_INQ_VARID (id_nc, trim(dim_name_meteo_nc(i)), var_id_nc)
            !write(*,*) id_nc, trim(dim_name_nc(i)), var_id_nc(i),dim_length_nc(i)
            var1d_nc_dp=0.
            !write(*,*) 'HERE',i,dim_start_nc(i),dim_length_nc(i)
            unit_dim_meteo_nc(i)=''
            if (status_nc .EQ. NF90_NOERR) then
                status_nc = nf90_get_att(id_nc, var_id_nc, "units", unit_dim_meteo_nc(i))
                status_nc = NF90_GET_VAR (id_nc, var_id_nc,var1d_nc_dp(1:dim_length_meteo_nc(i)),start=(/dim_start_meteo_nc(i)/),count=(/dim_length_meteo_nc(i)/));meteo_var1d_nc(1:dim_length_meteo_nc(i),i)=real(var1d_nc_dp(1:dim_length_meteo_nc(i)))  
                !Use the first file to give valid time stamps
                if (i_file.eq.3.and.i.eq.time_dim_nc_index) then
                    val_dim_meteo_nc(1:dim_length_meteo_nc(i),i)=real(var1d_nc_dp(1:dim_length_meteo_nc(i)))
                    valid_dim_length_meteo_nc(i)=dim_length_meteo_nc(i)
                endif  
                !Use first file to get the valid height of the wind measurements (height3)
                if (i_file.eq.3.and.i.eq.z_dim_nc_index) then
                    val_dim_meteo_nc(1:dim_length_meteo_nc(i),i)=real(var1d_nc_dp(1:dim_length_meteo_nc(i)))
                    valid_dim_length_meteo_nc(i)=dim_length_meteo_nc(i)
                endif            
            
                !Convert from meters to km for AROME data if necessary
                if ((i.eq.x_dim_nc_index.or.i.eq.y_dim_nc_index).and.trim(unit_dim_meteo_nc(i)).eq.'km') then
                    write(unit_logfile,'(A)') 'Units of x y data are in kilometres. Converting to metres'
                    val_dim_meteo_nc(1:dim_length_meteo_nc(i),i)=val_dim_meteo_nc(1:dim_length_meteo_nc(i),i)*1000.
                    meteo_var1d_nc(1:dim_length_meteo_nc(i),i)=meteo_var1d_nc(1:dim_length_meteo_nc(i),i)*1000.
                endif
            
                write(unit_logfile,'(3A,2es12.4)') ' ',trim(dim_name_meteo_nc(i)),' (min, max): ' &
                        ,minval(meteo_var1d_nc(1:dim_length_meteo_nc(i),i)),maxval(meteo_var1d_nc(1:dim_length_meteo_nc(i),i)) 
     
            else
                !meteo_var1d_nc(1:dim_length_meteo_nc(i),i)=0.
                !val_dim_meteo_nc(1:dim_length_meteo_nc(i),i)=0.
            endif
       enddo

        
        !i_conc=compound_index
        
        !Loop through the meteo_variables
        do i=1,num_var_meteo_nc
            
            !Identify the variable name and ID in the nc file
            var_name_nc_temp=var_name_meteo_nc(i)
            status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
            !write(*,*) 'Status1: ',status_nc,var_id_nc,trim(var_name_nc_temp),i_source
            
            !If a variable name is found in the file then go further
            if (status_nc.eq.NF90_NOERR) then
                scale_factor_nc=1.
                !Find the dimensions of the variable (temp_num_dims)
                status_nc = NF90_INQUIRE_VARIABLE(id_nc, var_id_nc, ndims = temp_num_dims)
                !write(*,*) temp_num_dims,status_nc
                if (temp_num_dims.eq.2.and.i_file.eq.3) then
                    !Read latitude and longitude data into a 2d grid if available. Only lat lon is 2d?
                    if (i.eq.lat_nc_index.or.i.eq.lon_nc_index) then
                    status_nc = NF90_GET_VAR (id_nc, var_id_nc, var2d_nc_dp);meteo_var2d_nc(:,:,i)=real(var2d_nc_dp)
                    write(unit_logfile,'(A,i3,A,2A,2f16.4)') ' Reading: ',temp_num_dims,' ',trim(var_name_nc_temp),' (min, max): ',minval(meteo_var2d_nc(:,:,i)),maxval(meteo_var2d_nc(:,:,i))
                    endif
                elseif (temp_num_dims.eq.3.and.i_file.eq.4) then
                    !Special case for z0 file as they are scaled integers and a single time file
                    !write(*,'(6i)') dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),dim_start_nc(time_dim_nc_index),dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)
                    status_nc = nf90_get_att(id_nc, var_id_nc, "scale_factor", scale_factor_nc)
                    if (status_nc.ne.NF90_NOERR) scale_factor_nc=1.
                    !write(*,*) 'scale_factor=',scale_factor_nc
                    status_nc = NF90_GET_VAR (id_nc, var_id_nc, meteo_var3d_nc(:,:,1,i),start=(/dim_start_meteo_nc(x_dim_nc_index),dim_start_meteo_nc(y_dim_nc_index),dim_start_meteo_nc(time_dim_nc_index)/),count=(/dim_length_meteo_nc(x_dim_nc_index),dim_length_meteo_nc(y_dim_nc_index),dim_length_meteo_nc(time_dim_nc_index)/))
                    meteo_var3d_nc(:,:,:,i)=real(meteo_var3d_nc(:,:,:,i)*scale_factor_nc)
                    !write(*,*) status_nc
                    write(unit_logfile,'(A,I,3A,2f16.4)') ' Reading: ',temp_num_dims,' ',trim(var_name_nc_temp),' (min, max): ',minval(meteo_var3d_nc(:,:,1:dim_length_meteo_nc(time_dim_nc_index),i)),maxval(meteo_var3d_nc(:,:,1:dim_length_meteo_nc(time_dim_nc_index),i))
                elseif (temp_num_dims.eq.4.and.i_file.eq.3) then
                    !write(*,*) dim_start_nc(z_dim_nc_index),dim_start_nc(z_dim_nc_index)+dim_length_nc(z_dim_nc_index)-1
                    !write(*,*) dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),dim_start_nc(z_dim_nc_index),dim_start_nc(time_dim_nc_index)
                    !write(*,*) dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index)
                    status_nc = NF90_GET_VAR (id_nc, var_id_nc, meteo_var4d_nc(:,:,dim_start_meteo_nc(z_dim_nc_index):dim_start_meteo_nc(z_dim_nc_index)+dim_length_meteo_nc(z_dim_nc_index)-1,:,i),start=(/dim_start_meteo_nc(x_dim_nc_index),dim_start_meteo_nc(y_dim_nc_index),dim_start_meteo_nc(z_dim_nc_index),dim_start_meteo_nc(time_dim_nc_index)-1/),count=(/dim_length_meteo_nc(x_dim_nc_index),dim_length_meteo_nc(y_dim_nc_index),dim_length_meteo_nc(z_dim_nc_index),dim_length_meteo_nc(time_dim_nc_index)+1/))
                    !status_nc = NF90_GET_VAR (id_nc, var_id_nc, meteo_var4d_nc(:,:,dim_start_meteo_nc(z_dim_nc_index):dim_start_meteo_nc(z_dim_nc_index)+dim_length_meteo_nc(z_dim_nc_index)-1,dim_start_meteo_nc(time_dim_nc_index)-1:dim_start_meteo_nc(time_dim_nc_index)+dim_length_meteo_nc(time_dim_nc_index)-1,i),start=(/dim_start_meteo_nc(x_dim_nc_index),dim_start_meteo_nc(y_dim_nc_index),dim_start_meteo_nc(z_dim_nc_index),dim_start_meteo_nc(time_dim_nc_index)-1/),count=(/dim_length_meteo_nc(x_dim_nc_index),dim_length_meteo_nc(y_dim_nc_index),dim_length_meteo_nc(z_dim_nc_index),dim_length_meteo_nc(time_dim_nc_index)+1/))
                    !status_nc = NF90_GET_VAR (id_nc, var_id_nc, meteo_var4d_nc(:,:,1,dim_start_meteo_nc(time_dim_nc_index):dim_length_meteo_nc(time_dim_nc_index),i),start=(/dim_start_meteo_nc(x_dim_nc_index),dim_start_meteo_nc(y_dim_nc_index),1,dim_start_meteo_nc(time_dim_nc_index)/),count=(/dim_length_meteo_nc(x_dim_nc_index),dim_length_meteo_nc(y_dim_nc_index),1,dim_length_meteo_nc(time_dim_nc_index)/))
                    !status_nc = NF90_GET_VAR (id_nc, var_id_nc, temp_var4d_nc(:,:,:,:),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),dim_start_nc(z_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                    !var4d_nc(:,val_dim_nc:,:,:,i,i_source)=real(temp_var4d_nc(:,:,:,:))
                    write(unit_logfile,'(A,I,3A,2f16.4)') ' Reading: ',temp_num_dims,' ',trim(var_name_nc_temp),' (min, max): ',minval(meteo_var4d_nc(:,:,dim_start_meteo_nc(z_dim_nc_index):dim_start_meteo_nc(z_dim_nc_index)+dim_length_meteo_nc(z_dim_nc_index)-1,1:dim_length_meteo_nc(time_dim_nc_index),i)),maxval(meteo_var4d_nc(1:dim_length_meteo_nc(x_dim_nc_index),1:dim_length_meteo_nc(y_dim_nc_index),dim_start_meteo_nc(z_dim_nc_index):dim_start_meteo_nc(z_dim_nc_index)+dim_length_meteo_nc(z_dim_nc_index)-1,1:dim_length_meteo_nc(time_dim_nc_index),i))
                    !write(*,*) dim_start_meteo_nc(time_dim_nc_index)-1,dim_length_meteo_nc(time_dim_nc_index)+1,dim_start_meteo_nc(z_dim_nc_index),dim_start_meteo_nc(z_dim_nc_index)+dim_length_meteo_nc(z_dim_nc_index)-1
                elseif (temp_num_dims.eq.3.and.i_file.eq.3) then
                    !NBV meteo data
                    status_nc = nf90_get_att(id_nc, var_id_nc, "scale_factor", scale_factor_nc)
                    if (status_nc.ne.NF90_NOERR) scale_factor_nc=1.
                    !write(*,*) dim_start_nc(z_dim_nc_index),dim_start_nc(z_dim_nc_index)+dim_length_nc(z_dim_nc_index)-1
                    !write(*,*) dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),dim_start_nc(z_dim_nc_index),dim_start_nc(time_dim_nc_index)
                    !write(*,*) dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index)
                    !status_nc = NF90_GET_VAR (id_nc, var_id_nc, meteo_var4d_nc(:,:,1,:,i),start=(/dim_start_meteo_nc(x_dim_nc_index),dim_start_meteo_nc(y_dim_nc_index),dim_start_meteo_nc(time_dim_nc_index)-1/),count=(/dim_length_meteo_nc(x_dim_nc_index),dim_length_meteo_nc(y_dim_nc_index),dim_length_meteo_nc(time_dim_nc_index)+1/))
                    status_nc = NF90_GET_VAR (id_nc, var_id_nc, meteo_var4d_nc(:,:,1,dim_start_meteo_nc(time_dim_nc_index):dim_length_meteo_nc(time_dim_nc_index),i),start=(/dim_start_meteo_nc(x_dim_nc_index),dim_start_meteo_nc(y_dim_nc_index),dim_start_meteo_nc(time_dim_nc_index)/),count=(/dim_length_meteo_nc(x_dim_nc_index),dim_length_meteo_nc(y_dim_nc_index),dim_length_meteo_nc(time_dim_nc_index)/))
                    !status_nc = NF90_GET_VAR (id_nc, var_id_nc, temp_var4d_nc(:,:,:,:),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),dim_start_nc(z_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                    meteo_var4d_nc(:,:,:,:,i)=real(meteo_var4d_nc(:,:,:,:,i)*scale_factor_nc)
                    write(unit_logfile,'(A,I,3A,2f16.4)') ' Reading: ',temp_num_dims,' ',trim(var_name_nc_temp),' (min, max): ',minval(meteo_var4d_nc(:,:,1,1:dim_length_meteo_nc(time_dim_nc_index),i)),maxval(meteo_var4d_nc(1:dim_length_meteo_nc(x_dim_nc_index),1:dim_length_meteo_nc(y_dim_nc_index),1,1:dim_length_meteo_nc(time_dim_nc_index),i))
                    !write(*,*) dim_start_meteo_nc(time_dim_nc_index)-1,dim_length_meteo_nc(time_dim_nc_index)+1,dim_start_meteo_nc(z_dim_nc_index),dim_start_meteo_nc(z_dim_nc_index)+dim_length_meteo_nc(z_dim_nc_index)-1
                elseif (temp_num_dims.eq.5.and.i_file.eq.3) then
                    !This is the case when there is an ensemble member in the format
                    !write(*,*) dim_start_nc(z_dim_nc_index),dim_start_nc(z_dim_nc_index)+dim_length_nc(z_dim_nc_index)-1
                    !write(*,*) dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),dim_start_nc(z_dim_nc_index),dim_start_nc(time_dim_nc_index)
                    !write(*,*) dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index)
                    status_nc = NF90_GET_VAR (id_nc, var_id_nc, meteo_var4d_nc(:,:,dim_start_meteo_nc(z_dim_nc_index):dim_start_meteo_nc(z_dim_nc_index)+dim_length_meteo_nc(z_dim_nc_index)-1,:,i),start=(/dim_start_meteo_nc(x_dim_nc_index),dim_start_meteo_nc(y_dim_nc_index),1,dim_start_meteo_nc(z_dim_nc_index),dim_start_meteo_nc(time_dim_nc_index)-1/),count=(/dim_length_meteo_nc(x_dim_nc_index),dim_length_meteo_nc(y_dim_nc_index),1,dim_length_meteo_nc(z_dim_nc_index),dim_length_meteo_nc(time_dim_nc_index)+1/))
                    !status_nc = NF90_GET_VAR (id_nc, var_id_nc, temp_var4d_nc(:,:,:,:),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),dim_start_nc(z_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                    !var4d_nc(:,val_dim_nc:,:,:,i,i_source)=real(temp_var4d_nc(:,:,:,:))
                    write(unit_logfile,'(A,I,3A,2f16.4)') ' Reading: ',temp_num_dims,' ',trim(var_name_nc_temp),' (min, max): ',minval(meteo_var4d_nc(:,:,dim_start_meteo_nc(z_dim_nc_index):dim_start_meteo_nc(z_dim_nc_index)+dim_length_meteo_nc(z_dim_nc_index)-1,1:dim_length_meteo_nc(time_dim_nc_index),i)),maxval(meteo_var4d_nc(1:dim_length_meteo_nc(x_dim_nc_index),1:dim_length_meteo_nc(y_dim_nc_index),dim_start_meteo_nc(z_dim_nc_index):dim_start_meteo_nc(z_dim_nc_index)+dim_length_meteo_nc(z_dim_nc_index)-1,1:dim_length_meteo_nc(time_dim_nc_index),i))
                    !write(*,*) dim_start_meteo_nc(time_dim_nc_index)-1,dim_length_meteo_nc(time_dim_nc_index)+1,dim_start_meteo_nc(z_dim_nc_index),dim_start_meteo_nc(z_dim_nc_index)+dim_length_meteo_nc(z_dim_nc_index)-1
                else
                    write(unit_logfile,'(8A,8A)') ' Cannot find a correct dimension for: ',trim(var_name_nc_temp)
                endif    
                
            else
                 !write(unit_logfile,'(8A,8A)') ' Cannot read: ',trim(var_name_nc_temp)
            endif

        enddo
                    
        !Read in 2m temperature completely to get the daily average for home heating
        !if (use_RWC_emission_data.and.save_emissions_for_EMEP(heating_index).and.i_file.eq.3) then
        if (use_RWC_emission_data.and.i_file.eq.3) then
            DMT_start_time_nc_index=start_time_meteo_nc_index
            DMT_end_time_nc_index=end_time_meteo_nc_index
            !DMT_start_time_nc_index=save_emissions_start_index
            !DMT_end_time_nc_index=save_emissions_end_index
            
            DMT_dim_length_nc=DMT_end_time_nc_index-DMT_start_time_nc_index+1
            if (allocated(DMT_EMEP_grid_nc)) deallocate (DMT_EMEP_grid_nc)
            if (.not.allocated(DMT_EMEP_grid_nc)) allocate (DMT_EMEP_grid_nc(dim_length_meteo_nc(x_dim_nc_index),dim_length_meteo_nc(y_dim_nc_index),DMT_dim_length_nc))
            
            !write(*,*) DMT_start_time_nc_index,DMT_end_time_nc_index,DMT_dim_length_nc
            if (calculate_source(heating_index).and.i_file.eq.3) then
                var_name_nc_temp=var_name_meteo_nc(t2m_nc_index)
                status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                status_nc = NF90_INQUIRE_VARIABLE(id_nc, var_id_nc, ndims = temp_num_dims)
                if (temp_num_dims.eq.4) then
                    status_nc = NF90_GET_VAR (id_nc, var_id_nc, DMT_EMEP_grid_nc(:,:,:),start=(/dim_start_meteo_nc(x_dim_nc_index),dim_start_meteo_nc(y_dim_nc_index),1,DMT_start_time_nc_index/),count=(/dim_length_meteo_nc(x_dim_nc_index),dim_length_meteo_nc(y_dim_nc_index),1,DMT_dim_length_nc/))
                elseif (temp_num_dims.eq.3) then
                    !NBV meteo data
                    status_nc = NF90_GET_VAR (id_nc, var_id_nc, DMT_EMEP_grid_nc(:,:,:),start=(/dim_start_meteo_nc(x_dim_nc_index),dim_start_meteo_nc(y_dim_nc_index),DMT_start_time_nc_index/),count=(/dim_length_meteo_nc(x_dim_nc_index),dim_length_meteo_nc(y_dim_nc_index),DMT_dim_length_nc/))
                elseif (temp_num_dims.eq.5) then
                    status_nc = NF90_GET_VAR (id_nc, var_id_nc, DMT_EMEP_grid_nc(:,:,:),start=(/dim_start_meteo_nc(x_dim_nc_index),dim_start_meteo_nc(y_dim_nc_index),1,1,DMT_start_time_nc_index/),count=(/dim_length_meteo_nc(x_dim_nc_index),dim_length_meteo_nc(y_dim_nc_index),1,1,DMT_dim_length_nc/))
                else
                    write(unit_logfile,'(8A,8A)') ' Cannot find a correct dimension for: ',trim(var_name_nc_temp)
                endif                    
                
                write(unit_logfile,'(A,i,3A,2f16.4)') ' Reading: ',temp_num_dims,' ',trim(var_name_nc_temp),' (min, max): ',minval(DMT_EMEP_grid_nc),maxval(DMT_EMEP_grid_nc)
                DMT_EMEP_grid_nc(:,:,1)=sum(DMT_EMEP_grid_nc,3)/DMT_dim_length_nc-273.13
                write(unit_logfile,'(3A,2f16.4)') ' Calculating mean: ',trim('Daily mean temperature'),' (min, max): ',minval(DMT_EMEP_grid_nc(:,:,1)),maxval(DMT_EMEP_grid_nc(:,:,1))
            
            endif
        
        endif

        status_nc = NF90_CLOSE (id_nc)
                    
    enddo !End file loop
    
    !Set the correct time dimensions to the first file value
    dim_length_meteo_nc=valid_dim_length_meteo_nc
        
    !Set the grid spacing
    if (meteo_nc_projection_type.eq.LL_projection_index) then
        meteo_dgrid_nc(lon_nc_index)=meteo_var1d_nc(2,x_dim_nc_index)-meteo_var1d_nc(1,x_dim_nc_index)
        meteo_dgrid_nc(lat_nc_index)=meteo_var1d_nc(2,y_dim_nc_index)-meteo_var1d_nc(1,y_dim_nc_index)
        write(unit_logfile,'(A,2f16.4)') ' Grid spacing meteo (lon,lat): ',meteo_dgrid_nc(lon_nc_index),meteo_dgrid_nc(lat_nc_index)
    else       
        meteo_dgrid_nc(lon_nc_index)=meteo_var1d_nc(2,x_dim_nc_index)-meteo_var1d_nc(1,x_dim_nc_index)
        meteo_dgrid_nc(lat_nc_index)=meteo_var1d_nc(2,y_dim_nc_index)-meteo_var1d_nc(1,y_dim_nc_index)
        write(unit_logfile,'(A,2f16.4)') ' Grid spacing meteo (x,y) in meters: ',meteo_dgrid_nc(lon_nc_index),meteo_dgrid_nc(lat_nc_index)
    endif
        

        !Do manipulations when the additional meteorology is read in
        !Put everything in 3d data since it is all surface values
            write(unit_logfile,'(A)') ' Calculating alternative meteorological data'
            write(unit_logfile,'(A,4i)') ' Dimensions: ',dim_length_meteo_nc
            !logz0 is read in as z0 and must be converted to logz0 
            do t=1,dim_length_meteo_nc(time_dim_nc_index)
                meteo_var3d_nc(:,:,t,logz0_nc_index)=meteo_var3d_nc(:,:,1,logz0_nc_index)
                !write(*,*) t,sum(meteo_var3d_nc(:,:,t,logz0_nc_index))/dim_length_meteo_nc(x_dim_nc_index)/dim_length_meteo_nc(y_dim_nc_index)
            enddo
            where (meteo_var3d_nc(:,:,:,logz0_nc_index).lt.0.001 ) meteo_var3d_nc(:,:,:,logz0_nc_index)=0.001 
            meteo_var3d_nc(:,:,:,logz0_nc_index)=log(meteo_var3d_nc(:,:,:,logz0_nc_index))
            
            meteo_var3d_nc(:,:,:,t2m_nc_index)=meteo_var4d_nc(:,:,surface_level_nc,:,t2m_nc_index)
            !Assumes that the first step is 0 and data is hourly. So that the start time step for meteo must correspond to the second hour of any calculation
            !t=1
            !meteo_var3d_nc(:,:,t,Hflux_nc_index)=meteo_var4d_nc(:,:,surface_level_nc,t,Hflux_nc_index)/3600.
            !meteo_var3d_nc(:,:,t,uw_nc_index)=meteo_var4d_nc(:,:,surface_level_nc,t,uw_nc_index)/3600.
            !meteo_var3d_nc(:,:,t,vw_nc_index)=meteo_var4d_nc(:,:,surface_level_nc,t,vw_nc_index)/3600.
            if (index(alternative_meteorology_type,'nbv').gt.0) then
                !Flux is upward and not integrated. Stress is not integrated
                meteo_var3d_nc(:,:,:,Hflux_nc_index)=-meteo_var4d_nc(:,:,surface_level_nc,:,Hflux_nc_index)
                meteo_var3d_nc(:,:,:,uw_nc_index)=meteo_var4d_nc(:,:,surface_level_nc,:,uw_nc_index)
                meteo_var3d_nc(:,:,:,vw_nc_index)=meteo_var4d_nc(:,:,surface_level_nc,:,vw_nc_index)
                meteo_var3d_nc(:,:,:,precip_nc_index)=meteo_var4d_nc(:,:,surface_level_nc,:,precip_nc_index)
            else
                do t=dim_length_meteo_nc(time_dim_nc_index),1,-1
                meteo_var3d_nc(:,:,t,Hflux_nc_index)=(meteo_var4d_nc(:,:,surface_level_nc,t,Hflux_nc_index)-meteo_var4d_nc(:,:,surface_level_nc,t-1,Hflux_nc_index))/3600.
                meteo_var3d_nc(:,:,t,uw_nc_index)=(meteo_var4d_nc(:,:,surface_level_nc,t,uw_nc_index)-meteo_var4d_nc(:,:,surface_level_nc,t-1,uw_nc_index))/3600.
                meteo_var3d_nc(:,:,t,vw_nc_index)=(meteo_var4d_nc(:,:,surface_level_nc,t,vw_nc_index)-meteo_var4d_nc(:,:,surface_level_nc,t-1,vw_nc_index))/3600.
                meteo_var3d_nc(:,:,t,precip_nc_index)=(meteo_var4d_nc(:,:,surface_level_nc,t,precip_nc_index)-meteo_var4d_nc(:,:,surface_level_nc,t-1,precip_nc_index))
                !write(*,*) t,sum(meteo_var3d_nc(:,:,t,Hflux_nc_index))/dim_length_meteo_nc(x_dim_nc_index)/dim_length_meteo_nc(y_dim_nc_index) &
                !    ,sum(meteo_var4d_nc(:,:,surface_level_nc,t,hmix_nc_index))/dim_length_meteo_nc(x_dim_nc_index)/dim_length_meteo_nc(y_dim_nc_index)
                enddo
            endif
            
            !Approximate density of air used (1.2 kg/m^3 +/- 10%)
            meteo_var3d_nc(:,:,:,ustar_nc_index)=sqrt(sqrt(meteo_var3d_nc(:,:,:,uw_nc_index)**2+meteo_var3d_nc(:,:,:,vw_nc_index)**2)/1.2)
            where (meteo_var3d_nc(:,:,:,ustar_nc_index).lt.ustar_min) meteo_var3d_nc(:,:,:,ustar_nc_index)=ustar_min
            !do t=dim_length_meteo_nc(time_dim_nc_index),0,-1
            !     write(*,*) t,sum(meteo_var3d_nc(:,:,t,ustar_nc_index))/dim_length_meteo_nc(x_dim_nc_index)/dim_length_meteo_nc(y_dim_nc_index) &
            !        ,sum(meteo_var3d_nc(:,:,t,uw_nc_index))/dim_length_meteo_nc(x_dim_nc_index)/dim_length_meteo_nc(y_dim_nc_index) &
            !        ,sum(meteo_var3d_nc(:,:,t,vw_nc_index))/dim_length_meteo_nc(x_dim_nc_index)/dim_length_meteo_nc(y_dim_nc_index) 
            !enddo
             
            !Approximate temperature (273 +/- 10%)
            !Have inserted the correct temperature now
            meteo_var3d_nc(:,:,:,invL_nc_index)=meteo_var3d_nc(:,:,:,Hflux_nc_index)*0.4*9.8/1.2/1004./meteo_var3d_nc(:,:,:,t2m_nc_index)/meteo_var3d_nc(:,:,:,ustar_nc_index)**3
            !Limit stable L to lowest_stable_L and to lowest_unstable_L (negative number) for unstable.
            where (meteo_var3d_nc(:,:,:,invL_nc_index).lt.1.0/lowest_unstable_L) meteo_var3d_nc(:,:,:,invL_nc_index)=1.0/lowest_unstable_L
            where (meteo_var3d_nc(:,:,:,invL_nc_index).gt.1.0/lowest_stable_L) meteo_var3d_nc(:,:,:,invL_nc_index)=1.0/lowest_stable_L
            
            !Put the 10 m wind vectors as the lowest grid level
            !H_meteo=val_dim_meteo_nc(surface_level_nc,z_dim_nc_index)
            H_meteo=10.
            write(unit_logfile,'(A,f8.2)')' Alternative meteo: setting lowest meteo grid height = ',H_meteo
            meteo_var3d_nc(:,:,:,ugrid_nc_index)=meteo_var4d_nc(:,:,surface_level_nc,:,u10_nc_index)
            meteo_var3d_nc(:,:,:,vgrid_nc_index)=meteo_var4d_nc(:,:,surface_level_nc,:,v10_nc_index)
            meteo_var3d_nc(:,:,:,FF10_nc_index)=meteo_var4d_nc(:,:,surface_level_nc,:,FF10_nc_index)
            if (sum(abs(meteo_var3d_nc(:,:,:,FF10_nc_index))).eq.0) then
                !Calculate wind speed if it can't read it
                meteo_var3d_nc(:,:,:,FF10_nc_index)=sqrt(meteo_var3d_nc(:,:,:,ugrid_nc_index)**2+meteo_var3d_nc(:,:,:,vgrid_nc_index)**2)
            else
                meteo_var3d_nc(:,:,:,FF10_nc_index)=sqrt(meteo_var3d_nc(:,:,:,u10_nc_index)**2+meteo_var3d_nc(:,:,:,v10_nc_index)**2)
            endif
            
            !Smooth the boundary layer height (running mean) and set minimum 
            meteo_var3d_nc(:,:,:,hmix_nc_index)=0.
            do j=2,dim_length_meteo_nc(y_dim_nc_index)-1
            do i=2,dim_length_meteo_nc(x_dim_nc_index)-1
                do jj=-1,1
                do ii=-1,1
                    meteo_var3d_nc(i,j,:,hmix_nc_index)=meteo_var3d_nc(i,j,:,hmix_nc_index)+meteo_var4d_nc(i+ii,j+jj,surface_level_nc,:,hmix_nc_index)/9.
                    !if (ii.ne.0.and.jj.ne.0) then
                    !    var3d_nc(i,j,:,hmix_nc_index,allsource_index)=var3d_nc(i,j,:,hmix_nc_index,allsource_index)-var3d_nc(i+ii,j+jj,:,hmix_nc_index,traffic_index)/8.
                    !endif
                enddo
                enddo
                 
            enddo
            enddo
            where (meteo_var3d_nc(:,:,:,hmix_nc_index).lt.hmix_min) meteo_var3d_nc(:,:,:,hmix_nc_index)=hmix_min
            where (meteo_var3d_nc(:,:,:,hmix_nc_index).gt.hmix_max) meteo_var3d_nc(:,:,:,hmix_nc_index)=hmix_max
            where (meteo_var3d_nc(:,:,:,ustar_nc_index).lt.ustar_min) meteo_var3d_nc(:,:,:,ustar_nc_index)=ustar_min
            !do t=dim_length_nc(time_dim_nc_index),2,-1
            !    write(*,*) t,sum(var3d_nc(:,:,t,hmix_nc_index,allsource_nc_index))/dim_length_nc(x_dim_nc_index)/dim_length_nc(y_dim_nc_index)
            !enddo
            
            write(unit_logfile,'(3A,2f16.4)') ' Alternative meteo: ',trim(var_name_meteo_nc(logz0_nc_index)),' (min, max): ',minval(meteo_var3d_nc(:,:,1:dim_length_meteo_nc(time_dim_nc_index),logz0_nc_index)),maxval(meteo_var3d_nc(:,:,1:dim_length_meteo_nc(time_dim_nc_index),logz0_nc_index))
            write(unit_logfile,'(3A,2f16.4)') ' Alternative meteo: ',trim(var_name_meteo_nc(ustar_nc_index)),' (min, max): ',minval(meteo_var3d_nc(:,:,1:dim_length_meteo_nc(time_dim_nc_index),ustar_nc_index)),maxval(meteo_var3d_nc(:,:,1:dim_length_meteo_nc(time_dim_nc_index),ustar_nc_index))
            write(unit_logfile,'(3A,2f16.4)') ' Alternative meteo: ',trim(var_name_meteo_nc(Hflux_nc_index)),' (min, max): ',minval(meteo_var3d_nc(:,:,1:dim_length_meteo_nc(time_dim_nc_index),Hflux_nc_index)),maxval(meteo_var3d_nc(:,:,1:dim_length_meteo_nc(time_dim_nc_index),Hflux_nc_index))
            write(unit_logfile,'(3A,2f16.4)') ' Alternative meteo: ',trim(var_name_meteo_nc(invL_nc_index)),' (min, max): ',minval(meteo_var3d_nc(:,:,1:dim_length_meteo_nc(time_dim_nc_index),invL_nc_index)),maxval(meteo_var3d_nc(:,:,1:dim_length_meteo_nc(time_dim_nc_index),invL_nc_index))
            write(unit_logfile,'(3A,2f16.4)') ' Alternative meteo: ',trim(var_name_meteo_nc(ugrid_nc_index)),' (min, max): ',minval(meteo_var3d_nc(:,:,1:dim_length_meteo_nc(time_dim_nc_index),ugrid_nc_index)),maxval(meteo_var3d_nc(:,:,1:dim_length_meteo_nc(time_dim_nc_index),ugrid_nc_index))
            write(unit_logfile,'(3A,2f16.4)') ' Alternative meteo: ',trim(var_name_meteo_nc(vgrid_nc_index)),' (min, max): ',minval(meteo_var3d_nc(:,:,1:dim_length_meteo_nc(time_dim_nc_index),vgrid_nc_index)),maxval(meteo_var3d_nc(:,:,1:dim_length_meteo_nc(time_dim_nc_index),vgrid_nc_index))
            write(unit_logfile,'(3A,2f16.4)') ' Alternative meteo: ',trim(var_name_meteo_nc(FF10_nc_index)),' (min, max): ',minval(meteo_var3d_nc(:,:,1:dim_length_meteo_nc(time_dim_nc_index),FF10_nc_index)),maxval(meteo_var3d_nc(:,:,1:dim_length_meteo_nc(time_dim_nc_index),FF10_nc_index))
            write(unit_logfile,'(3A,2f16.4)') ' Alternative meteo: ',trim(var_name_meteo_nc(hmix_nc_index)),' (min, max): ',minval(meteo_var3d_nc(:,:,1:dim_length_meteo_nc(time_dim_nc_index),hmix_nc_index)),maxval(meteo_var3d_nc(:,:,1:dim_length_meteo_nc(time_dim_nc_index),hmix_nc_index))
            write(unit_logfile,'(3A,2f16.4)') ' Alternative meteo: ',trim(var_name_meteo_nc(t2m_nc_index)),' (min, max): ',minval(meteo_var3d_nc(:,:,1:dim_length_meteo_nc(time_dim_nc_index),t2m_nc_index)),maxval(meteo_var3d_nc(:,:,1:dim_length_meteo_nc(time_dim_nc_index),t2m_nc_index))
            
            
        !If no logz0 available. Set to log(0.1)
        !For urban areas a value of 0.3 is used
        !where (var3d_nc(:,:,:,logz0_nc_index,:).eq.0.0) var3d_nc(:,:,:,logz0_nc_index,:)=log(0.3)
        if (replace_z0.ne.NODATA_value) then
            write(unit_logfile,'(A,f8.4)') ' Replacing z0 everywhere with: ',replace_z0
            meteo_var3d_nc(:,:,:,logz0_nc_index)=log(replace_z0)                
        endif
        if (replace_invL.ne.NODATA_value) then
            write(unit_logfile,'(A,f8.4)') ' Replacing inverse L everywhere with: ',replace_invL
            meteo_var3d_nc(:,:,:,invL_nc_index)=replace_invL                
        endif
        if (replace_hmix.ne.NODATA_value) then
            write(unit_logfile,'(A,f8.4)') ' Replacing HMIX everywhere with: ',replace_hmix
            meteo_var3d_nc(:,:,:,hmix_nc_index)=replace_hmix
        endif
        if (FF_scale.ne.NODATA_value) then
            write(unit_logfile,'(A,f8.4)') ' Rescaling wind fields everywhere with factor: ',FF_scale
            meteo_var3d_nc(:,:,:,ustar_nc_index)=meteo_var3d_nc(:,:,:,ustar_nc_index)*FF_scale
            meteo_var3d_nc(:,:,:,FF10_nc_index)=meteo_var3d_nc(:,:,:,FF10_nc_index)*FF_scale
            meteo_var3d_nc(:,:,:,inv_FF10_nc_index)=meteo_var3d_nc(:,:,:,inv_FF10_nc_index)/FF_scale
            meteo_var3d_nc(:,:,:,ugrid_nc_index)=meteo_var3d_nc(:,:,:,ugrid_nc_index)*FF_scale
            meteo_var3d_nc(:,:,:,vgrid_nc_index)=meteo_var3d_nc(:,:,:,vgrid_nc_index)*FF_scale
            meteo_var3d_nc(:,:,:,inv_FFgrid_nc_index)=meteo_var3d_nc(:,:,:,inv_FFgrid_nc_index)/FF_scale
            meteo_var3d_nc(:,:,:,ugrid_nc_index)=meteo_var3d_nc(:,:,:,ugrid_nc_index)*FF_scale
            meteo_var3d_nc(:,:,:,vgrid_nc_index)=meteo_var3d_nc(:,:,:,vgrid_nc_index)*FF_scale
        endif
        if (FF10_offset.ne.NODATA_value) then
            write(unit_logfile,'(A,f8.4)') ' Offsetting 10 m wind fields everywhere with a value: ',FF10_offset
            meteo_var3d_nc(:,:,:,FF10_nc_index)=meteo_var3d_nc(:,:,:,FF10_nc_index)+FF10_offset
        endif
        if (DD_offset.ne.NODATA_value) then
            write(unit_logfile,'(A,f8.4)') ' Rotating wind fields everywhere with a value: ',DD_offset

            if (.not.allocated(temp_meteo_var3d_nc)) allocate (temp_meteo_var3d_nc(valid_dim_length_meteo_nc(x_dim_nc_index),valid_dim_length_meteo_nc(y_dim_nc_index),valid_dim_length_meteo_nc(time_dim_nc_index),2))
            temp_meteo_var3d_nc(:,:,:,1) = meteo_var3d_nc(:,:,:,ugrid_nc_index)*cos(DD_offset/180.*3.14159)+meteo_var3d_nc(:,:,:,vgrid_nc_index)*sin(DD_offset/180.*3.14159)                                       
            temp_meteo_var3d_nc(:,:,:,2) =-meteo_var3d_nc(:,:,:,ugrid_nc_index)*sin(DD_offset/180.*3.14159)+meteo_var3d_nc(:,:,:,vgrid_nc_index)*cos(DD_offset/180.*3.14159)                                       
            meteo_var3d_nc(:,:,:,ugrid_nc_index) = temp_meteo_var3d_nc(:,:,:,1)
            meteo_var3d_nc(:,:,:,vgrid_nc_index) = temp_meteo_var3d_nc(:,:,:,2)
        endif
       
        !Set the magnitude of the gridded wind fields. Should probably be done after subgridding?
        meteo_var3d_nc(:,:,:,FFgrid_nc_index)=sqrt(meteo_var3d_nc(:,:,:,ugrid_nc_index)**2+meteo_var3d_nc(:,:,:,vgrid_nc_index)**2)
  
        !Check meteo time which comes in seconds since 1970. Converts to days.
        if (use_alternative_meteorology_flag) then
            date_num_temp=val_dim_meteo_nc(1,time_dim_nc_index)/3600./24.+30./24./3600.
            call number_to_date(date_num_temp,date_array,ref_year_meteo)
            write(unit_logfile,'(a,i6)') ' Time dimension meteo: ',dim_length_meteo_nc(time_dim_nc_index)
            write(unit_logfile,'(a,6i6)') ' Date start meteo = ',date_array
            !date_num_temp=dble(ceiling(val_dim_meteo_nc(dim_length_meteo_nc(time_dim_nc_index),time_dim_nc_index)/3600.))/24.
            date_num_temp=val_dim_meteo_nc(dim_length_meteo_nc(time_dim_nc_index),time_dim_nc_index)/3600./24.+30./24./3600.
            call number_to_date(date_num_temp,date_array,ref_year_meteo)
            write(unit_logfile,'(a,6i6)') ' Date end meteo =   ',date_array
            !write(*,*) start_time_meteo_nc_index,valid_dim_length_meteo_nc(time_dim_nc_index)
        endif
        !do t=1,dim_length_meteo_nc(time_dim_nc_index)
        !    date_num_temp=val_dim_meteo_nc(t,time_dim_nc_index)/3600./24.+0.1/24.
        !    call number_to_date(date_num_temp,date_array,ref_year_meteo)
        !    write(unit_logfile,'(a,i4,6i6,f12.4)') ' Date end meteo =   ',t,date_array,date_num_temp
        !enddo
        !stop
        
        if (allocated(var1d_nc_dp)) deallocate (var1d_nc_dp)
        if (allocated(var2d_nc_dp)) deallocate (var2d_nc_dp)
        if (allocated(temp_meteo_var3d_nc)) deallocate (temp_meteo_var3d_nc)
        if (allocated(meteo_var4d_nc)) deallocate (meteo_var4d_nc) 
 
    end subroutine uEMEP_read_meteo_nc
    
    