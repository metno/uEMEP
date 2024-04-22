!uEMEP_read_EMEP
    
    subroutine uEMEP_read_EMEP
    
    use uEMEP_definitions
    use netcdf
    
    implicit none
    
    integer i,j,k,t
    integer ii,jj,iii,jjj,iiii,jjjj
    logical exists
    character(256) pathfilename_nc
    integer status_nc     !Error message
    integer id_nc
    integer dim_id_nc(num_dims_nc)
    character(256) dimname_temp,var_name_nc_temp,var_name_nc_temp2,unit_name_nc_temp
    integer var_id_nc
    real :: local_fraction_scaling=1.0
    integer i_file,i_source,i_conc,i_dim
    integer temp_frac_index,temp_file_index,temp_compound_index,temp_source_index
    integer temp_num_dims
    integer temp_start_time_nc_index,temp_end_time_nc_index
    integer i_loop
    integer valid_dim_length_nc(num_dims_nc) !dimensions of file 1
    integer numAtts_projection
    integer surface_level_nc_2
    
    real temp_lat(4),temp_lon(4)
    real temp_y(4),temp_x(4)
    real temp_x_min,temp_x_max,temp_y_min,temp_y_max
    integer i_temp_min,i_temp_max,j_temp_min,j_temp_max
    double precision temp_var1d_nc_dp(2,2)
    real temp_delta(2)
    real H_emep_temp
    
    integer n_file
    double precision date_num_temp,date_num_2000
    integer date_array(6)
    double precision scale_factor_nc
    
    integer i_pollutant,p_loop,p_loop_index
    
    integer DMT_start_time_nc_index,DMT_end_time_nc_index,DMT_dim_length_nc

    logical nonzero_wind_notfound

    integer i_sp,ii_sp,pmxx_sp_index
    
    integer i_depo
    
    !Temporary reading variables
    double precision, allocatable :: var1d_nc_dp(:)
    double precision, allocatable :: var2d_nc_dp(:,:)
    double precision, allocatable :: var3d_nc_dp(:,:,:)
    double precision, allocatable :: var4d_nc_dp(:,:,:,:)
    
    !Temporary files for rotating wind field and PM
    real, allocatable :: temp_var4d_nc(:,:,:,:,:)
    
    real, allocatable :: temp_var3d_nc(:,:,:)

    !Temporary PM arrays for reading in PM10
    real, allocatable :: pm_var4d_nc(:,:,:,:,:,:,:)
    real, allocatable :: pm_var3d_nc(:,:,:,:,:,:)
    real, allocatable :: pm_lc_var4d_nc(:,:,:,:,:,:,:,:,:)
    
    real, allocatable :: species_temp_var3d_nc(:,:,:)
    
    
    !NOTE: temporary for nh3 false is not on
    logical :: use_comp_temporary=.false.
    logical :: EMEP_region_outside_domain=.false.
    
    real EMEP_grid_interpolation_size_temp
    
    real mean_phi_temp,mean_invL_temp
    integer phi_count

    !Functions
    double precision date_to_number
    
    write(unit_logfile,'(A)') ''
    write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Reading EMEP data (uEMEP_read_EMEP)'
	write(unit_logfile,'(A)') '================================================================'

    
    !This if statement is already specified in uEMEP_define_subgrid and is not necessary here
    if (hourly_calculations) then
        temp_start_time_nc_index=start_time_nc_index
        temp_end_time_nc_index=end_time_nc_index
    else
        temp_start_time_nc_index=1
        temp_end_time_nc_index=1
    endif
     
    if (use_single_time_loop_flag) then
        temp_start_time_nc_index=start_time_nc_index+t_loop-1
        temp_end_time_nc_index=temp_start_time_nc_index
    endif

    !Presettng the surface level to 1. Valid when there is no inverting of layers
    surface_level_nc=EMEP_surface_level_nc
    surface_level_nc_2=EMEP_surface_level_nc_2
    write(unit_logfile,'(A,I)') ' Surface level base set to: ',surface_level_nc
    write(unit_logfile,'(A,I)') ' Surface level local_contribution set to: ',surface_level_nc_2

        if (allocated(val_dim_nc)) deallocate (val_dim_nc)
        if (allocated(unit_dim_nc)) deallocate (unit_dim_nc)
        if (allocated(var1d_nc)) deallocate (var1d_nc)
        if (allocated(var2d_nc)) deallocate (var2d_nc)
        if (allocated(var3d_nc)) deallocate (var3d_nc)
        if (allocated(var4d_nc)) deallocate (var4d_nc)
        if (allocated(comp_var3d_nc)) deallocate (comp_var3d_nc)
        if (allocated(comp_var4d_nc)) deallocate (comp_var4d_nc)
        if (allocated(var1d_nc_dp)) deallocate (var1d_nc_dp) 
        if (allocated(var2d_nc_dp)) deallocate (var2d_nc_dp)
        if (allocated(lc_var3d_nc)) deallocate (lc_var3d_nc)
        if (allocated(lc_var4d_nc)) deallocate (lc_var4d_nc)
        if (allocated(DMT_EMEP_grid_nc)) deallocate (DMT_EMEP_grid_nc) !Daily mean temperature
        if (allocated(species_var3d_nc)) deallocate (species_var3d_nc)
        if (allocated(depo_var3d_nc)) deallocate (depo_var3d_nc)
        if (allocated(temp_var3d_nc)) deallocate (temp_var3d_nc)
      
    !Loop through the EMEP files containing the data

    n_file=2
    do i_file=1,n_file
        

        !Temporary fix. Must remove
        !if (i_file.eq.2) dim_name_nc(z_dim_nc_index)='klevel'
    
        !Set the filename
        pathfilename_EMEP(i_file)=trim(pathname_EMEP(i_file))//trim(filename_EMEP(i_file))
     
        !Test existence of the filename. If does not exist then stop
        inquire(file=trim(pathfilename_EMEP(i_file)),exist=exists)
        if (.not.exists) then
            write(unit_logfile,'(A,A)') ' ERROR: Netcdf file does not exist: ', trim(pathfilename_EMEP(i_file))
            write(unit_logfile,'(A)') '  STOPPING'
            stop
        endif

        !Open the netcdf file for reading
        write(unit_logfile,'(2A)') ' Opening netcdf file: ',trim(pathfilename_EMEP(i_file))
        status_nc = NF90_OPEN (pathfilename_EMEP(i_file), nf90_nowrite, id_nc)
        if (status_nc .NE. NF90_NOERR) write(unit_logfile,'(A,I)') 'ERROR opening netcdf file: ',status_nc
    
        EMEP_projection_type=LL_projection_index

        !Find the projection. If no projection then in lat lon coordinates
        status_nc = NF90_INQ_VARID (id_nc,'projection_lambert',var_id_nc)
        
        if (status_nc.eq.NF90_NOERR) then
            !If there is a projection then read in the attributes. All these are doubles
            !status_nc = nf90_inquire_variable(id_nc, var_id_nc, natts = numAtts_projection)
            status_nc = nf90_get_att(id_nc, var_id_nc, 'standard_parallel', EMEP_projection_attributes(1:2))
            status_nc = nf90_get_att(id_nc, var_id_nc, 'longitude_of_central_meridian', EMEP_projection_attributes(3))
            status_nc = nf90_get_att(id_nc, var_id_nc, 'latitude_of_projection_origin', EMEP_projection_attributes(4))
            status_nc = nf90_get_att(id_nc, var_id_nc, 'earth_radius', EMEP_projection_attributes(5))
                
            EMEP_projection_type=LCC_projection_index

                !Reset names of the x,y coordinates
            dim_name_nc(x_dim_nc_index)='i'
            dim_name_nc(y_dim_nc_index)='j'
            var_name_nc(lon_nc_index,:,allsource_index)='lon'
            var_name_nc(lat_nc_index,:,allsource_index)='lat'
            
            write(unit_logfile,'(A,5f12.2)') 'Reading lambert_conformal_conic projection. ',EMEP_projection_attributes(1:5)
            if (EMEP_projection_attributes(1).ne.EMEP_projection_attributes(4).or.EMEP_projection_attributes(2).ne.EMEP_projection_attributes(4)) then
                use_alternative_LCC_projection_flag=.true.
                write(unit_logfile,'(A,l)') 'Using alternative lambert_conformal_conic projection: ',use_alternative_LCC_projection_flag
            else
                use_alternative_LCC_projection_flag=.false.                
            endif
            !Always set to true
            use_alternative_LCC_projection_flag=.true.
        endif
        
        !Find the projection. If no projection then in lat lon coordinates
        status_nc = NF90_INQ_VARID (id_nc,'Polar_Stereographic',var_id_nc)
        
        if (status_nc.eq.NF90_NOERR) then
            !If there is a projection then read in the attributes. All these are doubles
            !status_nc = nf90_inquire_variable(id_nc, var_id_nc, natts = numAtts_projection)
            EMEP_projection_attributes=0.
            status_nc = nf90_get_att(id_nc, var_id_nc, 'straight_vertical_longitude_from_pole', EMEP_projection_attributes(1))
            status_nc = nf90_get_att(id_nc, var_id_nc, 'latitude_of_projection_origin', EMEP_projection_attributes(2))
            status_nc = nf90_get_att(id_nc, var_id_nc, 'false_easting', EMEP_projection_attributes(3))
            if (status_nc.ne.NF90_NOERR) EMEP_projection_attributes(3)=0.
            status_nc = nf90_get_att(id_nc, var_id_nc, 'false_northing', EMEP_projection_attributes(4))
            if (status_nc.ne.NF90_NOERR) EMEP_projection_attributes(4)=0.
            status_nc = nf90_get_att(id_nc, var_id_nc, 'earth_radius', EMEP_projection_attributes(5))
            if (status_nc.ne.NF90_NOERR) EMEP_projection_attributes(5)=6.370e6
            status_nc = nf90_get_att(id_nc, var_id_nc, 'scale_factor_at_projection_origin', EMEP_projection_attributes(6))
            if (status_nc.ne.NF90_NOERR) EMEP_projection_attributes(6)=1.
                
            EMEP_projection_type=PS_projection_index

            !Reset names of the x,y coordinates
            dim_name_nc(x_dim_nc_index)='i'
            dim_name_nc(y_dim_nc_index)='j'
            var_name_nc(lon_nc_index,:,allsource_index)='lon'
            var_name_nc(lat_nc_index,:,allsource_index)='lat'
            
            write(unit_logfile,'(A,6f12.2)') 'Reading Polar_Stereographic: ',EMEP_projection_attributes(1:6)
        endif

        !Find the (x,y,z,time,xdist,ydist) dimensions of the file
        do i_dim=1,num_dims_nc
            status_nc = NF90_INQ_DIMID (id_nc,dim_name_nc(i_dim),dim_id_nc(i_dim))
            status_nc = NF90_INQUIRE_DIMENSION (id_nc,dim_id_nc(i_dim),dimname_temp,dim_length_nc(i_dim))
            if (status_nc .NE. NF90_NOERR) then
                write(unit_logfile,'(A,A,A,I)') 'No dimension information available for ',trim(dim_name_nc(i_dim)),' Setting to 1 with status: ',status_nc
                dim_length_nc(i_dim)=1
            endif
        enddo
            
        if (subgrid_dim(t_dim_index).gt.dim_length_nc(time_dim_nc_index)) then
            write(unit_logfile,'(A,2I)') 'ERROR: Specified time dimensions are greater than EMEP netcdf dimensions. Stopping ',subgrid_dim(t_dim_index),dim_length_nc(time_dim_nc_index)
            stop
        endif 
              
        write(unit_logfile,'(A,6I)') ' Size of dimensions (x,y,z,t,xdist,ydist): ',dim_length_nc
        dim_length_EMEP_nc=dim_length_nc
        
        dim_start_nc(time_dim_nc_index)=temp_start_time_nc_index
        dim_length_nc(time_dim_nc_index)=min(dim_length_nc(time_dim_nc_index),subgrid_dim(t_dim_index))
       
        write(unit_logfile,'(A,6I)') ' New size of dimensions (x,y,z,t,xdist,ydist): ',dim_length_nc

        if (mod(dim_length_nc(xdist_dim_nc_index),2).ne.1.or.mod(dim_length_nc(ydist_dim_nc_index),2).ne.1) then
            write(unit_logfile,'(A,2I)') ' ERROR: Even sized dimensions for local contribution. Must be odd: ',dim_length_nc(xdist_dim_nc_index),dim_length_nc(ydist_dim_nc_index)
            stop
        endif
            
        if (i_file.eq.2) then
            xdist_centre_nc=1+dim_length_nc(xdist_dim_nc_index)/2
            ydist_centre_nc=1+dim_length_nc(ydist_dim_nc_index)/2
            write(unit_logfile,'(A,2I)') ' Centre index of local contribution dimensions: ',xdist_centre_nc,ydist_centre_nc
        endif
        
        !Calculate the necessary extent of the EMEP grid region and only read these grids
        if (reduce_EMEP_region_flag) then
            !EMEP_grid_interpolation_size_temp=max(EMEP_grid_interpolation_size*local_fraction_grid_size_scaling,EMEP_additional_grid_interpolation_size_original*local_fraction_additional_grid_size_scaling)
            EMEP_grid_interpolation_size_temp=EMEP_grid_interpolation_size*local_fraction_grid_size_scaling
            
            write(unit_logfile,'(A,f12.2)') 'Reducing EMEP domain. EMEP grid interpolation size is now = ',EMEP_grid_interpolation_size_temp
            
            !Determine the LL cordinates of the target grid
            !if (EMEP_projection_type.eq.LCC_projection_index) then
                !Retrieve the four corners of the target grid in lat and lon
            call PROJ2LL(init_subgrid_min(x_dim_index),init_subgrid_min(y_dim_index),temp_lon(1),temp_lat(1),projection_attributes,projection_type)
            call PROJ2LL(init_subgrid_max(x_dim_index),init_subgrid_max(y_dim_index),temp_lon(2),temp_lat(2),projection_attributes,projection_type)
            call PROJ2LL(init_subgrid_min(x_dim_index),init_subgrid_max(y_dim_index),temp_lon(3),temp_lat(3),projection_attributes,projection_type)
            call PROJ2LL(init_subgrid_max(x_dim_index),init_subgrid_min(y_dim_index),temp_lon(4),temp_lat(4),projection_attributes,projection_type)
            
            
            !if (projection_type.eq.RDM_projection_index) then
            !    call RDM2LL(init_subgrid_min(y_dim_index),init_subgrid_min(x_dim_index),temp_lat(1),temp_lon(1))
            !    call RDM2LL(init_subgrid_max(y_dim_index),init_subgrid_max(x_dim_index),temp_lat(2),temp_lon(2))
            !    call RDM2LL(init_subgrid_max(y_dim_index),init_subgrid_min(x_dim_index),temp_lat(3),temp_lon(3))
            !    call RDM2LL(init_subgrid_min(y_dim_index),init_subgrid_max(x_dim_index),temp_lat(4),temp_lon(4))
            !elseif (projection_type.eq.UTM_projection_index) then
            !    call UTM2LL(utm_zone,init_subgrid_min(y_dim_index),init_subgrid_min(x_dim_index),temp_lat(1),temp_lon(1))
            !    call UTM2LL(utm_zone,init_subgrid_max(y_dim_index),init_subgrid_max(x_dim_index),temp_lat(2),temp_lon(2))
            !    call UTM2LL(utm_zone,init_subgrid_max(y_dim_index),init_subgrid_min(x_dim_index),temp_lat(3),temp_lon(3))
            !    call UTM2LL(utm_zone,init_subgrid_min(y_dim_index),init_subgrid_max(x_dim_index),temp_lat(4),temp_lon(4))
            !endif   
                
                !This did not work because it was almost all of the grid and because the min and max lat lon did not cover all stations
                !if (read_EMEP_only_once_flag.and.use_multiple_receptor_grids_flag) then
                !    temp_lat(1)=minval(lat_receptor(1:n_receptor_in));temp_lon(1)=minval(lon_receptor(1:n_receptor_in))
                !    temp_lat(2)=maxval(lat_receptor(1:n_receptor_in));temp_lon(2)=maxval(lon_receptor(1:n_receptor_in))
                !    temp_lat(3)=maxval(lat_receptor(1:n_receptor_in));temp_lon(3)=minval(lon_receptor(1:n_receptor_in))
                !    temp_lat(4)=minval(lat_receptor(1:n_receptor_in));temp_lon(4)=maxval(lon_receptor(1:n_receptor_in))
                !    write(*,*) temp_lat
                !    write(*,*) temp_lon
                !endif
                
                temp_x_min=1.e32;temp_y_min=1.e32
                temp_x_max=-1.e32;temp_y_max=-1.e32

                    if (EMEP_projection_type.eq.LCC_projection_index) then
                        !Convert lat lon corners to lambert
                        do i=1,4
                            !if (use_alternative_LCC_projection_flag) then
                                call lb2lambert2_uEMEP(temp_x(i),temp_y(i),temp_lon(i),temp_lat(i),EMEP_projection_attributes)
                            !else
                            !    call lb2lambert_uEMEP(temp_x(i),temp_y(i),temp_lon(i),temp_lat(i),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                            !endif
                            !call lb2lambert_uEMEP(temp_x(i),temp_y(i),temp_lon(i),temp_lat(i),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                        enddo            
                            !write(*,*) temp_x
                            !write(*,*) temp_y
                    elseif (EMEP_projection_type.eq.PS_projection_index) then
                        !Convert lat lon corners to Polar Stereo
                        do i=1,4
                            call LL2PS_spherical(temp_x(i),temp_y(i),temp_lon(i),temp_lat(i),EMEP_projection_attributes)
                        enddo            
                    elseif (EMEP_projection_type.eq.LL_projection_index) then
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
                    !write(*,*) temp_x(i),temp_y(i)
                    if (temp_x(i).lt.temp_x_min) temp_x_min=temp_x(i)
                    if (temp_y(i).lt.temp_y_min) temp_y_min=temp_y(i)
                    if (temp_x(i).gt.temp_x_max) temp_x_max=temp_x(i)
                    if (temp_y(i).gt.temp_y_max) temp_y_max=temp_y(i)
                enddo
            
                !Read in the first 2 x and y position values from the nc file to get min values and delta values
                !write(*,*) temp_x_min,temp_x_max,temp_y_min,temp_y_max
                
                !Save these values, min and max extent of the calculation grid in emep coordinates, for possible use later
                if (limit_emep_grid_interpolation_region_to_calculation_region) then
                subgrid_proj_min(y_dim_index)=temp_y_min
                subgrid_proj_max(y_dim_index)=temp_y_max
                subgrid_proj_min(x_dim_index)=temp_x_min
                subgrid_proj_max(x_dim_index)=temp_x_max
                endif
                
                status_nc = NF90_INQ_VARID (id_nc, trim(dim_name_nc(x_dim_nc_index)), var_id_nc)
                status_nc = NF90_GET_VAR (id_nc, var_id_nc,temp_var1d_nc_dp(1,1:2),start=(/1/),count=(/2/))
                status_nc = NF90_INQ_VARID (id_nc, trim(dim_name_nc(y_dim_nc_index)), var_id_nc)
                status_nc = NF90_GET_VAR (id_nc, var_id_nc,temp_var1d_nc_dp(2,1:2),start=(/1/),count=(/2/))
                status_nc = nf90_get_att(id_nc, var_id_nc, "units", unit_name_nc_temp)
                if (trim(unit_name_nc_temp).eq.'km') then
                    write(unit_logfile,'(A)') 'Units of x y data are in kilometres. Converting to metres'
                    temp_var1d_nc_dp=temp_var1d_nc_dp*1000.
                endif
                
                !write(*,*) temp_var1d_nc_dp
                temp_delta(1)=temp_var1d_nc_dp(1,2)-temp_var1d_nc_dp(1,1)
                temp_delta(2)=temp_var1d_nc_dp(2,2)-temp_var1d_nc_dp(2,1)
                !write(*,*) temp_delta
                !Find grid position of the max and min coordinates and add2 grids*EMEP_grid_interpolation_size_temp
                i_temp_min=1+floor((temp_x_min-temp_var1d_nc_dp(1,1))/temp_delta(1)+0.5)
                i_temp_max=1+floor((temp_x_max-temp_var1d_nc_dp(1,1))/temp_delta(1)+0.5)
                j_temp_min=1+floor((temp_y_min-temp_var1d_nc_dp(2,1))/temp_delta(2)+0.5)
                j_temp_max=1+floor((temp_y_max-temp_var1d_nc_dp(2,1))/temp_delta(2)+0.5)
                i_temp_min=1+floor((temp_x_min-temp_var1d_nc_dp(1,1))/temp_delta(1)+0.5)
                i_temp_max=1+ceiling((temp_x_max-temp_var1d_nc_dp(1,1))/temp_delta(1)+0.5)
                j_temp_min=1+floor((temp_y_min-temp_var1d_nc_dp(2,1))/temp_delta(2)+0.5)
                j_temp_max=1+ceiling((temp_y_max-temp_var1d_nc_dp(2,1))/temp_delta(2)+0.5)
                !write(unit_logfile,'(A,2I)') ' Reading EMEP i grids: ',i_temp_min,i_temp_max
                !write(unit_logfile,'(A,2I)') ' Reading EMEP j grids: ',j_temp_min,j_temp_max
                i_temp_min=max(1,i_temp_min-1-ceiling(1.*EMEP_grid_interpolation_size_temp))
                i_temp_max=min(dim_length_nc(x_dim_nc_index),i_temp_max+1+ceiling(1.*EMEP_grid_interpolation_size_temp))
                j_temp_min=max(1,j_temp_min-1-ceiling(1.*EMEP_grid_interpolation_size_temp))
                j_temp_max=min(dim_length_nc(y_dim_nc_index),j_temp_max+1+ceiling(1.*EMEP_grid_interpolation_size_temp))
                dim_length_nc(x_dim_nc_index)=i_temp_max-i_temp_min+1
                dim_length_nc(y_dim_nc_index)=j_temp_max-j_temp_min+1
                dim_start_nc(x_dim_nc_index)=i_temp_min
                dim_start_nc(y_dim_nc_index)=j_temp_min
                write(unit_logfile,'(A,3I)') ' Reading EMEP i grids: ',i_temp_min,i_temp_max,dim_length_nc(x_dim_nc_index)
                write(unit_logfile,'(A,3I)') ' Reading EMEP j grids: ',j_temp_min,j_temp_max,dim_length_nc(y_dim_nc_index)
                dim_start_EMEP_nc=dim_start_nc
            !endif

        endif
        
        if (dim_length_nc(x_dim_nc_index).lt.1.or.dim_length_nc(y_dim_nc_index).lt.1) then
            write(unit_logfile,'(A,2I)') ' WARNING: Selected EMEP region dimensions are less than 1 (i,j): ',dim_length_nc
            write(unit_logfile,'(A)') ' Setting to 1 but this selected region is invalid and should not be calculated '
            write(unit_logfile,'(A)') ' This can happen if a receptor is outside the EMEP domain'
            dim_length_nc=1
            dim_start_nc(x_dim_nc_index)=1
            dim_start_nc(y_dim_nc_index)=1
            EMEP_region_outside_domain=.true.
        endif

        !Allocate the nc arrays for reading
        if (.not.allocated(val_dim_nc)) allocate (val_dim_nc(maxval(dim_length_nc),num_dims_nc)) !x, y, z and time dimension values
        if (.not.allocated(unit_dim_nc)) allocate (unit_dim_nc(num_dims_nc)) !x, y, z and time dimension values
        if (.not.allocated(var1d_nc)) allocate (var1d_nc(maxval(dim_length_nc),num_dims_nc)) !x, y, z and time maximum dimensions
        if (.not.allocated(var2d_nc)) allocate (var2d_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),2)) !Lat and lon
        if (.not.allocated(var3d_nc)) then
            allocate (var3d_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index),num_var_nc,n_source_nc_index,n_pollutant_loop))
            var3d_nc=0.
        endif
        if (.not.allocated(var4d_nc)) then
            allocate (var4d_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index),num_var_nc,n_source_nc_index,n_pollutant_loop))
            var4d_nc=0.
        endif
        if (.not.allocated(comp_var3d_nc)) then
            allocate (comp_var3d_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index),n_compound_nc_index))
            comp_var3d_nc=0.
        endif
        if (.not.allocated(comp_var4d_nc)) then
            allocate (comp_var4d_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index),n_compound_nc_index))
            comp_var4d_nc=0.
        endif
        if (.not.allocated(var1d_nc_dp)) allocate (var1d_nc_dp(maxval(dim_length_nc))) 
        if (.not.allocated(var2d_nc_dp)) allocate (var2d_nc_dp(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index))) !Lat and lon
        if (.not.allocated(temp_var3d_nc)) allocate (temp_var3d_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index))) !Lat and lon
        !if (.not.allocated(var3d_nc_dp)) allocate (var3d_nc_dp(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)))
        !allocate (var4d_nc_dp(dim_length_nc(x_index),dim_length_nc(y_index),1,dim_length_nc(time_index)))
        
        
        if (calculate_deposition_flag) then
        if (.not.allocated(depo_var3d_nc)) then
            allocate (depo_var3d_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index),n_landuse_index,n_pollutant_loop))
            depo_var3d_nc=0
        endif
        endif

        if (save_emep_species.or.save_seasalt) then
        if (.not.allocated(species_var3d_nc)) then
            allocate (species_var3d_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index),n_pmxx_sp_index,n_species_loop_index))
            species_var3d_nc=0
        endif
        if (.not.allocated(species_temp_var3d_nc)) allocate (species_temp_var3d_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)))        
        endif
        
        if (i_file.eq.2.and..not.allocated(lc_var3d_nc)) then
            allocate (lc_var3d_nc(dim_length_nc(xdist_dim_nc_index),dim_length_nc(ydist_dim_nc_index),dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index),num_lc_var_nc,n_source_nc_index,n_pollutant_loop))
            lc_var3d_nc=0.
        endif
        if (i_file.eq.2.and..not.allocated(lc_var4d_nc)) then
            allocate (lc_var4d_nc(dim_length_nc(xdist_dim_nc_index),dim_length_nc(ydist_dim_nc_index),dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),1,dim_length_nc(time_dim_nc_index),num_lc_var_nc,n_source_nc_index,n_pollutant_loop))
            lc_var4d_nc=0.
        endif
        if (i_file.eq.2.and..not.allocated(pm_lc_var4d_nc)) then
            allocate (pm_lc_var4d_nc(dim_length_nc(xdist_dim_nc_index),dim_length_nc(ydist_dim_nc_index),dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),1,dim_length_nc(time_dim_nc_index),num_lc_var_nc,n_source_nc_index,2))
            pm_lc_var4d_nc=0.
        endif
        if (.not.allocated(pm_var4d_nc)) then
            allocate (pm_var4d_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index),num_var_nc,n_source_nc_index,2))
            pm_var4d_nc=0.
        endif
        if (.not.allocated(pm_var3d_nc)) then
            allocate (pm_var3d_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index),num_var_nc,n_source_nc_index,2))
            pm_var3d_nc=0.
        endif
        if (.not.allocated(temp_var4d_nc).and.i_file.eq.1) then
            allocate (temp_var4d_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index),2))
            temp_var4d_nc=0.
        endif
        
        !write(*,*) x_dim_nc_index,y_dim_nc_index
        !write(*,*) shape(var1d_nc_dp)
        !write(*,*) dim_length_nc
        !Read in the dimensions and check values of the dimensions. Not necessary but diagnostic
        do i=1,num_dims_nc
            status_nc = NF90_INQ_VARID (id_nc, trim(dim_name_nc(i)), var_id_nc)
            !write(*,*) id_nc, trim(dim_name_nc(i)), var_id_nc(i),dim_length_nc(i)
            var1d_nc_dp=0.
            !write(*,*) 'HERE',i,dim_start_nc(i),dim_length_nc(i)
            unit_dim_nc(i)=''
            if (status_nc .EQ. NF90_NOERR) then
                status_nc = nf90_get_att(id_nc, var_id_nc, "units", unit_dim_nc(i))
                status_nc = NF90_GET_VAR (id_nc, var_id_nc,var1d_nc_dp(1:dim_length_nc(i)),start=(/dim_start_nc(i)/),count=(/dim_length_nc(i)/));var1d_nc(1:dim_length_nc(i),i)=real(var1d_nc_dp(1:dim_length_nc(i)))  
                !write(*,*) id_nc, trim(dim_name_nc(i)), var_id_nc,dim_length_nc(i),status_nc
                !Use the first file to give valid time stamps
                if (i_file.eq.1.and.i.eq.time_dim_nc_index) then
                    val_dim_nc(1:dim_length_nc(i),i)=(var1d_nc_dp(1:dim_length_nc(i)))
                    valid_dim_length_nc(i)=dim_length_nc(i)
                elseif (i_file.ne.1.and.i.ne.time_dim_nc_index) then 
                    val_dim_nc(1:dim_length_nc(i),i)=(var1d_nc_dp(1:dim_length_nc(i)))
                    valid_dim_length_nc(i)=dim_length_nc(i)
                endif
            !write(*,*) val_dim_nc(1:dim_length_nc(i),i),trim(unit_dim_nc(i))
            else
                var1d_nc(1:dim_length_nc(i),i)=0.
                val_dim_nc(1:dim_length_nc(i),i)=0.
            endif
            
            !Convert from meters to km for AROME data if necessary
            if ((i.eq.x_dim_nc_index.or.i.eq.y_dim_nc_index).and.trim(unit_dim_nc(i)).eq.'km') then
                write(unit_logfile,'(A)') 'Units of x y data are in kilometres. Converting to metres'
                val_dim_nc(1:dim_length_nc(i),i)=val_dim_nc(1:dim_length_nc(i),i)*1000.
                var1d_nc(1:dim_length_nc(i),i)=var1d_nc(1:dim_length_nc(i),i)*1000.
            endif
            
            if (i.eq.time_dim_nc_index) then
                !write(unit_logfile,'(3A,2i12)') ' ',trim(dim_name_nc(i)),' (min, max in hours): ' &
                !    ,minval(int((var1d_nc(1:dim_length_nc(i),i)-var1d_nc(dim_start_nc(i),i))/3600.+.5)+1) &
                !    ,maxval(int((var1d_nc(1:divar1d_ncm_length_nc(i),i)-var1d_nc(dim_start_nc(i),i))/3600.+.5)+1)                     
            else
                write(unit_logfile,'(3A,2f12.2)') ' ',trim(dim_name_nc(i)),' (min, max): ' &
                    ,minval(var1d_nc(1:dim_length_nc(i),i)),maxval(var1d_nc(1:dim_length_nc(i),i)) 
            endif       
        enddo

     
        
        !Loop through the pollutants.
        do p_loop=1,n_emep_pollutant_loop+1

        !Set the compound index for this pollutant
        if (p_loop.le.n_emep_pollutant_loop) then
            i_pollutant=pollutant_loop_index(p_loop)
            p_loop_index=p_loop
        else
            !Special case to read in meteo data that is defined as the all_nc_index
            i_pollutant=all_nc_index
            !Put the meteo data in the first pollutant index always
            p_loop_index=meteo_p_loop_index
        endif
        
        !Loop through the sources
        do i_source=1,n_source_nc_index
        !write(*,*) i_source,trim(source_file_str(i_source)),uEMEP_to_EMEP_sector(i_source),calculate_source(i_source),calculate_EMEP_source(i_source)
        if (calculate_source(i_source).or.calculate_EMEP_source(i_source).or.save_EMEP_source(i_source).or.i_source.eq.allsource_index.or.(i_source.eq.extrasource_nc_index.and.use_alternative_ppm_variable_for_lf)) then
        !var_name_nc(num_var_nc,n_compound_nc_index,n_source_nc_index)
            
        !Loop through the variables
        do i=1,num_var_nc
            !write(*,*) i,trim(var_name_nc(i))
            !if (i.eq.frac_nc_index) var_name_nc_temp=var_name_nc(i,i_pollutant,i_source)
            !if (i.eq.conc_nc_index) var_name_nc_temp=var_name_nc(i,i_pollutant,i_source)
            
            !Identify the variable name and ID in the nc file
            var_name_nc_temp=var_name_nc(i,i_pollutant,i_source)
            status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
            !write(*,*) 'Status1: ',status_nc,var_id_nc,trim(var_name_nc_temp),i_source
            !Exception for pm10
            !if ((i.eq.conc_nc_index.or.i.eq.frac_nc_index.or.i.eq.emis_nc_index).and.i_pollutant.eq.pm10_nc_index) then
            if ((i.eq.conc_nc_index.or.(i.ge.min_frac_nc_loop_index.and.i.le.max_frac_nc_loop_index).or.i.eq.emis_nc_index).and.i_pollutant.eq.pm10_nc_index) then
                var_name_nc_temp=var_name_nc(i,pmco_nc_index,i_source)
                status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                !write(*,*) '-------------------',trim(var_name_nc_temp),status_nc
            endif
            
            !If a variable name is found in the file then go further
            if (status_nc.eq.NF90_NOERR) then
                scale_factor_nc=1.
                !Find the dimensions of the variable (temp_num_dims)
                status_nc = NF90_INQUIRE_VARIABLE(id_nc, var_id_nc, ndims = temp_num_dims)
                !write(*,*) temp_num_dims,status_nc
                if (temp_num_dims.eq.2.and.i_file.eq.1) then
                    !Read latitude and longitude data into a 2d grid if available. Only lat lon is 2d? Only read for file 1 assuming file 2 is the same. Problem if not
                    if (i.eq.lat_nc_index.or.i.eq.lon_nc_index) then
                    status_nc = NF90_GET_VAR (id_nc, var_id_nc, var2d_nc_dp);var2d_nc(:,:,i)=real(var2d_nc_dp)
                    write(unit_logfile,'(A,i3,A,2A,2f16.4)') ' Reading: ',temp_num_dims,' ',trim(var_name_nc_temp),' (min, max): ',minval(var2d_nc(:,:,i)),maxval(var2d_nc(:,:,i))
                    endif
                elseif (temp_num_dims.eq.3.and.i_file.eq.1.and.i.eq.emis_nc_index.and.i_pollutant.eq.pm10_nc_index.and.i_source.ne.extrasource_nc_index) then
                        var_name_nc_temp2=var_name_nc(i,pmco_nc_index,i_source)
                        status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp2), var_id_nc)
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, pm_var3d_nc(:,:,:,i,i_source,1),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),dim_start_nc(time_dim_nc_index)/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,I,3A,2f16.4)') ' Reading: ',temp_num_dims,' ',trim(var_name_nc_temp2),' (min, max): ',minval(pm_var3d_nc(:,:,:,i,i_source,1)),maxval(pm_var3d_nc(:,:,:,i,i_source,1))
                        write(unit_logfile,'(2A,f16.4)') ' Average of: ',trim(var_name_nc_temp2),sum(pm_var3d_nc(:,:,:,i,i_source,1))/(size(pm_var3d_nc,1)*size(pm_var3d_nc,2)*size(pm_var3d_nc,4))
                        var_name_nc_temp2=var_name_nc(i,pm25_nc_index,i_source)
                        status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp2), var_id_nc)
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, pm_var3d_nc(:,:,:,i,i_source,2),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),dim_start_nc(time_dim_nc_index)/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,I,3A,2f16.4)') ' Reading: ',temp_num_dims,' ',trim(var_name_nc_temp2),' (min, max): ',minval(pm_var3d_nc(:,:,:,i,i_source,2)),maxval(pm_var3d_nc(:,:,:,i,i_source,2))
                        write(unit_logfile,'(2A,f16.4,3i)') ' Average of: ',trim(var_name_nc_temp2),sum(pm_var3d_nc(:,:,:,i,i_source,2))/(size(pm_var3d_nc,1)*size(pm_var3d_nc,2)*size(pm_var3d_nc,4)),size(pm_var3d_nc,1),size(pm_var3d_nc,2),size(pm_var3d_nc,4)
                        var3d_nc(:,:,:,i,i_source,p_loop_index)=pm_var3d_nc(:,:,:,i,i_source,1)+pm_var3d_nc(:,:,:,i,i_source,2)
                        write(unit_logfile,'(2A,f16.4,3i)') ' Average of: ',trim(var_name_nc(i,pm10_nc_index,i_source)),sum(var3d_nc(:,:,:,i,i_source,p_loop_index))/(size(var3d_nc,1)*size(var3d_nc,2)*size(var3d_nc,4)),size(var3d_nc,1),size(var3d_nc,2),size(var3d_nc,4)
                elseif (temp_num_dims.eq.3.and.i_file.eq.1) then
                    !write(*,'(6i)') dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),temp_start_time_nc_index,dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)
                    status_nc = NF90_GET_VAR (id_nc, var_id_nc, var3d_nc(:,:,:,i,i_source,p_loop_index),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),dim_start_nc(time_dim_nc_index)/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                    !write(*,*) status_nc
                    write(unit_logfile,'(A,I,A,I,3A,2f16.4)') ' Reading: ',temp_num_dims,' p_loop:',p_loop_index,' ',trim(var_name_nc_temp),' (min, max): ',minval(var3d_nc(:,:,:,i,i_source,p_loop_index)),maxval(var3d_nc(:,:,:,i,i_source,p_loop_index))
                    !if (i.eq.emis_nc_index) write(*,*) 'HERE',sum(var3d_nc(:,:,:,i,i_source,p_loop_index)),i_source,p_loop_index
                elseif (temp_num_dims.eq.4) then
                    if (i_file.eq.2.and.i.eq.ZTOP_nc_index) then
                        !Don't try to read
                    elseif (i.eq.conc_nc_index.and.i_pollutant.eq.pm10_nc_index.and.i_source.ne.extrasource_nc_index) then
                        var_name_nc_temp2=var_name_nc(i,pmco_nc_index,i_source)
                        status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp2), var_id_nc)
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, pm_var4d_nc(:,:,dim_start_nc(z_dim_nc_index):dim_start_nc(z_dim_nc_index)+dim_length_nc(z_dim_nc_index)-1,:,i,i_source,1),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),dim_start_nc(z_dim_nc_index),dim_start_nc(time_dim_nc_index)/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,I,A,I,3A,2f16.4)') ' Reading: ',temp_num_dims,' p_loop:',p_loop_index,' ',trim(var_name_nc_temp2),' (min, max): ',minval(pm_var4d_nc(:,:,:,:,i,i_source,1)),maxval(pm_var4d_nc(:,:,:,:,i,i_source,1))
                        write(unit_logfile,'(2A,f16.4)') ' Average of: ',trim(var_name_nc_temp2),sum(pm_var4d_nc(:,:,1,:,i,i_source,1))/(size(pm_var4d_nc,1)*size(pm_var4d_nc,2)*size(pm_var4d_nc,4))
                        var_name_nc_temp2=var_name_nc(i,pm25_nc_index,i_source)
                        status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp2), var_id_nc)
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, pm_var4d_nc(:,:,dim_start_nc(z_dim_nc_index):dim_start_nc(z_dim_nc_index)+dim_length_nc(z_dim_nc_index)-1,:,i,i_source,2),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),dim_start_nc(z_dim_nc_index),dim_start_nc(time_dim_nc_index)/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,I,A,I,3A,2f16.4)') ' Reading: ',temp_num_dims,' p_loop:',p_loop_index,' ',trim(var_name_nc_temp2),' (min, max): ',minval(pm_var4d_nc(:,:,:,:,i,i_source,2)),maxval(pm_var4d_nc(:,:,:,:,i,i_source,2))
                        write(unit_logfile,'(2A,f16.4,3i)') ' Average of: ',trim(var_name_nc_temp2),sum(pm_var4d_nc(:,:,1,:,i,i_source,2))/(size(pm_var4d_nc,1)*size(pm_var4d_nc,2)*size(pm_var4d_nc,4)),size(pm_var4d_nc,1),size(pm_var4d_nc,2),size(pm_var4d_nc,4)
                        var4d_nc(:,:,:,:,i,i_source,p_loop_index)=pm_var4d_nc(:,:,:,:,i,i_source,1)+pm_var4d_nc(:,:,:,:,i,i_source,2)
                        write(unit_logfile,'(2A,f16.4,3i)') ' Average of: ',trim(var_name_nc(i,pm10_nc_index,i_source)),sum(var4d_nc(:,:,:,:,i,i_source,p_loop_index))/(size(var4d_nc,1)*size(var4d_nc,2)*size(var4d_nc,4)),size(var4d_nc,1),size(var4d_nc,2),size(var4d_nc,4)
                    else
                    status_nc = NF90_GET_VAR (id_nc, var_id_nc, var4d_nc(:,:,dim_start_nc(z_dim_nc_index):dim_start_nc(z_dim_nc_index)+dim_length_nc(z_dim_nc_index)-1,:,i,i_source,p_loop_index),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),dim_start_nc(z_dim_nc_index),dim_start_nc(time_dim_nc_index)/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                    !status_nc = NF90_GET_VAR (id_nc, var_id_nc, temp_var4d_nc(:,:,:,:),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),dim_start_nc(z_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                    !var4d_nc(:,:,:,:,i,i_source)=real(temp_var4d_nc(:,:,:,:))
                    write(unit_logfile,'(A,I,A,I,3A,2f16.4)') ' Reading: ',temp_num_dims,' p_loop:',p_loop_index,' ',trim(var_name_nc_temp),' (min, max): ',minval(var4d_nc(:,:,:,:,i,i_source,p_loop_index)),maxval(var4d_nc(:,:,:,:,i,i_source,p_loop_index))
                    !if (i.eq.precip_nc_index) then 
                    !    write(*,*) maxval(var4d_nc(:,:,:,:,i,i_source,p_loop_index))
                    !    stop
                    !endif
                    
                    endif
                    !write(*,*) shape(var4d_nc)
                    !write(*,*) dim_start_nc(z_dim_nc_index),dim_length_nc(z_dim_nc_index)
                    !write(*,*) maxval(var4d_nc(:,:,1,1,i,i_source)),maxval(var4d_nc(:,:,1,2,i,i_source))
                elseif (temp_num_dims.eq.6.and.i_file.eq.2) then
                    !if (i.eq.frac_nc_index.and.i_pollutant.eq.pm10_nc_index) then
                    lc_frac_nc_index=convert_frac_to_lc_frac_loop_index(i)
                    !write(*,*) i,lc_frac_nc_index
                    if (i.ge.min_frac_nc_loop_index.and.i.le.max_frac_nc_loop_index.and.i_pollutant.eq.pm10_nc_index) then
                        var_name_nc_temp2=var_name_nc(i,pmco_nc_index,i_source)
                        status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp2), var_id_nc)
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, pm_lc_var4d_nc(:,:,:,:,:,:,lc_frac_nc_index,i_source,1),start=(/1,1,dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),dim_start_nc(z_dim_nc_index),dim_start_nc(time_dim_nc_index)/),count=(/dim_length_nc(xdist_dim_nc_index),dim_length_nc(ydist_dim_nc_index),dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,2I,3A,2f16.4)') ' Reading: ',i,temp_num_dims,' ',trim(var_name_nc_temp2),' (min, max): ',minval(pm_lc_var4d_nc(:,:,:,:,:,:,lc_frac_nc_index,i_source,1)),maxval(pm_lc_var4d_nc(:,:,:,:,:,:,lc_frac_nc_index,i_source,1))
                        var_name_nc_temp2=var_name_nc(i,pm25_nc_index,i_source)
                        status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp2), var_id_nc)
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, pm_lc_var4d_nc(:,:,:,:,:,:,lc_frac_nc_index,i_source,2),start=(/1,1,dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),dim_start_nc(z_dim_nc_index),dim_start_nc(time_dim_nc_index)/),count=(/dim_length_nc(xdist_dim_nc_index),dim_length_nc(ydist_dim_nc_index),dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,2I,3A,2f16.4)') ' Reading: ',i,temp_num_dims,' ',trim(var_name_nc_temp2),' (min, max): ',minval(pm_lc_var4d_nc(:,:,:,:,:,:,lc_frac_nc_index,i_source,2)),maxval(pm_lc_var4d_nc(:,:,:,:,:,:,lc_frac_nc_index,i_source,2))
                        !Not used but calculated for writing
                        !lc_var4d_nc(:,:,:,:,:,:,lc_frac_nc_index,i_source,p_loop_index)=pm_lc_var4d_nc(:,:,:,:,:,:,lc_frac_nc_index,i_source,1)+pm_lc_var4d_nc(:,:,:,:,:,:,lc_frac_nc_index,i_source,2)
                    else
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, lc_var4d_nc(:,:,:,:,:,:,lc_frac_nc_index,i_source,p_loop_index),start=(/1,1,dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),dim_start_nc(z_dim_nc_index),dim_start_nc(time_dim_nc_index)/),count=(/dim_length_nc(xdist_dim_nc_index),dim_length_nc(ydist_dim_nc_index),dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,2I,3A,2f16.4)') ' Reading: ',i,temp_num_dims,' ',trim(var_name_nc_temp),' (min, max): ',minval(lc_var4d_nc(:,:,:,:,:,:,lc_frac_nc_index,i_source,p_loop_index)),maxval(lc_var4d_nc(:,:,:,:,:,:,lc_frac_nc_index,i_source,p_loop_index))
                    endif
                    !write(*,*) shape(lc_var4d_nc)
                    !write(*,*) maxval(lc_var4d_nc(3,3,:,:,:,1,lc_frac_nc_index,i_source)),maxval(lc_var4d_nc(3,3,:,:,:,2,lc_frac_nc_index,i_source))

                else
                    write(unit_logfile,'(8A,8A)') ' Cannot find a correct dimension for: ',trim(var_name_nc_temp)
                endif    
                
            else
                 !write(unit_logfile,'(8A,8A)') ' Cannot read: ',trim(var_name_nc_temp)
            endif

        enddo !Variable loop
        
        
        endif
        enddo !Source loop
        
        !Loop through the additional compounds that are in the base file
        !if (i_file.eq.1) then
        i_source=allsource_index
        !i=conc_nc_index
        do i_loop=1,n_pollutant_compound_loop(p_loop)
            i_conc=pollutant_compound_loop_index(p_loop,i_loop)
            var_name_nc_temp=comp_name_nc(i_conc)
            !write(*,*) p_loop,i_conc,i_source,trim(var_name_nc_temp)
            status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
            !write(*,*) 'Status1: ',status_nc,id_nc,var_id_nc,trim(var_name_nc_temp)
            
            !If a variable name is found in the file then go further
            if (status_nc.eq.NF90_NOERR) then

                !Find the dimensions of the variable (temp_num_dims)
                status_nc = NF90_INQUIRE_VARIABLE(id_nc, var_id_nc, ndims = temp_num_dims)

                if (temp_num_dims.eq.3) then
                    if (i_file.eq.1) then
                    status_nc = NF90_GET_VAR (id_nc, var_id_nc, temp_var3d_nc(:,:,:),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                    comp_var4d_nc(:,:,surface_level_nc,:,i_conc)=temp_var3d_nc(:,:,:)*comp_scale_nc(i_conc)
                    write(unit_logfile,'(A,I,3A,2f16.4)') ' Reading compound file 1: ',temp_num_dims,' ',trim(var_name_nc_temp),' (min, max): ',minval(comp_var4d_nc(:,:,surface_level_nc,:,i_conc)),maxval(comp_var4d_nc(:,:,surface_level_nc,:,i_conc))
                    !write(*,*) comp_var4d_nc(:,:,1,:,i_conc)
                    elseif (i_file.eq.2) then
                    !In case the comp data is in the uEMEP file then read it here with no vertical extent
                    status_nc = NF90_GET_VAR (id_nc, var_id_nc, temp_var3d_nc(:,:,:),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                    comp_var4d_nc(:,:,surface_level_nc,:,i_conc)=temp_var3d_nc(:,:,:)*comp_scale_nc(i_conc)
                    write(unit_logfile,'(A,I,3A,2f16.4)') ' Reading compound file 2: ',temp_num_dims,' ',trim(var_name_nc_temp),' (min, max): ',minval(comp_var4d_nc(:,:,surface_level_nc,:,i_conc)),maxval(comp_var4d_nc(:,:,surface_level_nc,:,i_conc))
                    endif
                endif
                if (temp_num_dims.eq.4) then
                    if (i_file.eq.1) then
                        !write(*,'(4i)') dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),dim_start_nc(z_dim_nc_index),temp_start_time_nc_index
                        !write(*,'(4i)') dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index)
                        
                    status_nc = NF90_GET_VAR (id_nc, var_id_nc, comp_var4d_nc(:,:,:,:,i_conc),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),dim_start_nc(z_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                    comp_var4d_nc(:,:,:,:,i_conc)=comp_var4d_nc(:,:,:,:,i_conc)*comp_scale_nc(i_conc)
                    write(unit_logfile,'(A,I,3A,2f16.4)') ' Reading compound file 1: ',temp_num_dims,' ',trim(var_name_nc_temp),' (min, max): ',minval(comp_var4d_nc(:,:,:,:,i_conc)),maxval(comp_var4d_nc(:,:,:,:,i_conc))
                   
                    elseif (i_file.eq.2) then
                    !In case the comp data is in the uEMEP file then read it here with no vertical extent
                    status_nc = NF90_GET_VAR (id_nc, var_id_nc, comp_var4d_nc(:,:,1,:,i_conc),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),1,temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),1,dim_length_nc(time_dim_nc_index)/))
                    comp_var4d_nc(:,:,1,:,i_conc)=comp_var4d_nc(:,:,1,:,i_conc)*comp_scale_nc(i_conc)
                    write(unit_logfile,'(A,I,3A,2f16.4)') ' Reading compound file 2: ',temp_num_dims,' ',trim(var_name_nc_temp),' (min, max): ',minval(comp_var4d_nc(:,:,1,:,i_conc)),maxval(comp_var4d_nc(:,:,1,:,i_conc))
                    endif
                endif              
                
            else
                 write(unit_logfile,'(8A,8A)') ' Cannot read compound: ',trim(var_name_nc_temp)
            endif    
                
        enddo !pollutant
        enddo !compound loop
        !endif

        if (calculate_deposition_flag.and.i_file.eq.1) then
            write(unit_logfile,'(A)') ' Reading deposition velocity data from EMEP (cm/s) and converting to m/s: '
            do i_depo=1,n_landuse_index
            do p_loop=1,n_emep_pollutant_loop

                i_pollutant=pollutant_loop_index(p_loop)
                var_name_nc_temp=deposition_name_nc(i_depo,i_pollutant)
                    
                status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                if (status_nc.eq.NF90_NOERR) then
                    status_nc = NF90_GET_VAR (id_nc, var_id_nc, depo_var3d_nc(:,:,:,i_depo,p_loop),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                    write(unit_logfile,'(A,2A,2f16.4)') ' Reading deposition velocity: ',trim(var_name_nc_temp),' (min, max): ',minval(depo_var3d_nc(:,:,:,i_depo,p_loop)),maxval(depo_var3d_nc(:,:,:,i_depo,p_loop))
                    !write(unit_logfile,'(2A,f16.4)') ' Average of: ',trim(var_name_nc_temp),sum(species_temp_var3d_nc(:,:,:))/(size(species_temp_var3d_nc,1)*size(species_temp_var3d_nc,2)*size(species_temp_var3d_nc,3))
                else
                    !write(unit_logfile,'(8A,8A)') ' Cannot read deposition velocity: ',trim(var_name_nc_temp)
                endif    
            enddo    
            enddo
            
            !Converting to m/s
            depo_var3d_nc=depo_var3d_nc/100.
            
        endif
        
        !Loop through the species
        if ((save_emep_species.or.save_seasalt).and.i_file.eq.1) then
            write(unit_logfile,'(A)') ' Reading species data from EMEP: '

            do i_sp=1,n_species_loop_index
                ii_sp=species_loop_index(i_sp)
                if (ii_sp.eq.sp_soa_index) then
                    !Read a and b soa
                    species_temp_var3d_nc=0.
                    pmxx_sp_index=pm25_sp_index
                    var_name_nc_temp=species_name_nc(pmxx_sp_index,sp_asoa_in_index)
                    status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                    if (status_nc.eq.NF90_NOERR) then
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, species_temp_var3d_nc(:,:,:),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,2A,2f16.4)') ' Reading species: ',trim(var_name_nc_temp),' (min, max): ',minval(species_temp_var3d_nc(:,:,:)),maxval(species_temp_var3d_nc(:,:,:))
                        species_var3d_nc(:,:,:,pmxx_sp_index,sp_asoa_index)=species_temp_var3d_nc(:,:,:)
                        species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)=species_temp_var3d_nc(:,:,:)
                    !write(unit_logfile,'(2A,f16.4)') ' Average of: ',trim(var_name_nc_temp),sum(species_temp_var3d_nc(:,:,:))/(size(species_temp_var3d_nc,1)*size(species_temp_var3d_nc,2)*size(species_temp_var3d_nc,3))
                    else
                         write(unit_logfile,'(8A,8A)') ' Cannot read species: ',trim(var_name_nc_temp)
                    endif    
                    var_name_nc_temp=species_name_nc(pmxx_sp_index,sp_bsoa_in_index)
                    status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                    if (status_nc.eq.NF90_NOERR) then
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, species_temp_var3d_nc(:,:,:),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,2A,2f16.4)') ' Reading species: ',trim(var_name_nc_temp),' (min, max): ',minval(species_temp_var3d_nc(:,:,:)),maxval(species_temp_var3d_nc(:,:,:))
                        species_var3d_nc(:,:,:,pmxx_sp_index,sp_bsoa_index)=species_temp_var3d_nc(:,:,:)
                        species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)=species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)+species_temp_var3d_nc(:,:,:)
                        write(unit_logfile,'(A,2A,2f16.4)') ' Adding species: ',trim(species_name_nc(pmxx_sp_index,ii_sp)),' (min, max): ',minval(species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)),maxval(species_var3d_nc(:,:,:,pm25_sp_index,i_sp))
                    !write(unit_logfile,'(2A,f16.4)') ' Average of: ',trim(var_name_nc_temp),sum(species_temp_var3d_nc(:,:,:))/(size(species_temp_var3d_nc,1)*size(species_temp_var3d_nc,2)*size(species_temp_var3d_nc,3))
                    else
                         write(unit_logfile,'(8A,8A)') ' Cannot read species: ',trim(var_name_nc_temp)
                    endif    
                    species_var3d_nc(:,:,:,pm10_sp_index,i_sp)=species_var3d_nc(:,:,:,pm25_sp_index,i_sp)
                
                endif
                
                if (ii_sp.eq.sp_sia_index) then
                    !Read pm10 sia but currently not using. Use the sum of PM2.5 + 0.73*NO3_course instead. This is done at the end of the sia reading
                    pmxx_sp_index=pm10_sp_index
                    var_name_nc_temp=species_name_nc(pmxx_sp_index,sp_sia_in_index)
                    status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                    if (status_nc.eq.NF90_NOERR) then
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, species_temp_var3d_nc(:,:,:),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,2A,2f16.4)') ' Reading species: ',trim(var_name_nc_temp),' (min, max): ',minval(species_temp_var3d_nc(:,:,:)),maxval(species_temp_var3d_nc(:,:,:))
                        species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)=species_temp_var3d_nc(:,:,:)
                        !write(unit_logfile,'(2A,f16.4)') ' Average of: ',trim(var_name_nc_temp),sum(species_temp_var3d_nc(:,:,:))/(size(species_temp_var3d_nc,1)*size(species_temp_var3d_nc,2)*size(species_temp_var3d_nc,3))
                    else
                         write(unit_logfile,'(8A,8A)') ' Cannot read species: ',trim(var_name_nc_temp)
                    endif    

                    !Read course NO3
                    pmxx_sp_index=pmco_sp_index
                    var_name_nc_temp=species_name_nc(pmxx_sp_index,sp_no3_index)
                    status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                    if (status_nc.eq.NF90_NOERR) then
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, species_temp_var3d_nc(:,:,:),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,2A,2f16.4)') ' Reading species: ',trim(var_name_nc_temp),' (min, max): ',minval(species_temp_var3d_nc(:,:,:)),maxval(species_temp_var3d_nc(:,:,:))
                        species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)=species_temp_var3d_nc(:,:,:)
                        !write(unit_logfile,'(2A,f16.4)') ' Average of: ',trim(var_name_nc_temp),sum(species_temp_var3d_nc(:,:,:))/(size(species_temp_var3d_nc,1)*size(species_temp_var3d_nc,2)*size(species_temp_var3d_nc,3))
                    else
                         write(unit_logfile,'(8A,8A)') ' Cannot read species: ',trim(var_name_nc_temp)
                    endif    

                    !Set to 0
                    species_var3d_nc(:,:,:,pm25_sp_index,sp_sia_index)=0.
                    !Read pm25(fine) no3
                    pmxx_sp_index=pm25_sp_index
                    var_name_nc_temp=species_name_nc(pmxx_sp_index,sp_no3_index)
                    status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                    if (status_nc.eq.NF90_NOERR) then
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, species_temp_var3d_nc(:,:,:),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,2A,2f16.4)') ' Reading species: ',trim(var_name_nc_temp),' (min, max): ',minval(species_temp_var3d_nc(:,:,:)),maxval(species_temp_var3d_nc(:,:,:))
                        species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)=species_temp_var3d_nc(:,:,:)
                        !write(unit_logfile,'(2A,f16.4)') ' Average of: ',trim(var_name_nc_temp),sum(species_temp_var3d_nc(:,:,:))/(size(species_temp_var3d_nc,1)*size(species_temp_var3d_nc,2)*size(species_temp_var3d_nc,3))
                    else
                         write(unit_logfile,'(8A,8A)') ' Cannot read species: ',trim(var_name_nc_temp)
                    endif
                    
                    !Read so4 and add to no3
                    pmxx_sp_index=pm25_sp_index
                    var_name_nc_temp=species_name_nc(pmxx_sp_index,sp_so4_index)
                    status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                    if (status_nc.eq.NF90_NOERR) then
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, species_temp_var3d_nc(:,:,:),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,2A,2f16.4)') ' Reading species: ',trim(var_name_nc_temp),' (min, max): ',minval(species_temp_var3d_nc(:,:,:)),maxval(species_temp_var3d_nc(:,:,:))
                        species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)=species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)+species_temp_var3d_nc(:,:,:)
                        !write(unit_logfile,'(2A,f16.4)') ' Average of: ',trim(var_name_nc_temp),sum(species_temp_var3d_nc(:,:,:))/(size(species_temp_var3d_nc,1)*size(species_temp_var3d_nc,2)*size(species_temp_var3d_nc,3))
                    else
                         write(unit_logfile,'(8A,8A)') ' Cannot read species: ',trim(var_name_nc_temp)
                    endif  
                    
                    !Read nh4 and add to so4 and no3
                    pmxx_sp_index=pm25_sp_index
                    var_name_nc_temp=species_name_nc(pmxx_sp_index,sp_nh4_index)
                    status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                    if (status_nc.eq.NF90_NOERR) then
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, species_temp_var3d_nc(:,:,:),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,2A,2f16.4)') ' Reading species: ',trim(var_name_nc_temp),' (min, max): ',minval(species_temp_var3d_nc(:,:,:)),maxval(species_temp_var3d_nc(:,:,:))
                        species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)=species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)+species_temp_var3d_nc(:,:,:)
                        !write(unit_logfile,'(2A,f16.4)') ' Average of: ',trim(var_name_nc_temp),sum(species_temp_var3d_nc(:,:,:))/(size(species_temp_var3d_nc,1)*size(species_temp_var3d_nc,2)*size(species_temp_var3d_nc,3))
                   else
                         write(unit_logfile,'(8A,8A)') ' Cannot read species: ',trim(var_name_nc_temp)
                    endif    
                    
                    !Add up the species
                    species_var3d_nc(:,:,:,pm25_sp_index,i_sp)=species_var3d_nc(:,:,:,pm25_sp_index,i_sp) !+ 0.27*species_var3d_nc(:,:,:,pmco_sp_index,i_sp) !Not PM2.5 but PMfine
                    species_var3d_nc(:,:,:,pm10_sp_index,i_sp)=species_var3d_nc(:,:,:,pm25_sp_index,i_sp) +  species_var3d_nc(:,:,:,pmco_sp_index,i_sp) !0.73*species_var3d_nc(:,:,:,pmco_sp_index,i_sp)
                    write(unit_logfile,'(A,2A,2f16.4)') ' Adding species: ',trim(species_name_nc(pm25_sp_index,ii_sp)),' (min, max): ',minval(species_var3d_nc(:,:,:,pm25_sp_index,i_sp)),maxval(species_var3d_nc(:,:,:,pm25_sp_index,i_sp))
                    write(unit_logfile,'(A,2A,2f16.4)') ' Adding species: ',trim(species_name_nc(pm10_sp_index,ii_sp)),' (min, max): ',minval(species_var3d_nc(:,:,:,pm10_sp_index,i_sp)),maxval(species_var3d_nc(:,:,:,pm10_sp_index,i_sp))
                    !write(unit_logfile,'(2A,f16.4)') ' Average of: ',trim(species_name_nc(pm25_sp_index,i_sp)),sum(species_var3d_nc(:,:,:,pm25_sp_index,i_sp))/(size(species_temp_var3d_nc,1)*size(species_temp_var3d_nc,2)*size(species_temp_var3d_nc,3))
                    !write(unit_logfile,'(2A,f16.4)') ' Average of: ',trim(species_name_nc(pm10_sp_index,i_sp)),sum(species_var3d_nc(:,:,:,pm10_sp_index,i_sp))/(size(species_temp_var3d_nc,1)*size(species_temp_var3d_nc,2)*size(species_temp_var3d_nc,3))
                endif

            
               if (ii_sp.eq.sp_dust_index) then
                    !Read dust
                    do pmxx_sp_index=1,n_pmxx_sp_index
                    if (pmxx_sp_index.eq.pm25_sp_index.or.pmco_sp_index.eq.pmxx_sp_index) then
                    var_name_nc_temp=species_name_nc(pmxx_sp_index,sp_dust_sah_index)
                    status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                    if (status_nc.eq.NF90_NOERR) then
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, species_temp_var3d_nc(:,:,:),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,2A,2f16.4)') ' Reading species: ',trim(var_name_nc_temp),' (min, max): ',minval(species_temp_var3d_nc(:,:,:)),maxval(species_temp_var3d_nc(:,:,:))
                        species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)=species_temp_var3d_nc(:,:,:)
                    else
                         write(unit_logfile,'(8A,8A)') ' Cannot read species: ',trim(var_name_nc_temp)
                    endif    
                    var_name_nc_temp=species_name_nc(pmxx_sp_index,sp_dust_wb_index)
                    status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                    if (status_nc.eq.NF90_NOERR) then
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, species_temp_var3d_nc(:,:,:),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,2A,2f16.4)') ' Reading species: ',trim(var_name_nc_temp),' (min, max): ',minval(species_temp_var3d_nc(:,:,:)),maxval(species_temp_var3d_nc(:,:,:))
                        species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)=species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)+species_temp_var3d_nc(:,:,:)
                        write(unit_logfile,'(A,2A,2f16.4)') ' Adding species: ',trim(species_name_nc(pmxx_sp_index,ii_sp)),' (min, max): ',minval(species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)),maxval(species_var3d_nc(:,:,:,pmxx_sp_index,i_sp))
                    else
                         write(unit_logfile,'(8A,8A)') ' Cannot read species: ',trim(var_name_nc_temp)
                    endif
                    endif
                    enddo
                    species_var3d_nc(:,:,:,pm10_sp_index,i_sp)=species_var3d_nc(:,:,:,pm25_sp_index,i_sp)+species_var3d_nc(:,:,:,pmco_sp_index,i_sp)
                    write(unit_logfile,'(A,2A,2f16.4)') ' Adding species: ',trim(species_name_nc(pm10_sp_index,ii_sp)),' (min, max): ',minval(species_var3d_nc(:,:,:,pm10_sp_index,i_sp)),maxval(species_var3d_nc(:,:,:,pm10_sp_index,i_sp))
                    
               endif
  
               if (ii_sp.eq.sp_seasalt_index) then
                    !Read seasalt
                    do pmxx_sp_index=1,n_pmxx_sp_index
                    if (pmxx_sp_index.eq.pm25_sp_index.or.pmco_sp_index.eq.pmxx_sp_index) then
                    var_name_nc_temp=species_name_nc(pmxx_sp_index,sp_seasalt_in_index)
                    status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                    if (status_nc.eq.NF90_NOERR) then
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, species_temp_var3d_nc(:,:,:),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,2A,2f16.4)') ' Reading species: ',trim(var_name_nc_temp),' (min, max): ',minval(species_temp_var3d_nc(:,:,:)),maxval(species_temp_var3d_nc(:,:,:))
                        species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)=species_temp_var3d_nc(:,:,:)
                    else
                         write(unit_logfile,'(8A,8A)') ' Cannot read species: ',trim(var_name_nc_temp)
                    endif    
                    endif
                    enddo
                    species_var3d_nc(:,:,:,pm10_sp_index,i_sp)=species_var3d_nc(:,:,:,pm25_sp_index,i_sp)+species_var3d_nc(:,:,:,pmco_sp_index,i_sp)
                    write(unit_logfile,'(A,2A,2f16.4)') ' Adding species: ',trim(species_name_nc(pm10_sp_index,ii_sp)),' (min, max): ',minval(species_var3d_nc(:,:,:,pm10_sp_index,i_sp)),maxval(species_var3d_nc(:,:,:,pm10_sp_index,i_sp))
                    
               endif
            
               if (ii_sp.eq.sp_ffire_index) then
                    !Read fire
                    pmxx_sp_index=pm25_sp_index
                    var_name_nc_temp=species_name_nc(pmxx_sp_index,sp_ffire_bc_index)
                    status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                    if (status_nc.eq.NF90_NOERR) then
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, species_temp_var3d_nc(:,:,:),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,2A,2f16.4)') ' Reading species: ',trim(var_name_nc_temp),' (min, max): ',minval(species_temp_var3d_nc(:,:,:)),maxval(species_temp_var3d_nc(:,:,:))
                        species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)=species_temp_var3d_nc(:,:,:)
                    else
                         write(unit_logfile,'(8A,8A)') ' Cannot read species: ',trim(var_name_nc_temp)
                    endif    
                    var_name_nc_temp=species_name_nc(pmxx_sp_index,sp_ffire_rem_index)
                    status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                    if (status_nc.eq.NF90_NOERR) then
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, species_temp_var3d_nc(:,:,:),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,2A,2f16.4)') ' Reading species: ',trim(var_name_nc_temp),' (min, max): ',minval(species_temp_var3d_nc(:,:,:)),maxval(species_temp_var3d_nc(:,:,:))
                        species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)=species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)+species_temp_var3d_nc(:,:,:)
                        write(unit_logfile,'(A,2A,2f16.4)') ' Adding species: ',trim(species_name_nc(pmxx_sp_index,ii_sp)),' (min, max): ',minval(species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)),maxval(species_var3d_nc(:,:,:,pmxx_sp_index,i_sp))
                    else
                         write(unit_logfile,'(8A,8A)') ' Cannot read species: ',trim(var_name_nc_temp)
                    endif
                    species_var3d_nc(:,:,:,pm10_sp_index,i_sp)=species_var3d_nc(:,:,:,pm25_sp_index,i_sp)
                    
               endif

               if (ii_sp.eq.sp_ppm_index) then
                    !Read ppm
                    do pmxx_sp_index=1,n_pmxx_sp_index
                    if (pmxx_sp_index.eq.pm25_sp_index.or.pm10_sp_index.eq.pmxx_sp_index) then
                    var_name_nc_temp=species_name_nc(pmxx_sp_index,sp_ppm_in_index)
                    status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                    if (status_nc.eq.NF90_NOERR) then
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, species_temp_var3d_nc(:,:,:),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,2A,2f16.4)') ' Reading species: ',trim(var_name_nc_temp),' (min, max): ',minval(species_temp_var3d_nc(:,:,:)),maxval(species_temp_var3d_nc(:,:,:))
                        species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)=species_temp_var3d_nc(:,:,:)
                        !write(unit_logfile,'(2A,f16.4)') ' Average of: ',trim(var_name_nc_temp),sum(species_temp_var3d_nc(:,:,:))/(size(species_temp_var3d_nc,1)*size(species_temp_var3d_nc,2)*size(species_temp_var3d_nc,3))
                    else
                         write(unit_logfile,'(8A,8A)') ' Cannot read species: ',trim(var_name_nc_temp)
                    endif    
                    endif
                    enddo
                    
               endif
               
               if (ii_sp.eq.sp_water_index) then
                    !Read water pm25 only
                    pmxx_sp_index=pm25_sp_index
                    var_name_nc_temp=species_name_nc(pmxx_sp_index,sp_water_in_index)
                    status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                    if (status_nc.eq.NF90_NOERR) then
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, species_temp_var3d_nc(:,:,:),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,2A,2f16.4)') ' Reading species: ',trim(var_name_nc_temp),' (min, max): ',minval(species_temp_var3d_nc(:,:,:)),maxval(species_temp_var3d_nc(:,:,:))
                        species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)=species_temp_var3d_nc(:,:,:)
                    else
                         write(unit_logfile,'(8A,8A)') ' Cannot read species: ',trim(var_name_nc_temp)
                    endif
                    !Set water in pm10 the same as pm25
                    species_var3d_nc(:,:,:,pm10_sp_index,i_sp)=species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)
               endif

               !Should there be 0.27*NO3 here for PM2.5? Check naming!!
               if (ii_sp.eq.sp_pm_index) then
                    !Read total pm
                    do pmxx_sp_index=1,n_pmxx_sp_index
                    !if (pmxx_sp_index.eq.pm25_sp_index.or.pm10_sp_index.eq.pmxx_sp_index) then
                    var_name_nc_temp=species_name_nc(pmxx_sp_index,sp_pm_in_index)
                    status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                    if (status_nc.eq.NF90_NOERR) then
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, species_temp_var3d_nc(:,:,:),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,2A,2f16.4)') ' Reading species: ',trim(var_name_nc_temp),' (min, max): ',minval(species_temp_var3d_nc(:,:,:)),maxval(species_temp_var3d_nc(:,:,:))
                        species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)=species_temp_var3d_nc(:,:,:)
                        !write(unit_logfile,'(2A,f16.4)') ' Average of: ',trim(var_name_nc_temp),sum(species_temp_var3d_nc(:,:,:))/(size(species_temp_var3d_nc,1)*size(species_temp_var3d_nc,2)*size(species_temp_var3d_nc,3))
                    else
                         write(unit_logfile,'(8A,8A)') ' Cannot read species: ',trim(var_name_nc_temp)
                    endif    
                    !endif
                    enddo
                    
               endif

            if (save_emep_OP_species) then
               
                if (ii_sp.eq.sp_BBOA_index) then
                    !Read total pm
                    do pmxx_sp_index=1,n_pmxx_sp_index
                    !if (pmxx_sp_index.eq.pm25_sp_index.or.pm10_sp_index.eq.pmxx_sp_index) then
                    var_name_nc_temp=species_name_nc(pmxx_sp_index,sp_FFIRE_OM_in_index)
                    status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                    if (status_nc.eq.NF90_NOERR) then
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, species_temp_var3d_nc(:,:,:),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,2A,2f16.4)') ' Reading species: ',trim(var_name_nc_temp),' (min, max): ',minval(species_temp_var3d_nc(:,:,:)),maxval(species_temp_var3d_nc(:,:,:))
                        species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)=species_temp_var3d_nc(:,:,:)
                        !write(unit_logfile,'(2A,f16.4)') ' Average of: ',trim(var_name_nc_temp),sum(species_temp_var3d_nc(:,:,:))/(size(species_temp_var3d_nc,1)*size(species_temp_var3d_nc,2)*size(species_temp_var3d_nc,3))
                    else
                         write(unit_logfile,'(8A,8A)') ' Cannot read species: ',trim(var_name_nc_temp)
                    endif    
                    !endif
                    enddo
                    do pmxx_sp_index=1,n_pmxx_sp_index
                    !if (pmxx_sp_index.eq.pm25_sp_index.or.pm10_sp_index.eq.pmxx_sp_index) then
                    var_name_nc_temp=species_name_nc(pmxx_sp_index,sp_FFIRE_BC_in_index)
                    status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                    if (status_nc.eq.NF90_NOERR) then
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, species_temp_var3d_nc(:,:,:),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,2A,2f16.4)') ' Reading species: ',trim(var_name_nc_temp),' (min, max): ',minval(species_temp_var3d_nc(:,:,:)),maxval(species_temp_var3d_nc(:,:,:))
                        species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)=species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)+species_temp_var3d_nc(:,:,:)
                        !write(unit_logfile,'(2A,f16.4)') ' Average of: ',trim(var_name_nc_temp),sum(species_temp_var3d_nc(:,:,:))/(size(species_temp_var3d_nc,1)*size(species_temp_var3d_nc,2)*size(species_temp_var3d_nc,3))
                    else
                         write(unit_logfile,'(8A,8A)') ' Cannot read species: ',trim(var_name_nc_temp)
                    endif    
                    !endif
                    enddo
                    do pmxx_sp_index=1,n_pmxx_sp_index
                    !if (pmxx_sp_index.eq.pm25_sp_index.or.pm10_sp_index.eq.pmxx_sp_index) then
                    var_name_nc_temp=species_name_nc(pmxx_sp_index,sp_FFIRE_REM_in_index)
                    status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                    if (status_nc.eq.NF90_NOERR) then
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, species_temp_var3d_nc(:,:,:),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,2A,2f16.4)') ' Reading species: ',trim(var_name_nc_temp),' (min, max): ',minval(species_temp_var3d_nc(:,:,:)),maxval(species_temp_var3d_nc(:,:,:))
                        species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)=species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)+species_temp_var3d_nc(:,:,:)
                        !write(unit_logfile,'(2A,f16.4)') ' Average of: ',trim(var_name_nc_temp),sum(species_temp_var3d_nc(:,:,:))/(size(species_temp_var3d_nc,1)*size(species_temp_var3d_nc,2)*size(species_temp_var3d_nc,3))
                    else
                         write(unit_logfile,'(8A,8A)') ' Cannot read species: ',trim(var_name_nc_temp)
                    endif    
                    !endif
                    enddo
                    species_var3d_nc(:,:,:,pm10_sp_index,i_sp)=species_var3d_nc(:,:,:,pm25_sp_index,i_sp)+species_var3d_nc(:,:,:,pmco_sp_index,i_sp)
                    write(unit_logfile,'(A,2A,2f16.4)') ' Adding species: ',trim(species_name_nc(pm10_sp_index,ii_sp)),' (min, max): ',minval(species_var3d_nc(:,:,:,pm10_sp_index,i_sp)),maxval(species_var3d_nc(:,:,:,pm10_sp_index,i_sp))
                endif

                !write(*,*) i_sp,ii_sp,n_species_loop_index,size(species_var3d_nc,5)
                
                if (ii_sp.eq.sp_BBOA_RES_index) then
                    !Read total pm
                    do pmxx_sp_index=1,n_pmxx_sp_index
                    !if (pmxx_sp_index.eq.pm25_sp_index.or.pm10_sp_index.eq.pmxx_sp_index) then
                    var_name_nc_temp=species_name_nc(pmxx_sp_index,sp_POM_RES_in_index)
                    status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                    if (status_nc.eq.NF90_NOERR) then
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, species_temp_var3d_nc(:,:,:),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,2A,2f16.4)') ' Reading species: ',trim(var_name_nc_temp),' (min, max): ',minval(species_temp_var3d_nc(:,:,:)),maxval(species_temp_var3d_nc(:,:,:))
                        species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)=species_temp_var3d_nc(:,:,:)
                        !write(unit_logfile,'(2A,f16.4)') ' Average of: ',trim(var_name_nc_temp),sum(species_temp_var3d_nc(:,:,:))/(size(species_temp_var3d_nc,1)*size(species_temp_var3d_nc,2)*size(species_temp_var3d_nc,3))
                    else
                         write(unit_logfile,'(8A,8A)') ' Cannot read species: ',trim(var_name_nc_temp)
                    endif    
                    !endif
                    enddo
                    do pmxx_sp_index=1,n_pmxx_sp_index
                    !if (pmxx_sp_index.eq.pm25_sp_index.or.pm10_sp_index.eq.pmxx_sp_index) then
                    var_name_nc_temp=species_name_nc(pmxx_sp_index,sp_EC_RES_NEW_in_index)
                    status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                    if (status_nc.eq.NF90_NOERR) then
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, species_temp_var3d_nc(:,:,:),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,2A,2f16.4)') ' Reading species: ',trim(var_name_nc_temp),' (min, max): ',minval(species_temp_var3d_nc(:,:,:)),maxval(species_temp_var3d_nc(:,:,:))
                        species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)=species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)+species_temp_var3d_nc(:,:,:)
                        !write(unit_logfile,'(2A,f16.4)') ' Average of: ',trim(var_name_nc_temp),sum(species_temp_var3d_nc(:,:,:))/(size(species_temp_var3d_nc,1)*size(species_temp_var3d_nc,2)*size(species_temp_var3d_nc,3))
                    else
                         write(unit_logfile,'(8A,8A)') ' Cannot read species: ',trim(var_name_nc_temp)
                    endif    
                    !endif
                    enddo
                    do pmxx_sp_index=1,n_pmxx_sp_index
                    !if (pmxx_sp_index.eq.pm25_sp_index.or.pm10_sp_index.eq.pmxx_sp_index) then
                    var_name_nc_temp=species_name_nc(pmxx_sp_index,sp_EC_RES_AGE_in_index)
                    status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                    if (status_nc.eq.NF90_NOERR) then
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, species_temp_var3d_nc(:,:,:),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,2A,2f16.4)') ' Reading species: ',trim(var_name_nc_temp),' (min, max): ',minval(species_temp_var3d_nc(:,:,:)),maxval(species_temp_var3d_nc(:,:,:))
                        species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)=species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)+species_temp_var3d_nc(:,:,:)
                        !write(unit_logfile,'(2A,f16.4)') ' Average of: ',trim(var_name_nc_temp),sum(species_temp_var3d_nc(:,:,:))/(size(species_temp_var3d_nc,1)*size(species_temp_var3d_nc,2)*size(species_temp_var3d_nc,3))
                    else
                         write(unit_logfile,'(8A,8A)') ' Cannot read species: ',trim(var_name_nc_temp)
                    endif    
                    !endif
                    enddo
                    do pmxx_sp_index=1,n_pmxx_sp_index
                    !if (pmxx_sp_index.eq.pm25_sp_index.or.pm10_sp_index.eq.pmxx_sp_index) then
                    var_name_nc_temp=species_name_nc(pmxx_sp_index,sp_REM_RES_in_index)
                    status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                    if (status_nc.eq.NF90_NOERR) then
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, species_temp_var3d_nc(:,:,:),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,2A,2f16.4)') ' Reading species: ',trim(var_name_nc_temp),' (min, max): ',minval(species_temp_var3d_nc(:,:,:)),maxval(species_temp_var3d_nc(:,:,:))
                        species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)=species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)+species_temp_var3d_nc(:,:,:)
                        !write(unit_logfile,'(2A,f16.4)') ' Average of: ',trim(var_name_nc_temp),sum(species_temp_var3d_nc(:,:,:))/(size(species_temp_var3d_nc,1)*size(species_temp_var3d_nc,2)*size(species_temp_var3d_nc,3))
                    else
                         write(unit_logfile,'(8A,8A)') ' Cannot read species: ',trim(var_name_nc_temp)
                    endif    
                    !endif
                    enddo
                     do pmxx_sp_index=1,n_pmxx_sp_index
                    !if (pmxx_sp_index.eq.pm25_sp_index.or.pm10_sp_index.eq.pmxx_sp_index) then
                    var_name_nc_temp=species_name_nc(pmxx_sp_index,sp_EC_RES_in_index)
                    status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                    if (status_nc.eq.NF90_NOERR) then
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, species_temp_var3d_nc(:,:,:),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,2A,2f16.4)') ' Reading species: ',trim(var_name_nc_temp),' (min, max): ',minval(species_temp_var3d_nc(:,:,:)),maxval(species_temp_var3d_nc(:,:,:))
                        species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)=species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)+species_temp_var3d_nc(:,:,:)
                        !write(unit_logfile,'(2A,f16.4)') ' Average of: ',trim(var_name_nc_temp),sum(species_temp_var3d_nc(:,:,:))/(size(species_temp_var3d_nc,1)*size(species_temp_var3d_nc,2)*size(species_temp_var3d_nc,3))
                    else
                         write(unit_logfile,'(8A,8A)') ' Cannot read species: ',trim(var_name_nc_temp)
                    endif    
                    !endif
                    enddo
                    do pmxx_sp_index=1,n_pmxx_sp_index
                    !if (pmxx_sp_index.eq.pm25_sp_index.or.pm10_sp_index.eq.pmxx_sp_index) then
                    var_name_nc_temp=species_name_nc(pmxx_sp_index,sp_EC_RES_in_index)
                    status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                    if (status_nc.eq.NF90_NOERR) then
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, species_temp_var3d_nc(:,:,:),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,2A,2f16.4)') ' Reading species: ',trim(var_name_nc_temp),' (min, max): ',minval(species_temp_var3d_nc(:,:,:)),maxval(species_temp_var3d_nc(:,:,:))
                        species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)=species_var3d_nc(:,:,:,pmxx_sp_index,i_sp)+species_temp_var3d_nc(:,:,:)
                        !write(unit_logfile,'(2A,f16.4)') ' Average of: ',trim(var_name_nc_temp),sum(species_temp_var3d_nc(:,:,:))/(size(species_temp_var3d_nc,1)*size(species_temp_var3d_nc,2)*size(species_temp_var3d_nc,3))
                    else
                         write(unit_logfile,'(8A,8A)') ' Cannot read species: ',trim(var_name_nc_temp)
                    endif    
                    !endif
                    enddo
                   
                     species_var3d_nc(:,:,:,pm10_sp_index,i_sp)=species_var3d_nc(:,:,:,pm25_sp_index,i_sp)+species_var3d_nc(:,:,:,pmco_sp_index,i_sp)
                     write(unit_logfile,'(A,2A,2f16.4)') ' Adding species: ',trim(species_name_nc(pm10_sp_index,ii_sp)),' (min, max): ',minval(species_var3d_nc(:,:,:,pm10_sp_index,i_sp)),maxval(species_var3d_nc(:,:,:,pm10_sp_index,i_sp))
                   
               endif

                
            endif
            
        !species_name_nc(pm10_sp_index,sp_BBOA_index)='pm10_EMEP_BBOA'
        !species_name_nc(pm25_sp_index,sp_BBOA_index)='pm25_EMEP_BBOA'
        !species_name_nc(pmco_sp_index,sp_BBOA_index)='pmco_EMEP_BBOA'
        !species_name_nc(pm10_sp_index,sp_BBOA_RES_index)='pm10_EMEP_BBOA_RES'
        !species_name_nc(pm25_sp_index,sp_BBOA_RES_index)='pm25_EMEP_BBOA_RES'
        !species_name_nc(pmco_sp_index,sp_BBOA_RES_index)='pmco_EMEP_BBOA_RES'
        
        !species_name_nc(pm25_sp_index,sp_POM_RES_in_index)='SURF_ug_POM_F_RES'
        !species_name_nc(pm25_sp_index,sp_EC_RES_NEW_in_index)='SURF_ug_EC_F_RES_NEW'
        !species_name_nc(pm25_sp_index,sp_EC_RES_AGE_in_index)='SURF_ug_EC_F_RES_AGE'
        !species_name_nc(pm25_sp_index,sp_REM_RES_in_index)='SURF_ug_REMPPM25_RES'
        !species_name_nc(pm25_sp_index,sp_FFIRE_OM_in_index)='SURF_ug_FFIRE_OM'
        !species_name_nc(pm25_sp_index,sp_FFIRE_BC_in_index)='SURF_ug_FFIRE_BC'
        !species_name_nc(pm25_sp_index,sp_FFIRE_REM_in_index)='SURF_ug_FFIRE_REMPPM25'
        
        !species_name_nc(pmco_sp_index,sp_EC_RES_in_index)='SURF_ug_EC_C_RES'
        !species_name_nc(pmco_sp_index,sp_POM_RES_in_index)='SURF_ug_POM_C_RES'
        !species_name_nc(pmco_sp_index,sp_REM_RES_in_index)='SURF_ug_REMPPM_C_RES'
        !species_name_nc(pmco_sp_index,sp_FFIRE_in_index)='SURF_ug_FFIRE_C'

            enddo !sp_index
            
            !Derive SOA from the other species. Based on using PMFINE as the pm25 total 
            if (save_emep_species.and.derive_SOA_from_other_species) then
                species_var3d_nc(:,:,:,pm25_sp_index,sp_soa_index)= &
                     species_var3d_nc(:,:,:,pm25_sp_index,sp_pm_index) &
                    -species_var3d_nc(:,:,:,pm25_sp_index,sp_sia_index) &
                    -species_var3d_nc(:,:,:,pm25_sp_index,sp_seasalt_index) &
                    -species_var3d_nc(:,:,:,pm25_sp_index,sp_dust_index) &
                    -species_var3d_nc(:,:,:,pm25_sp_index,sp_ffire_index) &
                    -species_var3d_nc(:,:,:,pm25_sp_index,sp_ppm_index)
                
                species_var3d_nc(:,:,:,pm10_sp_index,sp_soa_index)=species_var3d_nc(:,:,:,pm25_sp_index,sp_soa_index)
                species_var3d_nc(:,:,:,pm10_sp_index,sp_asoa_index)=species_var3d_nc(:,:,:,pm25_sp_index,sp_asoa_index)
                species_var3d_nc(:,:,:,pm10_sp_index,sp_bsoa_index)=species_var3d_nc(:,:,:,pm25_sp_index,sp_bsoa_index)
                
                var_name_nc_temp=species_name_nc(pm25_sp_index,sp_soa_index)

                write(unit_logfile,'(A,2A,2f16.4)') ' Calculating species: ',trim(var_name_nc_temp),' (min, max): ',minval(species_var3d_nc(:,:,:,pm25_sp_index,sp_soa_index)),maxval(species_var3d_nc(:,:,:,pm25_sp_index,sp_soa_index))
                
                where (species_var3d_nc(:,:,:,:,sp_soa_index).le.0) species_var3d_nc(:,:,:,:,sp_soa_index)=0
                
                write(unit_logfile,'(A,2A,2f16.4)') ' Limitting species: ',trim(var_name_nc_temp),' (min, max): ',minval(species_var3d_nc(:,:,:,pm25_sp_index,sp_soa_index)),maxval(species_var3d_nc(:,:,:,pm25_sp_index,sp_soa_index))
               
                
                !SURF_ug_SOA=SURF_ug_PMFINE- SURF_ug_SO4- SURF_ug_NO3_F- SURF_ug_NH4_F- SURF_ug_SEASALT_F- SURF_ug_DUST_SAH_F- SURF_ug_DUST_WB_F- SURF_ug_FFIRE_REMPPM25- SURF_ug_FFIRE_BC- SURF_ug_PPM2.5                
            endif
            
            !Calculate the total of the parts by adding all and subtracting the total
            !Do not so this when only using sea salt
            if (save_emep_species) then

                write(unit_logfile,'(A)') ' Cross check average of species'
                do i_sp=1,n_species_loop_index 
                    ii_sp=species_loop_index(i_sp)
                    if (i_sp.ne.sp_pm_index) then
                    !Do not include the total here
                        write(unit_logfile,'(A,2A,2f16.4)') ' Average of species: ',trim(species_name_nc(pm10_sp_index,ii_sp)),' (pm25, pm10): ', &
                        sum(species_var3d_nc(:,:,:,pm25_sp_index,i_sp))/(size(species_temp_var3d_nc,1)*size(species_temp_var3d_nc,2)*size(species_temp_var3d_nc,3)), &
                        sum(species_var3d_nc(:,:,:,pm10_sp_index,i_sp))/(size(species_temp_var3d_nc,1)*size(species_temp_var3d_nc,2)*size(species_temp_var3d_nc,3))
                    endif            
                enddo
            
                 species_temp_var3d_nc(:,:,:)=sum(species_var3d_nc(:,:,:,pm25_sp_index,1:sp_ppm_index),4)
                
                 write(unit_logfile,'(A,2A,4f16.4)') ' Average of: ','pm25',' (sum species, total pm, D3, SURF): ', &
                    sum(species_temp_var3d_nc(:,:,:))/(size(species_temp_var3d_nc,1)*size(species_temp_var3d_nc,2)*size(species_temp_var3d_nc,3)), &
                    sum(species_var3d_nc(:,:,:,pm25_sp_index,sp_pm_index))/(size(species_temp_var3d_nc,1)*size(species_temp_var3d_nc,2)*size(species_temp_var3d_nc,3)), &
                    sum(comp_var4d_nc(:,:,surface_level_nc,:,pm25_nc_index))/(size(species_temp_var3d_nc,1)*size(species_temp_var3d_nc,2)*size(species_temp_var3d_nc,3)), &
                    sum(var3d_nc(:,:,:,conc_nc_index,allsource_index,pollutant_loop_back_index(pm25_nc_index)))/(size(species_temp_var3d_nc,1)*size(species_temp_var3d_nc,2)*size(species_temp_var3d_nc,3))
                 
                 species_temp_var3d_nc(:,:,:)=sum(species_var3d_nc(:,:,:,pm10_sp_index,1:sp_ppm_index),4)
                 
                 write(unit_logfile,'(A,2A,4f16.4)') ' Average of: ','pm10',' (sum species, total pm, D3, SURF): ', &
                    sum(species_temp_var3d_nc(:,:,:))/(size(species_temp_var3d_nc,1)*size(species_temp_var3d_nc,2)*size(species_temp_var3d_nc,3)), &
                    sum(species_var3d_nc(:,:,:,pm10_sp_index,sp_pm_index))/(size(species_temp_var3d_nc,1)*size(species_temp_var3d_nc,2)*size(species_temp_var3d_nc,3)), &
                    sum(comp_var4d_nc(:,:,surface_level_nc,:,pm10_nc_index))/(size(species_temp_var3d_nc,1)*size(species_temp_var3d_nc,2)*size(species_temp_var3d_nc,3)), &
                    sum(var3d_nc(:,:,:,conc_nc_index,allsource_index,pollutant_loop_back_index(pm10_nc_index)))/(size(species_temp_var3d_nc,1)*size(species_temp_var3d_nc,2)*size(species_temp_var3d_nc,3))
            
            endif
            
            
        endif
                  
        !Read in 2m temperature over the whole period to get the daily average for home heating
        if (use_RWC_emission_data) then
            DMT_start_time_nc_index=start_time_nc_index
            DMT_end_time_nc_index=end_time_nc_index
            DMT_dim_length_nc=DMT_end_time_nc_index-DMT_start_time_nc_index+1
            if (.not.allocated(DMT_EMEP_grid_nc)) allocate (DMT_EMEP_grid_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),DMT_dim_length_nc))
            
            if (calculate_source(heating_index).and.i_file.eq.1) then
                var_name_nc_temp=var_name_nc(t2m_nc_index,all_nc_index,allsource_nc_index)
                status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                status_nc = NF90_GET_VAR (id_nc, var_id_nc, DMT_EMEP_grid_nc,start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),DMT_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),DMT_dim_length_nc/))
                write(unit_logfile,'(3A,2f16.4)') ' Reading: ',trim(var_name_nc_temp),' (min, max): ',minval(DMT_EMEP_grid_nc),maxval(DMT_EMEP_grid_nc)
                DMT_EMEP_grid_nc(:,:,1)=sum(DMT_EMEP_grid_nc,3)/DMT_dim_length_nc-273.13
                write(unit_logfile,'(3A,2f16.4,a,i)') ' Calculating mean: ',trim(var_name_nc_temp),' (min, max): ',minval(DMT_EMEP_grid_nc(:,:,1)),maxval(DMT_EMEP_grid_nc(:,:,1)),' over this number of time steps: ',DMT_dim_length_nc
            
            endif
        
        endif

        
        status_nc = NF90_CLOSE (id_nc)
                
        
    enddo !End file loop
    
    !Set the correct time dimensions to the first file value
    dim_length_nc(time_dim_nc_index)=valid_dim_length_nc(time_dim_nc_index)
        
    !Set the grid spacing
    if (EMEP_projection_type.eq.LL_projection_index) then
        dgrid_nc(lon_nc_index)=var1d_nc(2,x_dim_nc_index)-var1d_nc(1,x_dim_nc_index)
        dgrid_nc(lat_nc_index)=var1d_nc(2,y_dim_nc_index)-var1d_nc(1,y_dim_nc_index)
        write(unit_logfile,'(A,2f16.4)') ' Grid spacing (lon,lat): ',dgrid_nc(lon_nc_index),dgrid_nc(lat_nc_index)
    else       
        dgrid_nc(lon_nc_index)=var1d_nc(2,x_dim_nc_index)-var1d_nc(1,x_dim_nc_index)
        dgrid_nc(lat_nc_index)=var1d_nc(2,y_dim_nc_index)-var1d_nc(1,y_dim_nc_index)
        write(unit_logfile,'(A,2f16.4)') ' Grid spacing (x,y) in meters: ',dgrid_nc(lon_nc_index),dgrid_nc(lat_nc_index)
        !
    endif
    
    
    !EMEP emissions for traffic do not differentiate between exhaust and nonexhaust
    !Place the read PM2.5 emissions also into the exhaust emissions. If PM10 is higher then these are then the non-exhaust emissions
    !This only needs to be done when the EMEP emissions are used and distributed using a proxy but we do it in general anyway
    !if (local_subgrid_method_flag.eq.3) then
        !write(unit_logfile,'(A,2es12.2)') ' In: Filling the EMEP exhaust emissions with PM25 emissions ',sum(var3d_nc(:,:,:,emis_nc_index,traffic_index,pollutant_loop_back_index(pmex_nc_index))),sum(var3d_nc(:,:,:,emis_nc_index,traffic_index,pollutant_loop_back_index(pm25_nc_index)))      
        !var3d_nc(:,:,:,emis_nc_index,traffic_index,pollutant_loop_back_index(pmex_nc_index))=var3d_nc(:,:,:,emis_nc_index,traffic_index,pollutant_loop_back_index(pm25_nc_index))      
        !write(unit_logfile,'(A,2es12.2)') ' Out: Filling the EMEP exhaust emissions with PM25 emissions ',sum(var3d_nc(:,:,:,emis_nc_index,traffic_index,pollutant_loop_back_index(pmex_nc_index))),sum(var3d_nc(:,:,:,emis_nc_index,traffic_index,pollutant_loop_back_index(pm25_nc_index)))      
        !var3d_nc(:,:,:,conc_nc_index,traffic_index,pollutant_loop_back_index(pmex_nc_index))=var3d_nc(:,:,:,conc_nc_index,traffic_index,pollutant_loop_back_index(pm25_nc_index))      
        !var4d_nc(:,:,:,:,conc_nc_index,traffic_index,pollutant_loop_back_index(pmex_nc_index))=var4d_nc(:,:,:,:,conc_nc_index,traffic_index,pollutant_loop_back_index(pm25_nc_index))      
        !write(unit_logfile,'(A,2es12.2)') ' 3D Filling the EMEP exhaust concentrations with PM25 emissions ',sum(var3d_nc(:,:,:,conc_nc_index,traffic_index,pollutant_loop_back_index(pmex_nc_index))),sum(var3d_nc(:,:,:,conc_nc_index,traffic_index,pollutant_loop_back_index(pm25_nc_index)))      
        !write(unit_logfile,'(A,2es12.2)') ' 4D Filling the EMEP exhaust emissions with PM25 emissions ',sum(var4d_nc(:,:,:,:,conc_nc_index,traffic_index,pollutant_loop_back_index(pmex_nc_index))),sum(var4d_nc(:,:,:,:,,traffic_index,pollutant_loop_back_index(pm25_nc_index)))      
    !endif
    if (use_alternative_ppm_variable_for_lf) then
        write(unit_logfile,'(4A)')' Replacing variable ',trim(var_name_nc(conc_nc_index,pm25_nc_index,allsource_nc_index)),' with ',trim(var_name_nc(conc_nc_index,pm25_nc_index,extrasource_nc_index))
        write(unit_logfile,'(4A)')' Replacing variable ',trim(var_name_nc(conc_nc_index,pmco_nc_index,allsource_nc_index)),' with ',trim(var_name_nc(conc_nc_index,pmco_nc_index,extrasource_nc_index))
        write(unit_logfile,'(A,f12.3)')' LF PPM2.5 before correction variables are used: ',sum(pm_var4d_nc(:,:,surface_level_nc_2,:,conc_nc_index,allsource_nc_index,2))/size(pm_var4d_nc,1)/size(pm_var4d_nc,2)/size(pm_var4d_nc,4)
        write(unit_logfile,'(A,f12.3)')' LF PPMCO  before correction variables are used: ',sum(pm_var4d_nc(:,:,surface_level_nc_2,:,conc_nc_index,allsource_nc_index,1))/size(pm_var4d_nc,1)/size(pm_var4d_nc,2)/size(pm_var4d_nc,4)
        if (alternative_ppm_variable_for_lf_dim.eq.4) then
            var4d_nc(:,:,surface_level_nc_2,:,conc_nc_index,allsource_nc_index,pollutant_loop_back_index(pm25_nc_index))=var4d_nc(:,:,surface_level_nc,:,conc_nc_index,extrasource_nc_index,pollutant_loop_back_index(pm25_nc_index))
            
            !pm_var4d_nc(:,:,surface_level_nc_2,:,conc_nc_index,allsource_nc_index,2)=var4d_nc(:,:,surface_level_nc_2,:,conc_nc_index,allsource_nc_index,pollutant_loop_back_index(pm25_nc_index))
            !pm_var4d_nc(:,:,surface_level_nc_2,:,conc_nc_index,allsource_nc_index,1)=var4d_nc(:,:,surface_level_nc,:,conc_nc_index,extrasource_nc_index,pollutant_loop_back_index(pm10_nc_index))
           
            !var4d_nc(:,:,:,:,conc_nc_index,allsource_nc_index,pmco_nc_index)=var4d_nc(:,:,:,:,conc_nc_index,extrasource_nc_index,pmco_nc_index)
           ! write(unit_logfile,'(4A)')' Replacing variable ',trim(var_name_nc(conc_nc_index,pm10_nc_index,allsource_nc_index)),' with ',trim(var_name_nc(conc_nc_index,pm10_nc_index,extrasource_nc_index))
            var4d_nc(:,:,surface_level_nc_2,:,conc_nc_index,allsource_nc_index,pollutant_loop_back_index(pm10_nc_index))=var4d_nc(:,:,surface_level_nc_2,:,conc_nc_index,allsource_nc_index,pollutant_loop_back_index(pm25_nc_index))+var4d_nc(:,:,surface_level_nc,:,conc_nc_index,extrasource_nc_index,pollutant_loop_back_index(pm10_nc_index))
        elseif (alternative_ppm_variable_for_lf_dim.eq.3) then
            var4d_nc(:,:,surface_level_nc_2,:,conc_nc_index,allsource_nc_index,pollutant_loop_back_index(pm25_nc_index))=var3d_nc(:,:,:,conc_nc_index,extrasource_nc_index,pollutant_loop_back_index(pm25_nc_index))
            
            !pm_var4d_nc(:,:,surface_level_nc_2,:,conc_nc_index,allsource_nc_index,2)=var4d_nc(:,:,surface_level_nc_2,:,conc_nc_index,allsource_nc_index,pollutant_loop_back_index(pm25_nc_index))
            !pm_var4d_nc(:,:,surface_level_nc_2,:,conc_nc_index,allsource_nc_index,1)=var3d_nc(:,:,:,conc_nc_index,extrasource_nc_index,pollutant_loop_back_index(pm10_nc_index))
            
            !var4d_nc(:,:,:,:,conc_nc_index,allsource_nc_index,pmco_nc_index)=var4d_nc(:,:,:,:,conc_nc_index,extrasource_nc_index,pmco_nc_index)
           ! write(unit_logfile,'(4A)')' Replacing variable ',trim(var_name_nc(conc_nc_index,pm10_nc_index,allsource_nc_index)),' with ',trim(var_name_nc(conc_nc_index,pm10_nc_index,extrasource_nc_index))
            !var4d_nc(:,:,surface_level_nc,:,conc_nc_index,allsource_nc_index,pollutant_loop_back_index(pm10_nc_index))=var3d_nc(:,:,:,conc_nc_index,allsource_nc_index,pollutant_loop_back_index(pm25_nc_index))+var3d_nc(:,:,:,conc_nc_index,extrasource_nc_index,pollutant_loop_back_index(pm10_nc_index))            
            !var4d_nc(:,:,surface_level_nc_2,:,conc_nc_index,allsource_nc_index,pollutant_loop_back_index(pm10_nc_index))=pm_var4d_nc(:,:,surface_level_nc_2,:,conc_nc_index,allsource_nc_index,2)+pm_var4d_nc(:,:,surface_level_nc_2,:,conc_nc_index,allsource_nc_index,1)           
            var4d_nc(:,:,surface_level_nc_2,:,conc_nc_index,allsource_nc_index,pollutant_loop_back_index(pm10_nc_index))=var3d_nc(:,:,:,conc_nc_index,extrasource_nc_index,pollutant_loop_back_index(pm10_nc_index))+var3d_nc(:,:,:,conc_nc_index,extrasource_nc_index,pollutant_loop_back_index(pm25_nc_index))
       else
            write(unit_logfile,'(4A)')' Not replacing variable ',trim(var_name_nc(conc_nc_index,pmco_nc_index,allsource_nc_index)),' with ',trim(var_name_nc(conc_nc_index,pmco_nc_index,extrasource_nc_index))            
        endif
        write(unit_logfile,'(A,f12.3)')' LF PPM2.5 after correction variables are used: ',sum(pm_var4d_nc(:,:,surface_level_nc_2,:,conc_nc_index,allsource_nc_index,2))/size(pm_var4d_nc,1)/size(pm_var4d_nc,2)/size(pm_var4d_nc,4)
        write(unit_logfile,'(A,f12.3)')' LF PPMCO  after correction variables are used: ',sum(pm_var4d_nc(:,:,surface_level_nc_2,:,conc_nc_index,allsource_nc_index,1))/size(pm_var4d_nc,1)/size(pm_var4d_nc,2)/size(pm_var4d_nc,4)
    endif
    
    !Transfer all source values to all the sources for use in source looping later
    do i_source=1,n_source_nc_index
        !if (calculate_source(i_source).and.i_source.ne.allsource_nc_index) then
        if (i_source.ne.allsource_nc_index) then
        var3d_nc(:,:,:,conc_nc_index,i_source,:)=var3d_nc(:,:,:,conc_nc_index,allsource_nc_index,:)
        var4d_nc(:,:,:,:,conc_nc_index,i_source,:)=var4d_nc(:,:,:,:,conc_nc_index,allsource_nc_index,:)
        pm_var4d_nc(:,:,:,:,conc_nc_index,i_source,:)=pm_var4d_nc(:,:,:,:,conc_nc_index,allsource_nc_index,:)
        endif
    enddo
        
        !Transfer local contribution 4d to 3d since this is the only one currently used
        !write(*,*) shape(var4d_nc)
        !write(*,*) surface_level_nc
        lc_var3d_nc=lc_var4d_nc(:,:,:,:,surface_level_nc_2,:,:,:,:)
        !var3d_nc(:,:,:,conc_nc_index,:)=var4d_nc(:,:,surface_level_nc,:,conc_nc_index,:)
        !Don't do this because if 3d is read then it is obliterated
        !if (sum(var3d_nc(:,:,:,conc_nc_index,:,:)).eq.0) then
        var3d_nc(:,:,:,conc_nc_index,:,:)=var4d_nc(:,:,surface_level_nc_2,:,conc_nc_index,:,:) !Changed tis from surface_level_nc_2 to surface_level_nc under lf changes
        !endif
        !if (sum(comp_var3d_nc(:,:,:,:)).eq.0) then
        comp_var3d_nc(:,:,:,:)=comp_var4d_nc(:,:,surface_level_nc,:,:)
        !endif
        
        if (allocated(lc_var4d_nc)) deallocate(lc_var4d_nc)
        if (allocated(comp_var4d_nc)) deallocate(comp_var4d_nc)

        !Adjust the depositions to be total depositions and not just, for example, N
        do i_pollutant=1,n_emep_pollutant_loop
        var3d_nc(:,:,:,drydepo_nc_index,:,i_pollutant)=var3d_nc(:,:,:,drydepo_nc_index,:,i_pollutant)*depo_scale_nc(pollutant_loop_index(i_pollutant))
        var3d_nc(:,:,:,wetdepo_nc_index,:,i_pollutant)=var3d_nc(:,:,:,wetdepo_nc_index,:,i_pollutant)*depo_scale_nc(pollutant_loop_index(i_pollutant))
        enddo
        
        !Use average lowest level grid thickness
        H_emep_temp=sum(var4d_nc(:,:,surface_level_nc,:,ZTOP_nc_index,allsource_index,meteo_p_loop_index))/dim_length_nc(x_dim_nc_index)/dim_length_nc(y_dim_nc_index)/dim_length_nc(time_dim_nc_index)
        if (H_emep_temp.ne.0) then
            H_emep=H_emep_temp
            write(unit_logfile,'(A,f8.2)')' Using model depth info. Setting lowest level depth = ',H_emep
        else
            write(unit_logfile,'(A,f8.2)')' No model depth info available. Setting lowest level depth = ',H_emep
        endif
    
        H_meteo=H_emep/2.
        write(unit_logfile,'(A,f8.2)')' Setting lowest meteo grid height = ',H_meteo
        
        !var3d_nc(:,:,:,conc_nc_index,:,:)=var4d_nc(:,:,surface_level_nc,:,conc_nc_index,:,:)
        
        !For nh3 only use the comp value not the value in the local fraction file
        !NOTE: TEMPORARY
        if (use_comp_temporary) then
            write(*,*) 'WARNING: TEMPORARALLY USING COMP FOR LF SOURCE SCALING FOR NH3',pollutant_loop_back_index(nh3_nc_index),nh3_nc_index
            do i_source=1,n_source_index
            if (calculate_source(i_source).or.calculate_EMEP_source(i_source).or.i_source.eq.allsource_index) then
            var3d_nc(:,:,:,conc_nc_index,i_source,pollutant_loop_back_index(nh3_nc_index))=comp_var3d_nc(:,:,:,nh3_nc_index)               
            endif
            enddo
        endif

        !At the moment the local contribution based on fraction. Convert to local contributions here
        !Remove this if we read local contributions in a later version
        do j=1,dim_length_nc(ydist_dim_nc_index)
        do i=1,dim_length_nc(xdist_dim_nc_index)
        
        do p_loop=1,n_pollutant_loop
            !do lc_local_nc_index=minval(lc_local_nc_loop_index),maxval(lc_local_nc_loop_index)
            do ii=1,n_local_fraction_grids
                lc_local_nc_index=lc_local_nc_loop_index(ii)
                lc_frac_nc_index=lc_frac_nc_loop_index(ii)
            !lc_frac_nc_index=convert_local_to_fraction_loop_index(lc_local_nc_index)
            !write(*,*) lc_local_nc_index,lc_frac_nc_index
            if (pollutant_loop_index(p_loop).eq.pm10_nc_index) then
                lc_var3d_nc(i,j,:,:,:,lc_local_nc_index,:,p_loop)=pm_var4d_nc(:,:,surface_level_nc_2,:,conc_nc_index,:,1)*pm_lc_var4d_nc(i,j,:,:,surface_level_nc_2,:,lc_frac_nc_index,:,1) &
                    +pm_var4d_nc(:,:,surface_level_nc_2,:,conc_nc_index,:,2)*pm_lc_var4d_nc(i,j,:,:,surface_level_nc_2,:,lc_frac_nc_index,:,2)
            elseif (pollutant_loop_index(p_loop).eq.pm25_nc_index) then
                lc_var3d_nc(i,j,:,:,:,lc_local_nc_index,:,p_loop)=pm_var4d_nc(:,:,surface_level_nc_2,:,conc_nc_index,:,2)*pm_lc_var4d_nc(i,j,:,:,surface_level_nc_2,:,lc_frac_nc_index,:,2)
            else
                lc_var3d_nc(i,j,:,:,:,lc_local_nc_index,:,p_loop)=var4d_nc(:,:,surface_level_nc_2,:,conc_nc_index,:,p_loop)*lc_var3d_nc(i,j,:,:,:,lc_frac_nc_index,:,p_loop) !This was before, so it used the D3 value, but not the surface value. D3 is not always read
                !lc_var3d_nc(i,j,:,:,:,lc_local_nc_index,:,p_loop)=max(lc_var3d_nc(i,j,:,:,:,lc_local_nc_index,:,p_loop),var3d_nc(:,:,:,conc_nc_index,:,p_loop)*lc_var3d_nc(i,j,:,:,:,lc_frac_nc_index,:,p_loop)) !Choose the max of the 3d and 4d values, so chooses one or the other
                !lc_var3d_nc(i,j,:,:,:,lc_local_nc_index,:,p_loop)=pm_var4d_nc(:,:,surface_level_nc_2,:,conc_nc_index,:,2)*lc_var3d_nc(i,j,:,:,:,lc_frac_nc_index,:,p_loop)
                !write(*,*) sum(lc_var3d_nc(i,j,:,:,:,lc_local_nc_index,:,p_loop)),sum(var3d_nc(:,:,:,conc_nc_index,:,p_loop)),sum(lc_var3d_nc(i,j,:,:,:,lc_frac_nc_index,:,p_loop))          
            endif
            enddo
        enddo
                 
        enddo
        enddo

    !Take account of the fact that the GNFR emissions can be version 13 or 19
    !In which case traffic is split into 4 and power is split into 2
    !This is valid only for the local fraction sources
    if (use_GNFR19_emissions_from_EMEP_flag) then
        write(unit_logfile,'(3A,f16.4)') ' Aggregating GNFR19 to GNFR13 or GNFR14'
        
        !Reset these flags so they are no longer calculated as EMEP contributions after being read
        !calculate_EMEP_source(traffic_exhaust_nc_index)=.false.
        !calculate_EMEP_source(traffic_nonexhaust_nc_index)=.false.

        do lc_local_nc_index=minval(lc_local_nc_loop_index),maxval(lc_local_nc_loop_index)
            !No allocation of exhaust emissions to PM10 because only course is read
            !Is this correct? Removing. PM10 should have been set earlier
            lc_var3d_nc(:,:,:,:,:,lc_local_nc_index,traffic_gasoline_nc_index,pollutant_loop_back_index(pm10_nc_index))=lc_var3d_nc(:,:,:,:,:,lc_local_nc_index,traffic_gasoline_nc_index,pollutant_loop_back_index(pm25_nc_index))
            lc_var3d_nc(:,:,:,:,:,lc_local_nc_index,traffic_diesel_nc_index,pollutant_loop_back_index(pm10_nc_index))=lc_var3d_nc(:,:,:,:,:,lc_local_nc_index,traffic_diesel_nc_index,pollutant_loop_back_index(pm25_nc_index))
            lc_var3d_nc(:,:,:,:,:,lc_local_nc_index,traffic_gas_nc_index,pollutant_loop_back_index(pm10_nc_index))=lc_var3d_nc(:,:,:,:,:,lc_local_nc_index,traffic_gas_nc_index,pollutant_loop_back_index(pm25_nc_index))
            
            !Put all source into traffic
            lc_var3d_nc(:,:,:,:,:,lc_local_nc_index,traffic_nc_index,:)=lc_var3d_nc(:,:,:,:,:,lc_local_nc_index,traffic_nc_index,:) &
                                                                        +lc_var3d_nc(:,:,:,:,:,lc_local_nc_index,traffic_gasoline_nc_index,:) &
                                                                        +lc_var3d_nc(:,:,:,:,:,lc_local_nc_index,traffic_diesel_nc_index,:) &
                                                                        +lc_var3d_nc(:,:,:,:,:,lc_local_nc_index,traffic_gas_nc_index,:) &
                                                                        +lc_var3d_nc(:,:,:,:,:,lc_local_nc_index,traffic_nonexhaust_nc_index,:)
            !Put the three exhaust into traffic. If exhaust not included then it will be 0
            !lc_var3d_nc(:,:,:,:,:,lc_local_nc_index,traffic_exhaust_nc_index,pollutant_loop_back_index(pm25_nc_index))=lc_var3d_nc(:,:,:,:,:,lc_local_nc_index,traffic_gasoline_nc_index,pollutant_loop_back_index(pm25_nc_index)) &
            !                                                            +lc_var3d_nc(:,:,:,:,:,lc_local_nc_index,traffic_diesel_nc_index,pollutant_loop_back_index(pm25_nc_index)) &
            !                                                            +lc_var3d_nc(:,:,:,:,:,lc_local_nc_index,traffic_gas_nc_index,pollutant_loop_back_index(pm25_nc_index))
            !lc_var3d_nc(:,:,:,:,:,lc_local_nc_index,traffic_exhaust_nc_index,pollutant_loop_back_index(pm10_nc_index))=lc_var3d_nc(:,:,:,:,:,lc_local_nc_index,traffic_exhaust_nc_index,pollutant_loop_back_index(pm25_nc_index))
            lc_var3d_nc(:,:,:,:,:,lc_local_nc_index,traffic_exhaust_nc_index,:)=lc_var3d_nc(:,:,:,:,:,lc_local_nc_index,traffic_gasoline_nc_index,:) &
                                                                        +lc_var3d_nc(:,:,:,:,:,lc_local_nc_index,traffic_diesel_nc_index,:) &
                                                                        +lc_var3d_nc(:,:,:,:,:,lc_local_nc_index,traffic_gas_nc_index,:)
            !Special case for PM10 because traffic exhaust is not read for PM10
            !lc_var3d_nc(:,:,:,:,:,lc_local_nc_index,traffic_exhaust_nc_index,pollutant_loop_back_index(pm10_nc_index))=lc_var3d_nc(:,:,:,:,:,lc_local_nc_index,traffic_exhaust_nc_index,pollutant_loop_back_index(pm25_nc_index))
            !Aggregate the two public powers
            lc_var3d_nc(:,:,:,:,:,lc_local_nc_index,publicpower_nc_index,:)=lc_var3d_nc(:,:,:,:,:,lc_local_nc_index,publicpower_nc_index,:) &
                                                                        +lc_var3d_nc(:,:,:,:,:,lc_local_nc_index,publicpower_point_nc_index,:) &
                                                                        +lc_var3d_nc(:,:,:,:,:,lc_local_nc_index,publicpower_area_nc_index,:)
        enddo
        do lc_local_nc_index=minval(lc_local_nc_loop_index),maxval(lc_local_nc_loop_index)
            write(unit_logfile,'(A,i,a,f16.4)') 'Mean exhaust PM2.5 EMEP contribution centre lc grid: ',lc_local_nc_index,' ' &
                ,sum(lc_var3d_nc(xdist_centre_nc,ydist_centre_nc,:,:,:,lc_local_nc_index,traffic_exhaust_nc_index,pollutant_loop_back_index(pm25_nc_index)))/(size(lc_var3d_nc,3)*size(lc_var3d_nc,4)*size(lc_var3d_nc,5))
            write(unit_logfile,'(A,i,a,f16.4)') 'Mean nonexhaust PM2.5 EMEP contribution centre lc grid: ',lc_local_nc_index,' ' &
                ,sum(lc_var3d_nc(xdist_centre_nc,ydist_centre_nc,:,:,:,lc_local_nc_index,traffic_nonexhaust_nc_index,pollutant_loop_back_index(pm25_nc_index)))/(size(lc_var3d_nc,3)*size(lc_var3d_nc,4)*size(lc_var3d_nc,5))
        enddo
    !else
    !    do lc_local_nc_index=minval(lc_local_nc_loop_index),maxval(lc_local_nc_loop_index)
    !    lc_var3d_nc(:,:,:,:,:,lc_local_nc_index,traffic_nc_index,pollutant_loop_back_index(pmex_nc_index))=0
    !    enddo
    endif
    
    !With GNFR19 add up all the traffic emissions and put them in the traffic emissions, according to the use of all or all_totals
    if (use_GNFR19_emissions_from_EMEP_flag) then
        if (pollutant_index.eq.all_totals_nc_index) write(unit_logfile,'(A)') 'Aggregating exhaust and non-exhaust traffic emissions when using GNFR19: '
        if (pollutant_index.eq.all_nc_index) write(unit_logfile,'(A)') 'Aggregating exhaust and non-exhaust traffic emissions seperately when using GNFR19: '
        
        !Aggregate exhaust emissions
        var3d_nc(:,:,:,emis_nc_index,traffic_exhaust_nc_index,:)=var3d_nc(:,:,:,emis_nc_index,traffic_gasoline_nc_index,:) &
                                                                +var3d_nc(:,:,:,emis_nc_index,traffic_diesel_nc_index,:) &
                                                                +var3d_nc(:,:,:,emis_nc_index,traffic_gas_nc_index,:)
        !Aggregating total emissions
        var3d_nc(:,:,:,emis_nc_index,traffic_nc_index,:)=var3d_nc(:,:,:,emis_nc_index,traffic_nc_index,:) &
                                                        +var3d_nc(:,:,:,emis_nc_index,traffic_exhaust_nc_index,:) &
                                                        +var3d_nc(:,:,:,emis_nc_index,traffic_nonexhaust_nc_index,:)
        !Aggregating public power emissions
        var3d_nc(:,:,:,emis_nc_index,publicpower_nc_index,:)=var3d_nc(:,:,:,emis_nc_index,publicpower_nc_index,:) &
                                                        +var3d_nc(:,:,:,emis_nc_index,publicpower_point_nc_index,:) &
                                                        +var3d_nc(:,:,:,emis_nc_index,publicpower_area_nc_index,:)
    endif
    
    !Scale the read in voc emissions so they are split into benzene
    if (extract_benzene_from_voc_emissions) then
        write(unit_logfile,'(A)') 'Converting VOC emissions to Benzene emissions: '
        do i_source=1,n_source_index
            if (calculate_source(i_source).or.calculate_EMEP_source(i_source).or.save_EMEP_source(i_source)) then
            var3d_nc(:,:,:,emis_nc_index,i_source,pollutant_loop_back_index(c6h6_nc_index))=benzene_split_voc_in_GNFR_sectors(uEMEP_to_EMEP_sector(i_source))*var3d_nc(:,:,:,emis_nc_index,i_source,pollutant_loop_back_index(c6h6_nc_index))
            write(unit_logfile,'(A,i4,a,a,a,f12.3)') 'GNFR Source=',uEMEP_to_EMEP_sector(i_source),' , uEMEP source= ',trim(source_file_str(i_source)), ' , split(%)= ',benzene_split_voc_in_GNFR_sectors(uEMEP_to_EMEP_sector(i_source))*100.
            endif
        enddo       
    endif
    
        !Check output for PM10. Not needed anymore
        !do i_source=1,n_source_index
        !    if (calculate_source(i_source).or.calculate_EMEP_source(i_source).or.save_EMEP_source(i_source)) then
        !    do lc_local_nc_index=minval(lc_local_nc_loop_index),maxval(lc_local_nc_loop_index)
        !        p_loop=pollutant_loop_back_index(pm10_nc_index)
        !        write(unit_logfile,'(3A,f16.4)') ' Average local contribution of: ',trim(var_name_nc(conc_nc_index,pm10_nc_index,allsource_index)),' '//trim(source_file_str(i_source)), &
        !        sum(lc_var3d_nc(xdist_centre_nc,ydist_centre_nc,:,:,:,lc_local_nc_index,i_source,p_loop))/(size(lc_var3d_nc,3)*size(lc_var3d_nc,4)*size(lc_var3d_nc,5))
        !    enddo
        !    endif
        !enddo

        if (allocated(pm_lc_var4d_nc)) deallocate(pm_lc_var4d_nc)
        if (allocated(pm_var4d_nc)) deallocate(pm_var4d_nc)
        if (allocated(pm_var3d_nc)) deallocate(pm_var3d_nc)
        
       !Set the local grid contribution for the individual grid
        !write(*,*) shape(lc_var4d_nc)
        !write(*,*) shape(var4d_nc)
        !Commented out as these are not used?

        
        !Reset these variables to what they would be with 1 EMEP grid
        !Not necessary once the rest of the code is adapted to the loop data
        frac_nc_index=num_var_nc_start+1
        lc_frac_nc_index=1
        local_nc_index=num_var_nc_start+n_local_fraction_grids+1
        lc_local_nc_index=n_local_fraction_grids+1
        
        !This is no longer in use, but just in case set the variables to the first EMEP LC grid
        var3d_nc(:,:,:,frac_nc_index,:,:)=lc_var3d_nc(xdist_centre_nc,ydist_centre_nc,:,:,:,lc_frac_nc_index,:,:)
        var3d_nc(:,:,:,local_nc_index,:,:)=var3d_nc(:,:,:,conc_nc_index,:,:)*var3d_nc(:,:,:,frac_nc_index,:,:)
                
        !write(*,*) minval(var4d_nc(:,:,:,:,local_nc_index,:)),maxval(var4d_nc(:,:,:,:,local_nc_index,:))
        !write(*,*) minval(var4d_nc(:,:,:,:,frac_nc_index,:)),maxval(var4d_nc(:,:,:,:,frac_nc_index,:))
        !write(*,*) minval(var3d_nc(:,:,:,frac_nc_index,:)),maxval(var3d_nc(:,:,:,frac_nc_index,:))
        !write(*,*) minval(var4d_nc(:,:,:,:,conc_nc_index,:)),maxval(var4d_nc(:,:,:,:,conc_nc_index,:))
        !write(*,*) minval(var3d_nc(:,:,:,conc_nc_index,:)),maxval(var3d_nc(:,:,:,conc_nc_index,:))
        !write(*,*) minval(var3d_nc(:,:,:,inv_FF10_nc_index,allsource_index)),maxval(var3d_nc(:,:,:,inv_FF10_nc_index,allsource_index))
        !write(*,*) minval(var3d_nc(:,:,:,FF10_nc_index,allsource_index)),maxval(var3d_nc(:,:,:,FF10_nc_index,allsource_index))
        
        
        where (var3d_nc(:,:,:,ustar_nc_index,:,meteo_p_loop_index).lt.ustar_min) var3d_nc(:,:,:,ustar_nc_index,:,meteo_p_loop_index)=ustar_min

        !Test for 0 wind speed components if valid
        if (.not.EMEP_region_outside_domain) then
        do j=1,dim_length_nc(y_dim_nc_index)
        do i=1,dim_length_nc(x_dim_nc_index)
        do t=1,dim_length_nc(time_dim_nc_index)
            if (var4d_nc(i,j,surface_level_nc,t,ugrid_nc_index,allsource_nc_index,meteo_p_loop_index).eq.0..and.var4d_nc(i,j,surface_level_nc,t,vgrid_nc_index,allsource_nc_index,meteo_p_loop_index).eq.0.) then
                write(unit_logfile,'(a,3i)') 'Zero wind fields at (i,j,t): ',i,j,t
                !Search for the nearest non double zero value
                k=0
                nonzero_wind_notfound=.true.
                do while (k.lt.20.and.nonzero_wind_notfound)
                k=k+1
                    do ii=-k,k
                    do jj =-k,k
                        iii=i+ii
                        jjj=j+jj
                        iii=min(max(1,iii),dim_length_nc(x_dim_nc_index))
                        jjj=min(max(1,jjj),dim_length_nc(y_dim_nc_index))
                        if (var4d_nc(iii,jjj,surface_level_nc,t,ugrid_nc_index,allsource_nc_index,meteo_p_loop_index).ne.0. &
                            .or.var4d_nc(iii,jjj,surface_level_nc,t,vgrid_nc_index,allsource_nc_index,meteo_p_loop_index).ne.0.) then
                            nonzero_wind_notfound=.false.
                            var4d_nc(i,j,surface_level_nc,t,ugrid_nc_index,allsource_nc_index,meteo_p_loop_index)=var4d_nc(iii,jjj,surface_level_nc,t,ugrid_nc_index,allsource_nc_index,meteo_p_loop_index)
                            var4d_nc(i,j,surface_level_nc,t,vgrid_nc_index,allsource_nc_index,meteo_p_loop_index)=var4d_nc(iii,jjj,surface_level_nc,t,vgrid_nc_index,allsource_nc_index,meteo_p_loop_index)
                            !write(unit_logfile,'(a,4i,2f10.2)') 'Wind found for (i,j,t): ',k,i,j,t,var4d_nc(i,j,surface_level_nc,t,ugrid_nc_index,allsource_nc_index,meteo_p_loop_index),var4d_nc(i,j,surface_level_nc,t,vgrid_nc_index,allsource_nc_index,meteo_p_loop_index)
                        endif
                    enddo
                    enddo
                        
                enddo
                if (nonzero_wind_notfound) then
                    write(unit_logfile,'(a,4i,2f10.2)') 'Wind not found for (loop,i,j,t): ',k,i,j,t,var4d_nc(i,j,surface_level_nc,t,ugrid_nc_index,allsource_nc_index,meteo_p_loop_index),var4d_nc(i,j,surface_level_nc,t,vgrid_nc_index,allsource_nc_index,meteo_p_loop_index)
                else
                    write(unit_logfile,'(a,4i,2f10.2)') 'Wind found for (loop,i,j,t): ',k,i,j,t,var4d_nc(i,j,surface_level_nc,t,ugrid_nc_index,allsource_nc_index,meteo_p_loop_index),var4d_nc(i,j,surface_level_nc,t,vgrid_nc_index,allsource_nc_index,meteo_p_loop_index)
                endif
            endif
        enddo
        enddo
        enddo

        do j=1,dim_length_nc(y_dim_nc_index)
        do i=1,dim_length_nc(x_dim_nc_index)
        do t=1,dim_length_nc(time_dim_nc_index)
            if (var4d_nc(i,j,surface_level_nc,t,ugrid_nc_index,allsource_nc_index,meteo_p_loop_index).eq.0..and.var4d_nc(i,j,surface_level_nc,t,vgrid_nc_index,allsource_nc_index,meteo_p_loop_index).eq.0.) then
                write(unit_logfile,'(a,3i)') 'ERROR: Found zero wind fields in both components (i,j,t). Stopping: ',i,j,t
                stop
            endif
        enddo
        enddo
        enddo
            
        endif
        
        !If no logz0 available. Set to log(0.1)
        !For urban areas a value of 0.3 is used
        where (var3d_nc(:,:,:,logz0_nc_index,:,meteo_p_loop_index).eq.0.0) var3d_nc(:,:,:,logz0_nc_index,:,meteo_p_loop_index)=log(0.3)
        if (replace_z0.ne.NODATA_value) then
            write(unit_logfile,'(A,f8.4)') ' Replacing z0 everywhere with: ',replace_z0
            var3d_nc(:,:,:,logz0_nc_index,:,meteo_p_loop_index)=log(replace_z0)
        endif
        if (replace_invL.ne.NODATA_value) then
            write(unit_logfile,'(A,f8.4)') ' Replacing inverse L everywhere with: ',replace_invL
            var3d_nc(:,:,:,invL_nc_index,:,meteo_p_loop_index)=replace_invL
        endif
        if (replace_hmix.ne.NODATA_value) then
            write(unit_logfile,'(A,f8.4)') ' Replacing HMIX everywhere with: ',replace_hmix
            var3d_nc(:,:,:,hmix_nc_index,:,meteo_p_loop_index)=replace_hmix
        endif
        
        if (use_phi_for_invL) then
            phi_count=0
            mean_phi_temp=0
            mean_invL_temp=0
            do j=1,dim_length_nc(y_dim_nc_index)
            do i=1,dim_length_nc(x_dim_nc_index)
            do t=1,dim_length_nc(time_dim_nc_index)
                call TROENKz_invL_from_phi(z_invL,var3d_nc(i,j,t,phi_nc_index,allsource_nc_index,meteo_p_loop_index),var3d_nc(i,j,t,invL_nc_index,allsource_nc_index,meteo_p_loop_index))
                phi_count=phi_count+1
                mean_phi_temp=mean_phi_temp+var3d_nc(i,j,t,phi_nc_index,allsource_nc_index,meteo_p_loop_index)
                mean_invL_temp=mean_invL_temp+var3d_nc(i,j,t,invL_nc_index,allsource_nc_index,meteo_p_loop_index)

            enddo
            enddo
            enddo
            
            mean_phi_temp=mean_phi_temp/phi_count
            mean_invL_temp=mean_invL_temp/phi_count

            write(unit_logfile,'(A,2f8.4)') ' Using phi instead of invL. Mean phi and invL: ',mean_phi_temp,mean_invL_temp
            
        endif
        
        !Limit stable L to lowest_stable_L and to lowest_unstable_L (negative number) for unstable.
        where (var3d_nc(:,:,:,invL_nc_index,:,meteo_p_loop_index).lt.1.0/lowest_unstable_L) var3d_nc(:,:,:,invL_nc_index,:,meteo_p_loop_index)=1.0/lowest_unstable_L
        where (var3d_nc(:,:,:,invL_nc_index,:,meteo_p_loop_index).gt.1.0/lowest_stable_L) var3d_nc(:,:,:,invL_nc_index,:,meteo_p_loop_index)=1.0/lowest_stable_L

        where (var3d_nc(:,:,:,hmix_nc_index,:,meteo_p_loop_index).lt.hmix_min) var3d_nc(:,:,:,hmix_nc_index,:,meteo_p_loop_index)=hmix_min

        !Limit Jd as well in case there is something wrong
        where (var4d_nc(:,:,:,:,J_nc_index,:,:).lt.0) var4d_nc(:,:,:,:,J_nc_index,:,:)=0.
        
        if (J_scale.ne.NODATA_value) then
            write(unit_logfile,'(A,f8.4)') ' Scaling J(NO2) everywhere with: ',J_scale
            var4d_nc(:,:,:,:,J_nc_index,:,:)=var4d_nc(:,:,:,:,J_nc_index,:,:)*J_scale
        endif

        !Correct the inverse of the wind speed for the factor 0.2 used to create it in EMEP
        write(unit_logfile,'(A)') ' Correcting inverse wind speed to account for the 0.2 m/s offset: '
        var3d_nc(:,:,:,inv_FF10_nc_index,:,:)=var3d_nc(:,:,:,inv_FF10_nc_index,:,:)/(1.-0.2*var3d_nc(:,:,:,inv_FF10_nc_index,:,:))
        var3d_nc(:,:,:,inv_FFgrid_nc_index,:,:)=var3d_nc(:,:,:,inv_FFgrid_nc_index,:,:)/(1.-0.2*var3d_nc(:,:,:,inv_FFgrid_nc_index,:,:))
        var4d_nc(:,:,:,:,inv_FFgrid_nc_index,:,:)=var4d_nc(:,:,:,:,inv_FFgrid_nc_index,:,:)/(1.-0.2*var4d_nc(:,:,:,:,inv_FFgrid_nc_index,:,:))
        where (var3d_nc(:,:,:,inv_FF10_nc_index,:,:).gt.4.5) var3d_nc(:,:,:,inv_FF10_nc_index,:,:)=4.5 !Set the limit so that FF can not be less than 0.2222
        where (var3d_nc(:,:,:,inv_FFgrid_nc_index,:,:).gt.4.5) var3d_nc(:,:,:,inv_FFgrid_nc_index,:,:)=4.5 !Set the limit so that FF can not be less than 0.2222
        where (var4d_nc(:,:,:,:,inv_FFgrid_nc_index,:,:).gt.4.5) var4d_nc(:,:,:,:,inv_FFgrid_nc_index,:,:)=4.5 !Set the limit so that FF can not be less than 0.2222
        
        if (FF_scale.ne.NODATA_value) then
            write(unit_logfile,'(A,f8.4)') ' Rescaling wind fields everywhere with factor: ',FF_scale
            var3d_nc(:,:,:,ustar_nc_index,:,meteo_p_loop_index)=var3d_nc(:,:,:,ustar_nc_index,:,meteo_p_loop_index)*FF_scale
            var3d_nc(:,:,:,FF10_nc_index,:,meteo_p_loop_index)=var3d_nc(:,:,:,FF10_nc_index,:,meteo_p_loop_index)*FF_scale
            var3d_nc(:,:,:,inv_FF10_nc_index,:,meteo_p_loop_index)=var3d_nc(:,:,:,inv_FF10_nc_index,:,meteo_p_loop_index)/FF_scale
            var3d_nc(:,:,:,ugrid_nc_index,:,meteo_p_loop_index)=var3d_nc(:,:,:,ugrid_nc_index,:,meteo_p_loop_index)*FF_scale
            var3d_nc(:,:,:,vgrid_nc_index,:,meteo_p_loop_index)=var3d_nc(:,:,:,vgrid_nc_index,:,meteo_p_loop_index)*FF_scale
            var3d_nc(:,:,:,u10_nc_index,:,meteo_p_loop_index)=var3d_nc(:,:,:,u10_nc_index,:,meteo_p_loop_index)*FF_scale
            var3d_nc(:,:,:,v10_nc_index,:,meteo_p_loop_index)=var3d_nc(:,:,:,v10_nc_index,:,meteo_p_loop_index)*FF_scale
            var3d_nc(:,:,:,inv_FFgrid_nc_index,:,meteo_p_loop_index)=var3d_nc(:,:,:,inv_FFgrid_nc_index,:,meteo_p_loop_index)/FF_scale
            var4d_nc(:,:,:,:,ugrid_nc_index,:,meteo_p_loop_index)=var4d_nc(:,:,:,:,ugrid_nc_index,:,meteo_p_loop_index)*FF_scale
            var4d_nc(:,:,:,:,vgrid_nc_index,:,meteo_p_loop_index)=var4d_nc(:,:,:,:,vgrid_nc_index,:,meteo_p_loop_index)*FF_scale
            var4d_nc(:,:,:,:,inv_FFgrid_nc_index,:,meteo_p_loop_index)=var4d_nc(:,:,:,:,inv_FFgrid_nc_index,:,meteo_p_loop_index)/FF_scale
        endif
        if (FF10_offset.ne.NODATA_value) then
            write(unit_logfile,'(A,f8.4)') ' Offsetting 10 m wind fields everywhere with a value: ',FF10_offset
            var3d_nc(:,:,:,FF10_nc_index,:,meteo_p_loop_index)=var3d_nc(:,:,:,FF10_nc_index,:,meteo_p_loop_index)+FF10_offset
         endif
        if (DD_offset.ne.NODATA_value) then
            write(unit_logfile,'(A,f8.4)') ' Rotating wind fields everywhere with a value: ',DD_offset

            !Make use of the spare source index parts of the array for the conversion
            temp_var4d_nc=0
            temp_var4d_nc(:,:,:,:,1) = var4d_nc(:,:,:,:,ugrid_nc_index,allsource_index,meteo_p_loop_index)*cos(DD_offset/180.*3.14159)+var4d_nc(:,:,:,:,vgrid_nc_index,allsource_index,meteo_p_loop_index)*sin(DD_offset/180.*3.14159)                                       
            temp_var4d_nc(:,:,:,:,2) =-var4d_nc(:,:,:,:,ugrid_nc_index,allsource_index,meteo_p_loop_index)*sin(DD_offset/180.*3.14159)+var4d_nc(:,:,:,:,vgrid_nc_index,allsource_index,meteo_p_loop_index)*cos(DD_offset/180.*3.14159)                                       
            var4d_nc(:,:,:,:,ugrid_nc_index,allsource_nc_index,meteo_p_loop_index) = temp_var4d_nc(:,:,:,:,1)
            var4d_nc(:,:,:,:,vgrid_nc_index,allsource_nc_index,meteo_p_loop_index) = temp_var4d_nc(:,:,:,:,2)
            temp_var4d_nc=0
            temp_var4d_nc(:,:,:,1,1) = var3d_nc(:,:,:,u10_nc_index,allsource_index,meteo_p_loop_index)*cos(DD_offset/180.*3.14159)+var3d_nc(:,:,:,vgrid_nc_index,allsource_index,meteo_p_loop_index)*sin(DD_offset/180.*3.14159)                                       
            temp_var4d_nc(:,:,:,1,2) =-var3d_nc(:,:,:,u10_nc_index,allsource_index,meteo_p_loop_index)*sin(DD_offset/180.*3.14159)+var3d_nc(:,:,:,vgrid_nc_index,allsource_index,meteo_p_loop_index)*cos(DD_offset/180.*3.14159)                                       
            var3d_nc(:,:,:,u10_nc_index,allsource_nc_index,meteo_p_loop_index) = temp_var4d_nc(:,:,:,1,1)
            var3d_nc(:,:,:,v10_nc_index,allsource_nc_index,meteo_p_loop_index) = temp_var4d_nc(:,:,:,1,2)
        endif
       
        !Set the magnitude of the gridded wind fields. Should probably be done after subgridding?
        var4d_nc(:,:,:,:,FFgrid_nc_index,allsource_index,meteo_p_loop_index)=sqrt(var4d_nc(:,:,:,:,ugrid_nc_index,allsource_index,meteo_p_loop_index)**2+var4d_nc(:,:,:,:,vgrid_nc_index,allsource_index,meteo_p_loop_index)**2)
        !Will override the read in FF10 if it is available
        if (sum(abs(var3d_nc(:,:,:,u10_nc_index,allsource_index,meteo_p_loop_index))).ne.0.and.sum(abs(var3d_nc(:,:,:,v10_nc_index,allsource_index,meteo_p_loop_index))).ne.0) then
            wind_vectors_10m_available=.true.
            var3d_nc(:,:,:,FF10_nc_index,allsource_index,meteo_p_loop_index)=sqrt(var3d_nc(:,:,:,u10_nc_index,allsource_index,meteo_p_loop_index)**2+var3d_nc(:,:,:,v10_nc_index,allsource_index,meteo_p_loop_index)**2)
        endif
        write(unit_logfile,'(A,L)') ' 10 m wind vectors available: ',wind_vectors_10m_available           
            
  
        !Check EMEP time
        date_num_temp=dble(ceiling(val_dim_nc(1,time_dim_nc_index)*24.))/24.
        call number_to_date(date_num_temp,date_array,ref_year_EMEP)
        write(unit_logfile,'(a,i6)') ' Time dimension EMEP:  ',dim_length_nc(time_dim_nc_index)
        write(unit_logfile,'(a,6i6)') ' Date start EMEP =  ',date_array
        date_num_temp=dble(ceiling(val_dim_nc(dim_length_nc(time_dim_nc_index),time_dim_nc_index)*24.))/24.
        call number_to_date(date_num_temp,date_array,ref_year_EMEP)
        write(unit_logfile,'(a,6i6)') ' Date end EMEP =    ',date_array
        
        !Test and correct dates
        if (1.eq.1) then
        if (.not.allocated(time_seconds_output)) allocate(time_seconds_output(dim_length_nc(time_dim_nc_index)))
        do t=1,dim_length_nc(time_dim_nc_index)
            date_num_temp=dble(ceiling(val_dim_nc(t,time_dim_nc_index)*24.))/24.
            date_num_temp=val_dim_nc(t,time_dim_nc_index)+0.55/24. !Add a bit over half an hour to compensate for average of time
            call number_to_date(date_num_temp,date_array,ref_year_EMEP)
            !write(unit_logfile,'(a,i4,6i6,d)') ' Date EMEP =   ',t,date_array,date_num_temp
            
            date_array(5:6)=0 !Set minutes and hours to 0
            date_num_temp=date_to_number(date_array,ref_year_EMEP)
            call number_to_date(date_num_temp,date_array,ref_year_EMEP)
            !write(unit_logfile,'(a,i4,6i6,d)') ' Date EMEP =   ',t,date_array,date_num_temp
            !val_dim_nc(t,time_dim_nc_index)=ceiling(date_num_temp*dble(24.)*dble(3600.))/dble(24.)/dble(3600.)
            date_num_temp=date_num_temp+dble(0.01)/dble(24.)/dble(3600.) !Add 0.01 of a second to avoid any rounding off errors
            call number_to_date(date_num_temp,date_array,ref_year_EMEP)
            !write(unit_logfile,'(a,i4,6i6,d)') ' Date EMEP =   ',t,date_array,date_num_temp
            val_dim_nc(t,time_dim_nc_index)=date_num_temp
            
            !Convert to seconds since 2000
            date_array=0
            date_array(1)=2000;date_array(2)=1;date_array(3)=1
            date_num_2000=date_to_number(date_array,ref_year_EMEP)
            time_seconds_output(t)=int4((date_num_temp-date_num_2000)*24*3600)
            unit_dim_nc(time_dim_nc_index)="seconds since 2000-1-1 0:0:0";;

        enddo
        !stop
        endif
        
        !Test for Ceclius or Kelvin
        if (maxval(var3d_nc(:,:,:,t2m_nc_index,:,meteo_p_loop_index)).lt.150) then
            write(unit_logfile,'(a,f12.2)')'WARNING: Temperature appears to be in Celcius. Converting to Kelvin. Max temperature is: ',maxval(var3d_nc(:,:,:,T2m_nc_index,:,meteo_p_loop_index))
            var3d_nc(:,:,:,t2m_nc_index,:,meteo_p_loop_index)=var3d_nc(:,:,:,T2m_nc_index,:,meteo_p_loop_index)+273.13
        endif
    
        !Deallocate temporary arrays
        if (allocated(var1d_nc_dp)) deallocate (var1d_nc_dp)
        if (allocated(var2d_nc_dp)) deallocate (var2d_nc_dp)
        if (allocated(var3d_nc_dp)) deallocate (var3d_nc_dp)
        if (allocated(var4d_nc_dp)) deallocate (var4d_nc_dp)
        if (allocated(temp_var4d_nc)) deallocate (temp_var4d_nc)
        if (allocated(temp_var3d_nc)) deallocate (temp_var3d_nc)
        if (allocated(species_temp_var3d_nc)) deallocate (species_temp_var3d_nc)
 
        !Shift the EMEP grid to the west by 0.1 degrees. Portugal test
        !var1d_nc(:,x_dim_nc_index)=var1d_nc(:,x_dim_nc_index)-0.1
        

    end subroutine uEMEP_read_EMEP
    
    