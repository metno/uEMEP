!uEMEP_read_EMEP
    
    subroutine uEMEP_read_EMEP
    
    use uEMEP_definitions
    use netcdf
    
    implicit none
    
    integer i,j,k,t
    integer ii,jj
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
    integer xdist_centre_nc,ydist_centre_nc
    integer temp_start_time_nc_index,temp_end_time_nc_index
    integer i_loop
    integer valid_dim_length_nc(num_dims_nc) !dimmensions of file 1
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
    
    integer n_file
    double precision date_num_temp
    integer date_array(6)
    double precision scale_factor_nc

    !Temporary reading variables
    double precision, allocatable :: var1d_nc_dp(:)
    double precision, allocatable :: var2d_nc_dp(:,:)
    double precision, allocatable :: var3d_nc_dp(:,:,:)
    double precision, allocatable :: var4d_nc_dp(:,:,:,:)
    real, allocatable :: swop_var4d_nc(:,:,:,:,:,:)
    real, allocatable :: swop_comp_var4d_nc(:,:,:,:,:)
    
    !Temporary files for roatating wind field and PM
    real, allocatable :: temp_var4d_nc(:,:,:,:,:)
    
    !Temporary PM arrays for reading in PM10
    real, allocatable :: pm_var4d_nc(:,:,:,:,:,:,:)
    real, allocatable :: pm_lc_var4d_nc(:,:,:,:,:,:,:,:,:)
    
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
    surface_level_nc=1
    surface_level_nc_2=1
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
        else
            EMEP_projection_type=LL_projection_index
        endif
        
        !Find the (x,y,z,time,xdist,ydist) dimmensions of the file
        do i_dim=1,num_dims_nc
            status_nc = NF90_INQ_DIMID (id_nc,dim_name_nc(i_dim),dim_id_nc(i_dim))
            status_nc = NF90_INQUIRE_DIMENSION (id_nc,dim_id_nc(i_dim),dimname_temp,dim_length_nc(i_dim))
            if (status_nc .NE. NF90_NOERR) then
                write(unit_logfile,'(A,A,A,I)') 'No dimension information available for ',trim(dim_name_nc(i_dim)),' Setting to 1 with status: ',status_nc
                dim_length_nc(i_dim)=1
            endif
        enddo
        
        
        
        if (subgrid_dim(t_dim_index).gt.dim_length_nc(time_dim_nc_index)) then
            write(unit_logfile,'(A,2I)') 'ERROR: Specified time dimensions are greater than EMEP netcdf dimmensions. Stopping ',subgrid_dim(t_dim_index),dim_length_nc(time_dim_nc_index)
            stop
        endif 

               
        write(unit_logfile,'(A,6I)') ' Size of dimensions (x,y,z,t,xdist,ydist): ',dim_length_nc
        

        dim_start_nc(time_dim_nc_index)=temp_start_time_nc_index
        dim_length_nc(time_dim_nc_index)=min(dim_length_nc(time_dim_nc_index),subgrid_dim(t_dim_index))

        
        write(unit_logfile,'(A,6I)') ' New size of dimensions (x,y,z,t,xdist,ydist): ',dim_length_nc

        if (mod(dim_length_nc(xdist_dim_nc_index),2).ne.1.or.mod(dim_length_nc(ydist_dim_nc_index),2).ne.1) then
            write(unit_logfile,'(A,2I)') ' ERROR: Even sized dimmensions for local contribution. Must be odd: ',dim_length_nc(xdist_dim_nc_index),dim_length_nc(ydist_dim_nc_index)
            stop
        endif
            
        if (i_file.eq.2) then
            xdist_centre_nc=1+dim_length_nc(xdist_dim_nc_index)/2
            ydist_centre_nc=1+dim_length_nc(ydist_dim_nc_index)/2
            write(unit_logfile,'(A,2I)') ' Centre index of local contribution dimmensions: ',xdist_centre_nc,ydist_centre_nc
        endif
        
        !Set the last vertical dimension value as the surface layer index
        if (i_file.eq.1.and.invert_levels_flag) then
            surface_level_nc=dim_length_nc(z_dim_nc_index)
            write(unit_logfile,'(A,I)') ' Surface level set to number of vertical layers: ',surface_level_nc
        endif

        !Calculate the necessary extent of the EMEP grid region and only read these grids
        if (reduce_EMEP_region_flag) then
            !Determine the LL cordinates of the target grid
            !if (EMEP_projection_type.eq.LCC_projection_index) then
                !Retrieve the four corners of the target grid in lat and lon
                call UTM2LL(utm_zone,init_subgrid_min(y_dim_index),init_subgrid_min(x_dim_index),temp_lat(1),temp_lon(1))
                call UTM2LL(utm_zone,init_subgrid_max(y_dim_index),init_subgrid_max(x_dim_index),temp_lat(2),temp_lon(2))
                call UTM2LL(utm_zone,init_subgrid_max(y_dim_index),init_subgrid_min(x_dim_index),temp_lat(3),temp_lon(3))
                call UTM2LL(utm_zone,init_subgrid_min(y_dim_index),init_subgrid_max(x_dim_index),temp_lat(4),temp_lon(4))
                !write(*,*) temp_lat
                !write(*,*) temp_lon
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
                !Find grid position of the max and min coordinates and add2 grids*EMEP_grid_interpolation_size
                i_temp_min=1+floor((temp_x_min-temp_var1d_nc_dp(1,1))/temp_delta(1)+0.5)
                i_temp_max=1+floor((temp_x_max-temp_var1d_nc_dp(1,1))/temp_delta(1)+0.5)
                j_temp_min=1+floor((temp_y_min-temp_var1d_nc_dp(2,1))/temp_delta(2)+0.5)
                j_temp_max=1+floor((temp_y_max-temp_var1d_nc_dp(2,1))/temp_delta(2)+0.5)
                !write(unit_logfile,'(A,2I)') ' Reading EMEP i grids: ',i_temp_min,i_temp_max
                !write(unit_logfile,'(A,2I)') ' Reading EMEP j grids: ',j_temp_min,j_temp_max
                i_temp_min=max(1,i_temp_min-1-ceiling(1.*EMEP_grid_interpolation_size))
                i_temp_max=min(dim_length_nc(x_dim_nc_index),i_temp_max+1+ceiling(1.*EMEP_grid_interpolation_size))
                j_temp_min=max(1,j_temp_min-1-ceiling(1.*EMEP_grid_interpolation_size))
                j_temp_max=min(dim_length_nc(y_dim_nc_index),j_temp_max+1+ceiling(1.*EMEP_grid_interpolation_size))
                dim_length_nc(x_dim_nc_index)=i_temp_max-i_temp_min+1
                dim_length_nc(y_dim_nc_index)=j_temp_max-j_temp_min+1
                dim_start_nc(x_dim_nc_index)=i_temp_min
                dim_start_nc(y_dim_nc_index)=j_temp_min
                write(unit_logfile,'(A,3I)') ' Reading EMEP i grids: ',i_temp_min,i_temp_max,dim_length_nc(x_dim_nc_index)
                write(unit_logfile,'(A,3I)') ' Reading EMEP j grids: ',j_temp_min,j_temp_max,dim_length_nc(y_dim_nc_index)
            !endif

        endif
        

        !Allocate the nc arrays for reading
        if (.not.allocated(val_dim_nc)) allocate (val_dim_nc(maxval(dim_length_nc),num_dims_nc)) !x, y, z and time dimmension values
        if (.not.allocated(unit_dim_nc)) allocate (unit_dim_nc(num_dims_nc)) !x, y, z and time dimmension values
        if (.not.allocated(var1d_nc)) allocate (var1d_nc(maxval(dim_length_nc),num_dims_nc)) !x, y, z and time maximum dimmensions
        if (.not.allocated(var2d_nc)) allocate (var2d_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),2)) !Lat and lon
        if (.not.allocated(var3d_nc)) allocate (var3d_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index),num_var_nc,n_source_nc_index))
        if (.not.allocated(var4d_nc)) allocate (var4d_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index),num_var_nc,n_source_nc_index))
        if (.not.allocated(comp_var3d_nc)) allocate (comp_var3d_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index),n_compound_nc_index))
        if (.not.allocated(comp_var4d_nc)) allocate (comp_var4d_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index),n_compound_nc_index))
        if (.not.allocated(var1d_nc_dp)) allocate (var1d_nc_dp(maxval(dim_length_nc))) 
        if (.not.allocated(var2d_nc_dp)) allocate (var2d_nc_dp(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index))) !Lat and lon
        !if (.not.allocated(var3d_nc_dp)) allocate (var3d_nc_dp(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)))
        !allocate (var4d_nc_dp(dim_length_nc(x_index),dim_length_nc(y_index),1,dim_length_nc(time_index)))
        if (i_file.eq.2.and..not.allocated(lc_var3d_nc)) allocate (lc_var3d_nc(dim_length_nc(xdist_dim_nc_index),dim_length_nc(ydist_dim_nc_index),dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index),num_lc_var_nc,n_source_nc_index))
        if (i_file.eq.2.and..not.allocated(lc_var4d_nc)) allocate (lc_var4d_nc(dim_length_nc(xdist_dim_nc_index),dim_length_nc(ydist_dim_nc_index),dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),1,dim_length_nc(time_dim_nc_index),num_lc_var_nc,n_source_nc_index))
        if (i_file.eq.2.and..not.allocated(pm_lc_var4d_nc)) allocate (pm_lc_var4d_nc(dim_length_nc(xdist_dim_nc_index),dim_length_nc(ydist_dim_nc_index),dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),1,dim_length_nc(time_dim_nc_index),num_lc_var_nc,n_source_nc_index,2))
        if (.not.allocated(pm_var4d_nc)) allocate (pm_var4d_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index),num_var_nc,n_source_nc_index,2))

        !if (.not.allocated(temp_var4d_nc)) allocate (temp_var4d_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index)))

        !write(*,*) x_dim_nc_index,y_dim_nc_index
        !write(*,*) shape(var1d_nc_dp)
        !write(*,*) dim_length_nc
        !Read in the dimensions and check values of the dimmensions. Not necessary but diagnostic
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
                !    ,maxval(int((var1d_nc(1:dim_length_nc(i),i)-var1d_nc(dim_start_nc(i),i))/3600.+.5)+1)                     
            else
                write(unit_logfile,'(3A,2f12.2)') ' ',trim(dim_name_nc(i)),' (min, max): ' &
                    ,minval(var1d_nc(1:dim_length_nc(i),i)),maxval(var1d_nc(1:dim_length_nc(i),i)) 
            endif       
        enddo

        !Set the compound index (could also be looped if required)
        if (compound_index.eq.nox_nc_index) then
            n_compound_loop=3
            compound_loop_index(1)=nox_nc_index
            compound_loop_index(2)=no2_nc_index
            compound_loop_index(3)=o3_nc_index
        else
            n_compound_loop=1
            compound_loop_index(1)=compound_index
        endif
        
        i_conc=compound_index
        
        !Loop through the sources
        do i_source=1,n_source_nc_index
        if (calculate_source(i_source).or.i_source.eq.allsource_index) then
        !var_name_nc(num_var_nc,n_compound_nc_index,n_source_nc_index)
            
        !Loop through the variables
        do i=1,num_var_nc
            !write(*,*) i,trim(var_name_nc(i))
            !if (i.eq.frac_nc_index) var_name_nc_temp=var_name_nc(i,i_conc,i_source)
            !if (i.eq.conc_nc_index) var_name_nc_temp=var_name_nc(i,i_conc,i_source)
            
            !Identify the variable name and ID in the nc file
            var_name_nc_temp=var_name_nc(i,i_conc,i_source)
            status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
            !write(*,*) 'Status1: ',status_nc,var_id_nc,trim(var_name_nc_temp),i_source
            !Exception for pm10
            if ((i.eq.conc_nc_index.or.i.eq.frac_nc_index).and.i_conc.eq.pm10_nc_index) then
                var_name_nc_temp=var_name_nc(i,pmco_nc_index,i_source)
                status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                !write(*,*) '-------------------',trim(var_name_nc_temp),status_nc
            endif
            
            !If a variable name is found in the file then go further
            if (status_nc.eq.NF90_NOERR) then
                scale_factor_nc=1.
                !Find the dimmensions of the variable (temp_num_dims)
                status_nc = NF90_INQUIRE_VARIABLE(id_nc, var_id_nc, ndims = temp_num_dims)
                !write(*,*) temp_num_dims,status_nc
                if (temp_num_dims.eq.2.and.i_file.eq.1) then
                    !Read latitude and longitude data into a 2d grid if available. Only lat lon is 2d?
                    if (i.eq.lat_nc_index.or.i.eq.lon_nc_index) then
                    status_nc = NF90_GET_VAR (id_nc, var_id_nc, var2d_nc_dp);var2d_nc(:,:,i)=real(var2d_nc_dp)
                    write(unit_logfile,'(A,i3,A,2A,2f16.4)') ' Reading: ',temp_num_dims,' ',trim(var_name_nc_temp),' (min, max): ',minval(var2d_nc(:,:,i)),maxval(var2d_nc(:,:,i))
                    endif
                elseif (temp_num_dims.eq.3.and.i_file.eq.1) then
                    !write(*,'(6i)') dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),temp_start_time_nc_index,dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)
                    status_nc = NF90_GET_VAR (id_nc, var_id_nc, var3d_nc(:,:,:,i,i_source),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),dim_start_nc(time_dim_nc_index)/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                    !write(*,*) status_nc
                    write(unit_logfile,'(A,I,3A,2f16.4)') ' Reading: ',temp_num_dims,' ',trim(var_name_nc_temp),' (min, max): ',minval(var3d_nc(:,:,:,i,i_source)),maxval(var3d_nc(:,:,:,i,i_source))
                elseif (temp_num_dims.eq.4) then
                    if (i_file.eq.2.and.i.eq.ZTOP_nc_index) then
                        !Don't try to read
                    elseif (i.eq.conc_nc_index.and.i_conc.eq.pm10_nc_index) then
                        var_name_nc_temp2=var_name_nc(i,pmco_nc_index,i_source)
                        status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp2), var_id_nc)
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, pm_var4d_nc(:,:,dim_start_nc(z_dim_nc_index):dim_start_nc(z_dim_nc_index)+dim_length_nc(z_dim_nc_index)-1,:,i,i_source,1),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),dim_start_nc(z_dim_nc_index),dim_start_nc(time_dim_nc_index)/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,I,3A,2f16.4)') ' Reading: ',temp_num_dims,' ',trim(var_name_nc_temp2),' (min, max): ',minval(pm_var4d_nc(:,:,:,:,i,i_source,1)),maxval(pm_var4d_nc(:,:,:,:,i,i_source,1))
                        var_name_nc_temp2=var_name_nc(i,pm25_nc_index,i_source)
                        status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp2), var_id_nc)
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, pm_var4d_nc(:,:,dim_start_nc(z_dim_nc_index):dim_start_nc(z_dim_nc_index)+dim_length_nc(z_dim_nc_index)-1,:,i,i_source,2),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),dim_start_nc(z_dim_nc_index),dim_start_nc(time_dim_nc_index)/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,I,3A,2f16.4)') ' Reading: ',temp_num_dims,' ',trim(var_name_nc_temp2),' (min, max): ',minval(pm_var4d_nc(:,:,:,:,i,i_source,2)),maxval(pm_var4d_nc(:,:,:,:,i,i_source,2))
                        var4d_nc(:,:,:,:,i,i_source)=pm_var4d_nc(:,:,:,:,i,i_source,1)+pm_var4d_nc(:,:,:,:,i,i_source,2)
                    else
                    status_nc = NF90_GET_VAR (id_nc, var_id_nc, var4d_nc(:,:,dim_start_nc(z_dim_nc_index):dim_start_nc(z_dim_nc_index)+dim_length_nc(z_dim_nc_index)-1,:,i,i_source),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),dim_start_nc(z_dim_nc_index),dim_start_nc(time_dim_nc_index)/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                    !status_nc = NF90_GET_VAR (id_nc, var_id_nc, temp_var4d_nc(:,:,:,:),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),dim_start_nc(z_dim_nc_index),temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                    !var4d_nc(:,:,:,:,i,i_source)=real(temp_var4d_nc(:,:,:,:))
                    write(unit_logfile,'(A,I,3A,2f16.4)') ' Reading: ',temp_num_dims,' ',trim(var_name_nc_temp),' (min, max): ',minval(var4d_nc(:,:,:,:,i,i_source)),maxval(var4d_nc(:,:,:,:,i,i_source))
                    endif
                    !write(*,*) shape(var4d_nc)
                    !write(*,*) dim_start_nc(z_dim_nc_index),dim_length_nc(z_dim_nc_index)
                    !write(*,*) maxval(var4d_nc(:,:,1,1,i,i_source)),maxval(var4d_nc(:,:,1,2,i,i_source))
                elseif (temp_num_dims.eq.6.and.i_file.eq.2) then
                    if (i.eq.frac_nc_index.and.i_conc.eq.pm10_nc_index) then
                        var_name_nc_temp2=var_name_nc(i,pmco_nc_index,i_source)
                        status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp2), var_id_nc)
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, pm_lc_var4d_nc(:,:,:,:,:,:,lc_frac_nc_index,i_source,1),start=(/1,1,dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),dim_start_nc(z_dim_nc_index),dim_start_nc(time_dim_nc_index)/),count=(/dim_length_nc(xdist_dim_nc_index),dim_length_nc(ydist_dim_nc_index),dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,I,3A,2f16.4)') ' Reading: ',temp_num_dims,' ',trim(var_name_nc_temp2),' (min, max): ',minval(pm_lc_var4d_nc(:,:,:,:,:,:,lc_frac_nc_index,i_source,1)),maxval(pm_lc_var4d_nc(:,:,:,:,:,:,lc_frac_nc_index,i_source,1))
                        var_name_nc_temp2=var_name_nc(i,pm25_nc_index,i_source)
                        status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp2), var_id_nc)
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, pm_lc_var4d_nc(:,:,:,:,:,:,lc_frac_nc_index,i_source,2),start=(/1,1,dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),dim_start_nc(z_dim_nc_index),dim_start_nc(time_dim_nc_index)/),count=(/dim_length_nc(xdist_dim_nc_index),dim_length_nc(ydist_dim_nc_index),dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,I,3A,2f16.4)') ' Reading: ',temp_num_dims,' ',trim(var_name_nc_temp2),' (min, max): ',minval(pm_lc_var4d_nc(:,:,:,:,:,:,lc_frac_nc_index,i_source,2)),maxval(pm_lc_var4d_nc(:,:,:,:,:,:,lc_frac_nc_index,i_source,2))
                        !Not used but calculated for writing
                        lc_var4d_nc(:,:,:,:,:,:,lc_frac_nc_index,i_source)=pm_lc_var4d_nc(:,:,:,:,:,:,lc_frac_nc_index,i_source,1)+pm_lc_var4d_nc(:,:,:,:,:,:,lc_frac_nc_index,i_source,2)
                    else
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, lc_var4d_nc(:,:,:,:,:,:,lc_frac_nc_index,i_source),start=(/1,1,dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),dim_start_nc(z_dim_nc_index),dim_start_nc(time_dim_nc_index)/),count=(/dim_length_nc(xdist_dim_nc_index),dim_length_nc(ydist_dim_nc_index),dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index)/))
                        write(unit_logfile,'(A,I,3A,2f16.4)') ' Reading: ',temp_num_dims,' ',trim(var_name_nc_temp),' (min, max): ',minval(lc_var4d_nc(:,:,:,:,:,:,lc_frac_nc_index,i_source)),maxval(lc_var4d_nc(:,:,:,:,:,:,lc_frac_nc_index,i_source))
                    endif
                    !write(*,*) shape(lc_var4d_nc)
                    !write(*,*) maxval(lc_var4d_nc(3,3,:,:,:,1,lc_frac_nc_index,i_source)),maxval(lc_var4d_nc(3,3,:,:,:,2,lc_frac_nc_index,i_source))

                else
                    write(unit_logfile,'(8A,8A)') ' Cannot find a correct dimmension for: ',trim(var_name_nc_temp)
                endif    
                
            else
                 !write(unit_logfile,'(8A,8A)') ' Cannot read: ',trim(var_name_nc_temp)
            endif

        enddo
        
        
        endif
        enddo
        
        !Loop through the additional compounds that are in the base file
        !if (i_file.eq.1) then
        i_source=allsource_index
        i=conc_nc_index
        do i_loop=1,n_compound_loop
            i_conc=compound_loop_index(i_loop)
            var_name_nc_temp=comp_name_nc(i_conc)
            status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
            !write(*,*) 'Status1: ',status_nc,id_nc,var_id_nc,trim(var_name_nc_temp)
            
            !If a variable name is found in the file then go further
            if (status_nc.eq.NF90_NOERR) then

                !Find the dimmensions of the variable (temp_num_dims)
                status_nc = NF90_INQUIRE_VARIABLE(id_nc, var_id_nc, ndims = temp_num_dims)

                if (temp_num_dims.eq.3) then
                    if (i_file.eq.1) then
                    status_nc = NF90_GET_VAR (id_nc, var_id_nc, comp_var4d_nc(:,:,1,:,i_conc),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),1,temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),1,dim_length_nc(time_dim_nc_index)/))
                    comp_var4d_nc(:,:,1,:,i_conc)=comp_var4d_nc(:,:,1,:,i_conc)*comp_scale_nc(i_conc)
                    write(unit_logfile,'(A,I,3A,2f16.4)') ' Reading compound file 1: ',temp_num_dims,' ',trim(var_name_nc_temp),' (min, max): ',minval(comp_var4d_nc(:,:,1,:,i_conc)),maxval(comp_var4d_nc(:,:,1,:,i_conc))
                    elseif (i_file.eq.2) then
                    !In case the comp data is in the uEMEP file then read it here with no vertical extent
                    status_nc = NF90_GET_VAR (id_nc, var_id_nc, comp_var4d_nc(:,:,1,:,i_conc),start=(/dim_start_nc(x_dim_nc_index),dim_start_nc(y_dim_nc_index),1,temp_start_time_nc_index/),count=(/dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),1,dim_length_nc(time_dim_nc_index)/))
                    comp_var4d_nc(:,:,1,:,i_conc)=comp_var4d_nc(:,:,1,:,i_conc)*comp_scale_nc(i_conc)
                    write(unit_logfile,'(A,I,3A,2f16.4)') ' Reading compound file 2: ',temp_num_dims,' ',trim(var_name_nc_temp),' (min, max): ',minval(comp_var4d_nc(:,:,1,:,i_conc)),maxval(comp_var4d_nc(:,:,1,:,i_conc))
                    endif
                endif
                if (temp_num_dims.eq.4) then
                    if (i_file.eq.1) then
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
                
        enddo !compound loop
        !endif
            
        status_nc = NF90_CLOSE (id_nc)
        
        !Invert the base file arrays in the z direction to be compatible with the uEMEP files
        !Not implemented any more (1.eq.2)
        if (i_file.eq.1.and.invert_levels_flag) then
            write(unit_logfile,'(a)') 'WARNING: Inverting layers in base file. Temporary measure since the base file layers are ordered differently to the uEMEP layers'
            if (.not.allocated(swop_var4d_nc)) allocate (swop_var4d_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index),num_var_nc,n_source_nc_index))
            do i=1,dim_length_nc(z_dim_nc_index)
                swop_var4d_nc(:,:,i,:,:,:)=var4d_nc(:,:,dim_length_nc(z_dim_nc_index)-i+1,:,:,:)
            enddo
            var4d_nc=swop_var4d_nc
            if (allocated(swop_var4d_nc)) deallocate (swop_var4d_nc)
            if (.not.allocated(swop_comp_var4d_nc)) allocate (swop_comp_var4d_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index),n_compound_nc_index))
            do i=1,dim_length_nc(z_dim_nc_index)
                swop_comp_var4d_nc(:,:,i,:,:)=comp_var4d_nc(:,:,dim_length_nc(z_dim_nc_index)-i+1,:,:)
            enddo
            comp_var4d_nc=swop_comp_var4d_nc
            if (allocated(swop_comp_var4d_nc)) deallocate (swop_comp_var4d_nc)
        endif
        
            
    enddo !End file loop
    
    !Set the correct time dimmensions tot he first file value
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
    
       
        !Transfer all source values to all the sources for use in source looping later
    do i_source=1,n_source_nc_index
        if (calculate_source(i_source).and.i_source.ne.allsource_nc_index) then
        var3d_nc(:,:,:,conc_nc_index,i_source)=var3d_nc(:,:,:,conc_nc_index,allsource_nc_index)
        var4d_nc(:,:,:,:,conc_nc_index,i_source)=var4d_nc(:,:,:,:,conc_nc_index,allsource_nc_index)
        endif
    enddo
        
        !Transfer local contribution 4d to 3d since this is the only one currently used
        !write(*,*) shape(var4d_nc)
        !write(*,*) surface_level_nc
        lc_var3d_nc=lc_var4d_nc(:,:,:,:,surface_level_nc_2,:,:,:)
        !var3d_nc(:,:,:,conc_nc_index,:)=var4d_nc(:,:,surface_level_nc,:,conc_nc_index,:)
        var3d_nc(:,:,:,conc_nc_index,:)=var4d_nc(:,:,surface_level_nc_2,:,conc_nc_index,:)
        comp_var3d_nc(:,:,:,:)=comp_var4d_nc(:,:,surface_level_nc,:,:)
        
        if (allocated(lc_var4d_nc)) deallocate(lc_var4d_nc)
        if (allocated(comp_var4d_nc)) deallocate(comp_var4d_nc)

        !Use average lowest level grid thickness
        H_emep_temp=sum(var4d_nc(:,:,surface_level_nc_2,:,ZTOP_nc_index,allsource_index))/dim_length_nc(x_dim_nc_index)/dim_length_nc(y_dim_nc_index)/dim_length_nc(time_dim_nc_index)
        if (H_emep_temp.ne.0) then
            H_emep=H_emep_temp
            write(unit_logfile,'(A,f8.2)')' Using model depth info. Setting lowest level depth = ',H_emep
        else
            write(unit_logfile,'(A,f8.2)')' No model depth info available. Setting lowest level depth = ',H_emep
        endif
    
        H_meteo=H_emep/2.
        write(unit_logfile,'(A,f8.2)')' Setting lowest meteo grid height = ',H_meteo
        
        var3d_nc(:,:,:,conc_nc_index,:)=var4d_nc(:,:,surface_level_nc_2,:,conc_nc_index,:)
        
        !At the moment the local contribution based on fraction. Convert to local contributions here
        !Remove this if we read local contributions in a later version
        do j=1,dim_length_nc(ydist_dim_nc_index)
        do i=1,dim_length_nc(xdist_dim_nc_index)
            if (compound_index.eq.pm10_nc_index) then
                lc_var3d_nc(i,j,:,:,:,lc_local_nc_index,:)=pm_var4d_nc(:,:,surface_level_nc_2,:,conc_nc_index,:,1)*pm_lc_var4d_nc(i,j,:,:,surface_level_nc_2,:,lc_frac_nc_index,:,1) &
                    +pm_var4d_nc(:,:,surface_level_nc_2,:,conc_nc_index,:,2)*pm_lc_var4d_nc(i,j,:,:,surface_level_nc_2,:,lc_frac_nc_index,:,2)
            else
            !lc_var4d_nc(i,j,:,:,:,:,lc_local_nc_index,:)=var4d_nc(:,:,:,:,conc_nc_index,:)*lc_var4d_nc(i,j,:,:,:,:,lc_frac_nc_index,:)
            lc_var3d_nc(i,j,:,:,:,lc_local_nc_index,:)=var3d_nc(:,:,:,conc_nc_index,:)*lc_var3d_nc(i,j,:,:,:,lc_frac_nc_index,:)
            endif
        enddo
        enddo
        
        if (allocated(pm_lc_var4d_nc)) deallocate(pm_lc_var4d_nc)
        if (allocated(pm_var4d_nc)) deallocate(pm_var4d_nc)
        
       !Set the local grid contribution for the individual grid
        !write(*,*) shape(lc_var4d_nc)
        !write(*,*) shape(var4d_nc)
        !Commented out as these are not used?
        !var4d_nc(:,:,:,:,frac_nc_index,:)=lc_var4d_nc(xdist_centre_nc,ydist_centre_nc,:,:,:,:,lc_frac_nc_index,:)
        var3d_nc(:,:,:,frac_nc_index,:)=lc_var3d_nc(xdist_centre_nc,ydist_centre_nc,:,:,:,lc_frac_nc_index,:)
        !var4d_nc(:,:,:,:,local_nc_index,:)=var4d_nc(:,:,:,:,conc_nc_index,:)*var4d_nc(:,:,:,:,frac_nc_index,:)
        var3d_nc(:,:,:,local_nc_index,:)=var3d_nc(:,:,:,conc_nc_index,:)*var3d_nc(:,:,:,frac_nc_index,:)
        
        !write(*,*) minval(var4d_nc(:,:,:,:,local_nc_index,:)),maxval(var4d_nc(:,:,:,:,local_nc_index,:))
        !write(*,*) minval(var4d_nc(:,:,:,:,frac_nc_index,:)),maxval(var4d_nc(:,:,:,:,frac_nc_index,:))
        !write(*,*) minval(var3d_nc(:,:,:,frac_nc_index,:)),maxval(var3d_nc(:,:,:,frac_nc_index,:))
        !write(*,*) minval(var4d_nc(:,:,:,:,conc_nc_index,:)),maxval(var4d_nc(:,:,:,:,conc_nc_index,:))
        !write(*,*) minval(var3d_nc(:,:,:,conc_nc_index,:)),maxval(var3d_nc(:,:,:,conc_nc_index,:))
        !write(*,*) minval(var3d_nc(:,:,:,inv_FF10_nc_index,allsource_index)),maxval(var3d_nc(:,:,:,inv_FF10_nc_index,allsource_index))
        !write(*,*) minval(var3d_nc(:,:,:,FF10_nc_index,allsource_index)),maxval(var3d_nc(:,:,:,FF10_nc_index,allsource_index))
        
            
        !If no logz0 available. Set to log(0.1)
        !For urban areas a value of 0.3 is used
        where (var3d_nc(:,:,:,logz0_nc_index,:).eq.0.0) var3d_nc(:,:,:,logz0_nc_index,:)=log(0.3)
        if (replace_z0.ne.NODATA_value) then
            write(unit_logfile,'(A,f8.4)') ' Replacing z0 everywhere with: ',replace_z0
            var3d_nc(:,:,:,logz0_nc_index,:)=log(replace_z0)
        endif
        if (replace_invL.ne.NODATA_value) then
            write(unit_logfile,'(A,f8.4)') ' Replacing inverse L everywhere with: ',replace_invL
            var3d_nc(:,:,:,invL_nc_index,:)=replace_invL
        endif
        if (replace_hmix.ne.NODATA_value) then
            write(unit_logfile,'(A,f8.4)') ' Replacing HMIX everywhere with: ',replace_hmix
            var3d_nc(:,:,:,hmix_nc_index,:)=replace_hmix
        endif
        if (FF_scale.ne.NODATA_value) then
            write(unit_logfile,'(A,f8.4)') ' Rescaling wind fields everywhere with factor: ',FF_scale
            var3d_nc(:,:,:,ustar_nc_index,:)=var3d_nc(:,:,:,ustar_nc_index,:)*FF_scale
            var3d_nc(:,:,:,FF10_nc_index,:)=var3d_nc(:,:,:,FF10_nc_index,:)*FF_scale
            var3d_nc(:,:,:,inv_FF10_nc_index,:)=var3d_nc(:,:,:,inv_FF10_nc_index,:)/FF_scale
            var3d_nc(:,:,:,ugrid_nc_index,:)=var3d_nc(:,:,:,ugrid_nc_index,:)*FF_scale
            var3d_nc(:,:,:,vgrid_nc_index,:)=var3d_nc(:,:,:,vgrid_nc_index,:)*FF_scale
            var3d_nc(:,:,:,inv_FFgrid_nc_index,:)=var3d_nc(:,:,:,inv_FFgrid_nc_index,:)/FF_scale
            var4d_nc(:,:,:,:,ugrid_nc_index,:)=var4d_nc(:,:,:,:,ugrid_nc_index,:)*FF_scale
            var4d_nc(:,:,:,:,vgrid_nc_index,:)=var4d_nc(:,:,:,:,vgrid_nc_index,:)*FF_scale
            var4d_nc(:,:,:,:,inv_FFgrid_nc_index,:)=var4d_nc(:,:,:,:,inv_FFgrid_nc_index,:)/FF_scale
        endif
        if (FF10_offset.ne.NODATA_value) then
            write(unit_logfile,'(A,f8.4)') ' Offsetting 10 m wind fields everywhere with a value: ',FF10_offset
            var3d_nc(:,:,:,FF10_nc_index,:)=var3d_nc(:,:,:,FF10_nc_index,:)+FF10_offset
         endif
        if (DD_offset.ne.NODATA_value) then
            write(unit_logfile,'(A,f8.4)') ' Rotating wind fields everywhere with a value: ',DD_offset
            if (.not.allocated(temp_var4d_nc)) allocate (temp_var4d_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(z_dim_nc_index),dim_length_nc(time_dim_nc_index),2))

            !Make use of the spare source index parts of the array for the conversion
            temp_var4d_nc(:,:,:,:,1) = var4d_nc(:,:,:,:,ugrid_nc_index,allsource_index)*cos(DD_offset/180.*3.14159)+var4d_nc(:,:,:,:,vgrid_nc_index,allsource_index)*sin(DD_offset/180.*3.14159)                                       
            temp_var4d_nc(:,:,:,:,2) =-var4d_nc(:,:,:,:,ugrid_nc_index,allsource_index)*sin(DD_offset/180.*3.14159)+var4d_nc(:,:,:,:,vgrid_nc_index,allsource_index)*cos(DD_offset/180.*3.14159)                                       
            var4d_nc(:,:,:,:,ugrid_nc_index,allsource_nc_index) = temp_var4d_nc(:,:,:,:,1)
            var4d_nc(:,:,:,:,vgrid_nc_index,allsource_nc_index) = temp_var4d_nc(:,:,:,:,2)
        endif
       
        !Set the magnitude of the gridded wind fields. Should probably be done after subgridding?
        var4d_nc(:,:,:,:,FFgrid_nc_index,allsource_index)=sqrt(var4d_nc(:,:,:,:,ugrid_nc_index,allsource_index)**2+var4d_nc(:,:,:,:,vgrid_nc_index,allsource_index)**2)
  
        !Check EMEP time
        date_num_temp=dble(ceiling(val_dim_nc(1,time_dim_nc_index)*24.))/24.
        call number_to_date(date_num_temp,date_array,ref_year_EMEP)
        write(unit_logfile,'(a,i6)') ' Time dimmension EMEP:  ',dim_length_nc(time_dim_nc_index)
        write(unit_logfile,'(a,6i6)') ' Date start EMEP =  ',date_array
        date_num_temp=dble(ceiling(val_dim_nc(dim_length_nc(time_dim_nc_index),time_dim_nc_index)*24.))/24.
        call number_to_date(date_num_temp,date_array,ref_year_EMEP)
        write(unit_logfile,'(a,6i6)') ' Date end EMEP =    ',date_array
        
        !Test and correct dates
        if (1.eq.1) then
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
        enddo
        !stop
        endif
        
        if (allocated(var1d_nc_dp)) deallocate (var1d_nc_dp)
        if (allocated(var2d_nc_dp)) deallocate (var2d_nc_dp)
        if (allocated(var3d_nc_dp)) deallocate (var3d_nc_dp)
        if (allocated(var4d_nc_dp)) deallocate (var4d_nc_dp)
        if (allocated(swop_var4d_nc)) deallocate (swop_var4d_nc)
        if (allocated(swop_comp_var4d_nc)) deallocate (swop_comp_var4d_nc)
        if (allocated(temp_var4d_nc)) deallocate (temp_var4d_nc)
 
    end subroutine uEMEP_read_EMEP
    
    