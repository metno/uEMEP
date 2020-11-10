!Saves receptor data in netcdf format
    
    subroutine uEMEP_save_netcdf_receptor_file(unit_logfile_in,filename_netcdf,nx,ny,nt_in,val_array_in,x_array,y_array,lon_array,lat_array,name_array,unit_array,title_str,create_file,valid_min &
        ,x_rec,y_rec,lon_rec,lat_rec,height_rec,name_rec_in,nr,variable_type,scale_factor)
    
    use uEMEP_definitions
    use netcdf
    
    implicit none
    
    character(256) filename_netcdf,name_array,unit_array,title_str,temp_name,temp_name3(3)
    integer unit_logfile_in
    integer nx,ny,nt_in,nr
    real val_array(nx,ny,nt_in),val_array_in(nx,ny,nt_in)!,val_array_temp(nx,ny,nt_in)
    real x_array(nx,ny)
    real y_array(nx,ny)
    real lon_array(nx,ny)
    real lat_array(nx,ny)!,lat_array_temp(nx,ny)
    !real time_array(nt_in)
    !real x_vector(nx)
    !real y_vector(ny)
    logical create_file
    real valid_min
    character(256) variable_type
    real scale_factor
    
    integer ncid
    integer station_dimid,lat_dimid,lon_dimid,val_dimid,time_dimid,charlen_dimid
    integer station_varid,station_name_varid,lat_varid,lon_varid,val_varid,time_varid,proj_varid,x_varid,y_varid,height_varid
    integer dimids3(3),dimids2(2)
    integer n_dims_length(3),n_dims_start(3)
    integer status
    integer tr,rr
    real x_rec(nr),y_rec(nr),height_rec(nr)
    real lon_rec(nr),lat_rec(nr)
    character(256) name_rec_in(nr)
    character(256) temp_char
    integer n_char
    !parameter (n_char=7)
    parameter (n_char=64)
    character(1) name_rec(n_char,nr)
    integer n_time_total
    real val_rec(nr,nt_in)
    real delta(2)
    real area_weighted_interpolation_function
    integer id_rec(nr)
    integer nf90_type
    integer nt
    integer(8) time_seconds_output_nc(nt_in)
    integer tr_0

    !Do not save if no receptor position data is available
    if (nr.eq.0) then
        return
    endif
    
    if (trim(variable_type).eq.'byte') nf90_type=NF90_BYTE
    if (trim(variable_type).eq.'float') nf90_type=NF90_FLOAT
    if (trim(variable_type).eq.'double') nf90_type=NF90_DOUBLE
    
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
                val_array(:,:,nt)=val_array_av(:,:,counter_av)/end_time_loop_index
                time_seconds_output_nc(nt)=time_seconds_output_av(counter_av)/end_time_loop_index
            endif
            !write(unit_logfile_in,*) 'Saving as average single time loop (nt,counter_av):',nt,counter_av
        else
            !write(unit_logfile_in,*) 'Saving as average multiple time loop (nt,counter_av):',nt,counter_av,time_seconds_output_nc(1),time_seconds_output_nc(nt)
            !write(*,*) time_seconds_output_nc(1:nt)
            val_array_av(:,:,counter_av)=sum(val_array(:,:,1:nt),3)/nt
            time_seconds_output_av(counter_av)=sum(time_seconds_output_nc(1:nt),1)/nt
            !write(*,*) sum(time_seconds_output_nc(1:nt),1)
            nt=1
            val_array(:,:,nt)=val_array_av(:,:,counter_av)
            time_seconds_output_nc(nt)=time_seconds_output_av(counter_av)
            !write(unit_logfile_in,*) 'Saving as average multiple time loop (nt,counter_av):',nt,counter_av,time_seconds_output_nc(nt),time_seconds_output_av(counter_av)
        endif
        
    endif

    !Interpolate to receptor position given the input array
    !write(unit_logfile,'(a)')' Interpolating to receptor point '
    !Assumes 2 elements to an array
    delta(1)=(x_array(2,1)-x_array(1,1))
    delta(2)=(y_array(1,2)-y_array(1,1))
    do rr=1,nr
    do tr=1,nt
        val_rec(rr,tr)=area_weighted_interpolation_function(x_array,y_array,val_array(:,:,tr),nx,ny,delta,x_rec(rr),y_rec(rr))
    enddo
    enddo
    
    !Make the receptor name to fit to netcdf requirements
    do rr=1,nr
    temp_char=name_rec_in(rr)
    do tr=1,n_char
        name_rec(tr,rr)=temp_char(tr:tr)
        !write(*,*) trim(name_rec_in(rr)),trim(name_rec(rr))
    enddo
        tr_0=len_trim(temp_char ) 
        name_rec(tr_0+1,rr)=char(0)
        id_rec(rr)=rr
    enddo
    
    if (save_netcdf_average_flag) then  
        n_time_total=1
    else
        n_time_total=end_time_nc_index-start_time_nc_index+1
    endif
        
    
    
    if (create_file) then
        !Create a netcdf file
        !call check(  nf90_create(filename_netcdf, nf90_clobber, ncid) )
        !call check(  nf90_create(filename_netcdf, NF90_HDF5, ncid) )
        call check(  nf90_create(filename_netcdf, IOR(NF90_HDF5, NF90_CLASSIC_MODEL), ncid) ) !New

        !Specify global attributes
        call check(  nf90_put_att(ncid, nf90_global, "Conventions", "CF-1.6" ) )
        call check(  nf90_put_att(ncid, nf90_global, "title", trim(title_str)) )
        call check(  nf90_put_att(ncid, nf90_global, "Model", "uEMEP" ) )        
        call check(  nf90_put_att(ncid, nf90_global, "featureType", "timeSeries" ) )        

        !Projection data
        call check(  nf90_def_var(ncid, "projection_utm", NF90_int, proj_varid) )
        call check(  nf90_put_att(ncid, proj_varid, "semi_major_axis", 6378137.0 ) )
        call check(  nf90_put_att(ncid, proj_varid, "inverse_flattening", 298.257222101 ) )

        call check(  nf90_put_att(ncid, proj_varid, "grid_mapping_name", "transverse_mercator" ) )
        call check(  nf90_put_att(ncid, proj_varid, "scale_factor_at_central_meridian", 0.9996 ) )
        call check(  nf90_put_att(ncid, proj_varid, "latitude_of_projection_origin", 0 ) )
        call check(  nf90_put_att(ncid, proj_varid, "false_easting", 500000. ) )
        call check(  nf90_put_att(ncid, proj_varid, "false_northing", 0. ) )
        call check(  nf90_put_att(ncid, proj_varid, "longitude_of_central_meridian", utm_lon0 ) )
        !call check(  nf90_put_att(ncid, proj_varid, "semi_major_axis", 6378140.0 ) )
        !call check(  nf90_put_att(ncid, proj_varid, "semi_minor_axis", 6356750.0 ) )
  
        !Define the dimensions for the entire dataset
        !write(*,*) 'n_valid_receptor_in',n_valid_receptor_in
        !write(*,*) 'n_valid_receptor',n_valid_receptor
        call check(  nf90_def_dim(ncid,"station_id",n_valid_receptor_in, station_dimid) )
        call check(  nf90_def_dim(ncid,"charlen",n_char, charlen_dimid) )
        
        !call check(  nf90_def_dim(ncid,"time",n_time_total, time_dimid) )
        !To have time as unlimittec (Heiko)
        call check(  nf90_def_dim(ncid,"time",NF90_UNLIMITED, time_dimid) )
        
     
        !Define the dimension variables
        call check(  nf90_def_var(ncid, "station_id", NF90_INT, station_dimid, station_varid) )
        !call check(  nf90_def_var(ncid, "time", NF90_DOUBLE, time_dimid, time_varid) )
        call check(  nf90_def_var(ncid, "time", NF90_INT, time_dimid, time_varid) )
    
        !Define the values
        dimids2 = (/ station_dimid, time_dimid /)
        !dimids1 = (/ station_dimid /)
        call check(  nf90_def_var(ncid, "lat", NF90_REAL, station_dimid, lat_varid) )
        call check(  nf90_def_var(ncid, "lon", NF90_REAL, station_dimid, lon_varid) )
        call check(  nf90_def_var(ncid, "y", NF90_REAL, station_dimid, y_varid) )
        call check(  nf90_def_var(ncid, "x", NF90_REAL, station_dimid, x_varid) )
        call check(  nf90_def_var(ncid, "station_name", NF90_CHAR, (/charlen_dimid,station_dimid/), station_name_varid) )
        call check(  nf90_def_var(ncid, "station_height", NF90_REAL, station_dimid, height_varid) )
        !call check(  nf90_def_var(ncid, "station_name", NF90_CHAR, (/station_dimid/), station_name_varid) )
        !call check(  nf90_def_var(ncid, "station_name", NF90_CHAR, (/charlen_dimid,station_dimid/), station_name_varid) )

        !Specify the units
        call check(  nf90_put_att(ncid, lat_varid, "units", "degrees_north") )
        call check(  nf90_put_att(ncid, lon_varid, "units", "degrees_east") )
        call check(  nf90_put_att(ncid, y_varid, "units", "m") )
        call check(  nf90_put_att(ncid, x_varid, "units", "m") )
        call check(  nf90_put_att(ncid, height_varid, "units", "m") )
        call check(  nf90_put_att(ncid, time_varid, "units", trim(unit_dim_nc(time_dim_nc_index))) )
        call check(  nf90_put_att(ncid, station_varid, "long_name", "station index" ) )
        call check(  nf90_put_att(ncid, station_name_varid, "long_name", "station name" ) )
        call check(  nf90_put_att(ncid, station_name_varid, "cf_role", "timeseries_id" ) )

        !Specify other dimension attributes
        !call check(  nf90_put_att(ncid, y_varid, "standard_name", "projection_y_coordinate") )
        !call check(  nf90_put_att(ncid, x_varid, "standard_name", "projection_x_coordinate") )
        !call check(  nf90_put_att(ncid, y_varid, "axis", "Y") )
        !call check(  nf90_put_att(ncid, x_varid, "axis", "X") )
  
        !Close the definitions
        call check( nf90_enddef(ncid) )

        !call check( nf90_put_var(ncid, time_varid, val_dim_nc(1:dim_length_nc(time_dim_nc_index),time_dim_nc_index)) )
        !call check( nf90_put_var(ncid, station_varid, name_rec(:,1:n_char) )
        !call check( nf90_put_var(ncid, lat_varid, lat_array) )
        !call check( nf90_put_var(ncid, lon_varid, lon_array) )
        
        !Put this in to fix unlimitted problem? If works on annual need to check it works other places as well!
        !call check( nf90_put_var(ncid, time_varid, time_seconds_output_nc(1:nt), start=(/1/), count=(/nt/)) )
        !call check( nf90_put_var(ncid, time_varid, time_seconds_output_nc(1:nt)) )
        
        call check( nf90_close(ncid) )
    
    endif
        
    !Add to the existing file
    call check( nf90_open(filename_netcdf, NF90_WRITE, ncid) )
    
    !Get the dimensions id from the existing file
    call check( nf90_inq_dimid(ncid,"time",time_dimid) )
    call check( nf90_inq_dimid(ncid, "station_id", station_dimid) )
   ! write(*,*) 'station_dimid ',station_dimid
   ! write(*,*) 'time_dimid ',time_dimid
    dimids2 = (/ station_dimid, time_dimid /)

    !Get the size of the dimensions
    call check( nf90_inquire_dimension(ncid, dimids2(1), temp_name, n_dims_length(1)) )
    call check( nf90_inquire_dimension(ncid, dimids2(2), temp_name, n_dims_length(2)) )
    !Set the starting point to 1
    n_dims_start(1:2)=1
    
    
    !Set time to full length in unlimitted case (Heiko)
    n_dims_length(2) = n_time_total

   ! write(*,*) 'n_dims_length(1) ',n_dims_length(1)
   ! write(*,*) 'n_dims_length(2) ',n_dims_length(2)

    status=nf90_inq_varid(ncid, trim(name_array), val_varid)
    if (status.ne.nf90_NoErr) then
        !if the variable does not exist then create a new one
        !write(*,*) 'Creating new: ',trim(name_array)
        call check( nf90_redef(ncid) )
        
        call check( nf90_def_var(ncid, trim(name_array), nf90_type, dimids2, val_varid) )
        call check( nf90_put_att(ncid, val_varid, "units", trim(unit_array)) )
    
        !Specify other variable attributes
        if (nf90_type.eq.NF90_byte) then
            call check(  nf90_put_att(ncid, val_varid, "missing_value", int1(NODATA_value) ) ) !New
            call check(  nf90_put_att(ncid, val_varid, "valid_min", int1(valid_min)) )
        else
            call check(  nf90_put_att(ncid, val_varid, "missing_value", NODATA_value ) ) !New
            call check(  nf90_put_att(ncid, val_varid, "valid_min", valid_min) )
        endif
        call check(  nf90_put_att(ncid, val_varid, "grid_mapping", "projection_utm") )
        call check(  nf90_put_att(ncid, val_varid, "coordinates", "station_name lat lon") )
        if (scale_factor.ne.1.) call check(  nf90_put_att(ncid, val_varid, "scale_factor", scale_factor) )
        
        !Close the definitions
        call check( nf90_enddef(ncid) )
    endif


    
    !Should not be used but can be
    if (use_single_time_loop_flag) then
        !Add time to the time dimension       
        !call check( nf90_inq_varid(ncid, "time", time_varid) )
        !call check( nf90_inquire_dimension(ncid, time_dimid, temp_name, n_dims(3)) )
        !n_dims(3)=n_dims(3)+1
    
        if (save_netcdf_average_flag) then
            n_dims_start(2)=1
            n_dims_length(2)=1
        else
            n_dims_start(2)=t_loop
            n_dims_length(2)=1
        endif
        !call check( nf90_put_var(ncid, time_varid, val_dim_nc(1,time_dim_nc_index), start = (/n_dims(2)/) ) )
        !write(*,*) n_dims(3),val_dim_nc(1,time_dim_nc_index)
        !write(*,*) n_dims
        
        !Add dimension and array to existing
        !call check( nf90_inq_varid(ncid, trim(name_array), val_varid) )       
        !call check( nf90_put_var(ncid, val_varid, val_rec, start=(/1,n_dims(2)/), count=(/n_dims(1),1/)) )

        !call check( nf90_put_var(ncid, time_varid, val_dim_nc(1:dim_length_nc(time_dim_nc_index),time_dim_nc_index), start=(/n_dims(2)/), count=(/1/)) )
        !call check( nf90_put_var(ncid, station_varid, name_rec(1,:), start=(/1,1/), count=(/n_dims(1),n_char/)) )

    elseif (use_multiple_receptor_grids_flag) then
        !Add time to the time dimension       
        !call check( nf90_inq_varid(ncid, "station", station_varid) )
        !call check( nf90_inquire_dimension(ncid, time_dimid, temp_name, n_dims(3)) )
        !n_dims(3)=n_dims(3)+1
        
        n_dims_start(1)=valid_receptor_inverse_index(g_loop)
        n_dims_length(1)=1
        id_rec(1)=valid_receptor_inverse_index(g_loop)
        !write(*,*) n_dims(3),val_dim_nc(1,time_dim_nc_index)
        !write(*,*) n_dims
        
        !Add dimension and array to existing
        !call check( nf90_inq_varid(ncid, trim(name_array), val_varid) )              
        !call check( nf90_put_var(ncid, val_varid, val_rec, start=(/n_dims(1),1/), count=(/1,n_dims(2)/)) )

        !call check( nf90_put_var(ncid, time_varid, val_dim_nc(1:dim_length_nc(time_dim_nc_index),time_dim_nc_index) , start=(/1/), count=(/n_dims(2)/)) )
        !call check( nf90_put_var(ncid, station_varid, name_rec(1,:), start = (/n_dims(1),1/), count=(/1,n_char/)) )
    endif
    
    !write(*,*) 'n_dims_start',n_dims_start
    !write(*,*) 'n_dims_length',n_dims_length
    
        !Fill in the complete dimension variables for time and receptor names
        call check( nf90_inq_varid(ncid, "station_id", station_varid) )
        call check( nf90_inq_varid(ncid, "time", time_varid) )
        call check( nf90_inq_varid(ncid, "station_name", station_name_varid) )
        call check( nf90_inq_varid(ncid, "x", x_varid) )
        call check( nf90_inq_varid(ncid, "y", y_varid) )
        call check( nf90_inq_varid(ncid, "lon", lon_varid) )
        call check( nf90_inq_varid(ncid, "lat", lat_varid) )
        call check( nf90_inq_varid(ncid, "station_height", height_varid) )
      
        !Write time to the file

        !!call check( nf90_put_var(ncid, time_varid, val_dim_nc(1:dim_length_nc(time_dim_nc_index),time_dim_nc_index), start=(/n_dims_start(2)/), count=(/n_dims_length(2)/)) )
        call check( nf90_put_var(ncid, time_varid, time_seconds_output_nc(1:nt), start=(/n_dims_start(2)/), count=(/n_dims_length(2)/)) )
        !!call check( nf90_put_var(ncid, station_varid, name_rec(:), start = (/1,1/), count=(/n_dims(1),n_char/)) )
        
        !Write station index and name
        call check( nf90_put_var(ncid, station_varid, id_rec, start = (/n_dims_start(1),1/), count=(/n_dims_length(1),n_char/)) )
        call check( nf90_put_var(ncid, station_name_varid, name_rec, start = (/1,n_dims_start(1)/), count=(/n_char,n_dims_length(1)/)) )
        
        !!call check( nf90_put_var(ncid, station_name_varid, name_rec, start = (/1,1/), count=(/n_char,n_dims(1)/)) )
        !!call check( nf90_put_var(ncid, station_name_varid, name_rec, start = (/1/), count=(/n_dims(1)/)) )
      
        !Write the variable to file
        if (nf90_type.eq.NF90_byte) then
            call check( nf90_put_var(ncid, val_varid, int1(val_rec(:,1:nt)), start = (/n_dims_start(1),n_dims_start(2)/), count=(/n_dims_length(1),n_dims_length(2)/)) )
        else
            call check( nf90_put_var(ncid, val_varid, val_rec(:,1:nt), start = (/n_dims_start(1),n_dims_start(2)/), count=(/n_dims_length(1),n_dims_length(2)/)) ) 
        endif

        !Write position data to the file
        call check( nf90_put_var(ncid, x_varid, x_rec, start = (/n_dims_start(1)/), count=(/n_dims_length(1)/)) )
        call check( nf90_put_var(ncid, y_varid, y_rec, start = (/n_dims_start(1)/), count=(/n_dims_length(1)/)) )
        call check( nf90_put_var(ncid, lon_varid, lon_rec, start = (/n_dims_start(1)/), count=(/n_dims_length(1)/)) )
        call check( nf90_put_var(ncid, lat_varid, lat_rec, start = (/n_dims_start(1)/), count=(/n_dims_length(1)/)) )
        call check( nf90_put_var(ncid, height_varid, height_rec, start = (/n_dims_start(1)/), count=(/n_dims_length(1)/)) )


    
        
    call check( nf90_close(ncid) )
    

    
    end subroutine uEMEP_save_netcdf_receptor_file
    
 