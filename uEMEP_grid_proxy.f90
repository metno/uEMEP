!uEMEP_grid_proxy.f90
    
!==========================================================================
!   uEMEP model save_gridded_proxy
!   Calculates proxy concentrations based on Gaussian distribution
!   Generic for all sources and subsources
!==========================================================================
    subroutine uEMEP_grid_proxy(source_index)
    
    use uEMEP_definitions

    implicit none
    
    integer         source_index
    character(256)  temp_name
    logical         exists
    integer         jj,ii,tt,tt_emis
    real            distance_subgrid
    real            dx_temp,dy_temp
    integer         i_start,i_end,j_start,j_end,t_start,t_end
    integer         i_cross,j_cross
    integer         i_cross_integral,j_cross_integral
    real            u_subgrid_loc,v_subgrid_loc,kz_subgrid_loc
    real            cos_subgrid_loc,sin_subgrid_loc
    integer         i_nc,j_nc
    integer         subsource_index
    real            ay_loc,by_loc,az_loc,bz_loc,sig_y_0_loc,sig_z_0_loc,h_emis_loc,z_rec_loc
    real            lon_limit,lat_limit
    real            time_weight(subgrid_dim(t_dim_index)),time_total(subgrid_dim(t_dim_index))
    integer         time_count(subgrid_dim(t_dim_index))
    real            x_pos,y_pos
    real            lon_area_min,lon_area_max,lat_area_min,lat_area_max
    real            distance_subgrid_min
    real xpos_subgrid,ypos_subgrid
    real xpos_emission_subgrid,ypos_emission_subgrid
    real temp_subgrid_internal
    
    !integer ii_temp,jj_temp
    real, allocatable :: temp_emission_subgrid(:,:)
    real, allocatable :: temp_subgrid(:,:)
    real, allocatable :: temp_FF_subgrid(:,:)
    real, allocatable :: trajectory_subgrid(:,:,:,:)
    real, allocatable :: angle_diff(:,:)
    !real, allocatable :: temp_lat_emission_subgrid(:,:)
    !real, allocatable :: temp_x_emission_subgrid(:,:)
    !real, allocatable :: temp_y_emission_subgrid(:,:)
    !real, allocatable :: temp_distance_subgrid(:,:)
    !real, allocatable :: temp_vector(:,:)
    !logical, allocatable :: logical_vector(:,:)
    !integer i_dim,j_dim
    !real temp_val

    !logical :: use_vector=.false.
    
    !real, allocatable :: r(:,:)
    !real, allocatable :: sig_th(:,:)
    !real, allocatable :: sig_z(:,:)
    !real, allocatable :: B(:,:)
    !real order_1,order_2

    !real, allocatable :: emission_subgrid_temp(:,:,:)
    integer traj_max_index
    logical valid_traj
    real traj_step_size,x_loc,y_loc
    real z0_temp,h_temp
    
    !functions
    real gauss_plume_second_order_rotated_func
    real gauss_plume_second_order_rotated_integral_func
    real gauss_plume_cartesian_func
    real gauss_plume_cartesian_integral_func
    real gauss_plume_cartesian_trajectory_func
    !real, allocatable :: gauss_plume_second_order_rotated_vector_func(:,:)
 
    write(unit_logfile,'(A)') ''
    write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Calculating dispersion of proxy (uEMEP_grid_proxy)'
	write(unit_logfile,'(A)') '================================================================'

    !First call the integral routine if it is needed. Only when using the concentration redistribution
    if (local_subgrid_method_flag.eq.1.or.EMEP_grid_interpolation_flag.eq.4) then
        call uEMEP_grid_proxy_integral(source_index)
    endif
    
    allocate (temp_emission_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index))) 
    allocate (temp_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index))) 
    allocate (temp_FF_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index))) 
    allocate (angle_diff(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index))) 

    temp_subgrid=0.
    temp_emission_subgrid=0.
    temp_FF_subgrid=0.
    
    !subgrid(:,:,:,proxy_integral_subgrid_index,source_index,:)=0.
    
        
    if (read_existing_grid_data(proxy_file_index(source_index))) then
        if (combine_emission_subsources_during_dispersion(source_index).and.n_subsource(source_index).gt.1) then
            n_subsource(source_index)=1
        endif
        do subsource_index=1,n_subsource(source_index)
            !Read in proxy file
            temp_name=trim(pathname_grid(proxy_file_index(source_index)))//trim(filename_grid(proxy_file_index(source_index)))//trim(subsource_str(subsource_index))//'_'//trim(file_tag)//'.asc'
            inquire(file=temp_name,exist=exists)
            if (.not.exists) then
                write(unit_logfile,*)'ERROR: '//trim(temp_name)//' does not exist.'
                return
            endif
            call read_esri_ascii_3d_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),subgrid_delta(1), &
                subgrid(:,:,:,proxy_subgrid_index,source_index,subsource_index),x_subgrid,y_subgrid)
             
        enddo
        return
    endif
 

    lon_limit=dgrid_nc(lon_nc_index)/2.*EMEP_grid_interpolation_size
    lat_limit=dgrid_nc(lat_nc_index)/2.*EMEP_grid_interpolation_size

    !Minimum distance for travel time calculation
    distance_subgrid_min=sqrt(subgrid_delta(x_dim_index)*subgrid_delta(x_dim_index)+subgrid_delta(y_dim_index)*subgrid_delta(y_dim_index))/2./sqrt(2.)*4./3.14159
    
    

    if (use_downwind_position_flag) traj_max_index=traj_max_index*2
    
    !write(*,*) traj_step_scale,traj_step_size,traj_max_index,integral_subgrid_loop_index(x_dim_index),integral_subgrid_loop_index(y_dim_index)
    
    do subsource_index=1,n_subsource(source_index)

        call uEMEP_set_dispersion_params_simple(source_index,subsource_index)
    
        !Set local dispersion parameters
        ay_loc=ay(source_index,subsource_index)
        by_loc=by(source_index,subsource_index)
        az_loc=az(source_index,subsource_index)
        bz_loc=bz(source_index,subsource_index)
        sig_y_0_loc=sig_y_0(source_index,subsource_index)
        sig_z_0_loc=sig_z_0(source_index,subsource_index)
        h_emis_loc=h_emis(source_index,subsource_index)
        z_rec_loc=z_rec(source_index,subsource_index)

        !write(*,*) z_rec_loc,ay_loc,by_loc,az_loc,bz_loc,sig_y_0_loc,sig_z_0_loc
    
        write(unit_logfile,'(a,i3)')'Calculating proxy concentration data for '//trim(source_file_str(source_index))//' with subsource index ',subsource_index

        !Set the start and end times of the loop
        t_start=1
        t_end=subgrid_dim(t_dim_index)
        
        !Copy emission_subgrid before looping for imporved performance
        !Not necessary with the time loop outside
        !do tt=t_start,t_end
        !    emission_subgrid_temp(tt,:,:)=emission_subgrid(:,:,tt,source_index,subsource_index)
        !enddo
                
        !Loop through the time
        do tt=t_start,t_end
        
        subgrid(:,:,tt,proxy_subgrid_index,source_index,:)=0.
        
        !Set a temporary emission array
        temp_emission_subgrid=emission_subgrid(:,:,tt,source_index,subsource_index)
        temp_subgrid=0.
        temp_FF_subgrid=0.
        
        if (.not.use_single_time_loop_flag) then
            if (tt.gt.t_start) then
                last_meteo_subgrid(:,:,:)=meteo_subgrid(:,:,tt-1,:)
            else
                last_meteo_subgrid(:,:,:)=meteo_subgrid(:,:,tt,:)      
            endif           
        endif

        !Precalculate information for the trajectory model
        !Maxium number of trajectory steps and size of steps based on the integral (meteorology) loop size
        if (use_trajectory_flag) then
            traj_step_size=min(integral_subgrid_delta(x_dim_index),integral_subgrid_delta(y_dim_index))*traj_step_scale
            traj_max_index=floor(max(integral_subgrid_loop_index(x_dim_index),integral_subgrid_loop_index(y_dim_index))/traj_step_scale)
            if (tt.eq.t_start) write(unit_logfile,'(a,f12.1,i)') 'Trajectory step (m) and dimensions: ', traj_step_size, traj_max_index
            if (.not.allocated(trajectory_subgrid)) allocate(trajectory_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),traj_max_index,2))
            
            trajectory_subgrid=NODATA_value
            
            !Loop through the emissions and create trajectories for all emissions source grids
            do j=1,emission_subgrid_dim(y_dim_index,source_index)
            do i=1,emission_subgrid_dim(x_dim_index,source_index)
            
                if (temp_emission_subgrid(i,j).ne.0) then
                    call uEMEP_calculate_all_trajectory(x_emission_subgrid(i,j,source_index),y_emission_subgrid(i,j,source_index),tt, &
                        traj_max_index,traj_step_size,trajectory_subgrid(i,j,:,x_dim_index),trajectory_subgrid(i,j,:,y_dim_index))
                    !write(*,'(50f12.1)') x_emission_subgrid(i,j,source_index),trajectory_subgrid(i,j,:,x_dim_index)
                    !write(*,'(50f12.1)') y_emission_subgrid(i,j,source_index),trajectory_subgrid(i,j,:,y_dim_index)
                endif
            
            enddo
            enddo
       
        endif

        !Create a temporary wind speed subgrid for each hour
        do j_cross=1,integral_subgrid_dim(y_dim_index)
        do i_cross=1,integral_subgrid_dim(x_dim_index)
            z0_temp=exp(meteo_subgrid(i_cross,j_cross,tt,logz0_subgrid_index))
            h_temp=h_emis(source_index,subsource_index)
            if (annual_calculations.and.wind_level_flag.eq.1) then
                temp_FF_subgrid(i_cross,j_cross)=1./meteo_subgrid(i_cross,j_cross,tt,inv_FFgrid_subgrid_index)
            elseif (annual_calculations.and.wind_level_flag.eq.2) then
                temp_FF_subgrid(i_cross,j_cross)=1./meteo_subgrid(i_cross,j_cross,tt,inv_FFgrid_subgrid_index)*(1.-(log((H_emep/2.+z0_temp)/z0_temp)-log((h_temp+z0_temp)/z0_temp))/log((H_emep/2.+z0_temp)/z0_temp))
            elseif (annual_calculations.and.wind_level_flag.eq.3) then
                temp_FF_subgrid(i_cross,j_cross)=1./meteo_subgrid(i_cross,j_cross,tt,inv_FF10_subgrid_index)
            elseif (annual_calculations.and.wind_level_flag.eq.4) then
                temp_FF_subgrid(i_cross,j_cross)=1./meteo_subgrid(i_cross,j_cross,tt,inv_FF10_subgrid_index)*(1.-(log((10.+z0_temp)/z0_temp)-log((h_temp+z0_temp)/z0_temp))/log((10.+z0_temp)/z0_temp))
            elseif (hourly_calculations.and.wind_level_flag.eq.1) then
                temp_FF_subgrid(i_cross,j_cross)=meteo_subgrid(i_cross,j_cross,tt,FFgrid_subgrid_index)
            elseif (hourly_calculations.and.wind_level_flag.eq.2) then
                temp_FF_subgrid(i_cross,j_cross)=meteo_subgrid(i_cross,j_cross,tt,FFgrid_subgrid_index)*(1.-(log((H_emep/2.+z0_temp)/z0_temp)-log((h_temp+z0_temp)/z0_temp))/log((H_emep/2.+z0_temp)/z0_temp))
            elseif (hourly_calculations.and.wind_level_flag.eq.3) then
                temp_FF_subgrid(i_cross,j_cross)=meteo_subgrid(i_cross,j_cross,tt,FF10_subgrid_index)
            elseif (hourly_calculations.and.wind_level_flag.eq.4) then
                temp_FF_subgrid(i_cross,j_cross)=meteo_subgrid(i_cross,j_cross,tt,FF10_subgrid_index)*(1.-(log((10.+z0_temp)/z0_temp)-log((h_temp+z0_temp)/z0_temp))/log((10.+z0_temp)/z0_temp))
            elseif (wind_level_flag.eq.0) then
                temp_FF_subgrid(i_cross,j_cross)=1.
            else
                write(unit_logfile,'(a)') 'No valid wind_level_flag selected. Stopping (uEMEP_grid_proxy)'
                stop
            endif
            
            !Setting a minimum value for wind for dispersion purposes (cannot be zero)
            temp_FF_subgrid(i_cross,j_cross)=sqrt(temp_FF_subgrid(i_cross,j_cross)*temp_FF_subgrid(i_cross,j_cross)+FF_min_dispersion*FF_min_dispersion)
            
            if (temp_FF_subgrid(i_cross,j_cross).eq.0) then
                write(unit_logfile,'(a,2i)') 'Zero wind speed at integral grid (stopping): ',i_cross,j_cross
                stop
            endif
            
            !Finds the angle difference between the current and last meteo field for dispersion.
            !Still need to centre the wind vector, which is not currently done.
            if (use_last_meteo_in_dispersion) then
                cos_subgrid_loc=meteo_subgrid(i_cross,j_cross,tt,cos_subgrid_index)
                sin_subgrid_loc=meteo_subgrid(i_cross,j_cross,tt,sin_subgrid_index)
                angle_diff(i_cross,j_cross)=abs(asin(meteo_subgrid(i_cross,j_cross,tt,sin_subgrid_index)*last_meteo_subgrid(i_cross,j_cross,cos_subgrid_index) &
                                           -meteo_subgrid(i_cross,j_cross,tt,cos_subgrid_index)*last_meteo_subgrid(i_cross,j_cross,sin_subgrid_index) ))/2.
                
                angle_diff(i_cross,j_cross)=min(angle_diff(i_cross,j_cross),3.14159/4.) !Maximum allowed is 45 degrees
                !write(*,*) i_cross,j_cross,angle_diff(i_cross,j_cross)*180/3.14159
                                        
            else
                angle_diff(i_cross,j_cross)=0.
            endif

            !Implement a meandering that is dependent on wind speed. Max 30 degrees at FF_min_dispersion (0.5) and around 10 degrees at 5*FF_min_dispersion (2.5)
            if (use_meandering_in_dispersion) then
                angle_diff(i_cross,j_cross)=angle_diff(i_cross,j_cross) &
                +(30.*3.14159/180.)*exp(-(temp_FF_subgrid(i_cross,j_cross)-FF_min_dispersion)/(FF_min_dispersion*2.))
            endif

        enddo
        enddo
        
        !Loop through the proxy grid
        do j=1,subgrid_dim(y_dim_index)
        do i=1,subgrid_dim(x_dim_index)
               
            !Only use those that are marked for use
            if (use_subgrid(i,j,source_index)) then    
             
                if (EMEP_projection_type.eq.LL_projection_index) then
                    xpos_subgrid=lon_subgrid(i,j)
                    ypos_subgrid=lat_subgrid(i,j)
                elseif (EMEP_projection_type.eq.LCC_projection_index) then
                    call lb2lambert_uEMEP(xpos_subgrid,ypos_subgrid,lon_subgrid(i,j),lat_subgrid(i,j),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                endif
                
                !Use the wind direction to move the target area downwind
                if (use_downwind_position_flag.and.hourly_calculations) then
                    
                    !Find the cross reference to the target grid from the emission grid
                    i_cross=crossreference_target_to_emission_subgrid(i,j,x_dim_index,source_index)
                    j_cross=crossreference_target_to_emission_subgrid(i,j,y_dim_index,source_index)
                    !Find the meteorological data
                    i_cross_integral=crossreference_target_to_integral_subgrid(i,j,x_dim_index)
                    j_cross_integral=crossreference_target_to_integral_subgrid(i,j,y_dim_index)                    
                    x_pos=max(-1.,min(1.,meteo_subgrid(i_cross_integral,j_cross_integral,tt,cos_subgrid_index)*sqrt(2.)))
                    y_pos=max(-1.,min(1.,meteo_subgrid(i_cross_integral,j_cross_integral,tt,sin_subgrid_index)*sqrt(2.)))                    
                    i_end=min(ceiling(i_cross+(1.-x_pos)*emission_subgrid_loop_index(x_dim_index,source_index)),emission_subgrid_dim(x_dim_index,source_index))
                    i_start=max(floor(i_cross-(1.+x_pos)*emission_subgrid_loop_index(x_dim_index,source_index)),1)
                    j_end=min(ceiling(j_cross+(1.-y_pos)*emission_subgrid_loop_index(y_dim_index,source_index)),emission_subgrid_dim(y_dim_index,source_index))
                    j_start=max(floor(j_cross-(1.+y_pos)*emission_subgrid_loop_index(y_dim_index,source_index)),1)
                                                    
                    !Set new lon and lat limits to twice their normal size to include the upwind source region
                    !lon_limit=dgrid_nc(lon_nc_index)*EMEP_grid_interpolation_size
                    !lat_limit=dgrid_nc(lat_nc_index)*EMEP_grid_interpolation_size
                    lon_area_max=xpos_subgrid+(1.-x_pos)*lon_limit
                    lon_area_min=xpos_subgrid-(1.+x_pos)*lon_limit
                    lat_area_max=ypos_subgrid+(1.-y_pos)*lat_limit
                    lat_area_min=ypos_subgrid-(1.+y_pos)*lat_limit

                else
                    
                    !Find the cross reference to the target grid from the emission grid
                    i_cross=crossreference_target_to_emission_subgrid(i,j,x_dim_index,source_index)
                    j_cross=crossreference_target_to_emission_subgrid(i,j,y_dim_index,source_index)
                    
                    !Set the size of the loop region around the target cell to be up to subgrid_loop_index
                    i_start=max(1,i_cross-emission_subgrid_loop_index(x_dim_index,source_index))
                    i_end=min(emission_subgrid_dim(x_dim_index,source_index),i_cross+emission_subgrid_loop_index(x_dim_index,source_index))
                    j_start=max(1,j_cross-emission_subgrid_loop_index(y_dim_index,source_index))
                    j_end=min(emission_subgrid_dim(y_dim_index,source_index),j_cross+emission_subgrid_loop_index(y_dim_index,source_index))
 
                    lon_area_max=xpos_subgrid+lon_limit
                    lon_area_min=xpos_subgrid-lon_limit
                    lat_area_max=ypos_subgrid+lat_limit
                    lat_area_min=ypos_subgrid-lat_limit

                endif
 
    

                time_weight(tt)=0.
                !time_count=0
                time_total(tt)=0.
                
                !i_start=1;i_end=emission_subgrid_dim(x_dim_index,source_index);j_start=1;j_end=emission_subgrid_dim(y_dim_index,source_index)

                !Loop through emission sub_grids in the nearby region
                !j_end=j_start;i_end=i_start
                do jj=j_start,j_end
                do ii=i_start,i_end
                    
                    !Limit distance to the size of the EMEP grid
                    if (temp_emission_subgrid(ii,jj).ne.0) then
                        !.and.abs(lon_emission_subgrid(ii,jj,source_index)-lon_subgrid(i,j)).le.lon_limit &
                        !.and.abs(lat_emission_subgrid(ii,jj,source_index)-lat_subgrid(i,j)).le.lat_limit) then
                            
                        if (EMEP_projection_type.eq.LL_projection_index) then
                            xpos_emission_subgrid=lon_emission_subgrid(ii,jj,source_index)
                            ypos_emission_subgrid=lat_emission_subgrid(ii,jj,source_index)
                        elseif (EMEP_projection_type.eq.LCC_projection_index) then
                            call lb2lambert_uEMEP(xpos_emission_subgrid,ypos_emission_subgrid,lon_emission_subgrid(ii,jj,source_index),lat_emission_subgrid(ii,jj,source_index),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                        endif

                        !write(*,*) xpos_emission_subgrid,lon_area_min,lon_area_max
                        !write(*,*) ypos_emission_subgrid,lat_area_min,lat_area_max

                        !if (xpos_emission_subgrid.gt.lon_area_min.and.xpos_emission_subgrid.lt.lon_area_max &
                        !    .and.ypos_emission_subgrid.gt.lat_area_min.and.ypos_emission_subgrid.lt.lat_area_max) then                  
                        if (xpos_emission_subgrid.ge.lon_area_min.and.xpos_emission_subgrid.le.lon_area_max &
                            .and.ypos_emission_subgrid.ge.lat_area_min.and.ypos_emission_subgrid.le.lat_area_max) then                  
                       
                            !write(*,*) xpos_emission_subgrid-lon_area_min,xpos_emission_subgrid-lon_area_max,ypos_emission_subgrid-lat_area_min,ypos_emission_subgrid-lat_area_max
                        !Loop through the time
                        !do tt=t_start,t_end
                        if (hourly_calculations) then 
                            
                            if (use_trajectory_flag) then
                                
                                !call uEMEP_local_trajectory(i,j,ii,jj,tt,traj_max_index,traj_step_size,source_index,x_loc,y_loc,valid_traj)
                                !call uEMEP_local_trajectory(x_subgrid(i,j),y_subgrid(i,j),x_emission_subgrid(ii,jj,source_index),y_emission_subgrid(ii,jj,source_index),tt,traj_max_index,traj_step_size,x_loc,y_loc,valid_traj)
                                call uEMEP_minimum_distance_trajectory(x_subgrid(i,j),y_subgrid(i,j),x_emission_subgrid(ii,jj,source_index),y_emission_subgrid(ii,jj,source_index),tt, &
                                traj_max_index,traj_step_size,trajectory_subgrid(ii,jj,:,x_dim_index),trajectory_subgrid(ii,jj,:,y_dim_index),x_loc,y_loc,valid_traj)
                               
                                if (valid_traj) then
                                 !subgrid(i,j,tt,proxy_subgrid_index,source_index,subsource_index)=subgrid(i,j,tt,proxy_subgrid_index,source_index,subsource_index) &
                                  !  +emission_subgrid(ii,jj,tt,source_index,subsource_index) &
                                  !  *gauss_plume_cartesian_trajectory_func(x_loc,y_loc,h_emis_loc,z_rec_loc,ay_loc,by_loc,az_loc,bz_loc,sig_y_0_loc,sig_z_0_loc)
                                 
                                    i_cross_integral=crossreference_emission_to_integral_subgrid(ii,jj,x_dim_index,source_index)
                                    j_cross_integral=crossreference_emission_to_integral_subgrid(ii,jj,y_dim_index,source_index)                            
                                    
                                    if (stability_scheme_flag.eq.1) call uEMEP_set_dispersion_params_simple(source_index,subsource_index)
                                    if (stability_scheme_flag.eq.2) call uEMEP_set_dispersion_params_PG(meteo_subgrid(i_cross_integral,j_cross_integral,tt,invL_subgrid_index),source_index,subsource_index)
                                                                        
                                    temp_subgrid_internal=temp_emission_subgrid(ii,jj) &
                                        *gauss_plume_cartesian_trajectory_func(x_loc,y_loc,h_emis_loc,z_rec_loc, &
                                        ay(source_index,subsource_index),by(source_index,subsource_index),az(source_index,subsource_index),bz(source_index,subsource_index), &
                                        sig_y_0(source_index,subsource_index),sig_z_0(source_index,subsource_index),angle_diff(i_cross_integral,j_cross_integral)) &
                                        /temp_FF_subgrid(i_cross_integral,j_cross_integral)
 
                                    temp_subgrid(i,j)=temp_subgrid(i,j)+temp_subgrid_internal

                                    !distance_subgrid=x_loc
                                    distance_subgrid=sqrt(x_loc*x_loc+y_loc*y_loc)
                                 
                                    !temp_subgrid(i,j)=temp_subgrid(i,j) &
                                    !    +temp_emission_subgrid(ii,jj) &
                                    !    *gauss_plume_cartesian_func(0.,0.,h_emis_loc,1.,0., &
                                    !    x_loc,y_loc,z_rec_loc, &
                                    !    ay(source_index,subsource_index),by(source_index,subsource_index),az(source_index,subsource_index),bz(source_index,subsource_index), &
                                    !    sig_y_0(source_index,subsource_index),sig_z_0(source_index,subsource_index)) &
                                    !    /temp_FF_subgrid(i_cross_integral,j_cross_integral)

                                    !write(*,*) x_loc, y_loc
                                    
                                    !if (x_loc.eq.0.and.y_loc.eq.0) write(*,*) x_subgrid(i,j)-x_emission_subgrid(ii,jj,source_index),y_subgrid(i,j)-y_emission_subgrid(ii,jj,source_index)
                                    !write(*,*) i,j,ii,jj,temp_subgrid(i,j)
                                else
                                    !write(*,*) 'Not a valid trajectory',i,j,ii,jj
                                    temp_subgrid_internal=0.
                                endif
                                
                            else
                                
                                !subgrid(i,j,tt,proxy_subgrid_index,source_index,subsource_index)=subgrid(i,j,tt,proxy_subgrid_index,source_index,subsource_index) &
                                !    +emission_subgrid(ii,jj,tt,source_index,subsource_index) &
                                !    *gauss_plume_cartesian_func(x_emission_subgrid(ii,jj,source_index),y_emission_subgrid(ii,jj,source_index),h_emis_loc,cos_subgrid_loc,sin_subgrid_loc, &
                                !    x_subgrid(i,j),y_subgrid(i,j),z_rec_loc,ay_loc,by_loc,az_loc,bz_loc,sig_y_0_loc,sig_z_0_loc)
                                
                                    i_cross_integral=crossreference_emission_to_integral_subgrid(ii,jj,x_dim_index,source_index)
                                    j_cross_integral=crossreference_emission_to_integral_subgrid(ii,jj,y_dim_index,source_index)                            
                                    cos_subgrid_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,cos_subgrid_index)
                                    sin_subgrid_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,sin_subgrid_index)                                    

                                    if (stability_scheme_flag.eq.1) call uEMEP_set_dispersion_params_simple(source_index,subsource_index)
                                    if (stability_scheme_flag.eq.2) call uEMEP_set_dispersion_params_PG(meteo_subgrid(i_cross_integral,j_cross_integral,tt,invL_subgrid_index),source_index,subsource_index)

                                    temp_subgrid_internal=temp_emission_subgrid(ii,jj) &
                                        *gauss_plume_cartesian_func(x_emission_subgrid(ii,jj,source_index),y_emission_subgrid(ii,jj,source_index),h_emis_loc,cos_subgrid_loc,sin_subgrid_loc, &
                                        x_subgrid(i,j),y_subgrid(i,j),z_rec_loc, &
                                        ay(source_index,subsource_index),by(source_index,subsource_index),az(source_index,subsource_index),bz(source_index,subsource_index), &
                                        sig_y_0(source_index,subsource_index),sig_z_0(source_index,subsource_index),angle_diff(i_cross_integral,j_cross_integral)) &
                                        /temp_FF_subgrid(i_cross_integral,j_cross_integral)
                            
                                    temp_subgrid(i,j)=temp_subgrid(i,j)+temp_subgrid_internal

                                    distance_subgrid=sqrt((x_emission_subgrid(ii,jj,source_index)-x_subgrid(i,j))*(x_emission_subgrid(ii,jj,source_index)-x_subgrid(i,j)) &
                                        +(y_emission_subgrid(ii,jj,source_index)-y_subgrid(i,j))*(y_emission_subgrid(ii,jj,source_index)-y_subgrid(i,j)))
                                
                                    !write(*,*) i,j,ii,jj,temp_subgrid(i,j)
                            
                            endif
                            
                            !Calculate weighted time based on the selected temp_FF_subgrid wind level
                            !Minimum distance is half of a grid diagonal weighted so the circle has the same area as the square with that diagonal
                            !if (subgrid(i,j,tt,proxy_subgrid_index,source_index,subsource_index).gt.0) then
                            if (temp_subgrid_internal.gt.0) then
                                    distance_subgrid=max(distance_subgrid,distance_subgrid_min)
                                    !Take weighted average (weighted by concentration) of the inverse of the time. Inverse because we want to weight the smaller values (I think)
                                    !time_weight(tt)=time_weight(tt)+meteo_subgrid(i_cross_integral,j_cross_integral,tt,FF10_subgrid_index)/distance_subgrid*subgrid(i,j,tt,proxy_subgrid_index,source_index,subsource_index)
                                    i_cross_integral=crossreference_emission_to_integral_subgrid(ii,jj,x_dim_index,source_index)
                                    j_cross_integral=crossreference_emission_to_integral_subgrid(ii,jj,y_dim_index,source_index)                            

                                    !time_weight(tt)=time_weight(tt)+meteo_subgrid(i_cross_integral,j_cross_integral,tt,FF10_subgrid_index)/distance_subgrid*temp_subgrid(i,j)
                                    !time_weight(tt)=time_weight(tt)+temp_FF_subgrid(i_cross_integral,j_cross_integral)/distance_subgrid*temp_subgrid_internal
                                    !Not inverse of the time
                                    time_weight(tt)=time_weight(tt)+distance_subgrid/temp_FF_subgrid(i_cross_integral,j_cross_integral)*temp_subgrid_internal
                                    
                                    !write(*,*) time_weight(tt),meteo_subgrid(i_cross_integral,j_cross_integral,tt,FF10_subgrid_index),distance_subgrid,temp_subgrid(i,j)
                                    !Calculate sum of the concentration for normalisation
                                    !time_total(tt)=time_total(tt)+subgrid(i,j,tt,proxy_subgrid_index,source_index,subsource_index)
                                    time_total(tt)=time_total(tt)+temp_subgrid_internal
                            endif

                        else
                            
                            distance_subgrid=sqrt((x_emission_subgrid(ii,jj,source_index)-x_subgrid(i,j))*(x_emission_subgrid(ii,jj,source_index)-x_subgrid(i,j)) &
                                +(y_emission_subgrid(ii,jj,source_index)-y_subgrid(i,j))*(y_emission_subgrid(ii,jj,source_index)-y_subgrid(i,j)))
                            
                            !subgrid(i,j,tt,proxy_subgrid_index,source_index,subsource_index)=subgrid(i,j,tt,proxy_subgrid_index,source_index,subsource_index) &
                            !    +emission_subgrid(ii,jj,tt,source_index,subsource_index) &
                            !    *gauss_plume_second_order_rotated_func(distance_subgrid,z_rec_loc,ay_loc,by_loc,az_loc,bz_loc,sig_y_0_loc,sig_z_0_loc,h_emis_loc)
                            
                            !Divide by wind speed at emission position
                            i_cross_integral=crossreference_emission_to_integral_subgrid(ii,jj,x_dim_index,source_index)
                            j_cross_integral=crossreference_emission_to_integral_subgrid(ii,jj,y_dim_index,source_index)                            

                            temp_subgrid(i,j)=temp_subgrid(i,j) &
                                +temp_emission_subgrid(ii,jj) &
                                *gauss_plume_second_order_rotated_func(distance_subgrid,z_rec_loc,ay_loc,by_loc,az_loc,bz_loc,sig_y_0_loc,sig_z_0_loc,h_emis_loc) &
                                /temp_FF_subgrid(i_cross_integral,j_cross_integral)

                            
                        endif
                        
                        !enddo !time loop
                        
                        endif
                    endif
                    
                    
                enddo
                enddo
                                
                !Add to the travel time
                !traveltime_subgrid(i,j,1:subgrid_dim(t_dim_index),1)=traveltime_subgrid(i,j,1:subgrid_dim(t_dim_index),1)+time_weight
                !traveltime_subgrid(i,j,1:subgrid_dim(t_dim_index),2)=traveltime_subgrid(i,j,1:subgrid_dim(t_dim_index),2)+time_total
                traveltime_subgrid(i,j,tt,1)=traveltime_subgrid(i,j,tt,1)+time_weight(tt)
                traveltime_subgrid(i,j,tt,2)=traveltime_subgrid(i,j,tt,2)+time_total(tt)
                !traveltime_subgrid(i,j,tt,1)=time_weight(tt)
                !traveltime_subgrid(i,j,tt,2)=time_total(tt)
                !time_weight=time_weight/time_total
                !Invert to get the time
                !time_weight=1./time_weight
                !Invert to get the time
                !time_weight=-log(time_weight)
                !Set none valid to 0. This happens if no emissions are available. Set to 12 hours (long time)
                !where (time_total.eq.0) time_weight=3600.*12.
                !comp_subgrid(i,j,:,traveltime_index)=time_weight
                !write(*,'(2i4,4f12.4)') i,j,temp_subgrid(i,j),time_weight(tt),traveltime_subgrid(i,j,tt,1),traveltime_subgrid(i,j,tt,2)
            
            else
                temp_subgrid(i,j)=NODATA_value
            endif
            
        enddo
        !if (mod(j,10).eq.0) write(*,'(3a,i5,a,i5,a,i3,a,i3)') 'Gridding ',trim(source_file_str(source_index)),' proxy',j,' of ',subgrid_dim(2),' and ',subsource_index,' of ',n_subsource(source_index)
        enddo
      
        if (mod(j,1).eq.0) write(*,'(3a,i5,a,i5,a,i3,a,i3)') 'Gridding ',trim(source_file_str(source_index)),' proxy for hour ',tt,' of ',subgrid_dim(t_dim_index),' and ',subsource_index,' of ',n_subsource(source_index)
        
        !Put the temporary subgrid back into the subgrid array
        subgrid(:,:,tt,proxy_subgrid_index,source_index,subsource_index)=temp_subgrid
        
        enddo !time loop
        
    enddo !subsource_index
               
        if (combine_emission_subsources_during_dispersion(source_index).and.n_subsource(source_index).gt.1) then
            do subsource_index=2,n_subsource(n_source_index)
                subgrid(:,:,:,proxy_subgrid_index,source_index,1)=subgrid(:,:,:,proxy_subgrid_index,source_index,1)+subgrid(:,:,:,proxy_subgrid_index,source_index,subsource_index)
                !subgrid(:,:,:,proxy_integral_subgrid_index,source_index,1)=subgrid(:,:,:,proxy_integral_subgrid_index,source_index,1)+subgrid(:,:,:,proxy_integral_subgrid_index,source_index,subsource_index)
            enddo
            n_subsource(source_index)=1
        endif
    
            
        if (save_intermediate_files) then
        if (.not.read_existing_grid_data(proxy_file_index(source_index))) then
            do subsource_index=1,n_subsource(n_source_index)
            temp_name=trim(pathname_grid(proxy_file_index(source_index)))//trim(filename_grid(proxy_file_index(source_index)))//trim(subsource_str(subsource_index))//'_'//trim(file_tag)//'.asc'
            write(unit_logfile,'(a)')'Writing to: '//trim(temp_name)
            call write_esri_ascii_3d_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index), &
                subgrid_delta(1),subgrid(:,:,:,proxy_subgrid_index,source_index,subsource_index),x_subgrid,y_subgrid)
            !temp_name=trim(pathname_grid(proxy_integral_file_index(source_index)))//trim(filename_grid(proxy_integral_file_index(source_index)))//trim(subsource_str(subsource_index))//'_'//trim(file_tag)//'.asc'
            !write(unit_logfile,'(a)')'Writing to: '//trim(temp_name)
            !call write_esri_ascii_3d_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index), &
            !    subgrid_delta(1),subgrid(:,:,:,proxy_integral_subgrid_index,source_index,subsource_index),x_subgrid,y_subgrid)
            enddo
        endif
        endif
  
    !deallocate (emission_subgrid_temp)
    !deallocate (temp_lon_emission_subgrid)
    !deallocate (temp_lat_emission_subgrid)
    !deallocate (temp_x_emission_subgrid)
    !deallocate (temp_y_emission_subgrid)
    if (allocated(trajectory_subgrid)) deallocate(trajectory_subgrid)
    if (allocated(temp_emission_subgrid)) deallocate(temp_emission_subgrid)
    if (allocated(temp_subgrid)) deallocate(temp_subgrid)
    if (allocated(temp_FF_subgrid)) deallocate(temp_FF_subgrid)

    end subroutine uEMEP_grid_proxy
    
    
!==========================================================================
!   uEMEP model grid_proxy_integral
!   Calculates the vertically integrated proxy concentrations based on Gaussian distribution
!   Generic for all sources and subsources
!   Seperate routine to the target grid as the integrated form may be of lower resolution
!==========================================================================
    subroutine uEMEP_grid_proxy_integral(source_index)
    
    use uEMEP_definitions

    implicit none
    
    integer         source_index
    character(256)  temp_name
    logical         exists
    integer         jj,ii,tt,tt_emis
    real            distance_subgrid
    real            dx_temp,dy_temp
    integer         i_start,i_end,j_start,j_end,t_start,t_end
    integer         i_cross,j_cross
    integer         i_cross_integral,j_cross_integral
    real            u_subgrid_loc,v_subgrid_loc,kz_subgrid_loc
    real            cos_subgrid_loc,sin_subgrid_loc
    integer         i_nc,j_nc
    integer         subsource_index
    real            ay_loc,by_loc,az_loc,bz_loc,sig_y_0_loc,sig_z_0_loc,h_emis_loc,z_rec_loc
    real            lon_limit,lat_limit
    real            x_pos,y_pos
    real            lon_area_min,lon_area_max,lat_area_min,lat_area_max
    real xpos_subgrid,ypos_subgrid
    real xpos_emission_subgrid,ypos_emission_subgrid
    real xpos_integral_subgrid,ypos_integral_subgrid

    integer traj_max_index
    logical valid_traj
    real traj_step_size,x_loc,y_loc
    real, allocatable :: temp_FF_subgrid(:,:)
    real, allocatable :: angle_diff(:,:)
    real z0_temp,h_temp
   
    !functions
    real gauss_plume_second_order_rotated_func
    real gauss_plume_second_order_rotated_integral_func
    real gauss_plume_cartesian_func
    real gauss_plume_cartesian_integral_func
    real gauss_plume_cartesian_trajectory_integral_func
    
    !Leave the subroutine if the integral is not required
    !if (local_subgrid_method_flag.ne.1) return
    
    lon_limit=dgrid_nc(lon_nc_index)/2.*EMEP_grid_interpolation_size
    lat_limit=dgrid_nc(lat_nc_index)/2.*EMEP_grid_interpolation_size
      
    !Read in previously saved data and leave the programme
    if (read_existing_grid_data(proxy_file_index(source_index))) then
        if (combine_emission_subsources_during_dispersion(source_index).and.n_subsource(source_index).gt.1) then
            n_subsource(source_index)=1
        endif
        
        do subsource_index=1,n_subsource(source_index)
            
            temp_name=trim(pathname_grid(proxy_integral_file_index(source_index)))//trim(filename_grid(proxy_integral_file_index(source_index)))//trim(subsource_str(subsource_index))//'_'//trim(file_tag)//'.asc'
            inquire(file=temp_name,exist=exists)
            if (.not.exists) then
                write(unit_logfile,*)'ERROR: '//trim(temp_name)//' does not exist.'
                return
            endif
            call read_esri_ascii_3d_file(unit_logfile,temp_name,integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),integral_subgrid_dim(t_dim_index),integral_subgrid_delta(1), &
                integral_subgrid(:,:,:,source_index,subsource_index),x_integral_subgrid,y_integral_subgrid)
            
        enddo
        
        return
        
    endif
    
    allocate (temp_FF_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index))) 
    allocate (angle_diff(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index))) 
    
    !Precalculate information for the trajectory model
    !Maxium number of trajectory steps and size of steps based on the integral (meteorology) loop size
    traj_step_size=min(integral_subgrid_delta(x_dim_index),integral_subgrid_delta(y_dim_index))*traj_step_scale
    traj_max_index=floor(max(integral_subgrid_loop_index(x_dim_index),integral_subgrid_loop_index(y_dim_index))/traj_step_scale)


    do subsource_index=1,n_subsource(source_index)

        call uEMEP_set_dispersion_params_simple(source_index,subsource_index)
    
        !Set local dispersion parameters prior to time loop for use in annual data
        ay_loc=ay(source_index,subsource_index)
        by_loc=by(source_index,subsource_index)
        az_loc=az(source_index,subsource_index)
        bz_loc=bz(source_index,subsource_index)
        sig_y_0_loc=sig_y_0(source_index,subsource_index)
        sig_z_0_loc=sig_z_0(source_index,subsource_index)
        h_emis_loc=h_emis(source_index,subsource_index)
        z_rec_loc=z_rec(source_index,subsource_index)

        !write(*,*) z_rec_loc,ay_loc,by_loc,az_loc,bz_loc,sig_y_0_loc,sig_z_0_loc
    
        write(unit_logfile,'(a,i3)')'Calculating proxy integral concentration data for '//trim(source_file_str(source_index))//' with subsource index ',subsource_index

        !Set the start and end times of the loop
        t_start=1
        t_end=integral_subgrid_dim(t_dim_index)
    
        !Loop through the time
        do tt=t_start,t_end
        
        integral_subgrid(:,:,tt,source_index,:)=0.
    
        !Create a temporary wind speed subgrid for each hour
        temp_FF_subgrid=0.
        do j_cross=1,integral_subgrid_dim(y_dim_index)
        do i_cross=1,integral_subgrid_dim(x_dim_index)
            z0_temp=exp(meteo_subgrid(i_cross,j_cross,tt,logz0_subgrid_index))
            h_temp=h_emis(source_index,subsource_index)
            if (annual_calculations.and.wind_level_integral_flag.eq.1) then
                temp_FF_subgrid(i_cross,j_cross)=1./meteo_subgrid(i_cross,j_cross,tt,inv_FFgrid_subgrid_index)
            elseif (annual_calculations.and.wind_level_integral_flag.eq.2) then
                temp_FF_subgrid(i_cross,j_cross)=1./meteo_subgrid(i_cross,j_cross,tt,inv_FFgrid_subgrid_index)*(1.-(log((H_emep/2.+z0_temp)/z0_temp)-log((h_temp+z0_temp)/z0_temp))/log((H_emep/2.+z0_temp)/z0_temp))
            elseif (annual_calculations.and.wind_level_integral_flag.eq.3) then
                temp_FF_subgrid(i_cross,j_cross)=1./meteo_subgrid(i_cross,j_cross,tt,inv_FF10_subgrid_index)
            elseif (annual_calculations.and.wind_level_integral_flag.eq.4) then
                temp_FF_subgrid(i_cross,j_cross)=1./meteo_subgrid(i_cross,j_cross,tt,inv_FF10_subgrid_index)*(1.-(log((10.+z0_temp)/z0_temp)-log((h_temp+z0_temp)/z0_temp))/log((10.+z0_temp)/z0_temp))
            elseif (hourly_calculations.and.wind_level_integral_flag.eq.1) then
                temp_FF_subgrid(i_cross,j_cross)=meteo_subgrid(i_cross,j_cross,tt,FFgrid_subgrid_index)
            elseif (hourly_calculations.and.wind_level_integral_flag.eq.2) then
                temp_FF_subgrid(i_cross,j_cross)=meteo_subgrid(i_cross,j_cross,tt,FFgrid_subgrid_index)*(1.-(log((H_emep/2.+z0_temp)/z0_temp)-log((h_temp+z0_temp)/z0_temp))/log((H_emep/2.+z0_temp)/z0_temp))
            elseif (hourly_calculations.and.wind_level_integral_flag.eq.3) then
                temp_FF_subgrid(i_cross,j_cross)=meteo_subgrid(i_cross,j_cross,tt,FF10_subgrid_index)
            elseif (hourly_calculations.and.wind_level_integral_flag.eq.4) then
                temp_FF_subgrid(i_cross,j_cross)=meteo_subgrid(i_cross,j_cross,tt,FF10_subgrid_index)*(1.-(log((10.+z0_temp)/z0_temp)-log((h_temp+z0_temp)/z0_temp))/log((10.+z0_temp)/z0_temp))
            elseif (wind_level_integral_flag.eq.0) then
                temp_FF_subgrid(i_cross,j_cross)=1.
            else
                write(unit_logfile,'(a)') 'No valid wind_level_integral_flag selected. Stopping (uEMEP_grid_proxy)'
                stop
            endif
            !Setting a minimum value for wind for dispersion purposes (cannot be zero)
            temp_FF_subgrid(i_cross,j_cross)=sqrt(temp_FF_subgrid(i_cross,j_cross)*temp_FF_subgrid(i_cross,j_cross)+FF_min_dispersion*FF_min_dispersion)

            if (use_last_meteo_in_dispersion) then
                cos_subgrid_loc=meteo_subgrid(i_cross,j_cross,tt,cos_subgrid_index)
                sin_subgrid_loc=meteo_subgrid(i_cross,j_cross,tt,sin_subgrid_index)
                angle_diff(i_cross,j_cross)=abs(asin(meteo_subgrid(i_cross,j_cross,tt,sin_subgrid_index)*last_meteo_subgrid(i_cross,j_cross,cos_subgrid_index) &
                                           -meteo_subgrid(i_cross,j_cross,tt,cos_subgrid_index)*last_meteo_subgrid(i_cross,j_cross,sin_subgrid_index) ))/2.
                
                angle_diff(i_cross,j_cross)=min(angle_diff(i_cross,j_cross),3.14159/4.) !Less than 45 degrees
                !write(*,*) i_cross,j_cross,angle_diff(i_cross,j_cross)*180/3.14159
                                        
            else
                angle_diff(i_cross,j_cross)=0.
            endif

        enddo
        enddo
                           
        !Loop through the proxy integral grid
        do j=1,integral_subgrid_dim(y_dim_index)
        do i=1,integral_subgrid_dim(x_dim_index)
                            

                if (EMEP_projection_type.eq.LL_projection_index) then
                    xpos_integral_subgrid=lon_integral_subgrid(i,j)
                    ypos_integral_subgrid=lat_integral_subgrid(i,j)
                elseif (EMEP_projection_type.eq.LCC_projection_index) then
                    call lb2lambert_uEMEP(xpos_integral_subgrid,ypos_integral_subgrid,lon_integral_subgrid(i,j),lat_integral_subgrid(i,j),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                endif

                !Find the cross reference to the integral grid from the emission grid
                i_cross=crossreference_integral_to_emission_subgrid(i,j,x_dim_index,source_index)
                j_cross=crossreference_integral_to_emission_subgrid(i,j,y_dim_index,source_index)

                if (use_downwind_position_flag.and.hourly_calculations) then
                    
                    i_cross_integral=i
                    j_cross_integral=j                   
                    
                    x_pos=max(-1.,min(1.,meteo_subgrid(i_cross_integral,j_cross_integral,tt,cos_subgrid_index)*sqrt(2.)))
                    y_pos=max(-1.,min(1.,meteo_subgrid(i_cross_integral,j_cross_integral,tt,sin_subgrid_index)*sqrt(2.)))                    
                    
                    i_end=min(ceiling(i_cross+(1.-x_pos)*emission_subgrid_loop_index(x_dim_index,source_index)),emission_subgrid_dim(x_dim_index,source_index))
                    i_start=max(floor(i_cross-(1.+x_pos)*emission_subgrid_loop_index(x_dim_index,source_index)),1)
                    j_end=min(ceiling(j_cross+(1.-y_pos)*emission_subgrid_loop_index(y_dim_index,source_index)),emission_subgrid_dim(y_dim_index,source_index))
                    j_start=max(floor(j_cross-(1.+y_pos)*emission_subgrid_loop_index(y_dim_index,source_index)),1)
                                                    
                    lon_area_max=xpos_integral_subgrid+(1.-x_pos)*lon_limit
                    lon_area_min=xpos_integral_subgrid-(1.+x_pos)*lon_limit
                    lat_area_max=ypos_integral_subgrid+(1.-y_pos)*lat_limit
                    lat_area_min=ypos_integral_subgrid-(1.+y_pos)*lat_limit

                else
                               
                    !Set the size of the loop region around the target cell to be up to integral_subgrid_loop_index
                    i_start=max(1,i_cross-emission_subgrid_loop_index(x_dim_index,source_index))
                    i_end=min(emission_subgrid_dim(x_dim_index,source_index),i_cross+emission_subgrid_loop_index(x_dim_index,source_index))
                    j_start=max(1,j_cross-emission_subgrid_loop_index(y_dim_index,source_index))
                    j_end=min(emission_subgrid_dim(y_dim_index,source_index),j_cross+emission_subgrid_loop_index(y_dim_index,source_index))
 

                    lon_area_max=xpos_integral_subgrid+lon_limit
                    lon_area_min=xpos_integral_subgrid-lon_limit
                    lat_area_max=ypos_integral_subgrid+lat_limit
                    lat_area_min=ypos_integral_subgrid-lat_limit

                endif
          
                !Loop through emission sub_grids in the nearby region
                do jj=j_start,j_end
                do ii=i_start,i_end
                            
                    if (emission_subgrid(ii,jj,tt,source_index,subsource_index).ne.0) then

                        if (EMEP_projection_type.eq.LL_projection_index) then
                            xpos_emission_subgrid=lon_emission_subgrid(ii,jj,source_index)
                            ypos_emission_subgrid=lat_emission_subgrid(ii,jj,source_index)
                        elseif (EMEP_projection_type.eq.LCC_projection_index) then
                            call lb2lambert_uEMEP(xpos_emission_subgrid,ypos_emission_subgrid,lon_emission_subgrid(ii,jj,source_index),lat_emission_subgrid(ii,jj,source_index),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                        endif

                        !Limit distance to the size of the EMEP grid
                        !if (abs(lon_emission_subgrid(ii,jj,source_index)-lon_integral_subgrid(i,j)).le.lon_limit &
                        !    .and.abs(lat_emission_subgrid(ii,jj,source_index)-lat_integral_subgrid(i,j)).le.lat_limit) then
                        !if (xpos_emission_subgrid.gt.lon_area_min.and.xpos_emission_subgrid.lt.lon_area_max &
                        !    .and.ypos_emission_subgrid.gt.lat_area_min.and.ypos_emission_subgrid.lt.lat_area_max) then                  
                        if (xpos_emission_subgrid.ge.lon_area_min.and.xpos_emission_subgrid.le.lon_area_max &
                            .and.ypos_emission_subgrid.ge.lat_area_min.and.ypos_emission_subgrid.le.lat_area_max) then                  
               
                            if (hourly_calculations) then
                            
                                !tt_emis=1   !Currently set to constant emissions
                            
                                !Determine meteorology at source position
                                i_cross_integral=crossreference_emission_to_integral_subgrid(ii,jj,x_dim_index,source_index)
                                j_cross_integral=crossreference_emission_to_integral_subgrid(ii,jj,y_dim_index,source_index)                          
                                kz_subgrid_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,kz_subgrid_index)
                                !u_subgrid_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,ugrid_subgrid_index)
                                !v_subgrid_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,vgrid_subgrid_index)
                                cos_subgrid_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,cos_subgrid_index)
                                sin_subgrid_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,sin_subgrid_index)
                            
                                !Calculate dispersion parameters based on meteo at source
                                if (stability_scheme_flag.eq.1) call uEMEP_set_dispersion_params_simple(source_index,subsource_index)
                                if (stability_scheme_flag.eq.2) call uEMEP_set_dispersion_params_PG(meteo_subgrid(i_cross_integral,j_cross_integral,tt,invL_subgrid_index),source_index,subsource_index)
    
                                !Set local dispersion parameters
                                ay_loc=ay(source_index,subsource_index)
                                by_loc=by(source_index,subsource_index)
                                az_loc=az(source_index,subsource_index)
                                bz_loc=bz(source_index,subsource_index)
                                sig_y_0_loc=sig_y_0(source_index,subsource_index)
                                sig_z_0_loc=sig_z_0(source_index,subsource_index)  
                            
                                if (use_trajectory_flag) then
                                
                                    !call uEMEP_local_trajectory(i,j,ii,jj,tt,traj_max_index,traj_step_size,source_index,x_loc,y_loc,valid_traj)
                                    call uEMEP_local_trajectory(x_integral_subgrid(i,j),y_integral_subgrid(i,j),x_emission_subgrid(ii,jj,source_index),y_emission_subgrid(ii,jj,source_index),tt,traj_max_index,traj_step_size,x_loc,y_loc,valid_traj)
                               
                                    if (valid_traj) then
                                        integral_subgrid(i,j,tt,source_index,subsource_index)=integral_subgrid(i,j,tt,source_index,subsource_index) &
                                            +emission_subgrid(ii,jj,tt,source_index,subsource_index) &
                                            *gauss_plume_cartesian_trajectory_integral_func(x_loc,y_loc,h_emis_loc,z_rec_loc,ay_loc,by_loc,az_loc,bz_loc,sig_y_0_loc,sig_z_0_loc,0.,H_emep,angle_diff(i_cross_integral,j_cross_integral)) &
                                            /temp_FF_subgrid(i_cross_integral,j_cross_integral)
                                    endif
                                
                                else
                                
                                    integral_subgrid(i,j,tt,source_index,subsource_index)=integral_subgrid(i,j,tt,source_index,subsource_index) &
                                        +emission_subgrid(ii,jj,tt,source_index,subsource_index) &
                                        *gauss_plume_cartesian_integral_func(x_emission_subgrid(ii,jj,source_index),y_emission_subgrid(ii,jj,source_index),h_emis_loc,cos_subgrid_loc,sin_subgrid_loc, &
                                        x_integral_subgrid(i,j),y_integral_subgrid(i,j),z_rec_loc,ay_loc,by_loc,az_loc,bz_loc,sig_y_0_loc,sig_z_0_loc,0.,H_emep,angle_diff(i_cross_integral,j_cross_integral)) &
                                        /temp_FF_subgrid(i_cross_integral,j_cross_integral)
                                                            
                                endif
                                
                            else
                            
                                distance_subgrid=sqrt((x_emission_subgrid(ii,jj,source_index)-x_integral_subgrid(i,j))*(x_emission_subgrid(ii,jj,source_index)-x_integral_subgrid(i,j)) &
                                    +(y_emission_subgrid(ii,jj,source_index)-y_integral_subgrid(i,j))*(y_emission_subgrid(ii,jj,source_index)-y_integral_subgrid(i,j)))
                            
                                !Divide by wind speed at emission position
                                i_cross_integral=crossreference_emission_to_integral_subgrid(ii,jj,x_dim_index,source_index)
                                j_cross_integral=crossreference_emission_to_integral_subgrid(ii,jj,y_dim_index,source_index)                          

                                integral_subgrid(i,j,tt,source_index,subsource_index)=integral_subgrid(i,j,tt,source_index,subsource_index) &
                                    +emission_subgrid(ii,jj,tt,source_index,subsource_index) &
                                    *gauss_plume_second_order_rotated_integral_func(distance_subgrid,z_rec_loc,ay_loc,by_loc,az_loc,bz_loc,sig_y_0_loc,sig_z_0_loc,h_emis_loc,0.,H_emep) &
                                    /temp_FF_subgrid(i_cross_integral,j_cross_integral)
        
                            endif
                        
                        endif
                    
                    endif
                    
                    enddo
                    enddo
            
        enddo
            !write(*,'(3a,i5,a,i5,a,i3,a,i3)') 'Gridding ',trim(source_file_str(source_index)),' integral proxy',j,' of ',integral_subgrid_dim(y_dim_index),' and ',subsource_index,' of ',n_subsource(source_index)
        enddo
        
        if (mod(j,1).eq.0) write(*,'(3a,i5,a,i5,a,i3,a,i3)') 'Integral gridding ',trim(source_file_str(source_index)),' proxy for hour ',tt,' of ',subgrid_dim(t_dim_index),' and ',subsource_index,' of ',n_subsource(source_index)
        enddo !time loop      

    enddo !subsource_index
               
        if (combine_emission_subsources_during_dispersion(source_index).and.n_subsource(source_index).gt.1) then
            do subsource_index=2,n_subsource(n_source_index)
                integral_subgrid(:,:,:,source_index,1)=integral_subgrid(:,:,:,source_index,1)+integral_subgrid(:,:,:,source_index,subsource_index)
            enddo
            n_subsource(source_index)=1
        endif
    
            
    if (save_intermediate_files) then
        if (.not.read_existing_grid_data(proxy_file_index(source_index))) then
            do subsource_index=1,n_subsource(n_source_index)
            temp_name=trim(pathname_grid(proxy_integral_file_index(source_index)))//trim(filename_grid(proxy_integral_file_index(source_index)))//trim(subsource_str(subsource_index))//'_'//trim(file_tag)//'.asc'
            write(unit_logfile,'(a)')'Writing to: '//trim(temp_name)
            call write_esri_ascii_3d_file(unit_logfile,temp_name,integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),integral_subgrid_dim(t_dim_index), &
                integral_subgrid_delta(1),integral_subgrid(:,:,:,source_index,subsource_index),x_integral_subgrid,y_integral_subgrid)
            enddo
        endif
    endif


    end subroutine uEMEP_grid_proxy_integral