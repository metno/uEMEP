!uEMEP_subgrid_dispersion.f90
    
!==========================================================================
!   uEMEP model save_gridded_proxy
!   Calculates proxy concentrations based on Gaussian distribution
!   Generic for all sources and subsources
!==========================================================================
    subroutine uEMEP_subgrid_dispersion(source_index)
    
    use uEMEP_definitions

    implicit none
    
    integer i,j,k
    integer         source_index
    character(256)  temp_name
    logical         exists
    integer         jj,ii,tt,tt_emis
    real            distance_subgrid
    real            dx_temp,dy_temp
    integer         i_start,i_end,j_start,j_end,t_start,t_end
    integer         i_cross,j_cross
    integer         i_cross_integral,j_cross_integral,i_cross_target_integral,j_cross_target_integral
    real            u_subgrid_loc,v_subgrid_loc,kz_subgrid_loc
    real            cos_subgrid_loc,sin_subgrid_loc,FF_loc,FF_zc_loc
    integer         i_nc,j_nc
    integer         subsource_index
    real            ay_loc,by_loc,az_loc,bz_loc,sig_y_0_loc,sig_z_0_loc,sig_y_00_loc,sig_z_00_loc,h_emis_loc,z_rec_loc,sig_z_loc,sig_y_loc,h_mix_loc
    real            xpos_limit,ypos_limit
    real            time_weight(subgrid_dim(t_dim_index),n_pollutant_loop),time_total(subgrid_dim(t_dim_index),n_pollutant_loop)
    integer         time_count(subgrid_dim(t_dim_index))
    real            x_downwind,y_downwind
    real            xpos_area_min,xpos_area_max,ypos_area_min,ypos_area_max
    real            distance_subgrid_min
    real            xpos_subgrid,ypos_subgrid
    real            xpos_emission_subgrid,ypos_emission_subgrid
    real            temp_subgrid_internal
    real            distance_emission_subgrid_min
    real            temp_sum_subgrid(n_pollutant_loop)
    integer         count
    
    real, allocatable :: temp_emission_subgrid(:,:,:)
    real, allocatable :: temp_subgrid(:,:,:)
    !real, allocatable :: diagnostic_subgrid(:,:,:)
    real, allocatable :: temp_FF_subgrid(:,:)
    real, allocatable :: temp_FF_emission_subgrid(:,:)
    real, allocatable :: trajectory_subgrid(:,:,:,:)
    real, allocatable :: angle_diff(:,:)

    integer traj_max_index
    logical valid_traj
    real traj_step_size,x_loc,y_loc,FFgrid_loc,logz0_loc,u_star0_loc,FF10_loc,zc_loc,invL_loc
    real z0_temp,h_temp
    
    real, allocatable :: temp_target_subgrid(:,:,:)
    real, allocatable :: x_target_subgrid(:,:)
    real, allocatable :: y_target_subgrid(:,:)
    real, allocatable :: traveltime_temp_target_subgrid(:,:,:,:)
    
    integer temp_target_subgrid_dim_min(2),temp_target_subgrid_dim_max(2)
    integer temp_target_subgrid_dim_length(2)
    real temp_target_subgrid_delta(2)
    logical :: use_target_subgrid=.true.
    integer i_target_start,i_target_end,j_target_start,j_target_end
    real area_weighted_interpolation_function
    logical temp_use_subgrid

    integer i_pollutant
    real temp_subgrid_internal_pollutant(n_pollutant_loop)
    
    !functions
    real gauss_plume_second_order_rotated_func
    !real gauss_plume_second_order_rotated_integral_func
    real gauss_plume_cartesian_func
    !real gauss_plume_cartesian_integral_func
    real gauss_plume_cartesian_trajectory_func
    real gauss_plume_cartesian_sigma_func
    real gauss_plume_second_order_rotated_reflected_func
    
    write(unit_logfile,'(A)') ''
    write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Calculating dispersion of proxy (uEMEP_subgrid_dispersion)'
	write(unit_logfile,'(A)') '================================================================'

    !First call the integral dispersion routine if it is needed. Only when using the concentration redistribution or the EMEP grid interpolation with proxy
    if (local_subgrid_method_flag.eq.1.or.EMEP_grid_interpolation_flag.eq.4) then
        call uEMEP_subgrid_dispersion_integral(source_index)
    endif
    
 
    allocate (temp_emission_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),n_pollutant_loop)) 
    allocate (temp_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),n_pollutant_loop)) 
    !allocate (diagnostic_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),n_pollutant_loop)) 
    allocate (temp_FF_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index))) 
    allocate (temp_FF_emission_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index))) 
    allocate (angle_diff(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index))) 
    
    temp_subgrid=0.
    temp_emission_subgrid=0.
    temp_FF_subgrid=0.
    temp_FF_emission_subgrid=0.
    angle_diff=0.
        
    !Set the x and y position limits to coincide to half the EMEP grid (refered here as lon and lat but can be also LCC projection) times the number of grids
    xpos_limit=dgrid_nc(lon_nc_index)/2.*EMEP_grid_interpolation_size
    ypos_limit=dgrid_nc(lat_nc_index)/2.*EMEP_grid_interpolation_size

    !Minimum distance for travel time calculation set to  half of a grid diagonal weighted so the circle has the same area as the square with that diagonal
    distance_subgrid_min=sqrt(subgrid_delta(x_dim_index)*subgrid_delta(x_dim_index)+subgrid_delta(y_dim_index)*subgrid_delta(y_dim_index))/2./sqrt(2.)*4./3.14159
    !Minimum distance for dispersion set to  half of an emission grid diagonal weighted so the circle has the same area as the square with that diagonal
    distance_emission_subgrid_min=sqrt(emission_subgrid_delta(x_dim_index,source_index)*emission_subgrid_delta(x_dim_index,source_index) &
        +emission_subgrid_delta(y_dim_index,source_index)*emission_subgrid_delta(y_dim_index,source_index))/2./sqrt(2.)*4./3.14159
    
    do subsource_index=1,n_subsource(source_index)

        !Do not use target subgrid if the grid is auto selected
        !if (emission_subgrid_delta(x_dim_index,source_index).le.subgrid_delta(x_dim_index).or. &
        if (use_emission_positions_for_auto_subgrid_flag(source_index)) then
            !If auto positions using emissions for that source then do not use the interpolation target grid
            use_target_subgrid=.false.
            write(unit_logfile,*) 'Using auto subgrid for source ',trim(source_file_str(source_index))
        elseif (emission_subgrid_delta(x_dim_index,source_index).le.subgrid_delta(x_dim_index)) then
            !If the subgrid emissions are less than or equal to the dispersion grid then do not use the target interpolation grid
            !No matter what auto subgrid is used
            use_target_subgrid=.false.
        else
            !Use the target subgrid even when the other auto subgrids are on. Slows it down but is necessary to get the right interpolation
            !(use_population_positions_for_auto_subgrid_flag.or.use_receptor_positions_for_auto_subgrid_flag)
            use_target_subgrid=.true.
        endif

        
        if (use_target_subgrid) then
            write(unit_logfile,*) 'Using emission subgrid with interpolation for source ',trim(source_file_str(source_index))
        else
            write(unit_logfile,*) 'Using normal subgrid with no interpolation for source ',trim(source_file_str(source_index))
        endif
        
        
        call uEMEP_set_dispersion_params_simple(source_index,subsource_index)
    
        !Set local dispersion parameters to be used only in the annual calculation, overwritten in the hourly files
        ay_loc=ay(source_index,subsource_index)
        by_loc=by(source_index,subsource_index)
        az_loc=az(source_index,subsource_index)
        bz_loc=bz(source_index,subsource_index)
        sig_y_00_loc=sig_y_00(source_index,subsource_index)
        sig_z_00_loc=sig_z_00(source_index,subsource_index)
        h_emis_loc=h_emis(source_index,subsource_index)
        z_rec_loc=z_rec(source_index,subsource_index)

    
        write(unit_logfile,'(a,i3)')'Calculating proxy concentration data for '//trim(source_file_str(source_index))//' with subsource index ',subsource_index

        !Set up a target grid that matches the emissions grid and is just slightly bigger than it
        !Find the grid index it belongs to
        if (use_target_subgrid) then
        temp_target_subgrid_dim_min(x_dim_index)=-1+1+floor((subgrid_min(x_dim_index)-emission_subgrid_min(x_dim_index,source_index))/emission_subgrid_delta(x_dim_index,source_index))
        temp_target_subgrid_dim_min(y_dim_index)=-1+1+floor((subgrid_min(y_dim_index)-emission_subgrid_min(y_dim_index,source_index))/emission_subgrid_delta(y_dim_index,source_index))
        temp_target_subgrid_dim_max(x_dim_index)=+1+1+ceiling((subgrid_max(x_dim_index)-emission_subgrid_min(x_dim_index,source_index))/emission_subgrid_delta(x_dim_index,source_index))
        temp_target_subgrid_dim_max(y_dim_index)=+1+1+ceiling((subgrid_max(y_dim_index)-emission_subgrid_min(y_dim_index,source_index))/emission_subgrid_delta(y_dim_index,source_index))

        temp_target_subgrid_dim_length(x_dim_index)=temp_target_subgrid_dim_max(x_dim_index)-temp_target_subgrid_dim_min(x_dim_index)+1
        temp_target_subgrid_dim_length(y_dim_index)=temp_target_subgrid_dim_max(y_dim_index)-temp_target_subgrid_dim_min(y_dim_index)+1
        temp_target_subgrid_delta(x_dim_index)=emission_subgrid_delta(x_dim_index,source_index)
        temp_target_subgrid_delta(y_dim_index)=emission_subgrid_delta(y_dim_index,source_index)

        !Reallocate internal target arrays for each source
        if (allocated(temp_target_subgrid)) deallocate (temp_target_subgrid)
        if (allocated(x_target_subgrid)) deallocate (x_target_subgrid)
        if (allocated(y_target_subgrid)) deallocate (y_target_subgrid)
        if (allocated(traveltime_temp_target_subgrid)) deallocate (traveltime_temp_target_subgrid)
        if (.not.allocated(temp_target_subgrid)) allocate (temp_target_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),n_pollutant_loop)) 
        if (.not.allocated(traveltime_temp_target_subgrid)) allocate (traveltime_temp_target_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),2,n_pollutant_loop)) 
        if (.not.allocated(x_target_subgrid)) allocate (x_target_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index))) 
        if (.not.allocated(y_target_subgrid)) allocate (y_target_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index))) 
        
        
        x_target_subgrid(:,:)=x_emission_subgrid(:,:,source_index)
        y_target_subgrid(:,:)=y_emission_subgrid(:,:,source_index)

        endif
        
        !Set the start and end times of the loop
        t_start=1
        t_end=subgrid_dim(t_dim_index)
        
        !Loop through the time
        do tt=t_start,t_end
        
        subgrid(:,:,tt,proxy_subgrid_index,source_index,:)=0.
        temp_target_subgrid=0.
        traveltime_temp_target_subgrid=0.
    
        !Set a temporary emission array
        temp_emission_subgrid=emission_subgrid(:,:,tt,source_index,:)
        temp_subgrid=0.
        temp_FF_subgrid=0.
        !diagnostic_subgrid=0.
        
        !Set the last meteo data subgrid in the case when the internal time loop is used
        if (.not.use_single_time_loop_flag) then
            if (tt.gt.t_start) then
                last_meteo_subgrid(:,:,:)=meteo_subgrid(:,:,tt-1,:)
            else
                last_meteo_subgrid(:,:,:)=meteo_subgrid(:,:,tt,:)      
            endif           
        endif

        !Precalculate information for the trajectory model
        !Maxium number of trajectory steps and size of steps based on the integral (meteorology) loop size
        if (use_trajectory_flag(source_index)) then
            
            traj_step_size=min(integral_subgrid_delta(x_dim_index),integral_subgrid_delta(y_dim_index))*traj_step_scale
            traj_max_index=floor(max(integral_subgrid_loop_index(x_dim_index),integral_subgrid_loop_index(y_dim_index))/traj_step_scale)

            if (tt.eq.t_start) write(unit_logfile,'(a,f12.1,i)') 'Trajectory step (m) and dimensions: ', traj_step_size, traj_max_index
            if (.not.allocated(trajectory_subgrid)) allocate(trajectory_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),traj_max_index,2))
            
            trajectory_subgrid=NODATA_value
            
            !Loop through the emissions and create trajectories for all emissions source grids
            do j=1,emission_subgrid_dim(y_dim_index,source_index)
            do i=1,emission_subgrid_dim(x_dim_index,source_index)
            
                if (sum(temp_emission_subgrid(i,j,:)).ne.0) then
                    call uEMEP_calculate_all_trajectory(x_emission_subgrid(i,j,source_index),y_emission_subgrid(i,j,source_index),tt, &
                        traj_max_index,traj_step_size,trajectory_subgrid(i,j,:,x_dim_index),trajectory_subgrid(i,j,:,y_dim_index))
                endif
            
            enddo
            enddo
       
        endif

        !Create a temporary wind speed subgrid for each hour
        temp_FF_subgrid=0.
        do j_cross=1,integral_subgrid_dim(y_dim_index)
        do i_cross=1,integral_subgrid_dim(x_dim_index)
            z0_temp=exp(meteo_subgrid(i_cross,j_cross,tt,logz0_subgrid_index))
            h_temp=h_emis(source_index,subsource_index)
            if (annual_calculations.and.wind_level_flag.eq.1) then
                temp_FF_subgrid(i_cross,j_cross)=1./meteo_subgrid(i_cross,j_cross,tt,inv_FFgrid_subgrid_index)
            elseif (annual_calculations.and.wind_level_flag.eq.2) then
                temp_FF_subgrid(i_cross,j_cross)=1./meteo_subgrid(i_cross,j_cross,tt,inv_FFgrid_subgrid_index)*(1.-(log((H_meteo+z0_temp)/z0_temp)-log((h_temp+z0_temp)/z0_temp))/log((H_meteo+z0_temp)/z0_temp))
            elseif (annual_calculations.and.wind_level_flag.eq.3) then
                temp_FF_subgrid(i_cross,j_cross)=1./meteo_subgrid(i_cross,j_cross,tt,inv_FF10_subgrid_index)
            elseif (annual_calculations.and.wind_level_flag.eq.4) then
                temp_FF_subgrid(i_cross,j_cross)=1./meteo_subgrid(i_cross,j_cross,tt,inv_FF10_subgrid_index)*(1.-(log((10.+z0_temp)/z0_temp)-log((h_temp+z0_temp)/z0_temp))/log((10.+z0_temp)/z0_temp))
            elseif (hourly_calculations.and.wind_level_flag.eq.1) then
                temp_FF_subgrid(i_cross,j_cross)=meteo_subgrid(i_cross,j_cross,tt,FFgrid_subgrid_index)
            elseif (hourly_calculations.and.wind_level_flag.eq.2) then
                temp_FF_subgrid(i_cross,j_cross)=meteo_subgrid(i_cross,j_cross,tt,FFgrid_subgrid_index)*(1.-(log((H_meteo+z0_temp)/z0_temp)-log((h_temp+z0_temp)/z0_temp))/log((H_meteo+z0_temp)/z0_temp))
            elseif (hourly_calculations.and.wind_level_flag.eq.3) then
                temp_FF_subgrid(i_cross,j_cross)=meteo_subgrid(i_cross,j_cross,tt,FF10_subgrid_index)
            elseif (hourly_calculations.and.wind_level_flag.eq.4) then
                temp_FF_subgrid(i_cross,j_cross)=meteo_subgrid(i_cross,j_cross,tt,FF10_subgrid_index)*(1.-(log((10.+z0_temp)/z0_temp)-log((h_temp+z0_temp)/z0_temp))/log((10.+z0_temp)/z0_temp))
            elseif (wind_level_flag.eq.0) then
                temp_FF_subgrid(i_cross,j_cross)=1.
            elseif (wind_level_flag.eq.5) then
                !Will set based on sigma z centre of mass
                temp_FF_subgrid(i_cross,j_cross)=1.
            elseif (wind_level_flag.eq.6) then
                !Will set based on sigma z centre of mass and emission height
                temp_FF_subgrid(i_cross,j_cross)=1.
            else
                
                write(unit_logfile,'(a)') 'No valid wind_level_flag selected. Stopping (uEMEP_subgrid_dispersion)'
                stop
            endif
            
            !Setting a minimum value for wind for dispersion purposes (cannot be zero)
            temp_FF_subgrid(i_cross,j_cross)=sqrt(temp_FF_subgrid(i_cross,j_cross)*temp_FF_subgrid(i_cross,j_cross)+FF_min_dispersion*FF_min_dispersion)
                    
            if (temp_FF_subgrid(i_cross,j_cross).eq.0) then
                write(unit_logfile,'(a,2i)') 'Zero wind speed at integral grid (stopping): ',i_cross,j_cross
                stop
            endif
            
            !Finds the angle difference between the current and last meteo field for dispersion and implements meandering if selected
            if (hourly_calculations) then
                call delta_wind_direction (i_cross,j_cross,tt,meteo_subgrid(i_cross,j_cross,tt,FF10_subgrid_index),angle_diff(i_cross,j_cross))
            else
                angle_diff(i_cross,j_cross)=0.
            endif
            

        enddo
        enddo

        !If wind level flag is set to 5, use of initial plume centre of mass, then set wind speed for each non-zero emission grid
        if (wind_level_flag.eq.5) then
            temp_FF_emission_subgrid=0.
            do jj=1,emission_subgrid_dim(y_dim_index,source_index)
            do ii=1,emission_subgrid_dim(x_dim_index,source_index)
            if (sum(temp_emission_subgrid(ii,jj,:)).ne.0) then
                                    
                !Set the integral meteorological grid position for the emission position
                i_cross_integral=crossreference_emission_to_integral_subgrid(ii,jj,x_dim_index,source_index)
                j_cross_integral=crossreference_emission_to_integral_subgrid(ii,jj,y_dim_index,source_index)
                
                !Set the local variables
                logz0_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,logz0_subgrid_index)
                FF10_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,FF10_subgrid_index)
                sig_y_00_loc=emission_properties_subgrid(ii,jj,emission_sigy00_index,source_index)
                sig_z_00_loc=emission_properties_subgrid(ii,jj,emission_sigz00_index,source_index)
                h_emis_loc=emission_properties_subgrid(ii,jj,emission_h_index,source_index)
                h_mix_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,hmix_subgrid_index)
                
                if (annual_calculations) then
                    FF10_loc=1./meteo_subgrid(i_cross_integral,j_cross_integral,tt,inv_FF10_subgrid_index)
                endif
                
                !Set sig_0's at the emission position
                x_loc=0.
                call uEMEP_set_dispersion_sigma_simple(sig_z_00_loc,sig_y_00_loc,sigy_0_subgid_width_scale,emission_subgrid_delta(:,source_index),angle_diff(i_cross_integral,j_cross_integral),x_loc,sig_z_loc,sig_y_loc,sig_z_0_loc,sig_y_0_loc)                                        

                !Use the initial plume centre of mass to determine wind advection height
                call z_centremass_gauss_func(sig_z_0_loc,h_emis_loc,h_mix_loc,zc_loc)
                call u_profile_neutral_val_func(zc_loc,FF10_loc,10.,h_mix_loc,exp(logz0_loc),FF_loc,u_star0_loc)
                
                !Set a minimum wind speed based on traffic (if use_traffic_for_minFF_flag=T)
!                FF_loc=sqrt(FF_loc*FF_loc+emission_properties_subgrid(ii,jj,emission_minFF_index,source_index,subsource_index)*emission_properties_subgrid(ii,jj,emission_minFF_index,source_index,subsource_index))
                
                !Set the minimum wind speed 
                FF_loc=sqrt(FF_loc*FF_loc+FF_min_dispersion*FF_min_dispersion)

                temp_FF_emission_subgrid(ii,jj)=FF_loc
                !write(*,*) FF10_loc,FF_loc,zc_loc,sig_z_0_loc
                
            
            endif
            enddo
            enddo
        endif

        !Loop through the target grid
        if (use_target_subgrid) then
            j_target_start=temp_target_subgrid_dim_min(y_dim_index);j_target_end=temp_target_subgrid_dim_max(y_dim_index)
            i_target_start=temp_target_subgrid_dim_min(x_dim_index);i_target_end=temp_target_subgrid_dim_max(x_dim_index)
        else
            j_target_start=1;j_target_end=subgrid_dim(y_dim_index)
            i_target_start=1;i_target_end=subgrid_dim(x_dim_index)
        endif
        !write(*,*) i_target_start,i_target_end,j_target_start,j_target_end
        
        do j=j_target_start,j_target_end
        do i=i_target_start,i_target_end
        
        !do j=1,subgrid_dim(y_dim_index)
        !do i=1,subgrid_dim(x_dim_index)
               
            !Only use those that are marked for use
            if (use_target_subgrid) then
                !Always use the grids because they cannot be tested
                temp_use_subgrid=.true.
            else
                temp_use_subgrid=use_subgrid(i,j,source_index)
            endif
            
            if (temp_use_subgrid) then    
             
                !Set the position of the target grid in terms of the EMEP projection
                if (use_target_subgrid) then
                    xpos_subgrid=xproj_emission_subgrid(i,j,source_index)
                    ypos_subgrid=yproj_emission_subgrid(i,j,source_index)
                else
                    xpos_subgrid=xproj_subgrid(i,j)
                    ypos_subgrid=yproj_subgrid(i,j)
                endif
                
                !Find the cross reference to the emission grid from the target grid
                if (use_target_subgrid) then
                    i_cross=i
                    j_cross=j
                else
                    i_cross=crossreference_target_to_emission_subgrid(i,j,x_dim_index,source_index)
                    j_cross=crossreference_target_to_emission_subgrid(i,j,y_dim_index,source_index)
                endif

                !Find the cross reference for the meteo grid at the target grid
                if (use_target_subgrid) then
                    i_cross_target_integral=crossreference_emission_to_integral_subgrid(i,j,x_dim_index,source_index)
                    j_cross_target_integral=crossreference_emission_to_integral_subgrid(i,j,y_dim_index,source_index)
                else
                    i_cross_target_integral=crossreference_target_to_integral_subgrid(i,j,x_dim_index)
                    j_cross_target_integral=crossreference_target_to_integral_subgrid(i,j,y_dim_index)
                endif
                

                !Set the travel time integral values to 0
                time_weight(tt,:)=0.
                time_total(tt,:)=0.
                
                !Use the wind direction to move the target area downwind. To reduce the search loop
                if (use_downwind_position_flag.and.hourly_calculations) then
                                    
                    !Set the emission grid loop region based on the downwind position
                    x_downwind=max(-1.,min(1.,meteo_subgrid(i_cross_target_integral,j_cross_target_integral,tt,cos_subgrid_index)*sqrt(2.)))
                    y_downwind=max(-1.,min(1.,meteo_subgrid(i_cross_target_integral,j_cross_target_integral,tt,sin_subgrid_index)*sqrt(2.)))
                    i_end=min(ceiling(i_cross+1+(1.-x_downwind)*emission_subgrid_loop_index(x_dim_index,source_index)),emission_subgrid_dim(x_dim_index,source_index))
                    i_start=max(floor(i_cross-1-(1.+x_downwind)*emission_subgrid_loop_index(x_dim_index,source_index)),1)
                    j_end=min(ceiling(j_cross+1+(1.-y_downwind)*emission_subgrid_loop_index(y_dim_index,source_index)),emission_subgrid_dim(y_dim_index,source_index))
                    j_start=max(floor(j_cross-1-(1.+y_downwind)*emission_subgrid_loop_index(y_dim_index,source_index)),1)
                                                    
                    !Set the EMEP projection limits to include the upwind source region
                    xpos_area_max=xpos_subgrid+(1.-x_downwind)*xpos_limit/2.+emission_subgrid_dim(x_dim_index,source_index)
                    xpos_area_min=xpos_subgrid-(1.+x_downwind)*xpos_limit/2.-emission_subgrid_dim(x_dim_index,source_index)
                    ypos_area_max=ypos_subgrid+(1.-y_downwind)*ypos_limit/2.+emission_subgrid_dim(y_dim_index,source_index)
                    ypos_area_min=ypos_subgrid-(1.+y_downwind)*ypos_limit/2.-emission_subgrid_dim(y_dim_index,source_index)

                else
                    
                    !Set the size of the loop region around the target cell to be up to subgrid_loop_index
                    i_start=max(1,i_cross-emission_subgrid_loop_index(x_dim_index,source_index))
                    i_end=min(emission_subgrid_dim(x_dim_index,source_index),i_cross+emission_subgrid_loop_index(x_dim_index,source_index))
                    j_start=max(1,j_cross-emission_subgrid_loop_index(y_dim_index,source_index))
                    j_end=min(emission_subgrid_dim(y_dim_index,source_index),j_cross+emission_subgrid_loop_index(y_dim_index,source_index))
 
                    !Set the emission limits (EMEP projection ) surrounding the target grid
                    xpos_area_max=xpos_subgrid+xpos_limit
                    xpos_area_min=xpos_subgrid-xpos_limit
                    ypos_area_max=ypos_subgrid+ypos_limit
                    ypos_area_min=ypos_subgrid-ypos_limit

                endif
    
                !Loop through emission sub_grids in the nearby region
                do jj=j_start,j_end
                do ii=i_start,i_end
                    
                    !Only non zero emissions to be calculated
                    if (sum(temp_emission_subgrid(ii,jj,:)).ne.0) then

                        !Set the EMEP projection position of the emission grid
                        xpos_emission_subgrid=xproj_emission_subgrid(ii,jj,source_index)
                        ypos_emission_subgrid=yproj_emission_subgrid(ii,jj,source_index)

                        !Select only emissions within the predefined region
                        if (xpos_emission_subgrid.ge.xpos_area_min.and.xpos_emission_subgrid.le.xpos_area_max &
                            .and.ypos_emission_subgrid.ge.ypos_area_min.and.ypos_emission_subgrid.le.ypos_area_max) then                  
                       
                            !Set the integral meteorological grid position for the emission position
                            i_cross_integral=crossreference_emission_to_integral_subgrid(ii,jj,x_dim_index,source_index)
                            j_cross_integral=crossreference_emission_to_integral_subgrid(ii,jj,y_dim_index,source_index)
                            i_cross_integral=min(max(1,i_cross_integral),integral_subgrid_dim(x_dim_index))
                            j_cross_integral=min(max(1,j_cross_integral),integral_subgrid_dim(y_dim_index))
                            
                            if (hourly_calculations) then 
                            
                                if (use_trajectory_flag(source_index)) then
                                
                                    !Calculate the minimum distance to the trajectory. Time consuming
                                    if (use_target_subgrid) then
                                    call uEMEP_minimum_distance_trajectory_fast(x_target_subgrid(i,j),y_target_subgrid(i,j),x_emission_subgrid(ii,jj,source_index),y_emission_subgrid(ii,jj,source_index), &
                                    traj_max_index,traj_step_size,trajectory_subgrid(ii,jj,:,x_dim_index),trajectory_subgrid(ii,jj,:,y_dim_index),x_loc,y_loc,valid_traj)
                                    else
                                    call uEMEP_minimum_distance_trajectory_fast(x_subgrid(i,j),y_subgrid(i,j),x_emission_subgrid(ii,jj,source_index),y_emission_subgrid(ii,jj,source_index), &
                                    traj_max_index,traj_step_size,trajectory_subgrid(ii,jj,:,x_dim_index),trajectory_subgrid(ii,jj,:,y_dim_index),x_loc,y_loc,valid_traj)
                                    endif
                                                                 
                                else
                                                                
                                    !Set the local wind cos and sin values
                                    cos_subgrid_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,cos_subgrid_index)
                                    sin_subgrid_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,sin_subgrid_index)                                    

                                    !Determine the rotated position along wind values
                                    if (use_target_subgrid) then
                                    x_loc=(x_target_subgrid(i,j)-x_emission_subgrid(ii,jj,source_index))*cos_subgrid_loc+(y_target_subgrid(i,j)-y_emission_subgrid(ii,jj,source_index))*sin_subgrid_loc
                                    y_loc=-(x_target_subgrid(i,j)-x_emission_subgrid(ii,jj,source_index))*sin_subgrid_loc+(y_target_subgrid(i,j)-y_emission_subgrid(ii,jj,source_index))*cos_subgrid_loc
                                    else
                                    x_loc=(x_subgrid(i,j)-x_emission_subgrid(ii,jj,source_index))*cos_subgrid_loc+(y_subgrid(i,j)-y_emission_subgrid(ii,jj,source_index))*sin_subgrid_loc
                                    y_loc=-(x_subgrid(i,j)-x_emission_subgrid(ii,jj,source_index))*sin_subgrid_loc+(y_subgrid(i,j)-y_emission_subgrid(ii,jj,source_index))*cos_subgrid_loc
                                    endif
                                    !write(*,*) x_loc,x_subgrid(i,j)-x_emission_subgrid(ii,jj,source_index),y_subgrid(i,j)-y_emission_subgrid(ii,jj,source_index)
                                                                       
                                    !If x is downwind then it is valid
                                    if (x_loc.ge.0.) then
                                        valid_traj=.true.
                                    else
                                        valid_traj=.false.
                                    endif
                                                            
                                endif
                                
                                !Calculate dispersion
                                if (valid_traj) then
                                    
                                    !Set the mixing height at the average of the emission and target position
                                    h_mix_loc=(meteo_subgrid(i_cross_integral,j_cross_integral,tt,hmix_subgrid_index)+meteo_subgrid(i_cross_target_integral,j_cross_target_integral,tt,hmix_subgrid_index))/2.
                                    !Set the local wind speed and other parameters at emission position
                                    FF_loc=temp_FF_subgrid(i_cross_integral,j_cross_integral)
                                    !L_loc=1./meteo_subgrid(i_cross_integral,j_cross_integral,tt,invL_subgrid_index)
                                    invL_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,invL_subgrid_index)
                                    FFgrid_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,FFgrid_subgrid_index)
                                    logz0_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,logz0_subgrid_index)
                                    u_star0_loc=max(meteo_subgrid(i_cross_integral,j_cross_integral,tt,ustar_subgrid_index),0.001)
                                    FF10_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,FF10_subgrid_index)
                                    sig_y_00_loc=emission_properties_subgrid(ii,jj,emission_sigy00_index,source_index)
                                    sig_z_00_loc=emission_properties_subgrid(ii,jj,emission_sigz00_index,source_index)
                                    h_emis_loc=emission_properties_subgrid(ii,jj,emission_h_index,source_index)
                                    
                                    if (wind_level_flag.eq.5.or.wind_level_flag.eq.6) then
                                        FF_loc=temp_FF_emission_subgrid(ii,jj)
                                    endif

                                    !Select method for assigning sigma
                                    if (stability_scheme_flag.eq.1) then
                                        call uEMEP_set_dispersion_sigma_simple(sig_z_00_loc,sig_y_00_loc,sigy_0_subgid_width_scale,emission_subgrid_delta(:,source_index),angle_diff(i_cross_integral,j_cross_integral),x_loc,sig_z_loc,sig_y_loc,sig_z_0_loc,sig_y_0_loc)                                        
                                    endif
                                    
                                    if (stability_scheme_flag.eq.2) then
                                        call uEMEP_set_dispersion_sigma_PG(invL_loc,logz0_loc,sig_z_00_loc,sig_y_00_loc,sigy_0_subgid_width_scale,emission_subgrid_delta(:,source_index),angle_diff(i_cross_integral,j_cross_integral),x_loc,sig_z_loc,sig_y_loc,sig_z_0_loc,sig_y_0_loc)
                                    endif
                                    
                                    if (stability_scheme_flag.eq.3) then
                                        !Set initial values for sigma. Initial sig_y is set here as well but is overridden by Kz dispersion
                                        call uEMEP_set_dispersion_sigma_simple(sig_z_00_loc,sig_y_00_loc,sigy_0_subgid_width_scale,emission_subgrid_delta(:,source_index),angle_diff(i_cross_integral,j_cross_integral),x_loc,sig_z_loc,sig_y_loc,sig_z_0_loc,sig_y_0_loc)
                                        !call uEMEP_set_dispersion_sigma_Kz_emulator(h_emis_loc,invL_loc,logz0_loc,h_mix_loc,sig_z_00_loc,sig_y_00_loc,emission_subgrid_delta(:,source_index),0.,x_loc,sig_z_loc,sig_y_loc,sig_z_0_loc,sig_y_0_loc)
                                        !write(*,*) 'H0: ',sig_z_loc
                                    
                                        call uEMEP_set_dispersion_sigma_Kz(x_loc,sig_z_00_loc,sig_y_00_loc,sigy_0_subgid_width_scale,sig_z_loc,h_emis_loc,h_mix_loc,invL_loc,FF10_loc,10.,logz0_loc,emission_subgrid_delta(:,source_index),u_star0_loc,average_zc_h_in_Kz_flag,n_kz_iterations,sig_z_loc,sig_y_loc,FF_zc_loc)
                                        !write(*,*) 'H1: ',sig_z_loc
                                        
                                        !Diagnostic
                                        if (2.eq.1.and.x_loc.lt.10.) then
                                        write(*,'(a,5f12.2,f12.4)') 'OUT KZ FULL: ',x_loc,sig_z_loc,sig_y_loc,h_mix_loc,exp(logz0_loc),invL_loc
                                        !call uEMEP_set_dispersion_sigma_Kz_emulator(h_emis_loc,invL_loc,logz0_loc,h_mix_loc,sig_z_00_loc,sig_y_00_loc,emission_subgrid_delta(:,source_index),angle_diff(i_cross_integral,j_cross_integral),x_loc,sig_z_loc,sig_y_loc,sig_z_0_loc,sig_y_0_loc)
                                        call uEMEP_set_dispersion_sigma_Kz_emulator(h_emis_loc,invL_loc,logz0_loc,h_mix_loc,sig_z_00_loc,sig_y_00_loc,sigy_0_subgid_width_scale,emission_subgrid_delta(:,source_index),0.,x_loc,sig_z_loc,sig_y_loc,sig_z_0_loc,sig_y_0_loc)
                                        write(*,'(A,5f12.2,f12.4)') 'OUT KZ EMUL: ',x_loc,sig_z_loc,sig_y_loc,h_mix_loc,exp(logz0_loc),invL_loc
                                        write(*,*)
                                        endif
                                        
                                        !Add the meandering and change in wind angle to the plume since not included in Kz calculation
                                        sig_y_loc=sig_y_loc+x_loc*angle_diff(i_cross_integral,j_cross_integral)
                                        
                                        !Use the average of the emision height and zc to determine wind speed. Is set to true if wind_level_flag=6
                                        if (wind_level_flag.eq.6) then
                                            !FF_loc=FF_zc_loc
                                            !Set the minimum wind speed 
                                            FF_loc=sqrt(FF_zc_loc*FF_zc_loc+FF_min_dispersion*FF_min_dispersion)
                                        endif

                                    endif
                                    
                                    if (stability_scheme_flag.eq.4) then
                                        call uEMEP_set_dispersion_sigma_Kz_emulator(h_emis_loc,invL_loc,logz0_loc,h_mix_loc,sig_z_00_loc,sig_y_00_loc,sigy_0_subgid_width_scale,emission_subgrid_delta(:,source_index),angle_diff(i_cross_integral,j_cross_integral),x_loc,sig_z_loc,sig_y_loc,sig_z_0_loc,sig_y_0_loc)
                                    endif
                                    
                                    !Adjust the height of the wind to the average of the emission and plume centre of mass height.
                                    !This is already the case in the Kz calculation so not repeated here.
                                    if (wind_level_flag.eq.6.and.stability_scheme_flag.ne.3) then
                                    !if (wind_level_flag.eq.6) then
                                        call z_centremass_gauss_func(sig_z_loc,h_emis_loc,h_mix_loc,zc_loc)
                                        zc_loc=(h_emis_loc+zc_loc)/2.
                                        call u_profile_neutral_val_func(zc_loc,FF10_loc,10.,h_mix_loc,exp(logz0_loc),FF_zc_loc,u_star0_loc)
                                        FF_loc=sqrt(FF_zc_loc*FF_zc_loc+FF_min_dispersion*FF_min_dispersion)
                                    endif
                                                                       
                                    !Calculate the dispersion
                                    temp_subgrid_internal=gauss_plume_cartesian_sigma_func(x_loc,y_loc,h_emis_loc,z_rec_loc,sig_z_loc,sig_y_loc,h_mix_loc,FF_loc)
                                    
                                    !if (ii.eq.i_cross.and.jj.eq.j_cross) write(*,'(9es12.2)') x_loc,y_loc,h_emis_loc,z_rec_loc,sig_z_loc,sig_y_loc,h_mix_loc,FF_loc,temp_subgrid_internal
                                    !if (tt.ge.18.and.tt.le.18.and.temp_subgrid_internal.gt.1.e-3) then
                                        !write(*,'(2i,12es12.2)') ii,jj,x_loc,y_loc,h_emis_loc,z_rec_loc,sig_z_loc,sig_y_loc,h_mix_loc,FF_loc,temp_emission_subgrid(ii,jj,1),temp_subgrid_internal,sin_subgrid_loc,cos_subgrid_loc
                                    !endif
                                    
                                                                        
                                    !For diagnostics only!!
                                    !diagnostic_subgrid(i,j,1)=diagnostic_subgrid(i,j,1)+temp_subgrid_internal
                                    !diagnostic_subgrid(i,j,2)=diagnostic_subgrid(i,j,2)+temp_emission_subgrid(ii,jj,1)
                                    
                                    do i_pollutant=1,n_pollutant_loop
                                    !Multiply by the emission factor
                                    temp_subgrid_internal_pollutant(i_pollutant)=temp_subgrid_internal*temp_emission_subgrid(ii,jj,i_pollutant)
                                    
                                    
                                    !Add to the receptor subgrid position
                                    if (use_target_subgrid) then
                                    temp_target_subgrid(i,j,i_pollutant)=temp_target_subgrid(i,j,i_pollutant)+temp_subgrid_internal_pollutant(i_pollutant)
                                    else
                                    temp_subgrid(i,j,i_pollutant)=temp_subgrid(i,j,i_pollutant)+temp_subgrid_internal_pollutant(i_pollutant)
                                    endif
                                    enddo
                                    
                                    !Determine the distance for the travel time calculation
                                    distance_subgrid=sqrt(x_loc*x_loc+y_loc*y_loc)
                                    
                                else
                                    
                                    temp_subgrid_internal=0.
                                    temp_subgrid_internal_pollutant=0.
                                    
                                endif
                                
                            
                                !Calculate weighted time based on the selected temp_FF_subgrid wind level
                                if (temp_subgrid_internal.gt.0) then
                                        
                                    distance_subgrid=max(distance_subgrid,distance_subgrid_min)

                                    !Take weighted average (weighted by concentration) of the time
                                    time_weight(tt,:)=time_weight(tt,:)+distance_subgrid/FF_loc*temp_subgrid_internal_pollutant
                                    
                                    !Calculate sum of the concentration for normalisation
                                    time_total(tt,:)=time_total(tt,:)+temp_subgrid_internal_pollutant
                                    
                                endif

                            else
                                                            
                                !If not hourly concentration then use the annual dispersion function
                                if (use_target_subgrid) then
                                    distance_subgrid=sqrt((x_emission_subgrid(ii,jj,source_index)-x_target_subgrid(i,j))*(x_emission_subgrid(ii,jj,source_index)-x_target_subgrid(i,j)) &
                                    +(y_emission_subgrid(ii,jj,source_index)-y_target_subgrid(i,j))*(y_emission_subgrid(ii,jj,source_index)-y_target_subgrid(i,j)))
                                else
                                    distance_subgrid=sqrt((x_emission_subgrid(ii,jj,source_index)-x_subgrid(i,j))*(x_emission_subgrid(ii,jj,source_index)-x_subgrid(i,j)) &
                                    +(y_emission_subgrid(ii,jj,source_index)-y_subgrid(i,j))*(y_emission_subgrid(ii,jj,source_index)-y_subgrid(i,j)))
                                endif
                                
                                !Set the simple as default
                                call uEMEP_set_dispersion_sigma_simple(sig_z_00_loc,sig_y_00_loc,sigy_0_subgid_width_scale,emission_subgrid_delta(:,source_index),angle_diff(i_cross_integral,j_cross_integral),x_loc,sig_z_loc,sig_y_loc,sig_z_0_loc,sig_y_0_loc)                                        

                                !Select method for assigning sigma
                                if (stability_scheme_flag.eq.1) then
                                    call uEMEP_set_dispersion_sigma_simple(sig_z_00_loc,sig_y_00_loc,sigy_0_subgid_width_scale,emission_subgrid_delta(:,source_index),angle_diff(i_cross_integral,j_cross_integral),x_loc,sig_z_loc,sig_y_loc,sig_z_0_loc,sig_y_0_loc)                                        
                                endif
                                    
                                if (stability_scheme_flag.eq.2) then
                                    call uEMEP_set_dispersion_sigma_PG(invL_loc,logz0_loc,sig_z_00_loc,sig_y_00_loc,sigy_0_subgid_width_scale,emission_subgrid_delta(:,source_index),angle_diff(i_cross_integral,j_cross_integral),x_loc,sig_z_loc,sig_y_loc,sig_z_0_loc,sig_y_0_loc)
                                endif

                                if (wind_level_flag.eq.5) then
                                    FF_loc=temp_FF_emission_subgrid(ii,jj)
                                else
                                    FF_loc=temp_FF_subgrid(i_cross_integral,j_cross_integral)
                                endif
                                
                                h_mix_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,hmix_subgrid_index)
                                invL_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,invL_subgrid_index)
                                logz0_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,logz0_subgrid_index)

                                
                                if (wind_level_flag.eq.6.and.stability_scheme_flag.ne.3) then
                                        FF10_loc=1./meteo_subgrid(i_cross_integral,j_cross_integral,tt,inv_FF10_subgrid_index)
                                        call z_centremass_gauss_func(sig_z_loc,h_emis_loc,h_mix_loc,zc_loc)
                                        zc_loc=(h_emis_loc+zc_loc)/2.
                                        call u_profile_neutral_val_func(zc_loc,FF10_loc,10.,h_mix_loc,exp(logz0_loc),FF_zc_loc,u_star0_loc)
                                        FF_loc=sqrt(FF_zc_loc*FF_zc_loc+FF_min_dispersion*FF_min_dispersion)
                                        !write(*,'(2i,7es12.2)') ii,jj,zc_loc,FF10_loc,10.,h_mix_loc,exp(logz0_loc),FF_zc_loc,u_star0_loc
                                endif
                                    
                                !write(*,*) h_mix_loc,invL_loc,logz0_loc,FF_loc
                          
                                !write(*,'(9es12.2)') distance_subgrid,z_rec_loc,ay_loc,by_loc,az_loc,bz_loc,sig_y_00_loc,sig_z_00_loc,h_emis_loc
                                !Divide by wind speed at emission position
                                !Changed from 00 to 0 in the sigmas
                                if (use_target_subgrid) then
                                temp_target_subgrid(i,j,:)=temp_target_subgrid(i,j,:) &
                                    +temp_emission_subgrid(ii,jj,:) &
                                    *gauss_plume_second_order_rotated_reflected_func(distance_subgrid,z_rec_loc,ay_loc,by_loc,az_loc,bz_loc,sig_y_0_loc,sig_z_0_loc,h_emis_loc,h_mix_loc) &
!                                    *gauss_plume_second_order_rotated_func(distance_subgrid,z_rec_loc,ay_loc,by_loc,az_loc,bz_loc,sig_y_0_loc,sig_z_0_loc,h_emis_loc) &
                                    /FF_loc
                                else
                                temp_subgrid(i,j,:)=temp_subgrid(i,j,:) &
                                    +temp_emission_subgrid(ii,jj,:) &
                                    *gauss_plume_second_order_rotated_reflected_func(distance_subgrid,z_rec_loc,ay_loc,by_loc,az_loc,bz_loc,sig_y_0_loc,sig_z_0_loc,h_emis_loc,h_mix_loc) &
!                                    *gauss_plume_second_order_rotated_func(distance_subgrid,z_rec_loc,ay_loc,by_loc,az_loc,bz_loc,sig_y_0_loc,sig_z_0_loc,h_emis_loc) &
                                    /FF_loc
                                endif
                                !write(*,'(4i5,2es12.2,4f12.3)') i,j,ii,jj,temp_subgrid(i,j,:), &
                                !    gauss_plume_second_order_rotated_func(distance_subgrid,z_rec_loc,ay_loc,by_loc,az_loc,bz_loc,sig_y_00_loc,sig_z_00_loc,h_emis_loc)/FF_loc &
                                !    ,distance_subgrid,az_loc,bz_loc,sig_z_00_loc
                            
                            endif
                        
                        endif
                    
                    endif
                    
                    
                enddo
                enddo
                                
                !Add to the travel time array
                if (use_target_subgrid) then
                    traveltime_temp_target_subgrid(i,j,1,:)=traveltime_temp_target_subgrid(i,j,1,:)+time_weight(tt,:)
                    traveltime_temp_target_subgrid(i,j,2,:)=traveltime_temp_target_subgrid(i,j,2,:)+time_total(tt,:)
                else
                    traveltime_subgrid(i,j,tt,1,:)=traveltime_subgrid(i,j,tt,1,:)+time_weight(tt,:)
                    traveltime_subgrid(i,j,tt,2,:)=traveltime_subgrid(i,j,tt,2,:)+time_total(tt,:)
                endif
                
            
            else
                !Set to nodata value for grids that should not be used for all pollutants
                temp_subgrid(i,j,:)=NODATA_value
            endif
            if (.not.use_target_subgrid) then             
                !write(*,'(3i,3es12.2)') tt,i,j,temp_subgrid(i,j,pollutant_loop_index(nox_index)),diagnostic_subgrid(i,j,1),diagnostic_subgrid(i,j,2)
            else
                !write(*,'(3i,3es12.2)') tt,i,j,temp_target_subgrid(i,j,pollutant_loop_index(nox_index)),diagnostic_subgrid(i,j,1),diagnostic_subgrid(i,j,2)
            endif
            
        enddo
        !if (mod(j,10).eq.0) write(*,'(3a,i5,a,i5,a,i3,a,i3)') 'Gridding ',trim(source_file_str(source_index)),' proxy',j,' of ',subgrid_dim(2),' and ',subsource_index,' of ',n_subsource(source_index)
        enddo
      
        if (mod(j,1).eq.0) write(*,'(3a,i5,a,i5,a,i3,a,i3)') 'Gridding ',trim(source_file_str(source_index)),' proxy for hour ',tt,' of ',subgrid_dim(t_dim_index),' and subsource ',subsource_index,' of ',n_subsource(source_index)
        
        !Put the temporary subgrid back into the subgrid array only for the selected grids
        if (use_target_subgrid) then
            !write(*,*) 'Mean temp traveltime target grid',tt,sum(traveltime_temp_target_subgrid(i_target_start:i_target_end,j_target_start:j_target_end,1))/temp_target_subgrid_dim_length(x_dim_index)/temp_target_subgrid_dim_length(y_dim_index)
            !write(*,*) 'Mean temp target grid',tt,sum(temp_target_subgrid(i_target_start:i_target_end,j_target_start:j_target_end,n_target_comp))/temp_target_subgrid_dim_length(x_dim_index)/temp_target_subgrid_dim_length(y_dim_index)
            !write(*,*) shape(traveltime_temp_target_subgrid),shape(traveltime_subgrid)
            do j=1,subgrid_dim(y_dim_index)
            do i=1,subgrid_dim(x_dim_index)
            if (use_subgrid(i,j,source_index)) then
                do i_pollutant=1,n_pollutant_loop
                temp_subgrid(i,j,i_pollutant)=area_weighted_interpolation_function(x_target_subgrid,y_target_subgrid,temp_target_subgrid(:,:,i_pollutant) &
                    ,emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),emission_subgrid_delta(:,source_index),x_subgrid(i,j),y_subgrid(i,j))
                traveltime_subgrid(i,j,tt,1,i_pollutant)=traveltime_subgrid(i,j,tt,1,i_pollutant) &
                    +area_weighted_interpolation_function(x_target_subgrid,y_target_subgrid,traveltime_temp_target_subgrid(:,:,1,i_pollutant) &
                    ,emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),emission_subgrid_delta(:,source_index),x_subgrid(i,j),y_subgrid(i,j))
                traveltime_subgrid(i,j,tt,2,i_pollutant)=traveltime_subgrid(i,j,tt,2,i_pollutant) &
                    +area_weighted_interpolation_function(x_target_subgrid,y_target_subgrid,traveltime_temp_target_subgrid(:,:,2,i_pollutant) &
                    ,emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),emission_subgrid_delta(:,source_index),x_subgrid(i,j),y_subgrid(i,j))
                !write(*,*) tt,i,j,temp_subgrid(i,j)
                enddo
            else
                temp_subgrid(i,j,:)=NODATA_value
                traveltime_subgrid(i,j,tt,:,:)=NODATA_value
            endif
            enddo
            enddo
        endif
        !write(unit_logfile,'(a,3f12.3)') 'Mean, min and max grid concentration: ',sum(temp_subgrid)/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index),minval(temp_subgrid),maxval(temp_subgrid)
        do i_pollutant=1,n_pollutant_loop
            temp_sum_subgrid(i_pollutant)=0.
            count=0
            do j=1,subgrid_dim(y_dim_index)
            do i=1,subgrid_dim(x_dim_index)
            if (use_subgrid(i,j,source_index)) then
                temp_sum_subgrid(i_pollutant)=temp_sum_subgrid(i_pollutant)+temp_subgrid(i,j,i_pollutant)
                count=count+1
            endif
            enddo
            enddo
            if (count.gt.0) then
                temp_sum_subgrid(i_pollutant)=temp_sum_subgrid(i_pollutant)/count
            else
                temp_sum_subgrid(i_pollutant)=0
            endif            
            write(unit_logfile,'(a,3f12.3)') 'Mean concentration '//trim(pollutant_file_str(pollutant_loop_index(i_pollutant)))//': ',temp_sum_subgrid(i_pollutant)
        enddo
        subgrid(:,:,tt,proxy_subgrid_index,source_index,:)=temp_subgrid

        
        enddo !time loop
        
    enddo !subsource_index
    
    !Combine the subsources in the dispersion if required
    !if (combine_emission_subsources_during_dispersion(source_index).and.n_subsource(source_index).gt.1) then
    !    do subsource_index=2,n_subsource(n_source_index)
    !        subgrid(:,:,:,proxy_subgrid_index,source_index,1)=subgrid(:,:,:,proxy_subgrid_index,source_index,1)+subgrid(:,:,:,proxy_subgrid_index,source_index,subsource_index)
    !    enddo
    !    n_subsource(source_index)=1
    !endif
      
    if (allocated(trajectory_subgrid)) deallocate(trajectory_subgrid)
    if (allocated(temp_emission_subgrid)) deallocate(temp_emission_subgrid)
    if (allocated(temp_subgrid)) deallocate(temp_subgrid)
    if (allocated(temp_FF_subgrid)) deallocate(temp_FF_subgrid)
    if (allocated(temp_FF_emission_subgrid)) deallocate(temp_FF_emission_subgrid)
    if (allocated(temp_subgrid)) deallocate(temp_subgrid)
    if (allocated(traveltime_temp_target_subgrid)) deallocate(traveltime_temp_target_subgrid)

    end subroutine uEMEP_subgrid_dispersion
    
    
!==========================================================================
!   uEMEP model grid_proxy_integral
!   Calculates the vertically integrated proxy concentrations based on Gaussian distribution
!   Generic for all sources and subsources
!   Seperate routine to the target grid as the integrated form may be of lower resolution
!==========================================================================
    subroutine uEMEP_subgrid_dispersion_integral(source_index)
    
    use uEMEP_definitions

    implicit none
    
    integer i,j,k
    integer         source_index
    character(256)  temp_name
    logical         exists
    integer         jj,ii,tt,tt_emis
    real            distance_subgrid
    real            dx_temp,dy_temp
    integer         i_start,i_end,j_start,j_end,t_start,t_end
    integer         i_cross,j_cross
    integer         i_cross_integral,j_cross_integral,i_cross_target_integral,j_cross_target_integral
    real            u_subgrid_loc,v_subgrid_loc,kz_subgrid_loc
    real            cos_subgrid_loc,sin_subgrid_loc,FF_loc,FF_zc_loc,zc_loc
    integer         i_nc,j_nc
    integer         subsource_index
    real            ay_loc,by_loc,az_loc,bz_loc,sig_y_0_loc,sig_z_0_loc,sig_y_00_loc,sig_z_00_loc,h_emis_loc,z_rec_loc,sig_z_loc,sig_y_loc,h_mix_loc
    real            xpos_limit,ypos_limit
    real            x_downwind,y_downwind
    real            xpos_area_min,xpos_area_max,ypos_area_min,ypos_area_max
    real            xpos_subgrid,ypos_subgrid
    real            xpos_emission_subgrid,ypos_emission_subgrid
    real            xpos_integral_subgrid,ypos_integral_subgrid
    real            distance_emission_subgrid_min

    integer         traj_max_index
    logical         valid_traj
    real            traj_step_size,x_loc,y_loc,invL_loc,FFgrid_loc,logz0_loc,u_star0_loc,FF10_loc
    real            z0_temp,h_temp
    
    integer i_pollutant
    
    !real, allocatable :: temp_emission_subgrid(:,:)
    !real, allocatable :: temp_subgrid(:,:)
    real, allocatable :: temp_FF_subgrid(:,:)
    real, allocatable :: temp_FF_emission_subgrid(:,:)
    real, allocatable :: trajectory_subgrid(:,:,:,:)
    real, allocatable :: angle_diff(:,:)

    !functions
    !real gauss_plume_second_order_rotated_func
    real gauss_plume_second_order_rotated_integral_func
    !real gauss_plume_cartesian_func
    real gauss_plume_cartesian_integral_func
    real gauss_plume_cartesian_trajectory_integral_func
    real gauss_plume_cartesian_sigma_integral_func
   
    
    allocate (temp_FF_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index))) 
    allocate (angle_diff(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index))) 
    allocate (temp_FF_emission_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index)))     
        
    temp_FF_subgrid=0.
    temp_FF_emission_subgrid=0.
    angle_diff=0.

    xpos_limit=dgrid_nc(lon_nc_index)/2.*EMEP_grid_interpolation_size
    ypos_limit=dgrid_nc(lat_nc_index)/2.*EMEP_grid_interpolation_size

    !Minimum distance for dispersion set to  half of an emission grid diagonal weighted so the circle has the same area as the square with that diagonal
    distance_emission_subgrid_min=sqrt(emission_subgrid_delta(x_dim_index,source_index)*emission_subgrid_delta(x_dim_index,source_index) &
        +emission_subgrid_delta(y_dim_index,source_index)*emission_subgrid_delta(y_dim_index,source_index))/2./sqrt(2.)*4./3.14159

    do subsource_index=1,n_subsource(source_index)

        call uEMEP_set_dispersion_params_simple(source_index,subsource_index)
    
        !Set local dispersion parameters prior to time loop for use in annual data
        ay_loc=ay(source_index,subsource_index)
        by_loc=by(source_index,subsource_index)
        az_loc=az(source_index,subsource_index)
        bz_loc=bz(source_index,subsource_index)
        sig_y_00_loc=sig_y_00(source_index,subsource_index)
        sig_z_00_loc=sig_z_00(source_index,subsource_index)
        h_emis_loc=h_emis(source_index,subsource_index)
        z_rec_loc=z_rec(source_index,subsource_index)

        !write(*,*) z_rec_loc,ay_loc,by_loc,az_loc,bz_loc,sig_y_0_loc,sig_z_0_loc
    
        write(unit_logfile,'(a,i3)')'Calculating proxy integral concentration data for '//trim(source_file_str(source_index))//' with subsource index ',subsource_index

        !Set the start and end times of the loop
        t_start=1
        t_end=integral_subgrid_dim(t_dim_index)
    
        !Loop through the time
        do tt=t_start,t_end
        
        integral_subgrid(:,:,tt,:,source_index,:)=0.
        temp_FF_subgrid=0.
        
         !Set the last meteo data subgrid in the case when the internal time loop is used
        if (.not.use_single_time_loop_flag) then
            if (tt.gt.t_start) then
                last_meteo_subgrid(:,:,:)=meteo_subgrid(:,:,tt-1,:)
            else
                last_meteo_subgrid(:,:,:)=meteo_subgrid(:,:,tt,:)      
            endif           
        endif

        !Precalculate information for the trajectory model
        !Maxium number of trajectory steps and size of steps based on the integral (meteorology) loop size
        if (use_trajectory_flag(source_index)) then

            traj_step_size=min(integral_subgrid_delta(x_dim_index),integral_subgrid_delta(y_dim_index))*traj_step_scale
            traj_max_index=floor(max(integral_subgrid_loop_index(x_dim_index),integral_subgrid_loop_index(y_dim_index))/traj_step_scale)
            !if (use_downwind_position_flag) traj_max_index=traj_max_index*2

            if (tt.eq.t_start) write(unit_logfile,'(a,f12.1,i)') 'Trajectory step (m) and dimensions: ', traj_step_size, traj_max_index
            if (.not.allocated(trajectory_subgrid)) allocate(trajectory_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),traj_max_index,2))
            
            trajectory_subgrid=NODATA_value
            
            !Loop through the emissions and create trajectories for all emissions source grids
            do j=1,emission_subgrid_dim(y_dim_index,source_index)
            do i=1,emission_subgrid_dim(x_dim_index,source_index)
            
                if (emission_subgrid(i,j,tt,source_index,subsource_index).ne.0) then
                    call uEMEP_calculate_all_trajectory(x_emission_subgrid(i,j,source_index),y_emission_subgrid(i,j,source_index),tt, &
                        traj_max_index,traj_step_size,trajectory_subgrid(i,j,:,x_dim_index),trajectory_subgrid(i,j,:,y_dim_index))
                endif
            
            enddo
            enddo
       
        endif

        !Create a temporary wind speed subgrid for each hour
        temp_FF_subgrid=0.
        do j_cross=1,integral_subgrid_dim(y_dim_index)
        do i_cross=1,integral_subgrid_dim(x_dim_index)
            z0_temp=exp(meteo_subgrid(i_cross,j_cross,tt,logz0_subgrid_index))
            h_temp=h_emis(source_index,subsource_index)
            if (annual_calculations.and.wind_level_integral_flag.eq.1) then
                temp_FF_subgrid(i_cross,j_cross)=1./meteo_subgrid(i_cross,j_cross,tt,inv_FFgrid_subgrid_index)
            elseif (annual_calculations.and.wind_level_integral_flag.eq.2) then
                temp_FF_subgrid(i_cross,j_cross)=1./meteo_subgrid(i_cross,j_cross,tt,inv_FFgrid_subgrid_index)*(1.-(log((H_meteo+z0_temp)/z0_temp)-log((h_temp+z0_temp)/z0_temp))/log((H_meteo+z0_temp)/z0_temp))
            elseif (annual_calculations.and.wind_level_integral_flag.eq.3) then
                temp_FF_subgrid(i_cross,j_cross)=1./meteo_subgrid(i_cross,j_cross,tt,inv_FF10_subgrid_index)
            elseif (annual_calculations.and.wind_level_integral_flag.eq.4) then
                temp_FF_subgrid(i_cross,j_cross)=1./meteo_subgrid(i_cross,j_cross,tt,inv_FF10_subgrid_index)*(1.-(log((10.+z0_temp)/z0_temp)-log((h_temp+z0_temp)/z0_temp))/log((10.+z0_temp)/z0_temp))
            elseif (hourly_calculations.and.wind_level_integral_flag.eq.1) then
                temp_FF_subgrid(i_cross,j_cross)=meteo_subgrid(i_cross,j_cross,tt,FFgrid_subgrid_index)
            elseif (hourly_calculations.and.wind_level_integral_flag.eq.2) then
                temp_FF_subgrid(i_cross,j_cross)=meteo_subgrid(i_cross,j_cross,tt,FFgrid_subgrid_index)*(1.-(log((H_meteo+z0_temp)/z0_temp)-log((h_temp+z0_temp)/z0_temp))/log((H_meteo+z0_temp)/z0_temp))
            elseif (hourly_calculations.and.wind_level_integral_flag.eq.3) then
                temp_FF_subgrid(i_cross,j_cross)=meteo_subgrid(i_cross,j_cross,tt,FF10_subgrid_index)
            elseif (hourly_calculations.and.wind_level_integral_flag.eq.4) then
                temp_FF_subgrid(i_cross,j_cross)=meteo_subgrid(i_cross,j_cross,tt,FF10_subgrid_index)*(1.-(log((10.+z0_temp)/z0_temp)-log((h_temp+z0_temp)/z0_temp))/log((10.+z0_temp)/z0_temp))
            elseif (wind_level_integral_flag.eq.0) then
                temp_FF_subgrid(i_cross,j_cross)=1.
            elseif (wind_level_integral_flag.eq.5) then
                !Will set later based on sigma z centre of mass
                temp_FF_subgrid(i_cross,j_cross)=1.
            elseif (wind_level_integral_flag.eq.6) then
                !Will set later based on average sigma z centre of mass and z_emis
                temp_FF_subgrid(i_cross,j_cross)=1.
            else
                write(unit_logfile,'(a)') 'No valid wind_level_integral_flag selected. Stopping (uEMEP_subgrid_dispersion)'
                stop
            endif
            
            !Setting a minimum value for wind for dispersion purposes (cannot be zero)
            temp_FF_subgrid(i_cross,j_cross)=sqrt(temp_FF_subgrid(i_cross,j_cross)*temp_FF_subgrid(i_cross,j_cross)+FF_min_dispersion*FF_min_dispersion)

            !Finds the angle difference between the current and last meteo field for dispersion and implements meandering if selected
            if (hourly_calculations) then
                call delta_wind_direction (i_cross,j_cross,tt,meteo_subgrid(i_cross,j_cross,tt,FF10_subgrid_index),angle_diff(i_cross,j_cross))
            else
                angle_diff(i_cross,j_cross)=0.
            endif

            
        enddo
        enddo
                           
        !If wind level flag is set to 5, use of initial plume centre of mass, then set wind speed for each non-zero emission grid
        if (wind_level_integral_flag.eq.5) then
            temp_FF_emission_subgrid=0.
            do jj=1,emission_subgrid_dim(y_dim_index,source_index)
            do ii=1,emission_subgrid_dim(x_dim_index,source_index)
            if (sum(emission_subgrid(ii,jj,tt,source_index,:)).ne.0) then
                                    
                !Set the integral meteorological grid position for the emission position
                i_cross_integral=crossreference_emission_to_integral_subgrid(ii,jj,x_dim_index,source_index)
                j_cross_integral=crossreference_emission_to_integral_subgrid(ii,jj,y_dim_index,source_index)
                
                !Set the local variables
                logz0_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,logz0_subgrid_index)
                FF10_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,FF10_subgrid_index)
                sig_y_00_loc=emission_properties_subgrid(ii,jj,emission_sigy00_index,source_index)
                sig_z_00_loc=emission_properties_subgrid(ii,jj,emission_sigz00_index,source_index)
                h_emis_loc=emission_properties_subgrid(ii,jj,emission_h_index,source_index)
                h_mix_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,hmix_subgrid_index)
                
                if (annual_calculations) then
                    FF10_loc=1./meteo_subgrid(i_cross_integral,j_cross_integral,tt,inv_FF10_subgrid_index)
                endif

                !Set sig_0's at the emission position
                x_loc=0.
                call uEMEP_set_dispersion_sigma_simple(sig_z_00_loc,sig_y_00_loc,sigy_0_subgid_width_scale,emission_subgrid_delta(:,source_index),angle_diff(i_cross_integral,j_cross_integral),x_loc,sig_z_loc,sig_y_loc,sig_z_0_loc,sig_y_0_loc)                                        

                !Use the initial plume centre of mass to determine wind advection height
                call z_centremass_gauss_func(sig_z_0_loc,h_emis_loc,h_mix_loc,zc_loc)
                call u_profile_neutral_val_func(zc_loc,FF10_loc,10.,h_mix_loc,exp(logz0_loc),FF_loc,u_star0_loc)
                
                !Set a minimum wind speed based on traffic (if use_traffic_for_minFF_flag=T)
!                FF_loc=sqrt(FF_loc*FF_loc+emission_properties_subgrid(ii,jj,emission_minFF_index,source_index,subsource_index)*emission_properties_subgrid(ii,jj,emission_minFF_index,source_index,subsource_index))
                
                !Set the minimum wind speed 
                FF_loc=sqrt(FF_loc*FF_loc+FF_min_dispersion*FF_min_dispersion)

                temp_FF_emission_subgrid(ii,jj)=FF_loc
                !write(*,*) FF10_loc,FF_loc,zc_loc,sig_z_0_loc
                
            
            endif
            enddo
            enddo
        endif

        !Loop through the proxy integral grid
        do j=1,integral_subgrid_dim(y_dim_index)
        do i=1,integral_subgrid_dim(x_dim_index)
                            
                xpos_integral_subgrid=xproj_integral_subgrid(i,j)
                ypos_integral_subgrid=yproj_integral_subgrid(i,j)

                !Find the cross reference to the integral grid from the emission grid
                i_cross=crossreference_integral_to_emission_subgrid(i,j,x_dim_index,source_index)
                j_cross=crossreference_integral_to_emission_subgrid(i,j,y_dim_index,source_index)

                i_cross_target_integral=i
                j_cross_target_integral=j
                
                if (use_downwind_position_flag.and.hourly_calculations) then
                    
                    x_downwind=max(-1.,min(1.,meteo_subgrid(i_cross_target_integral,j_cross_target_integral,tt,cos_subgrid_index)*sqrt(2.)))
                    y_downwind=max(-1.,min(1.,meteo_subgrid(i_cross_target_integral,j_cross_target_integral,tt,sin_subgrid_index)*sqrt(2.)))                                     
                    i_end=min(ceiling(i_cross+1+(1.-x_downwind)*emission_subgrid_loop_index(x_dim_index,source_index)),emission_subgrid_dim(x_dim_index,source_index))
                    i_start=max(floor(i_cross-1-(1.+x_downwind)*emission_subgrid_loop_index(x_dim_index,source_index)),1)
                    j_end=min(ceiling(j_cross+1+(1.-y_downwind)*emission_subgrid_loop_index(y_dim_index,source_index)),emission_subgrid_dim(y_dim_index,source_index))
                    j_start=max(floor(j_cross-1-(1.+y_downwind)*emission_subgrid_loop_index(y_dim_index,source_index)),1)
                                                    
                    !Set new lon and lat limits to include the upwind source region
                    xpos_area_max=xpos_integral_subgrid+(1.-x_downwind)*xpos_limit/2.+emission_subgrid_dim(x_dim_index,source_index)
                    xpos_area_min=xpos_integral_subgrid-(1.+x_downwind)*xpos_limit/2.-emission_subgrid_dim(x_dim_index,source_index)
                    ypos_area_max=ypos_integral_subgrid+(1.-y_downwind)*ypos_limit/2.+emission_subgrid_dim(y_dim_index,source_index)
                    ypos_area_min=ypos_integral_subgrid-(1.+y_downwind)*ypos_limit/2.-emission_subgrid_dim(y_dim_index,source_index)

                else
                               
                    !Set the size of the loop region around the target cell to be up to integral_subgrid_loop_index
                    i_start=max(1,i_cross-emission_subgrid_loop_index(x_dim_index,source_index))
                    i_end=min(emission_subgrid_dim(x_dim_index,source_index),i_cross+emission_subgrid_loop_index(x_dim_index,source_index))
                    j_start=max(1,j_cross-emission_subgrid_loop_index(y_dim_index,source_index))
                    j_end=min(emission_subgrid_dim(y_dim_index,source_index),j_cross+emission_subgrid_loop_index(y_dim_index,source_index))
 
                    xpos_area_max=xpos_integral_subgrid+xpos_limit
                    xpos_area_min=xpos_integral_subgrid-xpos_limit
                    ypos_area_max=ypos_integral_subgrid+ypos_limit
                    ypos_area_min=ypos_integral_subgrid-ypos_limit

                endif
          
                !Loop through emission sub_grids in the nearby region
                do jj=j_start,j_end
                do ii=i_start,i_end
                            
                    if (sum(emission_subgrid(ii,jj,tt,source_index,:)).ne.0) then
                
                        !
                        xpos_emission_subgrid=xproj_emission_subgrid(ii,jj,source_index)
                        ypos_emission_subgrid=yproj_emission_subgrid(ii,jj,source_index)
                                                
                        if (xpos_emission_subgrid.ge.xpos_area_min.and.xpos_emission_subgrid.le.xpos_area_max &
                            .and.ypos_emission_subgrid.ge.ypos_area_min.and.ypos_emission_subgrid.le.ypos_area_max) then                  
               
                            !Determine meteorology grid position
                            i_cross_integral=crossreference_emission_to_integral_subgrid(ii,jj,x_dim_index,source_index)
                            j_cross_integral=crossreference_emission_to_integral_subgrid(ii,jj,y_dim_index,source_index)                          

                            if (hourly_calculations) then
                            
                                if (use_trajectory_flag(source_index)) then
                                
                                    !Calculate the minimum distance to the trajectory. Time consuming
                                    call uEMEP_minimum_distance_trajectory_fast(x_integral_subgrid(i,j),y_integral_subgrid(i,j),x_emission_subgrid(ii,jj,source_index),y_emission_subgrid(ii,jj,source_index), &
                                    traj_max_index,traj_step_size,trajectory_subgrid(ii,jj,:,x_dim_index),trajectory_subgrid(ii,jj,:,y_dim_index),x_loc,y_loc,valid_traj)
                                
                                else
                                                                
                                    !Set the local wind cos and sin values
                                    cos_subgrid_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,cos_subgrid_index)
                                    sin_subgrid_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,sin_subgrid_index)                                    

                                    !Determine the rotated along wind values
                                    x_loc=(x_integral_subgrid(i,j)-x_emission_subgrid(ii,jj,source_index))*cos_subgrid_loc+(y_integral_subgrid(i,j)-y_emission_subgrid(ii,jj,source_index))*sin_subgrid_loc
                                    y_loc=-(x_integral_subgrid(i,j)-x_emission_subgrid(ii,jj,source_index))*sin_subgrid_loc+(y_integral_subgrid(i,j)-y_emission_subgrid(ii,jj,source_index))*cos_subgrid_loc
                                                                        
                                   !If x is downwind then it is valid
                                    if (x_loc.ge.0) then
                                        valid_traj=.true.
                                    else
                                        valid_traj=.false.
                                    endif
                                                            
                                endif

                                !Calculate dispersion
                                if (valid_traj) then
                                    
                                    !Set the mixing height at the average of the emission and target position
                                    h_mix_loc=(meteo_subgrid(i_cross_integral,j_cross_integral,tt,hmix_subgrid_index)+meteo_subgrid(i_cross_target_integral,j_cross_target_integral,tt,hmix_subgrid_index))/2.
                                    !Set the local wind speed and other parameters at emission position
                                    FF_loc=temp_FF_subgrid(i_cross_integral,j_cross_integral)
                                    !L_loc=1./meteo_subgrid(i_cross_integral,j_cross_integral,tt,invL_subgrid_index)
                                    invL_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,invL_subgrid_index)
                                    FFgrid_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,FFgrid_subgrid_index)
                                    logz0_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,logz0_subgrid_index)
                                    u_star0_loc=max(meteo_subgrid(i_cross_integral,j_cross_integral,tt,ustar_subgrid_index),0.001)
                                    FF10_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,FF10_subgrid_index)
                                    sig_y_00_loc=emission_properties_subgrid(ii,jj,emission_sigy00_index,source_index)
                                    sig_z_00_loc=emission_properties_subgrid(ii,jj,emission_sigz00_index,source_index)
                                    h_emis_loc=emission_properties_subgrid(ii,jj,emission_h_index,source_index)

                                    if (wind_level_integral_flag.eq.5) then
                                        FF_loc=temp_FF_emission_subgrid(ii,jj)
                                    endif

                                    !Select method for assigning sigma
                                    if (stability_scheme_flag.eq.1) then
                                        call uEMEP_set_dispersion_sigma_simple(sig_z_00_loc,sig_y_00_loc,sigy_0_subgid_width_scale,emission_subgrid_delta(:,source_index),angle_diff(i_cross_integral,j_cross_integral),x_loc,sig_z_loc,sig_y_loc,sig_z_0_loc,sig_y_0_loc)
                                    endif
                                    if (stability_scheme_flag.eq.2) then
                                        call uEMEP_set_dispersion_sigma_PG(invL_loc,logz0_loc,sig_z_00_loc,sig_y_00_loc,sigy_0_subgid_width_scale,emission_subgrid_delta(:,source_index),angle_diff(i_cross_integral,j_cross_integral),x_loc,sig_z_loc,sig_y_loc,sig_z_0_loc,sig_y_0_loc)
                                    endif
                                    
                                    if (stability_scheme_flag.eq.3) then
                                        
                                        !Set initial values for sigma. Sig_y is set here
                                        call uEMEP_set_dispersion_sigma_simple(sig_z_00_loc,sig_y_00_loc,sigy_0_subgid_width_scale,emission_subgrid_delta(:,source_index),angle_diff(i_cross_integral,j_cross_integral),x_loc,sig_z_loc,sig_y_loc,sig_z_0_loc,sig_y_0_loc)
                                    
                                        !write(*,*) 'IN:  ',x_loc,sig_z_loc,FF_loc
                                        call uEMEP_set_dispersion_sigma_Kz(x_loc,sig_z_0_loc,sig_y_0_loc,sigy_0_subgid_width_scale,sig_z_loc,h_emis_loc,h_mix_loc,invL_loc,FF10_loc,10.,logz0_loc,emission_subgrid_delta(:,source_index),u_star0_loc,average_zc_h_in_Kz_flag,n_kz_iterations,sig_z_loc,sig_y_loc,FF_zc_loc)
                                        !write(*,*) 'OUT: ',x_loc,sig_z_loc,FF_loc
                                        
                                        sig_y_loc=sig_y_loc+x_loc*abs(angle_diff(i_cross_integral,j_cross_integral))
                                        
                                        !Use the average of the emisiion height and zc to determine wind speed. Is set to true if wind_level_flag=6
                                        if (wind_level_integral_flag.eq.6) then
                                           FF_loc=sqrt(FF_zc_loc*FF_zc_loc+FF_min_dispersion*FF_min_dispersion)
                                        endif
                                   
                                    endif

                                    if (stability_scheme_flag.eq.4) then
                                        call uEMEP_set_dispersion_sigma_Kz_emulator(h_emis_loc,invL_loc,logz0_loc,h_mix_loc,sig_z_00_loc,sig_y_00_loc,sigy_0_subgid_width_scale,emission_subgrid_delta(:,source_index),angle_diff(i_cross_integral,j_cross_integral),x_loc,sig_z_loc,sig_y_loc,sig_z_0_loc,sig_y_0_loc)
                                    endif
                                    
                                    if (wind_level_integral_flag.eq.6.and.stability_scheme_flag.ne.3) then
                                        call z_centremass_gauss_func(sig_z_loc,h_emis_loc,h_mix_loc,zc_loc)
                                        zc_loc=(h_emis_loc+zc_loc)/2.
                                        call u_profile_neutral_val_func(zc_loc,FF10_loc,10.,h_mix_loc,exp(logz0_loc),FF_zc_loc,u_star0_loc)
                                        FF_loc=sqrt(FF_zc_loc*FF_zc_loc+FF_min_dispersion*FF_min_dispersion)
                                    endif
                                    
                                    !Calculate the dispersion
                                    integral_subgrid(i,j,tt,hsurf_integral_subgrid_index,source_index,:)=integral_subgrid(i,j,tt,hsurf_integral_subgrid_index,source_index,:) &
                                        + gauss_plume_cartesian_sigma_integral_func(x_loc,y_loc,h_emis_loc,z_rec_loc,sig_z_loc,sig_y_loc,h_mix_loc,FF_loc,0.,H_emep) &
                                        * emission_subgrid(ii,jj,tt,source_index,:)
                                                                        
                                endif
                                
                            else
                                !If annual calculations
                                distance_subgrid=sqrt((x_emission_subgrid(ii,jj,source_index)-x_integral_subgrid(i,j))*(x_emission_subgrid(ii,jj,source_index)-x_integral_subgrid(i,j)) &
                                    +(y_emission_subgrid(ii,jj,source_index)-y_integral_subgrid(i,j))*(y_emission_subgrid(ii,jj,source_index)-y_integral_subgrid(i,j)))

                                if (wind_level_integral_flag.eq.5) then
                                    FF_loc=temp_FF_emission_subgrid(ii,jj)
                                else
                                    FF_loc=temp_FF_subgrid(i_cross_integral,j_cross_integral)
                                endif

                                integral_subgrid(i,j,tt,hsurf_integral_subgrid_index,source_index,:)=integral_subgrid(i,j,hsurf_integral_subgrid_index,tt,source_index,:) &
                                    +emission_subgrid(ii,jj,tt,source_index,:) &
                                    *gauss_plume_second_order_rotated_integral_func(distance_subgrid,z_rec_loc,ay_loc,by_loc,az_loc,bz_loc,sig_y_0_loc,sig_z_0_loc,h_emis_loc,0.,H_emep) &
                                    /FF_loc
        
                            endif
                        
                        endif
                    
                    endif
                    
                    enddo
                    enddo
                !write(*,'(6i6,f)') i,j,i_start,i_end,j_start,j_end,integral_subgrid(i,j,tt,source_index,subsource_index)
        enddo
            !write(*,'(3a,i5,a,i5,a,i3,a,i3)') 'Gridding ',trim(source_file_str(source_index)),' integral proxy',j,' of ',integral_subgrid_dim(y_dim_index),' and ',subsource_index,' of ',n_subsource(source_index)
        enddo
        
        if (mod(j,1).eq.0) write(*,'(3a,i5,a,i5,a,i3,a,i3)') 'Integral gridding ',trim(source_file_str(source_index)),' proxy for hour ',tt,' of ',subgrid_dim(t_dim_index),' and ',subsource_index,' of ',n_subsource(source_index)
        enddo !time loop      

    enddo !subsource_index
               
    !if (combine_emission_subsources_during_dispersion(source_index).and.n_subsource(source_index).gt.1) then
    !    do subsource_index=2,n_subsource(n_source_index)
    !        integral_subgrid(:,:,:,:,source_index,1)=integral_subgrid(:,:,:,:,source_index,1)+integral_subgrid(:,:,:,:,source_index,subsource_index)
    !    enddo
    !    n_subsource(source_index)=1
    !endif


    if (allocated(trajectory_subgrid)) deallocate(trajectory_subgrid)
    !if (allocated(temp_emission_subgrid)) deallocate(temp_emission_subgrid)
    !if (allocated(temp_subgrid)) deallocate(temp_subgrid)
    if (allocated(temp_FF_subgrid)) deallocate(temp_FF_subgrid)

    end subroutine uEMEP_subgrid_dispersion_integral