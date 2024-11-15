module subgrid_deposition

    use uemep_configuration
    use set_dispersion_parameters, only: delta_wind_direction, &
        uEMEP_set_dispersion_sigma_PG, uEMEP_set_dispersion_sigma_simple, &
        uEMEP_set_dispersion_params_PG, uEMEP_set_dispersion_params_simple
    use local_trajectory, only: uEMEP_calculate_all_trajectory, &
        uEMEP_minimum_distance_trajectory_fast
    use dispersion_functions, only: gauss_plume_cartesian_sigma_integral_func, &
        gauss_plume_cartesian_sigma_func
    use kz_functions, only: z_centremass_gauss_func, u_profile_neutral_val_func, &
        uEMEP_set_dispersion_sigma_Kz
    use area_interpolation_functions, only: area_weighted_interpolation_function, &
        area_weighted_extended_interpolation_function
    use mod_rargsort, only: rargsort

    implicit none
    private

    public :: uEMEP_subgrid_deposition

contains

! uEMEP_subgrid_deposition.f90

!Calculations the dry deposition and source depletion of a plume
    subroutine uEMEP_subgrid_deposition(source_index)

        use uEMEP_definitions

        implicit none

        integer source_index

        integer i,j,ii,jj,tt,k
        integer i_pollutant

        !Define the target subgrid to be the same as the emission subgrid
        real, allocatable :: target_subgrid(:,:,:)
        real, allocatable :: target_deposition_subgrid(:,:,:,:)
        real, allocatable :: x_target_subgrid(:,:)
        real, allocatable :: y_target_subgrid(:,:)
        real, allocatable :: traveltime_target_subgrid(:,:,:,:)
        real :: target_subgrid_delta(2)
        integer target_subgrid_dim(n_dim_index)

        !define the temporary arrays for meteo
        real, allocatable :: temp_FF_subgrid(:,:)
        real, allocatable :: temp_FF_integral_subgrid(:,:)
        real, allocatable :: trajectory_vector(:,:)
        real, allocatable :: angle_diff(:,:)

        real temp_sum_subgrid(n_pollutant_loop)

        integer traj_max_index
        logical valid_traj
        real traj_step_size,x_loc,y_loc,FFgrid_loc,logz0_loc,u_star0_loc,FF10_loc,zc_loc,invL_loc
        real ay_loc,by_loc,az_loc,bz_loc,sig_y_0_loc,sig_z_0_loc,sig_y_00_loc,sig_z_00_loc,h_emis_loc,z_rec_loc,sig_z_loc,sig_y_loc,h_mix_loc
        real FF_loc,FF_zc_loc,precip_loc
        real FF_integral_loc

        integer n_plume_subgrid_max
        parameter (n_plume_subgrid_max=100000)
        integer plume_crossreference(n_plume_subgrid_max,2)
        real plume_distance(n_plume_subgrid_max,2)
        real temp_plume_distance(n_plume_subgrid_max)
        integer sorted_plume_index(n_plume_subgrid_max)
        real plume_source(n_pollutant_loop)
        integer plume_count,max_plume_count

        real xpos_limit,ypos_limit
        real distance_subgrid,distance_subgrid_min,distance_emission_subgrid_min

        integer i_target_start,i_target_end,j_target_start,j_target_end
        integer t_start,t_end

        integer i_cross,j_cross
        integer i_cross_integral,j_cross_integral
        integer i_cross_deposition,j_cross_deposition
        integer i_cross_target_integral,j_cross_target_integral

        integer subsource_index
        real subgrid_internal,subgrid_internal_integral
        real subgrid_internal_pollutant(n_pollutant_loop),drydepo_internal_pollutant(n_pollutant_loop),wetdepo_internal_pollutant(n_pollutant_loop)
        real drydepo_internal_integral
        real vertical_integral_internal,wetdepo_internal,wetdepo_internal_integral

        real xpos_emission_subgrid,ypos_emission_subgrid
        real xpos_area_max,xpos_area_min,ypos_area_max,ypos_area_min
        real xpos_target_subgrid,ypos_target_subgrid
        integer target_subgrid_dim_min(2),target_subgrid_dim_max(2)

        real plume_vertical_integral(n_integral_subgrid_index)
        real plume_vertical_integral_pollutant(n_integral_subgrid_index)
        real, allocatable :: target_vertical_integral_subgrid(:,:,:)
        integer i_integral


        integer count
        !Functions
        !integer rargsort(n_plume_subgrid_max)

        write(unit_logfile,'(A)') ''
        write(unit_logfile,'(A)') '================================================================'
        write(unit_logfile,'(A)') 'Calculating deposition and dispersion (uEMEP_subgrid_deposition)'
        write(unit_logfile,'(A)') '================================================================'

        !Set up the target subgrid to be the same as the emission subgrid for that source.
        target_subgrid_dim(:)=emission_subgrid_dim(:,source_index)
        target_subgrid_delta(:)=emission_subgrid_delta(:,source_index)

        !Deposition and plume depletion must be calculated also in the buffer zone
        if (calculate_source_depletion_flag) then
            !Need to calculate on all the grids to get the depletion
            target_subgrid_dim_min(:)=1;target_subgrid_dim_max(:)=target_subgrid_dim(1:2)
        else
            !Limit the target grid to slightly larger than concentration grid
            target_subgrid_dim_min(x_dim_index)=-1+1+floor((subgrid_min(x_dim_index)-emission_subgrid_min(x_dim_index,source_index))/emission_subgrid_delta(x_dim_index,source_index))
            target_subgrid_dim_min(y_dim_index)=-1+1+floor((subgrid_min(y_dim_index)-emission_subgrid_min(y_dim_index,source_index))/emission_subgrid_delta(y_dim_index,source_index))
            target_subgrid_dim_max(x_dim_index)=+1+1+ceiling((subgrid_max(x_dim_index)-emission_subgrid_min(x_dim_index,source_index))/emission_subgrid_delta(x_dim_index,source_index))
            target_subgrid_dim_max(y_dim_index)=+1+1+ceiling((subgrid_max(y_dim_index)-emission_subgrid_min(y_dim_index,source_index))/emission_subgrid_delta(y_dim_index,source_index))
        endif

        !Allocate the target subgrid to be the same as the emission subgrid but no time dimension
        if (allocated(target_subgrid)) deallocate (target_subgrid)
        if (allocated(target_deposition_subgrid)) deallocate (target_deposition_subgrid)
        if (allocated(x_target_subgrid)) deallocate (x_target_subgrid)
        if (allocated(y_target_subgrid)) deallocate (y_target_subgrid)
        if (allocated(traveltime_target_subgrid)) deallocate (traveltime_target_subgrid)
        if (.not.allocated(target_subgrid)) allocate (target_subgrid(target_subgrid_dim(x_dim_index),target_subgrid_dim(y_dim_index),n_pollutant_loop))
        if (.not.allocated(target_deposition_subgrid)) allocate (target_deposition_subgrid(target_subgrid_dim(x_dim_index),target_subgrid_dim(y_dim_index),n_deposition_index,n_pollutant_loop))
        if (.not.allocated(x_target_subgrid)) allocate (x_target_subgrid(target_subgrid_dim(x_dim_index),target_subgrid_dim(y_dim_index)))
        if (.not.allocated(y_target_subgrid)) allocate (y_target_subgrid(target_subgrid_dim(x_dim_index),target_subgrid_dim(y_dim_index)))
        if (.not.allocated(traveltime_target_subgrid)) allocate (traveltime_target_subgrid(target_subgrid_dim(x_dim_index),target_subgrid_dim(y_dim_index),2,n_pollutant_loop))

        if (adjust_wetdepo_integral_to_lowest_layer_flag.or.local_subgrid_method_flag.eq.1) then
            if (allocated(target_vertical_integral_subgrid)) deallocate (target_vertical_integral_subgrid)
            if (.not.allocated(target_vertical_integral_subgrid)) allocate (target_vertical_integral_subgrid(target_subgrid_dim(x_dim_index),target_subgrid_dim(y_dim_index),3))
        endif

        x_target_subgrid(:,:)=x_emission_subgrid(:,:,source_index)
        y_target_subgrid(:,:)=y_emission_subgrid(:,:,source_index)

        !Allocate temporary wind speed subgrid
        allocate (temp_FF_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index)))
        allocate (temp_FF_integral_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index)))
        allocate (angle_diff(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index)))

        !Set the x and y position limits to coincide to half the EMEP grid (refered here as lon and lat but can be also LCC projection) times the number of grids
        xpos_limit=dgrid_nc(lon_nc_index)/2.*EMEP_grid_interpolation_size*local_fraction_grid_size_scaling
        ypos_limit=dgrid_nc(lat_nc_index)/2.*EMEP_grid_interpolation_size*local_fraction_grid_size_scaling

        !Minimum distance for travel time calculation set to half of a grid diagonal weighted so the circle has the same area as the square with that diagonal
        distance_subgrid_min=sqrt(subgrid_delta(x_dim_index)*subgrid_delta(x_dim_index)+subgrid_delta(y_dim_index)*subgrid_delta(y_dim_index))/2./sqrt(2.)*4./3.14159
        !Minimum distance for dispersion set to  half of an emission grid diagonal weighted so the circle has the same area as the square with that diagonal
        distance_emission_subgrid_min=sqrt(emission_subgrid_delta(x_dim_index,source_index)*emission_subgrid_delta(x_dim_index,source_index) &
            +emission_subgrid_delta(y_dim_index,source_index)*emission_subgrid_delta(y_dim_index,source_index))/2./sqrt(2.)*4./3.14159

        !Set the subsource_index to 1 so no additional subsources in these routines
        subsource_index=1


        !Set local dispersion parameters to be used only in the annual calculation, overwritten in the hourly files
        !if (annual_calculations) then
        call uEMEP_set_dispersion_params_simple(source_index,subsource_index)
        ay_loc=ay(source_index,subsource_index)
        by_loc=by(source_index,subsource_index)
        az_loc=az(source_index,subsource_index)
        bz_loc=bz(source_index,subsource_index)
        sig_y_00_loc=sig_y_00(source_index,subsource_index)
        sig_z_00_loc=sig_z_00(source_index,subsource_index)
        h_emis_loc=h_emis(source_index,subsource_index)
        z_rec_loc=z_rec(source_index,subsource_index)
        !endif

        write(unit_logfile,'(a,i3)')'Calculating deposition and dispersion data for '//trim(source_file_str(source_index))

        !Set the target grid loop variables
        j_target_start=1;j_target_end=target_subgrid_dim(y_dim_index)
        i_target_start=1;i_target_end=target_subgrid_dim(x_dim_index)

        j_target_start=target_subgrid_dim_min(y_dim_index);j_target_end=target_subgrid_dim_max(y_dim_index)
        i_target_start=target_subgrid_dim_min(x_dim_index);i_target_end=target_subgrid_dim_max(x_dim_index)

        !write(*,*) j_target_start,j_target_end,i_target_start,i_target_end
        !stop
        !Set the start and end times of the loop
        t_start=1;t_end=subgrid_dim(t_dim_index)

        !Loop through the time
        do tt=t_start,t_end

            !Initialise the final grid to 0
            subgrid(:,:,tt,proxy_subgrid_index,source_index,:)=0.
            subgrid(:,:,tt,drydepo_local_subgrid_index,source_index,:)=0.
            subgrid(:,:,tt,wetdepo_local_subgrid_index,source_index,:)=0.
            integral_subgrid(:,:,tt,:,source_index,:)=0.
            !deposition_subgrid(:,:,tt,:,:)=0.

            !Initialise the target grid to 0
            target_subgrid=0.
            target_deposition_subgrid=0.
            target_vertical_integral_subgrid=0.

            !Set the last meteo data subgrid in the case when the internal time loop is used
            if (.not.use_single_time_loop_flag) then
                if (tt.gt.t_start) then
                    last_meteo_subgrid(:,:,:)=meteo_subgrid(:,:,tt-1,:)
                else
                    last_meteo_subgrid(:,:,:)=meteo_subgrid(:,:,tt,:)
                endif
            endif

            !Finds the angle difference between the current and last meteo field for dispersion and implements meandering if selected
            do j_cross=1,integral_subgrid_dim(y_dim_index)
                do i_cross=1,integral_subgrid_dim(x_dim_index)
                    if (hourly_calculations) then
                        call delta_wind_direction (i_cross,j_cross,tt,meteo_subgrid(i_cross,j_cross,tt,FF10_subgrid_index),angle_diff(i_cross,j_cross))
                    else
                        angle_diff(i_cross,j_cross)=0.
                    endif
                enddo
            enddo

            !Create a temporary wind speed subgrid for each hour
            call uEMEP_create_wind_field(temp_FF_subgrid,angle_diff,wind_level_flag,source_index,subsource_index,tt)

            !Create an integral wind speed subgrid
            call uEMEP_create_wind_field(temp_FF_integral_subgrid,angle_diff,wind_level_integral_flag,source_index,subsource_index,tt)


            !Define the trajectory and its length
            !Maxium number of trajectory steps and size of steps based on the integral (meteorology) loop size
            if (use_trajectory_flag(source_index)) then

                traj_step_size=min(integral_subgrid_delta(x_dim_index),integral_subgrid_delta(y_dim_index))*traj_step_scale
                traj_max_index=floor(max(integral_subgrid_loop_index(x_dim_index),integral_subgrid_loop_index(y_dim_index))/traj_step_scale)

                if (tt.eq.t_start) write(unit_logfile,'(a,f12.1,i)') 'Trajectory step (m) and dimensions: ', traj_step_size, traj_max_index
                if (.not.allocated(trajectory_vector)) allocate(trajectory_vector(traj_max_index,2))
            endif

            max_plume_count=0
            !Set up the deposition target subgrid based on the emission subgrid.
            !This works when the deposition subgrid is larger than the target emission subgrid
            !When the deposition subgrid is smaller than the target then we need to take the average of the grids
            !This should be area weighted but is not
            if (calculate_deposition_flag) then
                if (deposition_subgrid_delta(x_dim_index).lt.target_subgrid_delta(x_dim_index)) then
                    do jj=j_target_start,j_target_end
                        do ii=i_target_start,i_target_end
                            i_cross_deposition=crossreference_emission_to_deposition_subgrid(ii,jj,x_dim_index,source_index)
                            j_cross_deposition=crossreference_emission_to_deposition_subgrid(ii,jj,y_dim_index,source_index)
                            !i_loop_deposition=floor(target_subgrid_delta(x_dim_index)/deposition_subgrid_delta(x_dim_index)/2.)
                            !target_deposition_subgrid(ii,jj,:,:)=target_deposition_subgrid(ii,jj,:,:)+deposition_subgrid(i_cross_deposition,j_cross_deposition,tt,:,:)
                            !write(*,*) deposition_subgrid(i_cross_deposition,j_cross_deposition,tt,:,:)
                            do i_pollutant=1,n_pollutant_loop
                                target_deposition_subgrid(ii,jj,vd_index,i_pollutant)=target_deposition_subgrid(ii,jj,vd_index,i_pollutant)+ &
                                    area_weighted_extended_interpolation_function(x_deposition_subgrid,y_deposition_subgrid,deposition_subgrid(:,:,tt,vd_index,i_pollutant) &
                                    ,deposition_subgrid_dim(x_dim_index),deposition_subgrid_dim(y_dim_index),deposition_subgrid_delta(:),x_target_subgrid(ii,jj),y_target_subgrid(ii,jj),target_subgrid_delta)
                            enddo
                        enddo
                    enddo

                else

                    do jj=1,emission_subgrid_dim(y_dim_index,source_index)
                        do ii=1,emission_subgrid_dim(x_dim_index,source_index)
                            i_cross_deposition=crossreference_emission_to_deposition_subgrid(ii,jj,x_dim_index,source_index)
                            j_cross_deposition=crossreference_emission_to_deposition_subgrid(ii,jj,y_dim_index,source_index)
                            target_deposition_subgrid(ii,jj,:,:)=deposition_subgrid(i_cross_deposition,j_cross_deposition,tt,:,:)
                            !write(*,*) deposition_subgrid(i_cross_deposition,j_cross_deposition,tt,:,:)
                        enddo
                    enddo

                endif
            endif

            !Loop through all the emission subgrids
            do jj=1,emission_subgrid_dim(y_dim_index,source_index)
                do ii=1,emission_subgrid_dim(x_dim_index,source_index)

                    !Allocate temporary emission (plume source) for this grid and time
                    plume_source(:)=emission_subgrid(ii,jj,tt,source_index,:)

                    !Only calculate if there is an emissions

                    if (sum(plume_source(:)).ne.0) then

                        !Calculate the trajectory for this emission source
                        if (use_trajectory_flag(source_index)) then

                            trajectory_vector=NODATA_value

                            !Calculate the trajectory for this emission grid
                            call uEMEP_calculate_all_trajectory(x_emission_subgrid(ii,jj,source_index),y_emission_subgrid(ii,jj,source_index),tt, &
                                traj_max_index,traj_step_size,trajectory_vector(:,x_dim_index),trajectory_vector(:,y_dim_index))
                        else
                            !Create an artificial trajectory here for the straight line case by making the trajectory step size half the size of the integral grid
                            !Not working!!!!
                            traj_step_size=min(integral_subgrid_max(x_dim_index)-integral_subgrid_min(x_dim_index),integral_subgrid_max(y_dim_index)-integral_subgrid_min(y_dim_index))/2.
                            call uEMEP_calculate_all_trajectory(x_emission_subgrid(ii,jj,source_index),y_emission_subgrid(ii,jj,source_index),tt, &
                                traj_max_index,traj_step_size,trajectory_vector(:,x_dim_index),trajectory_vector(:,y_dim_index))

                        endif

                        !Set the integral meteorological grid position for the emission position
                        i_cross_integral=crossreference_emission_to_integral_subgrid(ii,jj,x_dim_index,source_index)
                        j_cross_integral=crossreference_emission_to_integral_subgrid(ii,jj,y_dim_index,source_index)
                        i_cross_integral=min(max(1,i_cross_integral),integral_subgrid_dim(x_dim_index))
                        j_cross_integral=min(max(1,j_cross_integral),integral_subgrid_dim(y_dim_index))

                        !Set the local wind speed and other parameters at emission position
                        FF_loc=temp_FF_subgrid(i_cross_integral,j_cross_integral)
                        FF_integral_loc=temp_FF_integral_subgrid(i_cross_integral,j_cross_integral)
                        !L_loc=1./meteo_subgrid(i_cross_integral,j_cross_integral,tt,invL_subgrid_index)
                        invL_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,invL_subgrid_index)
                        FFgrid_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,FFgrid_subgrid_index)
                        logz0_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,logz0_subgrid_index)
                        u_star0_loc=max(meteo_subgrid(i_cross_integral,j_cross_integral,tt,ustar_subgrid_index),ustar_min)
                        FF10_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt,FF10_subgrid_index)
                        sig_y_00_loc=emission_properties_subgrid(ii,jj,emission_sigy00_index,source_index)
                        sig_z_00_loc=emission_properties_subgrid(ii,jj,emission_sigz00_index,source_index)
                        h_emis_loc=emission_properties_subgrid(ii,jj,emission_h_index,source_index)


                        !MORE COMPLICATED BECAUSE THE PLUME DEPLETION MUST ALSO BE CALCULATED OUTSIDE THE TARGET REGION!!
                        !SO THE LOOP NEEDS TO BE THE EMISSION LOOP NOT JUST THE TARGET LOOP
                        !Loop through the target grids and find the distance of each target grid along the plume centre line

                        if (calculate_deposition_flag.and.calculate_source_depletion_flag) then
                            !Set the size of the loop region around the emission cell to be up to emission_subgrid_loop_index
                            i_target_start=max(1,ii-emission_subgrid_loop_index(x_dim_index,source_index))
                            i_target_end=min(emission_subgrid_dim(x_dim_index,source_index),ii+emission_subgrid_loop_index(x_dim_index,source_index))
                            j_target_start=max(1,jj-emission_subgrid_loop_index(y_dim_index,source_index))
                            j_target_end=min(emission_subgrid_dim(y_dim_index,source_index),jj+emission_subgrid_loop_index(y_dim_index,source_index))
                        else
                            !Usea smaller area
                            i_target_start=max(1,ii-emission_subgrid_loop_index(x_dim_index,source_index))
                            i_target_end=min(emission_subgrid_dim(x_dim_index,source_index),ii+emission_subgrid_loop_index(x_dim_index,source_index))
                            j_target_start=max(1,jj-emission_subgrid_loop_index(y_dim_index,source_index))
                            j_target_end=min(emission_subgrid_dim(y_dim_index,source_index),jj+emission_subgrid_loop_index(y_dim_index,source_index))

                        endif

                        !Set the emission limits (EMEP projection ) surrounding the target grid
                        xpos_emission_subgrid=xproj_emission_subgrid(ii,jj,source_index)
                        ypos_emission_subgrid=yproj_emission_subgrid(ii,jj,source_index)
                        xpos_area_max=xpos_emission_subgrid+xpos_limit
                        xpos_area_min=xpos_emission_subgrid-xpos_limit
                        ypos_area_max=ypos_emission_subgrid+ypos_limit
                        ypos_area_min=ypos_emission_subgrid-ypos_limit

                        plume_count=0
                        plume_crossreference=0
                        plume_distance=0

                        !write(*,*) i_target_end-i_target_start,j_target_end-j_target_start
                        do j=j_target_start,j_target_end
                            do i=i_target_start,i_target_end

                                !Only calculate if it is within the local region distance
                                !Set the EMEP projection position of the emission grid. This guarantees that it extentds to the right distance
                                xpos_target_subgrid=xproj_emission_subgrid(i,j,source_index)
                                ypos_target_subgrid=yproj_emission_subgrid(i,j,source_index)

                                !Select only target grids within the predefined region
                                if (xpos_target_subgrid.ge.xpos_area_min.and.xpos_target_subgrid.le.xpos_area_max &
                                    .and.ypos_target_subgrid.ge.ypos_area_min.and.ypos_target_subgrid.le.ypos_area_max) then

                                    !Find the integral index for the target grid
                                    i_cross_target_integral=crossreference_emission_to_integral_subgrid(i,j,x_dim_index,source_index)
                                    j_cross_target_integral=crossreference_emission_to_integral_subgrid(i,j,y_dim_index,source_index)
                                    i_cross_target_integral=min(max(1,i_cross_target_integral),integral_subgrid_dim(x_dim_index))
                                    j_cross_target_integral=min(max(1,j_cross_target_integral),integral_subgrid_dim(y_dim_index))

                                    !Set the mixing height at the average of the emission and target position
                                    h_mix_loc=(meteo_subgrid(i_cross_integral,j_cross_integral,tt,hmix_subgrid_index)+meteo_subgrid(i_cross_target_integral,j_cross_target_integral,tt,hmix_subgrid_index))/2.

                                    !Set the precipitation at the receptor grid
                                    precip_loc=meteo_subgrid(i_cross_target_integral,j_cross_target_integral,tt,precip_subgrid_index)

                                    !Find the minimum distance to the trajectory and check it is valid (downwind)
                                    call uEMEP_minimum_distance_trajectory_fast(x_target_subgrid(i,j),y_target_subgrid(i,j), &
                                        traj_max_index,traj_step_size,trajectory_vector(:,x_dim_index),trajectory_vector(:,y_dim_index),x_loc,y_loc,valid_traj)

                                    !valid_traj=.false.

                                    if (valid_traj) then
                                        !Check if y_loc is within 4*sigma region based on the PG stability functions at emission position
                                        call uEMEP_set_dispersion_sigma_PG(invL_loc,logz0_loc,sig_z_00_loc,sig_y_00_loc,sigy_0_subgid_width_scale,emission_subgrid_delta(:,source_index),angle_diff(i_cross_integral,j_cross_integral),x_loc,sig_z_loc,sig_y_loc,sig_z_0_loc,sig_y_0_loc)

                                        !Put the plume in the plume vector for later use
                                        if (abs(y_loc).lt.3.*sig_y_loc) then
                                            plume_count=plume_count+1

                                            !Check number of subgrids in plume does not exceed allowable
                                            if (plume_count.gt.n_plume_subgrid_max) then
                                                write(unit_logfile,'(a,i)') 'WARNING: Number of subgrids for a plume exceeds the maximum value of: ', n_plume_subgrid_max
                                                write(unit_logfile,'(a)')   '         Exiting the plume but will continue calculating'
                                                goto 10
                                            else
                                                plume_crossreference(plume_count,x_dim_index)=i
                                                plume_crossreference(plume_count,y_dim_index)=j
                                                plume_distance(plume_count,x_dim_index)=x_loc
                                                plume_distance(plume_count,y_dim_index)=y_loc
                                            endif

                                        endif
                                    endif

                                endif
                            enddo
                        enddo

                        max_plume_count=max(max_plume_count,plume_count)

                        !Sort the target grids from closest to furthest creating a crossreference index
                        !write(*,*) 'IN: ',plume_distance(1:plume_count,x_dim_index)
10                      temp_plume_distance(1:plume_count)=plume_distance(1:plume_count,x_dim_index)
                        call rargsort(temp_plume_distance(1:plume_count),sorted_plume_index(1:plume_count),plume_count)
                        !write(*,*) 'OUT:',sorted_plume_index(1:plume_count)

                        !Set the plume source that will be depleted with deposition
                        plume_source=emission_subgrid(ii,jj,tt,source_index,:)

                        !Loop through the sorted array starting at the closest (which should always be the emission grid)
                        !write(*,*) ii,jj,plume_count
                        do k=1,plume_count

                            x_loc=plume_distance(sorted_plume_index(k),x_dim_index)
                            y_loc=plume_distance(sorted_plume_index(k),y_dim_index)
                            i=plume_crossreference(sorted_plume_index(k),x_dim_index)
                            j=plume_crossreference(sorted_plume_index(k),y_dim_index)

                            !write(*,*) x_loc,y_loc

                            !Select method for assigning sigma
                            if (stability_scheme_flag.eq.1) then
                                call uEMEP_set_dispersion_sigma_simple(sig_z_00_loc,sig_y_00_loc,sigy_0_subgid_width_scale,emission_subgrid_delta(:,source_index),angle_diff(i_cross_integral,j_cross_integral),x_loc,sig_z_loc,sig_y_loc,sig_z_0_loc,sig_y_0_loc)
                            endif

                            if (stability_scheme_flag.eq.2) then
                                call uEMEP_set_dispersion_params_PG(invL_loc,source_index,subsource_index)
                                ay_loc=ay(source_index,subsource_index)
                                by_loc=by(source_index,subsource_index)
                                az_loc=az(source_index,subsource_index)
                                bz_loc=bz(source_index,subsource_index)
                                call uEMEP_set_dispersion_sigma_PG(invL_loc,logz0_loc,sig_z_00_loc,sig_y_00_loc,sigy_0_subgid_width_scale,emission_subgrid_delta(:,source_index),angle_diff(i_cross_integral,j_cross_integral),x_loc,sig_z_loc,sig_y_loc,sig_z_0_loc,sig_y_0_loc)
                            endif

                            if (stability_scheme_flag.eq.3) then
                                !Set initial values for sigma. Initial sig_y is set here as well but is overridden by Kz dispersion
                                call uEMEP_set_dispersion_sigma_simple(sig_z_00_loc,sig_y_00_loc,sigy_0_subgid_width_scale,emission_subgrid_delta(:,source_index),angle_diff(i_cross_integral,j_cross_integral),x_loc,sig_z_loc,sig_y_loc,sig_z_0_loc,sig_y_0_loc)

                                call uEMEP_set_dispersion_sigma_Kz(Kz_scheme,x_loc,sig_z_00_loc,sig_y_00_loc,sigy_0_subgid_width_scale,sig_z_loc,h_emis_loc,h_mix_loc,invL_loc,FF10_loc,10.,logz0_loc,emission_subgrid_delta(:,source_index),u_star0_loc,average_zc_h_in_Kz_flag,n_kz_iterations,sig_y_scaling_factor,sig_z_loc,sig_y_loc,FF_zc_loc)

                                !Add the meandering and change in wind angle to the plume since not included in Kz calculation
                                sig_y_loc=sig_y_loc+x_loc*angle_diff(i_cross_integral,j_cross_integral)

                                !Use the average of the emision height and zc to determine wind speed. Is set to true if wind_level_flag=6
                                !CHECK THIS
                                if (wind_level_flag.eq.6.or.wind_level_zc_flag) then
                                    !Set the minimum wind speed
                                    FF_loc=sqrt(FF_zc_loc*FF_zc_loc+FF_min_dispersion*FF_min_dispersion)
                                endif
                                if (wind_level_integral_flag.eq.6) then
                                    !Set the minimum wind speed
                                    FF_integral_loc=sqrt(FF_zc_loc*FF_zc_loc+FF_min_dispersion*FF_min_dispersion)
                                endif

                            endif

                            if (stability_scheme_flag.eq.4) then
                                write(unit_logfile,'(a,i)') 'Stability_scheme_flag=4 (Kz emulator) no longer an option. Stopping'
                                stop
                                !call uEMEP_set_dispersion_sigma_Kz_emulator(h_emis_loc,invL_loc,logz0_loc,h_mix_loc,sig_z_00_loc,sig_y_00_loc,sigy_0_subgid_width_scale,emission_subgrid_delta(:,source_index),angle_diff(i_cross_integral,j_cross_integral),x_loc,sig_z_loc,sig_y_loc,sig_z_0_loc,sig_y_0_loc)
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
                            if (wind_level_integral_flag.eq.6.and.stability_scheme_flag.ne.3) then
                                !if (wind_level_flag.eq.6) then
                                call z_centremass_gauss_func(sig_z_loc,h_emis_loc,h_mix_loc,zc_loc)
                                zc_loc=(h_emis_loc+zc_loc)/2.
                                call u_profile_neutral_val_func(zc_loc,FF10_loc,10.,h_mix_loc,exp(logz0_loc),FF_zc_loc,u_star0_loc)
                                FF_integral_loc=sqrt(FF_zc_loc*FF_zc_loc+FF_min_dispersion*FF_min_dispersion)
                            endif

                            !Calculate the dispersion based on the derived sigmas
                            subgrid_internal=gauss_plume_cartesian_sigma_func(x_loc,y_loc,h_emis_loc,z_rec_loc,sig_z_loc,sig_y_loc,h_mix_loc,FF_loc)
                            !write(*,'(9es12.2)') x_loc,y_loc,h_emis_loc,z_rec_loc,sig_z_loc,sig_y_loc,h_mix_loc,FF_loc,subgrid_internal

                            !Only use half of the source grid for deposition and depletion
                            if (k.gt.1) then
                                !s/m3 *m2=s/m
                                subgrid_internal_integral=(subgrid_internal)*target_subgrid_delta(x_dim_index)*target_subgrid_delta(y_dim_index)
                            else
                                subgrid_internal_integral=(subgrid_internal)*target_subgrid_delta(x_dim_index)*0.5*target_subgrid_delta(y_dim_index) !Half the grid
                            endif

                            !Calculate the vertically integrated mass of the plume (s/m2) up to the lowest level and up to the mixing height
                            if (adjust_wetdepo_integral_to_lowest_layer_flag.and.calculate_deposition_flag) then
                                plume_vertical_integral(hsurf_integral_subgrid_index)=gauss_plume_cartesian_sigma_integral_func(x_loc,y_loc,h_emis_loc,sig_z_loc,sig_y_loc,h_mix_loc,FF_loc,0.,H_emep)*H_emep
                                plume_vertical_integral(hmix_integral_subgrid_index)=gauss_plume_cartesian_sigma_integral_func(x_loc,y_loc,h_emis_loc,sig_z_loc,sig_y_loc,h_mix_loc,FF_loc,0.,h_mix_loc)*h_mix_loc
                                !write(*,*) h_mix_loc,x_loc,sig_z_loc,vertical_integral_internal/(sig_z_loc*sqrt(2.*3.14159)/(2.*3.14159*sig_y_loc*sig_z_loc)*exp(-0.5*(y_loc*y_loc)/(sig_y_loc*sig_y_loc))/FF_loc*plume_source(i_pollutant) )
                                !write(*,*) plume_vertical_integral(1)/plume_vertical_integral(2),H_emep/h_mix_loc,H_emep/h_mix_loc/(plume_vertical_integral(1)/plume_vertical_integral(2))
                                !Calculate the average concentration in the lowest layer

                                !These two give the same results but the second is quicker. Can put it in a seperate subroutine
                                !vertical_integral_internal=gauss_plume_cartesian_sigma_integral_func(x_loc,y_loc,h_emis_loc,z_rec_loc,sig_z_loc,sig_y_loc,h_mix_loc,FF_loc,0.,h_mix_loc)*h_mix_loc*plume_source(i_pollutant)
                                vertical_integral_internal=exp(-0.5*(y_loc*y_loc)/(sig_y_loc*sig_y_loc))/FF_loc/(sqrt(2.*3.14159)*sig_y_loc)
                                !write(*,*) h_mix_loc,x_loc,sig_z_loc,vertical_integral_internal/(exp(-0.5*(y_loc*y_loc)/(sig_y_loc*sig_y_loc))/FF_loc/(sqrt(2.*3.14159)*sig_y_loc)*plume_source(i_pollutant))

                            endif
                            if (local_subgrid_method_flag.eq.1) then
                                plume_vertical_integral(hsurf_average_subgrid_index)=gauss_plume_cartesian_sigma_integral_func(x_loc,y_loc,h_emis_loc,sig_z_loc,sig_y_loc,h_mix_loc,FF_integral_loc,0.,H_emep)
                            endif


                            if (subgrid_internal_integral.gt.0) then
                                if (calculate_deposition_flag) then
                                    do i_pollutant=1,n_pollutant_loop

                                        !Multiply by the emissions to get the concentration prior to depletion (ug/m3)
                                        subgrid_internal_pollutant(i_pollutant)=subgrid_internal*plume_source(i_pollutant)

                                        !Calculate the dry deposition flux (ug/m3 *m/s=ug/s/m2)
                                        drydepo_internal_pollutant(i_pollutant)=subgrid_internal_pollutant(i_pollutant)*target_deposition_subgrid(i,j,vd_index,i_pollutant)

                                        !Calculate the dimensionless integral for plume depletion (s/m *m/s=1)
                                        drydepo_internal_integral=subgrid_internal_integral*target_deposition_subgrid(i,j,vd_index,i_pollutant)

                                        !Set the scavenging (s/m2 /m *m/s = /m2). 1e-3/3600 converts mm/hr to m/s
                                        wetdepo_internal=vertical_integral_internal*wetdepo_scavanging_rate(pollutant_loop_index(i_pollutant))*(precip_loc/1000./3600.)

                                        !Calculate the wet deposition flux ( /m2 ug/s = ug/m2/s)
                                        wetdepo_internal_pollutant(i_pollutant)=wetdepo_internal*plume_source(i_pollutant)

                                        !Calculate the dimensionless integral for plume depletion (/m2 * m2 = 1)
                                        wetdepo_internal_integral=wetdepo_internal*target_subgrid_delta(x_dim_index)*target_subgrid_delta(y_dim_index)

                                        !Make integrated concentration values
                                        plume_vertical_integral_pollutant(:)=plume_vertical_integral(:)*plume_source(i_pollutant)

                                        !Calculate the plume depletion by dry and wet deposition (ug/s)
                                        if (calculate_source_depletion_flag) then
                                            plume_source(i_pollutant)=plume_source(i_pollutant)*exp(-drydepo_internal_integral-wetdepo_internal_integral)
                                        endif

                                        !write(*,'(f12.2,4es12.2,f12.3)') x_loc,subgrid_internal_pollutant(i_pollutant),subgrid_internal,plume_source(i_pollutant),drydepo_internal_pollutant(i_pollutant)*target_subgrid_delta(x_dim_index)*target_subgrid_delta(y_dim_index),plume_source(i_pollutant)/emission_subgrid(ii,jj,tt,source_index,i_pollutant)

                                        !Add to the dry deposition target grid
                                        target_deposition_subgrid(i,j,drydepo_index,i_pollutant)=target_deposition_subgrid(i,j,drydepo_index,i_pollutant)+drydepo_internal_pollutant(i_pollutant)
                                        !Add to the wet deposition target grid
                                        target_deposition_subgrid(i,j,wetdepo_index,i_pollutant)=target_deposition_subgrid(i,j,wetdepo_index,i_pollutant)+wetdepo_internal_pollutant(i_pollutant)

                                        !Add to the target concentration subgrid position
                                        target_subgrid(i,j,i_pollutant)=target_subgrid(i,j,i_pollutant)+subgrid_internal_pollutant(i_pollutant)

                                    enddo

                                else
                                    !No deposition
                                    do i_pollutant=1,n_pollutant_loop

                                        !Multiply by the emissions to get the concentration prior to depletion (ug/m3)
                                        subgrid_internal_pollutant(i_pollutant)=subgrid_internal*plume_source(i_pollutant)

                                        !Add to the target concentration subgrid position
                                        target_subgrid(i,j,i_pollutant)=target_subgrid(i,j,i_pollutant)+subgrid_internal_pollutant(i_pollutant)

                                        !Make integrated concentration values
                                        plume_vertical_integral_pollutant(:)=plume_vertical_integral(:)*plume_source(i_pollutant)
                                    enddo

                                endif

                                !Determine the distance for the travel time calculation
                                distance_subgrid=sqrt(x_loc*x_loc+y_loc*y_loc)
                                distance_subgrid=max(distance_subgrid,distance_subgrid_min)

                                !Add to the travel time array
                                traveltime_target_subgrid(i,j,1,:)=traveltime_target_subgrid(i,j,1,:)+distance_subgrid/FF_loc*subgrid_internal_pollutant
                                traveltime_target_subgrid(i,j,2,:)=traveltime_target_subgrid(i,j,2,:)+subgrid_internal_pollutant

                                !Calculate the vertically integrated mass of the plume (s/m2) up to the lowest level and up to the mixing height
                                if ((calculate_deposition_flag.and.adjust_wetdepo_integral_to_lowest_layer_flag).or.local_subgrid_method_flag.eq.1) then
                                    target_vertical_integral_subgrid(i,j,:)=target_vertical_integral_subgrid(i,j,:)+plume_vertical_integral_pollutant(:)
                                    !write(*,'(2i,4es12.2)') i,j,plume_vertical_integral(1),target_vertical_integral_subgrid(i,j,1),plume_vertical_integral(2),target_vertical_integral_subgrid(i,j,2)
                                endif

                            endif

                        enddo !plume loop

                    endif!end if emission not 0

                enddo!End of emission loop
            enddo

            write(unit_logfile,'(a,i)') 'Maximum number of subgrids per plume in this period and region = ',max_plume_count

            !This ratio should be 1 in a well mixed layer. It can be used to estimate the average column concentration based on the surface layer concentration
            !target_vertical_integral_subgrid(:,:,4)=target_vertical_integral_subgrid(:,:,1)/target_vertical_integral_subgrid(:,:,2)

            !Interpolate target grid to the concentration subgrid
            do j=1,subgrid_dim(y_dim_index)
                do i=1,subgrid_dim(x_dim_index)
                    if (use_subgrid(i,j,source_index)) then
                        do i_pollutant=1,n_pollutant_loop
                            if (calculate_deposition_flag) then
                                subgrid(i,j,tt,drydepo_local_subgrid_index,source_index,i_pollutant)=subgrid(i,j,tt,drydepo_local_subgrid_index,source_index,i_pollutant) &
                                    +area_weighted_interpolation_function(x_target_subgrid,y_target_subgrid,target_deposition_subgrid(:,:,drydepo_index,i_pollutant) &
                                    ,target_subgrid_dim(x_dim_index),target_subgrid_dim(y_dim_index),target_subgrid_delta(:),x_subgrid(i,j),y_subgrid(i,j))
                                subgrid(i,j,tt,wetdepo_local_subgrid_index,source_index,i_pollutant)=subgrid(i,j,tt,wetdepo_local_subgrid_index,source_index,i_pollutant) &
                                    +area_weighted_interpolation_function(x_target_subgrid,y_target_subgrid,target_deposition_subgrid(:,:,wetdepo_index,i_pollutant) &
                                    ,target_subgrid_dim(x_dim_index),target_subgrid_dim(y_dim_index),target_subgrid_delta(:),x_subgrid(i,j),y_subgrid(i,j))
                            endif
                            subgrid(i,j,tt,proxy_subgrid_index,source_index,i_pollutant)=subgrid(i,j,tt,proxy_subgrid_index,source_index,i_pollutant) &
                                +area_weighted_interpolation_function(x_target_subgrid,y_target_subgrid,target_subgrid(:,:,i_pollutant) &
                                ,target_subgrid_dim(x_dim_index),target_subgrid_dim(y_dim_index),target_subgrid_delta(:),x_subgrid(i,j),y_subgrid(i,j))
                            traveltime_subgrid(i,j,tt,1,i_pollutant)=traveltime_subgrid(i,j,tt,1,i_pollutant) &
                                +area_weighted_interpolation_function(x_target_subgrid,y_target_subgrid,traveltime_target_subgrid(:,:,1,i_pollutant) &
                                ,target_subgrid_dim(x_dim_index),target_subgrid_dim(y_dim_index),target_subgrid_delta(:),x_subgrid(i,j),y_subgrid(i,j))
                            traveltime_subgrid(i,j,tt,2,i_pollutant)=traveltime_subgrid(i,j,tt,2,i_pollutant) &
                                +area_weighted_interpolation_function(x_target_subgrid,y_target_subgrid,traveltime_target_subgrid(:,:,2,i_pollutant) &
                                ,target_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),target_subgrid_delta(:),x_subgrid(i,j),y_subgrid(i,j))
                        enddo
                    else
                        subgrid(i,j,tt,proxy_subgrid_index,source_index,:)=NODATA_value
                        traveltime_subgrid(i,j,tt,:,:)=NODATA_value
                    endif
                enddo
            enddo

            !Determine the final travel time
            traveltime_subgrid(:,:,tt,3,:)=traveltime_subgrid(:,:,tt,1,:)/traveltime_subgrid(:,:,tt,2,:)
            where (traveltime_subgrid(:,:,tt,2,:).eq.0) traveltime_subgrid(:,:,tt,3,:)=3600.*12.

            !Place the vertically integrated values in the integral subgrid from the target grid, that is the same as the emission grid
            !write(*,*) 'Interpolating vertically integrated values to the integral grid'
            if ((calculate_deposition_flag.and.adjust_wetdepo_integral_to_lowest_layer_flag).or.local_subgrid_method_flag.eq.1) then
                do j=1,integral_subgrid_dim(y_dim_index)
                    do i=1,integral_subgrid_dim(x_dim_index)

                        !i_cross_integral=crossreference_integral_to_emission_subgrid(i,j,x_dim_index,source_index)
                        !j_cross_integral=crossreference_integral_to_emission_subgrid(i,j,y_dim_index,source_index)
                        !i_cross_integral=max(1,i_cross_integral);i_cross_integral=min(integral_subgrid_dim(x_dim_index),i_cross_integral)
                        !j_cross_integral=max(1,j_cross_integral);j_cross_integral=min(integral_subgrid_dim(y_dim_index),j_cross_integral)

                        do i_integral=1,n_integral_subgrid_index
                            do i_pollutant=1,n_pollutant_loop
                                integral_subgrid(i,j,tt,i_integral,source_index,i_pollutant)=integral_subgrid(i,j,tt,i_integral,source_index,i_pollutant) &
                                    +area_weighted_interpolation_function(x_target_subgrid,y_target_subgrid,target_vertical_integral_subgrid(:,:,i_integral) &
                                    ,target_subgrid_dim(x_dim_index),target_subgrid_dim(y_dim_index),target_subgrid_delta(:),x_integral_subgrid(i,j),y_integral_subgrid(i,j))
                                !write(*,*) i,j,i_integral,integral_subgrid(i,j,tt,i_integral,source_index,i_pollutant)
                            enddo
                        enddo
                    enddo
                enddo

                !Add to allsource
                do i_integral=1,n_integral_subgrid_index
                    integral_subgrid(:,:,tt,i_integral,allsource_index,:)=integral_subgrid(:,:,tt,i_integral,allsource_index,:)+integral_subgrid(:,:,tt,i_integral,source_index,:)
                enddo

            endif

            !Show mean outputs for checking
            do i_pollutant=1,n_pollutant_loop
                temp_sum_subgrid(i_pollutant)=0.
                count=0
                do j=1,subgrid_dim(y_dim_index)
                    do i=1,subgrid_dim(x_dim_index)
                        if (use_subgrid(i,j,source_index)) then
                            temp_sum_subgrid(i_pollutant)=temp_sum_subgrid(i_pollutant)+subgrid(i,j,tt,proxy_subgrid_index,source_index,i_pollutant)
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


        enddo!End time loop

        if (allocated(trajectory_vector)) deallocate(trajectory_vector)
        if (allocated(target_subgrid)) deallocate(target_subgrid)
        if (allocated(target_deposition_subgrid)) deallocate(target_deposition_subgrid)
        if (allocated(temp_FF_subgrid)) deallocate(temp_FF_subgrid)
        if (allocated(temp_FF_integral_subgrid)) deallocate(temp_FF_integral_subgrid)
        if (allocated(traveltime_target_subgrid)) deallocate(traveltime_target_subgrid)
        if (allocated(target_vertical_integral_subgrid)) deallocate(target_vertical_integral_subgrid)

    end subroutine uEMEP_subgrid_deposition


    subroutine uEMEP_create_wind_field(temp_FF_subgrid,angle_diff,wind_level_flag_in,source_index_in,subsource_index_in,tt_in)

        use uEMEP_definitions

        implicit none

        integer, intent(in) ::  wind_level_flag_in,tt_in,source_index_in,subsource_index_in
        integer j_cross,i_cross
        integer i_cross_integral,j_cross_integral
        integer ii,jj
        real z0_temp,h_temp
        real sig_y_00_loc,h_emis_loc,h_mix_loc,FF10_loc,x_loc,sig_z_loc,sig_y_loc,ff_loc,logz0_loc
        real sig_z_00_loc,sig_y_0_loc,sig_z_0_loc,zc_loc,u_star0_loc

        !define the temporary arrays for meteo
        real, intent(out) :: temp_FF_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index))
        real, intent(in) :: angle_diff(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index))

        temp_FF_subgrid=0.
        if (wind_level_flag_in.ne.5) then
            do j_cross=1,integral_subgrid_dim(y_dim_index)
                do i_cross=1,integral_subgrid_dim(x_dim_index)
                    z0_temp=exp(meteo_subgrid(i_cross,j_cross,tt_in,logz0_subgrid_index))
                    h_temp=h_emis(source_index_in,subsource_index_in)
                    if (annual_calculations.and.wind_level_flag_in.eq.1) then
                        temp_FF_subgrid(i_cross,j_cross)=1./meteo_subgrid(i_cross,j_cross,tt_in,inv_FFgrid_subgrid_index)
                    elseif (annual_calculations.and.wind_level_flag_in.eq.2) then
                        temp_FF_subgrid(i_cross,j_cross)=1./meteo_subgrid(i_cross,j_cross,tt_in,inv_FFgrid_subgrid_index)*(1.-(log((H_meteo+z0_temp)/z0_temp)-log((h_temp+z0_temp)/z0_temp))/log((H_meteo+z0_temp)/z0_temp))
                    elseif (annual_calculations.and.wind_level_flag_in.eq.3) then
                        temp_FF_subgrid(i_cross,j_cross)=1./meteo_subgrid(i_cross,j_cross,tt_in,inv_FF10_subgrid_index)
                    elseif (annual_calculations.and.wind_level_flag_in.eq.4) then
                        temp_FF_subgrid(i_cross,j_cross)=1./meteo_subgrid(i_cross,j_cross,tt_in,inv_FF10_subgrid_index)*(1.-(log((10.+z0_temp)/z0_temp)-log((h_temp+z0_temp)/z0_temp))/log((10.+z0_temp)/z0_temp))
                    elseif (hourly_calculations.and.wind_level_flag_in.eq.1) then
                        temp_FF_subgrid(i_cross,j_cross)=meteo_subgrid(i_cross,j_cross,tt_in,FFgrid_subgrid_index)
                    elseif (hourly_calculations.and.wind_level_flag_in.eq.2) then
                        temp_FF_subgrid(i_cross,j_cross)=meteo_subgrid(i_cross,j_cross,tt_in,FFgrid_subgrid_index)*(1.-(log((H_meteo+z0_temp)/z0_temp)-log((h_temp+z0_temp)/z0_temp))/log((H_meteo+z0_temp)/z0_temp))
                    elseif (hourly_calculations.and.wind_level_flag_in.eq.3) then
                        temp_FF_subgrid(i_cross,j_cross)=meteo_subgrid(i_cross,j_cross,tt_in,FF10_subgrid_index)
                    elseif (hourly_calculations.and.wind_level_flag_in.eq.4) then
                        temp_FF_subgrid(i_cross,j_cross)=meteo_subgrid(i_cross,j_cross,tt_in,FF10_subgrid_index)*(1.-(log((10.+z0_temp)/z0_temp)-log((h_temp+z0_temp)/z0_temp))/log((10.+z0_temp)/z0_temp))
                    elseif (wind_level_flag_in.eq.0) then
                        temp_FF_subgrid(i_cross,j_cross)=1.
                    elseif (wind_level_flag_in.eq.5) then
                        !Will set based on sigma z centre of mass
                        temp_FF_subgrid(i_cross,j_cross)=1.
                    elseif (wind_level_flag_in.eq.6) then
                        !Will set based on sigma z centre of mass and emission height
                        temp_FF_subgrid(i_cross,j_cross)=1.
                    else
                        write(unit_logfile,'(a)') 'No valid wind_level_flag_in selected. Stopping (uEMEP_subgrid_dispersion)'
                        stop
                    endif

                    !Setting a minimum value for wind for dispersion purposes (cannot be zero)
                    temp_FF_subgrid(i_cross,j_cross)=sqrt(temp_FF_subgrid(i_cross,j_cross)*temp_FF_subgrid(i_cross,j_cross)+FF_min_dispersion*FF_min_dispersion)

                    if (temp_FF_subgrid(i_cross,j_cross).eq.0) then
                        write(unit_logfile,'(a,2i)') 'Zero wind speed at integral grid (stopping): ',i_cross,j_cross
                        stop
                    endif

                enddo
            enddo
        endif


        !If wind level flag is set to 5, use of initial plume centre of mass, then set wind speed for each non-zero emission grid
        if (wind_level_flag_in.eq.5) then
            temp_FF_subgrid=0.
            do jj=1,emission_subgrid_dim(y_dim_index,source_index_in)
                do ii=1,emission_subgrid_dim(x_dim_index,source_index_in)
                    if (sum(emission_subgrid(ii,jj,tt_in,source_index_in,:)).ne.0) then

                        !Set the integral meteorological grid position for the emission position
                        i_cross_integral=crossreference_emission_to_integral_subgrid(ii,jj,x_dim_index,source_index_in)
                        j_cross_integral=crossreference_emission_to_integral_subgrid(ii,jj,y_dim_index,source_index_in)

                        !Set the local variables
                        logz0_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt_in,logz0_subgrid_index)
                        FF10_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt_in,FF10_subgrid_index)
                        sig_y_00_loc=emission_properties_subgrid(ii,jj,emission_sigy00_index,source_index_in)
                        sig_z_00_loc=emission_properties_subgrid(ii,jj,emission_sigz00_index,source_index_in)
                        h_emis_loc=emission_properties_subgrid(ii,jj,emission_h_index,source_index_in)
                        h_mix_loc=meteo_subgrid(i_cross_integral,j_cross_integral,tt_in,hmix_subgrid_index)

                        if (annual_calculations) then
                            FF10_loc=1./meteo_subgrid(i_cross_integral,j_cross_integral,tt_in,inv_FF10_subgrid_index)
                        endif

                        !Set sig_0's at the emission position
                        x_loc=0.
                        call uEMEP_set_dispersion_sigma_simple(sig_z_00_loc,sig_y_00_loc,sigy_0_subgid_width_scale,emission_subgrid_delta(:,source_index_in),angle_diff(i_cross_integral,j_cross_integral),x_loc,sig_z_loc,sig_y_loc,sig_z_0_loc,sig_y_0_loc)

                        !Use the initial plume centre of mass to determine wind advection height
                        call z_centremass_gauss_func(sig_z_0_loc,h_emis_loc,h_mix_loc,zc_loc)
                        call u_profile_neutral_val_func(zc_loc,FF10_loc,10.,h_mix_loc,exp(logz0_loc),FF_loc,u_star0_loc)

                        !Set the minimum wind speed
                        FF_loc=sqrt(FF_loc*FF_loc+FF_min_dispersion*FF_min_dispersion)

                        temp_FF_subgrid(ii,jj)=FF_loc
                        !write(*,*) FF10_loc,FF_loc,zc_loc,sig_z_0_loc

                    endif
                enddo
            enddo
        endif

    end subroutine uEMEP_create_wind_field

end module subgrid_deposition

