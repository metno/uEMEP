!uEMEP_define_subgrid.f90

    subroutine uEMEP_define_subgrid
    
    use uEMEP_definitions
    
    implicit none
    
    integer i_source
    !integer ii,jj
    real dx_temp,dy_temp
    real lon_temp,lat_temp

    write(unit_logfile,'(A)') ''
    write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Define subgrids and buffer zones (uEMEP_define_subgrid)'
	write(unit_logfile,'(A)') '================================================================'

    !If annual calculations then always set time start and stop to 1
    if (annual_calculations) then
        start_time_nc_index=1
        end_time_nc_index=1
    endif
    
    !Set the time index to be the same as the EMEP time dimensions
    emission_subgrid_dim(t_dim_index,:)=subgrid_dim(t_dim_index)
    integral_subgrid_dim(t_dim_index)=subgrid_dim(t_dim_index)
    
    write(unit_logfile,'(A,I5)')'Number of external time steps:',end_time_loop_index-start_time_loop_index+1
    write(unit_logfile,'(A,I5)')'Number of internal time steps:',subgrid_dim(t_dim_index)
    write(unit_logfile,'(A,2I5)')'Number of target grids:',subgrid_dim(1:2)
    write(unit_logfile,'(A,2I5)')'Number of integral grids:',integral_subgrid_dim(1:2)
    write(unit_logfile,'(A,2I5)')'Max number of emission grids:',emission_max_subgrid_dim(1:2)
    do i_source=1,n_source_index
    if (calculate_source(i_source)) then
    write(unit_logfile,'(A,A,2I5)')'Number of emission grids: ',trim(source_file_str(i_source)),emission_subgrid_dim(1:2,i_source)        
    endif
    enddo

    !Allocate buffers and adjust the dimensions appropriately
    !Calculate the max loop size to cover the nearest EMEP grids. This avoids looping through all the grids
    loop_index_scale=1.5 !Was 1.5
    
    !Define the centre of the subgrid
    !ii=int(subgrid_dim(x_dim_index)/2)
    !jj=int(subgrid_dim(y_dim_index)/2)
    if (projection_type.eq.RDM_projection_index) then
        call RDM2LL((subgrid_min(y_dim_index)+subgrid_max(y_dim_index))/2.,(subgrid_min(x_dim_index)+subgrid_max(x_dim_index))/2.,lat_temp,lon_temp)
    elseif (projection_type.eq.UTM_projection_index) then
        call UTM2LL(utm_zone,(subgrid_min(y_dim_index)+subgrid_max(y_dim_index))/2.,(subgrid_min(x_dim_index)+subgrid_max(x_dim_index))/2.,lat_temp,lon_temp)
    endif   

    if (EMEP_projection_type.eq.LL_projection_index) then
        dx_temp=111000.*dgrid_nc(lon_nc_index)*cos(lat_temp*pi/180.)/2.
        dy_temp=111000.*dgrid_nc(lat_nc_index)/2.
    else
        dx_temp=dgrid_nc(lon_nc_index)/2.
        dy_temp=dgrid_nc(lat_nc_index)/2.        
    endif
    
    
    subgrid_loop_index(x_dim_index)=floor(dx_temp/subgrid_delta(x_dim_index)*loop_index_scale*EMEP_grid_interpolation_size)
    subgrid_loop_index(y_dim_index)=floor(dy_temp/subgrid_delta(y_dim_index)*loop_index_scale*EMEP_grid_interpolation_size)

    integral_subgrid_loop_index(x_dim_index)=floor(dx_temp/integral_subgrid_delta(x_dim_index)*loop_index_scale*EMEP_grid_interpolation_size)
    integral_subgrid_loop_index(y_dim_index)=floor(dy_temp/integral_subgrid_delta(y_dim_index)*loop_index_scale*EMEP_grid_interpolation_size)

    do i_source=1,n_source_index
    if (calculate_source(i_source)) then
        emission_subgrid_loop_index(x_dim_index,i_source)=floor(dx_temp/emission_subgrid_delta(x_dim_index,i_source)*loop_index_scale*EMEP_grid_interpolation_size)
        emission_subgrid_loop_index(y_dim_index,i_source)=floor(dy_temp/emission_subgrid_delta(y_dim_index,i_source)*loop_index_scale*EMEP_grid_interpolation_size)
    endif
    enddo

    !Set the buffer sizes according to these loops for emissions only
    !This will remove edge effects for dispersion but will only remove edge effects for moving window only when emissions are used for redistribution
    buffer_index_scale=1.5
    if (use_downwind_position_flag) buffer_index_scale=buffer_index_scale*2.
    
    if (use_buffer_zone) then
        buffer_index(x_dim_index)=floor(dx_temp/subgrid_delta(x_dim_index)*buffer_index_scale)
        buffer_index(y_dim_index)=floor(dy_temp/subgrid_delta(y_dim_index)*buffer_index_scale)

        integral_buffer_index(x_dim_index)=floor(dx_temp/integral_subgrid_delta(x_dim_index)*buffer_index_scale)
        integral_buffer_index(y_dim_index)=floor(dy_temp/integral_subgrid_delta(y_dim_index)*buffer_index_scale)
        do i_source=1,n_source_index
        if (calculate_source(i_source)) then
            emission_buffer_index(x_dim_index,i_source)=floor(dx_temp/emission_subgrid_delta(x_dim_index,i_source)*buffer_index_scale)
            emission_buffer_index(y_dim_index,i_source)=floor(dy_temp/emission_subgrid_delta(y_dim_index,i_source)*buffer_index_scale)       
        endif
        enddo
        !buffer_index=subgrid_loop_index
        !emission_buffer_index=emission_subgrid_loop_index
        !integral_buffer_index=integral_subgrid_loop_index
    else
        buffer_index=0
        emission_buffer_index=0
        integral_buffer_index=0
    endif
    buffer_size=buffer_index*subgrid_delta
    emission_buffer_size=emission_buffer_index*emission_subgrid_delta
    integral_buffer_size=integral_buffer_index*integral_subgrid_delta
    
    do i_source=1,n_source_index
    if (calculate_source(i_source)) then
        emission_subgrid_dim(1:2,i_source)=emission_subgrid_dim(1:2,i_source)+emission_buffer_index(1:2,i_source)*2
        emission_subgrid_min(1:2,i_source)=emission_subgrid_min(1:2,i_source)-emission_buffer_size(1:2,i_source)
        emission_subgrid_max(1:2,i_source)=emission_subgrid_max(1:2,i_source)+emission_buffer_size(1:2,i_source)
    endif
    enddo

    emission_max_subgrid_dim(1:2)=emission_max_subgrid_dim(1:2)+buffer_index(1:2)*2
    integral_subgrid_dim(1:2)=integral_subgrid_dim(1:2)+integral_buffer_index(1:2)*2
    integral_subgrid_min(1:2)=integral_subgrid_min(1:2)-integral_buffer_size(1:2)
    integral_subgrid_max(1:2)=integral_subgrid_max(1:2)+integral_buffer_size(1:2)

    write(unit_logfile,'(A,2I5)')'Number of target grids to be looped for each EMEP grid:',subgrid_loop_index(1:2)
    write(unit_logfile,'(A,2I5)')'Number of integral grids to be looped for each EMEP grid:',integral_subgrid_loop_index(1:2)
    write(unit_logfile,'(A,2I5)')'Size of integral grid buffer zone:',integral_buffer_index(1:2)
    do i_source=1,n_source_index
    if (calculate_source(i_source)) then
    write(unit_logfile,'(A,A,2I5)')'Size of emission grid buffer zone: ',trim(source_file_str(i_source)),emission_buffer_index(1:2,i_source)        
    endif
    enddo

    write(unit_logfile,'(A,2I5)')'Number of target grids:',subgrid_dim(1:2)
    write(unit_logfile,'(A,2I5)')'Number of integral grids:',integral_subgrid_dim(1:2)
    write(unit_logfile,'(A,2I5)')'Max number of emission grids:',emission_max_subgrid_dim(1:2)
    do i_source=1,n_source_index
    if (calculate_source(i_source)) then
    write(unit_logfile,'(A,A,2I5)')'Number of emission grids: ',trim(source_file_str(i_source)),emission_subgrid_dim(1:2,i_source)        
    endif
    enddo

    !Define target grid
    allocate (subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_subgrid_index,n_source_index,n_possible_subsource))
    allocate (x_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index)))
    allocate (y_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index)))
    allocate (lon_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index)))
    allocate (lat_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index)))
    allocate (traveltime_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),2)) !Last index 1 for weighted time, 2 for sum of weights
    traveltime_subgrid=0.

    !Define compound subgrid. Same as target in dimensions
    allocate (comp_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index))
    allocate (comp_EMEP_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index))
    allocate (orig_EMEP_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index))

    do j=1,subgrid_dim(y_dim_index)
    do i=1,subgrid_dim(x_dim_index)                   
        x_subgrid(i,j)=subgrid_min(x_dim_index)+subgrid_delta(x_dim_index)*(i-1)
        y_subgrid(i,j)=subgrid_min(y_dim_index)+subgrid_delta(y_dim_index)*(j-1)
        if (projection_type.eq.RDM_projection_index) then
            call RDM2LL(y_subgrid(i,j),x_subgrid(i,j),lat_subgrid(i,j),lon_subgrid(i,j))
        elseif (projection_type.eq.UTM_projection_index) then
            call UTM2LL(utm_zone,y_subgrid(i,j),x_subgrid(i,j),lat_subgrid(i,j),lon_subgrid(i,j))
        endif   
    enddo
    enddo

    !Define integral grid
    allocate (integral_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),integral_subgrid_dim(t_dim_index),n_source_index,n_possible_subsource))
    allocate (x_integral_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index)))
    allocate (y_integral_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index)))
    allocate (lon_integral_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index)))
    allocate (lat_integral_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index)))

    do j=1,integral_subgrid_dim(y_dim_index)
    do i=1,integral_subgrid_dim(x_dim_index)                 
        x_integral_subgrid(i,j)=integral_subgrid_min(x_dim_index)+integral_subgrid_delta(x_dim_index)*(i-1)
        y_integral_subgrid(i,j)=integral_subgrid_min(y_dim_index)+integral_subgrid_delta(y_dim_index)*(j-1)
        if (projection_type.eq.RDM_projection_index) then
            call RDM2LL(y_integral_subgrid(i,j),x_integral_subgrid(i,j),lat_integral_subgrid(i,j),lon_integral_subgrid(i,j))
        elseif (projection_type.eq.UTM_projection_index) then
            call UTM2LL(utm_zone,y_integral_subgrid(i,j),x_integral_subgrid(i,j),lat_integral_subgrid(i,j),lon_integral_subgrid(i,j))
        endif   
    enddo
    enddo

    !Define meteo grid, same positional coordinates as the integral grid (lower resolution)
    allocate (meteo_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),integral_subgrid_dim(t_dim_index),n_meteo_subgrid_index))   

    !Alocate the use_subgrid array and set to true for all subgrids
    allocate (use_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),n_source_index))
    use_subgrid=.true.
    
    !Allocate emission grids with same dimensions as target grid
    allocate (proxy_emission_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),n_source_index,n_possible_subsource)) 
    allocate (emission_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_source_index,n_possible_subsource)) 
    allocate (x_emission_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),n_source_index))
    allocate (y_emission_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),n_source_index))
    allocate (lon_emission_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),n_source_index))
    allocate (lat_emission_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),n_source_index))
    allocate (emission_time_profile_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_source_index,n_possible_subsource)) 

    !Define emission grids
    do i_source=1,n_source_index
    if (calculate_source(i_source)) then
    do j=1,emission_subgrid_dim(y_dim_index,i_source)
    do i=1,emission_subgrid_dim(x_dim_index,i_source)                  
        x_emission_subgrid(i,j,i_source)=emission_subgrid_min(x_dim_index,i_source)+emission_subgrid_delta(x_dim_index,i_source)*(i-1)
        y_emission_subgrid(i,j,i_source)=emission_subgrid_min(y_dim_index,i_source)+emission_subgrid_delta(y_dim_index,i_source)*(j-1)
        if (projection_type.eq.UTM_projection_index) then
            call UTM2LL(utm_zone,y_emission_subgrid(i,j,i_source),x_emission_subgrid(i,j,i_source), &
            lat_emission_subgrid(i,j,i_source),lon_emission_subgrid(i,j,i_source))
        elseif (projection_type.eq.RDM_projection_index) then
            call RDM2LL(y_emission_subgrid(i,j,i_source),x_emission_subgrid(i,j,i_source), &
            lat_emission_subgrid(i,j,i_source),lon_emission_subgrid(i,j,i_source))
        endif
    enddo
    enddo
    endif
    enddo

    !Define population grid
    allocate (population_subgrid(population_subgrid_dim(x_dim_index),population_subgrid_dim(y_dim_index),n_population_index))
    allocate (x_population_subgrid(population_subgrid_dim(x_dim_index),population_subgrid_dim(y_dim_index)))
    allocate (y_population_subgrid(population_subgrid_dim(x_dim_index),population_subgrid_dim(y_dim_index)))
    allocate (lon_population_subgrid(population_subgrid_dim(x_dim_index),population_subgrid_dim(y_dim_index)))
    allocate (lat_population_subgrid(population_subgrid_dim(x_dim_index),population_subgrid_dim(y_dim_index)))

    do j=1,population_subgrid_dim(y_dim_index)
    do i=1,population_subgrid_dim(x_dim_index)                 
        x_population_subgrid(i,j)=population_subgrid_min(x_dim_index)+population_subgrid_delta(x_dim_index)*(i-1)
        y_population_subgrid(i,j)=population_subgrid_min(y_dim_index)+population_subgrid_delta(y_dim_index)*(j-1)
        if (projection_type.eq.RDM_projection_index) then
            call RDM2LL(y_population_subgrid(i,j),x_population_subgrid(i,j),lat_population_subgrid(i,j),lon_population_subgrid(i,j))
        elseif (projection_type.eq.UTM_projection_index) then
            call UTM2LL(utm_zone,y_population_subgrid(i,j),x_population_subgrid(i,j),lat_population_subgrid(i,j),lon_population_subgrid(i,j))
        endif   
    enddo
    enddo

    !Define exposure subgrid
    allocate (exposure_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_source_index))

    end subroutine uEMEP_define_subgrid