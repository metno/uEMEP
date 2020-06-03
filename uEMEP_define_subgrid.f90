!uEMEP_define_subgrid.f90
    subroutine uEMEP_define_subgrid_extent
    
    use uEMEP_definitions
    
    implicit none
    
    integer i,j
    integer i_source
    !integer ii,jj
    real dx_temp,dy_temp
    real lon_temp,lat_temp
    integer :: subsource_index=1

    write(unit_logfile,'(A)') ''
    write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Define subgrid extent and buffer zones (uEMEP_define_subgrid_extent)'
	write(unit_logfile,'(A)') '================================================================'

    
    !If annual calculations then always set time start and stop to 1
    if (annual_calculations) then
        start_time_nc_index=1
        end_time_nc_index=1
        start_time_meteo_nc_index=1
        end_time_meteo_nc_index=1
    endif
    
    !Set the time index to be the same as the EMEP time dimensions
    emission_subgrid_dim(t_dim_index,:)=subgrid_dim(t_dim_index)
    init_emission_subgrid_dim(t_dim_index,:)=subgrid_dim(t_dim_index)
    integral_subgrid_dim(t_dim_index)=subgrid_dim(t_dim_index)
    
    write(unit_logfile,'(A,I6)')'Number of external time steps:',end_time_loop_index-start_time_loop_index+1
    write(unit_logfile,'(A,I6)')'Number of internal time steps:',subgrid_dim(t_dim_index)
    write(unit_logfile,'(A,2I6)')'Number of target grids:',subgrid_dim(1:2)
    write(unit_logfile,'(A,2I6)')'Number of integral grids:',integral_subgrid_dim(1:2)
    write(unit_logfile,'(A,2I6)')'Max number of emission grids:',emission_max_subgrid_dim(1:2)
    do i_source=1,n_source_index
    if (calculate_source(i_source)) then
    write(unit_logfile,'(A,A,2I6)')'Number of emission grids: ',trim(source_file_str(i_source)),emission_subgrid_dim(1:2,i_source)        
    write(unit_logfile,'(A,A,2I8)')'Number of initial emission grids: ',trim(source_file_str(i_source)),init_emission_subgrid_dim(1:2,i_source)        
    endif
    enddo

    !Allocate buffers and adjust the dimensions appropriately
    !Calculate the max loop size to cover the nearest EMEP grids. This avoids looping through all the grids
    loop_index_scale=1.2*EMEP_grid_interpolation_size/2. !Was 1.5
    
    !Define the centre of the subgrid
    !ii=int(subgrid_dim(x_dim_index)/2)
    !jj=int(subgrid_dim(y_dim_index)/2)
    
    call PROJ2LL((subgrid_min(x_dim_index)+subgrid_max(x_dim_index))/2.,(subgrid_min(y_dim_index)+subgrid_max(y_dim_index))/2.,lon_temp,lat_temp,projection_attributes,projection_type)
    
    !if (projection_type.eq.RDM_projection_index) then
    !    call RDM2LL((subgrid_min(y_dim_index)+subgrid_max(y_dim_index))/2.,(subgrid_min(x_dim_index)+subgrid_max(x_dim_index))/2.,lat_temp,lon_temp)
    !elseif (projection_type.eq.UTM_projection_index) then
    !    call UTM2LL(utm_zone,(subgrid_min(y_dim_index)+subgrid_max(y_dim_index))/2.,(subgrid_min(x_dim_index)+subgrid_max(x_dim_index))/2.,lat_temp,lon_temp)
    !endif   

    if (EMEP_projection_type.eq.LL_projection_index) then
        dx_temp=111000.*dgrid_nc(lon_nc_index)*cos(lat_temp*pi/180.)
        dy_temp=111000.*dgrid_nc(lat_nc_index)
    else
        !Assumed LCC
        dx_temp=dgrid_nc(lon_nc_index)
        dy_temp=dgrid_nc(lat_nc_index)
    endif
    
    subgrid_loop_index(x_dim_index)=floor(dx_temp/subgrid_delta(x_dim_index)*loop_index_scale)
    subgrid_loop_index(y_dim_index)=floor(dy_temp/subgrid_delta(y_dim_index)*loop_index_scale)

    integral_subgrid_loop_index(x_dim_index)=floor(dx_temp/integral_subgrid_delta(x_dim_index)*loop_index_scale)
    integral_subgrid_loop_index(y_dim_index)=floor(dy_temp/integral_subgrid_delta(y_dim_index)*loop_index_scale)

    do i_source=1,n_source_index
    if (calculate_source(i_source)) then
        emission_subgrid_loop_index(x_dim_index,i_source)=floor(dx_temp/emission_subgrid_delta(x_dim_index,i_source)*loop_index_scale)
        emission_subgrid_loop_index(y_dim_index,i_source)=floor(dy_temp/emission_subgrid_delta(y_dim_index,i_source)*loop_index_scale)
    endif
    enddo

    do i_source=1,n_source_index
    if (calculate_source(i_source)) then
        init_emission_subgrid_loop_index(x_dim_index,i_source)=floor(dx_temp/init_emission_subgrid_delta(x_dim_index,i_source)*loop_index_scale)
        init_emission_subgrid_loop_index(y_dim_index,i_source)=floor(dy_temp/init_emission_subgrid_delta(y_dim_index,i_source)*loop_index_scale)
    endif
    enddo

    if (calculate_deposition_flag) then
    deposition_subgrid_loop_index(x_dim_index)=floor(dx_temp/deposition_subgrid_delta(x_dim_index)*loop_index_scale)
    deposition_subgrid_loop_index(y_dim_index)=floor(dy_temp/deposition_subgrid_delta(y_dim_index)*loop_index_scale)
    endif
    if (read_landuse_flag) then
    landuse_subgrid_loop_index(x_dim_index)=floor(dx_temp/landuse_subgrid_delta(x_dim_index)*loop_index_scale)
    landuse_subgrid_loop_index(y_dim_index)=floor(dy_temp/landuse_subgrid_delta(y_dim_index)*loop_index_scale)
    endif
    
    !Set the buffer sizes according to these loops for emissions only
    !This will remove edge effects for dispersion but will only remove edge effects for moving window only when emissions are used for redistribution
    buffer_index_scale=loop_index_scale
    !if (use_downwind_position_flag) buffer_index_scale=buffer_index_scale*2.
    
    if (use_buffer_zone) then
        buffer_index(x_dim_index)=floor(dx_temp/subgrid_delta(x_dim_index)*buffer_index_scale)
        buffer_index(y_dim_index)=floor(dy_temp/subgrid_delta(y_dim_index)*buffer_index_scale)

        integral_buffer_index(x_dim_index)=floor(dx_temp/integral_subgrid_delta(x_dim_index)*buffer_index_scale)
        integral_buffer_index(y_dim_index)=floor(dy_temp/integral_subgrid_delta(y_dim_index)*buffer_index_scale)
        
        !Buffer index for emissions a half grid larger because of the possible use of the integral
        do i_source=1,n_source_index
        if (calculate_source(i_source)) then
            emission_buffer_index(x_dim_index,i_source)=floor(dx_temp/emission_subgrid_delta(x_dim_index,i_source)*(buffer_index_scale+0.5))
            emission_buffer_index(y_dim_index,i_source)=floor(dy_temp/emission_subgrid_delta(y_dim_index,i_source)*(buffer_index_scale+0.5))       
            if (local_subgrid_method_flag.eq.3) then
                !Extend the emission buffer zone 0.5 EMEP grids when redistributing emissions to gurantee all emissions are accounted for
                emission_buffer_index(x_dim_index,i_source)=floor(dx_temp/emission_subgrid_delta(x_dim_index,i_source)*(buffer_index_scale+1.0))
                emission_buffer_index(y_dim_index,i_source)=floor(dy_temp/emission_subgrid_delta(y_dim_index,i_source)*(buffer_index_scale+1.0))
            endif
        endif
        enddo

        do i_source=1,n_source_index
        if (calculate_source(i_source)) then
            init_emission_buffer_index(x_dim_index,i_source)=floor(dx_temp/init_emission_subgrid_delta(x_dim_index,i_source)*(buffer_index_scale+0.5))
            init_emission_buffer_index(y_dim_index,i_source)=floor(dy_temp/init_emission_subgrid_delta(y_dim_index,i_source)*(buffer_index_scale+0.5))       
            if (local_subgrid_method_flag.eq.3) then
                !Extend the emission buffer zone 0.5 EMEP grids when redistributing emissions to gurantee all emissions are accounted for
                init_emission_buffer_index(x_dim_index,i_source)=floor(dx_temp/init_emission_subgrid_delta(x_dim_index,i_source)*(buffer_index_scale+1.0))
                init_emission_buffer_index(y_dim_index,i_source)=floor(dy_temp/init_emission_subgrid_delta(y_dim_index,i_source)*(buffer_index_scale+1.0))       
            endif
        endif
        enddo

        if (calculate_deposition_flag) then
        deposition_buffer_index(x_dim_index)=floor(dx_temp/deposition_subgrid_delta(x_dim_index)*buffer_index_scale)
        deposition_buffer_index(y_dim_index)=floor(dy_temp/deposition_subgrid_delta(y_dim_index)*buffer_index_scale)
        endif
        if (read_landuse_flag) then
        landuse_buffer_index(x_dim_index)=floor(dx_temp/landuse_subgrid_delta(x_dim_index)*buffer_index_scale)
        landuse_buffer_index(y_dim_index)=floor(dy_temp/landuse_subgrid_delta(y_dim_index)*buffer_index_scale)
        endif
    else
        buffer_index=0
        emission_buffer_index=0
        init_emission_buffer_index=0
        integral_buffer_index=0
        deposition_buffer_index=0
        landuse_buffer_index=0
    endif
    buffer_size=buffer_index*subgrid_delta
    emission_buffer_size=emission_buffer_index*emission_subgrid_delta
    init_emission_buffer_size=emission_buffer_index*init_emission_subgrid_delta
    integral_buffer_size=integral_buffer_index*integral_subgrid_delta
    if (calculate_deposition_flag) then
        deposition_buffer_size=deposition_buffer_index*deposition_subgrid_delta
    endif
    if (read_landuse_flag) then
        landuse_buffer_size=landuse_buffer_index*landuse_subgrid_delta
    endif
    
    do i_source=1,n_source_index
    if (calculate_source(i_source)) then
        emission_subgrid_dim(1:2,i_source)=emission_subgrid_dim(1:2,i_source)+emission_buffer_index(1:2,i_source)*2
        emission_subgrid_min(1:2,i_source)=emission_subgrid_min(1:2,i_source)-emission_buffer_size(1:2,i_source)
        emission_subgrid_max(1:2,i_source)=emission_subgrid_max(1:2,i_source)+emission_buffer_size(1:2,i_source)
    endif
    enddo
    !emission_max_subgrid_dim(1:2)=emission_max_subgrid_dim(1:2)+buffer_index(1:2)*2
    emission_max_subgrid_dim(1:2)=maxval(emission_subgrid_dim(1:2,:),2)
    
    do i_source=1,n_source_index
    if (calculate_source(i_source)) then
        init_emission_subgrid_dim(1:2,i_source)=init_emission_subgrid_dim(1:2,i_source)+init_emission_buffer_index(1:2,i_source)*2
        init_emission_subgrid_min(1:2,i_source)=init_emission_subgrid_min(1:2,i_source)-init_emission_buffer_size(1:2,i_source)
        init_emission_subgrid_max(1:2,i_source)=init_emission_subgrid_max(1:2,i_source)+init_emission_buffer_size(1:2,i_source)
    endif
    enddo

    integral_subgrid_dim(1:2)=integral_subgrid_dim(1:2)+integral_buffer_index(1:2)*2
    integral_subgrid_min(1:2)=integral_subgrid_min(1:2)-integral_buffer_size(1:2)
    integral_subgrid_max(1:2)=integral_subgrid_max(1:2)+integral_buffer_size(1:2)

    if (calculate_deposition_flag) then
    deposition_subgrid_dim(1:2)=deposition_subgrid_dim(1:2)+deposition_buffer_index(1:2)*2
    deposition_subgrid_min(1:2)=deposition_subgrid_min(1:2)-deposition_buffer_size(1:2)
    deposition_subgrid_max(1:2)=deposition_subgrid_max(1:2)+deposition_buffer_size(1:2)
    endif
    if (read_landuse_flag) then
    landuse_subgrid_dim(1:2)=landuse_subgrid_dim(1:2)+landuse_buffer_index(1:2)*2
    landuse_subgrid_min(1:2)=landuse_subgrid_min(1:2)-landuse_buffer_size(1:2)
    landuse_subgrid_max(1:2)=landuse_subgrid_max(1:2)+landuse_buffer_size(1:2)
    endif
    
    write(unit_logfile,'(A,2I5)')'Number of target grids to be looped for each EMEP grid:',subgrid_loop_index(1:2)
    write(unit_logfile,'(A,2I5)')'Number of integral grids to be looped for each EMEP grid:',integral_subgrid_loop_index(1:2)
    write(unit_logfile,'(A,2I5)')'Size of integral grid buffer zone:',integral_buffer_index(1:2)
    if (calculate_deposition_flag) then
    write(unit_logfile,'(A,2I5)')'Size of deposition grid buffer zone:',deposition_buffer_index(1:2)
    endif
    if (read_landuse_flag) then
    write(unit_logfile,'(A,2I5)')'Size of landuse grid buffer zone:',landuse_buffer_index(1:2)
    endif
    
    do i_source=1,n_source_index
    if (calculate_source(i_source)) then
    write(unit_logfile,'(A,A,2I5)')'Size of emission grid buffer zone: ',trim(source_file_str(i_source)),emission_buffer_index(1:2,i_source)        
    write(unit_logfile,'(A,A,2I5)')'Size of initial emission grid buffer zone: ',trim(source_file_str(i_source)),emission_buffer_index(1:2,i_source)        
    endif
    enddo

    write(unit_logfile,'(A,2I6)')'Number of target grids:',subgrid_dim(1:2)
    write(unit_logfile,'(A,2I6)')'Number of integral grids:',integral_subgrid_dim(1:2)
    write(unit_logfile,'(A,2I6)')'Number of integral grids:',integral_subgrid_dim(1:2)
    if (calculate_deposition_flag) then
    write(unit_logfile,'(A,2I6)')'Number of deposition grids:',deposition_subgrid_dim(1:2)
    endif
    write(unit_logfile,'(A,2I6)')'Max number of emission grids:',emission_max_subgrid_dim(1:2)
    do i_source=1,n_source_index
    if (calculate_source(i_source)) then
    write(unit_logfile,'(A,A,2I6)')   'Number of emission grids:',trim(source_file_str(i_source)),emission_subgrid_dim(1:2,i_source)        
    write(unit_logfile,'(A,A,2f12.1)')'Min of emission grids:   ',trim(source_file_str(i_source)),emission_subgrid_min(1:2,i_source)     
    write(unit_logfile,'(A,A,2f12.1)')'Max of emission grids:   ',trim(source_file_str(i_source)),emission_subgrid_max(1:2,i_source)     
    write(unit_logfile,'(A,A,2f12.1)')'Delta of emission grids: ',trim(source_file_str(i_source)),emission_subgrid_delta(1:2,i_source)     
    endif
    if (calculate_source(i_source)) then
    write(unit_logfile,'(A,A,2I8)')   'Number of initial emission grids:',trim(source_file_str(i_source)),init_emission_subgrid_dim(1:2,i_source)        
    write(unit_logfile,'(A,A,2f12.1)')'Min of initial emission grids:   ',trim(source_file_str(i_source)),init_emission_subgrid_min(1:2,i_source)     
    write(unit_logfile,'(A,A,2f12.1)')'Max of initial emission grids:   ',trim(source_file_str(i_source)),init_emission_subgrid_max(1:2,i_source)     
    write(unit_logfile,'(A,A,2f12.1)')'Delta of initial emission grids: ',trim(source_file_str(i_source)),init_emission_subgrid_delta(1:2,i_source)     
    endif
    enddo

    end subroutine uEMEP_define_subgrid_extent
 
    subroutine uEMEP_define_subgrid
    
    use uEMEP_definitions
    
    implicit none
    
    integer i,j
    integer i_source
    !integer ii,jj
    real dx_temp,dy_temp
    real lon_temp,lat_temp
    integer :: subsource_index=1

    write(unit_logfile,'(A)') ''
    write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Define subgrid arrays (uEMEP_define_subgrid)'
	write(unit_logfile,'(A)') '================================================================'


    !Deallocate grids if they are already allocated. This will be in the case of the use_multiple_receptor_grids_flag=.true.
    if (allocated(subgrid)) deallocate (subgrid)
    if (allocated(x_subgrid)) deallocate (x_subgrid)
    if (allocated(y_subgrid)) deallocate (y_subgrid)
    if (allocated(lon_subgrid)) deallocate (lon_subgrid)
    if (allocated(lat_subgrid)) deallocate (lat_subgrid)
    if (allocated(xproj_subgrid)) deallocate (xproj_subgrid)
    if (allocated(yproj_subgrid)) deallocate (yproj_subgrid)
    if (allocated(traveltime_subgrid)) deallocate (traveltime_subgrid)
    
    !Define target grid
    if (.not.allocated(subgrid)) allocate (subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_subgrid_index,n_source_index,n_pollutant_loop))
    if (.not.allocated(x_subgrid)) allocate (x_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index)))
    if (.not.allocated(y_subgrid)) allocate (y_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index)))
    if (.not.allocated(lon_subgrid)) allocate (lon_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index)))
    if (.not.allocated(lat_subgrid)) allocate (lat_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index)))
    if (.not.allocated(xproj_subgrid)) allocate (xproj_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index)))
    if (.not.allocated(yproj_subgrid)) allocate (yproj_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index)))
    if (.not.allocated(traveltime_subgrid)) allocate (traveltime_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),2,n_pollutant_loop)) !Last index 1 for weighted time, 2 for sum of weights
    traveltime_subgrid=0.

    !Deallocate grids if they are already allocated. This will be in the case of the use_multiple_receptor_grids_flag=.true.
    if (allocated(comp_subgrid)) deallocate (comp_subgrid)
    if (allocated(comp_EMEP_subgrid)) deallocate (comp_EMEP_subgrid)
    if (allocated(orig_EMEP_subgrid)) deallocate (orig_EMEP_subgrid)
    if (allocated(species_EMEP_subgrid)) deallocate (species_EMEP_subgrid)
    
    !Define compound subgrid. Same as target in dimensions
    if (.not.allocated(comp_subgrid)) then
        allocate (comp_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index))
        comp_subgrid=0.
    endif
    if (.not.allocated(comp_EMEP_subgrid)) then
        allocate (comp_EMEP_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index))
        comp_EMEP_subgrid=0.
    endif
    if (.not.allocated(orig_EMEP_subgrid)) then
        allocate (orig_EMEP_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index))
        orig_EMEP_subgrid=0.
    endif
    if (.not.allocated(species_EMEP_subgrid).and.(save_emep_species.or.save_seasalt)) then
        allocate (species_EMEP_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_pmxx_sp_index,n_species_loop_index))
        species_EMEP_subgrid=0.
    endif

    
    
    do j=1,subgrid_dim(y_dim_index)
    do i=1,subgrid_dim(x_dim_index)                   
        x_subgrid(i,j)=subgrid_min(x_dim_index)+subgrid_delta(x_dim_index)*(i-0.5)
        y_subgrid(i,j)=subgrid_min(y_dim_index)+subgrid_delta(y_dim_index)*(j-0.5)
        
        call PROJ2LL(x_subgrid(i,j),y_subgrid(i,j),lon_subgrid(i,j),lat_subgrid(i,j),projection_attributes,projection_type)

        !if (projection_type.eq.RDM_projection_index) then
        !    call RDM2LL(y_subgrid(i,j),x_subgrid(i,j),lat_subgrid(i,j),lon_subgrid(i,j))
        !elseif (projection_type.eq.UTM_projection_index) then
        !    call UTM2LL(utm_zone,y_subgrid(i,j),x_subgrid(i,j),lat_subgrid(i,j),lon_subgrid(i,j))
        !endif
        
        !If the EMEP projection is lambert then set the proj coordinates to lambert, otherwise to lat-lon
        if (EMEP_projection_type.eq.LCC_projection_index) then
            call lb2lambert2_uEMEP(xproj_subgrid(i,j),yproj_subgrid(i,j),lon_subgrid(i,j),lat_subgrid(i,j),EMEP_projection_attributes)
        else
            xproj_subgrid(i,j)=lon_subgrid(i,j)
            yproj_subgrid(i,j)=lat_subgrid(i,j)            
        endif
    enddo
    enddo

    !Deallocate grids if they are already allocated. This will be in the case of the use_multiple_receptor_grids_flag=.true.
    if (allocated(integral_subgrid)) deallocate (integral_subgrid)
    if (allocated(x_integral_subgrid)) deallocate (x_integral_subgrid)
    if (allocated(y_integral_subgrid)) deallocate (y_integral_subgrid)
    if (allocated(lon_integral_subgrid)) deallocate (lon_integral_subgrid)
    if (allocated(lat_integral_subgrid)) deallocate (lat_integral_subgrid)
    if (allocated(xproj_integral_subgrid)) deallocate (xproj_integral_subgrid)
    if (allocated(yproj_integral_subgrid)) deallocate (yproj_integral_subgrid)
    if (allocated(meteo_nc_xproj_integral_subgrid)) deallocate (meteo_nc_xproj_integral_subgrid)
    if (allocated(meteo_nc_yproj_integral_subgrid)) deallocate (meteo_nc_yproj_integral_subgrid)

    !Define integral grid
    if (.not.allocated(integral_subgrid)) allocate (integral_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),integral_subgrid_dim(t_dim_index),n_integral_subgrid_index,n_source_index,n_pollutant_loop))
    if (.not.allocated(x_integral_subgrid)) allocate (x_integral_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index)))
    if (.not.allocated(y_integral_subgrid)) allocate (y_integral_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index)))
    if (.not.allocated(lon_integral_subgrid)) allocate (lon_integral_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index)))
    if (.not.allocated(lat_integral_subgrid)) allocate (lat_integral_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index)))
    if (.not.allocated(xproj_integral_subgrid)) allocate (xproj_integral_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index)))
    if (.not.allocated(yproj_integral_subgrid)) allocate (yproj_integral_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index)))
    if (use_alternative_meteorology_flag) then
        if (.not.allocated(meteo_nc_xproj_integral_subgrid)) allocate (meteo_nc_xproj_integral_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index)))
        if (.not.allocated(meteo_nc_yproj_integral_subgrid)) allocate (meteo_nc_yproj_integral_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index)))
    endif
    integral_subgrid=0.
    
    do j=1,integral_subgrid_dim(y_dim_index)
    do i=1,integral_subgrid_dim(x_dim_index)                 
        x_integral_subgrid(i,j)=integral_subgrid_min(x_dim_index)+integral_subgrid_delta(x_dim_index)*(i-0.5)
        y_integral_subgrid(i,j)=integral_subgrid_min(y_dim_index)+integral_subgrid_delta(y_dim_index)*(j-0.5)

        call PROJ2LL(x_integral_subgrid(i,j),y_integral_subgrid(i,j),lon_integral_subgrid(i,j),lat_integral_subgrid(i,j),projection_attributes,projection_type)
        !if (projection_type.eq.RDM_projection_index) then
        !    call RDM2LL(y_integral_subgrid(i,j),x_integral_subgrid(i,j),lat_integral_subgrid(i,j),lon_integral_subgrid(i,j))
        !elseif (projection_type.eq.UTM_projection_index) then
        !    call UTM2LL(utm_zone,y_integral_subgrid(i,j),x_integral_subgrid(i,j),lat_integral_subgrid(i,j),lon_integral_subgrid(i,j))
        !endif 
        
        !If the EMEP projection is lambert then set the proj coordinates to lambert, otherwise to lat-lon
        if (EMEP_projection_type.eq.LCC_projection_index) then
            call lb2lambert2_uEMEP(xproj_integral_subgrid(i,j),yproj_integral_subgrid(i,j),lon_integral_subgrid(i,j),lat_integral_subgrid(i,j),EMEP_projection_attributes)
        else
            xproj_integral_subgrid(i,j)=lon_integral_subgrid(i,j)
            yproj_integral_subgrid(i,j)=lat_integral_subgrid(i,j)            
        endif
        if (use_alternative_meteorology_flag) then
        if (meteo_nc_projection_type.eq.LCC_projection_index) then
            call lb2lambert2_uEMEP(meteo_nc_xproj_integral_subgrid(i,j),meteo_nc_yproj_integral_subgrid(i,j),lon_integral_subgrid(i,j),lat_integral_subgrid(i,j),meteo_nc_projection_attributes)
        else
            meteo_nc_xproj_integral_subgrid(i,j)=lon_integral_subgrid(i,j)
            meteo_nc_yproj_integral_subgrid(i,j)=lat_integral_subgrid(i,j)            
        endif
        endif
   enddo
    enddo

    !Deallocate grids if they are already allocated. This will be in the case of the use_multiple_receptor_grids_flag=.true.
    if (allocated(meteo_subgrid)) deallocate (meteo_subgrid)   
    if (allocated(last_meteo_subgrid)) deallocate (last_meteo_subgrid)   
    !Define meteo grid, same positional coordinates as the integral grid (lower resolution)
    if (.not.allocated(meteo_subgrid)) allocate (meteo_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),integral_subgrid_dim(t_dim_index),n_meteo_subgrid_index))   
    if (.not.allocated(last_meteo_subgrid)) allocate (last_meteo_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),n_meteo_subgrid_index))   

    !Deallocate grids if they are already allocated. This will be in the case of the use_multiple_receptor_grids_flag=.true.
    if (allocated(use_subgrid)) deallocate (use_subgrid)
    if (allocated(use_subgrid_val)) deallocate (use_subgrid_val)
    if (allocated(use_subgrid_interpolation_index)) deallocate (use_subgrid_interpolation_index)
    
    !Allocate the use_subgrid array and set to true for all subgrids
    if (.not.allocated(use_subgrid)) allocate (use_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),n_source_index))
    if (.not.allocated(use_subgrid_val)) allocate (use_subgrid_val(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),n_source_index))
    if (.not.allocated(use_subgrid_interpolation_index)) allocate (use_subgrid_interpolation_index(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),n_source_index))
    
    use_subgrid=.true.
    use_subgrid_val=1
    
    !Deallocate grids if they are already allocated. This will be in the case of the use_multiple_receptor_grids_flag=.true.
    if (allocated(proxy_emission_subgrid)) deallocate (proxy_emission_subgrid) 
    if (allocated(emission_subgrid)) deallocate (emission_subgrid) 
    if (allocated(x_emission_subgrid)) deallocate (x_emission_subgrid)
    if (allocated(y_emission_subgrid)) deallocate (y_emission_subgrid)
    if (allocated(lon_emission_subgrid)) deallocate (lon_emission_subgrid)
    if (allocated(lat_emission_subgrid)) deallocate (lat_emission_subgrid)
    if (allocated(xproj_emission_subgrid)) deallocate (xproj_emission_subgrid)
    if (allocated(yproj_emission_subgrid)) deallocate (yproj_emission_subgrid)
    if (allocated(emission_time_profile_subgrid)) deallocate (emission_time_profile_subgrid) 
    if (allocated(emission_properties_subgrid)) deallocate (emission_properties_subgrid) 
 
    !Allocate emission grids with same dimensions as the maximum emission subgrid
    if (.not.allocated(proxy_emission_subgrid)) allocate (proxy_emission_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),n_source_index,n_pollutant_loop)) 
    if (.not.allocated(emission_subgrid)) allocate (emission_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_source_index,n_pollutant_loop)) 
    if (.not.allocated(x_emission_subgrid)) allocate (x_emission_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),n_source_index))
    if (.not.allocated(y_emission_subgrid)) allocate (y_emission_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),n_source_index))
    if (.not.allocated(lon_emission_subgrid)) allocate (lon_emission_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),n_source_index))
    if (.not.allocated(lat_emission_subgrid)) allocate (lat_emission_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),n_source_index))
    if (.not.allocated(emission_time_profile_subgrid)) allocate (emission_time_profile_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_source_index,n_pollutant_loop)) 
    if (.not.allocated(xproj_emission_subgrid)) allocate (xproj_emission_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),n_source_index))
    if (.not.allocated(yproj_emission_subgrid)) allocate (yproj_emission_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),n_source_index))
    if (.not.allocated(emission_properties_subgrid)) allocate (emission_properties_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),n_emission_index,n_source_index)) 
    emission_time_profile_subgrid=1.
    proxy_emission_subgrid=0.
    emission_subgrid=0.
    emission_properties_subgrid=0.
    
    !Define emission grids
    do i_source=1,n_source_index
    if (calculate_source(i_source)) then
    do j=1,emission_subgrid_dim(y_dim_index,i_source)
    do i=1,emission_subgrid_dim(x_dim_index,i_source)
        
        x_emission_subgrid(i,j,i_source)=emission_subgrid_min(x_dim_index,i_source)+emission_subgrid_delta(x_dim_index,i_source)*(i-0.5)
        y_emission_subgrid(i,j,i_source)=emission_subgrid_min(y_dim_index,i_source)+emission_subgrid_delta(y_dim_index,i_source)*(j-0.5)

        call PROJ2LL(x_emission_subgrid(i,j,i_source),y_emission_subgrid(i,j,i_source),lon_emission_subgrid(i,j,i_source),lat_emission_subgrid(i,j,i_source),projection_attributes,projection_type)
        !if (projection_type.eq.UTM_projection_index) then
        !    call UTM2LL(utm_zone,y_emission_subgrid(i,j,i_source),x_emission_subgrid(i,j,i_source), &
        !    lat_emission_subgrid(i,j,i_source),lon_emission_subgrid(i,j,i_source))
        !elseif (projection_type.eq.RDM_projection_index) then
        !    call RDM2LL(y_emission_subgrid(i,j,i_source),x_emission_subgrid(i,j,i_source), &
        !    lat_emission_subgrid(i,j,i_source),lon_emission_subgrid(i,j,i_source))
        !endif
        
        !If the EMEP projection is lambert then set the proj coordinates to lambert, otherwise to lat-lon
        if (EMEP_projection_type.eq.LCC_projection_index) then
            call lb2lambert2_uEMEP(xproj_emission_subgrid(i,j,i_source),yproj_emission_subgrid(i,j,i_source),lon_emission_subgrid(i,j,i_source),lat_emission_subgrid(i,j,i_source),EMEP_projection_attributes)
        else
            xproj_emission_subgrid(i,j,i_source)=lon_emission_subgrid(i,j,i_source)
            yproj_emission_subgrid(i,j,i_source)=lat_emission_subgrid(i,j,i_source)            
        endif

    enddo
    enddo
    endif
    enddo

    !Deallocate grids if they are already allocated. This will be in the case of the use_multiple_receptor_grids_flag=.true.
    if (allocated(population_subgrid)) deallocate (population_subgrid)
    if (allocated(x_population_subgrid)) deallocate (x_population_subgrid)
    if (allocated(y_population_subgrid)) deallocate (y_population_subgrid)
    if (allocated(lon_population_subgrid)) deallocate (lon_population_subgrid)
    if (allocated(lat_population_subgrid)) deallocate (lat_population_subgrid)
    if (allocated(xproj_population_subgrid)) deallocate (xproj_population_subgrid)
    if (allocated(yproj_population_subgrid)) deallocate (yproj_population_subgrid)
    
    !Define population grid
    if (.not.allocated(population_subgrid)) allocate (population_subgrid(population_subgrid_dim(x_dim_index),population_subgrid_dim(y_dim_index),n_population_index))
    if (.not.allocated(x_population_subgrid)) allocate (x_population_subgrid(population_subgrid_dim(x_dim_index),population_subgrid_dim(y_dim_index)))
    if (.not.allocated(y_population_subgrid)) allocate (y_population_subgrid(population_subgrid_dim(x_dim_index),population_subgrid_dim(y_dim_index)))
    if (.not.allocated(lon_population_subgrid)) allocate (lon_population_subgrid(population_subgrid_dim(x_dim_index),population_subgrid_dim(y_dim_index)))
    if (.not.allocated(lat_population_subgrid)) allocate (lat_population_subgrid(population_subgrid_dim(x_dim_index),population_subgrid_dim(y_dim_index)))
    if (.not.allocated(xproj_population_subgrid)) allocate (xproj_population_subgrid(population_subgrid_dim(x_dim_index),population_subgrid_dim(y_dim_index)))
    if (.not.allocated(yproj_population_subgrid)) allocate (yproj_population_subgrid(population_subgrid_dim(x_dim_index),population_subgrid_dim(y_dim_index)))
    population_subgrid=0
    
    do j=1,population_subgrid_dim(y_dim_index)
    do i=1,population_subgrid_dim(x_dim_index)                 

        x_population_subgrid(i,j)=population_subgrid_min(x_dim_index)+population_subgrid_delta(x_dim_index)*(i-0.5)
        y_population_subgrid(i,j)=population_subgrid_min(y_dim_index)+population_subgrid_delta(y_dim_index)*(j-0.5)
        !Set the lat-lon coordinates of the population
        
        call PROJ2LL(x_population_subgrid(i,j),y_population_subgrid(i,j),lon_population_subgrid(i,j),lat_population_subgrid(i,j),projection_attributes,projection_type)
        !if (projection_type.eq.RDM_projection_index) then
        !    call RDM2LL(y_population_subgrid(i,j),x_population_subgrid(i,j),lat_population_subgrid(i,j),lon_population_subgrid(i,j))
        !elseif (projection_type.eq.UTM_projection_index) then
        !    call UTM2LL(utm_zone,y_population_subgrid(i,j),x_population_subgrid(i,j),lat_population_subgrid(i,j),lon_population_subgrid(i,j))
        !endif
        
        !If the EMEP projection is lambert then set the proj coordinates to lambert, otherwise to lat-lon
        if (EMEP_projection_type.eq.LCC_projection_index) then
            call lb2lambert2_uEMEP(xproj_population_subgrid(i,j),yproj_population_subgrid(i,j),lon_population_subgrid(i,j),lat_population_subgrid(i,j),EMEP_projection_attributes)
        else
            xproj_population_subgrid(i,j)=lon_population_subgrid(i,j)
            yproj_population_subgrid(i,j)=lat_population_subgrid(i,j)            
        endif

    enddo
    enddo

    !Place some properties in the emission properties subgrid
    do j=1,emission_max_subgrid_dim(y_dim_index)
    do i=1,emission_max_subgrid_dim(x_dim_index)
    
        emission_properties_subgrid(i,j,emission_h_index,:)=h_emis(:,subsource_index)
        emission_properties_subgrid(i,j,emission_sigz00_index,:)=sig_z_00(:,subsource_index)
        emission_properties_subgrid(i,j,emission_sigy00_index,:)=sig_y_00(:,subsource_index)
    
    enddo
    enddo

    !Deallocate grids if they are already allocated. This will be in the case of the use_multiple_receptor_grids_flag=.true.
    if (calculate_deposition_flag) then
    if (allocated(orig_EMEP_deposition_subgrid)) deallocate (orig_EMEP_deposition_subgrid)
    if (allocated(deposition_subgrid)) deallocate (deposition_subgrid)
    if (allocated(x_deposition_subgrid)) deallocate (x_deposition_subgrid)
    if (allocated(y_deposition_subgrid)) deallocate (y_deposition_subgrid)
    if (allocated(lon_deposition_subgrid)) deallocate (lon_deposition_subgrid)
    if (allocated(lat_deposition_subgrid)) deallocate (lat_deposition_subgrid)
    if (allocated(xproj_deposition_subgrid)) deallocate (xproj_deposition_subgrid)
    if (allocated(yproj_deposition_subgrid)) deallocate (yproj_deposition_subgrid)

    !Define deposition grid
    if (.not.allocated(orig_EMEP_deposition_subgrid)) allocate (orig_EMEP_deposition_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_deposition_index,n_compound_index))
    if (.not.allocated(deposition_subgrid)) allocate (deposition_subgrid(deposition_subgrid_dim(x_dim_index),deposition_subgrid_dim(y_dim_index),deposition_subgrid_dim(t_dim_index),n_deposition_index,n_pollutant_loop))
    if (.not.allocated(x_deposition_subgrid)) allocate (x_deposition_subgrid(deposition_subgrid_dim(x_dim_index),deposition_subgrid_dim(y_dim_index)))
    if (.not.allocated(y_deposition_subgrid)) allocate (y_deposition_subgrid(deposition_subgrid_dim(x_dim_index),deposition_subgrid_dim(y_dim_index)))
    if (.not.allocated(lon_deposition_subgrid)) allocate (lon_deposition_subgrid(deposition_subgrid_dim(x_dim_index),deposition_subgrid_dim(y_dim_index)))
    if (.not.allocated(lat_deposition_subgrid)) allocate (lat_deposition_subgrid(deposition_subgrid_dim(x_dim_index),deposition_subgrid_dim(y_dim_index)))
    if (.not.allocated(xproj_deposition_subgrid)) allocate (xproj_deposition_subgrid(deposition_subgrid_dim(x_dim_index),deposition_subgrid_dim(y_dim_index)))
    if (.not.allocated(yproj_deposition_subgrid)) allocate (yproj_deposition_subgrid(deposition_subgrid_dim(x_dim_index),deposition_subgrid_dim(y_dim_index)))
    orig_EMEP_deposition_subgrid=0
    deposition_subgrid=0
    
    do j=1,deposition_subgrid_dim(y_dim_index)
    do i=1,deposition_subgrid_dim(x_dim_index)                 

        x_deposition_subgrid(i,j)=deposition_subgrid_min(x_dim_index)+deposition_subgrid_delta(x_dim_index)*(i-0.5)
        y_deposition_subgrid(i,j)=deposition_subgrid_min(y_dim_index)+deposition_subgrid_delta(y_dim_index)*(j-0.5)
        !Set the lat-lon coordinates of the deposition

        call PROJ2LL(x_deposition_subgrid(i,j),y_deposition_subgrid(i,j),lon_deposition_subgrid(i,j),lat_deposition_subgrid(i,j),projection_attributes,projection_type)
        !if (projection_type.eq.RDM_projection_index) then
        !    call RDM2LL(y_deposition_subgrid(i,j),x_deposition_subgrid(i,j),lat_deposition_subgrid(i,j),lon_deposition_subgrid(i,j))
        !elseif (projection_type.eq.UTM_projection_index) then
        !    call UTM2LL(utm_zone,y_deposition_subgrid(i,j),x_deposition_subgrid(i,j),lat_deposition_subgrid(i,j),lon_deposition_subgrid(i,j))
        !endif
        
        !If the EMEP projection is lambert then set the proj coordinates to lambert, otherwise to lat-lon
        if (EMEP_projection_type.eq.LCC_projection_index) then
            call lb2lambert2_uEMEP(xproj_deposition_subgrid(i,j),yproj_deposition_subgrid(i,j),lon_deposition_subgrid(i,j),lat_deposition_subgrid(i,j),EMEP_projection_attributes)
        else
            xproj_deposition_subgrid(i,j)=lon_deposition_subgrid(i,j)
            yproj_deposition_subgrid(i,j)=lat_deposition_subgrid(i,j)            
        endif

    enddo
    enddo
    endif
    
    if (read_landuse_flag) then
    !Deallocate grids if they are already allocated. This will be in the case of the use_multiple_receptor_grids_flag=.true.
    if (allocated(landuse_subgrid)) deallocate (landuse_subgrid)
    if (allocated(x_landuse_subgrid)) deallocate (x_landuse_subgrid)
    if (allocated(y_landuse_subgrid)) deallocate (y_landuse_subgrid)
    if (allocated(lon_landuse_subgrid)) deallocate (lon_landuse_subgrid)
    if (allocated(lat_landuse_subgrid)) deallocate (lat_landuse_subgrid)
    if (allocated(xproj_landuse_subgrid)) deallocate (xproj_landuse_subgrid)
    if (allocated(yproj_landuse_subgrid)) deallocate (yproj_landuse_subgrid)
    
    !Define landuse grid
    if (.not.allocated(landuse_subgrid)) allocate (landuse_subgrid(landuse_subgrid_dim(x_dim_index),landuse_subgrid_dim(y_dim_index),n_landuse_index))
    if (.not.allocated(x_landuse_subgrid)) allocate (x_landuse_subgrid(landuse_subgrid_dim(x_dim_index),landuse_subgrid_dim(y_dim_index)))
    if (.not.allocated(y_landuse_subgrid)) allocate (y_landuse_subgrid(landuse_subgrid_dim(x_dim_index),landuse_subgrid_dim(y_dim_index)))
    if (.not.allocated(lon_landuse_subgrid)) allocate (lon_landuse_subgrid(landuse_subgrid_dim(x_dim_index),landuse_subgrid_dim(y_dim_index)))
    if (.not.allocated(lat_landuse_subgrid)) allocate (lat_landuse_subgrid(landuse_subgrid_dim(x_dim_index),landuse_subgrid_dim(y_dim_index)))
    if (.not.allocated(xproj_landuse_subgrid)) allocate (xproj_landuse_subgrid(landuse_subgrid_dim(x_dim_index),landuse_subgrid_dim(y_dim_index)))
    if (.not.allocated(yproj_landuse_subgrid)) allocate (yproj_landuse_subgrid(landuse_subgrid_dim(x_dim_index),landuse_subgrid_dim(y_dim_index)))
    landuse_subgrid=0
    
    do j=1,landuse_subgrid_dim(y_dim_index)
    do i=1,landuse_subgrid_dim(x_dim_index)                 

        x_landuse_subgrid(i,j)=landuse_subgrid_min(x_dim_index)+landuse_subgrid_delta(x_dim_index)*(i-0.5)
        y_landuse_subgrid(i,j)=landuse_subgrid_min(y_dim_index)+landuse_subgrid_delta(y_dim_index)*(j-0.5)
        !Set the lat-lon coordinates of the landuse

        call PROJ2LL(x_landuse_subgrid(i,j),y_landuse_subgrid(i,j),lon_landuse_subgrid(i,j),lat_landuse_subgrid(i,j),projection_attributes,projection_type)
        !if (projection_type.eq.RDM_projection_index) then
        !    call RDM2LL(y_landuse_subgrid(i,j),x_landuse_subgrid(i,j),lat_landuse_subgrid(i,j),lon_landuse_subgrid(i,j))
        !elseif (projection_type.eq.UTM_projection_index) then
        !    call UTM2LL(utm_zone,y_landuse_subgrid(i,j),x_landuse_subgrid(i,j),lat_landuse_subgrid(i,j),lon_landuse_subgrid(i,j))
        !endif

        !If the EMEP projection is lambert then set the proj coordinates to lambert, otherwise to lat-lon
        if (EMEP_projection_type.eq.LCC_projection_index) then
            call lb2lambert2_uEMEP(xproj_landuse_subgrid(i,j),yproj_landuse_subgrid(i,j),lon_landuse_subgrid(i,j),lat_landuse_subgrid(i,j),EMEP_projection_attributes)
        else
            xproj_landuse_subgrid(i,j)=lon_landuse_subgrid(i,j)
            yproj_landuse_subgrid(i,j)=lat_landuse_subgrid(i,j)            
        endif

    enddo
    enddo
    endif
    
    !Deallocate grids if they are already allocated. This will be in the case of the use_multiple_receptor_grids_flag=.true.
    if (allocated(exposure_subgrid)) deallocate (exposure_subgrid)
    !Define exposure subgrid
    if (.not.allocated(exposure_subgrid)) allocate (exposure_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_source_index,n_pollutant_loop))

    end subroutine uEMEP_define_subgrid