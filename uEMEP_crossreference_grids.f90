!uEMEP_crossreference_grids
    
    subroutine uEMEP_crossreference_grids

    use uEMEP_definitions
    
    implicit none
    
    integer i,j,k
    integer ii,jj
    integer i_source
    real x_temp,y_temp
    
    !Cross referencing must be done for each new grid when using multiple grids
    if (allocated(crossreference_target_to_emep_subgrid)) deallocate (crossreference_target_to_emep_subgrid)
    if (allocated(crossreference_integral_to_emep_subgrid)) deallocate (crossreference_integral_to_emep_subgrid)
    if (allocated(crossreference_target_to_integral_subgrid)) deallocate (crossreference_target_to_integral_subgrid)
    if (allocated(crossreference_target_to_emission_subgrid)) deallocate (crossreference_target_to_emission_subgrid)
    if (allocated(crossreference_emission_to_EMEP_subgrid)) deallocate (crossreference_emission_to_EMEP_subgrid)
    if (allocated(crossreference_integral_to_emission_subgrid)) deallocate (crossreference_integral_to_emission_subgrid)
    if (allocated(crossreference_emission_to_integral_subgrid)) deallocate (crossreference_emission_to_integral_subgrid)
    if (allocated(crossreference_target_to_population_subgrid)) deallocate (crossreference_target_to_population_subgrid)
    if (use_alternative_meteorology_flag) then
        if (allocated(crossreference_integral_to_meteo_nc_subgrid)) deallocate (crossreference_integral_to_meteo_nc_subgrid)
    endif
    
    allocate (crossreference_target_to_emep_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),2))
    allocate (crossreference_integral_to_emep_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),2))
    allocate (crossreference_target_to_integral_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),2))
    allocate (crossreference_target_to_emission_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),2,n_source_index))
    allocate (crossreference_emission_to_EMEP_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),2,n_source_index))
    allocate (crossreference_integral_to_emission_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),2,n_source_index))
    allocate (crossreference_emission_to_integral_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),2,n_source_index))
    allocate (crossreference_target_to_population_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),2))

    if (use_alternative_meteorology_flag) then
        allocate (crossreference_integral_to_meteo_nc_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),2))
    endif
    
    write(unit_logfile,'(A)')'Allocating EMEP grid index to subgrid index'
    !Loop through subgrid and find those subgrids within EMEP grids and allocate concentrations directly from EMEP grids.         
    do j=1,subgrid_dim(y_dim_index)
    do i=1,subgrid_dim(x_dim_index)
        if (EMEP_projection_type.eq.LL_projection_index) then
            ii=1+floor((lon_subgrid(i,j)-var1d_nc(1,lon_nc_index))/dgrid_nc(lon_nc_index)+0.5)
            jj=1+floor((lat_subgrid(i,j)-var1d_nc(1,lat_nc_index))/dgrid_nc(lat_nc_index)+0.5)     
            crossreference_target_to_emep_subgrid(i,j,x_dim_index)=ii
            crossreference_target_to_emep_subgrid(i,j,y_dim_index)=jj
        elseif (EMEP_projection_type.eq.LCC_projection_index) then
            !When EMEP is read as x,y projection then var1d_nc(:,lon/lat_nc_index) are the x, y projection indexes, actually
            !if (use_alternative_LCC_projection_flag) then
                call lb2lambert2_uEMEP(x_temp,y_temp,lon_subgrid(i,j),lat_subgrid(i,j),EMEP_projection_attributes)
            !else
            !    call lb2lambert_uEMEP(x_temp,y_temp,lon_subgrid(i,j),lat_subgrid(i,j),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
            !endif
            ii=1+floor((x_temp-var1d_nc(1,lon_nc_index))/dgrid_nc(lon_nc_index)+0.5)
            jj=1+floor((y_temp-var1d_nc(1,lat_nc_index))/dgrid_nc(lat_nc_index)+0.5)     
            crossreference_target_to_emep_subgrid(i,j,x_dim_index)=ii
            crossreference_target_to_emep_subgrid(i,j,y_dim_index)=jj
        else
            write(unit_logfile,'(A)')'No valid projection in use. Stopping'
            stop
        endif   
    enddo
    enddo
    write(unit_logfile,'(A)')'Allocating integral grid index to subgrid index'
    do j=1,subgrid_dim(y_dim_index)
    do i=1,subgrid_dim(x_dim_index)
        crossreference_target_to_integral_subgrid(i,j,x_dim_index)=1+floor((x_subgrid(i,j)-integral_subgrid_min(x_dim_index))/integral_subgrid_delta(x_dim_index))
        crossreference_target_to_integral_subgrid(i,j,y_dim_index)=1+floor((y_subgrid(i,j)-integral_subgrid_min(y_dim_index))/integral_subgrid_delta(y_dim_index))      
    enddo
    enddo
    write(unit_logfile,'(A)')'Allocating population grid index to subgrid index'
    do j=1,subgrid_dim(y_dim_index)
    do i=1,subgrid_dim(x_dim_index)
        crossreference_target_to_population_subgrid(i,j,x_dim_index)=1+floor((x_subgrid(i,j)-population_subgrid_min(x_dim_index))/population_subgrid_delta(x_dim_index))
        crossreference_target_to_population_subgrid(i,j,y_dim_index)=1+floor((y_subgrid(i,j)-population_subgrid_min(y_dim_index))/population_subgrid_delta(y_dim_index))      
    enddo
    enddo
    write(unit_logfile,'(A)')'Allocating EMEP grid index to integral subgrid index'
    do j=1,integral_subgrid_dim(y_dim_index)
    do i=1,integral_subgrid_dim(x_dim_index)
        if (EMEP_projection_type.eq.LL_projection_index) then
            ii=1+floor((lon_integral_subgrid(i,j)-var1d_nc(1,lon_nc_index))/dgrid_nc(lon_nc_index)+0.5)
            jj=1+floor((lat_integral_subgrid(i,j)-var1d_nc(1,lat_nc_index))/dgrid_nc(lat_nc_index)+0.5)     
            crossreference_integral_to_emep_subgrid(i,j,x_dim_index)=ii
            crossreference_integral_to_emep_subgrid(i,j,y_dim_index)=jj
        elseif (EMEP_projection_type.eq.LCC_projection_index) then
            !When EMEP is read as x,y projection then var1d_nc(:,lon/lat_nc_index) are the x, y projection indexes, actually
            !if (use_alternative_LCC_projection_flag) then
                call lb2lambert2_uEMEP(x_temp,y_temp,lon_integral_subgrid(i,j),lat_integral_subgrid(i,j),EMEP_projection_attributes)
            !else
            !    call lb2lambert_uEMEP(x_temp,y_temp,lon_integral_subgrid(i,j),lat_integral_subgrid(i,j),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
            !endif
            !call lb2lambert_uEMEP(x_temp,y_temp,lon_integral_subgrid(i,j),lat_integral_subgrid(i,j),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
            ii=1+floor((x_temp-var1d_nc(1,lon_nc_index))/dgrid_nc(lon_nc_index)+0.5)
            jj=1+floor((y_temp-var1d_nc(1,lat_nc_index))/dgrid_nc(lat_nc_index)+0.5)     
            crossreference_integral_to_emep_subgrid(i,j,x_dim_index)=ii
            crossreference_integral_to_emep_subgrid(i,j,y_dim_index)=jj            
        else
            write(unit_logfile,'(A)')'No valid projection in use. Stopping'
            stop
        endif
    enddo
    enddo

    if (use_alternative_meteorology_flag) then
    write(unit_logfile,'(A)')'Allocating alternative meteo nc grid index to integral subgrid index'
    do j=1,integral_subgrid_dim(y_dim_index)
    do i=1,integral_subgrid_dim(x_dim_index)
        if (meteo_nc_projection_type.eq.LL_projection_index) then
            ii=1+floor((lon_integral_subgrid(i,j)-meteo_var1d_nc(1,lon_nc_index))/meteo_dgrid_nc(lon_nc_index)+0.5)
            jj=1+floor((lat_integral_subgrid(i,j)-meteo_var1d_nc(1,lat_nc_index))/meteo_dgrid_nc(lat_nc_index)+0.5)     
            crossreference_integral_to_meteo_nc_subgrid(i,j,x_dim_index)=ii
            crossreference_integral_to_meteo_nc_subgrid(i,j,y_dim_index)=jj
        elseif (meteo_nc_projection_type.eq.LCC_projection_index) then
            !When EMEP is read as x,y projection then var1d_nc(:,lon/lat_nc_index) are the x, y projection indexes, actually
            !if (use_alternative_LCC_projection_flag) then
                call lb2lambert2_uEMEP(x_temp,y_temp,lon_integral_subgrid(i,j),lat_integral_subgrid(i,j),meteo_nc_projection_attributes)
            !else
            !    call lb2lambert_uEMEP(x_temp,y_temp,lon_integral_subgrid(i,j),lat_integral_subgrid(i,j),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
            !endif
            !call lb2lambert_uEMEP(x_temp,y_temp,lon_integral_subgrid(i,j),lat_integral_subgrid(i,j),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
            ii=1+floor((x_temp-meteo_var1d_nc(1,lon_nc_index))/meteo_dgrid_nc(lon_nc_index)+0.5)
            jj=1+floor((y_temp-meteo_var1d_nc(1,lat_nc_index))/meteo_dgrid_nc(lat_nc_index)+0.5)     
            crossreference_integral_to_meteo_nc_subgrid(i,j,x_dim_index)=ii
            crossreference_integral_to_meteo_nc_subgrid(i,j,y_dim_index)=jj            
        else
            write(unit_logfile,'(A)')'No valid projection in use. Stopping'
            stop
        endif
    enddo
    enddo
    endif
    
    do i_source=1,n_source_index
    if (calculate_source(i_source)) then
    !do j=1,emission_subgrid_dim(y_dim_index,i_source)
    !do i=1,emission_subgrid_dim(x_dim_index,i_source)
    !    crossreference_emission_to_target_subgrid(i,j,x_dim_index,i_source)=1+floor((x_emission_subgrid(i,j,i_source)-subgrid_min(x_dim_index))/subgrid_delta(x_dim_index)+0.5)
    !    crossreference_emission_to_target_subgrid(i,j,y_dim_index,i_source)=1+floor((y_emission_subgrid(i,j,i_source)-subgrid_min(y_dim_index))/subgrid_delta(y_dim_index)+0.5)      
    !enddo
    !enddo
    do j=1,emission_subgrid_dim(y_dim_index,i_source)
    do i=1,emission_subgrid_dim(x_dim_index,i_source)
        if (EMEP_projection_type.eq.LL_projection_index) then
            ii=1+floor((lon_emission_subgrid(i,j,i_source)-var1d_nc(1,lon_nc_index))/dgrid_nc(lon_nc_index)+0.5)
            jj=1+floor((lat_emission_subgrid(i,j,i_source)-var1d_nc(1,lat_nc_index))/dgrid_nc(lat_nc_index)+0.5)     
            crossreference_emission_to_emep_subgrid(i,j,x_dim_index,i_source)=ii
            crossreference_emission_to_emep_subgrid(i,j,y_dim_index,i_source)=jj
        elseif (EMEP_projection_type.eq.LCC_projection_index) then
            !When EMEP is read as x,y projection then var1d_nc(:,lon/lat_nc_index) are the x, y projection indexes, actually
            !if (use_alternative_LCC_projection_flag) then
                call lb2lambert2_uEMEP(x_temp,y_temp,lon_emission_subgrid(i,j,i_source),lat_emission_subgrid(i,j,i_source),EMEP_projection_attributes)
            !else
            !    call lb2lambert_uEMEP(x_temp,y_temp,lon_emission_subgrid(i,j,i_source),lat_emission_subgrid(i,j,i_source),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
            !endif
            !call lb2lambert_uEMEP(x_temp,y_temp,lon_emission_subgrid(i,j,i_source),lat_emission_subgrid(i,j,i_source),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
            ii=1+floor((x_temp-var1d_nc(1,lon_nc_index))/dgrid_nc(lon_nc_index)+0.5)
            jj=1+floor((y_temp-var1d_nc(1,lat_nc_index))/dgrid_nc(lat_nc_index)+0.5)
            !write(*,'(4i12,2es12.2)') i,j,ii,jj,x_temp,y_temp
            crossreference_emission_to_emep_subgrid(i,j,x_dim_index,i_source)=ii
            crossreference_emission_to_emep_subgrid(i,j,y_dim_index,i_source)=jj
        else
            write(unit_logfile,'(A)')'No valid projection in use. Stopping'
            stop
        endif
    enddo
    enddo
    do j=1,subgrid_dim(y_dim_index)
    do i=1,subgrid_dim(x_dim_index)
        crossreference_target_to_emission_subgrid(i,j,x_dim_index,i_source)=1+floor((x_subgrid(i,j)-emission_subgrid_min(x_dim_index,i_source))/emission_subgrid_delta(x_dim_index,i_source))
        crossreference_target_to_emission_subgrid(i,j,y_dim_index,i_source)=1+floor((y_subgrid(i,j)-emission_subgrid_min(y_dim_index,i_source))/emission_subgrid_delta(y_dim_index,i_source))      
    enddo
    enddo
    do j=1,emission_subgrid_dim(y_dim_index,i_source)
    do i=1,emission_subgrid_dim(x_dim_index,i_source)
        crossreference_emission_to_integral_subgrid(i,j,x_dim_index,i_source)=1+floor((x_emission_subgrid(i,j,i_source)-integral_subgrid_min(x_dim_index))/integral_subgrid_delta(x_dim_index))
        crossreference_emission_to_integral_subgrid(i,j,y_dim_index,i_source)=1+floor((y_emission_subgrid(i,j,i_source)-integral_subgrid_min(y_dim_index))/integral_subgrid_delta(y_dim_index))
        
        !At edge this can return negative distances due to the different sizes of emission and integral grids and buffer zones. Set the limits here. Should not be a problem 
        crossreference_emission_to_integral_subgrid(i,j,x_dim_index,i_source)=max(min(crossreference_emission_to_integral_subgrid(i,j,x_dim_index,i_source),integral_subgrid_dim(x_dim_index)),1)
        crossreference_emission_to_integral_subgrid(i,j,y_dim_index,i_source)=max(min(crossreference_emission_to_integral_subgrid(i,j,y_dim_index,i_source),integral_subgrid_dim(y_dim_index)),1)
        
        if (crossreference_emission_to_integral_subgrid(i,j,x_dim_index,i_source).lt.1.or.crossreference_emission_to_integral_subgrid(i,j,x_dim_index,i_source).gt.integral_subgrid_dim(x_dim_index) &
            .or.crossreference_emission_to_integral_subgrid(i,j,y_dim_index,i_source).lt.1.or.crossreference_emission_to_integral_subgrid(i,j,y_dim_index,i_source).gt.integral_subgrid_dim(y_dim_index)) then
            
            write(unit_logfile,'(A,4i,4f)') 'WARNING: crossreference_emission_to_integral_subgrid is out of bounds (i_emis,j_emis,i_integral,j_integral,x_emis,y_emis)',i,j, &
                crossreference_emission_to_integral_subgrid(i,j,x_dim_index,i_source),crossreference_emission_to_integral_subgrid(i,j,y_dim_index,i_source) &
                ,x_emission_subgrid(i,j,i_source)/1000,y_emission_subgrid(i,j,i_source)/1000 &
                ,(x_emission_subgrid(i,j,i_source)-integral_subgrid_min(x_dim_index))/integral_subgrid_delta(x_dim_index)+0.5,(y_emission_subgrid(i,j,i_source)-integral_subgrid_min(y_dim_index))/integral_subgrid_delta(y_dim_index)+0.5
        endif
            
    enddo
    enddo
    do j=1,integral_subgrid_dim(y_dim_index)
    do i=1,integral_subgrid_dim(x_dim_index)
        crossreference_integral_to_emission_subgrid(i,j,x_dim_index,i_source)=1+floor((x_integral_subgrid(i,j)-emission_subgrid_min(x_dim_index,i_source))/emission_subgrid_delta(x_dim_index,i_source))
        crossreference_integral_to_emission_subgrid(i,j,y_dim_index,i_source)=1+floor((y_integral_subgrid(i,j)-emission_subgrid_min(y_dim_index,i_source))/emission_subgrid_delta(y_dim_index,i_source))      
    enddo
    enddo
    
    endif
    enddo
    

    
    end subroutine uEMEP_crossreference_grids
   
    