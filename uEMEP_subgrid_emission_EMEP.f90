!==========================================================================
!   uEMEP model subgrid_emission_EMEP
!==========================================================================
    subroutine uEMEP_subgrid_emission_EMEP
       
    use uEMEP_definitions

    implicit none
    
    integer i,j,k
    character(256) temp_name
    logical exists
    integer ii,jj,tt,t
    integer i_temp,j_temp,i_file
    integer i_nc_temp,j_nc_temp
    real, allocatable :: weighting_nc(:,:),weighting_subgrid(:,:,:)
    real, allocatable :: total_weighting_nc(:,:,:),proxy_weighting_nc(:,:,:)
    real, allocatable :: total_proxy_emission_subgrid(:,:,:,:)
    real, allocatable :: total_proxy_subgrid_emission_in_EMEP_grid(:,:,:,:)
    integer, allocatable :: subgrid_count_nc(:,:)
    integer, allocatable :: subgrid_count_subgrid(:,:,:)
    integer i_nc_start,i_nc_end,j_nc_start,j_nc_end
    integer i_start,i_end,j_start,j_end,t_start,t_end
    real lon_min,lon_max,lat_min,lat_max
    integer i_nc,j_nc
    integer id,jd
    integer subsource_index
    integer i_source,i_subsource
    integer id_p,jd_p,id_m,jd_m,in_p,jn_p,in_m,jn_m
    integer ii_nc,jj_nc,ii_w,jj_w
    integer :: n_weight=3,ii_w0=2,jj_w0=2
    integer weighting_subgrid_dim(2,n_source_index)
    integer i_cross,j_cross
    !integer, allocatable :: crossreference_weighting_to_emep_subgrid(:,:,:,:)
    integer i_w_c,j_w_c
    integer i_nc_c,j_nc_c
    real sum_temp(n_pollutant_loop)
    real xpos_subgrid,ypos_subgrid
    real xpos_subgrid2,ypos_subgrid2
    real EMEP_emission_scaling
    integer i_pollutant
    
    !functions
    
    !Set the scaling factor for EMEP emissions depending on whether they are total emissions or average emissions
    !EMEP_emission_aggregation_period=365.*24.   !For aggregation over a year 
    !EMEP_emission_aggregation_period=1.         !For hourly average        
    
    if (local_subgrid_method_flag.ne.3.and.local_subgrid_method_flag.ne.4) return
    
    write(unit_logfile,'(A)') ''
    write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Distributing EMEP emission to subgrids (uEMEP_subgrid_emission_EMEP)'
	write(unit_logfile,'(A)') '================================================================'

 
    !Allocate and save the existing emission subgrid data
    allocate (total_proxy_emission_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),n_source_index,n_pollutant_loop)) 
    
    !There are no subsources in EMEP so this index is set to 1 and the subsource emissions are transferred to the proxy emission subgrid
    !subsource_index=1
    !do i_source=1,n_source_index
    !if (calculate_source(i_source)) then
    !    temp_proxy_emission_subgrid(:,:,i_source,:)=sum(proxy_emission_subgrid(:,:,i_source,:),3)
    !endif
    !enddo
        
    !Set the start and end times of the loop
    t_start=1
    t_end=subgrid_dim(t_dim_index)
    
    tt=1
    
    do i_source=1,n_source_index
        if (calculate_source(i_source)) then
            do i_pollutant=1,n_pollutant_loop
                !write(*,*) trim(source_file_str(i_source))
                !write(*,*) trim(var_name_nc(conc_nc_index,pollutant_loop_index(i_pollutant),allsource_nc_index)) !trim(pollutant_file_str(pollutant_loop_index(i_pollutant)))
                !write(*,*) sum(emission_subgrid(1:emission_subgrid_dim(x_dim_index,i_source),1:emission_subgrid_dim(y_dim_index,i_source),:,i_source,i_pollutant))
                !write(*,*) (t_end-t_start+1)
                
            !write(unit_logfile,'(A,A,A,A,ES10.2)') 'Emission source ',trim(source_file_str(i_source))//' ',trim(var_name_nc(conc_nc_index,pollutant_loop_index(i_pollutant),allsource_nc_index)),': Total hourly average emissions before use of EMEP (ug/s)=', &
            !    sum(emission_subgrid(1:emission_subgrid_dim(x_dim_index,i_source),1:emission_subgrid_dim(y_dim_index,i_source),:,i_source,i_pollutant))/(t_end-t_start+1)
            enddo
        endif
    enddo
                
    !Distribute the EMEP emissions evenly over the subgrids within an EMEP grid
    if (EMEP_emission_grid_interpolation_flag.eq.0.or.local_subgrid_method_flag.eq.4) then
    !if (EMEP_emission_grid_interpolation_flag.eq.0) then
        write(unit_logfile,'(A)')'Distributing EMEP emissions to all subgrids within an EMEP grid'

        tt=1
        allocate (total_proxy_subgrid_emission_in_EMEP_grid(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),n_source_index,n_pollutant_loop))
        allocate (subgrid_count_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),n_pollutant_loop))
        allocate (subgrid_count_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index)))
        
        !emission_subgrid=0.
        total_proxy_subgrid_emission_in_EMEP_grid=0.
        total_proxy_emission_subgrid=0.
        
        do i_source=1,n_source_index
        if (calculate_source(i_source)) then
                 
            emission_subgrid(:,:,:,i_source,:)=0.
            subgrid_count_subgrid=0
            subgrid_count_nc=0
            
           !Calculate total subgrid emissions and number of subgrids in each EMEP grid
            do j=1,emission_subgrid_dim(y_dim_index,i_source)
            do i=1,emission_subgrid_dim(x_dim_index,i_source)

                ii=crossreference_emission_to_emep_subgrid(i,j,x_dim_index,i_source)
                jj=crossreference_emission_to_emep_subgrid(i,j,y_dim_index,i_source)
                ii=max(min(ii,dim_length_nc(x_dim_nc_index)),1)
                jj=max(min(jj,dim_length_nc(y_dim_nc_index)),1)
                subgrid_count_nc(ii,jj)=subgrid_count_nc(ii,jj)+1
                total_proxy_subgrid_emission_in_EMEP_grid(ii,jj,i_source,:)=total_proxy_subgrid_emission_in_EMEP_grid(ii,jj,i_source,:)+proxy_emission_subgrid(i,j,i_source,:)
                emission_subgrid(i,j,:,i_source,:)=var3d_nc(ii,jj,:,emis_nc_index,i_source,:)
            
            enddo
            enddo
            
            !Transfer the total emissions in the EMEP grid to each subgrid within that grid
            !if (subgrid_emission_distribution_flag) then
            
            do j=1,emission_subgrid_dim(y_dim_index,i_source)
            do i=1,emission_subgrid_dim(x_dim_index,i_source)

                ii=crossreference_emission_to_emep_subgrid(i,j,x_dim_index,i_source)
                jj=crossreference_emission_to_emep_subgrid(i,j,y_dim_index,i_source)
                ii=max(min(ii,dim_length_nc(x_dim_nc_index)),1)
                jj=max(min(jj,dim_length_nc(y_dim_nc_index)),1)
                
                total_proxy_emission_subgrid(i,j,i_source,:)=total_proxy_subgrid_emission_in_EMEP_grid(ii,jj,i_source,:)  
                subgrid_count_subgrid(i,j,:)=subgrid_count_nc(ii,jj)  
                !write(*,*) subgrid_count_nc(ii,jj)
                !Converts from mg/m2/hour(year) to ug/s/subgrid assuming the original EMEP emissions are in mg/m2/hour(year)
                if (hourly_calculations) emission_subgrid(i,j,:,i_source,:)=emission_subgrid(i,j,:,i_source,:) &
                    *emission_subgrid_delta(x_dim_index,i_source)*emission_subgrid_delta(y_dim_index,i_source)*1000./3600.
                if (annual_calculations) emission_subgrid(i,j,:,i_source,:)=emission_subgrid(i,j,:,i_source,:) &
                    *emission_subgrid_delta(x_dim_index,i_source)*emission_subgrid_delta(y_dim_index,i_source)*1000./3600./EMEP_emission_aggregation_period
                    

                !write(*,*) emission_subgrid(i,j,1,i_source,subsource_index),var3d_nc(ii,jj,1,emis_nc_index,i_source)*emission_subgrid_delta(x_dim_index,i_source)*emission_subgrid_delta(y_dim_index,i_source)*1000./3600.
            enddo
            enddo
                 
            !endif
            
        
           
            !Determine subgrid normalised time profile per hour from EMEP grid emissions (average hourly emission conversion)
            !This is not quite right because the entire emission time profile is not available in a short period
            !Minimum of a day is needed, with the assumption that all days are the same
            if (local_subgrid_method_flag.eq.4) then
                write(unit_logfile,'(A)')'Calculating EMEP emission time profile'
                do j=1,emission_subgrid_dim(y_dim_index,i_source)
                do i=1,emission_subgrid_dim(x_dim_index,i_source)
                    if (hourly_calculations) then
                        sum_temp(:)=sum(emission_subgrid(i,j,:,i_source,:),1)
                        do i_pollutant=1,n_pollutant_loop
                            emission_time_profile_subgrid(i,j,:,i_source,i_pollutant)=emission_subgrid(i,j,:,i_source,i_pollutant)/sum_temp(i_pollutant)*emission_subgrid_dim(t_dim_index,i_source)
                            if (sum_temp(i_pollutant).eq.0.) emission_time_profile_subgrid(i,j,:,i_source,i_pollutant)=0.
                        enddo
                        
                    else
                        emission_time_profile_subgrid(i,j,:,i_source,:)=1.
                    endif
                    !write(*,'(<emission_subgrid_dim(t_dim_index,i_source)>f6.2)') emission_time_profile_subgrid(i,j,:,i_source,subsource_index)
                    !write(*,*) i,j,subgrid_count_subgrid(i,j,1,i_source),total_proxy_emission_subgrid(i,j,1,i_source,subsource_index),emission_subgrid(i,j,1,i_source,subsource_index)
                    
                    !Set emissions to 0 in the case when local_subgrid_method_flag.eq.4 since these are set later
                    !This way of doing things is not logical as it fills the grid unnecessarilly. Should be fixed and made logical
                    emission_subgrid(i,j,:,i_source,:)=0.
                enddo
                enddo
            endif
            
            !Distribute EMEP emissions to existing proxy subgrid emissions
            if (subgrid_emission_distribution_flag.and.local_subgrid_method_flag.ne.4) then
                if (local_subgrid_method_flag.eq.2) write(unit_logfile,'(2A)')'Distributing local emission data to subgrid emissions for: ',trim(source_file_str(i_source))
                if (local_subgrid_method_flag.eq.3) write(unit_logfile,'(2A)')'Distributing EMEP emissions to proxy subgrid emissions, no weighting used, for: ',trim(source_file_str(i_source))
                
                do t=t_start,t_end
                    emission_subgrid(:,:,t,i_source,:)=emission_subgrid(:,:,t,i_source,:)*subgrid_count_subgrid(:,:,:)*proxy_emission_subgrid(:,:,i_source,:)/total_proxy_emission_subgrid(:,:,i_source,:)              
                    where (total_proxy_emission_subgrid(:,:,i_source,:).eq.0.) emission_subgrid(:,:,t,i_source,:)=0.
                enddo
                
                !Fix round off negatives
                where (emission_subgrid(:,:,:,i_source,:).lt.0) emission_subgrid(:,:,:,i_source,:)=0.
                
            !write(*,*) sum(emission_subgrid(:,:,:,traffic_index,pollutant_loop_back_index(pmex_nc_index))),sum(emission_subgrid(:,:,:,traffic_index,pollutant_loop_back_index(pm25_nc_index)))
            !write(*,*) sum(var3d_nc(:,:,:,emis_nc_index,traffic_index,pollutant_loop_back_index(pmex_nc_index))),sum(var3d_nc(:,:,:,emis_nc_index,traffic_index,pollutant_loop_back_index(pm25_nc_index)))
            !write(*,*) pmex_index,pmex_nc_index
            !stop
            endif    

        endif
        enddo


    endif

    
    !Quick calculation of area weighting, no edge effects. Does not need to change with time
    !This is done also if there is moving window weighting later as it is used for the nonlocal contribution
    if (EMEP_emission_grid_interpolation_flag.eq.1.and.local_subgrid_method_flag.ne.4) then
        
        write(unit_logfile,'(A)')'Distributing EMEP emissions to proxy emission subgrids using area weighting of the EMEP grid'
    
        tt=1
        !tt=dim_length_nc(time_dim_nc_index)
        
        allocate (weighting_nc(n_weight,n_weight)) !EMEP grid weighting for interpolation. Does not need a source index for area weighting
        !allocate (subgrid_count_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),tt,n_source_index))
        allocate (subgrid_count_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),n_pollutant_loop))
        
        
        !allocate (emission_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_source_index,n_pollutant_loop)) 
        total_proxy_emission_subgrid=0.

        
        do i_source=1,n_source_index
        if (calculate_source(i_source)) then
            
            emission_subgrid(:,:,:,i_source,:)=0.
            subgrid_count_subgrid=0

            do j=1,emission_subgrid_dim(y_dim_index,i_source)
            do i=1,emission_subgrid_dim(x_dim_index,i_source)
        
                !Only calculate for valid emission subgrids when using the proxy emissions for distribution
                if (.not.subgrid_emission_distribution_flag.or.sum(proxy_emission_subgrid(i,j,i_source,:)).ne.0) then
                    
                !Assumes it is never on the edge of the EMEP grid
                i_nc=crossreference_emission_to_emep_subgrid(i,j,x_dim_index,i_source)
                j_nc=crossreference_emission_to_emep_subgrid(i,j,y_dim_index,i_source)
        
                weighting_nc=0.
            
                if (EMEP_projection_type.eq.LL_projection_index) then
                    xpos_subgrid=lon_emission_subgrid(i,j,i_source)
                    ypos_subgrid=lat_emission_subgrid(i,j,i_source)
                elseif (EMEP_projection_type.eq.LCC_projection_index) then
                    if (use_alternative_LCC_projection_flag) then
                        call lb2lambert2_uEMEP(xpos_subgrid,ypos_subgrid,lon_emission_subgrid(i,j,i_source),lat_emission_subgrid(i,j,i_source),EMEP_projection_attributes)
                    else
                        call lb2lambert_uEMEP(xpos_subgrid,ypos_subgrid,lon_emission_subgrid(i,j,i_source),lat_emission_subgrid(i,j,i_source),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                    endif
                    !call lb2lambert_uEMEP(xpos_subgrid,ypos_subgrid,lon_emission_subgrid(i,j,i_source),lat_emission_subgrid(i,j,i_source),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                elseif (EMEP_projection_type.eq.PS_projection_index) then
                    call LL2PS_spherical(xpos_subgrid,ypos_subgrid,lon_emission_subgrid(i,j,i_source),lat_emission_subgrid(i,j,i_source),EMEP_projection_attributes)
                endif

                !Calculate the area weighted EMEP grid emissions at each subgrid
                do jj=-1,+1
                do ii=-1,+1
                
                    ii_nc=ii+i_nc
                    jj_nc=jj+j_nc
                
                    ii_w=ii+2
                    jj_w=jj+2
                

                    lon_min=max(xpos_subgrid-dgrid_nc(lon_nc_index)/2.,var1d_nc(ii_nc,lon_nc_index)-dgrid_nc(lon_nc_index)/2.)
                    lon_max=min(xpos_subgrid+dgrid_nc(lon_nc_index)/2.,var1d_nc(ii_nc,lon_nc_index)+dgrid_nc(lon_nc_index)/2.)
                    lat_min=max(ypos_subgrid-dgrid_nc(lat_nc_index)/2.,var1d_nc(jj_nc,lat_nc_index)-dgrid_nc(lat_nc_index)/2.)
                    lat_max=min(ypos_subgrid+dgrid_nc(lat_nc_index)/2.,var1d_nc(jj_nc,lat_nc_index)+dgrid_nc(lat_nc_index)/2.)
            
                    if (lon_max.gt.lon_min.and.lat_max.gt.lat_min) then
                        weighting_nc(ii_w,jj_w)=(lat_max-lat_min)*(lon_max-lon_min)/dgrid_nc(lon_nc_index)/dgrid_nc(lat_nc_index)
                    else
                        weighting_nc(ii_w,jj_w)=0.
                    endif                

                    emission_subgrid(i,j,:,i_source,:)=emission_subgrid(i,j,:,i_source,:)+var3d_nc(ii_nc,jj_nc,:,emis_nc_index,i_source,:)*weighting_nc(ii_w,jj_w)
                    
                enddo
                enddo
                
                !write(*,*) sum(weighting_nc) !This is OK
                 
                !Calculate the total subgrid emissions in the EMEP grid region surrounding the subgrid
                i_start=max(1,i-emission_subgrid_loop_index(x_dim_index,i_source))
                i_end=min(emission_subgrid_dim(x_dim_index,i_source),i+emission_subgrid_loop_index(x_dim_index,i_source))
                j_start=max(1,j-emission_subgrid_loop_index(y_dim_index,i_source))
                j_end=min(emission_subgrid_dim(y_dim_index,i_source),j+emission_subgrid_loop_index(y_dim_index,i_source))
                do jj=j_start,j_end
                do ii=i_start,i_end

                    if (EMEP_projection_type.eq.LL_projection_index) then
                        xpos_subgrid2=lon_emission_subgrid(ii,jj,i_source)
                        ypos_subgrid2=lat_emission_subgrid(ii,jj,i_source)
                    elseif (EMEP_projection_type.eq.LCC_projection_index) then
                        if (use_alternative_LCC_projection_flag) then
                            call lb2lambert2_uEMEP(xpos_subgrid2,ypos_subgrid2,lon_emission_subgrid(ii,jj,i_source),lat_emission_subgrid(ii,jj,i_source),EMEP_projection_attributes)
                        else
                            call lb2lambert_uEMEP(xpos_subgrid2,ypos_subgrid2,lon_emission_subgrid(ii,jj,i_source),lat_emission_subgrid(ii,jj,i_source),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                        endif
                        !call lb2lambert_uEMEP(xpos_subgrid2,ypos_subgrid2,lon_emission_subgrid(ii,jj,i_source),lat_emission_subgrid(ii,jj,i_source),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                    elseif (EMEP_projection_type.eq.PS_projection_index) then
                        call LL2PS_spherical(xpos_subgrid2,ypos_subgrid2,lon_emission_subgrid(ii,jj,i_source),lat_emission_subgrid(ii,jj,i_source),EMEP_projection_attributes)
                    endif

                    if (abs(xpos_subgrid-xpos_subgrid2).le.dgrid_nc(lon_nc_index)/2. &
                        .and.abs(ypos_subgrid-ypos_subgrid2).le.dgrid_nc(lat_nc_index)/2.) then
                            
                        total_proxy_emission_subgrid(i,j,i_source,:)=total_proxy_emission_subgrid(i,j,i_source,:)+proxy_emission_subgrid(ii,jj,i_source,:)
                        subgrid_count_subgrid(i,j,:)=subgrid_count_subgrid(i,j,:)+1

                    endif                   

                enddo
                enddo
                                
                !Converts from mg/subgrid to ug/s/subgrid assuming the original EMEP emissions are in mg/m2/hour
                if (hourly_calculations) emission_subgrid(i,j,:,i_source,:)=emission_subgrid(i,j,:,i_source,:)*emission_subgrid_delta(x_dim_index,i_source)*emission_subgrid_delta(y_dim_index,i_source)*1000./3600.
                if (annual_calculations) emission_subgrid(i,j,:,i_source,:)=emission_subgrid(i,j,:,i_source,:)*emission_subgrid_delta(x_dim_index,i_source)*emission_subgrid_delta(y_dim_index,i_source)*1000./3600./EMEP_emission_aggregation_period

            endif !If valid emission point
                
            enddo
            if (mod(j,100).eq.0) write(*,'(A,A,A,i4,A,i4)') 'Subgrid emission EMEP interpolation for ',trim(source_file_str(i_source)),': ',j,' of ',emission_subgrid_dim(2,i_source)
            enddo
        
       
            !Distribute EMEP emissions to subgrid emissions
            if (subgrid_emission_distribution_flag) then
                write(unit_logfile,'(A)')'Distributing EMEP emissions to subgrid emissions within an area weighted EMEP grid'
                do t=t_start,t_end
                    emission_subgrid(:,:,t,i_source,:)=emission_subgrid(:,:,t,i_source,:)*subgrid_count_subgrid(:,:,:)*proxy_emission_subgrid(:,:,i_source,:)/total_proxy_emission_subgrid(:,:,i_source,:)              
                    where (total_proxy_emission_subgrid(:,:,i_source,:).eq.0.) emission_subgrid(:,:,t,i_source,:)=0.
                enddo
            endif    

        endif
        enddo !source loop
        
    endif

    !Loop through subgrid and carry out a subgrid weighted moving window interpolation
    if (EMEP_emission_grid_interpolation_flag.eq.2.and.local_subgrid_method_flag.ne.4) then

        write(unit_logfile,'(A)')'Distributing EMEP emissions to proxy subgrids using emission weighting of the EMEP grid'
        
        !NOTE: currently does the inteprolation for all time steps which is not ncessary. Only needs to do it for 1. Same with the area interpolation.
        !Only use time tt for weighting distribution to increase speed. It is possible this can change with time. Replace 'tt' with ':'
        tt=1
        !tt=dim_length_nc(time_dim_nc_index)
        
        allocate (total_weighting_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),n_pollutant_loop)) !EMEP grid weighting for interpolation
        allocate (proxy_weighting_nc(5,5,n_pollutant_loop)) !EMEP grid weighting for interpolation
        !allocate (subgrid_count_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),tt,n_source_index))
        allocate (subgrid_count_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),n_pollutant_loop))

        i_w_c=3
        j_w_c=3
        
        
        total_proxy_emission_subgrid=0.
        !emission_subgrid=0.
        
        !Set the start and end times of the loop
        !t_start=1
        !t_end=subgrid_dim(t_dim_index)
               
        allocate (weighting_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),n_pollutant_loop))
        weighting_subgrid_dim(:,:)=emission_subgrid_dim(1:2,:)
        
        !Calculate weighting sum for each EMEP grid.
        total_weighting_nc=0.
        do i_source=1,n_source_index
        if (calculate_source(i_source)) then
            
            emission_subgrid(:,:,:,i_source,:)=0.
            subgrid_count_subgrid=0
            
            weighting_subgrid(:,:,:)=proxy_emission_subgrid(:,:,i_source,:)
            
            !Calculate the total weighting (emission) in each emep grid
            do j=1,weighting_subgrid_dim(y_dim_index,i_source)
            do i=1,weighting_subgrid_dim(x_dim_index,i_source)

                i_nc=crossreference_emission_to_emep_subgrid(i,j,x_dim_index,i_source)
                j_nc=crossreference_emission_to_emep_subgrid(i,j,y_dim_index,i_source)
                total_weighting_nc(i_nc,j_nc,:)=total_weighting_nc(i_nc,j_nc,:)+weighting_subgrid(i,j,:)
                !write(*,*) i_source,i,j,i_nc,j_nc,weighting_subgrid(i,j,:,i_source)
                !write(*,*) i,j,i_nc,j_nc,i_nc-crossreference_weighting_to_emep_subgrid(i,j,x_dim_index,i_source),j_nc-crossreference_weighting_to_emep_subgrid(i,j,y_dim_index,i_source)
        
            enddo
            enddo
        
            !Calculate the proxy weighting in the nearest emep grids for each subgrid
        
            do j=1,emission_subgrid_dim(y_dim_index,i_source)
            do i=1,emission_subgrid_dim(x_dim_index,i_source)
            
                !Only calculate for valid emission subgrids when using the proxy for distribution
                if (.not.subgrid_emission_distribution_flag.or.sum(proxy_emission_subgrid(i,j,i_source,:)).ne.0) then

                proxy_weighting_nc=0.
                       
                i_cross=i !crossreference_target_to_emission_subgrid(i,j,x_dim_index,i_source)
                j_cross=j !crossreference_target_to_emission_subgrid(i,j,y_dim_index,i_source)                   
                i_nc_c=crossreference_emission_to_emep_subgrid(i_cross,j_cross,x_dim_index,i_source)
                j_nc_c=crossreference_emission_to_emep_subgrid(i_cross,j_cross,y_dim_index,i_source)                   
                !Limit the loop so that it doesn't go over more than necessary subgrids and does not go outside the domain
                i_start=max(1,i-emission_subgrid_loop_index(x_dim_index,i_source))
                i_end=min(emission_subgrid_dim(x_dim_index,i_source),i+emission_subgrid_loop_index(x_dim_index,i_source))
                j_start=max(1,j-emission_subgrid_loop_index(y_dim_index,i_source))
                j_end=min(emission_subgrid_dim(y_dim_index,i_source),j+emission_subgrid_loop_index(y_dim_index,i_source))
                !write(*,*) i,j,i_end-i_start,j_end-j_start
                
                if (EMEP_projection_type.eq.LL_projection_index) then
                    xpos_subgrid=lon_emission_subgrid(i,j,i_source)
                    ypos_subgrid=lat_emission_subgrid(i,j,i_source)
                elseif (EMEP_projection_type.eq.LCC_projection_index) then
                    if (use_alternative_LCC_projection_flag) then
                        call lb2lambert2_uEMEP(xpos_subgrid,ypos_subgrid,lon_emission_subgrid(i,j,i_source),lat_emission_subgrid(i,j,i_source),EMEP_projection_attributes)
                    else
                        call lb2lambert_uEMEP(xpos_subgrid,ypos_subgrid,lon_emission_subgrid(i,j,i_source),lat_emission_subgrid(i,j,i_source),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                    endif
                    !call lb2lambert_uEMEP(xpos_subgrid,ypos_subgrid,lon_emission_subgrid(i,j,i_source),lat_emission_subgrid(i,j,i_source),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                elseif (EMEP_projection_type.eq.PS_projection_index) then
                    call LL2PS_spherical(xpos_subgrid,ypos_subgrid,lon_emission_subgrid(i,j,i_source),lat_emission_subgrid(i,j,i_source),EMEP_projection_attributes)
                endif

                do jj=j_start,j_end
                do ii=i_start,i_end                
                    
                    if (EMEP_projection_type.eq.LL_projection_index) then
                        xpos_subgrid2=lon_emission_subgrid(ii,jj,i_source)
                        ypos_subgrid2=lat_emission_subgrid(ii,jj,i_source)
                    elseif (EMEP_projection_type.eq.LCC_projection_index) then
                        if (use_alternative_LCC_projection_flag) then
                            call lb2lambert2_uEMEP(xpos_subgrid2,ypos_subgrid2,lon_emission_subgrid(ii,jj,i_source),lat_emission_subgrid(ii,jj,i_source),EMEP_projection_attributes)
                        else
                            call lb2lambert_uEMEP(xpos_subgrid2,ypos_subgrid2,lon_emission_subgrid(ii,jj,i_source),lat_emission_subgrid(ii,jj,i_source),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                        endif
                        !call lb2lambert_uEMEP(xpos_subgrid2,ypos_subgrid2,lon_emission_subgrid(ii,jj,i_source),lat_emission_subgrid(ii,jj,i_source),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                    elseif (EMEP_projection_type.eq.PS_projection_index) then
                        call LL2PS_spherical(xpos_subgrid2,ypos_subgrid2,lon_emission_subgrid(ii,jj,i_source),lat_emission_subgrid(ii,jj,i_source),EMEP_projection_attributes)
                    endif

                    if (abs(xpos_subgrid-xpos_subgrid2).le.dgrid_nc(lon_nc_index)/2. &
                        .and.abs(ypos_subgrid-ypos_subgrid2).le.dgrid_nc(lat_nc_index)/2.) then                  
                        i_nc=crossreference_emission_to_emep_subgrid(ii,jj,x_dim_index,i_source)
                        j_nc=crossreference_emission_to_emep_subgrid(ii,jj,y_dim_index,i_source)                   
                        !proxy_weighting_nc(i_nc,j_nc,:,i_source)=proxy_weighting_nc(i_nc,j_nc,:,i_source)+weighting_subgrid(ii,jj,:,i_source)
                        proxy_weighting_nc(i_nc-i_nc_c+i_w_c,j_nc-j_nc_c+j_w_c,:)=proxy_weighting_nc(i_nc-i_nc_c+i_w_c,j_nc-j_nc_c+j_w_c,:)+weighting_subgrid(ii,jj,:)
                        
                        total_proxy_emission_subgrid(i,j,i_source,:)=total_proxy_emission_subgrid(i,j,i_source,:)+weighting_subgrid(ii,jj,:) !Weighting subgrid is the same as the existing proxy subgrid but without the subsource
                        subgrid_count_subgrid(i,j,:)=subgrid_count_subgrid(i,j,:)+1

                    endif                   
                enddo
                enddo
                
                i_cross=i !crossreference_target_to_emission_subgrid(i,j,x_dim_index,i_source)
                j_cross=j !crossreference_target_to_emission_subgrid(i,j,y_dim_index,i_source)    
                i_nc=crossreference_emission_to_emep_subgrid(i_cross,j_cross,x_dim_index,i_source)
                j_nc=crossreference_emission_to_emep_subgrid(i_cross,j_cross,y_dim_index,i_source)
                i_nc_start=max(1,i_nc-1)
                i_nc_end=min(dim_length_nc(x_dim_nc_index),i_nc+1)
                j_nc_start=max(1,j_nc-1)
                j_nc_end=min(dim_length_nc(y_dim_nc_index),j_nc+1)

                !write(*,*) i,j,i_nc,j_nc,i_nc_end-i_nc_start,j_nc_end-j_nc_start

                do jj=j_nc_start,j_nc_end
                do ii=i_nc_start,i_nc_end
            
                    proxy_weighting_nc(ii-i_nc+i_w_c,jj-j_nc+j_w_c,:)=proxy_weighting_nc(ii-i_nc+i_w_c,jj-j_nc+j_w_c,:)/total_weighting_nc(ii,jj,:)
                    !where (total_weighting_nc(ii,jj,tt,i_source).eq.0) proxy_weighting_nc(ii-i_nc+i_w_c,jj-j_nc+j_w_c,:,i_source)=0.
                    where (total_weighting_nc(ii,jj,:).eq.0) proxy_weighting_nc(ii-i_nc+i_w_c,jj-j_nc+j_w_c,:)=0.
                    
                enddo
                enddo
            
                                   
         
            !Add up the contributing weights
            !Note that only the subsource_index=1 can be determined since EMEP has no subsources
            !if (hourly_calculations) then
                !do i_subsource=1,n_subsource_index(i_source)
                do jj=j_nc_start,j_nc_end
                do ii=i_nc_start,i_nc_end
                    do i_pollutant=1,n_pollutant_loop
                    emission_subgrid(i,j,:,i_source,i_pollutant)=emission_subgrid(i,j,:,i_source,i_pollutant) &
                        +var3d_nc(ii,jj,:,emis_nc_index,i_source,i_pollutant)*proxy_weighting_nc(ii-i_nc+i_w_c,jj-j_nc+j_w_c,i_pollutant)
                    enddo
                    !if (weighting_nc(ii,jj,traffic_index).ne.0) then
                    !    write(*,*) i,j,ii,jj,weighting_nc(ii,jj,traffic_index),subgrid(i,j,traffic_index,emep_local_subgrid_index)
                    !endif
                    
                enddo
                enddo
                !write(*,*) i,j,proxy_weighting_nc(2,:,1,i_source)
                !write(*,*) i,j,proxy_weighting_nc(2,2,1,i_source),emission_subgrid(i,j,1,i_source,subsource_index)
                
                !write(*,*) sum(proxy_weighting_nc)
                
            !Converts from mg/subgrid to ug/s/subgrid assuming the original EMEP emissions are in mg/m2/hour
            if (hourly_calculations) emission_subgrid(i,j,:,i_source,:)=emission_subgrid(i,j,:,i_source,:)*emission_subgrid_delta(x_dim_index,i_source)*emission_subgrid_delta(y_dim_index,i_source)*1000./3600.
            if (annual_calculations) emission_subgrid(i,j,:,i_source,:)=emission_subgrid(i,j,:,i_source,:)*emission_subgrid_delta(x_dim_index,i_source)*emission_subgrid_delta(y_dim_index,i_source)*1000./3600./EMEP_emission_aggregation_period

            endif !If valid emission subgrid
                
        enddo  
        if (mod(j,100).eq.0) write(*,'(A,A,A,i4,A,i4)') 'Subgrid emission EMEP interpolation for ',trim(source_file_str(i_source)),': ',j,' of ',emission_subgrid_dim(2,i_source)
        enddo
        
            !Puts counts into the time array. If emission distributions change in time this needs to be replaced
            !do t=1,size(subgrid_count_subgrid,3)
            !    subgrid_count_subgrid(:,:,t,i_source)=subgrid_count_subgrid(:,:,tt,i_source)
            !enddo
            

        !Distribute EMEP emissions to existing subgrid emissions
        if (subgrid_emission_distribution_flag) then
            write(unit_logfile,'(A)')'Distributing EMEP emissions to subgrid emissions within an emission weighted EMEP grid'
            do t=t_start,t_end
                emission_subgrid(:,:,t,i_source,:)=emission_subgrid(:,:,t,i_source,:)*subgrid_count_subgrid(:,:,:)*proxy_emission_subgrid(:,:,i_source,:)/total_proxy_emission_subgrid(:,:,i_source,:)              
                where (total_proxy_emission_subgrid(:,:,i_source,:).eq.0.) emission_subgrid(:,:,t,i_source,:)=0.
            enddo
        endif    

        endif !End if calculate_source
        enddo !End source loop

    endif
    
    !if (subgrid_emission_distribution_flag) then
    do i_source=1,n_source_index
        if (calculate_source(i_source)) then
            do i_pollutant=1,n_pollutant_loop
            write(unit_logfile,'(A,A,A,ES10.2)') 'Emission source ',trim(source_file_str(i_source))//' '//trim(pollutant_file_str(pollutant_loop_index(i_pollutant))),': Total hourly average emissions after use of EMEP (ug/s)=', &
                sum(emission_subgrid(1:emission_subgrid_dim(x_dim_index,i_source),1:emission_subgrid_dim(y_dim_index,i_source),:,i_source,i_pollutant))/(t_end-t_start+1)
            enddo
        endif
    enddo
    !endif    


    if (allocated(weighting_nc)) deallocate(weighting_nc)
    if (allocated(total_weighting_nc)) deallocate(total_weighting_nc)
    if (allocated(proxy_weighting_nc)) deallocate(proxy_weighting_nc)
    if (allocated(weighting_subgrid)) deallocate(weighting_subgrid)
    if (allocated(total_proxy_emission_subgrid)) deallocate(total_proxy_emission_subgrid)
    if (allocated(total_proxy_subgrid_emission_in_EMEP_grid)) deallocate(total_proxy_subgrid_emission_in_EMEP_grid)
    if (allocated(subgrid_count_nc)) deallocate(subgrid_count_nc)
    if (allocated(subgrid_count_subgrid)) deallocate(subgrid_count_subgrid)
    
    end subroutine uEMEP_subgrid_emission_EMEP
