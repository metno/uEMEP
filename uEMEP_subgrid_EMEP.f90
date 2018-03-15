!==========================================================================
!   uEMEP model subgrid_EMEP
!==========================================================================
    subroutine uEMEP_subgrid_EMEP
       
    use uEMEP_definitions

    implicit none
    
    character(256) temp_name
    logical exists
    integer ii,jj,tt
    integer i_temp,j_temp,i_file
    integer i_nc_temp,j_nc_temp
    real, allocatable :: weighting_nc(:,:,:,:),weighting_subgrid(:,:,:,:)
    real, allocatable :: total_weighting_nc(:,:,:,:),proxy_weighting_nc(:,:,:,:)
    real, allocatable :: area_weighting_nc(:,:,:,:,:,:)
    integer i_nc_start,i_nc_end,j_nc_start,j_nc_end
    integer i_start,i_end,j_start,j_end,t_start,t_end
    real lon_min,lon_max,lat_min,lat_max
    real lon_area_min,lon_area_max,lat_area_min,lat_area_max
    integer i_nc,j_nc
    integer id,jd
    integer source_index,emep_subsource
    real, allocatable :: nonlocal_correction(:,:)
    real, allocatable :: nonlocal_correction_average(:)
    integer i_source,i_subsource
    integer id_p,jd_p,id_m,jd_m,in_p,jn_p,in_m,jn_m
    integer ii_nc,jj_nc,ii_w,jj_w
    integer :: n_weight=3,ii_w0=2,jj_w0=2
    integer weighting_subgrid_dim(2,n_source_index)
    integer i_cross,j_cross
    integer i_cross_integral,j_cross_integral
    integer, allocatable :: crossreference_weighting_to_emep_subgrid(:,:,:,:)
    integer i_w_c,j_w_c
    integer i_nc_c,j_nc_c
    integer i_comp
    integer count
    integer ii_start,ii_end,jj_start,jj_end
    real u_subgrid_loc,v_subgrid_loc
    real x_pos,y_pos
    real lon_limit,lat_limit
    real xpos_subgrid,ypos_subgrid
    real xpos_emission_subgrid,ypos_emission_subgrid
    real xpos_integral_subgrid,ypos_integral_subgrid
    real EMEP_grid_interpolation_size_sqr
    integer :: tt_dim=1

    
    !functions
    write(unit_logfile,'(A)') ''
    write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Distributing EMEP concentrations to subgrids  (uEMEP_subgrid_EMEP)'
	write(unit_logfile,'(A)') '================================================================'

    subgrid(:,:,:,emep_subgrid_index,:,:)=0
    subgrid(:,:,:,emep_frac_subgrid_index,:,:)=0
    subgrid(:,:,:,emep_local_subgrid_index,:,:)=0
    subgrid(:,:,:,emep_nonlocal_subgrid_index,:,:)=0
    comp_EMEP_subgrid(:,:,:,:)=0
    orig_EMEP_subgrid(:,:,:,:)=0
    
    EMEP_grid_interpolation_size_sqr=EMEP_grid_interpolation_size*EMEP_grid_interpolation_size
    
    tt_dim=1
    
    !CHECK
    !Need to check the time element of the nonlocal_correction when there is no time loop. Something strange
    if (.not.allocated(nonlocal_correction)) allocate (nonlocal_correction(tt_dim,n_source_index))
    if (.not.allocated(nonlocal_correction_average)) allocate (nonlocal_correction_average(n_source_index))

    !There are no subsources in EMEP
    emep_subsource=1

    !If EMEP gridded to subgrid files already exist then read them in and leave the subroutine
    do source_index=1,n_source_index
        if (read_existing_grid_data(emep_subgrid_file_index(source_index)).and.calculate_source(source_index)) then
           
            i_file=emep_subgrid_file_index(source_index)
            temp_name=trim(pathname_grid(i_file))//trim(filename_grid(i_file))//trim(subsource_str(emep_subsource))//'_'//trim(file_tag)//'.asc'
            inquire(file=temp_name,exist=exists)
            if (.not.exists) then
                write(unit_logfile,*)'WARNING: '//trim(temp_name)//' does not exist.'
                !return
            else
                call read_esri_ascii_3d_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),subgrid_delta(1),subgrid(:,:,:,emep_subgrid_index,source_index,emep_subsource),x_subgrid,y_subgrid)
            endif
            i_file=emep_subgrid_frac_file_index(source_index)
            temp_name=trim(pathname_grid(i_file))//trim(filename_grid(i_file))//trim(subsource_str(emep_subsource))//'_'//trim(file_tag)//'.asc'
            inquire(file=temp_name,exist=exists)
            if (.not.exists) then
                write(unit_logfile,*)'WARNING: '//trim(temp_name)//' does not exist.'
                !return
            else
                call read_esri_ascii_3d_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),subgrid_delta(1),subgrid(:,:,:,emep_frac_subgrid_index,source_index,emep_subsource),x_subgrid,y_subgrid)
            endif
            subgrid(:,:,:,emep_local_subgrid_index,source_index,emep_subsource)=subgrid(:,:,:,emep_subgrid_index,source_index,emep_subsource)*subgrid(:,:,:,emep_frac_subgrid_index,source_index,emep_subsource)
            subgrid(:,:,:,emep_nonlocal_subgrid_index,source_index,emep_subsource)=subgrid(:,:,:,emep_subgrid_index,source_index,emep_subsource)*(1.-subgrid(:,:,:,emep_frac_subgrid_index,source_index,emep_subsource))
            
        endif
    enddo
 
    if (read_existing_grid_data(emep_subgrid_file_index(allsource_index))) then
        write(unit_logfile,'(A,A)') 'Reading existing local and nonlocal files: ',trim(temp_name)
        return
    endif
    
    !tt=1
    
    nonlocal_correction_average=0.
    
 
    !Save the original EMEP directly to subgrid for comparison purposes
    
    write(unit_logfile,'(A)')'Calculating original EMEP subgrid using nearest neighbour interpolation'
        
        do j=1,subgrid_dim(y_dim_index)
        do i=1,subgrid_dim(x_dim_index)
            
            ii=crossreference_target_to_emep_subgrid(i,j,x_dim_index)
            jj=crossreference_target_to_emep_subgrid(i,j,y_dim_index)
         
            !Interpolate the other EMEP compounds as well to subgrid
            do i_comp=1,n_compound_loop
                orig_EMEP_subgrid(i,j,:,compound_loop_index(i_comp))=comp_var3d_nc(ii,jj,:,compound_loop_index(i_comp))
            enddo

        enddo
        enddo

    !Loop through subgrid and find those subgrids within EMEP grids and allocate concentrations directly from EMEP grids.

    if (EMEP_grid_interpolation_flag.eq.0) then
        
        if (tt.eq.t_start) write(unit_logfile,'(A)')'Calculating EMEP local subgrid contribution using nearest neighbour interpolation'
        
        do j=1,subgrid_dim(y_dim_index)
        do i=1,subgrid_dim(x_dim_index)
        if (use_subgrid(i,j,allsource_index)) then
            
            ii=crossreference_target_to_emep_subgrid(i,j,x_dim_index)
            jj=crossreference_target_to_emep_subgrid(i,j,y_dim_index)
        
            !if (hourly_calculations) then
            !    subgrid(i,j,:,emep_subgrid_index,:,subsource_index)=var4d_nc(ii,jj,surface_level_nc,:,conc_nc_index,:)
            !    subgrid(i,j,:,emep_frac_subgrid_index,:,subsource_index)=var4d_nc(ii,jj,surface_level_nc,:,frac_nc_index,:)
            !else
                subgrid(i,j,:,emep_subgrid_index,:,emep_subsource)=var3d_nc(ii,jj,:,conc_nc_index,:)
                subgrid(i,j,:,emep_local_subgrid_index,:,emep_subsource)=var3d_nc(ii,jj,:,local_nc_index,:)
            !endif
                    
                !Interpolate the other EMEP compounds as well to subgrid
                do i_comp=1,n_compound_loop
                    comp_EMEP_subgrid(i,j,:,compound_loop_index(i_comp))=comp_var3d_nc(ii,jj,:,compound_loop_index(i_comp))
                enddo

        endif
        enddo
        enddo
        
        subgrid(:,:,:,emep_nonlocal_subgrid_index,:,emep_subsource)=subgrid(:,:,:,emep_subgrid_index,:,emep_subsource)-subgrid(:,:,:,emep_local_subgrid_index,:,emep_subsource)
   endif

    
    !Set the start and end times of the loop
    t_start=1
    t_end=subgrid_dim(t_dim_index)

    !Loop through the time and the subgrids
    do tt=t_start,t_end

    !Quick calculation of area weighting, no edge effects. Does not need to change with time
    !This is done also if there is moving window weighting later as it is used for the nonlocal contribution
    if (EMEP_grid_interpolation_flag.ge.1) then
        
        if (tt.eq.t_start) write(unit_logfile,'(A)')'Calculating EMEP local subgrid contribution using area weighted interpolation'
        
        if (use_downwind_position_flag.and.hourly_calculations) then
            n_weight=1+2*floor(0.5*(EMEP_grid_interpolation_size-1.))
            n_weight=3+2*floor(EMEP_grid_interpolation_size*0.5)
            ii_w0=1+floor(n_weight*.5)
            jj_w0=1+floor(n_weight*.5)
        else
            n_weight=3+2*floor(EMEP_grid_interpolation_size*0.5)
            ii_w0=1+floor(n_weight*.5)
            jj_w0=1+floor(n_weight*.5)
        endif
        
        
        if (tt.eq.t_start) write(unit_logfile,*) 'Weighting grid dimensions and centres: ',n_weight,ii_w0,jj_w0
    
        if (.not.allocated(weighting_nc)) allocate (weighting_nc(n_weight,n_weight,tt_dim,n_source_index)) !EMEP grid weighting for interpolation. Does not need a source index for area weighting
        if (.not.allocated(area_weighting_nc)) allocate (area_weighting_nc(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),n_weight,n_weight,tt_dim,n_source_index)) !EMEP grid weighting for area interpolation

        subgrid(:,:,tt,emep_subgrid_index,:,:)=0
        subgrid(:,:,tt,emep_frac_subgrid_index,:,:)=0
        subgrid(:,:,tt,emep_local_subgrid_index,:,:)=0
        subgrid(:,:,tt,emep_nonlocal_subgrid_index,:,:)=0
        comp_EMEP_subgrid(:,:,tt,:)=0
        
        !Cover the search area necessary for the surorunding EMEP grids
        jj_start=-1-floor(0.5*(EMEP_grid_interpolation_size-1.))
        ii_start=-1-floor(0.5*(EMEP_grid_interpolation_size-1.))
        jj_end=1+floor(0.5*(EMEP_grid_interpolation_size-1.))
        ii_end=1+floor(0.5*(EMEP_grid_interpolation_size-1.))
        
        !Set the size of the region surorunding the target grid that is searched
        lon_limit=dgrid_nc(lon_nc_index)/2.*EMEP_grid_interpolation_size
        lat_limit=dgrid_nc(lat_nc_index)/2.*EMEP_grid_interpolation_size
        
        do j=1,subgrid_dim(y_dim_index)
        do i=1,subgrid_dim(x_dim_index)
        if (use_subgrid(i,j,allsource_index)) then
            
            !Assumes it is never on the edge of the EMEP grid as it is not limitted
            i_nc=crossreference_target_to_emep_subgrid(i,j,x_dim_index)
            j_nc=crossreference_target_to_emep_subgrid(i,j,y_dim_index)

            !write(*,*) i_nc,ii_start,ii_end
            !write(*,*) j_nc,jj_start,jj_end
            
            if (EMEP_projection_type.eq.LL_projection_index) then
                xpos_subgrid=lon_subgrid(i,j)
                ypos_subgrid=lat_subgrid(i,j)
            elseif (EMEP_projection_type.eq.LCC_projection_index) then
                if (use_alternative_LCC_projection_flag) then
                    call lb2lambert2_uEMEP(xpos_subgrid,ypos_subgrid,lon_subgrid(i,j),lat_subgrid(i,j),real(EMEP_projection_attributes(1)),real(EMEP_projection_attributes(2)),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                else
                    call lb2lambert_uEMEP(xpos_subgrid,ypos_subgrid,lon_subgrid(i,j),lat_subgrid(i,j),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                endif
                !call lb2lambert_uEMEP(xpos_subgrid,ypos_subgrid,lon_subgrid(i,j),lat_subgrid(i,j),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
            endif
            
            !Set the edges of the search area surounding the target grid
            lon_area_min=xpos_subgrid-lon_limit
            lon_area_max=xpos_subgrid+lon_limit
            lat_area_min=ypos_subgrid-lat_limit
            lat_area_max=ypos_subgrid+lat_limit
        
            !Use the wind direction to move the target area downwind
                if (use_downwind_position_flag.and.hourly_calculations) then
                    
                    i_cross_integral=crossreference_target_to_integral_subgrid(i,j,x_dim_index)
                    j_cross_integral=crossreference_target_to_integral_subgrid(i,j,y_dim_index)                   
                    x_pos=max(-1.,min(1.,meteo_subgrid(i_cross_integral,j_cross_integral,tt,cos_subgrid_index)*sqrt(2.)))
                    y_pos=max(-1.,min(1.,meteo_subgrid(i_cross_integral,j_cross_integral,tt,sin_subgrid_index)*sqrt(2.)))                    

                    lon_area_max=xpos_subgrid+(1.-x_pos)*lon_limit
                    lon_area_min=xpos_subgrid-(1.+x_pos)*lon_limit
                    lat_area_max=ypos_subgrid+(1.-y_pos)*lat_limit
                    lat_area_min=ypos_subgrid-(1.+y_pos)*lat_limit

                    ii_start=-1-floor(EMEP_grid_interpolation_size*0.5)
                    ii_end=1+floor(EMEP_grid_interpolation_size*0.5)
                    jj_start=-1-floor(EMEP_grid_interpolation_size*0.5)
                    jj_end=1+floor(EMEP_grid_interpolation_size*0.5)
                    
                endif


            !write(*,*) ii_start,ii_end,jj_start,jj_end,floor(0.5*(EMEP_grid_interpolation_size-1.))
        
            weighting_nc(:,:,tt_dim,:)=0.
            
            do jj=jj_start,jj_end
            do ii=ii_start,ii_end
                
                ii_nc=ii+i_nc
                jj_nc=jj+j_nc
                
                ii_w=ii+ii_w0
                jj_w=jj+jj_w0
                
                !Set the edges to an EMEP grid surorunding the EMEP grid being tested
                lon_min=max(lon_area_min,var1d_nc(ii_nc,lon_nc_index)-dgrid_nc(lon_nc_index)/2.)
                lon_max=min(lon_area_max,var1d_nc(ii_nc,lon_nc_index)+dgrid_nc(lon_nc_index)/2.)
                lat_min=max(lat_area_min,var1d_nc(jj_nc,lat_nc_index)-dgrid_nc(lat_nc_index)/2.)
                lat_max=min(lat_area_max,var1d_nc(jj_nc,lat_nc_index)+dgrid_nc(lat_nc_index)/2.)
            
                if (lon_max.gt.lon_min.and.lat_max.gt.lat_min) then
                    weighting_nc(ii_w,jj_w,tt_dim,:)=(lat_max-lat_min)*(lon_max-lon_min)/dgrid_nc(lon_nc_index)/dgrid_nc(lat_nc_index)/EMEP_grid_interpolation_size_sqr
                else
                    weighting_nc(ii_w,jj_w,tt_dim,:)=0.
                endif                

                !if (hourly_calculations) then !This if statement might slow down the loop
                !    subgrid(i,j,:,emep_local_subgrid_index,:,emep_subsource)=subgrid(i,j,:,emep_local_subgrid_index,:,emep_subsource) &
                !        +var4d_nc(ii_nc,jj_nc,surface_level_nc,:,local_nc_index,:)*weighting_nc(ii_w,jj_w,:,:)/EMEP_grid_interpolation_size_sqr
                !    subgrid(i,j,:,emep_nonlocal_subgrid_index,:,emep_subsource)=subgrid(i,j,:,emep_nonlocal_subgrid_index,:,emep_subsource) &
                !        +(var4d_nc(ii_nc,jj_nc,surface_level_nc,:,conc_nc_index,:)-var4d_nc(ii_nc,jj_nc,surface_level_nc,:,local_nc_index,:))*weighting_nc(ii_w,jj_w,:,:)/EMEP_grid_interpolation_size_sqr
                !else
                    subgrid(i,j,tt,emep_local_subgrid_index,:,emep_subsource)=subgrid(i,j,tt,emep_local_subgrid_index,:,emep_subsource) &
                        +var3d_nc(ii_nc,jj_nc,tt,local_nc_index,:)*weighting_nc(ii_w,jj_w,tt_dim,:)
                    subgrid(i,j,tt,emep_nonlocal_subgrid_index,:,emep_subsource)=subgrid(i,j,tt,emep_nonlocal_subgrid_index,:,emep_subsource) &
                        +(var3d_nc(ii_nc,jj_nc,tt,conc_nc_index,:)-var3d_nc(ii_nc,jj_nc,tt,local_nc_index,:))*weighting_nc(ii_w,jj_w,tt_dim,:)

                    !Interpolate the other EMEP compounds as well to subgrid
                    do i_comp=1,n_compound_loop
                        comp_EMEP_subgrid(i,j,tt,compound_loop_index(i_comp))=comp_EMEP_subgrid(i,j,tt,compound_loop_index(i_comp)) &
                        +comp_var3d_nc(ii_nc,jj_nc,tt,compound_loop_index(i_comp))*weighting_nc(ii_w,jj_w,tt_dim,allsource_index)
                    enddo
                    
                !endif

            enddo
            enddo
            !write(*,*) sum(weighting_nc(:,:,traffic_index,1)),sum(weighting_nc(:,:,shipping_index,1))
            !Subtract the additional local emissions from the nonlocal using the new scheme
            
            !Problems with the time and source index. Need to fix for hourly data
            nonlocal_correction(tt_dim,:)=0.
            do jj=jj_start,jj_end
            do ii=ii_start,ii_end
                
                ii_nc=ii+i_nc
                jj_nc=jj+j_nc
                
                ii_w=ii+ii_w0
                jj_w=jj+jj_w0

                !write(*,*) ii,jj,id,jd,id-ii,jd-jj
                if (jj.ne.0.or.ii.ne.0) then
                    !First weight is emission, the second is area
                    !nonlocal_correction(:,:)=nonlocal_correction(:,:) &
                    !    -lc_var3d_nc(ii_nc,jj_nc,:,NN_index(ii_w0,jj_w0),:)*weighting_nc(ii_w,jj_w,:,:)*weighting_nc(ii_w0,jj_w0,:,:) &
                    !    -lc_var3d_nc(i_nc,j_nc,:,NN_index(ii_w,jj_w),:)*weighting_nc(ii_w0,jj_w0,:,:)*weighting_nc(ii_w,jj_w,:,:)
                    nonlocal_correction(tt_dim,:)=nonlocal_correction(tt_dim,:) &
                        -lc_var3d_nc(ii_w0,jj_w0,ii_nc,jj_nc,tt,lc_local_nc_index,:)*weighting_nc(ii_w,jj_w,tt_dim,:)*weighting_nc(ii_w0,jj_w0,tt_dim,:) &
                        -lc_var3d_nc(ii_w,jj_w,i_nc,j_nc,tt,lc_local_nc_index,:)*weighting_nc(ii_w0,jj_w0,tt_dim,:)*weighting_nc(ii_w,jj_w,tt_dim,:)
                    
                    !write(*,'(i3,7es10.2,2f12.2)') NN_index(in_p,jn_p),lc_var3d_nc(i_nc,j_nc,tt,NN_index(2,2),traffic_index),subgrid(i,j,tt,emep_nonlocal_subgrid_index,traffic_index,emep_subsource), &
                    !    lc_var3d_nc(id_m,jd_m,tt,NN_index(2,2),traffic_index),lc_var3d_nc(i_nc,j_nc,tt,NN_index(in_m,jn_m),traffic_index), &
                    !    weighting_nc(id_m,jd_m,tt,1)*weighting_nc(i_nc,j_nc,tt,1),weighting_nc(i_nc,j_nc,tt,1)*weighting_nc(id_m,jd_m,tt,1), &
                    !    nonlocal_correction(tt,traffic_index),nonlocal_correction(tt,traffic_index)/lc_var3d_nc(i_nc,j_nc,tt,NN_index(2,2),traffic_index), &
                    !    nonlocal_correction(tt,traffic_index)/subgrid(i,j,tt,emep_nonlocal_subgrid_index,traffic_index,emep_subsource)
                    
                endif
                
            enddo
            enddo
            
            !Test
            !nonlocal_correction(tt,:)=0.
            
            !write(*,*) subgrid(i,j,tt,emep_local_subgrid_index,traffic_index,emep_subsource),lc_var3d_nc(i_nc,j_nc,tt,5,traffic_index),sum(lc_var3d_nc(i_nc,j_nc,tt,:,traffic_index))-lc_var3d_nc(i_nc,j_nc,tt,5,traffic_index),nonlocal_correction(traffic_index)
            
            !do source_index=1,n_source_index
            subgrid(i,j,tt,emep_nonlocal_subgrid_index,:,emep_subsource)=subgrid(i,j,tt,emep_nonlocal_subgrid_index,:,emep_subsource)+nonlocal_correction(tt_dim,:)
            subgrid(i,j,tt,emep_local_subgrid_index,:,emep_subsource)=subgrid(i,j,tt,emep_local_subgrid_index,:,emep_subsource)-nonlocal_correction(tt_dim,:)
            !enddo
            subgrid(i,j,tt,emep_subgrid_index,:,emep_subsource)=subgrid(i,j,tt,emep_nonlocal_subgrid_index,:,emep_subsource)+subgrid(i,j,tt,emep_local_subgrid_index,:,emep_subsource)

            area_weighting_nc(i,j,:,:,tt_dim,:)=weighting_nc(:,:,tt_dim,:)
            
            !For diagnostics only
            nonlocal_correction_average=nonlocal_correction_average+nonlocal_correction(tt_dim,:)
            
        
        endif
        enddo
            !write(*,*) 'Subgrid EMEP area interpolation: ',j,' of ',subgrid_dim(2)
        enddo
        
        
        if (tt.eq.t_end) then
            nonlocal_correction_average=nonlocal_correction_average/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index)
            write(unit_logfile,'(A,<n_source_index>f12.3)') 'Nonlocal correction for area weighting = ',nonlocal_correction_average
        endif
        
        
    endif

    !Loop through subgrid and carry out a subgrid weighted moving window interpolation
    if (EMEP_grid_interpolation_flag.gt.1) then

        nonlocal_correction_average=0.
        
        !This is already set in the previous call
        !n_weight=1+2*floor(EMEP_grid_interpolation_size)
   
        if (tt.eq.t_start) write(unit_logfile,'(A,2i)')'Calculating EMEP local subgrid contribution using moving window interpolation method ',EMEP_grid_interpolation_flag
        
        !if (.not.allocated(weighting_nc)) allocate (weighting_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index),n_source_index)) !EMEP grid weighting for interpolation
        if (.not.allocated(total_weighting_nc)) allocate (total_weighting_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),tt_dim,n_source_index)) !EMEP grid weighting for interpolation
        !allocate (proxy_weighting_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index),n_source_index)) !EMEP grid weighting for interpolation
        if (.not.allocated(proxy_weighting_nc)) allocate (proxy_weighting_nc(n_weight,n_weight,tt_dim,n_source_index)) !EMEP grid weighting for interpolation
        
        !Set the index offset for the local contribution based on a dimension of 5
        i_w_c=1+floor(n_weight*.5)
        j_w_c=1+floor(n_weight*.5)
        
        !subgrid(:,:,tt,emep_subgrid_index,:,:)=0
        !subgrid(:,:,tt,emep_frac_subgrid_index,:,:)=0
        subgrid(:,:,tt,emep_local_subgrid_index,:,:)=0
        
        !Set the start and end times of the loop
        !t_start=1
        !t_end=subgrid_dim(t_dim_index)
       
        !do tt=t_start,t_end

        !Emission moving window only works when the emission and subgrids are the same
            !This can be fixed here I think
        !Emission moving window adds all the subsource emissions since EMEP does not understand subsources
        !This means that combine_emission_subsources_during_dispersion must be set to true whenever using the redistribution of EMEP
        
        !Emission weighting
        if (EMEP_grid_interpolation_flag.eq.2) then
            if (.not.allocated(weighting_subgrid)) allocate (weighting_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),tt_dim,n_source_index))
            if (.not.allocated(crossreference_weighting_to_emep_subgrid)) allocate (crossreference_weighting_to_emep_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),2,n_source_index))
            weighting_subgrid(:,:,tt_dim,:)=sum(emission_subgrid(:,:,tt,:,:),4)
            weighting_subgrid_dim(:,:)=emission_subgrid_dim(1:2,:)
            crossreference_weighting_to_emep_subgrid=crossreference_emission_to_emep_subgrid
        endif
        
        if (EMEP_grid_interpolation_flag.eq.3) then

            !Aggregated emissions first on integral grid to increase speed
            if (.not.allocated(weighting_subgrid)) allocate (weighting_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),tt_dim,n_source_index))
            if (.not.allocated(crossreference_weighting_to_emep_subgrid)) allocate (crossreference_weighting_to_emep_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),2,n_source_index))
            do i_source=1,n_source_index
                crossreference_weighting_to_emep_subgrid(:,:,:,i_source)=crossreference_integral_to_emep_subgrid
                weighting_subgrid_dim(:,i_source)=integral_subgrid_dim(1:2)
            enddo
            weighting_subgrid(:,:,tt_dim,:)=0.
            do i_source=1,n_source_index
            if (calculate_source(i_source)) then
            do j=1,emission_subgrid_dim(y_dim_index,i_source)
            do i=1,emission_subgrid_dim(x_dim_index,i_source)
                i_cross=crossreference_emission_to_integral_subgrid(i,j,x_dim_index,i_source)
                j_cross=crossreference_emission_to_integral_subgrid(i,j,y_dim_index,i_source)                   
                weighting_subgrid(i_cross,j_cross,tt_dim,i_source)=weighting_subgrid(i_cross,j_cross,tt_dim,i_source)+sum(emission_subgrid(i,j,tt,i_source,:),1)
            enddo
            enddo
            endif
            enddo

        endif
        
        
        !Integral proxy weighting
        if (EMEP_grid_interpolation_flag.eq.4) then
            !Aggregated proxy on integral grid to increase speed
            if (.not.allocated(weighting_subgrid)) allocate (weighting_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),integral_subgrid_dim(t_dim_index),n_source_index))
            if (.not.allocated(crossreference_weighting_to_emep_subgrid)) allocate (crossreference_weighting_to_emep_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),2,n_source_index))
            do i_source=1,n_source_index
                crossreference_weighting_to_emep_subgrid(:,:,:,i_source)=crossreference_integral_to_emep_subgrid
                weighting_subgrid_dim(:,i_source)=integral_subgrid_dim(1:2)
            enddo
            weighting_subgrid(:,:,tt_dim,:)=sum(integral_subgrid(:,:,tt,:,:),4)
            
        endif
        
        !Calculate weighting sum for each EMEP grid.
        total_weighting_nc=0.
        do i_source=1,n_source_index
        if (calculate_source(i_source)) then
        do j=1,weighting_subgrid_dim(y_dim_index,i_source)
        do i=1,weighting_subgrid_dim(x_dim_index,i_source)
            i_nc=crossreference_weighting_to_emep_subgrid(i,j,x_dim_index,i_source)
            j_nc=crossreference_weighting_to_emep_subgrid(i,j,y_dim_index,i_source)
            total_weighting_nc(i_nc,j_nc,tt_dim,i_source)=total_weighting_nc(i_nc,j_nc,tt_dim,i_source)+weighting_subgrid(i,j,tt_dim,i_source)
            !write(*,*) i_source,i,j,i_nc,j_nc,weighting_subgrid(i,j,:,i_source)
        enddo
        enddo
        endif
        enddo

        nonlocal_correction_average=0.
        
        lon_limit=dgrid_nc(lon_nc_index)/2.*EMEP_grid_interpolation_size
        lat_limit=dgrid_nc(lat_nc_index)/2.*EMEP_grid_interpolation_size


        do i_source=1,n_source_index
        if (calculate_source(i_source)) then
        
        
        !Calculate the proxy weighting in the nearest emep grids for each subgrid
        
        do j=1,subgrid_dim(y_dim_index)
        do i=1,subgrid_dim(x_dim_index)
        if (use_subgrid(i,j,allsource_index)) then
            
            proxy_weighting_nc=0.
           
            if (EMEP_projection_type.eq.LL_projection_index) then
                xpos_subgrid=lon_subgrid(i,j)
                ypos_subgrid=lat_subgrid(i,j)
            elseif (EMEP_projection_type.eq.LCC_projection_index) then
                if (use_alternative_LCC_projection_flag) then
                    call lb2lambert2_uEMEP(xpos_subgrid,ypos_subgrid,lon_subgrid(i,j),lat_subgrid(i,j),real(EMEP_projection_attributes(1)),real(EMEP_projection_attributes(2)),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                else
                    call lb2lambert_uEMEP(xpos_subgrid,ypos_subgrid,lon_subgrid(i,j),lat_subgrid(i,j),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                endif
                !call lb2lambert_uEMEP(xpos_subgrid,ypos_subgrid,lon_subgrid(i,j),lat_subgrid(i,j),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
            endif
            
            !Set the edges of the search area surounding the target grid
            lon_area_min=xpos_subgrid-lon_limit
            lon_area_max=xpos_subgrid+lon_limit
            lat_area_min=ypos_subgrid-lat_limit
            lat_area_max=ypos_subgrid+lat_limit
            
            if (EMEP_grid_interpolation_flag.eq.2) then
                i_cross=crossreference_target_to_emission_subgrid(i,j,x_dim_index,i_source)
                j_cross=crossreference_target_to_emission_subgrid(i,j,y_dim_index,i_source)                   
                i_nc_c=crossreference_target_to_emep_subgrid(i,j,x_dim_index)
                j_nc_c=crossreference_target_to_emep_subgrid(i,j,y_dim_index)                   
                !Limit the loop so that it doesn't go over more than necessary subgrids and does not go outside the domain
                i_start=max(1,i_cross-emission_subgrid_loop_index(x_dim_index,i_source))
                i_end=min(emission_subgrid_dim(x_dim_index,i_source),i_cross+emission_subgrid_loop_index(x_dim_index,i_source))
                j_start=max(1,j_cross-emission_subgrid_loop_index(y_dim_index,i_source))
                j_end=min(emission_subgrid_dim(y_dim_index,i_source),j_cross+emission_subgrid_loop_index(y_dim_index,i_source))
                
                !Use the wind direction to move the target area downwind
                if (use_downwind_position_flag.and.hourly_calculations) then
                    
                    !Set new lon and lat limits to twice their normal size to include the upwind source region
                    !lon_limit=dgrid_nc(lon_nc_index)*EMEP_grid_interpolation_size
                    !lat_limit=dgrid_nc(lat_nc_index)*EMEP_grid_interpolation_size
                    
                    !do tt=t_start,t_end
    
                        i_cross_integral=crossreference_target_to_integral_subgrid(i,j,x_dim_index)
                        j_cross_integral=crossreference_target_to_integral_subgrid(i,j,y_dim_index)                    
                        x_pos=max(-1.,min(1.,meteo_subgrid(i_cross_integral,j_cross_integral,tt,cos_subgrid_index)*sqrt(2.)))
                        y_pos=max(-1.,min(1.,meteo_subgrid(i_cross_integral,j_cross_integral,tt,sin_subgrid_index)*sqrt(2.)))                    
                        i_end=min(ceiling(i_cross+(1.-x_pos)*emission_subgrid_loop_index(x_dim_index,i_source)),emission_subgrid_dim(x_dim_index,i_source))
                        i_start=max(floor(i_cross-(1.+x_pos)*emission_subgrid_loop_index(x_dim_index,i_source)),1)
                        j_end=min(ceiling(j_cross+(1.-y_pos)*emission_subgrid_loop_index(y_dim_index,i_source)),emission_subgrid_dim(y_dim_index,i_source))
                        j_start=max(floor(j_cross-(1.+y_pos)*emission_subgrid_loop_index(y_dim_index,i_source)),1)
                        lon_area_max=xpos_subgrid+(1.-x_pos)*lon_limit
                        lon_area_min=xpos_subgrid-(1.+x_pos)*lon_limit
                        lat_area_max=ypos_subgrid+(1.-y_pos)*lat_limit
                        lat_area_min=ypos_subgrid-(1.+y_pos)*lat_limit
                                                    
                        !write(*,*) i_start,i_end,j_start,j_end
                        
                        do jj=j_start,j_end
                        do ii=i_start,i_end                
 
                            if (EMEP_projection_type.eq.LL_projection_index) then
                                xpos_emission_subgrid=lon_emission_subgrid(ii,jj,i_source)
                                ypos_emission_subgrid=lat_emission_subgrid(ii,jj,i_source)
                            elseif (EMEP_projection_type.eq.LCC_projection_index) then
                                if (use_alternative_LCC_projection_flag) then
                                    call lb2lambert2_uEMEP(xpos_emission_subgrid,ypos_emission_subgrid,lon_emission_subgrid(ii,jj,i_source),lat_emission_subgrid(ii,jj,i_source),real(EMEP_projection_attributes(1)),real(EMEP_projection_attributes(2)),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                                else
                                    call lb2lambert_uEMEP(xpos_emission_subgrid,ypos_emission_subgrid,lon_emission_subgrid(ii,jj,i_source),lat_emission_subgrid(ii,jj,i_source),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                                endif
                                !call lb2lambert_uEMEP(xpos_emission_subgrid,ypos_emission_subgrid,lon_emission_subgrid(ii,jj,i_source),lat_emission_subgrid(ii,jj,i_source),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                            endif

                            !if (xpos_emission_subgrid.gt.lon_area_min.and.xpos_emission_subgrid.lt.lon_area_max &
                            !    .and.ypos_emission_subgrid.gt.lat_area_min.and.ypos_emission_subgrid.lt.lat_area_max) then                  
                            if (xpos_emission_subgrid.ge.lon_area_min.and.xpos_emission_subgrid.le.lon_area_max &
                                .and.ypos_emission_subgrid.ge.lat_area_min.and.ypos_emission_subgrid.le.lat_area_max) then                  
                                i_nc=crossreference_emission_to_emep_subgrid(ii,jj,x_dim_index,i_source)
                                j_nc=crossreference_emission_to_emep_subgrid(ii,jj,y_dim_index,i_source)                   
                                proxy_weighting_nc(i_nc-i_nc_c+i_w_c,j_nc-j_nc_c+j_w_c,tt_dim,i_source)=proxy_weighting_nc(i_nc-i_nc_c+i_w_c,j_nc-j_nc_c+j_w_c,tt_dim,i_source)+weighting_subgrid(ii,jj,tt_dim,i_source)
                                !write(*,*) ii,jj,proxy_weighting_nc(i_nc-i_nc_c+i_w_c,j_nc-j_nc_c+j_w_c,tt,i_source)
                            endif  
                        enddo
                        enddo
                        
                    !enddo
                
                else
                 
                    do jj=j_start,j_end
                    do ii=i_start,i_end                
 
                        if (EMEP_projection_type.eq.LL_projection_index) then
                            xpos_emission_subgrid=lon_emission_subgrid(ii,jj,i_source)
                            ypos_emission_subgrid=lat_emission_subgrid(ii,jj,i_source)
                        elseif (EMEP_projection_type.eq.LCC_projection_index) then
                            if (use_alternative_LCC_projection_flag) then
                                call lb2lambert2_uEMEP(xpos_emission_subgrid,ypos_emission_subgrid,lon_emission_subgrid(ii,jj,i_source),lat_emission_subgrid(ii,jj,i_source),real(EMEP_projection_attributes(1)),real(EMEP_projection_attributes(2)),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                            else
                                call lb2lambert_uEMEP(xpos_emission_subgrid,ypos_emission_subgrid,lon_emission_subgrid(ii,jj,i_source),lat_emission_subgrid(ii,jj,i_source),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                            endif
                            !call lb2lambert_uEMEP(xpos_emission_subgrid,ypos_emission_subgrid,lon_emission_subgrid(ii,jj,i_source),lat_emission_subgrid(ii,jj,i_source),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                        endif

                        if (abs(xpos_subgrid-xpos_emission_subgrid).le.lon_limit &
                            .and.abs(ypos_subgrid-ypos_emission_subgrid).le.lat_limit) then                  
                            i_nc=crossreference_emission_to_emep_subgrid(ii,jj,x_dim_index,i_source)
                            j_nc=crossreference_emission_to_emep_subgrid(ii,jj,y_dim_index,i_source)                   
                            !proxy_weighting_nc(i_nc,j_nc,:,i_source)=proxy_weighting_nc(i_nc,j_nc,:,i_source)+weighting_subgrid(ii,jj,:,i_source)
                            proxy_weighting_nc(i_nc-i_nc_c+i_w_c,j_nc-j_nc_c+j_w_c,tt_dim,i_source)=proxy_weighting_nc(i_nc-i_nc_c+i_w_c,j_nc-j_nc_c+j_w_c,tt_dim,i_source)+weighting_subgrid(ii,jj,tt_dim,i_source)
                            !write(*,*) tt, proxy_weighting_nc(i_nc-i_nc_c+i_w_c,j_nc-j_nc_c+j_w_c,tt_dim,i_source), weighting_subgrid(ii,jj,tt_dim,i_source)

                        endif                   
                    enddo
                    enddo
                
                endif
                
            elseif (EMEP_grid_interpolation_flag.eq.3.or.EMEP_grid_interpolation_flag.eq.4)  then
                !Find the cross reference to the integral grid from the target grid
                !do i_source=1,n_source_index
                !if (calculate_source(i_source)) then
                i_cross=crossreference_target_to_integral_subgrid(i,j,x_dim_index)
                j_cross=crossreference_target_to_integral_subgrid(i,j,y_dim_index)    
                i_nc_c=crossreference_target_to_emep_subgrid(i,j,x_dim_index)
                j_nc_c=crossreference_target_to_emep_subgrid(i,j,y_dim_index)                   
                i_start=max(1,i_cross-integral_subgrid_loop_index(x_dim_index))
                i_end=min(integral_subgrid_dim(x_dim_index),i_cross+integral_subgrid_loop_index(x_dim_index))
                j_start=max(1,j_cross-integral_subgrid_loop_index(y_dim_index))
                j_end=min(integral_subgrid_dim(y_dim_index),j_cross+integral_subgrid_loop_index(y_dim_index))
               
                if (use_downwind_position_flag.and.hourly_calculations) then
                    
                    !Set new lon and lat limits to twice their normal size to include the upwind source region
                    !lon_limit=dgrid_nc(lon_nc_index)*EMEP_grid_interpolation_size
                    !lat_limit=dgrid_nc(lat_nc_index)*EMEP_grid_interpolation_size
                    
                    !do tt=t_start,t_end
                   
                        i_cross_integral=i_cross
                        j_cross_integral=j_cross
                        x_pos=max(-1.,min(1.,meteo_subgrid(i_cross_integral,j_cross_integral,tt,cos_subgrid_index)*sqrt(2.)))
                        y_pos=max(-1.,min(1.,meteo_subgrid(i_cross_integral,j_cross_integral,tt,sin_subgrid_index)*sqrt(2.)))                    
                        i_end=min(ceiling(i_cross+(1.-x_pos)*integral_subgrid_loop_index(x_dim_index)),integral_subgrid_dim(x_dim_index))
                        i_start=max(floor(i_cross-(1.+x_pos)*integral_subgrid_loop_index(x_dim_index)),1)
                        j_end=min(ceiling(j_cross+(1.-y_pos)*integral_subgrid_loop_index(y_dim_index)),integral_subgrid_dim(y_dim_index))
                        j_start=max(floor(j_cross-(1.+y_pos)*integral_subgrid_loop_index(y_dim_index)),1)
                        lon_area_max=xpos_subgrid+(1.-x_pos)*lon_limit
                        lon_area_min=xpos_subgrid-(1.+x_pos)*lon_limit
                        lat_area_max=ypos_subgrid+(1.-y_pos)*lat_limit
                        lat_area_min=ypos_subgrid-(1.+y_pos)*lat_limit

                                                   
                        do jj=j_start,j_end
                        do ii=i_start,i_end
                            if (EMEP_projection_type.eq.LL_projection_index) then
                                xpos_integral_subgrid=lon_integral_subgrid(ii,jj)
                                ypos_integral_subgrid=lat_integral_subgrid(ii,jj)
                            elseif (EMEP_projection_type.eq.LCC_projection_index) then
                                if (use_alternative_LCC_projection_flag) then
                                    call lb2lambert2_uEMEP(xpos_integral_subgrid,ypos_integral_subgrid,lon_integral_subgrid(ii,jj),lat_integral_subgrid(ii,jj),real(EMEP_projection_attributes(1)),real(EMEP_projection_attributes(2)),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                                else
                                    call lb2lambert_uEMEP(xpos_integral_subgrid,ypos_integral_subgrid,lon_integral_subgrid(ii,jj),lat_integral_subgrid(ii,jj),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                                endif
                                !call lb2lambert_uEMEP(xpos_integral_subgrid,ypos_integral_subgrid,lon_integral_subgrid(ii,jj),lat_integral_subgrid(ii,jj),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                            endif
                            
                            !if (xpos_integral_subgrid.gt.lon_area_min.and.xpos_integral_subgrid.lt.lon_area_max &
                            !    .and.ypos_integral_subgrid.gt.lat_area_min.and.ypos_integral_subgrid.lt.lat_area_max) then                  
                            if (xpos_integral_subgrid.ge.lon_area_min.and.xpos_integral_subgrid.le.lon_area_max &
                                .and.ypos_integral_subgrid.ge.lat_area_min.and.ypos_integral_subgrid.le.lat_area_max) then                  
                                i_nc=crossreference_integral_to_emep_subgrid(ii,jj,x_dim_index)
                                j_nc=crossreference_integral_to_emep_subgrid(ii,jj,y_dim_index)                   
                                !proxy_weighting_nc(i_nc,j_nc,:,i_source)=proxy_weighting_nc(i_nc,j_nc,:,i_source)+weighting_subgrid(ii,jj,:,i_source)
                                !write(*,*) i_cross,j_cross,i_nc_c,j_nc_c,i_nc-i_nc_c+i_w_c,j_nc-j_nc_c+j_w_c
                                proxy_weighting_nc(i_nc-i_nc_c+i_w_c,j_nc-j_nc_c+j_w_c,tt_dim,i_source)=proxy_weighting_nc(i_nc-i_nc_c+i_w_c,j_nc-j_nc_c+j_w_c,tt_dim,i_source)+weighting_subgrid(ii,jj,tt_dim,i_source)
                                !write(*,*)tt,ii,jj,proxy_weighting_nc(i_nc-i_nc_c+i_w_c,j_nc-j_nc_c+j_w_c,tt,i_source)
                            endif                   
                        enddo
                        enddo

                    !enddo
                    
                else
                
                    do jj=j_start,j_end
                    do ii=i_start,i_end
                        if (EMEP_projection_type.eq.LL_projection_index) then
                            xpos_integral_subgrid=lon_integral_subgrid(ii,jj)
                            ypos_integral_subgrid=lat_integral_subgrid(ii,jj)
                        elseif (EMEP_projection_type.eq.LCC_projection_index) then
                            if (use_alternative_LCC_projection_flag) then
                                call lb2lambert2_uEMEP(xpos_integral_subgrid,ypos_integral_subgrid,lon_integral_subgrid(ii,jj),lat_integral_subgrid(ii,jj),real(EMEP_projection_attributes(1)),real(EMEP_projection_attributes(2)),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                            else
                                call lb2lambert_uEMEP(xpos_integral_subgrid,ypos_integral_subgrid,lon_integral_subgrid(ii,jj),lat_integral_subgrid(ii,jj),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                            endif
                            !call lb2lambert_uEMEP(xpos_integral_subgrid,ypos_integral_subgrid,lon_integral_subgrid(ii,jj),lat_integral_subgrid(ii,jj),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                        endif

                        if (abs(xpos_subgrid-xpos_integral_subgrid).le.lon_limit &
                            .and.abs(ypos_subgrid-ypos_integral_subgrid).le.lat_limit) then                  
                            i_nc=crossreference_integral_to_emep_subgrid(ii,jj,x_dim_index)
                            j_nc=crossreference_integral_to_emep_subgrid(ii,jj,y_dim_index)                   
                            !proxy_weighting_nc(i_nc,j_nc,:,i_source)=proxy_weighting_nc(i_nc,j_nc,:,i_source)+weighting_subgrid(ii,jj,:,i_source)
                            !write(*,*) i_cross,j_cross,i_nc_c,j_nc_c,i_nc-i_nc_c+i_w_c,j_nc-j_nc_c+j_w_c
                            proxy_weighting_nc(i_nc-i_nc_c+i_w_c,j_nc-j_nc_c+j_w_c,tt_dim,i_source)=proxy_weighting_nc(i_nc-i_nc_c+i_w_c,j_nc-j_nc_c+j_w_c,tt_dim,i_source)+weighting_subgrid(ii,jj,tt_dim,i_source)
                            !write(*,*)ii,jj,proxy_weighting_nc(i_nc-i_nc_c+i_w_c,j_nc-j_nc_c+j_w_c,1,i_source)
                        endif                   
                    enddo
                    enddo
                
                endif
                
                !endif
                !enddo
                !do ii=2,4
                !do jj=2,4
                    !write(*,*) i,j,ii,jj,proxy_weighting_nc(ii,jj,:,i_source)
                !enddo
                !enddo
                        
            endif
        
 
            
            !if (weighting_nc(ii,jj,traffic_index).ne.0) then
            !write(*,*) i,j,weighting_nc(ii,jj,traffic_index)
            !endif
            
            !Normalise weighting
            
            !Specify the nearest to the target grid
            !if (EMEP_grid_interpolation_flag.eq.2) then
            !    i_cross=crossreference_target_to_emission_subgrid(i,j,x_dim_index,i_source)
            !    j_cross=crossreference_target_to_emission_subgrid(i,j,y_dim_index,i_source)    
            !elseif (EMEP_grid_interpolation_flag.eq.3.or.EMEP_grid_interpolation_flag.eq.4)  then
            !    i_cross=crossreference_target_to_integral_subgrid(i,j,x_dim_index)
            !    j_cross=crossreference_target_to_integral_subgrid(i,j,y_dim_index)    
            !endif
            !i_nc=crossreference_weighting_to_emep_subgrid(i_cross,j_cross,x_dim_index,i_source)
            !j_nc=crossreference_weighting_to_emep_subgrid(i_cross,j_cross,y_dim_index,i_source)
            i_nc=crossreference_target_to_emep_subgrid(i,j,x_dim_index)
            j_nc=crossreference_target_to_emep_subgrid(i,j,y_dim_index)
            i_nc_start=max(1+i_nc-i_w_c,i_nc-1-floor(EMEP_grid_interpolation_size*0.5))
            i_nc_end=min(dim_length_nc(x_dim_nc_index)+i_nc-i_w_c,i_nc+1+floor(EMEP_grid_interpolation_size*0.5))
            j_nc_start=max(1+j_nc-j_w_c,j_nc-1-floor(EMEP_grid_interpolation_size*0.5))
            j_nc_end=min(dim_length_nc(y_dim_nc_index)+j_nc-j_w_c,j_nc+1+floor(EMEP_grid_interpolation_size*0.5))

            !weighting_nc=weighting_nc/total_weighting_nc
            !i_nc_start>1+i_nc-i_w_c
            !i_nc_end<dim_length_nc(x_dim_nc_index)+i_nc-i_w_c
            
            !if (2.eq.1) then

            do jj=j_nc_start,j_nc_end
            do ii=i_nc_start,i_nc_end
            
            !Can take out the time loop here if we do not have the if statement. Include a small value instead
            !do tt=t_start,t_end
            !do i_source=1,n_source_index
                if (total_weighting_nc(ii,jj,tt_dim,i_source).ne.0.) then
                    !proxy_weighting_nc(ii,jj,tt,i_source)=proxy_weighting_nc(ii,jj,tt,i_source)/total_weighting_nc(ii,jj,tt,i_source)
                    proxy_weighting_nc(ii-i_nc+i_w_c,jj-j_nc+j_w_c,tt_dim,i_source)=proxy_weighting_nc(ii-i_nc+i_w_c,jj-j_nc+j_w_c,tt_dim,i_source)/total_weighting_nc(ii,jj,tt_dim,i_source)/EMEP_grid_interpolation_size_sqr
                !write(*,*)  tt,proxy_weighting_nc(ii-i_nc+i_w_c,jj-j_nc+j_w_c,tt_dim,i_source)
                else
                    !proxy_weighting_nc(ii,jj,tt,i_source)=0.
                    proxy_weighting_nc(ii-i_nc+i_w_c,jj-j_nc+j_w_c,tt_dim,i_source)=0.
                endif
            !enddo
                if (proxy_weighting_nc(ii-i_nc+i_w_c,jj-j_nc+j_w_c,tt_dim,i_source).gt.1) write(*,'(A,7i6,f12.2)')'WEIGHTING>1: ',tt,i,j,ii,jj,ii-i_nc,jj-j_nc,proxy_weighting_nc(ii-i_nc+i_w_c,jj-j_nc+j_w_c,tt_dim,i_source)
            !enddo
            
            !write(*,*) i,j,ii-i_nc+i_w_c,jj-j_nc+j_w_c,proxy_weighting_nc(ii-i_nc+i_w_c,jj-j_nc+j_w_c,:,i_source)
            
            enddo
            enddo
                    
            !endif
            
            !Add up the contributing weights
            !Note that only the emep_subsource=1 can be determined since EMEP has no subsources
            !if (hourly_calculations) then
                !do i_subsource=1,n_subsource_index(i_source)
                do jj=j_nc_start,j_nc_end
                do ii=i_nc_start,i_nc_end
                    
                    subgrid(i,j,tt,emep_local_subgrid_index,i_source,emep_subsource)=subgrid(i,j,tt,emep_local_subgrid_index,i_source,emep_subsource) &
                        +var3d_nc(ii,jj,tt,local_nc_index,i_source)*proxy_weighting_nc(ii-i_nc+i_w_c,jj-j_nc+j_w_c,tt_dim,i_source)
                    !if (weighting_nc(ii,jj,traffic_index).ne.0) then
                        !write(*,'(4i6,2es12.2)') i,j,ii,jj,proxy_weighting_nc(ii-i_nc+i_w_c,jj-j_nc+j_w_c,tt_dim,i_source),subgrid(i,j,tt,emep_local_subgrid_index,i_source,emep_subsource)
                    !endif
                    
                enddo
                enddo
                !enddo
            
                !write(*,*) subgrid(i,j,1,emep_local_subgrid_index,i_source,emep_subsource),subgrid(i,j,2,emep_local_subgrid_index,i_source,emep_subsource),var3d_nc(i_nc,j_nc,1,local_nc_index,i_source),var3d_nc(i_nc,j_nc,2,local_nc_index,i_source),proxy_weighting_nc(i_nc-i_nc+i_w_c,j_nc-j_nc+j_w_c,1,i_source),proxy_weighting_nc(i_nc-i_nc+i_w_c,j_nc-j_nc+j_w_c,2,i_source)
                !write(*,*) tt,i,j,subgrid(i,j,tt,emep_local_subgrid_index,i_source,emep_subsource)
                
           ! else
                !Slightly quicker method, perhaps, for annual values
                !do i_source=1,n_source_index
                !    subgrid(i,j,1,emep_local_subgrid_index,i_source,emep_subsource)=sum(var3d_nc(i_nc_start:i_nc_end,j_nc_start:j_nc_end,1,local_nc_index,i_source) &
                !        *weighting_nc(i_nc_start:i_nc_end,j_nc_start:j_nc_end,i_source))
                !enddo
                !do i_subsource=1,n_subsource_index(i_source)
                !if (2.eq.1) then
             !   tt=1
                
             !   do jj=j_nc_start,j_nc_end
             !   do ii=i_nc_start,i_nc_end
             !       subgrid(i,j,:,emep_local_subgrid_index,i_source,emep_subsource)=subgrid(i,j,:,emep_local_subgrid_index,i_source,emep_subsource) &
             !           +var3d_nc(ii,jj,:,local_nc_index,i_source)*proxy_weighting_nc(ii-i_nc+i_w_c,jj-j_nc+j_w_c,:,i_source)/EMEP_grid_interpolation_size_sqr
             !   enddo
             !   enddo
                !enddo
                !endif
           ! endif

            !Subtract the additional local emissions from the nonlocal using the new scheme
            !i_nc=crossreference_target_to_emep_subgrid(i,j,x_dim_index)
            !j_nc=crossreference_target_to_emep_subgrid(i,j,y_dim_index)
            !Assumes never on edge of EMEP grid so that i_nc and j_nc +/-1 are within the bounds. Should set in a control for this
            nonlocal_correction(tt_dim,i_source)=0.
            !if (2.eq.1) then
            do jj=-1-floor(EMEP_grid_interpolation_size*0.5),+1+floor(EMEP_grid_interpolation_size*0.5)
            do ii=-1-floor(EMEP_grid_interpolation_size*0.5),+1+floor(EMEP_grid_interpolation_size*0.5)
                ii_nc=ii+i_nc
                jj_nc=jj+j_nc
                
                ii_w=ii+ii_w0
                jj_w=jj+jj_w0

                if (jj.ne.0.or.ii.ne.0) then
                    !First weight is emission, the second is area
                    nonlocal_correction(tt_dim,i_source)=nonlocal_correction(tt_dim,i_source) &
                        -lc_var3d_nc(ii_w0,jj_w0,ii_nc,jj_nc,tt,lc_local_nc_index,i_source)*proxy_weighting_nc(ii+i_w_c,jj+j_w_c,tt_dim,i_source)*area_weighting_nc(i,j,ii_w0,jj_w0,tt_dim,i_source) &
                        -lc_var3d_nc(ii_w,jj_w,i_nc,j_nc,tt,lc_local_nc_index,i_source)*proxy_weighting_nc(ii+i_w_c,jj+j_w_c,tt_dim,i_source)*area_weighting_nc(i,j,ii_w,jj_w,tt_dim,i_source)
                        !-lc_var3d_nc(ii_w,jj_w,i_nc,j_nc,:,lc_local_nc_index,i_source)*proxy_weighting_nc(i_w_c,j_w_c,:,i_source)*area_weighting_nc(i,j,ii_w,jj_w,:,i_source)
                endif
                        !-lc_var3d_nc(ii_w0,jj_w0,ii_nc,jj_nc,:,lc_local_nc_index,:)*weighting_nc(ii_w,jj_w,:,:)*weighting_nc(ii_w0,jj_w0,:,:) &
                        !-lc_var3d_nc(ii_w,jj_w,i_nc,j_nc,:,lc_local_nc_index,:)*weighting_nc(ii_w0,jj_w0,:,:)*weighting_nc(ii_w,jj_w,:,:)


            enddo
            enddo
            !endif
            !Adjust correction for size of area
            nonlocal_correction(tt_dim,i_source)=nonlocal_correction(tt_dim,i_source)
            
            !nonlocal_correction=0.
            
            subgrid(i,j,tt,emep_nonlocal_subgrid_index,i_source,emep_subsource)=subgrid(i,j,tt,emep_nonlocal_subgrid_index,i_source,emep_subsource)+nonlocal_correction(tt_dim,i_source)
            subgrid(i,j,tt,emep_local_subgrid_index,i_source,emep_subsource)=subgrid(i,j,tt,emep_local_subgrid_index,i_source,emep_subsource)-nonlocal_correction(tt_dim,i_source)            
            subgrid(i,j,tt,emep_subgrid_index,i_source,emep_subsource)=subgrid(i,j,tt,emep_nonlocal_subgrid_index,i_source,emep_subsource)+subgrid(i,j,tt,emep_local_subgrid_index,i_source,emep_subsource)
            
            !Averaged over time for diagnostic purposes only
            nonlocal_correction_average(i_source)=nonlocal_correction_average(i_source)+nonlocal_correction(tt_dim,i_source)
               
        !write(*,*) subgrid(i,j,shipping_index,emep_local_subgrid_index),sum(weighting_nc(:,:,shipping_index)),subgrid(i,j,shipping_index,proxy_integral_subgrid_index)
        !endif 
            !write (*,*) tt,i,j,subgrid(i,j,tt,emep_local_subgrid_index,i_source,emep_subsource) !Has a value
        endif !use subgrid
        
        enddo
        !i_source=allsource_index
        !if (mod(j,10).eq.0) write(*,'(A,A,A,i4,A,i4)') 'Subgrid EMEP with proxy interpolation for ',trim(source_file_str(i_source)),': ',j,' of ',subgrid_dim(2)
        !if (mod(i,10).eq.0) write(*,'(A,A,A,i4,A,i4)') 'Subgrid EMEP with proxy interpolation for ',trim(source_file_str(i_source)),': ',i,' of ',subgrid_dim(1)
        enddo
        
        
        endif !End if calculate_source
        enddo !End source loop

        !enddo !End time loop
        
        if (tt.eq.t_end) then
            nonlocal_correction_average=nonlocal_correction_average/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index)
            write(unit_logfile,'(A,<n_source_index>f12.3)') 'Nonlocal correction for proxy weighting = ',nonlocal_correction_average
        endif

    endif
         
    if (mod(j,1).eq.0) write(*,'(a,i5,a,i5)') 'Gridding EMEP for hour ',tt,' of ',subgrid_dim(t_dim_index)
   
    enddo   !End time loop

        !Create the all source version of the local and nonlocal contribution after calculating all the source contributions
        !The nonlocal contribution uses the difference between the local and total, here the total is based on the area interpolation. Is this correct?
        subgrid(:,:,:,emep_local_subgrid_index,allsource_index,emep_subsource)=0.
        subgrid(:,:,:,emep_subgrid_index,allsource_index,emep_subsource)=0.
        subgrid(:,:,:,emep_nonlocal_subgrid_index,allsource_index,emep_subsource)=0. !-subgrid(:,:,:,emep_subgrid_index,allsource_index,emep_subsource)
        count=0
        do i_source=1,n_source_index
        if (calculate_source(i_source)) then
                
            !Check values for local and totals for each source
            !write(*,*) trim(source_file_str(i_source))
            if (minval(subgrid(:,:,:,emep_nonlocal_subgrid_index,i_source,emep_subsource)).lt.0.0) then
                write(unit_logfile,'(A,A,f12.4,A)') 'WARNING: Min nonlocal source less than 0 for ',trim(source_file_str(i_source)),minval(subgrid(:,:,:,emep_nonlocal_subgrid_index,i_source,emep_subsource)),' Setting to 0 and adding to local'            
            endif
            !Set any negative nonlocal to 0 and add the value back into the local. Indicates a problem with the moving window method
            !Done in a loop because of stack overflow problem
                do j=1,subgrid_dim(y_dim_index)
                do i=1,subgrid_dim(x_dim_index)
                    where (subgrid(i,j,:,emep_nonlocal_subgrid_index,i_source,emep_subsource).lt.0.) 
                        subgrid(i,j,:,emep_nonlocal_subgrid_index,i_source,emep_subsource)=0.
                        subgrid(i,j,:,emep_local_subgrid_index,i_source,emep_subsource)=subgrid(i,j,:,emep_local_subgrid_index,i_source,emep_subsource)-subgrid(i,j,:,emep_nonlocal_subgrid_index,i_source,emep_subsource)
                    endwhere
                enddo
                enddo
                
                !Add the local subgrid sources together to get an allsource local contribution
                subgrid(:,:,:,emep_local_subgrid_index,allsource_index,emep_subsource)=subgrid(:,:,:,emep_local_subgrid_index,allsource_index,emep_subsource)+subgrid(:,:,:,emep_local_subgrid_index,i_source,emep_subsource)
                !subgrid(:,:,:,emep_nonlocal_subgrid_index,allsource_index,emep_subsource)=subgrid(:,:,:,emep_nonlocal_subgrid_index,allsource_index,emep_subsource)+subgrid(:,:,:,emep_nonlocal_subgrid_index,i_source,emep_subsource)
                subgrid(:,:,:,emep_subgrid_index,i_source,emep_subsource)=subgrid(:,:,:,emep_nonlocal_subgrid_index,i_source,emep_subsource)+subgrid(:,:,:,emep_local_subgrid_index,i_source,emep_subsource)
                count=count+1
                !Add up the total EMEP for all source (will be averaged with count)
                subgrid(:,:,:,emep_subgrid_index,allsource_index,emep_subsource)=subgrid(:,:,:,emep_subgrid_index,allsource_index,emep_subsource)+subgrid(:,:,:,emep_subgrid_index,i_source,emep_subsource)
                if (minval(subgrid(:,:,:,emep_subgrid_index,i_source,emep_subsource)).lt.0.0) then
                    write(unit_logfile,'(A,A,f12.4)') 'ERROR: Min total source less than 0 for ',trim(source_file_str(i_source)),minval(subgrid(:,:,:,emep_subgrid_index,i_source,emep_subsource))
                    stop
                endif
                
        
        endif
        enddo
        
        !NOW WORKING FOR WEIGHTING 0 BUT LOST THE NONLOCAL COMPONENT WITH AREA WEIGHTING 1!FIX!
        
        
        !Set the allsource nonlocal value to the average of the remainder. This can be negative
        !subgrid(:,:,:,emep_nonlocal_subgrid_index,allsource_index,emep_subsource)=(subgrid(:,:,:,emep_nonlocal_subgrid_index,allsource_index,emep_subsource)-subgrid(:,:,:,emep_local_subgrid_index,allsource_index,emep_subsource))/count
        !subgrid(:,:,:,emep_nonlocal_subgrid_index,allsource_index,emep_subsource)=(subgrid(:,:,:,emep_subgrid_index,allsource_index,emep_subsource)-subgrid(:,:,:,emep_local_subgrid_index,allsource_index,emep_subsource))/count
        subgrid(:,:,:,emep_nonlocal_subgrid_index,allsource_index,emep_subsource)=(subgrid(:,:,:,emep_subgrid_index,allsource_index,emep_subsource)/count-subgrid(:,:,:,emep_local_subgrid_index,allsource_index,emep_subsource))
        
        if (minval(subgrid(:,:,:,emep_nonlocal_subgrid_index,allsource_index,emep_subsource)).lt.0.0) then
            write(unit_logfile,'(A,f12.4,A)') 'WARNING: Min nonlocal allsource less than 0 with ',minval(subgrid(:,:,:,emep_nonlocal_subgrid_index,allsource_index,emep_subsource)),' Setting to 0'  
        endif
        !Remove any negative values. Loop to avoid stack overflow
        do j=1,subgrid_dim(y_dim_index)
        do i=1,subgrid_dim(x_dim_index)
            where (subgrid(i,j,:,emep_nonlocal_subgrid_index,allsource_index,emep_subsource).lt.0.) 
                subgrid(i,j,:,emep_nonlocal_subgrid_index,allsource_index,emep_subsource)=0.
            endwhere
        enddo
        enddo

        !Add up the sources and calculate fractions
        subgrid(:,:,:,emep_subgrid_index,:,emep_subsource)=subgrid(:,:,:,emep_nonlocal_subgrid_index,:,emep_subsource)+subgrid(:,:,:,emep_local_subgrid_index,:,emep_subsource)
        subgrid(:,:,:,emep_frac_subgrid_index,:,emep_subsource)=subgrid(:,:,:,emep_local_subgrid_index,:,emep_subsource)/subgrid(:,:,:,emep_subgrid_index,:,emep_subsource)
        if (minval(subgrid(:,:,:,emep_subgrid_index,allsource_index,emep_subsource)).lt.0.0) then
             write(unit_logfile,'(A,f12.4)')'ERROR: Minimum total allsource less than 0 with ',minval(subgrid(:,:,:,emep_subgrid_index,allsource_index,emep_subsource))
            stop
        endif
    
        !do tt=t_start,t_end
        !do j=1,subgrid_dim(y_dim_index)
        !do i=1,subgrid_dim(x_dim_index)
        !        write (*,*) tt,i,j,subgrid(i,j,tt,emep_local_subgrid_index,traffic_index,emep_subsource) !Fixed
        !enddo
        !enddo
        !enddo

    !Save files
    if (save_intermediate_files) then
    do source_index=1,n_source_index
        if (.not.read_existing_grid_data(emep_subgrid_file_index(source_index)).and.(calculate_source(source_index).or.source_index.eq.allsource_index)) then
            !do subsource_index=1,n_subsource(source_index)
            i_file=emep_subgrid_file_index(source_index)
            temp_name=trim(pathname_grid(i_file))//trim(filename_grid(i_file))//'_'//trim(file_tag)//'.asc'
            write(unit_logfile,'(a)')'Writing to: '//trim(temp_name)
            call write_esri_ascii_3d_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),subgrid_delta(1),subgrid(:,:,:,emep_subgrid_index,source_index,emep_subsource),x_subgrid,y_subgrid)
            !i_file=emep_subgrid_frac_file_index(source_index)
            !temp_name=trim(pathname_grid(i_file))//trim(filename_grid(i_file))//'_'//trim(file_tag)//'.asc'
            !write(unit_logfile,'(a)')'Writing to: '//trim(temp_name)
            !call write_esri_ascii_3d_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),subgrid_delta(1),subgrid(:,:,:,emep_frac_subgrid_index,source_index,emep_subsource),x_subgrid,y_subgrid)
            i_file=emep_subgrid_local_file_index(source_index)
            temp_name=trim(pathname_grid(i_file))//trim(filename_grid(i_file))//'_'//trim(file_tag)//'.asc'
            write(unit_logfile,'(a)')'Writing to: '//trim(temp_name)
            call write_esri_ascii_3d_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),subgrid_delta(1),subgrid(:,:,:,emep_local_subgrid_index,source_index,emep_subsource),x_subgrid,y_subgrid)
            i_file=emep_subgrid_nonlocal_file_index(source_index)
            temp_name=trim(pathname_grid(i_file))//trim(filename_grid(i_file))//'_'//trim(file_tag)//'.asc'
            write(unit_logfile,'(a)')'Writing to: '//trim(temp_name)
            call write_esri_ascii_3d_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),subgrid_delta(1),subgrid(:,:,:,emep_nonlocal_subgrid_index,source_index,emep_subsource),x_subgrid,y_subgrid)
        
                    !write(*,*) i_source
            !do jj=subgrid_dim(y_dim_index),1,-1
            !write(*,'(20es12.3)') (subgrid(ii,jj,1,i_source,emep_subgrid_index),ii=1,subgrid_dim(x_dim_index))
            !enddo
            !enddo
        endif
    enddo
    endif

    
    if (save_intermediate_files) then
         if (.not.read_existing_grid_data(emep_subgrid_file_index(allsource_index))) then
            !Save the other EMEP compounds used for nox chemistry as well
            do i_comp=1,n_compound_loop
                i_file=emep_subgrid_file_index(allsource_index)
                temp_name=trim(pathname_grid(i_file))//trim(filename_grid(i_file))//'_'//trim(var_name_nc(conc_nc_index,compound_loop_index(i_comp),allsource_index))//'_'//trim(file_tag)//'.asc'
                write(unit_logfile,'(a)')'Writing to: '//trim(temp_name)
                call write_esri_ascii_3d_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),subgrid_delta(1),comp_EMEP_subgrid(:,:,:,compound_loop_index(i_comp)),x_subgrid,y_subgrid)
            enddo
         endif
     endif    


    if (allocated(weighting_nc)) deallocate(weighting_nc)
    if (allocated(area_weighting_nc)) deallocate(area_weighting_nc)
    if (allocated(total_weighting_nc)) deallocate(total_weighting_nc)
    if (allocated(weighting_subgrid)) deallocate(weighting_subgrid)
  
  
    end subroutine uEMEP_subgrid_EMEP

