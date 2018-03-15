!==========================================================================
!   uEMEP model subgrid_meteo_EMEP
!==========================================================================
!   This routine interpolates EMEP meteo data to the integral subgrid 
!   using either nearest neighbour or area weighted interpolation
!==========================================================================
    
    subroutine uEMEP_subgrid_meteo_EMEP
      
    use uEMEP_definitions

    implicit none
    
    character(256) temp_name
    logical exists
    integer ii,jj,tt
    integer i_temp,j_temp,i_source,i_file
    integer i_nc_temp,j_nc_temp
    real, allocatable :: weighting_nc(:,:)
    integer i_nc_start,i_nc_end,j_nc_start,j_nc_end
    integer i_start,i_end,j_start,j_end,t_start,t_end
    real lon_min,lon_max,lat_min,lat_max
    integer i_nc,j_nc
    real angle_utm,angle_lcc,angle_utm2
    real, allocatable :: u_utm(:),v_utm(:),th(:),ff(:)
    real dlatx,dlaty
    real xpos_subgrid,ypos_subgrid
    
    if (.not.allocated(u_utm)) allocate (u_utm(dim_length_nc(time_dim_nc_index)))
    if (.not.allocated(v_utm)) allocate (v_utm(dim_length_nc(time_dim_nc_index)))
    if (.not.allocated(th)) allocate (th(dim_length_nc(time_dim_nc_index)))
    if (.not.allocated(ff)) allocate (ff(dim_length_nc(time_dim_nc_index)))

    if (use_single_time_loop_flag) then
        if (t_loop.gt.start_time_loop_index) then
            last_meteo_subgrid(:,:,:)=meteo_subgrid(:,:,1,:)        
        endif
    endif
    

    meteo_subgrid=0.
    
    !If EMEP meteo gridded to subgrid files already exist then read them in and leave the subroutine
    !Only tests ugrid_file_index
    if (read_existing_grid_data(subgrid_meteo_file_index)) then
            
        !Test existence of the uwind file only
        i_file=subgrid_ugrid_file_index
        temp_name=trim(pathname_grid(i_file))//trim(filename_grid(i_file))//'_'//trim(file_tag)//'.asc'
        inquire(file=temp_name,exist=exists)
        if (.not.exists) then
            write(unit_logfile,*)'WARNING: '//trim(temp_name)//' does not exist.'
                !return
        else
            call read_esri_ascii_3d_file(unit_logfile,temp_name,integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),integral_subgrid_dim(t_dim_index),integral_subgrid_delta(1),meteo_subgrid(:,:,:,ugrid_subgrid_index),x_integral_subgrid,y_integral_subgrid)
            i_file=subgrid_vgrid_file_index;temp_name=trim(pathname_grid(i_file))//trim(filename_grid(i_file))//'_'//trim(file_tag)//'.asc'
            call read_esri_ascii_3d_file(unit_logfile,temp_name,integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),integral_subgrid_dim(t_dim_index),integral_subgrid_delta(1),meteo_subgrid(:,:,:,vgrid_subgrid_index),x_integral_subgrid,y_integral_subgrid)
            i_file=subgrid_hmix_file_index;temp_name=trim(pathname_grid(i_file))//trim(filename_grid(i_file))//'_'//trim(file_tag)//'.asc'
            call read_esri_ascii_3d_file(unit_logfile,temp_name,integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),integral_subgrid_dim(t_dim_index),integral_subgrid_delta(1),meteo_subgrid(:,:,:,hmix_subgrid_index),x_integral_subgrid,y_integral_subgrid)
            i_file=subgrid_kz_file_index;temp_name=trim(pathname_grid(i_file))//trim(filename_grid(i_file))//'_'//trim(file_tag)//'.asc'
            call read_esri_ascii_3d_file(unit_logfile,temp_name,integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),integral_subgrid_dim(t_dim_index),integral_subgrid_delta(1),meteo_subgrid(:,:,:,kz_subgrid_index),x_integral_subgrid,y_integral_subgrid)
            i_file=subgrid_FFgrid_file_index;temp_name=trim(pathname_grid(i_file))//trim(filename_grid(i_file))//'_'//trim(file_tag)//'.asc'
            call read_esri_ascii_3d_file(unit_logfile,temp_name,integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),integral_subgrid_dim(t_dim_index),integral_subgrid_delta(1),meteo_subgrid(:,:,:,FFgrid_subgrid_index),x_integral_subgrid,y_integral_subgrid)
            i_file=subgrid_FF10_file_index;temp_name=trim(pathname_grid(i_file))//trim(filename_grid(i_file))//'_'//trim(file_tag)//'.asc'
            call read_esri_ascii_3d_file(unit_logfile,temp_name,integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),integral_subgrid_dim(t_dim_index),integral_subgrid_delta(1),meteo_subgrid(:,:,:,FF10_subgrid_index),x_integral_subgrid,y_integral_subgrid)
        endif
        
        return
    endif

    write(unit_logfile,'(A)')'Calculating EMEP subgrid meteo data'

    !Loop through the integral subgrid and find those subgrids within EMEP grids and allocate values directly from EMEP grids. Nearest neighbour
    if (EMEP_meteo_grid_interpolation_flag.eq.0) then
        do j=1,integral_subgrid_dim(y_dim_index)
        do i=1,integral_subgrid_dim(x_dim_index)
        
            !Assumes it is never on the edge of the EMEP grid, not limitted
            i_nc=crossreference_integral_to_emep_subgrid(i,j,x_dim_index)
            j_nc=crossreference_integral_to_emep_subgrid(i,j,y_dim_index)
        
            meteo_subgrid(i,j,:,ugrid_subgrid_index)=var4d_nc(i_nc,j_nc,surface_level_nc,:,ugrid_nc_index,allsource_index)
            meteo_subgrid(i,j,:,vgrid_subgrid_index)=var4d_nc(i_nc,j_nc,surface_level_nc,:,vgrid_nc_index,allsource_index)
            meteo_subgrid(i,j,:,FFgrid_subgrid_index)=var4d_nc(i_nc,j_nc,surface_level_nc,:,FFgrid_nc_index,allsource_index)
            meteo_subgrid(i,j,:,FF10_subgrid_index)=var3d_nc(i_nc,j_nc,:,FF10_nc_index,allsource_index)
            meteo_subgrid(i,j,:,kz_subgrid_index)=var4d_nc(i_nc,j_nc,surface_level_nc,:,kz_nc_index,allsource_index)
            meteo_subgrid(i,j,:,hmix_subgrid_index)=var3d_nc(i_nc,j_nc,:,hmix_nc_index,allsource_index)
            meteo_subgrid(i,j,:,logz0_subgrid_index)=var3d_nc(i_nc,j_nc,:,logz0_nc_index,allsource_index)
            meteo_subgrid(i,j,:,invL_subgrid_index)=var3d_nc(i_nc,j_nc,:,invL_nc_index,allsource_index)
            meteo_subgrid(i,j,:,inv_FFgrid_subgrid_index)=var4d_nc(i_nc,j_nc,surface_level_nc,:,inv_FFgrid_nc_index,allsource_index)
            meteo_subgrid(i,j,:,inv_FF10_subgrid_index)=var3d_nc(i_nc,j_nc,:,inv_FF10_nc_index,allsource_index)
            meteo_subgrid(i,j,:,ustar_subgrid_index)=var3d_nc(i_nc,j_nc,:,ustar_nc_index,allsource_index)
            meteo_subgrid(i,j,:,J_subgrid_index)=var4d_nc(i_nc,j_nc,surface_level_nc,:,J_nc_index,allsource_index)
                  
        enddo
        enddo
    endif
    
    !Area weighted interpolation of meteorology to integral grid
    if (EMEP_meteo_grid_interpolation_flag.ge.1) then
          
        allocate (weighting_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index))) !EMEP grid weighting for interpolation. Does not need a source index for area weighting
 
        do j=1,integral_subgrid_dim(y_dim_index)
        do i=1,integral_subgrid_dim(x_dim_index)
        
            !Assumes it is never on the edge of the EMEP grid, not limitted
            i_nc=crossreference_integral_to_emep_subgrid(i,j,x_dim_index)
            j_nc=crossreference_integral_to_emep_subgrid(i,j,y_dim_index)
            
            if (EMEP_projection_type.eq.LL_projection_index) then
                xpos_subgrid=lon_integral_subgrid(i,j)
                ypos_subgrid=lat_integral_subgrid(i,j)
            elseif (EMEP_projection_type.eq.LCC_projection_index) then
                if (use_alternative_LCC_projection_flag) then
                    call lb2lambert2_uEMEP(xpos_subgrid,ypos_subgrid,lon_integral_subgrid(i,j),lat_integral_subgrid(i,j),real(EMEP_projection_attributes(1)),real(EMEP_projection_attributes(2)),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                else
                    call lb2lambert_uEMEP(xpos_subgrid,ypos_subgrid,lon_integral_subgrid(i,j),lat_integral_subgrid(i,j),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
                endif
                !call lb2lambert_uEMEP(xpos_subgrid,ypos_subgrid,lon_integral_subgrid(i,j),lat_integral_subgrid(i,j),real(EMEP_projection_attributes(3)),real(EMEP_projection_attributes(4)))
            endif
            
        
            do jj=j_nc-1,j_nc+1
            do ii=i_nc-1,i_nc+1
                
                lon_min=max(xpos_subgrid-dgrid_nc(lon_nc_index)/2.,var1d_nc(ii,lon_nc_index)-dgrid_nc(lon_nc_index)/2.)
                lon_max=min(xpos_subgrid+dgrid_nc(lon_nc_index)/2.,var1d_nc(ii,lon_nc_index)+dgrid_nc(lon_nc_index)/2.)
                lat_min=max(ypos_subgrid-dgrid_nc(lat_nc_index)/2.,var1d_nc(jj,lat_nc_index)-dgrid_nc(lat_nc_index)/2.)
                lat_max=min(ypos_subgrid+dgrid_nc(lat_nc_index)/2.,var1d_nc(jj,lat_nc_index)+dgrid_nc(lat_nc_index)/2.)
            
                !write(*,*) lon_min,lon_max,lat_min,lat_max
                !write(*,*) xpos_subgrid,var1d_nc(ii,lon_nc_index),dgrid_nc(lon_nc_index)
                
                !Determine the area intersection of the EMEP grid and an EMEP grid size centred on the integral subgrid
                if (lon_max.gt.lon_min.and.lat_max.gt.lat_min) then
                    weighting_nc(ii,jj)=(lat_max-lat_min)*(lon_max-lon_min)/dgrid_nc(lon_nc_index)/dgrid_nc(lat_nc_index)
                else
                    weighting_nc(ii,jj)=0.
                endif                

                !write(*,*) i,j,ii-i_nc,jj-j_nc,weighting_nc(ii,jj)
                
                meteo_subgrid(i,j,:,ugrid_subgrid_index)=meteo_subgrid(i,j,:,ugrid_subgrid_index)+var4d_nc(ii,jj,surface_level_nc,:,ugrid_nc_index,allsource_index)*weighting_nc(ii,jj)
                meteo_subgrid(i,j,:,vgrid_subgrid_index)=meteo_subgrid(i,j,:,vgrid_subgrid_index)+var4d_nc(ii,jj,surface_level_nc,:,vgrid_nc_index,allsource_index)*weighting_nc(ii,jj)
                meteo_subgrid(i,j,:,FFgrid_subgrid_index)=meteo_subgrid(i,j,:,FFgrid_subgrid_index)+var4d_nc(ii,jj,surface_level_nc,:,FFgrid_nc_index,allsource_index)*weighting_nc(ii,jj)
                meteo_subgrid(i,j,:,kz_subgrid_index)=meteo_subgrid(i,j,:,kz_subgrid_index)+var4d_nc(ii,jj,surface_level_nc,:,kz_nc_index,allsource_index)*weighting_nc(ii,jj)
                meteo_subgrid(i,j,:,hmix_subgrid_index)=meteo_subgrid(i,j,:,hmix_subgrid_index)+var3d_nc(ii,jj,:,hmix_nc_index,allsource_index)*weighting_nc(ii,jj)
                meteo_subgrid(i,j,:,FF10_subgrid_index)=meteo_subgrid(i,j,:,FF10_subgrid_index)+var3d_nc(ii,jj,:,FF10_nc_index,allsource_index)*weighting_nc(ii,jj)
                meteo_subgrid(i,j,:,logz0_subgrid_index)=meteo_subgrid(i,j,:,logz0_subgrid_index)+var3d_nc(ii,jj,:,logz0_nc_index,allsource_index)*weighting_nc(ii,jj)
                meteo_subgrid(i,j,:,invL_subgrid_index)=meteo_subgrid(i,j,:,invL_subgrid_index)+var3d_nc(ii,jj,:,invL_nc_index,allsource_index)*weighting_nc(ii,jj)
                meteo_subgrid(i,j,:,inv_FFgrid_subgrid_index)=meteo_subgrid(i,j,:,inv_FFgrid_subgrid_index)+var4d_nc(ii,jj,surface_level_nc,:,inv_FFgrid_nc_index,allsource_index)*weighting_nc(ii,jj)
                meteo_subgrid(i,j,:,inv_FF10_subgrid_index)=meteo_subgrid(i,j,:,inv_FF10_subgrid_index)+var3d_nc(ii,jj,:,inv_FF10_nc_index,allsource_index)*weighting_nc(ii,jj)
                meteo_subgrid(i,j,:,ustar_subgrid_index)=meteo_subgrid(i,j,:,ustar_subgrid_index)+var3d_nc(ii,jj,:,ustar_nc_index,allsource_index)*weighting_nc(ii,jj)
                meteo_subgrid(i,j,:,J_subgrid_index)=meteo_subgrid(i,j,:,J_subgrid_index)+var4d_nc(ii,jj,surface_level_nc,:,J_nc_index,allsource_index)*weighting_nc(ii,jj)

            enddo
            enddo
            
        enddo
        enddo
        
    endif

    !Rotate wind fields if necessary and define cos and sin of the wind direction
        do j=1,integral_subgrid_dim(y_dim_index)
        do i=1,integral_subgrid_dim(x_dim_index)
        
            !Assumes it is never on the edge of the EMEP grid, not limitted
            i_nc=crossreference_integral_to_emep_subgrid(i,j,x_dim_index)
            j_nc=crossreference_integral_to_emep_subgrid(i,j,y_dim_index)
                        
            !Adjust wind direction to utm projection.
            !First determine rotation from LCC to latlon
            if (EMEP_projection_type.eq.LCC_projection_index) then
                dlatx=var2d_nc(i_nc+1,j_nc,lat_nc_index)-var2d_nc(i_nc-1,j_nc,lat_nc_index)
                dlaty=var2d_nc(i_nc,j_nc+1,lat_nc_index)-var2d_nc(i_nc,j_nc-1,lat_nc_index)
                angle_lcc=atan(dlatx/dlaty)              
            else
                angle_lcc=0.
            endif
            
            if (projection_type.eq.UTM_projection_index) then
                angle_utm2 = atan(tan((lon_integral_subgrid(i,j)-utm_lon0)/180.*pi)*sin(lat_integral_subgrid(i,j)/180.*pi)) !lon is longitude relative to middle of utm zone: lon=gl-Lambda0  
            else
                angle_utm2 = 0.
            endif
            
            angle_utm=angle_utm2+angle_lcc
                
            !write(*,*) i,j,angle_lcc*180./pi,angle_utm2*180./pi,angle_utm*180./pi
            
            u_utm = meteo_subgrid(i,j,:,ugrid_subgrid_index)*cos(angle_utm)+meteo_subgrid(i,j,:,vgrid_subgrid_index)*sin(angle_utm)                                       
            v_utm =-meteo_subgrid(i,j,:,ugrid_subgrid_index)*sin(angle_utm)+meteo_subgrid(i,j,:,vgrid_subgrid_index)*cos(angle_utm)                                       
            meteo_subgrid(i,j,:,ugrid_subgrid_index)=u_utm
            meteo_subgrid(i,j,:,vgrid_subgrid_index)=v_utm
            
            !Create cos and sin's of the lowest level wind direction  for efficient use in the dispersion equations
            ff=sqrt(meteo_subgrid(i,j,:,vgrid_subgrid_index)*meteo_subgrid(i,j,:,vgrid_subgrid_index)+meteo_subgrid(i,j,:,ugrid_subgrid_index)*+meteo_subgrid(i,j,:,ugrid_subgrid_index))
            meteo_subgrid(i,j,:,sin_subgrid_index)=meteo_subgrid(i,j,:,vgrid_subgrid_index)/ff
            meteo_subgrid(i,j,:,cos_subgrid_index)=meteo_subgrid(i,j,:,ugrid_subgrid_index)/ff
            where (ff.eq.0.)
                meteo_subgrid(i,j,:,sin_subgrid_index)=0.
                meteo_subgrid(i,j,:,cos_subgrid_index)=0.
            endwhere             
            
        enddo
        enddo

    if (use_single_time_loop_flag) then
        if (t_loop.eq.start_time_loop_index) then
            last_meteo_subgrid(:,:,:)=meteo_subgrid(:,:,1,:)        
        endif
    endif

    !Save files in ascii format
    if (save_intermediate_files) then
    if (.not.read_existing_grid_data(subgrid_meteo_file_index)) then
            i_file=subgrid_ugrid_file_index;temp_name=trim(pathname_grid(i_file))//trim(filename_grid(i_file))//'_'//trim(file_tag)//'.asc'
            write(unit_logfile,'(a)')'Writing to: '//trim(temp_name)
            call write_esri_ascii_3d_file(unit_logfile,temp_name,integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),integral_subgrid_dim(t_dim_index),integral_subgrid_delta(1),meteo_subgrid(:,:,:,ugrid_subgrid_index),x_integral_subgrid,y_integral_subgrid)
            
            i_file=subgrid_vgrid_file_index;temp_name=trim(pathname_grid(i_file))//trim(filename_grid(i_file))//'_'//trim(file_tag)//'.asc'
            write(unit_logfile,'(a)')'Writing to: '//trim(temp_name)
            call write_esri_ascii_3d_file(unit_logfile,temp_name,integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),integral_subgrid_dim(t_dim_index),integral_subgrid_delta(1),meteo_subgrid(:,:,:,vgrid_subgrid_index),x_integral_subgrid,y_integral_subgrid)
            
            i_file=subgrid_hmix_file_index;temp_name=trim(pathname_grid(i_file))//trim(filename_grid(i_file))//'_'//trim(file_tag)//'.asc'
            write(unit_logfile,'(a)')'Writing to: '//trim(temp_name)
            call write_esri_ascii_3d_file(unit_logfile,temp_name,integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),integral_subgrid_dim(t_dim_index),integral_subgrid_delta(1),meteo_subgrid(:,:,:,hmix_subgrid_index),x_integral_subgrid,y_integral_subgrid)
            
            i_file=subgrid_kz_file_index;temp_name=trim(pathname_grid(i_file))//trim(filename_grid(i_file))//'_'//trim(file_tag)//'.asc'
            write(unit_logfile,'(a)')'Writing to: '//trim(temp_name)
            call write_esri_ascii_3d_file(unit_logfile,temp_name,integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),integral_subgrid_dim(t_dim_index),integral_subgrid_delta(1),meteo_subgrid(:,:,:,kz_subgrid_index),x_integral_subgrid,y_integral_subgrid)

            i_file=subgrid_FFgrid_file_index;temp_name=trim(pathname_grid(i_file))//trim(filename_grid(i_file))//'_'//trim(file_tag)//'.asc'
            write(unit_logfile,'(a)')'Writing to: '//trim(temp_name)
            call write_esri_ascii_3d_file(unit_logfile,temp_name,integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),integral_subgrid_dim(t_dim_index),integral_subgrid_delta(1),meteo_subgrid(:,:,:,FFgrid_subgrid_index),x_integral_subgrid,y_integral_subgrid)

            i_file=subgrid_FF10_file_index;temp_name=trim(pathname_grid(i_file))//trim(filename_grid(i_file))//'_'//trim(file_tag)//'.asc'
            write(unit_logfile,'(a)')'Writing to: '//trim(temp_name)
            call write_esri_ascii_3d_file(unit_logfile,temp_name,integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),integral_subgrid_dim(t_dim_index),integral_subgrid_delta(1),meteo_subgrid(:,:,:,FF10_subgrid_index),x_integral_subgrid,y_integral_subgrid)

            i_file=subgrid_J_file_index;temp_name=trim(pathname_grid(i_file))//trim(filename_grid(i_file))//'_'//trim(file_tag)//'.asc'
            write(unit_logfile,'(a)')'Writing to: '//trim(temp_name)
            call write_esri_ascii_3d_file(unit_logfile,temp_name,integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),integral_subgrid_dim(t_dim_index),integral_subgrid_delta(1),meteo_subgrid(:,:,:,J_subgrid_index),x_integral_subgrid,y_integral_subgrid)
        !enddo
    endif
    endif
    
    if (allocated(u_utm)) deallocate(u_utm)
    if (allocated(v_utm)) deallocate(v_utm)
    if (allocated(th)) deallocate(th)
    if (allocated(ff)) deallocate(ff)
    if (allocated(weighting_nc)) deallocate(weighting_nc)

    end subroutine uEMEP_subgrid_meteo_EMEP
    