!uEMEP_read_deposition_EMEP.f90
    
    subroutine uEMEP_subgrid_deposition_EMEP
    
    use uEMEP_definitions

    implicit none

    integer i_source,i_pollutant
    integer i,j,i_cross_emep,j_cross_emep
    
    real xpos_min,xpos_max,ypos_min,ypos_max
    real xpos_area_min,xpos_area_max,ypos_area_min,ypos_area_max
    real xpos_limit,ypos_limit
    integer i_nc,j_nc
    integer ii,jj
    real weighting_nc(-1:1,-1:1)
    integer i_meteo,j_meteo
    real h_mix_loc(subgrid_dim(t_dim_index))
    real ratio(subgrid_dim(t_dim_index),n_source_index,n_pollutant_loop)
    real delta_area(2)
    real temp(subgrid_dim(t_dim_index),2)
    integer tt
    
    integer i_depo,i_loop
    
    !Functions
    real area_weighted_extended_interpolation_function

    write(unit_logfile,'(A)')'Distributing deposition to the subgrid'
    
        
        !Place the EMEP deposition velocities in the deposition subgrid, no interpolation. Now done using area interpolation
        do i_pollutant=1,n_emep_pollutant_loop
        
        !do j=1,deposition_subgrid_dim(y_dim_index)
        !do i=1,deposition_subgrid_dim(x_dim_index)
            
            !i_cross_emep=crossreference_deposition_to_emep_subgrid(i,j,x_dim_index)
            !j_cross_emep=crossreference_deposition_to_emep_subgrid(i,j,y_dim_index)
            !deposition_subgrid(i,j,:,vd_index,i_pollutant)=depo_var3d_nc(i_cross_emep,j_cross_emep,:,grid_index,i_pollutant)
            !write(*,*) i,j,deposition_subgrid(i,j,:,vd_index,i_pollutant)

        !enddo
        !enddo
        
        !Place the EMEP depositions in the concentration subgrid, no interpolation
        do j=1,subgrid_dim(y_dim_index)
        do i=1,subgrid_dim(x_dim_index)
            
            ii=crossreference_target_to_emep_subgrid(i,j,x_dim_index)
            jj=crossreference_target_to_emep_subgrid(i,j,y_dim_index)
         

            !Nearest neighbour interpolate the EMEP compounds to subgrid
            !do i_loop=1,n_pollutant_compound_loop(i_pollutant)
            !do i_depo=1,n_deposition_index
                !orig_EMEP_deposition_subgrid(i,j,:,i_depo,i_pollutant)=depo_var3d_nc(ii,jj,:,i_depo,pollutant_compound_loop_index(i_pollutant,i_loop))
                orig_EMEP_deposition_subgrid(i,j,:,drydepo_index,i_pollutant)=var3d_nc(ii,jj,:,drydepo_nc_index,allsource_index,i_pollutant)
                orig_EMEP_deposition_subgrid(i,j,:,wetdepo_index,i_pollutant)=var3d_nc(ii,jj,:,wetdepo_nc_index,allsource_index,i_pollutant)
                
            !enddo
            !enddo

            !This is overwritten in the weighted interpolation
            !write(*,*)sum(subgrid(i,j,:,emep_nonlocal_subgrid_index,:,:)),sum(subgrid(i,j,:,emep_subgrid_index,:,:))
            subgrid(i,j,:,drydepo_nonlocal_subgrid_index,:,:)=var3d_nc(ii,jj,:,drydepo_nc_index,:,:) &
                            *subgrid(i,j,:,emep_nonlocal_subgrid_index,:,:)/subgrid(i,j,:,emep_subgrid_index,:,:)
            subgrid(i,j,:,wetdepo_nonlocal_subgrid_index,:,:)=var3d_nc(ii,jj,:,wetdepo_nc_index,:,:) &
                            *subgrid(i,j,:,emep_nonlocal_subgrid_index,:,:)/subgrid(i,j,:,emep_subgrid_index,:,:)


        enddo
        enddo
        enddo

        !Place the EMEP nonlocal deposition velocities in the subgrid, area weighted interpolation
        subgrid(:,:,:,drydepo_nonlocal_subgrid_index,:,:)=0.
        subgrid(:,:,:,wetdepo_nonlocal_subgrid_index,:,:)=0.
        
        xpos_limit=dgrid_nc(lon_nc_index)/2.
        ypos_limit=dgrid_nc(lat_nc_index)/2.

        do j=1,subgrid_dim(y_dim_index)
        do i=1,subgrid_dim(x_dim_index)

            !Assumes it is never on the edge of the EMEP grid, not limitted
            i_nc=crossreference_target_to_emep_subgrid(i,j,x_dim_index)
            j_nc=crossreference_target_to_emep_subgrid(i,j,y_dim_index)

            if (EMEP_projection_type.eq.LL_projection_index) then
                delta_area(1)=xpos_limit*EMEP_grid_interpolation_size*110570.*cos(3.14159/180.*var1d_nc(j_nc,lat_nc_index))
                delta_area(2)=ypos_limit*EMEP_grid_interpolation_size*110570.                
            elseif (EMEP_projection_type.eq.LCC_projection_index) then
                delta_area(1)=xpos_limit*EMEP_grid_interpolation_size
                delta_area(2)=ypos_limit*EMEP_grid_interpolation_size
            endif
            
            !write(*,*) delta_area
            !The adjustment for the vertical profile is done here for wet deposition
            !The integral subgrid values are placed onto the target subgrid using the area weighted average over the lcoal region (xpos_limit*EMEP_grid_interpolation_size)
            if (adjust_wetdepo_integral_to_lowest_layer_flag) then
                i_meteo=crossreference_target_to_integral_subgrid(i,j,x_dim_index)
                j_meteo=crossreference_target_to_integral_subgrid(i,j,y_dim_index)
                h_mix_loc(:)=meteo_subgrid(i_meteo,j_meteo,:,hmix_subgrid_index)
                
                i_source=allsource_index
                do i_pollutant=1,n_pollutant_loop
                do tt=1,subgrid_dim(t_dim_index)
                    temp(tt,2)=area_weighted_extended_interpolation_function(x_integral_subgrid,y_integral_subgrid,integral_subgrid(:,:,tt,hmix_integral_subgrid_index,i_source,i_pollutant) &
                    ,integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),integral_subgrid_delta(x_dim_index),x_subgrid(i,j),y_subgrid(i,j),delta_area)
                    temp(tt,1)=area_weighted_extended_interpolation_function(x_integral_subgrid,y_integral_subgrid,integral_subgrid(:,:,tt,hsurf_integral_subgrid_index,i_source,i_pollutant) &
                    ,integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),integral_subgrid_delta(x_dim_index),x_subgrid(i,j),y_subgrid(i,j),delta_area)
                    
                    ratio(tt,i_source,i_pollutant)=temp(tt,2)/temp(tt,1)*H_emep/h_mix_loc(tt)
                    if (temp(tt,1).eq.0) ratio(tt,i_source,i_pollutant)=0.
                    ratio(tt,i_source,i_pollutant)=max(0.,ratio(tt,i_source,i_pollutant))
                    ratio(tt,i_source,i_pollutant)=min(1.,ratio(tt,i_source,i_pollutant))
                    !write(*,'(2i,3es12.2)') i,j,temp(tt,1),temp(tt,2),ratio(tt,i_source,i_pollutant)
                enddo
                enddo
                
            else
                ratio=1.
            endif
            

            xpos_area_max=xproj_subgrid(i,j)+xpos_limit
            xpos_area_min=xproj_subgrid(i,j)-xpos_limit
            ypos_area_max=yproj_subgrid(i,j)+ypos_limit
            ypos_area_min=yproj_subgrid(i,j)-ypos_limit
                  
            do jj=j_nc-1,j_nc+1
            do ii=i_nc-1,i_nc+1             
            
                xpos_min=max(xpos_area_min,var1d_nc(ii,lon_nc_index)-dgrid_nc(lon_nc_index)/2.)
                xpos_max=min(xpos_area_max,var1d_nc(ii,lon_nc_index)+dgrid_nc(lon_nc_index)/2.)
                ypos_min=max(ypos_area_min,var1d_nc(jj,lat_nc_index)-dgrid_nc(lat_nc_index)/2.)
                ypos_max=min(ypos_area_max,var1d_nc(jj,lat_nc_index)+dgrid_nc(lat_nc_index)/2.)
               
                !Determine the area intersection of the EMEP grid and an EMEP grid size centred on the integral subgrid
                if (xpos_max.gt.xpos_min.and.ypos_max.gt.ypos_min) then
                    weighting_nc(ii-i_nc,jj-j_nc)=(ypos_max-ypos_min)*(xpos_max-xpos_min)/dgrid_nc(lon_nc_index)/dgrid_nc(lat_nc_index)
                else
                    weighting_nc(ii-i_nc,jj-j_nc)=0.
                endif                

                i_source=allsource_index
                do i_pollutant=1,n_emep_pollutant_loop

                    subgrid(i,j,:,drydepo_nonlocal_subgrid_index,allsource_index,i_pollutant)=subgrid(i,j,:,drydepo_nonlocal_subgrid_index,allsource_index,i_pollutant) &
                            +var3d_nc(ii,jj,:,drydepo_nc_index,i_source,i_pollutant) &
                            *subgrid(i,j,:,emep_nonlocal_subgrid_index,allsource_index,i_pollutant)/subgrid(i,j,:,emep_subgrid_index,allsource_index,i_pollutant) &
                            *weighting_nc(ii-i_nc,jj-j_nc)
                    subgrid(i,j,:,wetdepo_nonlocal_subgrid_index,allsource_index,i_pollutant)=subgrid(i,j,:,wetdepo_nonlocal_subgrid_index,allsource_index,i_pollutant) &
                            +var3d_nc(ii,jj,:,wetdepo_nc_index,i_source,i_pollutant) &
                            *(1.-subgrid(i,j,:,emep_local_subgrid_index,allsource_index,i_pollutant)/subgrid(i,j,:,emep_subgrid_index,i_source,i_pollutant)*ratio(:,i_source,i_pollutant)) &
                            *weighting_nc(ii-i_nc,jj-j_nc)
                
                enddo
            
            enddo
            enddo
                    !write(*,*) i,j,ratio,subgrid(i,j,:,emep_local_subgrid_index,i_source,:)/subgrid(i,j,:,emep_subgrid_index,i_source,:)
            
        
        
            !write(*,*) i,j,deposition_subgrid(i,j,:,vd_index,:),depo_var3d_nc(i_nc,j_nc,:,grid_index,:)
        enddo
        enddo

    end subroutine uEMEP_subgrid_deposition_EMEP
    
    subroutine uEMEP_calculate_deposition
    
    use uEMEP_definitions
    
    implicit none

    character(256) temp_name
    integer subsource_index,source_index    
    integer i_pollutant,i_loop
    integer i,j
    !real sum_temp(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index))
    integer i_cross_deposition,j_cross_deposition
    
    write(unit_logfile,'(A)')'Combining local deposition sources'

    !Calculate redistributed subgrid allsource deposition
    !
    subgrid(:,:,:,drydepo_local_subgrid_index,allsource_index,:)=0.
    !subgrid(:,:,:,drydepo_nonlocal_subgrid_index,allsource_index,:)=0.
    subgrid(:,:,:,wetdepo_local_subgrid_index,allsource_index,:)=0.
    !subgrid(:,:,:,wetdepo_nonlocal_subgrid_index,allsource_index,:)=0.
    do source_index=1,n_source_index
    if (calculate_source(source_index)) then
        do j=1,subgrid_dim(y_dim_index)
        do i=1,subgrid_dim(x_dim_index)
            
            !i_cross_deposition=crossreference_target_to_deposition_subgrid(i,j,x_dim_index)
            !j_cross_deposition=crossreference_target_to_deposition_subgrid(i,j,y_dim_index)
            !subgrid(i,j,:,drydepo_nonlocal_subgrid_index,source_index,:)=subgrid(i,j,:,emep_nonlocal_subgrid_index,source_index,:) &
            !     *deposition_subgrid(i_cross_deposition,j_cross_deposition,:,vd_index,:)
            
            !subgrid(i,j,:,drydepo_nonlocal_subgrid_index,allsource_index,:)=subgrid(i,j,:,drydepo_nonlocal_subgrid_index,allsource_index,:)+subgrid(i,j,:,drydepo_nonlocal_subgrid_index,source_index,:)
            subgrid(i,j,:,drydepo_local_subgrid_index,allsource_index,:)=subgrid(i,j,:,drydepo_local_subgrid_index,allsource_index,:)+subgrid(i,j,:,drydepo_local_subgrid_index,source_index,:)
            !write(*,*) subgrid(i,j,:,emep_nonlocal_subgrid_index,source_index,1),deposition_subgrid(i_cross_deposition,j_cross_deposition,:,vd_index,1)
            !subgrid(i,j,:,wetdepo_nonlocal_subgrid_index,allsource_index,:)=subgrid(i,j,:,wetdepo_nonlocal_subgrid_index,allsource_index,:)+subgrid(i,j,:,wetdepo_nonlocal_subgrid_index,source_index,:)
            subgrid(i,j,:,wetdepo_local_subgrid_index,allsource_index,:)=subgrid(i,j,:,wetdepo_local_subgrid_index,allsource_index,:)+subgrid(i,j,:,wetdepo_local_subgrid_index,source_index,:)
            
            
        enddo
        enddo
    endif
    enddo
    
    !Convert the nonlocal depositions from mg/m2/hr to ug/m2/sec for compatability with local calculations
    subgrid(:,:,:,drydepo_nonlocal_subgrid_index,allsource_index,:)=subgrid(:,:,:,drydepo_nonlocal_subgrid_index,allsource_index,:)*1000./3600.
    subgrid(:,:,:,wetdepo_nonlocal_subgrid_index,allsource_index,:)=subgrid(:,:,:,wetdepo_nonlocal_subgrid_index,allsource_index,:)*1000./3600.

    
    
    end subroutine uEMEP_calculate_deposition
    

    
    subroutine uEMEP_set_deposition_velocities
    
    use uEMEP_definitions

    implicit none

    integer i_source,i_pollutant
    integer i,j,i_cross_emep,j_cross_emep
    
    real xpos_min,xpos_max,ypos_min,ypos_max
    real xpos_area_min,xpos_area_max,ypos_area_min,ypos_area_max
    real xpos_limit,ypos_limit
    integer i_nc,j_nc
    integer ii,jj
    real weighting_nc(-1:1,-1:1)
    
    integer i_landuse
    integer i_depo,i_loop

    write(unit_logfile,'(A)')'Placing EMEP deposition velocities on the deposition subgrid'
    
    !Fill in with a default value
    deposition_subgrid(:,:,:,:,1:n_pollutant_loop)=0.
    !write(*,*) size(deposition_subgrid)

    drydepo_vd_default(1:n_compound_index)=0.0 !Positive value in m/s
    drydepo_vd_default(nh3_index)=0.03 !Positive value in m/s
    wetdepo_scavanging_rate(1:n_compound_index)=0.0 !Positive value in m/s
    wetdepo_scavanging_rate(nh3_index)=0.5*1.e6/1000. !Dimensions of m-1 W/h
    
        do i_pollutant=1,n_pollutant_loop
            deposition_subgrid(:,:,:,vd_index,i_pollutant)=drydepo_vd_default(pollutant_loop_index(i_pollutant))
            !write(*,*) i_source,i_pollutant,sum(deposition_subgrid(:,:,:,vd_index,i_pollutant)),drydepo_vd_default(i_pollutant)
        enddo
    !stop
   
        !Place the EMEP deposition velocities in the deposition subgrid, area weighted interpolation
        deposition_subgrid(:,:,:,vd_index,:)=0.
        
        xpos_limit=dgrid_nc(lon_nc_index)/2.
        ypos_limit=dgrid_nc(lat_nc_index)/2.

        do j=1,deposition_subgrid_dim(y_dim_index)
        do i=1,deposition_subgrid_dim(x_dim_index)

            !Cross reference EMEP grid with limits
            i_nc=crossreference_deposition_to_emep_subgrid(i,j,x_dim_index)
            j_nc=crossreference_deposition_to_emep_subgrid(i,j,y_dim_index)
                        
            i_nc=min(max(2,i_nc),dim_length_nc(x_dim_nc_index)-1)
            j_nc=min(max(2,j_nc),dim_length_nc(y_dim_nc_index)-1)
            
            !write(*,*) i,j,i_nc,j_nc
            xpos_area_max=xproj_deposition_subgrid(i,j)+xpos_limit
            xpos_area_min=xproj_deposition_subgrid(i,j)-xpos_limit
            ypos_area_max=yproj_deposition_subgrid(i,j)+ypos_limit
            ypos_area_min=yproj_deposition_subgrid(i,j)-ypos_limit
                  
            do jj=j_nc-1,j_nc+1
            do ii=i_nc-1,i_nc+1             
            
                xpos_min=max(xpos_area_min,var1d_nc(ii,lon_nc_index)-dgrid_nc(lon_nc_index)/2.)
                xpos_max=min(xpos_area_max,var1d_nc(ii,lon_nc_index)+dgrid_nc(lon_nc_index)/2.)
                ypos_min=max(ypos_area_min,var1d_nc(jj,lat_nc_index)-dgrid_nc(lat_nc_index)/2.)
                ypos_max=min(ypos_area_max,var1d_nc(jj,lat_nc_index)+dgrid_nc(lat_nc_index)/2.)
               
                !Determine the area intersection of the EMEP grid and an EMEP grid size centred on the deposition subgrid
                if (xpos_max.gt.xpos_min.and.ypos_max.gt.ypos_min) then
                    weighting_nc(ii-i_nc,jj-j_nc)=(ypos_max-ypos_min)*(xpos_max-xpos_min)/dgrid_nc(lon_nc_index)/dgrid_nc(lat_nc_index)
                else
                    weighting_nc(ii-i_nc,jj-j_nc)=0.
                endif                

                
                if (read_landuse_flag) then
                    do i_landuse=1,n_landuse_index-1
                        deposition_subgrid(i,j,:,vd_index,:)=deposition_subgrid(i,j,:,vd_index,:)+depo_var3d_nc(ii,jj,:,i_landuse,:)*landuse_subgrid(i,j,i_landuse)*weighting_nc(ii-i_nc,jj-j_nc)
                    enddo
                else
                    deposition_subgrid(i,j,:,vd_index,:)=deposition_subgrid(i,j,:,vd_index,:)+depo_var3d_nc(ii,jj,:,grid_index,:)*weighting_nc(ii-i_nc,jj-j_nc)                    
                endif
                

            enddo
            enddo
            
            !Fill in zeros with the nearest EMEP value. 0 can appear because of non overlapping EMEP and landuse grids, I think
            where (deposition_subgrid(i,j,:,vd_index,:).eq.0) deposition_subgrid(i,j,:,vd_index,:)=depo_var3d_nc(i_nc,j_nc,:,grid_index,:)
            !write(*,'(2i,23es12.2,f)') i,j,deposition_subgrid(i,j,:,vd_index,:),depo_var3d_nc(i_nc,j_nc,:,grid_index,:),sum(depo_var3d_nc(ii,jj,:,1:n_landuse_index-1,:)),sum(landuse_subgrid(i,j,1:n_landuse_index-1))
        enddo
        enddo

    end subroutine uEMEP_set_deposition_velocities