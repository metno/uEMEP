module redistribute_data

    use save_netcdf_file, only: mean_mask

    implicit none
    private

    public :: uEMEP_redistribute_local_source, uEMEP_disperse_local_source, &
        uEMEP_combine_local_source

contains

!uEMEP_redistribute_local_source
!Same routine for all sources
    
    subroutine uEMEP_redistribute_local_source(source_index)
    
    use uEMEP_definitions
    
    implicit none

    integer i,j
    integer source_index
    integer ii,jj,tt
    integer integral_counter
    real sum_integral(n_pollutant_loop)
    integer i_start,i_end,j_start,j_end,t_start,t_end
    integer emep_subsource
    integer i_cross_integral,j_cross_integral
    real xpos_limit,ypos_limit
    real xpos_area_min,xpos_area_max,ypos_area_min,ypos_area_max
    real xpos_integral_subgrid,ypos_integral_subgrid
    integer i_pollutant

    !allocate (sum_integral(subgrid_dim(1),subgrid_dim(2))) !Can just be a scalar
    !allocate (scaling_factor_traffic_subgrid(subgrid_dim(1),subgrid_dim(2))) !Can just be a scalar
    !allocate (traffic_redistributed_local_subgrid(subgrid_dim(1),subgrid_dim(2))) !Can just be a scalar
    
    if (local_subgrid_method_flag.ne.1) return
    
 
    write(unit_logfile,'(A)') ''
    write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Redistribute local source using EMEP concentrations (uEMEP_redistribute_local_source)'
	write(unit_logfile,'(A)') '================================================================'

    !No subsources for the emep related arrays
    emep_subsource=1


    write(unit_logfile,'(2A)')'Calculating the scaling factor and local contribution at each subgrid for ',trim(source_file_str(source_index))

    xpos_limit=dgrid_nc(lon_nc_index)/2.*EMEP_grid_interpolation_size*local_fraction_grid_size_scaling
    ypos_limit=dgrid_nc(lat_nc_index)/2.*EMEP_grid_interpolation_size*local_fraction_grid_size_scaling

    !Set the start and end times of the loop
    t_start=1
    t_end=subgrid_dim(t_dim_index)

    do tt=t_start,t_end
        
    !Calculate the mean concentration of the integral values for each subgrid point
    do j=1,subgrid_dim(y_dim_index)
    do i=1,subgrid_dim(x_dim_index)
        
        if (use_subgrid(i,j,source_index)) then

            sum_integral=0.
            integral_counter=0
        
            i_cross_integral=crossreference_target_to_integral_subgrid(i,j,x_dim_index)
            j_cross_integral=crossreference_target_to_integral_subgrid(i,j,y_dim_index)

            xpos_integral_subgrid=xproj_subgrid(i,j)
            ypos_integral_subgrid=yproj_subgrid(i,j)


            !Use the wind direction to move the target area downwind
            if (use_downwind_position_flag.and.hourly_calculations) then
    
                !Use the same area as for non upwind
                i_start=max(1,i_cross_integral-integral_subgrid_loop_index(x_dim_index))
                i_end=min(integral_subgrid_dim(x_dim_index),i_cross_integral+integral_subgrid_loop_index(x_dim_index))
                j_start=max(1,j_cross_integral-integral_subgrid_loop_index(y_dim_index))
                j_end=min(integral_subgrid_dim(y_dim_index),j_cross_integral+integral_subgrid_loop_index(y_dim_index))

                xpos_area_max=xpos_integral_subgrid+xpos_limit
                xpos_area_min=xpos_integral_subgrid-xpos_limit
                ypos_area_max=ypos_integral_subgrid+ypos_limit
                ypos_area_min=ypos_integral_subgrid-ypos_limit

            else
       
                i_start=max(1,i_cross_integral-integral_subgrid_loop_index(x_dim_index))
                i_end=min(integral_subgrid_dim(x_dim_index),i_cross_integral+integral_subgrid_loop_index(x_dim_index))
                j_start=max(1,j_cross_integral-integral_subgrid_loop_index(y_dim_index))
                j_end=min(integral_subgrid_dim(y_dim_index),j_cross_integral+integral_subgrid_loop_index(y_dim_index))
 
                xpos_area_max=xpos_integral_subgrid+xpos_limit
                xpos_area_min=xpos_integral_subgrid-xpos_limit
                ypos_area_max=ypos_integral_subgrid+ypos_limit
                ypos_area_min=ypos_integral_subgrid-ypos_limit
                
            endif

            !write(*,*) i_start-i_cross_integral,i_end-i_cross_integral,j_start-j_cross_integral,j_end-j_cross_integral
                !Calculate the average grid concentration at each integral subgrid based on proxy integral
                do jj=j_start,j_end
                do ii=i_start,i_end
                    
                    xpos_integral_subgrid=xproj_integral_subgrid(ii,jj)
                    ypos_integral_subgrid=yproj_integral_subgrid(ii,jj)
                    
                    if (xpos_integral_subgrid.ge.xpos_area_min.and.xpos_integral_subgrid.le.xpos_area_max &
                        .and.ypos_integral_subgrid.ge.ypos_area_min.and.ypos_integral_subgrid.le.ypos_area_max) then                  

                        !do i_pollutant=1,n_pollutant_loop
                            sum_integral(:)=sum_integral(:)+integral_subgrid(ii,jj,tt,hsurf_average_subgrid_index,source_index,:)
                            !write(*,*) ii,jj,integral_subgrid(ii,jj,tt,hsurf_average_subgrid_index,source_index,1)
                        !enddo
                        
                        integral_counter=integral_counter+1
                        !write(*,*) i,j,ii-i_cross_integral,jj-j_cross_integral
                    endif  
                enddo
                enddo

            !Calculate scaling factor
            do i_pollutant=1,n_pollutant_loop
                if (sum_integral(i_pollutant).ne.0) then
                    subgrid(i,j,tt,scaling_factor_subgrid_index,source_index,i_pollutant)=subgrid(i,j,tt,proxy_subgrid_index,source_index,i_pollutant)/sum_integral(i_pollutant)*integral_counter
                else
                    subgrid(i,j,tt,scaling_factor_subgrid_index,source_index,i_pollutant)=0.
                endif
                !write(*,*) i_pollutant,subgrid(i,j,tt,proxy_subgrid_index,source_index,i_pollutant),sum_integral(i_pollutant)
                if (isnan(subgrid(i,j,tt,scaling_factor_subgrid_index,source_index,i_pollutant))) write(*,*) i,j,sum_integral(i_pollutant),integral_counter
                if (isnan(subgrid(i,j,tt,emep_local_subgrid_index,source_index,i_pollutant))) write(*,*) 'L',i,j,sum_integral(i_pollutant),integral_counter
            enddo
        
        endif
        
    enddo
        !write(*,*) 'Redistribution',j,' of ',subgrid_dim(2)
    enddo
        write(*,*) 'Redistribution time ',tt,' of ',subgrid_dim(t_dim_index)
    enddo

    do i_pollutant=1,n_pollutant_loop
    write(unit_logfile,'(i,A,f12.4)') i_pollutant,' Max scaling factor for pollutant  =',maxval(subgrid(:,:,:,scaling_factor_subgrid_index,source_index,i_pollutant))
    write(unit_logfile,'(i,A,f12.4)') i_pollutant,' Mean scaling factor for pollutant =',sum(subgrid(:,:,:,scaling_factor_subgrid_index,source_index,i_pollutant)) &
        /subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index)
    enddo
    
    !Calculate redistributed subgrid source concentrations
    subgrid(:,:,:,local_subgrid_index,source_index,:)=subgrid(:,:,:,scaling_factor_subgrid_index,source_index,:)*subgrid(:,:,:,emep_local_subgrid_index,source_index,:)
   
    
    end subroutine uEMEP_redistribute_local_source
    

!uEMEP_disperse_local_source
!Same routine for all sources
    !This routine uses the emission factors and meteorology to calculate the local concentrations
    !This is instead of redistributing using proxy data
    !Gives the same output as the uEMEP_redistribute_local_source
    
    subroutine uEMEP_disperse_local_source(source_index)
    
    use uEMEP_definitions
    
    implicit none

    integer source_index
    
    !allocate (sum_integral(subgrid_dim(1),subgrid_dim(2))) !Can just be a scalar
    !allocate (scaling_factor_traffic_subgrid(subgrid_dim(1),subgrid_dim(2))) !Can just be a scalar
    !allocate (traffic_redistributed_local_subgrid(subgrid_dim(1),subgrid_dim(2))) !Can just be a scalar

    if (local_subgrid_method_flag.ne.2.and.local_subgrid_method_flag.ne.3.and.local_subgrid_method_flag.ne.4) return
    
    write(unit_logfile,'(A)') ''
    write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Redistribute local source using dispersion (uEMEP_disperse_local_source)'
	write(unit_logfile,'(A)') '================================================================'

 
    !do subsource_index=1,n_subsource(source_index)

        !Calculate redistributed subgrid source concentrations
        subgrid(:,:,:,scaling_factor_subgrid_index,source_index,:)=1.
        subgrid(:,:,:,local_subgrid_index,source_index,:)=subgrid(:,:,:,scaling_factor_subgrid_index,source_index,:)*subgrid(:,:,:,proxy_subgrid_index,source_index,:)
        if (trace_emissions_from_in_region) then
        subgrid_from_in_region(:,:,:,local_subgrid_index,source_index,:)=subgrid(:,:,:,scaling_factor_subgrid_index,source_index,:)*subgrid_from_in_region(:,:,:,proxy_subgrid_index,source_index,:)
        !subgrid_from_in_region(:,:,:,local_subgrid_index,source_index,:)=0
        endif
        
    !enddo !Subsource loop
    
    end subroutine uEMEP_disperse_local_source
    
    
    !uEMEP_combine_local_source
    !Saves the combination of local sources for final results
    
    subroutine uEMEP_combine_local_source
    
    use uEMEP_definitions
    
    implicit none

    integer source_index    
    integer i_pollutant
    integer i,j
    integer i_sp
    real sum_temp(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index))
        
    real, allocatable :: subgrid_dummy(:,:,:,:,:,:)
    integer in_region_loop, n_in_region_loop

    !if (interpolate_subgrids_flag) then
        !write(unit_logfile,'(a)') 'Interpolate routines not currently active. Doing nothing'
        !call uEMEP_interpolate_auto_subgrid
        !return
        !stop
        !call uEMEP_interpolate_subgrids
        !call uEMEP_linear_interpolate_subgrids
        !call uEMEP_bilinear_interpolate_subgrids
        
        !Remember to reset the use_subgrids val and logical so that everything will be used in the end
    !endif
    
   !Loop over the normal subgrid and the in region subgrid by saving  to a dummy variable and pputting back when finished
    if (trace_emissions_from_in_region) then
        n_in_region_loop=2
        if (.not.allocated(subgrid_dummy)) allocate (subgrid_dummy(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_subgrid_index,n_source_index,n_pollutant_loop))
     else
        n_in_region_loop=1        
    endif
    
    do in_region_loop=1,n_in_region_loop
        
        write(unit_logfile,'(a)')'--------------------------'
        if (in_region_loop.eq.1) write(unit_logfile,'(a)')'Combining the local and nonlocal contributions at each subgrid'
        if (in_region_loop.eq.2) write(unit_logfile,'(a)')'Combining the local and nonlocal contributions at each subgrid from in region'
        write(unit_logfile,'(a)')'--------------------------'
        
   !Loop over the normal subgrid and the in region subgrid by saving  to a dummy variable and pputting back when finished
    if (trace_emissions_from_in_region.and.in_region_loop.eq.2) then
        subgrid_dummy=subgrid
        subgrid=subgrid_from_in_region
    endif
 
    !Calculate redistributed subgrid allsource concentrations
    subgrid(:,:,:,local_subgrid_index,allsource_index,:)=0.
    do source_index=1,n_source_index
        if (calculate_source(source_index).and.source_index.ne.allsource_index) then
             subgrid(:,:,:,local_subgrid_index,allsource_index,:)=subgrid(:,:,:,local_subgrid_index,allsource_index,:)+subgrid(:,:,:,local_subgrid_index,source_index,:)
        endif
        !Add the selected EMEP local sources to this as well, if they are not already included in the subgrid downscaling
        if (calculate_EMEP_source(source_index).and..not.calculate_source(source_index).and.source_index.ne.allsource_index) then
             subgrid(:,:,:,local_subgrid_index,allsource_index,:)=subgrid(:,:,:,local_subgrid_index,allsource_index,:)+subgrid(:,:,:,emep_local_subgrid_index,source_index,:)
            !write(*,*) source_index,sum(subgrid(:,:,:,emep_local_subgrid_index,source_index,:)),sum(subgrid(:,:,:,local_subgrid_index,allsource_index,:)),sum(subgrid(:,:,:,emep_nonlocal_subgrid_index,allsource_index,:))
        endif
    enddo
    
    if (EMEP_additional_grid_interpolation_size.gt.0) then    
    do i_pollutant=1,n_pollutant_loop
        !If the compound is PM2.5 or PM10 then add the non PPM part to the non-local
        if (pollutant_loop_index(i_pollutant).eq.pm25_index.or.pollutant_loop_index(i_pollutant).eq.pm10_index) then
            write(unit_logfile,'(A,A)') 'Pollutant: ',trim(var_name_nc(conc_nc_index,pollutant_loop_index(i_pollutant),allsource_nc_index))
            write(unit_logfile,'(A,f12.2)') 'MEAN PPM NONLOCAL ADDITIONAL: ',sum((subgrid(:,:,:,emep_additional_nonlocal_subgrid_index,allsource_index,i_pollutant)))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index)
            subgrid(:,:,:,emep_additional_nonlocal_subgrid_index,allsource_index,i_pollutant)=subgrid(:,:,:,emep_additional_nonlocal_subgrid_index,allsource_index,i_pollutant) &
            +(comp_EMEP_subgrid(:,:,:,pollutant_loop_index(i_pollutant))-subgrid(:,:,:,emep_subgrid_index,allsource_index,i_pollutant))
            write(unit_logfile,'(A,f12.2)') 'MEAN ADD REST NONLOCAL ADDITIONAL: ',sum((comp_EMEP_subgrid(:,:,:,pollutant_loop_index(i_pollutant))-subgrid(:,:,:,emep_subgrid_index,allsource_index,i_pollutant)))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index)
        endif
    
        !subgrid(:,:,:,total_subgrid_index,allsource_index,i_pollutant)=subgrid(:,:,:,local_subgrid_index,allsource_index,i_pollutant)+subgrid(:,:,:,emep_nonlocal_subgrid_index,allsource_index,i_pollutant)
     
    enddo
    endif
    
    do i_pollutant=1,n_pollutant_loop
        !If the compound is PM2.5 or PM10 then add the non PPM part to the non-local
        if (pollutant_loop_index(i_pollutant).eq.pm25_index.or.pollutant_loop_index(i_pollutant).eq.pm10_index) then
            write(unit_logfile,'(A,A)') 'Pollutant: ',trim(var_name_nc(conc_nc_index,pollutant_loop_index(i_pollutant),allsource_nc_index))
            write(unit_logfile,'(A,f12.2)') 'MEAN PPM NONLOCAL: ',mean_mask(subgrid(:,:,:,emep_nonlocal_subgrid_index,allsource_index,i_pollutant),use_subgrid(:,:,allsource_index),subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index)) !sum((subgrid(:,:,:,emep_nonlocal_subgrid_index,allsource_index,i_pollutant)))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index)
            write(unit_logfile,'(A,f12.2)') 'MEAN PPM LOCAL: ',mean_mask(subgrid(:,:,:,emep_local_subgrid_index,allsource_index,i_pollutant),use_subgrid(:,:,allsource_index),subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index)) !sum((subgrid(:,:,:,emep_local_subgrid_index,allsource_index,i_pollutant)))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index)
            write(unit_logfile,'(A,f12.2)') 'MEAN PPM TOTAL: ',mean_mask(subgrid(:,:,:,emep_subgrid_index,allsource_index,i_pollutant),use_subgrid(:,:,allsource_index),subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index)) !sum((subgrid(:,:,:,emep_subgrid_index,allsource_index,i_pollutant)))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index)
            write(unit_logfile,'(A,f12.2)') 'MEAN COMP PM TOTAL: ',mean_mask(comp_EMEP_subgrid(:,:,:,pollutant_loop_index(i_pollutant)),use_subgrid(:,:,allsource_index),subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index)) !sum((comp_EMEP_subgrid(:,:,:,pollutant_loop_index(i_pollutant))))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index)
            write(unit_logfile,'(A,f12.2)') 'MEAN COMP PM ORIGINAL: ',mean_mask(orig_EMEP_subgrid(:,:,:,pollutant_loop_index(i_pollutant)),use_subgrid(:,:,allsource_index),subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index)) !sum((orig_EMEP_subgrid(:,:,:,pollutant_loop_index(i_pollutant))))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index)
            subgrid(:,:,:,emep_nonlocal_subgrid_index,allsource_index,i_pollutant)=subgrid(:,:,:,emep_nonlocal_subgrid_index,allsource_index,i_pollutant) &
            +(comp_EMEP_subgrid(:,:,:,pollutant_loop_index(i_pollutant))-subgrid(:,:,:,emep_subgrid_index,allsource_index,i_pollutant))
            write(unit_logfile,'(A,f12.2)') 'MEAN ADD REST NONLOCAL: ',mean_mask(comp_EMEP_subgrid(:,:,:,pollutant_loop_index(i_pollutant))-subgrid(:,:,:,emep_subgrid_index,allsource_index,i_pollutant),use_subgrid(:,:,allsource_index),subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index)) !sum((comp_EMEP_subgrid(:,:,:,pollutant_loop_index(i_pollutant))-subgrid(:,:,:,emep_subgrid_index,allsource_index,i_pollutant)))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index)
            write(unit_logfile,'(A,f12.2)') 'MEAN NEW PM NONLOCAL: ',mean_mask(subgrid(:,:,:,emep_nonlocal_subgrid_index,allsource_index,i_pollutant),use_subgrid(:,:,allsource_index),subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index)) !sum((subgrid(:,:,:,emep_nonlocal_subgrid_index,allsource_index,i_pollutant)))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index)
            if (sum(comp_EMEP_subgrid(:,:,:,pollutant_loop_index(i_pollutant))-subgrid(:,:,:,emep_subgrid_index,allsource_index,i_pollutant)).lt.0) write(unit_logfile,'(A)') 'WARNING!!!: PPM EMEP is more than total EMEP PM. Negative non PPM contributions.'
        endif
    
        subgrid(:,:,:,total_subgrid_index,allsource_index,i_pollutant)=subgrid(:,:,:,local_subgrid_index,allsource_index,i_pollutant)+subgrid(:,:,:,emep_nonlocal_subgrid_index,allsource_index,i_pollutant)
    enddo

    !Replace back
    if (trace_emissions_from_in_region.and.in_region_loop.eq.2) then
        subgrid_from_in_region=subgrid
        subgrid=subgrid_dummy
    endif

    enddo !in_region_loop
    
    if (trace_emissions_from_in_region) then
        if (allocated(subgrid_dummy)) deallocate (subgrid_dummy)
    endif
 
    
     do i_pollutant=1,n_pollutant_loop
       !Place the results in the compound results
        !do i_loop=1,n_pollutant_compound_loop(i_pollutant)
            comp_subgrid(:,:,:,pollutant_loop_index(i_pollutant))=subgrid(:,:,:,total_subgrid_index,allsource_index,i_pollutant)
            if (trace_emissions_from_in_region) then
                comp_subgrid_from_in_region(:,:,:,pollutant_loop_index(i_pollutant))=subgrid_from_in_region(:,:,:,total_subgrid_index,allsource_index,i_pollutant)
            endif
        !enddo
        !write(*,'(2i,4f16.2)') i_pollutant,pollutant_loop_index(i_pollutant)&
        !    ,sum(subgrid(:,:,:,total_subgrid_index,allsource_index,i_pollutant))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index)&
        !    ,sum(subgrid(:,:,:,local_subgrid_index,allsource_index,i_pollutant))&
        !    ,sum(subgrid(:,:,:,emep_nonlocal_subgrid_index,allsource_index,i_pollutant))&
        !    ,sum(comp_subgrid(:,:,:,pollutant_loop_index(i_pollutant)))
        
     enddo

     !Replace the species PPM with the actual species PPM used in the local fraction and then replace this with the nonlocal part.
     !Only works for the complete group of species
     if (save_emep_species) then
         !Replace the primary species value with the one used in the calculations for consistency
         write(unit_logfile,'(A)') 'ppm read from surf and read emep. Difference should just be deposition'
         write(unit_logfile,'(A,2f12.2)') 'PPM25 (init_sp,comp)',sum(species_EMEP_subgrid(:,:,:,pm25_sp_index,sp_ppm_index))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index) &
             ,sum(subgrid(:,:,:,emep_subgrid_index,allsource_index,pollutant_loop_back_index(pm25_nc_index)))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index)
         write(unit_logfile,'(A,2f12.2)') 'PPM10 (init_sp,comp)',sum(species_EMEP_subgrid(:,:,:,pm10_sp_index,sp_ppm_index))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index) &
             ,sum(subgrid(:,:,:,emep_subgrid_index,allsource_index,pollutant_loop_back_index(pm10_nc_index)))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index)
         write(unit_logfile,'(A)') 'pm surf init species and summed species should be the same. comp should be different (difference between SURF and SURF_rh50)'
         write(unit_logfile,'(A,3f12.2)') 'PM25 (init_sp,sum_sp,comp)',sum(species_EMEP_subgrid(:,:,:,pm25_sp_index,sp_pm_index))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index) &
             ,sum(species_EMEP_subgrid(:,:,:,pm25_sp_index,1:sp_ppm_index))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index) &
             ,sum(comp_EMEP_subgrid(:,:,:,pm25_nc_index))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index)
         write(unit_logfile,'(A,3f12.2)') 'PM10 (init_sp,sum_sp,comp)',sum(species_EMEP_subgrid(:,:,:,pm10_sp_index,sp_pm_index))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index) &
             ,sum(species_EMEP_subgrid(:,:,:,pm10_sp_index,1:sp_ppm_index))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index) &
             ,sum(comp_EMEP_subgrid(:,:,:,pm10_nc_index))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index)
         !The remaining species will not, most likely, be normalised with the total. This is done here but the total ppm is kept
         !Something weird is going on here. Need to check
         
         !Normalise the species surface values with the total grid values to account for any difference du to deposition and water
         !Because t can have dimension 1 and 'sum' function does not think it is a dimension anymore then need to do the sum over the time loop. This turned out to not be the case so they are the same           
         if (use_single_time_loop_flag) then
             sum_temp(:,:,:)=sum(species_EMEP_subgrid(:,:,:,pm25_sp_index,1:sp_ppm_index),4)
         else
             sum_temp(:,:,:)=sum(species_EMEP_subgrid(:,:,:,pm25_sp_index,1:sp_ppm_index),4)
         endif
         do i_sp=1,sp_ppm_index
             species_EMEP_subgrid(:,:,:,pm25_sp_index,i_sp)=species_EMEP_subgrid(:,:,:,pm25_sp_index,i_sp) &
                *(comp_EMEP_subgrid(:,:,:,pm25_nc_index)) &
                /sum_temp
            where (sum_temp.le.0) species_EMEP_subgrid(:,:,:,pm25_sp_index,i_sp)=0
            where (isnan(species_EMEP_subgrid(:,:,:,pm25_sp_index,i_sp)).or.species_EMEP_subgrid(:,:,:,pm25_sp_index,i_sp).le.0) species_EMEP_subgrid(:,:,:,pm25_sp_index,i_sp)=0
         enddo
         if (use_single_time_loop_flag) then
             sum_temp(:,:,:)=sum(species_EMEP_subgrid(:,:,:,pm10_sp_index,1:sp_ppm_index),4)
         else
             sum_temp(:,:,:)=sum(species_EMEP_subgrid(:,:,:,pm10_sp_index,1:sp_ppm_index),4)
         endif
         do i_sp=1,sp_ppm_index
            species_EMEP_subgrid(:,:,:,pm10_sp_index,i_sp)=species_EMEP_subgrid(:,:,:,pm10_sp_index,i_sp) &
                *(comp_EMEP_subgrid(:,:,:,pm10_nc_index)) &
                /sum_temp
            where (sum_temp.le.0) species_EMEP_subgrid(:,:,:,pm10_sp_index,i_sp)=0
            where (isnan(species_EMEP_subgrid(:,:,:,pm10_sp_index,i_sp)).or.species_EMEP_subgrid(:,:,:,pm10_sp_index,i_sp).le.0) species_EMEP_subgrid(:,:,:,pm10_sp_index,i_sp)=0
         enddo
         
         write(unit_logfile,'(A)') 'sum_sp should be the same as comp after scaling, init should be unchanged'
         write(unit_logfile,'(A,3f12.2)') 'PM25 (init_sp,sum_sp,comp)',sum(species_EMEP_subgrid(:,:,:,pm25_sp_index,sp_pm_index))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index) &
             ,sum(species_EMEP_subgrid(:,:,:,pm25_sp_index,1:sp_ppm_index))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index) &
             ,sum(comp_EMEP_subgrid(:,:,:,pm25_nc_index))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index)
         write(unit_logfile,'(A,3f12.2)') 'PM10 (init_sp,sum_sp,comp)',sum(species_EMEP_subgrid(:,:,:,pm10_sp_index,sp_pm_index))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index) &
             ,sum(species_EMEP_subgrid(:,:,:,pm10_sp_index,1:sp_ppm_index))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index) &
             ,sum(comp_EMEP_subgrid(:,:,:,pm10_nc_index))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index)
         !Replace the ppm values in the species with the emep read version
         !This is not correct. emep_subgrid_index,allsource_index is the sum of all local and nonlocal sources, so a larger number than the primary.
         !This will lead to negative values set to 0. So comment out this line (21.03.2022)
         !species_EMEP_subgrid(:,:,:,pm25_sp_index,sp_ppm_index)=subgrid(:,:,:,emep_subgrid_index,allsource_index,pollutant_loop_back_index(pm25_nc_index))
         !species_EMEP_subgrid(:,:,:,pm10_sp_index,sp_ppm_index)=subgrid(:,:,:,emep_subgrid_index,allsource_index,pollutant_loop_back_index(pm10_nc_index))
         !Normalise the species other than ppm again after subtracting the ppm
         !This should not change anything
         if (use_single_time_loop_flag) then
            sum_temp(:,:,:)=sum(species_EMEP_subgrid(:,:,:,pm25_sp_index,1:sp_ppm_index-1),4)
         else
            sum_temp(:,:,:)=sum(species_EMEP_subgrid(:,:,:,pm25_sp_index,1:sp_ppm_index-1),4)
         endif
         do i_sp=1,sp_ppm_index-1
             species_EMEP_subgrid(:,:,:,pm25_sp_index,i_sp)=species_EMEP_subgrid(:,:,:,pm25_sp_index,i_sp) &
                *(comp_EMEP_subgrid(:,:,:,pm25_nc_index)-species_EMEP_subgrid(:,:,:,pm25_sp_index,sp_ppm_index)) &
                /sum_temp
            where (sum_temp.le.0) species_EMEP_subgrid(:,:,:,pm25_sp_index,i_sp)=0
            where (isnan(species_EMEP_subgrid(:,:,:,pm25_sp_index,i_sp)).or.species_EMEP_subgrid(:,:,:,pm25_sp_index,i_sp).le.0) species_EMEP_subgrid(:,:,:,pm25_sp_index,i_sp)=0
         enddo
         if (use_single_time_loop_flag) then
            sum_temp(:,:,:)=sum(species_EMEP_subgrid(:,:,:,pm10_sp_index,1:sp_ppm_index-1),4)
         else
            sum_temp(:,:,:)=sum(species_EMEP_subgrid(:,:,:,pm10_sp_index,1:sp_ppm_index-1),4)
         endif
         do i_sp=1,sp_ppm_index-1
            species_EMEP_subgrid(:,:,:,pm10_sp_index,i_sp)=species_EMEP_subgrid(:,:,:,pm10_sp_index,i_sp) &
                *(comp_EMEP_subgrid(:,:,:,pm10_nc_index)-species_EMEP_subgrid(:,:,:,pm10_sp_index,sp_ppm_index)) &
                /sum_temp
            where (sum_temp.le.0) species_EMEP_subgrid(:,:,:,pm10_sp_index,i_sp)=0
            where (isnan(species_EMEP_subgrid(:,:,:,pm10_sp_index,i_sp)).or.species_EMEP_subgrid(:,:,:,pm10_sp_index,i_sp).le.0) species_EMEP_subgrid(:,:,:,pm10_sp_index,i_sp)=0
         enddo
         
         
         !Remove the primary based on the local contribution.
         species_EMEP_subgrid(:,:,:,pm25_sp_index,sp_ppm_index)=species_EMEP_subgrid(:,:,:,pm25_sp_index,sp_ppm_index)-subgrid(:,:,:,emep_local_subgrid_index,allsource_index,pollutant_loop_back_index(pm25_nc_index))
         species_EMEP_subgrid(:,:,:,pm10_sp_index,sp_ppm_index)=species_EMEP_subgrid(:,:,:,pm10_sp_index,sp_ppm_index)-subgrid(:,:,:,emep_local_subgrid_index,allsource_index,pollutant_loop_back_index(pm10_nc_index))
         where (species_EMEP_subgrid(:,:,:,pm25_sp_index,sp_ppm_index).lt.0.or.isnan(species_EMEP_subgrid(:,:,:,pm25_sp_index,sp_ppm_index))) species_EMEP_subgrid(:,:,:,pm25_sp_index,sp_ppm_index)=0
         where (species_EMEP_subgrid(:,:,:,pm10_sp_index,sp_ppm_index).lt.0.or.isnan(species_EMEP_subgrid(:,:,:,pm10_sp_index,sp_ppm_index))) species_EMEP_subgrid(:,:,:,pm10_sp_index,sp_ppm_index)=0
         
         !Normalise again after the subtraction to the nonlocal contribution (21.03.2022)
         sum_temp(:,:,:)=sum(species_EMEP_subgrid(:,:,:,pm25_sp_index,1:sp_ppm_index),4)
         do i_sp=1,sp_ppm_index
             species_EMEP_subgrid(:,:,:,pm25_sp_index,i_sp)=species_EMEP_subgrid(:,:,:,pm25_sp_index,i_sp) &
                *(subgrid(:,:,:,emep_nonlocal_subgrid_index,allsource_index,pollutant_loop_back_index(pm25_nc_index))) &
                /sum_temp
            where (sum_temp.le.0) species_EMEP_subgrid(:,:,:,pm25_sp_index,i_sp)=0
            where (isnan(species_EMEP_subgrid(:,:,:,pm25_sp_index,i_sp)).or.species_EMEP_subgrid(:,:,:,pm25_sp_index,i_sp).le.0) species_EMEP_subgrid(:,:,:,pm25_sp_index,i_sp)=0
         enddo
         if (use_single_time_loop_flag) then
             sum_temp(:,:,:)=sum(species_EMEP_subgrid(:,:,:,pm10_sp_index,1:sp_ppm_index),4)
         else
             sum_temp(:,:,:)=sum(species_EMEP_subgrid(:,:,:,pm10_sp_index,1:sp_ppm_index),4)
         endif
         do i_sp=1,sp_ppm_index
            species_EMEP_subgrid(:,:,:,pm10_sp_index,i_sp)=species_EMEP_subgrid(:,:,:,pm10_sp_index,i_sp) &
                *(subgrid(:,:,:,emep_nonlocal_subgrid_index,allsource_index,pollutant_loop_back_index(pm10_nc_index))) &
                /sum_temp
            where (sum_temp.le.0) species_EMEP_subgrid(:,:,:,pm10_sp_index,i_sp)=0
            where (isnan(species_EMEP_subgrid(:,:,:,pm10_sp_index,i_sp)).or.species_EMEP_subgrid(:,:,:,pm10_sp_index,i_sp).le.0) species_EMEP_subgrid(:,:,:,pm10_sp_index,i_sp)=0
         enddo
         
         
         write(unit_logfile,'(A)') 'init_sp should be the same as before, sum_sp should be the same as nonlocal, comp should be larger'
         write(unit_logfile,'(A,4f12.2)') 'PM25 (init_sp,sum_sp,comp,nonlocal)',sum(species_EMEP_subgrid(:,:,:,pm25_sp_index,sp_pm_index))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index) &
             ,sum(species_EMEP_subgrid(:,:,:,pm25_sp_index,1:sp_ppm_index))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index) &
             ,sum(comp_EMEP_subgrid(:,:,:,pm25_nc_index))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index) &
             ,sum(subgrid(:,:,:,emep_nonlocal_subgrid_index,allsource_index,pollutant_loop_back_index(pm25_nc_index)))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index)
         write(unit_logfile,'(A,4f12.2)') 'PM10 (init_sp,sum_sp,comp,nonlocal)',sum(species_EMEP_subgrid(:,:,:,pm10_sp_index,sp_pm_index))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index) &
             ,sum(species_EMEP_subgrid(:,:,:,pm10_sp_index,1:sp_ppm_index))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index) &
             ,sum(comp_EMEP_subgrid(:,:,:,pm10_nc_index))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index) &
             ,sum(subgrid(:,:,:,emep_nonlocal_subgrid_index,allsource_index,pollutant_loop_back_index(pm10_nc_index)))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index)
    endif
     

    !Only show results where all the subgrids and all sources are valid
    !This is temporary
    if (interpolate_subgrids_flag) then
    do source_index=1,n_source_index
        if (calculate_source(source_index)) then
            do j=1,subgrid_dim(y_dim_index)
            do i=1,subgrid_dim(x_dim_index)
                if (.not.use_subgrid(i,j,source_index)) then
                    !subgrid(i,j,:,total_subgrid_index,allsource_index,:)=NODATA_value
                    !subgrid(i,j,:,total_subgrid_index,source_index,:)=NODATA_value
                    !comp_subgrid(i,j,:,:)=NODATA_value
                endif                
            enddo
            enddo
        endif
    enddo
    endif
        
    end subroutine uEMEP_combine_local_source

end module redistribute_data

