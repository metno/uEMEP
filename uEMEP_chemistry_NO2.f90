
    subroutine uEMEP_chemistry_control
    
    use uEMEP_definitions

    implicit none
    
    real, allocatable :: subgrid_dummy(:,:,:,:,:,:)
    real, allocatable :: comp_subgrid_dummy(:,:,:,:)
    real, allocatable :: comp_source_subgrid_dummy(:,:,:,:,:)
    real, allocatable :: comp_source_EMEP_subgrid_dummy(:,:,:,:,:)
    real, allocatable :: comp_source_EMEP_additional_subgrid_dummy(:,:,:,:,:)
    
    integer in_region_loop, n_in_region_loop

    !These are calculated in the Chemistry routine. Fist declared here. Are global variables
    if (.not.allocated(comp_source_subgrid)) allocate(comp_source_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index,n_source_index))
    if (.not.allocated(comp_source_EMEP_subgrid)) allocate(comp_source_EMEP_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index,n_source_index))
    if (.not.allocated(comp_source_EMEP_additional_subgrid)) allocate(comp_source_EMEP_additional_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index,n_source_index))

    if (trace_emissions_from_in_region) then
        n_in_region_loop=2
        if (.not.allocated(subgrid_dummy)) allocate (subgrid_dummy(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_subgrid_index,n_source_index,n_pollutant_loop))
        if (.not.allocated(comp_subgrid_dummy))allocate (comp_subgrid_dummy(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index))
        subgrid_dummy=0 
        comp_subgrid_dummy=0 
        if (.not.allocated(comp_source_subgrid_dummy)) allocate(comp_source_subgrid_dummy(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index,n_source_index))
        if (.not.allocated(comp_source_EMEP_subgrid_dummy)) allocate(comp_source_EMEP_subgrid_dummy(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index,n_source_index))
        if (.not.allocated(comp_source_EMEP_additional_subgrid_dummy)) allocate(comp_source_EMEP_additional_subgrid_dummy(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index,n_source_index))
        comp_source_subgrid_dummy=0
        comp_source_EMEP_subgrid_dummy=0
        comp_source_EMEP_additional_subgrid_dummy=0

    else
        n_in_region_loop=1        
    endif
    
    do in_region_loop=1,n_in_region_loop

        write(unit_logfile,'(a)')''
        write(unit_logfile,'(a)')'--------------------------'
        if (in_region_loop.eq.1) write(unit_logfile,'(a)')'Chemistry for all contributions'
        if (in_region_loop.eq.2) write(unit_logfile,'(a)')'Chemistry only for inside regional contributions'
        write(unit_logfile,'(a)')'--------------------------'
        
       !Loop over the normal subgrid and the in region subgrid by saving  to a dummy variable and pputting back when finished
        if (trace_emissions_from_in_region.and.in_region_loop.eq.2) then        
            subgrid_dummy=subgrid
            comp_subgrid_dummy=comp_subgrid
            
            subgrid=subgrid_from_in_region
            comp_subgrid=comp_subgrid_from_in_region
            !comp_source_subgrid=comp_source_subgrid_from_in_region
            !comp_source_EMEP_subgrid=comp_source_EMEP_subgrid_from_in_region
            !comp_source_EMEP_additional_subgrid=comp_source_EMEP_additional_subgrid_from_in_region
            
            !These are calculated in the chemistry routine
            comp_source_subgrid_dummy=comp_source_subgrid
            comp_source_EMEP_subgrid_dummy=comp_source_EMEP_subgrid
            comp_source_EMEP_additional_subgrid_dummy=comp_source_EMEP_additional_subgrid
        endif

        call uEMEP_chemistry
        
        if (trace_emissions_from_in_region.and.in_region_loop.eq.2) then
            subgrid_from_in_region=subgrid
            comp_subgrid_from_in_region=comp_subgrid
            comp_source_subgrid_from_in_region=comp_source_subgrid
            comp_source_EMEP_subgrid_from_in_region=comp_source_EMEP_subgrid
            comp_source_EMEP_additional_subgrid_from_in_region=comp_source_EMEP_additional_subgrid

            subgrid=subgrid_dummy
            comp_subgrid=comp_subgrid_dummy
            comp_source_subgrid=comp_source_subgrid_dummy
            comp_source_EMEP_subgrid=comp_source_EMEP_subgrid_dummy
            comp_source_EMEP_additional_subgrid=comp_source_EMEP_additional_subgrid_dummy
            
        endif

    enddo !from_in_region loop
    
    if (trace_emissions_from_in_region) then
        if (allocated(subgrid_dummy)) deallocate (subgrid_dummy)
        if (allocated(comp_subgrid_dummy)) deallocate (comp_subgrid_dummy)
        if (allocated(comp_source_subgrid_dummy)) deallocate(comp_source_subgrid_dummy)
        if (allocated(comp_source_EMEP_subgrid_dummy)) deallocate(comp_source_EMEP_subgrid_dummy)
        if (allocated(comp_source_EMEP_additional_subgrid_dummy)) deallocate(comp_source_EMEP_additional_subgrid_dummy)
    endif
    
    end subroutine uEMEP_chemistry_control
    
    
    
    subroutine uEMEP_chemistry
    !Routine for doing the chemistry calculations in uEMEP
    
    use uEMEP_definitions

    implicit none
    
    integer i,j
    real nox_bg,no2_bg,o3_bg,nox_loc,f_no2_loc,J_photo,temperature
    real nox_out,no2_out,o3_out,p_bg_out,p_out
    
    integer t,t_start,t_end
    integer i_source,i_subsource,emep_subsource
    integer i_pollutant
    logical :: nox_available=.false.
    !integer i_comp,i_file
    !character(256) temp_name
    integer i_integral,j_integral
    real FF_loc,distance_grid
    integer i_cross_integral,j_cross_integral,i_nc,j_nc
    real sum_p_bg_out,sum_p_out,count_p_out
    real max_p_bg_out,max_p_out,min_p_bg_out,min_p_out
    
    !NB. Additional is calculated but not necessarily saved!
    real nox_bg_additional,no2_bg_additional,o3_bg_additional
    
    !Search for nox in the pollutants
    do i_pollutant=1,n_pollutant_loop
        if (pollutant_loop_index(i_pollutant).eq.nox_nc_index) nox_available=.true.
    enddo
    
    !Leave the chemistry routine if nox is not available
    if (.not.nox_available) return  

    write(unit_logfile,'(A)') ''
    write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Calculating chemistry for NO2 (uEMEP_chemistry)'
	write(unit_logfile,'(A)') '================================================================'
    if (no2_chemistry_scheme_flag.eq.0) then
        write(unit_logfile,'(A)') 'No chemistry used'
    elseif (no2_chemistry_scheme_flag.eq.1) then
        write(unit_logfile,'(A)') 'Photostationary state used'
    elseif (no2_chemistry_scheme_flag.eq.2) then
        write(unit_logfile,'(A)') 'Photochemistry with time scale used'
    elseif (no2_chemistry_scheme_flag.eq.3) then
        write(unit_logfile,'(A)') 'Romberg parameterisation used'
    elseif (no2_chemistry_scheme_flag.eq.4) then
        write(unit_logfile,'(A)') 'SRM parameterisation used'
    elseif (no2_chemistry_scheme_flag.eq.5) then
        write(unit_logfile,'(A)') 'During parameterisation used'
    endif
    
    t_start=1
    t_end=subgrid_dim(t_dim_index)
    i_subsource=1
    emep_subsource=1
    comp_subgrid(:,:,:,no2_index)=0
    comp_subgrid(:,:,:,nox_index)=0
    comp_subgrid(:,:,:,o3_index)=0

    nox_bg=0.;no2_bg=0.;o3_bg=0.;nox_loc=0.;f_no2_loc=0.;J_photo=0.;temperature=0.;
        
    !Before calculating travel time then include the other EMEP sources not downscaled
    !Travel time is set to EMEP Grid_width/FFgrid
    do t=t_start,t_end
    do j=1,subgrid_dim(y_dim_index)
    do i=1,subgrid_dim(x_dim_index)
        i_cross_integral=crossreference_target_to_integral_subgrid(i,j,x_dim_index)
        j_cross_integral=crossreference_target_to_integral_subgrid(i,j,y_dim_index)
        FF_loc=1.
        if (hourly_calculations) then
            FF_loc=max(FF_min_dispersion,meteo_subgrid(i_cross_integral,j_cross_integral,t,FFgrid_subgrid_index))
        elseif (annual_calculations) then
            FF_loc=max(FF_min_dispersion,1./meteo_subgrid(i_cross_integral,j_cross_integral,t,inv_FFgrid_subgrid_index))
        endif
        
        i_nc=crossreference_target_to_emep_subgrid(i,j,x_dim_index)
        j_nc=crossreference_target_to_emep_subgrid(i,j,y_dim_index)

        if (EMEP_projection_type.eq.LL_projection_index) then
            distance_grid=111000.*sqrt(dgrid_nc(lon_nc_index)*cos(var1d_nc(j_nc,lat_nc_index)*pi/180.)*dgrid_nc(lat_nc_index))
        else
        !Assumed LCC or PS
            distance_grid=sqrt(dgrid_nc(lon_nc_index)*dgrid_nc(lat_nc_index))
        endif

        do i_source=1,n_source_index
        if (calculate_emep_source(i_source).and..not.calculate_source(i_source)) then
            !traveltime_subgrid(i,j,t,1,:)=traveltime_subgrid(i,j,t,1,:) &
            !    +distance_grid/FF_loc*subgrid(i,j,t,emep_local_subgrid_index,i_source,:)**traveltime_power
            !traveltime_subgrid(i,j,t,2,:)=traveltime_subgrid(i,j,t,2,:)+subgrid(i,j,t,emep_local_subgrid_index,i_source,:)**traveltime_power
        !write(*,'(3i,4f12.2)') i,j,i_source,distance_grid,FF_loc,distance_grid/FF_loc,subgrid(i,j,t,emep_local_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index))
        endif
        enddo
    enddo
    enddo
    enddo    
    
    !Calculate the weighted travel time from the totals calculated in uEMEP_subgrid_dispersion
    !do t=t_start,t_end
        !traveltime_subgrid(:,:,t,3,:)=traveltime_subgrid(:,:,t,1,:)/traveltime_subgrid(:,:,t,2,:)
        !Invert it to get the time scale
        !traveltime_subgrid(:,:,t,1)=1./traveltime_subgrid(:,:,t,1)
        !Set none valid to 12 hours (long time)
        !where (traveltime_subgrid(:,:,t,2,:).eq.0) traveltime_subgrid(:,:,t,3,:)=3600.*12.
        !write(*,*) t
        !write(*,*) traveltime_subgrid(:,:,t,1)
    !enddo
    
    do j=1,subgrid_dim(y_dim_index)
    do i=1,subgrid_dim(x_dim_index)
    !write(*,'(2i4,<subgrid_dim(t_dim_index)>f6.1)') i,j,traveltime_subgrid(i,j,:,1)/60.
    enddo
    enddo
    
    sum_p_bg_out=0.
    sum_p_out=0.
    count_p_out=0
    max_p_bg_out=-1000.;min_p_bg_out=1000.;max_p_out=-1000.;min_p_out=1000.
    do t=t_start,t_end
    do j=1,subgrid_dim(y_dim_index)
    do i=1,subgrid_dim(x_dim_index)
    if (use_subgrid(i,j,allsource_index)) then
        
        i_integral=crossreference_target_to_integral_subgrid(i,j,x_dim_index)
        j_integral=crossreference_target_to_integral_subgrid(i,j,y_dim_index)

        J_photo=meteo_subgrid(i_integral,j_integral,t,J_subgrid_index)
        temperature=meteo_subgrid(i_integral,j_integral,t,t2m_subgrid_index)

        nox_bg=subgrid(i,j,t,emep_nonlocal_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))
        
        !Do not do this because the EMEP contributions are already added into the local contribution
        !Add the additional non-downscaled EMEP local source to the background as these have not been used to calulate travel times
        !do i_source=1,n_source_index
        !if (calculate_EMEP_source(i_source).and..not.calculate_source(i_source).and.i_source.ne.allsource_index) then
        !    nox_bg=nox_bg+subgrid(i,j,t,emep_local_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index))
        !endif
        !enddo
        
        if (EMEP_additional_grid_interpolation_size.gt.0) then
        nox_bg_additional=subgrid(i,j,t,emep_additional_nonlocal_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))
        !do i_source=1,n_source_index
        !if (calculate_EMEP_source(i_source).and..not.calculate_source(i_source).and.i_source.ne.allsource_index) then
        !    nox_bg_additional=nox_bg_additional+subgrid(i,j,t,emep_additional_local_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index))
        !endif
        !enddo
        !subgrid(i,j,t,emep_additional_nonlocal_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))=nox_bg_additional
        endif
        
        !nox_bg=subgrid(i,j,t,emep_subgrid_index,allsource_index,emep_subsource)
        !nox_bg=comp_subgrid(i,j,t,nox_index)*(14.+16.*2.)/14.
        
        o3_bg=comp_EMEP_subgrid(i,j,t,o3_index)
        !o3_bg=comp_EMEP_subgrid(i,j,t,o3_index)+48./46.*comp_EMEP_subgrid(i,j,t,no2_index)*subgrid(i,j,t,emep_local_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))/subgrid(i,j,t,emep_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))
       
        f_no2_loc=0.
        nox_loc=0.
        
        do i_source=1,n_source_index
        if (calculate_source(i_source)) then
            do i_subsource=1,n_subsource(i_source)
            f_no2_loc=f_no2_loc+emission_factor(no2_index,i_source,i_subsource)/emission_factor(nox_index,i_source,i_subsource)*subgrid(i,j,t,local_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index))
            nox_loc=nox_loc+subgrid(i,j,t,local_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index))
            enddo         
        endif
        if (calculate_emep_source(i_source).and..not.calculate_source(i_source)) then
            do i_subsource=1,n_subsource(i_source)
            f_no2_loc=f_no2_loc+f_no2_emep*subgrid(i,j,t,emep_local_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index))
            nox_loc=nox_loc+subgrid(i,j,t,emep_local_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index))
            enddo         
        endif
        enddo
        if (nox_loc.ne.0.0) then
            f_no2_loc=f_no2_loc/nox_loc
        else
            f_no2_loc=0.
        endif
        !no2_bg=max(0.,comp_subgrid(i,j,t,no2_index)-f_no2_loc*subgrid(i,j,t,emep_local_subgrid_index,allsource_index,emep_subsource))
        !no2_bg=comp_EMEP_subgrid(i,j,t,no2_index)*subgrid(i,j,t,emep_nonlocal_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))/subgrid(i,j,t,emep_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))
        no2_bg=comp_EMEP_subgrid(i,j,t,no2_index)*nox_bg/subgrid(i,j,t,emep_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))
        o3_bg=max(0.0,comp_EMEP_subgrid(i,j,t,o3_index)+48./46.*(comp_EMEP_subgrid(i,j,t,no2_index)-no2_bg)) !Conserve Ox when removing NO2 in the background. Cannot be less than 0
        !no2_bg=comp_subgrid(i,j,t,no2_index)

        !Assume stationary state to derive no2 and o3 background
        if (no2_background_chemistry_scheme_flag.eq.1) then
            call uEMEP_nonlocal_NO2_O3(nox_bg,comp_EMEP_subgrid(i,j,t,nox_index),comp_EMEP_subgrid(i,j,t,no2_index),comp_EMEP_subgrid(i,j,t,o3_index),J_photo,temperature,f_no2_emep,no2_bg,o3_bg)
        endif
        
        if (EMEP_additional_grid_interpolation_size.gt.0) then
            no2_bg_additional=comp_EMEP_subgrid(i,j,t,no2_index)*nox_bg_additional/subgrid(i,j,t,emep_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))
            if (no2_background_chemistry_scheme_flag.eq.1) then             
                call uEMEP_nonlocal_NO2_O3(nox_bg_additional,comp_EMEP_subgrid(i,j,t,nox_index),comp_EMEP_subgrid(i,j,t,no2_index),comp_EMEP_subgrid(i,j,t,o3_index),J_photo,temperature,f_no2_emep,no2_bg_additional,o3_bg_additional)
            else
                !o3_bg_additional=comp_EMEP_subgrid(i,j,t,o3_index)
                o3_bg_additional=max(0.0,comp_EMEP_subgrid(i,j,t,o3_index)+48./46.*(comp_EMEP_subgrid(i,j,t,no2_index)-no2_bg_additional)) !Conserve Ox when removing NO2 in the background
            endif
            comp_source_EMEP_additional_subgrid(i,j,t,o3_index,allsource_index)=o3_bg_additional
            comp_source_EMEP_additional_subgrid(i,j,t,no2_index,allsource_index)=no2_bg_additional
        endif
        
        !Set the background O3 level. use all_source for the nonlocal.
        !subgrid(i,j,t,emep_nonlocal_subgrid_index,allsource_index,pollutant_loop_back_index(o3_nc_index))=o3_bg
        comp_source_EMEP_subgrid(i,j,t,o3_index,allsource_index)=o3_bg
        comp_source_EMEP_subgrid(i,j,t,no2_index,allsource_index)=no2_bg
        
        !write(*,*) nox_bg,subgrid(i,j,t,emep_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))
       !if (J_photo.ne.0) then
        !write(*,'(A,3I6,f5.2,6f12.3)') 'IN : ',i,j,t,no2_bg/nox_bg,nox_bg,no2_bg,o3_bg,nox_loc,f_no2_loc,J_photo
       !endif
        
        !if (nox_loc+nox_bg.lt.0.) write(*,*) t,nox_loc,nox_bg,nox_loc+nox_bg
       
        
        if (no2_chemistry_scheme_flag.eq.0) then
            nox_out=nox_bg+nox_loc
            no2_out=no2_bg+nox_loc*f_no2_loc
            o3_out=o3_bg
        elseif (no2_chemistry_scheme_flag.eq.1) then
            call uEMEP_photostationary_NO2(nox_bg,no2_bg,o3_bg,nox_loc,f_no2_loc,J_photo,temperature,nox_out,no2_out,o3_out,p_bg_out,p_out)
        elseif (no2_chemistry_scheme_flag.eq.2) then
            !write(*,'(7f8.2,f12.2,2f8.2)') nox_bg,no2_bg,o3_bg,nox_loc,f_no2_loc,J_photo,temperature,traveltime_subgrid(i,j,t,1,pollutant_loop_back_index(nox_index))
            call uEMEP_phototimescale_NO2(nox_bg,no2_bg,o3_bg,nox_loc,f_no2_loc,J_photo,temperature,traveltime_subgrid(i,j,t,3,pollutant_loop_back_index(nox_nc_index))*traveltime_scaling,nox_out,no2_out,o3_out,p_bg_out,p_out)
            !write(*,'(7f8.2,f12.2,2f8.2)') nox_bg,no2_bg,o3_bg,nox_loc,f_no2_loc,J_photo,temperature,traveltime_subgrid(i,j,t,1),no2_out/nox_out,o3_out/o3_bg
        elseif (no2_chemistry_scheme_flag.eq.3) then
            call uEMEP_Romberg_NO2(nox_bg,no2_bg,nox_loc,o3_bg,f_no2_loc,nox_out,no2_out,o3_out,romberg_parameters)
            !write(*,'(8f8.3)') nox_bg,no2_bg,nox_loc,o3_bg,f_no2_loc,nox_out,no2_out,o3_out
        elseif (no2_chemistry_scheme_flag.eq.4) then
            call uEMEP_SRM_NO2(nox_bg,no2_bg,nox_loc,o3_bg,f_no2_loc,nox_out,no2_out,o3_out,SRM_parameters)
        elseif (no2_chemistry_scheme_flag.eq.5) then
            call uEMEP_During_NO2(nox_bg,no2_bg,nox_loc,o3_bg,f_no2_loc,comp_EMEP_subgrid(i,j,t,nox_index),comp_EMEP_subgrid(i,j,t,no2_index),comp_EMEP_subgrid(i,j,t,o3_index),J_photo,temperature,nox_out,no2_out,o3_out,p_bg_out,p_out)
        endif
        
        sum_p_bg_out=sum_p_bg_out+p_bg_out
        sum_p_out=sum_p_out+p_out
        count_p_out=count_p_out+1
        max_p_bg_out=max(max_p_bg_out,p_bg_out);min_p_bg_out=min(min_p_bg_out,p_bg_out)
        max_p_out=max(max_p_out,p_out);min_p_out=min(min_p_out,p_out)
        
        !write(*,*) nox_out-subgrid(i,j,t,total_subgrid_index,allsource_index,1)
        
        comp_subgrid(i,j,t,o3_index)=o3_out
        comp_subgrid(i,j,t,no2_index)=no2_out
        comp_subgrid(i,j,t,nox_index)=nox_out
       
        
        !if (J_photo.ne.0) then
        !write(*,'(A,3I6,f5.2,6f12.3)') 'OUT: ',i,j,t,no2_out/nox_out,nox_out,no2_out,o3_out,p_bg_out,p_out
        !write(*,'(A,5I6,f5.2,3f12.3,1es12.2,2f12.3)') 'IN: ',i,j,i_integral,j_integral,t,no2_bg/nox_bg,nox_bg,no2_bg,o3_bg,J_photo,p_bg_out,p_out
       ! write(*,*)
        !endif
        !if (no2_out.gt.nox_out) write(*,*)no2_out,nox_out,nox_bg,nox_loc
    !else
    !    comp_subgrid(i,j,t,:)=0.
    !endif
    
    else
        comp_subgrid(i,j,t,o3_index)=NODATA_value
        comp_subgrid(i,j,t,no2_index)=NODATA_value
        comp_subgrid(i,j,t,nox_index)=NODATA_value
        
    endif
    
    enddo
    enddo
    enddo

    write(*,'(A,2f12.3)') 'P value (nonlocal,local): ',sum_p_bg_out/count_p_out,sum_p_out/count_p_out
    write(*,'(A,2f12.3)') 'P max (nonlocal,local): ',max_p_bg_out,max_p_out
    write(*,'(A,2f12.3)') 'P min (nonlocal,local): ',min_p_bg_out,min_p_out

    
    end subroutine uEMEP_chemistry
   
    subroutine uEMEP_source_fraction_chemistry
    !Special source allocation for no2 based on leaving out one source at a time in the chemistry calculation
    !This will always give a sum less, but not much less than, the total no2
    !This is normalised in order for it to be used
    !Vhemistry scheme must have been run prior to implementing this
    
    use uEMEP_definitions

    implicit none
    
    integer i,j
    real nox_bg,no2_bg,o3_bg,nox_loc,f_no2_loc,J_photo,temperature
    real nox_out,no2_out,o3_out,p_bg_out,p_out
    
    integer t,t_start,t_end
    integer i_source,i_subsource,emep_subsource
    integer i_pollutant
    logical :: nox_available=.false.
    !integer i_comp,i_file
    !character(256) temp_name
    integer i_integral,j_integral
    integer remove_source
    real sum_no2_source_subgrid,sum_o3_source_subgrid
    real, allocatable  :: comp_source_temp_subgrid(:,:,:,:,:)
    real, allocatable  :: comp_source_EMEP_temp_subgrid(:,:,:,:,:)
    
    !Search for nox in the pollutants
    do i_pollutant=1,n_pollutant_loop
        if (pollutant_loop_index(i_pollutant).eq.nox_nc_index) nox_available=.true.
    enddo
    
    !Leave the chemistry routine if nox is not available
    if (.not.nox_available) return  

    if (.not.allocated(comp_source_subgrid)) allocate(comp_source_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index,n_source_index))
    
    if (calculate_EMEP_additional_grid_flag) then
        if (.not.allocated(comp_source_additional_subgrid)) allocate(comp_source_additional_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index,n_source_index))
        !Temporary array for storing the comp_source_subgrid to avoid rewriting large parts of the routine when running the additional version
        if (.not.allocated(comp_source_temp_subgrid)) allocate(comp_source_temp_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index,n_source_index))
        if (.not.allocated(comp_source_EMEP_temp_subgrid)) allocate(comp_source_EMEP_temp_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index,n_source_index))
        comp_source_temp_subgrid=comp_source_subgrid
        comp_source_EMEP_temp_subgrid=comp_source_EMEP_subgrid
        comp_source_EMEP_subgrid=comp_source_EMEP_additional_subgrid
    endif
    ! Already allocated in chemistry call
    !if (.not.allocated(comp_source_EMEP_subgrid)) allocate(comp_source_EMEP_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index,n_source_index))

    !write(unit_logfile,'(A)') ''
    !write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Calculating chemistry source contribution for NO2 (uEMEP_source_fraction_chemistry)'
	!write(unit_logfile,'(A)') '================================================================'
    if (no2_chemistry_scheme_flag.eq.0) then
        write(unit_logfile,'(A)') 'No chemistry used'
    elseif (no2_chemistry_scheme_flag.eq.1) then
        write(unit_logfile,'(A)') 'Photostationary state used'
    elseif (no2_chemistry_scheme_flag.eq.2) then
        write(unit_logfile,'(A)') 'Photochemistry with time scale used'
    elseif (no2_chemistry_scheme_flag.eq.3) then
        write(unit_logfile,'(A)') 'Romberg parameterisation used'
    elseif (no2_chemistry_scheme_flag.eq.4) then
        write(unit_logfile,'(A)') 'SRM parameterisation used'
    elseif (no2_chemistry_scheme_flag.eq.5) then
        write(unit_logfile,'(A)') 'During parameterisation used'
    endif
    
    t_start=1
    t_end=subgrid_dim(t_dim_index)
    i_subsource=1
    emep_subsource=1

    nox_bg=0.;no2_bg=0.;o3_bg=0.;nox_loc=0.;f_no2_loc=0.;J_photo=0.;temperature=0.;
        
    !Weighted travel time already calculated
    
    do t=t_start,t_end
    do j=1,subgrid_dim(y_dim_index)
    do i=1,subgrid_dim(x_dim_index)
    if (use_subgrid(i,j,allsource_index)) then
        
        i_integral=crossreference_target_to_integral_subgrid(i,j,x_dim_index)
        j_integral=crossreference_target_to_integral_subgrid(i,j,y_dim_index)

        J_photo=meteo_subgrid(i_integral,j_integral,t,J_subgrid_index)
        temperature=meteo_subgrid(i_integral,j_integral,t,t2m_subgrid_index)
       
        
        if (calculate_EMEP_additional_grid_flag) then
            nox_bg=subgrid(i,j,t,emep_additional_nonlocal_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))
        else
            nox_bg=subgrid(i,j,t,emep_nonlocal_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))
        endif
        
        !Add the additional non-downscaled EMEP local source to the background as these have not been used to calulate travel times
        !do i_source=1,n_source_index
        !if (calculate_EMEP_source(i_source).and..not.calculate_source(i_source).and.i_source.ne.allsource_index) then
        !    nox_bg=nox_bg+subgrid(i,j,t,emep_local_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index))
        !endif
        !enddo
         
        o3_bg=comp_EMEP_subgrid(i,j,t,o3_index)
        !o3_bg=comp_EMEP_subgrid(i,j,t,o3_index)+48./46.*comp_EMEP_subgrid(i,j,t,no2_index)*subgrid(i,j,t,emep_local_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))/subgrid(i,j,t,emep_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))
            
        do remove_source=1,n_source_index
        if (calculate_source(remove_source).or.remove_source.eq.allsource_index.or.(calculate_emep_source(remove_source).and..not.calculate_source(remove_source))) then
        !if (calculate_source(remove_source).or.remove_source.eq.allsource_index.or.calculate_emep_source(remove_source)) then
        
            f_no2_loc=0.
            nox_loc=0.

            do i_source=1,n_source_index
            if (calculate_source(i_source)) then
                if (remove_source.ne.i_source) then
                    do i_subsource=1,n_subsource(i_source)
                    f_no2_loc=f_no2_loc+emission_factor(no2_index,i_source,i_subsource)/emission_factor(nox_index,i_source,i_subsource)*subgrid(i,j,t,local_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index))
                    nox_loc=nox_loc+subgrid(i,j,t,local_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index))
                    enddo       
                endif
            endif             

            !Include the local EMEP that are not being downscaled
            if (.not.calculate_EMEP_additional_grid_flag) then
            if (calculate_emep_source(i_source).and..not.calculate_source(i_source)) then
                if (remove_source.ne.i_source) then
                    do i_subsource=1,n_subsource(i_source)
                    f_no2_loc=f_no2_loc+f_no2_emep*subgrid(i,j,t,emep_local_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index))
                    nox_loc=nox_loc+subgrid(i,j,t,emep_local_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index))
                    enddo
                endif
            endif
            endif
            
            if (calculate_EMEP_additional_grid_flag) then
            !If calculating the additional region then use the additional local EMEP not being downscaled 
            if (calculate_emep_source(i_source).and..not.calculate_source(i_source)) then
                if (remove_source.ne.i_source) then
                    do i_subsource=1,n_subsource(i_source)
                    f_no2_loc=f_no2_loc+f_no2_emep*subgrid(i,j,t,emep_additional_local_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index))
                    nox_loc=nox_loc+subgrid(i,j,t,emep_additional_local_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index))
                    enddo
                endif
            endif
            !If calculating the additional region then include the difference BG-BG_additional to the local EMEP that is being downscaled
            if (calculate_source(i_source)) then
                if (remove_source.ne.i_source) then
                    do i_subsource=1,n_subsource(i_source)
                    f_no2_loc=f_no2_loc+f_no2_emep* &
                        (subgrid(i,j,t,emep_nonlocal_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index)) &
                        -subgrid(i,j,t,emep_additional_nonlocal_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index)))
                    nox_loc=nox_loc+subgrid(i,j,t,emep_nonlocal_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index)) &
                        -subgrid(i,j,t,emep_additional_nonlocal_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index))
                    enddo
                endif
            endif                
            endif
            
            enddo

            if (nox_loc.ne.0.0) then
                f_no2_loc=f_no2_loc/nox_loc
            else
                f_no2_loc=0.
            endif
 
            
            !Use the all source index to calculate the contribution from the background
            !This is done by removing all the sources, rather than the difference as done for the local sources
            !This is because the chemistry is disturbed when removing background nox and no2
            if (remove_source.ne.allsource_index) then
                !nox_bg=subgrid(i,j,t,emep_nonlocal_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))
                !no2_bg=comp_EMEP_subgrid(i,j,t,no2_index)*subgrid(i,j,t,emep_nonlocal_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))/subgrid(i,j,t,emep_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))       
                no2_bg=comp_EMEP_subgrid(i,j,t,no2_index)*nox_bg/subgrid(i,j,t,emep_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))
                o3_bg=max(0.0,comp_EMEP_subgrid(i,j,t,o3_index)+48./46.*(comp_EMEP_subgrid(i,j,t,no2_index)-no2_bg)) !Conserve Ox when removing NO2 in the background. Cannot be less than 0

            else
                !nox_bg=subgrid(i,j,t,emep_nonlocal_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))
                !no2_bg=comp_EMEP_subgrid(i,j,t,no2_index)*subgrid(i,j,t,emep_nonlocal_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))/subgrid(i,j,t,emep_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))       
                no2_bg=comp_EMEP_subgrid(i,j,t,no2_index)*nox_bg/subgrid(i,j,t,emep_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))
                o3_bg=max(0.0,comp_EMEP_subgrid(i,j,t,o3_index)+48./46.*(comp_EMEP_subgrid(i,j,t,no2_index)-no2_bg)) !Conserve Ox when removing NO2 in the background. Cannot be less than 0
                nox_loc=0.
                f_no2_loc=0.
            endif
            
        !Assume stationary state to derive no2 and o3 background. Overwrites the previous setting
        if (no2_background_chemistry_scheme_flag.eq.1) then
            call uEMEP_nonlocal_NO2_O3(nox_bg,comp_EMEP_subgrid(i,j,t,nox_index),comp_EMEP_subgrid(i,j,t,no2_index),comp_EMEP_subgrid(i,j,t,o3_index),J_photo,temperature,f_no2_emep,no2_bg,o3_bg)
        endif
        
            if (no2_chemistry_scheme_flag.eq.0) then
                nox_out=nox_bg+nox_loc
                no2_out=no2_bg+nox_loc*f_no2_loc
                o3_out=o3_bg
            elseif (no2_chemistry_scheme_flag.eq.1) then
                call uEMEP_photostationary_NO2(nox_bg,no2_bg,o3_bg,nox_loc,f_no2_loc,J_photo,temperature,nox_out,no2_out,o3_out,p_bg_out,p_out)
            elseif (no2_chemistry_scheme_flag.eq.2) then
                !write(*,'(i,7f8.2,f12.2,2f8.2)') remove_source,nox_bg,no2_bg,o3_bg,nox_loc,f_no2_loc,J_photo,temperature,traveltime_subgrid(i,j,t,1,pollutant_loop_back_index(nox_nc_index))
                call uEMEP_phototimescale_NO2(nox_bg,no2_bg,o3_bg,nox_loc,f_no2_loc,J_photo,temperature,traveltime_subgrid(i,j,t,3,pollutant_loop_back_index(nox_nc_index))*traveltime_scaling,nox_out,no2_out,o3_out,p_bg_out,p_out)
                !write(*,'(i,7f8.2,f12.2,3f8.4)') remove_source,nox_bg,no2_bg,o3_bg,nox_loc,f_no2_loc,J_photo,temperature,traveltime_subgrid(i,j,t,1,pollutant_loop_back_index(nox_nc_index)),nox_out,no2_out,o3_out
            elseif (no2_chemistry_scheme_flag.eq.3) then
                call uEMEP_Romberg_NO2(nox_bg,no2_bg,nox_loc,o3_bg,f_no2_loc,nox_out,no2_out,o3_out,romberg_parameters)   
            elseif (no2_chemistry_scheme_flag.eq.4) then
                call uEMEP_SRM_NO2(nox_bg,no2_bg,nox_loc,o3_bg,f_no2_loc,nox_out,no2_out,o3_out,SRM_parameters)
            elseif (no2_chemistry_scheme_flag.eq.5) then
                call uEMEP_During_NO2(nox_bg,no2_bg,nox_loc,o3_bg,f_no2_loc,comp_EMEP_subgrid(i,j,t,nox_index),comp_EMEP_subgrid(i,j,t,no2_index),comp_EMEP_subgrid(i,j,t,o3_index),J_photo,temperature,nox_out,no2_out,o3_out,p_bg_out,p_out)
            endif
        
            !write(*,*) nox_out-subgrid(i,j,t,total_subgrid_index,allsource_index,1)
        
            !For background just use the result without any sources.
            !There is a problem disturbing the chemistry by removing the background nox and no2 but not changing the o3
            if (remove_source.eq.allsource_index) then
                comp_source_subgrid(i,j,t,no2_index,remove_source)=no2_bg
                comp_source_subgrid(i,j,t,o3_index,remove_source)=o3_bg
            else
                !Avoid round off errors which can occur with small numbers
                comp_source_subgrid(i,j,t,no2_index,remove_source)=max(0.0,comp_subgrid(i,j,t,no2_index)-no2_out)
                !Can be negative and can be greater than 1 so do not limit
                comp_source_subgrid(i,j,t,o3_index,remove_source)=comp_subgrid(i,j,t,o3_index)-o3_out
                !comp_source_EMEP_subgrid(i,j,t,no2_index,remove_source)=max(0.0,comp_subgrid(i,j,t,no2_index)-no2_out)
                !Can be negative and can be greater than 1 so do not limit
                !comp_source_EMEP_subgrid(i,j,t,o3_index,remove_source)=comp_subgrid(i,j,t,o3_index)-o3_out
                !write(*,*) i,j,comp_source_subgrid(i,j,t,o3_index,remove_source)
            endif
      
        endif
        enddo
        
    else

        comp_source_subgrid(i,j,t,:,:)=NODATA_value
        
    endif
    
    !Normalise the contributions
    
    !Calculate the sum
            sum_no2_source_subgrid=0.
            sum_o3_source_subgrid=0.
            do i_source=1,n_source_index
            !if (calculate_source(i_source).or.i_source.eq.allsource_index) then
            !if (calculate_source(i_source)) then
            if (calculate_source(i_source).or.(calculate_emep_source(i_source).and..not.calculate_source(i_source))) then
                sum_no2_source_subgrid=sum_no2_source_subgrid+comp_source_subgrid(i,j,t,no2_index,i_source)
                sum_o3_source_subgrid=sum_o3_source_subgrid+comp_source_subgrid(i,j,t,o3_index,i_source)
            endif
            !if (calculate_emep_source(i_source).and..not.calculate_source(i_source)) then
            !    sum_no2_source_subgrid=sum_no2_source_subgrid+comp_source_EMEP_subgrid(i,j,t,no2_index,i_source)
            !    sum_o3_source_subgrid=sum_o3_source_subgrid+comp_source_EMEP_subgrid(i,j,t,o3_index,i_source)
            !endif
            enddo
            !Set the background fractions so they will not be adjusted with normalisation
            !comp_source_subgrid(i,j,t,no2_index,allsource_index)=comp_source_subgrid(i,j,t,no2_index,allsource_index)/comp_subgrid(i,j,t,no2_index)
            !comp_source_subgrid(i,j,t,o3_index,allsource_index)=comp_source_subgrid(i,j,t,o3_index,allsource_index)/comp_subgrid(i,j,t,o3_index)
            do i_source=1,n_source_index
            !if (calculate_source(i_source)) then
            if (calculate_source(i_source).or.(calculate_emep_source(i_source).and..not.calculate_source(i_source))) then
                !Adjust for the background and normalise
                if (sum_no2_source_subgrid.ne.0) then
                    comp_source_subgrid(i,j,t,no2_index,i_source)=comp_source_subgrid(i,j,t,no2_index,i_source)/sum_no2_source_subgrid &
                     *(comp_subgrid(i,j,t,no2_index)-comp_source_EMEP_subgrid(i,j,t,no2_index,allsource_index))
                else
                    comp_source_subgrid(i,j,t,no2_index,i_source)=0
                    !comp_source_subgrid(i,j,t,no2_index,allsource_index)=(comp_subgrid(i,j,t,no2_index)-comp_source_EMEP_subgrid(i,j,t,no2_index,allsource_index))
                endif
                !write(*,*) i,j,i_source,sum_o3_source_subgrid
                if (sum_o3_source_subgrid.ne.0) then
                    comp_source_subgrid(i,j,t,o3_index,i_source)=comp_source_subgrid(i,j,t,o3_index,i_source)/sum_o3_source_subgrid &
                     *(comp_subgrid(i,j,t,o3_index)-comp_source_EMEP_subgrid(i,j,t,o3_index,allsource_index))
                else
                    comp_source_subgrid(i,j,t,o3_index,i_source)=0
                    !comp_source_subgrid(i,j,t,o3_index,allsource_index)=(comp_subgrid(i,j,t,o3_index)-comp_source_EMEP_subgrid(i,j,t,o3_index,allsource_index))
                endif
                if (comp_subgrid(i,j,t,no2_index).le.0) comp_source_subgrid(i,j,t,no2_index,i_source)=0
                if (comp_subgrid(i,j,t,o3_index).le.0) comp_source_subgrid(i,j,t,o3_index,i_source)=0
                !if (comp_source_subgrid(i,j,t,no2_index,i_source).le.0) comp_source_subgrid(i,j,t,no2_index,i_source)=0.
                !if (comp_source_subgrid(i,j,t,o3_index,i_source).le.0) comp_source_subgrid(i,j,t,o3_index,i_source)=0
            
            endif
            enddo
            !write(*,'(2i4,6f12.6)') i,j,sum_no2_source_subgrid,sum_no2_source_subgrid/comp_subgrid(i,j,t,no2_index),comp_source_subgrid(i,j,t,no2_index,allsource_index) &
            !    ,comp_source_subgrid(i,j,t,no2_index,traffic_index),comp_source_subgrid(i,j,t,no2_index,shipping_index),comp_source_subgrid(i,j,t,no2_index,heating_index)
            !write(*,'(2i4,6f12.6)') i,j,sum_o3_source_subgrid,sum_o3_source_subgrid/comp_subgrid(i,j,t,o3_index),comp_source_subgrid(i,j,t,o3_index,allsource_index) &
            !    ,comp_source_subgrid(i,j,t,o3_index,traffic_index),comp_source_subgrid(i,j,t,o3_index,shipping_index),comp_source_subgrid(i,j,t,o3_index,heating_index)
            !if (comp_source_subgrid(i,j,t,nor_index,allsource_index).lt.0.or.comp_source_subgrid(i,j,t,no2_index,traffic_index).lt.0.or.no2_source_subgrid(i,j,t,shipping_index).lt.0.or.no2_source_fraction_subgrid(i,j,t,heating_index).lt.0) then
             !   write(*,*) 'Traffic value less than 0. comp_subgrid =',comp_subgrid(i,j,t,no2_index),comp_EMEP_subgrid(i,j,t,no2_index) &
             !       ,comp_EMEP_subgrid(i,j,t,no2_index)*subgrid(i,j,t,emep_nonlocal_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))/subgrid(i,j,t,emep_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))  &
             !       ,comp_EMEP_subgrid(i,j,t,o3_index)
            !    stop
            !endif
            
    enddo
    enddo
    enddo

    !Transfer the arrays to the right outputs
    if (calculate_EMEP_additional_grid_flag) then
        comp_source_additional_subgrid=comp_source_subgrid
        comp_source_subgrid=comp_source_temp_subgrid
        comp_source_EMEP_subgrid=comp_source_EMEP_temp_subgrid
        !EMEP_additional is unchanged
        if (allocated(comp_source_temp_subgrid)) deallocate(comp_source_temp_subgrid)
        if (allocated(comp_source_EMEP_temp_subgrid)) deallocate(comp_source_EMEP_temp_subgrid)        
    endif

    
    end subroutine uEMEP_source_fraction_chemistry
    
    
    subroutine uEMEP_photostationary_NO2(nox_bg,no2_bg,o3_bg,nox_loc,f_no2_loc,J_photo,temperature,nox_out,no2_out,o3_out,p_bg_out,p_out)

    !use uEMEP_definitions
    
    implicit none
    
    real, intent(in) :: nox_bg,no2_bg,o3_bg,nox_loc,f_no2_loc,J_photo,temperature
    real, intent(out) :: nox_out,no2_out,o3_out,p_bg_out,p_out
    integer no2_i,no_i,nox_i,o3_i,ox_i,nox_bg_i,no2_bg_i,n_i
    parameter (n_i=7)
    real Na,Na_fac,k1
    real mass(n_i)
    real mmass(n_i)
    real mol(n_i)
    real f_no2,f_ox,Jd,fac_sqrt
    real :: min_nox=1.0e-6

    DATA mmass /46.,30.,46.,48.,47.,46.,46./

    no2_i=1;no_i=2;nox_i=3;o3_i=4;ox_i=5;nox_bg_i=6;no2_bg_i=7
           

    Na=6.022e23        !(molecules/mol)
    Na_fac=Na/1.0e12   !Conversion from ug/m3 to molecules/cm3 included


    k1=1.4e-12*exp(-1310./temperature); !(cm^3/s) and temperature in Kelvin

    !mmass(1:n_i)=(/46.,30.,46.,48.,47.,46.,46.,46./)
    mass(1:n_i)=0.
    mol(1:n_i)=0.

    !Test for 0 NOx. If so leave the routine
    mass(nox_i)=nox_loc+nox_bg   
!    if (mass(nox_i).eq.0.) then
    if (mass(nox_i).le.min_nox) then
        nox_out=0.
        no2_out=0.
        o3_out=o3_bg
        return
    endif
    
    !Check the photostationary assumption for the input data
    mass(nox_i)=nox_bg
    mass(no2_i)=no2_bg
    mass(o3_i)=o3_bg
    mol=mass/mmass*Na_fac !(molecules per cm3)
    mol(ox_i)=mol(o3_i)+mol(no2_i)
    mol(no_i)=max(0.0,mol(nox_i)-mol(no2_i))
    !Test the photostationary state for the bg input data
    if (J_photo.ne.0.and.mol(no_i).ne.0..and.mol(o3_i).ne.0.) then
         p_bg_out=J_photo*mol(no2_i)/k1/mol(o3_i)/mol(no_i)
    else
        !p_bg_out=mol(no2_i)/mol(ox_i)
        p_bg_out=mol(no2_i)/(mol(ox_i)+mol(nox_i)-abs(mol(ox_i)-mol(nox_i)))*2.
    endif
   
    !if (J_photo.ne.0.) write(*,*) p_bg_out,J_photo,mol(no2_i),k1,mol(o3_i),mol(no_i)
    
    !Add the local contribution for calculation
    mass(nox_i)=nox_loc+nox_bg
    mass(no2_i)=f_no2_loc*nox_loc+no2_bg
    mass(o3_i)=o3_bg
    !mass(ox_i)=o3_bg+mass(no2_i)
    mol=mass/mmass*Na_fac !(molecules per cm3)

    mol(ox_i)=mol(o3_i)+mol(no2_i)
    mol(no_i)=max(0.0,mol(nox_i)-mol(no2_i))
    
    f_no2=mol(no2_i)/mol(nox_i)
    f_ox=mol(ox_i)/mol(nox_i)
    
    !Test the photostationary state for the input data. Will not be in equilibrium
    if (J_photo.ne.0.and.mol(no_i).ne.0..and.mol(o3_i).ne.0.) then
        p_out=J_photo*mol(no2_i)/k1/mol(o3_i)/mol(no_i)
    else
        !p_out=mol(no2_i)/mol(ox_i)
        p_out=mol(no2_i)/(mol(ox_i)+mol(nox_i)-abs(mol(ox_i)-mol(nox_i)))*2.
    endif

    !Set the photolysis rate
    Jd=J_photo/k1/mol(nox_i)
    
    !Calculate fraction of NO2 in photostationary state
    fac_sqrt=max(0.0,(1.+f_ox+Jd)**2-4.*f_ox)
    !if (J_photo.ne.0) then
        f_no2=0.5*((1.+f_ox+Jd)-sqrt(fac_sqrt))
    !else
    !    f_no2=min(1.0,f_ox)
    !endif
    
    !write(*,'(A,9ES12.1)') 'MOL:  ',mol(nox_i),mol(no2_i),mol(o3_i),mol(ox_i),f_no2,f_ox,Jd,f_no2,p_bg_out

    !Convert back to mass
    mol(no2_i)=f_no2*mol(nox_i);
    mol(o3_i)=max(0.,mol(ox_i)-mol(no2_i))  !Rounding errors possible
    mol(no_i)=max(0.,mol(nox_i)-mol(no2_i)) !Rounding errors possible
    mass=mol*mmass/Na_fac !(ug/m3)


    no2_out=mass(no2_i)
    nox_out=mass(nox_i)
    o3_out=mass(o3_i)

    !Check output
    if (J_photo.ne.0.and.mol(no_i).ne.0..and.mol(o3_i).ne.0.) then
        p_out=J_photo*mol(no2_i)/k1/mol(o3_i)/mol(no_i)
    else
        !p_out=mol(no2_i)/mol(ox_i)
        p_out=mol(no2_i)/(mol(ox_i)+mol(nox_i)-abs(mol(ox_i)-mol(nox_i)))*2.
    endif
    
    !write(*,'(A,9ES12.1)') 'MASS: ',mass(nox_i),mass(no2_i),mass(o3_i),mass(ox_i),f_no2,f_ox,Jd,f_no2,p_out

    end subroutine uEMEP_photostationary_NO2
    
    
    subroutine uEMEP_phototimescale_NO2(nox_bg,no2_bg,o3_bg,nox_loc,f_no2_loc,J_photo,temperature,time_scale,nox_out,no2_out,o3_out,p_bg_out,p_out)

    !use uEMEP_definitions
    
    implicit none
    
    real, intent(in) :: nox_bg,no2_bg,o3_bg,nox_loc,f_no2_loc,J_photo,temperature,time_scale
    real, intent(out) :: nox_out,no2_out,o3_out,p_bg_out,p_out
    integer no2_i,no_i,nox_i,o3_i,ox_i,nox_bg_i,no2_bg_i,n_i
    parameter (n_i=7)
    real Na,Na_fac,k1
    real mass(n_i)
    real mmass(n_i)
    real mol(n_i)
    real fac_sqrt
    real f_no2,f_ox,Jd,Jd_bg
    real :: min_nox=1.0e-6
    real c,b,BB,td,f_no2_0,f_no2_ps
    complex(4) AA
    real p_tot_out,f_ox_bg,f_no2_bg_ps,f_no2_bg
    DATA mmass /46.,30.,46.,48.,47.,46.,46./

    no2_i=1;no_i=2;nox_i=3;o3_i=4;ox_i=5;nox_bg_i=6;no2_bg_i=7
           

    Na=6.022e23        !(molecules/mol)
    Na_fac=Na/1.0e12   !Conversion from ug/m3 to molecules/cm3 included


    k1=1.4e-12*exp(-1310./temperature); !(cm^3/s) and temperature in Kelvin

    !mmass(1:n_i)=(/46.,30.,46.,48.,47.,46.,46.,46./)
    mass(1:n_i)=0.
    mol(1:n_i)=0.

    !Test for 0 NOx. If so leave the routine
    mass(nox_i)=nox_loc+nox_bg   
!    if (mass(nox_i).eq.0.) then
    if (mass(nox_i).le.min_nox) then
        nox_out=0.
        no2_out=0.
        o3_out=o3_bg
        return
    endif
    
    !Check the photostationary assumption for the input data
    mass(nox_i)=nox_bg
    mass(no2_i)=no2_bg
    mass(o3_i)=o3_bg
    mol=mass/mmass*Na_fac !(molecules per cm3)

    mol(ox_i)=mol(o3_i)+mol(no2_i)
    mol(no_i)=max(0.0,mol(nox_i)-mol(no2_i))
    f_ox_bg=mol(ox_i)/mol(nox_i)
    Jd_bg=J_photo/k1/mol(nox_i)
    f_no2_bg_ps=0.5*((1+f_ox_bg+Jd_bg)-sqrt((1+f_ox_bg+Jd_bg)**2-4.*f_ox_bg))
    f_no2_bg=mol(no2_i)/mol(nox_i)
    p_bg_out=f_no2_bg/f_no2_bg_ps
    
    !Check input
    if (J_photo.ne.0.and.mol(no_i).ne.0..and.mol(o3_i).ne.0.) then
        p_bg_out=J_photo*mol(no2_i)/k1/mol(o3_i)/mol(no_i)
    else
        !p_bg_out=mol(no2_i)/mol(ox_i)
        p_bg_out=mol(no2_i)/(mol(ox_i)+mol(nox_i)-abs(mol(ox_i)-mol(nox_i)))*2.
    endif
    
    !if (J_photo.ne.0.) write(*,*) p_bg_out,J_photo,mol(no2_i),k1,mol(o3_i),mol(no_i)
    
    !Add the local contribution for calculation
    mass(nox_i)=nox_loc+nox_bg
    mass(no2_i)=f_no2_loc*nox_loc+no2_bg
    mass(o3_i)=o3_bg
    !mass(ox_i)=o3_bg+mass(no2_i)
    mol=mass/mmass*Na_fac !(molecules per cm3)
    mol(ox_i)=mol(o3_i)+mol(no2_i)
    mol(no_i)=max(0.0,mol(nox_i)-mol(no2_i))
    
    f_no2=mol(no2_i)/mol(nox_i)
    f_ox=mol(ox_i)/mol(nox_i)

    !Set the photolysis rate
    Jd=J_photo/k1/mol(nox_i)

    !Calculate photostationary for total nox, ox
    fac_sqrt=max(0.,(1+f_ox+Jd)**2-4.*f_ox)
    f_no2_ps=0.5*((1+f_ox+Jd)-sqrt(fac_sqrt))
    p_tot_out=f_no2/f_no2_ps

    !Test the photostationary state for the input data
    !if (J_photo.ne.0.and.mol(no_i).ne.0..and.mol(o3_i).ne.0.) then
    !    p_out=J_photo*mol(no2_i)/k1/mol(o3_i)/mol(no_i)
    !else
    !    p_out=mol(no2_i)/mol(nox_i)
    !endif

    
    !Calculate fraction of NO2 based on the time scale
    !fac_sqrt=max(0.0,(1+f_ox+Jd)**2-4*f_ox)
    c=f_ox
    b=1+f_ox+Jd
    BB=sqrt(max(0.0,b**2-4.*c))!max avoids roundoff errors
    td=time_scale*k1*mol(nox_i)
    f_no2_0=f_no2

    !if (J_photo.ne.0) then
        !f_no2=0.5*((1+f_ox+Jd)-sqrt(fac_sqrt))
    !if (b.lt.100.)
    AA=clog(cmplx((BB+b-2.*f_no2_0)/(BB-b+2.*f_no2_0)))
    f_no2=real(-BB/2.*((exp(AA+BB*td)-1.)/(exp(AA+BB*td)+1.))+b/2.)
    !if (BB*td.gt.50.) f_no2=-BB/2.+b/2.
    if (isnan(f_no2)) f_no2=-BB/2.+b/2.

    !write(*,*) f_no2,AA,nox_loc,f_no2_loc,AA,BB,Jd,mol(nox_i),k1

    fac_sqrt=max(0.0,(1+f_ox+Jd)**2-4.*f_ox)
    f_no2_ps=0.5*((1+f_ox+Jd)-sqrt(fac_sqrt))
    p_out=f_no2/f_no2_ps
    !write(*,*) p_bg_out,p_tot_out,p_out,AA,BB,b,td,exp(AA+BB*td),f_ox
    !write(*,*) o3_bg,no2_bg,nox_bg
    !write(*,*) c,b,BB,AA,(BB/2.+b/2.-f_no2_0)/(BB/2.-b/2.+f_no2_0),f_no2
    
    !else
    !    f_no2=1.0
    !endif
    
    !write(*,'(A,9ES12.1)') 'MOL:  ',mol(nox_i),mol(no2_i),mol(o3_i),mol(ox_i),f_no2,f_ox,Jd,f_no2,p_bg_out

    !Convert back to mass
    mol(no2_i)=max(0.,f_no2*mol(nox_i))
    mol(o3_i)=max(0.,mol(ox_i)-mol(no2_i))  !Rounding errors possible
    mol(no_i)=max(0.,mol(nox_i)-mol(no2_i)) !Rounding errors possible
    mass=mol*mmass/Na_fac !(ug/m3)


    no2_out=mass(no2_i)
    nox_out=mass(nox_i)
    o3_out=mass(o3_i)
    
    if (isnan(no2_out)) then
        write(*,'(8a12)') 'nox_bg','no2_bg','o3_bg','nox_loc','f_no2_loc','J_photo','temperature','time_scale'
        write(*,'(8es12.2)') nox_bg,no2_bg,o3_bg,nox_loc,f_no2_loc,J_photo,temperature,time_scale
        write(*,'(4a12)') 'f_no2','BB','b','b**2-4.*c'
        write(*,'(4es12.2)') f_no2,BB,b,b**2-4.*c
        stop
    endif
    

    !Check output
    if (J_photo.ne.0.and.mol(no_i).ne.0..and.mol(o3_i).ne.0.) then
        p_out=J_photo*mol(no2_i)/k1/mol(o3_i)/mol(no_i)
    else
        !p_out=mol(no2_i)/mol(ox_i)
        p_out=mol(no2_i)/(mol(ox_i)+mol(nox_i)-abs(mol(ox_i)-mol(nox_i)))*2.
    endif

    !write(*,'(6f12.2)') 1./(k1*mol(nox_i)),1./J_photo,time_scale,f_no2/f_no2_ps,p_out,nox_loc/(nox_loc+nox_bg)


    !write(*,'(A,9ES12.1)') 'MASS: ',mass(nox_i),mass(no2_i),mass(o3_i),mass(ox_i),f_no2,f_ox,Jd,f_no2,p_out

    end subroutine uEMEP_phototimescale_NO2
    
    subroutine uEMEP_Romberg_NO2(nox_bg,no2_bg,nox_loc,o3_bg,f_no2_loc,nox_out,no2_out,o3_out,romberg_parameters)
    
    implicit none
    
    real, intent(in) :: nox_bg,no2_bg,nox_loc,o3_bg,f_no2_loc
    real, intent(in) :: romberg_parameters(3)
    real, intent(out) :: nox_out,no2_out,o3_out
    
    !From Norwegian obs fit
    !real :: a_rom=20
    !real :: b_rom=30
    !From model fit
    real :: a_rom=30
    real :: b_rom=35
    real :: c_rom=0.20
    real ox_init,no2_init,no2_equ
    real beta,F,K
    !Gral values 30 35 0.18
    !Bächlin and Bösinger (2008) 29 35 0.217
    
    if (romberg_parameters(1).ne.0) then
        a_rom=romberg_parameters(1)
        b_rom=romberg_parameters(2)
        c_rom=romberg_parameters(3)
    endif
    
    nox_out=nox_bg+nox_loc
    no2_equ=a_rom*nox_bg/(nox_bg+b_rom)+nox_bg*c_rom
    no2_out=a_rom*nox_out/(nox_out+b_rom)+nox_out*c_rom
    no2_out=no2_out-no2_equ+no2_bg
    no2_out=max(no2_bg,no2_out)

    no2_init=no2_bg+f_no2_loc*nox_loc
    
    !Small adjustments for molecular weights
    ox_init=no2_init*47./46.+o3_bg*47./48.
    o3_out=ox_init*48./47.-no2_out*48./46.
    
    !o3_out=o3_bg.+no2_bg*48./46.-no2_out*48./46.

    end subroutine uEMEP_Romberg_NO2

    subroutine uEMEP_SRM_NO2(nox_bg,no2_bg,nox_loc,o3_bg,f_no2_loc,nox_out,no2_out,o3_out,SRM_parameters)
    
    implicit none
    
    real, intent(in) :: nox_bg,no2_bg,nox_loc,o3_bg,f_no2_loc
    real, intent(in) :: SRM_parameters(3)
    real, intent(out) :: nox_out,no2_out,o3_out
    
    !From model fit
    real :: beta=0.45
    real :: K=30
    real :: F=0.2
    real ox_init,no2_init,no2_equ

    !From RIVM Briefrapport 2014-0109
    !beta=1
    !K=100
    !F=0.2
    
    !Reference
    !https://core.ac.uk/download/pdf/58774365.pdf
    
    if (SRM_parameters(1).ne.0) then
        beta=SRM_parameters(1)
        K=SRM_parameters(2)
        F=SRM_parameters(3)
    endif
    
    nox_out=nox_bg+nox_loc
    no2_out=no2_bg+beta*o3_bg*nox_loc/(nox_loc+K/(1-F))+F*nox_loc

    no2_init=no2_bg+f_no2_loc*nox_loc
    
    !Small adjustments for molecular weights
    ox_init=no2_init*47./46.+o3_bg*47./48.
    o3_out=ox_init*48./47.-no2_out*48./46.
    
    !o3_out=o3_bg.+no2_bg*48./46.-no2_out*48./46.
    
    end subroutine uEMEP_SRM_NO2

    subroutine uEMEP_During_NO2(nox_bg,no2_bg,nox_loc,o3_bg,f_no2_loc,nox_emep,no2_emep,o3_emep,J_photo,temperature,nox_out,no2_out,o3_out,p_bg_out,p_out)
    
    implicit none
    
    real, intent(in) :: nox_bg,no2_bg,nox_loc,o3_bg,f_no2_loc,J_photo,temperature
    real, intent(in) :: nox_emep,no2_emep,o3_emep
    real, intent(out) :: nox_out,no2_out,o3_out
    real, intent(out) :: p_bg_out,p_out
    
    real :: mol_nox_bg,mol_no2_bg,mol_nox_loc,mol_o3_bg,mol_no2_loc,mol_ox_loc,mol_no_bg,mol_ox_bg
    real :: mol_nox_out,mol_no2_out,mol_o3_out,mol_no_out
    real :: mol_nox_emep,mol_no2_emep,mol_o3_emep,mol_no_emep,mol_ox_emep,p_emep_out
    real :: b,d,r,c,k1
    real :: Na,Na_fac
    integer no2_i,no_i,nox_i,o3_i,ox_i,n_i
    parameter (n_i=5)
    real mmass(n_i)
    DATA mmass /46.,30.,46.,48.,47./

    no2_i=1;no_i=2;nox_i=3;o3_i=4;ox_i=5

    !Reference
    !Düring, I., Bächlin, W., Ketzel, M., Baum, A., Friedrich, U., Wurzler, S., 2011.
    !A new simplified NO/NO2 conversion model under consideration of direct NO2-emissions.
    !Meteorol. Zeitschrift 20, 67–73. doi:10.1127/0941-2948/2011/0491
    
    !Improved Methodologies for NO2 Exposure Assessment in the EU, page 53
    !https://ec.europa.eu/environment/air/pdf/NO2_Exposure_Final_Report.pdf
    
    k1=1.4e-12*exp(-1310./temperature); !(cm^3/s) and temperature in Kelvin

    Na=6.022e23        !(molecules/mol)
    Na_fac=Na/1.0e12   !Conversion from ug/m3 to molecules/cm3 included
    !Normally multiplied by *Na_fac but not necessary as it is just a scaling
        
    mol_no2_bg=no2_bg/mmass(no2_i)
    mol_no2_loc=f_no2_loc*nox_loc/mmass(no2_i)
    mol_nox_bg=nox_bg/mmass(nox_i)
    mol_nox_loc=nox_loc/mmass(nox_i)
    mol_o3_bg=o3_bg/mmass(o3_i)
    mol_no_bg=(nox_bg-no2_bg)/mmass(nox_i)
    mol_ox_bg=mol_o3_bg+mol_no2_bg
    
    mol_o3_emep=o3_emep/mmass(o3_i)
    mol_nox_emep=nox_emep/mmass(nox_i)
    mol_no2_emep=no2_emep/mmass(no2_i)
    mol_no_emep=max(0.,mol_nox_emep-mol_no2_emep)
    mol_ox_emep=mol_o3_emep+mol_no2_emep
    
    mol_ox_loc=mol_o3_bg+mol_no2_bg+mol_no2_loc

    if (mol_no2_emep.gt.0) then
        r=mol_o3_emep*mol_no_emep/mol_no2_emep
    else
        r=0.
    endif  
    
    b=mol_ox_loc+mol_nox_bg+mol_nox_loc+r
    c=max(0.,b**2-4.*mol_ox_loc*(mol_nox_bg+mol_nox_loc)) !Should never be less than 0 but can be -0
    d=sqrt(c)
    mol_no2_out=(b-d)/2.
    mol_o3_out=mol_ox_loc-mol_no2_out
    mol_no_out=mol_nox_bg+mol_nox_loc-mol_no2_out
    
    nox_out=nox_bg+nox_loc    
    no2_out=mol_no2_out*mmass(no2_i)
    o3_out=mol_o3_out*mmass(o3_i)

    !Not correct as it does not calculate the actual photostationary equation
    p_out=r !/Na_fac
    p_bg_out=mol_o3_bg*mol_no_bg/mol_no2_bg !/Na_fac
    !write(*,*) r,nox_out,no2_out,o3_out

        !Check output
    if (J_photo.ne.0.and.mol_no_out.ne.0..and.mol_o3_out.ne.0.) then
        p_out=J_photo*mol_no2_out/k1/mol_o3_out/mol_no_out/Na_fac
        p_emep_out=J_photo*mol_no2_emep/k1/mol_o3_emep/mol_no_emep/Na_fac
        p_bg_out=J_photo*mol_no2_bg/k1/mol_o3_bg/mol_no_bg/Na_fac
    else
        p_out=mol_no2_out/(mol_ox_loc+mol_nox_bg+mol_nox_loc-abs(mol_ox_loc-mol_nox_bg-mol_nox_loc))*2.
        p_emep_out=mol_no2_emep/(mol_ox_emep+mol_nox_emep-abs(mol_ox_emep-mol_nox_emep))*2.
        p_bg_out=mol_no2_bg/(mol_ox_bg+mol_nox_bg-abs(mol_ox_bg-mol_nox_bg))*2.
    endif
    !p_out and p_emep are the same, or should be. Tested that way
    !p_out=p_emep_out

    end subroutine uEMEP_During_NO2

    subroutine uEMEP_nonlocal_NO2_O3(nox_bg,nox_emep,no2_emep,o3_emep,J_photo,temperature,f_no2,no2_out,o3_out)
    
    implicit none
    
    real, intent(in) :: J_photo,temperature,f_no2
    real, intent(in) :: nox_bg
    real, intent(in) :: nox_emep,no2_emep,o3_emep
    real, intent(out) :: no2_out,o3_out
    
    real :: mol_nox_bg
    real :: mol_no2_out,mol_o3_out
    real :: mol_nox_emep,mol_no2_emep,mol_o3_emep,mol_ox_emep,mol_no_emep
    real :: b,d,r,c
    real :: Na,Na_fac,k1
    real :: p_phot,r_phot
    integer no2_i,no_i,nox_i,o3_i,ox_i,n_i
    parameter (n_i=5)
    real mmass(n_i)
    DATA mmass /46.,30.,46.,48.,47./

    no2_i=1;no_i=2;nox_i=3;o3_i=4;ox_i=5

    !Reference
    !
    Na=6.022e23        !(molecules/mol)
    Na_fac=Na/1.0e12   !Conversion from ug/m3 to molecules/cm3 included

    k1=1.4e-12*exp(-1310./temperature); !(cm^3/s) and temperature in Kelvin

    !Normally multiplied by *Na_fac but not necessary as it is just a scaling
    mol_o3_emep=o3_emep/mmass(o3_i)*Na_fac
    mol_no2_emep=no2_emep/mmass(no2_i)*Na_fac
    mol_nox_emep=nox_emep/mmass(nox_i)*Na_fac
    mol_ox_emep=mol_o3_emep+mol_no2_emep-f_no2*mol_nox_emep
    mol_nox_bg=nox_bg/mmass(nox_i)*Na_fac
    mol_no_emep=max(0.,mol_nox_emep-mol_no2_emep)
    
    if (mol_no2_emep.gt.0) then
        r=mol_o3_emep*mol_no_emep/mol_no2_emep
        p_phot=J_photo/k1*mol_no2_emep/mol_o3_emep/mol_no_emep
        r_phot=J_photo/k1/Na_fac
        r=r_phot*Na_fac
    
        b=mol_ox_emep+mol_nox_bg+r
        c=max(0.,b**2-4.*mol_ox_emep*mol_nox_bg) !Should never be less than 0 but can be -0.0
        d=sqrt(c)
        mol_no2_out=(b-d)/2.
        mol_o3_out=mol_ox_emep-mol_no2_out
    
        no2_out=max(0.0,mol_no2_out*mmass(no2_i)/Na_fac)
        o3_out=max(0.0,mol_o3_out*mmass(o3_i)/Na_fac)
    else
        no2_out=0.
        o3_out=o3_emep       
    endif
    
    !write(*,'(a,9f12.3)') 'BG(r,r_phot,p_phot,nox_emep,no2_emep,o3_emep,nox_bg,no2_out,o3_out): ',r/Na_fac,r_phot,p_phot,nox_emep,no2_emep,o3_emep,nox_bg,no2_out,o3_out
    
    end subroutine uEMEP_nonlocal_NO2_O3

    subroutine correct_annual_mean_chemistry
    
    use uEMEP_definitions

    implicit none

    integer i,j,t
    integer t_start,t_end
    integer i_integral,j_integral
    real o3_in,nox_in,no2_in,J_photo_in,temperature_in,lon_in,lat_in
    real ox_sigma_ratio_in,nox_sigma_ratio_in
    real o3_out,no2_out
    real sum_no2_in,sum_no2_out,mean_no2_in,mean_no2_out
    integer no2_count
    logical run_all_flag

    
    write(unit_logfile,'(A)') ''
    write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Correcting annual mean NO2 and O3 (correct_annual_mean_chemistry)'
	write(unit_logfile,'(A)') '================================================================'

    t_start=1
    t_end=subgrid_dim(t_dim_index)

    sum_no2_in=0
    sum_no2_out=0
    no2_count=0
    do t=t_start,t_end
        
        !Run the conversion routine once to get the Jd distribution which is saved. This is to save time as this is slow. Done in the centre of the grid
        if (quick_annual_mean_pdf_chemistry_correction) then
            run_all_flag=.false.
            i=subgrid_dim(x_dim_index)/2
            j=subgrid_dim(y_dim_index)/2
            o3_in=comp_subgrid(i,j,t,o3_index)
            no2_in=comp_subgrid(i,j,t,no2_index)
            nox_in=comp_subgrid(i,j,t,nox_index)
            i_integral=crossreference_target_to_integral_subgrid(i,j,x_dim_index)
            j_integral=crossreference_target_to_integral_subgrid(i,j,y_dim_index)
            J_photo_in=meteo_subgrid(i_integral,j_integral,t,J_subgrid_index)
            temperature_in=meteo_subgrid(i_integral,j_integral,t,t2m_subgrid_index)
            lon_in=lon_subgrid(i,j)
            lat_in=lat_subgrid(i,j)        
            ox_sigma_ratio_in=ox_sigma_ratio_pdf
            nox_sigma_ratio_in=nox_sigma_ratio_pdf
            call uEMEP_annual_mean_pdf_correction_NO2_O3(max_bin_pdf,log10_step_bin_pdf,.true.,no2_in,nox_in,o3_in,J_photo_in,temperature_in,ox_sigma_ratio_in,nox_sigma_ratio_in,lon_in,lat_in,no2_out,o3_out)
        else
            run_all_flag=.true.
        endif
        
    do j=1,subgrid_dim(y_dim_index)
    do i=1,subgrid_dim(x_dim_index)

        o3_in=comp_subgrid(i,j,t,o3_index)
        no2_in=comp_subgrid(i,j,t,no2_index)
        nox_in=comp_subgrid(i,j,t,nox_index)

        i_integral=crossreference_target_to_integral_subgrid(i,j,x_dim_index)
        j_integral=crossreference_target_to_integral_subgrid(i,j,y_dim_index)
        J_photo_in=meteo_subgrid(i_integral,j_integral,t,J_subgrid_index)
        temperature_in=meteo_subgrid(i_integral,j_integral,t,t2m_subgrid_index)

        lon_in=lon_subgrid(i,j)
        lat_in=lat_subgrid(i,j)
        
        ox_sigma_ratio_in=ox_sigma_ratio_pdf
        nox_sigma_ratio_in=nox_sigma_ratio_pdf
        
        if (o3_in.le.0.or.nox_in.le.0.or.no2_in.le.0) then
            o3_out=o3_in
            no2_out=no2_in
        else
            call uEMEP_annual_mean_pdf_correction_NO2_O3(max_bin_pdf,log10_step_bin_pdf,run_all_flag,no2_in,nox_in,o3_in,J_photo_in,temperature_in,ox_sigma_ratio_in,nox_sigma_ratio_in,lon_in,lat_in,no2_out,o3_out)
        endif

        comp_subgrid(i,j,t,o3_index)=o3_out
        comp_subgrid(i,j,t,no2_index)=no2_out

        !write(*,'(a)') 'no2_in,nox_in,o3_in,J_photo_in,temperature_in,ox_sigma_ratio_in,nox_sigma_ratio_in,lon_in,lat_in,no2_out,o3_out'
        !write(*,'(11es12.2)') no2_in,nox_in,o3_in,J_photo_in,temperature_in,ox_sigma_ratio_in,nox_sigma_ratio_in,lon_in,lat_in,no2_out,o3_out
        
        sum_no2_in=sum_no2_in+no2_in
        sum_no2_out=sum_no2_out+no2_out
        no2_count=no2_count+1
       if (isnan(no2_out)) then
            write(unit_logfile,'(a,2i,4f12.4)') 'NaN in pdf output. Stopping ',i,j,no2_in,no2_out,o3_in,o3_out
            stop
        endif
        
        
    enddo
    enddo
    enddo
    
    write(unit_logfile,'(a,f12.4)') 'Average NO2 scaling with pdf correction = ',sum_no2_out/sum_no2_in
    
    end subroutine correct_annual_mean_chemistry

    subroutine uEMEP_annual_mean_pdf_correction_NO2_O3(bin_max,delta_log10_bin,run_all,no2_in,nox_in,o3_in,J_photo_in,temperature_in,ox_sigma_ratio_in,nox_sigma_ratio_in,lon_in,lat_in,no2_out,o3_out)
    
    implicit none

    real, intent(in) :: bin_max,delta_log10_bin
    logical, intent(in) :: run_all
    real, intent(in) :: J_photo_in,temperature_in
    real, intent(in) :: no2_in,nox_in,o3_in
    real, intent(in) :: ox_sigma_ratio_in,nox_sigma_ratio_in
    real, intent(in) :: lon_in,lat_in
    real, intent(out) :: no2_out,o3_out

    real :: mol_nox_bg
    real :: mol_no2_out,mol_o3_out
    real :: mol_nox_emep,mol_no2_emep,mol_o3_emep,mol_ox_emep,mol_no_emep
    real :: Na,Na_fac,k1
    integer no2_i,no_i,nox_i,o3_i,ox_i,n_i
    parameter (no2_i=1,no_i=2,nox_i=3,o3_i=4,ox_i=5,n_i=5)
    real mmass(n_i)
    DATA mmass /46.,30.,46.,48.,47./

    real bin_min,log10_bin_max,log10_bin_min
    real, allocatable :: log10_bin(:),bin(:),delta_bin(:)
    integer n_bin
    real nox_sigma_ratio,ox_sigma_ratio
    real ox_sigma,ox_sig_sqr,ox_mu,nox_sigma,nox_sig_sqr,nox_mu
    real, allocatable :: y_nox(:),y_ox(:), y_Jd(:)
    integer nbin_Jd,bin_temp
    parameter (nbin_Jd=8760) !Hours in a year
    real, save :: y_Jd_acc(nbin_Jd)
    real mean_y_Jd_acc,y_Jd_acc_temp,y_Jd_acc_temp_log10
    real y_all_val,y_all_prob,y_all_sum,y_all_prob_sum,y_all,y_annual,y_scale
    real solar_net,azimuth_ang,zenith_ang
    real mol_nox,mol_no2,mol_ox,mol_o3,Jd
    double precision date_num
    integer i,j,k,l
    real :: pi=3.14159265
    integer date_array(6)
    real :: min_sigma_ratio=0.01
    real bin_temp2

    Na=6.022e23        !(molecules/mol)
    Na_fac=Na/1.0e12   !Conversion from ug/m3 to molecules/cm3 included

    k1=1.4e-12*exp(-1310./temperature_in) !(cm^3/s) and temperature in Kelvin
    mol_nox=nox_in*Na_fac/mmass(nox_i)
    mol_no2=no2_in*Na_fac/mmass(no2_i)
    mol_o3=o3_in*Na_fac/mmass(o3_i)
    mol_ox=mol_no2+mol_o3
    Jd=J_photo_in/k1

    no2_out=no2_in
    o3_out=o3_in
    
    !Create the bins for the pdf in (mol/cm3). ox, nox and J
    !bin_max=1000. !In ug/m^3
    !delta_log10_bin=0.05
    !delta_log10_bin=0.01
    bin_min=0.1
    log10_bin_max=log10(bin_max)
    log10_bin_min=log10(bin_min)
    n_bin=(log10_bin_max-log10_bin_min)/delta_log10_bin+1
    !Creates 80 bins with these settings
    
    allocate (log10_bin(n_bin))
    allocate (bin(n_bin))
    allocate (delta_bin(n_bin))
    
    do i=1,n_bin
        log10_bin(i)=log10_bin_min+(i-1)*delta_log10_bin
    enddo
    do i=1,n_bin
        bin(i)=(10.**log10_bin(i))*Na_fac/mmass(nox_i)
        delta_bin(i)=(10.**(log10_bin(i)+delta_log10_bin/2.)-10.**(log10_bin(i)-delta_log10_bin/2.))*Na_fac/mmass(nox_i)
    enddo
    
    !Distribute concentrations into the pdf bins
    !Minimum needed to avoid NaNs in the calculation
    nox_sigma_ratio=1.14
    ox_sigma_ratio=0.21
    if (nox_sigma_ratio_in.ne.0) nox_sigma_ratio=max(nox_sigma_ratio_in,min_sigma_ratio)
    if (ox_sigma_ratio_in.ne.0) ox_sigma_ratio=max(ox_sigma_ratio_in,min_sigma_ratio)
    
    ox_sigma=mol_ox*ox_sigma_ratio
    ox_sig_sqr=log(1+ox_sigma**2./mol_ox**2.)
    ox_mu=log(mol_ox**2./sqrt(mol_ox**2.+ox_sigma**2.))
    nox_sigma=mol_nox*nox_sigma_ratio
    nox_sig_sqr=log(1+nox_sigma**2./mol_nox**2.)
    nox_mu=log(mol_nox**2./sqrt(mol_nox**2.+nox_sigma**2.))

    allocate (y_ox(n_bin))
    allocate (y_nox(n_bin))
    if (.not.allocated(y_Jd)) allocate (y_Jd(n_bin))
    y_Jd=0
    y_ox=0
    y_nox=0
    
    do i=1,n_bin
        y_ox(i)=1./bin(i)/sqrt(ox_sig_sqr)/sqrt(2.*pi)*exp(-((log(bin(i))-ox_mu)**2)/2./ox_sig_sqr)*delta_bin(i)
        y_nox(i)=1./bin(i)/sqrt(nox_sig_sqr)/sqrt(2.*pi)*exp(-((log(bin(i))-nox_mu)**2)/2./nox_sig_sqr)*delta_bin(i)
        !write(*,*) i,bin(i),delta_bin(i),y_ox(i),y_nox(i)
    enddo
    
    !Create the Jd_acc distribution by looping through every hour in the year and extracting the zenith angle
    !Only do this if requested for the first time
    if (run_all) then
        
    y_Jd_acc=0
    do i=1,nbin_Jd
        date_num=1+i/24.
        !call number_to_date(date_num,date_array,2000)
        date_array=0
        zenith_ang=0
        call global_radiation_sub(lat_in,lon_in,date_array,date_num,0.,0.,0.,0.,solar_net,azimuth_ang,zenith_ang)    
        if (zenith_ang.ge.90) then
            y_Jd_acc(i)=0
        else
            y_Jd_acc(i)=((cosd(zenith_ang))**0.28)
        endif
        !write(*,*) i,zenith_ang,y_Jd_acc(i)
    enddo
    
    endif
    
    mean_y_Jd_acc=sum(y_Jd_acc)/nbin_Jd
    do i=1,nbin_Jd
        y_Jd_acc_temp=Jd*y_Jd_acc(i)/mean_y_Jd_acc   
        if (y_Jd_acc_temp/Na_fac*mmass(nox_i)<bin_min) then
            bin_temp=1
        else
            y_Jd_acc_temp_log10=log10(y_Jd_acc_temp/Na_fac*mmass(nox_i))
            bin_temp=floor((y_Jd_acc_temp_log10-log10_bin_min)/delta_log10_bin+0.5)+1
            bin_temp=min(max(bin_temp,1),n_bin)
        endif
        !write(*,'(2i,4es12.2)') i,bin_temp,Jd,y_Jd_acc_temp,y_Jd_acc_temp_log10,log10_bin_min
        y_Jd(bin_temp)= y_Jd(bin_temp)+1
    enddo

    
    !Test
    !y_Jd=0;y_Jd(1)=1
    !y_ox=0
    !y_nox=0
   
    !Normalise all distributions
    y_Jd=y_Jd/sum(y_Jd)
    y_ox=y_ox/sum(y_ox)
    y_nox=y_nox/sum(y_nox)
    
    !do i=1,n_bin
    !    write(*,'(i,5ES12.2)') i,bin(i),delta_bin(i),y_ox(i),y_nox(i),y_Jd(i)
    !enddo

    !allocate (y_all_val)
    !Calculate scaling based on photostationary assumption
    y_all_sum=0
    y_all_prob_sum=0
    
    do j=1,n_bin
    do k=1,n_bin
    do l=1,n_bin
        !y_all_val(j,k,l)=((bin(j)+bin(k)+bin(l)) - sqrt((bin(j)+bin(k)+bin(l))**2 - 4.*bin(j)*bin(l)))/2.
        !y_all_prob(j,k,l)=y_nox(j)*y_Jd(k)*y_ox(l)
        !Calculate weighting
        y_all_prob=y_nox(j)*y_Jd(k)*y_ox(l)
        if (y_all_prob.gt.0) then
        !Calculate NO2 value
        bin_temp2=bin(j)+bin(k)+bin(l)
        y_all_val=(bin_temp2 - sqrt(bin_temp2*bin_temp2 - 4.*bin(j)*bin(l)))/2.
        !Add weighted value
        y_all_sum=y_all_sum+y_all_val*y_all_prob
        !Calculate sum of weights for normalisation later
        y_all_prob_sum=y_all_prob_sum+y_all_prob
        endif
    enddo
    enddo
    enddo
    !y_all_prob=y_all_prob/sum(y_all_prob)
    !y_all=y_all_val*y_all_prob
    !y_all_sum=sum(y_all)
    y_all=y_all_sum/y_all_prob_sum
    
    !Calculate the mean value
    y_annual=((mol_nox+Jd+mol_ox) - sqrt((mol_nox+Jd+mol_ox)**2 - 4.*mol_nox*mol_ox))/2.
    y_scale=y_all/y_annual
    !write(*,*) y_scale,y_all,y_annual
    
    mol_no2_out=y_scale*mol_no2
    mol_o3_out=mol_ox-mol_no2_out
    
    no2_out=mol_no2_out/Na_fac*mmass(no2_i)
    o3_out=mol_o3_out/Na_fac*mmass(o3_i)
    
    deallocate (log10_bin)
    deallocate (bin)
    deallocate (delta_bin)
    deallocate (y_ox)
    deallocate (y_nox)
    deallocate (y_Jd)
    
    end subroutine uEMEP_annual_mean_pdf_correction_NO2_O3