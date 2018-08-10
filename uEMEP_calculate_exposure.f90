!uEMEP_calculate_exposure.f90

    subroutine uEMEP_calculate_exposure
    
    use uEMEP_definitions
    
    implicit none
    
    integer i,j,k
    integer i_source, i_subsource
    integer i_cross,j_cross
    real weighted_concentration(n_source_index)
    real subgrid_area_scaling
    real population_total
    real max_val
    integer i_max,j_max,i_cross_max,j_cross_max
    real val_limit
    real pop_over_limit(subgrid_dim(t_dim_index))
    real grids_over_limit(subgrid_dim(t_dim_index))
    integer t
    
    write(unit_logfile,'(A)') ''
    write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Calculating exposure (uEMEP_calculate_exposure)'
	write(unit_logfile,'(A)') '================================================================'

    !Exposure will be calculated on the target grid
    !The target grid will always be equal to or of a higher resolution
    !Exposure values are given for hourly data as people*concentration/period (period=hour or year)
    !Population weighted concentrations are also calculated
    
    !Loop through target grid and find the population
    exposure_subgrid=0.
    subgrid_area_scaling=1.
    !Scale population due to difference in grid size sizes only when population grid is used
    !Only works when the target subgrid is smaller than the population subgrid
    if (population_index.eq.2) then      
        !This does not work when target grid is larger than the population grid
        subgrid_area_scaling=(subgrid_delta(x_dim_index)*subgrid_delta(y_dim_index))/(population_subgrid_delta(x_dim_index)*population_subgrid_delta(y_dim_index))
    endif
    
    population_total=0
    
    max_val=-1.
    pop_over_limit=0.
    grids_over_limit=0.
    val_limit=40.
    i_max=0;j_max=0
    
    do j=1,subgrid_dim(y_dim_index)
    do i=1,subgrid_dim(x_dim_index)

        i_cross=crossreference_target_to_population_subgrid(i,j,x_dim_index)
        j_cross=crossreference_target_to_population_subgrid(i,j,y_dim_index)
            
        do i_source=1,n_source_index
        if (calculate_source(i_source).and.use_subgrid(i,j,i_source)) then
            exposure_subgrid(i,j,:,i_source)=sum(subgrid(i,j,:,local_subgrid_index,i_source,:),2)*population_subgrid(i_cross,j_cross,population_data_type)*subgrid_area_scaling         
        endif
        enddo

        !Calculate number over limit. Does nto work for hourly data properly yet
        i_source=allsource_index
        i_subsource=1
        do t=1,subgrid_dim(t_dim_index)
            
        if (use_subgrid(i,j,i_source)) then
             !write(*,*) i,j,t,comp_subgrid(i,j,t,no2_nc_index)
             if (comp_subgrid(i,j,t,no2_nc_index).gt.val_limit) then
                pop_over_limit(t)=pop_over_limit(t)+population_subgrid(i_cross,j_cross,population_data_type)*subgrid_area_scaling
                grids_over_limit(t)=grids_over_limit(t)+1
                !write(*,*) pop_over_limit(t)
             endif
             
        endif
        enddo
        

        if (use_subgrid(i,j,allsource_index)) then
            exposure_subgrid(i,j,:,allsource_index)=sum(subgrid(i,j,:,total_subgrid_index,allsource_index,:),2)*population_subgrid(i_cross,j_cross,population_data_type)*subgrid_area_scaling
            population_total=population_total+population_subgrid(i_cross,j_cross,population_data_type)*subgrid_area_scaling
            if (sum(subgrid(i,j,:,total_subgrid_index,allsource_index,:))/subgrid_dim(t_dim_index).gt.max_val) then
                max_val=sum(subgrid(i,j,:,total_subgrid_index,allsource_index,:))/subgrid_dim(t_dim_index)
                i_max=i;j_max=j
            endif
            if (sum(subgrid(i,j,:,total_subgrid_index,allsource_index,:))/subgrid_dim(t_dim_index).gt.max_val) then
                max_val=sum(subgrid(i,j,:,total_subgrid_index,allsource_index,:))/subgrid_dim(t_dim_index)
                i_max=i;j_max=j
            endif
        endif
        
    enddo
    enddo

    i_cross_max=crossreference_target_to_population_subgrid(i_max,j_max,x_dim_index)
    j_cross_max=crossreference_target_to_population_subgrid(i_max,j_max,y_dim_index)

    write(unit_logfile,'(A)') 'Population weighted concentration and max time average concentration (with population) by source per period for '//trim(input_comp_name)
    weighted_concentration(allsource_index)=sum(exposure_subgrid(:,:,:,allsource_index))/population_total/subgrid_dim(t_dim_index)
    write(unit_logfile,'(A24,2f12.2,2f12.0)') 'Total: ',weighted_concentration(allsource_index),max_val,population_total,population_subgrid(i_cross_max,j_cross_max,population_data_type)
    
    !Calculate population weighted values for each source
        do i_source=1,n_source_index
        if (calculate_source(i_source)) then
            weighted_concentration(i_source)=sum(exposure_subgrid(:,:,:,i_source))/population_total/subgrid_dim(t_dim_index)
            write(unit_logfile,'(A24,2f12.2)') trim(source_file_str(i_source))//': ',weighted_concentration(i_source),sum(subgrid(i_max,j_max,:,local_subgrid_index,i_source,:))/subgrid_dim(t_dim_index)
        endif
        enddo
    
        write(unit_logfile,'(A24,2f12.2)') 'nonlocal: ',weighted_concentration(allsource_index)-sum(weighted_concentration(:))+weighted_concentration(allsource_index) &
            ,2*sum(subgrid(i_max,j_max,:,total_subgrid_index,allsource_index,:))/subgrid_dim(t_dim_index)-sum(subgrid(i_max,j_max,:,local_subgrid_index,:,:))/subgrid_dim(t_dim_index)

    !In case of no2 recalculate the total and present it again as no2 and not nox
    if (compound_index.eq.nox_nc_index) then
    do j=1,subgrid_dim(y_dim_index)
    do i=1,subgrid_dim(x_dim_index)
        
        i_cross=crossreference_target_to_population_subgrid(i,j,x_dim_index)
        j_cross=crossreference_target_to_population_subgrid(i,j,y_dim_index)
            
        if (use_subgrid(i,j,allsource_index)) then
            exposure_subgrid(i,j,:,allsource_index)=comp_subgrid(i,j,:,no2_nc_index)*population_subgrid(i_cross,j_cross,population_data_type)*subgrid_area_scaling
        endif
        
    enddo
    enddo
       
    weighted_concentration(allsource_index)=sum(exposure_subgrid(:,:,:,allsource_index))/population_total/subgrid_dim(t_dim_index)
    write(unit_logfile,'(A24,2f12.2)') 'Total no2: ',weighted_concentration(allsource_index),sum(comp_subgrid(i_max,j_max,:,no2_nc_index))/subgrid_dim(t_dim_index)
    write(unit_logfile,'(A24,f12.2)') 'Population over limit: ',maxval(pop_over_limit)
    write(unit_logfile,'(A24,f12.2)') 'Grids over limit: ',maxval(grids_over_limit)
    
    endif
    
    
    end subroutine uEMEP_calculate_exposure
