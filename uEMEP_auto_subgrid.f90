!uEMEP_auto_subgrid.f90
    
!==========================================================================
!   uEMEP model uEMEP_auto_subgrid
!   Automatically creates a grid dependent on the distance to source
!   
!==========================================================================
    subroutine uEMEP_auto_subgrid
    
    use uEMEP_definitions

    implicit none
    
    integer i,j,k
    integer         i_source
    integer         t
    character(256)  temp_name
    logical         exists
    integer, allocatable :: use_subgrid_val(:,:,:)    
    
    real :: max_use_subgrid_size(n_source_index)
    real :: use_subgrid_range(n_source_index)
    integer n_use_subgrid_levels(n_source_index)
    integer use_subgrid_step
    integer use_emission_subgrid_space
    integer i_cross,j_cross
    integer i_start,i_end,j_start,j_end
    real sum_emission
    
        
    !Exit if emission positions are not to be used and if one of the other auto positions is to be used
    if (.not.use_emission_positions_for_auto_subgrid_flag.and.(use_population_positions_for_auto_subgrid_flag.or.use_receptor_positions_for_auto_subgrid_flag)) then
        return
    endif
    !Exit and fill subgrid if none of the auto subgrids are to be used
    if (.not.use_emission_positions_for_auto_subgrid_flag.and..not.use_population_positions_for_auto_subgrid_flag.and..not.use_receptor_positions_for_auto_subgrid_flag) then
        use_subgrid=.true.
        return
    endif
    

    allocate(use_subgrid_val(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),n_source_index))
    
    write(unit_logfile,'(A)') ''
    write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Setting subgrids automatically (uEMEP_auto_subgrid)'
	write(unit_logfile,'(A)') '================================================================'
   
    !Read in use subgrid files
    do i_source=1,n_source_index
    if (calculate_source(i_source)) then

    if (read_existing_grid_data(use_subgrid_file_index(i_source))) then
        if (calculate_source(i_source)) then
            temp_name=trim(pathname_grid(use_subgrid_file_index(i_source)))//trim(filename_grid(use_subgrid_file_index(i_source)))//'_'//trim(file_tag)//'.asc'
            inquire(file=temp_name,exist=exists)
            if (.not.exists) then
                write(unit_logfile,*)'ERROR: '//trim(temp_name)//' does not exist.'
                return
            endif
            call read_esri_ascii_3d_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),1,subgrid_delta(1), &
                real(use_subgrid_val(:,:,i_source)),x_subgrid,y_subgrid)
            
            !Convert values to logical
            do j=1,subgrid_dim(y_dim_index)
            do i=1,subgrid_dim(x_dim_index)
                if (use_subgrid_val(i,j,i_source).gt.0) then
                    use_subgrid(i,j,i_source)=.true.
                else
                    use_subgrid(i,j,i_source)=.true.
                endif      
            enddo
            enddo
            
        return
        endif
    endif
    endif
    enddo
    
    
    !Set all subgrids to do not use
    use_subgrid_val=0
    !Set time index used for emissions to 1
    t=1
    !Set the maximum grid size
    max_use_subgrid_size=2500.
    max_use_subgrid_size(traffic_index)=2500.
    max_use_subgrid_size(shipping_index)=2500.
    use_subgrid_range=8.
    use_subgrid_range(traffic_index)=6.
    use_subgrid_range(shipping_index)=12.
    


    !Set the number of levels to match this
    n_use_subgrid_levels=floor(log(max_use_subgrid_size/sqrt(subgrid_delta(x_dim_index)*subgrid_delta(y_dim_index)))/log(2.)+.5)
        

    do i_source=1,n_source_index
    if (calculate_source(i_source)) then
        write(unit_logfile,'(a,2f10.1,i10)') trim(source_file_str(i_source))//': maximum use subgrid size (m), grid range and number of levels: ',max_use_subgrid_size(i_source),use_subgrid_range(i_source),n_use_subgrid_levels(i_source)
        do k=0,n_use_subgrid_levels(i_source)
            use_subgrid_step=2**k
            use_emission_subgrid_space=floor(use_subgrid_step/sqrt(emission_subgrid_delta(x_dim_index,i_source)*emission_subgrid_delta(y_dim_index,i_source)/subgrid_delta(x_dim_index)/subgrid_delta(y_dim_index))/2.*use_subgrid_range(i_source))
            !use_emission_subgrid_space=floor(use_subgrid_step/2.*use_subgrid_range(i_source))
            use_emission_subgrid_space=max(1,use_emission_subgrid_space)
            !write(*,*) k,use_subgrid_step,use_emission_subgrid_space,use_subgrid_step*subgrid_delta(y_dim_index),use_emission_subgrid_space*emission_subgrid_delta(y_dim_index,i_source)
            
            do j=use_subgrid_step,subgrid_dim(y_dim_index),use_subgrid_step
            do i=use_subgrid_step,subgrid_dim(x_dim_index),use_subgrid_step
                
                i_cross=crossreference_target_to_emission_subgrid(i,j,x_dim_index,i_source)
                j_cross=crossreference_target_to_emission_subgrid(i,j,y_dim_index,i_source)

                !If the target grid is inside an emission grid with non zero emissions then include it
                !if (sum(emission_subgrid(i_cross,j_cross,t,i_source,:)).gt.0) use_subgrid_val(i,j,i_source)=1
                
                !Search in the neighbourhood
                i_start=max(1,i_cross-use_emission_subgrid_space)
                i_end=min(emission_subgrid_dim(x_dim_index,i_source),i_cross+use_emission_subgrid_space)
                j_start=max(1,j_cross-use_emission_subgrid_space)
                j_end=min(emission_subgrid_dim(y_dim_index,i_source),j_cross+use_emission_subgrid_space)

                sum_emission=sum(proxy_emission_subgrid(i_start:i_end,j_start:j_end,i_source,:))
                if (sum_emission.gt.0.) use_subgrid_val(i,j,i_source)=1
                !write(*,'(6i,f)') i,j,i_start,i_end,j_start,j_end,sum_emission
                
                !Always include the coarsest level everywhere
                if (k.eq.n_use_subgrid_levels(i_source)) use_subgrid_val(i,j,i_source)=1
                
            enddo
            enddo
        enddo

            !Put points around the domain edge at lowest resolution for proper interpolation
            k=n_use_subgrid_levels(i_source)
            use_subgrid_step=2**k
            do j=1,subgrid_dim(y_dim_index),use_subgrid_step
                use_subgrid_val(1,j,i_source)=1
                use_subgrid_val(subgrid_dim(x_dim_index),j,i_source)=1
            enddo
            do i=1,subgrid_dim(x_dim_index),use_subgrid_step
                use_subgrid_val(i,1,i_source)=1
                use_subgrid_val(i,subgrid_dim(y_dim_index),i_source)=1
            enddo            
            use_subgrid_val(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),i_source)=1

        
        use_subgrid_val(:,:,allsource_index)=use_subgrid_val(:,:,allsource_index)+use_subgrid_val(:,:,i_source)

    endif
    enddo
    
    use_subgrid_val=min(1,use_subgrid_val)
 
    
    !Convert values to logical
    do i_source=1,n_source_index 
    if (calculate_source(i_source).or.i_source.eq.allsource_index) then
        do j=1,subgrid_dim(y_dim_index)
        do i=1,subgrid_dim(x_dim_index)
            if (use_subgrid_val(i,j,i_source).gt.0) then
                use_subgrid(i,j,i_source)=.true.
            else
                use_subgrid(i,j,i_source)=.false.
            endif      
        enddo
        enddo
    endif
    enddo
    
    !If the grids are not to be interpolated then they must all be calculated at the same place for the different sources for interpolation later
    if (.not.interpolate_subgrids_flag) then
    do i_source=1,n_source_index 
        use_subgrid(:,:,i_source)=use_subgrid(:,:,allsource_index)
    enddo
    endif

    do i_source=1,n_source_index 
        if (.not.read_existing_grid_data(use_subgrid_file_index(i_source))) then
            if (calculate_source(i_source).or.i_source.eq.allsource_index) then
                write(unit_logfile,'(a,2i10,f6.1)') 'Number of calculation subgrids for '//trim(source_file_str(i_source))//' (number, total, percent):',sum(use_subgrid_val(:,:,i_source)),subgrid_dim(1)*subgrid_dim(2),sum(use_subgrid_val(:,:,i_source))*100./(subgrid_dim(1)*subgrid_dim(2))
            endif
        endif
    enddo

    if (save_intermediate_files) then
    do i_source=1,n_source_index 
        if (.not.read_existing_grid_data(use_subgrid_file_index(i_source))) then
            if (calculate_source(i_source).or.i_source.eq.allsource_index) then
                !write(unit_logfile,'(a,2i10,f6.1)') 'Number of calculation subgrids for '//trim(source_file_str(i_source))//' (number, total, fraction):',sum(use_subgrid_val(:,:,i_source)),subgrid_dim(1)*subgrid_dim(2),sum(use_subgrid_val(:,:,i_source))*1./(subgrid_dim(1)*subgrid_dim(2))
                temp_name=trim(pathname_grid(use_subgrid_file_index(i_source)))//trim(filename_grid(use_subgrid_file_index(i_source)))//'_'//trim(file_tag)//'.asc'
                write(unit_logfile,'(a)')'Writing to: '//trim(temp_name)
                call write_esri_ascii_3d_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),1,subgrid_delta(1), &
                     real(use_subgrid_val(:,:,i_source)),x_subgrid,y_subgrid)
            
            endif
        endif
    enddo
    endif
    
    end subroutine uEMEP_auto_subgrid
    