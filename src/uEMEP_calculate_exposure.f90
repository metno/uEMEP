module calculate_exposure

    use uemep_configuration
    use uEMEP_definitions

    implicit none
    private

    public :: uEMEP_calculate_exposure

contains

    subroutine uEMEP_calculate_exposure()
        integer :: i, j
        integer :: i_source
        integer :: i_cross, j_cross
        real :: weighted_concentration(n_source_index,n_pollutant_loop)
        real :: subgrid_area_scaling
        real :: population_total(n_pollutant_loop)
        real :: max_val(n_pollutant_loop)
        integer :: i_max(n_pollutant_loop), j_max(n_pollutant_loop), i_cross_max(n_pollutant_loop), j_cross_max(n_pollutant_loop)
        real :: val_limit(n_compound_nc_index)
        real :: pop_over_limit(subgrid_dim(t_dim_index),n_pollutant_loop)
        real :: grids_over_limit(subgrid_dim(t_dim_index),n_pollutant_loop)
        integer :: t
        integer :: i_pollutant

        write(unit_logfile,'(A)') ''
        write(unit_logfile,'(A)') '================================================================'
        write(unit_logfile,'(A)') 'Calculating exposure (uEMEP_calculate_exposure)'
        write(unit_logfile,'(A)') '================================================================'

        ! Loop through target grid and find the population
        exposure_subgrid = 0.0
        subgrid_area_scaling = 1.0
        
        ! Scale population due to difference in grid size sizes only when population grid is used
        ! Only works when the target subgrid is smaller than the population subgrid
        if (population_index .eq. 2) then
            ! This does not work when target grid is larger than the population grid
            subgrid_area_scaling = (subgrid_delta(x_dim_index)*subgrid_delta(y_dim_index))/(population_subgrid_delta(x_dim_index)*population_subgrid_delta(y_dim_index))
        end if

        population_total = 0
        max_val = -1.0
        pop_over_limit = 0.0
        grids_over_limit = 0.0
        val_limit(no2_index) = 100.0 ! TODO: Need to fix this later
        val_limit(pm10_index) = 50.0 ! TODO: Need to fix this later
        val_limit(pm25_index)=30.0 ! TODO: Need to fix this later
        i_max = 0; j_max = 0

        do j = 1, subgrid_dim(y_dim_index)
            do i = 1, subgrid_dim(x_dim_index)
                i_cross = crossreference_target_to_population_subgrid(i,j,x_dim_index)
                j_cross = crossreference_target_to_population_subgrid(i,j,y_dim_index)
                do i_source = 1, n_source_index
                    if (calculate_source(i_source) .and. use_subgrid(i,j,i_source)) then
                        exposure_subgrid(i,j,:,i_source,:) = subgrid(i,j,:,local_subgrid_index,i_source,:)*population_subgrid(i_cross,j_cross,population_data_type)*subgrid_area_scaling
                    end if
                end do

                ! Calculate number over limit. Does not work for hourly data properly yet
                i_source = allsource_index
                do i_pollutant = 1, n_pollutant_loop
                    do t = 1, subgrid_dim(t_dim_index)
                        if (use_subgrid(i,j,i_source)) then
                            if (pollutant_loop_index(i_pollutant) .eq. no2_nc_index) then
                                if (comp_subgrid(i,j,t,pollutant_loop_index(i_pollutant)) .gt. val_limit(pollutant_loop_index(i_pollutant))) then
                                    pop_over_limit(t,i_pollutant) = pop_over_limit(t,i_pollutant) + population_subgrid(i_cross,j_cross,population_data_type)*subgrid_area_scaling
                                    grids_over_limit(t,i_pollutant) = grids_over_limit(t,i_pollutant) + 1
                                end if
                            end if
                        end if
                    end do

                    if (use_subgrid(i,j,allsource_index)) then
                        exposure_subgrid(i,j,:,allsource_index,i_pollutant) = subgrid(i,j,:,total_subgrid_index,allsource_index,i_pollutant)*population_subgrid(i_cross,j_cross,population_data_type)*subgrid_area_scaling
                        population_total(i_pollutant) = population_total(i_pollutant) + population_subgrid(i_cross,j_cross,population_data_type)*subgrid_area_scaling
                        if (sum(subgrid(i,j,:,total_subgrid_index,allsource_index,i_pollutant))/subgrid_dim(t_dim_index) .gt. max_val(i_pollutant)) then
                            max_val(i_pollutant) = sum(subgrid(i,j,:,total_subgrid_index,allsource_index,i_pollutant))/subgrid_dim(t_dim_index)
                            i_max(i_pollutant) = i
                            j_max(i_pollutant) = j
                        end if
                        i_cross_max(i_pollutant) = crossreference_target_to_population_subgrid(i_max(i_pollutant),j_max(i_pollutant),x_dim_index)
                        j_cross_max(i_pollutant) = crossreference_target_to_population_subgrid(i_max(i_pollutant),j_max(i_pollutant),y_dim_index)
                    end if
                end do
            end do
        end do


        do i_pollutant = 1, n_pollutant_loop
            write(unit_logfile,'(A)') 'Population weighted concentration and max time average concentration (with population) by source per period for '//trim(input_comp_name)
            weighted_concentration(allsource_index,i_pollutant) = sum(exposure_subgrid(:,:,:,allsource_index,i_pollutant))/population_total(i_pollutant)/subgrid_dim(t_dim_index)
            write(unit_logfile,'(A24,2f12.2,2f12.0)') 'Total '//trim(pollutant_file_str(pollutant_loop_index(i_pollutant)))//': ', weighted_concentration(allsource_index,i_pollutant), max_val(i_pollutant), population_total(i_pollutant), population_subgrid(i_cross_max(i_pollutant), j_cross_max(i_pollutant), population_data_type)

            ! Calculate population weighted values for each source
            do i_source = 1, n_source_index
                if (calculate_source(i_source)) then
                    weighted_concentration(i_source,i_pollutant) = sum(exposure_subgrid(:,:,:,i_source,i_pollutant))/population_total(i_pollutant)/subgrid_dim(t_dim_index)
                    write(unit_logfile,'(A24,2f12.2)') trim(source_file_str(i_source))//': ',weighted_concentration(i_source,i_pollutant),sum(subgrid(i_max(i_pollutant),j_max(i_pollutant),:,local_subgrid_index,i_source,i_pollutant))/subgrid_dim(t_dim_index)
                end if
            end do

            write(unit_logfile,'(A24,2f12.2)') 'nonlocal: ', weighted_concentration(allsource_index,i_pollutant)-sum(weighted_concentration(:,i_pollutant))+weighted_concentration(allsource_index,i_pollutant), &
                2*sum(subgrid(i_max(i_pollutant),j_max(i_pollutant),:,total_subgrid_index,allsource_index,i_pollutant))/subgrid_dim(t_dim_index)-sum(subgrid(i_max(i_pollutant),j_max(i_pollutant),:,local_subgrid_index,:,i_pollutant))/subgrid_dim(t_dim_index)

            ! In case of no2 recalculate the total and present it again as no2 and not nox
            if (pollutant_loop_index(i_pollutant) .eq. nox_nc_index) then
                do j = 1, subgrid_dim(y_dim_index)
                    do i = 1, subgrid_dim(x_dim_index)
                        i_cross = crossreference_target_to_population_subgrid(i,j,x_dim_index)
                        j_cross = crossreference_target_to_population_subgrid(i,j,y_dim_index)
                        if (use_subgrid(i,j,allsource_index)) then
                            exposure_subgrid(i,j,:,allsource_index,i_pollutant) = comp_subgrid(i,j,:,no2_nc_index)*population_subgrid(i_cross,j_cross,population_data_type)*subgrid_area_scaling
                        end if
                    end do
                end do
                weighted_concentration(allsource_index,i_pollutant) = sum(exposure_subgrid(:,:,:,allsource_index,i_pollutant))/population_total(i_pollutant)/subgrid_dim(t_dim_index)
                write(unit_logfile,'(A24,2f12.2)') 'Total no2: ', weighted_concentration(allsource_index,i_pollutant), sum(comp_subgrid(i_max(i_pollutant),j_max(i_pollutant),:,no2_nc_index))/subgrid_dim(t_dim_index)
                write(unit_logfile,'(A24,f12.2)') 'Population over limit: ', maxval(pop_over_limit)
                write(unit_logfile,'(A24,f12.2)') 'Grids over limit: ', maxval(grids_over_limit)
            end if
        end do !pollutant loop

    end subroutine uEMEP_calculate_exposure

end module calculate_exposure

