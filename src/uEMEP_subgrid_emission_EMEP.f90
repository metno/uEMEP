module subgrid_emission_emep

    !! Distribute EMEP emissions onto uEMEP proxy and local subgrids
    !!
    !! Copyright (C) 2007 Free Software Foundation.
    !! License GNU LGPL-3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>.
    !! This is free software: you are free to change and redistribute it.
    !!
    !! Developed and maintained at the Norwegian Meteorological Institute.
    !! Contribute at: <https://github.com/metno/uEMEP>
    !!
    !! This module contains procedures to map EMEP grid emissions to proxy and local subgrids
    !!
    !! Supported algorithms are:
    !!   - distribute evenly
    !!   - distribute using area weighting
    !!   - distribute using emission weighting
    !!
    !! The procedures operate on global level arrays (e.g., emission_subgrid, proxy_emission_subgrid, var3d_nc) and
    !! apply unit conversion from EMEP source units (mg/m2/h or annual aggregates) to uEMEP units (ug/s/subgrid)

    use uemep_definitions
    use uemep_configuration
    use mod_lambert_projection, only: lb2lambert2_uEMEP, LL2PS_spherical, lb2lambert_uEMEP

    implicit none
    private

    public :: uEMEP_subgrid_emission_EMEP

contains

    subroutine uEMEP_subgrid_emission_EMEP()
        !! High-level driver for distributing EMEP emissions to subgrids
        !!
        !! Prepares working arrays, selects the distribution method according to configuration flags, 
        !! calls the lower-level distribution routines and applies optional GNFR scaling to emission_subgrid

        integer :: t
        real, allocatable :: weighting_nc(:,:), weighting_subgrid(:,:,:)
        real, allocatable :: total_weighting_nc(:,:,:), proxy_weighting_nc(:,:,:)
        real, allocatable :: total_proxy_emission_subgrid(:,:,:,:)
        real, allocatable :: total_proxy_subgrid_emission_in_EMEP_grid(:,:,:,:)
        integer, allocatable :: subgrid_count_nc(:,:)
        integer, allocatable :: subgrid_count_subgrid(:,:,:)
        integer :: t_start, t_end
        integer :: i_source, i_pollutant

        if (local_subgrid_method_flag /= 3 .and. local_subgrid_method_flag /= 4) then
            return
        end if

        write(unit_logfile,'(A)') ''
        write(unit_logfile,'(A)') '================================================================'
        write(unit_logfile,'(A)') 'Distributing EMEP emission to subgrids (uEMEP_subgrid_emission_EMEP)'
        write(unit_logfile,'(A)') '================================================================'

        ! Allocate and save the existing emission subgrid data (??)
        allocate (total_proxy_emission_subgrid(emission_max_subgrid_dim(x_dim_index), &
            emission_max_subgrid_dim(y_dim_index),n_source_index,n_pollutant_loop))

        ! Set the start and end times of the loop
        t_start = 1
        t_end = subgrid_dim(t_dim_index)

        ! Distribute the EMEP emissions evenly over the subgrids within an EMEP grid
        if (EMEP_emission_grid_interpolation_flag == 0 .or. local_subgrid_method_flag == 4) then
            write(unit_logfile,'(A)') 'Distributing EMEP emissions to all subgrids within an EMEP grid'
            call distribute_emep_evenly(t_start, t_end, total_proxy_emission_subgrid, &
                total_proxy_subgrid_emission_in_EMEP_grid, subgrid_count_subgrid, subgrid_count_nc)
        end if

        ! Distribute the EMEP emissions to proxy emission subgrids using area weighting of the EMEP grid
        ! This is done also if there is moving window weighting later as it is used for the nonlocal contribution
        if (EMEP_emission_grid_interpolation_flag == 1 .and. local_subgrid_method_flag /=4) then
            write(unit_logfile,'(A)') 'Distributing EMEP emissions to proxy emission subgrids using area weighting of the EMEP grid'
            call distribute_emep_to_proxy_by_area_weighting( &
                t_start, t_end, weighting_nc, subgrid_count_subgrid, total_proxy_emission_subgrid)
        end if

        ! Distribute EMEP emissions to proxy subgrids using emission weighting of the EMEP grid
        if (EMEP_emission_grid_interpolation_flag == 2 .and. local_subgrid_method_flag /= 4) then
            write(unit_logfile,'(A)') 'Distributing EMEP emissions to proxy subgrids using emission weighting of the EMEP grid'
            call distribute_emep_to_proxy_by_emission_weighting(t_start, t_end, total_weighting_nc, proxy_weighting_nc, &
                subgrid_count_subgrid, total_proxy_emission_subgrid, weighting_subgrid)
        endif

        ! Scale the subgrid emissions if GNFR is used
        do i_source = 1, n_source_index
            if (calculate_source(i_source)) then
                if (scale_GNFR_emission_source(i_source) /= 1.0) then
                    write(unit_logfile,'(2A,f12.2)') 'Scaling EMEP emissions for uEMEP source: ', &
                        trim(source_file_str(i_source)), scale_GNFR_emission_source(i_source)
                    do t = t_start, t_end
                        emission_subgrid(:,:,t,i_source,:) = emission_subgrid(:,:,t,i_source,:)*scale_GNFR_emission_source(i_source)
                    end do
                end if
            end if
        end do

        do i_source = 1, n_source_index
            if (calculate_source(i_source)) then
                do i_pollutant = 1, n_pollutant_loop
                    write(unit_logfile,'(3A,ES10.2)') 'Emission source ', &
                        trim(source_file_str(i_source))//' '//trim(pollutant_file_str(pollutant_loop_index(i_pollutant))), &
                        ': Total hourly average emissions after use of EMEP (ug/s)=', &
                        sum(emission_subgrid(1:emission_subgrid_dim(x_dim_index,i_source), &
                            1:emission_subgrid_dim(y_dim_index,i_source),:,i_source,i_pollutant))/(t_end - t_start + 1)
                end do
            end if
        end do

        if (allocated(weighting_nc)) deallocate(weighting_nc)
        if (allocated(total_weighting_nc)) deallocate(total_weighting_nc)
        if (allocated(proxy_weighting_nc)) deallocate(proxy_weighting_nc)
        if (allocated(weighting_subgrid)) deallocate(weighting_subgrid)
        if (allocated(total_proxy_emission_subgrid)) deallocate(total_proxy_emission_subgrid)
        if (allocated(total_proxy_subgrid_emission_in_EMEP_grid)) deallocate(total_proxy_subgrid_emission_in_EMEP_grid)
        if (allocated(subgrid_count_nc)) deallocate(subgrid_count_nc)
        if (allocated(subgrid_count_subgrid)) deallocate(subgrid_count_subgrid)

    end subroutine uEMEP_subgrid_emission_EMEP

    subroutine distribute_emep_evenly(t_start, t_end, tot_proxy_emis_sgrid, tot_proxy_sgrid_emis_in_emep_grid, &
        sgrid_nb_sgrid, sgrid_nb_nc)
        !! Evenly distribute EMEP grid emissions to subgrids inside each EMEP cell
        !!
        !! Fills emission_subgrid with EMEP values and, optionally, applies distribution to existing proxy_emission_subgrid
        !! Converts units from mg/m2/hour (or annual mg/m2/year) to ug/s/subgrid
        integer, intent(in) :: t_start
        integer, intent(in) :: t_end
        real, intent(inout) :: tot_proxy_emis_sgrid(:,:,:,:)
        real, allocatable, intent(inout) :: tot_proxy_sgrid_emis_in_emep_grid(:,:,:,:)
        integer, allocatable, intent(inout) :: sgrid_nb_sgrid(:,:,:)
        integer, allocatable, intent(inout) :: sgrid_nb_nc(:,:)

        integer :: i, j, ii, jj
        integer :: tt, t
        integer :: i_source, i_pollutant
        real :: sum_temp(n_pollutant_loop)

        tt = 1
        allocate(tot_proxy_sgrid_emis_in_emep_grid(dim_length_nc(x_dim_nc_index), &
            dim_length_nc(y_dim_nc_index),n_source_index,n_pollutant_loop))
        allocate(sgrid_nb_sgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),n_pollutant_loop))
        allocate(sgrid_nb_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index)))

        tot_proxy_sgrid_emis_in_emep_grid = 0.0
        tot_proxy_emis_sgrid = 0.0

        do i_source = 1, n_source_index
            if (calculate_source(i_source)) then

                emission_subgrid(:,:,:,i_source,:) = 0.0
                sgrid_nb_sgrid = 0
                sgrid_nb_nc = 0

                ! Calculate total subgrid emissions and number of subgrids in each EMEP grid
                do j = 1, emission_subgrid_dim(y_dim_index,i_source)
                    do i = 1, emission_subgrid_dim(x_dim_index,i_source)
                        ii = crossreference_emission_to_emep_subgrid(i,j,x_dim_index,i_source)
                        jj = crossreference_emission_to_emep_subgrid(i,j,y_dim_index,i_source)
                        ii = max(min(ii,dim_length_nc(x_dim_nc_index)),1)
                        jj = max(min(jj,dim_length_nc(y_dim_nc_index)),1)
                        sgrid_nb_nc(ii,jj) = sgrid_nb_nc(ii,jj) + 1
                        tot_proxy_sgrid_emis_in_emep_grid(ii,jj,i_source,:) = &
                            tot_proxy_sgrid_emis_in_emep_grid(ii,jj,i_source,:)+proxy_emission_subgrid(i,j,i_source,:)
                        emission_subgrid(i,j,:,i_source,:) = var3d_nc(ii,jj,:,emis_nc_index,i_source,:)
                    end do
                end do

                ! Transfer the total emissions in the EMEP grid to each subgrid within that grid
                do j = 1, emission_subgrid_dim(y_dim_index,i_source)
                    do i = 1, emission_subgrid_dim(x_dim_index,i_source)

                        ii = crossreference_emission_to_emep_subgrid(i,j,x_dim_index,i_source)
                        jj = crossreference_emission_to_emep_subgrid(i,j,y_dim_index,i_source)
                        ii = max(min(ii,dim_length_nc(x_dim_nc_index)),1)
                        jj = max(min(jj,dim_length_nc(y_dim_nc_index)),1)

                        tot_proxy_emis_sgrid(i,j,i_source,:) = tot_proxy_sgrid_emis_in_emep_grid(ii,jj,i_source,:)
                        sgrid_nb_sgrid(i,j,:)=sgrid_nb_nc(ii,jj)

                        ! Converts from mg/m2/hour(year) to ug/s/subgrid assuming the orig EMEP emissions are in mg/m2/hour(year)
                        if (hourly_calculations) then
                            emission_subgrid(i,j,:,i_source,:) = emission_subgrid(i,j,:,i_source,:) &
                                *emission_subgrid_delta(x_dim_index,i_source)*emission_subgrid_delta(y_dim_index,i_source) &
                                *1000.0/3600.0
                        end if
                        if (annual_calculations) then
                            emission_subgrid(i,j,:,i_source,:)=emission_subgrid(i,j,:,i_source,:) &
                                *emission_subgrid_delta(x_dim_index,i_source)*emission_subgrid_delta(y_dim_index,i_source) &
                                *1000.0/3600.0/EMEP_emission_aggregation_period
                        end if
                    end do
                end do

                ! Determine subgrid normalised time profile per hour from EMEP grid emissions (average hourly emission conversion)
                ! This is not quite right because the entire emission time profile is not available in a short period
                ! Minimum of a day is needed, with the assumption that all days are the same
                if (local_subgrid_method_flag == 4) then
                    write(unit_logfile,'(A)') 'Calculating EMEP emission time profile'
                    do j = 1, emission_subgrid_dim(y_dim_index,i_source)
                        do i = 1, emission_subgrid_dim(x_dim_index,i_source)
                            if (hourly_calculations) then
                                sum_temp(:) = sum(emission_subgrid(i,j,:,i_source,:),1)
                                do i_pollutant = 1, n_pollutant_loop
                                    emission_time_profile_subgrid(i,j,:,i_source,i_pollutant) = &
                                        emission_subgrid(i,j,:,i_source,i_pollutant)/sum_temp(i_pollutant) &
                                            *emission_subgrid_dim(t_dim_index,i_source)
                                    if (sum_temp(i_pollutant) == 0.0) then
                                        emission_time_profile_subgrid(i,j,:,i_source,i_pollutant) = 0.0
                                    end if
                                end do
                            else
                                emission_time_profile_subgrid(i,j,:,i_source,:) = 1.0
                            endif

                            ! Set emissions to 0 in the case when local_subgrid_method_flag.eq.4 since these are set later
                            ! This way of doing things is not logical as it fills the grid unnecessarilly. 
                            ! Should be fixed and made logical
                            emission_subgrid(i,j,:,i_source,:) = 0.0
                        end do
                    end do
                end if

                ! Distribute EMEP emissions to existing proxy subgrid emissions
                if (subgrid_emission_distribution_flag .and. local_subgrid_method_flag /= 4) then
                    if (local_subgrid_method_flag == 2) then
                        write(unit_logfile,'(2A)') 'Distributing local emission data to subgrid emissions for: ', &
                            trim(source_file_str(i_source))
                    end if
                    if (local_subgrid_method_flag == 3) then
                        write(unit_logfile,'(2A)') &
                            'Distributing EMEP emissions to proxy subgrid emissions, no weighting used, for: ', &
                            trim(source_file_str(i_source))
                    end if

                    do t = t_start, t_end
                        emission_subgrid(:,:,t,i_source,:) = emission_subgrid(:,:,t,i_source,:)*sgrid_nb_sgrid(:,:,:) &
                            *proxy_emission_subgrid(:,:,i_source,:)/tot_proxy_emis_sgrid(:,:,i_source,:)
                        where (tot_proxy_emis_sgrid(:,:,i_source,:) == 0.0) emission_subgrid(:,:,t,i_source,:) = 0.0
                    end do

                    ! Fix round off negatives
                    where (emission_subgrid(:,:,:,i_source,:) < 0) emission_subgrid(:,:,:,i_source,:) = 0.0
                endif
            endif
        enddo
    end subroutine distribute_emep_evenly

    subroutine distribute_emep_to_proxy_by_area_weighting(t_start, t_end, weighting_nc, subgrid_count_subgrid, &
        total_proxy_emission_subgrid)
        !! Area-weighted mapping of EMEP grid emissions to proxy subgrids
        !!
        !! Notes: 
        !!   - For each proxy subgrid cell the routine computes overlap area with nearby
        !!     EMEP cells (3×3 window) and accumulates weighted EMEP emissions
        !!   - Converts units to ug/s/subgrid
        !!   - Currently assumes neighbour indices exist
        !!   - Quick calculation of area weighting, no edge effects. Does not need to change with time
        integer, intent(in) :: t_start
        integer, intent(in) :: t_end
        real, allocatable, intent(inout) :: weighting_nc(:,:)
        integer, allocatable, intent(inout) :: subgrid_count_subgrid(:,:,:)
        real, allocatable, intent(inout) :: total_proxy_emission_subgrid(:,:,:,:)
        
        integer :: t, tt
        integer :: n_weight = 3
        integer :: i, j, ii, jj
        integer :: i_nc, j_nc, ii_nc, jj_nc, ii_w, jj_w
        integer :: i_start, j_start, i_end, j_end
        integer :: i_source
        real :: xpos_subgrid, ypos_subgrid
        real :: xpos_subgrid2, ypos_subgrid2
        real :: lon_min, lon_max, lat_min, lat_max

        tt = 1

        allocate (weighting_nc(n_weight,n_weight))
        allocate (subgrid_count_subgrid(emission_max_subgrid_dim(x_dim_index), &
            emission_max_subgrid_dim(y_dim_index),n_pollutant_loop))
        total_proxy_emission_subgrid = 0.0

        do i_source = 1, n_source_index
            if (calculate_source(i_source)) then
                
                emission_subgrid(:,:,:,i_source,:) = 0.0
                subgrid_count_subgrid = 0

                do j = 1, emission_subgrid_dim(y_dim_index,i_source)
                    do i = 1, emission_subgrid_dim(x_dim_index,i_source)

                        ! Only calculate for valid emission subgrids when using the proxy emissions for distribution
                        if (.not. subgrid_emission_distribution_flag .or. sum(proxy_emission_subgrid(i,j,i_source,:)) /= 0) then

                            ! Assumes it is never on the edge of the EMEP grid
                            i_nc = crossreference_emission_to_emep_subgrid(i,j,x_dim_index,i_source)
                            j_nc = crossreference_emission_to_emep_subgrid(i,j,y_dim_index,i_source)

                            weighting_nc = 0.0

                            ! Get subgrid x/y
                            call get_xy(i, j, i_source, xpos_subgrid, ypos_subgrid)

                            ! Calculate the area weighted EMEP grid emissions at each subgrid
                            do jj = -1, 1
                                do ii = -1, 1

                                    ii_nc = ii + i_nc
                                    jj_nc = jj + j_nc
                                    ii_w = ii + 2
                                    jj_w = jj + 2

                                    lon_min = max(xpos_subgrid - dgrid_nc(lon_nc_index)/2.0, &
                                        var1d_nc(ii_nc,lon_nc_index) - dgrid_nc(lon_nc_index)/2.0)
                                    lon_max = min(xpos_subgrid + dgrid_nc(lon_nc_index)/2.0, &
                                        var1d_nc(ii_nc,lon_nc_index) + dgrid_nc(lon_nc_index)/2.0)
                                    lat_min = max(ypos_subgrid - dgrid_nc(lat_nc_index)/2.0, &
                                        var1d_nc(jj_nc,lat_nc_index) - dgrid_nc(lat_nc_index)/2.0)
                                    lat_max = min(ypos_subgrid + dgrid_nc(lat_nc_index)/2.0, &
                                        var1d_nc(jj_nc,lat_nc_index) + dgrid_nc(lat_nc_index)/2.0)

                                    if (lon_max > lon_min .and. lat_max > lat_min) then
                                        weighting_nc(ii_w,jj_w) = &
                                            (lat_max-lat_min)*(lon_max-lon_min)/dgrid_nc(lon_nc_index)/dgrid_nc(lat_nc_index)
                                    else
                                        weighting_nc(ii_w,jj_w) = 0.0
                                    end if

                                    emission_subgrid(i,j,:,i_source,:) = emission_subgrid(i,j,:,i_source,:) &
                                        + var3d_nc(ii_nc,jj_nc,:,emis_nc_index,i_source,:)*weighting_nc(ii_w,jj_w)
                                end do
                            end do

                            ! Calculate the total subgrid emissions in the EMEP grid region surrounding the subgrid
                            i_start = max(1, i - emission_subgrid_loop_index(x_dim_index,i_source))
                            i_end = min(emission_subgrid_dim(x_dim_index,i_source), &
                                i + emission_subgrid_loop_index(x_dim_index,i_source))
                            j_start = max(1, j - emission_subgrid_loop_index(y_dim_index,i_source))
                            j_end = min(emission_subgrid_dim(y_dim_index,i_source), &
                                j + emission_subgrid_loop_index(y_dim_index,i_source))
                            
                            do jj = j_start, j_end
                                do ii = i_start, i_end

                                    call get_xy(ii, jj, i_source, xpos_subgrid2, ypos_subgrid2)

                                    if (abs(xpos_subgrid - xpos_subgrid2) <= dgrid_nc(lon_nc_index)/2.0 &
                                            .and. abs(ypos_subgrid - ypos_subgrid2) <= dgrid_nc(lat_nc_index)/2.0) then

                                        total_proxy_emission_subgrid(i,j,i_source,:) = &
                                            total_proxy_emission_subgrid(i,j,i_source,:) + proxy_emission_subgrid(ii,jj,i_source,:)
                                        subgrid_count_subgrid(i,j,:) = subgrid_count_subgrid(i,j,:) + 1
                                    end if
                                end do
                            end do

                            ! Converts from mg/subgrid to ug/s/subgrid assuming the original EMEP emissions are in mg/m2/hour
                            if (hourly_calculations) then
                                emission_subgrid(i,j,:,i_source,:) = &
                                    emission_subgrid(i,j,:,i_source,:)*emission_subgrid_delta(x_dim_index,i_source) &
                                    *emission_subgrid_delta(y_dim_index,i_source)*1000.0/3600.0
                            end if
                            if (annual_calculations) then
                                emission_subgrid(i,j,:,i_source,:) = &
                                    emission_subgrid(i,j,:,i_source,:)*emission_subgrid_delta(x_dim_index,i_source) &
                                    *emission_subgrid_delta(y_dim_index,i_source)*1000.0/3600.0/EMEP_emission_aggregation_period
                            end if
                        end if
                    end do

                    if (mod(j,100) == 0) then
                        write(*,'(3A,i4,A,i4)') 'Subgrid emission EMEP interpolation for ', &
                            trim(source_file_str(i_source)), ': ', j, ' of ', emission_subgrid_dim(2,i_source)
                    end if
                end do

                ! Distribute EMEP emissions to subgrid emissions
                if (subgrid_emission_distribution_flag) then
                    write(unit_logfile,'(A)') 'Distributing EMEP emissions to subgrid emissions within an area weighted EMEP grid'
                    do t = t_start, t_end
                        emission_subgrid(:,:,t,i_source,:) = emission_subgrid(:,:,t,i_source,:)*subgrid_count_subgrid(:,:,:) &
                            *proxy_emission_subgrid(:,:,i_source,:)/total_proxy_emission_subgrid(:,:,i_source,:)
                        where (total_proxy_emission_subgrid(:,:,i_source,:) == 0.0) emission_subgrid(:,:,t,i_source,:) = 0.0
                    end do
                end if

            end if
        end do
    end subroutine distribute_emep_to_proxy_by_area_weighting

    subroutine distribute_emep_to_proxy_by_emission_weighting(t_start, t_end, total_weighting_nc, proxy_weighting_nc, &
        subgrid_count_subgrid, total_proxy_emission_subgrid, weighting_subgrid)
        !! Emission-weighted moving-window interpolation from EMEP grid to proxy subgrids
        !!
        !! Notes
        !!   - Loop through subgrid and carry out a subgrid weighted moving window interpolation
        !!   - currently does the inteprolation for all time steps which is not ncessary. 
        !!     Only needs to do it for 1. Same with the area interpolation.
        integer, intent(in) :: t_start
        integer, intent(in) :: t_end
        real, allocatable, intent(inout) :: total_weighting_nc(:,:,:) !! EMEP grid weighting for interpolation
        real, allocatable, intent(inout) :: proxy_weighting_nc(:,:,:) !! EMEP grid weighting for interpolation
        integer, allocatable, intent(inout) :: subgrid_count_subgrid(:,:,:)
        real, allocatable, intent(inout) :: total_proxy_emission_subgrid(:,:,:,:)
        real, allocatable, intent(inout) :: weighting_subgrid(:,:,:)

        integer :: i, j, ii, jj
        integer :: i_nc, j_nc
        integer :: i_cross,j_cross
        integer :: t, tt
        integer :: i_w_c, j_w_c
        integer :: i_nc_start, i_nc_end, j_nc_start, j_nc_end
        integer :: i_nc_c, j_nc_c
        integer :: i_start, i_end, j_start, j_end
        integer :: i_source, i_pollutant
        integer :: weighting_subgrid_dim(2,n_source_index)
        real :: xpos_subgrid, ypos_subgrid
        real :: xpos_subgrid2, ypos_subgrid2

        ! Only use time tt for weighting distribution to increase speed. It is possible this can change with time. ??
        ! Replace 'tt' with ':' ??
        tt = 1

        allocate (total_weighting_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),n_pollutant_loop)) 
        allocate (proxy_weighting_nc(5,5,n_pollutant_loop)) 
        allocate (subgrid_count_subgrid(emission_max_subgrid_dim(x_dim_index), &
            emission_max_subgrid_dim(y_dim_index),n_pollutant_loop))

        i_w_c = 3
        j_w_c = 3

        total_proxy_emission_subgrid = 0.0

        allocate (weighting_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),n_pollutant_loop))
        weighting_subgrid_dim(:,:) = emission_subgrid_dim(1:2,:)

        ! Calculate weighting sum for each EMEP grid
        total_weighting_nc = 0.0
        do i_source = 1, n_source_index
            if (calculate_source(i_source)) then

                emission_subgrid(:,:,:,i_source,:) = 0.0
                subgrid_count_subgrid = 0
                weighting_subgrid(:,:,:) = proxy_emission_subgrid(:,:,i_source,:)

                ! Calculate the total weighting (emission) in each emep grid
                do j = 1, weighting_subgrid_dim(y_dim_index,i_source)
                    do i = 1, weighting_subgrid_dim(x_dim_index,i_source)
                        i_nc = crossreference_emission_to_emep_subgrid(i,j,x_dim_index,i_source)
                        j_nc = crossreference_emission_to_emep_subgrid(i,j,y_dim_index,i_source)
                        total_weighting_nc(i_nc,j_nc,:) = total_weighting_nc(i_nc,j_nc,:) + weighting_subgrid(i,j,:)
                    end do
                end do

                ! Calculate the proxy weighting in the nearest emep grids for each subgrid
                do j = 1, emission_subgrid_dim(y_dim_index,i_source)
                    do i = 1, emission_subgrid_dim(x_dim_index,i_source)

                        ! Only calculate for valid emission subgrids when using the proxy for distribution
                        if (.not. subgrid_emission_distribution_flag .or. sum(proxy_emission_subgrid(i,j,i_source,:)) /= 0) then

                            proxy_weighting_nc = 0.0

                            i_cross = i
                            j_cross = j
                            i_nc_c = crossreference_emission_to_emep_subgrid(i_cross,j_cross,x_dim_index,i_source)
                            j_nc_c = crossreference_emission_to_emep_subgrid(i_cross,j_cross,y_dim_index,i_source)
                            
                            ! Limit the loop so that it doesn't go over more than necessary subgrids and does not go outside the domain
                            i_start = max(1, i - emission_subgrid_loop_index(x_dim_index,i_source))
                            i_end = min(emission_subgrid_dim(x_dim_index,i_source), &
                                i + emission_subgrid_loop_index(x_dim_index,i_source))
                            j_start = max(1, j - emission_subgrid_loop_index(y_dim_index,i_source))
                            j_end = min(emission_subgrid_dim(y_dim_index,i_source), &
                                j + emission_subgrid_loop_index(y_dim_index,i_source))

                            call get_xy(i, j, i_source, xpos_subgrid, ypos_subgrid)

                            do jj = j_start, j_end
                                do ii = i_start, i_end

                                    call get_xy(ii, jj, i_source, xpos_subgrid2, ypos_subgrid2)

                                    if (abs(xpos_subgrid - xpos_subgrid2) <= dgrid_nc(lon_nc_index)/2.0 &
                                        .and. abs(ypos_subgrid - ypos_subgrid2) <= dgrid_nc(lat_nc_index)/2.0) then

                                        i_nc = crossreference_emission_to_emep_subgrid(ii,jj,x_dim_index,i_source)
                                        j_nc = crossreference_emission_to_emep_subgrid(ii,jj,y_dim_index,i_source)
                                        proxy_weighting_nc(i_nc-i_nc_c+i_w_c,j_nc-j_nc_c+j_w_c,:) = &
                                            proxy_weighting_nc(i_nc-i_nc_c+i_w_c,j_nc-j_nc_c+j_w_c,:) + weighting_subgrid(ii,jj,:)
                                        
                                        ! Weighting subgrid is the same as the existing proxy subgrid but without the subsource
                                        total_proxy_emission_subgrid(i,j,i_source,:) = &
                                            total_proxy_emission_subgrid(i,j,i_source,:) + weighting_subgrid(ii,jj,:) 
                                        subgrid_count_subgrid(i,j,:) = subgrid_count_subgrid(i,j,:) + 1
                                    end if
                                end do
                            end do

                            i_cross = i
                            j_cross = j
                            i_nc = crossreference_emission_to_emep_subgrid(i_cross,j_cross,x_dim_index,i_source)
                            j_nc = crossreference_emission_to_emep_subgrid(i_cross,j_cross,y_dim_index,i_source)
                            i_nc_start = max(1, i_nc-1)
                            i_nc_end = min(dim_length_nc(x_dim_nc_index), i_nc+1)
                            j_nc_start = max(1, j_nc-1)
                            j_nc_end = min(dim_length_nc(y_dim_nc_index), j_nc+1)

                            do jj = j_nc_start, j_nc_end
                                do ii = i_nc_start, i_nc_end
                                    proxy_weighting_nc(ii-i_nc+i_w_c,jj-j_nc+j_w_c,:) = &
                                        proxy_weighting_nc(ii-i_nc+i_w_c,jj-j_nc+j_w_c,:)/total_weighting_nc(ii,jj,:)
                                    where (total_weighting_nc(ii,jj,:) == 0.0) &
                                        proxy_weighting_nc(ii-i_nc+i_w_c,jj-j_nc+j_w_c,:) = 0.0
                                end do
                            end do

                            ! Add up the contributing weights
                            ! Note that only the subsource_index=1 can be determined since EMEP has no subsources
                            do jj = j_nc_start, j_nc_end
                                do ii = i_nc_start, i_nc_end
                                    do i_pollutant = 1, n_pollutant_loop
                                        emission_subgrid(i,j,:,i_source,i_pollutant) = &
                                            emission_subgrid(i,j,:,i_source,i_pollutant) &
                                            + var3d_nc(ii,jj,:,emis_nc_index,i_source,i_pollutant) &
                                            *proxy_weighting_nc(ii-i_nc+i_w_c,jj-j_nc+j_w_c,i_pollutant)
                                    end do
                                end do
                            end do

                            ! Converts from mg/subgrid to ug/s/subgrid assuming the original EMEP emissions are in mg/m2/hour
                            if (hourly_calculations) then
                                emission_subgrid(i,j,:,i_source,:) = emission_subgrid(i,j,:,i_source,:) &
                                    *emission_subgrid_delta(x_dim_index,i_source)*emission_subgrid_delta(y_dim_index,i_source) &
                                    *1000.0/3600.0
                            end if
                            if (annual_calculations) then
                                emission_subgrid(i,j,:,i_source,:) = emission_subgrid(i,j,:,i_source,:) &
                                    *emission_subgrid_delta(x_dim_index,i_source)*emission_subgrid_delta(y_dim_index,i_source) &
                                    *1000.0/3600.0/EMEP_emission_aggregation_period
                            end if
                        end if
                    end do
                    if (mod(j,100) == 0) then
                        write(*,'(3A,i4,A,i4)') 'Subgrid emission EMEP interpolation for ', &
                            trim(source_file_str(i_source)), ': ', j, ' of ', emission_subgrid_dim(2,i_source)
                    end if
                end do

                ! Distribute EMEP emissions to existing subgrid emissions
                if (subgrid_emission_distribution_flag) then
                    write(unit_logfile,'(A)') &
                        'Distributing EMEP emissions to subgrid emissions within an emission weighted EMEP grid'
                    do t = t_start, t_end
                        emission_subgrid(:,:,t,i_source,:) = emission_subgrid(:,:,t,i_source,:)*subgrid_count_subgrid(:,:,:) &
                            *proxy_emission_subgrid(:,:,i_source,:)/total_proxy_emission_subgrid(:,:,i_source,:)
                        where (total_proxy_emission_subgrid(:,:,i_source,:) == 0.0) emission_subgrid(:,:,t,i_source,:) = 0.0
                    end do
                end if

            end if
        end do
    end subroutine distribute_emep_to_proxy_by_emission_weighting

    subroutine get_xy(i, j, i_source, xpos, ypos)

        !! Convert emission subgrid lon/lat to projected coordinates used for overlap and distance calculations
        !!
        !! Notes
        !!   - Selects conversion based on EMEP_projection_type:
        !!     - LL : returns lon/lat directly
        !!     - LCC : calls lb2lambert_uEMEP / lb2lambert2_uEMEP
        !!     - PS : calls LL2PS_spherical
        !!   - Stops with an error if EMEP_projection_type is unrecognised

        integer, intent(in) :: i, j
        integer, intent(in) :: i_source
        real, intent(out) :: xpos, ypos

        select case(EMEP_projection_type)
        case(LL_projection_index)
            xpos = lon_emission_subgrid(i,j,i_source)
            ypos = lat_emission_subgrid(i,j,i_source)
        case(LCC_projection_index)
            if (use_alternative_LCC_projection_flag) then
                call lb2lambert2_uEMEP(xpos, ypos, lon_emission_subgrid(i,j,i_source), &
                    lat_emission_subgrid(i,j,i_source),EMEP_projection_attributes)
            else
                call lb2lambert_uEMEP(xpos, ypos, lon_emission_subgrid(i,j,i_source), &
                    lat_emission_subgrid(i,j,i_source), real(EMEP_projection_attributes(3)), &
                    real(EMEP_projection_attributes(4)))
            end if
        case(PS_projection_index)
            call LL2PS_spherical(xpos, ypos, lon_emission_subgrid(i,j,i_source), &
                lat_emission_subgrid(i,j,i_source), EMEP_projection_attributes)
        case default
            write(unit_logfile, "(a,i0)") "ERROR: Unknown projection index: ", EMEP_projection_type
            stop 1
        end select
    end subroutine get_xy

end module subgrid_emission_emep
