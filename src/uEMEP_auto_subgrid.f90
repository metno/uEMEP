module auto_subgrid

    use uemep_configuration
    use uEMEP_definitions
    use mod_area_interpolation, only: area_weighted_interpolation_function
    use mod_lambert_projection, only: LL2PROJ, PROJ2LL
    use netcdf

    implicit none
    private

    public :: uEMEP_auto_subgrid, uEMEP_region_mask, uEMEP_interpolate_auto_subgrid, nlreg_uEMEP_region_mask_new

contains

    subroutine uEMEP_auto_subgrid()
        !! uEMEP model uEMEP_auto_subgrid
        !! Automatically creates a grid dependent on the distance to source
        integer :: i,j,k
        integer :: i_source
        integer :: t
        real :: max_use_subgrid_size(n_source_index)
        real :: use_subgrid_range(n_source_index)
        integer :: use_subgrid_step
        integer :: use_emission_subgrid_space
        integer :: i_cross, j_cross
        integer :: i_start, i_end, j_start, j_end
        integer :: i_start_k, i_end_k, j_start_k, j_end_k
        real :: sum_emission

        ! Exit if emission positions are not to be used and if one of the other auto positions is to be used
        if ( .not. use_emission_positions_for_auto_subgrid_flag(allsource_index) .and. (use_population_positions_for_auto_subgrid_flag .or. use_receptor_positions_for_auto_subgrid_flag)) then
            return
        end if

        ! Exit and fill subgrid if none of the auto subgrids are to be used
        if ( .not. use_emission_positions_for_auto_subgrid_flag(allsource_index) .and. .not. use_population_positions_for_auto_subgrid_flag .and. .not. use_receptor_positions_for_auto_subgrid_flag) then
            use_subgrid = .true.
            return
        end if

        write(unit_logfile,'(A)') ''
        write(unit_logfile,'(A)') '================================================================'
        write(unit_logfile,'(A)') 'Setting subgrids automatically (uEMEP_auto_subgrid)'
        write(unit_logfile,'(A)') '================================================================'

        ! Set all subgrids to do not use
        use_subgrid_val = 0
        use_subgrid_interpolation_index = -1

        ! Set time index used for emissions to 1, so only tests for the first hour if there are emissions
        t = 1
        
        ! Set the maximum grid size
        max_use_subgrid_size = max_interpolation_subgrid_size
        use_subgrid_range = 8.0
        use_subgrid_range(traffic_index) = 4.0
        use_subgrid_range(shipping_index) = 12.0

        n_use_subgrid_levels = 0
        
        ! Sets the auto gridding
        if (subgrid_delta(x_dim_index) .eq. 100.0) then
            use_subgrid_step_delta(0) = 1
            use_subgrid_step_delta(1) = 2
            use_subgrid_step_delta(2) = 5
            use_subgrid_step_delta(3) = 10
            use_subgrid_step_delta(4) = 20
            if (max_interpolation_subgrid_size .eq. 250.0) n_use_subgrid_levels = 1
            if (max_interpolation_subgrid_size .eq. 500.0) n_use_subgrid_levels = 2
            if (max_interpolation_subgrid_size .eq. 1000.0) n_use_subgrid_levels = 3
            if (max_interpolation_subgrid_size .eq. 2000.0) n_use_subgrid_levels = 4 
            use_subgrid_range(traffic_index) = 6.0
        else if (subgrid_delta(x_dim_index) .eq. 25.0) then
            use_subgrid_step_delta(0) = 1
            use_subgrid_step_delta(1) = 2
            use_subgrid_step_delta(2) = 4
            use_subgrid_step_delta(3) = 10
            use_subgrid_step_delta(4) = 20
            use_subgrid_step_delta(5) = 40
            use_subgrid_step_delta(5) = 80
            if (max_interpolation_subgrid_size .eq. 250.0) n_use_subgrid_levels = 2
            if (max_interpolation_subgrid_size .eq. 500.0) n_use_subgrid_levels = 3
            if (max_interpolation_subgrid_size .eq. 1000.0) n_use_subgrid_levels = 4
            if (max_interpolation_subgrid_size .eq. 2000.0) n_use_subgrid_levels = 5
            use_subgrid_range(traffic_index) = 4.0
        else if (subgrid_delta(x_dim_index) .eq. 50.0) then
            use_subgrid_step_delta(0) = 1
            use_subgrid_step_delta(1) = 2
            use_subgrid_step_delta(2) = 4
            use_subgrid_step_delta(3) = 10
            use_subgrid_step_delta(4) = 20
            use_subgrid_step_delta(5) = 40
            if (max_interpolation_subgrid_size .eq. 250.0) n_use_subgrid_levels = 2
            if (max_interpolation_subgrid_size .eq. 500.0) n_use_subgrid_levels = 3
            if (max_interpolation_subgrid_size .eq. 1000.0) n_use_subgrid_levels = 4
            if (max_interpolation_subgrid_size .eq. 2000.0) n_use_subgrid_levels = 5
            use_subgrid_range(traffic_index) = 6.0
        else if (subgrid_delta(x_dim_index) .eq. 250.0) then
            use_subgrid_step_delta(0) = 1
            use_subgrid_step_delta(1) = 2
            use_subgrid_step_delta(2) = 4
            use_subgrid_step_delta(3) = 8
            if (max_interpolation_subgrid_size .eq. 250.) n_use_subgrid_levels = 0
            if (max_interpolation_subgrid_size .eq. 500.) n_use_subgrid_levels = 1
            if (max_interpolation_subgrid_size .eq. 1000.) n_use_subgrid_levels = 2
            if (max_interpolation_subgrid_size .eq. 2000.) n_use_subgrid_levels = 3
            use_subgrid_range(traffic_index) = 8.0
        else
            write(unit_logfile,'(a,f12.2)') 'When auto gridding target grid sizes must be either 25, 50, 100 or 250 metres. Stopping because target grid size is:', subgrid_delta(x_dim_index)
            stop
        end if

        if (n_use_subgrid_levels(allsource_index).eq.0) then
            write(unit_logfile,'(a,f12.2)') 'When auto gridding maximum grid sizes must be either 250, 500, 1000 or 2000 metres. Stopping because max grid size is:',max_interpolation_subgrid_size
            stop 1
        end if

        ! Set the number of levels to match this
        do i_source = 1, n_source_index
            if (calculate_source(i_source)) then
                write(*,*) 'Using auto subgrid for source ', trim(source_file_str(i_source)), use_emission_positions_for_auto_subgrid_flag(i_source)
            end if
        end do

        do i_source = 1, n_source_index
            if (calculate_source(i_source) .and. use_emission_positions_for_auto_subgrid_flag(i_source)) then
                write(unit_logfile,'(a,2f10.1,i10)') trim(source_file_str(i_source))//': maximum use subgrid size (m), grid range and number of levels: ', &
                    max_use_subgrid_size(i_source), use_subgrid_range(i_source), n_use_subgrid_levels(i_source)
                ! Fill the interpolation index with the highest level
                do k = n_use_subgrid_levels(i_source), 0, -1
                    use_subgrid_step = 2**k
                    use_subgrid_step = use_subgrid_step_delta(k)
                    use_emission_subgrid_space = floor(use_subgrid_step/sqrt(emission_subgrid_delta(x_dim_index,i_source)*emission_subgrid_delta(y_dim_index,i_source) &
                        /subgrid_delta(x_dim_index)/subgrid_delta(y_dim_index))/2.0*use_subgrid_range(i_source))
                    use_emission_subgrid_space = max(1, use_emission_subgrid_space)

                    do j = 1, subgrid_dim(y_dim_index), use_subgrid_step
                        do i = 1, subgrid_dim(x_dim_index), use_subgrid_step

                            i_cross=crossreference_target_to_emission_subgrid(i,j,x_dim_index,i_source)
                            j_cross=crossreference_target_to_emission_subgrid(i,j,y_dim_index,i_source)

                            ! Search in the neighbourhood and find the sum of the emissions
                            i_start = max(1, i_cross - use_emission_subgrid_space)
                            i_end = min(emission_subgrid_dim(x_dim_index,i_source), i_cross + use_emission_subgrid_space)
                            j_start = max(1, j_cross - use_emission_subgrid_space)
                            j_end = min(emission_subgrid_dim(y_dim_index,i_source), j_cross + use_emission_subgrid_space)
                            sum_emission = sum(proxy_emission_subgrid(i_start:i_end,j_start:j_end,i_source,:))

                            ! Select those with emission sums > 0 but always include the coarsest level everywhere
                            if (sum_emission .gt. 0. .or. k .eq. n_use_subgrid_levels(i_source)) then
                                use_subgrid_val(i,j,i_source) = 1
                                
                                ! Label the grids for interpolation later
                                i_start_k = max(1, i - use_subgrid_step/2)
                                i_end_k = min(subgrid_dim(x_dim_index), i + use_subgrid_step/2)
                                j_start_k = max(1, j - use_subgrid_step/2)
                                j_end_k = min(subgrid_dim(y_dim_index), j + use_subgrid_step/2)

                                use_subgrid_interpolation_index(i_start_k:i_end_k,j_start_k:j_end_k,i_source) = k
                            end if
                        end do
                    end do
                end do
                use_subgrid_val(:,:,allsource_index)=use_subgrid_val(:,:,allsource_index)+use_subgrid_val(:,:,i_source)

                ! Check
                do j = 1, subgrid_dim(y_dim_index)
                    do i = 1, subgrid_dim(x_dim_index)
                        if (use_subgrid_interpolation_index(i,j,i_source) .lt. 0) write(*,*) i, j, use_subgrid_interpolation_index(i,j,i_source)
                    end do
                end do
            end if
        end do

        use_subgrid_val = min(1, use_subgrid_val)

        ! Convert values to logical
        do i_source = 1, n_source_index
            if (use_emission_positions_for_auto_subgrid_flag(i_source)) then
                if (calculate_source(i_source) .or. i_source .eq. allsource_index) then
                    do j = 1, subgrid_dim(y_dim_index)
                        do i = 1, subgrid_dim(x_dim_index)
                            if (use_subgrid_val(i,j,i_source) .gt. 0) then
                                use_subgrid(i,j,i_source) = .true.
                            else
                                use_subgrid(i,j,i_source) = .false.
                            end if
                        end do
                    end do
                end if
            else
                ! If not to be auto gridded then set use to true
                use_subgrid(:,:,i_source) = .true.
            end if
        end do

        do i_source = 1, n_source_index
            if (calculate_source(i_source) .and. use_emission_positions_for_auto_subgrid_flag(i_source)) then
                write(unit_logfile,'(a,2i10,f6.1)') 'Number of calculation subgrids for '//trim(source_file_str(i_source))//' (number, total, percent):', &
                    sum(use_subgrid_val(:,:,i_source)), subgrid_dim(1)*subgrid_dim(2), sum(use_subgrid_val(:,:,i_source))*100.0/(subgrid_dim(1)*subgrid_dim(2))
            end if
        end do

    end subroutine uEMEP_auto_subgrid

    subroutine uEMEP_interpolate_auto_subgrid()
        ! This is the corresponding routine for interpolating the auto selected data
        integer :: xdim, ydim
        real :: delta(2)
        real :: xval, yval
        real :: xgrid(3,3), ygrid(3,3), zgrid(3,3)

        integer :: i, j, k, ii, jj
        integer :: i_source, i_pollutant, t
        integer :: i_in(3), j_in(3)
        integer :: use_subgrid_step

        write(unit_logfile,'(A)') ''
        write(unit_logfile,'(A)') '================================================================'
        write(unit_logfile,'(A)') 'Interpolating auto grid (uEMEP_interpolate_auto_subgrid)'
        write(unit_logfile,'(A)') '================================================================'

        xdim = 3
        ydim = 3

        do i_source = 1, n_source_index
            if (calculate_source(i_source) .and. i_source .ne. allsource_index .and. use_emission_positions_for_auto_subgrid_flag(i_source)) then
                ! Only interpolate for chosen sources that have been auto gridded based on emissions
                do k = n_use_subgrid_levels(i_source), 1, -1
                    use_subgrid_step = use_subgrid_step_delta(k)
                    delta = use_subgrid_step*subgrid_delta
                    do j = 1, subgrid_dim(y_dim_index)
                        do i = 1, subgrid_dim(x_dim_index)
                            xval = x_subgrid(i,j)
                            yval = y_subgrid(i,j)
                            ! Only do the interpolation if it is the right interpolation_index, it is inside the region and it is not a valid subgrid
                            ! TODO: use a logical flag here?
                            ! if (use_subgrid_interpolation_index(i,j,i_source).eq.k.and.use_subgrid_val(i,j,i_source).ne.outside_interpolation_region_index.and..not.use_subgrid_val(i,j,i_source)) then
                            ! Do it everywhere in the grid
                            if (use_subgrid_interpolation_index(i,j,i_source) .eq. k .and. use_subgrid_val(i,j,i_source) .ne. &
                                outside_interpolation_region_index .and. .not. use_subgrid(i,j,i_source)) then

                                i_in(2) = floor(real(i - 1)/use_subgrid_step + 0.5)*use_subgrid_step + 1
                                i_in(2) = min(subgrid_dim(x_dim_index), max(1, i_in(2)))
                                i_in(1) = i_in(2) - use_subgrid_step
                                i_in(3) = i_in(2) + use_subgrid_step

                                j_in(2) = floor(real(j - 1)/use_subgrid_step + 0.5)*use_subgrid_step + 1
                                j_in(2) = min(subgrid_dim(y_dim_index), max(1, j_in(2)))
                                j_in(1) = j_in(2) - use_subgrid_step
                                j_in(3) = j_in(2) + use_subgrid_step

                                if (i_in(1) .lt. 1) i_in(1) = i_in(2)
                                if (j_in(1) .lt. 1) j_in(1) = j_in(2)
                                if (i_in(3) .gt. subgrid_dim(x_dim_index)) i_in(3) = i_in(2)
                                if (j_in(3) .gt. subgrid_dim(y_dim_index)) j_in(3) = j_in(2)

                                do t = 1, subgrid_dim(t_dim_index)
                                    do i_pollutant = 1, n_pollutant_loop
                                        do jj = 1, 3
                                            do ii = 1, 3
                                                xgrid(ii,jj) = x_subgrid(i_in(2),j_in(2)) + (ii - 2)*delta(1)
                                                ygrid(ii,jj) = y_subgrid(i_in(2),j_in(2)) + (jj - 2)*delta(2)
                                                zgrid(ii,jj) = subgrid(i_in(ii),j_in(jj),t,proxy_subgrid_index,i_source,i_pollutant)
                                                if (i_in(ii) .lt. 1 .or. i_in(ii) .gt. subgrid_dim(x_dim_index) .or. j_in(jj) .lt. 1 .or. j_in(jj) .gt. subgrid_dim(y_dim_index)) then
                                                    zgrid(ii,jj) = subgrid(i_in(2),j_in(2),t,proxy_subgrid_index,i_source,i_pollutant)
                                                else if (use_subgrid_val(i_in(ii),j_in(jj),i_source).eq.outside_interpolation_region_index) then
                                                    zgrid(ii,jj) = subgrid(i_in(2),j_in(2),t,proxy_subgrid_index,i_source,i_pollutant)
                                                end if
                                            end do
                                        end do
                                        subgrid(i,j,t,proxy_subgrid_index,i_source,i_pollutant)=area_weighted_interpolation_function(xgrid,ygrid,zgrid,xdim,ydim,delta,xval,yval)

                                        ! Travel time interpolation as well
                                        do jj = 1, 3
                                            do ii = 1, 3
                                                zgrid(ii,jj) = traveltime_subgrid(i_in(ii),j_in(jj),t,1,i_pollutant)
                                                if (i_in(ii) .lt. 1 .or. i_in(ii) .gt. subgrid_dim(x_dim_index) .or. j_in(jj) .lt. 1 .or. j_in(jj) .gt. subgrid_dim(y_dim_index)) then
                                                    zgrid(ii,jj) = traveltime_subgrid(i_in(2),j_in(2),t,1,i_pollutant)
                                                else if (use_subgrid_val(i_in(ii),j_in(jj),i_source) .eq. outside_interpolation_region_index) then
                                                    zgrid(ii,jj) = traveltime_subgrid(i_in(2),j_in(2),t,1,i_pollutant)
                                                end if
                                            end do
                                        end do
                                        traveltime_subgrid(i,j,t,1,i_pollutant) = area_weighted_interpolation_function(xgrid, ygrid, zgrid, xdim, ydim, delta, xval, yval)

                                        do jj = 1, 3
                                            do ii = 1, 3
                                                zgrid(ii,jj) = traveltime_subgrid(i_in(ii),j_in(jj),t,2,i_pollutant)
                                                if (i_in(ii) .lt. 1 .or. i_in(ii) .gt. subgrid_dim(x_dim_index) .or. j_in(jj) .lt. 1 .or. j_in(jj) .gt. subgrid_dim(y_dim_index)) then
                                                    zgrid(ii,jj) = traveltime_subgrid(i_in(2),j_in(2),t,2,i_pollutant)
                                                else if (use_subgrid_val(i_in(ii),j_in(jj),i_source) .eq. outside_interpolation_region_index) then
                                                    zgrid(ii,jj) = traveltime_subgrid(i_in(2),j_in(2),t,2,i_pollutant)
                                                end if
                                            end do
                                        end do
                                        traveltime_subgrid(i,j,t,2,i_pollutant) = area_weighted_interpolation_function(xgrid, ygrid, zgrid, xdim, ydim, delta, xval, yval)
                                    end do
                                end do
                            end if
                        end do
                    end do
                end do
            end if
        end do

        ! Reset the use_subgrid values so chemistry and exposure happens everywhere but not outside the region
        use_subgrid(:,:,allsource_index) = .true.
        where (use_subgrid_val(:,:,allsource_index) .eq. outside_region_index) use_subgrid(:,:,allsource_index) = .false.

    end subroutine uEMEP_interpolate_auto_subgrid

    subroutine uEMEP_region_mask()
        ! This routine defines use_subgrid_val=2 for regions outside the selected region. This can be used to control the interpolation routine
        ! It also sets use_subgrid=.false. outside the region so that no calculations are made there either
        character(256) :: temp_name, temp_str, temp_str1
        integer, allocatable :: tile_municipality_subgrid(:,:,:)
        integer :: municipality_id
        real :: x_ssb, f_easting, ssb_dx, y_ssb, ssb_dy
        integer :: count
        integer(kind=8) :: ssb_id
        integer :: i, j, k, i_range, j_range, i_range_interp, j_range_interp, i_tile, j_tile
        logical :: exists
        integer :: unit_in
        integer :: index_val
        integer :: i_source
        character(256) :: region_number_str
        integer, parameter :: n_search = 5
        character(16) :: search_str(n_search)
        real :: search_delta(n_search)
        integer :: temp_search
        integer :: i_tile_region, j_tile_region
        integer :: i_range_region, j_range_region
        integer :: total_grids = 0
        integer :: source_count = 0
        integer :: io

        search_str = ['_1000m', '_500m ', '_250m ', '_100m ', '_50m  ']
        search_delta = [1000.0, 500.0, 250.0, 100.0, 50.0]

        write(unit_logfile,'(A)') ''
        write(unit_logfile,'(A)') '================================================================'
        write(unit_logfile,'(A)') 'Masking region (uEMEP_region_mask)'
        write(unit_logfile,'(A)') '================================================================'

        if ( .not. allocated(tile_municipality_subgrid)) allocate (tile_municipality_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),2))
        tile_municipality_subgrid = 0

        i_source = allsource_index

        ! Search file name to define the grid size
        ssb_dx = 0.0; ssb_dy = 0.0
        do k = 1, n_search
            temp_search = index(filename_population(municipality_index), trim(adjustl(search_str(k))))
            if (temp_search .ne. 0) then
                ssb_dx = search_delta(k)
                ssb_dy = search_delta(k)
                write(unit_logfile,'(i,A)') temp_search, ' Reading municipality data with resolution '//trim(adjustl(search_str(k)))
            end if
        end do

        if (ssb_dx .eq. 0) then
            write(unit_logfile,'(A)') 'Cannot find a valid SSB grid size. Stopping. '//trim(filename_population(municipality_index))
            stop 1
        end if

        region_number_str = ''
        write(region_number_str,*) region_index
        region_number_str = trim(region_number_str)//'_'

        ! Read in SSB file containing gridded municipality ids
        f_easting=2.0e6
        pathfilename_population(municipality_index) = trim(pathname_population(municipality_index))//trim(adjustl(region_number_str))//trim(filename_population(municipality_index))
        
        ! Test existence of the heating filename. If does not exist then use default
        inquire(file=trim(pathfilename_population(municipality_index)), exist=exists)
        if ( .not. exists) then
            write(unit_logfile,'(A,A)') ' ERROR: Masking region SSB file with municipality IDs does not exist: ', trim(pathfilename_population(municipality_index))
            stop 1
        end if
        temp_name = pathfilename_population(municipality_index)
        
        ! Open the file for reading
        unit_in = 20
        open(unit_in, file=temp_name, access='sequential', status='old', readonly)
        write(unit_logfile,'(a)') ' Opening SSB municipality file '//trim(temp_name)
        rewind(unit_in)
        
        ! Read header SSBID0250M;kommunenum
        read(unit_in,'(A)') temp_str
        write(unit_logfile,'(A)') 'Header: '//trim(temp_str)
        count = 0
        i_range = ceiling(ssb_dx/subgrid_delta(x_dim_index)/2.0) + 1
        j_range = ceiling(ssb_dy/subgrid_delta(y_dim_index)/2.0) + 1
        i_range_interp = ceiling(ssb_dx/subgrid_delta(x_dim_index)/2.0 + max_interpolation_subgrid_size/subgrid_delta(x_dim_index)) + 1
        j_range_interp = ceiling(ssb_dy/subgrid_delta(y_dim_index)/2.0 + max_interpolation_subgrid_size/subgrid_delta(y_dim_index)) + 1
        do
            ssb_id = 0; municipality_id = 0
            ! Read in file string
            read(unit_in,'(A)',iostat=io) temp_str
            if (io /= 0) exit
            index_val = index(temp_str,';', back=.false.)
            temp_str1 = temp_str(1:index_val-1)
            temp_str = temp_str(index_val+1:)
            if (index_val .gt. 1) read(temp_str1,*) ssb_id
            read(temp_str,*) municipality_id
            ! If this ssb municipality grid has the correct ID then find the target grid that matches it
            if (municipality_id .eq. region_id) then
                count = count + 1

                !Convert id to grid centre coordinates that are already in UTM33 for SSB data
                x_ssb = floor(ssb_id/10000000.0) - f_easting + ssb_dx/2.0
                y_ssb = mod(ssb_id, 10000000) + ssb_dy/2.0

                ! Find the tile this ssb grid is in
                i_tile = 1 + floor((x_ssb - subgrid_min(x_dim_index))/subgrid_delta(x_dim_index)) ! New definition
                j_tile = 1 + floor((y_ssb - subgrid_min(y_dim_index))/subgrid_delta(y_dim_index)) !New definition

                do j = j_tile-j_range, j_tile+j_range
                    do i = i_tile-i_range, i_tile+i_range
                        ! Make sure i and j are inside a valid range
                        if (i .ge. 1 .and. i .le. subgrid_dim(x_dim_index) .and. j .ge. 1 .and. j .le. subgrid_dim(y_dim_index)) then
                            ! Find the target subgrid within the ssb grid and the region
                            if (x_subgrid(i,j) .ge. x_ssb-ssb_dx/2.0 .and. x_subgrid(i,j) .lt. x_ssb + ssb_dx/2.0 .and. &
                                y_subgrid(i,j) .ge. y_ssb-ssb_dy/2.0 .and. y_subgrid(i,j) .lt. y_ssb+ssb_dy/2.0) then
                                !If there is a target subgrid within this range then allocate the mask value if it is not the correct region_id
                                tile_municipality_subgrid(i,j,1) = 1
                            end if
                        end if
                    end do
                end do

                do j = j_tile-j_range_interp, j_tile+j_range_interp
                    do i = i_tile-i_range_interp, i_tile+i_range_interp
                        ! Make sure i and j are inside a valid range
                        if (i .ge. 1 .and. i .le. subgrid_dim(x_dim_index) .and. j .ge. 1 .and. j .le. subgrid_dim(y_dim_index)) then
                            ! Find the target subgrid within the ssb grid and the region
                            if (x_subgrid(i,j) .ge. x_ssb-ssb_dx/2.0 - max_interpolation_subgrid_size .and. x_subgrid(i,j) .le. x_ssb + ssb_dx/2.0 + max_interpolation_subgrid_size .and. &
                                y_subgrid(i,j) .ge. y_ssb-ssb_dy/2.0 - max_interpolation_subgrid_size .and. y_subgrid(i,j) .le. y_ssb+ssb_dy/2.0 + max_interpolation_subgrid_size) then
                                ! If there is a target subgrid within this range then allocate the mask value if it is not the correct region_id
                                tile_municipality_subgrid(i,j,2) = 1
                            end if
                        end if
                    end do
                end do

                ! Specify the emissions grids that will be used for regional calculations
                if (trace_emissions_from_in_region) then
                    source_count = 0
                    do i_source = 1, n_source_index
                        if (calculate_source(i_source)) then
                            i_tile_region = 1 + floor((x_ssb-emission_subgrid_min(x_dim_index,i_source))/emission_subgrid_delta(x_dim_index,i_source))
                            j_tile_region = 1 + floor((y_ssb-emission_subgrid_min(y_dim_index,i_source))/emission_subgrid_delta(y_dim_index,i_source))
                            i_range_region = ceiling(ssb_dx/emission_subgrid_delta(x_dim_index,i_source)/2.0) + 1
                            j_range_region = ceiling(ssb_dy/emission_subgrid_delta(y_dim_index,i_source)/2.0) + 1

                            do j = j_tile_region-j_range_region, j_tile_region+j_range_region
                                do i = i_tile_region-i_range_region, i_tile_region+i_range_region
                                    ! Make sure i and j are inside a valid range
                                    if (i .ge. 1 .and. i .le. emission_subgrid_dim(x_dim_index,i_source) .and. j .ge. 1 .and. j .le. emission_subgrid_dim(y_dim_index,i_source)) then
                                        ! Find the target subgrid within the ssb grid and the region
                                        if (x_emission_subgrid(i,j,i_source) .ge. x_ssb-ssb_dx/2.0 .and. x_emission_subgrid(i,j,i_source) .lt. x_ssb+ssb_dx/2.0 .and. &
                                            y_emission_subgrid(i,j,i_source) .ge. y_ssb-ssb_dy/2.0 .and. y_emission_subgrid(i,j,i_source) .lt. y_ssb+ssb_dy/2.0) then
                                            ! Default is false
                                            use_subgrid_region(i,j,i_source) = .true.
                                            total_grids = total_grids + 1
                                        end if
                                    end if
                                end do
                            end do
                            source_count = source_count + 1
                        end if
                    end do
                end if
            end if
        end do
        close(unit_in)

        do j = 1, subgrid_dim(y_dim_index)
            do i = 1, subgrid_dim(x_dim_index)
                if (tile_municipality_subgrid(i,j,2) .eq. 0) then
                    use_subgrid_val(i,j,:) = outside_interpolation_region_index
                    use_subgrid(i,j,:) = .false.
                end if
                if (tile_municipality_subgrid(i,j,1) .eq. 0) then
                    use_subgrid_val(i,j,:) = outside_region_index
                end if
            end do
        end do

        if (allocated(tile_municipality_subgrid)) deallocate (tile_municipality_subgrid)

    end subroutine uEMEP_region_mask

    subroutine nlreg_uEMEP_region_mask_new()

        ! Variables used for reading the region mask netcdf file
        ! help-parameters for reading the file
        character(256) pathfilename_region_mask
        logical exists
        integer status_nc
        integer id_nc,var_id_nc,x_dim_id_nc,y_dim_id_nc
        integer temp_num_dims
        ! projection of the region mask
        integer region_mask_projection_type
        double precision region_mask_projection_attributes(10)
        ! the x and y coordinates
        real, allocatable :: x_values_regionmask(:)
        real, allocatable :: y_values_regionmask(:)
        ! the region ID data themselves
        integer, allocatable :: region_mask(:, :)
        ! assumed names for dimensions of the region mask file
        character(256) x_dim_name_regionmask,y_dim_name_regionmask,dimname_temp
        ! length of dimensions of the region mask file
        integer nx_regionmask,ny_regionmask
        ! resolution of the region mask file (to be verified is constant!)
        real dx_regionmask,dy_regionmask
        ! grid dimension looping indices
        integer i,j,ii,jj,i_sub,j_sub

        ! Variables used when applying the region mask to the EMEP and uEMEP grids
        ! x and y value of current location in the region mask projection
        real x_location,y_location
        ! index of closest grid-cell in the region mask to the current location
        integer x_index,y_index
        ! grid resolution in EMEP
        real dx_emep,dy_emep
        ! x and y value of midpoint of an EMEP grid and of subsamples of the EMEP grid
        real x_emepmid,y_emepmid
        real x_emepsub,y_emepsub
        ! corresponding lon and lat
        real lon_emepsub,lat_emepsub

        ! Additional variables used for creating list of regions and calculating fraction of EMEP cells in each region
        integer region_index
        integer, allocatable :: temp_region_ids(:)
        integer, allocatable :: temp_region_ids_dummy(:)
        integer counter
        integer current_region_id,previous_region_id
        integer list_allocated_size
        logical current_region_already_found
        ! number of elements in list of unique regions to make space for at a time
        integer, parameter :: n_regions_allocate = 20


        write(unit_logfile,'(A)') ''
        write(unit_logfile,'(A)') '================================================================'
        write(unit_logfile,'(A)') 'Creating region mask arrays (nlreg_uEMEP_region_mask_new)'
        write(unit_logfile,'(A)') '================================================================'

        ! Ensure arrays are not already allocated
        if (allocated(nlreg_subgrid_region_id)) then
            deallocate(nlreg_subgrid_region_id)
        end if
        if (allocated(nlreg_EMEP_subsample_region_id)) then
            deallocate(nlreg_EMEP_subsample_region_id)
        end if
        if (allocated(nlreg_regionfraction_per_EMEP_grid)) then
            deallocate(nlreg_regionfraction_per_EMEP_grid)
        end if
        if (allocated(nlreg_region_ids)) then
            deallocate(nlreg_region_ids)
        end if

        ! Read the region mask netcdf file (implementation based on uEMEP_read_EMEP)

        ! determine full filename
        pathfilename_region_mask = trim(nlreg_pathname_region_mask)//trim(nlreg_filename_region_mask)

        !Test existence of the region mask file. If does not exist then stop
        inquire(file=trim(pathfilename_region_mask),exist=exists)
        if (.not.exists) then
            write(unit_logfile,'(A,A)') ' ERROR: Netcdf file does not exist: ', trim(pathfilename_region_mask)
            write(unit_logfile,'(A)') '  STOPPING'
            stop
        end if

        !Open the netcdf file for reading
        write(unit_logfile,'(2A)') ' Opening netcdf file: ',trim(pathfilename_region_mask)
        status_nc = NF90_OPEN (pathfilename_region_mask, nf90_nowrite, id_nc)
        if (status_nc /= NF90_NOERR) write(unit_logfile,'(A,I0)') 'ERROR opening netcdf file: ',status_nc

        ! Initialize projection type to longitude-latitude
        region_mask_projection_type = LL_projection_index

        ! Check if it is a UTM projection
        status_nc = NF90_INQ_VARID (id_nc,'projection_utm',var_id_nc)
        if (status_nc == NF90_NOERR) then
            status_nc = nf90_get_att(id_nc, var_id_nc, 'longitude_of_central_meridian', region_mask_projection_attributes(2))
            ! calculate the UTM zone
            region_mask_projection_attributes(1)=(region_mask_projection_attributes(2)+180+3)/6
            ! verify this is an exact UTM zone
            if (abs(region_mask_projection_attributes(1) - floor(region_mask_projection_attributes(1)+0.1)) > 1e-6) then
                write(unit_logfile, '(3A,f,A)') ' ERROR during reading of attributes for "projection_utm" in file ',trim(pathfilename_region_mask),': "longitude_of_central_meridian" = ',region_mask_projection_attributes(2),', which does not correspond to a UTM zone.'
                ! IS THIS STRING FORMATTED CORRECTLY?
                stop
            end if
            region_mask_projection_type = UTM_projection_index
            x_dim_name_regionmask='x'
            y_dim_name_regionmask='y'
            write(unit_logfile,'(A,f3.0)') 'Region mask file has projection UTM zone ',region_mask_projection_attributes(1)
        end if

        if (region_mask_projection_type == LL_projection_index) then
            x_dim_name_regionmask='lon'
            y_dim_name_regionmask='lat'
            write(unit_logfile, '(A)') 'Assuming region mask file is in lon-lat, since it is variable projection_utm was not found. Other projections than UTM are not implemented'
        end if

        !Read dimensions of the file (x,y)
        ! x
        status_nc = NF90_INQ_DIMID (id_nc,x_dim_name_regionmask,x_dim_id_nc)
        status_nc = NF90_INQUIRE_DIMENSION (id_nc,x_dim_id_nc,dimname_temp,nx_regionmask)
        if (status_nc /= NF90_NOERR) then
            write(unit_logfile,'(A,I0)') 'Reading of x dimension failed with status: ',status_nc
            stop
        end if
        ! y
        status_nc = NF90_INQ_DIMID (id_nc,y_dim_name_regionmask,y_dim_id_nc)
        status_nc = NF90_INQUIRE_DIMENSION (id_nc,y_dim_id_nc,dimname_temp,ny_regionmask)
        if (status_nc /= NF90_NOERR) then
            write(unit_logfile,'(A,I0)') 'Reading of y dimension failed with status: ',status_nc
            stop
        end if
        write(unit_logfile,'(A,I0,A,I0)') ' Size of dimensions (x,y): ',nx_regionmask,' ',ny_regionmask

        ! Verify the dimensions are at least 2x2
        if (.not. (nx_regionmask > 1 .and. ny_regionmask > 1)) then
            write(unit_logfile, '(A)') 'Dimensions of regionmask is not at least 2x2'
            stop
        end if

        ! Allocate memory for the region mask and its coordinates
        allocate(region_mask(nx_regionmask, ny_regionmask))
        allocate(x_values_regionmask(nx_regionmask))
        allocate(y_values_regionmask(ny_regionmask))
        write(unit_logfile,'(A)') 'Allocation done.'
        ! NB: should I set these arrays to 0 before they are read from file?

        ! Read coordinate values
        ! x
        status_nc = NF90_INQ_VARID (id_nc, trim(x_dim_name_regionmask), var_id_nc)
        if (status_nc == NF90_NOERR) then
            status_nc = NF90_GET_VAR (id_nc,var_id_nc,x_values_regionmask)
        else
            write(unit_logfile,'(A)') 'Error while reading x values from region mask'
            stop
        end if
        ! y
        status_nc = NF90_INQ_VARID (id_nc, trim(y_dim_name_regionmask), var_id_nc)
        if (status_nc == NF90_NOERR) then
            status_nc = NF90_GET_VAR (id_nc,var_id_nc,y_values_regionmask)
        else
            write(unit_logfile,'(A)') 'Error while reading y values from region mask'
            stop
        end if

        write(unit_logfile,'(A)') 'Reading x and y values done.'

        ! Determine grid spacing and verify it is constant
        dx_regionmask = x_values_regionmask(2) - x_values_regionmask(1)
        do i = 2, nx_regionmask
            if (.not. (x_values_regionmask(i) - x_values_regionmask(i-1) == dx_regionmask)) then
                write(unit_logfile,'(A)') 'Not constant spacing in x coordinate in region mask'
                stop
            end if
        end do
        dy_regionmask = y_values_regionmask(2) - y_values_regionmask(1)
        do i = 2, ny_regionmask
            if (.not. (y_values_regionmask(i) - y_values_regionmask(i-1) == dy_regionmask)) then
                write(unit_logfile,'(A)') 'Not constant spacing in y coordinate in region mask'
                stop
            end if
        end do

        write(unit_logfile,'(A)') 'Reading mask itself'

        ! Read the mask itself
        status_nc = NF90_INQ_VARID (id_nc, trim(nlreg_varname_region_mask), var_id_nc)
        if (status_nc == NF90_NOERR) then
            status_nc = NF90_INQUIRE_VARIABLE(id_nc, var_id_nc, ndims = temp_num_dims)
            ! NB: the following line fails if region_mask is bigger than ca. 1 million elements when runnning interactively, but not as a qsub job (fix by increasing 'ulimit -s')
            status_nc = NF90_GET_VAR (id_nc, var_id_nc, region_mask)
            !write(unit_logfile, '(A,I0)') 'region_mask(1,1)= ',region_mask(1,1)
            !write(unit_logfile, '(A,I0)') 'region_mask(nx,ny)= ',region_mask(nx_regionmask,ny_regionmask)
            write(unit_logfile,'(A,i3,A,2A,2i16)') ' Reading: ',temp_num_dims,' ',trim(nlreg_varname_region_mask),' (min, max): ',minval(region_mask),maxval(region_mask)
        else
            write(unit_logfile,'(A)') 'Could not read region mask values from file'
            stop
        end if

        ! Close the region mask netcdf file
        status_nc = NF90_CLOSE (id_nc)

        ! Allocate arrays for region masks on EMEP and uEMEP grids
        allocate(nlreg_subgrid_region_id(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index)))
        allocate(nlreg_EMEP_subsample_region_id(nlreg_n_subsamples_per_EMEP_grid,nlreg_n_subsamples_per_EMEP_grid,dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index)))
        ! initalize to -1 (no-region).
        ! All -1 values will be replaced with positive values for the subgrid, but not necessarily for the EMEP subsamples
        nlreg_subgrid_region_id = -1
        nlreg_EMEP_subsample_region_id = -1

        write(unit_logfile, '(A)') 'Calculating region mask for the target grid'

        ! Set region ID of each target subgrid
        do i = 1, subgrid_dim(x_dim_index)
            do j = 1, subgrid_dim(y_dim_index)
                ! calculate x- and y- position in the region mask projection
                call LL2PROJ(lon_subgrid(i,j),lat_subgrid(i,j),x_location,y_location,region_mask_projection_attributes,region_mask_projection_type)
                ! Determine index in the region mask grid for this location
                x_index = nint(1+(x_location-x_values_regionmask(1))/dx_regionmask)
                y_index = nint(1+(y_location-y_values_regionmask(1))/dy_regionmask)
                
                ! Check if this location is inside the region mask grid
                if (x_index >= 1 .and. x_index <= nx_regionmask .and. y_index >= 1 .and. y_index <= ny_regionmask) then
                    nlreg_subgrid_region_id(i,j) = region_mask(x_index,y_index)
                    ! verify that the region ID is positive
                    if (.not. (nlreg_subgrid_region_id(i,j) > 0)) then
                        write(unit_logfile, '(A)') 'ERROR: Non-positive region ID found at a subgrid. This is not allowed'
                        stop
                    end if
                else
                    ! this receptor location is not within the region mask grid
                    write(unit_logfile, '(A)') 'ERROR: The target subgrid extends outside the given region mask. This is not allowed'
                    stop
                end if
                !write(unit_logfile,'(A,2i12,2f12.4,2f12.2,2i12,i4)') 'i,j,lon_subgrid(i,j),lat_subgrid(i,j),x_location,y_location,x_index,y_index,region_id = ',i,j,lon_subgrid(i,j),lat_subgrid(i,j),x_location,y_location,x_index,y_index,nlreg_subgrid_region_id(i,j)
            end do
        end do

        ! Set region ID of each subsample of the EMEP grid

        write(unit_logfile,'(A)') 'Calculating region mask for the EMEP grid'

        ! determine spacing in EMEP grid (NB: maybe this is alredy availabe somewhere?)
        ! NB: I will not verify it is constant, but I assume it is
        dx_emep = var1d_nc(2,x_dim_nc_index) - var1d_nc(1,x_dim_nc_index)
        dy_emep = var1d_nc(2,y_dim_nc_index) - var1d_nc(1,y_dim_nc_index)
        !write(unit_logfile, '(A,2I12,2f12.2)') 'EMEP nx,ny,dx,dy = ',dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dx_emep,dy_emep
        ! loop over all EMEP grid cells
        do ii = 1, dim_length_nc(x_dim_nc_index)
            do jj = 1, dim_length_nc(y_dim_nc_index)
                ! EMEP projection coordinate values at centre of this EMEP grid
                x_emepmid = var1d_nc(ii,x_dim_nc_index)
                y_emepmid = var1d_nc(jj,y_dim_nc_index)
                !write(unit_logfile,'(A,2i12,2f12.2)') 'ii,jj,x_emepmid,yemepmid = ',ii,jj,x_emepmid,y_emepmid
                ! go through all subsamples of this EMEP grid
                do i_sub = 1, nlreg_n_subsamples_per_EMEP_grid
                    do j_sub = 1, nlreg_n_subsamples_per_EMEP_grid
                        ! EMEP projection coordinate value at this subsample of the EMEP grid
                        x_emepsub = x_emepmid - dx_emep/2 + (i_sub-0.5)*dx_emep/nlreg_n_subsamples_per_EMEP_grid
                        y_emepsub = y_emepmid - dy_emep/2 + (j_sub-0.5)*dy_emep/nlreg_n_subsamples_per_EMEP_grid
                        ! calculate longitude and latitude from the EMEP projection
                        call PROJ2LL(x_emepsub,y_emepsub,lon_emepsub,lat_emepsub,EMEP_projection_attributes,EMEP_projection_type)
                        ! calculate projection coordinates in the region mask grid
                        call LL2PROJ(lon_emepsub,lat_emepsub,x_location,y_location,region_mask_projection_attributes,region_mask_projection_type)
                        ! Determine index in the region mask grid for this location
                        x_index = nint(1+(x_location-x_values_regionmask(1))/dx_regionmask)
                        y_index = nint(1+(y_location-y_values_regionmask(1))/dy_regionmask)
                        ! Check if this location is inside the region mask grid
                        if (x_index >= 1 .and. x_index <= nx_regionmask .and. y_index >= 1 .and. y_index <= ny_regionmask) then
                            nlreg_EMEP_subsample_region_id(i_sub,j_sub,ii,jj) = region_mask(x_index,y_index)
                        end if
                        ! NB: if it is outside the region ID, the 'no-region' value -1 is kept

                        !if (ii > 1 .and. ii < 4 .and. jj > 3 .and. jj < 6) !then
                        !    write(unit_logfile,'(A,4i4,2f12.2,2f12.4,2f12.2,2i12,i12)') 'ii,jj,i_sub,j_sub,x_sub,y_sub,lon_sub,lat_sub,x_loc,y_loc,x_ind,y_ind,reg = ',ii,jj,i_sub,j_sub,x_emepsub,y_emepsub,lon_emepsub,lat_emepsub,x_location,y_location,x_index,y_index,nlreg_EMEP_subsample_region_id(i_sub,j_sub,ii,jj)
                        !end if
                    end do
                end do
            end do
        end do

        ! Deallocate the region mask
        deallocate(region_mask)
        deallocate(x_values_regionmask)
        deallocate(y_values_regionmask)

        write(unit_logfile,'(A)') 'Finding the regions occurring in the target grid'

        ! Determine which regions occur in the target grid
        ! (NB: We don't really need to worry about regions occurring only in the EMEP grid!)
        ! make room for 'n_regions_allocate' regions in the list initially
        allocate(temp_region_ids(n_regions_allocate))
        temp_region_ids = -1
        !write(unit_logfile,'(20i4)') temp_region_ids
        list_allocated_size = size(temp_region_ids)
        counter = 0
        previous_region_id = -1
        do i = 1, subgrid_dim(x_dim_index)
            do j = 1, subgrid_dim(y_dim_index)
                current_region_id = nlreg_subgrid_region_id(i,j)
                !write(unit_logfile,'(A,4i4)') 'i,j,cur,prev,counter = ',i,j,current_region_id,previous_region_id,counter
                if (current_region_id > 0 .and. .not. current_region_id == previous_region_id) then
                    ! Region ID is different from previous subgrid: check if we already found it before
                    current_region_already_found = .false.
                    do region_index = 1, counter
                        ! NB: What happens if counter=0 ??????????????????????????????
                        if (temp_region_ids(region_index) == current_region_id) then
                            current_region_already_found = .true.
                            exit
                        end if
                    end do
                    !write(unit_logfile,'(A,l)') 'already_found? ',current_region_already_found
                    if (.not. current_region_already_found) then
                        ! new region ID found
                        ! NB: region_id <= 0 is ignored!
                        if (counter == list_allocated_size) then
                            ! We need to allocate more space in the list of regions
                            !write(unit_logfile,'(A)') 'We need to allocate more!'
                            allocate(temp_region_ids_dummy(list_allocated_size))
                            temp_region_ids_dummy = temp_region_ids
                            deallocate(temp_region_ids)
                            allocate(temp_region_ids(list_allocated_size+n_regions_allocate))
                            temp_region_ids(1:list_allocated_size) = temp_region_ids_dummy
                            deallocate(temp_region_ids_dummy)
                            list_allocated_size = size(temp_region_ids)
                        end if
                        counter = counter + 1
                        temp_region_ids(counter) = current_region_id
                        !write(unit_logfile,'(20i4)') temp_region_ids
                    end if
                    previous_region_id = current_region_id
                end if
            end do
        end do
        nlreg_n_regions = counter
        allocate(nlreg_region_ids(nlreg_n_regions))
        nlreg_region_ids = temp_region_ids(1:nlreg_n_regions)
        deallocate(temp_region_ids)
        write(unit_logfile,'(A,I0)') 'Number of regions within target grid: ',nlreg_n_regions
        write(unit_logfile,'(A,100I5)') 'ID of these regions are (printing max 100): ', nlreg_region_ids
        
        write(unit_logfile,'(A)') 'Calculating fraction of each EMEP cell in each region'

        ! Calculate fraction of each EMEP grid that is within each region, by counting the subsamples
        allocate(nlreg_regionfraction_per_EMEP_grid(dim_length_nc(x_dim_nc_index), dim_length_nc(y_dim_nc_index),nlreg_n_regions))
        do ii = 1, dim_length_nc(x_dim_nc_index)
            do jj = 1, dim_length_nc(y_dim_nc_index)
                do region_index = 1, nlreg_n_regions
                    current_region_id = nlreg_region_ids(region_index)
                    counter = 0
                    do i_sub = 1, nlreg_n_subsamples_per_EMEP_grid
                        do j_sub = 1, nlreg_n_subsamples_per_EMEP_grid
                            if (nlreg_EMEP_subsample_region_id(i_sub,j_sub,ii,jj) == current_region_id) then
                                counter = counter + 1
                            end if
                        end do
                    end do
                    nlreg_regionfraction_per_EMEP_grid(ii,jj,region_index) = counter*1.0/nlreg_n_subsamples_per_EMEP_grid**2
                    !write(unit_logfile,'(A,6i5,f12.4)') 'ii,jj,region_index,region_id,counter,nsubsamples,fraction =',ii,jj,region_index,current_region_id,counter,nlreg_n_subsamples_per_EMEP_grid,nlreg_regionfraction_per_EMEP_grid(ii,jj,region_index)
                end do
            end do
        end do

    end subroutine nlreg_uEMEP_region_mask_new

end module auto_subgrid

