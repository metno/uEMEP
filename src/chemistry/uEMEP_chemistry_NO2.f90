module chemistry_no2

    use uemep_constants, only: pi, epsilon0
    use uemep_configuration
    use uEMEP_definitions
    use time_functions, only: get_sun_angles

    implicit none
    private

    public :: uEMEP_chemistry, correct_annual_mean_chemistry, &
        uEMEP_source_fraction_chemistry

contains

    ! THIS SUBROUTINE IS NO LONGER USED. uEMEP_chemistry IS CALLED DIRECTLY INSTEAD
    subroutine uEMEP_chemistry_control()
        real, allocatable :: subgrid_dummy(:,:,:,:,:,:)
        real, allocatable :: comp_subgrid_dummy(:,:,:,:)
        real, allocatable :: comp_source_subgrid_dummy(:,:,:,:,:)
        real, allocatable :: comp_source_EMEP_subgrid_dummy(:,:,:,:,:)
        real, allocatable :: comp_source_EMEP_additional_subgrid_dummy(:,:,:,:,:)

        integer :: in_region_loop, n_in_region_loop

        ! These are calculated in the Chemistry routine. Fist declared here. Are global variables
        if ( .not. allocated(comp_source_subgrid)) then
            allocate(comp_source_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index,n_source_index))
        end if
        if ( .not. allocated(comp_source_EMEP_subgrid)) then
            allocate(comp_source_EMEP_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index,n_source_index))
        end if
        if ( .not. allocated(comp_source_EMEP_additional_subgrid)) then
            allocate(comp_source_EMEP_additional_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index,n_source_index))
        end if
        if (trace_emissions_from_in_region) then
            n_in_region_loop = 2
            if ( .not. allocated(subgrid_dummy)) then
                allocate (subgrid_dummy(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_subgrid_index,n_source_index,n_pollutant_loop))
            end if
            if ( .not. allocated(comp_subgrid_dummy)) then
                allocate (comp_subgrid_dummy(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index))
            end if
            subgrid_dummy = 0
            comp_subgrid_dummy = 0
            if ( .not. allocated(comp_source_subgrid_dummy)) then
                allocate(comp_source_subgrid_dummy(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index,n_source_index))
            end if
            if ( .not. allocated(comp_source_EMEP_subgrid_dummy)) then
                allocate(comp_source_EMEP_subgrid_dummy(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index,n_source_index))
            end if
            if ( .not. allocated(comp_source_EMEP_additional_subgrid_dummy)) then
                allocate(comp_source_EMEP_additional_subgrid_dummy(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index,n_source_index))
            end if
            comp_source_subgrid_dummy = 0
            comp_source_EMEP_subgrid_dummy = 0
            comp_source_EMEP_additional_subgrid_dummy = 0
        else
            n_in_region_loop = 1
        end if

        do in_region_loop = 1, n_in_region_loop
            write(unit_logfile,'(a)')''
            write(unit_logfile,'(a)')'--------------------------'
            if (in_region_loop .eq. 1) write(unit_logfile,'(a)') 'Chemistry for all contributions'
            if (in_region_loop .eq. 2) write(unit_logfile,'(a)') 'Chemistry only for inside regional contributions'
            write(unit_logfile,'(a)')'--------------------------'

            ! Loop over the normal subgrid and the in region subgrid by saving  to a dummy variable and pputting back when finished
            if (trace_emissions_from_in_region.and.in_region_loop.eq.2) then
                subgrid_dummy = subgrid
                comp_subgrid_dummy = comp_subgrid
                subgrid = subgrid_from_in_region
                comp_subgrid = comp_subgrid_from_in_region

                ! These are calculated in the chemistry routine
                comp_source_subgrid_dummy = comp_source_subgrid
                comp_source_EMEP_subgrid_dummy = comp_source_EMEP_subgrid
                comp_source_EMEP_additional_subgrid_dummy = comp_source_EMEP_additional_subgrid
            end if

            call uEMEP_chemistry()

            if (trace_emissions_from_in_region .and. in_region_loop .eq. 2) then
                subgrid_from_in_region = subgrid
                comp_subgrid_from_in_region = comp_subgrid
                comp_source_subgrid_from_in_region = comp_source_subgrid
                comp_source_EMEP_subgrid_from_in_region = comp_source_EMEP_subgrid
                comp_source_EMEP_additional_subgrid_from_in_region = comp_source_EMEP_additional_subgrid

                subgrid = subgrid_dummy
                comp_subgrid = comp_subgrid_dummy
                comp_source_subgrid = comp_source_subgrid_dummy
                comp_source_EMEP_subgrid = comp_source_EMEP_subgrid_dummy
                comp_source_EMEP_additional_subgrid = comp_source_EMEP_additional_subgrid_dummy
            end if
        end do ! from_in_region loop

        if (trace_emissions_from_in_region) then
            if (allocated(subgrid_dummy)) deallocate(subgrid_dummy)
            if (allocated(comp_subgrid_dummy)) deallocate(comp_subgrid_dummy)
            if (allocated(comp_source_subgrid_dummy)) deallocate(comp_source_subgrid_dummy)
            if (allocated(comp_source_EMEP_subgrid_dummy)) deallocate(comp_source_EMEP_subgrid_dummy)
            if (allocated(comp_source_EMEP_additional_subgrid_dummy)) deallocate(comp_source_EMEP_additional_subgrid_dummy)
        end if

    end subroutine uEMEP_chemistry_control

    subroutine uEMEP_chemistry()
        ! Routine for doing the chemistry calculations in uEMEP
        integer :: i, j
        real :: nox_bg, no2_bg, o3_bg, nox_loc, f_no2_loc, J_photo, temperature
        real :: nox_out, no2_out, o3_out, p_bg_out, p_out
        integer :: t, t_start, t_end
        integer :: i_source, i_subsource, emep_subsource
        integer :: i_pollutant
        logical :: nox_available = .false.
        integer :: i_integral, j_integral
        real :: FF_loc, distance_grid
        integer :: i_cross_integral, j_cross_integral, i_nc, j_nc
        real :: sum_p_bg_out, sum_p_out, count_p_out
        real :: max_p_bg_out, max_p_out, min_p_bg_out, min_p_out

        ! NB. Additional is calculated but not necessarily saved!
        real :: nox_bg_additional, no2_bg_additional, o3_bg_additional

        ! These are calculated in the Chemistry routine. Fist declared here. Are global variables
        if ( .not. allocated(comp_source_subgrid)) then
            allocate(comp_source_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index,n_source_index))
        end if
        if ( .not. allocated(comp_source_EMEP_subgrid)) then
            allocate(comp_source_EMEP_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index,n_source_index))
        end if
        if ( .not. allocated(comp_source_EMEP_additional_subgrid)) then
            allocate(comp_source_EMEP_additional_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index,n_source_index))
        end if

        ! Search for nox in the pollutants
        do i_pollutant = 1, n_pollutant_loop
            if (pollutant_loop_index(i_pollutant) .eq. nox_nc_index) nox_available = .true.
        end do

        ! Leave the chemistry routine if nox is not available
        if ( .not. nox_available) return

        write(unit_logfile,'(A)') ''
        write(unit_logfile,'(A)') '================================================================'
        write(unit_logfile,'(A)') 'Calculating chemistry for NO2 (uEMEP_chemistry)'
        write(unit_logfile,'(A)') '================================================================'

        if (no2_chemistry_scheme_flag .eq. 0) then
            write(unit_logfile,'(A)') 'No chemistry used'
        else if (no2_chemistry_scheme_flag .eq. 1) then
            write(unit_logfile,'(A)') 'Photostationary state used'
        else if (no2_chemistry_scheme_flag .eq. 2) then
            write(unit_logfile,'(A)') 'Photochemistry with time scale used'
        else if (no2_chemistry_scheme_flag .eq. 3) then
            write(unit_logfile,'(A)') 'Romberg parameterisation used'
        else if (no2_chemistry_scheme_flag .eq. 4) then
            write(unit_logfile,'(A)') 'SRM parameterisation used'
        else if (no2_chemistry_scheme_flag .eq. 5) then
            write(unit_logfile,'(A)') 'During parameterisation used'
        end if

        t_start = 1
        t_end = subgrid_dim(t_dim_index)
        i_subsource = 1
        emep_subsource = 1
        comp_subgrid(:,:,:,no2_index) = 0
        comp_subgrid(:,:,:,nox_index) = 0
        comp_subgrid(:,:,:,o3_index) = 0

        nox_bg = 0.0; no2_bg = 0.0; o3_bg = 0.0; nox_loc = 0.0; f_no2_loc = 0.0; J_photo = 0.0; temperature=0.0

        ! Before calculating travel time then include the other EMEP sources not downscaled
        ! Travel time is set to EMEP Grid_width/FFgrid
        do t = t_start, t_end
            do j = 1, subgrid_dim(y_dim_index)
                do i = 1, subgrid_dim(x_dim_index)
                    i_cross_integral = crossreference_target_to_integral_subgrid(i,j,x_dim_index)
                    j_cross_integral = crossreference_target_to_integral_subgrid(i,j,y_dim_index)
                    FF_loc = 1.0
                    if (hourly_calculations) then
                        FF_loc = max(FF_min_dispersion, meteo_subgrid(i_cross_integral,j_cross_integral,t,FFgrid_subgrid_index))
                    else if (annual_calculations) then
                        FF_loc = max(FF_min_dispersion, 1.0/meteo_subgrid(i_cross_integral,j_cross_integral,t,inv_FFgrid_subgrid_index))
                    end if
                    i_nc = crossreference_target_to_emep_subgrid(i,j,x_dim_index)
                    j_nc = crossreference_target_to_emep_subgrid(i,j,y_dim_index)
                    if (EMEP_projection_type .eq. LL_projection_index) then
                        distance_grid = 111000.0*sqrt(dgrid_nc(lon_nc_index)*cos(var1d_nc(j_nc,lat_nc_index)*pi/180.0)*dgrid_nc(lat_nc_index))
                    else
                        ! Assumed LCC or PS
                        distance_grid = sqrt(dgrid_nc(lon_nc_index)*dgrid_nc(lat_nc_index))
                    end if
                end do
            end do
        end do

        sum_p_bg_out = 0.0
        sum_p_out = 0.0
        count_p_out = 0
        max_p_bg_out = -1000.0; min_p_bg_out = 1000.0; max_p_out = -1000.0; min_p_out = 1000.0
        do t = t_start, t_end
            do j = 1, subgrid_dim(y_dim_index)
                do i = 1, subgrid_dim(x_dim_index)
                    if (use_subgrid(i,j,allsource_index)) then
                        i_integral = crossreference_target_to_integral_subgrid(i,j,x_dim_index)
                        j_integral = crossreference_target_to_integral_subgrid(i,j,y_dim_index)
                        J_photo = meteo_subgrid(i_integral,j_integral,t,J_subgrid_index)
                        temperature = meteo_subgrid(i_integral,j_integral,t,t2m_subgrid_index)
                        nox_bg = subgrid(i,j,t,emep_nonlocal_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))

                        if (EMEP_additional_grid_interpolation_size .gt. 0) then
                            nox_bg_additional=subgrid(i,j,t,emep_additional_nonlocal_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))
                        end if

                        o3_bg = comp_EMEP_subgrid(i,j,t,o3_index)
                        f_no2_loc = 0.0
                        nox_loc = 0.0
                        do i_source = 1, n_source_index
                            if (calculate_source(i_source)) then
                                do i_subsource = 1, n_subsource(i_source)
                                    f_no2_loc = f_no2_loc + emission_factor(no2_index,i_source,i_subsource)/emission_factor(nox_index,i_source,i_subsource)*subgrid(i,j,t,local_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index))
                                    nox_loc = nox_loc + subgrid(i,j,t,local_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index))
                                end do
                            end if
                            if (calculate_emep_source(i_source) .and. .not. calculate_source(i_source)) then
                                do i_subsource = 1, n_subsource(i_source)
                                    f_no2_loc = f_no2_loc + f_no2_emep*subgrid(i,j,t,emep_local_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index))
                                    nox_loc = nox_loc + subgrid(i,j,t,emep_local_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index))
                                end do
                            end if
                        end do
                        if (abs(nox_loc) > epsilon0) then
                            f_no2_loc = f_no2_loc/nox_loc
                        else
                            f_no2_loc = 0.0
                        end if
                        no2_bg = comp_EMEP_subgrid(i,j,t,no2_index)*nox_bg/subgrid(i,j,t,emep_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))
                        o3_bg = max(0.0, comp_EMEP_subgrid(i,j,t,o3_index) + 48.0/46.0*(comp_EMEP_subgrid(i,j,t,no2_index) - no2_bg)) ! Conserve Ox when removing NO2 in the background. Cannot be less than 0

                        ! Assume stationary state to derive no2 and o3 background
                        if (no2_background_chemistry_scheme_flag .eq. 1) then
                            call uEMEP_nonlocal_NO2_O3(nox_bg, comp_EMEP_subgrid(i,j,t,nox_index), comp_EMEP_subgrid(i,j,t,no2_index), comp_EMEP_subgrid(i,j,t,o3_index), J_photo, temperature, f_no2_emep, no2_bg, o3_bg)
                        end if

                        if (EMEP_additional_grid_interpolation_size .gt. 0) then
                            no2_bg_additional = comp_EMEP_subgrid(i,j,t,no2_index)*nox_bg_additional/subgrid(i,j,t,emep_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))
                            if (no2_background_chemistry_scheme_flag .eq. 1) then
                                call uEMEP_nonlocal_NO2_O3(nox_bg_additional, comp_EMEP_subgrid(i,j,t,nox_index), comp_EMEP_subgrid(i,j,t,no2_index), comp_EMEP_subgrid(i,j,t,o3_index), J_photo, temperature, f_no2_emep, no2_bg_additional, o3_bg_additional)
                            else
                                o3_bg_additional = max(0.0, comp_EMEP_subgrid(i,j,t,o3_index) + 48.0/46.0*(comp_EMEP_subgrid(i,j,t,no2_index) - no2_bg_additional)) ! Conserve Ox when removing NO2 in the background
                            end if
                            comp_source_EMEP_additional_subgrid(i,j,t,o3_index,allsource_index) = o3_bg_additional
                            comp_source_EMEP_additional_subgrid(i,j,t,no2_index,allsource_index) = no2_bg_additional
                        end if

                        ! Set the background O3 level. use all_source for the nonlocal.
                        comp_source_EMEP_subgrid(i,j,t,o3_index,allsource_index) = o3_bg
                        comp_source_EMEP_subgrid(i,j,t,no2_index,allsource_index) = no2_bg

                        if (no2_chemistry_scheme_flag .eq. 0) then
                            nox_out = nox_bg + nox_loc
                            no2_out = no2_bg + nox_loc*f_no2_loc
                            o3_out = o3_bg
                        else if (no2_chemistry_scheme_flag .eq. 1) then
                            call uEMEP_photostationary_NO2(nox_bg, no2_bg, o3_bg, nox_loc, f_no2_loc, J_photo, temperature, nox_out, no2_out, o3_out, p_bg_out, p_out)
                        else if (no2_chemistry_scheme_flag .eq. 2) then
                            call uEMEP_phototimescale_NO2(nox_bg, no2_bg, o3_bg, nox_loc, f_no2_loc, J_photo, temperature, traveltime_subgrid(i,j,t,3,pollutant_loop_back_index(nox_nc_index))*traveltime_scaling, nox_out, no2_out, o3_out, p_bg_out, p_out)
                        else if (no2_chemistry_scheme_flag .eq. 3) then
                            call uEMEP_Romberg_NO2(nox_bg, no2_bg, nox_loc, o3_bg, f_no2_loc, nox_out, no2_out, o3_out, romberg_parameters)
                        else if (no2_chemistry_scheme_flag .eq. 4) then
                            call uEMEP_SRM_NO2(nox_bg, no2_bg, nox_loc, o3_bg, f_no2_loc, nox_out, no2_out, o3_out, SRM_parameters)
                        else if (no2_chemistry_scheme_flag .eq. 5) then
                            call uEMEP_During_NO2(nox_bg, no2_bg, nox_loc, o3_bg, f_no2_loc, comp_EMEP_subgrid(i,j,t,nox_index), comp_EMEP_subgrid(i,j,t,no2_index), comp_EMEP_subgrid(i,j,t,o3_index), J_photo, temperature, nox_out, no2_out, o3_out, p_bg_out, p_out)
                        end if

                        sum_p_bg_out = sum_p_bg_out + p_bg_out
                        sum_p_out = sum_p_out + p_out
                        count_p_out = count_p_out + 1
                        max_p_bg_out = max(max_p_bg_out, p_bg_out); min_p_bg_out = min(min_p_bg_out, p_bg_out)
                        max_p_out = max(max_p_out, p_out); min_p_out = min(min_p_out, p_out)
                        comp_subgrid(i,j,t,o3_index) = o3_out
                        comp_subgrid(i,j,t,no2_index) = no2_out
                        comp_subgrid(i,j,t,nox_index) = nox_out
                    else
                        comp_subgrid(i,j,t,o3_index) = NODATA_value
                        comp_subgrid(i,j,t,no2_index) = NODATA_value
                        comp_subgrid(i,j,t,nox_index) = NODATA_value

                    end if

                end do
            end do
        end do

        write(*,'(A,2f12.3)') 'P value (nonlocal,local): ', sum_p_bg_out/count_p_out, sum_p_out/count_p_out
        write(*,'(A,2f12.3)') 'P max (nonlocal,local): ', max_p_bg_out, max_p_out
        write(*,'(A,2f12.3)') 'P min (nonlocal,local): ', min_p_bg_out, min_p_out

    end subroutine uEMEP_chemistry

    subroutine uEMEP_source_fraction_chemistry()
        ! Special source allocation for no2 based on leaving out one source at a time in the chemistry calculation
        ! This will always give a sum less, but not much less than, the total no2
        ! This is normalised in order for it to be used
        ! Vhemistry scheme must have been run prior to implementing this
        integer :: i, j
        real :: nox_bg, no2_bg, o3_bg, nox_loc, f_no2_loc, J_photo, temperature
        real :: nox_out, no2_out, o3_out, p_bg_out, p_out
        integer :: t, t_start, t_end
        integer :: i_source, i_subsource, emep_subsource
        integer :: i_pollutant
        logical :: nox_available = .false.
        integer :: i_integral, j_integral
        integer :: remove_source
        real :: sum_no2_source_subgrid, sum_o3_source_subgrid
        real, allocatable :: comp_source_temp_subgrid(:,:,:,:,:)
        real, allocatable :: comp_source_EMEP_temp_subgrid(:,:,:,:,:)

        ! additional delarations needed for the in-region calculations
        integer, parameter :: inregion_index = 1
        integer, parameter :: outregion_index = 2
        integer :: k
        integer subgrid_var_index
        real :: f_no2_isource, nox_isource
        real :: no2_inandout_region(2)
        real :: o3_inandout_region(2)
        real :: sum_no2_inregion_outregion, sum_o3_inregion_outregion

        if (trace_emissions_from_in_region .and. .not. calculate_EMEP_additional_grid_flag) then
            if (.not. allocated(comp_source_subgrid_from_in_region)) allocate(comp_source_subgrid_from_in_region(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index,n_source_index))
        end if

        ! Search for nox in the pollutants
        do i_pollutant = 1, n_pollutant_loop
            if (pollutant_loop_index(i_pollutant) .eq. nox_nc_index) nox_available = .true.
        end do

        ! Leave the chemistry routine if nox is not available
        if ( .not. nox_available) return

        if ( .not. allocated(comp_source_subgrid)) then
            allocate(comp_source_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index,n_source_index))
        end if

        if (calculate_EMEP_additional_grid_flag) then
            if ( .not. allocated(comp_source_additional_subgrid)) then
                allocate(comp_source_additional_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index,n_source_index))
            end if
            ! Temporary array for storing the comp_source_subgrid to avoid rewriting large parts of the routine when running the additional version
            if ( .not. allocated(comp_source_temp_subgrid)) then
                allocate(comp_source_temp_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index,n_source_index))
            end if
            if ( .not. allocated(comp_source_EMEP_temp_subgrid)) then
                allocate(comp_source_EMEP_temp_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index,n_source_index))
            end if
            comp_source_temp_subgrid = comp_source_subgrid
            comp_source_EMEP_temp_subgrid = comp_source_EMEP_subgrid
            comp_source_EMEP_subgrid = comp_source_EMEP_additional_subgrid
        end if
        
        write(unit_logfile,'(A)') 'Calculating chemistry source contribution for NO2 (uEMEP_source_fraction_chemistry)'
        if (no2_chemistry_scheme_flag .eq. 0) then
            write(unit_logfile,'(A)') 'No chemistry used'
        else if (no2_chemistry_scheme_flag .eq. 1) then
            write(unit_logfile,'(A)') 'Photostationary state used'
        else if (no2_chemistry_scheme_flag .eq. 2) then
            write(unit_logfile,'(A)') 'Photochemistry with time scale used'
        else if (no2_chemistry_scheme_flag .eq. 3) then
            write(unit_logfile,'(A)') 'Romberg parameterisation used'
        else if (no2_chemistry_scheme_flag .eq. 4) then
            write(unit_logfile,'(A)') 'SRM parameterisation used'
        else if (no2_chemistry_scheme_flag .eq. 5) then
            write(unit_logfile,'(A)') 'During parameterisation used'
        end if

        t_start = 1
        t_end = subgrid_dim(t_dim_index)
        i_subsource = 1
        emep_subsource = 1

        nox_bg=0.0; no2_bg = 0.0; o3_bg = 0.0; nox_loc = 0.0; f_no2_loc = 0.0; J_photo = 0.0; temperature=0.0

        ! Weighted travel time already calculated
        do t = t_start, t_end
            do j = 1,subgrid_dim(y_dim_index)
                do i = 1,subgrid_dim(x_dim_index)
                    if (use_subgrid(i,j,allsource_index)) then
                        i_integral = crossreference_target_to_integral_subgrid(i,j,x_dim_index)
                        j_integral = crossreference_target_to_integral_subgrid(i,j,y_dim_index)
                        J_photo = meteo_subgrid(i_integral,j_integral,t,J_subgrid_index)
                        temperature = meteo_subgrid(i_integral,j_integral,t,t2m_subgrid_index)

                        if (calculate_EMEP_additional_grid_flag) then
                            nox_bg = subgrid(i,j,t,emep_additional_nonlocal_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))
                        else
                            nox_bg = subgrid(i,j,t,emep_nonlocal_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))
                        end if
                        o3_bg = comp_EMEP_subgrid(i,j,t,o3_index)

                        do remove_source = 1, n_source_index
                            if (calculate_source(remove_source) .or. remove_source .eq. allsource_index .or. (calculate_emep_source(remove_source) .and. .not. calculate_source(remove_source))) then
                                f_no2_loc = 0.0
                                nox_loc = 0.0

                                do i_source = 1, n_source_index
                                    if (calculate_source(i_source)) then
                                        if (remove_source .ne. i_source) then
                                            do i_subsource = 1, n_subsource(i_source)
                                                f_no2_loc = f_no2_loc + emission_factor(no2_index,i_source,i_subsource)/emission_factor(nox_index,i_source,i_subsource)*subgrid(i,j,t,local_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index))
                                                nox_loc = nox_loc + subgrid(i,j,t,local_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index))
                                            end do
                                        end if
                                    end if

                                    ! Include the local EMEP that are not being downscaled
                                    if ( .not. calculate_EMEP_additional_grid_flag) then
                                        if (calculate_emep_source(i_source) .and. .not. calculate_source(i_source)) then
                                            if (remove_source .ne. i_source) then
                                                do i_subsource = 1, n_subsource(i_source)
                                                    f_no2_loc = f_no2_loc + f_no2_emep*subgrid(i,j,t,emep_local_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index))
                                                    nox_loc = nox_loc + subgrid(i,j,t,emep_local_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index))
                                                end do
                                            end if
                                        end if
                                    end if

                                    if (calculate_EMEP_additional_grid_flag) then
                                        ! If calculating the additional region then use the additional local EMEP not being downscaled
                                        if (calculate_emep_source(i_source) .and. .not. calculate_source(i_source)) then
                                            if (remove_source .ne. i_source) then
                                                do i_subsource = 1, n_subsource(i_source)
                                                    f_no2_loc = f_no2_loc + f_no2_emep*subgrid(i,j,t,emep_additional_local_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index))
                                                    nox_loc = nox_loc + subgrid(i,j,t,emep_additional_local_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index))
                                                end do
                                            end if
                                        end if

                                        !If calculating the additional region then include the difference BG-BG_additional to the local EMEP that is being downscaled
                                        if (calculate_source(i_source)) then
                                            if (remove_source .ne. i_source) then
                                                do i_subsource = 1, n_subsource(i_source)
                                                    f_no2_loc = f_no2_loc + f_no2_emep* &
                                                        (subgrid(i,j,t,emep_nonlocal_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index)) &
                                                        - subgrid(i,j,t,emep_additional_nonlocal_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index)))
                                                    nox_loc = nox_loc + subgrid(i,j,t,emep_nonlocal_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index)) &
                                                        - subgrid(i,j,t,emep_additional_nonlocal_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index))
                                                end do
                                            end if
                                        end if
                                    end if
                                end do

                                if (abs(nox_loc) > epsilon0) then
                                    f_no2_loc = f_no2_loc/nox_loc
                                else
                                    f_no2_loc = 0.0
                                end if


                                ! Use the all source index to calculate the contribution from the background
                                ! This is done by removing all the sources, rather than the difference as done for the local sources
                                ! This is because the chemistry is disturbed when removing background nox and no2
                                if (remove_source .ne. allsource_index) then
                                    no2_bg = comp_EMEP_subgrid(i,j,t,no2_index)*nox_bg/subgrid(i,j,t,emep_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))
                                    o3_bg = max(0.0, comp_EMEP_subgrid(i,j,t,o3_index) + 48.0/46.0*(comp_EMEP_subgrid(i,j,t,no2_index) - no2_bg)) ! Conserve Ox when removing NO2 in the background. Cannot be less than 0
                                else
                                    no2_bg = comp_EMEP_subgrid(i,j,t,no2_index)*nox_bg/subgrid(i,j,t,emep_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))
                                    o3_bg = max(0.0, comp_EMEP_subgrid(i,j,t,o3_index) + 48.0/46.0*(comp_EMEP_subgrid(i,j,t,no2_index) - no2_bg)) !Conserve Ox when removing NO2 in the background. Cannot be less than 0
                                    nox_loc = 0.0
                                    f_no2_loc = 0.0
                                end if

                                ! Assume stationary state to derive no2 and o3 background. Overwrites the previous setting
                                if (no2_background_chemistry_scheme_flag .eq. 1) then
                                    call uEMEP_nonlocal_NO2_O3(nox_bg, comp_EMEP_subgrid(i,j,t,nox_index), comp_EMEP_subgrid(i,j,t,no2_index), comp_EMEP_subgrid(i,j,t,o3_index), J_photo, temperature, f_no2_emep, no2_bg, o3_bg)
                                end if

                                if (no2_chemistry_scheme_flag .eq. 0) then
                                    nox_out = nox_bg + nox_loc
                                    no2_out = no2_bg + nox_loc*f_no2_loc
                                    o3_out = o3_bg
                                else if (no2_chemistry_scheme_flag .eq. 1) then
                                    call uEMEP_photostationary_NO2(nox_bg, no2_bg, o3_bg, nox_loc, f_no2_loc, J_photo, temperature, nox_out, no2_out, o3_out, p_bg_out, p_out)
                                else if (no2_chemistry_scheme_flag .eq. 2) then
                                    call uEMEP_phototimescale_NO2(nox_bg, no2_bg, o3_bg, nox_loc, f_no2_loc, J_photo, temperature, traveltime_subgrid(i,j,t,3,pollutant_loop_back_index(nox_nc_index))*traveltime_scaling, nox_out, no2_out, o3_out, p_bg_out, p_out)
                                else if (no2_chemistry_scheme_flag .eq. 3) then
                                    call uEMEP_Romberg_NO2(nox_bg, no2_bg, nox_loc, o3_bg, f_no2_loc, nox_out, no2_out, o3_out, romberg_parameters)
                                else if (no2_chemistry_scheme_flag .eq. 4) then
                                    call uEMEP_SRM_NO2(nox_bg, no2_bg, nox_loc, o3_bg, f_no2_loc, nox_out, no2_out, o3_out, SRM_parameters)
                                else if (no2_chemistry_scheme_flag .eq. 5) then
                                    call uEMEP_During_NO2(nox_bg, no2_bg, nox_loc, o3_bg, f_no2_loc, comp_EMEP_subgrid(i,j,t,nox_index), comp_EMEP_subgrid(i,j,t,no2_index), comp_EMEP_subgrid(i,j,t,o3_index), J_photo, temperature, nox_out, no2_out, o3_out, p_bg_out, p_out)
                                end if

                                ! For background just use the result without any sources.
                                ! There is a problem disturbing the chemistry by removing the background nox and no2 but not changing the o3
                                if (remove_source .eq. allsource_index) then
                                    comp_source_subgrid(i,j,t,no2_index,remove_source) = no2_bg
                                    comp_source_subgrid(i,j,t,o3_index,remove_source) = o3_bg
                                else
                                    !Avoid round off errors which can occur with small numbers
                                    comp_source_subgrid(i,j,t,no2_index,remove_source) = max(0.0,comp_subgrid(i,j,t,no2_index) - no2_out)
                                    
                                    !Can be negative and can be greater than 1 so do not limit
                                    comp_source_subgrid(i,j,t,o3_index,remove_source) = comp_subgrid(i,j,t,o3_index) - o3_out
                                end if
                            end if
                        end do

                        !Normalise the contributions
                        !Calculate the sum
                        sum_no2_source_subgrid = 0.0
                        sum_o3_source_subgrid = 0.0
                        do i_source = 1, n_source_index
                            if (calculate_source(i_source) .or. (calculate_emep_source(i_source) .and. .not. calculate_source(i_source))) then
                                sum_no2_source_subgrid = sum_no2_source_subgrid + comp_source_subgrid(i,j,t,no2_index,i_source)
                                sum_o3_source_subgrid = sum_o3_source_subgrid + comp_source_subgrid(i,j,t,o3_index,i_source)
                            end if
                        end do

                        ! Set the background fractions so they will not be adjusted with normalisation
                        do i_source = 1, n_source_index
                            if (calculate_source(i_source) .or. (calculate_emep_source(i_source) .and. .not. calculate_source(i_source))) then
                                ! Adjust for the background and normalise
                                if (abs(sum_no2_source_subgrid) > epsilon0) then
                                    comp_source_subgrid(i,j,t,no2_index,i_source) = comp_source_subgrid(i,j,t,no2_index,i_source)/sum_no2_source_subgrid &
                                        *(comp_subgrid(i,j,t,no2_index) - comp_source_EMEP_subgrid(i,j,t,no2_index,allsource_index))
                                else
                                    comp_source_subgrid(i,j,t,no2_index,i_source) = 0
                                end if
                                
                                if (abs(sum_o3_source_subgrid) > epsilon0) then
                                    comp_source_subgrid(i,j,t,o3_index,i_source) = comp_source_subgrid(i,j,t,o3_index,i_source)/sum_o3_source_subgrid &
                                        *(comp_subgrid(i,j,t,o3_index) - comp_source_EMEP_subgrid(i,j,t,o3_index,allsource_index))
                                else
                                    comp_source_subgrid(i,j,t,o3_index,i_source) = 0
                                end if
                                ! Setting local sources to 0 if total concentration is zero: No longer do this, because nonlocal might be non-zero even if total is zero
                                !if (comp_subgrid(i,j,t,no2_index) .le. 0) comp_source_subgrid(i,j,t,no2_index,i_source) = 0
                                !if (comp_subgrid(i,j,t,o3_index) .le. 0) comp_source_subgrid(i,j,t,o3_index,i_source) = 0
                            end if
                        end do

                        ! Calculate NO2 and O3 source contributions from-in-region
                        ! ********************************************************
                        if (trace_emissions_from_in_region .and. .not. calculate_EMEP_additional_grid_flag) then
                            do remove_source = 1, n_source_index
                                if (calculate_source(remove_source) .or. calculate_EMEP_source(remove_source)) then
                                    do k = 1,2 ! inregion and outregion
                                        ! add up all local sources to NOx, except 'remove_source' from either inregion or outregion
                                        f_no2_loc = 0
                                        nox_loc = 0
                                        do i_source = 1, n_source_index
                                            if (calculate_source(i_source) .or. calculate_EMEP_source(i_source)) then
                                                ! NB: loop over n_subsource is not logical, since we in practice weigh contributions from sources with 2 subsources more than sources with only 1, but I do it to follow Bruce's method above...
                                                do i_subsource = 1, n_subsource(i_source)
                                                    ! check whether to use downscaled or non-downscaled local contribution for this source
                                                    if (calculate_source(i_source)) then
                                                        ! downscaled contribution
                                                        f_no2_isource = emission_factor(no2_index,i_source,i_subsource)/emission_factor(nox_index,i_source,i_subsource)
                                                        subgrid_var_index = local_subgrid_index
                                                    else ! i.e. calculate_EMEP_source(i_source) .and. .not. calculate_source(i_source)
                                                        ! EMEP contribution
                                                        f_no2_isource = f_no2_emep
                                                        subgrid_var_index = emep_local_subgrid_index
                                                    end if
                                                    ! check if this source is the one we remove or not, to determine how to add it
                                                    if (i_source == remove_source) then
                                                        ! for the source to remove: We should then treat the inregion and outregion differently
                                                        if (k == inregion_index) then
                                                            ! in region: add only the local contribution from outside region
                                                            nox_isource = subgrid(i,j,t,subgrid_var_index,i_source,pollutant_loop_back_index(nox_nc_index)) - subgrid_from_in_region(i,j,t,subgrid_var_index,i_source,pollutant_loop_back_index(nox_nc_index))
                                                        else if (k == outregion_index) then
                                                            ! out region: add only the local contribution from inside region
                                                            nox_isource = subgrid_from_in_region(i,j,t,subgrid_var_index,i_source,pollutant_loop_back_index(nox_nc_index))
                                                        else
                                                            write(unit_logfile,'(A)') 'ERROR: value of k is not inregion_index or outregion_index'
                                                            stop
                                                        end if
                                                    else
                                                        ! for other sources, just add all the local contribution in both cases
                                                        nox_isource = subgrid(i,j,t,subgrid_var_index,i_source,pollutant_loop_back_index(nox_nc_index))
                                                    end if
                                                    f_no2_loc = f_no2_loc + f_no2_isource*nox_isource
                                                    nox_loc = nox_loc + nox_isource
                                                end do
                                            end if
                                        end do
                                        ! divide f_no2_loc by total NOx contribution, in both cases
                                        if (abs(nox_loc) > epsilon0) then
                                            f_no2_loc = f_no2_loc/nox_loc
                                        else
                                            f_no2_loc = 0.0
                                        end if

                                        ! Calculate background concentrations (following Bruce's approach above)
                                        no2_bg = comp_EMEP_subgrid(i,j,t,no2_index)*nox_bg/subgrid(i,j,t,emep_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))
                                        o3_bg = max(0.0, comp_EMEP_subgrid(i,j,t,o3_index) + 48.0/46.0*(comp_EMEP_subgrid(i,j,t,no2_index) - no2_bg)) !Conserve Ox when removing NO2 in the background. Cannot be less than 0

                                        ! Assume stationary state to derive no2 and o3 background. Overwrites the previous setting
                                        if (no2_background_chemistry_scheme_flag .eq. 1) then
                                            call uEMEP_nonlocal_NO2_O3(nox_bg, comp_EMEP_subgrid(i,j,t,nox_index), comp_EMEP_subgrid(i,j,t,no2_index), comp_EMEP_subgrid(i,j,t,o3_index), J_photo, temperature, f_no2_emep, no2_bg, o3_bg)
                                        end if

                                        ! Calculate NO2 and O3 with the chemistry scheme
                                        if (no2_chemistry_scheme_flag .eq. 0) then
                                            nox_out = nox_bg + nox_loc
                                            no2_out = no2_bg + nox_loc*f_no2_loc
                                            o3_out = o3_bg
                                        else if (no2_chemistry_scheme_flag .eq. 1) then
                                            call uEMEP_photostationary_NO2(nox_bg, no2_bg, o3_bg, nox_loc, f_no2_loc, J_photo, temperature, nox_out, no2_out, o3_out, p_bg_out, p_out)
                                        else if (no2_chemistry_scheme_flag .eq. 2) then
                                            call uEMEP_phototimescale_NO2(nox_bg, no2_bg, o3_bg, nox_loc, f_no2_loc, J_photo, temperature, traveltime_subgrid(i,j,t,3,pollutant_loop_back_index(nox_nc_index))*traveltime_scaling, nox_out, no2_out, o3_out, p_bg_out, p_out)
                                        else if (no2_chemistry_scheme_flag .eq. 3) then
                                            call uEMEP_Romberg_NO2(nox_bg, no2_bg, nox_loc, o3_bg, f_no2_loc, nox_out, no2_out, o3_out, romberg_parameters)
                                        else if (no2_chemistry_scheme_flag .eq. 4) then
                                            call uEMEP_SRM_NO2(nox_bg, no2_bg, nox_loc, o3_bg, f_no2_loc, nox_out, no2_out, o3_out, SRM_parameters)
                                        else if (no2_chemistry_scheme_flag .eq. 5) then
                                            call uEMEP_During_NO2(nox_bg, no2_bg, nox_loc, o3_bg, f_no2_loc, comp_EMEP_subgrid(i,j,t,nox_index), comp_EMEP_subgrid(i,j,t,no2_index), comp_EMEP_subgrid(i,j,t,o3_index), J_photo, temperature, nox_out, no2_out, o3_out, p_bg_out, p_out)
                                        end if

                                        !Avoid round off errors which can occur with small numbers
                                        no2_inandout_region(k) = max(0.0, comp_subgrid(i,j,t,no2_index) - no2_out)
                                        !Can be negative and can be greater than 1 so do not limit
                                        o3_inandout_region(k) = comp_subgrid(i,j,t,o3_index) - o3_out
                                    end do ! k=1,2

                                    ! scale the contributions so the sum equals the total contribution from the source
                                    sum_no2_inregion_outregion = no2_inandout_region(inregion_index) + no2_inandout_region(outregion_index)
                                    sum_o3_inregion_outregion = o3_inandout_region(inregion_index) + o3_inandout_region(outregion_index)

                                    if (abs(sum_no2_inregion_outregion) > epsilon0) then
                                        comp_source_subgrid_from_in_region(i,j,t,no2_index,remove_source) = comp_source_subgrid(i,j,t,no2_index,remove_source) * no2_inandout_region(inregion_index) / sum_no2_inregion_outregion
                                    else
                                        comp_source_subgrid_from_in_region(i,j,t,no2_index,remove_source) = 0
                                    end if
                                    if (abs(sum_o3_inregion_outregion) > epsilon0) then
                                        comp_source_subgrid_from_in_region(i,j,t,o3_index,remove_source) = comp_source_subgrid(i,j,t,o3_index,remove_source) * o3_inandout_region(inregion_index) / sum_o3_inregion_outregion
                                    else
                                        comp_source_subgrid_from_in_region(i,j,t,o3_index,remove_source) = 0
                                    end if

                                end if
                            end do ! remove_source = 1, n_source_index
                        end if !(trace_emissions_from_in_region .and. .not. calculate_EMEP_additional_grid_flag)
                        ! ***************************************************************
                        ! done calculating NO2 and O3 source contributions from-in-region

                    else ! i.e. if (.not. (use_subgrid(i,j,allsource_index)))
                        comp_source_subgrid(i,j,t,:,:) = NODATA_value
                        if (trace_emissions_from_in_region .and. .not. calculate_EMEP_additional_grid_flag) then
                            comp_source_subgrid_from_in_region(i,j,t,:,:) = NODATA_value
                        end if
                    end if
                end do
            end do
        end do

        ! Transfer the arrays to the right outputs
        if (calculate_EMEP_additional_grid_flag) then
            comp_source_additional_subgrid = comp_source_subgrid
            comp_source_subgrid = comp_source_temp_subgrid
            comp_source_EMEP_subgrid = comp_source_EMEP_temp_subgrid
            
            ! EMEP_additional is unchanged
            if (allocated(comp_source_temp_subgrid)) deallocate(comp_source_temp_subgrid)
            if (allocated(comp_source_EMEP_temp_subgrid)) deallocate(comp_source_EMEP_temp_subgrid)
        end if

    end subroutine uEMEP_source_fraction_chemistry

    subroutine uEMEP_photostationary_NO2(nox_bg, no2_bg, o3_bg, nox_loc, f_no2_loc, J_photo, temperature, nox_out, no2_out, o3_out, p_bg_out, p_out)
        real, intent(in) :: nox_bg, no2_bg, o3_bg, nox_loc, f_no2_loc, J_photo, temperature
        real, intent(out) :: nox_out, no2_out, o3_out, p_bg_out, p_out
        integer no2_i, no_i, nox_i, o3_i, ox_i, nox_bg_i, no2_bg_i
        integer, parameter :: n_i = 7
        real :: Na, Na_fac, k1
        real :: mass(n_i)
        real :: mmass(n_i) = [46.0, 30.0, 46.0, 48.0, 47.0, 46.0, 46.0]
        real :: mol(n_i)
        real :: f_no2, f_ox, Jd, fac_sqrt
        real :: min_nox = 1.0e-6


        no2_i = 1; no_i = 2; nox_i = 3; o3_i = 4; ox_i = 5; nox_bg_i = 6; no2_bg_i = 7

        Na = 6.022e23        ! (molecules/mol)
        Na_fac = Na/1.0e12   !Conversion from ug/m3 to molecules/cm3 included
        k1 = 1.4e-12*exp(-1310.0/temperature); !(cm^3/s) and temperature in Kelvin
        mass(1:n_i) = 0.0
        mol(1:n_i) = 0.0

        ! Test for 0 NOx. If so leave the routine
        mass(nox_i) = nox_loc + nox_bg
        if (mass(nox_i) .le. min_nox) then
            nox_out = 0.0
            no2_out = 0.0
            o3_out = o3_bg
            return
        end if

        ! Check the photostationary assumption for the input data
        mass(nox_i) = nox_bg
        mass(no2_i) = no2_bg
        mass(o3_i) = o3_bg
        mol = mass/mmass*Na_fac ! (molecules per cm3)
        mol(ox_i) = mol(o3_i) + mol(no2_i)
        mol(no_i) = max(0.0, mol(nox_i) - mol(no2_i))
        
        ! Test the photostationary state for the bg input data
        if (abs(J_photo) > epsilon0 .and. abs(mol(no_i)) > epsilon0 .and. abs(mol(o3_i)) > epsilon0) then
            p_bg_out = J_photo*mol(no2_i)/k1/mol(o3_i)/mol(no_i)
        else
            p_bg_out = mol(no2_i)/(mol(ox_i) + mol(nox_i) - abs(mol(ox_i) - mol(nox_i)))*2.0
        end if

        ! Add the local contribution for calculation
        mass(nox_i) = nox_loc + nox_bg
        mass(no2_i) = f_no2_loc*nox_loc + no2_bg
        mass(o3_i) = o3_bg
        mol = mass/mmass*Na_fac ! (molecules per cm3)

        mol(ox_i) = mol(o3_i) + mol(no2_i)
        mol(no_i) = max(0.0, mol(nox_i) - mol(no2_i))

        f_no2 = mol(no2_i)/mol(nox_i)
        f_ox = mol(ox_i)/mol(nox_i)

        ! Test the photostationary state for the input data. Will not be in equilibrium
        if (abs(J_photo) > epsilon0 .and. abs(mol(no_i)) > epsilon0 .and. abs(mol(o3_i)) > epsilon0) then
            p_out = J_photo*mol(no2_i)/k1/mol(o3_i)/mol(no_i)
        else
            p_out = mol(no2_i)/(mol(ox_i)+mol(nox_i) - abs(mol(ox_i) - mol(nox_i)))*2.0
        end if

        ! Set the photolysis rate
        Jd = J_photo/k1/mol(nox_i)

        ! Calculate fraction of NO2 in photostationary state
        fac_sqrt = max(0.0, (1.0 + f_ox + Jd)**2 - 4.0*f_ox)
        f_no2 = 0.5*((1.0 + f_ox + Jd) - sqrt(fac_sqrt))

        ! Convert back to mass
        mol(no2_i) = f_no2*mol(nox_i);
        mol(o3_i) = max(0.0, mol(ox_i) - mol(no2_i))  ! Rounding errors possible
        mol(no_i) = max(0.0, mol(nox_i) - mol(no2_i)) !Rounding errors possible
        mass = mol*mmass/Na_fac ! (ug/m3)
        no2_out = mass(no2_i)
        nox_out = mass(nox_i)
        o3_out = mass(o3_i)

        ! Check output
        if (abs(J_photo) > epsilon0 .and. abs(mol(no_i)) > epsilon0 .and. abs(mol(o3_i)) > epsilon0) then
            p_out = J_photo*mol(no2_i)/k1/mol(o3_i)/mol(no_i)
        else
            p_out = mol(no2_i)/(mol(ox_i) + mol(nox_i) - abs(mol(ox_i) - mol(nox_i)))*2.0
        end if

    end subroutine uEMEP_photostationary_NO2

    subroutine uEMEP_phototimescale_NO2(nox_bg, no2_bg, o3_bg, nox_loc, f_no2_loc, J_photo, temperature, time_scale, nox_out, no2_out, o3_out, p_bg_out, p_out)
        real, intent(in) :: nox_bg, no2_bg, o3_bg, nox_loc, f_no2_loc, J_photo, temperature, time_scale
        real, intent(out) :: nox_out, no2_out, o3_out, p_bg_out, p_out
        
        ! Local variables
        integer :: no2_i,no_i,nox_i,o3_i,ox_i,nox_bg_i,no2_bg_i
        integer, parameter :: n_i = 7
        double precision :: Na, Na_fac, k1
        double precision :: mass(n_i)
        double precision :: mmass(n_i) = [46.0, 30.0, 46.0, 48.0, 47.0, 46.0, 46.0]
        double precision :: mol(n_i)
        double precision :: fac_sqrt
        double precision :: f_no2,f_ox,Jd,Jd_bg
        double precision :: min_nox=1.0e-6
        double precision :: c, b, BB, td, f_no2_0, f_no2_ps
        complex(8) :: AA
        double precision :: p_tot_out, f_ox_bg, f_no2_bg_ps, f_no2_bg

        no2_i = 1; no_i = 2; nox_i = 3; o3_i = 4; ox_i = 5; nox_bg_i = 6; no2_bg_i = 7
        Na = 6.022e23        ! (molecules/mol)
        Na_fac = Na/1.0e12   ! Conversion from ug/m3 to molecules/cm3 included
        k1 = 1.4e-12*exp(-1310.0/temperature) ! (cm^3/s) and temperature in Kelvin
        mass(1:n_i) = 0.0
        mol(1:n_i) = 0.0

        ! Test for 0 NOx. If so leave the routine
        mass(nox_i) = nox_loc + nox_bg
        if (mass(nox_i) .le. min_nox) then
            nox_out = 0.0
            no2_out = 0.0
            o3_out = o3_bg
            return
        end if

        ! Check the photostationary assumption for the input data
        mass(nox_i) = nox_bg
        mass(no2_i) = no2_bg
        mass(o3_i) = o3_bg
        mol = mass/mmass*Na_fac ! (molecules per cm3)

        mol(ox_i) = mol(o3_i) + mol(no2_i)
        mol(no_i) = max(0.0, mol(nox_i) - mol(no2_i))
        f_ox_bg = mol(ox_i)/mol(nox_i)
        Jd_bg = J_photo/k1/mol(nox_i)
        f_no2_bg_ps = 0.5*((1 + f_ox_bg + Jd_bg) - sqrt((1 + f_ox_bg+Jd_bg)**2 - 4.0*f_ox_bg))
        f_no2_bg = mol(no2_i)/mol(nox_i)
        p_bg_out = f_no2_bg/f_no2_bg_ps

        ! Check input
        if (abs(J_photo) > epsilon0 .and. abs(mol(no_i)) > epsilon0 .and. abs(mol(o3_i)) > epsilon0) then
            p_bg_out = J_photo*mol(no2_i)/k1/mol(o3_i)/mol(no_i)
        else
            p_bg_out = mol(no2_i)/(mol(ox_i) + mol(nox_i) - abs(mol(ox_i) - mol(nox_i)))*2.0
        end if

        ! Add the local contribution for calculation
        mass(nox_i) = nox_loc + nox_bg
        mass(no2_i) = f_no2_loc*nox_loc + no2_bg
        mass(o3_i) = o3_bg
        mol = mass/mmass*Na_fac ! (molecules per cm3)
        mol(ox_i) = mol(o3_i) + mol(no2_i)
        mol(no_i) = max(0.0, mol(nox_i) - mol(no2_i))
        f_no2 = mol(no2_i)/mol(nox_i)
        f_ox = mol(ox_i)/mol(nox_i)

        ! Set the photolysis rate
        Jd = J_photo/k1/mol(nox_i)

        ! Calculate photostationary for total nox, ox
        fac_sqrt = max(0.0, (1 + f_ox + Jd)**2 - 4.0*f_ox)
        f_no2_ps = 0.5*((1 + f_ox + Jd) - sqrt(fac_sqrt))
        p_tot_out = f_no2/f_no2_ps
        
        ! Calculate fraction of NO2 based on the time scale
        c = f_ox
        b = 1 + f_ox + Jd
        BB = sqrt(max(0.0, b**2 - 4.0*c))! max avoids roundoff errors
        td = time_scale*k1*mol(nox_i)
        f_no2_0 = f_no2

        AA = clog(cmplx((BB + b - 2.0*f_no2_0)/(BB - b + 2.0*f_no2_0)))
        f_no2 = real(-BB/2.0*((exp(AA + BB*td) - 1.0)/(exp(AA + BB*td) + 1.0)) + b/2.0)
        if (isnan(f_no2)) f_no2 = -BB/2.0 + b/2.0

        fac_sqrt = max(0.0, (1 + f_ox + Jd)**2 - 4.0*f_ox)
        f_no2_ps = 0.5*((1.0 + f_ox + Jd) - sqrt(fac_sqrt))
        p_out = f_no2/f_no2_ps

        ! Convert back to mass
        mol(no2_i) = max(0.0, f_no2*mol(nox_i))
        mol(o3_i) = max(0.0, mol(ox_i) - mol(no2_i))  ! Rounding errors possible
        mol(no_i) = max(0.0, mol(nox_i) - mol(no2_i)) ! Rounding errors possible
        mass = mol*mmass/Na_fac ! (ug/m3)
        no2_out = mass(no2_i)
        nox_out = mass(nox_i)
        o3_out = mass(o3_i)

        if (isnan(no2_out)) then
            write(*,'(8a12)') 'nox_bg', 'no2_bg', 'o3_bg', 'nox_loc', 'f_no2_loc', 'J_photo', 'temperature', 'time_scale'
            write(*,'(8es12.2)') nox_bg, no2_bg, o3_bg, nox_loc, f_no2_loc, J_photo, temperature, time_scale
            write(*,'(4a12)') 'f_no2', 'BB', 'b', 'b**2-4.*c'
            write(*,'(4es12.2)') f_no2, BB, b, b**2 - 4.0*c
            stop 1
        end if

        ! Check output
        if (abs(J_photo) > epsilon0 .and. abs(mol(no_i)) > epsilon0 .and. abs(mol(o3_i)) > epsilon0) then
            p_out = J_photo*mol(no2_i)/k1/mol(o3_i)/mol(no_i)
        else
            p_out = mol(no2_i)/(mol(ox_i) + mol(nox_i) - abs(mol(ox_i) - mol(nox_i)))*2.0
        end if

    end subroutine uEMEP_phototimescale_NO2

    subroutine uEMEP_Romberg_NO2(nox_bg, no2_bg, nox_loc, o3_bg, f_no2_loc, nox_out, no2_out, o3_out, romberg_params)
        ! From Norwegian obs fit
        ! real :: a_rom=20
        ! real :: b_rom=30
        ! From model fit
        ! real :: a_rom=30
        ! real :: b_rom=35
        ! real :: c_rom=0.20
        ! Gral values 30 35 0.18
        ! B�chlin and B�singer (2008) 29 35 0.217
        real, intent(in) :: nox_bg, no2_bg, nox_loc, o3_bg, f_no2_loc
        real, intent(in) :: romberg_params(3)
        real, intent(out) :: nox_out, no2_out, o3_out

        ! Local variables
        real :: a_rom = 30.0
        real :: b_rom = 35.0
        real :: c_rom = 0.20
        real :: ox_init, no2_init, no2_equ
        
        ! If available, use custom parameter values
        if (abs(romberg_params(1)) > epsilon0) then
            a_rom = romberg_params(1)
            b_rom = romberg_params(2)
            c_rom = romberg_params(3)
        end if

        nox_out = nox_bg + nox_loc
        no2_equ = a_rom*nox_bg/(nox_bg+b_rom) + nox_bg*c_rom
        no2_out = a_rom*nox_out/(nox_out+b_rom) + nox_out*c_rom
        no2_out = no2_out - no2_equ + no2_bg
        no2_out = max(no2_bg, no2_out)
        no2_init = no2_bg + f_no2_loc*nox_loc

        ! Small adjustments for molecular weights
        ox_init = no2_init*47.0/46.0 + o3_bg*47.0/48.0
        o3_out = ox_init*48.0/47.0 - no2_out*48.0/46.0
    end subroutine uEMEP_Romberg_NO2

    subroutine uEMEP_SRM_NO2(nox_bg, no2_bg, nox_loc, o3_bg, f_no2_loc, nox_out, no2_out, o3_out, SRM_params)
        ! From model fit
        ! real :: beta = 0.45
        ! real :: K = 30.0
        ! real :: F = 0.2
        ! From RIVM Briefrapport 2014-0109
        ! beta=1
        ! K=100
        ! F=0.2
        !
        ! Reference
        ! https://core.ac.uk/download/pdf/58774365.pdf
        real, intent(in) :: nox_bg, no2_bg, nox_loc, o3_bg, f_no2_loc
        real, intent(in) :: SRM_params(3)
        real, intent(out) :: nox_out, no2_out, o3_out

        ! Local variables
        real :: beta = 0.45
        real :: K = 30.0
        real :: F = 0.2
        real :: ox_init, no2_init

        ! If available, use custom parameter values
        if (abs(SRM_params(1)) > epsilon0) then
            beta = SRM_params(1)
            K = SRM_params(2)
            F = SRM_params(3)
        end if

        nox_out = nox_bg + nox_loc
        no2_out = no2_bg + beta*o3_bg*nox_loc/(nox_loc + K/(1 - F)) + F*nox_loc
        no2_init = no2_bg + f_no2_loc*nox_loc

        ! Small adjustments for molecular weights
        ox_init = no2_init*47.0/46.0 + o3_bg*47.0/48.0
        o3_out = ox_init*48.0/47.0 - no2_out*48.0/46.0
    end subroutine uEMEP_SRM_NO2

    subroutine uEMEP_During_NO2(nox_bg, no2_bg, nox_loc, o3_bg, f_no2_loc, nox_emep, no2_emep, o3_emep, J_photo, temperature, nox_out, no2_out, o3_out, p_bg_out, p_out)
        ! Reference
        ! D�ring, I., B�chlin, W., Ketzel, M., Baum, A., Friedrich, U., Wurzler, S., 2011.
        ! A new simplified NO/NO2 conversion model under consideration of direct NO2-emissions.
        ! Meteorol. Zeitschrift 20, 67�73. doi:10.1127/0941-2948/2011/0491
        !
        ! Improved Methodologies for NO2 Exposure Assessment in the EU, page 53
        ! https://ec.europa.eu/environment/air/pdf/NO2_Exposure_Final_Report.pdf
        real, intent(in) :: nox_bg, no2_bg, nox_loc, o3_bg, f_no2_loc, J_photo, temperature
        real, intent(in) :: nox_emep, no2_emep, o3_emep
        real, intent(out) :: nox_out, no2_out, o3_out
        real, intent(out) :: p_bg_out, p_out

        ! Local variables
        real :: mol_nox_bg, mol_no2_bg, mol_nox_loc, mol_o3_bg, mol_no2_loc, mol_ox_loc, mol_no_bg, mol_ox_bg
        real :: mol_no2_out, mol_o3_out, mol_no_out
        real :: mol_nox_emep, mol_no2_emep, mol_o3_emep, mol_no_emep, mol_ox_emep, p_emep_out
        real :: b, d, r, c, k1
        real :: Na, Na_fac
        integer :: no2_i, no_i, nox_i, o3_i, ox_i
        integer, parameter :: n_i = 5
        real :: mmass(n_i) = [46.0, 30.0, 46.0, 48.0, 47.0]

        no2_i = 1; no_i = 2; nox_i = 3; o3_i = 4; ox_i = 5
        k1 = 1.4e-12*exp(-1310.0/temperature) ! (cm^3/s) and temperature in Kelvin
        Na = 6.022e23 ! (molecules/mol)
        Na_fac = Na/1.0e12 ! Conversion from ug/m3 to molecules/cm3 included
        
        ! Normally multiplied by *Na_fac but not necessary as it is just a scaling
        mol_no2_bg = no2_bg/mmass(no2_i)
        mol_no2_loc = f_no2_loc*nox_loc/mmass(no2_i)
        mol_nox_bg = nox_bg/mmass(nox_i)
        mol_nox_loc = nox_loc/mmass(nox_i)
        mol_o3_bg = o3_bg/mmass(o3_i)
        mol_no_bg = (nox_bg - no2_bg)/mmass(nox_i)
        mol_ox_bg = mol_o3_bg + mol_no2_bg

        mol_o3_emep = o3_emep/mmass(o3_i)
        mol_nox_emep = nox_emep/mmass(nox_i)
        mol_no2_emep = no2_emep/mmass(no2_i)
        mol_no_emep = max(0.0, mol_nox_emep - mol_no2_emep)
        mol_ox_emep = mol_o3_emep + mol_no2_emep
        mol_ox_loc = mol_o3_bg + mol_no2_bg + mol_no2_loc

        if (mol_no2_emep .gt. 0) then
            r = mol_o3_emep*mol_no_emep/mol_no2_emep
        else
            r = 0.0
        end if

        b = mol_ox_loc + mol_nox_bg + mol_nox_loc + r
        c = max(0.0, b**2 - 4.0*mol_ox_loc*(mol_nox_bg + mol_nox_loc)) ! Should never be less than 0 but can be -0
        d = sqrt(c)
        mol_no2_out = (b - d)/2.0
        mol_o3_out = mol_ox_loc - mol_no2_out
        mol_no_out = mol_nox_bg + mol_nox_loc - mol_no2_out

        nox_out = nox_bg + nox_loc
        no2_out = mol_no2_out*mmass(no2_i)
        o3_out = mol_o3_out*mmass(o3_i)

        ! Not correct as it does not calculate the actual photostationary equation
        p_out = r
        p_bg_out = mol_o3_bg*mol_no_bg/mol_no2_bg

        ! Check output
        if (abs(J_photo) > epsilon0 .and. abs(mol_no_out) > epsilon0 .and. abs(mol_o3_out) > epsilon0) then
            p_out = J_photo*mol_no2_out/k1/mol_o3_out/mol_no_out/Na_fac
            p_emep_out = J_photo*mol_no2_emep/k1/mol_o3_emep/mol_no_emep/Na_fac
            p_bg_out = J_photo*mol_no2_bg/k1/mol_o3_bg/mol_no_bg/Na_fac
        else
            p_out = mol_no2_out/(mol_ox_loc + mol_nox_bg + mol_nox_loc - abs(mol_ox_loc - mol_nox_bg - mol_nox_loc))*2.0
            p_emep_out = mol_no2_emep/(mol_ox_emep + mol_nox_emep - abs(mol_ox_emep - mol_nox_emep))*2.0
            p_bg_out = mol_no2_bg/(mol_ox_bg + mol_nox_bg - abs(mol_ox_bg - mol_nox_bg))*2.0
        end if
    end subroutine uEMEP_During_NO2

    subroutine uEMEP_nonlocal_NO2_O3(nox_bg, nox_emep, no2_emep, o3_emep, J_photo, temperature, f_no2, no2_out, o3_out)
        real, intent(in) :: J_photo, temperature, f_no2
        real, intent(in) :: nox_bg
        real, intent(in) :: nox_emep, no2_emep, o3_emep
        real, intent(out) :: no2_out, o3_out

        ! Local variables
        real :: mol_nox_bg
        real :: mol_no2_out, mol_o3_out
        real :: mol_nox_emep, mol_no2_emep, mol_o3_emep, mol_ox_emep, mol_no_emep
        real :: b, d, r, c
        real :: Na, Na_fac, k1
        real :: p_phot, r_phot
        integer no2_i, no_i, nox_i, o3_i, ox_i
        integer, parameter :: n_i = 5
        real :: mmass(n_i) = [46.0, 30.0, 46.0, 48.0, 47.0]

        no2_i = 1; no_i = 2; nox_i = 3; o3_i = 4; ox_i = 5
        Na = 6.022e23 ! (molecules/mol)
        Na_fac = Na/1.0e12 ! Conversion from ug/m3 to molecules/cm3 included
        k1 = 1.4e-12*exp(-1310.0/temperature) ! (cm^3/s) and temperature in Kelvin
        
        ! Normally multiplied by *Na_fac but not necessary as it is just a scaling
        mol_o3_emep = o3_emep/mmass(o3_i)*Na_fac
        mol_no2_emep = no2_emep/mmass(no2_i)*Na_fac
        mol_nox_emep = nox_emep/mmass(nox_i)*Na_fac
        mol_ox_emep = mol_o3_emep + mol_no2_emep - f_no2*mol_nox_emep
        mol_nox_bg = nox_bg/mmass(nox_i)*Na_fac
        mol_no_emep = max(0.0, mol_nox_emep - mol_no2_emep)

        if (mol_no2_emep .gt. 0) then
            r = mol_o3_emep*mol_no_emep/mol_no2_emep
            p_phot = J_photo/k1*mol_no2_emep/mol_o3_emep/mol_no_emep
            r_phot = J_photo/k1/Na_fac
            r = r_phot*Na_fac

            b = mol_ox_emep + mol_nox_bg + r
            c = max(0.0, b**2 - 4.0*mol_ox_emep*mol_nox_bg) ! Should never be less than 0 but can be -0.0
            d = sqrt(c)
            mol_no2_out = (b - d)/2.0
            mol_o3_out = mol_ox_emep - mol_no2_out

            no2_out = max(0.0, mol_no2_out*mmass(no2_i)/Na_fac)
            o3_out = max(0.0, mol_o3_out*mmass(o3_i)/Na_fac)
        else
            no2_out = 0.0
            o3_out = o3_emep
        end if
    end subroutine uEMEP_nonlocal_NO2_O3

    subroutine correct_annual_mean_chemistry()
        integer :: i, j, t
        integer :: t_start, t_end
        integer :: i_integral, j_integral
        real :: o3_in, nox_in, no2_in, J_photo_in, temperature_in, lon_in, lat_in
        real :: ox_sigma_ratio_in, nox_sigma_ratio_in
        real :: o3_out, no2_out
        real :: sum_no2_in, sum_no2_out
        integer :: no2_count
        logical :: run_all_flag

        write(unit_logfile,'(A)') ''
        write(unit_logfile,'(A)') '================================================================'
        write(unit_logfile,'(A)') 'Correcting annual mean NO2 and O3 (correct_annual_mean_chemistry)'
        write(unit_logfile,'(A)') '================================================================'

        t_start = 1
        t_end = subgrid_dim(t_dim_index)

        sum_no2_in = 0.0
        sum_no2_out = 0.0
        no2_count = 0
        do t = t_start, t_end

            !Run the conversion routine once to get the Jd distribution which is saved. This is to save time as this is slow. Done in the centre of the grid
            if (quick_annual_mean_pdf_chemistry_correction) then
                run_all_flag = .false.
                i = subgrid_dim(x_dim_index)/2
                j = subgrid_dim(y_dim_index)/2
                o3_in = comp_subgrid(i,j,t,o3_index)
                no2_in = comp_subgrid(i,j,t,no2_index)
                nox_in = comp_subgrid(i,j,t,nox_index)
                i_integral = crossreference_target_to_integral_subgrid(i,j,x_dim_index)
                j_integral = crossreference_target_to_integral_subgrid(i,j,y_dim_index)
                J_photo_in = meteo_subgrid(i_integral,j_integral,t,J_subgrid_index)
                temperature_in = meteo_subgrid(i_integral,j_integral,t,t2m_subgrid_index)
                lon_in = lon_subgrid(i,j)
                lat_in = lat_subgrid(i,j)
                ox_sigma_ratio_in = ox_sigma_ratio_pdf
                nox_sigma_ratio_in = nox_sigma_ratio_pdf
                call uEMEP_annual_mean_pdf_correction_NO2_O3(min_bin_pdf, max_bin_pdf, log10_step_bin_pdf, .true., no2_in, nox_in, o3_in, J_photo_in, &
                    temperature_in, ox_sigma_ratio_in, nox_sigma_ratio_in, lon_in, lat_in, no2_out, o3_out)
            else
                run_all_flag = .true.
            end if

            do j = 1, subgrid_dim(y_dim_index)
                do i = 1, subgrid_dim(x_dim_index)

                    o3_in = comp_subgrid(i,j,t,o3_index)
                    no2_in = comp_subgrid(i,j,t,no2_index)
                    nox_in = comp_subgrid(i,j,t,nox_index)
                    i_integral = crossreference_target_to_integral_subgrid(i,j,x_dim_index)
                    j_integral = crossreference_target_to_integral_subgrid(i,j,y_dim_index)
                    J_photo_in = meteo_subgrid(i_integral,j_integral,t,J_subgrid_index)
                    temperature_in = meteo_subgrid(i_integral,j_integral,t,t2m_subgrid_index)
                    lon_in = lon_subgrid(i,j)
                    lat_in = lat_subgrid(i,j)

                    ox_sigma_ratio_in = ox_sigma_ratio_pdf
                    nox_sigma_ratio_in = nox_sigma_ratio_pdf

                    if (o3_in .le. 0 .or. nox_in .le. 0 .or. no2_in .le. 0) then
                        o3_out = o3_in
                        no2_out = no2_in
                    else
                        call uEMEP_annual_mean_pdf_correction_NO2_O3(min_bin_pdf, max_bin_pdf, log10_step_bin_pdf, run_all_flag, no2_in, nox_in, o3_in, J_photo_in, &
                            temperature_in, ox_sigma_ratio_in, nox_sigma_ratio_in, lon_in, lat_in, no2_out, o3_out)
                    end if

                    comp_subgrid(i,j,t,o3_index) = o3_out
                    comp_subgrid(i,j,t,no2_index) = no2_out
                    sum_no2_in = sum_no2_in + no2_in
                    sum_no2_out = sum_no2_out + no2_out
                    no2_count = no2_count + 1

                    if (isnan(no2_out)) then
                        write(unit_logfile,'(a,2i5,4f12.4)') 'NaN in pdf output. Stopping ',i,j,no2_in,no2_out,o3_in,o3_out
                        stop 1
                    end if
                end do
            end do
        end do

        write(unit_logfile,'(a,f12.4)') 'Average NO2 scaling with pdf correction = ',sum_no2_out/sum_no2_in
    end subroutine correct_annual_mean_chemistry

    subroutine uEMEP_annual_mean_pdf_correction_NO2_O3(bin_min, bin_max, delta_log10_bin, run_all, no2_in, nox_in, o3_in, J_photo_in, temperature_in, ox_sigma_ratio_in, nox_sigma_ratio_in, lon_in, lat_in, no2_out, o3_out)
        real, intent(in) :: bin_min, bin_max, delta_log10_bin
        logical, intent(in) :: run_all
        real, intent(in) :: J_photo_in, temperature_in
        real, intent(in) :: no2_in, nox_in, o3_in
        real, intent(in) :: ox_sigma_ratio_in, nox_sigma_ratio_in
        real, intent(in) :: lon_in, lat_in
        real, intent(out) :: no2_out, o3_out

        ! Loca variables
        real :: mol_no2_out, mol_o3_out
        real :: Na, Na_fac, k1
        integer, parameter :: no2_i = 1
        !integer, parameter :: no_i = 2
        integer, parameter :: nox_i = 3
        integer, parameter :: o3_i = 4
        !integer, parameter :: ox_i = 5
        integer, parameter :: n_i = 5
        real :: mmass(n_i) = [46.0, 30.0, 46.0, 48.0, 47.0]
        real :: log10_bin_max, log10_bin_min
        real, allocatable :: log10_bin(:), bin(:), delta_bin(:)
        integer :: n_bin
        real :: nox_sigma_ratio, ox_sigma_ratio
        real :: ox_sigma, ox_sig_sqr, ox_mu, nox_sigma, nox_sig_sqr, nox_mu
        real, allocatable :: y_nox(:), y_ox(:), y_Jd(:)
        integer, parameter :: nbin_Jd = 8760 ! Hours in a year
        integer :: bin_temp
        real, save :: y_Jd_acc(nbin_Jd)
        real :: mean_y_Jd_acc, y_Jd_acc_temp, y_Jd_acc_temp_log10
        real :: y_all_val, y_all_prob, y_all_sum, y_all_prob_sum, y_all, y_annual, y_scale
        real :: azimuth_ang, zenith_ang
        real :: mol_nox, mol_no2, mol_ox, mol_o3, Jd
        double precision :: date_num
        integer :: i, j, k, l
        integer :: date_array(6)
        real :: min_sigma_ratio = 0.01
        real :: bin_temp2

        Na = 6.022e23 ! (molecules/mol)
        Na_fac = Na/1.0e12 ! Conversion from ug/m3 to molecules/cm3 included
        k1 = 1.4e-12*exp(-1310.0/temperature_in) ! (cm^3/s) and temperature in Kelvin
        mol_nox = nox_in*Na_fac/mmass(nox_i)
        mol_no2 = no2_in*Na_fac/mmass(no2_i)
        mol_o3 = o3_in*Na_fac/mmass(o3_i)
        mol_ox = mol_no2 + mol_o3
        Jd = J_photo_in/k1
        no2_out = no2_in
        o3_out = o3_in

        !Create the bins for the pdf in (mol/cm3). ox, nox and J
        log10_bin_max = log10(bin_max)
        log10_bin_min = log10(bin_min)
        n_bin = int((log10_bin_max - log10_bin_min)/delta_log10_bin) + 1
        ! Creates 80 bins with these settings

        allocate (log10_bin(n_bin))
        allocate (bin(n_bin))
        allocate (delta_bin(n_bin))

        do i = 1, n_bin
            log10_bin(i) = log10_bin_min + (i - 1)*delta_log10_bin
        end do
        do i = 1, n_bin
            bin(i) = (10.0**log10_bin(i))*Na_fac/mmass(nox_i)
            delta_bin(i) = (10.0**(log10_bin(i) + delta_log10_bin/2.0) - 10.0**(log10_bin(i) - delta_log10_bin/2.0))*Na_fac/mmass(nox_i)
        end do

        ! Distribute concentrations into the pdf bins
        ! Minimum needed to avoid NaNs in the calculation
        nox_sigma_ratio = 1.14
        ox_sigma_ratio = 0.21
        if (abs(nox_sigma_ratio_in) > epsilon0) nox_sigma_ratio = max(nox_sigma_ratio_in, min_sigma_ratio)
        if (abs(ox_sigma_ratio_in) > epsilon0) ox_sigma_ratio = max(ox_sigma_ratio_in, min_sigma_ratio)

        ox_sigma = mol_ox*ox_sigma_ratio
        ox_sig_sqr = log(1.0 + ox_sigma**2.0/mol_ox**2.0)
        ox_mu = log(mol_ox**2.0/sqrt(mol_ox**2.0 + ox_sigma**2.0))
        nox_sigma = mol_nox * nox_sigma_ratio
        nox_sig_sqr = log(1.0 + nox_sigma**2.0/mol_nox**2.0)
        nox_mu = log(mol_nox**2.0/sqrt(mol_nox**2.0 + nox_sigma**2.0))

        allocate (y_ox(n_bin))
        allocate (y_nox(n_bin))
        if ( .not. allocated(y_Jd)) allocate (y_Jd(n_bin))
        y_Jd = 0.0
        y_ox = 0.0
        y_nox = 0.0

        do i = 1, n_bin
            y_ox(i) = 1.0/bin(i)/sqrt(ox_sig_sqr)/sqrt(2.0*pi)*exp(-((log(bin(i)) - ox_mu)**2)/2.0/ox_sig_sqr)*delta_bin(i)
            y_nox(i) = 1.0/bin(i)/sqrt(nox_sig_sqr)/sqrt(2.0*pi)*exp(-((log(bin(i)) - nox_mu)**2)/2.0/nox_sig_sqr)*delta_bin(i)
        end do

        ! Create the Jd_acc distribution by looping through every hour in the year and extracting the zenith angle
        ! Only do this if requested for the first time
        if (run_all) then
            y_Jd_acc = 0
            do i = 1, nbin_Jd
                date_num = 1.0 + i/24.0
                date_array = 0
                zenith_ang = 0.0
                call get_sun_angles(lat_in, lon_in, date_array, date_num, 0.0, azimuth_ang, zenith_ang)
                if (zenith_ang .ge. 90.0) then
                    y_Jd_acc(i) = 0
                else
                    y_Jd_acc(i) = ((cosd(zenith_ang))**0.28)
                end if
            end do
        end if

        mean_y_Jd_acc = sum(y_Jd_acc)/nbin_Jd
        do i = 1, nbin_Jd
            y_Jd_acc_temp = Jd*y_Jd_acc(i)/mean_y_Jd_acc
            if (y_Jd_acc_temp/Na_fac*mmass(nox_i) < bin_min) then
                bin_temp = 1
            else
                y_Jd_acc_temp_log10 = log10(y_Jd_acc_temp/Na_fac*mmass(nox_i))
                bin_temp = floor((y_Jd_acc_temp_log10 - log10_bin_min)/delta_log10_bin + 0.5) + 1
                bin_temp = min(max(bin_temp, 1), n_bin)
            end if
            y_Jd(bin_temp) = y_Jd(bin_temp) + 1
        end do

        ! Normalise all distributions
        y_Jd = y_Jd/sum(y_Jd)
        y_ox = y_ox/sum(y_ox)
        y_nox = y_nox/sum(y_nox)

        !Calculate scaling based on photostationary assumption
        y_all_sum = 0
        y_all_prob_sum = 0

        do j = 1, n_bin
            do k = 1, n_bin
                do l = 1, n_bin
                    ! Calculate weighting
                    y_all_prob = y_nox(j)*y_Jd(k)*y_ox(l)
                    if (y_all_prob .gt. 0.0) then
                        ! Calculate NO2 value
                        bin_temp2 = bin(j) + bin(k) + bin(l)
                        y_all_val = (bin_temp2 - sqrt(bin_temp2*bin_temp2 - 4.0*bin(j)*bin(l)))/2.0
                        ! Add weighted value
                        y_all_sum = y_all_sum + y_all_val*y_all_prob
                        ! Calculate sum of weights for normalisation later
                        y_all_prob_sum = y_all_prob_sum + y_all_prob
                    end if
                end do
            end do
        end do
        y_all = y_all_sum/y_all_prob_sum

        ! Calculate the mean value
        y_annual = ((mol_nox + Jd + mol_ox) - sqrt((mol_nox + Jd + mol_ox)**2 - 4.0*mol_nox*mol_ox))/2.0
        y_scale = y_all/y_annual
        mol_no2_out = y_scale*mol_no2
        mol_o3_out = mol_ox - mol_no2_out
        no2_out = mol_no2_out/Na_fac*mmass(no2_i)
        o3_out = mol_o3_out/Na_fac*mmass(o3_i)

        deallocate (log10_bin)
        deallocate (bin)
        deallocate (delta_bin)
        deallocate (y_ox)
        deallocate (y_nox)
        deallocate (y_Jd)
    end subroutine uEMEP_annual_mean_pdf_correction_NO2_O3

end module chemistry_no2

