module read_time_profiles

    use time_functions, only: day_of_week, summer_time_europe, number_to_date

    implicit none
    private

    public :: uEMEP_read_time_profiles

contains

!uEMEP_read_time_profiles.f90

    subroutine uEMEP_read_time_profiles

        use uEMEP_definitions

        implicit none

        integer i,j
        character(256) temp_str
        integer unit_in
        integer exists
        integer week_day_temp
        double precision date_num_temp

        integer n_col
        !parameter (n_col=5)
        character(256), allocatable :: header_str(:)
        integer, allocatable :: source_index_in(:)
        integer temp_region_id
        integer n_hours_in_week,n_months_in_year
        integer, allocatable :: time_month_of_year_input(:),time_hour_of_week_input(:)
        real, allocatable :: val_month_of_year_input(:,:),val_hour_of_week_input(:,:)
        integer date_array(6)
        integer i_source,t,hour_of_week_index
        integer t_profile_loop
        integer emission_time_shift_temp

        integer i_cross,j_cross
        real hdd_temp
        integer col_val,index_val
        logical do_not_calculate_RWC_emissions

        !double precision date_to_number

        !emission_time_profile_subgrid=1.

        !Only read data if flag is correct
        if (local_subgrid_method_flag.ne.2) return

        write(unit_logfile,'(A)') ''
        write(unit_logfile,'(A)') '================================================================'
        write(unit_logfile,'(A)') 'Reading time profiles (uEMEP_read_time_profiles)'
        write(unit_logfile,'(A)') '================================================================'

        pathfilename_timeprofile=trim(pathname_timeprofile)//trim(filename_timeprofile)

        !Test existence of the filename. If does not exist then use default
        inquire(file=trim(pathfilename_timeprofile),exist=exists)
        if (.not.exists) then
            write(*,'(A,A)') ' ERROR: Time profile data file does not exist: ', trim(pathfilename_timeprofile)
            stop
        endif


        !Open the file for reading
        unit_in=20
        open(unit_in,file=pathfilename_timeprofile,access='sequential',status='old',readonly)
        write(unit_logfile,'(a)') ' Opening time profile file: '//trim(pathfilename_timeprofile)

        !Find the number of commas in the header to determine columns
        read(unit_in,'(a)') temp_str
        !write(*,*) trim(temp_str),index(temp_str,',')
        col_val=0
        do while (index(temp_str,',').ne.0)
            index_val=index(temp_str,',')
            temp_str=temp_str(index_val+1:)
            col_val=col_val+1
            !write(*,*) col_val,trim(temp_str)
        enddo

        n_col=col_val+1
        allocate(header_str(n_col))
        allocate(source_index_in(n_col-1))
        write(unit_logfile,'(a,i)') ' Number of columns read: ',n_col

        !stop

        rewind(unit_in)

        !Read source header string
        read(unit_in,*) header_str
        !write(unit_logfile,*) header_str
        !Read source index, correpsonding to uEMEP source indexes (change to SNAP or NFR later)
        read(unit_in,*) temp_str,source_index_in
        write(unit_logfile,*) trim(temp_str),source_index_in

        !Read region
        read(unit_in,*) temp_str,temp_region_id
        write(unit_logfile,'(a,i)') trim(temp_str),temp_region_id

        !Read Hour_of_week
        read(unit_in,*) temp_str,n_hours_in_week
        write(unit_logfile,'(a,i)') trim(temp_str),n_hours_in_week


        !write(*,*) num_week_traffic,days_in_week,hours_in_day,n_roadlinks
        if (.not.allocated(time_hour_of_week_input)) allocate (time_hour_of_week_input(n_hours_in_week))
        if (.not.allocated(val_hour_of_week_input)) allocate (val_hour_of_week_input(n_hours_in_week,n_col-1))

        do i=1,n_hours_in_week
            read(unit_in,*) time_hour_of_week_input(i),val_hour_of_week_input(i,1:n_col-1)
            !write(*,*) time_hour_of_week_input(i),val_hour_of_week_input(i,1:n_col-1)
        enddo

        !Read Hour_of_week
        read(unit_in,*) temp_str,n_months_in_year
        write(unit_logfile,'(a,i)') trim(temp_str),n_months_in_year


        !write(*,*) num_week_traffic,days_in_week,hours_in_day,n_roadlinks
        if (.not.allocated(time_month_of_year_input)) allocate (time_month_of_year_input(n_months_in_year))
        if (.not.allocated(val_month_of_year_input)) allocate (val_month_of_year_input(n_months_in_year,n_col-1))

        do i=1,n_months_in_year
            read(unit_in,*) time_month_of_year_input(i),val_month_of_year_input(i,1:n_col-1)
            !write(*,*) time_month_of_year_input(i),val_month_of_year_input(i,1:n_col-1)
        enddo

        close(unit_in,status='keep')
        !close(unit_in)

        !Normalise the data to be used with average emmissions to give an hourly average value and a monthly average value
        do i_source=1,n_col-1
            val_hour_of_week_input(:,i_source)=val_hour_of_week_input(:,i_source)/sum(val_hour_of_week_input(:,i_source))*n_hours_in_week
            val_month_of_year_input(:,i_source)=val_month_of_year_input(:,i_source)/sum(val_month_of_year_input(:,i_source))*n_months_in_year
        enddo
        !do i=1,n_hours_in_week
        !    write(*,*) time_hour_of_week_input(i),val_hour_of_week_input(i,1:n_col-1)
        !enddo
        !stop
        !Get time information for the current calculation
        !    if (use_single_time_loop_flag) then
        !        t_profile_loop=end_time_loop_index
        !    else
        !        t_profile_loop=dim_length_nc(time_dim_nc_index)
        !    endif

        t_profile_loop=dim_length_nc(time_dim_nc_index)

        do t=1,t_profile_loop
            !if (t.eq.t_loop) then
            !EMEP date is days since 1900
            !write(*,*) val_dim_nc(t,time_dim_nc_index)
            !Round up the hour. Not here, is done earlier now
            !date_num_temp=dble(ceiling(val_dim_nc(t,time_dim_nc_index)*24.))/24.
            date_num_temp=val_dim_nc(t,time_dim_nc_index)
            !write(*,*) real(ceiling(val_dim_nc(t,time_dim_nc_index)*24.)),date_num_temp

            call number_to_date(date_num_temp,date_array,ref_year_EMEP)
            if (t.eq.1) write(unit_logfile,'(a,6i6)') 'Date array start = ',date_array
            if (t.eq.t_profile_loop) write(unit_logfile,'(a,6i6)') 'Date array end = ',date_array
            week_day_temp= day_of_week(date_array)
            !write(unit_logfile,*) 'Day of week = ',week_day_temp
            !write(*,*) t,date_array

            if (summer_time_europe(date_array)) then
                if (auto_adjustment_for_summertime) then
                    emission_time_shift_temp=emission_timeprofile_hour_shift+1
                    if (t.eq.1) write(unit_logfile,'(a)') ' Emission profiles set to summer time. '
                else
                    emission_time_shift_temp=emission_timeprofile_hour_shift
                    if (t.eq.1) write(unit_logfile,'(a)') ' Emission profiles not adjusted for summer time. '
                endif
            else
                emission_time_shift_temp=emission_timeprofile_hour_shift
                if (t.eq.1) write(unit_logfile,'(a)') ' Emission profiles set to winter time. '
            endif

            hour_of_week_index=(week_day_temp-1)*24+date_array(4)+emission_time_shift_temp
            if (hour_of_week_index.gt.n_hours_in_week) hour_of_week_index=hour_of_week_index-n_hours_in_week
            if (hour_of_week_index.lt.1) hour_of_week_index=hour_of_week_index+n_hours_in_week

            do i_source=1,n_col-1
                if (calculate_source(source_index_in(i_source))) then
                    if (source_index_in(i_source).eq.shipping_index.and.read_monthly_and_daily_shipping_data_flag) then
                        !Do nothing as the time profile has already been set in uEMEP_read_monthly_and_daily_shipping_asi_data subroutine
                        !write(unit_logfile,'(a,i,es16.6)') 'Not resetting shipping time profiles: ',t,sum(emission_time_profile_subgrid(:,:,t,source_index_in(i_source),1))
                    else
                        !Set the time profile
                        emission_time_profile_subgrid(:,:,t,source_index_in(i_source),:)=val_hour_of_week_input(hour_of_week_index,i_source)*val_month_of_year_input(date_array(2),i_source)
                        !if (source_index_in(i_source).eq.industry_index) write(*,*) 'TIME PROFILE: ',val_hour_of_week_input(hour_of_week_index,i_source),val_month_of_year_input(date_array(2),i_source)
                    endif
                    !write(*,*) hour_of_week_index,val_hour_of_week_input(hour_of_week_index,i_source),val_month_of_year_input(date_array(2),i_source)
                    !write(*,*) emission_time_profile_subgrid(1,1,t,source_index_in(i_source),1)
                    !Will  only do this calculation if for heat and only if meteorology exists
                    !This is poorly poisitioned here because it can be called even when no meteorology is available, hence the allocation check
                    do_not_calculate_RWC_emissions=.false.
                    if (source_index_in(i_source).eq.heating_index.and.use_RWC_emission_data) then
                        do j=1,emission_subgrid_dim(y_dim_index,source_index_in(i_source))
                            do i=1,emission_subgrid_dim(x_dim_index,source_index_in(i_source))
                                if (save_emissions_for_EMEP(allsource_index)) then
                                    if (allocated(meteo_var1d_nc)) then
                                        !Need to cross reference the meteo grid to the emission grid as this is not done normally
                                        !Tricky using the two different emep grids?
                                        !Not certain if this 0.5 is correct. Not used else where
                                        i_cross=1+floor((x_emission_subgrid(i,j,source_index_in(i_source))-meteo_var1d_nc(1,lon_nc_index))/meteo_dgrid_nc(lon_nc_index)+0.5)
                                        j_cross=1+floor((y_emission_subgrid(i,j,source_index_in(i_source))-meteo_var1d_nc(1,lat_nc_index))/meteo_dgrid_nc(lat_nc_index)+0.5)
                                        !Because the meteo grid can be smaller than the EMEP grid then need to limit it
                                        !write(*,'(6i12)') i,j,i_cross,j_cross,dim_length_meteo_nc(x_dim_nc_index),dim_length_meteo_nc(y_dim_nc_index)
                                        i_cross=min(max(1,i_cross),dim_length_meteo_nc(x_dim_nc_index))
                                        j_cross=min(max(1,j_cross),dim_length_meteo_nc(y_dim_nc_index))
                                    else
                                        !Do not do this calculation until the meteo data is available
                                        !write(unit_logfile,'(a)') ' No meteo data available for RWC heating emission calculations. Stopping'
                                        !stop
                                        !i_cross=1
                                        !j_cross=1
                                        do_not_calculate_RWC_emissions=.true.
                                    endif
                                else
                                    !Use the EMEP meteorology
                                    i_cross=crossreference_emission_to_emep_subgrid(i,j,x_dim_index,source_index_in(i_source))
                                    j_cross=crossreference_emission_to_emep_subgrid(i,j,y_dim_index,source_index_in(i_source))
                                    i_cross=min(max(1,i_cross),dim_length_nc(x_dim_nc_index))
                                    j_cross=min(max(1,j_cross),dim_length_nc(y_dim_nc_index))
                                endif
                                if (.not.do_not_calculate_RWC_emissions) then
                                    hdd_temp=max(0.,HDD_threshold_value-max(DMT_min_value,DMT_EMEP_grid_nc(i_cross,j_cross,1)))
                                    !write (*,*) hdd_temp,DMT_EMEP_grid_nc(i_cross,j_cross,1),max(DMT_min_value,DMT_EMEP_grid_nc(i_cross,j_cross,1)),HDD_threshold_value-max(DMT_min_value,DMT_EMEP_grid_nc(i_cross,j_cross,1))

                                    !the heating emissions have already been normalised with RWC_grid_HDD in uEMEP_read_RWC_heating_data
                                    emission_time_profile_subgrid(i,j,t,source_index_in(i_source),:)=val_hour_of_week_input(hour_of_week_index,i_source)/24.*hdd_temp
                                    !if (i.eq.1.and.j.eq.1) write(*,*) t,hour_of_week_index,val_hour_of_week_input(hour_of_week_index,i_source)/24.,hdd_temp,HDD_threshold_value,DMT_EMEP_grid_nc(i_cross,j_cross,1)
                                    !if (i.eq.1.and.j.eq.1) write(*,*) t,hour_of_week_index,emission_time_profile_subgrid(i,j,t,source_index_in(i_source),1)
                                endif
                            enddo
                        enddo
                    endif

                endif
            enddo

            if (annual_calculations) then
                emission_time_profile_subgrid(:,:,t,:,:)=1.
            endif

        enddo

        deallocate (time_hour_of_week_input)
        deallocate (val_hour_of_week_input)
        deallocate (time_month_of_year_input)
        deallocate (val_month_of_year_input)


    end subroutine uEMEP_read_time_profiles

end module read_time_profiles

