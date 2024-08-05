module save_netcdf_file

    use uemep_configuration
    use chemistry_no2, only: uEMEP_source_fraction_chemistry
    use mod_read_esri_ascii_file, only: write_esri_ascii_file
    use mod_area_interpolation, only: area_weighted_interpolation_function

    implicit none
    private

    public :: uEMEP_save_netcdf_control, mean_mask, check

contains

!Saves data in netcdf format

    subroutine uEMEP_save_netcdf_control

        use uEMEP_definitions

        implicit none

        integer i,j,k,t,l
        integer i_comp,i_file,i_meteo
        character(256) temp_name,unit_str,title_str,title_str_rec,var_name_temp, temp_date_str,station_name_str,temp_name_rec,temp_compound_str
        logical create_file,create_file_rec
        integer i_source
        real :: valid_min=0.
        real, allocatable :: temp_subgrid(:,:,:)
        real, allocatable :: temp_integral_subgrid(:,:,:)
        real, allocatable :: exhaust_subgrid(:,:,:)
        real, allocatable :: aqi_subgrid(:,:,:,:)
        integer, allocatable :: aqi_responsible_pollutant_index(:,:,:)
        real, allocatable :: temp_subgrid_ascii(:,:)

        integer ii,jj,tt

        real aqi_limits_temp(n_compound_index,1:5)
        real max_aqi
        integer n_aqi_pollutant_index
        parameter (n_aqi_pollutant_index=4)
        integer aqi_pollutant_index(n_aqi_pollutant_index)
        integer i_pollutant,i_loop
        character(256) variable_type
        logical :: receptor_available=.true.
        real scale_factor
        integer n_save_aqi_pollutant_index
        real temp_sum_comp
        integer count
        character(256) filename_ascii

        integer i_sp,ii_sp

        real, allocatable :: subgrid_dummy(:,:,:,:,:,:)
        ! real, allocatable :: comp_subgrid_dummy(:,:,:,:)  ! THIS IS NO LONGER NEEDED, SINCE WE NO LONGER USE REGION LOOP FOR NO2 AND O3 CONTRIBUTION NOR SAVE IN-REGION VERSION OF TOTAL CONCENTRATIONS
        character(256) filename_append
        character(256),parameter :: inregion_suffix = '_from_in_region'   ! same function as filename_append, but does not change with in_region_loop
        integer in_region_loop, n_in_region_loop
        logical save_netcdf_fraction_as_contribution_flag_dummy

        ! THIS IS NO LONGER NEEDED, SINCE WE NO LONGER USE REGION LOOP FOR NO2 AND O3 CONTRIBUTION
        ! real, allocatable :: comp_source_subgrid_dummy(:,:,:,:,:)
        ! real, allocatable :: comp_source_EMEP_subgrid_dummy(:,:,:,:,:)
        ! real, allocatable :: comp_source_EMEP_additional_subgrid_dummy(:,:,:,:,:)

        if (include_o3_in_aqi_index) then
            n_save_aqi_pollutant_index=n_aqi_pollutant_index
        else
            n_save_aqi_pollutant_index=n_aqi_pollutant_index-1
        endif

        if (.not.allocated(temp_subgrid)) allocate(temp_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index)))
        if (.not.allocated(temp_integral_subgrid)) allocate(temp_integral_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),subgrid_dim(t_dim_index)))
        if (.not.allocated(exhaust_subgrid)) allocate(exhaust_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index)))
        if (.not.allocated(aqi_subgrid).and.save_aqi) allocate(aqi_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index))
        if (.not.allocated(aqi_responsible_pollutant_index).and.save_aqi) allocate(aqi_responsible_pollutant_index(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index)))
        if (.not.allocated(temp_subgrid_ascii).and.save_compounds_as_ascii) allocate(temp_subgrid_ascii(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index)))

        !Save subgrid calculations
        valid_min=0.
        unit_str="ug/m3"
        variable_type='byte'
        !variable_type='double'
        variable_type='float'
        scale_factor=1.

        !Save the data
        i_file=subgrid_total_file_index(allsource_index)
        !temp_name=trim(pathname_grid(i_file))//trim(filename_grid(i_file))//'_'//trim(file_tag)//'.nc'
        !if (len(config_date_str).gt.0) then
        !    temp_date_str='_'//trim(config_date_str)
        !else
        !    temp_date_str=''
        !endif

        if (len(filename_date_output_grid).gt.0) then
            temp_date_str='_'//filename_date_output_grid
        else
            temp_date_str=''
        endif

        if (use_multiple_receptor_grids_flag) then
            station_name_str=trim(name_receptor_in(g_loop,1))//'_';
        else
            station_name_str=''
        endif

        !Do not write 'all' in file name if all compounds are selected
        if (pollutant_index.eq.all_nc_index) then
            temp_compound_str=''
        else
            temp_compound_str='_'//trim(var_name_nc(conc_nc_index,compound_index,allsource_index))
        endif
        !Do not write anything about the compounds
        temp_compound_str=''

        !temp_name=trim(pathname_grid(i_file))//trim(station_name_str)//'uEMEP_'//trim(file_tag)//trim(temp_compound_str)//trim(temp_date_str)//'_'//trim(forecast_hour_str)//'.nc'
        !temp_name_rec=trim(pathname_grid(i_file))//'uEMEP_'//trim(file_tag)//'_station'//trim(temp_compound_str)//trim(temp_date_str)//'_'//trim(forecast_hour_str)//'.nc'
        !if (save_netcdf_average_flag) then
        !temp_name=trim(pathname_grid(i_file))//trim(station_name_str)//'uEMEP_'//trim(file_tag)//trim(temp_compound_str)//'_mean'//trim(temp_date_str)//'_'//trim(forecast_hour_str)//'.nc'
        !temp_name_rec=trim(pathname_grid(i_file))//'uEMEP_'//trim(file_tag)//'_station'//trim(temp_compound_str)//'_mean'//trim(temp_date_str)//'_'//trim(forecast_hour_str)//'.nc'
        !endif
        temp_name=trim(pathname_grid(i_file))//trim(station_name_str)//'uEMEP_'//trim(file_tag)//trim(temp_compound_str)//trim(temp_date_str)//'.nc'
        temp_name_rec=trim(pathname_grid(i_file))//'uEMEP_'//trim(file_tag)//'_station'//trim(temp_compound_str)//trim(temp_date_str)//'.nc'
        if (save_netcdf_average_flag) then
            temp_name=trim(pathname_grid(i_file))//trim(station_name_str)//'uEMEP_'//trim(file_tag)//trim(temp_compound_str)//'_mean'//trim(temp_date_str)//'.nc'
            temp_name_rec=trim(pathname_grid(i_file))//'uEMEP_'//trim(file_tag)//'_station'//trim(temp_compound_str)//'_mean'//trim(temp_date_str)//'.nc'
        endif

        finished_file=trim(pathname_grid(i_file))//trim(finished_subpath)//trim(station_name_str)//'uEMEP_'//trim(file_tag)//trim(temp_compound_str)//trim(temp_date_str)//'.'//trim(finished_filename)
        finished_file_rec=trim(pathname_grid(i_file))//trim(finished_subpath)//'uEMEP_'//trim(file_tag)//'_station'//trim(temp_compound_str)//trim(temp_date_str)//'.'//trim(finished_filename)
        if (save_netcdf_average_flag) then
            finished_file=trim(pathname_grid(i_file))//trim(finished_subpath)//trim(station_name_str)//'uEMEP_'//trim(file_tag)//trim(temp_compound_str)//'_mean'//trim(temp_date_str)//'.'//trim(finished_filename)
            finished_file_rec=trim(pathname_grid(i_file))//trim(finished_subpath)//'uEMEP_'//trim(file_tag)//'_station'//trim(temp_compound_str)//'_mean'//trim(temp_date_str)//'.'//trim(finished_filename)
        endif

        write(unit_logfile,'(A)') ''
        write(unit_logfile,'(A)') '================================================================'
        write(unit_logfile,'(A)') 'Saving netcdf data (uEMEP_save_netcdf_control)'
        write(unit_logfile,'(A)') '================================================================'

        i_comp=1
        if (save_netcdf_receptor_flag.and.n_valid_receptor.eq.0) then
            if (i_comp.eq.1.and.t_loop.eq.start_time_loop_index) then
                write(unit_logfile,'(a)')'No receptor positions available. Will not save receptor data.'
                receptor_available=.false.
            endif
        endif

        if (save_netcdf_average_flag) then
            if (.not.allocated(val_array_av)) then
                allocate(val_array_av(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),n_var_av))
                val_array_av=0.
            endif
            if (.not.allocated(time_seconds_output_av)) then
                allocate(time_seconds_output_av(n_var_av))
                time_seconds_output_av=0
            endif
            !Reset average file counter for every time step saving occurs
            counter_av=0
            write(unit_logfile,*) 'Saving as mean data'
        endif

        !Save the final result of the subgrid calculation in ascii format
        !This should only be used for annual mean calculations since it cannot have a time dimmension
        !Intended for FAIRMODE output
        !Makes a new file for each compound and averages over time, if there is a time dimension
        if (save_compounds_as_ascii) then
            variable_type='float'
            unit_str="ug/m3"
            do i_pollutant=1,n_pollutant_loop
                if (pollutant_loop_index(i_pollutant).ne.pmex_index &
                    .and.pollutant_loop_index(i_pollutant).ne.pm25_sand_index.and.pollutant_loop_index(i_pollutant).ne.pm10_sand_index &
                    .and.pollutant_loop_index(i_pollutant).ne.pm25_salt_index.and.pollutant_loop_index(i_pollutant).ne.pm10_salt_index) then
                    do i_loop=1,n_pollutant_compound_loop(i_pollutant)

                        i_comp=pollutant_compound_loop_index(i_pollutant,i_loop)
                        !write(*,*) i_comp

                        var_name_temp=trim(var_name_nc(conc_nc_index,i_comp,allsource_index))

                        title_str=trim(var_name_temp)//'_'//trim(file_tag)//trim(temp_date_str)
                        write(unit_logfile,'(a)')'Writing ascii data to: '//trim(title_str)

                        filename_ascii=trim(pathname_output_grid)//trim(title_str)//'.asc'

                        temp_subgrid_ascii(:,:)=sum(comp_subgrid(:,:,:,i_comp),3)/subgrid_dim(t_dim_index)


                        write(unit_logfile,'(a,f12.3)')'Writing ascii array variable: '//trim(var_name_temp),sum(comp_subgrid(:,:,:,i_comp))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index)

                        call write_esri_ascii_file(filename_ascii,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_delta(x_dim_index),temp_subgrid_ascii,x_subgrid,y_subgrid)

                    enddo
                endif
            enddo
        endif

        !Save the final result of the subgrid calculation
        if (save_compounds) then


            !Loop over the normal subgrid and the in region subgrid by saving  to a dummy variable and putting back when finished
            if (trace_emissions_from_in_region) then
                !original: n_in_region_loop=2
                ! NO LONGER INCLUDE SEPARATE TOTAL CONCENTRATION FOR FROM-IN-REGION. LATER REMOVE THE in_region_loop
                n_in_region_loop=1
                ! if (.not.allocated(comp_subgrid_dummy))allocate (comp_subgrid_dummy(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index))
            else
                n_in_region_loop=1
            endif

            do in_region_loop=1,n_in_region_loop

                write(unit_logfile,'(a)')'--------------------------'
                if (in_region_loop.eq.1) write(unit_logfile,'(a)')'Saving all total compounds'
                if (in_region_loop.eq.2) write(unit_logfile,'(a)')'Saving only regional compounds'
                write(unit_logfile,'(a)')'--------------------------'

                ! THIS IS NO LONGER NEEDED, SINCE WE NO LONGER SAVE IN_REGION VERSION OF TOTAL CONCENTRATIONS
                ! !Loop over the normal subgrid and the in region subgrid by saving  to a dummy variable and pputting back when finished
                ! if (trace_emissions_from_in_region.and.in_region_loop.eq.2) then
                !     filename_append='_from_in_region'
                !     comp_subgrid_dummy=comp_subgrid
                !     comp_subgrid=comp_subgrid_from_in_region
                ! else
                !     filename_append=''
                ! endif


                variable_type='float'
                unit_str="ug/m3"
                do i_pollutant=1,n_pollutant_loop
                    if (pollutant_loop_index(i_pollutant).ne.pmex_index &
                        .and.pollutant_loop_index(i_pollutant).ne.pm25_sand_index.and.pollutant_loop_index(i_pollutant).ne.pm10_sand_index &
                        .and.pollutant_loop_index(i_pollutant).ne.pm25_salt_index.and.pollutant_loop_index(i_pollutant).ne.pm10_salt_index) then
                        do i_loop=1,n_pollutant_compound_loop(i_pollutant)

                            i_comp=pollutant_compound_loop_index(i_pollutant,i_loop)
                            !write(*,*) i_comp

                            if (i_pollutant.eq.1.and.i_loop.eq.1.and.t_loop.eq.start_time_loop_index.and.save_netcdf_file_flag.and.in_region_loop.eq.1) then
                                create_file=.true.
                                title_str='uEMEP_concentration_'//trim(file_tag)//temp_date_str
                                write(unit_logfile,'(a)')'Writing to: '//trim(temp_name)
                            else
                                create_file=.false.
                            endif

                            if (i_pollutant.eq.1.and.i_loop.eq.1.and.t_loop.eq.start_time_loop_index.and.first_g_loop.and.save_netcdf_receptor_flag.and.in_region_loop.eq.1) then
                                create_file_rec=.true.
                                title_str_rec='uEMEP_receptor_'//trim(file_tag)//temp_date_str
                                if (receptor_available) write(unit_logfile,'(a)')'Writing to: '//trim(temp_name_rec)
                            else
                                create_file_rec=.false.
                            endif

                            var_name_temp=trim(var_name_nc(conc_nc_index,i_comp,allsource_index))//'_concentration'//filename_append
                            if (save_netcdf_file_flag) then
                                temp_sum_comp=0.
                                count=0
                                do j=1,subgrid_dim(y_dim_index)
                                    do i=1,subgrid_dim(x_dim_index)
                                        if (use_subgrid(i,j,allsource_index)) then
                                            temp_sum_comp=temp_sum_comp+sum(comp_subgrid(i,j,:,i_comp))/subgrid_dim(t_dim_index)
                                            count=count+1
                                        endif
                                    enddo
                                enddo
                                if (count.gt.0) then
                                    temp_sum_comp=temp_sum_comp/count
                                else
                                    temp_sum_comp=0
                                endif
                                !write(unit_logfile,'(a,f12.3)')'Writing netcdf array variable:    '//trim(var_name_temp),temp_sum_comp
                                write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_mask(comp_subgrid(:,:,:,i_comp),use_subgrid(:,:,allsource_index),size(comp_subgrid(:,:,:,i_comp),1),size(comp_subgrid(:,:,:,i_comp),2),size(comp_subgrid(:,:,:,i_comp),3))
                                call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                    ,comp_subgrid(:,:,:,i_comp),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                    ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                            endif
                            if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                                write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_mask(comp_subgrid(:,:,:,i_comp),use_subgrid(:,:,allsource_index),size(comp_subgrid(:,:,:,i_comp),1),size(comp_subgrid(:,:,:,i_comp),2),size(comp_subgrid(:,:,:,i_comp),3))
                                !write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(comp_subgrid(:,:,:,i_comp))/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index)
                                call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                    ,comp_subgrid(:,:,:,i_comp),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                    ,unit_str,title_str_rec,create_file_rec,valid_min &
                                    ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                    ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                    ,z_rec(allsource_index,1) &
                                    ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                            endif
                        enddo
                    endif
                enddo

                ! THIS IS NO LONGER NEEDED, SINCE WE NO LONGER SAVE IN_REGION VERSION OF TOTAL CONCENTRATIONS
                ! if (trace_emissions_from_in_region.and.in_region_loop.eq.2) then
                !     comp_subgrid_from_in_region=comp_subgrid
                !     comp_subgrid=comp_subgrid_dummy
                ! endif

            enddo !in region loop
        endif


        ! THIS IS NO LONGER NEEDED, SINCE WE NO LONGER SAVE IN_REGION VERSION OF TOTAL CONCENTRATIONS
        ! if (trace_emissions_from_in_region) then
        !     if (allocated(comp_subgrid_dummy)) deallocate (comp_subgrid_dummy)
        ! endif

        create_file=.false.
        create_file_rec=.false.

        !Loop over the normal subgrid and the in region subgrid by saving  to a dummy variable and pputting back when finished
        if (trace_emissions_from_in_region) then
            n_in_region_loop=2
            if (.not.allocated(subgrid_dummy)) allocate (subgrid_dummy(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_subgrid_index,n_source_index,n_pollutant_loop))
            subgrid_dummy=0
        else
            n_in_region_loop=1
        endif

        ! THIS IS NO LONGER NEEDED, SINCE WE NO LONGER USE REGION LOOP FOR NO2 AND O3 CONTRIBUTION
        ! if (trace_emissions_from_in_region.and.(save_no2_source_contributions.or.save_o3_source_contributions)) then
        !     if (.not.allocated(comp_subgrid_dummy))allocate (comp_subgrid_dummy(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index))
        !     subgrid_dummy=0
        !     comp_subgrid_dummy=0
        !     if (.not.allocated(comp_source_subgrid_dummy)) allocate(comp_source_subgrid_dummy(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index,n_source_index))
        !     if (.not.allocated(comp_source_EMEP_subgrid_dummy)) allocate(comp_source_EMEP_subgrid_dummy(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index,n_source_index))
        !     if (.not.allocated(comp_source_EMEP_additional_subgrid_dummy)) allocate(comp_source_EMEP_additional_subgrid_dummy(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_compound_index,n_source_index))
        !     comp_source_subgrid_dummy=0
        !     comp_source_EMEP_subgrid_dummy=0
        !     comp_source_EMEP_additional_subgrid_dummy=0

        ! else
        !     n_in_region_loop=1
        ! endif

        do in_region_loop=1,n_in_region_loop

            write(unit_logfile,'(a)')'--------------------------'
            if (in_region_loop.eq.1) write(unit_logfile,'(a)')'Saving all contributions'
            if (in_region_loop.eq.2) write(unit_logfile,'(a)')'Saving only regional contributions'
            write(unit_logfile,'(a)')'--------------------------'

            !Loop over the normal subgrid and the in region subgrid by saving  to a dummy variable and pputting back when finished
            if (trace_emissions_from_in_region.and.in_region_loop.eq.2) then
                filename_append='_from_in_region'
                subgrid_dummy=subgrid
                subgrid=subgrid_from_in_region
                save_netcdf_fraction_as_contribution_flag_dummy=save_netcdf_fraction_as_contribution_flag
                save_netcdf_fraction_as_contribution_flag=save_netcdf_fraction_as_contribution_from_in_region_flag
            else
                filename_append=''
            endif


            !Save the different local source contributions, not the total though
            if (save_source_contributions) then
                if (save_netcdf_fraction_as_contribution_flag) then
                    variable_type='float'
                    unit_str="ug/m3"
                else
                    variable_type='byte'
                    unit_str="%"
                endif

                do i_pollutant=1,n_pollutant_loop
                    do i_source=1,n_source_index
                        !if (calculate_source(i_source).or.i_source.eq.allsource_index) then
                        !Don't save any exhaust, sand or salt pollutant sources. Dealt with later
                        !(calculate_source(i_source).or.calculate_emep_source(i_source))
                        if (calculate_source(i_source).and.i_source.ne.allsource_index.and.pollutant_loop_index(i_pollutant).ne.pmex_index &
                            .and.pollutant_loop_index(i_pollutant).ne.pm25_sand_index.and.pollutant_loop_index(i_pollutant).ne.pm10_sand_index &
                            .and.pollutant_loop_index(i_pollutant).ne.pm25_salt_index.and.pollutant_loop_index(i_pollutant).ne.pm10_salt_index) then

                            i_file=subgrid_local_file_index(i_source)

                            !Only save nonexhaust pm and all the other sources
                            if ((pollutant_loop_index(i_pollutant).eq.nox_index.and.i_source.eq.traffic_index)) then
                                !if ((pollutant_loop_index(i_pollutant).ne.nox_index.or.i_source.ne.traffic_index)) then
                            else
                                if (i_source.eq.traffic_index.and.(pollutant_loop_index(i_pollutant).eq.pm10_index.or.pollutant_loop_index(i_pollutant).eq.pm25_index)) then
                                    if (pollutant_index.eq.all_totals_nc_index.or.pollutant_index.eq.aaqd_totals_nc_index.or.pollutant_index.eq.gp_totals_nc_index.or.pollutant_index.eq.op_totals_nc_index) then
                                        var_name_temp=trim(var_name_nc(conc_nc_index,pollutant_loop_index(i_pollutant),allsource_index))//'_'//trim(filename_grid(i_file))//trim(filename_append)
                                        if (save_netcdf_fraction_as_contribution_flag) then
                                            temp_subgrid=subgrid(:,:,:,local_subgrid_index,i_source,i_pollutant)
                                        else
                                            temp_subgrid=subgrid(:,:,:,local_subgrid_index,i_source,i_pollutant)/subgrid(:,:,:,total_subgrid_index,allsource_index,i_pollutant)*100.
                                        endif
                                    else
                                        var_name_temp=trim(var_name_nc(conc_nc_index,pollutant_loop_index(i_pollutant),allsource_index))//'_'//trim(filename_grid(i_file))//'_nonexhaust'//trim(filename_append)
                                        if (save_netcdf_fraction_as_contribution_flag) then
                                            temp_subgrid=(subgrid(:,:,:,local_subgrid_index,i_source,i_pollutant)-subgrid(:,:,:,local_subgrid_index,i_source,pollutant_loop_back_index(pmex_index)))
                                        else
                                            temp_subgrid=(subgrid(:,:,:,local_subgrid_index,i_source,i_pollutant)-subgrid(:,:,:,local_subgrid_index,i_source,pollutant_loop_back_index(pmex_index)))/subgrid(:,:,:,total_subgrid_index,allsource_index,i_pollutant)*100.
                                        endif
                                    endif
                                else
                                    var_name_temp=trim(var_name_nc(conc_nc_index,pollutant_loop_index(i_pollutant),allsource_index))//'_'//trim(filename_grid(i_file))//trim(filename_append)
                                    if (save_netcdf_fraction_as_contribution_flag) then
                                        temp_subgrid=subgrid(:,:,:,local_subgrid_index,i_source,i_pollutant)
                                    else
                                        temp_subgrid=subgrid(:,:,:,local_subgrid_index,i_source,i_pollutant)/subgrid(:,:,:,total_subgrid_index,allsource_index,i_pollutant)*100.
                                    endif
                                endif
                                !In case of any 0 concentrations when using the fractional contributions
                                where (subgrid(:,:,:,total_subgrid_index,allsource_index,i_pollutant).eq.0) temp_subgrid=0

                                if (save_netcdf_file_flag) then
                                    !write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                                    !write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_nodata(temp_subgrid,size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3),NODATA_value)
                                    write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_mask(temp_subgrid,use_subgrid(:,:,allsource_index),size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3))


                                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                        ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                                endif
                                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                                    write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_mask(temp_subgrid,use_subgrid(:,:,allsource_index),size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3))
                                    ! write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                        ,unit_str,title_str_rec,create_file_rec,valid_min &
                                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                        ,z_rec(allsource_index,1) &
                                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                                endif

                            endif

                            !Special case for exhaust as this must be given as a fraction for both PM2.5 and PM10.
                            !if ((pollutant_loop_index(i_pollutant).eq.pm10_index.or.pollutant_loop_index(i_pollutant).eq.pm25_index).and.i_source.eq.traffic_index) then
                            !Write all exhaust pollutants except the pmex
                            if (i_source.eq.traffic_index.and.(.not.pollutant_index.eq.all_totals_nc_index.and..not.pollutant_index.eq.aaqd_totals_nc_index.and..not.pollutant_index.eq.gp_totals_nc_index.and..not.pollutant_index.eq.op_totals_nc_index.or.pollutant_loop_index(i_pollutant).eq.nox_index)) then

                                !If PM then exhaust output is exhaust
                                if (pollutant_loop_index(i_pollutant).eq.pm10_index.or.pollutant_loop_index(i_pollutant).eq.pm25_index) then
                                    exhaust_subgrid=subgrid(:,:,:,local_subgrid_index,i_source,pollutant_loop_back_index(pmex_index))
                                else
                                    exhaust_subgrid=subgrid(:,:,:,local_subgrid_index,i_source,i_pollutant)
                                endif

                                if (save_netcdf_fraction_as_contribution_flag) then
                                    var_name_temp=trim(var_name_nc(conc_nc_index,pollutant_loop_index(i_pollutant),allsource_index))//'_'//'local_contribution_traffic_exhaust'//trim(filename_append)
                                    temp_subgrid=exhaust_subgrid
                                else
                                    var_name_temp=trim(var_name_nc(conc_nc_index,pollutant_loop_index(i_pollutant),allsource_index))//'_'//'local_fraction_traffic_exhaust'//trim(filename_append)
                                    temp_subgrid=exhaust_subgrid/subgrid(:,:,:,total_subgrid_index,allsource_index,i_pollutant)*100.
                                endif
                                where (subgrid(:,:,:,total_subgrid_index,allsource_index,i_pollutant).eq.0) temp_subgrid=0

                                if (save_netcdf_file_flag) then
                                    !write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                                    write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_mask(temp_subgrid,use_subgrid(:,:,allsource_index),size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3))
                                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                        ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                                endif
                                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                                    write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_mask(temp_subgrid,use_subgrid(:,:,allsource_index),size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3))
                                    !write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                        ,unit_str,title_str_rec,create_file_rec,valid_min &
                                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                        ,z_rec(allsource_index,1) &
                                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                                endif

                            endif
                        endif

                        !Special case for salt and sand
                        if (pollutant_loop_index(i_pollutant).eq.pm25_sand_index.or.pollutant_loop_index(i_pollutant).eq.pm10_sand_index &
                            .or.pollutant_loop_index(i_pollutant).eq.pm25_salt_index.or.pollutant_loop_index(i_pollutant).eq.pm10_salt_index) then
                            !Save the nonexhaust sand and salt
                            if (i_source.eq.traffic_index) then

                                i_file=subgrid_local_file_index(i_source)
                                var_name_temp=trim(var_name_nc(conc_nc_index,pollutant_loop_index(i_pollutant),allsource_index))//'_'//trim(filename_grid(i_file))//trim(filename_append)
                                if (save_netcdf_fraction_as_contribution_flag) then
                                    temp_subgrid=subgrid(:,:,:,local_subgrid_index,i_source,i_pollutant)
                                else
                                    if (pollutant_loop_index(i_pollutant).eq.pm10_sand_index.or.pollutant_loop_index(i_pollutant).eq.pm10_salt_index) then
                                        temp_subgrid=subgrid(:,:,:,local_subgrid_index,i_source,i_pollutant)/subgrid(:,:,:,total_subgrid_index,allsource_index,pollutant_loop_back_index(pm10_index))*100.
                                    else
                                        temp_subgrid=subgrid(:,:,:,local_subgrid_index,i_source,i_pollutant)/subgrid(:,:,:,total_subgrid_index,allsource_index,pollutant_loop_back_index(pm25_index))*100.
                                    endif

                                endif
                                where (subgrid(:,:,:,total_subgrid_index,allsource_index,i_pollutant).eq.0) temp_subgrid=0

                                if (save_netcdf_file_flag) then
                                    write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_mask(temp_subgrid,use_subgrid(:,:,allsource_index),size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3))
                                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                        ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                                endif
                                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                                    write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                        ,unit_str,title_str_rec,create_file_rec,valid_min &
                                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                        ,z_rec(allsource_index,1) &
                                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                                endif

                            endif
                        endif

                        if (pollutant_loop_index(i_pollutant).ne.pmex_index &
                            .and.pollutant_loop_index(i_pollutant).ne.pm25_sand_index.and.pollutant_loop_index(i_pollutant).ne.pm10_sand_index &
                            .and.pollutant_loop_index(i_pollutant).ne.pm25_salt_index.and.pollutant_loop_index(i_pollutant).ne.pm10_salt_index) then
                            !Save the total nonlocal part from EMEP
                            if (i_source.eq.allsource_index) then

                                i_file=emep_subgrid_nonlocal_file_index(i_source)
                                var_name_temp=trim(var_name_nc(conc_nc_index,pollutant_loop_index(i_pollutant),allsource_index))//'_'//trim(filename_grid(i_file))//trim(filename_append)
                                if (save_netcdf_fraction_as_contribution_flag) then
                                    temp_subgrid=subgrid(:,:,:,emep_nonlocal_subgrid_index,i_source,i_pollutant)
                                else
                                    temp_subgrid=subgrid(:,:,:,emep_nonlocal_subgrid_index,i_source,i_pollutant)/subgrid(:,:,:,total_subgrid_index,allsource_index,i_pollutant)*100.
                                endif
                                where (subgrid(:,:,:,total_subgrid_index,allsource_index,i_pollutant).eq.0) temp_subgrid=0

                                if (save_netcdf_file_flag) then
                                    write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_mask(temp_subgrid,use_subgrid(:,:,allsource_index),size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3))
                                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                        ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                                endif
                                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                                    write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                        ,unit_str,title_str_rec,create_file_rec,valid_min &
                                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                        ,z_rec(allsource_index,1) &
                                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                                endif

                                if (EMEP_additional_grid_interpolation_size.gt.0) then
                                    i_file=emep_additional_subgrid_nonlocal_file_index(i_source)
                                    var_name_temp=trim(var_name_nc(conc_nc_index,pollutant_loop_index(i_pollutant),allsource_index))//'_'//trim(filename_grid(i_file))//trim(filename_append)
                                    if (save_netcdf_fraction_as_contribution_flag) then
                                        temp_subgrid=subgrid(:,:,:,emep_additional_nonlocal_subgrid_index,i_source,i_pollutant)
                                    else
                                        temp_subgrid=subgrid(:,:,:,emep_additional_nonlocal_subgrid_index,i_source,i_pollutant)/subgrid(:,:,:,total_subgrid_index,allsource_index,i_pollutant)*100.
                                    endif
                                    where (subgrid(:,:,:,total_subgrid_index,allsource_index,i_pollutant).eq.0) temp_subgrid=0

                                    if (save_netcdf_file_flag) then
                                        write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_mask(temp_subgrid,use_subgrid(:,:,allsource_index),size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3))
                                        call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                            ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                            ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                                    endif
                                    if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                                        write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                                        call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                            ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                            ,unit_str,title_str_rec,create_file_rec,valid_min &
                                            ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                            ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                            ,z_rec(allsource_index,1) &
                                            ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                                    endif
                                endif

                            endif
                        endif
                    enddo
                enddo
            endif

            ! Save NO2 and O3 source contributions
            ! NB: This is no longer repeated in the second iteration of the in-region loop because
            ! the subroutine "uEMEP_source_fraction_chemistry" now calculates both normal and in_region contributions
            ! to NO2 and O3 in a single call. Note that in_region contributions are no longer defined for nonlocal
            ! nor for the additional domain
            if (in_region_loop == 1 .and. (save_no2_source_contributions.or.save_o3_source_contributions)) then

                write(unit_logfile,'(A)') '-------------------------------'
                write(unit_logfile,'(A)') 'Saving NO2 and O3 contributions'
                write(unit_logfile,'(A)') '-------------------------------'

                ! THIS IS NO LONGER NEEDED, SINCE WE NO LONGER USE REGION LOOP FOR NO2 AND O3 CONTRIBUTION
                ! if (trace_emissions_from_in_region.and.in_region_loop.eq.2) then
                !     !Save the normal variables in a dummy array
                !     comp_source_subgrid_dummy=comp_source_subgrid
                !     comp_source_EMEP_subgrid_dummy=comp_source_EMEP_subgrid
                !     comp_source_EMEP_additional_subgrid_dummy=comp_source_EMEP_additional_subgrid
                !     !Set in the in region variables
                !     comp_source_subgrid=comp_source_subgrid_from_in_region
                !     comp_source_EMEP_subgrid=comp_source_EMEP_subgrid_from_in_region
                !     comp_source_EMEP_additional_subgrid=comp_source_EMEP_additional_subgrid_from_in_region
                ! endif

                if (EMEP_additional_grid_interpolation_size.gt.0) then
                    calculate_EMEP_additional_grid_flag=.true.
                    call uEMEP_source_fraction_chemistry
                endif

                calculate_EMEP_additional_grid_flag=.false.
                call uEMEP_source_fraction_chemistry


                if (save_no2_source_contributions) then

                    valid_min=0.

                    if (save_netcdf_fraction_as_contribution_flag) then
                        variable_type='float'
                        unit_str="ug/m3"
                    else
                        variable_type='byte'
                        unit_str="%"
                    endif

                    do i_source=1,n_source_index
                        !if (calculate_source(i_source).or.i_source.eq.allsource_index) then
                        !if (calculate_source(i_source).or.i_source.eq.allsource_index) then
                        if (calculate_source(i_source).or.i_source.eq.allsource_index.or.(calculate_emep_source(i_source).and..not.calculate_source(i_source))) then
                            !Save all EMEP NO2 contributions
                            !if (calculate_source(i_source).or.i_source.eq.allsource_index.or.calculate_emep_source(i_source)) then
                            if (i_source.eq.traffic_nonexhaust_index) then
                                !Do not save nonexhaust for exhaust gas emissions
                            else

                                if (i_source.eq.allsource_index) then
                                    !if (i_source.eq.allsource_index.or.(calculate_emep_source(i_source).and..not.calculate_source(i_source))) then
                                    i_file=emep_subgrid_nonlocal_file_index(i_source)
                                else
                                    i_file=subgrid_local_file_index(i_source)
                                endif
                                if (i_source.ne.allsource_index.and.(calculate_emep_source(i_source).and..not.calculate_source(i_source)).and.save_emep_source_contributions) then
                                    i_file=emep_subgrid_local_file_index(i_source)
                                endif

                                if (i_source.eq.traffic_index) then
                                    var_name_temp=trim(var_name_nc(conc_nc_index,no2_nc_index,allsource_nc_index))//'_'//trim(filename_grid(i_file))//'_exhaust'//trim(filename_append)
                                else
                                    var_name_temp=trim(var_name_nc(conc_nc_index,no2_nc_index,allsource_nc_index))//'_'//trim(filename_grid(i_file))//trim(filename_append)
                                endif
                                if (save_netcdf_fraction_as_contribution_flag) then
                                    temp_subgrid=comp_source_subgrid(:,:,:,no2_index,i_source)
                                else
                                    temp_subgrid=comp_source_subgrid(:,:,:,no2_index,i_source)/comp_subgrid(:,:,:,no2_index)*100.
                                endif
                                if (i_source.eq.allsource_index) then
                                    if (save_netcdf_fraction_as_contribution_flag) then
                                        temp_subgrid=comp_source_EMEP_subgrid(:,:,:,no2_index,i_source)
                                    else
                                        !temp_subgrid=comp_source_EMEP_subgrid(:,:,:,no2_index,i_source)/comp_EMEP_subgrid(:,:,:,no2_index)*100.
                                        temp_subgrid=comp_source_EMEP_subgrid(:,:,:,no2_index,i_source)/comp_subgrid(:,:,:,no2_index)*100.
                                    endif
                                endif


                                if (save_netcdf_file_flag) then
                                    write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_mask(temp_subgrid,use_subgrid(:,:,allsource_index),size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3))
                                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                        ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                                endif
                                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                                    write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                        ,unit_str,title_str_rec,create_file_rec,valid_min &
                                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                        ,z_rec(allsource_index,1) &
                                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                                endif
                            endif
                        endif

                        !Save the additional EMEP nonlocal source contributions here
                        if (i_source.eq.allsource_index) then
                            if (EMEP_additional_grid_interpolation_size.gt.0) then
                                i_file=emep_additional_subgrid_nonlocal_file_index(i_source)
                                var_name_temp=trim(var_name_nc(conc_nc_index,no2_nc_index,allsource_nc_index))//'_'//trim(filename_grid(i_file))//trim(filename_append)
                                if (save_netcdf_fraction_as_contribution_flag) then
                                    temp_subgrid=comp_source_EMEP_additional_subgrid(:,:,:,no2_index,i_source)
                                else
                                    temp_subgrid=comp_source_EMEP_additional_subgrid(:,:,:,no2_index,i_source)/comp_EMEP_subgrid(:,:,:,no2_index)*100.
                                endif

                                if (save_netcdf_file_flag) then
                                    write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_mask(temp_subgrid,use_subgrid(:,:,allsource_index),size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3))
                                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                        ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                                endif
                                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                                    write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                        ,unit_str,title_str_rec,create_file_rec,valid_min &
                                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                        ,z_rec(allsource_index,1) &
                                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                                endif

                            endif
                        endif
                    enddo

                    ! Save in-region version of NO2 concentration for all sources, except for non-local (i.e. allsource_index is NOT included)
                    ! ********************************************************************************************************
                    if (trace_emissions_from_in_region) then
                        if (save_netcdf_fraction_as_contribution_from_in_region_flag) then
                            variable_type='float'
                            unit_str="ug/m3"
                        else
                            variable_type='byte'
                            unit_str="%"
                        endif

                        do i_source=1,n_source_index
                            ! Skip sources not to be calculated in uEMEP or EMEP (NB: also explicitly skip allsource)
                            if (i_source == allsource_index .or. .not. (calculate_source(i_source) .or. calculate_emep_source(i_source))) then
                                cycle
                            end if
                            !Save all EMEP NO2 contributions
                            !if (calculate_source(i_source).or.i_source.eq.allsource_index.or.calculate_emep_source(i_source)) then
                            if (i_source.eq.traffic_nonexhaust_index) then
                                !Do not save nonexhaust for exhaust gas emissions
                                cycle
                            end if
                            i_file=subgrid_local_file_index(i_source)
                            if ((calculate_emep_source(i_source) .and. .not. calculate_source(i_source))) then
                                ! NB: skipping check that save_emep_source_contributions==.true.
                                i_file=emep_subgrid_local_file_index(i_source)
                            end if

                            if (i_source.eq.traffic_index) then
                                var_name_temp=trim(var_name_nc(conc_nc_index,no2_nc_index,allsource_nc_index))//'_'//trim(filename_grid(i_file))//'_exhaust'//trim(inregion_suffix)
                            else
                                var_name_temp=trim(var_name_nc(conc_nc_index,no2_nc_index,allsource_nc_index))//'_'//trim(filename_grid(i_file))//trim(inregion_suffix)
                            end if
                            if (save_netcdf_fraction_as_contribution_flag) then
                                temp_subgrid=comp_source_subgrid_from_in_region(:,:,:,no2_index,i_source)
                            else
                                ! NB: normalizing with the total concentration, taken from 'comp_subgrid', not the in-region array
                                temp_subgrid=comp_source_subgrid_from_in_region(:,:,:,no2_index,i_source)/comp_subgrid(:,:,:,no2_index)*100.
                            end if

                            if (save_netcdf_file_flag) then
                                write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_mask(temp_subgrid,use_subgrid(:,:,allsource_index),size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3))
                                call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                    ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                    ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                            end if
                            if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                                write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                                call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                    ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                    ,unit_str,title_str_rec,create_file_rec,valid_min &
                                    ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                    ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                    ,z_rec(allsource_index,1) &
                                    ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                            end if
                        end do
                    end if
                    ! **************************************************
                    ! Done saving in-region version of NO2 contributions

                    valid_min=0.

                endif  !save_no2_source_contributions

                if (save_o3_source_contributions) then

                    if (save_netcdf_fraction_as_contribution_flag) then
                        variable_type='float'
                        unit_str="ug/m3"
                        valid_min=-1000.
                    else
                        variable_type='short'
                        unit_str="%"
                        valid_min=-100.
                    endif


                    do i_source=1,n_source_index
                        !if (calculate_source(i_source).or.i_source.eq.allsource_index) then
                        if (calculate_source(i_source).or.i_source.eq.allsource_index.or.(calculate_emep_source(i_source).and..not.calculate_source(i_source))) then
                            !if (calculate_source(i_source).or.i_source.eq.allsource_index.or.calculate_emep_source(i_source)) then
                            !if (calculate_source(i_source).or.i_source.eq.allsource_index.or.(calculate_emep_source(i_source).and..not.calculate_source(i_source))) then
                            !Only save the nonlocal part as 100% save_netcdf_fraction_as_contribution_flag=.false.
                            if ((i_source.eq.allsource_index.and..not.save_netcdf_fraction_as_contribution_flag).or.save_netcdf_fraction_as_contribution_flag) then
                                if (i_source.eq.traffic_nonexhaust_nc_index) then
                                    !Do not save nonexhaust for exhaust gas emissions
                                else

                                    if (i_source.eq.allsource_index) then
                                        i_file=emep_subgrid_nonlocal_file_index(i_source)
                                    else
                                        i_file=subgrid_local_file_index(i_source)
                                    endif
                                    if (i_source.ne.allsource_index.and.(calculate_emep_source(i_source).and..not.calculate_source(i_source)).and.save_emep_source_contributions) then
                                        i_file=emep_subgrid_local_file_index(i_source)
                                    endif

                                    var_name_temp=trim(var_name_nc(conc_nc_index,o3_nc_index,allsource_nc_index))//'_'//trim(filename_grid(i_file))//trim(filename_append)
                                    if (save_netcdf_fraction_as_contribution_flag) then
                                        !temp_subgrid=comp_source_EMEP_subgrid(:,:,:,o3_index,i_source)
                                        temp_subgrid=comp_source_subgrid(:,:,:,o3_index,i_source)
                                        !temp_subgrid=subgrid(:,:,:,emep_nonlocal_subgrid_index,i_source,pollutant_loop_back_index(o3_nc_index))
                                    else
                                        temp_subgrid=comp_source_subgrid(:,:,:,o3_index,i_source)/comp_subgrid(:,:,:,o3_index)*100.
                                        !Put some limits on the ozone because it can be very large if o3 is small
                                        temp_subgrid=min(temp_subgrid,1000.)
                                        temp_subgrid=max(temp_subgrid,-100.)

                                        !temp_subgrid=comp_source_EMEP_subgrid(:,:,:,o3_index,i_source)/comp_EMEP_subgrid(:,:,:,o3_index)*100.
                                        !temp_subgrid=100.*subgrid(:,:,:,emep_nonlocal_subgrid_index,i_source,pollutant_loop_back_index(o3_nc_index))/subgrid(:,:,:,emep_subgrid_index,i_source,pollutant_loop_back_index(o3_nc_index))
                                    endif
                                    if (i_source.eq.allsource_index) then
                                        if (save_netcdf_fraction_as_contribution_flag) then
                                            temp_subgrid=comp_source_EMEP_subgrid(:,:,:,o3_index,i_source)
                                        else
                                            !Note this divides by the EMEP value not the subgrid value because it has to be 100%. Do not use fractions for O3, asking for trouble
                                            temp_subgrid=comp_EMEP_subgrid(:,:,:,o3_index)/comp_EMEP_subgrid(:,:,:,o3_index)*100.
                                            !temp_subgrid=comp_source_EMEP_subgrid(:,:,:,o3_index,i_source)/comp_EMEP_subgrid(:,:,:,o3_index)*100.
                                            !temp_subgrid=comp_source_EMEP_subgrid(:,:,:,o3_index,i_source)/comp_subgrid(:,:,:,o3_index)*100.
                                            temp_subgrid=min(temp_subgrid,1000.)
                                            temp_subgrid=max(temp_subgrid,-100.)
                                        endif
                                    endif

                                    if (save_netcdf_file_flag) then
                                        write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_mask(temp_subgrid,use_subgrid(:,:,allsource_index),size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3))
                                        call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                            ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                            ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                                    endif
                                    if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                                        write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                                        call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                            ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                            ,unit_str,title_str_rec,create_file_rec,valid_min &
                                            ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                            ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                            ,z_rec(allsource_index,1) &
                                            ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                                    endif

                                    !Save the additional EMEP source contributions here, only background
                                    if (i_source.eq.allsource_index) then
                                        if (EMEP_additional_grid_interpolation_size.gt.0) then
                                            i_file=emep_additional_subgrid_nonlocal_file_index(i_source)
                                            var_name_temp=trim(var_name_nc(conc_nc_index,o3_nc_index,allsource_nc_index))//'_'//trim(filename_grid(i_file))//trim(filename_append)
                                            if (save_netcdf_fraction_as_contribution_flag) then
                                                temp_subgrid=comp_source_EMEP_additional_subgrid(:,:,:,o3_index,i_source)
                                                !temp_subgrid=subgrid(:,:,:,emep_nonlocal_subgrid_index,i_source,pollutant_loop_back_index(o3_nc_index))
                                            else
                                                temp_subgrid=comp_source_EMEP_additional_subgrid(:,:,:,o3_index,i_source)/comp_EMEP_subgrid(:,:,:,o3_index)*100.
                                                !temp_subgrid=100.*subgrid(:,:,:,emep_nonlocal_subgrid_index,i_source,pollutant_loop_back_index(o3_nc_index))/subgrid(:,:,:,emep_subgrid_index,i_source,pollutant_loop_back_index(o3_nc_index))
                                                temp_subgrid=min(temp_subgrid,1000.)
                                                temp_subgrid=max(temp_subgrid,-100.)
                                            endif

                                            if (save_netcdf_file_flag) then
                                                write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_mask(temp_subgrid,use_subgrid(:,:,allsource_index),size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3))
                                                call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                                    ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                                    ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                                            endif
                                            if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                                                write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                                                call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                                    ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                                    ,unit_str,title_str_rec,create_file_rec,valid_min &
                                                    ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                                    ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                                    ,z_rec(allsource_index,1) &
                                                    ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                                            endif

                                        endif
                                    endif

                                endif
                            endif
                        endif
                    enddo

                    ! Save in-region version of O3 concentration for all sources, except for non-local (i.e. allsource_index is NOT included)
                    ! ********************************************************************************************************
                    if (trace_emissions_from_in_region .and. save_netcdf_fraction_as_contribution_from_in_region_flag) then
                        ! NB: Not implementing the case of fraction-as-contribution, as is implemented but not used for the normal O3 source contributions above
                        variable_type='float'
                        unit_str="ug/m3"
                        valid_min=-1000.
                        do i_source=1,n_source_index
                            ! Skip sources not to be calculated in uEMEP or EMEP (NB: also explicitly skip allsource)
                            if (i_source == allsource_index .or. .not. (calculate_source(i_source) .or. calculate_emep_source(i_source))) then
                                cycle
                            end if
                            if (i_source.eq.traffic_nonexhaust_nc_index) then
                                !Do not save nonexhaust for exhaust gas emissions
                                cycle
                            end if
                            if (i_source.eq.allsource_index) then
                                i_file=emep_subgrid_nonlocal_file_index(i_source)
                            else
                                i_file=subgrid_local_file_index(i_source)
                            endif
                            if (calculate_emep_source(i_source) .and. .not. calculate_source(i_source)) then
                                ! NB: skipping check that save_emep_source_contributions==.true.
                                i_file=emep_subgrid_local_file_index(i_source)
                            endif

                            var_name_temp=trim(var_name_nc(conc_nc_index,o3_nc_index,allsource_nc_index))//'_'//trim(filename_grid(i_file))//trim(inregion_suffix)
                            temp_subgrid=comp_source_subgrid_from_in_region(:,:,:,o3_index,i_source)

                            if (save_netcdf_file_flag) then
                                write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_mask(temp_subgrid,use_subgrid(:,:,allsource_index),size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3))
                                call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                    ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                    ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                            endif
                            if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                                write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                                call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                    ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                    ,unit_str,title_str_rec,create_file_rec,valid_min &
                                    ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                    ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                    ,z_rec(allsource_index,1) &
                                    ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                            endif
                        enddo
                    end if
                    ! **************************************************
                    ! Done saving in-region version of O3 contributions

                    valid_min=0.

                endif  !save_o3_source_contributions

                ! THIS IS NO LONGER NEEDED, SINCE WE NO LONGER USE REGION LOOP FOR NO2 AND O3 CONTRIBUTION
                ! if (trace_emissions_from_in_region.and.in_region_loop.eq.2) then
                !     !Save the in region variables
                !     comp_source_subgrid_from_in_region=comp_source_subgrid
                !     comp_source_EMEP_subgrid_from_in_region=comp_source_EMEP_subgrid
                !     comp_source_EMEP_additional_subgrid_from_in_region=comp_source_EMEP_additional_subgrid
                !     !Restore the normal variables from the dummy array
                !     comp_source_subgrid=comp_source_subgrid_dummy
                !     comp_source_EMEP_subgrid=comp_source_EMEP_subgrid_dummy
                !     comp_source_EMEP_additional_subgrid=comp_source_EMEP_additional_subgrid_dummy
                ! endif

            endif  !(in_region_loop == 1 .and. (save_no2_source_contributions.or.save_o3_source_contributions))

            if (trace_emissions_from_in_region.and.in_region_loop.eq.2) then
                subgrid_from_in_region=subgrid
                subgrid=subgrid_dummy
                save_netcdf_fraction_as_contribution_flag=save_netcdf_fraction_as_contribution_flag_dummy
            endif

        enddo !in_region_loop

        if (trace_emissions_from_in_region) then
            if (allocated(subgrid_dummy)) deallocate (subgrid_dummy)

            ! THIS IS NO LONGER NEEDED, SINCE WE NO LONGER USE REGION LOOP FOR NO2 AND O3 CONTRIBUTION
            ! if (allocated(comp_subgrid_dummy)) deallocate (comp_subgrid_dummy)
            ! if (allocated(comp_source_subgrid_dummy)) deallocate(comp_source_subgrid_dummy)
            ! if (allocated(comp_source_EMEP_subgrid_dummy)) deallocate(comp_source_EMEP_subgrid_dummy)
            ! if (allocated(comp_source_EMEP_additional_subgrid_dummy)) deallocate(comp_source_EMEP_additional_subgrid_dummy)
        endif


        !Save the emissions interpolated to the target grid
        if (save_emissions) then
            write(unit_logfile,'(a)')'--------------------------'
            write(unit_logfile,'(a)')'Saving emissions'
            write(unit_logfile,'(a)')'--------------------------'
            variable_type='float'
            unit_str="g/s"
            do i_pollutant=1,n_pollutant_loop
                do i_source=1,n_source_index
                    if (calculate_source(i_source)) then

                        !Do not save emissions if the source is not traffic and the pollutant is road sand and salt or exhaust as these do not exist
                        if (i_source.ne.traffic_index.and.(pollutant_loop_index(i_pollutant).eq.pm25_sand_index.or.pollutant_loop_index(i_pollutant).eq.pm10_sand_index &
                            .or.pollutant_loop_index(i_pollutant).eq.pm25_salt_index.or.pollutant_loop_index(i_pollutant).eq.pm10_salt_index &
                            .or.pollutant_loop_index(i_pollutant).eq.pmex_index)) then
                            !Do nothing because these pollutants only exist for traffic
                        else

                            i_file=emission_file_index(i_source)
                            var_name_temp=trim(var_name_nc(conc_nc_index,pollutant_loop_index(i_pollutant),allsource_index))//'_'//trim(filename_grid(i_file))

                            !Calculate the emissions in the target grid
                            temp_subgrid=0.
                            do j=1,subgrid_dim(y_dim_index)
                                do i=1,subgrid_dim(x_dim_index)


                                    ii=crossreference_target_to_emission_subgrid(i,j,x_dim_index,i_source)
                                    jj=crossreference_target_to_emission_subgrid(i,j,y_dim_index,i_source)

                                    temp_subgrid(i,j,:)=emission_subgrid(ii,jj,:,i_source,i_pollutant)
                                    !Subgrid emissions, if relevant, are in units ug/sec/subgrid. Convert to g/s. Acount for the difference in subgrid sizes here
                                    temp_subgrid(i,j,:)=1.0e-6*temp_subgrid(i,j,:)*(subgrid_delta(y_dim_index)*subgrid_delta(x_dim_index)) &
                                        /(emission_subgrid_delta(y_dim_index,i_source)*emission_subgrid_delta(x_dim_index,i_source))

                                enddo
                            enddo

                            if (save_netcdf_file_flag) then
                                write(unit_logfile,'(a,es12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                                call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                    ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                            endif
                            if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                                write(unit_logfile,'(a,es12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                                call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                    ,unit_str,title_str_rec,create_file_rec,valid_min &
                                    ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                    ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                    ,z_rec(allsource_index,1) &
                                    ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                            endif

                        endif

                    endif
                enddo
            enddo
        endif

        !Save population interpolated to the target grid
        if (save_population) then
            variable_type='float'
            unit_str="inhabitants/grid"

            i_file=population_file_index(population_index)
            var_name_temp=trim(filename_grid(i_file))

            !Calculate the population in the target grid
            temp_subgrid=0.
            do j=1,subgrid_dim(y_dim_index)
                do i=1,subgrid_dim(x_dim_index)

                    ii=crossreference_target_to_population_subgrid(i,j,x_dim_index)
                    jj=crossreference_target_to_population_subgrid(i,j,y_dim_index)

                    temp_subgrid(i,j,:)=population_subgrid(ii,jj,population_index)
                    !Acount for the difference in subgrid sizes here
                    temp_subgrid(i,j,:)=temp_subgrid(i,j,:)*(subgrid_delta(y_dim_index)*subgrid_delta(x_dim_index)) &
                        /(population_subgrid_delta(y_dim_index)*population_subgrid_delta(x_dim_index))

                enddo
            enddo

            unit_str='inhabitants/grid'
            if (save_netcdf_file_flag) then
                write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_mask(temp_subgrid,use_subgrid(:,:,allsource_index),size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3))
                call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                    ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
            endif
            if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                    ,unit_str,title_str_rec,create_file_rec,valid_min &
                    ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,z_rec(allsource_index,1) &
                    ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
            endif

        endif

        !Save EMEP from in region distribution to the target grid
        !if (trace_emissions_from_in_region) save_emep_region_mask=.true.
        ! NO LONGER SAVE THIS (the arrays are outdated)
        save_emep_region_mask=.false.
        if (save_emep_region_mask) then
            variable_type='float'
            unit_str=""


            do k=1,2
                if (k.eq.1) then
                    var_name_temp=trim('EMEP_regional_fraction_mask')
                else
                    var_name_temp=trim('EMEP_additional_regional_fraction_mask')
                endif

                !Calculate the population in the target grid
                temp_subgrid=0.
                do j=1,subgrid_dim(y_dim_index)
                    do i=1,subgrid_dim(x_dim_index)

                        ii=crossreference_target_to_emep_subgrid(i,j,x_dim_index)
                        jj=crossreference_target_to_emep_subgrid(i,j,y_dim_index)

                        temp_subgrid(i,j,:)=EMEP_grid_fraction_in_region(ii,jj,allsource_index,k)
                        !temp_subgrid(i,j,:)=EMEP_grid_fraction_in_region(ii,jj,traffic_index,k)

                    enddo
                enddo

                unit_str='fraction'
                if (save_netcdf_file_flag) then
                    write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_mask(temp_subgrid,use_subgrid(:,:,allsource_index),size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3))
                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                endif
                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                    write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str_rec,create_file_rec,valid_min &
                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,z_rec(allsource_index,1) &
                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                endif
            enddo



        endif

        ! Save region mask at target subgrid
        if (trace_emissions_from_in_region) then  !! .or. use_region_select_and_mask_flag ??
            write(unit_logfile,'(A)') '-------------------------------'
            write(unit_logfile,'(A)') 'Saving region mask'
            write(unit_logfile,'(A)') '-------------------------------'

            variable_type = 'short'
            var_name_temp = 'region_id'
            unit_str='1'  !!!!!!! what to write here???
            temp_subgrid = 0.
            do tt = 1,subgrid_dim(t_dim_index)
                ! add 0.1 to all values to avoid them being rounded down to 1 less than original value
                temp_subgrid(:,:,tt) = nlreg_subgrid_region_id + 0.1
            end do
            if (save_netcdf_file_flag) then
                write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
                ! NB: hard-coding scale_factor=1 since it does not make sense to have scale_factor for an ID variable!
                call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                    ,unit_str,title_str,create_file,valid_min,variable_type,1.)
            endif
            if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(A)') 'Saving receptor version of region mask is not yet implemented!'
                stop
            endif
        end if

        ! Save contributions from nonlocal in-region
        if (save_emep_source_contributions .and. trace_emissions_from_in_region) then
            write(unit_logfile,'(A)') '-------------------------------'
            write(unit_logfile,'(A)') 'Saving EMEP semilocal contributions from-in-region'
            write(unit_logfile,'(A)') '-------------------------------'

            filename_append = '_from_in_region'
            do i_pollutant = 1,n_emep_pollutant_loop
                do i_source = 1,n_source_index
                    if (calculate_source(i_source).or.save_EMEP_source(i_source).or.i_source.eq.allsource_index) then
                        !Do not save for the additional GNFR sources
                        if (i_source == traffic_gasoline_nc_index .or. i_source == traffic_diesel_nc_index .or. i_source == traffic_gas_nc_index .or. i_source == publicpower_point_nc_index .or. i_source == publicpower_area_nc_index) cycle
                        !Do not save nonexhaust for nox
                        if (i_source == traffic_nonexhaust_index .and. (pollutant_loop_index(i_pollutant) == nox_index .or. pollutant_loop_index(i_pollutant) == no2_index)) cycle
                        ! loop over normal and additional version of the nonlocal in-region
                        do i_loop = 1,2
                            if (i_loop == 1) then
                                ! small domain
                                temp_subgrid = nlreg_subgrid_semilocal_from_in_region(:,:,:,i_source,i_pollutant)
                                i_file = nlreg_emep_subgrid_semilocal_file_index(i_source)
                            else
                                ! additional domain: add the additional increment
                                if (.not. (EMEP_additional_grid_interpolation_size > 0)) cycle
                                temp_subgrid = nlreg_subgrid_semilocal_from_in_region(:,:,:,i_source,i_pollutant) + nlreg_subgrid_semilocal_from_in_region_additional_increment(:,:,:,i_source,i_pollutant)
                                i_file = nlreg_emep_additional_subgrid_semilocal_file_index(i_source)
                            end if
                            variable_type='float'
                            unit_str='ug/m3'
                            var_name_temp=trim(var_name_nc(conc_nc_index,pollutant_loop_index(i_pollutant),allsource_index))//'_'//trim(filename_grid(i_file))//trim(filename_append)
                            if (save_netcdf_file_flag) then
                                write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(subgrid,1)/size(subgrid,2)/size(subgrid,3)
                                call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                    ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                    ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                            endif
                            if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                                write(unit_logfile,'(A)') 'Saving receptor version of semilocal not yet implemented!'
                                stop
                            endif
                        end do
                    end if
                end do
            end do
        end if

        !Save the EMEP data interpolated to the subgrid. These are based on the gridded concentrations
        ! Do not save additional contributions from-in-region
        if (save_emep_source_contributions) then

            if (trace_emissions_from_in_region) then
                n_in_region_loop=2
                if (.not.allocated(subgrid_dummy)) allocate (subgrid_dummy(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_subgrid_index,n_source_index,n_pollutant_loop))
                subgrid_dummy=0
            else
                n_in_region_loop=1
            endif

            do in_region_loop=1,n_in_region_loop

                write(unit_logfile,'(a)')'--------------------------'
                if (in_region_loop.eq.1) write(unit_logfile,'(a)')'Saving all EMEP contributions'
                if (in_region_loop.eq.2) write(unit_logfile,'(a)')'Saving only regional EMEP contributions'
                write(unit_logfile,'(a)')'--------------------------'

                !Loop over the normal subgrid and the in region subgrid by saving  to a dummy variable and pputting back when finished
                if (trace_emissions_from_in_region.and.in_region_loop.eq.2) then

                    filename_append='_from_in_region'
                    subgrid_dummy=subgrid
                    subgrid=subgrid_from_in_region
                    save_netcdf_fraction_as_contribution_flag_dummy=save_netcdf_fraction_as_contribution_flag
                    save_netcdf_fraction_as_contribution_flag=save_netcdf_fraction_as_contribution_from_in_region_flag
                else
                    filename_append=''
                endif

                do i_pollutant=1,n_emep_pollutant_loop
                    !EMEP does not have all pollutants so only save to n_emep_pollutant_loop
                    do i_source=1,n_source_index
                        if (calculate_source(i_source).or.save_EMEP_source(i_source).or.i_source.eq.allsource_index) then

                            !Only save the allsource value here 'EMEP_allsources'
                            if (i_source.eq.allsource_index) then
                                variable_type='float'
                                i_file=emep_subgrid_file_index(i_source)
                                var_name_temp=trim(var_name_nc(conc_nc_index,pollutant_loop_index(i_pollutant),allsource_index))//'_'//trim(filename_grid(i_file))//trim(filename_append)
                                unit_str='ug/m3'
                                if (save_netcdf_file_flag) then
                                    write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp),sum(subgrid(:,:,:,emep_subgrid_index,i_source,i_pollutant))/size(subgrid,1)/size(subgrid,2)/size(subgrid,3)
                                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                        ,subgrid(:,:,:,emep_subgrid_index,i_source,i_pollutant),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                        ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                                endif
                                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                                    write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(subgrid(:,:,:,emep_subgrid_index,i_source,i_pollutant))/size(subgrid,1)/size(subgrid,2)/size(subgrid,3)
                                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                        ,subgrid(:,:,:,emep_subgrid_index,i_source,i_pollutant),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                        ,unit_str,title_str_rec,create_file_rec,valid_min &
                                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                        ,z_rec(allsource_index,1) &
                                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                                endif

                            endif

                            !When pollutant is pmex then only save traffic
                            !if (pollutant_loop_index(i_pollutant).ne.pmex_nc_index.or.(pollutant_loop_index(i_pollutant).eq.pmex_nc_index.and.i_source.eq.traffic_nc_index)) then
                            !Do not save for the additional GNFR sources
                            if (i_source.ne.traffic_gasoline_nc_index.and.i_source.ne.traffic_diesel_nc_index.and.i_source.ne.traffic_gas_nc_index.and.i_source.ne.publicpower_point_nc_index.and.i_source.ne.publicpower_area_nc_index) then
                                !Do not save nonexhaust for no2 and nox
                                if (i_source.eq.traffic_nonexhaust_nc_index.and.(pollutant_loop_index(i_pollutant).eq.no2_index.or.pollutant_loop_index(i_pollutant).eq.nox_index.or.pollutant_loop_index(i_pollutant).eq.o3_index)) then
                                    !Do not save nonexhaust for exhaust gas emissions
                                else

                                    if (save_netcdf_fraction_as_contribution_flag) then
                                        variable_type='float'
                                        unit_str="ug/m3"
                                    else
                                        variable_type='byte'
                                        unit_str="%"
                                    endif

                                    i_file=emep_subgrid_local_file_index(i_source)
                                    var_name_temp=trim(var_name_nc(conc_nc_index,pollutant_loop_index(i_pollutant),allsource_index))//'_'//trim(filename_grid(i_file))//trim(filename_append)
                                    !if (i_source.eq.traffic_exhaust_nc_index) var_name_temp=var_name_temp//'_exhaust'
                                    !if (i_source.eq.traffic_nonexhaust_nc_index) var_name_temp=var_name_temp//'_nonexhaust'
                                    if (save_netcdf_fraction_as_contribution_flag) then
                                        temp_subgrid=subgrid(:,:,:,emep_local_subgrid_index,i_source,i_pollutant)
                                    else
                                        temp_subgrid=subgrid(:,:,:,emep_local_subgrid_index,i_source,i_pollutant)/subgrid(:,:,:,emep_subgrid_index,allsource_index,i_pollutant)*100.
                                        temp_subgrid=min(temp_subgrid,1000.)
                                        temp_subgrid=max(temp_subgrid,-1000.)
                                    endif
                                    if (save_netcdf_file_flag) then
                                        write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_mask(temp_subgrid,use_subgrid(:,:,allsource_index),size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3))
                                        call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                            ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                            ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                                    endif
                                    if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                                        write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                                        call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                            ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                            ,unit_str,title_str_rec,create_file_rec,valid_min &
                                            ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                            ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                            ,z_rec(allsource_index,1) &
                                            ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                                    endif

                                    !Save the additional EMEP source contributions here (not for from-in-region)
                                    if (EMEP_additional_grid_interpolation_size.gt.0 .and. in_region_loop == 1) then

                                        i_file=emep_additional_subgrid_local_file_index(i_source)
                                        var_name_temp=trim(var_name_nc(conc_nc_index,pollutant_loop_index(i_pollutant),allsource_index))//'_'//trim(filename_grid(i_file))//trim(filename_append)
                                        if (save_netcdf_fraction_as_contribution_flag) then
                                            temp_subgrid=subgrid(:,:,:,emep_additional_local_subgrid_index,i_source,i_pollutant)
                                        else
                                            temp_subgrid=subgrid(:,:,:,emep_additional_local_subgrid_index,i_source,i_pollutant)/subgrid(:,:,:,emep_subgrid_index,allsource_index,i_pollutant)*100.
                                            temp_subgrid=min(temp_subgrid,1000.)
                                            temp_subgrid=max(temp_subgrid,-1000.)
                                        endif
                                        if (save_netcdf_file_flag) then
                                            write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_mask(temp_subgrid,use_subgrid(:,:,allsource_index),size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3))
                                            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                                ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                                ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                                        endif
                                        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                                            write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                                            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                                ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                                ,unit_str,title_str_rec,create_file_rec,valid_min &
                                                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                                ,z_rec(allsource_index,1) &
                                                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                                        endif

                                    endif
                                endif
                            endif

                        endif
                    enddo
                enddo

                if (trace_emissions_from_in_region.and.in_region_loop.eq.2) then
                    subgrid_from_in_region=subgrid
                    subgrid=subgrid_dummy
                    save_netcdf_fraction_as_contribution_flag=save_netcdf_fraction_as_contribution_flag_dummy
                endif

            enddo !from_in_region loop
            if (trace_emissions_from_in_region) then
                if (allocated(subgrid_dummy)) deallocate (subgrid_dummy)
            endif
        endif


        !Save the other interpolated EMEP compounds used for nox chemistry as well. These are based on the surface comp values
        if (save_for_chemistry) then
            write(unit_logfile,'(a)')'--------------------------'
            write(unit_logfile,'(a)')'Saving for offline chemistry'
            write(unit_logfile,'(a)')'--------------------------'
            variable_type='float'
            do i_pollutant=1,n_emep_pollutant_loop
                if (pollutant_loop_index(i_pollutant).eq.nox_index) then
                    do i_loop=1,n_pollutant_compound_loop(i_pollutant)

                        i_comp=pollutant_compound_loop_index(i_pollutant,i_loop)

                        !Somo35 may be included here. Do not include it if it is.
                        if (i_comp.ne.somo35_nc_index) then

                            i_file=emep_subgrid_file_index(allsource_index)
                            var_name_temp=trim(var_name_nc(conc_nc_index,i_comp,allsource_index))//'_forchemistry'//'_'//trim(filename_grid(i_file))
                            temp_subgrid=comp_EMEP_subgrid(:,:,:,i_comp)
                            unit_str="ug/m3"
                            if (save_netcdf_file_flag) then
                                write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_mask(temp_subgrid,use_subgrid(:,:,allsource_index),size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3))
                                call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                    ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                            endif
                            if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                                write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(comp_EMEP_subgrid(:,:,:,i_comp))/size(subgrid,1)/size(subgrid,2)/size(subgrid,3)
                                call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                    ,unit_str,title_str_rec,create_file_rec,valid_min &
                                    ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                    ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                    ,z_rec(allsource_index,1) &
                                    ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                            endif

                        endif
                    enddo
                endif
            enddo
        endif

        !Save weighted travel time
        if (save_for_chemistry) then
            write(unit_logfile,'(a)')'--------------------------'
            write(unit_logfile,'(a)')'Saving weighted travel time'
            write(unit_logfile,'(a)')'--------------------------'
            variable_type='float'
            var_name_temp='Weighted_travel_time'
            unit_str='seconds'
            i_pollutant=pollutant_loop_back_index(nox_nc_index) !Only save the travel time for nox. Though this may be expanded for other compounds if necessary, like ammonia
            temp_subgrid=traveltime_subgrid(:,:,:,3,i_pollutant)
            if (save_netcdf_file_flag) then
                write(unit_logfile,'(a,es12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                !write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
                call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                    ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
            endif
            if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,es12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                !write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(traveltime_subgrid(:,:,:,1,i_pollutant))/size(subgrid,1)/size(subgrid,2)/size(subgrid,3)
                call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                    ,unit_str,title_str_rec,create_file_rec,valid_min &
                    ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,z_rec(allsource_index,1) &
                    ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
            endif
        endif

        if (save_for_chemistry) then
            variable_type='float'
            i_file=subgrid_J_file_index
            i_meteo=J_subgrid_index
            unit_str="1/s"
            !The same for all--------------------
            var_name_temp=trim(filename_grid(i_file))
            temp_subgrid=0.
            do j=1,subgrid_dim(y_dim_index)
                do i=1,subgrid_dim(x_dim_index)
                    ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
                enddo
            enddo
            if (save_netcdf_file_flag) then
                write(unit_logfile,'(a,es12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                !write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
                call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
            endif
            if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,es12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                !write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                    ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,z_rec(allsource_index,1) &
                    ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
            endif
            !The same for all--------------------

            !Save temperature if not also saving the meteo data.
            if (.not.save_other_meteo) then
                variable_type='float'
                i_file=subgrid_t2m_file_index
                i_meteo=t2m_subgrid_index
                unit_str="Celcius"
                !The same for all--------------------
                var_name_temp=trim(filename_grid(i_file))
                temp_subgrid=0.
                do j=1,subgrid_dim(y_dim_index)
                    do i=1,subgrid_dim(x_dim_index)
                        ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)-273.13
                    enddo
                enddo
                if (save_netcdf_file_flag) then
                    write(unit_logfile,'(a,es12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    !write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                endif
                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                    write(unit_logfile,'(a,es12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    !write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,z_rec(allsource_index,1) &
                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                endif
            endif
        endif

        if (save_emep_species) then
            write(unit_logfile,'(a)')'--------------------------'
            write(unit_logfile,'(a)')'Saving EMEP species'
            write(unit_logfile,'(a)')'--------------------------'
            variable_type='float'
            do i_sp=1,n_species_loop_index
                ii_sp=species_loop_index(i_sp)
                if (i_sp.ne.sp_pm_index) then
                    !Do not include the total
                    do i_loop=1,n_pmxx_sp_index
                        if (i_loop.eq.pm25_sp_index.or.i_loop.eq.pm10_sp_index) then

                            if (i_loop.eq.pm25_sp_index) i_pollutant=pollutant_loop_back_index(pm25_index)
                            if (i_loop.eq.pm10_sp_index) i_pollutant=pollutant_loop_back_index(pm10_index)

                            if (save_netcdf_fraction_as_contribution_flag) then
                                variable_type='float'
                                var_name_temp=trim(species_name_nc(i_loop,ii_sp))//'_nonlocal_contribution'
                                unit_str="ug/m3"
                                temp_subgrid=species_EMEP_subgrid(:,:,:,i_loop,i_sp)
                            else
                                variable_type='byte'
                                var_name_temp=trim(species_name_nc(i_loop,i_sp))//'_nonlocal_fraction'
                                unit_str="%"
                                temp_subgrid=species_EMEP_subgrid(:,:,:,i_loop,i_sp)/subgrid(:,:,:,total_subgrid_index,allsource_index,i_pollutant)*100.
                            endif

                            if (save_netcdf_file_flag) then
                                write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid(:,:,:))/size(subgrid,1)/size(subgrid,2)/size(subgrid,3)
                                call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                    ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                            endif
                            if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                                write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid(:,:,:))/size(subgrid,1)/size(subgrid,2)/size(subgrid,3)
                                call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                    ,unit_str,title_str_rec,create_file_rec,valid_min &
                                    ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                    ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                    ,z_rec(allsource_index,1) &
                                    ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                            endif
                        endif
                    enddo
                endif
            enddo
        endif

        !Only save sea salt in this way when emep species are not saved in the general way
        if (save_seasalt.and..not.save_emep_species) then
            write(unit_logfile,'(a)')'--------------------------'
            write(unit_logfile,'(a)')'Saving EMEP sea salt'
            write(unit_logfile,'(a)')'--------------------------'
            variable_type='float'
            i_sp=n_species_loop_index
            ii_sp=species_loop_index(i_sp)
            do i_loop=1,n_pmxx_sp_index
                if (i_loop.eq.pm25_sp_index.or.i_loop.eq.pm10_sp_index) then

                    if (i_loop.eq.pm25_sp_index) i_pollutant=pollutant_loop_back_index(pm25_index)
                    if (i_loop.eq.pm10_sp_index) i_pollutant=pollutant_loop_back_index(pm10_index)
                    if (save_netcdf_fraction_as_contribution_flag) then
                        variable_type='float'
                        var_name_temp=trim(species_name_nc(i_loop,ii_sp))//'_nonlocal_contribution'
                        unit_str="ug/m3"
                        temp_subgrid=species_EMEP_subgrid(:,:,:,i_loop,i_sp)
                    else
                        variable_type='byte'
                        var_name_temp=trim(species_name_nc(i_loop,ii_sp))//'_nonlocal_fraction'
                        unit_str="%"
                        temp_subgrid=species_EMEP_subgrid(:,:,:,i_loop,i_sp)/subgrid(:,:,:,total_subgrid_index,allsource_index,i_pollutant)*100.
                    endif

                    if (save_netcdf_file_flag) then
                        write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid(:,:,:))/size(subgrid,1)/size(subgrid,2)/size(subgrid,3)
                        call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                            ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                            ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                    endif
                    if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                        write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid(:,:,:))/size(subgrid,1)/size(subgrid,2)/size(subgrid,3)
                        call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                            ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                            ,unit_str,title_str_rec,create_file_rec,valid_min &
                            ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                            ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                            ,z_rec(allsource_index,1) &
                            ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                    endif
                endif
            enddo

        endif

        !Save the original EMEP compounds
        if (save_emep_original) then
            write(unit_logfile,'(a)')'--------------------------'
            write(unit_logfile,'(a)')'Saving original EMEP'
            write(unit_logfile,'(a)')'--------------------------'
            variable_type='float'
            do i_pollutant=1,n_emep_pollutant_loop
                if (pollutant_loop_index(i_pollutant).ne.pmex_index) then
                    do i_loop=1,n_pollutant_compound_loop(i_pollutant)

                        i_comp=pollutant_compound_loop_index(i_pollutant,i_loop)

                        i_file=emep_subgrid_file_index(allsource_index)
                        var_name_temp=trim(var_name_nc(conc_nc_index,i_comp,allsource_index))//'_original_EMEP_concentration'
                        unit_str="ug/m3"
                        if (i_comp.eq.somo35_index) unit_str="ppb�d"
                        if (save_netcdf_file_flag) then
                            !write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp),sum(orig_EMEP_subgrid(:,:,:,i_comp))/size(subgrid,1)/size(subgrid,2)/size(subgrid,3)
                            write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_mask(orig_EMEP_subgrid(:,:,:,i_comp),use_subgrid(:,:,allsource_index),size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3))
                            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                ,orig_EMEP_subgrid(:,:,:,i_comp),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                        endif
                        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                            ! write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(orig_EMEP_subgrid(:,:,:,i_comp))/size(subgrid,1)/size(subgrid,2)/size(subgrid,3)
                            write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_mask(orig_EMEP_subgrid(:,:,:,i_comp),use_subgrid(:,:,allsource_index),size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3))
                            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                ,orig_EMEP_subgrid(:,:,:,i_comp),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                ,unit_str,title_str_rec,create_file_rec,valid_min &
                                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                ,z_rec(allsource_index,1) &
                                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                        endif
                    enddo
                endif
            enddo
        endif

        !Save AQI
        if (save_aqi) then
            write(unit_logfile,'(a)')'--------------------------'
            write(unit_logfile,'(a)')'Saving AQI'
            write(unit_logfile,'(a)')'--------------------------'
            variable_type='byte'
            scale_factor=0.1
            !Hard coded AQI limits
            aqi_pollutant_index(1)=no2_index;aqi_pollutant_index(2)=pm10_index;aqi_pollutant_index(3)=pm25_index;aqi_pollutant_index(4)=o3_index
            aqi_limits_temp(:,2:4)=aqi_hourly_limits(:,1:3)
            aqi_limits_temp(:,1)=0.
            aqi_limits_temp(:,5)=2.*aqi_hourly_limits(:,3)
            aqi_subgrid=0.

            do t=1,subgrid_dim(t_dim_index)
                do j=1,subgrid_dim(y_dim_index)
                    do i=1,subgrid_dim(x_dim_index)
                        max_aqi=0.
                        do l=1,n_save_aqi_pollutant_index !pollutant (no2,pm10,pm2.5,o3)
                            do k=1,4 !level
                                if (comp_subgrid(i,j,t,aqi_pollutant_index(l)).ge.aqi_limits_temp(aqi_pollutant_index(l),k).and.comp_subgrid(i,j,t,aqi_pollutant_index(l)).lt.aqi_limits_temp(aqi_pollutant_index(l),k+1)) then
                                    aqi_subgrid(i,j,t,aqi_pollutant_index(l))=k+(comp_subgrid(i,j,t,aqi_pollutant_index(l))-aqi_limits_temp(aqi_pollutant_index(l),k))/(aqi_limits_temp(aqi_pollutant_index(l),k+1)-aqi_limits_temp(aqi_pollutant_index(l),k))
                                endif
                            enddo
                            aqi_subgrid(i,j,t,aqi_pollutant_index(l))=min(aqi_subgrid(i,j,t,aqi_pollutant_index(l)),4.99)
                            if (aqi_subgrid(i,j,t,aqi_pollutant_index(l)).gt.max_aqi) then
                                max_aqi=aqi_subgrid(i,j,t,aqi_pollutant_index(l))
                                aqi_responsible_pollutant_index(i,j,t)=l
                            endif
                            !write(*,*)  i,j,t,aqi_responsible_pollutant_index(i,j,t),aqi_subgrid(i,j,t,aqi_pollutant_index(l)),max_aqi
                        enddo
                    enddo
                enddo
            enddo

            do l=1,n_save_aqi_pollutant_index
                !write(unit_logfile,*)  'MAX AQI in time and space from '//trim(pollutant_file_str(aqi_pollutant_index(l)))//' = ',maxval(aqi_subgrid(:,:,:,aqi_pollutant_index(l)))
            enddo

            var_name_temp='AQI'
            unit_str='1'
            !Take the maximum of the pollutants
            temp_subgrid=maxval(aqi_subgrid(:,:,:,aqi_pollutant_index(1:n_save_aqi_pollutant_index)),4)/scale_factor
            if (save_netcdf_file_flag) then
                write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
                call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                    ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
            endif
            if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid*scale_factor)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                    ,unit_str,title_str_rec,create_file_rec,valid_min &
                    ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,z_rec(allsource_index,1) &
                    ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
            endif


            do l=1,n_save_aqi_pollutant_index

                i_comp=aqi_pollutant_index(l)
                var_name_temp='AQI_'//trim(var_name_nc(conc_nc_index,i_comp,allsource_index))
                unit_str='1'
                !Take the maximum of the pollutants
                temp_subgrid=aqi_subgrid(:,:,:,i_comp)/scale_factor
                if (save_netcdf_file_flag) then
                    write(unit_logfile,'(a)')'Writing netcdf variable: '//trim(var_name_temp)
                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                endif
                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                    write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid*scale_factor)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str_rec,create_file_rec,valid_min &
                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,z_rec(allsource_index,1) &
                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                endif

            enddo

            !Reset scale_factor
            scale_factor=1.

        endif

        !Save deposition
        if (save_deposition.and.calculate_deposition_flag) then
            write(unit_logfile,'(a)')'--------------------------'
            write(unit_logfile,'(a)')'Saving deposition'
            write(unit_logfile,'(a)')'--------------------------'
            variable_type='float'
            do i_pollutant=1,n_pollutant_loop
                do i_source=1,n_source_index
                    if (calculate_source(i_source)) then

                        i_comp=pollutant_loop_index(i_pollutant)

                        var_name_temp=trim(var_name_nc(conc_nc_index,i_comp,allsource_index))//'_local_dry_deposition_'//source_file_str(i_source)
                        unit_str="mg/m2/hr"
                        !Save in mg per hour, conversition from ug/m2/s
                        temp_subgrid=subgrid(:,:,:,drydepo_local_subgrid_index,i_source,i_pollutant)/1000.*3600.
                        if (save_netcdf_file_flag) then
                            write(unit_logfile,'(a,es12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                        endif
                        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                            write(unit_logfile,'(a,es12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                ,unit_str,title_str_rec,create_file_rec,valid_min &
                                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                ,z_rec(allsource_index,1) &
                                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                        endif

                        var_name_temp=trim(var_name_nc(conc_nc_index,i_comp,allsource_index))//'_local_wet_deposition_'//source_file_str(i_source)
                        unit_str="mg/m2/hr"
                        !Save in mg per hour, conversition from ug/m2/s
                        temp_subgrid=subgrid(:,:,:,wetdepo_local_subgrid_index,i_source,i_pollutant)/1000.*3600.
                        if (save_netcdf_file_flag) then
                            write(unit_logfile,'(a,es12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                            call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                        endif
                        if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                            write(unit_logfile,'(a,es12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                            call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                                ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                                ,unit_str,title_str_rec,create_file_rec,valid_min &
                                ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                                ,z_rec(allsource_index,1) &
                                ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                        endif

                    endif
                enddo

                var_name_temp=trim(var_name_nc(conc_nc_index,i_comp,allsource_index))//'_nonlocal_dry_deposition_'//source_file_str(allsource_index)
                unit_str="mg/m2/hr"
                !Save in mg per hour, conversition from ug/m2/s
                temp_subgrid=subgrid(:,:,:,drydepo_nonlocal_subgrid_index,allsource_index,i_pollutant)/1000.*3600.
                if (save_netcdf_file_flag) then
                    write(unit_logfile,'(a,es12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                endif
                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                    write(unit_logfile,'(a,es12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str_rec,create_file_rec,valid_min &
                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,z_rec(allsource_index,1) &
                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                endif

                var_name_temp=trim(var_name_nc(conc_nc_index,i_comp,allsource_index))//'_total_dry_deposition_'//source_file_str(allsource_index)
                unit_str="mg/m2/hr"
                !Save in mg per hour, conversition from ug/m2/s
                temp_subgrid=(subgrid(:,:,:,drydepo_nonlocal_subgrid_index,allsource_index,i_pollutant)+subgrid(:,:,:,drydepo_local_subgrid_index,allsource_index,i_pollutant))/1000.*3600.
                if (save_netcdf_file_flag) then
                    write(unit_logfile,'(a,es12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                endif
                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                    write(unit_logfile,'(a,es12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str_rec,create_file_rec,valid_min &
                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,z_rec(allsource_index,1) &
                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                endif

                var_name_temp=trim(var_name_nc(conc_nc_index,i_comp,allsource_index))//'_nonlocal_wet_deposition_'//source_file_str(allsource_index)
                unit_str="mg/m2/hr"
                !Save in mg per hour, conversition from ug/m2/s
                temp_subgrid=subgrid(:,:,:,wetdepo_nonlocal_subgrid_index,allsource_index,i_pollutant)/1000.*3600.
                if (save_netcdf_file_flag) then
                    write(unit_logfile,'(a,es12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                endif
                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                    write(unit_logfile,'(a,es12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str_rec,create_file_rec,valid_min &
                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,z_rec(allsource_index,1) &
                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                endif

                var_name_temp=trim(var_name_nc(conc_nc_index,i_comp,allsource_index))//'_total_wet_deposition_'//source_file_str(allsource_index)
                unit_str="mg/m2/hr"
                !Save in mg per hour, conversition from ug/m2/s
                temp_subgrid=(subgrid(:,:,:,wetdepo_nonlocal_subgrid_index,allsource_index,i_pollutant)+subgrid(:,:,:,wetdepo_local_subgrid_index,allsource_index,i_pollutant))/1000.*3600.
                if (save_netcdf_file_flag) then
                    write(unit_logfile,'(a,es12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                endif
                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                    write(unit_logfile,'(a,es12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str_rec,create_file_rec,valid_min &
                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,z_rec(allsource_index,1) &
                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                endif

                !Deposition velocity
                temp_subgrid=0.
                do j=1,subgrid_dim(y_dim_index)
                    do i=1,subgrid_dim(x_dim_index)
                        ii=crossreference_target_to_deposition_subgrid(i,j,x_dim_index);jj=crossreference_target_to_deposition_subgrid(i,j,y_dim_index)
                        temp_subgrid(i,j,:)=deposition_subgrid(ii,jj,:,vd_index,i_pollutant)
                    enddo
                enddo
                var_name_temp=trim(var_name_nc(conc_nc_index,i_comp,allsource_index))//'_deposition_velocity_'//source_file_str(allsource_index)
                unit_str="cm/s"
                !Save in cm/s, conversition from ug/m2/s
                temp_subgrid=temp_subgrid*100.
                if (save_netcdf_file_flag) then
                    write(unit_logfile,'(a,es12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                endif
                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                    write(unit_logfile,'(a,es12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str_rec,create_file_rec,valid_min &
                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,z_rec(allsource_index,1) &
                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                endif

                !Landuse category
                if (read_landuse_flag) then
                    temp_subgrid=0.
                    do j=1,subgrid_dim(y_dim_index)
                        do i=1,subgrid_dim(x_dim_index)
                            ii=crossreference_target_to_deposition_subgrid(i,j,x_dim_index);jj=crossreference_target_to_deposition_subgrid(i,j,y_dim_index)
                            temp_subgrid(i,j,:)=landuse_subgrid(ii,jj,grid_index)
                        enddo
                    enddo
                    var_name_temp=trim(var_name_nc(conc_nc_index,i_comp,allsource_index))//'_EMEP_landuse_category'
                    unit_str="1"
                    variable_type='byte'
                    if (save_netcdf_file_flag) then
                        write(unit_logfile,'(a,es12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                        call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                            ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                            ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                    endif
                    if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                        write(unit_logfile,'(a,es12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                        call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                            ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                            ,unit_str,title_str_rec,create_file_rec,valid_min &
                            ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                            ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                            ,z_rec(allsource_index,1) &
                            ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                    endif
                endif

                variable_type='float'

                !Save the original emep wet and dry deposition
                var_name_temp=trim(var_name_nc(conc_nc_index,i_comp,allsource_index))//'_original_EMEP_wet_deposition_'//source_file_str(allsource_index)
                unit_str="mg/m2/hr"
                !Save in mg per hour, conversition from ug/m2/s
                temp_subgrid=orig_emep_deposition_subgrid(:,:,:,wetdepo_index,i_pollutant)
                if (save_netcdf_file_flag) then
                    write(unit_logfile,'(a,es12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                endif
                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                    write(unit_logfile,'(a,es12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str_rec,create_file_rec,valid_min &
                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,z_rec(allsource_index,1) &
                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                endif
                var_name_temp=trim(var_name_nc(conc_nc_index,i_comp,allsource_index))//'_original_EMEP_dry_deposition_'//source_file_str(allsource_index)
                unit_str="mg/m2/hr"
                !Save in mg per hour, conversition from ug/m2/s
                temp_subgrid=orig_emep_deposition_subgrid(:,:,:,drydepo_index,i_pollutant)
                if (save_netcdf_file_flag) then
                    write(unit_logfile,'(a,es12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                endif
                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                    write(unit_logfile,'(a,es12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid,x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp &
                        ,unit_str,title_str_rec,create_file_rec,valid_min &
                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,z_rec(allsource_index,1) &
                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                endif

            enddo

        endif


        !Save the meteo interpolated to the target grid
        valid_min=-1.e24

        if (save_wind_vectors) then
            write(unit_logfile,'(a)')'--------------------------'
            write(unit_logfile,'(a)')'Saving wind vectors'
            write(unit_logfile,'(a)')'--------------------------'
            variable_type='float'
            if (wind_vectors_10m_available) then
                i_file=subgrid_u10_file_index
                i_meteo=u10_subgrid_index
            else
                i_file=subgrid_ugrid_file_index
                i_meteo=ugrid_subgrid_index
            endif

            unit_str="m/s"
            !The same for all--------------------
            var_name_temp=trim(filename_grid(i_file))
            temp_subgrid=0.
            do j=1,subgrid_dim(y_dim_index)
                do i=1,subgrid_dim(x_dim_index)
                    ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
                enddo
            enddo
            if (save_netcdf_file_flag) then
                write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp),mean_mask(temp_subgrid,use_subgrid(:,:,allsource_index),size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3))
                call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
            endif
            if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                    ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,z_rec(allsource_index,1) &
                    ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
            endif
            !The same for all--------------------

            if (wind_vectors_10m_available) then
                i_file=subgrid_v10_file_index
                i_meteo=v10_subgrid_index
            else
                i_file=subgrid_vgrid_file_index
                i_meteo=vgrid_subgrid_index
            endif
            unit_str="m/s"
            !The same for all--------------------
            var_name_temp=trim(filename_grid(i_file))
            temp_subgrid=0.
            do j=1,subgrid_dim(y_dim_index)
                do i=1,subgrid_dim(x_dim_index)
                    ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
                enddo
            enddo
            if (save_netcdf_file_flag) then
                write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp),mean_mask(temp_subgrid,use_subgrid(:,:,allsource_index),size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3))
                call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
            endif
            if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                    ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,z_rec(allsource_index,1) &
                    ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
            endif
            !The same for all--------------------
        endif


        if (save_other_meteo) then
            write(unit_logfile,'(a)')'--------------------------'
            write(unit_logfile,'(a)')'Saving other meteo data'
            write(unit_logfile,'(a)')'--------------------------'
            variable_type='float'
            i_file=subgrid_hmix_file_index
            i_meteo=hmix_subgrid_index
            unit_str="m"
            !The same for all--------------------
            var_name_temp=trim(filename_grid(i_file))
            temp_subgrid=0.
            do j=1,subgrid_dim(y_dim_index)
                do i=1,subgrid_dim(x_dim_index)
                    ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
                enddo
            enddo
            if (save_netcdf_file_flag) then
                write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_mask(temp_subgrid,use_subgrid(:,:,allsource_index),size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3))
                call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
            endif
            if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                    ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,z_rec(allsource_index,1) &
                    ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
            endif
            !The same for all--------------------

            variable_type='float'
            i_file=subgrid_t2m_file_index
            i_meteo=t2m_subgrid_index
            unit_str="Celcius"
            !The same for all--------------------
            var_name_temp=trim(filename_grid(i_file))
            temp_subgrid=0.
            do j=1,subgrid_dim(y_dim_index)
                do i=1,subgrid_dim(x_dim_index)
                    ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)-273.13
                enddo
            enddo
            if (save_netcdf_file_flag) then
                write(unit_logfile,'(a,f12.2)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
            endif
            if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.2)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                    ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,z_rec(allsource_index,1) &
                    ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
            endif
            !The same for all--------------------

            variable_type='float'
            i_file=subgrid_hmix_file_index
            i_meteo=precip_subgrid_index
            unit_str="mm/hr"
            !The same for all--------------------
            var_name_temp='precipitation'!trim(filename_grid(i_file))
            temp_subgrid=0.
            do j=1,subgrid_dim(y_dim_index)
                do i=1,subgrid_dim(x_dim_index)
                    ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
                enddo
            enddo
            if (save_netcdf_file_flag) then
                write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_mask(temp_subgrid,use_subgrid(:,:,allsource_index),size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3))
                call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
            endif
            if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                    ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,z_rec(allsource_index,1) &
                    ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
            endif
            !The same for all--------------------


            i_file=subgrid_kz_file_index
            i_meteo=kz_subgrid_index
            unit_str="m2/s"
            !The same for all--------------------
            var_name_temp=trim(filename_grid(i_file))
            temp_subgrid=0.
            do j=1,subgrid_dim(y_dim_index)
                do i=1,subgrid_dim(x_dim_index)
                    ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
                enddo
            enddo
            if (save_netcdf_file_flag) then
                write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_mask(temp_subgrid,use_subgrid(:,:,allsource_index),size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3))
                call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
            endif
            if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                    ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,z_rec(allsource_index,1) &
                    ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
            endif
            !The same for all--------------------

            if (hourly_calculations) then

                if (wind_vectors_10m_available) then
                    i_file=subgrid_FF10_file_index
                    i_meteo=FF10_subgrid_index
                else
                    i_file=subgrid_FFgrid_file_index
                    i_meteo=FFgrid_subgrid_index
                endif
                unit_str="m/s"
                !The same for all--------------------
                var_name_temp=trim(filename_grid(i_file))
                temp_subgrid=0.
                do j=1,subgrid_dim(y_dim_index)
                    do i=1,subgrid_dim(x_dim_index)
                        ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
                    enddo
                enddo
                if (save_netcdf_file_flag) then
                    write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_mask(temp_subgrid,use_subgrid(:,:,allsource_index),size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3))
                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                endif
                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                    write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,z_rec(allsource_index,1) &
                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                endif
                !The same for all--------------------

                if (wind_vectors_10m_available) then
                    i_file=subgrid_DD10_file_index
                else
                    i_file=subgrid_DDgrid_file_index
                endif
                !i_meteo=DDgrid_subgrid_index !i_meteo not used in this conversion to wind direction
                unit_str="degrees"
                do tt=1,subgrid_dim(t_dim_index)
                    do jj=1,integral_subgrid_dim(y_dim_index)
                        do ii=1,integral_subgrid_dim(x_dim_index)
                            if (wind_vectors_10m_available) then
                                temp_integral_subgrid(ii,jj,tt)=DIRECTION(meteo_subgrid(ii,jj,tt,u10_subgrid_index),meteo_subgrid(ii,jj,tt,v10_subgrid_index))
                            else
                                temp_integral_subgrid(ii,jj,tt)=DIRECTION(meteo_subgrid(ii,jj,tt,ugrid_subgrid_index),meteo_subgrid(ii,jj,tt,vgrid_subgrid_index))
                            endif

                        enddo
                    enddo
                enddo
                !This is not converted to real degrees, but subgrid degrees, so is wrong but close

                !The same for all--------------------
                var_name_temp=trim(filename_grid(i_file))
                temp_subgrid=0.
                do j=1,subgrid_dim(y_dim_index)
                    do i=1,subgrid_dim(x_dim_index)
                        ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);
                        !In this case this is not the same for all
                        temp_subgrid(i,j,:)=temp_integral_subgrid(ii,jj,:)
                    enddo
                enddo
                if (save_netcdf_file_flag) then
                    write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_mask(temp_subgrid,use_subgrid(:,:,allsource_index),size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3))
                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                endif
                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                    write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,z_rec(allsource_index,1) &
                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                endif
                !The same for all--------------------

            endif

            i_file=subgrid_invL_file_index
            i_meteo=invL_subgrid_index
            unit_str="1/m"
            !The same for all--------------------
            var_name_temp=trim(filename_grid(i_file))
            temp_subgrid=0.
            do j=1,subgrid_dim(y_dim_index)
                do i=1,subgrid_dim(x_dim_index)
                    ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
                enddo
            enddo
            if (save_netcdf_file_flag) then
                write(unit_logfile,'(a,f12.4)')'Writing netcdf variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
            endif
            if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.4)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                    ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,z_rec(allsource_index,1) &
                    ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
            endif
            !The same for all--------------------

            i_file=subgrid_logz0_file_index
            i_meteo=logz0_subgrid_index
            unit_str="log(m)"
            !The same for all--------------------
            var_name_temp=trim(filename_grid(i_file))
            temp_subgrid=0.
            do j=1,subgrid_dim(y_dim_index)
                do i=1,subgrid_dim(x_dim_index)
                    ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
                enddo
            enddo
            if (save_netcdf_file_flag) then
                write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_mask(temp_subgrid,use_subgrid(:,:,allsource_index),size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3))
                call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
            endif
            if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                    ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,z_rec(allsource_index,1) &
                    ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
            endif
            !The same for all--------------------

            i_file=subgrid_ustar_file_index
            i_meteo=ustar_subgrid_index
            unit_str="m/s"
            !The same for all--------------------
            var_name_temp=trim(filename_grid(i_file))
            temp_subgrid=0.
            do j=1,subgrid_dim(y_dim_index)
                do i=1,subgrid_dim(x_dim_index)
                    ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
                enddo
            enddo
            if (save_netcdf_file_flag) then
                write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_mask(temp_subgrid,use_subgrid(:,:,allsource_index),size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3))
                call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
            endif
            if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                    ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                    ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                    ,z_rec(allsource_index,1) &
                    ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
            endif
            !The same for all--------------------

            if (annual_calculations) then

                i_file=subgrid_invFFgrid_file_index
                i_meteo=inv_FFgrid_subgrid_index
                unit_str="s/m"
                !The same for all--------------------
                var_name_temp=trim(filename_grid(i_file))
                temp_subgrid=0.
                do j=1,subgrid_dim(y_dim_index)
                    do i=1,subgrid_dim(x_dim_index)
                        ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
                    enddo
                enddo
                if (save_netcdf_file_flag) then
                    write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_mask(temp_subgrid,use_subgrid(:,:,allsource_index),size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3))
                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                endif
                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                    write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,z_rec(allsource_index,1) &
                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                endif
                !The same for all--------------------

                i_file=subgrid_invFF10_file_index
                i_meteo=inv_FF10_subgrid_index
                unit_str="s/m"
                !The same for all--------------------
                var_name_temp=trim(filename_grid(i_file))
                temp_subgrid=0.
                do j=1,subgrid_dim(y_dim_index)
                    do i=1,subgrid_dim(x_dim_index)
                        ii=crossreference_target_to_integral_subgrid(i,j,x_dim_index);jj=crossreference_target_to_integral_subgrid(i,j,y_dim_index);temp_subgrid(i,j,:)=meteo_subgrid(ii,jj,:,i_meteo)
                    enddo
                enddo
                if (save_netcdf_file_flag) then
                    write(unit_logfile,'(a,f12.3)')'Writing netcdf variable: '//trim(var_name_temp), mean_mask(temp_subgrid,use_subgrid(:,:,allsource_index),size(temp_subgrid,1),size(temp_subgrid,2),size(temp_subgrid,3))
                    call uEMEP_save_netcdf_file(unit_logfile,temp_name,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str,create_file,valid_min,variable_type,scale_factor)
                endif
                if (save_netcdf_receptor_flag.and.n_valid_receptor.ne.0) then
                    write(unit_logfile,'(a,f12.3)')'Writing netcdf receptor variable: '//trim(var_name_temp),sum(temp_subgrid)/size(temp_subgrid,1)/size(temp_subgrid,2)/size(temp_subgrid,3)
                    call uEMEP_save_netcdf_receptor_file(unit_logfile,temp_name_rec,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index) &
                        ,temp_subgrid(:,:,:),x_subgrid,y_subgrid,lon_subgrid,lat_subgrid,var_name_temp,unit_str,title_str_rec,create_file_rec,valid_min &
                        ,x_receptor(valid_receptor_index(1:n_valid_receptor)),y_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,lon_receptor(valid_receptor_index(1:n_valid_receptor)),lat_receptor(valid_receptor_index(1:n_valid_receptor)) &
                        ,z_rec(allsource_index,1) &
                        ,name_receptor(valid_receptor_index(1:n_valid_receptor),1),n_valid_receptor,variable_type,scale_factor)
                endif
                !The same for all--------------------

            endif

        endif

        !Deallocate the averaging arrays with each receptor grid
        if (use_multiple_receptor_grids_flag.and..not.use_single_time_loop_flag.and.save_netcdf_average_flag) then
            if (allocated(val_array_av)) deallocate(val_array_av)
            if (allocated(time_seconds_output_av)) deallocate(time_seconds_output_av)
            counter_av=0
        endif


    end subroutine uEMEP_save_netcdf_control


    subroutine uEMEP_save_netcdf_file(unit_logfile_in,filename_netcdf,nx,ny,nt_in,val_array_in,x_array,y_array,lon_array,lat_array,name_array,unit_array,title_str,create_file,valid_min,variable_type,scale_factor)

        use uEMEP_definitions
        use netcdf

        implicit none

        character(256) filename_netcdf,name_array,unit_array,title_str,temp_name
        integer unit_logfile_in
        integer nx,ny,nt_in
        real val_array(nx,ny,nt_in),val_array_in(nx,ny,nt_in)!,val_array_temp(nx,ny,nt)
        real x_array(nx,ny)
        real y_array(nx,ny)
        real lon_array(nx,ny)
        real lat_array(nx,ny) !,lat_array_temp(nx,ny)
        !real time_array(nt)
        real x_vector(nx)
        real y_vector(ny)
        logical create_file
        real valid_min
        character(256) variable_type
        real scale_factor

        integer ncid
        integer y_dimid,x_dimid,time_dimid
        integer y_varid,x_varid,lat_varid,lon_varid,val_varid,time_varid,proj_varid
        integer dimids3(3),dimids2(2),chunks3(3)
        integer n_dims(3)
        integer status
        integer nf90_type
        integer t
        integer n_time_total
        integer nt
        integer(8) time_seconds_output_nc(nt_in) !Need integer 8 for the averaging
        integer i_source
        character(2) temp_str
        integer i
        character(256) temp_str2

        if (trim(variable_type).eq.'byte') nf90_type=NF90_BYTE
        if (trim(variable_type).eq.'short') nf90_type=NF90_SHORT
        if (trim(variable_type).eq.'float') nf90_type=NF90_FLOAT
        if (trim(variable_type).eq.'double') nf90_type=NF90_DOUBLE

        !Assumes x and y are the dimensions
        x_vector=x_array(:,1)
        y_vector=y_array(1,:)

        n_time_total=end_time_nc_index-start_time_nc_index+1
        nt=nt_in
        val_array=val_array_in
        time_seconds_output_nc=time_seconds_output

        !Save averages only
        if (save_netcdf_average_flag) then
            counter_av=counter_av+1
            if (counter_av.gt.n_var_av) then
                write(unit_logfile_in,*) 'ERROR: Array size for saving averages (n_var_av) not large enough. Stopping'
                stop
            endif
            if (use_single_time_loop_flag) then
                val_array_av(:,:,counter_av)=val_array_av(:,:,counter_av)+val_array(:,:,nt) !nt=1 in this case
                time_seconds_output_av(counter_av)=time_seconds_output_av(counter_av)+time_seconds_output_nc(nt)
                if (t_loop.eq.end_time_loop_index) then
                    val_array(:,:,nt)=val_array_av(:,:,counter_av)/n_time_total
                    time_seconds_output_nc(nt)=time_seconds_output_av(counter_av)/n_time_total
                endif
                !write(unit_logfile_in,'(a,3i)') 'Saving as average single time loop (nt,counter_av):',nt,counter_av,time_seconds_output_nc(nt)
            else
                val_array_av(:,:,counter_av)=sum(val_array,3)/size(val_array,3)
                time_seconds_output_av(counter_av)=sum(time_seconds_output_nc)/size(time_seconds_output_nc,1)
                nt=1
                val_array(:,:,nt)=val_array_av(:,:,counter_av)
                time_seconds_output_nc(nt)=time_seconds_output_av(counter_av)
                !write(unit_logfile_in,'(a,3i)') 'Saving as average multiple time loop (nt,counter_av):',nt,counter_av,time_seconds_output_nc(nt)

            endif

        endif

        !Mask the regions if required
        if (use_region_select_and_mask_flag) then
            do t=1,nt
                where (use_subgrid_val(:,:,allsource_index).eq.outside_region_index) val_array(:,:,t)=NODATA_value
                where (.not.use_subgrid(:,:,allsource_index)) val_array(:,:,t)=NODATA_value
            enddo
        endif

        if (create_file) then
            !Create a netcdf file
            !call check(  nf90_create(filename_netcdf, nf90_clobber, ncid) )
            !call check(  nf90_create(filename_netcdf, NF90_HDF5, ncid) )
            call check(  nf90_create(filename_netcdf, IOR(NF90_HDF5, NF90_CLASSIC_MODEL), ncid) ) !New

            !Specify global attributes
            call check(  nf90_put_att(ncid, nf90_global, "Conventions", "CF-1.6" ) )
            call check(  nf90_put_att(ncid, nf90_global, "title", trim(title_str)) )
            call check(  nf90_put_att(ncid, nf90_global, "Model", trim(model_version_str) ) )


            !Save some model flags for later reference
            call check(  nf90_put_att(ncid, nf90_global, "EMEP_grid_interpolation_flag", EMEP_grid_interpolation_flag ) )
            call check(  nf90_put_att(ncid, nf90_global, "no2_chemistry_scheme_flag", no2_chemistry_scheme_flag ) )
            call check(  nf90_put_att(ncid, nf90_global, "EMEP_grid_interpolation_size", EMEP_grid_interpolation_size ) )
            call check(  nf90_put_att(ncid, nf90_global, "EMEP_additional_grid_interpolation_size", EMEP_additional_grid_interpolation_size ) )
            call check(  nf90_put_att(ncid, nf90_global, "no2_background_chemistry_scheme_flag", no2_background_chemistry_scheme_flag ) )
            call check(  nf90_put_att(ncid, nf90_global, "local_subgrid_method_flag", local_subgrid_method_flag ) )
            call check(  nf90_put_att(ncid, nf90_global, "EMEP_emission_grid_interpolation_flag", EMEP_emission_grid_interpolation_flag ) )
            if (limit_emep_grid_interpolation_region_to_calculation_region) then
                call check(  nf90_put_att(ncid, nf90_global, "limit_emep_grid_interpolation_region_to_calculation_region", "true" ) )
            else
                call check(  nf90_put_att(ncid, nf90_global, "limit_emep_grid_interpolation_region_to_calculation_region", "false" ) )
            endif
            call check(  nf90_put_att(ncid, nf90_global, "n_local_fraction_grids", n_local_fraction_grids ) )
            do i=1,n_local_fraction_grids
                write (temp_str,'(i0)') i
                temp_str2="local_fraction_grid_size("//trim(temp_str)//")"
                call check(  nf90_put_att(ncid, nf90_global, trim(temp_str2), local_fraction_grid_size(i)) )
            enddo
            call check(  nf90_put_att(ncid, nf90_global, "local_fraction_grid_for_EMEP_grid_interpolation", local_fraction_grid_for_EMEP_grid_interpolation ) )
            call check(  nf90_put_att(ncid, nf90_global, "local_fraction_grid_for_EMEP_additional_grid_interpolation", local_fraction_grid_for_EMEP_additional_grid_interpolation ) )
            if (.not.use_GNFR_emissions_from_EMEP_flag) then
                call check(  nf90_put_att(ncid, nf90_global, "use_GNFR_emissions_from_EMEP_flag", "false" ) )
            endif
            if (use_GNFR_emissions_from_EMEP_flag.and..not.use_GNFR19_emissions_from_EMEP_flag) then
                call check(  nf90_put_att(ncid, nf90_global, "use_GNFR13_emissions_from_EMEP_flag", "true" ) )
            endif
            if (use_GNFR19_emissions_from_EMEP_flag) then
                call check(  nf90_put_att(ncid, nf90_global, "use_GNFR19_emissions_from_EMEP_flag", "true" ) )
            endif


            !Write out sources
            i=0
            do i_source=1,n_source_index
                if (calculate_source(i_source)) then
                    i=i+1
                    write (temp_str,'(i0)') i
                    temp_str2="uEMEP_source("//trim(temp_str)//")"
                    call check(  nf90_put_att(ncid, nf90_global, trim(temp_str2), trim(source_file_str(i_source)) ) )
                endif
            enddo
            call check(  nf90_put_att(ncid, nf90_global, "n_uEMEP_sources", i ))
            i=0
            do i_source=1,n_source_index
                if (calculate_EMEP_source(i_source)) then
                    i=i+1
                    write (temp_str,'(i0)') i
                    temp_str2="EMEP_source("//trim(temp_str)//")"
                    call check(  nf90_put_att(ncid, nf90_global, trim(temp_str2), trim(source_file_str(i_source)) ) )
                endif
            enddo
            call check(  nf90_put_att(ncid, nf90_global, "n_EMEP_sources", i ))

            !Projection data
            if (projection_type.eq.UTM_projection_index) then
                call check(  nf90_def_var(ncid, "projection_utm", NF90_int, proj_varid) )
                call check(  nf90_put_att(ncid, proj_varid, "semi_major_axis", 6378137.0 ) )
                call check(  nf90_put_att(ncid, proj_varid, "inverse_flattening", 298.257222101 ) )

                call check(  nf90_put_att(ncid, proj_varid, "grid_mapping_name", "transverse_mercator" ) )
                call check(  nf90_put_att(ncid, proj_varid, "scale_factor_at_central_meridian", 0.9996 ) )
                call check(  nf90_put_att(ncid, proj_varid, "latitude_of_projection_origin", 0 ) )
                call check(  nf90_put_att(ncid, proj_varid, "false_easting", 500000. ) )
                if (utm_zone.ge.0) then
                    call check(  nf90_put_att(ncid, proj_varid, "false_northing", 0. ) )
                else
                    call check(  nf90_put_att(ncid, proj_varid, "false_northing", 10000000. ) )
                endif
                call check(  nf90_put_att(ncid, proj_varid, "longitude_of_central_meridian", utm_lon0 ) )
            endif

            if (projection_type.eq.LTM_projection_index) then
                call check(  nf90_def_var(ncid, "projection_tm", NF90_int, proj_varid) )
                call check(  nf90_put_att(ncid, proj_varid, "semi_major_axis", 6378137.0 ) )
                call check(  nf90_put_att(ncid, proj_varid, "inverse_flattening", 298.257222101 ) )

                call check(  nf90_put_att(ncid, proj_varid, "grid_mapping_name", "transverse_mercator" ) )
                call check(  nf90_put_att(ncid, proj_varid, "scale_factor_at_central_meridian", 0.9996 ) )
                call check(  nf90_put_att(ncid, proj_varid, "latitude_of_projection_origin", 0 ) )
                call check(  nf90_put_att(ncid, proj_varid, "false_easting", 500000. ) )
                if (utm_zone.ge.0) then
                    call check(  nf90_put_att(ncid, proj_varid, "false_northing", 0. ) )
                else
                    call check(  nf90_put_att(ncid, proj_varid, "false_northing", 10000000. ) )
                endif
                call check(  nf90_put_att(ncid, proj_varid, "longitude_of_central_meridian", ltm_lon0 ) )

            endif

            if (projection_type.eq.RDM_projection_index) then
                !Not properly assigned, same as UTM
                call check(  nf90_def_var(ncid, "projection_RDM", NF90_int, proj_varid) )
                call check(  nf90_put_att(ncid, proj_varid, "semi_major_axis", 6378137.0 ) )
                call check(  nf90_put_att(ncid, proj_varid, "inverse_flattening", 298.257222101 ) )

                call check(  nf90_put_att(ncid, proj_varid, "grid_mapping_name", "transverse_mercator" ) )
                call check(  nf90_put_att(ncid, proj_varid, "scale_factor_at_central_meridian", 0.9996 ) )
                call check(  nf90_put_att(ncid, proj_varid, "latitude_of_projection_origin", 0 ) )
                call check(  nf90_put_att(ncid, proj_varid, "false_easting", 500000. ) )
                call check(  nf90_put_att(ncid, proj_varid, "false_northing", 0. ) )
                call check(  nf90_put_att(ncid, proj_varid, "longitude_of_central_meridian", utm_lon0 ) )
            endif

            if (projection_type.eq.LAEA_projection_index) then
                call check(  nf90_def_var(ncid, "projection_ETRS89_LAEA", NF90_int, proj_varid) )
                !call check(  nf90_put_att(ncid, proj_varid, "semi_major_axis", 6378137.0 ) )
                !call check(  nf90_put_att(ncid, proj_varid, "inverse_flattening", 298.257222101 ) )

                !https://github.com/mdsumner/rasterwise/blob/master/README.md
                !EPSG:3035
                !int ETRS89-LAEA ;
                !    ETRS89-LAEA:missing_value = -1. ;
                !   ETRS89-LAEA:grid_mapping_name = "lambert_azimuthal_equal_area" ;
                !   ETRS89-LAEA:longitude_of_projection_origin = 10. ;
                !   ETRS89-LAEA:latitude_of_projection_origin = 52. ;
                !   ETRS89-LAEA:false_easting = 4321000. ;
                !   ETRS89-LAEA:false_northing = 3210000. ;
                !   ETRS89-LAEA:inverse_flattening = 298.257222101 ;
                !   ETRS89-LAEA:semi_major_axis = 6378137. ;

                !call check(  nf90_put_att(ncid, proj_varid, "grid_mapping_name", "lambert_azimuthal_equal_area" ) )
                !call check(  nf90_put_att(ncid, proj_varid, "semi_major_axis", 6378137.0 ) )
                !call check(  nf90_put_att(ncid, proj_varid, "inverse_flattening", 298.257222101 ) )
                !call check(  nf90_put_att(ncid, proj_varid, "scale_factor_at_central_meridian", 0.9996 ) )
                !call check(  nf90_put_att(ncid, proj_varid, "latitude_of_projection_origin",  52. ) )
                !call check(  nf90_put_att(ncid, proj_varid, "false_easting", 4321000. ) )
                !call check(  nf90_put_att(ncid, proj_varid, "false_northing", 3210000. ) )
                !call check(  nf90_put_att(ncid, proj_varid, "longitude_of_projection_origin", 10.) )

                call check(  nf90_put_att(ncid, proj_varid, "grid_mapping_name", "lambert_azimuthal_equal_area" ) )
                call check(  nf90_put_att(ncid, proj_varid, "semi_major_axis", projection_attributes(5) ) )
                call check(  nf90_put_att(ncid, proj_varid, "inverse_flattening", projection_attributes(6) ) )
                call check(  nf90_put_att(ncid, proj_varid, "latitude_of_projection_origin",  projection_attributes(2) ) )
                call check(  nf90_put_att(ncid, proj_varid, "false_easting", projection_attributes(3) ) )
                call check(  nf90_put_att(ncid, proj_varid, "false_northing", projection_attributes(4) ) )
                call check(  nf90_put_att(ncid, proj_varid, "longitude_of_projection_origin", projection_attributes(1)) )
            endif

            !Define the dimensions
            call check(  nf90_def_dim(ncid,"time",NF90_UNLIMITED, time_dimid) )
            call check(  nf90_def_dim(ncid, "y", ny, y_dimid) )
            call check(  nf90_def_dim(ncid, "x", nx, x_dimid) )

            !Define the dimension variables
            !call check(  nf90_def_var(ncid, "time", NF90_DOUBLE, time_dimid, time_varid) )
            call check(  nf90_def_var(ncid, "time", NF90_INT, time_dimid, time_varid) )
            call check(  nf90_def_var(ncid, "y", NF90_REAL, y_dimid, y_varid) )
            call check(  nf90_def_var(ncid, "x", NF90_REAL, x_dimid, x_varid) )

            !Define the values
            dimids3 = (/ x_dimid, y_dimid, time_dimid /)
            dimids2 = (/ x_dimid, y_dimid /)
            call check(  nf90_def_var(ncid, "lat", NF90_REAL, dimids2, lat_varid) )
            call check(  nf90_def_var(ncid, "lon", NF90_REAL, dimids2, lon_varid) )

            !Specify the units
            call check(  nf90_put_att(ncid, lat_varid, "units", "degrees_north") )
            call check(  nf90_put_att(ncid, lon_varid, "units", "degrees_east") )
            call check(  nf90_put_att(ncid, y_varid, "units", "m") )
            call check(  nf90_put_att(ncid, x_varid, "units", "m") )
            call check(  nf90_put_att(ncid, time_varid, "units", trim(unit_dim_nc(time_dim_nc_index))) )

            !Specify other dimension attributes
            call check(  nf90_put_att(ncid, y_varid, "standard_name", "projection_y_coordinate") )
            call check(  nf90_put_att(ncid, x_varid, "standard_name", "projection_x_coordinate") )
            call check(  nf90_put_att(ncid, y_varid, "axis", "Y") )
            call check(  nf90_put_att(ncid, x_varid, "axis", "X") )

            !Close the definitions
            call check( nf90_enddef(ncid) )

            !call check( nf90_put_var(ncid, time_varid, val_dim_nc(1:dim_length_nc(time_dim_nc_index),time_dim_nc_index)) )
            !call check( nf90_put_var(ncid, time_varid, time_seconds_output(1:dim_length_nc(time_dim_nc_index))) )
            call check( nf90_put_var(ncid, time_varid, time_seconds_output_nc(1:nt)) )
            call check( nf90_put_var(ncid, y_varid, y_vector) )
            call check( nf90_put_var(ncid, x_varid, x_vector) )
            call check( nf90_put_var(ncid, lat_varid, lat_array) )
            call check( nf90_put_var(ncid, lon_varid, lon_array) )

            call check( nf90_close(ncid) )

        endif

        !Do not save any average data for single time loops until the last time step where it will save the average
        if (save_netcdf_average_flag.and.use_single_time_loop_flag.and.t_loop.ne.end_time_loop_index) then
            !call check( nf90_close(ncid) )
            return
        endif

        !Add to the existing file
        call check( nf90_open(filename_netcdf, NF90_WRITE, ncid) )

        !Get the dimensions id from the existing file
        call check( nf90_inq_dimid(ncid,"time",time_dimid) )
        call check( nf90_inq_dimid(ncid, "y", y_dimid) )
        call check( nf90_inq_dimid(ncid, "x", x_dimid) )
        dimids3 = (/ x_dimid, y_dimid, time_dimid /)
        chunks3 = (/ nx, ny, 1 /) !New
        call check( nf90_inquire_dimension(ncid, dimids3(1), temp_name, n_dims(1)) )
        call check( nf90_inquire_dimension(ncid, dimids3(2), temp_name, n_dims(2)) )
        call check( nf90_inquire_dimension(ncid, dimids3(3), temp_name, n_dims(3)) )


        status=nf90_inq_varid(ncid, trim(name_array), val_varid)
        if (status.ne.nf90_NoErr) then
            call check( nf90_redef(ncid) )
            !if the variable does not exist then create a new one
            !write(*,*) 'Creating new: ',trim(name_array)
            call check( nf90_def_var(ncid, trim(name_array), nf90_type, dimids3, val_varid) )
            ! gzip level 3 compression and shuffling
            ! optional _FillValue for values which never have been written, unpacked value
            call check( nf90_def_var_chunking(ncid, val_varid, NF90_CHUNKED, chunks3) ) !New
            call check( nf90_def_var_deflate(ncid, val_varid, 1, 1, 3) ) !New
            call check( nf90_put_att(ncid, val_varid, "units", trim(unit_array)) )

            !Specify other variable attributes
            if (nf90_type.eq.NF90_byte) then
                call check(  nf90_put_att(ncid, val_varid, "missing_value", int(NODATA_value,kind=1) ) ) !New
                call check(  nf90_put_att(ncid, val_varid, "valid_min", int(valid_min,kind=1)) )
            elseif (nf90_type.eq.NF90_short) then
                call check(  nf90_put_att(ncid, val_varid, "missing_value", int(NODATA_value,kind=2) ) ) !New
                call check(  nf90_put_att(ncid, val_varid, "valid_min", int(valid_min,kind=2)) )
            else
                call check(  nf90_put_att(ncid, val_varid, "missing_value", NODATA_value ) ) !New
                call check(  nf90_put_att(ncid, val_varid, "valid_min", valid_min) )
            endif
            if (projection_type.eq.UTM_projection_index) then
                call check(  nf90_put_att(ncid, val_varid, "grid_mapping", "projection_utm") )
            elseif (projection_type.eq.LTM_projection_index) then
                call check(  nf90_put_att(ncid, val_varid, "grid_mapping", "projection_tm") )
            elseif (projection_type.eq.LAEA_projection_index) then
                call check(  nf90_put_att(ncid, val_varid, "grid_mapping", "projection_ETRS89_LAEA") )
            elseif (projection_type.eq.RDM_projection_index) then
                call check(  nf90_put_att(ncid, val_varid, "grid_mapping", "projection_RDM") )
            endif

            call check(  nf90_put_att(ncid, val_varid, "coordinates", "lon lat") )
            if (scale_factor.ne.1.) call check(  nf90_put_att(ncid, val_varid, "scale_factor", scale_factor) )

            !Close the definitions
            call check( nf90_enddef(ncid) )

        endif


        if (use_single_time_loop_flag) then
            !Add time to the time dimension
            call check( nf90_inq_varid(ncid, "time", time_varid) )
            !call check( nf90_inquire_dimension(ncid, time_dimid, temp_name, n_dims(3)) )
            !n_dims(3)=n_dims(3)+1
            if (save_netcdf_average_flag) then
                n_dims(3)=nt
            else
                n_dims(3)=t_loop
            endif
            call check( nf90_put_var(ncid, time_varid, time_seconds_output_nc(1), start = (/n_dims(3)/) ) )
            !write(*,*) n_dims(3),val_dim_nc(1,time_dim_nc_index)
            !write(*,*) n_dims

            !Add dimension and array to existing
            call check( nf90_inq_varid(ncid, trim(name_array), val_varid) )
            if (nf90_type.eq.NF90_byte) then
                call check( nf90_put_var(ncid, val_varid, int(val_array(:,:,1:nt),kind=1), start=(/1,1,n_dims(3)/), count=(/n_dims(1),n_dims(2),1/)) )
            elseif (nf90_type.eq.NF90_short) then
                call check( nf90_put_var(ncid, val_varid, int(val_array(:,:,1:nt),kind=2), start=(/1,1,n_dims(3)/), count=(/n_dims(1),n_dims(2),1/)) )
            else
                call check( nf90_put_var(ncid, val_varid, val_array, start=(/1,1,n_dims(3)/), count=(/n_dims(1),n_dims(2),1/)) )
            endif

        else

            !Write the whole variable to file. Default is float
            if (nf90_type.eq.NF90_byte) then
                call check( nf90_put_var(ncid, val_varid, int(val_array(:,:,1:nt),kind=1)) )
            elseif (nf90_type.eq.NF90_short) then
                call check( nf90_put_var(ncid, val_varid, int(val_array(:,:,1:nt),kind=2)) )
            else
                call check( nf90_put_var(ncid, val_varid, val_array(:,:,1:nt)) )
            endif

        endif

        call check( nf90_close(ncid) )



    end subroutine uEMEP_save_netcdf_file


    subroutine check(status)

        use netcdf
        implicit none
        integer, intent ( in) :: status

        if(status /= nf90_noerr) then
            write(*,*) 'Stopping due to netcdf error: '//trim(nf90_strerror(status))
            error stop
        end if

    end subroutine check


! ######################################################################
    FUNCTION DIRECTION(UD,VD)
!	CALCULATES THE WIND DIRECTION
        !Taken from NBLM1
        IMPLICIT NONE
        REAL DIRECTION,UD,VD,PI
        PI=180./3.14159
        DIRECTION=0.
        IF (UD.GT.0.AND.VD.GE.0) DIRECTION=270.-ATAN(ABS(VD/UD))*PI
        IF (UD.LE.0.AND.VD.GT.0) DIRECTION=180.-ATAN(ABS(UD/VD))*PI
        IF (UD.LT.0.AND.VD.LE.0) DIRECTION=90.-ATAN(ABS(VD/UD))*PI
        IF (UD.GE.0.AND.VD.LT.0) DIRECTION=360.-ATAN(ABS(UD/VD))*PI
    END FUNCTION DIRECTION
! ######################################################################

    function mean_nodata(array,n1,n2,n3,nodata_num)

        !use uEMEP_definitions
        implicit none
        real :: array(n1,n2,n3)
        real :: nodata_num
        integer :: n1,n2,n3
        real :: mean_nodata
        integer i,j,t
        real :: count=0
        real :: sum_array=0

        do t=1,n3
            do j=1,n2
                do i=1,n1
                    if (array(i,j,t).ne.nodata_num) then
                        sum_array=sum_array+array(i,j,t)
                        count=count+1.
                    endif
                enddo
            enddo
        enddo

        if (count.gt.0) then
            mean_nodata=sum_array/count
        else
            mean_nodata=0
        endif


    end function mean_nodata

    function mean_mask(array,mask,n1,n2,n3)

        !use uEMEP_definitions
        implicit none
        real, intent (in) :: array(n1,n2,n3)
        logical, intent (in) :: mask(n1,n2)
        integer, intent (in) :: n1,n2,n3
        real :: mean_mask
        integer i,j,t
        real :: count=0
        real :: sum_array=0

        sum_array=0
        count=0

        do t=1,n3
            do j=1,n2
                do i=1,n1
                    if (mask(i,j)) then
                        sum_array=sum_array+array(i,j,t)
                        count=count+1.
                    endif
                enddo
            enddo
        enddo

        if (count.gt.0) then
            mean_mask=sum_array/count
        else
            mean_mask=0
        endif


    end function mean_mask


!Saves receptor data in netcdf format

    subroutine uEMEP_save_netcdf_receptor_file(unit_logfile_in,filename_netcdf,nx,ny,nt_in,val_array_in,x_array,y_array,lon_array,lat_array,name_array,unit_array,title_str,create_file,valid_min &
        ,x_rec,y_rec,lon_rec,lat_rec,height_rec,name_rec_in,nr,variable_type,scale_factor)

        use uEMEP_definitions
        use netcdf

        implicit none

        character(256) filename_netcdf,name_array,unit_array,title_str,temp_name
        integer unit_logfile_in
        integer nx,ny,nt_in,nr
        real val_array(nx,ny,nt_in),val_array_in(nx,ny,nt_in)!,val_array_temp(nx,ny,nt_in)
        real x_array(nx,ny)
        real y_array(nx,ny)
        real lon_array(nx,ny)
        real lat_array(nx,ny)!,lat_array_temp(nx,ny)
        !real time_array(nt_in)
        !real x_vector(nx)
        !real y_vector(ny)
        logical create_file
        real valid_min
        character(256) variable_type
        real scale_factor

        integer ncid
        integer station_dimid,time_dimid,charlen_dimid
        integer station_varid,station_name_varid,lat_varid,lon_varid,val_varid,time_varid,proj_varid,x_varid,y_varid,height_varid
        integer dimids2(2)
        integer n_dims_length(3),n_dims_start(3)
        integer status
        integer tr,rr
        real x_rec(nr),y_rec(nr),height_rec(nr)
        real lon_rec(nr),lat_rec(nr)
        character(256) name_rec_in(nr)
        character(256) temp_char
        integer n_char
        !parameter (n_char=7)
        parameter (n_char=64)
        character(1) name_rec(n_char,nr)
        integer n_time_total
        real val_rec(nr,nt_in)
        real delta(2)
        integer id_rec(nr)
        integer nf90_type
        integer nt
        integer(8) time_seconds_output_nc(nt_in)
        integer tr_0
        integer i_source
        character(2) temp_str
        integer i
        character(256) temp_str2

        !Do not save if no receptor position data is available
        if (nr.eq.0) then
            return
        endif

        if (trim(variable_type).eq.'byte') nf90_type=NF90_BYTE
        if (trim(variable_type).eq.'short') nf90_type=NF90_SHORT
        if (trim(variable_type).eq.'float') nf90_type=NF90_FLOAT
        if (trim(variable_type).eq.'double') nf90_type=NF90_DOUBLE

        nt=nt_in
        val_array=val_array_in
        time_seconds_output_nc=time_seconds_output

        !Save averages only
        if (save_netcdf_average_flag) then
            counter_av=counter_av+1
            if (counter_av.gt.n_var_av) then
                write(unit_logfile_in,*) 'ERROR: Array size for saving averages (n_var_av) not large enough. Stopping'
                stop
            endif
            if (use_single_time_loop_flag) then
                val_array_av(:,:,counter_av)=val_array_av(:,:,counter_av)+val_array(:,:,nt) !nt=1 in this case
                time_seconds_output_av(counter_av)=time_seconds_output_av(counter_av)+time_seconds_output_nc(nt)
                if (t_loop.eq.end_time_loop_index) then
                    val_array(:,:,nt)=val_array_av(:,:,counter_av)/end_time_loop_index
                    time_seconds_output_nc(nt)=time_seconds_output_av(counter_av)/end_time_loop_index
                endif
                !write(unit_logfile_in,*) 'Saving as average single time loop (nt,counter_av):',nt,counter_av
            else
                !write(unit_logfile_in,*) 'Saving as average multiple time loop (nt,counter_av):',nt,counter_av,time_seconds_output_nc(1),time_seconds_output_nc(nt)
                !write(*,*) time_seconds_output_nc(1:nt)
                val_array_av(:,:,counter_av)=sum(val_array(:,:,1:nt),3)/nt
                time_seconds_output_av(counter_av)=sum(time_seconds_output_nc(1:nt),1)/nt
                !write(*,*) sum(time_seconds_output_nc(1:nt),1)
                nt=1
                val_array(:,:,nt)=val_array_av(:,:,counter_av)
                time_seconds_output_nc(nt)=time_seconds_output_av(counter_av)
                !write(unit_logfile_in,*) 'Saving as average multiple time loop (nt,counter_av):',nt,counter_av,time_seconds_output_nc(nt),time_seconds_output_av(counter_av)
            endif

        endif

        !Interpolate to receptor position given the input array
        !write(unit_logfile,'(a)')' Interpolating to receptor point '
        !Assumes 2 elements to an array
        delta(1)=(x_array(2,1)-x_array(1,1))
        delta(2)=(y_array(1,2)-y_array(1,1))
        do rr=1,nr
            do tr=1,nt
                val_rec(rr,tr)=area_weighted_interpolation_function(x_array,y_array,val_array(:,:,tr),nx,ny,delta,x_rec(rr),y_rec(rr))
            enddo
        enddo

        !Make the receptor name to fit to netcdf requirements
        do rr=1,nr
            id_rec(rr)=rr
            temp_char=name_rec_in(rr)
            tr_0=len_trim(temp_char)
            do tr=1,tr_0
                name_rec(tr,rr)=temp_char(tr:tr)
                !write(*,*) trim(name_rec_in(rr)),trim(name_rec(rr))
            enddo
            do tr=tr_0+1,n_char
                name_rec(tr,rr)=char(0)
            enddo
        enddo


        if (save_netcdf_average_flag) then
            n_time_total=1
        else
            n_time_total=end_time_nc_index-start_time_nc_index+1
        endif



        if (create_file) then
            !Create a netcdf file
            !call check(  nf90_create(filename_netcdf, nf90_clobber, ncid) )
            !call check(  nf90_create(filename_netcdf, NF90_HDF5, ncid) )
            call check(  nf90_create(filename_netcdf, IOR(NF90_HDF5, NF90_CLASSIC_MODEL), ncid) ) !New

            !Specify global attributes
            call check(  nf90_put_att(ncid, nf90_global, "Conventions", "CF-1.6" ) )
            call check(  nf90_put_att(ncid, nf90_global, "title", trim(title_str)) )
            call check(  nf90_put_att(ncid, nf90_global, "Model", trim(model_version_str) ) )
            call check(  nf90_put_att(ncid, nf90_global, "featureType", "timeSeries" ) )

            !Save some model flags for later reference
            call check(  nf90_put_att(ncid, nf90_global, "EMEP_grid_interpolation_flag", EMEP_grid_interpolation_flag ) )
            call check(  nf90_put_att(ncid, nf90_global, "no2_chemistry_scheme_flag", no2_chemistry_scheme_flag ) )
            call check(  nf90_put_att(ncid, nf90_global, "EMEP_grid_interpolation_size", EMEP_grid_interpolation_size ) )
            call check(  nf90_put_att(ncid, nf90_global, "EMEP_additional_grid_interpolation_size", EMEP_additional_grid_interpolation_size ) )
            call check(  nf90_put_att(ncid, nf90_global, "no2_background_chemistry_scheme_flag", no2_background_chemistry_scheme_flag ) )
            call check(  nf90_put_att(ncid, nf90_global, "local_subgrid_method_flag", local_subgrid_method_flag ) )
            call check(  nf90_put_att(ncid, nf90_global, "EMEP_emission_grid_interpolation_flag", EMEP_emission_grid_interpolation_flag ) )
            if (limit_emep_grid_interpolation_region_to_calculation_region) then
                call check(  nf90_put_att(ncid, nf90_global, "limit_emep_grid_interpolation_region_to_calculation_region", "true" ) )
            else
                call check(  nf90_put_att(ncid, nf90_global, "limit_emep_grid_interpolation_region_to_calculation_region", "false" ) )
            endif
            call check(  nf90_put_att(ncid, nf90_global, "n_local_fraction_grids", n_local_fraction_grids ) )
            do i=1,n_local_fraction_grids
                write (temp_str,'(i0)') i
                temp_str2="local_fraction_grid_size("//trim(temp_str)//")"
                call check(  nf90_put_att(ncid, nf90_global, trim(temp_str2), local_fraction_grid_size(i)) )
            enddo
            call check(  nf90_put_att(ncid, nf90_global, "local_fraction_grid_for_EMEP_grid_interpolation", local_fraction_grid_for_EMEP_grid_interpolation ) )
            call check(  nf90_put_att(ncid, nf90_global, "local_fraction_grid_for_EMEP_additional_grid_interpolation", local_fraction_grid_for_EMEP_additional_grid_interpolation ) )
            if (.not.use_GNFR_emissions_from_EMEP_flag) then
                call check(  nf90_put_att(ncid, nf90_global, "use_GNFR_emissions_from_EMEP_flag", "false" ) )
            endif
            if (use_GNFR_emissions_from_EMEP_flag.and..not.use_GNFR19_emissions_from_EMEP_flag) then
                call check(  nf90_put_att(ncid, nf90_global, "use_GNFR13_emissions_from_EMEP_flag", "true" ) )
            endif
            if (use_GNFR19_emissions_from_EMEP_flag) then
                call check(  nf90_put_att(ncid, nf90_global, "use_GNFR19_emissions_from_EMEP_flag", "true" ) )
            endif

            !Write out sources
            i=0
            do i_source=1,n_source_index
                if (calculate_source(i_source)) then
                    i=i+1
                    write (temp_str,'(i0)') i
                    temp_str2="uEMEP_source("//trim(temp_str)//")"
                    call check(  nf90_put_att(ncid, nf90_global, trim(temp_str2), trim(source_file_str(i_source)) ) )
                endif
            enddo
            call check(  nf90_put_att(ncid, nf90_global, "n_uEMEP_sources", i ))
            i=0
            do i_source=1,n_source_index
                if (calculate_EMEP_source(i_source)) then
                    i=i+1
                    write (temp_str,'(i0)') i
                    temp_str2="EMEP_source("//trim(temp_str)//")"
                    call check(  nf90_put_att(ncid, nf90_global, trim(temp_str2), trim(source_file_str(i_source)) ) )
                endif
            enddo
            call check(  nf90_put_att(ncid, nf90_global, "n_EMEP_sources", i ))

            !Projection data
            if (projection_type.eq.UTM_projection_index) then
                call check(  nf90_def_var(ncid, "projection_utm", NF90_int, proj_varid) )
                call check(  nf90_put_att(ncid, proj_varid, "semi_major_axis", 6378137.0 ) )
                call check(  nf90_put_att(ncid, proj_varid, "inverse_flattening", 298.257222101 ) )

                call check(  nf90_put_att(ncid, proj_varid, "grid_mapping_name", "transverse_mercator" ) )
                call check(  nf90_put_att(ncid, proj_varid, "scale_factor_at_central_meridian", 0.9996 ) )
                call check(  nf90_put_att(ncid, proj_varid, "latitude_of_projection_origin", 0 ) )
                call check(  nf90_put_att(ncid, proj_varid, "false_easting", 500000. ) )
                call check(  nf90_put_att(ncid, proj_varid, "false_northing", 0. ) )
                call check(  nf90_put_att(ncid, proj_varid, "longitude_of_central_meridian", utm_lon0 ) )
                !call check(  nf90_put_att(ncid, proj_varid, "semi_major_axis", 6378140.0 ) )
                !call check(  nf90_put_att(ncid, proj_varid, "semi_minor_axis", 6356750.0 ) )
            endif

            if (projection_type.eq.LTM_projection_index) then
                call check(  nf90_def_var(ncid, "projection_tm", NF90_int, proj_varid) )
                call check(  nf90_put_att(ncid, proj_varid, "semi_major_axis", 6378137.0 ) )
                call check(  nf90_put_att(ncid, proj_varid, "inverse_flattening", 298.257222101 ) )

                call check(  nf90_put_att(ncid, proj_varid, "grid_mapping_name", "transverse_mercator" ) )
                call check(  nf90_put_att(ncid, proj_varid, "scale_factor_at_central_meridian", 0.9996 ) )
                call check(  nf90_put_att(ncid, proj_varid, "latitude_of_projection_origin", 0 ) )
                call check(  nf90_put_att(ncid, proj_varid, "false_easting", 500000. ) )
                call check(  nf90_put_att(ncid, proj_varid, "false_northing", 0. ) )
                call check(  nf90_put_att(ncid, proj_varid, "longitude_of_central_meridian", ltm_lon0 ) )
            endif

            if (projection_type.eq.RDM_projection_index) then
                !Not properly assigned, same as UTM
                call check(  nf90_def_var(ncid, "projection_RDM", NF90_int, proj_varid) )
                call check(  nf90_put_att(ncid, proj_varid, "semi_major_axis", 6378137.0 ) )
                call check(  nf90_put_att(ncid, proj_varid, "inverse_flattening", 298.257222101 ) )

                call check(  nf90_put_att(ncid, proj_varid, "grid_mapping_name", "transverse_mercator" ) )
                call check(  nf90_put_att(ncid, proj_varid, "scale_factor_at_central_meridian", 0.9996 ) )
                call check(  nf90_put_att(ncid, proj_varid, "latitude_of_projection_origin", 0 ) )
                call check(  nf90_put_att(ncid, proj_varid, "false_easting", 500000. ) )
                call check(  nf90_put_att(ncid, proj_varid, "false_northing", 0. ) )
                call check(  nf90_put_att(ncid, proj_varid, "longitude_of_central_meridian", utm_lon0 ) )
            endif

            if (projection_type.eq.LAEA_projection_index) then
                call check(  nf90_def_var(ncid, "projection_ETRS89_LAEA", NF90_int, proj_varid) )
                !call check(  nf90_put_att(ncid, proj_varid, "semi_major_axis", 6378137.0 ) )
                !call check(  nf90_put_att(ncid, proj_varid, "inverse_flattening", 298.257222101 ) )

                !https://github.com/mdsumner/rasterwise/blob/master/README.md
                !EPSG:3035
                call check(  nf90_put_att(ncid, proj_varid, "grid_mapping_name", "lambert_azimuthal_equal_area" ) )
                call check(  nf90_put_att(ncid, proj_varid, "semi_major_axis", projection_attributes(5) ) )
                call check(  nf90_put_att(ncid, proj_varid, "inverse_flattening", projection_attributes(6) ) )
                call check(  nf90_put_att(ncid, proj_varid, "latitude_of_projection_origin",  projection_attributes(2) ) )
                call check(  nf90_put_att(ncid, proj_varid, "false_easting", projection_attributes(3) ) )
                call check(  nf90_put_att(ncid, proj_varid, "false_northing", projection_attributes(4) ) )
                call check(  nf90_put_att(ncid, proj_varid, "longitude_of_projection_origin", projection_attributes(1)) )
            endif

            !Define the dimensions for the entire dataset
            !write(*,*) 'n_valid_receptor_in',n_valid_receptor_in
            !write(*,*) 'n_valid_receptor',n_valid_receptor
            call check(  nf90_def_dim(ncid,"station_id",n_valid_receptor_in, station_dimid) )
            call check(  nf90_def_dim(ncid,"charlen",n_char, charlen_dimid) )

            !call check(  nf90_def_dim(ncid,"time",n_time_total, time_dimid) )
            !To have time as unlimittec (Heiko)
            call check(  nf90_def_dim(ncid,"time",NF90_UNLIMITED, time_dimid) )


            !Define the dimension variables
            call check(  nf90_def_var(ncid, "station_id", NF90_INT, station_dimid, station_varid) )
            !call check(  nf90_def_var(ncid, "time", NF90_DOUBLE, time_dimid, time_varid) )
            call check(  nf90_def_var(ncid, "time", NF90_INT, time_dimid, time_varid) )

            !Define the values
            dimids2 = (/ station_dimid, time_dimid /)
            !dimids1 = (/ station_dimid /)
            call check(  nf90_def_var(ncid, "lat", NF90_REAL, station_dimid, lat_varid) )
            call check(  nf90_def_var(ncid, "lon", NF90_REAL, station_dimid, lon_varid) )
            call check(  nf90_def_var(ncid, "y", NF90_REAL, station_dimid, y_varid) )
            call check(  nf90_def_var(ncid, "x", NF90_REAL, station_dimid, x_varid) )
            call check(  nf90_def_var(ncid, "station_name", NF90_CHAR, (/charlen_dimid,station_dimid/), station_name_varid) )
            call check(  nf90_def_var(ncid, "station_height", NF90_REAL, station_dimid, height_varid) )
            !call check(  nf90_def_var(ncid, "station_name", NF90_CHAR, (/station_dimid/), station_name_varid) )
            !call check(  nf90_def_var(ncid, "station_name", NF90_CHAR, (/charlen_dimid,station_dimid/), station_name_varid) )

            !Specify the units
            call check(  nf90_put_att(ncid, lat_varid, "units", "degrees_north") )
            call check(  nf90_put_att(ncid, lon_varid, "units", "degrees_east") )
            call check(  nf90_put_att(ncid, y_varid, "units", "m") )
            call check(  nf90_put_att(ncid, x_varid, "units", "m") )
            call check(  nf90_put_att(ncid, height_varid, "units", "m") )
            call check(  nf90_put_att(ncid, time_varid, "units", trim(unit_dim_nc(time_dim_nc_index))) )
            call check(  nf90_put_att(ncid, station_varid, "long_name", "station index" ) )
            call check(  nf90_put_att(ncid, station_name_varid, "long_name", "station name" ) )
            call check(  nf90_put_att(ncid, station_name_varid, "cf_role", "timeseries_id" ) )

            !Specify other dimension attributes
            !call check(  nf90_put_att(ncid, y_varid, "standard_name", "projection_y_coordinate") )
            !call check(  nf90_put_att(ncid, x_varid, "standard_name", "projection_x_coordinate") )
            !call check(  nf90_put_att(ncid, y_varid, "axis", "Y") )
            !call check(  nf90_put_att(ncid, x_varid, "axis", "X") )

            !Close the definitions
            call check( nf90_enddef(ncid) )

            !call check( nf90_put_var(ncid, time_varid, val_dim_nc(1:dim_length_nc(time_dim_nc_index),time_dim_nc_index)) )
            !call check( nf90_put_var(ncid, station_varid, name_rec(:,1:n_char) )
            !call check( nf90_put_var(ncid, lat_varid, lat_array) )
            !call check( nf90_put_var(ncid, lon_varid, lon_array) )

            !Put this in to fix unlimitted problem? If works on annual need to check it works other places as well!
            !call check( nf90_put_var(ncid, time_varid, time_seconds_output_nc(1:nt), start=(/1/), count=(/nt/)) )
            !call check( nf90_put_var(ncid, time_varid, time_seconds_output_nc(1:nt)) )

            call check( nf90_close(ncid) )

        endif

        !Add to the existing file
        call check( nf90_open(filename_netcdf, NF90_WRITE, ncid) )

        !Get the dimensions id from the existing file
        call check( nf90_inq_dimid(ncid,"time",time_dimid) )
        call check( nf90_inq_dimid(ncid, "station_id", station_dimid) )
        ! write(*,*) 'station_dimid ',station_dimid
        ! write(*,*) 'time_dimid ',time_dimid
        dimids2 = (/ station_dimid, time_dimid /)

        !Get the size of the dimensions
        call check( nf90_inquire_dimension(ncid, dimids2(1), temp_name, n_dims_length(1)) )
        call check( nf90_inquire_dimension(ncid, dimids2(2), temp_name, n_dims_length(2)) )
        !Set the starting point to 1
        n_dims_start(1:2)=1


        !Set time to full length in unlimitted case (Heiko)
        n_dims_length(2) = n_time_total

        ! write(*,*) 'n_dims_length(1) ',n_dims_length(1)
        ! write(*,*) 'n_dims_length(2) ',n_dims_length(2)

        status=nf90_inq_varid(ncid, trim(name_array), val_varid)
        if (status.ne.nf90_NoErr) then
            !if the variable does not exist then create a new one
            !write(*,*) 'Creating new: ',trim(name_array)
            call check( nf90_redef(ncid) )

            call check( nf90_def_var(ncid, trim(name_array), nf90_type, dimids2, val_varid) )
            call check( nf90_put_att(ncid, val_varid, "units", trim(unit_array)) )

            !Specify other variable attributes
            if (nf90_type.eq.NF90_byte) then
                call check(  nf90_put_att(ncid, val_varid, "missing_value", int(NODATA_value,kind=1) ) ) !New
                call check(  nf90_put_att(ncid, val_varid, "valid_min", int(valid_min,kind=1)) )
            elseif (nf90_type.eq.NF90_short) then
                call check(  nf90_put_att(ncid, val_varid, "missing_value", int(NODATA_value,kind=2) ) ) !New
                call check(  nf90_put_att(ncid, val_varid, "valid_min", int(valid_min,kind=2)) )
            else
                call check(  nf90_put_att(ncid, val_varid, "missing_value", NODATA_value ) ) !New
                call check(  nf90_put_att(ncid, val_varid, "valid_min", valid_min) )
            endif
            call check(  nf90_put_att(ncid, val_varid, "grid_mapping", "projection_utm") )
            call check(  nf90_put_att(ncid, val_varid, "coordinates", "station_name lat lon") )
            if (scale_factor.ne.1.) call check(  nf90_put_att(ncid, val_varid, "scale_factor", scale_factor) )

            !Close the definitions
            call check( nf90_enddef(ncid) )
        endif



        !Should not be used but can be
        if (use_single_time_loop_flag) then
            !Add time to the time dimension
            !call check( nf90_inq_varid(ncid, "time", time_varid) )
            !call check( nf90_inquire_dimension(ncid, time_dimid, temp_name, n_dims(3)) )
            !n_dims(3)=n_dims(3)+1

            if (save_netcdf_average_flag) then
                n_dims_start(2)=1
                n_dims_length(2)=1
            else
                n_dims_start(2)=t_loop
                n_dims_length(2)=1
            endif
            !call check( nf90_put_var(ncid, time_varid, val_dim_nc(1,time_dim_nc_index), start = (/n_dims(2)/) ) )
            !write(*,*) n_dims(3),val_dim_nc(1,time_dim_nc_index)
            !write(*,*) n_dims

            !Add dimension and array to existing
            !call check( nf90_inq_varid(ncid, trim(name_array), val_varid) )
            !call check( nf90_put_var(ncid, val_varid, val_rec, start=(/1,n_dims(2)/), count=(/n_dims(1),1/)) )

            !call check( nf90_put_var(ncid, time_varid, val_dim_nc(1:dim_length_nc(time_dim_nc_index),time_dim_nc_index), start=(/n_dims(2)/), count=(/1/)) )
            !call check( nf90_put_var(ncid, station_varid, name_rec(1,:), start=(/1,1/), count=(/n_dims(1),n_char/)) )

        elseif (use_multiple_receptor_grids_flag) then
            !Add time to the time dimension
            !call check( nf90_inq_varid(ncid, "station", station_varid) )
            !call check( nf90_inquire_dimension(ncid, time_dimid, temp_name, n_dims(3)) )
            !n_dims(3)=n_dims(3)+1

            n_dims_start(1)=valid_receptor_inverse_index(g_loop)
            n_dims_length(1)=1
            id_rec(1)=valid_receptor_inverse_index(g_loop)
            !write(*,*) n_dims(3),val_dim_nc(1,time_dim_nc_index)
            !write(*,*) n_dims

            !Add dimension and array to existing
            !call check( nf90_inq_varid(ncid, trim(name_array), val_varid) )
            !call check( nf90_put_var(ncid, val_varid, val_rec, start=(/n_dims(1),1/), count=(/1,n_dims(2)/)) )

            !call check( nf90_put_var(ncid, time_varid, val_dim_nc(1:dim_length_nc(time_dim_nc_index),time_dim_nc_index) , start=(/1/), count=(/n_dims(2)/)) )
            !call check( nf90_put_var(ncid, station_varid, name_rec(1,:), start = (/n_dims(1),1/), count=(/1,n_char/)) )
        endif

        !write(*,*) 'n_dims_start',n_dims_start
        !write(*,*) 'n_dims_length',n_dims_length

        !Fill in the complete dimension variables for time and receptor names
        call check( nf90_inq_varid(ncid, "station_id", station_varid) )
        call check( nf90_inq_varid(ncid, "time", time_varid) )
        call check( nf90_inq_varid(ncid, "station_name", station_name_varid) )
        call check( nf90_inq_varid(ncid, "x", x_varid) )
        call check( nf90_inq_varid(ncid, "y", y_varid) )
        call check( nf90_inq_varid(ncid, "lon", lon_varid) )
        call check( nf90_inq_varid(ncid, "lat", lat_varid) )
        call check( nf90_inq_varid(ncid, "station_height", height_varid) )

        !Write time to the file

        !!call check( nf90_put_var(ncid, time_varid, val_dim_nc(1:dim_length_nc(time_dim_nc_index),time_dim_nc_index), start=(/n_dims_start(2)/), count=(/n_dims_length(2)/)) )
        call check( nf90_put_var(ncid, time_varid, time_seconds_output_nc(1:nt), start=(/n_dims_start(2)/), count=(/n_dims_length(2)/)) )
        !!call check( nf90_put_var(ncid, station_varid, name_rec(:), start = (/1,1/), count=(/n_dims(1),n_char/)) )

        !Write station index and name
        call check( nf90_put_var(ncid, station_varid, id_rec, start = (/n_dims_start(1),1/), count=(/n_dims_length(1),n_char/)) )
        call check( nf90_put_var(ncid, station_name_varid, name_rec, start = (/1,n_dims_start(1)/), count=(/n_char,n_dims_length(1)/)) )

        !!call check( nf90_put_var(ncid, station_name_varid, name_rec, start = (/1,1/), count=(/n_char,n_dims(1)/)) )
        !!call check( nf90_put_var(ncid, station_name_varid, name_rec, start = (/1/), count=(/n_dims(1)/)) )

        !Write the variable to file
        if (nf90_type.eq.NF90_byte) then
            call check( nf90_put_var(ncid, val_varid, int(val_rec(:,1:nt),kind=1), start = (/n_dims_start(1),n_dims_start(2)/), count=(/n_dims_length(1),n_dims_length(2)/)) )
        elseif (nf90_type.eq.NF90_short) then
            call check( nf90_put_var(ncid, val_varid, int(val_rec(:,1:nt),kind=2), start = (/n_dims_start(1),n_dims_start(2)/), count=(/n_dims_length(1),n_dims_length(2)/)) )
        else
            call check( nf90_put_var(ncid, val_varid, val_rec(:,1:nt), start = (/n_dims_start(1),n_dims_start(2)/), count=(/n_dims_length(1),n_dims_length(2)/)) )
        endif

        !Write position data to the file
        call check( nf90_put_var(ncid, x_varid, x_rec, start = (/n_dims_start(1)/), count=(/n_dims_length(1)/)) )
        call check( nf90_put_var(ncid, y_varid, y_rec, start = (/n_dims_start(1)/), count=(/n_dims_length(1)/)) )
        call check( nf90_put_var(ncid, lon_varid, lon_rec, start = (/n_dims_start(1)/), count=(/n_dims_length(1)/)) )
        call check( nf90_put_var(ncid, lat_varid, lat_rec, start = (/n_dims_start(1)/), count=(/n_dims_length(1)/)) )
        call check( nf90_put_var(ncid, height_varid, height_rec, start = (/n_dims_start(1)/), count=(/n_dims_length(1)/)) )




        call check( nf90_close(ncid) )



    end subroutine uEMEP_save_netcdf_receptor_file

end module save_netcdf_file

