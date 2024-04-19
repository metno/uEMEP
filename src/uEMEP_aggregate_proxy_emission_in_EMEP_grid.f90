module aggregate_proxy_emission_in_emep_grid
    ! NOTE This module is not currently in use. Should it be deleted?

    use uEMEP_definitions, only: x_dim_nc_index, dim_length_nc, y_dim_nc_index, &
        unit_conversion, n_pollutant_loop, pollutant_loop_index, n_source_index, &
        calculate_source, make_EMEP_grid_emission_data, unit_logfile, source_file_str, &
        pollutant_file_str, agriculture_index, agriculture_nc_index, traffic_index, &
        shipping_index, shipping_nc_index, traffic_nc_index, agriculture_nc_index, &
        heating_index, heating_nc_index, y_dim_index, x_dim_index, emission_subgrid_dim, &
        crossreference_emission_to_emep_subgrid, proxy_emission_subgrid, &
        compound_index, emission_factor_conversion, hourly_calculations, &
        emis_nc_index, allsource_index, var3d_nc, time_dim_nc_index, emission_subgrid_delta, &
        annual_calculations, var1d_nc, replace_EMEP_local_with_subgrid_local, &
        conc_nc_index, frac_nc_index, emission_file_index, pathname_grid, filename_grid, &
        file_tag, lat_nc_index, dgrid_nc
    use mod_read_esri_ascii_file, only: write_esri_ascii_file

    implicit none
    private

contains

!uEMEP_aggregate_proxy_emission_in_EMEP_grid.f90
    !This routine takes subgrid emissions and aggregates them in the EMEP grid
    !This is used for cross checking emissions
    
    subroutine uEMEP_aggregate_proxy_emission_in_EMEP_grid
    
    implicit none
    
    integer i,j
    integer i_source
    real, allocatable :: EMEP_aggregated_subgid_emission(:,:)
    real, allocatable :: EMEP_aggregated_emission(:,:)
    integer, allocatable :: EMEP_aggregated_subgid_emission_count(:,:)
    real, allocatable :: lon_array(:,:),lat_array(:,:)
    integer iii,jjj   
    character(256) temp_name
    real var3d_nc_local_temp
    integer t
    integer i_nc_source
    integer :: i_subsource=1
    
    integer i_pollutant,p_loop
       
    allocate (EMEP_aggregated_subgid_emission_count(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index)))
    allocate (EMEP_aggregated_subgid_emission(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index)))
    allocate (EMEP_aggregated_emission(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index)))
    allocate (lon_array(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index)))
    allocate (lat_array(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index)))
            
    EMEP_aggregated_subgid_emission=0.
    EMEP_aggregated_subgid_emission_count=0
    
    t=1
    unit_conversion=1.
    !unit_conversion(agriculture_index)=1.0e6   !Converts total emissions in kg to EMEPS mg/m^2
    !Units should be specified in the reading emission files
    
    !loop through the emission grid and aggregate data
    !Units are the same as EMEP mg/m2
    do i_pollutant=1,n_pollutant_loop
    p_loop=pollutant_loop_index(i_pollutant)
    
    do i_source=1,n_source_index
    if (calculate_source(i_source)) then
        if (make_EMEP_grid_emission_data(i_source)) then

            write(unit_logfile,'(A)') ''
            write(unit_logfile,'(A)') '================================================================'
	        write(unit_logfile,'(A)') 'Aggregating proxy emissions in EMEP grids for diagnostics (uEMEP_aggregate_proxy_emission_in_EMEP_grid)'
	        write(unit_logfile,'(A)') '================================================================'

            write(unit_logfile,'(A)') 'Saving proxy subgrid emissions in EMEP grid ('//trim(source_file_str(i_source))//'_'//trim(pollutant_file_str(p_loop))//')'
            if (i_source.eq.agriculture_index) then
                i_nc_source=agriculture_nc_index
            elseif (i_source.eq.traffic_index) then
                i_nc_source=traffic_nc_index
            elseif (i_source.eq.shipping_index) then
                i_nc_source=shipping_nc_index
            elseif (i_source.eq.heating_index) then
                i_nc_source=heating_nc_index
            else
                write(unit_logfile,'(A)') 'Undefined source in routine uEMEP_aggregate_proxy_emission_in_EMEP_grid. Stopping'
                stop
            endif
            

            EMEP_aggregated_subgid_emission=0.
            EMEP_aggregated_emission=0.
            EMEP_aggregated_subgid_emission_count=0
            do j=1,emission_subgrid_dim(y_dim_index,i_source)
            do i=1,emission_subgrid_dim(x_dim_index,i_source)
                  
                !ii=crossreference_emission_to_target_subgrid(i,j,x_dim_index,i_source)
                !jj=crossreference_emission_to_target_subgrid(i,j,y_dim_index,i_source)
                !iii=crossreference_target_to_emep_subgrid(ii,jj,x_dim_index)
                !jjj=crossreference_target_to_emep_subgrid(ii,jj,y_dim_index)
                iii=crossreference_emission_to_emep_subgrid(i,j,x_dim_index,i_source)
                jjj=crossreference_emission_to_emep_subgrid(i,j,y_dim_index,i_source)
            
                !EMEP_aggregated_subgid_emission(iii,jjj)=EMEP_aggregated_subgid_emission(iii,jjj)+sum(proxy_emission_subgrid(i,j,i_source,:)*emission_factor_conversion(compound_index,i_source,:)) &
                !    /(emission_subgrid_delta(y_dim_index,i_source)*emission_subgrid_delta(x_dim_index,i_source))*3600./1000. !Conversion from ug/sec/subgrid to mg/m2/hr (EMEP)
                !EMEP_aggregated_emission(iii,jjj)=EMEP_aggregated_emission(iii,jjj)+sum(var3d_nc(iii,jjj,:,emis_nc_index,i_nc_source))/dim_length_nc(time_dim_nc_index)
                EMEP_aggregated_subgid_emission_count(iii,jjj)=EMEP_aggregated_subgid_emission_count(iii,jjj)+1
                !Calculate as tonnes per year
                EMEP_aggregated_subgid_emission(iii,jjj)=EMEP_aggregated_subgid_emission(iii,jjj)+proxy_emission_subgrid(i,j,i_source,i_pollutant)*emission_factor_conversion(compound_index,i_source,i_subsource) &
                    *1.e-12*3600.*24.*365. !Conversion from ug/sec/subgrid to ton/year 
                if (hourly_calculations) then
                    EMEP_aggregated_emission(iii,jjj)=EMEP_aggregated_emission(iii,jjj)+sum(var3d_nc(iii,jjj,:,emis_nc_index,allsource_index,i_pollutant))/dim_length_nc(time_dim_nc_index) &
                    *(emission_subgrid_delta(y_dim_index,i_source)*emission_subgrid_delta(x_dim_index,i_source))*1.e-9*24.*365. !Conversion from mg/m2/hr to ton/year (EMEP)
                    !NOTE: Have put in allsource_index because on this occasion there were no sectors in the EMEP file for this. Beware for later comparisons!!!
                    !write(*,*) i_nc_source,i_pollutant,EMEP_aggregated_emission(iii,jjj),sum(var3d_nc(iii,jjj,:,emis_nc_index,i_nc_source,i_pollutant))
                endif
                if (annual_calculations) then
                    EMEP_aggregated_emission(iii,jjj)=EMEP_aggregated_emission(iii,jjj)+sum(var3d_nc(iii,jjj,:,emis_nc_index,i_nc_source,i_pollutant))/dim_length_nc(time_dim_nc_index) &
                    *(emission_subgrid_delta(y_dim_index,i_source)*emission_subgrid_delta(x_dim_index,i_source))*1.e-9 !Conversion from mg/m2/yr to ton/year (EMEP)
                endif
                
                
            enddo
            enddo
            
            do j=1,dim_length_nc(y_dim_nc_index)
            do i=1,dim_length_nc(x_dim_nc_index)
                if (EMEP_aggregated_subgid_emission_count(i,j).gt.0) then
                    EMEP_aggregated_subgid_emission(i,j)=EMEP_aggregated_subgid_emission(i,j) !/EMEP_aggregated_subgid_emission_count(i,j)
                    EMEP_aggregated_emission(i,j)=EMEP_aggregated_emission(i,j) !/EMEP_aggregated_subgid_emission_count(i,j)
                else
                    EMEP_aggregated_subgid_emission(i,j)=0.
                    EMEP_aggregated_emission(i,j)=0.
                endif
                lon_array(i,j)=var1d_nc(i,x_dim_nc_index)
                lat_array(i,j)=var1d_nc(j,y_dim_nc_index)
            enddo
            enddo
            
            !Check this. This is no longer valid and does not take into account neighbouring grids. Do not use
            if (replace_EMEP_local_with_subgrid_local(i_source)) then
                do j=1,dim_length_nc(y_dim_nc_index)
                do i=1,dim_length_nc(x_dim_nc_index)
                    if (EMEP_aggregated_subgid_emission(i,j).ge.0.and.var3d_nc(i,j,t,emis_nc_index,i_nc_source,i_pollutant).gt.0) then
                        var3d_nc_local_temp=var3d_nc(i,j,t,conc_nc_index,i_nc_source,i_pollutant)*(1.-var3d_nc(i,j,t,frac_nc_index,i_nc_source,i_pollutant)) !nonlocal contribution
                        var3d_nc(i,j,t,frac_nc_index,i_nc_source,i_pollutant)=EMEP_aggregated_subgid_emission(i,j)/var3d_nc(i,j,t,emis_nc_index,i_nc_source,i_pollutant)*var3d_nc(i,j,t,frac_nc_index,i_nc_source,i_pollutant) !New local fraction
                        var3d_nc(i,j,t,conc_nc_index,i_nc_source,i_pollutant)=var3d_nc(i,j,t,conc_nc_index,i_nc_source,i_pollutant)*var3d_nc(i,j,t,frac_nc_index,i_nc_source,i_pollutant)+var3d_nc_local_temp !new total contribution
                        var3d_nc(i,j,t,conc_nc_index,allsource_index,i_pollutant)=var3d_nc(i,j,t,conc_nc_index,i_nc_source,i_pollutant)
                    endif
                enddo
                enddo
            endif
                
    
            temp_name=trim(pathname_grid(emission_file_index(i_source)))//trim(filename_grid(emission_file_index(i_source)))//'_'//trim(pollutant_file_str(p_loop))//'_aggregated_proxy_EMEP_'//trim(file_tag)//'.asc'
            write(unit_logfile,'(a)')'Writing to: '//trim(temp_name)
            write(unit_logfile,'(a,f12.2)')'Total local emissions (ton/year): ',sum(EMEP_aggregated_subgid_emission)
            call write_esri_ascii_file(unit_logfile,temp_name,dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dgrid_nc(lat_nc_index), &
            EMEP_aggregated_subgid_emission(:,:),lon_array,lat_array)
            temp_name=trim(pathname_grid(emission_file_index(i_source)))//trim(filename_grid(emission_file_index(i_source)))//'_'//trim(pollutant_file_str(p_loop))//'_EMEP_'//trim(file_tag)//'.asc'
            write(unit_logfile,'(a)')'Writing to: '//trim(temp_name)
            write(unit_logfile,'(a,f12.2)')'Total EMEP emissions (ton/year): ',sum(EMEP_aggregated_emission)
            call write_esri_ascii_file(unit_logfile,temp_name,dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dgrid_nc(lat_nc_index), &
            EMEP_aggregated_emission(:,:),lon_array,lat_array)
            
            write(unit_logfile,'(a,f12.4)')'Ratio of local to EMEP emissions: ',sum(EMEP_aggregated_subgid_emission)/sum(EMEP_aggregated_emission)
            
        endif    
    endif
    enddo
    
    enddo
    
    if (allocated(EMEP_aggregated_subgid_emission_count)) deallocate (EMEP_aggregated_subgid_emission_count)
    if (allocated(EMEP_aggregated_subgid_emission)) deallocate (EMEP_aggregated_subgid_emission)
    if (allocated(EMEP_aggregated_emission)) deallocate (EMEP_aggregated_emission)
    if (allocated(lon_array)) deallocate (lon_array)
    if (allocated(lat_array)) deallocate (lat_array)
    
    end subroutine uEMEP_aggregate_proxy_emission_in_EMEP_grid

end module aggregate_proxy_emission_in_emep_grid

