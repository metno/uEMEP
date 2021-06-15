!uEMEP_set_emission_factors.f90
    
    subroutine uEMEP_set_emission_factors
    
    use uEMEP_definitions
    
    implicit none
    
    write(unit_logfile,'(A)') ''
    write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Setting emission factors (uEMEP_set_emission_factors)'
	write(unit_logfile,'(A)') '================================================================'

    !Converts the emission units of the input data to a standard ug/s/subgrid
    emission_factor_conversion(nox_index,traffic_index,:)=emission_factor(nox_index,traffic_index,:)*(1.e-3)*(1.e+6)/(3600.*24.) ![veh*m/day]*(g/km/veh)*(km/m)*(ug/g)*(day/sec)=ug/sec
    emission_factor_conversion(nox_index,shipping_index,:)=emission_factor(nox_index,shipping_index,:)*(1.e+12)/(3600.*24.*365./12.) ![tonne/month]*(ug/ton)*(month/sec)=ug/sec. Only valid for the current AIS data which is a total. Should change
    emission_factor_conversion(nox_index,heating_index,:)=emission_factor(nox_index,heating_index,:)*(1.e+9)/(3600.*24.*365.) ![dwellings]*(kg/dwelling/year)*(ug/kg)*(year/sec)=ug/sec

    emission_factor_conversion(no2_index,traffic_index,:)=emission_factor(no2_index,traffic_index,:)*(1.e-3)*(1.e+6)/(3600.*24.) 
    emission_factor_conversion(no2_index,shipping_index,:)=emission_factor(no2_index,shipping_index,:)*(1.e+12)/(3600.*24.*365./12.)
    emission_factor_conversion(no2_index,heating_index,:)=emission_factor(no2_index,heating_index,:)*(1.e+9)/(3600.*24.*365.) ![dwellings]*(kg/dwelling/year)*(ug/kg)*(year/sec)=ug/sec

    emission_factor_conversion(pm25_index,traffic_index,:)=emission_factor(pm25_index,traffic_index,:)*(1.e-3)*(1.e+6)/(3600.*24.) ![veh*m/day]*(g/km/veh)*(km/m)*(ug/g)*(day/sec)=ug/sec
    emission_factor_conversion(pm25_index,shipping_index,:)=emission_factor(pm25_index,shipping_index,:)*(1.e+12)/(3600.*24.*365./12.) ![tonne/month]*(ug/kg)*(month/sec)=ug/sec
    emission_factor_conversion(pm25_index,heating_index,:)=emission_factor(pm25_index,heating_index,:)*(1.e+9)/(3600.*24.*365.) ![dwellings]*(kg/dwelling/year)*(ug/kg)*(year/sec)=ug/sec

    emission_factor_conversion(pm10_index,traffic_index,:)=emission_factor(pm10_index,traffic_index,:)*(1.e-3)*(1.e+6)/(3600.*24.) ![veh*m/day]*(g/km/veh)*(km/m)*(ug/g)*(day/sec)=ug/sec
    emission_factor_conversion(pm10_index,shipping_index,:)=emission_factor(pm10_index,shipping_index,:)*(1.e+12)/(3600.*24.*365./12.) ![tonne/month]*(ug/kg)*(month/sec)=ug/sec
    emission_factor_conversion(pm10_index,heating_index,:)=emission_factor(pm10_index,heating_index,:)*(1.e+9)/(3600.*24.*365.) ![dwellings]*(kg/dwelling/year)*(ug/kg)*(year/sec)=ug/sec

    emission_factor_conversion(pmex_index,traffic_index,:)=emission_factor(pmex_index,traffic_index,:)*(1.e-3)*(1.e+6)/(3600.*24.) ![veh*m/day]*(g/km/veh)*(km/m)*(ug/g)*(day/sec)=ug/sec
    
    emission_factor_conversion(nh3_index,agriculture_index,:)=emission_factor(nh3_index,agriculture_index,:)*(1.e+9)/(3600.*24.*365.)   ![kg/yr]*(ug/kg)*(yr/sec)=ug/sec
    !Test for monthly data
    !emission_factor_conversion(nh3_index,agriculture_index,:)=emission_factor_conversion(nh3_index,agriculture_index,:)*1.06
    
    if (read_shipping_from_netcdf_flag) then
        emission_factor_conversion(:,shipping_index,:)=(1.e+6)/(3600.) ![g/hr]*(ug/sec)=ug/sec. 
    elseif (read_weekly_shipping_data_flag.or.read_monthly_and_daily_shipping_data_flag) then
        !Convert from g/hour to ug/s
        emission_factor_conversion(:,shipping_index,:)=(1.e+6)/(3600.) ![g/hr]*(ug/sec)=ug/sec. 
    endif

    if (use_RWC_emission_data) then
        !Convert from g/h/subgrid to ug/s/subgrid 
        emission_factor_conversion(:,heating_index,:)=1.e6/3600.
        if (annual_calculations) then
            !Convert from g/yr/subgrid to ug/s/subgrid. Temporal profile will be set to 1.
            !Should also be in shipping. No because the original data is read in as g/hr
            emission_factor_conversion(:,heating_index,:)=1.e6/3600./24./365.
        endif
        
    endif
    
    !Converts tonne per year to ug/s/subgrid
    emission_factor_conversion(:,industry_index,:)=1.e12/(3600.*24*365.)
    
    
    end subroutine uEMEP_set_emission_factors

!uEMEP_convert_proxy_to_emissions
    
    subroutine uEMEP_convert_proxy_to_emissions
    
    use uEMEP_definitions
    
    implicit none
    
    integer i_source,i_subsource
    integer tt
    real sum_emission_subgrid
    
    integer i_pollutant
    
    !Do not calculate emissions if EMEP emissions are to be used
    if (local_subgrid_method_flag.eq.3) return
    
    write(unit_logfile,'(A)') ''
    write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Converting proxy data to emissions (uEMEP_convert_proxy_to_emissions)'
	write(unit_logfile,'(A)') '================================================================'

    !Set all emissions to the same constant emission value with emissions in ug/sec for all sources
    do i_source=1,n_source_index
    if (calculate_source(i_source)) then
        i_subsource=1
        do i_pollutant=1,n_pollutant_loop
            !Do not calculate for traffic if use_NORTRIP_emission_data=.true. and  use_NORTRIP_emission_pollutant=.false.). This is done in uEMEP_grid_roads
            if (i_source.ne.traffic_index &
                .or.(i_source.eq.traffic_index.and..not.use_NORTRIP_emission_data) &
                .or.(i_source.eq.traffic_index.and.use_NORTRIP_emission_data.and..not.use_NORTRIP_emission_pollutant(pollutant_loop_index(i_pollutant)))) then               
                
                do tt=1,subgrid_dim(t_dim_index)
                    emission_subgrid(:,:,tt,i_source,i_pollutant)=proxy_emission_subgrid(:,:,i_source,i_pollutant) &
                        *emission_factor_conversion(pollutant_loop_index(i_pollutant),i_source,i_subsource) &
                        *emission_time_profile_subgrid(:,:,tt,i_source,i_pollutant)
                    
                    !write(*,*) sum(proxy_emission_subgrid(:,:,i_source,i_subsource)),emission_factor_conversion(compound_index,i_source,i_subsource),sum(emission_time_profile_subgrid(:,:,tt,i_source,i_subsource))
                enddo
            endif
            
            sum_emission_subgrid=sum(emission_subgrid(:,:,:,i_source,i_pollutant))
            write(unit_logfile,'(A,es12.2)') 'Sum of emissions for '//trim(source_file_str(i_source))//' '//trim(pollutant_file_str(pollutant_loop_index(i_pollutant)))//' over period (kg)=',sum_emission_subgrid*3600./1.e9
         enddo
    endif
    enddo
    
    !Check NOx shipping
    !i_source=shipping_index
    !i_pollutant=pollutant_loop_back_index(nox_nc_index)
    !do tt=1,subgrid_dim(t_dim_index)
    !    write(unit_logfile,'(A,es12.2)') 'Sum of emissions per hour '//trim(source_file_str(i_source))//' '//trim(pollutant_file_str(pollutant_loop_index(i_pollutant)))//' over period (kg)=',sum(emission_subgrid(:,:,tt,i_source,i_pollutant))*3600./1.e9
    !enddo
    
    end subroutine uEMEP_convert_proxy_to_emissions

    subroutine uEMEP_nox_emission_temperature
    !Routine for doing the chemistry calculations in uEMEP
    
    use uEMEP_definitions

    implicit none
    
    integer t_start,t_end
    integer i,j,t
    integer i_cross,j_cross
    real :: ref_temperature1=-10.
    real :: ref_temperature2=10.
    real :: ref_scaling1=3.
    real :: ref_scaling2=1.
    real :: emission_scaling
    real :: temperature
    real :: emission_scaling_average,temperature_average
    integer :: emission_scaling_average_count
    
    real :: weighting_loop_sum
    real :: temperature_loop
    real :: weighting_loop
    real :: emission_scaling_loop
    real :: sigma_loop
    integer :: k
    
    !Do not correct if no traffic is to be calculated
    if (.not.calculate_source(traffic_index)) then
        return
    endif
    
    write(unit_logfile,'(A)') ''
    write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Using temperature dependency of NOx traffic emissions (uEMEP_nox_emission_temperature)'
	write(unit_logfile,'(A)') '================================================================'

	if (annual_calculations) write(unit_logfile,'(A)') 'Using temperature distribution to calculate the emissions scaling.'

    t_start=1
    t_end=subgrid_dim(t_dim_index)
    ref_temperature1=traffic_nox_emission_temperature_ref_temperature(1)
    ref_temperature2=traffic_nox_emission_temperature_ref_temperature(2)
    ref_scaling1=traffic_nox_emission_temperature_ref_scaling(1)
    ref_scaling2=traffic_nox_emission_temperature_ref_scaling(2)

    !Use the meteo model input from EMEP or alternative directly
    !This is so it will work for saving emissions for EMEP as well since the meteo and cross reference subgrids are not calculated
    emission_scaling_average=0.
    emission_scaling_average_count=0
    temperature_average=0
    
    do t=t_start,t_end
    do j=1,emission_subgrid_dim(y_dim_index,traffic_index)
    do i=1,emission_subgrid_dim(x_dim_index,traffic_index)

        if (save_emissions_for_EMEP(allsource_index)) then
            if (allocated(meteo_var1d_nc)) then
                !Need to cross reference the meteo grid to the emission grid as this is not done normally
                i_cross=1+floor((x_emission_subgrid(i,j,traffic_index)-meteo_var1d_nc(1,lon_nc_index))/meteo_dgrid_nc(lon_nc_index)+0.5)
                j_cross=1+floor((y_emission_subgrid(i,j,traffic_index)-meteo_var1d_nc(1,lat_nc_index))/meteo_dgrid_nc(lat_nc_index)+0.5)     
                !Because the meteo grid can be smaller than the EMEP grid then need to limit it
                i_cross=min(max(1,i_cross),dim_length_meteo_nc(x_dim_nc_index))
                j_cross=min(max(1,j_cross),dim_length_meteo_nc(y_dim_nc_index))
            else
                write(unit_logfile,'(a)') 'No meteo data available for traffic temperature correction. Stopping.'
                stop
            endif
        else
            !Use the EMEP meteorology
            i_cross=crossreference_emission_to_emep_subgrid(i,j,x_dim_index,traffic_index)
            j_cross=crossreference_emission_to_emep_subgrid(i,j,y_dim_index,traffic_index)
            i_cross=min(max(1,i_cross),dim_length_nc(x_dim_nc_index))
            j_cross=min(max(1,j_cross),dim_length_nc(y_dim_nc_index))
        endif
        
        !i_integral=crossreference_emission_to_integral_subgrid(i,j,x_dim_index,traffic_index)
        !j_integral=crossreference_emission_to_integral_subgrid(i,j,y_dim_index,traffic_index)
 
        if (use_alternative_meteorology_flag) then
            temperature=meteo_var3d_nc(i_cross,j_cross,t,t2m_nc_index)-273.14
        else
            temperature=var3d_nc(i_cross,j_cross,t,t2m_nc_index,allsource_index,meteo_p_loop_index)-273.14
        endif
        !temperature=meteo_subgrid(i_integral,j_integral,t,t2m_subgrid_index)-273.14
        
        !If annual calculation then use a Gaussian temperature distribution based on Norwegian conditions with a sigma of 7 degrees
        !This is hardcoded awaiting the proper implementation in NORTRIP
        if (annual_calculations) then
            sigma_loop=7.
            weighting_loop_sum=0
            emission_scaling=0
            do k=-20,20
                temperature_loop=temperature+k
                weighting_loop=exp(-(temperature_loop-temperature)*(temperature_loop-temperature)/2./(sigma_loop*sigma_loop))
                weighting_loop_sum=weighting_loop_sum+weighting_loop
                emission_scaling_loop=ref_scaling1+(ref_scaling2-ref_scaling1)*(temperature_loop-ref_temperature1)/(ref_temperature2-ref_temperature1)
                emission_scaling_loop=max(min(emission_scaling_loop,ref_scaling1),ref_scaling2)
                emission_scaling=emission_scaling+emission_scaling_loop*weighting_loop
                !write(*,*) k,temperature_loop,weighting_loop,emission_scaling_loop
            enddo
            emission_scaling=emission_scaling/weighting_loop_sum
        else
            emission_scaling=ref_scaling1+(ref_scaling2-ref_scaling1)*(temperature-ref_temperature1)/(ref_temperature2-ref_temperature1)
            emission_scaling=max(min(emission_scaling,ref_scaling1),ref_scaling2)
        endif
            
        emission_scaling_average=emission_scaling_average+emission_scaling
        temperature_average=temperature_average+temperature
        emission_scaling_average_count=emission_scaling_average_count+1
        
        !subgrid(i,j,t,local_subgrid_index,traffic_index,pollutant_loop_back_index(nox_nc_index))=emission_scaling*subgrid(i,j,t,local_subgrid_index,traffic_index,pollutant_loop_back_index(nox_nc_index))
        emission_subgrid(i,j,t,traffic_index,pollutant_loop_back_index(nox_nc_index))=emission_scaling*emission_subgrid(i,j,t,traffic_index,pollutant_loop_back_index(nox_nc_index))
        !write(*,'(3i,2f12.2)') i,j,t,temperature,emission_scaling
    enddo
    enddo
    enddo
    

    write(unit_logfile,'(a,2f12.2)') 'Average emissions scaling factor and temperature for traffic NOx: ',emission_scaling_average/emission_scaling_average_count,temperature_average/emission_scaling_average_count
    write(unit_logfile,'(a,4f12.2)') 'Parameters ref_temperature1,ref_temperature2,ref_scaling1,ref_scaling2: ',ref_temperature1,ref_temperature2,ref_scaling1,ref_scaling2

    
    end subroutine uEMEP_nox_emission_temperature
