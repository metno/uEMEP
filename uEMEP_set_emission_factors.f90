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
    
    if (use_RWC_emission_data) then
        !Convert from g/h/subgrid to ug/s/subgrid 
        emission_factor_conversion(:,heating_index,:)=1.e6/3600.
    endif
    
    
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
    
    
    end subroutine uEMEP_convert_proxy_to_emissions
