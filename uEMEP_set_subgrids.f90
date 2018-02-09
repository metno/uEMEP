!uEMEP_set_subgrids.f90
    
    subroutine uEMEP_set_subgrids

    use uEMEP_definitions

    implicit none

    logical :: use_preset_grids=.false.
    
    
    !Set the subgrid step for integration. Reduces the time spent in performing the integral and the emission distribution
    !integral_subgrid_step=1
    !weighting step should be 1 until a routine is written that aggregates the data for the weighting
    !This is not used and can be removed
    weighting_step=1

    if (use_preset_grids) then
 
    write(unit_logfile,*)'Using preset grid parameters for: '//trim(file_tag)

    if (file_tag.eq.'Oslo_region_250m') then
        subgrid_delta(x_dim_index)=250.;subgrid_delta(y_dim_index)=250.
        subgrid_min(x_dim_index)=220000.;subgrid_min(y_dim_index)=6600000.
        subgrid_max(x_dim_index)=290000.;subgrid_max(y_dim_index)=6680000.
        filename_rl(1)='Road_data_ascii_ERFK_Norge.txt'
        if (annual_calculations) then
            pathname_EMEP='C:\uEMEP\EMEP_data\Hourly\'
            filename_EMEP(1)='Base_fullrun.nc'
            filename_EMEP(2)='uEMEP_full.nc'
        elseif (hourly_calculations) then
            pathname_EMEP='C:\uEMEP\EMEP_data\Hourly\'
            filename_EMEP(1)='Base_hour.nc'
            filename_EMEP(2)='uEMEP_hour.nc'
        endif
        
        calculate_source(traffic_index)=.true.
        calculate_source(shipping_index)=.true. 
        integral_subgrid_step=2
        h_emis(traffic_index,1)=2.
        h_emis(shipping_index,1)=25.
        sig_y_00(shipping_index,:)=5.
        sig_y_00(traffic_index,:)=1.
    endif

    if (file_tag.eq.'Oslo_region_100m') then
        subgrid_delta(x_dim_index)=100.;subgrid_delta(y_dim_index)=100.
        subgrid_min(x_dim_index)=220000.;subgrid_min(y_dim_index)=6600000.
        subgrid_max(x_dim_index)=290000.;subgrid_max(y_dim_index)=6680000.
        filename_rl(1)='Road_data_ascii_ERFK_Norge.txt'
        if (annual_calculations) then
            pathname_EMEP='C:\uEMEP\EMEP_data\Hourly\'
            filename_EMEP(1)='Base_fullrun.nc'
            filename_EMEP(2)='uEMEP_full.nc'
        elseif (hourly_calculations) then
            pathname_EMEP='C:\uEMEP\EMEP_data\Hourly\'
            filename_EMEP(1)='Base_hour.nc'
            filename_EMEP(2)='uEMEP_hour.nc'
        endif
        
        calculate_source(traffic_index)=.true.
        calculate_source(shipping_index)=.true. 
        integral_subgrid_step=2
        h_emis(traffic_index,1)=2.
        h_emis(shipping_index,1)=25.
        sig_y_00(shipping_index,:)=5.
        sig_y_00(traffic_index,:)=1.
     endif

    if (file_tag.eq.'Oslo_centre_25m') then
        subgrid_delta(x_dim_index)=25.;subgrid_delta(y_dim_index)=25.
        subgrid_min(x_dim_index)=250000.;subgrid_min(y_dim_index)=6645000.
        subgrid_max(x_dim_index)=270000.;subgrid_max(y_dim_index)=6660000.
        filename_rl(1)='Road_data_ascii_ERFK_large_Oslo.txt'
        integral_subgrid_step=10
        !weighting_step=4
    endif


    if (file_tag.eq.'Oslo_50m') then
        subgrid_delta(x_dim_index)=50.;subgrid_delta(y_dim_index)=50.
        subgrid_min(x_dim_index)=240000.;subgrid_min(y_dim_index)=6640000.
        subgrid_max(x_dim_index)=280000.;subgrid_max(y_dim_index)=6660000.
        filename_rl(1)='Road_data_ascii_ERFK_large_Oslo.txt'
        if (annual_calculations) then
            pathname_EMEP='C:\uEMEP\EMEP_data\Hourly\'
            filename_EMEP(1)='Base_fullrun.nc'
            filename_EMEP(2)='uEMEP_full.nc'
        elseif (hourly_calculations) then
            pathname_EMEP='C:\uEMEP\EMEP_data\Hourly\'
            filename_EMEP(1)='Base_hour.nc'
            filename_EMEP(2)='uEMEP_hour.nc'
        endif
        
        calculate_source(traffic_index)=.true.
        calculate_source(shipping_index)=.true. 
        integral_subgrid_step=10
        h_emis(traffic_index,1)=2.
        h_emis(shipping_index,1)=25.
        sig_y_00(shipping_index,:)=5.
        sig_y_00(traffic_index,:)=1.
    endif

    if (file_tag.eq.'Oslo_250m') then
        subgrid_delta(x_dim_index)=250.;subgrid_delta(y_dim_index)=250.
        subgrid_min(x_dim_index)=240000.;subgrid_min(y_dim_index)=6640000.
        subgrid_max(x_dim_index)=280000.;subgrid_max(y_dim_index)=6660000.
        filename_rl(1)='Road_data_ascii_ERFK_large_Oslo.txt'
        if (annual_calculations) then
            pathname_EMEP='C:\uEMEP\EMEP_data\Hourly\'
            filename_EMEP(1)='Base_fullrun.nc'
            filename_EMEP(2)='uEMEP_full.nc'
        elseif (hourly_calculations) then
            pathname_EMEP='C:\uEMEP\EMEP_data\Hourly\'
            filename_EMEP(1)='Base_hour.nc'
            filename_EMEP(2)='uEMEP_hour.nc'
        endif
        
        calculate_source(traffic_index)=.true.
        calculate_source(shipping_index)=.true. 
        calculate_source(heating_index)=.true. 
        integral_subgrid_step=2
        h_emis(traffic_index,1)=2.
        h_emis(shipping_index,1)=25.
        h_emis(heating_index,1)=25.
        sig_y_00(shipping_index,:)=5.
        sig_y_00(traffic_index,:)=1.
        sig_y_00(heating_index,:)=5.
   endif

    if (file_tag.eq.'Oslo_100m') then
        subgrid_delta(x_dim_index)=100.;subgrid_delta(y_dim_index)=100.
        subgrid_min(x_dim_index)=240000.;subgrid_min(y_dim_index)=6640000.
        subgrid_max(x_dim_index)=280000.;subgrid_max(y_dim_index)=6660000.
        filename_rl(1)='Road_data_ascii_ERFK_large_Oslo.txt'
        if (annual_calculations) then
            pathname_EMEP='C:\uEMEP\EMEP_data\Hourly\'
            filename_EMEP(1)='Base_fullrun.nc'
            filename_EMEP(2)='uEMEP_full.nc'
        elseif (hourly_calculations) then
            pathname_EMEP='C:\uEMEP\EMEP_data\Hourly\'
            filename_EMEP(1)='Base_hour.nc'
            filename_EMEP(2)='uEMEP_hour.nc'
        endif
        
        calculate_source(traffic_index)=.true.
        calculate_source(shipping_index)=.true. 
        integral_subgrid_step=5
        h_emis(traffic_index,1)=2.
        h_emis(shipping_index,1)=25.
        sig_y_00(shipping_index,:)=5.
        sig_y_00(traffic_index,:)=1.
    endif

     if (file_tag.eq.'Oslo_25m') then
        subgrid_delta(x_dim_index)=25.;subgrid_delta(y_dim_index)=25.
        subgrid_min(x_dim_index)=257000.;subgrid_min(y_dim_index)=6647000.
        subgrid_max(x_dim_index)=267000.;subgrid_max(y_dim_index)=6654000.
        filename_rl(1)='Road_data_ascii_ERFK_large_Oslo.txt'
        if (annual_calculations) then
            pathname_EMEP='C:\uEMEP\EMEP_data\Hourly\'
            filename_EMEP(1)='Base_fullrun.nc'
            filename_EMEP(2)='uEMEP_full.nc'
        elseif (hourly_calculations) then
            pathname_EMEP='C:\uEMEP\EMEP_data\Hourly\'
            filename_EMEP(1)='Base_hour.nc'
            filename_EMEP(2)='uEMEP_hour.nc'
        endif
        
        calculate_source(traffic_index)=.true.
        calculate_source(shipping_index)=.true. 
        integral_subgrid_step=20
        h_emis(traffic_index,1)=2.
        h_emis(shipping_index,1)=25.
        sig_y_00(shipping_index,:)=5.
        sig_y_00(traffic_index,:)=1.
     endif

     if (file_tag.eq.'Stavanger_250m') then
        subgrid_delta(x_dim_index)=250.;subgrid_delta(y_dim_index)=250.
        !subgrid_min(x_dim_index)=-60000.;subgrid_min(y_dim_index)=6520000.
        !subgrid_max(x_dim_index)=-5000.;subgrid_max(y_dim_index)=6600000.
        subgrid_min(x_dim_index)=-50000.;subgrid_min(y_dim_index)=6550000.
        subgrid_max(x_dim_index)=-25000.;subgrid_max(y_dim_index)=6600000.
        filename_rl(1)='Road_data_ascii_ERFK_Stavanger.txt'
        pathname_EMEP='C:\uEMEP\EMEP_data\Latest\'
        filename_EMEP='uEMEP_NOx_newmass_fullrun.nc'
        calculate_source(traffic_index)=.true.
        calculate_source(shipping_index)=.true. 
        calculate_source(agriculture_index)=.false.
        h_emis(traffic_index,1)=1.
        h_emis(shipping_index,1)=15.
        integral_subgrid_step=5
    endif

    if (file_tag.eq.'Stavanger_50m') then
        subgrid_delta(x_dim_index)=50.;subgrid_delta(y_dim_index)=50.
        subgrid_min(x_dim_index)=-50000.;subgrid_min(y_dim_index)=6550000.
        subgrid_max(x_dim_index)=-25000.;subgrid_max(y_dim_index)=6600000.
        filename_rl(1)='Road_data_ascii_ERFK_Stavanger.txt'
        pathname_EMEP='C:\uEMEP\EMEP_data\Latest\'
        filename_EMEP='uEMEP_NOx_newmass_fullrun.nc'
        calculate_source(traffic_index)=.true.
        calculate_source(shipping_index)=.true. 
        calculate_source(agriculture_index)=.false.
        h_emis(traffic_index,1)=1.
        h_emis(shipping_index,1)=25.
        integral_subgrid_step=1
        integral_subgrid_step=10
        !weighting_step=4        
   endif

    if (file_tag.eq.'Bergen_region_250m') then
        subgrid_delta(x_dim_index)=250.;subgrid_delta(y_dim_index)=250.
        subgrid_min(x_dim_index)=-60000.;subgrid_min(y_dim_index)=6700000.
        subgrid_max(x_dim_index)=-5000.;subgrid_max(y_dim_index)=6760000.
        filename_rl(1)='Road_data_ascii_ERFK_Bergen.txt'
        pathname_EMEP='C:\uEMEP\EMEP_data\Latest\'
        filename_EMEP='Base_fullrun.nc'
        calculate_source(traffic_index)=.true.
        calculate_source(shipping_index)=.true. 
        calculate_source(agriculture_index)=.false.
        h_emis(traffic_index,1)=1.
        h_emis(shipping_index,1)=15.
    endif

    if (file_tag.eq.'Bergen_50m') then
        subgrid_delta(x_dim_index)=50.;subgrid_delta(y_dim_index)=50.
        subgrid_min(x_dim_index)=-45000.;subgrid_min(y_dim_index)=6710000.
        subgrid_max(x_dim_index)=-20000.;subgrid_max(y_dim_index)=6750000.
        filename_rl(1)='Road_data_ascii_ERFK_Bergen.txt'
        if (annual_calculations) then
            pathname_EMEP='C:\uEMEP\EMEP_data\Hourly\'
            filename_EMEP(1)='Base_fullrun.nc'
            filename_EMEP(2)='uEMEP_full.nc'
        elseif (hourly_calculations) then
            pathname_EMEP='C:\uEMEP\EMEP_data\Hourly\'
            filename_EMEP(1)='Base_hour.nc'
            filename_EMEP(2)='uEMEP_hour.nc'
        endif      
        calculate_source(traffic_index)=.true.
        calculate_source(shipping_index)=.true. 
        integral_subgrid_step=10
        h_emis(traffic_index,1)=2.
        h_emis(shipping_index,1)=25.
        sig_y_00(shipping_index,:)=5.
        sig_y_00(traffic_index,:)=1.
    endif

    if (file_tag.eq.'Bergen_25m') then
        subgrid_delta(x_dim_index)=25.;subgrid_delta(y_dim_index)=25.
        subgrid_min(x_dim_index)=-45000.;subgrid_min(y_dim_index)=6710000.
        subgrid_max(x_dim_index)=-20000.;subgrid_max(y_dim_index)=6750000.
        filename_rl(1)='Road_data_ascii_ERFK_Bergen.txt'
         pathname_EMEP='C:\uEMEP\EMEP_data\Latest\'
        filename_EMEP='uEMEP_NOx_newmass_fullrun.nc'
        calculate_source(traffic_index)=.true.
        calculate_source(shipping_index)=.true. 
        calculate_source(agriculture_index)=.false.
        h_emis(traffic_index,1)=1.
        h_emis(shipping_index,1)=15.
        integral_subgrid_step=10
        sig_y_00(shipping_index,:)=5.
        sig_y_00(traffic_index,:)=1.
   endif

    if (file_tag.eq.'Southern_250m') then
        subgrid_delta(x_dim_index)=250.;subgrid_delta(y_dim_index)=250.
        !subgrid_min(x_dim_index)=180000.;subgrid_min(y_dim_index)=6520000.
        !subgrid_max(x_dim_index)=370000.;subgrid_max(y_dim_index)=6800000.
        !subgrid_min(x_dim_index)=-70000.;subgrid_min(y_dim_index)=6440000.
        !subgrid_max(x_dim_index)=350000.;subgrid_max(y_dim_index)=6830000.
        subgrid_min(x_dim_index)=-70000.;subgrid_min(y_dim_index)=6440000.
        subgrid_max(x_dim_index)=400000.;subgrid_max(y_dim_index)=7110000.
        filename_rl(1)='Road_data_ascii_ERFK_Norge.txt'
        if (annual_calculations) then
            pathname_EMEP='C:\uEMEP\EMEP_data\Hourly\'
            filename_EMEP(1)='Base_fullrun.nc'
            filename_EMEP(2)='uEMEP_full.nc'
        elseif (hourly_calculations) then
            pathname_EMEP='C:\uEMEP\EMEP_data\Hourly\'
            filename_EMEP(1)='Base_hour.nc'
            filename_EMEP(2)='uEMEP_hour.nc'
        endif
        
        calculate_source(traffic_index)=.true.
        calculate_source(shipping_index)=.true. 
        integral_subgrid_step=2
        h_emis(traffic_index,1)=2.
        h_emis(shipping_index,1)=25.
        sig_y_00(shipping_index,:)=5.
        sig_y_00(traffic_index,:)=1.
    endif
    if (file_tag.eq.'SouthEast_250m') then
        subgrid_delta(x_dim_index)=250.;subgrid_delta(y_dim_index)=250.
        subgrid_min(x_dim_index)=100000.;subgrid_min(y_dim_index)=6500000.
        subgrid_max(x_dim_index)=350000.;subgrid_max(y_dim_index)=6830000.
        filename_rl(1)='Road_data_ascii_ERFK_Norge.txt'
        if (annual_calculations) then
            pathname_EMEP='C:\uEMEP\EMEP_data\Hourly\'
            filename_EMEP(1)='Base_fullrun.nc'
            filename_EMEP(2)='uEMEP_full.nc'
        elseif (hourly_calculations) then
            pathname_EMEP='C:\uEMEP\EMEP_data\Hourly\'
            filename_EMEP(1)='Base_hour.nc'
            filename_EMEP(2)='uEMEP_hour.nc'
        endif
        
        calculate_source(traffic_index)=.true.
        calculate_source(shipping_index)=.true. 
        integral_subgrid_step=2
        h_emis(traffic_index,1)=2.
        h_emis(shipping_index,1)=25.
        sig_y_00(shipping_index,:)=5.
        sig_y_00(traffic_index,:)=1.
    endif
    if (file_tag.eq.'Narvik_50m') then
        subgrid_delta(x_dim_index)=50.;subgrid_delta(y_dim_index)=50.
        subgrid_min(x_dim_index)=600000.-20000.;subgrid_min(y_dim_index)=7594000.-20000.
        subgrid_max(x_dim_index)=600000.+20000.;subgrid_max(y_dim_index)=7594000.+20000.
        filename_rl(1)='Road_data_ascii_ERFK_Norge.txt'
        pathname_EMEP='C:\uEMEP\EMEP_data\Latest\'
        filename_EMEP='uEMEP_NOx_newmass_fullrun.nc'
        calculate_source(traffic_index)=.true.
        calculate_source(shipping_index)=.true. 
        calculate_source(agriculture_index)=.false.
        h_emis(traffic_index,1)=1.
        h_emis(shipping_index,1)=15.
        integral_subgrid_step=10
        sig_y_00(shipping_index,:)=5.
        sig_y_00(traffic_index,:)=1.
    endif
    if (file_tag.eq.'Nederland_500m') then
        hourly_calculations=.false.
        annual_calculations=.true.
        compound_index=nh3_nc_index
        compound_frac_index=nh3_nc_index

        subgrid_delta(x_dim_index)=500.;subgrid_delta(y_dim_index)=500.
        subgrid_min(x_dim_index)=13000.;subgrid_min(y_dim_index)=304000.
        subgrid_max(x_dim_index)=280000.;subgrid_max(y_dim_index)=618000.
        pathname_EMEP='C:\uEMEP\EMEP_data\nh3\'
        !filename_EMEP='uEMEP_nh3_RIVM_EMIS_fullrun.nc'
        filename_EMEP(1)='NH3_fullrun.nc'
        filename_EMEP(2)='uEMEP_full.nc'
        projection_type=RDM_projection_index
        calculate_source(traffic_index)=.false.
        calculate_source(shipping_index)=.false. 
        calculate_source(agriculture_index)=.true.
    endif
    if (file_tag.eq.'Nederland_250m') then
        subgrid_delta(x_dim_index)=250.;subgrid_delta(y_dim_index)=250.
        subgrid_min(x_dim_index)=13000.;subgrid_min(y_dim_index)=304000.
        subgrid_max(x_dim_index)=280000.;subgrid_max(y_dim_index)=618000.
        pathname_EMEP='C:\uEMEP\EMEP_data\nh3\'
        filename_EMEP='uEMEP_nh3_RIVM_EMIS_fullrun.nc'
        projection_type=RDM_projection_index
        calculate_source(traffic_index)=.false.
        calculate_source(shipping_index)=.false. 
        calculate_source(agriculture_index)=.true.
    endif
    if (file_tag.eq.'Nederland_100m') then
        subgrid_delta(x_dim_index)=100.;subgrid_delta(y_dim_index)=100.
        !subgrid_min(x_dim_index)=150000.;subgrid_min(y_dim_index)=430000.
        !subgrid_max(x_dim_index)=195000.;subgrid_max(y_dim_index)=485000.
        subgrid_min(x_dim_index)=13000.;subgrid_min(y_dim_index)=304000.
        subgrid_max(x_dim_index)=280000.;subgrid_max(y_dim_index)=618000.
        pathname_EMEP='C:\uEMEP\EMEP_data\nh3\'
        filename_EMEP='uEMEP_nh3_RIVM_EMIS_fullrun.nc'
        projection_type=RDM_projection_index
        calculate_source(traffic_index)=.false.
        calculate_source(shipping_index)=.false. 
        calculate_source(agriculture_index)=.true.
    endif
    if (file_tag.eq.'Nederland_50km_500m') then
        subgrid_delta(x_dim_index)=500.;subgrid_delta(y_dim_index)=500.
        subgrid_min(x_dim_index)=13000.;subgrid_min(y_dim_index)=304000.
        subgrid_max(x_dim_index)=280000.;subgrid_max(y_dim_index)=618000.
        pathname_EMEP='C:\uEMEP\EMEP_data\nh3\'
        filename_EMEP='uEMEP05_nh3_fullrun.nc'
        projection_type=RDM_projection_index
    endif

    endif
    
    !Reset min and max with the buffer and calculate dimensions
    !subgrid_min(x_dim_index)=subgrid_min(x_dim_index)-buffer(x_dim_index);subgrid_min(y_dim_index)=subgrid_min(y_dim_index)-buffer(y_dim_index)
    !subgrid_max(x_dim_index)=subgrid_max(x_dim_index)+buffer(x_dim_index);subgrid_max(y_dim_index)=subgrid_max(y_dim_index)+buffer(y_dim_index)
    subgrid_dim(x_dim_index)=floor((subgrid_max(x_dim_index)-subgrid_min(x_dim_index))/subgrid_delta(x_dim_index))
    subgrid_dim(y_dim_index)=floor((subgrid_max(y_dim_index)-subgrid_min(y_dim_index))/subgrid_delta(y_dim_index))
 
    !write(*,*) subgrid_dim(x_dim_index),subgrid_max(x_dim_index),subgrid_min(x_dim_index),subgrid_delta(x_dim_index)
    !write(*,*) subgrid_dim(y_dim_index),subgrid_max(y_dim_index),subgrid_min(y_dim_index),subgrid_delta(y_dim_index)
    
    !Set all integral subgrids relative to the target subgrid
    integral_subgrid_delta=subgrid_delta*integral_subgrid_step
    integral_subgrid_min=subgrid_min
    integral_subgrid_max=subgrid_max
    integral_subgrid_dim(x_dim_index)=floor((subgrid_max(x_dim_index)-subgrid_min(x_dim_index))/integral_subgrid_delta(x_dim_index))+1
    integral_subgrid_dim(y_dim_index)=floor((subgrid_max(y_dim_index)-subgrid_min(y_dim_index))/integral_subgrid_delta(y_dim_index))+1
    integral_subgrid_dim(t_dim_index)=subgrid_dim(t_dim_index)
    !Set the integral subgrid dimensions so they cannot be larger than the target subgrid
    integral_subgrid_dim(x_dim_index)=min(integral_subgrid_dim(x_dim_index),subgrid_dim(x_dim_index))
    integral_subgrid_dim(y_dim_index)=min(integral_subgrid_dim(y_dim_index),subgrid_dim(y_dim_index))

    !Set all population subgrids relative to the target subgrid
    if (population_data_type.eq.population_index) then
        !When 250 m population data is used then set this as a limit
        population_subgrid_delta(x_dim_index)=max(subgrid_delta(x_dim_index),limit_population_delta)
        population_subgrid_delta(y_dim_index)=max(subgrid_delta(y_dim_index),limit_population_delta)
    else
        !Allow the point population data to have the same grid as the target grid
        population_subgrid_delta(x_dim_index)=subgrid_delta(x_dim_index)
        population_subgrid_delta(y_dim_index)=subgrid_delta(y_dim_index)
    endif        
    
    population_subgrid_min=subgrid_min
    population_subgrid_max=subgrid_max
    population_subgrid_dim(x_dim_index)=floor((population_subgrid_max(x_dim_index)-population_subgrid_min(x_dim_index))/population_subgrid_delta(x_dim_index))+1
    population_subgrid_dim(y_dim_index)=floor((population_subgrid_max(y_dim_index)-population_subgrid_min(y_dim_index))/population_subgrid_delta(y_dim_index))+1
    !Set the population subgrid dimensions so they cannot be larger than the target subgrid. Not certain why I do this.
    population_subgrid_dim(x_dim_index)=min(population_subgrid_dim(x_dim_index),subgrid_dim(x_dim_index))
    population_subgrid_dim(y_dim_index)=min(population_subgrid_dim(y_dim_index),subgrid_dim(y_dim_index))

    !Set all emission subgrids to be the same as the target subgrid
    emission_max_subgrid_dim=subgrid_dim
    do i=1,n_source_index
        emission_subgrid_delta(:,i)=subgrid_delta
        emission_subgrid_min(:,i)=subgrid_min
        emission_subgrid_max(:,i)=subgrid_max
        emission_subgrid_dim(:,i)=subgrid_dim
    enddo

    !Set shipping data to a minimum vale for all sources (Cannot be smaller than the target subgrid)
    emission_subgrid_delta(x_dim_index,shipping_index)=max(subgrid_delta(x_dim_index),limit_shipping_delta)
    emission_subgrid_delta(y_dim_index,shipping_index)=max(subgrid_delta(y_dim_index),limit_shipping_delta)
    emission_subgrid_delta(x_dim_index,heating_index)=max(subgrid_delta(x_dim_index),limit_heating_delta)
    emission_subgrid_delta(y_dim_index,heating_index)=max(subgrid_delta(y_dim_index),limit_heating_delta)
    
    !Set all the emission subgrid dimmensions after changes
    do i=1,n_source_index
        emission_subgrid_dim(x_dim_index,i)=floor((emission_subgrid_max(x_dim_index,i)-emission_subgrid_min(x_dim_index,i))/emission_subgrid_delta(x_dim_index,i))
        emission_subgrid_dim(y_dim_index,i)=floor((emission_subgrid_max(y_dim_index,i)-emission_subgrid_min(y_dim_index,i))/emission_subgrid_delta(y_dim_index,i))
        write(unit_logfile,'(A,I6,A5,2I6)') 'Emission grid dimensions for source ',i,': ',emission_subgrid_dim(1:2,i)
    enddo
    
    write(unit_logfile,'(A,2I6)') 'Target grid dimensions: ',subgrid_dim(1:2)
    write(unit_logfile,'(A,2I6)') 'Integral grid dimensions: ',integral_subgrid_dim(1:2) 
    write(unit_logfile,'(A,2f10.1)') 'Target subgrid grid sizes: ',subgrid_delta
    write(unit_logfile,'(A,I6,2f10.1)') 'Integral subgrid step and grid sizes: ',integral_subgrid_step,integral_subgrid_delta
    
    end subroutine uEMEP_set_subgrids