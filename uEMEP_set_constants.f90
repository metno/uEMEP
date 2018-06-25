    subroutine uEMEP_set_constants
    
    use uEMEP_definitions
    
    implicit none
    
    integer i,j,k
    integer count
    character(1) str_temp
    
    
    !Initialise some array flags
    read_existing_grid_data=.false.
	replace_EMEP_local_with_subgrid_local=.false.
	calculate_source=.false.
    combine_emission_subsources_during_dispersion=.true.
    
    
    dim_name_nc(x_dim_nc_index)='lon'
    dim_name_nc(y_dim_nc_index)='lat'
    dim_name_nc(z_dim_nc_index)='lev'
    dim_name_nc(time_dim_nc_index)='time'
    dim_name_nc(xdist_dim_nc_index)='x_dist'
    dim_name_nc(ydist_dim_nc_index)='y_dist'


        !Concentrations
        var_name_nc=''
        var_name_nc(conc_nc_index,o3_nc_index,allsource_nc_index)='o3'
        var_name_nc(conc_nc_index,no2_nc_index,allsource_nc_index)='no2'
        var_name_nc(conc_nc_index,nox_nc_index,allsource_nc_index)='nox'
        var_name_nc(conc_nc_index,nh3_nc_index,allsource_nc_index)='nh3'
        var_name_nc(conc_nc_index,pm25_nc_index,allsource_nc_index)='pm25'
        var_name_nc(conc_nc_index,pm10_nc_index,allsource_nc_index)='pm10'
        
        !Local fractions
        var_name_nc(frac_nc_index,nox_nc_index,traffic_nc_index)='nox_sec07_local_fraction'
        var_name_nc(frac_nc_index,pm10_nc_index,traffic_nc_index)='pm10_sec07_local_fraction'
        var_name_nc(frac_nc_index,pm25_nc_index,traffic_nc_index)='pm25_sec07_local_fraction'
 
        var_name_nc(frac_nc_index,nox_nc_index,shipping_nc_index)='nox_sec08_local_fraction'
        var_name_nc(frac_nc_index,pm25_nc_index,shipping_nc_index)='pm25_sec08_local_fraction'
        var_name_nc(frac_nc_index,pm10_nc_index,shipping_nc_index)='pm10_sec08_local_fraction'

        var_name_nc(frac_nc_index,nh3_nc_index,agriculture_nc_index)='nh3_sec10_local_fraction'

        var_name_nc(frac_nc_index,nox_nc_index,heating_nc_index)='nox_sec02_local_fraction'
        var_name_nc(frac_nc_index,pm25_nc_index,heating_nc_index)='pm25_sec02_local_fraction'
        var_name_nc(frac_nc_index,pm10_nc_index,heating_nc_index)='pm10_sec02_local_fraction'

        !Total emissions
        var_name_nc(emis_nc_index,nh3_nc_index,allsource_nc_index)='Emis_mgm2_nh3'
        var_name_nc(emis_nc_index,nox_nc_index,allsource_nc_index)='Emis_mgm2_nox'
        var_name_nc(emis_nc_index,pm10_nc_index,allsource_nc_index)='Emis_mgm2_pm10'
        var_name_nc(emis_nc_index,pm25_nc_index,allsource_nc_index)='Emis_mgm2_pm25'
        
        !Sector emissions
        var_name_nc(emis_nc_index,nox_nc_index,traffic_nc_index)='Emis_mgm2_sec7nox'
        var_name_nc(emis_nc_index,pm25_nc_index,traffic_nc_index)='Emis_mgm2_sec7pm25'
        var_name_nc(emis_nc_index,pm10_nc_index,traffic_nc_index)='Emis_mgm2_sec7pm10'

        var_name_nc(emis_nc_index,nox_nc_index,shipping_nc_index)='Emis_mgm2_sec8nox'
        var_name_nc(emis_nc_index,pm25_nc_index,shipping_nc_index)='Emis_mgm2_sec8pm25'
        var_name_nc(emis_nc_index,pm10_nc_index,shipping_nc_index)='Emis_mgm2_sec8pm10'
        
        var_name_nc(emis_nc_index,nox_nc_index,heating_nc_index)='Emis_mgm2_sec2nox'
        var_name_nc(emis_nc_index,pm25_nc_index,heating_nc_index)='Emis_mgm2_sec2pm25'
        var_name_nc(emis_nc_index,pm10_nc_index,heating_nc_index)='Emis_mgm2_sec2pm10'
              
        !Meteorology
        var_name_nc(ugrid_nc_index,:,allsource_nc_index)='u_wind'
        var_name_nc(vgrid_nc_index,:,allsource_nc_index)='v_wind'
        var_name_nc(FFgrid_nc_index,:,allsource_nc_index)='wind_speed'
        var_name_nc(FF10_nc_index,:,allsource_nc_index)='ws10m'
        var_name_nc(inv_FFgrid_nc_index,:,allsource_nc_index)='inv_wind_speed'
        var_name_nc(inv_FF10_nc_index,:,allsource_nc_index)='inv_ws10m'
        var_name_nc(hmix_nc_index,:,allsource_nc_index)='HMIX'
        var_name_nc(kz_nc_index,:,allsource_nc_index)='Kz'
        var_name_nc(ustar_nc_index,:,allsource_nc_index)='met2d_ustar_nwp'
        var_name_nc(logz0_nc_index,:,allsource_nc_index)='av_logz0'
        var_name_nc(invL_nc_index,:,allsource_nc_index)='inv_L'
        var_name_nc(J_nc_index,:,allsource_nc_index)='J(NO2)'
        var_name_nc(ZTOP_nc_index,:,allsource_nc_index)='Z_TOP'
        var_name_nc(t2m_nc_index,:,allsource_nc_index)='met2d_t2m'
        surface_level_nc=7 !Will be reset as length of 'lev' dimension

        !Alternative meteorology.
        dim_name_meteo_nc(x_dim_nc_index)='x'
        dim_name_meteo_nc(y_dim_nc_index)='y'
        dim_name_meteo_nc(z_dim_nc_index)='height3'
        dim_name_meteo_nc(time_dim_nc_index)='time'

        var_name_meteo_nc=''
        var_name_meteo_nc(lon_nc_index)='lon'
        var_name_meteo_nc(lat_nc_index)='lat'
        var_name_meteo_nc(ugrid_nc_index)='u_wind' !!Doesn't exist but name is used
        var_name_meteo_nc(vgrid_nc_index)='v_wind' !!Doesn't exist but name is used
        var_name_meteo_nc(FFgrid_nc_index)='FF_wind' !Doesn't exist but name is used
        var_name_meteo_nc(FF10_nc_index)='wind_speed' 
        var_name_meteo_nc(inv_FFgrid_nc_index)='' 
        var_name_meteo_nc(inv_FF10_nc_index)='' 
        var_name_meteo_nc(hmix_nc_index)='atmosphere_boundary_layer_thickness'
        var_name_meteo_nc(kz_nc_index)=''
        var_name_meteo_nc(ustar_nc_index)='ustar' !Doesn't exist but name is used
        var_name_meteo_nc(logz0_nc_index)='Z0' !Needs to be converted to log(Z0)
        var_name_meteo_nc(invL_nc_index)='invL' !!Doesn't exist but name is used. Will be calculated
        var_name_meteo_nc(J_nc_index)=''
        var_name_meteo_nc(ZTOP_nc_index)=''
        var_name_meteo_nc(t2m_nc_index)='air_temperature_2m'
        !Additional     parameter (u10_nc_subgrid_index=22,v10_nc_subgrid_index=23,uw_nc_subgrid_index=24,vw_nc_subgrid_index=25,Hflux_nc_subgrid_index=26)
        var_name_meteo_nc(u10_nc_index)='x_wind_10m' !10 m wind not grid. Replaces ugrid. Used for direction
        var_name_meteo_nc(v10_nc_index)='y_wind_10m' !10 m wind not grid. Replaces vgrid. Used for direction
        var_name_meteo_nc(uw_nc_index)='downward_eastward_momentum_flux_in_air' !Will be used to determine ustar
        var_name_meteo_nc(vw_nc_index)='downward_northward_momentum_flux_in_air' !Will be used to determine ustar
        var_name_meteo_nc(Hflux_nc_index)='integral_of_surface_downward_sensible_heat_flux_wrt_time'
        
       
        !Additional compounds for chemistry and totals
        comp_name_nc(o3_nc_index)='D3_ug_O3'
        comp_name_nc(no2_nc_index)='D3_ug_NO2'
        comp_name_nc(nox_nc_index)='D3_ugN_NOX'
        comp_name_nc(nh3_nc_index)='D3_ug_NH3'
        comp_name_nc(pm25_nc_index)='pm25'
        comp_name_nc(pm10_nc_index)='pm10'
        comp_name_nc(pm25_nc_index)='D3_ug_PMFINE'
        !comp_name_nc(pm25_nc_index)='SURF_ug_PM25_rh50'

        comp_scale_nc(:)=1.
        comp_scale_nc(nox_nc_index)=(14.+2.*16.)/14. !Value read in is in ugN
        
    !Allocate the indexes for specifying compound and sources together. Never used and not correct either!!!
    count=0
    do i=1,size(compound_source_index,1)
        do j=1,size(compound_source_index,2)
            count=count+1
            compound_source_index(i,j)=count
        enddo
    enddo
    
     
    !Allocate source strings for writing to files
    source_file_str(allsource_index)='allsource'
    source_file_str(traffic_index)='traffic'
    source_file_str(shipping_index)='shipping'
    source_file_str(agriculture_index)='agriculture'
    source_file_str(heating_index)='heating'
    source_file_str(industry_index)='industry'

    do i=1,n_possible_subsource
        write(str_temp,'(i1)') i
        subsource_str(i)='_'//trim(str_temp)
    enddo
    
    
    !Set filename indexes for grids generated by uEMEP
    do i=1,n_source_index
        j=j+1;proxy_emission_file_index(i)=j
        j=j+1;emission_file_index(i)=j
        j=j+1;proxy_file_index(i)=j
        j=j+1;proxy_integral_file_index(i)=j
        j=j+1;emep_subgrid_file_index(i)=j
        j=j+1;emep_subgrid_nonlocal_file_index(i)=j
        j=j+1;emep_subgrid_local_file_index(i)=j
        j=j+1;emep_subgrid_frac_file_index(i)=j
        j=j+1;subgrid_local_file_index(i)=j
        j=j+1;subgrid_total_file_index(i)=j
        j=j+1;use_subgrid_file_index(i)=j
        j=j+1;emep_emission_subgrid_file_index(i)=j
    enddo
    do i=1,n_population_index
        j=j+1;population_file_index(i)=j
    enddo
    j=j+1;subgrid_ugrid_file_index=j
    j=j+1;subgrid_vgrid_file_index=j
    j=j+1;subgrid_hmix_file_index=j
    j=j+1;subgrid_kz_file_index=j
    j=j+1;subgrid_logz0_file_index=j
    j=j+1;subgrid_invL_file_index=j     
    j=j+1;subgrid_FFgrid_file_index=j     
    j=j+1;subgrid_FF10_file_index=j
    j=j+1;subgrid_invFFgrid_file_index=j     
    j=j+1;subgrid_invFF10_file_index=j
    j=j+1;subgrid_ustar_file_index=j
    j=j+1;subgrid_J_file_index=j
    j=j+1;subgrid_meteo_file_index=j
    
    
    !Set initial values  for the dispersion parameters
    sig_y_00=0.
    sig_z_00=0.
    h_emis=0.
    z_rec=2.
    
    h_emis(traffic_index,:)=2.
    h_emis(shipping_index,:)=25.
    h_emis(heating_index,:)=25.
    h_emis(agriculture_index,:)=1.
    h_emis(industry_index,:)=25.
    sig_y_00(shipping_index,:)=5.
    sig_y_00(traffic_index,:)=1.
    sig_y_00(heating_index,:)=5.
    sig_y_00(agriculture_index,:)=5.
    sig_y_00(industry_index,:)=5.
    sig_z_00(shipping_index,:)=5.
    sig_z_00(traffic_index,:)=1.
    sig_z_00(heating_index,:)=10.
    sig_z_00(agriculture_index,:)=10.
    sig_z_00(industry_index,:)=10.

    !Set all emission factors to unity
    emission_factor=1.
    
    !Preset all initial emission factors 
    emission_factor(nox_index,traffic_index,:)=0.4 !(g/km/veh)
    emission_factor(nox_index,shipping_index,:)=1. !Shipping data is in emissions [tonne/month]
    emission_factor(nox_index,heating_index,:)=3./15. !(kg/dwelling/year) Estimate only

    emission_factor(no2_index,traffic_index,:)=0.15*emission_factor(nox_index,traffic_index,:)
    emission_factor(no2_index,shipping_index,:)=0.10*emission_factor(nox_index,shipping_index,:)
    emission_factor(no2_index,heating_index,:)=0.1*emission_factor(nox_index,heating_index,:)!(kg/dwelling/year) Estimate only

    emission_factor(pm25_index,traffic_index,:)=0.01 !(g/km/veh)
    emission_factor(pm25_index,shipping_index,:)=1. !Shipping data is in emissions [tonne/month]
    emission_factor(pm25_index,heating_index,:)=3. !(kg/dwelling/year) SSB number is 6

    emission_factor(pm10_index,traffic_index,:)=0.01 !(g/km/veh)
    emission_factor(pm10_index,shipping_index,:)=1. !Shipping data is in emissions [tonne/month]
    emission_factor(pm10_index,heating_index,:)=3. !(kg/dwelling/year) SSB number is 6

    emission_factor(nh3_index,agriculture_index,:)=1. !Agriculture data is in emissions [kg/yr]

    ratio_truck_car_emission(nox_index)=4.86/.318 !From excel sheet for NOx
    ratio_truck_car_emission(no2_index)=4.86/.318 !Should perhaps be different
    ratio_truck_car_emission(pm25_index)=10.
    ratio_truck_car_emission(pm10_index)=10.
   
    end subroutine uEMEP_set_constants