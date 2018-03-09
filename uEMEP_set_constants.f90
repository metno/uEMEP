    subroutine uEMEP_set_constants
    
    use uEMEP_definitions
    
    implicit none
    
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
        var_name_nc(conc_nc_index,o3_nc_index,:)='o3'
        var_name_nc(conc_nc_index,no2_nc_index,:)='no2'
        var_name_nc(conc_nc_index,nox_nc_index,:)='nox'
        var_name_nc(conc_nc_index,nh3_nc_index,:)='nh3'
        var_name_nc(conc_nc_index,pm25_nc_index,:)='pm25'
        var_name_nc(conc_nc_index,pm10_nc_index,:)='pm10'
        
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
        var_name_nc(emis_nc_index,nh3_nc_index,:)='Emis_mgm2_nh3'
        var_name_nc(emis_nc_index,nox_nc_index,:)='Emis_mgm2_nox'
        var_name_nc(emis_nc_index,pm10_nc_index,:)='Emis_mgm2_pm10'
        var_name_nc(emis_nc_index,pm25_nc_index,:)='Emis_mgm2_pm25'
        
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
        var_name_nc(ugrid_nc_index,:,allsource_index)='u_wind'
        var_name_nc(vgrid_nc_index,:,allsource_index)='v_wind'
        var_name_nc(FFgrid_nc_index,:,allsource_index)='wind_speed'
        var_name_nc(FF10_nc_index,:,allsource_index)='ws10m'
        var_name_nc(inv_FFgrid_nc_index,:,allsource_index)='inv_wind_speed'
        var_name_nc(inv_FF10_nc_index,:,allsource_index)='inv_ws10m'
        var_name_nc(hmix_nc_index,:,allsource_index)='HMIX'
        var_name_nc(kz_nc_index,:,allsource_index)='Kz'
        var_name_nc(ustar_nc_index,:,allsource_index)='met2d_ustar_nwp'
        var_name_nc(logz0_nc_index,:,allsource_index)='av_logz0'
        var_name_nc(invL_nc_index,:,allsource_index)='inv_L'
        var_name_nc(J_nc_index,:,allsource_index)='J(NO2)'
        var_name_nc(ZTOP_nc_index,:,allsource_index)='Z_TOP'
        surface_level_nc=7 !Will be reset as length of 'lev' dimension

        !Additional compounds for chemistry and totals
        comp_name_nc(o3_nc_index)='D3_ug_O3'
        comp_name_nc(no2_nc_index)='D3_ug_NO2'
        comp_name_nc(nox_nc_index)='D3_ugN_NOX'
        comp_name_nc(nh3_nc_index)='D3_ug_NH3'
        comp_name_nc(pm25_nc_index)='pm25'
        comp_name_nc(pm10_nc_index)='pm10'

        comp_scale_nc(:)=1.
        comp_scale_nc(nox_nc_index)=(14+2.*16.)/14. !Value read in is in ugN
        
    !Allocate the indexes for specifying compound and sources
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
    

       
    end subroutine uEMEP_set_constants