module uEMEP_set_sources

    use uEMEP_definitions
    use uemep_configuration

    implicit none

contains

    subroutine set_sources()

        integer :: ic

        ic = 0

        ic = ic + 1; allsource_index = ic; allsource_nc_index = allsource_index
        ic = ic + 1; traffic_index = ic; traffic_nc_index = traffic_index
        ic = ic + 1; shipping_index = ic; shipping_nc_index = shipping_index
        ic = ic + 1; heating_index = ic; heating_nc_index = heating_index
        ic = ic + 1; agriculture_index = ic; agriculture_nc_index = agriculture_index
        ic = ic + 1; industry_index = ic; industry_nc_index = industry_index
        ic = ic + 1; publicpower_index = ic; publicpower_nc_index = publicpower_index
        ic = ic + 1; fugitive_index = ic; fugitive_nc_index = fugitive_index
        ic = ic + 1; solvents_index = ic; solvents_nc_index = solvents_index
        ic = ic + 1; aviation_index = ic; aviation_nc_index = aviation_index
        ic = ic + 1; offroad_index = ic; offroad_nc_index = offroad_index
        if (include_source_waste) ic = ic + 1; waste_index = ic; waste_nc_index = waste_index
        if (include_source_livestock) ic = ic + 1; livestock_index = ic; livestock_nc_index = livestock_index
        if (include_source_other) ic = ic + 1; other_index = ic; other_nc_index = other_index
        ic = ic + 1; traffic_exhaust_index = ic; traffic_exhaust_nc_index = traffic_exhaust_index
        ic = ic + 1; traffic_nonexhaust_index = ic; traffic_nonexhaust_nc_index = traffic_nonexhaust_index

        n_source_index = ic

        ic = ic + 1; traffic_gasoline_nc_index = ic
        ic = ic + 1; traffic_diesel_nc_index = ic
        ic = ic + 1; traffic_gas_nc_index = ic
        ic = ic + 1; publicpower_point_nc_index = ic
        ic = ic + 1; publicpower_area_nc_index = ic
        ic = ic + 1; extrasource_nc_index = ic

        n_source_nc_index = ic

        allocate(compound_source_index(n_compound_index, n_source_index))
        allocate(emission_subgrid_loop_index(2, n_source_index))
        allocate(init_emission_subgrid_loop_index(2, n_source_index))
        allocate(n_use_subgrid_levels(n_source_index))
        allocate(emission_subgrid_dim(n_dim_index, n_source_index))
        allocate(emission_subgrid_delta(2, n_source_index))
        allocate(emission_subgrid_min(2, n_source_index))
        allocate(emission_subgrid_max(2, n_source_index))
        allocate(init_emission_subgrid_dim(n_dim_index, n_source_index))
        allocate(init_emission_subgrid_delta(2, n_source_index))
        allocate(init_emission_subgrid_min(2, n_source_index))
        allocate(init_emission_subgrid_max(2, n_source_index))
        allocate(emission_buffer_index(2, n_source_index))
        allocate(emission_buffer_size(2, n_source_index))
        allocate(init_emission_buffer_index(2, n_source_index))
        allocate(init_emission_buffer_size(2, n_source_index))
        allocate(by(n_source_index, n_possible_subsource))
        allocate(az(n_source_index, n_possible_subsource))
        allocate(bz(n_source_index, n_possible_subsource))
        allocate(sig_y_0(n_source_index, n_possible_subsource))
        allocate(sig_z_0(n_source_index, n_possible_subsource))
        allocate(proxy_emission_file_index(n_source_index))
        allocate(emission_file_index(n_source_index))
        allocate(proxy_file_index(n_source_index))
        allocate(proxy_integral_file_index(n_source_index))
        allocate(emep_subgrid_file_index(n_source_index))
        allocate(emep_subgrid_nonlocal_file_index(n_source_index))
        allocate(emep_subgrid_local_file_index(n_source_index))
        allocate(emep_subgrid_frac_file_index(n_source_index))
        allocate(subgrid_local_file_index(n_source_index))
        allocate(subgrid_total_file_index(n_source_index))
        allocate(emep_additional_subgrid_nonlocal_file_index(n_source_index))
        allocate(emep_additional_subgrid_local_file_index(n_source_index))
        allocate(emep_subgrid_semilocal_file_index(n_source_index))
        allocate(subgrid_sourcetotal_inregion_file_index(n_source_index))
        allocate(subgrid_sourcetotal_file_index(n_source_index))
        allocate(use_subgrid_file_index(n_source_index))
        allocate(emep_emission_subgrid_file_index(n_source_index))
        allocate(unit_conversion(n_source_index))
        allocate(emission_factor_conversion(n_compound_nc_index, n_source_index, n_possible_subsource))
        allocate(local_fraction_grid_for_EMEP_grid_interpolation_source(n_source_index))
        allocate(EMEP_grid_interpolation_size_source(n_source_index))

        allocate(convert_GNFR_to_uEMEP_sector_index(n_source_nc_index))
        allocate(source_file_postfix(n_source_nc_index))
        allocate(save_EMEP_source(n_source_nc_index))
        allocate(source_file_str(n_source_nc_index))
        allocate(uEMEP_to_EMEP_sector(n_source_nc_index))
        allocate(uEMEP_to_EMEP_sector_str(n_source_nc_index))
        allocate(uEMEP_to_EMEP_emis_sector_str(n_source_nc_index))
        allocate(var_name_nc(num_var_nc_name, n_pollutant_nc_index, n_source_nc_index))
        allocate(calculate_source(n_source_nc_index))
        allocate(calculate_EMEP_source(n_source_nc_index))
        allocate(make_EMEP_grid_emission_data(n_source_nc_index))
        allocate(replace_EMEP_local_with_subgrid_local(n_source_nc_index))
        allocate(use_emission_positions_for_auto_subgrid_flag(n_source_index))
        allocate(use_trajectory_flag(n_source_index))
        allocate(save_emissions_for_EMEP(n_source_index))
        allocate(n_subsource(n_source_index))
        allocate(uEMEP_to_EMEP_replace_sector(n_source_nc_index))
        allocate(convert_uEMEP_to_GNFR_sector_index(n_source_nc_index))
        allocate(h_emis(n_source_index, n_possible_subsource))
        allocate(sig_y_00(n_source_index, n_possible_subsource))
        allocate(sig_z_00(n_source_index, n_possible_subsource))
        allocate(emission_factor(n_compound_nc_index, n_source_index, n_possible_subsource))
        allocate(z_rec(n_source_index, n_possible_subsource))
        allocate(ay(n_source_index, n_possible_subsource))
        allocate(landuse_proxy_weighting(n_source_index, n_clc_landuse_index))
        allocate(scale_GNFR_emission_source(n_source_index))

        unit_conversion = 1.0
        emission_factor_conversion = 0.0
        local_fraction_grid_for_EMEP_grid_interpolation_source = 1
        EMEP_grid_interpolation_size_source = 1.0
        
        save_EMEP_source = .false.
        source_file_str = ''
        uEMEP_to_EMEP_sector = 0
        uEMEP_to_EMEP_sector_str = ''
        uEMEP_to_EMEP_emis_sector_str = ''
        calculate_source = .false.
        calculate_EMEP_source = .false.
        make_EMEP_grid_emission_data = .false.
        replace_EMEP_local_with_subgrid_local = .false.
        use_emission_positions_for_auto_subgrid_flag = .false.
        use_trajectory_flag = .false.
        save_emissions_for_EMEP = .false.
        n_subsource = 1
        uEMEP_to_EMEP_replace_sector = 0
        emission_factor = 1.0
        landuse_proxy_weighting = 0.0
        scale_GNFR_emission_source = 1.0

    end subroutine set_sources

end module uEMEP_set_sources