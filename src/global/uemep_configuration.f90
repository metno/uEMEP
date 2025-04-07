module uemep_configuration

    use uEMEP_definitions, only: n_population_index, num_var_nc_name, n_pollutant_nc_index, n_source_nc_index, &
        n_compound_nc_index, num_var_population_nc, num_var_landuse_nc, n_source_index, UTM_projection_index, &
        LCC_projection_index, population_index, max_n_local_fraction_grids, n_dim_index, n_possible_subsource, &
        n_clc_landuse_index, num_pollen_nc, num_var_pollen_nc
    use uemep_constants, only: NODATA_value
    use uemep_logger

    implicit none

    ! Configuration file name entered in command line
    integer, parameter :: n_max_config_files = 10
    character(256) :: name_config_file(n_max_config_files) = ''
    character(256) :: emission_date_str = ''
    integer :: n_config_files = 0

    character(256) :: config_date_str = ''
    character(256) :: filename_log_file = 'uEMEP_log.txt'
    character(256) :: pathname_log_file = ''
    character(256) :: file_tag
    character(256) :: replacement_date_str = '<>'
    character(256) :: replacement_yesterday_date_str = '[]'
    character(256) :: replacement_hour_str = '<>'
    character(256) :: NORTRIP_replacement_hour_str = '<>'
    character(256) :: input_comp_name
    character(256) :: pathname_output_grid
    character(256) :: filename_date_output_grid = '<replace_date>_<replace_hour>'
    character(256) :: pathname_rl(2) ! Path name for input roadlink files
    character(256) :: filename_rl(2) ! File name for input roadlink files
    character(256) :: pathname_mrl(50) ! Path name for multiple road link (mrl) files
    character(256) :: filename_mrl(50) ! File name for multiple road link (mrl) files
    character(256) :: pathname_EMEP(4) ! Path name for input EMEP files
    character(256) :: filename_EMEP(4) ! File name for input EMEP files
    character(256) :: pathfilename_EMEP(4)  ! Combined path and file name for input EMEP files
    character(256) :: original_filename_EMEP(4)
    character(256) :: original_pathname_EMEP(4)
    character(256) :: filename_ship(2) ! File name for shipping ais files
    character(256) :: pathname_ship(2) ! Path name for shipping ais files
    character(256) :: pathfilename_ship(2) ! Combined path and file name for shipping ais files
    character(256) :: filename_agriculture(2) ! File name for input agriculture rivm files
    character(256) :: pathname_agriculture(2) ! Path name for input agriculture rivm files
    character(256) :: pathfilename_agriculture(2) ! Combined path and file name for input agriculture rivm files
    character(256) :: filename_emission_rivm(2) ! File name for input emission rivm files
    character(256) :: pathname_emission_rivm(2) ! Path name for input emission rivm files
    character(256) :: pathfilename_emission_rivm(2) ! Combined path and file name for input emission rivm files
    character(256) :: filename_industry(10) ! File name for input industry files
    character(256) :: pathname_industry(10) ! Path name for input industry files
    character(256) :: pathfilename_industry(10) ! Combined path and file name for input industry files
    character(256) :: filename_heating(10) ! File name for input heating files
    character(256) :: pathname_heating(10) ! Path name for input heating files
    character(256) :: pathfilename_heating(10) ! Combined path and file name for input heating files
    character(256) :: filename_population(n_population_index) ! File name for input population files
    character(256) :: pathname_population(n_population_index) ! Path name for input population files
    character(256) :: pathfilename_population(n_population_index) ! Combined path and file name for input population files
    character(256) :: filename_receptor ! File name for the receptor file
    character(256) :: pathname_receptor ! Path name for the receptor file
    character(256) :: pathfilename_receptor  ! Combined path and file name for the receptor file
    character(256) :: filename_timeprofile ! File name for the time profile file
    character(256) :: pathname_timeprofile ! Path name for the time profile file
    character(256) :: pathfilename_timeprofile  ! Combined path and file name for the time profile file
    character(256) :: filename_pollen(num_pollen_nc)
    character(256) :: pathname_pollen(num_pollen_nc)
    character(256) :: pathfilename_pollen(num_pollen_nc)
    ! character(256) :: filename_pollen_proxy(10) ! File name for pollen proxy input
    ! character(256) :: pathname_pollen_proxy(10) ! Path name for pollen proxy input
    character(256) :: alternative_meteorology_type = 'meps'
    character(256) :: pathname_region_id = ''
    character(256) :: filename_region_id = ''
    character(256) :: region_name = ''
    character(256) :: pathfilename_region_id = ''
    character(256) :: pathname_tiles = ''
    character(256) :: filename_tiles = ''
    character(256) :: tile_tag = ''
    character(256) :: save_tile_tag = ''
    character(256) :: inpath_region_heating_scaling = ''
    character(256) :: infile_region_heating_scaling = ''
    character(256) :: pathfilename_region_heating_scaling = ''
    character(256) :: pathname_rl_change = ''
    character(256) :: filename_rl_change = ''
    character(256) :: forecast_hour_str = '00' ! Forecast hour string for writing to files
    character(256) :: NORTRIP_hour_str = '01' ! NORTRIP hour string for writing to files
    character(256) :: pathname_emissions_for_EMEP = ''
    character(256) :: save_emissions_for_EMEP_projection = 'lambert'
    character(256) :: save_emissions_for_EMEP_region = 'NO'
    character(256) :: var_name_nc(num_var_nc_name, n_pollutant_nc_index, n_source_nc_index)
    character(256) :: filename_landuse = ''
    character(256) :: pathname_landuse = ''
    character(256) :: pathfilename_landuse = '' ! Combined path and filename
    ! character(256) :: pathfilename_pollen_proxy ! Combined path and filename for pollen proxy data
    character(256) :: emission_naming_template_str = 'Sec<n>_Emis_mgm2_'
    character(256) :: pathname_boundingbox = ''
    character(256) :: filename_boundingbox = ''
    character(256) :: pathfilename_boundingbox = ''
    character(256) :: select_country_by_name = ''
    character(256) :: comp_name_nc(n_compound_nc_index)
    character(256) :: var_name_population_nc(num_var_population_nc)
    character(256) :: var_name_pollen_nc(num_var_pollen_nc)
    character(256) :: local_fraction_naming_template_str = 'sec<n>_local_fraction'
    character(256) :: finished_filename = ''
    character(256) :: finished_subpath = 'finished/'
    character(256) :: var_name_landuse_nc(num_var_landuse_nc)
    ! Configurations needed for new way to read region mask
    character(256) :: pathname_region_mask = ''
    character(256) :: filename_region_mask = ''
    character(256) :: varname_region_mask = 'region_index'
    ! character(len=256) :: var_name_pollen_proxy_nc(2)

    logical :: hourly_calculations = .false.
    logical :: annual_calculations = .false.
    logical :: use_single_time_loop_flag = .false.
    logical :: reduce_EMEP_region_flag = .false.
    logical :: use_multiple_receptor_grids_flag = .false.
    logical :: reduce_roadlink_region_flag = .true.
    logical :: calculate_source(n_source_nc_index) = .false.
    logical :: calculate_EMEP_source(n_source_nc_index) = .false.
    logical :: make_EMEP_grid_emission_data(n_source_nc_index) = .false.
    logical :: replace_EMEP_local_with_subgrid_local(n_source_nc_index) = .false.
    logical :: subgrid_emission_distribution_flag = .false. ! If true then distributes the EMEP emissions to the existing emission subgrid
    logical :: EMEP_grid_interpolation_simple_flag = .false. ! not used?
    logical :: use_downwind_position_flag = .false. ! If true then searches the upwind EMEP grid position for emissions
    logical :: average_zc_h_in_Kz_flag = .true.
    logical :: use_emission_positions_for_auto_subgrid_flag(n_source_index) = .false.
    logical :: use_receptor_positions_for_auto_subgrid_flag = .false.
    logical :: use_population_positions_for_auto_subgrid_flag = .false.
    logical :: interpolate_subgrids_flag = .false.
    logical :: use_trajectory_flag(n_source_index) = .false.
    logical :: calculate_aggregated_shipping_emissions_flag = .false.
    logical :: use_aggregated_shipping_emissions_flag = .true.
    logical :: calculate_population_exposure_flag = .false.
    logical :: use_last_meteo_in_dispersion = .false.
    logical :: use_meandering_in_dispersion = .false.
    logical :: use_traffic_for_sigma0_flag = .false.
    logical :: use_alternative_meteorology_flag = .false.
    logical :: use_alternative_z0_flag = .false.
    logical :: save_netcdf_file_flag = .false.
    logical :: save_netcdf_receptor_flag = .false.
    logical :: save_netcdf_fraction_as_contribution_flag = .false.
    logical :: calculate_tiling_flag = .false.
    logical :: calculate_region_tiling_flag = .false.
    logical :: use_region_select_and_mask_flag = .false.
    logical :: use_NORTRIP_emission_data = .false.
    logical :: use_NORTRIP_emission_pollutant(n_pollutant_nc_index) = .true.
    logical :: use_RWC_emission_data = .false.
    logical :: include_o3_in_aqi_index = .false.
    logical :: read_weekly_shipping_data_flag = .false.
    logical :: read_monthly_and_daily_shipping_data_flag = .false.
    logical :: use_tunnel_emissions_flag = .true.
    logical :: use_tunnel_deposition_flag = .false.
    logical :: save_emissions_for_EMEP(n_source_index) = .false.
    logical :: save_compounds = .true. ! Output data saving flags
    logical :: save_source_contributions = .true. ! Output data saving flags
    logical :: save_emep_source_contributions = .false. ! Output data saving flags
    logical :: save_emep_additional_source_contributions = .false.
    logical :: save_total_source_contributions = .false.
    logical :: save_local_source_contributions_from_in_region = .false.
    logical :: save_semilocal_source_contributions_from_in_region = .false.
    logical :: save_total_source_contributions_from_in_region = .false.
    logical :: save_no2_source_contributions = .true. ! Output data saving flags
    logical :: save_o3_source_contributions = .true. ! Output data saving flags
    logical :: save_wind_vectors = .false. ! Output data saving flags
    logical :: save_other_meteo = .false. ! Output data saving flags
    logical :: save_emep_original = .true. ! Output data saving flags
    logical :: save_emissions = .false. ! Output data saving flags
    logical :: save_for_chemistry = .false. ! Output data saving flags
    logical :: save_population = .false. ! Output data saving flags
    logical :: save_pollen = .false.
    logical :: save_aqi = .true. ! Output data saving flags
    logical :: save_emep_species = .false. ! Output data saving flags
    logical :: save_deposition = .false. ! Output data saving flags
    logical :: save_seasalt = .false. ! Output data saving flags
    logical :: save_netcdf_average_flag = .false.
    logical :: use_traffic_nox_emission_temperature_dependency = .false.
    logical :: calculate_deposition_flag = .false.
    logical :: calculate_source_depletion_flag = .false.
    logical :: read_landuse_flag = .false.
    logical :: use_plume_dispersion_deposition_flag = .false.
    logical :: adjust_wetdepo_integral_to_lowest_layer_flag = .false.
    logical :: auto_adjustment_for_summertime = .true.
    logical :: use_EMEP_surface_ozone_flag = .false.
    logical :: use_EMEP_surface_compounds_flag = .false.
    logical :: use_water_in_EMEP_surface_pm_flag = .false.
    logical :: save_compounds_as_ascii = .false.
    logical :: use_GNFR_emissions_from_EMEP_flag = .false.
    logical :: use_GNFR19_emissions_from_EMEP_flag = .false.
    logical :: use_alphabetic_GNFR_emissions_from_EMEP_flag = .false.
    logical :: use_emission_naming_template_flag = .false.
    logical :: read_OSM_roadlink_data_flag = .false.
    logical :: no_header_roadlink_data_flag = .false.
    logical :: use_user_specified_sectors_flag = .false.
    logical :: read_population_from_netcdf_flag = .false.
    logical :: read_population_from_netcdf_local_flag = .false.
    logical :: auto_select_OSM_country_flag = .false.
    logical :: select_latlon_centre_domain_position_flag = .false.
    logical :: read_shipping_from_netcdf_flag = .false.
    logical :: read_RWC_file_with_extra_HDD = .false.
    logical :: read_RWC_file_with_extra_HDD_and_height = .false.
    logical :: use_alternative_traveltime_weighting = .false.
    logical :: use_straightline_traveltime_distance = .false.
    logical :: limit_emep_grid_interpolation_region_to_calculation_region = .false.
    logical :: use_local_fraction_naming_template_flag = .false.
    logical :: use_local_fraction_grid_size_in_template_flag = .false.
    logical :: save_traffic_emissions_for_EMEP_as_exhaust_nonexhaust_flag = .false.
    logical :: use_annual_mean_pdf_chemistry_correction = .false.
    logical :: quick_annual_mean_pdf_chemistry_correction = .true.
    logical :: use_landuse_as_proxy = .false.
    logical :: read_rivm_landuse_flag = .false.
    logical :: use_rivm_agricuture_emission_data = .false.
    logical :: read_subgrid_emission_data = .false.
    logical :: use_rivm_subgrid_emission_format = .false.
    logical :: save_EMEP_somo35 = .false.
    logical :: save_EMEP_comax = .false.
    logical :: save_EMEP_o3max = .false.
    logical :: save_EMEP_o3_26th = .false.
    logical :: save_EMEP_so2 = .false.
    logical :: derive_SOA_from_other_species = .false.
    logical :: use_phi_for_invL = .false.
    logical :: trace_emissions_from_in_region = .false.
    logical :: calc_grid_vertical_average_concentration_annual_flag = .false.
    logical :: wind_level_zc_flag = .false. ! This will use the centre of mass wind no matter what type of wind_level_flag is used
    logical :: use_alternative_ppm_variable_for_lf = .false.
    logical :: save_emep_OP_species = .false.


    integer :: start_time_nc_index = 1
    integer :: end_time_nc_index = 1
    integer :: start_time_meteo_nc_index = 1
    integer :: end_time_meteo_nc_index = 1
    integer :: use_receptor_region = 1 ! Surrounding grid region when just running for receptors
    integer :: projection_type = UTM_projection_index
    integer :: EMEP_projection_type = LCC_projection_index
    integer :: utm_zone = 33
    integer :: EMEP_grid_interpolation_flag = 0
    integer :: EMEP_meteo_grid_interpolation_flag = 1
    integer :: EMEP_emission_grid_interpolation_flag = 0
    integer :: local_subgrid_method_flag = 1
    integer :: stability_scheme_flag
    integer :: wind_level_flag = 3
    integer :: wind_level_integral_flag = 3
    integer :: no2_chemistry_scheme_flag = 1
    integer :: no2_background_chemistry_scheme_flag = 0
    integer :: integral_subgrid_step = 1
    integer :: n_subsource(n_source_index) = 1 !Initialise the number of actual emission subsources to 1 for all subsources
    integer :: num_multiple_roadlink_files = 0 ! Number of multiple road link files
    integer :: population_data_type = population_index
    integer :: emission_timeprofile_hour_shift = 1 ! Winter European time
    integer :: region_id = 0
    integer :: region_index = 0
    integer :: HDD_threshold_value = 15
    integer :: n_kz_iterations = 2
    integer :: save_emissions_start_index = 1
    integer :: save_emissions_end_index = 24
    integer :: EMEP_surface_level_nc = 1
    integer :: EMEP_surface_level_nc_2 = 1
    integer :: uEMEP_to_EMEP_replace_sector(n_source_nc_index) = 0
    integer :: local_fraction_grid_size(max_n_local_fraction_grids) = 1
    integer :: n_local_fraction_grids = 1
    integer :: local_fraction_grid_for_EMEP_grid_interpolation = 1
    integer :: local_fraction_grid_for_EMEP_additional_grid_interpolation = 1
    integer :: n_var_av = 100  ! Maximum number of variables to be saved as averages
    integer :: convert_uEMEP_to_GNFR_sector_index(n_source_nc_index)
    integer :: Kz_scheme = 2 ! 1 is O'Brian, 2 is Troen
    integer :: save_emission_subgrid_dim(n_dim_index)
    integer :: alternative_ppm_variable_for_lf_dim = 4


    real :: utm_lon0 = 15.
    real :: ltm_lon0 = 0.
    real :: EMEP_grid_interpolation_size = 1.
    real :: EMEP_additional_grid_interpolation_size = 0.
    real :: traj_step_scale = 2.0
    real :: subgrid_delta(2)
    real :: subgrid_min(2)
    real :: subgrid_max(2)
    real :: init_subgrid_delta(2)
    real :: init_subgrid_min(2)
    real :: init_subgrid_max(2)
    real :: deposition_subgrid_delta(2) = 0.0
    real :: landuse_subgrid_delta(2) = 0.0
    real :: h_emis(n_source_index, n_possible_subsource)
    real :: sig_y_00(n_source_index, n_possible_subsource)
    real :: sigy_0_subgid_width_scale = 0.25
    real :: sig_z_00(n_source_index, n_possible_subsource)
    real :: FF_min_dispersion = 0.1
    real :: ustar_min = 0.001
    real :: hmix_min = 25.0
    real :: hmix_max = 2000.0
    real :: emission_factor(n_compound_nc_index, n_source_index, n_possible_subsource) = 1.0
    real :: ratio_truck_car_emission(n_compound_nc_index) = 10.0
    real :: z_rec(n_source_index, n_possible_subsource) ! Pseudo dispersion parameters
    real :: ay(n_source_index, n_possible_subsource) ! Pseudo dispersion parameters
    real :: replace_invL = NODATA_value  ! Will not replace invL when it has a NODATA value
    real :: replace_hmix = NODATA_value  ! Will not replace mix when it has a NODATA value
    real :: FF_scale = NODATA_value
    real :: FF10_offset = NODATA_value
    real :: DD_offset = NODATA_value
    real :: J_scale = NODATA_value
    real :: replace_z0 = NODATA_value  ! Will not replace z0 when it has a NODATA value
    real :: region_subgrid_delta = 50.0
    real :: max_interpolation_subgrid_size = 1000.0
    real :: DMT_min_value = -20.0 !Minimum allowable daily mean temperature for heating degree day calculation
    real :: integral_subgrid_delta_ref = 0.
    real :: ventilation_factor = 1.0
    real :: min_ADT_ventilation_factor = 0.0
    real :: min_length_ventilation_factor = 0.0
    real :: windspeed_tunnel = 1.0
    real :: lowest_stable_L = 1.0e6
    real :: lowest_unstable_L = -10.0
    real :: tunnel_sig_z_00 = 5.0
    real :: bridge_h_emis = 10.0 ! Bridge height not in use yet
    real :: traffic_nox_emission_temperature_ref_temperature(2)
    real :: traffic_nox_emission_temperature_ref_scaling(2)
    real :: limit_industry_delta = 250.0
    real :: limit_shipping_delta = 250.0
    real :: limit_heating_delta = 250.0
    real :: limit_population_delta = 250.0
    real :: limit_pollen_delta = 250.0
    real :: EMEP_emission_aggregation_period = 1.0
    real :: select_lat_centre_position = 60.0
    real :: select_lon_centre_position = 11.0
    real :: select_domain_width_EW_km = 20.0
    real :: select_domain_height_NS_km = 20.0
    real :: osm_adt_power_scale = 1.0
    real :: romberg_parameters(3) = 0.0
    real :: SRM_parameters(3) = 0.0
    real :: sig_y_scaling_factor = 2.0
    real :: min_proxy_emission_shipping_value = 0.0
    real :: population_power_scale = 1.0
    real :: H_emep = 90.0 ! Height of lowest level in EMEP
    real :: comp_scale_nc(n_compound_nc_index)
    real :: traveltime_power = 1.
    real :: traveltime_scaling = 1.0
    real :: f_no2_emep = 0.1
    real :: ox_sigma_ratio_pdf = 0.0
    real :: nox_sigma_ratio_pdf = 0.0
    real :: max_bin_pdf = 1000.0
    real :: min_bin_pdf = 0.0001
    real :: log10_step_bin_pdf = 0.05
    real :: landuse_proxy_weighting(n_source_index, n_clc_landuse_index) = 0.0
    real :: scale_GNFR_emission_source(n_source_index) = 1.0
    real :: subgrid_receptor_offset(2) = 0.0
    real :: z_invL = 10.0
    real :: save_emission_subgrid_min(2)  !Only x and y
    real :: save_emission_subgrid_delta(2)


    double precision :: projection_attributes(10)
    double precision :: EMEP_projection_attributes(10)

end module uemep_configuration