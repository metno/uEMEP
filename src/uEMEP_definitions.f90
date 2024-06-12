module uEMEP_definitions
!UEMEP_definitions
!Defines variables and indexes used in the uEMEP routines
!Bruce Rolstad Denby 09.11.2016

    implicit none

    ! Directory seperator for linux (/) or windows (\)
    character(1) :: slash = '\'

    character(256) :: model_version_str = 'uEMEP_vx.x'

    logical :: read_EMEP_only_once_flag = .false. ! Note can lead to virtual memory overflow

    integer :: unit_logfile = 0
    integer :: unit_finishedfile = 20
    integer :: n_roadlinks = 0
    integer :: n_roadlinks_major = 0
    integer :: n_roadlinks_major_selected = 0
    logical :: calculate_EMEP_additional_grid_flag = .false.

    ! Single loop time variables
    integer :: t_loop
    integer :: start_time_loop_index = 1
    integer :: end_time_loop_index = 1

    ! Dimensions of the netcdf files that are read
    integer, parameter :: num_dims_nc = 6 ! (x,y,z,t,xpos,ypos)
    integer, parameter :: num_dims_meteo_nc = 4 ! (x,y,z,t)

    integer :: dim_length_nc(num_dims_nc)
    integer :: dim_start_nc(num_dims_nc) = [1, 1, 1, 1, 1, 1] ! start at first value
    integer :: dim_length_meteo_nc(num_dims_meteo_nc)
    integer :: dim_start_meteo_nc(num_dims_meteo_nc) = [1, 1, 1, 1] ! start at first value
    integer :: dim_length_EMEP_nc(num_dims_nc)
    integer :: dim_start_EMEP_nc(num_dims_nc)

    ! Dimensions of the netcdf files that are used
    integer :: end_dim_nc(num_dims_nc)
    integer :: start_dim_nc(num_dims_nc)

    integer, parameter :: lon_nc_index = 1
    integer, parameter :: lat_nc_index = 2
    integer, parameter :: ugrid_nc_index = 3
    integer, parameter :: vgrid_nc_index = 4
    integer, parameter :: FF10_nc_index = 5
    integer, parameter :: FFgrid_nc_index = 6
    integer, parameter :: inv_FFgrid_nc_index = 7
    integer, parameter :: inv_FF10_nc_index = 8
    integer, parameter :: hmix_nc_index = 9
    integer, parameter :: kz_nc_index = 10
    integer, parameter :: invL_nc_index = 11
    integer, parameter :: ustar_nc_index = 12
    integer, parameter :: logz0_nc_index = 13
    integer, parameter :: J_nc_index = 14
    integer, parameter :: conc_nc_index = 15
    ! integer, parameter :: frac_nc_index = 16 ! Set with the loop_index (multiple emep lf grids)
    ! integer, parameter :: local_nc_index = 17 ! Set with the loop_index (multiple emep lf grids)
    integer, parameter :: emis_nc_index = 18
    integer, parameter :: x_nc_index = 19
    integer, parameter :: y_nc_index = 20
    integer, parameter :: ZTOP_nc_index = 21
    integer, parameter :: u10_nc_index = 22
    integer, parameter :: v10_nc_index = 23
    integer, parameter :: uw_nc_index = 24
    integer, parameter :: vw_nc_index = 25
    integer, parameter :: Hflux_nc_index = 26
    integer, parameter :: t2m_nc_index = 27
    integer, parameter :: precip_nc_index = 28
    integer, parameter :: wetdepo_nc_index = 29
    integer, parameter :: drydepo_nc_index = 30
    integer, parameter :: phi_nc_index = 31
    integer :: num_var_nc  !This will be set later when the number local fraction grids is known
    integer, parameter :: num_var_nc_start = 31 ! number of readable variables
    integer, parameter :: num_var_nc_name = 51 ! number of possible variable names, add 20 to include any extra local fraction grids
    integer, parameter :: num_var_meteo_nc = num_var_nc_start

    integer, parameter :: num_lc_var_nc_start = 2 ! number of readable local contribution variables

    integer :: frac_nc_index = num_var_nc_start + 1
    integer :: local_nc_index = num_var_nc_start + 2
    integer :: lc_frac_nc_index = 1
    integer :: lc_local_nc_index = 2
    integer :: num_lc_var_nc = num_lc_var_nc_start

    integer :: compound_index

    integer, parameter :: no2_nc_index = 1
    integer, parameter :: nox_nc_index = 2
    integer, parameter :: pm25_nc_index = 3
    integer, parameter :: pm10_nc_index = 4
    integer, parameter :: nh3_nc_index = 5
    integer, parameter :: o3_nc_index = 6
    integer, parameter :: so2_nc_index = 7
    integer, parameter :: pmex_nc_index = 8

    ! Additional for nh3 compounds
    integer, parameter :: nh4_nc_index = 9

    ! Additional NORTRIP specific compounds
    integer, parameter :: pm25_sand_nc_index = 10
    integer, parameter :: pm10_sand_nc_index = 11
    integer, parameter :: pm25_salt_nc_index = 12
    integer, parameter :: pm10_salt_nc_index = 13

    ! Additional compounds. Need to change n_compound_index as well
    integer, parameter :: c6h6_nc_index = 14
    integer, parameter :: bap_nc_index = 15
    integer, parameter :: co_nc_index = 16
    integer, parameter :: somo35_nc_index = 17
    integer, parameter :: comax_nc_index = 18
    integer, parameter :: o3max_nc_index = 19
    integer, parameter :: o3_26th_nc_index = 20
    integer, parameter :: n_compound_nc_index = 20

    ! These are only used in names but need to change the variable n_pollutant_nc_index to fit these!
    integer, parameter :: pmco_nc_index = 21
    integer, parameter :: all_nc_index = 22
    integer, parameter :: pm_nc_index = 23
    integer, parameter :: all_sand_nc_index = 24
    integer, parameter :: all_salt_nc_index = 25
    integer, parameter :: all_sand_salt_nc_index = 26
    integer, parameter :: all_totals_nc_index = 27
    integer, parameter :: aaqd_totals_nc_index = 28
    integer, parameter :: gp_totals_nc_index = 29
    integer, parameter :: op_totals_nc_index = 30

    ! These must be the same as the subgrid source indexes. Should probably just use the one
    integer, parameter :: allsource_nc_index = 1
    integer, parameter :: traffic_nc_index = 2
    integer, parameter :: shipping_nc_index = 3
    integer, parameter :: heating_nc_index = 4
    integer, parameter :: agriculture_nc_index = 5
    integer, parameter :: industry_nc_index = 6

    ! All the other GNFR emissions
    integer, parameter :: publicpower_nc_index = 7
    integer, parameter :: fugitive_nc_index = 8
    integer, parameter :: solvents_nc_index = 9
    integer, parameter :: aviation_nc_index = 10
    integer, parameter :: offroad_nc_index = 11
    integer, parameter :: waste_nc_index = 12
    integer, parameter :: livestock_nc_index = 13
    integer, parameter :: other_nc_index = 14
    integer, parameter :: traffic_exhaust_nc_index = 15
    integer, parameter :: traffic_nonexhaust_nc_index = 16
    integer, parameter :: traffic_gasoline_nc_index = 17
    integer, parameter :: traffic_diesel_nc_index = 18
    integer, parameter :: traffic_gas_nc_index = 19
    integer, parameter :: publicpower_point_nc_index = 20
    integer, parameter :: publicpower_area_nc_index = 21
    integer, parameter :: extrasource_nc_index = 22
    integer, parameter :: n_source_nc_index = 22

    integer :: convert_GNFR_to_uEMEP_sector_index(n_source_nc_index)

    ! Loop for all pollutants to be calculated
    integer :: pollutant_index
    integer, parameter :: n_pollutant_nc_index = 30 ! Includes the addition naming indexes index
    integer :: n_pollutant_loop = 1
    integer :: n_emep_pollutant_loop = 1
    integer :: pollutant_loop_index(n_pollutant_nc_index)
    integer :: pollutant_loop_back_index(n_pollutant_nc_index)

    ! Set pollutant compound loop
    integer :: n_pollutant_compound_loop(n_pollutant_nc_index) = 0
    integer :: pollutant_compound_loop_index(n_pollutant_nc_index, n_compound_nc_index) = 0

    ! Set the pollutant reference index for meteo data. All meteo data read from EMEP uses this index
    integer :: meteo_p_loop_index = 1
    character(256) :: pollutant_file_str(n_pollutant_nc_index) = ''

    character(256) :: dim_name_nc(num_dims_nc)
    character(256) :: var_name_meteo_nc(num_var_meteo_nc)
    character(256) :: dim_name_meteo_nc(num_dims_meteo_nc)

    integer, parameter :: num_var_population_nc = 2
    integer, parameter :: num_dims_population_nc = 2
    integer, parameter :: population_nc_index = 1 ! population_nc_index=1 assumes population files in lat and lon
    integer, parameter :: dwelling_nc_index = 2
    character(256) :: dim_name_population_nc(num_dims_population_nc)

    integer, parameter :: num_var_shipping_nc = 1 ! assumes shipping files in lat and lon
    integer, parameter :: num_dims_shipping_nc = 2
    character(256) :: var_name_shipping_nc(num_var_shipping_nc)
    character(256) :: dim_name_shipping_nc(num_dims_shipping_nc)

    integer, parameter :: num_var_landuse_nc = 1 ! assumes landuse files in lat and lon
    integer, parameter :: num_dims_landuse_nc = 2
    character(256) :: dim_name_landuse_nc(num_dims_landuse_nc)
    integer :: dim_length_landuse_nc(num_dims_landuse_nc)
    integer :: dim_start_landuse_nc(num_dims_landuse_nc) = [1, 1] ! start at first value

    ! Dimension netcdf fields
    integer, parameter :: x_dim_nc_index = 1
    integer, parameter :: y_dim_nc_index = 2
    integer, parameter :: z_dim_nc_index = 3
    integer, parameter :: time_dim_nc_index = 4
    integer, parameter :: xdist_dim_nc_index = 5
    integer, parameter :: ydist_dim_nc_index = 6

    ! Declare netcdf files
    real, allocatable :: var1d_nc(:, :)
    real, allocatable :: var2d_nc(:, :, :)
    real, allocatable :: var3d_nc(:, :, :, :, :, :)
    real, allocatable :: var4d_nc(:, :, :, :, :, :, :)
    real, allocatable :: lc_var3d_nc(:, :, :, :, :, :, :, :) ! Netcdf local contribution array (idist,jdist,i,j,k,t,type,source)
    real, allocatable :: lc_var4d_nc(:, :, :, :, :, :, :, :, :) ! Netcdf local contribution array (idist,jdist,i,j,k,t,type,source)
    real, allocatable :: comp_var3d_nc(:, :, :, :) ! Netcdf additional compounds array (i,j,t,compound)
    real, allocatable :: comp_var4d_nc(:, :, :, :, :) ! Netcdf additional compounds array (i,j,k,t,compound)

    ! Alternative meteorological input data netcdf files. No source element
    real, allocatable :: meteo_var1d_nc(:, :)
    real, allocatable :: meteo_var2d_nc(:, :, :)
    real, allocatable :: meteo_var3d_nc(:, :, :, :)
    real, allocatable :: meteo_var4d_nc(:, :, :, :, :)

    real :: dgrid_nc(2)
    real :: meteo_dgrid_nc(2)
    double precision, allocatable :: val_dim_nc(:, :)
    double precision, allocatable :: val_dim_meteo_nc(:, :)
    character(256), allocatable :: unit_dim_nc(:)
    character(256), allocatable :: unit_dim_meteo_nc(:)

    integer :: surface_level_nc

    ! Indicates the centre of the EMEP local fraction distance array. Set in uEMEP_read_EMEP and used in uEMEP_subgrid_EMEP
    ! Source of a significant bug up to febuary 2021 when the EMEP lf grid used in the interpolation were not the same
    integer :: xdist_centre_nc = 0
    integer :: ydist_centre_nc = 0

    ! Declare road link data arrays
    real, allocatable :: inputdata_rl(:, :)
    integer, allocatable :: inputdata_int_rl(:, :)
    real, allocatable :: inputdata_rl_emissions(:, :, :)
    logical, allocatable :: valid_link_flag(:)

    ! Road link (rl) indexes
    integer, parameter :: x1_rl_index = 1
    integer, parameter :: x2_rl_index = 2
    integer, parameter :: y1_rl_index = 3
    integer, parameter :: y2_rl_index = 4
    integer, parameter :: x0_rl_index = 5
    integer, parameter :: y0_rl_index = 6
    integer, parameter :: lon1_rl_index = 7
    integer, parameter :: lon2_rl_index = 8
    integer, parameter :: lat1_rl_index = 9
    integer, parameter :: lat2_rl_index = 10
    integer, parameter :: lon0_rl_index = 11
    integer, parameter :: lat0_rl_index = 12
    integer, parameter :: length_rl_index = 13
    integer, parameter :: z0_rl_index = 14
    integer, parameter :: width_rl_index = 15
    integer, parameter :: adt_rl_index = 16
    integer, parameter :: speed_rl_index = 17
    integer, parameter :: hdv_rl_index = 18
    integer, parameter :: tunnel_length_rl_index = 19
    integer, parameter :: num_var_rl = 19

    integer, parameter :: id_rl_index = 1
    integer, parameter :: roadtype_rl_index = 2
    integer, parameter :: nlanes_rl_index = 3
    integer, parameter :: major_index_rl_index = 4
    integer, parameter :: num_int_rl = 4

    ! Declare file and path names for input roadlink files
    character(256) :: pathfilename_rl(2)

    ! Filenames for multiple road link (mrl) files
    character(256) :: pathfilename_mrl(50)

    ! Input shipping (ship) indexes
    integer, parameter :: lon_ship_index = 1 !! Ship longitude index
    integer, parameter :: lat_ship_index = 2 !! Ship latitude index
    integer, parameter :: nox_ship_index = 3 !! Ship total NOx emission index
    integer, parameter :: tpm_ship_index = 4 !! Ship total particulate matter emission
    integer, parameter :: num_var_ship = 4 !! Number of ship variables

    integer, parameter :: idloyds_ship_index = 1
    integer, parameter :: idnorwegian_ship_index = 2
    integer, parameter :: num_int_ship = 2

    integer, parameter :: date_ship_index = 1
    integer, parameter :: time_ship_index = 2
    integer :: num_char_ship ! TODO: Should this be set?
    integer, parameter :: num_char_rl = 2 ! TODO: Strange position, move to road links?

    ! Input agriculture
    integer, parameter :: lon_agriculture_index = 1
    integer, parameter :: lat_agriculture_index = 2
    integer, parameter :: nh3_agriculture_index = 3
    integer, parameter :: num_var_agriculture = 3

    ! Declaration for all filenames used to save gridded data produced by uEMEP. Up to 200 file names but of course can be more
    integer, parameter :: n_filenames_grid = 500
    character(256) :: filename_grid(n_filenames_grid)
    character(256) :: pathname_grid(n_filenames_grid)
    character(256) :: pathfilename_grid(n_filenames_grid)  ! Combined path and filename

    ! Declare subgrid variable indexes for 'subgrid' array
    integer, parameter :: proxy_subgrid_index = 1
    integer, parameter :: proxy_integral_subgrid_index = 2
    integer, parameter :: scaling_factor_subgrid_index = 3
    integer, parameter :: local_subgrid_index = 4
    integer, parameter :: nonlocal_subgrid_index = 5
    integer, parameter :: total_subgrid_index = 6
    integer, parameter :: emep_subgrid_index = 7
    integer, parameter :: emep_frac_subgrid_index = 8
    integer, parameter :: emep_local_subgrid_index = 9
    integer, parameter :: emep_nonlocal_subgrid_index = 10
    integer, parameter :: emep_additional_local_subgrid_index = 11
    integer, parameter :: emep_additional_nonlocal_subgrid_index = 12
    integer, parameter :: drydepo_local_subgrid_index = 13
    integer, parameter :: wetdepo_local_subgrid_index = 14
    integer, parameter :: drydepo_nonlocal_subgrid_index = 15
    integer, parameter :: wetdepo_nonlocal_subgrid_index = 16
    integer, parameter :: n_subgrid_index = 16

    ! Declare meteo subgrid variables. Does nothave to be the same as the nc version
    integer, parameter :: ugrid_subgrid_index = 1
    integer, parameter :: vgrid_subgrid_index = 2
    integer, parameter :: FF10_subgrid_index = 3
    integer, parameter :: FFgrid_subgrid_index = 4
    integer, parameter :: inv_FFgrid_subgrid_index = 5
    integer, parameter :: inv_FF10_subgrid_index = 6
    integer, parameter :: hmix_subgrid_index = 7
    integer, parameter :: kz_subgrid_index = 8
    integer, parameter :: invL_subgrid_index = 9
    integer, parameter :: ustar_subgrid_index = 10
    integer, parameter :: logz0_subgrid_index = 11
    integer, parameter :: J_subgrid_index = 12
    integer, parameter :: t2m_subgrid_index = 13
    integer, parameter :: cos_subgrid_index = 14
    integer, parameter :: sin_subgrid_index = 15
    integer, parameter :: precip_subgrid_index = 16
    integer, parameter :: u10_subgrid_index = 17
    integer, parameter :: v10_subgrid_index = 18
    integer, parameter :: phi_index = 19
    integer, parameter :: n_meteo_subgrid_index = 19

    ! Declare compound indexes for the subgrid. Same as nc_index values for compounds. Must be converted when necessary
    integer, parameter :: no2_index = 1
    integer, parameter :: nox_index = 2
    integer, parameter :: pm25_index = 3
    integer, parameter :: pm10_index = 4
    integer, parameter :: nh3_index = 5
    integer, parameter :: o3_index = 6
    integer, parameter :: so2_index = 7
    integer, parameter :: pmex_index = 8
    integer :: no_index !! TODO Does not seem to be defined or used anywhere, delete?

    ! Additional for nh3 compounds
    integer, parameter :: nh4_index = 9
    integer, parameter :: pm25_sand_index = 10
    integer, parameter :: pm10_sand_index = 11
    integer, parameter :: pm25_salt_index = 12
    integer, parameter :: pm10_salt_index = 13

    ! Additional compounds
    integer, parameter :: c6h6_index = 14
    integer, parameter :: bap_index = 15
    integer, parameter :: co_index = 16
    integer, parameter :: somo35_index = 17
    integer, parameter :: comax_index = 18
    integer, parameter :: o3max_index = 19
    integer, parameter :: o3_26th_index = 20
    integer, parameter :: n_compound_index = 20

    ! Declare source indexes (type_source) must be the same as source_nc_index
    integer, parameter :: allsource_index = 1
    integer, parameter :: traffic_index = 2
    integer, parameter :: shipping_index = 3
    integer, parameter :: heating_index = 4
    integer, parameter :: agriculture_index = 5
    integer, parameter :: industry_index = 6

    ! All the other GNFR emissions
    integer, parameter :: publicpower_index = 7
    integer, parameter :: fugitive_index = 8
    integer, parameter :: solvents_index = 9
    integer, parameter :: aviation_index = 10
    integer, parameter :: offroad_index = 11
    integer, parameter :: waste_index = 12
    integer, parameter :: livestock_index = 13
    integer, parameter :: other_index = 14
    integer, parameter :: traffic_exhaust_index = 15
    integer, parameter :: traffic_nonexhaust_index = 16
    integer, parameter :: n_source_index = 16
    integer, parameter :: n_source_calculate_index = 14
    integer :: compound_source_index(n_compound_index, n_source_index)

    character(256) :: source_file_postfix(n_source_nc_index)
    logical :: save_EMEP_source(n_source_nc_index) = .false.

    integer, parameter :: x_dim_index = 1
    integer, parameter :: y_dim_index = 2
    integer, parameter :: t_dim_index = 3
    integer, parameter :: n_dim_index = 3

    ! Target redistribution grid, the same for all sources and compounds
    integer :: subgrid_dim(n_dim_index)
    real :: subgrid_proj_min(2)
    real :: subgrid_proj_max(2)
    integer :: init_subgrid_dim(n_dim_index)
    real, allocatable :: subgrid(:, :, :, :, :, :) ! subgrid (i,j,t,type_subgrid,type_source,type_subsource)
    real, allocatable :: meteo_subgrid(:, :, :, :)
    real, allocatable :: last_meteo_subgrid(:, :, :)
    real, allocatable :: comp_subgrid(:, :, :, :)
    real, allocatable :: comp_EMEP_subgrid(:, :, :, :)
    real, allocatable :: orig_EMEP_subgrid(:, :, :, :)
    real, allocatable :: x_subgrid(:, :)
    real, allocatable :: y_subgrid(:, :)
    real, allocatable :: lon_subgrid(:, :)
    real, allocatable :: lat_subgrid(:, :)
    real, allocatable :: xproj_subgrid(:, :)
    real, allocatable :: yproj_subgrid(:, :)
    real, allocatable :: traveltime_subgrid(:, :, :, :, :)
    real, allocatable :: exposure_subgrid(:, :, :, :, :)
    integer :: subgrid_loop_index(2) ! Number of target subgrids to loop through, limitted by the size of the EMEP grid
    integer :: integral_subgrid_loop_index(2) ! Number of integral subgrids to loop through, limitted by the size of the EMEP grid
    integer :: emission_subgrid_loop_index(2, n_source_index) ! Number of emission subgrids to loop through, limitted by the size of the EMEP grid
    integer :: init_emission_subgrid_loop_index(2, n_source_index) ! Number of emission subgrids to loop through, limitted by the size of the EMEP grid
    logical, allocatable :: use_subgrid(:, :, :) ! Specifies if calculations are to be made at a particular set of target subgrids or not
    integer, allocatable :: use_subgrid_val(:, :, :) ! Same as use_subgrid but given a value to indicate if it is in the buffer zone of a region ore not
    integer, allocatable :: use_subgrid_interpolation_index(:, :, :) ! Assigns the resolution level for auto gridding to the target grid
    integer :: n_use_subgrid_levels(n_source_index)

    real :: loop_index_scale = 1.5
    real :: buffer_index_scale = 1.5

    ! Emission subgrid per source and subsource. Can be time dependent
    ! Each source type has its own x and y and dim
    ! Each source may be of lesser dimensions than the total array size (which is the same as the target grid)
    integer, parameter :: n_possible_subsource = 2
    character(2) :: subsource_str(n_possible_subsource)

    integer, parameter :: emission_h_index = 1
    integer :: emission ! TODO Does not appear to be defined or used anywhere, delete?
    integer, parameter :: emission_sigz00_index = 2
    integer, parameter :: emission_sigy00_index = 3
    integer, parameter :: emission_minFF_index = 4
    integer, parameter :: n_emission_index = 4

    integer :: emission_subgrid_dim(n_dim_index, n_source_index)
    integer :: emission_max_subgrid_dim(n_dim_index)
    real :: emission_subgrid_delta(2, n_source_index)
    real :: emission_subgrid_min(2, n_source_index)
    real :: emission_subgrid_max(2, n_source_index)
    integer :: init_emission_subgrid_dim(n_dim_index, n_source_index)
    real :: init_emission_subgrid_delta(2, n_source_index)
    real :: init_emission_subgrid_min(2, n_source_index)
    real :: init_emission_subgrid_max(2, n_source_index)

    real, allocatable :: emission_subgrid(:, :, :, :, :) ! emission_subgrid (i,j,t,n_source,n_subsource)
    real, allocatable :: proxy_emission_subgrid(:, :, :, :) ! No time dependence
    real, allocatable :: emission_time_profile_subgrid(:, :, :, :, :)
    real, allocatable :: emission_properties_subgrid(:, :, :, :) ! No time dependence and no pollutant dependence. Would then need to do the gaussian calculation for each pollutant

    real, allocatable :: x_emission_subgrid(:, :, :) ! x_emission_subgrid (i,j,n_source)
    real, allocatable :: y_emission_subgrid(:, :, :)
    real, allocatable :: lon_emission_subgrid(:, :, :)
    real, allocatable :: lat_emission_subgrid(:, :, :)
    real, allocatable :: xproj_emission_subgrid(:, :, :)
    real, allocatable :: yproj_emission_subgrid(:, :, :)

    logical :: use_buffer_zone = .true.
    integer :: buffer_index(2)
    real :: buffer_size(2)
    integer :: emission_buffer_index(2, n_source_index)
    real :: emission_buffer_size(2, n_source_index)
    integer :: init_emission_buffer_index(2, n_source_index)
    real :: init_emission_buffer_size(2, n_source_index)
    integer :: integral_buffer_index(2)
    real :: integral_buffer_size(2)

    integer, parameter :: hsurf_integral_subgrid_index = 1
    integer, parameter :: hmix_integral_subgrid_index = 2
    integer, parameter :: hsurf_average_subgrid_index = 3
    integer, parameter :: n_integral_subgrid_index = 3

    integer :: integral_subgrid_dim(n_dim_index)
    real :: integral_subgrid_delta(2) = 0.
    real :: integral_subgrid_min(2)
    real :: integral_subgrid_max(2)
    real, allocatable :: integral_subgrid(:, :, :, :, :, :) ! emission_subgrid (i,j,t,n_integral_type,n_source,n_pollutant)
    real, allocatable :: x_integral_subgrid(:, :)
    real, allocatable :: y_integral_subgrid(:, :)
    real, allocatable :: lon_integral_subgrid(:, :)
    real, allocatable :: lat_integral_subgrid(:, :)
    real, allocatable :: xproj_integral_subgrid(:, :)
    real, allocatable :: yproj_integral_subgrid(:, :)
    real, allocatable :: meteo_nc_xproj_integral_subgrid(:, :)
    real, allocatable :: meteo_nc_yproj_integral_subgrid(:, :)

    integer :: population_subgrid_dim(2)
    real :: population_subgrid_delta(2)
    real :: population_subgrid_min(2)
    real :: population_subgrid_max(2)
    real, allocatable :: population_subgrid(:, :, :)
    real, allocatable :: x_population_subgrid(:, :)
    real, allocatable :: y_population_subgrid(:, :)
    real, allocatable :: lon_population_subgrid(:, :)
    real, allocatable :: lat_population_subgrid(:, :)
    real, allocatable :: xproj_population_subgrid(:, :)
    real, allocatable :: yproj_population_subgrid(:, :)

    integer, allocatable :: crossreference_target_to_emep_subgrid(:, :, :) ! (i,j,dim)
    integer, allocatable :: crossreference_target_to_localfraction_subgrid(:, :, :, :) ! (i,j,dim,n_lf_grids)
    integer, allocatable :: crossreference_integral_to_emep_subgrid(:, :, :) ! (i,j,dim)
    integer, allocatable :: crossreference_integral_to_meteo_nc_subgrid(:, :, :) ! (i,j,dim)
    integer, allocatable :: crossreference_target_to_integral_subgrid(:, :, :) ! (i,j,dim)
    integer, allocatable :: crossreference_target_to_emission_subgrid(:, :, :, :) ! (i,j,dim,n_source)
    integer, allocatable :: crossreference_integral_to_emission_subgrid(:, :, :, :) ! (i,j,dim,n_source)
    integer, allocatable :: crossreference_emission_to_emep_subgrid(:, :, :, :) ! (i,j,dim,n_source)
    integer, allocatable :: crossreference_emission_to_integral_subgrid(:, :, :, :) ! (i,j,dim,n_source)
    integer, allocatable :: crossreference_target_to_population_subgrid(:, :, :) ! (i,j,dim)
    integer, allocatable :: crossreference_emission_to_deposition_subgrid(:, :, :, :) ! (i,j,dim,n_source)
    integer, allocatable :: crossreference_emission_to_landuse_subgrid(:, :, :, :) ! (i,j,dim,n_source)
    integer, allocatable :: crossreference_target_to_deposition_subgrid(:, :, :) ! (i,j,dim)
    integer, allocatable :: crossreference_deposition_to_emep_subgrid(:, :, :) ! (i,j,dim)

    real :: min_link_size = 50.0
    real :: min_adt = 1000.0
    real :: H_meteo = 45.0 ! Height of the gridded meteo values

    integer :: start_month_in_annual_calculations = 1
    integer :: end_month_in_annual_calculations = 12

    ! Pseudo dispersion parameters
    real :: by(n_source_index, n_possible_subsource)
    real :: az(n_source_index, n_possible_subsource)
    real :: bz(n_source_index, n_possible_subsource)
    real :: sig_y_0(n_source_index, n_possible_subsource)
    real :: sig_z_0(n_source_index, n_possible_subsource)

    integer, parameter :: UTM_projection_index = 1
    integer, parameter :: RDM_projection_index = 2
    integer, parameter :: LCC_projection_index = 3
    integer, parameter :: LL_projection_index = 4
    integer, parameter :: LAEA_projection_index = 5
    integer, parameter :: LTM_projection_index = 6
    integer, parameter :: PS_projection_index = 7
    double precision :: population_nc_projection_attributes(10)
    integer :: meteo_nc_projection_type = LCC_projection_index
    double precision :: meteo_nc_projection_attributes(10)
    logical :: use_alternative_LCC_projection_flag = .false.

    ! Filename index for files produced by uEMEP
    integer :: proxy_emission_file_index(n_source_index)
    integer :: emission_file_index(n_source_index)
    integer :: proxy_file_index(n_source_index)
    integer :: proxy_integral_file_index(n_source_index)
    integer :: emep_subgrid_file_index(n_source_index)
    integer :: emep_subgrid_nonlocal_file_index(n_source_index)
    integer :: emep_subgrid_local_file_index(n_source_index)
    integer :: emep_subgrid_frac_file_index(n_source_index)
    integer :: subgrid_local_file_index(n_source_index)
    integer :: subgrid_total_file_index(n_source_index)
    integer :: emep_additional_subgrid_nonlocal_file_index(n_source_index)
    integer :: emep_additional_subgrid_local_file_index(n_source_index)

    ! Filename index for meteorological parameters
    integer :: subgrid_ugrid_file_index
    integer :: subgrid_vgrid_file_index
    integer :: subgrid_u10_file_index
    integer :: subgrid_v10_file_index
    integer :: subgrid_hmix_file_index
    integer :: subgrid_kz_file_index
    integer :: subgrid_logz0_file_index
    integer :: subgrid_invL_file_index
    integer :: subgrid_FF10_file_index
    integer :: subgrid_FFgrid_file_index
    integer :: subgrid_invFF10_file_index
    integer :: subgrid_invFFgrid_file_index
    integer :: subgrid_ustar_file_index
    integer :: subgrid_J_file_index
    integer :: subgrid_meteo_file_index
    integer :: subgrid_DD10_file_index
    integer :: subgrid_DDgrid_file_index
    integer :: subgrid_t2m_file_index

    ! Filename index for grid auto grid parameters
    integer :: use_subgrid_file_index(n_source_index)
    integer :: emep_emission_subgrid_file_index(n_source_index)

    character(256) :: source_file_str(n_source_nc_index) = ''
    real :: unit_conversion(n_source_index) = 1.0

    real :: emission_factor_conversion(n_compound_nc_index, n_source_index, n_possible_subsource) = 0.0

    integer :: weighting_step = 1

    integer :: ref_year_EMEP = 1900
    integer :: ref_year_meteo = 1970

    integer :: n_receptor
    integer :: n_receptor_in
    integer, parameter :: n_receptor_max = 20000
    integer :: n_valid_receptor
    integer :: n_valid_receptor_in
    real :: lon_receptor(n_receptor_max)
    real :: lat_receptor(n_receptor_max)
    real :: x_receptor(n_receptor_max)
    real :: y_receptor(n_receptor_max)
    real :: height_receptor(n_receptor_max)
    real :: lon_receptor_in(n_receptor_max)
    real :: lat_receptor_in(n_receptor_max)
    real :: x_receptor_in(n_receptor_max)
    real :: y_receptor_in(n_receptor_max)
    real :: height_receptor_in(n_receptor_max)
    integer :: i_receptor_subgrid(n_receptor_max)
    integer :: j_receptor_subgrid(n_receptor_max)
    character(256) :: name_receptor(n_receptor_max, 2)
    character(256) :: name_receptor_in(n_receptor_max, 2)
    logical :: use_receptor(n_receptor_max) = .true.
    integer :: valid_receptor_index(n_receptor_max)
    integer :: valid_receptor_inverse_index(n_receptor_max)

    ! Indicies for SSB building and population data
    integer, parameter :: dwelling_index = 1
    integer, parameter :: population_index = 2
    integer, parameter :: establishment_index = 3
    integer, parameter :: school_index = 4
    integer, parameter :: kindergaten_index = 5
    integer, parameter :: home_index = 6
    integer, parameter :: municipality_index = 7
    integer, parameter :: RWC_heating_index = 8
    integer, parameter :: n_population_index = 8
    integer :: population_file_index(n_population_index)

    ! Is a temporary variable set when reading SSB data
    integer :: SSB_data_type = dwelling_index

    integer :: n_receptor_grid
    integer :: start_grid_loop_index
    integer :: end_grid_loop_index
    integer :: g_loop

    ! Correction output time array converting days 1900 to seconds 2000
    integer(kind=4), allocatable :: time_seconds_output(:)

    ! Residential wood combustion variables
    integer :: n_RWC_grids
    real, allocatable :: RWC_grid_emission(:, :)
    real, allocatable :: RWC_grid_HDD(:, :)
    integer(kind=8), allocatable :: RWC_grid_id(:)
    integer, allocatable :: RWC_region_id(:)
    real, allocatable :: RWC_grid_height(:, :)
    real, allocatable :: DMT_EMEP_grid_nc(:, :, :)

    ! Scenario calculator variables

    real :: aqi_hourly_limits(n_compound_index, 1:3)
    real :: aqi_daily_limits(n_compound_index, 1:3)
    real :: aqi_annual_limits(n_compound_index, 1:3)

    ! Special source allocation for no2 based on leaving out the source in the chemistry calculation
    real, allocatable :: comp_source_subgrid(:, :, :, :, :)
    real, allocatable :: comp_source_additional_subgrid(:, :, :, :, :)
    real, allocatable :: comp_source_EMEP_subgrid(:, :, :, :, :)
    real, allocatable :: comp_source_EMEP_additional_subgrid(:, :, :, :, :)

    ! Set the scaling factors for the auto gridding routine
    integer :: use_subgrid_step_delta(0:10)

    integer, parameter :: outside_region_index = -1
    integer, parameter :: outside_interpolation_region_index = -2
    integer, parameter :: inside_region_index = 0

    ! Variables for saving averages
    real, allocatable :: val_array_av(:, :, :)
    integer(8), allocatable :: time_seconds_output_av(:)
    integer :: counter_av = 0

    ! Species variables
    integer, parameter :: pm10_sp_index = 1
    integer, parameter :: pm25_sp_index = 2
    integer, parameter :: pmco_sp_index = 3 ! pmco_sp_index is just for reading
    integer, parameter :: n_pmxx_sp_index = 3

    integer, parameter :: sp_soa_index = 1
    integer, parameter :: sp_sia_index = 2
    integer, parameter :: sp_dust_index = 3
    integer, parameter :: sp_seasalt_index = 4
    integer, parameter :: sp_ffire_index = 5
    integer, parameter :: sp_ppm_index = 6
    integer, parameter :: sp_water_index = 7
    integer, parameter :: sp_pm_index = 8
    integer, parameter :: sp_asoa_index = 9
    integer, parameter :: sp_bsoa_index = 10
    integer, parameter :: n_sp_index = 10

    ! Additional BBOA species used in the OP calculations
    integer, parameter :: sp_BBOA_index = 11
    integer, parameter :: sp_BBOA_RES_index = 12
    integer, parameter :: n_sp_OP_index = 12

    ! These are used just for reading
    integer, parameter :: sp_no3_index = 13
    integer, parameter :: sp_so4_index = 14
    integer, parameter :: sp_nh4_index = 15
    integer, parameter :: sp_dust_sah_index = 16
    integer, parameter :: sp_dust_wb_index = 17
    integer, parameter :: sp_ffire_bc_index = 18
    integer, parameter :: sp_ffire_rem_index = 19

    ! Alternative input names so the other names are reserved for otuput
    integer, parameter :: sp_soa_in_index = 20
    integer, parameter :: sp_sia_in_index = 21
    integer, parameter :: sp_dust_in_index = 22
    integer, parameter :: sp_seasalt_in_index = 23
    integer, parameter :: sp_ffire_in_index = 24
    integer, parameter :: sp_ppm_in_index = 25
    integer, parameter :: sp_water_in_index = 26
    integer, parameter :: sp_pm_in_index = 27

    ! Alternative input names for OP so the other names are reserved for otuput
    integer, parameter :: sp_POM_RES_in_index = 28
    integer, parameter :: sp_EC_RES_NEW_in_index = 29
    integer, parameter :: sp_EC_RES_AGE_in_index = 30
    integer, parameter :: sp_REM_RES_in_index = 31
    integer, parameter :: sp_FFIRE_OM_in_index = 32
    integer, parameter :: sp_FFIRE_BC_in_index = 33
    integer, parameter :: sp_FFIRE_REM_in_index = 34
    integer, parameter :: sp_EC_RES_in_index = 35
    integer, parameter :: sp_asoa_in_index = 36
    integer, parameter :: sp_bsoa_in_index = 37
    integer, parameter :: n_sp_all_index = 37

    real, allocatable :: species_var3d_nc(:, :, :, :, :) ! (x,y,t,n_pmxx_sp_index,n_species_loop_index)
    real, allocatable :: species_EMEP_subgrid(:, :, :, :, :) ! (x,y,t,n_pmxx_sp_index,n_species_loop_index)
    integer :: species_loop_index(n_sp_all_index)
    integer :: n_species_loop_index = n_sp_index ! Variable length of species list, set in uEMEP_set_species_loop

    ! Dimension with the largest possible, including OP
    character(256) :: species_name_nc(n_pmxx_sp_index, n_sp_all_index)

    ! Deposition and land use
    real, allocatable :: orig_EMEP_deposition_subgrid(:, :, :, :, :)
    real, allocatable :: depo_var3d_nc(:, :, :, :, :)
    integer :: deposition_subgrid_dim(n_dim_index)
    real :: deposition_subgrid_min(2)
    real :: deposition_subgrid_max(2)

    integer, parameter :: vd_index = 1
    integer, parameter :: drydepo_index = 2
    integer, parameter :: wetdepo_index = 3
    integer, parameter :: n_deposition_index = 3
    real, allocatable :: deposition_subgrid(:, :, :, :, :) ! deposition_subgrid (i,j,t,n_deposition_index,n_pollutant_loop)
    real, allocatable :: x_deposition_subgrid(:, :)
    real, allocatable :: y_deposition_subgrid(:, :)
    real, allocatable :: lon_deposition_subgrid(:, :)
    real, allocatable :: lat_deposition_subgrid(:, :)
    real, allocatable :: xproj_deposition_subgrid(:, :)
    real, allocatable :: yproj_deposition_subgrid(:, :)
    integer :: deposition_subgrid_loop_index(2)
    integer :: deposition_buffer_index(2)
    real :: deposition_buffer_size(2)

    real :: wetdepo_scavanging_rate(n_compound_index)
    real :: drydepo_vd_default(n_compound_index)

    integer :: landuse_subgrid_dim(n_dim_index)
    real :: landuse_subgrid_min(2)
    real :: landuse_subgrid_max(2) ! Only x and y

    integer, parameter :: temp_conif_index = 1
    integer, parameter :: temp_decid_index = 2
    integer, parameter :: med_needle_index = 3
    integer, parameter :: med_broadleaf_index = 4
    integer, parameter :: temp_crop_index = 5
    integer, parameter :: med_crop_index = 6
    integer, parameter :: root_crop_index = 7
    integer, parameter :: moorland_index = 8
    integer, parameter :: grass_index = 9
    integer, parameter :: medscrub_index = 10
    integer, parameter :: wetlands_index = 11
    integer, parameter :: tundra_index = 12
    integer, parameter :: desert_index = 13 ! TODO Same as water below?
    integer, parameter :: water_index = 13 ! TODO Same as desert above?
    integer, parameter :: ice_index = 14
    integer, parameter :: urban_index = 15
    integer, parameter :: grid_index = 16
    integer, parameter :: clc_index = 17
    integer, parameter :: n_landuse_index = 17

    real, allocatable :: landuse_subgrid(:, :, :) ! landuse_subgrid (i,j,n_landuse_index) Fraction of landuse type
    real, allocatable :: x_landuse_subgrid(:, :)
    real, allocatable :: y_landuse_subgrid(:, :)
    real, allocatable :: lon_landuse_subgrid(:, :)
    real, allocatable :: lat_landuse_subgrid(:, :)
    real, allocatable :: xproj_landuse_subgrid(:, :)
    real, allocatable :: yproj_landuse_subgrid(:, :)
    integer :: landuse_subgrid_loop_index(2)
    integer :: landuse_buffer_index(2)
    real :: landuse_buffer_size(2)

    character(256) :: deposition_name_nc(n_landuse_index, n_compound_nc_index)

    real :: depo_scale_nc(n_compound_nc_index)

    logical :: first_g_loop = .true.

    ! Define the source sector match between uEMEP and EMEP
    integer :: uEMEP_to_EMEP_sector(n_source_nc_index) = 0
    character(2) :: uEMEP_to_EMEP_sector_str(n_source_nc_index) = ''
    character(2) :: uEMEP_to_EMEP_emis_sector_str(n_source_nc_index) = ''

    ! Define the aggregation period for EMEP emissions when these are to be used in calculations. Annual is 365*24=8760 or 8784 for leap years

    integer :: population_nc_projection_type = LL_projection_index

    ! Additional multiple local fraction variables
    integer, parameter :: max_n_local_fraction_grids = 3
    integer :: frac_nc_loop_index(max_n_local_fraction_grids)
    integer :: local_nc_loop_index(max_n_local_fraction_grids)
    integer :: lc_frac_nc_loop_index(max_n_local_fraction_grids)
    integer :: lc_local_nc_loop_index(max_n_local_fraction_grids)
    integer :: min_frac_nc_loop_index
    integer :: max_frac_nc_loop_index
    integer :: min_lc_frac_nc_loop_index
    integer :: max_lc_frac_nc_loop_index
    integer :: convert_frac_to_lc_frac_loop_index(num_var_nc_name)

    integer :: local_fraction_grid_for_EMEP_grid_interpolation_source(n_source_index) = 1
    real :: local_fraction_grid_size_scaling = 1.0
    real :: EMEP_grid_interpolation_size_original = 1.0
    real :: EMEP_grid_interpolation_size_source(n_source_index) = 1.0
    real :: local_fraction_additional_grid_size_scaling = 1.0
    real :: EMEP_additional_grid_interpolation_size_original = 0.0

    character(256) :: finished_file = ''
    character(256) :: finished_file_rec = ''

    ! Landuse proxy

    integer, parameter :: n_clc_landuse_index = 44

    ! Benzene split from VOC
    logical :: extract_benzene_from_voc_emissions = .false.
    real :: benzene_split_voc_in_GNFR_sectors(13) = &
        [0.0449, 0.0212, 0.0668, 0.0084, 0.0, 0.0266, 0.0226, 0.0214, 0.0223, 0.0362, 0.068, 0.0601, 0.068]

    ! Definitions of the grid when saving emissions

    real, allocatable :: subgrid_from_in_region(:, :, :, :, :, :)
    real, allocatable :: EMEP_grid_fraction_in_region(:, :, :, :)
    real, allocatable :: lf_EMEP_grid_fraction_in_region(:, :, :, :, :, :)
    logical, allocatable :: use_subgrid_region(:, :, :) ! Specifies the region emissions will be carried from for subgrid_from_in_region
    real, allocatable :: comp_subgrid_from_in_region(:, :, :, :)

    real, allocatable :: comp_source_subgrid_from_in_region(:, :, :, :, :)
    real, allocatable :: comp_source_additional_subgrid_from_in_region(:, :, :, :, :)
    real, allocatable :: comp_source_EMEP_subgrid_from_in_region(:, :, :, :, :)
    real, allocatable :: comp_source_EMEP_additional_subgrid_from_in_region(:, :, :, :, :)

    ! Setting this to true is for diagnostic puroses. Gived the integrated lowest grid average concentration instead of the receptor

    logical :: save_emep_region_mask = .false.
    logical :: wind_vectors_10m_available = .false.

    ! New definitions needed for the adjusted non-local from-in-region (nlreg) methodology

    ! nlreg_subgrid_region_id will contain the region ID of each uEMEP target subgrid
    integer, allocatable :: nlreg_subgrid_region_id(:, :)
    ! nlreg_EMEP_subsample_region_id will contain the region ID of each subsample of each EMEP grid
    integer, allocatable :: nlreg_EMEP_subsample_region_id(:, :, :, :)
    ! nlreg_regionfraction_per_EMEP_grid will give the fraction of each EMEP grid that is within each region
    real, allocatable :: nlreg_regionfraction_per_EMEP_grid(:, :, :)
    ! nlreg_region_ids specifies the region IDs along the region dimension and will get length nlreg_n_regions
    integer, allocatable :: nlreg_region_ids(:)
    ! nlreg_n_regions is the number of regions. Will be set after reading region mask file
    integer :: nlreg_n_regions = -1
    ! nlreg_n_subsamples_per_EMEP_grid specifies how many subsamples to use along each spatial dimension of an EMEP grid when creating a region mask for each grid
    integer :: nlreg_n_subsamples_per_EMEP_grid = 20

    ! The following arrays will contain the contributions from outside moving window but within region for each uEMEP target subgrid (x,y,t,source,pollutant)
    ! - from the small (1x1) LF grid data
    real, allocatable :: nlreg_subgrid_nonlocal_from_in_region(:, :, :, :, :)
    ! - from the big (additional) LF grid data, only the part outside the 1x1 LF domain
    ! (i.e. this array contains only the contribution NOT covered in the above array, so they can be added)
    real, allocatable :: nlreg_subgrid_nonlocal_from_in_region_additional_increment(:, :, :, :, :)

end module uEMEP_definitions

