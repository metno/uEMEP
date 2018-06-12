!UEMEP_definitions
!Defines variables and indexes used in the uEMEP routines
!Bruce Rolstad Denby 09.11.2016
   
    module uEMEP_definitions
    
    !Directory seperator for linux (/) or windows (\)
    character(1) :: slash='\'

    !Configuration file name entered in command line
    character(256) :: name_config_file(5)=''
    character(256) :: filename_log_file='uEMEP_log.txt'
    character(256) :: pathname_log_file=''
    character(256) :: config_date_str=''
    character(256) :: replacement_date_str='<>'
    integer :: n_config_files=0
    
    logical :: use_single_time_loop_flag=.false.
    logical :: reduce_EMEP_region_flag=.false.
    
    !Nodata value
    real NODATA_value
    parameter (NODATA_value=-999.)
    
    !often used
    !integer i,j,k
    
    integer :: unit_logfile=0
    integer :: n_roadlinks=0
    integer :: utm_zone=33
    real :: utm_lon0=15.
    integer :: EMEP_grid_interpolation_flag=0
    integer :: EMEP_meteo_grid_interpolation_flag=1
    real :: EMEP_grid_interpolation_size=1.
    integer :: local_subgrid_method_flag=1
    logical :: use_emission_positions_for_auto_subgrid_flag=.false.
    
    !Single loop time variables
    integer t_loop
    integer :: start_time_loop_index=1
    integer :: end_time_loop_index=1
    
    !Dimensions of the netcdf files that are read
    integer num_dims_nc 
    parameter (num_dims_nc=6) !(x,y,z,t,xpos,ypos)
    integer num_dims_meteo_nc 
    parameter (num_dims_meteo_nc=4) !(x,y,z,t)
    
    integer dim_length_nc(num_dims_nc)
    integer dim_start_nc(num_dims_nc)
    data dim_start_nc /1, 1, 1, 1, 1, 1/                 ! start at first value
    integer dim_length_meteo_nc(num_dims_meteo_nc)
    integer dim_start_meteo_nc(num_dims_meteo_nc)
    data dim_start_meteo_nc /1, 1, 1, 1/                 ! start at first value

    !Dimensions of the netcdf files that are used
    integer end_dim_nc(num_dims_nc)
    integer start_dim_nc(num_dims_nc)

    integer lat_nc_index,lon_nc_index
    integer ugrid_nc_index,vgrid_nc_index,FF10_nc_index,FFgrid_nc_index,inv_FFgrid_nc_index,inv_FF10_nc_index
    integer hmix_nc_index,kz_nc_index,invL_nc_index,ustar_nc_index,logz0_nc_index,J_nc_index
    integer conc_nc_index,frac_nc_index,local_nc_index,emis_nc_index
    integer x_nc_index,y_nc_index,ZTOP_nc_index
    integer u10_nc_index,v10_nc_index,uw_nc_index,vw_nc_index,Hflux_nc_index
    parameter (lon_nc_index=1,lat_nc_index=2)
    parameter (ugrid_nc_index=3,vgrid_nc_index=4,FF10_nc_index=5,FFgrid_nc_index=6,inv_FFgrid_nc_index=7,inv_FF10_nc_index=8)
    parameter (hmix_nc_index=9,kz_nc_index=10,invL_nc_index=11,ustar_nc_index=12,logz0_nc_index=13,J_nc_index=14)
    parameter (conc_nc_index=15,frac_nc_index=16,local_nc_index=17,emis_nc_index=18)
    parameter (x_nc_index=19,y_nc_index=20,ZTOP_nc_index=21)
    parameter (u10_nc_index=22,v10_nc_index=23,uw_nc_index=24,vw_nc_index=25,Hflux_nc_index=26)
    integer num_var_nc
    parameter (num_var_nc=26)                  ! number of readable variables
    integer num_var_meteo_nc
    parameter (num_var_meteo_nc=num_var_nc)
    
    integer lc_frac_nc_index,lc_local_nc_index
    parameter (lc_frac_nc_index=1,lc_local_nc_index=2)
    integer num_lc_var_nc
    parameter (num_lc_var_nc=2)                  ! number of readable local contribution variables
    
    integer compound_index,compound_frac_index

    integer no2_nc_index,nox_nc_index,pm25_nc_index,pm10_nc_index,nh3_nc_index,o3_nc_index
    parameter (no2_nc_index=1,nox_nc_index=2,pm25_nc_index=3,pm10_nc_index=4,nh3_nc_index=5,o3_nc_index=6)
    integer n_compound_nc_index
    parameter (n_compound_nc_index=6)
    !THese must be the same as the subgrid source indexes. Should probably just use the one
    integer allsource_nc_index,traffic_nc_index,shipping_nc_index,heating_nc_index,agriculture_nc_index,industry_nc_index
    parameter (allsource_nc_index=1,traffic_nc_index=2,shipping_nc_index=3,heating_nc_index=4,agriculture_nc_index=5,industry_nc_index=6)
    integer n_source_nc_index
    parameter (n_source_nc_index=6)
   
    character(256) var_name_nc(num_var_nc,n_compound_nc_index,n_source_nc_index)
    character(256) dim_name_nc(num_dims_nc)
    character(256) var_name_meteo_nc(num_var_meteo_nc)
    character(256) dim_name_meteo_nc(num_dims_meteo_nc)
    character(256) comp_name_nc(n_compound_nc_index)
    character(256) input_comp_name
    real comp_scale_nc(n_compound_nc_index)
    integer :: n_compound_loop = 1
    integer compound_loop_index(n_compound_nc_index)

    !dimension netcdf fields
    integer x_dim_nc_index,y_dim_nc_index,z_dim_nc_index,time_dim_nc_index,xdist_dim_nc_index,ydist_dim_nc_index
    parameter (x_dim_nc_index=1,y_dim_nc_index=2,z_dim_nc_index=3,time_dim_nc_index=4,xdist_dim_nc_index=5,ydist_dim_nc_index=6)

    !Declare netcdf files
    real, allocatable :: var1d_nc(:,:)
    real, allocatable :: var2d_nc(:,:,:)
    real, allocatable :: var3d_nc(:,:,:,:,:)
    real, allocatable :: var4d_nc(:,:,:,:,:,:)
    real, allocatable :: lc_var3d_nc(:,:,:,:,:,:,:)         !Netcdf local contribution array (idist,jdist,i,j,k,t,type,source)
    real, allocatable :: lc_var4d_nc(:,:,:,:,:,:,:,:)       !Netcdf local contribution array (idist,jdist,i,j,k,t,type,source)
    real, allocatable :: comp_var3d_nc(:,:,:,:)             !Netcdf additional compounds array (i,j,t,compound)
    real, allocatable :: comp_var4d_nc(:,:,:,:,:)           !Netcdf additional compounds array (i,j,k,t,compound)
    !Alternative meteorological input data netcdf files. No source element
    real, allocatable :: meteo_var1d_nc(:,:)
    real, allocatable :: meteo_var2d_nc(:,:,:)  
    real, allocatable :: meteo_var3d_nc(:,:,:,:)
    real, allocatable :: meteo_var4d_nc(:,:,:,:,:)
    
    real dgrid_nc(2)
    real meteo_dgrid_nc(2)
    double precision, allocatable :: val_dim_nc(:,:)
    double precision, allocatable :: val_dim_meteo_nc(:,:)
    character(256), allocatable :: unit_dim_nc(:)
    character(256), allocatable :: unit_dim_meteo_nc(:)
    
    integer surface_level_nc
    
    !Declare road link data arrays
    real, allocatable :: inputdata_rl(:,:)
    integer, allocatable :: inputdata_int_rl(:,:)
    
    !Road link (rl) indexes
    integer x1_rl_index,x2_rl_index,y1_rl_index,y2_rl_index,x0_rl_index,y0_rl_index
    parameter (x1_rl_index=1,x2_rl_index=2,y1_rl_index=3,y2_rl_index=4,x0_rl_index=5,y0_rl_index=6)
    integer lon1_rl_index,lon2_rl_index,lat1_rl_index,lat2_rl_index,lon0_rl_index,lat0_rl_index
    parameter (lon1_rl_index=7,lon2_rl_index=8,lat1_rl_index=9,lat2_rl_index=10,lon0_rl_index=11,lat0_rl_index=12)
    integer length_rl_index,z0_rl_index,width_rl_index
    parameter (length_rl_index=13,z0_rl_index=14,width_rl_index=15)
    integer adt_rl_index,speed_rl_index,hdv_rl_index
    parameter (adt_rl_index=16,speed_rl_index=17,hdv_rl_index=18)
    integer num_var_rl
    parameter(num_var_rl=18)

    integer id_rl_index,roadtype_rl_index,nlanes_rl_index
    parameter (id_rl_index=1,roadtype_rl_index=2,nlanes_rl_index=3)
    integer num_int_rl
    parameter(num_int_rl=3)

    !Declare file and path names for input roadlink files
    character(256) filename_rl(2)
    character(256) pathname_rl(2)
    character(256) pathfilename_rl(2)
    character(256) file_tag

    !Declare ais shiping data arrays
    !real, allocatable :: inputdata_ship(:,:)
    !integer, allocatable :: inputdata_int_ship(:,:)
    !character(8), allocatable :: inputdata_char_ship(:,:)
    
    !Input shipping (ship) indexes ddlatitude;ddlongitude;totalnoxemission;totalparticulatematteremission;fk_vessellloydstype;fk_ais_norwegianmainvesselcategory;date;time
    integer lon_ship_index,lat_ship_index,nox_ship_index,tpm_ship_index
    parameter (lon_ship_index=1,lat_ship_index=2,nox_ship_index=3,tpm_ship_index=4)
    integer num_var_ship
    parameter(num_var_ship=4)

    integer idloyds_ship_index,idnorwegian_ship_index
    parameter (idloyds_ship_index=1,idnorwegian_ship_index=2)
    integer num_int_ship
    parameter(num_int_ship=2)

    integer date_ship_index,time_ship_index
    parameter (date_ship_index=1,time_ship_index=2)
    integer num_char_ship
    parameter(num_char_rl=2)

    !Declare file and path names for shipping ais files
    character(256) filename_ship(2)
    character(256) pathname_ship(2)
    character(256) pathfilename_ship(2)
    
    !Input agriculture
    integer lon_agriculture_index,lat_agriculture_index,nh3_agriculture_index
    parameter (lon_agriculture_index=1,lat_agriculture_index=2,nh3_agriculture_index=3)
    integer num_var_agriculture
    parameter(num_var_agriculture=3)

    !Declare file and path names for input agriculture rivm files
    character(256) filename_agriculture(2)
    character(256) pathname_agriculture(2)
    character(256) pathfilename_agriculture(2)

    !Declare file and path names for SSB building and population files
    character(256) filename_heating(2)
    character(256) pathname_heating(2)
    character(256) pathfilename_heating(2)
    
    
    !integer traffic_emission_file_index,traffic_proxy_file_index,traffic_proxy_integral_file_index
    !parameter (traffic_emission_file_index=1,traffic_proxy_file_index=2,traffic_proxy_integral_file_index=3)
    
    !integer shipping_emission_file_index,shipping_proxy_file_index,shipping_proxy_integral_file_index
    !parameter (shipping_emission_file_index=4,shipping_proxy_file_index=5,shipping_proxy_integral_file_index=6)
    
    !integer emep_traffic_input_file_index,emep_traffic_subgrid_file_index,emep_traffic_subgrid_lf_file_index,traffic_subgrid_local_file_index
    !parameter(emep_traffic_input_file_index=7,emep_subgrid_file_index=8,emep_subgrid_frac_file_index=9,subgrid_local_file_index=10)
    
    !integer emep_shipping_input_file_index,emep_shipping_subgrid_file_index,emep_shipping_subgrid_lf_file_index,shipping_subgrid_local_file_index
    !parameter(emep_shipping_input_file_index=11,emep_shipping_subgrid_file_index=12,emep_shipping_subgrid_lf_file_index=13,shipping_subgrid_local_file_index=14)
    
    !integer subgrid_uwind_file_index,subgrid_vwind_file_index,subgrid_hmix_file_index,subgrid_kz_file_index
    !parameter (subgrid_uwind_file_index=15,subgrid_vwind_file_index=16,subgrid_hmix_file_index=17,subgrid_kz_file_index=18)
    
    !integer emep_agriculture_input_file_index,emep_agriculture_subgrid_file_index,emep_agriculture_subgrid_lf_file_index,agriculture_subgrid_local_file_index
    !parameter(emep_agriculture_input_file_index=19,emep_agriculture_subgrid_file_index=20,emep_agriculture_subgrid_lf_file_index=21,agriculture_subgrid_local_file_index=22)
 
    !integer agriculture_emission_file_index,agriculture_proxy_file_index,agriculture_proxy_integral_file_index
    !parameter (agriculture_emission_file_index=23,agriculture_proxy_file_index=24,agriculture_proxy_integral_file_index=25)

    !Definition of the EMEP file to be read. Two files. 1 is base file and 2 is uEMEP file
    character(256) filename_EMEP(4)
    character(256) pathname_EMEP(4)
    character(256) pathfilename_EMEP(4)  !Combined path and filename

    !Definition of the receptor file to be read.
    character(256) filename_receptor
    character(256) pathname_receptor
    character(256) pathfilename_receptor  !Combined path and filename

    !Definition of the timeprofilefile to be read.
    character(256) filename_timeprofile
    character(256) pathname_timeprofile
    character(256) pathfilename_timeprofile  !Combined path and filename

    !Declaration for all filenames used to save gridded data produced by uEMEP. Up to 200 file names but of course can be more
    integer n_filenames_grid
    parameter (n_filenames_grid=200)
    character(256) filename_grid(n_filenames_grid)
    character(256) pathname_grid(n_filenames_grid)
    character(256) pathfilename_grid(n_filenames_grid)  !Combined path and filename
    logical read_existing_grid_data(n_filenames_grid)
    character(256) pathname_output_grid
    
    logical :: save_intermediate_files=.false.
    
    !integer traffic_param_index,shipping_param_index,heating_param_index,agriculture_param_index
    !parameter (traffic_param_index=1,shipping_param_index=2,heating_param_index=3,agriculture_param_index=4)
    !integer num_param_index
    !parameter (num_param_index=4)
    
    !Declare subgrid variable indexes for 'subgrid' array
    integer proxy_subgrid_index,proxy_integral_subgrid_index
    integer scaling_factor_subgrid_index,local_subgrid_index,nonlocal_subgrid_index,total_subgrid_index
    integer emep_subgrid_index,emep_frac_subgrid_index,emep_local_subgrid_index,emep_nonlocal_subgrid_index,proxy_average_integral_subgrid_index
    parameter (proxy_subgrid_index=1,proxy_integral_subgrid_index=2)
    parameter (scaling_factor_subgrid_index=3,local_subgrid_index=4,nonlocal_subgrid_index=5,total_subgrid_index=6)
    parameter (emep_subgrid_index=7,emep_frac_subgrid_index=8,emep_local_subgrid_index=9,emep_nonlocal_subgrid_index=10,proxy_average_integral_subgrid_index=11)
    integer n_subgrid_index
    parameter (n_subgrid_index=11)

    !Declare meteo subgrid variables. Must be the same as the nc version
    integer ugrid_subgrid_index,vgrid_subgrid_index,FF10_subgrid_index,FFgrid_subgrid_index,inv_FFgrid_subgrid_index,inv_FF10_subgrid_index
    integer hmix_subgrid_index,kz_subgrid_index,invL_subgrid_index,ustar_subgrid_index,logz0_subgrid_index,J_subgrid_index,cos_subgrid_index,sin_subgrid_index
    parameter (ugrid_subgrid_index=1,vgrid_subgrid_index=2,FF10_subgrid_index=3,FFgrid_subgrid_index=4,inv_FFgrid_subgrid_index=5,inv_FF10_subgrid_index=6)
    parameter (hmix_subgrid_index=7,kz_subgrid_index=8,invL_subgrid_index=9,ustar_subgrid_index=10,logz0_subgrid_index=11,J_subgrid_index=12,cos_subgrid_index=13,sin_subgrid_index=14)
    integer n_meteo_subgrid_index
    parameter (n_meteo_subgrid_index=14)

    !Declare compund indexes
    integer no2_index,nox_index,pm25_index,pm10_index,nh3_index,o3_index,no_index,traveltime_index
    parameter (no2_index=1,nox_index=2,pm25_index=3,pm10_index=4,nh3_index=5,o3_index=6,no_index=7,traveltime_index=8)
    !Declare source indexes (type_source)
    integer allsource_index,traffic_index,shipping_index,heating_index,agriculture_index,industry_index
    parameter (allsource_index=1,traffic_index=2,shipping_index=3,heating_index=4,agriculture_index=5,industry_index=6)
    integer n_compound_index,n_source_index
    parameter (n_compound_index=8,n_source_index=6)
    integer compound_source_index(n_compound_index,n_source_index)
    
    character(256) source_file_postfix(n_source_index)
    logical calculate_source(n_source_index)
    logical make_EMEP_grid_emission_data(n_source_index)
    logical replace_EMEP_local_with_subgrid_local(n_source_index)
    logical combine_emission_subsources_during_dispersion(n_source_index)

    integer x_dim_index,y_dim_index,t_dim_index,n_dim_index
    parameter (x_dim_index=1,y_dim_index=2,t_dim_index=3,n_dim_index=3)
    
    !Target redistribution grid, the same for all sources and compounds
    !subgrid (i,j,t,type_subgrid,type_source,type_subsource)
    integer subgrid_dim(n_dim_index)
    real subgrid_delta(2),subgrid_min(2),subgrid_max(2)   !Only x and y
    integer init_subgrid_dim(n_dim_index)
    real init_subgrid_delta(2),init_subgrid_min(2),init_subgrid_max(2)   !Only x and y
    real, allocatable :: subgrid(:,:,:,:,:,:)
    real, allocatable :: meteo_subgrid(:,:,:,:)
    real, allocatable :: last_meteo_subgrid(:,:,:)
    real, allocatable :: comp_subgrid(:,:,:,:)
    real, allocatable :: comp_EMEP_subgrid(:,:,:,:)
    real, allocatable :: orig_EMEP_subgrid(:,:,:,:)
    real, allocatable :: x_subgrid(:,:)
    real, allocatable :: y_subgrid(:,:)
    real, allocatable :: lon_subgrid(:,:)
    real, allocatable :: lat_subgrid(:,:)
    real, allocatable :: xproj_subgrid(:,:)
    real, allocatable :: yproj_subgrid(:,:)
    real, allocatable :: traveltime_subgrid(:,:,:,:)
    real, allocatable :: exposure_subgrid(:,:,:,:)
    integer subgrid_loop_index(2)               !Number of target subgrids to loop through, limitted by the size of the EMEP grid
    integer integral_subgrid_loop_index(2)      !Number of integral subgrids to loop through, limitted by the size of the EMEP grid
    integer emission_subgrid_loop_index(2,n_source_index)      !Number of emission subgrids to loop through, limitted by the size of the EMEP grid
    logical, allocatable :: use_subgrid(:,:,:)    !Specifies if calculations are to be made at a particular set of target subgrids or not
    real :: loop_index_scale=1.5
    real :: buffer_index_scale=1.5
    
    !Emission subgrid per source and subsource. Can be time dependent
    !Each source type has its own x and y and dim
    !Each source may be of lesser dimmensions than the total array size (which is the same as the target grid)
    integer n_possible_subsource
    parameter (n_possible_subsource=2)
    integer :: n_subsource(n_source_index)=1 !Initialise the number of actual emission subsources to 1 for all subsources
    character(2) subsource_str(n_possible_subsource)
    
    integer emission_h_index,emission,emission_sigz00_index,emission_sigy00_index,emission_minFF_index,n_emission_index
    parameter (emission_h_index=1,emission_sigz00_index=2,emission_sigy00_index=3,emission_minFF_index=4,n_emission_index=4)
    
    integer emission_subgrid_dim(n_dim_index,n_source_index)
    integer emission_max_subgrid_dim(n_dim_index)
    real emission_subgrid_delta(2,n_source_index),emission_subgrid_min(2,n_source_index),emission_subgrid_max(2,n_source_index)  !Only x and y
    !emission_subgrid (i,j,t,n_source,n_subsource)
    real, allocatable :: emission_subgrid(:,:,:,:,:)
    real, allocatable :: proxy_emission_subgrid(:,:,:,:) !No time dependence
    real, allocatable :: emission_time_profile_subgrid(:,:,:,:,:)
    real, allocatable :: emission_properties_subgrid(:,:,:,:,:) !No time dependence
    !x_emission_subgrid (i,j,n_source)
    real, allocatable :: x_emission_subgrid(:,:,:)
    real, allocatable :: y_emission_subgrid(:,:,:)
    real, allocatable :: lon_emission_subgrid(:,:,:)
    real, allocatable :: lat_emission_subgrid(:,:,:)
    real, allocatable :: xproj_emission_subgrid(:,:,:)
    real, allocatable :: yproj_emission_subgrid(:,:,:)
    
    logical :: use_buffer_zone=.true.
    integer buffer_index(2)
    real    buffer_size(2)
    integer emission_buffer_index(2,n_source_index)
    real    emission_buffer_size(2,n_source_index)
    integer integral_buffer_index(2)
    real    integral_buffer_size(2)

    integer integral_subgrid_dim(n_dim_index)
    real integral_subgrid_delta(2),integral_subgrid_min(2),integral_subgrid_max(2)  !Only x and y
    !emission_subgrid (i,j,t,n_source,n_subsource)
    real, allocatable :: integral_subgrid(:,:,:,:,:)
    !x_emission_subgrid (i,j,n_source)
    real, allocatable :: x_integral_subgrid(:,:)
    real, allocatable :: y_integral_subgrid(:,:)
    real, allocatable :: lon_integral_subgrid(:,:)
    real, allocatable :: lat_integral_subgrid(:,:)
    real, allocatable :: xproj_integral_subgrid(:,:)
    real, allocatable :: yproj_integral_subgrid(:,:)
    real, allocatable :: meteo_nc_xproj_integral_subgrid(:,:)
    real, allocatable :: meteo_nc_yproj_integral_subgrid(:,:)

    integer population_subgrid_dim(2)
    real population_subgrid_delta(2),population_subgrid_min(2),population_subgrid_max(2)  !Only x and y
    real, allocatable :: population_subgrid(:,:,:)
    real, allocatable :: x_population_subgrid(:,:)
    real, allocatable :: y_population_subgrid(:,:)
    real, allocatable :: lon_population_subgrid(:,:)
    real, allocatable :: lat_population_subgrid(:,:)
    real, allocatable :: xproj_population_subgrid(:,:)
    real, allocatable :: yproj_population_subgrid(:,:)

    !Initial emission grid, can be different for different sources so specified individually
    !emission_subgrid(i,j,value_type)
    !integer x_emission_subgrid_index,y_emission_subgrid_index,lon_emission_subgrid_index,lat_emission_subgrid_index
    !parameter (x_emission_subgrid_index=1,y_emission_subgrid_index=2,lon_emission_subgrid_index=3,lat_emission_subgrid_index=4)
    !integer n_possible_emission_subgrids
    !parameter (n_possible_emission_subgrids=5)
    !integer proxy_emission_subgrid_index(n_possible_emission_subgrids) !5 possible emission subgrids for each emission source
    !integer :: n_emission(n_source_index)=1 !Initialise the number of possible emission subgrids to 1 for all subsources
    !integer n_emission_subgrid_index
    !parameter (n_emission_subgrid_index=9) !4+5
    
    !integer traffic_emission_subgrid_dim(2)
    !real traffic_emission_subgrid_delta(2),traffic_emission_subgrid_min(2),traffic_emission_subgrid_max(2)
    !real, allocatable :: traffic_emission_subgrid(:,:,:)
    !integer shipping_emission_subgrid_dim(2)
    !real shipping_emission_subgrid_delta(2),shipping_emission_subgrid_min(2),shipping_emission_subgrid_max(2)
    !real, allocatable :: shipping_emission_subgrid(:,:,:)
    !integer agriculture_emission_subgrid_dim(2)
    !real agriculture_emission_subgrid_delta(2),agriculture_emission_subgrid_min(2),agriculture_emission_subgrid_max(2)
    !real, allocatable :: agriculture_emission_subgrid(:,:,:)
    !integer heating_emission_subgrid_dim(2)
    !real heating_emission_subgrid_delta(2),heating_emission_subgrid_min(2),heating_emission_subgrid_max(2)
    !real, allocatable :: heating_emission_subgrid(:,:,:)
    
    !
    integer, allocatable :: crossreference_target_to_emep_subgrid(:,:,:) !(i,j,dim)
    integer, allocatable :: crossreference_integral_to_emep_subgrid(:,:,:) !(i,j,dim)
    integer, allocatable :: crossreference_integral_to_meteo_nc_subgrid(:,:,:) !(i,j,dim)
    integer, allocatable :: crossreference_target_to_integral_subgrid(:,:,:) !(i,j,dim)
    !integer, allocatable :: crossreference_emission_to_target_subgrid(:,:,:,:) !(i,j,dim,n_source)
    integer, allocatable :: crossreference_target_to_emission_subgrid(:,:,:,:) !(i,j,dim,n_source)
    integer, allocatable :: crossreference_integral_to_emission_subgrid(:,:,:,:) !(i,j,dim,n_source)
    integer, allocatable :: crossreference_emission_to_emep_subgrid(:,:,:,:) !(i,j,dim,n_source)
    integer, allocatable :: crossreference_emission_to_integral_subgrid(:,:,:,:) !(i,j,dim,n_source)
    integer, allocatable :: crossreference_target_to_population_subgrid(:,:,:) !(i,j,dim)
    !integer, allocatable :: crossreference_shipping_emission_subgrid(:,:,:)
    !integer, allocatable :: crossreference_agriculture_emission_subgrid(:,:,:)
    !integer, allocatable :: crossreference_heating_emission_subgrid(:,:,:)
    
    real :: min_link_size=50.
    real :: min_adt=1000.
    real :: H_emep=90. !Height of lowest level in EMEP
    real :: H_meteo=45. !Height of the gridded meteo values 
    
    logical :: hourly_calculations=.false.
    logical :: annual_calculations=.false.
    
    !Pseudo dispersion parameters
    real z_rec(n_source_index,n_possible_subsource)
    real ay(n_source_index,n_possible_subsource),by(n_source_index,n_possible_subsource),az(n_source_index,n_possible_subsource),bz(n_source_index,n_possible_subsource)
    real sig_y_0(n_source_index,n_possible_subsource),sig_z_0(n_source_index,n_possible_subsource)
    !To be set at input
    real sig_y_00(n_source_index,n_possible_subsource),sig_z_00(n_source_index,n_possible_subsource),h_emis(n_source_index,n_possible_subsource)
    integer stability_scheme_flag
    
    real :: FF_min_dispersion=0.1
    real :: emission_timeprofile_hour_shift=0.
    
    real pi
    parameter (pi=3.141592)

    !Near neighbour arrays
    !character(1) NNi_str(3),NNj_str(3)
    !character(256) NNpre_str
    !character(256) NNsource_str(n_source_index)
    !character(256) NN_name_nc(n_source_index,3,3)
    !integer NN_index(3,3)
    !integer num_NN_nc
    !parameter (num_NN_nc=9)                  ! number of readable near neighbours

    integer UTM_projection_index,RDM_projection_index,LCC_projection_index,LL_projection_index
    parameter (UTM_projection_index=1,RDM_projection_index=2,LCC_projection_index=3,LL_projection_index=4)
    integer :: projection_type=UTM_projection_index
    integer :: EMEP_projection_type=LL_projection_index
    double precision :: EMEP_projection_attributes(10)
    double precision :: projection_attributes(10)
    integer :: meteo_nc_projection_type=LCC_projection_index
    double precision :: meteo_nc_projection_attributes(10)
    logical :: use_alternative_LCC_projection_flag=.false.

    
    !Filename index for files produced by uEMEP
    integer proxy_emission_file_index(n_source_index),emission_file_index(n_source_index),proxy_file_index(n_source_index),proxy_integral_file_index(n_source_index)
    integer emep_subgrid_file_index(n_source_index),emep_subgrid_nonlocal_file_index(n_source_index),emep_subgrid_local_file_index(n_source_index),emep_subgrid_frac_file_index(n_source_index)
    integer subgrid_local_file_index(n_source_index),subgrid_total_file_index(n_source_index)
    !Filename index for meteorological parameters
    integer subgrid_ugrid_file_index,subgrid_vgrid_file_index,subgrid_hmix_file_index,subgrid_kz_file_index,subgrid_logz0_file_index,subgrid_invL_file_index,subgrid_FF10_file_index,subgrid_FFgrid_file_index
    integer subgrid_invFF10_file_index,subgrid_invFFgrid_file_index,subgrid_ustar_file_index,subgrid_J_file_index,subgrid_meteo_file_index
    !Filename index for grid auto grid parameters
    integer use_subgrid_file_index(n_source_index)
    integer emep_emission_subgrid_file_index(n_source_index)
    
    character(256) source_file_str(n_source_index)
    real :: unit_conversion(n_source_index)=1.
    
    real :: emission_factor_conversion(n_compound_nc_index,n_source_index,n_possible_subsource)=0.
    real :: emission_factor(n_compound_nc_index,n_source_index,n_possible_subsource)=1.
    real :: ratio_truck_car_emission(n_compound_nc_index)=10.
    
    integer :: integral_subgrid_step=1
    integer :: weighting_step=1

    real :: limit_shipping_delta=250.
    real :: limit_heating_delta=250.
    real :: limit_population_delta=250.
    real :: traj_step_scale=2.
    
    integer :: start_time_nc_index=1
    integer :: end_time_nc_index=1
    integer :: start_time_meteo_nc_index=1
    integer :: end_time_meteo_nc_index=1
    integer :: ref_year_EMEP=1900
    integer :: ref_year_meteo=1970

    integer :: EMEP_emission_grid_interpolation_flag=0
    logical :: subgrid_emission_distribution_flag=.false.  !If true then distributes the EMEP emissions to the existing emission subgrid
    logical :: use_downwind_position_flag=.false.           !If true then searches the upwind EMEP grid position for emissions
    logical :: use_trajectory_flag=.false.
    integer :: no2_chemistry_scheme_flag=1
    integer :: wind_level_flag=3
    integer :: wind_level_integral_flag=3
    logical :: use_receptor_positions_for_auto_subgrid_flag=.false.
    logical :: interpolate_subgrids_flag=.false.
    logical :: use_aggregated_shipping_emissions_flag=.true.
    logical :: calculate_aggregated_shipping_emissions_flag=.false.
    logical :: average_zc_h_in_Kz_flag=.true.

    
    integer n_receptor,n_receptor_in,n_receptor_max,n_valid_receptor,n_valid_receptor_in
    parameter (n_receptor_max=1000)
    real lon_receptor(n_receptor_max),lat_receptor(n_receptor_max),x_receptor(n_receptor_max),y_receptor(n_receptor_max)
    real lon_receptor_in(n_receptor_max),lat_receptor_in(n_receptor_max),x_receptor_in(n_receptor_max),y_receptor_in(n_receptor_max)
    integer i_receptor_subgrid(n_receptor_max),j_receptor_subgrid(n_receptor_max)
    character(256) name_receptor(n_receptor_max,2)
    character(256) name_receptor_in(n_receptor_max,2)
    logical :: use_receptor(n_receptor_max)=.true.
    integer :: use_receptor_region=1 !Surrounding grid region when just running for receptors
    integer valid_receptor_index(n_receptor_max)
    integer valid_receptor_inverse_index(n_receptor_max)
    
    !Indicies for SSB building and population data
    integer dwelling_index,population_index,establishment_index,school_index,home_index,n_population_index
    parameter(dwelling_index=1,population_index=2,establishment_index=3,school_index=4,kindergaten_index=5,home_index=6,n_population_index=6)
    integer population_file_index(n_population_index)
    
    character(256) filename_population(n_population_index),pathname_population(n_population_index),pathfilename_population(n_population_index)
    
    !Is a temporary variable set when reading SSB data
    integer :: SSB_data_type=dwelling_index
    
    logical :: calculate_population_exposure_flag=.false.
    logical :: use_population_positions_for_auto_subgrid_flag=.false.
    integer :: population_data_type=population_index
    
    logical :: use_multiple_receptor_grids_flag=.false.
    integer :: n_receptor_grid
    integer start_grid_loop_index,end_grid_loop_index
    integer g_loop
    
    logical :: use_last_meteo_in_dispersion=.false.
    logical :: use_meandering_in_dispersion=.false.
    logical :: use_traffic_for_sigma0_flag=.false.
!    logical :: use_traffic_for_minFF_flag=.false.
    logical :: use_emission_grid_gradient_flag=.false.
    logical :: use_alternative_meteorology_flag=.false.
    logical :: use_alternative_z0_flag=.false.
    logical :: 	save_netcdf_file_flag=.false.
    logical ::	save_netcdf_receptor_flag=.false.
    
    !Some testing and scaling values
    real :: replace_z0=NODATA_value  !Will not replace z0 when it has a NODATA value
    real :: replace_invL=NODATA_value  !Will not replace invL when it has a NODATA value
    real :: replace_hmix=NODATA_value  !Will not replace mix when it has a NODATA value
	real :: FF_scale=NODATA_value
	real :: FF10_offset=NODATA_value
	real :: DD_offset=NODATA_value
    
    real :: hmix_max=2000.
    real :: hmix_min=25.
    real :: ustar_min=0.001
    
    end module uEMEP_definitions
    
    
