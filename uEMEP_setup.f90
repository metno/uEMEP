!   Set up of uEMEP before starting calculations
    subroutine uEMEP_setup
    
    use uEMEP_definitions

    implicit none
    
    real read_name_real
    logical read_name_logical
    integer read_name_integer
    character(256) read_name_char
    
    integer :: unit_in=30
    
!==========================================================================
!   uEMEP model setup
!==========================================================================

    write(unit_logfile,'(A)') ''
    write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Setting model run characteristics (uEMEP_setup)'
	write(unit_logfile,'(A)') '================================================================'

    
    !Choose multiple hourly or single average calculations. Will be changed based on the time dimensions in EMEP
    hourly_calculations=.true.
    annual_calculations=.false.

    !Compound name
    input_comp_name='nox'
    
    !Set the compound that is to be redistributed and the local fraction compound that is to be used
    compound_index=nox_nc_index
    compound_frac_index=nox_nc_index
    !compound_index=nh3_nc_index
    !compound_frac_index=nh3_nc_index
    compound_index=pm25_nc_index
    compound_frac_index=pm25_nc_index
   
    !Specify which files to read from existing calculations
    read_existing_grid_data=.false.
    
    read_existing_grid_data(proxy_emission_file_index(traffic_index))=.false.
    read_existing_grid_data(proxy_emission_file_index(shipping_index))=.false.
    read_existing_grid_data(proxy_emission_file_index(agriculture_index))=.false.
    read_existing_grid_data(proxy_emission_file_index(heating_index))=.false.
    
    read_existing_grid_data(proxy_file_index(traffic_index))=.false.
    read_existing_grid_data(proxy_file_index(shipping_index))=.false.
    read_existing_grid_data(proxy_file_index(agriculture_index))=.false.
    
    read_existing_grid_data(emep_subgrid_file_index)=.false.
    !read_existing_grid_data(subgrid_local_file_index)=.false. !Not used

    read_existing_grid_data(subgrid_meteo_file_index)=.false.
    
    read_existing_grid_data(use_subgrid_file_index(traffic_index))=.false.
    read_existing_grid_data(use_subgrid_file_index(shipping_index))=.false.
    read_existing_grid_data(use_subgrid_file_index(agriculture_index))=.false.

    !Specify which source to calculate. Overridden in the set_subgrids routine, usually.
    calculate_source=.false.
    calculate_source(traffic_index)=.false.
    calculate_source(shipping_index)=.false. 
    calculate_source(agriculture_index)=.false.
    calculate_source(heating_index)=.true.
    
    !Specify if the proxy data is to be made into EMEP grids for cross checking
    !Will create both an EMEP and a proxy grid file per EMEP grid
    make_EMEP_grid_emission_data=.false.
    make_EMEP_grid_emission_data=.true.
    replace_EMEP_local_with_subgrid_local=.false.
    !replace_EMEP_local_with_subgrid_local(agriculture_index)=.true.
    
    !Set the utm zone being used and it's central position
    utm_zone=33;utm_lon0=15
    !utm_zone=32;utm_lon0=9
    
    file_tag='Oslo_region_100m'
    !file_tag='Oslo_centre_25m'
    file_tag='Oslo_250m'
    !file_tag='Oslo_100m'
    !file_tag='Oslo_50m'
    !file_tag='Oslo_25m'
    !file_tag='Bergen_region_250m'
    !file_tag='Bergen_50m'
    !file_tag='Stavanger_region_250m'
    !file_tag='Stavanger_50m'
    !file_tag='Stavanger_50m'
    !file_tag='Oslo_50m'
    !file_tag='Oslo_25m'
    !file_tag='Bergen_25m'
    !file_tag='Oslo_region_250m'
    !file_tag='Southern_250m'
    !file_tag='SouthEast_250m'
    !file_tag='Narvik_50m'
    !file_tag='Nederland_500m'
    !file_tag='Nederland_250m'
    !file_tag='Nederland_100m'
    !file_tag='Nederland_50km_500m'
    
    !Choose EMEP grid interpolation method 0=none, 1= area weighted, 2=emission weighted, 3=aggregated to emission integral weighted, 4=proxy integral weighted
    EMEP_grid_interpolation_flag=3
    
    !Choose EMEP emission grid interpolation method 0=none, 1=area weighted, 2=emission weighted(recommended)
    EMEP_emission_grid_interpolation_flag=2
    !If true then distributes the EMEP grid emissions to the existing emission subgrid using the proxy emission data
    !If false then distributes the EMEP grid emissions to all subgrids. Should be true when using proxy distributions
    subgrid_emission_distribution_flag=.true.
 
    !Choose the region for interpolation 1 is one EMEP grid 2 is 2 EMEP grids. Triple check of this method required
    EMEP_grid_interpolation_size=1.0
    
    !Use downwind positioning of subgrids. Only works properly when area weighting for EMEP_emission_grid_interpolation_flag is chosen
    use_downwind_position_flag=.true.
    
    !Choose the way the local contribution is calculated.
    !1 is redistribution of EMEP concentrations
    !2 is dispersion based on prescribed emission factors and time profiles (if available)
    !3 is dispersion based on redistribution of EMEP emissions
    !4 is dispersion based on prescribed emission factors but using EMEP time profiles
    local_subgrid_method_flag=4
    
    !Choose dispersion parameters and scheme
    !1 Neutral scheme
    !2 PG scheme using Klug parameters 
    stability_scheme_flag=2
    
    !Wind flag
    !1 Use gridded wind speed without interpolation
    !2 Use gridded wind speed with log interpolation to emission heights
    !3 Use 10 m wind speed without interpolation
    !4 Use 10 m wind speed with log interpolation to emission heights
    wind_level_flag=4
    wind_level_integral_flag=4
     
    !Choose no2 chemistry scheme
    !0 no chemistry
    !1 photostationary
    !2 time scale photo chemistry
    no2_chemistry_scheme_flag=2
    
    !Choose the use of auto subgridding
    use_emission_positions_for_auto_subgrid_flag=.false.
    use_receptor_positions_for_auto_subgrid_flag=.false.
    use_population_positions_for_auto_subgrid_flag=.true.

    interpolate_subgrids_flag=.false.
    
    !Use trajectory for disperion calculations.
    !traj_step_scale indicates the step size of the trajectories in regard to the meteorological subgrid.
    use_trajectory_flag=.true.
    traj_step_scale=2.
    
    !Set shipping emission aggregation flags
    use_aggregated_shipping_emissions_flag=.true.
    calculate_aggregated_shipping_emissions_flag=.false.
    
    !Calculate population exposure
    calculate_population_exposure_flag  =  .true.

    !Choose the time in the input file
    start_time_nc_index=1
    end_time_nc_index=8+start_time_nc_index-1
    
    !Set constants, indexes and names based on choice of annual or hourly
    call uEMEP_set_constants
    
    call uEMEP_set_filenames
    
    !call uEMEP_set_flags
    
    !call uEMEP_set_params
    
    call uEMEP_set_subgrids
    
    call uEMEP_set_emission_factors
    
    !Read this routine as a name file
    name_config_file='C:\uEMEP\Fortran\code\uEMEP_setup.f90'
    open(unit_in,file=name_config_file,access='sequential',status='old',readonly) 
    
    traj_step_scale=                                                    read_name_real('traj_step_scale',traj_step_scale,unit_in,unit_logfile)
    
    calculate_population_exposure_flag=                                 read_name_logical('calculate_population_exposure_flag',.false.,unit_in,unit_logfile)
    read_existing_grid_data(proxy_emission_file_index(traffic_index))=  read_name_logical('read_existing_grid_data(proxy_emission_file_index(traffic_index))',.false.,unit_in,unit_logfile)
    
    EMEP_grid_interpolation_flag=                                       read_name_integer('EMEP_grid_interpolation_flag',0,unit_in,unit_logfile)
    
    file_tag=                                                           read_name_char('file_tag','',unit_in,unit_logfile)
    
    !Find the correct compound index
    input_comp_name=                                                    read_name_char('input_comp_name','',unit_in,unit_logfile)

    do i=1,n_compound_nc_index
        if (trim(var_name_nc(conc_nc_index,i,allsource_index)).eq.trim(input_comp_name)) then
            compound_index=i
        endif
        if (trim(var_name_nc(frac_nc_index,i,allsource_index)).eq.trim(input_comp_name)) then
            compound_frac_index=i
        endif
    enddo

    close (unit_in)
    stop
    
    end subroutine uEMEP_setup
    