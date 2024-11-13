---
title: Tutorial
---

**Contents:**

- [Introduction](#introduction)
- [Getting started](#getting-started)
    - [Installing uEMEP](#installing-uemep)
    - [The configuration file](#the-configuration-file)
    - [Initial configuration](#initial-configuration)
- [Setting the spatial extent](#setting-the-spatial-extent)
- [EMEP forcing data](#emep-forcing-data)
- [Emission redistribution proxy data](#emission-redistribution-proxy-data)
    - [Traffic emissions](#traffic-emissions)
    - [Residential heating emissions](#residential-heating-emissions)
    - [Shipping emissions](#shipping-emissions)
- [uEMEP processes](#uemep-processes)
    - [Wind and dispersion](#wind-and-dispersion)
    - [NO2 chemistry](#no2-chemistry)
    - [Tunnels](#tunnels)
- [Finalize configuration](#finalize-configuration)
- [Running uEMEP](#running-uemep)
    - [Multiple configuration files](#multiple-configuration-files)
    - [Running in parallel](#running-in-parallel)
- [References](#references)


## Introduction <a name="introduction"></a>

In general, uEMEP can be configured to run using two different downscaling methods. Both of these methods use Gaussian dispersion modeling to recalculate air pollutant concentrations in high resolution. The major difference between the methods is the origin and handling of emission sources:

1. Emission redistribution of EMEP gridded emissions using proxy data
2. Independent emissions data

The Norwegian forecast system uses independent emission sources to deliver daily forcasts of air quality, covering the entire country in high resolution. However, independent emission data is not readily available across Europe, and instead redistribution of EMEP gridded emissions are used (Mu et al., 2022). 

For the scope of this tutorial, instructions will thus focus on method 2. Proxy data for the redistribution of emissions are based on openly available data sources and the tutorial allows uEMEP to be configured for downscaling pollutant anywhere within the EMEP domain. 

Disclaimer: The data and configuration provided in this tutorial is not necessarily the only or best way to configure uEMEP or downscale air quality in Europe in general. The purpose of this tutorial is to demonstrate the process of configuring uEMEP.

## Getting started <a name="getting-started"></a>

To follow this tutorial, it is assumed that the user has access to a relatively new Linux distribution with the necessary software installed or the rights to install the software. This guide has been tested using a virtual machine running Ubuntu 24.04 with 4 CPUs, 4GB memory and 40GB storage.

The uEMEP model requires the output from an existing EMEP simulation with calculations of the local fractions. In addition, suitable redistribution proxy data is required. How to generate and/or acquire these data sources from scratch is beyond the scope of this tutorial. 

The data for following this tutorial is available at: 

The download is 5.1GB which will expand to 12.5 when extracted. Download and extract the data as:

```bash
# Update and install necessary software
sudo apt update && sudo apt upgrade -y
sudo apt install wget unzip

# Download and extract the data
cd {local_dir}
unzip uemep_demo.zip
```

### Installing uEMEP <a name="installing-uemep"></a>

uEMEP is programmed using Fortran, and a compatible Fortran compiler is therefore necessary. In addition, as uEMEP reads and writes NetCDF files, a compatible NetCDF library should be installed. Currently, uEMEP supports the Intel (`ifort`) and GNU (`gfortran`) Fortran compilers.

In this tutorial, we will use `gfortran` and the NetCDF libraries from the default Ubuntu repository as this is the easiest to setup. Note however that uEMEP is developed with Intel `ifort` as the default compiler and this setup is better tested. In Ubuntu, install the necessary dependencies as:

```bash
sudo apt install build-essential git cmake gfortran libnetcdf-dev libnetcdff-dev
```

The uEMEP model code is open source, distributed under the LGPL-3.0 license. The latest release can be downloaded at github.com/metno/uEMEP/releases. Alternatively, the current stable version can be downloaded using git as:

```bash
cd {local_dir}/uemep_demo/
git clone https://github.com/metno/uEMEP.git
```

Since `ifort` is the default compiler, we need to change the compiler to `gfortran`:

```bash
# Open the CMakeLists.txt file
nano {local-dir}/uemep_demo/uEMEP/CMakeLists.txt

# Change the compiler to gfortran by modifying lines 10-11:
#set(CMAKE_Fortran_COMPILER "ifort")
set(CMAKE_Fortran_COMPILER "gfortran")
```

Then proceed to compile uEMEP as:

```bash
cd uEMEP
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make && make test
```

Verify that uEMEP is running, which should return the following information (the printed version may differ):

```bash
$ ./uemep --version
 uEMEP: Air quality dispersion model for high resolution downscaling of EMEP MSC-W
 
 Version: 7.0.1
 Copyright (C) 2007 Free Software Foundation.
 License GNU LGPL-3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>.
 This is free software: you are free to change and redistribute it.
 
 Developed and maintained at the Norwegian Meteorological Institute.
 Contribute at: <https://github.com/metno/uEMEP>
```

Lastly, symlink or copy the executable to the bin folder:

```bash
cd {local_dir}/uemep_demo/uEMEP/build
cp uemep ../../bin/
```

### The configuration file <a name="the-configuration-file"></a>

The configuration of uEMEP is flexible and thus has a large number of configurable parameters. While most of these parameters have default values, uEMEP cannot produce any useable output without additional user input. Which parameters that are necessary to define depends on the intended configuration and cannot easily be generalized.

In this tutorial, we will begin with an empty configuration file to demonstrate the various options requiring user input. A complete configuration file is provided with the downloaded data for reference. Create a new empty file as:

```bash
cd {local_dir}/uemep_demo/config_files
touch uemep_demo_config.txt
```

In the following sections, we will edit this file so uEMEP can run for any area in Europe.

### Initial configuration <a name="initial-configuration"></a>

uEMEP can run in two different time modes: hourly or annual. Since the EMEP output are annual averages, we set the uEMEP model to calculate annual average concentrations as well. In addition, we have to define which pollutants (i.e., compounds) we want to include:

```fortran
!---- Initial configuration ----!

! Set to calculate annual mean concentrations
annual_concentrations = .true.

! Set which compounds are calculated (gp_totals include NOx, NO2, PM25 and PM10)
input_comp_name = "gp_totals"
```

Next we specify some variables relating to how uEMEP handles path, files, and variable names.

```fortran
! Set the filename tag to be included in the output files
file_tag = "uemep_demo"

! Set the hour string to be included in the output files
forecast_hour_str = "00"

! Set some replacement strings which are used internally while reading and/or writing
replacement_date_str = "<replace_date>"
replacement_yesterday_date_str = "[replace_yesterday_date]"
replacement_hour_str = "<replace_hour>"

! Set date string for output files
filename_date_output_grid = "<yyyymmdd>_<replace_hour>"
```

Lastly we need to specify where uEMEP should write output, what kind of output it should produce:

```fortran
! Force uEMEP to write log output to the terminal
filename_log_file = ""

! Set the output directory (all output will go here)
pathname_output_grid = "{local_dir}/uemep_demo/model_output/"

! Set to save output as NetCDF (setting this to false will produce no output)
save_netcdf_file_flag = .true.

! Set to save downscaled and EMEP original pollutant concentrations
save_compounds = .true.
save_emep_original = .true.

! Set to save source contributions
save_source_contributions = .false.
save_no2_source_contributions = .false.
save_o3_source_contributions = .false.

! Set to turn off Norwegian air quality indicator output
save_aqi = .false.

! Set to save EMEP species (currently required to be set to true)
save_emep_species = .true.
```

## Setting the spatial extent <a name="setting-the-spatial-extent"></a>

uEMEP can run in two different spatial modes: tiles or stations. When running tiles, the output will be gridded concentrations, whereas when running stations, the output will be the concentration at specific geographical locations. In both cases, some information on the spatial setup is necessary. 

Here we will downscale tiles, so we set the size of the tile (area) as well as the resolution within the tile.

In this tutorial, the coordinates are defined to cover an area of 50x50 km in the Copenhagen-Malmö general area and with a horizontal resolution of 250 m. Coordinates uses the Lambert Azimuthal Equal Area (LAEA) projection [](https://epsg.io/3035) and have the unit meters. As an approximation, LAEA coordinates for Europe can be extracted from here: [](https://epsg.io/map#srs=3035-1149). Click on "Reproject Map" in the upper right corner. Clicking anywhere in this map will provide the LAEA coordinates. For example, for modeling the city of Paris the lower left coordinate would be 

![Paris example](../media/epsg_laea_paris_example.png "Paris example")

Remember however that increasing the tile size will increase the execution time and memory use (see [Running in parallel](#running-in-parallel) for modelling larger regions). 

```fortran
! Set the area to model
subgrid_min(x_dim_index) = 4465000.0
subgrid_min(y_dim_index) = 3595000.0
subgrid_max(x_dim_index) = 4515000.0
subgrid_max(y_dim_index) = 3645000.0

! Set the target resolution (in meters)
subgrid_delta(x_dim_index) = 250.0
subgrid_delta(y_dim_index) = 250.0
```

In addition, we provide information on the geographical projection that should be used. As mentioned above we use the Lambert Azimuthal Equal Area (LAEA) projection which is very useful because it is easy to define rectangular tiles which line up when a larger area is split into multiple smaller areas (tiles):

```fortran
! Set the projection type and attributes (5 = LAEA)
projection_type = 5
projection_attributes(1) = 10.0
projection_attributes(2) = 52.0
projection_attributes(3) = 4321000.0
projection_attributes(4) = 3210000.0
projection_attributes(5) = 6378137.0
projection_attributes(6) = 298.2572221
```

## EMEP forcing data <a name="emep-forcing-data"></a>

First set the path and filenames, and provide some general information about the EMEP data:

```fortran
!---- EMEP forcing data ----!

! Set path and filenames for the EMEP output
pathname_EMEP(1) = "{local_dir}/uemep_demo/model_input/emep/"
filename_EMEP(1) = "Base_fullrun.nc"
pathname_EMEP(2) = "{local_dir}/uemep_demo/model_input/emep/"
filename_EMEP(2) = "Base_LF_full.nc"

! Set the EMEP projection type
EMEP_projection_type = 4

! Set the EMEP aggregation period in hours (1 year = 8760 hours)
EMEP_emissions_aggregation_period = 8760

! Set the EMEP naming template to use the latest version
use_emission_naming_template_flag = .true.
emission_naming_template_str = "GNFR_<n>_Emis_mgm2_"

! Set to use GNFR emission sectors
use_GNFR_emissions_from_EMEP_flag = .true.
use_alphabetic_GNFR_emissions_from_EMEP_flag = .true.

! Set to use EMEP surface values instead of gridded to include deposition
use_EMEP_surface_ozone_flag = .true.
use_EMEP_surface_compounds_flag = .true.

! Set to use EMEP surface PM including water
use_water_in_EMEP_surface_pm_flag = .true.
```

We then set some limits to the mixing and boundary layer conditions provided by EMEP:

```fortran
! Set to contrain the minimum and maximum mixing heights
hmix_min = 100.0
hmix_max = 2000.0

! Set to contrain the minimum friction velocity
ustar_min = 0.01

! Set to contrain the lower bound for stable Monin-Obukhov length
lowest_stable_L = 25.0

! Set to contrain the upper bound for unstable Monin-Obukhob length
lowest_unstable_L = -10.0
```

Lastly, we provide information in the local fractions:

```fortran
! Set naming template used in the NetCDF file
use_local_fraction_naming_template_flag = .true.
use_local_fraction_grid_size_in_template_flag = .true.
local_fraction_naming_template_str = 'sec<n>_fraction'
```

## Emission redistribution proxy data <a name="emission-redistribution-proxy-data"></a>

In Mu et al. (2022), three GNFR sectors were downscaled including traffic, residential heating and shipping. Denby et al. (2024) extended this list to include aviation and off road emission sources. In this tutorial, we will use the simpler setup by Mu et al.

In general, emission redistribution proxy data refers to high resolution spatiotemporal data which can represent the emissions sources. For example, emissions emissions from traffic is high associated with location and size of roads within an EMEP grid cell. Thus, we assume that the sub-grid cell distribution of emissions can be proxied by roadlink information. In the following sections, we present suitable redistribution proxies for traffic, residential heating and shipping.

First we specify that local contribution should be calculated based on redistribution of EMEP emissions:

```fortran
! Set to distribute EMEP grid emissions to the existing emission subgrid using proxy emission data
subgrid_emission_distribution_flag = .true.

! Set to downscale EMEP emissions using proxy data
local_subgrid_method_flag = 3
```

### Traffic emissions <a name="traffic-emissions"></a>

As mentioned above, traffic emissions are associated with the location of roads, and to redistribute EMEP emissions, proxy data from Open Street Maps (OSM) is used. The implementation is discussed in Mu et al. (2022), and is based on the assumption that emissions are correlated with locations as well as traffic intensity (measured as average daily traffic, ADT). ADT is not provided in OSM, and instead a relative weighting scheme is used, ranging from the smallest (residential roads with a relative weight of 0.05) to the largest (highways with a relative weight of 2.0).

To include traffic emissions, we set the emission height as well as the dispersion parameters necessary for the Gaussian plume model:

```fortran
! Include traffic emission source in the calculation
calculate_source(traffic_index) = .true.

! Set the traffic emission height
h_emis(traffic_index,1) = 1.0

! Set dispersion parameters for the Gaussian plume model
sig_y_00(traffic_index,1) = 2.0
sig_z_00(traffic_index,1) = 1.0
```

Then we specify that OSM should be used to redistribute traffic emissions as well as where this information is found:

```fortran
! Set to use OSM data to redistribution traffic emissions
read_OSM_roadlink_data_flag = .true.

! Specify the header type of the OSM data
no_header_roadlink_data_flag = .false.

! Set pathname for the OSM data
pathname_rl(1) = "{local_dir}/uemep_demo/model_input/traffic/"

! Set to automatically select the OSM country based on grid position
auto_select_OSM_country_flag = .true.

! Set path and filename for the country bounding box data (used for auto country selection)
pathname_boundingbox = "{local_dir}/uemep_demo/model_input/traffic/"
filename_boundingbox = "country_bounding_box.txt"
```

Lastly, we set the conversion factors for converting proxies to emissions:

```fortran
! Set traffic emission conversion factors
emission_factor(nox_index,traffic_index,:) = 0.5
emission_factor(no2_index,traffic_index,:) = 0.075
emission_factor(pm25_index,traffic_index,:) = 0.01
emission_factor(pm10_index,traffic_index,:) = 0.01

! Set the truck-car emission ratio
ratio_truck_car_emission(nox_index) = 10.0
ratio_truck_car_emission(no2_index) = 10.0
ratio_truck_car_emission(pm25_index) = 10.0
ratio_truck_car_emission(pm10_index) = 10.0
```

### Residential heating emissions <a name="residential-heating-emissions"></a>

Emissions from residential heating is redistributed using data on building density from GHS (ref). 

To include residential heating emissions, we set the emission height as well as the dispersion parameters necessary for the Gaussian plume model:

```fortran
! Include residential heating emission source in the calculations
calculate_source(heating_index) = .true.

! Set the residential heating emission height
h_emis(heating_index) = 15.0

! Set dispersion parameters for the Gaussian plume model
sig_y_00(heating_index) = 5.0
sig_z_00(heating_index) = 10.0
```

We then specify the path and filename of the heating emissions redistribution proxy data:

```fortran
! Set path and filename for the heating proxy data
pathname_population(dwelling_index) = "{local_dir}/uemep_demo/model_input/heating/"
filename_population(dwelling_index) = "europe_buildings_popmask_latlon_demo.nc"

! Specify the proxy variable name in the NetCDF file
var_name_population_nc(dwelling_nc_index) = "N_buildings"
```

### Shipping emissions <a name="shipping-emissions"></a>

```fortran
! Include shipping emission source in the calculations
calculate_source(shipping_index) = .true.

! Set the shipping emission height
h_emis(shipping_index,1) = 45.0

! Set dispersion parameters for the Gaussian plume model
sig_y_00(shipping_index) = 5.0
sig_z_00(shipping_index) = 20.0
```

Then we specify the path and filename, as well as the type of shipping redistribution data:

```fortran
! Set to read the shipping proxy data as NetCDF
read_shipping_from_netcdf_flag = .true.

! Set data type to monthly
read_monthly_and_daily_shipping_data_flag = .true.

! Set path and filename for the shipping proxy data
pathname_ship(1) = "{local_dir}/uemep_demo/model_input/shipping/"
filename_ship(1) = "2017-2018_average_demo.nc"
```

## uEMEP processes <a name="uemep-processes"></a>

### Moving window, integration and interpolation

```fortran
! Set the moving windows size (EMEP grids)
EMEP_grid_interpolation_size = 2.0
EMEP_additional_grid_interpolation_size = 4.0

! Set the EMEP grid interpolation method for the non-local contribution in the moving window
EMEP_grid_interpolation_flag = 6

! Set the integral subgrid size scaling parameters
integral_subgrid_step = 4
integral_subgrid_delta_ref = 1000.0

! Set the number and size of local fraction grid cells
n_local_fraction_grids = 2
local_fraction_grid_size(1) = 1
local_fraction_grid_size(2) = 4

! Set the local fraction grid interpolation method
local_fraction_grid_for_EMEP_grid_interpolation = 1
local_fraction_grid_for_EMEP_additional_grid_interpolation = 2
```

### Wind and dispersion <a name="wind-and-dispersion"></a>

First we define the type of wind data used for dispersion. 

```
! Set the wind type used for dispersion
wind_level_flag = 6
wind_level_integral_flag = 1
```

Next we choose the dispersion scheme and set the parameters. Here we will use the default K scheme, using K from EMEP and variable height. 

```fortran
! Set to use the K dispersion scheme
stability_scheme_flag = 3

! Set the number of iterations in the K scheme method
n_kz_iterations = 2

! Set to use the average of the plume center of mass and emission height to calculate Kz
average_zc_h_in_Kz_flag = .true.

! Set to use changes in hourly wind direction to increase horizontal dispersion
use_last_meteo_in_dispersion = .true.

! Set to increase horizontal dispersion at low wind speeds, controlled by `FF_min_dispersion`:
use_meandering_in_dispersion = .true.
FF_min_dispersion = 0.5

! Set to scale initial horizontal standard deviation parameter to subgrid size
sigy_0_subgid_width_scale = 0.8
```

### NO2 chemistry <a name="no2-chemistry"></a>

```fortran
! Set the no2 chemistry scheme
no2_chemistry_scheme_flag = 2

! Set the no2 background chemistry scheme
no2_background_chemistry_scheme_flag = 1

! Set the Romberg parameters
romberg_parameters(1) = 30.0
romberg_parameters(2) = 35.0
romberg_parameters(3) = 0.23

! Set to use annual NO2 and O3 pdf correction
use_annual_mean_pdf_chemistry_correction = .true.
quick_annual_mean_pdf_chemistry_correction = .true.

! Set standard deviation for the pdf correction
nox_sigma_ratio_pdf = 1.14
ox_sigma_ratio_pdf = 0.21
```

### Tunnels <a name="tunnels"></a>

```fortran
! Set to use tunnel deposition scheme
use_tunnel_deposition_scheme = .true.
```

## Finalize configuration <a name="finalize-configuration"></a>

Firstly, we will set some options to increase performance and decrease memory usage:

```fortran
! Set to avoid reading EMEP data outside the region necessary for the subgrid
reduce_EMEP_region_flag = .true.
```

At this point, the primary configuration parameters are set. To finalize the configuration we need to set a few more parameters which we will not categorize, but are necessary for various technical reasons.

```fortran
! Set the region subgrid resolution
region_subgrid_delta = 250.0

! Set to adjust the wet deposition based on assumptions about the vertical distribution
adjust_wetdepo_integral_to_lowest_layer_flag = .true.

! Set the population data type
population_data_type = 2

! Set path and filenames for population density proxy data for calculation of exposure
read_population_from_netcdf_flag = .true.
pathname_population(population_index) = "{local_dir}/uemep_demo/model_input/heating/"
filename_population(population_index) = "europe_latlon_demo.nc"
var_name_population_nc(population_nc_index) = "Band1"
limit_population_delta = 250.0

! Set to derive secondary organic aerosols
derive_SOA_from_other_species = .true.

! Include EMEP sources
calculate_EMEP_source(traffic_index) = .true.
calculate_EMEP_source(shipping_index) = .true. 
calculate_EMEP_source(agriculture_index) = .true.
calculate_EMEP_source(heating_index) = .true.
calculate_EMEP_source(industry_index) = .true.
calculate_EMEP_source(publicpower_index) = .true.
calculate_EMEP_source(fugitive_index) = .true.
calculate_EMEP_source(solvents_index) = .true.
calculate_EMEP_source(aviation_index) = .true.
calculate_EMEP_source(offroad_index) = .true.
calculate_EMEP_source(waste_index) = .true.
calculate_EMEP_source(livestock_index) = .true.
calculate_EMEP_source(other_index) = .true.

use_phi_for_invL = .true.
var_name_nc(phi_nc_index) = 'phih_10m'
```

## Running uEMEP <a name="running-uemep"></a>

uEMEP is now configured to calculate downscaled concentrations of NOx, NO2, PM25 and PM10, based on EMEP emission sources from traffic, residential heating and shipping.

To run the uEMEP with this configuration, run the executable with the configuration file as the first argument, and a string for the intended date as the second string. In this case, the EMEP data are annual averages, and it is enough to provide the year followed by january 1st, i.e., `20220101`:

```bash
cd {local_dir}/uemep_demo/bin
./uemep ../config_files/uemep_demo_config.txt "20220101"
```

As uEMEP runs, it will output details to the terminal.

The output will be placed in `{local_dir}/uemep_demo/model_output/`. Further analyses of the output is beyond the scope of this tutorial, but `ncview` and `ncdump` are excellent applications to assess the content of the results:

```bash
# Go to the results folder
cd {local_dir}/uemep_demo/model_output

# Install ncview and ncdump
sudo apt install ncview netcdf-bin

# Print the content of the netcdf 
ncdump -h uEMEP_uemep_demo_20220101_00.nc

# Investigate the content of the netcdf using ncview
ncview uEMEP_uemep_demo_20220101_00.nc
```

### Multiple configuration files <a name="multiple-configuration-files"></a>

### Running in parallel <a name="running-in-parallel"></a>

## References

Denby, B.R., Gauss, M., Wind, P., Mu, Q., Wærsted, E.G., Fagerli, H., Valdebenito, A., Klein, H. 2020. Description of the uEMEP_v5 downscaling approach for the EMEP MSC-W chemistry transport model. Geoscientified Model Development 13, 6303-6323. https://doi.org/10.5194/gmd-13-6303-2020

Denby, B.R., Klimont, Z., Nyiri, A., Kiesewetter, G., Heyes, C., Fagerli, H. 2024. Future scenarios for air quality in Europe, the Western Balkans and EECCA countries: An assessment for the Gothenburg protocol review. Atmospheric Environment 333, 120602. https://doi.org/10.1016/j.atmosenv.2024.120602

Mu, Q., Denby, B.R., Wærsted, E.G., Fagerli, H. 2022. Downscaling of air pollutants in Europe using uEMEP_v6. Geoscientific Model Development 15, 449-465. https://doi.org/10.5194/gmd-15-449-2022