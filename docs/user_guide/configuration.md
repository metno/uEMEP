---
title: Configuration
---

uEMEP is flexible and has a large amount of configuration options. The syntax follows Fortran syntax, similar to the use of namelists.

## Contents:
- [File tag](#file-tag)
- [Temporal scale](#temporal-scale)
- [Emission sources](#emission-sources)
- [EMEP forcing data](#emep-forcing-data)

## File tag <a name="file-tag"></a>

uEMEP writes the model output as NetCDF format and the file tag is used to create file names associated with the specific calculation.

```fortran
file_tag = 'example'
```

## Temporal scale <a name="temporal-scale"></a>

uEMEP can be run on either the hourly or the annual scale. By default, both of these options are set to false, and one must be set to true for uEMEP to run.

```fortran
! Set calculations to be performed at the hourly scale
hourly_calculations = .true.

! Set calculations to be performed at the annual scale
annual_calculations = .true.
```

Note that setting both to true will not work and will result in the model crashing.

## Emission sources <a name="emission-sources"></a>

uEMEP can include a number of emission sources, which by default are set to false. One or more emission sources should be set to true for uEMEP to run. In addition to setting the inclusion of the emission sources to true, valid proxy redistribution data is required. Currently, uEMEP includes procedures to include the following emission sources:

```fortran
! Include emissions from traffic in the calculations
calculate_source(traffic_index) = .true.

! Include emissions from shipping in the calculations
calculate_source(shipping_index) = .true.

! Include emissions from agriculture in the calculations
calculate_source(agriculture_index) = .true.

! Include emissions from residential heating in the calculations
calculate_source(heating_index) = .true.

! Include emissions from industrial sources in the calculations
calculate_source(industry_index) = .true.
```
## EMEP forcing data <a name="emep-forcing-data"></a>

uEMEP will downscale output from the EMEP model, and thus requires an EMEP simulation with output from local fractions included to be present. The path and filename pointing to the EMEP model output should be provided for uEMEP to run.

```fortran
! Set the path and filename for the EMEP NetCDF files including meteorology
pathname_EMEP(1) = 'path_to_emep_files'
filename_EMEP(1) = 'emep_file.nc'

! Set the path and filename for the EMEP NetCDF files including local fractions
pathname_EMEP(2) = 'path_to_emep_lf_files'
filename_EMEP(2) = 'emep_lf_file.nc'
```

