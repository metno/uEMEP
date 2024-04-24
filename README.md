# uEMEP
Air quality dispersion model for high resolution downscaling of EMEP MSC-W

Comments, questions to brucerd@met.no

This github repository contains the fortran code for compiling the uEMEP model.

## This version
Version 6.0

## Installation

Download and compile the latest version:

```bash
git clone https://github.com/metno/uEMEP.git
cd uEMEP
mkdir -p build
cd build
cmake ..
make
```

Note: uEMEP supports multiple cores during the build process, e.g., `make -j 4`.

Compilation requires an Intel Fortran compiler and a compatible NetCDF installation.

## Testing

Tests are currently built by default when building uEMEP. To run the tests, simply run `make test` or `ctest` ib the build directory after running `make`.

## Implementation
The command line structure for uEMEP is as follows:
uEMEPvX.exe config_file_1 config_file_2 … config_file_10 yyyymmddHH

where X is the current version.
The file names config_file_n are up to 10 configuration files that can be read that specify the model calculation.

Each new configuration file will overwrite the previous values of the parameters specified in the new configuration file.

Parameters that are not specified will be unchanged.

The date string, required, ‘yyyymmddHH’ refers to the date string of the EMEP file to be read, specified in the configuration files.

uEMEP uses the time stamps provided by EMEP to specify the calculation times.

uEMEP requires EMEP output files, one of which contains the local fraction data, for implementation.



