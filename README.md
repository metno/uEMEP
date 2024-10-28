# uEMEP
Air quality dispersion model for high resolution downscaling of EMEP MSC-W

This Github repository contains the source code for the uEMEP model.

## This version
Version 7.0.0

## Installation

Download and compile the latest version:

```bash
git clone https://github.com/metno/uEMEP.git
cd uEMEP
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

Note: uEMEP supports building source files in parallel, e.g., `make -j 4`.

Compilation requires an Intel Fortran compiler and a compatible NetCDF installation.

## Testing

Tests are currently built by default when building uEMEP. 

To run the tests, simply run `make test` or `ctest` in the build directory after running `make`.

## Running

For help on running uEMEP, run the following from the `build` directory:

```bash
./uemep --help
```
