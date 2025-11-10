# uEMEP
Air quality dispersion model for high resolution downscaling of EMEP MSC-W

This Github repository contains the source code for the uEMEP model.

## Installation

Download and compile the latest version:

```bash
git clone https://github.com/metno/uEMEP.git
cd uEMEP
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=<build-type> ..
make
```

> Note: uEMEP supports building source files in parallel, e.g., `make -j 4`.

Compilation uses the Intel Fortran compiler (ifort) by default requires a compatible NetCDF installation. To use `gfortran` instead if `ifort`, the compiler can be prepended as:

```bash
FF=gfortran cmake -DCMAKE_BUILD_TYPE=<build-type> ..
```

There are currently three different build types, `Release`, `Debug` and `Coverage`. See the top level `CMakeLists.txt` file for the compiler flags used. The `Coverage` build type is a special build for assessing code coverage using `gcov` and currently requires the `gfortran` compiler version >= 14.

## Testing

Tests are currently built by default when building uEMEP. 

To run the tests, simply run `make test` or `ctest` in the build directory after running `make`.

## Running

For help on running uEMEP, run the following from the `build` directory:

```bash
./uemep --help
```
