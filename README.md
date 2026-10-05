# uEMEP
Air quality dispersion model for high resolution downscaling of EMEP MSC-W.

This repository contains the source code for the uEMEP model.

## Installation

Download and compile the latest version:

```bash
git clone https://github.com/metno/uEMEP.git
cd uEMEP
cmake -S . -B release
cmake --build release -j 8
```

This builds an optimized version of uEMEP. Intel Fortran (`ifort`) is used by
default.

To use GNU Fortran (`gfortran`) instead, configure a separate build directory:

```bash
FC=gfortran cmake -S . -B release-gnu
cmake --build release-gnu -j 8
```

Building requires a Fortran compiler and a NetCDF-Fortran installation built
with that compiler.

## Build types

| Build type | Purpose |
|------------|---------|
| `Release` | Optimized build for normal use. |
| `Check` | Build with extra checks to help find errors when testing. |
| `Debug` | Build for investigating problems with a debugger. |
| `Coverage` | Build for measuring code coverage; **requires GNU Fortran 14 or newer and `gcov`**. |

`Release` is selected by default. To choose another build type, add
`-DCMAKE_BUILD_TYPE=<type>` to the configure command. For example, replace the
Release configure command above with:

```bash
cmake -S . -B debug -DCMAKE_BUILD_TYPE=Debug
cmake --build debug -j 8
```

For a coverage build, use GNU Fortran and a separate build directory:

```bash
FC=gfortran cmake -S . -B coverage -DCMAKE_BUILD_TYPE=Coverage
cmake --build coverage -j 8
```

## Testing

Unit tests are included in the build. To run them with extra error checks:

```bash
cmake -S . -B check -DCMAKE_BUILD_TYPE=Check
cmake --build check -j 8
cd check
ctest --output-on-failure
```

Run the tests from inside the build directory (otherwise tests depending on relative paths will fail).

## Running

For help on running uEMEP, run this from the build directory:

```bash
./uemep --help
```
