# Important

This branch contains the legacy source code for uEMEP (version <= 6) and is no longer updated.

## uEMEP
Air quality dispersion model for high resolution downscaling of EMEP MSC-W

Comments, questions to brucerd@met.no

## This version
Version 5.0

## Installation
This github repository contains the fortran code and makefiles for compiling the uEMEP model.

Download/pull these files and adjust the makefile for your directory configuration.

Compilation requires an intel fortran compiler.

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

## Example files and configuration
A zip file 'uEMEP_demo_startup_files.zip' containing the data and configuration files for first implementing uEMEP.


