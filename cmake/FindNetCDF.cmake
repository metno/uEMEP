# Simple script for finding NetCDF using nf-config
#
# Note that the script will always report that the interface has been found

# Find NetCDF using nf-config
execute_process(COMMAND nf-config --includedir OUTPUT_VARIABLE NC_INC)
execute_process(COMMAND nf-config --flibs OUTPUT_VARIABLE NC_LIB)

# Strip trailing whitespace
string(STRIP ${NC_LIB} NC_LIB)

# Pretend we always have interface
set(NC_HAS_INTERFACE "YES")

# Handle arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NetCDF DEFAULT_MSG NC_LIB NC_INC NC_HAS_INTERFACE)
mark_as_advanced(NC_LIB NC_INC)

