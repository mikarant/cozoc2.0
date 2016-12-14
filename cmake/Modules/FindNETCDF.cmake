include (LibFindMacros)

libfind_pkg_check_modules (NETCDF_PKG_CONFIG netcdf)

find_path(NETCDF_INCLUDE_DIR
  NAMES netcdf.h
  PATHS ${NETCDF_PKG_CONFIG_INCLUDE_DIRS})

find_library(NETCDF_LIBRARY
  NAMES netcdf
  PATHS ${NETCDF_PKG_CONFIG_LIBRARY_DIRS})

set (NETCDF_PROCESS_INCLUDES
  NETCDF_INCLUDE_DIR)
set (NETCDF_PROCESS_LIBS
  NETCDF_LIBRARY)
libfind_process (NETCDF)
