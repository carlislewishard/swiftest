# Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 

# - Finds the NetCDF libraries 

find_path(NETCDF_INCLUDE_DIR NAMES netcdf.mod HINTS ENV NETCDF_FORTRAN_HOME ENV CPATH)
find_library(NETCDF_FORTRAN_LIBRARY NAMES netcdff HINTS ENV NETCDF_FORTRAN_HOME ENV LD_LIBRARY_PATH)
find_library(NETCDF_LIBRARY NAMES netcdf HINTS ENV NETCDF_FORTRAN_HOME ENV LD_LIBRARY_PATH)

set(NETCDF_FOUND TRUE)
set(NETCDF_INCLUDE_DIRS ${NETCDF_INCLUDE_DIR})
# Note for posterity: When building static libraries, NETCDF_FORTRAN_LIBRARY must come *before* NETCDF_LIBRARY. Otherwise you get a bunch of "undefined reference to" errors
set(NETCDF_LIBRARIES ${NETCDF_FORTRAN_LIBRARY} ${NETCDF_LIBRARY})

mark_as_advanced(NETCDF_LIBRARY NETCDF_FORTRAN_LIBRARY NETCDF_INCLUDE_DIR)