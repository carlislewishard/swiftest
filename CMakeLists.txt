# Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 

# CMake project file for FOO

##################################################
# Define the project and the depencies that it has
##################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5...3.20.1)
PROJECT(Swiftest Fortran)

# Set the Swiftest version
SET(VERSION 1.0.0)

# Add our local modlues to the module path
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Uncomment if it is required that Fortran 90 is supported
IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

# Set some options the user may choose
# Uncomment the below if you want the user to choose a parallelization library
OPTION(USE_MPI "Use the MPI library for parallelization" OFF)
OPTION(USE_OPENMP "Use OpenMP for parallelization" OFF)

# This INCLUDE statement executes code that sets the compile flags for DEBUG,
# RELEASE, and TESTING.  You should  review this file and make sure the flags 
# are to your liking.
INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake) 
# Locate and set parallelization libraries.  There are some CMake peculiarities
# taken care of here, such as the fact that the FindOpenMP routine doesn't know
# about Fortran.
INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)

INCLUDE(${CMAKE_MODULE_PATH}/SetUpNetCDF.cmake)

# There is an error in CMAKE with this flag for pgf90.  Unset it
GET_FILENAME_COMPONENT(FCNAME ${CMAKE_Fortran_COMPILER} NAME)
IF(FCNAME STREQUAL "pgf90")
    UNSET(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS)
ENDIF(FCNAME STREQUAL "pgf90")

############################################################
# Define the actual files and folders that make up the build
############################################################

# Define the executable name
SET(FOOEXE swiftest_driver)

# Define some directories
SET(SRC ${CMAKE_SOURCE_DIR}/src)
SET(LIB ${CMAKE_SOURCE_DIR}/lib)
SET(BIN ${CMAKE_SOURCE_DIR}/bin)

# Have the .mod files placed in the lib folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${LIB})

# The source for the FOO binary and have it placed in the bin folder
ADD_SUBDIRECTORY(${SRC} ${BIN})

# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean 
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/distclean.cmake
)