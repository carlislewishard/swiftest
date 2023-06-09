# Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 

######################################################
# Determine and set the Fortran compiler flags we want 
######################################################

####################################################################
# Make sure that the default build type is RELEASE if not specified.
####################################################################
INCLUDE(${CMAKE_MODULE_PATH}/SetCompileFlag.cmake)

# Make sure the build type is uppercase
STRING(TOUPPER "${CMAKE_BUILD_TYPE}" BT)

SET(BUILD_TYPE_MSG "Choose the type of build, options are DEBUG, RELEASE, PROFILE, or TESTING.")

IF(BT STREQUAL "RELEASE")
    SET(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      ${BUILD_TYPE_MSG}
      FORCE)
ELSEIF(BT STREQUAL "DEBUG")
    SET (CMAKE_BUILD_TYPE DEBUG CACHE STRING
      ${BUILD_TYPE_MSG}
      FORCE)
ELSEIF(BT STREQUAL "TESTING")
    SET (CMAKE_BUILD_TYPE TESTING CACHE STRING
      ${BUILD_TYPE_MSG}
      FORCE)
ELSEIF(BT STREQUAL "PROFILE")
    SET (CMAKE_BUILD_TYPE PROFILE CACHE STRING
      ${BUILD_TYPE_MSG}
      FORCE)
ELSEIF(NOT BT)
    SET(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      ${BUILD_TYPE_MSG}
      FORCE)
    MESSAGE(STATUS "CMAKE_BUILD_TYPE not given, defaulting to RELEASE")
ELSE()
    MESSAGE(FATAL_ERROR "CMAKE_BUILD_TYPE not valid! ${BUILD_TYPE_MSG}")
ENDIF(BT STREQUAL "RELEASE")

IF (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    IF (APPLE)
        SET(MACHINE_CODE_VALUE "tune=native" CACHE STRING "Tells the compiler which processor features it may target, including which instruction sets and optimizations it may generate.")
    ELSE ()
        SET(MACHINE_CODE_VALUE "arch=native" CACHE STRING "Tells the compiler which processor features it may target, including which instruction sets and optimizations it may generate.")
    ENDIF ()
ELSE ()
    SET(MACHINE_CODE_VALUE "host" CACHE STRING "Tells the compiler which processor features it may target, including which instruction sets and optimizations it may generate.")
ENDIF ()

#########################################################
# If the compiler flags have already been set, return now
#########################################################

IF(CMAKE_Fortran_FLAGS_RELEASE AND CMAKE_Fortran_FLAGS_TESTING AND CMAKE_Fortran_FLAGS_DEBUG AND CMAKE_Fortran_FLAGS_PROFILE )
    RETURN ()
ENDIF(CMAKE_Fortran_FLAGS_RELEASE AND CMAKE_Fortran_FLAGS_TESTING AND CMAKE_Fortran_FLAGS_DEBUG AND CMAKE_Fortran_FLAGS_PROFILE)

########################################################################
# Determine the appropriate flags for this compiler for each build type.
# For each option type, a list of possible flags is given that work
# for various compilers.  The first flag that works is chosen.
# If none of the flags work, nothing is added (unless the REQUIRED 
# flag is given in the call).  This way unknown compiles are supported.
#######################################################################

#####################
### GENERAL FLAGS ###
#####################

# Free form
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                Fortran "-ffree-form" # GNU
                ) 

# Don't add underscores in symbols for C-compatability
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                Fortran "-fno-underscoring" # GNU
                ) 
# Compile code assuming that IEEE signaling NaNs may generate user-visible traps during floating-point operations. 
# Setting this option disables optimizations that may change the number of exceptions visible with signaling NaNs. 
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                Fortran "-fsignaling-nans " # GNU
                ) 
               
# Allows for lines longer than 80 characters without truncation
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                Fortran "-ffree-line-length-none" # GNU (gfortran)
                )

# Disables right margin wrapping in list-directed output
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                Fortran  "-no-wrap-margin" # Intel
                         "/wrap-margin-"   # Intel Windows        
                )

# Aligns a variable to a specified boundary and offset
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                Fortran "-align all -align array64byte" # Intel
                        "/align:all /align:array64byte" # Intel Windows
                )

# Enables changing the variable and array memory layout
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                Fortran "-pad"  # Intel
                        "/Qpad" # Intel Windows

                )

IF (NOT BUILD_SHARED_LIBS)
        # Use static Intel libraries
        SET_COMPILE_FLAG(CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS}"
                Fortran "-static-intel"  # Intel
        )
        # Use static Intel MPI libraries
        SET_COMPILE_FLAG(CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS}"
                Fortran "-static_mpi"  # Intel
        )

        IF (USE_OPENMP)
                SET_COMPILE_FLAG(CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS}"
                        Fortran "-qopenmp-link=static"  # Intel
                )
        ENDIF (USE_OPENMP)
ENDIF ()

IF (USE_SIMD)
        # Enables OpenMP SIMD compilation when OpenMP parallelization is disabled. 
        IF (NOT USE_OPENMP)
                SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                        Fortran "-qno-openmp -qopenmp-simd" # Intel
                        Fortran "/Qopenmp- /Qopenmp-simd" # Intel Windows
                        )     
        ENDIF (NOT USE_OPENMP)

        # Optimize for an old enough processor that it should run on most computers
        SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                        Fortran "-x${MACHINE_CODE_VALUE}" # Intel
                                "/Qx${MACHINE_CODE_VALUE}" # Intel Windows
                                "-m${MACHINE_CODE_VALUE}"    # GNU
                        )

        # Generate an extended set of vector functions
        SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                        Fortran "-vecabi=cmdtarget" # Intel
                        Fortran "/Qvecabi:cmdtarget" # Intel Windows
                        )
ENDIF (USE_SIMD)                




###################
### DEBUG FLAGS ###
###################

# Disable optimizations
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran REQUIRED "-O0" # All compilers not on Windows
                                  "/Od" # Intel Windows
                                  "-Og" # GNU (gfortran)
                )

# Turn on all warnings 
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-warn all" # Intel
                         "/warn:all" # Intel Windows
                         "-Wall"     # GNU
                )
# This enables some extra warning flags that are not enabled by -Wall
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-Wextra" # GNU
                )

# Disable the warning that arrays may be uninitialized, which comes up due to a known bug in gfortran
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-Wno-maybe-uninitialized" # GNU
                )
# Disable the warning about unused dummy arguments. These primarily occur due to interface rules for type-bound procedures used in extendable types.
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-Wno-unused-dummy-argument" # GNU
                )

# Tells the compiler to issue compile-time messages for nonstandard language elements (Fortran 2018).                
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran "-stand f18"  # Intel
                        "/stand:f18"  # Intel Windows
                        "-fstd=f2018" # GNU
                )  

# Traceback
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-traceback"   # Intel Group
                         "/traceback"   # Intel Windows
                         "-fbacktrace"  # GNU (gfortran)
                )
# Sanitize
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran "-fsanitize=address, undefined"  # Gnu 
                )
                

# Check everything
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran "-check all"       # Intel
                         "/check:all"      # Intel Windows
                         "-fcheck=all" # GNU 
                )
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran "-fstack-check" # GNU 
                )
                

# Initializes matrices/arrays with NaN values
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran "-init=snan,arrays"  # Intel
                        "/Qinit:snan,arrays" # Intel Windows
                        "-finit-real=snan"   # GNU
                )

# Does not generate an interface block for each routine in a source file
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-nogen-interfaces" # Intel
                         "/nogen-interfaces" # Intel Windows
                )

# Does not generate aposition independent executable
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-no-pie" # Intel
                )

# Does not set denormal results from floating-point calculations to zero
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-no-ftz" # Intel
                         "/Qftz-"  # Intel Windows
                )

# Enables floating-point invalid, divide-by-zero, and overflow exceptions
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-fpe-all=0"                         # Intel
                         "/fpe-all:0"                         # Intel Windows
                         "-ffpe-trap=zero,overflow,underflow" # GNU
                )

# List of floating-point exceptions, whose flag status is printed to ERROR_UNIT when invoking STOP and ERROR STOP
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-ffpe-summary=all" # GNU
                )

# Enables floating-point invalid, divide-by-zero, and overflow exceptions
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-fpe0"  # Intel
                         "/fpe:0" # Intel Windows
                )

# Enables debug info
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-debug all" # Intel
                         "/debug:all" # Intel Windows
                )

# Disables additional interprocedural optimizations for a single file compilation
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-no-ip" # Intel
                         "/Qip-"  # Intel Windows
                )

# Disables prefetch insertion optimization
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran "-qno-opt-prefetch" # Intel
                        "/Qopt-prefetch-"   # Intel Windows
                )
                
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-fstack-check" # GNU 
                )   

#####################
### TESTING FLAGS ###
#####################

# Optimizations
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_TESTING "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran REQUIRED "-O3" # All compilers not on Windows
                                 "/O3" # Intel Windows
                )

#####################
### RELEASE FLAGS ###
#####################
# NOTE: agressive optimizations (-O3) are already turned on by default

# Unroll loops
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                Fortran "-unroll"        # Intel
                        "/unroll"        # Intel Windows
                        "-funroll-loops" # GNU
                )

# Inline functions
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                 Fortran "-inline"            # Intel
                         "/Qinline"           # Intel Windows
                         "-finline-functions" # GNU
                )

# Calls the Matrix Multiply library
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                Fortran "-qopt-matmul" # Intel
                        "/Qopt-matmul" # Intel Windows
                )

# Saves the compiler options and version number to the executable 
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                Fortran "-sox" # Intel
                )

# Aligns a variable to a specified boundary and offset
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                Fortran "-align all" # Intel
                        "/align:all" # Intel Windows
                )

# No floating-point exceptions
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                Fortran "-fp-model no-except" # Intel
                        "/fp:no-except"       # Intel Windows
                )

# Generate fused multiply-add instructions
 SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                Fortran "-fma"  # Intel
                        "/Qfma" # Intel Windows
                )

# Tells the compiler to link to certain libraries in the Intel oneAPI Math Kernel Library (oneMKL). 
 SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                Fortran "-qmkl=cluster" # Intel
                        "-qmkl"         # Intel
                        "/Qmkl:cluster" # Intel Windows
                        "/Qmkl"         # Intel Windows
                ) 

# Enables additional interprocedural optimizations for a single file compilation
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                Fortran "-ip"  # Intel
                        "/Qip" # Intel Windows
                )

 
#####################
### MATH FLAGS ###
#####################
# Some subroutines require more strict floating point operation optimizations for repeatability
SET_COMPILE_FLAG(STRICTMATH_FLAGS "${STRICTMATH_FLAGS}"
                Fortran "-fp-model=precise" # Intel
                         "/fp:precise" # Intel Windows 
                         "-fno-unsafe-math-optimizations" # GNU
                )
# Disable transformations and optimizations that assume default floating-point rounding behavior. 
SET_COMPILE_FLAG(STRICTMATH_FLAGS "${STRICTMATH_FLAGS}"
                Fortran "-frounding-math"
                )

SET_COMPILE_FLAG(STRICTMATH_FLAGS "${STRICTMATH_FLAGS}"
                Fortran "-prec-div"  # Intel
                        "/Qprec-div" # Intel Windows 
                ) 

SET_COMPILE_FLAG(STRICTMATH_FLAGS "${STRICTMATH_FLAGS}"
                Fortran "-prec-sqrt"   # Intel
                         "/Qprec-sqrt" # Intel Windows 
                )

SET_COMPILE_FLAG(STRICTMATH_FLAGS "${STRICTMATH_FLAGS}"
                Fortran "-assume protect-parens" # Intel
                        "/assume:protect-parens" # Intel Windows 
                ) 

# Improves floating-point precision and consistency
SET_COMPILE_FLAG(STRICTMATH_FLAGS "${STRICTMATH_FLAGS}"
                Fortran "-mp1"   # Intel
                        "/Qprec" # Intel Windows
                ) 

# Most subroutines can use aggressive optimization of floating point operations without problems.           
SET_COMPILE_FLAG(FASTMATH_FLAGS "${FASTMATH_FLAGS}"
                Fortran "-fp-model=fast" # Intel
                        "/fp:fast"       # Intel Windows
                        "-ffast-math"    # GNU
                )


SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}" 
                Fortran ${STRICTMATH_FLAGS}
                )

#####################
### PROFILE FLAGS ###
#####################
# Enables the optimization reports to be generated
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_RELEASE}"
                Fortran "-O2 -pg -qopt-report=5 -traceback -p -g3" # Intel
                        "/O2 /Qopt-report:5 /traceback /Z7"        # Intel Windows
                        "-O2 -pg -fbacktrace"                      # GNU
                )
