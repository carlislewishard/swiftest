#!/bin/bash
# This script will build the Swiftest package. It is assumed that compatible dependencies have been built already before this is run
# 
# Copyright 2023 - David Minton
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 
SCRIPT_DIR=$(realpath $(dirname $0))
ROOT_DIR=$(realpath ${SCRIPT_DIR}/..)
USTMT="Usage: $0 <{Intel}|GNU> [/path/to/nf-nc-hdf5]"
if [[ ( $@ == "--help") ||  $@ == "-h" ]]; then 
	echo $USTMT
	exit 0
fi
COMPILER=${1:-Intel}
DEPDIR={$2:-$(realpath ${ROOT_DIR}/build)}

case $COMPILER in
    Intel)
        if command -v ifx &> /dev/null; then
            export FC=$(command -v ifx)
            export CC=$(command -v icx)
            export CXX=$(command -v icpx)
        elif command -v ifort &> /dev/null; then
            export FC=$(command -v ifort) 
            export CC=$(command -v icc)
            export CXX=$(command -v icpc)
        else
            echo "Error. Cannot find valid Intel fortran compiler."
            exit 1
        fi
        export F77="${FC}"
        ;;
    GNU)
        export FC=$(command -v gfortran)
        export CC=$(command -v gcc)
        export CXX=$(command -v g++)
        ;;
    *)
        echo "Unknown compiler type: ${COMPILER}"
        echo $USTMT
        exit 1
        ;;
esac
export F77=${FC}
echo "Using $COMPILER compilers:\nFC: $FC\nCC: $CC\nCXX: $CXX\n"

export CPATH=$DEPDIR
export LD_LIBRARY_PATH="${CPATH}/lib:${LD_LIBRARY_PATH}"
export LIBS="-lhdf5_hl -lhdf5 -lz"
export LDFLAGS="-L${DEPDIR}/lib"
export CFLAGS="-fPIC"
export CMAKE_ARGS="-DBUILD_SHARED_LIBS=OFF"

if [ $COMPILER = "Intel" ]; then 
    export FCFLAGS="${CFLAGS} -standard-semantics"
    export FFLAGS=${CFLAGS}
else
    export FCFLAGS="${CFLAGS}"
    export FFLAGS="${CFLAGS}"
fi
cd $ROOT_DIR
pip install .
