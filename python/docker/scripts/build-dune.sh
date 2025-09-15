#!/bin/bash

# This script is called from docker/Dockerfile to build the dune-common library.

set -e

BUILD_JOBS=$1
LIBTYPE=$2

if [ "$LIBTYPE" == "static" ]; then
    shared=0
    static=1
    libext="a"
else
    shared=1
    static=0
    libext="so"
fi

export CMAKE_GENERATOR=Ninja

pushd /tmp/opm

# Build dune-common
git clone --depth 1 --branch releases/opm/2024.04 https://gitlab.dune-project.org/core/dune-common.git
pushd dune-common
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=$shared -DCMAKE_POSITION_INDEPENDENT_CODE=1 -DDUNE_ENABLE_PYTHONBINDINGS=0 -DBLA_STATIC=$static -DCMAKE_DISABLE_FIND_PACKAGE_QuadMath=1 -DLAPACK_LIBRARIES=/usr/lib64/liblapack.$libext -DBLAS_LIBRARIES=/usr/lib64/libblas.$libext -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON
cmake --build . -j ${BUILD_JOBS}
cmake --build . --target install
popd

# Build dune-geometry
git clone --depth 1 --branch v2.9.1 https://gitlab.dune-project.org/core/dune-geometry.git
pushd dune-geometry
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=$shared -DCMAKE_POSITION_INDEPENDENT_CODE=1 -DDUNE_ENABLE_PYTHONBINDINGS=0 -DBLA_STATIC=$static -DCMAKE_DISABLE_FIND_PACKAGE_QuadMath=1 -DLAPACK_LIBRARIES=/usr/lib64/liblapack.$libext -DBLAS_LIBRARIES=/usr/lib64/libblas.$libext -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON
cmake --build . -j ${BUILD_JOBS}
cmake --build . --target install
popd

# Build dune-istl
git clone --depth 1 --branch v2.9.0 https://gitlab.dune-project.org/core/dune-istl.git
pushd dune-istl
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=$shared -DCMAKE_POSITION_INDEPENDENT_CODE=1 -DDUNE_ENABLE_PYTHONBINDINGS=0 -DBLA_STATIC=$static -DCMAKE_DISABLE_FIND_PACKAGE_QuadMath=1 -DLAPACK_LIBRARIES=/usr/lib64/liblapack.$libext -DBLAS_LIBRARIES=/usr/lib64/libblas.$libext -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON
cmake --build . -j ${BUILD_JOBS}
cmake --build . --target install
popd

# Build dune-uggrid
# NOTE: The core dune-uggrid repository requires user authentication to clone, so we use the staging repository.
git clone --depth 1 --branch v2.9.1 https://gitlab.dune-project.org/staging/dune-uggrid.git
pushd dune-uggrid
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=$shared -DCMAKE_POSITION_INDEPENDENT_CODE=1 -DDUNE_ENABLE_PYTHONBINDINGS=0 -DBLA_STATIC=$static -DCMAKE_DISABLE_FIND_PACKAGE_QuadMath=1 -DLAPACK_LIBRARIES=/usr/lib64/liblapack.$libext -DBLAS_LIBRARIES=/usr/lib64/libblas.$libext -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON
cmake --build . -j ${BUILD_JOBS}
cmake --build . --target install
popd

# Build dune-grid
git clone --depth 1 --branch v2.9.1 https://gitlab.dune-project.org/core/dune-grid.git
pushd dune-grid
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=$shared -DCMAKE_POSITION_INDEPENDENT_CODE=1 -DDUNE_ENABLE_PYTHONBINDINGS=0 -DBLA_STATIC=$static -DCMAKE_DISABLE_FIND_PACKAGE_QuadMath=1 -DLAPACK_LIBRARIES=/usr/lib64/liblapack.$libext -DBLAS_LIBRARIES=/usr/lib64/libblas.$libext -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON
cmake --build . -j ${BUILD_JOBS}
cmake --build . --target install
popd

# Build dune-localfunctions
git clone --depth 1 --branch v2.9.1 https://gitlab.dune-project.org/core/dune-localfunctions.git
pushd dune-localfunctions
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=$shared -DCMAKE_POSITION_INDEPENDENT_CODE=1 -DDUNE_ENABLE_PYTHONBINDINGS=0 -DBLA_STATIC=$static -DCMAKE_DISABLE_FIND_PACKAGE_QuadMath=1 -DLAPACK_LIBRARIES=/usr/lib64/liblapack.$libext -DBLAS_LIBRARIES=/usr/lib64/libblas.$libext -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON
cmake --build . -j ${BUILD_JOBS}
cmake --build . --target install
popd

popd
