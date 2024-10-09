#!/bin/bash

set -e

BUILD_JOBS=$1

export CMAKE_GENERATOR=Ninja

pushd /tmp/opm

# Build boost
git clone --depth 1 --branch boost-1.84.0 https://github.com/boostorg/boost
pushd boost
git submodule init
git submodule update
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=0 -DCMAKE_POSITION_INDEPENDENT_CODE=1
cmake --build . -- -j${BUILD_JOBS}
cmake --build . --target install
popd

# Build dune-common
git clone --depth 1 --branch releases/opm/2024.04 https://gitlab.dune-project.org/core/dune-common.git
pushd dune-common
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=0 -DCMAKE_POSITION_INDEPENDENT_CODE=1 -DDUNE_ENABLE_PYTHONBINDINGS=0 -DBLA_STATIC=1 -DCMAKE_DISABLE_FIND_PACKAGE_QuadMath=1 -DBLAS_LIBRARIES=/usr/lib64/libblas.a -DLAPACK_LIBRARIES=/usr/lib64/liblapack.a
cmake --build . -- -j${BUILD_JOBS}
cmake --build . --target install
popd

# Build dune-geometry
git clone --depth 1 --branch v2.9.1 https://gitlab.dune-project.org/core/dune-geometry.git
pushd dune-geometry
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=0 -DCMAKE_POSITION_INDEPENDENT_CODE=1 -DDUNE_ENABLE_PYTHONBINDINGS=0 -DBLA_STATIC=1 -DCMAKE_DISABLE_FIND_PACKAGE_QuadMath=1  -DBLAS_LIBRARIES=/usr/lib64/libblas.a -DLAPACK_LIBRARIES=/usr/lib64/liblapack.a
cmake --build . -- -j${BUILD_JOBS}
cmake --build . --target install
popd

# Build dune-istl
git clone --depth 1 --branch releases/opm/2024.04 https://gitlab.dune-project.org/core/dune-istl.git
pushd dune-istl
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=0 -DCMAKE_POSITION_INDEPENDENT_CODE=1 -DDUNE_ENABLE_PYTHONBINDINGS=0 -DBLA_STATIC=1 -DCMAKE_DISABLE_FIND_PACKAGE_QuadMath=1  -DBLAS_LIBRARIES=/usr/lib64/libblas.a -DLAPACK_LIBRARIES=/usr/lib64/liblapack.a
cmake --build . -- -j${BUILD_JOBS}
cmake --build . --target install
popd

# Build dune-uggrid
git clone --depth 1 --branch v2.9.1 https://gitlab.dune-project.org/staging/dune-uggrid.git
pushd dune-uggrid
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=0 -DCMAKE_POSITION_INDEPENDENT_CODE=1 -DDUNE_ENABLE_PYTHONBINDINGS=0 -DBLA_STATIC=1 -DCMAKE_DISABLE_FIND_PACKAGE_QuadMath=1  -DBLAS_LIBRARIES=/usr/lib64/libblas.a -DLAPACK_LIBRARIES=/usr/lib64/liblapack.a
cmake --build . -- -j${BUILD_JOBS}
cmake --build . --target install
popd

# Build dune-grid
git clone --depth 1 --branch v2.9.1 https://gitlab.dune-project.org/core/dune-grid.git
pushd dune-grid
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=0 -DCMAKE_POSITION_INDEPENDENT_CODE=1 -DDUNE_ENABLE_PYTHONBINDINGS=0 -DBLA_STATIC=1 -DCMAKE_DISABLE_FIND_PACKAGE_QuadMath=1  -DBLAS_LIBRARIES=/usr/lib64/libblas.a -DLAPACK_LIBRARIES=/usr/lib64/liblapack.a
cmake --build . -- -j${BUILD_JOBS}
cmake --build . --target install
popd

# Build dune-localfunctions
git clone --depth 1 --branch v2.9.1 https://gitlab.dune-project.org/core/dune-localfunctions.git
pushd dune-localfunctions
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=0 -DCMAKE_POSITION_INDEPENDENT_CODE=1 -DDUNE_ENABLE_PYTHONBINDINGS=0 -DBLA_STATIC=1 -DCMAKE_DISABLE_FIND_PACKAGE_QuadMath=1  -DBLAS_LIBRARIES=/usr/lib64/libblas.a -DLAPACK_LIBRARIES=/usr/lib64/liblapack.a
cmake --build . -- -j${BUILD_JOBS}
cmake --build . --target install
popd

popd
