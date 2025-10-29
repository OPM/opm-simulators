#!/bin/bash

# Download source code

function usage() {
    echo "Usage: install.sh [options]"
    echo "Options:"
    echo "  --version=<version>    Version of Dune to install"
    echo "  --prefix=<prefix>      Prefix to install Dune to"
    echo "  --use-mpi=<yes|no>     Whether to use MPI"
    echo "  --build-type=<debug|release> Build type"
    echo "  --use-sudo=<yes|no>    Whether to use sudo (default: no)"
    echo "  --help                 Display this help message"
    echo
    echo "Example: install.sh --version=2.9.1 --prefix=/opt/dune-2.9.1 --use-mpi=no --build-type=debug --use-sudo=yes"
    exit
}

# Parse command line arguments
use_sudo=no
while [[ $# -gt 0 ]]; do
    case $1 in
        --use-sudo=*)
            use_sudo="${1#*=}"
            ;;
        --version=*)
            version="${1#*=}"
            ;;
        --prefix=*)
            prefix="${1#*=}"
            ;;
        --use-mpi=*)
            use_mpi="${1#*=}"
            ;;
        --build-type=*)
            build_type="${1#*=}"
            ;;
        --help)
            usage
            ;;
        *)
            echo "Error: Unknown option $1"
            # Display usage
            usage
            ;;
    esac
done
if [ -z "$version" ]; then
    echo "Error: Version is required"
    usage
fi
if [ -z "$prefix" ]; then
    echo "Error: Prefix is required"
    usage
fi
if [ -z "$use_mpi" ]; then
    echo "Error: Use MPI is required"
    usage
fi
if [ -z "$build_type" ]; then
    echo "Error: Build type is required"
    usage
fi
echo "Downloading source code for Dune version ${version}"
for module in common functions geometry grid istl localfunctions typetree uggrid ; do
    wget "https://dune-project.org/download/${version}/dune-${module}-${version}.tar.gz"
    tarball="dune-${module}-${version}.tar.gz"
    tar zxvf "$tarball"
    rm "$tarball"
done
echo
echo "Configuring Dune version ${version}"
sleep 1
if [ "$use_mpi" == "yes" ]; then
    disable_mpi="OFF"
else
    disable_mpi="ON"
fi
if [ "$build_type" == "debug" ]; then
    build_type="Debug"
else
    build_type="Release"
fi
dune-common-"$version"/bin/dunecontrol cmake \
   -DCMAKE_INSTALL_PREFIX=$prefix \
   -DCMAKE_BUILD_TYPE=$build_type \
   -DCMAKE_DISABLE_FIND_PACKAGE_MPI=$disable_mpi \
   -DUG_ENABLE_PARALLEL=$use_mpi
echo
echo "Building Dune version ${version}"
sleep 1
dune-common-"$version"/bin/dunecontrol make -j$(nproc)
echo
echo "Installing Dune version ${version}"
sleep 1
if [ "$use_sudo" == "yes" ]; then
    sudo dune-common-"$version"/bin/dunecontrol make install
else
    dune-common-"$version"/bin/dunecontrol make install
fi
echo
echo "Dune version ${version} installed to ${prefix}"
