#!/bin/bash

# This script is called from docker/Dockerfile to generate Python wheels for a set of python versions.
# For each python version it:
#  - first builds the python bindings for opm-common and opm-simulators (the .so files)
#  - then uses setup.py to build python package wheels for both opm-common and opm-simulators
#  - and finally uses auditwheel to repair the wheels and copy them to the wheelhouse directory
#
set -e

# Silence pip's root-user warning inside manylinux build containers. This is a
# controlled, ephemeral environment and we intentionally run as root.
export PIP_ROOT_USER_ACTION=ignore

VERSION_COMMON=${1:-"master"}
VERSION_GRID=${2:-"master"}
VERSION_SIMULATORS=${3:-"master"}
BUILD_JOBS=${4:-16}
VERSION_TAG=${5:-""}
LIBTYPE=${6:-"static"}
PYTHON_VERSIONS_ARG=${7:-""}
TARGET_COMMON=${8:-"opmcommon_python"}
TARGET_SIMULATORS=${9:-"simulators"}

# Read default Python versions and manylinux platform from JSON config
if command -v python3 >/dev/null 2>&1 && [ -f "docker/scripts/read_python_config.py" ]; then
    DEFAULT_VERSIONS=$(python3 docker/scripts/read_python_config.py default_versions 2>/dev/null)
    if [ $? -ne 0 ] || [ -z "$DEFAULT_VERSIONS" ]; then
        echo "Error: Could not read Python versions from docker/python_versions.json"
        exit 1
    fi
    MANYLINUX_PLATFORM=$(python3 docker/scripts/read_python_config.py manylinux_platform 2>/dev/null)
    if [ $? -ne 0 ] || [ -z "$MANYLINUX_PLATFORM" ]; then
        echo "Error: Could not read manylinux platform from docker/python_versions.json"
        exit 1
    fi
else
    echo "Error: python3 or docker/scripts/read_python_config.py not found"
    exit 1
fi
PYTHON_VERSIONS=${PYTHON_VERSIONS_ARG:-"$DEFAULT_VERSIONS"}
echo "Using repository versions - Common: $VERSION_COMMON, Grid: $VERSION_GRID, Simulators: $VERSION_SIMULATORS"
echo "Using build targets - Common: $TARGET_COMMON, Simulators: $TARGET_SIMULATORS"
echo "Requested Python versions: $PYTHON_VERSIONS"

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

# TODO: These version mappings could also be read from docker/python_versions.json in the future
# NOTE: These version mappings should match python_versions.json - run sync_versions.sh after JSON changes
declare -A all_python_versions
unset python_versions  # Clear any existing array
declare -A python_versions
all_python_versions[3.8]="cp38-cp38:/opt/python/cp38-cp38/bin/python"
all_python_versions[3.9]="cp39-cp39:/opt/python/cp39-cp39/bin/python"
all_python_versions[3.10]="cp310-cp310:/opt/python/cp310-cp310/bin/python"
all_python_versions[3.11]="cp311-cp311:/opt/python/cp311-cp311/bin/python"
all_python_versions[3.12]="cp312-cp312:/opt/python/cp312-cp312/bin/python"
all_python_versions[3.13]="cp313-cp313:/opt/python/cp313-cp313/bin/python"

# Build python_versions array from requested versions
IFS=',' read -ra VERSIONS_ARRAY <<< "$PYTHON_VERSIONS"
for version in "${VERSIONS_ARRAY[@]}"; do
    version=$(echo "$version" | tr -d ' ')  # Remove whitespace
    if [[ -n "${all_python_versions[$version]}" ]]; then
        IFS=':' read -r tag path <<< "${all_python_versions[$version]}"
        python_versions[$tag]=$path
        echo "Added Python $version (tag: $tag, path: $path)"
    else
        echo "Warning: Python version $version not supported. Supported: 3.8, 3.9, 3.10, 3.11, 3.12, 3.13"
    fi
done

if [[ ${#python_versions[@]} -eq 0 ]]; then
    echo "Error: No valid Python versions specified!"
    exit 1
fi

# Install Python packages
for python_bin in ${python_versions[*]}
do
  ${python_bin} -m pip install pip --upgrade
  # NOTE (2025-09-12): Removed pytest-runner; we no longer use setup_requires
  # and do not need it at build time. Keep the minimal toolchain for wheel
  # building and auditwheel repair.
  ${python_bin} -m pip install wheel setuptools twine auditwheel scikit-build cmake numpy
done

DIR=`pwd`

# Setup opm modules. Remove existing directories to avoid conflicts
rm -rf opm-common opm-grid opm-simulators opm-utilities
echo "Cloning repositories..."

# Function to clone repository with PR support
clone_repo() {
    local repo_name=$1
    local version=$2

    if [[ "$version" == pull/* ]]; then
        # Handle PR references
        echo "Cloning $repo_name and fetching PR: $version"
        git clone --depth=1 https://github.com/OPM/$repo_name
        cd $repo_name
        git fetch origin $version  # NOTE: implies --depth=1
        git checkout FETCH_HEAD
        cd ..
    else
        # Handle regular branches/tags
        echo "Cloning $repo_name at version: $version"
        git clone --depth=1 https://github.com/OPM/$repo_name -b $version
    fi
}

# Clone each repository
clone_repo "opm-common" "$VERSION_COMMON"
clone_repo "opm-grid" "$VERSION_GRID"
clone_repo "opm-simulators" "$VERSION_SIMULATORS"
git clone https://github.com/OPM/opm-utilities

ln -sf opm-utilities/opm-super/CMakeLists.txt CMakeLists.txt
sed -e 's/add_subdirectory(opm-upscaling)//' -e 's/add_dependencies(opmupscaling opmgrid)//g' -i CMakeLists.txt

mkdir -p /tmp/opm/wheelhouse

for tag in ${!python_versions[@]}
do
    # Delete the folder if it already exists
    if [ -d $tag ]; then
      rm -rf $tag
    fi
    mkdir $tag && pushd $tag

    cmake -DPYTHON_EXECUTABLE=${python_versions[$tag]} -DWITH_NATIVE=0 -DBoost_USE_STATIC_LIBS=$static \
    -DOPM_ENABLE_PYTHON=ON -DOPM_PYTHON_PACKAGE_VERSION_TAG=${VERSION_TAG} -DBLA_STATIC=$static \
    -DBLAS_LIBRARIES=/usr/lib64/libblas.${libext} -DSUITESPARSE_USE_STATIC=$static \
    -DCMAKE_DISABLE_FIND_PACKAGE_QuadMath=1 -DBUILD_SHARED_LIBS=$shared -DOPM_INTERPROCEDURAL_OPTIMIZATION_TYPE=GNU -DOPM_INTERPROCEDURAL_OPTIMIZATION_JOBS=${BUILD_JOBS} ..

    cmake --build . --target $TARGET_COMMON $TARGET_SIMULATORS --parallel ${BUILD_JOBS}

    # Package opm-common bindings
    cd opm-common/python
    # NOTE (2025-09-12): Switch to the modern PEP 517 build frontend.
    # - Avoids setuptools deprecation warnings from direct setup.py usage.
    # - Uses isolated builds honoring pyproject.toml build-system requires.
    # - We intentionally drop --plat-name/--python-tag; environment tags and
    #   auditwheel will set/normalize these appropriately.
    # Ensure the 'build' frontend is available for this interpreter.
    ${python_versions[$tag]} -m pip -q install --upgrade build
    ${python_versions[$tag]} -m build --sdist --wheel

    # Set LD_LIBRARY_PATH so auditwheel can find and bundle OPM libraries, important when using shared==1
    export LD_LIBRARY_PATH="${PWD}/../lib:${PWD}/../../opm-grid/lib:${PWD}/../../opm-simulators/lib:$LD_LIBRARY_PATH"
    ${python_versions[$tag]} -m auditwheel repair dist/*$tag*.whl
    cp wheelhouse/*$tag*.whl /tmp/opm/wheelhouse  # NOTE: auditwheel puts repaired wheels in wheelhouse/
    cd ../..

    # Package opm-simulators bindings
    cd opm-simulators/python
    # (build already installed for this interpreter above)
    ${python_versions[$tag]} -m build --sdist --wheel

    # Set LD_LIBRARY_PATH so auditwheel can find and bundle OPM libraries, important when using shared==1
    export LD_LIBRARY_PATH="${PWD}/../lib:${PWD}/../../opm-common/lib:${PWD}/../../opm-grid/lib:$LD_LIBRARY_PATH"
    ${python_versions[$tag]} -m auditwheel repair dist/*$tag*.whl
    cp wheelhouse/*$tag*.whl /tmp/opm/wheelhouse  # NOTE: auditwheel puts repaired wheels in wheelhouse/

    popd
done
