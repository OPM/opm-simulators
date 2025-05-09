#!/bin/bash

set -e

VERSION=${1:-"master"}
BUILD_JOBS=${2:-16}
VERSION_TAG=${3:-""}

export CMAKE_GENERATOR=Ninja

declare -A python_versions
python_versions[cp38-cp38]=/opt/python/cp38-cp38/bin/python
python_versions[cp39-cp39]=/opt/python/cp39-cp39/bin/python
python_versions[cp310-cp310]=/opt/python/cp310-cp310/bin/python
python_versions[cp311-cp311]=/opt/python/cp311-cp311/bin/python
python_versions[cp312-cp312]=/opt/python/cp312-cp312/bin/python
python_versions[cp313-cp313]=/opt/python/cp313-cp313/bin/python

for python_bin in ${python_versions[*]}
do
  ${python_bin} -m pip install pip --upgrade
  ${python_bin} -m pip install wheel setuptools twine pytest-runner auditwheel scikit-build cmake numpy
done

DIR=`pwd`

# Setup opm modules
git clone https://github.com/OPM/opm-common -b $VERSION
git clone https://github.com/OPM/opm-grid -b $VERSION
git clone https://github.com/OPM/opm-simulators -b $VERSION
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
    cmake -DPYTHON_EXECUTABLE=${python_versions[$tag]} -DWITH_NATIVE=0 -DBoost_USE_STATIC_LIBS=1 \
    -DOPM_ENABLE_PYTHON=ON -DOPM_PYTHON_PACKAGE_VERSION_TAG=${VERSION_TAG} -DBLA_STATIC=1 -DBLAS_LIBRARIES=/usr/lib64/libblas.a -DSUITESPARSE_USE_STATIC=1 -DCMAKE_DISABLE_FIND_PACKAGE_QuadMath=1 ..

    cmake --build . --target opmcommon_python simulators --parallel ${BUILD_JOBS}

    # Package opm-common bindings
    cd opm-common/python
    ${python_versions[$tag]} setup.py sdist bdist_wheel --plat-name manylinux_2_28_x86_64 --python-tag $tag
    ${python_versions[$tag]} -m auditwheel repair dist/*$tag*.whl
    cp dist/*$tag*.whl /tmp/opm/wheelhouse
    cd ../..

    cd opm-simulators/python
    ${python_versions[$tag]} setup.py sdist bdist_wheel --plat-name manylinux_2_28_x86_64 --python-tag $tag
    ${python_versions[$tag]} -m auditwheel repair dist/*$tag*.whl
    cp dist/*$tag*.whl /tmp/opm/wheelhouse

    popd
done
