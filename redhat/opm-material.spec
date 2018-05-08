#
# spec file for package opm-material
#

%define tag final

Name:           opm-material
Version:        2018.04
Release:        0
Summary:        Open Porous Media - thermodynamic framework library
License:        GPL-3.0
Group:          Development/Libraries/C and C++
Url:            http://www.opm-project.org/
Source0:        https://github.com/OPM/%{name}/archive/release/%{version}/%{tag}.tar.gz#/%{name}-%{version}.tar.gz
BuildRequires:  blas-devel lapack-devel dune-common-devel
BuildRequires:  git suitesparse-devel doxygen bc
BuildRequires:  tinyxml-devel dune-istl-devel ecl-devel zlib-devel
BuildRequires:  opm-common-devel opm-common-openmpi-devel devtoolset-6-toolchain openmpi-devel
%{?el6:BuildRequires:  cmake3 boost148-devel}
%{?!el6:BuildRequires:  cmake boost-devel}
BuildRoot:      %{_tmppath}/%{name}-%{version}-build

%description
The Open Porous Media (OPM) initiative provides a set of open-source tools centered around the simulation of flow and transport of fluids in porous media. The goal of the initiative is to establish a sustainable environment for the development of an efficient and well-maintained software suite.

%package devel
Summary:        Development and header files for opm-material
Group:          Development/Libraries/C and C++

%description devel
This package contains the development and header files for opm-material

%package openmpi-devel
Summary:        Development and header files for opm-material with openmpi
Group:          Development/Libraries/C and C++

%description openmpi-devel
This package contains the development and header files for opm-material with
openMPI.

%package doc
Summary:        Documentation files for opm-material
Group:          Documentation
BuildArch:	noarch

%description doc
This package contains the documentation files for opm-material

%prep
%setup -q -n %{name}-release-%{version}-%{tag}

# consider using -DUSE_VERSIONED_DIR=ON if backporting
%build
scl enable devtoolset-6 bash

mkdir serial
cd serial
%{?el6:cmake3} %{?!el6:cmake} -DENABLE_MPI=0 -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=RelWithDebInfo -DSTRIP_DEBUGGING_SYMBOLS=ON -DCMAKE_INSTALL_PREFIX=%{_prefix} -DCMAKE_INSTALL_DOCDIR=share/doc/%{name}-%{version} -DUSE_RUNPATH=OFF -DWITH_NATIVE=OFF -DCMAKE_CXX_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/g++ -DCMAKE_C_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/gcc -DCMAKE_Fortran_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/gfortran %{?el6:-DBOOST_LIBRARYDIR=%{_libdir}/boost148 -DBOOST_INCLUDEDIR=%{_includedir}/boost148} ..
make
make test
cd ..

mkdir openmpi
cd openmpi
%{?el6:module load openmpi-x86_64}
%{?!el6:module load mpi/openmpi-x86_64}
%{?el6:cmake3} %{?!el6:cmake} -DUSE_MPI=1 -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=RelWithDebInfo -DSTRIP_DEBUGGING_SYMBOLS=ON -DCMAKE_INSTALL_PREFIX=%{_prefix}/lib64/openmpi DUSE_RUNPATH=OFF -DWITH_NATIVE=OFF -DCMAKE_CXX_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/g++ -DCMAKE_C_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/gcc -DCMAKE_Fortran_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/gfortran %{?el6:-DBOOST_LIBRARYDIR=%{_libdir}/boost148 -DBOOST_INCLUDEDIR=%{_includedir}/boost148} -DCMAKE_INSTALL_INCLUDE_DIR=%{_prefix}/include/openmpi-x86_64 ..
make
make test

%install
cd serial
make install DESTDIR=${RPM_BUILD_ROOT}
make install-html DESTDIR=${RPM_BUILD_ROOT}
cd ..
cd openmpi
make install DESTDIR=${RPM_BUILD_ROOT}
mv ${RPM_BUILD_ROOT}/usr/lib64/openmpi/include/* ${RPM_BUILD_ROOT}/usr/include/openmpi-x86_64/

%clean
rm -rf %{buildroot}

%files doc
%{_docdir}/*

%files devel
%defattr(-,root,root,-)
/usr/lib/dunecontrol/*
/usr/lib/pkgconfig/*
%{_includedir}/*
%{_datadir}/cmake/*
%{_datadir}/opm/cmake/Modules/*

%files openmpi-devel
%defattr(-,root,root,-)
/usr/lib64/openmpi/lib/dunecontrol/*
/usr/lib64/openmpi/lib/pkgconfig/*
%{_includedir}/openmpi-x86_64/*
/usr/lib64/openmpi/share/cmake/*
/usr/lib64/openmpi/share/opm/cmake/Modules/*
