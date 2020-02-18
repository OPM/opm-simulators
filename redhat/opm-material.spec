#
# spec file for package opm-material
#

%define tag final

Name:           opm-material
Version:        2018.10
Release:        0
Summary:        Open Porous Media - thermodynamic framework library
License:        GPL-3.0
Group:          Development/Libraries/C and C++
Url:            http://www.opm-project.org/
Source0:        https://github.com/OPM/%{name}/archive/release/%{version}/%{tag}.tar.gz#/%{name}-%{version}.tar.gz
BuildRequires:  blas-devel lapack-devel dune-common-devel
BuildRequires:  git suitesparse-devel doxygen bc
BuildRequires:  tinyxml-devel dune-istl-devel
BuildRequires:  opm-common-devel opm-common-openmpi-devel openmpi-devel opm-common-mpich-devel mpich-devel
BuildRequires:  cmake3
%{?!el8:BuildRequires: devtoolset-8-toolchain boost148-devel}
%{?el8:BuildRequires: boost-devel}
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

%package mpich-devel
Summary:        Development and header files for opm-material with mpich
Group:          Development/Libraries/C and C++

%description mpich-devel
This package contains the development and header files for opm-material with
mpich.

%package doc
Summary:        Documentation files for opm-material
Group:          Documentation
BuildArch:	noarch

%description doc
This package contains the documentation files for opm-material

%global debug_package %{nil}

%prep
%setup -q -n %{name}-release-%{version}-%{tag}

# consider using -DUSE_VERSIONED_DIR=ON if backporting
%build
mkdir serial
cd serial
cmake3 -DENABLE_MPI=0 -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=RelWithDebInfo -DSTRIP_DEBUGGING_SYMBOLS=ON -DCMAKE_INSTALL_PREFIX=%{_prefix} -DCMAKE_INSTALL_DOCDIR=share/doc/%{name}-%{version} -DUSE_RUNPATH=OFF -DWITH_NATIVE=OFF %{?!el8:-DCMAKE_CXX_COMPILER=/opt/rh/devtoolset-8/root/usr/bin/g++ -DCMAKE_C_COMPILER=/opt/rh/devtoolset-8/root/usr/bin/gcc -DCMAKE_Fortran_COMPILER=/opt/rh/devtoolset-8/root/usr/bin/gfortran} %{?!el8:-DBOOST_LIBRARYDIR=%{_libdir}/boost148 -DBOOST_INCLUDEDIR=%{_includedir}/boost148} ..
make %{?_smp_mflags}
make test
cd ..

mkdir openmpi
cd openmpi
%{?el6:module load openmpi-x86_64}
%{?!el6:module load mpi/openmpi-x86_64}
cmake3 -DUSE_MPI=1 -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=RelWithDebInfo -DSTRIP_DEBUGGING_SYMBOLS=ON -DCMAKE_INSTALL_PREFIX=%{_prefix}/lib64/openmpi -DCMAKE_INSTALL_LIBDIR=lib -DUSE_RUNPATH=OFF -DWITH_NATIVE=OFF %{?!el8:-DCMAKE_CXX_COMPILER=/opt/rh/devtoolset-8/root/usr/bin/g++ -DCMAKE_C_COMPILER=/opt/rh/devtoolset-8/root/usr/bin/gcc -DCMAKE_Fortran_COMPILER=/opt/rh/devtoolset-8/root/usr/bin/gfortran} %{?!el8:-DBOOST_LIBRARYDIR=%{_libdir}/boost148 -DBOOST_INCLUDEDIR=%{_includedir}/boost148} -DCMAKE_INSTALL_INCLUDE_DIR=%{_prefix}/include/openmpi-x86_64 ..
make %{?_smp_mflags}
make test
cd ..

mkdir mpich
cd mpich
%{?el6:module rm openmpi-x86_64}
%{?el6:module load mpich-x86_64}
%{?!el6:module rm mpi/openmpi-x86_64}
%{?!el6:module load mpi/mpich-x86_64}
cmake3 -DUSE_MPI=1 -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=RelWithDebInfo -DSTRIP_DEBUGGING_SYMBOLS=ON -DCMAKE_INSTALL_PREFIX=%{_prefix}/lib64/mpich -DCMAKE_INSTALL_LIBDIR=lib -DUSE_RUNPATH=OFF -DWITH_NATIVE=OFF %{?!el8:-DCMAKE_CXX_COMPILER=/opt/rh/devtoolset-8/root/usr/bin/g++ -DCMAKE_C_COMPILER=/opt/rh/devtoolset-8/root/usr/bin/gcc -DCMAKE_Fortran_COMPILER=/opt/rh/devtoolset-8/root/usr/bin/gfortran} %{?!el8:-DBOOST_LIBRARYDIR=%{_libdir}/boost148 -DBOOST_INCLUDEDIR=%{_includedir}/boost148} -DCMAKE_INSTALL_INCLUDE_DIR=%{_prefix}/include/mpich-x86_64 ..
make %{?_smp_mflags}
make test

%install
cd serial
make install DESTDIR=${RPM_BUILD_ROOT}
make install-html DESTDIR=${RPM_BUILD_ROOT}
cd ..
cd openmpi
make install DESTDIR=${RPM_BUILD_ROOT}
mv ${RPM_BUILD_ROOT}/%{_libdir}/openmpi/include/* ${RPM_BUILD_ROOT}/usr/include/openmpi-x86_64/
cd ..
cd mpich
make install DESTDIR=${RPM_BUILD_ROOT}
mv ${RPM_BUILD_ROOT}/%{_libdir}/mpich/include/* ${RPM_BUILD_ROOT}/usr/include/mpich-x86_64/

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
%{_libdir}/openmpi/lib/dunecontrol/*
%{_libdir}/openmpi/lib/pkgconfig/*
%{_includedir}/openmpi-x86_64/*
%{_libdir}/openmpi/share/cmake/*
%{_libdir}/openmpi/share/opm/cmake/Modules/*

%files mpich-devel
%defattr(-,root,root,-)
%{_libdir}/mpich/lib/dunecontrol/*
%{_libdir}/mpich/lib/pkgconfig/*
%{_includedir}/mpich-x86_64/*
%{_libdir}/mpich/share/cmake/*
%{_libdir}/mpich/share/opm/cmake/Modules/*
