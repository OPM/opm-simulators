#
# spec file for package opm-material
#

Name:           opm-material
Version:        2013.09
Release:        0
Summary:        Open Porous Media - thermodynamic framework library
License:        GPL-3.0
Group:          Development/Libraries/C and C++
Url:            http://www.opm-project.org/
Source0:        %{name}-%{version}.tar.gz
BuildRequires:  blas-devel gcc-c++ gcc-gfortran lapack-devel dune-common-devel
BuildRequires:  boost-devel git suitesparse-devel cmake28 doxygen bc
BuildRequires:  tinyxml-devel dune-istl-devel
BuildRoot:      %{_tmppath}/%{name}-%{version}-build

%description
The Open Porous Media (OPM) initiative provides a set of open-source tools centered around the simulation of flow and transport of fluids in porous media. The goal of the initiative is to establish a sustainable environment for the development of an efficient and well-maintained software suite.

%package devel
Summary:        Development and header files for opm-material
Group:          Development/Libraries/C and C++
Requires:       %{name} = %{version}

%description devel
This package contains the development and header files for opm-material

%package doc
Summary:        Documentation files for opm-material
Group:          Documentation
BuildArch:	noarch

%description doc
This package contains the documentation files for opm-material

%prep
%setup -q

# consider using -DUSE_VERSIONED_DIR=ON if backporting
%build
cmake28 -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=%{_prefix} -DCMAKE_INSTALL_DOCDIR=share/doc/%{name}-%{version} -DUSE_RUNPATH=OFF
make

%install
make install DESTDIR=${RPM_BUILD_ROOT}
make install-html DESTDIR=${RPM_BUILD_ROOT}

%clean
rm -rf %{buildroot}

%files doc
%{_docdir}/*

%files devel
%defattr(-,root,root,-)
%{_libdir}/dunecontrol/*
%{_libdir}/pkgconfig/*
%{_includedir}/*
%{_datadir}/cmake/*
