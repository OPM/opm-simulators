#
# spec file for package opm-autodiff
#

%define tag rc1

Name:           opm-simulators
Version:        2017.10
Release:        0
Summary:        Open Porous Media - core library
License:        GPL-3.0
Group:          Development/Libraries/C and C++
Url:            http://www.opm-project.org/
Source0:        https://github.com/OPM/%{name}/archive/release/%{version}/%{tag}.tar.gz#/%{name}-%{version}.tar.gz
BuildRequires:  blas-devel lapack-devel dune-common-devel opm-output-devel
BuildRequires:  git suitesparse-devel doxygen bc devtoolset-6-toolchain 
BuildRequires:  opm-parser-devel opm-core-devel opm-grid-devel
BuildRequires:  tinyxml-devel dune-istl-devel eigen3-devel ecl-devel ewoms-devel
%{?el6:BuildRequires: cmake3 boost148-devel}
%{!?el6:BuildRequires: cmake boost-devel}
BuildRoot:      %{_tmppath}/%{name}-%{version}-build
Requires:       libopm-simulators1 = %{version}

%description
The Open Porous Media (OPM) initiative provides a set of open-source tools centered around the simulation of flow and transport of fluids in porous media. The goal of the initiative is to establish a sustainable environment for the development of an efficient and well-maintained software suite.

%package -n libopm-simulators1
Summary:        Open Porous Media - automatic differentiation library
Group:          System/Libraries

%description -n libopm-simulators1
The Open Porous Media (OPM) initiative provides a set of open-source tools centered around the simulation of flow and transport of fluids in porous media. The goal of the initiative is to establish a sustainable environment for the development of an efficient and well-maintained software suite.

%package devel
Summary:        Development and header files for opm-simulators
Group:          Development/Libraries/C and C++
Requires:       libopm-simulators1 = %{version}

%description devel
This package contains the development and header files for opm-simulators

%package doc
Summary:        Documentation files for opm-simulators
Group:          Documentation
BuildArch:	noarch

%description doc
This package contains the documentation files for opm-simulators

%package bin
Summary:        Applications in opm-simulators
Group:          Scientific
Requires:       libopm-simulators1 = %{version}

%description bin
This package contains the applications for opm-simulators

%prep
%setup -q -n %{name}-release-%{version}-%{tag}

%build
scl enable devtoolset-6 bash
%{?el6:cmake3} %{?!el6:cmake} -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=%{_prefix} -DCMAKE_INSTALL_DOCDIR=share/doc/%{name}-%{version} -DUSE_RUNPATH=OFF -DWITH_NATIVE=OFF -DCMAKE_CXX_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/g++ -DCMAKE_C_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/gcc %{?el6:-DBOOST_LIBRARYDIR=%{_libdir}/boost148 -DBOOST_INCLUDEDIR=%{_includedir}/boost148}
make

%install
make install DESTDIR=${RPM_BUILD_ROOT}
make install-html DESTDIR=${RPM_BUILD_ROOT}

%clean
rm -rf %{buildroot}

%post -n libopm-simulators1 -p /sbin/ldconfig

%postun -n libopm-simulators1 -p /sbin/ldconfig

%files doc
%{_docdir}/*

%files -n libopm-simulators1
%defattr(-,root,root,-)
%{_libdir}/*.so.*

%files devel
%defattr(-,root,root,-)
%{_libdir}/*.so
%{_libdir}/dunecontrol/*
%{_libdir}/pkgconfig/*
%{_includedir}/*
%{_datadir}/cmake/*

%files bin
%{_bindir}/*
