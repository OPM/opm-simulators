#
# spec file for package opm-simulators
#

%define tag final

Name:           opm-simulators
Version:        2018.10
Release:        0
Summary:        Open Porous Media - core library
License:        GPL-3.0
Group:          Development/Libraries/C and C++
Url:            http://www.opm-project.org/
Source0:        https://github.com/OPM/%{name}/archive/release/%{version}/%{tag}.tar.gz#/%{name}-%{version}.tar.gz
BuildRequires:  blas-devel lapack-devel
BuildRequires: dune-common-devel dune-geometry-devel dune-grid-devel dune-localfunctions-devel
BuildRequires:  git suitesparse-devel doxygen bc devtoolset-6-toolchain
BuildRequires:  opm-grid-devel opm-grid-openmpi-devel opm-grid-mpich-devel
BuildRequires:  ewoms-devel ewoms-openmpi-devel ewoms-mpich-devel
BuildRequires:  opm-material-devel opm-material-openmpi-devel opm-material-mpich-devel
BuildRequires:  tinyxml-devel dune-istl-devel ecl-devel zlib-devel
BuildRequires:  openmpi-devel trilinos-openmpi-devel ptscotch-openmpi-devel scotch-devel
BuildRequires:  mpich-devel trilinos-mpich-devel ptscotch-mpich-devel
BuildRequires:  opm-common-devel opm-common-openmpi-devel opm-common-mpich-devel
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

%package -n libopm-simulators1-openmpi
Summary:        Open Porous Media - automatic differentiation library
Group:          System/Libraries

%description -n libopm-simulators1-openmpi
The Open Porous Media (OPM) initiative provides a set of open-source tools centered around the simulation of flow and transport of fluids in porous media. The goal of the initiative is to establish a sustainable environment for the development of an efficient and well-maintained software suite.

%package -n libopm-simulators1-mpich
Summary:        Open Porous Media - automatic differentiation library
Group:          System/Libraries

%description -n libopm-simulators1-mpich
The Open Porous Media (OPM) initiative provides a set of open-source tools centered around the simulation of flow and transport of fluids in porous media. The goal of the initiative is to establish a sustainable environment for the development of an efficient and well-maintained software suite.

%package devel
Summary:        Development and header files for opm-simulators
Group:          Development/Libraries/C and C++
Requires:       libopm-simulators1 = %{version}

%description devel
This package contains the development and header files for opm-simulators

%package openmpi-devel
Summary:        Development and header files for opm-simulators
Group:          Development/Libraries/C and C++
Requires:       libopm-simulators1-openmpi = %{version}

%description openmpi-devel
This package contains the development and header files for opm-simulators

%package mpich-devel
Summary:        Development and header files for opm-simulators
Group:          Development/Libraries/C and C++
Requires:       libopm-simulators1-mpich = %{version}

%description mpich-devel
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

%package openmpi-bin
Summary:        Applications in opm-simulators
Group:          Scientific
Requires:       libopm-simulators1-openmpi = %{version}
Requires:       libopm-grid1-openmpi
Requires:       libopm-common1-openmpi
Requires:       trilinos-openmpi
Requires:       ptscotch-openmpi

%description openmpi-bin
This package contains the applications for opm-simulators

%package mpich-bin
Summary:        Applications in opm-simulators
Group:          Scientific
Requires:       libopm-simulators1-mpich = %{version}
Requires:       libopm-grid1-mpich
Requires:       libopm-common1-mpich
Requires:       trilinos-mpich
Requires:       ptscotch-mpich

%description mpich-bin
This package contains the applications for opm-simulators

%prep
%setup -q -n %{name}-release-%{version}-%{tag}

%build
scl enable devtoolset-6 bash
mkdir serial
cd serial
%{?el6:cmake3} %{?!el6:cmake} -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=%{_prefix} -DCMAKE_INSTALL_DOCDIR=share/doc/%{name}-%{version} -DUSE_RUNPATH=OFF -DWITH_NATIVE=OFF -DUSE_QUADMATH=0 -DCMAKE_CXX_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/g++ -DCMAKE_C_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/gcc %{?el6:-DBOOST_LIBRARYDIR=%{_libdir}/boost148 -DBOOST_INCLUDEDIR=%{_includedir}/boost148} -DCMAKE_INSTALL_SYSCONFDIR=/etc ..
make %{?_smp_mflags}
#make test
cd ..

mkdir openmpi
cd openmpi
%{?el6:module load openmpi-x86_64}
%{?!el6:module load mpi/openmpi-x86_64}
%{?el6:cmake3} %{?!el6:cmake} -DUSE_MPI=1 -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=%{_prefix}/lib64/openmpi -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_INSTALL_DOCDIR=share/doc/%{name}-%{version} -DUSE_RUNPATH=OFF -DWITH_NATIVE=OFF -DUSE_QUADMATH=0 -DCMAKE_CXX_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/g++ -DCMAKE_C_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/gcc %{?el6:-DBOOST_LIBRARYDIR=%{_libdir}/boost148 -DBOOST_INCLUDEDIR=%{_includedir}/boost148} -DZOLTAN_ROOT=/usr/lib64/openmpi -DCMAKE_CXX_FLAGS=-I/usr/include/openmpi-x86_64/trilinos -DZOLTAN_INCLUDE_DIRS=/usr/include/openmpi-x86_64/trilinos -DPTSCOTCH_ROOT=/usr/lib64/openmpi -DPTSCOTCH_INCLUDE_DIR=/usr/include/openmpi-x86_64 -DCMAKE_INSTALL_SYSCONFDIR=/etc ..
make %{?_smp_mflags}
#make test
cd ..

mkdir mpich
cd mpich
%{?el6:module rm openmpi-x86_64}
%{?el6:module load mpich-x86_64}
%{?!el6:module rm mpi/openmpi-x86_64}
%{?!el6:module load mpi/mpich-x86_64}
%{?el6:cmake3} %{?!el6:cmake} -DUSE_MPI=1 -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=%{_prefix}/lib64/mpich -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_INSTALL_DOCDIR=share/doc/%{name}-%{version} -DUSE_RUNPATH=OFF -DWITH_NATIVE=OFF -DUSE_QUADMATH=0 -DCMAKE_CXX_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/g++ -DCMAKE_C_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/gcc %{?el6:-DBOOST_LIBRARYDIR=%{_libdir}/boost148 -DBOOST_INCLUDEDIR=%{_includedir}/boost148} -DZOLTAN_ROOT=/usr/lib64/mpich -DCMAKE_CXX_FLAGS=-I/usr/include/mpich-x86_64/trilinos -DZOLTAN_INCLUDE_DIRS=/usr/include/mpich-x86_64/trilinos -DPTSCOTCH_ROOT=/usr/lib64/mpich -DPTSCOTCH_INCLUDE_DIR=/usr/include/mpich-x86_64 -DCMAKE_INSTALL_SYSCONFDIR=/etc ..
make %{?_smp_mflags}
#make test

%install
cd serial
make install DESTDIR=${RPM_BUILD_ROOT}
make install-html DESTDIR=${RPM_BUILD_ROOT}
cd ..

cd openmpi
make install DESTDIR=${RPM_BUILD_ROOT}
mv ${RPM_BUILD_ROOT}/usr/lib64/openmpi/include/* ${RPM_BUILD_ROOT}/usr/include/openmpi-x86_64/
cd ..

cd mpich
make install DESTDIR=${RPM_BUILD_ROOT}
mv ${RPM_BUILD_ROOT}/usr/lib64/mpich/include/* ${RPM_BUILD_ROOT}/usr/include/mpich-x86_64/

%clean
rm -rf %{buildroot}

%post -n libopm-simulators1 -p /sbin/ldconfig
%post -n libopm-simulators1-openmpi -p /sbin/ldconfig
%post -n libopm-simulators1-mpich -p /sbin/ldconfig

%postun -n libopm-simulators1 -p /sbin/ldconfig
%postun -n libopm-simulators1-openmpi -p /sbin/ldconfig
%postun -n libopm-simulators1-mpich -p /sbin/ldconfig

%files doc
%{_docdir}/*
/etc/bash_completion.d/*

%files -n libopm-simulators1
%defattr(-,root,root,-)
%{_libdir}/*.so.*

%files -n libopm-simulators1-openmpi
%defattr(-,root,root,-)
%{_libdir}/openmpi/lib/*.so.*

%files -n libopm-simulators1-mpich
%defattr(-,root,root,-)
%{_libdir}/mpich/lib/*.so.*

%files devel
%defattr(-,root,root,-)
%{_libdir}/*.so
/usr/lib/dunecontrol/*
%{_libdir}/pkgconfig/*
%{_includedir}/*
%{_datadir}/cmake/*
%{_datadir}/opm/cmake/Modules/*

%files openmpi-devel
%defattr(-,root,root,-)
%{_libdir}/openmpi/lib/*.so
%{_libdir}/openmpi/lib/dunecontrol/*
%{_libdir}/openmpi/lib/pkgconfig/*
%{_includedir}/openmpi-x86_64/*
%{_libdir}/openmpi/share/cmake/*
%{_libdir}/openmpi/share/opm/cmake/Modules/*

%files mpich-devel
%defattr(-,root,root,-)
%{_libdir}/mpich/lib/*.so
%{_libdir}/mpich/lib/dunecontrol/*
%{_libdir}/mpich/lib/pkgconfig/*
%{_includedir}/mpich-x86_64/*
%{_libdir}/mpich/share/cmake/*
%{_libdir}/mpich/share/opm/cmake/Modules/*

%files bin
%{_bindir}/*

%files openmpi-bin
%{_libdir}/openmpi/bin/*

%files mpich-bin
%{_libdir}/mpich/bin/*
