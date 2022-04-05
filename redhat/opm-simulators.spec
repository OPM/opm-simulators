#
# spec file for package opm-simulators
#

%define tag final
%define rtype release
%define toolset devtoolset-9
%define build_openmpi 1
%define build_openmpi3 1
%define build_mpich 1

Name:           opm-simulators
Version:        2018.10
Release:        0
Summary:        Open Porous Media - core library
License:        GPL-3.0
Group:          Development/Libraries/C and C++
Url:            http://www.opm-project.org/
Source0:        https://github.com/OPM/%{name}/archive/release/%{version}/%{tag}.tar.gz#/%{name}-%{version}.tar.gz
BuildRequires:  blas-devel lapack-devel openblas-devel
BuildRequires:  git suitesparse-devel doxygen bc graphviz
BuildRequires:  tinyxml-devel zlib-devel
BuildRequires: zoltan-devel
BuildRequires: cmake3
BuildRequires: %{toolset}-toolchain
BuildRequires: boost-devel python3-devel tbb-devel
BuildRequires: dune-common-devel
BuildRequires: dune-geometry-devel
BuildRequires: dune-uggrid-devel
BuildRequires: dune-grid-devel
BuildRequires: dune-localfunctions-devel
BuildRequires: dune-istl-devel
BuildRequires: opm-common-devel
BuildRequires: opm-material-devel
BuildRequires: opm-grid-devel
BuildRequires: opm-models-devel

%if %{build_openmpi}
BuildRequires: openmpi-devel
BuildRequires: zoltan-openmpi-devel
BuildRequires: dune-common-openmpi-devel
BuildRequires: dune-geometry-openmpi-devel
BuildRequires: dune-uggrid-openmpi-devel
BuildRequires: dune-grid-openmpi-devel
BuildRequires: dune-localfunctions-openmpi-devel
BuildRequires: dune-istl-openmpi-devel
BuildRequires: opm-common-openmpi-devel
BuildRequires: opm-material-openmpi-devel
BuildRequires: opm-grid-openmpi-devel
BuildRequires: opm-models-openmpi-devel
%endif

%if %{build_openmpi3}
BuildRequires: openmpi3-devel
BuildRequires: zoltan-openmpi3-devel
BuildRequires: dune-common-openmpi3-devel
BuildRequires: dune-geometry-openmpi3-devel
BuildRequires: dune-uggrid-openmpi3-devel
BuildRequires: dune-grid-openmpi3-devel
BuildRequires: dune-localfunctions-openmpi3-devel
BuildRequires: dune-istl-openmpi3-devel
BuildRequires: opm-common-openmpi3-devel
BuildRequires: opm-material-openmpi3-devel
BuildRequires: opm-grid-openmpi3-devel
BuildRequires: opm-models-openmpi3-devel
%endif

%if %{build_mpich}
BuildRequires: mpich-devel
BuildRequires: zoltan-mpich-devel
BuildRequires: dune-common-mpich-devel
BuildRequires: dune-geometry-mpich-devel
BuildRequires: dune-uggrid-mpich-devel
BuildRequires: dune-grid-mpich-devel
BuildRequires: dune-localfunctions-mpich-devel
BuildRequires: dune-istl-mpich-devel
BuildRequires: opm-common-mpich-devel
BuildRequires: opm-material-mpich-devel
BuildRequires: opm-grid-mpich-devel
BuildRequires: opm-models-mpich-devel
%endif

BuildRoot:      %{_tmppath}/%{name}-%{version}-build

%description
The Open Porous Media (OPM) initiative provides a set of open-source tools centered around the simulation of flow and transport of fluids in porous media. The goal of the initiative is to establish a sustainable environment for the development of an efficient and well-maintained software suite.

%package -n libopm-simulators%{opm_package_version}
Summary:        Open Porous Media - automatic differentiation library
Group:          System/Libraries

%description -n libopm-simulators%{opm_package_version}
The Open Porous Media (OPM) initiative provides a set of open-source tools centered around the simulation of flow and transport of fluids in porous media. The goal of the initiative is to establish a sustainable environment for the development of an efficient and well-maintained software suite.

%package devel
Summary:        Development and header files for opm-simulators
Group:          Development/Libraries/C and C++
Requires:       libopm-simulator%{opm_package_version}s = %{version}

%description devel
This package contains the development and header files for opm-simulators

%package doc
Summary:        Documentation files for opm-simulators
Group:          Documentation
BuildArch:	noarch

%description doc
This package contains the documentation files for opm-simulators

%package bin%{opm_package_version}
Summary:        Applications in opm-simulators
Group:          Scientific

%description bin%{opm_package_version}
This package contains the applications for opm-simulators

%if %{build_openmpi}
%package -n libopm-simulators-openmpi%{opm_package_version}
Summary:        Open Porous Media - automatic differentiation library
Group:          System/Libraries

%description -n libopm-simulators-openmpi%{opm_package_version}
The Open Porous Media (OPM) initiative provides a set of open-source tools centered around the simulation of flow and transport of fluids in porous media. The goal of the initiative is to establish a sustainable environment for the development of an efficient and well-maintained software suite.

%package openmpi-devel
Summary:        Development and header files for opm-simulators
Group:          Development/Libraries/C and C++
Requires:       libopm-simulators-openmpi%{opm_package_version} = %{version}

%description openmpi-devel
This package contains the development and header files for opm-simulators

%package openmpi-bin%{opm_package_version}
Summary:        Applications in opm-simulators
Group:          Scientific

%description openmpi-bin%{opm_package_version}
This package contains the applications for opm-simulators
%endif

%if %{build_openmpi3}
%package -n libopm-simulators-openmpi3%{opm_package_version}
Summary:        Open Porous Media - automatic differentiation library
Group:          System/Libraries

%description -n libopm-simulators-openmpi3%{opm_package_version}
The Open Porous Media (OPM) initiative provides a set of open-source tools centered around the simulation of flow and transport of fluids in porous media. The goal of the initiative is to establish a sustainable environment for the development of an efficient and well-maintained software suite.

%package openmpi3-devel
Summary:        Development and header files for opm-simulators
Group:          Development/Libraries/C and C++
Requires:       libopm-simulators-openmpi3%{opm_package_version} = %{version}

%description openmpi3-devel
This package contains the development and header files for opm-simulators

%package openmpi3-bin%{opm_package_version}
Summary:        Applications in opm-simulators
Group:          Scientific

%description openmpi3-bin%{opm_package_version}
This package contains the applications for opm-simulators
%endif

%if %{build_mpich}
%package -n libopm-simulators-mpich%{opm_package_version}
Summary:        Open Porous Media - automatic differentiation library
Group:          System/Libraries

%description -n libopm-simulators-mpich%{opm_package_version}
The Open Porous Media (OPM) initiative provides a set of open-source tools centered around the simulation of flow and transport of fluids in porous media. The goal of the initiative is to establish a sustainable environment for the development of an efficient and well-maintained software suite.

%package mpich-devel
Summary:        Development and header files for opm-simulators
Group:          Development/Libraries/C and C++
Requires:       libopm-simulators-mpich%{opm_package_version} = %{version}

%description mpich-devel
This package contains the development and header files for opm-simulators

%package mpich-bin%{opm_package_version}
Summary:        Applications in opm-simulators
Group:          Scientific

%description mpich-bin%{opm_package_version}
This package contains the applications for opm-simulators
%endif

%global debug_package %{nil}

%prep
%setup -q -n %{name}-%{rtype}-%{version}-%{tag}

%build
mkdir serial
pushd serial
scl enable %{toolset} 'cmake3 -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=%{_prefix} -DCMAKE_INSTALL_DOCDIR=share/doc/%{name}-%{version} -DUSE_RUNPATH=OFF -DWITH_NATIVE=OFF -DUSE_QUADMATH=0  -DCMAKE_INSTALL_SYSCONFDIR=/etc .. '
scl enable %{toolset} 'make %{?_smp_mflags}'
scl enable %{toolset} 'make test'
popd

%if %{build_openmpi}
mkdir openmpi
pushd openmpi
module load mpi/openmpi-x86_64
scl enable %{toolset} 'cmake3 -DUSE_MPI=1 -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=%{_prefix}/lib64/openmpi -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_INSTALL_DOCDIR=share/doc/%{name}-%{version} -DUSE_RUNPATH=OFF -DWITH_NATIVE=OFF -DUSE_QUADMATH=0 -DZOLTAN_INCLUDE_DIRS=/usr/include/openmpi-x86_64/zoltan -DCMAKE_INSTALL_SYSCONFDIR=/etc ..'
scl enable %{toolset} 'make %{?_smp_mflags}'
scl enable %{toolset} 'make test'
module unload mpi/openmpi-x86_64
popd
%endif

%if %{build_openmpi3}
mkdir openmpi3
pushd openmpi3
module load mpi/openmpi3-x86_64
scl enable %{toolset} 'cmake3 -DUSE_MPI=1 -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=%{_prefix}/lib64/openmpi3 -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_INSTALL_DOCDIR=share/doc/%{name}-%{version} -DUSE_RUNPATH=OFF -DWITH_NATIVE=OFF -DUSE_QUADMATH=0 -DZOLTAN_INCLUDE_DIRS=/usr/include/openmpi3-x86_64/zoltan -DCMAKE_INSTALL_SYSCONFDIR=/etc ..'
scl enable %{toolset} 'make %{?_smp_mflags}'
scl enable %{toolset} 'make test'
module unload mpi/openmpi3-x86_64
popd
%endif

%if %{build_mpich}
mkdir mpich
pushd mpich
module load mpi/mpich-x86_64
scl enable %{toolset} 'cmake3 -DUSE_MPI=1 -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=%{_prefix}/lib64/mpich -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_INSTALL_DOCDIR=share/doc/%{name}-%{version} -DUSE_RUNPATH=OFF -DWITH_NATIVE=OFF -DUSE_QUADMATH=0 -DZOLTAN_INCLUDE_DIRS=/usr/include/mpich-x86_64/zoltan -DCMAKE_INSTALL_SYSCONFDIR=/etc ..'
scl enable %{toolset} 'make %{?_smp_mflags}'
scl enable %{toolset} 'make test'
module unload mpi/openmpi3-x86_64
popd
%endif

%install
scl enable %{toolset} 'make install DESTDIR=${RPM_BUILD_ROOT} -C serial'
scl enable %{toolset} 'make install-html DESTDIR=${RPM_BUILD_ROOT} -C serial'
mv ${RPM_BUILD_ROOT}/usr/bin/flow ${RPM_BUILD_ROOT}/usr/bin/flow%{opm_package_version}

%if %{build_openmpi}
scl enable %{toolset} 'make install DESTDIR=${RPM_BUILD_ROOT} -C openmpi'
mkdir -p ${RPM_BUILD_ROOT}/usr/include/openmpi-x86_64/
mv ${RPM_BUILD_ROOT}/usr/lib64/openmpi/include/* ${RPM_BUILD_ROOT}/usr/include/openmpi-x86_64/
mv ${RPM_BUILD_ROOT}/usr/lib64/openmpi/bin/flow ${RPM_BUILD_ROOT}/usr/lib64/openmpi/bin/flow%{opm_package_version}
%endif

%if %{build_openmpi3}
scl enable %{toolset} 'make install DESTDIR=${RPM_BUILD_ROOT} -C openmpi3'
mkdir -p ${RPM_BUILD_ROOT}/usr/include/openmpi3-x86_64/
mv ${RPM_BUILD_ROOT}/usr/lib64/openmpi3/include/* ${RPM_BUILD_ROOT}/usr/include/openmpi3-x86_64/
mv ${RPM_BUILD_ROOT}/usr/lib64/openmpi3/bin/flow ${RPM_BUILD_ROOT}/usr/lib64/openmpi3/bin/flow%{opm_package_version}
%endif

%if %{build_mpich}
scl enable %{toolset} 'make install DESTDIR=${RPM_BUILD_ROOT} -C mpich'
mkdir -p ${RPM_BUILD_ROOT}/usr/include/mpich-x86_64/
mv ${RPM_BUILD_ROOT}/usr/lib64/mpich/include/* ${RPM_BUILD_ROOT}/usr/include/mpich-x86_64/
mv ${RPM_BUILD_ROOT}/usr/lib64/mpich/bin/flow ${RPM_BUILD_ROOT}/usr/lib64/mpich/bin/flow%{opm_package_version}
%endif

%clean
rm -rf %{buildroot}

%post -n libopm-simulators%{opm_package_version} -p /sbin/ldconfig
%postun -n libopm-simulators%{opm_package_version} -p /sbin/ldconfig

%if %{build_openmpi}
%post -n libopm-simulators-openmpi%{opm_package_version} -p /sbin/ldconfig
%postun -n libopm-simulators-openmpi%{opm_package_version} -p /sbin/ldconfig
%endif

%if %{build_openmpi3}
%post -n libopm-simulators-openmpi3%{opm_package_version} -p /sbin/ldconfig
%postun -n libopm-simulators-openmpi3%{opm_package_version} -p /sbin/ldconfig
%endif

%if %{build_mpich}
%post -n libopm-simulators-mpich%{opm_package_version} -p /sbin/ldconfig
%postun -n libopm-simulators-mpich%{opm_package_version} -p /sbin/ldconfig
%endif

%files doc
%{_docdir}/*

%files -n libopm-simulators%{opm_package_version}
%defattr(-,root,root,-)
%{_libdir}/*.so.*

%files devel
%defattr(-,root,root,-)
%{_libdir}/*.so
/usr/lib/dunecontrol/*
%{_libdir}/pkgconfig/*
%{_includedir}/*
%{_datadir}/cmake/*
%{_datadir}/opm/cmake/Modules/*

%files bin%{opm_package_version}
%{_bindir}/*
/etc/bash_completion.d/*
%{_datadir}/man/*

%if %{build_openmpi}
%files -n libopm-simulators-openmpi%{opm_package_version}
%defattr(-,root,root,-)
%{_libdir}/openmpi/lib/*.so.*

%files openmpi-devel
%defattr(-,root,root,-)
%{_libdir}/openmpi/lib/*.so
%{_libdir}/openmpi/lib/dunecontrol/*
%{_libdir}/openmpi/lib/pkgconfig/*
%{_includedir}/openmpi-x86_64/*
%{_libdir}/openmpi/share/cmake/*
%{_libdir}/openmpi/share/opm/cmake/Modules/*

%files openmpi-bin%{opm_package_version}
%{_libdir}/openmpi/bin/*
%{_libdir}/openmpi/share/man/*
%endif

%if %{build_openmpi3}
%files -n libopm-simulators-openmpi3%{opm_package_version}
%defattr(-,root,root,-)
%{_libdir}/openmpi3/lib/*.so.*

%files openmpi3-devel
%defattr(-,root,root,-)
%{_libdir}/openmpi3/lib/*.so
%{_libdir}/openmpi3/lib/dunecontrol/*
%{_libdir}/openmpi3/lib/pkgconfig/*
%{_includedir}/openmpi3-x86_64/*
%{_libdir}/openmpi3/share/cmake/*
%{_libdir}/openmpi3/share/opm/cmake/Modules/*

%files openmpi3-bin%{opm_package_version}
%{_libdir}/openmpi3/bin/*
%{_libdir}/openmpi3/share/man/*
%endif

%if %{build_mpich}
%files -n libopm-simulators-mpich%{opm_package_version}
%defattr(-,root,root,-)
%{_libdir}/mpich/lib/*.so.*

%files mpich-devel
%defattr(-,root,root,-)
%{_libdir}/mpich/lib/*.so
%{_libdir}/mpich/lib/dunecontrol/*
%{_libdir}/mpich/lib/pkgconfig/*
%{_includedir}/mpich-x86_64/*
%{_libdir}/mpich/share/cmake/*
%{_libdir}/mpich/share/opm/cmake/Modules/*

%files mpich-bin%{opm_package_version}
%{_libdir}/mpich/bin/*
%{_libdir}/mpich/share/man/*
%endif
