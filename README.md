The Open Porous Media Material Framework
========================================

CONTENT
-------

opm-material is an infrastructural OPM module for code related with
material properties like relative-permeability/capillary pressure
laws, thermodynamic relations, flash solvers, empirical heat
conduction laws, et cetera. It is a "library-less" module and only
requires the availability of the DUNE module "dune-common" as well as
the OPM modules "opm-core" and "opm-parser". Historically,
opm-material emerged as a spin-off of OPM's eWoms module [1], which in
turn is a heavily modified version of the Dumux [2] simulation toolkit
for flow and transport in porous media.

LICENSE
-------

This module is distributed under the GNU General Public License,
version 2 or later (GPLv2+). Note, that for these parts the DUNE
template exception does NOT apply.

If you need to know the license of an individual header file, refer to
the copyright statement at the beginning of each file. To avoid legal
issues, the Open Porous Media initiative recommends you to use a
license which is compatible with the GNU General Public license,
version 3 or later (GPLv3+) for code which you publish that uses any
module provided by the Open Porous Media initiative.

PLATFORMS
---------

The opm-material module is designed to run on Linux platforms. It is
also regularly run on Mac OS X. No efforts have been made to ensure
that the code will compile and run on windows platforms, but
contributions for this are certainly welcome.


DEPENDENCIES FOR DEBIAN BASED DISTRIBUTIONS (Debian Squeeze/Ubuntu Precise)
---------------------------------------------------------------------------

# packages necessary for building
sudo apt-get install -y build-essential gfortran cmake cmake-data util-linux

# packages necessary for documentation
sudo apt-get install -y doxygen ghostscript texlive-latex-recommended pgf

# packages necessary for version control
sudo apt-get install -y git-core

# basic libraries
sudo apt-get install -y libboost-all-dev

# required DUNE parts
sudo apt-get install libdune-common-dev libdune-istl-dev

DEPENDENCIES FOR SUSE BASED DISTRIBUTIONS
-----------------------------------------

# repository containing prerequisites
sudo zypper ar http://download.opensuse.org/repositories/science/openSUSE_12.3/science.repo

# utility libraries
sudo zypper in boost-devel

# tools necessary for building
sudo zypper in gcc gcc-c++ gcc-fortran cmake git doxygen

# DUNE libraries
sudo zypper in dune-common-devel dune-istl-devel

DEPENDENCIES FOR RHEL BASED DISTRIBUTIONS
-----------------------------------------

# packages necessary for building
sudo yum install make gcc-c++ gcc-gfortran cmake28 util-linux

# packages necessary for documentation
sudo yum install doxygen ghostscript texlive

# packages necessary for version control
sudo yum install git

# basic libraries
sudo yum install boost-devel

# DUNE libraries
sudo yum-config-manager --add-repo \
    http://www.opm-project.org/packages/current/redhat/6/opm.repo
sudo yum install dune-common-devel dune-istl-devel

DOWNLOADING FROM GIT
--------------------

The prerequisite "dune-common" module can be downloaded like this:

git clone git://github.com/dune-project/dune-common.git

The prerequisite OPM modules are available using the following commands:

git clone git://github.com/OPM/opm-parser.git
git clone git://github.com/OPM/opm-core.git

For a read-only download of the actual opm-material module use:

git clone git://github.com/OPM/opm-material.git

If you want to contribute to the opm-material development, fork
OPM/opm-material on github and open pull requests.

BUILDING
--------

There are two ways to build the opm-material module.

1. As a stand-alone module.
In this setup we recommend creating an entirely separate directory
outside the directory containing the source code and doing the build
from that separate directory (termed "the build directory").  This
configuration is sometimes referred to as an "out-of-source build".

As an example, consider the following layout in which "opm-material" refers
to the directory containing the package source code as downloaded from
GitHub

    workspace
      |
      +-- build
      |
      +-- opm-material
      |     |
      |     +-- ...
      |     |
      |     +-- opm
      |     |
      |     +-- ...

We will configure a release-type (optimised) build using traditional
Unix Makefiles within the "build" directory.  The following command
configures the build

    cd path/to/build
    cmake ../opm-material -DCMAKE_BUILD_TYPE=Release

If you want to debug the library you should specify the build type
"Debug" instead of "Release" in the command above. This will disable
optimizations and make it easier to step through the code.

Building the software then amounts to typing

    make

in the top-level "build" directory; i.e., the directory from which we
invoked the "cmake" utility.  On a multi-core computer system you may
want to build the software in parallel (make(1)'s "job-server" mode) in
order to reduce the total amount of time needed to complete the build.
To do so, replace the above "make" command with

    make -j N

or, possibly,

    nice -20 make -j N

in which "N" is an integer that should typically not exceed the number
of cores in the system.

Once the library has been built, it can be installed in a central,
system-wide location (often in "/usr/local") through the command

    sudo make install


2. As a dune module.
 - Put the opm-material directory in the same directory as the other
   dune modules to be built (e.g. dune-commmon, dune-grid). Note that
   for Ubuntu you can install Dune from the ppa as outlined above.
 - Run dunecontrol normally. For more information on the dune build
   system, see http://www.dune-project.org/doc/installation-notes.html


DOCUMENTATION
-------------

Efforts have been made to document the code with Doxygen.
In order to build the documentation, enter the command

 make doc

in the topmost directory.


REPORTING ISSUES
----------------

Issues can be reported in the Git issue tracker online at:

    http://github.com/OPM/opm-material/issues

To help diagnose build errors, please provide a link to a build log together
with the issue description.

You can capture such a log from the build using the `script' utility, e.g.:

    LOGFILE=$(date +%Y%m%d-%H%M-)build.log ;
        cmake -E cmake_echo_color --cyan --bold "Log file: $LOGFILE" ;
    script -q $LOGFILE -c 'cmake ../opm-material -DCMAKE_BUILD_TYPE=Debug' &&
    script -q $LOGFILE -a -c 'ionice nice make -j 4 -l 3' ||
    cat CMakeCache.txt CMakeFiles/CMake*.log >> $LOGFILE

The resulting file can be uploaded to for instance gist.github.com.

[1] http://opm-project.org/ewoms
[2] http://dumux.org
