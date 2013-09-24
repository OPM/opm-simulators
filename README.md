Open Porous Media Automatic Differentiation Library
===================================================


CONTENT
-------

opm-autodiff contains a small library for automatic differentiation
built on the Eigen linear algebra package. It also contains some
experimental solver and simulator prototypes demonstrating how it can
be used. The most important parts are:

* AutoDiffBlock.hpp (class for AD on vectorized data with sparse jacobians)
* AutoDiffHelpers.hpp (various utilities to simplify writing solvers)
* sim_fibo_ad.cpp (prototype fully implicit black-oil simulator)

LICENSE
-------

The library is distributed under the GNU General Public License,
version 3 or later (GPLv3+).


PLATFORMS
---------

The opm-autodiff module is designed to run on Linux platforms. It is
also regularly run on Mac OS X. No efforts have been made to ensure
that the code will compile and run on windows platforms.


REQUIREMENTS
------------

opm-autodiff requires opm-core, and all its requirements (see
opm-core/README). In addition, opm-autodiff requires Eigen, version
3.1 (has not been tested with later versions).


DOWNLOADING
-----------

For a read-only download:
git clone git://github.com/OPM/opm-autodiff.git

If you want to contribute, fork OPM/opm-autodiff on github.


BUILDING
--------

See build instructions in opm-core/README, or at
http://opm-project.org/download.php.


DOCUMENTATION
-------------

Efforts have been made to document the code with Doxygen.
In order to build the documentation, enter the command

 make doc

in the topmost directory. The class AutoDiffBlock is the most
important and most well-documented.


REPORTING ISSUES
----------------

Issues can be reported in the Git issue tracker online at:

    https://github.com/OPM/opm-autodiff/issues

To help diagnose build errors, please provide a link to a build log together
with the issue description.

You can capture such a log from the build using the `script' utility, e.g.:

    LOGFILE=$(date +%Y%m%d-%H%M-)build.log ;
    cmake -E cmake_echo_color --cyan --bold "Log file: $LOGFILE" ;
    script -q $LOGFILE -c 'cmake ../opm-core -DCMAKE_BUILD_TYPE=Debug' &&
    script -q $LOGFILE -a -c 'ionice nice make -j 4 -l 3' ||
    cat CMakeCache.txt CMakeFiles/CMake*.log >> $LOGFILE

The resulting file can be uploaded to for instance gist.github.com.
