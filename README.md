# Open Porous Media Simulators and Automatic Differentiation Library

CONTENT
-------

opm-simulators contains the OPM Flow reservoir simulator, which is a fully
implicit black-oil simulator that also supports CO2 storage, H2 storage, thermal,
solvent, and polymer options. It supports input and output in the common Eclipse formats,
allowing easy integration in existing workflows.
Its system assembly approach builds on automatic differentiation,
using the local AD class Evaluation from opm-common. The linear solver subsystem is based on
the dune-istl library.

For more information see https://opm-project.org, for the OPM manual see
https://opm-project.org/?page_id=955

LICENSE
-------

The library is distributed under the GNU General Public License,
version 3 or later (GPLv3+).


PLATFORMS
---------

The opm-simulators module is designed to run on Linux platforms. It is
also regularly run on Mac OS X. No efforts have been made to ensure
that the code will compile and run on windows platforms.


REQUIREMENTS
------------

opm-simulators requires several other OPM modules, see
http://opm-project.org/?page_id=274. In addition, opm-simulators
requires Dune and some other software to be available, for details
see https://opm-project.org/?page_id=239.


DOWNLOADING
-----------

For a read-only download:
git clone git://github.com/OPM/opm-simulators.git

If you want to contribute, fork OPM/opm-simulators on github.


BUILDING
--------

See build instructions at http://opm-project.org/?page_id=36


IN-CODE DOCUMENTATION
---------------------

Efforts have been made to document the code with Doxygen.
In order to build the documentation, enter the command

 make doc

in the topmost directory.


REPORTING ISSUES
----------------

The OPM mailing list can be used for questions and issue reporting,
see https://opm-project.org/?page_id=358

Issues can also be reported in the Git issue tracker online at:

    https://github.com/OPM/opm-simulators/issues

To help diagnose build errors, please provide a link to a build log together
with the issue description.

You can capture such a log from the build using the `script' utility, e.g.:

    LOGFILE=$(date +%Y%m%d-%H%M-)build.log ;
    cmake -E cmake_echo_color --cyan --bold "Log file: $LOGFILE" ;
    script -q $LOGFILE -c 'cmake ../opm-core -DCMAKE_BUILD_TYPE=Debug' &&
    script -q $LOGFILE -a -c 'ionice nice make -j 4 -l 3' ||
    cat CMakeCache.txt CMakeFiles/CMake*.log >> $LOGFILE

The resulting file can be uploaded to for instance gist.github.com.


CONTRIBUTING
------------

Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on how to contribute to this project.
