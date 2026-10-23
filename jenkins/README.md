# opm-simulators jenkins build scripts:

**build.sh**:
This script will build dependencies, then build opm-simulators and execute its tests.
It also inspects the $ghbPrBuildComment environmental variable and builds
downstreams if requested. It inspects the $ghbPrBuildComment
environmental variable to obtain a pull request to use for the modules.

It is intended for pre-merge builds of pull requests.

To specify a given pull request to use for upstreams and downstreams,
trigger line needs to contain &lt;module-name&gt;=&lt;pull request number&gt;.

To build with downstreams the trigger line needs to contain 'with downstreams'.

**static_analysis.sh**:
Run static analysis for opm-simulators

**run-spe.sh**:
This script will execute the SPE1, SPE3 and SPE9 cases, then compare against
OPM and Eclipse reference results. It is meant to be executed after a build.
The binary used is from the build identified by the configuration
environment variable.

**run-norne.sh**:
This script will execute the Norne case. It is meant to be executed after
a build. The binary used is from the build identified by the configuration
given as the first argument, and the case is run in parallel if the number
of processes is given as the second argument.

**plot-norne.sh**:
This script generates documents with plots of the Norne wells, after the
runs by run-norne.sh. They include the serial run and, if the number of
processes is given as the second argument, the parallel run.
