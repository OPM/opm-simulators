# Build and install dune from source

- Install the build dependencies, e.g. on Ubuntu:
```
$ sudo apt-get build-dep dune-common dune-geometry dune-istl dune-grid
```
See also:
- https://opm-project.org/?page_id=239
- https://dune-project.org/installation/installation-buildsrc/

# Determine the version of dune to install

Example: On Ubuntu, you can run the following command:

```
$ apt-cache showsrc dune-common | grep ^Version:
```

# Run the install script:

```
$ ./install.sh --version=2.10.0 --prefix=/opt/dune-2.10.0 --use-mpi=no --build-type=debug --use-sudo=yes
```

# Use with OPM flow

## Serial build (without MPI)

As an example, assuming the above `install.sh` command to first install dune in `/opt/dune-2.10.0`: You can give the follwing cmake flags when building opm-common, opm-grid, and opm-simulators to build without MPI support:
```
-DCMAKE_PREFIX_PATH=/opt/dune-2.10.0 -DUSE_MPI=OFF -DCMAKE_DISABLE_FIND_PACKAGE_MPI=ON
```
