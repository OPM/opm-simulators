# Python scripts for testing the generated wheels

## Installation
```
$ python -m venv .venv
$ source .venv/bin/activate
$ pip install .
```
or if you plan to change the source code, install for development like this:
```
$ uv venv
$ uv sync   # Install dependencies
$ source .venv/bin/activate
```

## Build the PyPI wheels

Build the wheels from the `python` folder of the `opm-simulators` repository:

```
$ docker build -t manylinux_2_28_opm:built . -f Dockerfile --output wheelhouse
```

## Build a test container

```
$ opm-wheels build-docker-image --docker-os="ubuntu:24.04"
```

## Run tests

After the wheels and testing containers have been built as shown above, you can run the
`opm-common` and `opm-simulators` unit tests for each wheel in `python/wheelhouse`
in a given container, for example:

```
$ opm-wheels run-docker-image --docker-os="ubuntu:24.04"
```

