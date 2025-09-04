# Python scripts for testing the generated wheels

## Installation

**Important:** This package requires access to the `python_versions.json` configuration file in the opm-simulators repository. Choose one of these installation methods:

### Option 1: Install within repository (Recommended)
```bash
# Create venv inside the opm-simulators repository
$ cd /path/to/opm-simulators/python/docker/test_wheels/scripts
$ python -m venv .venv
$ source .venv/bin/activate
$ pip install .    # Installs package and dependencies
```

### Option 2: Global installation with environment variable
```bash
# Install anywhere, but set environment variable
$ python -m venv .venv
$ source .venv/bin/activate
$ pip install .
$ export OPM_SIMULATORS_ROOT=/path/to/opm-simulators
```

### For development (editable install with locked dependency versions):
```bash
$ cd /path/to/opm-simulators/python/docker/test_wheels/scripts
$ uv venv
$ uv sync          # Installs package in editable mode with exact versions from uv.lock
$ source .venv/bin/activate
```

## Build the PyPI wheels

Build the wheels from the `python` folder of the `opm-simulators` repository using
the provided helper script for simplified building with automatic logging:

```bash
# Static build with wheel extraction
./docker/run-docker-build.sh static --output wheelhouse
# Creates: docker/wheelhouse/ with wheels

# Shared libraries build
./docker/run-docker-build.sh shared --libtype shared --output wheelhouse-shared
# Creates: docker/wheelhouse-shared/ with wheels
```
see `.docker/run-docker-build.sh` for more examples.

## Build a test container

Build a Docker container for testing the generated wheels:

```bash
# Basic usage (Ubuntu 24.04)
$ opm-wheels build-docker-image --docker-os="ubuntu:24.04"

# Build with custom repository and branch
$ opm-wheels build-docker-image \
    --docker-os="ubuntu:22.04" \
    --opm-simulators-branch="pr-branch" \
    --opm-common-branch="development" \
    --docker-progress

# Build for Debian
$ opm-wheels build-docker-image --docker-os="debian:11"
```

### Available options for `build-docker-image`:
- `--docker-os`, `-d` - OS for Docker image (default: ubuntu:22.04)
  - Supported: `3.8,3.9,3.10,3.11,3.12,3.13`
- `--opm-simulators-repo`, `-osr` - OPM simulators repository URL
- `--opm-common-repo`, `-ocr` - OPM common repository URL
- `--opm-simulators-branch`, `-osb` - Simulators branch (default: master)
- `--opm-common-branch`, `-ocb` - Common branch (default: master)
- `--docker-progress`, `-dp` - Show Docker build progress output

## Run tests

Run unit tests for `opm-common` and `opm-simulators` wheels in the test container:

```bash
# Test all Python versions with default wheel directory
$ opm-wheels run-docker-image --docker-os="ubuntu:24.04"

# Test specific Python versions
$ opm-wheels run-docker-image \
    --docker-os="ubuntu:22.04" \
    --python-versions="3.11,3.12,3.13"

# Test with custom wheel directory
$ opm-wheels run-docker-image \
    --docker-os="debian:11" \
    --wheel-dir="/path/to/wheelhouse-shared" \
    --python-versions="3.12"
```

### Available options for `run-docker-image`:
- `--docker-os`, `-d` - OS for Docker image (default: ubuntu:22.04)
- `--wheel-dir` - Directory containing wheel files (default: python/wheelhouse)
- `--python-versions`, `-p` - Python versions for testing
  - Supported: `3.8,3.9,3.10,3.11,3.12,3.13`
  - Format: comma-separated list without spaces

## Additional commands

### Show wheel contents
```bash
# List files inside a wheel
$ opm-wheels show-wheel-files opm-2025.10-cp312-cp312-manylinux_2_27_x86_64.whl

# Use custom wheel directory
$ opm-wheels show-wheel-files \
    --wheel-dir="wheelhouse-shared" \
    opm_simulators-2025.10-cp312-cp312-manylinux_2_27_x86_64.whl
```

## Troubleshooting

### FileNotFoundError: Could not locate python_versions.json

This error occurs when the CLI cannot find the configuration file. Solutions:

1. **Run from within the repository:**
   ```bash
   cd /path/to/opm-simulators
   cd python/docker/test_wheels/scripts
   source .venv/bin/activate
   opm-wheels run-docker-image --help
   ```

2. **Set environment variable for global installations:**
   ```bash
   export OPM_SIMULATORS_ROOT=/path/to/opm-simulators
   opm-wheels run-docker-image --help
   ```

3. **Verify the configuration file exists:**
   ```bash
   ls /path/to/opm-simulators/python/docker/python_versions.json
   ```

### RuntimeError: Not a valid Git repository

This error means you're not inside a Git repository. Either:
- Navigate to the opm-simulators repository directory, or  
- Use the `OPM_SIMULATORS_ROOT` environment variable method

